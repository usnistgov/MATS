from bisect import bisect
import re
import warnings

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from .hapi import ISO, PYTIPS2017, PYTIPS2011, PYTIPS2021, PYTIPS2025, pcqsdhc, PROFILE_LORENTZ
from .utilities import molecularMass, etalon, convolveSpectrumSame
from .codata import CONSTANTS
from .o2_cia_karman import O2_CIA_Karman_Model
from numba import jit, prange
from .mHT.profile import mHTprofile_vector
from .mHT.numba_mHTP import mHTprofile_vector_numba

@jit(nopython=True, parallel=True, fastmath=True, cache=True)
def fast_broadening_params(nlines, p, T, pref, Tref,
                           abundances, 
                           # Arrays for Gamma0
                           g0_arrs, n_g0_arrs,
                           # Arrays for Delta0
                           d0_arrs, n_d0_arrs,
                           # Arrays for Gamma2 (SD)
                           sd_gamma_arrs, n_gamma2_arrs,
                           # Arrays for Delta2 (SD)
                           sd_delta_arrs, n_delta2_arrs,
                           
                           # Flexible "Real" Part (nuVC or NuOptRe)
                           param_re_arrs, n_param_re_arrs,
                           # Flexible "Imag" Part (eta or NuOptIm)
                           param_im_arrs, n_param_im_arrs,
                           # Arrays for Y (Line Mixing)
                           y_arrs, n_y_arrs, 
                           lineprofile_mode):
    """
    Calculates weighted sums of broadening parameters for all lines in parallel.
    Avoids allocating temporary arrays for every intermediate math step.
    """
    
    # Initialize output arrays
    Gamma0 = np.zeros(nlines, dtype=np.float64)
    Delta0 = np.zeros(nlines, dtype=np.float64)
    Gamma2 = np.zeros(nlines, dtype=np.float64)
    Delta2 = np.zeros(nlines, dtype=np.float64)
    # Generic outputs (Re=Real, Im=Imaginary) # This is loose nomenclatrue 
    Out_Re = np.zeros(nlines, dtype=np.float64) 
    Out_Im = np.zeros(nlines, dtype=np.float64)
    #mHTP
    Y      = np.zeros(nlines, dtype=np.float64)
    
    # Pre-calculate T ratio
    T_ratio = Tref / T
    delta_T = T - Tref
    p_norm = p / pref
    
    # Loop over diluents (outer loop, usually small, e.g. 2 for Air)
    num_diluents = len(abundances)
    
    for d in range(num_diluents):
        abun = abundances[d]
        
        # Inner loop over lines (Parallelized)
        for i in prange(nlines):
            
            # Gamma0
            term_g0 = g0_arrs[d][i] * (T_ratio ** n_g0_arrs[d][i])
            Gamma0[i] += abun * term_g0 * p_norm
            
            # Delta0
            Delta0[i] += abun * (d0_arrs[d][i] + n_d0_arrs[d][i] * delta_T) * p_norm
            
            # Gamma2 (Speed Dependent)
            Gamma2[i] += abun * (sd_gamma_arrs[d][i] * g0_arrs[d][i] * p_norm * (T_ratio ** n_gamma2_arrs[d][i]))
            
            # Delta2 (Speed Dependent)
            Delta2[i] += abun * ((sd_delta_arrs[d][i] * d0_arrs[d][i] + n_delta2_arrs[d][i] * delta_T) * p_norm)

            # NuVC (Dicke)
            Out_Re[i] += abun * (param_re_arrs[d][i] * p_norm * (T_ratio ** n_param_re_arrs[d][i]))

            if lineprofile_mode == 0:
                # Eta
                Out_Im[i] += param_im_arrs[d][i] * abun
            else: 
                Out_Im[i] += abun * (param_im_arrs[d][i] * p_norm * (T_ratio ** n_param_im_arrs[d][i]))

            # Line Mixing
            Y[i] += abun * (y_arrs[d][i] * p_norm * (T_ratio ** n_y_arrs[d][i]))
            
    return Gamma0, Delta0, Gamma2, Delta2, Out_Re, Out_Im, Y

class Spectroscopic_model:
    #converts from pandas DF to numpy arrays
    # LBL calcultion
    def __init__(self, parameter_linelist, isotope_list = ISO, 
                 lineprofile = 'mHTP', numba_lineprofile = True,
                 beta_formalism = False):
        self.nlines = len(parameter_linelist)
        self.lineprofile = lineprofile
        if lineprofile != 'HTP':
            self.beta_formalism = beta_formalism
            self.numba_lineprofile = numba_lineprofile
        else:
            self.beta_formalism = False
            self.numba_lineprofile = False

        #Static arrays
        self.elower = parameter_linelist['elower'].to_numpy(dtype=np.float64)
        self.molec_id = parameter_linelist['molec_id'].to_numpy(dtype=np.int32)   
        self.local_iso_id = parameter_linelist['local_iso_id'].to_numpy(dtype=np.int32)
        self.sw_scale_factor = parameter_linelist['sw_scale_factor'].to_numpy(dtype=np.float64)

        self.global_to_local_idx = {g_idx: l_idx for l_idx, g_idx in enumerate(parameter_linelist.index)}

        if self.lineprofile == 'HTP':
            # HTP: Real = nuVC, Imag = eta (Eta has no T-dep in HTP)
            self.col_map = {
                'param_Re': 'nuVC_', 'n_param_Re': 'n_nuVC_',
                'param_Im': 'eta_',  'n_param_Im': 'n_eta_' # Map even if unused, will default to 0
            }
        else:
            # mHTP: Real = NuOptRe, Imag = NuOptIm (Both have T-dependence)
            self.col_map = {
                'param_Re': 'nuOptRe_', 'n_param_Re': 'n_nuOptRe_',
                'param_Im': 'nuOptIm_', 'n_param_Im': 'n_nuOptIm_'
            }


        #Dynamic
        self.dynamic_arrays = {}
        for col in parameter_linelist.columns:
            # Check if this column is one of our special mapped ones
            is_mapped = False
            for key in ['param_Re', 'n_param_Re', 'param_Im', 'n_param_Im']:
                prefix = self.col_map[key]
                if col.startswith(prefix):
                    diluent = col[len(prefix):] 
                    internal_name = f"{key}_{diluent}" # e.g. re_air
                    data_array = parameter_linelist[col].values.astype(np.float64)
                    self.dynamic_arrays[internal_name] = data_array
                    self.dynamic_arrays[col] = data_array
                    is_mapped = True
            
            # If not mapped, check if it is a standard parameter
            if not is_mapped:
                if col.startswith(('nu', 'sw', 'gamma0_', 'n_gamma0_', 'delta0_', 'n_delta0_', 
                                   'SD_gamma_', 'n_gamma2_', 'SD_delta_', 'n_delta2_', 
                                   'y_', 'n_y', 'BIA_slope_', 'BIA_collision_duration_')):
                    self.dynamic_arrays[col] = parameter_linelist[col].values.astype(np.float64)

        self.mass = np.zeros(self.nlines, dtype=np.float64)
        self.sigma_Tref = np.zeros(self.nlines, dtype=np.float64) # Tref is fixed at 296K
        #self.abundance_ratio = np.ones(self.nlines, dtype=np.float64) 
        self.unique_pairs = np.unique(np.column_stack((self.molec_id, self.local_iso_id)), axis=0)

        for m, i in self.unique_pairs:
            m, i = int(m), int(i)
            mask = (self.molec_id == m) & (self.local_iso_id == i)
            
            # Mass
            self.mass[mask] = molecularMass(m, i, isotope_list=isotope_list)

        #MAPPING CACHE (Filled by the fit_dataset later)
        self.param_index_map = []
        #Caching State
        # SMART CACHING STATE VARIABLES
        self.active_line_indices = np.array([], dtype=np.int32)
        self.static_line_indices = np.array([], dtype=np.int32)
        self.cached_static_spectrum = {}
        self.cache_state = {}

        
        
    def _resolve_array(self, param_base_name, spectrum_number):
        """
        The Waterfall Lookup.
        Tries to find specific '{param}_{id}', then falls back to '{param}'.
        
        Parameters
        ----------
        param_base_name : str
            The generic name of the parameter (e.g., 'sw', 'gamma0_air').
        spectrum_id : int or str
            The ID of the current spectrum being simulated (e.g., 5).
            
        Returns
        -------
        np.array or None
            The requested array if found, otherwise None.
        """
        # 1. Try Specific (e.g., 'sw_5' or 'gamma0_air_5')
        # This allows you to fit a parameter for just ONE spectrum.
        specific_key = f"{param_base_name}_{spectrum_number}"
        if specific_key in self.dynamic_arrays:
            return self.dynamic_arrays[specific_key]
        
        if param_base_name in self.dynamic_arrays:
            return self.dynamic_arrays[param_base_name]
        
        # FAIL-SAFE: Return zeros to satisfy Numba types
        return np.zeros(self.nlines, dtype=np.float64)

    def calculate_spectrum(self, waves, 
                            T, p, molefraction, abundance_ratios, Diluent, spectrum_number, spectrum_min, 
                            baseline_coeffs = None, 
                            etalon_dict = None, 
                            cia_config = None, 
                            ILS_function = None, ILS_parameters = None, ILS_wing = None,
                            interpolated_compressability_file=None,
                            #pass through for LBL_alpha
                             TIPS = PYTIPS2025, isotope_list = ISO, 
                            natural_abundance = True, abundance_ratio_MI = {}, 
                            BIA_slope=False, BIA_FW_LBL=False,
                            IntensityThreshold = 1e-30, 
                            wing_cutoff = 25, wing_wavenumbers = 25, wing_method = 'wing_wavenumbers',segment = None, 
                            pathlength = 0,
                            dataspace = 'alpha'):
        
        #Calculate LBL
        LBL_alpha = self.calculate_lbl_absorbance(waves, T, p, molefraction, abundance_ratios, Diluent, 
                                 spectrum_number, 
                                 interpolated_compressability_file = interpolated_compressability_file, TIPS = TIPS, isotope_list = isotope_list, 
                                 natural_abundance = natural_abundance, 
                                 BIA_slope=BIA_slope, BIA_FW_LBL=BIA_FW_LBL,
                                 IntensityThreshold = IntensityThreshold, 
                                 wing_cutoff = wing_cutoff, wing_wavenumbers = wing_wavenumbers, wing_method = wing_method,segment=segment)      
        
        if dataspace == 'alpha':
            LBL_alpha = LBL_alpha* 1e6

        # Calculate CIA
        cia_alpha = np.zeros_like(waves)
        if cia_config:
            if cia_config['model'] == 'Karman':
                    cia_alpha = cia_config['calculator'].calculate_cia(
                        waves, T, p, Diluent, 
                        **cia_config['params'])
            if dataspace == 'alpha':
                pass  
                    
            if cia_config['model'] == 'ad hoc':
                cia_alpha = cia_config['values']

        #Baseline
        baseline = np.zeros_like(waves)
        w_rel = waves - spectrum_min
        if baseline_coeffs is not None:
            baseline = np.polyval(baseline_coeffs, w_rel)

        #Etalon Calculation
        etalon_signal = np.zeros_like(waves)
        if etalon_dict is not None:
            for i, params in etalon_dict.items():
                etalon_signal += etalon(w_rel, params['amp'], params['period'], params['phase'])
        
        if dataspace == 'alpha':
            model_y = LBL_alpha + cia_alpha + baseline + etalon_signal
        elif dataspace == 'absorbance':
            model_y = (LBL_alpha + cia_alpha)*pathlength + (baseline + etalon_signal)
        elif (dataspace == 'absorption')  or (dataspace == 'transmittance'):
            OD = (LBL_alpha + cia_alpha)*pathlength
            model_y = np.exp(-OD)*(1+baseline + etalon_signal)


        if ILS_function is not None and ILS_parameters is not None:
            # Note: convolveSpectrumSame returns (waves, alpha, ...)
            # We ignore the extra returns for the model result
            res_to_pass = ILS_parameters
            if (isinstance(ILS_parameters, list) or isinstance(ILS_parameters, np.ndarray)) and len(ILS_parameters) == 1:
                res_to_pass = ILS_parameters[0]
                
            _, model_y, _, _, _ = convolveSpectrumSame(
                waves, model_y, 
                SlitFunction=ILS_function, 
                Resolution=res_to_pass, 
                AF_wing=ILS_wing)
        if dataspace == 'absorption':
            model_y = 1-model_y

        return model_y
    
    @staticmethod
    @jit(nopython=True, fastmath=True, cache=True)
    def calculate_lbl_numba_kernel(waves, 
                                   nu_arr, intensity_arr, 
                                   gamma0_arr, gamma2_arr, delta0_arr, delta2_arr, 
                                   re_arr, im_arr, y_arr, 
                                   gammaD_arr, alpha_arr, 
                                   idx_low_arr, idx_high_arr,
                                   mol_dens, out_spectrum):
        
        n_lines = len(nu_arr)
        
        # Allocate ONCE, reuse for all lines.
        MAX_POINTS = 20000 
        scratch_buffer = np.zeros(MAX_POINTS, dtype=np.float64)
        
        for i in range(n_lines):
            idx_start = idx_low_arr[i]
            idx_end = idx_high_arr[i]
            
            slice_len = idx_end - idx_start
            if slice_len <= 0: continue
            
            if slice_len > MAX_POINTS:
                line_buffer = np.zeros(slice_len, dtype=np.float64)
            else:
                line_buffer = scratch_buffer[:slice_len]
            
            wave_slice = waves[idx_start:idx_end]
            
            mHTprofile_vector_numba(
                nu_arr[i], gammaD_arr[i], 
                gamma0_arr[i], gamma2_arr[i], delta0_arr[i], delta2_arr[i],
                re_arr[i], im_arr[i], 
                wave_slice, 
                y_arr[i], 0.0, alpha_arr[i], 0.0,
                line_buffer 
            )
            
            factor = mol_dens * intensity_arr[i]
            
            for k in range(slice_len):
                # ADD IN-PLACE TO THE FED ARRAY
                out_spectrum[idx_start + k] += factor * line_buffer[k]
            
        return out_spectrum
    


    def calculate_lbl_absorbance(self, waves, T, p, molefraction, abundance_ratios, Diluent, spectrum_number, 
                                 interpolated_compressability_file = None, TIPS = PYTIPS2025, isotope_list = ISO, 
                                 natural_abundance = True, #abundance_ratio_MI = {}, 
                                 BIA_slope=False, BIA_FW_LBL=False,
                                 IntensityThreshold = 1e-30, 
                                 wing_cutoff = 25, wing_wavenumbers = 25, wing_method = 'wing_wavenumbers',segment = None):

        Tref = 296. 
        pref = 1. 
        mol_density_ref = (pref/ CONSTANTS['cpa_atm'])/(CONSTANTS['k']*273.15) 
        mol_dens = (p/ CONSTANTS['cpa_atm'])/(CONSTANTS['k']*T)
        if interpolated_compressability_file != None:
            mol_dens = mol_dens / interpolated_compressability_file([p, T])[0]
        density_amagat = mol_dens / mol_density_ref

        sigma_T = np.ones(self.nlines, dtype=np.float64)
        sigma_Tref = np.ones(self.nlines, dtype=np.float64)
        #abundance_ratio = np.ones(self.nlines, dtype=np.float64)
        for m, i in self.unique_pairs:
            m, i = int(m), int(i)
            mask = (self.molec_id == m) & (self.local_iso_id == i)
            try:
                sigma_T[mask] = TIPS(m, i, T)
                sigma_Tref[mask] = TIPS(m, i, Tref)
                #if (not natural_abundance) and (abundance_ratio_MI != {}):
                    #abundance_ratio[mask] = abundance_ratio_MI[m][i]
            except:
                pass

        sw_array = self._resolve_array('sw', spectrum_number)
        sw_array = sw_array*self.sw_scale_factor
        nu_array = self._resolve_array('nu', spectrum_number)
        GammaD = np.sqrt(2*CONSTANTS['k']*CONSTANTS['Na']*T*np.log(2)/(self.mass))*nu_array / CONSTANTS['c']
        
        line_intensity = sw_array * (sigma_Tref / sigma_T) * \
                    (np.exp(-CONSTANTS['c2'] * self.elower / T) * (1 - np.exp(-CONSTANTS['c2'] * nu_array / T))) / \
                    (np.exp(-CONSTANTS['c2'] * self.elower / Tref) * (1 - np.exp(-CONSTANTS['c2'] * nu_array / Tref)))
        
        abundances_list = []
        g0_list = []; n_g0_list = []; d0_list = []; n_d0_list = []
        sd_gamma_list = []; n_gamma2_list = []; sd_delta_list = []; n_delta2_list = []
        param_re_list = []; n_param_re_list = []; param_im_list = []; n_param_im_list = []
        y_list = []; n_y_list = []

        perturber_mass = 0
        for species, info in Diluent.items():
            abundances_list.append(info['composition'])
            perturber_mass += Diluent[species]['composition']*Diluent[species]['m']
            
            g0_list.append(self._resolve_array(f'gamma0_{species}', spectrum_number))
            n_g0_list.append(self._resolve_array(f'n_gamma0_{species}', spectrum_number))
            d0_list.append(self._resolve_array(f'delta0_{species}', spectrum_number))
            n_d0_list.append(self._resolve_array(f'n_delta0_{species}', spectrum_number))
            sd_gamma_list.append(self._resolve_array(f'SD_gamma_{species}', spectrum_number))
            n_gamma2_list.append(self._resolve_array(f'n_gamma2_{species}', spectrum_number))
            sd_delta_list.append(self._resolve_array(f'SD_delta_{species}', spectrum_number))
            n_delta2_list.append(self._resolve_array(f'n_delta2_{species}', spectrum_number))
            param_re_list.append(self._resolve_array(f'param_Re_{species}', spectrum_number))
            n_param_re_list.append(self._resolve_array(f'n_param_Re_{species}', spectrum_number))
            param_im_list.append(self._resolve_array(f'param_Im_{species}', spectrum_number))
            n_param_im_list.append(self._resolve_array(f'n_param_Im_{species}', spectrum_number))
            y_list.append(self._resolve_array(f'y_{species}', spectrum_number))
            n_y_list.append(self._resolve_array(f'n_y_{species}', spectrum_number))
        
        if self.beta_formalism:
            alpha = perturber_mass / self.mass
        else:
            alpha = np.ones(self.nlines, dtype=np.float64)*10.0
        
        lineprofile_mode = 0 if self.lineprofile == 'HTP' else 1

        Gamma0, Delta0, Gamma2, Delta2, ParamRe, ParamIm, Y = fast_broadening_params(
            self.nlines, p, T, pref, Tref,
            np.array(abundances_list, dtype=np.float64),
            np.array(g0_list, dtype=np.float64), np.array(n_g0_list, dtype=np.float64),
            np.array(d0_list, dtype=np.float64), np.array(n_d0_list, dtype=np.float64),
            np.array(sd_gamma_list, dtype=np.float64), np.array(n_gamma2_list, dtype=np.float64),
            np.array(sd_delta_list, dtype=np.float64), np.array(n_delta2_list, dtype=np.float64),
            np.array(param_re_list, dtype=np.float64), np.array(n_param_re_list, dtype=np.float64),
            np.array(param_im_list, dtype=np.float64), np.array(n_param_im_list, dtype=np.float64),
            np.array(y_list, dtype=np.float64), np.array(n_y_list, dtype=np.float64), lineprofile_mode)

        intensity_bia = np.zeros(self.nlines)
        bia_duration = np.zeros(self.nlines)
        if BIA_slope:
            for species, info in Diluent.items():
                abun = info['composition']
                bia_slope_arr = self._resolve_array(f'BIA_slope_{species}', spectrum_number)
                intensity_bia += abun*(line_intensity*(1-0.01*bia_slope_arr*(density_amagat)))
                if BIA_FW_LBL:
                    bia_dur_arr = self._resolve_array(f'BIA_collision_duration_{species}', spectrum_number)
                    bia_duration += abun*(bia_dur_arr)
        
        if wing_method == 'wing_wavenumbers':
            line_cutoffs = np.ones(self.nlines)*wing_wavenumbers
        elif wing_method == 'wing_cutoff':
            line_cutoffs = (0.5346 * Gamma0 + np.sqrt(0.2166 * Gamma0**2 + GammaD**2)) * wing_cutoff

        valid_indices = np.where(line_intensity >= IntensityThreshold)[0]
        if len(valid_indices) == 0:
            return np.zeros_like(waves)

        valid_nu = nu_array[valid_indices]
        valid_cut = line_cutoffs[valid_indices]
        idx_low_arr = np.searchsorted(waves, valid_nu - valid_cut)
        idx_high_arr = np.searchsorted(waves, valid_nu + valid_cut)

        cache_key = f"{spectrum_number}_{segment}" if segment is not None else spectrum_number
        has_static = (hasattr(self, 'static_line_indices') and len(self.static_line_indices) > 0)
        
        mf_state = tuple(sorted(molefraction.items()))
        abund_state = tuple(sorted(abundance_ratios.items()))
        current_state = {
            'P': p, 
            'T': T, 
            'waves_start': waves[0], 
            'waves_end': waves[-1], 
            'waves_len': len(waves), 
            'mf': mf_state,
            'abund': abund_state  
        }
        cache_valid = (cache_key in self.cache_state and self.cache_state[cache_key] == current_state)

        # =========================================================
        # FAST ARRAY VECTORIZATION (Fixes Dictionary Bug)
        # =========================================================
        mf_array = np.zeros(self.nlines, dtype=np.float64)
        abund_array = np.ones(self.nlines, dtype=np.float64)
        
        for m, i in self.unique_pairs:
            m_int, i_int = int(m), int(i)
            mask = (self.molec_id == m_int) & (self.local_iso_id == i_int)
            if m in molefraction:
                mf_array[mask] = molefraction[m_int]
            if (m_int, i_int) in abundance_ratios:
                abund_array[mask] = abundance_ratios[(m_int, i_int)]
                
        # Calculate final effective intensity for all lines at once
        eff_intensity = np.zeros(self.nlines, dtype=np.float64)
        if BIA_slope:
            eff_intensity[valid_indices] = intensity_bia[valid_indices] * mf_array[valid_indices] * abund_array[valid_indices]
        else:
            eff_intensity[valid_indices] = line_intensity[valid_indices] * mf_array[valid_indices] * abund_array[valid_indices]

        # =========================================================
        #  PATH A: NUMBA KERNEL (mHTP-Numba)
        # =========================================================
        if self.lineprofile != 'HTP' and self.numba_lineprofile:
            
            # 1. Update the Static Cache if needed
            if has_static and not cache_valid:
                static_to_calc = np.intersect1d(valid_indices, self.static_line_indices)
                static_spectrum = np.zeros_like(waves)
                
                if len(static_to_calc) > 0:
                    self.calculate_lbl_numba_kernel(
                        waves, nu_array[static_to_calc], eff_intensity[static_to_calc],
                        Gamma0[static_to_calc], Gamma2[static_to_calc], Delta0[static_to_calc], Delta2[static_to_calc],
                        ParamRe[static_to_calc], ParamIm[static_to_calc], Y[static_to_calc], GammaD[static_to_calc], alpha[static_to_calc],
                        idx_low_arr[np.searchsorted(valid_indices, static_to_calc)], 
                        idx_high_arr[np.searchsorted(valid_indices, static_to_calc)],
                        mol_dens, static_spectrum
                    )
                self.cached_static_spectrum[cache_key] = static_spectrum
                self.cache_state[cache_key] = current_state

            # 2. Feed the Cache (or zeros) into the Final Spectrum
            if has_static and cache_key in self.cached_static_spectrum:
                final_spectrum = self.cached_static_spectrum[cache_key].copy()
            else:
                final_spectrum = np.zeros_like(waves)

            # 3. Add Active Lines In-Place directly into final_spectrum
            active_to_calc = np.intersect1d(valid_indices, self.active_line_indices) if has_static else valid_indices
            
            if len(active_to_calc) > 0:
                self.calculate_lbl_numba_kernel(
                    waves, nu_array[active_to_calc], eff_intensity[active_to_calc],
                    Gamma0[active_to_calc], Gamma2[active_to_calc], Delta0[active_to_calc], Delta2[active_to_calc],
                    ParamRe[active_to_calc], ParamIm[active_to_calc], Y[active_to_calc], GammaD[active_to_calc], alpha[active_to_calc],
                    idx_low_arr[np.searchsorted(valid_indices, active_to_calc)], 
                    idx_high_arr[np.searchsorted(valid_indices, active_to_calc)],
                    mol_dens, final_spectrum # <-- Numba stacks active lines right on top!
                )
            
            # 4. Handle Far Wing Continuum
            if BIA_slope and BIA_FW_LBL:
                for k, i in enumerate(valid_indices):
                    if bia_duration[i] != 0:
                        cut_i = line_cutoffs[i]
                        BoundIndexLower_BIA = bisect(waves, nu_array[i] - 5*cut_i)
                        BoundIndexUpper_BIA = bisect(waves, nu_array[i] + 5*cut_i)
                        wave_slice_BIA = waves[BoundIndexLower_BIA:BoundIndexUpper_BIA]
                        BIA_profile = PROFILE_LORENTZ(nu_array[i], 1/(2*np.pi*CONSTANTS['c']*1e-12*bia_duration[i]), 0, wave_slice_BIA)
                        
                        val = (line_intensity[i] - intensity_bia[i])
                        final_spectrum[BoundIndexLower_BIA:BoundIndexUpper_BIA] += mol_dens * mf_array[i] * abund_array[i] * val * BIA_profile
            
            return final_spectrum

        # =========================================================
        #  PATH B: PYTHON LOOP (HTP/mHTP Standard)
        # =========================================================
        else:
            calc_groups = []
            if has_static and not cache_valid:
                calc_groups.append(('static', np.intersect1d(valid_indices, self.static_line_indices)))
            
            if has_static:
                calc_groups.append(('active', np.intersect1d(valid_indices, self.active_line_indices)))
            else:
                calc_groups.append(('active', valid_indices))

            results = {}
            for label, indices in calc_groups:
                group_xsect = np.zeros_like(waves)
                if len(indices) > 0:
                    grp_nu = nu_array[indices]
                    grp_cut = line_cutoffs[indices]
                    grp_low = np.searchsorted(waves, grp_nu - grp_cut)
                    grp_high = np.searchsorted(waves, grp_nu + grp_cut)
                    
                    for k, i in enumerate(indices):
                        idx_low = max(0, grp_low[k])
                        idx_high = min(len(waves), grp_high[k])
                        if idx_low >= idx_high: continue

                        wave_slice = waves[idx_low:idx_high]
                        
                        if self.lineprofile == 'HTP':
                            lineshape_PT, _ = pcqsdhc(nu_array[i], GammaD[i], Gamma0[i], Gamma2[i], Delta0[i], Delta2[i],
                                                      ParamRe[i], ParamIm[i], wave_slice, Ylm=Y[i])
                        else:
                            lineshape_PT = mHTprofile_vector(nu_array[i], GammaD[i], Gamma0[i], Gamma2[i], Delta0[i], Delta2[i], 
                                                             ParamRe[i], ParamIm[i], wave_slice, Y[i], 0, alpha[i])
                        
                        group_xsect[idx_low:idx_high] += mol_dens * eff_intensity[i] * lineshape_PT
                            
                results[label] = group_xsect

            if 'static' in results:
                self.cached_static_spectrum[cache_key] = results['static']
                self.cache_state[cache_key] = current_state
            
            active_spectrum = results.get('active', np.zeros_like(waves))
            
            if has_static and cache_key in self.cached_static_spectrum:
                final_spectrum = active_spectrum + self.cached_static_spectrum[cache_key]
            else:
                final_spectrum = active_spectrum

            if BIA_slope and BIA_FW_LBL:
                for i in valid_indices:
                    if bia_duration[i] != 0:
                        cut_i = line_cutoffs[i]
                        BoundIndexLower_BIA = bisect(waves, nu_array[i] - 5*cut_i)
                        BoundIndexUpper_BIA = bisect(waves, nu_array[i] + 5*cut_i)
                        wave_slice_BIA = waves[BoundIndexLower_BIA:BoundIndexUpper_BIA]
                        BIA_profile = PROFILE_LORENTZ(nu_array[i], 1/(2*np.pi*CONSTANTS['c']*1e-12*bia_duration[i]), 0, wave_slice_BIA)
                        
                        val = (line_intensity[i] - intensity_bia[i])
                        final_spectrum[BoundIndexLower_BIA:BoundIndexUpper_BIA] += mol_dens * mf_array[i] * abund_array[i] * val * BIA_profile

            return final_spectrum
    
    def configure_for_fitting(self, lmfit_params):
        """
        Populate param_index_map and identifying Active Lines.
        """
        self.param_index_map = []
        varying_indices = set()

        for name, param in lmfit_params.items():
            if param.vary and "_line_" in name:
                prefix, sep, suffix = name.rpartition('_line_')

                if prefix in self.dynamic_arrays:
                    try:
                        global_idx = int(suffix)
                        
                        local_idx = self.global_to_local_idx[global_idx]
                        
                        self.param_index_map.append((name, prefix, local_idx))
                        varying_indices.add(local_idx)
                    except ValueError:
                        pass
        
        # Sort and Store Active/Static Lists
        self.active_line_indices = np.sort(list(varying_indices)).astype(np.int32)
        all_indices = np.arange(self.nlines, dtype=np.int32)
        self.static_line_indices = np.setdiff1d(all_indices, self.active_line_indices)
        
        # Reset Cache
        self.cached_static_spectrum = {}
        self.cache_state = {}
    

    def update_from_lmfit(self, params):
        """
        The Fast Unpacker. Called every iteration of the fit.
        """
        for name, key, idx in self.param_index_map:
            if name in params:
                self.dynamic_arrays[key][idx] = params[name].value










 