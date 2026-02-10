from bisect import bisect
import re
import warnings

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from .hapi import ISO, PYTIPS2017, PYTIPS2011, PYTIPS2021, pcqsdhc, PROFILE_LORENTZ
from .utilities import molecularMass, etalon, convolveSpectrumSame
from .codata import CONSTANTS
from .o2_cia_karman import O2_CIA_Karman_Model
from numba import jit, prange

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
                 lineprofile = 'mHTP'):
        self.nlines = len(parameter_linelist)
        self.lineprofile = lineprofile

        #Static arrays
        self.elower = parameter_linelist['elower'].to_numpy(dtype=np.float64)
        self.molec_id = parameter_linelist['molec_id'].to_numpy(dtype=np.int32)   
        self.local_iso_id = parameter_linelist['local_iso_id'].to_numpy(dtype=np.int32)
        self.sw_scale_factor = parameter_linelist['sw_scale_factor'].to_numpy(dtype=np.float64)

        if self.lineprofile == 'HTP':
            # HTP: Real = nuVC, Imag = eta (Eta has no T-dep in HTP)
            self.col_map = {
                're': 'nuVC_', 'n_re': 'n_nuVC_',
                'im': 'eta_',  'n_im': 'n_eta_' # Map even if unused, will default to 0
            }
        else:
            # mHTP: Real = NuOptRe, Imag = NuOptIm (Both have T-dependence)
            self.col_map = {
                're': 'nuOptRe_', 'n_re': 'n_nuOptRe_',
                'im': 'nuOptIm_', 'n_im': 'n_nuOptIm_'
            }


        #Dynamic
        self.dynamic_arrays = {}
        for col in parameter_linelist.columns:
            # Check if this column is one of our special mapped ones
            is_mapped = False
            for key in ['re', 'n_re', 'im', 'n_im']:
                prefix = self.col_map[key]
                if col.startswith(prefix):
                    diluent = col[len(prefix):] 
                    internal_name = f"p_{key}_{diluent}" # e.g. p_re_air
                    self.dynamic_arrays[internal_name] = parameter_linelist[col].values.astype(np.float64)
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
                            T, p, molefraction, Diluent, spectrum_number, spectrum_min, 
                            baseline_coeffs = None, 
                            etalon_dict = None, 
                            cia_config = None, 
                            ILS_function = None, ILS_parameters = None, ILS_wing = None,
                            interpolated_compressability_file=None,
                            #pass through for LBL_alpha
                             TIPS = PYTIPS2021, isotope_list = ISO, 
                            natural_abundance = True, abundance_ratio_MI = {}, 
                            BIA_slope=False, BIA_FW_LBL=False,
                            IntensityThreshold = 1e-30, 
                            wing_cutoff = 25, wing_wavenumbers = 25, wing_method = 'wing_wavenumbers',):
        
        #Calculate LBL
        LBL_alpha = self.calculate_lbl_absorbance(waves, T, p, molefraction, Diluent, 
                                 spectrum_number, 
                                 interpolated_compressability_file = interpolated_compressability_file, TIPS = TIPS, isotope_list = isotope_list, 
                                 natural_abundance = natural_abundance, abundance_ratio_MI = abundance_ratio_MI, 
                                 BIA_slope=BIA_slope, BIA_FW_LBL=BIA_FW_LBL,
                                 IntensityThreshold = IntensityThreshold, 
                                 wing_cutoff = wing_cutoff, wing_wavenumbers = wing_wavenumbers, wing_method = wing_method,)
        LBL_alpha = LBL_alpha* 1e6
        
        # Calculate CIA
        alpha_cia = np.zeros_like(waves)
        if cia_config:
            if cia_config['model'] == 'Karman':
                    alpha_cia = cia_config['calculator'].calculate_cia(
                        waves, T, p, Diluent, 
                        **cia_config['params'])  
                    
            if cia_config['model'] == 'ad hoc':
                alpha_cia = cia_config['values']
                 
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

        total_alpha = baseline + etalon_signal + LBL_alpha + alpha_cia

        if ILS_function is not None and ILS_parameters is not None:
            # Note: convolveSpectrumSame returns (waves, alpha, ...)
            # We ignore the extra returns for the model result
            _, total_alpha, _, _, _ = convolveSpectrumSame(
                waves, total_alpha, 
                SlitFunction=ILS_function, 
                Resolution=ILS_parameters, 
                AF_wing=ILS_wing)
        
        return total_alpha

    def calculate_lbl_absorbance(self, waves, T, p, molefraction, Diluent, 
                                 spectrum_number, 
                                 interpolated_compressability_file = None, TIPS = PYTIPS2021, isotope_list = ISO, 
                                 natural_abundance = True, abundance_ratio_MI = {}, 
                                 BIA_slope=False, BIA_FW_LBL=False,
                                 IntensityThreshold = 1e-30, 
                                 wing_cutoff = 25, wing_wavenumbers = 25, wing_method = 'wing_wavenumbers',):


        #define reference temperature/pressure and calculate molecular density
        Tref = 296. # K
        pref = 1. # atm
        mol_density_ref = (pref/ CONSTANTS['cpa_atm'])/(CONSTANTS['k']*273.15) #density at 1 atm and 273.15 K
        mol_dens = (p/ CONSTANTS['cpa_atm'])/(CONSTANTS['k']*T)
        if interpolated_compressability_file != None:
            mol_dens = mol_dens / interpolated_compressability_file([p, T])[0]
        density_amagat = mol_dens / mol_density_ref

        #Vectorized
        sigma_T = np.ones(self.nlines, dtype=np.float64)
        sigma_Tref = np.ones(self.nlines, dtype=np.float64)
        #mass = np.ones(self.nlines, dtype=np.float64)
        abundance_ratio = np.ones(self.nlines, dtype=np.float64)
        for m, i in self.unique_pairs:
            m, i = int(m), int(i)
            mask = (self.molec_id == m) & (self.local_iso_id == i)
            try:
                sigma_T[mask] = TIPS(m, i, T)
                sigma_Tref[mask] = TIPS(m, i, Tref)
                if (not natural_abundance) and (abundance_ratio_MI != {}):
                    abundance_ratio[mask] = abundance_ratio_MI[m][i]
            except:
                pass


        #Calculate Line Intensity and Doppler Broadening
        
        sw_array = self._resolve_array('sw', spectrum_number)
        sw_array = sw_array*self.sw_scale_factor
        nu_array = self._resolve_array('nu', spectrum_number)
        GammaD = np.sqrt(2*CONSTANTS['k']*CONSTANTS['Na']*T*np.log(2)/(self.mass))*nu_array / CONSTANTS['c']

        line_intensity = sw_array * (sigma_Tref / sigma_T) * \
                    (np.exp(-CONSTANTS['c2'] * self.elower / T) * (1 - np.exp(-CONSTANTS['c2'] * nu_array / T))) / \
                    (np.exp(-CONSTANTS['c2'] * self.elower / Tref) * (1 - np.exp(-CONSTANTS['c2'] * nu_array / Tref)))
        
        # Calculated Line Parameters across Broadeners
        abundances_list = []
        
        # Accumulate arrays for Numba
        g0_list = []; n_g0_list = []
        d0_list = []; n_d0_list = []
        sd_gamma_list = []; n_gamma2_list = []
        sd_delta_list = []; n_delta2_list = []
        #nuvc_list = []; n_nuvc_list = []
        #eta_list = []
        param_re_list = []; n_param_re_list = []
        param_im_list = []; n_param_im_list = []
        y_list = []; n_y_list = []

        for species, info in Diluent.items():
            abundances_list.append(info['composition'])
            
            # Resolve all arrays (using _resolve_array which now defaults to zeros if missing)
            g0_list.append(self._resolve_array(f'gamma0_{species}', spectrum_number))
            n_g0_list.append(self._resolve_array(f'n_gamma0_{species}', spectrum_number))
            d0_list.append(self._resolve_array(f'delta0_{species}', spectrum_number))
            n_d0_list.append(self._resolve_array(f'n_delta0_{species}', spectrum_number))
            sd_gamma_list.append(self._resolve_array(f'SD_gamma_{species}', spectrum_number))
            n_gamma2_list.append(self._resolve_array(f'n_gamma2_{species}', spectrum_number))
            sd_delta_list.append(self._resolve_array(f'SD_delta_{species}', spectrum_number))
            n_delta2_list.append(self._resolve_array(f'n_delta2_{species}', spectrum_number))
            #nuvc_list.append(self._resolve_array(f'nuVC_{species}', spectrum_number))
            #n_nuvc_list.append(self._resolve_array(f'n_nuVC_{species}', spectrum_number))
            #eta_list.append(self._resolve_array(f'eta_{species}', spectrum_number))
            param_re_list.append(self._resolve_array(f'p_re_{species}', spectrum_number))
            n_param_re_list.append(self._resolve_array(f'p_n_re_{species}', spectrum_number))
            
            param_im_list.append(self._resolve_array(f'p_im_{species}', spectrum_number))
            n_param_im_list.append(self._resolve_array(f'p_n_im_{species}', spectrum_number))

            y_list.append(self._resolve_array(f'y_{species}', spectrum_number))
            n_y_list.append(self._resolve_array(f'n_y_{species}', spectrum_number))
            
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

        # BIA Logic (Keep in NumPy for now as it's conditional)
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
        
        #LBL Loop
        Xsect = np.zeros_like(waves)
        if wing_method == 'wing_wavenumbers':
            line_cutoffs = np.ones(self.nlines)*wing_wavenumbers
        elif wing_method == 'wing_cutoff':
            line_cutoffs = (0.5346 * Gamma0 + np.sqrt(0.2166 * Gamma0**2 + GammaD**2)) * wing_cutoff


        valid_indices = np.where(line_intensity >= IntensityThreshold)[0]

        for i in valid_indices:
            nu_i = nu_array[i]
            cut_i = line_cutoffs[i]

            
            BoundIndexLower = bisect(waves, nu_i - cut_i)
            BoundIndexUpper = bisect(waves, nu_i + cut_i)
            if BoundIndexLower >= len(waves) or BoundIndexUpper <= 0:
                continue
            
            wave_slice = waves[BoundIndexLower:BoundIndexUpper]

            if self.lineprofile == 'HTP':
                lineshape_PT, lineshape_vals_imag = pcqsdhc(
                    nu_i, GammaD[i], Gamma0[i], Gamma2[i], Delta0[i], Delta2[i],
                    ParamRe[i], ParamIm[i], wave_slice, Ylm=Y[i])  # This now includes the line-mixing component
            else:
                pass
            
            mf = molefraction[self.molec_id[i]]
            line_abundance_ratio = abundance_ratio[i]
            if BIA_slope:
                 Xsect[BoundIndexLower:BoundIndexUpper] += mol_dens  * \
                                                            mf * line_abundance_ratio * \
                                                            intensity_bia[i] * ( lineshape_PT)
                 if BIA_FW_LBL and bia_duration[i]!=0:
                    BoundIndexLower_BIA = bisect(waves, nu_i - 5*cut_i)
                    BoundIndexUpper_BIA = bisect(waves, nu_i + 5*cut_i)
                    wave_slice_BIA = waves[BoundIndexLower_BIA:BoundIndexUpper_BIA]
                    BIA_profile = PROFILE_LORENTZ(nu_i, 1/(2*np.pi*CONSTANTS['c']*1e-12*bia_duration[i]), 0, wave_slice_BIA)
                    Xsect[BoundIndexLower_BIA:BoundIndexUpper_BIA] += mol_dens  * \
                                                                mf * line_abundance_ratio * \
                                                                (line_intensity[i] - intensity_bia[i]) * BIA_profile
                     
            else:
                Xsect[BoundIndexLower:BoundIndexUpper] += mol_dens  * \
                                                            mf * line_abundance_ratio * \
                                                            line_intensity[i] * ( lineshape_PT)
        
        return (np.asarray(Xsect))
    
    def configure_for_fitting(self, lmfit_params):
        """
       Populate param_index_map from lmfit_params 
        """
        for name, param in lmfit_params.items():
            if param.vary and "_line_" in name:
                prefix, sep, suffix = name.rpartition('_line_')

                if prefix in self.dynamic_arrays:
                    try:
                        line_idx = int(suffix)
                        # Store tuple: (lmfit_name, array_key, array_index)
                        self.param_index_map.append((name, prefix, line_idx))
                    except ValueError:
                        pass # Suffix wasn't an integer, ignore
    
    def update_from_lmfit(self, params):
        """
        The Fast Unpacker. Called every iteration of the fit.
        """
        for name, key, idx in self.param_index_map:
            if name in params:
                self.dynamic_arrays[key][idx] = params[name].value










 