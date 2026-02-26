#Import Packages
# from .Utilities import *
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
from .spectroscopic_model import Spectroscopic_model
from .mHT.profile import beta

from lmfit import Minimizer,  Parameters, conf_interval, printfuncs

# lmfit generates warnings from the uncertainties module about params with zero uncertainty
# this is expected behavior and we should be able to safely ignore them
warnings.filterwarnings("ignore", category=UserWarning, module="uncertainties")

def convert_int_to_float(df, exclude_cols=None):
    mask = (df.drop(columns=exclude_cols, axis=1) if exclude_cols else df).select_dtypes(int)
    df[mask.columns] = mask.astype(float)
    return df


class Fit_DataSet:
    """Provides the fitting functionality for a Dataset."""

    def __init__(self, dataset, base_linelist_file, param_linelist_file, CIA_linelist_file = None,
                 lineprofile = 'mHTP', numba_lineprofile = True,
                 minimum_parameter_fit_intensity = 1e-30, minimum_simulation_intensity=1e-30,
                 weight_spectra = False,
                 
                 pressure_bounds = None,
                 temperature_bounds = None,
                 abundance_ratio_bounds = None,
                 etalon_amp_bounds = None,  
                 etalon_period_bounds = None, 
                 etalon_phase_bounds = None, 
                 x_shift_bounds = None,

                 nu_relative_bounds = None, 
                 sw_bound_scalars = None, 

                 gamma0_bounds = None, n_gamma0_bounds = None,
                 delta0_bounds = None, n_delta0_bounds = None,
                 SD_gamma_bounds  = None, n_gamma2_bounds  = None,
                 SD_delta_bounds  = None, n_delta2_bounds  = None,
                 nuVC_bounds  = None, n_nuVC_bounds = None,
                 eta_bounds  = None,
                 nuOptRe_bounds = None, n_nuOptRe_bounds = None, 
                 nuOptIm_bounds = None, n_nuOptIm_bounds = None,
                 linemixing_bounds = None, n_linemixing_bounds = None, 

                 BIA_intensity_depletion_bounds = None, 
                 BIA_collision_duration_bounds = None,

                 beta_formalism = False, additional_columns = []):
        
        self.dataset = dataset
        #baseline linelist ingest
        self.base_linelist_file = base_linelist_file
        base_int_cols = ['Spectrum Number', 'Segment Number']
        self.baseline_list = pd.read_csv(self.base_linelist_file + '.csv')
        self.baseline_list = convert_int_to_float(self.baseline_list, exclude_cols=base_int_cols)
        
        #parameter linelist ingest
        int_cols = additional_columns.copy()
        int_cols += ['molec_id', 'local_iso_id']
        self.param_linelist_file = param_linelist_file
        raw_df = pd.read_csv(param_linelist_file + ".csv", index_col = 0)
        if 'nu' in raw_df.columns:
            raw_df.sort_values('nu', inplace=True)
        #raw_df.reset_index(drop=True, inplace=True)
        #raw_df = raw_df.loc[:, ~raw_df.columns.str.contains('Unnamed:')]

        self.lineparam_list = convert_int_to_float(raw_df, exclude_cols=int_cols)
        self.lineprofile = lineprofile
        self.numba_lineprofile = numba_lineprofile
        self.engine = Spectroscopic_model(self.lineparam_list, lineprofile = self.lineprofile, numba_lineprofile = self.numba_lineprofile,
                                          isotope_list = self.dataset.isotope_list, beta_formalism=beta_formalism)
        
        #CIA linelist ingest
        self.CIA_linelist_file = CIA_linelist_file
        if self.CIA_linelist_file == None:
            self.CIAparam_list = None
        else:
            self.CIAparam_list = pd.read_csv(self.CIA_linelist_file + '.csv')
        
        
        self.minimum_parameter_fit_intensity = minimum_parameter_fit_intensity # Minimum fit intensity
        self.minimum_simulation_intensity = minimum_simulation_intensity # Minimum simulation intensity
        self.weight_spectra = weight_spectra #Spectrum Weighting boolean

        #Limits!
        def init_bounds(b): return b if b is not None else [-np.inf, np.inf]

        self.pressure_bounds = init_bounds(pressure_bounds)
        self.temperature_bounds = init_bounds(temperature_bounds)
        self.x_shift_bounds = init_bounds(x_shift_bounds)
        self.abundance_ratio_bounds = init_bounds(abundance_ratio_bounds)

        self.etalon_amp_bounds = init_bounds(etalon_amp_bounds)
        self.etalon_period_bounds = init_bounds(etalon_period_bounds)
        self.etalon_phase_bounds = init_bounds(etalon_phase_bounds)

        self.nu_relative_bounds = init_bounds(nu_relative_bounds) 
        self.sw_bound_scalars = init_bounds(sw_bound_scalars)  

        self.gamma0_bounds = init_bounds(gamma0_bounds)
        self.n_gamma0_bounds = init_bounds(n_gamma0_bounds)
        self.delta0_bounds = init_bounds(delta0_bounds)
        self.n_delta0_bounds = init_bounds(n_delta0_bounds)
        self.SD_gamma_bounds = init_bounds(SD_gamma_bounds)
        self.n_gamma2_bounds = init_bounds(n_gamma2_bounds)
        self.SD_delta_bounds = init_bounds(SD_delta_bounds)
        self.n_delta2_bounds = init_bounds(n_delta2_bounds)
        
        if lineprofile == 'HTP':
            self.paramRe_bounds = init_bounds(nuVC_bounds)
            self.n_paramRe_bounds = init_bounds(n_nuVC_bounds)
            self.paramIm_bounds = init_bounds(eta_bounds)
            self.n_paramIm_bounds = [-np.inf, np.inf]
        else:
            self.paramRe_bounds = init_bounds(nuOptRe_bounds)
            self.n_paramRe_bounds = init_bounds(n_nuOptRe_bounds)
            self.paramIm_bounds = init_bounds(nuOptIm_bounds)
            self.n_paramIm_bounds = init_bounds(n_nuOptIm_bounds)
        
        self.linemixing_bounds = init_bounds(linemixing_bounds)
        self.n_linemixing_bounds = init_bounds(n_linemixing_bounds)
        self.BIA_intensity_depletion_bounds = init_bounds(BIA_intensity_depletion_bounds)
        self.BIA_collision_duration_bounds = init_bounds(BIA_collision_duration_bounds)

        self.beta_formalism = beta_formalism #Beta formalism

        self.spec_attrs = self.prep_sim()
        
        # Optimization Cache
        self.fit_data_cache = []

    def _extract_baseline_coeffs(self, params, spec_num, seg_num):
        """Reconstructs polynomial coeffs list [c_n, ..., c_0] from params."""
        baseline_coeffs = []
        for i in range(self.dataset.baseline_order + 1):
            char = chr(97 + i) # 0->'a', 1->'b'
            p_name = f"baseline_{char}_{spec_num}_{seg_num}"
            if p_name in params:
                baseline_coeffs.append(params[p_name].value)
            else:
                baseline_coeffs.append(0.0)
        
        return baseline_coeffs[::-1]
    
    def _extract_etalon_dict(self, params, spec_num, seg_num):
        """Reconstructs nested etalon dictionary {1: {'amp':...}, ...}"""
        etalon_dict = {}
        for i in range(1, 11): 
            amp_name = f"etalon_{i}_amp_{spec_num}_{seg_num}"
            if amp_name in params:
                period_name = f"etalon_{i}_period_{spec_num}_{seg_num}"
                phase_name = f"etalon_{i}_phase_{spec_num}_{seg_num}"
                
                etalon_dict[i] = {
                    'amp': params[amp_name].value,
                    'period': params[period_name].value if period_name in params else 1.0,
                    'phase': params[phase_name].value if phase_name in params else 0.0
                }
        return etalon_dict if etalon_dict else None
    
    def _extract_ils_resolution(self, params, spectrum, segment):
        """Extracts ILS resolution parameters list."""
        if spectrum.ILS_function is None: return None
            
        func_name = spectrum.ILS_function.__name__
        num_params = self.dataset.ILS_function_dict.get(func_name, 0)
        
        if num_params == 0: return []
            
        resolution_params = []
        for i in range(num_params):
            p_name = f"{func_name}_res_{i}_{spectrum.spectrum_number}_{segment}"
            if p_name in params:
                resolution_params.append(params[p_name].value)
            else:
                resolution_params.append(0.0) 
        return resolution_params

    def _extract_cia_config(self, params):
        """Packs CIA model and parameters into a config dict."""
        if self.dataset.CIA_model['model'] == 'Karman':
            cia_params = {}
            param_map = {
                'S_SO_O2_O2': 'SO_O2', 'S_SO_O2_N2': 'SO_N2', 'S_EXCH_O2_O2': 'EXCH_O2',
                'EXCH_b_O2_O2': 'EXCH_b', 'EXCH_c_O2_O2': 'EXCH_c',
                'SO_b_O2_O2': 'SO_b_O2_O2', 'SO_c_O2_O2': 'SO_c_O2_O2',
                'SO_b_O2_N2': 'SO_b_O2_N2', 'SO_c_O2_N2': 'SO_c_O2_N2',
                'SO_shift_O2_O2': 'SO_shift_O2_O2', 'SO_shift_O2_N2': 'SO_shift_O2_N2',
                'EXCH_shift_O2_N2': 'EXCH_shift' 
            }
            
            for fit_name, arg_name in param_map.items():
                if fit_name in params:
                    cia_params[arg_name] = params[fit_name].value
                else:
                    cia_params[arg_name] = 0.0
                    
            return {
                'model': 'Karman',
                'calculator': self.spec_attrs['Dataset']['CIA model'],
                'params': cia_params
            }
        else:
            return None

    def prep_sim(self):
        spectrum_attributes = {"Compressability Factor": None, 'CIA model': None} 
        spectra_numbers = self.dataset.get_list_spectrum_numbers()
        spectra_numbers.append('Dataset')

        spectra_attribute_dict = {spec_num: spectrum_attributes.copy() for spec_num in spectra_numbers}        

        for spectrum in self.dataset.spectra:
            if spectrum.compressability_file != None:
                comp_factor = pd.read_csv(spectrum.compressability_file + '.csv')
                pressures = np.asarray(comp_factor['Pressure (MPa)'].values*1e6/101325).astype(float)
                temperatures = list(comp_factor)
                temperatures.remove('Pressure (MPa)')
                temperatures = np.asarray(temperatures).astype(float)
                comp_factor.drop('Pressure (MPa)', inplace=True, axis=1) 
                comp_factor_array = comp_factor.to_numpy()
                spectra_attribute_dict[spectrum.spectrum_number]['Compressability Factor'] = RegularGridInterpolator(points = [pressures, temperatures], values = comp_factor_array)
            else:
                spectra_attribute_dict[spectrum.spectrum_number]['Compressability Factor'] = None
        
        if self.dataset.CIA_model['model'] == 'Karman':
            spectra_attribute_dict['Dataset']['CIA model'] =O2_CIA_Karman_Model(self.dataset.CIA_model['band'])

        return spectra_attribute_dict
 
    def generate_params(self):
        """Generates the lmfit parameter object that will be used in fitting."""
       
        params = Parameters()
        baseline_parameters = []
        for base_param in list(self.baseline_list):
            if ('_vary' not in base_param) and ('_err' not in base_param) and ('Spectrum Number' not in base_param) and ('Segment Number' not in base_param):
                baseline_parameters.append(base_param)
        for index in self.baseline_list.index.values:
            spec_num = int(self.baseline_list.iloc[index]['Spectrum Number'])
            seg_num = int(self.baseline_list.iloc[index]['Segment Number'])
            
            for base_param in baseline_parameters:
                val = self.baseline_list.loc[index][base_param]
                vary = self.baseline_list.loc[index][base_param + '_vary']
                lmfit_name = f"{base_param}_{spec_num}_{seg_num}"
                if ('Pressure' in base_param) and self.pressure_bounds != [-np.inf, np.inf]:
                    params.add(lmfit_name, val, vary, min = self.pressure_bounds[0], max = self.pressure_bounds[1])
                elif ('Temperature' in base_param) and self.temperature_bounds != [-np.inf, np.inf]:
                    params.add(lmfit_name, val, vary, min = self.temperature_bounds[0], max = self.temperature_bounds[1])
                elif ('etalon_' in base_param) and self.etalon_amp_bounds and ('amp' in base_param):
                    params.add(lmfit_name, val, vary, min = self.etalon_amp_bounds[0], max = self.etalon_amp_bounds[1])
                elif ('etalon_' in base_param) and self.etalon_amp_bounds and ('period' in base_param):
                    params.add(lmfit_name, val, vary, min = self.etalon_period_bounds[0], max = self.etalon_period_bounds[1])
                elif ('etalon_' in base_param) and self.etalon_amp_bounds and ('phase' in base_param):
                    params.add(lmfit_name, val, vary, min = self.etalon_phase_bounds[0], max = self.etalon_phase_bounds[1])
                elif ('x_shift' in base_param) and self.x_shift_bounds != [-np.inf, np.inf]:
                    params.add(lmfit_name, val, vary, min = self.x_shift_bounds[0], max = self.x_shift_bounds[1])
                elif ('abundance_ratio_' in base_param) and self.abundance_ratio_bounds != [-np.inf, np.inf]:
                     params.add(lmfit_name, val, vary, min = self.abundance_ratio_bounds[0], max = self.abundance_ratio_bounds[1])
                else:
                    params.add(lmfit_name, val, vary)
        
        if self.lineprofile == 'HTP':
            param_Re = 'nuVC'; param_Im = 'eta'            
        else:
            param_Re = 'nuOptRe'; param_Im = 'nuOptIm'
        
        linelist_params = []
        num_nominal_temps = self.dataset.get_number_nominal_temperatures()[0]
        valid_column_prefix = ['nu', 'sw', 'gamma0', 'delta0', 'SD_gamma', 'SD_delta', param_Re, param_Im, 'y']
        if num_nominal_temps > 1: valid_column_prefix.append('n_')
        if self.dataset.BIA_model['sw_depletion']:
            valid_column_prefix.append('BIA_slope')
            if self.dataset.BIA_model['farwing_continuum'] == 'LBL':
                valid_column_prefix.append('BIA_collision_duration')

        column_suffix_to_ignore = ['_err', '_vary']
        excluded_cols = ['molec_id', 'local_iso_id', 'elower', 'sw_scale_factor']

        candidate_params = [
            col for col in self.lineparam_list.columns
            if any(col.startswith(p) for p in valid_column_prefix) and not any(col.endswith(s) for s in column_suffix_to_ignore) and col not in excluded_cols
        ]
        self.constrain_dictionary = {'nu':True,'sw':True, 'gamma0': True, 'delta0': True, 'SD_gamma': True, 'SD_delta': True, param_Re: True, param_Im: True, 'y': True }
        sorted_constrain_keys = sorted(self.constrain_dictionary.keys(), key=len, reverse=True)
        for col in candidate_params:
            spectrum_specific = False
            for other_col in candidate_params:
                if other_col.startswith(col + '_') and other_col != col:
                    suffix = other_col[len(col)+1:]
                    if suffix.isdigit():
                        spectrum_specific = True
                        break
            if spectrum_specific:
                for key in sorted_constrain_keys:
                    if col.startswith(key):
                        self.constrain_dictionary[key] = False
                        break
            else:
                linelist_params.append(col)
        

        for spec_line in self.lineparam_list.index.values:
            sw_scaled = self.lineparam_list.loc[spec_line]['sw'] * self.lineparam_list.loc[spec_line]['sw_scale_factor']
            if sw_scaled < self.minimum_parameter_fit_intensity: continue
            for line_param in linelist_params:
                val = self.lineparam_list.loc[spec_line][line_param]
                vary = self.lineparam_list.loc[spec_line][line_param + '_vary']
                lmfit_name = f"{line_param}_line_{spec_line}"

                if (line_param == 'nu') or (('nu_' in line_param) and (param_Re not in line_param) and (param_Im not in line_param)):
                    params.add(lmfit_name, val, vary, min = val + self.nu_relative_bounds[0], max = val + self.nu_relative_bounds[1])
                elif ('sw' in line_param):
                    params.add(lmfit_name, val, vary, min = val * self.sw_bound_scalars[0], max = val * self.sw_bound_scalars[1])
                elif ('gamma0_' in line_param) and ('n_' not in line_param): 
                    params.add(lmfit_name, val, vary, min =self.gamma0_bounds[0], max = self.gamma0_bounds[1])
                elif ('n_gamma0_' in line_param):
                    params.add(lmfit_name, val, vary, min = self.n_gamma0_bounds[0], max = self.n_gamma0_bounds[1])
                elif ('delta0_' in line_param) and ('n_' not in line_param):
                    params.add(lmfit_name, val, vary, min = self.delta0_bounds[0], max = self.delta0_bounds[1])
                elif ('n_delta0' in line_param):
                    params.add(lmfit_name, val, vary, min =self.n_delta0_bounds[0], max = self.n_delta0_bounds[1])
                elif ('SD_gamma' in line_param): 
                    params.add(lmfit_name, val, vary, min = self.SD_gamma_bounds[0], max = self.SD_gamma_bounds[1])
                elif ('n_gamma2' in line_param):
                    params.add(lmfit_name, val, vary, min = self.n_gamma2_bounds[0], max = self.n_gamma2_bounds[1])
                elif ('SD_delta' in line_param):
                    params.add(lmfit_name, val, vary, min = self.SD_delta_bounds[0], max = self.SD_delta_bounds[1])
                elif ('n_delta2' in line_param):
                    params.add(lmfit_name, val, vary, min = self.SD_delta_bounds[0], max = self.SD_delta_bounds[1])
                elif (param_Re in line_param) :
                    params.add(lmfit_name, val, vary, min = self.paramRe_bounds[0], max = self.paramRe_bounds[1])
                elif ('n_' + param_Re in line_param):
                    params.add(lmfit_name, val, vary, min = self.n_paramRe_bounds[0], max = self.n_paramRe_bounds[1])
                elif (param_Im in line_param) and ('n_'  + param_Im not in line_param):
                    params.add(lmfit_name, val, vary, min = self.paramIm_bounds[0], max = self.paramIm_bounds[1])
                elif ('n_' + param_Im in line_param):
                    params.add(lmfit_name, val, vary, min = self.n_paramIm_bounds[0], max = self.n_paramIm_bounds[1])
                elif ('y_' in line_param) and ('n_' not in line_param):
                    params.add(lmfit_name, val, vary, min = self.linemixing_bounds[0], max = self.linemixing_bounds[1])
                elif ('n_y_' in line_param):
                    params.add(lmfit_name, val, vary, min = self.n_linemixing_bounds[0], max = self.n_linemixing_bounds[1])
                elif ('BIA_slope_' in line_param) and (self.constrain_dictionary['sw']):
                    params.add(lmfit_name, val, vary, min = self.BIA_intensity_depletion_bounds[0], max = self.BIA_intensity_depletion_bounds[1])
                elif ('BIA_collision_duration_' in line_param) and (self.constrain_dictionary['sw']):
                    params.add(lmfit_name, val, vary, min = self.BIA_collision_duration_bounds[0], max = self.BIA_collision_duration_bounds[1])
                else:
                    print ('yup' + line_param)

        if self.dataset.CIA_model['model'] == "Karman":
            cia_parameters = []
            for cia_param in list(self.CIAparam_list):
                if ('_vary' not in cia_param) and ('_err' not in cia_param) and ('CIA Pair' not in cia_param):
                    cia_parameters.append(cia_param)
            for cia_pair in self.CIAparam_list['CIA Pair']:
                index = self.CIAparam_list[self.CIAparam_list['CIA Pair']==cia_pair].index.values[0]   
                for cia_param in cia_parameters:
                    if 'SO' in cia_param:
                        params.add(cia_param + '_'+ cia_pair, self.CIAparam_list.loc[index][cia_param], self.CIAparam_list.loc[index][cia_param + '_vary'])
                    elif 'EXCH' in cia_param:
                        if cia_pair == 'O2_O2':
                            params.add(cia_param + '_'+ cia_pair, self.CIAparam_list.loc[index][cia_param], self.CIAparam_list.loc[index][cia_param + '_vary'])
                        elif cia_pair =='O2_N2':
                            params.add(cia_param + '_'+ cia_pair, 0, False)

        self.engine.configure_for_fitting(params)
        return (params)

    def constrained_baseline(self, params, baseline_segment_constrained = True, xshift_segment_constrained = True, molefraction_segment_constrained = True,
                                    etalon_amp_segment_constrained = True, etalon_period_segment_constrained = True, etalon_phase_segment_constrained = True,
                                    pressure_segment_constrained = True, temperature_segment_constrained = True, abundance_ratio_segment_constrained = True):

        spectrum_segment_min = {}
        for spectrum in self.dataset.spectra:
            spectrum_segment_min[spectrum.spectrum_number] = np.min(list(set(spectrum.segments)))

        for param in params:
            if ('Pressure' in param) and pressure_segment_constrained:
                indices = [m.start() for m in re.finditer('_', param)]
                spectrum_num = int(param[indices[0]+1:indices[1]])
                segment_num = int(param[indices[1]+1:])
            elif ('Temperature' in param) and temperature_segment_constrained:
                indices = [m.start() for m in re.finditer('_', param)]
                spectrum_num = int(param[indices[0]+1:indices[1]])
                segment_num = int(param[indices[1]+1:])

            elif ('baseline' in param) and baseline_segment_constrained:
                indices = [m.start() for m in re.finditer('_', param)]
                spectrum_num = int(param[indices[1]+1:indices[2]])
                segment_num = int(param[indices[2]+1:])
                if segment_num != spectrum_segment_min[spectrum_num]:
                    params[param].set(expr = param[:indices[1]+1] + str(spectrum_num) + '_' + str(spectrum_segment_min[spectrum_num]))
            elif ('x_shift' in param) and xshift_segment_constrained:
                indices = [m.start() for m in re.finditer('_', param)]
                spectrum_num = int(param[indices[1]+1:indices[2]])
                segment_num = int(param[indices[2]+1:])
                if segment_num != spectrum_segment_min[spectrum_num]:
                    params[param].set(expr = param[:indices[1]+1] + str(spectrum_num) + '_' + str(spectrum_segment_min[spectrum_num]))
            elif ('molefraction' in param) and molefraction_segment_constrained:
                indices = [m.start() for m in re.finditer('_', param)]
                spectrum_num = int(param[indices[1]+1:indices[2]])
                segment_num = int(param[indices[2]+1:])
                if segment_num != spectrum_segment_min[spectrum_num]:
                    params[param].set(expr = param[:indices[1]+1] + str(spectrum_num) + '_' + str(spectrum_segment_min[spectrum_num]))
            elif ('abundance_ratio' in param) and abundance_ratio_segment_constrained: #
                parts = param.split('_')
                spectrum_num = int(parts[-2])
                segment_num = int(parts[-1])
                if segment_num != spectrum_segment_min[spectrum_num]:
                    base_name = '_'.join(parts[:-2])
                    params[param].set(expr = f"{base_name}_{spectrum_num}_{spectrum_segment_min[spectrum_num]}")
                    
                
            elif ('etalon' in param):
                indices = [m.start() for m in re.finditer('_', param)]
                spectrum_num = int(param[indices[2]+1:indices[3]])
                segment_num = int(param[indices[3]+1:])
                if 'amp' in param and etalon_amp_segment_constrained:
                    if segment_num != spectrum_segment_min[spectrum_num]:
                        params[param].set(expr = param[:indices[2]+1] + str(spectrum_num) + '_' + str(spectrum_segment_min[spectrum_num]))
                elif 'period' in param and etalon_period_segment_constrained:
                    if segment_num != spectrum_segment_min[spectrum_num]:
                        params[param].set(expr = param[:indices[2]+1] + str(spectrum_num) + '_' + str(spectrum_segment_min[spectrum_num]))
                elif 'phase' in param and etalon_phase_segment_constrained:
                    if segment_num != spectrum_segment_min[spectrum_num]:
                        params[param].set(expr = param[:indices[2]+1] + str(spectrum_num) + '_' + str(spectrum_segment_min[spectrum_num]))
        return params
    
    def constrained_CIA(self, params, S_temperature_dependence_constrained = True, shift_constrained = True):
        if self.dataset.CIA_model['model'] == 'Karman':
            for param in params:
                if (('SO_b' in param) or ('SO_c' in param)) and S_temperature_dependence_constrained:
                    if 'O2_O2' not in param:
                        params[param].set(expr = param[:5] + 'O2_O2')
                if ('SO_shift' in param) and shift_constrained:
                    if 'O2_O2' not in param:
                        params[param].set(expr = param[:9] + 'O2_O2')
        return params

    def _precompute_fitting_data(self):
        """Pre-computes and caches data segments and weights to avoid inner-loop overhead."""
        self.fit_data_cache = []
        for spectrum in self.dataset.spectra:
            wavenumber_segments, alpha_segments, indices_segments = spectrum.segment_wave_alpha()
            for segment in set(spectrum.segments):
                # Pre-calculate what we can (base waves, data, indices)
                idx_min = np.min(indices_segments[segment])
                idx_max = np.max(indices_segments[segment])
                
                # Pre-calc Weighting Array
                w_array = None
                if self.weight_spectra:
                    if spectrum.tau_stats.all() == 0:
                        w_array = spectrum.weight
                    else:
                        seg_stats = spectrum.tau_stats[idx_min : idx_max + 1]
                        w_array = spectrum.weight * (1.0 / seg_stats)
                
                # Pre-calc CIA Slice indices if needed
                cia_slice_indices = None
                if self.dataset.CIA_model['model'] == 'ad hoc':
                    cia_slice_indices = (idx_min, idx_max + 1)

                self.fit_data_cache.append({
                    'spectrum': spectrum,
                    'segment': segment,
                    'base_waves': wavenumber_segments[segment],
                    'data_y': alpha_segments[segment],
                    'weights': w_array,
                    'cia_slice_indices': cia_slice_indices
                })

    def objective_function(self, params, wing_cutoff = 25, wing_wavenumbers=25, wing_method='wing_wavenumbers'):
        """
        Calculates residuals for the Minimizer using pre-computed data structures.
        """
        self.engine.update_from_lmfit(params)
        total_residuals = []
        
        cia_config = None
        if self.dataset.CIA_model['model'] == 'Karman':
            cia_config = self._extract_cia_config(params)

        # Iterate over pre-computed segments instead of nested object loops
        for item in self.fit_data_cache:
            spectrum = item['spectrum']
            segment = item['segment']
            
            # Apply X-Shift to base waves (this floats, so we calculate here)
            x_shift_name = f'x_shift_{spectrum.spectrum_number}_{segment}'
            x_shift = params[x_shift_name].value if x_shift_name in params else 0.0
            waves = item['base_waves'] + x_shift

            # Extract Parameters
            T = params[f'Temperature_{spectrum.spectrum_number}_{segment}'].value
            p = params[f'Pressure_{spectrum.spectrum_number}_{segment}'].value

            # Extract Molefraction
            mf = spectrum.molefraction.copy()
            for molec_id in mf:
                iso_name = self.dataset.isotope_list.get((molec_id, 1), [0,0,0,0,'Unknown'])[4]
                mf_param_name = f'molefraction_{iso_name}_{spectrum.spectrum_number}_{segment}'
                if mf_param_name in params:
                    mf[molec_id] = params[mf_param_name].value
            
            #Extract abundance ratios
            abundance_ratios = {}
            for m, iso_dict in spectrum.abundance_ratio_MI.items():
                for i, val in iso_dict.items():
                    abundance_ratios[(int(m), int(i))] = val

            unique_pairs = np.unique(np.column_stack((self.lineparam_list['molec_id'], self.lineparam_list['local_iso_id'])), axis=0)
            for m, i in unique_pairs:
                m, i = int(m), int(i)
                mol_name = self.dataset.isotope_list.get((m, 1), [0,0,0,0,'Unknown'])[4]
                abund_param_name = f'abundance_ratio_{mol_name}_{i}_{spectrum.spectrum_number}_{segment}'
                if abund_param_name in params:
                    abundance_ratios[(m, i)] = params[abund_param_name].value


            baseline_coeffs = self._extract_baseline_coeffs(params, spectrum.spectrum_number, segment)
            etalon_dict = self._extract_etalon_dict(params, spectrum.spectrum_number, segment)
            ils_res = self._extract_ils_resolution(params, spectrum, segment)

            # Handle Ad Hoc CIA using pre-calc indices
            if item['cia_slice_indices']:
                start, end = item['cia_slice_indices']
                cia_config = {'model': 'ad hoc', 'values': spectrum.cia[start:end]}

            model_y = self.engine.calculate_spectrum(
                waves=waves,
                T=T, p=p, molefraction=mf, abundance_ratios=abundance_ratios,
                Diluent=spectrum.Diluent,
                spectrum_number=spectrum.spectrum_number,
                spectrum_min=np.min(spectrum.wavenumber),
                segment = segment,
                baseline_coeffs=baseline_coeffs,
                etalon_dict=etalon_dict,
                cia_config=cia_config,
                ILS_function=spectrum.ILS_function,
                ILS_parameters=ils_res,
                ILS_wing=spectrum.ILS_wing,
                interpolated_compressability_file=self.spec_attrs[spectrum.spectrum_number]['Compressability Factor'],
                BIA_slope=self.dataset.BIA_model['sw_depletion'],
                BIA_FW_LBL=(self.dataset.BIA_model['farwing_continuum'] == 'LBL'),
                IntensityThreshold=self.minimum_simulation_intensity,
                wing_cutoff=wing_cutoff,
                wing_wavenumbers=wing_wavenumbers,
                wing_method=wing_method
            )

            resid = item['data_y'] - model_y
            
            # Apply pre-calculated weights
            if self.weight_spectra and item['weights'] is not None:
                resid = resid * item['weights']

            total_residuals.append(resid)
        
        return np.concatenate(total_residuals)

    def fit_data(self, params, wing_cutoff = 25, wing_wavenumbers = 25, wing_method = 'wing_wavenumbers', xtol = 1e-10, maxfev = 10000, ftol = 1e-10, 
                 method = 'least_squares'):
        """Uses the lmfit minimizer to do the fitting through the simulation model function."""

        fcn_args = (wing_cutoff, wing_wavenumbers, wing_method)

        # 1. Initialize the cache/config for the spectroscopic engine (Optimizes calculate_spectrum)
        self.engine.configure_for_fitting(params)
        
        # 2. Precompute data segments to avoid inner-loop overhead (Optimizes fit_dataset loop)
        self._precompute_fitting_data()

        floated_parameters = any(p.vary for p in params.values())
        if not floated_parameters and method == 'least_squares':
            method = 'leastsq'

        kws = {
            'xtol': xtol,
            'ftol': ftol,
        }

        if method == 'least_squares':
            kws['gtol'] = ftol
            kws['x_scale'] = 'jac'

        self.minner = Minimizer(self.objective_function, params, max_nfev=maxfev,
                           fcn_args = fcn_args, **kws)
        
        
        result = self.minner.minimize(method = method)
        return result
        

    def residual_analysis(self, result, indv_resid_plot = False):
        """Updates the model and residual arrays in each spectrum object."""
        residual_array = result.residual

        for spectrum in self.dataset.spectra:
            spectrum_residual, residual_array = np.split(residual_array, [len(spectrum.wavenumber)])

            if self.weight_spectra:
                if spectrum.tau_stats.all() == 0:
                    weights = len(spectrum_residual)*[spectrum.weight]
                else:
                    pt_by_pt_weights= 1 / (spectrum.tau_stats)
                    weights = spectrum.weight * pt_by_pt_weights
                if spectrum.weight == 0:
                    spectrum_residual = np.asarray(len(spectrum_residual)*[0])
                else:
                    spectrum_residual  = spectrum_residual / weights

            spectrum.set_residuals(spectrum_residual)
            spectrum.set_model(spectrum.alpha - spectrum_residual)
            if indv_resid_plot:
                spectrum.plot_model_residuals()       


    def update_params(self, result, base_linelist_update_file = None , param_linelist_update_file = None, 
                      CIA_linelist_update_file = None):
        """Updates the baseline and line parameter files based on fit results."""
        # [This function body remains the same as in your input, ensuring correct indentation]
        if base_linelist_update_file == None:
            base_linelist_update_file = self.base_linelist_file
        if param_linelist_update_file == None:
            param_linelist_update_file = self.param_linelist_file
        if CIA_linelist_update_file == None:
            CIA_linelist_update_file = self.CIA_linelist_file


        for key, par in result.params.items():
            #Baseline
            if ('Pressure' in par.name) or ('Temperature' in par.name):
                indices = [m.start() for m in re.finditer('_', par.name)]
                parameter = (par.name[:indices[0]])
                spectrum = int(par.name[indices[0] + 1:indices[1]])
                segment = int(par.name[indices[1] + 1:])
                mask = (self.baseline_list['Segment Number'] == segment) & (self.baseline_list['Spectrum Number'] == spectrum)
                self.baseline_list.loc[mask, parameter] = par.value
                if par.vary:
                    self.baseline_list.loc[mask, parameter + '_err'] = par.stderr
            elif ('molefraction' in par.name) or ('baseline' in par.name) or ('x_shift' in par.name) or ('abundance_ratio' in par.name):
                parts = par.name.split('_')
                segment = int(parts[-1])
                spectrum = int(parts[-2])

                parameter = '_'.join(parts[:-2])
                mask = (self.baseline_list['Segment Number'] == segment) & (self.baseline_list['Spectrum Number'] == spectrum)
                self.baseline_list.loc[mask, parameter] = par.value
                if par.vary:
                    self.baseline_list.loc[mask, parameter + '_err'] = par.stderr


                '''
                indices = [m.start() for m in re.finditer('_', par.name)]
                parameter = (par.name[:indices[1]])
                spectrum = int(par.name[indices[1] + 1:indices[2]])
                segment = int(par.name[indices[2] + 1:])
                mask = (self.baseline_list['Segment Number'] == segment) & (self.baseline_list['Spectrum Number'] == spectrum)
                self.baseline_list.loc[mask, parameter] = par.value
                if par.vary:
                    self.baseline_list.loc[mask, parameter + '_err'] = par.stderr
                '''

            elif ('etalon' in par.name):
                indices = [m.start() for m in re.finditer('_', par.name)]
                parameter = par.name[:indices[2]]
                spectrum = int(par.name[indices[2]+1:indices[3]])
                segment = int(par.name[indices[3]+1:])
                mask = (self.baseline_list['Segment Number'] == segment) & (self.baseline_list['Spectrum Number'] == spectrum)
                self.baseline_list.loc[mask, parameter] = par.value
                if par.vary:
                    self.baseline_list.loc[mask, parameter + '_err'] = par.stderr

            elif ('_res_' in par.name):
                indices = [m.start() for m in re.finditer('_', par.name[par.name.find('_res_') + 5:])]
                spectrum = par.name[par.name.find('_res_') + 5:][indices[0]+1:indices[1]]
                segment = par.name[par.name.find('_res_') + 5:][indices[1]+1:]
                parameter = par.name[:par.name.find('_res_')] + '_res_' + par.name[par.name.find('_res_') + 5:][:indices[0]]
                mask = (self.baseline_list['Segment Number'] == segment) & (self.baseline_list['Spectrum Number'] == spectrum)
                self.baseline_list.loc[mask, parameter] = par.value
                if par.vary:
                    self.baseline_list.loc[mask, parameter + '_err'] = par.stderr
            #CIA
            elif (self.dataset.CIA_model['model']=='Karman') and ('O2_O2' in par.name):
                indices = [m.start() for m in re.finditer('_', par.name)]
                parameter = par.name[:indices[1]]
                mask = (self.CIAparam_list['CIA Pair'] == 'O2_O2')
                self.CIAparam_list.loc[mask, parameter] = par.value

                if par.vary:
                    self.CIAparam_list.loc[mask, parameter + '_err'] = par.stderr

            elif (self.dataset.CIA_model['model']=='Karman') and ('O2_N2' in par.name):
                indices = [m.start() for m in re.finditer('_', par.name)]
                parameter = par.name[:indices[1]]
                mask = (self.CIAparam_list['CIA Pair'] == 'O2_N2')
                self.CIAparam_list.loc[mask, parameter] = par.value
                if par.vary:
                    self.CIAparam_list.loc[mask, parameter + '_err'] = par.stderr


            #Line shape Parameters
            else:
                parameter = par.name[:par.name.find('_line')]
                line = int(par.name[par.name.find('_line')+6:])
                self.lineparam_list.loc[line, parameter] = par.value
                if par.vary:
                    self.lineparam_list.loc[line, parameter + '_err'] = par.stderr



        self.baseline_list.to_csv(base_linelist_update_file + '.csv', index = False)
        self.lineparam_list.to_csv(param_linelist_update_file + '.csv')
        if self.dataset.CIA_model['model'] == "Karman":
            self.CIAparam_list.to_csv(CIA_linelist_update_file + '.csv', index = False)

        #Calculate Baseline + Etalons and add to the Baseline term for each spectra
        for spectrum in self.dataset.spectra:
            wavenumber_segments, alpha_segments, indices_segments = spectrum.segment_wave_alpha()
            baseline = len(spectrum.wavenumber)*[0]
            for segment in set(spectrum.segments):
                waves = wavenumber_segments[segment]
                bound_min = np.min(indices_segments[segment])
                bound_max = np.max(indices_segments[segment])
                wave_rel = waves - np.min(spectrum.wavenumber)
                baseline_param_array = [0]*(self.dataset.baseline_order+1)
                fit_etalon_parameters = {}
                for i in range(1, len(spectrum.etalons)+1):
                    fit_etalon_parameters[i] = {'amp': 0, 'period':0, 'phase':0}
                for key, par in result.params.items():
                    if ('baseline' in par.name) and ((str(spectrum.spectrum_number) + '_' + str(segment)) in par.name):
                        baseline_param_array[ord(par.name[9:par.name.find('_',9)])-97] = float(par.value)
                    elif ('etalon' in par.name) and ((str(spectrum.spectrum_number) + '_' + str(segment)) in par.name[par.name.find('_', 7):]):
                        etalon_num = int(par.name[par.name.find('_')+1: par.name.find('_', par.name.find('_')+1)])
                        if ('amp' in par.name) and (str(etalon_num) in par.name):
                            fit_etalon_parameters[etalon_num]['amp'] = par.value
                        if ('period' in par.name) and (str(etalon_num) in par.name):
                            fit_etalon_parameters[etalon_num]['period'] = par.value
                        if ('phase' in par.name) and (str(etalon_num) in par.name):
                            fit_etalon_parameters[etalon_num]['phase'] = par.value
                baseline_param_array = baseline_param_array[::-1] 
                baseline[bound_min: bound_max + 1] += np.polyval(baseline_param_array, wave_rel)
                for i in range(1, len(spectrum.etalons)+1):
                    baseline[bound_min: bound_max +1] += etalon(wave_rel, fit_etalon_parameters[i]['amp'], fit_etalon_parameters[i]['period'], fit_etalon_parameters[i]['phase'])
            spectrum.set_background(baseline)
        
        #Calculate CIA
        if self.dataset.CIA_model['model']!= None:
            if self.dataset.CIA_model['model']=='Karman':
                for spectrum in self.dataset.spectra:
                    CIA = self.spec_attrs['Dataset']['CIA model'].calculate_cia(spectrum.wavenumber, spectrum.get_temperature(), spectrum.get_pressure(), spectrum.get_Diluent(),
                        float(result.params['S_SO_O2_O2'].value),
                        float(result.params['S_SO_O2_N2'].value),
                        float(result.params['S_EXCH_O2_O2'].value),
                        float(result.params['EXCH_b_O2_O2'].value),
                        float(result.params['EXCH_c_O2_O2'].value),
                        float(result.params['SO_b_O2_O2'].value),
                        float(result.params['SO_c_O2_O2'].value),
                        float(result.params['SO_b_O2_N2'].value),
                        float(result.params['SO_c_O2_N2'].value),
                        float(result.params['SO_shift_O2_O2'].value),
                        float(result.params['SO_shift_O2_N2'].value),
                        float(result.params['EXCH_shift_O2_N2'].value),)
                    spectrum.set_cia(CIA)
                    
    def generate_beta_output_file(self, beta_summary_filename = None ):
        """ Generates a file that summarizes the beta values used in the fitting in the case that beta was used to correct the Dicke narrowing term (beta_formalism = True).


        Parameters
        ----------
        beta_summary_filename : str, optional
            Filename to save the beta information. The default is Beta Summary File.

        Returns
        -------
        The generated file has the beta values for each line and spectra that were used in the fitting.  Access to this information is critical as the relationship between beta and nuVC is what generates a spoecific nuVC.

        """
        if beta_summary_filename == None:
            beta_summary_filename = 'Beta Summary File'
        if self.beta_formalism == True:
            beta_summary_list = self.lineparam_list.copy()

            for spectrum in self.dataset.spectra:
                
                mass = np.zeros(len(beta_summary_list), dtype=np.float64)
                unique_pairs = np.unique(np.column_stack((beta_summary_list['molec_id'], beta_summary_list['local_iso_id'])), axis=0)
                for m, i in unique_pairs:
                    m, i = int(m), int(i)
                    mask = (beta_summary_list['molec_id'] == m) & (beta_summary_list['local_iso_id'] == i)
                    mass[mask] = molecularMass(m, i, isotope_list=spectrum.isotope_list)
                beta_summary_list['mass'] = mass
                
                for segment in list(set(spectrum.segments)):
                    perturber_mass = 0
                    nuOptRe = np.zeros(len(beta_summary_list), dtype=np.float64)
                    beta_values = np.zeros(len(beta_summary_list), dtype=np.float64)

                    wavenumber_segments, alpha_segments, indices_segments = spectrum.segment_wave_alpha()
                    p = self.baseline_list[(self.baseline_list['Spectrum Number'] == spectrum.spectrum_number) & (self.baseline_list['Segment Number'] == segment)]['Pressure'].values[0]
                    T = self.baseline_list[(self.baseline_list['Spectrum Number'] == spectrum.spectrum_number) & (self.baseline_list['Segment Number'] == segment)]['Temperature'].values[0]

                    wave_min = np.min(wavenumber_segments[segment])
                    wave_max = np.max(wavenumber_segments[segment])

                    for species, info in spectrum.Diluent.items():
                        perturber_mass += spectrum.Diluent[species]['composition']*spectrum.Diluent[species]['m']
                        if self.constrain_dictionary['nuOptRe']:
                            nuOptRe += spectrum.Diluent[species]['composition']*(beta_summary_list['nuOptRe_%s'%species]*(p/1)*((296/T)**(beta_summary_list['n_nuOptRe_%s'%species])))
                        else:
                            nuOptRe += spectrum.Diluent[species]['composition']*(beta_summary_list['nuOptRe_%s_%s'%(species,str(spectrum.spectrum_number))]*(p/1)*((296/T)**(beta_summary_list['n_nuVC_%s'%species])))
                    alpha = perturber_mass / mass

                    if self.constrain_dictionary['nu']:
                        GammaD = np.sqrt(2*CONSTANTS['k']*CONSTANTS['Na']*T*np.log(2)/(beta_summary_list['mass'].values))*beta_summary_list['nu'] / CONSTANTS['c']
                    else:
                        GammaD = np.sqrt(2*CONSTANTS['k']*CONSTANTS['Na']*T*np.log(2)/(beta_summary_list['m_mass'].values))*beta_summary_list['nu' + '_' + str(spectrum.spectrum_number)] / CONSTANTS['c']
                    
                    mask = (beta_summary_list['nu'].values >= wave_min) & (beta_summary_list['nu'].values <= wave_max)
                    valid_indices = np.where(mask)[0]
                    beta_segment_values = np.ones(len(valid_indices), dtype=np.float64)

                    for k, ii in enumerate(valid_indices):
                        if nuOptRe[ii] > 0.0:
                            beta_segment_values[k] = beta(GammaD[ii], nuOptRe[ii], alpha[ii])
                    
                    beta_summary_list.loc[mask, 'Beta_' + str(spectrum.spectrum_number)] = beta_segment_values
                    
            select_columns = ['molec_id', 'local_iso_id', 'nu']
            for param in beta_summary_list.columns:
                if ('nuOptRe' in param) & ('n_nuOptRe' not in param):
                    select_columns.append(param)
                if ('Beta' in param):
                    select_columns.append(param)
            beta_summary_list = beta_summary_list[select_columns]
            beta_summary_list.to_csv(beta_summary_filename + '.csv')
