#Import Packages
# from .Utilities import *
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
from .spectroscopic_model import Spectroscopic_model

from lmfit import Minimizer,  Parameters

# lmfit generates warnings from the uncertainties module about params with zero uncertainty
# this is expected behavior and we should be able to safely ignore them
warnings.filterwarnings("ignore", category=UserWarning, module="uncertainties")

def convert_int_to_float(df, exclude_cols=None):
    mask = (df.drop(columns=exclude_cols, axis=1) if exclude_cols else df).select_dtypes(int)
    df[mask.columns] = mask.astype(float)
    return df


class Fit_DataSet:
    """Provides the fitting functionality for a Dataset.


    Parameters
    ----------
    dataset : object
        Dataset Object.
    base_linelist_file : str
        filename for file containing baseline parameters.
    param_linelist_file : str
        filename for file containing parmeter parameters.
    CIA_linelist_file : str, optional
        Future Feature: filename for file constraining CIA parameters
    minimum_parameter_fit_intensity : float, optional
        minimum intensity for parameters to be generated for fitting. NOTE: Even if a value is floated in the param_linelist if it is below this threshold then it won't be a floated.. The default is 1e-27.
    weight_spectra : boolean
        If True, then the pt by pt percent uncertainty for each spectrum and the spectrum weighting will be used in the calculation of the residuals.  Default is False.
    baseline_limit : bool, optional
        If True, then impose min/max limits on baseline parameter solutions. The default is False.
    baseline_limit_factor : float, optional
        The factor variable describes the multiplicative factor that the value can vary by min = 1/factor * init_guess, max = factor* init_guess. NOTE: If the init_guess for a parameter is equal to zero, then the limits aren't imposed because 1) then it would constrain the fit to 0 and 2) LMfit won't let you set min = max.. The default is 10.
    pressure_limit : bool, optional
        If True, then impose min/max limits on pressure solutions. The default is False.
    pressure_limit_factor : float, optional
        The factor variable describes the multiplicative factor the value can vary by min = 1/factor * init_guess, max = factor* init_guess. NOTE: If the init_guess for a parameter is equal to zero, then the limits aren't imposed because 1) then it would constrain the fit to 0 and 2) LMfit won't let you set min = max.. The default is 10.
    temperature_limit : bool, optional
        If True, then impose min/max limits on temperature solutions. The default is False.
    temperature_limit_factor : float, optional
        The factor variable describes the multiplicative factor the value can vary by min = 1/factor * init_guess, max = factor* init_guess. NOTE: If the init_guess for a parameter is equal to zero, then the limits aren't imposed because 1) then it would constrain the fit to 0 and 2) LMfit won't let you set min = max.. The default is 10.
    molefraction_limit : bool, optional
        If True, then impose min/max limits on mole fraction solutions. The default is False.
    molefraction_limit_factor : float, optional
        DESCRIPTION. The default is 10.
    etalon_limit : bool, optional
        If True, then impose min/max limits on etalon solutions. The default is False.
    etalon_limit_factor : float, optional
        The factor variable describes the multiplicative factor the value can vary by min = 1/factor * init_guess, max = factor* init_guess. NOTE: If the init_guess for a parameter is equal to zero, then the limits aren't imposed because 1) then it would constrain the fit to 0 and 2) LMfit won't let you set min = max.. The default is 50. #phase is constrained to +/- 2pi
    x_shift_limit : bool, optional
        If True, then impose min/max limits on x-shift solutions. The default is False.
    x_shift_limit_magnitude : float, optional
        The magnitude variables set the +/- range of the variable in cm-1. The default is 0.1.
    nu_limit : bool, optional
        If True, then impose min/max limits on line center solutions. The default is False.
    nu_limit_magnitude : float, optional
        The magnitude variables set the +/- range of the variable in cm-1. The default is 0.1.
    sw_limit : bool, optional
        If True, then impose min/max limits on line intensity solutions. The default is False.
    sw_limit_factor : float, optional
        The factor variable describes the multiplicative factor the value can vary by min = 1/factor * init_guess, max = factor* init_guess. NOTE: If the init_guess for a parameter is equal to zero, then the limits aren't imposed because 1) then it would constrain the fit to 0 and 2) LMfit won't let you set min = max.. The default is 10.
    gamma0_limit : bool, optional
        If True, then impose min/max limits on collisional half-width solutions. The default is False.
    gamma0_limit_factor : float, optional
        The factor variable describes the multiplicative factor the value can vary by min = 1/factor * init_guess, max = factor* init_guess. NOTE: If the init_guess for a parameter is equal to zero, then the limits aren't imposed because 1) then it would constrain the fit to 0 and 2) LMfit won't let you set min = max.. The default is 10.
    n_gamma0_limit : bool, optional
        DESCIf True, then impose min/max limits on temperature exponent for half width solutions. The default is True.
    n_gamma0_limit_factor : float, optional
        The factor variable describes the multiplicative factor the value can vary by min = 1/factor * init_guess, max = factor* init_guess. NOTE: If the init_guess for a parameter is equal to zero, then the limits aren't imposed because 1) then it would constrain the fit to 0 and 2) LMfit won't let you set min = max.. The default is 10.
    delta0_limit : bool, optional
        If True, then impose min/max limits on collisional shift solutions. The default is False.
    delta0_limit_factor : float, optional
        The factor variable describes the multiplicative factor the value can vary by min = 1/factor * init_guess, max = factor* init_guess. NOTE: If the init_guess for a parameter is equal to zero, then the limits aren't imposed because 1) then it would constrain the fit to 0 and 2) LMfit won't let you set min = max.. The default is 10.
    n_delta0_limit : bool, optional
        If True, then impose min/max limits on temperature exponent of the collisional shift solutions. The default is True.
    n_delta0_limit_factor : float, optional
        DESCRIPTION. The default is 10.
    SD_gamma_limit : bool, optional
        If True, then impose min/max limits on the aw solutions. The default is False.
    SD_gamma_limit_factor : float, optional
        The factor variable describes the multiplicative factor the value can vary by min = 1/factor * init_guess, max = factor* init_guess. NOTE: If the init_guess for a parameter is equal to zero, then the limits aren't imposed because 1) then it would constrain the fit to 0 and 2) LMfit won't let you set min = max.. The default is 10.
    n_gamma2_limit : bool, optional
        If True, then impose min/max limits on temperature exponent of the speed-dependent width solutions. The default is True.
    n_gamma2_limit_factor : float, optional
        The factor variable describes the multiplicative factor the value can vary by min = 1/factor * init_guess, max = factor* init_guess. NOTE: If the init_guess for a parameter is equal to zero, then the limits aren't imposed because 1) then it would constrain the fit to 0 and 2) LMfit won't let you set min = max.. The default is 10.
    SD_delta_limit : bool, optional
        If True, then impose min/max limits on as solutions. The default is True.
    SD_delta_limit_factor : float, optional
        The factor variable describes the multiplicative factor the value can vary by min = 1/factor * init_guess, max = factor* init_guess. NOTE: If the init_guess for a parameter is equal to zero, then the limits aren't imposed because 1) then it would constrain the fit to 0 and 2) LMfit won't let you set min = max.. The default is 10.
    n_delta2_limit : bool, optional
        If True, then impose min/max limits on temperature exponent of the speed-dependent shift solutions. The default is True.
    n_delta2_limit_factor : float, optional
        The factor variable describes the multiplicative factor the value can vary by min = 1/factor * init_guess, max = factor* init_guess. NOTE: If the init_guess for a parameter is equal to zero, then the limits aren't imposed because 1) then it would constrain the fit to 0 and 2) LMfit won't let you set min = max.. The default is 10.
    nuVC_limit : bool, optional
        If True, then impose min/max limits on dicke narrowing solutions. The default is False.
    nuVC_limit_factor : float, optional
        The factor variable describes the multiplicative factor the value can vary by min = 1/factor * init_guess, max = factor* init_guess. NOTE: If the init_guess for a parameter is equal to zero, then the limits aren't imposed because 1) then it would constrain the fit to 0 and 2) LMfit won't let you set min = max.. The default is 10.
    n_nuVC_limit : bool, optional
        If True, then impose min/max limits on temperature exponent of dicke narrowing solutions. The default is True.
    n_nuVC_limit_factor : float, optional
        The factor variable describes the multiplicative factor the value can vary by min = 1/factor * init_guess, max = factor* init_guess. NOTE: If the init_guess for a parameter is equal to zero, then the limits aren't imposed because 1) then it would constrain the fit to 0 and 2) LMfit won't let you set min = max.. The default is 10.
    eta_limit : bool, optional
        If True, then impose min/max limits on correlation parameter solutions.. The default is True.
    eta_limit_factor : float, optional
         The factor variable describes the multiplicative factor the value can vary by min = 1/factor * init_guess, max = factor* init_guess. NOTE: If the init_guess for a parameter is equal to zero, then the limits aren't imposed because 1) then it would constrain the fit to 0 and 2) LMfit won't let you set min = max.The default is 10.
    linemixing_limit : bool, optional
        The factor variable describes the multiplicative factor the value can vary by min = 1/factor * init_guess, max = factor* init_guess. NOTE: If the init_guess for a parameter is equal to zero, then the limits aren't imposed because 1) then it would constrain the fit to 0 and 2) LMfit won't let you set min = max.. The default is False.
    linemixing_limit_factor : float, optional
        If True, then impose min/max limits on line-mixing solutions.. The default is 10.
    beta_formalism : boolean, optional
        If True, then the beta correction on the Dicke narrowing is used in the simulation model.
    """

    def __init__(self, dataset, base_linelist_file, param_linelist_file, CIA_linelist_file = None,
                 lineprofile = 'mHTP',
                minimum_parameter_fit_intensity = 1e-30, minimum_simulation_intensity=1e-30,
                weight_spectra = False,
                baseline_limit = False, baseline_limit_factor = 10,
                pressure_limit = False, pressure_limit_factor = 10,
                temperature_limit = False, temperature_limit_factor = 10,
                molefraction_limit = False, molefraction_limit_factor = 10,
                etalon_limit = False, etalon_limit_factor = 50, #phase is constrained to +/- 2pi,
                x_shift_limit = False, x_shift_limit_magnitude = 0.1,
                nu_limit = False, nu_limit_magnitude = 0.1,

                sw_limit = False, sw_limit_factor = 10,
                gamma0_limit = False, gamma0_limit_factor = 10, n_gamma0_limit= False, n_gamma0_limit_factor = 10,
                delta0_limit = False, delta0_limit_factor = 10, n_delta0_limit = False, n_delta0_limit_factor = 10,
                SD_gamma_limit = False, SD_gamma_limit_factor  = 10, n_gamma2_limit = False, n_gamma2_limit_factor  = 10,
                SD_delta_limit = False, SD_delta_limit_factor  = 10, n_delta2_limit = False, n_delta2_limit_factor  = 10,
                nuVC_limit = False, nuVC_limit_factor  = 10, n_nuVC_limit = False, n_nuVC_limit_factor = 10,
                eta_limit = False, eta_limit_factor  = 10,

                nuOptRe_limit = False, nuOptRe_limit_factor = 10, n_nuOptRe_limit = False, n_nuOptRe_limit_factor = 10, 
                nuOptIm_limit = False, nuOptIm_limit_factor = 10, n_nuOptIm_limit = False, n_nuOptIm_limit_factor = 10,

                linemixing_limit = False, linemixing_limit_factor  = 10, n_linemixing_limit = False, n_linemixing_limit_factor = 10,
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
        raw_df = pd.read_csv(param_linelist_file + ".csv")
        if 'nu' in raw_df.columns:
            raw_df.sort_values('nu', inplace=True)
        raw_df.reset_index(drop=True, inplace=True)
        raw_df = raw_df.loc[:, ~raw_df.columns.str.contains('Unnamed:')]
        
        self.lineparam_list = convert_int_to_float(raw_df, exclude_cols=int_cols)
        self.lineprofile = lineprofile
        self.engine = Spectroscopic_model(self.lineparam_list, lineprofile = self.lineprofile, 
                                          isotope_list = self.dataset.isotope_list)
        
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
        self.baseline_limit = baseline_limit
        self.baseline_limit_factor  = baseline_limit_factor
        self.pressure_limit = pressure_limit
        self.pressure_limit_factor = pressure_limit_factor
        self.temperature_limit = temperature_limit
        self.temperature_limit_factor = temperature_limit_factor
        self.etalon_limit = etalon_limit
        self.etalon_limit_factor = etalon_limit_factor
        self.molefraction_limit = molefraction_limit
        self.molefraction_limit_factor = molefraction_limit_factor
        self.etalon_limit = etalon_limit
        self.etalon_limit_factor = etalon_limit_factor
        self.x_shift_limit = x_shift_limit
        self.x_shift_limit_magnitude = x_shift_limit_magnitude
        self.nu_limit = nu_limit
        self.nu_limit_magnitude = nu_limit_magnitude
        self.sw_limit = sw_limit
        self.sw_limit_factor = sw_limit_factor
        self.gamma0_limit = gamma0_limit
        self.gamma0_limit_factor = gamma0_limit_factor
        self.n_gamma0_limit = n_gamma0_limit
        self.n_gamma0_limit_factor = n_gamma0_limit_factor
        self.delta0_limit = delta0_limit
        self.delta0_limit_factor= delta0_limit_factor
        self.n_delta0_limit = n_delta0_limit
        self.n_delta0_limit_factor = n_delta0_limit_factor
        self.SD_gamma_limit = SD_gamma_limit
        self.SD_gamma_limit_factor = SD_gamma_limit_factor
        self.n_gamma2_limit = n_gamma2_limit
        self.n_gamma2_limit_factor = n_gamma2_limit_factor
        self.SD_delta_limit = SD_delta_limit
        self.SD_delta_limit_factor = SD_delta_limit_factor
        self.n_delta2_limit = n_delta2_limit
        self.n_delta2_limit_factor = n_delta2_limit_factor
        

        if lineprofile == 'HTP':
            self.paramRe_limit = nuVC_limit
            self.paramRe_limit_factor = nuVC_limit_factor
            self.n_paramRe_limit = n_nuVC_limit
            self.n_paramRe_limit_factor = n_nuVC_limit_factor
            self.paramIm_limit = eta_limit
            self.paramIm_limit_factor = eta_limit_factor
            self.n_paramIm_limit = False
            self.n_paramIm_limit_factor = None

        else:
            self.paramRe_limit = nuOptRe_limit
            self.paramRe_limit_factor = nuOptRe_limit_factor
            self.n_paramRe_limit = n_nuOptRe_limit
            self.n_paramRe_limit_factor = n_nuOptRe_limit_factor
            self.paramIm_limit = nuOptIm_limit
            self.paramIm_limit_factor = nuOptIm_limit_factor
            self.n_paramIm_limit = n_nuOptIm_limit
            self.n_paramIm_limit_factor = n_nuOptIm_limit_factor

        self.linemixing_limit = linemixing_limit
        self.linemixing_limit_factor = linemixing_limit_factor
        self.n_linemixing_limit = n_linemixing_limit
        self.n_linemixing_limit_factor = n_linemixing_limit_factor

        
        self.beta_formalism = beta_formalism #Beta formalism

        self.spec_attrs = self.prep_sim()

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
        # Safely scan for up to 10 etalons per spectrum
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
            # Naming convention: ILSFunc_res_0_Spec_Seg
            p_name = f"{func_name}_res_{i}_{spectrum.spectrum_number}_{segment}"
            if p_name in params:
                resolution_params.append(params[p_name].value)
            else:
                resolution_params.append(0.0) 
        return resolution_params

    def _extract_cia_config(self, params):
        """Packs CIA model and parameters into a config dict."""
        if self.dataset.CIA_model['model'] == 'Karman':
            # Extract the specific Karman parameters
            cia_params = {}
            param_map = {
                # --- Scalar Magnitudes ---
                'S_SO_O2_O2':   'SO_O2',       # Fit param S_SO_O2_O2 -> Arg SO_O2
                'S_SO_O2_N2':   'SO_N2',       # Fit param S_SO_O2_N2 -> Arg SO_N2
                'S_EXCH_O2_O2': 'EXCH_O2',     # Fit param S_EXCH_O2_O2 -> Arg EXCH_O2
                
                # --- Shape Parameters (Global) ---
                'EXCH_b_O2_O2': 'EXCH_b',
                'EXCH_c_O2_O2': 'EXCH_c',
                
                # --- Specific Shape Parameters ---
                'SO_b_O2_O2':   'SO_b_O2_O2',
                'SO_c_O2_O2':   'SO_c_O2_O2',
                'SO_b_O2_N2':   'SO_b_O2_N2',
                'SO_c_O2_N2':   'SO_c_O2_N2',
                
                # --- Shifts ---
                'SO_shift_O2_O2':   'SO_shift_O2_O2',
                'SO_shift_O2_N2':   'SO_shift_O2_N2',
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
        spectrum_attributes = {"Compressability Factor": None, 'CIA model': None} # can add all potential pre-calculated parts
        spectra_numbers = self.dataset.get_list_spectrum_numbers()
        spectra_numbers.append('Dataset')

        spectra_attribute_dict = {spec_num: spectrum_attributes.copy() for spec_num in spectra_numbers}        

        for spectrum in self.dataset.spectra:
            if spectrum.compressability_file != None:
                comp_factor = pd.read_csv(spectrum.compressability_file + '.csv')
                pressures = np.asarray(comp_factor['Pressure (MPa)'].values*1e6/101325)
                pressures = pressures.astype(float)
                temperatures = list(comp_factor)
                temperatures.remove('Pressure (MPa)')
                temperatures = np.asarray(temperatures)
                temperatures = temperatures.astype(float)
                comp_factor.drop('Pressure (MPa)', inplace=True, axis=1) 
                comp_factor_array = comp_factor.to_numpy()
                spectra_attribute_dict[spectrum.spectrum_number]['Compressability Factor'] = RegularGridInterpolator(points = [pressures, temperatures], values = comp_factor_array)
            else:
                spectra_attribute_dict[spectrum.spectrum_number]['Compressability Factor'] = None
        
        if self.dataset.CIA_model['model'] == 'Karman':
            spectra_attribute_dict['Dataset']['CIA model'] =O2_CIA_Karman_Model(self.dataset.CIA_model['band'])

        return spectra_attribute_dict
 

    def generate_params(self):
        """Generates the lmfit parameter object that will be used in fitting.


        Returns
        -------
        params : lmfit parameter object
            the parameter object is a dictionary comprised of all parameters translated from dataframes into a dictionary format compatible with lmfit

        """

        params = Parameters()
        #Baseline Parameters
        baseline_parameters = []
        for base_param in list(self.baseline_list):
            if ('_vary' not in base_param) and ('_err' not in base_param) and ('Spectrum Number' not in base_param) and ('Segment Number' not in base_param):
                baseline_parameters.append(base_param)
        for index in self.baseline_list.index.values:
            spec_num = self.baseline_list.iloc[index]['Spectrum Number']
            seg_num = self.baseline_list.iloc[index]['Segment Number']
            for base_param in baseline_parameters:
                if self.baseline_list.loc[index][base_param] == 0:
                    params.add(base_param + '_'+str(int(spec_num))+'_'+ str(int(seg_num)), self.baseline_list.loc[index][base_param], self.baseline_list.loc[index][base_param + '_vary'])
                elif ('Pressure' in base_param) and self.pressure_limit:
                    params.add(base_param + '_'+str(int(spec_num))+'_'+ str(int(seg_num)), self.baseline_list.loc[index][base_param], self.baseline_list.loc[index][base_param + '_vary'],
                              min = (1 / self.pressure_limit_factor)*self.baseline_list.loc[index][base_param],
                              max = self.pressure_limit_factor*self.baseline_list.loc[index][base_param])
                elif ('Temperature' in base_param) and self.temperature_limit:
                    params.add(base_param + '_'+str(int(spec_num))+'_'+ str(int(seg_num)), self.baseline_list.loc[index][base_param], self.baseline_list.loc[index][base_param + '_vary'],
                              min = (1 / self.temperature_limit_factor)*self.baseline_list.loc[index][base_param],
                              max = self.temperature_limit_factor*self.baseline_list.loc[index][base_param])
                elif ('molefraction' in base_param) and self.molefraction_limit:
                    params.add(base_param + '_'+str(int(spec_num))+'_'+ str(int(seg_num)), self.baseline_list.loc[index][base_param], self.baseline_list.loc[index][base_param + '_vary'],
                              min = (1 / self.molefraction_limit_factor)*self.baseline_list.loc[index][base_param],
                              max = self.molefraction_limit_factor*self.baseline_list.loc[index][base_param])
                elif ('baseline' in base_param) and self.baseline_limit:
                    params.add(base_param + '_'+str(int(spec_num))+'_'+ str(int(seg_num)), self.baseline_list.loc[index][base_param], self.baseline_list.loc[index][base_param + '_vary'],
                              min = (1 / self.baseline_limit_factor)*self.baseline_list.loc[index][base_param],
                              max = self.baseline_limit_factor *self.baseline_list.loc[index][base_param])
                elif ('etalon_' in base_param) and self.etalon_limit and ('phase' not in base_param):
                    params.add(base_param + '_'+str(int(spec_num))+'_'+ str(int(seg_num)), self.baseline_list.loc[index][base_param], self.baseline_list.loc[index][base_param + '_vary'],
                              min = (1 / self.etalon_limit_factor )*self.baseline_list.loc[index][base_param],
                              max = self.etalon_limit_factor *self.baseline_list.loc[index][base_param])
                elif ('etalon_' in base_param) and self.etalon_limit and ('phase' in base_param):
                    params.add(base_param + '_'+str(int(spec_num))+'_'+ str(int(seg_num)), self.baseline_list.loc[index][base_param], self.baseline_list.loc[index][base_param + '_vary'],
                              min = -2*np.pi, max = 2*np.pi)
                elif ('x_shift' in base_param) and self.x_shift_limit:
                    params.add(base_param + '_'+str(int(spec_num))+'_'+ str(int(seg_num)), self.baseline_list.loc[index][base_param], self.baseline_list.loc[index][base_param + '_vary'],
                              min = (self.baseline_list.loc[index][base_param] - self.x_shift_limit_magnitude),
                              max = self.x_shift_limit_magnitude + self.baseline_list.loc[index][base_param])
                else:
                    params.add(base_param + '_'+str(int(spec_num))+'_'+ str(int(seg_num)), self.baseline_list.loc[index][base_param], self.baseline_list.loc[index][base_param + '_vary'])
        
        if self.lineprofile == 'HTP':
            param_Re = 'nuVC'
            param_Im = 'eta'            
        else: #all other uses mHTP defintiion.
            param_Re = 'nuOptRe'
            param_Im = 'nuOptIm'
        
        diluent_list = []
        for spectrum in self.dataset.spectra:
            for diluent in spectrum.Diluent:
                if diluent not in diluent_list:
                    diluent_list.append(diluent)
        n_dil = len(diluent_list)

        #Lineshape parameters
        linelist_params = []
        num_nominal_temps = self.dataset.get_number_nominal_temperatures()[0]

        for line_param in list(self.lineparam_list):
            print ()
            if num_nominal_temps == 1:
                if self.dataset.BIA_model['sw_depletion']:
                    if self.dataset.BIA_model['farwing_continuum'] == 'LBL':
                        if ('_vary' not in line_param) and ('_err' not in line_param) and (line_param != 'molec_id') and (line_param != 'local_iso_id') and (line_param != 'elower') and ('n_' not in line_param):
                            linelist_params.append(line_param)
                    else:
                        if ('_vary' not in line_param) and ('_err' not in line_param) and (line_param != 'molec_id') and (line_param != 'local_iso_id') and (line_param != 'elower') and ('n_' not in line_param) and ('BIA_collision_duration' not in line_param):
                            linelist_params.append(line_param)
                else:
                    if ('_vary' not in line_param) and ('_err' not in line_param) and (line_param != 'molec_id') and (line_param != 'local_iso_id') and (line_param != 'elower') and ('n_' not in line_param) and ('BIA' not in line_param):
                        linelist_params.append(line_param)
            else:
                if self.dataset.BIA_model['sw_depletion']:
                    if self.dataset.BIA_model['farwing_continuum'] == 'LBL':
                        if ('_vary' not in line_param) and ('_err' not in line_param) and (line_param != 'molec_id') and (line_param != 'local_iso_id') and (line_param != 'elower'):
                            linelist_params.append(line_param)
                    else:
                        if ('_vary' not in line_param) and ('_err' not in line_param) and (line_param != 'molec_id') and (line_param != 'local_iso_id') and (line_param != 'elower') and ('BIA_collision_duration' not in line_param):
                            linelist_params.append(line_param)
                else:
                    if ('_vary' not in line_param) and ('_err' not in line_param) and (line_param != 'molec_id') and (line_param != 'local_iso_id') and (line_param != 'elower') and ('BIA' not in line_param):
                        linelist_params.append(line_param)


        def count_cols(key): return sum(key in p for p in linelist_params)

        nu_constrain = True
        sw_constrain = True
        if ((sum(('nu' in param) & (param_Re not in param) & (param_Im not in param) & ('n_' not in param) for param in linelist_params))) > 1:
            nu_constrain = False
        if (sum(('sw' in param) & ('sw_scale_factor' not in param) for param in linelist_params)) > 1:
            sw_constrain = False

        gamma0_constrain = True
        delta0_constrain = True
        SD_gamma_constrain = True
        SD_delta_constrain = True
        paramRe_constrain = True
        paramIm_constrain = True
        linemix_constrain = True

        limit_g0 = n_dil * 2 if num_nominal_temps > 1 else n_dil
        if count_cols('gamma0_') > limit_g0: gamma0_constrain = False
        limit_d0 = n_dil * 2 if num_nominal_temps > 1 else n_dil
        if count_cols('delta0_') > limit_d0: delta0_constrain = False

        if count_cols('SD_gamma_') > n_dil: SD_gamma_constrain = False # Should this be +1 for >1 temperature
        if count_cols('SD_delta_') > n_dil: SD_delta_constrain = False # Should this be +1 for >1 temperature

        limit_re = n_dil * 2 if num_nominal_temps > 1 else n_dil
        if count_cols(param_Re) > limit_re: paramRe_constrain = False

        limit_im = n_dil * 2 if (self.lineprofile != 'HTP' and num_nominal_temps > 1) else n_dil
        if count_cols(param_Im) > limit_im: paramIm_constrain = False

        limit_y = n_dil * 2 if num_nominal_temps > 1 else n_dil
        if count_cols('y_') > limit_y: linemix_constrain = False
        
        for spec_line in self.lineparam_list.index.values:
            sw_scaled = self.lineparam_list.loc[spec_line]['sw'] * self.lineparam_list.loc[spec_line]['sw_scale_factor']
            if sw_scaled < self.minimum_parameter_fit_intensity:
                continue
            for line_param in linelist_params:
                print (linelist_params)
                val = self.lineparam_list.loc[spec_line][line_param]
                print (line_param, val)
                vary = self.lineparam_list.loc[spec_line][line_param + '_vary']
                lmfit_name = f"{line_param}_line_{spec_line}"


                #indices = [m.start() for m in re.finditer('_', line_param)]
                #index_length = len(indices)
                #NU
                if line_param == 'nu' and nu_constrain:
                    if self.nu_limit:
                        params.add(lmfit_name, val, vary, min = val - self.nu_limit_magnitude, max = val + self.nu_limit_magnitude)
                    else:
                        params.add(lmfit_name, val, vary)
                elif (line_param != 'nu') and ('nu_' in line_param) and (param_Re not in line_param) and (param_Im not in line_param) and (not nu_constrain):
                    if self.nu_limit:
                        params.add(lmfit_name, val, vary, min = val - self.nu_limit_magnitude, max = val + self.nu_limit_magnitude)
                    else:
                        params.add(lmfit_name, val, vary)
                #SW
                elif line_param == 'sw' and sw_constrain:
                    if self.sw_limit:
                        params.add(lmfit_name, val, vary, min = (1 / self.sw_limit_factor)* val, max = self.sw_limit_factor* val)
                    else:
                        params.add(lmfit_name, val, vary)
                elif (line_param != 'sw') and ('sw' in line_param) and (not sw_constrain) and (line_param != 'sw_scale_factor'):
                    if self.sw_limit:
                        params.add(lmfit_name, val, vary,min =  (1 / self.sw_limit_factor)* val, max = self.sw_limit_factor*val)
                    else:
                        params.add(lmfit_name, val, vary)
                #GAMMA0
                elif ('gamma0_' in line_param) and ('n_' not in line_param) and (gamma0_constrain) : #and (index_length==1)
                    if self.gamma0_limit and val != 0:
                        params.add(lmfit_name, val, vary, min = (1 / self.gamma0_limit_factor)*val, max = self.gamma0_limit_factor*val)
                    else:
                        params.add(lmfit_name, val, vary)
                elif ('gamma0_' in line_param) and ('n_' not in line_param) and (not gamma0_constrain): #and (index_length>1)
                    if self.gamma0_limit and val != 0:
                            params.add(lmfit_name, val, vary, min = (1 / self.gamma0_limit_factor)*val,max = self.gamma0_limit_factor*val)
                    else:
                            params.add(lmfit_name, val, vary)
                elif ('n_gamma0' in line_param):
                    if self.n_gamma0_limit and val != 0:
                        params.add(lmfit_name, val, vary, min = (1 / self.n_gamma0_limit_factor) *val, max = self.n_gamma0_limit_factor*val)
                    else:
                        params.add(lmfit_name, val, vary)
                #DELTA0
                elif ('delta0' in line_param) and ('n_' not in line_param) and (delta0_constrain) : #and (index_length==1)
                    if self.delta0_limit and val != 0:
                        params.add(lmfit_name, val, vary, min = (1 / self.delta0_limit_factor )*val, max = self.delta0_limit_factor*val)
                    else:
                        params.add(lmfit_name, val, vary)
                elif ('delta0_' in line_param) and ('n_' not in line_param) and (not delta0_constrain) : #and (index_length>1)
                    if self.delta0_limit and val != 0:
                            params.add(lmfit_name, val, vary, min = (1 / self.delta0_limit_factor)*val, max = self.delta0_limit_factor*val)
                    else:
                            params.add(lmfit_name, val, vary)
                elif ('n_delta0' in line_param):
                    if self.n_delta0_limit and val != 0:
                        params.add(lmfit_name, val, vary, min = (1 / self.n_delta0_limit_factor)*val, max = self.n_delta0_limit_factor / 100*val)
                    else:
                        params.add(lmfit_name, val, vary)
                #SD Gamma
                elif ('SD_gamma' in line_param) and (SD_gamma_constrain): #and (index_length==2)
                    if self.SD_gamma_limit and val != 0:
                        params.add(lmfit_name, val, vary, min = (1 / self.SD_gamma_limit_factor) *val, max = self.SD_gamma_limit_factor*val)
                    else:
                        params.add(lmfit_name, val, vary)
                elif ('SD_gamma' in line_param) and (not SD_gamma_constrain): #and (index_length>2)
                    if self.SD_gamma_limit and val != 0:
                            params.add(lmfit_name, val, vary, min = (1 / self.SD_gamma_limit_factor)*val, max = self.SD_gamma_limit_factor*val)
                    else:
                            params.add(lmfit_name, val, vary)
                elif ('n_gamma2' in line_param):
                    if self.n_gamma2_limit and val != 0:
                        params.add(lmfit_name, val, vary, min = (1 / self.n_gamma2_limit_factor)*val, max = (self.n_gamma2_limit_factor / 100)*val)
                    else:
                        params.add(lmfit_name, val, vary)
                #SD Delta
                elif ('SD_delta' in line_param) and (SD_delta_constrain): #and (index_length==2)
                    if self.SD_delta_limit and val != 0:
                        params.add(lmfit_name, val, vary, min = (1 / self.SD_delta_limit_factor )*val, max = self.SD_delta_limit_factor*val)
                    else:
                        params.add(line_param + '_' + 'line_' + str(spec_line), val, vary)
                elif ('SD_delta' in line_param) and (not SD_delta_constrain): #and (index_length>2)
                    if self.SD_delta_limit and val != 0:
                            params.add(lmfit_name, val, vary,
                                min = (1 / self.SD_delta_limit_factor )*val,
                                max =self.SD_delta_limit_factor * val)
                    else:
                            params.add(lmfit_name, val, vary)
                elif ('n_delta2' in line_param):
                    if self.n_delta2_limit and val != 0:
                        params.add(lmfit_name, val, vary, min = (1 / self.n_delta2_limit_factor )*val, max = self.n_delta2_limit_factor*val)
                    else:
                        params.add(lmfit_name, val, vary)
                #nuVC
                elif (param_Re in line_param) and ('n_'+ param_Re not in line_param) and (paramRe_constrain): #and (index_length==1)
                    if self.paramRe_limit and val!= 0:
                        params.add(lmfit_name, val, vary, min = (1 /self.paramRe_limit_factor)*val, max = self.paramRe_limit_factor*val)
                    else:
                        params.add(lmfit_name, val, vary)
                elif (param_Re in line_param) and ('n_' + param_Re not in line_param) and (not paramRe_constrain): #(index_length>1)
                    if self.paramRe_limit and val != 0:
                            params.add(lmfit_name, val, vary, min = (1 / self.paramRe_limit_factor)*val, max = self.paramRe_limit_factor*val)  
                    else:
                        params.add(lmfit_name, val, vary)

                elif ('n_' + param_Re in line_param):
                    if self.n_paramRe_limit and val != 0:
                        params.add(lmfit_name, val, vary, min = (1 / self.n_paramRe_limit_factor )*val, max = self.n_paramRe_limit_factor*val)
                    else:
                        params.add(lmfit_name, val, vary)                
                #eta
                elif (param_Im in line_param) and ('n_'  + param_Im not in line_param) and (paramIm_constrain): # and (index_length==1)
                    if self.paramIm_limit and val != 0:
                        params.add(lmfit_name, val, vary, min = (1 / self.paramIm_limit_factor)*val, max = (self.paramIm_limit_factor)*val)
                    else:
                        params.add(lmfit_name, val, vary)
                elif (param_Im in line_param) and ('n_'  + param_Im not in line_param) and (not paramIm_constrain): # and (index_length>1)
                    if self.paramIm_limit and val != 0:
                            params.add(lmfit_name, val, vary, min = (1 / self.paramIm_limit)*val, max = self.paramIm_limit*val)
                    else:
                            params.add(lmfit_name, val, vary)

                elif ('n_' + param_Im in line_param):
                    if self.n_paramIm_limit and val != 0:
                        params.add(lmfit_name, val, vary, min = (1 / self.n_paramIm_limit_factor )*val, max = self.n_paramIm_limit_factor*val)
                    else:
                        params.add(lmfit_name, val, vary)
                
                # linemixing
                
                elif ('y_' in line_param) and ('n_' not in line_param) and (linemix_constrain): # and (index_length==1)
                    if self.linemixing_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                        params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'],
                                min = (1 / self.linemixing_limit_factor)*self.lineparam_list.loc[int(spec_line)][line_param] ,
                                max = self.linemixing_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])
                    else:
                        params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                elif ('y_' in line_param) and ('n_' not in line_param) and (not linemix_constrain): # and (index_length>1)
                    if self.linemixing_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'],
                                min = (1 / self.linemixing_limit_factor)*self.lineparam_list.loc[int(spec_line)][line_param],
                                max = self.linemixing_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])
                    else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                elif ('n_y_' in line_param):
                    if self.n_linemixing_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                        params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'],
                                min = (1 / self.n_linemixing_limit_factor) *self.lineparam_list.loc[int(spec_line)][line_param],
                                max = self.n_linemixing_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])
                    else:
                        params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param],self.lineparam_list.loc[spec_line][line_param + '_vary'])
                
                #BIA
                elif ('BIA_slope_' in line_param) and (sw_constrain):
                    params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                #BIA farwing
                elif ('BIA_collision_duration_' in line_param) and (sw_constrain):
                    params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])

        #CIA Parameters (O2 Karman Model)
        if self.dataset.CIA_model['model'] == "Karman":
            cia_parameters = []
            for cia_param in list(self.CIAparam_list):
                if ('_vary' not in cia_param) and ('_err' not in cia_param) and ('CIA Pair' not in cia_param):
                    cia_parameters.append(cia_param)
            for cia_pair in self.CIAparam_list['CIA Pair']:
                index = self.CIAparam_list[self.CIAparam_list['CIA Pair']==cia_pair].index.values[0]   
                for cia_param in cia_parameters:
                    if 'SO' in cia_param:
                        params.add(cia_param + '_'+ cia_pair, self.CIAparam_list.loc[index][cia_param], 
                                   self.CIAparam_list.loc[index][cia_param + '_vary'])
                    elif 'EXCH' in cia_param:
                        if cia_pair == 'O2_O2':
                            params.add(cia_param + '_'+ cia_pair, self.CIAparam_list.loc[index][cia_param], 
                                   self.CIAparam_list.loc[index][cia_param + '_vary'])
                        elif cia_pair =='O2_N2':
                            params.add(cia_param + '_'+ cia_pair, 0, False)

        self.engine.configure_for_fitting(params)
        return (params)

    def constrained_baseline(self, params, baseline_segment_constrained = True, xshift_segment_constrained = True, molefraction_segment_constrained = True,
                                    etalon_amp_segment_constrained = True, etalon_period_segment_constrained = True, etalon_phase_segment_constrained = True,
                                    pressure_segment_constrained = True, temperature_segment_constrained = True):
        """Imposes baseline constraints when using multiple segments per spectrum, ie all baseline parameters can be the same across the entire spectrum except for the etalon phase, which is allowed to vary per segment.


        Parameters
        ----------
        params : lmfit parameter object
            the params object is a dictionary comprised of all parameters translated from dataframes into a dictionary format compatible with lmfit.
        baseline_segment_constrained : bool, optional
            True means the baseline terms are constrained across each spectrum. The default is True.
        xshift_segment_constrained : bool, optional
            True means the x_shift terms are constrained across each spectrum. The default is True.
        molefraction_segment_constrained : bool, optional
            True means the mole fraction for that molecule is constrained across each spectrum. The default is True.
        etalon_amp_segment_constrained : bool, optional
            True means the etalon amplitude is constrained across each spectrum. The default is True.
        etalon_period_segment_constrained : bool, optional
            True means the etalon period is constrained across each spectrum. The default is True.
        etalon_phase_segment_constrained : bool, optional
            True means the etalon phase is constrained across each spectrum. The default is True.
        pressure_segment_constrained : bool, optional
            True means the pressure is constrained across each spectrum. The default is True.
        temperature_segment_constrained : bool, optional
            True means the temperature is constrained across each spectrum. The default is True.

        Returns
        -------
        params : lmfit parameter object
            the params object is a dictionary comprised of all parameters translated from dataframes into a dictionary format compatible with lmfit.

        """

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
        '''All user to set the temperature dependence of the SO mechanism scalar and/or the shift of the SO to be equal for O2-O2 and O2_N2

        Parameters
        ----------
        params : lmfit parameter object
            the params object is a dictionary comprised of all parameters translated from dataframes into a dictionary format compatible with lmfit.
        S_temperature_dependence_constrained : boolean, optional
            Constrains the temperature dependence of the SO mechanism scalar to be the same for O2-N2 and O2-O2. The default is True.
        shift_constrained : boolean, optional
            Constrains the SO shift to be the same for O2-N2 and O2-O2. The default is True.

        Returns
        -------
        params : lmfit parameter object
            the params object is a dictionary comprised of all parameters translated from dataframes into a dictionary format compatible with lmfit.

        '''
        if self.dataset.CIA_model['model'] == 'Karman':
        
            for param in params:
                if (('SO_b' in param) or ('SO_c' in param)) and S_temperature_dependence_constrained:
                    if 'O2_O2' not in param:
                        params[param].set(expr = param[:5] + 'O2_O2')
                if ('SO_shift' in param) and shift_constrained:
                    if 'O2_O2' not in param:
                        params[param].set(expr = param[:9] + 'O2_O2')
        return params

    def objective_function(self, params, wing_cutoff = 25, wing_wavenumbers=25, wing_method='wing_wavenumbers'):
        """
        Calculates residuals for the Minimizer.
        """

        #Update arrays in spectroscopic_model
        self.engine.update_from_lmfit(params)
        total_residuals = []
        cia_config = None
        if self.dataset.CIA_model['model'] == 'Karman':
            cia_config = self._extract_cia_config(params) #Get CIA parameters
      

        #Loop over spectra
        for spectrum in self.dataset.spectra:
            wavenumber_segments, alpha_segments, indices_segments = spectrum.segment_wave_alpha()
            #Loop over segments
            for segment in set(spectrum.segments):

                #x_shift
                x_shift_name = f'x_shift_{spectrum.spectrum_number}_{segment}'
                x_shift = params[x_shift_name].value if x_shift_name in params else 0.0
                waves = wavenumber_segments[segment] + x_shift

                #Environmental Parameters
                T = params[f'Temperature_{spectrum.spectrum_number}_{segment}'].value
                p = params[f'Pressure_{spectrum.spectrum_number}_{segment}'].value

                mf = spectrum.molefraction.copy()
                for molec_id in mf:
                    iso_name = self.dataset.isotope_list.get((molec_id, 1), [0,0,0,0,'Unknown'])[4]
                    mf_param_name = f'molefraction_{iso_name}_{spectrum.spectrum_number}_{segment}'
                    if mf_param_name in params:
                        mf[molec_id] = params[mf_param_name].value

                baseline_coeffs = self._extract_baseline_coeffs(params, spectrum.spectrum_number, segment)
                etalon_dict = self._extract_etalon_dict(params, spectrum.spectrum_number, segment)
                ils_res = self._extract_ils_resolution(params, spectrum, segment)

                if self.dataset.CIA_model['model'] == 'ad hoc':
                    idx_min = np.min(indices_segments[segment])
                    idx_max = np.max(indices_segments[segment])
                    cia_slice = spectrum.cia[idx_min : idx_max + 1]
                    cia_config = {
                        'model': 'ad hoc',
                        'values': cia_slice}

                model_y = self.engine.calculate_spectrum(
                    waves=waves,
                    T=T, 
                    p=p, 
                    molefraction=mf,
                    Diluent=spectrum.Diluent,
                    spectrum_number=spectrum.spectrum_number,
                    spectrum_min=np.min(spectrum.wavenumber), # For baseline relative shift logic
                    
                    baseline_coeffs=baseline_coeffs,
                    etalon_dict=etalon_dict,
                    cia_config=cia_config,
                    
                    ILS_function=spectrum.ILS_function,
                    ILS_parameters=ils_res,
                    ILS_wing=spectrum.ILS_wing,
                    
                    # Engine pass-throughs
                    interpolated_compressability_file=self.spec_attrs[spectrum.spectrum_number]['Compressability Factor'],
                    BIA_slope=self.dataset.BIA_model['sw_depletion'],
                    BIA_FW_LBL=(self.dataset.BIA_model['farwing_continuum'] == 'LBL'),
                    IntensityThreshold=self.minimum_simulation_intensity, # Use class attribute
                    wing_cutoff=wing_cutoff,
                    wing_wavenumbers=wing_wavenumbers,
                    wing_method=wing_method
                )

                data_y = alpha_segments[segment]
                resid = data_y - model_y  #Earlier versions of MATS did model - data, we are switching this to the more excepted obs - calc

                #Weighting
                if self.weight_spectra:
                    idx_min = np.min(indices_segments[segment])
                    idx_max = np.max(indices_segments[segment])
                    if spectrum.tau_stats.all() == 0:
                        w = spectrum.weight
                    else:
                        seg_stats = spectrum.tau_stats[idx_min : idx_max + 1]
                        w = spectrum.weight * (1.0 / seg_stats)
                    
                    resid = resid * w

                total_residuals.append(resid)
        
        return np.concatenate(total_residuals)
                    


    def fit_data(self, params, wing_cutoff = 25, wing_wavenumbers = 25, wing_method = 'wing_wavenumbers', xtol = 1e-7, maxfev = 2000, ftol = 1e-7, 
                 method = 'least_squares'):
        """Uses the lmfit minimizer to do the fitting through the simulation model function.


        Parameters
        ----------
        params : lmfit parameter object
            the params object is a dictionary comprised of all parameters translated from dataframes into a dictionary format compatible with lmfit.
        wing_cutoff : float, optional
            number of voigt half-widths to simulate on either side of each line. The default is 50.
        wing_wavenumbers : float, optional
            number of wavenumbers to simulate on either side of each line. The default is 50.
        wing_method : str, optional
            Provides choice between the wing_cutoff and wing_wavenumbers line cut-off options. The default is 'wing_wavenumbers'.
        xtol : float, optional
             Absolute error in xopt between iterations that is acceptable for convergence. The default is 1e-7.
        maxfev : float, optional
            DESCRIPTION. The default is 2000.
        ftol : The maximum number of calls to the function., optional
            Absolute error in func(xopt) between iterations that is acceptable for convergence.. The default is 1e-7.
        method : string, optional
            Defines the minimization method from the options in LMFIT (not all will work).  Has been tested on the 'leastsq' which uses the Levenberg-Marquardt algorithm and 'least_squares' which uses the Trust Region Reflective method.

        Returns
        -------
        result : LMFit result Object
            contains all fit results as LMFit results object.

        """

        fcn_args = (wing_cutoff, wing_wavenumbers, wing_method)

        minner = Minimizer(self.objective_function, params, 
                           xtol = xtol, max_nfev = maxfev, ftol = ftol, 
                           fcn_args = fcn_args)
        
        floated_parameters = any(p.vary for p in params.values())
        if not floated_parameters and method == 'least_squares':
            method = 'leastsq'
        result = minner.minimize(method = method)
        return result
        


    def residual_analysis(self, result, indv_resid_plot = False):
        """Updates the model and residual arrays in each spectrum object with the results of the fit and gives the option of generating the combined absorption and residual plot for each spectrum.


        Parameters
        ----------
        result : LMFit result Object
            contains all fit results as LMFit results object.
        indv_resid_plot : bool, optional
            True if you want to show residual plots for each spectra.. The default is False.


        """
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
            spectrum.set_model(spectrum_residual + spectrum.alpha)
            if indv_resid_plot:
                spectrum.plot_model_residuals()
    def update_params(self, result, base_linelist_update_file = None , param_linelist_update_file = None, 
                      CIA_linelist_update_file = None):
        """Updates the baseline and line parameter files based on fit results with the option to write over the file (default) or save as a new file and updates baseline values in the spectrum objects.


        Parameters
        ----------
        result : LMFit result Object
            contains all fit results as LMFit results object.
        base_linelist_update_file : str, optional
            Name of file to save the updated baseline parameters. Default is to override the input. The default is None.
        param_linelist_update_file : str, optional
            Name of file to save the updated line parameters. Default is to override the input. The default is None.
        cia_linelist_update_file : str, optional
            Name of file to save the updated CIA parameters.  Default is to override the input.   

        """

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
                self.baseline_list.loc[(self.baseline_list['Segment Number'] == segment) & (self.baseline_list['Spectrum Number'] == spectrum), parameter] = par.value
                if par.vary:
                    self.baseline_list.loc[(self.baseline_list['Segment Number'] == segment) & (self.baseline_list['Spectrum Number'] == spectrum), parameter + '_err'] = par.stderr

            elif ('molefraction' in par.name) or ('baseline' in par.name) or ('x_shift' in par.name):
                indices = [m.start() for m in re.finditer('_', par.name)]
                parameter = (par.name[:indices[1]])
                spectrum = int(par.name[indices[1] + 1:indices[2]])
                segment = int(par.name[indices[2] + 1:])
                self.baseline_list.loc[(self.baseline_list['Segment Number'] == segment) & (self.baseline_list['Spectrum Number'] == spectrum), parameter] = par.value
                if par.vary:
                    self.baseline_list.loc[(self.baseline_list['Segment Number'] == segment) & (self.baseline_list['Spectrum Number'] == spectrum), parameter + '_err'] = par.stderr
            elif ('etalon' in par.name):
                indices = [m.start() for m in re.finditer('_', par.name)]
                parameter = par.name[:indices[2]]
                spectrum = int(par.name[indices[2]+1:indices[3]])
                segment = int(par.name[indices[3]+1:])
                self.baseline_list.loc[(self.baseline_list['Segment Number'] == segment) & (self.baseline_list['Spectrum Number'] == spectrum), parameter] = par.value
                if par.vary:
                    self.baseline_list.loc[(self.baseline_list['Segment Number'] == segment) & (self.baseline_list['Spectrum Number'] == spectrum), parameter + '_err'] = par.stderr
            elif ('_res_' in par.name):
                indices = [m.start() for m in re.finditer('_', par.name[par.name.find('_res_') + 5:])]
                spectrum = par.name[par.name.find('_res_') + 5:][indices[0]+1:indices[1]]
                segment = par.name[par.name.find('_res_') + 5:][indices[1]+1:]
                parameter = par.name[:par.name.find('_res_')] + '_res_' + par.name[par.name.find('_res_') + 5:][:indices[0]]
                self.baseline_list.loc[(self.baseline_list['Segment Number'] == segment) & (self.baseline_list['Spectrum Number'] == spectrum), parameter] = par.value
                if par.vary:
                    self.baseline_list.loc[(self.baseline_list['Segment Number'] == segment) & (self.baseline_list['Spectrum Number'] == spectrum), parameter + '_err'] = par.stderr
            #CIA
            #Karman CIA
            elif (self.dataset.CIA_model['model']=='Karman') and ('O2_O2' in par.name):
                indices = [m.start() for m in re.finditer('_', par.name)]
                parameter = par.name[:indices[1]]
                self.CIAparam_list.loc[(self.CIAparam_list['CIA Pair'] == 'O2_O2'), parameter] = par.value
                if par.vary:
                    self.CIAparam_list.loc[(self.CIAparam_list['CIA Pair'] == 'O2_O2'), parameter + '_err'] = par.stderr
            elif (self.dataset.CIA_model['model']=='Karman') and ('O2_N2' in par.name):
                indices = [m.start() for m in re.finditer('_', par.name)]
                parameter = par.name[:indices[1]]
                self.CIAparam_list.loc[(self.CIAparam_list['CIA Pair'] == 'O2_N2'), parameter] = par.value
                if par.vary:
                    self.CIAparam_list.loc[(self.CIAparam_list['CIA Pair'] == 'O2_N2'), parameter + '_err'] = par.stderr 
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
                baseline_param_array = baseline_param_array[::-1] # reverses array to be used for polyval
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
            #List of all Diluents in Dataset
            diluent_list = []
            for spectrum in self.dataset.spectra:
                for diluent in spectrum.Diluent:
                    if diluent not in diluent_list:
                        diluent_list.append(diluent)
            #Determine if line center and nuVC terms are constrained across Dataset or vary for each spectrum
            nu_constrain = True
            nuVC_constrain = True


            if ((sum(('nu' in param) & ('nuVC' not in param) &( '_err' not in param) & ('_vary' not in param ) for param in beta_summary_list))) > 1:
                nu_constrain = False
            if (self.dataset.get_number_nominal_temperatures()[0]) == 1:
                if ((sum(('nuVC' in param) &( '_err' not in param) & ('_vary' not in param ) & ('n_nuVC' not in param) for param in beta_summary_list))) > len(diluent_list):
                    nuVC_constrain = False
                    print (sum(('nuVC' in param) &( '_err' not in param) & ('_vary' not in param ) & ('_vary' not in param ) & ('n_nuVC' not in param) for param in beta_summary_list))
            else:
                if ((sum(('nuVC' in param) &( '_err' not in param) & ('_vary' not in param ) & ('n_nuVC' not in param) for param in beta_summary_list))) > 2*len(diluent_list):
                    nuVC_constrain = False
            #Add Column for mass
            for molec in beta_summary_list['molec_id'].unique():
                for iso in beta_summary_list ['local_iso_id'].unique():
                    beta_summary_list.loc[(beta_summary_list['molec_id']==molec) & (beta_summary_list['local_iso_id']==iso), 'm_mass'] = molecularMass(molec,iso, isotope_list = self.dataset.isotope_list, )

            #Single or MS for nu and nuVC
            for spectrum in self.dataset.spectra:
                wavenumber_segments, alpha_segments, indices_segments = spectrum.segment_wave_alpha()
                mp = 0
                for diluent in spectrum.Diluent:
                    mp += spectrum.Diluent[diluent]['composition']*spectrum.Diluent[diluent]['m']

                for segment in list(set(spectrum.segments)):
                    p = self.baseline_list[(self.baseline_list['Spectrum Number'] == spectrum.spectrum_number) & (self.baseline_list['Segment Number'] == segment)]['Pressure'].values[0]
                    T = self.baseline_list[(self.baseline_list['Spectrum Number'] == spectrum.spectrum_number) & (self.baseline_list['Segment Number'] == segment)]['Temperature'].values[0]
                    wave_min = np.min(wavenumber_segments[segment])
                    wave_max = np.max(wavenumber_segments[segment])

                    beta_summary_list['alpha'] = mp / beta_summary_list['m_mass']
                    if nu_constrain:
                        GammaD = np.sqrt(2*CONSTANTS['k']*CONSTANTS['Na']*T*np.log(2)/(beta_summary_list['m_mass'].values))*beta_summary_list['nu'] / CONSTANTS['c']#change with nu
                    else:
                        GammaD = np.sqrt(2*CONSTANTS['k']*CONSTANTS['Na']*T*np.log(2)/(beta_summary_list['m_mass'].values))*beta_summary_list['nu' + '_' + str(spectrum.spectrum_number)] / CONSTANTS['c'] #change with nu
                    nuVC = len(GammaD)*[0]
                    for diluent in spectrum.Diluent:
                        abun = spectrum.Diluent[diluent]['composition']
                        if nuVC_constrain:
                            nuVC += abun*(beta_summary_list['nuVC_%s'%diluent]*(p/1)*((296/T)**(beta_summary_list['n_nuVC_%s'%diluent])))
                        else:
                            nuVC += abun*(beta_summary_list['nuVC_%s'%diluent]*(p/1)*((296/T)**(beta_summary_list['n_nuVC_%s_%s'%(diluent,str(spectrum.spectrum_number))])))
                    Chi = nuVC/ GammaD
                    A = 0.0534 + 0.1585*np.exp(-0.451*beta_summary_list['alpha'].values)
                    B = 1.9595 - 0.1258*beta_summary_list['alpha'].values + 0.0056*beta_summary_list['alpha'].values**2 + 0.0050*beta_summary_list['alpha'].values**3
                    C = -0.0546 + 0.0672*beta_summary_list['alpha'].values - 0.0125*beta_summary_list['alpha'].values**2+0.0003*beta_summary_list['alpha'].values**3
                    D = 0.9466 - 0.1585*np.exp(-0.4510*beta_summary_list['alpha'].values)
                    beta_summary_list['beta'] =  A*np.tanh(B * np.log10(Chi) + C) + D
                    beta_summary_list.loc[(beta_summary_list['nu'] >= wave_min) & (beta_summary_list['nu'] <= wave_max),'Beta_' + str(spectrum.spectrum_number) ] = beta_summary_list[(beta_summary_list['nu'] >= wave_min) & (beta_summary_list['nu'] <= wave_max)]['beta'].values
                    beta_summary_list.loc[(beta_summary_list['nu'] >= wave_min) & (beta_summary_list['nu'] <= wave_max),'Chi_' + str(spectrum.spectrum_number) ] = beta_summary_list[(beta_summary_list['nu'] >= wave_min) & (beta_summary_list['nu'] <= wave_max)]['beta'].values

            select_columns = ['molec_id', 'local_iso_id', 'nu']
            for param in beta_summary_list.columns:
                if ('nuVC' in param) & ('n_nuVC' not in param):
                    select_columns.append(param)
                if ('Beta' in param):
                    select_columns.append(param)
            beta_summary_list = beta_summary_list[select_columns]
            #beta_summary_list  = beta_summary_list[(beta_summary_list['nu'] >= wave_min) & (beta_summary_list['nu'] <= wave_max)]
            beta_summary_list.to_csv(beta_summary_filename + '.csv')
