import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.stats import norm
import re

from .spectrum import simulate_spectrum
from .dataset import Dataset

from .utilities import (
    isotope_list_molecules_isotopes,
)

from .hapi import ISO

def sigma_to_alpha(sigma_level):
        alpha = 2 * norm.sf(sigma_level)
        
        return alpha

class MonteCarlo:

    def __init__(self, dataset, base_linelist_file, param_linelist_file, directory = None):
        #Assumes parameters have been saved as .csvs
        # Currently written for one dataset type
        # Currently assumes that you aren't using the segment feature or the ILS features
        # This also currently won't adjust CIA parameters



        self.dataset = dataset
        self.dataset_name= self.dataset.dataset_name



        self.directory = directory

        self.base_linelist_file = base_linelist_file
        self.param_linelist_file = param_linelist_file

        if self.directory == None:

            self.base_linelist = pd.read_csv(base_linelist_file + '.csv')
            self.param_linelist = pd.read_csv(param_linelist_file + '.csv')  
            self.fit_summary = pd.read_csv(self.dataset_name + '.csv') #assumes only one dataset type
        
            
        else:
            self.base_linelist = pd.read_csv(os.path.join(self.directory, base_linelist_file + '.csv'))
            self.param_linelist = pd.read_csv(os.path.join(self.directory, param_linelist_file + '.csv'))
            self.fit_summary = pd.read_csv(os.path.join(self.directory, self.dataset_name  + '.csv'))

        self.iterations_param_linelist = pd.DataFrame()
        self.iterations_base_linelist = pd.DataFrame()


    def _extract_baseline_params(self, spec_base_linelist):
        """Reconstructs polynomial coeffs list [c_0, ..., c_n] from baseline_params."""
        baseline_params= []
        for i in range(self.dataset.baseline_order + 1):
            char = chr(97 + i) # 0->'a', 1->'b'
            p_name = f"baseline_{char}"
            baseline_params.append(spec_base_linelist[p_name].values[0])
        return np.asarray(baseline_params)
    
    def _extract_molefraction_parms(self, spec_base_linelist, spectrum):
        molefraction_dict = spectrum.molefraction
        for molec_id in molefraction_dict:
            mol_name = self.dataset.isotope_list.get((molec_id, 1), [0,0,0,0,'Unknown'])[4]
            mf_param_name = f'molefraction_{mol_name}'
            molefraction_dict[molec_id] = spec_base_linelist[mf_param_name].values[0]
        return molefraction_dict
    
    def _extract_isoabund_parms(self, spec_base_linelist, spectrum):
        abundance_ratio_dict = spectrum.abundance_ratio_MI
        for molec_id, isotope_dict in abundance_ratio_dict.items():
            mol_name = self.dataset.isotope_list.get((molec_id, 1), [0,0,0,0,'Unknown'])[4]
            for isotope_id in isotope_dict:
                isotope_param_name = f'abundance_ratio_{mol_name}_{isotope_id}' 
                abundance_ratio_dict[molec_id][isotope_id] = spec_base_linelist[isotope_param_name].values[0]

        return abundance_ratio_dict
    
    def _extract_etalon_params(self, spec_base_linelist, spectrum):
        etalon_dict = spectrum.etalons
        if etalon_dict != {}:
            for etalon in etalon_dict:
                etalon_prefix =  f'etalon_{etalon}_'
                etalon_dict[etalon] = [spec_base_linelist[etalon_prefix + '_amp'].values[0], 
                                       spec_base_linelist[etalon_prefix + '_period'].values[1]]
        return etalon_dict

  
    def simulate_dataset(self, at_SNR = True, from_parametric_bs = False, iteration = 0,
                         beta_formalism = False,
                         wave_error = 0.0, temperature_err =  {'bias': 0, 'function': None, 'params': {}}, 
                         pressure_err = {'per_bias': 0, 'function': None, 'params': {}},
                         molefraction_err = {},
                         wing_method = 'wing_wavenumbers', wing_cutoff = 25, wing_wavenumbers = 25, sim_window = 5):
    

        list_spec = []

        if from_parametric_bs:
            baseline_linelist_for_sim = self.iterations_base_linelist[self.iterations_base_linelist['Iteration'] == iteration].copy()
            parameter_linelist_for_sim = self.iterations_param_linelist[self.iterations_param_linelist['Iteration'] == iteration].copy()
        
        else:
            baseline_linelist_for_sim = self.base_linelist.copy()
            parameter_linelist_for_sim = self.param_linelist.copy()


 
        for spectrum in self.dataset.spectra:
            spec_num = spectrum.spectrum_number
            spec_fit_summary = self.fit_summary[self.fit_summary['Spectrum Number']==spec_num].copy()

            if at_SNR:
                SNR = spec_fit_summary['QF'].mean()
            else:
                SNR = None

            spec_base_linelist = baseline_linelist_for_sim[baseline_linelist_for_sim['Spectrum Number']==spec_num].copy()
            baseline_params = self._extract_baseline_params(spec_base_linelist)
            if spectrum.dataspace == 'alpha':
                pathlength = 0
            else:
                pathlength = spec_base_linelist['pathlength'].values[0]
            

            molefraction_dict = self._extract_molefraction_parms(spec_base_linelist, spectrum)
            isotopic_abund_dict = self._extract_isoabund_parms(spec_base_linelist, spectrum)
            etalon_dict = self._extract_etalon_params(spec_base_linelist, spectrum)

            #generate simulated param_list for parameter_linelist this is where to pertub the parameters
            #Regenerate sw_based on sw_scale
            spec_param_linelist = parameter_linelist_for_sim.copy()

            lineshape_parameters = ['nu', 'sw', 'gamma0', 'delta0', 'SD_gamma', 'SD_delta', 'nuVC', 'nuOptRe', 'eta', 'nuOptIm', 'y']
            broadener_suffix = self.dataset.get_broadener_list()   
            column_suffix_to_ignore = ['_err', '_vary']

            suffix = '|'.join(map(re.escape, column_suffix_to_ignore))

            spec_param_linelist = spec_param_linelist.drop(columns=spec_param_linelist.columns[spec_param_linelist.columns.str.contains(suffix, case=False)])
            for ls in lineshape_parameters:
                if (ls == 'nu') | (ls == 'sw'):
                    if ls + '_' + str(spec_num) in list(spec_param_linelist):
                        spec_param_linelist[ls] = spec_param_linelist[ls + '_' + str(spec_num)].values # replaces with spectrum value if not constrained
                    if ls == 'sw':
                        spec_param_linelist[ls] *= spec_param_linelist['sw_scale_factor']          
                else:
                    for broadener in broadener_suffix:
                        if ls + '_' + broadener + '_' + str(spec_num) in list(spec_param_linelist):
                            spec_param_linelist[ls + '_' + broadener] = spec_param_linelist[ls + '_' + broadener + '_' + str(spec_num)].values # replaces with spectrum value if not constrained


            #If fittable parameter need to pull from Baseline and Parameter linelist files

            sim_spec = simulate_spectrum(spec_param_linelist, 
                                         lineprofile = 'mHTP', numba_lineprofile = False,
                                         wavenumbers = spec_fit_summary['Wavenumber (cm-1)'].values,  wave_error = wave_error,
                                         SNR = SNR, 
                                         dataspace = spectrum.dataspace, 
                                         temperature = spec_base_linelist['Temperature'].values[0], temperature_err =temperature_err,
                                         pressure = spec_base_linelist['Pressure'].values[0],   pressure_err = pressure_err, 


                                         baseline_terms = baseline_params, 
                                         pathlength = pathlength, 
                                         molefraction = molefraction_dict, molefraction_err = molefraction_err,
                                         abundance_ratio_MI = isotopic_abund_dict,
                                         Diluent= spectrum.Diluent,
                                         nominal_temperature = spectrum.nominal_temperature, 
                                         etalons = etalon_dict, x_shift = spec_base_linelist['x_shift'].values[0], 
                                         beta_formalism = beta_formalism, TIPS = spectrum.TIPS, compressability_file= spectrum.compressability_file,
                                         BIA_model = self.dataset.BIA_model, CIA_model = self.dataset.CIA_model,
                                         
                                         wing_cutoff = wing_cutoff, wing_wavenumbers = wing_wavenumbers, wing_method = wing_method, sim_window = sim_window, 
                                         filename = 'temp', isotope_list = self.dataset.isotope_list, natural_abundance=spectrum.natural_abundance,)

            list_spec.append(sim_spec)
        return Dataset(list_spec, self.dataset_name + '_MC', self.dataset.param_linelist, BIA_model=self.dataset.BIA_model)


    def update_fit_files(self, initial_at_optimized = False):

        parameter_linelist = self.param_linelist.copy()
        baseline_linelist = self.base_linelist.copy()

        for col in list(parameter_linelist):
            if '_vary' in col:
                param = col.replace('_vary', "")
                if not initial_at_optimized:
                    parameter_linelist[param] += np.random.normal(loc = 0.0, scale = parameter_linelist[param + '_err'], size = len(parameter_linelist) )
                parameter_linelist[param + '_vary'] = False
                parameter_linelist.loc[parameter_linelist[param + '_err']!=0, param + '_vary'] = True
        
        for col in list(baseline_linelist):
            if '_vary' in col:
                param = col.replace('_vary', "")
                if not initial_at_optimized:
                    baseline_linelist[param] += np.random.normal(loc = 0.0, scale = baseline_linelist[param + '_err'], size = len(baseline_linelist) )
                baseline_linelist[param + '_vary'] = False
                baseline_linelist.loc[baseline_linelist[param + '_err']!=0, param + '_vary'] = True

        if self.directory == None:
            parameter_linelist.to_csv(self.param_linelist_file + '_MC.csv')
            baseline_linelist.to_csv(self.base_linelist_file + '_MC.csv')
        else:
            baseline_linelist.to_csv(os.path.join(self.directory,self.base_linelist_file + '_MC.csv'), index = False)
            parameter_linelist.to_csv(os.path.join(self.directory,self.param_linelist_file + '_MC.csv'))
    
    def add_iteration_results(self, fit_data_obj, iteration):

        lineparam_list = fit_data_obj.lineparam_list.copy()
        baseline_list = fit_data_obj.baseline_list.copy()

        lineparam_list['Iteration'] = iteration
        baseline_list['Iteration'] = iteration

        for spectrum in fit_data_obj.dataset.spectra:
            baseline_list.loc[ baseline_list['Spectrum Number'] == spectrum.spectrum_number, 'QF'] = spectrum.calculate_QF()

        self.iterations_param_linelist = pd.concat([self.iterations_param_linelist, lineparam_list], ignore_index=True)
        self.iterations_base_linelist = pd.concat([self.iterations_base_linelist, baseline_list], ignore_index=True)

    def calculate_kappa(self):
        param_drops = [c for c in self.param_linelist.columns if c.endswith('_kappa') or c.endswith('_experr')]
        if param_drops:
            self.param_linelist = self.param_linelist.drop(columns=param_drops)
            
        base_drops = [c for c in self.base_linelist.columns if c.endswith('_kappa') or c.endswith('_experr')]
        if base_drops:
            self.base_linelist = self.base_linelist.drop(columns=base_drops)
        
        param_names = [col.replace("_vary", "") for col in self.iterations_param_linelist.columns if "_vary" in col]
        new_param_cols = {f"{p}_kappa": 1.0 for p in param_names}
        df_new_param = pd.DataFrame(new_param_cols, index=self.param_linelist.index)

        # Calculate Kappa values
        for transition in self.iterations_param_linelist['Unnamed: 0'].unique():
            select = self.iterations_param_linelist[self.iterations_param_linelist['Unnamed: 0'] == transition]
            for param_name in param_names:
                vary_col = param_name + '_vary'
                if vary_col in select.columns and select[vary_col].values[0]:
                    kappa = select[param_name].std() / select[param_name + '_err'].mean()
                    if kappa < 1.0 or np.isnan(kappa):
                        kappa = 1.0
                    
                    # Update our pre-allocated structures in memory
                    df_new_param.loc[self.param_linelist['Unnamed: 0'] == transition, param_name + '_kappa'] = kappa

        self.param_linelist = pd.concat([self.param_linelist, df_new_param], axis=1)
        new_experr_cols = {}
        for param in param_names:
            new_experr_cols[param + '_experr'] = self.param_linelist[param + '_kappa'].values * self.param_linelist[param + '_err'].values
        self.param_linelist = pd.concat([self.param_linelist, pd.DataFrame(new_experr_cols, index=self.param_linelist.index)], axis=1)

        base_params = [col.replace("_vary", "") for col in self.iterations_base_linelist.columns if "_vary" in col]
        new_base_cols = {f"{p}_kappa": 1.0 for p in base_params}
        df_new_base = pd.DataFrame(new_base_cols, index=self.base_linelist.index)
            
        for spec_num in self.iterations_base_linelist['Spectrum Number'].unique():
            select = self.iterations_base_linelist[self.iterations_base_linelist['Spectrum Number'] == spec_num]
            
            for param_name in base_params:
                vary_col = param_name + '_vary'
                if vary_col in select.columns and select[vary_col].values[0]:
                    kappa = select[param_name].std() / select[param_name + '_err'].mean()
                    if kappa < 1.0 or np.isnan(kappa):
                        kappa = 1.0
                    df_new_base.loc[self.base_linelist['Spectrum Number'] == spec_num, param_name + '_kappa'] = kappa

        # Concatenate baseline kappas all at once
        self.base_linelist = pd.concat([self.base_linelist, df_new_base], axis=1)

        # Calculate baseline experimental errors all at once
        new_base_experr_cols = {}
        for param in base_params:
            new_base_experr_cols[param + '_experr'] = self.base_linelist[param + '_kappa'].values * self.base_linelist[param + '_err'].values
        self.base_linelist = pd.concat([self.base_linelist, pd.DataFrame(new_base_experr_cols, index=self.base_linelist.index)], axis=1)
    
    

    def calculate_ci(self, sigma = 1):

        alpha = sigma_to_alpha(sigma)
        lower_q = (alpha) /2
        upper_q = 1 - lower_q

        param_drops = [c for c in self.param_linelist.columns if c.endswith('_lower_ci') or c.endswith('_upper_ci')]
        if param_drops:
            self.param_linelist = self.param_linelist.drop(columns=param_drops)
                
        base_drops = [c for c in self.base_linelist.columns if  c.endswith('_lower_ci') or c.endswith('_upper_ci')]
        if base_drops:
            self.base_linelist = self.base_linelist.drop(columns=base_drops)
            
        param_names = [col.replace("_vary", "") for col in self.iterations_param_linelist.columns if "_vary" in col]
        new_param_cols_lower = {f"{p}_lower_ci": 0.0 for p in param_names}
        new_param_cols_upper = {f"{p}_upper_ci": 0.0 for p in param_names}
        new_param_cols = new_param_cols_lower | new_param_cols_upper
        df_new_param = pd.DataFrame(new_param_cols, index=self.param_linelist.index)

        for transition in self.iterations_param_linelist['Unnamed: 0'].unique():
            select = self.iterations_param_linelist[self.iterations_param_linelist['Unnamed: 0'] == transition]
            for param_name in param_names:
                
                vary_col = param_name + '_vary'
                if vary_col in select.columns and select[vary_col].values[0]:
                    lower_bound = select[param_name].quantile(lower_q)
                    upper_bound = select[param_name].quantile(upper_q)
                    exp_fit_value = self.param_linelist[self.param_linelist['Unnamed: 0'] == transition][param_name].values[0]
                    
                    
                    df_new_param.loc[self.param_linelist['Unnamed: 0'] == transition, param_name + '_lower_ci'] = lower_bound - exp_fit_value 
                    df_new_param.loc[self.param_linelist['Unnamed: 0'] == transition, param_name + '_upper_ci'] = upper_bound - exp_fit_value

        self.param_linelist = pd.concat([self.param_linelist, df_new_param], axis=1)
        
        base_params = [col.replace("_vary", "") for col in self.iterations_base_linelist.columns if "_vary" in col]
        
        new_base_cols_lower = {f"{p}_lower_ci": 0.0 for p in base_params}
        new_base_cols_upper = {f"{p}_upper_ci": 0.0 for p in base_params}
        new_base_cols = new_base_cols_lower | new_base_cols_upper
        df_new_base = pd.DataFrame(new_base_cols, index=self.base_linelist.index)

        for spec_num in self.iterations_base_linelist['Spectrum Number'].unique():
            select = self.iterations_base_linelist[self.iterations_base_linelist['Spectrum Number'] == spec_num]
            for param_name in base_params:
                vary_col = param_name + '_vary'
                if vary_col in select.columns and select[vary_col].values[0]:
                    lower_bound = select[param_name].quantile(lower_q)
                    upper_bound = select[param_name].quantile(upper_q)
                    exp_fit_value = self.base_linelist[self.base_linelist['Spectrum Number'] == spec_num][param_name].values[0]

                    df_new_base.loc[self.base_linelist['Spectrum Number'] == spec_num, param_name + '_lower_ci'] = lower_bound - exp_fit_value 
                    df_new_base.loc[self.base_linelist['Spectrum Number'] == spec_num, param_name + '_upper_ci'] = upper_bound - exp_fit_value
        self.base_linelist = pd.concat([self.base_linelist, df_new_base], axis=1)
    
    def parametric_bootstrap_analysis(self, kappa = True, CI = True, sigma = 1):

        if kappa:
            self.calculate_kappa()
        
        if CI:
            self.calculate_ci(sigma = sigma)

        if self.directory == None:

            self.base_linelist.to_csv(self.base_linelist_file + '_pbs.csv')
            self.param_linelist.to_csv(self.param_linelist_file + '_pbs.csv')         
        else:
            self.base_linelist.to_csv(os.path.join(self.directory, self.base_linelist_file + '_pbs.csv'))
            self.param_linelist.to_csv(os.path.join(self.directory, self.param_linelist_file + '_pbs.csv'))










        
        

        



            
                



