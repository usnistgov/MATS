import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import re

from .spectrum import simulate_spectrum
from .dataset import Dataset

from .utilities import (
    isotope_list_molecules_isotopes,
)

from .hapi import ISO

class MonteCarlo:

    def __init__(self, dataset, base_linelist_file, param_linelist_file):
        #Assumes parameters have been saved as .csvs
        # Currently written for one dataset type
        # Currently assumes that you aren't using the segment feature or the ILS features
        # This also currently won't adjust CIA parameters



        self.dataset = dataset
        self.dataset_name= self.dataset.dataset_name

        self.base_linelist = pd.read_csv(base_linelist_file + '.csv')
        self.param_linelist = pd.read_csv(param_linelist_file + '.csv')

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

  
   
    def simulate_dataset(self, at_SNR = True,
                         wave_error = 0.0, temperature_err =  {'bias': 0, 'function': None, 'params': {}}, 
                         pressure_err = {'per_bias': 0, 'function': None, 'params': {}},
                         molefraction_err = {},
                         wing_method = 'wing_wavenumbers', wing_cutoff = 25, wing_wavenumbers = 25, sim_window = 5):
    


        fit_summary = pd.read_csv(self.dataset_name + '.csv') #assumes only one dataset type

        list_spec = []
 
        for spectrum in self.spectra:
            spec_num = spectrum.spectrum_number
            spec_fit_summary = fit_summary[fit_summary['Spectrum Number']==spec_num]

            if at_SNR:
                SNR = spec_fit_summary['QF'].mean()
            else:
                SNR = None

            spec_base_linelist = self.base_linelist[self.base_linelist['Spectrum Number']==spec_num]
            baseline_params = self._extract_baseline_params(spec_base_linelist)
            if spectrum.dataspace == 'alpha':
                pathlength = 0
            else:
                pathlength = spec_base_linelist['pathlength'].values[0]
            

            molefraction_dict = self._extract_molefraction_parms(spec_base_linelist, spectrum)
            isotopic_abund_dict = self._extract_isoabund_parms(spec_base_linelist, spectrum)
            etalon_dict = self._extract_etalon_params(self, spec_base_linelist, spectrum)

            #generate simulated param_list for parameter_linelist this is where to pertub the parameters
            #Regenerate sw_based on sw_scale
            spec_param_linelist = self.param_linelist.copy()

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
                                         beta_formalism = spectrum.beta_formalism, TIPS = spectrum.TIPS, compressability_file= spectrum.compressability_file,
                                         BIA_model = spectrum.BIA_model, CIA_model = spectrum.CIA_model,
                                         
                                         wing_cutoff = wing_cutoff, wing_wavenumbers = wing_wavenumbers, wing_method = wing_method, sim_window = sim_window, 
                                         filename = 'temp', isotope_list = self.dataset.isotope_list, natural_abundance=spectrum.natural_abundance,)

            list_spec.append(sim_spec)
        return Dataset(list_spec, self.dataset_name + 'MC', self.dataset.param_linelist, BIA_model=self.dataset.BIA_model)


