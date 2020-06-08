#Import Packages

from Utilities import *




class Dataset:
    """Combines spectrum objects into a Dataset object to enable multi-spectrum fitting.
    
    Parameters
    ----------
    spectra : list
        list of spectrum objects to be included in the Dataset. Example [spectrum_1, spectrum_2, . . .]
    dataset_name : str
        Used to provide a name for the Dataset to use when saving files
    param_linelist : pandas dataframe
        Reads in the  parameter linelist used in fitting.  This enables for consistency checks between the input spectra and the parameter line list. 
    baseline_order : int
        sets the baseline order for all spectra in the dataset.  This will automatically be set to the maximum baseline order across all spectrum included in the Dataset.
        
    """

    
    
    def __init__(self, spectra, dataset_name, param_linelist, baseline_order = 1, CIA_model = None):
        self.spectra = spectra
        self.dataset_name = dataset_name
        self.param_linelist = param_linelist
        self.baseline_order = baseline_order
        self.CIA_model = CIA_model
        self.renumber_spectra()
        self.molecule_list = self.correct_component_list()
        self.correct_etalon_list()
        self.max_baseline_order()
        self.broadener_list = self.get_broadener_list()
        
    def renumber_spectra(self):
        """renumbers the spectra to be sequential starting at 1 (called in the initialization of the class).
         """

        count = 1
        for spectrum in self.spectra:
            spectrum.set_spectrum_number(count)
            count+=1

    def max_baseline_order(self):
        """ sets the baseline order to be equal to the maximum in any of the included spectra.
        """
        
        baseline_order_list = []
        for spectrum in self.spectra:
            baseline_order_list.append(spectrum.baseline_order)
        self.baseline_order = max(baseline_order_list)

    def correct_component_list(self):
        """Corrects so that all spectra and the parameter line list share the same molecules, but the mole fraction is fixed to zero where molecules are not present (called in the initialization of the class). 
        """
        
        dataset_molecule_list = []
        for spectrum in self.spectra:
            dataset_molecule_list += (spectrum.molefraction.keys())
        molecules_in_paramlist = self.param_linelist['molec_id'].unique()
        for i in range(0, len(molecules_in_paramlist)):
            dataset_molecule_list.append(molecules_in_paramlist[i]) 
        dataset_molecule_list = list(set(dataset_molecule_list))
        
        for spectrum in self.spectra:
            spectrum_molefraction_dictionary = spectrum.get_molefraction()
            for molecule in dataset_molecule_list:
                if molecule not in spectrum_molefraction_dictionary:
                    spectrum_molefraction_dictionary[molecule] = 0
            spectrum.set_molefraction(spectrum_molefraction_dictionary)
        return dataset_molecule_list
    
    def get_broadener_list(self):
        """Provides a list of all broadeners in the dataset.
        

        Returns
        -------
        dataset_broadener_list : list
            list of all broadeners in the dataset.

        """
        dataset_broadener_list = []
        for spectrum in self.spectra:
            dataset_broadener_list += spectrum.Diluent.keys()
        dataset_broadener_list = list(set(dataset_broadener_list))
        return dataset_broadener_list
            


    def correct_etalon_list(self):
        """Corrects so that all spectrum share the same number of etalons, but the amplitude and period are fixed to zero where appropriate(called in the initialization of the class). 
        """

        dataset_etalon_list = []
        for spectrum in self.spectra:
            dataset_etalon_list += spectrum.etalons.keys()
        dataset_etalon_list = list(set(dataset_etalon_list))
        for spectrum in self.spectra:
            spectrum_etalon_dictionary = spectrum.get_etalons()
            for etalon_number in dataset_etalon_list:
                if etalon_number not in spectrum_etalon_dictionary:
                    spectrum_etalon_dictionary[etalon_number] = [0,0]
            spectrum.set_etalons(spectrum_etalon_dictionary)


    def get_etalons(self):
        """ Get list of number of etalons for spectra.
        

        Returns
        -------
        dataset_etalon_list : list
            etalon keys across spectra

        """
        dataset_etalon_list = []
        for spectrum in self.spectra:
            dataset_etalon_list += spectrum.etalons.keys()
        dataset_etalon_list = list(set(dataset_etalon_list))
        return dataset_etalon_list

    def get_molecules(self):
        """ Get list of molecules in spectra.
        

        Returns
        -------
        dataset_molecule_list : list
            list of molecules in spectra

        """

        dataset_molecule_list = []
        for spectrum in self.spectra:
            dataset_molecule_list += (spectrum.molefraction.keys())
        dataset_molecule_list = list(set(dataset_molecule_list)) 
        return dataset_molecule_list
       
    def get_spectra(self):
        return list(self.spectra)
    def get_dataset_name(self):
        return self.dataset_name
    def get_baseline_order(self):
        return self.baseline_order
    def set_dataset_name(self, new_dataset_name):
        self.dataset_name = new_dataset_name
    def set_baseline_order(self, new_baseline_order):
        self.baseline_order = new_baseline_order
    def set_spectra(self, new_spectra):
        self.spectra = new_spectra
        self.renumber_spectra()
    def get_number_spectra(self):
        return len(self.spectra)

    def get_spectrum_filename(self, spectrum_num):
        """ Gets spectrum filename for spectrum in Dataset.
        Parameters
        ----------
        spectrum_num : int
            Integer assigned to a given spectrum.

        Returns
        -------
        filename : str
            if spectrum_num in Dataset then the filename for that the spectrum_number is returned

        """

        for spectrum in self.spectra:
            if spectrum.spectrum_number == spectrum_num:
                return (spectrum.get_filename())
        return None

    def get_spectrum_pressure(self, spectrum_num):
        """Gets spectrum pressure for spectrum in Dataset.
        

        Parameters
        ----------
        spectrum_num : int
            Integer assigned to a given spectrum.

        Returns
        -------
        pressure : float
            if spectrum_num in Dataset then the pressure (torr) for that the spectrum_number is returned.

        """
        
        for spectrum in self.spectra:
            if spectrum.spectrum_number == spectrum_num:
                return (spectrum.get_pressure_torr())
        return None

    def get_spectrum_temperature(self, spectrum_num):
        """Gets spectrum temperature for spectrum in Dataset.
        

        Parameters
        ----------
        spectrum_num : int
            Integer assigned to a given spectrum.

        Returns
        -------
        temperature : float
            if spectrum_num in Dataset then the temperature (K) for that the spectrum_number is returned.

        """

        for spectrum in self.spectra:
            if spectrum.spectrum_number == spectrum_num:
                return (spectrum.get_temperature())
        return None

    def get_spectra_extremes(self):
        for spectrum in self.spectra:
            if spectrum.get_spectrum_number() == 1:
                wave_min = np.min(spectrum.wavenumber)
                wave_max = np.max(spectrum.wavenumber)
            else:
                if np.min(spectrum.wavenumber) < wave_min:
                    wave_min = np.min(spectrum.wavenumber)
                if np.max(spectrum.wavenumber) > wave_max:
                    wave_max = np.max(spectrum.wavenumber)
        return wave_min, wave_max

    def get_spectrum_extremes(self):
        """Gets the minimum and maximum wavenumber for the entire Dataset.
        

        Returns
        -------
        wave_min : float
            The minimum wavenumber in all spectra in the Dataset.
        wave_max : float
            The maximum wavenumber in all spectra in the Dataset.

        """

        extreme_dictionary = {}
        for spectrum in self.spectra:
            extreme_dictionary[spectrum.get_spectrum_number()] = [np.min(spectrum.wavenumber), np.max(spectrum.wavenumber)]
        return extreme_dictionary
            
    def get_number_nominal_temperatures(self):
        """ Get the number of nominal temperatures in the .
        

        Returns
        -------
        num_nominal_temperatures : int
            number of nominal temperatures in the Dataset
        nominal_temperatures : list
            list of all nominal temperatures listed in the input spectra

        """

        nominal_temperatures = []
        for spectrum in self.spectra:
            if spectrum.nominal_temperature not in nominal_temperatures:
                nominal_temperatures.append(spectrum.nominal_temperature)
        return len(nominal_temperatures), nominal_temperatures
         
    def average_QF(self):
        """Calculates the Average QF from all spectra.
        

        Returns
        -------
        average_QF : float
            average QF of all spectra in dataset

        """

        sum_ = 0
        for spectrum in self.spectra:
            sum_ += spectrum.calculate_QF()
        return sum_ / self.get_number_spectra()

    def get_list_spectrum_numbers(self):
        """ Generates a list of all spectrum_numbers.
        

        Returns
        -------
        spec_num_list : list
            list of all spectrum numbers used in the dataset.  

        """        
        spec_num_list = []
        for spectrum in self.spectra:
            spec_num_list.append(spectrum.spectrum_number)
        return spec_num_list    
        
    def generate_baseline_paramlist(self):
        """Generates a csv file called dataset_name + _baseline_paramlist, which will be used to generate another csv file that is used for fitting spectrum dependent parameters with columns for 
        spectrum number, segment number, x_shift, concentration for each molecule in the dataset, baseline terms (a = 0th order term, b = 1st order, etc), and etalon terms (set an amplitude, period, and phase for the number of etalons listed for each spectrum in the Dataset).

        Returns
        -------
        baseline_paramlist : dataframe
            dataframe containing information describing baseline parameters.  Also saves to .csv.  Either file can be edited before making the baseline parameter list used for fitting.  If editting the .csv file will need to regenerate dataframe from .csv.

        """
        
        baseline_paramlist = pd.DataFrame()
        for spectrum in self.spectra:
            for segment in list(set(spectrum.segments)):
                line = {}
                line['Spectrum Number'] = spectrum.spectrum_number
                line['Segment Number'] = segment
                line['Baseline Order'] = spectrum.baseline_order
                line['Pressure'] = spectrum.get_pressure()
                line['Temperature'] = spectrum.get_temperature()
                line['x_shift'] = spectrum.x_shift
                for molecule in spectrum.molefraction:
                    line['molefraction_' + (ISO[(molecule, 1)][4])] = (spectrum.molefraction[molecule])
                for i in range(0, self.baseline_order + 1):
                    if chr(i+97) == 'a':
                        line['baseline_' + chr(i+97)] = spectrum.alpha[0]
                    else:
                        line['baseline_' + chr(i+97)] = 0
                for etalon_name in spectrum.etalons:
                    line['etalon_' + str(etalon_name) + '_amp'] = spectrum.etalons[etalon_name][0]
                    line['etalon_' + str(etalon_name) + '_period'] = spectrum.etalons[etalon_name][1]
                    line['etalon_' + str(etalon_name) + '_phase'] = 0
                baseline_paramlist  = baseline_paramlist.append(line, ignore_index=True)
        baseline_paramlist = baseline_paramlist.set_index('Spectrum Number')
        baseline_paramlist.to_csv(self.dataset_name + '_baseline_paramlist.csv', index = True)
        return baseline_paramlist
    def generate_CIA_paramlist(self):
        """Generates a csv file called dataset_name + _CIA_paramlist, which will be used to generate another csv file that is used for fitting the broadband CIA that is common across all spectra, where the columns will be dependent on the CIA model used.  
        

        Returns
        -------
        CIA_paramlist : pandas dataframe
            dataframe containing information decribing the CIA parameters based on the CIA model chosen.  This dataframe is also saved to a dataframe.  Either file can be edited before making the CIA parameter list used for fitting.  If editting the .csv file will need to regenerate dataframe from .csv.

        """
        if self.CIA_model == None:
            return None
        elif self.CIA_model == 'Karman':
            CIA_paramlist = pd.DataFrame()
            CIA_list = []
            for molecule in self.molecule_list:
                molecule_name = ISO[(molecule, 1)][4]
                for broadener in self.broadener_list:
                    if broadener == 'self':
                        if molecule_name + '_' + molecule_name not in CIA_list:
                            CIA_list.append(molecule_name + '_' + molecule_name)
                    else:
                        if (molecule_name + '_' + broadener not in CIA_list) & (broadener + '_' + molecule_name not in CIA_list):
                            CIA_list.append(molecule_name + '_' + broadener)
                        
            CIA_paramlist['CIA Pair'] = CIA_list
            CIA_paramlist['EXCH_scalar'] = [0]*len(CIA_list)
            CIA_paramlist['EXCH_gamma'] = [3]*len(CIA_list)
            CIA_paramlist['EXCH_l'] = [2]*len(CIA_list)
            CIA_paramlist['SO_scalar'] = [0]*len(CIA_list)
            CIA_paramlist['SO_ahard'] = [7]*len(CIA_list)
            CIA_paramlist['SO_l'] = [2]*len(CIA_list)
            CIA_paramlist['bandcenter'] = [13122]*len(CIA_list)
            CIA_paramlist['Nmax'] = [31]*len(CIA_list)
            #index definition?            
            CIA_paramlist.to_csv(self.dataset_name + '_CIA_paramlist.csv', index = False) 
            return CIA_paramlist
        else:
            return None

    def generate_summary_file(self, save_file = False):
        """ Generates a summary file combining spectral information from all spectra in the Dataset.
        

        Parameters
        ----------
        save_file : bool, optional
            If True, then a .csv is generated in addition to the dataframe. The default is False.

        Returns
        -------
        summary_file : dataframe
            Summary dataframe comprised of spectral information inculding model and residuals for all spectra in Dataset.

        """

        summary_file = pd.DataFrame()
        for spectrum in self.spectra:
            spectrum_data = spectrum.save_spectrum_info(save_file = False)
            summary_file = summary_file.append(spectrum_data)
        if save_file:
            summary_file.to_csv(self.dataset_name + '.csv', index = False)
        return summary_file

    def plot_model_residuals(self):
        """ Generates a plot showing both the model and experimental data as a function of wavenumber in the main plot with a subplot showing the residuals as function of wavenumber.
        """

        fig = plt.figure(figsize = (16,10))
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
        ax0 = plt.subplot(gs[0])
        ax0.set_ylabel('$\\alpha (\\frac{ppm}{cm})$')
        ax1 = plt.subplot(gs[1])
        ax1.set_xlabel('Wavenumbers ($cm^{-1}$)')
        ax1.set_ylabel('Residuals $(\\frac{ppm}{cm})$')
        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
        ax0.ticklabel_format(useOffset=False)
        ax1.ticklabel_format(useOffset=False)
        for spectrum in self.spectra:
            plot_color = colors[((spectrum.spectrum_number % 7)-1)]
            ax0.plot(spectrum.wavenumber,spectrum.model, plot_color+'-')
            ax0.plot(spectrum.wavenumber, spectrum.alpha, plot_color+'.', label = spectrum.filename)
            ax1.plot(spectrum.wavenumber,spectrum.residuals, plot_color+"-")
        ax0.legend(bbox_to_anchor=(1, 1))
        plt.show()          