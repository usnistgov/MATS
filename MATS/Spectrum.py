#Import Packages
from Utilities import *
from Fit_Dataset import HTP_from_DF_select, HTP_wBeta_from_DF_select





class Spectrum:
    """Spectrum class provides all information describing experimental or simulated spectrum.
    
    Parameters
    ----------
    filename : str
        file containing spectrum with .csv extension. File extension is not included in the name
    molefraction : dict
        mole fraction of each molecule in spectra in the format {molec_id: mole fraction (out of 1), molec_id: molefraction, . . . }
    natural_abundance : bool, optional
        flag for if the molecular species in the spectrum are at natural abundance
    abundance_ratio_MI : dict, optional
        if not at natural abundance sets the enhancement factor for each molecule and isotope in the following format {molec_id:{iso_id: enhancement, iso_id: enhancement}, . . . }
    isotope_list : dict, optional
        provides opportunity to specify the isotope look-up table.  Default is ISO, which is from HAPI.  If not using ISO, then must use this format and suggested you use function to add to ISO
    diluent : str, optional
        sets the diluent for the sample.  Default = 'air'
    Diluent : dict, optional
        sets the diluent for the sample if there are a combination of several diluents. Format {'he': 0.5, 'air': 0.5). NOTE: the line parameter file must have parameters that correspond to the diluent (ie gamma0_he, and gamma0_air). Additionally, the contribution from all diluents must sum to 1.
    spectrum_number : int, optional
        sets a number for the spectrum that will correspond to fit parameters. This is set in the Dataset class, so should not need to be defined manually
    input_freq : bool, optional
        True indicates that the frequency_column is in MHz. False indicates that the frequency column is in wavenumbers.
    input_tau : bool, optional
        True indicates that the tau_column is in us. False indicates that the tau column is in 1/c*tau (ppm/cm).
    pressure_column : str, optional
        Name of the pressure column in input file. Column should be in torr. Default is Cavity Pressure /Torr.
    temperature_column : str, optional
        Name of the temperature column in input file. Column should be in celsius. Default is Cavity Temperature Side 2 /C.
    frequency_column  : str, optional
        Name of the frequency column in input file. Column should be in MHz or Wavenumbers (with appropriate flag for input_freq). NOTE: This is not a detuning axis. Default is Total Frequency /MHz.
    tau_column : str, optional
        Name of the tau column in input file. Column should be in us or in 1/ctau (ppm/cm) (with appropriate flag for input_tau). Default is Mean tau/us.
    tau_stats_column  : str, optional
        Name of column containing the pt by pt error as a percent of the y-axis. Default is 'tau rel. std. dev./%
    segment_column  : str, optional
        Name of column with segment numbers in input file. This column allows spectrum background parameters to be treated in segments across the spectrum. If None then the entire spectrum will share the same segment.
    etalons : dict, optional
        Allows you to define etalons by an amplitude and period (1/ frequency in cm-1). Default is no etalons. Input is dictionary with keys being a number and the value being an array with the first index being the amplitude and the second being the period.
    nominal_temperature : int, optional
        nominal temperature indicates the approximate temperature of the spectrum. When combining spectra into a dataset this acts as a flag to whether temperature dependence terms should be parameters that can be fit or whether they should act as constants.  Default = 296
    x_shift : float, optional
        value in wavenumbers of the x shift for the spectrum axis. This is a fittable parameter. Be careful in using this parameter as floating multiple parameters with similar effects cause fits to not converge (ie. unconstrained line centers + x_shift or fits with line center, shifts, and x_shifts terms without enough lines per spectrum to isolate the different effects). Floating this term works best if you have a good initial guess. Default = 0
    baseline_order : int, optional
        sets the baseline order for this spectrum. Allows you to change the baseline order across the broader dataset.()
    weight : float, optional,
        set the weighting to used for the spectrum if using weighted fits.  This can be used to weight a whole spectrum in addition to the pt by pt weighting using the stats column.  
    ILS_function: string, optional
        Default is None and means that no instrument line shape is used in the fitting.  
        Function can be: SLIT_MICHELSON, SLIT_DIFFRACTION, SLIT_COSINUS, SLIT_DISPERSION, SLIT_GAUSSIAN, SLIT_TRIANGULAR, SLIT_RECTANGULAR corresponding to the ILS functions defined in HAPI or a user defined function.
    ILS_resolution: float/array, optional   
        Resolution is a float or array of ILS resolutions in wavenumbers.  The SlitFunctions defined in HAPI have 1 resolution, but this opens the option for the user defined function to be comprised of several functions with varying resolutions.  Default is 0.1 cm-1.
    ILS_wing: float/array, optional
        AF_wing is the a float or array consisting of the range the ILS is calculted over in cm-1.  This is a single value fort he HAPI slit functions, but could be an array of multiple values for user-defined functions.  Default is 10 cm-1 
    
    """


    
    def __init__(self, filename, molefraction = {}, natural_abundance = True, isotope_list = ISO, diluent = 'air', Diluent = {}, abundance_ratio_MI = {}, spectrum_number = 1, 
                    input_freq = True, input_tau = True, 
                    pressure_column = 'Cavity Pressure /Torr', temperature_column = 'Cavity Temperature Side 2 /C', frequency_column = 'Total Frequency /MHz', 
                    tau_column = 'Mean tau/us', tau_stats_column = None, segment_column = None, 
                    etalons = {}, nominal_temperature = 296, x_shift = 0, baseline_order = 1, weight = 1, 
                    ILS_function = None, ILS_resolution = 0.1, ILS_wing = 10):
        self.filename = filename
        self.molefraction = molefraction
        self.natural_abundance = natural_abundance
        self.abundance_ratio_MI = abundance_ratio_MI
        self.isotope_list = isotope_list
        self.diluent = diluent
        if Diluent == {}: #if Diluent was not set as the dictionary of various broadeners, then define dictionary with all of the broadening contribution coming from the diluent broadener
            if self.diluent == 'air':
                self.Diluent = {self.diluent: {'composition':1, 'm': 28.95734}}
            elif self.diluent == 'self': 
                self.Diluent = {self.diluent: {'composition':1, 'm': 0}}
                #mass will be set during HTP_wBeta_from_DF_select if necessary 
                print ("If using the self broadening term, then consider explicitly labeling the broadener (ie in an oxygen spectra use 'O2' instead of self).  This may avoid confusion in multiple species fits. ")
            else:
                print ('If using the HTP_wBeta_from_DF_select then you need to go back and use the Diluent{diluent:{"composition": 1, "m": mass}} format')
                self.Diluent = {self.diluent: {'composition':1, 'm': 0}}
                
                
        else:
            self.Diluent = Diluent
            if 'self' in self.Diluent:
                print ("You are using the 'self' term, consider explicitly labeling the broadener (ie in an oxygen spectra use 'O2' instead of 'self').  This may avoid confusion in multiple species fits. For single species fits it should not matter.")  
                print ("Double check that you did not include the equivalent of the self term explicitly (ie in an oxygen spectra having both 'O2' and 'self').")
        self.spectrum_number = spectrum_number
        self.pressure_column = pressure_column
        self.temperature_column = temperature_column
        self.frequency_column = frequency_column
        self.tau_column = tau_column
        self.tau_stats_column = tau_stats_column
        self.segment_column = segment_column
        self.input_freq = input_freq
        self.input_tau = input_tau
        self.etalons = etalons
        self.nominal_temperature = nominal_temperature
        self.x_shift = x_shift
        self.baseline_order = baseline_order
        self.weight = weight
        if self.weight ==0:
            print ('Change the weight to a non-zero value or remove the spectrum from the dataset.  If the weight is 0, then the residuals returned for that spectrum will be 0. ')
        self.ILS_function = ILS_function
        self.ILS_resolution = ILS_resolution
        self.ILS_wing = ILS_wing
        

        self.diluent_sum_check() # Makes sure that the diluent contributions sum to 1
        
        #Defined from contents of file
        file_contents = pd.read_csv(self.filename + '.csv',float_precision = 'high')
        self.pressure = file_contents[self.pressure_column].mean() / 760
        self.temperature = file_contents[self.temperature_column].mean() + 273.15
        if self.input_freq:
            self.frequency = file_contents[self.frequency_column].values
            self.wavenumber = self.frequency*10**6 / c
        else:
            self.wavenumber = file_contents[self.frequency_column].values
            self.frequency = self.wavenumber*c / 10**6
        if self.input_tau:
            self.tau = file_contents[self.tau_column].values
            self.alpha = (self.tau*c / 1e12)**-1
        else:
            self.alpha = file_contents[self.tau_column].values
            self.tau = (self.alpha*c / 1e12)**-1
            
        if self.tau_stats_column != None:
            stats = file_contents[self.tau_stats_column].values
            stats= np.nan_to_num(stats)
            median_tau_stats = np.median(stats [stats  > 0])
            stats[stats <= 0] = median_tau_stats
            self.tau_stats = stats
        else:
            self.tau_stats  = np.asarray(len(file_contents)*[0])
        if self.segment_column != None:
            self.segments = file_contents[self.segment_column].values
        else:
            self.segments = len(file_contents)*[1] 
        self.model = len(self.alpha)*[0]
        self.residuals = self.alpha - self.model
        self.background = len(self.alpha)*[0]
        self.cia = len(self.alpha)*[0]
    
    def diluent_sum_check(self):
        """Checks that if multiple broadeners are used that the contributions sum to one.
               

        Returns
        -------
        str
            Warning if the diluents don't sum to one

        """

        diluent_sum = 0
        for dil in self.Diluent:
            diluent_sum+=self.Diluent[dil]['composition']
        if diluent_sum != 1:
            print ("YOUR DILUENTS DO NOT SUM TO ONE!  They sum to " + str(diluent_sum))

    def segment_wave_alpha(self):
        """Defines the wavenumber, alpha, and indices of spectrum that correspond to a given spectrum segment.
        

        Returns
        -------
        wavenumber_segments : dict
            dictionary where the key corresponds to a segment number and the values correspond to the wavenumbers for that segment
        alpha_segments : dict
            dictionary where the key corresponds to a segment number and the values correspond to the alpha values for that segment.
        indices_segments : dict
            dictionary where the key corresponds to a segment number and the values correspond to the array indices for that segment.

        """
        
        wavenumber_segments = {}
        alpha_segments = {}
        indices_segments = {}
        for segment in list(set(self.segments)):
            indices = [i for i, x in enumerate(self.segments) if x == segment]
            indices_segments[segment] = indices
            wavenumber_segments[segment] = self.wavenumber[indices]
            alpha_segments[segment] = self.alpha[indices]
        return wavenumber_segments, alpha_segments, indices_segments
 
    ## GETTERS    
    def get_filename(self):
        return self.filename
    def get_molefraction(self):
        return self.molefraction
    def get_natural_abundance(self):
        return self.natural_abundance
    def get_abundance_ratio_MI(self):
        return self.abundance_ratio_MI
    def get_diluent(self):
        return self.diluent
    def get_Diluent(self):
        return self.Diluent
    def get_spectrum_number(self):
        return self.spectrum_number
    def get_pressure(self):
        return self.pressure
    def get_temperature(self):
        return self.temperature
    def get_pressure_torr(self):
        return self.pressure *760
    def get_temperature_C(self):
        return self.temperature - 273.15
    def get_frequency(self):
        return self.frequency
    def get_tau(self):
        return self.tau
    def get_tau_stats(self):
        return self.tau_stats
    def get_wavenumber(self):
        return self.wavenumber
    def get_alpha(self):
        return self.alpha
    def get_etalons(self):
        return self.etalons
    def get_model(self):
        return self.model
    def get_residuals(self):
        return self.residuals
    def get_background(self):
        return self.background
    def get_cia(self):
        return self.cia
    def get_nominal_temperature(self):
        return self.nominal_temperature

    ##SETTERS 
    def set_weight(self, new_weight):
        self.weight = new_weight
        if new_weight == 0:
            print ('Change the weight to a non-zero value or remove the spectrum from the dataset.  If the weight is 0, then the residuals returned for that spectrum will be 0. ')
    def set_molefraction(self, new_molefraction):
        self.molefraction = new_molefraction
    def set_natural_abundance(self, new_natural_abundance):
        self.natural_abundance = new_natural_abundance
    def set_abundance_ration_MI(self, new_abundance_ratio_MI):
        self.abundance_ratio_MI = new_abundance_ratio_MI
    def set_diluent(self, new_diluent):
        self.diluent = new_diluent
        if self.diluent == 'air':
            self.Diluent = {self.diluent: {'composition':1, 'm': 28.95734}}
        elif self.diluent == 'self': 
            self.Diluent = {self.diluent: {'composition':1, 'm': 0}}
                #mass will be set during HTP_wBeta_from_DF_select if necessary 
        else:
            print ('If using the HTP_wBeta_from_DF_select then you need to go back and use the Diluent{diluent:{"composition": 1, "m": mass}} format')
            self.Diluent = {self.diluent: {'composition':1, 'm': 0}}
    def set_Diluent(self, new_Diluent):
        self.Diluent = new_Diluent
    def set_spectrum_number(self, new_spectrum_number):
        self.spectrum_number = new_spectrum_number
    def set_pressure_column(self, new_pressure_column):
        self.pressure_column = new_pressure_column
        file_contents = pd.read_csv(self.filename + '.csv')
        self.pressure = file_contents[self.pressure_column].mean() / 760
    def set_temperature_column(self, new_temperature_column):
        self.temperature_colum = new_temperature_column
        file_contents = pd.read_csv(self.filename + '.csv')
        self.temperature = file_contents[self.temperature_column].mean() + 273.15
    def set_frequency_column(self, new_frequency_column):
        self.frequency_column = new_frequency_column
        file_contents = pd.read_csv(self.filename + '.csv')
        self.frequency = file_contents[self.frequency_column].values
        self.wavenumber = self.frequency*10**6 / c
    def set_tau_column(self, new_tau_column):
        self.tau_column = new_tau_column
        file_contents = pd.read_csv(self.filename + '.csv')
        self.tau = file_contents[self.tau_column].values
        self.alpha = (self.tau*c*1e-12)**-1
    def set_tau_stats_column(self, new_tau_stats_column):
        self.tau_stats_column = new_tau_stats_column
        file_contents = pd.read_csv(self.filename + '.csv')
        stats = file_contents[self.tau_stats_column].values
        stats= np.nan_to_num(stats)
        median_tau_stats = np.median(stats [stats  > 0])
        stats[stats <= 0] = median_tau_stats
        self.tau_stats = stats
        
    def set_etalons(self, new_etalons):
        self.etalons = new_etalons
    def set_model(self, new_model):
        self.model = new_model
    def set_residuals(self, new_residuals):
        self.residuals = new_residuals     
    def set_background(self, new_background):
        self.background = new_background
    def set_cia(self, new_cia):
        self.cia = new_cia
    def set_nominal_temperature(self, new_nominal_temperature):
        self.nominal_temperature = new_nominal_temperature 

    ##Other Functions
    def plot_freq_tau(self):
        """Generates a plot of tau (us) as a function of frequency (MHz).
        """
        
        plt.plot(self.frequency, self.tau)
        plt.xlabel('Frequency (MHz)')
        plt.ylabel('$\\tau (\mu s)$')
        plt.show()

    def plot_wave_alpha(self):
        """Generates a plot of alpha (ppm/cm) as a function of wavenumber (cm-1).
        """
        
        plt.plot(self.wavenumber, self.alpha)
        plt.xlabel('Wavenumber ($cm^{-1}$)')
        plt.ylabel('$\\alpha (\\frac{ppm}{cm})$')
        plt.show()

    def calculate_QF(self):
        """Calculates the quality of fit factor (QF) for a spectrum - QF = (maximum alpha - minimum alpha) / std(residuals).

        Returns
        -------
        float
            QF.

        """

        return np.around((self.alpha.max() - self.alpha.min()) / self.residuals.std(),0)

    def plot_model_residuals(self):
        """Generates a plot of the alpha and model (ppm/cm) as a function of wavenumber (cm-1) and on lower plot shows the residuals (ppm/cm) as a function of wavenumber (cm-1).
        """
        
        fig = plt.figure(figsize = (16,10))
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
        QF = self.calculate_QF()
        ax0 = plt.subplot(gs[0])
        ax0.plot(self.wavenumber,self.model, 'r-' )
        ax0.plot(self.wavenumber, self.alpha, 'k.')
        ax0.set_ylabel('$\\alpha (\\frac{ppm}{cm})$')
        ax0.ticklabel_format(useOffset=False)
        ax0.text(0.25, 0.95,'QF: ' + str(QF), horizontalalignment='center', verticalalignment='center', transform = ax0.transAxes)
        ax0.set_title(str(self.spectrum_number) +': ' + self.filename)
        ax1 = plt.subplot(gs[1])
        ax1.ticklabel_format(useOffset=False)
        ax1.plot(self.wavenumber,self.residuals, "r-")
        ax1.set_xlabel('Wavenumbers ($cm^{-1}$)')
        ax1.set_ylabel('Residuals $(\\frac{ppm}{cm})$')
        plt.show()

    def save_spectrum_info(self, save_file = False):
        """Saves spectrum information to a pandas dataframe with option to also save as as a csv file.
        

        Parameters
        ----------
        save_file : bool, optional
           If False, then only a dataframe is created. If True, then a csv file will be generated with the name filename + '_saved.csv'. The default is False.

        Returns
        -------
        new_file : dataframe
            returns a pandas dataframe with columns related to the spectrum information

        """

        file_contents = pd.read_csv(self.filename + '.csv')       
        new_file = pd.DataFrame()
        new_file['Spectrum Number'] = [self.spectrum_number]*len(self.alpha)
        new_file['Spectrum Name'] = [self.filename]*len(self.alpha)
        new_file['Frequency (MHz)'] = self.frequency
        new_file['Wavenumber (cm-1)'] = self.wavenumber
        new_file['Pressure (Torr)'] = file_contents[self.pressure_column].values
        new_file['Temperature (C)'] = file_contents[self.temperature_column].values
        new_file['Tau (us)'] = self.tau
        new_file['Tau Error (%)'] = self.tau_stats
        new_file['Alpha (ppm/cm)'] = self.alpha
        new_file['Model (ppm/cm)'] = self.model
        new_file['Residuals (ppm/cm)'] = self.residuals
        new_file['QF'] = [self.calculate_QF()]*len(new_file)
        new_file['Background'] = self.background
        new_file['CIA (ppm/cm)'] = self.cia 
        if save_file:
            new_file.to_csv(self.filename + '_saved.csv', index = False)
        return (new_file)

    def fft_spectrum(self):
        
        """Takes the FFT of the residuals of the spectrum, generates a plot of frequency (cm-1) versus amplitude (ppm/cm), and prints a dataframe with the 20 highest amplitude frequencies with the FFT frequency (period), amplitude, FFT phase, and frequency (cm-1).  
     
        """

        
        wave = self.wavenumber
        y = self.residuals
        wave_step = wave[1] - wave[0]
        A = np.fft.rfft(y)
        fft_freq = np.fft.rfftfreq(wave.shape[-1], wave_step)
        fft_amplitude = np.sqrt(A.real**2 + A.imag**2) / (len(A))
        fft_phase = np.arctan2(A.imag, A.real)
        FFT = pd.DataFrame()
        FFT['Period (cm)'] = fft_freq
        FFT['Amplitude'] = fft_amplitude
        FFT['Phase'] = fft_phase
        FFT['Freq (cm-1)'] = 1 / fft_freq
        fft_ = FFT.replace([np.inf, -np.inf], np.nan).dropna(how = 'any')
        fft_ = (fft_[fft_['Amplitude'] > 1e-5].sort_values(['Amplitude'], ascending = [0]).reset_index(drop = True))
        print (fft_.loc[0:20])
        plt.subplot(111)
        plt.plot(1 / fft_freq, fft_amplitude, '-')
        plt.ylabel('Amplitude')
        plt.xlabel('Experimental Frequency ($cm^{-1}$)')
        plt.ylabel('Amplitude (ppm/cm')
        plt.show()

def simulate_spectrum(parameter_linelist, 
                        wave_min=None, wave_max= None, wave_space=None, wavenumbers = [],  wave_error = 0, 
                        SNR = None, baseline_terms = [0], temperature = 25, temperature_err = {'bias': 0, 'function': None, 'params': {}}, pressure = 760, 
                        pressure_err = {'per_bias': 0, 'function': None, 'params': {}}, 
                        wing_cutoff = 25, wing_wavenumbers = 25, wing_method = 'wing_cutoff', filename = 'temp', molefraction = {}, molefraction_err = {},
                        isotope_list = ISO, natural_abundance = True, abundance_ratio_MI = {},diluent = 'air', Diluent = {}, 
                        nominal_temperature = 296, etalons = {}, x_shift = 0, IntensityThreshold = 1e-30, num_segments = 1, beta_formalism = False, 
                        ILS_function = None, ILS_resolution = 0.1, ILS_wing = 10):
    """Generates a synthetic spectrum, where the output is a spectrum object that can be used in MATS classes.
    

    Parameters
    ----------
    parameter_linelist : dataframe
        linelist following the convention of the linelists used for the HTP_from_DF_select.  Note that there will need to be a linemixing column for each nominal temperature, which you will have to do manually (ie y_air_296, y_self_296).
    wavenumbers : array of floats, optional
        array of wavenumbers for the simulation (cm-1).  If provided, then this axis will be used.  If wavenumbers = None, then the wave_min, wave_max, and wave_space will be used to calculate wavenumber grid.
    wave_min : float, optional
         minimum wavenumber for the simulation (cm-1)
    wave_max : float, optional
        maximum wavenumber for the simulation (cm-1).
    wave_space : float, optional
        wavenumber spacing for the simulation (cm-1).
    wave_error : float, optional
        absolute error on the wavenumber axis (cm-1) to include in simulations. The default is 0.
    SNR : float, optional
        Signal to noise ratio to impose on the simulated spectrum. The default is None.  If SNR is none there is no noise on the simulation.
    baseline_terms : list, optional
        polynomial baseline coefficients where the index is equal to the coefficient order, ie. [0, 1, 2] would correspond to baseline = 0 + 1*(wavenumber - minimum wavenumber) + 2*(wavenumber - minimum wavenumber)^2. The default is [0].
    temperature : float, optional
         temperature for simulation in celsius. The default is 25.
    temperature_err : dict, optional
        possible keys include 'bias', 'function', and 'params'. The bias indicates the absolute bias in Celsius of the temperature reading, which will be added to the input temperature. Function can be 'linear' with params 'm' and 'b' or 'sine' with parameters 'amp', 'phase', and 'phase'. These define a function that is added to both the bias and set temperature as a function of the wavenumber. Note: if 'function' key is not equal to None, then there also needs to be a params key to define the function.. The default is {'bias': 0, 'function': None, 'params': {}}.
    pressure : float, optional
        pressure for simulation in torr. The default is 760.
    pressure_err : dict, optional
        possible keys include bias, function, and params. The bias indicates the percent bias in of the pressure reading, which will be added to the input pressure. Function can be 'linear' with params 'm' and 'b' or 'sine' with parameters 'amp', 'phase', and 'phase'. These define a function that is added to both the bias and set pressure as a function of the wavenumber. Note: if 'function' key is not equal to None, then there also needs to be a params key to define the function.. The default is {'per_bias': 0, 'function': None, 'params': {}}.
    wing_cutoff : float, optional
        number of voigt half-widths to simulate on either side of each line. The default is 25.
    wing_wavenumbers : float, optional
        number of wavenumbers to simulate on either side of each line. The default is 25.
    wing_method : str, optional
        Provides choice between the wing_cutoff and wing_wavenumbers line cut-off options. The default is 'wing_cutoff'.
    filename : str, optional
        allows you to pick the output filename for the simulated spectra. The default is 'temp'.
    molefraction : dict, optional
        mole fraction of each molecule to be simulated in the spectrum in the format {molec_id: mole fraction (out of 1), molec_id: molefraction, . . . }. The default is {}.
    molefraction_err : dict, optional
        percent error in the mole fraction of each molecule to be simulated in the spectrum in the format {molec_id: percent error in mole fraction, molec_id: percent error in mole fraction, . . . }. The default is {}.
    natural_abundance : bool, optional
        flag for if the spectrum contains data at natural abundance. The default is True.
    abundance_ratio_MI : dict, optional
        if not at natural abundance sets the enhancement factor for each molecule and isotope in the following format {molec_id:{iso_id: enhancement, iso_id: enhancement}, . . . }. The default is {}.
    isotope_list : dict, optional
        provides opportunity to specify the isotope look-up table.  Default is ISO, which is from HAPI.  If not using ISO, then must use this format and suggested you use function to add to ISO
    diluent : str, optional
        sets the diluent for the sample if only using one broadener. The default is 'air'.
    Diluent : dict, optional
        sets the diluent for the sample if there are a combination of several. Format {'he': 0.5, 'air': 0.5). NOTE: the line parameter file must have parameters that correspond to the diluent (ie gamma0_he, and gamma0_air). Additionally, the contribution from all diluents must sum to 1.. The default is {}.
    nominal_temperature : int, optional
        nominal temperature indicates the approximate temperature of the spectrum. When combining spectra into a dataset this acts as a flag to whether temperature dependence terms should be parameters that can be fit or whether they should act as constants.. The default is 296.
    etalons : dict, optional
        Allows you to define etalons by an amplitude and period (1/ frequency in cm-1). Default is no etalons. Input is dictionary with keys being a number and the value being an array with the first index being the amplitude and the second being the period.. The default is {}.
    x_shift : float, optional
        value in wavenumbers of the x shift for the spectrum axis.. The default is 0.
    IntensityThreshold : float, optional
        minimum line intensity to use in the simulation. The default is 1e-30.
    num_segments : int, optional
        Number of segments in the file, which is implemented labeling the segment column into equal sequential se . The default is 10.
    beta_formalism : boolean, optional
        Indicates whether the beta correction for Dicke Narrowing should be used.  The default is False.
    ILS_function: string, optional
        Default is None and means that no instrument line shape is used in the fitting.  
        Function can be: SLIT_MICHELSON, SLIT_DIFFRACTION, SLIT_COSINUS, SLIT_DISPERSION, SLIT_GAUSSIAN, SLIT_TRIANGULAR, SLIT_RECTANGULAR corresponding to the ILS functions defined in HAPI or a user defined function.
    ILS_resolution: float/array, optional   
        Resolution is a float or array of ILS resolutions in wavenumbers.  The SlitFunctions defined in HAPI have 1 resolution, but this opens the option for the user defined function to be comprised of several functions with varying resolutions.  Default is 0.1 cm-1.
    ILS_wing: float, optional
        AF_wing is the a float consisting of the range the ILS is calculted over in cm-1. Default is 10 cm-1 
    
    Returns
    -------
    spectrum_file : .csv
        File that contains the simulated wavenumber axis, noisy wavenumber axis, absorbance data, noisy absorbance data, percent noise, pressure (torr), and temperature (C). The filename will correspond to the filename parameter, which has a default value of temp. The pressure and temperature columns will include whatever functional change there is to the pressure or temperature, but not the bias offset. This is coded to match how this error would manifest in experiments.
    spectrum_object : object
        Outputs a Spectrum class object. This makes it so the you can easily switch between reading in an experimental spectrum and simulated a synthetic spectrum by simply switching out whether the spectrum object is defined through the class definition or through the simulate_spectrum function.

    """

    #Checks to make a Diluent dictionary has been assigned    
    if not Diluent:
         if diluent == 'air':
             Diluent = {diluent: {'composition':1, 'm':28.95734}}
         elif diluent == 'self':
            Diluent = {diluent: {'composition':1, 'm':0}}
         else:
            Diluent = {diluent: {'composition':1, 'm':0}}
            print ('THIS IS GOING TO BREAK WITH A DIVISION ERROR IF YOU USE THE BETA VERSION')

    #Generates a linemixing column for each Diluent
    for dil in Diluent:
        parameter_linelist['y_' + dil] = parameter_linelist['y_' + dil + '_' + str(nominal_temperature)]
    #Set-Up Parameters
    baseline_terms = np.flip(baseline_terms)
    temperature_K = temperature + 273.15
    pressure_atm = pressure / 760
    if wavenumbers == []:
        wavenumbers = np.arange(wave_min, wave_max + wave_space, wave_space)

    wavenumbers_err = wavenumbers + wave_error*np.random.normal(loc = 0, scale =1, size = len(wavenumbers))
    #molefraction error
    #check that all moleules included in the parameter list are included in spectrum molefraction
    dataset_molecule_list = list(molefraction.keys())
    molecules_in_paramlist = parameter_linelist['molec_id'].unique()
    for i in range(0, len(molecules_in_paramlist)):
        dataset_molecule_list.append(molecules_in_paramlist[i]) 
    dataset_molecule_list = list(set(dataset_molecule_list))
    for molecule in dataset_molecule_list:
        if molecule not in molefraction:
            molefraction[molecule] = 0

    
    molefraction_w_error = {}        
    for species in molefraction:
        if molefraction_err == {}:
            molefraction_err[species] = 0
        elif species not in molefraction_err:
            molefraction_err[species] = 0
        molefraction_w_error[species] = molefraction[species] + molefraction[species]*(molefraction_err[species]/100)
    #pressure error
    pressure_w_error = pressure_atm + pressure_atm*(pressure_err['per_bias']/100)#adds the pressure bias based on the percent bias of the pressure measaurement.  Can loop this or set this to be constant for all spectra
    pressure_w_error = len(wavenumbers)*[pressure_w_error]   
    if pressure_err['function'] == 'linear':
        if 'params' in pressure_err:
            pressure_w_error += pressure_err['params']['m']*(wavenumbers-np.min(wavenumbers)) + pressure_err['params']['b']
    elif pressure_err['function'] == 'sine':
        if 'params' in pressure_err:
            pressure_w_error += etalon((wavenumbers-np.min(wavenumbers)), pressure_err['params']['amp'], pressure_err['params']['period'], pressure_err['params']['phase'])
    #temperature error  
    temperature_w_error = temperature_K + temperature_err['bias']
    temperature_w_error = len(wavenumbers)*[temperature_w_error] 
    if temperature_err['function'] == 'linear':
        if 'params' in pressure_err:
            temperature_w_error += temperature_err['params']['m']*(wavenumbers-np.min(wavenumbers)) + temperature_err['params']['b']
    elif pressure_err['function'] == 'sine':
        if 'params' in pressure_err:
            temperature_w_error += etalon((wavenumbers-np.min(wavenumbers)), temperature_err['params']['amp'], temperature_err['params']['period'], temperature_err['params']['phase'])

    #Define Segments
    seg_number = np.arange(len(wavenumbers))
    seg_number = np.abs(seg_number// (len(wavenumbers)/num_segments)).astype(int)

    alpha_array = len(wavenumbers)*[0]
    pressure_array = len(wavenumbers)*[0]
    temperature_array = len(wavenumbers)*[0]
    for seg in range(0, num_segments):

        segment_array = (np.where(seg_number == seg)[0])
        waves = np.take(wavenumbers, segment_array)
        segment_pressure = np.mean(np.take(pressure_w_error, segment_array))
        segment_temperature = np.mean(np.take(temperature_w_error, segment_array))

        if beta_formalism:
            waves, alpha = HTP_wBeta_from_DF_select(parameter_linelist,waves , wing_cutoff, wing_wavenumbers, wing_method,
                                p = segment_pressure, T = segment_temperature,  molefraction = molefraction_w_error, isotope_list = isotope_list,
                                natural_abundance = natural_abundance, abundance_ratio_MI = abundance_ratio_MI,  
                                Diluent = Diluent, diluent = diluent, IntensityThreshold = IntensityThreshold)
        else:
            waves, alpha = HTP_from_DF_select(parameter_linelist,waves , wing_cutoff, wing_wavenumbers, wing_method,
                    p = segment_pressure, T = segment_temperature,  molefraction = molefraction_w_error, isotope_list = isotope_list,
                    natural_abundance = natural_abundance, abundance_ratio_MI = abundance_ratio_MI,  
                    Diluent = Diluent, diluent = diluent, IntensityThreshold = IntensityThreshold)
        alpha_array[np.min(segment_array): np.max(segment_array)+1] = alpha * 1e6

        pressure_array[np.min(segment_array): np.max(segment_array)+1] = len(alpha)*[segment_pressure]
        temperature_array[np.min(segment_array): np.max(segment_array)+1] = len(alpha)*[segment_temperature]
    pressure_array = np.asarray(pressure_array) - pressure_atm*(pressure_err['per_bias'] / 100)
    temperature_array = np.asarray(temperature_array) - temperature_err['bias']
    
    #Calculate Baseline
    baseline = np.polyval(baseline_terms, wavenumbers -np.min(wavenumbers) )
    # Calculate Etalons
    etalon_model = len(wavenumbers)*[0]
    for r in range(1, len(etalons)+1):
        amp = etalons[r][0]
        period = etalons[r][1]
        phase = np.random.rand()
        x = wavenumbers - np.min(wavenumbers)
        etalon_model += amp*np.sin((2*np.pi * period)*x+ phase) 
    alpha_array += (baseline + etalon_model)
    if ILS_function != None:
        wavenumbers, alpha_array, i1, i2m, slit = convolveSpectrumSame(wavenumbers, alpha_array, SlitFunction = ILS_function, Resolution = ILS_resolution ,AF_wing=ILS_wing)  
    #Calculate Noisy Spectrum
    if SNR == None:
        alpha_noise = alpha_array
    else:   
        alpha_noise = alpha_array + np.max(alpha_array)*np.random.normal(loc = 0, scale =1, size = len(alpha_array))*1/SNR
 
    
    #Generate and save Simulated Spectrum File
    spectrum = pd.DataFrame()
    spectrum['Segment Number'] = seg_number
    spectrum['Wavenumber (cm-1)'] = wavenumbers
    spectrum['Wavenumber + Noise (cm-1)'] = wavenumbers_err
    spectrum['Alpha (ppm/cm)'] = alpha_array 
    spectrum['Alpha + Noise (ppm/cm)'] = alpha_noise
    spectrum['Noise (%)'] = 100 *(alpha_noise - (alpha_array + baseline + etalon_model)) / np.max(alpha_noise)
    spectrum['Pressure (Torr)'] = pressure_array*760
    spectrum['Temperature (C)'] = temperature_array - 273.15
    spectrum.to_csv(filename + '.csv', index = False)
    # Returns a spectrum class object for facile integration into the fitting workflow
    return Spectrum(filename, molefraction = molefraction, natural_abundance = natural_abundance, diluent = diluent, Diluent = Diluent, abundance_ratio_MI = abundance_ratio_MI, isotope_list = isotope_list,
                    spectrum_number = 1, input_freq = False, input_tau = False, 
                pressure_column = 'Pressure (Torr)', temperature_column = 'Temperature (C)', frequency_column = 'Wavenumber + Noise (cm-1)', 
                tau_column = 'Alpha + Noise (ppm/cm)', tau_stats_column = 'Noise (%)', segment_column = 'Segment Number',
                etalons = etalons, nominal_temperature = nominal_temperature, x_shift = x_shift, baseline_order = len(baseline_terms)-1, weight = 1, 
                ILS_function = ILS_function, ILS_resolution = ILS_resolution ,ILS_wing=ILS_wing)