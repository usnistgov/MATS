#Import Packages
import numpy as np
import pandas as pd
from bisect import bisect
import sys
import matplotlib.pyplot as plt
from matplotlib import gridspec
import re
import qgrid
from lmfit import Parameters, Minimizer
#import seaborn and set display properties
import seaborn as sns
sns.set_style("whitegrid")
sns.set_style("ticks")
sns.set_context("poster")
from hapi import EnvironmentDependency_Intensity, PYTIPS2017, molecularMass, pcqsdhc, ISO


def HTP_from_DF_select(linelist, waves, wing_cutoff = 50, wing_wavenumbers = 50, wing_method = 'wing_cutoff',
                p = 1, T = 296, molefraction = {}, 
                natural_abundance = True, abundance_ratio_MI = {},  Diluent = {}, diluent = 'air', IntensityThreshold = 1e-30):
    """Simulates the absorbance based on input line list, wavenumbers, and sample parameters.
    
    Outline
    
    1.  Uses provided wavenumber axis
    
    2.  Calculates the molecular density from pressure and temperature
    
    3.  Set-up Diluent dictionary if not given as input
    
    4.  Calculate line intensity and doppler width at temperature for all lines
    
    5.  Loop through each line in the line list and loop through each diluent and generate a line parameter 
        that is the appropriate ratio of each diluent species corrected for pressure and temperature.  For each line then simulate the line for the given simulation cutoffs and add to cross section
    
    6.  Return wavenumber and cross section arrays 
    

    Parameters
    ----------
    linelist : dataframe
        Pandas dataframe with the following column headers, where species corresponds to each diluent in the spectrum objects included in the dataset and nominal temperature corresponds to the nominal temperatures included in the dataset:
            nu = wavenumber of the spectral line transition (cm-1) in vacuum
            
            sw = The spectral line intensity (cm−1/(molecule⋅cm−2)) at Tref=296K
            
            elower = The lower-state energy of the transition (cm-1)
            
            molec_id = HITRAN molecular ID number
            
            local_iso_id = HITRAN isotopologue ID number
            
            gamma_0_species = half width at half maximum (HWHM) (cm−1/atm) at Tref=296K and reference pressure pref=1atm for a given diluent (air, self, CO2, etc.)
            
            n_gamma0_species = The coefficient of the temperature dependence of the half width
            
            delta_0_species = The pressure shift (cm−1/atm) at Tref=296K and pref=1atm of the line position with respect to the vacuum transition wavenumber νij
            
            n_delta0_species = the coefficient of the temperature dependence of the pressure shift
            
            SD_gamma_species = the ratio of the speed dependent width to the half-width at reference temperature and pressure
            
            n_gamma2_species = the coefficient of the temperature dependence of the speed dependent width NOTE: This is the temperature dependence of the speed dependent width not the ratio of the speed dependence to the half-width
            
            SD_delta_species = the ratio of the speed dependent shift to the collisional shift at reference temperature and pressure
            
            n_delta2_species = the coefficient of the temperature dependence of the speed dependent shift NOTE: This is the temperature dependence of the speed dependent shift not the ratio of the speed dependence to the shift
            
            nuVC_species = dicke narrowing term at reference temperature
            
            n_nuVC_species = coefficient of the temperature dependence of the dicke narrowing term
            
            eta_species = the correlation parameter for the VC and SD effects
            
            y_species_nominaltemperature = linemixing term (as currently written this doesn't have a temperature dependence, so read in a different column for each nominal temperature)
    waves : array
        1-D array comprised of wavenumbers (cm-1) to use as x-axis for simulation.
    wing_cutoff : float, optional
        number of half-widths for each line to be calculated for. The default is 50.
    wing_wavenumbers : float, optional
        number of wavenumber for each line to be calculated. The default is 50.
    wing_method : str, optional
        defines which wing cut-off option to use.  Options are wing_cutoff or wing_wavenumbers The default is 'wing_cutoff'.
    p : float, optional
        pressure for simulation in atmospheres. The default is 1.
    T : float, optional
        temperature for simulation in kelvin. The default is 296.
    molefraction : dict, optional
        takes the form {molecule_id : mole fraction, molecule_id: mole fraction, . . .}. The default is {}.
    natural_abundance : bool, optional
        True indicates the sample is at natural abundance. The default is True.
    abundance_ratio_MI : dictionary, optional
        If sample is not at natural abundance, then the natural abundance defines the enrichment factor compared to natural abundance(ie a sample where the main isotope is the only isotope would have a 1/natural abundance as the enrichment factor). Defined as {M:{I:enrichment factor, I: enrichment factor, I: enrichment factor}, . . . }. The default is {}.
    Diluent : dict, optional
        contains the species as the key with the value being the abundance of that diluent in the sample, ie {'He':0.5, 'self':0.5}. The default is {}.
    diluent : str, optional
        If Diluent = {}, then this value will be used to set the only diluent to be equal to this diluent. The default is 'air'.
    IntensityThreshold : float, optional
        minimum line intensity that will be simulated. The default is 1e-30.

    Returns
    -------
    wavenumbers : array
        wavenumber axis that should match the input waves
    xsect : array
        simulated cross section as a function of wavenumbers (ppm/cm)
        

    """
    #Generate X-axis for simulation
    wavenumbers = waves
        
    #Set Omegas to X-values
    Xsect = [0]*len(wavenumbers)
    
    #define reference temperature/pressure and calculate molecular density
    Tref = 296. # K
    pref = 1. # atm
    mol_dens = (p/9.869233e-7)/(1.380648813E-16*T) 
    
    #Sets-up the  Diluent (currently limited to air or self, unless manual input in Diluent)
    if not Diluent:
        Diluent = {diluent:1.}

    #Calculate line intensity
    linelist['SigmaT'] = 0
    linelist['SigmaTref'] = 0
    linelist['GammaD'] = 0
    linelist['m'] = 0
    linelist['abun_ratio'] = 1
    
    for molec in linelist['molec_id'].unique():
        for iso in linelist ['local_iso_id'].unique():
            linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'SigmaT'] = PYTIPS2017(molec,iso,T)
            linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'SigmaTref'] = PYTIPS2017(molec,iso,Tref)
            linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'm'] = molecularMass(molec,iso) * 1.66053873e-27 * 1000 #cmassmol and kg conversion
            if ( natural_abundance == False) and abundance_ratio_MI != {}:
                linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'abun_ratio'] = abundance_ratio_MI[molec][iso]
    linelist['LineIntensity'] = EnvironmentDependency_Intensity(linelist['sw'],T,Tref,linelist['SigmaT'],linelist['SigmaTref'],linelist['elower'],linelist['nu'])
    
    #Calculate Doppler Broadening
    linelist['GammaD'] = np.sqrt(2*1.380648813E-16*T*np.log(2)/linelist['m']/2.99792458e10**2)*linelist['nu']
    # Calculated Line Parameters across Broadeners
    linelist['Gamma0'] = 0
    linelist['Shift0'] = 0
    linelist['Gamma2'] = 0
    linelist['Shift2'] = 0
    linelist['NuVC'] = 0
    linelist['Eta'] = 0
    linelist['Y'] = 0
    for species in Diluent:
        abun = Diluent[species]
        #Gamma0: pressure broadening coefficient HWHM
        linelist['Gamma0'] += abun*(linelist['gamma0_%s'%species]*(p/pref)*((Tref/T)**linelist['n_gamma0_%s'%species]))
        #Delta0
        linelist['Shift0'] += abun*((linelist['delta0_%s'%species] + linelist['n_delta0_%s'%species]*(T-Tref))*p/pref)
        #Gamma2
        linelist['Gamma2'] += abun*(linelist['SD_gamma_%s'%species]*linelist['gamma0_%s'%species]*(p/pref)*((Tref/T)**linelist['n_gamma2_%s'%species]))
        #Delta2
        linelist['Shift2'] += abun*((linelist['SD_delta_%s'%species]*linelist['delta0_%s'%species] + linelist['n_delta2_%s'%species]*(T-Tref))*p/pref)
        #nuVC
        linelist['NuVC'] += abun*(linelist['nuVC_%s'%species]*(p/pref)*((Tref/T)**(linelist['n_nuVC_%s'%species])))
        #eta
        linelist['Eta'] += abun*linelist['eta_%s'%species] 
        #Line mixing
        linelist['Y'] += abun*(linelist['y_%s'%species]*(p/pref))
    
    #Line profile simulation cut-off determination    
    if wing_method == 'wing_cutoff':
        linelist['line cut-off'] = (0.5346*linelist['Gamma0'] + (0.2166*linelist['Gamma0']**2 + linelist['GammaD']**2)**0.5)*wing_cutoff
    else:
        linelist['line cut-off'] = wing_wavenumbers
    
    #Enforce  Line Intensity Simulation Threshold
    linelist = linelist[linelist['LineIntensity']>= IntensityThreshold] 
    
    #For each line calculate spectrum and add to the global spectrum
    for index, line in linelist.iterrows():
        BoundIndexLower = bisect(wavenumbers, line['nu'] - line['line cut-off'])
        BoundIndexUpper = bisect(wavenumbers, line['nu'] + line['line cut-off'])
        lineshape_vals_real, lineshape_vals_imag = pcqsdhc(line['nu'],line['GammaD'],line['Gamma0'],line['Gamma2'],line['Shift0'],line['Shift2'],
                                                    line['NuVC'],line['Eta'],wavenumbers[BoundIndexLower:BoundIndexUpper])    
        Xsect[BoundIndexLower:BoundIndexUpper] += mol_dens  * \
                                                    molefraction[line['molec_id']] * line['abun_ratio'] * \
                                                    line['LineIntensity'] * (lineshape_vals_real + line['Y']*lineshape_vals_imag) 
    
    # Return two arrays corresponding to the wavenumber axis and the calculated cross-section                                                                                                      
    return (wavenumbers, np.asarray(Xsect))  


class Spectrum:
    """Spectrum class provides all information describing experimental or simulated spectrum necessary for fitting.
    
    Parameters
    ----------
    filename : str
        file containing spectrum with .csv extension. File extension is not included in the name
    molefraction : dict
        mole fraction of each molecule in spectra in the format {molec_id: mole fraction (out of 1), molec_id: molefraction, . . . }
    natural_abundance : bool, optional
        flag for if the spectrum contains data at natural abundance
    abundance_ratio_MI : dict, optional
        if not at natural abundance sets the enhancement factor for each molecule and isotope in the following format {molec_id:{iso_id: enhancement, iso_id: enhancement}, . . . }
    diluent : str, optional
        sets the diluent for the sample.  Default = 'air'
    Diluent : dict, optional
        sets the diluent for the sample if there are a combination of several. Format {'he': 0.5, 'air': 0.5). NOTE: the line parameter file must have parameters that correspond to the diluent (ie gamma0_he, and gamma0_air). Additionally, the contribution from all diluents must sum to 1.
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
        Name of the tau stats column in input file. Default is 'tau rel. std. dev./%
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
    """

    def __init__(self, filename, molefraction = {}, natural_abundance = True, diluent = 'air', Diluent = {}, abundance_ratio_MI = {}, spectrum_number = 1, 
                    input_freq = True, input_tau = True, 
                    pressure_column = 'Cavity Pressure /Torr', temperature_column = 'Cavity Temperature Side 2 /C', frequency_column = 'Total Frequency /MHz', 
                    tau_column = 'Mean tau/us', tau_stats_column = None, segment_column = None, 
                    etalons = {}, nominal_temperature = 296, x_shift = 0, baseline_order = 1):
        self.filename = filename
        self.molefraction = molefraction
        self.natural_abundance = natural_abundance
        self.abundance_ratio_MI = abundance_ratio_MI
        self.diluent = diluent
        if Diluent == {}: #if Diluent was not set as the dictionary of various broadeners, then define dictionary with all of the broadening contribution coming from the diluent broadener
            self.Diluent = {self.diluent: 1}
        else:
            self.Diluent = Diluent
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
        self.diluent_sum_check() # Makes sure that the diluent contributions sum to 1
        
        #Defined from contents of file
        file_contents = pd.read_csv(self.filename + '.csv',float_precision = 'High')
        self.pressure = file_contents[self.pressure_column].mean() / 760
        self.temperature = file_contents[self.temperature_column].mean() + 273.15
        if self.input_freq:
            self.frequency = file_contents[self.frequency_column].values
            self.wavenumber = self.frequency*10**6 / 29979245800
        else:
            self.wavenumber = file_contents[self.frequency_column].values
            self.frequency = self.wavenumber*29979245800 / 10**6
        if self.input_tau:
            self.tau = file_contents[self.tau_column].values
            self.alpha = (self.tau*0.0299792458)**-1
        else:
            self.alpha = file_contents[self.tau_column].values
            self.tau = (self.alpha*0.0299792458)**-1
            
        if self.tau_stats_column != None:
            self.tau_stats = file_contents[self.tau_stats_column].values
        else:
            self.tau_stats = None
        if self.segment_column != None:
            self.segments = file_contents[self.segment_column].values
        else:
            self.segments = len(file_contents)*[1] 
        self.model = len(self.alpha)*[0]
        self.residuals = self.alpha - self.model
        self.background = len(self.alpha)*[0]
    
    def diluent_sum_check(self):
        """Checks that if multiple broadeners are used that the contributions sum to one
               

        Returns
        -------
        str
            Warning if the diluents don't sum to one

        """
        diluent_sum = 0
        for dil in self.Diluent:
            diluent_sum+=self.Diluent[dil]
        if diluent_sum != 1:
            print ("YOUR DILUENTS DO NOT SUM TO ONE!")

    def segment_wave_alpha(self):
        """Defines the wavenumber, alpha, and indices of spectrum that correspond to a given spectrum segment
        

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
    def get_nominal_temperature(self):
        return self.nominal_temperature

    ##SETTERS 
    def set_molefraction(self, new_molefraction):
        self.molefraction = new_molefraction
    def set_natural_abundance(self, new_natural_abundance):
        self.natural_abundance = new_natural_abundance
    def set_abundance_ration_MI(self, new_abundance_ratio_MI):
        self.abundance_ratio_MI = new_abundance_ratio_MI
    def set_diluent(self, new_diluent):
        self.diluent = new_diluent
        self.Diluent = {new_diluent : 1}
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
        self.wavenumber = self.frequency*10**6 / 29979245800
    def set_tau_column(self, new_tau_column):
        self.tau_column = new_tau_column
        file_contents = pd.read_csv(self.filename + '.csv')
        self.tau = file_contents[self.tau_column].values
        self.alpha = (self.tau*0.0299792458)**-1
    def set_tau_stats_column(self, new_tau_stats_column):
        self.tau_stats_column = new_tau_stats_column
        file_contents = pd.read_csv(self.filename + '.csv')
        self.tau_stats = file_contents[self.tau_stats_column].values
    def set_etalons(self, new_etalons):
        self.etalons = new_etalons
    def set_model(self, new_model):
        self.model = new_model
    def set_residuals(self, new_residuals):
        self.residuals = new_residuals     
    def set_background(self, new_background):
        self.background = new_background
    def set_nominal_temperature(self, new_nominal_temperature):
        self.nominal_temperature = new_nominal_temperature 

    ##Other Functions
    def plot_freq_tau(self):
        """Generates a plot of Tau (us) as a function of Frequency (MHz)
        """
        plt.plot(self.frequency, self.tau)
        plt.xlabel('Frequency (MHz)')
        plt.ylabel('$\\tau (\mu s)$')
        plt.show()

    def plot_wave_alpha(self):
        """Generates a plot of alpha (ppm/cm) as a function of Wavenumber (cm-1)
        """
        plt.plot(self.wavenumber, self.alpha)
        plt.xlabel('Wavenumber ($cm^{-1}$)')
        plt.ylabel('$\\alpha (\\frac{ppm}{cm})$')
        plt.show()

    def calculate_QF(self):
        """Calculates the quality of fit factor (QF) for spectrum
        
        QF = (maximum alpha - minimum alpha) / std(residuals)
        

        Returns
        -------
        float
            QF.

        """
        return np.around((self.alpha.max() - self.alpha.min()) / self.residuals.std(),0)

    def plot_model_residuals(self):
        """Generates a plot of the alpha and model (ppm/cm) as a function of Wavenumber (cm-1) and on lower plot shows the residuals (ppm/cm) as a function of wavenumber (cm-1)
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
        """Saves spectrum information to a pandas dataframe with option to also save as .csv
        

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
        if save_file:
            new_file.to_csv(self.filename + '_saved.csv', index = False)
        return (new_file)

    def fft_spectrum(self):
        """Takes the FFT of the residuals of the spectrum
        Prints a dataframe with the 20 highest amplitude frequencies with the FFT frequency (period), amplitude, FFT phase, and frequency (cm-1).
        Also generates a plot of frequency (cm-1) versus amplitude (ppm/cm)
        """
        wave = self.wavenumber
        y = self.residuals
        wave_step = wave[1] - wave[0]
        A = np.fft.rfft(y)
        fft_freq = np.fft.rfftfreq(wave.shape[-1], wave_step)
        fft_amplitude = np.sqrt(A.real**2 + A.imag**2) / (len(A))
        fft_phase = np.arctan2(A.imag, A.real)
        FFT = pd.DataFrame()
        FFT['Frequency'] = fft_freq
        FFT['Amplitude'] = fft_amplitude
        FFT['Phase'] = fft_phase
        FFT['Freq (cm-1)'] = 1 / fft_freq
        fft_ = FFT.replace([np.inf, -np.inf], np.nan).dropna(how = 'any')
        fft_ = (fft_[fft_['Amplitude'] > 1e-5].sort_values(['Amplitude'], ascending = [0]).reset_index(drop = True))
        print (fft_.loc[0:20])
        plt.subplot(111)
        plt.plot(1 / fft_freq, fft_amplitude, '-')
        plt.ylabel('Amplitude')
        plt.xlabel('Experimental Wavenumber ($cm^{-1}$)')
        plt.ylabel('Amplitude (ppm/cm')
        plt.show()

class Dataset:
    """Combines spectrum objects into a Dataset object to enable multi-spectrum fitting
    
    Parameters
    ----------
    spectra : list
        list of spectrum objects to be included in the Dataset. Example [spectrum_1, spectrum_2, . . .]
    dataset_name : str
        Used to provide a name for the Dataset to use when saving files
    baseline_order : int
        sets the baseline order for all spectra in the dataset.  This will automaticall be set to the maximum baseline order across all spectrum included in the Dataset.
        
    """
    def __init__(self, spectra, dataset_name, baseline_order = 1):
        self.spectra = spectra
        self.dataset_name = dataset_name
        self.baseline_order = baseline_order
        self.renumber_spectra()
        self.correct_component_list()
        self.correct_etalon_list()
        self.max_baseline_order()
        
    def renumber_spectra(self):
        """renumbers the spectra to be sequential starting at 1. It is called in the initialization of the class
         """
        count = 1
        for spectrum in self.spectra:
            spectrum.set_spectrum_number(count)
            count+=1

    def max_baseline_order(self):
        """ sets the baseline order to be equal to the maximum set in any of the included spectra       
        """
        baseline_order_list = []
        for spectrum in self.spectra:
            baseline_order_list.append(spectrum.baseline_order)
        self.baseline_order = max(baseline_order_list)

    def correct_component_list(self):
        """Corrects so that all spectrum share the same molecules, but the mole fraction is fixed to zero where molecule is not present.  Called at the initialization of the class
        """
        dataset_molecule_list = []
        for spectrum in self.spectra:
            dataset_molecule_list += (spectrum.molefraction.keys())
        dataset_molecule_list = list(set(dataset_molecule_list))
        for spectrum in self.spectra:
            spectrum_molefraction_dictionary = spectrum.get_molefraction()
            for molecule in dataset_molecule_list:
                if molecule not in spectrum_molefraction_dictionary:
                    spectrum_molefraction_dictionary[molecule] = 0
            spectrum.set_molefraction(spectrum_molefraction_dictionary)

    def correct_etalon_list(self):
        """Corrects so that all spectrum share the same number of etalons, but the amplitude and period are fixed to zero where appropriate.  Called at the initialization of the class
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
        """ Get list of number of etalons for spectra
        

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
        """ Get list of molecules in spectra
        

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
        """ Gets spectrum filename for spectrum in Dataset
        

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
        """Gets spectrum pressure for spectrum in Dataset
        

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
        """Gets spectrum temperature for spectrum in Dataset
        

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
        """Gets the minimum and maximum wavenumber for the entire Dataset
        

        Returns
        -------
        wave_min : float
            The minimum wavenumber in all spectra in the Dataset.
        wave_max : float
            The maximum wavenumber in all spectra in the Dataset.

        """
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
        """Gets the minimum and maximum wavenumber for each spectrum in the Dataset
        

        Returns
        -------
        extreme_dictionary : dict
            Dictionary where the key corresponds to the spectrum_num with the value equal to a list where list[0] = minimum wavenumber in spectrum and list[1] = maximum wavenumber in spectrum
        """
        extreme_dictionary = {}
        for spectrum in self.spectra:
            extreme_dictionary[spectrum.get_spectrum_number()] = [np.min(spectrum.wavenumber), np.max(spectrum.wavenumber)]
        return extreme_dictionary
            
    def get_number_nominal_temperatures(self):
        """ Get the number of nominal temperatures in the Dataset
        

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
        """Calculates the Average QF from all spectrum
        

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
        """ Generates a list of all spectrum_numbers
        

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
        """Generates a .csv file called dataset_name + _baseline_paramlist.csv. This file will be used to generate another .csv that is used for fitting these spectrum dependent parameters. With columns for 
        spectrum number, segment number, x_shift, concentration for each molecule in the dataset, baseline terms (a = 0th order term, b = 1st order , . . .), and etalon terms (set an amplitude, period, and phase for the number of etalons listed for each spectrum in the Dataset)

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
                    line['etalon_' + str(etalon_name) + '_freq'] = spectrum.etalons[etalon_name][1]
                    line['etalon_' + str(etalon_name) + '_phase'] = 0
                baseline_paramlist  = baseline_paramlist.append(line, ignore_index=True)
        baseline_paramlist = baseline_paramlist.set_index('Spectrum Number')
        baseline_paramlist.to_csv(self.dataset_name + '_baseline_paramlist.csv')
        return baseline_paramlist

    def generate_summary_file(self, save_file = False):
        """ Generates a summary file combining spectral information from all spectra in the Dataset
        

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
        """ Generates a plot showing both the model and experimental data as a function of wavenumber in the main plot. A subplot shows the residuals as function of wavenumber.  Each spectra in the Dataset will have a different color.  
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
           
def max_iter(pars, iter, resid, *args, **kws):
    """Used to stop fitting if the number of iterations is greater that 2500
    """
    if iter > 2500:
        return True
    else:
        return False
                   
def etalon(x, amp, freq, phase):
    """Etalon definition
    

    Parameters
    ----------
    x : array
        array of floats used to define the x-axis
    amp : float
        amplitude of the etalon.
    freq : float
        period of the etalon.
    phase : float
        phase of the etalon.

    Returns
    -------
    etalon : array
        etalon as a function of input x-axis, amplitude, period, and phase.

    """
    return amp*np.sin((2*np.pi * freq)*x+ phase)   

    
        
def simulate_spectrum(parameter_linelist, wave_min, wave_max, wave_space, wave_error = 0, 
                        SNR = None, baseline_terms = [0], temperature = 25, temperature_err = {'bias': 0, 'function': None, 'params': {}}, pressure = 760, 
                        pressure_err = {'per_bias': 0, 'function': None, 'params': {}}, 
                        wing_cutoff = 25, wing_wavenumbers = 25, wing_method = 'wing_cutoff', filename = 'temp', molefraction = {}, molefraction_err = {},
                        natural_abundance = True, abundance_ratio_MI = {},diluent = 'air', Diluent = {}, 
                        nominal_temperature = 296, etalons = {}, x_shift = 0, IntensityThreshold = 1e-30, num_segments = 10):
    """allows for a synthetic spectrum to be generated, where the output is a spectrum object that can be used in MATS classes.
    

    Parameters
    ----------
    parameter_linelist : dataframe
        linelist following the convention of the linelists used for the HTP_from_DF_select.  Note that there will need to be a linemixing column for each nominal temperature, which you will have to do manually (ie y_air_296, y_self_296).
    wave_min : float
         minimum wavenumber for the simulation (cm-1)
    wave_max : float
        maximum wavenumber for the simulation (cm-1).
    wave_space : float
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
        possible keys include 'bias', 'function', and 'params'. The bias indicates the absolute bias in Celsius of the temperature reading, which will be added to the input temperature. Function can be 'linear' with params 'm' and 'b' or 'sine' with parameters 'amp', 'freq', and 'phase'. These define a function that is added to both the bias and set temperature as a function of the wavenumber. Note: if 'function' key is not equal to None, then there also needs to be a params key to define the function.. The default is {'bias': 0, 'function': None, 'params': {}}.
    pressure : float, optional
        pressure for simulation in torr. The default is 760.
    pressure_err : dict, optional
        possible keys include bias, function, and params. The bias indicates the percent bias in of the pressure reading, which will be added to the input pressure. Function can be 'linear' with params 'm' and 'b' or 'sine' with parameters 'amp', 'freq', and 'phase'. These define a function that is added to both the bias and set pressure as a function of the wavenumber. Note: if 'function' key is not equal to None, then there also needs to be a params key to define the function.. The default is {'per_bias': 0, 'function': None, 'params': {}}.
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
    num_segments : TYPE, optional
        Number of segments in the file, which is implemented labeling the segment column into equal sequential se . The default is 10.

    Returns
    -------
    spectrum_file : .csv
        File that contains the simulated wavenumber axis, noisy wavenumber axis, absorbance data, noisy absorbance data, pressure (torr), and temperature (C). The filename will correspond to the filename parameter, which has a default value of temp. The pressure and temperature columns will include whatever functional change there is to the pressure or temperature, but not the bias offset. This is coded to match how this error would manifest in experiments.
    spectrum_object : object
        Outputs a Spectrum class object. This makes it so the you can easily switch between reading in an experimental spectrum and simulated a synthetic spectrum by simply switching out whether the spectrum object is defined through the class definition or through the simulate_spectrum function.

    """
    #Checks to make a Diluent dictionary has been assigned    
    if not Diluent:
        Diluent = {diluent:1.}
    #Generates a linemixing column for each Diluent
    for dil in Diluent:
        parameter_linelist['y_' + dil] = parameter_linelist['y_' + dil + '_' + str(nominal_temperature)]
    #Set-Up Parameters
    baseline_terms = np.flip(baseline_terms)
    temperature_K = temperature + 273.15
    pressure_atm = pressure / 760
    wavenumbers = np.arange(wave_min, wave_max + wave_space, wave_space)

    wavenumbers_err = wavenumbers + wave_error*np.random.normal(loc = 0, scale =1, size = len(wavenumbers))
    #molefraction error
    molefraction_w_error = {}        
    for species in molefraction:
        if molefraction_err == {}:
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
            pressure_w_error += etalon((wavenumbers-np.min(wavenumbers)), pressure_err['params']['amp'], pressure_err['params']['freq'], pressure_err['params']['phase'])
    #temperature error  
    temperature_w_error = temperature_K + temperature_err['bias']
    temperature_w_error = len(wavenumbers)*[temperature_w_error] 
    if temperature_err['function'] == 'linear':
        if 'params' in pressure_err:
            temperature_w_error += temperature_err['params']['m']*(wavenumbers-np.min(wavenumbers)) + temperature_err['params']['b']
    elif pressure_err['function'] == 'sine':
        if 'params' in pressure_err:
            temperature_w_error += etalon((wavenumbers-np.min(wavenumbers)), temperature_err['params']['amp'], temperature_err['params']['freq'], temperature_err['params']['phase'])

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

  
        waves, alpha = HTP_from_DF_select(parameter_linelist,waves , wing_cutoff, wing_wavenumbers, wing_method,
                            p = segment_pressure, T = segment_temperature,  molefraction = molefraction_w_error, 
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
        freq = etalons[r][1]
        phase = np.random.rand()
        x = wavenumbers - np.min(wavenumbers)
        etalon_model += amp*np.sin((2*np.pi * freq)*x+ phase) 
    #Calculate Noisy Spectrum
    if SNR == None:
        alpha_noise = alpha_array
    else:   
        alpha_noise = alpha_array + np.max(alpha_array)*np.random.normal(loc = 0, scale =1, size = len(alpha_array))*1/SNR
    alpha_noise += (baseline + etalon_model) 
    #Generate and save Simulated Spectrum File
    spectrum = pd.DataFrame()
    spectrum['Segment Number'] = seg_number
    spectrum['Wavenumber (cm-1)'] = wavenumbers
    spectrum['Wavenumber + Noise (cm-1)'] = wavenumbers_err
    spectrum['Alpha (ppm/cm)'] = alpha_array + baseline + etalon_model
    spectrum['Alpha + Noise (ppm/cm)'] = alpha_noise
    spectrum['Pressure (Torr)'] = pressure_array*760
    spectrum['Temperature (C)'] = temperature_array - 273.15
    spectrum.to_csv(filename + '.csv', index = False)
    # Returns a spectrum class object for facile integration into the fitting workflow
    return Spectrum(filename, molefraction = molefraction, natural_abundance = natural_abundance, diluent = diluent, Diluent = Diluent, abundance_ratio_MI = abundance_ratio_MI, spectrum_number = 1, 
                input_freq = False, input_tau = False, 
                pressure_column = 'Pressure (Torr)', temperature_column = 'Temperature (C)', frequency_column = 'Wavenumber + Noise (cm-1)', 
                tau_column = 'Alpha + Noise (ppm/cm)', tau_stats_column = None, segment_column = 'Segment Number',
                etalons = etalons, nominal_temperature = nominal_temperature, x_shift = x_shift, baseline_order = len(baseline_terms)-1)
        
class Generate_FitParam_File:
    """Class generates the parameter files used in fitting and sets up fit settings.
    
    Parameters
    ----------
    dataset : object
        Dataset object to be used in the fits
    param_linelist : dataframe
        parameter linelist dataframe name, where the dataframe has the appropriate column headers
    base_linelist : dataframe
        baseline parameter dataframe name generated from the dataset.generate_baseline_paramlist() function.
    lineprofile : str
        lineprofile to use for the simulation. This sets values in the parameter linelist to 0 and forces them to not vary unless manually changed. Default values is VP, voigt profile. Options are VP, SDVP, NGP, SDNGP, HTP.
    linemixing : bool
        If False, then all linemixing parameters are set to 0 and not able to float.  Default is False.
    threshold_intensity : float
        This is the minimum line intensity that will be simulated. Default value is 1e-30.
    fit_intensity  : float
        This is the minimum line intensity for which parameters will be set to float. This can be overwritten manually. Default value is 1e-26.
    fit_window : float
        currently not used
    sim_window : float
        This is the region outside of the wavelength region of the dataset that will be simulated for analysis. Default value is 5 cm-1.
    base_linelist_savename : str
        filename that the baseline linelist will be saved as. Default is Baseline_LineList
    param_linelist_savename : str
        filename that the parameter linelist will be saved as. Default is Parameter_LineList.
    nu_constrain : bool
        True indicates that the line centers will be a global variable across all spectra. False generates a new value for each spectrum in the dataset.
    sw_constrain : bool
        True indicates that the line intensities will be a global variable across all spectra. False generates a new value for each spectrum in the dataset.
    gamma0_constrain : bool
        True indicates that the collisional width will be a global variable across all spectra. False generates a new value for each spectrum in the dataset.
    delta0_constrain : bool
        True indicates that the shift will be a global variable across all spectra. False generates a new value for each spectrum in the dataset.
    aw_constrain : bool
        True indicates that the speed dependent width will be a global variable across all spectra. False generates a new value for each spectrum in the dataset.
    as_constrain : bool
        True indicates that the speed dependent shift will be a global variable across all spectra. False generates a new value for each spectrum in the dataset.
    nuVC_constrain : bool
        True indicates that the dicke narrowing term will be a global variable across all spectra. False generates a new value for each spectrum in the dataset.
    eta_constrain : bool
        True indicates that the correlation parameter will be a global variable across all spectra. False generates a new value for each spectrum in the dataset.
    linemixing_constrain : bool
        True indicates that the first-order linemixing term will be a global variable across all spectra. False generates a new value for each spectrum in the datas
    """
    def __init__ (self, dataset, param_linelist, base_linelist, 
                  lineprofile = 'VP', linemixing = False, threshold_intensity = 1e-30, fit_intensity = 1e-26, fit_window = 1.5, sim_window = 5, 
                  param_linelist_savename = 'Parameter_LineList', base_linelist_savename = 'Baseline_LineList', 
                 nu_constrain = True, sw_constrain = True, gamma0_constrain = True, delta0_constrain = True, aw_constrain = True, as_constrain = True, 
                 nuVC_constrain = True, eta_constrain =True, linemixing_constrain = True):
        self.dataset = dataset
        self.param_linelist = param_linelist
        self.base_linelist = base_linelist
        self.lineprofile = lineprofile
        self.linemixing = linemixing
        self.threshold_intensity = threshold_intensity
        self.fit_intensity = fit_intensity
        self.fit_window = fit_window
        self.sim_window = sim_window
        self.base_linelist_savename = base_linelist_savename
        self.param_linelist_savename = param_linelist_savename
        self.nu_constrain = nu_constrain
        self.sw_constrain = sw_constrain
        self.gamma0_constrain = gamma0_constrain
        self.delta0_constrain = delta0_constrain
        self.aw_constrain = aw_constrain
        self.as_constrain = as_constrain
        self.nuVC_constrain = nuVC_constrain
        self.eta_constrain = eta_constrain
        self.linemixing_constrain = linemixing_constrain
    def get_dataset(self):
        return self.dataset
    def get_param_linelist(self):
        return self.param_linelist
    def get_base_linelist(self):
        return self.base_linelist
    def generate_fit_param_linelist_from_linelist(self, vary_nu = {}, vary_sw = {},
                                   vary_gamma0 = {}, vary_n_gamma0 = {}, 
                                   vary_delta0 = {}, vary_n_delta0 = {}, 
                                   vary_aw = {}, vary_n_gamma2 = {}, 
                                   vary_as = {}, vary_n_delta2 = {}, 
                                   vary_nuVC = {}, vary_n_nuVC = {},
                                   vary_eta = {}, vary_linemixing = {}):
        """Generates the actual parameter line list used and updated in fitting.
               
        

        Parameters
        ----------
        vary_nu : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope line centers should be floated.  Each key in the dictionary corresponds to a molecule, where the values are a dictionary corresponding to local_iso_ids as keys and boolean values acting as boolean flags. True means float the nu for that molecule and isotope. Example vary_nu = {7:{1:True, 2: False, 3: False}, 1:{1:False, 2: False, 3: False}}, would indicate that line centers of all main oxygen isotopes would be floated (if the line intensity is greater than the fit intensity), where the line centers for all other O2 isotopes and all water isotopes would not be varied. If these value are left blank, then all variable would have to be manually switched to float. The default is {}.
        vary_sw : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope line intensities should be floated.  Follows nu_vary example.  The default is {}.
        vary_gamma0 : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope collisional half-width should be floated.  Follows nu_vary example.  . The default is {}.
        vary_n_gamma0 : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope coefficient of the temperature dependence of the half width  should be floated.  Follows nu_vary example.  . The default is {}.
        vary_delta0 : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope pressure shift should be floated.  Follows nu_vary example.  . The default is {}.
        vary_n_delta0 : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope coefficient of the temperature dependence of the pressure shift should be floated.  Follows nu_vary example.  . The default is {}.
        vary_aw : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope ratio of the speed dependent width to the half-width should be floated.  Follows nu_vary example.  . The default is {}.
        vary_n_gamma2 : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope coefficient of the temperature dependence of the speed dependent width should be floated.  Follows nu_vary example.  . The default is {}.
        vary_as : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope the ratio of the speed dependent shift to the shift should be floated.  Follows nu_vary example.  . The default is {}.
        vary_n_delta2 : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope coefficient of the temperature dependence of the speed dependent shift should be floated.  Follows nu_vary example.  . The default is {}.
        vary_nuVC : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope dicke narrowing should be floated.  Follows nu_vary example.  . The default is {}.
        vary_n_nuVC : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope  coefficient of the temperature dependence of the dicke narrowing should be floated.  Follows nu_vary example.  . The default is {}.
        vary_eta : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope correlation parameter for the VC and SD effects should be floated.  Follows nu_vary example.  . The default is {}.
        vary_linemixing : bool, optional
            Dictionary of dictionaries setting whether the molecule and isotope first-order line-mixing should be floated.  Follows nu_vary example.  . The default is {}.

        Returns
        -------
        param_linelist_df : dataframe
            returns dataframe based on parameter line list with addition of a vary and err column for every floatable parameter.  The vary columns are defined by the inputs and the fit_intensity value.  The err columns will be populated from fit results.  The dataframe is also saved as a .csv file.  line intensity will be normalized by the fit_intensity (set to the sw_scale_factor). The term 'sw' is now equal to the normalized value, such that in the simulation 'sw'*sw_scale_factor is used for the line intensity. Because line intensities are so small it is difficult to fit the value without normalization. 

        """
        
        param_linelist_df = self.get_param_linelist().copy()
        #Intensity Thresholding
        param_linelist_df = param_linelist_df[param_linelist_df['sw'] > self.threshold_intensity] #intensity thresholding
        #Cutdown linelist to parameters within specified simulation window of minimum and maximum of dataset
        dataset_min, dataset_max = (self.dataset.get_spectra_extremes())
        extreme_dictionary = self.dataset.get_spectrum_extremes()
        param_linelist_df = param_linelist_df[param_linelist_df['nu'] < (dataset_max + self.sim_window)]
        param_linelist_df = param_linelist_df[param_linelist_df['nu'] > (dataset_min - self.sim_window)]
        #delete parameters not partaining to species (where relevant)
        diluent_list = []
        for spectrum in self.dataset.spectra:
            for diluent in spectrum.Diluent:
                if diluent not in diluent_list:
                    diluent_list.append(diluent)
        num_nominal_temps, list_nominal_temps = self.dataset.get_number_nominal_temperatures()
        column_list = ['molec_id', 'local_iso_id','elower', 'nu', 'sw']
        for diluent in diluent_list:
            column_list.append('gamma0_' + diluent)
            column_list.append('n_gamma0_' + diluent)
            column_list.append('delta0_' + diluent)
            column_list.append('n_delta0_' + diluent)
            column_list.append('SD_gamma_' + diluent)
            column_list.append('n_gamma2_' + diluent)
            column_list.append('SD_delta_' + diluent)
            column_list.append('n_delta2_' + diluent)
            column_list.append('nuVC_' + diluent)
            column_list.append('n_nuVC_' + diluent)
            column_list.append('eta_' + diluent)
            for nominal_temperature in list_nominal_temps:
                column_list.append('y_' + diluent + '_' + str(nominal_temperature)) ## Fix this so it has the nominal temperatures
        param_linelist_df = param_linelist_df[column_list]
        param_linelist_df = param_linelist_df.reset_index(drop = True)
        #Re-defines the Line intensity as sw*sw_scale_factor
        param_linelist_df['sw'] = param_linelist_df['sw'] / self.fit_intensity
        param_linelist_df['sw_scale_factor'] = [self.fit_intensity]*len(param_linelist_df)
        
        # Defines the Linecenter parameters in the event that the line center is held constant across all samples and is not
        ## Starting point is equal to the inital value
        order_nu = ['nu']
        param_linelist_df['nu_vary'] = len(param_linelist_df)*[False]
        param_linelist_df['nu_err'] = len(param_linelist_df)*[0]
        if self.nu_constrain:
            if vary_nu != {}:
                for molecule in vary_nu:
                    for isotope in vary_nu[molecule]:
                        param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) & (param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'nu_vary'] = (vary_nu[molecule][isotope])
        else:
            for spec in self.dataset.get_list_spectrum_numbers():
                param_linelist_df['nu_' + str(spec)] = param_linelist_df['nu'].values
                order_nu.append('nu_' + str(spec))
                param_linelist_df['nu_' + str(spec) + '_vary'] = len(param_linelist_df)*[False]
                param_linelist_df['nu_' + str(spec) + '_err'] = len(param_linelist_df)*[0]
                if vary_nu != {}:
                    for molecule in vary_nu:
                        for isotope in vary_nu[molecule]:
                            param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) & (param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'nu_' + str(spec) + '_vary'] = (vary_nu[molecule][isotope])             
        
        # Defines the linestrength in the event that the sw is held constant across all samples and not
        ## Starting point is equal to the initial value
        order_sw = ['sw', 'sw_scale_factor']
        param_linelist_df['sw_vary'] = len(param_linelist_df)*[False]
        param_linelist_df['sw_err'] = len(param_linelist_df)*[0]
        if self.sw_constrain:
            if vary_sw != {}:
                 for molecule in vary_sw:
                    for isotope in vary_sw[molecule]:
                        param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'sw_vary'] = (vary_sw[molecule][isotope])
        else:
            for spec in self.dataset.get_list_spectrum_numbers():
                param_linelist_df['sw_' + str(spec)] = param_linelist_df['sw'].values
                order_sw.append('sw_' + str(spec))
                param_linelist_df['sw_' + str(spec) + '_vary'] = len(param_linelist_df)*[False]
                param_linelist_df['sw_' + str(spec) + '_err'] = len(param_linelist_df)*[0]
                if vary_sw != {}:
                    for molecule in vary_sw:
                        for isotope in vary_sw[molecule]:
                            param_linelist_df.loc[(param_linelist_df['nu'] >= extreme_dictionary[spec][0])&(param_linelist_df['nu'] <= extreme_dictionary[spec][1])&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'sw_' + str(spec) + '_vary'] = (vary_sw[molecule][isotope])
        
        #Loop through other parameters and then set things to 0 based on lineshape 
        order_gamma0 = []
        order_delta0 = []
        order_SD_gamma = []
        order_SD_delta = []
        order_nuVC = []
        order_eta = []
        order_linemixing = []
        for diluent in diluent_list:
            #Gamma0 option for constrain and not constrained
            order_gamma0.append('gamma0_' + diluent)
            param_linelist_df['gamma0_' + diluent + '_vary'] = len(param_linelist_df)*[False]
            param_linelist_df['gamma0_' + diluent + '_err'] = len(param_linelist_df)*[0]
            if self.gamma0_constrain:
                if vary_gamma0 != {}:
                    for molecule in vary_gamma0:
                        for isotope in vary_gamma0[molecule]:
                            param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'gamma0_' +diluent + '_vary'] = (vary_gamma0[molecule][isotope])
            else:
                for spec in self.dataset.get_list_spectrum_numbers():
                    order_gamma0.append('gamma0_' +diluent + '_' +str(spec))
                    param_linelist_df['gamma0_' +diluent + '_' +str(spec)] = (param_linelist_df['gamma0_' + diluent].values)
                    param_linelist_df['gamma0_' + diluent + '_'+str(spec) + '_vary'] = len(param_linelist_df)*[False]
                    param_linelist_df['gamma0_'+ diluent + '_'+ str(spec) + '_err'] = len(param_linelist_df)*[0]
                    if vary_gamma0 != {}:
                        for molecule in vary_gamma0:
                            for isotope in vary_gamma0[molecule]:
                                param_linelist_df.loc[(param_linelist_df['nu'] >= extreme_dictionary[spec][0])&(param_linelist_df['nu'] <= extreme_dictionary[spec][1])&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'gamma0_' +diluent +'_'+str(spec) + '_vary'] = (vary_gamma0[molecule][isotope])
            order_gamma0.append('n_gamma0_' +diluent )

            #Delta0 option for constrain and not constrained
            order_delta0.append('delta0_' + diluent)
            param_linelist_df['delta0_' + diluent + '_vary'] = len(param_linelist_df)*[False]
            param_linelist_df['delta0_' + diluent + '_err'] = len(param_linelist_df)*[0]
            if self.delta0_constrain:
                if vary_delta0 != {}:
                    for molecule in vary_delta0:
                        for isotope in vary_delta0[molecule]:
                            param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'delta0_' +diluent + '_vary'] = (vary_delta0[molecule][isotope])
            else:
                for spec in self.dataset.get_list_spectrum_numbers():
                    order_delta0.append('delta0_' +diluent + '_' +str(spec))
                    param_linelist_df['delta0_' +diluent + '_' +str(spec)] = (param_linelist_df['delta0_' + diluent].values)
                    param_linelist_df['delta0_' +  diluent + '_'+str(spec) + '_vary'] = len(param_linelist_df)*[0]
                    param_linelist_df['delta0_' + diluent + '_'+ str(spec) + '_err'] = len(param_linelist_df)*[0]
                    if vary_delta0 != {}:
                        for molecule in vary_delta0:
                            for isotope in vary_delta0[molecule]:
                                param_linelist_df.loc[(param_linelist_df['nu'] >= extreme_dictionary[spec][0])&(param_linelist_df['nu'] <= extreme_dictionary[spec][1])&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'delta0_' +diluent +'_'+str(spec) + '_vary'] = (vary_delta0[molecule][isotope])
            order_delta0.append('n_delta0_' +diluent )
 
            #SD Gamma option for constrain and not constrained
            order_SD_gamma.append('SD_gamma_' + diluent )
            param_linelist_df['SD_gamma_' + diluent + '_vary'] = len(param_linelist_df)*[False]
            param_linelist_df['SD_gamma_' + diluent + '_err'] = len(param_linelist_df)*[0]
            if self.aw_constrain:
                if (self.lineprofile == 'VP') or (self.lineprofile == 'NGP'):
                    param_linelist_df.loc[:, 'SD_gamma_' + diluent] = 0
                else:
                    if vary_aw != {}:
                        for molecule in vary_aw:
                            for isotope in vary_aw[molecule]:
                                param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'SD_gamma_' +diluent + '_vary'] = (vary_aw[molecule][isotope])
            else:
                for spec in self.dataset.get_list_spectrum_numbers():
                    order_SD_gamma.append('SD_gamma_' +diluent + '_' +str(spec))
                    param_linelist_df['SD_gamma_' + diluent + '_'+ str(spec) + '_vary'] = len(param_linelist_df)*[False]
                    param_linelist_df['SD_gamma_' + diluent + '_'+ str(spec) + '_err'] = len(param_linelist_df)*[0]
                    if (self.lineprofile == 'VP') or (self.lineprofile == 'NGP'):
                        param_linelist_df.loc[:, 'SD_gamma_' +diluent + '_' +str(spec)] = 0
                    else:
                        param_linelist_df['SD_gamma_' +diluent + '_' +str(spec)] = (param_linelist_df['SD_gamma_' + diluent].values)
                        if vary_aw != {}:
                            for molecule in vary_aw:
                                for isotope in vary_aw[molecule]:
                                    param_linelist_df.loc[(param_linelist_df['nu'] >= extreme_dictionary[spec][0])&(param_linelist_df['nu'] <= extreme_dictionary[spec][1])&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'SD_gamma_' +diluent +'_'+str(spec) + '_vary'] = (vary_aw[molecule][isotope])
            order_SD_gamma.append('n_gamma2_' +diluent )

            #SD Delta option for constrain and not constrained
            order_SD_delta.append('SD_delta_' + diluent)
            param_linelist_df['SD_delta_' + diluent + '_vary'] = len(param_linelist_df)*[False]
            param_linelist_df['SD_delta_' + diluent + '_err'] = len(param_linelist_df)*[0]
            if self.as_constrain:
                if (self.lineprofile == 'VP') or (self.lineprofile == 'NGP'):
                    param_linelist_df.loc[:, 'SD_delta_' + diluent] = 0
                else:
                    if vary_as != {}:
                        for molecule in vary_as:
                            for isotope in vary_as[molecule]:
                                param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'SD_delta_' +diluent + '_vary'] = (vary_as[molecule][isotope])
            else:
                for spec in self.dataset.get_list_spectrum_numbers():
                    order_SD_delta.append('SD_delta_' +diluent + '_' +str(spec))
                    param_linelist_df['SD_delta_' + diluent + '_'+ str(spec) + '_vary'] = len(param_linelist_df)*[False]
                    param_linelist_df['SD_delta_' + diluent + '_'+ str(spec) + '_err'] = len(param_linelist_df)*[0]
                    if (self.lineprofile == 'VP') or (self.lineprofile == 'NGP'):
                        param_linelist_df.loc[:, 'SD_delta_' +diluent + '_' +str(spec)] = 0
                    else:
                        param_linelist_df['SD_delta_' +diluent + '_' +str(spec)] = (param_linelist_df['SD_deltaamma_' + diluent].values)

                        if vary_as != {}:
                            for molecule in vary_as:
                                for isotope in vary_as[molecule]:
                                    param_linelist_df.loc[(param_linelist_df['nu'] >= extreme_dictionary[spec][0])&(param_linelist_df['nu'] <= extreme_dictionary[spec][1])&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'SD_delta_' +diluent +'_'+str(spec) + '_vary'] = (vary_as[molecule][isotope])           
            order_SD_delta.append('n_delta2_' +diluent )

            #nuVC option for constrain and not constrained
            order_nuVC.append('nuVC_' + diluent)
            param_linelist_df['nuVC_' + diluent + '_vary'] = len(param_linelist_df)*[False]
            param_linelist_df['nuVC_' + diluent + '_err'] = len(param_linelist_df)*[0]
            if self.nuVC_constrain:
                if (self.lineprofile == 'VP') or (self.lineprofile == 'SDVP'):
                    param_linelist_df.loc[:, 'nuVC_' + diluent ] = 0
                else:
                    if vary_nuVC != {}:
                        for molecule in vary_nuVC:
                            for isotope in vary_nuVC[molecule]:
                                param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'nuVC_' +diluent + '_vary'] = (vary_nuVC[molecule][isotope])
            else:
                for spec in self.dataset.get_list_spectrum_numbers():
                    order_nuVC.append('nuVC_' + diluent + '_'+str(spec))
                    param_linelist_df['nuVC_' + diluent + '_'+str(spec) + '_vary'] = len(param_linelist_df)*[False]
                    param_linelist_df['nuVC_' + diluent + '_'+ str(spec) + '_err'] = len(param_linelist_df)*[0]
                    if (self.lineprofile == 'VP') or (self.lineprofile == 'SDVP'):
                        param_linelist_df.loc[:, 'nuVC_' + diluent + '_'+str(spec)] = 0                       
                    else:
                        param_linelist_df['nuVC_' +diluent + '_' +str(spec)] = (param_linelist_df['nuVC_' + diluent].values)

                        if vary_nuVC != {}:
                            for molecule in vary_nuVC:
                                for isotope in vary_nuVC[molecule]:
                                    param_linelist_df.loc[(param_linelist_df['nu'] >= extreme_dictionary[spec][0])&(param_linelist_df['nu'] <= extreme_dictionary[spec][1])&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'nuVC_' +diluent +'_'+str(spec) + '_vary'] = (vary_nuVC[molecule][isotope])
            order_nuVC.append('n_nuVC_' +diluent )

            #eta option for constrain and not constrained
            order_eta.append('eta_' + diluent)
            param_linelist_df['eta_' + diluent + '_vary'] = len(param_linelist_df)*[False]
            param_linelist_df['eta_' + diluent + '_err'] = len(param_linelist_df)*[0]
            if self.eta_constrain:
                if (self.lineprofile == 'VP') or (self.lineprofile == 'SDVP') or (self.lineprofile == 'NGP') or (self.lineprofile == 'SDNGP'):
                    param_linelist_df.loc[:, 'eta_' + diluent] = 0
                else:
                    if vary_eta != {}:
                        for molecule in vary_eta:
                            for isotope in vary_eta[molecule]:
                                param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'eta_' +diluent + '_vary'] = (vary_eta[molecule][isotope])
            else:
                for spec in self.dataset.get_list_spectrum_numbers():
                    order_eta.append('eta_' +diluent + '_' +str(spec))
                    param_linelist_df['eta_' + str(spec) + '_vary'] = len(param_linelist_df)*[False]
                    param_linelist_df['eta_' + str(spec) + '_err'] = len(param_linelist_df)*[0]
                    if (self.lineprofile == 'VP') or (self.lineprofile == 'SDVP') or (self.lineprofile == 'NGP') or (self.lineprofile == 'SDNGP'):
                        param_linelist_df.loc[:, 'eta_' +diluent + '_' +str(spec)] = 0                       
                    else:
                        param_linelist_df['eta_' +diluent + '_' +str(spec)] = (param_linelist_df['eta_' + diluent].values)
                        if vary_eta != {}:
                            for molecule in vary_eta:
                                for isotope in vary_eta[molecule]:
                                    param_linelist_df.loc[(param_linelist_df['nu'] >= extreme_dictionary[spec][0])&(param_linelist_df['nu'] <= extreme_dictionary[spec][1])&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'eta_' +diluent +'_'+str(spec) + '_vary'] = (vary_eta[molecule][isotope])
            
            # Linemixing
            for nominal_temp in list_nominal_temps:
                order_linemixing.append('y_' + diluent + '_'+ str(nominal_temp))
                param_linelist_df['y_' + diluent + '_'+ str(nominal_temp) + '_vary'] = len(param_linelist_df)*[False]
                param_linelist_df['y_' + diluent + '_'+ str(nominal_temp)+ '_err'] = len(param_linelist_df)*[0]
                if self.linemixing_constrain:
                    if not self.linemixing:
                        param_linelist_df.loc[:,'y_' + diluent + '_' + str(nominal_temp)] = 0
                    else:
                        if vary_linemixing != {}:
                            for molecule in vary_linemixing:
                                for isotope in vary_linemixing[molecule]:
                                    param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'y_' +diluent + '_'+ str(nominal_temp)+ '_vary'] = (vary_linemixing[molecule][isotope])
                else:
                    for spec in self.dataset.get_list_spectrum_numbers():
                        order_linemixing.append('y_' + diluent + '_'+ str(nominal_temp)+ '_'+ str(spec))
                        param_linelist_df['y_' +diluent + '_'+ str(nominal_temp) + '_'+ str(spec) + '_vary'] = len(param_linelist_df)*[False]
                        param_linelist_df['y_' +diluent + '_'+ str(nominal_temp) + '_'+ str(spec) + '_err'] = len(param_linelist_df)*[0]
                        if not self.linemixing:
                            param_linelist_df.loc[:,'y_' + diluent + '_'+str(nominal_temp) +'_' + str(spec)] = 0
                        else:
                            param_linelist_df['y_' +diluent + '_'+ str(nominal_temp)+'_' +str(spec)] = (param_linelist_df['y_' + diluent].values)
                            if vary_linemixing != {}:
                                for molecule in vary_linemixing:
                                    for isotope in vary_linemixing[molecule]:
                                        param_linelist_df.loc[(param_linelist_df['nu'] >= extreme_dictionary[spec][0])&(param_linelist_df['nu'] <= extreme_dictionary[spec][1])&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'y_' +diluent+ '_'+ str(nominal_temp) +'_'+str(spec) + '_vary'] = (vary_linemixing[molecule][isotope])

            #Temperature Dependence
            if num_nominal_temps > 1:
                param_linelist_df['n_gamma0_'+diluent+'_vary'] = len(param_linelist_df)*[False]
                param_linelist_df['n_gamma0_'+diluent+'_err'] = len(param_linelist_df)*[0]
                param_linelist_df['n_delta0_'+diluent+'_vary'] = len(param_linelist_df)*[False]
                param_linelist_df['n_delta0_'+diluent+'_err'] = len(param_linelist_df)*[0]
                param_linelist_df['n_gamma2_'+diluent+'_vary'] = len(param_linelist_df)*[False]
                param_linelist_df['n_gamma2_'+diluent+'_err'] = len(param_linelist_df)*[0]
                param_linelist_df['n_delta2_'+diluent+'_vary'] = len(param_linelist_df)*[False]
                param_linelist_df['n_delta2_'+diluent+'_err'] = len(param_linelist_df)*[0]
                param_linelist_df['n_nuVC_'+diluent+'_vary'] = len(param_linelist_df)*[False]
                param_linelist_df['n_nuVC_'+diluent+'_err'] = len(param_linelist_df)*[0]
                #n_Gamma0
                if vary_n_gamma0 != {}:
                    for molecule in vary_n_gamma0:
                        for isotope in vary_n_gamma0[molecule]:
                            param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'n_gamma0_' +diluent + '_vary'] = (vary_n_gamma0[molecule][isotope])
                #n_Delta0
                if vary_n_delta0 != {}:
                    for molecule in vary_n_delta0:
                        for isotope in vary_n_delta0[molecule]:
                            param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'n_delta0_' +diluent + '_vary'] = (vary_n_delta0[molecule][isotope])
                #n_Gamma2
                if not (self.lineprofile == 'VP') or  not (self.lineprofile == 'NGP') :
                    if vary_n_gamma2 != {}:
                        for molecule in vary_n_gamma2:
                            for isotope in vary_n_gamma2[molecule]:
                                param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'n_gamma2_' +diluent + '_vary'] = (vary_n_gamma2[molecule][isotope])
                #n_Delta2
                if not (self.lineprofile == 'VP') or  not (self.lineprofile == 'NGP') :
                    if vary_n_delta2 != {}:
                        for molecule in vary_n_delta2:
                            for isotope in vary_n_delta2[molecule]:
                                param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'n_delta2_' +diluent + '_vary'] = (vary_n_delta2[molecule][isotope])
                #n_nuVC
                if not (self.lineprofile == 'VP') or  not (self.lineprofile == 'SDVP') :
                    if vary_n_nuVC != {}:
                        for molecule in vary_n_nuVC:
                            for isotope in vary_n_nuVC[molecule]:
                                param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'n_nuVC_' +diluent + '_vary'] = (vary_n_nuVC[molecule][isotope])
        ordered_list = ['molec_id', 'local_iso_id','elower']
        
        for item in order_nu:
            ordered_list.append(item)
            ordered_list.append(item + '_err')
            ordered_list.append(item + '_vary')
        for item in order_sw:
            ordered_list.append(item)
            if item != 'sw_scale_factor':
                ordered_list.append(item + '_err')
                ordered_list.append(item + '_vary')
        for item in order_gamma0:
            ordered_list.append(item)
            if num_nominal_temps > 1:
                ordered_list.append(item + '_err')
                ordered_list.append(item + '_vary')
            else:
                if 'n_' != item[:2]:
                    ordered_list.append(item + '_err')
                    ordered_list.append(item + '_vary')
        for item in order_delta0:
            ordered_list.append(item)
            if num_nominal_temps > 1:
                ordered_list.append(item + '_err')
                ordered_list.append(item + '_vary')
            else:
                if 'n_' != item[:2]:
                    ordered_list.append(item + '_err')
                    ordered_list.append(item + '_vary')
        for item in order_nuVC:
            ordered_list.append(item)
            if num_nominal_temps > 1:
                ordered_list.append(item + '_err')
                ordered_list.append(item + '_vary')
            else:
                if 'n_' != item[:2]:
                    ordered_list.append(item + '_err')
                    ordered_list.append(item + '_vary')
        for item in order_SD_gamma:
            ordered_list.append(item)
            if num_nominal_temps > 1:
                ordered_list.append(item + '_err')
                ordered_list.append(item + '_vary')
            else:
                if 'n_' != item[:2]:
                    ordered_list.append(item + '_err')
                    ordered_list.append(item + '_vary')
        for item in order_SD_delta:
            ordered_list.append(item)
            if num_nominal_temps > 1:
                ordered_list.append(item + '_err')
                ordered_list.append(item + '_vary')
            else:
                if 'n_' != item[:2]:
                    ordered_list.append(item + '_err')
                    ordered_list.append(item + '_vary')
        for item in order_eta:
            ordered_list.append(item)
            ordered_list.append(item + '_err')
            ordered_list.append(item + '_vary')
        for item in order_linemixing:
            ordered_list.append(item)
            ordered_list.append(item + '_err')
            ordered_list.append(item + '_vary')
        param_linelist_df = param_linelist_df[ordered_list]
        param_linelist_df.to_csv(self.param_linelist_savename + '.csv')
        return param_linelist_df 

    def generate_fit_baseline_linelist(self, vary_baseline = True, vary_pressure = False, vary_temperature = False,vary_molefraction = {}, vary_xshift = False, 
                                      vary_etalon_amp= False, vary_etalon_freq= False, vary_etalon_phase= False):
        """Generates the actual baseline line list used and updated in fitting.  
        

        Parameters
        ----------
        vary_baseline : bool, optional
            If True, then sets all baseline parameters for all spectra to float. The default is True.
        vary_pressure : bool, optional
            If True, then the pressures for all spectra are floated.  This should be used with caution as the impact these parameters have on other floated parameters might lead to an unstable solution. The default is False.
        vary_temperature : bool, optional
            If True, then the temperatures for all spectra are floated. This should be used with caution as the impact these parameters have on other floated parameters might lead to an unstable solution.The default is False.
        vary_molefraction : dict, optional
            keys to dictionary correspond to molecule id where the value is boolean flag, which dictates whether to float the mole fraction. The default is {}. Example: {7: True, 1: False}
        vary_xshift : bool, optional
            If True, then sets x-shift parameters for all spectra to float. The default is False.
        vary_etalon_amp : bool, optional
            If True, then sets etalon amplitude parameters for all spectra to float. The default is False.
        vary_etalon_freq : bool, optional
            If True, then sets etalon period parameters for all spectra to float. . The default is False.
        vary_etalon_phase : bool, optional
            If True, then sets etalon phase parameters for all spectra to float.. The default is False.

        Returns
        -------
        base_linelist_df : dataframe
            returns dataframe based on baseline line list with addition of a vary and err column for every floatable parameter.  The vary columns are defined by the inputs.  The err columns will be populated from fit results.  The dataframe is also saved as a .csv file..



        """
        base_linelist_df = self.get_base_linelist().copy()
        parameters =  (list(base_linelist_df))
        for param in parameters:
            if ('Baseline Order' != param) and ('Segment Number' != param):
                base_linelist_df[param + '_err'] = 0
                base_linelist_df[param + '_vary']= False
            if 'Pressure' in param:
                base_linelist_df[param + '_vary'] = len(base_linelist_df)*[(vary_pressure)]
                if (vary_pressure):
                    print ('USE CAUTION WHEN FLOATING PRESSURES')
            if 'Temperature' in param:
                base_linelist_df[param + '_vary'] = len(base_linelist_df)*[(vary_temperature)]
                if (vary_temperature):
                    print ('USE CAUTION WHEN FLOATING TEMPERATURES')
            if 'x_shift' in param:
                base_linelist_df[param + '_vary'] = len(base_linelist_df)*[(vary_xshift)]
            if 'baseline' in param:
                order = ord(param.replace('baseline_', '')) - 97
                base_linelist_df.loc[base_linelist_df['Baseline Order']>= order, param + '_vary'] = vary_baseline         
            if 'molefraction' in param:
                for molecule in vary_molefraction:
                    if (ISO[(molecule, 1)][4]) in param:  
                        base_linelist_df.loc[base_linelist_df[param]!=0, param + '_vary'] = (vary_molefraction[molecule])                                
            if 'amp' in param:
                base_linelist_df.loc[base_linelist_df[param]!=0, param + '_vary'] = (vary_etalon_amp) 
            if 'freq' in param:
                base_linelist_df.loc[base_linelist_df[param]!=0, param + '_vary'] = (vary_etalon_freq) 
            if 'phase' in param:
                base_linelist_df.loc[base_linelist_df[param.replace("phase", "freq")]!=0, param + '_vary'] = (vary_etalon_phase)
        base_linelist_df.drop(['Baseline Order'], axis=1, inplace = True)
        base_linelist_df = base_linelist_df.reindex(sorted(base_linelist_df.columns), axis=1)
        base_linelist_df.to_csv(self.base_linelist_savename + '.csv')
        return base_linelist_df
    
class Edit_Fit_Param_Files:
    """Allows for the baseline and parameter files to be editted in jupyter notebook.  Can also edit everything in .csv
    
    Parameters
    ----------
    base_linelist_file : str
        name of the .csv file containing the baseline parameters generated in format specified in the generate fit parameters class.
    param_linelist_file : str
        name of the .csv file containing the linelist parameters generated in format specified in the generate fit parameters class.
    new_base_linelist_file : str
        name of the to save the baseline param list as after editing. IF the value is None (default) then it will edit the input.
    new_param_linelist_file  : str
        name of the to save the linelist param list as after editing. IF the value is None (default) then it will edit the input.
    """
    def __init__(self, base_linelist_file, param_linelist_file, new_base_linelist_file = None, new_param_linelist_file = None):
        self.base_linelist_file = base_linelist_file
        self.param_linelist_file = param_linelist_file
        if new_base_linelist_file == None:
            self.new_base_linelist_file = base_linelist_file
        else:
            self.new_base_linelist_file = new_base_linelist_file
        if new_param_linelist_file == None:
            self.new_param_linelist_file = param_linelist_file
        else:
            self.new_param_linelist_file = new_param_linelist_file
    def edit_generated_baselist(self):
        """Allows editting of baseline linelist in notebook        
        """
        base_linelist_df = pd.read_csv(self.base_linelist_file + '.csv', index_col = 0)
        baseline_widget = qgrid.show_grid(base_linelist_df, grid_options={'forceFitColumns': False, 'defaultColumnWidth': 200})
        return baseline_widget
    def save_editted_baselist(self, baseline_widget):
        """Saves edits of baseline linelist in notebook        
        """
        base_linelist_df = baseline_widget.get_changed_df()
        base_linelist_df.to_csv(self.new_base_linelist_file + '.csv')
        return base_linelist_df
    def edit_generated_paramlist(self):
        """Allows editting of parameter linelist in notebook        
        """
        param_linelist_df = pd.read_csv(self.param_linelist_file + '.csv', index_col = 0)
        param_widget = qgrid.show_grid(param_linelist_df, grid_options={'forceFitColumns': False, 'defaultColumnWidth': 200})
        return param_widget
    def save_editted_paramlist(self, param_widget):
        """Saves edits of parameter linelist in notebook        
        """
        param_linelist_df = param_widget.get_changed_df()
        param_linelist_df.to_csv(self.new_param_linelist_file + '.csv') 
        return param_linelist_df
     
def hasNumbers(inputString):
    """Determines whether there are numbers in a string    

    Parameters
    ----------
    inputString : str
        string for analysis

    Returns
    -------
    bool
        DESCRIPTION.

    """
    for char in inputString:
        if char.isdigit():
            return True
    return False

class Fit_DataSet:
    """Provides the fitting functionality
        

        Parameters
        ----------
        dataset : object
            Dataset Object.
        base_linelist_file : str
            filename for file containing baseline parameters.
        param_linelist_file : TYPE
            filename for file containing parmeter parameters..
        minimum_parameter_fit_intensity : float, optional
            minimum intensity for parameters to be generated for fitting. NOTE: Even if a value is floated in the param_linelist if it is below this threshold then it won't be a floated.. The default is 1e-27.
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

        """
    def __init__(self, dataset, base_linelist_file, param_linelist_file, minimum_parameter_fit_intensity = 1e-27,
                baseline_limit = False, baseline_limit_factor = 10, 
                pressure_limit = False, pressure_limit_factor = 10,
                temperature_limit = False, temperature_limit_factor = 10,
                molefraction_limit = False, molefraction_limit_factor = 10, 
                etalon_limit = False, etalon_limit_factor = 50, #phase is constrained to +/- 2pi, 
                x_shift_limit = False, x_shift_limit_magnitude = 0.1, 
                nu_limit = False, nu_limit_magnitude = 0.1, 
                sw_limit = False, sw_limit_factor = 10, 
                gamma0_limit = False, gamma0_limit_factor = 10, n_gamma0_limit= True, n_gamma0_limit_factor = 10, 
                delta0_limit = False, delta0_limit_factor = 10, n_delta0_limit = True, n_delta0_limit_factor = 10, 
                SD_gamma_limit = False, SD_gamma_limit_factor  = 10, n_gamma2_limit = True, n_gamma2_limit_factor  = 10, 
                SD_delta_limit = True, SD_delta_limit_factor  = 10, n_delta2_limit = True, n_delta2_limit_factor  = 10, 
                nuVC_limit = False, nuVC_limit_factor  = 10, n_nuVC_limit = True, n_nuVC_limit_factor = 10, 
                eta_limit = True, eta_limit_factor  = 10, linemixing_limit = False, linemixing_limit_factor  = 10):
            
        self.dataset = dataset
        self.base_linelist_file = base_linelist_file
        self.baseline_list = pd.read_csv(self.base_linelist_file + '.csv')#, index_col = 0
        self.param_linelist_file = param_linelist_file
        self.lineparam_list = pd.read_csv(self.param_linelist_file + '.csv', index_col = 0)
        self.minimum_parameter_fit_intensity = minimum_parameter_fit_intensity
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
        self.SD_delta_limit = SD_gamma_limit
        self.SD_delta_limit_factor = SD_delta_limit_factor
        self.n_delta2_limit = n_delta2_limit
        self.n_delta2_limit_factor = n_delta2_limit_factor
        self.nuVC_limit = nuVC_limit
        self.nuVC_limit_factor = nuVC_limit_factor
        self.n_nuVC_limit = n_nuVC_limit
        self.n_nuVC_limit_factor = n_nuVC_limit_factor
        self.eta_limit = eta_limit
        self.eta_limit_factor = eta_limit_factor
        self.linemixing_limit = linemixing_limit
        self.linemixing_limit_factor = linemixing_limit_factor
        
    def generate_params(self):
        """Generates the lmfit params object that will be used in fitting
        

        Returns
        -------
        params : lmfit parameter object
            the params object is a dictionary comprised of all parameters translated from dataframes into a dictionary format compatible with lmfit

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
        
        #Lineshape parameters
        linelist_params = []
        for line_param in list(self.lineparam_list):
            if (self.dataset.get_number_nominal_temperatures()[0]) == 1:
                if ('_vary' not in line_param) and ('_err' not in line_param) and (line_param != 'molec_id') and (line_param != 'local_iso_id') and (line_param != 'elower') and ('n_' not in line_param):
                    linelist_params.append(line_param)
            else:
                if ('_vary' not in line_param) and ('_err' not in line_param) and (line_param != 'molec_id') and (line_param != 'local_iso_id') and (line_param != 'elower'):
                    linelist_params.append(line_param)
        diluent_list = []
        for spectrum in self.dataset.spectra:
            for diluent in spectrum.Diluent:
                if diluent not in diluent_list:
                    diluent_list.append(diluent)
        nu_constrain = True
        sw_constrain = True
        gamma0_constrain = True
        delta0_constrain = True
        SD_gamma_constrain = True
        SD_delta_constrain = True
        nuVC_constrain = True
        eta_constrain = True
        linemix_constrain = True
        if ((sum(('nu' in param) & ('nuVC' not in param) for param in linelist_params))) > 1:
            nu_constrain = False
        if (sum('sw' in param for param in linelist_params)) > 2:
            sw_constrain = False
        if (self.dataset.get_number_nominal_temperatures()[0]) == 1:
            if (sum('gamma0' in param for param in linelist_params)) >  len(diluent_list):
                gamma0_constrain = False
            if (sum('delta0' in param for param in linelist_params)) >  len(diluent_list):
                delta0_constrain = False
            if (sum('SD_gamma' in param for param in linelist_params)) >  len(diluent_list):
                SD_gamma_constrain = False
            if (sum('SD_delta' in param for param in linelist_params)) >  len(diluent_list):
                SD_delta_constrain = False
            if (sum('nuVC' in param for param in linelist_params)) >  len(diluent_list):
                nuVC_constrain = False
            if (sum('eta' in param for param in linelist_params)) >  len(diluent_list):
                eta_constrain = False 
        else:
            if (sum('gamma0' in param for param in linelist_params)) >  2*len(diluent_list):
                gamma0_constrain = False
            if (sum('delta0' in param for param in linelist_params)) >  2*len(diluent_list):
                delta0_constrain = False
            if (sum('SD_gamma' in param for param in linelist_params)) >  len(diluent_list) + 1:
                SD_gamma_constrain = False
            if (sum('SD_delta' in param for param in linelist_params)) >  len(diluent_list) + 1:
                SD_delta_constrain = False
            if (sum('nuVC' in param for param in linelist_params)) >  2*len(diluent_list) :
                nuVC_constrain = False
            if (sum('eta' in param for param in linelist_params)) >  len(diluent_list):
                eta_constrain = False 
        if (sum('y' in param for param in linelist_params)) >  len(diluent_list)*(self.dataset.get_number_nominal_temperatures()[0]):
            linemix_constrain = False
        linemix_terms_constrained = []
        for diluent in diluent_list:
            for temperature in (self.dataset.get_number_nominal_temperatures()[1]):
                linemix_terms_constrained.append('y_'+diluent + '_' + str(temperature))
        for spec_line in self.lineparam_list.index.values:
            if self.lineparam_list.loc[spec_line]['sw'] >= self.minimum_parameter_fit_intensity / self.lineparam_list.loc[spec_line]['sw_scale_factor']:# bigger than 1 because fit_intensity / fit_intensity
                for line_param in linelist_params:
                    indices = [m.start() for m in re.finditer('_', line_param)]
                    index_length = len(indices)
                    #NU
                    if line_param == 'nu' and nu_constrain:
                        if self.nu_limit:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                      min = self.lineparam_list.loc[spec_line][line_param] - self.nu_limit_magnitude, 
                                      max = self.lineparam_list.loc[spec_line][line_param] + self.nu_limit_magnitude)
                        else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    elif (line_param != 'nu') and ('nu' in line_param) and ('nuVC' not in line_param) and (not nu_constrain):
                        if self.nu_limit:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                      min = self.lineparam_list.loc[spec_line][line_param] - self.nu_limit_magnitude, 
                                      max = self.lineparam_list.loc[spec_line][line_param] + self.nu_limit_magnitude)
                        else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    #SW
                    elif line_param == 'sw' and sw_constrain:
                        if self.sw_limit:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.sw_limit_factor)* self.lineparam_list.loc[int(spec_line)][line_param], 
                                  max = self.sw_limit_factor* self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    elif (line_param != 'sw') and ('sw' in line_param) and (not sw_constrain) and (line_param != 'sw_scale_factor'):
                        if self.sw_limit:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min =  (1 / self.sw_limit_factor)* self.lineparam_list.loc[int(spec_line)][line_param],
                                  max = self.sw_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    #GAMMA0
                    elif ('gamma0_' in line_param) and ('n_' not in line_param) and (gamma0_constrain) and (index_length==1):
                        if self.gamma0_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.gamma0_limit_factor)*self.lineparam_list.loc[int(spec_line)][line_param] ,#(1 - (self.gamma0_limit_percent / 100))*lineparam_list.loc[int(spec_line)][line_param], 
                                  max = self.gamma0_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])#(1 + (self.gamma0_limit_percent / 100))*lineparam_list.loc[int(spec_line)][line_param])
                                
                        else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    elif ('gamma0_' in line_param) and ('n_' not in line_param) and (not gamma0_constrain) and (index_length>1):
                        if self.gamma0_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                             params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.gamma0_limit_factor)*self.lineparam_list.loc[int(spec_line)][line_param],#(1 - (self.gamma0_limit_percent / 100))*lineparam_list.loc[int(spec_line)][line_param], 
                                  max = self.gamma0_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])#(1 + (self.gamma0_limit_percent / 100))*lineparam_list.loc[int(spec_line)][line_param])
                        else:
                             params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    elif ('n_gamma0' in line_param):
                        if self.n_gamma0_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.n_gamma0_limit_factor) *self.lineparam_list.loc[int(spec_line)][line_param], 
                                  max = self.n_gamma0_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param],self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    #DELTA0
                    elif ('delta0' in line_param) and ('n_' not in line_param) and (delta0_constrain) and (index_length==1):
                        if self.delta0_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.delta0_limit_factor )*self.lineparam_list.loc[int(spec_line)][line_param], 
                                  max = self.delta0_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    elif ('delta0_' in line_param) and ('n_' not in line_param) and (not delta0_constrain) and (index_length>1):
                        if self.delta0_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                             params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.delta0_limit_factor)*self.lineparam_list.loc[int(spec_line)][line_param], 
                                  max = self.delta0_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                             params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    elif ('n_delta0' in line_param):
                        if self.n_delta0_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param],self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.n_delta0_limit_factor)*self.lineparam_list.loc[int(spec_line)][line_param], 
                                  max = self.n_delta0_limit_factor / 100*self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    #SD Gamma
                    elif ('SD_gamma' in line_param) and (SD_gamma_constrain) and (index_length==2):
                        if self.SD_gamma_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.SD_gamma_limit_factor) *self.lineparam_list.loc[int(spec_line)][line_param], 
                                  max = self.SD_gamma_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'],
                            min = 0)
                    elif ('SD_gamma' in line_param) and (not SD_gamma_constrain) and (index_length>2):
                        if self.SD_gamma_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                             params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.SD_gamma_limit_factor)*self.lineparam_list.loc[int(spec_line)][line_param], 
                                  max = self.SD_gamma_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                             params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    elif ('n_gamma2' in line_param):
                    
                        if self.n_gamma2_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                min = (1 / self.n_gamma2_limit_factor)*self.lineparam_list.loc[int(spec_line)][line_param], 
                                max = (self.n_gamma2_limit_factor / 100)*self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    #SD Delta
                    elif ('SD_delta' in line_param) and (SD_delta_constrain) and (index_length==2):
                        if self.SD_delta_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.SD_delta_limit_factor )*self.lineparam_list.loc[int(spec_line)][line_param], 
                                  max = self.SD_delta_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    elif ('SD_delta' in line_param) and (not SD_delta_constrain) and (index_length>2):
                        if self.SD_delta_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                             params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.SD_delta_limit_factor )*self.lineparam_list.loc[int(spec_line)][line_param], 
                                  max =self.SD_delta_limit_factor / 100*self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                             params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    elif ('n_delta2' in line_param):

                        if self.n_delta2_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                min = (1 / self.n_delta2_limit_factor )*self.lineparam_list.loc[int(spec_line)][line_param], 
                                max = self.n_delta2_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    #nuVC
                    elif ('nuVC' in line_param) and ('n_nuVC_' not in line_param) and (nuVC_constrain) and (index_length==1):
                        if self.nuVC_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 /self.nuVC_limit_factor)*self.lineparam_list.loc[int(spec_line)][line_param], 
                                  max = self.nuVC_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    elif ('nuVC' in line_param) and ('n_nuVC' not in line_param) and (not nuVC_constrain) and (index_length>1):
                        if self.nuVC_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                             params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.nuVC_limit_factor)*self.lineparam_list.loc[int(spec_line)][line_param], 
                                  max = self.nuVC_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                             params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    elif ('n_nuVC' in line_param):
                        if self.n_nuVC_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.n_nuVC_limit_factor )*self.lineparam_list.loc[int(spec_line)][line_param], 
                                  max = self.n_nuVC_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])                    
                    #eta
                    elif ('eta_' in line_param) and (eta_constrain) and (index_length==1):
                        if self.eta_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.eta_limit_factor)*self.lineparam_list.loc[int(spec_line)][line_param], 
                                  max = (self.eta_limit_factor)*self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    elif ('eta_' in line_param) and (not eta_constrain) and (index_length>1):
                        if self.eta_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                             params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.eta_limit_factor)*self.lineparam_list.loc[int(spec_line)][line_param], 
                                  max = self.eta_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                             params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary']) 
                    # linemixing
                    elif ('y_' in line_param) and (linemix_constrain) and (line_param in linemix_terms_constrained):
                        if self.linemixing_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.linemixing_limit_factor)*self.lineparam_list.loc[int(spec_line)][line_param], 
                                  max = self.linemixing_limit_factor *self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    elif ('y_' in line_param) and (not linemix_constrain) and (not line_param in linemix_terms_constrained):
                        if self.linemixing_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                             params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.linemixing_limit_factor)*self.lineparam_list.loc[int(spec_line)][line_param], 
                                  max = self.linemixing_limit_factor / 100*self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                             params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary']) 
        return (params)

    def constrained_baseline(self, params, baseline_segment_constrained = True, xshift_segment_constrained = True, molefraction_segment_constrained = True,
                                    etalon_amp_segment_constrained = True, etalon_freq_segment_constrained = True, etalon_phase_segment_constrained = True, 
                                    pressure_segment_constrained = True, temperature_segment_constrained = True):
        """Allows you to impose baseline constraints when using multiple segments per spectrum, ie. all baseline parameters can be the same across the entire spectrum except for the etalon phase, which is allowed to vary per segment
        

        Parameters
        ----------
        params : lmfit parameter object
            the params object is a dictionary comprised of all parameters translated from dataframes into a dictionary format compatible with lmfit.
        baseline_segment_constrained : bool, optional
            True means the baseline terms are constrained across each spectrum.. The default is True.
        xshift_segment_constrained : bool, optional
            True means the x_shift terms are constrained across each spectrum.. The default is True.
        molefraction_segment_constrained : bool, optional
            True means the mole fraction for that molecule is constrained across each spectrum.. The default is True.
        etalon_amp_segment_constrained : bool, optional
            True means the etalon amplitude is constrained across each spectrum.. The default is True.
        etalon_freq_segment_constrained : bool, optional
            True means the etalon frequency is constrained across each spectrum.. The default is True.
        etalon_phase_segment_constrained : bool, optional
            True means the etalon phase is constrained across each spectrum.. The default is True.
        pressure_segment_constrained : bool, optional
            True means the pressure is constrained across each spectrum.. The default is True.
        temperature_segment_constrained : bool, optional
            True means the temperature is constrained across each spectrum.. The default is True.

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
                elif 'freq' in param and etalon_freq_segment_constrained:
                    if segment_num != spectrum_segment_min[spectrum_num]:
                        params[param].set(expr = param[:indices[2]+1] + str(spectrum_num) + '_' + str(spectrum_segment_min[spectrum_num]))
                elif 'phase' in param and etalon_phase_segment_constrained:
                    if segment_num != spectrum_segment_min[spectrum_num]:
                        params[param].set(expr = param[:indices[2]+1] + str(spectrum_num) + '_' + str(spectrum_segment_min[spectrum_num]))
        return params

    def simulation_model(self, params, wing_cutoff = 50, wing_wavenumbers = 50, wing_method = 'wing_cutoff'):
        """This is the model that is used for fitting.  It combines the baseline terms (baseline, etalons) and resonant absorption models
        

        Parameters
        ----------
        params : lmfit parameter object
            the params object is a dictionary comprised of all parameters translated from dataframes into a dictionary format compatible with lmfit.
        wing_cutoff : float, optional
            number of voigt half-widths to simulate on either side of each line. The default is 50.
        wing_wavenumbers : float, optional
            number of wavenumbers to simulate on either side of each line. The default is 50.
        wing_method : TYPE, optional
            Provides choice between the wing_cutoff and wing_wavenumbers line cut-off options. The default is 'wing_cutoff'.

        Returns
        -------
        total_residuals : array
            residuals for all spectra in Dataset.

        """
        total_simulated = []
        total_residuals = []
        baseline_params = []
        linelist_params = []
        
        # Set-up Baseline Parameters
        for param in (list(params.valuesdict().keys())):
            if ('molefraction' in param) or ('baseline' in param) or ('etalon' in param) or ('x_shift' in param) or ('Pressure' in param) or ('Temperature' in param):
                baseline_params.append(param)
            else:
                linelist_params.append(param)

        for spectrum in self.dataset.spectra:
            simulated_spectra = len(spectrum.wavenumber)*[0]    
            residuals = len(spectrum.alpha)*[0]        
            wavenumber_segments, alpha_segments, indices_segments = spectrum.segment_wave_alpha()
            Diluent = spectrum.Diluent
            spectrum_number = spectrum.spectrum_number
            nominal_temp = spectrum.nominal_temperature          
            columns = ['molec_id', 'local_iso_id', 'elower', 'nu', 'sw', 'sw_scale_factor']
            for species in Diluent:
                columns.append('gamma0_' + species)
                columns.append('n_gamma0_'+ species)
                columns.append('delta0_' + species)
                columns.append('n_delta0_'+ species)
                columns.append('SD_gamma_' + species)
                columns.append('n_gamma2_'+ species )
                columns.append('SD_delta_' + species)
                columns.append('n_delta2_' + species)
                columns.append('nuVC_' + species)
                columns.append('n_nuVC_' + species)
                columns.append('eta_' + species)
                columns.append('y_' + species + '_' + str(nominal_temp))
            rename_dictionary = {}
            for column in self.lineparam_list:
                if ('vary' not in column) and ('err' not in column):
                    if ((column + '_' + str(spectrum_number)) in self.lineparam_list) and ('y_' not in column):
                        columns = [(column + '_'  + str(spectrum_number))if x==column else x for x in columns]
                        rename_dictionary[(column + '_'  + str(spectrum_number))] = column
                    elif ('y_' in column):
                        if ((column + '_' + str(spectrum_number)) in self.lineparam_list):
                            columns = [(column + '_'  +str(spectrum_number))if x==column else x for x in columns]
                            rename_dictionary[(column) + '_' + str(spectrum_number)] = (column[:column.find(str(nominal_temp))-1])
                        elif column.count('_') < 3:
                            rename_dictionary[(column)] = (column[:column.find(str(nominal_temp))-1])
            linelist_for_sim = self.lineparam_list[columns].copy()
           
            # Replaces the relevant linelist locations with the 
            for parameter in linelist_params:
                line = int(parameter[parameter.find('_line_') + 6:])
                param = parameter[:parameter.find('_line_')]
                if param in columns:
                    if 'sw' in param:
                        #print (linelist_for_sim[line]['sw_scale_factor'])
                        linelist_for_sim.loc[line, param] = np.float(params[parameter])
                    else:
                        linelist_for_sim.loc[line, param] = np.float(params[parameter])
            #Renames columns to generic (no scan number)
            linelist_for_sim=linelist_for_sim.rename(columns = rename_dictionary) 
            linelist_for_sim['sw'] = linelist_for_sim['sw']*linelist_for_sim['sw_scale_factor']
            for segment in list(set(spectrum.segments)):
                wavenumbers = wavenumber_segments[segment]
                wavenumbers_relative = wavenumbers - np.min(spectrum.wavenumber)
                x_shift = np.float(params['x_shift_' + str(spectrum_number) + '_' + str(segment)])
                #linelist_for_sim['nu'] = linelist_for_sim['nu'] + x_shift # Q
                wavenumbers += x_shift
                wavenumbers_relative+= x_shift
                #Set-up MoleFraction for Fitting
                fit_molefraction = spectrum.molefraction
                for molecule in spectrum.molefraction:
                    if ('molefraction_'+ ISO[(molecule, 1)][4]) + '_' + str(spectrum_number) + '_' + str(segment) in baseline_params:
                        fit_molefraction[molecule] = np.float(params[('molefraction_'+ ISO[(molecule, 1)][4]) + '_' + str(spectrum_number) + '_' + str(segment)])
                #Get Environmental Parameters
                p = np.float(params['Pressure_' + str(spectrum_number) + '_' + str(segment)])
                T = np.float(params['Temperature_' + str(spectrum_number) + '_' + str(segment)])               
                #Simulate Spectra
                fit_nu, fit_coef = HTP_from_DF_select(linelist_for_sim, wavenumbers, wing_cutoff = wing_cutoff, wing_wavenumbers = wing_wavenumbers, wing_method = wing_method,
                        p = p, T = T, molefraction = fit_molefraction, 
                        natural_abundance = spectrum.natural_abundance, abundance_ratio_MI = spectrum.abundance_ratio_MI,  Diluent = Diluent)
                fit_coef = fit_coef * 1e6
                ## Baseline Calculation
                baseline_param_array = [0]*(self.dataset.baseline_order+1)
                for param in baseline_params:
                    if ('baseline' in param):
                        indices = [m.start() for m in re.finditer('_', param)]
                        spectrum_num = int(param[indices[1]+1:indices[2]])
                        segment_num = int(param[indices[2]+1:])
                        if (spectrum_num == spectrum_number) and (segment_num == segment): 
                            baseline_param_array[ord(param[9:param.find('_',9)])-97] = np.float(params[param])
                baseline_param_array = baseline_param_array[::-1] # reverses array to be used for polyval
                baseline = np.polyval(baseline_param_array, wavenumbers_relative)
                #Etalon Calculation
                fit_etalon_parameters = {}
                for i in range(1, len(spectrum.etalons)+1):
                    fit_etalon_parameters[i] = {'amp': 0, 'freq':1, 'phase':0}
                for param in baseline_params:
                    if ('etalon' in param) and (str(spectrum_number) in param[param.find('_', 7):]):
                        etalon_num = int(param[param.find('_')+1: param.find('_', param.find('_')+1)])
                        if param == 'etalon_' + str(etalon_num) + '_amp_' + str(spectrum_number) + '_' +str(segment):#('amp' in param) and (str(etalon_num) in param):
                            fit_etalon_parameters[etalon_num]['amp'] = np.float(params[param])
                        if param == 'etalon_' + str(etalon_num) + '_freq_' + str(spectrum_number) + '_' +str(segment):#('freq' in param) and (str(etalon_num) in param):
                            fit_etalon_parameters[etalon_num]['freq'] = np.float(params[param])
                        if param == 'etalon_' + str(etalon_num) + '_phase_' + str(spectrum_number) +'_' + str(segment):#('phase' in param) and (str(etalon_num) in param):
                            fit_etalon_parameters[etalon_num]['phase'] = np.float(params[param])
                etalons = len(wavenumbers)*[0]
                for i in range(1, len(spectrum.etalons)+1):
                    etalons += etalon(wavenumbers_relative, fit_etalon_parameters[i]['amp'], fit_etalon_parameters[i]['freq'], fit_etalon_parameters[i]['phase']) 
                simulated_spectra[np.min(indices_segments[segment]): np.max(indices_segments[segment])+1] = (baseline + etalons + fit_coef)
                residuals[np.min(indices_segments[segment]): np.max(indices_segments[segment])+1]  = (baseline + etalons + fit_coef) - alpha_segments[segment]
            total_simulated = np.append(total_simulated, simulated_spectra)
            total_residuals = np.append(total_residuals, residuals)
        total_residuals = np.asarray(total_residuals)
        total_simulated = np.asarray(total_simulated)
        return total_residuals
    def fit_data(self, params, wing_cutoff = 50, wing_wavenumbers = 50, wing_method = 'wing_cutoff', xtol = 1e-7, maxfev = 2000, ftol = 1e-7):
        """Uses the LMFit minimizer to do the fitting through the simulation model function
        

        Parameters
        ----------
        params : lmfit parameter object
            the params object is a dictionary comprised of all parameters translated from dataframes into a dictionary format compatible with lmfit.
        wing_cutoff : float, optional
            number of voigt half-widths to simulate on either side of each line. The default is 50.
        wing_wavenumbers : float, optional
            number of wavenumbers to simulate on either side of each line. The default is 50.
        wing_method : str, optional
            Provides choice between the wing_cutoff and wing_wavenumbers line cut-off options. The default is 'wing_cutoff'.
        xtol : float, optional
             Absolute error in xopt between iterations that is acceptable for convergence. The default is 1e-7.
        maxfev : float, optional
            DESCRIPTION. The default is 2000.
        ftol : The maximum number of calls to the function., optional
            Absolute error in func(xopt) between iterations that is acceptable for convergence.. The default is 1e-7.

        Returns
        -------
        result : LMFit result Object
            contains all fit results as LMFit results object.

        """
        minner = Minimizer(self.simulation_model, params, xtol =xtol, maxfev =  maxfev, ftol = ftol, fcn_args=(wing_cutoff, wing_wavenumbers, wing_method))
        result = minner.minimize(method = 'leastsq')#'
        return result
    def residual_analysis(self, result, indv_resid_plot = False):
        """Updates the model and residual arrays in each spectrum object with the results of the fit and gives the option of plotting each spectra with residuals.
        

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
            spectrum.set_residuals(spectrum_residual)
            spectrum.set_model(spectrum_residual + spectrum.alpha)
            if indv_resid_plot:
                spectrum.plot_model_residuals()
    def update_params(self, result, base_linelist_update_file = None , param_linelist_update_file = None):
        """Updates the baseline and line parameter files based on fit results with the option to write over the file (default) or save as a new file. It also calculates and add the baseline and etalon terms into a baseline array for the spectrum.
        

        Parameters
        ----------
        result : LMFit result Object
            contains all fit results as LMFit results object.
        base_linelist_update_file : str, optional
            Name of file to save the updated baseline parameters. Default is to override the input. The default is None.
        param_linelist_update_file : str, optional
            Name of file to save the updated line parameters. Default is to override the input. The default is None.

        """
        if base_linelist_update_file == None:
            base_linelist_update_file = self.base_linelist_file
        if param_linelist_update_file == None:
            param_linelist_update_file = self.param_linelist_file
        for key, par in result.params.items():
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
            else:
                parameter = par.name[:par.name.find('_line')]
                line = int(par.name[par.name.find('_line')+6:])
                self.lineparam_list.loc[line, parameter] = par.value
                if par.vary:
                    self.lineparam_list.loc[line, parameter + '_err'] = par.stderr
        self.baseline_list.to_csv(base_linelist_update_file + '.csv', index = False)
        self.lineparam_list.to_csv(param_linelist_update_file + '.csv')
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
                    fit_etalon_parameters[i] = {'amp': 0, 'freq':1, 'phase':0}
                for key, par in result.params.items():
                    if ('baseline' in par.name) and (str(spectrum.spectrum_number) in par.name):
                        baseline_param_array[ord(par.name[9:par.name.find('_',9)])-97] = np.float(par.value)
                    elif ('etalon' in par.name) and (str(spectrum.spectrum_number) in par.name[par.name.find('_', 7):]):
                        etalon_num = int(par.name[par.name.find('_')+1: par.name.find('_', par.name.find('_')+1)])
                        if ('amp' in par.name) and (str(etalon_num) in par.name):
                            fit_etalon_parameters[etalon_num]['amp'] = par.value
                        if ('freq' in par.name) and (str(etalon_num) in par.name):
                            fit_etalon_parameters[etalon_num]['freq'] = par.value
                        if ('phase' in par.name) and (str(etalon_num) in par.name):
                            fit_etalon_parameters[etalon_num]['phase'] = par.value
                baseline_param_array = baseline_param_array[::-1] # reverses array to be used for polyval
                baseline[bound_min: bound_max + 1] += np.polyval(baseline_param_array, wave_rel)
                for i in range(1, len(spectrum.etalons)+1):
                    baseline[bound_min: bound_max +1] += etalon(wave_rel, fit_etalon_parameters[i]['amp'], fit_etalon_parameters[i]['freq'], fit_etalon_parameters[i]['phase'])
            spectrum.set_background(baseline)