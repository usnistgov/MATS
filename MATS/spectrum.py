#Import Packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.interpolate import RegularGridInterpolator

from .utilities import etalon, convolveSpectrumSame

from .hapi import ISO, PYTIPS2011, PYTIPS2017, PYTIPS2021, PYTIPS2025
from .codata import CONSTANTS
from .spectroscopic_model import Spectroscopic_model
from .o2_cia_karman import O2_CIA_Karman_Model




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
        if not at natural abundance sets the enhancement factor for each molecule and isotope in the following format {molec_id:{iso_id: enhancement, iso_id: enhancement}, . . . }.  The enhancement is the ratio of the new abundance to the natural abundance.
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
    TIPS : definition, optional
        selects the HAPI provided TIPS version to use for the partition function
    compressability_file : str, optional
        compressability_file is the name of a file that contains the compressability information, generated by the NIST Refprop program with the rows containing compressability (used to correct the ideal gas law) as a function pressure (MPa) and the columns as a function of temperature (K).  Pressure and temperature are header/index in the given axis.  
    """

    ALLOWED_Y_UNITS = {
        'alpha': ['cm-1', '10-6 cm-1', '10^-6 cm^-1', 'ppm/cm', 'us'],
        'absorbance': ['unitless', 'none', 'a.u.', '', 'ppm', '10-6'],
        'absorption': ['unitless', 'none', 'a.u.', ''],
        'transmittance': ['unitless', 'none', 'a.u.', '', "%"], 
    }

    ALLOWED_Y_UNC_UNITS = {
        'alpha': ['cm-1', '10-6 cm-1', '10^-6 cm^-1','ppm/cm', 'us', '%', 'rel'],
        'absorbance': ['unitless', 'none', 'a.u.', '', '%', 'ppm', '10-6', 'rel'],
        'absorption': ['unitless', 'none', 'a.u.', '', '%', 'rel'],
        'transmittance': ['unitless', 'none', 'a.u.', '', '%', 'rel']
    }

    CIA_CONVERSION = {'cm-1': 1, 
                      'ppm/cm': 1e-6,
                      '10-6 cm-1': 1e-6}

    # 2. X-Axis Conversion (Target: cm-1)
    # 1 cm-1 = 29979.2458 MHz
    X_CONVERSION = {
        'cm-1': 1.0,
        'mhz': 10**6 / CONSTANTS['c'],
        'ghz': 10**9 / CONSTANTS['c']
    }

    PRESSURE_CONVERSION = {
        'atm': 1.0,
        'torr': 1 / 760,
        'mmHg': 1 / 760,
        'bar': 100000/101325,
        'mbar':100 / 101325,
        'pa': 1/101325,
        'kPa': 1/101.325, 

    }

    TEMPERATURE_CONVERSION = {
        'k': 0, 
        'c': 273.15, 
    }


    def __init__(self, filename, molefraction = {}, natural_abundance = True, isotope_list = ISO, diluent = 'air', Diluent = {}, abundance_ratio_MI = {}, spectrum_number = 1,
                    x_column = 'Total Frequency /MHz', x_input_units = 'cm-1', 
                    y_column = 'Mean tau/us',  y_input_units = 'ppm/cm', y_unc_column = None, y_unc_input_units = '',
                    pressure_column = 'Cavity Pressure /Torr', pressure_input_units = 'torr',
                    temperature_column = 'Cavity Temperature Side 2 /C', temperature_input_units = 'C',
                    dataspace = 'alpha', pathlength = 0., 
                     
                    segment_column = None,
                    etalons = {}, nominal_temperature = 296, x_shift = 0.0, baseline_order = 1, weight = 1,
                    ILS_function = None, ILS_resolution = 0.1, ILS_wing = 10, TIPS = PYTIPS2025, 
                    compressability_file = None,
                    cia = None, cia_input_units = None):
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
                print ("If using the self broadening term, then consider explicitly labeling the broadener (ie in an oxygen spectra use 'O2' instead of self).  This may avoid confusion in multiple species fits. ")
            else:
                print ('If using the beta formalism use the Diluent{diluent:{"composition": 1, "m": mass}} format')
                self.Diluent = {self.diluent: {'composition':1, 'm': 0}}
        else:
            self.Diluent = Diluent
            if 'self' in self.Diluent:
                print ("You are using the 'self' term, consider explicitly labeling the broadener (ie in an oxygen spectra use 'O2' instead of 'self').  This may avoid confusion in multiple species fits. For single species fits it should not matter.")
                print ("Double check that you did not include the equivalent of the self term explicitly (ie in an oxygen spectra having both 'O2' and 'self').")
        self.spectrum_number = spectrum_number

        self.dataspace = dataspace.lower()
        if self.dataspace == 'alpha':
            self.pathlength = 0.
        else:
            self.pathlength = pathlength
        
        #Validate Units
        self.x_input_units = x_input_units.lower()
        if self.x_input_units not in self.X_CONVERSION:
            raise ValueError(f"Unknown x_units '{self.x_input_units}'. Add it to X_CONVERSION dict.")

        self.y_input_units = y_input_units.lower()
        if self.y_input_units not in self.ALLOWED_Y_UNITS.get(self.dataspace, []):
            raise ValueError(
                f"Invalid y_units '{self.y_input_units}' for dataspace '{self.dataspace}'. "
                f"Allowed units: {self.ALLOWED_Y_UNITS[self.dataspace]}"
            )
        
        self.y_unc_input_units = y_unc_input_units.lower()
        if y_unc_column is not None:
            if self.y_unc_input_units not in self.ALLOWED_Y_UNC_UNITS.get(self.dataspace, []):
                raise ValueError(
                    f"Invalid y_units '{self.y_unc_input_units}' for dataspace '{self.dataspace}'. "
                    f"Allowed units: {self.ALLOWED_Y_UNC_UNITS[self.dataspace]}"
                )

        self.pressure_input_units = pressure_input_units.lower()
        if self.pressure_input_units not in self.PRESSURE_CONVERSION:
            raise ValueError(f"Unknown x_units '{self.pressure_input_units}'. Add it to PRESSURE_CONVERSION dict.")
        self.temperature_input_units = temperature_input_units.lower()
        if self.temperature_input_units not in self.TEMPERATURE_CONVERSION:
            raise ValueError(f"Unknown x_units '{self.temperature_input_units}'. Add it to TEMPERATURE_CONVERSION dict.")
        
        #Dataframe
        file_contents = pd.read_csv(self.filename + '.csv',float_precision = 'high')
        raw_x = file_contents[x_column].values
        raw_y = file_contents[y_column].values
        raw_pressure = file_contents[pressure_column].values
        raw_temperature =file_contents[temperature_column].values

        if y_unc_column is not None:
            raw_y_unc = file_contents[y_unc_column].values if y_unc_column in file_contents.columns else None
        else:
            raw_y_unc = None
        
        if segment_column is not None:
            self.segments = file_contents[segment_column].values.astype(int)
        else:
            self.segments = len(file_contents)*[1]

        #Convert to DF inputs to Desired Units
        self.wavenumber = raw_x * self.X_CONVERSION[self.x_input_units]
        self.y_data, self.y_unc, self.y_output_units, self.y_output_label = self._convert_y_data(raw_y, raw_y_unc)

        self.pressure_array = raw_pressure * self.PRESSURE_CONVERSION[self.pressure_input_units]
        self.pressure = np.mean(self.pressure_array)
        self.temperature_array = raw_temperature + self.TEMPERATURE_CONVERSION[self.temperature_input_units]
        self.temperature = np.mean(self.temperature_array)

        self.pressure_output_units = 'atm'
        self.temperature_output_units = 'K'
        self.x_output_units = 'cm$^{-1}$'
                

        #Initialize Model, Residual, Background, and CIA
        spectrum_length = len(self.y_data)
        self.model = np.zeros(spectrum_length)
        self.residuals = self.y_data - self.model
        self.background = np.zeros(spectrum_length)
        if cia is None:
            self.cia = np.zeros(spectrum_length)
        else:
            self.cia = cia
            self.cia_input_units = cia_input_units.lower()
            if self.cia_input_units not in self.CIA_CONVERSION:
                raise ValueError(f"Unknown x_units '{self.cia_input_units}'. Add it to CIA_CONVERSION dict.")
            self.cia *= self.CIA_CONVERSION[self.cia_input_units]

        
        #Other Spectrum Inputs
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
        self.TIPS = TIPS
        self.compressability_file = compressability_file

        #Initial Functions and Checks
        self.diluent_sum_check() # Makes sure that the diluent contributions sum to 1

   
    def _convert_y_data(self, y_raw, y_unc_raw):
        y_converted = y_raw.copy()
        y_unc_converted = y_unc_raw.copy() if y_unc_raw is not None else None

        if self.dataspace == 'alpha':
            y_output_units = {'text': '10^-6 cm^-1', 'formatted': ' (10$^{-6}$ cm$^{-1}$)'}
            y_output_label = {'text': 'alpha', 'formatted' : '$\\alpha$'}
            if self.y_input_units == 'cm-1':
                y_converted *= 1e6 
            elif self.y_input_units == 'us':
                y_converted = (CONSTANTS['c']*y_raw / 1e12)**-1

            if y_unc_raw is not None:
                if self.y_unc_input_units == 'cm-1':
                    y_unc_converted *=1e6
                elif self.y_unc_input_units == 'us':
                    y_unc_converted = (CONSTANTS['c']*y_unc_raw / 1e12)**-1
                elif self.y_unc_input_units == '%':
                    y_unc_converted *= (y_converted/100)
                elif self.y_unc_input_units == 'rel':
                    y_unc_converted *= (y_converted)            
            
        elif (self.dataspace == 'absorbance'):
            y_output_units = {'text': '', 'formatted': ''}
            y_output_label = {'text': 'absorbance', 'formatted' : 'absorbance'}
            if (self.y_input_units == 'ppm') or (self.y_input_units == '10-6'):
                y_converted /= 1e6 

            if y_unc_raw is not None:
                if (self.y_unc_input_units == 'ppm') or (self.y_unc_input_units == '10-6'):
                    y_unc_converted /=1e6
                elif self.y_unc_input_units == '%':
                    y_unc_converted *= (y_converted/100)   
                elif self.y_unc_input_units == 'rel':
                    y_unc_converted *= (y_converted)   
   
        elif (self.dataspace == 'absorption'):
            y_output_units = {'text': '', 'formatted': ''}
            y_output_label = {'text': 'absorption', 'formatted' : 'absorption'}
            if y_unc_raw is not None:
                if self.y_unc_input_units == '%':
                    y_unc_converted *= (y_converted/100)   
                elif self.y_unc_input_units == 'rel':
                    y_unc_converted *= (y_converted)   
        elif (self.dataspace == 'transmittance'):
            y_output_units = {'text': '', 'formatted': ''}
            y_output_label = {'text': 'transmittance', 'formatted' : 'transmittance'}
            if (self.y_input_units == '%'):
                y_converted /= 100
            if y_unc_raw is not None:
                if (self.y_unc_input_units == '%') & (self.y_input_units == '%'):
                    print ('Ambiguous meaning of percent in uncertainty.  Assume meaning is transmission units')
                elif (self.y_unc_input_units == "%"):
                    y_unc_converted *= (y_converted/100)
                elif self.y_unc_input_units == 'rel':
                    y_unc_converted *= (y_converted)  
        if y_unc_converted is None:
            y_unc_converted = np.zeros(len(y_converted))


        return y_converted, y_unc_converted, y_output_units, y_output_label
    

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

    def segment_wave_y(self):
        """Defines the wavenumber, alpha, and indices of spectrum that correspond to a given spectrum segment.


        Returns
        -------
        wavenumber_segments : dict
            dictionary where the key corresponds to a segment number and the values correspond to the wavenumbers for that segment
        y_segments : dict
            dictionary where the key corresponds to a segment number and the values correspond to the y-axis values for that segment.
        indices_segments : dict
            dictionary where the key corresponds to a segment number and the values correspond to the array indices for that segment.

        """

        wavenumber_segments = {}
        y_segments = {}
        indices_segments = {}
        for segment in list(set(self.segments)):
            indices = [i for i, x in enumerate(self.segments) if x == segment]
            indices_segments[segment] = indices
            wavenumber_segments[segment] = self.wavenumber[indices]
            y_segments[segment] = self.y_data[indices]
        return wavenumber_segments, y_segments, indices_segments

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
    
    def get_temperature_at_segment(self, segment_number):
        mask = (np.asarray(self.segments) == segment_number)
        if not np.any(mask):
            raise ValueError(f"Segment {segment_number} not found in this spectrum.")
        return np.mean(self.temperature_array[mask])
        
    def get_pressure_at_segment(self, segment_number):
        mask = (np.asarray(self.segments) == segment_number)
        if not np.any(mask):
            raise ValueError(f"Segment {segment_number} not found in this spectrum.")
        return np.mean(self.pressure_array[mask])

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
            self.Diluent = {self.diluent: {'composition':1, 'm': 0.0}}
                #mass will be set during HTP_wBeta_from_DF_select if necessary
        else:
            print ('If using the HTP_wBeta_from_DF_select then you need to go back and use the Diluent{diluent:{"composition": 1, "m": mass}} format')
            self.Diluent = {self.diluent: {'composition':1, 'm': 0.0}}
    def set_Diluent(self, new_Diluent):
        self.Diluent = new_Diluent
    def set_spectrum_number(self, new_spectrum_number):
        self.spectrum_number = new_spectrum_number
    

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

    def plot_wave_y(self):
        """Generates a plot of y (for dataspace) as a function of wavenumber (cm-1).
        """
        fig, ax = plt.subplots()
        ax.plot(self.wavenumber, self.y_data)
        ax.set_xlabel('Wavenumber (cm$^{-1}$)')
        ax.set_ylabel(self.y_output_label['formatted'] + self.y_output_units['formatted'])
        ax.ticklabel_format(axis='both', style='plain', useOffset=False)
       
        plt.show()


    def calculate_QF(self):
        """Calculates the quality of fit factor (QF) for a spectrum - QF = abs(maximum y - minimum y) / std(residuals).  Absolute value accounts for the sign in tranmission where the minimum is the peak height in Transmission

        Returns
        -------
        float
            QF.

        """
        return np.around(np.abs(self.y_data.max() - self.y_data.min()) / self.residuals.std(),0)

    def plot_model_residuals(self):
        """Generates a plot of the alpha and model (ppm/cm) as a function of wavenumber (cm-1) and on lower plot shows the residuals (ppm/cm) as a function of wavenumber (cm-1).
        """

        fig = plt.figure(figsize = (16,10))
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
        QF = self.calculate_QF()
        ax0 = plt.subplot(gs[0])
        ax0.plot(self.wavenumber,self.model, 'r-' )
        ax0.plot(self.wavenumber, self.y_data, 'k.')
        ax0.ticklabel_format(useOffset=False)
        ax0.text(0.25, 0.95,'QF: ' + str(QF), horizontalalignment='center', verticalalignment='center', transform = ax0.transAxes)
        ax0.set_title(str(self.spectrum_number) +': ' + self.filename)
        ax1 = plt.subplot(gs[1])
        ax1.ticklabel_format(useOffset=False)
        ax1.plot(self.wavenumber,self.residuals, "r-")
        ax1.set_xlabel('Wavenumber (cm$^{-1}$)')  

        ax0.set_ylabel(self.y_output_label['formatted'] + self.y_output_units['formatted'])
        ax1.set_ylabel('Residuals ' + + self.y_output_units['formatted'])
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
        new_file = pd.DataFrame()
        new_file['Spectrum Number'] = [self.spectrum_number]*len(self.y_data)
        new_file['Spectrum Name'] = [self.filename]*len(self.y_data)
        new_file['Wavenumber (cm-1)'] = self.wavenumber     
        new_file['Pressure (atm)'] = self.pressure_array
        new_file['Temperature (K)'] = self.temperature_array

        new_file[self.y_output_label['text'] + self.y_output_units['text']] = self.y_data
        new_file['Model' + self.y_output_units['text'] ] = self.model
        new_file['Residuals' + self.y_output_units['text']] = self.residuals

        new_file['QF'] = [self.calculate_QF()]*len(new_file)
        new_file['Background'] = self.background
        new_file['CIA (cm-1)' ] = self.cia

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


def simulate_spectrum(parameter_linelist, lineprofile = 'mHTP', numba_lineprofile = True,
                        wave_min=None, wave_max= None, wave_space=None, wavenumbers = [],  wave_error = 0.0,
                        SNR = None, baseline_terms = [0.0], 
                        dataspace = 'alpha', pathlength = 0,
                        temperature = 296, temperature_err = {'bias': 0, 'function': None, 'params': {}}, 
                        pressure = 1,   pressure_err = {'per_bias': 0, 'function': None, 'params': {}},
                        wing_cutoff = 25, wing_wavenumbers = 25, wing_method = 'wing_wavenumbers', sim_window = 5, 
                        filename = 'temp', 
                        molefraction = {}, molefraction_err = {},
                        isotope_list = ISO, 
                        natural_abundance = True, abundance_ratio_MI = {},
                        diluent = 'air', Diluent = {},
                        nominal_temperature = 296, etalons = {}, x_shift = 0.0, 
                        IntensityThreshold = 1e-30, 
                        num_segments = 1, beta_formalism = False,
                        ILS_function = None, ILS_resolution = 0.1, ILS_wing = 10, 
                        TIPS = PYTIPS2025, 
                        compressability_file = None, 
                        BIA_model = {'sw_depletion': False, 'farwing_continuum': None}, 
                        CIA_model = {'model': None, 'params': None}) : #ad hoc CIA model should be in cm-1
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
         temperature for simulation in K. The default is 296.
    temperature_err : dict, optional
        possible keys include 'bias', 'function', and 'params'. The bias indicates the absolute bias in Celsius of the temperature reading, which will be added to the input temperature. Function can be 'linear' with params 'm' and 'b' or 'sine' with parameters 'amp', 'phase', and 'phase'. These define a function that is added to both the bias and set temperature as a function of the wavenumber. Note: if 'function' key is not equal to None, then there also needs to be a params key to define the function.. The default is {'bias': 0, 'function': None, 'params': {}}.
    pressure : float, optional
        pressure for simulation in atm.  The default is 1
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
        if not at natural abundance sets the enhancement factor for each molecule and isotope in the following format {molec_id:{iso_id: enhancement, iso_id: enhancement}, . . . }. The default is {}.  The enhancement is the ratio of the new abundance to the natural abundance.
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
    TIPS : definition, optional
        selects the HAPI provided TIPS version to use for the partition function
    
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
            Diluent = {diluent: {'composition':1, 'm':0.0}}
         else:
            Diluent = {diluent: {'composition':1, 'm':0.0}}
            print ('THIS IS GOING TO BREAK WITH A DIVISION ERROR IF YOU USE THE BETA VERSION')

    if len(wavenumbers) == 0:
        wavenumbers = np.arange(wave_min, wave_max + wave_space, wave_space)

    linelist_for_sim = parameter_linelist.copy()

    linelist_for_sim = linelist_for_sim[
        (linelist_for_sim['nu'] >= np.min(wavenumbers) - sim_window) & 
        (linelist_for_sim['nu'] <= np.max(wavenumbers) + sim_window)]

    if 'nu' in linelist_for_sim.columns:
        linelist_for_sim.sort_values('nu', inplace=True)
    linelist_for_sim.reset_index(drop=True, inplace=True)

    linelist_for_sim['sw'] /= IntensityThreshold
    linelist_for_sim['sw_scale_factor'] = IntensityThreshold

    if lineprofile != 'HTP' and lineprofile != 'mHTP':
        for col in linelist_for_sim.columns:
            if col.startswith(('eta_', 'nuOptIm_', 'n_nuOptIm_')):
                linelist_for_sim[col] = 0
            if lineprofile == 'VP' or lineprofile == 'SDVP':
                if col.startswith(('nuVC_', 'n_nuVC_', 'nuOptRe_', 'n_nuOptRe_')):
                    linelist_for_sim[col] = 0
            if lineprofile == 'VP' or lineprofile == 'NGP':
                if col.startswith(('SD_gamma_', 'n_gamma2_', 'SD_delta_', 'n_delta2_')):
                    linelist_for_sim[col] = 0

    engine = Spectroscopic_model(linelist_for_sim, lineprofile=lineprofile, numba_lineprofile = numba_lineprofile,
                                 isotope_list = isotope_list, beta_formalism=beta_formalism)

    #Frequency axis
    
    wavenumbers_err = wavenumbers + wave_error*np.random.normal(loc = 0, scale =1, size = len(wavenumbers))

    baseline_coeffs = np.flip(baseline_terms)

    engine_etalon_dict = {}
    if etalons:
        for i, vals in etalons.items():
            # vals is [Amplitude, Period]
            # Generate random phase (0 to 2pi)
            phi = np.random.rand() * 2 * np.pi 
            engine_etalon_dict[i] = {'amp': vals[0], 'period': vals[1], 'phase': phi}

    # ILS Resolution: Ensure it is a list for the engine
    if not isinstance(ILS_resolution, list) and not isinstance(ILS_resolution, np.ndarray):
        engine_ils_res = [ILS_resolution]
    else:
        engine_ils_res = ILS_resolution
    
    #Temperature
    temperature_w_error = np.full(len(wavenumbers), temperature + temperature_err['bias'])
    if temperature_err['function'] == 'linear' and 'params' in temperature_err:
        temperature_w_error += temperature_err['params']['m']*(wavenumbers-np.min(wavenumbers)) + temperature_err['params']['b']
    elif temperature_err['function'] == 'sine' and 'params' in temperature_err:
         temperature_w_error += etalon((wavenumbers-np.min(wavenumbers)), temperature_err['params']['amp'], temperature_err['params']['period'], temperature_err['params']['phase'])

    #Pressure
    pressure_w_error = np.full(len(wavenumbers), pressure * (1 + pressure_err['per_bias']/100))
    if pressure_err['function'] == 'linear' and 'params' in pressure_err:
         pressure_w_error += pressure_err['params']['m']*(wavenumbers-np.min(wavenumbers)) + pressure_err['params']['b']
    elif pressure_err['function'] == 'sine' and 'params' in pressure_err:
         pressure_w_error += etalon((wavenumbers-np.min(wavenumbers)), pressure_err['params']['amp'], pressure_err['params']['period'], pressure_err['params']['phase'])

    #Molefraction
    dataset_molecule_list = list(molefraction.keys())
    molecules_in_paramlist = linelist_for_sim['molec_id'].unique()
    for m in molecules_in_paramlist:
        if m not in dataset_molecule_list: dataset_molecule_list.append(m)
    
    for molecule in dataset_molecule_list:
        if molecule not in molefraction: molefraction[molecule] = 0

    molefraction_w_error = {}
    for species in molefraction:
        err = molefraction_err.get(species, 0.0)
        molefraction_w_error[species] = molefraction[species] * (1 + err/100.0)

    # Compressibility Interpolator
    interp_comp_factor = None
    if compressability_file is not None:
        comp_factor = pd.read_csv(compressability_file + '.csv')
        pressures = np.asarray(comp_factor['Pressure (MPa)'].values*1e6/101325).astype(float)
        temperatures = np.asarray([x for x in list(comp_factor) if 'Pressure' not in x]).astype(float)
        comp_factor.drop('Pressure (MPa)', inplace=True, axis=1) 
        interp_comp_factor = RegularGridInterpolator(points = [pressures, temperatures], values = comp_factor.to_numpy())

    BIA_FW_LBL = (BIA_model['farwing_continuum'] == 'LBL')

    #CIA
    cia_config = None
    cia_array = np.zeros_like(wavenumbers)
    if CIA_model.get('model') == 'ad hoc':
        cia_array = CIA_model.get('values') #Should be in cm-1

    elif CIA_model.get('model') == 'Karman':
        cia_calc = O2_CIA_Karman_Model(band = CIA_model.get('band'))
    
        cia_config = {
                'model': 'Karman',
                'calculator': cia_calc,
                'params': CIA_model.get('parameters')}
        
        cia_array = cia_calc.calculate_cia(
            wavenumbers, 
            temperature, 
            pressure, 
            Diluent, 
            **CIA_model.get('parameters'))
        


    flat_abundance_ratios = {}
    unique_pairs = np.unique(np.column_stack((linelist_for_sim['molec_id'], linelist_for_sim['local_iso_id'])), axis=0)
    for m, i in unique_pairs:
        flat_abundance_ratios[(int(m), int(i))] = 1.0
    if not natural_abundance and abundance_ratio_MI:
        for m, iso_dict in abundance_ratio_MI.items():
            for i, val in iso_dict.items():
                if (int(m), int(i)) in flat_abundance_ratios:
                    flat_abundance_ratios[(int(m), int(i))] = val



    seg_number = np.arange(len(wavenumbers))
    seg_number = np.abs(seg_number// (len(wavenumbers)/num_segments)).astype(int)

    y_array = np.zeros_like(wavenumbers)
    final_pressure_array = np.zeros_like(wavenumbers)
    final_temp_array = np.zeros_like(wavenumbers)

    for seg in range(num_segments):
        idx = np.where(seg_number == seg)[0]
        if len(idx) == 0: continue
            
        # Segment conditions
        waves_seg = wavenumbers[idx] + x_shift
        seg_P = np.mean(pressure_w_error[idx])
        seg_T = np.mean(temperature_w_error[idx])
        
        # Store for output file
        final_pressure_array[idx] = seg_P
        final_temp_array[idx] = seg_T

        # --- CALL ENGINE (Handles LBL + Baseline + Etalon + ILS) ---
        y_seg = engine.calculate_spectrum(
            waves=waves_seg,
            T=seg_T,
            p=seg_P,
            molefraction=molefraction_w_error,
            abundance_ratios=flat_abundance_ratios,
            Diluent=Diluent,
            spectrum_number=1,
            spectrum_min=np.min(wavenumbers), # Ensure relative baseline x-axis is consistent
            
            baseline_coeffs=baseline_coeffs,
            etalon_dict=engine_etalon_dict,
            
            ILS_function=ILS_function,
            ILS_parameters=engine_ils_res,
            ILS_wing=ILS_wing,
            
            interpolated_compressability_file=interp_comp_factor,
            TIPS=TIPS,
            isotope_list=isotope_list,
            natural_abundance=natural_abundance,
            abundance_ratio_MI=abundance_ratio_MI,
            BIA_slope=BIA_model['sw_depletion'],
            BIA_FW_LBL=BIA_FW_LBL,
            cia_config=cia_config,
            IntensityThreshold=IntensityThreshold,
            wing_cutoff=wing_cutoff,
            wing_wavenumbers=wing_wavenumbers,
            wing_method=wing_method, 
            pathlength = pathlength,
            dataspace = dataspace
        )

        y_array[idx] = y_seg

    if SNR is None:
        y_noise = y_array
    else:
        # Scale noise by max signal
        y_noise = y_array + np.max(y_array) * np.random.normal(0, 1, len(y_array)) * (1/SNR)

    if dataspace == 'alpha':
        y_output_units = {'text': '10^-6 cm^-1', 'formatted': ' (10$^{-6}$ cm$^{-1}$)'}
        y_output_label = {'text': 'alpha', 'formatted' : '$\\alpha$'} 
    elif dataspace == 'absorbance':
        y_output_units = {'text': '', 'formatted': ''}
        y_output_label = {'text': 'absorbance', 'formatted' : 'absorbance'} 
    elif dataspace == 'absorption':
        y_output_units = {'text': '', 'formatted': ''}
        y_output_label = {'text': 'absorption', 'formatted' : 'absorption'} 
    elif dataspace == 'transmittance':
        y_output_units = {'text': '', 'formatted': ''}
        y_output_label = {'text': 'transmittance', 'formatted' : 'transmittance'} 

    spectrum = pd.DataFrame()
    spectrum['Segment Number'] = seg_number
    spectrum['Wavenumber (cm-1)'] = wavenumbers
    spectrum['Wavenumber + Noise (cm-1)'] = wavenumbers_err

    
    spectrum[y_output_label['text'] + y_output_units['text']] = y_array
    spectrum[y_output_label['text'] + ' + Noise' + y_output_units['text']] = y_noise
    # Calculate Noise % (Data - Clean_Model)
    # Note: alpha_array ALREADY includes Baseline+Etalon from the engine
    spectrum['CIA (cm-1)'] = cia_array
    spectrum['Noise (%)'] = 100 * (y_noise - y_array) / np.max(y_noise) if np.max(y_noise) != 0 else 0
    spectrum['Pressure (atm)'] = final_pressure_array
    spectrum['Temperature (K)'] = final_temp_array
    spectrum.to_csv(filename + '.csv', index=False)

    return Spectrum(filename, molefraction=molefraction, natural_abundance=natural_abundance, abundance_ratio_MI=abundance_ratio_MI, isotope_list=isotope_list,
                    diluent=diluent, Diluent=Diluent, 
                    spectrum_number=1, 
                    x_column = 'Wavenumber + Noise (cm-1)', x_input_units = 'cm-1', 
                    y_column = y_output_label['text'] + ' + Noise' + y_output_units['text'],  y_input_units = y_output_units['text'], 
                    y_unc_column = 'Noise (%)', y_unc_input_units = '%',
                    pressure_column = 'Pressure (atm)', pressure_input_units = 'atm',
                    temperature_column = 'Temperature (K)', temperature_input_units = 'K',
                    dataspace=dataspace, pathlength = pathlength,
                    etalons=etalons, nominal_temperature=nominal_temperature, x_shift=x_shift, 
                    baseline_order=len(baseline_terms)-1, weight=1,
                    ILS_function=ILS_function, ILS_resolution=ILS_resolution, ILS_wing=ILS_wing, 
                    TIPS=TIPS, compressability_file=compressability_file, cia = cia_array, cia_input_units = 'cm-1')



    