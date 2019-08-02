from numpy import sqrt,abs,exp,pi,log,sin,cos,tan
import numpy as np

import pandas as pd
from bisect import bisect
import os, sys
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy import fftpack
from numba import jit


sys.path.append(r'C:\Users\ema3\Documents\Python Scripts\HAPI')#Add hapi.py folder location to system path
sys.path.append(r'C:\Users\ema3\Documents\Cold Cavity - O2 A Band')# set location of HAPI.py module
from hapi import EnvironmentDependency_Intensity, PYTIPS2017, molecularMass, PROFILE_HT, ISO, volumeConcentration

from matplotlib import gridspec
import pandas as pd
import qgrid

from lmfit import minimize, Parameters, report_fit, Model, Minimizer
#from lmfit.models import GaussianModel, LinearModel, Model, VoigtModel



import seaborn as sns
sns.set_style("whitegrid")
sns.set_style("ticks")
sns.set_context("poster")


# proposed inputs
#Update to take out all hapi calls and just work from the limited 

'''
linelist = pd.read_csv('HITRAN_2016_linelist.csv')
wavenumbers = np.arange(13040, 13045, 0.001)
wing_cutoff = 50 # wingcut-off will be x halfwidths or end of range which ever is later
wing_wavenumbers = 50
pressure = 1#atm
temperature = 296 #K
concentration ={ 7 : 0.2} # Concentration is a dictionary with a concentration (out of 1) with the Molecule ID as the key
natural_abundance = True # if false then need to correct calculations
abundance_ratio_MI = {7 : {1: 1, 2: 1/0.01, 3: 1/0.01}}
diluent = 'air' # can define also as dictionary that it loops through with the value being the abun for the calc of each param
Diluent = {}
'''

@jit
def HTP_from_DF_select(linelist, waves, wing_cutoff = 50, wing_wavenumbers = 50, 
                pressure = 1, temperature = 296, concentration = {}, 
                natural_abundance = True, abundance_ratio_MI = {},  Diluent = {}, diluent = 'air', IntensityThreshold = 1e-30):
    
    #Generate X-axis for simulation
    wavenumbers = waves
        
    
    #Set Omegas to X-values
    Xsect = [0]*len(wavenumbers)
    
    #define reference temperature and pressure
    Tref = 296. # K
    pref = 1. # atm
    # define actual temperature and pressure
    T = temperature # K
    p = pressure # atm
   
    mol_dens = volumeConcentration(p,T)
    
    
    #Sets-up the  Diluent currently limited to air or self
    if not Diluent:
        Diluent = {diluent:1.}
  
    #Iterate through lines in linelist
    
    #Calculate line intensity
    
    num_lines = len(linelist)
    linelist['SigmaT'] = num_lines*[0]
    linelist['SigmaTref'] = num_lines*[0]
    linelist['GammaD'] = num_lines*[0]
    linelist['m'] = num_lines*[0]
    
    for molec in linelist['molec_id'].unique():
        for iso in linelist ['local_iso_id'].unique():
            linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'SigmaT'] = PYTIPS2017(molec,iso,T)
            linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'SigmaTref'] = PYTIPS2017(molec,iso,Tref)
            linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'm'] = molecularMass(molec,iso) * 1.66053873e-27 * 1000
    
    linelist['LineIntensity'] = EnvironmentDependency_Intensity(linelist['sw'],T,Tref,linelist['SigmaT'],linelist['SigmaTref'],linelist['elower'],linelist['nu'])

    
    ##Calculate Doppler Broadening
    linelist['GammaD'] = sqrt(2*1.380648813E-16*T*log(2)/linelist['m']/2.99792458e10**2)*linelist['nu']
    #GammaD = sqrt(2*cBolts*T*log(2)/m/cc**2)*LineCenterDB
    

    for index, line in linelist.iterrows():

        #Get Base line parameters
        LineCenterDB = line['nu']
        #LineIntensityDB = line['sw']
        LineIntensity = line['LineIntensity']
        #LowerStateEnergyDB = line['elower']
        MoleculeNumberDB = line['molec_id']
        IsoNumberDB = line['local_iso_id']
        GammaD = line['GammaD']
   
        
        #Calculate partition function
        #SigmaT = PYTIPS2017(MoleculeNumberDB,IsoNumberDB,T) #Partition Function T using TIPS2017
        #SigmaTref = PYTIPS2017(MoleculeNumberDB,IsoNumberDB,Tref) #Partition Function T using TIPS2017
        
        #Calculate Line Intensity
        #LineIntensity = EnvironmentDependency_Intensity(LineIntensityDB,T,Tref,SigmaT,SigmaTref, LowerStateEnergyDB,LineCenterDB)
        
        if LineIntensity < IntensityThreshold: continue
        
        #Isotopic Abundance Calculations
        abun_ratio = 1
        
        '''
        if (natural_abundance == True):
            abun_ratio *= abundance(MoleculeNumberDB,IsoNumberDB)
        '''
        
        
        if ( natural_abundance == False) and abundance_ratio_MI != {}:
            abun_ratio = abundance_ratio_MI[MoleculeNumberDB][IsoNumberDB]

        
        #Set values for parameter summation across diluents
        Gamma0 = 0.; Shift0 = 0.; Gamma2 = 0.; Shift2 = 0.; NuVC = 0.; EtaNumer = 0.; Y = 0;
        for species in Diluent:
            abun = Diluent[species]
            
            
            
            #Gamma0: pressure broadening coefficient HWHM
            Gamma0DB = line['gamma0_%s'%species]
            TempRatioPowerDB_Gamma0 =line['n_gamma0_%s'%species]
            Gamma0T = Gamma0DB*(p/pref)*((Tref/T)**TempRatioPowerDB_Gamma0)
            Gamma0 += abun*Gamma0T
            
            #Delta0
            Shift0DB = line['delta0_%s'%species]
            deltap = line['n_delta0_%s'%species]  
            Shift0T = (Shift0DB + deltap*(T-Tref))*p/pref
            Shift0 += abun*Shift0T
            
            
            
            #Gamma2
            SDDB = line['SD_gamma_%s'%species]
            Gamma2DB = SDDB*Gamma0DB
            TempRatioPowerDB_Gamma2 =line['n_gamma2_%s'%species]
            Gamma2T = Gamma2DB*p/pref*(Tref/T)**TempRatioPowerDB_Gamma2
            Gamma2 += abun*Gamma2T
        
            #Delta 2
            SDshiftDB = line['SD_delta_%s'%species]
            Shift2DB = SDshiftDB*Shift0DB
            delta2p =line['n_delta2_%s'%species]
            Shift2T = (Shift2DB + delta2p*(T-Tref))*p/pref
            Shift2 += abun*Shift2T
                      
            #nuVC
            NuVCDB = line['nuVC_%s'%species]
            KappaDB = line['n_nuVC_%s'%species]
            NuVCT = NuVCDB*(p/pref)*(Tref/T)**(KappaDB)
            NuVC += abun*NuVCT
            
            #Eta
            EtaDB = line['eta_%s'%species]        
            EtaNumer += EtaDB*abun*(Gamma0T+1j*Shift0T)
            
            # Line mixing

            YDB = line['y_%s'%species]
            # What does temperature look like here
            Y += abun*YDB

        Eta = EtaNumer/(Gamma0 + 1j*Shift0)

        
        
        
        appx_voigt_width = 0.5346*Gamma0 + (0.2166*Gamma0**2 + GammaD**2)**0.5
        OmegaWingF = max(wing_wavenumbers,wing_cutoff*appx_voigt_width)

        
        
        
        BoundIndexLower = bisect(wavenumbers,LineCenterDB-OmegaWingF)
        BoundIndexUpper = bisect(wavenumbers,LineCenterDB+OmegaWingF)
        


        lineshape_vals_real, lineshape_vals_imag = PROFILE_HT(LineCenterDB,GammaD,Gamma0,Gamma2,Shift0,Shift2,NuVC,Eta,wavenumbers[BoundIndexLower:BoundIndexUpper])#[BoundIndexLower:BoundIndexUpper])

    
        Xsect[BoundIndexLower:BoundIndexUpper] += mol_dens  * \
                                                    concentration[MoleculeNumberDB] * abun_ratio * \
                                                    LineIntensity * (lineshape_vals_real + Y*lineshape_vals_imag) 

                                                        
    return (wavenumbers, np.asarray(Xsect))         

'''
nu, alpha = HTP_from_DF_select(linelist, wavenumbers, wing_cutoff = 50, wing_wavenumbers = 50, 
                pressure = pressure, temperature = temperature, concentration =concentration, 
                natural_abundance = natural_abundance, abundance_ratio_MI = abundance_ratio_MI, Diluent = {}, diluent = 'air') 

                                                              
plt.plot(nu, alpha)
plt.show()
'''
class Spectrum:

    def __init__(self, filename, concentration = {}, natural_abundance = True, diluent = 'air', Diluent = {}, abundance_ratio_MI = {}, spectrum_number = 1, 
                    input_freq = True, input_tau = True, 
                    pressure_column = 'Cavity Pressure /Torr', temperature_column = 'Cavity Temperature Side 2 /C', frequency_column = 'Total Frequency /MHz', 
                    tau_column = 'Mean tau/us', tau_stats_column = 'tau rel. std. dev./%', 
                    etalons = {}, nominal_temperature = 296, x_shift = 0):
        self.filename = filename
        self.concentration = concentration
        self.natural_abundance = natural_abundance
        self.abundance_ratio_MI = abundance_ratio_MI
        self.diluent = diluent
        if Diluent == {}:
            self.Diluent = {self.diluent: 1}
        else:
            self.Diluent = Diluent
        self.spectrum_number = spectrum_number
        self.pressure_column = pressure_column
        self.temperature_column = temperature_column
        self.frequency_column = frequency_column
        self.tau_column = tau_column
        self.tau_stats_column = tau_stats_column
        self.input_freq = input_freq
        self.input_tau = input_tau
        self.etalons = etalons
        self.nominal_temperature = nominal_temperature
        self.x_shift = x_shift
        
        
        #Defined from contents of file
        file_contents = pd.read_csv(self.filename + '.csv')
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
        
        
        
        self.model = len(self.alpha)*[0]
        self.residuals = self.alpha - self.model
        self.background = len(self.alpha)*[0]
        
    ## GETTERS    
    def get_filename(self):
        return self.filename
    def get_concentration(self):
        return self.concentration
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
    def set_concentration(self, new_concentration):
        self.concentration = new_concentration
   
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
        plt.plot(self.frequency, self.tau)
        plt.xlabel('Frequency (MHz)')
        plt.ylabel('$\\tau (\mu s)$')
        plt.show()
    def plot_wave_alpha(self):
        plt.plot(self.wavenumber, self.alpha)
        plt.xlabel('Wavenumber ($cm^{-1}$)')
        plt.ylabel('$\\alpha (\\frac{ppm}{cm})$')
        plt.show()
    def calculate_QF(self):
        return np.around((self.alpha.max() - self.alpha.min()) / self.residuals.std(),0)
    def plot_model_residuals(self):
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
        new_file['Background'] = self.background
        if save_file:
            new_file.to_csv(self.filename + '_saved.csv', index = False)
        return (new_file)
    def fft_spectrum(self):
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
        #plt.xlim(0, 100)
        plt.xlabel('Experimental Wavenumber ($cm^{-1}$)')
        plt.ylabel('Amplitude (ppm/cm')
        plt.show()
                

        


class Dataset:
    def __init__(self, spectra, dataset_name, baseline_order = 1):
        self.spectra = spectra
        self.dataset_name = dataset_name
        self.baseline_order = baseline_order
        self.renumber_spectra()
    def renumber_spectra(self):
        count = 1
        for spectrum in self.spectra:
            spectrum.set_spectrum_number(count)
            count+=1
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
        for spectrum in self.spectra:
            if spectrum.spectrum_number == spectrum_num:
                return (spectrum.get_filename())
        return None
    def get_spectrum_pressure(self, spectrum_num):
        for spectrum in self.spectra:
            if spectrum.spectrum_number == spectrum_num:
                return (spectrum.get_pressure_torr())
        return None
    def get_spectrum_temperature(self, spectrum_num):
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
    def get_number_nominal_temperatures(self):
        nominal_temperatures = []
        for spectrum in self.spectra:
            if spectrum.nominal_temperature not in nominal_temperatures:
                nominal_temperatures.append(spectrum.nominal_temperature)
        return len(nominal_temperatures), nominal_temperatures
            
            
    def average_QF(self):
        sum_ = 0
        for spectrum in self.spectra:
            sum_ += spectrum.calculate_QF()
        return sum_ / self.get_number_spectra()
    def get_list_spectrum_numbers(self):
        spec_num_list = []
        for spectrum in self.spectra:
            spec_num_list.append(spectrum.spectrum_number)
        return spec_num_list            
    def generate_baseline_paramlist(self):
        baseline_paramlist = pd.DataFrame()
        
        for spectrum in self.spectra:
            line = {}
            line['Spectrum Number'] = spectrum.spectrum_number
            line['x_shift'] = spectrum.x_shift
            for molecule in spectrum.concentration:
                line['Concentration_' + (ISO[(molecule, 1)][4])] = (spectrum.concentration[molecule])
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
        
        summary_file = pd.DataFrame()
        for spectrum in self.spectra:
            spectrum_data = spectrum.save_spectrum_info(save_file = False)
            summary_file = summary_file.append(spectrum_data)
        if save_file:
            summary_file.to_csv(self.dataset_name + '.csv', index = False)
        return summary_file
    def plot_model_residuals(self):
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
            ax1.plot(spectrum.wavenumber,spectrum.residuals, plot_color+".")
        ax0.legend(bbox_to_anchor=(1, 1))
        
        plt.show()
            
            
           
def max_iter(pars, iter, resid, *args, **kws):
        if iter > 2500:
            return True
        else:
            return False
            
        
def etalon(x, amp, freq, phase):
    return amp*np.sin((2*np.pi * freq)*x+ phase) 
    
    
class Generate_FitParam_File:
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
    # Do I need getters and setters or will this just be re-run if needed
    def get_dataset(self):
        return self.dataset
    def get_param_linelist(self):
        return self.param_linelist
    def get_base_linelist(self):
        return self.base_linelist
    def generate_fit_param_linelist_from_linelist(self, vary_nu = {7:{1:True, 2:False, 3:False}, 1:{1:False}}, vary_sw = {7:{1:True, 2:False, 3:False}},
                                   vary_gamma0 = {7:{1: True, 2:False, 3: False}, 1:{1:False}}, vary_n_gamma0 = {}, 
                                   vary_delta0 = {7:{1: True, 2:False, 3: False}, 1:{1:False}}, vary_n_delta0 = {}, 
                                   vary_aw = {7:{1: True, 2:False, 3: False}, 1:{1:False}}, vary_n_gamma2 = {}, 
                                   vary_as = {}, vary_n_delta2 = {}, 
                                   vary_nuVC = {}, vary_n_nuVC = {},
                                   vary_eta = {}, vary_linemixing = {}):
        param_linelist_df = self.get_param_linelist().copy()
        #Intensity Thresholding
        param_linelist_df = param_linelist_df[param_linelist_df['sw'] > self.threshold_intensity] #intensity thresholding
        #Cutdown linelist to parameters within specified simulation window of minimum and maximum of dataset
        dataset_min, dataset_max = (self.dataset.get_spectra_extremes())
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
                            param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'sw_' + str(spec) + '_vary'] = (vary_sw[molecule][isotope])
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
                                param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'gamma0_' +diluent +'_'+str(spec) + '_vary'] = (vary_gamma0[molecule][isotope])
            
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
                                param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'delta0_' +diluent +'_'+str(spec) + '_vary'] = (vary_delta0[molecule][isotope])
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
                                    param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'SD_gamma_' +diluent +'_'+str(spec) + '_vary'] = (vary_aw[molecule][isotope])
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
                                    param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'SD_delta_' +diluent +'_'+str(spec) + '_vary'] = (vary_as[molecule][isotope])           
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
                                    param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'nuVC_' +diluent +'_'+str(spec) + '_vary'] = (vary_nuVC[molecule][isotope])
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
                                    param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'eta_' +diluent +'_'+str(spec) + '_vary'] = (vary_eta[molecule][isotope])
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
                                        param_linelist_df.loc[(param_linelist_df['nu'] >= dataset_min)&(param_linelist_df['nu'] <= dataset_max)&(param_linelist_df['sw'] > 1) &(param_linelist_df['molec_id'] == molecule) & (param_linelist_df['local_iso_id'] == isotope), 'y_' +diluent+ '_'+ str(nominal_temp) +'_'+str(spec) + '_vary'] = (vary_linemixing[molecule][isotope])


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
                
    
        #param_linelist_df = param_linelist_df.reindex(sorted(param_linelist_df.columns), axis=1)
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
        
        #print (list(param_linelist_df))
            
        
        param_linelist_df = param_linelist_df[ordered_list]
        param_linelist_df.to_csv(self.param_linelist_savename + '.csv')
        return param_linelist_df 
 

    def generate_fit_baseline_linelist(self, vary_baseline = True, vary_concentration = {7:True, 1:False}, vary_xshift = False, 
                                      vary_etalon_amp= False, vary_etalon_freq= False, vary_etalon_phase= False):
        base_linelist_df = self.get_base_linelist().copy()
        parameters =  (list(base_linelist_df))
        for param in parameters:
            base_linelist_df[param + '_err'] = len(base_linelist_df)*[0]
            if 'baseline' in param:
                base_linelist_df[param + '_vary'] = len(base_linelist_df)*[vary_baseline]
            if 'Concentration' in param:
                for molecule in vary_concentration:
                    if (ISO[(molecule, 1)][4]) in param:
                        base_linelist_df[param + '_vary'] = len(base_linelist_df)*[(vary_concentration[molecule])]
                        base_linelist_df[param + '_err'] = len(base_linelist_df)*[0]
            if 'x_shift' in param:
                base_linelist_df[param + '_vary'] = len(base_linelist_df)*[(vary_xshift)]
                base_linelist_df[param + '_err'] = len(base_linelist_df)*[0]
            if 'amp' in param:
                base_linelist_df[param + '_vary'] = len(base_linelist_df)*[(vary_etalon_amp)]
                base_linelist_df[param + '_err'] = len(base_linelist_df)*[0]
            if 'freq' in param:
                base_linelist_df[param + '_vary'] = len(base_linelist_df)*[(vary_etalon_freq)]
                base_linelist_df[param + '_err'] = len(base_linelist_df)*[0]
            if 'phase' in param:
                base_linelist_df[param + '_vary'] = len(base_linelist_df)*[(vary_etalon_phase)]  
                base_linelist_df[param + '_err'] = len(base_linelist_df)*[0] 
        base_linelist_df = base_linelist_df.reindex(sorted(base_linelist_df.columns), axis=1)
        base_linelist_df.to_csv(self.base_linelist_savename + '.csv')
        return base_linelist_df

    
class Edit_Fit_Param_Files:
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
        base_linelist_df = pd.read_csv(self.base_linelist_file + '.csv', index_col = 0)
        baseline_widget = qgrid.show_grid(base_linelist_df, grid_options={'forceFitColumns': False, 'defaultColumnWidth': 200})
        return baseline_widget
    def save_editted_baselist(self, baseline_widget):
        base_linelist_df = baseline_widget.get_changed_df()
        base_linelist_df.to_csv(self.new_base_linelist_file + '.csv')
        return base_linelist_df
    def edit_generated_paramlist(self):
        param_linelist_df = pd.read_csv(self.param_linelist_file + '.csv', index_col = 0)
        param_widget = qgrid.show_grid(param_linelist_df, grid_options={'forceFitColumns': False, 'defaultColumnWidth': 200})
        return param_widget
    def save_editted_paramlist(self, param_widget):
        param_linelist_df = param_widget.get_changed_df()
        param_linelist_df.to_csv(self.new_param_linelist_file + '.csv') 
        return param_linelist_df
     
def hasNumbers(inputString):
    has_digits = False
    for char in inputString:
        if char.isdigit():
            return True
    return False

class Fit_DataSet:
    def __init__(self, dataset, base_linelist_file, param_linelist_file, fit_intensity = 1e-27,
                baseline_limit = False, baseline_limit_factor = 10, 
                concentration_limit = False, concentration_limit_factor = 10, 
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
        self.baseline_list = pd.read_csv(self.base_linelist_file + '.csv', index_col = 0)
        
        self.param_linelist_file = param_linelist_file
        self.lineparam_list = pd.read_csv(self.param_linelist_file + '.csv', index_col = 0)
        
        self.fit_intensity = fit_intensity
        self.baseline_limit = baseline_limit
        self.baseline_limit_factor  = baseline_limit_factor
        self.etalon_limit = etalon_limit
        self.etalon_limit_factor = etalon_limit_factor
        self.concentration_limit = concentration_limit
        self.concentration_limit_factor = concentration_limit_factor
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
        params = Parameters()
        
        #Baseline Parameters
        #baseline_list = pd.read_csv(self.base_linelist_file + '.csv', index_col = 0)
        baseline_parameters = []
        for base_param in list(self.baseline_list):
            if ('_vary' not in base_param) and ('_err' not in base_param):
                baseline_parameters.append(base_param)
        for spec_num in self.baseline_list.index.values:
            for base_param in baseline_parameters:
                #print (base_param)
                if self.baseline_list.loc[int(spec_num)][base_param] == 0:
                    params.add(base_param + '_'+str(int(spec_num)), self.baseline_list.loc[int(spec_num)][base_param], self.baseline_list.loc[int(spec_num)][base_param + '_vary'])
                elif ('Concentration' in base_param) and self.concentration_limit:
                    params.add(base_param + '_'+str(int(spec_num)), self.baseline_list.loc[int(spec_num)][base_param], self.baseline_list.loc[int(spec_num)][base_param + '_vary'], 
                              min = (1 / self.concentration_limit_factor)*self.baseline_list.loc[int(spec_num)][base_param], 
                              max = self.concentration_limit_factor*self.baseline_list.loc[int(spec_num)][base_param])
                elif ('baseline' in base_param) and self.baseline_limit:
                    params.add(base_param + '_'+str(int(spec_num)), self.baseline_list.loc[int(spec_num)][base_param], self.baseline_list.loc[int(spec_num)][base_param + '_vary'], 
                              min = (1 / self.baseline_limit_factor)*self.baseline_list.loc[int(spec_num)][base_param], 
                              max = self.baseline_limit_factor *self.baseline_list.loc[int(spec_num)][base_param])
                elif ('etalon_' in base_param) and self.etalon_limit and ('phase' not in base_param):
                    params.add(base_param + '_'+str(int(spec_num)), self.baseline_list.loc[int(spec_num)][base_param], self.baseline_list.loc[int(spec_num)][base_param + '_vary'], 
                              min = (1 / self.etalon_limit_factor )*self.baseline_list.loc[int(spec_num)][base_param], 
                              max = self.etalon_limit_factor *self.baseline_list.loc[int(spec_num)][base_param])
                elif ('etalon_' in base_param) and self.etalon_limit and ('phase' in base_param):
                    params.add(base_param + '_'+str(int(spec_num)), self.baseline_list.loc[int(spec_num)][base_param], self.baseline_list.loc[int(spec_num)][base_param + '_vary'], 
                              min = -2*np.pi, max = 2*np.pi)
                elif ('x_shift' in base_param) and self.x_shift_limit:
                    params.add(base_param + '_'+str(int(spec_num)), self.baseline_list.loc[int(spec_num)][base_param], self.baseline_list.loc[int(spec_num)][base_param + '_vary'], 
                              min = (self.baseline_list.loc[int(spec_num)][base_param] - self.x_shift_limit_magnitude), 
                              max = self.x_shift_limit_magnitude + self.baseline_list.loc[int(spec_num)][base_param])   
                else:
                    params.add(base_param + '_'+str(int(spec_num)), self.baseline_list.loc[int(spec_num)][base_param], self.baseline_list.loc[int(spec_num)][base_param + '_vary'])
        
        #Lineshape parameters
        #lineparam_list = pd.read_csv(self.param_linelist_file + '.csv', index_col = 0)
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
        
        if ((sum('nu' in param for param in linelist_params)) - (sum('nuVC' in param for param in linelist_params))) > 1:
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
            if (sum('gamma0' in param for param in linelist_params)) >  len(diluent_list) + 1:
                gamma0_constrain = False
            if (sum('delta0' in param for param in linelist_params)) >  len(diluent_list) + 1:
                delta0_constrain = False
            if (sum('SD_gamma' in param for param in linelist_params)) >  len(diluent_list) + 1:
                SD_gamma_constrain = False
            if (sum('SD_delta' in param for param in linelist_params)) >  len(diluent_list) + 1:
                SD_delta_constrain = False
            if (sum('nuVC' in param for param in linelist_params)) >  len(diluent_list) + 1:
                nuVC_constrain = False
            if (sum('eta' in param for param in linelist_params)) >  len(diluent_list) + 1:
                eta_constrain = False 
        if (sum('y' in param for param in linelist_params)) >  len(diluent_list)*(self.dataset.get_number_nominal_temperatures()[0]):
            linemix_constrain = False
        linemix_terms_constrained = []
        for diluent in diluent_list:
            for temperature in (self.dataset.get_number_nominal_temperatures()[1]):
                linemix_terms_constrained.append('y_'+diluent + '_' + str(temperature))
        
        for spec_line in self.lineparam_list.index.values:
            if self.lineparam_list.loc[spec_line]['sw'] >= self.fit_intensity / self.lineparam_list.loc[spec_line]['sw_scale_factor']:# bigger than 1 because fit_intensity / fit_intensity
                for line_param in linelist_params:
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
                    elif ('gamma0_' in line_param) and ('n_' not in line_param) and (gamma0_constrain) and (not hasNumbers(line_param.replace("gamma0_", ''))):
                        if self.gamma0_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.gamma0_limit_factor)*self.lineparam_list.loc[int(spec_line)][line_param] ,#(1 - (self.gamma0_limit_percent / 100))*lineparam_list.loc[int(spec_line)][line_param], 
                                  max = self.gamma0_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])#(1 + (self.gamma0_limit_percent / 100))*lineparam_list.loc[int(spec_line)][line_param])
                                
                        else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    elif ('gamma0_' in line_param) and ('n_' not in line_param) and (not gamma0_constrain) and (hasNumbers(line_param.replace("gamma0_", ''))):
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
                    elif ('delta0' in line_param) and ('n_' not in line_param) and (delta0_constrain) and (not hasNumbers(line_param.replace("delta0", ''))):
                        if self.delta0_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.delta0_limit_factor )*self.lineparam_list.loc[int(spec_line)][line_param], 
                                  max = self.delta0_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    elif ('delta0_' in line_param) and ('n_' not in line_param) and (not delta0_constrain) and (hasNumbers(line_param.replace("delta0_", ''))):
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
                    elif ('SD_gamma' in line_param) and (SD_gamma_constrain) and (not hasNumbers(line_param.replace("SD_gamma", ''))):
                        if self.SD_gamma_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.SD_gamma_limit_factor) *self.lineparam_list.loc[int(spec_line)][line_param], 
                                  max = self.SD_gamma_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'],
                            min = 1e-4)
                    elif ('SD_gamma' in line_param) and (not SD_gamma_constrain) and (hasNumbers(line_param.replace("SD_gamma", ''))):
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
                    elif ('SD_delta' in line_param) and (SD_delta_constrain) and (not hasNumbers(line_param.replace("SD_delta", ''))):
                        if self.SD_delta_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.SD_delta_limit_factor )*self.lineparam_list.loc[int(spec_line)][line_param], 
                                  max = self.SD_delta_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    elif ('SD_delta' in line_param) and (not SD_delta_constrain) and (hasNumbers(line_param.replace("SD_delta", ''))):
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
                    elif ('nuVC' in line_param) and ('n_' not in line_param) and (nuVC_constrain) and (not hasNumbers(line_param.replace("nuVC", ''))):
                        if self.nuVC_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 /self.nuVC_limit_factor)*self.lineparam_list.loc[int(spec_line)][line_param], 
                                  max = self.nuVC_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    elif ('nuVC_' in line_param) and ('n_' not in line_param) and (not nuVC_constrain) and (hasNumbers(line_param.replace("nuVC_", ''))):
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
                    elif ('eta_' in line_param) and (eta_constrain) and (not hasNumbers(line_param.replace("eta", ''))):
                        if self.eta_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.eta_limit_factor)*self.lineparam_list.loc[int(spec_line)][line_param], 
                                  max = (self.eta_limit_factor)*self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    elif ('eta_' in line_param) and (not eta_constrain) and (hasNumbers(line_param.replace("eta_", ''))):
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
    def add_expr_params(self, params, param_name, expr):
        pass
    def simulation_model(self, params, wing_cutoff = 50, wing_wavenumbers = 50):
        total_simulated = []
        total_residuals = []
        baseline_params = []
        linelist_params = []
        
        # Set-up Baseline Parameters
        for param in (list(params.valuesdict().keys())):
            if ('Concentration' in param) or ('baseline' in param) or ('etalon' in param) or ('x_shift' in param):
                baseline_params.append(param)
            else:
                linelist_params.append(param)
        #Lineparameter fitting initialization       
        #lineparam_list = pd.read_csv(self.param_linelist_file + '.csv', index_col = 0)
        for spectrum in self.dataset.spectra:
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
            #print (linelist_for_sim)
            
            #Set-up Concentration for Fitting
            fit_concentration = spectrum.concentration
            for molecule in spectrum.concentration:
                if ('Concentration_'+ ISO[(molecule, 1)][4]) + '_' + str(spectrum_number) in baseline_params:
                    fit_concentration[molecule] = np.float(params[('Concentration_'+ ISO[(molecule, 1)][4]) + '_' + str(spectrum_number)])
            
            #Simulate Spectra

            x_shift = np.float(params['x_shift_' + str(spectrum_number)])
            linelist_for_sim['nu'] = linelist_for_sim['nu'] + x_shift
            
            fit_nu, fit_coef = HTP_from_DF_select(linelist_for_sim, spectrum.wavenumber, wing_cutoff = wing_cutoff, wing_wavenumbers = wing_wavenumbers, 
                    pressure = spectrum.pressure, temperature = spectrum.temperature, concentration = fit_concentration, 
                    natural_abundance = spectrum.natural_abundance, abundance_ratio_MI = spectrum.abundance_ratio_MI,  Diluent = Diluent)
            fit_coef *= 1e6

            
            
            
            ## Baseline Calculation
            baseline_param_array = [0]*(self.dataset.baseline_order+1)
            for param in baseline_params:
                if ('baseline' in param) and ((spectrum_number) == int(param[11:])):
                    baseline_param_array[ord(param[9:param.find('_',9)])-97] = np.float(params[param])
            
            baseline_param_array = baseline_param_array[::-1] # reverses array to be used for polyval
            baseline = np.polyval(baseline_param_array, spectrum.wavenumber - np.min(spectrum.wavenumber))
            
            
            #Etalon Calculation
            #print ((spectrum.etalons))
            fit_etalon_parameters = {}
            for i in range(1, len(spectrum.etalons)+1):
                fit_etalon_parameters[i] = {'amp': 0, 'freq':1, 'phase':0}
            
            for param in baseline_params:
                if ('etalon' in param) and (str(spectrum_number) in param[param.find('_', 7):]):
                    etalon_num = int(param[param.find('_')+1: param.find('_', param.find('_')+1)])
                    if ('amp' in param) and (str(etalon_num) in param):
                        fit_etalon_parameters[etalon_num]['amp'] = np.float(params[param])
                    if ('freq' in param) and (str(etalon_num) in param):
                        fit_etalon_parameters[etalon_num]['freq'] = np.float(params[param])
                    if ('phase' in param) and (str(etalon_num) in param):
                        fit_etalon_parameters[etalon_num]['phase'] = np.float(params[param])
            etalons = len(spectrum.wavenumber)*[0]
            for i in range(1, len(spectrum.etalons)+1):
                etalons += etalon(spectrum.wavenumber - np.min(spectrum.wavenumber), fit_etalon_parameters[i]['amp'], fit_etalon_parameters[i]['freq'], fit_etalon_parameters[i]['phase'])
                
            simulated_spectra = baseline + etalons + fit_coef
            residuals = simulated_spectra - spectrum.alpha
       
            total_simulated = np.append(total_simulated, simulated_spectra)
            total_residuals = np.append(total_residuals, residuals)
        
        total_residuals = np.asarray(total_residuals)
        total_simulated = np.asarray(total_simulated)
        return total_residuals
    def fit_data(self, params, wing_cutoff = 50, wing_wavenumbers = 50):
        minner = Minimizer(self.simulation_model, params, iter_cb = max_iter, fcn_args=(wing_cutoff, wing_wavenumbers))
        result = minner.minimize(method = 'leastsq')#'
        return result
    def residual_analysis(self, result, indv_resid_plot = False):
        residual_array = result.residual
        for spectrum in self.dataset.spectra:
            spectrum_residual, residual_array = np.split(residual_array, [len(spectrum.wavenumber)])
            spectrum.set_residuals(spectrum_residual)
            spectrum.set_model(spectrum_residual + spectrum.alpha)
            if indv_resid_plot:
                spectrum.plot_model_residuals()
    def update_params(self, result, base_linelist_update_file = None , param_linelist_update_file = None):
        if base_linelist_update_file == None:
            base_linelist_update_file = self.base_linelist_file
        if param_linelist_update_file == None:
            param_linelist_update_file = self.param_linelist_file
        #baseline_list = pd.read_csv(self.base_linelist_file + '.csv', index_col = 0) # from the input
        #lineparam_list = pd.read_csv(self.param_linelist_file + '.csv', index_col = 0)
        
        
        for key, par in result.params.items():
            if ('Concentration' in par.name) or ('baseline' in par.name) or ('x_shift' in par.name):
                parameter = par.name[:(par.name.find('_', par.name.find('_')+1))]
                line = int(par.name[(par.name.find('_', par.name.find('_')+1))+1:])
                self.baseline_list.loc[line, parameter] = par.value
                if par.vary:
                    self.baseline_list.loc[line, parameter + '_err'] = par.stderr
            elif ('etalon' in par.name):
                parameter =  (par.name[:par.name.find('_',par.name.find('_', par.name.find('_')+1)+1)])
                line=  int(par.name[par.name.find('_',par.name.find('_', par.name.find('_')+1)+1)+1:])
                self.baseline_list.loc[line, parameter] = par.value
                if par.vary:
                    self.baseline_list.loc[line, parameter + '_err'] = par.stderr
            else:
                parameter = par.name[:par.name.find('_line')]
                line = int(par.name[par.name.find('_line')+6:])
                self.lineparam_list.loc[line, parameter] = par.value
                if par.vary:
                    self.lineparam_list.loc[line, parameter + '_err'] = par.stderr
        self.baseline_list.to_csv(base_linelist_update_file + '.csv')
        self.lineparam_list.to_csv(param_linelist_update_file + '.csv')
        
        #Calculate Baseline + Etalons and add to the Baseline term for each spectra
        for spectrum in self.dataset.spectra:
            baseline = len(spectrum.wavenumber)*[0]
            
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
            baseline += np.polyval(baseline_param_array, spectrum.wavenumber - np.min(spectrum.wavenumber))
            for i in range(1, len(spectrum.etalons)+1):
                baseline += etalon(spectrum.wavenumber - np.min(spectrum.wavenumber), fit_etalon_parameters[i]['amp'], fit_etalon_parameters[i]['freq'], fit_etalon_parameters[i]['phase'])
            spectrum.set_background(baseline)

                     

                
            

        
'''
spectrum_80 = Spectrum('15Feb19_80Torr Air_ oxygen Computer _new config_corr', 
                       concentration = { 7 : 0.2072}, natural_abundance = True, diluent = 'air', 
                       etalons = {1: [0.001455, 0.83629422058], 2: [ 6.2325e-04, 1.67258844116]})
spectrum_100 = Spectrum('15Feb19_100Torr Air_ oxygen Computer _new config_corr', 
                        concentration = { 7 : 0.2072}, natural_abundance = True, diluent = 'air', 
                        etalons = {1: [0.001455, 0.83629422058], 2: [ 6.2325e-04, 1.67258844116]})
spectrum_120 = Spectrum('15Feb19_120Torr Air_ Oxygen Computer _new config_corr', 
                        concentration = { 7 : 0.2072}, natural_abundance = True, diluent = 'air',
                        etalons = {1: [0.001455, 0.83629422058], 2: [ 6.2325e-04, 1.67258844116]})

#Add all spectrum to a Dataset object
SPECTRA = Dataset([spectrum_80, spectrum_100, spectrum_120], 'New Configuration Oxygen A-Band')  
print (SPECTRA.plot_model_residuals())
'''





                
