import numpy as np
import pandas as pd
from bisect import bisect
import os, sys
import matplotlib.pyplot as plt
sys.path.append(r'C:\Users\ema3\Documents\Python Scripts\HAPI')#Add hapi.py folder location to system path
from hapi import *
from matplotlib import gridspec

# proposed inputs

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


def HTP_from_DF_select(linelist, wave_min, wave_max,wave_space = 0.01, wing_cutoff = 50, wing_wavenumbers = 50, 
                pressure = 1, temperature = 296, concentration = {}, 
                natural_abundance = True, abundance_ratio_MI = {},  Diluent = {}, diluent = 'air', IntensityThreshold = 1e-30):
    
    #Generate X-axis for simulation
    wavenumbers = np.arange(wave_min, wave_max + wave_space, wave_space)
        
    
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
        Diluent == {diluent:1.}
        '''
        if diluent == 'air': 
            Diluent = {'air':1.}
        elif diluent == 'self':
            Diluent = {'self':1.}
        else:
            raise Exception('Unknown GammaL value: %s' % GammaL)
        '''
              
    #Iterate through lines in linelist
    for index, line in linelist.iterrows():
        #Get Base line parameters
        LineCenterDB = line['nu']
        LineIntensityDB = line['sw']
        LowerStateEnergyDB = line['elower']
        MoleculeNumberDB = line['molec_id']
        IsoNumberDB = line['local_iso_id']
   
        
        #Calculate partition function
        SigmaT = PYTIPS2017(MoleculeNumberDB,IsoNumberDB,T) #Partition Function T using TIPS2017
        SigmaTref = PYTIPS2017(MoleculeNumberDB,IsoNumberDB,Tref) #Partition Function T using TIPS2017
        
        #Calculate Line Intensity
        LineIntensity = EnvironmentDependency_Intensity(LineIntensityDB,T,Tref,SigmaT,SigmaTref, LowerStateEnergyDB,LineCenterDB)
        
        if LineIntensity < IntensityThreshold: continue
        
        #Isotopic Abundance Calculations
        abun_ratio = 1
        '''
        if (natural_abundance == True):
            abun_ratio = 1 / abundance(MoleculeNumberDB,IsoNumberDB)
        '''
        
        if ( natural_abundance == False) and abundance_ratio_MI != {}:
            abun_ratio = abundance_ratio_MI[MoleculeNumberDB][IsoNumberDB]
            
        
        
        ##Calculate Doppler Broadening
        cMassMol = 1.66053873e-27 # hapi
        m = molecularMass(MoleculeNumberDB,IsoNumberDB) * cMassMol * 1000
        GammaD = sqrt(2*cBolts*T*log(2)/m/cc**2)*LineCenterDB
        
        #Set values for parameter summation across diluents
        Gamma0 = 0.; Shift0 = 0.; Gamma2 = 0.; Shift2 = 0.; NuVC = 0.; EtaNumer = 0.; Y = 0;
        for species in Diluent:
            abun = Diluent[species]
            
            
            
            #Gamma0: pressure broadening coefficient HWHM
            Gamma0DB = line['gamma0_%s'%species]
            TempRatioPowerDB_Gamma0 =line['n_gamma0_%s'%species]
            Gamma0T = Gamma0DB*p/pref*(Tref/T)**TempRatioPowerDB_Gamma0
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
        
        
        
        #/ abundance(MoleculeNumberDB,IsoNumberDB)
        
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
                    pressure_column = 'Cavity Pressure /Torr', temperature_column = 'Cavity Temperature Side 2 /C', frequency_column = 'Total Frequency /MHz', 
                    tau_column = 'Mean tau/us', tau_stats_column = 'tau rel. std. dev./%', 
                    etalons = {}, nominal_temperature = 296):
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
        self.etalons = etalons
        self.nominal_temperature = nominal_temperature
        
        
        #Defined from contents of file
        file_contents = pd.read_csv(self.filename + '.csv')
        self.pressure = file_contents[self.pressure_column].mean() / 760
        self.temperature = file_contents[self.temperature_column].mean() + 273.15
        self.frequency = file_contents[self.frequency_column].values
        self.tau = file_contents[self.tau_column].values
        self.tau_stats = file_contents[self.tau_stats_column].values
        self.wavenumber = self.frequency*10**6 / 29979245800
        self.alpha = (self.tau*0.0299792458)**-1
        
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
        
        ax0.plot(self.wavenumber,self.model )
        ax0.plot(self.wavenumber, self.alpha, 'o')
        ax0.text(0.25, 0.95,'QF: ' + str(QF), horizontalalignment='center', verticalalignment='center', transform = ax0.transAxes)
        
        ax1 = plt.subplot(gs[1])
        ax1.plot(self.wavenumber,self.residuals, "-")
        plt.show()
    def save_spectrum_info(self, save_file = False):
        file_contents = pd.read_csv(self.filename + '.csv')       
        new_file = pd.DataFrame()
        new_file['Spectrum Number'] = [self.spectrum_number]*len(self.alpha)
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
            line['x-shift'] = 0
            for molecule in spectrum.concentration:
                line['Concentration_' + (ISO[(molecule, 1)][4])] = (spectrum.concentration[molecule])
            for i in range(0, self.baseline_order + 1):
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
            
           
def max_iter(pars, iter, resid, *args, **kws):
        if iter > 2500:
            return True
        else:
            return False
        
def etalon(x, amp, freq, phase):
    return amp*np.sin((2*np.pi * freq)*x+ phase) 
     
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
print (SPECTRA.get_number_nominal_temperatures())
'''




                
