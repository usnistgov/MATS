#Import Packages
import numpy as np
import pandas as pd
from bisect import bisect
#import sys
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
#Files that should be in the same folder
from hapi import PYTIPS2017, molecularMass, pcqsdhc, ISO
from Karman_CIA import Karman_CIA_Model


#Constants
h = 6.62607015e-27 #erg s https://physics.nist.gov/cgi-bin/cuu/Value?h|search_for=h as of 5/21/2020
c = 29979245800 #cm/s # https://physics.nist.gov/cgi-bin/cuu/Value?c|search_for=c as of 5/21/2020
k = 1.380649e-16 # erg / K https://physics.nist.gov/cgi-bin/cuu/Value?k as of 5/21/2020     
Na = 6.02214076e23 # mol-1 https://physics.nist.gov/cgi-bin/cuu/Value?na as of 5/21/2020
cpa_atm = (10*101325)**-1 #convert from cpa to atm  https://physics.nist.gov/cgi-bin/cuu/Value?stdatm|search_for=atmosphere as of 5/21/2020
c2 =  (h*c)/k


def HTP_from_DF_select(linelist, waves, wing_cutoff = 50, wing_wavenumbers = 50, wing_method = 'wing_cutoff',
                p = 1, T = 296, molefraction = {}, 
                natural_abundance = True, abundance_ratio_MI = {},  Diluent = {}, diluent = 'air', IntensityThreshold = 1e-30):
    
    #Generate X-axis for simulation
    wavenumbers = waves
        
    #Set Omegas to X-values
    Xsect = [0]*len(wavenumbers)
    
    #define reference temperature/pressure and calculate molecular density
    Tref = 296. # K
    pref = 1. # atm
    mol_dens = (p/cpa_atm)/(k*T) 
    
    #Sets-up the  Diluent (currently limited to air or self, unless manual input in Diluent)
    #Sets-up the  Diluent (currently limited to air or self, unless manual input in Diluent)
    if not Diluent:
         #Diluent = {diluent:1.}
         if diluent == 'air':
             Diluent = {diluent: {'composition':1, 'm':28.95734}}
         elif diluent == 'self':
            Diluent = {diluent: {'composition':1, 'm':0}}

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
            linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'm'] = molecularMass(molec,iso) #* 1.66053873e-27 * 1000 #cmassmol and kg conversion 
            if ( natural_abundance == False) and abundance_ratio_MI != {}:
                linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'abun_ratio'] = abundance_ratio_MI[molec][iso]
    
    #linelist['LineIntensity'] = EnvironmentDependency_Intensity(linelist['sw'],T,Tref,linelist['SigmaT'],linelist['SigmaTref'],linelist['elower'],linelist['nu'])

    linelist['LineIntensity'] = linelist['sw']*linelist['SigmaTref']/linelist['SigmaT']*(np.exp(-c2*linelist['elower']/T)*(1-np.exp(-c2*linelist['nu']/T)))/(np.exp(-c2*linelist['elower']/Tref)*(1-np.exp(-c2*linelist['nu']/Tref)))

    
    #Calculate Doppler Broadening
    linelist['GammaD'] = np.sqrt(2*k*Na*T*np.log(2)/(linelist['m'].values))*linelist['nu'] / c
    # Calculated Line Parameters across Broadeners
    linelist['Gamma0'] = 0
    linelist['Shift0'] = 0
    linelist['Gamma2'] = 0
    linelist['Shift2'] = 0
    linelist['NuVC'] = 0
    linelist['Eta'] = 0
    linelist['Y'] = 0
    for species in Diluent:
        abun = Diluent[species]['composition']
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
        linelist['Eta'] += linelist['eta_%s'%species] *abun
        #Line mixing
        linelist['Y'] += abun*(linelist['y_%s'%species] *(p/pref))
    
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

def HTP_wBeta_from_DF_select(linelist, waves, wing_cutoff = 50, wing_wavenumbers = 50, wing_method = 'wing_cutoff',
                p = 1, T = 296, molefraction = {}, 
                natural_abundance = True, abundance_ratio_MI = {},  Diluent = {}, diluent = 'air', IntensityThreshold = 1e-30):
    
    #Generate X-axis for simulation
    wavenumbers = waves
        
    #Set Omegas to X-values
    Xsect = [0]*len(wavenumbers)
    
    #define reference temperature/pressure and calculate molecular density
    Tref = 296. # K
    pref = 1. # atm
    mol_dens = (p/cpa_atm)/(k*T) 
    
    #Sets-up the  Diluent (currently limited to air or self, unless manual input in Diluent)
    if not Diluent:
         #Diluent = {diluent:1.}
         if diluent == 'air':
             Diluent = {diluent: {'composition':1, 'm':28.95734}}
         elif diluent == 'self':
            Diluent = {diluent: {'composition':1, 'm':0}}
         else:
            Diluent = {diluent: {'composition':1, 'm':0}}
            print ('THIS IS GOING TO BREAK WITH A DIVISION ERROR.  GO BACK AND USE DILUENT FORMALISM FOR THE SPECTRUM DEFINITION')
            
            

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
            #linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'ma'] = molecularMass(molec,iso)
            linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'm'] = molecularMass(molec,iso) #* 1.66053873e-27 * 1000 #cmassmol and kg conversion
            if (len(Diluent) == 1) & ('self' in Diluent):
                Diluent['self']['mp'] = molecularMass(molec,iso)
            if ( natural_abundance == False) and abundance_ratio_MI != {}:
                linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'abun_ratio'] = abundance_ratio_MI[molec][iso]
    # Calculate mp
    mp = 0
    for diluent in Diluent:
        mp += Diluent[diluent]['composition']*Diluent[diluent]['m']
     
    # Get Line Intensity
    linelist['LineIntensity'] = linelist['sw']*linelist['SigmaTref']/linelist['SigmaT']*(np.exp(-c2*linelist['elower']/T)*(1-np.exp(-c2*linelist['nu']/T)))/(np.exp(-c2*linelist['elower']/Tref)*(1-np.exp(-c2*linelist['nu']/Tref)))

    #Calculate Doppler Broadening
    #linelist['GammaD'] = np.sqrt(2*1.380648813E-16*T*np.log(2)/linelist['m']/2.99792458e10**2)*linelist['nu']
    linelist['GammaD'] = np.sqrt(2*k*Na*T*np.log(2)/(linelist['m'].values))*linelist['nu'] / c
    
    # Calculated Line Parameters across Broadeners
    linelist['Gamma0'] = 0
    linelist['Shift0'] = 0
    linelist['Gamma2'] = 0
    linelist['Shift2'] = 0
    linelist['NuVC'] = 0
    linelist['Eta'] = 0
    linelist['Y'] = 0
    for species in Diluent:
        abun = Diluent[species]['composition']
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
        linelist['Eta'] += linelist['eta_%s'%species] *abun
        #Line mixing
        linelist['Y'] += abun*(linelist['y_%s'%species]*(p/pref))
        
    linelist['alpha'] =  mp / linelist['m']
    linelist['Chi'] = linelist['NuVC'] / linelist['GammaD']
    linelist['A'] = 0.0534 + 0.1585*np.exp(-0.451*linelist['alpha'].values)
    linelist['B'] = 1.9595 - 0.1258*linelist['alpha'].values + 0.0056*linelist['alpha'].values**2 + 0.0050*linelist['alpha'].values**3
    linelist['C'] = -0.0546 + 0.0672*linelist['alpha'].values - 0.0125*linelist['alpha'].values**2+0.0003*linelist['alpha'].values**3
    linelist['D'] = 0.9466 - 0.1585*np.exp(-0.4510*linelist['alpha'].values)
    linelist['Beta'] = linelist['A'].values*np.tanh(linelist['B'].values * np.log10(linelist['Chi'].values) + linelist['C'].values) + linelist['D'].values
   
    
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
                                                    line['NuVC']*line['Beta'],line['Eta'],wavenumbers[BoundIndexLower:BoundIndexUpper])    
        Xsect[BoundIndexLower:BoundIndexUpper] += mol_dens  * \
                                                    molefraction[line['molec_id']] * line['abun_ratio'] * \
                                                    line['LineIntensity'] * (lineshape_vals_real + line['Y']*lineshape_vals_imag) 
    
    # Return two arrays corresponding to the wavenumber axis and the calculated cross-section                                                                                                      
    return (wavenumbers, np.asarray(Xsect)) 

class Spectrum:
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
        self.cia = len(self.alpha)*[0]
    
    def diluent_sum_check(self):
        diluent_sum = 0
        for dil in self.Diluent:
            diluent_sum+=self.Diluent[dil]['composition']
        if diluent_sum != 1:
            print ("YOUR DILUENTS DO NOT SUM TO ONE!  They sum to " + str(diluent_sum))

    def segment_wave_alpha(self):
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
        self.alpha = (self.tau*c*1e-10)**-1
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
    def set_cia(self, new_cia):
        self.cia = new_cia
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
        new_file['QF'] = [self.calculate_QF()]*len(new_file)
        new_file['Background'] = self.background
        new_file['CIA (ppm/cm)'] = self.cia 
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

class Dataset:
    def __init__(self, spectra, dataset_name, baseline_order = 1, CIA_model = None):
        self.spectra = spectra
        self.dataset_name = dataset_name
        self.baseline_order = baseline_order
        self.CIA_model = CIA_model
        self.renumber_spectra()
        self.molecule_list = self.correct_component_list()
        self.correct_etalon_list()
        self.max_baseline_order()
        self.broadener_list = self.get_broadener_list()
        
    def renumber_spectra(self):
        count = 1
        for spectrum in self.spectra:
            spectrum.set_spectrum_number(count)
            count+=1

    def max_baseline_order(self):
        baseline_order_list = []
        for spectrum in self.spectra:
            baseline_order_list.append(spectrum.baseline_order)
        self.baseline_order = max(baseline_order_list)

    def correct_component_list(self):
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
        return dataset_molecule_list
    
    def get_broadener_list(self):
        dataset_broadener_list = []
        for spectrum in self.spectra:
            dataset_broadener_list += spectrum.Diluent.keys()
        dataset_broadener_list = list(set(dataset_broadener_list))
        return dataset_broadener_list
            


    def correct_etalon_list(self):
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
        dataset_etalon_list = []
        for spectrum in self.spectra:
            dataset_etalon_list += spectrum.etalons.keys()
        dataset_etalon_list = list(set(dataset_etalon_list))
        return dataset_etalon_list

    def get_molecules(self):
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

    def get_spectrum_extremes(self):
        extreme_dictionary = {}
        for spectrum in self.spectra:
            extreme_dictionary[spectrum.get_spectrum_number()] = [np.min(spectrum.wavenumber), np.max(spectrum.wavenumber)]
        return extreme_dictionary
            
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
            ax1.plot(spectrum.wavenumber,spectrum.residuals, plot_color+"-")
        ax0.legend(bbox_to_anchor=(1, 1))
        plt.show()          
           
def max_iter(pars, iter, resid, *args, **kws):
        if iter > 2500:
            return True
        else:
            return False
                   
def etalon(x, amp, period, phase):
    return amp*np.sin((2*np.pi * period)*x+ phase)   

   
        
def simulate_spectrum(parameter_linelist, wave_min, wave_max, wave_space, wave_error = 0, 
                        SNR = None, baseline_terms = [0], temperature = 25, temperature_err = {'bias': 0, 'function': None, 'params': {}}, pressure = 760, 
                        pressure_err = {'per_bias': 0, 'function': None, 'params': {}}, 
                        wing_cutoff = 25, wing_wavenumbers = 25, wing_method = 'wing_cutoff', filename = 'temp', molefraction = {}, molefraction_err = {},
                        natural_abundance = True, abundance_ratio_MI = {},diluent = 'air', Diluent = {}, 
                        nominal_temperature = 296, etalons = {}, x_shift = 0, IntensityThreshold = 1e-30, num_segments = 1, beta_formalism = False):
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
                                p = segment_pressure, T = segment_temperature,  molefraction = molefraction_w_error, 
                                natural_abundance = natural_abundance, abundance_ratio_MI = abundance_ratio_MI,  
                                Diluent = Diluent, diluent = diluent, IntensityThreshold = IntensityThreshold)
        else:
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
        period = etalons[r][1]
        phase = np.random.rand()
        x = wavenumbers - np.min(wavenumbers)
        etalon_model += amp*np.sin((2*np.pi * period)*x+ phase) 
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
    def __init__ (self, dataset, param_linelist, base_linelist, CIA_linelist = None, 
                  lineprofile = 'VP', linemixing = False, threshold_intensity = 1e-30, fit_intensity = 1e-26, fit_window = 1.5, sim_window = 5, 
                  param_linelist_savename = 'Parameter_LineList', base_linelist_savename = 'Baseline_LineList', CIA_linelist_savename = 'CIA_LineList', 
                 nu_constrain = True, sw_constrain = True, gamma0_constrain = True, delta0_constrain = True, aw_constrain = True, as_constrain = True, 
                 nuVC_constrain = True, eta_constrain =True, linemixing_constrain = True):
        self.dataset = dataset
        self.param_linelist = param_linelist
        self.base_linelist = base_linelist
        if self.dataset.CIA_model == None:
            self.CIA_linelist = None
            self.CIA_linelist_savename = None
        else:
            self.CIA_linelist = CIA_linelist
            self.CIA_linelist_savename = CIA_linelist_savename

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
    def get_CIA_linelist(self):
        return self.CIA_linelist
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
                            param_linelist_df['y_' +diluent + '_'+ str(nominal_temp)+'_' +str(spec)] = (param_linelist_df['y_' + diluent+ '_' +str(nominal_temp)].values)
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
        param_linelist_df.to_csv(self.param_linelist_savename + '.csv') #
        
        #Warnings of SD_gamma, SD_delta, nuVC, eta are zero and floated
        for param in ordered_list:
            if ('SD_gamma' in param) & ('vary' not in param) & ('err' not in param):
                if len(param_linelist_df[(param_linelist_df[param]==0) & (param_linelist_df[param + '_vary']==True)])!=0:
                    print (param, ': floating SD_gamma terms when the initial guess is zero can make it difficult for the solution to converge. Consider setting initial guess to a non-zero value.')
            if ('SD_delta' in param) & ('vary' not in param) & ('err' not in param):
                if len(param_linelist_df[(param_linelist_df[param]==0) & (param_linelist_df[param + '_vary']==True)])!=0:
                    print (param, ': floating SD_delta terms when the initial guess is zero can make it difficult for the solution to converge. Consider setting initial guess to a non-zero value.')
            if ('nuVC' in param) & ('vary' not in param) & ('err' not in param) & ('n_nuVC' not in param):
                if len(param_linelist_df[(param_linelist_df[param]==0) & (param_linelist_df[param + '_vary']==True)])!=0:
                    print (param, ': floating nuVC terms when the initial guess is zero can make it difficult for the solution to converge. Consider setting initial guess to a non-zero value.  This is generally less of an issue than the case where SD_gamma is floated when the initial guess is zero.')
                    
        return param_linelist_df 

    def generate_fit_baseline_linelist(self, vary_baseline = True, vary_pressure = False, vary_temperature = False,vary_molefraction = {7:True, 1:False}, vary_xshift = False, 
                                      vary_etalon_amp= False, vary_etalon_period= False, vary_etalon_phase= False):
        base_linelist_df = self.get_base_linelist().copy()
        parameters =  (list(base_linelist_df))

        #Generate Fit Baseline file
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
            if 'period' in param:
                base_linelist_df.loc[base_linelist_df[param]!=0, param + '_vary'] = (vary_etalon_period) 
            if 'phase' in param:
                base_linelist_df.loc[base_linelist_df[param.replace("phase", "period")]!=0, param + '_vary'] = (vary_etalon_phase)
        base_linelist_df.drop(['Baseline Order'], axis=1, inplace = True)
        #base_linelist_df = base_linelist_df.reindex(sorted(base_linelist_df.columns), axis=1)
        base_linelist_df.to_csv(self.base_linelist_savename + '.csv', index = True)
        return base_linelist_df
        
    def generate_fit_KarmanCIA_linelist(self, vary_EXCH_scalar = False, vary_EXCH_gamma = False, vary_EXCH_l = False, 
                                  vary_SO_scalar = False, vary_SO_ahard = False, vary_SO_l = False, vary_bandcenter = False, vary_Nmax = False, wave_step = 5):
        CIA_linelist_df = self.get_CIA_linelist().copy()
        parameters =  (list(CIA_linelist_df))
        #Initial Calculation of the Karman CIA based on initial guesses and application to spectra
        if self.dataset.CIA_model == 'Karman':
            CIA_pairs = CIA_linelist_df['CIA Pair'].unique()
            for spectrum in self.dataset.spectra:
                CIA = len(spectrum.wavenumber)*[0]
                for CIA_pair in CIA_pairs:
                    for molecule in self.dataset.molecule_list:
                        molecule_name = ISO[(molecule, 1)][4]
                        for broadener in self.dataset.broadener_list:
                            try:
                                broadener_composition = spectrum.Diluent[broadener]['composition']
                            except:
                                broadener_composition = 0
                            if broadener == 'self':
                                broadener = molecule_name
                            if (molecule_name in CIA_pair) and (broadener in CIA_pair):
                                CIA_select = CIA_linelist_df[CIA_linelist_df['CIA Pair'] == CIA_pair]
                                CIA += broadener_composition*Karman_CIA_Model(spectrum.wavenumber, spectrum.get_pressure_torr(), spectrum.get_temperature_C(), wave_step = wave_step,
                                     EXCH_scalar = CIA_select['EXCH_scalar'].values[0], EXCH_gamma = CIA_select['EXCH_gamma'].values[0], EXCH_l = CIA_select['EXCH_l'].values[0],
                                     SO_scalar = CIA_select['SO_scalar'].values[0], SO_ahard = CIA_select['SO_ahard'].values[0], SO_l = CIA_select['SO_l'].values[0],
                                     bandcenter = CIA_select['bandcenter'].values[0], Nmax = CIA_select['Nmax'].values[0])                  
                spectrum.set_cia(CIA)

            # Set parameter floats
            for param in parameters:
                if ('CIA Pair' != param):
                    CIA_linelist_df[param + '_err'] = 0
                    CIA_linelist_df[param + '_vary']= False
                if 'EXCH_scalar' in param:
                    CIA_linelist_df[param + '_vary']= len(CIA_linelist_df)*[(vary_EXCH_scalar)]
                if 'EXCH_gamma' in param:
                    CIA_linelist_df[param + '_vary']= len(CIA_linelist_df)*[(vary_EXCH_gamma)]
                    if vary_EXCH_gamma:
                        print ('USE CAUTION WHEN FLOATING EXCH_GAMMA')
                if 'EXCH_l' in param:
                    CIA_linelist_df[param + '_vary']= len(CIA_linelist_df)*[(vary_EXCH_l)]
                    if vary_EXCH_l:
                        print ('USE CAUTION WHEN FLOATING EXCH_L')
                if 'SO_scalar' in param:
                    CIA_linelist_df[param + '_vary']= len(CIA_linelist_df)*[(vary_SO_scalar)]
                if 'SO_ahard' in param:
                    CIA_linelist_df[param + '_vary']= len(CIA_linelist_df)*[(vary_SO_ahard)]
                    if vary_SO_ahard:
                        print ('USE CAUTION WHEN FLOATING SO_AHARD')
                if 'SO_l' in param:
                    CIA_linelist_df[param + '_vary']= len(CIA_linelist_df)*[(vary_SO_l)] 
                    if vary_SO_l:
                        print ('USE CAUTION WHEN FLOATING SO_L')
                if 'bandcenter' in param:
                    CIA_linelist_df[param + '_vary']= len(CIA_linelist_df)*[(vary_bandcenter)]   
                if 'Nmax' in param:
                    CIA_linelist_df[param + '_vary']= len(CIA_linelist_df)*[(vary_Nmax)] 
                    if vary_Nmax:
                        print ('USE CAUTION WHEN FLOATING NMAX')
            CIA_linelist_df = CIA_linelist_df[['CIA Pair', 
                      'EXCH_scalar', 'EXCH_scalar_err', 'EXCH_scalar_vary', 
                      'EXCH_gamma', 'EXCH_gamma_err', 'EXCH_gamma_vary',
                      'EXCH_l', 'EXCH_l_err', 'EXCH_l_vary', 
                      'SO_scalar', 'SO_scalar_err', 'SO_scalar_vary', 
                      'SO_ahard', 'SO_ahard_err', 'SO_ahard_vary',
                      'SO_l', 'SO_l_err','SO_l_vary',
                      'bandcenter','bandcenter_err', 'bandcenter_vary', 
                      'Nmax', 'Nmax_err', 'Nmax_vary']]
            CIA_linelist_df.to_csv(self.CIA_linelist_savename + '.csv', index = False)
            return CIA_linelist_df
        else:
            print ('Generate_fit_KarmanCIA_linelist only applies to the CIA files generated for the Karman CIA model')
            return None
            

    
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
        base_linelist_df.to_csv(self.new_base_linelist_file + '.csv', index = False)
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
    for char in inputString:
        if char.isdigit():
            return True
    return False

class Fit_DataSet:
    def __init__(self, dataset, base_linelist_file, param_linelist_file, CIA_linelist_file = None,
                minimum_parameter_fit_intensity = 1e-27,
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
                eta_limit = True, eta_limit_factor  = 10, 
                linemixing_limit = False, linemixing_limit_factor  = 10, 
                beta_formalism = False, 
                CIA_calculate = False, CIA_model = None, CIA_wavestep = 5):
        self.dataset = dataset
        self.base_linelist_file = base_linelist_file
        self.baseline_list = pd.read_csv(self.base_linelist_file + '.csv')#, index_col = 0
        self.param_linelist_file = param_linelist_file
        self.lineparam_list = pd.read_csv(self.param_linelist_file + '.csv', index_col = 0)
        self.CIA_linelist_file = CIA_linelist_file
        if self.CIA_linelist_file != None:
            self.CIA_linelist_file = CIA_linelist_file
            self.CIAparam_list = pd.read_csv(self.CIA_linelist_file + '.csv')
        else:
            self.CIAparam_list = None
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
        self.beta_formalism = beta_formalism
        self.CIA_calculate = CIA_calculate
        self.CIA_model = CIA_model
        self.CIA_wavestep = CIA_wavestep
        
    def generate_params(self):
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
        # CIA parameters
        if self.CIA_calculate & (self.CIA_linelist_file != None):
            if (self.CIA_model != None):
                CIA_parameters = []
                for CIA_param in list(self.CIAparam_list):
                    if ('_vary' not in CIA_param) and ('_err' not in CIA_param) and ('CIA Pair' not in CIA_param):
                        CIA_parameters.append(CIA_param)
                for index in self.CIAparam_list.index.values:
                    CIA_pair = self.CIAparam_list.iloc[index]['CIA Pair']
                    for CIA_param in CIA_parameters:
                        params.add(CIA_param + '_'+ CIA_pair, self.CIAparam_list.loc[index][CIA_param], self.CIAparam_list.loc[index][CIA_param + '_vary'])

            
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
                                  min = (1 / self.gamma0_limit_factor)*self.lineparam_list.loc[int(spec_line)][line_param] ,
                                  max = self.gamma0_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])
                                
                        else:
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    elif ('gamma0_' in line_param) and ('n_' not in line_param) and (not gamma0_constrain) and (index_length>1):
                        if self.gamma0_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                             params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.gamma0_limit_factor)*self.lineparam_list.loc[int(spec_line)][line_param]
                                  max = self.gamma0_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])
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
                            params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
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
                            if self.beta_formalism:
                                params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], min = 0)
                            else:
                                params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'])
                    elif ('nuVC' in line_param) and ('n_nuVC' not in line_param) and (not nuVC_constrain) and (index_length>1):
                        if self.nuVC_limit and self.lineparam_list.loc[spec_line][line_param] != 0:
                             params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], 
                                  min = (1 / self.nuVC_limit_factor)*self.lineparam_list.loc[int(spec_line)][line_param], 
                                  max = self.nuVC_limit_factor*self.lineparam_list.loc[int(spec_line)][line_param])
                        else:
                            if self.beta_formalism:
                                params.add(line_param + '_' + 'line_' + str(spec_line), self.lineparam_list.loc[spec_line][line_param], self.lineparam_list.loc[spec_line][line_param + '_vary'], min = 0)
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
                                    etalon_amp_segment_constrained = True, etalon_period_segment_constrained = True, etalon_phase_segment_constrained = True, 
                                    pressure_segment_constrained = True, temperature_segment_constrained = True):
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

    def simulation_model(self, params, wing_cutoff = 50, wing_wavenumbers = 50, wing_method = 'wing_cutoff'):
        total_simulated = []
        total_residuals = []
        baseline_params = []
        linelist_params = []
        CIA_params = []
        
        # Set-up Baseline Parameters
        for param in (list(params.valuesdict().keys())):
            if ('molefraction' in param) or ('baseline' in param) or ('etalon' in param) or ('x_shift' in param) or ('Pressure' in param) or ('Temperature' in param):
                baseline_params.append(param)
            elif ('EXCH_scalar' in param) or ('EXCH_gamma' in param) or ('EXCH_l' in param) or ('SO_scalar' in param) or ('SO_ahard' in param) or ('SO_ahard' in param) or ('SO_l' in param) or ('bandcenter' in param) or ('Nmax' in param):
                CIA_params.append(param)
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
                
                if self.beta_formalism == True:
                    fit_nu, fit_coef = HTP_wBeta_from_DF_select(linelist_for_sim, wavenumbers, wing_cutoff = wing_cutoff, wing_wavenumbers = wing_wavenumbers, wing_method = wing_method,
                            p = p, T = T, molefraction = fit_molefraction, 
                            natural_abundance = spectrum.natural_abundance, abundance_ratio_MI = spectrum.abundance_ratio_MI,  Diluent = Diluent)
                else:
                    fit_nu, fit_coef = HTP_from_DF_select(linelist_for_sim, wavenumbers, wing_cutoff = wing_cutoff, wing_wavenumbers = wing_wavenumbers, wing_method = wing_method,
                            p = p, T = T, molefraction = fit_molefraction, 
                            natural_abundance = spectrum.natural_abundance, abundance_ratio_MI = spectrum.abundance_ratio_MI,  Diluent = Diluent)
                fit_coef = fit_coef * 1e6
                ## CIA Calculation
                if self.CIA_calculate and self.CIA_model == 'Karman':
                    CIA_pairs = []
                    for param in CIA_params:
                        if 'EXCH_scalar' in param:
                            CIA_pairs.append(param.replace('EXCH_scalar_', ''))
                    CIA = len(wavenumbers)*[0]
                    for CIA_pair in CIA_pairs:
                        for molecule in self.dataset.molecule_list:
                            molecule_name = ISO[(molecule, 1)][4]
                            for broadener in self.dataset.broadener_list:
                                try:
                                    broadener_composition = spectrum.Diluent[broadener]['composition']
                                except:
                                    broadener_composition = 0
                                if broadener == 'self':
                                    broadener = molecule_name
                            if (molecule_name in CIA_pair) and (broadener in CIA_pair):
                                CIA += broadener_composition*Karman_CIA_Model(wavenumbers, p*760, T-273.15, wave_step = self.CIA_wavestep,
                                     EXCH_scalar = params['EXCH_scalar' + '_' + CIA_pair], EXCH_gamma = params['EXCH_gamma' + '_' + CIA_pair], EXCH_l = params['EXCH_l' + '_' + CIA_pair],
                                     SO_scalar = params['SO_scalar' + '_' + CIA_pair], SO_ahard = params['SO_ahard' + '_' + CIA_pair], SO_l = params['SO_l' + '_' + CIA_pair],
                                     bandcenter = params['bandcenter' + '_' + CIA_pair], Nmax = params['Nmax' + '_' + CIA_pair])
                else:
                    CIA = spectrum.cia[np.min(indices_segments[segment]): np.max(indices_segments[segment])+1]

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
                    fit_etalon_parameters[i] = {'amp': 0, 'period':1, 'phase':0}
                for param in baseline_params:
                    if ('etalon' in param) and (str(spectrum_number) in param[param.find('_', 7):]):
                        etalon_num = int(param[param.find('_')+1: param.find('_', param.find('_')+1)])
                        if param == 'etalon_' + str(etalon_num) + '_amp_' + str(spectrum_number) + '_' +str(segment):#('amp' in param) and (str(etalon_num) in param):
                            fit_etalon_parameters[etalon_num]['amp'] = np.float(params[param])
                        if param == 'etalon_' + str(etalon_num) + '_period_' + str(spectrum_number) + '_' +str(segment):#('period' in param) and (str(etalon_num) in param):
                            fit_etalon_parameters[etalon_num]['period'] = np.float(params[param])
                        if param == 'etalon_' + str(etalon_num) + '_phase_' + str(spectrum_number) +'_' + str(segment):#('phase' in param) and (str(etalon_num) in param):
                            fit_etalon_parameters[etalon_num]['phase'] = np.float(params[param])
                etalons = len(wavenumbers)*[0]
                for i in range(1, len(spectrum.etalons)+1):
                    etalons += etalon(wavenumbers_relative, fit_etalon_parameters[i]['amp'], fit_etalon_parameters[i]['period'], fit_etalon_parameters[i]['phase']) 
                simulated_spectra[np.min(indices_segments[segment]): np.max(indices_segments[segment])+1] = (baseline + etalons + fit_coef + CIA)
                residuals[np.min(indices_segments[segment]): np.max(indices_segments[segment])+1]  = (baseline + etalons + fit_coef + CIA) - alpha_segments[segment]
            total_simulated = np.append(total_simulated, simulated_spectra)
            total_residuals = np.append(total_residuals, residuals)
        total_residuals = np.asarray(total_residuals)
        total_simulated = np.asarray(total_simulated)
        return total_residuals
    def fit_data(self, params, wing_cutoff = 50, wing_wavenumbers = 50, wing_method = 'wing_cutoff', xtol = 1e-7, maxfev = 2000, ftol = 1e-7):
        minner = Minimizer(self.simulation_model, params, xtol =xtol, maxfev =  maxfev, ftol = ftol, fcn_args=(wing_cutoff, wing_wavenumbers, wing_method))
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
    def update_params(self, result, base_linelist_update_file = None , param_linelist_update_file = None, CIA_linelist_update_file = None):
        if base_linelist_update_file == None:
            base_linelist_update_file = self.base_linelist_file
        if param_linelist_update_file == None:
            param_linelist_update_file = self.param_linelist_file
        if CIA_linelist_update_file == None:
            CIA_linelist_update_file = self.CIA_linelist_file
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
            elif ('EXCH_scalar' in par.name) or ('EXCH_gamma' in par.name) or ('EXCH_l' in par.name) or ('SO_scalar' in par.name) or ('SO_ahard' in par.name) or ('SO_l' in par.name):
                if self.CIA_calculate and self.CIA_model == 'Karman':
                    indices = [m.start() for m in re.finditer('_', par.name)]
                    parameter = (par.name[:indices[1]])
                    CIA_pair = (par.name[indices[1] + 1:])
                    self.CIAparam_list.loc[(self.CIAparam_list['CIA Pair'] == CIA_pair), parameter] = par.value
                    if par.vary:
                        self.CIAparam_list.loc[(self.CIAparam_list['CIA Pair'] == CIA_pair), parameter + '_err'] = par.stderr
            elif ('bandcenter' in par.name) or ('Nmax' in par.name):
                if self.CIA_calculate and self.CIA_model == 'Karman':
                    indices = [m.start() for m in re.finditer('_', par.name)]
                    parameter = (par.name[:indices[0]])
                    CIA_pair = (par.name[indices[0] + 1:indices[1]])
                    self.CIAparam_list.loc[(self.CIAparam_list['CIA Pair'] == CIA_pair), parameter] = par.value
                    if par.vary:
                        self.CIAparam_list.loc[(self.CIAparam_list['CIA Pair'] == CIA_pair), parameter + '_err'] = par.stderr
            else:
                parameter = par.name[:par.name.find('_line')]
                line = int(par.name[par.name.find('_line')+6:])
                self.lineparam_list.loc[line, parameter] = par.value
                if par.vary:
                    self.lineparam_list.loc[line, parameter + '_err'] = par.stderr
        self.baseline_list.to_csv(base_linelist_update_file + '.csv', index = False)
        if CIA_linelist_update_file != None:
            self.CIAparam_list.to_csv(CIA_linelist_update_file + '.csv', index = False)
        self.lineparam_list.to_csv(param_linelist_update_file + '.csv')
        #Calculated CIA and add to the CIA term for each spectra
        if self.CIA_calculate and self.CIA_model == 'Karman':
            for spectrum in self.dataset.spectra:
                CIA_pairs = self.CIAparam_list['CIA Pair'].unique()
                for spectrum in self.dataset.spectra:
                    CIA = len(spectrum.wavenumber)*[0]
                    for CIA_pair in CIA_pairs:
                        for molecule in self.dataset.molecule_list:
                            molecule_name = ISO[(molecule, 1)][4]
                            for broadener in self.dataset.broadener_list:
                                try:
                                    broadener_composition = spectrum.Diluent[broadener]['composition']
                                except:
                                    broadener_composition = 0
                                if broadener == 'self':
                                    broadener = molecule_name
                                if (molecule_name in CIA_pair) and (broadener in CIA_pair):
                                    CIA_select = self.CIAparam_list[self.CIAparam_list['CIA Pair'] == CIA_pair]
                                    CIA += broadener_composition*Karman_CIA_Model(spectrum.wavenumber, spectrum.get_pressure_torr(), spectrum.get_temperature_C(), wave_step = self.CIA_wavestep,
                                         EXCH_scalar = CIA_select['EXCH_scalar'].values[0], EXCH_gamma = CIA_select['EXCH_gamma'].values[0], EXCH_l = CIA_select['EXCH_l'].values[0],
                                         SO_scalar = CIA_select['SO_scalar'].values[0], SO_ahard = CIA_select['SO_ahard'].values[0], SO_l = CIA_select['SO_l'].values[0],
                                         bandcenter = CIA_select['bandcenter'].values[0], Nmax = CIA_select['Nmax'].values[0])                  
                    spectrum.set_cia(CIA)
                
        
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
                    if ('baseline' in par.name) and (str(spectrum.spectrum_number) in par.name):
                        baseline_param_array[ord(par.name[9:par.name.find('_',9)])-97] = np.float(par.value)
                    elif ('etalon' in par.name) and (str(spectrum.spectrum_number) in par.name[par.name.find('_', 7):]):
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
    def generate_beta_output_file(self, beta_summary_filename = None ):
        if beta_summary_filename == None:
            beta_summary_filename = 'Beta Summary File'
        if self.beta_formalism == True:
            beta_summary_list = lineparam_list.copy()
            #List of all Diluents in Dataset
            diluent_list = []
            for spectrum in self.dataset.spectra:
                for diluent in spectrum.Diluent:
                    if diluent not in diluent_list:
                        diluent_list.append(diluent)
            #Determine if line center and nuVC terms are constrained across Dataset or vary for each spectrum
            nu_constrain = True
            nuVC_constrain = True

            
            if ((sum(('nu' in param) & ('nuVC' not in param) for param in linelist_params))) > 1:
                nu_constrain = False
            if (self.dataset.get_number_nominal_temperatures()[0]) == 1:
                if (sum('nuVC' in param for param in linelist_params)) >  len(diluent_list):
                    nuVC_constrain = False
            else:
                if (sum('nuVC' in param for param in linelist_params)) >  2*len(diluent_list) :
                    nuVC_constrain = False
                        
            #Add Column for mass
            for molec in beta_summary_list['molec_id'].unique():
                for iso in beta_summary_list ['local_iso_id'].unique():
                    beta_summary_list.loc[(beta_summary_list['molec_id']==molec) & (beta_summary_list['local_iso_id']==iso), 'm'] = molecularMass(molec,iso) 

            #Single or MS for nu and nuVC
            for spectrum in self.dataset.spectra:
                wavenumber_segments, alpha_segments, indices_segments = spectrum.segment_wave_alpha()
                mp = 0
                for diluent in spectrum.Diluent:
                    mp += spectrum.Diluent[diluent]['composition']*spectrum.Diluent[diluent]['m']
                
                for segment in list(set(spectrum.segments)):
                    p = self.baseline_list[(baseline_list['Spectrum Number'] == spectrum.spectrum_number) & (baseline_list['Segment Number'] == segment))]['Pressure'].values[0]
                    T = self.baseline_list[(baseline_list['Spectrum Number'] == spectrum.spectrum_number) & (baseline_list['Segment Number'] == segment))]['Temperature'].values[0]
                    wave_min = np.min(wavenumber_segments[segment])
                    wave_max = np.max(wavenumber_segments[segment])

                    alpha = mp / beta_summary_list['m']
                    if nu_constrain:
                        GammaD = np.sqrt(2*k*Na*T*np.log(2)/(beta_summary_list['m'].values))*beta_summary_list['nu'] / c #change with nu
                    else:
                        GammaD = np.sqrt(2*k*Na*T*np.log(2)/(beta_summary_list['m'].values))*beta_summary_list['nu' + '_' + str(spectrum.spectrum_number)] / c #change with nu
                    nuVC = len(GammaD)*[0]
                    for diluent in spectrum.Diluent:
                        abun = spectrum.Diluent[species]['composition']
                        if nuVC_constrain:
                            nuVC += abun*(beta_summary_list['nuVC_%s'%diluent]*(p/1)*((296/T)**(linelist['n_nuVC_%s'%diluent]))) 
                        else:
                            nuVC += abun*(beta_summary_list['nuVC_%s'%diluent]*(p/1)*((296/T)**(linelist['n_nuVC_%s_%s'%(diluent,str(spectrum.spectrum_number))]))) 
                    Chi = nuVC/ GammaD
                    A = 0.0534 + 0.1585*np.exp(-0.451*linelist['alpha'].values)
                    B = 1.9595 - 0.1258*linelist['alpha'].values + 0.0056*linelist['alpha'].values**2 + 0.0050*linelist['alpha'].values**3
                    C = -0.0546 + 0.0672*linelist['alpha'].values - 0.0125*linelist['alpha'].values**2+0.0003*linelist['alpha'].values**3
                    D = 0.9466 - 0.1585*np.exp(-0.4510*linelist['alpha'].values)
                    beta_summary_list['beta'] =  A*np.tanh(B * np.log10(Chi) + C) + D
                    beta_summary_list.loc[(beta_summary_list['nu'] >= wave_min) & (beta_summary_list['nu'] <= wave_max),'Beta_' + str(spectrum.spectrum_number) ] =beta_summary_list['(beta_summary_list['nu'] >= wave_min) & (beta_summary_list['nu'] <= wave_max)']['beta'].values
                    

            select_columns = ['molec_id', 'local_iso_id', 'nu']
            for param in beta_summary_list.columns:
                if ('nuVC' in param) & ('n_nuVC' not in param):
                    select_columns.append(param)
                if ('Beta' in param):
                    select_columns.append(param)
            beta_summary_list = beta_summary_list[select_columns]
            beta_summary_list.to_csv(beta_summary_filename + '.csv')


       