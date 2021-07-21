#Import Packages
# from .Utilities import *
from bisect import bisect
import re

import numpy as np
import pandas as pd
from .hapi import ISO, PYTIPS2017, pcqsdhc
from .utilities import molecularMass, etalon, convolveSpectrumSame
from .codata import CONSTANTS

from lmfit import Minimizer,  Parameters


def HTP_from_DF_select(linelist, waves, wing_cutoff = 25, wing_wavenumbers = 25, wing_method = 'wing_cutoff',
                p = 1, T = 296, molefraction = {}, isotope_list = ISO,
                natural_abundance = True, abundance_ratio_MI = {},  Diluent = {}, diluent = 'air', IntensityThreshold = 1e-30):
    """Calculates the absorbance (ppm/cm) based on input line list, wavenumbers, and spectrum environmental parameters.

    Outline

    1.  Uses provided wavenumber axis

    2.  Calculates the molecular density from pressure and temperature

    3.  Sets up Diluent dictionary if not given as input

    4.  Calculates line intensity and doppler width at temperature for all lines

    5.  Loops through each line in the line list and loops through each diluent, generating a line parameter at experimental conditions
        that is the appropriate ratio of each diluent species corrected for pressure and temperature.  For each line, simulate the line for the given simulation cutoffs and add to cross section

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
    isotope_list : dict, optional
        provides opportunity to specify the isotope look-up table.  Default is ISO, which is from HAPI.  If not using ISO, then must use this format and suggested you use function to add to ISO
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
    mol_dens = (p/ CONSTANTS['cpa_atm'])/(CONSTANTS['k']*T)

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
            try:
                linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'SigmaT'] = PYTIPS2017(molec,iso,T)
            except:
                linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'SigmaT'] = 1
            try:
                linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'SigmaTref'] = PYTIPS2017(molec,iso,Tref)
            except:
                linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'SigmaTref'] = 1
            linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'm'] = molecularMass(molec,iso, isotope_list = isotope_list) #* 1.66053873e-27 * 1000 #cmassmol and kg conversion
            if ( natural_abundance == False) and abundance_ratio_MI != {}:
                linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'abun_ratio'] = abundance_ratio_MI[molec][iso]

    #linelist['LineIntensity'] = EnvironmentDependency_Intensity(linelist['sw'],T,Tref,linelist['SigmaT'],linelist['SigmaTref'],linelist['elower'],linelist['nu'])

    linelist['LineIntensity'] = linelist['sw']*linelist['SigmaTref']/linelist['SigmaT']*(np.exp(-CONSTANTS['c2']*linelist['elower']/T)*(1-np.exp(-CONSTANTS['c2']*linelist['nu']/T)))/(np.exp(-CONSTANTS['c2']*linelist['elower']/Tref)*(1-np.exp(-CONSTANTS['c2']*linelist['nu']/Tref)))
    if isotope_list != ISO:
        linelist.loc[(linelist['SigmaT'] == 1) & (linelist['SigmaTref'] == 1), 'LineIntensity'] = linelist[(linelist['SigmaT'] == 1) & (linelist['SigmaTref'] == 1)]['sw'].values

    #Calculate Doppler Broadening
    linelist['GammaD'] = np.sqrt(2*CONSTANTS['k']*CONSTANTS['Na']*T*np.log(2)/(linelist['m'].values))*linelist['nu'] / CONSTANTS['c']
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

def HTP_wBeta_from_DF_select(linelist, waves, wing_cutoff = 25, wing_wavenumbers = 25, wing_method = 'wing_cutoff',
                p = 1, T = 296, molefraction = {}, isotope_list = ISO,
                natural_abundance = True, abundance_ratio_MI = {},  Diluent = {}, diluent = 'air', IntensityThreshold = 1e-30):
    """Calculates the absorbance (ppm/cm) based on input line list, wavenumbers, and spectrum environmental parameters with capability of incorporating the beta correction to the Dicke Narrowing proposed in Analytical-function correction to the Hartmann–Tran profile for more reliable representation of the Dicke-narrowed molecular spectra.

    Outline

    1.  Uses provided wavenumber axis

    2.  Calculates the molecular density from pressure and temperature

    3.  Sets up Diluent dictionary if not given as input

    4.  Calculates line intensity and doppler width at temperature for all lines

    5.  Loops through each line in the line list and loops through each diluent, generating a line parameter at experimental conditions
        that is the appropriate ratio of each diluent species corrected for pressure and temperature.  For each line, simulate the line for the given simulation cutoffs and add to cross section

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
    isotope_list : dict, optional
        provides opportunity to specify the isotope look-up table.  Default is ISO, which is from HAPI.  If not using ISO, then must use this format and suggested you use function to add to ISO
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
    mol_dens = (p/CONSTANTS['cpa_atm'])/(CONSTANTS['k']*T)

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
            try:
                linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'SigmaT'] = PYTIPS2017(molec,iso,T)
            except:
                linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'SigmaT'] = 1
            try:
                linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'SigmaTref'] = PYTIPS2017(molec,iso,Tref)
            except:
                linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'SigmaTref'] = 1


            linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'm'] = molecularMass(molec,iso, isotope_list = isotope_list) #* 1.66053873e-27 * 1000 #cmassmol and kg conversion
            if (len(Diluent) == 1) & ('self' in Diluent):
                Diluent['self']['mp'] = molecularMass(molec,iso, isotope_list = isotope_list)
            if ( natural_abundance == False) and abundance_ratio_MI != {}:
                linelist.loc[(linelist['molec_id']==molec) & (linelist['local_iso_id']==iso), 'abun_ratio'] = abundance_ratio_MI[molec][iso]
    # Calculate mp
    mp = 0
    for diluent in Diluent:
        mp += Diluent[diluent]['composition']*Diluent[diluent]['m']

    # Get Line Intensity
    linelist['LineIntensity'] = linelist['sw']*linelist['SigmaTref']/linelist['SigmaT']*(np.exp(-CONSTANTS['c2']*linelist['elower']/T)*(1-np.exp(-CONSTANTS['c2']*linelist['nu']/T)))/(np.exp(-CONSTANTS['c2']*linelist['elower']/Tref)*(1-np.exp(-CONSTANTS['c2']*linelist['nu']/Tref)))
    if isotope_list != ISO:
        linelist.loc[(linelist['SigmaT' == 1]) & (linelist['SigmaTref' == 1])] = linelist[(linelist['SigmaT' == 1]) & (linelist['SigmaTref' == 1])]['sw'].values

    #Calculate Doppler Broadening
    #linelist['GammaD'] = np.sqrt(2*1.380648813E-16*T*np.log(2)/linelist['m']/2.99792458e10**2)*linelist['nu']
    linelist['GammaD'] = np.sqrt(2*CONSTANTS['k']*CONSTANTS['Na']*T*np.log(2)/(linelist['m'].values))*linelist['nu'] / CONSTANTS['c']

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
                minimum_parameter_fit_intensity = 1e-27, weight_spectra = False,
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
                beta_formalism = False):


        self.dataset = dataset
        self.base_linelist_file = base_linelist_file
        self.baseline_list = pd.read_csv(self.base_linelist_file + '.csv')#, index_col = 0
        self.param_linelist_file = param_linelist_file
        self.lineparam_list = pd.read_csv(self.param_linelist_file + '.csv', index_col = 0)
        self.CIA_linelist_file = None
        self.CIAparam_list = None
        self.minimum_parameter_fit_intensity = minimum_parameter_fit_intensity
        self.weight_spectra = weight_spectra
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
                                  min = (1 / self.gamma0_limit_factor)*self.lineparam_list.loc[int(spec_line)][line_param],
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

    def simulation_model(self, params, wing_cutoff = 25, wing_wavenumbers = 25, wing_method = 'wing_cutoff'):
        """This is the model used for fitting that includes baseline, resonant absorption, and CIA models.


        Parameters
        ----------
        params : lmfit parameter object
            the params object is a dictionary comprised of all parameters translated from dataframes into a dictionary format compatible with lmfit.
        wing_cutoff : float, optional
            number of voigt half-widths to simulate on either side of each line. The default is 25.
        wing_wavenumbers : float, optional
            number of wavenumbers to simulate on either side of each line. The default is 25
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
            if ('molefraction' in param) or ('baseline' in param) or ('etalon' in param) or ('x_shift' in param) or ('Pressure' in param) or ('Temperature' in param) or ('_res_' in param):
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
                    if ('molefraction_'+ self.dataset.isotope_list[(molecule, 1)][4]) + '_' + str(spectrum_number) + '_' + str(segment) in baseline_params:
                        fit_molefraction[molecule] = np.float(params[('molefraction_'+ self.dataset.isotope_list[(molecule, 1)][4]) + '_' + str(spectrum_number) + '_' + str(segment)])
                #Get Environmental Parameters
                p = np.float(params['Pressure_' + str(spectrum_number) + '_' + str(segment)])
                T = np.float(params['Temperature_' + str(spectrum_number) + '_' + str(segment)])
                #Simulate Spectra

                if self.beta_formalism == True:
                    fit_nu, fit_coef = HTP_wBeta_from_DF_select(linelist_for_sim, wavenumbers, wing_cutoff = wing_cutoff, wing_wavenumbers = wing_wavenumbers, wing_method = wing_method,
                            p = p, T = T, molefraction = fit_molefraction, isotope_list = self.dataset.isotope_list,
                            natural_abundance = spectrum.natural_abundance, abundance_ratio_MI = spectrum.abundance_ratio_MI,  Diluent = Diluent)
                else:
                    fit_nu, fit_coef = HTP_from_DF_select(linelist_for_sim, wavenumbers, wing_cutoff = wing_cutoff, wing_wavenumbers = wing_wavenumbers, wing_method = wing_method,
                            p = p, T = T, molefraction = fit_molefraction, isotope_list = self.dataset.isotope_list,
                            natural_abundance = spectrum.natural_abundance, abundance_ratio_MI = spectrum.abundance_ratio_MI,  Diluent = Diluent)
                fit_coef = fit_coef * 1e6
                ## CIA Calculation
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
                segment_alpha = baseline + etalons + fit_coef + CIA
                #ILS_Function
                if spectrum.ILS_function != None:
                    if self.dataset.ILS_function_dict[spectrum.ILS_function.__name__] ==1:
                        spec_seg_ILS_resolution  = np.float(params[spectrum.ILS_function.__name__ + '_res_0_' + str(spectrum_number) + '_' +str(segment)])
                    else:
                        spec_seg_ILS_resolution  = []
                        for res_param in range(0, self.dataset.ILS_function_dict[spectrum.ILS_function.__name__]):
                            spec_seg_ILS_resolution  += np.float(params[spectrum.ILS_function.__name__ + '_res_' + str(res_param) +'_' + str(spectrum_number) + '_' +str(segment)])
                    wavenumbers, segment_alpha, i1, i2m, slit = convolveSpectrumSame(wavenumbers, segment_alpha, SlitFunction = spectrum.ILS_function, Resolution = spec_seg_ILS_resolution ,AF_wing=spectrum.ILS_wing)
                simulated_spectra[np.min(indices_segments[segment]): np.max(indices_segments[segment])+1] = (segment_alpha)
                #Weighted Spectra
                if self.weight_spectra:
                    if spectrum.tau_stats.all() == 0:
                        weights = len(alpha_segments[segment])*[spectrum.weight]
                    else:
                        pt_by_pt_weights= 1 / (spectrum.tau_stats[np.min(indices_segments[segment]): np.max(indices_segments[segment])+1])
                        weights = spectrum.weight * pt_by_pt_weights
                    residuals[np.min(indices_segments[segment]): np.max(indices_segments[segment])+1]  = ((segment_alpha) - alpha_segments[segment])*weights
                else:
                    residuals[np.min(indices_segments[segment]): np.max(indices_segments[segment])+1]  = (segment_alpha) - alpha_segments[segment]

            total_simulated = np.append(total_simulated, simulated_spectra)
            total_residuals = np.append(total_residuals, residuals)
        total_residuals = np.asarray(total_residuals)
        total_simulated = np.asarray(total_simulated)
        return total_residuals
    def fit_data(self, params, wing_cutoff = 25, wing_wavenumbers = 25, wing_method = 'wing_cutoff', xtol = 1e-7, maxfev = 2000, ftol = 1e-7):
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

        minner = Minimizer(self.simulation_model, params, xtol =xtol, max_nfev =  maxfev, ftol = ftol, fcn_args=(wing_cutoff, wing_wavenumbers, wing_method))
        result = minner.minimize(method = 'leastsq')#'
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
    def update_params(self, result, base_linelist_update_file = None , param_linelist_update_file = None):
        """Updates the baseline and line parameter files based on fit results with the option to write over the file (default) or save as a new file and updates baseline values in the spectrum objects.


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
            elif ('_res_' in par.name):
                indices = [m.start() for m in re.finditer('_', par.name[par.name.find('_res_') + 5:])]
                spec_num = par.name[par.name.find('_res_') + 5:][indices[0]+1:indices[1]]
                seg_num = par.name[par.name.find('_res_') + 5:][indices[1]+1:]
                parameter = par.name[:par.name.find('_res_')] + '_res_' + par.name[par.name.find('_res_') + 5:][:indices[0]]
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
                    beta_summary_list.loc[(beta_summary_list['molec_id']==molec) & (beta_summary_list['local_iso_id']==iso), 'm'] = molecularMass(molec,iso, isotope_list = self.dataset.isotope_list, )

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

                    beta_summary_list['alpha'] = mp / beta_summary_list['m']
                    if nu_constrain:
                        GammaD = np.sqrt(2*CONSTANTS['k']*CONSTANTS['Na']*T*np.log(2)/(beta_summary_list['m'].values))*beta_summary_list['nu'] / CONSTANTS['c']#change with nu
                    else:
                        GammaD = np.sqrt(2*CONSTANTS['k']*CONSTANTS['Na']*T*np.log(2)/(beta_summary_list['m'].values))*beta_summary_list['nu' + '_' + str(spectrum.spectrum_number)] / CONSTANTS['c'] #change with nu
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
