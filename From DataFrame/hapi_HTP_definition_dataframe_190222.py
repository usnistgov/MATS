import numpy as np
import pandas as pd
from bisect import bisect
import os, sys
import matplotlib.pyplot as plt
sys.path.append(r'C:\Users\ema3\Documents\Python Scripts\HAPI')#Add hapi.py folder location to system path
from hapi import *

# proposed inputs

'''
linelist = pd.read_csv('TestTable_initguess.csv')
wavenumbers = np.arange(13040, 13045, 0.001)
wing_cutoff = 50 # wingcut-off will be x halfwidths or end of range which ever is later
wing_wavenumbers = 50
pressure = 1#atm
temperature = 296 #K
concentration ={ 7 : 0.2} # Concentration is a dictionary with a concentration (out of 1) with the Molecule ID as the key
natural_abundance = True # if false then need to correct calculations
abundance_ratio_MI = {7 : {1: 1 /0.995262 , 2: 0, 3: 0}}
diluent = 'air' # can define also as dictionary that it loops through with the value being the abun for the calc of each param
Diluent = {}
'''



def HTP_from_DF_select(linelist, wavenumbers, wing_cutoff = 50, wing_wavenumbers = 50, 
                pressure = 1, temperature = 296, concentration = [], 
                natural_abundance = True, abundance_ratio_MI = {},  Diluent = {}, diluent = 'air', IntensityThreshold = 1e-30):


    #Set Omegas to X-values
    Xsect = [0]*len(wavenumbers)
    
    #define reference temperature and pressure
    Tref = 296. # K
    pref = 1. # atm
    # define actual temperature and pressure
    T = temperature # K
    p = pressure # atm
   
    mol_dens = volumeConcentration(p,T)
    
    
    #Sets-up the  Diluent
    if not Diluent:
        if diluent == 'air':
            Diluent = {'air':1.}
        elif diluent == 'self':
            Diluent = {'self':1.}
        else:
            raise Exception('Unknown GammaL value: %s' % GammaL)
              
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
        if (not natural_abundance) and abundance_ratio_MI != {}:
            abun_ratio = abundance_ratio_MI[MoleculeNumberDB][IsoNumberDB]
            print (abun_ratio)

        
        
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
        
        
        
        
        Xsect[BoundIndexLower:BoundIndexUpper] += mol_dens / abundance(MoleculeNumberDB,IsoNumberDB) * \
                                                    concentration[MoleculeNumberDB] * abun_ratio * \
                                                    LineIntensity * (lineshape_vals_real + Y*lineshape_vals_imag) 
    return (wavenumbers, Xsect)        

'''
nu, alpha = HTP_from_DF_select(linelist, wavenumbers, wing_cutoff = 50, wing_wavenumbers = 50, 
                pressure = pressure, temperature = temperature, concentration =concentration, 
                natural_abundance = False, abundance_ratio_MI = abundance_ratio_MI, Diluent = {}, diluent = 'air') 

                                                              
plt.plot(nu, alpha)
plt.show()
'''

    
    

    


                
