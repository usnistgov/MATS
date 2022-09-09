#Import Packages
import numpy as np
import pandas as pd
import scipy.interpolate as interpolate


def o2_cia_karman_model(wavenumbers, T, P,molefraction,
                        SO_O2, SO_N2, EXCH_O2, 
                        EXCH_b, EXCH_c, 
                        SO_b, SO_c, 
                        SO_shift = 0, EXCH_shift = 0,
                        band = 'singlet_delta'):
    
    if band == 'singlet_delta':
        df_Karman = pd.read_csv('Singlet_Delta_Karman.csv')
    elif band == 'a_band':
        df_Karman = pd.read_csv('A_Band_Karman.csv')
    model_wavenumber = df_Karman['Wavenumber (cm-1)'].values
    SO = (df_Karman['SO Temp Dep - d'].values*(T-296)**3 +df_Karman['SO Temp Dep - c'].values*(T-296)**2 + df_Karman['SO Temp Dep - b'].values*(T-296) + df_Karman['SO Temp Dep - a'].values)*df_Karman['Spin Orbit (cm-1 amagat-2) Normalized'].values
    EXCH = (df_Karman['EXCH Temp Dep - d'].values*(T-296)**3 +df_Karman['EXCH Temp Dep - c'].values*(T-296)**2 + df_Karman['EXCH Temp Dep - b'].values*(T-296) + df_Karman['EXCH Temp Dep - a'].values)*df_Karman['Exchange (cm-1 amagat-2) Normalized'].values
    
    f = interpolate.interp1d(wavenumbers, SO, fill_value = 'extrapolate', bounds_error = False, kind = 'slinear')
    SO_ = f(x_ - SO_shift)
    g = interpolate.interp1d(wavenumbers, EXCH, fill_value = 'extrapolate', bounds_error = False, kind = 'slinear')
    EXCH_ = g(x_ - EXCH_shift)
    
    model_O2_O2 = EXCH_O2*(1 + EXCH_b*(T-296) + EXCH_c*(T-296)**2)*EXCH_ + SO_O2*(1 + SO_b*(T-296) + SO_c*(T-296)**2)*SO_
    model_O2_N2 = SO_N2*(1 + SO_b*(T-296) + SO_c*(T-296)**2)*SO_

    CIA_model = molefraction['O2']*model_O2_O2 + molefraction['N2']*model_O2_N2
    
    return CIA_model
  
    
