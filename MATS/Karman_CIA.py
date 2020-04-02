import numpy as np
from scipy.special import factorial
import mpmath
from scipy import interpolate



def clebsch(j1, j2, j3, m1, m2, m3):
    #Provided by Tijs Karman
    """Calculates the Clebsch-Gordon coefficient
    for coupling (j1,m1) and (j2,m2) to give (j3,m3).

    Parameters
    ----------
    j1 : float
        Total angular momentum 1.

    j2 : float
        Total angular momentum 2.

    j3 : float
        Total angular momentum 3.

    m1 : float
        z-component of angular momentum 1.

    m2 : float
        z-component of angular momentum 2.

    m3 : float
        z-component of angular momentum 3.

    Returns
    -------
    cg_coeff : float
        Requested Clebsch-Gordan coefficient.

    """
    if m3 != m1 + m2:
        return 0
    vmin = int(np.max([-j1 + j2 + m3, -j1 + m1, 0]))
    vmax = int(np.min([j2 + j3 + m1, j3 - j1 + j2, j3 + m3]))

    C = np.sqrt((2.0 * j3 + 1.0) * factorial(j3 + j1 - j2) *
                factorial(j3 - j1 + j2) * factorial(j1 + j2 - j3) *
                factorial(j3 + m3) * factorial(j3 - m3) /
                (factorial(j1 + j2 + j3 + 1) *
                factorial(j1 - m1) * factorial(j1 + m1) *
                factorial(j2 - m2) * factorial(j2 + m2)))
    S = 0
    for v in range(vmin, vmax + 1):
        S += (-1.0) ** (v + j2 + m2) / factorial(v) * \
            factorial(j2 + j3 + m1 - v) * factorial(j1 - m1 + v) / \
            factorial(j3 - j1 + j2 - v) / factorial(j3 + m3 - v) / \
            factorial(v + j1 - j2 - m3)
    C = C * S
    return C



def SO(omegas, ahard=7 , lso=2, Nmax=31, Temp=296):
    #Provided by Tijs Karman
    # Constants
    cm1=4.5563353e-6
    mO=2.9156946e+4
    kboltz=3.1668152e-6
    # O2 constants
    mumass  = mO;
    BROTX   = 1.4376766*cm1;
    #partition function
    Zrot=0.0;
    for n in range(1,Nmax+2,2) :
        Zrot=Zrot+(2.0*n+1.0)*np.exp(-BROTX*n*(n+1.0)/(kboltz*Temp))
        #print(Zrot)

    VGSO=np.zeros(len(omegas))
    k=np.sqrt(2.0*mumass*kboltz*Temp)
    for n1 in range(1,Nmax+2,2) :
        E1=BROTX*n1*(n1+1.0) 
        P1=np.exp(-E1/(kboltz*Temp))/Zrot
        for n2 in range(abs(n1-lso),n1+lso+2,2) :
            fc=clebsch(n1,lso,n2,0,0,0)
            E2=BROTX*n2*(n2+1.0)
            dE=E2-E1
            fc=fc*fc*(2*n1+1)*P1
            for iw in range(len(omegas)):
                kp=np.sqrt(2.0*mumass*(kboltz*Temp+abs(omegas[iw]*cm1-dE)))
                mga=(ahard*(k-kp)/2)**2
                if abs(mga) < 260 :
                    mgf=mpmath.meijerg([[0],[]], [[0,1.5,2],[]], abs(mga)+1e-300)/(k*kp)
                else :
                    mgf=0.0
                mgf=mgf*(1.0+np.exp((omegas[iw]*cm1-dE)/(kboltz*Temp)))/(2.0*np.cosh((omegas[iw]*cm1-dE)/(kboltz*Temp)))
                VGSO[iw]=VGSO[iw] + fc*mgf
    return VGSO/max(VGSO)


def EXCH(omegas, gamma=3, lexch=2, Nmax=31, Temp=296):
    #Provided by Tijs Karman
    # Constants
    cm1=4.5563353e-6
    mO=2.9156946e+4
    kboltz=3.1668152e-6
    # O2 constants
    mumass  = mO;
    BROTX   = 1.4376766*cm1;
    #partition function
    Zrot=0.0;
    for n in range(1,Nmax+2,2) :
        Zrot=Zrot+(2.0*n+1.0)*np.exp(-BROTX*n*(n+1.0)/(kboltz*Temp))
        #print(Zrot)

    VGEXCH=np.zeros(len(omegas))
    k=np.sqrt(2.0*mumass*kboltz*Temp)
    for n1 in range(1,Nmax+2,2) :
        E1=BROTX*n1*(n1+1.0) 
        P1=np.exp(-E1/(kboltz*Temp))/Zrot
        for n2 in range(abs(n1-lexch),n1+lexch+2,2) :
            fc=clebsch(n1,lexch,n2,0,0,0)
            E2=BROTX*n2*(n2+1.0)
            dE=E2-E1
            fc=fc*fc*(2*n1+1)*P1
            for iw in range(len(omegas)):
                kp=np.sqrt(2.0*mumass*(kboltz*Temp+abs(omegas[iw]*cm1-dE)))
                mgf=gamma**2*k*kp/( ( (gamma**2+k**2)**2+2*kp**2*(gamma**2-k**2) + kp**4 )**2 )
                mgf=mgf*(1.0+np.exp((omegas[iw]*cm1-dE)/(kboltz*Temp)))/(2.0*np.cosh((omegas[iw]*cm1-dE)/(kboltz*Temp)))
                VGEXCH[iw]=VGEXCH[iw] + fc*mgf
    return VGEXCH/max(VGEXCH)



def Karman_CIA_Model(wavenumbers, pressure, temperature, wave_step = 5,
                     EXCH_scalar = 1, EXCH_gamma = 3, EXCH_l = 2,
                     SO_scalar = 1, SO_ahard = 7, SO_l = 2,
                     band_center = 13122, Nmax = 31):
    #X-Axis
    wave_min = np.min(wavenumbers)
    wave_max = np.max(wavenumbers)
    sim_wave = np.arange(wave_min, wave_max + wave_step, wave_step)
    omegas = sim_wave - band_center
    #Density Correction
    amagat = (760 / (62.36357736*273.15)) #mol / L
    at_P_T = (pressure/ (62.36357736*(temperature+273.15)))
    density_correction = ((at_P_T / amagat)**2)
    #Calculate Exchange and Spin Orbit Components
    if EXCH_scalar == 0:
        VGEXCH = len(omegas)*[0]
    else:
        VGEXCH = EXCH_scalar*EXCH(omegas, gamma=EXCH_gamma, lexch=EXCH_l, Nmax=Nmax, Temp=temperature + 273.15)
    if SO_scalar == 0:
        VGSO = len(omegas)*[0]
    else:
        VGSO = SO_scalar*SO(omegas, ahard=SO_ahard , lso=SO_l, Nmax=Nmax, Temp=temperature + 273.15)
    #Add the Exchange and SO components and correct for density
    CIA_model_sim = (VGEXCH + VGSO)*density_correction
    f = interpolate.interp1d(sim_wave, CIA_model_sim)
    return f(wavenumbers)
    

