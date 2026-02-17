from math import tanh as math_tanh
from math import log10 as math_log10
from cmath import sqrt as cmath_sqrt
try: 
  from numpy import array as numpy_array
  from numpy import median as numpy_median
except ImportError as msg: 
  raise SystemExit (str(msg) + '\nprofile.py:  Numpy not found. Numpy module is needed to run the code!')

# -----------------------------------  
# Choice of CPF, comment one of below
# -----------------------------------
#from mHT.CPF import cpf_accurate_vector as cpf_vector, cpf_accurate as cpf
# from mHT.CPF import cpf_fast_vector as cpf_vector, cpf_fast as cpf
from .CPF import cpf_accurate_vector as cpf_vector, cpf_accurate as cpf

def mHTprofile(nu0: float, GammaD: float, Gamma0: float, Gamma2: float, Delta0: float, Delta2: float, NuOptRe: float, NuOptIm: float, nu: float, *args) -> float:
  """ Modified Hartman Tran profile
  =====
  Subroutine to compute the complex normalized spectral-line shape using modified Hartman Tran profile model.
  
  Parameters
  ----------
  nu0 : float
    Unperturbed line position in cm-1.
  GammaD : float
    Doppler broadening in cm-1.
  Gamma0 : float
    Speed-averaged line-width in cm-1.  
  Gamma2 : float
    Unperturbed line position in cm-1.
  Delta0 : float
    Doppler broadening in cm-1.
  Delta2 : float
    Speed-averaged line-width in cm-1.  
  NuOptRe : float
    Unperturbed line position in cm-1.
  NuOptIm : float
    Doppler broadening in cm-1.
  nu : float
    Speed-averaged line-width in cm-1. 
  Ylm : float [optional, default=0.0]
    Imaginary part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless.
  Xlm : float [optional, default=0.0]
    Real part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless.
  alpha : float [optional, default=10.0]
    Mass ratio in the molecule for calculating beta-correction, applicable up to alpha=5, dimensionless.
  disp : boolean [optional, default=False]
    Boolean trigger for including dispersion profile in the output. 
  
  Returns
  -------
  float
    Real or imaginary (depending on disp value) part of the normalized spectral shape in cm.
  """
  nuD = 1.2011224087864498*GammaD

  match len(args):
    case 4:
      disp  = args[3]
      nuR   = NuOptRe*beta(GammaD,NuOptRe,args[2])
      LM    = args[1] + 1.0 - args[0]*1j
    case 3:
      disp  = False
      nuR   = NuOptRe*beta(GammaD,NuOptRe,args[2])
      LM    = args[1] + 1.0 - args[0]*1j
    case 2:
      disp  = False
      nuR   = NuOptRe
      LM    = args[1] + 1.0 - args[0]*1j
    case 1:
      disp  = False
      nuR   = NuOptRe
      LM    = 1.0 - args[0]*1j
    case _:
      disp  = False
      nuR   = NuOptRe
      LM    = 1.0   

  c2  = Gamma2+Delta2*1j
  c0  = Gamma0+Delta0*1j-1.5*c2+nuR+NuOptIm*1j  
  if abs(c2) > 1.0e-9: # 1.0e-9 - limit where speed dependence impact is lower than numerical noise level
    X    = ((nu0-nu)*1j+c0)/c2
    Y    = 0.25*(nuD/c2)**2.0
    if abs(Y)>abs(X)*1.0e-15: # 1.0e-15 - numerical zero  
      csqY = cmath_sqrt(Y)
      z2   = (X+Y)**0.5+csqY   
      z1   = z2-2*csqY if abs(X)>abs(Y)*3e-8 else ((nu0-nu)*1j+c0)/nuD    
      w1   = cpf(-z1.imag,z1.real)
      w2   = cpf(-z2.imag,z2.real)
      A    = 1.772453850905516/nuD*(w1-w2)
    else:
      rX = X**0.5
      if abs(rX) < 4.0e3: # 4.0e3 - numerical infinity             
        wX = cpf(-rX.imag,rX.real)
        A  = (2-3.5449077018110318*rX*wX)/c2
      else: A = (1-1.5/X)/X/c2
  else:
    z = ((nu0-nu)*1j+c0)/nuD
    w = cpf(-z.imag,z.real)
    A = 1.772453850905516*w/nuD
  I = 0.3183098861837907*LM/(1/A-(nuR+NuOptIm*1j))
    
  match disp:
    case True:
      return -I.imag
    case _: # If disp variable was defined as not boolean its equivalent to not(True)
      return I.real          

def mHTprofile_vector(nu0: float, GammaD: float, Gamma0: float, Gamma2: float, Delta0: float, Delta2: float, NuOptRe: float, NuOptIm: float, nu, *args) -> float:
  """ Modified Hartman Tran profile (vectorized version)
  =====
  Subroutine to compute the complex normalized spectral-line shape using modified Hartman Tran profile model.
  
  Parameters
  ----------
  nu0 : float
    Unperturbed line position in cm-1.
  GammaD : float
    Doppler broadening in cm-1.
  Gamma0 : float
    Speed-averaged line-width in cm-1.  
  Gamma2 : float
    Unperturbed line position in cm-1.
  Delta0 : float
    Doppler broadening in cm-1.
  Delta2 : float
    Speed-averaged line-width in cm-1.  
  NuOptRe : float
    Unperturbed line position in cm-1.
  NuOptIm : float
    Doppler broadening in cm-1.
  nu : float
    Speed-averaged line-width in cm-1. 
  Ylm : float [optional, default=0.0]
    Imaginary part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless.
  Xlm : float [optional, default=0.0]
    Real part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless.
  alpha : float [optional, default=10.0]
    Mass ratio in the molecule for calculating beta-correction, applicable up to alpha=5, dimensionless.
  disp : boolean [optional, default=False]
    Boolean trigger for including dispersion profile in the output. 
  
  Returns
  -------
  float
    Real or imaginary (depending on disp value) part of the normalized spectral shape in cm.
  """
  nuD = 1.2011224087864498*GammaD
  #nu  = numpy_array(nu,dtype=float).flatten()
  
  match len(args):
    case 4:
      disp  = args[3]
      nuR   = NuOptRe*beta(GammaD,NuOptRe,args[2])
      LM    = args[1] + 1.0 - args[0]*1j
    case 3:
      disp  = False
      nuR   = NuOptRe*beta(GammaD,NuOptRe,args[2])
      LM    = args[1] + 1.0 - args[0]*1j
    case 2:
      disp  = False
      nuR   = NuOptRe
      LM    = args[1] + 1.0 - args[0]*1j
    case 1:
      disp  = False
      nuR   = NuOptRe
      LM    = 1.0 - args[0]*1j
    case _:
      disp  = False
      nuR   = NuOptRe
      LM    = 1.0
  
  c2  = Gamma2+Delta2*1j
  c0  = Gamma0+Delta0*1j-1.5*c2+nuR+NuOptIm*1j
  
  if abs(c2) > 1.0e-9: # 1.0e-9 - limit where speed dependence impact is lower than numerical noise level
    X    = ((nu0-nu)*1j+c0)/c2
    Y    = 0.25*(nuD/c2)**2.0
    if c0.real>1e-9: vecabsX = min(abs(X));
    else:            vecabsX = numpy_median(abs(X));
    if abs(Y)>vecabsX*1.0e-15: # 1.0e-15 - numerical zero  
      csqY = cmath_sqrt(Y)
      z2   = (X+Y)**0.5+csqY   
      z1   = z2-2*csqY if vecabsX>abs(Y)*3e-8 else ((nu0-nu)*1j+c0)/nuD    
      w1   = cpf_vector(-z1.imag,z1.real)
      w2   = cpf_vector(-z2.imag,z2.real)
      A    = 1.772453850905516/nuD*(w1-w2)
    else:
      rX = X**0.5
      if c0.real>1e-9: vecabssqX = min(abs(rX))
      else:            vecabssqX = numpy_median(abs(rX))
      if vecabssqX < 4.0e3: # 4.0e3 - numerical infinity             
        wX = cpf_vector(-rX.imag,rX.real)
        A  = (2-3.5449077018110318*rX*wX)/c2
      else: A = (1-1.5/X)/X/c2
  else:
    z = ((nu0-nu)*1j+c0)/nuD
    w = cpf_vector(-z.imag,z.real)
    A = 1.772453850905516*w/nuD
  I = 0.3183098861837907*LM/(1/A-(nuR+NuOptIm*1j))
  
  match disp:
    case True:
      return -I.imag
    case _: # If disp variable was defined as not boolean its equivalent to not(True)
      return I.real 

def beta(GammaD: float, NuOptRe: float, alpha: float) -> float:
  """ Beta-Correction Function
  =====
  Subroutine to compute beta-correction used for hard-collision based line-shape profiles. To correct NuOptRe value in the profile . Applicable up to alpha = 5.0, for higher alpha values correction neglected. Source: 10.1016/j.jqsrt.2019.106784.

  Parameters
  ----------
  GammaD : float
    Doppler broadening in cm-1.
  NuOptRe : float
    Real part of the Dicke parameter in cm-1.
  alpha : float
    Mass ratio in the molecule, applicable up to alpha = 5.0, dimensionless.

  Returns
  -------
  float
    Value of the beta correction, dimensionless.
  """
  if alpha<5.0: # the mass ratio up to which the beta correction is applicable
    # <a> *math_tanh( <b> *math_log10(NuOptRe/GammaD)+ <c> )+ <d>
    return (0.0534+0.1585*2.718281828459045**(-0.4510*alpha))\
           *math_tanh(\
           (1.9595+alpha*(-0.1258+alpha*(0.0056+alpha*0.0050)))\
           *math_log10(NuOptRe/GammaD)+\
           (-0.0546+alpha*(0.0672+alpha*(-0.0125+alpha*0.0003)))\
           )+\
           0.9466-0.1585*2.718281828459045**(-0.4510*alpha)
  else:
    return 1.0