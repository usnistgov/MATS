try: 
  from numba import jit        as numba_jit
  from numba import njit       as numba_njit
  from numba import float64    as numba_f8
  from numba import complex128 as numba_c16
except ImportError as msg: 
  raise SystemExit (str(msg) + '\nCPF.py:  Numba not found. Numba module is needed to run the code!')
try: 
  from numpy import ndarray    as numpy_ndarray
  from numpy import empty_like as numpy_empty_like
  from numpy import complex128 as numpy_complex128
except ImportError as msg: 
  raise SystemExit (str(msg) + '\nCPF.py:  Numpy not found. Numpy module is needed to run the code!')

@numba_jit(numba_c16(numba_f8,numba_f8), nopython=True, cache=True)
def cpf_accurate(x: float, y: float) -> complex:
  """ Accurate CPF algorithm
  =====
  Computes the complex probability function using a rational series with 42 terms. It is assumed that Im(z) > 0 or Im(z) = 0 (Source: jstor.org/stable/2158232). A series was simplified to 37 terms introducing less than 1e-17 deviations on mHT profile.

  Parameters
  ----------
  x : float
    Real part of input complex parameter
  y : float
    Imaginary part of input complex parameter

  Returns
  -------
  complex
    Complex probability function
  """
  z = -y + x*1j
  Z = (5.449631621480024+z)/(5.449631621480024-z)
  return (2*(+2.975931371735470E+00+Z*(+2.697763665856064E+00+Z*(+2.288734169675538E+00+Z*(+1.814714451499866E+00+Z*(+1.342044484596932E+00+Z*(+9.230959991941070E-01+\
          Z*(+5.882708203344523E-01+Z*(+3.455278077566057E-01+Z*(+1.857036333535562E-01+Z*(+9.038744880336540E-02+Z*(+3.922970169744468E-02+Z*(+1.480296368764821E-02+\
          Z*(+4.631075611097791E-03+Z*(+1.070131083157417E-03+Z*(+1.044072210002090E-04+Z*(-4.816473680511106E-05+Z*(-2.914574364851397E-05+Z*(-6.617492208403963E-06+\
          Z*(+4.936972061734341E-07+Z*(+8.511689670641750E-07+Z*(+2.206733163926054E-07+Z*(-2.752070035718561E-08+Z*(-3.294773119114329E-08+Z*(-5.468780625369738E-09+\
          Z*(+2.636162411919059E-09+Z*(+1.282896083944607E-09+Z*(-6.353807951660892E-11+Z*(-1.799169607159564E-10+Z*(-2.392389527320517E-11+Z*(+2.063579921011804E-11+\
          Z*(+6.430136110306704E-12+Z*(-2.069163661083667E-12+Z*(-1.184560208678836E-12+Z*(+1.790586243645278E-13+Z*(+1.951777029849348E-13+Z*(-1.188364999909099E-14-\
          Z*3.129493160727961E-14))))))))))))))))))))))))))))))))))))/(5.449631621480024-z)+0.5641895835477563)/(5.449631621480024-z)
      
@numba_jit(numba_c16[:](numba_f8[:], numba_f8[:]), nopython=True, cache=True)
def cpf_accurate_vector(x, y):
  """ Accurate CPF algorithm (vectorized version)
  =====
  Computes the complex probability function using a rational series with 42 terms. It is assumed that Im(z) > 0 or Im(z) = 0 (Source: jstor.org/stable/2158232). A series was simplified to 37 terms introducing less than 1e-17 deviations on mHT profile.

  Parameters
  ----------
  x : numpy.ndarray
    Real part of input complex parameter
  y : numpy.ndarray
    Imaginary part of input complex parameter

  Returns
  -------
  numpy.ndarray
    Complex probability function
  """
  z = -y + x*1j
  Z = (5.449631621480024+z)/(5.449631621480024-z)
  return (2*(+2.975931371735470E+00+Z*(+2.697763665856064E+00+Z*(+2.288734169675538E+00+Z*(+1.814714451499866E+00+Z*(+1.342044484596932E+00+Z*(+9.230959991941070E-01+\
          Z*(+5.882708203344523E-01+Z*(+3.455278077566057E-01+Z*(+1.857036333535562E-01+Z*(+9.038744880336540E-02+Z*(+3.922970169744468E-02+Z*(+1.480296368764821E-02+\
          Z*(+4.631075611097791E-03+Z*(+1.070131083157417E-03+Z*(+1.044072210002090E-04+Z*(-4.816473680511106E-05+Z*(-2.914574364851397E-05+Z*(-6.617492208403963E-06+\
          Z*(+4.936972061734341E-07+Z*(+8.511689670641750E-07+Z*(+2.206733163926054E-07+Z*(-2.752070035718561E-08+Z*(-3.294773119114329E-08+Z*(-5.468780625369738E-09+\
          Z*(+2.636162411919059E-09+Z*(+1.282896083944607E-09+Z*(-6.353807951660892E-11+Z*(-1.799169607159564E-10+Z*(-2.392389527320517E-11+Z*(+2.063579921011804E-11+\
          Z*(+6.430136110306704E-12+Z*(-2.069163661083667E-12+Z*(-1.184560208678836E-12+Z*(+1.790586243645278E-13+Z*(+1.951777029849348E-13+Z*(-1.188364999909099E-14-\
          Z*3.129493160727961E-14))))))))))))))))))))))))))))))))))))/(5.449631621480024-z)+0.5641895835477563)/(5.449631621480024-z)

@numba_jit(numba_c16(numba_f8,numba_f8), nopython=True, cache=True)
def cpf_fast(x: float, y: float) -> complex:
  """ Fast CPF algorithm
  =====
  Computes the complex probability function using Humlicek's algorithm in its first subregion (Source: 10.1016/0022-4073(82)90078-4) and using a rational series with 24 terms in other subregions (Source: jstor.org/stable/2158232).

  Parameters
  ----------
  x : float
    Real part of input complex parameter
  y : float
    Imaginary part of input complex parameter

  Returns
  -------
  complex
    Complex probability function
  """
  if abs(x)+y>15.0: # The border of the first Humlicek's region
    t=y-x*1j
    return t*0.5641895835477563/(0.5+t*t)
  else:
    z=-y+x*1j
    Z=(4.119534287814236+z)/(4.119534287814236-z)
    return (2*(+2.197858936531542E+00+Z*(+1.856286499205540E+00+Z*(+1.394819673379119E+00+Z*(+9.257087138588670E-01+Z*(+5.361139535729116E-01+Z*(+2.654963959880772E-01+\
            Z*(+1.083872348456673E-01+Z*(+3.372336685531603E-02+Z*(+6.215006362949147E-03+Z*(-4.936426901286291E-04+Z*(-7.816642995626165E-04+Z*(-2.074843151143828E-04+\
            Z*(+2.433141546207148E-05+Z*(+3.047106608295325E-05+Z*(+4.139461724429617E-06+Z*(-3.038893184366094E-06+Z*(-1.085647579417637E-06+Z*(+2.568264135399530E-07+\
            Z*(+1.873834346505099E-07+Z*(-1.912225887484805E-08+Z*(-3.008282344381996E-08+Z*(+1.331045329581992E-09+Z*(+4.904820407381768E-09-\
	    Z*1.513747622620502E-10)))))))))))))))))))))))/(4.119534287814236-z)+0.5641895835477563)/(4.119534287814236-z)        

@numba_njit(cache=True)
def cpf_fast_vector(x: numpy_ndarray, y: numpy_ndarray) -> numpy_ndarray:
  """ Fast CPF algorithm (vectorized version)
  =====
  Computes the complex probability function using Humlicek's algorithm in its first subregion (Source: 10.1016/0022-4073(82)90078-4) and using a rational series with 24 terms in other subregions (Source: jstor.org/stable/2158232).
  
  Parameters
  ----------
  x : numpy.ndarray
    Real part of input complex parameter
  y : numpy.ndarray
    Imaginary part of input complex parameter
  
  Returns
  -------
  numpy.ndarray
    Complex probability function
  """
  res = numpy_empty_like(x, dtype=numpy_complex128)
  for i in range(x.size):
    xi = x[i]
    yi = y[i]
    if abs(xi) + yi > 15.0:
      t      = yi - 1j*xi
      res[i] = t*0.5641895835477563/(0.5+t*t)
    else:
      z      = -yi + 1j*xi
      Z      = (4.119534287814236 + z) / (4.119534287814236 - z)
      res[i] = (2*(+2.197858936531542E+00+Z*(+1.856286499205540E+00+Z*(+1.394819673379119E+00+Z*(+9.257087138588670E-01+Z*(+5.361139535729116E-01+Z*(+2.654963959880772E-01+\
                Z*(+1.083872348456673E-01+Z*(+3.372336685531603E-02+Z*(+6.215006362949147E-03+Z*(-4.936426901286291E-04+Z*(-7.816642995626165E-04+Z*(-2.074843151143828E-04+\
                Z*(+2.433141546207148E-05+Z*(+3.047106608295325E-05+Z*(+4.139461724429617E-06+Z*(-3.038893184366094E-06+Z*(-1.085647579417637E-06+Z*(+2.568264135399530E-07+\
                Z*(+1.873834346505099E-07+Z*(-1.912225887484805E-08+Z*(-3.008282344381996E-08+Z*(+1.331045329581992E-09+Z*(+4.904820407381768E-09-\
      	        Z*1.513747622620502E-10)))))))))))))))))))))))/(4.119534287814236-z)+0.5641895835477563)/(4.119534287814236-z)  
  return res