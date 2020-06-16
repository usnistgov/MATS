import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from hapi import PYTIPS2017, molecularMass, pcqsdhc, ISO
import qgrid
from bisect import bisect
import re
from lmfit import Parameters, Minimizer
from Karman_CIA import Karman_CIA_Model


#Constants
h = 6.62607015e-27 #erg s https://physics.nist.gov/cgi-bin/cuu/Value?h|search_for=h as of 5/21/2020
c = 29979245800 #cm/s # https://physics.nist.gov/cgi-bin/cuu/Value?c|search_for=c as of 5/21/2020
k = 1.380649e-16 # erg / K https://physics.nist.gov/cgi-bin/cuu/Value?k as of 5/21/2020     
Na = 6.02214076e23 # mol-1 https://physics.nist.gov/cgi-bin/cuu/Value?na as of 5/21/2020
cpa_atm = (10*101325)**-1 #convert from cpa to atm  https://physics.nist.gov/cgi-bin/cuu/Value?stdatm|search_for=atmosphere as of 5/21/2020
c2 =  (h*c)/k


           
def max_iter(pars, iter, resid, *args, **kws):
        if iter > 2500:
            return True
        else:
            return False
                   
def etalon(x, amp, period, phase):
    """Etalon definition
    

    Parameters
    ----------
    x : array
        array of floats used to define the x-axis
    amp : float
        amplitude of the etalon.
    period : float
        period of the etalon.
    phase : float
        phase of the etalon.

    Returns
    -------
    etalon : array
        etalon as a function of input x-axis, amplitude, period, and phase.

    """

    return amp*np.sin((2*np.pi * period)*x+ phase)   
 
def hasNumbers(inputString):
    """Determines whether there are numbers in a string    

    Parameters
    ----------
    inputString : str
        string for analysis

    Returns
    -------
    bool
        Returns True if the there are numbers in a string

    """

    for char in inputString:
        if char.isdigit():
            return True
    return False