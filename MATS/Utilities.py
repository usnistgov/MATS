import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from hapi import PYTIPS2017, pcqsdhc, ISO, ISO_INDEX
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

def molecularMass(M,I, isotope_list = ISO):
    """ molecular mass look-up based on the HAPI definition adapted to allow used to specify ISO list
    INPUT PARAMETERS: 
        M: HITRAN molecule number
        I: HITRAN isotopologue number
    OUTPUT PARAMETERS: 
        MolMass: molecular mass
    ---
    DESCRIPTION:
        Return molecular mass of HITRAN isotolopogue.
    ---
    EXAMPLE OF USAGE:
        mass = molecularMass(1,1) # H2O
    ---
    """
    return isotope_list[(M,I)][ISO_INDEX['mass']]

def isotope_list_molecules_isotopes(isotope_list = ISO):
    ''' The HITRAN style isotope list in the format (M,I), this function creates a dictionary from this with M as the keys and lists of I as values.
    

    Parameters
    ----------
    isotope_list : dictionary, optional
        HITRAN style isotope list. The default is ISO.

    Returns
    -------
    molecule_isotope_dictionary : dictionary
        Dictionary with the molecules in an isotope list with the available I values.

    '''
    M_I_in_isotope_list = (isotope_list.keys())
    molecule_isotope_dictionary = {}
    for MI in M_I_in_isotope_list:
        if MI[0] not in molecule_isotope_dictionary.keys():
            molecule_isotope_dictionary[MI[0]] = [MI[1]]
        else:
            isotope_list = molecule_isotope_dictionary[MI[0]] 
            isotope_list.append(MI[1])
            molecule_isotope_dictionary[MI[0]] = isotope_list    
    return molecule_isotope_dictionary

def add_to_HITRANstyle_isotope_list(input_isotope_list = ISO, molec_id = 100, local_iso_id = 1, global_isotope_id = 200, iso_name = '', 
                                    abundance = 1, mass = 1, mol_name = ''):
    '''  Allows for used to add to an existing isotope line list in the HITRAN format
    

    Parameters
    ----------
    input_isotope_list : dictionary, optional
        Isotope list dictionary in HAPI format. The default is ISO.
    molec_id : int, optional
        molec_id for new isotope entry. The default is 100.
    local_iso_id : int, optional
        local isotope id number for new isotope entry. The default is 1.
    global_isotope_id : int, optional
        global isotope id for new istope entry. The default is 200.
    iso_name : str, optional
        isotoope name for new isotope entry. The default is ''.
    abundance : float, optional
        relative abundance of new isotope entry. The default is 1.
    mass : float, optional
        Mass of isotope in g. The default is 1.
    mol_name : str, optional
        name of the molecule for new istope entry. The default is ''.

    Returns
    -------
    output_isotope_list : dictionary
        Isotope list with additionall entry

    '''
    # check that there is not another entry
    
    if (molec_id, local_iso_id) in input_isotope_list:
        print ('Already entry with that molec_id and local_iso_id.  This result will write over that entry.')
    elif (molec_id, 1) in  input_isotope_list:
        mol_name_ = input_isotope_list[(molec_id, 1)][ISO_INDEX['mol_name']]
        print ('This is being added as a new isotope of ' + mol_name_ + '.  Proceed if that was the intention.')
    global_iso_id_list = []
    for entry in input_isotope_list:
        global_iso_id_list.append(input_isotope_list[entry][ISO_INDEX['id']])
    if global_isotope_id in global_iso_id_list:
        print ('There is another entry with this global isotope id.  Consider changing the value for consistency')
    output_isotope_list = input_isotope_list.copy()
    output_isotope_list[(molec_id, local_iso_id)] = [global_isotope_id, iso_name, abundance, mass, mol_name]
    
    
    
    return output_isotope_list
        