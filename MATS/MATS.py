"""
MATS: Multi-spectrum Analysis Tool for Spectroscopy
===================================================

The purpose of the MATS project is to develop a NIST-based multi-spectrum fitting 
and analysis tool for spectroscopic data that allows the flexibility to test and 
adapt to experimental and data-driven needs. 

This software allows for the use of several commonly used spectroscopic line profiles 
(Voigt, Nelkin-Ghatak, speed-dependent Voigt, speed-dependent Nelkin-Ghatak, and Hartmann-Tran) 
and allows for pressure, temperature, and sample composition constraints to be imposed during fits. 

In addition to fitting experimental spectra, MATS can generate simulated spectra, 
which allows for its use as an error analysis tool.
"""




from Spectrum import *
from Dataset import *
from Generate_FitParam_File import *
from Fit_Dataset import *




       