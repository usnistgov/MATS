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


from .spectrum import Spectrum, simulate_spectrum
from .dataset import Dataset
from .generate_fitparam_file import Generate_FitParam_File
from .fit_dataset import Fit_DataSet
from .utilities import add_to_HITRANstyle_isotope_list
from .linelistdata import linelistdata
from .o2_cia_karman import o2_cia_karman_model

try:
    import pkg_resources

    __version__ = pkg_resources.get_distribution("MATS").version
except Exception:
    # Local copy or not installed with setuptools.
    # Disable minimum version checks on downstream libraries.
    __version__ = "999"
