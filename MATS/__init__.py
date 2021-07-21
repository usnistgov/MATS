from .spectrum import Spectrum, simulate_spectrum
from .dataset import Dataset
from .generate_fitparam_file import Generate_FitParam_File
from .fit_dataset import Fit_DataSet
from .utilities import add_to_HITRANstyle_isotope_list

from .linelistdata import linelistdata

try:
    import pkg_resources

    __version__ = pkg_resources.get_distribution("MATS").version
except Exception:
    # Local copy or not installed with setuptools.
    # Disable minimum version checks on downstream libraries.
    __version__ = "999"
