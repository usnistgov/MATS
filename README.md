**The purpose of this project is to develop a NIST based multi-spectrum fitting tool that allows the flexibility to test and adapt to experimental/data-driven needs.  The tool uses HAPI and LMFit as the engine for spectral fitting.**  

The package is based on spectrum objects, which read in x, y column of the wavenumber and absorption from a file.  Additionally, the spectrum object contains information on the pressure, temperature (and nominal temperature), etalons, and sample composition, and baseline behavior.  Spectrum object also allows for segmentation of the spectrum through a segment column.  The spectrum objects are bundled together to form a dataset object, which is the collection of spectra that are being analyzed.  

There are two files that contain parameters that are fit in this model, one for spectrum dependent parameters (baseline parameters, etalons, sample composition, x-shift) and the other for linelist parameters that are common across all spectrum.  These files are saved locally as .csv files with a column for each parameter and with rows corresponding to either the spectrum number or spectral line number.  In the event that there are numerous segments in a spectrum file, the the background file contains rows corresponding to each spectrum and each spectrum within the segment.  In addition to the columns for fit parameter, there are two additional columns for each parameter called param_vary and param_err.  The param_vary column is a boolean flag that is toggled to indicate whether a given parameter will be varied in the fit.  The param_err column will be set to zero initially and replaced with the fit error value, if it is floated after fitting.  Calls of the generate fit parameter file class not only make these input files, but also set the lineshape (HTP derivatives) and define initial conditions for whether a parameter should be varied and whether the parameter should be constrained to multi-spectrum fit or allowed to vary by spectrum.  The edit fit parameter class allows for edits to the parameters files in the ipython interface rather than editing in the .csv file. 

The fitting class fits the data and allows you to impose constraints on the parameters (min and max values), which are defined as a factor of the absolute value (for most cases the x-shift and line center terms are absolute).  The fitting class also allows you to define the wing cut-off and fitting convergence criteria. 

In addition to fitting real spectra, there is a simulate_spectrum function, which returns a spectrum object that match input simulation parameters.  This is useful for performing error analysis in the same framework as primary data analysis.

**The process of fitting spectrum can be broken down into these basic components:**
1.  [Set Up and Load Packages](Set-Up-and-Load-Packages)
2.  Load Spectra by generating [Spectrum](Spectrum) class objects for each spectrum, which can done by loading spectrum .csv files or by using the [Simulate Spectrum](Simulate-Spectrum) function 
3.  Create a dataset by generating a [Dataset](Dataset) class object combining all loaded spectra together (or desired subset).
4.  Set-up fit parameters using [Generate Fit Parameters File](Generate-Fit-Parameters-File) class.  This will need a spectral parameter line list.  Provided along with the package is a separate jupyter notebook that generates one from HITRAN with the option to then manually update with the literature values of your choice.  [Generating Panda Linelists from HITRAN](Generating-Pandas-Linelists)
5.  Use the [Edit Fit Parameters File](Edit-Fit-Parameters-File) class to edit the fit parameters files in the jupyter notebook.  Additionally, fit parameter files can be edited in the generated .csv files.
6.  Use the [Fit Dataset](Fit-Dataset) files to perform fits.  Iterating through steps 5 and 6 will optimize fits.  


**Code documentation can be found here for each class/definition**
*  [HTP definition from Dataframe](HTP-definition-from-Dataframe)
*  [Spectrum](Spectrum)
*  [Simulate Spectrum](Simulate-Spectrum) 
*  [Dataset](Dataset)
*  [Generate Fit Parameters File](Generate-Fit-Parameters-File)
*  [Edit Fit Parameters File](Edit-Fit-Parameters-File)
*  [Fit Dataset](Fit-Dataset)

[**Examples**](https://gitlab.nist.gov/gitlab/ema3/HAPI-spectral-fitting/tree/master/From%20DataFrame/Examples)


[**Linelists**](https://gitlab.nist.gov/gitlab/ema3/HAPI-spectral-fitting/tree/master/From%20DataFrame/Linelists)


