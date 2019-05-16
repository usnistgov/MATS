**The purpose of this project is to develop a NIST based multi-spectrum fitting tool that allows the flexibility to test and adapt to experimental/data-driven needs.  The tool uses HAPI and LMFit as the engine for spectral fitting.**  

The package is based on spectrum objects, which read in x, y column of the wavenumber and absorption from a file.  Additionally, the spectrum object contains information on the pressure, temperature (and nominal temperature), etalons, and sample composition.  The spectrum objects are bundled together to form a dataset object, which is the collection of spectra that are being analyzed.  Additionally, the dataset defines the polynomial order of the baseline for all of spectra in the dataset.  

There are two files that contain parameters that are fit in this model, one for spectrum dependendent parameters (baseline parameters, etalons, sample composition, x-shift) and the other for linelist parameters that are common across spectrum.  These files are saved locally as .csv files with a column for each parameter and with rows corresponding to either the spectrum number or spectral line number.  In addition to the columns for fit parameter, there are two additional columns for each parameter that called param_vary and param_err.  The param_vary column is a boolean flag that is toggled to indicate whether a given parameter will be varied in the fit.  The param_err column will be set to zero initially and replaced with the fit error value, if it is floated after fitting.  Calls of the generate fit parameter file class not only make these input files, but also set the lineshape (HTP derivatives) and define initial conditions for 
whether a parameter should be varied (line intensity minimum and parameter), and whether the parameter should be constrained to multi-spectrum fit or allowed to vary by spectrum.  The edit fit parameter class allows for edits to the parameters files in the ipython interface rather than editing in the .csv file. 

The fitting class imposes fits the data and allows you to impose constraints on the parameters (min and max values), which are defined as a factor of the absolute value (for most cases the x-shift and line center terms are absolute).  The fitting class also allows you to define the wing cut-off and spacing for simulation.  After fitting is done the residuals parameter for each spectrum will be updated to reflect the current model and the parameter files will be updated.  

**The process of fitting spectrum can be broken down into these basic components:**
1.  [Set Up and Load Packages](https://gitlab.nist.gov/gitlab/ema3/HAPI-spectral-fitting/wikis/Set-Up-and-Load-Packages)
2.  Load Spectra by generating [Spectrum](https://gitlab.nist.gov/gitlab/ema3/HAPI-spectral-fitting/wikis/Spectrum) class objects for each spectrum
3.  Create a dataset by generating a [Dataset](https://gitlab.nist.gov/gitlab/ema3/HAPI-spectral-fitting/wikis/Dataset) class object combining all of the loaded spectra together (or desired subset.
4.  Set-up fit parameters using [Generate Fit Parameters File](https://gitlab.nist.gov/gitlab/ema3/HAPI-spectral-fitting/wikis/Generate-Fit-Parameters-File) class.  This will need a spectral parameter line list.  Provided along with the package is a separate jupyter notebook that generates one from HITRAN with the option to then manually update with the literature values of your choice.  [Generating Panda Linelists from HITRAN](https://gitlab.nist.gov/gitlab/ema3/HAPI-spectral-fitting/wikis/Generating-Pandas-Linelists)
5.  Use the [Edit Fit Parameters File](https://gitlab.nist.gov/gitlab/ema3/HAPI-spectral-fitting/wikis/Edit-Fit-Parameters-File) class to edit the fit parameters files in the jupyter notebook.  Additionally, fit parameter files can be edited in the generated .csv files.
6.  Use the [Fit Dataset](https://gitlab.nist.gov/gitlab/ema3/HAPI-spectral-fitting/wikis/Fit-Dataset) files to perform fits.  Iterating through steps 5 and 6 will optimize fits.  


**Code documentation can be found here for each class/definition**
*  [HTP definition from Dataframe](https://gitlab.nist.gov/gitlab/ema3/HAPI-spectral-fitting/wikis/HTP-definition-from-Dataframe)
*  [Spectrum](https://gitlab.nist.gov/gitlab/ema3/HAPI-spectral-fitting/wikis/Spectrum)
*  [Dataset](https://gitlab.nist.gov/gitlab/ema3/HAPI-spectral-fitting/wikis/Dataset)
*  [Generate Fit Parameters File](https://gitlab.nist.gov/gitlab/ema3/HAPI-spectral-fitting/wikis/Generate-Fit-Parameters-File)
*  [Edit Fit Parameters File](https://gitlab.nist.gov/gitlab/ema3/HAPI-spectral-fitting/wikis/Edit-Fit-Parameters-File)
*  [Fit Dataset](https://gitlab.nist.gov/gitlab/ema3/HAPI-spectral-fitting/wikis/Fit-Dataset)

[**Examples**](https://gitlab.nist.gov/gitlab/ema3/HAPI-spectral-fitting/tree/master/From%20DataFrame/Examples)


[**Linelists**](https://gitlab.nist.gov/gitlab/ema3/HAPI-spectral-fitting/tree/master/From%20DataFrame/Linelists)
