MATS Version Summary
====================

MATS Version 2
++++++++++++++
published on 6/10/2021

NEW FEATURES
* Beta correction to the Dicke Narrowing accounting for the hardness of coliisions based on broadener and perturber.  
* Added CIA functionality.  Currently, can manually enter CIA for each file
* Added capability pass non-fitting columns through the fit definition.  This will help for future addition of bandwide functions or for ease of transforming fit outputs to reported results.
* Add ability to weight the spectra
* Added beyond HITRAN molecule capabilities, including update of isotope list.
* Added preliminary instrument line shape functionality.


BUG FIXES
* Changed the structure of MATS, so that each class is in its own python file
* Changed all constant values to be consistent with CODATA values and hardwired the values to definition in the utility script
* Added warning for floating parameters with an initial guess equal to 0.
* Changed the etalon_freq variables to etalon_period, since that is how it is coded.
* Added checks for molecule consistency between dataset and parameter line list.
* Ability to simulate at infinite SNR
* Indexing error for the spectrum number in baseline parameters
* Changed initial baseline guess to 0

MATS Version 1
++++++++++++++
published on 12/27/19


