MATS Version Summary
====================

MATS Version 3.0
++++++++++++++++
published on 7/17/2023

**New Features**

- Power law temperature dependence for line-mixing
- Oxygen Theoretical CIA Model based on work of Karman et al
- Added ability to select TIPS version
- Updated to HAPI v1.2.2.0
- Added Example for treating static HITRAN LBL molecules
- Added ability to correct ideal gas law w/ compressability factor as calculated with NIST RefProp

MATS Version 2.1
++++++++++++++++
published on 12/03/2021

**Focus of version update is cleanup to make package git installable**

- renamed python files to lowercase, as this is the preferred python style.
- added __init__.py to include the most important functionality at the top level
- reformatted constants.  Added codata.py file with includes CODATA dictionary (with values and metadata) and CONSTANTS dictionary (with only values)
- Added linelistdata.py file to autoload csv files from MATS/LineList.
- Got rid of all `from package import *` commands.  Use explicit imports only
- Reworked example notebooks to reflect updates.
- Added INSTALL.md file to explain the install options.  


MATS Version 2
++++++++++++++
published on 6/10/2021

**New Features**

- Beta correction to the Dicke Narrowing accounting for the hardness of coliisions based on broadener and perturber.  
- Added CIA functionality.  Currently, can manually enter CIA for each file
- Added capability pass non-fitting columns through the fit definition.  This will help for future addition of bandwide functions or for ease of transforming fit outputs to reported results.
- Add ability to weight the spectra
- Added beyond HITRAN molecule capabilities, including update of isotope list.
- Added preliminary instrument line shape functionality.


**Bug Fixes**

- Changed the structure of MATS, so that each class is in its own python file
- Changed all constant values to be consistent with CODATA values and hardwired the values to definition in the utility script
- Added warning for floating parameters with an initial guess equal to 0.
- Changed the etalon_freq variables to etalon_period, since that is how it is coded.
- Added checks for molecule consistency between dataset and parameter line list.
- Ability to simulate at infinite SNR
- Indexing error for the spectrum number in baseline parameters
- Changed initial baseline guess to 0

MATS Version 1
++++++++++++++
published on 12/27/19


