Generating Parameter Line lists 
===============================

Parameter Line list Overview
++++++++++++++++++++++++++++

The MATS program uses a spectroscopic line list that varies from the traditional HITRAN and HAPI format to accomodate the full HTP parameterization including temperature dependences for all parmeters (other than the correlation parameter and nearest neighbor line-mixing).  The linelist format has rows that correspond to each line and columns that correspond to the various line parameters.  The necessary column headers are below, where there needs to be a column for each parameter and each diluent being fit (ie. air, self) and the line-mixing term needs to include the nominal temperature as a suffix.  
	* molec_id: HITRAN molecule id number
	* local_iso_id: HITRAN local isotope id number
	* nu: line center (cm-1)
	* sw : spectral line intensity (cm−1/(molecule⋅cm−2)) at Tref=296K
	* elower: The lower-state energy of the transition (cm-1)
	* gamma0_diluent: collisional half-width (cm−1/atm) at Tref=296K and reference pressure pref=1atm for a given diluent
	* n_gamma0_diluent: coefficient of the temperature dependence of the half width
	* delta0_diluent: pressure shift (cm−1/atm) at Tref=296K and pref=1atm of the line position with respect to line center
	* n_delta0_diluent:  the coefficient of the temperature dependence of the pressure shift
	* SD_gamma_diluent: the ratio of the speed dependent width to the collisional half-width at reference temperature and pressure
	* n_gamma2_diluent: the coefficient of the temperature dependence of the speed dependent width NOTE: This is the temperature dependence of the speed dependent width not the ratio of the speed dependence to the half-width
	* SD_delta_diluent:  the ratio of the speed dependent shift to the collisional shift at reference temperature and pressure
	* n_delta2_diluent: he coefficient of the temperature dependence of the speed dependent shift NOTE: This is the temperature dependence of the speed dependent shift not the ratio of the speed dependence to the shift
	* nuVC_diluent: dicke narrowing term at reference temperature and pressure
	* n__nuVC_diluent:  coefficient of the temperature dependence of the dicke narrowing term
	* eta_diluent:  the correlation parameter for the VC and SD effects
	* y_diluent_nominaltemp: linemixing term (as currently written this doesn't have a temperature dependence, so different column for each nominal temperature)

Generating Parameter Line list from HITRAN
++++++++++++++++++++++++++++++++++++++++++

The parameter line list .csv file can be generated manually, but `this code <https://github.com/usnistgov/MATS/blob/master/MATS/HITRAN_to_Dataframe.ipynb>`_ generates a parameter list using `HITRAN Application Programming Interface (HAPI) <https://hitran.org/hapi/>`_.  The overview below walks through this notebook.  

Imports
-------
Import the numpy, pandas, os, and sys packages.  Set pandas dataframe to show all of the rows.

.. code:: ipython3

   import numpy as np
   import pandas as pd
   pd.set_option("display.max_rows", 101)
   import os, sys

Optional imports include matplotlib and seaborn for plotting.   

.. code:: ipython3

   import matplotlib.pyplot as plt
   import seaborn as sns
   sns.set_style("whitegrid")
   sns.set_style("ticks")
   sns.set_context("poster")
   
Add location of hapi.py to system path

.. code:: ipython3

   sys.path.append(r'C:\Users\Documents\Python Scripts\HAPI')#Add hapi.py folder location to system path
   from hapi import *
   
Set Working File Location

.. code:: ipython3

   spec_path = r'C:\Users\Documents\Python Scripts\HAPI' # Location of the Summary Data File
   os.chdir(spec_path)

Molecule and Isotope Information
---------------------------------

HAPI provids a dictionary called ISO_ID, where the global isotope id acts as the key, with a sub dictionary containing the molecule id, local isotope id, isotope name, relative abundance, mass, and molecule name.  The following command will generate the dictionary as the output.  This provides the necessary information for interfacing with HAPI to pull-down HITRAN infomation.

.. code:: ipython3

   print_iso_id()
 
.. parsed-literal:: 


   The dictionary "ISO_ID" contains information on "global" IDs of isotopologues in HITRAN

      id            M    I                    iso_name       abundance       mass        mol_name
       1     :      1    1                     H2(16O)    0.9973170000  18.010565             H2O
       2     :      1    2                     H2(18O)    0.0019998300  20.014811             H2O
       3     :      1    3                     H2(17O)    0.0003720000  19.014780             H2O
       4     :      1    4                     HD(16O)    0.0003106900  19.016740             H2O
       5     :      1    5                     HD(18O)    0.0000006230  21.020985             H2O
       6     :      1    6                     HD(17O)    0.0000001160  20.020956             H2O
     129     :      1    7                     D2(16O)    0.0000000242  20.022915             H2O
       7     :      2    1                 (12C)(16O)2    0.9842040000  43.989830             CO2
       8     :      2    2                 (13C)(16O)2    0.0110570000  44.993185             CO2
       9     :      2    3             (16O)(12C)(18O)    0.0039471000  45.994076             CO2
      10     :      2    4             (16O)(12C)(17O)    0.0007340000  44.994045             CO2
      11     :      2    5             (16O)(13C)(18O)    0.0000443400  46.997431             CO2
      12     :      2    6             (16O)(13C)(17O)    0.0000082500  45.997400             CO2
      13     :      2    7                 (12C)(18O)2    0.0000039573  47.998322             CO2
      14     :      2    8             (17O)(12C)(18O)    0.0000014700  46.998291             CO2
     121     :      2    9                 (12C)(17O)2    0.0000001368  45.998262             CO2
      15     :      2   10                 (13C)(18O)2    0.0000000450  49.001675             CO2
     120     :      2   11             (18O)(13C)(17O)    0.0000000165  48.001650             CO2
     122     :      2   12                 (13C)(17O)2    0.0000000015  47.001618             CO2
      16     :      3    1                      (16O)3    0.9929010000  47.984745              O3
      17     :      3    2             (16O)(16O)(18O)    0.0039819400  49.988991              O3
      18     :      3    3             (16O)(18O)(16O)    0.0019909700  49.988991              O3
      19     :      3    4             (16O)(16O)(17O)    0.0007400000  48.988960              O3
      20     :      3    5             (16O)(17O)(16O)    0.0003700000  48.988960              O3
      21     :      4    1                 (14N)2(16O)    0.9903330000  44.001062             N2O
      22     :      4    2             (14N)(15N)(16O)    0.0036409000  44.998096             N2O
      23     :      4    3             (15N)(14N)(16O)    0.0036409000  44.998096             N2O
      24     :      4    4                 (14N)2(18O)    0.0019858200  46.005308             N2O
      25     :      4    5                 (14N)2(17O)    0.0003690000  45.005278             N2O
      26     :      5    1                  (12C)(16O)    0.9865400000  27.994915              CO
      27     :      5    2                  (13C)(16O)    0.0110800000  28.998270              CO
      28     :      5    3                  (12C)(18O)    0.0019782000  29.999161              CO
      29     :      5    4                  (12C)(17O)    0.0003680000  28.999130              CO
      30     :      5    5                  (13C)(18O)    0.0000222200  31.002516              CO
      31     :      5    6                  (13C)(17O)    0.0000041300  30.002485              CO
      32     :      6    1                     (12C)H4    0.9882700000  16.031300             CH4
      33     :      6    2                     (13C)H4    0.0111000000  17.034655             CH4
      34     :      6    3                    (12C)H3D    0.0006157500  17.037475             CH4
      35     :      6    4                    (13C)H3D    0.0000049203  18.040830             CH4
      36     :      7    1                      (16O)2    0.9952620000  31.989830              O2
      37     :      7    2                  (16O)(18O)    0.0039914100  33.994076              O2
      38     :      7    3                  (16O)(17O)    0.0007420000  32.994045              O2
      39     :      8    1                  (14N)(16O)    0.9939740000  29.997989              NO
      40     :      8    2                  (15N)(16O)    0.0036543000  30.995023              NO
      41     :      8    3                  (14N)(18O)    0.0019931200  32.002234              NO
      42     :      9    1                 (32S)(16O)2    0.9456800000  63.961901             SO2
      43     :      9    2                 (34S)(16O)2    0.0419500000  65.957695             SO2
      44     :     10    1                 (14N)(16O)2    0.9916160000  45.992904             NO2
      45     :     11    1                     (14N)H3    0.9958715000  17.026549             NH3
      46     :     11    2                     (15N)H3    0.0036613000  18.023583             NH3
      47     :     12    1                H(14N)(16O)3    0.9891100000  62.995644            HNO3
     117     :     12    2                H(15N)(16O)3    0.0036360000  63.992680            HNO3
      48     :     13    1                      (16O)H    0.9974730000  17.002740              OH
      49     :     13    2                      (18O)H    0.0020001400  19.006986              OH
      50     :     13    3                      (16O)D    0.0001553700  18.008915              OH
      51     :     14    1                      H(19F)    0.9998442500  20.006229              HF
     110     :     14    2                      D(19F)    0.0001150000  21.012505              HF
      52     :     15    1                     H(35Cl)    0.7575870000  35.976678             HCl
      53     :     15    2                     H(37Cl)    0.2422570000  37.973729             HCl
     107     :     15    3                     D(35Cl)    0.0001180050  36.982954             HCl
     108     :     15    4                     D(37Cl)    0.0000377350  38.980004             HCl
      54     :     16    1                     H(79Br)    0.5067800000  79.926160             HBr
      55     :     16    2                     H(81Br)    0.4930600000  81.924115             HBr
     111     :     16    3                     D(79Br)    0.0000582935  80.932439             HBr
     112     :     16    4                     D(81Br)    0.0000567065  82.930392             HBr
      56     :     17    1                     H(127I)    0.9998442500 127.912297              HI
     113     :     17    2                     D(127I)    0.0001150000 128.918575              HI
      57     :     18    1                 (35Cl)(16O)    0.7559100000  50.963768             ClO
      58     :     18    2                 (37Cl)(16O)    0.2417200000  52.960819             ClO
      59     :     19    1             (16O)(12C)(32S)    0.9373900000  59.966986             OCS
      60     :     19    2             (16O)(12C)(34S)    0.0415800000  61.962780             OCS
      61     :     19    3             (16O)(13C)(32S)    0.0105300000  60.970341             OCS
      62     :     19    4             (16O)(12C)(33S)    0.0105300000  60.966371             OCS
      63     :     19    5             (18O)(12C)(32S)    0.0018800000  61.971231             OCS
      64     :     20    1                H2(12C)(16O)    0.9862400000  30.010565            H2CO
      65     :     20    2                H2(13C)(16O)    0.0110800000  31.013920            H2CO
      66     :     20    3                H2(12C)(18O)    0.0019776000  32.014811            H2CO
      67     :     21    1                H(16O)(35Cl)    0.7557900000  51.971593            HOCl
      68     :     21    2                H(16O)(37Cl)    0.2416800000  53.968644            HOCl
      69     :     22    1                      (14N)2    0.9926874000  28.006147              N2
     118     :     22    2                  (14N)(15N)    0.0072535000  29.997989              N2
      70     :     23    1                 H(12C)(14N)    0.9851100000  27.010899             HCN
      71     :     23    2                 H(13C)(14N)    0.0110700000  28.014254             HCN
      72     :     23    3                 H(12C)(15N)    0.0036217000  28.007933             HCN
      73     :     24    1               (12C)H3(35Cl)    0.7489400000  49.992328           CH3Cl
      74     :     24    2               (12C)H3(37Cl)    0.2394900000  51.989379           CH3Cl
      75     :     25    1                    H2(16O)2    0.9949520000  34.005480            H2O2
      76     :     26    1                    (12C)2H2    0.9776000000  26.015650            C2H2
      77     :     26    2                (12C)(13C)H2    0.0219700000  27.019005            C2H2
     105     :     26    3                    (12C)2HD    0.0003045500  27.021825            C2H2
      78     :     27    1                    (12C)2H6    0.9769900000  30.046950            C2H6
     106     :     27    2              (12C)H3(13C)H3    0.0219526110  31.050305            C2H6
      79     :     28    1                     (31P)H3    0.9995328300  33.997238             PH3
      80     :     29    1            (12C)(16O)(19F)2    0.9865400000  65.991722            COF2
     119     :     29    2            (13C)(16O)(19F)2    0.0110834000  66.995083            COF2
     126     :     30    1                 (32S)(19F)6    0.9501800000 145.962492             SF6
      81     :     31    1                     H2(32S)    0.9498800000  33.987721             H2S
      82     :     31    2                     H2(34S)    0.0421400000  35.983515             H2S
      83     :     31    3                     H2(33S)    0.0074980000  34.987105             H2S
      84     :     32    1           H(12C)(16O)(16O)H    0.9838980000  46.005480           HCOOH
      85     :     33    1                     H(16O)2    0.9951070000  32.997655             HO2
      86     :     34    1                       (16O)    0.9976280000  15.994915               O
      87     :     36    1                 (14N)(16O)+    0.9939740000  29.997989             NOp
      88     :     37    1                H(16O)(79Br)    0.5056000000  95.921076            HOBr
      89     :     37    2                H(16O)(81Br)    0.4919000000  97.919027            HOBr
      90     :     38    1                    (12C)2H4    0.9773000000  28.031300            C2H4
      91     :     38    2              (12C)H2(13C)H2    0.0219600000  29.034655            C2H4
      92     :     39    1               (12C)H3(16O)H    0.9859300000  32.026215           CH3OH
      93     :     40    1               (12C)H3(79Br)    0.5013000000  93.941811           CH3Br
      94     :     40    2               (12C)H3(81Br)    0.4876600000  95.939764           CH3Br
      95     :     41    1           (12C)H3(12C)(14N)    0.9748200000  41.026549           CH3CN
      96     :     42    1                 (12C)(19F)4    0.9893000000  87.993616             CF4
     116     :     43    1                    (12C)4H2    0.9559980000  50.015650            C4H2
     109     :     44    1                H(12C)3(14N)    0.9646069000  51.010899            HC3N
     103     :     45    1                          H2    0.9996880000   2.015650              H2
     115     :     45    2                          HD    0.0003114320   3.021825              H2
      97     :     46    1                  (12C)(32S)    0.9396240000  43.971036              CS
      98     :     46    2                  (12C)(34S)    0.0416817000  45.966787              CS
      99     :     46    3                  (13C)(32S)    0.0105565000  44.974368              CS
     100     :     46    4                  (12C)(33S)    0.0074166800  44.970399              CS
     114     :     47    1                 (32S)(16O)3    0.9423964000  79.956820             SO3
     123     :     48    1                (12C)2(14N)2    0.9707524330  52.006148            C2N2
     124     :     49    1           (12C)(16O)(35Cl)2    0.5663917610  97.932620           COCl2
     125     :     49    2      (12C)(16O)(35Cl)(37Cl)    0.3622352780  99.929670           COCl2


Set-Up Variables for Line list
------------------------------
To generate a line list you will need to provide a tablename (str), an array containing the global isotope numbers of the molecules/isotopes that you are interested in, the minimum and maximum wavenumbers, and the minimum line intensity of interested. The example below would generate a HITRAN table named 'CO' that contains the all CO isotopes and the most abundant CO2 isotope in the spectral region between 6200 and 6500 cm-1 that have line intensities greater than 1e-30. 

.. code:: ipython3

   tablename = 'CO'
   global_isotopes = [26, 27, 28, 29, 30,31,7]
   wave_min = 6200 #7903.5#cm-1
   wave_max = 6500 #7904.5 #cm-1
   intensity_cutoff = 1e-30


Generate HITRAN and Initial Guess Line lists from HAPI Call
-----------------------------------------------------------
The next section of the example contains a function and function call teh output a MATS compatible line list.  NOTE: The line mixing term that the fitting script wants has a subscript with the nominal temperatures included in the dataset. I have been adding these columns by hand to the .csv by copying and pasting. Code can be updated to do this manually.


.. code-block:: python

   def HITRANlinelist_to_csv(isotopes, minimum_wavenumber, maximum_wavenumber, tablename = 'tmp', filename = tablename, temperature = 296): 
		
		"""Generates two .csv files generated information available from HTIRAN.  The first line list matches the information available from HITRAN (_HITRAN.csv) and the second supplements the HITRAN information with theoretical values and translates into MATS input format (_initguess.csv)

		Outline

		1. Gets a line list from HITRAN and saves all available parameters to filename_HITRAN.csv
		2. Goes through the data provided from HITRAN and collects the highest order line shape information.
		3.  Where there is missing information for the complete HTP linelist set to 0 or make the following substitutions
			- for missing diluent information fill values with air
			- set missing shift temperature dependences equal to 0 (linear temperature dependence)
			- calculate the SD_gamma based on theory 
			- set the gamma_2 temperature exponent equal to the gamma0 temperature exponent
			- set the delta_2 temperature exponent equal to the delta0 temperature exponent
			- set the dicke narrowing temperature exponent to 1
		4. Save the supplemented and MATS formatted HITRAN information as filename_initguess.csv


		Parameters
		----------
		isotopes : list
			list of the HITRAN global isotope numbers to include in the HAPI call
		minimum_wavenumber : float
			minimum line center (cm-1) to include in the HAPI call.
		maximum_wavenumber : float
			maximum line center (cm-1) to include in the HAPI call.
		tablename : str, optional
			desired name for table generated from HAPI call. The default is 'tmp'.
		filename : str, optional
			acts as a base filename for the .csv files generated. The default is tablename.
		temperature : float, optional
			Nominal temperature of interest.  HITRAN breaks-up the HTP line parameters into temperature regimes.  This allows for selection of the most approriate parameter information. The default is 296.

		Returns
		-------
		linelist_select : dataframe
			pandas dataframe corresponding to the HITRAN information supplemented by theoretical values/assumptions.
		filename_HITRAN.csv : .csv file
			file corresponding to available HITRAN information
		filename_initguess.csv : .csv file
			file corresponding to available HITRAN information supplemented by theory and assumptions in MATS format

		"""

