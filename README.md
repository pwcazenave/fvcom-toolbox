fvcom-toolbox
=============

The fvcom-toolbox is a collection of MATLAB and some FORTRAN90 scripts for the purpose of preparing and postprocessing data from the Finite Volume Community Ocean Model (FVCOM). These include:

1. Scripts for preparing input files for FVCOM (wind forcing, open boundary forcing, river forcing, etc.)
2. Scripts for converting meshes from SMS to FVCOM
3. Scripts for postprocessing FVCOM data using MATLAB
4. Scripts for preparing data for the unstructured SWAN model

Notes:

(0) Some third-party MATLAB toolboxes are required for some functions:

* The Tidal Model Driver available at http://polaris.esr.org/ptm_index.html.
* The air-sea toolbox available at http://woodshole.er.usgs.gov/operations/sea-mat/air_sea-html/index.html.
* The OPeNDAP toolbox (for versions of MATLAB older than 2011b) available at http://www.opendap.org/pub/contributed/source/ml-toolbox/.

(1) The html based documentation is generated using m2html and is available with the download (see doc/index.html)

(2) The code was originally maintained at a Google Code repository (http://code.google.com/p/fvcom-toolbox/). This repository was used between Sep, 2010 (initial commit) and July, 2013 when it was moved to github (https://github.com/GeoffCowles/fvcom-toolbox). Commit history was not maintained during the move as substantial revisions had been made to the code by Plymouth Marine Laboratory members outside of version control. The github trunk includes most of these changes noted in the file headers and noted in the file PML_ChangeLog.txt. Although the development from the version included in FVCOM (github.com/GeoffCowles/fvcom-toolbox) is discontinuous from the version worked on by PML (and NOC and others), a commit history is available at http://gitlab.ecosystem-modelling.pml.ac.uk/pica/fvcom-toolbox or https://github.com/pwcazenave/fvcom-toolbox.

The PML version of the toolbox includes tagged releases, which can be downloaded as standalone (and thus relatively stable) versions. See the PML_ChangeLog.txt for details. Links to the direct downloads are:

- v20150319: https://github.com/pwcazenave/fvcom-toolbox/releases/tag/20150319
- v20141017: https://github.com/pwcazenave/fvcom-toolbox/releases/tag/20141017
- v20140728: https://github.com/pwcazenave/fvcom-toolbox/releases/tag/20140728
- v20140423: https://github.com/pwcazenave/fvcom-toolbox/releases/tag/20140423
- v20140131: https://github.com/pwcazenave/fvcom-toolbox/releases/tag/20140131
- v20130917: https://github.com/pwcazenave/fvcom-toolbox/releases/tag/20130917
- v20130719: https://github.com/pwcazenave/fvcom-toolbox/releases/tag/20130719
- v20130521: https://github.com/pwcazenave/fvcom-toolbox/releases/tag/20130521
- v20130403: https://github.com/pwcazenave/fvcom-toolbox/releases/tag/20130403
- v20130204: https://github.com/pwcazenave/fvcom-toolbox/releases/tag/20130204

To download the version included in the FVCOM source code, go to https://github.com/GeoffCowles/fvcom-toolbox.

Todo:

If you are stuck for something to get started with on the toolbox, there are some parts which are in need of some attention:

(0) The original MATLAB code used a third-party netCDF toolbox to write FVCOM input files. This library is largely redundant as recent versions of MATLAB include this functionality by default. Some functions in the toolbox still use the old toolbox:
- fvcom_prepro:
    * add_sigma_forDT.m
    * add_var_FVCOM_river.m
    * example_FVCOM_tsobc.m
    * example_FVCOM_wind_ts.m
    * example_FVCOM_wind_ts_speed.m
    * example_init_lag.m
    * wrf2fvcom_U10V10.m
    * write_FVCOM_bedflag.m
- swan_scripts:
    * calc_tauw.m
    * swan2netcdf.m
- utilities:
    * calc_tauw.m
    * gridvecs.m

(1) We need more tests in the tests subdirectory!

(2) Port the few FORTRAN codes to MATLAB to make this a more portable toolbox.
