
# License
  - This software is licensed by ECMWF with an APACHE license. Please check file LICENSE file in main directory.

# External codes included/used:

  - Currently, the code uses the following Fortran codes:
    - **Mie-Code**: originally developed by Olivier Boucher. This code has been modified to be modular 
                and a ISO-C-binnding has been added. With the new modular structure tests for the 
                series of coefficients an and bn, as well as, single sphere has been added. 
    - **T-Matrix**: code developed by Dr. Michael I. Mishchenko. This code was programmed in F77. It has been also modified
                to Fortran 90. The COMMON and GOTO statements has been removed.

# Contact Info

  - Purpose:  New code to create IFS aerosol optical information
  - Authors:  Ramiro Checa-Garcia
              contributions of other authors vie Mie Scattering Codes.
  - Contacts: ramiro.checa-garcia@ecmwf.int
  - TODO:     stored at TODO.org file
  - History:  open file docs.html for further information (after make build)
  - **Thanks!**:   Tim Stockdale and Robin Hogan for testing/comments.

```
     v0.8     Translated from Julia version
     v0.9     Updated the format for config files, added makefile and testing framework.
     v1.0     Added IFS-CY49R1 with updated optical models from S. Remy
     v1.1     Added optical automatic plot functionality per specie of optical properties 
              in the setting files. 
     v1.16    Several fixes of issues reported by Tim.
```

# How to install

To have a working version of this tool the steps are:

```
git clone ssh://git@git.ecmwf.int/~parc/ecaeropt.git  # This download the master branch (stable version)
make build
```

Note that need python3 with the standard library and following external libraries, numpy, toml and netcdf4 plus nco tools. In the server ATOS it is needed to use:
```
module load python3
module load nco
```

We recommend to run one the test to check everthing is working properly (this is just calculating for few wavelengths 3 cases: an externally mixed aerosol with 3 components, an soluble aerosol and an insolube aerosol.

```
./ecaeropt -i  # this show simple checks of mie-scatt code for monodisperse case.
make test      # tests including 3 kind of aerosols (about 10 minutes)
```

# How to run

## Create standard IFS model netcdf files

To build ifs standard files you may just type (script using SBATCH is implemented (openmp is a work in progress for mieBB),

```
make ifs-CY46R1
make ifs-CY48R1
make ifs-CY49R1     # which same optical models that CY48R1 but metadata compatible with CY49R1
make ifs-CY49R1_v2  # like CY49R12 but this has additional optical modes 
```
in both cases the log files with be stored in logs/ directory. The final netcdf will be stored at: `outputnc/store_ifsnc/`

**Important**: these `make` will clean (remove) ancillary files in `tmp` and `outputnc/*.nc` (but not `output/store_ifsnc/`)

## Custom calculations

For this you can use directly the tool, you can access to the description of the run options with

```
./ecaeropt -h 
```

remember that:

-  a configuration file (CONFIG-FILE) describes a calculation of optical properties of a single aerosol species
-  a setting files (SETTING-FILE) describes a more complex calculation: mixed aerosols, several single species and storage all in a single netcdf.
-  you can create also setting files to perform automatic tests. Everthing will be described in documentation.

# Documentation

You can have access to a more detailed documentation, by default `make build` will create a documentation and a link in the main folder named docs.html which you can open with any browser (ex. ` > firefox docs.html`)
You can also build documentation directly with `make docum`.

# How to contribute

You can use the code and report any issue or propose any new feature. Please, check also
the documentation you can open the file **docs/build/html/index.html**

If you want to contribute at coding level, there is a dev branch, where changes and implemented and evetually when full tested they are integrated in master branch. When a full set of features is included then we add a tag with the version. 

