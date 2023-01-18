
# License
  - This software is licensed by ECMWF with an APACHE license. Please check file LICENSE file in main directory.

# Contact Info

  - Purpose:  New code to create IFS aerosol optical information
  - Authors:  Ramiro Checa-Garcia
              contributions of other authors vie Mie Scattering Codes.
  - Contacts: ramiro.checa-garcia@ecmwf.int
  - TODO:     stored at TODO.org file
  - History:  open file docs.html for further information (after make build)

```
     v0.8     Translated from Julia version
     v0.9     Updated the format for config files, added makefile and testing framework.
```

# How to install

To have a working version of this tool the steps are:

```
git clone ssh://git@git.ecmwf.int/~parc/ecaeropt.git
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



