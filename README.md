 

        ______ _____             ______ _____         ____  _____ _______   
       |  ____/ ____|      /\   |  ____|  __ \       / __ \|  __ \__   __|  
       | |__ | |   ______ /  \  | |__  | |__) |_____| |  | | |__) | | |     
       |  __|| |  |______/ /\ \ |  __| |  _  /______| |  | |  ___/  | |     
       | |___| |____    / ____ \| |____| | \ \      | |__| | |      | |     
       |______\_____|  /_/    \_\______|_|  \_\      \____/|_|      |_|      
     

# Description

Tool to calculate optical properties of atmospheric particles and store as netcdf files. It is designed to support different methodologies and input information for each specie uses human-readable files. The tool also can produce input netcdf files for the aerosols optical models used by IFS numerical weather prediction model and IFS/CAMs.

# License
  - This software is licensed by ECMWF with an APACHE license. Please check file LICENSE file in main directory.

# Contact Info

  - Authors:      Ramiro Checa-Garcia
                  contributions of other authors via Mie Scattering Codes.
  - Contacts:     ramiro.checa-garcia@ecmwf.int
  - Contributors:
      - Tims Stockfale: for testing/comments.
      - Robin Hogan: for testing/comments.
      - Samuel Remy: info of several optical models.

# External codes included/used:

  - Currently, the code uses the following Fortran codes:
    - **Mie-Code**: originally developed by Olivier Boucher. This code has been modified to be modular 
                and a ISO-C-binnding has been added. With the new modular structure tests for the 
                series of coefficients an and bn, as well as, single sphere has been added. 

# History

```
     v0.8     Translated from Julia version
     v0.9     Updated the format for config files, added makefile and testing framework.
     v1.0     Added IFS-CY49R1 with updated optical models from S. Remy
     v1.1     Added optical automatic plot functionality per specie of optical properties 
              in the setting files. 
     v1.16    Several fixes of issues reported by Tim.
```

  Check docs.html for further information (after make build)


# How to install

To have a working version of this tool the steps are:

```
git clone ssh://git@git.ecmwf.int/~parc/ecaeropt.git  # This download the master branch (stable version)
make build
```

Note that need python3 with the standard library and following external libraries, **numpy**, **toml** and **netcdf4** plus **nco tools**. In the server ATOS it is needed to use:
```
module load python3
module load nco
```

We recommend to run once the test to check everthing is working properly (this is just calculating for few wavelengths 3 cases: an externally mixed aerosol with 3 components, an soluble aerosol and an insolube aerosol.

```
./ecaeropt -i  # this show simple checks of mie-scatt code for monodisperse case.
make test      # tests including 3 kind of aerosols (about 10 minutes)
```

# How to run

## Create standard IFS model netcdf files

To build ifs standard files you may just type (script using SBATCH is implemented),

```
make ifs-CY46R1
make ifs-CY48R1
make ifs-CY49R1    
```

in both cases the log files with be stored in logs/ directory. The final netcdf will be stored at: `outputnc/store_ifsnc/`

**Important**: these `make` will clean (remove) ancillary files in `tmp` and `outputnc/*.nc` (but not `output/store_ifsnc/`)

## Custom calculations

For this you can use directly the tool, you can access to the description of the run options with

```
./ecaeropt -h 
```

remember that:

-  a configuration file (CONFIG-FILE) describes a calculation of optical properties of a *single aerosol species*
-  a setting files (SETTING-FILE) describes a more complex calculation: mixed aerosols, several single species and storage all in a single netcdf.
-  you can create also setting files to perform automatic tests. Everthing will be described in documentation.

# Documentation

You can have access to a more detailed documentation, by default `make build` will create a documentation and a link in the main folder named docs.html which you can open with any browser (ex. ` > firefox docs.html`)
You can also build documentation directly with `make docum`. Note that you need to have **sphinx** installed (the python library).

# How to contribute

You can use the code and report any issue or propose any new feature. Please, check also
the documentation you can open the file **docs/build/html/index.html**

If you want to contribute at coding level, there is a dev branch, where changes and implemented and evetually when full tested they are integrated in master branch. When a full set of features is included then we add a tag with the version. 


