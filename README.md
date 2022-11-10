


# License
  - Probably the final license will be MIT or GPL. At this moment it is not
    decided. 


# Contact Info

  - Purpose:  New code to create IFS aerosol optical information

  - Authors:  Ramiro Checa-Garcia
              contributions of other authors vie Mie Scattering Codes.

  - Contacts: ramiro.checa-garcia@ecmwf.int

  - TODO:     stored at TODO.org file

  - History:  check directory docs/ for further information 

```
     v0.8     Translated from Julia version
     v0.9     Updated the format for config files
```

You need python3 with the following libraries:

- numpy
- toml
- netcdf4


in the server ATOS it is needed to use:
```
module load python3
```
other modules are used but they are part of standard library.


# How to run

You can access to the description of the run options with
```
./ecaeropt -h 
```

# Create documentation and how to contribute

You can use the code and report any issue or propose any new feature. Please, check also
the documentation:

```
cd docs
make html
```

and you can open the file **docs/build/html/index.html**


# Description

Tool to calculate optical properties of different kind of aerosols
and store the results as a netcdf file. It allows a calculation for single
and mixed aerosols for size distributions and bins.

## Approach for generating IFS aerosol optics

The tool is composed of several elements:
    - Read and parse configuration files and refractive index information.
    - Interface to several engines that calculate optical properties based
    on different codes (currently in fortran but extensible)
    - Storage of calculations in netcdfs with full metadata
    - A testing environment 

## Folder Structure
```
.
├── README.md             => this file
├── TODO.org              => todo items to improve code
├── ecaeropt              => main program executable [executable script]
├── build.sh              => exectuable script to build the engine libraries
│
├── data                  => data-files
│   ├── config_toml       => translation to TOML structured files
│   ├── config_toml_tests => new config TOML files for quick tests
│   ├── refr_idx          => refractive index files
│   └── wavelengths       => wavelengths files
│
├── docs                  => documentation using sphinx
│
├── engines
│   └── mie_Boucher_Bozzo => mie engine in Fortran code
│
├── outputnc              => default folder to store results
│   │
│   ├── dust_Dubovik_optics_IFS_2022-09-20.nc
│   ├── ...
│   └── sulfate_optics_IFS_2022-09-20.nc
│
├── settings              => complete settings files for calculations
│   ├── test_ifs.toml
│   └── test_single_mixing.toml
│
├── aeropt                => source code of tool (not engines)
│  
├── tmp                   => temporary files (testing)
│  
├── logs                  => log files (debug, building)
│  
└── tests                 => folder to store tests
    ├── config_toml_tests => new config TOML files for quick tests
    └── references        => folder to have references of tested netcdf
                             for future developments 

```
