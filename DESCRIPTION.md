


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
