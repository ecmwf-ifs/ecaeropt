[comment]: ################################################################################
[comment]: # DESCRIPTION.md
[comment]: #
[comment]: #   (C) Copyright 2022- ECMWF.
[comment]: #  
[comment]: #   This software is licensed under the terms of the Apache Licence Version 2.0
[comment]: #   which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
[comment]: # 
[comment]: #   In applying this licence, ECMWF does not waive the privileges and immunities
[comment]: #   granted to it by virtue of its status as an intergovernmental organisation
[comment]: #   nor does it submit to any jurisdiction.
[comment]: #
[comment]: #  Author:
[comment]: #     Ramiro Checa-Garcia. ECMWF
[comment]: # 
[comment]: #  Modifications:
[comment]: #     10-Dec-2022   Ramiro Checa-Garcia    1st. version
[comment]: #
[comment]: #################################################################################


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
├── README.md              => GENERAL INFO
├── DESCRIPTION.md         => THIS file: list of original directories 
├── TODO.org               => todo items to improve code
├── ecaeropt               => main program executable [executable script]
├── Makefile               => easy acces to scripts to build engines and doc, 
│                             and run several cases.
├── data                   => data-files
│   ├── config_toml        => translation to TOML structured files
│   ├── config_toml_CY46R1 => config TOML files for CY48R1 cycle
│   ├── config_toml_CY48R1 => new TOML files for CY48R1 cycle
│   ├── config_toml_CY49R1 => testing TOML files for next cycle
│   ├── config_toml_new    => new config TOML files for quick test
│   ├── refr_idx           => refractive index files
│   ├── non_sphere_scaling => scaling factor for non-spherical particles
│   └── wavelengths        => wavelengths files
│
├── docs                   => documentation using sphinx
│
├── engines
│   └── mie_Boucher_Bozzo  => mie engine in Fortran code
│
├── outputnc               => default folder to store single netcdf results
│   ├── dust_Dubovik_optics_IFS_2022-09-20.nc
│   ├── ...
│   └── sulfate_optics_IFS_2022-09-20.nc
│
├── outputplt              => default folder to store plots
│   ├── ...
│   └── ....png
│
├── settings              => complete settings files for calculations
│   ├── IFS_CY46R1.toml
│   ├── ...
│   └── IFS_CY49R1_v3.toml
│
├── scripts               => script used for gnu make system
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
