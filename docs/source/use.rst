.. docs/source/use.rst 

   (C) Copyright 2022- ECMWF.
  
   This software is licensed under the terms of the Apache Licence Version 2.0
   which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 
   In applying this licence, ECMWF does not waive the privileges and immunities
   granted to it by virtue of its status as an intergovernmental organisation
   nor does it submit to any jurisdiction.

  Author:
     Ramiro Checa-Garcia. ECMWF
 
  Modifications:
     10-Dec-2022   Ramiro Checa-Garcia    1st. version


How to use
==========

The tool is an executable script **./ecaeropt**, you see the options with `./ecaeropt -h`. Around this tool we have also
created a set of automatic calculations using *make*.

Automatic calculations
----------------------

After build the tool using make build you can use make for specific prescribed calculations,
this calculations are using the submission `sbatch` for ECMWF-ATOS supercomputer. User of
other environments can easily run the equivalent scripts included in the scripts folder as:
`bash scripts/create_XXXX.sh`. Remember that many environments allow the use of tab after make 
to see all the options.

.. code-block:: bash

   make IFS-example    # runs a short calculation of IFS file with only two wavelengths for testing purposes.
   make IFS-CY46R1     # creates the standard IFS netcdf for IFS model cycle CY46R1
   make IFS-CY48R1     # creates the standard IFS netcdf for IFS model cycle CY48R1
   make IFS-CY49R1     # creates the standard IFS netcdf for IFS model cycle CY49R1
   make IFS-CY49R1_v2  # creates the standard IFS netcdf for IFS model cycle CY49R1 with new optical models

.. warning::
   We recommend NOT run two of these make cases as the same time. Currently, they will both store single species
   on the same folder, *outputnc*, and files can be overwritten. You can safely run several instances at the same
   time by changing output folder entry in the setting files manually if needed (before use make). The correspoding
   setting files for each make are easy to figure out but also there is table in section **IFS Calculations** of 
   this documentation.
 
   
Custom Calculations
-------------------
   
To prepare custom calculations usign ecaeropt you need to understand the use of **ecaeropt** you should be familiar with 3 concepts:

- configuration file
- setting file
- engine

Configuration file
++++++++++++++++++

It's a file in TOML which describes all the information needed to calculate the optical properties of an aerosols specie (with a given size distribution on a given size intervals). The TOML file has fields where external files for refractive index and wavelengths have to be provided, all other information should be explicit in the TOML file.

In *ecaeropt* the configuration file plus the referred files included on it is translated into an object that has all the information for the calculation.

Setting file
++++++++++++

It is another kind of TOML file, in this case it contains further information about how to create a set of netcdf files of optical properties for a set of configuration files. The final output can be also indicated as an IFS netcdf file used by ecrad.

The setting file file, when we include information about reference and testing files can be use also to compare calcuations (between versions, different engines, ...). Which is very useful in providing consistent versions of the code and ancillary files.


Engines
+++++++

Engines are the part of the code that according to configuration/setting files perform the actual calculations. They are expected to be wrappers to code in Fortran/C for faster calculations.




