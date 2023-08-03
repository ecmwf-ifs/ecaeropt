.. docs/source/news_info.rst 

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

Description
===========

**ecaeropt** is an tool to calculate aerosol optical properties with a broad set of options, configuration and methods. 
The main goal is to produce netcdf files with the estimated properties to be used by radiative transfers models, in particular *ecRad*.
Additional goals of the tool are reproducibily and full documentation of the datasets produced.


Last news anc history
*********************

**ecaeropt** python version is ready for use. Please report any issue to developers.


**2023-Jan-25**  (Plotting)
    It has been added an optional plot functionality to setting files to plot optical properties
    and refractive index of each optical model to be integrated on a given IFS netcdf file.
**2023-Jan-10**  (CY49R1)
    Added support to create CY49R1 files for ecRad 1.5 with new varibles in netcdf indicating
    code-string, optical_model and bins [currently under testing]
**2022-Nov-22**  (Documentation)
    Added API of the code to the documentation
**2022-Nov-21**  (ATOS, feature)
    Added make instance for IFS model using sbatch to be run at ECMWF-ATOS HPC
**2022-Nov-18** (data)
    Added data to create IFS-CY48R1 consistent with last CAMS added aerosols.
    TODO: validate the final output



