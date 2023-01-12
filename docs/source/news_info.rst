

Description
===========

*ecaeropt* is an tool to calculate aerosol optical properties with a broad set of options, configuration and methods. 
The main goal is to produce netcdf files with the estimated properties to be used by radiative transfers models, in particular ecrad.
Additional goals of the tool are reproducibily and full documentation of the datasets produced.


Last news!
**********

*ecaeropt* python version is ready for use. Please report any issue to developers.


**2023-Jan-10**  (CY49R1)
    Added support to create CY49R1 files for ecRad 1.5 with new varibles in netcdf indicating
    code-string, optical_model and bins.
**2022-Nov-22**  (Documentation)
    Added API of the code to the documentation
**2022-Nov-21**  (ATOS, feature)
    Added make instance for IFS model using sbatch to be run at ECMWF-ATOS HPC
**2022-Nov-18** (data)
    Added data to create IFS-CY48R1 consistent with last CAMS added aerosols.
    TODO: validate the final output



