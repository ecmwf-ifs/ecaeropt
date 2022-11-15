


IFS calculations
================

For each IFS cycle we are preparing setting files based on configuration files to create a reproducible
enviroment for aerosol optical properties in the IFS model. Here it is described the main information
contained in setting files and configuration files per cycle.





+------------------------+--------------------+-------------------------+-----------------+
| IFS cycle              | Setting file       | Config files            | Notes           |
+========================+====================+=========================+=================+
| CY46R1                 | settings/IFS_46R1  | data/config_toml_CY46R1 | make IFS-CY46R1 |
+------------------------+--------------------+-------------------------+-----------------+
| CY47R1                 |                    |                         |                 |
+------------------------+--------------------+-------------------------+-----------------+
| CY48R1                 |                    |                         |                 |
+------------------------+--------------------+-------------------------+-----------------+


IFS CY46R1
----------

The directory **data/config_toml_CY46R1** has all the configuration files for the different aerosols. Here
we summarise the main properties


+------------------------+--------------------+-------------------------+-----------------+
| Aerosol                | Modes              | Output bins             | Refractive Index|
+========================+====================+=========================+=================+
| Mineral Dust           |                    | 3 bins                  |  3 models       |
+------------------------+--------------------+-------------------------+-----------------+
| Organic Aerosols       |                    |                         |                 |
+------------------------+--------------------+-------------------------+-----------------+
|                        |                    |                         |                 |
+------------------------+--------------------+-------------------------+-----------------+

IFS CY48R1
----------

Setting file
~~~~~~~~~~~~
Configuration files
~~~~~~~~~~~~~~~~~~~
