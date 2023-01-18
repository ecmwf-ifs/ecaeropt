


IFS calculations
================

For each IFS cycle we are preparing **setting files** based on **configuration files** to create a 
reproducible enviroment for aerosol optical properties in the IFS model. Here we describe the main
information contained in setting files and configuration files per cycle. The generated optical
properties file have informative important metadata.

+-------------------+----------------------------+-------------------------+----------------------+
| IFS cycle         | Setting file               | Config files            | Notes                |
+===================+============================+=========================+======================+
| CY46R1            | settings/IFS_CY46R1.toml   | data/config_toml_CY46R1 | make IFS-CY46R1      |
+-------------------+----------------------------+-------------------------+----------------------+
| CY47R1            |                            |                         |                      |
+-------------------+----------------------------+-------------------------+----------------------+
| CY48R1            | settings/IFS_CY48R1.toml   | data/config_toml_CY48R1 | make IFS-CY48R1      |
+-------------------+----------------------------+-------------------------+----------------------+
| CY49R1            | settings/IFS_CY49R1.toml   | data/config_toml_CY48R1 | make IFS-CY49R1      |
+-------------------+----------------------------+-------------------------+----------------------+
| CY49R1 _v2        | settings/IFS_CY49R1_v2.toml| data/config_toml_CY48R1 | make IFS-CY49R1-v2   |
+-------------------+----------------------------+-------------------------+----------------------+


IFS CY46R1
----------

The directory **data/config_toml_CY46R1** has all the configuration files for the different aerosols. Here
we summarise the main properties

Summary table
+++++++++++++

+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+
| Aerosol                  | Modes            | Output bins     | Refrac. Ind.        | Soluble   | Mixed    | Config-files             |
+==========================+==================+=================+=====================+===========+==========+==========================+
| Mineral Dust             | 1 logn modes     | 3 bins          |  3 optical models   | No        | No       | config_du1.toml          |
|                          |                  |                 |                     |           |          | config_du2.toml          |
|                          |                  |                 |                     |           |          | config_du3.toml          |
+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+
| Organic Aerosols         | mixture          |                 |                     | Yes       | External | config_organic_inso.toml |
|                          |                  |                 |                     |           |          | config_organic_waso.toml |
|                          |                  |                 |                     |           |          | config_organic_soot.toml |
+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+
| Black carbon             | 1 logn mode      | 1 bins          |  3 optical models   |           | No       | config_bc1.toml          |
|                          |                  |                 |                     |           |          | config_bc2.toml          |
|                          |                  |                 |                     |           |          | config_bc3.toml          |
+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+
| Sea Salt                 | 2 logn modes     | 3 bins          |  1 optical models   | Yes       | No       | config_ss.toml           |
+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+
| Sulfate                  | 1 logn mode      | 1 bins          |  1 optical model    | Yes       | No       | config_su.toml           |
+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+
| Secondary Organic A      | 1 logn mode      | 1 bins          |  1 optical model    |           | No       | config_SOA1.toml         |
+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+
| Secondary Organic B      | 1 logn mode      | 1 bins          |  1 optical model    |           | No       | config_SOA2.toml         |
+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+
| Nitrate fine             | 1 logn mode      | 1 bins          |  1 optical model    | Yes       | No       | config_ni_fine.toml      |
+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+
| Nitrate coarse           | 1 logn mode      | 1 bins          |  1 optical models   | Yes       | No       | config_ni_coarse.toml    |
+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+
| NH3                      | 1 logn mode      | 1 bin           |  1 optical model    | Yes       | External | config_ammonia_fine.toml |
+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+


Mineral Dust
++++++++++++

The mineral dust is defined by 1 lognormal mode and then translated into 3 bins for the calculation of optical
properties.

+------------+--------------------+---------------------+----------------+--------------+------------------+ 
| Modes      | N                  | :math:`\sigma_{g}`  | :math:`r_{0}`  | Density      | Refrac. Ind.     |
+============+====================+=====================+================+==============+==================+
| single     | 1.0                | 2.0                 | 0.29           | 2610 kg/m3   | 3 models         |
+------------+--------------------+---------------------+----------------+--------------+------------------+


+------------+--------------------+----------------------+
| Bins       | :math:`r_{min}`    | :math:`r_{max}`      | 
+============+====================+======================+
| small      | 0.03               | 0.55                 |
+------------+--------------------+----------------------+
| medium     | 0.55               | 0.9                  |
+------------+--------------------+----------------------+
| large      | 0.9                | 20.0                 |
+------------+--------------------+----------------------+


IFS CY47R1
----------

No differences with CY46R1


IFS CY48R1
----------

During this cycle there are new optical models, these optical models are not included by the setting file for CY48R1 althought the
configuration files are developed during this cycle. These new optical models are included in the aerosols optical properties `CY49R1_v2`.
For this reason we summarize the new optical models in the next IFS CY49R1 but here we are considering no differences with previous
versions.

IFS CY49R1 and IFS CY49R1_v2
----------------------------

The CY49R1 is equivalent to CY48R1 in terms of optical models, but the netcdf files have additional metadata.

- CY49R1    -> identical to CY48R1 but with new metadata
- CY49R1_v2 -> new metadata and also new optical models


Summary Table
+++++++++++++

+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+
| Aerosol                  | Modes            | Output          | Refrac. Ind.        | Soluble   | Mixed    | Config-files             |
+==========================+==================+=================+=====================+===========+==========+==========================+
| Mineral Dust             | 1 logn           | 3 bins          |  3 optical models   | No        | No       | du1                      |
|                          |                  |                 |                     |           |          | du2                      |
|                          |                  |                 |                     |           |          | config_du3.toml          |
|                          | 4 logn           | 3 bins          |  1 optical model    |           |          | config_du4.toml          |
+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+
| Organic Aerosols         | mixture          |                 |                     | Yes       | External | config_organic_inso.toml |
|                          |                  |                 |                     |           |          | config_organic_waso.toml |
|                          |                  |                 |                     |           |          | config_organic_soot.toml |
+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+
| Organic Matter           | 1 logn           | 1 bin           |                     | Yes       | No       | config_om.toml           |
+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+
| Black carbon             | 1 logn           | 1 bins          |  3 optical models   |           | No       | config_bc1.toml          |
|                          |                  |                 |                     |           |          | config_bc2.toml          |
|                          |                  |                 |                     |           |          | config_bc3.toml          |
+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+
| Sea Salt                 | 2 logn           | 3 bins          |  1 optical models   | Yes       | No       | config_ss.toml           |
+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+
| Sulfate  I               | 1 logn mode      | 1 bins          |  1 optical model    | Yes       | No       | config_su.toml           |
+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+
| Sulfate  II              | 1 logn mode      | 1 bins          |  1 optical model    | Yes       | No       | config_su1.toml          |
+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+
| Sulfate  III             | 1 logn mode      | 1 bins          |  1 optical model    | Yes       | No       | config_su2.toml          |
+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+
| Secondary Organic A      | 1 logn mode      | 1 bins          |  1 optical model    |           | No       | config_SOA1.toml         |
+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+
| Secondary Organic B      | 1 logn mode      | 1 bins          |  1 optical model    |           | No       | config_SOA2.toml         |
+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+
| Nitrate fine             | 1 logn mode      | 1 bins          |  1 optical model    | Yes       | No       | config_ni_fine.toml      |
+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+
| Nitrate coarse           | 1 logn mode      | 1 bins          |  1 optical models   | Yes       | No       | config_ni_coarse.toml    |
+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+
| NH3                      | 1 logn mode      | 1 bin           |  1 optical model    | Yes       | External | config_ammonia_fine.toml |
+--------------------------+------------------+-----------------+---------------------+-----------+----------+--------------------------+


Mineral Dust New
++++++++++++++++


The mineral dust is defined by 1 lognormal mode and then translated into 3 bins for the calculation of optical
properties. The refractive index is a combination of different sources: Remy, Balkanski-2007 and Di Biaggio-2017,
the size distribution is derived from Ryder et al. (which is mostly derived from air-craft measurements).

+------------+--------------------+---------------------+------------------+------------+---------------+ 
| Modes      | N                  | :math:`\sigma_{g}`  | :math:`r_{0}`    | Density    | Refrac. Ind.  |
+============+====================+=====================+==================+============+===============+
| Fine       | 391.0              |  2.0                | 0.05             | 2610 kg/m3 | composite     |
+------------+--------------------+---------------------+------------------+------------+---------------+
| Medium     | 8.390              |  1.18               | 0.42             | 2610 kg/m3 | composite     |
+------------+--------------------+---------------------+------------------+------------+---------------+
| Coarse     | 11.6               |  1.93               | 0.79             | 2610 kg/m3 | composite     |
+------------+--------------------+---------------------+------------------+------------+---------------+
| Coarse     | 0.000138           |  1.53               | 16.2             | 2610 kg/m3 | composite     |
+------------+--------------------+---------------------+------------------+------------+---------------+

Organic Matter
++++++++++++++

Described in `config_om.toml` is using the reference Brown et al 2018 for the refractive index, and it is not anymore a mixture of aerosols species like organics aerosols before The optical model is named `Brown2018`

Sulfate II
++++++++++

Sulfate II which is described in `config_su1.toml` is like Sulfate I (`config_su.toml`) but without scaling of extinction (in other words, scaling=1.0)
The optical model is named `GACP-noscaling`.

Sulfate III
+++++++++++

This new sulfate has an updated distribution but same refractive index. Note that scaling is also equal to 1.0. The optical model is named `GACP-newPSD`


List of Optical Models
----------------------

+-------------------+----------------------------+-------------------------+----------------------+
| IFS cycle         | Config file                | Optical model           | Notes                |
+===================+============================+=========================+======================+
|CY48R1             | config_bc1.toml            | "OPAC"                  | OPAC==Hess1988(?)    |
+-------------------+----------------------------+-------------------------+----------------------+
|CY48R1             | config_bc2.toml            | "Stier2007"             |                      |
+-------------------+----------------------------+-------------------------+----------------------+
|CY48R1             | config_bc3.toml            | "Bond2006"              |                      |
+-------------------+----------------------------+-------------------------+----------------------+
|CY48R1             | config_du1.toml            | "Dubovik2002"           |                      |
+-------------------+----------------------------+-------------------------+----------------------+
|CY48R1             | config_ni_coarse.toml      | "GLOMAP"                | Default CAMS         |
+-------------------+----------------------------+-------------------------+----------------------+
|CY48R1             | config_su1.toml            | "GACP-NoScaling"        |                      |
+-------------------+----------------------------+-------------------------+----------------------+
|CY48R1             | config_du3.toml            | "Fouquart1987"          |                      |
+-------------------+----------------------------+-------------------------+----------------------+
|CY48R1             | config_ni_fine.toml        | "GLOMAP"                | Default CAMS         |
+-------------------+----------------------------+-------------------------+----------------------+
|CY48R1             | config_organic_inso.toml   | "OPAC"                  | NWP Climat.          |
+-------------------+----------------------------+-------------------------+----------------------+
|CY48R1             | config_organic_soot.toml   | "OPAC"                  | NWP Climat.          |
+-------------------+----------------------------+-------------------------+----------------------+
|CY48R1             | config_organic_waso.toml   | "OPAC"                  | NWP Climat.          |
+-------------------+----------------------------+-------------------------+----------------------+
|CY48R1             | config_SOA1.toml           | "Moise2015"             | Default CAMS         |
+-------------------+----------------------------+-------------------------+----------------------+
|CY48R1             | config_SOA2.toml           | "Moise2015"             | Default CAMS         |
+-------------------+----------------------------+-------------------------+----------------------+
|CY48R1             | config_ss.toml             | "OPAC"                  | Default CAMS         |
+-------------------+----------------------------+-------------------------+----------------------+
|CY48R1             | config_su.toml             | "GACP"                  |                      |
+-------------------+----------------------------+-------------------------+----------------------+
|CY48R1             | config_su2.toml            | "GACP-NewPSD"           | Default CAMS         |
+-------------------+----------------------------+-------------------------+----------------------+
|CY48R1             | config_ammonia_fine.toml   | "GACP"                  | Default CAMS         |
+-------------------+----------------------------+-------------------------+----------------------+
|CY48R1             | config_du4.toml            | "Composite"             | Default CAMS         |
+-------------------+----------------------------+-------------------------+----------------------+
|CY48R1             | config_om.toml             | "Brown2018"             | Default CAMS         |
+-------------------+----------------------------+-------------------------+----------------------+
|CY48R1             | config_du1_test.toml       | "Dubovik2002"           |                      |
+-------------------+----------------------------+-------------------------+----------------------+
|CY48R1             | config_du2.toml            | "Woodward2001"          | NWP Climat.          |
+-------------------+----------------------------+-------------------------+----------------------+

