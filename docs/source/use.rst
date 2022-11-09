


How to use
==========


The tool is an executable script **./ecaeropt**, you see the options with **./ecaeropt -h**


Concepts
--------

To understand the use of **ecaeropt** you should be familiar with 3 concepts:

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






