


How to install
==============

Dependencies
------------

Because **ecaeropt** is a python program which uses calculation engines or backends in other languages, mainly in Fortran. Therefore it's needed:

- python3 (currently tested with python3.8 to python3.10) with libraries
  + numpy
  + toml
  + netcdf4
  + f2py3    (usually installed with numpy)

- nco  (it's used for comparing netcdf between versions)
- fortran 90 compiler (currently only tested/used gfortran)


Installation
------------

The most easy way to install is use:

.. code-block:: bash

   make build

which create all needed directories, compile engine libraries and prepare documentation in html. Log files about the process are stored in logs folder.

Clean
-----


.. code-block:: bash

   make clean

It delete files in `tmp`, `logs` directories and the compiled libraries as well.



