.. docs/source/install.rst 

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



