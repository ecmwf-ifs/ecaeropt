


How to contribute
=================

There are several ways to help to *ecaeropt*:

- The most easy way is just to run the code and report problems or request features, this can be done open an issue in the repository.
- You can also work in code improvements. Below you can read how we organise the code development:


Code Style
----------

* As much as possible standard PEP style for python. 
* Add docstrings to new methods/function/classes using the google docstrings syntax for python to be processed by the documentation system.
* For engines in fortran we recommend Fortran 90 or later, with similar ideas to: https://www.fortran90.org/src/best-practices.html 

Engines
-------

The engines need an interface that translates the current aerosol configuration into the details needed by their computational methods.
There are two scenarios, a wrapper that allows to call a procedure of Fortran, C etc, inside python (for example with f2py) or a python
code that creates the input for an external tool. The outputs have to be also introduced in the aerosol optical properties object.

.. IMPORTANT::
  
   We recommend/request to create tests example to validate the engine inside the tool, and ideally incorporate to the ``info
   mode`` which informs about the engines. For example, the mie scattering code is testing that a monodisperse distribution of spheres
   is providing right values (for absorbing and not absorbing spheres).

.. TODO::
  
   We are currently working in:
   (a) parallel and more structured version of Mie Code
   (b) new engine based on T-matrix code
   (c) Check and improve the legendre expansion calculations

Methods
-------


