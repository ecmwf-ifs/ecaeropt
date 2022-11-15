

Engines Implemented
*******************


Here we describe the use and physical grounds of the engines implemented. We aim here to include also those coding and mathematical details useful to understand the used engines.



Design
======

Interface
---------

The interface translate the language of aerosol configuration to the parameters needed by the actual engine. It is also required a flag that activated a debug.

.. code-block:: python

   def interface_new(aerconf, debug=False):

        # (1) transform aerconf to information for engine
        # (2) calling to engine 
        # (3) transform engine output into aeropt object

        return aeropt

Internals
---------

Most easy implementation is using **f2py**, however the mie-Boucher-Bozzo engine has a new wrapper that shows also how to create an interface using iso-c-bindings for fortran. This kind of code allow two things:
- use other implementations different that f2py.
- figure out how to create a wrapper to a C procedure using iso-c-bindings to use f2py to call C methods.


Common Physical Grounds
=======================




Terminology and units
---------------------

* Extinction, Scattering and Absorption coefficients will be noted by  :math:`\beta_{e}`,  :math:`\beta_{s}` and  :math:`\beta_{a}` respectively, and all they have dimensions :math:`[Lˆ{-1}]`.
  The single scattering albedo (:math:`\tilde{\omega}`) is defined by the relation :math:`\tilde{\omega}\beta_{e}=\beta_{s}`.
* Mass extinction coefficent (:math:`\kappa_{e}`) is defined by :math:`\kappa_{e}\rho=\beta_{e}` and the units are :math:`[Lˆ{2}Mˆ{-1}]`  (area per unit mass). One reason this ammount is useful is because we can express the optical depth of an homogeneous cell (of heigth H and area A) at top by :math:`\tau_{e}=\kappa_{e}M/A` where M is the mass (with :math:`\rho=M/(AH)`). So having the mass and the area of a grid-cell and :math:`\kappa_{e}` we can estimate the optical depth.
* If we use the concentration X rather than density, then :math:`\beta_{e}=\sigma_{e}X`, and :math:`\kappa_{e}m=\sigma_{e}` with m the mass per particle (for monodisperse distribution), and in this context the extinction efficiency (:math:`Q_{e}`) is defined by :math:`Q_{e}A=\sigma_{e}` 

The definition of :math:`Q_{e}` and previous amounts can be also applied to :math:`Q_{s}` and :math:`Q_{a}`, with the useful relations for polydisperse distributions:

.. math::

   \beta_{e} = \int_{0}^{\infty} n(r)Q_{e}(r)\pi rˆ{2}dr \\
   \beta_{s} = \int_{0}^{\infty} n(r)Q_{s}(r)\pi rˆ{2}dr \\
   \beta_{a} = \int_{0}^{\infty} n(r)Q_{a}(r)\pi rˆ{2}dr 

but even:

.. math::

    P(cos \theta)=\beta_{s}^{-1}\int_{0}^{\infty}n(r)Q_{s}(r)\pi rˆ{2} p(cos \theta, r)dr

and

.. math::

   g = \beta_{s}^{-1}\int_{0}^{\infty} n(r)Q_{s}(r)\pi rˆ{2}g(r)dr


Log-normal distribution
"""""""""""""""""""""""


Size parameter
""""""""""""""

:math:`x=\frac{2\pi r}{\lambda}=\frac{D\pi}{\lambda}`




