.. docs/source/engines.rst 

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

* use other implementations different that f2py.
* figure out how to create a wrapper to a C procedure using iso-c-bindings to use f2py to call C methods.


Common Physical Grounds
=======================


Terminology and units
---------------------

* **Extinction**, **Scattering** and **Absorption coefficients** will be noted by  :math:`\beta_{e}`,  :math:`\beta_{s}` and  :math:`\beta_{a}` respectively, and all they have dimensions :math:`[L^{-1}]`.
* The **single scattering albedo** (:math:`\tilde{\omega}`) is defined by the relation :math:`\tilde{\omega}\beta_{e}=\beta_{s}`.
* **Mass extinction coefficent** (:math:`\kappa_{e}`) is defined by :math:`\kappa_{e}\rho=\beta_{e}` and the units are :math:`[L^{2}M^{-1}]`  (area per unit mass). One reason this ammount is useful is because we can express the optical depth of an homogeneous cell (of heigth H and area A) at top by :math:`\tau_{e}=\kappa_{e}M/A` where M is the mass (with :math:`\rho=M/(AH)`). So having the mass, the area of a grid-cell and :math:`\kappa_{e}` we can estimate the optical depth.
* If we use the concentration :math:`\chi` rather than density, then :math:`\beta_{e}=\sigma_{e}\chi`, and :math:`\kappa_{e}m=\sigma_{e}` with `m` the mass per particle (for monodisperse distribution), and in this context the extinction efficiency (:math:`Q_{e}`) is defined by :math:`Q_{e}A=\sigma_{e}` 

The definition of :math:`Q_{e}` and previous amounts can be also applied to :math:`Q_{s}` and :math:`Q_{a}`, with the useful relations for **polydisperse distributions**:

.. math::

   \beta_{e} = \int_{0}^{\infty} n(r)Q_{e}(r)\pi r^{2}dr \\
   \beta_{s} = \int_{0}^{\infty} n(r)Q_{s}(r)\pi r^{2}dr \\
   \beta_{a} = \int_{0}^{\infty} n(r)Q_{a}(r)\pi r^{2}dr 

but even:

.. math::

    P(cos \theta)=\beta_{s}^{-1}\int_{0}^{\infty}n(r)Q_{s}(r)\pi r^{2} p(cos \theta, r)dr

and

.. math::

   g = \beta_{s}^{-1}\int_{0}^{\infty} n(r)Q_{s}(r)\pi r^{2}g(r)dr


Log-normal distribution
"""""""""""""""""""""""

Currently, this is the only distribution supported by our Mie-Scattering engine because it is the most widely used. Based on a given optical model with a lognormal distribution the properties for a given bin or size interval are calculated as (note that size bins limits of config file are multiplied by 1.e-6, also the value of r0):

.. math::

   n(r)d(lnr)=N_{tot}\frac{1}{\sqrt{2\pi}ln(\sigma_{g})}\exp\left[-\frac{1}{2}\left(\frac{ln(r)-ln(r_{0})}{\sigma_{g}}\right)^{2}\right]

It is important the d(lnr). Note that the effective radius would be:

.. math::

   r_{eff}=\frac{\int_{a}^{b}r(\pi r^{2})n(r)dr}{\int_{a}^{b}\pi r^{2}n(r)dr}

in the case of the lognormal distribution and for a=0 and :math:`b\arrow\infty`, we have:  :math:`r_{eff}=r_{0}e^{2.5(ln\sigma_{g})^{2}}`
 
* For each wavelength we proceed to an integration of the optical properties (see equations above) in a specific size interval [a,b]
    * For each bin we estimate 3 parameters of the distribution using size-bin limits, Nbin_points, r_0, sigma_g, Ntot and Ndis
    * Mass is also estimated based on these parameters and hygrospopic growth
    * Following equations above, we have the integrals below to be solved for each size interval [a,b] with a discrete quadrature (trapezoid method). 

.. math::

   \beta_{e}^{\lambda} = \int_{a}^{b} n(r)Q_{e}(r)\pi r^{2}dr \\
   \beta_{s}^{\lambda} = \int_{a}^{b} n(r)Q_{s}(r)\pi r^{2}dr \\
   \beta_{a}^{\lambda} = \int_{a}^{b} n(r)Q_{a}(r)\pi r^{2}dr 

The mie-core calculates :math:`Q_{e,s,a}` (as well as asymmetry parameter) therefore other methods that provide these quantities can be integrated in the code. Note that
in the mie-core the important parameter is the mie-size-parameter: :math:`x=\frac{2\pi r}{\lambda}=\frac{D\pi}{\lambda}`


In the case of a multimodal lognormal distribution with m modes:

.. math::

   n(r)d(lnr)=\sum_{i=1}^{m} N_{i,tot}\frac{1}{\sqrt{2\pi}ln(\sigma_{i,g})}\exp\left[-\frac{1}{2}\left(\frac{ln(r)-ln(r_{i,0})}{\sigma_{i,g}}\right)^{2}\right]

.. note:: How this is dscribed in the code
   
   In each configuration file which described a single aerosol species (either multimodal or unimodal) the parameters of the above equation are:
   
   lognormal.r0       
   lognormal.sigma_g 
   lognormal.Ntot

   if they are a single value, then we have a unimodal distribution, if they are a list of values then we have a multimodal distribution. For example, in the case
   of black carbon models we have only one single mode with values:

   lognormal.r0 = 0.0118
   lognormal.sigma_g = 2.0
   lognormal.Ntot = 1.0
   
   then we are describing: 
   
     
   :math:`1.0\frac{1}{\sqrt{2\pi}ln(2.0)}\exp\left[-\frac{1}{2}\left(\frac{ln(r)-ln(0.0118)}{2}\right)^{2}\right]`


The code also calculates the phase function, and the scattering matrix, or more specifically the Phase matrix. For spherical particles the Mueller matrix has 4 different components
F11, F12, F33 and F34. It is a 4x4 matrix for each angle, but F44=F33, F43=-F34, F21=F12 and F32=F23=F13=F31=F41=F42=0. The element M11 is basically the phase function, and F12/F11 is
related to the linear polarization. 

.. math::

   F11 = \frac{1}{2}\left( S_{1}\dot\S_{1}^{*}+S_{2}\dot\S_{2}^{*} \right)

.. math::

   F12 = \frac{1}{2}\left( S_{1}\dot\S_{1}^{*}-S_{2}\dot\S_{2}^{*} \right)

.. math::

   F33 = \frac{1}{2}\left( S_{1}\dot\S_{2}^{*}+S_{2}\dot\S_{1}^{*} \right)

.. math::

   F34 = \frac{i}{2}\left( S_{1}\dot\S_{2}^{*}-S_{2}\dot\S_{1}^{*} \right)

The normalization of the final matrix P11, P12, P33, P34 is done in such way that P11 has the same normalization than the phase function. Regarding the 
calculation for a size distribution, we perform an integration for these Pij like it is done for the phase function.
