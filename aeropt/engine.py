

##########################################################################################
# aeropt/aer.py
#
# (C) Copyright 2022- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#
#
# Author:
#    Ramiro Checa-Garcia. ECMWF
#
# Modifications:
#    05-Nov-2022   Ramiro Checa-Garcia    1st. version
#
#                                                                                         
# Info: 
#      This file is likely to be extended with new functions, one per inferface/wrapper   
#      to external libraries.                                                              
#                                                                                         
#      FUNCTIONS   
#       * test_approx   : Implemented for future tests 
#       * test_engine   : Test engine (mie-spherical) for single spheres
#       * interface_mie : Interface for the external mie code engine.
###########################################################################################


import sys
import numpy as np
import pytest

import aeropt.aer as aer
from engines.mie.libs.mie import mie as mie


def test_approx(val, ref, tol, info, passed, failed):

    if val == pytest.approx(ref, tol):
        print("             Passed     : ",info)
        passed=passed+1
    else: 
        print("             Failed     : ",info)
        print("              calculated:", val) 
        print("              reference :", ref)
        failed=failed+1
    
    return passed, failed

def test_engine(enginename, tol=1.e-7):
    """
    Function to tests a given interface for simple cases:
      - single sphere non-absorbing
      - single sphere absorbing

    Args:

    """
    ok  =0
    fail=0

    fact = 1.e06
    lambda_out    = np.array([0.6328])
    bins_min      = np.array([4.9646*lambda_out[0]/(2.0*np.pi)])*fact
    bins_max      = np.array([bins_min])
    sigma_g       = np.array([0.0])
    Ntot          = np.array([1.0])
    r0            = bins_min
    rho           = 1000.0
    rh_tab        = np.array([1.0])
    rh_growth     = np.array([1.0])
    znr_tab       = np.array([1.50])
    zni_tab       = np.array([0.0])
    lambda_tab    = np.array([0.6328])
    angles        = np.array([0.0,90.0,180.0])
    cutoff_radius = 190.5

   
    print("          Non-absorbing sphere...(tolerance +-",tol,")")
    out_noabs = mieBB.main_mie_aerosols(lambda_out, bins_min,
                                  bins_max, sigma_g,
                                  r0, Ntot, rho,
                                  rh_tab,rh_growth,
                                  cutoff_radius, 3, lambda_tab,
                                  znr_tab, zni_tab,
                                  angles)
    ok, fail = test_approx(out_noabs[8][0], 3.89618111, tol, "Qext", ok, fail)
    ok, fail = test_approx(out_noabs[8][1], 3.89618111, tol, "Qsca", ok, fail)
    ok, fail = test_approx(out_noabs[8][2], 0.0       , tol, "Qabs", ok, fail)
    ok, fail = test_approx(out_noabs[8][3], 0.70765370, tol, "gpar", ok, fail)

    print("          Absorbing sphere...    (tolerence +-",tol,")")
    zni_tab_abs  = np.array([+0.01])
    out_abs = mieBB.main_mie_aerosols(lambda_out, bins_min,
                                  bins_max, sigma_g,
                                  r0, Ntot, rho,
                                  rh_tab,rh_growth,
                                  cutoff_radius, 3, lambda_tab,
                                  znr_tab, zni_tab_abs,
                                  angles)

    ok, fail = test_approx(out_abs[8][0], 3.80065414, tol, "Qext", ok, fail)
    ok, fail = test_approx(out_abs[8][1], 3.54171034, tol, "Qsca", ok, fail)
    ok, fail = test_approx(out_abs[8][2], 0.25894380, tol, "Qabs", ok, fail)
    ok, fail = test_approx(out_abs[8][3], 0.73124641, tol, "gpar", ok, fail)

    print("          Passed ", ok  ," calculated variables")
    print("          Failed ", fail," calculated variables")

    return

def interface_mie(aerconf, logfile, debug=False, mix=0, verbose=0):
    """
    Interface to calculate optical properties given an object of aerosol configurations

    Args:
        aerconf (aer) : object of aerosol configuration.
        logfile (str) : string log filename in case it's used to debug
        debug   (bool): True if debug is activated, False otherwise
        mix     (int) : not used
        verbose (int) : not used

    Return:
        aeropt  (aeropt) : object of optical properties

    TODO: implement properly the degub info output inside fortran code.
    """
    # This is now to create a shared library (for Julia, and probably cpython)
    # gfortran -shared -fPIC -o mie_Boucher_Boss.so main_mie_aersols.f90 interp_ri.f90 parkind1.F90

    # This part of the code is unclear in the original code, why this is done. So cutoff_radius=7.5 is
    # hard-coded
    if mix >= 3:
        cutoff_radius = 7.5
    else:
        cutoff_radius = np.array(aerconf.bins_max) + 100.0  # at this point I don't understand well this statement

    #cutoff_radius = 7.5

    print("                     calculation of optical properties")
    print("                     engine: [Boucher-Bozzo code]")
    print("                     ...calculating optical properties...")

    if debug==False:
        idebug = 0
    else:
        idebug = 2

    logdate = int(logfile.split("/")[1].split(".")[0][0:6])
    logtime = int(logfile.split("/")[1].split(".")[0][6:13])

    # With this two number is possible to reconstruct the filename in Fortran.
    # Rather than add to the default log file we can create
    # tmp/datetime_mieBB.log

    sys.stdout.flush()
    out = mie.main_mie_aerosols(aerconf.lambda_out, aerconf.bins_min,
                                aerconf.bins_max, aerconf.sigma_g,
                                aerconf.r0, aerconf.Ntot, aerconf.rho,
                                aerconf.rh_tab, aerconf.rh_growth,
                                cutoff_radius, idebug, aerconf.ri_lambdatab,
                                aerconf.znr_tab, aerconf.zni_tab,
                                aerconf.angles )

    aeropt = aer.mie_to_aeropt(aerconf, out, "mie_Boucher_Bozzo")
    
    return aeropt
