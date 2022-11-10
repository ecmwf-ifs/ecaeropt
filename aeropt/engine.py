
 ###########################################################################################
 #                                                                                         #
 # aeropt/engine.py                                                                        #
 #                                                                                         #
 # author: Ramiro Checa-Garcia                                                             #
 # email:  ramiro.checa-garcia@ecmwf.int                                                   #
 #                                                                                         #
 # history:                                                                                #
 #    - Nov-2022 [Ramiro Checa-Garcia]      1st tested version                             #
 #                                                                                         #
 # info:                                                                                   #
 #      This file is likely to be extended with new functions, one per inferface/wrapper   #
 #      to external libraries.                                                             # 
 #                                                                                         #
 #      FUNCTIONS                                                                          #
 #        * readconf        : creates an aerosol object based in configuration file        #
 #        * mie_to_aeropt   : creates an aeropt object from outputs of a mie code          #
 #        * mixing          : calculates an external mixing aerosol given the aeropt of    #
 #                            the components                                               #
 ###########################################################################################


import sys
import numpy as np
import aeropt.aer as aer
from engines.mie_Boucher_Bozzo.libs.mie_BB import mie_boucher_bozzo as mieBB

def interface_mie_Boucher_Bozzo(aerconf, logfile, debug=False, mix=0, verbose=0):
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
    if mix > 3:
        cutoff_radius = 7.5
    else:
        cutoff_radius = np.array(aerconf.bins_max) + 100.0  # at this point I don't understand well this statement

    cutoff_radius = 7.5

    print("                     calculation of optical properties")
    print("                     engine: [Boucher-Bozzo code]")
    print("                     ...calculating optical properties...")

    if debug==False:
        idebug = 0
    else:
        idebug = 2

    #print(logfile)
    logdate = int(logfile.split("/")[1].split(".")[0][0:6])
    logtime = int(logfile.split("/")[1].split(".")[0][6:13])

    #print(logdate)
    #print(logtime)

    # With this two number is possible to reconstruct the filename in Fortran.
    # Rather than add to the default log file we can create
    # tmp/datetime_mieBB.log

    #print('\n\n\n*********************')
    #print(aerconf.ri_lambdatab)
    #print(aerconf.znr_tab)
    #print('\n\n\n*********************')

    sys.stdout.flush()
    out = mieBB.main_mie_aerosols(aerconf.lambda_out, aerconf.bins_min,
                                  aerconf.bins_max, aerconf.sigma_g,
                                  aerconf.r0, aerconf.Ntot, aerconf.rho,
                                  aerconf.rh_tab, aerconf.rh_growth,
                                  cutoff_radius, idebug, aerconf.ri_lambdatab,
                                  aerconf.znr_tab, aerconf.zni_tab,
                                  aerconf.angles )

    aeropt = aer.mie_to_aeropt(aerconf, out, "mie_Boucher_Bozzo")
    
    return aeropt
