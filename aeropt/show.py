
#  +----------------------------------------------------------------------------------------+
#  | aeropt/show.py                                                                         |
#  |                                                                                        |
#  | (C) Copyright 2022- ECMWF.                                                             |
#  |                                                                                        |
#  | This software is licensed under the terms of the Apache Licence Version 2.0            |
#  | which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.                   |
#  |                                                                                        |
#  | In applying this licence, ECMWF does not waive the privileges and immunities           |
#  | granted to it by virtue of its status as an intergovernmental organisation             |
#  | nor does it submit to any jurisdiction.                                                |
#  |                                                                                        |
#  |                                                                                        |
#  | Author:                                                                                |
#  |    Ramiro Checa-Garcia. ECMWF                                                          |
#  |                                                                                        |
#  | Modifications:                                                                         |
#  |    30-Oct-2022   Ramiro Checa-Garcia    Added quick-tests to mie-info                  |
#  |                                                                                        |
#  |                                                                                        |
#  | Info:                                                                                  |
#  |      Note: function run_header depends on global variables.                            |
#  |                                                                                        |
#  |      CLASSES                                                                           |
#  |        * bcolor                                                                        |
#  |                                                                                        |
#  |      FUNCTIONS                                                                         |
#  |        * add_header   :                                                                |
#  |        * add_footer   :                                                                |
#  |        * add_info     :                                                                |
#  |                                                                                        |
#  +----------------------------------------------------------------------------------------+


import aeropt.engine as engine 

def add_footer():

    print("\n Run completed. Have a nice day ")
    print("=====================================================================================\n")

    return

def add_info():

    mie="""
     - Mie-Spherical Distributions:
     Engine based in Fortran 90 code. It assumed log-normal distribution. It estimates
     the aerosol opt. per bin based on an integral with 999 intervals between bin-min
     and bin-max (and function is the log-normal distribution). The number of maximum
     terms for the sum of mie scattering is set as: ...
     """
    print("\n Engines implemented:\n", mie)
    print("     Quick tests of Mie Code")
    
    engine.test_engine("Mie")


    return


def add_header(rinfo, fsetting): # this function depends on global variables



    print("\n=====================================================================================")
    print("       ______ _____             ______ _____         ____  _____ _______   ")
    print("      |  ____/ ____|      /\   |  ____|  __ \       / __ \|  __ \__   __|  ")
    print("      | |__ | |   ______ /  \  | |__  | |__) |_____| |  | | |__) | | |     ")
    print("      |  __|| |  |______/ /\ \ |  __| |  _  /______| |  | |  ___/  | |     ")
    print("      | |___| |____    / ____ \| |____| | \ \      | |__| | |      | |     ")
    print("      |______\_____|  /_/    \_\______|_|  \_\      \____/|_|      |_|     ") 
    print("\n")
    print("   ec offline aerosol optics        : ", rinfo.version+"  (in code)")
    print("   architecture                     : ", rinfo.machine)
    print("   current path                     : ", rinfo.runpath)
    print("   by the user                      : ", rinfo.user   )
    print("   date                             : ", rinfo.date, "\n" )
    if rinfo.debug==True:
        print("   --> run debug mode, logfile      : ", rinfo.logfile) 
    print(" * Processing ...                   : ", fsetting     )

    return


class bcolor:
    HEADER= '\033[95m'
    BLUE  = '\033[94m'
    CYAN  = '\033[96m'
    GREEN = '\033[92m'
    WARN  = '\033[93m'
    FAIL  = '\033[91m'
    ENDC  = '\033[0m'
    BOLD  = '\033[1m'
    UNDER = '\033[4m'

