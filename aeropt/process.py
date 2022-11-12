


 #########################################################################################################
 #                                                                                                       #
 # process.jl                                                                                            #
 #                                                                                                       #
 # author: Ramiro Checa-Garcia                                                                           #
 # email:  ramiro.checa-garcia@ecmwf.int                                                                 #
 #                                                                                                       #
 # history:                                                                                              #
 #                                                                                                       #
 #     26-Oct-2022     Translated from Julia to Python                                                   #
 #                                                                                                       #
 # info:                                                                                                 #
 #    There are two functions that process a given aerosol either single or mixed                        #
 #                                                                                                       #
 #                                                                                                       #
 #########################################################################################################


import sys
import os
import os.path
import logging

import aeropt.engine   as engine
import aeropt.store_nc as store
import aeropt.aer      as aer


import numpy as np
from netCDF4 import Dataset as NCDataset
from datetime import datetime

def single(rinfo, path_conf_file, ncoutname, angles, aerengine="mie_Boucher_Bozzo", wl_out="none"):          

    if os.path.isfile(ncoutname):
        print("********** File with name: ", ncoutname )
        print("           is already on your disk.")
        print("           Please change your settings or options.")
        print("           ...  it is finished now the processing.")
        print("================================================================================\n")
        sys.exit()
    
    aer_conf = aer.readconf(path_conf_file, angles, wl_out=wl_out)

    print("                     config: ", path_conf_file)
    print("                     nϕ, nλ: ", aer_conf.nmumax, ", ", aer_conf.nb_lambda)
    print("                     rh, nD: ", aer_conf.rh_int, ", ", aer_conf.size_bins)
    
    if rinfo.debug==True:
       logging.debug("\n== Debug information block ================\n")
       logging.debug("\n\n Showing aerosol configuration loaded from: "+path_conf_file)
       logging.debug("\n\n")
       logging.debug(aer_conf)
    
    if aerengine=="mie_Boucher_Bozzo":
        aer_opt = engine.interface_mie_Boucher_Bozzo(aer_conf, rinfo.logfile, debug=rinfo.debug)
    else:
        print("---- ERROR ---- ", aerengine , " not yet implemented")
        sys.exit()
   
    print("                     =>   optical prop. calculated")
    print("               Storing results as netcdf with path:")
    print("               ",ncoutname," \n")

    store.store_nc_single(aer_conf, aer_opt, rinfo, ncname=ncoutname)

    return ncoutname


def mixture(rinfo, path_conf_files, ncoutname, nangle, laerengine, wl_out="none"):


    print("\n   Processing mixing aerosols.........\n")

    # Reading first component to fix array sizes

    mix_aer_obj = []
    mix_aer_opt = []
    mix_ri_ltab = []
    mix_ri_rtab = []
    mix_ri_itab = []

    # We have to refactor this, here in process there should be a loop that
    # and store all the aer_opt for each component, then the mixing is calculated
    # in a subroutine. It will use more memory but it is more structured.

    mix_component=0
    num_component=len(path_conf_files)

    for component_path, component_engine in zip(path_conf_files, laerengine):

        mix_component=mix_component+1

        print("         * Aerosol component ",mix_component, "/", num_component,":")

        aer_conf = aer.readconf(component_path, nangle, wl_out=wl_out)

        print("                     config: ", component_path)
        print("                     nϕ, nλ: ", aer_conf.nmumax, ", ", aer_conf.nb_lambda)
        print("                     rh, nD: ", aer_conf.rh_int, ", ", aer_conf.size_bins)

        if rinfo.debug==True:
            logging.debug("\n== Debug information block ================\n")
            logging.debug("\n\n Showing aerosol configuration loaded from: "+path_conf_file)
            logging.debug("\n\n")
            logging.debug(aer_conf)

        if component_engine=="mie_Boucher_Bozzo":
            aer_opt = engine.interface_mie_Boucher_Bozzo(aer_conf, rinfo.logfile)
        else:
            print("---- ERROR ---- ", component_engine , " not yet implemented")
            sys.exit()

        print("                     => component calculated")

        mix_aer_obj.append(aer_conf)
        mix_aer_opt.append(aer_opt)
        mix_ri_ltab.append(aer_conf.ri_lambdatab)
        mix_ri_rtab.append(aer_conf.znr_tab)
        mix_ri_itab.append(aer_conf.zni_tab)

    
    print("\n         All components calculated: estimating mixture ")

    mix_aer, aer_mix_opt = aer.mixing(mix_aer_obj, mix_aer_opt,
                                      mix_ri_ltab, mix_ri_rtab,
                                      mix_ri_itab, num_component)

    print("\n         Storing results as netcdf with path:")
    print("         ", ncoutname )

    store.store_nc_mixture(mix_aer, aer_mix_opt, rinfo, ncname=ncoutname)

    #print("\n         Total Cpu Time used: ")

    return ncoutname





