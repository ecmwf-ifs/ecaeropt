#  +----------------------------------------------------------------------------------------+
#  | aeropt/modes.py                                                                        |
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
#  |    26-Nov-2022   Ramiro Checa-Garcia    Added check_vars_nc                            |
#  |                                                                                        |
#  |                                                                                        |
#  | Info:                                                                                  |
#  |      Provides the FUNCTIONS:                                                           |
#  |                                                                                        |
#  |      FUNCTIONS                                                                         |
#  |         * config_file_mode   :                                                         |
#  |         * info_mode          :                                                         |
#  |         * setting_file_mode  :                                                         |
#  |         * test_mode          :                                                         |
#  |         * compare_ifsnc      :                                                         |
#  |         * check_vars_nc      :                                                         |
#  |                                                                                        |
#  +----------------------------------------------------------------------------------------+


try:
    import tomllib as toml
except ModuleNotFoundError:
    import toml


import os.path
import subprocess 
import numpy as np
from netCDF4 import Dataset as NCDataset

import aeropt.show       as show
import aeropt.process    as run
import aeropt.parsetoml  as parse
import aeropt.ifs        as ifs
import aeropt.plt        as plt

def config_file_mode(rinfo, config_file, outname):
    """
    Perform a calculation with the information of a configuration file
    (that provides the information of an aerosol configuration).

    Args:
        rinfo       (object) : run info of calculation
        config_file (str)    : configuration file path
        outname     (str)    : filename for output netcdf

    Return:
        None

    """


    show.add_header(rinfo, config_file)
    print("\n   RUN in CONFIG FILE MODE: \n")
    angles = np.linspace(180.0,0.0,180)
    run.single(rinfo, config_file, outname, angles )
    show.add_footer()

    return


def info_mode(rinfo):
    """
    Show info about the current engines installed and
    perform basic tests about their consistency.

    Args:
        rinfo       (object) : run info of calculation

    Return:
        None

    """

    show.add_header(rinfo, "info")
    show.add_info()
    show.add_footer()

    return


def setting_file_mode(rinfo, fsetting, test=False):
   """
   Perfom a set of calculations described in a setting file

   Args:
       rinfo    (object): run info object 
       fsetting (str)   : path to setting file
       test     (bool)  : flag to know if called for testing.

   Return:
       dicnc_iaer (dic) : dictionary with reference and calculated netcdf for each aerosol
       runset           : extra-info

   """

   cwarn = (show.bcolor.WARN +show.bcolor.BOLD,show.bcolor.ENDC)
   show.add_header(rinfo, fsetting)

   if os.path.isfile(fsetting):
       runset = toml.load(fsetting)     # This should be load only once
   else:
       runset = None

   if test==False:
       print("\n   RUN in SETTING FILE MODE\n")
   else:

       print("\n   RUN in "+cwarn[0]+"TEST MODE"+cwarn[1]+"\n")
   
   
   dicnc_iaer = {}
   if os.path.isfile(fsetting):
        ifs_flag = "ifs" in runset.keys()
        if "plt" in runset["process"].keys():
            plt_flag = runset["process"]["plt"]
        else:
            plt_flag = False
        if "wavelengths" in runset["default"].keys():
            wl_out_default = runset["default"]["wavelengths"]
        else:
            wl_out_default = "none"

        if parse.check_toml(runset):
            # Typically a process with have several steps: 
            # mixing aerosols, single aerosols,...
            steps = runset["process"] 
            
            if "mixing" in steps.keys():
                aermix = parse.toml_mixing(runset)
                if steps["skip"]==True:
                    nc_mix = aermix["ncname"]
                else:
                    nc_mix = run.mixture( rinfo
                                        , aermix["lconf"] , aermix["ncname"]
                                        , aermix["nangles"], aermix["lengine"]
                                        , wl_out=wl_out_default, opt_model=aermix["opt_model"])
                if test==False:
                    dicnc_iaer[aermix["mixture"]]=nc_mix
                else:
                    nc_ref = runset["process"]["mixing"][aermix["mixture"]]["ref"] 
 
                    dicnc_iaer[aermix["mixture"]]={"nc_test":nc_mix 
                                                  ,"nc_ref" :nc_ref}

            if "single" in steps.keys():
                dic_aersingle, angles = parse.toml_single(runset)

                n_single = len(dic_aersingle.keys())
                print("\n   Processing ",n_single," single aerosols ..........\n")
                for iaer, aer in enumerate(dic_aersingle.keys()):
                    print("         * Single aerosol ", iaer+1, "/", n_single, "........")
                    if steps["skip"]==True:
                       nc_iaer = dic_aersingle[aer]["ncname"]
                    else:
                       # here we can add process to a thread pool
                       nc_iaer = run.single(rinfo,
                                            dic_aersingle[aer]["config"], 
                                            dic_aersingle[aer]["ncname"], 
                                            angles, 
                                            aerengine=dic_aersingle[aer]["engine"],
                                            wl_out=wl_out_default)
                    if test==False:
                        dicnc_iaer[aer]=nc_iaer
                    else:
                        #print(runset["species"][aer]["ref"])
                        dicnc_iaer[aer]={ "nc_test":nc_iaer, 
                                          "nc_ref": runset["species"][aer]["ref"]}

            if plt_flag:

                print("\n Plotting refractive index and optical properties ..........\n")
                for specie in dicnc_iaer.keys():
                    print("\n      Aerosol Specie netcdf: ",specie)
                    plt.aerplt("refindex", dicnc_iaer[specie])
                    plt.aerplt("optical",  dicnc_iaer[specie])
                    plt.aerplt("phasefunction",  dicnc_iaer[specie])

            if ifs_flag:
                ifs_vers = runset["ifs"]["ifs_cycle"]
                print("\n Creating a IFS netcdf file with all aerosols ..........\n")
                if ifs_vers < "49R1":
                   ifs.process_ifs(dicnc_iaer, runset, fsetting, rinfo)
                else:
                   ifs.process_ifs_49R1(dicnc_iaer, runset, fsetting, rinfo)
            elif test==False:
                print(show.bcolor.WARN+show.bcolor.BOLD+"   WARNING: "+show.bcolor.ENDC)
                print(" IFS section not in setting file. Code ends here. \n")
 
   else:
        print(show.bcolor.FAIL+show.bcolor.BOLD+"   ERROR: "+show.bcolor.ENDC)
        print(" Setting file ",fsetting,"not found in filesystem. \n")

   if test==False:

       if "nc_ref" in runset["ifs"].keys():
           if runset["ifs"]["nc_ref"]!=None:
              compare_ifsnc(runset, "IFS final netcdf")
              show.add_footer()
           else: 
              show.add_footer()
       else:
           show.add_footer()
   else:
        print("\n   Calculations with ecaerrad done. Proceeding with comparison. \n")
        #show.add_footer()
   return dicnc_iaer, runset
                  


def test_mode(rinfo, fsetting):
    """ 
    Perfom a calculation of a setting file prepared for a test mode
          Note: after calling the setting mode, it creates temporal netcdf
          with the difference between calculation and reference netcdf
          using nco cli-tool

    Args:
        rinfo  (object): run info object  
        fsetting  (str): path to settting file

    Return:
        None

    """


    list_nc, dict_fset = setting_file_mode(rinfo, fsetting, test=True)
    for iaer, aer in list_nc.items():
           pathnctest = aer["nc_test"]
           pathncref  = aer["nc_ref"]
           nctest     = os.path.basename(pathnctest)
           ncref      = os.path.basename(pathncref)
           diffcmd    = "aeropt/compare_secure.sh "+pathnctest+" "+pathncref 
           subprocess.Popen("aeropt/compare_secure.sh "+pathnctest+" "+pathncref, shell=True).wait()
           ncdiff = "tmp/"+nctest.replace(".nc","")+"_minus_"+ncref
           check_vars_nc(ncdiff, iaer, pathnctest, pathncref)

    runset = toml.load(fsetting)
    if "ifs" in runset.keys():

       if "nc_ref" in runset["ifs"].keys():
          if runset["ifs"]["nc_ref"] != None:
             compare_ifsnc(runset,"IFS-netcdf")

    show.add_footer()
    
    return


def compare_ifsnc(runset, iaer):
    pathnctest = runset["ifs"]["netcdfname"]
    pathncref  = runset["ifs"]["nc_ref"]
    nctest     = os.path.basename(pathnctest)
    ncref      = os.path.basename(pathncref)
    subprocess.Popen("aeropt/compare_secure.sh "+pathnctest+" "+pathncref, shell=True).wait()
    ncdiff = "tmp/"+nctest.replace(".nc","")+"_minus_"+ncref
    check_vars_nc(ncdiff, iaer, pathnctest, pathncref)
    return

def check_vars_nc(ncdiff, iaer, pathnctest, pathncref):
    """
    Perform comparison between variables in two netcdf files.

    Args: 
        ncdiff     (str): path to nc difference between calc and ref.
        iaer       (int): just an integer of the aerosol compared
        pathnctest (str): path to the calculated netcdf
        pathncref  (str): path to the reference netcdf
    
    Return:
        passing   (bool): indicates if test has been passed or not
                          a test is passed if all variables are the 
                          same within a tolerance value.
    
    """


    cpass = (show.bcolor.GREEN+show.bcolor.BOLD,show.bcolor.ENDC)
    cwarn = (show.bcolor.WARN +show.bcolor.BOLD,show.bcolor.ENDC)
    cfail = (show.bcolor.FAIL +show.bcolor.BOLD,show.bcolor.ENDC)
    cbold = (show.bcolor.BOLD, show.bcolor.ENDC) 
    ds    = NCDataset(ncdiff, "r")
    dsref = NCDataset(pathncref,"r") 
    print("\n      Checking difference for ",iaer, " between: \n")
    print("            netcdf to test   : ", pathnctest)
    print("            netcdf reference : ", pathncref,"\n")
    print("")
    str_vars = ["phase_function","code_hydrophilic", "code_hydrophobic", "optical_model_hydrophilic", "optical_model_hydrophobic"]
    passing=True
    for varname in ds.variables:
        if varname not in str_vars:
           if varname in ds.dimensions:
              mydim=varname
           else:
              vmin = np.amin(ds[varname][:])
              vmax = np.amax(ds[varname][:])
              rmin = np.amin(ds[varname][:]/dsref[varname][:])
              rmax = np.amax(ds[varname][:]/dsref[varname][:])
              print(cbold[0]+"       "+varname+cbold[1]+"\n")
              print(         "             minimum difference: ", vmin, rmin)
              print(         "             maximum difference: ", vmax, rmax)
              if abs(vmin)<=1e-10 and abs(vmax)<= 1e-10:
                 print(cpass[0]+"             Test passed absolute diff."+cpass[1]+" (threshold 1.e-10) \n")
              elif abs(rmin)<=1.e-5 and abs(rmax)<= 1.e-5:
                 print(cpass[0]+"             Test passed relative diff."+cpass[1]+" (threshold 1.e-5)  \n")
              else:
                 if varname in ["rel_hum_growth", "ref_idx_real", "ref_idx_img"]:
                    print(cwarn[0]+"             Test unclear "+cwarn[1]+"\n")
                    print(    "                    ... for externally mixed aerosols this variable ", varname)
                    print(    "                    ... because new version stores all components not only 1st")
                 else:
                    print(cfail[0]+"             Test not passed "+cfail[1]+"  ... keep testing...\n")
                    ashape=ds[varname].shape
                    for idim in range(ashape[0]):
                            print("           specie number: ",idim+1, "\n")
                            if len(ashape)>1:
                               pvmin = np.amin(ds[varname][idim,:])
                               pvmax = np.amax(ds[varname][idim,:])
                               prmin = np.amin(ds[varname][idim,:]/dsref[varname][idim,:])
                               prmax = np.amax(ds[varname][idim,:]/dsref[varname][idim,:])

                               print(         "             minimum difference: ", pvmin, prmin)
                               print(         "             maximum difference: ", pvmax, prmax)
                    passing=False
               
    return passing

