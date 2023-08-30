#  +----------------------------------------------------------------------------------------+
#  | aeropt/parsetoml.py                                                                    |
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
#  |    06-Dec-2022   Ramiro Checa-Garcia    Added check_toml                               |
#  |                                                                                        |
#  |                                                                                        |
#  | Info:                                                                                  |
#  |      Provides CLASSES and FUNCTIONS:                                                   |
#  |                                                                                        |
#  |                                                                                        |
#  |      FUNCTIONS                                                                         |
#  |        * nangle         :                                                              |
#  |        * ncname         :                                                              |
#  |        * select_ncname  :                                                              |
#  |        * check_toml     :                                                              |
#  |        * toml_single    :                                                              |
#  |        * each_mixture   :                                                              |
#  |        * toml_mixing    :                                                              |
#  |                                                                                        |
#  +----------------------------------------------------------------------------------------+


import os
import os.path
import numpy as np

try:
    import tomllib as toml
except ModuleNotFoundError:
    import toml

from datetime import datetime

def nangle(nangle):
     """
     Creates a numpy array of angles based on input

     Args:
         nangle (int) : create a linear spaced array of nangles between 0 and 180.
         nangle (str) : read the file with path nangle with the given angles.

     Return:
         angles (numpy array) with angles for phase function calculation.

     """

     if isinstance(nangle, int) and nangle==0:
        angles = np.array([])
     
     if isinstance(nangle, int) and nangle>0:
        angles = np.linspace(0.0, 180.0, nangle) 
     
     if isinstance(nangle, list) and len(nangle)==3:
         angles= np.linspace(nangle[0],nangle[1],nangle[2])
     else:
         print("\n    ERROR: ****** Inconsistency in nangle specification in input file *******")
         exit()

     if isinstance(nangle, str): 
        angles = np.loadtext(nangle, skiprows=1)
        # The format of the nangle file is one description line
        # and n rows each with one angle in degrees.
        #return 
     
     return angles
 

def ncname(aer_conf, ncname,  outncdir,  mix=False, skip=False):
    """
    Ascertain the ncname given an list with the rules.

    Args:
        aer_conf
        ncname
        outncdir
        mix
        skip

    Return:
       outncname (str): path with outncname 

    """
    for sname in ncname:
        if sname=="C:today":
            ncname = [a for a in ncname if a != sname] 
            ncname.append("C:today")
        
    
    if mix==False:
        aer_type = toml.load(aer_conf)["tags"]["aer_type"]
    else:
        aer_type = str(aer_conf)
    
    today = datetime.today().strftime('%Y-%m-%d')
    strncname = "_".join(ncname)
    strncname = strncname.replace( "C:aer_type" , aer_type )
    strncname = strncname.replace( "C:today"    , today )
    if strncname[-3:-1] !=".nc":
        strncname = strncname+".nc"
    

    ncname = os.path.join(outncdir, strncname)
    if os.path.isfile(ncname) and skip==False:
        print("********** File with name: ", ncname )
        print("           is already on your disk.")
        print("           Please change your settings or options.")
        print("           ...  it is finished now the processing.")
        print("================================================================================\n")
        exit()
    
    return os.path.join(outncdir, strncname)


def select_ncname(iaer, aerconfig, defaultnc, outncdir, mix=False, skip=False):
    """
    Create the output ncname based on inputs and rules

    Args:


    Return:
       outncname (string) : filename for output netcdf
    """

    if "ncname" in iaer.keys():
         if iaer["ncname"]=="default":
             outncname = ncname(aerconfig, defaultnc,      outncdir, mix=mix, skip=skip)
         else:
             outncname = ncname(aerconfig, iaer["ncname"], outncdir, mix=mix, skip=skip)
         
    else:
         outncname = ncname(aerconfig, defaultnc , outncdir, mix=mix, skip=skip )
    
    return outncname



def check_toml(runset, goal="settings"):
    """
    Perform a basic check of setting file

    Args:


    Return:


    """
    tomlkeys = runset.keys()
    if goal == "settings":
        if "main" in tomlkeys and len(tomlkeys)<=2:
           print("********** settings with name: ", ftoml )
           print("           At command line you aim to use a setting toml file")
           print("           but it seems that probably you used a config file ")
           print("           rather than a setting file")
           print("           ...  it is finished now the processing.")
           print("================================================================================\n")
           exit()
       
        if set(["process", "species", "default"]).issubset(set(tomlkeys)):
           if "info" in tomlkeys:
               return True
           else:
               print("   [WARNING]   It's recommented to add info section in settings file")
               return True
           
        else:
           print("********** settings with name: ", ftoml )
           print("           At command line you aim to use a setting toml file")
           print("           but it might be incomplete: we recomm sections: ")
           print("           default => with your default settings")
           print("           species => with your list of available config files.")
           print("           process => with your request of processing.")
           print("           info    => (optional) with additional metadata for netcdf.")
           print("           ...  it is finished now the processing.")
           print("================================================================================\n")
           exit()
        
    return


def toml_single(runset):
    """
    
    
    """

    skip = runset["process"]["skip"]
    defaultnc= runset["default"]["ncname"]
    outncdir = runset["default"]["outncdir"]
    nangles  = nangle(runset["default"]["nangle"])

    dict_aer = {}
    for aerconf in runset["process"]["single"]:
        species_aer = runset["species"][aerconf]
        if "config" in species_aer.keys():
            dict_aer[aerconf]={}
            dict_aer[aerconf]["config"]= species_aer["config"]

            if "engine" in species_aer.keys() and species_aer["engine"]!="default":
                 dict_aer[aerconf]["engine"]=species_aer["engine"]
            else:
                 dict_aer[aerconf]["engine"]=runset["default"]["engine"]
             
            dict_aer[aerconf]["ncname"] = select_ncname( species_aer,
                                                          dict_aer[aerconf]["config"],
                                                          defaultnc,
                                                          outncdir, mix=False, skip=skip )
        else:
            print(" It is mandatory to add a config file to the species section ")
            print(" ... Now the code will finish prematurely")
            exit()
         
    return dict_aer, nangles
    

 


def each_mixture(species_mix, l_aer_config, l_aer_engine, default_engine):

    for component in species_mix.keys():
        mix_component = species_mix[component]
        if "config" in mix_component.keys():
            l_aer_config.append( mix_component["config"])

            if "engine" in mix_component.keys() and mix_component["engine"]!="default":
                l_aer_iengine.append(species_mix["engine"])
            else:
                l_aer_engine.append(default_engine)
            
        else:
            print(" It is mandatory to add a config file to the species section ")
            print(" ... Now the code will finish prematurely")
            exit()
        
    return l_aer_config, l_aer_engine



def toml_mixing(runset):

    skip         = runset["process"]["skip"]
    defaultnc    = runset["default"]["ncname"]
    outncdir     = runset["default"]["outncdir"]
    nangles      = nangle(runset["default"]["nangle"])

    l_aer_config = []
    l_aer_engine = []
    outncname    = ""
    aerconf_mix  = []
    l_opt_model  = []

    dic_run_mixing = runset["process"]["mixing"]
    dic_species    = runset["species"]

    for mixture in dic_run_mixing.keys():
        species_mix = dic_species[mixture]
        l_aer_config, l_aer_engine = each_mixture(species_mix,
                                                  l_aer_config,
                                                  l_aer_engine,
                                                  runset["default"]["engine"])
        outncname = select_ncname( dic_run_mixing[mixture],
                                   mixture, defaultnc, 
                                   outncdir,  mix=True, 
                                   skip=skip )
        aerconf_mix.append(mixture)
        if "opt_model" in dic_run_mixing[mixture].keys():
            l_opt_model.append(dic_run_mixing[mixture]["opt_model"])
        else:
            l_opt_model.append(None)

    # Here there is a kind of inconsistency as this is outside of the loop.
    # this code will not work with more than one mixing case.

    # Here we also want to check that opt_model is present before this


    aermix = { "lconf":l_aer_config, 
               "lengine":l_aer_engine,
               "ncname":outncname,
               "nangles":nangles, 
               "mixture":aerconf_mix[0],
               "opt_model": l_opt_model[0] }

    return aermix # l_aer_config, l_aer_engine, outncname, nangles, aerconf_mix[0]
    


