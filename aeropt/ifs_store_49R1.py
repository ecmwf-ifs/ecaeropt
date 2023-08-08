
#############################################################################################################
# aeropt/ifs.py
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
#    10-Nov-2022   R. Checa-Garcia   1st. tested version
#    12-Nov-2022   R. Checa-Garcia   Added the ifs process to create single netcdf with all species
#                                    this is split in several functions to have also consistency test   
#    25-Jan-2023   T. Stockdale      Checked to introduce volcanic aerosols
#    28-Jan-2023   R. Checa-Garcia   Created function with 49R1 ifs netcdf final format
#
#                                                                                         
# Info: 
#    There are functions to create a single netcdf with a set of aerosol species   
#
#      FUNCTIONS                                                                        
#        * ifs_testdim           : perform few tests of consistency of rules to create ifs file            
#        * create_hydro_xxx_dict : create dictionaries to include aerosols netcdfs into ifs new file       
#        * process_ifs           : main function to create ifs netcdf file                                                                                                    
#        * process_ifs_49R1      : main function to create ifs netcdf file with 49R1 format           
###########################################################################################################

from pprint import pprint
import numpy as np
from netCDF4 import Dataset as NCDataset
from netCDF4 import stringtochar as to_char
from netCDF4 import stringtoarr as to_arr
import sys
from datetime import datetime


def process_ifs_49R1(dic_nciaer, runset, fsetting, rinfo, ncformat="NETCDF3_CLASSIC"):
    # Here we will store in a single netcdf all the nciaer aerosols using dict_settings

    ifs = runset["ifs"]
    dim_rh, dim_wl, ifsphobic, ifsphilic, rev_species = ifs_testdim(dic_nciaer, runset)
    str_today = datetime.today().strftime('%Y-%m-%d')
  

    if ifs["netcdfname"]=="standard":
        outifsnc="outputnc/store_ifsnc/aerosol_ifs_CY"+ifs["ifs_cycle"]+"_"+str_today.replace("-","")+".nc"
    else:
        outifsnc=ifs["netcdfname"]

    print("\n === Storing all the processed aerosols in a single netcdf \n")

    # Opening NC dataset  ======================================================
    ds = NCDataset(outifsnc,"w", format=ncformat)

    dim_code     = 2
    dim_optmodel = 15

    # Defining dimensions ======================================================
    ds.createDimension("wavenumber"        , dim_wl)
    ds.createDimension("relative_humidity" , dim_rh)
    ds.createDimension("hydrophobic"       , ifsphobic)
    ds.createDimension("hydrophilic"       , ifsphilic)
    ds.createDimension("code_str"          , dim_code)
    ds.createDimension("optical_model_str" , dim_optmodel)

    # Global attributes   ======================================================
    # at the end

    # Tuples with dimensions to define variables later on ======================

    s1rh = ("relative_humidity",)
    s2phi= ("hydrophilic", "relative_humidity",)
    s3phi= ("hydrophilic","relative_humidity","wavenumber",)
    s2pho =("hydrophobic","wavenumber",)

    # Dimension variables  =====================================================

    ν             = ds.createVariable("wavenumber", 'f4', ("wavenumber",))
    ν.units       = "cm-1",
    ν.long_name   = "wavenumber"

    rh            = ds.createVariable("relative_humidity", 'f4', s1rh)
    rh.units      = "1"
    rh.long_name  = "Central relative humidity"

    #phi           = ds.createVariable("hydrophilic", 'f4', ("hydrophilic",))
    #phi.units     = "-"
    #phi.long_name = "hydrophilic index"

    #pho           = ds.createVariable("hydrophobic", 'f4', ("hydrophobic",))
    #pho.units     = "-"
    #pho.long_name = "hydrophobic index"

    # Defining variables  ======================================================

    λ             = ds.createVariable("wavelength", 'f4', ("wavenumber",))
    λ.units       = "m"
    λ.long_name   = "wavelength"

    rh1           = ds.createVariable("relative_humidity1", 'f4', s1rh )
    rh1.units     = "1"
    rh1.long_name = "Lower bounds of relative humidity bins"

    rh2           = ds.createVariable("relative_humidity2", 'f4', s1rh )
    rh2.units     = "1"
    rh2.long_name = "Upper bounds of relative humidity bins"


    # Hydrophilic Variables ====================================================

    bin_hydrophilic           = ds.createVariable("bin_hydrophilic", 'i4', ("hydrophilic",))
    bin_hydrophilic.long_name = "Hydrophilic aerosol size bin"
    bin_hydrophilic.comment   = "A value of zero indicates that this aerosol type is not partitioned into bins."
    
    code_philic               = ds.createVariable("code_hydrophilic",'S1', ("hydrophilic","code_str",))
    code_philic.longname      = "Hydrophilic aerosol code"
    code_philic.description   = """SS: Sea salt \nOM: Hydrophilic organic matter \nSU: Sulfate  \nOB: Secondary organic biogenic \nOA: Secondary organic anthropogenic \nAM: Fine-mode ammonium sulfate \nNI: Nitrate """
    optmod_philic             = ds.createVariable("optical_model_hydrophilic",'S1', ("hydrophilic","optical_model_str",))
    optmod_philic.long_name   = "Hydrophilic aerosol optical model"
    
    min_rad_phi            = ds.createVariable( "min_radius_hydrophilic", 'f4', ("hydrophilic",))
    min_rad_phi.units      = "m",
    min_rad_phi.long_name  = "Minimum radius of distribution of hydrophilic aerosols"

    max_rad_phi            = ds.createVariable( "max_radius_hydrophilic", 'f4', ("hydrophilic",))
    max_rad_phi.units      = "m",
    max_rad_phi.long_name  = "Maximum radius of distribution of hydrophilic aerosols"

    growth_f_phi           = ds.createVariable( "growth_factor_hydrophilic", 'f4', s2phi)
    growth_f_phi.units     = "m",
    growth_f_phi.long_name = "Relativity humidity growth factor of hydrophilic aerosols"

    ri_real_phi            = ds.createVariable( "ref_index_real_hydrophilic", 'f4', s3phi)
    ri_real_phi.units      = "1",
    ri_real_phi.long_name  = "Real part of refractive index of hydrophilic aerosols"
    ri_imag_phi            = ds.createVariable( "ref_index_imag_hydrophilic", 'f4', s3phi)
    ri_imag_phi.units      = "1",
    ri_imag_phi.long_name  = "Imaginary part of refractive index of hydrophilic aerosols"

    mme_phi                = ds.createVariable( "mass_ext_hydrophilic", 'f4', s3phi)
    mme_phi.units          = "m2 kg-1",
    mme_phi.long_name      = "Mass-extinction coefficient of hydrophilic aerosols"

    ssa_phi                = ds.createVariable( "ssa_hydrophilic", 'f4', s3phi)
    ssa_phi.units          = "1",
    ssa_phi.long_name      = "Single scattering albedo of hydrophilic aerosols"

    asy_phi                = ds.createVariable( "asymmetry_hydrophilic", 'f4', s3phi)
    asy_phi.units          = "1",
    asy_phi.long_name      = "Asymmetry factor of hydrophilic aerosols"

    lid_phi                = ds.createVariable( "lidar_ratio_hydrophilic", 'f4', s3phi)
    lid_phi.units          = "sr",
    lid_phi.long_name      = "Lidar ratio of hydrophilic aerosols"

    # Hydrophobic Variables ====================================================

    bin_hydrophobic           = ds.createVariable("bin_hydrophobic", 'i4', ("hydrophobic",))
    bin_hydrophobic.long_name = "Hydrophobic aerosol size bin" ;
    bin_hydrophobic.comment   = "A value of zero indicates that this aerosol type is not partitioned into bins." 

    code_phobic               = ds.createVariable("code_hydrophobic",'S1', ("hydrophobic","code_str",))
    code_phobic.longname      = "Hydrophobic aerosol code"
    code_phobic.description   = """DD: Desert dust \nOM: Hydrophobic organic matter  \nBC: Black carbon \nSU: Stratospheric sulfate """

    optmod_phobic             = ds.createVariable("optical_model_hydrophobic",'S1', ("hydrophobic","optical_model_str",))
    optmod_phobic.long_name   = "Hydrophobic aerosol optical model"
    
    min_rad_pho           = ds.createVariable( "min_radius_hydrophobic", 'f4', ("hydrophobic",))
    min_rad_pho.units     = "m",
    min_rad_pho.long_name = "Minimum radius of distribution of hydrophobic aerosols"

    max_rad_pho           = ds.createVariable( "max_radius_hydrophobic", 'f4', ("hydrophobic",))
    max_rad_pho.units     = "m",
    max_rad_pho.long_name = "Maximum radius of distribution of hydrophobic aerosols"

    ri_real_pho           = ds.createVariable( "ref_index_real_hydrophobic", 'f4', s2pho)
    ri_real_pho.units     = "1",
    ri_real_pho.long_name = "Real part of refractive index of hydrophobic aerosols"

    ri_imag_pho           = ds.createVariable( "ref_index_imag_hydrophobic", 'f4', s2pho)
    ri_imag_pho.units     = "1",
    ri_imag_pho.long_name = "Imaginary part of refractive index of hydrophobic aerosols"

    mme_pho               = ds.createVariable( "mass_ext_hydrophobic", 'f4', s2pho)
    mme_pho.units         = "m2 kg-1",
    mme_pho.long_name     = "Mass-extinction coefficient of hydrophobic aerosols"

    ssa_pho               = ds.createVariable( "ssa_hydrophobic", 'f4',  s2pho)
    ssa_pho.units         = "1",
    ssa_pho.long_name     = "Single scattering albedo of hydrophobic aerosols"

    asy_pho               = ds.createVariable( "asymmetry_hydrophobic", 'f4',  s2pho)
    asy_pho.units         = "1",
    asy_pho.long_name     = "Asymmetry factor of hydrophobic aerosols"

    lid_pho               = ds.createVariable( "lidar_ratio_hydrophobic", 'f4',  s2pho)
    lid_pho.units         = "sr",
    lid_pho.long_name     = "Lidar ratio of hydrophobic aerosols"


    dd_pho = create_hydro_cdf_dict(ifsphobic, "phobic", rev_species, dic_nciaer, ifs)
    dd_phi = create_hydro_cdf_dict(ifsphilic, "philic", rev_species, dic_nciaer, ifs)
    
    # Setting dimensions values:
    λ[:]   = dd_phi[1][0]["wavelength"][:]
    rh1[:] = dd_phi[1][0]["rel_hum"][:]*0.01
    rh2[:] = np.concatenate((rh1[1::], [1.0]), axis=None)
    rh[:]  = (rh1[:]+rh2[:])*0.5*0.01
    ν[:]   = 0.01 / dd_phi[1][0]["wavelength"][:]

    species_nc = dic_nciaer.keys()

    asy_pho[:,:]   = np.zeros( (ifsphobic, dim_wl))
    mme_pho[:,:]   = np.zeros( (ifsphobic, dim_wl))
    ssa_pho[:,:]   = np.zeros( (ifsphobic, dim_wl))
    lid_pho[:,:]   = np.zeros( (ifsphobic, dim_wl))
    min_rad_phi[:] = np.zeros( ifsphilic)
    max_rad_phi[:] = np.zeros( ifsphilic)
    min_rad_pho[:] = np.zeros( ifsphobic)
    max_rad_pho[:] = np.zeros( ifsphobic)



    for ipho in range(ifsphobic):

        if dd_pho[ipho][4] > 0:
            rhphobic=dd_pho[ipho][4]
        else:
            rhphobic=0
        #print(ipho, "=> ", rhphobic)
        #print(dd_pho[ipho][0]["ref_idx_img" ][:].shape)
        #print(dd_pho[ipho][0]["ref_idx_img" ].dimensions) 
        optmod_phobic[ipho] = to_arr(dd_pho[ipho][0].getncattr("optical_model"),15,'S')

        if dd_pho[ipho][3]=="nobinned":
           bin_hydrophobic[ipho]=dd_pho[ipho][5]
        else:
           bin_hydrophobic[ipho]=dd_pho[ipho][1]+1

        code_phobic[ipho]=to_arr(dd_pho[ipho][2],2,'S')
        if "component" in dd_pho[ipho][0]["size_bin_min"].dimensions:
           #min_rad_pho[ipho] = np.mean(dd_pho[ipho][0]["size_bin_min"][dd_pho[ipho][1], :])*1.e-6
           #max_rad_pho[ipho] = np.mean(dd_pho[ipho][0]["size_bin_max"][dd_pho[ipho][1], :])*1.e-6 
           min_rad_pho[ipho] = dd_pho[ipho][0]["size_bin_min"][dd_pho[ipho][1], -1]*1.e-6
           max_rad_pho[ipho] = dd_pho[ipho][0]["size_bin_max"][dd_pho[ipho][1], -1]*1.e-6
       
        else:
           min_rad_pho[ipho] = dd_pho[ipho][0]["size_bin_min"][dd_pho[ipho][1]]*1.e-6
           max_rad_pho[ipho] = dd_pho[ipho][0]["size_bin_max"][dd_pho[ipho][1]]*1.e-6
    

        if dd_pho[ipho][0]["rel_hum"][:].size==1:
            if "component" in dd_pho[ipho][0]["ref_idx_real"]:
                #ri_real_pho[ipho, :]  = np.mean(dd_pho[ipho][0]["ref_idx_real"][:, 0, :], axis=1)
                #ri_imag_pho[ipho, :]  = np.mean(dd_pho[ipho][0]["ref_idx_img" ][:, 0, :], axis=1)
                ri_real_pho[ipho, :]  = dd_pho[ipho][0]["ref_idx_real"][:, 0, -1]
                ri_imag_pho[ipho, :]  = dd_pho[ipho][0]["ref_idx_img" ][:, 0, -1]
            else:
                ri_real_pho[ipho, :]  = dd_pho[ipho][0]["ref_idx_real"][:, 0]
                ri_imag_pho[ipho, :]  = dd_pho[ipho][0]["ref_idx_img" ][:, 0]

            asy_pho[ipho, :] = dd_pho[ipho][0]["asymmetry_factor"     ][:, 0, dd_pho[ipho][1]]
            mme_pho[ipho, :] = dd_pho[ipho][0]["extinction"           ][:, 0, dd_pho[ipho][1]]
            ssa_pho[ipho, :] = dd_pho[ipho][0]["single_scatter_albedo"][:, 0, dd_pho[ipho][1]]
            lid_pho[ipho, :] = dd_pho[ipho][0]["lidar_ratio"          ][:, 0, dd_pho[ipho][1]]
        else:
            if "component" in dd_pho[ipho][0]["ref_idx_real"].dimensions:
                # Mixing.... mean in last dimension which is components!!
                #ri_real_pho[ipho, :]  = np.mean(dd_pho[ipho][0]["ref_idx_real"][:,rhphobic,:], axis=1)
                #ri_imag_pho[ipho, :]  = np.mean(dd_pho[ipho][0]["ref_idx_img" ][:,rhphobic,:], axis=1)
                ri_real_pho[ipho, :]  = dd_pho[ipho][0]["ref_idx_real"][:,rhphobic,-1]
                ri_imag_pho[ipho, :]  = dd_pho[ipho][0]["ref_idx_img" ][:,rhphobic,-1]
            else:
                # rh is last dimension should be rhphobic
                ri_real_pho[ipho, :]  = dd_pho[ipho][0]["ref_idx_real"][:,rhphobic]
                ri_imag_pho[ipho, :]  = dd_pho[ipho][0]["ref_idx_img" ][:,rhphobic]

            asy_pho[ipho, :] = dd_pho[ipho][0]["asymmetry_factor"     ][:, rhphobic, dd_pho[ipho][1]]
            mme_pho[ipho, :] = dd_pho[ipho][0]["extinction"           ][:, rhphobic, dd_pho[ipho][1]]
            ssa_pho[ipho, :] = dd_pho[ipho][0]["single_scatter_albedo"][:, rhphobic, dd_pho[ipho][1]]
            lid_pho[ipho, :] = dd_pho[ipho][0]["lidar_ratio"          ][:, rhphobic, dd_pho[ipho][1]]

    asy_phi[:,:,:] = np.zeros( (ifsphilic, dim_rh, dim_wl))
    mme_phi[:,:,:] = np.zeros( (ifsphilic, dim_rh, dim_wl))
    ssa_phi[:,:,:] = np.zeros( (ifsphilic, dim_rh, dim_wl))
    lid_phi[:,:,:] = np.zeros( (ifsphilic, dim_rh, dim_wl))
    for iphi in range(ifsphilic):
       
        if dd_phi[iphi][3]=="nobinned":
           bin_hydrophilic[iphi]=dd_phi[iphi][5]
        else:
           bin_hydrophilic[iphi]=dd_phi[iphi][1]+1

        optmod_philic[iphi] = to_arr(dd_phi[iphi][0].getncattr("optical_model"),15,'S') 
        code_philic[iphi]=to_arr(dd_phi[iphi][2], 2, 'S')
        if "component" in dd_phi[iphi][0]["size_bin_min"].dimensions:
            #min_rad_phi[iphi]  = np.mean(dd_phi[iphi][0]["size_bin_min"][dd_phi[iphi][1], :])*1.e-6
            #max_rad_phi[iphi]  = np.mean(dd_phi[iphi][0]["size_bin_max"][dd_phi[iphi][1], :])*1.e-6
            #growth_f_phi[iphi, :] = np.mean(dd_phi[iphi][0]["rel_hum_growth"][:,:], axis=1)

            min_rad_phi[iphi]  = dd_phi[iphi][0]["size_bin_min"][dd_phi[iphi][1], -1]*1.e-6
            max_rad_phi[iphi]  = dd_phi[iphi][0]["size_bin_max"][dd_phi[iphi][1], -1]*1.e-6
            growth_f_phi[iphi, :] = dd_phi[iphi][0]["rel_hum_growth"][:,-1]

        else:
            min_rad_phi[iphi]  = dd_phi[iphi][0]["size_bin_min"][dd_phi[iphi][1]]*1.e-6
            max_rad_phi[iphi]  = dd_phi[iphi][0]["size_bin_max"][dd_phi[iphi][1]]*1.e-6
            growth_f_phi[iphi, :] = dd_phi[iphi][0]["rel_hum_growth"][:]


        if "component" in dd_phi[iphi][0]["ref_idx_real"].dimensions:
            #ri_real_phi[iphi, :, :] = np.transpose(np.mean(dd_phi[iphi][0]["ref_idx_real"][:,:,:], axis=2))
            #ri_imag_phi[iphi, :, :] = np.transpose(np.mean(dd_phi[iphi][0]["ref_idx_img" ][:,:,:], axis=2))
            ri_real_phi[iphi, :, :] = np.transpose(dd_phi[iphi][0]["ref_idx_real"][:,:,-1])
            ri_imag_phi[iphi, :, :] = np.transpose(dd_phi[iphi][0]["ref_idx_img" ][:,:,-1])
        else:
            ri_real_phi[iphi, :, :] = np.transpose(dd_phi[iphi][0]["ref_idx_real"][:,:])
            ri_imag_phi[iphi, :, :] = np.transpose(dd_phi[iphi][0]["ref_idx_img" ][:,:])

        asy_phi[iphi, :, :] = np.transpose(dd_phi[iphi][0]["asymmetry_factor"     ][:, :, dd_phi[iphi][1]])
        mme_phi[iphi, :, :] = np.transpose(dd_phi[iphi][0]["extinction"           ][:, :, dd_phi[iphi][1]])
        ssa_phi[iphi, :, :] = np.transpose(dd_phi[iphi][0]["single_scatter_albedo"][:, :, dd_phi[iphi][1]])
        lid_phi[iphi, :, :] = np.transpose(dd_phi[iphi][0]["lidar_ratio"          ][:, :, dd_phi[iphi][1]])


    metainfo                   = runset["info"]
    ds.history                 = "Created "+str_today+". With ecaerrad v"+rinfo.version+ " and setting file at path "+fsetting+". Run by "+rinfo.user+"."
    ds.contact                 = metainfo["contact"]
    ds.institution             = metainfo["institution"]
    ds.comment                 = metainfo["comment_string"]+metainfo["product_version"]
    description_hydrophilic = []
    for ihydro in range(ifsphilic):
        description_hydrophilic.append(metainfo["description_hydrophilic"][str(ihydro+1)])
    ds.description_hydrophilic = "\n ".join(description_hydrophilic)

    description_hydrophobic = []
    for ihydro in range(ifsphobic):
        description_hydrophobic.append(metainfo["description_hydrophobic"][str(ihydro+1)])
    ds.description_hydrophobic = "\n ".join(description_hydrophobic)
    ds.close()

    print("\n   IFS aerosol optical file stored at ", outifsnc)

    return


