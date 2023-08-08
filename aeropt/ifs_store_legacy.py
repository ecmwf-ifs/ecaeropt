
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

def process_ifs(dic_nciaer, runset, fsetting, rinfo, ncformat="NETCDF4"):
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

    # Defining dimensions ======================================================
    ds.createDimension("wavenumber"        , dim_wl)
    ds.createDimension("relative_humidity" , dim_rh)
    ds.createDimension("hydrophobic"       , ifsphobic)
    ds.createDimension("hydrophilic"       , ifsphilic)

    # Global attributes   ======================================================
    # at the end

    # Tuples with dimensions to define variables later on ======================

    s1rh = ("relative_humidity",)
    s2phi= ("relative_humidity","hydrophilic",)
    s3phi= ("wavenumber","relative_humidity","hydrophilic",)

    # Dimension variables  =====================================================

    ν             = ds.createVariable("wavenumber", 'f4', ("wavenumber",))
    ν.units       = "cm-1",
    ν.long_name   = "wavenumber"

    rh            = ds.createVariable("relative_humidity", 'f4', s1rh)
    rh.units      = "1"
    rh.long_name  = "Central relative humidity"

    phi           = ds.createVariable("hydrophilic", 'f4', ("hydrophilic",))
    phi.units     = "-"
    phi.long_name = "hydrophilic index"

    pho           = ds.createVariable("hydrophobic", 'f4', ("hydrophobic",))
    pho.units     = "-"
    pho.long_name = "hydrophobic index"

    # Defining variables  ======================================================

    λ             = ds.createVariable("wavelength", 'f4', ("wavenumber",))
    λ.units       = "m"
    λ.long_name   = "wavelength"

    rh1           = ds.createVariable("relative_humidity1", 'f4', s1rh )
    rh1.units     = "1"
    rh1.long_name = "Lower bounds of relative humidity"

    rh2           = ds.createVariable("relative_humidity2", 'f4', s1rh )
    rh2.units     = "1"
    rh2.long_name = "Upper bounds of relative humidity"


    # Hydrophilic Variables ====================================================
    
    min_rad_phi            = ds.createVariable( "min_rad_hydrophilic", 'f4', ("hydrophilic",))
    min_rad_phi.units      = "m",
    min_rad_phi.long_name  = "Minimum radius of distribution of each hydrophilic aerosol"

    max_rad_phi            = ds.createVariable( "max_rad_hydrophilic", 'f4', ("hydrophilic",))
    max_rad_phi.units      = "m",
    max_rad_phi.long_name  = "Maximum radiations of distribution of each hydrophilic aerosol"

    growth_f_phi           = ds.createVariable( "growth_factor_hydrophilic", 'f4', s2phi)
    growth_f_phi.units     = "m",
    growth_f_phi.long_name = "Relativity Humidity growth factor of hydrophilic aerosols"

    ri_real_phi            = ds.createVariable( "ri_real_hydrophilic", 'f4', s3phi)
    ri_real_phi.units      = "1",
    ri_real_phi.long_name  = "Real part of refractive index of each hydrophilic aerosols"
    ri_imag_phi            = ds.createVariable( "ri_imag_hydrophilic", 'f4', s3phi)
    ri_imag_phi.units      = "1",
    ri_imag_phi.long_name  = "Imaginary part of refractive index of each hydrophilic aerosols"

    mme_phi                = ds.createVariable( "mass_ext_hydrophilic", 'f4', s3phi)
    mme_phi.units          = "m2 kg-1",
    mme_phi.long_name      = "Mass-extinction coefficient of each hydrophilic aerosols"

    ssa_phi                = ds.createVariable( "ssa_hydrophilic", 'f4', s3phi)
    ssa_phi.units          = "1",
    ssa_phi.long_name      = "Single scattering of each hydrophilic aerosols"

    asy_phi                = ds.createVariable( "asymetry_hydrophilic", 'f4', s3phi)
    asy_phi.units          = "1",
    asy_phi.long_name      = "Asymmetry factor of each hydrophilic aerosols"

    lid_phi                = ds.createVariable( "lidar_ratio_hydrophilic", 'f4', s3phi)
    lid_phi.units          = "sr",
    lid_phi.long_name      = "Lidar ratio of each hydrophilic aerosols"

    # Hydrophobic Variables ====================================================

    s2pho = ("wavenumber","hydrophobic",)

    min_rad_pho           = ds.createVariable( "min_rad_hydrophobic", 'f4', ("hydrophobic",))
    min_rad_pho.units     = "m",
    min_rad_pho.long_name = "Minimum radiations of distribution of each hydrophobic aerosol"

    max_rad_pho           = ds.createVariable( "max_rad_hydrophobic", 'f4', ("hydrophobic",))
    max_rad_pho.units     = "m",
    max_rad_pho.long_name = "Maximum radiations of distribution of each hydrophobic aerosol"

    ri_real_pho           = ds.createVariable( "ri_real_hydrophobic", 'f4', s2pho)
    ri_real_pho.units     = "1",
    ri_real_pho.long_name = "Real part of refractive index of each hydrophobic aerosols"

    ri_imag_pho           = ds.createVariable( "ri_imag_hydrophobic", 'f4', s2pho)
    ri_imag_pho.units     = "1",
    ri_imag_pho.long_name = "Imaginary part of refractive index of each hydrophobic aerosols"

    mme_pho               = ds.createVariable( "mass_ext_hydrophobic", 'f4', s2pho)
    mme_pho.units         = "m2 kg-1",
    mme_pho.long_name     = "Mass-extinction coefficient of each hydrophobic aerosols"

    ssa_pho               = ds.createVariable( "ssa_hydrophobic", 'f4',  s2pho)
    ssa_pho.units         = "1",
    ssa_pho.long_name     = "Single scattering of each hydrophobic aerosols"

    asy_pho               = ds.createVariable( "asymetry_hydrophobic", 'f4',  s2pho)
    asy_pho.units         = "1",
    asy_pho.long_name     = "Asymmetry factor of each hydrophobic aerosols"

    lid_pho               = ds.createVariable( "lidar_ratio_hydrophobic", 'f4',  s2pho)
    lid_pho.units         = "sr",
    lid_pho.long_name     = "Lidar ratio of each hydrophobic aerosols"


    dd_pho = create_hydro_cdf_dict(ifsphobic, "phobic", rev_species, dic_nciaer, ifs)
    dd_phi = create_hydro_cdf_dict(ifsphilic, "philic", rev_species, dic_nciaer, ifs)
    
    # Setting dimensions values:
    λ[:]   = dd_phi[1][0]["wavelength"][:]
    rh1[:] = dd_phi[1][0]["rel_hum"][:]
    rh2[:] = np.concatenate((rh1[1::], [1.0]), axis=None)
    rh[:]  = (rh1[:]+rh2[:])*0.5
    ν[:]   = 0.01 / dd_phi[1][0]["wavelength"][:]

    species_nc = dic_nciaer.keys()


    asy_pho[:,:]   = np.zeros( (dim_wl, ifsphobic))
    mme_pho[:,:]   = np.zeros( (dim_wl, ifsphobic))
    ssa_pho[:,:]   = np.zeros( (dim_wl, ifsphobic))
    lid_pho[:,:]   = np.zeros( (dim_wl, ifsphobic))
    min_rad_phi[:] = np.zeros( ifsphilic)
    max_rad_phi[:] = np.zeros( ifsphilic)
    min_rad_pho[:] = np.zeros( ifsphobic)
    max_rad_pho[:] = np.zeros( ifsphobic)

    rhphobic=ifs["rhphobic"]

    for ipho in range(ifsphobic):
        if "component" in dd_pho[ipho][0]["size_bin_min"].dimensions:
           min_rad_pho[ipho] = np.mean(dd_pho[ipho][0]["size_bin_min"][dd_pho[ipho][1], :])
           max_rad_pho[ipho] = np.mean(dd_pho[ipho][0]["size_bin_max"][dd_pho[ipho][1], :])
        else:
           min_rad_pho[ipho] = dd_pho[ipho][0]["size_bin_min"][dd_pho[ipho][1]]
           max_rad_pho[ipho] = dd_pho[ipho][0]["size_bin_max"][dd_pho[ipho][1]]
    

        if dd_pho[ipho][0]["rel_hum"][:].size==1:
            if "component" in dd_pho[ipho][0]["ref_idx_real"]:
                ri_real_pho[:, ipho]  = np.mean(dd_pho[ipho][0]["ref_idx_real"][:, 0, :], axis=1)
                ri_imag_pho[:, ipho]  = np.mean(dd_pho[ipho][0]["ref_idx_img" ][:, 0, :], axis=1)

            else:
                ri_real_pho[:, ipho]  = dd_pho[ipho][0]["ref_idx_real"][:, 0]
                ri_imag_pho[:, ipho]  = dd_pho[ipho][0]["ref_idx_img" ][:, 0]

            asy_pho[:, ipho] = dd_pho[ipho][0]["asymmetry_factor"     ][:, 0, dd_pho[ipho][1]]
            mme_pho[:, ipho] = dd_pho[ipho][0]["extinction"           ][:, 0, dd_pho[ipho][1]]
            ssa_pho[:, ipho] = dd_pho[ipho][0]["single_scatter_albedo"][:, 0, dd_pho[ipho][1]]
            lid_pho[:, ipho] = dd_pho[ipho][0]["lidar_ratio"          ][:, 0, dd_pho[ipho][1]]
        else:

            if "component" in dd_pho[ipho][0]["ref_idx_real"].dimensions:
                ri_real_pho[:, ipho]  = np.mean(dd_pho[ipho][0]["ref_idx_real"][:,rhphobic,:], axis=1)
                ri_imag_pho[:, ipho]  = np.mean(dd_pho[ipho][0]["ref_idx_img" ][:,rhphobic,:], axis=1)
            else:
                ri_real_pho[:, ipho]  = dd_pho[ipho][0]["ref_idx_real"][:,1]
                ri_imag_pho[:, ipho]  = dd_pho[ipho][0]["ref_idx_img" ][:,1]

            asy_pho[:, ipho] = dd_pho[ipho][0]["asymmetry_factor"     ][:, rhphobic, dd_pho[ipho][1]]
            mme_pho[:, ipho] = dd_pho[ipho][0]["extinction"           ][:, rhphobic, dd_pho[ipho][1]]
            ssa_pho[:, ipho] = dd_pho[ipho][0]["single_scatter_albedo"][:, rhphobic, dd_pho[ipho][1]]
            lid_pho[:, ipho] = dd_pho[ipho][0]["lidar_ratio"          ][:, rhphobic, dd_pho[ipho][1]]

    asy_phi[:,:,:] = np.zeros( (dim_wl, dim_rh, ifsphilic))
    mme_phi[:,:,:] = np.zeros( (dim_wl, dim_rh, ifsphilic))
    ssa_phi[:,:,:] = np.zeros( (dim_wl, dim_rh, ifsphilic))
    lid_phi[:,:,:] = np.zeros( (dim_wl, dim_rh, ifsphilic))

    for iphi in range(ifsphilic):

        if "component" in dd_phi[iphi][0]["size_bin_min"].dimensions:
            min_rad_phi[iphi]  = np.mean(dd_phi[iphi][0]["size_bin_min"][dd_phi[iphi][1], :])
            max_rad_phi[iphi]  = np.mean(dd_phi[iphi][0]["size_bin_max"][dd_phi[iphi][1], :])
            growth_f_phi[:, iphi] = np.mean(dd_phi[iphi][0]["rel_hum_growth"][:,:], axis=1)
        else:
            min_rad_phi[iphi]  = dd_phi[iphi][0]["size_bin_min"][dd_phi[iphi][1]]
            max_rad_phi[iphi]  = dd_phi[iphi][0]["size_bin_max"][dd_phi[iphi][1]]
            growth_f_phi[:, iphi] = dd_phi[iphi][0]["rel_hum_growth"][:]


        if "component" in dd_phi[iphi][0]["ref_idx_real"].dimensions:
            ri_real_phi[:, :, iphi] = np.mean(dd_phi[iphi][0]["ref_idx_real"][:,:,:], axis=2)
            ri_imag_phi[:, :, iphi] = np.mean(dd_phi[iphi][0]["ref_idx_img" ][:,:,:], axis=2)
        else:
            ri_real_phi[:, :, iphi] = dd_phi[iphi][0]["ref_idx_real"][:,:]
            ri_imag_phi[:, :, iphi] = dd_phi[iphi][0]["ref_idx_img" ][:,:]

        asy_phi[:, :, iphi] = dd_phi[iphi][0]["asymmetry_factor"     ][:, :, dd_phi[iphi][1]]
        mme_phi[:, :, iphi] = dd_phi[iphi][0]["extinction"           ][:, :, dd_phi[iphi][1]]
        ssa_phi[:, :, iphi] = dd_phi[iphi][0]["single_scatter_albedo"][:, :, dd_phi[iphi][1]]
        lid_phi[:, :, iphi] = dd_phi[iphi][0]["lidar_ratio"          ][:, :, dd_phi[iphi][1]]


    metainfo                   = runset["info"]
    ds.history                 = "Created "+str_today+". With ecaerrad v"+rinfo.version+ " and setting file at path "+fsetting+". Run by "+rinfo.user+"."
    ds.contact                 = metainfo["contact"]
    ds.institution             = metainfo["institution"]
    ds.comment                 = metainfo["comment_string"]+metainfo["product_version"]
    description_hydrophilic = []
    for ihydro in range(ifsphilic):
        description_hydrophilic.append(metainfo["description_hydrophilic"][str(ihydro+1)])
    ds.description_hydrophilic = " ; ".join(description_hydrophilic)

    description_hydrophobic = []
    for ihydro in range(ifsphobic):
        description_hydrophobic.append(metainfo["description_hydrophobic"][str(ihydro+1)])
    ds.description_hydrophobic = " ; ".join(description_hydrophobic) 

    ds.close()


    print("\n   IFS aerosol optical file stored at ", outifsnc)

    return

