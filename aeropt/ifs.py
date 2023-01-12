
 #########################################################################################################
 #                                                                                                       #
 # create_ifs.jl                                                                                         #
 #                                                                                                       #
 # author: Ramiro Checa-Garcia                                                                           #
 # email:  ramiro.checa-garcia@ecmwf.int                                                                 #
 #                                                                                                       #
 # history:                                                                                              #
 #                                                                                                       #
 #     08-Oct-2022     Added the ifs process to create single netcdf with all species                    #
 #                     this is split in several functions to have also consistency test                  #
 #                                                                                                       #
 # info:                                                                                                 #
 #    There are functions to create a single netcdf with a set of aerosol species                        #
 #    - ifs_testdim           => perform few tests of consistency of rules to create ifs file            #
 #    - create_hydro_xxx_dict => create dictionaries to include aerosols netcdfs into ifs new file       #
 #    - create_ifs            => main function to create ifs netcdf file                                 #
 #                                                                                                       #
 #########################################################################################################

import numpy as np
from netCDF4 import Dataset as NCDataset
from netCDF4 import stringtochar as to_char
from netCDF4 import stringtoarr as to_arr
import sys
from datetime import datetime

def ifs_testdim(dic_nciaer, runset):
    """ Evaluates the consistency between a ifs setting and a set of netcdf files used to build it.

    Args:
       list_nciaer (dict): Dictionary species => ncname
       runset      (dict): Dictionary from toml setting file.


    Return:
       dim_rh      (int) : Size dimension rel_hum.
       dim_wl      (int) : Size dimension waveleght
       ifsphobic   (int) : Number of phobic species
       ifsphilic   (int) : Number of philic species
       dic_species (dict): Dict with species

    Auxiliary functions:

       alleq(x)  =True if all elements in x are equal
       phipho(x) = 0   if x="philic" otherwise 1

    """

    def alleq(x):
        """Check if all elements in array/list are equal
        """
        return all(y == x[0] for y in x)

    def phipho(x):
        """0 if x=="philic" 1 otherwise
        """
        out=0 if x=="philic" else 1
        return out


    #=========================================================================
    
    ifs = runset["ifs"]

    # From the dimensions of netcdf files are created the dimensions of nc_ifs
    # -> there are 3 categories, those that should be equal
    # -> those in which ifs should be the maximum between all netcdfs
    # -> those not present in nc_ifs
    #
    dim_names   = ["wavelength", "angle", "size_bin", "rel_hum"]
    dim_equal   = ["wavelength", "angle"]
    dim_maxim   = ["rel_hum"]
    dim_notifs  = ["size_bin"]


    list_nc     = []            # List with ncnames
    dims_nc     = []
    type_nc     = []
    posi_nc     = []
    dim_rh      = 0
    dim_wl      = 0
    species_nc  = dic_nciaer.keys()
    dic_species = {}

    #=========================================================================
    # ifs["aer"] has the species names under a dictionary structure
    # -> here we create two dictionaries:
    #    dic_species: from ifs["aer"]   => species name in species section
    #    rev_sepcies:  from species name => ifs["aer"]
    #

    for ispe in ifs["aer"].keys():
        dic_species[ispe]=ifs["aer"][ispe]["species"]
    
    rev_species = {}
    for k in dic_species.keys():
        rev_species[dic_species[k]]=k


    iaer=0
    for species, ncfile in dic_nciaer.items():
        list_nc.append( NCDataset(ncfile, "r") )
        list_dm = list(list_nc[iaer].dimensions.keys())
        pair_dm = []
        for item in list_dm:
            pair_dm.append( (list_nc[iaer].dimensions[item].name,
                             list_nc[iaer].dimensions[item].size) )
        dims_nc.append( dict(pair_dm) )
        posi_nc.append( ifs["aer"][rev_species[species]]["position"] )
        type_nc.append( phipho(ifs["aer"][rev_species[species]]["position"]) )
        iaer=iaer+1
        

    λ  = []   # Temporal list of all wavelengths dimension size of ncfiles
    ϕ  = []   # Temporal list of all angles      dimension size of ncfiles
    rh = []   # Temporal list of all rel.hum.    dimension size of ncfiles
    for ncdim in dims_nc:
        λ.append(  ncdim["wavelength"])
        ϕ.append(  ncdim["angle"])
        rh.append( ncdim["rel_hum"])

    # Here we perform the tests and inform to the user
    if alleq(λ)==True and alleq(ϕ)==True:
         shared_dim = True
         print("    * [OK] Shared dimensions wavelength and angle consistent \
                           in all species")
    else:
         print("    Shared dimensions wavelength and angle INconsistent in all\
                    species")
         print("    Please, double check the species ncfiles")
         sys.exit()
        
    dim_rh = max(rh)
    if set(rh)==set([dim_rh,1]):
         print("    * [OK] Dimensions relative humidity is reasonable")
    else:
         print("    Aerosols species hydrophilic do not have consistent \
                    values of relative humidity " )
         print("    Revise configuration files of species")
         sys.exit()
        
    pos_rh     = np.argmax(np.array(rh))
    dim_wl     = λ[1]

    ifsphobic  = ifs["number"]["phobic"]
    ifsphilic  = ifs["number"]["philic"]
    phobic     = 0
    philic     = 0
    pos_pho = np.zeros(ifsphobic, dtype = int)
    pos_phi = np.zeros(ifsphilic, dtype = int)
    for namespecies, species in ifs["aer"].items():
        if species["type"]=="phobic":
            phobic = phobic+species["bin"]
            for i in range(species["bin"]):
                pos_pho[species["position"]+i-1]=i+species["position"]-1
            
        elif species["type"]=="philic":
            philic = philic+species["bin"]
            for i in range(species["bin"]):
                pos_phi[species["position"]+i-1]=i+species["position"]-1
        else:
            print("    We only expect type philic or phobic but not :"+species["type"])

    if ifsphobic==phobic:
        print("    * [OK] Number of phobic binned species ",phobic, " consistent with:",ifsphobic)
    else:
        print("    Number of phobic binned species ", phobic," INconsistent with:", ifsphobic)
        print("    Please, double check the species ncfiles")
        sys.exit()

    if ifsphilic==philic:
        print("    * [OK] Number of philic binned species ", philic, " consistent with:" ,ifsphilic)
    else:
        print("    Number of philic binned species ", philic, " INconsistent with:", ifsphilic)
        print("    Please, double check the species ncfiles")
        sys.exit()

    if (pos_pho == np.arange(ifsphobic)).all():
        print("    * [OK] Position of phobic species is consistent")
    else:
        print("    Position of phobic species is INconsistent")
        print("    Please, double check the ifs section in settings file")
        sys.exit()

    if (pos_phi == np.arange(ifsphilic)).all():
        print("    * [OK] Position of philic species is consistent")
    else:
        print("    Position of philic species is INconsistent")
        print("    Please, double check the ifs section in settings file")
        sys.exit()
    
    #TODO: we can check that the sum of all size_bins for philic and phobic should be consistent
    #      with inputs.

    for ii in list_nc:
        ii.close()
    
    return dim_rh, dim_wl, ifsphobic, ifsphilic, dic_species

def create_hydro_cdf_dict(nhydro, hydro, rev_species, dic_nciaer, ifs):
    """Return a dictionary with the hydro species in the right other and the
    corresponding ncfiles opened with python-netcdf4.


    """

    dic_code = {'dd1': 'DD', 'dd2': 'DD', 'dd3': 'DD', 'dd4': 'DD'
               , 'org_dry': 'OM', 'om_dry': 'OM'
               , 'bc1': 'BC', 'bc2': 'BC', 'bc3': 'BC'
               , 'su_dry': 'SU', 'ss': 'SS', 'org': 'OM', 'om': 'OM'
               , 'su': 'SU', 'su1': 'SU', 'su2': 'SU'
               , 'soab': 'OB', 'soaa': 'OA'
               , 'ammf': 'AM', 'nif': 'NI', 'nic': 'NI'}

    d_cdf_hydro={}

    for ii in range(nhydro+1):
        for ikey in ifs["aer"].keys():
            if ifs["aer"][ikey]["position"]==ii and ifs["aer"][ikey]["type"]==hydro:
               if ifs["aer"][ikey]["bin"]==1:
                   d_cdf_hydro[ii-1]=(NCDataset(dic_nciaer[rev_species[ikey]]), 0, dic_code[ikey])
               else:
                   for jj in range(ifs["aer"][ikey]["bin"]): 
                       n = ii+jj-1
                       d_cdf_hydro[n]=(NCDataset(dic_nciaer[rev_species[ikey]]), jj, dic_code[ikey])
                  
    return d_cdf_hydro


def create_hydro_ncname_dict(nhydro, hydro, rev_species, dic_nciaer, ifs):

    dic_code = {'dd1': 'DD', 'dd2': 'DD', 'dd3': 'DD', 'dd4': 'DD'
               , 'org_dry': 'OM', 'om_dry': 'OM'
               , 'bc1': 'BC', 'bc2': 'BC', 'bc3': 'BC'
               , 'su_dry': 'SU', 'ss': 'SS', 'org': 'OM', 'om': 'OM'
               , 'su': 'SU', 'su1': 'SU', 'su2': 'SU'
               , 'soab': 'OB', 'soaa': 'OA'
               , 'ammf': 'AM', 'nif': 'NI', 'nic': 'NI'}

    d_nc_hydro={}
    for ii in range(nhydro):
        for ikey in keys(ifs["aer"]):
            if ifs["aer"][ikey]["position"]==ii and ifs["aer"][ikey]["type"]==hydro:
                if ifs["aer"][ikey]["bin"]==1:
                    d_nc_hydro[ii-1]=(dic_nciaer[rev_species[ikey]],0, dic_code[ikey])
                else:
                    for jj in range(ifs["aer"][ikey]["bin"]):
                       n = ii+jj-1
                       d_nc_hydro[n]=(dic_nciaer[rev_species[ikey]], jj, dic_code[ikey])

    return d_nc_hydro


def process_ifs(dic_nciaer, runset, fsetting, rinfo, ncformat="NETCDF4"):
    # Here we will store in a single netcdf all the nciaer aerosols using dict_settings

    ifs = runset["ifs"]
    dim_rh, dim_wl, ifsphobic, ifsphilic, rev_species = ifs_testdim(dic_nciaer, runset)
    str_today = datetime.today().strftime('%Y-%m-%d')

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

    for ipho in range(ifsphobic):
        print(ipho, ifsphobic)
        if "component" in dd_pho[ipho][0]["size_bin_min"].dimensions:
           min_rad_pho[ipho] = np.mean(dd_pho[ipho][0]["size_bin_min"][dd_pho[ipho][1], :])
           max_rad_pho[ipho] = np.mean(dd_pho[ipho][0]["size_bin_max"][dd_pho[ipho][1], :])
        else:
           print(ipho, min_rad_pho)
           print(dd_pho.keys())
           print(dd_pho[ipho][0]["size_bin_min"][:])
           print(dd_pho[ipho][1])
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
                ri_real_pho[:, ipho]  = np.mean(dd_pho[ipho][0]["ref_idx_real"][:,3,:], axis=1)
                ri_imag_pho[:, ipho]  = np.mean(dd_pho[ipho][0]["ref_idx_img" ][:,3,:], axis=1)
            else:
                ri_real_pho[:, ipho]  = dd_pho[ipho][0]["ref_idx_real"][:,1]
                ri_imag_pho[:, ipho]  = dd_pho[ipho][0]["ref_idx_img" ][:,1]

            asy_pho[:, ipho] = dd_pho[ipho][0]["asymmetry_factor"     ][:, 3, dd_pho[ipho][1]]
            mme_pho[:, ipho] = dd_pho[ipho][0]["extinction"           ][:, 3, dd_pho[ipho][1]]
            ssa_pho[:, ipho] = dd_pho[ipho][0]["single_scatter_albedo"][:, 3, dd_pho[ipho][1]]
            lid_pho[:, ipho] = dd_pho[ipho][0]["lidar_ratio"          ][:, 3, dd_pho[ipho][1]]

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



def process_ifs_49R1(dic_nciaer, runset, fsetting, rinfo, ncformat="NETCDF4"):
    # Here we will store in a single netcdf all the nciaer aerosols using dict_settings

    # dic_code = {'dd1': 'du1', 'dd2': 'du2', 'dd3': 'du3', 'dd4': 'du4', 'org_dry': 'organic', 'om_dry': 'om', 'bc1': 'bc1', 'bc2': 'bc2', 'bc3': 'bc3', 'su_dry': 'su', 'ss': 'ss', 'org': 'organic', 'om': 'om', 'su': 'su', 'su1': 'su1', 'su2': 'su2', 'soab': 'SOA1', 'soaa': 'SOA2', 'ammf': 'NH3_fine', 'nif': 'ni_fine', 'nic': 'ni_coarse'}
    ifs = runset["ifs"]
    dim_rh, dim_wl, ifsphobic, ifsphilic, rev_species = ifs_testdim(dic_nciaer, runset)
    str_today = datetime.today().strftime('%Y-%m-%d')
   
    print(ifsphobic)
    print(ifsphilic)
    print(rev_species)
    
    
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
    #s2phi= ("relative_humidity","hydrophilic",)
    s2phi= ("hydrophilic", "relative_humidity",)
    #s3phi= ("wavenumber","relative_humidity","hydrophilic",)
    s3phi= ("hydrophilic","relative_humidity","wavenumber",)
    s2pho = ("hydrophobic","wavenumber",)

    # Dimension variables  =====================================================

    #code          = ds.createVariable("code_str", 'S2', ("code_str",))

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

    bin_hydrophilic           = ds.createVariable("bin_hydrophilic", 'f4', ("hydrophilic",))
    bin_hydrophilic.long_name = "Hydrophilic aerosol size bin"
    bin_hydrophilic.comment   = "A value of zero indicates that this aerosol type is not partitioned into bins."
    
    code_philic               = ds.createVariable("code_hydrophilic",'S1', ("hydrophilic","code_str",))
    code_philic.longname      = "Hydrophilic aerosol code"
    code_philic.description   = """SS: Sea salt \n OM: Hydrophilic organic matter \n SU: Sulfate  \n OB: Secondary organic biogenic \n OA: Secondary organic anthropogenic \n AM: Fine-mode ammonium sulfate \n NI: Nitrate """
    optmod_philic             = ds.createVariable("optical_model_hydrophilic",'S1', ("hydrophilic","optical_model_str",))
    optmod_philic.long_name   = "Hydrophilic aerosol optical model"
    
    min_rad_phi            = ds.createVariable( "min_radius_hydrophilic", 'f4', ("hydrophilic",))
    min_rad_phi.units      = "m",
    min_rad_phi.long_name  = "Minimum radius of distribution of each hydrophilic aerosol"

    max_rad_phi            = ds.createVariable( "max_radius_hydrophilic", 'f4', ("hydrophilic",))
    max_rad_phi.units      = "m",
    max_rad_phi.long_name  = "Maximum radiations of distribution of each hydrophilic aerosol"

    growth_f_phi           = ds.createVariable( "growth_factor_hydrophilic", 'f4', s2phi)
    growth_f_phi.units     = "m",
    growth_f_phi.long_name = "Relativity Humidity growth factor of hydrophilic aerosols"

    ri_real_phi            = ds.createVariable( "ref_index_real_hydrophilic", 'f4', s3phi)
    ri_real_phi.units      = "1",
    ri_real_phi.long_name  = "Real part of refractive index of each hydrophilic aerosols"
    ri_imag_phi            = ds.createVariable( "ref_index_imag_hydrophilic", 'f4', s3phi)
    ri_imag_phi.units      = "1",
    ri_imag_phi.long_name  = "Imaginary part of refractive index of each hydrophilic aerosols"

    mme_phi                = ds.createVariable( "mass_ext_hydrophilic", 'f4', s3phi)
    mme_phi.units          = "m2 kg-1",
    mme_phi.long_name      = "Mass-extinction coefficient of each hydrophilic aerosols"

    ssa_phi                = ds.createVariable( "ssa_hydrophilic", 'f4', s3phi)
    ssa_phi.units          = "1",
    ssa_phi.long_name      = "Single scattering of each hydrophilic aerosols"

    asy_phi                = ds.createVariable( "asymmetry_hydrophilic", 'f4', s3phi)
    asy_phi.units          = "1",
    asy_phi.long_name      = "Asymmetry factor of each hydrophilic aerosols"

    lid_phi                = ds.createVariable( "lidar_ratio_hydrophilic", 'f4', s3phi)
    lid_phi.units          = "sr",
    lid_phi.long_name      = "Lidar ratio of each hydrophilic aerosols"

    # Hydrophobic Variables ====================================================

    bin_hydrophobic           = ds.createVariable("bin_hydrophobic", 'f4', ("hydrophobic",))
    bin_hydrophobic.long_name = "Hydrophobic aerosol size bin" ;
    bin_hydrophobic.comment   = "A value of zero indicates that this aerosol type is not partitioned into bins." 

    code_phobic               = ds.createVariable("code_hydrophobic",'S1', ("hydrophobic","code_str",))
    code_phobic.longname      = "Hydrophobic aerosol code"
    code_phobic.description   = """DD: Desert dust \n OM: Hydrophobic organic matter  \n BC: Black carbon \n SU: Stratospheric sulfate """

    optmod_phobic             = ds.createVariable("optical_model_hydrophobic",'S1', ("hydrophobic","optical_model_str",))
    optmod_phobic.long_name   = "Hydrophobic aerosol optical model"
    
    min_rad_pho           = ds.createVariable( "min_radius_hydrophobic", 'f4', ("hydrophobic",))
    min_rad_pho.units     = "m",
    min_rad_pho.long_name = "Minimum radiations of distribution of each hydrophobic aerosol"

    max_rad_pho           = ds.createVariable( "max_radius_hydrophobic", 'f4', ("hydrophobic",))
    max_rad_pho.units     = "m",
    max_rad_pho.long_name = "Maximum radiations of distribution of each hydrophobic aerosol"

    ri_real_pho           = ds.createVariable( "ref_index_real_hydrophobic", 'f4', s2pho)
    ri_real_pho.units     = "1",
    ri_real_pho.long_name = "Real part of refractive index of each hydrophobic aerosols"

    ri_imag_pho           = ds.createVariable( "ref_index_imag_hydrophobic", 'f4', s2pho)
    ri_imag_pho.units     = "1",
    ri_imag_pho.long_name = "Imaginary part of refractive index of each hydrophobic aerosols"

    mme_pho               = ds.createVariable( "mass_ext_hydrophobic", 'f4', s2pho)
    mme_pho.units         = "m2 kg-1",
    mme_pho.long_name     = "Mass-extinction coefficient of each hydrophobic aerosols"

    ssa_pho               = ds.createVariable( "ssa_hydrophobic", 'f4',  s2pho)
    ssa_pho.units         = "1",
    ssa_pho.long_name     = "Single scattering of each hydrophobic aerosols"

    asy_pho               = ds.createVariable( "asymmetry_hydrophobic", 'f4',  s2pho)
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
    rh2[:] = np.concatenate((rh1[1::]*0.01, [1.0]), axis=None)
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
        print(dd_pho[ipho][0])
        optmod_phobic[ipho] = to_arr(dd_pho[ipho][0].getncattr("optical_model"),15,'S') 

        code_phobic[ipho]=to_arr(dd_pho[ipho][2],2,'S')
        if "component" in dd_pho[ipho][0]["size_bin_min"].dimensions:
           min_rad_pho[ipho] = np.mean(dd_pho[ipho][0]["size_bin_min"][dd_pho[ipho][1], :])
           max_rad_pho[ipho] = np.mean(dd_pho[ipho][0]["size_bin_max"][dd_pho[ipho][1], :])
        else:
           min_rad_pho[ipho] = dd_pho[ipho][0]["size_bin_min"][dd_pho[ipho][1]]
           max_rad_pho[ipho] = dd_pho[ipho][0]["size_bin_max"][dd_pho[ipho][1]]
    

        if dd_pho[ipho][0]["rel_hum"][:].size==1:
            if "component" in dd_pho[ipho][0]["ref_idx_real"]:
                ri_real_pho[ipho, :]  = np.mean(dd_pho[ipho][0]["ref_idx_real"][:, 0, :], axis=1)
                ri_imag_pho[ipho, :]  = np.mean(dd_pho[ipho][0]["ref_idx_img" ][:, 0, :], axis=1)

            else:
                ri_real_pho[ipho, :]  = dd_pho[ipho][0]["ref_idx_real"][:, 0]
                ri_imag_pho[ipho, :]  = dd_pho[ipho][0]["ref_idx_img" ][:, 0]

            asy_pho[ipho, :] = dd_pho[ipho][0]["asymmetry_factor"     ][:, 0, dd_pho[ipho][1]]
            mme_pho[ipho, :] = dd_pho[ipho][0]["extinction"           ][:, 0, dd_pho[ipho][1]]
            ssa_pho[ipho, :] = dd_pho[ipho][0]["single_scatter_albedo"][:, 0, dd_pho[ipho][1]]
            lid_pho[ipho, :] = dd_pho[ipho][0]["lidar_ratio"          ][:, 0, dd_pho[ipho][1]]
        else:
            if "component" in dd_pho[ipho][0]["ref_idx_real"].dimensions:
                ri_real_pho[ipho, :]  = np.mean(dd_pho[ipho][0]["ref_idx_real"][:,3,:], axis=1)
                ri_imag_pho[ipho, :]  = np.mean(dd_pho[ipho][0]["ref_idx_img" ][:,3,:], axis=1)
            else:
                ri_real_pho[ipho, :]  = dd_pho[ipho][0]["ref_idx_real"][:,1]
                ri_imag_pho[ipho, :]  = dd_pho[ipho][0]["ref_idx_img" ][:,1]

            asy_pho[ipho, :] = dd_pho[ipho][0]["asymmetry_factor"     ][:, 3, dd_pho[ipho][1]]
            mme_pho[ipho, :] = dd_pho[ipho][0]["extinction"           ][:, 3, dd_pho[ipho][1]]
            ssa_pho[ipho, :] = dd_pho[ipho][0]["single_scatter_albedo"][:, 3, dd_pho[ipho][1]]
            lid_pho[ipho, :] = dd_pho[ipho][0]["lidar_ratio"          ][:, 3, dd_pho[ipho][1]]

    asy_phi[:,:,:] = np.zeros( (ifsphilic, dim_rh, dim_wl))
    mme_phi[:,:,:] = np.zeros( (ifsphilic, dim_rh, dim_wl))
    ssa_phi[:,:,:] = np.zeros( (ifsphilic, dim_rh, dim_wl))
    lid_phi[:,:,:] = np.zeros( (ifsphilic, dim_rh, dim_wl))

    for iphi in range(ifsphilic):
       
        optmod_philic[iphi] = to_arr(dd_phi[iphi][0].getncattr("optical_model"),15,'S') 
        code_philic[iphi]=to_arr(dd_phi[iphi][2], 2, 'S')
        if "component" in dd_phi[iphi][0]["size_bin_min"].dimensions:
            min_rad_phi[iphi]  = np.mean(dd_phi[iphi][0]["size_bin_min"][dd_phi[iphi][1], :])
            max_rad_phi[iphi]  = np.mean(dd_phi[iphi][0]["size_bin_max"][dd_phi[iphi][1], :])
            growth_f_phi[iphi, :] = np.mean(dd_phi[iphi][0]["rel_hum_growth"][:,:], axis=1)
        else:
            min_rad_phi[iphi]  = dd_phi[iphi][0]["size_bin_min"][dd_phi[iphi][1]]
            max_rad_phi[iphi]  = dd_phi[iphi][0]["size_bin_max"][dd_phi[iphi][1]]
            growth_f_phi[iphi, :] = dd_phi[iphi][0]["rel_hum_growth"][:]


        if "component" in dd_phi[iphi][0]["ref_idx_real"].dimensions:
            ri_real_phi[iphi, :, :] = np.transpose(np.mean(dd_phi[iphi][0]["ref_idx_real"][:,:,:], axis=2))
            ri_imag_phi[iphi, :, :] = np.transpose(np.mean(dd_phi[iphi][0]["ref_idx_img" ][:,:,:], axis=2))
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
    ds.description_hydrophilic = " \n ".join(description_hydrophilic)

    description_hydrophobic = []
    for ihydro in range(ifsphobic):
        description_hydrophobic.append(metainfo["description_hydrophobic"][str(ihydro+1)])
    ds.description_hydrophobic = " \n ".join(description_hydrophobic)
    ds.close()


    print("\n   IFS aerosol optical file stored at ", outifsnc)

    return



