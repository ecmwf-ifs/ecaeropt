#  +----------------------------------------------------------------------------------------+
#  | aeropt/store_nc.py                                                                     |
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
#  |    30-Oct-2022   Ramiro Checa-Garcia    Added documentation to functions               |
#  |                                                                                        |
#  |                                                                                        |
#  | Info:                                                                                  |
#  |                                                                                        |
#  |      FUNCTIONS                                                                         |
#  |        * store_nc_single   :                                                           |
#  |        * store_nc_mixture  :                                                           |
#  |                                                                                        |
#  +----------------------------------------------------------------------------------------+


from netCDF4 import Dataset as NCDataset
from datetime import datetime

try:
    import tomllib as toml
except ModuleNotFoundError:
    import toml

def store_nc_single(aer, aeropt, rinfo, ncname="output_test.nc", ncformat="NETCDF4"):
    """
    Creates a netcdf with the optical properties of a single aerosol specie based on
    aer and aeropt objects (usually derived from config files)

    Properties of the netcdf
    Dimensions:
       wavelength -> nb_lambda
       size_bin   -> size_bins
       rel_hum    -> ri_nhr
       angle      -> mux
       leg_coeff  -> nmom

     Variables
       ref_ind_real
       ref_ind_img
       lidar_radio
       asymmetry_factor
       single_scatter_albedo
       extinction
       phase_function
       (legendre_coeff)
       wavelength
       rel_hum
       rel_hum_growth
       size_bin_min
       size_bim_max
       angle
    """

    str_today = rinfo.date
    
    # Opening NC dataset  =================================================
    ds = NCDataset(ncname,"w", format=ncformat)

    # Defining dimensions =================================================

    ds.createDimension("wavelength", aer.nb_lambda)
    ds.createDimension("size_bin"  , aer.size_bins)
    ds.createDimension("rel_hum"   , aer.ri_nrh)
    ds.createDimension("angle"     , aeropt.nmux)

    # Global attributes   =================================================

    opticalmodel=toml.load(aer.config_file, _dict=dict)["tags"]["opt_model"]

    str_config = "Aerosol Config.  file: "+aer.config_file
    str_refind = "Refractive index file: "+aer.ri_file

    str_radii="Mode radius: ("+",".join([str(a) for a in aer.r0])+") microns "
    str_sigma="; geometric standard deviation: (" + ",".join([str(a) for a in aer.sigma_g]) +")"
    str_densi="; density "+str(aer.rho)+"kg/m3 "
    str_ntot ="; total number concentration: ("+",".join([str(a) for a in aer.Ntot])+")  m-3 "


    ds.description               = "Single aerosol netcdf file created with ecaeropt code from ECMWF"
    ds.history                   = "Created " + str_today+" with the ecaeropt tool v1.0"

    ds.source_config             = str_config
    ds.source_refidx             = str_refind
    ds.source_size_distribution  = str_radii+str_sigma+str_densi+str_ntot
    ds.optical_model             = opticalmodel

    ds.institution   = "European Centre for Medium-Range Weather Forecasts (ECMWF)."
    ds.contact       = "ramiro.checa-garcia@ecmwf.int"

    s2 = ("wavelength","rel_hum",)
    s3 = ("wavelength","rel_hum","size_bin",)
    s4 = ("wavelength","rel_hum","size_bin","angle",)

    # Defining variables  ====================================== ref idx ==

    ref_idx_real             = ds.createVariable("ref_idx_real", 'f4', s2)
    ref_idx_real[:]          = aeropt.ri_r
    ref_idx_real.units       = "-"
    ref_idx_real.long_name   = "real part of the refractive index"

    ref_idx_img              = ds.createVariable("ref_idx_img", 'f4',  s2)
    ref_idx_img[:]           = aeropt.ri_i
    ref_idx_img.units        = "-"
    ref_idx_img.long_name    = "imaginary part of the refractive index"

    # Defining variables  ====================================== Lidar ====

    lidar_ratio              = ds.createVariable("lidar_ratio", 'f4', s3)
    lidar_ratio[:]           = aeropt.lidar
    lidar_ratio.units        = "sr"
    lidar_ratio.long_name    = "lidar ratio"


    # Defining variables  ================================== Bulk Optical ==

    asymmetry_fac            = ds.createVariable("asymmetry_factor", 'f4', s3)
    asymmetry_fac.units      = "-"
    asymmetry_fac.long_name  = "asymmetry factor"
    asymmetry_fac[:]         = aeropt.asy

    single_sca_alb           = ds.createVariable("single_scatter_albedo", 'f4', s3)
    single_sca_alb.units     = "-"
    single_sca_alb.long_name = "single scatter albedo"
    single_sca_alb[:]        = aeropt.ssa

    extinction               = ds.createVariable("extinction", 'f4', s3)
    extinction.units         = "m^2/kg"
    extinction.long_name     = "mass extinction coefficient"
    extinction[:]            = aeropt.ext * aer.ext_scaling

    # Defining variables  =============================== Phase Function ==

    phase_func               = ds.createVariable("phase_function", 'f4', s4)
    phase_func.units         = "-"
    phase_func.long_name     = "scattering phase function"
    phase_func.description   = ("normalized: pfun(cos(ang))*d(cos(ang)) "
                              + "between 0 and pi is =2 an the integral "
                              + "on the full solid angle is 4*PI")
    phase_func[:]            = aeropt.pfun

    # Defining variables  ============================== Humidity growth ==

    rel_hum_growth           = ds.createVariable("rel_hum_growth", 'f4', ("rel_hum",))
    rel_hum_growth.units     = "-"
    rel_hum_growth.long_name = "relative humidity growth factor"
    rel_hum_growth[:]        = aer.rh_growth


    # Defining variables  =============================== Mueller P11    ==

    mueller_p11              = ds.createVariable("p11", 'f4', s4)
    mueller_p11.units         = "-"
    mueller_p11.long_name     = "Mueller Matrix p11"
    mueller_p11.description   = "Mueller Matrix p11. Still we need to check normalization"
    mueller_p11[:]            = aeropt.p11

    # Defining variables  =============================== Mueller P11    ==

    mueller_p12              = ds.createVariable("p12", 'f4', s4)
    mueller_p12.units         = "-"
    mueller_p12.long_name     = "Mueller Matrix p12"
    mueller_p12.description   = "Mueller Matrix p12. Still we need to check normalization"
    mueller_p12[:]            = aeropt.p12

    # Defining variables  =============================== Mueller P11    ==

    mueller_p33              = ds.createVariable("p33", 'f4', s4)
    mueller_p33.units         = "-"
    mueller_p33.long_name     = "Mueller Matrix p33"
    mueller_p33.description   = "Mueller Matrix p33. Still we need to check normalization"
    mueller_p33[:]            = aeropt.p33

    # Defining variables  =============================== Mueller P11    ==

    mueller_p34              = ds.createVariable("p34", 'f4', s4)
    mueller_p34.units         = "-"
    mueller_p34.long_name     = "Mueller Matrix p34"
    mueller_p34.description   = "Mueller Matrix p34. Still we need to check normalization"
    mueller_p34[:]            = aeropt.p34


    # Defining variables  =================================== Dimensions ==

    wavelength               = ds.createVariable("wavelength", 'f4', ("wavelength",))
    wavelength.units         = "m"
    wavelength.long_name     = "monochromatic wavelength"
    wavelength[:]            = aer.lambda_out

    rel_hum                  = ds.createVariable("rel_hum", 'f4', ("rel_hum",))
    rel_hum.units            = "%"
    rel_hum.long_name        = "ambient relative humidity"
    rel_hum[:]               = aer.rh_tab

    size_bin_min             = ds.createVariable("size_bin_min", 'f4', ("size_bin",))
    size_bin_min.units       = "m"
    size_bin_min.long_name   = "lower limit of the distribution size bin"
    size_bin_min[:]          = aer.bins_min

    size_bin_max             = ds.createVariable("size_bin_max", 'f4', ("size_bin",))
    size_bin_max.units       = "m"
    size_bin_max.long_name   = "upper limit of the distribution size bin"
    size_bin_max[:]          = aer.bins_max

    angle                    = ds.createVariable("angle", 'f4', ("angle",))
    angle.unit               = "degrees"
    angle.long_name          = "phase function angle"
    angle[:]                 = aeropt.pfun_ang

    # Closing netcdf ======================================================
    ds.close()




def store_nc_mixture(laer_conf, aeropt, rinfo, ncname="output_test.nc", ncformat="NETCDF4", opt_model=None):
    """
    Creates a netcdf with the optical properties of a mixture aerosol species based on
    aer and aeropt objects (usually derived from config files)

    Properties of the netcdf
    Dimensions:
       wavelength -> nb_lambda
       size_bin   -> size_bins
       rel_hum    -> ri_nhr
       angle      -> mux
       (leg_coeff  -> nmom)

     Variables
       ref_ind_real
       ref_ind_img
       lidar_radio
       asymmetry_factor
       single_scatter_albedo
       extinction
       phase_function
       (legendre_coeff)
       wavelength
       rel_hum
       rel_hum_growth
       size_bin_min
       size_bim_max
       angle
    """
    def _alleq(lst):
        leq = [False]*len(lst)
        for i,x in enumerate(lst):
            if x == lst[0]:
                leq[i]=True
        if False in leq:
            return False
        else:
            return True

        return
    aer = laer_conf[0]


    num_components = len(laer_conf)
    str_today = rinfo.date
    # Opening NC dataset  =================================================
    ds = NCDataset(ncname,"w", format=ncformat)

    # Defining dimensions =================================================
    ds.createDimension("wavelength", aer.nb_lambda)
    ds.createDimension("size_bin"  , aer.size_bins)
    ds.createDimension("rel_hum"   , aer.ri_nrh)
    ds.createDimension("angle"     , aeropt.nmux)
    ds.createDimension("component" , num_components)

    # Define a global attribute
    ds.description = ("Aerosol Optical Properties of a mixed aerosol with "
                     + str(num_components)+" components")
    ds.history     = "Created "+str_today+" with the ecaeropt tool v1.0"

    if opt_model == None:
        lopticalmodel = [None]*len(laer_conf)
        for iaer, aerconf in enumerate(laer_conf):    
            lopticalmodel[iaer]=toml.load(aerconf.config_file, _dict=dict)["tags"]["opt_model"]

        if _alleq(lopticalmodel):
            ds.optical_model = lopticalmodel[0]
        else:
            print(" ==> INCONSISTENCY in OPTICAL MODEL OF AEROSOL-EXTERNAL MIXTURE")
            print(" ==> Revise configuration files: ")
            for ii, optmodel in zip([a.config_file for a in laer_conf], lopticalmodel):
                print(" -> File: ",ii, optmodel)
            exit()
    else: 
        ds.optical_model = opt_model


    str_config       = "Configuration files   : "+",".join([aer.config_file for aer in laer_conf])
    str_refidx       = "Refractive index files: "+",".join([aer.ri_file for aer in laer_conf])

    ds.configuration = str_config
    ds.refr_index    = str_refidx
    ds.institution   = "European Centre for Medium-Range Weather Forecasts (ECMWF)."
    ds.contact       = "ramiro.checa-garcia@ecmwf.int"
    #str_radii="Mode radius: ("*join([string(a) for a in aer.r0],",")*") microns "
    #str_sigma="; geometric standard deviation: ("*join([string(a) for a in aer.σ_g],",")*")"
    #str_densi="; density "*string(aer.ρ)*"kg/m3 "
    #str_ntot ="; total number concentration: (" *join([string(a) for a in aer.Ntot],",")*")  m-3 "

    #ds.parameters"]  = str_radii*str_sigma*str_densi*str_ntot

    # Define the variablesg

    s3idx = ("wavelength", "rel_hum",   "component", )
    s3bin = ("wavelength", "rel_hum",   "size_bin",  )
    s4bin = ("wavelength", "rel_hum",   "size_bin",  "angle", )
    s2com = ("size_bin",   "component", )
    s2hum = ("rel_hum", "component",)

    ref_idx_real             = ds.createVariable("ref_idx_real", 'f4', s3idx)
    ref_idx_real[:]          = aeropt.ri_r
    ref_idx_real.units       = "-"
    ref_idx_real.long_name   = "real part of the refractive index"

    ref_idx_img              = ds.createVariable("ref_idx_img", 'f4', s3idx)
    ref_idx_img[:]           = aeropt.ri_i
    ref_idx_img.units        = "-"
    ref_idx_img.long_name    = "imaginary part of the refractive index"

    lidar_ratio              = ds.createVariable( "lidar_ratio", 'f4', s3bin)
    lidar_ratio[:]           = aeropt.lidar
    lidar_ratio.units        = "sr"
    lidar_ratio.long_name    = "lidar ratio"

    asymmetry_fact           = ds.createVariable("asymmetry_factor", 'f4', s3bin)
    asymmetry_fact.units     = "-"
    asymmetry_fact.long_name = "asymmetry factor"
    asymmetry_fact[:]        = aeropt.asy

    single_sca_alb           = ds.createVariable("single_scatter_albedo", 'f4', s3bin)
    single_sca_alb.units     = "-"
    single_sca_alb.long_name = "single scatter albedo"
    single_sca_alb[:]        = aeropt.ssa

    extinction               = ds.createVariable("extinction", 'f4', s3bin)
    extinction.units         = "m^2/kg"
    extinction.long_name     = "mass extinction coefficient (not scaled)"
    extinction[:]            = aeropt.ext # .* aer.ext_scaling

    phase_fun                = ds.createVariable("phase_function", 'f4', s4bin)
    phase_fun.units          = "-"
    phase_fun.long_name      = "scattering phase function"
    phase_fun.description    = ("normalized so that pfun(cos(ang))*d(cos(ang)) between 0 "
                              +"and pi is =2 an the integral on the full solid angle is 4*PI")
    phase_fun[:]             = aeropt.pfun

    wavelength               = ds.createVariable("wavelength", 'f4', ("wavelength",))
    wavelength.units         = "m"
    wavelength.long_name     = "monochromatic wavelength"
    wavelength[:]            = aer.lambda_out
  
    rel_hum                  = ds.createVariable("rel_hum", 'f4', ("rel_hum",))
    rel_hum.units            = "%"
    rel_hum.long_name        = "ambient relative humidity"
    rel_hum[:]               = aer.rh_tab # note that this is asummed to be
                                          # identical for all components in the mixture.

    rel_hum_growth           = ds.createVariable("rel_hum_growth", 'f4', s2hum)
    rel_hum_growth.units     = "-"
    rel_hum_growth.long_name = "relative humidity growth factor"
    for i in range(num_components):
        rel_hum_growth[:,i] = laer_conf[i].rh_growth
    
    size_bin_min             = ds.createVariable("size_bin_min", 'f4',  s2com)
    size_bin_min.units       = "m"
    size_bin_min.long_name   = "lower limit of the distribution size bin"
    for i in range(num_components):
        size_bin_min[:,i] = laer_conf[i].bins_min
    

    size_bin_max             = ds.createVariable("size_bin_max", 'f4', s2com)
    size_bin_max.units       = "m"
    size_bin_max.long_name   = "upper limit of the distribution size bin"
    for i in range(num_components):
        size_bin_max[:,i] = laer_conf[i].bins_max
    
    angle           = ds.createVariable("angle", 'f4', ("angle",))
    angle.units     = "degrees"
    angle.long_name = "phase function angle"
    angle[:]        = aeropt.pfun_ang

    ds.close()

    return 

