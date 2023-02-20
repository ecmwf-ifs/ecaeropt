

 ###########################################################################################
 #                                                                                         #
 # aeropt/aer.py                                                                           #
 #                                                                                         #
 # author: Ramiro Checa-Garcia                                                             #
 # email:  ramiro.checa-garcia@ecmwf.int                                                   #
 #                                                                                         #
 # history:                                                                                #
 #    - Oct-2022 [Ramiro Checa-Garcia]      1st tested version                             #
 #                                                                                         #
 # info:                                                                                   #
 #        CLASSES                                                                          #
 #        * aerosol         : one has all the information of the aerosol configuration.    #
 #        * aeropt          : encapsulate the single aerosol optical properties            #
 #                                                                                         #
 #        FUNCTIONS                                                                        #
 #        * readconf        : creates an aerosol object based in configuration file        #
 #        * mie_to_aeropt   : creates an aeropt object from outputs of a mie code          #                
 #        * mixing          : calculates an external mixing aerosol given the aeropt of    #
 #                            the components                                               #
 ###########################################################################################


import numpy as np
import toml


class aeropt:
    """Class to store aerosol optical properties

    Attributes:
       kind  (str)     : lognormal, mixing-logn 
       ri_r  (np.array): real      part refractive index 
       ri_i  (np.array): imaginary part refractive index 
       ext   (np.array): extinction 
       ssa   (np.array): single scattering albedo
       asy   (np.array): asysmetry factor
       mass  (np.array): 
       lidar (np.array):
       pfun  (np.array): phase function 
       angle (np.array): angles used for phase function
       nmux  (int)     :

    """
    def __init__(self, kind, ri_r, ri_i, ext, ssa, asy, mass, lidar,
                 pfun, pfun_ang, nmux, engine):

        self.kind     = kind  # lognormal, or mixing-logn
        self.ri_r     = ri_r
        self.ri_i     = ri_i
        self.ext      = ext
        self.ssa      = ssa
        self.asy      = asy
        self.mass     = mass
        self.lidar    = lidar
        self.pfun     = pfun
        self.pfun_ang = pfun_ang
        self.nmux     = nmux
        self.engine   = engine


def mie_to_aeropt(aerconf, out, engine):
    """
    Given the output of the mie-scatteting engine, and the aerosols configuration
    it returns a aerosol optical properties object.

    Args:
      aerconf (object) : aerosol configuration obj.
      out     (tuple)  : output of mie-scattering engine is a tuple of different 
                        arrays.
      engine   (string): string with the name of the tuple.

    Returns:

      aeropt  (object) : returns an aerosols optical properties object.

    """
    
    def _check_scaling(lbins):
       nlbin = lbins[-1]
       if len(lbins)%nlbin==0:
           return True
       else:
           return False
       return


    znr, zni, ext, omg, asy, lidar, mass, ph, qval = out 

    if aerconf.scaling_file != None:
        scaling = np.loadtxt(aerconf.scaling_file).T
        if _check_scaling(scaling[4]):
           bins    = scaling[4] 
           nbin    = int(bins[-1])
           nwl     = int(len(bins)/nbin)
           sca_wl  = scaling[0].reshape(nwl,nbin) * 1.e-6
           sca_ext = scaling[1].reshape(nwl,nbin)
           sca_ssa = scaling[2].reshape(nwl,nbin)
           sca_asy = scaling[3].reshape(nwl,nbin)
           wl      = aerconf.lambda_out
           print(wl)
           s_ext   = np.zeros((len(wl), nbin))
           s_ssa   = np.zeros((len(wl), nbin))
           s_asy   = np.zeros((len(wl), nbin))
           for nb in range(nbin):
               s_ext[:,nb]  = np.interp(wl, sca_wl[:,nb], sca_ext[:,nb])
               s_ssa[:,nb]  = np.interp(wl, sca_wl[:,nb], sca_ssa[:,nb])
               s_asy[:,nb]  = np.interp(wl, sca_wl[:,nb], sca_asy[:,nb])
           print(s_ext.shape, ext.shape)
           print(ext)
           
           if len(ext.shape) == 3:
              for nrh in range(ext.shape[1]):
                  ext[:,nrh,:] = ext[:,nrh,:] * s_ext[:,:]
                  omg[:,nrh,:] = omg[:,nrh,:] * s_ssa[:,:]
                  asy[:,nrh,:] = asy[:,nrh,:] * s_asy[:,:]
           else:
               ext[:,:] = ext[:,:] * s_ext[:,:] 
               omg[:,:] = omy[:,:] * s_ssa[:,:]
               asy[:,:] = asy[:,:] * s_asy[:,:]

    return aeropt( aerconf.kind, znr, zni, ext, omg, asy, mass, lidar,
                   ph, aerconf.angles, aerconf.nmumax, engine)

class aerosol:
    """Class to store a complete aerosol configuration

    Attributes:

       kind         (string)  : kind or family distribution (at this moment only lognormal)
       r0           (array)   : r0 : mean radius such as mu=ln r0 (one value per mode)
       sigma_g      (array)   : sigma_g : standard deviation of distribution (one value per mode)
       Ntot         (array)   : Ntot : Total number concentration (one value per mode)
       rho          (float)   : rho : density of particles in [kg m-3] (same for all modes)
       bins_min     (array)   : bins_min : lower bounds of bins 
       bins_max     (array)   : bins_max : upper bounds of bins
       rh_tab       (list)    : rh_tab   : relative humidity values
       rh_growth    (list)    : rh_growth : growth of particles for the relative humidity values (rh_tab)
       ri_file      (string)  : ri_file : plain text file (csv with spaces as delimiters) with refractive index per wave lenght.
       aer_type     (string)  : aer_type : aerosol type
       ext_scaling  (float)   : ext_scaling : scaling introduced in extinction after calculation
       lambda_out   (array)   : output wavelenghts : wavelength for optical properties (can be different from ri file)
       ri_nrh       (       ) : ri_nrh
       NInp         (       ) : NInp
       Ntot_mix     (       ) : Ntot_mix
       nb_lambda    (integer) : number of output wavelenghts (it would be the length of lambda_out)
       Ndis         (       ) : Ndis (number of modes in the lognormal distribution)
       size_bins    (float)   : size_bins (number of size bins)
       rh_int       (integer) : size relative humidity arrays 
       config_file  (string)  : filepath of config file : where info and rules for optical properties are detailed
       wl_file      (string)  : filepath of output wavelenghts : lambda_out and nb_lambda are taken from this file.
       nmumax       (      )  : nmumax
       angles       (array)   : angles for phase function calculations
       ri_lambdatab (array)   : tabulated wavelengths of refr.i. : input wavelenghts for species 
       znr_tab      (array)   : tabulated real  refractive index : input real refractive index for species
       zni_tab      (array)   : tabulated imag. refractive index : input imaginary refractive index for species.



    """

    def __init__(self, kind, r0 , sigma_g, Ntot, rho, bins_min, bins_max, rh_tab, rh_growth, 
                 ri_file  , aer_type, ext_scaling, lambda_out, ri_nrh, NInp, Ntot_mix,
                 nb_lambda, Ndis, size_bins, rh_int, config_file, wl_file, nmumax, angles,
                 ri_lambdatb, znr_tab, zni_tab, scaling_file ):

       # It is good to group these also in classes: distribution, bins, properties,....
       self.kind        = kind
       self.r0          = r0
       self.sigma_g     = sigma_g
       self.Ntot        = Ntot
       self.rho         = rho
       self.bins_min    = bins_min
       self.bins_max    = bins_max
       self.rh_tab      = rh_tab
       self.rh_growth   = rh_growth
       self.ri_file     = ri_file
       self.aer_type    = aer_type
       self.ext_scaling = ext_scaling
       self.lambda_out  = lambda_out
       self.ri_nrh      = ri_nrh
       self.NInp        = NInp
       self.Ntot_mix    = Ntot_mix
       self.nb_lambda   = nb_lambda
       self.Ndis        = Ndis
       self.size_bins   = size_bins
       self.rh_int      = rh_int
       self.config_file = config_file
       self.wl_file     = wl_file
       self.nmumax      = nmumax
       # It is reasonable to add this information to aerosol-class
       self.angles       = angles
       self.ri_lambdatab= ri_lambdatb
       self.znr_tab     = znr_tab
       self.zni_tab     = zni_tab
       # Here we can add a new item: self.scaling => which is one file per bin.
       # this will not be used until we estimate optical properties
       # so it is after run the engine.
       self.scaling_file= scaling_file

    def __str__(self):
       a = []

       a.append("\n ===== AEROSOL CONFIGURATION ==================")
       a.append("\n **** DISTRIBUTION: ")
       a.append("        KIND            : "+str(self.kind        ))
       a.append("        (r0)            : "+str(self.r0          ))
       a.append("        (sigma_g)       : "+str(self.sigma_g     ))
       a.append("        (Ntot)          : "+str(self.Ntot        )) 
       a.append("\n **** BINDS ")
       a.append("        min bins array  : "+str(self.bins_min  ))
       a.append("        max bins array  : "+str(self.bins_max  ))
       a.append("\n **** PROPERTIES ")
       a.append("        rho             : "+str(self.rho         )) 
       a.append("        rh.tab array    : "+str(self.rh_tab      )) 
       a.append("        rh.growth array : "+str(self.rh_growth   )) 
       a.append("\n **** OPTICAL ")
       a.append("        ref. index file : "+str(self.ri_file     )) 
       a.append("        extinction sca. : "+str(self.ext_scaling )) 
       a.append("        lambda out [0]. : "+str(self.lambda_out[0])) 
       a.append("        lambda out [1]. : "+str(self.lambda_out[1])) 
       a.append("        lambda out dim. : "+str(self.lambda_out.shape))
       a.append("        ri_nrh          : "+str(self.ri_nrh)) 
       a.append("        NInp            : "+str(self.NInp        )) 
       a.append("        nb_lambda       : "+str(self.nb_lambda   )) 

       a.append("\n **** TAGS ")
       a.append("        aer_type        : "+ str(self.aer_type    ))


       a.append("\n **** ANGLES ")
       a.append("        min angle       : "+str(np.min(self.angles)      )) 
       a.append("        max angle       : "+str(np.max(self.angles)      )) 
       a.append("        number angles   : "+str(self.angles.shape        )) 
       

       a.append("\n **** OTHERS ")

       a.append("        Ntot_mix        : "+ str(self.Ntot_mix    ))  
       a.append("        Ndis            : "+ str(self.Ndis        )) 

       a.append("        config-file     : "+str(self.config_file )) 
       a.append("        wl-file         : "+str(self.wl_file )) 
       a.append("\n **** SIZES ")
       a.append("        num. bins       : "+str(self.size_bins   )) 
       a.append("        num. rh         : "+str(self.rh_int      )) 
       a.append("        n-mu-max        : "+str(self.nmumax      )) 

       a.append("\n **** REFRACTIVE INDEX ")
       a.append("        ri lambda  dims : " + str(self.ri_lambdatab.shape)) 
       a.append("        ri real    dims : " + str(self.znr_tab.shape     )) 
       a.append("        ri imag    dims : " + str(self.zni_tab.shape     )) 
       if self.scaling != None:
            a.append("\n **** SCALING PRESENT ")
       return "\n".join(a)



def readconf(config_file, angles, wl_out="none", debug=False):
    """
     Return an aerosol configuration object given a configuration file,
     scattering angles and list of wavelengths.

     Args:
        config_file (string): file path with configuration toml file.
        angles      (array) : array of angles or file path with them.
        wl_out      (string): name of wl_out file with wavelengths
        debug       (bool)  : True shows debug info.

     Return:
        aerosol     (object): aerosol configuration object

    """


    def as_farray(value):

        if isinstance(value, float) or isinstance(value, int):
            arval = np.array([value], dtype=float)
        else:
            arval = np.array(value, dtype=float)

        return arval


    conf = toml.load(config_file, _dict=dict)

    if "scaling_file" in conf["optical"]:
        scaling_file = conf["optical"]["scaling_file"]
    else:
        scaling_file = None

    # Distribution ====================================
    dist = conf["distribution"]["lognormal"]
    kind = "lognormal"

    r0      = as_farray(dist["r0"])
    sigma_g = as_farray(dist["sigma_g"])
    Ntot    = as_farray(dist["Ntot"])

    # Bins ============================================
    bins = conf["bins"]

    bins_min = as_farray(bins["min"])
    bins_max = as_farray(bins["max"])

    # Physical properties =============================
    prop = conf["properties"]

    rho       = prop["rho"]
    rh_tab    = as_farray(prop["rh"]["tab"])
    rh_growth = as_farray(prop["rh"]["growth"])

    # Optical properties ==============================
    opti = conf["optical"]

    ext_scaling = as_farray(opti["ext_scaling"])
    ri_file     = opti["ri_file"]

    # Tags ============================================
    tags = conf["tags"]

    aer_type = tags["aer_type"]
    config_file = config_file

    # Mixture =========================================
    if "mixture" in conf.keys():
        mixt = conf["mixture"]
        Ntot_mix = as_farray(mixt["Ntot_mix"])
    else:
        Ntot_mix = np.array([0.0])


    # Processing linked files of configuration

    if wl_out == 'none':
        wl_file = opti["wl_out_file"]
    else:
        wl_file = wl_out

    lambda_out = np.loadtxt(wl_file, skiprows=1)

    nb_lambda  = len(lambda_out)

    ri_nrh=len(rh_tab)

    # INFO:
    #
    # Note the file with refractive index repeat the wavelength as many
    # times as we have ri_nrh this is controlled by the variable ninp
    # which is the number of input lambdas in the refractive index file
    #
    # - soa_a_ri_table_mix
    #   wc -l = 1898 => 1898-2 = 1896
    # - rh_tab in config file -> 12
    #   1896/12 = 158 with is the ninp => number wavelengths in input file

    ni_data = np.loadtxt(ri_file, skiprows=2).T

    if ni_data.ndim > 1:
        if ni_data.shape[1]>3: # we have rh included
            NInp=int(len(ni_data[0,:])/ri_nrh)
        else:
            NInp=1
        nwl = int(len(ni_data[0,0:NInp]))
        ri_lambdatab  = ni_data[0,0:nwl]
        znr_tab  = ni_data[1,:].reshape(ri_nrh, nwl).T
        zni_tab  = ni_data[2,:].reshape(ri_nrh, nwl).T
    else:
        NInp=1
        nwl = 1
        ri_lambdatab  = ni_data[0]
        znr_tab  = ni_data[1]
        zni_tab  = ni_data[2]

    nmumax = len(angles)
    rh_int = len(rh_tab)

    Ndis   = len(r0)
    size_bins = len(bins_min)
    


    return aerosol(kind, r0 , sigma_g, Ntot, rho, bins_min, bins_max, rh_tab,
                    rh_growth, ri_file  , aer_type, ext_scaling, lambda_out,
                    ri_nrh, NInp, Ntot_mix, nb_lambda, Ndis, size_bins, rh_int,
                    config_file, wl_file, nmumax, np.array(angles), ri_lambdatab,
                    znr_tab, zni_tab, scaling_file)




def mixing(mix_aer, mix_opt, mix_lbtab, mix_ri_rtab, mix_ri_itab, num_components):
    """
    Calculates the mixture optical properties given the optical properties of the components

    Args:
       mix_aer
       mix_opt
       mix_lbtab
       mix_ri_rtab
       mix_ri_itab    ()
       num_components (integer)


    Return:
       aer_mix_opt (aeropt) : aeropt object with optical properties of the mixture

    """


    aer = mix_aer[1] # This aerconf object

    s3   = ( aer.nb_lambda, aer.rh_int, aer.size_bins )
    s2   = (                aer.rh_int, aer.size_bins )
    s4   = ( aer.nb_lambda, aer.rh_int, aer.size_bins, aer.nmumax )
    smix = ( aer.nb_lambda, aer.rh_int, num_components)


    # accumulating arrays -------
    ext_acc   = np.zeros(s3)
    sca_acc   = np.zeros(s3)
    asy_acc   = np.zeros(s3)
    mass_acc  = np.zeros(s2)
    lidr_acc  = np.zeros(s3) 
    pfun_acc  = np.zeros(s4)
    pfun_ang_acc = np.zeros(aer.nmumax)

    # final mixing arrays --------
    ext_mix   = np.zeros(s3)
    sca_mix   = np.zeros(s3)
    asy_mix   = np.zeros(s3)
    mass_mix  = np.zeros(s2)
    lidr_mix  = np.zeros(s3) 
    pfun_mix  = np.zeros(s4)
    pfun_ang_mix = np.zeros(aer.nmumax)
    nmux_mix = 0.0

    ri_r_mix = np.zeros(smix)
    ri_i_mix = np.zeros(smix) 

    # Per component arrays ---------

    ext_com   = np.zeros(s3)
    sca_com   = np.zeros(s3)
    asy_com   = np.zeros(s3)
    mass_com  = np.zeros(s2)
    lidr_com  = np.zeros(s3) 
    pfun_com  = np.zeros(s4)
    nmux_com = 0.0

    # Revise these loops and write with python index and numpy broadcasting
    for i_com in range(num_components):

        sca_com  = mix_opt[i_com].ext  * mix_opt[i_com].ssa
        ext_com  = mix_opt[i_com].ext
        ssa_com  = mix_opt[i_com].ssa
        asy_com  = mix_opt[i_com].asy
        mass_com = mix_opt[i_com].mass * mix_aer[i_com].Ntot_mix
        lidr_com = mix_opt[i_com].lidar
        pfun_com = mix_opt[i_com].pfun


        for j in range(aer.nb_lambda):
            ext_acc[ j,:,:] = ext_acc[ j,:,:] + ext_com[ j,:,:]  * mass_com[ :,:]
            sca_acc[ j,:,:] = sca_acc[ j,:,:] + sca_com[ j,:,:]  * mass_com[ :,:]
            asy_acc[ j,:,:] = asy_acc[ j,:,:] + asy_com[ j,:,:]  * sca_com[j,:,:] * mass_com[:,:]
            lidr_acc[j,:,:] = lidr_acc[j,:,:] + lidr_com[j,:,:]  * sca_com[j,:,:] * mass_com[:,:]

        mass_acc[:,:] = mass_acc[:,:] + mass_com[:,:]

        pfun_ang_mix = mix_opt[i_com].pfun_ang
        nmux_mix     = mix_opt[i_com].nmux

        ri_r_mix[:,:, i_com] = mix_opt[i_com].ri_r
        ri_i_mix[:,:, i_com] = mix_opt[i_com].ri_i

        for j in range(aer.nb_lambda):
            for l in range(len(pfun_ang_mix)):
                pfun_acc[j,:,:,l] = pfun_acc[j,:,:,l] + pfun_com[j,:,:,l] * sca_com[j,:,:] * mass_com[:,:]

        for j in range(aer.nb_lambda):
            ext_mix[ j, :, :] = ext_acc[ j, :, :] / mass_acc[ :, :]
            for l in range(len(pfun_ang_mix)):
                  pfun_mix[ j, :, :, l] = pfun_acc[ j, :, :, l] / sca_acc[ j, :, :]

    ssa_mix   = sca_acc  / ext_acc
    asy_mix   = asy_acc  / sca_acc
    lidr_mix  = lidr_acc / sca_acc

    mass_mix  = mass_acc # think on that not in original code.


    #print("\n         Calculated mixed aerosol\n")

    # Here we should create aer_mix_opt as an aerosol optics object

    aer_mix_opt = aeropt("mixing-logn", ri_r_mix, ri_i_mix,
                         ext_mix, ssa_mix, asy_mix, mass_mix,
                         lidr_mix, pfun_mix, pfun_ang_mix, nmux_mix, mix_opt[0].engine)

    return mix_aer, aer_mix_opt

