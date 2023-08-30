

module mie_code

! +----------------------------------------------------------------------------------------+
! | engines/mie_code/src/mie_aerosols_struct.f90                                           |
! |                                                                                        |
! | (C) Copyright 1996-2023 CNRS (Author: Olivier Boucher)                                 |
! | (C) Copyright 2022- ECMWF                                                              |
! |                                                                                        |
! | This software is licensed under the terms of the Apache Licence Version 2.0            |
! | which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.                   |
! |                                                                                        |
! | In applying this licence, ECMWF does not waive the privileges and immunities           |
! | granted to it by virtue of its status as an intergovernmental organisation             |
! | nor does it submit to any jurisdiction.                                                |
! |                                                                                        |
! | Author:                                                                                |
! |    Olivier Boucher.      CNRS    Original code for Mie homogenous spheres              |
! |    Alessio Bozzo.       ECMWF    1st implementation in ECMWF                           |
! |    Ramiro Checa-Garcia. ECMWF    Changes/improvements for ecaeropt                     |
! |                                                                                        |
! | Modifications:                                                                         |
! |    10-Oct-2022   R. Checa-Garcia  The original code of Boucher/Bozzo has been changed  |
! |                                   to a module and it is created a wrapper of the       |
! |                                   fortran subroutine using iso_c_bindings which allows |
! |                                   the access of the mie_aerosols subroutine from       |
! |                                   different languages                                  |
! |    26-Nov-2022   R. Checa-Garcia  using iso_fortran_env / code tests also in Julia.    |
! |    28-Nov-2022   R. Checa-Garcia  Refactor to an structured code. Tests for an,bn done |
! |                                   Test for single sphere. Distribution methods         |
! |                                   programmed as functions. Modif. refract. index input |
! |                                   (now array)                                          |
! |                                                                                        |
! +----------------------------------------------------------------------------------------+
  use iso_c_binding
  use iso_fortran_env

  use parkind1, only : JPIM, JPRB 
  implicit none
  
  contains


  subroutine mie_aerosols( &
          & Ninp, Nmumax, nb_lambda, lambda_out, Ndis, size_bins, &   ! 6
          & size_bins_min, size_bins_max,                         &   ! 2
          & sigma_g, r_0_inp, Ntot, rho,          ri_nrh,         &   ! 5
          & rh_int, rh_tab, rh_growth, cutoff_radius, verbose,    &   ! 5
          & zlambtab, z_n_r_tab, z_n_i_tab,                       &   ! 3 => 22 in
          & znr,zni,                                              &   ! 2
          & ext_stored, omg_stored, asy_stored,                   &   ! 3 + 4 below => 9
          & lidar_ratio_stored, mass_stored, ph_stored, ph_ang,   &
          & test_single_sphere)   bind(c, name="mie_aerosols")


        implicit none

        integer(KIND=c_int) , intent(IN) :: size_bins
        integer(KIND=c_int) , intent(IN) :: rh_int
        integer(KIND=c_int) , intent(IN) :: Ndis
        integer(KIND=c_int) , intent(IN) :: nb_lambda
        integer(KIND=c_int) , intent(IN) :: ri_nrh
        integer(KIND=c_int) , intent(IN) :: NInp, nmumax, verbose
        real(KIND=c_double) , intent(IN) :: lambda_out(nb_lambda)
        real(KIND=c_double) , intent(IN) :: size_bins_min(size_bins)
        real(KIND=c_double) , intent(IN) :: size_bins_max(size_bins)
        real(KIND=c_double) , intent(IN) :: sigma_g(Ndis)
        real(KIND=c_double) , intent(IN) :: r_0_inp(Ndis)
        real(KIND=c_double) , intent(IN) :: Ntot(Ndis)
        real(KIND=c_double) , intent(IN) :: rho
        real(KIND=c_double) , intent(IN) :: rh_tab(rh_int)
        real(KIND=c_double) , intent(IN) :: rh_growth(rh_int)
        real(KIND=c_double) , intent(IN) :: cutoff_radius  !mass cut-off radius in microns
        real(KIND=c_double) , intent(IN) :: zlambtab(NInp)
        real(KIND=c_double) , intent(IN) :: z_n_r_tab(NInp,ri_nrh), z_n_i_tab(NInp,ri_nrh)

        ! output quantities (per specie size bin ; per wavelength ; per rh(if present))
        real(KIND=c_double) , intent(OUT)  :: znr(nb_lambda,rh_int),zni(nb_lambda,rh_int)
        real(KIND=c_double) , intent(OUT)  :: ext_stored(nb_lambda,rh_int,size_bins)
        real(KIND=c_double) , intent(OUT)  :: omg_stored(nb_lambda,rh_int,size_bins)
        real(KIND=c_double) , intent(OUT)  :: asy_stored(nb_lambda,rh_int,size_bins)
        real(KIND=c_double) , intent(OUT)  :: mass_stored(rh_int,size_bins), test_single_sphere(4)
        real(KIND=c_double) , intent(OUT)  :: lidar_ratio_stored(nb_lambda,rh_int,size_bins)
        real(KIND=c_double) , intent(OUT)  :: ph_stored(nb_lambda,rh_int,size_bins,Nmumax)
        real(KIND=c_double) , intent(INOUT):: ph_ang(Nmumax) ! this has been changed and we have to include
                                                              ! input angles in the code

        ! verbose ==2 print some information to check size and values between
        !                   external code and called subrotines
        if (verbose.eq.2) then
            print *, "=== FORTRAN-WRAPPER ==========================================="
            print *, " testing size arrays inside F90 sub."
            print *, ' znr                ', shape(znr)
            print *, ' zni                ', shape(zni)
            print *, ' ext_stored         ', shape(ext_stored)
            print *, ' omg_stored         ', shape(omg_stored)
            print *, ' asy_stored         ', shape(asy_stored)
            print *, ' mass_stored        ', shape(mass_stored)
            print *, ' lidar_ratio_stored ', shape(lidar_ratio_stored)
            print *, ' ph_stored          ', shape(ph_stored)
            print *, ' ph_ang             ', shape(ph_ang)
            print *, "=============================================================="
        endif


        if (verbose.eq.2) then
            print *, "=== FORTRAN-WRAPPER =========================================="
            print *, " testing some variables inside F90 sub."
            print *, ' NInp       ', NInp
            print *, ' Nmumax     ', Nmumax
            print *, ' nb_lambda  ', nb_lambda
            print *, ' lambda_out ', lambda_out(:)
            print *, ' Ndis       ', Ndis
            print *, ' size_bins  ', size_bins
            print *, ' ri_nrh     ', ri_nrh
            print *, ' ri_lamb_tab', zlambtab(:)
            print *, ' rho        ', rho
            print *, ' size_bins  ', size_bins
            print *, "============================================================="
        endif


        call main_mie_aerosols( &
            & Ninp, Nmumax, nb_lambda, lambda_out, Ndis, size_bins, &   ! 6
            & size_bins_min, size_bins_max,                         &   ! 2
            & sigma_g, r_0_inp, Ntot, rho,          ri_nrh,         &   ! 5
            & rh_int, rh_tab, rh_growth, cutoff_radius, verbose,    &   ! 5
            & zlambtab, z_n_r_tab, z_n_i_tab,                       &   ! 3 => 21 in
            & znr,zni,                                              &   ! 2
            & ext_stored, omg_stored, asy_stored,                   &   ! 3
            & lidar_ratio_stored, mass_stored, ph_stored, ph_ang,   &
            & test_single_sphere)       ! 4 => 9  out)



        if (verbose.eq.2) then
            print *, "COMPLETED (fortran wrapper) ====================="
            print *, '- Fortran wrapper -------------------------------'
            print *, ' znr (1,1)'     , znr(1,1)
            print *, ' zni (1,1)'     , zni(1,1)
            print *, ' ext (1,1,1)'   , ext_stored(1,1,1)
            print *, ' omg (1,1,1)'   , omg_stored(1,1,1)
            print *, ' asy (1,1,1)'   , asy_stored(1,1,1)
            print *, ' lidar (1,1,1)' , lidar_ratio_stored(1,1,1)
            print *, ' pfun (1,1,1,1)', ph_stored(1,1,1,1)
            print *, '-------------------------------------------------'
        endif

  end subroutine mie_aerosols


    
  subroutine main_mie_aerosols( &
      & Ninp, Nmumax, nb_lambda, lambda_out, Ndis, size_bins, &   ! 6
      & size_bins_min, size_bins_max,                         &   ! 2
      & sigma_g, r_0_inp, Ntot, rho,          ri_nrh,         &   ! 5
      & rh_int, RH_tab, rh_growth, cutoff_radius, verbose,    &   ! 5
      & zlambtab, z_n_r_tab, z_n_i_tab,                       &   ! 3 => 21 in
      & znr,zni,                                              &   ! 2
      & ext_stored, omg_stored, asy_stored,                   &   ! 3
      & lidar_ratio_stored, mass_stored, ph_stored, ph_ang,   &
      & test_single_sphere)       ! 4 => 9  out


      USE PARKIND1 , ONLY : JPIM, JPRB

      IMPLICIT NONE

  !-------Mie computations for a size distribution of homogeneous spheres---------------------------
  !
  !-------Ref : Toon and Ackerman, Applied Optics, 1981
  !             Stephens, CSIRO, 1979
  !
  !------------------------------------------------------------------------------------------------
  !
  ! Input/configuration parameters:
  !
  !  NInp      : number of reference wavelengths in the refractive index file
  !  nb_lambda : number of output  wavelengths
  !  Ndis      : number of distributions considered (mono modal or bi-modal)
  !  size_bins : number of size bins considered for the specie
  !
  !----------------------------------------------------------------------------------------------
  ! SIZE distribution properties              => it considers a log-normal distribution
  !  sigma_g   : geometric standard deviation
  !  r_0       : geometric number mean radius (um)/modal radius (units um)!!
  !  Ntot      : total concentration in m-3
  !
  !  rho       : density in kgm-3
  !----------------------------------------------------------------------------------------------
  !
  ! OUTPUT
  !
  ! znr               (nb_lambda,rh_int)
  ! zni               (nb_lambda,rh_int)
  ! ext_stored        (nb_lambda,rh_int,size_bins)            units:
  ! omg_stored        (nb_lambda,rh_int,size_bins)            units:
  ! asy_stored        (nb_lambda,rh_int,size_bins)            units:
  ! lidar_ratio_stored(nb_lambda,rh_int,size_bins)            units:
  ! mass_stored       (rh_int,size_bins)                      units:
  ! ph_stored         (nb_lambda,rh_int,size_bins,0:Nmumax)
  !
  !----------------------------------------------------------------------------------------------

    INTEGER(4), INTENT(IN) :: size_bins
    INTEGER(4), INTENT(IN) :: rh_int
    INTEGER(4), INTENT(IN) :: Ndis
    INTEGER(4), INTENT(IN) :: nb_lambda
    INTEGER(4), INTENT(IN) :: ri_nrh
    INTEGER(4), INTENT(IN) :: NInp, nmumax, verbose
    REAL(8)   , INTENT(IN) :: lambda_out(nb_lambda)
    REAL(8)   , INTENT(IN) :: size_bins_min(size_bins)
    REAL(8)   , INTENT(IN) :: size_bins_max(size_bins)
    REAL(8)   , INTENT(IN) :: sigma_g(Ndis)
    REAL(8)   , INTENT(IN) :: r_0_inp(Ndis)
    REAL(8)   , INTENT(IN) :: Ntot(Ndis)
    REAL(8)   , INTENT(IN) :: rho
    REAL(8)   , INTENT(IN) :: RH_tab(rh_int)
    REAL(8)   , INTENT(IN) :: rh_growth(rh_int)
    REAL(8)   , INTENT(IN) :: cutoff_radius                   !mass cut-off radius in microns
    REAL(8)   , INTENT(IN) :: zlambtab(NInp) 
    REAL(8)   , INTENT(IN) :: z_n_r_tab(NInp,ri_nrh), z_n_i_tab(NInp,ri_nrh)
    ! intent(in) 19 objects

    ! output quantities (per specie size bin, per wavelength, per rh(if present))
    REAL(8),    INTENT(OUT)    :: znr(nb_lambda,rh_int),zni(nb_lambda,rh_int)
    REAL(8),    INTENT(OUT)    :: ext_stored(nb_lambda,rh_int,size_bins)
    REAL(8),    INTENT(OUT)    :: omg_stored(nb_lambda,rh_int,size_bins)
    REAL(8),    INTENT(OUT)    :: asy_stored(nb_lambda,rh_int,size_bins)
    REAL(8),    INTENT(OUT)    :: mass_stored(rh_int,size_bins)
    REAL(8),    INTENT(OUT)    :: lidar_ratio_stored(nb_lambda,rh_int,size_bins)
    REAL(8),    INTENT(OUT)    :: ph_stored(nb_lambda,rh_int,size_bins,Nmumax)
    REAL(8),    INTENT(OUT)    :: test_single_sphere(4)
    REAL(8),    INTENT(INOUT)  :: ph_ang(Nmumax)

    ! INTERNAL VARIABLES
    INTEGER(4)             :: Nbin_points
    REAL(8)                :: rmin,rmax,r_0(Ndis)

    INTEGER(4) , PARAMETER :: Nsize =210000  !--- reserved number for complex iterative functions  
    INTEGER(4) , PARAMETER :: Nsize2=210100  !    additinal space for specific arrays

    COMPLEX(KIND=JPRB)     :: m              !--- refractive index m=n_r-i*n_i
    INTEGER(4)             :: Nmax           !--- number of iterations for the Wiscombe algorithm
    COMPLEX(KIND=JPRB)     :: a(1:Nsize), b(1:Nsize) 

    !-- additional variables and arrays for computation of the
    !-- phase function and lidar ratio

    INTEGER(4) :: Nmu   !, Nmumax
    REAL(8)    :: volume, surface, rgv, rgn
    REAL(8)    :: tot_num_concentration, radius, radius2, mass
    REAL(8)    :: sigma_abs,   Q_abs
    REAL(8)    :: PP(1:Nmumax), phase(1:Nmumax)
    REAL(8)    :: Q_ext, Q_sca, g, omega   !--parameters for radius r
    REAL(8)    :: x  !--size parameter
    REAL(8)    :: r  !--radius
    REAL(8)    :: sigma_sca, sigma_ext, omegatot,  gtot !--averaged parameters
    REAL(8)    :: numb, deltar
    !-- interpolation of refractive indices on output wavelengths
    REAL(8)    :: n_r(nb_lambda)   , n_i(nb_lambda)

    ! intermediate values
    REAL(8) :: final_a(nb_lambda), final_g(nb_lambda), final_w(nb_lambda)
    REAL(8) :: final_lratio(nb_lambda) !, final_PH(nb_lambda,0:Nmumax)

    !general looping indices and utility parameters
    INTEGER(4) :: Nwv, IRH, dis, bin,dbin 
    REAL(8)    :: pi,  NR, counter

    !----------------------------------------------------------------------------------------------
    !-- START OF COMPUTATIONS

    !total number of iterations (bins*RH*wl)
    NR=real(nb_lambda*rh_int*size_bins)
    counter=0


    Nbin_points=999
    test_single_sphere(:)=0.0

    !! refractve index input file
    !! expected 3 columns: wl(m), real part, im part
    !! 2 lines of header

    if (verbose.eq.2) then
        print *, "=== FORTRAN ===================================================="
        print *, " testing size arrays inside F90 sub."
        print *, ' znr                ', shape(znr)
        print *, ' zni                ', shape(zni)
        print *, ' ext_stored         ', shape(ext_stored)
        print *, ' omg_stored         ', shape(omg_stored)
        print *, ' asy_stored         ', shape(asy_stored)
        print *, ' mass_stored        ', shape(mass_stored)
        print *, ' lidar_ratio_stored ', shape(lidar_ratio_stored)
        print *, ' ph_stored          ', shape(ph_stored)
        print *, ' ph_ang             ', shape(ph_ang)
        print *, "==============================================================="
    endif


    if (verbose.eq.2) then
        print *, "=== FORTRAN ==================================================="
        print *, " testing some variables inside F90 sub."
        print *, ' NInp       ', NInp
        print *, ' Nmumax     ', Nmumax
        print *, ' nb_lambda  ', nb_lambda
        print *, ' lambda_out ', lambda_out(:)
        print *, ' Ndis       ', Ndis
        print *, ' size_bins  ', size_bins
        print *, ' z_n_r_tab  ', shape(z_n_r_tab)
        print *, ' ri_nrh     ', ri_nrh
        print *, ' ri_lamb_tab', zlambtab(1), zlambtab(2) 
        print *, ' rho        ', rho
        print *, ' size_bins  ', size_bins
        print *, "==============================================================="
    endif

    !OMP_set_num_thread(4)

  DO IRH=1,rh_int
     if (NInp.gt.1) then
        call interp_ri(NInp,zlambtab, z_n_r_tab(:,IRH), z_n_i_tab(:,IRH), nb_lambda,lambda_out,znr(:,IRH),zni(:,IRH),verbose)
     else
        znr(:,IRH)=z_n_r_tab(1,IRH)
        zni(:,IRH)=z_n_i_tab(1,IRH)
     endif
  ENDDO


  DO bin=1,size_bins

     DO IRH=1,rh_int

        if (rh_int .eq. 1) then
           DO Nwv=1, nb_lambda
              n_r(Nwv)=znr(Nwv,1)
              n_i(Nwv)=zni(Nwv,1)

              rmin=size_bins_min(bin)*1.E-6_JPRB
              rmax=size_bins_max(bin)*1.E-6_JPRB
             DO dis=1,Ndis
                 r_0(dis)=r_0_inp(dis)*1.E-6_JPRB
              ENDDO
           ENDDO
        else
           DO Nwv=1,nb_lambda
              n_r(Nwv)=znr(Nwv,IRH)
              n_i(Nwv)=zni(Nwv,IRH)
              
              rmin=size_bins_min(bin)*1.E-6_JPRB*rh_growth(IRH)
              rmax=size_bins_max(bin)*1.E-6_JPRB*rh_growth(IRH)
              
              DO dis=1,Ndis
                 r_0(dis)=r_0_inp(dis)*1.E-6_JPRB*rh_growth(IRH)
              ENDDO
           ENDDO
        endif


        !
        !----------Into the MIE code-----------
        !
        DO Nwv=1,nb_lambda

           m=CMPLX(n_r(Nwv),-n_i(Nwv))
          
           pi=4._JPRB*ATAN(1.)

           tot_num_concentration = 0.0_JPRB         
           sigma_sca = 0.0_JPRB
           sigma_ext = 0.0_JPRB
           gtot      = 0.0_JPRB
           omegatot  = 0.0_JPRB
           mass      = 0.0_JPRB
           volume    = 0.0_JPRB
           surface   = 0.0_JPRB
           sigma_abs = 0.0_JPRB
           radius    = 0.0_JPRB
           radius2   = 0.0_JPRB
           rgv       = 0.0_JPRB
           rgn       = 0.0_JPRB
           PP(:)     = 0.0_JPRB


           if(rmin.eq.rmax) Nbin_points=0
           
              if (verbose==3.and.Nwv==1) then
                ! This is a short test just for single spheres monodisperse
                ! we moved outside dbin to avoid be included in parallel do
                r = r_0(1)
                x = 2._JPRB*pi*r/lambda_out(Nwv)
                call mie_core_anbn(x, m, Nmax, a, b, Nsize)
                call mie_core_q(x, m, Nmax, a, b, Nsize, Q_ext, Q_sca, Q_abs, omega, g)
                test_single_sphere(1)=Q_ext
                test_single_sphere(2)=Q_sca
                test_single_sphere(3)=Q_abs
                test_single_sphere(4)=g
              endif

           DO dbin=0, Nbin_points 
              
              IF(rmin.eq.rmax) then
                 r=r_0(1)
                 deltar=1.0
                 numb=Ntot(1)
              else
                 r      = logn_r( rmin, rmax, dbin, Nbin_points)
                 deltar = logn_dr(rmin, rmax, Nbin_points)
                 numb   = logn_nt(Ntot, sigma_g, r, r_0, pi, Ndis, Nbin_points, dbin)
                 ! inside logn_nt we could have rho(dis) to allow for different densities
                 ! in the multi-mode distributions
                 ! mass=mass+4._JPRB/3._JPRB*pi*r**3*numberdis*deltar*rho(dis)        !--kg/m3*m3*1/m3

              endif

              mass=mass+logn_mass(rh_int, cutoff_radius, numb, deltar, rho, r, rh_growth(IRH), pi)

              !size parameter (assuming refr. index in vacuum)
              x=2._JPRB*pi*r/lambda_out(Nwv)

              !-- relevant for LR --------------------------------------
              tot_num_concentration=tot_num_concentration+numb*deltar
              radius  = radius +  r*numb*deltar
              radius2 = radius2+r*r*numb*deltar
              volume  = volume +4._JPRB/3._JPRB*pi*(r**3._JPRB)*numb*deltar
              surface = surface+4._JPRB   *pi*(r**2._JPRB)*numb*deltar
              rgv     = rgv+log(r)*(r**3._JPRB)*numb*deltar
              rgn     = rgn+log(r)       *numb*deltar
              !-----------------------------------------------------
              
              call mie_core_anbn(x, m , Nmax, a, b, 210000)
              call mie_core_q(x, m, Nmax, a, b, 210000, Q_ext, Q_sca, Q_abs, omega, g)
              call mie_core_phase(Nsize, Nmax, Nmumax, x, a, b, ph_ang, Q_sca, phase, pi)
            
              sigma_sca=sigma_sca+r**2*Q_sca*numb*deltar
              sigma_ext=sigma_ext+r**2*Q_ext*numb*deltar
              !-- for LR------------------------------------
              sigma_abs=sigma_abs+r**2*Q_abs*numb*deltar
              !---------------------------------------------------
              omegatot=omegatot+r**2*Q_ext*omega*numb*deltar
              gtot    =gtot+r**2*Q_sca*g*numb*deltar

              !-- for computation of LR:  phase function / lidar ratio -----------
              DO Nmu=1, Nmumax 
                 PP(Nmu)=PP(Nmu)+Q_sca*phase(Nmu)*numb*(r**2)*deltar
              ENDDO
           
           !-- end for LR--------------------------------------------------


           ENDDO   !---integration bin

           sigma_sca=pi*sigma_sca
           sigma_ext=pi*sigma_ext
           gtot=pi*gtot/sigma_sca
           omegatot=pi*omegatot/sigma_ext

           !-- for LR
           rgv=exp(rgv/(volume*3./4./pi))
           rgn=exp(rgn/tot_num_concentration)
           sigma_abs=pi*sigma_abs
           radius=radius  /tot_num_concentration
           radius2=radius2/tot_num_concentration
           
           final_g(Nwv)=gtot
           final_w(Nwv)=omegatot 
           final_a(Nwv)=sigma_ext/mass  !--ext coeff / mass or 80% RH mass
           !         print*,' mass =', mass  ,' g(aerosol) m-3(air); ext: ', sigma_ext, ' m-1'
           !         print*,' vol =', volume  ,' g(aerosol) m-3(air); ext: ', sigma_ext, ' m-1'

           DO Nmu=1, Nmumax
              ! here not it is not Nmu=0, Nmumax, but Nmu=1,Nmumax
              PP(Nmu)=pi/sigma_sca*PP(Nmu)
              !ph_ang(Nmu+1)=ACOS(Muk(Nmu))*180./pi
              ph_stored(Nwv,IRH,bin,Nmu)=PP(Nmu)
           ENDDO

           ! now PP(n) begins in PP(1)
           final_lratio(Nwv)=4.*pi/PP(1)/omegatot


           mass_stored(IRH,bin)=mass

           ext_stored(Nwv,IRH,bin)=final_a(Nwv)
           omg_stored(Nwv,IRH,bin)=final_w(Nwv)
           asy_stored(Nwv,IRH,bin)=final_g(Nwv)

           lidar_ratio_stored(Nwv,IRH,bin)=final_lratio(Nwv) 

           !progress status
           !counter=counter+1
           !if (verbose.eq.0.or.verbose.eq.2) then
           !   write(*,FMT="(A1,A,t42,F6.0,A)",ADVANCE="NO") achar(13), &
           !        & "                     Percent Complete:  ", (counter/NR)*100.0, "%"
           !   CALL sleep(2) !give a delay in sec to see the output 
           !end if


        ENDDO  !--loop on wavelength

     ENDDO  !----rh

  ENDDO   !-------size_bins



  if (verbose.eq.2) then
      print *, ''
      print *, '-FORTRAN subroutine -----------------------------------------'
      print *, ' znr (1,1)     ', znr(1,1)
      print *, ' zni (1,1)     ', zni(1,1)
      print *, ' ext (1,1,1)   ', ext_stored(1,1,1)
      print *, ' omg (1,1,1)   ', omg_stored(1,1,1)
      print *, ' asy (1,1,1)   ', asy_stored(1,1,1)
      print *, ' lidar (1,1,1) ', lidar_ratio_stored(1,1,1)
      print *, ' pfun (1,1,1,1)', ph_stored(1,1,1,1)
      print *, '-------------------------------------------------------------'
      print *, ' znr (2,1)     ', znr(2,1)
      print *, ' zni (2,1)     ', zni(2,1)
      print *, ' ext (2,1,1)   ', ext_stored(2,1,1)
      print *, ' omg (2,1,1)   ', omg_stored(2,1,1)
      print *, ' asy (2,1,1)   ', asy_stored(2,1,1)
      print *, ' lidar (2,1,1) ', lidar_ratio_stored(2,1,1)
      print *, ' pfun (2,1,1,1)', ph_stored(2,1,1,1)
      print *, '-------------------------------------------------------------'
      print *, ' znr (1,2)     ', znr(1,2)
      print *, ' zni (1,2)     ', zni(1,2)
      print *, ' ext (1,2,1)   ', ext_stored(1,2,1)
      print *, ' omg (1,2,1)   ', omg_stored(1,2,1)
      print *, ' asy (1,2,1)   ', asy_stored(1,2,1)
      print *, ' lidar (1,2,1) ', lidar_ratio_stored(1,2,1)
      print *, ' pfun (1,2,1,1)', ph_stored(1,2,1,1)
      print *, '-------------------------------------------------------------'
      print *, "*** Completed Fortran subroutine ***"
  endif

  contains


     real(8) function logn_r(rmin, rmax, dbin, Nbin_points)
         implicit none
         real(8)    , INTENT(IN)  :: rmin, rmax
         integer(4) , INTENT(IN)  :: Nbin_points, dbin

         logn_r=exp(log(rmin)+DBLE(dbin)/DBLE(Nbin_points)*(log(rmax)-log(rmin)))
     
     endfunction

     real(8) function logn_dr(rmin, rmax, Nbin_points)
         implicit none
         real(8)    , INTENT(IN)  :: rmin, rmax
         integer(4) , INTENT(IN)  :: Nbin_points

         logn_dr = 1._JPRB/DBLE(Nbin_points)*(log(rmax)-log(rmin))

     endfunction 

     real(8) function logn_nt(Ntot, sigma_g, r, r_0, pi, Ndis, Nbin_points, dbin)
         implicit none
         real(8)    , INTENT(IN)                     :: r, pi
         real(8)    , INTENT(IN), dimension(1:Ndis)  :: r_0, sigma_g, Ntot
         integer(4) , INTENT(IN)                     :: dbin, Ndis, Nbin_points
         real(8)                                     :: fbin, fexp
         if (dbin.eq.0.or.dbin.eq.Nbin_points) then
            fbin = 0.5
         else
            fbin = 1.0
         endif
         logn_nt=0.0
         do dis=1,Ndis
            fexp = exp(-0.5_JPRB*(log(r/r_0(dis))/log(sigma_g(dis)))**2)
            logn_nt=logn_nt + fbin*fexp*Ntot(dis)/SQRT(2._JPRB*pi)/log(sigma_g(dis))
         enddo
     endfunction


     real(8) function logn_mass(rh_int, cutoff_radius, numb, deltar, rho, r, rh_growth, pi)
         implicit none
         real(8),    INTENT(IN)  :: pi,cutoff_radius, numb, deltar, rho, r, rh_growth
         integer(4), INTENT(IN)  :: rh_int 
         real(8)                 :: rnew

         if (r.le.(cutoff_radius*1.e-6)) then
              if (rh_int .eq. 1) then
                 rnew = r
              else
                 rnew = r/rh_growth
              endif
              logn_mass = 4._JPRB/3._JPRB*pi*rnew**3*numb*deltar*rho        !--kg/m3*m3*1/m3
         else
              logn_mass = 0._JPRB 
         endif

     endfunction

  end subroutine



  subroutine mie_core_phase(Nsize, Nmax, Nmumax,x, a,b, ph_ang, Q_sca, phase, pi)
    
    implicit none 
    INTEGER(4)    , INTENT(IN)     :: Nsize, Nmax
    REAL(8)       , INTENT(IN)     :: pi, Q_sca, x
    INTEGER(4)    , INTENT(IN)     :: Nmumax
    COMPLEX(KIND=JPRB), INTENT(IN) :: a(1:Nsize), b(1:Nsize) 
    REAL(8)       , INTENT(INOUT)  :: ph_ang(Nmumax)
    REAL(8)       , INTENT(INOUT)  :: phase(1:Nmumax)
    REAL(8)                        :: Muk(1:Nmumax)
    REAL(8)                        :: ang, mu, n_float
    REAL(8)                        :: ss, tt, rho1, rho2
    REAL(8)                        :: wist(1:Nsize), wisp(0:Nsize)
    COMPLEX(KIND=JPRB)             :: s1, s2, nn
    INTEGER(4)                     :: Nmu, n  



    DO Nmu=1, Nmumax
       ! before we had Num=0, Nmumax because the eq. for the
       !angles uses Nmu=0, but now when we pass an array
       ! we can keep everything from 1. 
       !ang=pi-pi*DBLE(Nmu)/DBLE(Nmumax)
       
       ang = pi*ph_ang(Nmu)/180.0
       mu=COS(ang)
       Muk(Nmu)=mu

       !-- Wiscombe method: it is faster than usual method to calculate
       !   S1 and S2 matrix by calculating ss1 = (S1+S2) and ss2=(S1-S2)
       !   but eventually the original method is used. However current compilers
       !   could be equality efficent in both methods.
       !   Here we could have a more structured calculation by separating the 
       !   wist and wisp functions to actually evaluate this.
       wisp(0)=0.
       wisp(1)=1.
       DO n=1,Nmax-1
          n_float   = DBLE(n)
          ss        = mu*wisp(n)
          tt        = ss-wisp(n-1)
          wist(n)   = n_float*tt-wisp(n-1)
          wisp(n+1) = ss+((n_float+1.)/n_float)*tt
       ENDDO
       n       = Nmax
       n_float = DBLE(n)
       ss      = mu*wisp(n)
       tt      = ss-wisp(n-1)
       wist(n) = n_float*tt-wisp(n-1)

       s1=0.0
       s2=0.0
       
       DO n=1,Nmax
          nn=CMPLX(DBLE(n),0.0)
          s1=s1+(2.*nn+1.)/nn/(nn+1.)*(a(n)*wisp(n)+b(n)*wist(n))
          s2=s2+(2.*nn+1.)/nn/(nn+1.)*(b(n)*wisp(n)+a(n)*wist(n))
       ENDDO

       s1=s1*DCONJG(s1)
       s2=s2*DCONJG(s2)

       rho1=4./x**2 * DBLE(s1)/Q_sca
       rho2=4./x**2 * DBLE(s2)/Q_sca
       phase(Nmu)=0.5*(rho1+rho2)
       ! here the output should be the Fmatrix
    ENDDO

  endsubroutine 

  subroutine mie_core_q(x, m, Nmax, a, b, Nsize, Qext, Qsca, Qabs, omega, g)

    use parkind1, only : JPIM, JPRB 
    implicit none
    COMPLEX(KIND=JPRB), INTENT(IN) :: m    
    INTEGER(4)        , INTENT(IN) :: Nmax, Nsize
    COMPLEX(KIND=JPRB), INTENT(IN) :: a(1:Nsize), b(1:Nsize) 
    REAL(8)           , INTENT(OUT):: Qabs,  Qext, Qsca, g, omega, x 
    INTEGER(4)                     :: n
    REAL(8)                        :: n_float 
    
    Qext=0.0_JPRB
    Qsca=0.0_JPRB
    g=0.0_JPRB
    DO n=Nmax-1,1,-1
       n_float=DBLE(n)
       Qext=Qext+ (2._JPRB*n_float+1._JPRB) * DBLE( a(n)+b(n) )
       Qsca=Qsca+ (2._JPRB*n_float+1._JPRB) * &
            &           DBLE( a(n)*DCONJG(a(n)) + b(n)*DCONJG(b(n)) )
       g=g + n_float*(n_float+2._JPRB)/(n_float+1._JPRB) * &
            &           DBLE( a(n)*DCONJG(a(n+1))+b(n)*DCONJG(b(n+1)) )  + &
            &           (2._JPRB*n_float+1._JPRB)/n_float/(n_float+1._JPRB) * DBLE(a(n)*DCONJG(b(n)))
    ENDDO
    Qext=2._JPRB/x**2 * Qext
    Qsca=2._JPRB/x**2 * Qsca
    !-- for LR----------------------
    Qabs=Qext-Qsca
    if (AIMAG(m).EQ.0.0) Qabs=0.0
    !-------------------------------------
    omega=Qsca/Qext
    g=g*4./x**2/Qsca

  endsubroutine 
 

  subroutine mie_core_anbn(x, m, Nmax, a, b, Nsize)
    
    use parkind1, only : JPIM, JPRB 
    implicit none

    REAL(8)           , INTENT(IN)    :: x
    COMPLEX(KIND=JPRB), INTENT(IN)    :: m    !--- refractive index m=n_r-i*n_i
    COMPLEX(KIND=JPRB), INTENT(INOUT) :: a(1:Nsize), b(1:Nsize)
    INTEGER(4)        , INTENT(OUT)   :: Nmax
    INTEGER(4)        , INTENT(IN)    :: Nsize
    INTEGER(4)         :: Nstart    !--- number of iterations for the Wiscombe algorithm
    INTEGER(4)         :: n
    COMPLEX(KIND=JPRB) :: k2, k3, z1, z2
    COMPLEX(KIND=JPRB) :: u1,u5,u6,u8,nn,I
    COMPLEX(KIND=JPRB) :: ksiz2(-1:Nsize), psiz2(1:Nsize)
    COMPLEX(KIND=JPRB) :: nu1z1(1:Nsize+100), nu1z2(1:Nsize+100)
    COMPLEX(KIND=JPRB) :: nu3z2(0:Nsize)


    I=CMPLX(0.,1.)
    k2=m
    k3=CMPLX(1.0,0.0)
    z2=CMPLX(x,0.0)
    z1=m*z2





    IF (0.0.LE.x.AND.x.LE.8.) THEN
       Nmax=INT(x+4*x**(1./3.)+1.)+2
    ELSEIF (8..LT.x.AND.x.LT.4200.) THEN
       Nmax=INT(x+4.05*x**(1./3.)+2.)+1
    ELSEIF (4200..LE.x.AND.x.LE.20000.) THEN
       Nmax=INT(x+4*x**(1./3.)+2.)+1
    ELSE
       WRITE(*,*) 'x out of bound, x=', x
       STOP
    ENDIF

    !this has considerable impact on the accuracy of results
    !especially the scattering cross section for large size parameters.
    !and low absorption. Using Nmax+10 works for size parameters up to ~50
    !Nmax+100 works for size parameters up to about 100 (for m_i>~1e-3). 
    !I haven't check for larger values.
    !Nstart controls the number of terms to use in the downward recurrence loop
    !to compute the derivative of the spherical Riccati-Bessel functions necessary in the
    !computation of the a(n) and b(n). Equations explaining the following loops are
    !in the 6SV documentation for the MIE subroutine (part 2, pag 108-125). nu1z1 and
    !nu1z2 correspond respectively to Dn(mx) and Dn(x) while nu3z2 corresponds to Gm(x)
    !Further discussion of this downward recurrence is given in Hong Du 2004 Appl. Optics.

    Nstart=Nmax+100

    !-----------loop for nu1z1, nu1z2

    nu1z1(Nstart)=CMPLX(0.0,0.0)
    nu1z2(Nstart)=CMPLX(0.0,0.0)
    DO n=Nstart-1, 1 , -1
       nn=CMPLX(DBLE(n),0.0)
       nu1z1(n)=(nn+1.)/z1 - 1./( (nn+1.)/z1 + nu1z1(n+1) )
       nu1z2(n)=(nn+1.)/z2 - 1./( (nn+1.)/z2 + nu1z2(n+1) )
    ENDDO

    !------------loop for nu3z2

    nu3z2(0)=-I     

    DO n=1, Nmax
       nn=CMPLX(DBLE(n),0.0)
       nu3z2(n)=-nn/z2 + 1./ (nn/z2 - nu3z2(n-1) )
    ENDDO

    !-----------loop for psiz2 and ksiz2 (z2)
    ksiz2(-1)=COS(DBLE(z2))-I*SIN(DBLE(z2))
    ksiz2(0)=SIN(DBLE(z2))+I*COS(DBLE(z2))
    DO n=1,Nmax
       nn=CMPLX(DBLE(n),0.0)
       ksiz2(n)=(2.*nn-1.)/z2 * ksiz2(n-1) - ksiz2(n-2)
       psiz2(n)=CMPLX(DBLE(ksiz2(n)),0.0)
    ENDDO

    !-----------loop for a(n) and b(n)

    DO n=1, Nmax
       u1=k3*nu1z1(n) - k2*nu1z2(n)
       u5=k3*nu1z1(n) - k2*nu3z2(n)
       u6=k2*nu1z1(n) - k3*nu1z2(n)
       u8=k2*nu1z1(n) - k3*nu3z2(n)
       a(n)=psiz2(n)/ksiz2(n) * u1/u5
       b(n)=psiz2(n)/ksiz2(n) * u6/u8
    ENDDO

 end subroutine
end module
