




module mie_Boucher_Bozzo

  !====================================================================================
  ! author:  Ramiro Checa-Garcia                                                      !
  ! email:   ramiro.checa-garcia@ecmwf.int                                            !
  ! history:                                                                          !
  !          Sep-2022 -------> First implmentation and testing                        !
  !                                                                                   !
  ! sources: it is using a subroutine from O.Boucher and A.Bozzo with                 !
  !          several modifications                                                    !
  !                                                                                   !
  ! info:    The original code of Boucher and Bozzo has been changed to a module      !
  !          and it is created a wrapper of the fortran subroutine using              !
  !          iso_c_bindings which allows the access of the mie_aerosols subroutine    !
  !          from different languages.                                                !
  !                                                                                   !
  !                                                                                   !
  ! Also rather than use the KIND defined in parkind for the in and out               !
  ! variables it is real(8) and int(4)                                                !
  !                                                                                   !
  !                                                                                   !
  ! TODO:                                                                             !
  !    - [DONE] evaluate the impact of real(8) and int(4) in the code by comparing    !
  !      with the original results. The difference between both codes gives 0.0       !
  !      in all the variables.                                                        !
  !                                                                                   !
  !====================================================================================

  use iso_c_binding
  use iso_fortran_env

  implicit none
  
  contains

    subroutine mie_aerosols( &
          & Ninp, Nmumax, nb_lambda, lambda_out, Ndis, size_bins, &   ! 6
          & size_bins_min, size_bins_max,                         &   ! 2
          & sigma_g, r_0_inp, Ntot, rho,          ri_nrh,         &   ! 5
          & rh_int, rh_tab, rh_growth, cutoff_radius, verbose,    &   ! 5
          & zlambtab0, z_n_r_tab, z_n_i_tab,                      &   ! 3 => 22 in
          & znr,zni,                                              &   ! 2
          & ext_stored, omg_stored, asy_stored,                   &   ! 3 + 4 below => 9
          & lidar_ratio_stored, mass_stored, ph_stored, ph_ang)   bind(c, name="mie_aerosols")


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
        real(KIND=c_double) , intent(IN) :: zlambtab0(NInp)
        real(KIND=c_double) , intent(IN) :: z_n_r_tab(NInp,ri_nrh), z_n_i_tab(NInp,ri_nrh)

        ! output quantities (per specie size bin ; per wavelength ; per rh(if present))
        real(KIND=c_double) , intent(OUT)  :: znr(nb_lambda,rh_int),zni(nb_lambda,rh_int)
        real(KIND=c_double) , intent(OUT)  :: ext_stored(nb_lambda,rh_int,size_bins)
        real(KIND=c_double) , intent(OUT)  :: omg_stored(nb_lambda,rh_int,size_bins)
        real(KIND=c_double) , intent(OUT)  :: asy_stored(nb_lambda,rh_int,size_bins)
        real(KIND=c_double) , intent(OUT)  :: mass_stored(rh_int,size_bins)
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
            print *, ' ri_lamb_tab', zlambtab0(:)
            print *, ' rho        ', rho
            print *, ' size_bins  ', size_bins
            print *, "============================================================="
        endif


        call main_mie_aerosols( &
            & Ninp, Nmumax, nb_lambda, lambda_out, Ndis, size_bins, &   ! 6
            & size_bins_min, size_bins_max,                         &   ! 2
            & sigma_g, r_0_inp, Ntot, rho,          ri_nrh,         &   ! 5
            & rh_int, rh_tab, rh_growth, cutoff_radius, verbose,    &   ! 5
            & zlambtab0, z_n_r_tab, z_n_i_tab,                      &   ! 3 => 21 in
            & znr,zni,                                              &   ! 2
            & ext_stored, omg_stored, asy_stored,                   &   ! 3
            & lidar_ratio_stored, mass_stored, ph_stored, ph_ang)       ! 4 => 9  out)



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
          & zlambtab0, z_n_r_tab, z_n_i_tab,                      &   ! 3 => 21 in
          & znr,zni,                                              &   ! 2
          & ext_stored, omg_stored, asy_stored,                   &   ! 3
          & lidar_ratio_stored, mass_stored, ph_stored, ph_ang)       ! 4 => 9  out


          USE PARKIND1 , ONLY : JPIM, JPRB

          IMPLICIT NONE


!-------Mie computations for a size distribution of homogeneous spheres---------------------------
!
!-------Ref : Toon and Ackerman, Applied Optics, 1981
!             Stephens, CSIRO, 1979
!       Attention : surdimensionement des tableaux
!       to be compiled with double precision option
!
! AUTHOR: Olivier Boucher
!
! modifs: Alessio Bozzo 2017
!         Ramiro Checa-Garcia 2022 => code changes, added wrapper, refrac index info
!                                     is passed by arrays not by reading files.
!                                     changes also about KIND fortran specification.
!------------------------------------------------------------------------------------------------
!
! Input/configuration parameters:
!
!  NInp      : number of reference wavelengths in the refractive index file
!  nb_lambda : number of output  wavelengths
!  Ndis      : number of distributions considered (mono modal or bi-modal)
!  size_bins : number of size bins considered for the specie
!
!-------limits min/max of each size bin. Example
! size_bins=3
!   size_bins_min = (/ 0.03, 0.55, 0.90 /)
!   size_bins_max = (/ 0.55, 0.90, 20.0 /)
! size_bins=9
!   size_bins_min = (/ 0.03,0.05,0.1,0.2,0.5,1.,2.,5.,10. /)
!   size_bins_max = (/ 0.05,0.1,0.2,0.5,1.,2.,5.,10.,20. /)
!
!----------------------------------------------------------------------------------------------
!-------SIZE distribution properties---------------- => it considers a log-normal distribution
!  sigma_g : geometric standard deviation
!  r_0     : geometric number mean radius (um)/modal radius
!  Ntot    : total concentration in m-3
!
!  rho       : density in kgm-3
!----------------------------------------------------------------------------------------------
!
! OUTPUT
!
! znr(nb_lambda,rh_int)
! zni(nb_lambda,rh_int)
! ext_stored(nb_lambda,rh_int,size_bins)
! omg_stored(nb_lambda,rh_int,size_bins)
! asy_stored(nb_lambda,rh_int,size_bins)
! lidar_ratio_stored(nb_lambda,rh_int,size_bins)
! mass_stored(rh_int,size_bins)
! ph_stored(nb_lambda,rh_int,size_bins,0:Nmumax)
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
    REAL(8)   , INTENT(IN) :: cutoff_radius  !mass cut-off radius in microns
    REAL(8)   , INTENT(IN) :: zlambtab0(NInp)
    REAL(8)   , INTENT(IN) :: z_n_r_tab(NInp,ri_nrh), z_n_i_tab(NInp,ri_nrh)
    !CHARACTER(LEN=200), INTENT(IN) :: RI_file

    ! intent(in) 19 objects


    ! output quantities (per specie size bin, per wavelength, per rh(if present))
    REAL(8),INTENT(OUT)   :: znr(nb_lambda,rh_int),zni(nb_lambda,rh_int)
    REAL(8),INTENT(OUT)   :: ext_stored(nb_lambda,rh_int,size_bins)
    REAL(8),INTENT(OUT)   :: omg_stored(nb_lambda,rh_int,size_bins)
    REAL(8),INTENT(OUT)   :: asy_stored(nb_lambda,rh_int,size_bins)
    REAL(8),INTENT(OUT)   :: mass_stored(rh_int,size_bins)
    REAL(8),INTENT(OUT)   :: lidar_ratio_stored(nb_lambda,rh_int,size_bins)
    REAL(8),INTENT(OUT)   :: ph_stored(nb_lambda,rh_int,size_bins,Nmumax)
    REAL(8),INTENT(INOUT) :: ph_ang(Nmumax)

    INTEGER(4) :: Nbin_points,Nlimit
    REAL(8) :: rmin,rmax,r_0(Ndis)


    COMPLEX(KIND=JPRB) :: m           !--- refractive index m=n_r-i*n_i
    INTEGER(4) :: Nmax,Nstart !--- number of iterations for the Wiscombe algorithm
    COMPLEX(KIND=JPRB) :: k2, k3, z1, z2
    COMPLEX(KIND=JPRB) :: u1,u5,u6,u8,nn,I
    COMPLEX(KIND=JPRB) :: a(1:210000), b(1:210000)
    COMPLEX(KIND=JPRB) :: ksiz2(-1:210000), psiz2(1:210000)
    COMPLEX(KIND=JPRB) :: nu1z1(1:210100), nu1z2(1:210100)
    COMPLEX(KIND=JPRB) :: nu3z2(0:210000)

    !-- additional variables and arrays for computation of the
    !-- phase function and lidar ratio
    REAL(8) :: volume, surface, rgv, rgn, tot_num_concentration, radius, radius2, mass, pfun_angle,ang
    REAL(8) :: sigma_abs, ss, tt, rho1, rho2, Q_abs
    REAL(8) :: wist(1:210000), wisp(0:210000)
    COMPLEX(KIND=JPRB) :: s1, s2
    INTEGER(4) :: Nmu!, Nmumax
    !PARAMETER (Nmumax=1200)
    REAL(8) :: PP(0:Nmumax), phase(0:Nmumax), Muk(0:Nmumax)


    REAL(8) :: Q_ext, Q_sca, g, omega   !--parameters for radius r
    REAL(8) :: x  !--size parameter
    REAL(8) :: r  !--radius
    REAL(8) :: sigma_sca, sigma_ext, omegatot,  gtot !--averaged parameters
    REAL(8) :: number, deltar

    !-- interpolation of refractive indices on output wavelengths

    REAL(8)    :: zlambtab(NInp)
    !z_n_r_tab(NInp,ri_nrh), z_n_i_tab(NInp,ri_nrh)
    !REAL(8)    :: znr(nb_lambda,rh_int)   , zni(nb_lambda,rh_int)
    REAL(8)    :: n_r(nb_lambda)   , n_i(nb_lambda)
    REAL(8)    :: zint


    ! intermediate values
    REAL(8) :: final_a(nb_lambda), final_g(nb_lambda), final_w(nb_lambda)
    REAL(8) :: final_lratio(nb_lambda) !, final_PH(nb_lambda,0:Nmumax)

    !general looping indices and utility parameters
    INTEGER(4) :: jwv, jwl, n, nb, Nwv, IRH, dis, bin,dbin,dum
    REAL(8)    :: pi, n_float,numberdis, NR, counter


    !-------------------OUTPUT FILES---------------------------------------------------------------
    !OPEN(unit=50,form='formatted',file='size_distribution_ss')

    !----------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------
    !-- START OF COMPUTATIONS
   !   ---------------------

    !total number of iterations (bins*RH*wl)
    NR=real(nb_lambda*rh_int*size_bins)
    counter=0


    Nbin_points=999


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
        print *, ' ri_lamb_tab', zlambtab0(1), zlambtab0(2) 
        print *, ' rho        ', rho
        print *, ' size_bins  ', size_bins
        print *, "==============================================================="
    endif


    !open(10,file=RI_file)
    !read(10,*) cblabla
    !read(10,*) cblabla

    if(verbose.eq.1) then
    print*,'Input refr. idx on N wavelengths: ',NInp
    endif

    !if (ri_nrh.gt.1) then
    !   DO irh=1,ri_nrh
    !      DO jwv=1,NInp
    !         read(10,*) zlambtab(jwv), z_n_r_tab(jwv,irh), z_n_i_tab(jwv,irh), dum, dum
    !         if(verbose.eq.1) then
    !            PRINT 9002,jwv,zlambtab(jwv), z_n_r_tab(jwv,irh),z_n_i_tab(jwv,irh)
    !         end if
    !      ENDDO
    !   ENDDO
    !else
    !   DO jwv=1,NInp
    !      read(10,*) zlambtab(jwv), z_n_r_tab(jwv,1), z_n_i_tab(jwv,1)
    !      if(verbose.eq.1) then
    !         PRINT 9002,jwv,zlambtab(jwv), z_n_r_tab(jwv,1),z_n_i_tab(jwv,1)
    !      end if
    !   ENDDO
    !end if
    !9002 format(1x,I3,E10.3,2F10.6)

    !check if lambda in input ri file are in m
    do jwv=1,NInp
        if (zlambtab0(jwv) .gt. 100e-6) then
            !print *, ' index ', jwv, ' value :', zlambtab0(jwv)
            zlambtab(jwv)=zlambtab0(jwv)*1.e-6
            !print *,'Input wavelength for refractive index not in m! converting to m.'
            !print *,'CHECK Input wavelength in refractive index file!'
        else
            zlambtab(jwv)=zlambtab0(jwv)
        end if
    enddo

if(verbose.eq.1) then
   PRINT *,'Wavelengths at which Mie computations are actually performed'
   PRINT *,''

   DO nb=1, nb_lambda
     print 9999,nb,lambda_out(nb)
9999  format(1x,I3,E12.5)
   ENDDO
end if
! external subroutine performing the interpolation of the input refr index
! onto the working wavelengths. input wavelength always given in m


!print *, NInp, rh_int

!DO IRH=1,rh_int
!   print *, z_n_r_tab(1:10,IRH)
!ENDDO

DO IRH=1,rh_int
   if (NInp.gt.1) then
      call interp_ri(NInp,zlambtab, z_n_r_tab(:,IRH), z_n_i_tab(:,IRH), nb_lambda,lambda_out,znr(:,IRH),zni(:,IRH),verbose)
   else
      znr(:,IRH)=z_n_r_tab(1,IRH)
      zni(:,IRH)=z_n_i_tab(1,IRH)
   endif
ENDDO

!print *, ''
!print *, 'INTERPOLATED'
!print *, shape(znr)
!DO IRH=1,rh_int
!   print *, IRH, znr(:,IRH)
!   print *, IRH, lambda_out
!ENDDO
!print *,'---------------'

!-- loop on Bins
if(verbose.eq.1) then
   print *,'Loops on bins/RH/wavelengths start'
end if


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
         !print *,'Enter the Mie code with NWv=',NWv
         if(verbose.eq.1) then
            print *,'bin, IRH, Nwv= ', bin, IRH, Nwv, '(', size_bins, rh_int, nb_lambda, ')'
         end if

         m=CMPLX(n_r(Nwv),-n_i(Nwv))
        
         pi=4._JPRB*ATAN(1.)


         I=CMPLX(0.,1.)

         sigma_sca=0.0_JPRB
         sigma_ext=0.0_JPRB
         gtot=0.0_JPRB
         omegatot=0.0_JPRB
         mass = 0.0_JPRB

!-- following quantities relevanty for Lidar ratio computations --------
         volume=0.0_JPRB
         surface=0.0_JPRB
         sigma_abs=0.0_JPRB
         tot_num_concentration=0.0_JPRB
         radius=0.0_JPRB
         radius2=0.0_JPRB
         rgv=0.0_JPRB
         rgn=0.0_JPRB
         PP(:)=0.0
!         DO Nmu=0,Nmumax
!            PP(Nmu)=0.0
!         ENDDO
!-----------------------
!
!---integration over size distribution(s)
!---loop on size points within a size bin of the specie
!---if in input rmin=rmax the computations will be for a
!---single particle

         if(rmin.eq.rmax) Nbin_points=0
         
         DO dbin=0, Nbin_points 
            
            if(rmin.eq.rmax) then
               r=r_0(1)
               deltar=1.0
               number=Ntot(1)
            else
               r=exp(log(rmin)+DBLE(dbin)/DBLE(Nbin_points)*(log(rmax)-log(rmin)))
               deltar=1._JPRB/DBLE(Nbin_points)*(log(rmax)-log(rmin))
               !          print*,r
               
               number=0._JPRB
               DO dis=1, Ndis
                  numberdis=Ntot(dis)/SQRT(2._JPRB*pi)/log(sigma_g(dis))* &
                       &           exp(-0.5_JPRB*(log(r/r_0(dis))/log(sigma_g(dis)))**2)
                  IF (dbin.EQ.0.OR.dbin.EQ.Nbin_points) numberdis=numberdis/2.0_JPRB
                  number=number+numberdis
!here to allow for different densities in the multi-mode distributions
!             mass=mass+4._JPRB/3._JPRB*pi*r**3*numberdis*deltar*rho(dis)        !--kg/m3*m3*1/m3
               ENDDO

            endif

!       if (IRH.eq.9) then
!          !write the size distribution at 80%
!          WRITE(50,*) r*1.e6,number*deltar
!       endif
       !cut-off radius for the computation of mass
       !this is used in the OPAC aerosol types
            if (r.le.(cutoff_radius*1.e-6)) then

               if (rh_int .eq. 1) then
                  mass=mass+4._JPRB/3._JPRB*pi*r**3*number*deltar*rho        !--kg/m3*m3*1/m3
               else
                  !mass relative to fixed RH (i.e. dry mass in case rh_growth=1 at 0% RH)
                  mass=mass+4._JPRB/3._JPRB*pi*((r/rh_growth(IRH))**3)*number*  &
                       &                    deltar*rho          
               endif
          
            endif

            !size parameter (assuming refr. index in vacuum)
            x=2._JPRB*pi*r/lambda_out(Nwv)


!-- relevant for LR --------------------------------------
            tot_num_concentration=tot_num_concentration+number*deltar
            !print*,tot_num_concentration
            radius =radius +  r*number*deltar
            radius2=radius2+r*r*number*deltar
            volume =volume +4._JPRB/3._JPRB*pi*(r**3._JPRB)*number*deltar
            surface=surface+4._JPRB   *pi*(r**2._JPRB)*number*deltar
            rgv=rgv+log(r)*(r**3._JPRB)*number*deltar
            rgn=rgn+log(r)       *number*deltar
!-----------------------------------------------------
            
            k2=m
            k3=CMPLX(1.0,0.0)
            !k3=CMPLX(1.00029,0.0)

            z2=CMPLX(x,0.0)
            z1=m*z2

            IF (0.0.LE.x.AND.x.LE.8.) THEN
               Nmax=INT(x+4*x**(1./3.)+1.)+2
            ELSEIF (8..LT.x.AND.x.LT.4200.) THEN
               Nmax=INT(x+4.05*x**(1./3.)+2.)+1
            ELSEIF (4200..LE.x.AND.x.LE.20000.) THEN
               Nmax=INT(x+4*x**(1./3.)+2.)+1
            ELSE
               !        print*,'x out of bound, x=', x
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

!-----------------final loop--------------
            Q_ext=0.0_JPRB
            Q_sca=0.0_JPRB
            g=0.0_JPRB
            DO n=Nmax-1,1,-1
               n_float=DBLE(n)
               Q_ext=Q_ext+ (2._JPRB*n_float+1._JPRB) * DBLE( a(n)+b(n) )
               Q_sca=Q_sca+ (2._JPRB*n_float+1._JPRB) * &
                    &           DBLE( a(n)*DCONJG(a(n)) + b(n)*DCONJG(b(n)) )
               g=g + n_float*(n_float+2._JPRB)/(n_float+1._JPRB) * &
                    &           DBLE( a(n)*DCONJG(a(n+1))+b(n)*DCONJG(b(n+1)) )  + &
                    &           (2._JPRB*n_float+1._JPRB)/n_float/(n_float+1._JPRB) * DBLE(a(n)*DCONJG(b(n)))
            ENDDO
            Q_ext=2._JPRB/x**2 * Q_ext
            Q_sca=2._JPRB/x**2 * Q_sca
!-- for LR----------------------
            Q_abs=Q_ext-Q_sca
            if (AIMAG(m).EQ.0.0) Q_abs=0.0
!-------------------------------------
            omega=Q_sca/Q_ext
            g=g*4./x**2/Q_sca

            sigma_sca=sigma_sca+r**2*Q_sca*number*deltar
            sigma_ext=sigma_ext+r**2*Q_ext*number*deltar
!-- for LR------------------------------------
            sigma_abs=sigma_abs+r**2*Q_abs*number*deltar
!---------------------------------------------------
            omegatot=omegatot+r**2*Q_ext*omega*number*deltar
            gtot    =gtot+r**2*Q_sca*g*number*deltar
            !      print*,2.*r,pi*r**2*Q_ext
            !      print*,2.*r,g
            !      print*,2.*r,pi*r**2*Q_sca



!-- for computation of LR:  phase function / lidar ratio -----------
!   Ramiro Checa-Garcia. Here this has been changed to have angles as
!                        input.
            DO Nmu=1, Nmumax
               ! before we had Num=0, Nmumax because the eq. for the
               ! angles uses Nmu=0, but now when we pass an array
               ! we can keep everything frim 1. 
         
               !ang=pi-pi*DBLE(Nmu)/DBLE(Nmumax)
               
               ang = pi*ph_ang(Nmu)/180.0
               !print *, Nmu, ph_ang(Nmu), ang
               !print *, Nmu, ph_ang(Nmu+1), angname::AbstractString,
               mu=COS(ang)
               Muk(Nmu)=mu

!-- Algorithme de Wiscombe
               wisp(0)=0.
               wisp(1)=1.
               DO n=1,Nmax-1
                  n_float=DBLE(n)
                  ss=mu*wisp(n)
                  tt=ss-wisp(n-1)
                  wist(n)=n_float*tt-wisp(n-1)
                  wisp(n+1)=ss+((n_float+1.)/n_float)*tt
               ENDDO
               n=Nmax
               n_float=DBLE(n)
               ss=mu*wisp(n)
               tt=ss-wisp(n-1)
               wist(n)=n_float*tt-wisp(n-1)

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

               PP(Nmu)=PP(Nmu)+Q_sca*phase(Nmu)*number*(r**2)*deltar
            ENDDO
!-- end for LR--------------------------------------------------


         ENDDO   !---integration bin
!------------------------------------------------------------------

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
         counter=counter+1
         if (verbose.eq.0.or.verbose.eq.2) then
            write(*,FMT="(A1,A,t42,F6.0,A)",ADVANCE="NO") achar(13), &
                 & "                     Percent Complete:  ", (counter/NR)*100.0, "%"
            CALL sleep(2) !give a delay in sec to see the output 
         end if


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


end subroutine


end module
