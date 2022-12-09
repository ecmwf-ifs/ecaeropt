module tmat_parameters

    !implicit none

    integer, parameter :: NPN1=100, NPNG1=300
    integer, parameter :: NPNG2=2*NPNG1, NPN2=2*NPN1
    integer, parameter :: NPL=NPN2+1, NPN3=NPN1+1
    integer, parameter :: NPN4=80, NPN5=2*NPN4, NPN6=NPN4+1, NPL1=NPN5+1

end module tmat_parameters

module tmat_ctt
  
    use tmat_parameters

    real*8 ::    QR(NPN2,NPN2),QI(NPN2,NPN2)
    real*8 ::    RGQR(NPN2,NPN2),RGQI(NPN2,NPN2)

end module tmat_ctt


module tmat_tmat

    use tmat_parameters
    REAL*8 ::  TR11(NPN6,NPN4,NPN4),TR12(NPN6,NPN4,NPN4)
    REAL*8 ::  TR21(NPN6,NPN4,NPN4),TR22(NPN6,NPN4,NPN4)
    REAL*8 ::  IT11(NPN6,NPN4,NPN4),IT12(NPN6,NPN4,NPN4)
    REAL*8 ::  IT21(NPN6,NPN4,NPN4),IT22(NPN6,NPN4,NPN4)

end module tmat_tmat

module tmat_ct

    use tmat_parameters

    real*8 ::    TR1(NPN2,NPN2)
    real*8 ::    TI1(NPN2,NPN2)

end module tmat_ct
