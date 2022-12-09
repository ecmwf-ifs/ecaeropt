
program tmatrix_example

      use tmatrixlib

      IMPLICIT NONE

      integer :: NDISTR, NPNAX, NKMAX, NPNA, NDGS, NP
      REAL*8  :: RAT, AXMAX, B, GAM, EPS, LAM, MRR, MRI, DDELT

 
      RAT=0.5 
      NDISTR=3
      AXMAX=1D0
      NPNAX=2
      B=0.1D0
      GAM=0.5D0
      NKMAX=5
      EPS=2D0
      NP=-1
      LAM=0.5D0 
      MRR=1.53
      MRI=0.008D0 
      DDELT=0.001D0 
      NPNA=19 
      NDGS=2


      call tmatrix( RAT, NDISTR, AXMAX, NPNAX, B, GAM, NKMAX, EPS, NP, &
                  & LAM, MRR, MRI, DDELT, NPNA, NDGS)



end program

