
module tmat_CGcoeff

  use iso_fortran_env
  use tmat_parameters
  
  contains

!****************************************************************
!
!   CALCULATION OF THE QUANTITIES F(N+1)=0.5*LN(N!)
!   0.LE.N.LE.899
 
      SUBROUTINE FACT(F)
      REAL*8 F(900)
      !COMMON /FAC/ F
      F(1)=0D0
      F(2)=0D0
      DO 2 I=3,900
         I1=I-1
         F(I)=F(I1)+0.5D0*DLOG(DFLOAT(I1))
    2 CONTINUE
      RETURN
      END
 
!************************************************************
!
!   CALCULATION OF THE ARRAY SSIGN(N+1)=SIGN(N)
!   0.LE.N.LE.899
 
      SUBROUTINE SIGNUM(SSIGN)
      REAL*8 SSIGN(900)
      !COMMON /SS/ SSIGN
      SSIGN(1)=1D0
      DO 2 N=2,899 
         SSIGN(N)=-SSIGN(N-1)
    2 CONTINUE
      RETURN
      END
 
!******************************************************************
!
!   CALCULATION OF CLEBSCH-GORDAN COEFFICIENTS
!   (N,M:N1,M1/NN,MM)
!   FOR GIVEN N AND N1. M1=MM-M, INDEX MM IS FOUND FROM M AS
!   MM=M*K1+K2
!
!   INPUT PARAMETERS :  N,N1,NMAX,K1,K2
!                               N.LE.NMAX
!                               N.GE.1
!                               N1.GE.0
!                               N1.LE.N+NMAX
!   OUTPUT PARAMETERS : GG(M+NPN6,NN+1) - ARRAY OF THE CORRESPONDING
!                                       COEFFICIENTS
!                               /M/.LE.N
!                               /M1/=/M*(K1-1)+K2/.LE.N1
!                               NN.LE.MIN(N+N1,NMAX)
!                               NN.GE.MAX(/MM/,/N-N1/)
!   IF K1=1 AND K2=0, THEN 0.LE.M.LE.N
 
 
      SUBROUTINE CCG(N,N1,NMAX,K1,K2,GG)

      use tmat_parameters
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 GG(NPL1,NPN6),CD(0:NPN5),CU(0:NPN5)
      REAL*8 FAC(900), SSIGN(900)

      IF(NMAX.LE.NPN4.AND.0.LE.N1.AND.N1.LE.NMAX+N.AND.N.GE.1.AND.N.LE.NMAX) GO TO 1
      PRINT 5001
      STOP
 5001 FORMAT(' ERROR IN SUBROUTINE CCG')
    1 NNF=MIN0(N+N1,NMAX)
      MIN=NPN6-N
      MF=NPN6+N
      CALL FACT(FAC)
      CALL SIGNUM(SSIGN)
      IF(K1.EQ.1.AND.K2.EQ.0) MIN=NPN6
      DO 100 MIND=MIN,MF
         M=MIND-NPN6
         MM=M*K1+K2
         M1=MM-M
         IF(IABS(M1).GT.N1) GO TO 90
         NNL=MAX0(IABS(MM),IABS(N-N1))
         IF(NNL.GT.NNF) GO TO 90
         NNU=N+N1
         NNM=(NNU+NNL)*0.5D0
         IF (NNU.EQ.NNL) NNM=NNL
         
         CALL CCGIN(N,N1,M,MM,C, SSIGN, FAC)  !!!!!!!!!!!!!!!!!!
!         print *, 'INSIDE CCG', N, N1, M, MM, C
         CU(NNL)=C  
         IF (NNL.EQ.NNF) GO TO 50
         C2=0D0
         C1=C
         DO 7 NN=NNL+1,MIN0(NNM,NNF)
            A=DFLOAT((NN+MM)*(NN-MM)*(N1-N+NN))
            A=A*DFLOAT((N-N1+NN)*(N+N1-NN+1)*(N+N1+NN+1))
            A=DFLOAT(4*NN*NN)/A
            A=A*DFLOAT((2*NN+1)*(2*NN-1))
            A=DSQRT(A)
            B=0.5D0*DFLOAT(M-M1)
            D=0D0
            IF(NN.EQ.1) GO TO 5
            B=DFLOAT(2*NN*(NN-1))
            B=DFLOAT((2*M-MM)*NN*(NN-1)-MM*N*(N+1)+ MM*N1*(N1+1))/B
            D=DFLOAT(4*(NN-1)*(NN-1))
            D=D*DFLOAT((2*NN-3)*(2*NN-1))
            D=DFLOAT((NN-MM-1)*(NN+MM-1)*(N1-N+NN-1))/D
            D=D*DFLOAT((N-N1+NN-1)*(N+N1-NN+2)*(N+N1+NN))
            D=DSQRT(D)
    5       C=A*(B*C1-D*C2)
            C2=C1
            C1=C
            CU(NN)=C
    7    CONTINUE
         IF (NNF.LE.NNM) GO TO 50

         CALL DIRECT(N,M,N1,M1,NNU,MM,C,FAC) !!!!!!!!!!!!!!!!!!!!!!

         CD(NNU)=C
         IF (NNU.EQ.NNM+1) GO TO 50
         C2=0D0
         C1=C
         DO 12 NN=NNU-1,NNM+1,-1
            A=DFLOAT((NN-MM+1)*(NN+MM+1)*(N1-N+NN+1))
            A=A*DFLOAT((N-N1+NN+1)*(N+N1-NN)*(N+N1+NN+2))
            A=DFLOAT(4*(NN+1)*(NN+1))/A
            A=A*DFLOAT((2*NN+1)*(2*NN+3))
            A=DSQRT(A)
            B=DFLOAT(2*(NN+2)*(NN+1))
            B=DFLOAT((2*M-MM)*(NN+2)*(NN+1)-MM*N*(N+1) + MM*N1*(N1+1))/B
            D=DFLOAT(4*(NN+2)*(NN+2))
            D=D*DFLOAT((2*NN+5)*(2*NN+3))
            D=DFLOAT((NN+MM+2)*(NN-MM+2)*(N1-N+NN+2))/D
            D=D*DFLOAT((N-N1+NN+2)*(N+N1-NN-1)*(N+N1+NN+3))
            D=DSQRT(D)
            C=A*(B*C1-D*C2)
            C2=C1
            C1=C
            CD(NN)=C
   12    CONTINUE
   50    DO 9 NN=NNL,NNF
            IF (NN.LE.NNM) GG(MIND,NN+1)=CU(NN)
            IF (NN.GT.NNM) GG(MIND,NN+1)=CD(NN)
!           WRITE (6,*) N,M,N1,M1,NN,MM,GG(MIND,NN+1)
    9    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
 
!*********************************************************************
 
      SUBROUTINE DIRECT (N,M,N1,M1,NN,MM,C, FAC) !!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 FAC(900)
      !COMMON /FAC/ F
      C=FAC(2*N+1)+FAC(2*N1+1)+FAC(N+N1+M+M1+1)+FAC(N+N1-M-M1+1)    
      C=C-FAC(2*(N+N1)+1)-FAC(N+M+1)-FAC(N-M+1)-FAC(N1+M1+1)-FAC(N1-M1+1)
      C=DEXP(C)
      RETURN
      END
 
!*********************************************************************
!
!   CALCULATION OF THE CLEBCSH-GORDAN COEFFICIENTS
!   G=(N,M:N1,MM-M/NN,MM)
!   FOR GIVEN N,N1,M,MM, WHERE NN=MAX(/MM/,/N-N1/)
!                               /M/.LE.N
!                               /MM-M/.LE.N1
!                               /MM/.LE.N+N1
 
      SUBROUTINE CCGIN(N,N1,M,MM,G,SSIGN,F) !!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8, intent(in)            :: F(900),SSIGN(900)
      !COMMON /SS/ SSIGN
      !COMMON /FAC/ F
      M1=MM-M
      IF(N.GE.IABS(M).AND.N1.GE.IABS(M1).AND.IABS(MM).LE.(N+N1)) GO TO 1
      PRINT 5001
      STOP
 5001 FORMAT(' ERROR IN SUBROUTINE CCGIN')
    1 IF (IABS(MM).GT.IABS(N-N1)) GO TO 100
      L1=N
      L2=N1
      L3=M
      IF(N1.LE.N) GO TO 50
      K=N
      N=N1
      N1=K
      K=M
      M=M1
      M1=K
   50 N2=N*2
      M2=M*2
      N12=N1*2
      M12=M1*2
!      print *, 'SSIGN', N1,M1,SSIGN(N1+M1+1)
      G=SSIGN(N1+M1+1) &
     & *DEXP(F(N+M+1)+F(N-M+1)+F(N12+1)+F(N2-N12+2)-F(N2+2)  &
     &       -F(N1+M1+1)-F(N1-M1+1)-F(N-N1+MM+1)-F(N-N1-MM+1))
      N=L1
      N1=L2
      M=L3
      RETURN
  100 A=1D0
      L1=M
      L2=MM
      IF(MM.GE.0) GO TO 150
      MM=-MM
      M=-M
      M1=-M1
      A=SSIGN(MM+N+N1+1)
  150 G=A*SSIGN(N+M+1) &
     &   *DEXP(F(2*MM+2)+F(N+N1-MM+1)+F(N+M+1)+F(N1+M1+1)  &
     &        -F(N+N1+MM+2)-F(N-N1+MM+1)-F(-N+N1+MM+1)-F(N-M+1) &
     &        -F(N1-M1+1))
      M=L1
      MM=L2
      RETURN
      END

end module
