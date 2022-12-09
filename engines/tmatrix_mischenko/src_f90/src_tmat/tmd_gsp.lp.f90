

module tmat_gsp

contains

!********************************************************************
!                                                                   *
!   CALCULATION OF THE EXPANSION COEFFICIENTS FOR (I,Q,U,V) -       *
!   REPRESENTATION.                                                 *
!                                                                   *
!   INPUT PARAMETERS:                                               *
!                                                                   *
!      LAM - WAVELENGTH OF LIGHT                                    *
!      CSCA - SCATTERING CROSS SECTION                              *
!      TR AND TI - ELEMENTS OF THE T-MATRIX. TRANSFERRED THROUGH    *
!                  COMMON /CTM/                                     *
!      NMAX - DIMENSION OF T(M)-MATRICES                            *
!                                                                   *
!   OUTPUT INFORTMATION:                                            *
!                                                                   *
!      ALF1,...,ALF4,BET1,BET2 - EXPANSION COEFFICIENTS             *
!      LMAX - NUMBER OF COEFFICIENTS MINUS 1                        *
!                                                                   *
!********************************************************************
 
      SUBROUTINE GSP(NMAX,CSCA,LAM,ALF1,ALF2,ALF3,ALF4,BET1,BET2,LMAX)
      use tmat_parameters
      use tmat_CGcoeff
      use tmat_tmat


      IMPLICIT REAL*8 (A-B,D-H,O-Z),COMPLEX*16 (C)

      REAL*8 :: LAM, SSIGN(900), F(900)
      REAL*8 :: CSCA,SSI(NPL),SSJ(NPN1)
      REAL*8 :: ALF1(NPL),ALF2(NPL),ALF3(NPL)
      REAL*8 :: ALF4(NPL),BET1(NPL),BET2(NPL)
      REAL*8 :: TR1(NPL1,NPN4),TR2(NPL1,NPN4)
      REAL*8 :: TI1(NPL1,NPN4),TI2(NPL1,NPN4)
      REAL*8 :: G1(NPL1,NPN6),G2(NPL1,NPN6)
      REAL*8 :: AR1(NPN4),AR2(NPN4),AI1(NPN4),AI2(NPN4)
      REAL*8 :: FR(NPN4,NPN4),FI(NPN4,NPN4),FF(NPN4,NPN4)
      REAL*8 :: B1R(NPL1,NPL1,NPN4),B1I(NPL1,NPL1,NPN4)
      REAL*8 :: B2R(NPL1,NPL1,NPN4),B2I(NPL1,NPL1,NPN4)
      REAL*8 :: D1(NPL1,NPN4,NPN4),D2(NPL1,NPN4,NPN4)
      REAl*8 :: D3(NPL1,NPN4,NPN4),D4(NPL1,NPN4,NPN4)
      REAL*8 :: D5R(NPL1,NPN4,NPN4),D5I(NPL1,NPN4,NPN4)
      REAL*8 :: PLUS1(NPN6*NPN4*NPN4*8)         
!      REAL*4
!     &     TR11(NPN6,NPN4,NPN4),TR12(NPN6,NPN4,NPN4),
!     &     TR21(NPN6,NPN4,NPN4),TR22(NPN6,NPN4,NPN4),
!     &     TI11(NPN6,NPN4,NPN4),TI12(NPN6,NPN4,NPN4),
!     &     TI21(NPN6,NPN4,NPN4),TI22(NPN6,NPN4,NPN4)
      COMPLEX*16 CIM(NPN1)

!      COMMON /TMAT/ TR11,TR12,TR21,TR22,TI11,TI12,TI21,TI22
!      COMMON /CBESS/ B1R,B1I,B2R,B2I    
!      COMMON /SS/ SSIGN
!      EQUIVALENCE (PLUS1(1),TR11(1,1,1) )
!      EQUIVALENCE (D1(1,1,1),PLUS1(1)),        
!     &            (D2(1,1,1),PLUS1(NPL1*NPN4*NPN4+1)),
!     &            (D3(1,1,1),PLUS1(NPL1*NPN4*NPN4*2+1)),
!     &            (D4(1,1,1),PLUS1(NPL1*NPN4*NPN4*3+1)), 
!     &            (D5R(1,1,1),PLUS1(NPL1*NPN4*NPN4*4+1)) 

!      print *, 'CALL-FACT'
      CALL FACT(F)
!      print *, 'CALL-SIGNUM'
      CALL SIGNUM(SSIGN)
      LMAX=2*NMAX
      L1MAX=LMAX+1
!      print *, L1MAX, NMAX
      CI=(0D0,1D0)
      CIM(1)=CI
      DO 2 I=2,NMAX
         CIM(I)=CIM(I-1)*CI
    2 CONTINUE
      SSI(1)=1D0
      DO 3 I=1,LMAX
         I1=I+1
         SI=DFLOAT(2*I+1)
         SSI(I1)=SI
         IF(I.LE.NMAX) SSJ(I)=DSQRT(SI)
    3 CONTINUE
      CI=-CI
      DO 5 I=1,NMAX
         SI=SSJ(I)
         CCI=CIM(I)
         DO 4 J=1,NMAX
            SJ=1D0/SSJ(J)
            CCJ=CIM(J)*SJ/CCI
            FR(J,I)=CCJ
            FI(J,I)=CCJ*CI
            FF(J,I)=SI*SJ
    4    CONTINUE
    5 CONTINUE
      NMAX1=NMAX+1
 
! *****  CALCULATION OF THE ARRAYS B1 AND B2  *****
 
      K1=1
      K2=0
      K3=0
      K4=1
      K5=1
      K6=2
 
!     PRINT 3300, B1,B2
 3300 FORMAT (' B1 AND B2')
      DO 100 N=1,NMAX
 
! *****  CALCULATION OF THE ARRAYS T1 AND T2  *****
 
 
         DO 10 NN=1,NMAX
            M1MAX=MIN0(N,NN)+1
            DO 6 M1=1,M1MAX
               M=M1-1
               L1=NPN6+M
               TT1=TR11(M1,N,NN)
               TT2=TR12(M1,N,NN)
               TT3=TR21(M1,N,NN)
               TT4=TR22(M1,N,NN)
               TT5=IT11(M1,N,NN)
               TT6=IT12(M1,N,NN)
               TT7=IT21(M1,N,NN)
               TT8=IT22(M1,N,NN)
               T1=TT1+TT2
               T2=TT3+TT4
               T3=TT5+TT6
               T4=TT7+TT8
               TR1(L1,NN)=T1+T2
               TR2(L1,NN)=T1-T2
               TI1(L1,NN)=T3+T4
               TI2(L1,NN)=T3-T4
               IF(M.EQ.0) GO TO 6
               L1=NPN6-M
               T1=TT1-TT2
               T2=TT3-TT4
               T3=TT5-TT6
               T4=TT7-TT8
               TR1(L1,NN)=T1-T2
               TR2(L1,NN)=T1+T2
               TI1(L1,NN)=T3-T4
               TI2(L1,NN)=T3+T4
    6       CONTINUE
   10    CONTINUE
 
!  *****  END OF THE CALCULATION OF THE ARRAYS T1 AND T2  *****
 
         NN1MAX=NMAX1+N
         DO 40 NN1=1,NN1MAX
            N1=NN1-1
 
!  *****  CALCULATION OF THE ARRAYS A1 AND A2  *****
!            print *, 'CALL-CCG1'
            CALL CCG(N,N1,NMAX,K1,K2,G1)
            NNMAX=MIN0(NMAX,N1+N)
            NNMIN=MAX0(1,IABS(N-N1))
            KN=N+NN1
            DO 15 NN=NNMIN,NNMAX
               NNN=NN+1
               SIG=SSIGN(KN+NN)
               M1MAX=MIN0(N,NN)+NPN6
               AAR1=0D0
               AAR2=0D0
               AAI1=0D0
               AAI2=0D0
               DO 13 M1=NPN6,M1MAX
                  M=M1-NPN6
                  SSS=G1(M1,NNN)
                  RR1=TR1(M1,NN)
                  RI1=TI1(M1,NN)
                  RR2=TR2(M1,NN)
                  RI2=TI2(M1,NN)
                  IF(M.EQ.0) GO TO 12
                  M2=NPN6-M
                  RR1=RR1+TR1(M2,NN)*SIG
                  RI1=RI1+TI1(M2,NN)*SIG
                  RR2=RR2+TR2(M2,NN)*SIG
                  RI2=RI2+TI2(M2,NN)*SIG
   12             AAR1=AAR1+SSS*RR1
                  AAI1=AAI1+SSS*RI1
                  AAR2=AAR2+SSS*RR2
                  AAI2=AAI2+SSS*RI2
   13          CONTINUE
               XR=FR(NN,N)
               XI=FI(NN,N)
               AR1(NN)=AAR1*XR-AAI1*XI
               AI1(NN)=AAR1*XI+AAI1*XR
               AR2(NN)=AAR2*XR-AAI2*XI
               AI2(NN)=AAR2*XI+AAI2*XR
   15       CONTINUE
 
!  *****  END OF THE CALCULATION OF THE ARRAYS A1 AND A2 ****
 
            CALL CCG(N,N1,NMAX,K3,K4,G2)
!            print *,'CALL-CCG2', N, N1, NMAX, K3, K4
            M1=MAX0(-N1+1,-N)
            M2=MIN0(N1+1,N)
            M1MAX=M2+NPN6
            M1MIN=M1+NPN6
            DO 30 M1=M1MIN,M1MAX
               BBR1=0D0
               BBI1=0D0
               BBR2=0D0
               BBI2=0D0
               DO 25 NN=NNMIN,NNMAX
                  NNN=NN+1
                  SSS=G2(M1,NNN)
                  BBR1=BBR1+SSS*AR1(NN)
                  BBI1=BBI1+SSS*AI1(NN)
                  BBR2=BBR2+SSS*AR2(NN)
                  BBI2=BBI2+SSS*AI2(NN)
   25          CONTINUE
               B1R(NN1,M1,N)=BBR1
               B1I(NN1,M1,N)=BBI1
               B2R(NN1,M1,N)=BBR2
               B2I(NN1,M1,N)=BBI2
   30       CONTINUE
   40    CONTINUE
  100 CONTINUE
 
!  *****  END OF THE CALCULATION OF THE ARRAYS B1 AND B2 ****
 
!  *****  CALCULATION OF THE ARRAYS D1,D2,D3,D4, AND D5  *****
 
!     PRINT 3301
! 3301 FORMAT(' D1, D2, ...')
      DO 200 N=1,NMAX
         DO 190 NN=1,NMAX
            M1=MIN0(N,NN)
            M1MAX=NPN6+M1
            M1MIN=NPN6-M1
            NN1MAX=NMAX1+MIN0(N,NN)
            DO 180 M1=M1MIN,M1MAX
               M=M1-NPN6
               NN1MIN=IABS(M-1)+1
               DD1=0D0
               DD2=0D0
               DO 150 NN1=NN1MIN,NN1MAX
                  XX=SSI(NN1)
                  X1=B1R(NN1,M1,N)
                  X2=B1I(NN1,M1,N)
                  X3=B1R(NN1,M1,NN)
                  X4=B1I(NN1,M1,NN)
                  X5=B2R(NN1,M1,N)
                  X6=B2I(NN1,M1,N)
                  X7=B2R(NN1,M1,NN)
                  X8=B2I(NN1,M1,NN)
                  DD1=DD1+XX*(X1*X3+X2*X4)
                  DD2=DD2+XX*(X5*X7+X6*X8)
  150          CONTINUE
               D1(M1,NN,N)=DD1
               D2(M1,NN,N)=DD2
  180       CONTINUE
            MMAX=MIN0(N,NN+2)
            MMIN=MAX0(-N,-NN+2)
            M1MAX=NPN6+MMAX
            M1MIN=NPN6+MMIN
            DO 186 M1=M1MIN,M1MAX
               M=M1-NPN6
               NN1MIN=IABS(M-1)+1
               DD3=0D0
               DD4=0D0
               DD5R=0D0
               DD5I=0D0
               M2=-M+2+NPN6
               DO 183 NN1=NN1MIN,NN1MAX
                  XX=SSI(NN1)
                  X1=B1R(NN1,M1,N)
                  X2=B1I(NN1,M1,N)
                  X3=B2R(NN1,M1,N)
                  X4=B2I(NN1,M1,N)
                  X5=B1R(NN1,M2,NN)
                  X6=B1I(NN1,M2,NN)
                  X7=B2R(NN1,M2,NN)
                  X8=B2I(NN1,M2,NN)
                  DD3=DD3+XX*(X1*X5+X2*X6)
                  DD4=DD4+XX*(X3*X7+X4*X8)
                  DD5R=DD5R+XX*(X3*X5+X4*X6)
                  DD5I=DD5I+XX*(X4*X5-X3*X6)
  183          CONTINUE
               D3(M1,NN,N)=DD3
               D4(M1,NN,N)=DD4
               D5R(M1,NN,N)=DD5R
               D5I(M1,NN,N)=DD5I
  186       CONTINUE
  190    CONTINUE
  200 CONTINUE
 
!  *****  END OF THE CALCULATION OF THE D-ARRAYS *****
 
!  *****  CALCULATION OF THE EXPANSION COEFFICIENTS *****
 
!     PRINT 3303
! 3303 FORMAT (' G1, G2, ...')
 
      DK=LAM*LAM/(4D0*CSCA*DACOS(-1D0))
      DO 300 L1=1,L1MAX
         G1L=0D0
         G2L=0D0
         G3L=0D0
         G4L=0D0
         G5LR=0D0
         G5LI=0D0
!         print *, '(300)L1L1MAX',L1, L1MAX
         L=L1-1
         SL=SSI(L1)*DK
         DO 290 N=1,NMAX
            NNMIN=MAX0(1,IABS(N-L))
            NNMAX=MIN0(NMAX,N+L)
            IF(NNMAX.LT.NNMIN) GO TO 290
            CALL CCG(N,L,NMAX,K1,K2,G1)
!            print *, 'CALL-CCG3', N,L,NMAX,K1,K2
!            print *, 'CALL-CCG3G1', G1(NPN6+1,NNN), NPN6,NNN 
            IF(L.GE.2) THEN
!              print *, 'CALL-CCG4', L
              CALL CCG(N,L,NMAX,K5,K6,G2)
            endif
            NL=N+L
            DO 280  NN=NNMIN,NNMAX
               NNN=NN+1
               MMAX=MIN0(N,NN)
               M1MIN=NPN6-MMAX
               M1MAX=NPN6+MMAX
               SI=SSIGN(NL+NNN)
               DM1=0D0
               DM2=0D0
               DO 270 M1=M1MIN,M1MAX
                  M=M1-NPN6
                  IF(M.GE.0) SSS1=G1(M1,NNN)
                  IF(M.LT.0) SSS1=G1(NPN6-M,NNN)*SI
                  DM1=DM1+SSS1*D1(M1,NN,N)
                  DM2=DM2+SSS1*D2(M1,NN,N)
  270          CONTINUE
               FFN=FF(NN,N)
               SSS=G1(NPN6+1,NNN)*FFN
!               print *, 'G1, FFN', G1(NPN6+1,NNN),FFN,NNN
!               print *,'G1L,SSS,DM1',G1L,SSS,DM1
               G1L=G1L+SSS*DM1
               G2L=G2L+SSS*DM2*SI
               IF(L.LT.2) GO TO 280
               DM3=0D0
               DM4=0D0
               DM5R=0D0
               DM5I=0D0
               MMAX=MIN0(N,NN+2)
               MMIN=MAX0(-N,-NN+2)
               M1MAX=NPN6+MMAX
               M1MIN=NPN6+MMIN
               DO 275 M1=M1MIN,M1MAX
                  M=M1-NPN6
                  SSS1=G2(NPN6-M,NNN)
                  DM3=DM3+SSS1*D3(M1,NN,N)
                  DM4=DM4+SSS1*D4(M1,NN,N)
                  DM5R=DM5R+SSS1*D5R(M1,NN,N)
                  DM5I=DM5I+SSS1*D5I(M1,NN,N)
  275          CONTINUE
               G5LR=G5LR-SSS*DM5R
               G5LI=G5LI-SSS*DM5I
               SSS=G2(NPN4,NNN)*FFN
               G3L=G3L+SSS*DM3
               G4L=G4L+SSS*DM4*SI
  280       CONTINUE
  290    CONTINUE
!         print *,'G1L, SL',G1L, SL
         G1L=G1L*SL
         G2L=G2L*SL
         G3L=G3L*SL
         G4L=G4L*SL
         G5LR=G5LR*SL
         G5LI=G5LI*SL
         ALF1(L1)=G1L+G2L
         ALF2(L1)=G3L+G4L
         ALF3(L1)=G3L-G4L
         ALF4(L1)=G1L-G2L
         BET1(L1)=G5LR*2D0
         BET2(L1)=G5LI*2D0
         LMAX=L
!         print *, 'DABS(G1L)',DABS(G1L)
         IF(DABS(G1L).LT.1D-6) GO TO 500
  300 CONTINUE
  500 CONTINUE
      RETURN
      END


end module tmat_gsp
