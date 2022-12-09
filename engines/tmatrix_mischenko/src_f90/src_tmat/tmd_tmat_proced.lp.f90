
module tmat_shared_tmat2

      use tmat_parameters

      REAL*8  :: R11(NPN1,NPN1),R12(NPN1,NPN1)
      REAL*8  :: R21(NPN1,NPN1),R22(NPN1,NPN1)
      REAL*8  :: I11(NPN1,NPN1),I12(NPN1,NPN1)
      REAL*8  :: I21(NPN1,NPN1),I22(NPN1,NPN1)
      REAL*8  :: RG11(NPN1,NPN1),RG12(NPN1,NPN1)
      REAL*8  :: RG21(NPN1,NPN1),RG22(NPN1,NPN1)
      REAL*8  :: IG11(NPN1,NPN1),IG12(NPN1,NPN1)
      REAL*8  :: IG21(NPN1,NPN1),IG22(NPN1,NPN1)

end module tmat_shared_tmat2

      
module tmat_procedures


  contains

!**********************************************************************
 
      SUBROUTINE TMATR0 ( NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,DDR, &
                        & DRR,DRI,NMAX,NCHECK)
      use tmat_parameters
      use tmat_shared_bessel
      use tmat_ctt
      use tmat_ct
      use tmat_shared_tmat2

      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  :: X(NPNG2),W(NPNG2),AN(NPN1),S(NPNG2),SS(NPNG2)
      REAL*8  :: R(NPNG2),DR(NPNG2),SIG(NPN2)
      REAL*8  :: DDR(NPNG2),DRR(NPNG2)
      REAL*8  :: D1(NPNG2,NPN1),D2(NPNG2,NPN1)
      REAL*8  :: DRI(NPNG2),DS(NPNG2),DSS(NPNG2),RR(NPNG2)
      REAL*8  :: DV1(NPN1),DV2(NPN1)
 
      !REAL*8  :: R11(NPN1,NPN1),R12(NPN1,NPN1)
      !REAL*8  :: R21(NPN1,NPN1),R22(NPN1,NPN1)
      !REAL*8  :: I11(NPN1,NPN1),I12(NPN1,NPN1)
      !REAL*8  :: I21(NPN1,NPN1),I22(NPN1,NPN1)
      !REAL*8  :: RG11(NPN1,NPN1),RG12(NPN1,NPN1)
      !REAL*8  :: RG21(NPN1,NPN1),RG22(NPN1,NPN1)
      !REAL*8  :: IG11(NPN1,NPN1),IG12(NPN1,NPN1)
      !REAL*8  :: IG21(NPN1,NPN1),IG22(NPN1,NPN1)
      
      REAL*8  :: ANN(NPN1,NPN1)
      REAL*8  :: TQR(NPN2,NPN2),TQI(NPN2,NPN2)
      REAL*8  :: TRGQR(NPN2,NPN2),TRGQI(NPN2,NPN2)
      !REAL*8  :: PLUS(NPN6*NPN4*NPN4*8)
     ! COMMON /TMAT2/ PLUS,
     !&            R11,R12,R21,R22,I11,I12,I21,I22,RG11,RG12,RG21,RG22,
     !&            IG11,IG12,IG21,IG22
      MM1=1
      NNMAX=NMAX+NMAX
      NG=2*NGAUSS
      NGSS=NG
      FACTOR=1D0
      IF (NCHECK.EQ.1) THEN
            NGSS=NGAUSS
            FACTOR=2D0
         ELSE
            CONTINUE
      ENDIF
      SI=1D0
      DO 5 N=1,NNMAX
           SI=-SI
           SIG(N)=SI
    5 CONTINUE
   20 DO 25 I=1,NGAUSS
         I1=NGAUSS+I
         I2=NGAUSS-I+1
         CALL VIG ( X(I1), NMAX, 0, DV1, DV2)
         DO 25 N=1,NMAX
            SI=SIG(N)
            DD1=DV1(N)
            DD2=DV2(N)
            D1(I1,N)=DD1
            D2(I1,N)=DD2
            D1(I2,N)=DD1*SI
            D2(I2,N)=-DD2*SI
   25 CONTINUE
   30 DO 40 I=1,NGSS
           RR(I)=W(I)*R(I)
   40 CONTINUE
 
      DO 300  N1=MM1,NMAX
           AN1=AN(N1)
           DO 300 N2=MM1,NMAX
                AN2=AN(N2)
                AR12=0D0
                AR21=0D0
                AI12=0D0
                AI21=0D0
                GR12=0D0
                GR21=0D0
                GI12=0D0
                GI21=0D0
                IF (NCHECK.EQ.1.AND.SIG(N1+N2).LT.0D0) GO TO 205

                DO 200 I=1,NGSS
                    D1N1=D1(I,N1)
                    D2N1=D2(I,N1)
                    D1N2=D1(I,N2)
                    D2N2=D2(I,N2)
                    A12=D1N1*D2N2
                    A21=D2N1*D1N2
                    A22=D2N1*D2N2
                    AA1=A12+A21
 
                    QJ1=J(I,N1)
                    QY1=Y(I,N1)
                    QJR2=JR(I,N2)
                    QJI2=JI(I,N2)
                    QDJR2=DJR(I,N2)
                    QDJI2=DJI(I,N2)
                    QDJ1=DJ(I,N1)
                    QDY1=DY(I,N1)
 
                    C1R=QJR2*QJ1
                    C1I=QJI2*QJ1
                    B1R=C1R-QJI2*QY1
                    B1I=C1I+QJR2*QY1
 
                    C2R=QJR2*QDJ1
                    C2I=QJI2*QDJ1
                    B2R=C2R-QJI2*QDY1
                    B2I=C2I+QJR2*QDY1
 
                    DDRI=DDR(I)
                    C3R=DDRI*C1R
                    C3I=DDRI*C1I
                    B3R=DDRI*B1R
                    B3I=DDRI*B1I
 
                    C4R=QDJR2*QJ1
                    C4I=QDJI2*QJ1
                    B4R=C4R-QDJI2*QY1
                    B4I=C4I+QDJR2*QY1
 
                    DRRI=DRR(I)
                    DRII=DRI(I)
                    C5R=C1R*DRRI-C1I*DRII
                    C5I=C1I*DRRI+C1R*DRII
                    B5R=B1R*DRRI-B1I*DRII
                    B5I=B1I*DRRI+B1R*DRII
 
                    URI=DR(I)
                    RRI=RR(I)
 
                    F1=RRI*A22
                    F2=RRI*URI*AN1*A12
                    AR12=AR12+F1*B2R+F2*B3R
                    AI12=AI12+F1*B2I+F2*B3I
                    GR12=GR12+F1*C2R+F2*C3R
                    GI12=GI12+F1*C2I+F2*C3I
 
                    F2=RRI*URI*AN2*A21
                    AR21=AR21+F1*B4R+F2*B5R
                    AI21=AI21+F1*B4I+F2*B5I
                    GR21=GR21+F1*C4R+F2*C5R
                    GI21=GI21+F1*C4I+F2*C5I
  200           CONTINUE
 
  205           AN12=ANN(N1,N2)*FACTOR
                R12(N1,N2)=AR12*AN12
                R21(N1,N2)=AR21*AN12
                I12(N1,N2)=AI12*AN12
                I21(N1,N2)=AI21*AN12
                RG12(N1,N2)=GR12*AN12
                RG21(N1,N2)=GR21*AN12
                IG12(N1,N2)=GI12*AN12
                IG21(N1,N2)=GI21*AN12
  300 CONTINUE
 
      TPIR=PIR
      TPII=PII
      TPPI=PPI
 
      NM=NMAX
      DO 310 N1=MM1,NMAX
           K1=N1-MM1+1
           KK1=K1+NM
           DO 310 N2=MM1,NMAX
                K2=N2-MM1+1
                KK2=K2+NM
 
                TAR12= I12(N1,N2)
                TAI12=-R12(N1,N2)
                TGR12= IG12(N1,N2)
                TGI12=-RG12(N1,N2)
 
                TAR21=-I21(N1,N2)
                TAI21= R21(N1,N2)
                TGR21=-IG21(N1,N2)
                TGI21= RG21(N1,N2)
 
                TQR(K1,K2)=TPIR*TAR21-TPII*TAI21+TPPI*TAR12
                TQI(K1,K2)=TPIR*TAI21+TPII*TAR21+TPPI*TAI12
                TRGQR(K1,K2)=TPIR*TGR21-TPII*TGI21+TPPI*TGR12
                TRGQI(K1,K2)=TPIR*TGI21+TPII*TGR21+TPPI*TGI12
 
                TQR(K1,KK2)=0D0
                TQI(K1,KK2)=0D0
                TRGQR(K1,KK2)=0D0
                TRGQI(K1,KK2)=0D0
 
                TQR(KK1,K2)=0D0
                TQI(KK1,K2)=0D0
                TRGQR(KK1,K2)=0D0
                TRGQI(KK1,K2)=0D0
 
                TQR(KK1,KK2)=TPIR*TAR12-TPII*TAI12+TPPI*TAR21
                TQI(KK1,KK2)=TPIR*TAI12+TPII*TAR12+TPPI*TAI21
                TRGQR(KK1,KK2)=TPIR*TGR12-TPII*TGI12+TPPI*TGR21
                TRGQI(KK1,KK2)=TPIR*TGI12+TPII*TGR12+TPPI*TGI21
  310 CONTINUE
 
      NNMAX=2*NM
      DO 320 N1=1,NNMAX
           DO 320 N2=1,NNMAX
                QR(N1,N2)=TQR(N1,N2)
                QI(N1,N2)=TQI(N1,N2)
                RGQR(N1,N2)=TRGQR(N1,N2)
                RGQI(N1,N2)=TRGQI(N1,N2)
  320 CONTINUE
      
      CALL TT(NMAX,NCHECK)
      RETURN
      END
 
!**********************************************************************
 
      SUBROUTINE TMATR ( M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,DDR, &
                       & DRR,DRI,NMAX,NCHECK)
      use tmat_parameters
      use tmat_shared_bessel
      use tmat_ctt
      use tmat_ct
      use tmat_shared_tmat2

      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 :: X(NPNG2),W(NPNG2),AN(NPN1),S(NPNG2),SS(NPNG2)
      REAL*8 :: R(NPNG2),DR(NPNG2),SIG(NPN2)
      REAL*8 :: DDR(NPNG2),DRR(NPNG2)
      REAL*8 :: D1(NPNG2,NPN1),D2(NPNG2,NPN1)
      REAL*8 :: DRI(NPNG2),DS(NPNG2),DSS(NPNG2),RR(NPNG2)
      REAL*8 :: DV1(NPN1),DV2(NPN1)
      !REAL*8 :: R11(NPN1,NPN1),R12(NPN1,NPN1)
      !REAL*8 :: R21(NPN1,NPN1),R22(NPN1,NPN1)
      !REAL*8 :: I11(NPN1,NPN1),I12(NPN1,NPN1)
      !REAL*8 :: I21(NPN1,NPN1),I22(NPN1,NPN1)
      !REAL*8 :: RG11(NPN1,NPN1),RG12(NPN1,NPN1)
      !REAL*8 :: RG21(NPN1,NPN1),RG22(NPN1,NPN1)
      !REAL*8 :: IG11(NPN1,NPN1),IG12(NPN1,NPN1)
      !REAL*8 :: IG21(NPN1,NPN1),IG22(NPN1,NPN1)
      REAL*8 :: ANN(NPN1,NPN1)
      REAL*8 :: TQR(NPN2,NPN2),TQI(NPN2,NPN2)
      REAL*8 :: TRGQR(NPN2,NPN2),TRGQI(NPN2,NPN2)
      !REAL*8 :: TR1(NPN2,NPN2),TI1(NPN2,NPN2)
      !REAL*8 :: PLUS(NPN6*NPN4*NPN4*8)
      !COMMON /TMAT2/ PLUS,
      !&            R11,R12,R21,R22,I11,I12,I21,I22,RG11,RG12,RG21,RG22,
      !&            IG11,IG12,IG21,IG22
      MM1=M
      QM=DFLOAT(M)
      QMM=QM*QM
      NG=2*NGAUSS
      NGSS=NG
      FACTOR=1D0
      IF (NCHECK.EQ.1) THEN
            NGSS=NGAUSS
            FACTOR=2D0
         ELSE
            CONTINUE
      ENDIF
      SI=1D0
      NM=NMAX+NMAX
      DO 5 N=1,NM
           SI=-SI
           SIG(N)=SI
    5 CONTINUE
   20 DO 25 I=1,NGAUSS
         I1=NGAUSS+I
         I2=NGAUSS-I+1
         CALL VIG (X(I1),NMAX,M,DV1,DV2)
         DO 25 N=1,NMAX
            SI=SIG(N)
            DD1=DV1(N)
            DD2=DV2(N)
            D1(I1,N)=DD1
            D2(I1,N)=DD2
            D1(I2,N)=DD1*SI
            D2(I2,N)=-DD2*SI
   25 CONTINUE
   30 DO 40 I=1,NGSS
           WR=W(I)*R(I)
           DS(I)=S(I)*QM*WR
           DSS(I)=SS(I)*QMM
           RR(I)=WR
   40 CONTINUE
 
      DO 300  N1=MM1,NMAX
           AN1=AN(N1)
           DO 300 N2=MM1,NMAX
                AN2=AN(N2)
                AR11=0D0
                AR12=0D0
                AR21=0D0
                AR22=0D0
                AI11=0D0
                AI12=0D0
                AI21=0D0
                AI22=0D0
                GR11=0D0
                GR12=0D0
                GR21=0D0
                GR22=0D0
                GI11=0D0
                GI12=0D0
                GI21=0D0
                GI22=0D0
                SI=SIG(N1+N2)
 
                DO 200 I=1,NGSS
                    D1N1=D1(I,N1)
                    D2N1=D2(I,N1)
                    D1N2=D1(I,N2)
                    D2N2=D2(I,N2)
                    A11=D1N1*D1N2
                    A12=D1N1*D2N2
                    A21=D2N1*D1N2
                    A22=D2N1*D2N2
                    AA1=A12+A21
                    AA2=A11*DSS(I)+A22
                    QJ1=J(I,N1)
                    QY1=Y(I,N1)
                    QJR2=JR(I,N2)
                    QJI2=JI(I,N2)
                    QDJR2=DJR(I,N2)
                    QDJI2=DJI(I,N2)
                    QDJ1=DJ(I,N1)
                    QDY1=DY(I,N1)
 
                    C1R=QJR2*QJ1
                    C1I=QJI2*QJ1
                    B1R=C1R-QJI2*QY1
                    B1I=C1I+QJR2*QY1
 
                    C2R=QJR2*QDJ1
                    C2I=QJI2*QDJ1
                    B2R=C2R-QJI2*QDY1
                    B2I=C2I+QJR2*QDY1
 
                    DDRI=DDR(I)
                    C3R=DDRI*C1R
                    C3I=DDRI*C1I
                    B3R=DDRI*B1R
                    B3I=DDRI*B1I
 
                    C4R=QDJR2*QJ1
                    C4I=QDJI2*QJ1
                    B4R=C4R-QDJI2*QY1
                    B4I=C4I+QDJR2*QY1
 
                    DRRI=DRR(I)
                    DRII=DRI(I)
                    C5R=C1R*DRRI-C1I*DRII
                    C5I=C1I*DRRI+C1R*DRII
                    B5R=B1R*DRRI-B1I*DRII
                    B5I=B1I*DRRI+B1R*DRII
 
                    C6R=QDJR2*QDJ1
                    C6I=QDJI2*QDJ1
                    B6R=C6R-QDJI2*QDY1
                    B6I=C6I+QDJR2*QDY1
 
                    C7R=C4R*DDRI
                    C7I=C4I*DDRI
                    B7R=B4R*DDRI
                    B7I=B4I*DDRI
 
                    C8R=C2R*DRRI-C2I*DRII
                    C8I=C2I*DRRI+C2R*DRII
                    B8R=B2R*DRRI-B2I*DRII
                    B8I=B2I*DRRI+B2R*DRII
 
                    URI=DR(I)
                    DSI=DS(I)
                    DSSI=DSS(I)
                    RRI=RR(I)
 
                    IF (NCHECK.EQ.1.AND.SI.GT.0D0) GO TO 150
 
                    E1=DSI*AA1
                    AR11=AR11+E1*B1R
                    AI11=AI11+E1*B1I
                    GR11=GR11+E1*C1R
                    GI11=GI11+E1*C1I
                    IF (NCHECK.EQ.1) GO TO 160
 
  150               F1=RRI*AA2
                    F2=RRI*URI*AN1*A12
                    AR12=AR12+F1*B2R+F2*B3R
                    AI12=AI12+F1*B2I+F2*B3I
                    GR12=GR12+F1*C2R+F2*C3R
                    GI12=GI12+F1*C2I+F2*C3I
 
                    F2=RRI*URI*AN2*A21
                    AR21=AR21+F1*B4R+F2*B5R
                    AI21=AI21+F1*B4I+F2*B5I
                    GR21=GR21+F1*C4R+F2*C5R
                    GI21=GI21+F1*C4I+F2*C5I
                    IF (NCHECK.EQ.1) GO TO 200
 
  160               E2=DSI*URI*A11
                    E3=E2*AN2
                    E2=E2*AN1
                    AR22=AR22+E1*B6R+E2*B7R+E3*B8R
                    AI22=AI22+E1*B6I+E2*B7I+E3*B8I
                    GR22=GR22+E1*C6R+E2*C7R+E3*C8R
                    GI22=GI22+E1*C6I+E2*C7I+E3*C8I
  200           CONTINUE
                AN12=ANN(N1,N2)*FACTOR
                R11(N1,N2)=AR11*AN12
                R12(N1,N2)=AR12*AN12
                R21(N1,N2)=AR21*AN12
                R22(N1,N2)=AR22*AN12
                I11(N1,N2)=AI11*AN12
                I12(N1,N2)=AI12*AN12
                I21(N1,N2)=AI21*AN12
                I22(N1,N2)=AI22*AN12
                RG11(N1,N2)=GR11*AN12
                RG12(N1,N2)=GR12*AN12
                RG21(N1,N2)=GR21*AN12
                RG22(N1,N2)=GR22*AN12
                IG11(N1,N2)=GI11*AN12
                IG12(N1,N2)=GI12*AN12
                IG21(N1,N2)=GI21*AN12
                IG22(N1,N2)=GI22*AN12
 
  300 CONTINUE
      TPIR=PIR
      TPII=PII
      TPPI=PPI
      NM=NMAX-MM1+1
      DO 310 N1=MM1,NMAX
           K1=N1-MM1+1
           KK1=K1+NM
           DO 310 N2=MM1,NMAX
                K2=N2-MM1+1
                KK2=K2+NM
 
                TAR11=-R11(N1,N2)
                TAI11=-I11(N1,N2)
                TGR11=-RG11(N1,N2)
                TGI11=-IG11(N1,N2)
 
                TAR12= I12(N1,N2)
                TAI12=-R12(N1,N2)
                TGR12= IG12(N1,N2)
                TGI12=-RG12(N1,N2)
 
                TAR21=-I21(N1,N2)
                TAI21= R21(N1,N2)
                TGR21=-IG21(N1,N2)
                TGI21= RG21(N1,N2)
 
                TAR22=-R22(N1,N2)
                TAI22=-I22(N1,N2)
                TGR22=-RG22(N1,N2)
                TGI22=-IG22(N1,N2)
 
                TQR(K1,K2)=TPIR*TAR21-TPII*TAI21+TPPI*TAR12
                TQI(K1,K2)=TPIR*TAI21+TPII*TAR21+TPPI*TAI12
                TRGQR(K1,K2)=TPIR*TGR21-TPII*TGI21+TPPI*TGR12
                TRGQI(K1,K2)=TPIR*TGI21+TPII*TGR21+TPPI*TGI12
 
                TQR(K1,KK2)=TPIR*TAR11-TPII*TAI11+TPPI*TAR22
                TQI(K1,KK2)=TPIR*TAI11+TPII*TAR11+TPPI*TAI22
                TRGQR(K1,KK2)=TPIR*TGR11-TPII*TGI11+TPPI*TGR22
                TRGQI(K1,KK2)=TPIR*TGI11+TPII*TGR11+TPPI*TGI22
 
                TQR(KK1,K2)=TPIR*TAR22-TPII*TAI22+TPPI*TAR11
                TQI(KK1,K2)=TPIR*TAI22+TPII*TAR22+TPPI*TAI11
                TRGQR(KK1,K2)=TPIR*TGR22-TPII*TGI22+TPPI*TGR11
                TRGQI(KK1,K2)=TPIR*TGI22+TPII*TGR22+TPPI*TGI11
 
                TQR(KK1,KK2)=TPIR*TAR12-TPII*TAI12+TPPI*TAR21
                TQI(KK1,KK2)=TPIR*TAI12+TPII*TAR12+TPPI*TAI21
                TRGQR(KK1,KK2)=TPIR*TGR12-TPII*TGI12+TPPI*TGR21
                TRGQI(KK1,KK2)=TPIR*TGI12+TPII*TGR12+TPPI*TGI21
  310 CONTINUE
 
      NNMAX=2*NM
      DO 320 N1=1,NNMAX
           DO 320 N2=1,NNMAX
                QR(N1,N2)=TQR(N1,N2)
                QI(N1,N2)=TQI(N1,N2)
                RGQR(N1,N2)=TRGQR(N1,N2)
                RGQI(N1,N2)=TRGQI(N1,N2)
  320 CONTINUE
 
!      print *,'CALL-TT',NM, NCHECK
      CALL TT(NM,NCHECK)
 
      RETURN
      END
 
!*****************************************************************
 
      SUBROUTINE VIG (X, NMAX, M, DV1, DV2)
      use tmat_parameters
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 DV1(NPN1),DV2(NPN1)
 
      A=1D0
      QS=DSQRT(1D0-X*X)
      QS1=1D0/QS
      DO N=1,NMAX
         DV1(N)=0D0
         DV2(N)=0D0
      ENDDO   
      IF (M.NE.0) GO TO 20
      D1=1D0
      D2=X  
      DO N=1,NMAX  
         QN=DFLOAT(N)
         QN1=DFLOAT(N+1)
         QN2=DFLOAT(2*N+1)
         D3=(QN2*X*D2-QN*D1)/QN1 
         DER=QS1*(QN1*QN/QN2)*(-D1+D3)
         DV1(N)=D2
         DV2(N)=DER
         D1=D2
         D2=D3
      ENDDO   
      RETURN
   20 QMM=DFLOAT(M*M)
      DO I=1,M
         I2=I*2
         A=A*DSQRT(DFLOAT(I2-1)/DFLOAT(I2))*QS
      ENDDO   
      D1=0D0
      D2=A 
      DO N=M,NMAX
         QN=DFLOAT(N)
         QN2=DFLOAT(2*N+1)
         QN1=DFLOAT(N+1)
         QNM=DSQRT(QN*QN-QMM)
         QNM1=DSQRT(QN1*QN1-QMM)
         D3=(QN2*X*D2-QNM*D1)/QNM1
         DER=QS1*(-QN1*QNM*D1+QN*QNM1*D3)/QN2
         DV1(N)=D2
         DV2(N)=DER
         D1=D2
         D2=D3
      ENDDO   
      RETURN
      END 
 
!**********************************************************************
!                                                                     *
!   CALCULATION OF THE MATRIX    T = - RG(Q) * (Q**(-1))              *
!                                                                     *
!   INPUT INFORTMATION IS IN COMMON /CTT/                             *
!   OUTPUT INFORMATION IS IN COMMON /CT/                              *
!                                                                     *
!**********************************************************************
 
      SUBROUTINE TT(NMAX,NCHECK)

      use tmat_parameters
      use tmat_ctt
      use tmat_ct

      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 :: F(NPN2,NPN2),B(NPN2),WORK(NPN2)
!     *       QR(NPN2,NPN2),QI(NPN2,NPN2),
!     *       RGQR(NPN2,NPN2),RGQI(NPN2,NPN2),
       REAL*8 :: A(NPN2,NPN2),C(NPN2,NPN2),D(NPN2,NPN2),E(NPN2,NPN2)
!      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)
      COMPLEX*16 ZQ(NPN2,NPN2),ZW(NPN2)
      INTEGER IPIV(NPN2),IPVT(NPN2)
!      COMMON /CT/ TR1,TI1
!      COMMON /CTT/ QR,QI,RGQR,RGQI
      NDIM=NPN2
      NNMAX=2*NMAX
 
!     Matrix inversion from LAPACK 
 
!      print *,'TT',INFO, NCHECK, NMAX
      DO I=1,NNMAX
       DO J=1,NNMAX
         ZQ(I,J)=DCMPLX(QR(I,J),QI(I,J))
       ENDDO
      ENDDO

      INFO=0
      CALL ZGETRF(NNMAX,NNMAX,ZQ,NPN2,IPIV,INFO)
!      print *,'TT2', NNMAX, NPN2, INFO
      IF (INFO.NE.0) WRITE (6,1100) INFO, NCHECK, NMAX
      CALL ZGETRI(NNMAX,ZQ,NPN2,IPIV,ZW,NPN2,INFO)
!      print *,'TT3', NNMAX, NPN2, INFO
      IF (INFO.NE.0) WRITE (6,1100) INFO, NCHECK, NMAX

 1100 FORMAT ('WARNING (TT subroutine):  info=', I2, I2, I2)
      DO I=1,NNMAX
         DO J=1,NNMAX
            TR=0D0
            TI=0D0
         DO K=1,NNMAX
                 ARR=RGQR(I,K)
                 ARI=RGQI(I,K)
                 AR=ZQ(K,J)
                 AI=DIMAG(ZQ(K,J))
                 TR=TR-ARR*AR+ARI*AI
                 TI=TI-ARR*AI-ARI*AR
         ENDDO
         TR1(I,J)=TR
         TI1(I,J)=TI
         ENDDO
      ENDDO
      RETURN
      END


!**********************************************************************
 
      SUBROUTINE CONST (NGAUSS,NMAX,MMAX,P,X,W,AN,ANN,S,SS,NP,EPS)
     
      use tmat_parameters
      use tmat_aux

      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 :: X(NPNG2),W(NPNG2),X1(NPNG1),W1(NPNG1)
      REAL*8 :: X2(NPNG1),W2(NPNG1)
      REAL*8 :: S(NPNG2),SS(NPNG2) 
      REAL*8 :: AN(NPN1),ANN(NPN1,NPN1),DD(NPN1)
 
      DO 10 N=1,NMAX
           NN=N*(N+1)
           AN(N)=DFLOAT(NN)
           D=DSQRT(DFLOAT(2*N+1)/DFLOAT(NN))
           DD(N)=D
           DO 10 N1=1,N
                DDD=D*DD(N1)*0.5D0
                ANN(N,N1)=DDD
                ANN(N1,N)=DDD
   10 CONTINUE
      NG=2*NGAUSS
      IF (NP.EQ.-2) GO  TO 11
      CALL GAUSS(NG,0,0,X,W)
      GO TO 19
   11 NG1=DFLOAT(NGAUSS)/2D0
      NG2=NGAUSS-NG1
      XX=-DCOS(DATAN(EPS))
      CALL GAUSS(NG1,0,0,X1,W1)
      CALL GAUSS(NG2,0,0,X2,W2)
      DO 12 I=1,NG1
         W(I)=0.5D0*(XX+1D0)*W1(I)
         X(I)=0.5D0*(XX+1D0)*X1(I)+0.5D0*(XX-1D0)
   12 CONTINUE
      DO 14 I=1,NG2
         W(I+NG1)=-0.5D0*XX*W2(I)
         X(I+NG1)=-0.5D0*XX*X2(I)+0.5D0*XX
   14 CONTINUE
      DO 16 I=1,NGAUSS
         W(NG-I+1)=W(I)
         X(NG-I+1)=-X(I)
   16 CONTINUE
   19 DO 20 I=1,NGAUSS
           Y=X(I)
           Y=1D0/(1D0-Y*Y)
           SS(I)=Y
           SS(NG-I+1)=Y
           Y=DSQRT(Y)
           S(I)=Y
           S(NG-I+1)=Y
   20 CONTINUE
      RETURN
      END

end module

