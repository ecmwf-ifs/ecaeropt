
module tmatrixlib

      contains

      subroutine tmatrix( RAT, NDISTR, AXMAX, NPMAX, B, GAM, NKMAX, EPS, NP, LAM,    &
                & MRR, MRI, DDELT, NPNA, NDGS)

      use tmat_parameters
      use tmat_aux
      use tmat_CGcoeff
      use tmat_bessel
      use tmat_ct
      use tmat_tmat
      use tmat_gsp
      use tmat_procedures



      IMPLICIT REAL*8 (A-H,O-Z)

      REAL*8  :: LAM,MRR,MRI,X(NPNG2),W(NPNG2),S(NPNG2),SS(NPNG2)
      REAL*8  :: AN(NPN1),R(NPNG2),DR(NPNG2)
      REAL*8  :: DDR(NPNG2),DRR(NPNG2),DRI(NPNG2),ANN(NPN1,NPN1)
      REAL*8  :: XG(1000),WG(1000)
      REAL*8  :: ALPH1(NPL),ALPH2(NPL),ALPH3(NPL),ALPH4(NPL),BET1(NPL)
      REAL*8  :: BET2(NPL),XG1(2000),WG1(2000)
      REAL*8  :: AL1(NPL),AL2(NPL),AL3(NPL),AL4(NPL),BE1(NPL),BE2(NPL)

      P=DACOS(-1D0)
 
!  OPEN FILES *******************************************************
 
      OPEN (6,FILE='test')
      OPEN (10,FILE='tmatr.write')
 
!  INPUT DATA ********************************************************
 
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
 
      NCHECK=0
      IF (NP.EQ.-1.OR.NP.EQ.-2) NCHECK=1
      IF (NP.GT.0.AND.(-1)**NP.EQ.1) NCHECK=1
      WRITE (6,5454) NCHECK
 5454 FORMAT ('NCHECK=',I1)
      DAX=AXMAX/NPNAX
      !write (*,*) 'RAT',RAT
      IF (DABS(RAT-1D0).GT.1D-8.AND.NP.EQ.-1) CALL SAREA (EPS,RAT)
      if (DABS(RAT-1D0).GT.1D-8.AND.NP.GE.0)  CALL SURFCH(NP,EPS,RAT)
      IF (DABS(RAT-1D0).GT.1D-8.AND.NP.EQ.-2) CALL SAREAC (EPS,RAT)
      !PRINT 8000, RAT
 8000 FORMAT ('RAT=',F8.6)
      IF(NP.EQ.-1.AND.EPS.GE.1D0) PRINT 7000,EPS
      IF(NP.EQ.-1.AND.EPS.LT.1D0) PRINT 7001,EPS
      IF(NP.GE.0) PRINT 7100,NP,EPS
      IF(NP.EQ.-2.AND.EPS.GE.1D0) PRINT 7150,EPS
      IF(NP.EQ.-2.AND.EPS.LT.1D0) PRINT 7151,EPS
      PRINT 7400, LAM,MRR,MRI
      PRINT 7200, DDELT
 7000 FORMAT('RANDOMLY ORIENTED OBLATE SPHEROIDS, A/B=',F11.7)
 7001 FORMAT('RANDOMLY ORIENTED PROLATE SPHEROIDS, A/B=',F11.7)
 7100 FORMAT('RANDOMLY ORIENTED CHEBYSHEV PARTICLES, T',I1,'(',F5.2,')')
 7150 FORMAT('RANDOMLY ORIENTED OBLATE CYLINDERS, D/L=',F11.7)
 7151 FORMAT('RANDOMLY ORIENTED PROLATE CYLINDERS, D/L=',F11.7)
 7200 FORMAT ('ACCURACY OF COMPUTATIONS DDELT = ',D8.2)
 7400 FORMAT('LAM=',F10.6,3X,'MRR=',D10.4,3X,'MRI=',D10.4)
      DDELT=0.1D0*DDELT
      DO 600 IAX=1,NPNAX
         AXI=AXMAX-DAX*DFLOAT(IAX-1)
         R1=0.89031D0*AXI
         R2=1.56538D0*AXI
         NK=IDINT(AXI*NKMAX/AXMAX+2)                       
         IF (NK.GT.1000) PRINT 8001,NK
         IF (NK.GT.1000) STOP
         IF (NDISTR.EQ.3) CALL POWER (AXI,B,R1,R2)
 8001    FORMAT ('NK=',I4,' I.E., IS GREATER THAN 1000. ','EXECUTION TERMINATED.')
         CALL GAUSS (NK,0,0,XG,WG)
         Z1=(R2-R1)*0.5D0
         Z2=(R1+R2)*0.5D0
         Z3=R1*0.5D0
         IF (NDISTR.EQ.5) GO TO 3
         DO I=1,NK
            XG1(I)=Z1*XG(I)+Z2
            WG1(I)=WG(I)*Z1
         ENDDO
         GO TO 4
    3    DO I=1,NK
            XG1(I)=Z3*XG(I)+Z3
            WG1(I)=WG(I)*Z3
         ENDDO
         DO I=NK+1,2*NK
            II=I-NK
            XG1(I)=Z1*XG(II)+Z2
            WG1(I)=WG(II)*Z1
         ENDDO
         NK=NK*2
    4    CALL DISTRB (NK,XG1,WG1,NDISTR,AXI,B,GAM,R1,R2,REFF,VEFF,P)
         PRINT 8002,R1,R2
 8002    FORMAT('R1=',F10.6,'   R2=',F10.6)
         IF (DABS(RAT-1D0).LE.1D-6) PRINT 8003, REFF,VEFF
         IF (DABS(RAT-1D0).GT.1D-6) PRINT 8004, REFF,VEFF
 8003    FORMAT('EQUAL-VOLUME-SPHERE REFF=',F8.4,'   VEFF=',F7.4)
 8004    FORMAT('EQUAL-SURFACE-AREA-SPHERE REFF=',F8.4,'   VEFF=',F7.4)
         PRINT 7250,NK
 7250    FORMAT('NUMBER OF GAUSSIAN QUADRATURE POINTS ','IN SIZE AVERAGING =',I4)
         DO I=1,NPL
            ALPH1(I)=0D0
            ALPH2(I)=0D0
            ALPH3(I)=0D0
            ALPH4(I)=0D0
            BET1(I)=0D0
            BET2(I)=0D0
         ENDDO      
         CSCAT=0D0
         CEXTIN=0D0
         L1MAX=0
         DO 500 INK=1,NK
            I=NK-INK+1
            A=RAT*XG1(I)
            XEV=2D0*P*A/LAM
            IXXX=XEV+4.05D0*XEV**0.333333D0
            INM1=MAX0(4,IXXX)
            IF (INM1.GE.NPN1) PRINT 7333, NPN1
            IF (INM1.GE.NPN1) STOP
 7333 FORMAT('CONVERGENCE IS NOT OBTAINED FOR NPN1=',I3,'.  EXECUTION TERMINATED')
            QEXT1=0D0
            QSCA1=0D0
            DO 50 NMA=INM1,NPN1
               NMAX=NMA
               MMAX=1
               NGAUSS=NMAX*NDGS
               IF (NGAUSS.GT.NPNG1) PRINT 7340, NGAUSS
               IF (NGAUSS.GT.NPNG1) STOP
 7340          FORMAT('NGAUSS =',I3,' I.E. IS GREATER THAN NPNG1.','  EXECUTION TERMINATED')
 7334          FORMAT(' NMAX =', I3,'  DC2=',D8.2,'   DC1=',D8.2)
 7335 FORMAT('                              NMAX1 =', I3,'  DC2=',D8.2,'  DC1=',D8.2)
               CALL CONST(NGAUSS,NMAX,MMAX,P,X,W,AN,ANN,S,SS,NP,EPS)
               CALL VARY( LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,PPI,PIR,PII,R, &
                        & DR,DDR,DRR,DRI,NMAX)
               CALL TMATR0 ( NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,    &
                           &  DDR,DRR,DRI,NMAX,NCHECK)
               QEXT=0D0
               QSCA=0D0
               DO N=1,NMAX
                  N1=N+NMAX
                  TR1NN=TR1(N,N)
                  TI1NN=TI1(N,N)
                  TR1NN1=TR1(N1,N1)
                  TI1NN1=TI1(N1,N1)
                  DN1=DFLOAT(2*N+1)
                  QSCA=QSCA+DN1*(TR1NN*TR1NN+TI1NN*TI1NN+TR1NN1*TR1NN1+TI1NN1*TI1NN1)
                  QEXT=QEXT+(TR1NN+TR1NN1)*DN1
               ENDDO
               DSCA=DABS((QSCA1-QSCA)/QSCA)
               DEXT=DABS((QEXT1-QEXT)/QEXT)
!              PRINT 7334, NMAX,DSCA,DEXT
               QEXT1=QEXT
               QSCA1=QSCA
               NMIN=DFLOAT(NMAX)/2D0+1D0
               DO 10 N=NMIN,NMAX
                  N1=N+NMAX
                  TR1NN=TR1(N,N)
                  TI1NN=TI1(N,N)
                  TR1NN1=TR1(N1,N1)
                  TI1NN1=TI1(N1,N1)
                  DN1=DFLOAT(2*N+1)
                  DQSCA=DN1*(TR1NN*TR1NN+TI1NN*TI1NN+TR1NN1*TR1NN1+TI1NN1*TI1NN1)
                  DQEXT=(TR1NN+TR1NN1)*DN1
                  DQSCA=DABS(DQSCA/QSCA)
                  DQEXT=DABS(DQEXT/QEXT)
                  NMAX1=N
                  IF (DQSCA.LE.DDELT.AND.DQEXT.LE.DDELT) GO TO 12
   10          CONTINUE
   12          CONTINUE
!              PRINT 7335, NMAX1,DQSCA,DQEXT
               IF(DSCA.LE.DDELT.AND.DEXT.LE.DDELT) GO TO 55
               IF (NMA.EQ.NPN1) PRINT 7333, NPN1
               IF (NMA.EQ.NPN1) STOP      
   50       CONTINUE
   55       NNNGGG=NGAUSS+1
            IF (NGAUSS.EQ.NPNG1) PRINT 7336
            MMAX=NMAX1
            DO 150 NGAUS=NNNGGG,NPNG1
               NGAUSS=NGAUS
               NGGG=2*NGAUSS
 7336          FORMAT('WARNING: NGAUSS=NPNG1')
 7337          FORMAT(' NG=',I3,'  DC2=',D8.2,'   DC1=',D8.2)
               CALL CONST(NGAUSS,NMAX,MMAX,P,X,W,AN,ANN,S,SS,NP,EPS)
               CALL VARY( LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,PPI,PIR,PII,R, &
                        & DR,DDR,DRR,DRI,NMAX)
               CALL TMATR0 ( NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,  &
                           & DDR,DRR,DRI,NMAX,NCHECK)
               QEXT=0D0
               QSCA=0D0
               DO 104 N=1,NMAX
                  N1=N+NMAX
                  TR1NN=TR1(N,N)
                  TI1NN=TI1(N,N)
                  TR1NN1=TR1(N1,N1)
                  TI1NN1=TI1(N1,N1)
                  DN1=DFLOAT(2*N+1)
                  QSCA=QSCA+DN1*(TR1NN*TR1NN+TI1NN*TI1NN+TR1NN1*TR1NN1+TI1NN1*TI1NN1)
                  QEXT=QEXT+(TR1NN+TR1NN1)*DN1
  104          CONTINUE
               DSCA=DABS((QSCA1-QSCA)/QSCA)
               DEXT=DABS((QEXT1-QEXT)/QEXT)
!              PRINT 7337, NGGG,DSCA,DEXT
               QEXT1=QEXT
               QSCA1=QSCA
               IF(DSCA.LE.DDELT.AND.DEXT.LE.DDELT) GO TO 155
               IF (NGAUS.EQ.NPNG1) PRINT 7336
  150       CONTINUE
  155       CONTINUE
            QSCA=0D0
            QEXT=0D0
            NNM=NMAX*2
            DO 204 N=1,NNM
               QEXT=QEXT+TR1(N,N)
  204       CONTINUE
            IF (NMAX1.GT.NPN4) PRINT 7550, NMAX1
 7550       FORMAT ('NMAX1 = ',I3, ', i.e. greater than NPN4.',' Execution terminated')
            IF (NMAX1.GT.NPN4) STOP              
            DO 213 N2=1,NMAX1
               NN2=N2+NMAX
               DO 213 N1=1,NMAX1
                  NN1=N1+NMAX
                  ZZ1=TR1(N1,N2)
                  TR11(1,N1,N2)=ZZ1
                  ZZ2=TI1(N1,N2)
                  IT11(1,N1,N2)=ZZ2
                  ZZ3=TR1(N1,NN2)
                  TR12(1,N1,N2)=ZZ3
                  ZZ4=TI1(N1,NN2)
                  IT12(1,N1,N2)=ZZ4
                  ZZ5=TR1(NN1,N2)
                  TR21(1,N1,N2)=ZZ5
                  ZZ6=TI1(NN1,N2)
                  IT21(1,N1,N2)=ZZ6
                  ZZ7=TR1(NN1,NN2)
                  TR22(1,N1,N2)=ZZ7
                  ZZ8=TI1(NN1,NN2)
                  IT22(1,N1,N2)=ZZ8
                  QSCA=QSCA+ZZ1*ZZ1+ZZ2*ZZ2+ZZ3*ZZ3+ZZ4*ZZ4+ZZ5*ZZ5+ZZ6*ZZ6+ZZ7*ZZ7+ZZ8*ZZ8
  213       CONTINUE
!           PRINT 7800,0,DABS(QEXT),QSCA,NMAX
            DO 220 M=1,NMAX1
               CALL TMATR(M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,DDR,DRR,DRI,NMAX,NCHECK)
               NM=NMAX-M+1
               NM1=NMAX1-M+1
               M1=M+1
               QSC=0D0
               DO 214 N2=1,NM1
                  NN2=N2+M-1
                  N22=N2+NM
                  DO 214 N1=1,NM1
                     NN1=N1+M-1
                     N11=N1+NM
                     ZZ1=TR1(N1,N2)
                     TR11(M1,NN1,NN2)=ZZ1
                     ZZ2=TI1(N1,N2)
                     IT11(M1,NN1,NN2)=ZZ2
                     ZZ3=TR1(N1,N22)
                     TR12(M1,NN1,NN2)=ZZ3
                     ZZ4=TI1(N1,N22)
                     IT12(M1,NN1,NN2)=ZZ4
                     ZZ5=TR1(N11,N2)
                     TR21(M1,NN1,NN2)=ZZ5
                     ZZ6=TI1(N11,N2)
                     IT21(M1,NN1,NN2)=ZZ6
                     ZZ7=TR1(N11,N22)
                     TR22(M1,NN1,NN2)=ZZ7
                     ZZ8=TI1(N11,N22)
                     IT22(M1,NN1,NN2)=ZZ8
                     QSC=QSC+(ZZ1*ZZ1+ZZ2*ZZ2+ZZ3*ZZ3+ZZ4*ZZ4 + ZZ5*ZZ5+ZZ6*ZZ6+ZZ7*ZZ7+ZZ8*ZZ8)*2D0
  214          CONTINUE
               NNM=2*NM
               QXT=0D0
               DO 215 N=1,NNM
                  QXT=QXT+TR1(N,N)*2D0
  215          CONTINUE
               QSCA=QSCA+QSC
               QEXT=QEXT+QXT
!              PRINT 7800,M,DABS(QXT),QSC,NMAX
 7800          FORMAT(' m=',I3,'  qxt=',d12.6,'  qsc=',d12.6,'  nmax=',I3)
  220       CONTINUE
            COEFF1=LAM*LAM*0.5D0/P
            CSCA=QSCA*COEFF1
            CEXT=-QEXT*COEFF1
!           PRINT 7880, NMAX,NMAX1
 7880       FORMAT ('nmax=',I3,'   nmax1=',I3)
            CALL GSP (NMAX1,CSCA,LAM,AL1,AL2,AL3,AL4,BE1,BE2,LMAX)
            L1M=LMAX+1
            L1MAX=MAX(L1MAX,L1M)
            WGII=WG1(I)
            WGI=WGII*CSCA
            DO 250 L1=1,L1M
               ALPH1(L1)=ALPH1(L1)+AL1(L1)*WGI
               ALPH2(L1)=ALPH2(L1)+AL2(L1)*WGI
               ALPH3(L1)=ALPH3(L1)+AL3(L1)*WGI
               ALPH4(L1)=ALPH4(L1)+AL4(L1)*WGI
               BET1(L1)=BET1(L1)+BE1(L1)*WGI
               BET2(L1)=BET2(L1)+BE2(L1)*WGI
  250       CONTINUE
            CSCAT=CSCAT+WGI
            CEXTIN=CEXTIN+CEXT*WGII
!           PRINT 6070, I,NMAX,NMAX1,NGAUSS
 6070       FORMAT(4I6)
  500    CONTINUE
         DO 510 L1=1,L1MAX
            ALPH1(L1)=ALPH1(L1)/CSCAT
            ALPH2(L1)=ALPH2(L1)/CSCAT
            ALPH3(L1)=ALPH3(L1)/CSCAT
            ALPH4(L1)=ALPH4(L1)/CSCAT
            BET1(L1)=BET1(L1)/CSCAT
            BET2(L1)=BET2(L1)/CSCAT
  510    CONTINUE
         WALB=CSCAT/CEXTIN
         CALL HOVENR(L1MAX,ALPH1,ALPH2,ALPH3,ALPH4,BET1,BET2)
         ASYMM=ALPH1(2)/3D0
         PRINT 9100,CEXTIN,CSCAT,WALB,ASYMM
 9100    FORMAT('CEXT=',D12.6,2X,'CSCA=',D12.6,2X,2X,'W=',D12.6,2X,'<COS>=',D12.6)
         IF (WALB.GT.1D0) PRINT 9111
 9111    FORMAT ('WARNING: W IS GREATER THAN 1')
         WRITE (10,580) WALB,L1MAX
         DO L=1,L1MAX
            WRITE (10,575) ALPH1(L),ALPH2(L),ALPH3(L),ALPH4(L),BET1(L),BET2(L)           
         ENDDO   
  575    FORMAT(6D14.7)
  580    FORMAT(D14.8,I8)
         LMAX=L1MAX-1
         CALL MATR (ALPH1,ALPH2,ALPH3,ALPH4,BET1,BET2,LMAX,NPNA)
  600 CONTINUE
      ITIME=MCLOCK()
      TIME=DFLOAT(ITIME)/6000D0
      PRINT 1001,TIME
 1001 FORMAT (' time =',F8.2,' min')
      STOP
      END
 

end module

