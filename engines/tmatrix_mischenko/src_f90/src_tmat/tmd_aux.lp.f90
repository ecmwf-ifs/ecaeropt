module tmat_aux

  use iso_fortran_env
  use tmat_parameters
  implicit none

  contains

!*****************************************************************
 
      SUBROUTINE SAREA (D,RAT)

      implicit none
      real(kind=real64)     :: D
      real(kind=real64)     :: RAT
      real(kind=real64)     :: E, R, R0, R1

      !write (*,*) 'SAREA'
      IF (D.GE.1) then
         E=DSQRT(1D0-1D0/(D*D))
         R0 = D**(-4D0/3D0)*DLOG((1D0+E)/(1D0-E))
         R1 = 0.25D0*(2D0*D**(2D0/3D0) + R0/E)
         R=DSQRT(R1)
         RAT=1D0/R
      else
         E=DSQRT(1D0-D*D)
         R=0.5D0*(D**(2D0/3D0) + D**(-1D0/3D0)*DASIN(E)/E)
         R=DSQRT(R)
         RAT=1D0/R 
      endif

      END subroutine SAREA
 
!****************************************************************
 
      SUBROUTINE SURFCH(N,E,RAT)
      
      implicit none
      real(kind=real64)     :: E
      integer(kind=int32)   :: N
      real(kind=real64)     :: RAT


      integer, parameter    :: NG=60
      real(kind=real64)     :: DN, E2, EN, S, V, XI, DX, DXN
      real(kind=real64)     :: DS, DSN, DCN, A, A2, ENS, RS, RV
      real(kind=real64)     :: X(NG), W(NG)
      integer               :: I


      !write (*,*) 'SURFCH'
      
      DN=DFLOAT(N)
      E2=E*E
      EN=E*DN
      !NG=60
      CALL GAUSS (NG,0,0,X,W)
      S=0D0
      V=0D0
      DO I=1,NG
         XI=X(I)
         DX=DACOS(XI)
         DXN=DN*DX
         DS=DSIN(DX)
         DSN=DSIN(DXN)
         DCN=DCOS(DXN)
         A=1D0+E*DCN
         A2=A*A
         ENS=EN*DSN
         S=S+W(I)*A*DSQRT(A2+ENS*ENS)
         V=V+W(I)*(DS*A+XI*ENS)*DS*A2
      ENDDO
      RS=DSQRT(S*0.5D0)
      RV=(V*3D0/4D0)**(1D0/3D0)
      RAT=RV/RS
      
      END subroutine SURFCH

!********************************************************************
 
      SUBROUTINE SAREAC (EPS,RAT)

      implicit none
      real(kind=real64)     :: EPS
      real(kind=real64)     :: RAT
      real(kind=real64)     :: E, R, R0, R1
     
      RAT=(1.5D0/EPS)**(1D0/3D0)
      RAT=RAT/DSQRT( (EPS+2D0)/(2D0*EPS) )


      !write (*,*) 'SAREAC', RAT
      END subroutine SAREAC
 
!********************************************************************
!  COMPUTATION OF R1 AND R2 FOR A POWER LAW SIZE DISTRIBUTION WITH
!  EFFECTIVE RADIUS A AND EFFECTIVE VARIANCE B

      SUBROUTINE POWER(A,B,R1,R2)
      
      implicit none
      real(kind=real64), intent(in)  :: A, B
      real(kind=real64), intent(out) :: R1, R2
      real(kind=real64)              :: AX, BX, AA, BB
      
      !interface
      !      PURE function FUN(R, A, B) RESULT(Rout)
      !          use iso_fortran_env
      !          implicit none
      !          real(kind=real64) , intent(in) :: R, A, B
      !          real(kind=real64)              :: Rout
      !      end function FUN
      !end interface 

      !write(*,*) '--> POWER'
      AA=A
      BB=B
      AX=1D-5
      BX=A-1D-5
      R1=ZEROIN(FUN, AX,BX,0D0, AA, BB)
      R2=(1D0+B)*2D0*A-R1

 
      END

!***********************************************************************
! 
      PURE FUNCTION FUN(R1, A, B) RESULT(R3)

      implicit none
      real(kind=real64), intent(in) :: R1, A, B
      real(kind=real64)             :: R2, R3
      R2=(1D0+B)*2D0*A-R1
      R3=(R2-R1)/DLOG(R2/R1)-A

      END FUNCTION

!***********************************************************************
 
      REAL(KIND=real64) FUNCTION ZEROIN(FUN, AX,BX,TOL, AA, BB)
 
      implicit none
      real(kind=real64), intent(in) :: AX, BX, TOL, AA, BB
      real(kind=real64)             :: R,EPS, TOL1, A, B, FA, FB, C, FC, D, E, XM, S, P, Q
      
      interface
            PURE function FUN(R, A, B) RESULT(Rout)
                use iso_fortran_env
                implicit none
                real(kind=real64) , intent(in) :: R, A, B
                real(kind=real64)              :: Rout
            end function FUN
      end interface 

      !write(*,*) 'ZEROIN'
      EPS=1D0
   10 EPS=0.5D0*EPS
      TOL1=1D0+EPS
      IF (TOL1.GT.1D0) GO TO 10
   15 A=AX
      B=BX
      FA=FUN(A, AA, BB)
      FB=FUN(B, AA, BB)
      !write(*,*) 'FA, FB',FA,FB
   20 C=A
      FC=FA
      D=B-A
      E=D
   30 IF (DABS(FC).GE.DABS(FB)) GO TO 40
   35 A=B
      B=C
      C=A
      FA=FB
      FB=FC
      FC=FA
   40 TOL1=2D0*EPS*DABS(B)+0.5D0*TOL
      XM=0.5D0*(C-B)
      IF (DABS(XM).LE.TOL1) GO TO 90
   44 IF (FB.EQ.0D0) GO TO 90
   45 IF (DABS(E).LT.TOL1) GO TO 70
   46 IF (DABS(FA).LE.DABS(FB)) GO TO 70
   47 IF (A.NE.C) GO TO 50
   48 S=FB/FA
      P=2D0*XM*S
      Q=1D0-S
      GO TO 60
   50 Q=FA/FC
      R=FB/FC
      S=FB/FA
      P=S*(2D0*XM*Q*(Q-R)-(B-A)*(R-1D0))
      Q=(Q-1D0)*(R-1D0)*(S-1D0)
   60 IF (P.GT.0D0) Q=-Q
      P=DABS(P)
      IF ((2D0*P).GE.(3D0*XM*Q-DABS(TOL1*Q))) GO TO 70
   64 IF (P.GE.DABS(0.5D0*E*Q)) GO TO 70
   65 E=D
      D=P/Q
      GO TO 80
   70 D=XM
      E=D
   80 A=B
      FA=FB
      IF (DABS(D).GT.TOL1) B=B+D
      IF (DABS(D).LE.TOL1) B=B+DSIGN(TOL1,XM)
      FB=FUN(B, AA, BB)
      IF ((FB*(FC/DABS(FC))).GT.0D0) GO TO 20
   85 GO TO 30
   90 ZEROIN=B
      RETURN
      END
 
!**********************************************************************
!    CALCULATION OF POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE         *
!    FORMULA. IF IND1 = 0 - ON INTERVAL (-1,1), IF IND1 = 1 - ON      *
!    INTERVAL  (0,1). IF  IND2 = 1 RESULTS ARE PRINTED.               *
!    N - NUMBER OF POINTS                                             *
!    Z - DIVISION POINTS                                              *
!    W - WEIGHTS                                                      *
!**********************************************************************
 
      SUBROUTINE GAUSS ( N,IND1,IND2,Z,W )

      implicit none
      integer            :: N, IND1, IND2
      real(kind=real64)  :: Z(N), W(N)
      real(kind=real64)  :: A,B,C,F,CHECK, PB, PC, PA, DJ, X, ZZ
      integer            :: I,IND, M, NITER, K, J

      A=1D0
      B=2D0
      C=3D0
      IND=MOD(N,2)
      K=N/2+IND
      F=DFLOAT(N)
      DO 100 I=1,K
        M=N+1-I
        IF(I.EQ.1) X=A-B/((F+A)*F)
        IF(I.EQ.2) X=(Z(N)-A)*4D0+Z(N)
        IF(I.EQ.3) X=(Z(N-1)-Z(N))*1.6D0+Z(N-1)
        IF(I.GT.3) X=(Z(M+1)-Z(M+2))*C+Z(M+3)
        IF(I.EQ.K.AND.IND.EQ.1) X=0D0
        NITER=0
        CHECK=1D-16
  10    PB=1D0
        NITER=NITER+1
        IF (NITER.LE.100) GO TO 15
        CHECK=CHECK*10D0
  15    PC=X
        DJ=A
        DO J=2,N
            DJ=DJ+A
            PA=PB
            PB=PC
            PC=X*PB+(X*PB-PA)*(DJ-A)/DJ
        ENDDO
        PA=A/((PB-X*PC)*F)
        PB=PA*PC*(A-X*X)
        X=X-PB
        IF(DABS(PB).GT.CHECK*DABS(X)) GO TO 10
        Z(M)=X
        W(M)=PA*PA*(A-X*X)
        IF(IND1.EQ.0) W(M)=B*W(M)
        IF(I.EQ.K.AND.IND.EQ.1) GO TO 100
        Z(I)=-Z(M)
        W(I)=W(M)
  100 CONTINUE
      IF(IND2.NE.1) GO TO 110
!      PRINT 1100,N
! 1100 FORMAT(' ***  POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE FORMULA',
!     * ' OF ',I4,'-TH ORDER')
      DO  I=1,K
          ZZ=-Z(I)
 !         PRINT 1200,I,ZZ,I,W(I)
      ENDDO
! 1200 FORMAT(' ',4X,'X(',I4,') = ',F17.14,5X,'W(',I4,') = ',F17.14)
      GO TO 115
  110 CONTINUE
!C     PRINT 1300,N
! 1300 FORMAT(' GAUSSIAN QUADRATURE FORMULA OF ',I4,'-TH ORDER IS USED')
  115 CONTINUE
      IF(IND1.EQ.0) GO TO 140
      DO  I=1,N
          Z(I)=(A+Z(I))/B
      ENDDO
  140 CONTINUE
      RETURN
      END
 
 
 !*************************************************************************
 
      SUBROUTINE DISTRB (NNK,YY,WY,NDISTR,AA,BB,GAM,R1,R2,REFF,VEFF,PI)                                               
 
      implicit none
      integer, intent(in) :: NDISTR, NNK
      real(kind=real64),  intent(in) :: AA, BB, GAM, R1, R2
      real(kind=real64),  intent(out):: REFF, VEFF, PI
      real(kind=real64),  intent(in) :: YY(NNK)
      real(kind=real64),  intent(out):: WY(NNK)

      integer            :: I
      real(kind=real64)  :: A2, DB, X, Y, DA, B2, DAB, G, wysum, XI

      IF (NDISTR.EQ.1) THEN
        WY = modgamma(NNK, YY, AA, BB, GAM, WY)
      ENDIF
      IF (NDISTR.EQ.2) THEN
        WY = lognormal(NNK, YY, AA, BB, WY)
      ENDIF
      IF (NDISTR.EQ.3) THEN
        WY = hansentravis(NNK, YY, WY)
      ENDIF
      IF (NDISTR.EQ.4) THEN 
        WY = gammadis(NNK, YY, AA, BB, WY)
      ENDIF
      IF (NDISTR.EQ.5) THEN
        WY = modpowerlaw(NNK, YY, BB, R1, WY)
      ENDIF
                                                                 
      wysum=0D0
      DO I=1,NNK
         wysum=wysum+WY(I)
      enddo
      DO I=1,NNK
         WY(I)=WY(I)/wysum
      enddo
      G=0D0
      DO I=1,NNK
         X=YY(I)
         G=G+X*X*WY(I)
      enddo
      REFF=0D0
      DO I=1,NNK
         X=YY(I)
         REFF=REFF+X*X*X*WY(I)
      enddo
      REFF=REFF/G
      VEFF=0D0
      DO I=1,NNK
         X=YY(I)
         XI=X-REFF
         VEFF=VEFF+XI*XI*X*X*WY(I)
      enddo
      VEFF=VEFF/(G*REFF*REFF)
      RETURN


      ENDSUBROUTINE

      function modpowerlaw(NNK, YY, BB, R1, WY)

        implicit none
        integer, intent(in) :: NNK
        real(kind=real64),  intent(in) :: BB, R1
        real(kind=real64),  intent(in) :: YY(NNK), WY(NNK)
        real(kind=real64)              :: modpowerlaw(NNK)

        real(kind=real64)   :: X
        integer  :: I
      

        write(*, '(A,D10.4)') 'MODIFIED POWER LAW DISTRIBUTION,  alpha=',BB
        DO I=1,NNK
          X=YY(I)
          IF (X.LE.R1) THEN
            modpowerlaw(I)=WY(I)
          ELSE !(X.GT.R1) 
            modpowerlaw(I)=WY(I)*(X/R1)**BB
          ENDIF
        ENDDO

      end function

      function modgamma(NNK, YY, AA, BB, GAM, WY)

        implicit none
        integer, intent(in) :: NNK
        real(kind=real64),  intent(in) :: AA, BB, GAM
        real(kind=real64),  intent(in) :: YY(NNK), WY(NNK)
        real(kind=real64)              :: modgamma(NNK)

        real(kind=real64)   :: A2, DB, X, Y
        integer  :: I 
        
        write(*, '(A,F6.4,A,F6.4,A,F6.4)') 'MODIFIED GAMMA DISTRIBUTION, alpha=',AA, &
          & '  r_c=',BB,'  gamma=',GAM                                                   
        A2=AA/GAM                                                                 
        DB=1D0/BB
        DO I=1,NNK                                                             
          X=YY(I)                                                             
          Y=X**AA                                                                
          X=X*DB
          Y=Y*DEXP(-A2*(X**GAM))                                                 
          modgamma(I)=WY(I)*Y                                                       
        enddo

      endfunction


      function hansentravis(NNK, YY, WY)

        implicit none
        integer, intent(in) :: NNK
        real(kind=real64),  intent(in) :: YY(NNK), WY(NNK)
        real(kind=real64)              :: hansentravis(NNK)

        real(kind=real64)   :: DA, X
        integer  :: I
        
        write(*, '(A)') 'Hansen and Travis 1974 distribution'
      
        DO I=1,NNK                                                            
          X=YY(I)                                                                
          hansentravis(I)=WY(I)/(X*X*X)                                                 
        enddo                                                                  

      endfunction


      function gammadis(NNK, YY, AA, BB, WY)

        implicit none
        integer, intent(in) :: NNK
        real(kind=real64),  intent(in) :: AA, BB
        real(kind=real64),  intent(in) :: YY(NNK), WY(NNK)
        real(kind=real64)              :: gammadis(NNK)

        real(kind=real64)   :: DAB, B2, X
        integer  :: I
      

        write(*, '(A,F6.4,A,F6.4)') 'Gamma Distribution, a=',AA,'b=',BB                                                   
        B2=(1D0-3D0*BB)/BB                                                        
        DAB=1D0/(AA*BB)                                                          
        DO I=1,NNK                                                            
          X=YY(I)                                                                
          X=(X**B2)*DEXP(-X*DAB)                                                 
          gammadis(I)=WY(I)*X                                                       
        ENDDO 

      end function

      function lognormal(NNK, YY, AA, BB, WY)

        implicit none
        integer, intent(in) :: NNK
        real(kind=real64),  intent(in) :: AA, BB
        real(kind=real64),  intent(in) :: YY(NNK), WY(NNK)
        real(kind=real64)              :: lognormal(NNK)

        real(kind=real64)   :: DA, X, Y
        integer  :: I
        
        write(*, '(A,F6.4,A,F6.4)') 'Lognormal Distribution, r_g=',AA,'[ln(sigma_g)]**2=',BB                                                   
        DA=1D0/AA                                                                 
        DO I=1,NNK                                                            
          X=YY(I)                                                                
          Y=DLOG(X*DA)                                                          
          Y=DEXP(-Y*Y*0.5D0/BB)/X                                             
          lognormal(I)=WY(I)*Y                                                          
        enddo

      endfunction
                                                       
 

!*************************************************************
 
      SUBROUTINE HOVENR(L1,A1,A2,A3,A4,B1,B2)
      
      implicit none

      integer, intent(in) :: L1
      real(kind=real64),  intent(in) :: A1(L1), A2(L1), A3(L1)
      real(kind=real64),  intent(in) :: A4(L1), B1(L1), B2(L1)
      integer                        :: L, KONTR, LL, I
      real(kind=real64)              :: DL, DDL, AA1, AA2, AA3, AA4, BB1, BB2
      real(kind=real64)              :: CC, D, C, C1, C2, C3

      DO L=1,L1
         KONTR=1
         LL=L-1
         DL=DFLOAT(LL)*2D0+1D0
         DDL=DL*0.48D0
         AA1=A1(L)
         AA2=A2(L)
         AA3=A3(L)
         AA4=A4(L)
         BB1=B1(L)
         BB2=B2(L)
         IF(LL.GE.1.AND.DABS(AA1).GE.DL) KONTR=2
         IF(DABS(AA2).GE.DL) KONTR=2
         IF(DABS(AA3).GE.DL) KONTR=2
         IF(DABS(AA4).GE.DL) KONTR=2
         IF(DABS(BB1).GE.DDL) KONTR=2
         IF(DABS(BB2).GE.DDL) KONTR=2
         IF(KONTR.EQ.2) PRINT 3000,LL
         C=-0.1D0
         DO I=1,11
            C=C+0.1D0
            CC=C*C
            C1=CC*BB2*BB2
            C2=C*AA4
            C3=C*AA3
            IF((DL-C*AA1)*(DL-C*AA2)-CC*BB1*BB1.LE.-1D-4) KONTR=2
            IF((DL-C2)*(DL-C3)+C1.LE.-1D-4) KONTR=2
            IF((DL+C2)*(DL-C3)-C1.LE.-1D-4) KONTR=2
            IF((DL-C2)*(DL+C3)-C1.LE.-1D-4) KONTR=2
            IF(KONTR.EQ.2) PRINT 4000,LL,C
         ENDDO
      ENDDO

      IF(KONTR.EQ.1) PRINT 2000
 2000 FORMAT('TEST OF VAN DER MEE & HOVENIER IS SATISFIED')
 3000 FORMAT('TEST OF VAN DER MEE & HOVENIER IS NOT SATISFIED, L=',I3)
 4000 FORMAT('TEST OF VAN DER MEE & HOVENIER IS NOT SATISFIED, L=',I3,'   A=',D9.2)
      RETURN
      END
 


!****************************************************************
!
!    CALCULATION OF THE SCATTERING MATRIX FOR GIVEN EXPANSION
!    COEFFICIENTS
!
!    A1,...,B2 - EXPANSION COEFFICIENTS
!    LMAX - NUMBER OF COEFFICIENTS MINUS 1
!    N - NUMBER OF SCATTERING ANGLES
!        THE CORRESPONDING SCATTERING ANGLES ARE GIVEN BY
!        180*(I-1)/(N-1) (DEGREES), WHERE I NUMBERS THE ANGLES
 
      SUBROUTINE MATR(A1,A2,A3,A4,B1,B2,LMAX,NPNA)

      use tmat_parameters
      implicit none

      real(kind=real64)      :: A1(NPL),A2(NPL),A3(NPL),A4(NPL),B1(NPL),B2(NPL)
      integer                :: LMAX, NPNA
      integer                :: N, L1MAX, L, L1, I1
      real(kind=real64)      :: DN, DA, DB, TB, TAA, D6, U, DL, DL1
      real(kind=real64)      :: F11, F2, F3, F44, F12, F34, F22, F33
      real(kind=real64)      :: P1, P2, P3, P4, PP1, PP2, PP3, PP4
      real(kind=real64)      :: P, PL1, PL2, PL3, PL4

      character(len=90)      :: fmt1, fmt2

      fmt1='(3X,A,6X,A,6X,A,6X,A,6X,A,6X,A,6X,A)'
      fmt2='(5X,A,8X,A,8X,A,8X,A,8X,A,8X,A,8X,A)'
      
      N=NPNA
      DN=1D0/DFLOAT(N-1)
      DA=DACOS(-1D0)*DN
      DB=180D0*DN
      L1MAX=LMAX+1
      write(*,*) ' '
      write(*,fmt1) 'S','ALPHA1','ALPHA2','ALPHA3','ALPHA4','BETA1','BETA2'
      DO L1=1,L1MAX
         L=L1-1
         write(*,'(I3,6F12.5)') L,A1(L1),A2(L1),A3(L1),A4(L1),B1(L1),B2(L1)
      ENDDO
      TB=-DB
      TAA=-DA
      write(*,*)
      write(*, fmt2) '<','F11','F22','F33','F44','F12','F34'
      D6=DSQRT(6D0)*0.25D0
      
      !N=NPNA
      !DN=1D0/DFLOAT(N-1)
      !DA=DACOS(-1D0)*DN
      !DB=180D0*DN
      !L1MAX=LMAX+1
      !PRINT 1000
      !1000 FORMAT(' ')
      !PRINT 1001
      !1001 FORMAT(' ',2X,'S',6X,'ALPHA1',6X,'ALPHA2',6X,'ALPHA3',6X,'ALPHA4',7X,'BETA1',7X,'BETA2')
      !DO 10 L1=1,L1MAX
      !   L=L1-1
      !   PRINT 1002,L,A1(L1),A2(L1),A3(L1),A4(L1),B1(L1),B2(L1)
      !10 CONTINUE
      !1002 FORMAT(' ',I3,6F12.5)
      !TB=-DB
      !TAA=-DA
      !PRINT 1000
      !PRINT 1003
      !1003 FORMAT(' ',5X,'<',8X,'F11',8X,'F22',8X,'F33',8X,'F44', 8X,'F12',8X,'F34')

      D6=DSQRT(6D0)*0.25D0
      DO I1=1,N
         TAA=TAA+DA
         TB=TB+DB
         U=DCOS(TAA)
         F11=0D0
         F2=0D0
         F3=0D0
         F44=0D0
         F12=0D0
         F34=0D0
         P1=0D0
         P2=0D0
         P3=0D0
         P4=0D0
         PP1=1D0
         PP2=0.25D0*(1D0+U)*(1D0+U)
         PP3=0.25D0*(1D0-U)*(1D0-U)
         PP4=D6*(U*U-1D0)
         DO L1=1,L1MAX
            L=L1-1
            DL=DFLOAT(L)
            DL1=DFLOAT(L1)
            F11=F11+A1(L1)*PP1
            F44=F44+A4(L1)*PP1
            IF(L.NE.LMAX) then
               PL1=DFLOAT(2*L+1)
               P=(PL1*U*PP1-DL*P1)/DL1
               P1=PP1
               PP1=P
            ENDIF 
            IF(L.LT.2) CYCLE

            F2=F2+(A2(L1)+A3(L1))*PP2
            F3=F3+(A2(L1)-A3(L1))*PP3
            F12=F12+B1(L1)*PP4
            F34=F34+B2(L1)*PP4
            
            IF(L.EQ.LMAX) CYCLE
            
            PL2=DFLOAT(L*L1)*U
            PL3=DFLOAT(L1*(L*L-4))
            PL4=1D0/DFLOAT(L*(L1*L1-4))
            P=(PL1*(PL2-4D0)*PP2-PL3*P2)*PL4
            P2=PP2
            PP2=P
            P=(PL1*(PL2+4D0)*PP3-PL3*P3)*PL4
            P3=PP3
            PP3=P
            P=(PL1*U*PP4-DSQRT(DFLOAT(L*L-4))*P4)/DSQRT(DFLOAT(L1*L1-4))
            P4=PP4
            PP4=P
         ENDDO
         F22=(F2+F3)*0.5D0
         F33=(F2-F3)*0.5D0
!        F22=F22/F11
!        F33=F33/F11
!        F44=F44/F11
!        F12=-F12/F11
!        F34=F34/F11
         write(*, '(F6.2, 6F11.4)') TB,F11,F22,F33,F44,F12,F34
      ENDDO
      RETURN
      END


end module
