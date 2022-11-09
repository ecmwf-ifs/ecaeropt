subroutine interp_ri(nin,zlambtab, z_n_r_tab, z_n_i_tab, nout, lambda_int,znr,zni,verbose)

USE PARKIND1 , ONLY : JPIM, JPRB

implicit none

!-- interpolation of refractive indices on output wavelengths

REAL(KIND=JPRB),INTENT(IN)     :: zlambtab(nin), z_n_r_tab(nin), z_n_i_tab(nin)
REAL(KIND=JPRB),INTENT(IN)     :: lambda_int(nout)
INTEGER(KIND=JPIM),INTENT(IN)  :: nin, nout, verbose

REAL(KIND=JPRB),INTENT(OUT) :: znr(nout), zni(nout)
REAL(KIND=JPRB)    :: zint
INTEGER(KIND=JPIM) :: jwl,jwv


!-- interpolate the refractive index
if (verbose.eq.1) then
   PRINT *,'Interpolating the refractive index '
end if
DO jwv=1, nout
  IF (lambda_int(jwv) <  zlambtab(1)) THEN
    znr(jwv)=z_n_r_tab(1)
    zni(jwv)=z_n_i_tab(1)
  ELSEIF (lambda_int(jwv) >  zlambtab(nin)) THEN
    znr(jwv)=z_n_r_tab(nin)
    zni(jwv)=z_n_i_tab(nin)
  ELSEIF (lambda_int(jwv) >= zlambtab(1) .AND. lambda_int(jwv) <= zlambtab(nin)) THEN
    DO jwl=1,nin-1
      IF (lambda_int(jwv) >= zlambtab(jwl) .AND. lambda_int(jwv) <= zlambtab(jwl+1)) THEN
        zint=(lambda_int(jwv)-zlambtab(jwl))/ (zlambtab(jwl+1)-zlambtab(jwl))
        znr(jwv)=z_n_r_tab(jwl)+ zint * (z_n_r_tab(jwl+1)-z_n_r_tab(jwl))
        zni(jwv)=z_n_i_tab(jwl)+ zint * (z_n_i_tab(jwl+1)-z_n_i_tab(jwl))
       ENDIF
     ENDDO
   ENDIF
   if (verbose.eq.1) then
      PRINT 9002,jwv,lambda_int(jwv)*1.E+06, znr(jwv),zni(jwv)
   end if
9002 format(1x,I3,2F10.6,F10.6)
ENDDO


end subroutine
