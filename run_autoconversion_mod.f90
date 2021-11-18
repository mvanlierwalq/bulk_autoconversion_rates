module run_autoconversion_mod

!implicit none
PRIVATE
PUBLIC    :: kk2000_autoconversion, sb2001_autoconversion, berry_reinhardt_autoconversion, &
             gammln, wgamma
!real :: qc,nc,qauto,nauto,qr,mvd

!real :: qauto_kk, nauto_kk, qauto_sb, nauto_sb, qauto_br, nauto_br

CONTAINS

!..................................................................................
! KK2000

  subroutine kk2000_autoconversion(qc,nc,auto,nauto,ncauto)

    implicit none
    real :: qc,nc,auto,nauto,ncauto
    real :: pi,rhow

    pi=3.14159
    rhow=1000.

    IF (QC.GE.1.E-6) THEN ! limit as applied in MORR kinematic code

        auto=1350.*QC**2.47*  &
           (NC/1.e6)**(-1.79) ! there should be a factor of air density rho here, but
! we can safely neglect this (implicit assumption rho = 1)

        nauto = auto/(4./3.*PI*RHOW*(25.E-6)**3)
        ncauto = auto/(qc/nc)
    end if

 end subroutine kk2000_autoconversion

!.................................................................................. 
! SB2001 autoconversion

  subroutine sb2001_autoconversion(qc,nc,qr,auto,nauto,ncauto)

    implicit none

    real :: qc,nc,auto,nauto,qr

    real, dimension(16) :: dnu
    real :: nu,pgam,dum,dum1,ncauto

    integer :: dumii

         dnu(1) = -0.557
         dnu(2) = -0.557
         dnu(3) = -0.430
         dnu(4) = -0.307
         dnu(5) = -0.186
         dnu(6) = -0.067
         dnu(7) = 0.050
         dnu(8) = 0.167
         dnu(9) = 0.282
         dnu(10) = 0.397
         dnu(11) = 0.512
         dnu(12) = 0.626
         dnu(13) = 0.739
         dnu(14) = 0.853
         dnu(15) = 0.966
         dnu(16) = 0.966

         PGAM=0.0005714*(NC/1.E6)+0.2714
         PGAM=1./(PGAM**2)-1.
         PGAM=MAX(PGAM,2.)
         PGAM=MIN(PGAM,10.)

          dumii=int(pgam)
         nu=dnu(dumii)+(dnu(dumii+1)-dnu(dumii))* &
               (pgam-real(dumii))

    IF (QC.GE.1.E-6) THEN ! limit as applied in MORR kinematic code

        dum = 1.-qc/(qc+qr)
        dum1 = 600.*dum**0.68*(1.-dum**0.68)**3

        auto = 9.44e9/(20.*2.6e-7)* &
        (nu+2.)*(nu+4.)/(nu+1.)**2* &
        (qc/1000.)**4/(nc/1.e6)**2* &
        (1.+dum1/(1.-dum)**2)*1000.    ! again, there should be factors of rho here which
! we will neglect, implicit assumption is rho = 1
        ncauto = auto*2./2.6e-7*1000.
        nauto = 0.5*ncauto

     end if

 end subroutine sb2001_autoconversion

!..................................................................................
! Berry and Reinhardt 1974 autoconversion
! Inputs are nc and rc (bulk number and mass)
! Note: rc is technically a mass concentration in the Thompson scheme (air density times mixing ratio) but
! we can ignore that here and just assume mass mixing ratio. Same with number concentration/number mixing ratio.

! Thus, we can use inputs of number and mass mixing ratio for nc and rc, units of kg^-1 and kg/kg, respectively.
! NOTE the latter is NOT M3, but bulk mass (so, factor of pi/6*rhow different from M3).

  subroutine berry_reinhardt_autoconversion(rc,nc,prr_wau,pnr_wau,pnc_wau)

    implicit none

    real, intent(in) :: rc,nc
    real, intent(out) :: prr_wau,pnr_wau, pnc_wau
    real :: pi,rho_w,am_r,bm_r,r1,obmr,d0r,d0c
    real :: nu_c,lamc,xdc
    real :: dc_g,dc_b,zeta1,zeta,taud,tau
    integer :: n,nu_ci

    !real :: wgamma,gammln

        REAL, DIMENSION(3,15) :: cce, ccg
        REAL, DIMENSION(15) ::  ocg1, ocg2

         PI = 3.1415926536
         rho_w = 1000.0
         am_r = PI*rho_w/6.0
         bm_r = 3.0


         r1 = 1.e-20 ! a very small number
 
         D0c = 1.E-6
         D0r = 50.E-6

         obmr = 1./bm_r

         do n=1,15
         cce(1,n) = n + 1.
         cce(2,n) = bm_r + n + 1.
         cce(3,n) = bm_r + n + 4.

         ccg(1,n) = WGAMMA(cce(1,n)) 
         ccg(2,n) = WGAMMA(cce(2,n))
         ccg(3,n) = WGAMMA(cce(3,n))

         ocg1(n) = 1./ccg(1,n)
         ocg2(n) = 1./ccg(2,n)
         end do

!         nu_c = MIN(15, NINT(1000.E6/nc) + 2)                                                                                                                 
!         lamc = (nc*am_r*ccg(2,nu_c)*ocg1(nu_c)/rc)**obmr
!         xDc = (bm_r + nu_c + 1.) / lamc
!            if (xDc.lt. D0c) then
!             lamc = cce(2,nu_c)/D0c
!            elseif (xDc.gt. D0r*2.) then
!             lamc = cce(2,nu_c)/(D0r*2.)
!            endif


          nu_c = MIN(15, NINT(1000.E6/nc) + 2)
          nu_ci = INT(nu_c)
          xDc = MAX(D0c*1.E6, ((rc/(am_r*nc))**obmr) * 1.E6)
          lamc = (nc*am_r* ccg(2,nu_ci) * ocg1(nu_ci) / rc)**obmr

          mvd_c = (3.0 + nu_c + 0.672)/lamc

! Now do actual auto conversion calculation

         if (rc.gt. 0.01e-3) then                                                                                                                                   
          Dc_g = ((ccg(3,nu_ci)*ocg2(nu_ci))**obmr / lamc) * 1.E6                                                                                                        
          Dc_b = (xDc*xDc*xDc*Dc_g*Dc_g*Dc_g - xDc*xDc*xDc*xDc*xDc*xDc) & 
                 **(1./6.)                                                                                                                            

          zeta1 = 0.5*((6.25E-6*xDc*Dc_b*Dc_b*Dc_b - 0.4) &                                                                                                            
                     + abs(6.25E-6*xDc*Dc_b*Dc_b*Dc_b - 0.4))                                                                                                  
          zeta = 0.027*rc*zeta1                                                                                                                                     
          taud = 0.5*((0.5*Dc_b - 7.5) + abs(0.5*Dc_b - 7.5)) + R1                                                                                                     
          tau  = 3.72/(rc*taud)
          prr_wau = zeta/tau
! Note: min statement below is only to conserve water within 1 time step of Thompson scheme
! This is not needed for our off-line box model tests
!          prr_wau = MIN(DBLE(rc*odts), prr_wau)
          pnr_wau = prr_wau / (am_r*nu_c*D0r*D0r*D0r)              ! RAIN2M       

          pnc_wau = prr_wau/(am_r*mvd_c*mvd_c*mvd_c)    !..ncloud auto

         endif

         end subroutine berry_reinhardt_autoconversion

      REAL FUNCTION WGAMMA(y)

      IMPLICIT NONE
      REAL, INTENT(IN):: y
      !real gammln

      WGAMMA = EXP(GAMMLN(y))

      END FUNCTION WGAMMA

!+---+-----------------------------------------------------------------+                                                                                                          
      REAL FUNCTION GAMMLN(XX)
!     --- RETURNS THE VALUE LN(GAMMA(XX)) FOR XX > 0.                                                                                                                             
      IMPLICIT NONE
      REAL, INTENT(IN):: XX
      DOUBLE PRECISION, PARAMETER:: STP = 2.5066282746310005D0
      DOUBLE PRECISION, DIMENSION(6), PARAMETER:: &
               COF = (/76.18009172947146D0, -86.50532032941677D0, &
                       24.01409824083091D0, -1.231739572450155D0, &
                      .1208650973866179D-2, -.5395239384953D-5/)
      DOUBLE PRECISION:: SER,TMP,X,Y
      INTEGER:: J

      X=XX
      Y=X
      TMP=X+5.5D0
      TMP=(X+0.5D0)*LOG(TMP)-TMP
      SER=1.000000000190015D0
      DO 11 J=1,6
        Y=Y+1.D0
        SER=SER+COF(J)/Y
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER/X)
      END FUNCTION GAMMLN
!  (C) Copr. 1986-92 Numerical Recipes Software 2.02 

end module run_autoconversion_mod
