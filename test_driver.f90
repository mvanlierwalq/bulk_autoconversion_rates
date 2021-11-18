PROGRAM test_driver
USE run_autoconversion_mod
IMPLICIT NONE

REAL :: qc, nc, qr, qauto, nauto, mvd

! auto and nauto rates below are mass and number rates for RAIN (not cloud, mass rate is equal in
! magnitude but opposite in sign for cloud and rain but magnitude of number rates differs in both
! sign and magnitude)

! specify cloud mass and number mixing ratio and rain mass mixing ratio (used in SB2001 only)
qc=1.e-3 ! mass mixing ratio in kg/kg
nc=100.e6 ! number mixing ratio in kg^-1

! calculate mean size, use this to filter out implausible qc and nc combinations
mvd = (6.*qc/(nc*1000.*3.14))**0.333333
print*,'input qc,nc,mean volume diam=',qc,nc,mvd

qr=0.1e-3 ! mass mixing ratio for rain (kg/kg) --? needed for SB2001 only
!qr=0.1e-3 ! mass mixing ratio for rain (kg/kg) --? needed for SB2001 only

! call berry and reinhardt autoconversion
call berry_reinhardt_autoconversion(qc,nc,qauto,nauto)

print*,'BR74',qauto,nauto

call kk2000_autoconversion(qc,nc,qauto,nauto)

print*,'KK2000',qauto,nauto

call sb2001_autoconversion(qc,nc,qr,qauto,nauto)

print*,'SB2001',qauto,nauto






end program test_driver
