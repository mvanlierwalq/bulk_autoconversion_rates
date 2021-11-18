PROGRAM driver
USE run_autoconversion_mod
USE netcdf
IMPLICIT NONE

!REAL :: qc, nc, qr, qauto, nauto, mvd
REAL, ALLOCATABLE, DIMENSION(:)      :: qc, nc
REAL, ALLOCATABLE, DIMENSION(:,:)  :: mvd, &
                                      qauto_br, nauto_br, &
                                      qauto_kk, nauto_kk, &
                                      qauto_sb, nauto_sb

REAL        :: qc_max, qc_min, nc_max, nc_min, qr, qstep, nstep
!..
INTEGER                   :: mcid, qdimid, ndimid, &
                             qc_varid, nc_varid, &
                             mvd_varid, &
                             qauto_br_varid, nauto_br_varid, &
                             qauto_sb_varid, nauto_sb_varid, &
                             qauto_kk_varid, nauto_kk_varid, &
                             i,j

INTEGER, PARAMETER :: nresq=200, nresn=200, spacy=1

ALLOCATE( qc(nresq) )
ALLOCATE( nc(nresn) )
ALLOCATE( qauto_sb(nresn,nresq) )
ALLOCATE( nauto_sb(nresn,nresq) )
ALLOCATE( qauto_kk(nresn,nresq) )
ALLOCATE( nauto_kk(nresn,nresq) )
ALLOCATE( qauto_br(nresn,nresq) )
ALLOCATE( nauto_br(nresn,nresq) )
ALLOCATE( mvd(nresn,nresq)   )

! auto and nauto rates below are mass and number rates for RAIN (not cloud, mass rate is equal in
! magnitude but opposite in sign for cloud and rain but magnitude of number rates differs in both
! sign and magnitude)

! specify cloud mass and number mixing ratio and rain mass mixing ratio (used in SB2001 only)
qc_max = 1.e-2 ! mass mixing ratio in kg/kg
qc_min = 1.e-4 ! mass mixing ratio in kg/kg
nc_max = 100.e8 ! number mixing ratio in kg^-1
nc_min = 100.e5 ! number mixing ratio in kg^-1

qr=0.1e-3 ! mass mixing ratio for rain (kg/kg) --? needed for SB2001 only
!qr=0.1e-3 ! mass mixing ratio for rain (kg/kg) --? needed for SB2001 only


!..If we're doing linear spacing
IF (spacy == 0) THEN
  qstep = (qc_max-qc_min)/REAL(nresq)
  nstep = (nc_max-nc_min)/REAL(nresn)
  DO i = 1,nresq
    qc(i) = qc_min + (i-1)*(qc_max-qc_min)/REAL(nresq)
  ENDDO
  DO i = 1,nresn
    nc(i) = nc_min + (i+1)*(nc_max-nc_min)/REAL(nresn)
  ENDDO
!..If we're doing logarithmic spacing
ELSEIF (spacy == 1) THEN
  qstep = (LOG10(qc_max) - LOG10(qc_min))/REAL(nresq)
  nstep = (LOG10(nc_max) - LOG10(nc_min))/REAL(nresn)
  DO i = 1,nresq
    qc(i) = qc_min*10.**(qstep*(i-1))
  ENDDO
  DO i = 1,nresn
    nc(i) = nc_min*10.**(nstep*(i-1))
  ENDDO
ENDIF


DO i = 1,nresq   !..q loop
  DO j = 1,nresn   !..n loop
    ! calculate mean size, use this to filter out implausible qc and nc combinations
    mvd(j,i) = (6.*qc(i)/(nc(j)*1000.*3.14))**0.333333
    CALL berry_reinhardt_autoconversion(qc(i),nc(j),qauto_br(j,i),nauto_br(j,i))
    CALL kk2000_autoconversion(qc(i), nc(j), qauto_kk(j,i),nauto_kk(j,i))
    CALL sb2001_autoconversion(qc(i), nc(j), qr, qauto_sb(j,i), nauto_sb(j,i))
  ENDDO
ENDDO



!..Output to NetCDF
CALL check(nf90_create( path='autoconversion_grid_01.nc', cmode=nf90_clobber, ncid=mcid))
!..
CALL check(nf90_def_dim( mcid, 'qc', nresq, qdimid))
CALL check(nf90_def_dim( mcid, 'nc', nresn, ndimid))
!..
CALL check(nf90_def_var( mcid, 'qc', NF90_FLOAT, qdimid, qc_varid))
CALL check(nf90_def_var( mcid, 'nc', NF90_FLOAT, ndimid, nc_varid))
CALL check(nf90_def_var( mcid, 'qauto_br', NF90_FLOAT, [ndimid, qdimid], qauto_br_varid))
CALL check(nf90_def_var( mcid, 'nauto_br', NF90_FLOAT, [ndimid, qdimid], nauto_br_varid))
CALL check(nf90_def_var( mcid, 'qauto_sb', NF90_FLOAT, [ndimid, qdimid], qauto_sb_varid))
CALL check(nf90_def_var( mcid, 'nauto_sb', NF90_FLOAT, [ndimid, qdimid], nauto_sb_varid))
CALL check(nf90_def_var( mcid, 'qauto_kk', NF90_FLOAT, [ndimid, qdimid], qauto_kk_varid))
CALL check(nf90_def_var( mcid, 'nauto_kk', NF90_FLOAT, [ndimid, qdimid], nauto_kk_varid))
CALL check(nf90_def_var( mcid, 'mvd', NF90_FLOAT, [ndimid,qdimid], mvd_varid))
!..
CALL check(nf90_enddef( mcid ))
!..
CALL check(nf90_put_var( mcid, qc_varid, qc))
CALL check(nf90_put_var( mcid, nc_varid, nc))
CALL check(nf90_put_var( mcid, qauto_br_varid, qauto_br) )
CALL check(nf90_put_var( mcid, nauto_br_varid, nauto_br) )
CALL check(nf90_put_var( mcid, qauto_kk_varid, qauto_kk) )
CALL check(nf90_put_var( mcid, nauto_kk_varid, nauto_kk) )
CALL check(nf90_put_var( mcid, qauto_sb_varid, qauto_sb) )
CALL check(nf90_put_var( mcid, nauto_sb_varid, nauto_sb) )
CALL check(nf90_put_var( mcid, mvd_varid, mvd) )
!..
CALL check(nf90_close( mcid ))


DEALLOCATE( qc )
DEALLOCATE( nc )
DEALLOCATE( qauto_kk )
DEALLOCATE( nauto_kk )
DEALLOCATE( qauto_sb )
DEALLOCATE( nauto_sb )
DEALLOCATE( qauto_br )
DEALLOCATE( nauto_br )
DEALLOCATE( mvd)

CONTAINS
!....................................................
SUBROUTINE check(status)
USE netcdf
INTEGER, INTENT(IN) :: status

IF (status /= nf90_noerr) THEN
  PRINT*, TRIM(nf90_strerror(status))
  STOP "Stopped"
END IF
END SUBROUTINE check


end program driver
