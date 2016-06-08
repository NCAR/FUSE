SUBROUTINE HANDLE_ERR(IERR)
! Used to print our error statements from NetCDF calls and stop
USE nrtype
INTEGER(I4B)                             :: IERR    ! error code
include 'netcdf.inc'
IF (IERR.NE.NF_NOERR) THEN
 PRINT *, NF_STRERROR(IERR)
 STOP
ENDIF
END SUBROUTINE HANDLE_ERR
