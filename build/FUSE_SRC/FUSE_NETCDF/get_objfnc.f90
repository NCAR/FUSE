MODULE GET_OBJFNC_MODULE
IMPLICIT NONE
CONTAINS
! ---------------------------------------------------------------------------------------
SUBROUTINE GET_OBJFNC(NETCDF_FILE,OF_NAME,IMOD,IPARSET,OF,XOPT)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Read data in variable "OF_NAME" from file "NETCDF_FILE"
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE multiparam, ONLY: LPARAM, NUMPAR                  ! parameter names
IMPLICIT NONE
! input
CHARACTER(LEN=*), INTENT(IN)           :: NETCDF_FILE ! NetCDF file name
CHARACTER(LEN=*), INTENT(IN)           :: OF_NAME     ! Objective function name
INTEGER(I4B), INTENT(IN)               :: IMOD        ! Model index
INTEGER(I4B), INTENT(IN)               :: IPARSET     ! index of the parameter set
! internal
LOGICAL(LGT)                           :: LEXIST      ! .TRUE. if NetCDF file exists
INTEGER(I4B)                           :: IERR        ! error code
INTEGER(I4B)                           :: NCID        ! NetCDF file ID
INTEGER(I4B)                           :: IDIMID      ! NetCDF dimension ID
INTEGER(I4B)                           :: IVARID      ! NetCDF variable ID
INTEGER(I4B)                           :: IPAR        ! loop through model parameters
REAL(MSP)                              :: OF_VAL      ! objective function value (single precision)
REAL(MSP)                              :: APAR        ! parameter value (single precision)
! output
REAL(SP), INTENT(OUT)                  :: OF          ! objective function value (whatever precision SP is)
REAL(SP), DIMENSION(:), INTENT(OUT)    :: XOPT        ! optimal parameter set
include 'netcdf.inc'                                  ! use netCDF libraries
! ---------------------------------------------------------------------------------------
! check that the file exists
INQUIRE(FILE=TRIM(NETCDF_FILE),EXIST=LEXIST)
IF (.NOT.LEXIST) THEN
 print *, ' NetCDF file defining the desired model does not exist '
 print *, '   File = ', TRIM(NETCDF_FILE)
 STOP
ENDIF
! open file
IERR = NF_OPEN(TRIM(NETCDF_FILE),NF_NOWRITE,NCID); CALL HANDLE_ERR(IERR)
 ! read objective function value
 IERR = NF_INQ_VARID(NCID,TRIM(OF_NAME),IVARID); CALL HANDLE_ERR(IERR)
 IERR = NF_GET_VAR1_REAL(NCID,IVARID,(/IMOD,IPARSET/),OF_VAL); CALL HANDLE_ERR(IERR)
 OF = OF_VAL  ! return objective function value IN whatever precision SP is
 ! loop through parameters
 DO IPAR=1,NUMPAR
  ! get parameter value
  IERR = NF_INQ_VARID(NCID,TRIM(LPARAM(IPAR)%PARNAME),IVARID); CALL HANDLE_ERR(IERR)
  IERR = NF_GET_VAR1_REAL(NCID,IVARID,(/IMOD,IPARSET/),APAR); CALL HANDLE_ERR(IERR)
  ! put parameter value in the output vector
  XOPT(IPAR) = APAR
 END DO
! close NetCDF file
IERR = NF_CLOSE(NCID)
! ---------------------------------------------------------------------------------------
END SUBROUTINE GET_OBJFNC
END MODULE GET_OBJFNC_MODULE
