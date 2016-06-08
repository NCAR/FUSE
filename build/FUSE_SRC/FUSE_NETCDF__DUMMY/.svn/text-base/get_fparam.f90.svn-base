SUBROUTINE GET_FPARAM(NETCDF_FILE,IMOD,MPAR,XPAR)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Read parameters in LPARAM from the last parameter set in the specified NetCDF file
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE ddirectory                                        ! defines data directory
USE multiparam, ONLY: LPARAM, NUMPAR                  ! parameter names
IMPLICIT NONE
! input
CHARACTER(LEN=*), INTENT(IN)           :: NETCDF_FILE ! NetCDF file name
INTEGER(I4B), INTENT(IN)               :: IMOD        ! model index
INTEGER(I4B), INTENT(IN)               :: MPAR        ! number of model parameters
! internal
LOGICAL(LGT)                           :: LEXIST      ! .TRUE. if NetCDF file exists
INTEGER(I4B)                           :: IERR        ! error code
INTEGER(I4B)                           :: NCID        ! NetCDF file ID
INTEGER(I4B)                           :: IDIMID      ! NetCDF dimension ID
INTEGER(I4B)                           :: IVARID      ! NetCDF variable ID
INTEGER(I4B)                           :: IPAR        ! loop through model parameters
INTEGER(I4B)                           :: NPAR        ! number of parameter sets in output file
REAL(MSP)                              :: APAR        ! parameter value (single precision)
! output
REAL(SP), DIMENSION(MPAR), INTENT(OUT) :: XPAR        ! parameter value (whatever precision SP is)
include 'netcdf.inc'                                  ! use netCDF libraries
! ---------------------------------------------------------------------------------------
! check that the file exists
INQUIRE(FILE=TRIM(OUTPUT_PATH)//TRIM(NETCDF_FILE),EXIST=LEXIST)
IF (.NOT.LEXIST) THEN
 print *, ' NetCDF file defining the desired model does not exist '
 print *, '   File = ', TRIM(OUTPUT_PATH)//TRIM(NETCDF_FILE)
 STOP
ENDIF
! open file
IERR = NF_OPEN(TRIM(OUTPUT_PATH)//TRIM(NETCDF_FILE),NF_NOWRITE,NCID); CALL HANDLE_ERR(IERR)
 ! get number of parameter sets
 IERR = NF_INQ_DIMID(NCID,'par',IDIMID); CALL HANDLE_ERR(IERR)
 IERR = NF_INQ_DIMLEN(NCID,IDIMID,NPAR); CALL HANDLE_ERR(IERR)
 ! loop through parameters
 DO IPAR=1,NUMPAR
  ! get parameter value
  IERR = NF_INQ_VARID(NCID,TRIM(LPARAM(IPAR)%PARNAME),IVARID); CALL HANDLE_ERR(IERR)
  IERR = NF_GET_VAR1_REAL(NCID,IVARID,(/IMOD,NPAR/),APAR); CALL HANDLE_ERR(IERR)
  ! put parameter value in the output vector
  XPAR(IPAR) = APAR
 END DO
! close NetCDF file
IERR = NF_CLOSE(NCID)
! ---------------------------------------------------------------------------------------
END SUBROUTINE GET_FPARAM
