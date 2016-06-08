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
USE fuse_fileManager,only:INPUT_PATH                  ! defines data directory
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
! ---------------------------------------------------------------------------------------
! CONTENT REMOVED FOR COPYRIGHT VIOLATION
! ---------------------------------------------------------------------------------------
END SUBROUTINE GET_OBJFNC
END MODULE GET_OBJFNC_MODULE
