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
USE fuse_fileManager, only : OUTPUT_PATH              ! define output path
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
INTEGER(I4B)                           :: I_RAW_RMSE  ! NetCDF RMSE ID
INTEGER(I4B), DIMENSION(1)             :: I_OPT_PARA  ! index of the optimum parameter set (e.g. lowest RSME)
REAL(DP), DIMENSION(:),ALLOCATABLE     :: RAW_RMSE    ! RMSE for each parameter set
INTEGER(I4B)                           :: IPAR        ! loop through model parameters
INTEGER(I4B)                           :: NPAR        ! number of parameter sets in output file
REAL(DP)                               :: APAR        ! parameter value (single precision)
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

 ! extract RMSE for each parameter set
 ALLOCATE(RAW_RMSE(NPAR),STAT=IERR); IF(IERR.NE.0) STOP ' problem allocating space for RAW_RMSE '
 
 IERR = NF_INQ_VARID(NCID,'raw_rmse',I_RAW_RMSE); CALL HANDLE_ERR(IERR)
 IERR = NF_GET_VAR_DOUBLE(NCID,I_RAW_RMSE,RAW_RMSE); CALL HANDLE_ERR(IERR)

 I_OPT_PARA = MINLOC(RAW_RMSE)

 print *, 'INDEX OF LOWEST RMSE', I_OPT_PARA,'RMSE',RAW_RMSE(I_OPT_PARA)
 
  DEALLOCATE(RAW_RMSE,STAT=IERR); IF (IERR.NE.0) STOP ' problem deallocating ATIME/TDATA '

 ! loop through parameters
 DO IPAR=1,NUMPAR

  ! get parameter id
  IERR = NF_INQ_VARID(NCID,TRIM(LPARAM(IPAR)%PARNAME),IVARID); CALL HANDLE_ERR(IERR)

  ! get parameter value for the optimal parameter set
  IERR = NF_GET_VAR1_DOUBLE(NCID,IVARID,(/IMOD,I_OPT_PARA/),APAR); CALL HANDLE_ERR(IERR) 
  
  ! put parameter value in the output vector
  XPAR(IPAR) = APAR

 END DO

! close NetCDF file
IERR = NF_CLOSE(NCID)
! ---------------------------------------------------------------------------------------
END SUBROUTINE GET_FPARAM
