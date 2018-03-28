module get_fparam_module
USE nrtype
USE netcdf
implicit none
private
public::GET_PRE_PARAM,GET_SCE_PARAM

contains

SUBROUTINE GET_PRE_PARAM(NETCDF_FILE,ISET,IMOD,MPAR,XPAR)

! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Nans Addor, 2017 - Based on GET_SCE_PARAM
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Load a pre-defined set of parameter values from a NetCDF file
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE fuse_fileManager, only : OUTPUT_PATH              ! define output path
USE multiparam, ONLY: LPARAM, NUMPAR                  ! parameter names
IMPLICIT NONE
! input
CHARACTER(LEN=*), INTENT(IN)           :: NETCDF_FILE ! NetCDF file name
INTEGER(I4B), INTENT(IN)               :: ISET        ! indice of parameter set to extract
INTEGER(I4B), INTENT(IN)               :: IMOD        ! model index
INTEGER(I4B), INTENT(IN)               :: MPAR        ! number of model parameters
! internal
INTEGER(I4B), DIMENSION(1)             :: INDX        ! indices for parameter extraction
LOGICAL(LGT)                           :: LEXIST      ! .TRUE. if NetCDF file exists
INTEGER(I4B)                           :: IERR        ! error code
INTEGER(I4B)                           :: NCID        ! NetCDF file ID
INTEGER(I4B)                           :: IDIMID      ! NetCDF dimension ID
INTEGER(I4B)                           :: IVARID      ! NetCDF variable ID
INTEGER(I4B)                           :: IPAR        ! loop through model parameters
INTEGER(I4B)                           :: NPAR        ! number of parameter sets in output file
REAL(DP)                               :: APAR        ! parameter value (single precision)
! output
REAL(SP), DIMENSION(MPAR), INTENT(OUT) :: XPAR        ! parameter value (whatever precision SP is)
include 'netcdf.inc'                                  ! use netCDF libraries
! ---------------------------------------------------------------------------------------
! check that the file exists
INQUIRE(FILE=TRIM(NETCDF_FILE),EXIST=LEXIST)
IF (.NOT.LEXIST) THEN
 print *, 'The NetCDF file containing the predefined parameter set does not exist:'
 print *, TRIM(NETCDF_FILE)
 STOP
ENDIF

print *, 'Opening parameter file:', TRIM(NETCDF_FILE)

! open file
IERR = NF_OPEN(TRIM(NETCDF_FILE),NF_NOWRITE,NCID); CALL HANDLE_ERR(IERR)

! get number of parameter sets
IERR = NF_INQ_DIMID(NCID,'par',IDIMID); CALL HANDLE_ERR(IERR)
IERR = NF_INQ_DIMLEN(NCID,IDIMID,NPAR); CALL HANDLE_ERR(IERR)

print *, 'NPAR - total number of parameter sets in parameter file:', NPAR

IF (ISET.GT.NPAR) THEN
 print *, 'Impossible to extract parameter set', ISET, 'since there are only', NPAR, 'parameter sets'
 STOP
ENDIF

print *, 'Extracting parameter set', ISET

! loop through parameters
DO IPAR=1,NUMPAR

  ! get parameter id
  IERR = NF_INQ_VARID(NCID,TRIM(LPARAM(IPAR)%PARNAME),IVARID); CALL HANDLE_ERR(IERR)

  ! get parameter value for the selected parameter set
  INDX = (/ISET/)
  IERR = NF_GET_VAR1_DOUBLE(NCID,IVARID,INDX,APAR); CALL HANDLE_ERR(IERR)

  ! put parameter value in the output vector
  XPAR(IPAR) = APAR

  print *, 'PARAM VALUES:',LPARAM(IPAR)%PARNAME, '->', APAR

END DO

PRINT *, 'Predefined parameter set loaded into XPAR!'

! close NetCDF file
IERR = NF_CLOSE(NCID)
! ---------------------------------------------------------------------------------------
END SUBROUTINE GET_PRE_PARAM

SUBROUTINE GET_SCE_PARAM(NETCDF_FILE,IMOD,MPAR,XPAR)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! Modified by Nans Addor - previously the parameter set extracted was the last one,
! now it is the one associated with the lowest RMSE
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Read parameters in LPARAM from the parameter set with the lowest RMSE in the specified NetCDF file
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
INTEGER(I4B), DIMENSION(1)             :: INDX        ! indices for parameter extraction
LOGICAL(LGT)                           :: LEXIST      ! .TRUE. if NetCDF file exists
INTEGER(I4B)                           :: IERR        ! error code
INTEGER(I4B)                           :: NCID        ! NetCDF file ID
INTEGER(I4B)                           :: IDIMID      ! NetCDF dimension ID
INTEGER(I4B)                           :: IVARID      ! NetCDF variable ID
INTEGER(I4B)                           :: I_RAW_RMSE  ! NetCDF RMSE ID
INTEGER(I4B), DIMENSION(1)             :: I_OPT_PARA  ! index of the optimum parameter set (e.g. lowest RSME) - MUST BE DIMENSIONED
REAL(DP), DIMENSION(:),ALLOCATABLE     :: RAW_RMSE    ! RMSE for each parameter set
REAL(DP), DIMENSION(1)                 :: LOWEST_RAW_RMSE    ! LOWEST RAW RMSE
INTEGER(I4B)                           :: IPAR        ! loop through model parameters
INTEGER(I4B)                           :: NPAR        ! number of parameter sets in output file
REAL(DP)                               :: APAR        ! parameter value (single precision)
! output
REAL(SP), DIMENSION(MPAR), INTENT(OUT) :: XPAR        ! parameter value (whatever precision SP is)
include 'netcdf.inc'                                  ! use netCDF libraries
! ---------------------------------------------------------------------------------------
! check that the file exists
INQUIRE(FILE=TRIM(NETCDF_FILE),EXIST=LEXIST)
IF (.NOT.LEXIST) THEN
 print *, 'The NetCDF file containing the SCE parameter sets does not exist:'
 print *, TRIM(NETCDF_FILE)
 STOP
ELSE

  print *, 'Loading SCE parameters from'
  print *, TRIM(NETCDF_FILE)

ENDIF

! open file
IERR = NF_OPEN(TRIM(NETCDF_FILE),NF_NOWRITE,NCID); CALL HANDLE_ERR(IERR)

 ! get number of parameter sets
 IERR = NF_INQ_DIMID(NCID,'par',IDIMID); CALL HANDLE_ERR(IERR)
 IERR = NF_INQ_DIMLEN(NCID,IDIMID,NPAR); CALL HANDLE_ERR(IERR)

 ! extract RMSE for each parameter set
 print *, 'Length of the par dimension (the number of parameter sets produced by SCE is lower)', NPAR

 ALLOCATE(RAW_RMSE(NPAR),STAT=IERR); IF(IERR.NE.0) STOP ' problem allocating space for RAW_RMSE '

 IERR = NF_INQ_VARID(NCID,'raw_rmse',I_RAW_RMSE); CALL HANDLE_ERR(IERR)
 IERR = NF_GET_VAR_DOUBLE(NCID,I_RAW_RMSE,RAW_RMSE); CALL HANDLE_ERR(IERR)

 I_OPT_PARA = MINLOC(RAW_RMSE,DIM=1) !TODO: use argument MASK to find best parameter set for each of the SCE run
 LOWEST_RAW_RMSE=RAW_RMSE(I_OPT_PARA)
 print *, 'Index of parameter set with lowest RMSE =',I_OPT_PARA
 print *, 'Lowest RMSE =',LOWEST_RAW_RMSE

 DEALLOCATE(RAW_RMSE,STAT=IERR); IF (IERR.NE.0) STOP ' problem deallocating ATIME/TDATA '

 PRINT *, 'Reading from NetCDF file parameter values for best parameter set:'

 ! loop through parameters
 DO IPAR=1,NUMPAR

  ! get parameter id
  IERR = NF_INQ_VARID(NCID,TRIM(LPARAM(IPAR)%PARNAME),IVARID); CALL HANDLE_ERR(IERR)

  ! get parameter value for the optimal parameter set
  INDX = (/I_OPT_PARA/)
  IERR = NF_GET_VAR1_DOUBLE(NCID,IVARID,INDX,APAR); CALL HANDLE_ERR(IERR)

  ! put parameter value in the output vector
  XPAR(IPAR) = APAR

  print *, 'PARAM VALUES:',LPARAM(IPAR)%PARNAME, '->', APAR

 END DO

 PRINT *, 'Best parameter set loaded into XPAR!'

! close NetCDF file
IERR = NF_CLOSE(NCID)
! ---------------------------------------------------------------------------------------
END SUBROUTINE GET_SCE_PARAM

end module get_fparam_module
