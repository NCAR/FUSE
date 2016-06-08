SUBROUTINE GET_SMODEL(NETCDF_FILE,IMOD)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Read model decisions from a NetCDF output file
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE model_defn -- populate structure SMODL
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! data types, etc.
USE fuse_fileManager, only : OUTPUT_PATH              ! define output path
USE model_defn                                        ! model definition structures
IMPLICIT NONE
! input
CHARACTER(LEN=*), INTENT(IN)           :: NETCDF_FILE ! NetCDF file name
INTEGER(I4B), INTENT(IN)               :: IMOD        ! model index
! internal
LOGICAL(LGT)                           :: LEXIST      ! .TRUE. if NetCDF file exists
INTEGER(I4B)                           :: IERR        ! error code
INTEGER(I4B)                           :: NCID        ! NetCDF file ID
INTEGER(I4B)                           :: IDIMID      ! NetCDF dimension ID
INTEGER(I4B)                           :: IVARID      ! NetCDF variable ID
INTEGER(I4B)                           :: NDESC       ! number of model descriptors
INTEGER(I4B)                           :: NCHAR       ! length of model descriptors
INTEGER(I4B)                           :: IDESC       ! loop thru model descriptors
INTEGER(I4B), DIMENSION(3)             :: ISTART      ! start indices for data read
INTEGER(I4B), DIMENSION(3)             :: ICOUNT      ! number of elements read in each dimension
CHARACTER(LEN=50)                      :: TXTVEC      ! text vector
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
 ! get number of model differences
 IERR = NF_INQ_DIMID(NCID,'model_differences',IDIMID); CALL HANDLE_ERR(IERR)
 IERR = NF_INQ_DIMLEN(NCID,IDIMID,NDESC); CALL HANDLE_ERR(IERR)
 ! get length of model description
 IERR = NF_INQ_DIMID(NCID,'model_name_length',IDIMID); CALL HANDLE_ERR(IERR)
 IERR = NF_INQ_DIMLEN(NCID,IDIMID,NCHAR); CALL HANDLE_ERR(IERR)
 IF (LEN(TXTVEC).LT.NCHAR) STOP ' text vector is not long enough to hold model description '
 ! get variable name
 IERR = NF_INQ_VARID(NCID,'model_description',IVARID); CALL HANDLE_ERR(IERR)
 ! loop through descriptors
 DO IDESC=1,NDESC
  ! define indices
  ISTART = (/    1,IDESC,IMOD/)   ! starting position of array
  ICOUNT = (/NCHAR,    1,   1/)   ! number of array elements (one descriptor, one model)
  ! extract text string
  IERR = NF_GET_VARA_TEXT(NCID,IVARID,ISTART,ICOUNT,TXTVEC); CALL HANDLE_ERR(IERR)
  ! put text string in structure
  IF (IDESC.EQ.1) SMODL%RFERR(1:NCHAR) = TXTVEC(1:NCHAR)
  IF (IDESC.EQ.2) SMODL%ARCH1(1:NCHAR) = TXTVEC(1:NCHAR) 
  IF (IDESC.EQ.3) SMODL%ARCH2(1:NCHAR) = TXTVEC(1:NCHAR) 
  IF (IDESC.EQ.4) SMODL%QSURF(1:NCHAR) = TXTVEC(1:NCHAR) 
  IF (IDESC.EQ.5) SMODL%QPERC(1:NCHAR) = TXTVEC(1:NCHAR) 
  IF (IDESC.EQ.6) SMODL%ESOIL(1:NCHAR) = TXTVEC(1:NCHAR) 
  IF (IDESC.EQ.7) SMODL%QINTF(1:NCHAR) = TXTVEC(1:NCHAR) 
  IF (IDESC.EQ.8) SMODL%Q_TDH(1:NCHAR) = TXTVEC(1:NCHAR) 
  print *, TXTVEC(1:NCHAR) 
 END DO
IERR = NF_CLOSE(NCID); CALL HANDLE_ERR(IERR)
! ---------------------------------------------------------------------------------------
END SUBROUTINE GET_SMODEL
