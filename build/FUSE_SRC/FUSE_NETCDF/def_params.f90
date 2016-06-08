SUBROUTINE DEF_PARAMS(NMOD)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Define NetCDF output files -- parameter variables
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition (includes filename)
USE metaparams                                        ! metadata for all model parameters
USE multistats, ONLY:MSTATS                           ! model statistics structure
IMPLICIT NONE
! input
INTEGER(I4B), INTENT(IN)               :: NMOD        ! number of models
! internal
INTEGER(I4B)                           :: IERR,NCID   ! error code; NetCDF ID
INTEGER(I4B)                           :: NPAR_DIM    ! number of parameter sets
INTEGER(I4B)                           :: NMOD_DIM    ! number of models
INTEGER(I4B)                           :: NDIF_DIM    ! differences in models
INTEGER(I4B)                           :: NAME_DIM    ! length of string defining models
INTEGER(I4B)                           :: ERRM_DIM    ! length of string defining error message
INTEGER(I4B), DIMENSION(2)             :: FVAR        ! fixed dimensions
INTEGER(I4B), DIMENSION(3)             :: SVAR        ! model descriptor dimensions
INTEGER(I4B), DIMENSION(3)             :: EVAR        ! error message dimensions
INTEGER(I4B)                           :: IVAR        ! loop through variables
INTEGER(I4B)                           :: IVAR_ID     ! variable ID
include 'netcdf.inc'                                  ! use netCDF libraries
! ---------------------------------------------------------------------------------------
CALL PARDESCRIBE()               ! get list of parameter descriptions
! ---------------------------------------------------------------------------------------
! open file
IERR = NF_CREATE(TRIM(FNAME_NETCDF),NF_CLOBBER,NCID); CALL HANDLE_ERR(IERR)
 ! define dimensions
 IERR = NF_DEF_DIM(NCID,'mod',NMOD,NMOD_DIM); CALL HANDLE_ERR(IERR)
 IERR = NF_DEF_DIM(NCID,'par',NF_UNLIMITED,NPAR_DIM); CALL HANDLE_ERR(IERR)
 IERR = NF_DEF_DIM(NCID,'model_differences',8,NDIF_DIM); CALL HANDLE_ERR(IERR)
 IERR = NF_DEF_DIM(NCID,'model_name_length',10,NAME_DIM); CALL HANDLE_ERR(IERR)
 IERR = NF_DEF_DIM(NCID,'error_message_length',LEN(MSTATS%ERR_MESSAGE),ERRM_DIM)
 CALL HANDLE_ERR(IERR)
 ! assign dimensions to indices
 FVAR = (/NMOD_DIM,NPAR_DIM/)            ! dimensions for fixed output (parameters)
 SVAR = (/NAME_DIM,NDIF_DIM,NMOD_DIM/)   ! dimensions for model names
 EVAR = (/ERRM_DIM,NMOD_DIM,NPAR_DIM/)   ! dimensions for error messages
 ! define fixed output variables
 DO IVAR=1,NOUTPAR
  IERR = NF_DEF_VAR(NCID,TRIM(PNAME(IVAR)),NF_REAL,2,FVAR,IVAR_ID); CALL HANDLE_ERR(IERR)
  IERR = NF_PUT_ATT_TEXT(NCID,IVAR_ID,'long_name',LEN_TRIM(PDESC(IVAR)),TRIM(PDESC(IVAR)))
  CALL HANDLE_ERR(IERR)
  IERR = NF_PUT_ATT_TEXT(NCID,IVAR_ID,'units',LEN_TRIM(PUNIT(IVAR)),TRIM(PUNIT(IVAR)))
  CALL HANDLE_ERR(IERR)
  IERR = NF_PUT_ATT_REAL(NCID,IVAR_ID,'_FillValue',NF_REAL,1,-9999.); CALL HANDLE_ERR(IERR)
 END DO  ! ivar
 ! define model definitions
 IERR = NF_DEF_VAR(NCID,'model_description',NF_CHAR,3,SVAR,IVAR_ID); CALL HANDLE_ERR(IERR)
 ! define error messages
 IERR = NF_DEF_VAR(NCID,'error_message',NF_CHAR,3,EVAR,IVAR_ID); CALL HANDLE_ERR(IERR)
! end definitions and close file
IERR = NF_ENDDEF(NCID)
IERR = NF_CLOSE(NCID)
! ---------------------------------------------------------------------------------------
END SUBROUTINE DEF_PARAMS
