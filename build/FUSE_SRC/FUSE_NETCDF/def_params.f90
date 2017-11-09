SUBROUTINE DEF_PARAMS(NPAR)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! Modified by Nans Addor to include snow routine
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Define NetCDF output files -- parameter variables
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition (includes filename)
USE metaparams                                        ! metadata for all model parameters
USE multistats, ONLY: MSTATS                          ! model statistics structure
USE multistate, only: ncid_out                        ! NetCDF output file ID
IMPLICIT NONE
! input
INTEGER(I4B), INTENT(IN)               :: NPAR        ! number of parameter sets
! internal
INTEGER(I4B)                           :: IERR        ! error code
INTEGER(I4B)                           :: NPAR_DIM    ! number of parameter sets
INTEGER(I4B)                           :: NMOD_DIM    ! number of models
INTEGER(I4B)                           :: NDIF_DIM    ! differences in models
INTEGER(I4B)                           :: NAME_DIM    ! length of string defining models
INTEGER(I4B)                           :: ERRM_DIM    ! length of string defining error message
INTEGER(I4B), DIMENSION(1)             :: FVAR        ! fixed dimensions
INTEGER(I4B), DIMENSION(3)             :: SVAR        ! model descriptor dimensions
INTEGER(I4B), DIMENSION(3)             :: EVAR        ! error message dimensions
INTEGER(I4B)                           :: IVAR        ! loop through variables
INTEGER(I4B)                           :: IVAR_ID     ! variable ID
include 'netcdf.inc'                                  ! use netCDF libraries
! ---------------------------------------------------------------------------------------
CALL PARDESCRIBE()               ! get list of parameter descriptions
! ---------------------------------------------------------------------------------------
PRINT *, 'Define NetCDF output files - parameter variables = ', TRIM(FNAME_NETCDF_PARA)
! Create file
IERR = NF_CREATE(TRIM(FNAME_NETCDF_PARA),NF_CLOBBER,ncid_out); CALL HANDLE_ERR(IERR)
 ! define dimensions
 ! IERR = NF_DEF_DIM(ncid_out,'mod',NMOD,NMOD_DIM); CALL HANDLE_ERR(IERR)
! IERR = NF_DEF_DIM(ncid_out,'par',NF_UNLIMITED,NPAR_DIM); CALL HANDLE_ERR(IERR)
 IERR = NF_DEF_DIM(ncid_out,'par',NPAR,NPAR_DIM); CALL HANDLE_ERR(IERR) ! TODO : max number of parameter - should not be hard-coded
 !IERR = NF_DEF_DIM(ncid_out,'model_differences',9,NDIF_DIM); CALL HANDLE_ERR(IERR) !TODO: this should not be hard-coded
 !IERR = NF_DEF_DIM(ncid_out,'model_name_length',10,NAME_DIM); CALL HANDLE_ERR(IERR)
 !IERR = NF_DEF_DIM(ncid_out,'error_message_length',LEN(MSTATS%ERR_MESSAGE),ERRM_DIM)
 ! assign dimensions to indices
 FVAR = (/NPAR_DIM/)            ! dimensions for fixed output (parameters)
 !SVAR = (/NAME_DIM,NDIF_DIM,NMOD_DIM/)   ! dimensions for model names
 !EVAR = (/ERRM_DIM,NMOD_DIM,NPAR_DIM/)   ! dimensions for error messages
 ! define fixed output variables
 DO IVAR=1,NOUTPAR
  IERR = NF_DEF_VAR(ncid_out,TRIM(PNAME(IVAR)),NF_REAL,1,FVAR,IVAR_ID); CALL HANDLE_ERR(IERR)
  IERR = NF_PUT_ATT_TEXT(ncid_out,IVAR_ID,'long_name',LEN_TRIM(PDESC(IVAR)),TRIM(PDESC(IVAR)))
  CALL HANDLE_ERR(IERR)
  IERR = NF_PUT_ATT_TEXT(ncid_out,IVAR_ID,'units',LEN_TRIM(PUNIT(IVAR)),TRIM(PUNIT(IVAR)))
  CALL HANDLE_ERR(IERR)
  IERR = NF_PUT_ATT_REAL(ncid_out,IVAR_ID,'_FillValue',NF_REAL,1,-9999.); CALL HANDLE_ERR(IERR)
 END DO  ! ivar
  ! define model definitions
 !IERR = NF_DEF_VAR(ncid_out,'model_description',NF_CHAR,3,SVAR,IVAR_ID); CALL HANDLE_ERR(IERR)
  ! define error messages
 !IERR = NF_DEF_VAR(ncid_out,'error_message',NF_CHAR,3,EVAR,IVAR_ID); CALL HANDLE_ERR(IERR)
! end definitions and close file
IERR = NF_ENDDEF(ncid_out)
IERR = NF_CLOSE(ncid_out)
! ---------------------------------------------------------------------------------------
END SUBROUTINE DEF_PARAMS
