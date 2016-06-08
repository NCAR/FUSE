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
! ---------------------------------------------------------------------------------------
! CONTENT REMOVED FOR COPYRIGHT VIOLATION
! ---------------------------------------------------------------------------------------
END SUBROUTINE DEF_PARAMS
