SUBROUTINE PUT_PARAMS(IPAR,IMOD)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! write NetCDF output files  -- model parameters
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition structures (includes filename)
USE metaparams                                        ! metadata for model parameters
USE multistats, ONLY:MSTATS                           ! provide access to error message
USE parextract_module                                 ! extract parameters
IMPLICIT NONE
! input
INTEGER(I4B), INTENT(IN)               :: IPAR        ! parameter set index
INTEGER(I4B), INTENT(IN)               :: IMOD        ! model index
! internal
INTEGER(I4B)                           :: IERR,NCID   ! error code; NetCDF ID
INTEGER(I4B), DIMENSION(2)             :: INDX        ! indices for parameter write
INTEGER(I4B)                           :: IVAR        ! loop through parameters
REAL(SP)                               :: XPAR        ! desired parameter
REAL(MSP)                              :: APAR        ! convert to SP (need for SP write)
INTEGER(I4B)                           :: IVAR_ID     ! variable ID
INTEGER(I4B), PARAMETER                :: NDESC=8     ! number of model descriptors
INTEGER(I4B), PARAMETER                :: NCHAR=10    ! length of model descriptors
INTEGER(I4B), DIMENSION(3)             :: ISTART      ! starting position for array write
INTEGER(I4B), DIMENSION(3)             :: ICOUNT      ! count for array write
CHARACTER(LEN=10)                      :: TXTVEC      ! single model descriptor
! ---------------------------------------------------------------------------------------
! CONTENT REMOVED FOR COPYRIGHT VIOLATION
! ---------------------------------------------------------------------------------------
END SUBROUTINE PUT_PARAMS
