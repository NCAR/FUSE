SUBROUTINE PUT_SSTATS(IPAR,IMOD)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! write NetCDF output files -- summary statistics
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition structures (includes filename)
USE meta_stats                                        ! metadata for summary statistics
USE multistats                                        ! model summary statistics
USE model_numerix                                     ! model numerix parameters and arrays
USE sumextract_module                                 ! module to extract summary statistics
IMPLICIT NONE
! input
INTEGER(I4B), INTENT(IN)               :: IPAR        ! parameter set index
INTEGER(I4B), INTENT(IN)               :: IMOD        ! model index
! internal
INTEGER(I4B)                           :: IERR,NCID   ! error code; NetCDF ID
INTEGER(I4B), DIMENSION(2)             :: INDX        ! indices for parameter write
INTEGER(I4B)                           :: IVAR        ! loop through parameters
REAL(SP)                               :: XPAR        ! desired parameter (SP may not be SP)
REAL(MSP)                              :: APAR        ! desired parameter (...but MSP is SP)
INTEGER(I4B)                           :: IVAR_ID     ! variable ID
! ---------------------------------------------------------------------------------------
! CONTENT REMOVED FOR COPYRIGHT VIOLATION
! ---------------------------------------------------------------------------------------
END SUBROUTINE PUT_SSTATS
