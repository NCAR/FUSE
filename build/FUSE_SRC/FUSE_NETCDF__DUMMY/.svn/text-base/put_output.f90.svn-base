SUBROUTINE PUT_OUTPUT(IPAR,IMOD,ITIM)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! write NetCDF output files
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition (includes filename)
USE metaoutput                                        ! metadata for time-varying model output
USE varextract_module                                 ! interface for the function to extract variables
IMPLICIT NONE
! input
INTEGER(I4B), INTENT(IN)               :: IPAR        ! parameter set index
INTEGER(I4B), INTENT(IN)               :: IMOD        ! model index
INTEGER(I4B), INTENT(IN)               :: ITIM        ! time step index
! internal
INTEGER(I4B)                           :: IERR,NCID   ! error code; NetCDF ID
INTEGER(I4B), DIMENSION(3)             :: INDX        ! indices for time series write
INTEGER(I4B)                           :: IVAR        ! loop through variables
REAL(SP)                               :: XVAR        ! desired variable (SP NOT NECESSARILY SP)
REAL(MSP)                              :: AVAR        ! desired variable (SINGLE PRECISION)
INTEGER(I4B)                           :: IVAR_ID     ! variable ID
! ---------------------------------------------------------------------------------------
! CONTENT REMOVED FOR COPYRIGHT VIOLATION
! ---------------------------------------------------------------------------------------
END SUBROUTINE PUT_OUTPUT
