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
LOGICAL(LGT)                           :: WRITE_VAR   ! used to denote if the variable is written
INTEGER(I4B)                           :: IERR,NCID   ! error code; NetCDF ID
INTEGER(I4B), DIMENSION(3)             :: INDX        ! indices for time series write
INTEGER(I4B)                           :: IVAR        ! loop through variables
REAL(SP)                               :: XVAR        ! desired variable (SP NOT NECESSARILY SP)
REAL(MSP)                              :: AVAR        ! desired variable (SINGLE PRECISION)
INTEGER(I4B)                           :: IVAR_ID     ! variable ID
include 'netcdf.inc'                                  ! use netCDF libraries
! ---------------------------------------------------------------------------------------
! open file
IERR = NF_OPEN(TRIM(FNAME_NETCDF),NF_WRITE,NCID); CALL HANDLE_ERR(IERR)

  ! define indices for model output
 INDX = (/ITIM,IMOD,IPAR/)

 ! loop through time-varying model output
 DO IVAR=1,NOUTVAR
  ! check if there is a need to write the variable
  IF (Q_ONLY) THEN
   WRITE_VAR=.FALSE.
   IF (TRIM(VNAME(IVAR)).EQ.'q_routed') WRITE_VAR=.TRUE.
   IF (TRIM(VNAME(IVAR)).EQ.'watr_1')   WRITE_VAR=.TRUE.
   IF (TRIM(VNAME(IVAR)).EQ.'watr_2')   WRITE_VAR=.TRUE.
   IF (.NOT.WRITE_VAR) CYCLE
  ENDIF
  ! write the variable
  XVAR = VAREXTRACT(VNAME(IVAR)); AVAR=XVAR                                  ! get variable ivar
  IERR = NF_INQ_VARID(NCID,TRIM(VNAME(IVAR)),IVAR_ID); CALL HANDLE_ERR(IERR) ! get variable ID
  IERR = NF_PUT_VAR1_REAL(NCID,IVAR_ID,INDX,AVAR); CALL HANDLE_ERR(IERR)     ! write data
 END DO  ! (ivar)
! close NetCDF file
IERR = NF_CLOSE(NCID)
! ---------------------------------------------------------------------------------------
END SUBROUTINE PUT_OUTPUT
