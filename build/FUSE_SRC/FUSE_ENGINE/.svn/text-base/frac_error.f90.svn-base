module frac_error_mod
implicit none
contains
FUNCTION FRAC_ERROR(X_END1,X_END2)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Calculates the fractional error in each state (relative to state capacity)
!  for one-step and two step implicit solutions
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definitions
USE model_defnames
USE multiparam                                        ! model parameters
IMPLICIT NONE
! input/output
REAL(SP), DIMENSION(:), INTENT(IN)     :: X_END1      ! one-step solution
REAL(SP), DIMENSION(:), INTENT(IN)     :: X_END2      ! two-step solution
REAL(SP), DIMENSION(SIZE(X_END1))      :: FRAC_ERROR  ! fractional error
! internal
INTEGER(I4B)                           :: ISTT        ! loop through model states
! ---------------------------------------------------------------------------------------
DO ISTT=1,NSTATE
 IF (CSTATE(ISTT)%iSNAME.EQ.iopt_TENS1A) FRAC_ERROR(ISTT) = ABS(X_END1(ISTT)-X_END2(ISTT))/DPARAM%MAXTENS_1A
 IF (CSTATE(ISTT)%iSNAME.EQ.iopt_TENS1B) FRAC_ERROR(ISTT) = ABS(X_END1(ISTT)-X_END2(ISTT))/DPARAM%MAXTENS_1B 
 IF (CSTATE(ISTT)%iSNAME.EQ.iopt_TENS_1) FRAC_ERROR(ISTT) = ABS(X_END1(ISTT)-X_END2(ISTT))/DPARAM%MAXTENS_1
 IF (CSTATE(ISTT)%iSNAME.EQ.iopt_FREE_1) FRAC_ERROR(ISTT) = ABS(X_END1(ISTT)-X_END2(ISTT))/DPARAM%MAXFREE_1
 IF (CSTATE(ISTT)%iSNAME.EQ.iopt_WATR_1) FRAC_ERROR(ISTT) = ABS(X_END1(ISTT)-X_END2(ISTT))/MPARAM%MAXWATR_1
 IF (CSTATE(ISTT)%iSNAME.EQ.iopt_TENS_2) FRAC_ERROR(ISTT) = ABS(X_END1(ISTT)-X_END2(ISTT))/DPARAM%MAXTENS_2
 IF (CSTATE(ISTT)%iSNAME.EQ.iopt_FREE2A) FRAC_ERROR(ISTT) = ABS(X_END1(ISTT)-X_END2(ISTT))/DPARAM%MAXFREE_2A
 IF (CSTATE(ISTT)%iSNAME.EQ.iopt_FREE2B) FRAC_ERROR(ISTT) = ABS(X_END1(ISTT)-X_END2(ISTT))/DPARAM%MAXFREE_2B
 IF (CSTATE(ISTT)%iSNAME.EQ.iopt_WATR_2) FRAC_ERROR(ISTT) = ABS(X_END1(ISTT)-X_END2(ISTT))/MPARAM%MAXWATR_2
END DO
! ---------------------------------------------------------------------------------------
END FUNCTION FRAC_ERROR
endmodule frac_error_mod
