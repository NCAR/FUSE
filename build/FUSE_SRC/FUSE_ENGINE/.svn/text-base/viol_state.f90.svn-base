MODULE VIOL_STATE_MODULE
IMPLICIT NONE
CONTAINS
FUNCTION VIOL_STATE(X_TRY)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Identifies if the state vector is feasible
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition structures
USE model_defnames
USE multiparam                                        ! model parameters
USE multistate                                        ! model states (USE NSTATE)
USE model_numerix                                     ! model numerix
IMPLICIT NONE
! input
REAL(SP), DIMENSION(:), INTENT(IN)     :: X_TRY       ! vector of model states
! internal
REAL(SP)                               :: XMIN        ! very small number
INTEGER(I4B)                           :: ISTT        ! loop through model states
! output
LOGICAL(LGT)                           :: VIOL_STATE  ! .TRUE. if the states violate constraints
! ---------------------------------------------------------------------------------------
XMIN=FRACSTATE_MIN ! used to avoid zero derivatives
! ---------------------------------------------------------------------------------------
VIOL_STATE=.FALSE.
! loop through model states
DO ISTT=1,NSTATE
 SELECT CASE(CSTATE(ISTT)%iSNAME)
  ! upper tanks
  CASE (iopt_TENS1A)
   IF(X_TRY(ISTT).LT.XMIN*DPARAM%MAXTENS_1A) VIOL_STATE=.TRUE.
   IF(X_TRY(ISTT).GT.     DPARAM%MAXTENS_1A) VIOL_STATE=.TRUE. 
  CASE (iopt_TENS1B)
   IF(X_TRY(ISTT).LT.XMIN*DPARAM%MAXTENS_1B) VIOL_STATE=.TRUE.
   IF(X_TRY(ISTT).GT.     DPARAM%MAXTENS_1B) VIOL_STATE=.TRUE.
  CASE (iopt_TENS_1)
   IF(X_TRY(ISTT).LT.XMIN*DPARAM%MAXTENS_1)  VIOL_STATE=.TRUE.
   IF(X_TRY(ISTT).GT.     DPARAM%MAXTENS_1)  VIOL_STATE=.TRUE.
   print *, VIOL_STATE, desc_int2str(CSTATE(ISTT)%iSNAME), X_TRY(ISTT), XMIN*DPARAM%MAXTENS_1, DPARAM%MAXTENS_1
  CASE (iopt_FREE_1)
   IF(X_TRY(ISTT).LT.XMIN*DPARAM%MAXFREE_1)  VIOL_STATE=.TRUE.
   IF(X_TRY(ISTT).GT.     DPARAM%MAXFREE_1)  VIOL_STATE=.TRUE.
   print *, VIOL_STATE, desc_int2str(CSTATE(ISTT)%iSNAME), X_TRY(ISTT), XMIN*DPARAM%MAXFREE_1, DPARAM%MAXFREE_1
  CASE (iopt_WATR_1)
   IF(X_TRY(ISTT).LT.XMIN*MPARAM%MAXWATR_1)  VIOL_STATE=.TRUE.
   IF(X_TRY(ISTT).GT.     MPARAM%MAXWATR_1)  VIOL_STATE=.TRUE.
   print *, VIOL_STATE, desc_int2str(CSTATE(ISTT)%iSNAME), X_TRY(ISTT), XMIN*MPARAM%MAXWATR_1, MPARAM%MAXWATR_1
  ! lower tanks
  CASE (iopt_TENS_2)
   IF(X_TRY(ISTT).LT.XMIN*DPARAM%MAXTENS_2)  VIOL_STATE=.TRUE.
   IF(X_TRY(ISTT).GT.     DPARAM%MAXTENS_2)  VIOL_STATE=.TRUE.
   print *, VIOL_STATE, desc_int2str(CSTATE(ISTT)%iSNAME), X_TRY(ISTT), XMIN*DPARAM%MAXTENS_2, DPARAM%MAXTENS_2
  CASE (iopt_FREE2A)
   IF(X_TRY(ISTT).LT.XMIN*DPARAM%MAXFREE_2A) VIOL_STATE=.TRUE.
   IF(X_TRY(ISTT).GT.     DPARAM%MAXFREE_2A) VIOL_STATE=.TRUE.
   print *, VIOL_STATE, desc_int2str(CSTATE(ISTT)%iSNAME), X_TRY(ISTT), XMIN*DPARAM%MAXFREE_2A, DPARAM%MAXFREE_2A
  CASE (iopt_FREE2B)
   IF(X_TRY(ISTT).LT.XMIN*DPARAM%MAXFREE_2B) VIOL_STATE=.TRUE.
   IF(X_TRY(ISTT).GT.     DPARAM%MAXFREE_2B) VIOL_STATE=.TRUE.
   print *, VIOL_STATE, desc_int2str(CSTATE(ISTT)%iSNAME), X_TRY(ISTT), XMIN*DPARAM%MAXFREE_2B, DPARAM%MAXFREE_2B
  CASE (iopt_WATR_2)
   ! *** SET LOWER LIMITS ***
   IF (X_TRY(ISTT).LT.XMIN*MPARAM%MAXWATR_2) VIOL_STATE=.TRUE.
   ! *** SET UPPER LIMITS ***
   IF (SMODL%iARCH2.EQ.iopt_tens2pll_2 .OR. SMODL%iARCH2.EQ.iopt_fixedsiz_2) THEN
    ! cannot exceed capacity
    IF (X_TRY(ISTT).GT.MPARAM%MAXWATR_2)     VIOL_STATE=.TRUE.
   ELSE
    ! unlimited storage, but make sure the values are still sensible
    IF (X_TRY(ISTT).GT.MPARAM%MAXWATR_2*10._sp) VIOL_STATE=.TRUE.
   ENDIF
 END SELECT
END DO ! (loop through states)
! ---------------------------------------------------------------------------------------
END FUNCTION VIOL_STATE
END MODULE VIOL_STATE_MODULE
