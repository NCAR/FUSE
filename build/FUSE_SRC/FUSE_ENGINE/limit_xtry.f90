MODULE LIMIT_XTRY_MODULE
IMPLICIT NONE
CONTAINS
SUBROUTINE LIMIT_XTRY(X_TRY)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Imposes constraints on the vector of model states
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition structures
USE model_defnames
USE multiparam                                        ! model parameters
USE multistate                                        ! model states (USE NSTATE)
USE model_numerix                                     ! model numerix
IMPLICIT NONE
! input/output
REAL(SP), DIMENSION(:), INTENT(INOUT)  :: X_TRY       ! vector of model states
! internal
REAL(SP)                               :: XMIN        ! very small number
INTEGER(I4B)                           :: ISTT        ! loop through model states
! ---------------------------------------------------------------------------------------
XMIN=FRACSTATE_MIN ! used to avoid zero derivatives
! ---------------------------------------------------------------------------------------
! loop through model states
DO ISTT=1,NSTATE
 SELECT CASE(CSTATE(ISTT)%iSNAME)
  ! upper tanks
  CASE (iopt_TENS1A)
    IF(X_TRY(ISTT).LT.XMIN*DPARAM%MAXTENS_1A) X_TRY(ISTT) = XMIN*DPARAM%MAXTENS_1A
    IF(X_TRY(ISTT).GT.     DPARAM%MAXTENS_1A) X_TRY(ISTT) =      DPARAM%MAXTENS_1A
  CASE (iopt_TENS1B)
    IF(X_TRY(ISTT).LT.XMIN*DPARAM%MAXTENS_1B) X_TRY(ISTT) = XMIN*DPARAM%MAXTENS_1B
    IF(X_TRY(ISTT).GT.     DPARAM%MAXTENS_1B) X_TRY(ISTT) =      DPARAM%MAXTENS_1B
  CASE (iopt_TENS_1)
    IF(X_TRY(ISTT).LT.XMIN*DPARAM%MAXTENS_1)  X_TRY(ISTT) = XMIN*DPARAM%MAXTENS_1
    IF(X_TRY(ISTT).GT.     DPARAM%MAXTENS_1)  X_TRY(ISTT) =      DPARAM%MAXTENS_1
  CASE (iopt_FREE_1)
    IF(X_TRY(ISTT).LT.XMIN*DPARAM%MAXFREE_1)  X_TRY(ISTT) = XMIN*DPARAM%MAXFREE_1
    IF(X_TRY(ISTT).GT.     DPARAM%MAXFREE_1)  X_TRY(ISTT) =      DPARAM%MAXFREE_1
  CASE (iopt_WATR_1)
    IF(X_TRY(ISTT).LT.XMIN*MPARAM%MAXWATR_1)  X_TRY(ISTT) = XMIN*MPARAM%MAXWATR_1
    IF(X_TRY(ISTT).GT.     MPARAM%MAXWATR_1)  X_TRY(ISTT) =      MPARAM%MAXWATR_1
  ! lower tanks
  CASE (iopt_TENS_2)
    IF(X_TRY(ISTT).LT.XMIN*DPARAM%MAXTENS_2)  X_TRY(ISTT) = XMIN*DPARAM%MAXTENS_2
    IF(X_TRY(ISTT).GT.     DPARAM%MAXTENS_2)  X_TRY(ISTT) =      DPARAM%MAXTENS_2
  CASE (iopt_FREE2A)
    IF(X_TRY(ISTT).LT.XMIN*DPARAM%MAXFREE_2A) X_TRY(ISTT) = XMIN*DPARAM%MAXFREE_2A
    IF(X_TRY(ISTT).GT.     DPARAM%MAXFREE_2A) X_TRY(ISTT) =      DPARAM%MAXFREE_2A
  CASE (iopt_FREE2B)
    IF(X_TRY(ISTT).LT.XMIN*DPARAM%MAXFREE_2B) X_TRY(ISTT) = XMIN*DPARAM%MAXFREE_2B
    IF(X_TRY(ISTT).GT.     DPARAM%MAXFREE_2B) X_TRY(ISTT) =      DPARAM%MAXFREE_2B
  CASE (iopt_WATR_2)
    ! *** SET LOWER LIMITS ***
    IF (SMODL%iARCH2.NE.iopt_topmdexp_2) THEN
     ! enforce lower limit
     IF (X_TRY(ISTT).LT.XMIN*MPARAM%MAXWATR_2) X_TRY(ISTT) = XMIN*MPARAM%MAXWATR_2
    ELSE
     ! MPARAM%MAXWATR_2 is just a scaling parameter, but don't allow stupid values
     IF (X_TRY(ISTT).LT.-MPARAM%MAXWATR_2*10._sp) X_TRY(ISTT) = -MPARAM%MAXWATR_2*10._sp
    ENDIF
    ! *** SET UPPER LIMITS ***
    IF (SMODL%iARCH2.EQ.iopt_tens2pll_2 .OR. SMODL%iARCH2.EQ.iopt_fixedsiz_2) THEN
     ! cannot exceed capacity
     IF (X_TRY(ISTT).GT.MPARAM%MAXWATR_2) X_TRY(ISTT) = MPARAM%MAXWATR_2
    ELSE
     ! unlimited storage, but make sure the values are still sensible
     !IF (X_TRY(ISTT).GT.MPARAM%MAXWATR_2*100._sp) X_TRY(ISTT) = MPARAM%MAXWATR_2*100._sp
    ENDIF
 END SELECT
END DO ! (loop through states)
! ---------------------------------------------------------------------------------------
END SUBROUTINE LIMIT_XTRY
END MODULE LIMIT_XTRY_MODULE
