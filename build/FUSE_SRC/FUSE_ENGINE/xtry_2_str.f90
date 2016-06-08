MODULE XTRY_2_STR_MODULE
IMPLICIT NONE
CONTAINS
SUBROUTINE XTRY_2_STR(X_TRY,TMPSTR)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Populate the temporary state structure with values of X_TRY
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! Temporary model states updated in MODULE multistate
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! Numerical Recipes data types
USE model_defn, ONLY: CSTATE,NSTATE,SMODL             ! model definitions
USE model_defnames
USE multistate, ONLY: STATEV                          ! model states
USE multiparam, ONLY: DPARAM                          ! model parameters
IMPLICIT NONE
! input
REAL(SP), DIMENSION(:), INTENT(IN)     :: X_TRY       ! vector of model states
! output
TYPE(STATEV), INTENT(OUT)              :: TMPSTR      ! temporary state structure
! internal
INTEGER(I4B)                           :: ISTT        ! loop through model states
REAL(SP),PARAMETER::missingValue=-9999._sp
! ---------------------------------------------------------------------------------------
! (A) POPULATE THE TEMPORARY STATE STRUCTURE WITH VALUES OF XSTATE
! ---------------------------------------------------------------------------------------
DO ISTT=1,NSTATE
 SELECT CASE(CSTATE(ISTT)%iSNAME)
  CASE (iopt_TENS1A); TMPSTR%TENS_1A = X_TRY(ISTT)
  CASE (iopt_TENS1B); TMPSTR%TENS_1B = X_TRY(ISTT)
  CASE (iopt_TENS_1); TMPSTR%TENS_1  = X_TRY(ISTT)
  CASE (iopt_FREE_1); TMPSTR%FREE_1  = X_TRY(ISTT)
  CASE (iopt_WATR_1); TMPSTR%WATR_1  = X_TRY(ISTT)
  CASE (iopt_TENS_2); TMPSTR%TENS_2  = X_TRY(ISTT)
  CASE (iopt_FREE2A); TMPSTR%FREE_2A = X_TRY(ISTT)
  CASE (iopt_FREE2B); TMPSTR%FREE_2B = X_TRY(ISTT)
  CASE (iopt_WATR_2); TMPSTR%WATR_2  = X_TRY(ISTT)
 END SELECT
END DO  ! istt
! ---------------------------------------------------------------------------------------
! (B) ESTIMATE THE "STATES" THAT ARE NOT REALLY STATES
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iARCH1)  ! (upper layer architecture)
 CASE(iopt_onestate_1) ! upper layer defined by a single state variable
  TMPSTR%TENS_1A = missingValue                              ! 1st tension store (undefined)
  TMPSTR%TENS_1B = missingValue                              ! 2nd tension store (undefined)
  TMPSTR%TENS_1  = MIN(TMPSTR%WATR_1, DPARAM%MAXTENS_1)      ! tension storage
  TMPSTR%FREE_1  = MAX(0._sp, TMPSTR%WATR_1 - DPARAM%MAXTENS_1) ! free storage
 CASE(iopt_tension1_1) ! upper layer broken up into tension and free storage
  TMPSTR%TENS_1A = missingValue                              ! 1st tension store (undefined)
  TMPSTR%TENS_1B = missingValue                              ! 2nd tension store (undefined)
  TMPSTR%WATR_1  = TMPSTR%TENS_1 + TMPSTR%FREE_1             ! total storage
 CASE(iopt_tension2_1) ! tension storage sub-divided into recharge and excess
  TMPSTR%TENS_1  = TMPSTR%TENS_1A + TMPSTR%TENS_1B           ! tension storage
  TMPSTR%WATR_1  = TMPSTR%TENS_1  + TMPSTR%FREE_1            ! total storage
 CASE DEFAULT       ! (error check)
  print *, "MDEFN(IMOD)%iARCH1 must be iopt_tension2_1, iopt_tension1_1, or iopt_onestate_1"
  STOP
END SELECT  ! (upper layer architechure)
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iARCH2)  ! (lower layer architecture)
 CASE(iopt_tens2pll_2) ! tension reservoir plus two parallel tanks
  TMPSTR%FREE_2  = TMPSTR%FREE_2A + TMPSTR%FREE_2B           ! free storage
  TMPSTR%WATR_2  = TMPSTR%TENS_2  + TMPSTR%FREE_2            ! total storage
 CASE(iopt_unlimfrc_2,iopt_unlimpow_2,iopt_topmdexp_2,iopt_fixedsiz_2) ! single baseflow reservoir
  TMPSTR%TENS_2  = MIN(TMPSTR%WATR_2, DPARAM%MAXTENS_2)      ! tension storage
  TMPSTR%FREE_2  = MAX(0._sp, TMPSTR%WATR_2 - DPARAM%MAXTENS_2) ! free storage
  TMPSTR%FREE_2A = missingValue                              ! primary reservoir (undefined) 
  TMPSTR%FREE_2A = missingValue                              ! secondary reservoir (undefined) 
 CASE DEFAULT       ! (error check)
  print *, "MDEFN(IMOD)%ARCH2 must be 'tens2pll_2', 'unlimfrc_2', 'unlimpow_2'"
  print *, "  'topmdexp_2', or 'fixedsiz_2'"
  STOP
END SELECT
! ---------------------------------------------------------------------------------------
END SUBROUTINE XTRY_2_STR
END MODULE XTRY_2_STR_MODULE
