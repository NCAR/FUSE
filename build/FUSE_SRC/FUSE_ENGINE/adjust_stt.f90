SUBROUTINE ADJUST_STT()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2008
! --------
! Modified by Dmitri Kavetski, 5 June 2013 AD (EAWAG) to replace IF with SELECTCASE
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Ensure that states are consistent with parameter values (needed for the special case of
!  stochastic parameters)
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! Model states updated in MODULE multistate
! ---------------------------------------------------------------------------------------
USE model_defn                                        ! model definitions
USE model_defnames
USE multistate                                        ! model states
USE multiparam                                        ! model parameters
IMPLICIT NONE
! internal
INTEGER(I4B)                           :: ISTT        ! loop through model states
! ----------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------
DO ISTT=1,NSTATE  ! NSTATE is in module model_defn
 SELECTCASE(CSTATE(ISTT)%iSNAME)
 ! ---------------------------------------------------------------------------------------
 ! states in the upper layer
 ! ---------------------------------------------------------------------------------------
 CASE (iopt_TENS1A)  ! tension 1a
  IF (MSTATE%TENS_1A .GT. DPARAM%MAXTENS_1A) MSTATE%TENS_1A=DPARAM%MAXTENS_1A
 CASE (iopt_TENS1B)  ! tension 1b
  IF (MSTATE%TENS_1B .GT. DPARAM%MAXTENS_1B) MSTATE%TENS_1B=DPARAM%MAXTENS_1B
 CASE (iopt_TENS_1)  ! tension 1
  IF (MSTATE%TENS_1  .GT. DPARAM%MAXTENS_1)  MSTATE%TENS_1 =DPARAM%MAXTENS_1
 CASE (iopt_FREE_1)  ! free 1
  IF (MSTATE%FREE_1  .GT. DPARAM%MAXFREE_1)  MSTATE%FREE_1 =DPARAM%MAXFREE_1
 CASE (iopt_WATR_1)  ! total 1
  IF (MSTATE%WATR_1  .GT. MPARAM%MAXWATR_1)  MSTATE%WATR_1 =MPARAM%MAXWATR_1
 ! ---------------------------------------------------------------------------------------
 ! states in the lower layer
 ! ---------------------------------------------------------------------------------------
 CASE (iopt_TENS_2)  ! tension 2
  IF (MSTATE%TENS_2  .GT. DPARAM%MAXTENS_2)  MSTATE%TENS_2 =DPARAM%MAXTENS_2
 CASE (iopt_FREE2A)  ! free 2a
  IF (MSTATE%FREE_2A .GT. DPARAM%MAXFREE_2A) MSTATE%FREE_2A=DPARAM%MAXFREE_2A
 CASE (iopt_FREE2B)  ! free 2b
  IF (MSTATE%FREE_2B .GT. DPARAM%MAXFREE_2B) MSTATE%FREE_2B=DPARAM%MAXFREE_2B
 CASE (iopt_WATR_2)  ! total 2
  IF (MSTATE%WATR_2  .GT. MPARAM%MAXWATR_2)  MSTATE%WATR_2 =MPARAM%MAXWATR_2
 END SELECT
END DO  ! (loop through model states)
! ----------------------------------------------------------------------------------------
END SUBROUTINE ADJUST_STT