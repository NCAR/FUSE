SUBROUTINE COMP_STATS()
! ----------------------------------------------------------------------------------------
! Creator:
!   Martyn Clark, 2009
!
! ----------------------------------------------------------------------------------------
! Purpose:
!   Used to compute summary statistics of model output
!
! ----------------------------------------------------------------------------------------
! Future revisions:
!
!   (add other summary statistics)
!
! ----------------------------------------------------------------------------------------
USE nrtype                                                 ! variable types (DP, I4B, etc.)
USE multistats
USE model_numerix
IMPLICIT NONE
! ----------------------------------------------------------------------------------------
! compute numerical stats
MSTATS%NUM_FUNCS     = MSTATS%NUM_FUNCS     + REAL(NUM_FUNCS, KIND(SP))     ! number of function calls
MSTATS%NUM_JACOBIAN  = MSTATS%NUM_JACOBIAN  + REAL(NUM_JACOBIAN, KIND(SP))  ! number of times Jacobian is calculated
MSTATS%NUMSUB_ACCEPT = MSTATS%NUMSUB_ACCEPT + REAL(NUMSUB_ACCEPT, KIND(SP)) ! number of sub-steps accepted (taken)
MSTATS%NUMSUB_REJECT = MSTATS%NUMSUB_REJECT + REAL(NUMSUB_REJECT, KIND(SP)) ! number of sub-steps tried but rejected
MSTATS%NUMSUB_NOCONV = MSTATS%NUMSUB_NOCONV + REAL(NUMSUB_NOCONV, KIND(SP)) ! number of sub-steps tried that did not converge
! compute maximum number of iterations
IF (MAXNUM_ITERNS > MSTATS%MAXNUM_ITERNS) MSTATS%MAXNUM_ITERNS = MAXNUM_ITERNS
! compute probability distributions
WHERE(ORD_NSUBS.GE.NUMSUB_ACCEPT) PRB_NSUBS = PRB_NSUBS + 1
! ----------------------------------------------------------------------------------------
END SUBROUTINE COMP_STATS
