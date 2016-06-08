SUBROUTINE INIT_STATS()
! ----------------------------------------------------------------------------------------
! Creator:
!   Martyn Clark, 2009
!
! ----------------------------------------------------------------------------------------
! Purpose:
!   Used to initialize summary statistics
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
! initialize numerical statistics
MSTATS%NUM_FUNCS     = 0 
MSTATS%NUM_JACOBIAN  = 0
MSTATS%NUMSUB_ACCEPT = 0
MSTATS%NUMSUB_REJECT = 0
MSTATS%NUMSUB_NOCONV = 0
! initialize probability distributions
PRB_NSUBS(:) = 0
! ----------------------------------------------------------------------------------------
END SUBROUTINE INIT_STATS
