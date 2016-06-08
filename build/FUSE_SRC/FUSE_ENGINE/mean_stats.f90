SUBROUTINE MEAN_STATS()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Computes summary statistics from model simulations
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multistats -- summary statistics stored in MODULE multistats
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
! FUSE modules
USE multiforce                                        ! model forcing data (obs streamflow)
USE multiroute                                        ! routed runoff
USE multistats                                        ! summary statistics
USE model_numerix                                     ! model numerix parameters and data
IMPLICIT NONE
! internal
INTEGER(I4B)                           :: I           ! looping
INTEGER(I4B)                           :: NS          ! number of samples
INTEGER(I4B)                           :: IERR        ! error code for allocate/deallocate statements
REAL(SP), DIMENSION(:), ALLOCATABLE    :: QOBS        ! observed runoff
REAL(SP), DIMENSION(:), ALLOCATABLE    :: QSIM        ! simulated runoff
REAL(SP), DIMENSION(:), ALLOCATABLE    :: DOBS        ! observed runoff anomalies
REAL(SP), DIMENSION(:), ALLOCATABLE    :: DSIM        ! simulated runoff anomalies
REAL(SP), DIMENSION(:), ALLOCATABLE    :: RAWD        ! observed-simulated differences in flow
REAL(SP), DIMENSION(:), ALLOCATABLE    :: LOGD        ! observed-simulated differences in LOG flow
REAL(SP)                               :: XB_OBS      ! mean observed runoff
REAL(SP)                               :: XB_SIM      ! mean simulated runoff
REAL(SP)                               :: SS_OBS      ! sum of squared observed runoff anomalies
REAL(SP)                               :: SS_SIM      ! sum of squared simulated runoff anomalies
REAL(SP)                               :: SS_LOBS     ! sum of squared lagged differences in observed runoff
REAL(SP)                               :: SS_LSIM     ! sum of squared lagged differences in simulated runoff
REAL(SP)                               :: SS_RAW      ! sum of squared differences in observed - simulated
REAL(SP)                               :: SS_LOG      ! sum of squared differences in LOG observed - LOG simulated
REAL(SP), PARAMETER                    :: NO_ZERO=1.E-20  ! avoid divide by zero
! ---------------------------------------------------------------------------------------
! (1) PRELIMINARIES
! ---------------------------------------------------------------------------------------
! define sample size
NS = (NUMTIM-ISTART) + 1    ! (ISTART is shared in MODULE multiforce)
! allocate space for observed and simulated runoff
ALLOCATE(QOBS(NS),QSIM(NS),DOBS(NS),DSIM(NS),RAWD(NS),LOGD(NS),STAT=IERR)
IF (IERR.NE.0) STOP ' PROBLEM ALLOCATING SPACE IN MEAN_STATS.F90 '
! extract vectors from data structures
QOBS = AFORCE(ISTART:NUMTIM)%OBSQ
QSIM = AROUTE(ISTART:NUMTIM)%Q_ROUTED
! compute mean
XB_OBS = SUM(QOBS(:)) / REAL(NS, KIND(SP))
XB_SIM = SUM(QSIM(:)) / REAL(NS, KIND(SP))
! compute the sum of squares of simulated and observed vectors
DOBS(:) = QOBS(:) - XB_OBS
DSIM(:) = QSIM(:) - XB_SIM
SS_OBS  = DOT_PRODUCT(DOBS,DOBS)  ! = SUM( DOBS(:)*DOBS(:) )
SS_SIM  = DOT_PRODUCT(DSIM,DSIM)  ! = SUM( DSIM(:)*DSIM(:) )
! compute the sum of squares of lagged differences
SS_LOBS = DOT_PRODUCT(DOBS(2:NS),DOBS(1:NS-1))
SS_LSIM = DOT_PRODUCT(DSIM(2:NS),DSIM(1:NS-1))
! compute sum of squared differences between model and observations
RAWD(:) = QSIM(:) - QOBS(:)
LOGD(:) = LOG(QSIM(:)+TINY(QSIM)) - LOG(QOBS(:)+TINY(QOBS)) ! TINY ADDED TO AVOID TROUBLES WHEN Q=0
SS_RAW  = DOT_PRODUCT(RAWD,RAWD)  ! = SUM( RAWD(:)*RAWD(:) )
SS_LOG  = DOT_PRODUCT(LOGD,LOGD)  ! = SUM( LOGD(:)*LOGD(:) )
! ---------------------------------------------------------------------------------------
! (2) COMPUTE ERROR STATISTICS
! ---------------------------------------------------------------------------------------
! compute the mean
MSTATS%QOBS_MEAN = XB_OBS
MSTATS%QSIM_MEAN = XB_SIM
! compute the coefficient of variation 
MSTATS%QOBS_CVAR = SQRT( SS_OBS / REAL(NS-1, KIND(SP)) ) / (XB_OBS+NO_ZERO)
MSTATS%QSIM_CVAR = SQRT( SS_SIM / REAL(NS-1, KIND(SP)) ) / (XB_SIM+NO_ZERO)
! compute the lag-1 correlation coefficient
MSTATS%QOBS_LAG1 = SS_LOBS / (SQRT(SS_OBS*SS_OBS)+NO_ZERO)
MSTATS%QSIM_LAG1 = SS_LSIM / (SQRT(SS_SIM*SS_SIM)+NO_ZERO)
! compute the root-mean-squared-error of flow
MSTATS%RAW_RMSE  = SQRT( SS_RAW / REAL(NS, KIND(SP)) )
! compute the root-mean-squared-error of LOG flow
MSTATS%LOG_RMSE  = SQRT( SS_LOG / REAL(NS, KIND(SP)) )
! compute the Nash-Sutcliffe score
MSTATS%NASH_SUTT = 1. - SS_RAW/(SS_OBS+NO_ZERO)
! ---------------------------------------------------------------------------------------
! (4) COMPUTE STATISTICS ON NUMERICAL ACCURACY AND EFFICIENCY
! ---------------------------------------------------------------------------------------
! compute RMSE between "more accurate" and "less accurate" solutions
QOBS = AROUTE(ISTART:NUMTIM)%Q_ACCURATE
RAWD(:) = QSIM(:) - QOBS(:); SS_RAW  = DOT_PRODUCT(RAWD,RAWD)    ! = SUM( RAWD(:)*RAWD(:) )
MSTATS%NUM_RMSE = SQRT( SS_RAW / REAL(NS, KIND(SP)) )
! compute summary statistics for efficiency
MSTATS%NUM_FUNCS     = MSTATS%NUM_FUNCS     / REAL(NUMTIM, KIND(SP)) ! number of function calls
MSTATS%NUM_JACOBIAN  = MSTATS%NUM_JACOBIAN  / REAL(NUMTIM, KIND(SP)) ! number of times Jacobian is calculated
MSTATS%NUMSUB_ACCEPT = MSTATS%NUMSUB_ACCEPT / REAL(NUMTIM, KIND(SP)) ! number of sub-steps accepted (taken)
MSTATS%NUMSUB_REJECT = MSTATS%NUMSUB_REJECT / REAL(NUMTIM, KIND(SP)) ! number of sub-steps tried but rejected
MSTATS%NUMSUB_NOCONV = MSTATS%NUMSUB_NOCONV / REAL(NUMTIM, KIND(SP)) ! number of sub-steps tried that did not converge
! compute cumulative probability distributions
MSTATS%NUMSUB_PROB   = REAL(PRB_NSUBS(:), KIND(SP)) / REAL(NUMTIM, KIND(SP))
! ---------------------------------------------------------------------------------------
DEALLOCATE(QOBS,QSIM,DOBS,DSIM,RAWD,LOGD,STAT=IERR)
IF (IERR.NE.0) STOP ' PROBLEM DEALLOCATING SPACE IN MEAN_STATS.F90 '
! ---------------------------------------------------------------------------------------
END SUBROUTINE MEAN_STATS
