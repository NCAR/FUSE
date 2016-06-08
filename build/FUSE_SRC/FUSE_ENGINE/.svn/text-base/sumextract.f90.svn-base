MODULE SUMEXTRACT_MODULE
IMPLICIT NONE
CONTAINS
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
PURE FUNCTION SUMEXTRACT(STATNAME)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Extracts variable "VNAME(IVAR)" from relevant data structures
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE multistats                                        ! summary statistics
IMPLICIT NONE
! input
CHARACTER(*), INTENT(IN)               :: STATNAME    ! variable name
! internal
REAL(SP)                               :: XVAR        ! variable
! output
REAL(SP)                               :: SUMEXTRACT  ! FUNCTION name
! ---------------------------------------------------------------------------------------
! initialize XVAR
XVAR=-9999._sp
! DMSL diagnostix
IF (TRIM(STATNAME).EQ.'var_residul') XVAR = MSTATS%VAR_RESIDUL
IF (TRIM(STATNAME).EQ.'logp_simuln') XVAR = MSTATS%LOGP_SIMULN
IF (TRIM(STATNAME).EQ.'jump_taken')  XVAR = MSTATS%JUMP_TAKEN
! extract summary statistics
IF (TRIM(STATNAME).EQ.'qobs_mean')   XVAR = MSTATS%QOBS_MEAN
IF (TRIM(STATNAME).EQ.'qsim_mean')   XVAR = MSTATS%QSIM_MEAN 
IF (TRIM(STATNAME).EQ.'qobs_cvar')   XVAR = MSTATS%QOBS_CVAR
IF (TRIM(STATNAME).EQ.'qsim_cvar')   XVAR = MSTATS%QSIM_CVAR
IF (TRIM(STATNAME).EQ.'qobs_lag1')   XVAR = MSTATS%QOBS_LAG1
IF (TRIM(STATNAME).EQ.'qsim_lag1')   XVAR = MSTATS%QSIM_LAG1
IF (TRIM(STATNAME).EQ.'raw_rmse')    XVAR = MSTATS%RAW_RMSE
IF (TRIM(STATNAME).EQ.'log_rmse')    XVAR = MSTATS%LOG_RMSE
IF (TRIM(STATNAME).EQ.'nash_sutt')   XVAR = MSTATS%NASH_SUTT
! extract numerix stats
IF (TRIM(STATNAME).EQ.'numerx_rmse') XVAR = MSTATS%NUM_RMSE      
IF (TRIM(STATNAME).EQ.'mean_nfuncs') XVAR = MSTATS%NUM_FUNCS
IF (TRIM(STATNAME).EQ.'mean_njacob') XVAR = MSTATS%NUM_JACOBIAN
IF (TRIM(STATNAME).EQ.'mean_accept') XVAR = MSTATS%NUMSUB_ACCEPT
IF (TRIM(STATNAME).EQ.'mean_reject') XVAR = MSTATS%NUMSUB_REJECT
IF (TRIM(STATNAME).EQ.'mean_noconv') XVAR = MSTATS%NUMSUB_NOCONV
IF (TRIM(STATNAME).EQ.'maxnum_iter') XVAR = REAL(MSTATS%MAXNUM_ITERNS, KIND(SP))
! and, save the output
SUMEXTRACT = XVAR
! ---------------------------------------------------------------------------------------
END FUNCTION SUMEXTRACT
END MODULE SUMEXTRACT_MODULE
