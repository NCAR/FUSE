MODULE multistats
 USE nrtype
 TYPE SUMMARY
  ! DMSL diagnostix
  REAL(SP)                             :: VAR_RESIDUL   ! variance of the model residuals
  REAL(SP)                             :: LOGP_SIMULN   ! log density of the model simulation
  REAL(SP)                             :: JUMP_TAKEN    ! defines a jump in the MCMC production run
  ! comparisons between model output and observations
  REAL(SP)                             :: QOBS_MEAN     ! mean observed runoff (mm day-1)
  REAL(SP)                             :: QSIM_MEAN     ! mean simulated runoff (mm day-1)
  REAL(SP)                             :: QOBS_CVAR     ! coefficient of variation of observed runoff (-)
  REAL(SP)                             :: QSIM_CVAR     ! coefficient of variation of simulated runoff (-)
  REAL(SP)                             :: QOBS_LAG1     ! lag-1 correlation of observed runoff (-)
  REAL(SP)                             :: QSIM_LAG1     ! lag-1 correlation of simulated runoff (-)
  REAL(SP)                             :: RAW_RMSE      ! root-mean-squared-error of flow (mm day-1)
  REAL(SP)                             :: LOG_RMSE      ! root-mean-squared-error of LOG flow (mm day-1)
  REAL(SP)                             :: NASH_SUTT     ! Nash-Sutcliffe score
  ! attributes of model output
  REAL(SP)                             :: NUM_RMSE      ! error of the approximate solution
  REAL(SP)                             :: NUM_FUNCS     ! number of function calls
  REAL(SP)                             :: NUM_JACOBIAN  ! number of times Jacobian is calculated
  REAL(SP)                             :: NUMSUB_ACCEPT ! number of sub-steps taken
  REAL(SP)                             :: NUMSUB_REJECT ! number of sub-steps taken
  REAL(SP)                             :: NUMSUB_NOCONV ! number of sub-steps tried that did not converge
  INTEGER(I4B)                         :: MAXNUM_ITERNS ! maximum number of iterations in implicit scheme
  REAL(SP), DIMENSION(20)              :: NUMSUB_PROB   ! probability distribution for number of sub-steps
  ! error checking
  CHARACTER(LEN=1024)                  :: ERR_MESSAGE   ! error message
 ENDTYPE SUMMARY
 ! final data structures
 TYPE(SUMMARY)                         :: MSTATS        ! (model summary statistics)
 INTEGER(I4B)                          :: MOD_IX=1      ! (model index)
 INTEGER(I4B)                          :: PCOUNT        ! (number of parameter sets in model output files)
 INTEGER(I4B)                          :: FCOUNT        ! (number of model simulations)
END MODULE multistats
