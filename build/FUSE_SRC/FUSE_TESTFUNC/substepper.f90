SUBROUTINE SUBSTEPPER()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Run a given model and model parameter set for one time step, with adaptive sub-steps.
!
! The implicit solution is computed in the routine NEWTONITER, which finds the state vector "X_TRY"
! so that
!  X_TRY(:) = X_NEW(:)
!  X_NEW(:) = X_OLD(:) + DYDX(:) * HSTATE%STEP, with DYDX(:) evaluated at X_TRY(:)
!
! The "business=end" of the model is all within NEWTONITER (in the FUNCTION funcv) which computes
! model derivatives (DYDX) and model states (X_NEW) for a given state vector X_TRY(:)
!
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable definitions, etc.
USE newtoniter_mod, ONLY : newtoniter                 ! interface block for NEWTONITER
USE model_defn                                        ! model definitions
USE multiforce                                        ! model forcing data
USE multi_flux                                        ! model fluxes
USE multistate                                        ! model states
USE multiparam                                        ! model parameters
USE xtry_2_str_module                                 ! puts state vector into structure in multistate
USE str_2_xtry_module                                 ! gets state vector from structure in multistate
use model_numerix                                     ! define method/parameters used for numerical solution
IMPLICIT NONE
! internal variables
REAL(SP)                               :: STEP        ! new step size
REAL(SP)                               :: ETIME       ! part of the time step completed
REAL(SP)                               :: PREVSTEP    ! save pen-ultimate step size so small steps not carried over
REAL(SP), DIMENSION(:), ALLOCATABLE    :: X_START     ! state vector at start of time interval
REAL(SP), DIMENSION(:), ALLOCATABLE    :: X0_DYDT     ! derivative at X0 (start)
REAL(SP), DIMENSION(:), ALLOCATABLE    :: XM_DYDT     ! derivative at XM (middle)
REAL(SP), DIMENSION(:), ALLOCATABLE    :: X1_DYDT     ! derivative at X1 (end)
REAL(SP), DIMENSION(:), ALLOCATABLE    :: XC_DYDT     ! corrected derivatives
REAL(SP), DIMENSION(:), ALLOCATABLE    :: X_END0      ! explicit one-step solution, end of time interval
REAL(SP), DIMENSION(:), ALLOCATABLE    :: X_END1      ! implicit one-step solution, end of time interval
REAL(SP), DIMENSION(:), ALLOCATABLE    :: X_MID2      ! implicit two-step solution, middle of time interval
REAL(SP), DIMENSION(:), ALLOCATABLE    :: X_END2      ! implicit two-step solution, end of time interval
REAL(SP), DIMENSION(:), ALLOCATABLE    :: EVEC        ! error estimate for each state
REAL(SP), DIMENSION(:), ALLOCATABLE    :: TVEC        ! error threshold for each state
REAL(SP)                               :: DT          ! time step used in explicit euler
LOGICAL(LGT)                           :: ERROR_FLAG  ! .TRUE. if extrapolation error
LOGICAL(LGT)                           :: NEW_DERIVS  ! .TRUE. if need to calculate new derivatives
REAL(SP)                               :: STEPSAVE    ! save the time step (HSTATE%STEP altered for two-step solution)
REAL(SP)                               :: MULT        ! multiplier for new step size
REAL(SP), PARAMETER                    :: EPS=1.E-10_SP ! machine constant to prevent floating point errors
INTEGER(I4B)                           :: IERR        ! error code for allocate/deallocate
INTEGER(I4B), DIMENSION(1)             :: IMAX        ! index of maximum error
LOGICAL(LGT)                           :: CHECK       ! convergence check in SUBROUTINE newtoniter
INTEGER(I4B)                           :: NITER       ! number of iterations in newtoniter
REAL(SP)                               :: TEMPSTEP    ! suggested new time step, for case of non-convergence
REAL(SP)                               :: FTIM        ! fraction of model time interval to advance states
LOGICAL(LGT)                           :: NEWSTEP     ! FLAG to determine if a new Jacobian is needed
LOGICAL(LGT)                           :: STEP_INCREASE ! FLAG to determine if the end time step has been increased
INTEGER(I4B)                           :: I           ! looping variable
! interface blocks
INTERFACE
 SUBROUTINE limit_xtry(x)
 USE nrtype
 IMPLICIT NONE
 REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
 END SUBROUTINE limit_xtry
END INTERFACE
! ---------------------------------------------------------------------------------------
! (0) INITIALIZATION
! ---------------------------------------------------------------------------------------
ALLOCATE(X_START(NSTATE),X0_DYDT(NSTATE),XM_DYDT(NSTATE),X1_DYDT(NSTATE),XC_DYDT(NSTATE),&
         X_END0(NSTATE),X_END1(NSTATE),X_MID2(NSTATE),X_END2(NSTATE),&
         EVEC(NSTATE),TVEC(NSTATE), STAT=IERR)
IF (IERR.NE.0) STOP ' PROBLEM ALLOCATING SPACE IN MODEL1STEP '
ETIME         = 0._sp                        ! part of the time step completed
ASTATE        = FSTATE                       ! save model states at the start of the full step
newStep       = .true.                       ! initialize newstep (force re-calculation of Jacobian)
PREVSTEP      = HSTATE%STEP                  ! initialize the previous time step (used in next iteration)
STEP_INCREASE = .FALSE.                      ! used to check if the final sub-step has been increased
NUM_FUNCS     = 0                            ! number of function calls
NUM_JACOBIAN  = 0                            ! number of times Jacobian is calculated
NUMSUB_ACCEPT = 0                            ! number of sub-steps accepted (taken)
NUMSUB_REJECT = 0                            ! number of sub-steps tried but rejected
NUMSUB_NOCONV = 0                            ! number of sub-steps tried that did not converge
MAXNUM_ITERNS = 0                            ! maximum number of iterations taken in the newton method
! ---------------------------------------------------------------------------------------
! ensure time step is within bounds (can be out of bounds when processing remainder of last sub-step)
HSTATE%STEP = MIN( MAX(MIN_TSTEP,HSTATE%STEP), MAX_TSTEP) ! (MIN_TSTEP < stepsize < MAX_TSTEP)

SUBSTEPS: DO  ! continuous (recursive) loop over sub-steps

 ! ---------------------------------------------------------------------------------------
 ! (0) SAVE VECTOR OF STATES AND DERIVATIVES AT THE START OF THE SUB-STEP
 ! ---------------------------------------------------------------------------------------
 MSTATE = FSTATE                             ! model states at start of sub-step
 TSTATE = FSTATE; CALL STR_2_XTRY(X_START)   ! copy states (here TSTATE) to X_START
 ! determine if there is a need to calculate derivatives
 NEW_DERIVS=.FALSE.
 IF (ETIME.EQ.0._SP) THEN
  NEW_DERIVS=.TRUE.
 ELSE
  IF (SOLUTION_METHOD.EQ.EXPLICIT_EULER) NEW_DERIVS=.TRUE.
  IF (SOLUTION_METHOD.EQ.IMPLICIT_EULER .AND. & ! test for Crank-Nicholson
      TRUNCATION_ERROR.EQ.EMBEDDED_ERR  .AND. ORDER_ACCEPT.EQ.HIGHER_ORDER) NEW_DERIVS=.TRUE.
 ENDIF
 ! calculate new derivatives
 IF (NEW_DERIVS) THEN
  CALL MOD_DERIVS()                          ! model derivatives at start of sub-step
  FLUX_0 = M_FLUX                            ! save fluxes from explicit solution
  TSTATE=DY_DT; CALL STR_2_XTRY(X0_DYDT)     ! copy derivatives (here TSTATE) to X0_DYDT
 ELSE
  TSTATE=DYDT_OLD; CALL STR_2_XTRY(X0_DYDT)  ! copy derivatives (here TSTATE) to X0_DYDT
 ENDIF
 ! select solution method
 SELECT CASE(SOLUTION_METHOD)

  ! ---------------------------------------------------------------------------------------
  ! (1A) EXPLICIT ONE-STEP SOLUTION
  ! ---------------------------------------------------------------------------------------
  CASE (EXPLICIT_EULER)
   DT     = HSTATE%STEP                      ! define time step
   X_END1 = X_START + X0_DYDT*DT             ! explicit solution (can be out of range, but OK for error control)
   IF (ORDER_ACCEPT.EQ.LOWER_ORDER) THEN
    ! modify fluxes (M_FLUX) so that states (X_END1) are within bounds
    CALL XTRY_2_STR(X_START); BSTATE=TSTATE  ! populate state structure BSTATE with values of X_START
    CALL XTRY_2_STR(X_END1) ; ESTATE=TSTATE  ! populate state structure ESTATE with values of X_END1
    CALL FIX_STATES(DT,ERROR_FLAG)           ! ensure states are in bounds and disaggregate fluxes
    TSTATE=ESTATE; CALL STR_2_XTRY(X_END1)   ! copy states (here TSTATE) to X_END1
    ! EXIT here if there are no adaptive sub-steps
    IF (TEMPORAL_ERROR_CONTROL.EQ.TS_FIXED) THEN
     CALL WGT_FLUXES()                        ! just use W_FLUX=M_FLUX if no adaptive time steps
     EXIT SUBSTEPS                            ! EXIT the sub-steps loop
    ENDIF
    ! save M_FLUX, because modified below in MOD_DERIVS()
    FLUX_0 = M_FLUX  ! NOTE: unmodified FLUX_0 for higher-order solution is saved above
   ENDIF

   ! -------------------------------------------------------------------------------------
   ! (1B) EXPLICIT ERROR ESTIMATE
   ! -------------------------------------------------------------------------------------
   DT = HSTATE%STEP/2.0_SP                   ! define the time step
   X_MID2 = X_START + X0_DYDT*DT             ! explicit solution at the mid-point
   CALL XTRY_2_STR(X_START); BSTATE=TSTATE   ! populate state structure BSTATE with values of X_START
   ! ensure states are within range, and (if HIGHER_ORDER STEP_HALVING) make appropriate modifications
   IF (TRUNCATION_ERROR.EQ.STEP_HALVING .AND. ORDER_ACCEPT.EQ.HIGHER_ORDER) THEN
    CALL XTRY_2_STR(X_MID2) ; ESTATE=TSTATE  ! populate state structure ESTATE with values of X_MID2
    CALL FIX_STATES(DT,ERROR_FLAG)           ! ensure states are in bounds and disaggregate fluxes
    BSTATE=ESTATE                            ! set end state to start state
    TSTATE=ESTATE; CALL STR_2_XTRY(X_MID2)   ! copy states (here TSTATE) to X_MID2
    FLUX_1 = M_FLUX                          ! save fluxes
   ELSE
    CALL LIMIT_XTRY(X_MID2)                  ! ensure states are in bounds (no need to disagg fluxes)
    CALL XTRY_2_STR(X_MID2)                  ! populate state structure TSTATE with values of X_MID2
   ENDIF
   ! calculate derivative at the mid-point (TSTATE set above, TSTATE=ESTATE, or XTRY_2_STR(X_MID2))
   CALL MOD_DERIVS()                         ! evaluate dxdt for state vector TSTATE
   TSTATE = DY_DT; CALL STR_2_XTRY(XM_DYDT)  ! copy derivatives (here TSTATE) to XM_DYDT
   ! calculate different estimates of X_END2
   SELECT CASE(TRUNCATION_ERROR)
    CASE (STEP_HALVING)
     DT     = HSTATE%STEP/2.0_SP
     X_END2 = X_MID2  + XM_DYDT*DT           ! two-step method
    CASE (EMBEDDED_ERR)
     DT     = HSTATE%STEP
     X_END2 = X_START + XM_DYDT*DT           ! mid-point method
    CASE DEFAULT; STOP ' TRUNCATION_ERROR methods must be either STEP_HALVING or EMBEDDED_ERR '
   END SELECT ! select method for estimating temporal truncation error
   ! ensure states are within range, and make appropriate modifications (modifies M_FLUX)
   IF (ORDER_ACCEPT.EQ.HIGHER_ORDER) THEN
    CALL XTRY_2_STR(X_END2) ; ESTATE=TSTATE  ! populate state structure ESTATE with values of X_END2
    CALL FIX_STATES(DT,ERROR_FLAG)           ! ensure states are in bounds and disaggregate fluxes
    TSTATE=ESTATE; CALL STR_2_XTRY(X_END2)   ! copy states (here TSTATE) to X_END2
   ELSE
    M_FLUX = FLUX_0                          ! solution over the full step (saved earlier)
   ENDIF
   ! average fluxes for the two-step solution
   IF (TRUNCATION_ERROR.EQ.STEP_HALVING .AND. ORDER_ACCEPT.EQ.HIGHER_ORDER) THEN
    FLUX_2 = M_FLUX                          ! solution for the second half of the time step
    CALL MEANFLUXES()                        ! M_FLUX = FLUX_1 (first half) + FLUX_2 (second half)
   ENDIF

  ! --------------------------------------------------------------------------------------
  ! (2A) IMPLICIT ONE-STEP SOLUTION 
  ! --------------------------------------------------------------------------------------
  CASE (IMPLICIT_EULER)
   ! if use embedded error control, the "higher-order" solution is Crank-Nicholson, so
   ! need fluxes at the start of the current sub-step (calculated above)
   IF (TRUNCATION_ERROR.EQ.EMBEDDED_ERR) FLUX_1 = FLUX_0
   ! estimate the initial conditions used in the Newton scheme
   SELECT CASE (INITIAL_NEWTON) 
    CASE (STATE_OLD);     X_END1 = X_START
    CASE (EXPLICIT_MID);  X_END1 = X_START + X0_DYDT*HSTATE%STEP/2.0_SP ! estimate at mid-point
    CASE (EXPLICIT_FULL); X_END1 = X_START + X0_DYDT*HSTATE%STEP        ! estimate at end
   END SELECT
   ! estimate state vector at end of time step
   CALL NEWTONITER(X_END1,newStep,CHECK,NITER)   ! try different values of X until converge
   IF (NITER > MAXNUM_ITERNS) MAXNUM_ITERNS=NITER
   newStep=.false.
   ! just use this solution if no adaptive time steps
   IF (TEMPORAL_ERROR_CONTROL.EQ.TS_FIXED .AND. ORDER_ACCEPT.EQ.LOWER_ORDER) THEN
    CALL WGT_FLUXES()                        ! just use this solution if no adaptive time steps
    EXIT SUBSTEPS                            ! EXIT the sub-steps loop
   ENDIF
   ! save fluxes, if using lower-order solution
   IF (ORDER_ACCEPT.EQ.LOWER_ORDER) FLUX_0 = M_FLUX
   ! save fluxes at end of the current time step
   IF (TRUNCATION_ERROR.EQ.EMBEDDED_ERR) FLUX_2 = M_FLUX
   ! check for non-convergence
   IF (CHECK) THEN
    NUMSUB_NOCONV = NUMSUB_NOCONV + 1
    STEP = MAX(MIN_TSTEP, HSTATE%STEP*RMIN)  ! (avoid stepsize < MIN_TSTEP)
    TEMPSTEP = REVISE_STEP()  ! avoid small steps at the end of a time interval
    IF (TEMPSTEP.NE.STEP) THEN; PREVSTEP=STEP; ELSE; PREVSTEP=TEMPSTEP; ENDIF
    ! avoid the case of a continuous do loop where TEMPSTEP is at a minimum
    IF (TEMPSTEP.LT.HSTATE%STEP) THEN ! TEMPSTEP may equal HSTATE%STEP (MIN_TSTEP, or end of interval)
     newStep=.true.
     HSTATE%STEP=TEMPSTEP
     CYCLE SUBSTEPS
    ENDIF
    pause ' did not converge, and unable to make steps small enough '
   ENDIF
   ! -------------------------------------------------------------------------------------
   ! (2B) IMPLICIT ERROR ESTIMATE
   ! -------------------------------------------------------------------------------------
   SELECT CASE(TRUNCATION_ERROR)
    ! ------------------------------------------------------------------------------------
    ! temporal truncation error estimate = step halving
    CASE (STEP_HALVING)
     STEPSAVE=HSTATE%STEP                         ! need to alter HSTATE%STEP because used in FUNCV
     HSTATE%STEP = HSTATE%STEP/2._sp              ! new HSTATE%STEP for use in FUNCV
     ! implicit solution over the first half of the sub-step
     MSTATE = FSTATE                              ! model states at start of sub-step
     X_MID2 = X_START + X0_DYDT*HSTATE%STEP       ! explicit solution
     CALL NEWTONITER(X_MID2,newStep,CHECK,NITER)  ! solve for X_MID
     IF (NITER > MAXNUM_ITERNS) MAXNUM_ITERNS=NITER
     IF (NITER.GT.NITER_TOTAL) pause ' did not converge, two-step solution) '
     FLUX_1 = M_FLUX                              ! save fluxes over the first half of the time step
     ! implicit solution over the second half of the sub-step
     MSTATE  = TSTATE                             ! model states at start of next sub-step (TSTATE = X_MID2)
     TSTATE  = DY_DT                              ! temporarily populate TSTATE with derivatives  
     CALL STR_2_XTRY(XM_DYDT)                     ! copy derivatives (here TSTATE) to XM_DYDT
     X_END2  = X_MID2 + XM_DYDT*HSTATE%STEP       ! explicit solution
     CALL NEWTONITER(X_END2,newStep,CHECK,NITER)  ! try different values of X_END2 until converge
     IF (NITER > MAXNUM_ITERNS) MAXNUM_ITERNS=NITER
     IF (NITER.GT.NITER_TOTAL) pause ' did not converge, two-step solution '
     FLUX_2 = M_FLUX                              ! save fluxes over the second half of the time step
     ! calculate fluxes used in WGT_FLUXES()
     IF (ORDER_ACCEPT.EQ.HIGHER_ORDER) THEN
      CALL MEANFLUXES()                           ! M_FLUX = average explicit (FLUX_1) + implicit (FLUX_2) solution
     ELSE
      M_FLUX    = FLUX_0                          ! just use implicit one-step solution
      DYDT_OLD  = DY_DT                           ! save derivatives
     ENDIF
     HSTATE%STEP = STEPSAVE                       ! re-set time step again (used in FUNCV)
    ! -------------------------------------------------------------------------------------
    ! temporal truncation error estimate = embedded error estimate
    CASE (EMBEDDED_ERR)
     ! get derivative vector at the end of the time step (NOTE: don't enter two-step case)
     TSTATE  = DY_DT                              ! temporarily populate TSTATE with derivatives
     CALL STR_2_XTRY(X1_DYDT)                     ! copy derivatives (here TSTATE) to X1_DYDT
     ! alternative solution
     DT     = HSTATE%STEP
     X_END2 = X_START + 0.5_SP*(X0_DYDT+X1_DYDT)*DT
     ! ensure states are within range, and make appropriate modifications
     IF (ORDER_ACCEPT.EQ.HIGHER_ORDER) THEN
      CALL MEANFLUXES()                           ! M_FLUX = average explicit (FLUX_1) + implicit (FLUX_2) solution
      CALL XTRY_2_STR(X_START); BSTATE=TSTATE     ! populate state structure BSTATE with values of X_START
      CALL XTRY_2_STR(X_END2) ; ESTATE=TSTATE     ! populate state structure ESTATE with values of X_END2
      CALL FIX_STATES(DT,ERROR_FLAG)              ! ensure states are in bounds and disaggregate fluxes (M_FLUX)
      TSTATE=ESTATE; CALL STR_2_XTRY(X_END2)      ! copy states (here TSTATE) to X_END2
     ENDIF
    CASE DEFAULT
     STOP ' TRUNCATION_ERROR methods must be either STEP_HALVING or EMBEDDED_ERR '
   END SELECT ! select method for estimating temporal truncation error
  CASE DEFAULT; STOP ' SOLUTION_METHOD must be either EXPLICIT_EULER or IMPLICIT_EULER '
 END SELECT  ! select method for numerical solution  

 ! --------------------------------------------------------------------------------------
 ! (4) CALCULATE ERROR, CHECK IF ACCEPT/REJECT THE CURRENT STEP, AND NEW STEP SIZE
 ! --------------------------------------------------------------------------------------
 ! calculate the maximum error over all states
 EVEC = ABS(X_END2-X_END1)                        ! error estimate
 TVEC = ERR_TRUNC_REL*ABS(X_END2) + ERR_TRUNC_ABS ! error thresholds
 IMAX = MAXLOC(EVEC - TVEC)                       ! index of maximum error
 !WRITE(*,'(A10,1X,10(E12.5,1X))') 'X_START', ETIME, HSTATE%STEP, X_START
 !WRITE(*,'(A10,1X,10(E12.5,1X))')  'X_END0', ETIME, HSTATE%STEP, X_END0
 !WRITE(*,'(A10,1X,10(E12.5,1X))')  'X_MID2', ETIME, HSTATE%STEP, X_MID2
 !WRITE(*,'(A10,1X,10(E12.5,1X))')  'X_END1', ETIME, HSTATE%STEP, X_END1
 !WRITE(*,'(A10,1X,10(E12.5,1X))')  'X_END2', ETIME, HSTATE%STEP, X_END2
 !WRITE(*,'(A10,1X,10(E12.5,1X))')    'EVEC', ETIME, HSTATE%STEP, EVEC
 !WRITE(*,'(A10,1X,10(E12.5,1X))')    'TVEC', ETIME, HSTATE%STEP, TVEC
 ! --------------------------------------------------------------------------------------
 ! check to accept time step
 IF (TEMPORAL_ERROR_CONTROL.EQ.TS_FIXED .OR. & ! (accept if using fixed time steps)
      EVEC(IMAX(1)) < TVEC(IMAX(1)) .OR. &     ! (accept if error is less than critical threshold)
      HSTATE%STEP <= MIN_TSTEP) THEN           ! (accept if time step is already minimum allowable)
  ! accept step -- calculate new (increased) step size
  ! NOTE: step size not necessarily increased because of the safety factor
  IF (TEMPORAL_ERROR_CONTROL.EQ.TS_ADAPT) THEN
   MULT = SAFETY * SQRT( TVEC(IMAX(1)) / MAX(EVEC(IMAX(1)),EPS) )
   STEP = MIN( MAX(MIN_TSTEP, HSTATE%STEP * MIN(MULT,RMAX) ), MAX_TSTEP) ! (MIN_TSTEP < stepsize < MAX_TSTEP)
  ELSE
   STEP = MAX_TSTEP
  ENDIF
  ! add contribution of sub-step flux to the timestep-average flux
  !print *, 'm_flux%qbase_2a = ', m_flux%qbase_2a
  CALL WGT_FLUXES()                       ! add M_FLUX to W_FLUX
  ! save states at the end of the sub-step
  SELECT CASE (ORDER_ACCEPT)
   CASE  (LOWER_ORDER); CALL XTRY_2_STR(X_END1) ! populate TSTATE with X_END1
   CASE (HIGHER_ORDER); CALL XTRY_2_STR(X_END2) ! populate TSTATE with X_END2
  END SELECT
  FSTATE    = TSTATE                            ! states at the end of the sub-step
  ! save derivatives at the end of the sub-step
  ! NOTE: explicit euler solution calculated at start of SUBSTEP loop (no need to save derivatives)
  IF (SOLUTION_METHOD.EQ.IMPLICIT_EULER) THEN
   ! NOTE: derivatives for implicit one-step solution saved earlier 
   IF (TRUNCATION_ERROR.EQ.STEP_HALVING .AND. ORDER_ACCEPT.EQ.HIGHER_ORDER) DYDT_OLD  = DY_DT
   ! NOTE: implicit Crank-Nicholson solution also calculated at start of SUBSTEP loop
   IF (TRUNCATION_ERROR.EQ.EMBEDDED_ERR .AND. ORDER_ACCEPT.EQ. LOWER_ORDER) DYDT_OLD  = DY_DT
  ENDIF
  ! keep track of the number of sub-steps taken
  NUMSUB_ACCEPT = NUMSUB_ACCEPT + 1
  !print *, 'accept step ', numsub_accept
  ! compute fraction of big step that is finished, and check for exit criteria
  ETIME = ETIME + HSTATE%STEP             ! identify position within the time step
  IF (ETIME.GE.DELTIM) EXIT
  ! revise the length of time steps to avoid small steps at the end of a time interval
  HSTATE%STEP = REVISE_STEP()  ! avoid small steps at the end of a time interval
  IF (HSTATE%STEP.NE.STEP) THEN; PREVSTEP=STEP; ELSE; PREVSTEP=HSTATE%STEP; ENDIF
 ! --------------------------------------------------------------------------------------
 ELSE  ! REJECT STEP AND DECREASE STEP SIZE
  ! calculate new (decreased) step size
  !print *, 'reject step ', numsub_reject
  NUMSUB_REJECT = NUMSUB_REJECT + 1
  MULT = SAFETY * SQRT( TVEC(IMAX(1)) / MAX(EVEC(IMAX(1)),EPS) )
  STEP = MAX(MIN_TSTEP, HSTATE%STEP * MAX(MULT,RMIN) ) ! (avoid stepsize < MIN_TSTEP)
  HSTATE%STEP = REVISE_STEP()  ! avoid small steps at the end of a time interval
  IF (HSTATE%STEP.NE.STEP) THEN; PREVSTEP=STEP; ELSE; PREVSTEP=HSTATE%STEP; ENDIF
 ENDIF
 !print *, prevstep, step
 !IF (NUMSUB_REJECT.GT.10000) PAUSE
 ! (keep looping)

END DO SUBSTEPS ! continuous (recursive) do loop

! ---------------------------------------------------------------------------------------
! (9) RE-COMPUTE STATES AT THE END OF THE FULL STEP
! ---------------------------------------------------------------------------------------
! The implicit solution is not exact.  To conserve mass, we uses the weighted average of
! model fluxes throughout the time step to re-compute states at the end of the time step
! ---------------------------------------------------------------------------------------
!WRITE(*,'(A15,1X,5(E15.8,1X))') 'FINAL FLUXES = ', &
! W_FLUX%QSURF, W_FLUX%OFLOW_1, W_FLUX%QINTF_1, W_FLUX%OFLOW_2, W_FLUX%QBASE_2
! update model states
FTIM   = DELTIM                            ! fraction of time step in subroutine updatstate  
M_FLUX = W_FLUX                            ! SUBROUTINE mstate_eqn uses M_FLUX
FSTATE = ASTATE                            ! state at the start of the time step
CALL MSTATE_EQN()                          ! use time-step-average fluxes to compute model derivatives 
CALL UPDATSTATE(FTIM)                      ! update model states
W_FLUX%SATAREA = W_FLUX%SATAREA/DELTIM     ! normalize saturated area (weighted sum over sub-steps)
HSTATE%STEP=PREVSTEP                       ! ensure stepsize is not equal to the small remainder
! ---------------------------------------------------------------------------------------
DEALLOCATE(X_START,X0_DYDT,XM_DYDT,X1_DYDT,XC_DYDT,X_END0,X_END1,X_MID2,X_END2,EVEC,TVEC, &
           STAT=IERR); IF (IERR.NE.0) STOP ' PROBLEM DEALLOCATING SPACE IN MODEL1STEP '
! ---------------------------------------------------------------------------------------
CONTAINS
 FUNCTION REVISE_STEP()
 REAL(SP)    :: REVISE_STEP
 REAL(SP)    :: T_MGN
 SELECT CASE(SMALL_ENDSTEP)
  ! -------------------------------------------------------------------------------------
  CASE(STEP_TRUNC)  ! truncate the time step if near the end
   IF (ETIME + STEP .GE. DELTIM) REVISE_STEP = DELTIM - ETIME
   IF (ETIME + STEP .LT. DELTIM) REVISE_STEP = STEP 
  ! -------------------------------------------------------------------------------------
  CASE(LOOK_AHEAD)  ! the look-ahead method of Shampine (1994)
   IF (ETIME + STEP .GE. DELTIM) THEN
    REVISE_STEP = DELTIM - ETIME
   ELSE
    IF (ETIME + STEP*2._SP .GE. DELTIM) THEN
     REVISE_STEP = (DELTIM - ETIME)/2._SP
    ELSE
     REVISE_STEP = STEP
    ENDIF
   ENDIF
  ! -------------------------------------------------------------------------------------
  CASE(STEP_ABSORB) ! the step-absorption method
   IF (STEP_INCREASE) THEN  ! only try and increase step size once
    REVISE_STEP = STEP
   ELSE
    T_MGN = STEP/SAFETY - STEP  ! margin of error
    IF (ETIME + STEP + T_MGN .GE. DELTIM) THEN
     REVISE_STEP   = DELTIM - ETIME
     STEP_INCREASE = .TRUE. 
    ELSE
     IF (ETIME + STEP + T_MGN*2._SP .GE. DELTIM) THEN
      REVISE_STEP   = STEP + T_MGN*(1._SP - (DELTIM-(ETIME+STEP))/T_MGN)
      STEP_INCREASE = .TRUE.
     ELSE
      REVISE_STEP = STEP
     ENDIF
    ENDIF
   ENDIF
  CASE DEFAULT; STOP ' must use the STEP_TRUNC, LOOK_AHEAD, or STEP_ABSORB methods '
 END SELECT
 END FUNCTION REVISE_STEP
END SUBROUTINE SUBSTEPPER
