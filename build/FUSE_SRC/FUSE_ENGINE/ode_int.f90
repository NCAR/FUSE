SUBROUTINE ODE_INT(MODL_SOLVE,STATE_START,STATE_END,DT_SUB,DT_FULL,IERR,MESSAGE)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
!
! Used for the temporal integration of ordinary differential equations, using different
!  numerical methods
!
! Based on the FUSE "sub-stepper" routine, but all FUSE-specific data structures have
!  been stripped out to call a simple test function
!
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable definitions, etc.
USE model_numerix                                     ! define method/parameters used for numerical solution
IMPLICIT NONE
! input/output variables
REAL(SP), DIMENSION(:), INTENT(IN)     :: STATE_START ! state vector at the start of the full step
REAL(SP), DIMENSION(:), INTENT(OUT)    :: STATE_END   ! state vector at the end of the full step
REAL(SP), INTENT(INOUT)                :: DT_SUB      ! length of the sub-step
REAL(SP), INTENT(IN)                   :: DT_FULL     ! length of the full step
INTEGER(I4B), INTENT(OUT)              :: IERR        ! error code
CHARACTER(LEN=*), INTENT(OUT)          :: MESSAGE     ! error message
! internal variables
REAL(SP)                               :: STEP        ! new step size
REAL(SP)                               :: ETIME       ! part of the time step completed
REAL(SP)                               :: PREVSTEP    ! save pen-ultimate step size so small steps not carried over
LOGICAL(LGT)                           :: NEWSTEP     ! .TRUE. if new step (determine if a new Jacobian is needed)
LOGICAL(LGT)                           :: NEW_SUBSTEP ! .TRUE. if new sub-step (determine if need to calculate derivatives)
LOGICAL(LGT)                           :: STEP_INCREASE ! FLAG to determine if the end time step has been increased
REAL(SP), DIMENSION(SIZE(STATE_START)) :: STATE0      ! state vector at the start of the sub-step
REAL(SP), DIMENSION(SIZE(STATE_START)) :: STATE1_LO   ! state vector at the end of the sub-step (lower-order solution)
REAL(SP), DIMENSION(SIZE(STATE_START)) :: STATE1_HI   ! state vector at the end of the sub-step (higher-order solution)
REAL(SP), DIMENSION(SIZE(STATE_START)) :: STATE1_LO_S ! safeguarded explicit Euler solution, also used in explicit Heun
REAL(SP), DIMENSION(SIZE(STATE_START)) :: STATE1_HI_S ! safeguarded explicit Heun and implicit Heun solutions
REAL(SP), DIMENSION(SIZE(STATE_START)) :: STATE1_INIT ! initial state vector used in the implicit solution
REAL(SP), DIMENSION(SIZE(STATE_START)) :: STATE1_RETAIN ! states retained at the end of the sub-step
REAL(SP), DIMENSION(SIZE(STATE_START)) :: DYDT_0      ! model derivatives at the start of the sub-step
REAL(SP), DIMENSION(SIZE(STATE_START)) :: DYDT_1      ! model derivatives at the end of the sub-step
REAL(SP), DIMENSION(SIZE(STATE_START)) :: EVEC        ! error estimate for each state
REAL(SP), DIMENSION(SIZE(STATE_START)) :: TVEC        ! error threshold for each state
REAL(SP)                               :: MULT        ! multiplier for new step size
REAL(SP), PARAMETER                    :: EPS=1.E-10_SP ! machine constant to prevent floating point errors
INTEGER(I4B), DIMENSION(1)             :: IMAX        ! index of maximum error
INTEGER(I4B)                           :: NITER       ! number of iterations in newtoniter
LOGICAL(LGT)                           :: CHECK       ! convergence check in SUBROUTINE newtoniter
LOGICAL(LGT)                           :: FEXCESS     ! FLAG to denote if states are corrected for excessive extrapolation
REAL(SP)                               :: TEMPSTEP    ! suggested new time step, for case of non-convergence
! -------------------------------------------------------------------------------------------------
INTERFACE
 SUBROUTINE MODL_SOLVE(CALCDSDT,IE_SOLVE,SI_SOLVE,B_IMPOSE,AVG_FLUX,ADD_FLUX,NEWSTATE, & ! define functionality of the routine
                       DT,S0,S1,DSDT,NEWSTEP,CONVCHECK,NITER,SOLUTION,HBOUND, &          ! input/output
                       IERR,MESSAGE)                                                     ! error control
 USE nrtype                                                   ! variable definitions, etc.
 IMPLICIT NONE
 LOGICAL(LGT), INTENT(IN),OPTIONAL             :: CALCDSDT    ! FLAG to calculate derivatives at S0
 LOGICAL(LGT), INTENT(IN),OPTIONAL             :: IE_SOLVE    ! FLAG to compute the implicit Euler solution
 LOGICAL(LGT), INTENT(IN),OPTIONAL             :: SI_SOLVE    ! FLAG to compute the semi-implicit Euler solution
 LOGICAL(LGT), INTENT(IN),OPTIONAL             :: B_IMPOSE    ! FLAG to impose bounds on model states
 LOGICAL(LGT), INTENT(IN),OPTIONAL             :: AVG_FLUX    ! FLAG to average fluxes from start & end states
 LOGICAL(LGT), INTENT(IN),OPTIONAL             :: ADD_FLUX    ! FLAG to add accepted fluxes to the total flux
 LOGICAL(LGT), INTENT(IN),OPTIONAL             :: NEWSTATE    ! FLAG to use weighted fluxes to compute end state
 REAL(SP), INTENT(IN), OPTIONAL                :: DT          ! length of the sub-step
 REAL(SP), DIMENSION(:),INTENT(IN), OPTIONAL   :: S0          ! input state vector
 REAL(SP), DIMENSION(:), INTENT(OUT),OPTIONAL  :: S1          ! state vector from the implicit euler solution
 REAL(SP), DIMENSION(:),INTENT(INOUT),OPTIONAL :: DSDT        ! state derivatives
 LOGICAL(LGT), INTENT(IN),OPTIONAL             :: NEWSTEP     ! FLAG to denote a new model time step
 LOGICAL(LGT), INTENT(IN),OPTIONAL             :: CONVCHECK   ! FLAG to check for convergence of the implicit scheme
 INTEGER(I4B), INTENT(OUT), OPTIONAL           :: NITER       ! number of iterations
 INTEGER(I4B), INTENT(IN), OPTIONAL            :: SOLUTION    ! solution is at start (0) or end (1) of sub-step
 LOGICAL(LGT), INTENT(IN),OPTIONAL             :: HBOUND      ! FLAG to denote if the states were out of bounds
 INTEGER(I4B), INTENT(OUT)                     :: IERR        ! error code
 CHARACTER(LEN=*), INTENT(OUT)                 :: MESSAGE     ! error message
 END SUBROUTINE MODL_SOLVE
END INTERFACE
! ---------------------------------------------------------------------------------------
! (0) INITIALIZATION
! ---------------------------------------------------------------------------------------
! intilize states and counters
ETIME         = 0._sp                        ! part of the time step completed
CHECK         = .FALSE.
STATE0        = STATE_START                  ! save model states at the start of the full step
STATE1_RETAIN = STATE_START                  ! initial state (needed for rejected steps)
newStep       = .true.                       ! initialize newstep (force re-calculation of Jacobian)
NEW_SUBSTEP   = .TRUE.                       ! initialize new sub-step (check if need new derivatives)
! initialize diagnostix
NUM_FUNCS     = 0                            ! number of function calls
NUM_JACOBIAN  = 0                            ! number of times Jacobian is calculated
NUMSUB_ACCEPT = 0                            ! number of sub-steps accepted (taken)
NUMSUB_REJECT = 0                            ! number of sub-steps tried but rejected
NUMSUB_NOCONV = 0                            ! number of sub-steps tried that did not converge
MAXNUM_ITERNS = 0                            ! maximum number of iterations taken in the newton method
! ---------------------------------------------------------------------------------------
! DT_SUB (sub-step length) is carried over from previous step; ensure that it is in bounds
DT_SUB        = MIN( MAX(MIN_TSTEP,DT_SUB), MAX_TSTEP) ! (MIN_TSTEP < stepsize < MAX_TSTEP)
PREVSTEP      = DT_SUB    ! initialize the previous time step (tracked to avoid using small interval at end of step)
STEP_INCREASE = .FALSE.   ! used to check if the final sub-step has been increased

SUBSTEPS: DO  ! continuous (recursive) loop over sub-steps

 ! ---------------------------------------------------------------------------------------
 ! (0) SAVE VECTOR OF STATES AND DERIVATIVES AT THE START OF THE SUB-STEP
 ! ---------------------------------------------------------------------------------------

 ! refresh model states at the start of the sub-step
 IF (NEW_SUBSTEP .AND. .NOT.newStep) STATE0 = STATE1_RETAIN

 ! calculate new derivatives
 IF (NEW_SUBSTEP) THEN
  CALL MODL_SOLVE(CALCDSDT=.TRUE.,S0=STATE0,DT=DT_SUB,DSDT=DYDT_0,SOLUTION=0,IERR=IERR,MESSAGE=MESSAGE)
  IF (IERR.NE.0) THEN; PRINT *, IERR, MESSAGE; STOP; ENDIF
 ELSE
  CALL MODL_SOLVE(CALCDSDT=.FALSE.,SOLUTION=0,IERR=IERR,MESSAGE=MESSAGE)
  IF (IERR.NE.0) THEN; PRINT *, IERR, MESSAGE; STOP; ENDIF
 ENDIF

 ! select solution method
 SELECT CASE(SOLUTION_METHOD)

  ! ---------------------------------------------------------------------------------------
  ! (1) CALCULATE EXPLICIT EULER SOLUTION -- NOTE, NO ERROR CONTROL
  ! ---------------------------------------------------------------------------------------
  CASE (EXPLICIT_EULER)
   ! calculate explicit Euler solution
   STATE1_HI = STATE0 + DYDT_0*DT_SUB        ! explicit solution (can be out of range, but OK for error control)
   ! get a safegaurded solution to account for excessive extrapolation (includes flux disaggregation)
   CALL MODL_SOLVE(B_IMPOSE=.TRUE.,S0=STATE1_HI,S1=STATE1_HI_S,DT=DT_SUB,HBOUND=FEXCESS,&
                   IERR=IERR,MESSAGE=MESSAGE)
   IF (IERR.NE.0) THEN; PRINT *, IERR, MESSAGE; STOP; ENDIF
   newStep=.false.

  ! --------------------------------------------------------------------------------------
  ! (2) CALCULATE IMPLICIT EULER SOLUTION -- NOTE, NO ERROR CONTROL
  ! --------------------------------------------------------------------------------------
  CASE (IMPLICIT_EULER)
   ! estimate the initial conditions used in the Newton scheme
   SELECT CASE (INITIAL_NEWTON)
    CASE (STATE_OLD);     STATE1_INIT = STATE0
    CASE (EXPLICIT_MID);  STATE1_INIT = STATE0 + DYDT_0*DT_SUB/2.0_SP ! estimate at mid-point
    CASE (EXPLICIT_FULL); STATE1_INIT = STATE0 + DYDT_0*DT_SUB        ! estimate at end
   END SELECT
   ! estimate state vector at end of time step
   CALL MODL_SOLVE(IE_SOLVE=.TRUE.,S0=STATE1_INIT,S1=STATE1_HI,DSDT=DYDT_1,DT=DT_SUB,&
                   NEWSTEP=newStep,CONVCHECK=CHECK,NITER=NITER,&
                   IERR=IERR,MESSAGE=MESSAGE)
   IF (IERR.NE.0) THEN; PRINT *, IERR, TRIM(MESSAGE); STOP; ENDIF
   IF (NITER > MAXNUM_ITERNS) MAXNUM_ITERNS=NITER
   newStep=.false.
   ! re-compute state vector at the end of the sub-step (needed for non-convergence)
   STATE1_HI_S = STATE0 + DYDT_1*DT_SUB

  ! ---------------------------------------------------------------------------------------
  ! (3) CALCULATE EXPLICIT HEUN SOLUTION
  ! ---------------------------------------------------------------------------------------
  CASE (EXPLICIT_HEUN)
   ! calculate explicit Euler solution
   STATE1_LO = STATE0 + DYDT_0*DT_SUB        ! explicit solution (can be out of range, but OK for error control)
   ! get a safegaurded solution to account for excessive extrapolation (includes flux disaggregation)
   CALL MODL_SOLVE(B_IMPOSE=.TRUE.,S0=STATE1_LO,S1=STATE1_LO_S,DT=DT_SUB,HBOUND=FEXCESS,&
                   IERR=IERR,MESSAGE=MESSAGE)
   IF (IERR.NE.0) THEN; PRINT *, IERR, TRIM(MESSAGE); STOP; ENDIF
   newStep=.false.
   ! calculate explicit Heun solution (NOTE: using safeguarded states)
   CALL MODL_SOLVE(CALCDSDT=.TRUE.,S0=STATE1_LO_S,DT=DT_SUB,DSDT=DYDT_1,SOLUTION=1,IERR=IERR,MESSAGE=MESSAGE)
   IF (IERR.NE.0) THEN; PRINT *, IERR, MESSAGE; STOP; ENDIF
   STATE1_HI = STATE0 + 0.5_SP*(DYDT_0+DYDT_1)*DT_SUB
   ! average fluxes (average fluxes before imposing bounds)
   CALL MODL_SOLVE(AVG_FLUX=.TRUE.,IERR=IERR,MESSAGE=MESSAGE)
   IF (IERR.NE.0) THEN; PRINT *, IERR, TRIM(MESSAGE); STOP; ENDIF
   ! get a safegaurded solution to account for excessive extrapolation (includes flux disaggregation)
   CALL MODL_SOLVE(B_IMPOSE=.TRUE.,S0=STATE1_HI,S1=STATE1_HI_S,DT=DT_SUB,HBOUND=FEXCESS,&
                   IERR=IERR,MESSAGE=MESSAGE)
   IF (IERR.NE.0) THEN; PRINT *, IERR, TRIM(MESSAGE); STOP; ENDIF

  ! --------------------------------------------------------------------------------------
  ! (4) CALCULATE IMPLICIT HEUN SOLUTION
  ! --------------------------------------------------------------------------------------
  CASE (IMPLICIT_HEUN)
   ! estimate the initial conditions used in the Newton scheme
   SELECT CASE (INITIAL_NEWTON)
    CASE (STATE_OLD);     STATE1_INIT = STATE0
    CASE (EXPLICIT_MID);  STATE1_INIT = STATE0 + DYDT_0*DT_SUB/2.0_SP ! estimate at mid-point
    CASE (EXPLICIT_FULL); STATE1_INIT = STATE0 + DYDT_0*DT_SUB        ! estimate at end
   END SELECT
   ! estimate state vector at end of sub-step
   CALL MODL_SOLVE(IE_SOLVE=.TRUE.,S0=STATE1_INIT,S1=STATE1_HI,DSDT=DYDT_1,DT=DT_SUB,&
                   NEWSTEP=newStep,CONVCHECK=CHECK,NITER=NITER,&
                   IERR=IERR,MESSAGE=MESSAGE)
   IF (IERR.NE.0) THEN; PRINT *, IERR, TRIM(MESSAGE); STOP; ENDIF
   IF (NITER > MAXNUM_ITERNS) MAXNUM_ITERNS=NITER
   newStep=.false.
   ! re-compute state vector at the end of the sub-step (needed for non-convergence)
   STATE1_HI = STATE0 + 0.5_SP*(DYDT_0+DYDT_1)*DT_SUB
   ! average fluxes (average fluxes before imposing bounds)
   CALL MODL_SOLVE(AVG_FLUX=.TRUE.,IERR=IERR,MESSAGE=MESSAGE)
   IF (IERR.NE.0) THEN; PRINT *, IERR, TRIM(MESSAGE); STOP; ENDIF
   ! get a safegaurded solution to account for excessive extrapolation (includes flux disaggregation)
   CALL MODL_SOLVE(B_IMPOSE=.TRUE.,S0=STATE1_HI,S1=STATE1_HI_S,DT=DT_SUB,HBOUND=FEXCESS,&
                   IERR=IERR,MESSAGE=MESSAGE)
   IF (IERR.NE.0) THEN; PRINT *, IERR, TRIM(MESSAGE); STOP; ENDIF
   ! check for non-convergence
   IF (TEMPORAL_ERROR_CONTROL.EQ.TS_ADAPT) THEN
    IF (CHECK) THEN
     NUMSUB_NOCONV = NUMSUB_NOCONV + 1
     STEP = MAX(MIN_TSTEP, DT_SUB*RMIN)  ! (avoid stepsize < MIN_TSTEP)
     TEMPSTEP = REVISE_STEP()  ! avoid small steps at the end of a time interval
     IF (TEMPSTEP.NE.STEP) THEN; PREVSTEP=STEP; ELSE; PREVSTEP=TEMPSTEP; ENDIF
     ! avoid the case of a continuous do loop where TEMPSTEP is at a minimum
     IF (TEMPSTEP.LT.DT_SUB) THEN ! TEMPSTEP may equal DT_SUB (MIN_TSTEP, or end of interval)
      newStep = .true.
      DT_SUB  = TEMPSTEP
      CYCLE SUBSTEPS
     ENDIF
     IERR=10; MESSAGE='newton did not converge, and unable to make steps small enough'; RETURN
    ENDIF
   ENDIF
   ! compute auxillary lower-order solution
   STATE1_LO = STATE0 + DYDT_1*DT_SUB

  ! --------------------------------------------------------------------------------------
  ! (5) CALCULATE SEMI-IMPLICIT EULER SOLUTION
  ! --------------------------------------------------------------------------------------
  CASE (SEMI_IMPLICIT)
   ! use explicit Euler for lower-order solution
   STATE1_LO = STATE0 + DYDT_0*DT_SUB
   ! estimate state vector at end of time step
   CALL MODL_SOLVE(SI_SOLVE=.TRUE.,S0=STATE0,S1=STATE1_HI,DSDT=DYDT_0,DT=DT_SUB,&
                   IERR=IERR,MESSAGE=MESSAGE)
   IF (IERR.LT.0) PRINT *, IERR, TRIM(MESSAGE)
   IF (IERR.GT.0) THEN; PRINT *, IERR, TRIM(MESSAGE); STOP; ENDIF
   newStep=.false.
   ! get a safegaurded solution to account for excessive extrapolation (includes flux disaggregation)
   CALL MODL_SOLVE(B_IMPOSE=.TRUE.,S0=STATE1_HI,S1=STATE1_HI_S,DT=DT_SUB,HBOUND=FEXCESS,&
                   IERR=IERR,MESSAGE=MESSAGE)


  ! check that the solution method is OK
  CASE DEFAULT
   IERR=20; MESSAGE='ode_int: unknown solution method'; RETURN

 END SELECT

 ! --------------------------------------------------------------------------------------
 ! (3) CALCULATE ERROR, CHECK IF ACCEPT/REJECT THE CURRENT STEP, AND NEW STEP SIZE
 ! --------------------------------------------------------------------------------------
 ! calculate the maximum error over all states
 NEW_SUBSTEP = .FALSE.
 IF (TEMPORAL_ERROR_CONTROL.EQ.TS_ADAPT) THEN
  EVEC = ABS(STATE1_HI_S - STATE1_LO)                   ! error estimate
  TVEC = ERR_TRUNC_REL*ABS(STATE1_HI_S) + ERR_TRUNC_ABS ! error thresholds
  IMAX = MAXLOC(EVEC - TVEC)                            ! index of maximum error
  IF (EVEC(IMAX(1)) < TVEC(IMAX(1)) .OR. &     ! (accept if error is less than critical threshold)
      DT_SUB <= MIN_TSTEP) THEN
    NEW_SUBSTEP = .TRUE.
  ENDIF
 ELSE
  EVEC = 0._SP; TVEC = 0._SP; IMAX = 0
  NEW_SUBSTEP = .TRUE. ! (accept if fixed time steps)
 ENDIF
 ! --------------------------------------------------------------------------------------
 IF (NEW_SUBSTEP) THEN                ! (accept if time step is already minimum allowable)
  !WRITE(*,'(I1,1X,2(F8.5,1X),I1,1X,20(F8.3,1X))') 0, ETIME, DT_SUB, IMAX, STATE1_HI_S, EVEC, TVEC
  ! accept step -- calculate new (increased) step size
  ! NOTE: step size not necessarily increased because of the safety factor
  IF (TEMPORAL_ERROR_CONTROL.EQ.TS_ADAPT) THEN
   MULT = SAFETY * SQRT( TVEC(IMAX(1)) / MAX(EVEC(IMAX(1)),EPS) )
   STEP = MIN( MAX(MIN_TSTEP, DT_SUB * MIN(MULT,RMAX) ), MAX_TSTEP) ! (MIN_TSTEP < stepsize < MAX_TSTEP)
  ELSE
   STEP = MAX_TSTEP
  ENDIF
  ! add contribution of sub-step flux to the timestep-average flux
  CALL MODL_SOLVE(ADD_FLUX=.TRUE.,S1=STATE1_HI_S,DT=DT_SUB,IERR=IERR,MESSAGE=MESSAGE)
  IF (IERR.NE.0) THEN; PRINT *, IERR, TRIM(MESSAGE); STOP; ENDIF
  STATE1_RETAIN = STATE1_HI_S

  NUMSUB_ACCEPT = NUMSUB_ACCEPT + 1
  ! compute fraction of big step that is finished, and check for exit criteria
  ETIME = ETIME + DT_SUB      ! identify position within the time step
  IF (ETIME.GE.DT_FULL) THEN
   EXIT SUBSTEPS              ! exit the substeps loop
  ENDIF
  ! revise the length of time steps to avoid small steps at the end of a time interval
  DT_SUB = REVISE_STEP()      ! avoid small steps at the end of a time interval
  IF (DT_SUB.NE.STEP) THEN; PREVSTEP=STEP; ELSE; PREVSTEP=DT_SUB; ENDIF
 ! --------------------------------------------------------------------------------------
 ELSE  ! REJECT STEP AND DECREASE STEP SIZE
  NEW_SUBSTEP = .FALSE.
  !WRITE(*,'(I1,1X,2(F8.5,1X),I1,1X,20(F8.3,1X))') 1, ETIME, DT_SUB, IMAX, STATE1_HI_S, EVEC, TVEC
  ! calculate new (decreased) step size
  NUMSUB_REJECT = NUMSUB_REJECT + 1
  MULT = SAFETY * SQRT( TVEC(IMAX(1)) / MAX(EVEC(IMAX(1)),EPS) )
  STEP = MAX(MIN_TSTEP, DT_SUB * MAX(MULT,RMIN) ) ! (avoid stepsize < MIN_TSTEP)
  DT_SUB = REVISE_STEP()  ! avoid small steps at the end of a time interval
  IF (DT_SUB.NE.STEP) THEN; PREVSTEP=STEP; ELSE; PREVSTEP=DT_SUB; ENDIF
 ENDIF

 ! (keep looping)
END DO SUBSTEPS ! continuous (recursive) do loop
!print *, 'num_funcs = ', num_funcs

! ---------------------------------------------------------------------------------------
! (9) RE-COMPUTE STATES AT THE END OF THE FULL STEP
! ---------------------------------------------------------------------------------------
! The implicit solution is not exact.  To conserve mass, we uses the weighted average of
! model fluxes throughout the time step to re-compute states at the end of the time step
! ---------------------------------------------------------------------------------------
! update model states (note use of DT_FULL)
CALL MODL_SOLVE(NEWSTATE=.TRUE.,S1=STATE_END,DT=DT_FULL,IERR=IERR,MESSAGE=MESSAGE)
IF (IERR.NE.0) THEN; PRINT *, IERR, TRIM(MESSAGE); STOP; ENDIF
! NOTE: may need to modify diagnostic variables that do not have time units, e.g., satarea = satarea/dt_full
DT_SUB=PREVSTEP                       ! ensure stepsize is not equal to the small remainder

! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
CONTAINS
 FUNCTION REVISE_STEP()
 REAL(SP)    :: REVISE_STEP
 REAL(SP)    :: T_MGN
 SELECT CASE(SMALL_ENDSTEP)
  ! -------------------------------------------------------------------------------------
  CASE(STEP_TRUNC)  ! truncate the time step if near the end
   IF (ETIME + STEP .GE. DT_FULL) REVISE_STEP = DT_FULL - ETIME
   IF (ETIME + STEP .LT. DT_FULL) REVISE_STEP = STEP
  ! -------------------------------------------------------------------------------------
  CASE(LOOK_AHEAD)  ! the look-ahead method of Shampine (1994)
   IF (ETIME + STEP .GE. DT_FULL) THEN
    REVISE_STEP = DT_FULL - ETIME
   ELSE
    IF (ETIME + STEP*2._SP .GE. DT_FULL) THEN
     REVISE_STEP = (DT_FULL - ETIME)/2._SP
    ELSE
     REVISE_STEP = STEP
    ENDIF
   ENDIF
  ! -------------------------------------------------------------------------------------
  CASE(STEP_ABSORB) ! the step-absorption method
   IF (STEP_INCREASE) THEN  ! only try and increase step size once
    IF (ETIME + STEP .GE. DT_FULL) REVISE_STEP = DT_FULL - ETIME
    IF (ETIME + STEP .LT. DT_FULL) REVISE_STEP = STEP
   ELSE
    T_MGN = STEP/SAFETY - STEP  ! margin of error
    IF (ETIME + STEP + T_MGN .GE. DT_FULL) THEN
     REVISE_STEP   = DT_FULL - ETIME
     STEP_INCREASE = .TRUE.
    ELSE
     IF (ETIME + STEP + T_MGN*2._SP .GE. DT_FULL) THEN
      REVISE_STEP   = STEP + T_MGN*(T_MGN/(DT_FULL-(ETIME+STEP)))
      STEP_INCREASE = .TRUE.
     ELSE
      REVISE_STEP = STEP
     ENDIF
    ENDIF
   ENDIF
  CASE DEFAULT; STOP ' must use the STEP_TRUNC, LOOK_AHEAD, or STEP_ABSORB methods '
 END SELECT
 END FUNCTION REVISE_STEP
END SUBROUTINE ODE_INT
