SUBROUTINE FUSE_SOLVE(CALCDSDT,IE_SOLVE,SI_SOLVE,B_IMPOSE,AVG_FLUX,ADD_FLUX,NEWSTATE, & ! define functionality of the routine
                      DT,S0,S1,DSDT,NEWSTEP,CONVCHECK,NITER,SOLUTION,HBOUND, &          ! input/output
                      IERR,MESSAGE)                                                     ! error control
! Used to
!  (1) calculate dS/dt for the input vector S0
!  (2) solve for S using the implicit Euler method
!  (3) solve for S using the semi-implicit Euler method
!  (4) average fluxes from the start and end of the sub-step
!  (5) impose bounds on model states (and disaggregate fluxes) 
!  (6) add fluxes from accepted sub-steps to the total timestep flux
!  (7) estimate state at end of a full step, based on sum of fluxes 
USE nrtype                                                   ! variable definitions, etc.
USE multi_flux, ONLY: M_FLUX,FLUX_0,FLUX_1,W_FLUX,&          ! model fluxes
                       CURRENT_DT                            ! model fluxes (continued)
USE multistate, ONLY: FSTATE,MSTATE,BSTATE,ESTATE,&          ! model states
                       DY_DT,DYDT_0,DYDT_1,HSTATE            ! model states (continued)
USE fminln, ONLY: fmin_x0p,fmin_dtp,fmin_dt2p,fmin_dseep     ! variables used for residual vector in IE 
USE xtry_2_str_module                                        ! provide access to xtry_2_str
USE str_2_xtry_module                                        ! provide access to str_2_xtry
USE fuse_deriv_module                                        ! provide access to derivatives
USE fuse_sieul_module                                        ! provide access to the semi-implicit Euler function
USE newtoniter_mod                                           ! provide access to newtoniter
IMPLICIT NONE
! input/output variables
LOGICAL(LGT), INTENT(IN),OPTIONAL             :: CALCDSDT    ! FLAG to calculate derivatives at S0
LOGICAL(LGT), INTENT(IN),OPTIONAL             :: IE_SOLVE    ! FLAG to compute the implicit Euler solution
LOGICAL(LGT), INTENT(IN),OPTIONAL             :: SI_SOLVE    ! FLAG to compute the semi-implicit Euler solution
LOGICAL(LGT), INTENT(IN),OPTIONAL             :: B_IMPOSE    ! FLAG to impose bounds on model state
LOGICAL(LGT), INTENT(IN),OPTIONAL             :: AVG_FLUX    ! FLAG to average fluxes from start & end states
LOGICAL(LGT), INTENT(IN),OPTIONAL             :: ADD_FLUX    ! FLAG to add accepted fluxes to the total flux
LOGICAL(LGT), INTENT(IN),OPTIONAL             :: NEWSTATE    ! FLAG to use weighted fluxes to compute end state
REAL(SP), INTENT(IN), OPTIONAL                :: DT          ! length of the sub-step
REAL(SP), DIMENSION(:),INTENT(IN), OPTIONAL   :: S0          ! input state vector
REAL(SP), DIMENSION(:), INTENT(OUT),OPTIONAL  :: S1          ! state vector from the implicit euler solution
REAL(SP), DIMENSION(:),INTENT(INOUT),OPTIONAL :: DSDT        ! state derivatives
LOGICAL(LGT), INTENT(IN),OPTIONAL             :: NEWSTEP     ! FLAG to denote a new model time step
LOGICAL(LGT), INTENT(OUT),OPTIONAL            :: CONVCHECK   ! FLAG to check for convergence of the implicit scheme
INTEGER(I4B), INTENT(OUT), OPTIONAL           :: NITER       ! number of iterations
INTEGER(I4B), INTENT(IN), OPTIONAL            :: SOLUTION    ! solution is at start (0) or end (1) of sub-step
LOGICAL(LGT), INTENT(OUT),OPTIONAL            :: HBOUND      ! FLAG to denote if the states were out of bounds
INTEGER(I4B), INTENT(OUT)                     :: IERR        ! error code
CHARACTER(LEN=*), INTENT(OUT)                 :: MESSAGE     ! error message
! internal variables
REAL(SP), PARAMETER                           :: XACC=1.E-10 ! accuracy of implicit estimate 
LOGICAL(LGT)                                  :: ERROR_FLAG  ! FLAG to denote if violated constraints
REAL(SP), TARGET                              :: DT1         ! full time step
REAL(SP), TARGET                              :: DT2         ! half time step
REAL(SP), DIMENSION(:), ALLOCATABLE, TARGET   :: XI          ! initial state vector
REAL(SP), DIMENSION(:), ALLOCATABLE, TARGET   :: DSEE        ! change in state by explicit euler
REAL(SP), DIMENSION(:), ALLOCATABLE           :: DSDT0       ! state derivative at start of step
REAL(SP), DIMENSION(:), ALLOCATABLE           :: DSDT_SIE    ! state derivative for semi-implicit euler
! ---------------------------------------------------------------------------------------
IERR=0; MESSAGE='fuse_solve, just started'
! ---------------------------------------------------------------------------------------
! (1) CALCULATE DERIVATIVES
! ---------------------------------------------------------------------------------------
IF (PRESENT(CALCDSDT)) THEN
 IF (CALCDSDT) THEN
  ! check that we have passed what we need
  IF (.NOT.PRESENT(S0) .OR. .NOT.PRESENT(DT) .OR. .NOT.PRESENT(DSDT) .OR. .NOT.PRESENT(SOLUTION) ) THEN
   IF (.NOT.PRESENT(S0))       MESSAGE='need S0       to calculate model derivatives'
   IF (.NOT.PRESENT(DT))       MESSAGE='need DT       to calculate model derivatives'
   IF (.NOT.PRESENT(DSDT))     MESSAGE='need DSDT     to calculate model derivatives'
   IF (.NOT.PRESENT(SOLUTION)) MESSAGE='need SOLUTION to calculate model derivatives'
   IERR=20; RETURN
  ENDIF
  ! put DT into model flux structures
  CURRENT_DT = DT
  ! calculate derivatives
  DSDT = FUSE_DERIV(S0)   ! calculate derivatives
  ! save information in model structures
  SELECT CASE(SOLUTION)
   CASE(0)
    FLUX_0 = M_FLUX       ! save fluxes at the start of the sub-step
    DYDT_0 = DY_DT        ! save derivatives at the start of the sub-step
   CASE(1)
    FLUX_1 = M_FLUX       ! save fluxes at the end of the sub-step
    DYDT_1 = DY_DT        ! save derivatives at the start of the sub-step
  END SELECT
 ELSE
  ! check that we have passed what we need
  IF (.NOT.PRESENT(SOLUTION)) THEN
   MESSAGE='need SOLUTION to calculate model derivatives'; IERR=20; RETURN
  ENDIF
  ! extract information from model structures
  SELECT CASE(SOLUTION)
   CASE(0)
    M_FLUX = FLUX_0     ! extract fluxes from the start of the sub-step
    DY_DT  = DYDT_0     ! extract derivatives from the start of the sub-step
   CASE(1)
    M_FLUX = FLUX_1     ! extract fluxes from the end of the sub-step
    DY_DT  = DYDT_1     ! extract derivatives from the start of the sub-step
  END SELECT
 ENDIF
ENDIF 
! ---------------------------------------------------------------------------------------
! (2) ESTIMATE NEW VECTOR OF STATES USING THE IMPLICIT EULER/HEUN METHOD
! ---------------------------------------------------------------------------------------
IF (PRESENT(IE_SOLVE)) THEN
 IF (IE_SOLVE) THEN
  ! check that we have passed what we need
  IF (.NOT.PRESENT(S0) .OR. .NOT.PRESENT(S1) .OR. .NOT.PRESENT(DSDT) .OR. .NOT.PRESENT(DT) .OR. &
      .NOT.PRESENT(NEWSTEP) .OR. .NOT.PRESENT(CONVCHECK) .OR. .NOT.PRESENT(NITER)) THEN
   IF (.NOT.PRESENT(S0))        MESSAGE='need S0 for the implicit euler solution'
   IF (.NOT.PRESENT(S1))        MESSAGE='need S1 for the implicit euler solution'
   IF (.NOT.PRESENT(DT))        MESSAGE='need DT for the implicit euler solution'
   IF (.NOT.PRESENT(DSDT))      MESSAGE='need DYDT for the implicit euler solution'
   IF (.NOT.PRESENT(NEWSTEP))   MESSAGE='need NEWSTEP for the implicit euler solution'
   IF (.NOT.PRESENT(CONVCHECK)) MESSAGE='need CONVCHECK for the implicit euler solution'
   IF (.NOT.PRESENT(NITER))     MESSAGE='need NITER for the implicit euler solution'
   IERR=20; RETURN
  ENDIF
  ! alolocate space for pointer targets
  allocate(xi(size(s0)),dsee(size(s0)),dsdt0(size(s0)), stat=ierr)
  if (ierr.ne.0) then; ierr=20; message='fuse_solve: problem allocating space'; endif
  ! make pointer assignments for initial state and time steps (used in fminln for calc residual vector)
  fmin_x0p  =>xi                ! provide access to the initial state used in fmin
  fmin_dtp  =>dt1               ! provide access to the time step used in fmin
  fmin_dt2p =>dt2               ! provide access to the half time step used in fmin
  fmin_dseep=>dsee              ! provide access to the vector of change in state by explicit euler
  ! put DT into model flux structures
  CURRENT_DT = DT
  ! populate targets
  DT1=DT                        ! full sub-step
  DT2=DT/2._SP                  ! half sub-step
  CALL STR_2_XTRY(MSTATE,XI)    ! retrieve state at the start of the sub-step
  CALL STR_2_XTRY(DYDT_0,DSDT0) ! retrieve derivatives at the start of the sub-step
  DSEE = DSDT0*DT2              ! calculate explicit euler component of Heun solution
  ! compute the IE solution
  S1     = S0                                 ! S1 over-written on output
  CALL NEWTONITER(S1,NEWSTEP,CONVCHECK,NITER) ! try different values of X until converge
  FLUX_1 = M_FLUX                             ! save fluxes at end of sub-step (save in model structure)
  DYDT_1 = DY_DT                              ! save derivs at end of sub-step (save in model structure)
  CALL STR_2_XTRY(DY_DT,DSDT)                 ! extract derivatives from model structure, and return to ODE_INT
  ! deallocate space for pointer targets
  deallocate(xi,dsee,dsdt0, stat=ierr)
  if (ierr.ne.0) then; ierr=20; message='fuse_solve: problem deallocating space'; endif
 ENDIF
ENDIF
! ---------------------------------------------------------------------------------------
! (3) ESTIMATE NEW VECTOR OF STATES USING THE SEMI-IMPLICIT EULER METHOD
! ---------------------------------------------------------------------------------------
IF (PRESENT(SI_SOLVE)) THEN
 IF (SI_SOLVE) THEN
  ! check that we have passed what we need
  IF (.NOT.PRESENT(S0) .OR. .NOT.PRESENT(S1) .OR. .NOT.PRESENT(DSDT) .OR. .NOT.PRESENT(DT)) THEN
   IF (.NOT.PRESENT(S0))        MESSAGE='need S0   for the semi-implicit euler solution'
   IF (.NOT.PRESENT(S1))        MESSAGE='need S1   for the semi-implicit euler solution'
   IF (.NOT.PRESENT(DSDT))      MESSAGE='need DSDT for the semi-implicit euler solution'
   IF (.NOT.PRESENT(DT))        MESSAGE='need DT   for the semi-implicit euler solution'
   IERR=20; RETURN
  ENDIF
  ! allocate space
  ALLOCATE(DSDT_SIE(SIZE(S0)), STAT=IERR)
  IF (IERR.NE.0) THEN; IERR=20; MESSAGE='fuse_solve: problem allocating space'; ENDIF
  ! put DT into model flux structures
  CURRENT_DT = DT
  ! estimate new derivatives using the semi-implicit method
  CALL FUSE_SIEUL(S0,DSDT,DT,IERR,MESSAGE)   ! somewhat FUSE-specific
  CALL STR_2_XTRY(DY_DT,DSDT_SIE)  ! extract derivatives from the FUSE data structures
  ! compute new state
  S1 = S0 + DSDT_SIE*DT
  ! deallocate space
  DEALLOCATE(DSDT_SIE, STAT=IERR)
  IF (IERR.NE.0) THEN; IERR=20; MESSAGE='fuse_solve: problem deallocating space'; ENDIF
 ENDIF
ENDIF 
! ---------------------------------------------------------------------------------------
! (4) AVERAGE FLUXES FROM START & END OF STEP (NECESSARY IF ACCEPT HIGHER ORDER SOLUTION)
! ---------------------------------------------------------------------------------------
IF (PRESENT(AVG_FLUX)) THEN
 IF (AVG_FLUX) THEN   ! Case 1: Higher-order solution accepted
  ! average fluxes and derivatives from the start and end of the step
  CALL MEANFLUXES()
 ELSE                 ! Case 2: Lower-order solution accepted
  ! check that the solution argument is present
  IF (.NOT.PRESENT(SOLUTION)) THEN
   MESSAGE='need SOLUTION to assign fluxes and derivatives'; IERR=20; RETURN
  ENDIF
  ! assign fluxes from the appropriate solution
  SELECT CASE(SOLUTION)
   CASE(0) ! explicit euler: save fluxes and derivatives at start of sub-step
    M_FLUX = FLUX_0
    DY_DT  = DYDT_0
   CASE(1) ! implicit euler: save fluxes and derivatives at end of sub-step
    M_FLUX = FLUX_1
    DY_DT  = DYDT_1
  END SELECT
 ENDIF
ENDIF
! ---------------------------------------------------------------------------------------
! (5) IMPOSE BOUNDS ON MODEL STATES (AND DISAGGREGATE FLUXES)
! --------------------------------------------------------------------------------------- 
IF (PRESENT(B_IMPOSE)) THEN
 IF (B_IMPOSE) THEN
  ! check that we have passed what we need
  IF (.NOT.PRESENT(S0) .OR. .NOT.PRESENT(S1) .OR. .NOT.PRESENT(DT) .OR. &
      .NOT.PRESENT(HBOUND)) THEN
   IF (.NOT.PRESENT(S0))     MESSAGE='need S0     to impose bounds on model states'
   IF (.NOT.PRESENT(S1))     MESSAGE='need S1     to impose bounds on model states'
   IF (.NOT.PRESENT(DT))     MESSAGE='need DT     to impose bounds on model states'
   IF (.NOT.PRESENT(HBOUND)) MESSAGE='need HBOUND to impose bounds on model states'
   IERR=20; RETURN
  ENDIF
  ! put the model states in the appropriate structures
  BSTATE = MSTATE                 ! state at the start of the sub-step
  CALL XTRY_2_STR(S0,ESTATE)      ! extrapolated state at the end of the sub-step
  ! constrain bounds
  CALL FIX_STATES(DT,ERROR_FLAG)  ! ERROR_FLAG is a logical flag to denote if hit bound
  HBOUND=ERROR_FLAG
  ! extract states from the model structure
  CALL STR_2_XTRY(ESTATE,S1)      ! corrected state at the end of the sub-step
 ENDIF
ENDIF
! ---------------------------------------------------------------------------------------
! (6) ADD FLUXES FROM ACCEPTED SUB-STEPS TO THE TOTAL TIMESTEP FLUX
! ---------------------------------------------------------------------------------------
IF (PRESENT(ADD_FLUX)) THEN
 IF (ADD_FLUX) THEN
  ! check that S1 and DT are present
  IF (.NOT.PRESENT(S1) .OR. .NOT.PRESENT(DT)) THEN
   IF (.NOT.PRESENT(S1)) MESSAGE='need S1 to aggregate fluxes and save states'
   IF (.NOT.PRESENT(DT)) MESSAGE='need DT to aggregate fluxes and save states'
   IERR=20; RETURN
  ENDIF
  ! aggregate fluxes and save states
  HSTATE%STEP = DT            ! insert the time interval into the data structures
  CALL WGT_FLUXES()           ! compute the contribution of the flux over the time interval DT
  CALL XTRY_2_STR(S1,MSTATE)  ! update MSTATE
 ENDIF
ENDIF
! ---------------------------------------------------------------------------------------
! (7) COMPUTE STATE AT THE END OF THE TIME INTERVAL
! ---------------------------------------------------------------------------------------
IF (PRESENT(NEWSTATE)) THEN
 ! check that S1 and DT are present
 IF (.NOT.PRESENT(S1) .OR. .NOT.PRESENT(DT)) THEN
  IF (.NOT.PRESENT(S1)) MESSAGE='need S1 to aggregate fluxes and save states'
  IF (.NOT.PRESENT(DT)) MESSAGE='need DT to aggregate fluxes and save states'
  IERR=20; RETURN
 ENDIF
 ! update state
 IF (NEWSTATE) THEN
  M_FLUX = W_FLUX; CALL MSTATE_EQN()  ! compute model derivatives using aggregated fluxes
  CALL UPDATSTATE(DT)                 ! compute new value of FSTATE
  MSTATE = FSTATE                     ! update MSTATE
  CALL STR_2_XTRY(FSTATE,S1)          ! extract state vector
 ENDIF
ENDIF
! ---------------------------------------------------------------------------------------
END SUBROUTINE FUSE_SOLVE
