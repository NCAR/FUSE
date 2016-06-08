!MODULE TEST_SOLVE__MODULE
!IMPLICIT NONE
!CONTAINS
! ---------------------------------------------------------------------------------------
SUBROUTINE TEST_SOLVE(CALCDSDT,IE_SOLVE,B_IMPOSE,AVG_FLUX,ADD_FLUX,NEWSTATE, & ! define functionality of the routine
                      DT,S0,S1,DSDT,NEWSTEP,CONVCHECK,NITER,SOLUTION,HBOUND, & ! input/output
                      IERR,MESSAGE)                                            ! error control
! Used to
!  (1) calculate dS/dt for the input vector S0
!  (2) solve for S using the implicit Euler method
!  (3) add fluxes from accepted sub-steps to the total timestep flux
USE nrtype                                                   ! variable definitions, etc.
USE test_modvar, ONLY : DT_SUB,MS_MIN,MS_MAX,TSTATE,CSTATE,& ! model variables
                        M_FLUX,FLUX_0,FLUX_1,W_FLUX,&        ! model variables (continued)
                        MDS_DT,DSDT_0,DSDT_1,&               ! model variables (continued)
                        MSTATE,FSTATE                        ! model variables (continued)
USE test_deriv__module                                       ! provide access to derivatives
IMPLICIT NONE
! input/output variables
LOGICAL(LGT), INTENT(IN),OPTIONAL             :: CALCDSDT    ! FLAG to calculate derivatives at S0
LOGICAL(LGT), INTENT(IN),OPTIONAL             :: IE_SOLVE    ! FLAG to compute the implicit Euler solution
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
REAL(SP)                                      :: ERROR_LOSS  ! extrapolation error
REAL(SP)                                      :: TOTAL_FLUX  ! total fluxes involved in extrapolation
! ---------------------------------------------------------------------------------------
INTERFACE
 SUBROUTINE IMPL_ERROR(S,F,DF)
 ! Calculates the error for the implicit scheme (used in RTNEWT_SUB)
 !  S(n+1) = S(n) + dS(n+1)/dt * delT
 !  F = S(try) - (S(n) + dS(try)/dt * delT)
 USE nrtype                                          ! numerical recipes data types
 IMPLICIT NONE
 REAL(SP), INTENT(IN)                :: S            ! storage
 REAL(SP), INTENT(OUT)               :: F            ! function value
 REAL(SP), INTENT(OUT)               :: DF           ! function derivative
 END SUBROUTINE IMPL_ERROR
END INTERFACE
! ---------------------------------------------------------------------------------------
IERR=0; MESSAGE='test_solve, just started'
! ---------------------------------------------------------------------------------------
! (1) CALCULATE DERIVATIVES
! ---------------------------------------------------------------------------------------
IF (PRESENT(CALCDSDT)) THEN
 ! check that we have passed what we need
 IF (.NOT.PRESENT(S0) .OR. .NOT.PRESENT(DSDT) .OR. .NOT.PRESENT(SOLUTION) ) THEN
  IF (.NOT.PRESENT(S0))       MESSAGE='need S0       to calculate model derivatives'
  IF (.NOT.PRESENT(DSDT))     MESSAGE='need DSDT     to calculate model derivatives'
  IF (.NOT.PRESENT(SOLUTION)) MESSAGE='need SOLUTION to calculate model derivatives'
  IERR=20; RETURN
 ENDIF
 ! calculate derivatives
 IF (CALCDSDT) DSDT = TEST_DERIV(S0)   ! calculate derivatives
 ! save information in model structures
 SELECT CASE(SOLUTION)
  CASE(0)
   FLUX_0 = M_FLUX    ! save fluxes at the start of the sub-step
   DSDT_0 = MDS_DT    ! save derivatives at the start of the sub-step
  CASE(1)
   FLUX_1 = M_FLUX    ! save fluxes at the end of the sub-step
   DSDT_1 = MDS_DT    ! save derivatives at the start of the sub-step
 END SELECT
ENDIF 
! ---------------------------------------------------------------------------------------
! (2) ESTIMATE NEW VECTOR OF STATES USING THE IMPLICIT EULER METHOD
! ---------------------------------------------------------------------------------------
IF (PRESENT(IE_SOLVE)) THEN
 IF (IE_SOLVE) THEN
  ! check that we have passed what we need
  IF (.NOT.PRESENT(S1) .OR. .NOT.PRESENT(DSDT) .OR. .NOT.PRESENT(DT) .OR. &
      .NOT.PRESENT(NEWSTEP) .OR. .NOT.PRESENT(CONVCHECK) .OR. .NOT.PRESENT(NITER)) THEN
   IF (.NOT.PRESENT(S1))        MESSAGE='need S1 for the implicit euler solution'
   IF (.NOT.PRESENT(DT))        MESSAGE='need DT for the implicit euler solution'
   IF (.NOT.PRESENT(DSDT))      MESSAGE='need DYDT for the implicit euler solution'
   IF (.NOT.PRESENT(NEWSTEP))   MESSAGE='need NEWSTEP for the implicit euler solution'
   IF (.NOT.PRESENT(CONVCHECK)) MESSAGE='need CONVCHECK for the implicit euler solution'
   IF (.NOT.PRESENT(NITER))     MESSAGE='need NITER for the implicit euler solution'
   IERR=20; RETURN
  ENDIF
  ! compute the IE solution
  DT_SUB    = DT            ! DT_SUB is stored in MODULE test_modvar
  CALL RTNEWT_SUB(IMPL_ERROR,S0(1),MS_MIN%WATR_1,MS_MAX%WATR_1,XACC,S1(1),NITER)
  FLUX_1    = M_FLUX        ! save fluxes at the end of the sub-step (save in model structure)
  DSDT_1    = MDS_DT        ! save derivatives at the end of the sub-step (save in model structure)
  DSDT(1)   = MDS_DT%WATR_1 ! extract derivatives from model structure (return to ODE_INT routine)
  CONVCHECK = .FALSE.       ! no check for convergence
 ENDIF
ENDIF
! ---------------------------------------------------------------------------------------
! (3) AVERAGE FLUXES FROM START & END OF STEP (NECESSARY IF ACCEPT HIGHER ORDER SOLUTION)
! ---------------------------------------------------------------------------------------
IF (PRESENT(AVG_FLUX)) THEN
 IF (AVG_FLUX) THEN   ! Case 1: Higher-order solution accepted
  ! average fluxes and derivatives from the start and end of the step
  M_FLUX%DRAINAGE = (FLUX_0%DRAINAGE + FLUX_1%DRAINAGE)/2._SP
  MDS_DT%WATR_1   = (DSDT_0%WATR_1 + DSDT_1%WATR_1)/2._SP
 ELSE                 ! Case 2: Lower-order solution accepted
  ! check that the solution argument is present
  IF (.NOT.PRESENT(SOLUTION)) THEN
   MESSAGE='need SOLUTION to assign fluxes and derivatives'; IERR=20; RETURN
  ENDIF
  ! assign fluxes from the appropriate solution
  SELECT CASE(SOLUTION)
   CASE(0) ! explicit euler: save fluxes and derivatives at start of sub-step
    M_FLUX = FLUX_0
    MDS_DT = DSDT_0
   CASE(1) ! implicit euler: save fluxes and derivatives at end of sub-step
    M_FLUX = FLUX_1
    MDS_DT = DSDT_1
  END SELECT
 ENDIF
ENDIF
! ---------------------------------------------------------------------------------------
! (4) IMPOSE BOUNDS ON MODEL STATES (AND DISAGGREGATE FLUXES)
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
  S1     = S0        ! get copy of S0
  HBOUND = .FALSE.   ! initialize bounds
  ! only need to constrain minimum
  IF (S1(1).LT.MS_MIN%WATR_1) THEN
   ERROR_LOSS      = (S1(1) - MS_MIN%WATR_1)/DT   ! error (L/T)
   TOTAL_FLUX      = M_FLUX%DRAINAGE              ! total flux (L/T)
   M_FLUX%DRAINAGE = M_FLUX%DRAINAGE + (M_FLUX%DRAINAGE/TOTAL_FLUX)*ERROR_LOSS
   S1(1)           = MS_MIN%WATR_1
   HBOUND          = .TRUE.
   print *, 'dude, hit zee bounds'
  ENDIF
 ENDIF
ENDIF
! ---------------------------------------------------------------------------------------
! (5) ADD FLUXES FROM ACCEPTED SUB-STEPS TO THE TOTAL TIMESTEP FLUX
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
  W_FLUX%DRAINAGE = W_FLUX%DRAINAGE + M_FLUX%DRAINAGE*DT
  W_FLUX%CHECKTIM = W_FLUX%CHECKTIM + DT
  MSTATE%WATR_1   = S1(1)
 ENDIF
ENDIF
! ---------------------------------------------------------------------------------------
! (6) COMPUTE STATE AT THE END OF THE TIME INTERVAL
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
  MDS_DT%WATR_1   = -W_FLUX%DRAINAGE
  FSTATE%WATR_1   = FSTATE%WATR_1 + MDS_DT%WATR_1*DT 
  MSTATE%WATR_1   = FSTATE%WATR_1
  S1(1)           = FSTATE%WATR_1
  print *, 'newstate ', S1(1), W_FLUX%CHECKTIM, W_FLUX%DRAINAGE, DT
 ENDIF
ENDIF
! ---------------------------------------------------------------------------------------
END SUBROUTINE TEST_SOLVE
!END MODULE TEST_SOLVE__MODULE
