PROGRAM driver_testfunc
! Used to test the routine for temporal integration of ordinary differential equations
!  with the test function dS/dt = -sqrt(S)
USE nrtype                                        ! numerical recipes data types
USE interfaceb, ONLY: ODE_INT,TEST_SOLVE
USE test_modvar, ONLY: MS_MIN,MS_MAX,MSTATE,&     ! model variables
                       FSTATE,W_FLUX              ! model variables (continued)
USE model_numerix, ONLY: SOLUTION_METHOD,ERR_TRUNC_ABS,ERR_TRUNC_REL
IMPLICIT NONE
REAL(SP), DIMENSION(1)             :: STATE0      ! initial state
REAL(SP), DIMENSION(1)             :: STATE1      ! final state
REAL(SP)                           :: DT_SUB      ! length of sub-step
REAL(SP)                           :: DT_FULL     ! length of full step
INTEGER(I4B)                       :: IERR        ! error code
CHARACTER(LEN=1024)                :: MESSAGE     ! error message
INTEGER(I4B)                       :: ITIME       ! time index (loop through time steps)
INTEGER(I4B), PARAMETER            :: NTIME=5     ! number of time steps
INTEGER(I4B)                       :: ITNC        ! loop through truncation errors
INTEGER(I4B)                       :: ITYP        ! loop through solution methods
CHARACTER(LEN=4)                   :: CH          ! character string for output
CHARACTER(LEN=14)                  :: CTYP        ! character string for output
! --------------------------------------------------------------------------------------
! define numerical solution methods
CALL DEFAULT_NUMERIX()
DO ITYP=0,1
 SOLUTION_METHOD=ITYP
 IF (ITYP.EQ.0) CTYP='EXPLICIT_EULER'
 IF (ITYP.EQ.1) CTYP='IMPLICIT_EULER'
 DO ITNC=1,4
  IF (ITNC.EQ.1) THEN; ERR_TRUNC_ABS = 1.e-1; ERR_TRUNC_REL = 1.e-1; CH='_e-1'; ENDIF
  IF (ITNC.EQ.2) THEN; ERR_TRUNC_ABS = 1.e-2; ERR_TRUNC_REL = 1.e-2; CH='_e-2'; ENDIF
  IF (ITNC.EQ.3) THEN; ERR_TRUNC_ABS = 1.e-3; ERR_TRUNC_REL = 1.e-3; CH='_e-3'; ENDIF
  IF (ITNC.EQ.4) THEN; ERR_TRUNC_ABS = 1.e-4; ERR_TRUNC_REL = 1.e-4; CH='_e-4'; ENDIF
  ! initialize variables
  STATE0         =     1._SP   ! state at the start of the time step
  STATE1         = -9999._SP   ! state at the end of the time step
  MS_MIN%WATR_1  = 1.E-10_SP   ! minimum values of model states (shared in MODULE test_modvar)
  MS_MAX%WATR_1  = 1.E+00_SP   ! maximum values of model states (shared in MODULE test_modvar)
  FSTATE%WATR_1  =     1._SP   ! initial value of model states (shared in MODULE test_modvar)
  DT_SUB         =     1._SP   ! length of sub-step
  DT_FULL        =     1._SP   ! length of full step
  print *, '*********************************************************************************'
  print *, '*********************************************************************************'
  print *, 'in driver', state0, ityp, itnc
  DO ITIME=1,NTIME
   ! open files
   if (itime.eq.1) open(21,file=CTYP//'1'//CH//'.dat',status='unknown')
   if (itime.eq.2) open(21,file=CTYP//'2'//CH//'.dat',status='unknown')
   if (itime.eq.3) open(21,file=CTYP//'3'//CH//'.dat',status='unknown')
   if (itime.eq.4) open(21,file=CTYP//'4'//CH//'.dat',status='unknown')
   if (itime.eq.5) open(21,file=CTYP//'5'//CH//'.dat',status='unknown')
   ! initialize states and fluxes
   MSTATE%WATR_1   = FSTATE%WATR_1
   W_FLUX%DRAINAGE = 0._SP 
   W_FLUX%CHECKTIM = 0._SP 
   ! temporally integrate the ode
   CALL ODE_INT(TEST_SOLVE,STATE0,STATE1,DT_SUB,DT_FULL,IERR,MESSAGE)



   STATE0          = STATE1
   print *, '***** in driver *****', itime
   close(21)
  END DO
 END DO
END DO
STOP
END PROGRAM DRIVER_TESTFUNC
! --------------------------------------------------------------------------------------
SUBROUTINE DEFAULT_NUMERIX()
USE model_numerix
SOLUTION_METHOD        = IMPLICIT_EULER        ! implicit euler solution
TEMPORAL_ERROR_CONTROL = TS_ADAPT              ! adaptive time steps
TRUNCATION_ERROR       = EMBEDDED_ERR          ! embedded error control
ORDER_ACCEPT           = HIGHER_ORDER          ! accept higher-order solutions
INITIAL_NEWTON         = EXPLICIT_FULL         ! initial conditions for Newton
JAC_RECOMPUTE          = FULLYVARIABLE         ! fully variable Jacobian
CHECK_OVERSHOOT        = LINE_SEARCH           ! use line search to trap/fix overshoot problems
ERR_TRUNC_ABS          = 1.e-3                 ! absolute temporal truncation error tolerance
ERR_TRUNC_REL          = 1.e-3                 ! relative temporal truncation error tolerance
ERR_ITER_FUNC          = 1.e-9                 ! iteration convergence tolerance for function values
ERR_ITER_DX            = 1.e-9                 ! iteration convergence tolerance for dx
FRACSTATE_MIN          = 1.e-9                 ! fractional minimum value of state (for non-zero derivatives)
SAFETY                 = 0.9_sp                ! safety factor in step-size equation
RMIN                   = 0.1_sp                ! minimum step size multiplier
RMAX                   = 4.0_sp                ! maximum step size multiplier
NITER_TOTAL            = 100                   ! total number of iterations used in the implicit scheme
MIN_TSTEP              = 0.01_sp/60._sp/24._sp ! minimum time step length (minutes --> days)
MAX_TSTEP              = 1440._sp/60._sp/24._sp ! maximum time step length (minutes --> days)
END SUBROUTINE DEFAULT_NUMERIX
