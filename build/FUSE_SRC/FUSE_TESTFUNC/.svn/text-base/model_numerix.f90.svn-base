!******************************************************************
MODULE model_numerix
! Purpose: To define method/parameters used for numerical solution
! Programmer: Dmitri Kavetski and Martyn Clark
! Last modified:
! Comments:
USE nrtype
implicit none
! ---------------------------------------------------------------------------------------
! (A) METHODS
! ---------------------------------------------------------------------------------------
! 1. Solution technique
INTEGER(I4B), PARAMETER    :: EXPLICIT_EULER=0, IMPLICIT_EULER=1
INTEGER(I4B)               :: SOLUTION_METHOD
! 2. Temporal error control
INTEGER(I4B), PARAMETER    :: TS_FIXED=0, TS_ADAPT=1
INTEGER(I4B)               :: TEMPORAL_ERROR_CONTROL
! 3. Method used to estimate temporal truncation error
INTEGER(I4B), PARAMETER    :: STEP_HALVING=0, EMBEDDED_ERR=1
INTEGER(I4B)               :: TRUNCATION_ERROR
! 4. Order of solution that is accepted
INTEGER(I4B), PARAMETER    :: HIGHER_ORDER=0, LOWER_ORDER=1
INTEGER(I4B)               :: ORDER_ACCEPT
! 5. Method used to estimate the initial conditions for the Newton scheme
INTEGER(I4B), PARAMETER    :: STATE_OLD=0, EXPLICIT_MID=1, EXPLICIT_FULL=2
INTEGER(I4B)               :: INITIAL_NEWTON 
! 6. Jacobian re-evaluation strategy
INTEGER(I4B), PARAMETER    :: FULLYVARIABLE=0, CONST_SUBSTEP=1, CONSTFULLSTEP=2
INTEGER(I4B)               :: JAC_RECOMPUTE
REAL(SP), ALLOCATABLE      :: fjacDCMP(:,:), fjacCOPY(:,:), fjacINDX(:) ! (temporary arrays)
! 7. Method used to trap/fix errors in Newton
INTEGER(I4B), PARAMETER    :: FULL_NEWTON=0, LINE_SEARCH=1
INTEGER(I4B)               :: CHECK_OVERSHOOT
! 8. Method used to process the small interval at the end of a time step
INTEGER(I4B), PARAMETER    :: STEP_TRUNC=0, LOOK_AHEAD=1, STEP_ABSORB=2
INTEGER(I4B)               :: SMALL_ENDSTEP
! ---------------------------------------------------------------------------------------
! (B) PARAMETERS
! ---------------------------------------------------------------------------------------
REAL(SP)                   :: ERR_TRUNC_ABS  ! Absolute temporal truncation error tolerance
REAL(SP)                   :: ERR_TRUNC_REL  ! Relative temporal truncation error tolerance
REAL(SP)                   :: ERR_ITER_FUNC  ! Iteration convergence tolerance for function values
REAL(SP)                   :: ERR_ITER_DX    ! Iteration convergence tolerance for dx
REAL(SP)                   :: FRACSTATE_MIN  ! Fractional minimum value of state (for non-zero derivatives)
REAL(SP)                   :: SAFETY         ! Safety factor in step-size equation
REAL(SP)                   :: RMIN           ! Minimum step size multiplier
REAL(SP)                   :: RMAX           ! Maximum step size multiplier
INTEGER(I4B)               :: NITER_TOTAL    ! Total number of iterations used in the implicit scheme
REAL(SP)                   :: MIN_TSTEP      ! Minimum time step length
REAL(SP)                   :: MAX_TSTEP      ! Maximum time step length
! ---------------------------------------------------------------------------------------
! (C) DIAGNOSTIX
! ---------------------------------------------------------------------------------------
INTEGER(I4B)               :: NUM_FUNCS      ! number of function calls
INTEGER(I4B)               :: NUM_JACOBIAN   ! number of times Jacobian is calculated
INTEGER(I4B)               :: NUMSUB_ACCEPT  ! number of sub-steps accepted (taken)
INTEGER(I4B)               :: NUMSUB_REJECT  ! number of sub-steps tried but rejected
INTEGER(I4B)               :: NUMSUB_NOCONV  ! number of sub-steps tried that did not converge
INTEGER(I4B)               :: MAXNUM_ITERNS  ! maximum number of iterations in the implicit scheme
INTEGER(I4B),DIMENSION(20) :: ORD_NSUBS = (/  1,  2,  5,  10,  20,  30,   50,   75,  100,   200, &
                                            300,500,750,1000,2000,5000,10000,20000,50000,100000/)
INTEGER(I4B),DIMENSION(20) :: PRB_NSUBS      ! cumulative probability for number of substeps taken
! ---------------------------------------------------------------------------------------
END MODULE MODEL_NUMERIX
