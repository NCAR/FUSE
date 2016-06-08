PROGRAM NMX_DRIVER
! ---------------------------------------------------------------------------------------
! Creator:
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! Driver program for model numerix tests
! ---------------------------------------------------------------------------------------
USE nrtype                                        ! variable types, etc.
USE ddirectory                                    ! directory for data files
! data modules
USE model_defn,nstateFUSE=>nstate                 ! model definition structures
USE multiforce, ONLY: DELTIM                      ! data interval = maximum model time step
USE multiparam, ONLY: LPARAM, PARATT, NUMPAR      ! parameter metadata structures
USE multiroute                                    ! model routing structures
USE multistats                                    ! model statistics structures
! informational modules
USE selectmodl_module                             ! reads model control file
USE getpar_str_module                             ! extracts parameter metadata
USE par_insert_module                             ! inserts model parameters
! model numerix
USE model_numerix
IMPLICIT NONE
! get forcing data
CHARACTER(LEN=8)                       :: CBASID  ! basin id
INTEGER(I4B)                           :: NTIM    ! number of time steps
INTEGER(I4B)                           :: INFERN_START  ! start of inference period
! get model setup
INTEGER(I4B)                           :: I,J,K   ! looping
INTEGER(I4B)                           :: NMOD    ! number of models
INTEGER(I4B)                           :: ERR     ! error code
CHARACTER(LEN=256)                     :: MESSAGE ! error message
! define output files
INTEGER(I4B)                           :: ONEMOD  ! index for defining output file (one file per model)
LOGICAL(LGT)                           :: OUTPUT_FLAG  ! .TRUE. if desire time series output
LOGICAL(LGT)                           :: SSTATS_FLAG  ! .TRUE. if desire summary statistics
! get command-line arguments
CHARACTER(LEN=11)                      :: NUM_EXPERIMENT  ! name of numerical experiment
CHARACTER(LEN=11)                      :: PARNAM  ! parameter name
CHARACTER(LEN=11)                      :: CRANGE  ! range for parameter cut
CHARACTER(LEN=11)                      :: PAR_IDX ! index of parameter set
! loop through different model parameters
INTEGER(I4B)                           :: IPAR    ! looping variable
INTEGER(I4B)                           :: JPAR    ! looping variable
INTEGER(I4B)                           :: IPARSET ! looping variable
INTEGER(I4B)                           :: NCUT    ! number of parameter values in the "cut"
TYPE(PARATT)                           :: PARAM_META ! parameter metadata (model parameters)
REAL(SP)                               :: XDEF    ! default parameter value
REAL(SP)                               :: XLOW    ! lower parameter bound
REAL(SP)                               :: XUPP    ! upper parameter bound
REAL(SP), DIMENSION(:), ALLOCATABLE    :: BL      ! vector of lower parameter bounds
REAL(SP), DIMENSION(:), ALLOCATABLE    :: BU      ! vector of upper parameter bounds
REAL(SP), DIMENSION(:), ALLOCATABLE    :: APAR    ! model parameter set
INTEGER(KIND=4)                        :: ISEED   ! seed for the random sequence
REAL(KIND=4),DIMENSION(:), ALLOCATABLE :: URAND   ! vector of quasi-random numbers U[0,1]
REAL(SP)                               :: XFRC    ! fractional range for the parameter cut
REAL(SP)                               :: XRNG    ! range for the parameter cut
REAL(SP)                               :: XINC    ! parameter increment
REAL(SP)                               :: XPAR    ! parameter value
! loop through different parameter sets
INTEGER(I4B)                           :: ITRY    ! (looping)
INTEGER(I4B)                           :: JTRY    ! (looping)
INTEGER(I4B)                           :: KTRY    ! (looping)
INTEGER(I4B)                           :: MTRY    ! (looping)
INTEGER(I4B)                           :: NTRY    ! (looping)
! ---------------------------------------------------------------------------------------
! (0) RETRIEVE COMMAND-LINE ARGUMENTS
! ---------------------------------------------------------------------------------------
! get name of numerical experiment
CALL GETARG(1,NUM_EXPERIMENT)
IF (LEN_TRIM(NUM_EXPERIMENT).EQ.0) &
 STOP ' need name of numerical experiment as 1st command-line argument '
! ---------------------------------------------------------------------------------------
! get parameters to diagnose smoothing
IF (TRIM(NUM_EXPERIMENT).EQ.'DIAG_SMOOTH') THEN
 ! get parameter name
 CALL GETARG(2,PARNAM)
 IF (LEN_TRIM(PARNAM).EQ.0) STOP ' need parameter name as 2nd command-line argument '
 ! get range for cut
 CALL GETARG(3,CRANGE)
 IF (LEN_TRIM(CRANGE).EQ.0) STOP ' need range for cut as 3rd command-line argument '
 READ(CRANGE,*) XFRC   ! convert range to to a real number
ENDIF
! ---------------------------------------------------------------------------------------
! get index of parameter set in the sobol sequence
IF (TRIM(NUM_EXPERIMENT).EQ.'ADAPT_STEPS') THEN
 CALL GETARG(2,PAR_IDX)
 IF (LEN_TRIM(PAR_IDX).EQ.0) STOP ' need index for parameter set as 2nd command-line argument '
 READ(PAR_IDX,*) ISEED   ! convert index to an integer
ENDIF
! ---------------------------------------------------------------------------------------
! (1) GET MODEL FORCING DATA AND STORE IN MEMORY
! ---------------------------------------------------------------------------------------
! Read data from the "BATEA-compliant" ASCII files
CALL GETFORCING(INFERN_START,NTIM)
! ---------------------------------------------------------------------------------------
! (2) GET MODEL SETUP -- MODEL DEFINITION, AND PARAMETER AND VARIABLE INFO FOR ALL MODELS
! ---------------------------------------------------------------------------------------
! Define model attributes (valid for all models)
CALL UNIQUEMODL(NMOD)    ! get nmod unique models
CALL GETNUMERIX()        ! defines method/parameters used for numerical solution
CALL GETPARMETA()        ! read parameter metadata (parameter bounds etc.)
! Identify a single model (read control file ../DataFiles/m_decisions.txt)
CALL SELECTMODL(ERR,MESSAGE)
IF (ERR.NE.0) THEN
 PRINT *, TRIM(MESSAGE); STOP
ENDIF
! Define list of states and parameters for the current model
CALL ASSIGN_STT()        ! state definitions are stored in module model_defn
CALL ASSIGN_FLX()        ! flux definitions are stored in module model_defn
CALL ASSIGN_PAR()        ! parameter defintions are stored in module multiparam
! compute derived model parameters (bucket sizes, etc.)
CALL PAR_DERIVE()
! --------------------------------------------------------------------------------------
! (3) DEFINE NETCDF OUTPUT FILES 
! --------------------------------------------------------------------------------------
! Define output file names (shared in MODULE model_defn)
SELECT CASE(TRIM(NUM_EXPERIMENT))
 CASE('DIAG_SMOOTH')
  FNAME_NETCDF = TRIM(OUTPUT_PATH)//TRIM(SMODL%MNAME)//'__'//TRIM(PARNAM)//'__'//TRIM(CRANGE)//'.nc'
 CASE('EVAL_JACOBN')
  FNAME_NETCDF = TRIM(OUTPUT_PATH)//TRIM(SMODL%MNAME)//'_eval_jacobn.nc'
 CASE('CONV_PARAMS')
  FNAME_NETCDF = TRIM(OUTPUT_PATH)//TRIM(SMODL%MNAME)//'_conv_params.nc'
 CASE('LIMIT_ITERS')
  FNAME_NETCDF = TRIM(OUTPUT_PATH)//TRIM(SMODL%MNAME)//'_limit_iters.nc'
 CASE('FIXED_STEPS')
  FNAME_NETCDF = TRIM(OUTPUT_PATH)//TRIM(SMODL%MNAME)//'_fixed-steps.nc'
 CASE('ADAPT_STEPS')
  FNAME_NETCDF = TRIM(OUTPUT_PATH)//TRIM(SMODL%MNAME)//'_adapt-steps_'//TRIM(PAR_IDX)//'.nc'
 CASE DEFAULT
  STOP ' 1st command line argument must be DIAG_SMOOTH, EVAL_JACOBN, CONV_PARAMS, LIMIT_ITERS, FIXED_STEPS, or ADAPT_STEPS '
END SELECT
! Define NetCDF output files (only write parameters and summary statistics)
ONEMOD=1                 ! one file per model (i.e., model dimension = 1)
PCOUNT=0                 ! counter for parameter sets evaluated (shared in MODULE multistats)
CALL DEF_PARAMS(ONEMOD)  ! define model parameters (initial CREATE)
SELECT CASE(TRIM(NUM_EXPERIMENT))
 CASE('DIAG_SMOOTH'); OUTPUT_FLAG = .FALSE.  ! .TRUE. if desire time series output
 CASE('EVAL_JACOBN'); OUTPUT_FLAG = .TRUE.   ! .TRUE. if desire time series output
 CASE('CONV_PARAMS'); OUTPUT_FLAG = .TRUE.   ! .TRUE. if desire time series output
 CASE('LIMIT_ITERS'); OUTPUT_FLAG = .TRUE.   ! .TRUE. if desire time series output
 CASE('FIXED_STEPS'); OUTPUT_FLAG = .TRUE.   ! .TRUE. if desire time series output
 CASE('ADAPT_STEPS'); OUTPUT_FLAG = .FALSE.  ! .TRUE. if desire time series output
END SELECT
SSTATS_FLAG = .TRUE.     ! .TRUE. if desire summary statistics
IF (OUTPUT_FLAG) CALL DEF_OUTPUT(NTIM)    ! define model output (REDEF)
IF (SSTATS_FLAG) CALL DEF_SSTATS()        ! define summary statistics (REDEF)
! --------------------------------------------------------------------------------------
! (4) TRY DIFFERENT NUMERICAL METHODS/CONSTANTS
! --------------------------------------------------------------------------------------
SELECT CASE(TRIM(NUM_EXPERIMENT))
 CASE('DIAG_SMOOTH')
  ! get parameter bounds and the parameter default values
  CALL GETPAR_STR(TRIM(PARNAM),PARAM_META) 
  XLOW = PARAM_META%PARLOW
  XUPP = PARAM_META%PARUPP
  XDEF = PARAM_META%PARDEF
  ! re-set parameters
  CALL DEFAULT_NUMERIX()            ! get default numerix parameters
  NCUT      = 100                   ! number of parameter sets in the "cut"                      
  MAX_TSTEP =   1.                  ! max step length = 1 day
  TEMPORAL_ERROR_CONTROL = TS_ADAPT ! adaptive time steps 
  ! loop through different numerical methods
  DO SOLUTION_METHOD=1,1
   DO TRUNCATION_ERROR=1,1
    DO ORDER_ACCEPT=1,1
     ! evaluate different parameters for step-size control
     DO ITRY=0,4  ! play with different ERR_TRUNC_ABS parameters
      ERR_TRUNC_ABS = 1. * 10.**-REAL(ITRY, KIND(SP))
      DO JTRY=0,4    ! play with different ERR_TRUNC_REL parameters
       ERR_TRUNC_REL = 1. * 10.**-REAL(JTRY, KIND(SP))
       ! evaluate different error parameters for connvergence of implicit solution
       DO MTRY=0,9,2    !  play with different ERR_ITER_FUNC parameters
        ERR_ITER_FUNC = 1. * 10.**-REAL(MTRY, KIND(SP))
        DO NTRY=0,9,2    ! play with different ERR_ITER_DX parameters
         ERR_ITER_DX = 1. * 10.**-REAL(NTRY, KIND(SP))
         ! get NCUT increments
         XRNG = (XUPP-XLOW)*XFRC
         XINC = XRNG/REAL(NCUT,KIND(SP))
         DO IPAR=0,NCUT
          ! modify parameter value
          XPAR = (XDEF - XRNG/2._SP) + REAL(IPAR,KIND(SP))*XINC
          IF (XPAR.LT.XLOW) XPAR=XLOW
          IF (XPAR.GT.XUPP) XPAR=XUPP
          CALL PAR_INSERT(XPAR,PARNAM)
          CALL NMODEL_RUN(OUTPUT_FLAG,SSTATS_FLAG)
          ! write summary statistics
          IF (SSTATS_FLAG) THEN
           CALL MEAN_STATS()               ! compute summary statistics
           CALL PUT_SSTATS(PCOUNT,1)       ! 1 = just one model for numerix test
          ENDIF
         END DO  ! ipar
         WRITE(*,'(2(A11,1X),7(I2,1X))') PARNAM, CRANGE, &
          solution_method, truncation_error, order_accept, itry, jtry, mtry, ntry
        END DO  ! ntry
       END DO  ! mtry
      END DO  ! jtry
     END DO  ! itry
    END DO  ! order_accept
   END DO  ! truncation_error
  END DO  ! solution_method
 ! -------------------------------------------------------------------------------------
 CASE('EVAL_JACOBN')
  ! assess different Jacobian re-evaluation strategies
  CALL DEFAULT_NUMERIX()                   ! get default numerix parameters
  MAX_TSTEP =   1.                         ! max step length = 1 day
  SOLUTION_METHOD        = IMPLICIT_EULER  ! implicit euler solution
  TEMPORAL_ERROR_CONTROL = TS_ADAPT        ! adaptive time steps
  ! loop through different numerical methods
  DO TRUNCATION_ERROR=0,1
   DO ORDER_ACCEPT=0,1
    ! evaluate different parameters for step-size control
    DO ITRY=0,4  ! play with different ERR_TRUNC_ABS parameters
     ERR_TRUNC_ABS = 1. * 10.**-REAL(ITRY, KIND(SP))
     DO JTRY=0,4    ! play with different ERR_TRUNC_REL parameters
      ERR_TRUNC_REL = 1. * 10.**-REAL(JTRY, KIND(SP))
      DO KTRY=0,2
       JAC_RECOMPUTE=KTRY
       IF (JAC_RECOMPUTE.EQ.CONSTFULLSTEP) &
        ALLOCATE(fjacCOPY(nstateFUSE,nstateFUSE),fjacDCMP(nstateFUSE,nstateFUSE),&
                 fjacINDX(nstateFUSE) )
       ! evaluate different error parameters for connvergence of implicit solution
       DO MTRY=0,9,2    !  play with different ERR_ITER_FUNC parameters
        ERR_ITER_FUNC = 1. * 10.**-REAL(MTRY, KIND(SP))
        DO NTRY=0,9,2    ! play with different ERR_ITER_DX parameters
         ERR_ITER_DX = 1. * 10.**-REAL(NTRY, KIND(SP))
         write(*,'(7(I2,1X))') TRUNCATION_ERROR, ORDER_ACCEPT, ITRY, JTRY, KTRY, MTRY, NTRY
         CALL NMODEL_RUN(OUTPUT_FLAG,SSTATS_FLAG)
        END DO  ! ntry
       END DO  ! mtry
       IF (JAC_RECOMPUTE.EQ.CONSTFULLSTEP) DEALLOCATE(fjacCOPY,fjacDCMP,fjacINDX)
      END DO  ! ktry
     END DO  ! jtry
    END DO  ! itry
   END DO  ! order_accept
  END DO  ! truncation_error
 ! -------------------------------------------------------------------------------------
 CASE('CONV_PARAMS')
  ! evaluate impact of convergence parameters in the implicit scheme
  CALL DEFAULT_NUMERIX()                   ! get default numerix parameters
  ! modify numerix parameters
  MAX_TSTEP              = 1.              ! max step length = 1 day
  SOLUTION_METHOD        = IMPLICIT_EULER  ! implicit euler solution
  TEMPORAL_ERROR_CONTROL = TS_ADAPT        ! adaptive time steps
  DO TRUNCATION_ERROR=0,1
   DO ORDER_ACCEPT=0,1
    ! evaluate different parameters for step-size control
    DO ITRY=0,9,2    ! play with different ERR_TRUNC_ABS parameters
     ERR_TRUNC_ABS = 1. * 10.**-REAL(ITRY, KIND(SP))
     DO JTRY=0,9,2    ! play with different ERR_TRUNC_REL parameters
      ERR_TRUNC_REL = 1. * 10.**-REAL(JTRY, KIND(SP))
      ! evaluate different error parameters for connvergence of implicit solution
      DO MTRY=0,9,2    !  play with different ERR_ITER_FUNC parameters
       ERR_ITER_FUNC = 1. * 10.**-REAL(MTRY, KIND(SP))
       DO NTRY=0,9,2    ! play with different ERR_ITER_DX parameters
        ERR_ITER_DX = 1. * 10.**-REAL(NTRY, KIND(SP))
        ! run zee model
        CALL NMODEL_RUN(OUTPUT_FLAG,SSTATS_FLAG)
        print *, TRUNCATION_ERROR, ORDER_ACCEPT, ITRY, JTRY, MTRY, NTRY
       END DO  ! ntry
      END DO  ! mtry
     END DO  ! jtry
    END DO  ! itry
   END DO  ! order_accept
  END DO  ! truncation_error
 ! -------------------------------------------------------------------------------------
 CASE('LIMIT_ITERS')
  ! limit the number of iterations in the implicit scheme
  CALL DEFAULT_NUMERIX()            ! get default numerix parameters
  ! modify numerix parameters
  MAX_TSTEP              = 1.       ! max step length = 1 day
  TEMPORAL_ERROR_CONTROL = TS_ADAPT ! adaptive time steps
  DO TRUNCATION_ERROR=0,1
   DO ORDER_ACCEPT=0,1
    ! evaluate different parameters for step-size control
    DO ITRY=0,9,2    ! play with different ERR_TRUNC_ABS parameters
     ERR_TRUNC_ABS = 1. * 10.**-REAL(ITRY, KIND(SP))
     DO JTRY=0,9,2    ! play with different ERR_TRUNC_REL parameters
      ERR_TRUNC_REL = 1. * 10.**-REAL(JTRY, KIND(SP))
      DO KTRY=0,1
       ! modify minimum step-size multiplier
       IF (KTRY.EQ.0) RMIN = 0.1_sp
       IF (KTRY.EQ.1) RMIN = 0.5_sp
       ! loop through different number of iterations
       DO NITER_TOTAL=1,10      
        ! run zee model
        CALL NMODEL_RUN(OUTPUT_FLAG,SSTATS_FLAG)
       END DO
      END DO  ! ktry 
     END DO  ! jtry
    END DO  ! itry
   END DO
  END DO
 ! -------------------------------------------------------------------------------------
 CASE('FIXED_STEPS')
  ! fixed time steps, different solution methods and error control
  SSTATS_FLAG=.FALSE.                      ! don't compute statistics
  !CALL DEFAULT_NUMERIX()                   ! get default numerix parameters
  !CALL NMODEL_RUN(OUTPUT_FLAG,SSTATS_FLAG) ! base run with default parameters
  ! save solution for subsequent testing
  AROUTE(:)%Q_ACCURATE = AROUTE(:)%Q_ROUTED
  ! modify numerix parameters
  SSTATS_FLAG            = .TRUE.          ! compute summary statistics 
  MAX_TSTEP              = 1.              ! max step length = 1 day
  TEMPORAL_ERROR_CONTROL = TS_FIXED        ! fixed time steps 
  ! loop through different numerical methods
  DO SOLUTION_METHOD=0,1
   DO TRUNCATION_ERROR=0,1
    DO ORDER_ACCEPT=0,1
     ! run zee model
     CALL NMODEL_RUN(OUTPUT_FLAG,SSTATS_FLAG)
     ! write summary statistics
     IF (SSTATS_FLAG) THEN
      CALL MEAN_STATS()               ! compute summary statistics
      CALL PUT_SSTATS(PCOUNT,1)       ! 1 = just one model for numerix test
     ENDIF
    END DO
   END DO
  END DO
 ! -------------------------------------------------------------------------------------
 CASE('ADAPT_STEPS')  ! adaptive time steps for multiple parameter sets
  ! get parameter bounds and random numbers
  ALLOCATE(APAR(NUMPAR),BL(NUMPAR),BU(NUMPAR),URAND(NUMPAR))
  DO IPAR=1,NUMPAR
   CALL GETPAR_STR(LPARAM(IPAR)%PARNAME,PARAM_META) 
   BL(IPAR) = PARAM_META%PARLOW
   BU(IPAR) = PARAM_META%PARUPP
  END DO
  ! get new parameter sets
  CALL I4_SOBOL(NUMPAR,ISEED,URAND)
  WRITE(*,'(I4,1X,12(E10.2,1X))') ISEED-1, URAND
  APAR = BL + URAND*(BU-BL)
  CALL PUT_PARSET(APAR) 
  ! create the exact solution
  SSTATS_FLAG=.FALSE.                      ! don't compute statistics
  CALL DEFAULT_NUMERIX()                   ! get default numerix parameters
  CALL NMODEL_RUN(OUTPUT_FLAG,SSTATS_FLAG) ! base run with default parameters
  ! save solution for subsequent testing
  AROUTE(:)%Q_ACCURATE = AROUTE(:)%Q_ROUTED
  ! modify numerix parameters
  SSTATS_FLAG            = .TRUE.          ! compute summary statistics
  TEMPORAL_ERROR_CONTROL = TS_ADAPT        ! adaptive time steps 
  MAX_TSTEP              = DELTIM          ! max step length = data interval
  ! loop through different numerical methods
  DO SOLUTION_METHOD=0,1
   ! evaluate different parameters for step-size control
   DO ITRY=3,9,3    ! play with different ERR_TRUNC_ABS parameters
    ERR_TRUNC_ABS = 1. * 10.**-REAL(ITRY, KIND(SP))
    DO JTRY=1,9      ! play with different ERR_TRUNC_REL parameters
     ERR_TRUNC_REL = 1. * 10.**-REAL(JTRY, KIND(SP))
     ! run zee model
     CALL NMODEL_RUN(OUTPUT_FLAG,SSTATS_FLAG)
     WRITE(*,'(I4,1X,F9.4,1X,5(I1,1X))') &
      ISEED-1, DELTIM, SOLUTION_METHOD, TRUNCATION_ERROR, ITRY, JTRY
     ! compute and write summary statistics
     IF (SSTATS_FLAG) THEN
      CALL MEAN_STATS()               ! compute summary statistics
      CALL PUT_SSTATS(PCOUNT,1)       ! 1 = just one model for numerix test
     ENDIF
    END DO  ! (loop through different numerix parameter combinations)
   END DO  ! (loop through different numerix parameter combinations)

  END DO
  ! for reference, include the fixed-step implicit euler method
  SOLUTION_METHOD        = IMPLICIT_EULER        ! implicit euler solution
  TEMPORAL_ERROR_CONTROL = TS_FIXED              ! fixed time steps
  ORDER_ACCEPT           = LOWER_ORDER           ! accept lower-order solutions
  CALL NMODEL_RUN(OUTPUT_FLAG,SSTATS_FLAG)
  ! compute and write summary statistics
  IF (SSTATS_FLAG) THEN
   CALL MEAN_STATS()               ! compute summary statistics
   CALL PUT_SSTATS(PCOUNT,1)       ! 1 = just one model for numerix test
  ENDIF
 ! -------------------------------------------------------------------------------------
 CASE DEFAULT
 STOP ' 1st command line argument must be DIAG_SMOOTH, EVAL_JACOBN, CONV_PARAMS, LIMIT_ITERS, FIXED_STEPS, or ADAPT_STEPS '
 ! -------------------------------------------------------------------------------------
END SELECT
STOP
END PROGRAM NMX_DRIVER
! --------------------------------------------------------------------------------------
SUBROUTINE DEFAULT_NUMERIX()
USE model_numerix
SOLUTION_METHOD        = IMPLICIT_EULER        ! implicit euler solution
TEMPORAL_ERROR_CONTROL = TS_ADAPT              ! adaptive time steps
INITIAL_NEWTON         = EXPLICIT_FULL         ! initial conditions for Newton
JAC_RECOMPUTE          = FULLYVARIABLE         ! fully variable Jacobian
CHECK_OVERSHOOT        = LINE_SEARCH           ! use line search to trap/fix overshoot problems
ERR_TRUNC_ABS          = 1.e-9                 ! absolute temporal truncation error tolerance
ERR_TRUNC_REL          = 1.e-9                 ! relative temporal truncation error tolerance
ERR_ITER_FUNC          = 1.e-9                 ! iteration convergence tolerance for function values
ERR_ITER_DX            = 1.e-9                 ! iteration convergence tolerance for dx
FRACSTATE_MIN          = 1.e-9                 ! fractional minimum value of state (for non-zero derivatives)
SAFETY                 = 0.9_sp                ! safety factor in step-size equation
RMIN                   = 0.1_sp                ! minimum step size multiplier
RMAX                   = 4.0_sp                ! maximum step size multiplier
NITER_TOTAL            = 100                   ! total number of iterations used in the implicit scheme
MIN_TSTEP              = 0.01_sp/60._sp/24._sp ! minimum time step length (minutes --> days)
MAX_TSTEP              = 60.0_sp/60._sp/24._sp ! maximum time step length (minutes --> days)
END SUBROUTINE DEFAULT_NUMERIX
