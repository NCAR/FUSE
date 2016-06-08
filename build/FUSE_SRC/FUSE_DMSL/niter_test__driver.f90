PROGRAM NITER_TEST__DRIVER
! ---------------------------------------------------------------------------------------
! Creator:
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! Driver program to assess the number of function evaluations taken in four different
!  configurations of the Newton scheme:
!   1) Newton's method with line searches (fixed Jacobian)
!   2) Newton's method without line searches (fixed Jacobian)
!   3) Newton's method with line searches (variable Jacobian)
!   4) Newton's method without line searches (variable Jacobian)
! ---------------------------------------------------------------------------------------
USE nrtype                                                ! variable types, etc.
USE ddirectory                                            ! directory for data files
! data modules
USE model_defn,nstateFUSE=>nstate                         ! model definition structures
USE multiforce, ONLY: AFORCE, DELTIM, NUMTIM              ! data interval = maximum model time step
USE multiparam, ONLY: LPARAM, PARATT, NUMPAR              ! parameter metadata structures
USE multiroute, ONLY: AROUTE                              ! model routing structures
USE multistats                                            ! model statistics structures
! informational modules
USE selectmodl_module                                     ! reads model control file
USE getpar_str_module                                     ! extracts parameter metadata
USE par_insert_module                                     ! inserts model parameters
! model numerix
USE model_numerix                                         ! defines decisions on model numerix
! access to model simulation modules
USE fuse_rmse_module                                      ! run model and compute the root mean squared error
IMPLICIT NONE
! ---------------------------------------------------------------------------------------
! (0) GET COMMAND-LINE ARGUMENTS...
! ---------------------------------------------------------------------------------------
CHARACTER(LEN=6)                       :: FMODEL_ID='      ' ! integer defining FUSE model
CHARACTER(LEN=6)                       :: NUMPARSET='      ' ! number of model parameter sets
INTEGER(I4B)                           :: FUSE_ID            ! integer definining FUSE model
INTEGER(I4B)                           :: NUMPSET            ! number of model parameter sets
! ---------------------------------------------------------------------------------------
! (1) SETUP MODELS FOR SIMULATION -- POPULATE DATA STRUCTURES
! ---------------------------------------------------------------------------------------
! get model forcing data
INTEGER(I4B)                           :: NTIM            ! number of time steps
INTEGER(I4B)                           :: INFERN_START    ! start of inference period
! get model setup
INTEGER(I4B)                           :: NMOD            ! number of models
INTEGER(I4B)                           :: ERR             ! error code
CHARACTER(LEN=256)                     :: MESSAGE         ! error message
! define model output
LOGICAL(LGT)                           :: OUTPUT_FLAG     ! .TRUE. = write time series output
INTEGER(I4B)                           :: ONEMOD=1        ! just specify one model
! ---------------------------------------------------------------------------------------
! (2) RUN MODEL FOR GIVEN PARAMETER SET AND DIFFERENT NUMERIX CONFIGURATIONS
! ---------------------------------------------------------------------------------------
INTEGER(I4B)                           :: IPAR    ! looping variable
TYPE(PARATT)                           :: PARAM_META ! parameter metadata (model parameters)
REAL(SP), DIMENSION(:), ALLOCATABLE    :: BL      ! vector of lower parameter bounds
REAL(SP), DIMENSION(:), ALLOCATABLE    :: BU      ! vector of upper parameter bounds
REAL(SP), DIMENSION(:), ALLOCATABLE    :: APAR    ! model parameter set
INTEGER(KIND=4)                        :: ISEED   ! seed for the random sequence
REAL(KIND=4),DIMENSION(:), ALLOCATABLE :: URAND   ! vector of quasi-random numbers U[0,1]
INTEGER(I4B)                           :: IPSET   ! (looping)
INTEGER(I4B)                           :: IJAC    ! (looping)
INTEGER(I4B)                           :: ISCH    ! (looping)
REAL(SP)                               :: RMSE    ! error from the simulation
! ---------------------------------------------------------------------------------------
! (0) READ COMMAND LINE ARGUMENTS
! ---------------------------------------------------------------------------------------
! read command-line arguments
CALL GETARG(1,FMODEL_ID)  ! integer defining FUSE model
CALL GETARG(2,NUMPARSET)  ! number of model parameter sets
! check command-line arguments
IF (LEN_TRIM(FMODEL_ID).EQ.0) STOP '1st command-line argument is missing (FMODEL_ID)'
IF (LEN_TRIM(NUMPARSET).EQ.0) STOP '2nd command-line argument is missing (NUMPARSET)'
! process command-line arguments
CALL GETNUMERIX()         ! defines method/parameters used for numerical solution
READ(FMODEL_ID,*) FUSE_ID ! integer definining FUSE model
READ(NUMPARSET,*) NUMPSET ! number of model parameter sets
! ---------------------------------------------------------------------------------------
! (1) GET MODEL SETUP -- MODEL DEFINITION, AND PARAMETER AND VARIABLE INFO FOR ALL MODELS
! ---------------------------------------------------------------------------------------
! Read data from the "BATEA-compliant" ASCII files
CALL GETFORCING(INFERN_START,NTIM)
! Define model attributes (valid for all models)
CALL UNIQUEMODL(NMOD)    ! get nmod unique models
CALL GETPARMETA()        ! read parameter metadata (parameter bounds etc.)
! Identify a single model
CALL SELECTMODL(FUSE_ID,ISTATUS=ERR,MESSAGE=MESSAGE)
IF (ERR.NE.0) THEN
 PRINT *, TRIM(MESSAGE); STOP
ENDIF
! Define list of states and parameters for the current model
CALL ASSIGN_STT()        ! state definitions are stored in module model_defn
CALL ASSIGN_FLX()        ! flux definitions are stored in module model_defn
CALL ASSIGN_PAR()        ! parameter defintions are stored in module multiparam
! Allocate space for the constant Jacobians
ALLOCATE(fjacCOPY(nstateFUSE,nstateFUSE),fjacDCMP(nstateFUSE,nstateFUSE),fjacINDX(nstateFUSE))
! compute derived model parameters (bucket sizes, etc.)
CALL PAR_DERIVE()
! Define output file names (shared in MODULE model_defn)
FNAME_NETCDF = TRIM(OUTPUT_PATH)//TRIM(SMODL%MNAME)//'_niter-test.nc'
! Define NetCDF output files (only write parameters and summary statistics)
ONEMOD=1                 ! one file per model (i.e., model dimension = 1)
PCOUNT=0                 ! counter for parameter sets evaluated (shared in MODULE multistats)
OUTPUT_FLAG = .TRUE.     ! .TRUE. if desire time series output
CALL DEF_PARAMS(ONEMOD)  ! define model parameters (initial CREATE)
CALL DEF_SSTATS()        ! define summary statistics (REDEF)
IF (OUTPUT_FLAG) CALL DEF_OUTPUT(NTIM)    ! define model time series (REDEF)
! ---------------------------------------------------------------------------------------
! (2) RUN MODEL FOR THE CURRENT PARAMETER SET WITH DIFFERENT NUMERIX OPTIONS
! ---------------------------------------------------------------------------------------
! get default numerix parameters
CALL DEFAULT_NUMERIX()
! get parameter bounds and random numbers
ALLOCATE(APAR(NUMPAR),BL(NUMPAR),BU(NUMPAR),URAND(NUMPAR))
DO IPAR=1,NUMPAR
 CALL GETPAR_STR(LPARAM(IPAR)%PARNAME,PARAM_META)
 BL(IPAR) = PARAM_META%PARLOW
 BU(IPAR) = PARAM_META%PARUPP
END DO
! loop through parameter sets
DO IPSET=963,963+NUMPSET-1
 ! get new parameter sets
 ISEED=IPSET; CALL I4_SOBOL(NUMPAR,ISEED,URAND)
 WRITE(*,'(I4,1X,12(E10.2,1X))') ISEED-1, URAND
 APAR = BL + URAND*(BU-BL)
 ! base run with fully-variable Jacobian
 JAC_RECOMPUTE = FULLYVARIABLE
 CALL FUSE_RMSE(APAR,RMSE,OUTPUT_FLAG)
 ! try freezing the Jacobian once we get "sufficiently close" to the solution
 JAC_RECOMPUTE = PERIOD_FREEZE
 DO IJAC=0,10,2
  THRESH_FRZE = 1. * 10.**-REAL(IJAC, KIND(SP))
  !CALL FUSE_RMSE(APAR,RMSE,OUTPUT_FLAG)
 END DO  ! (loop through different numerix parameter combinations)
 print *, '**********'
 ! try only re-computing Jacobian if don't get sufficiently large decrease
 ! in the norm of the residual vector
 JAC_RECOMPUTE = SMALL_F_RATIO
 DO IJAC=10,10,2
  THRESH_FRZE = REAL(IJAC, KIND(SP))/10._sp
  print *, THRESH_FRZE
  CALL FUSE_RMSE(APAR,RMSE,OUTPUT_FLAG)
 END DO  ! (loop through different numerix parameter combinations)
END DO  ! (loop through different parameter sets)
! and, deallocate space
DEALLOCATE(APAR,BL,BU,URAND)
STOP
END PROGRAM NITER_TEST__DRIVER
! --------------------------------------------------------------------------------------------------------------
SUBROUTINE DEFAULT_NUMERIX()
USE model_numerix
SOLUTION_METHOD        = IMPLICIT_EULER        ! implicit euler solution
TEMPORAL_ERROR_CONTROL = TS_FIXED              ! fixed time steps
INITIAL_NEWTON         = STATE_OLD             ! initial conditions for Newton
JAC_RECOMPUTE          = FULLYVARIABLE         ! fully variable Jacobian
CHECK_OVERSHOOT        = LINE_SEARCH           ! use line search to trap/fix overshoot problems
ERR_TRUNC_ABS          = 1.e-9                 ! absolute temporal truncation error tolerance
ERR_TRUNC_REL          = 1.e-9                 ! relative temporal truncation error tolerance
ERR_ITER_FUNC          = 1.e-9                 ! iteration convergence tolerance for function values
ERR_ITER_DX            = 1.e-9                 ! iteration convergence tolerance for dx
THRESH_FRZE            = 1.e-9                 ! Threshold for freezing the Jacobian
FRACSTATE_MIN          = 1.e-9                 ! fractional minimum value of state (for non-zero derivatives)
SAFETY                 = 0.9_sp                ! safety factor in step-size equation
RMIN                   = 0.1_sp                ! minimum step size multiplier
RMAX                   = 4.0_sp                ! maximum step size multiplier
NITER_TOTAL            = 100                   ! total number of iterations used in the implicit scheme
MIN_TSTEP              = 0.01_sp/60._sp/24._sp ! minimum time step length (minutes --> days)
MAX_TSTEP              = 1440.0_sp/60._sp/24._sp ! maximum time step length (minutes --> days)
END SUBROUTINE DEFAULT_NUMERIX
