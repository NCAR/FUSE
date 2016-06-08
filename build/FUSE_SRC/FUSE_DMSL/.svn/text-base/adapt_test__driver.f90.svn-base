PROGRAM ADAPT_TEST__DRIVER
! ---------------------------------------------------------------------------------------
! Creator:
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! Driver program to evaluate the accuracy and efficiency of adaptive sub-stepping routines
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
CHARACTER(LEN=6)                       :: NSOLUTION='      ' ! numerical solution (0=implicit, 1=explicit)
CHARACTER(LEN=6)                       :: PAR_IDX  ='      ' ! index of parameter set
! ---------------------------------------------------------------------------------------
! (1) SETUP MODELS FOR SIMULATION -- POPULATE DATA STRUCTURES
! ---------------------------------------------------------------------------------------
! get model forcing data
INTEGER(I4B)                           :: NTIM            ! number of time steps
INTEGER(I4B)                           :: INFERN_START    ! start of inference period
! get model setup
INTEGER(I4B)                           :: FUSE_ID         ! integer defining FUSE model
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
INTEGER(I4B)                           :: ITRY    ! (looping)
INTEGER(I4B)                           :: JTRY    ! (looping)
REAL(SP)                               :: RMSE    ! error from the simulation
! ---------------------------------------------------------------------------------------
! (0) READ COMMAND LINE ARGUMENTS
! ---------------------------------------------------------------------------------------
! read command-line arguments
CALL GETARG(1,FMODEL_ID)  ! integer defining FUSE model
CALL GETARG(2,NSOLUTION)  ! numerical solution
CALL GETARG(3,PAR_IDX)    ! index in the Sobol sequence
! check command-line arguments
IF (LEN_TRIM(NSOLUTION).EQ.0) STOP '1st command-line argument is missing (FMODEL_ID)'
IF (LEN_TRIM(NSOLUTION).EQ.0) STOP '2nd command-line argument is missing (NSOLUTION)'
IF (LEN_TRIM(PAR_IDX)  .EQ.0) STOP '3rd command-line argument is missing (PAR_IDX)'
! convert command-line arguments to integer flags and real numbers
CALL GETNUMERIX()                         ! defines method/parameters used for numerical solution
READ(FMODEL_ID,*) FUSE_ID                 ! integer definining FUSE model
READ(NSOLUTION,*) SOLUTION_METHOD         ! numerical solution (0=EE, 1=EH, 2=IE, 3=IH, 4=SI)
READ(PAR_IDX,*)   ISEED                   ! convert index to an integer
! check solution method
SELECT CASE(SOLUTION_METHOD); CASE(EXPLICIT_EULER,EXPLICIT_HEUN,IMPLICIT_EULER,IMPLICIT_HEUN,SEMI_IMPLICIT)
CASE DEFAULT;
 PRINT *, 'solution method (2nd command line argument) must equal 0 (explicit_euler), 1 (explicit heun),'//&
          ' 2 (implicit_euler), 3 (implicit_heun), or 4 (semi_implicit)'
 STOP
END SELECT
! ---------------------------------------------------------------------------------------
! (1) GET MODEL SETUP -- MODEL DEFINITION, AND PARAMETER AND VARIABLE INFO FOR ALL MODELS
! ---------------------------------------------------------------------------------------
! Read data from the "BATEA-compliant" ASCII files
CALL GETFORCING(INFERN_START,NTIM)
! Define model attributes (valid for all models)
CALL UNIQUEMODL(NMOD)    ! get nmod unique models
CALL GETPARMETA()        ! read parameter metadata (parameter bounds etc.)
! Identify a single model (read control file ../DataFiles/m_decisions.txt)
CALL SELECTMODL(FUSE_ID,ISTATUS=ERR,MESSAGE=MESSAGE)
IF (ERR.NE.0) THEN
 PRINT *, TRIM(MESSAGE); STOP
ENDIF
! Define list of states and parameters for the current model
CALL ASSIGN_STT()        ! state definitions are stored in module model_defn
CALL ASSIGN_FLX()        ! flux definitions are stored in module model_defn
CALL ASSIGN_PAR()        ! parameter defintions are stored in module multiparam
! compute derived model parameters (bucket sizes, etc.)
CALL PAR_DERIVE()
! Define output file names (shared in MODULE model_defn)
FNAME_NETCDF = TRIM(OUTPUT_PATH)//TRIM(SMODL%MNAME)//'_adapt-steps_'//&
                TRIM(PAR_IDX)//'_'//TRIM(NSOLUTION)//'.nc'
! Define NetCDF output files (only write parameters and summary statistics)
ONEMOD=1                 ! one file per model (i.e., model dimension = 1)
PCOUNT=0                 ! counter for parameter sets evaluated (shared in MODULE multistats)
OUTPUT_FLAG = .FALSE.    ! .TRUE. if desire time series output
CALL DEF_PARAMS(ONEMOD)  ! define model parameters (initial CREATE)
CALL DEF_SSTATS()        ! define summary statistics (REDEF)
IF (OUTPUT_FLAG) CALL DEF_OUTPUT(NTIM)    ! define model time series (REDEF)
! ---------------------------------------------------------------------------------------
! (2) RUN MODEL FOR THE CURRENT PARAMETER SET WITH DIFFERENT NUMERIX OPTIONS
! ---------------------------------------------------------------------------------------
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
!DO IPAR=1,NUMPAR; WRITE(*,'(A11,1X,F9.3)') LPARAM(IPAR)%PARNAME, APAR(IPAR); END DO
! create the exact solution
TEMPORAL_ERROR_CONTROL = TS_ADAPT              ! adaptive time steps
ERR_TRUNC_ABS          = 1.e-9                 ! absolute temporal truncation error tolerance
ERR_TRUNC_REL          = 1.e-9                 ! relative temporal truncation error tolerance
MIN_TSTEP              = 0.01_sp/60._sp/24._sp ! minimum time step length (minutes --> days)
MAX_TSTEP              = 10.0_sp/60._sp/24._sp ! maximum time step length (minutes --> days)
! run model (parameters and statistics are written in FUSE_RMSE)
CALL FUSE_RMSE(APAR,RMSE,OUTPUT_FLAG)
! save solution for subsequent testing
AROUTE(:)%Q_ACCURATE = AROUTE(:)%Q_ROUTED
! modify numerix parameters
MAX_TSTEP              = DELTIM                ! max step length = data interval
! evaluate different parameters for step-size control
DO ITRY=3,9,3    ! play with different ERR_TRUNC_ABS parameters
 ERR_TRUNC_ABS = 1. * 10.**-REAL(ITRY, KIND(SP))
 DO JTRY=1,9      ! play with different ERR_TRUNC_REL parameters
  ERR_TRUNC_REL = 1. * 10.**-REAL(JTRY, KIND(SP))
  ! run zee model
  write(*,'(2(E15.7,1X))') ERR_TRUNC_ABS, ERR_TRUNC_REL
  CALL FUSE_RMSE(APAR,RMSE,OUTPUT_FLAG)
 END DO  ! (loop through different numerix parameter combinations)
END DO  ! (loop through different numerix parameter combinations)
! for reference, include the fixed-step method
TEMPORAL_ERROR_CONTROL = TS_FIXED              ! fixed time steps
CALL FUSE_RMSE(APAR,RMSE,OUTPUT_FLAG)          ! run zee model
! and, deallocate space
DEALLOCATE(APAR,BL,BU,URAND)
STOP
END PROGRAM ADAPT_TEST__DRIVER
