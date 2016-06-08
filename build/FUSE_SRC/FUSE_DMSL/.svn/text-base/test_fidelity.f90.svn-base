PROGRAM TEST_FIDELITY
! ---------------------------------------------------------------------------------------
! Creator:
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! Driver program to test the fidelity of the different numerical methods
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
CHARACTER(LEN=6)                       :: FADAPTIVE='      ' ! identifier for adaptive sub-steps (0=fixed, 1=adaptive)
CHARACTER(LEN=6)                       :: TRUNC_ABS='      ' ! absolute temporal truncation error tolerance
CHARACTER(LEN=6)                       :: TRUNC_REL='      ' ! relative temporal truncation error tolerance
! ---------------------------------------------------------------------------------------
! (1) SETUP MODELS FOR SIMULATION -- POPULATE DATA STRUCTURES
! ---------------------------------------------------------------------------------------
! get model forcing data
INTEGER(I4B)                           :: NTIM            ! number of time steps
INTEGER(I4B)                           :: INFERN_START    ! start of inference period
! get model setup
INTEGER(I4B)                           :: FUSE_ID         ! integer defining FUSE model
INTEGER(I4B)                           :: I               ! looping
INTEGER(I4B)                           :: NMOD            ! number of models
INTEGER(I4B)                           :: ERR             ! error code
CHARACTER(LEN=256)                     :: MESSAGE         ! error message
! define model output
LOGICAL(LGT)                           :: OUTPUT_FLAG     ! .TRUE. = write time series output
INTEGER(I4B)                           :: ONEMOD          ! index of the model used (=1)
! ---------------------------------------------------------------------------------------
! (2) TEST FIDELITY
! ---------------------------------------------------------------------------------------
! Define error code for I/O
INTEGER(I4B)                           :: IERR            ! error code for I/O
! Identify index of the parameter set
INTEGER(I4B)                           :: IPARSET         ! parameter set index
CHARACTER(LEN=4)                       :: CPARSET         ! convert parameter set index to a string
! Parameter vectors
REAL(SP),DIMENSION(:),ALLOCATABLE      :: XDF             ! default parameter vector
REAL(SP)                               :: FPAR            ! function value for parameter set
! Loop through different time steps
INTEGER(I4B)                           :: IDEL            ! loop through different time steps
! ---------------------------------------------------------------------------------------
! (0) READ COMMAND LINE ARGUMENTS
! ---------------------------------------------------------------------------------------
! read command-line arguments
CALL GETARG(1,FMODEL_ID)  ! integer defining FUSE model
CALL GETARG(2,NSOLUTION)  ! numerical solution (0=explicit, 1=implicit)
CALL GETARG(3,FADAPTIVE)  ! identifier for adaptive sub-steps (0=fixed, 1=adaptive)
CALL GETARG(4,TRUNC_ABS)  ! absolute temporal truncation error tolerance
CALL GETARG(5,TRUNC_REL)  ! relative temporal truncation error tolerance
! check command-line arguments
IF (LEN_TRIM(FMODEL_ID).EQ.0) STOP '1st command-line argument is missing (FMODEL_ID)'
IF (LEN_TRIM(NSOLUTION).EQ.0) STOP '2nd command-line argument is missing (NSOLUTION)'
IF (LEN_TRIM(FADAPTIVE).EQ.0) STOP '3rd command-line argument is missing (FADAPTIVE)'
IF (LEN_TRIM(TRUNC_ABS).EQ.0) STOP '4th command-line argument is missing (TRUNC_ABS)'
IF (LEN_TRIM(TRUNC_REL).EQ.0) STOP '5th command-line argument is missing (TRUNC_REL)'
! convert command-line arguments to integer flags and real numbers
CALL GETNUMERIX()                         ! defines method/parameters used for numerical solution
READ(FMODEL_ID,*) FUSE_ID                 ! integer definining FUSE model
READ(NSOLUTION,*) SOLUTION_METHOD         ! numerical solution (0=EE, 1=EH, 2=IE, 3=IH)
READ(FADAPTIVE,*) TEMPORAL_ERROR_CONTROL  ! identifier for adaptive sub-steps (0=fixed, 1=adaptive)
READ(TRUNC_ABS,*) ERR_TRUNC_ABS           ! absolute temporal truncation error tolerance
READ(TRUNC_REL,*) ERR_TRUNC_REL           ! relative temporal truncation error tolerance
! additional checks
SELECT CASE(SOLUTION_METHOD); CASE(EXPLICIT_EULER,EXPLICIT_HEUN,IMPLICIT_EULER,IMPLICIT_HEUN,SEMI_IMPLICIT)
CASE DEFAULT;
 STOP 'solution method (1st command line argument) must equal 0 (explicit_euler), 1 (explicit heun), &
      &2 (implicit_euler), 3 (implicit_heun), or 4 (semi_implicit)'
END SELECT
SELECT CASE(TEMPORAL_ERROR_CONTROL); CASE(TS_FIXED,TS_ADAPT); CASE DEFAULT;
 STOP 'temporal error control (2nd command line argument) must equal 0 (fixed steps) or 1 (adaptive steps)'
END SELECT
write(*,'(A5,1X,2(I1,1X),2(E12.5,1X),I6,1X,A11,1X,I6)') 'FUSE ', &
SOLUTION_METHOD, TEMPORAL_ERROR_CONTROL, ERR_TRUNC_ABS, ERR_TRUNC_REL
! ---------------------------------------------------------------------------------------
! (1) GET MODEL SETUP -- MODEL DEFINITION, AND PARAMETER AND VARIABLE INFO FOR ALL MODELS
! ---------------------------------------------------------------------------------------
! Just assign data
INFERN_START=1; NTIM=1; NUMTIM=NTIM; DELTIM=1._SP
ALLOCATE(AFORCE(NTIM),AROUTE(NTIM))  ! (shared in module multiroute)
AFORCE(INFERN_START:NTIM)%PPT   = (/50./)
AFORCE(INFERN_START:NTIM)%PET   = (/ 5./)
! Define model attributes (valid for all models)
CALL UNIQUEMODL(NMOD)    ! get nmod unique models
CALL GETPARMETA()        ! read parameter metadata (parameter bounds etc.)
! Identify a single model (use command-line argument)
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
! ---------------------------------------------------------------------------------------
! (2) TEST METHOD FIDELITY
! ---------------------------------------------------------------------------------------
! allocate arrays
ALLOCATE(XDF(NUMPAR), STAT=IERR)
IF (IERR.NE.0) STOP ' problem allocating space for parameter arrays '
IPARSET = 0
! loop through example parameter sets
OPEN(21,FILE=TRIM(DATA_PATH)//'param_fidelity.dat')
DO
 ! read parameter set
 READ(21,*,IOSTAT=IERR) XDF; IF (IERR.NE.0) EXIT
 WRITE(*,'(20(A,1X))') LPARAM(1:NUMPAR); WRITE(*,'(20(F9.3,1X))') XDF
 ! increment counter
 IPARSET = IPARSET + 1
 ! convert counter to a character string
 CPARSET='    '; WRITE(CPARSET,'(I4)') IPARSET; CPARSET=ADJUSTR(CPARSET)
 FORALL(I=1:LEN(CPARSET)-LEN_TRIM(ADJUSTL(CPARSET))) CPARSET(I:I)='0'
 ! define NetCDF files (filename shared in MODULE model_defn)
 FNAME_NETCDF = TRIM(OUTPUT_PATH)//TRIM(SMODL%MNAME)//'__'//CPARSET//'__'//&
                TRIM(NSOLUTION)//'-'//TRIM(FADAPTIVE)//'__fidelity.nc'
 write(*,'(a)') trim(fname_netcdf)
 ONEMOD=1                 ! one file per model (i.e., model dimension = 1)
 PCOUNT=0                 ! counter for parameter sets in output file (shared in MODULE multistats)
 FCOUNT=0                 ! counter for parameter sets evaluated (shared in MODULE multistats)
 OUTPUT_FLAG = .TRUE.     ! write model output
 CALL DEF_PARAMS(ONEMOD)  ! define model parameters (initial CREATE)
 IF (OUTPUT_FLAG) CALL DEF_OUTPUT(NTIM)    ! define model time series (REDEF)
 CALL DEF_SSTATS()        ! define summary statistics (REDEF)
 DO IDEL=1,100
  DELTIM = REAL(IDEL,KIND(SP))/100._SP
  ! run model with example parameter sets
  CALL FUSE_RMSE(XDF,FPAR,OUTPUT_FLAG)
 END DO
END DO  ! looping through example parameter sets
DEALLOCATE(XDF, STAT=IERR)
IF (IERR.NE.0) STOP ' problem deallocating space for parameter arrays '
STOP
END PROGRAM TEST_FIDELITY
! --------------------------------------------------------------------------------------
