PROGRAM URS_DRIVER
! ---------------------------------------------------------------------------------------
! Creator:
! Martyn Clark, 2011
! ---------------------------------------------------------------------------------------
! Purpose:
! Driver program to perform multiple runs of a model by uniform random sampling from the
!  feasible parameter space.
! ---------------------------------------------------------------------------------------
USE nrtype                                                ! variable types, etc.
USE fuse_fileManager,only:fuse_SetDirsUndPhiles,&         ! sets directories and filenames
     OUTPUT_PATH,FORCINGINFO
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
CHARACTER(LEN=80)                      :: MUSTRFILE='                                                                                ' ! path/name of muster file
CHARACTER(LEN=12)                      :: MBASIN_ID='            ' ! MOPEX basin ID
CHARACTER(LEN=6)                       :: FMODEL_ID='      ' ! integer defining FUSE model
CHARACTER(LEN=6)                       :: NSOLUTION='      ' ! numerical solution (0=explicit Euler; 1=explicit Heun; 2=implicit Euler; 3=implicit Heun, 4=semi-implicit)
CHARACTER(LEN=6)                       :: FADAPTIVE='      ' ! identifier for adaptive sub-steps (0=fixed, 1=adaptive)
CHARACTER(LEN=6)                       :: TRUNC_ABS='      ' ! absolute temporal truncation error tolerance
CHARACTER(LEN=6)                       :: TRUNC_REL='      ' ! relative temporal truncation error tolerance
CHARACTER(LEN=12)                      :: TSTEP_LEN='            ' ! maximum length of the time step (days)
CHARACTER(LEN=6)                       :: NUMPARSET='      ' ! number of parameter sets
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
! (2) RUN MODEL FOR DIFFERENT PARAMETER SETS
! ---------------------------------------------------------------------------------------
INTEGER(I4B)                           :: IPAR    ! loop thru model parameters
INTEGER(I4B)                           :: IPSET   ! loop thru model parameter sets
INTEGER(I4B)                           :: NUMPSET ! number of parameter sets
TYPE(PARATT)                           :: PARAM_META ! parameter metadata (model parameters)
REAL(SP), DIMENSION(:), ALLOCATABLE    :: BL      ! vector of lower parameter bounds
REAL(SP), DIMENSION(:), ALLOCATABLE    :: BU      ! vector of upper parameter bounds
REAL(SP), DIMENSION(:), ALLOCATABLE    :: APAR    ! model parameter set
INTEGER(KIND=4)                        :: ISEED   ! seed for the random sequence
REAL(KIND=4),DIMENSION(:), ALLOCATABLE :: URAND   ! vector of quasi-random numbers U[0,1]
REAL(SP)                               :: RMSE    ! error from the simulation
! ---------------------------------------------------------------------------------------
! (0) READ COMMAND LINE ARGUMENTS
! ---------------------------------------------------------------------------------------
! read command-line arguments
CALL GETARG(1,MUSTRFILE)  ! path/name of muster file
CALL GETARG(2,MBASIN_ID)  ! MOPEX basin ID
CALL GETARG(3,FMODEL_ID)  ! integer defining FUSE model
CALL GETARG(4,NSOLUTION)  ! numerical solution (0=explicit, 1=implicit)
CALL GETARG(5,FADAPTIVE)  ! identifier for adaptive sub-steps (0=fixed, 1=adaptive)
CALL GETARG(6,TRUNC_ABS)  ! absolute temporal truncation error tolerance
CALL GETARG(7,TRUNC_REL)  ! relative temporal truncation error tolerance
CALL GETARG(8,TSTEP_LEN)  ! maximum length of the time step (days)
CALL GETARG(9,NUMPARSET)  ! number of parameter sets
! check command-line arguments
IF (LEN_TRIM(MUSTRFILE).EQ.0) STOP '1st command-line argument is missing (MUSTRFILE)'
IF (LEN_TRIM(MBASIN_ID).EQ.0) STOP '2nd command-line argument is missing (MBASIN_ID)'
IF (LEN_TRIM(FMODEL_ID).EQ.0) STOP '3rd command-line argument is missing (FMODEL_ID)'
IF (LEN_TRIM(NSOLUTION).EQ.0) STOP '4th command-line argument is missing (NSOLUTION)'
IF (LEN_TRIM(FADAPTIVE).EQ.0) STOP '5th command-line argument is missing (FADAPTIVE)'
IF (LEN_TRIM(TRUNC_ABS).EQ.0) STOP '6th command-line argument is missing (TRUNC_ABS)'
IF (LEN_TRIM(TRUNC_REL).EQ.0) STOP '7th command-line argument is missing (TRUNC_REL)'
IF (LEN_TRIM(TSTEP_LEN).EQ.0) STOP '8th command-line argument is missing (TSTEP_LEN)'
IF (LEN_TRIM(NUMPARSET).EQ.0) STOP '9th command-line argument is missing (NUMPARSET)'
! get directories and filenames for control files
call fuse_SetDirsUndPhiles(trim(MUSTRFILE),err=err,message=message)
if (err.ne.0) write(*,*) trim(message); if (err.gt.0) stop
! define basin desired
FORCINGINFO = 'forcinginfo.'//TRIM(MBASIN_ID)//'.txt'
! convert command-line arguments to integer flags and real numbers
CALL GETNUMERIX(err,message)              ! defines method/parameters used for numerical solution
if (err.ne.0) write(*,*) trim(message); if (err.gt.0) stop
READ(FMODEL_ID,*) FUSE_ID                 ! integer definining FUSE model
READ(NSOLUTION,*) SOLUTION_METHOD         ! numerical solution (0=implicit, 1=explicit)
READ(FADAPTIVE,*) TEMPORAL_ERROR_CONTROL  ! identifier for adaptive sub-steps (0=fixed, 1=adaptive)
READ(TRUNC_ABS,*) ERR_TRUNC_ABS           ! absolute temporal truncation error tolerance
READ(TRUNC_REL,*) ERR_TRUNC_REL           ! relative temporal truncation error tolerance
READ(TSTEP_LEN,*) MAX_TSTEP               ! maximum length of the time step (days)
READ(NUMPARSET,*) NUMPSET                 ! number of parameter sets
! additional checks
SELECT CASE(SOLUTION_METHOD); CASE(EXPLICIT_EULER,EXPLICIT_HEUN,IMPLICIT_EULER,IMPLICIT_HEUN,SEMI_IMPLICIT)
CASE DEFAULT
 PRINT *, 'solution method (1st command line argument) must equal 0 (explicit_euler), 1 (explicit heun), '//&
          '2 (implicit_euler), 3 (implicit_heun), or 4 (semi_implicit)'
 STOP
END SELECT
SELECT CASE(TEMPORAL_ERROR_CONTROL); CASE(TS_FIXED,TS_ADAPT); CASE DEFAULT;
 STOP 'temporal error control (2nd command line argument) must equal 0 (fixed steps) or 1 (adaptive steps)'
END SELECT
write(*,'(A5,1X,2(I1,1X),2(E12.5,1X),I6,1X,A11,1X,I6)') 'FUSE ', &
SOLUTION_METHOD, TEMPORAL_ERROR_CONTROL, ERR_TRUNC_ABS, ERR_TRUNC_REL
! ---------------------------------------------------------------------------------------
! (1) GET MODEL SETUP -- MODEL DEFINITION, AND PARAMETER AND VARIABLE INFO FOR ALL MODELS
! ---------------------------------------------------------------------------------------
! Read data from the "BATEA-compliant" ASCII files
CALL GETFORCING(INFERN_START,NTIM)
! Define model attributes (valid for all models)
CALL UNIQUEMODL(NMOD,err,message)   ! get nmod unique models
IF (ERR.NE.0) WRITE(*,*) TRIM(MESSAGE); IF (ERR.GT.0) STOP
CALL GETPARMETA(ERR,MESSAGE)        ! read parameter metadata (parameter bounds etc.)
IF (ERR.NE.0) WRITE(*,*) TRIM(MESSAGE); IF (ERR.GT.0) STOP
! Identify a single model (read control file ../DataFiles/m_decisions.txt)
CALL SELECTMODL(FUSE_ID,ISTATUS=ERR,MESSAGE=MESSAGE)
IF (ERR.NE.0) WRITE(*,*) TRIM(MESSAGE); IF (ERR.GT.0) STOP
! Define list of states and parameters for the current model
CALL ASSIGN_STT()        ! state definitions are stored in module model_defn
CALL ASSIGN_FLX()        ! flux definitions are stored in module model_defn
CALL ASSIGN_PAR()        ! parameter defintions are stored in module multiparam
! compute derived model parameters (bucket sizes, etc.)
CALL PAR_DERIVE(ERR,MESSAGE)
IF (ERR.NE.0) WRITE(*,*) TRIM(MESSAGE); IF (ERR.GT.0) STOP
! Define output file names (shared in MODULE model_defn)
FNAME_NETCDF = TRIM(OUTPUT_PATH)//'URS_MOPEX__'//TRIM(MBASIN_ID)//'__'//TRIM(SMODL%MNAME)//'__'//&
               TRIM(NSOLUTION)//'-'//TRIM(FADAPTIVE)//'__'//&
               TRIM(TRUNC_ABS)//'-'//TRIM(TRUNC_REL)//'__'//&
               TRIM(TSTEP_LEN)//'.nc'
write(*,'(a)') trim(fname_netcdf)
! Define NetCDF output files (only write parameters and summary statistics)
ONEMOD=1                 ! one file per model (i.e., model dimension = 1)
PCOUNT=0                 ! counter for parameter sets evaluated (shared in MODULE multistats)
OUTPUT_FLAG = .TRUE.    ! .TRUE. if desire time series output
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
! loop through parameter sets
DO IPSET=1,NUMPSET
 ! get new parameter sets
 ISEED=IPSET; CALL I4_SOBOL(NUMPAR,ISEED,URAND)
 WRITE(*,'(I4,1X,12(E10.2,1X))') ISEED-1, URAND
 APAR = BL + URAND*(BU-BL)
 ! run zee model
 CALL FUSE_RMSE(APAR,RMSE,OUTPUT_FLAG)
END DO
! and, deallocate space
DEALLOCATE(APAR,BL,BU,URAND)
STOP
END PROGRAM URS_DRIVER
