PROGRAM SOBOL_DRIVER
! ---------------------------------------------------------------------------------------
! Creator:
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! Driver program to loop through example parameter sets
! ---------------------------------------------------------------------------------------
USE nrtype                                                ! variable types, etc.
USE ddirectory                                            ! directory for data files
! data modules
USE model_defn,nstateFUSE=>nstate                         ! model definition structures
USE multiforce, ONLY: DELTIM                              ! data interval = maximum model time step
USE multiparam, ONLY: LPARAM, PARATT, NUMPAR, SOBOL_INDX  ! parameter metadata structures
USE multiroute                                            ! model routing structures
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
CHARACTER(LEN=12)                      :: MBASIN_ID='      ' ! MOPEX basin ID
CHARACTER(LEN=6)                       :: FMODEL_ID='      ' ! integer defining FUSE model
CHARACTER(LEN=6)                       :: NSOLUTION='      ' ! numerical solution (0=implicit, 1=explicit)
CHARACTER(LEN=6)                       :: FADAPTIVE='      ' ! identifier for adaptive sub-steps (0=fixed, 1=adaptive)
CHARACTER(LEN=6)                       :: TRUNC_ABS='      ' ! absolute temporal truncation error tolerance
CHARACTER(LEN=6)                       :: TRUNC_REL='      ' ! relative temporal truncation error tolerance
CHARACTER(LEN=256)                     :: PARAMFILE          ! filename with list of parameter sets 
! ---------------------------------------------------------------------------------------
! (1) SETUP MODELS FOR SIMULATION -- POPULATE DATA STRUCTURES
! ---------------------------------------------------------------------------------------
! get model forcing data
INTEGER(I4B)                           :: NTIM            ! number of time steps
INTEGER(I4B)                           :: INFERN_START    ! start of inference period
! get model setup
INTEGER(I4B)                           :: FUSE_ID         ! integer defining FUSE model
INTEGER(I4B)                           :: I,J,K           ! looping
INTEGER(I4B)                           :: NMOD            ! number of models
INTEGER(I4B)                           :: ERR             ! error code
CHARACTER(LEN=256)                     :: MESSAGE         ! error message
! define model output
LOGICAL(LGT)                           :: OUTPUT_FLAG     ! .TRUE. = write time series output
! ---------------------------------------------------------------------------------------
! (2) CREATE PARAMETER GRID
! ---------------------------------------------------------------------------------------
! Identify existence of the parameter file
LOGICAL(LGT)                           :: LEXIST          ! .TRUE. if the parameter file exists
INTEGER(I4B), PARAMETER                :: IN_UNIT=21      ! file unit for parameter fie
! Define error code for I/O
INTEGER(I4B)                           :: IERR            ! error code for I/O
! Define parameters
CHARACTER(LEN=512)                     :: CHEAD           ! header text
INTEGER(I4B)                           :: NUM_ALLPAR      ! number of all possible parameters
TYPE PAR_TXT
 CHARACTER(LEN=11)                     :: PARAM_NAME      ! parameter name
ENDTYPE PAR_TXT
TYPE(PAR_TXT),DIMENSION(:),ALLOCATABLE :: PARNAMES_ALL    ! list of all possible parameter names 
INTEGER(I4B)                           :: IPOS,JPOS,KPOS  ! position in header string
INTEGER(I4B)                           :: IPAR_ALL        ! loop through all possible model parameters
! Index and values of parameters
REAL(SP),DIMENSION(:),ALLOCATABLE      :: ALLPARS         ! vector of model all parameters
REAL(SP),DIMENSION(:),ALLOCATABLE      :: TRYPARS         ! vector of model parameters to trial
INTEGER(I4B)                           :: IPAR_MOD        ! loop through parameters of the current model
INTEGER(I4B)                           :: ONEMOD          ! index of the model used (=1)
REAL(SP)                               :: FPAR            ! function value for parameter set
! ---------------------------------------------------------------------------------------
! (0) READ COMMAND LINE ARGUMENTS
! ---------------------------------------------------------------------------------------
! read command-line arguments
CALL GETARG(1,MBASIN_ID)  ! MOPEX basin ID 
CALL GETARG(2,FMODEL_ID)  ! integer defining FUSE model
CALL GETARG(3,NSOLUTION)  ! numerical solution (0=explicit, 1=implicit)
CALL GETARG(4,FADAPTIVE)  ! identifier for adaptive sub-steps (0=fixed, 1=adaptive)
CALL GETARG(5,TRUNC_ABS)  ! absolute temporal truncation error tolerance
CALL GETARG(6,TRUNC_REL)  ! relative temporal truncation error tolerance
CALL GETARG(7,PARAMFILE)  ! filename of the parameter sets
! check command-line arguments
IF (LEN_TRIM(MBASIN_ID).EQ.0) STOP '1st command-line argument is missing (MBASIN_ID)'
IF (LEN_TRIM(FMODEL_ID).EQ.0) STOP '2nd command-line argument is missing (FMODEL_ID)'
IF (LEN_TRIM(NSOLUTION).EQ.0) STOP '3rd command-line argument is missing (NSOLUTION)'
IF (LEN_TRIM(FADAPTIVE).EQ.0) STOP '4th command-line argument is missing (FADAPTIVE)'
IF (LEN_TRIM(TRUNC_ABS).EQ.0) STOP '5th command-line argument is missing (TRUNC_ABS)'
IF (LEN_TRIM(TRUNC_REL).EQ.0) STOP '6th command-line argument is missing (TRUNC_REL)'
IF (LEN_TRIM(PARAMFILE).EQ.0) STOP '7th command-line argument is missing (PARAMFILE)'
! define basin desired
FORCINGINFO = 'forcinginfo.'//TRIM(MBASIN_ID)//'.txt'
! convert command-line arguments to integer flags and real numbers
CALL GETNUMERIX()                         ! defines method/parameters used for numerical solution
READ(FMODEL_ID,*) FUSE_ID                 ! integer definining FUSE model
READ(NSOLUTION,*) SOLUTION_METHOD         ! numerical solution (0=implicit, 1=explicit)
READ(FADAPTIVE,*) TEMPORAL_ERROR_CONTROL  ! identifier for adaptive sub-steps (0=fixed, 1=adaptive)
READ(TRUNC_ABS,*) ERR_TRUNC_ABS           ! absolute temporal truncation error tolerance
READ(TRUNC_REL,*) ERR_TRUNC_REL           ! relative temporal truncation error tolerance
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
CALL UNIQUEMODL(NMOD)    ! get nmod unique models
CALL GETPARMETA()        ! read parameter metadata (parameter bounds etc.)
! Identify a single model (read control file ../DataFiles/m_decisions.txt)
CALL SELECTMODL(FUSE_ID,ISTATUS=ERR,MESSAGE=MESSAGE)
IF (ERR.NE.0) THEN
 PRINT *, TRIM(MESSAGE); STOP
ENDIF
! Define list of states and parameters for the current model
CALL ASSIGN_STT()        ! state definitions are stored in module model_defn
CALL ASSIGN_PAR()        ! parameter defintions are stored in module multiparam
CALL ASSIGN_FLX()        ! flux definitions stored in module model_defn
! compute derived model parameters (bucket sizes, etc.)
CALL PAR_DERIVE()
! define NetCDF files (filename shared in MODULE model_defn)
FNAME_NETCDF = TRIM(OUTPUT_PATH)//'DMSL_'//TRIM(MBASIN_ID)//'__'//TRIM(SMODL%MNAME)//'__'//&
               TRIM(NSOLUTION)//'-'//TRIM(FADAPTIVE)//'__'//TRIM(PARAMFILE)//'.nc'
write(*,'(a)') trim(fname_netcdf)
ONEMOD=1                 ! one file per model (i.e., model dimension = 1)
PCOUNT=0                 ! counter for parameter sets in output file (shared in MODULE multistats)
FCOUNT=0                 ! counter for parameter sets evaluated (shared in MODULE multistats)
OUTPUT_FLAG = .FALSE.    ! write model output
CALL DEF_PARAMS(ONEMOD)  ! define model parameters (initial CREATE)
IF (OUTPUT_FLAG) CALL DEF_OUTPUT(NTIM)    ! define model time series (REDEF)
CALL DEF_SSTATS()        ! define summary statistics (REDEF)
! ---------------------------------------------------------------------------------------
! (2) LOOP THROUGH EXAMPLE PARAMETER SETS
! ---------------------------------------------------------------------------------------
! check that the file exists
write(*,'(a)') TRIM(DATA_PATH)//TRIM(PARAMFILE)//'.dat'
INQUIRE(FILE=TRIM(DATA_PATH)//TRIM(PARAMFILE)//'.dat',EXIST=LEXIST)
IF (.NOT.LEXIST) STOP ' parameter file does not exist '
! open file
OPEN(IN_UNIT,FILE=TRIM(DATA_PATH)//TRIM(PARAMFILE)//'.dat',STATUS='old')
 NUM_ALLPAR=0
 ! read header
 DO
  READ(IN_UNIT,'(A512)') CHEAD      ! read header line
  IF (CHEAD(1:7).EQ.'ParFlag') EXIT ! title line identified by 'ParFlag'
  NUM_ALLPAR=NUM_ALLPAR+1           ! increment number total parameters
 END DO
 ! strip out the parameter names
 ALLOCATE(PARNAMES_ALL(NUM_ALLPAR)); IPOS=8
 IPAR_ALL=0
 DO
  ! get param index
  IPAR_ALL=IPAR_ALL+1
  IF (IPAR_ALL.GT.NUM_ALLPAR) EXIT
  ! extract a "word"
  JPOS=INDEX(CHEAD(IPOS:LEN_TRIM(CHEAD)),' ')
  KPOS=INDEX(CHEAD(JPOS:LEN_TRIM(CHEAD)),' ')
  ! add the parameter name to the structure (and fill with white space)
  PARNAMES_ALL(IPAR_ALL)%PARAM_NAME(JPOS:KPOS+1) = CHEAD(IPOS+JPOS:IPOS+JPOS+KPOS)
  IF (KPOS+1.LT.LEN(PARNAMES_ALL(IPAR_ALL)%PARAM_NAME)) &
   FORALL(I=KPOS+2:LEN(PARNAMES_ALL(IPAR_ALL)%PARAM_NAME)) PARNAMES_ALL(IPAR_ALL)%PARAM_NAME(I:I)=' '
  ! move to the next word
  IPOS=IPOS+JPOS+KPOS
  DO
   IF (CHEAD(IPOS+1:IPOS+1).NE.' ') EXIT
   IPOS=IPOS+1
  END DO
  ! check exit criteria
  IF (IPOS.GT.LEN_TRIM(CHEAD)) EXIT
 END DO
 ! allocate vector for the parameters
 ALLOCATE(ALLPARS(NUM_ALLPAR),TRYPARS(NUMPAR))
 ! loop through parameters
 DO
  ! read a line of parameters (SOBOL_INDX is stored in module multiparam)
  READ(IN_UNIT,*,IOSTAT=IERR) SOBOL_INDX, ALLPARS  
  IF (IERR.NE.0) EXIT
  ! extract the parameters that we need
  DO IPAR_MOD=1,NUMPAR
   DO IPAR_ALL=1,NUM_ALLPAR 
    IF (TRIM(LPARAM(IPAR_MOD)%PARNAME).EQ.TRIM(PARNAMES_ALL(IPAR_ALL)%PARAM_NAME)) THEN
     TRYPARS(IPAR_MOD) = ALLPARS(IPAR_ALL)
     !WRITE(*,'(A11,1X,F9.3,1X)') TRIM(LPARAM(IPAR_MOD)%PARNAME), TRYPARS(IPAR_MOD)
    ENDIF
   END DO
  END DO
  ! run model (parameters and statistics are written in FUSE_RMSE)
  CALL FUSE_RMSE(TRYPARS,FPAR,OUTPUT_FLAG)
 END DO
CLOSE(IN_UNIT)
STOP
END PROGRAM SOBOL_DRIVER
! --------------------------------------------------------------------------------------
