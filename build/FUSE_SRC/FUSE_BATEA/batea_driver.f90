PROGRAM BATEA_DRIVER
! ---------------------------------------------------------------------------------------
! Creator:
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! Driver program for BATEA simulations
! ---------------------------------------------------------------------------------------
USE nrtype                                        ! variable types, etc.
USE ddirectory                                    ! directory for data files
! data modules
USE model_defn,nstateFUSE=>nstate                 ! model definition structures
USE multiparam, ONLY: PARATT, NUMPAR              ! parameter metadata structures
USE multiforce                                    ! model forcing structure
USE multiroute                                    ! model routing structures
USE multistats                                    ! model statistics structures
! informational modules
USE selectmodl_module                             ! reads model control file
USE getpar_str_module                             ! extracts parameter metadata
USE par_insert_module                             ! inserts model parameters
IMPLICIT NONE
! get forcing data
CHARACTER(LEN=8)                       :: CBASID  ! basin id
INTEGER(I4B)                           :: NTIM    ! number of time steps
INTEGER(I4B)                           :: WARM_START   ! index of start of warm-up period
INTEGER(I4B)                           :: INFERN_START ! index of start of inference period
INTEGER(I4B)                           :: INFERN_END   ! index of end of inference period
! get model setup
INTEGER(I4B)                           :: I,J,K   ! looping
INTEGER(I4B)                           :: NMOD    ! number of models
INTEGER(I4B)                           :: ERR     ! error code
CHARACTER(LEN=256)                     :: MESSAGE ! error message
! get command line arguments
CHARACTER(LEN=256)                     :: BATEA_PATH   ! path of BATEA file
CHARACTER(LEN=256)                     :: BATEA_NAME   ! name of BATEA file
CHARACTER(LEN=256)                     :: FN_EXTENSION ! filename extension
CHARACTER(LEN=256)                     :: STN_NUMBER   ! station number
INTEGER(I4B)                           :: ISTN         ! station index
! define output files
CHARACTER(LEN=256)                     :: FNAME_ASCII  ! ascii output file name
INTEGER(I4B)                           :: ONEMOD  ! index for defining output file (one file per model)
LOGICAL(LGT)                           :: OUTPUT_FLAG  ! .TRUE. if desire time series output
LOGICAL(LGT)                           :: SSTATS_FLAG  ! .TRUE. if desire summary statistics
! model parameters
INTEGER(I4B), PARAMETER                :: IUNIT=21 ! file unit
INTEGER(I4B)                           :: IHEAD   ! loop through header lines
INTEGER(I4B)                           :: IO_ERR  ! error code for I/O
CHARACTER(LEN=256),DIMENSION(4)        :: TMPTXT  ! temporary text
CHARACTER(LEN=256)                     :: PARINFO ! parameter information string
INTEGER(I4B)                           :: IPOS    ! position of sub-string in string (identify multipliers)
INTEGER(I4B)                           :: JPOS    ! position of sub-string in string (identify multipliers)
INTEGER(I4B)                           :: INDX    ! index of rainfall multipliers
INTEGER(I4B), PARAMETER                :: NM=100  ! number of rainfall multipliers
INTEGER(I4B)                           :: N_MULT  ! actual number of rainfall multipliers
INTEGER(I4B), DIMENSION(NM)            :: IM_FILE ! index of rainfall multipliers in file
INTEGER(I4B), DIMENSION(NM)            :: IM_INF  ! index of rainfall multipliers in inference period
INTEGER(I4B)                           :: N_PARS  ! number of parameter values
INTEGER(I4B)                           :: M_PARS  ! number of deterministic parameter values
INTEGER(I4B), PARAMETER                :: PMAX=20 ! maximum number of deterministic parameters
CHARACTER(LEN=256),DIMENSION(PMAX)     :: PARNAME ! names of model parameters
INTEGER(I4B),DIMENSION(PMAX)           :: IP_FILE ! index of deterministic parameters in file
! loop through different model parameters
INTEGER(I4B), PARAMETER                :: N_SETS=5001 ! number of parameter sets
REAL(SP)                               :: APAR    ! parameter value
REAL(SP), DIMENSION(:), ALLOCATABLE    :: PAR_VEC    ! parameter vector
INTEGER(I4B)                           :: IERR    ! error code for allocate statement
INTEGER(I4B)                           :: IMULT   ! looping variable
INTEGER(I4B)                           :: IPARS   ! looping variable
INTEGER(I4B)                           :: JPARS   ! looping variable
INTEGER(I4B)                           :: KPARS   ! looping variable
INTEGER(I4B)                           :: IPARSET ! looping variable
INTEGER(I4B)                           :: ISTART  ! start index used to compute summary statistics
! ---------------------------------------------------------------------------------------
! (1) GET COMMAND LINE ARGUMENTS
! ---------------------------------------------------------------------------------------
! get path for BATEA file
CALL GETARG(1,BATEA_PATH)
IF (LEN_TRIM(BATEA_PATH).EQ.0) STOP ' need path of BATEA file as 1st command-line argument '
! get name of BATEA file
CALL GETARG(2,BATEA_NAME)
IF (LEN_TRIM(BATEA_NAME).EQ.0) STOP ' need name of BATEA file as 2nd command-line argument '
! get filename extension
CALL GETARG(3,FN_EXTENSION)
IF (LEN_TRIM(FN_EXTENSION).EQ.0) STOP ' need filename extension as 3rd command-line argument '
! get station number
CALL GETARG(4,STN_NUMBER)
IF (LEN_TRIM(STN_NUMBER).EQ.0) STOP ' need station number as 4th command-line argument '
CALL SYSTEM('cp ../input/forcinginfo_default.txt ../input/forcinginfo.txt')
CALL SYSTEM('sed -i "s/mpc_column_index/'//TRIM(STN_NUMBER)//'/" ../input/forcinginfo.txt')
! ---------------------------------------------------------------------------------------
! (2) GET MODEL SETUP -- MODEL DEFINITION, AND PARAMETER AND VARIABLE INFO FOR ALL MODELS
! ---------------------------------------------------------------------------------------
! Read data from the "BATEA-compliant" ASCII files (and make a copy of it)
CALL GETFORCING(INFERN_START,NTIM)
ALLOCATE(CFORCE(NTIM),STAT=IERR); IF (IERR.NE.0) STOP ' problem allocating CFORCE '
CFORCE=AFORCE
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
CALL ASSIGN_PAR()        ! parameter defintions are stored in module multiparam
! compute derived model parameters (bucket sizes, etc.)
CALL PAR_DERIVE()
! --------------------------------------------------------------------------------------
! (3) DEFINE NETCDF OUTPUT FILES 
! --------------------------------------------------------------------------------------
! Define output file names (shared in MODULE model_defn)
FNAME_NETCDF = TRIM(OUTPUT_PATH)//TRIM(SMODL%MNAME)//TRIM(FN_EXTENSION)//'.nc'
! Define NetCDF output files (only write parameters and summary statistics)
ONEMOD=1                 ! one file per model (i.e., model dimension = 1)
PCOUNT=0                 ! counter for parameter sets evaluated (shared in MODULE multistats)
CALL DEF_PARAMS(ONEMOD)  ! define model parameters (initial CREATE)
OUTPUT_FLAG = .TRUE.     ! .TRUE. if desire model output
SSTATS_FLAG = .TRUE.     ! .TRUE. if desire summary statistics
IF (OUTPUT_FLAG) CALL DEF_OUTPUT(NTIM)    ! define model output (REDEF)
IF (SSTATS_FLAG) CALL DEF_SSTATS()        ! define summary statistics (REDEF)
! --------------------------------------------------------------------------------------
! (4) GET MODEL PARAMETERS
! --------------------------------------------------------------------------------------
! open parameter file
N_PARS=0  ! initialize number of model parameters
N_MULT=0  ! initialize number of rainfall multipliers
M_PARS=0  ! initialize number of deterministic model parameters
print *, TRIM(BATEA_NAME)
IF (INDEX(TRIM(BATEA_NAME),'sls').GT.0) THEN
 OPEN(IUNIT,FILE=TRIM(BATEA_PATH)//TRIM(BATEA_NAME)//'.MODAL',STATUS='old')
  ! read header
  DO IHEAD=1,3; READ(IUNIT,*,IOSTAT=IO_ERR) TMPTXT(1); END DO  ! read header
  DO
   ! read data
   READ(IUNIT,*,IOSTAT=IO_ERR) TMPTXT
   IF (IO_ERR.NE.0) EXIT
   ! insert multipliers
   READ(TMPTXT(1),*) APAR

   print *, trim(tmptxt(1)), trim(tmptxt(3))
  END DO
 CLOSE(IUNIT)
ENDIF
pause
OPEN(IUNIT,FILE=TRIM(BATEA_PATH)//TRIM(BATEA_NAME)//'.parfile',STATUS='old')
 ! read header
 DO IHEAD=1,3; READ(IUNIT,*,IOSTAT=IO_ERR) TMPTXT(1); END DO  ! read header
 DO
  ! read data
  READ(IUNIT,*,IOSTAT=IO_ERR) TMPTXT; PARINFO=TMPTXT(4)
  IF (IO_ERR.NE.0 .OR. INDEX(TMPTXT(1),'-----').GT.0) EXIT
  N_PARS=N_PARS+1
  ! identify the index of rainfall multipliers in file and in inference period
  IPOS = INDEX(PARINFO,'RFERR_MLT'); JPOS=INDEX(PARINFO,'lat')
  IF (IPOS.GT.0 .AND. JPOS.GT.0) THEN
   N_MULT=N_MULT+1; IF (N_MULT.GT.NM) STOP ' insufficient number of rainfall multipliers '
   READ(PARINFO(INDEX(PARINFO,'=')+1:INDEX(PARINFO,']')-1),*) INDX
   IM_FILE(N_MULT)=N_PARS
   IM_INF(N_MULT) =INDX
  ENDIF
  ! identify the index of deterministic multipliers in file (and retrieve their names)
  IF (IPOS.LE.0) THEN
   M_PARS=M_PARS+1
   PARNAME(M_PARS) = PARINFO
   IP_FILE(M_PARS) = N_PARS
  ENDIF
 END DO
CLOSE(IUNIT)
! --------------------------------------------------------------------------------------
! (5) EXTRACT PARAMETERS AND PRIME MODEL
! --------------------------------------------------------------------------------------
KPARS=0
! open file with parameter multipliers
OPEN(IUNIT,FILE=TRIM(BATEA_PATH)//TRIM(BATEA_NAME)//'[prodkt].sam',STATUS='old')
 ALLOCATE(PAR_VEC(N_PARS), STAT=IERR); IF (IERR.NE.0) STOP ' problem allocating parameter vectors '
 ! loop through parameter sets
 DO IPARS=1,N_SETS
  ! get parameters
  READ(IUNIT,*) PAR_VEC
  ! only process every 10th sample
  KPARS=KPARS+1; IF (KPARS.EQ.10) KPARS=0
  IF (KPARS.NE.0) CYCLE
  ! refresh forcing data (structures corrupted with precipitation multipliers)
  AFORCE=CFORCE
  ! assign parameter multipliers
  DO IMULT=1,N_MULT
   AFORCE(INFERN_START+IM_INF(IMULT)-1)%PPT = AFORCE(INFERN_START+IM_INF(IMULT)-1)%PPT*PAR_VEC(IM_FILE(IMULT))
   !print *, IM_INF(IMULT), IM_FILE(IMULT), PAR_VEC(IM_FILE(IMULT))
   !print *, INFERN_START+IM_INF(IMULT)-1
  END DO
  print *, '********** PARAMETER SET ', IPARS, ' **********'
  WRITE(*,'(25(F4.2,1X))') PAR_VEC(IM_FILE(:))
  ! assign model parameters
  DO JPARS=1,M_PARS
   CALL PAR_INSERT(PAR_VEC(IP_FILE(JPARS)),TRIM(PARNAME(JPARS)))
   WRITE(*,'(A11,1X,F12.5)') TRIM(PARNAME(JPARS)), PAR_VEC(IP_FILE(JPARS))
  END DO
  ! compute derived model parameters (bucket sizes, etc.)
  CALL PAR_DERIVE()
  ! --------------------------------------------------------------------------------------
  ! (5) RUN MODEL
  ! --------------------------------------------------------------------------------------
  CALL BMODEL_RUN(OUTPUT_FLAG,SSTATS_FLAG)
  ! write summary statistics
  IF (SSTATS_FLAG) THEN
   ISTART=INFERN_START       ! model warm-up period
   CALL MEAN_STATS()         ! compute summary statistics
   CALL PUT_SSTATS(PCOUNT,1) ! 1 = just one model for numerix test
  ENDIF
 END DO  ! (looping through parameter sets)
CLOSE(IUNIT)
STOP
END PROGRAM BATEA_DRIVER
