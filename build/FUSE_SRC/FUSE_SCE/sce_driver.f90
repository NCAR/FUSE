PROGRAM sce_driver
! ---------------------------------------------------------------------------------------
! Creator:
! Martyn Clark, 2008
! ---------------------------------------------------------------------------------------
! Purpose:
! Driver program for SCE
! ---------------------------------------------------------------------------------------
USE nrtype                                        ! variable types, etc.
!USE ddirectory                                   ! directory for data files - commented because couldn't be found
USE fuse_fileManager,only:fuse_SetDirsUndPhiles,& ! use fuse_fileManager instead - sets directories and filenames
     OUTPUT_PATH,FORCINGINFO
! data modules
USE model_defn                                    ! model definition structures
USE multiparam, ONLY: PARATT, LPARAM, NUMPAR      ! parameter metadata structures
USE multistats                                    ! model statistics structures
USE model_numerix                                 ! model numerix structures
! informational modules
USE selectmodl_module                             ! reads model control file
USE getpar_str_module                             ! extracts parameter metadata
IMPLICIT NONE
! command-line arguments
CHARACTER(LEN=6)                       :: FMODEL_ID='      ' ! integer defining FUSE model
CHARACTER(LEN=6)                       :: NSOLUTION='      ' ! numerical solution (0=implicit, 1=explicit)
CHARACTER(LEN=6)                       :: FADAPTIVE='      ' ! identifier for adaptive sub-steps (0=fixed, 1=adaptive)
CHARACTER(LEN=6)                       :: TRUNC_ABS='      ' ! absolute temporal truncation error tolerance
CHARACTER(LEN=6)                       :: TRUNC_REL='      ' ! relative temporal truncation error tolerance
! forcing data
INTEGER(I4B)                           :: INFERN_START    ! start of inference period
INTEGER(I4B)                           :: NTIM    ! number of time steps
! model setup
INTEGER(I4B)                           :: FUSE_ID ! integer definining FUSE model
INTEGER(I4B)                           :: I,J,K   ! looping
INTEGER(I4B)                           :: NMOD    ! number of models
INTEGER(I4B)                           :: ERR     ! error code
CHARACTER(LEN=256)                     :: MESSAGE ! error message
TYPE(PARATT)                           :: PARAM_META ! parameter metadata
! define output files
INTEGER(I4B)                           :: ONEMOD  ! index for defining output file (one file per model)
! SCE variables
REAL(MSP), DIMENSION(16)               :: A       ! parameter set
REAL(MSP)                              :: AF      ! objective function value
REAL(MSP), DIMENSION(16)               :: BL      ! lower bound of model parameters
REAL(MSP), DIMENSION(16)               :: BU      ! upper bound of model parameters
INTEGER(I4B)                           :: NOPT    ! number of parameters to be optimized
INTEGER(I4B)                           :: MAXN    ! maximum number of trials before optimization is terminated
INTEGER(I4B)                           :: KSTOP   ! number of shuffling loops the value must change by PCENTO
REAL(MSP)                              :: PCENTO  ! the percentage
INTEGER(I4B)                           :: ISEED   ! starting seed for the random sequence
CHARACTER(LEN=3)                       :: CSEED   ! starting seed converted to a character
INTEGER(I4B)                           :: NGS     ! # complexes in the initial population
INTEGER(I4B)                           :: NPG     ! # points in each complex
INTEGER(I4B)                           :: NPS     ! # points in a sub-complex
INTEGER(I4B)                           :: NSPL    ! # evolution steps allowed for each complex before shuffling
INTEGER(I4B)                           :: MINGS   ! minimum number of complexes required
INTEGER(I4B)                           :: INIFLG  ! 1 = include initial point in the population
INTEGER(I4B)                           :: IPRINT  ! 0 = supress printing
INTEGER(I4B)                           :: ISCE    ! unit number for SCE write
REAL(MSP)                              :: FUNCTN  ! function name for the model run
! ---------------------------------------------------------------------------------------
! (1) GET COMMAND-LINE ARGUMENTS...
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
! read model numerix parameters
CALL GETNUMERIX()         ! defines method/parameters used for numerical solution
! process command-line arguments
READ(FMODEL_ID,*) FUSE_ID ! integer definining FUSE model
READ(NSOLUTION,*) SOLUTION_METHOD         ! numerical solution (0=implicit, 1=explicit)
READ(FADAPTIVE,*) TEMPORAL_ERROR_CONTROL  ! identifier for adaptive sub-steps (0=fixed, 1=adaptive)
READ(TRUNC_ABS,*) ERR_TRUNC_ABS           ! absolute temporal truncation error tolerance
READ(TRUNC_REL,*) ERR_TRUNC_REL           ! relative temporal truncation error tolerance
! ---------------------------------------------------------------------------------------
! (2) GET MODEL SETUP -- MODEL DEFINITION, AND PARAMETER AND VARIABLE INFO FOR ALL MODELS
! ---------------------------------------------------------------------------------------
! Read data from the "BATEA-compliant" ASCII files
CALL GETFORCING(INFERN_START,NTIM)
! Define model attributes (valid for all models)
CALL UNIQUEMODL(NMOD)    ! get nmod unique models
CALL GETPARMETA()        ! parameter meta data (parameter bounds, etc.)
! Identify a single model
CALL SELECTMODL(FUSE_ID,ISTATUS=ERR,MESSAGE=MESSAGE)
IF (ERR.NE.0) THEN
 PRINT *, TRIM(MESSAGE); STOP
ENDIF
! Define list of states and parameters for the current model
CALL ASSIGN_STT()        ! state definitions are stored in module model_defn
CALL ASSIGN_FLX()        ! flux definitions are stored in module model_defn
CALL ASSIGN_PAR()        ! parameter defintions are stored in module multiparam
! Get parameter bounds and a default parameter set
IF (NUMPAR.GT.16) STOP ' NUMBER OF PARAMETERS MUST NOT EXCEED 16 IN SCE '
DO I=1,NUMPAR
 CALL GETPAR_STR(TRIM(LPARAM(I)%PARNAME),PARAM_META)
 BL(I) = PARAM_META%PARLOW
 BU(I) = PARAM_META%PARUPP
 A(I)  = PARAM_META%PARDEF
END DO
! --------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------
! loop through different starting seeds
DO ISEED=10,100,10
 ! get the seed as a character string
 WRITE(CSEED,'(i3.3)') ISEED
 ! --------------------------------------------------------------------------------------
 ! (3) DEFINE NETCDF OUTPUT FILES 
 ! --------------------------------------------------------------------------------------
 ! Define output file names
 FNAME_NETCDF = TRIM(OUTPUT_PATH)//TRIM(SMODL%MNAME)//'__'//&
                TRIM(NSOLUTION)//'-'//TRIM(FADAPTIVE)//'_SCE_'//CSEED//'.nc'  ! shared in MODULE model_defn
 FNAME_ASCII  = TRIM(OUTPUT_PATH)//TRIM(SMODL%MNAME)//'__'//&
                TRIM(NSOLUTION)//'-'//TRIM(FADAPTIVE)//'_SCE_'//CSEED//'.dat' ! shared in MODULE model_defn
 ! Define NetCDF output files (only write parameters and summary statistics)
 ONEMOD=1                 ! one file per model (i.e., model dimension = 1)
 PCOUNT=0                 ! counter for parameter sets evaluated (shared in MODULE multistats)
 CALL DEF_PARAMS(ONEMOD)  ! define model parameters (initial CREATE)
 !CALL DEF_OUTPUT(NTIM)    ! define model output (REDEF)
 CALL DEF_SSTATS()        ! define summary statistics (REDEF)
 ! --------------------------------------------------------------------------------------
 ! (4) SCE WRAPPER
 ! --------------------------------------------------------------------------------------
 ! assign algorithmic control parameters for SCE
 NOPT   =  NUMPAR         ! number of parameters to be optimized (NUMPAR in module multiparam)
 MAXN   =   1000          ! maximum number of trials before optimization is terminated
 KSTOP  =      9          ! number of shuffling loops the value must change by PCENTO (MAX=9)
 PCENTO =      0.001      ! the percentage
 NGS    =     10          ! number of complexes in the initial population
 NPG    =  2*NOPT + 1     ! number of points in each complex
 NPS    =    NOPT + 1     ! number of points in a sub-complex
 NSPL   =  2*NOPT + 1     ! number of evolution steps allowed for each complex before shuffling
 MINGS  =  NGS            ! minimum number of complexes required
 INIFLG =  1              ! 1 = include initial point in the population
 IPRINT =  1              ! 0 = supress printing
 ! open up ASCII output file
 ISCE = 96; OPEN(ISCE,FILE=TRIM(FNAME_ASCII))
 ! optimize (returns A and AF)
 CALL SCEUA(A,AF,BL,BU,NOPT,MAXN,KSTOP,PCENTO,ISEED,&
            NGS,NPG,NPS,NSPL,MINGS,INIFLG,IPRINT,ISCE)
 ! close ASCII output file
 CLOSE(ISCE)
 ! call the function again with the optimized parameter set (to ensure the last parameter set is the optimum(
 AF = FUNCTN(NOPT,A) 
 ! --------------------------------------------------------------------------------------
END DO  ! looping through seeds
! ---------------------------------------------------------------------------------------
STOP
END
