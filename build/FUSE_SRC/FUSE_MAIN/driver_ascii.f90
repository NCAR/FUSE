PROGRAM DRIVER_ASCII
! ---------------------------------------------------------------------------------------
! Creator:
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! Simple driver program for FUSE (output ASCII files)
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
! get command-line arguments
CHARACTER(LEN=11)                      :: PAR_IDX ! start index of parameter set
CHARACTER(LEN=11)                      :: PAR_JDX ! end index of parameter set
INTEGER(I4B)                           :: IPAR1   ! start index of parameter set
INTEGER(I4B)                           :: IPAR2   ! end index of parameter set
! get forcing data
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
! generate a new parameter set
INTEGER(I4B)                           :: IPAR    ! loop through model parameters
INTEGER(I4B)                           :: JPAR    ! loop through model parameters
TYPE(PARATT)                           :: PARAM_META ! parameter metadata (model parameters)
REAL(SP), DIMENSION(:), ALLOCATABLE    :: BL      ! vector of lower parameter bounds
REAL(SP), DIMENSION(:), ALLOCATABLE    :: BU      ! vector of upper parameter bounds
REAL(SP), DIMENSION(:), ALLOCATABLE    :: APAR    ! model parameter set
INTEGER(KIND=4)                        :: ISEED   ! seed for the random sequence
REAL(KIND=4),DIMENSION(:), ALLOCATABLE :: URAND   ! vector of quasi-random numbers U[0,1]
! ---------------------------------------------------------------------------------------
! (0) RETRIEVE COMMAND-LINE ARGUMENTS 
! ---------------------------------------------------------------------------------------
! get start index for parameter set
CALL GETARG(1,PAR_IDX)
IF (LEN_TRIM(PAR_IDX).EQ.0) STOP ' need start index for parameter set as 1st command-line argument '
READ(PAR_IDX,*) IPAR1   ! convert index to an integer
! get end index for parameter set
CALL GETARG(2,PAR_JDX)
IF (LEN_TRIM(PAR_JDX).EQ.0) STOP ' need end index for parameter set as 2nd command-line argument '
READ(PAR_JDX,*) IPAR2   ! convert index to an integer
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
CALL ASSIGN_PAR()        ! parameter defintions are stored in module multiparam
! compute derived model parameters (bucket sizes, etc.)
CALL PAR_DERIVE()
! --------------------------------------------------------------------------------------
! (3) DEFINE NETCDF OUTPUT FILES 
! --------------------------------------------------------------------------------------
! Define output file names (shared in MODULE model_defn)
FNAME_NETCDF = TRIM(OUTPUT_PATH)//TRIM(SMODL%MNAME)//'_'//TRIM(PAR_IDX)//'.nc'
FNAME_ASCII  = TRIM(OUTPUT_PATH)//TRIM(SMODL%MNAME)//'_'//TRIM(PAR_IDX)//'.dat'
! Define indices and flags for model output
ONEMOD=1                 ! one file per model (i.e., model dimension = 1)
PCOUNT=0                 ! counter for parameter sets evaluated (shared in MODULE multistats)
OUTPUT_FLAG = .FALSE.    ! .TRUE. if desire NetCDF time series output
SSTATS_FLAG = .FALSE.    ! .TRUE. if desire NETCDF summary statistics
! open output file
OPEN(UNIT=OUTFILE_UNIT,NAME=TRIM(FNAME_ASCII),STATUS='unknown')
! --------------------------------------------------------------------------------------
! (4) RUN MODEL
! --------------------------------------------------------------------------------------
! allocate space for parameter vectors
ALLOCATE(APAR(NUMPAR),BL(NUMPAR),BU(NUMPAR),URAND(NUMPAR))
! get parameter bounds
DO IPAR=1,NUMPAR
 CALL GETPAR_STR(LPARAM(IPAR)%PARNAME,PARAM_META) 
 BL(IPAR) = PARAM_META%PARLOW
 BU(IPAR) = PARAM_META%PARUPP
END DO
! loop through parameter sets
DO IPAR=IPAR1,IPAR2
 ISEED=IPAR
 ! get new parameter sets
 CALL I4_SOBOL(NUMPAR,ISEED,URAND)
 WRITE(*,'(I4,1X,12(E10.2,1X))') ISEED-1, URAND
 APAR = BL + URAND*(BU-BL)
 CALL PUT_PARSET(APAR)
 ! write parameter set to the file
 WRITE(OUTFILE_UNIT,'(20(A9,1X))')   (TRIM(LPARAM(JPAR)%PARNAME),JPAR=1,NUMPAR)
 WRITE(OUTFILE_UNIT,'(20(F9.3,1X))') (APAR(JPAR),JPAR=1,NUMPAR)
 ! run zee model
 CALL FMODEL_RUN_ASCII()
END DO
! close the output file
CLOSE(OUTFILE_UNIT)
STOP
END PROGRAM DRIVER_ASCII
