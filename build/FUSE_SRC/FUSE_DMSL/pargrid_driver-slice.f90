PROGRAM PARGRID_DRIVER
! ---------------------------------------------------------------------------------------
! Creator:
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! Driver program for a parameter grid
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
CHARACTER(LEN=6)                       :: NSOLUTION='      ' ! numerical solution (0=implicit, 1=explicit)
CHARACTER(LEN=6)                       :: FADAPTIVE='      ' ! identifier for adaptive sub-steps (0=fixed, 1=adaptive)
CHARACTER(LEN=6)                       :: TRUNC_ABS='      ' ! absolute temporal truncation error tolerance
CHARACTER(LEN=6)                       :: TRUNC_REL='      ' ! relative temporal truncation error tolerance
CHARACTER(LEN=11)                      :: PAR_NAME1='           ' ! name of the 1st parameter in the grid
CHARACTER(LEN=11)                      :: PAR_NAME2='           ' ! name of the 2nd parameter in the grid
! ---------------------------------------------------------------------------------------
! (1) SETUP MODELS FOR SIMULATION -- POPULATE DATA STRUCTURES
! ---------------------------------------------------------------------------------------
! get model forcing data
INTEGER(I4B)                           :: NTIM            ! number of time steps
INTEGER(I4B)                           :: INFERN_START    ! start of inference period
! get model setup
INTEGER(I4B)                           :: I               ! looping
INTEGER(I4B)                           :: NMOD            ! number of models
INTEGER(I4B)                           :: ERR             ! error code
CHARACTER(LEN=256)                     :: MESSAGE         ! error message
! define model output
LOGICAL(LGT)                           :: OUTPUT_FLAG     ! .TRUE. = write time series output
! ---------------------------------------------------------------------------------------
! (2) CREATE PARAMETER GRID
! ---------------------------------------------------------------------------------------
! Define error code for I/O
INTEGER(I4B)                           :: IERR            ! error code for I/O
! Identify index of the parameter set
INTEGER(I4B)                           :: IPARSET         ! parameter set index
CHARACTER(LEN=4)                       :: CPARSET         ! convert parameter set index to a string
! Define the number of points in each direction
INTEGER(I4B),PARAMETER                 :: NGRID=1001       ! number of samples across a single parameter dimension
! Looping variables
INTEGER(I4B)                           :: IPAR            ! index of 1st model parameter
INTEGER(I4B)                           :: JPAR            ! index of 2nd model parameter
INTEGER(I4B)                           :: KPAR            ! loop through model parameters
INTEGER(I4B)                           :: MPAR            ! loop through model parameter values
INTEGER(I4B)                           :: NPAR            ! loop through model parameter values
! Identify the initial parameter set
TYPE(PARATT)                           :: PARAM_META      ! parameter metadata (model parameters)
INTEGER(I4B)                           :: ONEMOD          ! index of the model used (=1)
! Parameter vectors
REAL(SP),DIMENSION(:),ALLOCATABLE      :: X0I             ! parameter vector
REAL(SP),DIMENSION(:),ALLOCATABLE      :: XLO             ! lower bound on solution, either none or both bounds must be present
REAL(SP),DIMENSION(:),ALLOCATABLE      :: XHI             ! upper bound on solution, either none or both bounds must be present
REAL(SP),DIMENSION(:),ALLOCATABLE      :: XDF             ! default parameter vector
REAL(SP)                               :: FPAR            ! function value for parameter set
! ---------------------------------------------------------------------------------------
! (0) READ COMMAND LINE ARGUMENTS
! ---------------------------------------------------------------------------------------
! read command-line arguments
CALL GETARG(1,NSOLUTION)  ! numerical solution (0=explicit, 1=implicit)
CALL GETARG(2,FADAPTIVE)  ! identifier for adaptive sub-steps (0=fixed, 1=adaptive)
CALL GETARG(3,TRUNC_ABS)  ! absolute temporal truncation error tolerance
CALL GETARG(4,TRUNC_REL)  ! relative temporal truncation error tolerance
CALL GETARG(5,PAR_NAME1)  ! name of the 1st parameter in the grid
CALL GETARG(6,PAR_NAME2)  ! name of the 2nd parameter in the grid
! check command-line arguments
IF (LEN_TRIM(NSOLUTION).EQ.0) STOP '1st command-line argument is missing (NSOLUTION)'
IF (LEN_TRIM(FADAPTIVE).EQ.0) STOP '2nd command-line argument is missing (FADAPTIVE)'
IF (LEN_TRIM(TRUNC_ABS).EQ.0) STOP '3rd command-line argument is missing (TRUNC_ABS)'
IF (LEN_TRIM(TRUNC_REL).EQ.0) STOP '4th command-line argument is missing (TRUNC_REL)'
IF (LEN_TRIM(PAR_NAME1).EQ.0) STOP '5th command-line argument is missing (PAR_NAME1)'
IF (LEN_TRIM(PAR_NAME2).EQ.0) STOP '6th command-line argument is missing (PAR_NAME2)'
! convert command-line arguments to integer flags and real numbers
CALL GETNUMERIX()                         ! defines method/parameters used for numerical solution
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
! Read data from the "BATEA-compliant" ASCII files
CALL GETFORCING(INFERN_START,NTIM)
!INFERN_START=1; NTIM=20; NUMTIM=NTIM; DELTIM=1._SP
!ALLOCATE(AFORCE(NTIM),AROUTE(NTIM))  ! (shared in module multiroute)
!AFORCE(INFERN_START:NTIM)%PPT   = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,50.,50.,50.,50.,50.,0.,0.,0.,0.,0./)
!AFORCE(INFERN_START:NTIM)%PET   = (/5.,5.,5.,5.,5.,5.,5.,5.,5.,5., 5., 5., 5., 5., 5.,5.,5.,5.,5.,5./)
! Define model attributes (valid for all models)
CALL UNIQUEMODL(NMOD)    ! get nmod unique models
CALL GETPARMETA()        ! read parameter metadata (parameter bounds etc.)
! Identify a single model (read control file ../DataFiles/m_decisions.txt)
CALL SELECTMODL(ISTATUS=ERR,MESSAGE=MESSAGE)
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
! (2) MAKE A PARAMETER GRID
! ---------------------------------------------------------------------------------------
! allocate arrays
ALLOCATE(X0I(NUMPAR),XLO(NUMPAR),XHI(NUMPAR),XDF(NUMPAR), STAT=IERR)
IF (IERR.NE.0) STOP ' problem allocating space for parameter arrays '
! get parameter bounds
DO IPAR=1,NUMPAR
 CALL GETPAR_STR(LPARAM(IPAR)%PARNAME,PARAM_META)
 XLO(IPAR) = PARAM_META%PARLOW
 XHI(IPAR) = PARAM_META%PARUPP
END DO
IPARSET = 0
! loop through example parameter sets
OPEN(21,FILE=TRIM(DATA_PATH)//'param_example.dat')
DO
 ! read parameter set
 READ(21,*,IOSTAT=IERR) XDF
 WRITE(*,'(20(A,1X))') LPARAM(1:NUMPAR); WRITE(*,'(20(F9.3,1X))') XDF
 IF (IERR.NE.0) EXIT
 ! increment counter
 IPARSET = IPARSET + 1
 ! convert counter to a character string
 CPARSET='    '; WRITE(CPARSET,'(I4)') IPARSET; CPARSET=ADJUSTR(CPARSET)
 FORALL(I=1:LEN(CPARSET)-LEN_TRIM(ADJUSTL(CPARSET))) CPARSET(I:I)='0'
 ! define NetCDF files (filename shared in MODULE model_defn)
 FNAME_NETCDF = TRIM(OUTPUT_PATH)//TRIM(SMODL%MNAME)//'__'//CPARSET//'__'//&
                TRIM(NSOLUTION)//'-'//TRIM(FADAPTIVE)//'__'//&
                TRIM(TRUNC_ABS)//'-'//TRIM(TRUNC_REL)//'__'//&
                TRIM(PAR_NAME1)//'-'//TRIM(PAR_NAME2)//'__parslice.nc'
 write(*,'(a)') trim(fname_netcdf)
 ONEMOD=1                 ! one file per model (i.e., model dimension = 1)
 PCOUNT=0                 ! counter for parameter sets in output file (shared in MODULE multistats)
 FCOUNT=0                 ! counter for parameter sets evaluated (shared in MODULE multistats)
 OUTPUT_FLAG = .FALSE.    ! write model output
 CALL DEF_PARAMS(ONEMOD)  ! define model parameters (initial CREATE)
 IF (OUTPUT_FLAG) CALL DEF_OUTPUT(NTIM)    ! define model time series (REDEF)
 CALL DEF_SSTATS()        ! define summary statistics (REDEF)
 X0I = XDF  ! set parameters to their default value
 ! initial run with default parameter sets
 !CALL FUSE_RMSE(X0I,FPAR,OUTPUT_FLAG)
 !PAUSE 
 ! identify IPAR and JPAR
 DO KPAR=1,NUMPAR
  IF (TRIM(LPARAM(KPAR)%PARNAME).EQ.TRIM(PAR_NAME1)) IPAR = KPAR
  IF (TRIM(LPARAM(KPAR)%PARNAME).EQ.TRIM(PAR_NAME2)) JPAR = KPAR
 END DO
 ! loop through parameter perturbations
 !DO MPAR=1,NGRID
  DO NPAR=1,NGRID
   ! perturb parameters
   !X0I(IPAR) = XLO(IPAR) + REAL(MPAR-1,KIND(SP))/REAL(NGRID-1,KIND(SP)) * (XHI(IPAR)-XLO(IPAR))
   X0I(JPAR) = XLO(JPAR) + REAL(NPAR-1,KIND(SP))/REAL(NGRID-1,KIND(SP)) * (XHI(JPAR)-XLO(JPAR))
   ! run model (parameters and statistics are written in FUSE_RMSE)
   CALL FUSE_RMSE(X0I,FPAR,OUTPUT_FLAG)
   write(*,'(i6,1x,20(f9.3,1x))') PCOUNT, X0I
  END DO  ! npar
 !END DO  ! mpar
END DO  ! looping through example parameter sets
DEALLOCATE(X0I,XLO,XHI,XDF, STAT=IERR)
IF (IERR.NE.0) STOP ' problem deallocating space for parameter arrays '
STOP
END PROGRAM PARGRID_DRIVER
! --------------------------------------------------------------------------------------
