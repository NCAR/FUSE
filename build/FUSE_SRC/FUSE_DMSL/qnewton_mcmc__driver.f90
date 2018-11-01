
!******************************************************************
!module softwareData
! Purpose:
! Programmer: Dmitri Kavetski
! Last modified:
! Comments:
!use DFWIN,only:GetCurrentProcessId
!implicit none
! type definitions
! variable definitions
!integer(INT_PTR_KIND())::myProcID
!----------------------------------------------------
!contains
!----------------------------------------------------
!function getMyProcID()
! Purpose: Returns the processID of the callling process.
! Programmer: Dmitri Kavetski
! Last modified:
! Performance
! IN:
! OUT:
! Comments:
!use DFWIN,only:GetCurrentProcessId
!implicit none
! dummies
!integer(INT_PTR_KIND())::getMyProcID
! Start procedure here
!getMyProcID=GetCurrentProcessId()
! End procedure here
!endfunction getMyProcID
!----------------------------------------------------
!endmodule softwareData
!******************************************************************



PROGRAM QNEW_MCMC__DRIVER
! ---------------------------------------------------------------------------------------
! Creator:
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! Driver program for multi-start quasi-newton optimization
! ---------------------------------------------------------------------------------------
USE nrtype                                                ! variable types, etc.
USE ddirectory                                            ! directory for data files
! data modules
USE model_defn,nstateFUSE=>nstate                         ! model definition structures
USE multiforce, ONLY: DELTIM                              ! data interval = maximum model time step
USE multiparam, ONLY: LPARAM, PARATT, NUMPAR              ! parameter metadata structures
USE multiroute                                            ! model routing structures
USE multistats                                            ! model statistics structures
! informational modules
USE selectmodl_module                                     ! reads model control file
USE getpar_str_module                                     ! extracts parameter metadata
USE par_insert_module                                     ! inserts model parameters
USE get_objfnc_module                                     ! wrapper to get objective function from NetCDF output files
USE fuse_fileManager,only: Q_ONLY                         ! only write streamflow to output file?
! model numerix
USE model_numerix                                         ! defines decisions on model numerix
! access to qnewton and model simulation modules
USE dmsl_wrapper_module                                   ! wrapper for dmsl
USE fuse_rmse_module                                      ! run model and compute the root mean squared error
! software settings (Windows only)
!use softwareData
IMPLICIT NONE
! ---------------------------------------------------------------------------------------
! (0) GET COMMAND-LINE ARGUMENTS...
! ---------------------------------------------------------------------------------------
LOGICAL(LGT)                           :: READ_ARG           ! .true. to read command-line arguments
CHARACTER(LEN=12)                      :: MBASIN_ID='            ' ! MOPEX basin ID
CHARACTER(LEN=6)                       :: FMODEL_ID='      ' ! integer defining FUSE model
CHARACTER(LEN=6)                       :: NSOLUTION='      ' ! numerical solution (0=implicit, 1=explicit)
CHARACTER(LEN=6)                       :: FADAPTIVE='      ' ! identifier for adaptive sub-steps (0=fixed, 1=adaptive)
CHARACTER(LEN=6)                       :: TRUNC_ABS='      ' ! absolute temporal truncation error tolerance
CHARACTER(LEN=6)                       :: TRUNC_REL='      ' ! relative temporal truncation error tolerance
CHARACTER(LEN=6)                       :: NUM_MULTI='      ' ! number of multiple re-starts
CHARACTER(LEN=6)                       :: SOBOLSEED='      ' ! starting seed in the Sobol sequence
CHARACTER(LEN=6)                       :: NUMDIGITS='      ' ! number of reliable digits in function evaluation
CHARACTER(LEN=6)                       :: DO_QNEWTN='      ' ! T means do the quasi-Newton
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
! (2) MULTI-START QUASI-NETWON OPTIMIZATION
! ---------------------------------------------------------------------------------------
! Check if there is a need to run the multi-start qNewton method
LOGICAL(LGT)                           :: QNEW_FLAG       ! .TRUE. means run multi-start qNewton
CHARACTER(LEN=32)                      :: OF_NAME         ! name of the desired objective function
REAL(SP), DIMENSION(:), ALLOCATABLE    :: OF_VALS         ! objective functioni values
! Control of the multi-start method
INTEGER(I4B)                           :: NMULTI          ! number of multiple re-starts
INTEGER(I4B)                           :: IBEGIN          ! starting seed in the Sobol sequence
! Define file unit
INTEGER(I4B), PARAMETER                :: UOUT_QNEW=21    ! output unit for run-time information (quasi-newton)
! Looping variables
INTEGER(I4B)                           :: ISEED           ! loop through seeds in the Sobol sequence
INTEGER(I4B)                           :: IPAR            ! loop through model parameters
! Identify the initial parameter set
INTEGER(KIND=4)                        :: JSEED           ! index in the Sobol sequence
REAL(KIND=4),DIMENSION(:),ALLOCATABLE  :: URAND           ! vector of uniform random numbers (from the Sobol sequence)
TYPE(PARATT)                           :: PARAM_META      ! parameter metadata (model parameters)
REAL(SP),PARAMETER                     :: PSELECT=0.9_SP  ! fraction of parameter space to select initial seed
INTEGER(I4B)                           :: ONEMOD          ! index of the model used (=1)
! Input to qNewton
REAL(SP),DIMENSION(:),ALLOCATABLE      :: X0I             ! initial estimate of solution
REAL(SP),DIMENSION(:),ALLOCATABLE      :: XLO             ! lower bound on solution, either none or both bounds must be present
REAL(SP),DIMENSION(:),ALLOCATABLE      :: XHI             ! upper bound on solution, either none or both bounds must be present
REAL(SP),DIMENSION(:),ALLOCATABLE      :: XSC             ! typical scale of the parameters
INTEGER(I4B)                           :: FDIGITS         ! number of reliable digits in function evaluation
!*****                                                    ! (-2=estimate,-1=full machine precision)
! Approximate optimal solution
REAL(SP),DIMENSION(:),ALLOCATABLE      :: XOPT            ! optimum value of "x", for which f(x) takes its minimum value
REAL(SP)                               :: FOPT            ! function value at optimum
! Computational cost report
INTEGER(I4B)                           :: ITER            ! number of steps (iterations)
INTEGER(I4B)                           :: FCALLS          ! number of function calls
INTEGER(I4B)                           :: GCALLS          ! number of gradient calls
INTEGER(I4B)                           :: HCALLS          ! number of Hessian calls
! ---------------------------------------------------------------------------------------
! (2) MONTE-CARLO MARKOV CHAINS
! ---------------------------------------------------------------------------------------
! Define initial sample and diagonal of the covariance matrix
real(mrk),dimension(:),allocatable     :: sample0,covDiag0
! Define files
CHARACTER(LEN=256)                     :: FNAME_PRODKT    ! name of MCMC production file
LOGICAL(LGT)                           :: LEXIST          ! logical flag if the file exists
INTEGER(I4B), PARAMETER                :: UIN_MCMC=31     ! input unit for MCMC production files
LOGICAL(LGT)                           :: JUMP_FLAG       ! flag to denote a jump in MCMC
INTEGER(I4B), DIMENSION(1)             :: IMIN            ! index of seed with highest OF value
REAL(SP),DIMENSION(:,:),ALLOCATABLE    :: XPAR            ! parameter sets for all local optima
INTEGER(I4B)                           :: IWANT           ! used to skip parameter sets



!MyProcID=GetCurrentProcessId()
! ---------------------------------------------------------------------------------------
! (1) READ COMMAND LINE ARGUMENTS
! ---------------------------------------------------------------------------------------
! read command-line arguments
read_arg=.true.
if(read_arg)then
  CALL GETARG( 1,MBASIN_ID)  ! MOPEX basin ID
  CALL GETARG( 2,FMODEL_ID)  ! integer defining FUSE model
  CALL GETARG( 3,NSOLUTION)  ! numerical solution (0=explicit, 1=implicit)
  CALL GETARG( 4,FADAPTIVE)  ! identifier for adaptive sub-steps (0=fixed, 1=adaptive)
  CALL GETARG( 5,TRUNC_ABS)  ! absolute temporal truncation error tolerance
  CALL GETARG( 6,TRUNC_REL)  ! relative temporal truncation error tolerance
  CALL GETARG( 7,NUM_MULTI)  ! number of re-starts
  CALL GETARG( 8,SOBOLSEED)  ! starting seed in the Sobol sequence
  CALL GETARG( 9,NUMDIGITS)  ! number of reliable digits in function evaluation
  CALL GETARG(10,DO_QNEWTN)  ! T = run multi-start quasi-Newton
else
  MBASIN_ID="mahurangi"
  FMODEL_ID="070"
  NSOLUTION="2"  ! implicit Euler
  FADAPTIVE="0"  ! fixed-step
  TRUNC_ABS="1.e-2"
  TRUNC_REL="1.e-2"
  NUM_MULTI="2"
  SOBOLSEED="1"
  NUMDIGITS="10"
  DO_QNEWTN="T"
endif
! check command-line arguments
IF (LEN_TRIM(MBASIN_ID).EQ.0) STOP '1st command-line argument is missing (MBASIN_ID)'
IF (LEN_TRIM(FMODEL_ID).EQ.0) STOP '2nd command-line argument is missing (FMODEL_ID)'
IF (LEN_TRIM(NSOLUTION).EQ.0) STOP '3rd command-line argument is missing (NSOLUTION)'
IF (LEN_TRIM(FADAPTIVE).EQ.0) STOP '4th command-line argument is missing (FADAPTIVE)'
IF (LEN_TRIM(TRUNC_ABS).EQ.0) STOP '5th command-line argument is missing (TRUNC_ABS)'
IF (LEN_TRIM(TRUNC_REL).EQ.0) STOP '6th command-line argument is missing (TRUNC_REL)'
IF (LEN_TRIM(NUM_MULTI).EQ.0) STOP '7th command-line argument is missing (NUM_MULTI)'
IF (LEN_TRIM(SOBOLSEED).EQ.0) STOP '8th command-line argument is missing (SOBOLSEED)'
IF (LEN_TRIM(NUMDIGITS).EQ.0) STOP '9th command-line argument is missing (NUMDIGITS)'
! define basin desired
FORCINGINFO = 'forcinginfo.'//TRIM(MBASIN_ID)//'.all'
! convert command-line arguments to integer flags and real numbers
CALL GETNUMERIX()                         ! defines method/parameters used for numerical solution
READ(FMODEL_ID,*) FUSE_ID                 ! integer definining FUSE model
READ(NSOLUTION,*) SOLUTION_METHOD         ! numerical solution (0=implicit, 1=explicit)
READ(FADAPTIVE,*) TEMPORAL_ERROR_CONTROL  ! identifier for adaptive sub-steps (0=fixed, 1=adaptive)
READ(TRUNC_ABS,*) ERR_TRUNC_ABS           ! absolute temporal truncation error tolerance
READ(TRUNC_REL,*) ERR_TRUNC_REL           ! relative temporal truncation error tolerance
READ(NUM_MULTI,*) NMULTI                  ! define the number of re-starts
READ(SOBOLSEED,*) IBEGIN                  ! starting seed in the Sobol sequence
READ(NUMDIGITS,*) FDIGITS                 ! number of reliable digits in function evaluation
! check if there is a need to run the multi-start quasi-Newton method
QNEW_FLAG=.FALSE.
IF (LEN_TRIM(DO_QNEWTN).EQ.1) THEN
 IF (DO_QNEWTN.EQ.'T') QNEW_FLAG=.TRUE.
ENDIF
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
IF (NMULTI.LE.0) STOP 'number of re-starts (6th command line argument) must be > 0'
IF (IBEGIN.LE.0) STOP 'starting seed in the Sobol sequence must be greater > 0'
write(*,'(A5,1X,2(I1,1X),2(E12.5,1X),I6,1X,A11,1X,2(I6,1X))') 'FUSE ', &
SOLUTION_METHOD, TEMPORAL_ERROR_CONTROL, ERR_TRUNC_ABS, ERR_TRUNC_REL, &
NMULTI, TRIM(SOBOLSEED), IBEGIN, FDIGITS
! ---------------------------------------------------------------------------------------
! (1) GET MODEL SETUP -- MODEL DEFINITION, AND PARAMETER AND VARIABLE INFO FOR ALL MODELS
! ---------------------------------------------------------------------------------------
! Read data from the "BATEA-compliant" ASCII files
CALL GETFORCING(INFERN_START,NTIM)
! Define model attributes (valid for all models)
CALL UNIQUEMODL(NMOD)    ! get nmod unique models
CALL GETPARMETA()        ! read parameter metadata (parameter bounds etc.)
! Identify a single model (use FUSE_ID instead of reading ../input/m_decisions.txt)
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
! allocate arrays for quasi-Newton
ALLOCATE(X0I(NUMPAR),XLO(NUMPAR),XHI(NUMPAR),XSC(NUMPAR),URAND(NUMPAR),XOPT(NUMPAR))
! allocate arrays for MCMC
allocate(sample0(0:numpar),covDiag0(0:numpar))
! get parameter bounds and random numbers
DO IPAR=1,NUMPAR
 CALL GETPAR_STR(LPARAM(IPAR)%PARNAME,PARAM_META)
 XLO(IPAR) = PARAM_META%PARLOW  ! lower bound
 XHI(IPAR) = PARAM_META%PARUPP  ! upper bound
 XSC(IPAR) = PARAM_META%PARSCL  ! typical scale
END DO
! --------------------------------------------------------------------------------------
! (2) MULTI START QUASI-NEWTON...
! --------------------------------------------------------------------------------------
! define the desired objective function and allocate space for the objective function values
OF_NAME = 'raw_rmse'; ALLOCATE(OF_VALS(NMULTI),XPAR(NUMPAR,NMULTI))
! loop through different starting positions (use the Sobol sequence)
DO ISEED=IBEGIN,(IBEGIN+NMULTI)-1
 ! get the seed as a character string
 WRITE(SOBOLSEED,'(i3.3)') ISEED
 ! define file prefix (add seeds)
 FNAME_PREFIX = TRIM(OUTPUT_PATH)//'DMSL_'//TRIM(MBASIN_ID)//'__'//TRIM(SMODL%MNAME)//'__'//TRIM(SOBOLSEED)//'__'//&
                TRIM(NSOLUTION)//'-'//TRIM(FADAPTIVE)//'__'//TRIM(NUMDIGITS)//'__'//&
                TRIM(TRUNC_ABS)//'__'//TRIM(TRUNC_REL)
 ! define NetCDF files (filename shared in MODULE model_defn)
 FNAME_NETCDF = TRIM(FNAME_PREFIX)//'__qnewton.nc'
 ONEMOD=1                  ! one file per model (i.e., model dimension = 1)

 ! check if there is a need to run quasi-Newton
 IF (QNEW_FLAG) THEN  ! need to run quasi-Newton
  PCOUNT=0                 ! counter for parameter sets in output file (shared in MODULE multistats)
  FCOUNT=0                 ! counter for parameter sets evaluated (shared in MODULE multistats)
  OUTPUT_FLAG = .TRUE.     ! write model time series
  Q_ONLY      = .TRUE.     ! restrict output time series to simulated runofff
  CALL DEF_PARAMS(ONEMOD)  ! define model parameters (initial CREATE)
  IF (OUTPUT_FLAG) CALL DEF_OUTPUT(NTIM)    ! define model time series (REDEF)
  CALL DEF_SSTATS()        ! define summary statistics (REDEF)
  ! define ASCII files (filename shared in MODULE model_defn)
  FNAME_ASCII  = TRIM(FNAME_PREFIX)//'__qnewton.txt'
  ! open ASCII file (unit 21)
  OPEN(UOUT_QNEW,FILE=FNAME_ASCII, STATUS='unknown')
   ! get new parameter sets
   JSEED=ISEED; CALL I4_SOBOL(NUMPAR,JSEED,URAND)
   WRITE(*,'(2(I4,1X),20(E10.2,1X))') ISEED, JSEED-1, URAND
   X0I = XLO + ((1._SP - PSELECT)/2._SP)*(XHI-XLO) + (PSELECT*REAL(URAND,KIND(SP)))*(XHI-XLO)
   ! find local optimum in the vicinity of the starting point
   CALL QNEWTON_WRAPPER(X0I,XLO,XHI,XSC,FDIGITS,UOUT_QNEW,XOPT,FOPT,ITER,FCALLS,GCALLS,HCALLS,&
                        ERR,MESSAGE)
   IF (ERR.NE.0) PRINT *, TRIM(MESSAGE)
   WRITE(*,'(5(I6,1X),20(F9.3,1X))') FCOUNT,ITER,FCALLS,GCALLS,HCALLS,FOPT,XOPT
   ! run model again with optimum parameter set (to populate structures and write model output)
   CALL FUSE_RMSE(XOPT,FOPT,OUTPUT_FLAG)
   ! write model parameters and summary statistics
   CALL PUT_PARAMS(PCOUNT,ONEMOD)  ! PCOUNT = index for parameter set; ONEMOD=1 (just one model structure)
   CALL PUT_SSTATS(PCOUNT,ONEMOD)
  CLOSE(UOUT_QNEW)
 ENDIF
 ! get objective function value for the first parameter set
 PCOUNT=1; CALL GET_OBJFNC(FNAME_NETCDF,OF_NAME,ONEMOD,PCOUNT,FOPT,XOPT)
 OF_VALS(ISEED) = FOPT
 XPAR(:,ISEED)  = XOPT(:)
 !write(*,'(I4,1X,20(f12.6,1x))') ISEED, OF_VALS(ISEED), XPAR(:,ISEED)
END DO
! --------------------------------------------------------------------------------------
! (3a) MCMC...
! --------------------------------------------------------------------------------------
! identify the maximum seed and retrieve model parameter set
IMIN = MINLOC(OF_VALS)
FOPT = OF_VALS(IMIN(1))
XOPT(:) = XPAR(:,IMIN(1))
! ensure the parameter set is within bounds
DO IPAR=1,NUMPAR
 CALL GETPAR_STR(LPARAM(IPAR)%PARNAME,PARAM_META)
 XOPT(IPAR) = MAX(XOPT(IPAR),PARAM_META%PARLOW)
 XOPT(IPAR) = MIN(XOPT(IPAR),PARAM_META%PARUPP)
END DO
!write(*,'(i3,1x,20(f12.6,1x))') IMIN(1), FOPT, XOPT
! define write parameters for model output
PCOUNT=0                 ! counter for parameter sets in output file (shared in MODULE multistats)
FCOUNT=0                 ! counter for parameter sets evaluated (shared in MODULE multistats)
OUTPUT_FLAG = .TRUE.     ! write model time series
Q_ONLY      = .TRUE.     ! restrict output time series to simulated runofff
! define file prefix (no seeds in the filename)
FNAME_PREFIX = TRIM(OUTPUT_PATH)//'DMSL_'//TRIM(MBASIN_ID)//'__'//TRIM(SMODL%MNAME)//'__'//&
               TRIM(NSOLUTION)//'-'//TRIM(FADAPTIVE)//'__'//TRIM(NUMDIGITS)//'__'//&
               TRIM(TRUNC_ABS)//'__'//TRIM(TRUNC_REL)
! define NetCDF files (filename shared in MODULE model_defn)
FNAME_NETCDF = TRIM(FNAME_PREFIX)//'__predict.nc'
CALL DEF_PARAMS(ONEMOD)  ! define model parameters (initial CREATE)
IF (OUTPUT_FLAG) CALL DEF_OUTPUT(NTIM)    ! define model time series (REDEF)
CALL DEF_SSTATS()        ! define summary statistics (REDEF)
! get MCMC samples
sample0(0)  = log10(FOPT**2)  ! variance
sample0(1:) = XOPT(1:)
!write(*,'(es22.14e3)') sample0
!write(*,'(Z16.16)') sample0
covDiag0(0)   = 0.1_mrk*max(abs(sample0(0)),1._mrk)   ! assume typical RMSE = 1 mm/day for "optimal" param set
covDiag0(1:)  = 0.1_mrk*max(abs(sample0(1:)),xsc(1:))
!CALL MCMC_WRAPPER(sample0,covDiag0,err,message)
IF (ERR.NE.0) PRINT *, TRIM(MESSAGE)
! read the ASCII production MCMC output file and re-run for each of the parameter sets
FNAME_PRODKT = TRIM(FNAME_PREFIX)//'__prodkt.sam'
INQUIRE(FILE=TRIM(FNAME_PRODKT),EXIST=LEXIST)
IF (.NOT.LEXIST) STOP ' PRODKT FILE DOES NOT EXIST '
OPEN(UIN_MCMC,FILE=TRIM(FNAME_PRODKT),IOSTAT=ERR)
IF (ERR.NE.0) THEN; PRINT *, ERR, ' PROBLEM OPENING FILE '; STOP; ENDIF
 IWANT = 1
 DO ! continuous do loop with exit clause
  ! read a parameter set
  READ(UIN_MCMC,*,IOSTAT=ERR) sample0, MSTATS%LOGP_SIMULN, JUMP_FLAG; IF (ERR.NE.0) EXIT
  MSTATS%JUMP_TAKEN = 0._SP; IF (JUMP_FLAG) MSTATS%JUMP_TAKEN = 1._SP ! (convert flag to real)
  IF (IWANT.EQ.50) THEN
   ! run FUSE
   CALL FUSE_RMSE(sample0(1:),FOPT,OUTPUT_FLAG)
   ! write model parameters and summary statistics
   CALL PUT_PARAMS(PCOUNT,ONEMOD)  ! PCOUNT = index for parameter set; ONEMOD=1 (just one model structure)
   CALL PUT_SSTATS(PCOUNT,ONEMOD)
   IWANT = 1
  ELSE
   IWANT = IWANT+1
  ENDIF
 END DO
CLOSE(UIN_MCMC)
! --------------------------------------------------------------------------------------
! deallocate parameter vectors
DEALLOCATE(X0I,XLO,XHI,XSC,URAND,XOPT,sample0,covDiag0)
STOP
END PROGRAM QNEW_MCMC__DRIVER
! --------------------------------------------------------------------------------------
