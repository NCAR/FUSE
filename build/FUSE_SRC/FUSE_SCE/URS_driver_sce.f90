PROGRAM URS_DRIVER
! ---------------------------------------------------------------------------------------
! Creator:
! Martyn Clark, 2011
! Modified by Brian Henn to include snow model, 6/2013
! Modified by Nans Addor to include SCE
! ---------------------------------------------------------------------------------------
! Purpose:
! Driver program to perform multiple runs of a model by uniform random sampling from the
!  feasible parameter space.
!
! Current purpose:
! Calibrate FUSE using SCE
! ---------------------------------------------------------------------------------------
USE nrtype                                                ! variable types, etc.
USE fuse_fileManager,only:fuse_SetDirsUndPhiles,&         ! sets directories and filenames
     SETNGS_PATH,OUTPUT_PATH,FORCINGINFO,MBANDS_INFO
! data modules
USE model_defn,nstateFUSE=>nstate                         ! model definition structures
USE model_defnames                                        ! defines the integer model options
USE multiforce, ONLY: AFORCE, DELTIM, NUMTIM              ! data interval = maximum model time step
USE multibands                                            ! basin band stuctures
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
CHARACTER(LEN=1024)                    :: FFMFILE='        ' ! name of fuse_file_manager file
CHARACTER(LEN=12)                      :: MBASIN_ID='      ' ! MOPEX basin ID
CHARACTER(LEN=6)                       :: FMODEL_ID='      ' ! integer defining FUSE model
CHARACTER(LEN=6)                       :: NSOLUTION='      ' ! numerical solution (0=explicit Euler; 1=explicit Heun; 2=implicit Euler; 3=implicit Heun, 4=semi-implicit)
CHARACTER(LEN=6)                       :: FADAPTIVE='      ' ! identifier for adaptive sub-steps (0=fixed, 1=adaptive)
CHARACTER(LEN=6)                       :: TRUNC_ABS='      ' ! absolute temporal truncation error tolerance
CHARACTER(LEN=6)                       :: TRUNC_REL='      ' ! relative temporal truncation error tolerance
CHARACTER(LEN=12)                      :: TSTEP_LEN='      ' ! maximum length of the time step (days)
CHARACTER(LEN=6)                       :: NUMPARSET='      ' ! number of parameter sets
CHARACTER(LEN=10)                      :: MAX_T='          ' ! maximum number of trials before optimization is terminated
CHARACTER(LEN=6)                       :: MAX_SL='         ' ! number of shuffling loops the objective function must change by PCENTO (max=9)
CHARACTER(LEN=10)                      :: PERC='           ' ! percentage PCENTO (1 is 1%)
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
! (3) SCE VARIABLES
! ---------------------------------------------------------------------------------------
REAL(MSP)                              :: AF      ! objective function value
INTEGER(I4B)                           :: NOPT    ! number of parameters to be optimized
INTEGER(I4B)                           :: KSTOP   ! number of shuffling loops the value must change by PCENTO
INTEGER(I4B)                           :: MAXN    ! maximum number of trials before optimization is terminated
REAL(MSP)                              :: PCENTO  ! the percentage
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
! (0) READ COMMAND LINE ARGUMENTS
! ---------------------------------------------------------------------------------------
! read command-line arguments
CALL GETARG(1,FFMFILE)    ! name of fuse_file_manager file
CALL GETARG(2,MBASIN_ID)  ! basin ID
CALL GETARG(3,FMODEL_ID)  ! integer defining FUSE model
CALL GETARG(4,NSOLUTION)  ! numerical solution (0=explicit, 1=implicit)
CALL GETARG(5,FADAPTIVE)  ! identifier for adaptive sub-steps (0=fixed, 1=adaptive)
CALL GETARG(6,TRUNC_ABS)  ! absolute temporal truncation error tolerance
CALL GETARG(7,TRUNC_REL)  ! relative temporal truncation error tolerance
CALL GETARG(8,TSTEP_LEN)  ! maximum length of the time step (days)
CALL GETARG(9,NUMPARSET)  ! number of parameter sets
CALL GETARG(10,MAX_T)     ! maximum number of trials before optimization is terminated
CALL GETARG(11,MAX_SL)    ! number of shuffling loops the objective function must change by PCENTO (max=9)
CALL GETARG(12,PERC)      ! percentage PCENTO (1 is 1%)

! check command-line arguments
IF (LEN_TRIM(FFMFILE).EQ.0)   STOP '1st command-line argument is missing (FFMFILE)'
IF (LEN_TRIM(MBASIN_ID).EQ.0) STOP '2nd command-line argument is missing (MBASIN_ID)'
IF (LEN_TRIM(FMODEL_ID).EQ.0) STOP '3rd command-line argument is missing (FMODEL_ID)'
IF (LEN_TRIM(NSOLUTION).EQ.0) STOP '4th command-line argument is missing (NSOLUTION)'
IF (LEN_TRIM(FADAPTIVE).EQ.0) STOP '5th command-line argument is missing (FADAPTIVE)'
IF (LEN_TRIM(TRUNC_ABS).EQ.0) STOP '6th command-line argument is missing (TRUNC_ABS)'
IF (LEN_TRIM(TRUNC_REL).EQ.0) STOP '7th command-line argument is missing (TRUNC_REL)'
IF (LEN_TRIM(TSTEP_LEN).EQ.0) STOP '8th command-line argument is missing (TSTEP_LEN)'
IF (LEN_TRIM(NUMPARSET).EQ.0) STOP '9th command-line argument is missing (NUMPARSET)'
IF (LEN_TRIM(MAX_T).EQ.0)  stop 	'10th command-line argument is missing (MAX_T)'
IF (LEN_TRIM(MAX_SL).EQ.0) STOP 	'11th command-line argument is missing (MAX_SL)'
IF (LEN_TRIM(PERC).EQ.0) STOP 		'12th command-line argument is missing (PERC)'

! set path to fuse_file_manager
!FFMFILE=TRIM(SETNGS_PATH)//TRIM(MBASIN_ID)//'_fuse_file_manager.txt'
print *, 'fuse_file_manager:', TRIM(FFMFILE)

! get directories and filenames for control files
call fuse_SetDirsUndPhiles(fuseFileManagerIn=FFMFILE,err=err,message=message)

if (err.ne.0) write(*,*) trim(message); if (err.gt.0) stop
CALL GETNUMERIX(ERR,MESSAGE)              ! defines method/parameters used for numerical solution
if (err.ne.0) write(*,*) trim(message); if (err.gt.0) stop
! define basin desired
! convert command-line arguments to integer flags and real numbers
FORCINGINFO = TRIM(MBASIN_ID)//'_input_info_calib.txt'
MBANDS_INFO = TRIM(MBASIN_ID)//'_elev_bands_info.txt'
READ(FMODEL_ID,*) FUSE_ID                 ! integer defining FUSE model
READ(NSOLUTION,*) SOLUTION_METHOD         ! numerical solution (0=implicit, 1=explicit)
READ(FADAPTIVE,*) TEMPORAL_ERROR_CONTROL  ! identifier for adaptive sub-steps (0=fixed, 1=adaptive)
READ(TRUNC_ABS,*) ERR_TRUNC_ABS           ! absolute temporal truncation error tolerance
READ(TRUNC_REL,*) ERR_TRUNC_REL           ! relative temporal truncation error tolerance
READ(TSTEP_LEN,*) MAX_TSTEP               ! maximum length of the time step (days)
READ(NUMPARSET,*) NUMPSET                 ! number of parameter sets
READ(MAX_T,*) MAXN                        ! maximum number of trials before optimization is terminated
READ(MAX_SL,*) KSTOP                      ! number of shuffling loops the objective function must change by PCENTO (max=9)
READ(PERC,*) PCENTO                       ! percentage PCENTO (1 is 1%)

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
! Define model attributes (valid for all models)
CALL UNIQUEMODL(NMOD)           ! get nmod unique models
CALL GETPARMETA(ERR,MESSAGE)    ! read parameter metadata (parameter bounds etc.)

IF (ERR.NE.0) WRITE(*,*) TRIM(MESSAGE); IF (ERR.GT.0) STOP
! Select model: fuse_zDecisions.txt is ignored if FUSE_ID is provided, otherwise model defined by fuse_zDecisions.txt
CALL SELECTMODL(FUSE_ID,ERR=ERR,MESSAGE=MESSAGE)
IF (ERR.NE.0) WRITE(*,*) TRIM(MESSAGE); IF (ERR.GT.0) STOP
! Define list of states and parameters for the current model
! Read data from the "BATEA-compliant" ASCII files
CALL GETFORCING(INFERN_START,NTIM) ! read forcing data
IF (SMODL%iSNOWM.EQ.iopt_temp_index) CALL GET_MBANDS(err,message) ! read band data if snow model 
CALL ASSIGN_STT()        ! state definitions are stored in module model_defn
CALL ASSIGN_FLX()        ! flux definitions are stored in module model_defn
CALL ASSIGN_PAR()        ! parameter definitions are stored in module multiparam
! compute derived model parameters (bucket sizes, etc.)
CALL PAR_DERIVE(ERR,MESSAGE)
IF (ERR.NE.0) WRITE(*,*) TRIM(MESSAGE); IF (ERR.GT.0) STOP
! Define output file names (shared in MODULE model_defn)
FNAME_NETCDF = TRIM(OUTPUT_PATH)//'sce_params_'//TRIM(MBASIN_ID)//'_'//TRIM(SMODL%MNAME)//'_'//&
               TRIM(NSOLUTION)//'_'//TRIM(FADAPTIVE)//'_'//&
               TRIM(TRUNC_ABS)//'_'//TRIM(TRUNC_REL)//'_'//&
               TRIM(TSTEP_LEN)//'_SCE_'//TRIM(NUMPARSET)//'_'//&
			   TRIM(MAX_T)//'_'//TRIM(MAX_SL)//'_'//TRIM(PERC)//'.nc' 

write(*,'(a)') trim(fname_netcdf)
! Define NetCDF output files (only write parameters and summary statistics)
ONEMOD=1                 ! one file per model (i.e., model dimension = 1)
PCOUNT=0                 ! counter for parameter sets evaluated (shared in MODULE multistats)
OUTPUT_FLAG = .FALSE.    ! .TRUE. if desire time series output
CALL DEF_PARAMS(ONEMOD)  ! define model parameters (initial CREATE)
CALL DEF_SSTATS()        ! define summary statistics (REDEF)
IF (OUTPUT_FLAG) CALL DEF_OUTPUT(NTIM)    ! define model time series (REDEF)
! --------------------------------------------------------------------------------------
! (2) SCE WRAPPER
! --------------------------------------------------------------------------------------
! assign algorithmic control parameters for SCE
NOPT   =  NUMPAR         ! number of parameters to be optimized (NUMPAR in module multiparam)
PRINT *, 'NUMBER OF PARAM PASSED TO SCE WRAPPER', NUMPAR
!MAXN   = 5000			 ! maximum number of trials before optimization is terminated
!KSTOP  =      3          ! number of shuffling loops the value must change by PCENTO (MAX=9)
!PCENTO =      0.001      ! the percentage
NGS    =     10          ! number of complexes in the initial population
NPG    =  2*NOPT + 1     ! number of points in each complex
NPS    =    NOPT + 1     ! number of points in a sub-complex
NSPL   =  2*NOPT + 1     ! number of evolution steps allowed for each complex before shuffling
MINGS  =  NGS            ! minimum number of complexes required
INIFLG =  1              ! 1 = include initial point in the population
IPRINT =  1              ! 0 = suppress printing

! get parameter bounds and random numbers
ALLOCATE(APAR(NUMPAR),BL(NUMPAR),BU(NUMPAR),URAND(NUMPAR))
DO IPAR=1,NUMPAR ! loop through parameters and retrieve 
	
 CALL GETPAR_STR(LPARAM(IPAR)%PARNAME,PARAM_META)
 
 BL(IPAR) = PARAM_META%PARLOW
 BU(IPAR) = PARAM_META%PARUPP

END DO
! loop through parameter sets
DO IPSET=1,NUMPSET
	! assign new seed
	ISEED=IPSET
	
	! get new parameter sets using seed
	CALL I4_SOBOL(NUMPAR,ISEED,URAND)
	WRITE(*,'(I4,1X,12(E10.2,1X))') ISEED-1, URAND
	APAR = BL + URAND*(BU-BL)

	 ! get the seed as a character string
	WRITE(CSEED,'(i3.3)') ISEED-1

	FNAME_ASCII = TRIM(OUTPUT_PATH)//'sce_output_'//TRIM(MBASIN_ID)//'_'//TRIM(SMODL%MNAME)//'_'//& ! shared in MODULE model_defn
	           TRIM(NSOLUTION)//'_'//TRIM(FADAPTIVE)//'_'//&
	           TRIM(TRUNC_ABS)//'_'//TRIM(TRUNC_REL)//'_'//&
	           TRIM(TSTEP_LEN)//'_SCE_'//TRIM(NUMPARSET)//'_'//&
			   TRIM(MAX_T)//'_'//TRIM(MAX_SL)//'_'//TRIM(PERC)//'_'//CSEED//'.dat'

	write(*,'(a)') trim(FNAME_ASCII)

	print *, 'NUMPAR:', NUMPAR
	print *, 'NOPT:', NOPT

	! open up ASCII output file
	ISCE = 96; OPEN(ISCE,FILE=TRIM(FNAME_ASCII))
	! optimize (returns A and AF)

	! run SCE, which will repeatedly call FUNCTN
	CALL SCEUA(APAR,AF,BL,BU,NOPT,MAXN,KSTOP,PCENTO,ISEED,& 
	        NGS,NPG,NPS,NSPL,MINGS,INIFLG,IPRINT,ISCE)
	! close ASCII output file
	CLOSE(ISCE)
	! call the function again with the optimized parameter set (to ensure the last parameter set is the optimum
	AF = FUNCTN(NOPT,APAR) 

END DO
! and, deallocate space
DEALLOCATE(APAR,BL,BU,URAND)
STOP
END PROGRAM URS_DRIVER
