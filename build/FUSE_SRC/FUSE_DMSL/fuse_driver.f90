PROGRAM DISTRIBUTED_DRIVER
! ---------------------------------------------------------------------------------------
! Creators:
! Martyn Clark, 2011
! Modified by Brian Henn to include snow model, 6/2013
! Modified by Nans Addor to include distributed modeling, 9/2016
! Modified by Nans Addor to re-enable catchment-scale modeling, 4/2017
! ---------------------------------------------------------------------------------------
! Purpose:
! Driver program to run FUSE with a snow module as either at the catchment-scale or
! at the grid-scale
! ---------------------------------------------------------------------------------------
USE nrtype                                                ! variable types, etc.
USE netcdf                                                ! NetCDF library
USE fuse_fileManager,only:fuse_SetDirsUndPhiles,&         ! sets directories and filenames
          SETNGS_PATH,MBANDS_INFO,MBANDS_NC, &
          OUTPUT_PATH,FORCINGINFO,INPUT_PATH

! data modules
USE model_defn,nstateFUSE=>nstate                         ! model definition structures
USE model_defnames                                        ! defines the integer model options
USE multiforce, ONLY: forcefile,vname_aprecip             ! model forcing structures
USE multiforce, ONLY: AFORCE, aValid                      ! time series of lumped forcing/response data
USE multiforce, ONLY: nspat1, nspat2                      ! grid dimensions
USE multiforce, ONLY: GFORCE, GFORCE_3d                   ! spatial arrays of gridded forcing data
USE multiforce, only: ancilF, ancilF_3d                   ! ancillary forcing data
USE multiforce, ONLY: valDat                              ! response data
USE multiforce, only: DELTIM
USE multiforce, only: ISTART                              ! index for start of inference
USE multiforce, only: numtim_in, itim_in                  ! length of input time series and associated index
USE multiforce, only: numtim_sim, itim_sim                ! length of simulated time series and associated index
USE multiforce, only: numtim_sub, itim_sub                ! length of subperiod time series and associated index
USE multiforce,only:  warmup_beg,infern_beg,infern_end    ! timestep indices
USE multiforce,only:  longrun_beg,longrun_end             ! timestep indices
USE multiforce,only:  sim_beg,sim_end                     ! timestep indices
USE multiforce, only: ncid_forc                           ! NetCDF forcing file ID
USE multiforce, only: ncid_var                            ! NetCDF forcing variable ID
USE multistate, only: ncid_out                            ! NetCDF output file ID
USE multiforce, only: NA_VALUE                            ! NA_VALUE for the forcing

USE multibands                                            ! basin band stuctures
USE multiparam, ONLY: LPARAM, PARATT, NUMPAR              ! parameter metadata structures
USE multistate, only: gState                              ! gridded state variables
USE multistate, only: gState_3d                           ! gridded state variables with a time dimension
USE multiroute, ONLY: AROUTE                              ! model routing structures
USE multiroute, ONLY: AROUTE_3d                           ! model routing structures with a time dimension
USE multistats                                            ! model statistics structures

! informational modules
USE selectmodl_module                                     ! reads model control file
USE getpar_str_module                                     ! extracts parameter metadata
USE par_insert_module                                     ! inserts model parameters
USE force_info_module,only:force_info                     ! get forcing info for NetCDF files
USE get_gforce_module,only:read_ginfo                     ! get dimension lengths from the NetCDF file
USE get_gforce_module,only:get_varid                      ! get netCDF ID for forcing variables
USE get_gforce_module,only:get_gforce_3d                  ! get forcing
USE get_mbands_module,only:get_mbands, GET_MBANDS_INFO    ! get elevation bands for snow modeling
USE get_fparam_module                                     ! get SCE parameters from NetCDF file
USE getf_ascii_module,only:prelim_asc                     ! get preliminary data from the ASCII file
USE getf_ascii_module,only:close_file                     ! close ASCII file
USE getf_ascii_module,only:read_ascii                     ! read ascii forcing data for a given time step

! model numerix
USE model_numerix                                         ! defines decisions on model numerix

! access to model simulation modules
USE fuse_rmse_module                                      ! run model and compute the root mean squared error
IMPLICIT NONE

! ---------------------------------------------------------------------------------------
! GET COMMAND-LINE ARGUMENTS...
! ---------------------------------------------------------------------------------------
CHARACTER(LEN=64)                      :: DatString          ! string defining forcing data
CHARACTER(LEN=64)                      :: dom_id             ! ID of the domain
CHARACTER(LEN=6)                       :: FMODEL_ID='      ' ! integer defining FUSE model
CHARACTER(LEN=10)                      :: fuse_mode='      ' ! fuse execution mode (run_def, run_best, run_pre, calib_sce)
CHARACTER(LEN=64)                      :: file_para_list     ! txt file containing list of parameter sets

! ---------------------------------------------------------------------------------------
! SETUP MODELS FOR SIMULATION -- POPULATE DATA STRUCTURES
! ---------------------------------------------------------------------------------------
! fuse_file_manager
CHARACTER(LEN=1024)                    :: FFMFILE      	  ! name of fuse_file_manager file
CHARACTER(LEN=1024)                    :: ELEV_BANDS_NC	  ! name of NetCDF file for elevation bands
! get model forcing data
INTEGER(I4B)                           :: NTIM            ! number of time steps - still needed ?
INTEGER(I4B)                           :: INFERN_START    ! start of inference period - still needed?
! get model setup
INTEGER(I4B)                           :: FUSE_ID         ! integer defining FUSE model
INTEGER(I4B)                           :: NMOD            ! number of models
INTEGER(I4B)                           :: ERR             ! error code
CHARACTER(LEN=1024)                    :: MESSAGE         ! error message
! get spatial option
CHARACTER(LEN=6)                       :: SPATIAL_OPTION  ! spatial option (catch or grid)
LOGICAL(LGT)                           :: GRID_FLAG       ! spatial flag .true. if grid
INTEGER(I4B),PARAMETER                 :: LUMPED=0        ! named variable for lumped simulations
INTEGER(I4B),PARAMETER                 :: DISTRIBUTED=1   ! named variable for distributed simulations
! define model output
LOGICAL(LGT)                           :: OUTPUT_FLAG     ! .TRUE. = write time series output
INTEGER(I4B)                           :: ONEMOD=1        ! just specify one model
! timers
INTEGER(I4B)                           :: T_start_import_forcing ! system clock
INTEGER(I4B)                           :: T_end_import_forcing   ! system clock

! ---------------------------------------------------------------------------------------
! RUN MODEL FOR DIFFERENT PARAMETER SETS
! ---------------------------------------------------------------------------------------
INTEGER(I4B)                           :: ITIM    ! loop thru time steps
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
! SCE VARIABLES
! ---------------------------------------------------------------------------------------
REAL(MSP)                              :: AF_MSP    ! objective function value
REAL(MSP), DIMENSION(:), ALLOCATABLE   :: APAR_MSP  ! ! lower bound of model parameters
REAL(MSP), DIMENSION(:), ALLOCATABLE   :: BL_MSP    ! ! lower bound of model parameters
REAL(MSP), DIMENSION(:), ALLOCATABLE   :: BU_MSP    ! ! upper bound of model parameters
REAL(MSP), DIMENSION(:), ALLOCATABLE   :: URAND_MSP   ! vector of quasi-random numbers U[0,1]
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
! READ COMMAND LINE ARGUMENTS
! ---------------------------------------------------------------------------------------
! read command-line arguments
CALL GETARG(1,DatString)  ! string defining forcinginfo file
CALL GETARG(2,dom_id)     ! ID of the domain
CALL GETARG(3,FMODEL_ID)  ! ID of FUSE model
CALL GETARG(4,fuse_mode)  ! fuse execution mode (run_def, run_best, calib_sce)
IF(TRIM(fuse_mode).EQ.'run_pre')  CALL GETARG(5,file_para_list)  ! fuse execution mode txt file containing list of parameter sets

! check command-line arguments
IF (LEN_TRIM(DatString).EQ.0) STOP '1st command-line argument is missing (DatString)'
IF (LEN_TRIM(dom_id).EQ.0) STOP '2nd command-line argument is missing (dom_id)'
IF (LEN_TRIM(FMODEL_ID).EQ.0) STOP '3rd command-line argument is missing (fuse_mode)'
IF (LEN_TRIM(fuse_mode).EQ.0) STOP '4th command-line argument is missing (fuse_mode)'
IF(TRIM(fuse_mode).EQ.'run_pre')THEN
  IF(LEN_TRIM(file_para_list).EQ.0)  STOP '5th command-line argument is missing (file_para_list) and is required in mode run_pre'
ENDIF

! print command-line arguments
print*, '1st command-line argument (DatString) = ', trim(DatString)
print*, '2nd command-line argument (dom_id) = ', trim(dom_id)
print*, '3rd command-line argument (FMODEL_ID) = ', FMODEL_ID
print*, '4th command-line argument (fuse_mode) = ', fuse_mode
IF(TRIM(fuse_mode).EQ.'run_pre')THEN
  print*, '5th command-line argument (file_para_list) = ', file_para_list
ENDIF

! set path to fuse_file_manager
FFMFILE=DatString ! must be in bin folder and you must be in bin to run FUSE - TODO read argument to FFMFILE directly

! set directories and filenames for control files
call fuse_SetDirsUndPhiles(fuseFileManagerIn=FFMFILE,err=err,message=message)
if (err.ne.0) write(*,*) trim(message); if (err.gt.0) stop

PRINT *, 'SETNGS_PATH:', TRIM(SETNGS_PATH)
PRINT *, 'INPUT_PATH:', TRIM(INPUT_PATH)
PRINT *, 'OUTPUT_PATH:', TRIM(OUTPUT_PATH)

! defines method/parameters used for numerical solution
CALL GETNUMERIX(ERR,MESSAGE)

! define basin desired - ?overwrite fuse_SetDirsUndPhiles?
FORCINGINFO = TRIM(dom_id)//'_input_info.txt'
MBANDS_INFO = TRIM(dom_id)//'_elev_bands_info.txt'
ELEV_BANDS_NC = TRIM(dom_id)//'_'//MBANDS_NC

! convert command-line arguments to integer flags and real numbers
READ(FMODEL_ID,*) FUSE_ID                 ! integer defining FUSE model

! ---------------------------------------------------------------------------------------
! GET MODEL SETUP -- MODEL DEFINITION, AND PARAMETER AND VARIABLE INFO FOR ALL MODELS
! ---------------------------------------------------------------------------------------

! get forcing info from the txt file, including NA_VALUE
call force_info(fuse_mode,err,message)
if(err/=0)then; write(*,*) trim(message); stop; endif

! allocate space for the basin-average time series
allocate(aForce(numtim_sub),aRoute(numtim_sub),stat=err)
!allocate(aForce(numtim_sub),aRoute(numtim_sub),aValid(numtim_sub),stat=err)
if(err/=0)then; write(*,*) 'unable to allocate space for basin-average time series [aForce,aRoute]'; stop; endif

print *, 'Open forcing file:', trim(INPUT_PATH)//trim(forcefile)

! open NetCDF forcing file
err = nf90_open(trim(INPUT_PATH)//trim(forcefile), nf90_nowrite, ncid_forc)
if (err.ne.0) write(*,*) trim(message); if (err.gt.0) stop
PRINT *, 'NCID_FORC is', ncid_forc

! get the grid info (spatial and temporal dimensions) from the NetCDF file
call read_ginfo(ncid_forc,err,message)
if(err/=0)then; write(*,*) trim(message); stop; endif

! allocate space for the forcing grid and states
allocate(ancilF(nspat1,nspat2), gForce(nspat1,nspat2), gState(nspat1,nspat2), stat=err)
if(err/=0)then; write(*,*) 'unable to allocate space for forcing grid GFORCE'; stop; endif

! allocate space for the forcing grid and states with a time dimension - only for subperiod
allocate(AROUTE_3d(nspat1,nspat2,numtim_sub), gState_3d(nspat1,nspat2,numtim_sub+1),gForce_3d(nspat1,nspat2,numtim_sub),aValid(nspat1,nspat2,numtim_sub),stat=err)
if(err/=0)then; write(*,*) 'unable to allocate space for 3d structure'; stop; endif

! get elevation band info, in particular N_BANDS
CALL GET_MBANDS_INFO(ELEV_BANDS_NC,err,message) ! read band data from NetCDF file - for a 2D grid

! allocate space for elevation bands
allocate(MBANDS_VAR_4d(nspat1,nspat2,N_BANDS,numtim_sub+1),stat=err)
if(err/=0)then; write(*,*) 'unable to allocate space for elevation bands'; stop; endif

! get variable ID from the NetCDF file
call get_varID(ncid_forc,err,message)
if(err/=0)then; write(*,*) 'unable to get NetCDF variables ID'; stop; endif

! define the spatial flag (.true. is distributed)
if(nSpat1.GT.1) THEN

  GRID_FLAG=.TRUE.

ELSE

  GRID_FLAG=.FALSE.

ENDIF

print*, 'spatial dimensions = ', nSpat1, nSpat2
print*, 'netCDF ID for forcing file', ncid_forc
print*, 'NA_VALUE = ', NA_VALUE
print*, 'GRID_FLAG = ', GRID_FLAG

! Define model attributes (valid for all models)
CALL UNIQUEMODL(NMOD)           ! get nmod unique models
CALL GETPARMETA(ERR,MESSAGE)    ! read parameter metadata (parameter bounds etc.)

IF (ERR.NE.0) WRITE(*,*) TRIM(MESSAGE); IF (ERR.GT.0) STOP

! Identify a single model
CALL SELECTMODL(FMODEL_ID,ERR=ERR,MESSAGE=MESSAGE)
IF (ERR.NE.0) WRITE(*,*) TRIM(MESSAGE); IF (ERR.GT.0) STOP

! Define list of states and parameters for the current model
CALL ASSIGN_STT()        ! state definitions are stored in module model_defn
CALL ASSIGN_FLX()        ! flux definitions are stored in module model_defn
CALL ASSIGN_PAR()        ! parameter definitions are stored in module multiparam

! Compute derived model parameters (bucket sizes, etc.)
CALL PAR_DERIVE(ERR,MESSAGE)
IF (ERR.NE.0) WRITE(*,*) TRIM(MESSAGE); IF (ERR.GT.0) STOP

! Define output and parameter files
ONEMOD=1                 ! one file per model (i.e., model dimension = 1)
PCOUNT=0                 ! counter for parameter sets evaluated (shared in MODULE multistats)

IF(fuse_mode == 'run_def')THEN ! run FUSE with default parameter values

  ! files to which model run and parameter set will be saved
  FNAME_NETCDF_RUNS = TRIM(OUTPUT_PATH)//TRIM(dom_id)//'_'//TRIM(FMODEL_ID)//'_runs_def.nc'
  FNAME_NETCDF_PARA = TRIM(OUTPUT_PATH)//TRIM(dom_id)//'_'//TRIM(FMODEL_ID)//'_para_def.nc'

ELSE IF(fuse_mode == 'run_pre')THEN  ! run FUSE with pre-defined parameter values

  ! files to which model run and parameter set will be saved
  FNAME_NETCDF_RUNS = TRIM(OUTPUT_PATH)//TRIM(dom_id)//'_'//TRIM(FMODEL_ID)//'_runs_pre.nc'
  FNAME_NETCDF_PARA = TRIM(OUTPUT_PATH)//TRIM(dom_id)//'_'//TRIM(FMODEL_ID)//'_para_pre_out.nc'

ELSE IF(fuse_mode == 'calib_sce')THEN ! calibrate FUSE using SCE

  ! files to which model run and parameter set will be saved
  FNAME_NETCDF_RUNS = TRIM(OUTPUT_PATH)//TRIM(dom_id)//'_'//TRIM(FMODEL_ID)//'_runs_sce.nc'
  FNAME_NETCDF_PARA = TRIM(OUTPUT_PATH)//TRIM(dom_id)//'_'//TRIM(FMODEL_ID)//'_para_sce.nc'

ELSE IF(fuse_mode == 'run_best')THEN  ! run FUSE with best (highest RMSE) parameter set from a previous SCE calibration

  ! file from which SCE parameters will be loaded - same as FNAME_NETCDF_PARA above
  FNAME_NETCDF_PARA_SCE = TRIM(OUTPUT_PATH)//TRIM(dom_id)//'_'//TRIM(FMODEL_ID)//'_para_sce.nc'

  ! files to which "best" SCE model run and parameter set will be saved
  FNAME_NETCDF_RUNS = TRIM(OUTPUT_PATH)//TRIM(dom_id)//'_'//TRIM(FMODEL_ID)//'_runs_best.nc'
  FNAME_NETCDF_PARA = TRIM(OUTPUT_PATH)//TRIM(dom_id)//'_'//TRIM(FMODEL_ID)//'_para_best.nc'

ELSE

  print *, 'Unexpected fuse_mode!'

ENDIF

CALL DEF_PARAMS(ONEMOD)  ! define model parameters (initial CREATE)
CALL DEF_SSTATS()        ! define summary statistics (REDEF)
CALL DEF_OUTPUT(numtim_sim,nSpat1,nSpat2)    ! define model time series (REDEF)

! ---------------------------------------------------------------------------------------
! RUN FUSE IN DESIRED MODE
! ---------------------------------------------------------------------------------------

! get parameter bounds and random numbers
ALLOCATE(APAR(NUMPAR),BL(NUMPAR),BU(NUMPAR),URAND(NUMPAR))

print *, 'NUMPAR = ', NUMPAR

print *, 'Using default parameter values:'

DO IPAR=1,NUMPAR
 CALL GETPAR_STR(LPARAM(IPAR)%PARNAME,PARAM_META)
 BL(IPAR)   = PARAM_META%PARLOW
 BU(IPAR)   = PARAM_META%PARUPP
 APAR(IPAR) = PARAM_META%PARDEF ! using default parameter values
 if(PARAM_META%PARFIT) print*, LPARAM(IPAR)%PARNAME, PARAM_META%PARDEF
END DO

IF(fuse_mode == 'run_def')THEN ! run FUSE with default parameter values

  OUTPUT_FLAG=.TRUE.

  print *, 'Running FUSE with default parameter values'
  CALL FUSE_RMSE(APAR,GRID_FLAG,NCID_FORC,RMSE,OUTPUT_FLAG)
  print *, 'Done running FUSE with default parameter values'

ELSE IF(fuse_mode == 'run_pre')THEN ! run FUSE with pre-defined parameter values

  OUTPUT_FLAG=.TRUE.

  OPEN(21,FILE=TRIM(file_para_list))
    DO   ! loop through parameter files

      ! get output filename
      !READ(21,*,IOSTAT=IERR) FNAME_NETCDF_PARA_PRE
      !IF (IERR.NE.0) EXIT

      ! file from which parameters will be loaded
      !FNAME_NETCDF_PARA_PRE = TRIM(OUTPUT_PATH)//TRIM(FNAME_NETCDF_PARA_PRE)
      !PRINT *, 'Loading parameter set from', $FNAME_NETCDF_PARA_PRE

      !CALL GET_PRE_PARAM(FNAME_NETCDF_PARA_PRE,1,ONEMOD,NUMPAR,APAR) ! load specific parameter set

      !print *, 'Running FUSE with pre-defined parameter set'
      !CALL FUSE_RMSE(APAR,GRID_FLAG,NCID_FORC,RMSE,OUTPUT_FLAG)
      !print *, 'Done running FUSE with pre-defined parameter set'

    END DO ! (looping through output files)

  CLOSE(21)

ELSE IF(fuse_mode == 'calib_sce')THEN ! calibrate FUSE using SCE

  ! Calibrate FUSE with SCE
  OUTPUT_FLAG=.FALSE.

  ! assign algorithmic control parameters for SCE
  NOPT   =  NUMPAR         ! number of parameters to be optimized (NUMPAR in module multiparam)
  MAXN   =     10000 			 ! maximum number of trials before optimization is terminated
  KSTOP  =      3          ! number of shuffling loops the value must change by PCENTO (MAX=9)
  PCENTO =      0.001      ! the percentage
  NGS    =     10          ! number of complexes in the initial population
  NPG    =  2*NOPT + 1     ! number of points in each complex
  NPS    =    NOPT + 1     ! number of points in a sub-complex
  NSPL   =  2*NOPT + 1     ! number of evolution steps allowed for each complex before shuffling
  MINGS  =  NGS            ! minimum number of complexes required
  INIFLG =  1              ! 1 = include initial point in the population
  IPRINT =  1              ! 0 = supress printing
  FNAME_ASCII = TRIM(OUTPUT_PATH)//TRIM(dom_id)//'_'//TRIM(FMODEL_ID)//'_sce_output.txt'

  ! convert from SP used in FUSE to MSP used in SCE
  ALLOCATE(APAR_MSP(NUMPAR),BL_MSP(NUMPAR),BU_MSP(NUMPAR),URAND_MSP(NUMPAR))

  APAR_MSP=APAR
  BL_MSP=BL
  BU_MSP=BU
  URAND_MSP=URAND

  ! open up ASCII output file
  print *, 'Creating SCE output file:', FNAME_ASCII
  ISCE = 96; OPEN(ISCE,FILE=TRIM(FNAME_ASCII))

  ! optimize (returns A and AF)
  ! note that SCE requires the kind of APAR, BL, BU to be MSP
  CALL SCEUA(APAR_MSP,AF_MSP,BL_MSP,BU_MSP,NOPT,MAXN,KSTOP,PCENTO,ISEED,&
          NGS,NPG,NPS,NSPL,MINGS,INIFLG,IPRINT,ISCE)

  ! close ASCII output file
  CLOSE(ISCE)

  PRINT *, 'Done running SCE!'

  ! call the function again with the optimized parameter set (to ensure the last parameter set is the optimum)
  !AF_MSP = FUNCTN(NOPT,AF_MSP)

  !PRINT *, 'Done calling the function again with the optimized parameter set!'

ELSE IF(fuse_mode == 'run_best')THEN ! run FUSE with best (highest RMSE) parameter set from a previous SCE calibration

  OUTPUT_FLAG=.TRUE.

  ! load best SCE parameter set from NetCDF file into APAR
  CALL GET_SCE_PARAM(FNAME_NETCDF_PARA_SCE,ONEMOD,NUMPAR,APAR)

  print *, 'Running FUSE with best SCE parameter set'
  CALL FUSE_RMSE(APAR,GRID_FLAG,NCID_FORC,RMSE,OUTPUT_FLAG)
  print *, 'Done running FUSE with best SCE parameter set'

ELSE

print *, 'Unexpected fuse_mode!'
stop

ENDIF

! deallocate space
DEALLOCATE(APAR,BL,BU,URAND)
IF(SPATIAL_OPTION == 'CATCH')THEN
  DEALLOCATE(aForce,aRoute,aValid)
  if(err/=0)then; write(*,*) 'unable to deallocate space for catchment modeling'; stop; endif

ELSE
  DEALLOCATE(gForce, gState)
  !DEALLOCATE(ancilF_3d, gForce_3d, gState_3d,AROUTE_3d)
  DEALLOCATE(gForce_3d, gState_3d,AROUTE_3d)
  if(err/=0)then; write(*,*) 'unable to deallocate space for grid modeling'; stop; endif

ENDIF

! close NetCDF files
IF(GRID_FLAG)THEN
  PRINT *, 'Closing forcing file'
  err = nf90_close(ncid_forc)
  if (err.ne.0) write(*,*) trim(message); if (err.gt.0) stop
ENDIF

PRINT *, 'Closing output file'
err = nf90_close(ncid_out)
if (err.ne.0) write(*,*) trim(message); if (err.gt.0) stop
PRINT *, 'Done'

STOP
END PROGRAM DISTRIBUTED_DRIVER
