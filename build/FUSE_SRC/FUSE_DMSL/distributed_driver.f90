PROGRAM DISTRIBUTED_DRIVER
! ---------------------------------------------------------------------------------------
! Creators:
! Martyn Clark, 2011
! Modified by Brian Henn to include snow model, 6/2013
! Modified by Nans Addor to include distributed modeling, 9/2016
! ---------------------------------------------------------------------------------------
! Purpose:
! Driver program to run FUSE as a distributed model with a snow module
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
CHARACTER(LEN=6)                       :: FMODEL_ID='      ' ! integer defining FUSE model
CHARACTER(LEN=6)                       :: F_SPATIAL='      ' ! spatial option (0=lumped, 1= distributed)

! ---------------------------------------------------------------------------------------
! SETUP MODELS FOR SIMULATION -- POPULATE DATA STRUCTURES
! ---------------------------------------------------------------------------------------
! fuse_file_manager
CHARACTER(LEN=1024)                    :: FFMFILE      	  ! name of fuse_file_manager file - still needed?
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
INTEGER(I4B)                           :: SPATIAL_OPTION  ! spatial option (0=lumped, 1= distributed)
INTEGER(I4B),PARAMETER                 :: LUMPED=0        ! named variable for lumped simulations
INTEGER(I4B),PARAMETER                 :: DISTRIBUTED=1   ! named variable for distributed simulations
LOGICAL(LGT)                           :: SPATIAL_FLAG    ! spatial flag .true. if distributed
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
! READ COMMAND LINE ARGUMENTS
! ---------------------------------------------------------------------------------------
! read command-line arguments
CALL GETARG(1,DatString)  ! string defining forcinginfo file
CALL GETARG(2,FMODEL_ID)  ! integer defining FUSE model
CALL GETARG(3,F_SPATIAL)  ! spatial option (0=lumped, 1= distributed)

! check command-line arguments
IF (LEN_TRIM(DatString).EQ.0) STOP '1st command-line argument is missing (DatString)'
IF (LEN_TRIM(FMODEL_ID).EQ.0) STOP '2nd command-line argument is missing (FMODEL_ID)'
IF (LEN_TRIM(F_SPATIAL).EQ.0) STOP '3rd command-line argument is missing (F_SPATIAL)'

! print command-line arguments
print*, '1st command-line argument (DatString) = ', trim(DatString)
print*, '2nd command-line argument (FMODEL_ID) = ', FMODEL_ID
print*, '3rd command-line argument (F_SPATIAL) = ', F_SPATIAL

! set path to fuse_file_manager
FFMFILE=TRIM(SETNGS_PATH)//TRIM(DatString)//'_fuse_file_manager.txt'

! set directories and filenames for control files
call fuse_SetDirsUndPhiles(fuseFileManagerIn=FFMFILE,err=err,message=message)
if (err.ne.0) write(*,*) trim(message); if (err.gt.0) stop

! defines method/parameters used for numerical solution
CALL GETNUMERIX(ERR,MESSAGE)

! define basin desired - ?overwrite fuse_SetDirsUndPhiles?
FORCINGINFO = TRIM(DatString)//'_input_info.txt'
MBANDS_INFO = TRIM(DatString)//'_elev_bands_info.txt'
ELEV_BANDS_NC = TRIM(DatString)//'_'//MBANDS_NC

! convert command-line arguments to integer flags and real numbers
READ(FMODEL_ID,*) FUSE_ID                 ! integer defining FUSE model
READ(F_SPATIAL,*) SPATIAL_OPTION          ! spatial option (0=lumped, 1= distributed)

! ---------------------------------------------------------------------------------------
! GET MODEL SETUP -- MODEL DEFINITION, AND PARAMETER AND VARIABLE INFO FOR ALL MODELS
! ---------------------------------------------------------------------------------------

! define the spatial flag (.true. is distributed)
SPATIAL_FLAG=.TRUE.; IF(SPATIAL_OPTION == LUMPED) SPATIAL_FLAG=.FALSE.
! get forcing info from the txt file, including NA_VALUE

call force_info(err,message)
if(err/=0)then; write(*,*) trim(message); stop; endif

print*, 'Number of timesteps per subperiod (numtim_sub) = ', numtim_sub

! allocate space for the basin-average time series
allocate(aForce(numtim_sub),aRoute(numtim_sub),aValid(numtim_sub),stat=err)
if(err/=0)then; write(*,*) 'unable to allocate space for basin-average time series [aForce,aRoute]'; stop; endif

! get dimensions of the grid
IF(SPATIAL_OPTION == LUMPED)THEN
	print *, 'Running FUSE as a lumped model'
 ! specify as a 1x1 grid
 nspat1=1; nspat2=1
 ! allocate space for the forcing grid and states
 allocate(gForce(nspat1,nspat2), gState(nspat1,nspat2), stat=err)
 if(err/=0)then; write(*,*) 'unable to allocate space for forcing grid GFORCE'; stop; endif
 ! get the column indices from the ascii file (also read to the start of the file)
 call prelim_asc(err,message); if(err/=0)then; write(*,*) trim(message); stop; endif
 ! get the ASCII forcing data for a given time step
 do iTim=1,numtim_sim
  ! get the ASCII forcing data for a given time step
  call read_ascii(err,message)
  if(err/=0)then; print*, trim(message); stop; endif
  ! save forcing/response data in the structure
  aForce(iTim) = gForce(1,1)
  aValid(iTim) = valDat
 end do  ! itim
 ! close file
 call close_file(err,message)
 if(err/=0)then; print*,trim(message); read(*,*); endif
 print*, istart,numtim_sim
ELSE
 print *, 'Running FUSE as a distributed model'
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
 allocate(AROUTE_3d(nspat1,nspat2,numtim_sub+1), gState_3d(nspat1,nspat2,numtim_sub+1),gForce_3d(nspat1,nspat2,numtim_sub), stat=err)
 if(err/=0)then; write(*,*) 'unable to allocate space for 3d structure'; stop; endif

 ! allocate space for elevation bands
 allocate(MBANDS_VAR_4d(nspat1,nspat2,6,numtim_sub+1),stat=err) ! TODO: replace 6 by N_BANDS
 if(err/=0)then; write(*,*) 'unable to allocate space for elevation bands'; stop; endif

 ! get variable ID from the NetCDF file
 call get_varID(ncid_forc,err,message)
 if(err/=0)then; write(*,*) 'unable to get NetCDF variables ID'; stop; endif

ENDIF

print*, 'spatial dimensions = ', nSpat1, nSpat2
print*, 'netCDF ID for forcing file', ncid_forc
print*, 'NA_VALUE = ', NA_VALUE

! Define model attributes (valid for all models)
CALL UNIQUEMODL(NMOD)           ! get nmod unique models
CALL GETPARMETA(ERR,MESSAGE)    ! read parameter metadata (parameter bounds etc.)

IF (ERR.NE.0) WRITE(*,*) TRIM(MESSAGE); IF (ERR.GT.0) STOP

! Identify a single model (read control file ../fuse_zDecisions.txt)
CALL SELECTMODL(ERR=ERR,MESSAGE=MESSAGE)
IF (ERR.NE.0) WRITE(*,*) TRIM(MESSAGE); IF (ERR.GT.0) STOP

! Create elevations bands if snow module is on
IF (SMODL%iSNOWM.EQ.iopt_temp_index) THEN

  IF(SPATIAL_OPTION == LUMPED)THEN

   CALL GET_MBANDS(err,message) ! read band data from ASCII file - for one basin

  ELSE

   CALL GET_MBANDS_INFO(ELEV_BANDS_NC,err,message) ! read band data from NetCDF file - for a 2D grid

  ENDIF
ENDIF

! Define list of states and parameters for the current model
CALL ASSIGN_STT()        ! state definitions are stored in module model_defn
CALL ASSIGN_FLX()        ! flux definitions are stored in module model_defn
CALL ASSIGN_PAR()        ! parameter definitions are stored in module multiparam

! Compute derived model parameters (bucket sizes, etc.)
CALL PAR_DERIVE(ERR,MESSAGE)
IF (ERR.NE.0) WRITE(*,*) TRIM(MESSAGE); IF (ERR.GT.0) STOP

! Define output file and create it
FNAME_NETCDF = TRIM(OUTPUT_PATH)//TRIM(DatString)//'__'//TRIM(SMODL%MNAME)//'.nc'
ONEMOD=1                 ! one file per model (i.e., model dimension = 1)
PCOUNT=0                 ! counter for parameter sets evaluated (shared in MODULE multistats)
OUTPUT_FLAG = .TRUE.     ! .TRUE. if desire time series output
CALL DEF_PARAMS(ONEMOD)  ! define model parameters (initial CREATE)
CALL DEF_SSTATS()        ! define summary statistics (REDEF)
CALL DEF_OUTPUT(numtim_sim,nSpat1,nSpat2)    ! define model time series (REDEF)

! ---------------------------------------------------------------------------------------
! RUN MODEL FOR THE CURRENT PARAMETER SET
! ---------------------------------------------------------------------------------------

! get parameter bounds and random numbers
ALLOCATE(APAR(NUMPAR),BL(NUMPAR),BU(NUMPAR),URAND(NUMPAR))

print *, 'Using default parameter values:'

DO IPAR=1,NUMPAR
 CALL GETPAR_STR(LPARAM(IPAR)%PARNAME,PARAM_META)
 BL(IPAR)   = PARAM_META%PARLOW
 BU(IPAR)   = PARAM_META%PARUPP
 APAR(IPAR) = PARAM_META%PARDEF ! Using default parameter values
 if(PARAM_META%PARFIT) print*, LPARAM(IPAR)%PARNAME, PARAM_META%PARDEF
END DO

! run zee model
print *, 'Entering FUSE_RMSE'
!CALL FUSE_RMSE(APAR,SPATIAL_FLAG,NCID_FORC,RMSE,OUTPUT_FLAG)
CALL FUSE_RMSE(APAR,SPATIAL_FLAG,NCID_FORC,OUTPUT_FLAG)
print *, 'Done with FUSE_RMSE'

! deallocate space
DEALLOCATE(APAR,BL,BU,URAND)
IF(SPATIAL_OPTION == LUMPED)THEN
  DEALLOCATE(aForce,aRoute,aValid)
  if(err/=0)then; write(*,*) 'unable to deallocate space for lumped modeling'; stop; endif

ELSE
  DEALLOCATE(gForce, gState)
  !DEALLOCATE(ancilF_3d, gForce_3d, gState_3d,AROUTE_3d)
  DEALLOCATE(gForce_3d, gState_3d,AROUTE_3d)
  if(err/=0)then; write(*,*) 'unable to deallocate space for distributed modeling'; stop; endif

ENDIF

! close NetCDF files
IF(SPATIAL_OPTION /= LUMPED)THEN
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
