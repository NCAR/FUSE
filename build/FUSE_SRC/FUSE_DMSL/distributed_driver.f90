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
          SETNGS_PATH,MBANDS_INFO, &
          OUTPUT_PATH,FORCINGINFO,INPUT_PATH

! data modules
USE model_defn,nstateFUSE=>nstate                         ! model definition structures
USE model_defnames                                        ! defines the integer model options
USE multiforce, ONLY: forcefile,vname_aprecip             ! model forcing structures
USE multiforce, ONLY: AFORCE, aValid                      ! time series of lumped forcing/response data
USE multiforce, ONLY: GFORCE, nspat1, nspat2              ! spatial array of gridded forcing data
USE multiforce, only: ancilF                              ! ancillary forcing data
USE multiforce, ONLY: valDat                              ! response data
USE multiforce, only: DELTIM
USE multiforce, only: ISTART, NUMTIM                      ! index for start of inference, number of data steps
USE multiforce, only: ncid_forc                           ! NetCDF forcing file ID
USE multiforce, only: ncid_var                            ! NetCDF forcing variable ID

!USE multiforce, ONLY: AFORCE, DELTIM, NUMTIM              ! data interval = maximum model time step - now redundant
USE multibands                                            ! basin band stuctures
USE multiparam, ONLY: LPARAM, PARATT, NUMPAR              ! parameter metadata structures
USE multistate, only: gState                              ! gridded state variables
USE multiroute, ONLY: AROUTE                              ! model routing structures
USE multistats                                            ! model statistics structures

! informational modules
USE selectmodl_module                                     ! reads model control file
USE getpar_str_module                                     ! extracts parameter metadata
USE par_insert_module                                     ! inserts model parameters
USE force_info_module,only:force_info                     ! get forcing info for NetCDF files
USE get_gforce_module,only:read_ginfo                     ! get dimension lengths from the NetCDF file
USE get_gforce_module,only:get_varid                      ! get netCDF ID for forcing variables
USE getf_ascii_module,only:prelim_asc                     ! get preliminary data from the ASCII file
USE getf_ascii_module,only:close_file                     ! close ASCII file
USE getf_ascii_module,only:read_ascii                     ! read ascii forcing data for a given time step

! model numerix
USE model_numerix                                         ! defines decisions on model numerix

! access to model simulation modules
USE fuse_rmse_module                                      ! run model and compute the root mean squared error
IMPLICIT NONE
! ---------------------------------------------------------------------------------------
! (0) GET COMMAND-LINE ARGUMENTS...
! ---------------------------------------------------------------------------------------
CHARACTER(LEN=64)                      :: DatString          ! string defining forcing data
CHARACTER(LEN=6)                       :: FMODEL_ID='      ' ! integer defining FUSE model
CHARACTER(LEN=6)                       :: F_SPATIAL='      ' ! spatial option (0=lumped, 1= distributed)
!CHARACTER(LEN=12)                      :: MBASIN_ID='      ' ! MOPEX basin ID
!CHARACTER(LEN=6)                       :: FMODEL_ID='      ' ! integer defining FUSE model
!CHARACTER(LEN=6)                       :: NSOLUTION='      ' ! numerical solution (0=explicit Euler; 1=explicit Heun; 2=implicit Euler; 3=implicit Heun, 4=semi-implicit)
!CHARACTER(LEN=6)                       :: FADAPTIVE='      ' ! identifier for adaptive sub-steps (0=fixed, 1=adaptive)
!CHARACTER(LEN=6)                       :: TRUNC_ABS='      ' ! absolute temporal truncation error tolerance
!CHARACTER(LEN=6)                       :: TRUNC_REL='      ' ! relative temporal truncation error tolerance
!CHARACTER(LEN=12)                      :: TSTEP_LEN='            ' ! maximum length of the time step (days)
!CHARACTER(LEN=6)                       :: NUMPARSET='      ' ! number of parameter sets
! ---------------------------------------------------------------------------------------
! (1) SETUP MODELS FOR SIMULATION -- POPULATE DATA STRUCTURES
! ---------------------------------------------------------------------------------------
! fuse_file_manager
CHARACTER(LEN=1024)                    :: FFMFILE      	  ! name of fuse_file_manager file - still needed?
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
! ---------------------------------------------------------------------------------------
! (2) RUN MODEL FOR DIFFERENT PARAMETER SETS
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
! (0) READ COMMAND LINE ARGUMENTS
! ---------------------------------------------------------------------------------------
! read command-line arguments
CALL GETARG(1,DatString)  ! string defining forcinginfo file
CALL GETARG(2,FMODEL_ID)  ! integer defining FUSE model
CALL GETARG(3,F_SPATIAL)  ! spatial option (0=lumped, 1= distributed)
!CALL GETARG(1,MBASIN_ID)  ! MOPEX basin ID
!CALL GETARG(2,FMODEL_ID)  ! integer defining FUSE model
!CALL GETARG(3,NSOLUTION)  ! numerical solution (0=explicit, 1=implicit)
!CALL GETARG(4,FADAPTIVE)  ! identifier for adaptive sub-steps (0=fixed, 1=adaptive)
!CALL GETARG(5,TRUNC_ABS)  ! absolute temporal truncation error tolerance
!CALL GETARG(6,TRUNC_REL)  ! relative temporal truncation error tolerance
!CALL GETARG(7,TSTEP_LEN)  ! maximum length of the time step (days)
!CALL GETARG(8,NUMPARSET)  ! number of parameter sets
! check command-line arguments
IF (LEN_TRIM(DatString).EQ.0) STOP '1st command-line argument is missing (DatString)'
IF (LEN_TRIM(FMODEL_ID).EQ.0) STOP '2nd command-line argument is missing (FMODEL_ID)'
IF (LEN_TRIM(F_SPATIAL).EQ.0) STOP '3rd command-line argument is missing (F_SPATIAL)'

!IF (LEN_TRIM(MBASIN_ID).EQ.0) STOP '1st command-line argument is missing (MBASIN_ID)'
!IF (LEN_TRIM(FMODEL_ID).EQ.0) STOP '2nd command-line argument is missing (FMODEL_ID)'
!IF (LEN_TRIM(NSOLUTION).EQ.0) STOP '3rd command-line argument is missing (NSOLUTION)'
!IF (LEN_TRIM(FADAPTIVE).EQ.0) STOP '4th command-line argument is missing (FADAPTIVE)'
!IF (LEN_TRIM(TRUNC_ABS).EQ.0) STOP '5th command-line argument is missing (TRUNC_ABS)'
!IF (LEN_TRIM(TRUNC_REL).EQ.0) STOP '6th command-line argument is missing (TRUNC_REL)'
!IF (LEN_TRIM(TSTEP_LEN).EQ.0) STOP '7th command-line argument is missing (TSTEP_LEN)'
!IF (LEN_TRIM(NUMPARSET).EQ.0) STOP '8th command-line argument is missing (NUMPARSET)'

! print command-line arguments
print*, '1st command-line argument (DatString) = ', trim(DatString)
print*, '2nd command-line argument (FMODEL_ID) = ', FMODEL_ID
print*, '3rd command-line argument (F_SPATIAL) = ', F_SPATIAL

! set path to fuse_file_manager - TODO check which FileManager is used
FFMFILE=TRIM(SETNGS_PATH)//TRIM(DatString)//'_fuse_file_manager.txt' ! still needed?
print *, 'fuse_file_manager:', TRIM(FFMFILE) ! still needed?

! get directories and filenames for control files
call fuse_SetDirsUndPhiles(fuseFileManagerIn=FFMFILE,err=err,message=message)
! call fuse_SetDirsUndPhiles(err=err,message=message) ! TODO: check which FileManager is used

if (err.ne.0) write(*,*) trim(message); if (err.gt.0) stop
CALL GETNUMERIX(ERR,MESSAGE)              ! defines method/parameters used for numerical solution
! define basin desired
FORCINGINFO = TRIM(DatString)//'_input_info.txt'
MBANDS_INFO = TRIM(DatString)//'_elev_bands_info.txt'

! convert command-line arguments to integer flags and real numbers
READ(FMODEL_ID,*) FUSE_ID                 ! integer defining FUSE model
READ(F_SPATIAL,*) SPATIAL_OPTION          ! spatial option (0=lumped, 1= distributed)
!READ(NSOLUTION,*) SOLUTION_METHOD         ! numerical solution (0=implicit, 1=explicit)
!READ(FADAPTIVE,*) TEMPORAL_ERROR_CONTROL  ! identifier for adaptive sub-steps (0=fixed, 1=adaptive)
!READ(TRUNC_ABS,*) ERR_TRUNC_ABS           ! absolute temporal truncation error tolerance
!READ(TRUNC_REL,*) ERR_TRUNC_REL           ! relative temporal truncation error tolerance
!READ(TSTEP_LEN,*) MAX_TSTEP               ! maximum length of the time step (days)
!READ(NUMPARSET,*) NUMPSET                 ! number of parameter sets

! ---------------------------------------------------------------------------------------
! (1) GET MODEL SETUP -- MODEL DEFINITION, AND PARAMETER AND VARIABLE INFO FOR ALL MODELS
! ---------------------------------------------------------------------------------------
! define the spatial flag (.true. is distributed)
SPATIAL_FLAG=.TRUE.; IF(SPATIAL_OPTION == LUMPED) SPATIAL_FLAG=.FALSE.
! get forcing info from the txt file
call force_info(err,message)
if(err/=0)then; write(*,*) trim(message); stop; endif
! allocate space for the basin-average time series
allocate(aForce(numtim),aRoute(numtim),aValid(numtim),stat=err)
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
 do iTim=1,numtim
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
 print*, istart,numtim
ELSE
 print *, 'Running FUSE as a distributed model'
 ! open NetCDF forcing file
 err = nf90_open(trim(INPUT_PATH)//trim(forcefile), nf90_nowrite, ncid_forc)
 if (err.ne.0) write(*,*) trim(message); if (err.gt.0) stop
 PRINT *, 'NCID is', ncid_forc

 ! get the grid info (dimensions) from the NetCDF file
 call read_ginfo(ncid_forc,err,message)
 if(err/=0)then; write(*,*) trim(message); stop; endif

 ! allocate space for the forcing grid and states
 allocate(ancilF(nspat1,nspat2), gForce(nspat1,nspat2), gState(nspat1,nspat2), stat=err)
 if(err/=0)then; write(*,*) 'unable to allocate space for forcing grid GFORCE'; stop; endif

 ! get variable ID from the NetCDF file
 call get_varID(ncid_forc,err,message)
 if(err/=0)then; write(*,*) 'unable to get NetCDF variables ID'; stop; endif

ENDIF

print*, 'spatial dimensions = ', nSpat1, nSpat2
print*, 'indices for start of inference and number of time steps', istart, numtim
print*, 'netCDF ID for forcing file', ncid_forc
print *, 'netCDF ID for second forcing variable:', ncid_var(2)

! Define model attributes (valid for all models)
CALL UNIQUEMODL(NMOD)           ! get nmod unique models
CALL GETPARMETA(ERR,MESSAGE)    ! read parameter metadata (parameter bounds etc.)
IF (ERR.NE.0) WRITE(*,*) TRIM(MESSAGE); IF (ERR.GT.0) STOP

! Identify a single model (read control file ../fuse_zDecisions.txt)
CALL SELECTMODL(ERR=ERR,MESSAGE=MESSAGE)
IF (ERR.NE.0) WRITE(*,*) TRIM(MESSAGE); IF (ERR.GT.0) STOP

! Define list of states and parameters for the current model
! Read data from the "BATEA-compliant" ASCII files
! CALL GETFORCING(INFERN_START,NTIM) ! read forcing data - not needed anymore
IF (SMODL%iSNOWM.EQ.iopt_temp_index) CALL GET_MBANDS(err,message) ! read band data if snow model
CALL ASSIGN_STT()        ! state definitions are stored in module model_defn
CALL ASSIGN_FLX()        ! flux definitions are stored in module model_defn
CALL ASSIGN_PAR()        ! parameter definitions are stored in module multiparam

! compute derived model parameters (bucket sizes, etc.)
CALL PAR_DERIVE(ERR,MESSAGE)
IF (ERR.NE.0) WRITE(*,*) TRIM(MESSAGE); IF (ERR.GT.0) STOP

! Define output file names (shared in MODULE model_defn)
FNAME_NETCDF = TRIM(OUTPUT_PATH)//TRIM(DatString)//'__'//TRIM(SMODL%MNAME)//'.nc'
write(*,'(a)') trim(fname_netcdf)

! Define NetCDF output files (only write parameters and summary statistics)
ONEMOD=1                 ! one file per model (i.e., model dimension = 1)
PCOUNT=0                 ! counter for parameter sets evaluated (shared in MODULE multistats)
OUTPUT_FLAG = .TRUE.    ! .TRUE. if desire time series output
CALL DEF_PARAMS(ONEMOD)  ! define model parameters (initial CREATE)
CALL DEF_SSTATS()        ! define summary statistics (REDEF)

IF (OUTPUT_FLAG) CALL DEF_OUTPUT(NUMTIM,nSpat1,nSpat2)    ! define model time series (REDEF)

! ---------------------------------------------------------------------------------------
! (2) RUN MODEL FOR THE CURRENT PARAMETER SET WITH DIFFERENT NUMERIX OPTIONS
! ---------------------------------------------------------------------------------------
! get parameter bounds and random numbers
ALLOCATE(APAR(NUMPAR),BL(NUMPAR),BU(NUMPAR),URAND(NUMPAR))
DO IPAR=1,NUMPAR
 CALL GETPAR_STR(LPARAM(IPAR)%PARNAME,PARAM_META)
 BL(IPAR)   = PARAM_META%PARLOW
 BU(IPAR)   = PARAM_META%PARUPP
 APAR(IPAR) = PARAM_META%PARDEF
 if(PARAM_META%PARFIT) print*, LPARAM(IPAR)%PARNAME, PARAM_META%PARDEF
END DO

! run zee model
print *, 'Entering FUSE_RMSE'
CALL FUSE_RMSE(APAR,SPATIAL_FLAG,NCID_FORC,RMSE,OUTPUT_FLAG)
print *, 'Done with FUSE_RMSE'

! and, deallocate space
DEALLOCATE(APAR,BL,BU,URAND)

! close the NetCDF file
IF(SPATIAL_OPTION /= LUMPED)THEN
  PRINT *, 'Closing forcing file'
  err = nf90_close(ncid_forc)
  if (err.ne.0) write(*,*) trim(message); if (err.gt.0) stop
ENDIF

STOP
END PROGRAM DISTRIBUTED_DRIVER
