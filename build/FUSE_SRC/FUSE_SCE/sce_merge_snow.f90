PROGRAM SCE_MERGE
! ---------------------------------------------------------------------------------------
! Creator:
! Martyn Clark, 2009
! Modified by Nans Addor to include snow module
! ---------------------------------------------------------------------------------------
! Purpose:
! Driver program to merge SCE runs from multiple models
! ---------------------------------------------------------------------------------------
USE nrtype                                                ! variable types, etc.
! USE ddirectory                                          ! directory for data files commented because couldn't be found
USE fuse_fileManager,only:fuse_SetDirsUndPhiles,&         ! sets directories and filenames - added this line because program looking for OUTPUT_PATH, that probably was in ddirectory
     OUTPUT_PATH,FORCINGINFO
! data modules
USE model_defn,nstateFUSE=>nstate                         ! model definition structures
USE multiparam, ONLY: LPARAM, PARATT, NUMPAR              ! parameter metadata structures
USE multistats, ONLY: PCOUNT, MOD_IX                      ! parameter set / model counters
USE MODEL_DEFNAMES, ONLY: iopt_temp_index                 ! indice corresponding to snow module
! access to model simulation modules
USE fuse_rmse_module                                      ! run model and compute the root mean squared error
IMPLICIT NONE
! ---------------------------------------------------------------------------------------
! (0) GET COMMAND-LINE ARGUMENTS
! ---------------------------------------------------------------------------------------
CHARACTER(LEN=1024)                    :: FFMFILE='        ' ! name of fuse_file_manager file
CHARACTER(LEN=12)                      :: MBASIN_ID='      ' ! MOPEX basin ID
CHARACTER(LEN=1024)                    :: NC_FILE='        ' ! path/name of muster file
! ---------------------------------------------------------------------------------------
! (1) PRELIMINARIES... GET DATA AND NUMERIX DECISIONS
! ---------------------------------------------------------------------------------------
! get model forcing data
INTEGER(I4B)                           :: NTIM            ! number of time steps
INTEGER(I4B)                           :: INFERN_START    ! start of inference period
! ---------------------------------------------------------------------------------------
! (2) READ LIST OF OUTPUT FILES, AND RUN MODEL FOR BEST PARAMETER SET IN EACH ONE
! ---------------------------------------------------------------------------------------
INTEGER(I4B)                           :: I,J,K           ! looping   
INTEGER(I4B)                           :: IERR            ! error code for reading input files
INTEGER(I4B)                           :: ERR             ! error code
CHARACTER(LEN=256)                     :: MESSAGE         ! error message
LOGICAL(LGT)                           :: LEXIST          ! .TRUE. if the file exists
INTEGER(I4B)                           :: NMODEL=1        ! number of models in the file list
CHARACTER(LEN=1024)                    :: FILE_NAME       ! name of single NetCDF output file
INTEGER(I4B)                           :: ONEMOD=1        ! just one model in output file
LOGICAL(LGT)                           :: OUTPUT_FLAG     ! switch off/on model output
INTEGER(I4B)                           :: MPAR            ! number of model parameters
REAL(SP), DIMENSION(:), ALLOCATABLE    :: XPAR            ! model parameters
REAL(SP)                               :: RMSE            ! root mean squared error
! ---------------------------------------------------------------------------------------
! (0) GET COMMAND-LINE ARGUMENTS
! ---------------------------------------------------------------------------------------
! pad FILE_LIST with blanks
!DO I=1,LEN(FILE_LIST); FILE_LIST(I:I)=' '; END DO
! read input filename from the command line
!CALL GETARG(1,FILE_LIST)
!IF (LEN_TRIM(FILE_LIST).EQ.0) STOP '1st command-line argument is missing (FILE_LIST)'
CALL GETARG(1,FFMFILE)    ! name of fuse_file_manager file
CALL GETARG(2,MBASIN_ID)  ! basin ID
CALL GETARG(3,NC_FILE)    ! file containing list of SCE output files

IF (LEN_TRIM(FFMFILE).EQ.0) STOP '1st command-line argument is missing (MUSTRFILE)'
IF (LEN_TRIM(MBASIN_ID).EQ.0) STOP '2nd command-line argument is missing (MBASIN_ID)'
IF (LEN_TRIM(NC_FILE).EQ.0) STOP '3rd command-line argument is missing (NC_FILE)'

! get directories and filenames for control files
call fuse_SetDirsUndPhiles(fuseFileManagerIn=trim(FFMFILE),err=err,message=message)
if (err.ne.0) write(*,*) trim(message); if (err.gt.0) stop
! define basin desired
!FORCINGINFO = 'forcinginfo.'//TRIM(MBASIN_ID)//'.txt'

! ---------------------------------------------------------------------------------------
! (1) PRELIMINARIES... GET DATA AND NUMERIX DECISIONS
! ---------------------------------------------------------------------------------------
! Define  method/parameters used for numerical solution
CALL GETNUMERIX(ERR,MESSAGE)         
PRINT *, 'Method/parameters used for numerical solution defined!'

! Read data from the "BATEA-compliant" ASCII files
CALL GETFORCING(INFERN_START,NTIM)
PRINT *, 'Forcing loaded!'

! Read parameter metadata (parameter bounds etc.) for all models
CALL GETPARMETA(ERR,MESSAGE)         
PRINT *, 'Got parameters!'

! ---------------------------------------------------------------------------------------
! (2) READ LIST OF OUTPUT FILES, AND RUN MODEL FOR BEST PARAMETER SET IN EACH ONE
! ---------------------------------------------------------------------------------------
! check that the file containing list of SCE output files exists
PRINT *, 'File containing list of SCE output files: ', TRIM(NC_FILE)

INQUIRE(FILE=TRIM(NC_FILE),EXIST=LEXIST)
IF (.NOT.LEXIST) STOP 'file containing list of SCE output files does not exist'
! get number of output files (models) to process
NMODEL = 0
OPEN(21,FILE=TRIM(NC_FILE))
 DO; READ(21,*,IOSTAT=IERR) FILE_NAME; IF (IERR.NE.0) EXIT; NMODEL=NMODEL+1; END DO
CLOSE(21)

! Define structure of NetCDF file to store the simulation using optimum param sets
FNAME_NETCDF = TRIM(OUTPUT_PATH)//'sce_modrun_'//FILE_NAME(12:LEN_TRIM(FILE_NAME)) ! WARNING MIGHT WANT TO USE NC_FILE INSTEAD 
print*, 'NetCDF file to store the simulation using optimum param sets:',FNAME_NETCDF
OUTPUT_FLAG = .TRUE.     ! .TRUE. if desire time series output
CALL DEF_PARAMS(NMODEL)  
print*, 'NetCDF parameters ready!'

CALL DEF_SSTATS()        ! define summary statistics (REDEF)
print*, 'NetCDF summary statistics ready!'
!IF (OUTPUT_FLAG) CALL DEF_OUTPUT(NTIM)    ! define model time series (REDEF) - moved below because SMODL not populated at this stage

print *, 'Loading model information from NetCDF files'
! initialize the model index (stared in module multistats)
MOD_IX = 0
! loop through output files
OPEN(21,FILE=TRIM(NC_FILE))
 DO   ! loop through output files
  ! get output filename
  READ(21,*,IOSTAT=IERR) FILE_NAME
  IF (IERR.NE.0) EXIT

  ! identify model (populate SMODL)
  CALL GET_SMODEL(FILE_NAME,ONEMOD)
  print *, 'Done populating SMODL!'

  ! read band data if snow model 
  IF (SMODL%iSNOWM.EQ.iopt_temp_index) CALL GET_MBANDS(err,message) 
  
  IF (OUTPUT_FLAG) CALL DEF_OUTPUT(NTIM)    ! define model time series (REDEF) - moved here
  
  ! Define list of states and parameters for the current model
  CALL ASSIGN_STT()        ! state definitions are stored in module model_defn
  CALL ASSIGN_FLX()        ! flux definitions are stored in module model_defn
  
  print *, 'Assigning parameters necessary for this model structure'
  CALL ASSIGN_PAR()        ! parameter definitions are stored in module multiparam
  ! get final parameter set
  MPAR = NUMPAR            ! (number of model parameter sets)

  ALLOCATE(XPAR(MPAR),STAT=IERR); IF (IERR.NE.0) STOP ' problem allocating XPAR '

  print *, 'Determining parameter set with lowest RSME:'
  CALL GET_FPARAM(FILE_NAME,ONEMOD,MPAR,XPAR)

  WRITE(*,'(20(A11,1X))') LPARAM(1:NUMPAR)
  WRITE(*,'(20(F11.3,1X))') XPAR(1:NUMPAR)
  ! compute derived model parameters (bucket sizes, etc.)

  CALL PAR_DERIVE(ERR,MESSAGE) ! added err,message (from URS_driver_sce.f90) because PAR_DERIVE requires arguments
  ! define indices for data write
  PCOUNT=0          ! ensure the parameter counter is set to zero (incremented in fuse_rmse)
  MOD_IX=MOD_IX + 1 ! increment the model index
  ! run zee model
  PRINT *, XPAR

  print *, 'Running model with best parameter set'
  CALL FUSE_RMSE(XPAR,RMSE,OUTPUT_FLAG,ERR=ERR,MESSAGE=MESSAGE)
  print *, 'Simulation completed'

  ! deallocate space
  DEALLOCATE(XPAR, STAT=IERR); IF (IERR.NE.0) STOP ' problem deallocating XPAR '
 END DO ! (looping through output files)
CLOSE(21)
STOP
END PROGRAM SCE_MERGE
