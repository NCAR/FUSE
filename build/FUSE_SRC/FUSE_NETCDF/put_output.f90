SUBROUTINE PUT_OUTPUT(iSpat1,iSpat2,ITIM,IMOD,IPAR)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! write NetCDF output files
  ! ---------------------------------------------------------------------------------------
  USE nrtype                                            ! variable types, etc.
  USE model_defn                                        ! model definition (includes filename)
  USE metaoutput                                        ! metadata for time-varying model output
  USE varextract_module                                 ! interface for the function to extract variables
  USE multiforce,ONLY:timDat                            ! time data
  USE multistate, only: ncid_out                        ! NetCDF output file ID

  IMPLICIT NONE
  ! input
  INTEGER(I4B), INTENT(IN)               :: iSpat1      ! index of 1st spatial dimension
  INTEGER(I4B), INTENT(IN)               :: iSpat2      ! index of 2nd spatial dimension
  INTEGER(I4B), INTENT(IN)               :: ITIM        ! time step index
  INTEGER(I4B), INTENT(IN)               :: IMOD        ! model index
  INTEGER(I4B), INTENT(IN)               :: IPAR        ! parameter set index
  ! internal
  LOGICAL(LGT)                           :: WRITE_VAR   ! used to denote if the variable is written
  INTEGER(I4B)                           :: IERR        ! error code
  !INTEGER(I4B), DIMENSION(5)             :: INDX        ! indices for time series write
  INTEGER(I4B), DIMENSION(3)             :: INDX        ! indices for time series write
  INTEGER(I4B)                           :: IVAR        ! loop through variables
  REAL(SP)                               :: XVAR        ! desired variable (SP NOT NECESSARILY SP)
  REAL(MSP)                              :: AVAR        ! desired variable (SINGLE PRECISION)
  REAL(MSP)                              :: tDat        ! time data
  INTEGER(I4B)                           :: IVAR_ID     ! variable ID
  INCLUDE 'netcdf.inc'                                  ! use netCDF libraries
  ! ---------------------------------------------------------------------------------------
  ! open file
  IERR = NF_OPEN(TRIM(FNAME_NETCDF),NF_WRITE,ncid_out); CALL HANDLE_ERR(IERR)

  ! define indices for model output
  INDX = (/iSpat1,iSpat2,ITIM/)

  ! loop through time-varying model output
  DO IVAR=1,NOUTVAR

     ! check if there is a need to write the variable - see also def_output
     IF (Q_ONLY) THEN
        WRITE_VAR=.FALSE.
        IF (TRIM(VNAME(IVAR)).EQ.'ppt') WRITE_VAR=.TRUE.
        IF (TRIM(VNAME(IVAR)).EQ.'pet') WRITE_VAR=.TRUE.
        IF (TRIM(VNAME(IVAR)).EQ.'evap_1') WRITE_VAR=.TRUE.
        IF (TRIM(VNAME(IVAR)).EQ.'evap_2') WRITE_VAR=.TRUE.
        IF (TRIM(VNAME(IVAR)).EQ.'q_instnt') WRITE_VAR=.TRUE.
        !IF (TRIM(VNAME(IVAR)).EQ.'q_routed') WRITE_VAR=.TRUE.
        IF (TRIM(VNAME(IVAR)).EQ.'watr_1')   WRITE_VAR=.TRUE.
        IF (TRIM(VNAME(IVAR)).EQ.'watr_2')   WRITE_VAR=.TRUE.
        IF (TRIM(VNAME(IVAR)).EQ.'swe_tot')   WRITE_VAR=.TRUE.
        IF (.NOT.WRITE_VAR) CYCLE
     ENDIF

     ! write the variable
     XVAR = VAREXTRACT(VNAME(IVAR)); AVAR=XVAR                                  ! get variable ivar
     IERR = NF_INQ_VARID(ncid_out,TRIM(VNAME(IVAR)),IVAR_ID); CALL HANDLE_ERR(IERR) ! get variable ID
     IERR = NF_PUT_VAR1_REAL(ncid_out,IVAR_ID,INDX,AVAR); CALL HANDLE_ERR(IERR)     ! write data

  END DO  ! (ivar)

  ! write the time
  tDat = timDat%dtime ! convert to actual single precision
  ierr = nf_inq_varid(ncid_out,'time',ivar_id); CALL handle_err(ierr)        ! get variable ID for time
  ierr = nf_put_var1_real(ncid_out,ivar_id,itim,tDat); CALL handle_err(ierr) ! write time variable

  ! close NetCDF file
  IERR = NF_CLOSE(ncid_out)

END SUBROUTINE PUT_OUTPUT

SUBROUTINE PUT_GOUTPUT_3D(istart_sim,istart_in,numtim)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Nans Addor, based on Martyn Clark's 2007 PUT_OUTPUT
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! write a 3D data structure to the NetCDF output file
  ! ---------------------------------------------------------------------------------------
  USE nrtype                                            ! variable types, etc.
  USE model_defn                                        ! model definition (includes filename)
  USE metaoutput                                        ! metadata for time-varying model output
  USE varextract_module                                 ! interface for the function to extract variables
  USE multiforce, ONLY: timDat,time_steps               ! time data
  USE multistate, only: ncid_out                        ! NetCDF output file ID
  USE multiforce, ONLY: nspat1,nspat2                   ! spatial dimensions
  USE multiforce, ONLY: gForce_3d                       ! test only
  USE multiforce, only: NUMTIM                          ! number of data steps

  IMPLICIT NONE
  ! input
  INTEGER(I4B), INTENT(IN)               :: istart_sim  ! index start time step relative to numtim_sim
  INTEGER(I4B), INTENT(IN)               :: istart_in   ! index start time step relative to numtim_in - for time dimension
  INTEGER(I4B), INTENT(IN)               :: numtim      ! number of time steps to write
  ! internal
  LOGICAL(LGT)                           :: WRITE_VAR   ! used to denote if the variable is written
  INTEGER(I4B)                           :: IERR        ! error code
  INTEGER(I4B), DIMENSION(3)             :: IND_START   ! start indices
  INTEGER(I4B), DIMENSION(3)             :: IND_COUNT   ! count indices
  INTEGER(I4B)                           :: IVAR        ! loop through variables
  REAL(SP)                               :: XVAR        ! desired variable (SP NOT NECESSARILY SP)
  REAL(MSP)                              :: AVAR        ! desired variable (SINGLE PRECISION)
  REAL(SP), DIMENSION(nspat1,nspat2,numtim)    :: XVAR_3d        ! desired variable (SINGLE PRECISION)
  REAL(MSP), DIMENSION(nspat1,nspat2,numtim)   :: AVAR_3d        ! desired variable (SINGLE PRECISION)
  REAL(MSP), DIMENSION(:), ALLOCATABLE   :: tDat            ! time data
  REAL(SP), DIMENSION(:), ALLOCATABLE    :: time_steps_sub  ! time data
  INTEGER(I4B)                           :: IVAR_ID     ! variable ID
  INCLUDE 'netcdf.inc'                                  ! use netCDF libraries

  ! open file
  IERR = NF_OPEN(TRIM(FNAME_NETCDF),NF_WRITE,ncid_out); CALL HANDLE_ERR(IERR)

  ! define indices for model output
  IND_START = (/1,1,istart_sim/) ! The indices are relative to 1, i.e. the first data value of a variable would have index (1, 1, ..., 1)
  IND_COUNT = (/nspat1,nspat2,numtim/)

  ! loop through time-varying model output
  DO IVAR=1,NOUTVAR

    ! check if there is a need to write the variable - see also def_output
    IF (Q_ONLY) THEN
       WRITE_VAR=.FALSE.
       IF (TRIM(VNAME(IVAR)).EQ.'ppt') WRITE_VAR=.TRUE.
       IF (TRIM(VNAME(IVAR)).EQ.'pet') WRITE_VAR=.TRUE.
       IF (TRIM(VNAME(IVAR)).EQ.'evap_1') WRITE_VAR=.TRUE.
       IF (TRIM(VNAME(IVAR)).EQ.'evap_2') WRITE_VAR=.TRUE.
       IF (TRIM(VNAME(IVAR)).EQ.'q_instnt') WRITE_VAR=.TRUE.
       !IF (TRIM(VNAME(IVAR)).EQ.'q_routed') WRITE_VAR=.TRUE.
       IF (TRIM(VNAME(IVAR)).EQ.'watr_1')   WRITE_VAR=.TRUE.
       IF (TRIM(VNAME(IVAR)).EQ.'watr_2')   WRITE_VAR=.TRUE.
       IF (TRIM(VNAME(IVAR)).EQ.'swe_tot')   WRITE_VAR=.TRUE.
       !IF (TRIM(VNAME(IVAR)).EQ.'qsurf') WRITE_VAR=.TRUE.
       !IF (TRIM(VNAME(IVAR)).EQ.'oflow_1') WRITE_VAR=.TRUE.
       !IF (TRIM(VNAME(IVAR)).EQ.'qintf_1') WRITE_VAR=.TRUE.
       !IF (TRIM(VNAME(IVAR)).EQ.'oflow_2') WRITE_VAR=.TRUE.
       !IF (TRIM(VNAME(IVAR)).EQ.'qbase_2') WRITE_VAR=.TRUE.
       IF (.NOT.WRITE_VAR) CYCLE
    ENDIF

     ! write the variable
     XVAR_3d = VAREXTRACT_3d(VNAME(IVAR),numtim)   ! get variable
     AVAR_3d = XVAR_3d                             ! convert format

     IERR = NF_INQ_VARID(ncid_out,TRIM(VNAME(IVAR)),IVAR_ID); CALL HANDLE_ERR(IERR) ! get variable ID

     IERR = NF_PUT_VARA_REAL(ncid_out,IVAR_ID,IND_START,IND_COUNT,AVAR_3d); CALL HANDLE_ERR(IERR) ! write data

     IF(ierr/=0)THEN; PRINT*, TRIM('Problem while writing data to output file'); STOP; ENDIF

  END DO  ! (ivar)

  ! write the time
  allocate(tDat(numtim),time_steps_sub(numtim))

  time_steps_sub = time_steps(istart_in:(istart_in+numtim-1)) ! extract time for subperiod
  tDat = time_steps_sub ! convert to actual single precision
  ierr = nf_inq_varid(ncid_out,'time',ivar_id); CALL handle_err(ierr)             ! get variable ID for time
  ierr = nf_put_vara_real(ncid_out,ivar_id,istart_sim,numtim,tDat); CALL handle_err(ierr)  ! write time variable

  deallocate(tDat,time_steps_sub)

  ! close NetCDF file
  IERR = NF_CLOSE(ncid_out)
  IERR = NF_CLOSE(ncid_out)

END SUBROUTINE PUT_GOUTPUT_3D
