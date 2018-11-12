SUBROUTINE DEF_OUTPUT(nSpat1,nSpat2,NPSET,NTIM)

  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Define NetCDF output files -- time-varying model output
  ! ---------------------------------------------------------------------------------------

  USE nrtype                                            ! variable types, etc.
  USE model_defn                                        ! model definition (includes filename)
  USE metaoutput                                        ! metadata for all model variables
  USE fuse_fileManager,only: Q_ONLY                     ! only write streamflow to output file?
  USE multiforce, only: latitude,longitude              ! dimension arrays
  USE multiforce, only: name_psets,time_steps           ! dimension arrays
  USE multiforce, only: latUnits,lonUnits               ! units string
  USE multiforce, only: timeUnits                       ! units string
  USE multistate, only: ncid_out                        ! NetCDF output file ID

  IMPLICIT NONE

  ! input
  INTEGER(I4B), INTENT(IN)               :: NTIM           ! number of time steps
  INTEGER(I4B), INTENT(IN)               :: nSpat1,nSpat2  ! length of spatial dimensions
  INTEGER(I4B), INTENT(IN)               :: NPSET           ! number of parameter sets

  ! internal
  REAL(MSP),DIMENSION(nspat1)            :: longitude_msp        ! desired variable (SINGLE PRECISION)
  REAL(MSP),DIMENSION(nspat2)            :: latitude_msp         ! desired variable (SINGLE PRECISION)
  REAL(SP),parameter                     :: NA_VALUE_OUT= -9999. ! NA_VALUE for output file
  REAL(MSP)                              :: NA_VALUE_OUT_MSP     ! NA_VALUE for output file

  LOGICAL(LGT)                           :: WRITE_VAR   ! used to denote if the variable is written
  INTEGER(I4B)                           :: IERR        ! error code
  INTEGER(I4B)                           :: NTIM_DIM    ! time
  INTEGER(I4B)                           :: lon_dim     ! 1st spatial dimension
  INTEGER(I4B)                           :: lat_dim     ! 2nd spatial dimension
  INTEGER(I4B)                           :: param_dim   ! parameter set dimension
  INTEGER(I4B)                           :: NMOD_DIM    ! number of models
  INTEGER(I4B), DIMENSION(4)             :: TVAR        ! time-varying dimensions
  INTEGER(I4B)                           :: IVAR        ! loop through variables
  INTEGER(I4B)                           :: IVAR_ID     ! variable ID

  INTEGER(I4B)                           :: CHID          ! char position dimension id
  INTEGER(I4B),parameter                 :: TDIMS=2       ! char position dimension id
  INTEGER(I4B)                           :: TXDIMS(TDIMS) ! variable shape
  INTEGER(I4B)                           :: TSTART(TDIMS), TCOUNT(TDIMS)

  include 'netcdf.inc'                                  ! use netCDF libraries

  ! ---------------------------------------------------------------------------------------
  CALL VARDESCRIBE()  ! get list of variable descriptions
  ! ---------------------------------------------------------------------------------------
! put file in define mode
  print *, 'Create NetCDF file for runs:'
  PRINT *, FNAME_NETCDF_RUNS

  IERR = NF_CREATE(TRIM(FNAME_NETCDF_RUNS),NF_CLOBBER,ncid_out); CALL HANDLE_ERR(IERR)
  !IERR = NF_OPEN(TRIM(FNAME_NETCDF_RUNS),NF_WRITE,ncid_out); CALL HANDLE_ERR(IERR)
  !IERR = NF_REDEF(ncid_out); CALL HANDLE_ERR(IERR)

  ! define dimensions
  IERR = NF_DEF_DIM(ncid_out,'time',NF_UNLIMITED,NTIM_DIM); CALL HANDLE_ERR(IERR) !record dimension (unlimited length)
  IERR = NF_DEF_DIM(ncid_out,'longitude',nSpat1,lon_dim); CALL HANDLE_ERR(IERR)
  IERR = NF_DEF_DIM(ncid_out,'latitude',nSpat2,lat_dim); CALL HANDLE_ERR(IERR)
  IERR = NF_DEF_DIM(ncid_out,'param_set',NPSET,param_dim); CALL HANDLE_ERR(IERR)

  ! define character-position dimension for strings of max length 40
  !IERR = NF_DEF_DIM(ncid_out, "chid", 40, CHID); CALL HANDLE_ERR(IERR)

  ! define a character-string variable
  ! TXDIMS(1) = CHID   ! character-position dimension first
  ! TXDIMS(2) = NTIM_DIM ! record dimension ID
  ! IERR = NF_DEF_VAR(ncid_out, 'param_set',NF_CHAR, TDIMS, TXDIMS, param_dim); CALL HANDLE_ERR(IERR)

  ! retrieve ID for the model and parameter dimensions
  !IERR = NF_INQ_DIMID(ncid_out,'par',NPAR_DIM); CALL HANDLE_ERR(IERR)
  !IERR = NF_INQ_DIMID(ncid_out,'mod',NMOD_DIM); CALL HANDLE_ERR(IERR)

  ! assign dimensions to indices: for efficiency reasons, param_dim should be
  ! last, because it varies the slowest, but the NetCDF standard imposes
  ! the unlimited dimension to be last.
  TVAR = (/lon_dim,lat_dim,param_dim,NTIM_DIM/)

  ! define time-varying output variables
  DO IVAR=1,NOUTVAR

    ! check if there is a need to write the variable - see also put_output
    ! uncomment variables that should be written to output file
    IF (Q_ONLY) THEN
      WRITE_VAR=.FALSE.
      !IF (TRIM(VNAME(IVAR)).EQ.'ppt')      WRITE_VAR=.TRUE.
      !IF (TRIM(VNAME(IVAR)).EQ.'pet')      WRITE_VAR=.TRUE.
      !IF (TRIM(VNAME(IVAR)).EQ.'obsq')     WRITE_VAR=.TRUE.
      IF (TRIM(VNAME(IVAR)).EQ.'evap_1')   WRITE_VAR=.TRUE.
      IF (TRIM(VNAME(IVAR)).EQ.'evap_2')   WRITE_VAR=.TRUE.
      IF (TRIM(VNAME(IVAR)).EQ.'q_instnt') WRITE_VAR=.TRUE.
      !IF (TRIM(VNAME(IVAR)).EQ.'q_routed') WRITE_VAR=.TRUE.
      IF (TRIM(VNAME(IVAR)).EQ.'watr_1')   WRITE_VAR=.TRUE.
      IF (TRIM(VNAME(IVAR)).EQ.'watr_2')   WRITE_VAR=.TRUE.
      IF (TRIM(VNAME(IVAR)).EQ.'swe_tot')  WRITE_VAR=.TRUE.
      !IF (TRIM(VNAME(IVAR)).EQ.'qsurf')   WRITE_VAR=.TRUE.
      !IF (TRIM(VNAME(IVAR)).EQ.'oflow_1') WRITE_VAR=.TRUE.
      !IF (TRIM(VNAME(IVAR)).EQ.'qintf_1') WRITE_VAR=.TRUE.
      !IF (TRIM(VNAME(IVAR)).EQ.'oflow_2') WRITE_VAR=.TRUE.
      !IF (TRIM(VNAME(IVAR)).EQ.'qbase_2') WRITE_VAR=.TRUE.
      IF (.NOT.WRITE_VAR) CYCLE ! start new iteration of do loop, i.e. skip writting variable
    ENDIF

    ! write the variable
    IERR = NF_DEF_VAR(ncid_out,TRIM(VNAME(IVAR)),NF_REAL,4,TVAR,IVAR_ID); CALL HANDLE_ERR(IERR)

    IERR = NF_PUT_ATT_TEXT(ncid_out,IVAR_ID,'long_name',LEN_TRIM(LNAME(IVAR)),TRIM(LNAME(IVAR)))
    CALL HANDLE_ERR(IERR)
    IERR = NF_PUT_ATT_TEXT(ncid_out,IVAR_ID,'units',LEN_TRIM(VUNIT(IVAR)),TRIM(VUNIT(IVAR)))
    CALL HANDLE_ERR(IERR)
    !IERR = NF_DEF_VAR_FILL(ncid_out,IVAR_ID,0,NA_VALUE) ! define _FillValue for NetCDF4 files only
    NA_VALUE_OUT_MSP=NA_VALUE_OUT
    IERR = NF_PUT_ATT_REAL(ncid_out,IVAR_ID,'_FillValue',NF_FLOAT,1,NA_VALUE_OUT_MSP)
    CALL HANDLE_ERR(IERR)

  END DO  ! ivar

  ! define the time variable
  ierr = nf_def_var(ncid_out,'time',nf_real,1,ntim_dim,ivar_id); call handle_err(ierr)
  ierr = nf_put_att_text(ncid_out,ivar_id,'units',len_trim(timeUnits),trim(timeUnits))
  call handle_err(ierr)

  ! define the latitude variable
  ierr = nf_def_var(ncid_out,'latitude',nf_real,1,lat_dim,ivar_id); call handle_err(ierr)
  ierr = nf_put_att_text(ncid_out,ivar_id,'units',8,'degreesN'); call handle_err(ierr)
  ierr = nf_put_att_text(ncid_out,ivar_id,'axis',1,'Y'); call handle_err(ierr)

  ! define the longitude variable
  ierr = nf_def_var(ncid_out,'longitude',nf_real,1,lon_dim,ivar_id); call handle_err(ierr)
  ierr = nf_put_att_text(ncid_out,ivar_id,'units',8,'degreesE'); call handle_err(ierr)
  ierr = nf_put_att_text(ncid_out,ivar_id,'axis',1,'X'); call handle_err(ierr)

  ! define the param_set variable
  ierr = nf_def_var(ncid_out,'param_set',nf_char,1,param_dim,ivar_id); call handle_err(ierr)
  ierr = nf_put_att_text(ncid_out,ivar_id,'units',1,'-'); call handle_err(ierr)

  ! end definitions
  IERR = NF_ENDDEF(ncid_out); call handle_err(ierr)

  !IERR = NF_OPEN(TRIM(FNAME_NETCDF),NF_WRITE,ncid_out); CALL HANDLE_ERR(IERR)
  latitude_msp=latitude ! convert to actual single precision
  IERR = NF_INQ_VARID(ncid_out,'latitude',IVAR_ID); CALL HANDLE_ERR(IERR) ! get variable ID
  IERR = NF_PUT_VARA_REAL(ncid_out,IVAR_ID,1,nspat2,latitude_msp); CALL HANDLE_ERR(IERR) ! write data

  longitude_msp=longitude ! convert to actual single precision
  IERR = NF_INQ_VARID(ncid_out,'longitude',IVAR_ID); CALL HANDLE_ERR(IERR) ! get variable ID
  IERR = NF_PUT_VARA_REAL(ncid_out,IVAR_ID,1,nspat1,longitude_msp); CALL HANDLE_ERR(IERR) ! write data

  !TSTART(1) = 1      ! start at beginning of variable
  !TSTART(2) = 1      ! record number to write
  !TCOUNT(1) = 20     ! number of chars to write
  !TCOUNT(2) = 1      ! only write one record

  !IERR = NF_INQ_VARID(ncid_out,'param_set',IVAR_ID); CALL HANDLE_ERR(IERR) ! get variable ID
  !IERR = NF_PUT_VARA_TEXT(ncid_out,IVAR_ID,1,NPSET,name_psets); CALL HANDLE_ERR(IERR) ! write data
  !IERR = NF_PUT_VARA_TEXT(ncid_out,IVAR_ID,TSTART,TCOUNT,name_psets); CALL HANDLE_ERR(IERR) ! write data

  PRINT *, 'NetCDF file for model runs defined with dimensions', nSpat1 , nSpat2, NPSET, NTIM

  IERR = NF_ENDDEF(ncid_out)
  IERR = NF_CLOSE(ncid_out)

! ---------------------------------------------------------------------------------------
END SUBROUTINE DEF_OUTPUT
