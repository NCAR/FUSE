SUBROUTINE DEF_OUTPUT(NTIM,nSpat1,nSpat2)
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
USE multiforce,only:latitude,longitude,time_steps     ! dimension arrays
USE multiforce,only:latUnits,lonUnits,timeUnits       ! units string for time
USE multiforce, only:timeUnits                        ! units string for time
USE multistate, only: ncid_out                        ! NetCDF output file ID

IMPLICIT NONE
! input
INTEGER(I4B), INTENT(IN)               :: NTIM        ! number of time steps
INTEGER(I4B), INTENT(IN)               :: nSpat1,nSpat2 ! length of spatial dimensions
! internal
REAL(MSP),DIMENSION(nspat1)            :: longitude_msp        ! desired variable (SINGLE PRECISION)
REAL(MSP),DIMENSION(nspat2)            :: latitude_msp         ! desired variable (SINGLE PRECISION)

LOGICAL(LGT)                           :: WRITE_VAR   ! used to denote if the variable is written
INTEGER(I4B)                           :: IERR        ! error code
INTEGER(I4B)                           :: NTIM_DIM    ! time
INTEGER(I4B)                           :: lon_dim    ! 1st spatial dimension
INTEGER(I4B)                           :: lat_dim    ! 2nd spatial dimension
INTEGER(I4B)                           :: NPAR_DIM    ! number of parameter sets
INTEGER(I4B)                           :: NMOD_DIM    ! number of models
!INTEGER(I4B), DIMENSION(5)             :: TVAR        ! time-varying dimensions
INTEGER(I4B), DIMENSION(3)             :: TVAR        ! time-varying dimensions
INTEGER(I4B)                           :: IVAR        ! loop through variables
INTEGER(I4B)                           :: IVAR_ID     ! variable ID
include 'netcdf.inc'                                  ! use netCDF libraries
! ---------------------------------------------------------------------------------------
CALL VARDESCRIBE()  ! get list of variable descriptions
! ---------------------------------------------------------------------------------------
! put file in define mode
IERR = NF_OPEN(TRIM(FNAME_NETCDF),NF_WRITE,ncid_out); CALL HANDLE_ERR(IERR)
print *, 'ncid_out for outputfile in def mode - start:',ncid_out

IERR = NF_REDEF(ncid_out); CALL HANDLE_ERR(IERR)
 ! define time dimension
 IERR = NF_DEF_DIM(ncid_out,'time',NTIM,NTIM_DIM); CALL HANDLE_ERR(IERR)
 ! define spatial dimensions
 IERR = NF_DEF_DIM(ncid_out,'longitude',nSpat1,lon_dim); CALL HANDLE_ERR(IERR)
 IERR = NF_DEF_DIM(ncid_out,'latitude',nSpat2,lat_dim); CALL HANDLE_ERR(IERR)
 ! retrieve ID for the model and parameter dimensions
 IERR = NF_INQ_DIMID(ncid_out,'par',NPAR_DIM); CALL HANDLE_ERR(IERR)
 IERR = NF_INQ_DIMID(ncid_out,'mod',NMOD_DIM); CALL HANDLE_ERR(IERR)
 ! assign dimensions to indices - note that this specific dimension order
 ! was selected to optimize access to data
 TVAR = (/lon_dim,lat_dim,NTIM_DIM/) ! dimensions for time-varying output

 ! define time-varying output variables
 DO IVAR=1,NOUTVAR
  ! check if there is a need to write the variable - see also put_output
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
  IERR = NF_DEF_VAR(ncid_out,TRIM(VNAME(IVAR)),NF_REAL,3,TVAR,IVAR_ID); CALL HANDLE_ERR(IERR)

  IERR = NF_PUT_ATT_TEXT(ncid_out,IVAR_ID,'long_name',LEN_TRIM(LNAME(IVAR)),TRIM(LNAME(IVAR)))
  CALL HANDLE_ERR(IERR)
  IERR = NF_PUT_ATT_TEXT(ncid_out,IVAR_ID,'units',LEN_TRIM(VUNIT(IVAR)),TRIM(VUNIT(IVAR)))
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
 print *, 'latitude IVAR_ID in def mode:',IVAR_ID

 ! define the longitude variable
 ierr = nf_def_var(ncid_out,'longitude',nf_real,1,lon_dim,ivar_id); call handle_err(ierr)
 ierr = nf_put_att_text(ncid_out,ivar_id,'units',8,'degreesE'); call handle_err(ierr)
 ierr = nf_put_att_text(ncid_out,ivar_id,'axis',1,'X'); call handle_err(ierr)
 print *, 'longitude IVAR_ID in def mode:',IVAR_ID

 print *, 'ncid_out for outputfile in def mode - end :',ncid_out

! end definitions
IERR = NF_ENDDEF(ncid_out); call handle_err(ierr)

!IERR = NF_OPEN(TRIM(FNAME_NETCDF),NF_WRITE,ncid_out); CALL HANDLE_ERR(IERR)

print *, 'ncid_out for outputfile in write mode:',ncid_out

print *, 'latitude', latitude
latitude_msp=latitude ! convert to actual single precision
IERR = NF_INQ_VARID(ncid_out,'latitude',IVAR_ID); CALL HANDLE_ERR(IERR) ! get variable ID
print *, IVAR_ID
IERR = NF_PUT_VARA_REAL(ncid_out,IVAR_ID,1,nspat2,latitude_msp); CALL HANDLE_ERR(IERR) ! write data

print *, 'longitude', longitude
longitude_msp=longitude ! convert to actual single precision
IERR = NF_INQ_VARID(ncid_out,'longitude',IVAR_ID); CALL HANDLE_ERR(IERR) ! get variable ID
print *, IVAR_ID
IERR = NF_PUT_VARA_REAL(ncid_out,IVAR_ID,1,nspat1,longitude_msp); CALL HANDLE_ERR(IERR) ! write data

PRINT *, 'NetCDF file defined with dimensions', nSpat1 , nSpat2, NTIM

IERR = NF_CLOSE(ncid_out) ! closing it seems to be necessary to write dimensions

! ---------------------------------------------------------------------------------------
END SUBROUTINE DEF_OUTPUT
