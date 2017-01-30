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
USE multiforce,only:timeUnits                         ! units string for time
IMPLICIT NONE
! input
INTEGER(I4B), INTENT(IN)               :: NTIM        ! number of time steps
INTEGER(I4B), INTENT(IN)               :: nSpat1,nSpat2 ! length of spatial dimensions
! internal
LOGICAL(LGT)                           :: WRITE_VAR   ! used to denote if the variable is written
INTEGER(I4B)                           :: IERR,NCID   ! error code; NetCDF ID
INTEGER(I4B)                           :: NTIM_DIM    ! time
INTEGER(I4B)                           :: nSp1_DIM    ! 1st spatial dimension
INTEGER(I4B)                           :: nSp2_DIM    ! 2nd spatial dimension
INTEGER(I4B)                           :: NPAR_DIM    ! number of parameter sets
INTEGER(I4B)                           :: NMOD_DIM    ! number of models
INTEGER(I4B), DIMENSION(5)             :: TVAR        ! time-varying dimensions
INTEGER(I4B)                           :: IVAR        ! loop through variables
INTEGER(I4B)                           :: IVAR_ID     ! variable ID
include 'netcdf.inc'                                  ! use netCDF libraries
! ---------------------------------------------------------------------------------------
CALL VARDESCRIBE()  ! get list of variable descriptions
! ---------------------------------------------------------------------------------------
! open file and put in define mode
IERR = NF_OPEN(TRIM(FNAME_NETCDF),NF_WRITE,NCID); CALL HANDLE_ERR(IERR)
IERR = NF_REDEF(NCID); CALL HANDLE_ERR(IERR)
 ! define time dimension
 IERR = NF_DEF_DIM(NCID,'time',NTIM,NTIM_DIM); CALL HANDLE_ERR(IERR)
 ! define spatial dimensions
 IERR = NF_DEF_DIM(NCID,'sp1',nSpat1,nSp1_DIM); CALL HANDLE_ERR(IERR)
 IERR = NF_DEF_DIM(NCID,'sp2',nSpat2,nSp2_DIM); CALL HANDLE_ERR(IERR)
 ! retrieve ID for the model and parameter dimensions
 IERR = NF_INQ_DIMID(NCID,'par',NPAR_DIM); CALL HANDLE_ERR(IERR)
 IERR = NF_INQ_DIMID(NCID,'mod',NMOD_DIM); CALL HANDLE_ERR(IERR)
 ! assign dimensions to indices
 TVAR = (/NTIM_DIM,nSp2_DIM,nSp1_DIM,NMOD_DIM,NPAR_DIM/)   ! dimensions for time-varying output
 ! define time-varying output variables

 print *, 'Defining NetCDF output file'

 DO IVAR=1,NOUTVAR
  ! check if there is a need to write the variable
  IF (Q_ONLY) THEN
   WRITE_VAR=.FALSE.
   IF (TRIM(VNAME(IVAR)).EQ.'q_routed') WRITE_VAR=.TRUE.
   IF (TRIM(VNAME(IVAR)).EQ.'watr_1')   WRITE_VAR=.TRUE.
   IF (TRIM(VNAME(IVAR)).EQ.'watr_2')   WRITE_VAR=.TRUE.
   IF (.NOT.WRITE_VAR) CYCLE
  ENDIF
  ! write the variable
  IERR = NF_DEF_VAR(NCID,TRIM(VNAME(IVAR)),NF_REAL,5,TVAR,IVAR_ID); CALL HANDLE_ERR(IERR)
  IERR = NF_PUT_ATT_TEXT(NCID,IVAR_ID,'long_name',LEN_TRIM(LNAME(IVAR)),TRIM(LNAME(IVAR)))
  CALL HANDLE_ERR(IERR)
  IERR = NF_PUT_ATT_TEXT(NCID,IVAR_ID,'units',LEN_TRIM(VUNIT(IVAR)),TRIM(VUNIT(IVAR)))
  CALL HANDLE_ERR(IERR)
 END DO  ! ivar
 ! define the time variable
 ierr = nf_def_var(ncid,'time',nf_real,1,ntim_dim,ivar_id); call handle_err(ierr)
 ierr = nf_put_att_text(ncid,ivar_id,'units',len_trim(timeUnits),trim(timeUnits))
 call handle_err(ierr)
! end definitions and close file
IERR = NF_ENDDEF(NCID)
IERR = NF_CLOSE(NCID)
! ---------------------------------------------------------------------------------------
END SUBROUTINE DEF_OUTPUT
