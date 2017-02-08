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
USE multiforce, only:timeUnits                        ! units string for time
USE multistate, only: ncid_out                        ! NetCDF output file ID

IMPLICIT NONE
! input
INTEGER(I4B), INTENT(IN)               :: NTIM        ! number of time steps
INTEGER(I4B), INTENT(IN)               :: nSpat1,nSpat2 ! length of spatial dimensions
! internal
LOGICAL(LGT)                           :: WRITE_VAR   ! used to denote if the variable is written
INTEGER(I4B)                           :: IERR        ! error code
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
! put file in define mode
IERR = NF_OPEN(TRIM(FNAME_NETCDF),NF_WRITE,ncid_out); CALL HANDLE_ERR(IERR)
IERR = NF_REDEF(ncid_out); CALL HANDLE_ERR(IERR)
 ! define time dimension
 IERR = NF_DEF_DIM(ncid_out,'time',NTIM,NTIM_DIM); CALL HANDLE_ERR(IERR)
 ! define spatial dimensions
 IERR = NF_DEF_DIM(ncid_out,'sp1',nSpat1,nSp1_DIM); CALL HANDLE_ERR(IERR)
 IERR = NF_DEF_DIM(ncid_out,'sp2',nSpat2,nSp2_DIM); CALL HANDLE_ERR(IERR)
 ! retrieve ID for the model and parameter dimensions
 IERR = NF_INQ_DIMID(ncid_out,'par',NPAR_DIM); CALL HANDLE_ERR(IERR)
 IERR = NF_INQ_DIMID(ncid_out,'mod',NMOD_DIM); CALL HANDLE_ERR(IERR)
 ! assign dimensions to indices - note that this specific dimension order
 ! was selected to optimize access to data
 !TVAR = (/NTIM_DIM,nSp2_DIM,nSp1_DIM,NMOD_DIM,NPAR_DIM/) ! dimensions for time-varying output
 TVAR = (/nSp1_DIM,nSp2_DIM,NTIM_DIM,NMOD_DIM,NPAR_DIM/) ! dimensions for time-varying output

 ! define time-varying output variables
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
  IERR = NF_DEF_VAR(ncid_out,TRIM(VNAME(IVAR)),NF_REAL,5,TVAR,IVAR_ID); CALL HANDLE_ERR(IERR)

  IERR = NF_PUT_ATT_TEXT(ncid_out,IVAR_ID,'long_name',LEN_TRIM(LNAME(IVAR)),TRIM(LNAME(IVAR)))
  CALL HANDLE_ERR(IERR)
  IERR = NF_PUT_ATT_TEXT(ncid_out,IVAR_ID,'units',LEN_TRIM(VUNIT(IVAR)),TRIM(VUNIT(IVAR)))
  CALL HANDLE_ERR(IERR)
 END DO  ! ivar
 ! define the time variable
 ierr = nf_def_var(ncid_out,'time',nf_real,1,ntim_dim,ivar_id); call handle_err(ierr)
 ierr = nf_put_att_text(ncid_out,ivar_id,'units',len_trim(timeUnits),trim(timeUnits))
 call handle_err(ierr)
! end definitions
IERR = NF_ENDDEF(ncid_out)

PRINT *, 'NCID_OUT is', ncid_out

!IERR = NF_CLOSE(ncid_out)
PRINT *, 'Keeping NetCDF output file open'

! ---------------------------------------------------------------------------------------
END SUBROUTINE DEF_OUTPUT
