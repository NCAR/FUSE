SUBROUTINE DEF_SSTATS()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Define NetCDF output files -- summary statistics
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition (includes filename)
USE meta_stats                                        ! metadata for summary statistics
USE model_numerix                                     ! model numerix decisions
USE multistate, only: ncid_out                        ! NetCDF output file ID
IMPLICIT NONE
! internal
INTEGER(I4B)                           :: IERR        ! error code; NetCDF ID
INTEGER(I4B)                           :: NPAR_DIM    ! number of parameter sets
INTEGER(I4B)                           :: NMOD_DIM    ! number of models
!INTEGER(I4B)                           :: NORD_DIM    ! number of ordinates in prob distn
INTEGER(I4B), DIMENSION(1)             :: FVAR        ! dimensions for summary statistics
INTEGER(I4B), DIMENSION(2)             :: PVAR        ! dimensions for probability distributions
INTEGER(I4B)                           :: IVAR        ! loop through variables
INTEGER(I4B)                           :: IVAR_ID     ! variable ID
!INTEGER(I4B)                           :: IORD_ID     ! ordinates ID
!real(MSP), dimension(size(ORD_NSUBS))  :: rORD        ! ordinates of the prob dist (real numbers)
include 'netcdf.inc'                                  ! use netCDF libraries
! ---------------------------------------------------------------------------------------
CALL SUMDESCRIBE()  ! get list of summary statistics
! ---------------------------------------------------------------------------------------
! open file and put in define mode
IERR = NF_OPEN(TRIM(FNAME_NETCDF_PARA),NF_WRITE,ncid_out); CALL HANDLE_ERR(IERR)
IERR = NF_REDEF(ncid_out); CALL HANDLE_ERR(IERR)
 ! retrieve ID for the model and parameter dimensions
 IERR = NF_INQ_DIMID(ncid_out,'par',NPAR_DIM); CALL HANDLE_ERR(IERR)
 !IERR = NF_INQ_DIMID(ncid_out,'mod',NMOD_DIM); CALL HANDLE_ERR(IERR)

 ! define ord dimension
 !IERR = NF_DEF_DIM(ncid_out,'ord',SIZE(ORD_NSUBS),NORD_DIM); CALL HANDLE_ERR(IERR)

 ! define variables
 FVAR = (/NPAR_DIM/)            ! dimensions for fixed output (parameters)
 DO IVAR=1,NSUMVAR
  IERR = NF_DEF_VAR(ncid_out,TRIM(XNAME(IVAR)),NF_REAL,1,FVAR,IVAR_ID); CALL HANDLE_ERR(IERR)
  IERR = NF_PUT_ATT_TEXT(ncid_out,IVAR_ID,'long_name',LEN_TRIM(XDESC(IVAR)),TRIM(XDESC(IVAR)))

  CALL HANDLE_ERR(IERR)
  IERR = NF_PUT_ATT_TEXT(ncid_out,IVAR_ID,'units',LEN_TRIM(XUNIT(IVAR)),TRIM(XUNIT(IVAR)))
  CALL HANDLE_ERR(IERR)
 END DO  ! ivar

 ! define ordinates of probability distributions
 ! IERR = NF_DEF_VAR(ncid_out,'ordinates',NF_REAL,1,NORD_DIM,IORD_ID); CALL HANDLE_ERR(IERR)
 ! IERR = NF_PUT_ATT_TEXT(ncid_out,IORD_ID,'long_name',37,'ordinates of probability distribution')
 ! CALL HANDLE_ERR(IERR)

 ! IERR = NF_PUT_ATT_TEXT(ncid_out,IORD_ID,'units',1,'-'); CALL HANDLE_ERR(IERR)

 ! define probability distributions
 ! PVAR = (/NPAR_DIM,NORD_DIM/)   ! dimensions for probability distributions
 ! IERR = NF_DEF_VAR(ncid_out,'probability',NF_REAL,2,PVAR,IVAR_ID); CALL HANDLE_ERR(IERR)
 ! IERR = NF_PUT_ATT_TEXT(ncid_out,IVAR_ID,'long_name',44,'cumulative probability of number of substeps'); CALL HANDLE_ERR(IERR)
 ! IERR = NF_PUT_ATT_TEXT(ncid_out,IVAR_ID,'units',1,'-'); CALL HANDLE_ERR(IERR)

! end definitions and close file
IERR = NF_ENDDEF(ncid_out)
! write the ordinates of the probability distribution
!rORD = real(ORD_NSUBS,kind(MSP))
!IERR = NF_PUT_VAR_REAL(ncid_out,IORD_ID,rORD); CALL HANDLE_ERR(IERR)     ! write data
IERR = NF_CLOSE(ncid_out)

! ---------------------------------------------------------------------------------------
END SUBROUTINE DEF_SSTATS
