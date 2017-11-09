SUBROUTINE PUT_SSTATS(IPAR)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! write NetCDF output files -- summary statistics
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition structures (includes filename)
USE meta_stats                                        ! metadata for summary statistics
USE multistats                                        ! model summary statistics
USE model_numerix                                     ! model numerix parameters and arrays
USE sumextract_module                                 ! module to extract summary statistics
IMPLICIT NONE
! input
INTEGER(I4B), INTENT(IN)               :: IPAR        ! parameter set index
! internal
INTEGER(I4B)                           :: IERR,NCID   ! error code; NetCDF ID
INTEGER(I4B), DIMENSION(1)             :: INDX        ! indices for parameter write
INTEGER(I4B)                           :: IVAR        ! loop through parameters
REAL(SP)                               :: XPAR        ! desired parameter (SP may not be SP)
REAL(MSP)                              :: APAR        ! desired parameter (...but MSP is SP)
INTEGER(I4B)                           :: IVAR_ID     ! variable ID
!INTEGER(I4B), DIMENSION(2)             :: iBeg        ! start index
!INTEGER(I4B), DIMENSION(2)             :: iCnt        ! count
!REAL(MSP),dimension(size(MSTATS%NUMSUB_PROB))  :: TVEC      ! temporary vector
include 'netcdf.inc'                                  ! use netCDF libraries
! ---------------------------------------------------------------------------------------
! open NetCDF parameter file
IERR = NF_OPEN(TRIM(FNAME_NETCDF_PARA),NF_WRITE,NCID); CALL HANDLE_ERR(IERR)

 ! define indices for model output
 INDX = (/IPAR/)

 ! loop through summary statistics
 DO IVAR=1,NSUMVAR
  XPAR = SUMEXTRACT(XNAME(IVAR)); APAR=XPAR                                  ! get parameter ivar
  IERR = NF_INQ_VARID(NCID,TRIM(XNAME(IVAR)),IVAR_ID); CALL HANDLE_ERR(IERR) ! get variable ID
  IERR = NF_PUT_VAR1_REAL(NCID,IVAR_ID,INDX,APAR); CALL HANDLE_ERR(IERR)     ! write data
 END DO  ! (ivar)

 ! write probability distribution of number of substeps
 ! iBeg = (/INDX(1),1/)                          ! start index
 ! iCnt = (/1,SIZE(PRB_NSUBS)/)                  ! count
 ! TVEC = REAL(MSTATS%NUMSUB_PROB, KIND(MSP))    ! data

 ! IERR = NF_INQ_VARID(NCID,'probability',IVAR_ID); CALL HANDLE_ERR(IERR)  ! get variable ID
 ! IERR = NF_PUT_VARA_REAL(NCID,IVAR_ID,iBeg,iCnt,TVEC);  CALL HANDLE_ERR(IERR)  ! write data

IERR = NF_CLOSE(NCID)
! ---------------------------------------------------------------------------------------
END SUBROUTINE PUT_SSTATS
