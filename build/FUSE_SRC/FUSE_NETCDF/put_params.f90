SUBROUTINE PUT_PARAMS(IPAR)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! Modified by Nans Addor to include snow module
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! write NetCDF output files  -- model parameters
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition structures (includes filename)
USE model_defnames                                    ! define variable names
USE metaparams                                        ! metadata for model parameters
USE multistats, ONLY:MSTATS                           ! provide access to error message
USE parextract_module                                 ! extract parameters
IMPLICIT NONE
! input
INTEGER(I4B), INTENT(IN)               :: IPAR        ! parameter set index
! internal
INTEGER(I4B)                           :: IERR,NCID   ! error code; NetCDF ID
INTEGER(I4B), DIMENSION(1)             :: INDX        ! indices for parameter write
INTEGER(I4B)                           :: IVAR        ! loop through parameters
REAL(SP)                               :: XPAR        ! desired parameter
REAL(MSP)                              :: APAR        ! convert to SP (need for SP write)
INTEGER(I4B)                           :: IVAR_ID     ! variable ID
INTEGER(I4B), PARAMETER                :: NDESC=9     ! number of model descriptors - TODO: THIS SHOULDN'T BE HARD-CODED
INTEGER(I4B), PARAMETER                :: NCHAR=10    ! length of model descriptors - TODO: THIS SHOULDN'T BE HARD-CODED
INTEGER(I4B), DIMENSION(3)             :: ISTART      ! starting position for array write
INTEGER(I4B), DIMENSION(3)             :: ICOUNT      ! count for array write
CHARACTER(LEN=10)                      :: TXTVEC      ! single model descriptor
include 'netcdf.inc'                                  ! use netCDF libraries
! ---------------------------------------------------------------------------------------

! open file
IERR = NF_OPEN(TRIM(FNAME_NETCDF_PARA),NF_WRITE,NCID); CALL HANDLE_ERR(IERR)

 ! define indices for model output
 INDX = (/IPAR/)

 ! loop through model parameters
 DO IVAR=1,NOUTPAR  ! NOUTPAR is stored in module metaparams

  XPAR = PAREXTRACT(PNAME(IVAR)); APAR=XPAR                                  ! get parameter PNAME(IVAR)
  IERR = NF_INQ_VARID(NCID,TRIM(PNAME(IVAR)),IVAR_ID); CALL HANDLE_ERR(IERR) ! get variable ID
  IERR = NF_PUT_VAR1_REAL(NCID,IVAR_ID,INDX,APAR); CALL HANDLE_ERR(IERR)     ! write data

 END DO  ! (ivar)

 ! put model description
 !IERR = NF_INQ_VARID(NCID,'model_description',IVAR_ID); CALL HANDLE_ERR(IERR)

 ! print *, 'Writing model decisions to this NetCDF file:', TRIM(FNAME_NETCDF)
 !
 ! DO IVAR=1,NDESC
 !  ! extract text string
 !  IF (IVAR.EQ.1) TXTVEC = desc_int2str(SMODL%iRFERR)
 !  IF (IVAR.EQ.2) TXTVEC = desc_int2str(SMODL%iARCH1)
 !  IF (IVAR.EQ.3) TXTVEC = desc_int2str(SMODL%iARCH2)
 !  IF (IVAR.EQ.4) TXTVEC = desc_int2str(SMODL%iQSURF)
 !  IF (IVAR.EQ.5) TXTVEC = desc_int2str(SMODL%iQPERC)
 !  IF (IVAR.EQ.6) TXTVEC = desc_int2str(SMODL%iESOIL)
 !  IF (IVAR.EQ.7) TXTVEC = desc_int2str(SMODL%iQINTF)
 !  IF (IVAR.EQ.8) TXTVEC = desc_int2str(SMODL%iQ_TDH)
 !  IF (IVAR.EQ.9) TXTVEC = desc_int2str(SMODL%iSNOWM)
 !
 !  ISTART = (/    1,IVAR,IMOD/)   ! starting position of array
 !  ICOUNT = (/NCHAR,   1,   1/)   ! number of array elements (one descriptor, one model)
 !  IERR = NF_PUT_VARA_TEXT(NCID,IVAR_ID,ISTART,ICOUNT,TXTVEC); CALL HANDLE_ERR(IERR)
 ! END DO
 ! put error message
 !ISTART = (/                      1,IMOD,IPAR/)   ! starting position of array
 !ICOUNT = (/LEN(MSTATS%ERR_MESSAGE),   1,   1/)   ! number of array elements (one descriptor, one model)
 !IERR = NF_INQ_VARID(NCID,'error_message',IVAR_ID); CALL HANDLE_ERR(IERR)
 !IERR = NF_PUT_VARA_TEXT(NCID,IVAR_ID,ISTART,ICOUNT,MSTATS%ERR_MESSAGE); CALL HANDLE_ERR(IERR)
! close NetCDF file
IERR = NF_CLOSE(NCID)
! ---------------------------------------------------------------------------------------
END SUBROUTINE PUT_PARAMS
