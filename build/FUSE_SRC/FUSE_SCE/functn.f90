FUNCTION FUNCTN(NOPT,A)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Wrapper for SCE (used to compute the objective function)
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE FUSE_RMSE_MODULE                           		  ! run model and compute the root mean squared error
IMPLICIT NONE
! input
INTEGER(I4B)                           :: NOPT        ! number of parameters
REAL(SP), DIMENSION(100), INTENT(IN)  :: A            ! model parameter set - can be bumped up to 100 elements

! internal
REAL(SP), DIMENSION(:), ALLOCATABLE    :: SCE_PAR     ! sce parameter set
INTEGER(I4B)                           :: IERR        ! error code for allocate/deallocate
INTEGER(I4B)                           :: ERR         ! error code for fuse_rmse
CHARACTER(LEN=256)                     :: MESSAGE     ! error message for fuse_rmse
LOGICAL(LGT)                           :: OUTPUT_FLAG ! .TRUE. = write model time series
REAL(SP)                               :: RMSE        ! root mean squared error
! output
REAL(SP)                              :: FUNCTN      ! objective function value

! ---------------------------------------------------------------------------------------
! get SCE parameter set
ALLOCATE(SCE_PAR(NOPT), STAT=IERR); IF (IERR.NE.0) STOP ' problem allocating space '
SCE_PAR(1:NOPT) = A(1:NOPT)  ! convert from MSP used in SCE to SP used in FUSE

!print *, 'A_MSP', A(1:3)
!print *, 'SCE_PAR', SCE_PAR(1:3)

! compute RMSE
OUTPUT_FLAG=.true.          ! .TRUE. = write model time series
CALL FUSE_RMSE(SCE_PAR,RMSE,OUTPUT_FLAG,ERR=ERR,MESSAGE=MESSAGE)

IF (ERR.NE.0) THEN
  print *, MESSAGE
  return
ENDIF

! deallocate parameter set
DEALLOCATE(SCE_PAR, STAT=IERR); IF (IERR.NE.0) STOP ' problem deallocating space '
print *, 'RMSE =', RMSE
! save objective function value
FUNCTN = RMSE 
! ---------------------------------------------------------------------------------------
END FUNCTION FUNCTN
