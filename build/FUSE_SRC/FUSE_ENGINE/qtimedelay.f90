SUBROUTINE QTIMEDELAY(err,message)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Computes the fraction of runoff in future time steps
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multiparam -- runoff fractions stored in DPARAM%FRAC_FUTURE(:)
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE nr, ONLY : gammp                                  ! interface for the incomplete gamma function
USE model_defn                                        ! model definition structure
USE model_defnames
USE multiforce                                        ! model forcing (need DELTIM)
USE multiparam                                        ! model parameters
IMPLICIT NONE
! dummies
integer(i4b),intent(out)::err
character(*),intent(out)::message
! locals
INTEGER(I4B)                           :: NTDH        ! maximum number of future time steps
REAL(SP)                               :: ALPHA       ! shape parameter
REAL(SP)                               :: ALAMB       ! scale parameter
INTEGER(I4B)                           :: JTIM        ! (loop through future time steps)
REAL(SP)                               :: TFUTURE     ! future time (units of days)
REAL(SP)                               :: CUMPROB     ! cumulative probability at JTIM
REAL(SP)                               :: PSAVE       ! cumulative probability at JTIM-1
! ---------------------------------------------------------------------------------------
err=0
SELECT CASE(SMODL%iQ_TDH)
 CASE(iopt_rout_gamma) ! use a Gamma distribution with shape parameter = 2.5
  ALPHA = 2.5_SP                                             ! shape parameter

  !PRINT *, 'MPARAM= ', MPARAM

  ALAMB = ALPHA/MPARAM%TIMEDELAY                             ! scale parameter
  PSAVE = 0._SP                                              ! cumulative probability at JTIM-1
  NTDH = SIZE(DPARAM%FRAC_FUTURE)                            ! maximum number of future time steps
  ! loop through time steps and compute the fraction of runoff in future time steps
  DO JTIM=1,NTDH
   TFUTURE                   = REAL(JTIM,SP)*DELTIM          ! future time (units of days)
   CUMPROB                   = GAMMP(ALPHA,ALAMB*TFUTURE)    ! cumulative probability at JTIM
   DPARAM%FRAC_FUTURE(JTIM)  = MAX(0._SP, CUMPROB-PSAVE)     ! probability between JTIM-1 and JTIM
   PSAVE                     = CUMPROB                       ! cumulative probability at JTIM-1
   !WRITE(*,'(3(F11.5))') TFUTURE, DPARAM%FRAC_FUTURE(JTIM), CUMPROB
   IF(DPARAM%FRAC_FUTURE(JTIM)<EPSILON(1._SP))EXIT
  END DO
  DPARAM%NTDH_NEED = MIN(JTIM,NTDH)
  DPARAM%FRAC_FUTURE(DPARAM%NTDH_NEED+1:)=0._SP
  ! check there are enough bins
  IF (CUMPROB.LT.0.99_SP) THEN
   err=100; message='f-QTIMEDELAY/not enough bins in dparam%frac_future'
   return
  ENDIF
  ! ensure that the fractions sum to 1.0 (account for rounding errors, and not enough bins)
  DPARAM%FRAC_FUTURE(:) = DPARAM%FRAC_FUTURE(:) / SUM(DPARAM%FRAC_FUTURE(:))
 CASE(iopt_no_routing) ! no routing
  NTDH                       = SIZE(DPARAM%FRAC_FUTURE)
  DPARAM%NTDH_NEED           = 2
  DPARAM%FRAC_FUTURE(1)      = 1._SP
  DPARAM%FRAC_FUTURE(2:NTDH) = 0._SP
 CASE DEFAULT       ! check for errors
  err=100; message="f-QTIMEDELAY/SMODL%iQ_TDH must be either iopt_rout_gamma or iopt_no_routing"
  return
END SELECT
! ---------------------------------------------------------------------------------------
END SUBROUTINE QTIMEDELAY
