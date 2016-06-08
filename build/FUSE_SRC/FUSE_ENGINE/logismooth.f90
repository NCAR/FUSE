PURE FUNCTION LOGISMOOTH(STATE,STATE_MAX,PSMOOTH)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Uses a logistic function to smooth the threshold at the top of a bucket
! ---------------------------------------------------------------------------------------
USE nrtype
IMPLICIT NONE
REAL(SP), INTENT(IN)                   :: STATE       ! model state
REAL(SP), INTENT(IN)                   :: STATE_MAX   ! maximum model state
REAL(SP), INTENT(IN)                   :: PSMOOTH     ! smoothing parameter (fraction of state)
REAL(SP)                               :: ASMOOTH     ! actual smoothing
REAL(SP)                               :: LOGISMOOTH  ! FUNCTION name
! ---------------------------------------------------------------------------------------
ASMOOTH = PSMOOTH*STATE_MAX                           ! actual smoothing
LOGISMOOTH = 1._SP / ( 1._SP + EXP(-(STATE - (STATE_MAX - ASMOOTH*5._SP) ) / ASMOOTH) )
! ---------------------------------------------------------------------------------------
END FUNCTION LOGISMOOTH
