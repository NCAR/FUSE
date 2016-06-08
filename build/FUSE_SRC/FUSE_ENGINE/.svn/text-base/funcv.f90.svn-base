module funcv_mod
implicit none
contains
FUNCTION FUNCV(X_TRY)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Returns a vector of errors from the state equations, evaluated at X_TRY
!
! That is,
!  X_NEW(1) = X_OLD(1) + DYDX( X_TRY(1) ) * delT
!  X_NEW(2) = X_OLD(2) + DYDX( X_TRY(2) ) * delT
!  ...
!  X_NEW(N) = X_OLD(N) + DYDX( X_TRY(N) ) * delT
!
! So...
!  FUNCV(1) = X_NEW(1) - X_TRY(1)
!  FUNCV(2) = X_NEW(2) - X_TRY(2)
!
!  FUNCV(N) = X_NEW(N) - X_TRY(N)
!
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn, ONLY:CSTATE,NSTATE                    ! model definition structures
USE model_defnames
USE multistate, ONLY:TSTATE,MSTATE,DY_DT,HSTATE       ! model states
USE xtry_2_str_module                                 ! puts state vector into structure in multistate
IMPLICIT NONE
! input/output
REAL(SP), DIMENSION(:), INTENT(IN)     :: X_TRY       ! vector of model states
REAL(SP), DIMENSION(SIZE(X_TRY))       :: FUNCV       ! function evaluations
! internal
INTEGER(I4B)                           :: ISTT        ! loop through model states
! ---------------------------------------------------------------------------------------
! (1) COMPUTE MODEL DERIVATIVES
! ---------------------------------------------------------------------------------------
CALL XTRY_2_STR(X_TRY,TSTATE)  ! populate state structure TSTATE with values of X
CALL MOD_DERIVS()              ! evaluate dxdt for state vector X_TRY
! ---------------------------------------------------------------------------------------
! (2) COMPUTE FUNCTION VALUES
! ---------------------------------------------------------------------------------------
DO ISTT=1,NSTATE
 SELECT CASE(CSTATE(ISTT)%iSNAME)
  CASE (iopt_TENS1A); FUNCV(ISTT) = MSTATE%TENS_1A + DY_DT%TENS_1A*HSTATE%STEP - X_TRY(ISTT)
  CASE (iopt_TENS1B); FUNCV(ISTT) = MSTATE%TENS_1B + DY_DT%TENS_1B*HSTATE%STEP - X_TRY(ISTT)
  CASE (iopt_TENS_1); FUNCV(ISTT) = MSTATE%TENS_1  + DY_DT%TENS_1 *HSTATE%STEP - X_TRY(ISTT)
  CASE (iopt_FREE_1); FUNCV(ISTT) = MSTATE%FREE_1  + DY_DT%FREE_1 *HSTATE%STEP - X_TRY(ISTT)
  CASE (iopt_WATR_1); FUNCV(ISTT) = MSTATE%WATR_1  + DY_DT%WATR_1 *HSTATE%STEP - X_TRY(ISTT)
  CASE (iopt_TENS_2); FUNCV(ISTT) = MSTATE%TENS_2  + DY_DT%TENS_2 *HSTATE%STEP - X_TRY(ISTT)
  CASE (iopt_FREE2A); FUNCV(ISTT) = MSTATE%FREE_2A + DY_DT%FREE_2A*HSTATE%STEP - X_TRY(ISTT)
  CASE (iopt_FREE2B); FUNCV(ISTT) = MSTATE%FREE_2B + DY_DT%FREE_2B*HSTATE%STEP - X_TRY(ISTT)
  CASE (iopt_WATR_2); FUNCV(ISTT) = MSTATE%WATR_2  + DY_DT%WATR_2 *HSTATE%STEP - X_TRY(ISTT)
  CASE DEFAULT; STOP 'fatal error: cannot identify the state variable'
 END SELECT
 print *, desc_int2str(CSTATE(ISTT)%iSNAME), FUNCV(ISTT), HSTATE%STEP
END DO
! ---------------------------------------------------------------------------------------
END FUNCTION FUNCV
endmodule funcv_mod
