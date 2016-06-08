SUBROUTINE IMPL_ERROR(S,F,DF)
! Used to calculate the error for the implicit scheme
!  S(n+1) = S(n) + dS(n+1)/dt * delT
!  F = S(try) - (S(n) + dS(try)/dt * delT)
USE nrtype                                          ! numerical recipes data types
USE model_numerix, ONLY: NUM_JACOBIAN               ! number of times calculate the derivative
USE test_modvar, ONLY: MSTATE,MDS_DT,DT_SUB         ! model variables
USE test_deriv__module                              ! provide access to model derivatives function
IMPLICIT NONE
! input/output
REAL(SP), INTENT(IN)                :: S            ! storage
REAL(SP), INTENT(OUT)               :: F            ! function value
REAL(SP), INTENT(OUT)               :: DF           ! function derivative
! internal
REAL(SP)                            :: S0           ! state at the start of the sub-step
REAL(SP), PARAMETER                 :: RH=1.e-4_sp  ! relative step size for finite difference
REAL(SP)                            :: H            ! step size for finite difference
REAL(SP)                            :: SPH          ! perturbed state
REAL(SP), DIMENSION(1)              :: DSDT         ! state derivative (NOTE, pass as vector)
REAL(SP)                            :: FTRY         ! perturbed function value
! keep track of the number of times calculate the derivative
NUM_JACOBIAN = NUM_JACOBIAN + 1
! extract state at the start of the time step
S0   = MSTATE%WATR_1
! calculate perturbed function value
H    = RH*S                                         ! step size
SPH  = S+H                                          ! perturbed state
H    = SPH-S                                        ! actual step size (trick to account for roundoff errors)
DSDT = TEST_DERIV((/SPH/))                          ! calculate state derivative (NOTE, pass as vector)
FTRY = SPH - (S0 + DSDT(1)*DT_SUB)                  ! perturbed function value
! calculate function value
DSDT = TEST_DERIV((/S/))                            ! calculate state derivative (NOTE, pass as vector)
F    = S - (S0 + DSDT(1)*DT_SUB)                    ! calculate function value
MDS_DT%WATR_1 = DSDT(1)                             ! save state derivative
! calculate function derivative
DF   = (FTRY-F)/H                                   ! function derivative
END SUBROUTINE IMPL_ERROR
