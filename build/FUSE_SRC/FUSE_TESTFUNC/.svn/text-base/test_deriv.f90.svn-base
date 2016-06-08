MODULE TEST_DERIV__MODULE
IMPLICIT NONE
CONTAINS
! ---------------------------------------------------------------------------------------
FUNCTION TEST_DERIV(S)
! Used to calculate derivatives using the simple test function
!  dS/dt = -sqrt(S)
! For generality, includes
!  (1) Put state vector in model data structures
!  (2) Compute fluxes
!  (3) Compute derivatives
!  (4) Extract derivatives from model structure
USE nrtype                                        ! numerical recipes data types 
USE test_modvar, ONLY: TSTATE,M_FLUX,MDS_DT       ! model data structures
USE model_numerix, ONLY: NUM_FUNCS                ! number of function calls
IMPLICIT NONE
REAL(SP), DIMENSION(:), INTENT(IN)  :: S          ! storage
REAL(SP), DIMENSION(SIZE(S))        :: TEST_DERIV ! FUNCTION name
NUM_FUNCS       = NUM_FUNCS + 1        ! (0) Keep track of the number of function calls
TSTATE%WATR_1   = S(1)                 ! (1) Put state vector in model data structures
M_FLUX%DRAINAGE = SQRT(TSTATE%WATR_1)  ! (2) Compute fluxes
MDS_DT%WATR_1   = -M_FLUX%DRAINAGE     ! (3) Compute derivatives
TEST_DERIV(1)   = MDS_DT%WATR_1        ! (4) Extract derivatives from model structure
END FUNCTION TEST_DERIV
! ---------------------------------------------------------------------------------------
END MODULE TEST_DERIV__MODULE
