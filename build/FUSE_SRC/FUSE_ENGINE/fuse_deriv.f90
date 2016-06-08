MODULE FUSE_DERIV_MODULE
IMPLICIT NONE
CONTAINS
! ---------------------------------------------------------------------------------------
FUNCTION FUSE_DERIV(S)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Used to calculate derivatives from a specified FUSE model, includes
!  (1) Put state vector in model data structures
!  (2) Compute fluxes and derivatives
!  (3) Extract derivatives from model structure
! ---------------------------------------------------------------------------------------
USE nrtype                                        ! numerical recipes data types 
USE multistate, ONLY:TSTATE,DY_DT                 ! model data structures
USE str_2_xtry_module                             ! provide access to str_2_xtry
USE xtry_2_str_module                             ! provide access to xtry_2_str
IMPLICIT NONE
REAL(SP), DIMENSION(:), INTENT(IN)  :: S          ! storage
REAL(SP), DIMENSION(SIZE(S))        :: FUSE_DERIV ! FUNCTION name
CALL XTRY_2_STR(S,TSTATE)              ! (1) Put state vector in model data structures
CALL MOD_DERIVS()                      ! (2) Compute fluxes and derivatives
CALL STR_2_XTRY(DY_DT,FUSE_DERIV)      ! (3) Extract derivatives from model structure
END FUNCTION FUSE_DERIV
! ---------------------------------------------------------------------------------------
END MODULE FUSE_DERIV_MODULE
