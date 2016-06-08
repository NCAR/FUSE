MODULE fminln
  USE nrtype; USE nrutil, ONLY : nrerror
  REAL(SP), POINTER               :: fmin_dtp    ! time step
  REAL(SP), POINTER               :: fmin_dt2p   ! half time step
  REAL(SP), DIMENSION(:), POINTER :: fmin_x0p    ! initial state
  REAL(SP), DIMENSION(:), POINTER :: fmin_dseep  ! change in state by explicit euler
  REAL(SP), DIMENSION(:), POINTER :: fmin_dsdtp  ! state derivatives
  REAL(SP), DIMENSION(:), POINTER :: fmin_fvecp  ! residuals of the discrete system
CONTAINS
!BL
  FUNCTION fmin(x)
  USE model_numerix      ! provide access to the model numerix decisions
  USE fuse_deriv_module  ! provide access to the function to compute model derivatives
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: x
  REAL(SP) :: fmin
  if (.not.associated(fmin_x0p) .or. .not.associated(fmin_dtp) .or. &
      .not.associated(fmin_dt2p) .or. .not.associated(fmin_dseep) .or. &
      .not.associated(fmin_dsdtp) .or. .not.associated(fmin_fvecp) ) &
    call nrerror('fmin: problem with pointer for returned values')
  fmin_dsdtp = fuse_deriv(x)      ! calculate derivatives
  SELECT CASE(SOLUTION_METHOD)
   CASE(IMPLICIT_EULER); fmin_fvecp = x - (fmin_x0p + fmin_dtp*fmin_dsdtp)
   !print *, 'in fmin, x = ', x
   !print *, 'in fmin, x0 + dt * dsdt = ', fmin_x0p + fmin_dtp*fmin_dsdtp
   !print *, 'in fmin, fvec = ', fmin_fvecp
   CASE(IMPLICIT_HEUN);  fmin_fvecp = x - (fmin_x0p + fmin_dt2p*fmin_dsdtp + fmin_dseep)
   CASE DEFAULT; call nrerror('fmin: solution method must be either implicit euler or implicit heun')
  END SELECT
  fmin=0.5_sp*dot_product(fmin_fvecp,fmin_fvecp)
  END FUNCTION fmin 
END MODULE fminln
