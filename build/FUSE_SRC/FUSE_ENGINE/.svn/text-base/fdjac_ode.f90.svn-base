MODULE fdjac_ode_module
IMPLICIT NONE
CONTAINS
SUBROUTINE fdjac_ode(x,dsdt,df,simeth)
USE nrtype; USE nrutil, ONLY : assert_eq
USE model_numerix, ONLY : num_jacobian
USE fuse_deriv_module
! Used to compute Jacobian of the ODE, based on the NR routine fdjac
IMPLICIT NONE
! input/output
REAL(SP), DIMENSION(:), INTENT(IN)    :: dsdt       ! state derivative
REAL(SP), DIMENSION(:), INTENT(INOUT) :: x          ! trial state vector
REAL(SP), DIMENSION(:,:), INTENT(OUT) :: df         ! Jacobian
LOGICAL(LGT), INTENT(IN), OPTIONAL    :: simeth     ! flag for semi-implicit Euler method
! internal
LOGICAL(LGT)                          :: fdflux     ! flag to compute flux derivatives
REAL(SP), PARAMETER                   :: EPS=-1.0e-4_sp  ! relative state change, NOTE force h to be negative
INTEGER(I4B)                          :: j,n        ! loop through statesm number of states
REAL(SP)                              :: dx         ! relative change in state
REAL(SP), DIMENSION(size(x))          :: xsav,xph,h ! perturbed states and change in states
! check size of input argumets
n=assert_eq(size(x),size(dsdt),size(df,1),size(df,2),'fdjac')
! if semi-implicit Euler method, then compute flux derivatives
fdflux=.false.; if (present(simeth)) fdflux=.true.
! save input x value
xsav=x
! compute step size
dx = EPS                                     ! relative state change
!DK: dx-determination can be improved using the characteristic scale of state variables.
!    current approach ok for the moment.
h  = dx*abs(xsav)                            ! state change
where (h == 0.0) h=dx                        ! force state change to be non-zero
xph= xsav+h                                  ! perturbed state
h  = xph-xsav                                ! size of perturbation (trick to avoid rounding errors)
! compute Jacobian (and, if desired, compute the derivatives of the fluxes)
do j=1,n
 x(j)=xph(j)                                 ! perturb state
 !print *, 'computing jacobian, j, x = '; write(*,'(i3,1x,10(e20.8,1x))') j, x
 df(:,j)=(fuse_deriv(x)-dsdt(:))/h(j)        ! compute row of the Jacobian
 !print *, 'jac result '; write(*,'(10(e20.8,1x))') df(:,j)
 if (fdflux) call flux_deriv(j,h(j))         ! compute flux derivatives for state j
 x(j)=xsav(j)                                ! set state back to original value
end do
! keep track of the number of times computing the Jacobian
num_jacobian = num_jacobian + 1              ! num_jacobian shared in module model_numerix
END SUBROUTINE fdjac_ode
END MODULE fdjac_ode_module
