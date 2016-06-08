SUBROUTINE fdjac(x,fvec,df)
USE nrtype; USE nrutil, ONLY : assert_eq
use funcv_mod
use model_numerix, ONLY : num_jacobian
IMPLICIT NONE
REAL(SP), DIMENSION(:), INTENT(IN) :: fvec
REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
REAL(SP), DIMENSION(:,:), INTENT(OUT) :: df
!INTERFACE
! FUNCTION funcv(xtry)
! USE nrtype
! IMPLICIT NONE
! REAL(SP), DIMENSION(:), INTENT(IN) :: xtry
! REAL(SP), DIMENSION(size(xtry)) :: funcv
! END FUNCTION funcv
!END INTERFACE
REAL(SP), PARAMETER :: EPS=-1.0e-4_sp  ! NOTE force h to be negative
INTEGER(I4B) :: j,n
REAL(SP), DIMENSION(size(x)) :: xsav,xph,h
REAL(SP), DIMENSION(size(df,1)) :: vv
n=assert_eq(size(x),size(fvec),size(df,1),size(df,2),'fdjac')
xsav=x
h=EPS*abs(xsav)
where (h == 0.0) h=EPS
xph=xsav+h
h=xph-xsav
do j=1,n
 x(j)=xph(j)
 df(:,j)=(funcv(x)-fvec(:))/h(j)
 x(j)=xsav(j)
end do
! MPC check for zero derivative
vv=maxval(abs(df),dim=2)
if (any(vv == 0.0)) then
 do j=1,n; write(*,'(10(e12.5,1x))') df(:,j); end do
 stop ' fatal error: zero derivative in Jacobian '
endif
! keep track of the number of times computing the Jacobian
num_jacobian = num_jacobian + 1  ! num_jacobian shared in module model_numerix
END SUBROUTINE fdjac
