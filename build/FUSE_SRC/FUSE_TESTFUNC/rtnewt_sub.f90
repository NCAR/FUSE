SUBROUTINE rtnewt_sub(funcd,xold,x1,x2,xacc,xnew,niter)
! From Numerical Recipes, but converted from a function to a subroutine
USE nrtype; USE nrutil, ONLY : nrerror
IMPLICIT NONE
REAL(SP), INTENT(IN)      :: xold,x1,x2,xacc
REAL(SP), INTENT(OUT)     :: xnew
INTEGER(I4B), INTENT(OUT) :: niter
INTERFACE
 SUBROUTINE funcd(x,fval,fderiv)
 USE nrtype
 IMPLICIT NONE
 REAL(SP), INTENT(IN)  :: x
 REAL(SP), INTENT(OUT) :: fval,fderiv
 END SUBROUTINE funcd
END INTERFACE
INTEGER(I4B), PARAMETER :: MAXIT=20
INTEGER(I4B) :: j
REAL(SP) :: df,dx,f,xsave
xnew = xold
if (xnew < x1) xnew=x1
if (xnew > x2) xnew=x2
do j=1,MAXIT
 call funcd(xnew,f,df)    ! calculate function and derivative
 dx   =f/df               ! calculate dx
 xsave=xnew               ! save last trial value
 xnew =xsave-dx           ! calculate next trial value
 if (xnew < x1) xnew=x1   ! check > minimum
 if (xnew > x2) xnew=x2   ! check < maximum
 if (abs(xnew-xsave) < xacc .or. abs(dx) < xacc) then ! check for convergence
  niter=j                 ! save number of iterations
  RETURN
 endif
end do
call nrerror('rtnewt exceeded maximum iterations')
END SUBROUTINE rtnewt_sub
