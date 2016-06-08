MODULE newtoniter_mod
IMPLICIT NONE
CONTAINS
  ! ---------------------------------------------------------------------------------------
  SUBROUTINE newtoniter(X,newJacIn,check,niter)
  USE nrtype; USE nrutil, ONLY : nrerror,diagadd,vabs
  USE nr, ONLY : fdjac,lnsrch,lubksb,ludcmp
  USE fminln, ONLY : fmin,fmin_dsdtp,fmin_fvecp,fmin_dtp,fmin_dt2p
  USE limit_xtry_module
  USE fdjac_ode_module
  USE model_numerix
  ! Purpose: finds the state vector "X_NEW", so that
  !  X_NEW(:) = X_OLD(:) + DYDX(:) * HSTATE%STEP, with DYDX(:) evaluated at X_NEW(:)
  ! (based loosely on the Numerical Recipes routine newt.f90)
  ! Programmers: Dmitri Kavetski and Martyn Clark
  IMPLICIT NONE
  ! dummies
  REAL(SP), DIMENSION(:), INTENT(INOUT)    :: x         ! state vector
  LOGICAL(LGT), INTENT(IN)                 :: newJacIn  ! .TRUE. if new Jacobian required
  LOGICAL(LGT), INTENT(OUT)                :: check     ! .TRUE. if spurious minimum
  INTEGER(I4B), INTENT(OUT)                :: niter     ! number of iterations
  ! algorithmic control parameters (most passed through MODULE model_numerix)
  REAL(SP), PARAMETER                      :: TOLMIN=1.0e-10_sp ! check for spurious minima
  REAL(SP), PARAMETER                      :: STPMX=100.0_sp    ! maximum step in lnsrch
  ! locals
  INTEGER(I4B)                             :: i,j,k     ! looping (test)
  INTEGER(I4B)                             :: its       ! iteration counter
  INTEGER(I4B), DIMENSION(size(x))         :: indx      ! used in ludcmp
  REAL(SP)                                 :: d         ! used in ludcmp
  REAL(SP)                                 :: f,fold    ! function values
  REAL(SP)                                 :: absf_old  ! absolute value of the residual vector (last iter)
  REAL(SP)                                 :: absf_new  ! absolute value of the residual vector (current iter)
  REAL(SP)                                 :: stpmax    ! step size for lnsrch
  REAL(SP), DIMENSION(size(x))             :: g         ! gradient used in lnsrch
  REAL(SP), DIMENSION(size(x))             :: p,dx      ! p = newton step, dx = actual step
  REAL(SP), DIMENSION(size(x))             :: xold      ! old state vector
  REAL(SP), DIMENSION(size(x)), TARGET     :: dsdt      ! model derivatives 
  REAL(SP), DIMENSION(size(x)), TARGET     :: fvec      ! model residuals
  REAL(SP), DIMENSION(size(x),size(x))     :: jac_ode,fjac,fjacSave ! Jacobian matrices
  LOGICAL(LGT)                             :: newjac    ! .TRUE. if calculate a new Jacobian matrix

  ! ---------------------------------------------------------------------------------------
  ! (0) INITIALIZATION
  ! ---------------------------------------------------------------------------------------
  NITER=0          ! initialize number of iterations (intent=out)
  CHECK=.FALSE.    ! initialize check on convergence (intent=out)
  fmin_dsdtp=>dsdt ! provide access to the vector of derivatives used in fmin
  fmin_fvecp=>fvec ! provide access to the vector of residuals used in fmin

  ! ---------------------------------------------------------------------------------------
  ! (1) TEST FOR THE INITIAL GUESS BEING A ROOT (MORE STRINGENT TEST THAN SIMPLY ERR_ITER_DX)
  ! ---------------------------------------------------------------------------------------
  CALL LIMIT_XTRY(X)         ! ensure that the values of X are physically reasonable
  F=FMIN(X)                  ! compute function evaluation (populates vectors DSDT and FVEC)
  !write(*,'(10(f20.10,1x))') x
  ABSF_OLD=MAXVAL(ABS(FVEC)) ! initial norm of the residual vector
  IF (ABSF_OLD < 0.01_SP*ERR_ITER_DX) THEN
   CHECK=.FALSE.
   RETURN
  ENDIF

  ! ---------------------------------------------------------------------------------------
  ! (2) ITERATE TO EITHER NITER_TOTAL OR CONVERGENCE
  ! ---------------------------------------------------------------------------------------
  ! compute maximum step size used in line searches
  IF (CHECK_OVERSHOOT.EQ.LINE_SEARCH) STPMAX = STPMX*MAX(VABS(X),REAL(SIZE(X),SP))
  DO ITS=1,NITER_TOTAL
   NITER = ITS

   !print *, '***** new iteration *****', its, check, ABSF_OLD
   ! ---------------------------------------------------------------------------------------
   ! (2A) CHECK IF WE NEED A NEW JACOBIAN, AND, IF SO, COMPUTE IT
   ! ---------------------------------------------------------------------------------------
   SELECT CASE(JAC_RECOMPUTE)
    CASE(FULLYVARIABLE)
     NEWJAC=.TRUE.                ! always re-compute Jacobian
    CASE(CONST_SUBSTEP)
     NEWJAC=(ITS==1)              ! only recompute Jacobian on the first iteration
    CASE(CONSTFULLSTEP,PERIOD_FREEZE,SMALL_F_RATIO)
     IF (JAC_RECOMPUTE==CONSTFULLSTEP) THEN
      NEWJAC=newJacIn              ! only recompute Jacobian at start of full step (defined by input flag)
      IF (ITS>1) NEWJAC=.FALSE.
     ENDIF
     IF (JAC_RECOMPUTE==PERIOD_FREEZE) THEN
      NEWJAC=(MAXVAL(ABS(FVEC)) > THRESH_FRZE)
      IF (ITS==1) NEWJAC=newJacIn  ! always recompute Jacobian at start of full step (defined by input flag)
     ENDIF
     IF (JAC_RECOMPUTE==SMALL_F_RATIO) THEN
      IF (ITS==1) THEN
       NEWJAC=.TRUE.
      ELSE
       NEWJAC=(ABSF_NEW/ABSF_OLD > THRESH_FRZE)   
       ABSF_OLD=ABSF_NEW
      ENDIF
     ENDIF
     IF (.NOT.NEWJAC) THEN
      if (.not.allocated(fjacCOPY) .or. .not.allocated(fjacDCMP) .or. .not.allocated(fjacINDX)) &
       stop ' constant Jacobian copies not allocated '
      fjacSave=fjacCOPY      ! (used to compute the gradient, for use in lnsrch)
      FJAC=fjacDCMP          ! (used to compute p=dx in lubksb)
      INDX=fjacINDX          ! (used to compute p=dx in lubksb)
     ENDIF
   END SELECT
   !print *, 'newjac = ', newjac
   !print *, 'X = ', X
   IF (NEWJAC) THEN
    ! compute new jacobian matrix
    CALL FDJAC_ODE(X,DSDT,JAC_ODE)             ! calculate Jacobian of the ODE
    SELECT CASE(SOLUTION_METHOD)
     CASE(IMPLICIT_EULER); FJAC=-fmin_dtp*JAC_ODE    ! working towards (I - DT dg/dS), identity matrix added later
     CASE(IMPLICIT_HEUN);  FJAC=-fmin_dt2p*JAC_ODE   ! working towards (I - DT2 dg/dS), identity matrix added later
     CASE DEFAULT; STOP ' solution method muct be either implicit_euler or implicit heun '
    END SELECT
    CALL DIAGADD(FJAC,1._SP)  ! add identify matrix
    !print *, 'fjac = '; DO I=1,SIZE(X); WRITE(*,'(10(E12.5,1X))') FJAC(:,I); END DO
    !print *, 'fvec = ';  WRITE(*,'(10(E12.5,1X))') FVEC(:)
    IF (CHECK_OVERSHOOT==LINE_SEARCH) fjacSave=FJAC  ! need because FJAC overwritten in LUDCMP
    IF (JAC_RECOMPUTE==CONSTFULLSTEP .OR. JAC_RECOMPUTE==PERIOD_FREEZE .OR. JAC_RECOMPUTE==SMALL_F_RATIO) THEN
     IF (.NOT.ALLOCATED(fjacCOPY)) STOP ' constant Jacobian copies not allocated '
     fjacCOPY=FJAC  ! stored in MODULE model_numerix and re-used
    ENDIF
   ENDIF
   IF (CHECK_OVERSHOOT==LINE_SEARCH) THEN
    G=MATMUL(FVEC,fjacSave)
   ENDIF

   ! ---------------------------------------------------------------------------------------
   ! (2B) DECOMPOSE THE JACOBIAN MATRIX AND ESTIMATE DX (DX=P)
   ! ---------------------------------------------------------------------------------------
   XOLD=X
   FOLD=F
   IF (NEWJAC) THEN
    CALL LUDCMP(FJAC,INDX,D)
    IF (JAC_RECOMPUTE==CONSTFULLSTEP .OR. JAC_RECOMPUTE==PERIOD_FREEZE .OR. JAC_RECOMPUTE==SMALL_F_RATIO) THEN
     if (.not.allocated(fjacDCMP) .or. .not.allocated(fjacINDX)) &
       stop ' constant Jacobian copies not allocated '
     fjacDCMP=FJAC
     fjacINDX=INDX
    ENDIF
   ENDIF
   P=-FVEC
   CALL LUBKSB(FJAC,INDX,P)
   !print *, 'p = ', p 
 
   ! ---------------------------------------------------------------------------------------
   ! (2C) CHECK FOR OVERSHOOT AND FIX
   ! ---------------------------------------------------------------------------------------
   IF (CHECK_OVERSHOOT.EQ.LINE_SEARCH) THEN
    ! undertake line search
    CALL LNSRCH(XOLD,FOLD,G,P,X,F,STPMAX,CHECK,FMIN)
    ABSF_NEW = MAXVAL(ABS(FVEC))
    !print *, 'fvec = ', fvec, absf_new, check
    IF (ABSF_NEW < ERR_ITER_FUNC) THEN  ! test for convergence on function values
     CHECK=.FALSE.
     EXIT
    ENDIF
    IF (CHECK) THEN  ! test for a gradient of f zero (i.e., spurious convergence)
     CHECK=(MAXVAL( ABS(G)*MAX(ABS(X),1.0_SP) / MAX(F,0.5_SP*SIZE(X)) ) < TOLMIN)
     !print *, 'in check ', MAXVAL( ABS(G)*MAX(ABS(X),1.0_SP) / MAX(F,0.5_SP*SIZE(X)) ), check
     EXIT
    ENDIF
    DX = X-XOLD  ! done to account for constraints in LIMIT_XTRY (i.e., dx ne newton step)
   ELSE
    ! take full newton step
    X = XOLD + P
    CALL LIMIT_XTRY(X) ! ensure that the values of X are physically reasonable
    F = FMIN(X)        ! compute function evaluation (also populates DSDT and FVEC)
    ! test for convergence on function values
    ABSF_NEW = MAXVAL(ABS(FVEC))
    IF (ABSF_NEW < ERR_ITER_FUNC) THEN
     CHECK=.FALSE.
     EXIT
    ENDIF
    ! test that the function decreased
    IF (F.GE.FOLD) THEN
     X=XOLD
     CHECK=.TRUE.
     EXIT
    ENDIF 
    DX = X-XOLD  ! done to account for constraints in LIMIT_XTRY (i.e., dx ne newton step) 
   ENDIF
   !WRITE(*,'(I4,1X,10(E15.8,1X))') NITER, F, X
   !WRITE(*,'(I4,1X,10(E15.8,1X))') NITER, F, DX
   !WRITE(*,'(I4,1X,10(E15.8,1X))') NITER, F, ABS(FVEC)
   ! ---------------------------------------------------------------------------------------
   ! (2D) CHECK FOR CONVERGENCE
   ! ---------------------------------------------------------------------------------------
   ! check for convergence on dx
   IF (MAXVAL( ABS(DX) / MAX(ABS(X),1.0_SP) ) < ERR_ITER_DX) THEN
    CHECK=.FALSE.
    EXIT
   ENDIF
   ! check for non-convergence
   IF (ITS.EQ.NITER_TOTAL) CHECK=.TRUE.
   ! ---------------------------------------------------------------------------------------
  END DO  ! iteration loop
  ! ----------------------------------------------------------------------------------------
  END SUBROUTINE newtoniter
END MODULE newtoniter_mod
