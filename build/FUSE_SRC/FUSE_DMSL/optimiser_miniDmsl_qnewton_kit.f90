!******************************************************************
! (C) Copyright 2000-2008  ---  Dmitri Kavetski  ---  All rights reserved
!******************************************************************
module optimiser_dmsl_kit
! Purpose: Contains advanced numerical optimisation methods in Fortran-95.
! Programmer: Dmitri Kavetski, Created September 2003.
! Last modified 17 January 2005.
! This code is part of the DMSL library and is subject to use restrictions
! (see DMSL readme file for details).
! ---
! Primary References:
! * NW2000:  Nocedal,J. and Wright,S.J.(2000) Numerical Optimization, Springer.
! * F1996:   Fletcher,R.(1996) Practical Methods of Optimization,2nd Ed,Wiley.
! * DS1996:  Dennis Jr,J.E. and Schnabel,R.B.(1996) Numerical Methods for
!            Unconstrained Optimization and Nonlinear Equations, SIAM reprint.
! * GMW1981: Gill,P.E.,Murray,W. and Wright,M.H.(1981) Practical Optimization,
!            Academic Press.
! * GW1976:  Gill,P.E. and Murray,W.(1976) Minimization subject to bounds on
!            variables, NPL report NAC72.
! * P1992:   Press et al.(1992) Numerical Recipes in F-77, 2nd ed, Cambridge Press.
! ---
! Notes:
! * The module follows fairly closely the classical modular system presented by
!   Dennis and Schnabel (DS1996) also utilising material from more recent references
!   (NW2000). However, exploiting Fortran-95, the code is much more compact, at the
!   expense of some possible memory and efficiency losses (array operations).
!   These 'deficiencies' are (i) usually negligible, assuming the function
!   being optimised is expensive to evaluate, but (ii) greatly simplify developing,
!   debugging and modifying the code.
! * In addition to the 'proper' Newton-type methods (classic Newton, quasi-Newton),
!   this code also includes limited-memory quasi-Newton and conjugate-gradient methods.
!   This makes the code suitable for superdimensional problems.
! * If you function is very cheap to evaluate (i.e., cheaper than the linear algebra
!   in Newton-type equations), this code may become inefficient relative to carefull
!   micro-optimized codes, particularly for high-dimensional bounded problems.
! ---
use kinds_dmsl_kit    ! kind definitions
implicit none
private
!public::multiStart
public::qnewton,qnewtonUnwise_type
!public::LBFGS
!-----------------------------
public::QN_DMSL_mometh,LBFGS_mometh
!-----------------------------
! * Parameterised external settings
! Generic indicator: user supplied evaluator
integer(mik),parameter::user_meth=      0     ! user-provided (Hessian,gradient,etc.)
! Method selection for multistart optimisation
integer(mik),parameter::none_mometh=    1,&   ! random search
                        QN_DMSL_mometh= 2,&   ! native DMSL quasi-Newton (best for middle-D problems)
                        QN_IMSL_mometh= 3,&   ! IMSL-based quasi-Newton
                        LBFGS_mometh=   4,&   ! LBFGS scheme (best for huge-D problems)
                        SCE_mometh=     5     ! SCE search (multistart version probably redundant)
! Iteration globalisation methods
integer(mik),parameter::null_imeth=     0,&   ! No globalisation (usually 4 testing only)
                        armijo_imeth=   1,&   ! Armijo backtracking linesearch
                        wolfe_imeth=    2,&   ! Wolfe condition linesearch
                        stwolfe_imeth=  3,&   ! Strong Wolfe condition linesearch
                        brentmin_imeth= 4,&   ! Brent line minimisation
                        trustEx_imeth=  5,&   ! Near-exact trust method ("hookstep")
                        dogLeg_imeth=   6     ! Generalized dogleg trust method
! Gradient computation method
integer(mik),parameter::fd_gmeth=       1,&   ! Forward difference gradient
                        cd_gmeth=       2     ! Central difference gradient
! Hessian computation methods
integer(mik),parameter::fdg_hmeth=      1,&   ! Newton, Hessian by fwd differencing gradient
                        cdg_hmeth=      2,&   ! Newton, Hessian by cntrl differencing gradient
                        fdf_hmeth=      3,&   ! Newton, Hessian by fwd differencing function
                        cdf_hmeth=      4,&   ! Newton, Hessian by cntrl differencing function
                        bfgsInv_hmeth=  5,&   ! Quasi-Newton, BFGS update of inverse Hessian
                        bfgsUnfac_hmeth=6,&   ! Quasi-Newton, BFGS update of unfactored Hessian
                        bfgsFac_hmeth=  7,&   ! Quasi-Newton, BFGS update of factored Hessian
                        SR1unFac_hmeth= 8,&   ! Quasi-Newton, SR1 update of unfactored Hessian
                        NCG_FR_hmeth=   9,&   ! Conjugate-gradient method, Fletcher-Reeves
                        NCG_PR_hmeth=   10,&  ! Conjugate-gradient method, Polak-Ribiere
                        NCG_PPR_hmeth=  11    ! Conjugate-gradient method, Positive Polak-Ribiere
! quasi-Hessian initialisation method (ignored for non-quasi-Newton methods)
integer(mik),parameter::unt_himeth=     1,&   ! Unit matrix
                        untcnd1_himeth= 2,&   ! Unit matrix with conditioning on 1st step
                        scld_himeth=    3,&   ! Scaled matrix
                        scldcnd1_himeth=4,&   ! Scaled matrix with conditioning on 1st step
                        d2fdx2_himeth=  5,&   ! Diagonal of estimated d2f/dx2
                        hessX0_himeth=  6     ! Hessian at initial point (approx)
!-----------------------------
! * Parameterised internal settings
! Termination test values
integer(mik),parameter::no_con=         0,&   ! No convergence yet
                        grad_con=       1,&   ! Gradient criterion satisfied
                        search_con=     2,&   ! Search tolerance satisfied
                        fred_con=       3,&   ! Function reduction criterion satisfied
                        switchCD_con=   4,&   ! Switch to central difference gradient
                        srchBadGrad_con=5,&   ! Search convergent but grad still large
                        fredBadGrad_con=6     ! Function convergent but grad still large
! Globalisation return codes
integer(mik),parameter::badFunc_glob=  -3,&   ! Function evaluation returned error
                        unfeas_glob=   -2,&   ! Unfeasible points along search direction
                        badDir_glob=   -1,&   ! Bad direction provided (eg, not descent)
                        failed_glob=    1,&   ! Failed to achieve globalisation objective
                        success_glob=   0,&   ! Globalisation objective achieved
                        fconv_glob=     2     ! Globalisation converges to function precision
! Trust region iteration return codes
integer(mik),parameter::unfeas_tr=     -2,&   ! unfeasible region inside trust region
                        failed_tr=     -1,&   ! trust region did not achieve required f(x) reduction
                        suceed_tr=      0,&   ! trust region successful
                        collapsed_tr=   1,&   ! trust region collapsed to stol
                        blown_tr=       2,&   ! trust region blown up to stepmax
                        dxTiny_tr=      3,&   ! restricted step negligible
                        fconExpObs_tr=  4,&   ! expected reduction within machine precision
                        goBig_tr=       5,&   ! step to be retaken with larger trust radius
                        expRedNonP_tr=  6     ! expected function reduction nonpositive  
! Trust region subproblem outcomes
integer(mik),parameter::onTrustBound=   0,&   ! trust region step constrained by bound
                        insideTrust=    1,&   ! trust region step well inside trust radius
                        hardCase=       2,&   ! trust region encountered More's "hard case"
                        failed2Solve=  -1     ! failed to solve the trust subproblem
! Status codes for finite difference gradient estimation
integer(mik),parameter::fresh_hx=       0,&   ! freshly estimated stepsize
                        old_hx=         1     ! stepsize estimated at different point
! Active set values
integer(mik),parameter::freeVar_as=     0,&   ! free variable inside search bounds
                        loVar_as=      -1,&   ! variable fixed at lower bound
                        hiVar_as=      +1,&   ! variable fixed at upper bound
                        freeLoVar_as=  -2,&   ! variable on lower bound, Lgrng mult < 0
                        freeHiVar_as=  +2     ! variable on upper bound, Lgrng mult < 0
! Miscellaneous return codes
integer(mik),parameter::okAlg=          0,&   ! algorithmic sucess
                        failAlg=        1,&   ! algorithmic failure (not bug)
                        bugFail=       +100   ! failure due to apparent bug
! Iteration info summary
integer(mik),parameter::iterNfo_no=     0,&   ! no iter info
                        iterNfo_summ=   1,&   ! iteration summary
                        iterNfo_var=    2     ! summary and variables after each iteration
! Gradient check strategy 
integer(mik),parameter::chkG_neva=     -1,&   ! never check gradient
                        chkG_fail=      0,&   ! fast check when failed to globalise
                        chkG_f2g=       1,&   ! full check when failed to globalise
                        chkG_dxstp=     2,&   ! fast check every step with dx
                        chkG_hxstp=     3,&   ! fast check every step with hx
                        chkG_full=      4     ! full check every step
! Hessian checking strategy
integer(mik),parameter::chkHess_no=     0,&   ! no Hessian checking
                        chkHess_f2g=    1,&   ! full check when failed to globalise
                        chkHess_full=   2     ! full check every step
! Finite difference 'reliable' gradient estimation method
integer(mik),parameter::gradFD_gill=    1,&   ! Gill et al. method
                        gradFD_sw1=     2,&   ! Stepleman and Winarsky method,O(h)
                        gradFD_sw2=     3     ! Stepleman and Winarsky method,O(h2)
! Implementation of Strong Wolfe linesearch
integer(mik),parameter::strongwolfe_more=    1,&  ! Fairly sophisticated Strong Wolfe linesearch
                        strongwolfe_fletcher=2    ! Brute force "bisection"-style beast
! Elliptical scaling of Hessian
integer(mik),parameter::xscaleH_sphere= 0,&   ! Spherical Hessian
                        xscaleH_user=   1,&   ! User-supplied ellipticity (xscale)
                        xscaleH_hdiag=  2     ! Adaptive ellipticity based on Hessian diagonal
! Modified Hessian factorization
integer(mik),parameter::schnab_facmeth= 0,&   ! revised modified Cholesky-Gershgorin of Schnabel/Eskew
                        dennis_facmeth= 1     ! perturbed Cholesky-Gershgorin of Dennis/Schnabel
! Hybrid FD<->CD gradient (hybridFDCD)
integer(mik),parameter::useFDCDhybrid= -12    ! allows use of hybrid FD<->CD gradient
!                       fd_gmeth(1) and cd_gmeth(2) force strict O(1) and O(2) methods
! Types of Cauchy steps
integer(mik),parameter::cauchyInside=   0,&   ! Cauchy step inside trust region
                        cauchyOnBound=  1,&   ! Cauchy step constrained by trust
                        cauchyInfin=    2,&   ! Cauchy step wants infinity
                        cauchyZeroGrad= 3     ! Cauchy collapses because grad~0
! Active set diagonal fixing option
integer(mik),parameter::setUnit_fixDiag= 0,&  ! Set fixed diagonals to unity
                        keepDiag_fixDiag=1    ! Keep fixed diagonals as is
! Types of bounds for L-BFGS [DK: defined in LBFGS engine]
integer(mik),parameter::no_btype=        0,&  ! No bounds
                        lo_btype=        1,&  ! lower bound only
                        lh_btype=        2,&  ! low & high bounds
                        hi_btype=        3    ! high bound only
!-----------------------------
! Global constants
character(*),parameter::unknownMethodChar="UNKNOWN METHOD (USER INPUT ERROR)" ! text of error message
logical(mlk),parameter::bfgsInvNR=.true.,bfgsInvUt=bfgsInvNR  ! NR-based method for inverse quasi-Hessian update
!-----------------------------
! * External bundle
! The type below parameterises esoteric settings which are preset to default values.
! Only those users with some idea of the method (and the code!) should mess with them...
type qnewtonUnwise_type
! Initial point analysis
  real(mrk)::gtol0fac=1.e-3_mrk   ! reduction in gtol for initial point analysis
! Linesearch settings
  real(mrk)::alpha_ls=1.e-4_mrk     ! Wolfe criterion
  real(mrk)::beta_ls=0.9_mrk        ! Wolfe criterion (linesearch for Newton methods)
  real(mrk)::beta_ls_CG=1.e-3_mrk   ! Wolfe criterion (line minimisation for CG)
  integer(mik)::LNSstrongwolfe=strongwolfe_more ! strong Wolfe linesearch algorithm
  logical(mlk)::useDirDer=.false.               ! allows cheap directional derivatives (Wolfe)
  integer(mik)::linmin_ometh=2            ! line minimisation method (0=golden,1=Brent,2=dBrent)
  real(mrk)::linmin_tol=1.e-8_mrk         ! tolerance in line minimisation
  integer(mik)::linmin_itmax=1000         ! max number of iterations in line minimisation
! Trust region settings
  real(mrk)::acceptRatio_tr=1.e-4_mrk     ! acceptable fred ratio (obs/pred)
  real(mrk)::roDown_tr=0.10_mrk           ! below this fred ratio trust is decreased
  real(mrk)::radDown_tr=0.25_mrk          ! trust reduction factor
  real(mrk)::roUp_tr=0.7_mrk              ! above this fred ratio trust can be increased
  real(mrk)::stepOtrustUp_tr=0.5_mrk      ! if stepLen/trustRad>stepOtrust increase trust
  real(mrk)::radUp_tr=2.00_mrk            ! trust increase factor
  real(mrk)::trustOstepMax_tr=1.e2_mrk    ! if trustRad/stepLen>trustOstepMax truncate trust
  real(mrk)::roUpNow_tr=0.8_mrk           ! "increase trust now!" fred threshold
  integer(mik)::niter_tr=20               ! max outer iterations of trust region
  integer(mik)::ncholMax_tr=200            ! max Cholesky decomposition per trust solution
  real(mrk)::SR1forceUpdt=-1.e1_mrk       ! if SR1 perform below this ratio, force update
  logical(mlk)::pivotCholTrust=.true.     ! true for pivoted Cholesky in trust region (can be over-ruled)
  real(mrk)::dogNewtBias=0.8_mrk          ! Dogleg bias towards Newton (0=single dogleg)
  real(mrk)::boundFrac=1.0_mrk            ! prevents small trust expansions constrained by bounds
! Quasi-Hessian update settings
  logical(mlk)::skipQNupdtClassic=.false. ! forces "classic" update-skip condition in QN methods
  logical(mlk)::allowQHreset=.false.      ! reset quasi-Hessian to identity when failing
  logical(mlk)::maxSR1update=.false.      ! force frequent SR1 updates
  logical(mlk)::facBFGS_useR2=.false.     ! requests rank-2 BFGS updates (QR method)
  logical(mlk)::facBFGS_getLLt=.false.    ! DEBUG: requests backup unfactored BFGS Hessian
  logical(mlk)::dampedBFGS=.true.         ! requests damped BFGS updating (better than classic skips)
  real(mrk)::dampFac=0.2_mrk              ! BFGS damping factor
! Hessian scaling method
  integer(mik)::xscaleHmeth=xscaleH_user  ! ellipticity of Hessian
! Function roundoff estimation
  real(mrk)::Hscale=1._mrk            ! scale for roundoff estimation in f(x)
  real(mrk)::hammPow=1._mrk/3._mrk    ! power of epsRe in "h" for Hammings analysis
! Performance output
  integer(mik)::iterNfo=iterNfo_var   ! iteration info option
! Active set bound constraints handling
  real(mrk)::tolGfree_bnd=1.e-1_mrk   ! tolerance on gradient (Lgrng mult) for fast release (>1.0=>ignore)
  real(mrk)::tolOptSlack_bnd=1.e3_mrk ! slack factor on "stol&ftol" to release vars (<1.0=>ignore)
  real(mrk)::tolGfree2_bnd=1.e-1_mrk  ! tolerance for standard release (>1.0=>4 1 @ a time del)
  integer(mik)::fixDiagOption=keepDiag_fixDiag ! what to do with diagonals of fixed variables
! False convergence analysis
  real(mrk)::tolFalseDx=1.e3_mrk*epsRe    ! false convergence tolerance on dx
  integer(mik)::nFalseDxMax=100           ! max consecutive steps satisfying false tol
  integer(mik)::nFalseRfrshDxMax=20       ! max consecutive steps satisfying false tol for refresh
! Gradient checking
  integer(mik)::chkGrd=chkG_fail      ! gradient checking option
  integer(mik)::chkGrd_gmeth=fd_gmeth ! gradient checking method
  real(mrk)::chkGrd_tG=1.e-2_mrk      ! gradient check tolerance on g(x) agreement
  real(mrk)::chkGrd_tGdf=1.e-4_mrk    ! gradient check tolerance on df
  real(mrk)::chkGrd_tF=1.e-2_mrk      ! gradient check tolerance on f(x) vals
  real(mrk)::chkGrd_h=1.e0_mrk        ! h-value (scale) in gradient check
! Hessian checking
  integer(mik)::chkHess=chkHess_no        ! Hessian checking option
  integer(mik)::chkHess_hmeth=fdg_hmeth   ! Hessian checking method
  logical(mlk)::ignoreBadHess=.true.      ! no action taken on bad Hessians
! Finite difference gradient approximation
  real(mrk)::FDscale=1._mrk       ! scale for finite difference gradient (GMW,p345)
  logical(mlk)::useHxDef=.true.   ! forces default finite difference stepsize
  logical(mlk)::hybridFDCD=.false.    ! mixed FD/CD componentwise gradient approximation
  integer(mik)::dfdx0meth=gradFD_gill ! initial dfdx estimator method
  logical(mlk)::allowFDCD=.false.     ! allows enhanced switches FD<->CD gradient
  real(mrk)::tolFDCD=1.e-2_mrk    ! truncation error tolerance for FD->CD (enhanced)
  real(mrk)::fracFDCD=0.3_mrk     ! critical fraction for FD->CD switch   (enhanced)
  real(mrk)::tolCDFD=1.e+1_mrk    ! truncation error tolerance for CD->FD (enhanced)
  real(mrk)::fracCDFD=0.5_mrk     ! critical fraction for CD->FD switch   (enhanced)
  real(mrk)::tolGradFDCD=1.e-1_mrk  ! gradient tolerance for FD->CD switch
  real(mrk)::tolGradCDFD=1.e+1_mrk  ! gradient tolerance for CD->FD switch
  real(mrk)::tolDxFDCD=0.e-6_mrk    ! step tolerance for FD->CD switch
  logical(mlk)::adaptFDhX=.false.   ! adapt FD hx using Hessian diagonal
  logical(mlk)::adaptCDhX=.false.   ! adapt CD hx using Hessian diagonal
! Modified Hessian factorization settings
  integer(mik)::facmeth=schnab_facmeth ! modified factorization method
  real(mrk)::tau=undefRN          ! (schnab) these values indicate default initial e^1/3
  real(mrk)::tauBar=undefRN       ! (schnab) e^2/3. but F-95 precludes initialisation here
  real(mrk)::mu=0.1_mrk           ! (schnab) 
  real(mrk)::maxHessCond=undefRN        ! (dennis) bound on condition of modified Hessian
  logical(mlk)::controlHessCond=.false. ! Hessian condition control
endtype qnewtonUnwise_type
!-----------------------------
! * Internal data bundles
!---
type gmethBundle_type       ! - Gradient evaluation bundle
  integer(mik),pointer::gmeth_now
  logical(mlk)::useHxDef
  real(mrk)::FDscale
  real(mrk),pointer::hx(:)=>null()
  logical(mlk)::hybridFDCD
  real(mrk)::tolGradFDCD
  logical(mlk)::useDirDer
endtype gmethBundle_type
!---
type trustBundle_type       ! - Trust region bundle
  real(mrk)::acceptRatio_tr
  real(mrk)::roDown_tr
  real(mrk)::radDown_tr
  real(mrk)::roUp_tr
  real(mrk)::stepOtrustUp_tr
  real(mrk)::radUp_tr
  real(mrk)::roUpNow_tr
  real(mrk)::trustOstepMax_tr
  integer(mik)::niter_tr
  integer(mik)::ncholMax_tr
  real(mrk)::trustMax
  real(mrk)::trustMin
  real(mrk)::SR1skipTol
  real(mrk)::SR1forceUpdt
  logical(mlk)::pivotCholTrust
  real(mrk)::dogNewtBias
  real(mrk)::boundFrac
endtype trustBundle_type
!---
type objFuncBundle_type     ! - Function properties bundle
  real(mrk)::epsF
  real(mrk)::Hscale
endtype objFuncBundle_type
!---
type hessFacBundle_type     ! - Hessian factorization bundle
  integer(mik)::facmeth
  real(mrk)::tau
  real(mrk)::tauBar
  real(mrk)::mu
  real(mrk)::maxHessCond
endtype hessFacBundle_type
!-----------------------------
contains
!----------------------------------------------------
subroutine qnewton(             &
  evalFunc,dataIN,dataOUT,      & ! Objective function to be minimised
  x0,                           & ! Initial estimate of optimum
  xLo,xHi,activeSet,            & ! Upper and lower bounds on solution, active set
  gtol,stol,ftol,               & ! Termination tolerances
  xscale,fscale,fdigits,stpmax, & ! Scaling settings
  imeth,gmeth,hmeth,            & ! Computational algorithms
  himeth,trustRad,              & ! Initialisation settings
  maxIter,maxFev,               & ! Termination due to excessive effort
  uout,                         & ! Output unit for runtime information
  xopt,fopt,                    & ! Approximated optimal solution
  gradOpt,hessOpt,              & ! Useful diagnostics and information
  iter,fcalls,gcalls,hcalls,    & ! Computational cost report
  memHess2,                     & ! Memory footprint
  qnewtonUnwise,                & ! Esoteric settings that shouldnt be touched
  err,message)                    ! Return codes and runtime messages
!---------------
! Purpose: Implements Newton-type and conjugate-gradient methods for optimisation
! (minimisation) of differentiable (C2 class) functions f(x) evaluated by subroutine
! "evalFunc". Analytical derivatives (gradient and Hessian) need not be known, but
! can be approximated internally or externally.
! But if the objective function is genuinely non-differentiable, the
! Newton-type methods in this code will have (major, fatal) difficulties.
! Methods may still work on C1 functions (eg, singular Hessian), but with loss of
! computational efficiency. For C2 functions with reasonably accurate gradient estimates
! the methods generally converge at least superlinearly in the vicinity of the solution.
! Non-differentiable functions should be handled using simplex-type or Monte Carlo methods.
! ---
! Programmer: Dmitri Kavetski
! 17 January 2005.
! ---
! * For difficult problems there may be large variations in efficiency depending on
!   algorithmic settings (factors of 10-1000 not too unusual even for simple test
!   functions such as Rosenbrock with poor initial guesses/scaling). Hence, if a
!   problem is taking too long to solve, experimenting with different methods may prove
!   beneficial (may also boost confidence in results).
! * When solution bounds xLo and xHi not provided solves unconstrained problems.
!   Unconstrained algorithm may still work on "softly" constrained problems, ie,
!   when the optimum of a (possibly nonlinearly) constrained problem is "well" away
!   from the constrains. Method "should" also work on box-constrained problems
!   where the descent direction at the boundary points inwards (the box can be of
!   arbitrary shape). For truly contrained problems, supplying the optional arguments
!   "xLo,xHi,activeSet" invokes the active set strategy, which is the correct way
!   to solve bound-constrained problems. The conjugate-gradient method here also uses
!   the active set scheme, but is probably less suited to it than Newton-type methods.
!   NB: either none or both bounds (xLo and xHi) need to be supplied. If bounds are
!   supplied, then activeSet must also be supplied.
! * When solution bounds "xLo" and xHi" provided, solves bound-constrained problems
!   using the classic active set strategy. If a variable hits a bound, it is fixed
!   and its quadratic information discarded. The variable is then released only when
!   special conditions are met (based largely on its Lagrange multiplier). These
!   settings can be tuned to minimise zigzagging.
! * The Fortran-95 code here is designed for balanced clarity/efficiency. However,
!   no huge effort is taken to minimise storage (eg, 1D array storage of matrices),
!   since this negatively impacts on code readability/maintenance. In this code,
!   reduction of computations takes precedence over reduction of memory, with the
!   rationale that modern computing operates in memory-rich environment.
! * It is assumed that the primary computational cost of the optimisation is the
!   evaluation of functions. Hence no huge effort directed to optimise "small-scale"
!   arithmetic. In addition, note that the near-exact (hookstep) trust region method
!   may require several Cholesky decompositions per step. If your function is
!   considerably cheaper than a Cholesky decomposition then a dogleg-type trust region
!   (or linesearch) algorithm may be more efficient in computing time. Conjugate-gradient
!   methods are particularly fast per-step, but may be somewhat less robust.
! * Code testing and QA: The code was written from scratch by DK, and tested
!   "moderately" both internally (checking intermediate results) and externally
!   (checking performance on known test problems and also on fairly difficult problems
!   in water engineering where solutions have already been obtained using alternative
!   methods).
!   The current version has been tested using the following functions,
!   with satisfactory results in both constrained and unconstrained conditions
!     - Quadratic functions (n=2)
!     - Rosenbrock function (n=2...10,000)
!     - Powell Singular function (n=4...20) (which has singular Hessian at optimum)
!     - Trigonometric function (which has multiple 'global' optima)
!     - Helical valley function (n=3)
!     - Wood function (n=4)
!     - Water engineering objective functions, mixture of strong and mild
!       nonlinearities (n=20-80)
!   In all cases the DMSL code performed comparatively the same as equivalent IMSL
!   code (both have the same basic pseudocode, DS96). In general, it is always
!   recommended to verify numerical approximations using independent methods and codes.
!   Having a version of IMSL is often quite handy psychologically if another
!   code is having difficulties and one suspects its computer implementation.
! * Algorithm selection (Exact=exact Hessian)
!     - Typically superior options
!       . Trust region Exact/BFGS/SR1 method                (imeth=5    &  hmeth=0,6,8)
!       . Strong-Wolfe linesearch Exact/BFGS method         (imeth=3    &  hmeth=0,6-7)
!       . Strong-Wolfe linesearch with PR method (large N)  (imeth=3    &  hmeth=10,11)
!       . Brent line minimisation with PR method (large N)  (imeth=4    &  hmeth=10,11)
!     - Intrinsically incompatible options
!       . Inverse  BFGS Hessian and hookstep trust region   (imeth=5    &  hmeth=5)
!       . Factored BFGS Hessian and hookstep trust region   (imeth=5    &  hmeth=7)
!       . Conjugate gradient method and trust regions       (imeth=5-6  &  hmeth=9,10,11)
!     - Currently incompatible options
!       . Inverse  BFGS Hessian and active set method       (hmeth=5    &  xLo/xHi)
!       . Brent line minimisation and active set method     (imeth=4    &  xLo/xHi)
!     - Poor combined performance likely
!       . SR1 Hessian and linesearch methods                (imeth=1-4  &  hmeth=8)
!       . Conjugate gradient method and crude linesearches  (imeth=0-2  &  hmeth=9,10,11)
! * For problems poorly solved (or not solved at all) by the default algorithms,
!   changing the method often helps, sometimes dramatically. E.g., the trust region
!   method is generally solid, but in some cases linesearches are more efficient.
!   In very rare cases, the finite difference gradient can perform better than
!   analytical gradients (eg, on a narrow plateau analytical gradient can be
!   small, but the finite difference stepsize may indicate nearby regions of
!   larger gradient). However, the more derivative information is used, the better
!   performance is to be expected. In DK's experience, the best option is
!   exact gradient/exact Hessian with trust region. If exact gradient available,
!   estimating the Hessian from the gradient is typically the next reliable option,
!   since quasi-Newton method can be fooled by rapid changes in curvature.
! * For quasi-Newton methods, a good initial Hessian is often very beneficial,
!   and often has more influence on the final results than other settings.
! * For difficult problems the default esoteric settings can be inefficient and even
!   prevent success. The optional 'qnewtonUnwise' argument gives the user access
!   to virtually all algorithmic parameters of this code. But use them with care!
! * Whenever the exact gradient is supplied, the algorithm performance will typically
!   improve, sometimes dramatically. Note that finite difference gradients inevitably 
!   become inaccurate near optima, which limits the accuracy to which the solution
!   can be approximated. If some elements of the gradient can be calculated analytically
!   but others cannot, it may be worthwhile to use gmeth=0 option and approximate
!   the remaining elements numerically inside the user-supplied routine "evalFunc".
! * Whenever the exact Hessian is supplied, the algorithm performance will improve,
!   but often not as dramatically as when approximate gradients are replaced
!   by exact gradients. This is because the gradient determines the entirety of
!   descent directions, whereas the Hessian merely selects the optimal member of this
!   set of directions. When the gradient is approximate, the set of descent directions
!   will also be inaccurate, which more fundamentally degrades the algorithm.
! * Whenever supplying any derivatives analytically, it is strongly recommended to
!   check them as thoroughly as possible. It has been claimed (and the author
!   fully believes this) that the majority of failures when using optimisation
!   methods is due to inaccurate calculation of derivatives by the user. The code has
!   options to check the supplied gradients and Hessians. Since such checking usually
!   finds any errors very quickly, but is often expensive in terms of function calls,
!   it is usually sufficient to run the check for a few iterations only.
!   The most informative gradient checking option is "chkG_full", which checks every
!   gradient component at every step. This is very expensive and not recommended
!   except in initial stages of a project where code verification is necessary.
!   It is more efficient (but less reliable) to use the directional derivative to
!   check the gradient. Option "chkG_hxstp" check the gradient at every step using
!   a much cheaper method (2 function calls per check). By default, however, the
!   gradient is checked only when failing to globalise.
! * The code can be (and has been) used as a 'black-box' for special optimisation,
!   e.g., nonlinear least squares (NLS). However, NLS problems not only possess
!   extra structure that can be exploited for efficiency, but are also often subject
!   to strong ill-conditioning. In these cases supplying the analytical SS Hessian
!   has the negative effect of squaring the condition number of the problem.
!   For tough NLS problems a dedicated solver is often essential, based on QR or SVD
!   decomposition of Gauss-Newton Hessian equations.
! * For superdimensional problems Newton-type methods will fail due to O(N2) memory
!   growth and the O(N2)-O(N3) linear algebra on the Hessians. Instead, the
!   conjugate gradient methods should be used, which are only O(N) in both memory
!   (no Hessian stored) and cost-per-step (no linear alegbra). But there is some
!   reduction in robustness of the code, so accurate initial estimates and good
!   scaling become particularly important.
!---------------
! * Available selection of algorithms:
! 1. Iteration schemes:
!       Newton methods: Classical and Discrete
!       Quasi-Newton methods: BFGS (factored/unfactored) and SR1
!       Conjugate gradient methods: Fletcher-Reeves, Polak-Ribiere and PR+
! 2. Globalisation methods:
!       Linesearch methods
!         - Armijo, Wolfe and Strong Wolfe search conditions
!         - Brent line minimisation with/without derivatives
!       Trust region methods
!         - Near-exact (hookstep) solutions
!         - Generalized dogleg (2D subspace minimization) solutions
! 3. Active set method for bounded optimisation
! 4. Optional "smart" semi-adaptive finite difference gradient approximations:
!       Forward differences (adaptive dynamic/static stepsize)
!       Central differences (adaptive static stepsize)
! 5. Optional default-stepsize finite-difference Hessian approximations:
!       Forward differences of gradient
!       Forward differences of function
!       Central differences of function
! 6. Static/adaptive step scaling
!       Static scaling based on user xscale (Dennis and Schnabel)
!       Dynamic scaling based on Hessian ellipticity (Nocedal)
! 10.Auxiliary tools:
!       Hamming's empirical determination of function evaluation precision
!       Fast/full gradient checking
!       Full Hessian checking
!---------------
! * Code peculiarities
! 1. The statement "goto 1000 !return" is effectively a return statement - 
!   it directs the code to the cleanupMem routine to free the heap space ("hessScaled").
!---------------
! **** Troubleshooting Newton optimisation ****
! 1.  Do not expect to solve a hard problem with a default algorithm. Many problems
!     have some specific structure that makes them difficult to some methods.
!     Special insight and experimentation is often required for such problems.
!     This code offers many algorithm selections and, provided the problem has
!     a reasonably well-posed solution, it is likely one of the methods will succeed
!     and be efficient. In particular, having the analytical gradient puts at our
!     disposal far more reliable algorithms than when finite difference gradients
!     are used. Hence a day's effort in programming reliable derivatives is often
!     repaid by more consistent and efficient performance (and allows higher final accuracy).
! 2.  Algorithm returns with a warning (err<0), with the solution NOT satisfying the
!     prescribed tolerances. This typically occurs if the tolerances are too stringent and
!     iterations are stopped due to lack of progress. This warning is often innocuous and
!     in some codes is actually a sucessful return condition (suggesting convergence).
!     This code is more cautious in reporting sucess, so looking at the log file is useful
!     to establish whether a satisfactory solution is obtained. More generally, realistic
!     tolerances depend on whether analytical or approximate gradients are used. Exceedingly
!     stringent tolerances hamper efficiency. E.g., the algorithm can start thrashing around
!     at the end, since FD gradients (and often even analytical) are highly inaccurate
!     near stationary points. Even analytical gradients can be inaccurate for such points,
!     due to scaling constraints in floating point computation. Finally, Taylor series
!     analysis shows its generally impossible to zero the gradient to full machine precision.
!     See recommened ('default') tolerances below, but be prepared to modify (relax) them!
! 3.  Recommended values for primary tolerances:
!     With exact grads:  gtol=epsRe**1/2;               stol=aE*epsRe; ftol=bE*epsRe
!     With approx grads: gtol=epsRe**1/3 or epsRe**1/4; stol=aA*epsRe; ftol=bA*epsRe
!     where aE=aA=bE=bA ~ 10-100. Alternatively try something like stol=ftol ~ epsRe**2/3
!     See good discussions in GMW1981 and DS1996.
!     NB: for small-residual least-squares problems ftol can be set to epsRe**2 due
!     to special properties of these problems (see GMW1981 & DS1996).
! 4.  Finite difference gradient poorly approximating actual gradient.
!     - This algorithm can use default stepsize ("useHxDef=.true.") which assumes that
!       the function being optimised is well-scaled. This is fast and typically
!       satisfactory (at least good enough for IMSL).
!     - Poorly scaled functions may require the adaptive stepsize. "useHxDef=.false."
!       deploys stepsize adaption at the initial point and whenever slow progress
!       is taking place. This option is usually more robust, but can also be rather costly
!       in terms of function calls, particularly the SW (Stepleman/Winarsky) option.
!       A potential weakness of the "useHxDef=.false." option is that the gradient is
!       optimised for the selected points only. If the optimal stepsize varies
!       significantly, the mechanism does not always recognise that slow progress is
!       being made and hence perseveres with potentially poor gradients, slowing
!       the whole thing down. In these cases, the slow-progress diagnostics 'tolDxFDCD',
!       and 'tolFalseDx' may help.
!     - Often the Quasi-Newton Hessian gives useful order-of-magnitude estimates
!       of d2f/dx2 that can be used to cheaply optimise the gradient stepsize at each
!       Newton step ("adaptFDhX=.true.").
!     - Since forward differences become progressively inaccurate as the iterations
!       converge to an optimum, a timely switch to central differences can boost
!       efficiency by preventing many iterations at the limiting accuracy of the
!       gradient approximation, where very little progress is made at each step. This
!       algorithm will never terminate on a forward difference gradient alone, but
!       setting "allowFDCD=.true." enables additional switches to central differences,
!       which is often beneficial but sometimes wasteful if premature switches occur.
!     - The tolerance "tolGradFDCD" allows forcing central differences whenever
!       the scaled gradient is below this tolerance. "tolGradCDFD" regulates the switch
!       back to central differences, which is often beneficial when the function has
!       plateaus that must be negotiated using central differences before the area
!       of faster function variation is reached.
!       This switch is standard and is attempted regardless of allowFDCD (enhanced switches)
!     - "tolFDCD/CDFD" and "fracFDCD/CDFD" allow switches FD<->CD governed by truncation
!       error analysis based on the Hessian diagonal. This switch is non-standard and
!       seems rarely needed (a similar/more reliable effect given "tolGradFDCD/tolGradCDFD".
!       Set allowFDCD=.false. to ignore this option.
!     - "tolDxFDCD" requests a switch to central differences if the scaled step is small.
!     - "dfdx0meth" selects between the Gill et al. and Stepleman/Winarsky finite
!       differencing. The former is usually cheaper, the latter usually more robust.
!     - In general, use of finite differences gradients can degrade the convergence
!       and final accuracy of the solution. Hence it is often beneficial to take care
!       to provide the best accurate derivatives possible, preferably analytically.
!       Sometimes the user may have a better idea of stepsize selection, which can
!       certainly be very helpful. In this case, set "gmeth=user" and implement the
!       optimal finite differencing yourself.
! 5.  Zigzagging on/off bounds (this can be diagnosed from the iteration info file)
!     - For functions with many active constraints, undesirable zigzagging may occur
!       where a variable is repeatedly fixed/released, causing small steps to be taken
!       at each iteration.
!     - setting "tolOptSlack_bnd" to a low value will increase the tolerance to which
!       the function is optimised before changes in the active set are allowed.
!       Set "tolOptSlack_bnd<1" to optimise the function to full extent before releasing
!       any variables (NB: gtol is not affected by tolOptSlack_bnd, only stol&ftol).
!     - "tolGfree_bnd" specifies the critical gradient for fast active set release.
!       if a fixed variable has a positive scaled gradient in excess of the max gradient
!       times "tolGfree_bnd" then the variable is released immediately. To possibly limit
!       zigzagging, set "tolGfree_bnd>1" to disable fast release.
!     - "tolGfree2_bnd" allows to release more than 1 variable at a time when the
!       active set is forced to change. Set "tolGfree2_bnd>1" to ensure vars are released
!       one at a time. In some cases it is best to release several variables at a time.
!     - Even with optimal settings, zigzagging is still possible whenever the
!       unconstrained optimum is close to a bound. Sometimes a change of variables
!       can be beneficial to eliminate the bound or 'move' it away (eg, in log-space).
! 6.  Crazy/Small steps followed by failure/overflow/underflow
!     - Usually indicative of serious user error, eg, incorrect gradient/Hessian.
!     - It is not unheard of peoples trying to maximise instead of minimise:
!       incorrect implementation of the subroutine "evalFunc" is something to look
!       out for. If you want to maximise f(x), minimise fm=-f(x), with grad[fm]=-grad[f]
!       and hessian[fm]=-hessian[f].
!     - Small steps followed by failure are sometimes indicative of the function
!       being non-smooth, in which case the gradient may not be defined or be fairly
!       useless. Functions with higher numerical noise can sometimes be treated by
!       setting "fdigits" to less than machine precision. For truly rough functions use
!       simplex (polytope) optimisation, such as Nelder-Mead or (if fearing multimodality) SCE.
! 7.  Long cyclical iterations
!     - If very stringent accuracy is requested, the iterations may end up cycling
!       endlessly around the optimum, governed more by roundoff than by the genuine
!       function behaviour. In general set ftol~epsRe,stol~epsRe, but gtol~sqrt(epsRe)
!       and even lower whenever inexact derivatives are used.
!     - Endless cycling used to occur near bounds in older code versions.
!       This is possibly a bug in the code, so inform the author for fixing.
! 8.  Very slow progress of quasi-Newton method on unconstrained problems
!     - If no bounds present, zigzagging cannot be blamed. However, this code can
!       handle bounds 'implicitly' (do not provide xLo and xHi but use 'feas' argument
!       in evalfunc). This option should only be used if the bounds are purely safeguards
!       and the optimum is certainly to be expected well in the interior of the feasible
!       domain. Generally if you use 'feas' then it is best to supply xLo and xHi.
!     - The Newton method works best when the Hessian is reasonably-conditioned along
!       the Newton trajectory. If Hessian is ill-conditioned, two thing can occur:
!       i)  The underlying model is divergent, so that large steps are suggested and then
!           curbed by the linesearch and trust region globalisation. Although in principle
!           convergence will still occur in the limit as iter->infinity, this can become
!           impractical.
!       ii) The underlying model suggests small steps that are accepted by the
!           globalisation methods, but overall very small progress is made towards the solution.
!           Numerical evidence (Rosenbrock, 2D, x0=f*{-1.2,1},f>1000) suggests the optimisers
!           can be particularly affected by such problems when quasi-Newton Hessians are used.
!           One would expect a quasi-Newton Hessian such as BFGS or (more preferably, SR1) to
!           become progressively unreliable as the true Hessian becomes ill-conditioned.
!           Linesearches based on strong Wolfe conditions sometimes alleviate the problem,
!           since the search is allowed to expand the step length.
!     - The spherical/elliptical trust region is hard to orient in directions
!       un-aligned with coordinate axis (unless eigen-decompositions are used). Therefore
!       for rapidly curving high-dimensional valleys trust region methods may not work
!       as reliably as expected. These cases are often best handled by transforming the
!       solution space by the user so that the scaling of the variables is consistent.
!       Trying a linesearch algorithm is sometimes helpful, since it is more robust
!       wrt (affine) rescaling. Conversely, trust regions can be sensitive to poor scaling.
! 9.  Each iteration is painfully slow, even though the function is cheap to compute
!     - For very high-dimensional functions, the cost of Newton iterations can become
!       dominated by the linear algebra of Cholesky decompositions (trust region and
!       unfactored Hessians). Use of factored quasi-Hessians can then be highly beneficial
!       since it reduces the cost of each iteration from O(N3) to O(N2). Note that exact
!       trust regions cannot then be used, only linesearches (which is usually A-OK).
!     - For extremely high-dimensional functions even the O(N2) memory requirement
!       for the Hessian can be burdensome. These cases are best handled using
!       conjugate-gradient methods, which require only O(N) storage, but are typically
!       inferior to Newton methods in terms of function evaluations. Truncated Newton
!       and limited memory quasi-Newton methods are another alternatives. This code
!       implements a family of conjugate gradient methods.
! 10. Poor function scaling
!     The characteristic scale of the function, including its independent variables,
!     is critical in optimization. The scaling vector "xscale" and "fscale" should
!     be used whenever the variables are not uniformly scaled.
!     For example, if x1~1.e-10->5.e-10 but x2~1.e4->5.e4, then it is best to reflect
!     this a priori scaling information in xscale. Otherwise "mis-scaling" may occur,
!     which will be manifested in badly-conditioned Hessians and very slow progress of
!     the algorithms, especially conjugate gradient and maybe also quasi-Newton schemes.
!     Vector xscale affects the following components of all algorithms:
!     - Termination tests.
!     - Conditioning perturbation of the Hessian matrix in linesearch globalisation.
!       (except inverse BFGS updating, which is not monitored for conditioning)
!     - Trust region steps will be progressively affected by scaling, since the
!       trust region is always spherical in the scaled coordinates.
!     - Accuracy to which the Newton equations can be solved (since bad scaling
!       leads to poor conditining).
!     In addition, quasi-Newton (secant) methods are further affected as follows
!     - Initial quasi-Hessian estimate.
!     - Skipping conditions in quasi-Hessian updates.
!     In addition, finite difference versions of the algorithms are affected as follows
!     - Default stepsize estimation becomes progressively wrong.
!     Moderate mis-scaling is tolerable, but in general it is best to scale the problem
!     so that the variables are approximately of order unity in the optimal regions.
!     If variables range over many orders of magnitude nonlinear (eg, exponential)
!     transformations may be necessary.
!     By default, the Hessian 'ellipticity' scaling coincides with the user-provided
!     xscale, so that xscaleHmeth=xscaleH_user.
!     In some cases, a spherical trust region may be appropriate even though variables
!     are disparate in scale. This occurs for curving ellipsoidal objective functions
!     where the long axis is curving and an elliptical trust region is thus detrimental.
!     Then set xscaleHmeth=xscaleH_sphere. Note that this will affect Hessian scaling
!     only - all other algorithm components will still use xscale.
!     In other cases, it is possible to determine a "dimensionless" Hessian by
!     scaling by the square root of diagonal elements. This type of scaling is attractive
!     due to its invariance properties and works well for near-Gaussian objective functions,
!     where it is equivalent to converting a covariance matrix into a correlation matrix.
!     Set xscaleHmeth=xscaleH_hdiag to request this type of Hessian scaling. Again,
!     all other algorithm components will use xscale. Note that if this option is used
!     in a quasi-Newton method, default initial Hessians will not be scaled consistently,
!     which may degrade efficiency or even result in failure.
! 11. Stack/Heap/Memory overflow inside the qnewton code when using Hessian-based methods
!     This code does not exploit any sparsity of the objective function Hessian, and hence will 
!     become memory-bound for huge problems . Generally, the quasi-Newton storage constraints are
!     ndim^2, due to the storage of the Hessian matrix.
!     By design, the memory footprint of qnewton is 2*(ndim,ndim)+O(ndim). Whereas most other
!     quasi-Newton algorithms use 1*(ndim,ndim)+O(ndim), doing so significantly (IMO) increases
!     the code complexity (particularly if direction-scaling is implemented) and requires increased
!     arithmetic (since intermediate storage is unavailable).
!     This code therefore requires the storage of an additional [ndim,ndim] scaled Hessian matrix
!     ("hessScaled"), effectively halving the maximum problem memory size that could be solved if
!     all the arithmetic was rolled into the single user-provided array.
!     Note that halving the problem memory size reduces the maximum number of variables only
!     by sqrt(2) ~ 1.3, so that the actual penalty is quite small. If you really have a
!     huge problem its likely you may need to use conjugate-gradient-type algorithms.
!     Note that the [ndim,ndim] matrix is now allocated on the heap, rather than (automatically)
!     on the stack. Therefore any stack overflows are likely due to [ndim] vectors, all
!     of which are automatically allocated to the stack. Setting the stack to ~80% of available
!     RAM (~1-1.5GB on WinXP-based machines, but potentially higher for linux machines) should
!     avoid stack overflows. Of course if [ndim] vectors are too big (billions of variables!), then
!     forget about using any quasi-Newton code, and even conjugate gradients are likely
!     unfeasible. Therefore your problem would appear unsolvable using current computing
!     resources. Note that heap storage is no larger than stack storage provided the
!     later is set large enough using compiler switches.
! 12. Newton-type schemes becoming extremely inefficient for cheap super-dimensional functions.
!     Conjugate-gradient methods are ideal for very large problems since the memory footprint is
!     O(N) and computational cost per step is also O(N) - very favourable compared to Newton-type
!     methods which are O(N^2) and O(N^2)-O(N^3). Therefore, only tera-dimensional problems are
!     too big for CG methods. However CG methods can be sensitive to scaling (they are
!     related to the steepest descent scheme) and thus are not as robust for strongly nonlinear
!     problems. They work best when the Hessian eigenvalues are clustered into a few groups.
!---------------
use types_dmsl_kit,only:data_ricz_type
use utilities_dmsl_kit,only:&
  zero,half,one,&
  assertEq,assertEqLog,checkBounds,checkOnBounds,putDiag,&
  fmatmul_mv,&
  iFirstTrueLoc,getdiag,norm2,&
  getRelHxFromHx,getHxFromRelHx,&
  epsF_to_epsA,&
  getFDCDgrad,getCDgrad,&
  getHessDiagFromFunc,getHessFromGrad,getHessFromFunc
use linalg_dmsl_kit,only:choles_fwbw
implicit none
! dummies
type(data_ricz_type),intent(in),optional::dataIN
type(data_ricz_type),intent(inout),optional::dataOUT
real(mrk),intent(in)::x0(:)       ! initial estimate of solution
real(mrk),optional,intent(in)::xLo(:)             ! lower bound on solution, either none or both bounds must be present
real(mrk),optional,intent(in)::xHi(:)             ! upper bound on solution, either none or both bounds must be present
integer(mik),optional,intent(inout)::activeSet(:) ! active set (-1=lo,0=free,+1=hi), must be present if using xLo and xHi
real(mrk),intent(in)::gtol        ! scaled gradient tolerance
real(mrk),intent(in)::stol        ! scaled step tolerance
real(mrk),intent(in)::ftol        ! scaled function tolerance
real(mrk),intent(in)::xscale(:)   ! scale of independent variables
real(mrk),intent(in)::fscale      ! scale of function
integer(mik),intent(in)::fdigits  ! number of reliable digits in function evaluation (-2=estimate,-1=full machine precision)
real(mrk),intent(in)::stpmax      ! maximum scaled stepsize/trust radius (set<0 for default)
! recommended: stpmax=stmax*max(norm2(x0/xscale),norm2(one/xscale)),stmax=1.e2_mrk
integer(mik),intent(in)::maxIter  ! maximum number of iterations
integer(mik),intent(in)::maxFev   ! maximum number of function calls
integer(mik),intent(in)::imeth    ! iteration globalisation method
integer(mik),intent(in)::gmeth    ! gradient evaluation method
integer(mik),intent(in)::hmeth    ! Hessian evaluation method
integer(mik),optional,intent(in)::himeth    ! Hessian initialisation method
real(mrk),optional,intent(inout)::trustRad  ! initial scaled trust region radius (set<0 for internal default)
real(mrk),intent(out)::xopt(:)    ! optimum value of "x", for which f(x) takes its minimum value.
real(mrk),intent(out)::fopt       ! function value at optimum
integer(mik),intent(out)::iter    ! number of steps (iterations)
integer(mik),intent(out)::fcalls  ! number of function calls
integer(mik),intent(out)::gcalls  ! number of gradient calls
integer(mik),intent(out)::hcalls  ! number of Hessian calls
real(mrk),intent(inout)::gradOpt(:)             ! gradient at the optimum
real(mrk),intent(inout),optional::hessOpt(:,:)  ! Hessian at optimum
integer(mik),intent(in)::uout     ! output unit for runtime info
real(mrk),intent(out),optional::memHess2   ! additional memory necessary for allocating internal Hessian storage
type(qnewtonUnwise_type),intent(in),optional::qnewtonUnwise ! esoteric settings (use with care)
integer(mik),intent(out)::err     ! error diagnostic, err=0->ok,<0=warning,>0=error
character(*),intent(out)::message ! status description
! * user-provided function to be minimised ("objective function")
interface
  subroutine evalFunc(dataIN,dataOUT,x,feas,fx,gradFx,hessFx,err,message)
  use kinds_dmsl_kit
  use types_dmsl_kit,only:data_ricz_type
  implicit none
  type(data_ricz_type),intent(in),optional::dataIN
  type(data_ricz_type),intent(inout),optional::dataOUT
  real(mrk),intent(in)::x(:)
  logical(mlk),intent(out)::feas
  real(mrk),intent(out),optional::fx,gradFx(:),hessFx(:,:)
  integer(mik),intent(out)::err
  character(*),intent(out)::message
  endsubroutine evalfunc
endinterface
! locals
integer(mik)::ndim              ! dimensionality of objective function
real(mrk)::epsF                 ! relative function accuracy
! locals for Newton iterations
real(mrk)::stepmax              ! generic maximum step
real(mrk)::stepToBound          ! maximum step to nearest bound
real(mrk)::stepmaxL             ! "local" maximum step
real(mrk)::dx(size(x0))         ! search / shift vector
real(mrk),allocatable::hessScaled(:,:)            ! scratch Hessian workspace: heap
real(mrk)::xold(size(x0)),gradold(size(x0)),fold  ! previous iteration data
real(mrk)::fredExp,fredAct      ! expected and actual reduction in function values
! locals for modified Hessian factorization methods
real(mrk)::logdet,condEst,Einf  ! Hessian properties
integer(mik)::nfacstats(2)      ! Hessian factorization cost during step
! locals for linesearch
real(mrk)::lambda               ! steplength in linesearch
! locals for trust region
real(mrk)::trustRadTemp         ! temp copy of trust radius
logical(mlk)::trustDidGradHess  ! 
! locals for factored BFGS updating
real(mrk)::Ld(size(x0))         ! diagonal of lower Cholesky factor in factored BFGS
! locals for FD gradient
real(mrk),target::hx(size(x0))  ! stepsize for finite difference gradient
real(mrk)::d2fdx2(size(x0))     ! d2f/dx2 estimates
integer(mik),target::gmeth_now  ! current gradient evaluation method
integer(mik)::gradHx            ! FD stepsize status
real(mrk)::gOh_fdcd             ! fraction of FD gradient affected by truncation error
integer(mik)::errj(size(x0))
character(100)::messagej(size(x0))
! locals for cost reporting
integer(mik)::addFcalls,addGcalls
! locals for initial Hessian
logical(mlk)::sclHess1it
! locals for gradient checking
integer(mik)::gradCheckAnalysis
! locals for bounded search
logical(mlk)::boundedSearch     ! true if bounds supplied
logical(mlk)::hitBound          ! indicates that step truncated due to hitting bound
integer(mik)::nfree,nfree0,nfix,nthawn  ! number of variables in different status
logical(mlk)::skipDxDfCheck     ! skip convergence check (when releasing variable)
logical(mlk)::delCon            ! indicates that a constraint is to be deleted
! locals for false convergence detection
integer(mik)::nFalseDx          ! consecutive steps satisfying false convergence tolerance
! other locals
logical(mlk)::ok                ! general purpose logical
integer(mik)::retcode,globcode  ! algorithm status indicator
logical(mlk)::freshHess         ! indicates quasi-Hessian reset
! conjugate gradient variables
real(mrk)::dgg,gam,gg
! algorithm indicators
logical(mlk)::useConjGrad       ! true if conjugate gradient method in use
logical(mlk)::useQuasiHessian   ! true if quasi-Newton method in use
logical(mlk)::useTrust          ! true if trust region method in use
!----------
! 'Esoteric' parameters (can be user-defined via qnewtonUnwise_type)
! Initial point analysis
real(mrk)::gtol0fac             ! reduction in gtol for initial point analysis
! Linesearch settings
real(mrk)::alpha_ls             ! Wolfe criterion
real(mrk)::beta_ls              ! Wolfe criterion
integer(mik)::LNSstrongwolfe    ! Implementatation of Strong-Wolfe linesearch
logical(mlk)::useDirDer         ! allows use of cheap directional derivatives (Wolfe)
integer(mik)::linmin_ometh      ! line minimisation method (0=golden,1=Brent,2=dBrent)
real(mrk)::linmin_tol           ! tolerance in line minimisation
integer(mik)::linmin_itmax      ! max number of iterations in line minimisation
! Trust region settings
real(mrk)::acceptRatio_tr       ! trust region settings
real(mrk)::roDown_tr            ! below this fred ratio trust is decreased
real(mrk)::radDown_tr           ! trust reduction factor
real(mrk)::roUp_tr              ! above this fred ratio trust can be increased
real(mrk)::stepOtrustUp_tr      ! if stepLen/trustRad>stepOtrust increase trust
real(mrk)::radUp_tr             ! trust increase factor
real(mrk)::trustOstepMax_tr     ! if trustRad/stepLen>trustOstepMax truncate trust
real(mrk)::roUpNow_tr           ! "increase trust now!" threshold
integer(mik)::niter_tr          ! trust region max iterations
integer(mik)::ncholMax_tr       ! max Cholesky decomposition per trust solver
real(mrk)::SR1forceUpdt         ! if SR1 perform below this ratio, force update
logical(mlk)::pivotCholTrust    ! true for pivoted Cholesky in trust region 
real(mrk)::dogNewtBias          ! Dogleg bias towards Newton (0=single dogleg)
real(mrk)::boundFrac            ! prevents small trust expansions constrained by bounds
! Quasi-Hessian update settings
logical(mlk)::skipQNupdtClassic ! update skip condition in QN methods
logical(mlk)::allowQHreset      ! reset quasi-Hessian to identity when failing
logical(mlk)::maxSR1update      ! force frequent SR1 updates
logical(mlk)::facBFGS_useR2     ! requests rank-2 BFGS updates
logical(mlk)::facBFGS_getLLt    ! DEBUG: requests backup unfactored BFGS Hessian
logical(mlk)::dampedBFGS        ! requests damped BFGS updating
real(mrk)::dampFac              ! BFGS damping factor
! Hessian scaling method
integer(mik)::xscaleHmeth       ! ellipticity of Hessian
! Function roundoff estimation
real(mrk)::Hscale               ! scale for roundoff estimation in f(x)
real(mrk)::hammPow              ! power of epsRe in "h" for Hammings analysis
! Performance output
integer(mik)::iterNfo           ! iteration info option
! Active set bound constraints handling
real(mrk)::tolGfree_bnd         ! tolerance on gradient (Lgrng mult) for fast release
real(mrk)::tolOptSlack_bnd      ! slack factor on main "tol" to release vars
real(mrk)::tolGfree2_bnd        ! tolerance for standard release
integer(mik)::fixDiagOption     ! what to do with diagonals of fixed variables
! False convergence analysis
real(mrk)::tolFalseDx           ! false convergence tolerance on dx
integer(mik)::nFalseDxMax       ! max consecutive steps with false tol
integer(mik)::nFalseRfrshDxMax  ! max consecutive steps with false tol for refresh
! Gradient checking
integer(mik)::chkGrd            ! gradient checking option
integer(mik)::chkGrd_gmeth      ! gradient checking method
real(mrk)::chkGrd_tG            ! gradient check tolerance on g(x) agreement
real(mrk)::chkGrd_tGdf          ! gradient check tolerance on df
real(mrk)::chkGrd_tF            ! gradient check tolerance on f(x) vals
real(mrk)::chkGrd_h             ! h-value (scale) in gradient check
! Hessian checking
integer(mik)::chkHess           ! Hessian checking option
integer(mik)::chkHess_hmeth     ! Hessian checking method
logical(mlk)::ignoreBadHess     ! no action taken on bad Hessians
! Finite difference gradient approximation
real(mrk)::FDscale              ! scale for finite difference gradient (p345,GMW)
logical(mlk)::useHxDef          ! forces use of default finite difference stepsize
logical(mlk)::hybridFDCD        ! mixed hybrid FDCD componentwise gradient evaluation
integer(mik)::dfdx0meth         ! initial dfdx estimator
logical(mlk)::allowFDCD         ! allows enhanced switch FD->CD gradient
real(mrk)::tolFDCD              ! truncation error tolerance for FD->CD (enhanced)
real(mrk)::fracFDCD             ! critical fraction for FD->CD switch   (enhanced)
real(mrk)::tolCDFD              ! truncation error tolerance for CD->FD (enhanced)
real(mrk)::fracCDFD             ! critical fraction for CD->FD switch   (enhanced)
real(mrk)::tolGradFDCD          ! gradient tolerance for FD->CD switch
real(mrk)::tolGradCDFD          ! gradient tolerance for CD->FD switch
real(mrk)::tolDxFDCD            ! step tolerance for FD->CD switch
logical(mlk)::adaptFDhX         ! adapt FD hx using Hessian diagonal
logical(mlk)::adaptCDhX         ! adapt CD hx using Hessian diagonal
! Modified Hessian factorization settings
integer(mik)::facmeth           ! modified factorization method
real(mrk)::tau                  ! (schnab) these values indicate default initial e^1/3
real(mrk)::tauBar               ! (schnab) e^2/3. but F-95 precludes initialisation here
real(mrk)::mu                   ! (schnab) 
real(mrk)::maxHessCond          ! (dennis) bound on condition of modified Hessian
logical(mlk)::controlHessCond   ! Hessian condition control
!----------
! Bundles
type(gmethBundle_type)::  gmethBundle
type(trustBundle_type)::  trustBundle
type(objFuncBundle_type)::objFuncBundle
type(hessFacBundle_type)::hessFacBundle
! Default parameters hardly worth changing
real(mrk),parameter::stmax=1.e2_mrk ! used in default stepmax computation, DS96,IMSL
character(*),parameter::epsFtitle="Quasi-Newton-associated epsF estimation"
logical(mlk),parameter::useHxDefIni=.true.  ! sets stepsize in finite difference Hessians
! Debugging parameters
character(3)::facBFGS_typeH ! factored BFGS Hessian matrix storage scheme
! General purpose locals
!----------
! Start procedure here
err=0; message="qnewton/justStarted"
useConjGrad=useConjGrad_inq(hmeth)  ! method categories
useQuasiHessian=useQuasiHessian_inq(hmeth); useTrust=useTrust_inq(imeth)
! * Initialise and check input, dimensioning, etc.
if(useConjGrad)then
  call assertEq(size(x0),size(xopt),size(xscale),size(gradopt),ok,ndim)
else
  call assertEq(size(x0),size(xopt),size(xscale),size(gradopt),&
                size(hessopt,1),size(hessopt,2),ok,ndim)
endif
if(.not.ok)then                                 ! dimension error
  err=10;message="f-qnewton/dimError[mainVars]"
  call write_exitInfo(skipDetailedExitInfo=.true.)
  goto 1000 !return
endif 
call checkBounds(xLo=zero,x=xscale,chkLeq=.true.,err=err,message=message)
if(err/=0)then                                  ! check xscale is positive
  err=20;message="f-qnewton/inError/illegal[xscale<=0]/&"//message
  call write_exitInfo(skipDetailedExitInfo=.true.)
  goto 1000 !return
endif
if(fscale<=zero)then                            ! fscale must be positive
  err=21;message="f-qnewton/inError/illegal[fscale<=0]"
  call write_exitInfo(skipDetailedExitInfo=.true.)
  goto 1000 !return
elseif(gtol<zero.or.stol<zero.or.ftol<zero)then ! tolerances must be non-negative
  err=30;message="f-qnewton/inError/illegal[gtol,stol,ftol<0]"
  call write_exitInfo(skipDetailedExitInfo=.true.)
  goto 1000 !return
elseif(gtol==zero.and.stol==zero.and.ftol==zero)then ! all tolerances are zero (bad idea)
  err=31;message="f-qnewton/inError/badIdea[gtol,stol,ftol=0]"
  call write_exitInfo(skipDetailedExitInfo=.true.)
  goto 1000 !return
elseif(maxIter<1.or.maxFev<1)then               ! max cost must be positive
  err=40;message="f-qnewton/inError/illegal[maxIter,maxFev<1]"
  call write_exitInfo(skipDetailedExitInfo=.true.)
  goto 1000 !return
else
  err=0;message="qnewton/started[notFinishedYet]"
endif
call startupMem() ! generate memory space for the large array "hessScaled"
if(err/=0)return  ! memory issue - quit in a civilised manner rather than on stack overflow
call processUnwiseSettings()  ! process esoteric algorithm settings
xopt=x0     ! current (initial) point
gmeth_now=gmeth
! * check whether bounds supplied
if(present(activeSet).and.present(xLo).and.present(xHi))then ! bounds supplied correctly
  if(.not.assertEqLog(ndim,size(xLo),size(xHi),size(activeSet)))then
    err=50;message="f-qnewton/dimError[boundVars]"
    call write_exitInfo(skipDetailedExitInfo=.true.)
    goto 1000 !return
  endif                 ! check bounds for inconsistencies and violations
  call checkBounds(xLo=xLo,xHi=xHi,x=xopt,chkBnd=.true.,err=err,message=message)
  if(err/=0)then
    err=60; message="f-qnewton/&"//message
    call write_exitInfo(skipDetailedExitInfo=.true.)
    goto 1000 !return
  endif
  boundedSearch=.true.
elseif(present(activeSet).or.present(xLo).or.present(xHi))then
  err=14;message="f-qnewton/missingArg[activeSet,xLo,xHi]:checkUsage[needAll]";return
  call write_exitInfo(skipDetailedExitInfo=.true.)
  goto 1000 !return
else  ! correct specification of unbounded search
  boundedSearch=.false.
endif
call writeSettings()          ! report input settings
! * initialise active set scheme
hitBound=.false.
if(boundedSearch)then ! check for variables sitting exactly on their bounds
  call checkOnBounds(xLo=xLo,xHi=xHi,x=xopt,& !epsL=,epsH=,&
    symbI=freeVar_as,symbL=loVar_as,symbH=hiVar_as,actS=activeSet,err=err,message=message)
  if(err/=0)then
    err=70; message="f-qnewton/&"//message
    call write_exitInfo(skipDetailedExitInfo=.true.)
    goto 1000 !return
  endif
endif
fcalls=0;gcalls=0;hcalls=0;err=0;nfacstats=0;trustDidGradHess=.false.
condEst=-hugeRe; logDet=-hugeRe; Einf=-hugeRe ! flag quantities unknown b4 iterations
gOh_fdcd=-1._mrk
! * check validity of input
! - meaningless settings
if(gmeth/=user_meth)then
  selectcase(hmeth)
  case(fdg_hmeth,cdg_hmeth) ! makes no sense to finite difference approximate gradient
    err=20;message="f-qnewton/incompatible1[gmeth,hmeth]"
    call write_exitInfo(skipDetailedExitInfo=.true.)
    goto 1000 !return
  endselect
endif
! - incompatible settings
selectcase(imeth)
case(trustEx_imeth) ! hookstep pretty picky in its subsidiaries
  selectcase(hmeth)
  case(bfgsFac_hmeth) ! - can't handle factored Hessians
    err=20;message="f-qnewton/incompatible2[imeth,hmeth]"
    call write_exitInfo(skipDetailedExitInfo=.true.)
    goto 1000 !return
  case(bfgsInv_hmeth) ! - nor inverse Hessians
    err=20;message="f-qnewton/incompatible3[imeth,hmeth]"
    call write_exitInfo(skipDetailedExitInfo=.true.)
    goto 1000 !return
  endselect
endselect
if(useTrust.and.useConjGrad)then
    err=30;message="f-qnewton/incompatible3[useTrust,useConjGrad]"
    call write_exitInfo(skipDetailedExitInfo=.true.)
    goto 1000 !return
endif
! * convert number of accurate digits in function to function precision
call fdigits2epsF(fdigits,evalFunc,dataIN,dataOUT,x0,xLo,xHi,xscale,fscale,Hscale,hammPow,&
                  uout,epsFtitle,epsF,addFcalls,err,message)
fcalls=fcalls+addFcalls
if(err/=0)then
  err=20;message="f-qnewton/&"//message
  call write_exitInfo(skipDetailedExitInfo=.true.)
  goto 1000 !return
endif
! * prepare for Hessian initialisation
selectcase(hmeth)
case(bfgsInv_hmeth,bfgsUnfac_hmeth,bfgsFac_hmeth,SR1unFac_hmeth)
  sclHess1it=(himeth==untcnd1_himeth).or.(himeth==scldcnd1_himeth)
endselect
! *** Part I: Initialise function, gradient and (quasi) Hessian
selectcase(gmeth)
case(user_meth)     ! ** Analytical gradient available (best situation)
  gradHx=fresh_hx
  selectcase(hmeth)
  case(user_meth)           ! * user-supplied Hessian
    call evalFuncMacro(xx=xopt,ff=fopt,gg=gradopt,hh=hessopt,xxIsX0=.true.)
    if(err/=0)goto 1000 !return
  case(fdg_hmeth,cdg_hmeth) ! * finite difference Hessian from gradient
    call evalFuncMacro(xx=xopt,ff=fopt,gg=gradopt,xxIsX0=.true.)
    if(err/=0)goto 1000 !return
    call getHessFromGrad(evalFunc,dataIN,dataOUT,xopt,gradopt,xscale,epsF,useHxDef=useHxDefIni,&
      hmeth=hmeth,hessfd=hessopt,gcalls=addGcalls,err=retcode,message=message)
    gcalls=gcalls+addGcalls
    if(retcode/=okAlg)then
      err=-10;message="f-qnewton/&"//message
      call write_exitInfo(skipDetailedExitInfo=.true.)
      goto 1000 !return
    endif
  case(fdf_hmeth,cdf_hmeth) ! * finite difference Hessian from function
    call evalFuncMacro(xx=xopt,ff=fopt,gg=gradopt,xxIsX0=.true.)
    if(err/=0)goto 1000 !return
    call getHessFromFunc(evalFunc,dataIN,dataOUT,xopt,fopt,xscale,epsF,useHxDef=useHxDefIni,&
      hmeth=hmeth-2,hessfd=hessopt,fcalls=addFcalls,err=retcode,message=message)
    fcalls=fcalls+addFcalls
    if(retcode/=okAlg)then
      err=-10;message="f-qnewton/&"//message
      call write_exitInfo(skipDetailedExitInfo=.true.)
      goto 1000 !return
    endif
  case(bfgsInv_hmeth,bfgsUnfac_hmeth,bfgsFac_hmeth,SR1unFac_hmeth)
    call evalFuncMacro(xx=xopt,ff=fopt,gg=gradopt,xxIsX0=.true.)
    if(err/=0)goto 1000 !return
! Initialisation of quasi-Hessian. The equivalent block for FD/CD gradients
! is dispersed though the code...
    selectcase(himeth)
    case(unt_himeth,untcnd1_himeth,scld_himeth,scldcnd1_himeth)
      selectcase(hmeth)
      case(bfgsInv_hmeth)       ! * inverse-Hessian BFGS quasi-Newton
        call initQHess_inv(fopt,fscale,xscale,himeth,hessopt)
      case(bfgsUnfac_hmeth)     ! * unfactored-Hessian BFGS quasi-Newton
        call initQHess_unfac(fopt,fscale,xscale,himeth,hessopt)
      case(bfgsFac_hmeth)       ! * factored (Cholesky) Hessian BFGS quasi-Newton
        call initQHess_fac(fopt,fscale,xscale,himeth,hessopt,Ld,facBFGS_getLLt)
      case(SR1unFac_hmeth)      ! * unfactored SR1 quasi-Newton
        call initQHess_unfac(fopt,fscale,xscale,himeth,hessopt)
      endselect
    case(d2fdx2_himeth) ! * generate Hessian diagonal from function
! When the gradient is obtained from finite differences this section is not needed
      call getHessDiagFromFunc(evalFunc,dataIN,dataOUT,xopt,fopt,xscale,epsF,useHxDef=useHxDefIni,&
        hmeth=1,hessDiag=d2fdx2,fcalls=addFcalls,err=retcode,message=message)
      fcalls=fcalls+addFcalls
      if(retcode/=okAlg)then
        err=-10;message="f-qnewton/&"//message
        call write_exitInfo(skipDetailedExitInfo=.true.)
        goto 1000 !return
      endif
      call makeGoodHessDiag(d2fdx2,controlHessCond,maxHessCond)
      if(hmeth==bfgsFac_hmeth)Ld=d2fdx2
      selectcase(hmeth)
      case(bfgsInv_hmeth)       ! * inverse-Hessian BFGS quasi-Newton
        call initQHess_inv(fopt,fscale,xscale,himeth,hessopt,d2fdx2)
      case(bfgsUnfac_hmeth)     ! * unfactored-Hessian BFGS quasi-Newton
        call initQHess_unfac(fopt,fscale,xscale,himeth,hessopt,d2fdx2)
      case(bfgsFac_hmeth)       ! * factored (Cholesky) Hessian BFGS quasi-Newton
        call initQHess_fac(fopt,fscale,xscale,himeth,hessopt,Ld,facBFGS_getLLt)
      case(SR1unFac_hmeth)      ! * unfactored SR1 quasi-Newton
        call initQHess_unfac(fopt,fscale,xscale,himeth,hessopt,d2fdx2)
      endselect
    case(hessX0_himeth)  ! * initialise Hessian by finite differencing the gradient
      call getHessFromGrad(evalFunc,dataIN,dataOUT,xopt,gradopt,xscale,epsF,useHxDef=useHxDefIni,&
        hmeth=1,hessfd=hessopt,gcalls=addGcalls,err=retcode,message=message)
      gcalls=gcalls+addGcalls
      if(retcode/=okAlg)then
        err=-10;message="f-qnewton/&"//message
        call write_exitInfo(skipDetailedExitInfo=.true.)
        goto 1000 !return
      endif
      selectcase(hmeth)
      case(bfgsInv_hmeth,bfgsFac_hmeth) ! currently unsupported since requires early Cholesky
        err=-10;message="f-qnewton/unsupported[himeth=hessX0_himeth]"
        call write_exitInfo(skipDetailedExitInfo=.true.)
        goto 1000 !return
      endselect
    endselect
  case(NCG_FR_hmeth,NCG_PR_hmeth,NCG_PPR_hmeth) ! * Conjugate gradient methods
    call evalFuncMacro(xx=xopt,ff=fopt,gg=gradopt,xxIsX0=.true.)
    if(err/=0)goto 1000 !return
  case default              ! * unrecognized Hessian method
    err=+10;message="f-qnewton/unknown[hmeth]"
    call write_exitInfo(skipDetailedExitInfo=.true.)
    goto 1000 !return
  endselect
case(fd_gmeth,cd_gmeth)    ! ** Finite difference gradient
  hx=sqrt(epsF)*max(abs(x0),xscale)  ! default stepsize
  selectcase(hmeth)
  case(user_meth)               ! * user-supplied Hessian
    call evalFuncMacro(xx=xopt,ff=fopt,hh=hessopt,xxIsX0=.true.)
    if(err/=0)goto 1000 !return
    call getHx_macro()  ! optimise FD gradient stepsize
    if(err/=0)goto 1000 !return
  case(fdg_hmeth,cdg_hmeth)     ! * finite difference Hessian from gradient
    err=+100; message="f-qnewton/bug/invalidInputNotCaught[fdg_hmeth,cdg_hmeth]"
    call write_exitInfo(skipDetailedExitInfo=.true.)
    goto 1000 !return
  case(fdf_hmeth,cdf_hmeth)     ! * finite difference Hessian from function
    call evalFuncMacro(xx=xopt,ff=fopt,xxIsX0=.true.)
    if(err/=0)goto 1000 !return
    call getHessFromFunc(evalFunc,dataIN,dataOUT,xopt,fopt,xscale,epsF,useHxDef=useHxDefIni,&
      hmeth=hmeth-2,hessfd=hessopt,fcalls=addFcalls,err=retcode,message=message)
    fcalls=fcalls+addFcalls
    if(retcode/=okAlg)then
      err=-10;message="f-qnewton/&"//message
      call write_exitInfo(skipDetailedExitInfo=.true.)
      goto 1000 !return
    endif
    call getHx_macro()  ! optimise FD gradient stepsize
    if(err/=0)goto 1000 !return
  case(bfgsInv_hmeth)           ! * inverse-Hessian BFGS quasi-Newton
    call evalFuncMacro(xx=xopt,ff=fopt,xxIsX0=.true.)
    if(err/=0)goto 1000 !return
    call getHx_macro()  ! optimise FD gradient stepsize
    if(err/=0)goto 1000 !return
    selectcase(himeth)  ! initialise inverse Hessian to "standard" options
    case(unt_himeth,untcnd1_himeth,scld_himeth,scldcnd1_himeth)
      call initQHess_inv(fopt,fscale,xscale,himeth,hessopt)
    case(d2fdx2_himeth) ! use estimated d2f/dx2 for a diagonal Hessian
      call makeGoodHessDiag(d2fdx2,controlHessCond,maxHessCond)
      call initQHess_inv(fopt,fscale,xscale,himeth,hessopt,d2fdx2)
    case(hessX0_himeth) ! initialise Hessian by differencing function
      err=-10;message="f-qnewton/unsupported[himeth=hessX0_himeth]"
      call write_exitInfo(skipDetailedExitInfo=.true.) ! not supported since requires inversion of dense matrix
      goto 1000 !return
    endselect
  case(bfgsUnfac_hmeth)         ! * unfactored-Hessian BFGS quasi-Newton
    call evalFuncMacro(xx=xopt,ff=fopt,xxIsX0=.true.)
    if(err/=0)goto 1000 !return
    call getHx_macro()  ! optimise FD gradient stepsize
    if(err/=0)goto 1000 !return
    selectcase(himeth)  ! initialise Hessian to "standard" options
    case(unt_himeth,untcnd1_himeth,scld_himeth,scldcnd1_himeth)
      call initQHess_unfac(fopt,fscale,xscale,himeth,hessopt)
    case(d2fdx2_himeth) ! use estimated d2f/dx2 for a diagonal Hessian
      call makeGoodHessDiag(d2fdx2,controlHessCond,maxHessCond)
      call initQHess_unfac(fopt,fscale,xscale,himeth,hessopt,d2fdx2)
    case(hessX0_himeth) ! initialise Hessian by differencing function
      call getHessFromFunc(evalFunc,dataIN,dataOUT,xopt,fopt,xscale,epsF,useHxDef=useHxDefIni,&
        hmeth=1,hessfd=hessopt,fcalls=addFcalls,err=retcode,message=message)
      fcalls=fcalls+addFcalls
      if(retcode/=okAlg)then
        err=-10;message="f-qnewton/&"//message
        call write_exitInfo(skipDetailedExitInfo=.true.)
        goto 1000 !return
      endif
    endselect
  case(bfgsFac_hmeth)           ! * factored (Cholesky) Hessian BFGS quasi-Newton
    call evalFuncMacro(xx=xopt,ff=fopt,xxIsX0=.true.)
    if(err/=0)goto 1000 !return
    call getHx_macro()
    if(err/=0)goto 1000 !return
    selectcase(himeth)  ! initialise Hessian to "standard" options
    case(unt_himeth,untcnd1_himeth,scld_himeth,scldcnd1_himeth)
      call initQHess_fac(fopt,fscale,xscale,himeth,hessopt,Ld,facBFGS_getLLt)
    case(d2fdx2_himeth) ! use estimated d2f/dx2 for a diagonal Hessian
      call makeGoodHessDiag(d2fdx2,controlHessCond,maxHessCond)
      Ld=d2fdx2
      call initQHess_fac(fopt,fscale,xscale,himeth,hessopt,Ld,facBFGS_getLLt)
    case(hessX0_himeth) ! initialise Hessian by differencing function
      err=-10;message="f-qnewton/unsupported[himeth=hessX0_himeth]"
      call write_exitInfo(skipDetailedExitInfo=.true.) ! currently unsupported since requires early Cholesky
      goto 1000 !return
    endselect
  case(SR1unFac_hmeth)          ! * unfactored SR1 quasi-Newton
    call evalFuncMacro(xx=xopt,ff=fopt,xxIsX0=.true.)
    if(err/=0)goto 1000 !return
    call getHx_macro()
    if(err/=0)goto 1000 !return
    selectcase(himeth)  ! initialise Hessian to "standard" options
    case(unt_himeth,untcnd1_himeth,scld_himeth,scldcnd1_himeth)
      call initQHess_unfac(fopt,fscale,xscale,himeth,hessopt)
    case(d2fdx2_himeth) ! use estimated d2f/dx2 for a diagonal Hessian
      call makeGoodHessDiag(d2fdx2,controlHessCond,maxHessCond)
      call initQHess_unfac(fopt,fscale,xscale,himeth,hessopt,d2fdx2)
    case(hessX0_himeth) ! initialise Hessian by differencing function
      call getHessFromFunc(evalFunc,dataIN,dataOUT,xopt,fopt,xscale,epsF,useHxDef=useHxDefIni,&
        hmeth=1,hessfd=hessopt,fcalls=addFcalls,err=retcode,message=message)
      fcalls=fcalls+addFcalls
      if(retcode/=okAlg)then
        err=-10;message="f-qnewton/&"//message
        call write_exitInfo(skipDetailedExitInfo=.true.)
        goto 1000 !return
      endif
    endselect
  case(NCG_FR_hmeth,NCG_PR_hmeth,NCG_PPR_hmeth) ! * Conjugate gradient methods
    call evalFuncMacro(xx=xopt,ff=fopt,xxIsX0=.true.)
    if(err/=0)goto 1000 !return
    call getHx_macro()
    if(err/=0)goto 1000 !return
  case default
    err=+10;message="f-qnewton/unknown[hmeth]"
    call write_exitInfo(skipDetailedExitInfo=.true.)
    goto 1000 !return
  endselect
! convert stepsize into relative interval
  hx=getRelHxFromHx(hx,xopt,xscale,FDscale)
  gradHx=fresh_hx ! freshly estimated stepsize
case default
  err=+20;message="f-qnewton/unknown[gmeth]"
  call write_exitInfo(skipDetailedExitInfo=.true.)
  goto 1000 !return
endselect
if(useConjGrad)then ! Initialisation of conjugate gradient methods
  Ld=-gradOpt !/xscale ! start with steepest descent direction
endif
!------------
! check stepmax
call checkStepmax_macro()
! Initialise globalisation strategy
if(useTrust)then ! trust-region methods require trust radius
  if(present(trustRad))then
    if(trustRad<=zero)then  ! - initialise trust region radius to Cauchy step
      call setTrustToCauchy(hess=hessopt,hessScaled=hessScaled,grad=gradopt,&
        xscale=xscale,trustRad=trustRad,steepStepLen=trustRadTemp,&
        err=err,message=message)
      if(err/=0)then
        err=+30;message="f-qnewton/BUG?/&"//message
        call write_exitInfo(skipDetailedExitInfo=.true.)
        goto 1000 !return
      elseif(trustRad>stepmax)then
        trustRad=stepmax    ! limit radius by stepmax
      endif
    endif ! else            ! - user-supplied (positive) value
  else
    err=+10;message="f-qnewton/missing[trust]"
    call write_exitInfo(skipDetailedExitInfo=.true.)
    goto 1000 !return
  endif
endif
! * All algorithmic constant initialization should be finished
! - prepare bundles
call makeGmethBundle()                        ! bundle for gradient evaluation
if(useTrust)call makeTrustBundle()            ! bundle for trust evaluation
call makeObjFuncBundle()                      ! bundle for function properties
if(.not.useConjGrad)call makeHessFacBundle()  ! bundle of Hessian factorization settings
iter=0  ! define 'iter' for output of routines below
! check gradient?
selectcase(chkGrd)
case(chkG_hxstp,chkG_full)
  call checkGrad_macro()
  if(err>0)then
    call write_exitInfo(skipDetailedExitInfo=.true.)
    goto 1000 !return
  endif
endselect
! check Hessian?
selectcase(chkHess)
case(chkHess_full)
  call checkHess_macro()
  if(err>0)then
    call write_exitInfo(skipDetailedExitInfo=.true.)
    goto 1000 !return
  endif
endselect
! check whether initial point is on any bound
nfree=ndim;nfree0=ndim;nfix=0;nthawn=0
if(boundedSearch)then
  call checkActiveSet(xopt,xscale,activeSet,fixDiagOption,hmeth,gradOpt,hessOpt,&
                      Ld,xLo,xHi,hitBound,nfree,nfix,nthawn)
  if(nfix==ndim)then    ! all variables fixed
    err=0;message="qnewton/problem/onBound&fixed[x0]"
    call write_exitInfo(skipDetailedExitInfo=.true.)
    goto 1000 !return
  else                  ! check for constraints to be deleted
    call checkReleaseActiveSet(xopt,xscale,activeSet,gradopt,hessopt,Ld,&
       tolGfree_bnd,nfree==0,tolGfree2_bnd,nfree,nfix,nthawn)
  endif
endif
iter=0;dx=zero;fredAct=zero;fredExp=zero;freshHess=.true.
retcode=0;globcode=0;nFalseDx=0;delCon=.false.
call write_iterationInfo()   ! write initial info to output file as required
! check initial point using stringent convergence
call checkConvergence0(xopt,fopt,gradOpt,activeSet,xscale,fscale,gtol0fac*gtol,retcode)
selectcase(retcode)
case(no_con)
case(grad_con)
  err=0;message="qnewton/grad[x0]~0"
  call write_exitInfo(skipDetailedExitInfo=.true.)
  goto 1000 !return
endselect
! *** Part II. Iteration loop
do iter=1, maxiter
! * Get optimal local descent direction using quasi-newton eqn
  do  ! allow several attempts when using finite difference gradient
    selectcase(imeth)
    case(null_imeth,armijo_imeth,wolfe_imeth,stwolfe_imeth,brentmin_imeth)
      selectcase(hmeth)             ! * linesearch methods
      case(user_meth,&    ! construct and solve modified Newton equations
           fdg_hmeth,cdg_hmeth,fdf_hmeth,cdf_hmeth,&
           bfgsUnfac_hmeth,SR1unFac_hmeth)
        call solveModNewtHess(hess=hessopt,hessScaled=hessScaled,Ld=Ld,grad=gradopt,&
          hessFacBundle=hessFacBundle,xscaleHmeth=xscaleHmeth,xscale=xscale,fscale=fscale,&
          activeset=activeset,dx=dx,ncholstats=nfacstats,&
          logdet=logdet,condest=condest,Einf=Einf,err=err,message=message)
        if(err/=okAlg)then
          message="f-qnewton/&"//message
          call write_exitInfo()
          goto 1000 !return
        endif
      case(bfgsFac_hmeth) ! positive definite Cholesky decomposition already in-place
        call choles_fwbw(a=HessOpt,Ld=Ld,b=gradOpt,usePivot=.false.,x=dx,err=retcode,message=message)
        if(retcode/=okAlg)then
          err=-30;message="f-qnewton/&"//message
          call write_exitInfo()
          goto 1000 !return
        endif
        dx=-dx
      case(bfgsInv_hmeth) ! inverse Hessian available directly
        if(bfgsInvUt)then   ! work with upper triangle only
          dx=-fmatmul_mv(m=HessOpt,v=gradOpt,typeMV="SUV")
        else                ! full matmul
          dx=-matmul(HessOpt,gradOpt)
        endif
      case(NCG_FR_hmeth,NCG_PR_hmeth,NCG_PPR_hmeth) ! * conjugate gradient methods
        dx=Ld
      endselect
    case(trustEx_imeth,dogLeg_imeth) ! * trust methods (implemented in dedicated sub)
    case default
      err=10;message="f-qnewton/unknownIMETH"
      call write_exitInfo()
      goto 1000 !return
    endselect
! * Standard ("full") Newton step determined.
    xold=xopt; fold=fopt; gradold=gradopt ! save old location (nb: conjugate gradient uses this info)
    selectcase(imeth)                 ! - linesearch methods: check dx for bounds...
    case(null_imeth,armijo_imeth,wolfe_imeth,stwolfe_imeth,brentmin_imeth)
      if(boundedSearch)then ! ... and truncate if violating 'em.
        call checkStepBounds(xopt,xLo,xHi,activeSet,dx,stepToBound)
        stepmaxL=getStepLen2(dx*stepToBound,xscale)
        stepmaxL=min(stepmax,stepmaxL)
      endif
    case(trustEx_imeth,dogLeg_imeth)  ! - trust region checks for bounds internally
    endselect
! * Globalisation strategy (largely independent of Hessian method)
    selectcase(imeth)
    case(null_imeth)                  ! * Testing only: no globalisation
      xopt=xopt+dx; globcode=success_glob
      call getfredExp() ! expected function reduction based on step dx
      call evalFunc(dataIN,dataOUT,xopt,ok,fopt,err=err,message=message)
      fcalls=fcalls+1
      if(err/=0)then
        write(message,'(a,i0,a)')"f-qnewton/userErr[iter=",iter,"]/&"//trim(message)
        globcode=badFunc_glob
      elseif(.not.ok)then
        write(message,'(a,i0,a)')"f-qnewton/userUnfeas[iter=",iter,"]/&"//trim(message)
        globcode=unfeas_glob
      endif
      fredAct=fold-fopt   !   actual reduction in function value
    case(armijo_imeth)                ! * Armijo backtracking
      lambda=one                      !   always start with natural Newton step
      call getfredExp() ! expected function reduction based on step dx
      call linesearch_armijo(evalFunc,dataIN,dataOUT,xold,fold,gradopt,dx,xscale,&
        stol,alpha_ls,stepmaxL,xopt,fopt,fredAct,lambda,addFcalls,globcode,message)
      fcalls=fcalls+addFcalls
    case(wolfe_imeth)                 ! * Wolfe linesearch
      lambda=one                      !   always start with natural Newton step
      call getfredExp() ! expected function reduction based on step dx
      call linesearch_wolfe(evalFunc,dataIN,dataOUT,xold,fold,gradold,gmethBundle,objFuncBundle,&
        dx,xscale,fscale,stol,alpha_ls,beta_ls,stepmaxL,xopt,fopt,gradopt,fredAct,&
        lambda,addFcalls,addGcalls,globcode,message)
      fcalls=fcalls+addFcalls; gcalls=gcalls+addGcalls
    case(stwolfe_imeth)               ! * Strong Wolfe linesearch
      lambda=one                      !   always start with natural Newton step
      call getfredExp() ! expected function reduction based on step dx
      call linesearch_strongwolfe(evalFunc,dataIN,dataOUT,xold,fold,gradold,gmethBundle,objFuncBundle,&
        dx,xscale,fscale,stol,alpha_ls,beta_ls,stepmaxL,LNSstrongwolfe,&
        xopt,fopt,gradopt,fredAct,lambda,addFcalls,addGcalls,globcode,message)
      fcalls=fcalls+addFcalls; gcalls=gcalls+addGcalls
    case(brentmin_imeth)              ! * Brent line minimisation
      xopt=xold; lambda=one
      call getfredExp() ! expected function reduction based on step dx
      call brentmin(evalFunc,dataIN,dataOUT,linmin_ometh,xopt,fold,dx,stepmaxL,stol,&
        linmin_tol,linmin_itmax,xscale,fopt,lambda,addFcalls,addGcalls,globcode,message)
      fcalls=fcalls+addFcalls; gcalls=gcalls+addGcalls; fredAct=fold-fopt
    case(trustEx_imeth,dogLeg_imeth)  ! * Trust region globalization
      trustRadTemp=trustRad   ! save trust in case step needs retaken with better grad
      call trustDriver(evalFunc,dataIN,dataOUT,x0=xold,fx0=fold,grad0=gradold,&
        hess0=hessopt,Ld0=Ld,hessScaled=hessScaled,&
        boundedSearch=boundedSearch,xLo=xLo,xHi=xHi,activeSet=activeSet,&
        imeth=imeth,hmeth=hmeth,&
        quadTypeH=facBFGS_typeH,maxSR1update=maxSR1update,gmethBundle=gmethBundle,&
        xscaleHmeth=xscaleHmeth,xscale=xscale,fscale=fscale,&
        trustBundle=trustBundle,objFuncBundle=objFuncBundle,&
        hessFacBundle=hessFacBundle,didGradNewHess=trustDidGradHess,&
        x=xopt,fx=fopt,gradx=gradopt,dx=dx,trustRad=trustRad,&
        fredExp=fredExp,fredAct=fredAct,&
        fcalls=addFcalls,gcalls=addGcalls,ncholstats=nfacstats,&
        logdet=logdet,condest=condest,Einf=Einf,err=globcode,message=message)
      fcalls=fcalls+addFcalls; gcalls=gcalls+addGcalls
    case default
      err=+10;write(message,'(a,i0,a)')"f-qnewton/unknown[imeth=",imeth,"]"
      call write_exitInfo()
      goto 1000 !return
    endselect
!    call write_iterationInfo()
    if(fredAct<zero.and.imeth/=null_imeth)then
      err=100;message="f-qnewton/BUG/fredAct<0/&"//message
      call write_exitInfo()
      goto 1000 !return
    endif
! * Take action depending on globalisation strategy result
    selectcase(globcode)
    case(bugFail)       ! - found apparent bug in code
      err=100;message="f-qnewton/BUG/&"//message
      call write_exitInfo()
      goto 1000 !return
    case(badFunc_glob)  ! - error in user function
      err=150;message="f-qnewton/USERERROR/&"//message
      call write_exitInfo()
      goto 1000 !return
    case(unfeas_glob)   ! - could not find globalised feasible point
      err=20            !   satisfying necessary conditions
      selectcase(imeth)
      case(null_imeth,armijo_imeth,wolfe_imeth,stwolfe_imeth,brentmin_imeth)
        message="f-qnewton/searchDirUnfeas/&"//message
      case(trustEx_imeth,dogLeg_imeth)
        message="f-qnewton/couldNotFindFeasSearchDir/&"//message
      endselect
      call write_exitInfo()
      goto 1000 !return
    case(badDir_glob)   ! - search direction not descent
      err=30;
      if(useConjGrad)then
        message="f-qnewton/searchDirBad/checkConvergence/&"//message
      else
        message="f-qnewton/searchDirBad/bug??/&"//message
      endif
      call write_exitInfo()
      goto 1000 !return
    case(failed_glob)   ! - failed globalisation
      if(uout>0)write(uout,'(a,i7,a)')"iter=",iter,"  Globalization FAILED ..."//trim(message)
      call checkGrad_macro()  ! check gradient accuracy
      if(err>0)then
        call write_exitInfo()
        goto 1000 !return
      endif
      call checkHess_macro()  ! check gradient accuracy
      if(err>0)then
        call write_exitInfo()
        goto 1000 !return
      endif
      selectcase(gmeth_now) ! * forward difference gradient
      case(fd_gmeth)        ! try refreshing stepsize
        if(useHxDef)gradHx=fresh_hx  ! botch to switch to CD immediately
        selectcase(gradHx)
        case(old_hx)
          hx=getHxFromRelHx(hx,xopt,xscale,FDscale)
          call getHx_macro()    ! re-optimise stepsize
          if(err/=0)goto 1000 !return
          hx=getRelHxFromHx(hx,xopt,xscale,FDscale)
          if(any(errj/=0))then  ! perceived errors not unusual when grad[fx] is small
            call write_FDCDswitchInfo(fd_to_cd=.true.)
            gmeth_now=cd_gmeth  ! switch to central differences immediately
          endif
          gradHx=fresh_hx
          trustRad=trustRadTemp ! restore trust region
        case(fresh_hx)      ! switch to central differences
          gmeth_now=cd_gmeth
          hx=getHxFromRelHx(hx,xopt,xscale,FDscale)
          call getHx_macro()    ! re-optimise stepsize
          if(err/=0)goto 1000 !return
          hx=getRelHxFromHx(hx,xopt,xscale,FDscale)
          call write_FDCDswitchInfo(fd_to_cd=.true.)
          gradHx=fresh_hx
          trustRad=trustRadTemp ! restore trust region
        case default
          exit
        endselect ! refreshed gradient may alter the active set
        if(boundedSearch)then
          call checkActiveSet(xopt,xscale,activeSet,fixDiagOption,hmeth,&
                              gradopt,hessopt,Ld,xLo,xHi,&
            hitBound,nfree,nfix,nthawn)
          if(nfix==ndim)then    ! all dimensions fixed
            exit
          else                  ! check for constraints to be deleted
            call checkReleaseActiveSet(xopt,xscale,activeSet,gradopt,hessopt,Ld,&
               tolGfree_bnd,nfree==0,tolGfree2_bnd,nfree,nfix,nthawn)
          endif
        endif
        if(useQuasiHessian.and.allowQHreset)call setUnitQhess_macro()
      case(cd_gmeth)        ! * central difference gradient
        if(useHxDef)gradHx=fresh_hx  ! botch to exit immediately
        selectcase(gradHx)
        case(old_hx)
          hx=getHxFromRelHx(hx,xopt,xscale,FDscale)
          call getHx_macro()    ! re-optimise stepsize
          if(err/=0)goto 1000 !return
          hx=getRelHxFromHx(hx,xopt,xscale,FDscale)
          gradHx=fresh_hx
          trustRad=trustRadTemp ! restore trust region
        case(fresh_hx)      ! already tried refreshing central differences
          exit              ! not much else can be done...
        case default
          exit
        endselect ! refreshed gradient may alter the active set
        if(boundedSearch)then
          call checkActiveSet(xopt,xscale,activeSet,fixDiagOption,hmeth,&
                              gradopt,hessopt,Ld,xLo,xHi,&
            hitBound,nfree,nfix,nthawn)
          if(nfix==ndim)then    ! all dimensions fixed
            exit
          else                  ! check for constraints to be deleted
            call checkReleaseActiveSet(xopt,xscale,activeSet,gradopt,hessopt,Ld,&
               tolGfree_bnd,nfree==0,tolGfree2_bnd,nfree,nfix,nthawn)
          endif
        endif
        if(useQuasiHessian.and.allowQHreset)call setUnitQhess_macro()
      case(user_meth)       ! * analytical gradient was used... cant do much else
        if(freshHess)then
          exit
        elseif(useQuasiHessian.and.allowQHreset)then   ! ..except try resetting quasi-Hessian
          call setUnitQhess_macro()
          trustRad=trustRadTemp ! restore trust region
          freshHess=.true.
        else
          exit
        endif
      endselect
    case(success_glob)  ! - succeeded in globalisation strategy
!      gradHx=old_hx; freshHess=.false. ! these now set after convergence test
      exit
    case default
      err=10;message="f-qnewton/unknown/BUG?/&"//message
      call write_exitInfo()
      goto 1000 !return
    endselect
  enddo
! * Evaluate required derivatives (gradient/Hessian) at "globalised" point
  selectcase(hmeth)
  case(user_meth)             ! ** user-supplied Hessian
    selectcase(imeth)
    case(null_imeth,armijo_imeth,brentmin_imeth,trustEx_imeth,dogLeg_imeth)  ! need gradient/Hessian
      selectcase(gmeth_now)
      case(user_meth)     ! - analytical gradient
        call evalFuncMacro(xx=xopt,gg=gradopt,hh=hessopt,xxIsX0=.false.)
        if(err/=0)goto 1000 !return
      case(fd_gmeth)      ! - FD gradient
        call evalFuncMacro(xx=xopt,hh=hessopt,xxIsX0=.false.)
        if(err/=0)goto 1000 !return
        if(adaptFDhX)then
          call adaptFDgradHx(hx,xopt,xscale,FDscale,&
                             epsF_to_epsA(epsF,fopt,fscale,Hscale),hessopt)
        endif
        call getFDCDmacro()
      case(cd_gmeth)      ! - CD gradient
        call evalFuncMacro(xx=xopt,hh=hessopt,xxIsX0=.false.)
        if(err/=0)goto 1000 !return
        call getCDmacro()
      endselect
    case(wolfe_imeth,stwolfe_imeth)  ! just need Hessian
        call evalFuncMacro(xx=xopt,hh=hessopt,xxIsX0=.false.)
        if(err/=0)goto 1000 !return
    endselect
  case(fdg_hmeth,cdg_hmeth)   ! ** finite difference Hessian from gradient
    selectcase(imeth)
    case(null_imeth,armijo_imeth,brentmin_imeth,trustEx_imeth,dogLeg_imeth)  ! need gradient
      selectcase(gmeth_now)
      case(user_meth)     ! - analytical gradient
        call evalFuncMacro(xx=xopt,gg=gradopt,xxIsX0=.false.)
        if(err/=0)goto 1000 !return
      case default
        err=200;message="f-qnewton/BUG/shouldntBeHere:hmeth=fdg&gmeth/=user"
        goto 1000 !return
      endselect
    case(wolfe_imeth,stwolfe_imeth)  ! gradient already available
    endselect
    call getHessFromGrad(evalFunc,dataIN,dataOUT,xopt,gradopt,xscale,epsF,useHxDef=useHxDefIni,&
      hmeth=hmeth,hessfd=hessopt,gcalls=addGcalls,err=retcode,message=message)
    gcalls=gcalls+addGcalls
    if(retcode/=okAlg)then
      err=+10;message="f-qnewton/&"//message
      call write_exitInfo()
      goto 1000 !return
    endif
  case(fdf_hmeth,cdf_hmeth)   ! ** finite difference Hessian from function
    selectcase(imeth)
    case(null_imeth,armijo_imeth,brentmin_imeth,trustEx_imeth,dogLeg_imeth)  ! need gradient
      selectcase(gmeth_now)
      case(user_meth)     ! - analytical gradient
        call evalFuncMacro(xx=xopt,gg=gradopt,xxIsX0=.false.)
        if(err/=0)goto 1000 !return
      case(fd_gmeth)      ! - FD gradient
        if(adaptFDhX)then
          call adaptFDgradHx(hx,xopt,xscale,FDscale,&
                             epsF_to_epsA(epsF,fopt,fscale,Hscale),hessopt)
        endif
        call getFDCDmacro()
      case(cd_gmeth)      ! - CD gradient
        call getCDmacro()
      endselect
    case(wolfe_imeth,stwolfe_imeth)  ! gradient already available
    endselect
    call getHessFromFunc(evalFunc,dataIN,dataOUT,xopt,fopt,xscale,epsF,useHxDef=useHxDefIni,&
      hmeth=hmeth-2,hessfd=hessopt,fcalls=addFcalls,err=retcode,message=message)
    fcalls=fcalls+addFcalls
    if(retcode/=okAlg)then
      err=+10;message="f-qnewton/&"//message
      call write_exitInfo()
      goto 1000 !return
    endif
  case(bfgsInv_hmeth)         ! ** BFGS Hessian (inverse)
    selectcase(imeth)
    case(null_imeth,armijo_imeth,brentmin_imeth,trustEx_imeth,dogLeg_imeth)  ! need gradient
      selectcase(gmeth_now)
      case(user_meth)     ! - analytical gradient
        call evalFuncMacro(xx=xopt,gg=gradopt,xxIsX0=.false.)
        if(err/=0)goto 1000 !return
      case(fd_gmeth)      ! - FD gradient
        if(adaptFDhX)then
          call adaptFDgradHx(hx,xopt,xscale,FDscale,&
                             epsF_to_epsA(epsF,fopt,fscale,Hscale),hessopt)
        endif
        call getFDCDmacro()
      case(cd_gmeth)      ! - CD gradient
        call getCDmacro()
      endselect
    case(wolfe_imeth,stwolfe_imeth)  ! gradient already available
    endselect
    if(bfgsInvNR)then ! NR-based approach (efficient)
      call bfgsInv_update1(dx,xscale,activeSet,gradopt,gradold,hessopt,&
                           rescale=(iter==1.and.sclHess1it))
    else              ! Nocedal-based inverse quasi-Hessian (v inefficient)
      call bfgsInv_update2(dx,xscale,activeSet,gradopt,gradold,hessopt,hessScaled,&
                           rescale=(iter==1.and.sclHess1it))
    endif
  case(bfgsUnfac_hmeth)       ! ** BFGS Hessian (unfactored)
    selectcase(imeth)
    case(null_imeth,armijo_imeth,brentmin_imeth,trustEx_imeth,dogLeg_imeth)  ! need gradient
      selectcase(gmeth_now)
      case(user_meth)     ! - analytical gradient
        call evalFuncMacro(xx=xopt,gg=gradopt,xxIsX0=.false.)
        if(err/=0)goto 1000 !return
      case(fd_gmeth)      ! - FD gradient
        if(adaptFDhX)then
          call adaptFDgradHx(hx,xopt,xscale,FDscale,&
                             epsF_to_epsA(epsF,fopt,fscale,Hscale),hessopt)
        endif
        call getFDCDmacro()
      case(cd_gmeth)      ! - CD gradient
        call getCDmacro()
      endselect
    case(wolfe_imeth,stwolfe_imeth)  ! gradient already available
    endselect
    call bfgsUnfac_update(dx,xscale,activeSet,gradopt,gradold,hessopt,&
                          merge(epsF,sqrt(epsF),gmeth==user_meth),&
                          skipQNupdtClassic,dampedBFGS=dampedBFGS,dampFac=dampFac,&
                          rescale=(iter==1.and.sclHess1it),err=err,message=message)
  case(bfgsFac_hmeth)         ! ** BFGS Hessian (factored)
    selectcase(imeth)
    case(null_imeth,armijo_imeth,brentmin_imeth,trustEx_imeth,dogLeg_imeth)  ! need gradient
      selectcase(gmeth_now)
      case(user_meth)     ! - analytical gradient
        call evalFuncMacro(xx=xopt,gg=gradopt,xxIsX0=.false.)
        if(err/=0)goto 1000 !return
      case(fd_gmeth)      ! - FD gradient
        if(adaptFDhX)then
          call adaptFDgradHx(hx,xopt,xscale,FDscale,&
                             epsF_to_epsA(epsF,fopt,fscale,Hscale),hessopt)
        endif
        call getFDCDmacro()
      case(cd_gmeth)      ! - CD gradient
        call getCDmacro()
      endselect
    case(wolfe_imeth,stwolfe_imeth)  ! gradient already available
    endselect
    call bfgsFac_update(dx,xscale,activeSet,gradopt,gradold,hessopt,Ld,&
      tol=merge(epsF,sqrt(epsF),gmeth==user_meth),&
      controlHessCond=controlHessCond,maxHessCond=hessFacBundle%maxHessCond,&
      logDet=logdet,condest=condest,&
      facBFGS_useR2=facBFGS_useR2,facBFGS_getLLt=facBFGS_getLLt,&
      skipClassic=skipQNupdtClassic,dampedBFGS=dampedBFGS,dampFac=dampFac,&
      rescale=(iter==1.and.sclHess1it),err=err,message=message)
    if(err/=0)then  ! probably a bug
      call write_exitInfo()
      goto 1000 !return
    endif
  case(SR1unfac_hmeth)        ! ** SR1 Hessian (unfactored)
    selectcase(imeth)
    case(null_imeth,armijo_imeth,brentmin_imeth,trustEx_imeth,dogLeg_imeth)  ! need gradient
      selectcase(gmeth_now)
      case(user_meth)     ! - analytical gradient
        if(trustDidGradHess)then  ! already computed
        else
          call evalFuncMacro(xx=xopt,gg=gradopt,xxIsX0=.false.)
          if(err/=0)goto 1000 !return
        endif
      case(fd_gmeth)      ! - FD gradient
        if(trustDidGradHess)then   ! - already done (note that this precludes adaption)
        else
          if(adaptFDhX)then
            call adaptFDgradHx(hx,xopt,xscale,FDscale,&
                               epsF_to_epsA(epsF,fopt,fscale,Hscale),hessopt)
          endif
          call getFDCDmacro()
        endif
      case(cd_gmeth)      ! - CD gradient
        if(trustDidGradHess)then   ! - already done (note that this precludes adaption)
        else
          call getCDmacro()
        endif
      endselect
    case(wolfe_imeth,stwolfe_imeth)  ! gradient already available
    endselect
    if(err/=0)then
      message="f-qnewton/unknown/&"//message
      call write_exitInfo()
      goto 1000 !return
    endif
    if(.not.trustDidGradHess)then ! avoid updating twice
      call SR1unFac_update(dx,xscale,activeSet,gradopt,gradold,&
        trustBundle%SR1skipTol,hessopt,rescale=(iter==1.and.sclHess1it),&
        err=err,message=message)
    endif
  case(NCG_FR_hmeth,NCG_PR_hmeth,NCG_PPR_hmeth) ! ** Conjugate gradient methods
    selectcase(imeth)
    case(null_imeth,armijo_imeth,brentmin_imeth)  ! need gradient
      selectcase(gmeth_now)
      case(user_meth)     ! - analytical gradient
        call evalFuncMacro(xx=xopt,gg=gradopt,xxIsX0=.false.)
        if(err/=0)goto 1000 !return
      case(fd_gmeth)      ! - FD gradient
        if(adaptFDhX)then
          call adaptFDgradHx(hx,xopt,xscale,FDscale,&
                             epsF_to_epsA(epsF,fopt,fscale,Hscale),hessopt)
        endif
        call getFDCDmacro()
      case(cd_gmeth)      ! - CD gradient
        call getCDmacro()
      endselect
    case(wolfe_imeth,stwolfe_imeth)  ! gradient already available
    endselect
    gg=dot_product(gradOld,gradOld)
    selectcase(hmeth)
    case(NCG_FR_hmeth)              ! - Fletcher-Reeves algorithm
      dgg=dot_product(gradOpt,gradOpt)
    case(NCG_PR_hmeth)              ! - Polak-Ribiere (somewhat more graceful - NR-77)
      dgg=dot_product(gradOpt-gradOld,gradOpt)
    case(NCG_PPR_hmeth)             ! - PR+ (Nocedal and Wright)
      dgg=dot_product(gradOpt-gradOld,gradOpt)  !   more resistant to inexact linesearches than PR
      if(dgg<zero)dgg=zero
    endselect
    gam=dgg/gg; Ld=-gradopt+gam*Ld
  case default
    err=+10;message="f-qnewton/unknownHessMethod"
    call write_exitInfo()
    goto 1000 !return
  endselect
  selectcase(chkGrd)    ! * check gradient at every step?
  case(chkG_dxstp,chkG_hxstp,chkG_full)
    call checkGrad_macro()
    if(err>0)then
      call write_exitInfo()
      goto 1000 !return
    endif
  endselect
  selectcase(chkHess)   ! * check Hessian at every step?
  case(chkHess_full)
    call checkHess_macro()
    if(err>0)then
      call write_exitInfo()
      goto 1000 !return
    endif
  endselect
! * Check active set after each iteration.
  skipDxDfCheck=.false.
  if(boundedSearch.and.globcode/=failed_glob)then
!  if(boundedSearch)then
! - Safeguard rare case when failing to globalise in bounded optimisation creates
!   an infinite loop with skipping convergence test. Failure to globalise is final!!!
    call checkActiveSet(xopt,xscale,activeSet,fixDiagOption,hmeth,gradopt,hessopt,&
                        Ld,xLo,xHi,hitBound,nfree0,nfix,nthawn)
! - check whether bound-deletion convergence satisfied
    call checkConvergence(xopt,dx,fopt,gradOpt,&  ! this checks loose convergence
      merge(freeVar_as,loVar_as,activeSet==freeVar_as),&
      user_meth,fredExp,fredAct,xscale,fscale,&       
      gtol,tolOptSlack_bnd*stol,tolOptSlack_bnd*ftol,&  ! gtol not relaxed by tolOptSlack_bnd
!      tolOptSlack_bnd*gtol,tolOptSlack_bnd*stol,tolOptSlack_bnd*ftol,&
      hitBound,retcode)
    selectcase(retcode) ! establish whether to force constraint deletion
    case(grad_con,search_con,fred_con)
      delCon=.true.
    case default
      delCon=.false.
    endselect
    call checkReleaseActiveSet(xopt,xscale,activeSet,gradopt,hessopt,Ld,&
         tolGfree_bnd,nfree==0.or.delCon,tolGfree2_bnd,nfree,nfix,nthawn)
    if(nFree>nfree0)then
      skipDxDfCheck=.true. ! skip convergence check if variable released
      retcode=no_con
    endif
  elseif(delCon)then  ! do no skip convergence test after releasing constraint
    hitBound=.false.  ! (otherwise infinite loop can occur release/fix)
  elseif(globcode==failed_glob)then ! failed to globalize
    skipDxDfCheck=.true.
  endif
  call write_iterationInfo()
! * Check convergence criteria (gradient/search/function tolerance)
  call checkConvergence(xopt,dx,fopt,gradOpt,activeSet,gmeth_now,fredExp,fredAct,&
                        xscale,fscale,gtol,stol,ftol,hitBound.or.skipDxDfCheck,retcode)
  selectcase(gmeth) ! some extra logic when gradient is approximated
  case(fd_gmeth,cd_gmeth)
    if(useHxDef.and.hybridFDCD)then
! - if using default FD gradient stepsize with hybrid component adaption,
!   do not bother switching, since gradient would have already been estimated
!   using CD (provided tolGradFDCD ~ 0.1 or so, as it should be)
      selectcase(retcode)
      case(-grad_con,-search_con,-fred_con)
        retcode=abs(retcode)
      endselect
    elseif(.not.useHxDef.and.gradHx/=fresh_hx)then  ! insist on refreshing stepsize 
      selectcase(retcode)                           ! before termination
      case(grad_con,search_con,fred_con)
        retcode=switchCD_con
      endselect
    endif
  endselect
  selectcase(retcode)
  case(no_con)              ! no convergence yet
    gradHx=old_hx; freshHess=.false.  ! indicate that gradient stepsize unrefreshed
  case(grad_con)            ! gradient converged
    err=0;message="qnewton/ok/grad[x]~0"
    call write_exitInfo()
    goto 1000 !return
  case(search_con)          ! search converged
    err=0;message="w-qnewton/prob.ok/||dx||~0"
    call write_exitInfo()
    goto 1000 !return
  case(fred_con)            ! function converged
    err=0;message="w-qnewton/prob.ok/df~0"
    call write_exitInfo()
    goto 1000 !return
  case(srchBadGrad_con)     ! search converged but grad still large
    err=-10;message="w-qnewton/sus/||dx||~0,grad[x]>>0"
    call write_exitInfo()
    goto 1000 !return
  case(fredBadGrad_con)     ! function converged but grad still large
    err=-20;message="w-qnewton/sus/df~0,grad[x]>>0"
    call write_exitInfo()
    goto 1000 !return
  case(switchCD_con,&
       -grad_con,-search_con,-fred_con,&
       -srchBadGrad_con,-fredBadGrad_con)
! need to switch to central differences to continue progress
    retcode=switchCD_con
    call write_FDCDswitchInfo(fd_to_cd=.true.)
    gmeth_now=cd_gmeth
    hx=getHxFromRelHx(hx,xopt,xscale,FDscale)
    call getHx_macro()    ! re-optimise stepsize
    if(err/=0)goto 1000 !return
    hx=getRelHxFromHx(hx,xopt,xscale,FDscale)
    gradHx=fresh_hx
  endselect
  if(fcalls>maxFev)then
    err=-100;message="f-qnewton/maxFevalExceeded"
    call write_exitInfo()
    goto 1000 !return
  endif
  selectcase(globcode)  ! check global convergence
  case(success_glob)    ! ... keep going
  case default          ! ... failed to globalize but convergence not verified
    err=-50;message="f-qnewton/warning/lastStepFailedNoTermination/checkResults..."
    call write_exitInfo()
    goto 1000 !return
  endselect
!-----
! Check conditions for switch between forward <-> central differences
  if(gradHx/=fresh_hx)then
    if(.not.(hybridFDCD.and.useHxDef))then
      call checkFDCDswitch_macro()
      if(err/=0)goto 1000 !return
!    elseif(gmeth_now==cd_gmeth)then ! encourage switches CD->FD
!      call checkFDCDswitch_macro()
!      if(err/=0)goto 1000 !return
    endif
  endif
! check for false convergence (succesive scaled step lengths too small)
  if(scaledStepLen(dx,xopt,xscale)<=tolFalseDx)then
    nFalseDx=nFalseDx+1
  else    ! reset false convergence counter
    nFalseDx=0
  endif
  if(nFalseDx>nFalseDxMax)then
    err=30; message="f-qnewton/falseConvergence(nonCriticalPoint)/slowProgress"
    call write_exitInfo()
    goto 1000 !return
  elseif(mod(nFalseDx+1,nFalseRfrshDxMax)==0)then ! refresh stepsize?
    if(gmeth_now/=user_meth.and..not.useHxDef)then
      if(uout>0)write(uout,'(a)')"False convergence / slow progress detected"
      hx=getHxFromRelHx(hx,xopt,xscale,FDscale)
      call getHx_macro()    ! re-optimise stepsize
      if(err/=0)goto 1000 !return
      hx=getRelHxFromHx(hx,xopt,xscale,FDscale)
      gradHx=fresh_hx
    endif
    if(useQuasiHessian.and.allowQHreset)call setUnitQhess_macro()
  endif
enddo ! ... proceed to next Newton iteration
err=20; message="f-qnewton/maxIterExceeded"
call write_exitInfo()
! just before return: clean up heap memory
1000 call cleanupMem()
! End main procedure here
contains
!-----
subroutine startupMem()   ! macro to allocate heap memory
use utilities_dmsl_kit,only:BperMB
implicit none
integer(mik)::memHess2temp
if(useConjGrad)then   ! no Hessian storge needed
  memHess2temp=0      ! (big advantage of conjugate-gradient methods for large problems)
else
  memHess2temp=ndim**2*mrkBy/BperMB ! extra MB's needed for Hessian allocation on the heap
  allocate(hessScaled(ndim,ndim),stat=err)
  if(err/=0)then  ! no need to preserve err values: its always err=0 when this routine is called
    err=101
    write(message,'(a,i0,a)')"f-startupMem/allocError:hessScaled[",size(x0),"]/tryConjugateGradientMethods"
  endif
endif
if(present(memHess2))memHess2=memHess2temp
endsubroutine startupMem
!-----
subroutine cleanupMem()   ! macro to clean up the memory space and avoid memory leaks
implicit none             ! cannot rely on the compiler to do so
integer(mik)::jerr        ! local var to avoid obliterating the "err" return status
if(allocated(hessScaled))then
  deallocate(hessScaled,stat=jerr)
else
  jerr=0
endif
if(jerr/=0)then
  jerr=102;if(err==0)err=jerr      ! but if this happened on a normal return then flag a problem
  write(message,'(a,i0,a)')"f-cleanupMem/deAllocError:hessScaled[",size(x0),"]/strange(bug?)"
endif
endsubroutine cleanupMem
!-----
subroutine processUnwiseSettings()  ! macro to handle unwise parameters
use utilities_dmsl_kit,only:oneThird,twoThirds
implicit none
! locals
type(qnewtonUnwise_type)::qnewtonUnwiseLoc ! default settings held here
! Start procedure here
if(present(qnewtonUnwise))then  ! overwrite default settings
  qnewtonUnwiseLoc=qnewtonUnwise
  if(uout>0)write(uout,'(a)')"Warning:qnewtonUnwise prescribed"
endif
! Initial point analysis
gtol0fac        =qnewtonUnwiseLoc%gtol0fac
! Linesearch settings
alpha_ls        =qnewtonUnwiseLoc%alpha_ls 
beta_ls         =merge(qnewtonUnwiseLoc%beta_ls_CG,qnewtonUnwiseLoc%beta_ls,useConjGrad)
LNSstrongwolfe  =qnewtonUnwiseLoc%LNSstrongwolfe
useDirDer       =qnewtonUnwiseLoc%useDirDer
linmin_ometh    =qnewtonUnwiseLoc%linmin_ometh
linmin_tol      =qnewtonUnwiseLoc%linmin_tol
linmin_itmax    =qnewtonUnwiseLoc%linmin_itmax
! Trust region settings
acceptRatio_tr  =qnewtonUnwiseLoc%acceptRatio_tr
roDown_tr       =qnewtonUnwiseLoc%roDown_tr
radDown_tr      =qnewtonUnwiseLoc%radDown_tr
roUp_tr         =qnewtonUnwiseLoc%roUp_tr
stepOtrustUp_tr =qnewtonUnwiseLoc%stepOtrustUp_tr
radUp_tr        =qnewtonUnwiseLoc%radUp_tr
roUpNow_tr      =qnewtonUnwiseLoc%roUpNow_tr
trustOstepMax_tr=qnewtonUnwiseLoc%trustOstepMax_tr
niter_tr        =qnewtonUnwiseLoc%niter_tr
ncholMax_tr     =qnewtonUnwiseLoc%ncholMax_tr
SR1forceUpdt    =qnewtonUnwiseLoc%SR1forceUpdt
pivotCholTrust  =qnewtonUnwiseLoc%pivotCholTrust
! - do not use pivoting with factored BFGS quasi-Newton
if(imeth==dogLeg_imeth.and.hmeth==bfgsFac_hmeth)pivotCholTrust=.false.
dogNewtBias     =qnewtonUnwiseLoc%dogNewtBias
boundFrac       =qnewtonUnwiseLoc%boundFrac
! Quasi-Hessian update settings
skipQNupdtClassic=qnewtonUnwiseLoc%skipQNupdtClassic
allowQHreset    =qnewtonUnwiseLoc%allowQHreset
maxSR1update    =qnewtonUnwiseLoc%maxSR1update
facBFGS_useR2   =qnewtonUnwiseLoc%facBFGS_useR2
facBFGS_getLLt  =qnewtonUnwiseLoc%facBFGS_getLLt
facBFGS_typeH   =merge("SCL","SU ",hmeth==bfgsFac_hmeth.and..not.facBFGS_getLLt)
dampedBFGS      =qnewtonUnwiseLoc%dampedBFGS
dampFac         =qnewtonUnwiseLoc%dampFac
! Hessian scaling method
xscaleHmeth     =qnewtonUnwiseLoc%xscaleHmeth
! Function roundoff estimation
Hscale          =qnewtonUnwiseLoc%Hscale
hammPow         =qnewtonUnwiseLoc%hammPow
! Performance output
iterNfo         =qnewtonUnwiseLoc%iterNfo
! Active set bound constraints handling
tolGfree_bnd    =qnewtonUnwiseLoc%tolGfree_bnd
tolOptSlack_bnd =qnewtonUnwiseLoc%tolOptSlack_bnd
tolGfree2_bnd   =qnewtonUnwiseLoc%tolGfree2_bnd
fixDiagOption   =qnewtonUnwiseLoc%fixDiagOption
! False convergence analysis
tolFalseDx      =qnewtonUnwiseLoc%tolFalseDx
nFalseDxMax     =qnewtonUnwiseLoc%nFalseDxMax
nFalseRfrshDxMax=qnewtonUnwiseLoc%nFalseRfrshDxMax
! Gradient checking
chkGrd          =qnewtonUnwiseLoc%chkGrd   
chkGrd_gmeth    =qnewtonUnwiseLoc%chkGrd_gmeth
chkGrd_tG       =qnewtonUnwiseLoc%chkGrd_tG
chkGrd_tGdf     =qnewtonUnwiseLoc%chkGrd_tGdf
chkGrd_tF       =qnewtonUnwiseLoc%chkGrd_tF
chkGrd_h        =qnewtonUnwiseLoc%chkGrd_h
! Hessian checking
chkHess         =qnewtonUnwiseLoc%chkHess
chkHess_hmeth   =qnewtonUnwiseLoc%chkHess_hmeth
ignoreBadHess   =qnewtonUnwiseLoc%ignoreBadHess
! Finite difference gradient approximation
FDscale         =qnewtonUnwiseLoc%FDscale
useHxDef        =qnewtonUnwiseLoc%useHxDef
hybridFDCD      =qnewtonUnwiseLoc%hybridFDCD
dfdx0meth       =qnewtonUnwiseLoc%dfdx0meth
allowFDCD       =qnewtonUnwiseLoc%allowFDCD
tolFDCD         =qnewtonUnwiseLoc%tolFDCD
fracFDCD        =qnewtonUnwiseLoc%fracFDCD
tolCDFD         =qnewtonUnwiseLoc%tolCDFD
fracCDFD        =qnewtonUnwiseLoc%fracCDFD
tolGradFDCD     =qnewtonUnwiseLoc%tolGradFDCD
tolGradCDFD     =qnewtonUnwiseLoc%tolGradCDFD
tolDxFDCD       =qnewtonUnwiseLoc%tolDxFDCD
adaptFDhX       =qnewtonUnwiseLoc%adaptFDhX
adaptCDhX       =qnewtonUnwiseLoc%adaptCDhX
! Modified Hessian factorization settings
facmeth         =qnewtonUnwiseLoc%facmeth
! - no pivoting with Cholesky-Gershgorin factorization
if(facmeth==dennis_facmeth)pivotCholTrust=.false.
tau             =qnewtonUnwiseLoc%tau    
tau             =merge(epsRe**oneThird,tau,tau<zero) ! default, but F-95 precludes
tauBar          =qnewtonUnwiseLoc%tauBar             ! automatic initialisation
tauBar          =merge(epsRe**twoThirds,tauBar,tauBar<zero)
mu              =qnewtonUnwiseLoc%mu     
maxHessCond     =qnewtonUnwiseLoc%maxHessCond ! * check for +ve condition bound
maxHessCond     =merge(one/sqrt(epsRe),maxHessCond,maxHessCond<zero)
controlHessCond =qnewtonUnwiseLoc%controlHessCond
! End procedure here
endsubroutine processUnwiseSettings
!-----
subroutine makeGmethBundle() ! macro to collect gradient-evaluation bundle
implicit none
! Start procedure here
gmethBundle%gmeth_now         =>gmeth_now
gmethBundle%useHxDef          = useHxDef
gmethBundle%FDscale           = FDscale
gmethBundle%hx                =>hx
gmethBundle%hybridFDCD        = hybridFDCD
gmethBundle%tolGradFDCD       = tolGradFDCD
gmethBundle%useDirDer         = useDirDer
! End procedure here
endsubroutine makeGmethBundle
!-----
subroutine makeTrustBundle() ! macro to collect trust-evaluation bundle
use utilities_dmsl_kit,only:oneThird
implicit none
! Start procedure here
trustBundle%acceptRatio_tr    = acceptRatio_tr  
trustBundle%roDown_tr         = roDown_tr       
trustBundle%radDown_tr        = radDown_tr      
trustBundle%roUp_tr           = roUp_tr         
trustBundle%stepOtrustUp_tr   = stepOtrustUp_tr 
trustBundle%radUp_tr          = radUp_tr        
trustBundle%roUpNow_tr        = roUpNow_tr
trustBundle%trustOstepMax_tr  = trustOstepMax_tr
trustBundle%niter_tr          = niter_tr        
trustBundle%ncholMax_tr       = ncholMax_tr     
trustBundle%trustMax          = stepmax
trustBundle%trustMin          = stol
!trustBundle%SR1skipTol        = merge(sqrt(epsF),sqrt(sqrt(epsF)),gmeth==user_meth)
trustBundle%SR1skipTol        = merge(sqrt(epsF),epsF**oneThird,gmeth==user_meth)
trustBundle%SR1forceUpdt      = SR1forceUpdt
trustBundle%pivotCholTrust    = pivotCholTrust
trustBundle%dogNewtBias       = dogNewtBias
trustBundle%boundFrac         = boundFrac
! End procedure here
endsubroutine makeTrustBundle
!-----
subroutine makeObjFuncBundle() ! macro to collect function-evaluation bundle
implicit none
! Start procedure here
objFuncBundle%epsF   = epsF
objFuncBundle%Hscale = Hscale
! End procedure here
endsubroutine makeObjFuncBundle
!-----
subroutine makeHessFacBundle() ! macro to collect Hessian factorization bundle
implicit none
! Start procedure here
hessFacBundle%facmeth = facmeth
hessFacBundle%tau     = tau
hessFacBundle%tauBar  = tauBar
hessFacBundle%mu      = mu
hessFacBundle%maxHessCond = maxHessCond
! End procedure here
endsubroutine makeHessFacBundle
!-----
elemental function useConjGrad_inq(hmeth)     ! returns true if conjugate gradient used
implicit none
! dummies
integer(mik),intent(in)::hmeth
logical(mlk)::useConjGrad_inq
! Start procedure here
selectcase(hmeth)
case(NCG_FR_hmeth,NCG_PR_hmeth,NCG_PPR_hmeth)
  useConjGrad_inq=.true.
case default
  useConjGrad_inq=.false.
endselect
! End procedure here
endfunction useConjGrad_inq
!-----
elemental function useQuasiHessian_inq(hmeth) ! returns true if quasi-Newton method used
implicit none
! dummies
integer(mik),intent(in)::hmeth
logical(mlk)::useQuasiHessian_inq
! Start procedure here
selectcase(hmeth)
case(bfgsInv_hmeth,bfgsUnfac_hmeth,bfgsFac_hmeth,SR1unFac_hmeth)
  useQuasiHessian_inq=.true.
case default
  useQuasiHessian_inq=.false.
endselect
! End procedure here
endfunction useQuasiHessian_inq
!-----
elemental function useTrust_inq(imeth) ! returns true if trust region method used
implicit none
! dummies
integer(mik),intent(in)::imeth
logical(mlk)::useTrust_inq
! Start procedure here
selectcase(imeth)
case(trustEx_imeth,dogLeg_imeth)
  useTrust_inq=.true.
case default
  useTrust_inq=.false.
endselect
! End procedure here
endfunction useTrust_inq
!-----
subroutine getFDCDmacro() ! macro to evaluate FDCD gradient
use utilities_dmsl_kit,only:getFDCDgrad
implicit none
! Start procedure here
call getFDCDgrad(evalFunc,dataIN,dataOUT,xopt,fopt,xscale,fscale,epsF,&
  getHxFromRelHx(hx,xopt,xscale,FDscale),useHxDef,&
  merge(useFDCDhybrid,fd_gmeth,hybridFDCD),tolGradFDCD,&
  gradopt,addFcalls,err,message)
fcalls=fcalls+addFcalls
! End procedure here
endsubroutine getFDCDmacro
!-----
subroutine getCDmacro() ! macro to evaluate CD gradient
use utilities_dmsl_kit,only:getCDgrad
implicit none
! Start procedure here
call getCDgrad(evalFunc,dataIN,dataOUT,xopt,fopt,xscale,epsF,&
  getHxFromRelHx(hx,xopt,xscale,FDscale),useHxDef,&
  gradopt,addFcalls,err,message)
fcalls=fcalls+addFcalls
! End procedure here
endsubroutine getCDmacro
!-----
subroutine evalFuncMacro(xx,ff,gg,hh,xxIsX0)  ! macro to evaluate function
implicit none
! dummies
real(mrk),intent(in)::xx(:)
real(mrk),intent(out),optional::ff,gg(:),hh(:,:)
logical(mlk),intent(in)::xxIsX0
! locals
logical(mlk)::ok
! Start procedure here
call evalFunc(dataIN,dataOUT,xx,ok,ff,gg,hh,err=err,message=message)
if(err/=0)then
  err=100; write(message,'(a,i0,a)')"f-qnewton/userErr[iter=",iter,"]/&"//trim(message)
elseif(.not.ok)then
  if(xxIsX0)then
    err=10;message="f-qnewton/x0unfeas"
  else
    err=20;write(message,'(a,i0,a)')"f-qnewton/bug/x(accepted)Unfeas[iter=",iter,"]"
  endif
  call write_exitInfo()
else
  err=0
  if(present(ff))fcalls=fcalls+1
  if(present(gg))gcalls=gcalls+1
  if(present(hh))hcalls=hcalls+1
endif
! End procedure here
endsubroutine evalFuncMacro
!-----
subroutine checkStepmax_macro() ! macro to set stepmax
implicit none
! Start procedure here
stepmax=stpmax
if(stepmax<=zero)then ! default value for stepmax (DS96,IMSL)
  stepmax=stmax*max(norm2(x0/xscale),norm2(one/xscale))
endif
stepmaxL=stepmax
! End procedure here
endsubroutine checkStepmax_macro
!-----
subroutine checkFDCDswitch_macro()  ! macro to switch FD<->CD
implicit none
! Start procedure here
selectcase(gmeth_now)
case(fd_gmeth)
  if(allowFDCD.and.fracFDCD<=one)then
! Designed by DK: compare FD gradient dfdx_fd with curvature d2f/dx2.
! Since as df/dx(true)->0, dfdx_fd->0.5*h*d2f/dx2
! the quasi-Hessian diagonal gives possibly useful order-of-magnitude
! estimate of d2f/dx2 and hence can be used to construct a switch condition.
    selectcase(hmeth)  ! unfactored quasi-Hessian approximations
    case(user_meth,fdg_hmeth,cdg_hmeth,fdf_hmeth,cdf_hmeth,&
         bfgsUnfac_hmeth,SR1unFac_hmeth)
      d2fdx2=abs(getdiag(hessopt))
    case(bfgsFac_hmeth)
      d2fdx2=abs(getdiag(hessopt))
    case(bfgsInv_hmeth)
      d2fdx2=one/abs(getdiag(hessopt))  ! this is wrong but probably still works...
    endselect ! ...otherwise use recursion to get d2f/dx2 from inverse Hessian
! fraction of variables where gradient is "small" relative to curvature.
    gOh_fdcd=get_gOh_fdcd(gradopt,d2fdx2,&
                          getHxFromRelHx(hx,xopt,xscale,FDscale),activeset,tolFDCD)
    if(gOh_fdcd>=fracFDCD)gmeth_now=cd_gmeth ! switch to central diffs
    gOh_fdcd=-1._mrk  ! indicate on potential later switches that value not fresh
  endif
  if(scaledGrad(gradopt,xopt,fopt,xscale,fscale,activeSet)<tolGradFDCD)then
    gmeth_now=cd_gmeth  ! gradient small enough to switch to CD
  elseif(.not.hitbound.and.&  ! poor progress possibly due to inaccurate derivatives
      scaledStepLen(dx,xopt,xscale)<tolDxFDCD)then
    gmeth_now=cd_gmeth
  endif
  selectcase(gmeth_now)
  case(cd_gmeth)    ! switch to central differences
    call write_FDCDswitchInfo(fd_to_cd=.true.)
    hx=getHxFromRelHx(hx,xopt,xscale,FDscale)
    call getHx_macro()    ! re-optimise stepsize
    if(err/=0)return
    hx=getRelHxFromHx(hx,xopt,xscale,FDscale)
    gradHx=fresh_hx
  endselect
case(cd_gmeth)
  selectcase(gmeth)
  case(fd_gmeth) ! "economy" quasi-Newton using CD gradient
    if(allowFDCD.and.fracCDFD<=one)then
      selectcase(hmeth) ! unfactored quasi-Hessian approximations
      case(user_meth,fdg_hmeth,cdg_hmeth,fdf_hmeth,cdf_hmeth,&
           bfgsUnfac_hmeth,SR1unFac_hmeth)
        d2fdx2=abs(getdiag(hessopt))
      case(bfgsFac_hmeth)
        d2fdx2=abs(getdiag(hessopt))
      case(bfgsInv_hmeth)
        d2fdx2=one/abs(getdiag(hessopt))  ! this maybe a bit optimistique (see above)
      endselect
      gOh_fdcd=get_gOh_fdcd(gradopt,d2fdx2,&
                            getHxFromRelHx(hx,xopt,xscale,FDscale),activeset,tolCDFD)
      if(gOh_fdcd>=one-fracCDFD)gmeth_now=fd_gmeth  ! switch back to forward differences
      gOh_fdcd=-1._mrk  ! indicate on potential later switches that value not fresh
    endif
    if(scaledGrad(gradopt,xopt,fopt,xscale,fscale,activeSet)>tolGradCDFD)then
      gmeth_now=fd_gmeth
    endif
    selectcase(gmeth_now)
    case(fd_gmeth) ! switch from FD -> CD but do not re-optimize stepsize
      call write_FDCDswitchInfo(fd_to_cd=.false.)
    endselect
  endselect
endselect
! End procedure here
endsubroutine checkFDCDswitch_macro
!-----
subroutine setUnitQhess_macro() ! macro to reset quasi-Hessian to unit matrix
implicit none
! Start procedure here
selectcase(hmeth) ! reset quasi-Hessian to identity
case(bfgsInv_hmeth)
  call initQHess_inv(fopt,fscale,xscale,unt_himeth,hessopt)
case(bfgsUnfac_hmeth)
  call initQHess_unfac(fopt,fscale,xscale,unt_himeth,hessopt)
case(bfgsFac_hmeth)
  call initQHess_fac(fopt,fscale,xscale,unt_himeth,hessopt,Ld,facBFGS_getLLt)
case(SR1unFac_hmeth)
  call initQHess_unfac(fopt,fscale,xscale,unt_himeth,hessopt)
endselect
! End procedure here
endsubroutine setUnitQhess_macro
!-----
subroutine getfredExp()      ! macro to evaluate predicted function reduction
use utilities_dmsl_kit,only:quadDf
implicit none
! Start procedure here
selectcase(hmeth)
case(NCG_FR_hmeth,NCG_PR_hmeth,NCG_PPR_hmeth) ! * Conjugate gradient methods
  fredExp=-dot_product(gradOpt,dx)            !   ignore Hessian portion
case(bfgsInv_hmeth)                           ! * Inverse BFGS Hessian
  fredExp=-dot_product(gradOpt,dx)            !   ignore Hessian portion
case default  ! * All other Hessians can be used in full quadratic form
  fredExp=-quadDf(dx=dx,dfdx=gradOpt,d2fdx2=hessOpt,typeH=facBFGS_typeH)
endselect
! End procedure here
endsubroutine getfredExp
!-----
subroutine writeSettings()  ! macro to write algorithm settings
implicit none
! locals
character(200)::infoString !,infoStringB
character(*),parameter::fmtIS='(3x,a,a)',fmtSN='(3x,a,es7.1)'
! Start procedure here
if(uout<=0)return
write(uout,'(a)')"*********************************************************************"
write(uout,'(a)')"------- ALGORITHMIC SETTINGS: DMSL NEWTON OPTIMIZATION MODULE -------"
write(uout,'(a)')"List of major settings and some (not all) 'esoteric' settings..."
! * Some problems specs
write(uout,'(a)')"0. PRIMARY PROBLEM CHARACTERISTIX"
write(uout,'(3x,a,i0)')"Number of variables (ndim) is ",ndim
write(uout,'(a)')"1. PRIMARY ALGORITHM SELECTION"
! * Globalization method
selectcase(imeth)
case(null_imeth)
  infoString="Null: pure Newton iterations (debugging only)"
case(armijo_imeth)
  infoString="Armijo backtracking linesearch"
case(wolfe_imeth)
  infoString="Wolfe linesearch"
case(stwolfe_imeth)
  infoString="Strong Wolfe linesearch"
case(brentmin_imeth)
  infoString="Brent line minimisation"
case(trustEx_imeth)
  infoString="Hookstep (near-exact) trust region"
case(dogLeg_imeth)
  infoString="Dogleg trust region"
case default
  infoString=unknownMethodChar
endselect
write(uout,fmtIS)"Globalization method:  ",trim(infoString)
! * Gradient method
selectcase(gmeth)
case(user_meth)
  infoString="User-supplied"
case(fd_gmeth)
  infoString="Forward difference approximation ("//&
              trim(merge("default ","adaptive",useHxDef))//" stepsize)"
case(cd_gmeth)
  infoString="Central difference approximation ("//&
              trim(merge("default ","adaptive",useHxDef))//" stepsize)"
case default
  infoString=unknownMethodChar
endselect
write(uout,fmtIS)"Gradient method:       ",trim(infoString)
! * Hessian method
selectcase(hmeth)
case(user_meth)
  infoString="User-supplied"
case(fdg_hmeth)
  infoString="Forward difference from gradient"
case(cdg_hmeth)
  infoString="Central difference from gradient"
case(fdf_hmeth)
  infoString="Forward difference from function"
case(cdf_hmeth)
  infoString="Central difference from function"
case(bfgsInv_hmeth)
  infoString="BFGS Quasi-Newton: inverse Hessian"
case(bfgsUnfac_hmeth)
  infoString="BFGS Quasi-Newton: unfactored Hessian"
case(bfgsFac_hmeth)
  infoString="BFGS Quasi-Newton: factored Hessian"
case(SR1unFac_hmeth)
  infoString="SR1 Quasi-Newton: unfactored Hessian"
case(NCG_FR_hmeth)
  infoString="Conjugate-gradient method, Fletcher-Reeves"
case(NCG_PR_hmeth)
  infoString="Conjugate-gradient method, Polak-Ribiere"
case(NCG_PPR_hmeth)
  infoString="Conjugate-gradient method, Positive Polak-Ribiere"
case default
  infoString=unknownMethodChar
endselect
write(uout,fmtIS)"Hessian method:        ",trim(infoString)
! * Initial Hessian if appropriate
if(useQuasiHessian)then
  selectcase(himeth)
  case(user_meth)
    infoString="User-supplied"
  case(unt_himeth)
    infoString="Unit matrix"
  case(untcnd1_himeth)
    infoString="Conditioned unit matrix"
  case(scld_himeth)
    infoString="Scaled unit matrix"
  case(scldcnd1_himeth)
    infoString="Conditioned scaled unit matrix"
  case(d2fdx2_himeth)
    infoString="Approximate Hessian diagonal (can be expensive)"
  case(hessX0_himeth)
    infoString="Full Hessian (very expensive)"
  case default
    infoString=unknownMethodChar
  endselect
else
    infoString="Non-quasi-Newton method: 'himeth' ignored"
endif
write(uout,fmtIS)"Initial quasi-Hessian: ",trim(infoString)
! - Quasi-Hessian update
if(useQuasiHessian)then
  selectcase(hmeth)
  case(bfgsInv_hmeth,bfgsUnfac_hmeth,bfgsFac_hmeth)
    if(skipQNupdtClassic)then
      infoString="Classic BFGS update skipping"
    elseif(dampedBFGS)then
      infoString="Damped BFGS update skipping"
    else
      infoString="DK-modified BFGS update (maybe not positive-definite)"
    endif
    write(uout,fmtIS)"BFGS updating:     ",trim(infoString)
  endselect
  if(allowQHreset)then
    infoString="Reset Quasi-Hessian to unit matrix when failing"
  else
    infoString="Do not reset Quasi-Hessian to unit matrix when failing"
  endif
  write(uout,fmtIS)"Hessian resetting: ",trim(infoString)
endif
! * Convergence tolerance information
write(uout,'(a)')"2. TERMINATION CRITERIA"
write(uout,fmtSN)"Scaled gradient tolerance (gtol): ", gtol
write(uout,fmtSN)"Scaled step tolerance (stol):     ",     stol
write(uout,fmtSN)"Scaled function tolerance (ftol): ", ftol
write(uout,fmtSN)"False convergence tolerance (tolFalseDx): ",tolFalseDx
! * Active set (bound) information
write(uout,'(a)')"3. ACTIVE-SET INFORMATION"
if(boundedSearch)then
  write(uout,fmtIS)"Bounds supplied by user: bound-constrained minimization",""
  write(uout,fmtSN)"Gradient tolerance for fast release     (>1.0=ignore):      ",tolGfree_bnd
  write(uout,fmtSN)"Slack factor for variable release       (<1.0=ignore):      ",tolOptSlack_bnd
  write(uout,fmtSN)"Gradient tolerance for standard release (>1.0=one-at-time): ",tolGfree2_bnd
else
  write(uout,fmtIS)"No bounds supplied by user: unconstrained minimization",""
endif
! * Additional information
if(useConjGrad)then  ! method
else
  write(uout,'(a)')"4. SECONDARY ALGORITHM SELECTION"
! - Hessian inversion method
  selectcase(facmeth)
  case(schnab_facmeth)
    infoString="Revised modified Cholesky of Schnabel and Eskew"
  case(dennis_facmeth)
    infoString="Robust Cholesky-Gershgorin of Dennis and Schnabel (Gill et al)"
  case default
    infoString=unknownMethodChar
  endselect
  write(uout,fmtIS)"Hessian Cholesky method: ",trim(infoString)
! - Cholesky pivoting information
  selectcase(imeth)
  case(armijo_imeth,wolfe_imeth,stwolfe_imeth,brentmin_imeth)
    selectcase(facmeth)
    case(schnab_facmeth)
      infoString="Pivoting enabled (recommended)"
    case(dennis_facmeth)
      infoString="Pivotion disabled (not recommended)"
    endselect
  case(trustEx_imeth)
    if(pivotCholTrust)then
      infoString="Pivoting enabled (probably un-necessarily)"
    else
      infoString="Pivotion disabled"
    endif
  case(dogLeg_imeth)
    if(pivotCholTrust)then
      infoString="Pivoting enabled"
    else
      infoString="Pivotion disabled (not recommended)"
    endif
  endselect
  write(uout,fmtIS)"Cholesky pivoting: ",trim(infoString)
! - Hessian ellipticity scaling
  selectcase(xscaleHmeth)
  case(xscaleH_sphere)
    infoString="Non-scaled Hessian"
  case(xscaleH_user)
    infoString="Scaled Hessian (user scale: xscale)"
  case(xscaleH_hdiag)
    infoString="Scaled Hessian (based on Hessian diagonal)"
  case default
    infoString=unknownMethodChar
  endselect
  write(uout,fmtIS)"Hessian 'ellipticity' scaling: ",trim(infoString)
endif
! - Finite difference gradient
if(gmeth/=user_meth)then
  write(uout,'(a)')"5. FINITE DIFFERENCE GRADIENT"
  write(uout,fmtIS)"Hybrid FD/CD gradient: ",merge("enabled ","disabled",hybridFDCD)
  selectcase(imeth)
  case(wolfe_imeth,stwolfe_imeth)
    if(useDirDer)then
      infoString="Fast (1 func eval) method"
    else
      infoString="Slow (n func eval) method"
    endif
    write(uout,fmtIS)"Directional derivative: ",trim(infoString)
  endselect
  write(uout,fmtIS)"Enhanced FD<->CD switches: ",merge("enabled ","disabled",allowFDCD)
endif
! - 
write(uout,'(a)')"----- END ALGORITHMIC SETTINGS: DMSL NEWTON OPTIMIZATION MODULE -----"
write(uout,'(a)')"*********************************************************************"
! End main procedure here
endsubroutine writeSettings
!-----
subroutine write_iterationInfo(exitInfo,skipDetailedExitInfo)  ! macro to write iteration info
use utilities_dmsl_kit,only:quickif
implicit none
! dummies
logical(mlk),intent(in),optional::exitInfo,skipDetailedExitInfo
! locals
character(200)::infoString
logical(mlk),parameter::exitInfoDef=.false.,skipDetailedExitInfoDef=.false.
character(*),parameter::&
  fmtM ="(a,i7,             &
          &2x,a,es22.14e3,  &
          &4(2x,a,es12.4e3),&
          &2x,a,es11.4,     &
          &2x,a,i2,         &
          &2x,a,i2,         &
          &2x,a,i3,a,i2,a,f4.1,a,&
          &3(2x,a,es11.3e3))",  &  ! main format
  fmtR ="(a,es22.14e3)",        &
  fmtI ="(a,i0)",               &
  fmtC ="(a)"
! Start procedure here
if(uout<=0)return
if(quickif(exitInfo,exitInfoDef))then ! exit information
  write(uout,fmtC)   "---------------------------"
  selectcase(err)
  case(0)   ! succesful exit
    write(uout,fmtC) "ALGORITHM EXIT SUCCESSFUL"
  case(:-1) ! warning
    write(uout,fmtC) "ALGORITHM EXIT WITH WARNING (SOLUTION MAY STILL BE ACCURATE)"
  case(1:)  ! error
    write(uout,fmtC) "ALGORITHM EXIT WITH ERROR"
  endselect
  write(uout,fmtI) "Error code: ",err
  write(uout,fmtC) "Message: "//trim(message)
  if(quickif(skipDetailedExitInfo,skipDetailedExitInfoDef))then
    write(uout,fmtC) "---------------------------"
    write(uout,fmtC) "Algorithm did not fully initialise: no further exit details available"
  else
    write(uout,fmtC) "---------------------------"
    write(uout,fmtC) "I.   Termination information at the minimum"
    write(uout,fmtR) "Function value at minimum:  fopt=    ",fopt
    write(uout,fmtR) "Scaled gradient at minimum: grad[f]= ",scaledGrad(gradopt,xopt,fopt,xscale,fscale,activeSet,.true.)
    selectcase(gmeth_now)
    case(user_meth)
      infoString="User-supplied"
    case(fd_gmeth)
      infoString="Forward difference approximation ("//&
                  trim(merge("default ","adaptive",useHxDef))//" stepsize)"
    case(cd_gmeth)
      infoString="Central difference approximation ("//&
                  trim(merge("default ","adaptive",useHxDef))//" stepsize)"
    case default
      infoString=unknownMethodChar
    endselect
    write(uout,fmtC) "Gradient method at termination: "//trim(infoString)
    write(uout,fmtR) "Final scaled step: dx= ",scaledStepLen(dx,xopt,xscale)
    write(uout,fmtR) "Final observed scaled function reduction: dfObs= ",scaledFred(fredAct,fopt,fscale)
    write(uout,fmtR) "Final expected scaled function reduction: dfExp= ",scaledFred(fredExp,fopt,fscale)
    write(uout,fmtI) "Total number of function calls: fcalls= ",fcalls
    write(uout,fmtC) "---------------------------"
    write(uout,fmtC) "II.  Hessian information at the minimum (estimated)"
    write(uout,fmtR) "Condition number of Hessian:      condH=   ",condEst
    write(uout,fmtR) "Log(e)-determinant of Hessian:    logDetH= ",logDet
    write(uout,fmtR) "Bound on the magnitude of the most negative Hessian eigenvalue: |Einf|= ",Einf
    write(uout,fmtC) "---------------------------"
    write(uout,fmtC) "III. Constraint information at the minimum"
    if(boundedSearch)then
      write(uout,fmtC) "Hit bound on last step: "//merge("yes","no ",hitBound)
      write(uout,fmtI) "Number of free variables:  nfree= ",nfree
      write(uout,fmtI) "Number of fixed variables: nfix=  ",nfix
      write(uout,fmtI) "Number of thawed variables (Lagrange<0): nthawn= ",nthawn
    else
      write(uout,fmtC) "Unconstrained optimisation (No constraints were specified)"
    endif
    write(uout,fmtC) "---------------------------"
    write(uout,fmtC) "IV.  Function information at the minimum"
    write(uout,fmtR) "Function minimum: fopt= ",fopt
    call write_iterationInfo_aux1(uout=uout,vecLabel="xOpt:",vecLabel_i="x",&
      vecA1=xOpt,   prec=8,boundedSearch=boundedSearch,activeSet=activeSet)
    call write_iterationInfo_aux1(uout=uout,vecLabel="gOpt:",vecLabel_i="g",&
      vecA1=gradOpt,prec=8,boundedSearch=boundedSearch,activeSet=activeSet)
    call write_iterationInfo_aux1(uout=uout,vecLabel="hOpt:",vecLabel_i="h",&
      vecA2=hessOpt,prec=8,boundedSearch=boundedSearch,activeSet=activeSet)
  endif
  write(uout,fmtC) "---------------------------"
  write(uout,fmtC) "Thnx 4 uzing d'z q-nutn, ketch u L8er ..."
  write(uout,fmtC) "---------------------------"
else
  selectcase(iterNfo)
  case(iterNfo_no)                ! no iteration info
  case(iterNfo_summ,iterNfo_var)  ! iteration summary w/wo variables
! standard summary
    write(uout,fmtM,advance="no")&
      "iter"//"["//merge("U","C",nfree==ndim)//"]=",iter,& ! iteration number
! (U=unconstrained iteration (inside domain) vs C=constrained iteration (sliding/hitting bounds)
      "fx=",        fopt,                             & ! function value
      "grad[f]=",   scaledGrad(gradopt,xopt,fopt,xscale,fscale,activeSet,.true.),& ! scaled gradient
      "dx=",        scaledStepLen(dx,xopt,xscale),    & ! scaled step
      "dfObs=",     scaledFred(fredAct,fopt,fscale),  & ! observed scaled function reduction
      "dfExp=",     scaledFred(fredExp,fopt,fscale),  & ! expected scaled function reduction
      "fcalls=",    real(fcalls,4),                   & ! function calls
      "glob_code=", globcode,                         & ! globalisation return code
      "gmeth_now=", gmeth_now,                        & ! gradient method
      "nfac=",nfacstats(1),"(",nfacstats(2)," x",     & ! number of Hessian factorizations
              real(nfacstats(1))/real(max(nfacstats(2),1)),")",&
      "condH=",     condEst,                          & ! condition number of Hessian (estimated)
      "logdetH=",   logDet,                           & ! log-det[Hessian], estimated
      "EinfH=",     Einf                                ! |Einf| (magnitude of most negative Hessian eigenvalue), estimated
! trust region info
    selectcase(imeth)
    case(trustEx_imeth,dogLeg_imeth)
      write(uout,'(2x,a,es11.3e3)',advance="no")      &
      "trustRad=",  trustRad                            ! trust radius
    endselect
! active set info
    if(boundedSearch)then
      write(uout,'(2x,a,3(2x,a,i6))',advance="no")    &
      "hitBound= "//merge("yes","no ",hitBound),      & ! indicates whether step hit bound
      "nfree=",     nfree,                            & ! free variables
      "nfix=",      nfix,                             & ! fixed variables
      "nthaw=",     nthawn                              ! thawed variables (Lagrange<0)
    endif
! prepare for possible dump of variable's info
    write(uout,'(a)',advance=merge("yes","no ",iterNfo==iterNfo_summ)) " "
    selectcase(iterNfo)
    case(iterNfo_var)  ! append current optimum to line
      call write_iterationInfo_aux1(uout=uout,vecLabel="  xOpt:  ",vecLabel_i="x",&
        vecA1=xopt,prec=4,boundedSearch=boundedSearch,activeSet=activeSet)
    endselect
  endselect
endif
! End procedure here
endsubroutine write_iterationInfo
!-----
subroutine write_iterationInfo_aux1(uout,vecLabel,vecLabel_i,vecA1,vecA2,prec,boundedSearch,activeSet)
! Purpose: writes a standard inline vector (eg, optimum, gradient, Hessian diagonal)
use utilities_dmsl_kit,only:quickif
implicit none
! dummies
integer(mik),intent(in)::uout,activeSet(:)
character(*),intent(in),optional::vecLabel
character(*),intent(in)::vecLabel_i
real(mrk),intent(in),optional::vecA1(:)
real(mrk),intent(in),optional::vecA2(:,:)
integer(mik),intent(in)::prec ! output precision: 4 (single) -> 8 (double)
logical(mlk),intent(in)::boundedSearch
! locals
integer(mik)::i
integer(mik),parameter::var_len=100
character(var_len)::fmtSS,fmtSP,fmtU,fmtV,vecLabel0
character(*),parameter::&
  fmtSS4="(a,i0,a,ss,i2,s,a,es14.6e3, 2x)", & ! single-precision format for activeSet with no "+"
  fmtSS8="(a,i0,a,ss,i2,s,a,es22.14e3,2x)", & ! double-precision format for activeSet with no "+"
  fmtSP4="(a,i0,a,sp,i2,s,a,es14.6e3, 2x)", & ! single-precision format for activeSet with "+"
  fmtSP8="(a,i0,a,sp,i2,s,a,es22.14e3,2x)", & ! double-precision format for activeSet with "+"
  fmtU4= "(a,i0,          a,es14.6e3, 2x)", & ! single-precision format for variable (no active set)
  fmtU8= "(a,i0,          a,es22.14e3,2x)", & ! double-precision format for variable (no active set)
  fmtLA='a,2x',fmtLA1='('//fmtLA//')',fmtLA2='('//fmtLA//',', &  ! format for vecLabel_i 
  fmtLB='a',   fmtLB1='('//fmtLB//')',fmtLB2='('//fmtLB//','     ! format for vecLabel_i 
integer(mik)::vecType
integer(mik),parameter::vectType_arr1=1,vectType_arr2=2
logical(mlk)::pres1,pres2
! Start procedure here
selectcase(prec)
case(4)
  fmtSS=fmtSS4;fmtSP=fmtSP4;fmtU=fmtU4
case(8)
  fmtSS=fmtSS8;fmtSP=fmtSP8;fmtU=fmtU8
case(16)
!case default
  write(uout,'(a)')"BUGERRO:f-write_iterationInfo_aux1/badIN:prec/={4,8}"
endselect
vecLabel0=quickif(vecLabel," ",var_len)
pres1=present(vecA1); pres2=present(vecA2)
if(pres1.and.pres2)then
  write(uout,'(a)')"BUGERRO:f-write_iterationInfo_aux1/badIN:pres1.and.pres2"
elseif(pres1)then
  vecType=vectType_arr1
elseif(pres2)then
  vecType=vectType_arr2
else
  return
endif
if(boundedSearch)then ! include active set info using pretty notation
  if(len_trim(vecLabel0)>0)then ! acrobatics to get right spacing before 1st entry
    write(uout,fmtLA1,advance="no")trim(vecLabel0)
  else
    write(uout,fmtLB1,advance="no")trim(vecLabel0)
  endif
  do i=1,ndim
    fmtV=merge(fmtSS,fmtSP,activeSet(i)==0)
    selectcase(vecType)
    case(vectType_arr1)
      write(uout,fmtV,advance="no")vecLabel_i,i,"(",activeSet(i),")=",vecA1(i)
    case(vectType_arr2)
      write(uout,fmtV,advance="no")vecLabel_i,i,"(",activeSet(i),")=",vecA2(i,i)
    endselect
  enddo
  write(uout,'(a)',advance="yes")" " ! terminate line
else
  if(len_trim(vecLabel0)>0)then ! acrobatics to get right spacing before 1st entry
    write(fmtV,'(a,i0,a,a)')fmtLA2,ndim,trim(fmtU),')'
  else
    write(fmtV,'(a,i0,a,a)')fmtLB2,ndim,trim(fmtU),')'
  endif
  selectcase(vecType)
  case(vectType_arr1)
    write(uout,fmtV)trim(vecLabel0),(vecLabel_i,i,"=",vecA1(i),i=1,ndim)
  case(vectType_arr2)
    write(uout,fmtV)trim(vecLabel0),(vecLabel_i,i,"=",vecA2(i,i),i=1,ndim)
  endselect
endif
! End procedure here
endsubroutine write_iterationInfo_aux1
!-----
subroutine write_exitInfo(skipDetailedExitInfo)  ! macro to write exit info
implicit none
! dummies
logical(mlk),intent(in),optional::skipDetailedExitInfo
! locals
logical(mlk),parameter::exitInfo=.true.
! Start procedure here
call write_iterationInfo(exitInfo,skipDetailedExitInfo)
! End procedure here
endsubroutine write_exitInfo
!-----
subroutine getHx_macro()  ! macro to compact hx-estimation code
implicit none
! Start procedure here
! optimise finite difference interval and compute derivatives at initial point
if(useHxDef)then
  selectcase(gmeth_now)
  case(fd_gmeth)
    call getFDCDgrad(evalFunc,dataIN,dataOUT,xopt,fopt,xscale,fscale,epsF,&
                     hx,useHxDef,&
                     merge(useFDCDhybrid,fd_gmeth,hybridFDCD),tolGradFDCD,&
                     gradopt,addFcalls,err,message)
  case(cd_gmeth)
    call getCDgrad(evalFunc,dataIN,dataOUT,xopt,fopt,xscale,epsF,&
                   hx,useHxDef,gradopt,addFcalls,err,message)
  endselect
  fcalls=fcalls+addFcalls
  selectcase(himeth)
  case(d2fdx2_himeth) ! independent evaluation of d2f/dx2
    call getHessDiagFromFunc(evalFunc,dataIN,dataOUT,xopt,fopt,xscale,epsF,useHxDef=useHxDefIni,&
      hmeth=1,hessDiag=d2fdx2,fcalls=addFcalls,err=err,message=message)
    fcalls=fcalls+addFcalls
  endselect
  if(uout>0)write(uout,'(a)',advance="no")"automatic stepsize in getHx_macro..."
  if(err/=0)then
    err=+20;message="f-getHx_macro/&"//message
    if(uout>0)then
      write(uout,'(a)')"NOT.OK"
      write(uout,'(a)')"Message: "//trim(message)
    endif
  else
    if(uout>0)write(uout,'(a)')"OK"
  endif
else
  if(uout>0)write(uout,'(a)',advance="no")"Estimating FD stepsize... "
  call getFDgradHx(evalFunc,dataIN,dataOUT,xopt,fopt,xscale,fscale,activeSet,epsF,Hscale,&
                   merge(dfdx0meth,gradFD_sw2,gmeth_now==fd_gmeth),&
                   hx,xLo,xHi,gradopt,d2fdx2,addFcalls,err,errj,messagej)
  fcalls=fcalls+addFcalls
  if(err/=0)then  ! some kind of error
    message="f-getHx_macro/BUG?/&"//messagej(1)
    call write_exitInfo()
    return
  elseif(any(errj/=0))then
    if(uout>0)write(uout,'(a)')"NOT.OK"
    err=merge(+30,0,any(errj>0))
    message="f-getHx_macro/&"//messagej(iFirstTrueLoc(errj/=0))
    call write_getFDgradHx_error()
    where(errj/=0)d2fdx2=one  ! go for safety: if df/dx poor then d2f/dx2 also
  else                        ! likely to be poor
    if(uout>0)write(uout,'(a)')"OK"
  endif
endif
! End procedure here
endsubroutine getHx_macro
!-----
subroutine write_getFDgradHx_error()  ! macro to write error from getFDgradHx
implicit none
integer(mik)::j
! Start procedure here
if(uout>0)then ! write warning to output file
  write(uout,'(a)')"WARNING IN GETFDGRADHX"
  do j=1,size(errj) ! report troublesome components
    if(errj(j)/=0)then
      write(uout,'(a,i6,a,i4,a,2(a,es19.12))')&
          "Var",j,";  err=",errj(j),"; message="//trim(messagej(j)),&
          "; df/dx=",gradopt(j),"; d2f/dx2=",d2fdx2(j)
    endif
  enddo
endif
! End procedure here
endsubroutine write_getFDgradHx_error
!-----
subroutine write_FDCDswitchInfo(fd_to_cd) ! macro to write switch info
implicit none
logical(mlk),intent(in)::fd_to_cd
character(21)::string
! Start procedure here
if(uout>0)then
  string=merge("SWITCH:FD->CD (ratio=","SWITCH:CD->FD (ratio=",fd_to_cd)
  write(uout,'(a,f5.2,a,i7,a)')string,gOh_fdcd,") (iter=",iter,")"
endif
! End procedure here
endsubroutine write_FDCDswitchInfo
!-----
subroutine checkGrad_macro()  ! macro to check gradient accuracy
use numerix_dmsl_kit,only:checkGradFast,checkGradFast0,checkGrad
implicit none 
! locals
real(mrk)::dfObs,dfObsB,dfPred
integer(mik)::afcalls,nfigOK(ndim),i,nfigOKmax
real(mrk)::hh(ndim)
logical(mlk)::gradOK
character(100)::status
! local pars
integer(mik),parameter::checkgrad_meth_old=0,checkgrad_meth_new=1
integer(mik),parameter::checkgrad_meth=checkgrad_meth_new
character(*),parameter::fmtNfig="es12.3e3" ! "es23.14e3" ! "es17.8e3" ! 
! Start procedure here
selectcase(chkGrd)
case(chkG_dxstp,chkG_hxstp,chkG_fail)
  selectcase(chkGrd)
  case(chkG_dxstp)            ! * in direction of last step (this becomes rather
      hh=dx                   !   un-informative when sliding along boundaries
  case(chkG_hxstp,chkG_fail)  ! * in default direction
    selectcase(gmeth)
    case(user_meth) ! analytical gradient (default estimate of "hx" for checking)
      hh=sqrt(epsF)*max(abs(xopt),xscale)
    case default    ! use existing finite difference perturbations
      hh=getHxFromRelHx(hx,xopt,xscale,FDscale)
    endselect
    if(boundedSearch)then
      where(xopt+hh<xLo.or.xopt+hh>xHi)hh=-hh
    endif
  endselect
  selectcase(checkgrad_meth)
  case(checkgrad_meth_old)  ! old fast method
    call checkGradFast0(evalFunc,dataIN,dataOUT,x=xopt,fx=fopt,grad=gradopt,&
      xdir=hh,h=chkGrd_h,fscale=fscale,scalingAnalysis=.true.,&
      tolG=chkGrd_tG,tolGdf=chkGrd_tGdf,tolF=chkGrd_tF,&
      gradAnalysis=gradCheckAnalysis,&
      dfObs=dfObs,dfObsB=dfObsB,dfPred=dfPred,fcalls=afcalls,&
      err=err,message=message)
  case(checkgrad_meth_new)  ! new fast method (recommended)
    call checkGradFast(evalFunc,dataIN,dataOUT,x=xopt,fx=fopt,grad=gradopt,&
      xdir=hh,h=chkGrd_h,fscale=fscale,&
      tolG=chkGrd_tG,tolGdf=chkGrd_tGdf,tolF=chkGrd_tF,&
      gradAnalysis=gradCheckAnalysis,&
      dfA=dfObs,dfB=dfObsB,dfPred=dfPred,fcalls=afcalls,&
      err=err,message=message)
  endselect
  fcalls=fcalls+afcalls
  if(uout>0)then
    if(err/=0)then;write(status,'(a,i0)')"err:",err
    else;          status="OK";endif
    write(uout,'(a,i7,a,a,i3,3(a,1x,es10.3),a)')&
      "iter=",iter,&
      "  chk Grad... "//trim(status),&
      "; result:",gradCheckAnalysis,&
      "; dfPred=",dfPred,"; dfObs=",dfObs,"; dfObsB=",dfObsB,&
      "; message="//trim(message)
  endif
case(chkG_full,chkG_f2g)  ! full gradient check
  call checkGrad(evalFunc,dataIN,dataOUT,xopt,fopt,gradopt,hx,chkGrd_gmeth,&
                  xscale,epsF,chkGrd_tG,hh,nfigOK,gradOK,afcalls,err,message)
  fcalls=fcalls+afcalls
  if(uout>0)then
    if(err/=0)then;write(status,'(a,i0)')"err:",err
    else;          status="OK";endif
    write(uout,'(a,i7,a,a,i4,2(a,1x,i4),a)')&
      "iter=",iter,&
      "  chk Grad... "//trim(status),&
      "; bestAgree:",maxval(nfigOK),&
      "; worstAgree:",minval(nfigOK),&
      "; message="//trim(message)
    selectcase(iterNfo) ! possibly print entire gradient analysis
    case(iterNfo_var)
      write(uout,'(a)') "-----"
      write(uout,'(a)') "Gradient analysis using "&
                        //merge("Forward O(1) Diffs",&
                                "Central O(2) Diffs",chkGrd_gmeth==fd_gmeth)&
                        //" method"
      do i=1,ndim
        write(uout,'(a,i6,a,2(a,'//fmtNfig//'),a,i4)')&
          "var",i,": ",&              ! variable
          "gradOpt=",  gradopt(i),&   ! supplied gradient
          "; gradFD=",hh(i),&         ! estimated gradient
          "; nFigOK=",nfigOK(i)       ! decimal figures agreement
      enddo
      write(uout,'(a)') "-----"
    endselect
  endif
case default  ! no gradient check
  nfigOK=10
endselect
if(err/=0)then
  message="f-checkGrad_macro/&"//message; return
endif
nfigOKmax=maxval(nfigOK)
if(nfigOKmax<1)then
  err=-10
  write(message,'(a,i0,a)')"f-checkGrad_macro/veryBadGradAccuracy[nfigOKmax=",nfigOKmax,"]"
else
  err=0
endif
! Start procedure here
endsubroutine checkGrad_macro
!-----
subroutine checkHess_macro()  ! macro to check Hessian accuracy
use utilities_dmsl_kit,only:ns=>number_string,write_matrix,arthsi,trimv=>trim,flip_UtoL
use numerix_dmsl_kit,only:checkHess
implicit none 
! locals
integer(mik)::afcalls,agcalls,nfigOK(ndim,ndim),nfigOKmax
real(mrk)::hfd(ndim,ndim)
logical(mlk)::hessOK
character(100)::hmethChar,status
! local pars
! Start procedure here
selectcase(chkHess)
case(chkHess_full,chkHess_f2g)
  selectcase(hmeth)
  case(bfgsFac_hmeth) ! currently cannot check Cholesky decomposition of Hessian
! (need to compute full Hessian, or, if using cheap method, projected Hessian)
    err=100;message="f-checkHess_macro/CholeskyHessian:checkNotSupported"
    return
  endselect
  call flip_UtoL(hessopt)  ! make Hessian symmetric
  call checkHess(evalFunc,dataIN,dataOUT,xopt,fopt,gradopt,hessopt,chkHess_hmeth,&
                 xscale,epsF,chkGrd_tG,hfd,nfigOK,hessOK,afcalls,agcalls,err,message)
  fcalls=fcalls+afcalls; gcalls=gcalls+agcalls
  if(uout>0)then
    if(err/=0)then;write(status,'(a,i0)')"err:",err
    else;          status="OK";endif
    write(uout,'(a,i7,a,a,i3,2(a,1x,i4),a)')&
      "iter=",iter,&
      "  chk Hess... "//trim(status),&
      "; bestAgree:",maxval(nfigOK),&
      "; worstAgree:",minval(nfigOK),&
      "; message="//trim(message)
    selectcase(iterNfo) ! possibly print entire Hessian analysis
    case(iterNfo_var)
      selectcase(chkHess_hmeth)
      case(fdg_hmeth)
        hmethChar="'Gradient differencing, one-sided, O(1)'"
      case(cdg_hmeth)
        hmethChar="'Gradient differencing, central, O(2)'"
      case(fdf_hmeth)
        hmethChar="'Function differencing, one-sided, O(1)'"
      case(cdf_hmeth)
        hmethChar="'Function differencing, central, O(2)'"
      case default
        hmethChar="'chkHess_hmeth=Unknown' option, probably user input error"
      endselect
      write(uout,'(a)') "-----"
      write(uout,'(a)') "Hessian analysis using "//trim(hmethChar)//" method"
      call write_matrix(unt=uout,header="Supplied Hessian",&
        m=hessopt,nfig=4,display=-1,vLabel="var"//trimv(ns(arthsi(ndim)),3),&
        err=err,message=message)
      write(uout,'(a)') " "
      call write_matrix(unt=uout,header="Approxim Hessian",&
        m=hfd,nfig=4,display=-1,vLabel="var"//trimv(ns(arthsi(ndim)),3),&
        err=err,message=message)
      write(uout,'(a)') " "
      call write_matrix(unt=uout,header="Decimal digits of agreement",&
        m=nfigOK,display=-1,vLabel="var"//trimv(ns(arthsi(ndim)),3),&
        err=err,message=message)
      write(uout,'(a)') "-----"
    endselect
  endif
case default  ! no Hessian check
  err=0;message="w-checkHess_macro/hessCheckNotCarriedOut"
  nfigOK=10
endselect
if(err/=0)then
  message="f-checkHess_macro/&"//message; return
endif
nfigOKmax=maxval(nfigOK)
if(nfigOKmax<1)then
  err=merge(0,-10,ignoreBadHess)
  write(message,'(a,i0,a)')"f-checkHess_macro/veryBadHessAccuracy[nfigOKmax=",nfigOKmax,"]"
else
  err=0
endif
! Start procedure here
endsubroutine checkHess_macro
!-----
endsubroutine qnewton
!----------------------------------------------------
pure subroutine checkStepBounds(x,xLo,xHi,activeSet,dx,stepToBound,hitBound)
! When carrying out box-constrained optimisation, the natural suggested Newton
! step may be too long and must be truncated. However, it is too early to freeze
! any variables since the final step may be even shorter.
! Optionally returns largest step to nearest bound (stepToBound)
! Comments
! - Method may fail (overflow) on numerical condition if x~xBound~0~dx.
!   Need scale to handle this safely.
use utilities_dmsl_kit,only:zero,one
implicit none
! dummies
real(mrk),intent(in)::x(:),xLo(:),xHi(:)
integer(mik),intent(in)::activeSet(:)
real(mrk),intent(inout)::dx(:)
real(mrk),intent(out),optional::stepToBound
logical(mlk),optional,intent(out)::hitBound
! locals
real(mrk)::boundLmax,boundLj,dxx
integer(mik)::j,jMax
! Start procedure here
boundLmax=hugeRe;jMax=0
do j=1,size(dx)
  selectcase(activeSet(j))
  case(freeVar_as)                ! * free variable
    if(dx(j)>zero)then      !   - check upper bound
      dxx=max(epsRe*max(abs(xHi(j)),abs(x(j))),abs(dx(j)))
      boundLj=(xHi(j)-x(j))/dxx
      if(boundLj<boundLmax)then
        boundLmax=boundLj; jMax=j
      endif
    elseif(dx(j)<zero)then  !   - check lower bound
      dxx=max(epsRe*max(abs(xLo(j)),abs(x(j))),abs(dx(j)))
      boundLj=(x(j)-xLo(j))/dxx
      if(boundLj<boundLmax)then
        boundLmax=boundLj; jMax=j
      endif
    endif
  case(loVar_as,hiVar_as,&        ! * fixed variable
       freeLoVar_as,freeHiVar_as) ! * fixed variable with negative Lagrange multiplier
    dx(j)=zero  ! zero increment for frozen variables (should normally be zero already)
  endselect
enddo
if(boundLmax<=one)then  ! step crashes into bound and should be truncated
  dx=dx*boundLmax
  if(present(hitBound))hitBound=.true.
  if(present(stepToBound))stepToBound=one
else                    ! step does not reach bound
  if(present(hitBound))hitBound=.false.
  if(present(stepToBound))stepToBound=boundLmax
endif
! End procedure here
endsubroutine checkStepBounds
!----------------------------------------------------
pure subroutine truncStepCompsBounds(x,xLo,xHi,activeSet,dx,newActiveSet)
! When carrying out box-constrained optimisation, the natural step may be too long
! and must be truncated. In some cases, only the "illegal" values need to be zeroed.
! This sub is currently unused.
use utilities_dmsl_kit,only:zero,one
implicit none
! dummies
real(mrk),intent(in)::x(:),xLo(:),xHi(:)
integer(mik),intent(in)::activeSet(:)
integer(mik),intent(out),optional::newactiveSet(:)
real(mrk),intent(inout)::dx(:)
! locals
integer(mik)::j
! Start procedure here
if(present(newActiveSet))newActiveSet=activeSet
do j=1,size(dx)
  selectcase(activeSet(j))
  case(freeVar_as)                ! * free variable
    if(x(j)+dx(j)>xHi(j))then       !   - check upper bound
      dx(j)=max(zero,xHi(j)-x(j))   !     (guarding against roundoff)
      if(present(newActiveSet))newActiveSet(j)=hiVar_as
    elseif(x(j)+dx(j)<xLo(j))then   !   - check lower bound
      dx(j)=min(zero,xLo(j)-x(j))
      if(present(newActiveSet))newActiveSet(j)=loVar_as
    endif
  case(loVar_as,hiVar_as,&        ! * fixed variable
       freeLoVar_as,freeHiVar_as) ! * fixed variable with negative Lagrange multiplier
    dx(j)=zero  ! zero increment for frozen variables (should normally be zero already)
  endselect
enddo
! End procedure here
endsubroutine truncStepCompsBounds
!----------------------------------------------------
pure subroutine checkActiveSet(x,xscale,activeSet,fixDiagOption,&
    hmeth,grad,hess,Ld,xLo,xHi,onBound,nfree,nfix,nthawn)
! Purpose: check active set, freezing/thawing variables and adjusting Hessian.
! The projected Hessian is formed by zeroing the rows and columns corresponding
! to the frozen/thawn variable (except for the diagonal, which remains non-zero).
! fixDiagOption chooses what to do with Hessian diagonals of fixed/thawed variables:
! setting them to unity or keeping existing value (for better Hessian scaling)
! Comments:
! A. Moving onto a bound when using factored Hessians is slightly tricky, since
!    it is WRONG just to zero a portion of the Cholesky factor. Instead need to
!    use the updating method of Gill and Murray (1976) NPL Report NAC72.
! B. The implementation of the factored BFGS Hessian requires the diagonal of
!    matrix hess to contain the Cholesky diagonal Ld (duplicating the latter).
!    In later parts of the code, when computing expected model reduction
!    the call to quadDf allows supplying a Cholesky factor in place of the Hessian,
!    but it must be contained within the matrix (not in separate array).
!    Therefore diagonal of hess must always stay synchronized with Ld.
!    All Cholesky routines in DMSL handle the Cholesky in the lower triangle
!    and do not modify the upper triangle including diagonal (this allows
!    storing another symmetric or triangular matrix there).
use utilities_dmsl_kit,only:zero,one,replaceRowColMat,putDiag
use linalg_dmsl_kit,only:choles_update
implicit none
! dummies
real(mrk),intent(in)::x(:),xscale(:),grad(:),xLo(:),xHi(:)
integer(mik),intent(in)::hmeth
integer(mik),intent(in)::fixDiagOption
real(mrk),intent(inout)::hess(:,:),Ld(:)
integer(mik),intent(inout)::activeSet(:)
logical(mlk),intent(out)::onBound
integer(mik),intent(out)::nfree,nfix,nthawn
! locals
integer(mik)::j,jerr
character(100)::jmsg
real(mrk),parameter::safeEps=10._mrk*epsRe
logical(mlk)::actv(size(Ld))
! Start procedure here
onBound=.false.; actv=activeSet==freeVar_as
do j=1,size(x)
  selectcase(activeSet(j))
  case(freeVar_as)                ! * free variable
    if    (x(j)<=xLo(j)+safeEps*max(abs(xLo(j)),xscale(j)))then ! freeze variable
      onBound=.true.    ! now at lower bound, zero Hessian entries
      selectcase(hmeth) ! update Cholesky factors
      case(bfgsFac_hmeth)
        call choles_update(L=hess,Ld=Ld,useLDL=.false.,actvrc=actv,&
                           irc=j,err=jerr,message=jmsg)
        call putDiag(hess,Ld)       ! See comment B/C
      endselect
      selectcase(fixDiagOption)
      case(setUnit_fixDiag)                     ! set scaled unit diagonal
        call replaceRowColMat(hess,j,newDiag=one/xscale(j)**2)
        selectcase(hmeth)
        case(bfgsFac_hmeth)
          Ld(j)=one/xscale(j); hess(j,j)=Ld(j)  ! See comment B/C
        endselect
      case(keepDiag_fixDiag)
        call replaceRowColMat(hess,j)
      endselect
      activeSet(j)=merge(loVar_as,freeLoVar_as,grad(j)>zero)
      actv(j)=.false.
    elseif(x(j)>=xHi(j)-safeEps*max(abs(xHi(j)),xscale(j)))then ! freeze variable
      onBound=.true.    ! now at upper bound, zero Hessian entries
      selectcase(hmeth) ! update Cholesky factors
      case(bfgsFac_hmeth)
        call choles_update(L=hess,Ld=Ld,useLDL=.false.,actvrc=actv,&
                           irc=j,err=jerr,message=jmsg)
        call putDiag(hess,Ld)       ! See comment B/C
      endselect
      selectcase(fixDiagOption)
      case(setUnit_fixDiag)
        call replaceRowColMat(hess,j,newDiag=one/xscale(j)**2)
        selectcase(hmeth)
        case(bfgsFac_hmeth)
          Ld(j)=one/xscale(j); hess(j,j)=Ld(j)  ! See comment B/C
        endselect
      case(keepDiag_fixDiag)
        call replaceRowColMat(hess,j)
      endselect
      activeSet(j)=merge(hiVar_as,freeHiVar_as,grad(j)<zero)
      actv(j)=.false.
    endif
  case(loVar_as)                    ! * variable at lower bound
    if(grad(j)<=zero)&
      activeSet(j)=freeLoVar_as     ! ... thaw variable and
    selectcase(fixDiagOption)       ! ... ensure Hessian entries remained zeroed
    case(setUnit_fixDiag)
      call replaceRowColMat(hess,j,newDiag=one/xscale(j)**2)
      selectcase(hmeth)
      case(bfgsFac_hmeth)
        Ld(j)=one/xscale(j); hess(j,j)=Ld(j)  ! See comment B/C
      endselect
    case(keepDiag_fixDiag)
      call replaceRowColMat(hess,j)
    endselect
  case(freeLoVar_as)                ! * semi-free variable at lower bound
    if(grad(j)>zero)&
      activeSet(j)=loVar_as         ! ... fix variable and
    selectcase(fixDiagOption)       ! ... ensure Hessian entries remained zeroed
    case(setUnit_fixDiag)
      call replaceRowColMat(hess,j,newDiag=one/xscale(j)**2)
      selectcase(hmeth)
      case(bfgsFac_hmeth)
        Ld(j)=one/xscale(j); hess(j,j)=Ld(j)  ! See comment B/C
      endselect
    case(keepDiag_fixDiag)
      call replaceRowColMat(hess,j)
    endselect
  case(hiVar_as)                    ! * variable at higher bound
    if(grad(j)>=zero)&
      activeSet(j)=freeHiVar_as     ! ... thaw variable and
    selectcase(fixDiagOption)       ! ... ensure Hessian entries remained zeroed
    case(setUnit_fixDiag)
      call replaceRowColMat(hess,j,newDiag=one/xscale(j)**2)
      selectcase(hmeth)
      case(bfgsFac_hmeth)
        Ld(j)=one/xscale(j); hess(j,j)=Ld(j)  ! See comment B/C
      endselect
    case(keepDiag_fixDiag)
      call replaceRowColMat(hess,j)
    endselect
  case(freeHiVar_as)                ! * semi-free variable at higher bound
    if(grad(j)<zero)&
      activeSet(j)=hiVar_as         ! ... fix variable and
    selectcase(fixDiagOption)       ! ... ensure Hessian entries remained zeroed
    case(setUnit_fixDiag)
      call replaceRowColMat(hess,j,newDiag=one/xscale(j)**2)
      selectcase(hmeth)
      case(bfgsFac_hmeth)
        Ld(j)=one/xscale(j); hess(j,j)=Ld(j)  ! See comment B/C
      endselect
    case(keepDiag_fixDiag)
      call replaceRowColMat(hess,j)
    endselect
  endselect
enddo
nfree= count(activeset==freeVar_as)
nfix=  count(activeset==loVar_as.or.activeset==hiVar_as)
nthawn=count(activeset==freeLoVar_as.or.activeset==freeHiVar_as)
! End procedure here
endsubroutine checkActiveSet
!----------------------------------------------------
subroutine checkReleaseActiveSet(x,xscale,activeSet,grad,hess,Ld,&
  tolFast,forceRel,tolForce,nfree,nfix,nthawn)
! Purpose: check conditions for release of thawn variables from bounds
! Conditions:
! * Immediate release of var(i) if
!                                  grad(i) >= tolFast  * ||grad(free)||
! * If forceRel=.true. (ie, must release at least one var), release var(i) where
!                                  grad(i) >= tolForce * ||grad(thawn)||
! Actions:
! * Sets the status of variable
! * Could perhaps adjust Hessian components corresponding to the fixed variables
!  (currently retains Hessian as is). The questions is really of scaling.
!  See also usage of "fixDiagOption", which controls whether to overwrite fixed
!  diagonals with unity.
use utilities_dmsl_kit,only:imaxloc
implicit none
! dummies
real(mrk),intent(in)::x(:),xscale(:),grad(:)
logical(mlk),intent(in)::forceRel
real(mrk),intent(in)::tolFast,tolForce
real(mrk),intent(inout)::hess(:,:),Ld(:)
integer(mik),intent(inout)::activeSet(:)
integer(mik),intent(out)::nfree,nfix,nthawn
! locals
logical(mlk)::thawn(size(x))
integer(mik)::imax
real(mrk)::gradMaxFree,gradMaxThawn
! Start procedure here
thawn=activeset==freeLoVar_as.or.activeset==freeHiVar_as
if(count(thawn)>0)then ! something can be released
! * Immediate release of variables
  if(any(activeset==freeVar_as))then
    gradMaxFree=maxval(abs(grad)*max(abs(x),xscale),mask=activeset==freeVar_as)
    where(thawn.and.abs(grad)*max(abs(x),xscale)>=tolFast*gradMaxFree)&
      activeSet=freeVar_as
  endif
  if(forceRel)then
! * Time to release at least one variable
    imax=imaxloc(abs(grad)*max(abs(x),xscale),mask=thawn)
    activeSet(imax)=freeVar_as
    gradMaxThawn=abs(grad(imax))*max(abs(x(imax)),xscale(imax))
    where(thawn.and.abs(grad)*max(abs(x),xscale)>=tolForce*gradMaxThawn)&
      activeSet=freeVar_as  ! allow more than one variable to the released
  endif
endif
nfree= count(activeset==freeVar_as)
nfix=  count(activeset==loVar_as.or.activeset==hiVar_as)
nthawn=count(activeset==freeLoVar_as.or.activeset==freeHiVar_as)
! End procedure here
endsubroutine checkReleaseActiveSet
!----------------------------------------------------
subroutine fdigits2epsF(fdigits,evalFunc,dataIN,dataOUT,x,xLo,xHi,xscale,fscale,Hscale,hammPow,&
                        uout,uoutTitle,epsF,fcalls,err,message)
! Purpose: Converts number of reliable digits to function evaluation precision.
use types_dmsl_kit,only:data_ricz_type
use utilities_dmsl_kit,only:estimateEpsF
implicit none
! dummies
type(data_ricz_type),intent(in),optional::dataIN
type(data_ricz_type),intent(inout),optional::dataOUT
integer(mik),intent(in)::fdigits
real(mrk),intent(in)::x(:),xscale(:),fscale,Hscale,hammPow
real(mrk),optional,intent(in)::xLo(:),xHi(:)
integer(mik),intent(in)::uout
character(*),intent(in),optional::uoutTitle
real(mrk),intent(out)::epsF
integer(mik),intent(out)::fcalls
integer(mik),intent(out)::err
character(*),intent(out)::message
interface
  subroutine evalFunc(dataIN,dataOUT,x,feas,fx,gradFx,hessFx,err,message)
  use kinds_dmsl_kit
  use types_dmsl_kit,only:data_ricz_type
  implicit none
  type(data_ricz_type),intent(in),optional::dataIN
  type(data_ricz_type),intent(inout),optional::dataOUT
  real(mrk),intent(in)::x(:)
  logical(mlk),intent(out)::feas
  real(mrk),intent(out),optional::fx,gradFx(:),hessFx(:,:)
  integer(mik),intent(out)::err
  character(*),intent(out)::message
  endsubroutine evalfunc
endinterface
! locals
logical(mlk)::ok
real(mrk)::eA,eAfast,epsFfast
logical(mlk),parameter::fastOnly=.false.,useFast=.false.
! Start procedure here
fcalls=0;err=0
selectcase(fdigits)
case(-2)      ! use method of Hamming to estimate fdigits (function accuracy)
  call estimateEpsF(evalFunc,dataIN,dataOUT,x,xLo,xHi,xscale,fscale,Hscale,hammPow,fastOnly,&
    uout,uoutTitle,eA=eA,epsF=epsF,eAfast=eAfast,epsFfast=epsFfast,&
    feas=ok,fcalls=fcalls,err=err,message=message)
  if(.not.ok.or.err/=0)then
    err=20;   message="f-fdigits2epsF/&"//message
  else
    err=0;    message="fdigits2epsF/ok"
  endif
  if(useFast)epsF=epsFfast
case(-1)      ! full precision
  epsF=epsRe; message="fdigits2epsF/fullPrecision"
case(0:2)     ! virtually unworkable precision for optimisation methods
  err=30;     message="f-fdigits2epsF/fdigitsTooLow"
case(3:)      ! user-supplied relative precision
  epsF=10._mrk**(-fdigits); message="fdigits2epsF/ok(base10)"
case default  ! unknown specification
  err=-100;   message="f-fdigits2epsF/fdigits:unknown"
endselect
! End procedure here
endsubroutine fdigits2epsF
!----------------------------------------------------
pure function get_gOh_fdcd(gradFD,d2fdx2,h,activeset,tol)
! Purpose: calculates fraction where the estimated truncation error
! of a forward difference derivative exceeds the threshold
use utilities_dmsl_kit,only:half,zero
implicit none
! dummies
real(mrk),intent(in)::gradFD(:),d2fdx2(:),h(:),tol
integer(mik),optional,intent(in)::activeset(:)
real(mrk)::get_gOh_fdcd
! locals
integer(mik)::ndim
logical(mlk)::active(size(gradFD))
! Start procedure here
if(present(activeSet))then
  active=(activeSet==freeVar_as)
  ndim=count(active)
  if(ndim>0)then  ! * at least one active variable
    get_gOh_fdcd=&
      real(count(half*abs(h)*abs(d2fdx2)>tol*abs(gradFD).and.active),mrk)/&
           real(ndim,mrk)
  else            ! * all variables fixed
    get_gOh_fdcd=zero
  endif
else              ! * all variables active
    ndim=size(gradFD)
    get_gOh_fdcd=&
      real(count(half*abs(h)*abs(d2fdx2)>tol*abs(gradFD)),mrk)/real(ndim,mrk)
endif
! End procedure here
endfunction get_gOh_fdcd
!----------------------------------------------------
pure function getStepLen2(dx,xscale)
! Purpose: Computes scaled norm-2 steplength. Used in assessing stepmax
use utilities_dmsl_kit,only:norm2
implicit none
! dummies
real(mrk),intent(in)::dx(:),xscale(:)
real(mrk)::getStepLen2
! Start procedure here
getStepLen2=norm2(dx/xscale)
! End procedure here
endfunction getStepLen2
!----------------------------------------------------
pure function scaledGrad(grad,x,fx,xscale,fscale,activeSet,incAllFree)
! Purpose: Computes scaled gradient at point "x".
! Optional activeSet and incAllFree allows to regulate which variables
! included in analysis.
use utilities_dmsl_kit,only:zero
implicit none
! dummies
real(mrk),intent(in)::grad(:),x(:),fx,xscale(:),fscale
integer(mik),intent(in),optional::activeSet(:)
logical(mlk),intent(in),optional::incAllFree
real(mrk)::scaledGrad
! locals
logical(mlk)::active(size(grad))
! Start procedure here
if(present(activeSet))then
  if(present(incAllFree))then
    if(incAllFree)then
      active= activeSet==freeVar_as   .or.&
              activeSet==freeLoVar_as .or.&
              activeSet==freeHiVar_as
    else
      active= activeSet==freeVar_as
    endif
  else
      active= activeSet==freeVar_as
  endif
  if(any(active))then   ! * at least one active variable
    scaledGrad=&
          maxval(abs(grad)*max(abs(x),xscale)/max(abs(fx),fscale),&
          mask=active)
  else                  ! * all variables fixed
    scaledGrad=zero
  endif
else                    ! * all variables active
    scaledGrad=&
          maxval(abs(grad)*max(abs(x),xscale)/max(abs(fx),fscale))
endif
! End procedure here
endfunction scaledGrad
!----------------------------------------------------
pure function scaledStepLen(dx,x,xscale)
! Purpose: Computes scaled steplength
implicit none
! dummies
real(mrk),intent(in)::dx(:),x(:),xscale(:)
real(mrk)::scaledStepLen
! Start procedure here
scaledStepLen=maxval(abs(dx)/max(abs(x),xscale))
! End procedure here
endfunction scaledStepLen
!----------------------------------------------------
pure function scaledFred(fred,fx,fscale)
! Purpose: Computes scaled function reduction
implicit none
! dummies
real(mrk),intent(in)::fred,fx,fscale
real(mrk)::scaledFred
! Start procedure here
scaledFred=fred/max(abs(fx),fscale)
! End procedure here
endfunction scaledFred
!----------------------------------------------------
pure subroutine checkConvergence0(x,fx,gradFx,activeSet,xscale,fscale,gtol,termcode)
! Purpose: check convergence of optimisation algorithm at initial point using
! (a) max-norm of scaled gradient;
! Comments:
! * The gradient tolerance gtol supplied to this procedure should be very stringent
!   to avoid spurious termination on the startinhg point.
! * Note this procedure does not request CD gradient approximations. This can
!   be requested by the calling program
use utilities_dmsl_kit,only:one
implicit none
! dummies
real(mrk),intent(in)::x(:),fx,gradFx(:),xscale(:),fscale,gtol
integer(mik),intent(in),optional::activeSet(:)   ! current active set
integer(mik),intent(out)::termcode
! Start procedure here
termcode=no_con
if(scaledGrad(gradFx,x,fx,xscale,fscale,activeSet,.true.)<=gtol)then
  termcode=grad_con
endif
! End procedure here
endsubroutine checkConvergence0
!----------------------------------------------------
pure subroutine checkConvergence(x,dx,fx,gradFx,activeSet,gmeth,fredExp,fredAct,&
  xscale,fscale,gtol,stol,ftol,skipDxDfCheck,termcode)
! Purpose: check convergence of optimisation algorithm using
! (a) max-norm of scaled gradient;
! (b) scaled step tolerance
! (c) expected and predicted function reduction
! Will not terminate happily unless gradient is no larger than gtolMin
!    (i) user-provided; or
!   (ii) approximated using central differences;
use utilities_dmsl_kit,only:one
implicit none
! dummies
real(mrk),intent(in)::x(:),dx(:),fx,gradFx(:) ! current point properties
integer(mik),intent(in),optional::activeSet(:)! current active set
integer(mik),intent(in)::gmeth                ! method of gradient evaluation
real(mrk),intent(in)::fredExp,fredAct         ! expected actual and reduction
real(mrk),intent(in)::xscale(:),fscale        ! scale modifiers
real(mrk),intent(in)::gtol,stol,ftol          ! convergence tolerances
logical(mlk),intent(in)::skipDxDfCheck        ! requests skipping dx and df checks
integer(mik),intent(out)::termcode            ! termination code
! locals
real(mrk)::scaledG,scaledS,scaledFobs,scaledFexp
real(mrk),parameter::gtolMin=0.1_mrk  ! minimal gradient for succesful termination
! Start procedure here
termcode=no_con
scaledG=scaledGrad(gradFx,x,fx,xscale,fscale,activeSet,.true.)
if(scaledG<=gtol)then
! * check gradient convergence
  selectcase(gmeth)
  case(fd_gmeth)
!    termcode=switchCD_con
    termcode=-grad_con
  case default
    termcode=grad_con         ! gradient criterion satisfied
  endselect
elseif(.not.skipDxDfCheck)then
! check steplength tolerance and function convergence tolerance.
! these should be skipped if a bound has been hit (since steplength can be very small)
! or if releasing variables from the freezing set (in which case the only reliable
! measure of convergence is the gradient.
  scaledS=scaledStepLen(dx,x,xscale)
  scaledFobs=scaledFred(abs(fredAct),fx,fscale)
  scaledFexp=scaledFred(abs(fredExp),fx,fscale)
  if(scaledS<=stol.and.scaledG<=gtolMin)then
! * check step convergence
    selectcase(gmeth)
    case(fd_gmeth)
!      termcode=switchCD_con
      termcode=-search_con
    case default
      termcode=search_con
    endselect
  elseif(scaledFobs<=ftol.and.scaledFexp<=ftol.and.scaledG<=gtolMin)then
! * check function convergence (expected and actual)
    selectcase(gmeth)
    case(fd_gmeth)
!      termcode=switchCD_con
      termcode=-fred_con
    case default
      termcode=fred_con
    endselect
  elseif(scaledS<=stol)then   ! iterates converged but gradient too large
    selectcase(gmeth)
    case(fd_gmeth)
!      termcode=switchCD_con
      termcode=-srchBadGrad_con
    case default
      termcode=srchBadGrad_con
    endselect
  elseif(scaledFobs<=ftol.and.scaledFexp<=ftol)then ! function converged
!  elseif(scaledFobs<=ftol)then ! function converged but gradient large
    selectcase(gmeth)
    case(fd_gmeth)
!      termcode=switchCD_con
      termcode=-fredBadGrad_con
    case default
      termcode=fredBadGrad_con
    endselect
  endif
endif
! End procedure here
endsubroutine checkConvergence
!----------------------------------------------------
pure subroutine makeGoodHessDiag(hdiag,controlHessCond,maxHessCond)
! Purpose: Makes a good Hessian diagonal suitable for use in the quasi-Newton
! method. Must be positive non-singular with moderate conditioning
use utilities_dmsl_kit,only:zero,one
implicit none
! dummies
real(mrk),intent(inout)::hdiag(:)
logical(mlk),intent(in)::controlHessCond
real(mrk),intent(in)::maxHessCond
! locals
real(mrk)::dMax
! Start procedure here
hdiag=abs(hdiag); dMax=maxval(hdiag)
if(dMax==zero)then        ! handle zero case
  hdiag=one
elseif(controlHessCond)then   ! control conditioning
  where(hdiag<maxHessCond*dMax)hdiag=maxHessCond*dMax
endif
! End procedure here
endsubroutine makeGoodHessDiag
!----------------------------------------------------
pure subroutine initQHess_inv(fx,fscale,xscale,himeth,hessinv,hdiag)
! Purpose: Initialisation of inverse quasi-Hessian.
! Multiple options available, including inv[b*I], where
! I is the identity matrix and b is a scale factor chosen following
! Algorithm 9.4.3 in Dennis and Schnabel.
use utilities_dmsl_kit,only:zero,one,putDiag
implicit none
! dummies
real(mrk),intent(in)::fx,fscale,xscale(:)
integer(mik),intent(in)::himeth
real(mrk),intent(in),optional::hdiag(:)
real(mrk),intent(out)::hessinv(:,:)
! locals
real(mrk)::fscale0,hmax,hmin
real(mrk),parameter::safe=1.e2_mrk
! Start procedure here
selectcase(himeth)
case(unt_himeth,untcnd1_himeth)
  hessinv=zero
  call putDiag(hessinv,xscale**2)
case(d2fdx2_himeth)
  hessinv=zero
  hmax=maxval(abs(hdiag))  ! safeguard division by zero
  if(hmax==zero)hmax=one; hmin=hmax*safe*epsRe
  call putDiag(hessinv,one/sign(max(abs(hdiag),hmin),hdiag))
case(scld_himeth,scldcnd1_himeth)
  hessinv=zero; fscale0=max(abs(fx),fscale)
  call putDiag(hessinv,xscale**2/fscale0)
endselect
! End procedure here
endsubroutine initQHess_inv
!----------------------------------------------------
pure subroutine initQHess_unfac(fx,fscale,xscale,himeth,hess,hdiag)
! Purpose: Initialisation of unfactored quasi-Hessian.
! Multiple options available, including b*I, where
! I is the identity matrix and b is a scale factor chosen following
! Algorithm 9.4.3 in Dennis and Schnabel.
use utilities_dmsl_kit,only:zero,one,putDiag
implicit none
! dummies
real(mrk),intent(in)::fx,fscale,xscale(:)
integer(mik),intent(in)::himeth
real(mrk),intent(in),optional::hdiag(:)
real(mrk),intent(out)::hess(:,:)
! locals
real(mrk)::fscale0
! Start procedure here
selectcase(himeth)
case(unt_himeth,untcnd1_himeth)
  hess=zero
  call putDiag(hess,one/xscale**2)
case(d2fdx2_himeth)
  hess=zero
  call putDiag(hess,hdiag)
case(scld_himeth,scldcnd1_himeth)
  hess=zero; fscale0=max(abs(fx),fscale)
  call putDiag(hess,fscale0/xscale**2)
endselect
! End procedure here
endsubroutine initQHess_unfac
!----------------------------------------------------
pure subroutine initQHess_fac(fx,fscale,xscale,himeth,hess,Ld,facBFGS_getLLt)
! Purpose: Initialisation of factored quasi-Hessian.
! Multiple options available, including b*I, where
! I is the identity matrix and b is a scale factor chosen following
! Algorithm 9.4.4 in Dennis and Schnabel.
use utilities_dmsl_kit,only:zero,one,putDiag
implicit none
! dummies
real(mrk),intent(in)::fx,fscale,xscale(:)
integer(mik),intent(in)::himeth
real(mrk),intent(out)::hess(:,:),Ld(:)
logical(mlk),intent(in)::facBFGS_getLLt
! locals
real(mrk)::fscale0
! Start procedure here
selectcase(himeth)
case(unt_himeth,untcnd1_himeth)   ! unit matrix
  hess=zero; Ld=one/xscale
  call putDiag(hess,Ld**2)
case(d2fdx2_himeth)               ! given diagonal
  hess=zero; Ld=abs(Ld)
  call putDiag(hess,Ld)
  Ld=sqrt(Ld)
case(scld_himeth,scldcnd1_himeth) ! scaled unit method of Dennis and Schnabel
  hess=zero; fscale0=sqrt(max(abs(fx),fscale))
  Ld=fscale0/xscale
  call putDiag(hess,Ld**2)
endselect
if(.not.facBFGS_getLLt)then ! insert Cholesky diagonal to correctly compute
  call putDiag(hess,Ld)     ! quadratic model reduction from Cholesky factors.
endif
! End procedure here
endsubroutine initQHess_fac
!----------------------------------------------------
pure subroutine setTrustToCauchy(hess,hessScaled,grad,xscale,&
  trustRad,steepStepLen,err,message)
! Purpose: Initialises the trust region to the length of the Cauchy step
! constrained by max length ('trust') equal to scaled gradient length.
! IN: unscaled hessian and gradient
! OUT:trust radius set to scaled Cauchy step
! Programmer: Dmitri Kavetski, 17 January 2004.
use utilities_dmsl_kit,only:norm2,quadform
implicit none
! dummies
real(mrk),intent(in)::hess(:,:),grad(:),xscale(:)
real(mrk),intent(inout)::hessScaled(:,:)
real(mrk),intent(out)::trustRad,steepStepLen
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
integer(mik)::i,j,n,cauchyStepType
real(mrk)::gBg,gradScaled(size(grad)),trustRadTemp
logical(mlk),parameter::useL=.false.,haveGBG=.true.
! Start procedure here
n=size(hess,1)        ! construct scaled Hessian and gradient
forall(i=1:n,j=1:n,j>=i)hessScaled(i,j)=hess(i,j)*xscale(i)*xscale(j)
gradScaled=grad*xscale; steepStepLen=norm2(gradScaled)
! if(steepStepLen<tinyRe)then ! point with zero gradient
!   trustRad=one
!   err=0; message="w-setTrustToCauchy/unexpected[zeroGrad]"; return
! endif
trustRadTemp=steepStepLen    ! very initial trust - scaled gradient step
gBg=quadform(v=gradScaled,mm=hessScaled,typeM="SU")
call getCauchyStep(hess=hessScaled,useL=useL,grad=gradScaled,&
  trustRad=trustRadTemp,normG=steepStepLen,haveGBG=haveGBG,gBg=gBg,&
  cauchyStep=gradScaled,cauchyLen=trustRad,cauchyStepType=cauchyStepType)
err=0
selectcase(cauchyStepType)
case(cauchyInside)    ! - Cauchy step less than steepest descent step
  message="setTrustToCauchy/ok/trustLessThanGrad"
case(cauchyOnBound)   ! - Cauchy step same as steepest descent step (constrained by trust)
  message="setTrustToCauchy/ok/trustSameGrad"
case(cauchyInfin)     ! - Cauchy step infinite, truncated at gradient length
  message="w-setTrustToCauchy/ok?/cauchyInfin(negCurv?)"
case(cauchyZeroGrad)  ! - Cauchy step collapses because grad ~ 0
  message="w-setTrustToCauchy/unusual/cauchyZero(grad~0)"
endselect
! End procedure here
endsubroutine setTrustToCauchy
!----------------------------------------------------
pure subroutine adaptFDgradHx(hx,x,xscale,FDscale,epsFa,hess)
! Purpose: adapts forward difference gradient perturbation using modified
! Gill et al. / Stewarts / DK approach.
! Given d2f/dx2 and epsFa, the optimal stepsize is given by Gill et al (1983) as
! hF=2*sqrt(epsFa/d2fdx2). Earlier, Stewart derived a similar formula and assumed
! the quasi-Hessian provided order-of-magnitude-estimates of d2f/dx2.
! Gill et al noted that their approach can be combined with Stewarts, here DK
! modifies the approach further to avoid excessive hx-variation from step to step
! The new stepsize is the geometric average of GMSW/S and the old stepsize.
! Ref:
! * Gill,P.E.,W.Murray,M.A.Saunders,and M.H.Wright.1983.Computing forward-difference
!   intervals for numerical optimization. SIAM Journal on Scientific and Statistical
!   Computing 4:310-321.
! * Stewart,G.W.1967.A modification of Davidon's minimization method to accept
!   difference approximations of derivatives. Journal of the Association for
!   Computing Machinery 14:72-83.
use utilities_dmsl_kit,only:two,getdiag,getHxFromRelHx,getRelHxFromHx
implicit none
! dummies
real(mrk),intent(inout)::hx(:)
real(mrk),intent(in)::epsFa,hess(:,:),x(:),xscale(:),FDscale
! Start procedure here
hx=getHxFromRelHx(hx,x,xscale,FDscale)
hx=sign(sqrt(abs(hx)*two*sqrt(epsFa/abs(getdiag(hess)))),hx)
hx=getRelHxFromHx(hx,x,xscale,FDscale)
! End procedure here
endsubroutine adaptFDgradHx
!----------------------------------------------------
pure subroutine adaptCDgradHx(hx,x,xscale,FDscale,eTa,eCa)
! Purpose: adapts central difference gradient perturbation using modified
! Curtis and Reid approach.
! IN: hx = old stepsize (relative)
!     x,xscale,FD = used to obtain absolute stepsize
!     eTa = estimated absolute truncation error
!     eCa = estimated absolute condition/rounding error
! OUT:
!     hx = new relative stepsize
use utilities_dmsl_kit,only:two,zero,getHxFromRelHx,getRelHxFromHx
implicit none
! dummies
real(mrk),intent(inout)::hx(:)
real(mrk),intent(in)::x(:),xscale(:),FDscale,eTa,eCa
! locals
real(mrk),parameter::umin=1.e1_mrk,umax=1.e3_mrk,uaim=1.e2_mrk
! Start procedure here
hx=getHxFromRelHx(hx,x,xscale,FDscale)
if(eTa<eCa)then             ! increase to the max
  hx=hx*sqrt(uaim)
elseif(abs(eCa)==zero)then  ! condition error zero, when dfdx==0
  hx=hx*sqrt(uaim)
else                        ! increase using Curtis/Reid approach
  hx=hx*sqrt(uaim*eCa/eTa)
endif
hx=getRelHxFromHx(hx,x,xscale,FDscale)
! End procedure here
endsubroutine adaptCDgradHx
!----------------------------------------------------
subroutine getFDgradHx(evalFunc,dataIN,dataOUT,x,fx,xscale,fscale,activeSet,epsF,Hscale,&
  gradFD_imethod,hx,xLo,xHi,dfdx,d2fdx2,fcalls,retcode,err,message)
! Purpose: Return finite difference interval optimised for forward differences.
! also returns the central difference estimate.
! Programmer: Dmitri Kavetski
! Method: uses Gill and Murray's stepsize optimisation method bason on
! minimising the sum of truncation and condition errors. The engine procedure
! returns central, forward and backward estimates, which are by-products of the
! stepsize selection process. Since central differences are typically more accurate,
! it would make sense to use them when available. Hence this routine returns
! the central difference estimate, rather than forward or backward estimate.
! In any case, the difference between them is very small since the numerical error
! in them is minimised by this procedure.
use types_dmsl_kit,only:data_ricz_type
use utilities_dmsl_kit,only:half,zero,one,getMaxMagn,getHessDiagFromFunc,&
  epsF_to_epsA,ifirstTrueLoc
use numerix_dmsl_kit,only:dfdx_gill,dfdx_sw
implicit none
! dummies
type(data_ricz_type),intent(in),optional::dataIN
type(data_ricz_type),intent(inout),optional::dataOUT
real(mrk),intent(in)::x(:),fx,epsF,xscale(:),fscale,Hscale
integer(mik),intent(in),optional::activeSet(:)
integer(mik),intent(in)::gradFD_imethod
real(mrk),intent(inout)::hx(:)
real(mrk),intent(out)::dfdx(:)
real(mrk),optional,intent(in)::xlo(:),xhi(:)
real(mrk),optional,intent(out)::d2fdx2(:)
integer(mik),intent(out)::fcalls
integer(mik),intent(out)::err(:),retcode
character(*),intent(out)::message(:)
! user-provided function
interface
  subroutine evalFunc(dataIN,dataOUT,x,feas,fx,gradFx,hessFx,err,message)
  use kinds_dmsl_kit
  use types_dmsl_kit,only:data_ricz_type
  implicit none
  type(data_ricz_type),intent(in),optional::dataIN
  type(data_ricz_type),intent(inout),optional::dataOUT
  real(mrk),intent(in)::x(:)
  logical(mlk),intent(out)::feas
  real(mrk),intent(out),optional::fx,gradFx(:),hessFx(:,:)
  integer(mik),intent(out)::err
  character(*),intent(out)::message
  endsubroutine evalfunc
endinterface
! locals
real(mrk)::epsFa
integer(mik),dimension(size(x))::addFcalls,varStatus,whatdfdx
real(mrk),dimension(size(x))::Edfdx,d2fdx2loc,dfdxF,dfdxB,dfdxC,h0,hmax
integer(mik)::lerr,ibad
character(100)::lmessage
! local pars
real(mrk),parameter::hmaxFac=1.e2_mrk,hmaxSafe=0.1_mrk
! local for Gill
integer(mik),parameter::dfdxF_go=1,dfdxB_go=2,dfdxC_go=3,dfdxCFB_go=4
integer(mik),parameter::dfdx_go=dfdxC_go
integer(mik),parameter::Kmax_gill=10          ! max number of Gill calls to optimise hx(i)
real(mrk),parameter::hxMultMax_gill=1.e3_mrk  ! estimated max hx increase in Gill
real(mrk),parameter::h0safety=100._mrk        ! safety factor to avoid overstepping boundaries
! locals for SW
integer(mik),parameter::fcmax_sw=100 ! max number of SW calls to optimise hx(i)
real(mrk),parameter::beta_sw=3._mrk  ! step reduction in SW method
! registered values of whatdfdx (copied from numerix_dmsl_kit,SW codes
integer(mik),parameter::dfdxC2=0,dfdxF1=1,dfdxB1=-1,dfdxF2=2,dfdxB2=-2,&
                        d2fdx2C2=3,d2fdx2F1=4,d2fdx2B1=-4
! other local pars
logical(mlk),parameter::useHxDef_d2fdx2=.true.
! Start procedure here
err=0; h0=hx; fcalls=0; message="getFDgradHx/ok"; retcode=0
epsFa=epsF_to_epsA(epsF,fx,fscale,Hscale) ! absolute function accuracy
if(present(activeSet))then
  varStatus=activeSet
else
  varStatus=freeVar_as
endif
selectcase(gradFD_imethod)
case(gradFD_gill)
! * use Gill et al.'s method
  if(present(xLo).and.present(xHi))then  ! check internal vars for bound proximity
    where    (varStatus==hiVar_as.or.varStatus==freeHiVar_as) ! backward at upper bound
      hmax=half*(xLo-x)
    elsewhere(varStatus==loVar_as.or.varStatus==freeLoVar_as) ! forward at lower bound
      hmax=half*(xHi-x)
    elsewhere(varStatus==freeVar_as)  ! internal: check either side
      hmax=min(x-xLo,xhi-x)
    endwhere
    where(varStatus==freeVar_as.and.abs(hmax)<=abs(hx)*hxMultMax_gill)
      hmax=getMaxMagn(xLo-x,xHi-x)
      varStatus=merge(hiVar_as,loVar_as,hmax<zero)
    endwhere
    hmax=hmaxSafe*hmax
  else
    hmax=hmaxFac*max(abs(x),xscale)
  endif
  call dfdx_gill(evalFunc,dataIN,dataOUT,x,varStatus,fx,xscale,fscale,Kmax_gill,h0,hmax,epsFa,hx,dfdxF,Edfdx,&
                 d2fdx2loc,dfdxB,dfdxC,addFcalls,err,message)
  fcalls=fcalls+sum(addFcalls)
  if(any(err/=0))then
    ibad=ifirstTrueLoc(err/=0); err(1)=err(ibad)
    retcode=bugFail; message(1)="f-getFDgradHx/&"//message(ibad); return
  endif
  selectcase(dfdx_go)
  case(dfdxF_go)   ! forward (original method)
    dfdx=dfdxF
  case(dfdxB_go)   ! backward
    dfdx=dfdxB
  case(dfdxC_go)   ! return central difference since it is available
    dfdx=dfdxC
  case(dfdxCFB_go) ! return an 'improved' central difference
    dfdx=half*(dfdxF+dfdxB)
  endselect
  if(present(d2fdx2))d2fdx2=d2fdx2loc
case(gradFD_sw1)
! * Stepleman and Winarsky method,O(h) forward/backward diffs
  where(varStatus==freeVar_as)                              ! internal
    whatdfdx=dfdxF1 ! central approximation
  elsewhere(varStatus==loVar_as.or.varStatus==freeLoVar_as) ! lower bound
    whatdfdx=dfdxF1 ! forward app.
  elsewhere(varStatus==hiVar_as.or.varStatus==freeHiVar_as) ! upper bound
    whatdfdx=dfdxB1 ! backward app.
  endwhere
  if(present(xLo).and.present(xHi))then
    where    (varStatus==hiVar_as.or.varStatus==freeHiVar_as)  ! backward at upper bound
      hmax=(xLo-x)
    elsewhere(varStatus==loVar_as.or.varStatus==freeLoVar_as)  ! forward at lower bound
      hmax=(xHi-x)
    elsewhere(varStatus==freeVar_as)  ! internal: pick widest separation
      hmax=getMaxMagn(xLo-x,xHi-x)
      whatdfdx=merge(dfdxF1,dfdxB1,hmax>zero)
    endwhere
    hmax=hmaxSafe*hmax                ! and safeguard just in case
  else
    hmax=hmaxFac*max(abs(x),xscale)
  endif
  call dfdx_sw(evalFunc,dataIN,dataOUT,x,whatdfdx,fxin=fx,epsF=epsFa,fcallsmax=fcmax_sw,h0in=h0,&
    betain=spread(beta_sw,1,size(x)),xscale=xscale,fscale=fscale,hmax=hmax,&
    dfdx=dfdx,Edfdx=Edfdx,hopt=hx,&
    fcalls=addFcalls,err=err,message=message)
  fcalls=fcalls+sum(addFcalls)
  if(any(err/=0))then
    ibad=ifirstTrueLoc(err/=0); err(1)=err(ibad)
    retcode=bugFail; message(1)="f-getFDgradHx/&"//message(ibad); return
  endif
  if(present(d2fdx2))then
    call getHessDiagFromFunc(evalFunc,dataIN,dataOUT,x,fx,xscale,epsF,useHxDef=useHxDef_d2fdx2,&
      hmeth=1,hessDiag=d2fdx2,fcalls=addFcalls(1),err=lerr,message=lmessage)
    fcalls=fcalls+addFcalls(1)
    if(lerr/=0)then
      retcode=bugFail; err(1)=10; message(1)=trim(lmessage) !//message(1)
    endif
  endif
case(gradFD_sw2)
! * Stepleman and Winarsky method,O(h2) analysis
  where(varStatus==freeVar_as)                              ! internal
    whatdfdx=dfdxC2 ! central approximation
  elsewhere(varStatus==loVar_as.or.varStatus==freeLoVar_as) ! lower bound
    whatdfdx=dfdxF2 ! forward app.
  elsewhere(varStatus==hiVar_as.or.varStatus==freeHiVar_as) ! upper bound
    whatdfdx=dfdxB2 ! backward app.
  endwhere
  if(present(xLo).and.present(xHi))then
    where    (whatdfdx==dfdxB2)  ! backward at upper bound
      hmax=half*(xLo-x)
    elsewhere(whatdfdx==dfdxF2)  ! forward at lower bound
      hmax=half*(xHi-x)
    elsewhere(whatdfdx==dfdxC2)  ! internal: check either side
      hmax=min(x-xLo,xhi-x)
    endwhere
    where(whatdfdx==dfdxC2.and.abs(hmax)<h0safety*abs(h0)) ! ensure hmax does not
      hmax=getMaxMagn(xLo-x,xHi-x)                         ! reach any bound
      whatdfdx=merge(dfdxF2,dfdxB2,hmax>zero) ! switch to forward/backward method
      hmax=hmaxSafe*hmax                      ! and safeguard just in case
    endwhere
  else
    hmax=hmaxFac*max(abs(x),xscale)
  endif
  call dfdx_sw(evalFunc,dataIN,dataOUT,x,whatdfdx,fxin=fx,epsF=epsFa,fcallsmax=fcmax_sw,h0in=h0,&
    betain=spread(beta_sw,1,size(x)),xscale=xscale,fscale=fscale,hmax=hmax,&
    dfdx=dfdx,Edfdx=Edfdx,dfdxFree=d2fdx2loc,hopt=hx,&
    fcalls=addFcalls,err=err,message=message)
  fcalls=fcalls+sum(addFcalls)
  if(any(err/=0))then
    ibad=ifirstTrueLoc(err/=0); err(1)=err(ibad)
    retcode=bugFail; message(1)="f-getFDgradHx/&"//message(ibad); return
  endif
  if(present(d2fdx2))d2fdx2=d2fdx2loc
case default  ! bug: unknown method
  retcode=bugFail; err=bugFail; message="f-getFDgradHx/unknownMethod"
endselect
! End procedure here
endsubroutine getFDgradHx
!----------------------------------------------------
pure subroutine bfgsInv_update1(dx,xscale,activeSet,grad,gradold,qhessinv,rescale)
! Purpose: BFGS update of inverse quasi-Hessian (NR-based method).
! * Update is skipped if "fac" is not sufficiently positive.
! * Skipping condition 2 (when change in dx expected to be below noise)
!   is not implemented, since quasi-Hessian itself is unavailable.
! * Classic skipping condition requires 'fac>0' to ensure +ve definite q-Hessian.
!   Modified conditions (BFGS damping) not implemented for the inverse-updating.
! * The implementation below requires far fewer matrix multiplies than
!   "bfgsInv_update2" and takes advantage of symmetry.
! * Option available to rescale the initial diagonal Hessian after first
!   iteration but before first update using eqn (8.20) in Nocedal.
!   This can improve the scaling of Hessian for subsequent updates.
! * Routine can work with upper Hessian only. However, for some compilers,
!   the matmul is so fast that it could be preferred over DMSL's symmetric mamtul...
use utilities_dmsl_kit,only:zero,one,norm2,rank1updt,fmatmul_mv,flip_UtoL
implicit none
! dummies
real(mrk),intent(in)::dx(:),xscale(:),grad(:),gradold(:)
integer(mik),intent(in),optional::activeSet(:)
real(mrk),intent(inout)::qhessinv(:,:)
!logical(mlk),intent(in)::skipClassic
logical(mlk),intent(in)::rescale
! locals
real(mrk)::dg(size(dx)),hdg(size(dx)),fac,fad,fae
! Start procedure here
dg=grad-gradold
if(present(activeSet))then  ! need to zero dg for fixed variables
  where(activeSet/=freeVar_as)dg=zero
endif
fac=dot_product(dg,dx)
if(rescale)then ! rescale Hessian during first iteration
  call qhessRescale(qhess=qhessinv,invrs=.true.,dg=dg,fac=fac)
endif
!if(skipClassic)then ! skip update (condition 1), classic, ensures +ve def update
  if(fac<=sqrt(epsRe)*norm2(dx/xscale)*norm2(dg*xscale))return
!else                ! skip update (condition 1), new, can yield indefinite updates
!  if(abs(fac)<=sqrt(epsRe)*norm2(dx/xscale)*norm2(dg*xscale))return
!endif
if(bfgsInvUt)then   ! avoid accessing lower triangle
  hdg=fmatmul_mv(m=qhessinv,v=dg,typeMV="SUV")
else                ! full matmul
  hdg=matmul(qhessinv,dg)
endif
fae=dot_product(dg,hdg)
fac=one/fac; fad=one/fae; dg=fac*dx-fad*hdg  ! vector that makes BFGS different from DFP
call rank1updt(a=qhessinv,facX=fac,x=dx,facY=-fad,y=hdg,facZ=fae,z=dg,symm="U")
if(.not.bfgsInvUt)call flip_UtoL(qhessinv)      ! make symmetric if requested
! End procedure here
endsubroutine bfgsInv_update1
!----------------------------------------------------
pure subroutine bfgsInv_update2(dx,xscale,activeSet,grad,gradold,qhessinv,mtemp,rescale)
! Purpose: BFGS update of inverse quasi-Hessian (eqn 8.16 in Nocedal).
! * Update is skipped if "fac" is not sufficiently positive.
! * Skipping condition 2 (when change in dx expected to be below noise)
!   is not implemented, since quasi-Hessian itself is unavailable.
! * Classic skipping condition requires 'fac>0' to ensure +ve definite q-Hessian.
!   Modified conditions (BFGS damping) not implemented for the inverse-updating.
! * Option available to rescale the initial diagonal Hessian after first
!   iteration but before first update using eqn (8.20) in Nocedal.
!   This can improve the scaling of Hessian for subsequent updates.
! * This routine should be used as backup only - it implements the BFGS
!   equations in a rather cumbersome inefficient manner.
! * Routine works with entire matrix. This makes it not quite compatible
!   with "bfgsInv_update1", which works solely with upper triangle.
use utilities_dmsl_kit,only:zero,one,norm2,rank1updt,addDiag,outerprod
implicit none
! dummies
real(mrk),intent(in)::dx(:),xscale(:),grad(:),gradold(:)
integer(mik),intent(in),optional::activeSet(:)
real(mrk),intent(inout)::qhessinv(:,:),mtemp(:,:)
!logical(mlk),intent(in)::skipClassic
logical(mlk),intent(in)::rescale
! locals
real(mrk)::dg(size(dx)),fac
! Start procedure here
dg=grad-gradold
if(present(activeSet))then  ! need to zero dg for fixed variables
  where(activeSet/=freeVar_as)dg=zero
endif
fac=dot_product(dg,dx)
if(rescale)then ! rescale Hessian during first iteration
  call qhessRescale(qhess=qhessinv,invrs=.true.,dg=dg,fac=fac)
endif
!if(skipClassic)then ! skip update (condition 1), classic, ensures +ve def update
  if(fac<=sqrt(epsRe)*norm2(dx/xscale)*norm2(dg*xscale))return
!else                ! skip update (condition 1), new, can yield indefinite updates
!  if(abs(fac)<=sqrt(epsRe)*norm2(dx/xscale)*norm2(dg*xscale))return
!endif
fac=one/fac  ! this construction is far less efficient than "bfgsInv_update1"
mtemp=-fac*outerprod(dx,dg); call addDiag(mtemp,one) ! cos it needs several full matmul's
qhessinv=matmul(mtemp,qhessinv); qhessinv=matmul(qhessinv,transpose(mtemp))
call rank1updt(a=qhessinv,fac=fac,x=dx,symm="N")
! End procedure here
endsubroutine bfgsInv_update2
!----------------------------------------------------
pure subroutine bfgsUnfac_update(dx,xscale,activeSet,grad,gradold,qhess,tol,&
  skipClassic,dampedBFGS,dampFac,rescale,err,message)
! Purpose: BFGS update of unfactored quasi-Hessian. This allows monitoring the
! condition number of Hessian and ensuring "sufficient" positive definiteness.
! This naive implementation in this procedure leads to O(N3) cost since the
! Cholesky decomposition of the quasi-Hessian needs to be performed at each iteration.
! Comments:
! * Skipping conditions 1 and 2 implemented, to ensure positive definiteness
!   and prevent numerical noise from degrading the quasi-Hessian
! * Allows BFGS damping as described by Nocedal and Wright 1999,p.201&540,
!   to improve the performance in difficult regions where Hessian not +ve definite.
! * Classic skipping condition ensures BFGS Hessian remains positive
!   definite by skipping updates when 'fac~0'. Nocedal and Wright experience
!   (as well as DK's!) suggests that in some cases this forces excessive
!   skipping and inhibits the methods to the point of failure.
!   Damped BFGS handles 'fac~0' in a different way, still ensuring +ve
!   definite Hessians. A more drastic DK change is merely guard overflow and
!   accept indefinite Hessians. Indeed,when using trust-region methods,
!   indefinite q-Hessians can be OK, indeed, desirable. In this case use
!   skipClassic=.false. and dampedBFGS=.false.
! * Option available to rescale the initial diagonal Hessian after first
!   iteration but before first update using eqn (8.20) in Nocedal.
!   This can improve the scaling of Hessian for subsequent updates.
! * Routine works with upper triangle of Hessian only.
use utilities_dmsl_kit,only:zero,one,norm2,fmatmul_mv,rank1updt
implicit none
! dummies
real(mrk),intent(in)::dx(:),xscale(:),grad(:),gradold(:),tol
integer(mik),intent(in),optional::activeSet(:)
real(mrk),intent(inout)::qhess(:,:)
logical(mlk),intent(in)::skipClassic
logical(mlk),intent(in)::dampedBFGS
real(mrk),intent(in)::dampFac
logical(mlk),intent(in)::rescale
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
real(mrk)::dg(size(dx)),hdx(size(dx)),fac,fad
! BFGS damping
real(mrk)::dampTheta
! Start procedure here
err=0;message="bfgsUnFac_update/ok"; dg=grad-gradold
if(present(activeSet))then  ! need to zero dg for fixed variables
  where(activeSet/=freeVar_as)dg=zero
endif
fac=dot_product(dg,dx)
if(rescale)then ! rescale Hessian during first iteration
  call qhessRescale(qhess=qhess,dg=dg,invrs=.false.,fac=fac)
endif
!hdx=matmul(qhess,dx)
hdx=fmatmul_mv(m=qhess,v=dx,typeMV="SUV")   ! avoid accessing lower triangle
fad=dot_product(dx,hdx)
if(all(abs(dg-Hdx)<=tol*max(abs(grad),abs(gradold))))then   ! prevents noisy updates
  err=0;message="bfgsUnFac_update/skipCond2(noisyUpdate)"
  return
endif
if(skipClassic)then ! - classical (Dennis and Schnabel) skipping conditions
  if(fac<=sqrt(epsRe)*norm2(dx/xscale)*norm2(dg*xscale))then  ! prevents indefinite updates
    err=0;message="bfgsFac_update/skipCond1(posDefUpdate)(classic)"
    return
  endif
elseif(dampedBFGS.and.fac<dampFac*fad)then  ! - BFGS damping
  dampTheta=(one-dampFac)*fad/(fad-fac)
  dg=dampTheta*dg+(one-dampTheta)*Hdx       ! r in Nocedal eq,(18.23)
  fac=dot_product(dg,dx)
elseif(.not.dampedBFGS)then  ! - minimal skipping condition required for stability
  if(abs(fac)<=sqrt(epsRe)*norm2(dx/xscale)*norm2(dg*xscale))then
    err=0;message="bfgsFac_update/skipCond1(minimal)"
    return
  endif
endif
fac=one/fac; fad=one/fad
call rank1updt(a=qhess,facX=-fad,x=hdx,facY=fac,y=dg,symm="U")
! End procedure here
endsubroutine bfgsUnfac_update
!----------------------------------------------------
pure subroutine bfgsFac_update(dx,xscale,activeSet,grad,gradold,Lqhess,Ld,tol,&
  controlHessCond,maxHessCond,logdet,condest,facBFGS_useR2,facBFGS_getLLt,&
  skipClassic,dampedBFGS,dampFac,rescale,err,message)
! Purpose: BFGS update of factored quasi-Hessian.
! The Cholesky factor of H, input as Lqhess and diagonal Ld, is updated
! into the Cholesky factor of the updated Hessian using either:
! (a) a single rank-2 update using QR machinery;
! (b) two consecutive rank-1 updates using Cholesky updates;
! Comments:
! * Skipping conditions 1 and 2 implemented, to ensure positive definiteness
!   and prevent numerical noise from degrading the quasi-Hessian.
! * Allows BFGS damping as described by Nocedal and Wright 1999,p.201&540,
!   to improve the performance in difficult regions where Hessian not +ve definite.
! * Classic skipping condition ensures BFGS Hessian remains positive
!   definite by skipping updates when 'fac~0'. Nocedal and Wright experience
!   (as well as DK's!) suggests that in some cases this forces excessive
!   skipping and inhibits the methods to the point of failure.
!   Damped BFGS handles 'fac~0' in a different way, still ensuring +ve
!   definite Hessians. A more drastic DK change is merely guard overflow and
!   accept indefinite Hessians. Indeed,when using trust-region methods,
!   indefinite q-Hessians can be OK, indeed, desirable. In this case use
!   skipClassic=.false. and dampedBFGS=.false.
! * Option available to rescale the initial diagonal Hessian after first
!   iteration but before first update using eqn (8.20) in Nocedal.
!   This can improve the scaling of Hessian for subsequent updates.
! * Routine works with:
!   - lower triangle Lqhess + diagonal Ld contain Cholesky factor of Hessian.
!   - if(facBFGS_getLLt) upper triangle Lqhess contains Hessian itself (debugging).
! * Since Cholesky factors tend to demonstrate the condition of the matrix
!   (ratio of largest/smallest Cholesky diagonals), it is often beneficial to bound
!   the ratio using maxHessCond. This (trivially) tends to avoid quasi-Newton steps that
!   are nearly-orthogonal to the gradient (and hence delay global convergence)
!   This conditioning is done (implicitly) using the robust Cholesky factorization,
!   but the latter is not needed when Cholesky factors of the BFGS Hessian are evolved.
! Ref
! Dennis and Schnabel (1996) Numerical methods for unconstrained optimization
! and nonlinear equations, Algorithm A.9.4.2, p.356.
use utilities_dmsl_kit,only:zero,one,norm2,fmatmul_mv,putdiag,triang_minEig
use linalg_dmsl_kit,only:choles_update
implicit none
! dummies
real(mrk),intent(in)::dx(:),xscale(:),grad(:),gradold(:),tol
logical(mlk),intent(in)::controlHessCond
real(mrk),intent(in)::maxHessCond
integer(mik),intent(in),optional::activeSet(:)
real(mrk),intent(inout)::Lqhess(:,:),Ld(:),logdet,condest
logical(mlk),intent(in)::facBFGS_getLLt,facBFGS_useR2,skipClassic
logical(mlk),intent(in)::dampedBFGS
real(mrk),intent(in)::dampFac
logical(mlk),intent(in)::rescale
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
real(mrk)::dg(size(dx)),Hdx(size(dx)),u(size(dx)),fac,fad,alpha,beta,Ldmin
character(10)::jmsg
! BFGS damping
real(mrk)::dampTheta
! local pars
logical(mlk),parameter::forceC2inA=.true. ! robust (slow) Cholesky update
! Start procedure here
err=0;message="bfgsFac_update/ok"
dg=grad-gradold
if(present(activeSet))then  ! need to zero dg for fixed variables
  where(activeSet/=freeVar_as)dg=zero
endif
fac=dot_product(dg,dx) !DK: in some cases seems to become -ve, which can affect rescaling
if(rescale)then ! rescale Hessian during first iteration
  call qhessRescale(qhess=Lqhess,Ld=Ld,invrs=.false.,dg=dg,fac=fac)
endif
if(facBFGS_getLLt.and.facBFGS_useR2)then  ! - flag error as u is not formed
  err=10;message="f-bfgsFac_update/unsupported:facBFGS_getLLt.and.facBFGS_useR2"
  return
elseif(facBFGS_getLLt)then  ! - test case: unfactored Hessian in upper triangle
  hdx=fmatmul_mv(m=Lqhess,v=dx,typeMV="SUV")
else                    ! - simple but effective: get dx(t).L.L(t).dx without
  call putDiag(Lqhess,Ld) ! constructing L.L(t), saving O(3) matmul flops
  u  =fmatmul_mv(m=Lqhess,v=dx,typeMV="LtV")
  Hdx=fmatmul_mv(m=Lqhess,v=u, typeMV="LV")
endif
if(facBFGS_getLLt)then
  fad=dot_product(dx,hdx)
else
  fad=dot_product(u,u)
endif
if(all(abs(dg-Hdx)<=tol*max(abs(grad),abs(gradold))))then   ! prevents noisy updates
  err=0;message="bfgsFac_update/skipCond2(noisyUpdate)"
  return
endif
if(skipClassic)then ! - classical (Dennis and Schnabel) skipping conditions
  if(fac<=sqrt(epsRe)*norm2(dx/xscale)*norm2(dg*xscale))then  ! prevents indefinite updates
    err=0;message="bfgsFac_update/skipCond1(posDefUpdate)(classic)"
    return
  endif
elseif(dampedBFGS.and.fac<dampFac*fad)then  ! - BFGS damping
  dampTheta=(one-dampFac)*fad/(fad-fac)
  dg=dampTheta*dg+(one-dampTheta)*Hdx       ! r in Nocedal eq,(18.23)
  fac=dot_product(dg,dx)
elseif(.not.dampedBFGS)then  ! - minimal skipping condition required for stability
  if(abs(fac)<=sqrt(epsRe)*norm2(dx/xscale)*norm2(dg*xscale))then
    err=0;message="bfgsFac_update/skipCond1(minimal)"
    return
  endif
endif
if(facBFGS_useR2)then ! * do single rank-2 QR update
  alpha=sqrt(fac/fad); beta=sqrt(fac*fad)
  dg=dg-alpha*Hdx  ! t in DS96
  u=u/beta            ! u in DS96
  call choles_update(L=Lqhess,Ld=Ld,z=u,v=dg,keepLup=.false.,&
    getLLt=facBFGS_getLLt,logdet=logdet,condest=condest,err=err,message=jmsg)
  if(err/=0)then
    message="f-bfgsFac_update/BUG?/r2/&"//jmsg
    return
  endif
else                  ! * do two consecutive updates
  fac=one/fac; fad=one/fad
! - first rank-1 update, update preserves positive definiteness
  call choles_update(L=Lqhess,Ld=Ld,assumeLL=.true.,z=dg,alpha=fac,forceC2=forceC2inA,&
    logdet=logdet,condest=condest,err=err,message=jmsg)
  if(err/=0)then
    message="f-bfgsFac_update/BUG?/r1a/&"//jmsg
    return
  endif
! - second rank-1 update, needs to be robust to preserve +ve definiteness
  call choles_update(L=Lqhess,Ld=Ld,assumeLL=.true.,z=hdx,alpha=-fad,forceC2=.true.,&
    logdet=logdet,condest=condest,getLLt=facBFGS_getLLt,err=err,message=jmsg)
  if(err/=0)then
    message="f-bfgsFac_update/BUG?/r1b/&"//jmsg
    return
  endif
endif
if(controlHessCond.and.maxHessCond>zero)then  ! (crudely) maintain reasonable conditioning
  Ldmin=triang_minEig(Ld=Ld,condMax=maxHessCond,cholLd=.true.)
  where(Ld<Ldmin)Ld=Ldmin     ! by bumping up small diagonal elements
  condest=(maxval(Ld)/minval(Ld))**2
endif
if(.not.facBFGS_getLLt)then ! insert Cholesky diagonal to correctly compute
  call putDiag(Lqhess,Ld)   ! quadratic model reduction from Cholesky factors.
endif
! End procedure here
endsubroutine bfgsFac_update
!----------------------------------------------------
pure subroutine SR1unfac_update(dx,xscale,activeSet,grad,gradold,tol,qhess,rescale,err,message)
! Purpose: SR1 update of unfactored quasi-Hessian. Although the update is
! not guaranteed to be positive definite, it often becomes a better Hessian
! approximation than BFGS and the Broyden family updates. In these circumstances
! a trust region globalisation can become very effective.
! * Skipping conditions implemented as described by Nocedal,
!   to ensure the SR1 quasi-Hessian captures as much information as possible
! * Option available to rescale the initial diagonal Hessian after first
!   iteration but before first update using eqn (8.20) in Nocedal.
!   This can improve the scaling of Hessian for subsequent updates.
! * Routine works with upper triangle of Hessian only
use utilities_dmsl_kit,only:zero,one,norm2,fmatmul_mv,rank1updt
implicit none
! dummies
real(mrk),intent(in)::dx(:),grad(:),gradold(:),tol,xscale(:)
integer(mik),intent(in),optional::activeSet(:)
real(mrk),intent(inout)::qhess(:,:)
logical(mlk),intent(in)::rescale
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
real(mrk)::dg(size(dx)),vec(size(dx)),fac
! Start procedure here
err=0;message="SR1unfac_update/ok"; dg=grad-gradold
if(present(activeSet))then  ! need to zero dg for fixed variables
  where(activeSet/=freeVar_as)dg=zero
endif
vec=dg-fmatmul_mv(m=qhess,v=dx,typeMV="SUV")
fac=dot_product(vec,dx)
if(rescale)then ! rescale Hessian during first iteration
  call qhessRescale(qhess=qhess,invrs=.false.,dg=dg,fac=fac)
endif
if(abs(fac)<=tol*norm2(dx/xscale)*norm2(vec*xscale))then ! skip update
  err=0;message="SR1unfac_update/skipCond(NocedalWright)"
  return
endif
fac=one/fac
call rank1updt(a=qhess,fac=fac,x=vec,symm="U")
! End procedure here
endsubroutine SR1unfac_update
!----------------------------------------------------
pure subroutine qhessRescale(qhess,Ld,invrs,dg,fac)
! Purpose: Rescales the quasi-Newton Hessian diagonal after
! 1st iteration to improve its eigenvalue distribution.
! The update is a weak version of the secant requirement.
! This procedure handles unfactored, factored and inverse Hessians
! Ref: Nocedal and Wright, Dennis and Schnabel.
! Comments:
! 1. Works with diagonal of Hessian only;
! 2. I am not entirely sure whether fac is legally allowed to
!    become negative. For some cases I have seen -ve facs, so either
!    a bug, or legitimate occurence or some numerical precision issues.
use utilities_dmsl_kit,only:zero,putdiag
implicit none
! dummies
real(mrk),intent(out)::qhess(:,:)     ! Hessian (assumed diagonal)
real(mrk),intent(out),optional::Ld(:) ! Choelsky diagonal
logical(mlk),intent(in)::invrs        ! requests inverse update
real(mrk),intent(in)::dg(:),fac       ! required quantities
! locals
real(mrk)::dgdg
integer(mik)::i
! Start procedure here
dgdg=dot_product(dg,dg)
if(invrs)then ! inverse
  if(dgdg/=zero)call putdiag(qhess,fac/dgdg)  ! inverse Hessian diagonal
else
  if(fac/=zero) call putdiag(qhess,dgdg/fac)  ! Hessian diagonal
endif
if(present(Ld))then ! take some special care to avoid -ve sqrts
  if(fac<zero)forall(i=1:size(Ld))qhess(i,i)=abs(qhess(i,i))
  forall(i=1:size(Ld))Ld(i)=sqrt(qhess(i,i))  ! Cholesky diagonal
endif
! End procedure here
endsubroutine qhessRescale
!----------------------------------------------------
pure subroutine getXscaleH(xscaleHmeth,hess,xscale,fscale,xscaleH)
! Purpose: constructs Hessian ellipticity scaling.
! Supported scaling options:
! * Unit spherical
! * Elliptical - based on user xscale
! * Elliptical - based on Hessian (transforming to correlation-type matrix)
! ----
! Algorithm:
! When constructing scaling based on Hessian diagonal, it is worthwhile
! to bound the scaling, since d2f/dx2=0 would cause overflow and other problems.
! Currently the check is enabled only for those rows/columns where d2fdx2 seems
! too small. The rationale is that these elements may be simply numerical noise
! from cancellative computations, etc.
! The safeguarding of the factors is such that a multiplicative range of about
! safe*epsRe is given around the "default" scaling of xscale/sqrt(fscale).
! This heuristic may need revision in the future, so that sufficient freedom
! is given to perform the scaling as intended in principle, but safeguarding
! division by zero, where the definition of "zero" should reflect the a priori scaling
! of the Hessian, given by xscale and fscale. The current combination of this
! information is such that it reduces the initial Hessian estimate of Dennis and
! Schnabel (Algorithm 9.4.3) to the unit matrix.
use utilities_dmsl_kit,only:zero,one,getDiag,minmax
implicit none
! dummies
integer(mik),intent(in)::xscaleHmeth
real(mrk),intent(in)::hess(:,:),xscale(:),fscale
real(mrk),intent(out)::xscaleH(:)
! locals
!logical(mlk)::constrain(size(xscale))
real(mrk)::xscaleHmin(size(xscale)),xscaleHmax(size(xscale))
real(mrk),parameter::safe=1.e2_mrk,xscaleHmin0=safe*epsRe,xscaleHmax0=one/xscaleHmin0
! Start procedure here
selectcase(xscaleHmeth)
case(xscaleH_sphere)  ! - assume spherical contours
  xscaleH=one
case(xscaleH_user)    ! - user-prescribed (elliptical) scaling
  xscaleH=xscale
case(xscaleH_hdiag)   ! - ellipticity based on Hessian diagonal
  xscaleH=sqrt(fscale)/xscale     ! construct a priori bounds for scaling
  xscaleHmin=xscaleHmin0*xscaleH  ! factors to ensure stability and avoid
  xscaleHmax=xscaleHmax0*xscaleH  ! division by zero during scaling
  xscaleH=getDiag(hess)
!  constrain(:)=xscaleH(:)<=xscaleHmin**2 ! flag insufficiently positive diagonals
  xscaleH=sqrt(abs(xscaleH))  ! Hessian diagonal gives measure of sensitivity
!  where(constrain)xscaleH=minmax(xscaleH,xscaleHmin,xscaleHmax) ! of function to x-args
  xscaleH=minmax(xscaleH,xscaleHmin,xscaleHmax) ! of function to x-args
  xscaleH=one/xscaleH
endselect
! End procedure here
endsubroutine getXscaleH
!----------------------------------------------------
pure subroutine xscaleNewt(xscale,hess0,hessS,L0,LS,Ld0,LdS,grad0,gradS,p0,pS)
! Purpose: scales Hessian, L-factors, gradient and step
! using a diagonal scaling matrix xscale.
! NB: note step scaling is inverse!
implicit none
! dummies
real(mrk),intent(in)::xscale(:)
real(mrk),intent(in),   optional::hess0(:,:),L0(:,:),Ld0(:),grad0(:),p0(:)
real(mrk),intent(inout),optional::hessS(:,:),LS(:,:),LdS(:),gradS(:),pS(:)
! locals
integer(mik)::i,j,n
! Start procedure here
n=size(xscale)
if(present(hessS).and.present(hess0))then ! upper triangle of Hessian
  forall(i=1:n,j=1:n,i<=j)hessS(i,j)=hess0(i,j)*xscale(i)*xscale(j)
elseif(present(hessS))then
  forall(i=1:n,j=1:n,i<=j)hessS(i,j)=hessS(i,j)*xscale(i)*xscale(j)
endif
if(present(LS).and.present(L0))then       ! lower triangular factor L of Hessian
  forall(i=1:n,j=1:n,i>j)LS(i,j)=L0(i,j)*xscale(i)
elseif(present(LS))then
  forall(i=1:n,j=1:n,i>j)LS(i,j)=LS(i,j)*xscale(i)
endif
if(present(LdS).and.present(Ld0))then     ! diagonal of L-factor
  LdS=Ld0*xscale
elseif(present(LdS))then
  LdS=LdS*xscale
endif
if(present(gradS).and.present(grad0))then ! gradient
  gradS=grad0*xscale
elseif(present(gradS))then
  gradS=gradS*xscale
endif
if(present(pS).and.present(p0))then       ! step
  pS=p0/xscale
elseif(present(pS))then
  pS=pS/xscale
endif
! End procedure here
endsubroutine xscaleNewt
!----------------------------------------------------
pure subroutine unXscaleNewt(xscale,hess0,hessS,L0,LS,Ld0,LdS,grad0,gradS,p0,pS)
! Purpose: unscales Hessian, L-factors, gradient and step
! using a diagonal scaling matrix xscale.
! NB: note step scaling is inverse!
implicit none
! dummies
real(mrk),intent(in)::xscale(:)
real(mrk),intent(inout),optional::hess0(:,:),L0(:,:),Ld0(:),grad0(:),p0(:)
real(mrk),intent(in),   optional::hessS(:,:),LS(:,:),LdS(:),gradS(:),pS(:)
! locals
integer(mik)::i,j,n
! Start procedure here
n=size(xscale)
if(present(hess0).and.present(hessS))then ! upper triangle of Hessian
  forall(i=1:n,j=1:n,i<=j)hess0(i,j)=hessS(i,j)/xscale(i)/xscale(j)
elseif(present(hess0))then
  forall(i=1:n,j=1:n,i<=j)hess0(i,j)=hess0(i,j)/xscale(i)/xscale(j)
endif
if(present(L0).and.present(LS))then       ! lower triangular factor L of Hessian
  forall(i=1:n,j=1:n,i>j)L0(i,j)=LS(i,j)/xscale(i)
elseif(present(L0))then
  forall(i=1:n,j=1:n,i>j)L0(i,j)=L0(i,j)/xscale(i)
endif
if(present(Ld0).and.present(LdS))then     ! diagonal of L-factor
  Ld0=LdS/xscale
elseif(present(Ld0))then
  Ld0=Ld0/xscale
endif
if(present(grad0).and.present(gradS))then ! gradient
  grad0=gradS/xscale
elseif(present(grad0))then
  grad0=grad0/xscale
endif
if(present(p0).and.present(pS))then       ! step
  p0=pS*xscale
elseif(present(p0))then
  p0=p0*xscale
endif
! End procedure here
endsubroutine unXscaleNewt
!----------------------------------------------------
subroutine solveModNewtHess(hess,hessScaled,Ld,grad,hessFacBundle,&
  xscaleHmeth,xscale,fscale,activeset,dx,ncholstats,logdet,condest,Einf,err,message)
! Purpose: Processes the model Hessian equations for Newton-type optimisation
! using a modified factorization guaranteed to produce a positive definite
! matrix and hence descent direction.
! INPUT:
! hess        = full raw (unscaled) Hessian (may be indefinite, singular, etc.)
! hessScaled  = work array for Hessian decomposition
! grad        = full raw (unscaled) gradient
! hessFacBundle     = modified factorization bundle (settings etc.)
! xscaleHmeth = Hessian scaling method
! xscale      = user-provided xscale
! activeSet   = active set
! OUTPUT
! dx          = full solution of modified Hessian equations (Newton step)
! logdet      = log-determinant of modified Hessian
! condest     = condition estimate of modified Hessian
! Einf        = estimated most negative eigenvalue of input Hessian
! err         = error status
! message     = description of problems.
! Currently implemented factorization methods
! - modified Cholesky-Gershgorin w/wo pivoting, which perturb the Hessian diagonal
!   to achieve +ve definiteness and improve conditioning.
!---------
! Algorithm flowchart:
! Input:  Raw Hessian and gradient for Newton step
! Output: Full modified Newton solution
!
!         Raw Hessian (may be indefinite) ->
!         Active Hessian (excludes constrained variables) ->
!         Scaled active Hessian (accounting for diagonal scaling of vars) ->
!         Pivoted modified scaled active Hessian for Cholesky solution ->
!         Pivoted scaled active solution (to the modified problem) ->
!         Scaled active solution ->
!         Active solution ->
!         Full solution
!---------
! Comments
! * The work array hessScaled greatly simplifies memory management and reduces
!   arithmetic load in constructing the active Hessian, scaling and permuting it.
!   If this extra array is memory-busting, you should not be using dense Newton
!   in the first place - try conjugate gradient or truncated / limited memory Newton.
! * If the Cholesky algebra here is too much (but memory OK), can use factored
!   BFGS approximations which do not require explicit factorizations.
use utilities_dmsl_kit,only:zero,one,arthsi,terminateRowColMat
use linalg_dmsl_kit,only:choles_dcmp,choles_fwbw
implicit none
! dummies
real(mrk),intent(in)::hess(:,:),grad(:),xscale(:),fscale
real(mrk),intent(inout)::hessScaled(:,:) ! scratch Hessian
real(mrk),intent(out)::Ld(:)
type(hessFacBundle_type),intent(in)::hessFacBundle
integer(mik),intent(in)::xscaleHmeth
integer(mik),intent(in),optional::activeset(:)
real(mrk),intent(out)::dx(:)
integer(mik),intent(out)::ncholstats(:)
real(mrk),intent(out)::logdet,condest,Einf
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals for cholesky
integer(mik)::ndim,nchol
logical(mlk)::ok
!--locals for bounds
integer(mik)::lerr,nact
logical(mlk)::active(size(dx))
real(mrk)::xscaleH(size(dx))
integer(mik)::avar(size(dx))
character(100)::lmessage
! Cholesky pivoting
logical(mlk),parameter::doPivot=.true.
integer(mik)::indx(size(dx))
real(mrk)::gradScaled(size(dx))
! Start procedure here
ndim=size(dx);ncholstats=0
err=0;message="solveModNewtHess/ok"
if(present(activeset))then   ! * bound-contrained optimisation
  active=(activeSet==freeVar_as)
  nact=count(active)
else
  nact=ndim
endif
call getXscaleH(xscaleHmeth,hess,xscale,fscale,xscaleH)
if(nact==0)then         ! * all variables on bounds
  dx=zero; ok=.true.; logdet=zero; condest=zero; Einf=zero
  err=okAlg; message="w-solveModNewtHess/allVarsFixed"
else
  if(nact==ndim)then    ! * effectively unconstrained
    avar=arthsi(ndim)
    call xscaleNewt(xscaleH,hess0=hess,hessS=hessScaled)
  else                  ! * active constraints present: need to muck around
    avar=pack(arthsi(ndim),active) ! index of active variables
! pack active hessian, removing rows/columns corresponding to fixed variables
    call terminateRowColMat(hess,hessScaled(1:nact,1:nact),active,lerr,lmessage)
    call xscaleNewt(xscaleH(avar(1:nact)),hessS=hessScaled(1:nact,1:nact))
  endif ! also scale gradient
  call xscaleNewt(xscaleH(avar(1:nact)),grad0=grad(avar(1:nact)),gradS=gradScaled(1:nact))
  selectcase(hessFacBundle%facmeth)   ! select modified Hessian factorization method
  case(schnab_facmeth)  ! - revised modified Cholesky-Gershgorin of Schnabel and Eskew
    indx=-1             !  (indicate no pre-scrambling)
    call choles_dcmp(A=hessScaled(1:nact,1:nact),Ld=Ld(1:nact),&
      tau=hessFacBundle%tau,tauBar=hessFacBundle%tauBar,mu=hessFacBundle%mu,&
      doPivot=doPivot,indx=indx(1:nact),&
      posDefinite=ok,logDet=logdet,condest=condEst,Einf=Einf,err=err,message=lmessage)
    nchol=1 ! note if pivoting used then everything will be scrambled (and scaled)
  case(dennis_facmeth)  ! - perturbed Cholesky-Gershgorin of Dennis and Schnabel
    call choles_dcmp(A=hessScaled(1:nact,1:nact),Ld=Ld(1:nact),&
      maxCond=hessFacBundle%maxHessCond,&
      posDefinite=ok,nchol=nchol,logDet=logdet,condest=condEst,&
      Einf=Einf,err=err,message=lmessage)
  endselect
  ncholstats(1)=ncholstats(1)+nchol ! number of O(3) Cholesky factorizations
  ncholstats(2)=ncholstats(2)+1     ! number of internal iterations
  if(err/=okAlg)then                ! (usually 1 for linesearch, >1 for trusts)
    err=-20;message="f-solveModNewtHess/&"//lmessage
    return
  endif ! solve the scaled permuted Newton equations. the solution is unscrambled ...
  call choles_fwbw(a=hessScaled(1:nact,1:nact),Ld=Ld(1:nact),indx=indx(1:nact),&
                   usePivot=hessFacBundle%facmeth==schnab_facmeth.and.doPivot,&
                   b=gradScaled(1:nact),x=dx(1:nact),err=err,message=lmessage)
  call unXscaleNewt(xscaleH(avar(1:nact)),p0=dx(1:nact))  ! ... and now unscaled
  if(err/=okAlg)then
    err=-30;message="f-solveModNewtHess/&"//lmessage
    return
  endif
  if(nact<ndim)dx=unpack(dx(1:nact),active,zero) ! and expanded with fixed variables
  dx=-dx
endif
! End procedure here
endsubroutine solveModNewtHess
!----------------------------------------------------
subroutine linesearch_armijo(evalFunc,dataIN,dataOUT,x0,fx0,grad0,sdir,xscale,stol,alpha,&
                             stepmax,x,fx,fredAct,lambda,fcalls,retcode,message)
! Purpose: Linesearch using the Armijo condition (backtracking only).
! sdir is local search direction (typically Newton-derived).
! Note this linesearch does not impose the Wolfe conditions and hence
! may degrade the performance of quasi-Newton methods.
! Algorithm A6.3.1 in Dennis and Schnabel.
use types_dmsl_kit,only:data_ricz_type
use utilities_dmsl_kit,only:half,zero,one,three
implicit none
! dummies
type(data_ricz_type),intent(in),optional::dataIN
type(data_ricz_type),intent(inout),optional::dataOUT
real(mrk),intent(in)::x0(:),fx0,grad0(:),xscale(:),stol,alpha,stepmax
real(mrk),intent(inout)::lambda,sdir(:)
real(mrk),intent(out)::x(:),fx,fredAct
integer(mik),intent(out)::fcalls,retcode
character(*),intent(out)::message
! user-provided function
interface
  subroutine evalFunc(dataIN,dataOUT,x,feas,fx,gradFx,hessFx,err,message)
  use kinds_dmsl_kit
  use types_dmsl_kit,only:data_ricz_type
  implicit none
  type(data_ricz_type),intent(in),optional::dataIN
  type(data_ricz_type),intent(inout),optional::dataOUT
  real(mrk),intent(in)::x(:)
  logical(mlk),intent(out)::feas
  real(mrk),intent(out),optional::fx,gradFx(:),hessFx(:,:)
  integer(mik),intent(out)::err
  character(*),intent(out)::message
  endsubroutine evalfunc
endinterface
! locals
integer(mik)::err
real(mrk)::slope0,relLen,lambdaMin,lambdaTemp,lambdaPrev,fprev,stepLen
real(mrk)::m(2,2),c(2),disc
logical(mlk)::feas,firstRed
real(mrk),parameter::lambdaRedMax=0.1_mrk,lambdaRedMin=half
! Start procedure here
message="linesearch_armijo/ok"; fcalls=0; stepLen=getStepLen2(sdir,xscale)
if(stepLen>stepmax)then ! scale down large steps
  sdir=sdir*stepmax/stepLen; stepLen=stepmax
endif
slope0=dot_product(grad0,sdir) ! initial slope
if(slope0>zero)then       ! search direction is uphill
  x=x0; fx=fx0; retcode=badDir_glob; return
elseif(slope0==zero)then  ! no perceived slope in search direction
  x=x0; fx=fx0; fredAct=zero; retcode=failed_glob; return
endif
relLen=scaledStepLen(sdir,x0,xscale) ! scaled step length (used in termination test)
lambdaMin=stol/relLen ! minimum allowable steplength, lambda must be above noise
lambdaMin=max(minval(epsRe*max(abs(x0),xscale)/max(abs(sdir),xscale)),lambdaMin)
firstRed=.true.; retcode=failed_glob
do  ! loop to compute lambda that satisfies Armijo condition
  x=x0+lambda*sdir
  call evalFunc(dataIN,dataOUT,x,feas,fx,err=err,message=message)
  fcalls=fcalls+1
  if(err/=0)then
    retcode=badFunc_glob; message="f-linesearch_armijo/userErr/&"//message; return
  elseif(.not.feas)then ! reduce lambda and try again
    lambda=lambda*lambdaRedMax
    if(lambda<lambdaMin)then
      retcode=unfeas_glob; x=x0; fx=fx0
      exit
    endif
    cycle
  endif
  if(fx<fx0+alpha*lambda*slope0)then ! Armijo sufficient decrease condition ok
    retcode=success_glob; exit
  elseif(lambda<lambdaMin)then  ! no satisfactory x distinct from x0 in direction sdir
    retcode=failed_glob         ! was found (this may indicate convergence)
    if(fx>fx0)then  ! return original point if it was the best found
      x=x0; fx=fx0
    endif
    exit
  endif
  if(firstRed)then  ! * quadratic interpolation on backtrack
    firstRed=.false.; lambdaTemp=-slope0*half/(fx-fx0-slope0)
  else              ! * cubic interpolation on subsequent backtrack
    c(1:2)=  (/ fx-   fx0-lambda    *slope0,& ! enjoy some Fortran side-tracks
                fprev-fx0-lambdaPrev*slope0 /)
    m(1,1:2)=(/         one/lambda**2,   -one/lambdaPrev**2 /)
    m(2,1:2)=(/ -lambdaPrev/lambda**2, lambda/lambdaPrev**2 /)
    c=matmul(m,c)/(lambda-lambdaPrev)  ! c1*L^3+c2*L^2+slope0*L+fx0
    if(c(1)==zero)then  ! cubic is quadratic
      lambdaTemp=-slope0*half/c(2)
    else                ! legitimate cubic
      disc=c(2)**2-three*c(1)*slope0  ! discriminant of cubic
      if(disc<zero)then ! use minimum reduction when D<0
        lambdaTemp=lambdaRedMin*lambda
      else              ! cubic interpolation for minimum
        lambdaTemp=(-c(2)+sqrt(disc))/(three*c(1))
      endif
    endif
  endif
  lambdaPrev=lambda; fprev=fx ! safeguards to avoid extreme changes in lambda
  if(lambdaTemp>lambdaRedMin*lambda)then  ! ensure lambda is sufficiently decreased
    lambdaTemp=lambda*lambdaRedMin        ! to avoid stagnation
  elseif(lambdaTemp<lambdaRedMax*lambda)then  ! yet ensure lambda does not decrease
    lambdaTemp=lambda*lambdaRedMax            ! too fast to avoid stepsize collapse
  endif
  lambda=lambdaTemp
enddo
sdir=x-x0  ! get shift vector (must equal sdir*lambda)
fredAct=fx0-fx      ! actual reduction in function value
! End procedure here
endsubroutine linesearch_armijo
!----------------------------------------------------
subroutine linesearch_wolfe(evalFunc,dataIN,dataOUT,x0,fx0,grad0,gmethBundle,objFuncBundle,sdir,&
            xscale,fscale,stol,alpha,beta,stepmax,x,fx,gradFx,&
            fredAct,lambda,fcalls,gcalls,retcode,message)
! Purpose: Linesearch using the Wolfe conditions A (alpha) and B (beta).
! Condition A (sufficient decrease condition)
! Condition B (gradient condition)
! sdir is local search direction (typically Newton-derived).
! Algorithm A6.3.1mod in Dennis and Schnabel.
use types_dmsl_kit,only:data_ricz_type
use utilities_dmsl_kit,only:half,zero,one,three,two
implicit none
! dummies
type(data_ricz_type),intent(in),optional::dataIN
type(data_ricz_type),intent(inout),optional::dataOUT
real(mrk),intent(in)::x0(:),fx0,grad0(:),xscale(:),fscale,stepmax
real(mrk),intent(in)::stol,alpha,beta
type(gmethBundle_type),intent(in)::gmethBundle
type(objFuncBundle_type),intent(in)::objFuncBundle
real(mrk),intent(inout)::lambda,sdir(:)
real(mrk),intent(out)::x(:),fx,gradFx(:),fredAct
integer(mik),intent(out)::fcalls,gcalls,retcode
character(*),intent(out)::message
! user-provided function
interface
  subroutine evalFunc(dataIN,dataOUT,x,feas,fx,gradFx,hessFx,err,message)
  use kinds_dmsl_kit
  use types_dmsl_kit,only:data_ricz_type
  implicit none
  type(data_ricz_type),intent(in),optional::dataIN
  type(data_ricz_type),intent(inout),optional::dataOUT
  real(mrk),intent(in)::x(:)
  logical(mlk),intent(out)::feas
  real(mrk),intent(out),optional::fx,gradFx(:),hessFx(:,:)
  integer(mik),intent(out)::err
  character(*),intent(out)::message
  endsubroutine evalfunc
endinterface
! locals
real(mrk)::slope0,relLen,lambdaMin,lambdaTemp,lambdaPrev,fprev,stepLen
real(mrk)::m(2,2),c(2),disc
logical(mlk)::feas,firstRed,firstSearch,finalGrad
real(mrk),parameter::lambdaRedMax=0.1_mrk,lambdaRedMin=half
real(mrk)::slopeX,lambdaMax,Llo,Ldiff,Lincr,flo,fhi
integer(mik)::addFcalls,err
! Start procedure here
err=0; stepLen=getStepLen2(sdir,xscale)
if(stepLen>stepmax)then ! scale down large steps
  sdir=sdir*stepmax/stepLen; stepLen=stepmax
endif
slope0=dot_product(grad0,sdir) ! initial slope
if(slope0>=zero)then
  retcode=badDir_glob
  x=x0;fx=fx0;gradFx=grad0
  return
endif
relLen=scaledStepLen(sdir,x0,xscale) ! scaled step length (used in termination test)
lambdaMin=stol/relLen ! minimum allowable steplength, lambda must be above noise
lambdaMin=max(minval(epsRe*max(abs(x0),xscale)/max(abs(sdir),xscale)),lambdaMin)
lambdaMax=stepMax/stepLen
firstRed=.true.; fcalls=0; gcalls=0; firstSearch=.true.; retcode=failed_glob
finalGrad=.false.
outer_loop: do  ! loop to compute lambda that satisfies Wolfe conditions
  x=x0+lambda*sdir
  call evalFunc(dataIN,dataOUT,x,feas,fx,err=err,message=message)
  fcalls=fcalls+1
  if(err/=0)then
    retcode=badFunc_glob; message="f-linesearch_wolfe/userErr1/&"//message; return
  elseif(.not.feas)then ! reduce lambda and try again
    lambda=lambda*lambdaRedMax
    if(lambda<lambdaMin)then
      retcode=unfeas_glob
      x=x0;fx=fx0;gradFx=grad0
      exit
    endif
    cycle
  endif
  if(fx<fx0+alpha*lambda*slope0)then ! condition A ok
! ---- start modification to enforce Wolfe condition B (gradient condition)
    call getGradSlopeX()
    if(err/=0)return
    if3a4:if(slopeX<beta*slope0)then  ! violation of Wolfe condition B
      if41: if(firstSearch)then ! 10.3a.4.1
        firstSearch=.false.
        do  ! 10.3a.4.1.2
          lambdaprev=lambda; fprev=fx
          lambda=min(two*lambda,lambdaMax)  ! increase Lambda
          x=x0+lambda*sdir
          call evalFunc(dataIN,dataOUT,x,feas,fx,err=err,message=message)
          fcalls=fcalls+1
          if(err/=0)then
            retcode=badFunc_glob; message="f-linesearch_wolfe/userErr2/&"//message; return
          elseif(.not.feas)then ! unfeasible point encountered
            retcode=unfeas_glob
            x=x0;fx=fx0;gradFx=grad0
            exit outer_loop
          endif
          if(fx<fx0+alpha*lambda*slope0)then ! condition A satisfied
            call getGradSlopeX()
            if(err/=0)return
          endif
          if(fx>=fx0+alpha*lambda*slope0.or. & ! cond A now violated
             slopeX>=beta*slope0.or.        &  ! cond B now satisfied
             lambda>=lambdaMax)             &  ! ran out of lambda's
            exit  ! time to pop-out 10.3a.4.1.2U
        enddo
      endif if41
      if42: if(lambda<one.or.(lambda>one.and.fx>=fx0+alpha*lambda*slope0))then
        Llo=min(lambda,lambdaPrev)
        Ldiff=abs(LambdaPrev-Lambda)
        if(lambda<lambdaPrev)then
          flo=fx; fhi=fprev
        else
          flo=fprev; fhi=fx
        endif
        do  ! 10.3a.4.2.4
          Lincr=(-slopeX*Ldiff**2)*half/(fhi-(flo+slopeX*Ldiff))
          if(Lincr<0.2_mrk*Ldiff)Lincr=0.2_mrk*Ldiff
          lambda=Llo+Lincr
          x=x0+lambda*sdir
          call evalFunc(dataIN,dataOUT,x,feas,fx,err=err,message=message)
          fcalls=fcalls+1
          if(err/=0)then
            retcode=badFunc_glob; message="f-linesearch_wolfe/userErr3/&"//message; return
          elseif(.not.feas)then ! unfeasible point encountered
            retcode=unfeas_glob
            x=x0;fx=fx0;gradFx=grad0
            exit outer_loop
          endif
          if(fx>=fx0+alpha*lambda*slope0)then  ! 10.3a.4.2.4.6
            Ldiff=Lincr; fhi=fx
          else
            call getGradSlopeX()
            if(err/=0)return
            if(slopeX<beta*slope0)then
              Llo=lambda; Ldiff=Ldiff-Lincr; flo=fx
            endif
          endif
          if(slopeX>=beta*slope0.or.Ldiff<LambdaMin)exit  ! 10.3a.4.2.4U
        enddo
        if(slopeX<beta*slope0)then  ! could not satisfy condition B
          fx=flo; x=x0+Llo*sdir
        endif
      endif if42
    endif if3a4
    retcode=success_glob
    exit outer_loop
! ---- end modification to enforce Wolfe condition B
  elseif(lambda<lambdaMin)then  ! no satisfactory x distinct from x0 in direction sdir
    retcode=failed_glob         ! was found (this may indicate convergence)
    if(fx>fx0)then  ! return original point if it was the best found
      x=x0; fx=fx0; gradFx=grad0
    endif
    exit outer_loop
  endif
  if(firstRed)then  ! * quadratic interpolation on backtrack
    firstRed=.false.
    lambdaTemp=-slope0*half/(fx-fx0-slope0)
  else              ! * cubic interpolation on subsequent backtrack
    c(1:2)=  (/ fx-   fx0-lambda    *slope0,& ! enjoy some Fortran side-tracks
                fprev-fx0-lambdaPrev*slope0 /)
    m(1,1:2)=(/         one/lambda**2,   -one/lambdaPrev**2 /)
    m(2,1:2)=(/ -lambdaPrev/lambda**2, lambda/lambdaPrev**2 /)
    c=matmul(m,c)/(lambda-lambdaPrev)  ! c1*L^3+c2*L^2+slope0*L+fx0
    if(c(1)==zero)then  ! cubic is quadratic
      lambdaTemp=-slope0*half/c(2)
    else                ! legitimate cubic
      disc=c(2)**2-three*c(1)*slope0  ! discriminant of cubic
      if(disc<zero)then ! use minimum reduction when D<0
        lambdaTemp=lambdaRedMin*lambda
      else              ! cubic interpolation for minimum
        lambdaTemp=(-c(2)+sqrt(disc))/(three*c(1))
      endif
    endif
    if(lambdaTemp>lambdaRedMin*lambda)&     ! ensure lambda is sufficiently decreased
          lambdaTemp=lambdaRedMin*lambda    ! to avoid stagnation
  endif
  lambdaPrev=lambda; fprev=fx   ! safeguards to avoid spurious changes in lambda
  if(lambdaTemp<lambdaRedMax*lambda)&       ! yet ensure lambda does not decrease
          lambdaTemp=lambdaRedMax*lambda    ! too fast to avoid stepsize collapse
  lambda=lambdaTemp
enddo outer_loop
sdir=x-x0 ! get shift vector (must equal sdir*lambda)
fredAct=fx0-fx     ! actual reduction in function value
if(gmethBundle%useDirDer)then ! final call to get full gradient
  finalGrad=.true.; call getGradSlopeX()
endif
! End procedure here
contains
!--
subroutine getGradSlopeX()  ! macro to compute directional derivative
use utilities_dmsl_kit,only:getHxFromRelHx,getFDCDgrad,getCDgrad,getFDdirDer
implicit none
! dummies
! local registered settings
integer(mik),parameter::scal_smeth=0,imax_smeth=1,ave_smeth=2,wei_smeth=3
! Start procedure here
selectcase(gmethBundle%gmeth_now)
case(user_meth)     ! * user-provided gradient
  call evalFunc(dataIN,dataOUT,x,feas=feas,gradFx=gradFx,err=err,message=message); gcalls=gcalls+1
  if(err/=0)then        ! strange, since point x has already been trialled with err=0
    message="f-linesearch_wolfe/getGradSlopeX/userBug/userErr/&"//message
    retcode=bugFail;  err=bugFail; return
  elseif(.not.feas)then ! strange, since point x already tested ok...
    message="f-linesearch_wolfe/getGradSlopeX/userBug/userUnfeas/&"//message
    retcode=bugFail;  err=bugFail; return
  endif
  slopeX=dot_product(gradFx,sdir) ! new slope
case(fd_gmeth)      ! * FD gradient
  if(.not.finalGrad.and.gmethBundle%useDirDer)then ! use cheap directional derivative
    call getFDdirDer(evalFunc,dataIN,dataOUT,x=x,p=sdir,fx=fx,xscale=xscale,fscale=fscale,&
      epsF=objFuncBundle%epsF,&
      hx=getHxFromRelHx(gmethBundle%hx,x,xscale,gmethBundle%FDscale),&
      useHxDef=gmethBundle%useHxDef,&
      dmeth=merge(useFDCDhybrid,fd_gmeth,gmethBundle%hybridFDCD),&
      dFDCD=gmethBundle%tolGradFDCD,&
      smeth=scal_smeth,normalize=.false.,&
      fdDirDer=slopeX,fcalls=addFcalls,err=err,message=message)
  else
    call getFDCDgrad(evalFunc,dataIN,dataOUT,x,fx,xscale,fscale,objFuncBundle%epsF,&
      getHxFromRelHx(gmethBundle%hx,x,xscale,gmethBundle%FDscale),&
      gmethBundle%useHxDef,&
      merge(useFDCDhybrid,fd_gmeth,gmethBundle%hybridFDCD),gmethBundle%tolGradFDCD,&
      gradFx,addFcalls,err,message)
    slopeX=dot_product(gradFx,sdir) ! new slope
  endif
  fcalls=fcalls+addFcalls
case(cd_gmeth)      ! * CD gradient
  if(.not.finalGrad.and.gmethBundle%useDirDer)then ! use cheap directional derivative
    call getFDdirDer(evalFunc,dataIN,dataOUT,x=x,p=sdir,fx=fx,xscale=xscale,fscale=fscale,&
      epsF=objFuncBundle%epsF,&
      hx=getHxFromRelHx(gmethBundle%hx,x,xscale,gmethBundle%FDscale),&
      useHxDef=gmethBundle%useHxDef,&
      dmeth=cd_gmeth,&
      dFDCD=gmethBundle%tolGradFDCD,&
      smeth=scal_smeth,normalize=.false.,&
      fdDirDer=slopeX,fcalls=addFcalls,err=err,message=message)
  else
    call getCDgrad(evalFunc,dataIN,dataOUT,x,fx,xscale,objFuncBundle%epsF,&
      getHxFromRelHx(gmethBundle%hx,x,xscale,gmethBundle%FDscale),&
      gmethBundle%useHxDef,gradFx,addFcalls,err,message)
    fcalls=fcalls+addFcalls
    slopeX=dot_product(gradFx,sdir) ! new slope
  endif
  fcalls=fcalls+addFcalls
endselect
if(err/=0)then  ! either a big or rather unfeasible problem
  retcode=unfeas_glob; message="f-linesearch_wolfe/getGradSlopeX/&"//message
endif
! End procedure here
endsubroutine getGradSlopeX
!--
endsubroutine linesearch_wolfe
!----------------------------------------------------
subroutine linesearch_strongwolfe(evalFunc,dataIN,dataOUT,x0,fx0,grad0,gmethBundle,objFuncBundle,&
            sdir,xscale,fscale,stol,alpha,beta,stepmax,strongwolfe,x,fx,gradFx,&
            fredAct,lambda,fcalls,gcalls,retcode,message)
! Purpose: Linesearch using the strong Wolfe conditions A (alpha) and B (beta).
! Condition A (sufficient decrease condition)
! Condition B (absolute gradient condition)
! sdir is local search direction (typically Newton-derived).
! Programmer: Dmitri Kavetski
use types_dmsl_kit,only:data_ricz_type
implicit none
! dummies
type(data_ricz_type),intent(in),optional::dataIN
type(data_ricz_type),intent(inout),optional::dataOUT
real(mrk),intent(in)::x0(:),fx0,grad0(:),xscale(:),fscale,stepmax,stol,alpha,beta
type(gmethBundle_type),intent(in)::gmethBundle
type(objFuncBundle_type),intent(in)::objFuncBundle
real(mrk),intent(inout)::lambda,sdir(:)
integer(mik),intent(in)::strongwolfe
real(mrk),intent(out)::x(:),fx,gradFx(:),fredAct
integer(mik),intent(out)::fcalls,gcalls,retcode
character(*),intent(out)::message
! user-provided function
interface
  subroutine evalFunc(dataIN,dataOUT,x,feas,fx,gradFx,hessFx,err,message)
  use kinds_dmsl_kit
  use types_dmsl_kit,only:data_ricz_type
  implicit none
  type(data_ricz_type),intent(in),optional::dataIN
  type(data_ricz_type),intent(inout),optional::dataOUT
  real(mrk),intent(in)::x(:)
  logical(mlk),intent(out)::feas
  real(mrk),intent(out),optional::fx,gradFx(:),hessFx(:,:)
  integer(mik),intent(out)::err
  character(*),intent(out)::message
  endsubroutine evalfunc
endinterface
! locals
! Start procedure here
selectcase(strongwolfe)
case(strongwolfe_more)      ! Fairly sophisticated Strong Wolfe linesearch
  call linesearch_more(evalFunc,dataIN,dataOUT,x0,fx0,grad0,gmethBundle,objFuncBundle,sdir,&
                       xscale,fscale,stol,alpha,beta,stepmax,x,fx,gradFx,&
                       fredAct,lambda,fcalls,gcalls,retcode,message)
case(strongwolfe_fletcher)  ! Brute force "bisection"-style beast
  call linesearch_fletcher(evalFunc,dataIN,dataOUT,x0,fx0,grad0,gmethBundle,objFuncBundle,sdir,&
                       xscale,fscale,stol,alpha,beta,stepmax,x,fx,gradFx,&
                       fredAct,lambda,fcalls,gcalls,retcode,message)
endselect
! End procedure here
endsubroutine linesearch_strongwolfe
!----------------------------------------------------
subroutine linesearch_more(evalFunc,dataIN,dataOUT,x0,fx0,grad0,gmethBundle,objFuncBundle,&
            sdir,xscale,fscale,stol,alpha,beta,stepmax,x,fx,gradFx,&
            fredAct,lambda,fcalls,gcalls,retcode,message)
! Purpose: Linesearch using the strong Wolfe conditions A (alpha) and B (beta).
! Condition A (sufficient decrease condition)
! Condition B (absolute gradient condition)
! sdir is local search direction (typically Newton-derived).
! Programmer: Dmitri Kavetski
! Ref: More, J.J. and Thuente, D.J. (1994) Line search algorithms with
!      guaranteed sufficient decrease, ACM Transactions on Mathematical
!      software, vol. 20(3), p.286-307.
! Method: clever concoction of safeguarded quadratic/cubic extrapolation.
! Note: by reducing beta the linesearch becomes a line-minimisation.
! More's algorithm is probably not the best univariate optimisation method,
! (instead is is geared more towards "loose" linesearches that appear
! best for quasi-Newton methods). If beta is small then the final convergence
! can be rather finicky due to finite preciison arithmetic. It is rather
! surprising how problematic is polynomial interpolation of closely spaces
! abscissas. You can always fall back on the brute force bisection-type linesearch
! using Fletcher's method (linesearch_fletcher).
use types_dmsl_kit,only:data_ricz_type
use utilities_dmsl_kit,only:zero,one,two,half,average,&
  quadFitStation,cubiqFitStation
implicit none
! dummies
type(data_ricz_type),intent(in),optional::dataIN
type(data_ricz_type),intent(inout),optional::dataOUT
real(mrk),intent(in)::x0(:),fx0,grad0(:),xscale(:),fscale,stepmax
real(mrk),intent(in)::stol,alpha,beta
type(gmethBundle_type),intent(in)::gmethBundle
type(objFuncBundle_type),intent(in)::objFuncBundle
real(mrk),intent(inout)::lambda,sdir(:)
real(mrk),intent(out)::x(:),fx,gradFx(:),fredAct
integer(mik),intent(out)::fcalls,gcalls,retcode
character(*),intent(out)::message
! user-provided function
interface
  subroutine evalFunc(dataIN,dataOUT,x,feas,fx,gradFx,hessFx,err,message)
  use kinds_dmsl_kit
  use types_dmsl_kit,only:data_ricz_type
  implicit none
  type(data_ricz_type),intent(in),optional::dataIN
  type(data_ricz_type),intent(inout),optional::dataOUT
  real(mrk),intent(in)::x(:)
  logical(mlk),intent(out)::feas
  real(mrk),intent(out),optional::fx,gradFx(:),hessFx(:,:)
  integer(mik),intent(out)::err
  character(*),intent(out)::message
  endsubroutine evalfunc
endinterface
! local parameters
real(mrk),parameter::lambdaRedMax=0.1_mrk,lambdaRedUnfeas=0.2_mrk
real(mrk),parameter::lambdaUpFac=two
real(mrk),parameter::delSafe=0.66_mrk,safeEps=10._mrk
integer(mik),parameter::itsearchmax=100
! locals
integer(mik)::itsearch
real(mrk)::slope0,slopeX,stepLen,relLen,lambdaMin,lambdaMax,lambdaTemp,temp
real(mrk)::Llo,Lhi,Flo,Fhi,Glo,Ghi,Lcube,Lquad,Lquad2,LcubeMax
integer(mik)::typeQ1,typeQ2,typeC1,typeC2,err,addFcalls
logical(mlk)::feas,haveHi,useMod,finalGrad
! Start procedure here
err=0; message="linesearch_more/ok"; stepLen=getStepLen2(sdir,xscale)
if(stepLen>stepmax)then ! scale down large steps
  sdir=sdir*stepmax/stepLen; stepLen=stepmax
endif
slope0=dot_product(grad0,sdir) ! initial slope
if(slope0>=zero)then
  retcode=badDir_glob;message="f-linesearch_more/badDir"
  x=x0;fx=fx0;gradFx=grad0
  return
endif
relLen=scaledStepLen(sdir,x0,xscale) ! relative steplength (used in termination test)
lambdaMin=stol/relLen ! minimum allowable steplength, lambda must be above noise
lambdaMin=max(minval(epsRe*max(abs(x0),xscale)/max(abs(sdir),xscale)),lambdaMin)
lambdaMax=stepMax/stepLen         ! maximum allowable steplength
Llo=zero; Flo=fx0; Glo=slope0
Flo=Flo-fx0-alpha*slope0*Llo; Glo=Glo-alpha*slope0
Lhi=lambdaMax; haveHi=.false.; finalGrad=.false.
fcalls=0; gcalls=0; retcode=failed_glob; useMod=.false.
finalGrad=.false.
do itsearch=1, itsearchmax
! * evaluate trial lambda
  x=x0+lambda*sdir
  call getGradSlopeX()  ! evaluate function and directional derivative
  if(err/=0)then
    return
  elseif(.not.feas)then ! reduce lambda and try again
    Lhi=half*lambda; haveHi=.false.; lambdaMax=Lhi
    lambda=max(lambda*lambdaRedMax,Llo*(one-lambdaRedUnfeas)+lambdaRedUnfeas*Lhi)
    if(lambda<=lambdaMin)then  ! nothing feasible in requested direction
      if((useMod.and.fx>fx0).or.&
         (.not.useMod.and.fx+fx0+alpha*slope0*lambda>fx0))then
        x=x0;fx=fx0;gradFx=grad0
      endif
      retcode=unfeas_glob; exit
    endif
    cycle
  endif
! * check Wolfe conditions
  if(fx<fx0+alpha*lambda*slope0.and.abs(slopeX)<=-beta*slope0)then
    retcode=success_glob; exit  ! Wolfe conditions OK
  elseif(abs(Llo-Lhi)<=safeEps*epsRe*max(Llo,Lhi))then
    retcode=success_glob; exit  ! bracket collapsed - mite as well quit...
  elseif(average(n1=Llo,n2=Lhi)<lambdaMin)then
    retcode=failed_glob; x=x0; fx=fx0; gradFx=grad0
    exit  ! bracket collapsed to initial point ...
  endif
! * check for switch to modified updating algorithm
  if(.not.useMod)then
    temp=fx-fx0-alpha*slope0*lambda
    if(temp<=zero.and.slopeX>=zero)then
      useMod=.true. ! use modified algorithm from now on
      Flo=Flo+fx0+alpha*slope0*Llo; Glo=Glo+alpha*slope0 ! recover lower bracket
    else
      fx=temp; slopeX=slopeX-alpha*slope0
    endif
  endif
! * generate safeguarded trial value
!  write(*,*)"fx=",fx,"flo=",flo,"glo=",glo,"slopeX=",slopeX
  if(fx>Flo)then                            ! ** case 1 (interpolation)
    call quadFitStation(xa=Llo,xb=lambda,fA=Flo,fB=fx,dfA=Glo,&
                        xs=Lquad,ts=typeQ1,err=err,message=message)
    if(err/=0)then
      typeQ1=10;err=0
!      retcode=bugFail;message="f-linesearch_more/bug?/A/&"//message
!      return
    endif
    call cubiqFitStation(xa=Llo,xb=lambda,fA=Flo,fB=fx,dfA=Glo,dfB=slopeX,&
                         xs1=Lcube,xs2=LcubeMax,ts1=typeC1,ts2=typeC2,&
                         err=err,message=message)
    if(err/=0)then
      typeC1=10;err=0
!      retcode=bugFail;message="f-linesearch_more/bug?/B/&"//message
!      return
    endif
    if((typeQ1==-1.and.typeC1==-1.and.abs(Lcube-Llo)<abs(Lquad-Llo)).or.&
       (typeQ1/=-1.and.typeC1==-1))then
      LambdaTemp=Lcube
    elseif(typeQ1==-1.and.typeC1==-1)then
      LambdaTemp=half*(Lquad+Lcube)
    else
      LambdaTemp=average(n1=Llo,n2=Lhi) ! safeguard degeneration
    endif
  elseif(fx<=Flo.and.slopeX*Glo<zero)then   ! ** case 2 (interpolation)
    call quadFitStation(xa=Llo,xb=lambda,fA=Flo,dfA=Glo,dfB=slopeX,&
                        xs=Lquad2,ts=typeQ2,err=err,message=message)
    if(err/=0)then
      typeQ2=10;err=0
!      retcode=bugFail;message="f-linesearch_more/bug?/C/&"//message
!      return
    endif
    call cubiqFitStation(xa=Llo,xb=lambda,fA=Flo,fB=fx,dfA=Glo,dfB=slopeX,&
                         xs1=Lcube,xs2=LcubeMax,ts1=typeC1,ts2=typeC2,&
                         err=err,message=message)
    if(err/=0)then
      typeC1=10;err=0
!      retcode=bugFail;message="f-linesearch_more/bug?/D/&"//message
!      return
    endif
    if((typeQ2==-1.and.typeC1==-1.and.abs(Lcube-lambda)>=abs(Lquad2-lambda)).or.&
       (typeQ2/=-1.and.typeC1==-1))then
      LambdaTemp=Lcube
    elseif(typeQ2==-1)then
      LambdaTemp=Lquad2
    else
      LambdaTemp=average(n1=Llo,n2=Lhi)
    endif
  elseif(fx<=Flo.and.slopeX*Glo>=zero.and.abs(slopeX)<=abs(Glo))then
                                            ! ** case 3 (extrapolation)
    call quadFitStation(xa=Llo,xb=lambda,fA=Flo,dfB=slopeX,dfA=Glo,&
                        xs=Lquad2,ts=typeQ2,err=err,message=message)
    if(err/=0)then
      typeQ2=10;err=0
!      retcode=bugFail;message="f-linesearch_more/bug?/E/&"//message
!      return
    endif
    call cubiqFitStation(xa=Llo,xb=lambda,fA=Flo,fB=fx,dfA=Glo,dfB=slopeX,&
                         xs1=Lcube,xs2=LcubeMax,ts1=typeC1,ts2=typeC2,&
                         err=err,message=message)
! - in this case err/=0 often indicates convergence and roundoff error
    if(err/=0)then
      typeC1=10;err=0
!      retcode=bugFail;message="f-linesearch_more/bug?/F/&"//message
!      return
    endif
    if(typeC1==-1.and.Lcube>lambda)then
      if(abs(Lcube-lambda)<abs(Lquad2-lambda).or.typeQ2/=-1)then  ! cubic step valid
        lambdaTemp=Lcube
      elseif(typeQ2==-1)then  ! secant step valid
        lambdaTemp=Lquad2
      else
        lambdaTemp=min(lambda*lambdaUpFac,max(Lhi,Llo))
      endif
    elseif(typeQ2==-1)then  ! secant step valid
      lambdaTemp=Lquad2
    else                    ! no valid steps: signal apparent convergence
      lambdaTemp=lambda
    endif
    if(abs(LambdaTemp-lambda)<safeEps*epsRe*(LambdaTemp+lambda))then
! Extrapolation converged to machine precision
      if(.not.useMod)then  ! recover function value
        fx=fx+fx0+alpha*slope0*Llo
      endif
      retcode=success_glob
      exit
    elseif(LambdaTemp>=lambda.and.Lambda>=lambdaMax)then ! maximum step taken
      if(.not.useMod)then  ! recover function value
        fx=fx+fx0+alpha*slope0*Llo
      endif
      retcode=success_glob
      exit
    elseif(LambdaTemp>lambdaMax.and.Lambda<lambdaMax)then ! stick to maximum step
      lambdaTemp=lambdaMax
    endif
    if(lambdaTemp>Llo)then  ! employ safeguards
      LambdaTemp=min(lambda+delSafe*(Lhi-lambda),lambdaTemp)
    else
      LambdaTemp=max(lambda+delSafe*(Lhi-lambda),lambdaTemp)
    endif
  else                                      ! ** case 4 (interpolation)
    if(haveHi)then
      call cubiqFitStation(xa=Lhi,xb=lambda,fA=Fhi,fB=fx,dfA=Ghi,dfB=slopeX,&
                           xs1=Lcube,xs2=LcubeMax,ts1=typeC1,ts2=typeC2,&
                           err=err,message=message)
      if(err/=0)then
        typeC1=10;err=0
!        retcode=bugFail;message="f-linesearch_more/bug?/G/&"//message
!        return
      endif
      if(typeC1==-1)then  ! cubic minimum "safely" computed
        LambdaTemp=Lcube
      else
        LambdaTemp=average(n1=Llo,n2=Lhi) ! safeguard degeneration
      endif
    else
      LambdaTemp=min(lambda*lambdaUpFac,max(Lhi,Llo)) ! Lhi and Llo may not be ordered
    endif
  endif
! * update brackets using "modified updating algorithm", page 297
  if(fx>Flo)then
    Lhi=lambda; Fhi=fx; Ghi=slopeX; haveHi=.true.
  else
    if(slopeX*(Llo-lambda)<zero)then
      Lhi=Llo; Fhi=Flo; Ghi=Glo; haveHi=.true.
    endif
    Llo=lambda; Flo=fx; Glo=slopeX
  endif
  Lambda=LambdaTemp
!  if(Lambda<Llo.or.Lambda>Lhi)then ! safeguard against bugs just in case.
!    retcode=bugFail;message="f-linesearch_more/bug?/lambdaOutOfBrackets"
!    return
!    Lambda=average(n1=Llo,n2=Lhi)
!  endif
enddo
if(itsearch>=itsearchmax+1)then
  x=x0; fx=fx0; gradFx=grad0; sdir=zero; fredAct=zero
  retcode=failed_glob
else    ! * some additional postcalculations
  sdir=x-x0 ! get shift vector (must equal sdir*lambda)
  fredAct=fx0-fx     ! actual reduction in function value
  if(gmethBundle%useDirDer)then ! final call to get full gradient
    finalGrad=.true.; call getGradSlopeX()
  endif
endif
! End procedure here
contains
!--
subroutine getGradSlopeX()  ! macro to get directional derivative
use utilities_dmsl_kit,only:getFDCDgrad,getCDgrad,getHxFromRelHx,getFDdirDer
implicit none
! dummies
! locals
! local registered settings
integer(mik),parameter::scal_smeth=0,imax_smeth=1,ave_smeth=2,wei_smeth=3
! Start procedure here
selectcase(gmethBundle%gmeth_now)
case(user_meth)   ! analytical derivatives available
  call evalFunc(dataIN,dataOUT,x,feas,fx,gradFx,err=err,message=message);fcalls=fcalls+1;gcalls=gcalls+1
  if(err/=0)then
    message="f-linesearch_wolfe/getGradSlopeX/userErrA/&"//message
    retcode=badFunc_glob; return
  elseif(.not.feas)then
    message="f-linesearch_wolfe/getGradSlopeX/userUnfeasA/&"//message; return
  endif
  slopeX=dot_product(gradFx,sdir) ! new slope
case(fd_gmeth)    ! forward difference gradient
  if(.not.finalGrad)then
    call evalFunc(dataIN,dataOUT,x,feas,fx,err=err,message=message); fcalls=fcalls+1
    if(err/=0)then
      retcode=badFunc_glob; message="f-linesearch_more/getGradSlopeX/userErrFD/&"//message; return
    elseif(.not.feas)then  ! do not bother with gradient if unfeasible
      message="f-linesearch_wolfe/getGradSlopeX/userUnfeasFD/&"//message; return
    endif
  endif
  if(.not.finalGrad.and.gmethBundle%useDirDer)then ! use cheap directional derivative
    call getFDdirDer(evalFunc,dataIN,dataOUT,x=x,p=sdir,fx=fx,xscale=xscale,fscale=fscale,&
      epsF=objFuncBundle%epsF,&
      hx=getHxFromRelHx(gmethBundle%hx,x,xscale,gmethBundle%FDscale),&
      useHxDef=gmethBundle%useHxDef,&
      dmeth=merge(useFDCDhybrid,fd_gmeth,gmethBundle%hybridFDCD),&
      dFDCD=gmethBundle%tolGradFDCD,&
      smeth=scal_smeth,normalize=.false.,&
      fdDirDer=slopeX,fcalls=addFcalls,err=err,message=message)
  else
    call getFDCDgrad(evalFunc,dataIN,dataOUT,x,fx,xscale,fscale,objFuncBundle%epsF,&
      getHxFromRelHx(gmethBundle%hx,x,xscale,gmethBundle%FDscale),&
      gmethBundle%useHxDef,&
      merge(useFDCDhybrid,fd_gmeth,gmethBundle%hybridFDCD),gmethBundle%tolGradFDCD,&
      gradFx,addFcalls,err,message)
    slopeX=dot_product(gradFx,sdir)
  endif
  fcalls=fcalls+addFcalls
case(cd_gmeth)    ! central difference gradient
  if(.not.finalGrad)then
    call evalFunc(dataIN,dataOUT,x,feas,fx,err=err,message=message); fcalls=fcalls+1
    if(err/=0)then
      retcode=badFunc_glob; message="f-linesearch_more/getGradSlopeX/userErrCD/&"//message; return
    elseif(.not.feas)then  ! do not bother with gradient if unfeasible
      message="f-linesearch_wolfe/getGradSlopeX/userUnfeasCD/&"//message; return
    endif
  endif
  if(.not.finalGrad.and.gmethBundle%useDirDer)then ! use cheap directional derivative
    call getFDdirDer(evalFunc,dataIN,dataOUT,x=x,p=sdir,fx=fx,xscale=xscale,fscale=fscale,&
      epsF=objFuncBundle%epsF,&
      hx=getHxFromRelHx(gmethBundle%hx,x,xscale,gmethBundle%FDscale),&
      useHxDef=gmethBundle%useHxDef,&
      dmeth=merge(useFDCDhybrid,fd_gmeth,gmethBundle%hybridFDCD),&
      dFDCD=gmethBundle%tolGradFDCD,&
      smeth=scal_smeth,normalize=.false.,&
      fdDirDer=slopeX,fcalls=addFcalls,err=err,message=message)
  else
    call getCDgrad(evalFunc,dataIN,dataOUT,x,fx,xscale,objFuncBundle%epsF,&
      getHxFromRelHx(gmethBundle%hx,x,xscale,gmethBundle%FDscale),&
      gmethBundle%useHxDef,gradFx,addFcalls,err,message)
    slopeX=dot_product(gradFx,sdir)
  endif
  fcalls=fcalls+addFcalls
endselect
if(err/=0)then
  retcode=unfeas_glob; message="f-linesearch_more/getGradSlopeX/&"//message
endif
! End procedure here
endsubroutine getGradSlopeX
!--
endsubroutine linesearch_more
!----------------------------------------------------
subroutine linesearch_fletcher(evalFunc,dataIN,dataOUT,x0,fx0,grad0,gmethBundle,objFuncBundle,&
            sdir,xscale,fscale,stol,alpha,beta,stepmax,x,fx,gradFx,&
            fredAct,lambda,fcalls,gcalls,retcode,message)
! Purpose: Linesearch using the strong Wolfe conditions A (alpha) and B (beta).
! Condition A (sufficient decrease condition)
! Condition B (absolute gradient condition)
! sdir is local search direction (typically Newton-derived).
! Programmer: Dmitri Kavetski
! Ref: * Fletcher,R.(1996) Practical Methods of Optimization,2nd Ed,Wiley.
!      * Nocedal,J. and Wright,S.J.(2000) Numerical Optimization, Springer.
!      * More, J.J. and Thuente, D.J. (1994) Line search algorithms with
!        guaranteed sufficient decrease, ACM Transactions on Mathematical
!        software, vol. 20(3), p.286-307.
use types_dmsl_kit,only:data_ricz_type
use utilities_dmsl_kit,only:zero,one,two,half
implicit none
! dummies
type(data_ricz_type),intent(in),optional::dataIN
type(data_ricz_type),intent(inout),optional::dataOUT
real(mrk),intent(in)::x0(:),fx0,grad0(:),xscale(:),fscale,stepmax
real(mrk),intent(in)::stol,alpha,beta
type(gmethBundle_type),intent(in)::gmethBundle
type(objFuncBundle_type),intent(in)::objFuncBundle
real(mrk),intent(inout)::lambda,sdir(:)
real(mrk),intent(out)::x(:),fx,gradFx(:),fredAct
integer(mik),intent(out)::fcalls,gcalls,retcode
character(*),intent(out)::message
! user-provided function
interface
  subroutine evalFunc(dataIN,dataOUT,x,feas,fx,gradFx,hessFx,err,message)
  use kinds_dmsl_kit
  use types_dmsl_kit,only:data_ricz_type
  implicit none
  type(data_ricz_type),intent(in),optional::dataIN
  type(data_ricz_type),intent(inout),optional::dataOUT
  real(mrk),intent(in)::x(:)
  logical(mlk),intent(out)::feas
  real(mrk),intent(out),optional::fx,gradFx(:),hessFx(:,:)
  integer(mik),intent(out)::err
  character(*),intent(out)::message
  endsubroutine evalfunc
endinterface
! local parameters
real(mrk),parameter::lambdaRedMax=0.1_mrk,lambdaRedUnfeas=0.2_mrk
real(mrk),parameter::tau1=9._mrk,tau2def=0.1_mrk,tau3=half
real(mrk),parameter::lambdaUpCoarse=2._mrk,closeToMax=0.99
real(mrk),parameter::safeEps=10._mrk
integer(mik)::status
integer(mik),parameter::keepgoing=0,dontbother=1
! locals
real(mrk)::slope0,stepLen,relLen,lambdaMin,lambdaMax,tau2
real(mrk)::slopeX,lambdaTemp,lambdaPrev,fprev,Gprev
real(mrk)::Llo,Lhi,flo,fhi,Glo,Ghi
integer(mik)::err,addFcalls
logical(mlk)::feas,ipp,finalGrad
! Start procedure here
err=0; message="linesearch_fletcher/ok"; stepLen=getStepLen2(sdir,xscale)
if(stepLen>stepmax)then ! scale down large steps
  sdir=sdir*stepmax/stepLen
  stepLen=stepmax
endif
slope0=dot_product(grad0,sdir) ! initial slope
if(slope0>=zero)then
  retcode=badDir_glob
  x=x0;fx=fx0;gradFx=grad0
  return
endif
relLen=scaledStepLen(sdir,x0,xscale) ! relative steplength (used in termination test)
lambdaMin=stol/relLen ! minimum allowable steplength, lambda must be above noise
lambdaMin=max(minval(epsRe*max(abs(x0),xscale)/max(abs(sdir),xscale)),lambdaMin)
lambdaMax=stepMax/stepLen         ! maximum allowable steplength
lambdaPrev=zero; fprev=fx0; Gprev=slope0; Llo=zero; Flo=fx0; Glo=slope0; Lhi=lambdaMax
fcalls=0;gcalls=0; retcode=failed_glob; status=keepgoing; ipp=.false.
finalGrad=.false.
! ** bracket lambda that satisfies the strong Wolfe conditions
do
  x=x0+lambda*sdir
  call getGradSlopeX()
  if(err/=0)then
    return
  elseif(.not.feas)then ! reduce lambda and try again
    Lhi=half*lambda
    lambda=max(lambda*lambdaRedMax,&
               lambdaPrev*(one-lambdaRedUnfeas)+lambdaRedUnfeas*Lhi)
    if(lambda<lambdaMin)then
      status=dontbother; retcode=unfeas_glob
      x=x0; fx=fx0; gradFx=grad0
      exit
    endif
    cycle
  endif
  if(fx>=fx0+alpha*lambda*slope0.or.(fx>=fprev.and.ipp))then ! function increasing
    Llo=lambdaPrev; flo=fprev; Glo=Gprev; Lhi=lambda; fhi=fx; Ghi=slopeX
    exit
  endif
  if(.not.ipp)ipp=.true.
  if(abs(slopeX)<=-beta*slope0)then ! strong Wolfe conditions satisfied
    status=dontbother; retcode=success_glob; exit
  elseif(slopeX>=zero)then  ! positive slope (function increasing)
    Llo=lambda; flo=fx; Glo=slopeX; Lhi=lambdaPrev; fhi=fprev; Ghi=Gprev
    exit
  elseif(lambda>=lambdaMax)then ! maximum step not big enough, it seems
    status=dontbother; retcode=success_glob; exit
  endif
  lambdaTemp=lambda*lambdaUpCoarse    ! function still decreasing: increase step ...
  lambdaTemp=max(lambdaTemp,two*lambda-lambdaPrev) ! ... using safeguards
  lambdaTemp=min(lambdaTemp,lambda+tau1*(lambda-lambdaPrev))
  if(lambdaTemp>=closeToMax*Lhi)then  ! very close to max step
    status=dontbother; retcode=success_glob; exit
  endif
  lambdaprev=lambda; fprev=fx ! keep previous evaluation
  lambda=lambdatemp           ! and update steplength
  if(lambda<lambdaMin)then    ! flag tiny steps getting nowhere
    status=dontbother; retcode=failed_glob
    if(fx>fx0)then  ! return original point if it was the best found
      x=x0; fx=fx0; gradFx=grad0
    endif
    exit
  elseif(lambda>lambdaMax)then  ! flag large steps
    lambda=lambdaMax; Lhi=lambdaMax
  endif
enddo
! "zoom", or bracket contraction
tau2=min(tau2def,beta)
do
  if(status==dontbother)exit
  lambdaprev=lambda; fprev=fx; Gprev=slopeX
  lambdaTemp=half*(Llo+Lhi)                     ! bisection
  lambdaTemp=max(lambdaTemp,Llo+tau2*(Lhi-Llo)) ! use safeguards
  lambda=min(lambdaTemp,Lhi-tau3*(Lhi-Llo))
! evaluate function
  x=x0+lambda*sdir
  call getGradSlopeX()
  if(err/=0)then
    return
  elseif(.not.feas)then ! unfeasible inside bracket: too hard basket
    retcode=unfeas_glob
    x=x0; fx=fx0; gradFx=grad0
    exit
  endif
  if(fx>=fx0+alpha*lambda*slope0.or.fx>=Flo)then ! function increasing
    Lhi=lambda; fhi=fx
  else
    if(abs(slopeX)<=-beta*slope0)then ! satisfaction of strong gradient condition
      retcode=success_glob; exit
    elseif(slopeX*(Llo-lambda)<zero)then
      Lhi=Llo; Fhi=Flo
    endif
    Llo=lambda;Flo=fx
  endif
! - check for convergence
  if(abs(Lhi-Llo)<=safeEps*epsRe*max(Lhi,Llo))then   ! brackets merged
    retcode=success_glob; exit      ! (not much point persevering...)
  elseif(half*(Llo+Lhi)<lambdaMin)then
    if(fx>fx0)then  ! restore original point if all trials worse...
      x=x0;fx=fx0;gradFx=grad0
    endif
    retcode=failed_glob;  exit      ! bracket collapsed to initial point ...
  endif
enddo
sdir=x-x0 ! get shift vector (must equal sdir*lambda)
fredAct=fx0-fx     ! actual reduction in function value
if(gmethBundle%useDirDer)then ! final call to get full gradient
  finalGrad=.true.; call getGradSlopeX()
endif
! End procedure here
contains
!--
subroutine getGradSlopeX()  ! macro to get directional derivative
use utilities_dmsl_kit,only:getFDCDgrad,getCDgrad,getHxFromRelHx,getFDdirDer
implicit none
! dummies
! locals
! local registered settings
integer(mik),parameter::scal_smeth=0,imax_smeth=1,ave_smeth=2,wei_smeth=3
! Start procedure here
selectcase(gmethBundle%gmeth_now)
case(user_meth)   ! analytical derivatives available
  call evalFunc(dataIN,dataOUT,x,feas,fx,gradFx,err=err,message=message);fcalls=fcalls+1;gcalls=gcalls+1
  if(err/=0)then
    message="f-linesearch_fletcher/getGradSlopeX/userErrA/&"//message
    retcode=badFunc_glob; return
  elseif(.not.feas)then
    message="f-linesearch_fletcher/getGradSlopeX/userUnfeas/&"//message; return
  endif
  slopeX=dot_product(gradFx,sdir) ! new slope
case(fd_gmeth)    ! forward difference gradient
  if(.not.finalGrad)then
    call evalFunc(dataIN,dataOUT,x,feas,fx,err=err,message=message); fcalls=fcalls+1
    if(err/=0)then
      message="f-linesearch_fletcher/getGradSlopeX/userErrFD/&"//message
      retcode=badFunc_glob; return
    elseif(.not.feas)then
      message="f-linesearch_fletcher/getGradSlopeX/userUnfeasFD/&"//message; return
    endif
  endif
  if(.not.finalGrad.and.gmethBundle%useDirDer)then ! use cheap directional derivative
    call getFDdirDer(evalFunc,dataIN,dataOUT,x=x,p=sdir,fx=fx,xscale=xscale,fscale=fscale,&
      epsF=objFuncBundle%epsF,&
      hx=getHxFromRelHx(gmethBundle%hx,x,xscale,gmethBundle%FDscale),&
      useHxDef=gmethBundle%useHxDef,&
      dmeth=merge(useFDCDhybrid,fd_gmeth,gmethBundle%hybridFDCD),&
      dFDCD=gmethBundle%tolGradFDCD,&
      smeth=scal_smeth,normalize=.false.,&
      fdDirDer=slopeX,fcalls=addFcalls,err=err,message=message)
  else
    call getFDCDgrad(evalFunc,dataIN,dataOUT,x,fx,xscale,fscale,objFuncBundle%epsF,&
      getHxFromRelHx(gmethBundle%hx,x,xscale,gmethBundle%FDscale),&
      gmethBundle%useHxDef,&
      merge(useFDCDhybrid,fd_gmeth,gmethBundle%hybridFDCD),gmethBundle%tolGradFDCD,&
      gradFx,addFcalls,err,message)
    slopeX=dot_product(gradFx,sdir)
  endif
  fcalls=fcalls+addFcalls
case(cd_gmeth)    ! central difference gradient
  if(.not.finalGrad)then
    call evalFunc(dataIN,dataOUT,x,feas,fx,err=err,message=message); fcalls=fcalls+1
    if(err/=0)then
      message="f-linesearch_fletcher/getGradSlopeX/userErrCD/&"//message
      retcode=badFunc_glob; return
    elseif(.not.feas)then
      message="f-linesearch_fletcher/getGradSlopeX/userUnfeasCD/&"//message; return
    endif
  endif
  if(.not.finalGrad.and.gmethBundle%useDirDer)then ! use cheap directional derivative
    call getFDdirDer(evalFunc,dataIN,dataOUT,x=x,p=sdir,fx=fx,xscale=xscale,fscale=fscale,&
      epsF=objFuncBundle%epsF,&
      hx=getHxFromRelHx(gmethBundle%hx,x,xscale,gmethBundle%FDscale),&
      useHxDef=gmethBundle%useHxDef,&
      dmeth=merge(useFDCDhybrid,fd_gmeth,gmethBundle%hybridFDCD),&
      dFDCD=gmethBundle%tolGradFDCD,&
      smeth=scal_smeth,normalize=.false.,&
      fdDirDer=slopeX,fcalls=addFcalls,err=err,message=message)
  else
    call getCDgrad(evalFunc,dataIN,dataOUT,x,fx,xscale,objFuncBundle%epsF,&
      getHxFromRelHx(gmethBundle%hx,x,xscale,gmethBundle%FDscale),&
      gmethBundle%useHxDef,gradFx,addFcalls,err,message)
    slopeX=dot_product(gradFx,sdir)
  endif
  fcalls=fcalls+addFcalls
endselect
if(err/=0)then
  retcode=unfeas_glob; message="f-linesearch_fletcher/getGradSlopeX/&"//message
endif
! End procedure here
endsubroutine getGradSlopeX
!--
endsubroutine linesearch_fletcher
!----------------------------------------------------
subroutine brentmin(evalFunc,dataIN,dataOUT,linmin_ometh,xopt,fold,sdir,stpmax,stol,Ltol,itmax,xscale,&
  fopt,lambda,fcalls,gcalls,retcode,message)
! Purpose: Brent line minimisation: search from xopt in direction sdir.
use types_dmsl_kit,only:data_ricz_type
use utilities_dmsl_kit,only:zero,one,two,assertEq
use numerix_dmsl_kit,only:linmin
implicit none
! Dummies
type(data_ricz_type),intent(in),optional::dataIN
type(data_ricz_type),intent(inout),optional::dataOUT
real(mrk),intent(inout)::xopt(:),sdir(:),lambda
real(mrk),intent(in)::fold
integer(mik),intent(in)::linmin_ometh
real(mrk),intent(in)::stpmax,stol,Ltol,xscale(:)
integer(mik),intent(in)::itmax
real(mrk),intent(inout)::fopt
integer(mik),intent(out)::fcalls,gcalls
integer(mik),intent(out)::retcode
character(*),intent(out)::message
interface
  subroutine evalFunc(dataIN,dataOUT,x,feas,fx,gradFx,hessFx,err,message)
  use kinds_dmsl_kit
  use types_dmsl_kit,only:data_ricz_type
  implicit none
  type(data_ricz_type),intent(in),optional::dataIN
  type(data_ricz_type),intent(inout),optional::dataOUT
  real(mrk),intent(in)::x(:)
  logical(mlk),intent(out)::feas
  real(mrk),intent(out),optional::fx,gradFx(:),hessFx(:,:)
  integer(mik),intent(out)::err
  character(*),intent(out)::message
  endsubroutine evalfunc
endinterface
! Locals
integer(mik)::err,ndum
real(mrk)::lambdaMin,lambdaMin1,lambdaMin2,relLen
logical(mlk),parameter::useFold=.true.
logical(mlk)::ok
! Start procedure here
call assertEq(size(xopt),size(sdir),size(xscale),ok,ndum)
if(ok)then
  retcode=success_glob; message="brentmin/ok"
else
  retcode=bugFail; message="f-brentmin/dimError"; return
endif
relLen=scaledStepLen(sdir,xopt,xscale) ! relative steplength (used in termination test)
if(relLen<=stol)then
  retcode=success_glob; message="w-brentmin/zeroLen[sdir<=stol]"
  fopt=fold; fcalls=0; gcalls=0; return
endif
lambdaMin1=stol/relLen ! minimum allowable steplength based on convergence test
lambdaMin2=epsRe*minval(max(abs(xopt),xscale)/max(abs(sdir),xscale)) ! lambda must be above noise
lambdaMin=max(lambdaMin1,lambdaMin2); fopt=fold
call linmin(evalFunc,dataIN,dataOUT,linmin_ometh,xopt,sdir,-stpmax,+stpmax,Ltol,xscale,&
            itmax,lambda,fopt,useFold,fcalls,gcalls,err,message)
if(err/=0)then
  retcode=failed_glob; message="f-brentmin/&"//message; return
elseif(lambda<zero)then       ! search went backwards ... hmmm ...
  retcode=bugFail; message="f-brentmin/weird[lambda<0]/&"//message; return
elseif(lambda<=lambdaMin)then ! search back-collapsed to initial point ...
  retcode=failed_glob; message="f-brentmin/bracketCollapsed/&"//message; return
endif
! End procedure here
endsubroutine brentmin
!----------------------------------------------------
subroutine trustDriver(evalFunc,dataIN,dataOUT,x0,fx0,grad0,hess0,Ld0,hessScaled,&
  boundedSearch,xLo,xHi,activeSet,&
  imeth,hmeth,quadTypeH,maxSR1update,gmethBundle,xscaleHmeth,xscale,fscale,&
  trustBundle,objFuncBundle,hessFacBundle,didGradNewHess,&
  x,fx,gradx,dx,trustRad,fredExp,fredAct,fcalls,gcalls,ncholstats,&
  logdet,condest,Einf,err,message)
! --
! Purpose: Implements a complete trust region step:
! 1. Solve scaled trust region subproblem
! 2. Accept or reject trust region step
! 3. Update trust region.
! --
! Available methods
! * Near-exact hookstep of More and Sorensen: ~2-5 Cholesky per hook,
!   handles arbitrary Hessians, but makes most sense working with exact Hessian
!   (or SR1) and does not work with factored quasi-Hessians.
! * Generalized dogleg step (2D subspace minimization): <=1 Cholesky per dog,
!   handles arbitrary Hessians, can handle factored/unfactored Hessians.
!   approximate trust-region solution, best for quasi-Newton methods since
!   then it makes little sense to solve the trust problem to great accuracy.
! --
! Programmer: Dmitri Kavetski, 1 April 2004.
! Ref: Nocedal and Wright, 1999.
!      Dennis and Schanbel,1996.
! --
! Comments:
! * For bound-constrained optimisation this routine constructs the masked
!   active Hessian on each call, which can be expensive if the problem is
!   very highly dimensional and only a few constraints are active
! * Trust updating is normally done within the subsidiary routine updateTrust.
!   However, this routine does some basic operations accounting for active
!   set constraints (truncating steps that reach beyond bounds, etc.).
! * For extremely expensive functions it may be preferable to conduct a
!   near-exact curvilinear trust region solution search, particularly near
!   bounds. The cost of this strategy is a potentially large number of Cholesky
!   solutions per outer iteration. The "internal doubling" strategy implemented
!   here is probably nearly as good in the unconstrained case.
! * When working near bounds, the trust region is not altered, since reducing
!   the radius of a spherical trust region may unnecesarily constrain steps
!   away from bounds. A box-step trust region may be more appropriate if
!   spherical trust becomes too inefficient (behaves more like linesearch).
use types_dmsl_kit,only:data_ricz_type
use utilities_dmsl_kit,only:zero,one,terminateRowColMat,arthsi,quadDf
implicit none
! dummies
type(data_ricz_type),intent(in),optional::dataIN
type(data_ricz_type),intent(inout),optional::dataOUT
real(mrk),intent(in)::x0(:),fx0,grad0(:)
real(mrk),intent(inout)::hess0(:,:),Ld0(:)
real(mrk),intent(inout)::hessScaled(:,:) ! scratch Hessian
type(trustBundle_type),intent(in)::trustBundle
type(gmethBundle_type),intent(in)::gmethBundle
type(objFuncBundle_type),intent(in)::objFuncBundle
type(hessFacBundle_type),intent(in)::hessFacBundle
integer(mik),intent(in)::imeth,hmeth
character(*),intent(in)::quadTypeH
logical(mlk),intent(in)::maxSR1update
logical(mlk),intent(out)::didGradNewHess  ! only used for SR1 updating
integer(mik),intent(in)::xscaleHmeth
logical(mlk),intent(in)::boundedSearch
integer(mik),optional,intent(in)::activeSet(:)
real(mrk),optional,intent(in)::xLo(:),xHi(:)
real(mrk),intent(in)::xscale(:),fscale
real(mrk),intent(inout)::trustRad
real(mrk),intent(out)::x(:),dx(:),gradx(:),fx,fredExp,fredAct
real(mrk),intent(out)::logdet,condest,Einf
integer(mik),intent(out)::fcalls,gcalls,ncholstats(:)
integer(mik),intent(out)::err
character(*),intent(out)::message
! user-provided function
interface
  subroutine evalFunc(dataIN,dataOUT,x,feas,fx,gradFx,hessFx,err,message)
  use kinds_dmsl_kit
  use types_dmsl_kit,only:data_ricz_type
  implicit none
  type(data_ricz_type),intent(in),optional::dataIN
  type(data_ricz_type),intent(inout),optional::dataOUT
  real(mrk),intent(in)::x(:)
  logical(mlk),intent(out)::feas
  real(mrk),intent(out),optional::fx,gradFx(:),hessFx(:,:)
  integer(mik),intent(out)::err
  character(*),intent(out)::message
  endsubroutine evalfunc
endinterface
!--locals for trust driver algorithm
integer(mik)::ndim,nchol,addFcalls,i,itot,lerr,stepResult
real(mrk)::xtry(size(x)),fxtry,stepLenB4Big
logical(mlk)::firstTrustIter,reducedTrust,expandingTrust,redoHess
character(100)::lmessage
!--common locals for trust update
real(mrk)::redRatio,redExp
!--common locals for trust solvers
real(mrk)::newtStep(size(x)),newtLen,negStep(size(x))
real(mrk)::stepLen,logdetTemp,condestTemp,EinfTemp
!--locals for bounds
integer(mik)::nact,avar(size(x))
logical(mlk)::active(size(x)),hitBound
real(mrk)::gradScaled(size(x)),xscaleH(size(x)),gradTry(size(x))
!--locals for hook solver
type hook_type
  real(mrk)::lambdaPD   = zero
  real(mrk)::negLen     = undefRN
endtype hook_type
type(hook_type)::hookNew,hook
!--locals for dogleg solver
real(mrk)::LdScaled(size(x))
type doglegVars_type
  logical(mlk)::haveFac =.false.
  logical(mlk)::haveNewt=.false.
  logical(mlk)::haveNeg =.false.
  logical(mlk)::haveGBG =.false.
  logical(mlk)::posDef  =.false.
  real(mrk)::normG      = undefRN
  real(mrk)::gBg        = undefRN
  real(mrk)::absGdotNewt= undefRN
endtype doglegVars_type
type(doglegVars_type)::dogNew,dog
logical(mlk)::useL
integer(mik)::indx(size(x))
! local pars
integer(mik),parameter::itotMax=100 ! prevents (?) bugs in infinite looop.
! Start procedure here
err=failed_glob;message="trustDriver/started(notFinished)"
ndim=size(x);fcalls=0;gcalls=0;i=0;itot=0;ncholstats=0
firstTrustIter=.true.;reducedTrust=.false.;expandingTrust=.false.;didGradNewHess=.false.
newtLen=min(-one,undefRN);stepLenB4Big=zero;hitBound=.false.
useL=imeth==dogLeg_imeth.and.hmeth==bfgsFac_hmeth
!logdetTemp=logdet;condestTemp=condest;EinfTemp=Einf
hook=hookNew; dog=dogNew; indx=-20  ! initialize solver vars
do      ! this loop is not fixed length since it is possible to expand trust region
  i=i+1; itot=itot+1  ! if function behaviour agreement is good
  if(i>trustBundle%niter_tr)then
    err=failed_glob;message="trustDriver/tooManyBadTrustTries"
    exit
  elseif(itot>itotMax)then
    err=bugFail;message="f-trustDriver/stuckInLoop?itotMaxExceeded"
    return
  endif
! 1. solve trust region subproblem
  if(present(activeSet))then  ! * bound-contrained optimisation
    active=activeSet==freeVar_as
    nact=count(active)
  else
    nact=ndim
  endif
  if(nact==0)then         ! * all variables fixed
    dx=zero;redExp=zero;logdetTemp=zero;condestTemp=zero;nchol=zero
    err=bugFail;message="f-trustDriver/allVarsFixed"
    return
  endif
  selectcase(hmeth)
  case(bfgsFac_hmeth)         ! ** Factored Hessian (BFGS only)
    if(imeth==trustEx_imeth)then            ! - hook step unsupported
      err=bugFail;message="f-trustDriver/invalidIN:factoredHookstep"
      return
    elseif(trustBundle%pivotCholTrust)then  ! - pivoting unsupported
      err=bugFail;message="f-trustDriver/invalidIN:pivotedFactoredDog"
      return
    elseif(xscaleHmeth==xscaleH_hdiag)then  ! - diagonal Hessian scaling unsupported
      err=bugFail;message="f-trustDriver/invalidIN:factoredDog:xscaleH=hdiag"
      return
    endif
    redoHess=(itot==1.or.didGradNewHess)  ! scale (factored) scratch Hessian
    if(redoHess)then  ! (currently didGradNewHess===false for BFGS)
      call getXscaleH(xscaleHmeth,hess0,xscale,fscale,xscaleH)
      if(nact==ndim)then  ! - effectively unconstrained
        avar=arthsi(ndim)
        call xscaleNewt(xscaleH,L0=hess0,LS=hessScaled,Ld0=Ld0,LdS=LdScaled)
      else                ! - active constraints present: need to muck around
        avar=pack(arthsi(ndim),active) ! index of active variables
        call terminateRowColMat(hess0,hessScaled(1:nact,1:nact),active,lerr,lmessage)
        LdScaled=pack(Ld0,activeSet==freeVar_as)
        call xscaleNewt(xscaleH(avar(1:nact)),&
                        LS=hessScaled(1:nact,1:nact),LdS=LdScaled(1:nact))
      endif
      call xscaleNewt(xscaleH(avar(1:nact)),&
                      grad0=grad0(avar(1:nact)),gradS=gradScaled(1:nact))
      dog%haveFac=.true.  ! factored (BFGS) Hessian available immediately
      dog%posDef=.true.   ! and is always positive definite
    endif
  case default                ! ** Unfactored Hessians (all others)
    redoHess=(itot==1.or.didGradNewHess)  ! scale Hessian on first iteration or if refreshed
    if(redoHess)then      ! reconstruct active scaled Hessian
      call getXscaleH(xscaleHmeth,hess0,xscale,fscale,xscaleH)
      if(nact==ndim)then  ! - effectively unconstrained
        avar=arthsi(ndim)
        call xscaleNewt(xscaleH,hess0=hess0,hessS=hessScaled)
      else                ! - active constraints present: need to muck around
        avar=pack(arthsi(ndim),active) ! index of active variables
        call terminateRowColMat(hess0,hessScaled(1:nact,1:nact),active,lerr,lmessage)
        call xscaleNewt(xscaleH(avar(1:nact)),hessS=hessScaled(1:nact,1:nact))
      endif
      call xscaleNewt(xscaleH(avar(1:nact)),&
                      grad0=grad0(avar(1:nact)),gradS=gradScaled(1:nact))
      dog=dogNew;hook=hookNew    ! reset all dogs and hooks since Hessian is new
    endif
    selectcase(imeth)
    case(trustEx_imeth) ! * Near-exact hook step trust only handles unfactored Hessians
      call solveTrustHook(B=hessScaled(1:nact,1:nact),grad=gradScaled(1:nact),& ! data
        doPivot=trustBundle%pivotCholTrust,trustRad=trustRad,&  ! pivoting
        ncholMax=trustBundle%ncholMax_tr,lambda=hook%lambdaPD,& ! trust settings
        psol=dx(1:nact),pnorm=stepLen,stepResult=stepResult,&   ! hookstep and its length
        logdet=logdetTemp,condest=condestTemp,Einf=EinfTemp,&   ! properties of Hessian
        firstHook=redoHess,newtStep=newtStep(1:nact),newtLen=newtLen,&   ! full Newton step on first call
        negStep=negStep(1:nact),negLen=hook%negLen,&
        nchol=nchol,err=lerr,message=lmessage)
    endselect
  endselect
! - if dogleg requested compute it from scaled factored/unfactored Hessians
  selectcase(imeth)
  case(dogLeg_imeth)    ! * Dogleg trust: factored/unfactored input Hessian
      call solveGeneralDogTrust(&
        B=hessScaled(1:nact,1:nact),Ld=LdScaled(1:nact),&
        indx=indx(1:nact),grad=gradScaled(1:nact),&
        trustRad=trustRad,dogNewtBias=trustBundle%dogNewtBias,&
        haveFac=dog%haveFac,haveNewt=dog%haveNewt,haveNeg=dog%haveNeg,&
        doPivot=trustBundle%pivotCholTrust,hessFacBundle=hessFacBundle,&
        haveGBG=dog%haveGBG,useL=useL,posDef=dog%posDef,&
        newtStep=newtStep(1:nact),newtLen=newtLen,&
        negEigen=EinfTemp,negStep=negStep(1:nact),&
        normG=dog%normG,gBg=dog%gBg,absGdotNewt=dog%absGdotNewt,&
        logdet=logdetTemp,condest=condestTemp,Einf=EinfTemp,&
        psol=dx(1:nact),pLen=stepLen,stepResult=stepResult,&
        nchol=nchol,err=lerr,message=lmessage)
  case(trustEx_imeth)   ! (-) Hook step already computed
  case default
    err=bugFail;message="trustDriver/BUG/unknownIMETH"
    return
  endselect
! - basic check of trust solver
  selectcase(lerr)
  case(okAlg)   ! - sucesful completion
    ncholstats(1)=ncholstats(1)+nchol ! total number of Cholesky factorizations
    ncholstats(2)=ncholstats(2)+1     ! number of "minor" trust region iterations
    call unXscaleNewt(xscaleH(avar(1:nact)),p0=dx(1:nact))  ! unscale trust step
    if(nact<ndim)dx=unpack(dx(1:nact),active,zero)
  case default  ! - some kind of error: most likely a bug
    err=bugFail;message="trustDriver/BUG/&"//lmessage
    return
  endselect
! - reduce trust radius (4 safety) if natural (Newton) step is well inside radius
  selectcase(stepResult)
  case(insideTrust)           ! prevents trust radius from creeping up
    if(trustRad>stepLen*trustBundle%trustOstepMax_tr)then
     trustRad=stepLen*trustBundle%trustOstepMax_tr
     if(expandingTrust)reducedTrust=.true.
    endif
  case(onTrustBound,hardCase) ! keep trust region intact (for now...)
  case default
    err=bugFail;message="trustDriver/BUG/unknownStepRes/&"//lmessage
    return
  endselect
! - check for bound violation
  if(boundedSearch)then ! may need 2 truncate step when colliding with bounds
    call checkStepBounds(x0,xLo,xHi,activeset,dx,hitBound=hitBound)
  elseif(present(xLo).or.present(xHi).or.present(activeset))then
    err=10;message="trustDriver/inError/bug/bothBoundsMustBePresent";return
  endif
! - safeguarded scaled step length and expected reduction in quadratic model
  stepLen=getStepLen2(dx,xscaleH)
  redExp=-quadDf(dx=dx,dfdx=grad0,d2fdx2=hess0,typeH=quadTypeH)
  if(firstTrustIter)then
! save original expected function reduction given initial trust:
! this is used when assessing convergence of the globalisation strategy
    fredExp=redExp
  endif
  if(firstTrustIter)then ! store Hessian properties (later values will be affected by
    logdet=logdetTemp;condest=condestTemp;Einf=EinfTemp ! the trust region solution)
  endif
  firstTrustIter=.false.
! bounds grossly interfering with trust expansion: dont bother checking new point
  if(expandingTrust.and.hitBound.and.stepLen<=trustBundle%boundFrac*stepLenB4Big)then
    x=xtry; fx=fxtry; gradx=gradTry; dx=xtry-x0 ! fall back unto previous trust iteration
    err=success_glob; message="trustDriver/ok/expansionFailed(bounds)"
    exit
  endif
! 2. Accept / Reject solution and Update trust region
  if(nact==0)then     ! all vars fixed
    x=x0; fx=fx0; redRatio=zero; addFcalls=0
    err=success_glob; message="w-trustDriver/ok/&"//lmessage
    exit
  endif
  call updateTrust(evalFunc,dataIN,dataOUT,x0,fx0,grad0,dx,stepLen,stepResult,redExp,&
                   xscale,fscale,reducedTrust,objFuncBundle,trustBundle,&
                   x,fx,trustRad,redRatio,addFcalls,lerr,lmessage)
  fcalls=fcalls+addFcalls
  selectcase(lerr)
  case(suceed_tr)     ! succesful iteration
    call checkSR1updt()
    err=success_glob;message="trustDriver/ok/&"//lmessage
    exit
  case(goBig_tr)      ! succesful iteration, but re-take step with larger trust
    call checkSR1updt()
    xtry=x; fxtry=fx; gradTry=gradx; expandingTrust=.true.
    stepLenB4Big=stepLen  ! steplength before go-big step
    i=i-1 ! do not count trust expansion as a trust "try" since it is a good thing!
  case(collapsed_tr)  ! collapsed trust region
    if(redRatio<=zero)then ! reset current point of trial point worse than current
      x=x0; fx=fx0; gradx=grad0; dx=zero
    endif
    err=failed_glob;message="trustDriver/&"//lmessage
    exit
  case(fconExpObs_tr) ! exp/obs reduction within function precision
    didGradNewHess=.false.;err=success_glob;message="trustDriver/&"//lmessage
    if(redRatio<=zero)then ! reset current point of trial point worse than current
      x=x0; fx=fx0; gradx=grad0; dx=zero
    endif
    exit
  case(blown_tr)      ! trust region blown
    call checkSR1updt()
    err=success_glob;message="trustDriver/&"//lmessage
    exit
  case(unfeas_tr)     ! unfeasible iteration: keep going...
    if(expandingTrust)then  ! ...fall back on pre-expanded results
      x=xtry; fx=fxtry; gradx=gradTry; dx=xtry-x0 ! and do not update SR1 Hessian
      err=success_glob; message="trustDriver/ok/expansionFailed(unfeas)"
      exit
    else
      dx=zero; err=unfeas_glob
    endif
  case(failed_tr)     ! failed iteration:
    call checkSR1updt()
    if(expandingTrust)then  ! ...fall back on pre-expanded results
      x=xtry; fx=fxtry; gradx=gradTry; dx=xtry-x0
      err=success_glob; message="trustDriver/ok/expansionFailed(normal)"
      exit
    else                            ! ...or keep going with reduced trust
      dx=zero;reducedTrust=.true.   !    (preventing agressive go-big steps)
      if(redRatio<=zero)then ! reset current point if trial point worse than current
        x=x0; fx=fx0; gradx=grad0; dx=zero
      endif
      err=failed_glob; message="trustDriver/trustFailed/&"//lmessage
    endif
  case(dxTiny_tr)     ! negligible step suggested: exit with ...
    didGradNewHess=.false.;message=trim(message)//"/&"//lmessage  ! ... previous code (either unfeas or failed)
    exit
  case(expRedNonP_tr) ! expected reduction nonpositive: error?
    selectcase(stepResult)
    case(hardCase)  ! hard case a bit dubious: force "steeper-descent" move
      if(expandingTrust)then  ! ...fall back on pre-expanded results
        x=xtry; fx=fxtry; gradx=gradTry; dx=xtry-x0
        err=success_glob; message="trustDriver/ok/expansionFailed(HCfail)"
        exit
      else ! reset current point if trial point worse than current
        reducedTrust=.true.   !    (prevent agressive go-big steps)
        dx=zero; x=x0; fx=fx0; gradx=grad0; dx=zero
        err=failed_glob
        trustRad=trustRad*trustBundle%radDown_tr
      endif
!    if(present(activeSet).and.&   ! hard-case possibly interfering: try reducing trust
!       stepResult==hardCase)then  ! which forces a more "steepest-descent" move
!      err=failed_glob
!      trustRad=trustRad*trustBundle%radDown_tr
    case default
      err=bugFail; message="trustDriver/BUG?/&"//lmessage
      dx=zero
      exit
    endselect
  case default
    dx=zero; err=bugFail
!    err=failed_glob
    write(message,'(a,i0,a)')"trustDriver/unknown/[lerr=",lerr,"]/&"//trim(lmessage)
    exit
  endselect
enddo
fredAct=fx0-fx       ! actual reduction in function value
! End procedure here
contains
!----
subroutine checkSR1updt() ! macro to check internal SR1 update
use utilities_dmsl_kit,only:getHxFromRelHx,getFDCDgrad,getCDgrad
implicit none
! locals
logical(mlk)::feas
integer(mik)::jerr
! Start procedure here
feas=.true.; jerr=0
didGradNewHess=hmeth==SR1unFac_hmeth.and.&  ! SR1 quasi-Hessian
        (maxSR1update.or.&  ! - maximal-frequency updating requested or model way off
        ((fx0-fx)-redExp)<trustBundle%SR1forceUpdt*max(abs(fx0),fscale))
if(didGradNewHess)then
  selectcase(gmethBundle%gmeth_now)
  case(user_meth)     !   - user-supplied gradient
    call evalFunc(dataIN,dataOUT,x,feas,gradFx=gradx,err=err,message=lmessage); gcalls=gcalls+1
    if(err/=0)then        ! strange, since point x has already been trialled with err=0
      message="f-trustDriver/getGrad/userBug/userErr/&"//lmessage
      err=bugFail; return
    elseif(.not.feas)then ! strange, since point x already tested ok...
      message="f-trustDriver/getGrad/userBug/userUnfeas/&"//lmessage
      err=bugFail; return
    endif
  case(fd_gmeth)      !   - FD/CD method
    call getFDCDgrad(evalFunc,dataIN,dataOUT,x,fx,xscale,fscale,objFuncBundle%epsF,&
      getHxFromRelHx(gmethBundle%hx,x,xscale,gmethBundle%FDscale),gmethBundle%useHxDef,&
      merge(useFDCDhybrid,fd_gmeth,gmethBundle%hybridFDCD),gmethBundle%tolGradFDCD,&
      gradx,addFcalls,jerr,lmessage)
    fcalls=fcalls+addFcalls
  case(cd_gmeth)      !   - CD method
    call getCDgrad(evalFunc,dataIN,dataOUT,x,fx,xscale,objFuncBundle%epsF,&
      getHxFromRelHx(gmethBundle%hx,x,xscale,gmethBundle%FDscale),gmethBundle%useHxDef,&
      gradx,addFcalls,jerr,lmessage)
    fcalls=fcalls+addFcalls
  endselect
  if(.not.feas.or.jerr/=0)then
    err=bugFail
    write(message,'(a,i0,a)')"f-trustDriver/getGrad/unfeas??/BUG?/",jerr,"/&"//trim(lmessage)
  endif
  call SR1unFac_update(dx,xscale,activeSet,gradx,grad0,trustBundle%SR1skipTol,&
                       hess0,rescale=.false.,err=err,message=message)
endif
! End procedure here
endsubroutine checkSR1updt
!----
endsubroutine trustDriver
!----------------------------------------------------
pure subroutine solveTrustHook(B,grad,doPivot,trustRad,&
    ncholMax,lambda,psol,Pnorm,stepResult,&
    logdet,condest,Einf,firstHook,newtStep,newtLen,negStep,negLen,&
    nchol,err,message)
! Purpose: Solves the scaled trust region problem using the exact "hookstep" method.
! Algorithm 3.14 of More and Sorensen,1983 (see also Nocedal and Wright, 2000).
! For bounded problems expects the active Hessian submatrix only is supplied.
!---
! Input:  B       = scaled Hessian (or quasi-Hessian): full usage
!         grad    = gradient vector
!         trustRad= size of trust region
!         lambda  = initial diagonal increment at solution
! Output:
!         lambda  = final diagonal increment at solution
!         psol    = solution of the trust region: optimal constrained step
!         firstHook = (a) indicates first hook attempt with lambda=0, in this case
!                     return full Newton step newtStep. If B not +ve definite
!                     newtStep is undefined and newtLen set <0. In the hard case
!                     returns the negative curvature step.
!                     (b) on subsequent calls (expanding/contracting iterations)
!                     will check whether this information can be re-used.
!         Pnorm   = length of constrained step (could be <trustRad in the Newton case)
!         nchol   = number of Cholesky factorisations performed
!         err     = output status
! If trust region problem cannot be solved in ncholmax Cholesky iterations
! then procedure exits with error message.
! Ref:
! * More and Sorensen (1983) Computing a trust region step,
!   SIAM Journal of Scientific and Statistical Computing,4(3),553-572.
! Notes:
! * The algorithm works in feasible space only, and has no knowledge of any
!   bounds and constrains. These should be processed after the solution is
!   computed. In some cases a box-step may be more appropriate, but requires
!   a different algorithm.
! * The routine also returns estimates 'condest','logdet' of the input 
!   Hessian plus initial modification. On the first call to solveTrustHook()
!   the modification is zero and hence condest and logdet refer to the current
!   Hessian estimate. When solving bound-constrained problems, the reduced Hessian
!   is supplied, so that condest and logdet refer to the active sub-Hessian only.
! * Einf estimates the most negative eigenvalue of the input matrix.
! * The "hard-case" arise in the vicinity to saddle points or where the gradient
!   is small or near points with negative curvature. This code computes the
!   eigenvector of the smallest Hessian eigenvalue and upscales it to make up
!   the trust region step.
! * The method is iterative and each iteration requires a Cholesky factorization
!   (adding diagonal lambda). Although typically 1-3 iterations are needed, for
!   large problems this may be unacceptable. For such problems use
!   (i)  dogleg step with unfactored matrix, single Cholesky factorization,
!        which retains most if not all advantages of trust regions (provided
!        indefinite Hessians are generated near saddle points).
!   (ii) dogleg step with factored quasi-Newton (BFGS) Hessian, no
!        Cholesky factorizations but looses robustness near saddle-points,
!        since negative curvature information is not available.
! * doPivot requests pivoted Cholesky decompositions (provided allowPivot=.true.)
!   In the hook trust region algorithm indefinite Cholesky are handled using
!   diagonal addition until positive definite result achieved. Since pivoting
!   is not needed for positive-definite Cholesky methods, it is unlikely that
!   the hook step benefits from pivoting (unlike the dogstep, where pivoting
!   may help determine the eigenvalues and negative curvature). The disadvantage
!   of pivoting is that it may require lots of data-moving. Overall, not recommended.
use utilities_dmsl_kit,only:zero,half,one,two,norm2,getDiag,getKnorm,&
  lower_rsolv,putDiag,addDiag,flip_UtoL
use linalg_dmsl_kit,only:choles_dcmp,choles_negEigVec,choles_fwbw
implicit none
! dummies
real(mrk),intent(inout)::B(:,:),lambda
integer(mik),intent(inout)::stepResult
logical(mlk),intent(in)::doPivot
real(mrk),intent(in)::grad(:),trustRad
integer(mik),intent(in)::ncholMax
real(mrk),intent(out)::psol(:),Pnorm
real(mrk),intent(out)::logdet,condest,Einf
logical(mlk),intent(in)::firstHook
real(mrk),intent(inout)::newtStep(:),newtLen,negStep(:),negLen
integer(mik),intent(out)::nchol
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
integer(mik)::indx(size(grad)),ibad
real(mrk)::Hdiag(size(grad)),q(size(grad)),Ld(size(grad)),zHat(size(grad))
real(mrk)::L_lo,L_hi,L_s,Gnorm,Bnorm1,Qnorm,tau,RZ2,RP2,eigval
real(mrk)::logdetTemp,condestTemp
logical(mlk)::pd
integer(mik)::jerr
character(100)::jmessage
! local parameters
real(mrk),parameter::L_s_red=1.e-3_mrk,tolGnorm=1.e2_mrk*epsRe
real(mrk),parameter::powLhi=one,rootLhi=one/(one+powLhi)  ! powLhi>1 biases faster L increase
real(mrk),parameter::sigma1=0.1_mrk,sigma2=0._mrk
logical(mlk),parameter::normW=.false.
real(mrk),parameter::cholTau=zero,cholTauBar=zero,cholMu=zero ! robust Cholesky
integer(mik),parameter::largeNorm=1 ! method for eigenvector estimation
logical(mlk),parameter::allowNewtReuse=.true.,allowNegReuse=.false.
real(mrk),parameter::sigma3reuseHard=0.3
logical(mlk),parameter::usePivotDef=.false.,allowPivot=.true.
logical(mlk)::usePivot
! Start procedure here
! 0. Initial brackets for lambda
usePivot=merge(doPivot,usePivotDef,allowPivot)
if(trustRad<=zero)then  ! this caused DK a lot of grief once...
  err=failAlg;psol=zero;Pnorm=zero;stepResult=failed2Solve;return
elseif(.not.firstHook)then  ! internal trust iteration:
! check if computed results can be reused to avoid redundant factorizations.
! this option is not fully tested unsupported
  if(allowNewtReuse.and.stepResult==onTrustBound.and.newtLen<(one+sigma1)*trustRad)then
! * use available Newton step on this expansion since it is already known to be inside trust
    if(newtLen>zero)then
      stepResult=insideTrust;psol=newtStep;Pnorm=newtLen;nchol=0
      err=0;message="solveTrustHook/reusedNewtonStep"
      return
    endif
  elseif(allowNegReuse.and.stepResult==hardCase.and.&
    (negLen<(one+sigma1)*trustRad.or.negLen>sigma3reuseHard*trustRad))then
! * previous trust iteration was the "hard case" and now the trust has now
!   (a) increased. Simply upscale negative curvature since this step would have been
!       even "harder" and results are fairly predictable: a multiple of the negative
!       curvature eigenvector
!   (b) decreased not too much. Then downscale negative curvature direction since
!       the current direction is probably reasonable.
    if(negLen>zero)then
      stepResult=hardCase;psol=(trustRad/negLen)*negStep;Pnorm=trustRad;nchol=0
      err=0;message="solveTrustHook/reusedUpscaledNegStep"
      return
    else
! - this should never occur and suggests a bug, since the hard case on the previous
!   iteration must have used a scaled negative curvature step...
      err=200;message="f-solveTrustHook/BUG/(negLen<0)"
      return
    endif
  endif
endif
call flip_UtoL(B) ! ugly but currently necessary to compute Bnorm1 below
Hdiag=getdiag(B);Gnorm=norm2(grad);L_s=maxval(-Hdiag);Bnorm1=getKnorm(B,1)
L_lo=max(zero,L_s,Gnorm/trustRad-Bnorm1);Einf=hugeRe;stepResult=failed2Solve
L_hi=Gnorm/trustRad+Bnorm1; RZ2=zero; err=okAlg; message="solveTrustHook/ok"
do nchol = 1, ncholMax
! 1. Safeguard Lambda
! NB: inputting lambda=0 does not guarantee first solution will use lambda=0.
  lambda=max(lambda,L_lo); lambda=min(lambda,L_hi)
  if(lambda<=L_s)lambda=max(L_s_red*L_hi,(L_lo*L_hi**powLhi)**rootLhi)
! 2. Check positive definiteness of perturbed Hessian
  call addDiag(B,lambda)   ! augment Hessian diagonal
  if(usePivot)then ! - requests pivoted (robust) Cholesky: not really necessary here
    indx=-10  ! no prescrambling
    call choles_dcmp(a=B,Ld=Ld,posDefinite=pd,&
      doPivot=usePivot,indx=indx,skipRobust=.true.,&
      tau=cholTau,tauBar=cholTauBar,mu=cholMu,&
      ibad=ibad,logDet=logdetTemp,condest=condestTemp,&
      err=jerr,message=jmessage)
  else            ! - standard un-pivoted Cholesky (satisfactory in trust methods)
    indx=-20  ! to trap any possible bugs: indx should not be used.
    call choles_dcmp(a=B,Ld=Ld,posDefinite=pd,&
      ibad=ibad,logDet=logdetTemp,condest=condestTemp,&
      err=jerr,message=jmessage)
  endif
  call putDiag(B,Hdiag) ! restore Hessian diagonal
  if(nchol==1)then  ! keep original properties (not always useful, though)
    if(pd)then      ! - input Hessian positive definite
      logdet=logdetTemp;condest=condestTemp
    else            ! - input Hessian indefinite
      logdet=-hugeRe;condest=-hugeRe
    endif
  endif
  if(pd)then
! 3. Solve for restricted step direction using forward/backward substitution
    if(lambda<Einf)Einf=lambda
    q=-grad   ! solve Newton equations w/wo pivoting
    call choles_fwbw(A=B,Ld=Ld,usePivot=usePivot,indx=indx,b=q,x=psol,&
                     ynorm=RP2,err=jerr,message=jmessage)
    RP2=RP2**2      ! RP2 is ||y||^2, where y is L.y=q
! 4. Compute terms possibly needed for the upcoming Newton iteration
    Pnorm=norm2(psol)
    if(nchol==1.and.firstHook.and.lambda==zero)then
! store full Newton step (even if out of current trust bounds)
! if original B positive definite. This may be re-used later.
      newtStep=psol; newtLen=Pnorm
    else
      newtLen=min(-one,undefRN)
    endif
    if(Pnorm<(one-sigma1)*trustRad)then ! possibly the hard case: get eigenvector
      call choles_negEigVec(B=B,Ld=Ld,eigmeth=largeNorm,& ! of negative
        usePivot=usePivot,indx=indx,normW=normW,iBad=-10,& ! (or least positive)
        eigvec=zHat,eigval=eigval,RZ2=RZ2,err=jerr,message=jmessage) ! curvature
      if(jerr/=0)then ! this error shouldnt really occur
        err=failAlg;message="f-solveTrustHook/BUG?/&"//jmessage
        psol=hugeRe;pnorm=-hugeRe
        return
      endif ! compute scaling tau to construct final step with length trustRad
      call solveTrustHardCase(newtStep=psol,newtLen=pnorm,negStep=zHat,&
                              trustRad=trustRad,tau=tau)
    endif !NB: neg. curvature eigenvector zHat assumed to be normalized
  else
    newtLen=min(-one,undefRN) ! flag that Newton step undefined (B not +ve definite)
  endif
! 5. Update safeguards
  if(pd)then
    if(lambda>L_s.and.Pnorm<trustRad)then
      L_hi=min(L_hi,lambda)
    else
      L_lo=max(L_lo,lambda)
    endif
    L_s=max(L_s,lambda-RZ2)
  else
    L_lo=max(L_lo,lambda)
    L_s=L_lo ! DK: seems to better handle indefinite original B
  endif
  L_lo=max(L_lo,L_s)
! 6. Check convergence criteria
  if(pd)then  ! check condition 3.10 (using definition of trust region)
    if(abs(trustRad-Pnorm)<=sigma1*trustRad)then  ! constrained step on trust boundary
      stepResult=onTrustBound
      exit
    elseif(Pnorm<=trustRad.and.lambda==zero)then  ! full step inside trust region
      stepResult=insideTrust
      exit
    elseif(Pnorm<=trustRad)then     ! possibly the "hard case". Condition 3.11
      if(RZ2*tau**2<=sigma1*(two-sigma1)*max(sigma2,RP2+lambda*trustRad**2))then
        psol=psol+tau*zHat ! in some cases even non-hard cases are
        pnorm=norm2(psol)   ! terminated using this test, which checks that Newton
        stepResult=hardCase ! plus scaled least-positive curvature close to bound
        negStep=zHat;negLen=tau  ! - store eigenvector for possible reuse
        exit
      endif
    endif
  endif
! 7. Update lambda using safeguarded Newton iteration
  if(pd.and.Gnorm>tolGnorm)then  ! Newton update of lambda
! algorithm (3.2)-step 3: auxiliary vector
    if(usePivot)psol=psol(indx) ! account for pivoting scrambling
    call lower_rsolv(a=B,d=Ld,b=psol,x=q,transp=.false.,err=jerr)
    Qnorm=norm2(q)  ! Qnorm is residual for Newton-Hebden iteration
! algorithm (3.2)-step 4: Newton iteration on lambda (Hebden iteration)
    lambda=lambda+(Pnorm/Qnorm)**2*(Pnorm-trustRad)/trustRad
  elseif(pd.and.Pnorm>trustRad.and.L_s==zero.and.L_lo>=L_s.and.lambda==L_lo)then
    lambda=max(2._mrk*lambda,L_s_red*L_hi,sqrt(L_lo*L_hi))  ! prevent rare cycling
  else        ! safeguarded update when Newton iteration fails
    lambda=L_s
  endif
enddo
if(nchol>ncholmax)then  ! - failed to solve the trust region problem
  err=failAlg
  write(message,'(a,i0)')"solveTrustHook/nchol>ncholmax:",nchol
  psol=zero;Pnorm=zero;stepResult=failed2Solve
else
  err=okAlg; message="solveTrustHook/ok"
!  call putDiag(B,Hdiag) ! may want the original Hessian in the calling routine
endif
! End procedure here
endsubroutine solveTrustHook
!----------------------------------------------------
pure subroutine solveGeneralDogTrust(B,Ld,indx,grad,trustRad,dogNewtBias,&
  haveFac,haveNewt,haveNeg,doPivot,hessFacBundle,haveGBG,useL,&
  posDef,newtStep,newtLen,negEigen,negStep,&
  normG,gBg,absGdotNewt,&
  logdet,condest,Einf,psol,pLen,stepResult,nchol,err,message)
! ---
! Purpose: Implements the generalized dogleg trust region solution,
! handling Hessians with arbitrary convexity properties at the cost
! of a single (modified) Cholesky decomposition (except in the hard case).
! The method may not produce as accurate trust solutions as the near-exact
! approach of More and Sorensen, but can be much cheaper (for large problems) and
! not too shabby either.
! Recommended for large-scaled problems where memory is not an issue
! but where repeated Cholesky inversions are becoming onerous.
! The generalized dogleg is often referred to as the 2D subspace minimisation
! solution of the trust region subproblem. The implementation here is DK's
! concoction of Dennis and Schnabel, Nocedal and Wright, Schultz et al.,
! Gill et al. concepts, resulting in a simple 1-factorization algorithm
! that often would require 1 Cholesky per _outer_ trust region iteration.
! In addition, factored quasi-Newton updating allows using the dogleg
! with no Cholesky factorizations, ie, solving the trust subproblem
! in O(2) cost. not shabby... (but then cannot handle the hard-case).
! NB:
! * Dogleg requires single Cholesky except in the hard case where inverse
! iteration eigenvector polish requested, in which case 4-5 Cholesky factorizations
! are usually sufficient. NB: hard case only arises for indefinite functions
! near saddle-points, so BFGS and Gauss-Newton Hessians 'should' never invoke
! the hard case code. SR1 can sometimes be problematic since it can become strongly
! indefinite with very large negative diagonals. Indeed, my experimentation suggests
! SR1 quasi-Newton can be expensive in the hard case.
! * In addition, the use of robust Cholesky in the dogleg method can flag marginal
! positive-definite matrices as indefinite. This can occur for unfactored BFGS
! and Gauss-Newton Hessians. Code then uses its "hard-case" algorithm to proceed.
! ---
! Programmer: Dmitri Kavetski. 17 January 2004
! ---
! INPUT:
! B         = depending on doChol, either
!              (i)  raw scaled active Hessian or
!              (ii) permuted modified scaled active Hessian L factor (Ld=diag of L)
! grad      = gradient
! haveFac   = true if B already decomposed with L diagonal in Ld
!             (otherwise instructs to carry out (robust) Cholesky decomposition
! indx      = permutation vector
! dogNewtBias  = dogleg bias towards Newton,0=single dogleg,1=scaledNewton (~0.8)
! INPUT/OUTPUT (input if haveFac,output if .not.haveFac)
! posDef    = true if original B is positive definite
! newtStep  = Newton step (posdef) or modified Newton step (.not.posdef)
! cauchyStep= Cauchy step
! negEigen  = most negative eigenvalue of B
! negStep   = step of negative curvature (eigenvector)
! logdet,condest,Einf = Hessian properties
! OUTPUT
! psol      = solution of trust region problem
! pLen      = length of trust region solution
! err       = status completion
! message   = description of performance
! ---
! USAGE
! * Trust region optimization has two iteration loops
!   - Inner iteration ("trust acceptance") loop, where the function is
!     trialled along the trust region trajectory until sufficient decrease
!     obtained (not unlike a curvilinear search).
!   - Outer iteration ("step") loop, 
! * Normally, call this routine with all "have" vars set to false and
!   {normG,gDotNewt}<0 and the routine will calculate (and return) whatever is
!   needed for the trust region solution.
! * If have robust Cholesky factors (eg, factored quasi-Newton) then set
!  "haveFac=true" and supply the lower triangle of L and diagonal Ld.
!   Also need to supply posdef, if .false. then will assume it received
!   a modified Newton step and will solve the "hard case" if it is too short.
!   Note that pivoting complicates the use of this subroutine, so take care!
! * If lots of information is known apriori (eg, Newton steps), set
!   corresponding "have"'s to .true. and the routine will use the supplied data
!   (with no checks!, so be sure u no what u dooing...).
! * If the dogleg step failed to achieve sufficient decrease, call dogleg
!   again with decreased trust region but do not alter any "have" variables.
! ---
! Algorithm:
! * The generalized dogleg step is a subset of the 2D subspace minimization
!   strategy of solving the trust region subproblem.
!   - If the Hessian is positive definite and Newton step inside trust region,
!     simply take Newton step.
!   - Whenever the (possibly modified) Newton step exceeds the trust region, the
!     "exact" hookstep curve is replaced by piecewise linear intervals, connecting
!     point A = current point
!     point B = Cauchy point (constrained minimizer of linear model)
!     (point C) = for double dogleg step, point along BD which biases towards Newton.
!     point D = (modified) Newton point (unconstrained minimizer of quadratic model)
!   - If the modified Newton step is shorter than the trust radius, we are dealing
!     with the "hard case" and the modified Newton step is pumped all the way to
!     the trust bound along a direction of negative curvature. This allows the
!     trust region optimizer to escape from saddle points along directions of
!     negative curvature (ie, eigenvectors of negative eigenvalues).
!   - The "2D" bit in "2D subspace minimization" comes from the 2D subspace
!     obtained by joining the Cauchy and Newton steps. The actual trust
!     solution is a curve (hook), so the sucess of doglegs depends on
!     whether curvature is important and whether the computed Hessian
!     is reliable source of this information.
! * Dogleg method comes in two sub-flavours: single-dogleg and double-dogleg
!   the double-dogleg is biased towards the Newton step even if the latter is
!   outside the trust region, whereas the single-dogleg step (originally by Powell)
!   simply connects the Cauchy and Newton points. Set dogNewtBias=0 for single dogleg
!   or dogNewtBias=0.8 for standard double-dogleg (Dennis and Schnabel).
!   dogNewtBias=1 will give total bias to Newton step and will simply scale it to the
!   trust radius, discarding the Cauchy step, which seems a bit extreme and kind of
!   contrary to the spirit of trust regions (interpolating steepest descent and Newton).
! * Three options for the hard case, set by eigmeth
!   - fastChol: Fast O(N2) eigenvector estimation from the modified Cholesky factors using
!     the method of Gill et al. Method generally reliable, but requires the
!     reliable identification of offending rows of the Hessian. When pivoting
!     enabled, this is usually accurately determined.
!   - largeNorm: More accurate but more expensive O(N2) method based on large-norm method.
!   - Both these methods are approximate, require the single robust Cholesky but
!     can break down in some cases (since robust Cholesky usually overadds).
!   - invIter: uses inverse iteration to polish up eigenvectors, which often
!     saves otherwise ruined estimates. This may require additional O(N3) Choleskying,
!     which may be unavoidable to reliably compute eigenvectors in the hard case.
!   NB: The hook step code is particularly robust in the hard case.
! ---
! Refs:
! * More and Sorensen (1983) Computing a trust region step,
!   SIAM Journal of Scientific and Statistical Computing,4(3),pp.553-572.
! * Dennis and Schnabel (1996) Numerical methods for unconstrained
!   optimization and nonlinear equations. text and pseudocode.
! * Nocedal and Wright (1999) Numerical Optimization (dogleg chapters)
! * Schultz,G.A., Schnabel,R.B. and Byrd,R.H.(1985) A family of
!   trust-region-based algorithms for unconstrained minimization with
!   strong global convergence properties, SIAM Journal on Numerical
!   Analysis,V.22(1),Feb.1985,pp.47-67.
! ---
use utilities_dmsl_kit,only:zero,twoThirds,one,norm2,assertEq,arthsi,&
  getDiag,putDiag,addDiag,triang_minEig
use linalg_dmsl_kit,only:choles_dcmp,choles_fwbw,choles_negEigVec
implicit none
! dummies
logical(mlk),intent(inout)::haveFac,haveNewt,haveNeg,haveGBG
logical(mlk),intent(in)::doPivot,useL
type(hessFacBundle_type),intent(in)::hessFacBundle
real(mrk),intent(inout)::B(:,:),Ld(:)
real(mrk),intent(inout)::normG,gBg,absGdotNewt
real(mrk),intent(in)::grad(:),trustRad,dogNewtBias
integer(mik),intent(inout)::indx(:)
logical(mlk),intent(inout)::posDef
real(mrk),intent(inout)::newtStep(:),newtLen
real(mrk),intent(inout)::negStep(:),negEigen
real(mrk),intent(inout)::psol(:),pLen
real(mrk),intent(inout)::logdet,condest,Einf
integer(mik),intent(out)::stepResult,nchol,err
character(*),intent(out)::message
! locals
integer(mik)::n,job,cauchyStepType,iBad,ncholHC
real(mrk)::E(size(grad)),Gersh(size(grad)),cauchyStep(size(grad)),cauchyLen
real(mrk)::tempv(size(grad)),negeigChol,negeigHard
real(mrk)::tau,gamma,tempA,tempB,nu,eigZero
logical(mlk)::usePivot
! auxiliary
logical(mlk)::ok
integer(mik)::jerr
character(100)::jmsg
! algorithm parameters
real(mrk),parameter::sigma1=0.2_mrk     ! tolerance on step and trust agreement
real(mrk),parameter::eigtol=1.e-1_mrk   ! tolerance on negative curvature eigenvalue
integer(mik),parameter::ncholmaxHC=100  ! max Cholesky factorizations in the hard case
integer(mik),parameter::useNewt=0,doDog=1,hardCase=2
integer(mik),parameter::fastChol=0,largeNorm=1,inviter=2
integer(mik),parameter::eigmeth=invIter ! largeNorm ! fastChol ! ! hard-case eigenmethod
!integer(mik)::eigmeth ! hard-case eigenmethod
integer(mik)::cmeth ! eigenmethod Cholesky option
!logical(mlk),parameter::checkEig=.false.  ! forces Cauchy step for ill-conditioned Hessians
logical(mlk)::checkEig  ! forces Cauchy step for ill-conditioned Hessians
logical(mlk),parameter::normW=.true.
! Start procedure here
checkEig=.false.
call assertEq(size(B,1),size(B,2),size(grad),size(indx),size(psol),&
              size(newtStep),size(cauchyStep),size(negStep),ok,n)
if(.not.ok)then
  err=100;message="f-solveGeneralDogTrust/dimError"
  return
endif
! * Process Hessian matrix
if(haveFac)then
! - Assumes B already decomposed and uses {Ld,indx,posDef,negEigen,negCurv}
  usePivot=doPivot
  nchol=0;err=0;message="solveGeneralDogTrust/usingInputLfac"
else
! - Robust Cholesky decomposition of B to establish whether it is positive
! definite and perturb if not (estimating negative eigenvalue)
  selectcase(hessFacBundle%facmeth)   ! select modified Hessian factorization method
  case(schnab_facmeth)  ! - revised modified Cholesky-Gershgorin of Schnabel and Eskew
    usePivot=doPivot; indx=-1 !  (indicate no pre-scrambling)
    call choles_dcmp(A=B,Ld=Ld,iBad=iBad,&
      tau=hessFacBundle%tau,tauBar=hessFacBundle%tauBar,mu=hessFacBundle%mu,&
      doPivot=usePivot,indx=indx,posDefinite=posDef,logDet=logdet,condest=condEst,&
      Einf=negeigChol,E=E,Gout=Gersh,err=jerr,message=jmsg)
    nchol=1 ! note if pivoting used then everything will be scrambled (and scaled)
  case(dennis_facmeth)  ! - perturbed Cholesky-Gershgorin of Dennis and Schnabel
    usePivot=.false.; indx=-100 ! no pivoting here
    call choles_dcmp(A=B,Ld=Ld,&
      maxCond=hessFacBundle%maxHessCond,&
      posDefinite=posDef,nchol=nchol,logDet=logdet,condest=condEst,&
      Einf=Einf,err=jerr,message=jmsg)
  endselect
  haveFac=.true.
  Einf=negeigChol ! store estimated most negative eigenvalue
!  call addDiag(B,E) ! explicitly construct modified matrix in upper triangle
endif
! * Compute generalized Newton (if(posDef)=>classic Newton, else=modified)
if(.not.haveNewt)then   ! perform Cholesky forward/backward substitution
  if(.not.usePivot)then
    indx=arthsi(n)   ! assume unpivoted solution
  elseif(any(indx<1.or.indx>n))then ! this catches evident errors, but is not
    err=10;message="f-solveGeneralDogTrust/indxContentError"  ! bombproof...
    return
  endif
  call choles_fwbw(a=B,Ld=Ld,indx=indx,usePivot=usePivot,&
                   b=grad,x=newtStep,err=jerr,message=jmsg)
  newtStep=-newtStep; haveNewt=.true.
  newtLen=norm2(newtStep) ! length of (modified) Newton step
endif
! * Decide whether to (i) dogleg or (ii) use negative curvature
if(posDef)then
! - original Hessian positive definite: use dogleg step if Newton too big.
  job=merge(doDog,useNewt,newtLen>(one+sigma1)*trustRad)
elseif(newtLen>(one-sigma1)*trustRad)then
! - modified Newton at least long enough, use dogleg if too big.
  job=merge(doDog,useNewt,newtLen>(one+sigma1)*trustRad)
else
! - original Hessian indefinite and we are faced with the "hard case".
  job=hardCase  ! (since modified step too short)
endif
! * Carry out requested trust job: either dogleg or negative curvature
selectcase(job)
case(useNewt)   ! - simply return Newton step as it is inside trust region
  psol=newtStep; pLen=newtLen; stepResult=insideTrust
  err=0;  message="f-solveGeneralDogTrust/ok/usedNewtonInsideTrust"
case(doDog)     ! - standard dogleg step if Newton or modified Newton too long
! - Compute Cauchy point  ! ** this code is copied to below
  if(useL)then                    ! - put Cholesky diagonal into B for a sec...
    E=getDiag(B); call putDiag(B,Ld)
  endif
  if(useL.and.usePivot)then       ! - account for pivoting (for gBg computation)
    tempv=grad(indx)
  else
    tempv=grad
  endif
! (NB: if normG,gBg known then this cheap call merely scales the Cauchy step)
  call getCauchyStep(hess=B,useL=useL,grad=tempv,trustRad=trustRad,&
                     normG=normG,haveGBG=haveGBG,gBg=gBg,&
                     cauchyStepType=cauchyStepType,&
                     cauchyStep=cauchyStep,cauchyLen=cauchyLen)
  if(useL)call putDiag(B,E)                         ! - restore diagonal of B
  if(useL.and.usePivot)cauchyStep(indx)=cauchyStep  ! - unscramble Cauchy
  if(dogNewtBias>zero)then  ! double dogleg step (biased towards Newton)
    if(absGdotNewt<zero)absGdotNewt=abs(dot_product(grad,newtStep))
    gamma=normG**4/(gBg*absGdotNewt)
    nu=(one-dogNewtBias)+dogNewtBias*gamma
  else              ! single dogleg step
    nu=one
  endif
  if(nu*newtLen<=trustRad)then
! - in interval CD: use scaled Newton step
    psol=(trustRad/newtLen)*newtStep; pLen=trustRad
    err=0;  message="solveGeneralDogTrust/ok/usedDoglegScaledNewton(CD)"
  elseif(cauchyStepType/=cauchyInside)then
! - scaled Cauchy step if hitting trust bound
    psol=cauchyStep; pLen=cauchyLen
    err=0;  message="solveGeneralDogTrust/ok/usedDoglegCauchy(CauchOnTrust)"
  else
! - in interval BC: linear segment connecting Cauchy and Newton.
! note that in this branch Cauchy is always inside trust but
! Newton is always outside trust region. The solution is analytical
! and is obtained by solving a quadratic eqns (picking positive root)
    psol=nu*newtStep-cauchyStep
    tempA=dot_product(psol,psol); tempB=dot_product(psol,cauchyStep)
    tau=(cauchyLen-trustRad)*(cauchyLen+trustRad)
    tau=(-tempB+sqrt(tempB**2-tempA*tau))/tempA
    psol=cauchyStep+tau*psol; pLen=trustRad
    err=0;  message="solveGeneralDogTrust/ok/usedDoglegPath(BC)"
  endif
  stepResult=onTrustBound
case(hardCase)  ! - hard case needs pumping up with negative curvature:
! need to pump an exceedlingly short modified Newton step with a good multiple
! of negative curvature. This occurs near saddle points where a small gradient
! suppresses natural Newton steps (even for modified Hessians).
  if(.not.haveNeg)then  ! - compute eigenvector afresh
    iBad=min(ibad,n)    ! ... usually but not always ibad~n, particularly with pivoting
!    eigmeth=merge(invIter,largeNorm,hessFacBundle%facmeth==schnab_facmeth)
    selectcase(hessFacBundle%facmeth)
    case(schnab_facmeth)
!      cmeth=merge(1,2,usePivot) ! single Schnabel Cholesky with pivoting good enough
      cmeth=2 ! always multiple Cholesky to polish eigenvalues
    case default                ! need more expensive polish
      cmeth=2
    endselect
! - force inverse iteration polish when revised robust Cholesky used, since
! its diagonal perturbation E is not constant and hence does not just shift
! eigenvalues, but also modifies eigenvectors.
    selectcase(eigmeth)
    case(fastChol,largeNorm)
! - use modified Cholesky factors to estimate negative curvature
      call choles_negEigVec(B=B,Ld=Ld,indx=indx,eigmeth=eigmeth,&
        usePivot=usePivot,iBad=iBad,normW=normW,&
        eigvec=negStep,eigval=negeigHard,err=jerr,message=jmsg)
    case(invIter)
! - use inverse iteration polish to ensure more reliable results
      negeigChol=-Einf
      call choles_negEigVec(A=B,doRob=.false.,Lrob=B,Ld=Ld,eigRob=negeigChol,&
        doPivot=usePivot,cmeth=cmeth,indx=indx,ibad=ibad,&
        nitermax=ncholmaxHC,Etol=eigtol,Evtol=eigtol**2,&
        posDef=ok,eigval=negeigHard,eigvec=negStep,&
        nchol=ncholHC,err=jerr,message=jmsg)
      nchol=nchol+ncholHC ! - requires additional Cholesky decompositions
    endselect
    if(jerr/=0)then
      err=20;message="f-solveGeneralDogTrust/bug?badSet?/&"//jmsg
      return
    endif
    haveNeg=.true.; negEigen=negeigHard
  endif
  eigZero=triang_minEig(Ld=Ld,condMax=one/epsRe**twoThirds,cholLd=.true.)
  if(checkEig.and.negEigen>-eigZero)then
! - pathological positive semi-definite (near-singular) Hessian
! DK's experimentation suggests taking Cauchy step for ill-conditioned Hessians
! is not always efficient, particularly with SR1 quasi-Hessians.
! take Cauchy step---
    if(useL)then  ! ** this code is copied from Cauchy above
      E=getDiag(B); call putDiag(B,Ld)
    endif
    if(useL.and.usePivot)then
      tempv=grad(indx)
    else
      tempv=grad
    endif
    call getCauchyStep(hess=B,useL=useL,grad=grad,trustRad=trustRad,&
                       normG=normG,haveGBG=haveGBG,gBg=gBg,&
                       cauchyStepType=cauchyStepType,&
                       cauchyStep=cauchyStep,cauchyLen=cauchyLen)
    if(useL)call putDiag(B,E)
    if(useL.and.usePivot)cauchyStep(indx)=cauchyStep
    stepResult=merge(insideTrust,onTrustBound,cauchyStepType==cauchyInside)
    psol=cauchyStep; pLen=cauchyLen
! endtake Cauchy step---
    err=0;  message="w-solveGeneralDogTrust/ok/usedCauchy(pathosB)"
  else                  ! - yep.. hard case
! employs More and Sorensen algorithm for scaling the negative curvature
! eigenvector up to the trust bound, using the solution of the corresponding
! quadratic to preserve as much Newton direction as possible.
! This seems efficient even for near-singular Hessian matrices, since in this
! case the step comprises the direction where the function curves-up least.
    call solveTrustHardCase(newtStep,newtLen,negStep,trustRad,psol=psol)
    stepResult=hardCase
    err=0;  message="f-solveGeneralDogTrust/ok/usedHardCase"
  endif
case default
  err=200;message="f-solveGeneralDogTrust/BUG/unknownJob"
endselect
! End procedure here
endsubroutine solveGeneralDogTrust
!----------------------------------------------------
pure subroutine getCauchyStep(hess,useL,grad,trustRad,normG,haveGBG,gBg,&
  cauchyStep,cauchyLen,cauchyStepType)
! Purpose: Compute scaled Cauchy step to the minimizer of the linear model.
! All data is assumed to be prescaled (for efficiency).
! INPUT
!   hess            = scaled (possibly modified) Hessian.
!                     not used if haveGBG=true, in which case user supplies
!                     gBg = grad(tranpose).dot.hess.dot.grad
!   grad            = scaled gradient. Not used if normG=||grad||>0
!   trustRad        = scaled trust radius
!   useL            = instructs to use Cholesky L in lowerTriangle of hess
! OUTPUT
!   cauchyStep      = step to Cauchy point CP
!   cauchyLen       = length of Cauchy step
!   cauchyStepType  = type of Cauchy step (eg, if reaches trust boundary, etc.)
!---
! See eqn 4.7-4.8 in Nocedal.
! Programmer: Dmitri Kavetski, 17 January 2004.
use utilities_dmsl_kit,only:zero,one,norm2,quadform,fmatmul_mv
implicit none
! dummies
real(mrk),intent(in)::hess(:,:),grad(:),trustRad
logical(mlk),intent(in)::useL,haveGBG
real(mrk),intent(inout)::normG,gBg
real(mrk),intent(out)::cauchyStep(:),cauchyLen
integer(mik),intent(out)::cauchyStepType
! locals
real(mrk)::tau
real(mrk),parameter::normGminFac=1.e3_mrk,normGmin=tinyRe*normGminFac 
! Start procedure here
if(.not.haveGBG)then  ! compute {grad(t).B.grad}
  if(useL)then  ! - use (possibly modified) Cholesky factor L
    cauchyStep=fmatmul_mv(m=hess,v=grad,typeMV="LTV")
    gBg=norm2(cauchyStep)**2
  else          ! - use (possibly modified) Hessian matrix
    gBg=quadform(v=grad,mm=hess,typeM="SU")
  endif
endif
if(normG<zero)normG=norm2(grad) ! length of steepest descent step
if(gBg<=zero)then ! - linear model decreases monotonically ('infinite' Cauchy step)
  cauchyLen=trustRad;       cauchyStepType=cauchyInfin
else              ! - convex quadratic in tau, possibly constrained by trust bound
  tau=normG**3/(trustRad*gBg)
  if(tau>=one)then  ! -- yep, constrained
    cauchyLen=trustRad;     cauchyStepType=cauchyOnBound
  else              ! -- inside trust
    cauchyLen=tau*trustRad; cauchyStepType=cauchyInside
  endif
endif
if(normG>normGmin)then  ! safeguard division by zero
  cauchyStep=-cauchyLen*grad/normG
else
  cauchyStep=-grad/normGmin; cauchyStepType=cauchyZeroGrad
endif
! End procedure here
endsubroutine getCauchyStep
!----------------------------------------------------
pure subroutine solveTrustHardCase(newtStep,newtLen,negStep,trustRad,psol,tau)
! Purpose: Given modified Newton step and a direction of negative curvature,
! construct a globally convergent trust region step using a weighted
! combination of the modified Newton step and the negative curvature direction.
! This allows trust region methods to escape saddle point regions of attraction.
! * More and Sorensen (1983) Computing a trust region step,
!   SIAM Journal of Scientific and Statistical Computing,4(3),553-572.
implicit none
! dummies
real(mrk),intent(in)::newtStep(:),negStep(:),newtLen,trustRad
real(mrk),intent(inout),optional::tau,psol(:)
! locals
real(mrk)::dotPZ,tau0
! Start procedure here
dotPZ=dot_product(newtStep,negStep) ! multiple of eigenvector is calculated so that
tau0=(trustRad-newtLen)*(trustRad+newtLen) ! it does not cancel out the Newton
tau0=tau0/(dotPZ+sign(sqrt(dotPZ**2+tau0),dotPZ))  ! component of the step
if(present(tau))tau=tau0  ! (and thus keeps at least some Newton...)
if(present(psol))psol=newtStep+tau0*negStep
! End procedure here
endsubroutine solveTrustHardCase
!----------------------------------------------------
subroutine updateTrust(evalFunc,dataIN,dataOUT,xold,fold,gold,dx,stepLen,stepResult,redExp,&
  xscale,fscale,reducedTrust,objFuncBundle,trustBundle,&
  x,fx,trustRad,redRatio,fcalls,retcode,message)
! Purpose: Given a proposal step dx obtained by solving the trust region
! subproblem, accept or reject step dx and update the trust region for the
! next iteration or step.
!
! Input:  current point xold with function value fold
!         expected reduction in function value
!         trial (constrained) step dx of length stepLen
! Output: updated point x with function value fx
!         updated trust region and reduction ratio
!         status diagnostix
! Method: A modification of the basic trust update algorithm described by
!         Nocedal and Wright (1999) and Fletcher (1996).
!         Also includes Dennis and Schnabel details.
!
! Programmer: Dmitri Kavetski
! Algorithm:
! * The method assumes the input step has been truncated accounting for any bounds.
!   It is a bit awkward to handle elliptical trust regions near rectangular bounds
!   (see,eg.Nocedal). In these cases a box-step method may be more appropriate.
!   Note that reducing the trust region to the distance to nearest bound is
!   generally unsatisfactory since it will often unnecesarily truncate steps _away_
!   from 'em...
! * Care must be taken since the step may have been truncated by solution bounds
!   In this case the trust region may be accurate but truncated step may be far shorted.
!   In this case do not alter trust region.
! * If quadratic model is good but the trust region constrained the step,
!   the step is re-attempted with larger trust radius. This may achieve a greater
!   function reduction without additional gradient/Hessian calls, which becomes
!   particularly beneficial whenever the dimensionality of the objective function is
!   high and if derivatives are approximated by finite differences.
! * If trust region was close to constraining the step and the quadratic model
!   was good, increase trust region for next step.
! * If quadratic model poor, reduce trust region (but usually accept step)
! * Alternative strategies may include a (curvilinear linesearch in the trust region
!   direction, this would save those extra Cholesky decompositions when expanding
!   trust region.
! * For SR1 quasi-Newton Hessian updating, it may be preferable to update the
!   Hessian even along failed directions, in order to incorporate as much curvature
!   information as possible into the quadratic model.
use types_dmsl_kit,only:data_ricz_type
use utilities_dmsl_kit,only:one,zero,half,minmax
implicit none
! dummies
type(data_ricz_type),intent(in),optional::dataIN
type(data_ricz_type),intent(inout),optional::dataOUT
real(mrk),intent(in)::xold(:),gold(:),fold,dx(:),stepLen,redExp
integer(mik),intent(in)::stepResult
type(trustBundle_type),intent(in)::trustBundle
type(objFuncBundle_type),intent(in)::objFuncBundle
logical(mlk),intent(in)::reducedTrust
real(mrk),intent(in)::xscale(:),fscale
real(mrk),intent(out)::x(:),fx,redRatio
real(mrk),intent(inout)::trustRad
integer(mik),intent(out)::fcalls
integer(mik),intent(out)::retcode
character(*),intent(out)::message
! user-provided function
interface
  subroutine evalFunc(dataIN,dataOUT,x,feas,fx,gradFx,hessFx,err,message)
  use kinds_dmsl_kit
  use types_dmsl_kit,only:data_ricz_type
  implicit none
  type(data_ricz_type),intent(in),optional::dataIN
  type(data_ricz_type),intent(inout),optional::dataOUT
  real(mrk),intent(in)::x(:)
  logical(mlk),intent(out)::feas
  real(mrk),intent(out),optional::fx,gradFx(:),hessFx(:,:)
  integer(mik),intent(out)::err
  character(*),intent(out)::message
  endsubroutine evalfunc
endinterface
! locals
integer(mik)::err
logical(mlk)::feas
real(mrk)::redObs,dirDer
! local parameters
logical(mlk),parameter::backtrackTrust=.false.  ! trust reduction using backtracking
real(mrk),parameter::radDownMin=0.1_mrk,radDownMax=0.5_mrk  ! safeguarded by [min,max]
real(mrk),parameter::safeEps=10._mrk*epsRe
logical(mlk),parameter::checkExpObs=.false.
! Start procedure here
retcode=failed_tr;fcalls=0;x=xold;fx=fold;redRatio=zero
if(all(abs(dx)<=safeEps*max(abs(xold),xscale)))then  ! step a bit too small
  retcode=dxTiny_tr; message="updateTrust/dxTiny"
  return
elseif(redExp<=zero)then  ! expected reduction non-positive: return now
  redRatio=-one; retcode=expRedNonP_tr
  write(message,'(a,sp,es13.6,s)')"f-updateTrust/expFredNonPosit(BUG?):",redExp
  return
endif
x=xold+dx
! ** evaluate function at trial point
call evalFunc(dataIN,dataOUT,x,feas,fx,err=err,message=message); fcalls=fcalls+1
if(err/=0)then        ! strange, since point x has already been trialled with err=0
  message="f-updateTrust/userErr1/&"//message
  retcode=badFunc_glob; return
elseif(.not.feas)then         ! unfeasible point encountered
  retcode=unfeas_tr; message="updateTrust/unfeasX"
  trustRad=trustBundle%radDown_tr*stepLen
  if(trustRad<trustBundle%trustMin)then   ! check for collapsed trust
    retcode=unfeas_glob; message="updateTrust/unfeas->trust->0"
  endif
  return
endif
! ** calculate observed reductions in the function
redObs=fold-fx; redRatio=redObs/redExp
if(checkExpObs.and.fx<fold.and.&  ! small reduction around noise levels
   abs(redObs)<=objFuncBundle%epsF*max(abs(fold),fscale).and. &
   redExp     <=objFuncBundle%epsF*max(abs(fold),fscale))then
! converged to function precision
  if(fx>fold)then ! return initial point if (marginally) better
    x=xold;fx=fold;redRatio=zero
  endif
  retcode=fconExpObs_tr;message="updateTrust/red[f](obs&exp)~epsF";return
endif
! Adjust trust radius if necessary
if(redRatio<trustBundle%acceptRatio_tr)then
! ** Step does not satisfy sufficient decrease condition: reject step
!  x(:)=xold(:);fx=fold ! point and function may be needed by callee
  trustRad=trustBundle%radDown_tr*stepLen
  retcode=failed_tr; message="updateTrust/fredTooLow/trustGoingDown"
  if(trustRad<=trustBundle%trustMin)then      ! * check for collapsed trust
    retcode=collapsed_tr;message="updateTrust/collapsedTrust"
  endif
elseif(redRatio<trustBundle%roDown_tr)then
!     .and.redExp>safeEps*max(abs(fold),fscale))then
! ** Step satisfies sufficient decrease condition but agreement with quadratic model
!    is poor. Accept step but deflate trust region for next step
  if(backTrackTrust)then  ! - use linesearch back-tracking to controllably reduce trust
    dirDer=dot_product(gold,dx)
    trustRad=-half*dirDer/(fx-fold-dirDer)
    trustRad=minmax(trustRad,radDownMin,radDownMax)
    trustRad=trustRad*stepLen
  else                    ! - simple reduction (tends to work better - fewer assumptions)
    trustRad=trustBundle%radDown_tr*stepLen
  endif
  retcode=suceed_tr; message="updateTrust/ok/trustGoingDown"
elseif(redRatio>trustBundle%roUpNow_tr.and.(    &
!       stepLen>trustRad*trustBundle%radUp_tr.or.&
       stepResult==onTrustBound.or.             &
       stepResult==hardCase).and.               &
       .not.reducedTrust)then
! ** Step in good agreement with quadratic model and was near-constrained by the
!    trust bound (likely unneccesarily). Re-attempt step with larger trust
  trustRad=trustBundle%radUp_tr*trustRad
  retcode=goBig_tr; message="updateTrust/goBig/trustGoingUp"
  if(trustRad>trustBundle%trustMax)then    ! * check for blown trust
    trustRad=trustBundle%trustMax; retcode=blown_tr
    message="updateTrust/trustRegionWantsBig"
  endif
elseif(redRatio>trustBundle%roUp_tr.and.&
       stepLen>trustBundle%stepOtrustUp_tr*trustRad.and.&
      .not.reducedTrust)then
! ** Step in good agreement with quadratic model and close to trust radius
!    Avoid possible future interference by pumping the trust up.
  trustRad=trustBundle%radUp_tr*trustRad
  retcode=suceed_tr; message="updateTrust/ok/trustGoingUp"
  if(trustRad>trustBundle%trustMax)then    ! * check for blown trust
    trustRad=trustBundle%trustMax; retcode=blown_tr
    message="updateTrust/trustRegionWantsBig"
  endif
else
! ** (i)  sufficient (but not great) decrease achieved or
!    (ii) good agreement with quadratic model but trust region non-interfering
! accept step but do not alter trust radius
  retcode=suceed_tr; message="updateTrust/ok/keepTrust"
endif
!print *, 'in updateTrust, fold    = ', fold
!print *, 'in updateTrust, xold    = ', xold
!print *, 'in updateTrust, dx      = ', dx
!print *, 'in updateTrust, xold+dx = ', xold+dx
!print *, 'in updateTrust, fx      = ', fx
!print *, 'in updateTrust, x       = ', x
! End procedure here
endsubroutine updateTrust
!----------------------------------------------------
endmodule optimiser_dmsl_kit
!******************************************************************
