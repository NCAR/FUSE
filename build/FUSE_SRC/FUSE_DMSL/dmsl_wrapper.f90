MODULE DMSL_WRAPPER_MODULE
USE kinds_dmsl_kit
IMPLICIT NONE
PRIVATE
PUBLIC::QNEWTON_WRAPPER,MCMC_WRAPPER,OBJFUNC_WRAPPER_OPTI,OBJFUNC_WRAPPER_MCMC
CONTAINS
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! (A) QUASI-NEWTON OPTIMIZATION
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
SUBROUTINE QNEWTON_WRAPPER(X0I,XLO,XHI,XSCALE,FDIGITS,UOUT, &     ! input
                           XOPT,FOPT,ITER,FCALLS,GCALLS,HCALLS, & ! output
                           IERR,MESSAGE)                          ! error handling
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Wrapper for the DMSL quasi-Newton optimization
! ---------------------------------------------------------------------------------------
USE nrtype                                                ! data types
USE multistats, ONLY:MSTATS                               ! provide access to error message
USE optimiser_dmsl_kit, ONLY:QNEWTON                      ! provide access to qnewton
IMPLICIT NONE
! ---------------------------------------------------------------------------------------
! (1) DUMMIES
! ---------------------------------------------------------------------------------------
! input
REAL(SP),DIMENSION(:),INTENT(IN)        :: X0I             ! initial estimate of solution
REAL(SP),DIMENSION(:),INTENT(IN)        :: XLO             ! lower bound on solution, either none or both bounds must be present
REAL(SP),DIMENSION(:),INTENT(IN)        :: XHI             ! upper bound on solution, either none or both bounds must be present
REAL(SP),DIMENSION(:),INTENT(IN)        :: XSCALE          ! typical scale of parameters
INTEGER(I4B),INTENT(IN)                 :: FDIGITS         ! number of reliable digits in function evaluation
!                                                          ! (-2=estimate,-1=full machine precision)
INTEGER(I4B),INTENT(IN)                 :: UOUT            ! output unit for run-time information 
! output
REAL(SP),DIMENSION(:),INTENT(OUT)       :: XOPT            ! optimum value of "x", for which f(x) takes its minimum value
REAL(SP),INTENT(OUT)                    :: FOPT            ! function value at optimum
INTEGER(I4B),INTENT(OUT)                :: ITER            ! number of steps (iterations)
INTEGER(I4B),INTENT(OUT)                :: FCALLS          ! number of function calls
INTEGER(I4B),INTENT(OUT)                :: GCALLS          ! number of gradient calls
INTEGER(I4B),INTENT(OUT)                :: HCALLS          ! number of Hessian calls
! error handling
INTEGER(I4B),INTENT(OUT)                :: IERR            ! error code
CHARACTER(*),INTENT(OUT)                :: MESSAGE         ! error message
! ---------------------------------------------------------------------------------------
! (2) LOCALS
! ---------------------------------------------------------------------------------------
! Active set (to identify parameters on bounds)
INTEGER(I4B),DIMENSION(SIZE(X0I))       :: ACTIVESET       ! active set (-1=lo,0=free,+1=hi), must be present if using xLo and xHi
! Define termination tolerances
REAL(SP)                                :: EPSF            ! desired precision
REAL(SP)                                :: GTOL            ! scaled gradient tolerance
REAL(SP)                                :: STOL            ! scaled step tolerance
REAL(SP)                                :: FTOL            ! scaled function tolerance
! Define scaling settings
REAL(SP),PARAMETER                      :: FSCALE=1._SP    ! scale of function
REAL(SP)                                :: STPMAX          ! maximum scaled stepsize/trust radius (set<0 for default)
! Define computational algorithms used in qnewton
INTEGER(I4B),PARAMETER                  :: IMETH=5         ! iteration globalisation method; 5=Near-exact trust method ("hookstep")
INTEGER(I4B),PARAMETER                  :: GMETH=1         ! gradient evaluation method; 1=Forward difference gradient
INTEGER(I4B),PARAMETER                  :: HMETH=6         ! Hessian evaluation method; 6=BFGS update of unfactored Hessian
! Define initialization settings
INTEGER(I4B),PARAMETER                  :: HIMETH=5        ! Diagonal of estimated d2f/dx2
REAL(SP)                                :: TRUSTRAD        ! initial scaled trust region radius (set<0 for internal default)
! Define maximum effort expended before termination
INTEGER(I4B),PARAMETER                  :: MAXITER=5000   ! Maximum number of iterations
INTEGER(I4B),PARAMETER                  :: MAXFEV=500     ! Maximum number of function calls
! Useful diagnostics and information
REAL(SP),DIMENSION(SIZE(X0I))           :: GRADOPT         ! gradient at the optimum
REAL(SP),DIMENSION(SIZE(X0I),SIZE(X0I)) :: HESSOPT         ! Hessian at optimum
! Memory footprint
REAL(SP)                                :: MEMHESS2        ! additional memory necessary for allocating internal Hessian storage
! Return codes and runtime messages
INTEGER(I4B)                            :: ERR_QN          ! error diagnostic, err=0->ok,<0=warning,>0=error
CHARACTER(LEN=256)                      :: MESSAGE_QN      ! status description
INTEGER(I4B)                            :: ILEN            ! length of error message
INTEGER(I4B)                            :: I               ! looping variable
! ---------------------------------------------------------------------------------------
! initialize variables
ACTIVESET(:) =  0      ! define active set (-1=lo,0=free,+1=hi)
TRUSTRAD     = -1._SP  ! use internal default for trust region radius
STPMAX       = -1._SP  ! use internal default for maximum scaled stepsize/trust radius
MSTATS%ERR_MESSAGE(1:31)='searching for the local optimum'
FORALL(I=32:LEN(MSTATS%ERR_MESSAGE)) MSTATS%ERR_MESSAGE(I:I)=' '
! define termination tolerances
EPSF = 10._SP**(-FDIGITS)  ! desired precision
GTOL = SQRT(EPSF)          ! scaled gradient tolerance
STOL = EPSF                ! scaled step tolerance
FTOL = EPSF                ! scaled function tolerance
! find local optimum in the vicinity of the starting point
CALL QNEWTON(OBJFUNC_WRAPPER_OPTI,                                      & ! Objective function to be minimised
             x0=x0i,                                                    & ! Initial estimate of optimum
             xLo=xlo,xHi=xhi,activeSet=activeSet,                       & ! Upper and lower bounds on solution, active set
             gtol=gtol,stol=stol,ftol=ftol,                             & ! Termination tolerances
             xscale=xscale,fscale=fscale,fdigits=fdigits,stpmax=stpmax, & ! Scaling settings
             imeth=imeth,gmeth=gmeth,hmeth=hmeth,                       & ! Computational algorithms
             himeth=himeth,trustRad=trustRad,                           & ! Initialisation settings
             maxIter=maxIter,maxFev=maxFev,                             & ! Termination due to excessive effort
             uout=uout,                                                 & ! Output unit for runtime information
             xopt=xopt,fopt=fopt,                                       & ! Approximated optimal solution
             gradOpt=gradOpt,hessOpt=hessOpt,                           & ! Useful diagnostics and information
             iter=iter,fcalls=fcalls,gcalls=gcalls,hcalls=hcalls,       & ! Computational cost report
             memHess2=memHess2,                                         & ! Memory footprint
             err=err_qn,message=message_qn)                               ! Return codes and runtime messages
! save errors
MSTATS%ERR_MESSAGE = MESSAGE_QN
IERR=ERR_QN; MESSAGE=MESSAGE_QN
!WRITE(*,'(4(I6,1X),20(F9.3,1X))') ITER,FCALLS,GCALLS,HCALLS,FOPT,XOPT
END SUBROUTINE QNEWTON_WRAPPER
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
SUBROUTINE MCMC_WRAPPER(sample0,sdevDiag0,ierr,message)                   ! initial values for samples
! ---------------------------------------------------------------------------------------
! Creators:
! ---------
! Martyn Clark and Dmitri Kavetski, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Wrapper for the DMSL MCMC routines
! ---------------------------------------------------------------------------------------
USE mcmc_dmsl_kit,ONLY:mbrSettings_type,mbrOut_type,metropolis_RK  ! MCMC data types
USE model_defn, ONLY:FNAME_PREFIX,FNAME_TEMPRY             ! prefix for filenames
USE multiparam, ONLY:LPARAM,NUMPAR                         ! list of model parameters
IMPLICIT NONE
! input
REAL(mrk),DIMENSION(:),INTENT(IN)       :: sample0         ! initial sample
REAL(mrk),DIMENSION(:),INTENT(IN)       :: sdevDiag0       ! initial diagonal of the covariance matrix
! output
integer(mik)                            :: ierr            ! error code
character(*)                            :: message         ! error message
! local
integer(mik), parameter                 :: text_len=256    ! string length
integer(mik)                            :: ipar            ! loop through model parameters
character(len=text_len),dimension(:),allocatable :: parNames ! parameter names
type(mbrSettings_type)                  :: mbrSettings     ! Algorithmic control parameters
type(mbrOut_type)                       :: mbrOut          ! Performance diagnostix
character(len=text_len)                 :: lineFmtIn       ! user-specified formatting for output
character(len=text_len)                 :: lineFmtOut      ! actual format used
! ---------------------------------------------------------------------------------------
! initialize errors
ierr=0; message='start of mcmc_wrapper, everything is a-ok'
! populate parameter names
allocate(parNames(0:NUMPAR), stat=ierr)
if (ierr.ne.0) then; message='mcmc_wrapper: problem allocating parNames'; stop; endif
parNames(0) = 'Variance'
DO IPAR=1,NUMPAR
 parNames(ipar) = LPARAM(IPAR)%PARNAME
END DO
! set filenames in mbrSettings
mbrSettings%samfiles = TRIM(FNAME_PREFIX)//'__'//mbrSettings%samfiles
! open up run time diagnostix
open(mbrSettings%uInfo,file=TRIM(FNAME_PREFIX)//'__'//'mcmc_info.txt',status='unknown')
CALL metropolis_RK(OBJFUNC_WRAPPER_MCMC,                                & ! Objective function to be minimised
                   title="FUSE MCMC",                                   &
                   varNames=parnames,                                   & ! Parameter names
                   mbrSettings=mbrSettings,                             & ! Algorithmic control parameters
                   sample0=sample0,sdevDiag0=sdevDiag0,                 & ! Initial values for samples
                   mbrOut=mbrOut,                                       & ! Performance diagnostix
!                    lineFmtIn=lineFmtIn,&!lineFmtOut,&                   & ! Formatting for the sample
                   err=ierr,message=message)                               ! Return codes and runtime messages
! Phase 1 - OnerPerTime - Use one-variable-per-time Metropolis
! Phase 2 - Scaling     - Compute covariance and adjust its scale using single-block Metropolis
! Phase 3 - Burnin      - Sample using fixed covariance matrix
! Phase 4 - Production  - Production samples

if (ierr.ne.0) message='mcmc_wrapper: '//trim(message)

! deallocate parNames
deallocate(parNames, stat=ierr)
if (ierr.ne.0) message='mcmc_wrapper: problem deallocating parNames'


END SUBROUTINE MCMC_WRAPPER
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! (B) OBJECTIVE FUNCTION WRAPPaz
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
SUBROUTINE OBJFUNC_WRAPPER_OPTI(dataIN,dataOUT,argInf,&
                                feas,objFuncM,gradObjFuncM,hessObjFuncM,err,message)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Wrapper for the objective function used in DMSL optimization routines, based on the
!  bateauDK_objFunc_opt wrapper coded by Dmitri Kavetski
! Calls the SUBROUTINE fuse_rmse.f90 to calculate the RMSE for a given
!  FUSE model and parameter set
! ---------------------------------------------------------------------------------------
use kinds_dmsl_kit                       ! numeric kind definitions
use types_dmsl_kit,only:data_ricz_type   ! data types (dataIN,dataOUT; not actually used)
use fuse_rmse_module,only:fuse_rmse      ! provide access to fuse_rmse (run model)
use multiparam,only:lparam,paratt,numpar ! provide access to the FUSE model parameter structures
use multistats,only:fcount               ! provide access to the number of function evaluations
use getpar_str_module                    ! provide access to getpar_str (get parameter metadata)
implicit none
! dummies
type(data_ricz_type),intent(in),optional::dataIN
type(data_ricz_type),intent(inout),optional::dataOUT
real(mrk),intent(in)::argInf(:)
logical(mlk),intent(out)::feas
real(mrk),intent(out),optional::objFuncM,gradObjFuncM(:),hessObjFuncM(:,:)
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
integer(mik)::ipar                    ! loop through model parameters
type(paratt)::param_meta              ! parameter metadata
logical(mlk)::output_flag             ! switch to write model output
logical(mlk)::mparam_flag             ! switch to turn off writing of parameters and statistics
! default error code and message
err=0; message='no error checking'
! define flags to write model output and compute summary statistics
output_flag = .false. 
mparam_flag = .false.
! check for the feasability of the parameters
feas=.true. ! initialize feasability flag
do ipar=1,numpar
 call getpar_str(lparam(ipar)%parname,param_meta)    ! get parameter metadata structure
 if (argInf(ipar).lt.param_meta%parlow) feas=.false. ! check above lower limit
 if (argInf(ipar).gt.param_meta%parupp) feas=.false. ! check below upper limit
 !write(*,'(a11,1x,3(f12.6,1x),l1)') &
 ! lparam(ipar)%parname,argInf(ipar), param_meta%parlow, param_meta%parupp, feas
end do  ! looping through parameters
! calculate objective function and increment counter
if (present(objFuncM) .and. feas) then
 call fuse_rmse(argInf,objFuncM,output_flag,mparam_flag)
endif
if (present(objFuncM)) fcount = fcount+1
!if (present(objFuncM)) write(*,'(i8,1x,20(f9.3,1x))') fcount,objFuncM,argInf
! populate un-used output with missing values
if(present(gradObjFuncM))gradObjFuncM=undefRN
if(present(hessObjFuncM))hessObjFuncM=undefRN
END SUBROUTINE OBJFUNC_WRAPPER_OPTI
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
SUBROUTINE OBJFUNC_WRAPPER_MCMC(dataIN,dataOUT,x,&
                                feas,logp,faux,gradLogP,hessLogP,err,message)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Wrapper for the objective function used in DMSL MCMC routines, based on the
!  bateauDK_objFunc_opt wrapper coded by Dmitri Kavetski
! Calls the SUBROUTINE fuse_rmse.f90 to calculate the RMSE for a given
!  FUSE model and parameter set
! ---------------------------------------------------------------------------------------
! FUSE modules
use fuse_rmse_module,only:fuse_rmse      ! provide access to fuse_rmse (run model)
use multiforce,only:istart,numtim,aforce ! start+count of the calibration period; forcing data
use multiparam,only:lparam,paratt,numpar ! provide access to the FUSE model parameter structures
use multiroute,only:aroute               ! provide access to the FUSE simulated runoff
use multistats,only:fcount               ! provide access to the number of function evaluations
use getpar_str_module                    ! provide access to getpar_str (get parameter metadata)
! DMSL modules
use kinds_dmsl_kit                       ! numeric kind definitions
use types_dmsl_kit,only:data_ricz_type   ! data types (dataIN,dataOUT; not actually used)
USE numerix_dmsl_kit,only:normal_logp    ! log-density of a normal deviate
implicit none
! dummies
type(data_ricz_type),intent(in),optional::dataIN
type(data_ricz_type),intent(inout),optional::dataOUT
real(mrk),intent(in)::x(0:)
logical(mlk),intent(out)::feas
real(mrk),intent(out),optional::logp,faux(:),gradLogP(:),hessLogP(:,:)
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
integer(mik)::ipar                    ! loop through model parameters
integer(mik)::itim                    ! loop through calibration period
type(paratt)::param_meta              ! parameter metadata
logical(mlk)::output_flag             ! switch to write model output
logical(mlk)::mparam_flag             ! switch to turn off writing of parameters and statistics
real(mrk)   ::rmse                    ! root mean squared error
real(mrk),dimension(:),allocatable :: resd ! individual residuals
real(mrk),dimension(:),allocatable :: dens ! log-density of individual residuals
real(mrk)::VAR
! default error code and message
err=0; message='start of fuse wrapper'
! define flags to write model output and compute summary statistics
output_flag = .false. 
mparam_flag = .false.
! check for the feasability of the parameters
feas=.true. ! initialize feasability flag
do ipar=1,numpar
 call getpar_str(lparam(ipar)%parname,param_meta)    ! get parameter metadata structure
 if (x(ipar).lt.param_meta%parlow) feas=.false. ! check above lower limit
 if (x(ipar).gt.param_meta%parupp) feas=.false. ! check below upper limit
 !write(*,'(a11,1x,3(f12.6,1x),l1)') &
 ! lparam(ipar)%parname,x(ipar), param_meta%parlow, param_meta%parupp, feas
end do  ! looping through parameters
! add error checking
if (.not.feas) then
 message='parameter set is infeasible'
 return
endif
! calculate objective function and increment counter
if (present(logp) .and. feas) then
 VAR=10._mrk**x(0)
 call fuse_rmse(x(1:),rmse,output_flag,mparam_flag)
 ! allocate space for log-density of individual residuals
 allocate(resd(istart:numtim),dens(istart:numtim), stat=err)
 if (err.ne.0) then; err=-20; message='problem allocating space for dens'; return; endif
 ! loop thru time steps to get log-density of individual residuals
 do itim=istart,numtim
  resd(itim) = AROUTE(itim)%Q_ROUTED - AFORCE(itim)%OBSQ 
  dens(itim) = normal_logp(x=resd(itim),mean=0._mrk,var=VAR)
 end do
 logp = sum(dens)  ! log density of the simulation
 ! deallocate space for log-density of individual residuals
 deallocate(resd,dens, stat=err)
 if (err.ne.0) then; err=-30; message='problem deallocating space for dens'; return; endif
endif
if (present(logp)) fcount = fcount+1
!if (present(logp)) write(*,'(i8,1x,20(f9.3,1x))') fcount,logp,x
! populate un-used output with missing values
if(present(gradLogP))gradLogP=undefRN
if(present(hessLogP))hessLogP=undefRN
END SUBROUTINE OBJFUNC_WRAPPER_MCMC
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
END MODULE DMSL_WRAPPER_MODULE
