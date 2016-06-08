MODULE meta_stats
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Describe all summary statistics (used to define NetCDF output files, etc.)
! ---------------------------------------------------------------------------------------
! variable definitions
USE nrtype
IMPLICIT NONE
CHARACTER(LEN=11), DIMENSION(100)      :: XNAME       ! variable names
CHARACTER(LEN=52), DIMENSION(100)      :: XDESC       ! variable long names (descrition of variable)
CHARACTER(LEN=13), DIMENSION(100)      :: XUNIT       ! variable units
INTEGER(I4B)                           :: I           ! loop through variables
INTEGER(I4B)                           :: NSUMVAR     ! number of summary variables
CONTAINS
! ---------------------------------------------------------------------------------------
SUBROUTINE SUMDESCRIBE()
I=0  ! initialize counter
! DMSL diagnostix
I=I+1; XNAME(I)='var_residul'; XDESC(I)='variance of the model residuals, used in MCMC      '; XUNIT(I)='mm**2        '
I=I+1; XNAME(I)='logp_simuln'; XDESC(I)='log density of the simulation                      '; XUNIT(I)='problem_depnt'
I=I+1; XNAME(I)='jump_taken '; XDESC(I)='MCMC jump diagnostix; 0 = no jump; 1 = jumping     '; XUNIT(I)='-            '
! comparisons between model output and observations
I=I+1; XNAME(I)='qobs_mean  '; XDESC(I)='mean observed runoff                               '; XUNIT(I)='mm timestep-1'
I=I+1; XNAME(I)='qsim_mean  '; XDESC(I)='mean simulated runoff                              '; XUNIT(I)='mm timestep-1'
I=I+1; XNAME(I)='qobs_cvar  '; XDESC(I)='coefficient of variation of observed runoff        '; XUNIT(I)='-            '
I=I+1; XNAME(I)='qsim_cvar  '; XDESC(I)='coefficient of variation of simulated runoff       '; XUNIT(I)='-            '
I=I+1; XNAME(I)='qobs_lag1  '; XDESC(I)='lag-1 correlation of observed runoff               '; XUNIT(I)='-            '
I=I+1; XNAME(I)='qsim_lag1  '; XDESC(I)='lag-1 correlation of simulated runoff              '; XUNIT(I)='-            '
I=I+1; XNAME(I)='raw_rmse   '; XDESC(I)='root-mean-squared-error of flow                    '; XUNIT(I)='mm timestep-1'
I=I+1; XNAME(I)='log_rmse   '; XDESC(I)='root-mean-squared-error of LOG flow                '; XUNIT(I)='mm timestep-1'
I=I+1; XNAME(I)='nash_sutt  '; XDESC(I)='Nash-Sutcliffe score                               '; XUNIT(I)='-            '
! attributes of model output
I=I+1; XNAME(I)='numerx_rmse'; XDESC(I)='RMSE between exact and approximate solution        '; XUNIT(I)='mm timestep-1'
I=I+1; XNAME(I)='mean_nfuncs'; XDESC(I)='mean number function evaluations                   '; XUNIT(I)='-            '
I=I+1; XNAME(I)='mean_njacob'; XDESC(I)='mean number jacobian evaluations                   '; XUNIT(I)='-            '
I=I+1; XNAME(I)='mean_accept'; XDESC(I)='mean number sub-steps accepted (taken)             '; XUNIT(I)='-            '
I=I+1; XNAME(I)='mean_reject'; XDESC(I)='mean number sub-steps tried but rejected           '; XUNIT(I)='-            '
I=I+1; XNAME(I)='mean_noconv'; XDESC(I)='mean number sub-steps tried that did not converge  '; XUNIT(I)='-            '
I=I+1; XNAME(I)='maxnum_iter'; XDESC(I)='maximum number of iterations in the implicit scheme'; XUNIT(I)='-            '
NSUMVAR=I
END SUBROUTINE SUMDESCRIBE
END MODULE meta_stats
