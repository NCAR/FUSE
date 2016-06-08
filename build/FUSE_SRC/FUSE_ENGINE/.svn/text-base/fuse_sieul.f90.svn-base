MODULE FUSE_SIEUL_MODULE
IMPLICIT NONE
CONTAINS
! ---------------------------------------------------------------------------------------
SUBROUTINE FUSE_SIEUL(SINI,DSDT0,DT,IERR,MESSAGE)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! A FUSE-specific routine for the temporal integration of ordinary differential equations
!  using the semi-implicit Euler method
! ---------------------------------------------------------------------------------------
! Modules Modified:
! --------
! Populates M_FLUX in module multi_flux (in the routine disaggflux)
! ---------------------------------------------------------------------------------------
USE nrtype                                                   ! data types
USE nrutil, ONLY: diagadd                                    ! utility to add identity matrix
USE nr, ONLY: ludcmp,lubksb                                  ! provide access to the LU solver
USE fdjac_ode_module                                         ! provide access to fdjac_ode
USE disaggflux_module                                        ! provide access to disaggflux
IMPLICIT NONE
! input
REAL(SP), DIMENSION(:), INTENT(IN)            :: SINI        ! initial state vector
REAL(SP), DIMENSION(:), INTENT(IN)            :: DSDT0       ! initial state derivatives
REAL(SP), INTENT(IN)                          :: DT          ! time step
! internal
INTEGER(I4B)                                  :: ISTT        ! looping through states
REAL(SP), DIMENSION(SIZE(SINI))               :: STRY        ! trial state vector, used in FDJAC_ODE
REAL(SP), DIMENSION(SIZE(SINI),SIZE(SINI))    :: JAC_ODE     ! Jacobian of the ODE
REAL(SP), DIMENSION(SIZE(SINI),SIZE(SINI))    :: FJAC        ! Jacobian matrix
INTEGER(I4B), DIMENSION(SIZE(SINI))           :: INDX        ! Row permutations from partial pivoting (LUDCMP)
REAL(SP)                                      :: D           ! Denotes the number of row interchanges (LUDCMP)
REAL(SP), DIMENSION(SIZE(SINI))               :: DELS        ! Change in state variables
LOGICAL(LGT)                                  :: EFLAG       ! Error flag
! output -- note: derivatives stored in the FUSE data structures)
INTEGER(I4B), INTENT(OUT)                     :: IERR        ! error code
CHARACTER(*), INTENT(OUT)                     :: MESSAGE     ! error message
! ---------------------------------------------------------------------------------------
! initialize errors
IERR=0; MESSAGE='fuse_sieul: just started'
! calculate Jacobian at S(n) -- and also calculate flux derivatives (SIMETH=.true.)
STRY=SINI ! need to calculate Jacobian at S(n), but want to preserve SINI
CALL FDJAC_ODE(STRY,DSDT0,JAC_ODE,SIMETH=.TRUE.) ! calculate Jacobian of the ODE
FJAC=-DT*JAC_ODE; CALL DIAGADD(FJAC,1._SP)       ! compute (I - DT dg/dS)
! preliminaries before solving linear system
DELS=DT*DSDT0                     ! set up RHS of the linear system
! solve linear system delS = Jac**-1 dt*dSdt
CALL LUDCMP(FJAC,INDX,D)          ! decompose Jacobian
CALL LUBKSB(FJAC,INDX,DELS)       ! solve for delS
! disaggregate fluxes
CALL DISAGGFLUX(DELS,EFLAG)       ! disaggregate fluxes (store in structure M_FLUX)
! process warning (negative error code)
IF (EFLAG) THEN; IERR=-20; MESSAGE='fuse_sieul: unusual flux calculation; truncated'; ENDIF
! re-compute derivatives (use structure M_FLUX populated in disaggflux)
CALL MSTATE_EQN()
! ---------------------------------------------------------------------------------------
END SUBROUTINE FUSE_SIEUL
END MODULE FUSE_SIEUL_MODULE
