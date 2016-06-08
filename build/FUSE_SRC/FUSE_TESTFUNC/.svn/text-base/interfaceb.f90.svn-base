MODULE INTERFACEB
! -------------------------------------------------------------------------------------------------
INTERFACE
 SUBROUTINE ODE_INT(MODL_SOLVE,STATE_START,STATE_END,DT_SUB,DT_FULL,IERR,MESSAGE)
 USE nrtype                                            ! variable definitions, etc.
 IMPLICIT NONE
 REAL(SP), DIMENSION(:), INTENT(IN)     :: STATE_START ! state vector at the start of the full step
 REAL(SP), DIMENSION(:), INTENT(OUT)    :: STATE_END   ! state vector at the end of the full step
 REAL(SP), INTENT(INOUT)                :: DT_SUB      ! length of the sub-step
 REAL(SP), INTENT(IN)                   :: DT_FULL     ! length of the full step
 INTEGER(I4B), INTENT(OUT)              :: IERR        ! error code
 CHARACTER(LEN=*), INTENT(OUT)          :: MESSAGE     ! error message
 INTERFACE
  SUBROUTINE MODL_SOLVE(CALCDSDT,IE_SOLVE,B_IMPOSE,AVG_FLUX,ADD_FLUX,NEWSTATE, & ! define functionality of the routine
                        DT,S0,S1,DSDT,NEWSTEP,CONVCHECK,NITER,SOLUTION,HBOUND, & ! input/output
                        IERR,MESSAGE)                                            ! error control
  USE nrtype                                                   ! variable definitions, etc.
  IMPLICIT NONE
  LOGICAL(LGT), INTENT(IN),OPTIONAL             :: CALCDSDT    ! FLAG to calculate derivatives at S0
  LOGICAL(LGT), INTENT(IN),OPTIONAL             :: IE_SOLVE    ! FLAG to compute the implicit Euler solution
  LOGICAL(LGT), INTENT(IN),OPTIONAL             :: B_IMPOSE    ! FLAG to impose bounds on model states
  LOGICAL(LGT), INTENT(IN),OPTIONAL             :: AVG_FLUX    ! FLAG to average fluxes from start & end states 
  LOGICAL(LGT), INTENT(IN),OPTIONAL             :: ADD_FLUX    ! FLAG to add accepted fluxes to the total flux
  LOGICAL(LGT), INTENT(IN),OPTIONAL             :: NEWSTATE    ! FLAG to use weighted fluxes to compute end state
  REAL(SP), INTENT(IN), OPTIONAL                :: DT          ! length of the sub-step
  REAL(SP), DIMENSION(:),INTENT(IN), OPTIONAL   :: S0          ! input state vector
  REAL(SP), DIMENSION(:), INTENT(OUT),OPTIONAL  :: S1          ! state vector from the implicit euler solution
  REAL(SP), DIMENSION(:),INTENT(INOUT),OPTIONAL :: DSDT        ! state derivatives
  LOGICAL(LGT), INTENT(IN),OPTIONAL             :: NEWSTEP     ! FLAG to denote a new model time step
  LOGICAL(LGT), INTENT(IN),OPTIONAL             :: CONVCHECK   ! FLAG to check for convergence of the implicit scheme
  INTEGER(I4B), INTENT(OUT), OPTIONAL           :: NITER       ! number of iterations
  INTEGER(I4B), INTENT(IN), OPTIONAL            :: SOLUTION    ! solution is at start (0) or end (1) of sub-step
  LOGICAL(LGT), INTENT(OUT),OPTIONAL            :: HBOUND      ! FLAG to denote if the states were out of bounds 
  INTEGER(I4B), INTENT(OUT)                     :: IERR        ! error code
  CHARACTER(LEN=*), INTENT(OUT)                 :: MESSAGE     ! error message
  END SUBROUTINE MODL_SOLVE
 END INTERFACE
 END SUBROUTINE ODE_INT
END INTERFACE
! -------------------------------------------------------------------------------------------------
INTERFACE
 SUBROUTINE TEST_SOLVE(CALCDSDT,IE_SOLVE,B_IMPOSE,AVG_FLUX,ADD_FLUX,NEWSTATE, & ! define functionality of the routine
                       DT,S0,S1,DSDT,NEWSTEP,CONVCHECK,NITER,SOLUTION,HBOUND, & ! input/output
                       IERR,MESSAGE)                                            ! error control
 USE nrtype                                                   ! variable definitions, etc.
 IMPLICIT NONE
 LOGICAL(LGT), INTENT(IN),OPTIONAL             :: CALCDSDT    ! FLAG to calculate derivatives at S0
 LOGICAL(LGT), INTENT(IN),OPTIONAL             :: IE_SOLVE    ! FLAG to compute the implicit Euler solution
 LOGICAL(LGT), INTENT(IN),OPTIONAL             :: B_IMPOSE    ! FLAG to impose bounds on model states
 LOGICAL(LGT), INTENT(IN),OPTIONAL             :: AVG_FLUX    ! FLAG to average fluxes from start & end states 
 LOGICAL(LGT), INTENT(IN),OPTIONAL             :: ADD_FLUX    ! FLAG to add accepted fluxes to the total flux
 LOGICAL(LGT), INTENT(IN),OPTIONAL             :: NEWSTATE    ! FLAG to use weighted fluxes to compute end state
 REAL(SP), INTENT(IN), OPTIONAL                :: DT          ! length of the sub-step
 REAL(SP), DIMENSION(:),INTENT(IN), OPTIONAL   :: S0          ! input state vector
 REAL(SP), DIMENSION(:), INTENT(OUT),OPTIONAL  :: S1          ! state vector from the implicit euler solution
 REAL(SP), DIMENSION(:),INTENT(INOUT),OPTIONAL :: DSDT        ! state derivatives
 LOGICAL(LGT), INTENT(IN),OPTIONAL             :: NEWSTEP     ! FLAG to denote a new model time step
 LOGICAL(LGT), INTENT(IN),OPTIONAL             :: CONVCHECK   ! FLAG to check for convergence of the implicit scheme
 INTEGER(I4B), INTENT(OUT), OPTIONAL           :: NITER       ! number of iterations
 INTEGER(I4B), INTENT(IN), OPTIONAL            :: SOLUTION    ! solution is at start (0) or end (1) of sub-step
 LOGICAL(LGT), INTENT(OUT),OPTIONAL            :: HBOUND      ! FLAG to denote if the states were out of bounds 
 INTEGER(I4B), INTENT(OUT)                     :: IERR        ! error code
 CHARACTER(LEN=*), INTENT(OUT)                 :: MESSAGE     ! error message
 END SUBROUTINE TEST_SOLVE
END INTERFACE
! -------------------------------------------------------------------------------------------------
END MODULE INTERFACEB
