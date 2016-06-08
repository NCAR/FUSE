SUBROUTINE GETNUMERIX(err,message)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Reads decisions/parameters that defines the numerical scheme 
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE model_numerix -- model parameters stored in MODULE multiparam
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE fuse_fileManager,only:SETNGS_PATH,MOD_NUMERIX     ! defines data directory                               
USE model_numerix                                     ! defines numerix decisions
IMPLICIT NONE
! dummies
integer(I4B),intent(out)::err
character(*),intent(out)::message
! locals
INTEGER(I4B)                           :: IUNIT       ! file unit
integer(i4b),parameter::lenPath=1024 !DK211008: allows longer file paths
CHARACTER(LEN=lenPath)                 :: CFILE       ! name of constraints file
LOGICAL(LGT)                           :: LEXIST      ! .TRUE. if file exists
! ---------------------------------------------------------------------------------------
! read in control file
err=0; message="GETNUMERIX/ok"
IUNIT = 21  ! file unit
CFILE = TRIM(SETNGS_PATH)//TRIM(MOD_NUMERIX)      ! control file info shared in MODULE ddirectory
INQUIRE(FILE=CFILE,EXIST=LEXIST)  ! check that control file exists
IF (LEXIST) THEN
 ! open up model numerix file
 OPEN(IUNIT,FILE=CFILE,STATUS='old')
  READ(IUNIT,*) SOLUTION_METHOD         ! Method used to solve state equations (explicit vs implicit)
  READ(IUNIT,*) TEMPORAL_ERROR_CONTROL  ! Method used for temporal error control (adaptive time steps)
  READ(IUNIT,*) INITIAL_NEWTON          ! Method used to estimate the initial conditions for the Newton scheme
  READ(IUNIT,*) JAC_RECOMPUTE           ! Jacobian re-evaluation strategy
  READ(IUNIT,*) CHECK_OVERSHOOT         ! Method used to trap/fix errors in Newton
  READ(IUNIT,*) SMALL_ENDSTEP           ! Method used to process the small time interval at the end of a time step
  READ(IUNIT,*) ERR_TRUNC_ABS           ! Absolute temporal truncation error tolerance
  READ(IUNIT,*) ERR_TRUNC_REL           ! Relative temporal truncation error tolerance
  READ(IUNIT,*) ERR_ITER_FUNC           ! Iteration convergence tolerance for function values
  READ(IUNIT,*) ERR_ITER_DX             ! Iteration convergence tolerance for dx
  READ(IUNIT,*) THRESH_FRZE             ! Threshold for freezing the Jacobian
  READ(IUNIT,*) FRACSTATE_MIN           ! Fractional minimum value of state (used so that derivatives are non-zero) 
  READ(IUNIT,*) SAFETY                  ! Safety factor in step-size equation
  READ(IUNIT,*) RMIN                    ! Minimum step size multiplier
  READ(IUNIT,*) RMAX                    ! Maximum step size multiplier
  READ(IUNIT,*) NITER_TOTAL             ! Total number of iterations used in the implicit scheme
  READ(IUNIT,*) MIN_TSTEP               ! Minimum time step length (minutes)
  READ(IUNIT,*) MAX_TSTEP               ! Maximum time step length (minutes)
 CLOSE(IUNIT)
 MIN_TSTEP = MIN_TSTEP/(24._SP*60._SP)  ! Convert from minutes to days
 MAX_TSTEP = MAX_TSTEP/(24._SP*60._SP)  ! Convert from minutes to days
ELSE
 message="f-GETNUMERIX/model numerix file '"//trim(CFILE)//"' does not exist"
 err=100; return
ENDIF 
END SUBROUTINE GETNUMERIX
