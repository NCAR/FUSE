MODULE FUSE_RMSE_MODULE  ! have as a module because of dynamic arrays
IMPLICIT NONE
CONTAINS
SUBROUTINE FUSE_RMSE(XPAR,RMSE,OUTPUT_FLAG,MPARAM_FLAG,ERR,MESSAGE)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! Modified by Brian Henn to include snow model, 6/2013
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Calculate the RMSE for single FUSE model and single parameter set
!   input: model parameter set
!  output: root mean squared error
! ---------------------------------------------------------------------------------------
USE nrtype,only:SP,LGT                                 ! variable types, etc.
! data modules
USE model_defn, ONLY:NSTATE,SMODL                      ! number of state variables
USE model_defnames                                     ! integer model definitions
USE multiparam, ONLY:LPARAM,NUMPAR,MPARAM              ! list of model parameters
USE multiforce, ONLY:MFORCE,AFORCE,DELTIM,ISTART,&     ! model forcing data
                     NUMTIM                            ! model forcing data (continued)
USE multistate, ONLY:fracstate0,TSTATE,MSTATE,FSTATE,& ! model states
                     HSTATE                            ! model states (continued)
USE multiroute, ONLY:MROUTE,AROUTE                     ! routed runoff
USE multistats, ONLY:MSTATS,PCOUNT,MOD_IX              ! access model statistics; counter for param set
! informational modules
USE par_insert_module,only:PUT_PARSET                  ! insert parameters into data structures                          
USE str_2_xtry_module,only:STR_2_XTRY                  ! provide access to the routine str_2_xtry
! interface blocks
USE interfaceb, ONLY:ode_int,fuse_solve                ! provide access to FUSE_SOLVE through ODE_INT
! model numerix structures
!USE model_numerix
!USE fuse_deriv_module
!USE fdjac_ode_module
IMPLICIT NONE
! input
REAL(SP),DIMENSION(:),INTENT(IN)       :: XPAR         ! model parameter set
LOGICAL(LGT), INTENT(IN)               :: OUTPUT_FLAG  ! .TRUE. if desire time series output
LOGICAL(LGT), INTENT(IN), OPTIONAL     :: MPARAM_FLAG  ! .FALSE. (used to turn off writing statistics)
! output
REAL(SP),INTENT(OUT)                   :: RMSE         ! root mean squared error
INTEGER(I4B),INTENT(OUT)               :: ERR          ! error indicator
CHARACTER(LEN=256),INTENT(OUT)         :: MESSAGE      ! error message
! internal
REAL(SP)                               :: T1,T2        ! CPU time
INTEGER(I4B)                           :: ITIM         ! loop through time series
INTEGER(I4B)                           :: IPAR         ! loop through model parameters
REAL(SP)                               :: DT_SUB       ! length of sub-step
REAL(SP)                               :: DT_FULL      ! length of time step
REAL(SP), DIMENSION(:), ALLOCATABLE    :: STATE0       ! vector of model states at the start of the time step
REAL(SP), DIMENSION(:), ALLOCATABLE    :: STATE1       ! vector of model states at the end of the time step
REAL(SP), DIMENSION(:,:), ALLOCATABLE  :: J            ! used to compute the Jacobian (just as a test)
REAL(SP), DIMENSION(:), ALLOCATABLE    :: DSDT         ! used to compute the ODE (just as a test)
INTEGER(I4B)                           :: ITEST,JTEST  ! used to compute a grid of residuals
REAL(SP)                               :: TEST_A,TEST_B ! used to compute a grid of residuals

! ---------------------------------------------------------------------------------------
ERR=0

! allocate state vectors
ALLOCATE(STATE0(NSTATE),STATE1(NSTATE),STAT=ERR)

IF (ERR.NE.0) THEN
  print *, 'f-fuse_rmse/problem allocating state vectors'
  !return
ENDIF
! increment parameter counter for model output (shared in module MULTISTATS)
IF (.NOT.PRESENT(MPARAM_FLAG)) THEN
 PCOUNT = PCOUNT + 1
ELSE
 IF (MPARAM_FLAG) PCOUNT = PCOUNT + 1
ENDIF
! add parameter set to the data structure
CALL PUT_PARSET(XPAR)
!DO IPAR=1,NUMPAR; WRITE(*,'(A11,1X,F9.3)') LPARAM(IPAR), XPAR(IPAR); END DO
! compute derived model parameters (bucket sizes, etc.)
CALL PAR_DERIVE(ERR,message)
IF (ERR.NE.0) then
 message= 'f-FUSE_RMSE/problem allocating state vectors'
 return
endif
! initialize model states and model time step
CALL INIT_STATE(fracState0)            ! fracState0 is shared in MODULE multistate
CALL STR_2_XTRY(FSTATE,STATE0)         ! get the vector of states from the FSTATE structure
DT_SUB  = DELTIM                       ! init stepsize to full step (DELTIM shared in module multiforce)
DT_FULL = DELTIM                       ! init stepsize to full step (DELTIM shared in module multiforce)
! initialize summary statistics
CALL INIT_STATS()
CALL CPU_TIME(T1)
! loop through time
DO ITIM=1,NUMTIM            ! (NUMTIM is shared in MODULE multiforce)
 ! PRINT *, "FUSE_RMSE, ITIM = ",ITIM
 ! run model for one time step
 MFORCE = AFORCE(ITIM)      ! assign model forcing data
 MSTATE = FSTATE            ! refresh model states 
 CALL INITFLUXES()          ! set weighted sum of fluxes to zero
 ! if snow model, call UPDATE_SWE first to calculate snow fluxes and update snow bands 
 ! using explicit Euler approach; if not, call QRAINERROR
 SELECT CASE(SMODL%iSNOWM)
  CASE(iopt_temp_index)
   CALL UPDATE_SWE(DELTIM)
  CASE(iopt_no_snowmod)
   CALL QRAINERROR()
  CASE DEFAULT
   message="f-fuse_rmse/SMODL%iSNOWM must be either iopt_temp_index or iopt_no_snowmod"
   return
 END SELECT
 ! temporally integrate the ordinary differential equations
 CALL ODE_INT(FUSE_SOLVE,STATE0,STATE1,DT_SUB,DT_FULL,ERR,MESSAGE)
 IF (ERR.NE.0) THEN
  message="f-fuse_rmse/&"//message
  return
 ENDIF
 ! perform overland flow routing
 CALL Q_OVERLAND()
 ! save state
 STATE0=STATE1
 ! save instantaneous and routed runoff
 AROUTE(ITIM)%Q_INSTNT = MROUTE%Q_INSTNT  ! save instantaneous runoff
 AROUTE(ITIM)%Q_ROUTED = MROUTE%Q_ROUTED  ! save routed runoff
 IF (AROUTE(ITIM)%Q_ROUTED.LT.0._sp) STOP ' Q_ROUTED is less than zero '
 IF (AROUTE(ITIM)%Q_ROUTED.GT.1000._sp) STOP ' Q_ROUTED is enormous '
 ! compute summary statistics
 CALL COMP_STATS()
 ! write model output
 IF (OUTPUT_FLAG) THEN
  !PRINT *, 'OUTPUT_FLAG is TRUE, now writing model output'
  CALL PUT_OUTPUT(PCOUNT,MOD_IX,ITIM)
  !WRITE(*,'(I10,1X,2(F15.8,1X))') ITIM, FSTATE%WATR_1, FSTATE%WATR_2
  !WRITE(*,'(I10,1X,I4,1X,4(I2,1X),F9.3,1X,F20.1,1X,4(F11.3,1X))') ITIM, AFORCE(ITIM), AROUTE(ITIM)%Q_ROUTED
 ENDIF
END DO  ! (itim)
CALL CPU_TIME(T2)
WRITE(*,*) "TIME ELAPSED = ", t2-t1
! calculate mean summary statistics
CALL MEAN_STATS() 
RMSE = MSTATS%RAW_RMSE
! WRITE(unt,'(2(I6,1X),3(F20.15,1X))') MOD_IX, PCOUNT, MSTATS%RAW_RMSE, MSTATS%NASH_SUTT, MSTATS%NUM_FUNCS
! write model parameters and summary statistics
IF (.NOT.PRESENT(MPARAM_FLAG)) THEN
 CALL PUT_PARAMS(PCOUNT,MOD_IX)  ! PCOUNT = index for parameter set; ONEMOD=1 (just one model structure)
 CALL PUT_SSTATS(PCOUNT,MOD_IX)
ELSE
 IF (MPARAM_FLAG) THEN
  CALL PUT_PARAMS(PCOUNT,MOD_IX)  ! PCOUNT = index for parameter set; ONEMOD=1 (just one model structure)
  CALL PUT_SSTATS(PCOUNT,MOD_IX)
 ENDIF
ENDIF
! deallocate state vectors
DEALLOCATE(STATE0,STATE1,STAT=ERR)
IF (ERR.NE.0)THEN
 message='f-fuse_rmse/problem deallocating state vectors'
 return
ENDIF
! ---------------------------------------------------------------------------------------
END SUBROUTINE FUSE_RMSE
END MODULE FUSE_RMSE_MODULE
