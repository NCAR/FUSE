MODULE FUSE_RMSE_MODULE  ! have as a module because of dynamic arrays
IMPLICIT NONE
CONTAINS
SUBROUTINE FUSE_RMSE(XPAR,DISTRIBUTED,RMSE,OUTPUT_FLAG,MPARAM_FLAG)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! Modified by Brian Henn to include snow model, 6/2013
! Modified by Nans Addor to enable distributed modeling, 9/2016
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Calculate the RMSE for single FUSE model and single parameter set
!   input: model parameter set
!  output: root mean squared error
! ---------------------------------------------------------------------------------------
USE nrtype                                               ! variable types, etc.
! data modules
USE model_defn, ONLY:NSTATE,SMODL                        ! number of state variables
USE model_defnames                                       ! integer model definitions
USE multiparam, ONLY:LPARAM,NUMPAR,MPARAM                ! list of model parameters
USE multiforce, ONLY:MFORCE,AFORCE,DELTIM,ISTART,&       ! model forcing data
                     NUMTIM,warmup_beg                   ! model forcing data (continued)
USE multiforce, only:nspat1,nspat2                       ! spatial dimensions
USE multiforce, ONLY:gForce                              ! gridded forcing data
USE multistate, ONLY:fracstate0,TSTATE,MSTATE,FSTATE,&   ! model states
                     HSTATE                              ! model states (continued)
USE multistate, ONLY:gState                              ! gridded state variables
USE multiroute, ONLY:MROUTE,AROUTE                       ! routed runoff
USE multistats, ONLY:MSTATS,PCOUNT,MOD_IX                ! access model statistics; counter for param set
! code modules
USE get_gforce_module,only:get_modtim                    ! get model time for a given time step
USE get_gforce_module,only:get_gforce                    ! get gridded forcing data for a given time step
USE getPETgrid_module,only:getPETgrid                    ! get gridded PET
USE par_insert_module                                    ! insert parameters into data structures                          
USE str_2_xtry_module                                    ! provide access to the routine str_2_xtry
USE xtry_2_str_module                                    ! provide access to the routine xtry_2_str
! interface blocks
USE interfaceb, ONLY:ode_int,fuse_solve                  ! provide access to FUSE_SOLVE through ODE_INT
! model numerix structures
USE model_numerix
USE fuse_deriv_module
USE fdjac_ode_module
IMPLICIT NONE
! input
REAL(SP),DIMENSION(:),INTENT(IN)       :: XPAR           ! model parameter set
LOGICAL(LGT), INTENT(IN)               :: DISTRIBUTED    ! .TRUE. if doing distributed simulations
LOGICAL(LGT), INTENT(IN)               :: OUTPUT_FLAG    ! .TRUE. if desire time series output
LOGICAL(LGT), INTENT(IN), OPTIONAL     :: MPARAM_FLAG    ! .FALSE. (used to turn off writing statistics)
! output
REAL(SP),INTENT(OUT)                   :: RMSE           ! root mean squared error
! internal  
logical(lgt),parameter                 :: computePET=.false. ! flag to compute PET
REAL(SP)                               :: T1,T2          ! CPU time
INTEGER(I4B)                           :: ITIM           ! loop through time series
INTEGER(I4B)                           :: iSpat1,iSpat2  ! loop through spatial dimensions
INTEGER(I4B)                           :: IPAR           ! loop through model parameters
REAL(SP)                               :: DT_SUB         ! length of sub-step
REAL(SP)                               :: DT_FULL        ! length of time step
REAL(SP), DIMENSION(:), ALLOCATABLE    :: STATE0         ! vector of model states at the start of the time step
REAL(SP), DIMENSION(:), ALLOCATABLE    :: STATE1         ! vector of model states at the end of the time step
REAL(SP), DIMENSION(:,:), ALLOCATABLE  :: J              ! used to compute the Jacobian (just as a test)
REAL(SP), DIMENSION(:), ALLOCATABLE    :: DSDT           ! used to compute the ODE (just as a test)
INTEGER(I4B)                           :: ITEST,JTEST    ! used to compute a grid of residuals
REAL(SP)                               :: TEST_A,TEST_B  ! used to compute a grid of residuals
INTEGER(I4B)                           :: IERR           ! error code
INTEGER(I4B), PARAMETER                :: CLEN=1024      ! length of character string
INTEGER(I4B)                           :: ERR            ! error code
CHARACTER(LEN=CLEN)                    :: MESSAGE        ! error message
CHARACTER(LEN=CLEN)                    :: CMESSAGE       ! error message of downwind routine
INTEGER(I4B),PARAMETER::UNT=6  !1701 ! 6
! ---------------------------------------------------------------------------------------
! allocate state vectors
ALLOCATE(STATE0(NSTATE),STATE1(NSTATE),STAT=IERR)
IF (IERR.NE.0) STOP ' problem allocating space for state vectors in fuse_rmse '
! increment parameter counter for model output (shared in module MULTISTATS)
IF (.NOT.PRESENT(MPARAM_FLAG)) THEN
 PCOUNT = PCOUNT + 1
ELSE
 IF (MPARAM_FLAG) PCOUNT = PCOUNT + 1
ENDIF
! add parameter set to the data structure
CALL PUT_PARSET(XPAR)
print *, 'Parameter set added to data structure'
!DO IPAR=1,NUMPAR; WRITE(*,'(A11,1X,F9.3)') LPARAM(IPAR), XPAR(IPAR); END DO
! compute derived model parameters (bucket sizes, etc.)
CALL PAR_DERIVE(ERR,MESSAGE)
IF (ERR.NE.0) WRITE(*,*) TRIM(MESSAGE); IF (ERR.GT.0) STOP
! initialize model states
do iSpat1=1,nSpat1
 do iSpat2=1,nSpat2
  CALL INIT_STATE(fracState0)          ! fracState0 is shared in MODULE multistate
  gState(iSpat1,iSpat2) = FSTATE       ! put the state into the 2-d structure
 end do
end do
print *, 'Model state initialized'

! initialize model time step
DT_SUB  = DELTIM                       ! init stepsize to full step (DELTIM shared in module multiforce)
DT_FULL = DELTIM                       ! init stepsize to full step (DELTIM shared in module multiforce)
! initialize summary statistics
CALL INIT_STATS()
CALL CPU_TIME(T1)
! --------------------------------------------------------------------------------------------------------
! loop through time

print *, 'Running the model'
DO ITIM=1,NUMTIM            ! (NUMTIM is shared in MODULE multiforce)
 ! if not distributed (i.e., lumped)
 if(.not.distributed)then
  ! retrieve data from memory
  gForce(1,1) = aForce(iTim)
 else ! distributed 

  ! get the model time
  call get_modtim(warmup_beg+itim,ierr,message)
  if(ierr/=0)then; print*, trim(cmessage); stop; endif
  ! get the gridded model forcing data

  call get_gforce(warmup_beg+itim,ierr,cmessage)
  if(ierr/=0)then; print*, trim(cmessage); stop; endif
  ! compute potential ET
  if(computePET) call getPETgrid(ierr,cmessage)
  if(ierr/=0)then; print*, trim(cmessage); stop; endif
 endif

 ! -------------------------------------------------------------------------------------------------------
 ! loop through grid points, and run the model for one time step
 do iSpat1=1,nSpat1
  do iSpat2=1,nSpat2
   ! extract forcing data
   MFORCE = gForce(iSpat1,iSpat2)      ! assign model forcing data
   ! extract model states
   MSTATE = gState(iSpat1,iSpat2)      ! refresh model states
   FSTATE = gState(iSpat1,iSpat2)      ! refresh model states
   ! get the vector of model states from the structure
   CALL STR_2_XTRY(FSTATE,STATE0)      ! get the vector of states from the FSTATE structure
   ! initialize model fluxes
   CALL INITFLUXES()                   ! set weighted sum of fluxes to zero

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
   CALL ODE_INT(FUSE_SOLVE,STATE0,STATE1,DT_SUB,DT_FULL,IERR,MESSAGE)
   IF (IERR.NE.0) THEN; PRINT *, TRIM(MESSAGE); PAUSE; ENDIF
   ! perform overland flow routing
   CALL Q_OVERLAND()

   ! save the state
   CALL XTRY_2_STR(STATE1,FSTATE)      ! update FSTATE
   gState(iSpat1,iSpat2) = FSTATE      ! put the state into the 2-d structure
   ! save forcing data
   if(distributed)then
    aForce(iTim)%ppt = sum(gForce(:,:)%ppt)/real(size(gForce), kind(sp))
    aForce(iTim)%pet = sum(gForce(:,:)%pet)/real(size(gForce), kind(sp))
   endif

   ! save instantaneous and routed runoff
   AROUTE(ITIM)%Q_INSTNT = MROUTE%Q_INSTNT  ! save instantaneous runoff
   AROUTE(ITIM)%Q_ROUTED = MROUTE%Q_ROUTED  ! save routed runoff
   ! sanity check
   IF (AROUTE(ITIM)%Q_ROUTED.LT.0._sp) STOP ' Q_ROUTED is less than zero '
   IF (AROUTE(ITIM)%Q_ROUTED.GT.1000._sp) STOP ' Q_ROUTED is enormous '

   ! compute summary statistics
   CALL COMP_STATS()
   ! write model output
   IF (OUTPUT_FLAG) THEN
    CALL PUT_OUTPUT(PCOUNT,MOD_IX,iSpat1,iSpat2,ITIM)
    !WRITE(*,'(I10,1X,2(F15.8,1X))') ITIM, FSTATE%WATR_1, FSTATE%WATR_2
    !WRITE(*,'(I10,1X,I4,1X,4(I2,1X),F9.3,1X,F20.1,1X,4(F11.3,1X))') ITIM, AFORCE(ITIM), AROUTE(ITIM)%Q_ROUTED
   ENDIF
  end do  ! (looping thru 2nd spatial dimension)
 end do  ! (looping thru 1st spatial dimension)
END DO  ! (itim)

! get timing information
CALL CPU_TIME(T2)
WRITE(*,*) "TIME ELAPSED = ", t2-t1
! calculate mean summary statistics

CALL MEAN_STATS() 
RMSE = MSTATS%RAW_RMSE
WRITE(unt,'(2(I6,1X),3(F20.15,1X))') MOD_IX, PCOUNT, MSTATS%RAW_RMSE, MSTATS%NASH_SUTT, MSTATS%NUM_FUNCS
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
DEALLOCATE(STATE0,STATE1,STAT=IERR); IF (IERR.NE.0) STOP ' problem deallocating state vectors in fuse_rmse '
! ---------------------------------------------------------------------------------------
END SUBROUTINE FUSE_RMSE
END MODULE FUSE_RMSE_MODULE
