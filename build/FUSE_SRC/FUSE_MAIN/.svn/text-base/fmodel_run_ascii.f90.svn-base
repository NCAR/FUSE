SUBROUTINE FMODEL_RUN_ASCII()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Run a single model with one parameter set (ASCII output)
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
! data modules
USE model_defn, ONLY: OUTFILE_UNIT                    ! file unit for ASCII output
USE multiforce                                        ! model forcing data
USE multiparam                                        ! model parameters
USE multi_flux                                        ! model fluxes
USE multistate                                        ! model states
USE multiroute                                        ! routed runoff
USE multistats                                        ! summary statistics
! informational modules
USE par_insert_module                                 ! insert parameters into data structures                          
IMPLICIT NONE
! internal
INTEGER(I4B)                           :: ITIM        ! loop through time series
! ---------------------------------------------------------------------------------------
! increment parameter counter
PCOUNT = PCOUNT + 1
! initialize model states and model time step
CALL INIT_STATE(fracState0) ! fracState0 is shared in MODULE multistate
HSTATE%STEP = DELTIM        ! deltim is shared in module multiforce.
! write header for time series output
WRITE(OUTFILE_UNIT,'(A4,1X,3(A2,1X),8(A12,1X))') &
 'YEAR','MM','DD','HH','PPT','EFF_PPT','PET','WATR_1','WATR_2','Q_INSTNT','Q_ROUTED','OBSQ'
! loop through time
DO ITIM=1,NUMTIM            ! (NUMTIM is shared in MODULE multiforce)
 ! run model for one time step
 MFORCE = AFORCE(ITIM)      ! assign model forcing data
 CALL INITFLUXES()          ! set weighted sum of fluxes to zero
 CALL SUBSTEPPER()          ! run model for one time step using implicit solution with variable sub-steps
 CALL Q_OVERLAND()          ! overland flow routing
 ! save instantaneous and routed runoff
 AROUTE(ITIM)%Q_INSTNT = MROUTE%Q_INSTNT  ! save instantaneous runoff
 AROUTE(ITIM)%Q_ROUTED = MROUTE%Q_ROUTED  ! save routed runoff
 ! write model output to ASCII output file
 WRITE(OUTFILE_UNIT,'(I4,1X,3(I2,1X),8(ES12.5,1X))') &
  MFORCE%IY,MFORCE%IM,MFORCE%ID,MFORCE%IH, &
  MFORCE%PPT,W_FLUX%EFF_PPT,MFORCE%PET, &
  FSTATE%WATR_1,FSTATE%WATR_2, &
  MROUTE%Q_INSTNT,MROUTE%Q_ROUTED, &
  MFORCE%OBSQ
END DO  ! (itim)
! ---------------------------------------------------------------------------------------
END SUBROUTINE FMODEL_RUN_ASCII
