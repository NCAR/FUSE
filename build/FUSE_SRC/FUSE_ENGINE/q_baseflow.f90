SUBROUTINE Q_BASEFLOW()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Computes the baseflow from the lower soil layer
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multi_flux -- baseflow stored in  MODULE multi_flux
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition structure
USE model_defnames
USE multiparam                                        ! model parameters
USE multistate                                        ! model states
USE multi_flux                                        ! model fluxes
IMPLICIT NONE
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iARCH2)
 ! --------------------------------------------------------------------------------------
 CASE(iopt_tens2pll_2) ! tension reservoir plus two parallel tanks
  M_FLUX%QBASE_2A = MPARAM%QBRATE_2A * TSTATE%FREE_2A    ! qbrate_2a is a fraction (T-1)
  M_FLUX%QBASE_2B = MPARAM%QBRATE_2B * TSTATE%FREE_2B    ! qbrate_2b is a fraction (T-1)
  M_FLUX%QBASE_2  = M_FLUX%QBASE_2A + M_FLUX%QBASE_2B    ! total baseflow
  !WRITE(*,'(3(F9.3,1X))') MPARAM%QBRATE_2A, TSTATE%FREE_2A, M_FLUX%QBASE_2A
  !WRITE(*,'(4(F9.3,1X))') MPARAM%QBRATE_2B, TSTATE%FREE_2B, M_FLUX%QBASE_2B, M_FLUX%QBASE_2
 ! --------------------------------------------------------------------------------------
 CASE(iopt_unlimfrc_2) ! baseflow resvr of unlimited size (0-HUGE), frac rate
  M_FLUX%QBASE_2  = MPARAM%QB_PRMS * TSTATE%WATR_2       ! qb_prms is a fraction (T-1)
 ! --------------------------------------------------------------------------------------
 CASE(iopt_unlimpow_2) ! baseflow resvr of unlimited size (0-HUGE), power recession
  M_FLUX%QBASE_2  = DPARAM%QBSAT * (TSTATE%WATR_2/MPARAM%MAXWATR_2)**MPARAM%QB_POWR
 ! --------------------------------------------------------------------------------------
 CASE(iopt_topmdexp_2) ! topmodel exponential reservoir (-HUGE to HUGE)
  M_FLUX%QBASE_2  = DPARAM%QBSAT * EXP( -(1. - TSTATE%WATR_2/MPARAM%MAXWATR_2) )
 ! --------------------------------------------------------------------------------------
 CASE(iopt_fixedsiz_2) ! baseflow reservoir of fixed size
  M_FLUX%QBASE_2  = MPARAM%BASERTE * (TSTATE%WATR_2/MPARAM%MAXWATR_2)**MPARAM%QB_POWR
 ! --------------------------------------------------------------------------------------
 CASE DEFAULT
  print *, "SMODL%iARCH2 must be iopt_tens2pll_2, iopt_unlimfrc_2, iopt_unlimpow_2"
  print *, "  iopt_topmdexp_2, or iopt_fixedsiz_2"
  STOP
 ! --------------------------------------------------------------------------------------
END SELECT
! ---------------------------------------------------------------------------------------
END SUBROUTINE Q_BASEFLOW
