SUBROUTINE MSTATE_EQN()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Computes derivatives of all states for all model combinations
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multistate -- populates the MODULE multistate with derivatives DY_DT%(*)
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition structure
USE model_defnames
USE multiparam                                        ! model parameters
USE multiforce                                        ! model forcing data
USE multistate                                        ! model states
USE multi_flux                                        ! model fluxes
! ---------------------------------------------------------------------------------------
! (1) COMPUTE DERIVATIVES FOR STATES IN THE UPPER LAYER
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iARCH1)
 CASE(iopt_tension2_1) ! tension storage sub-divided into recharge and excess
  DY_DT%TENS_1A = M_FLUX%EFF_PPT - M_FLUX%QSURF - M_FLUX%EVAP_1A - M_FLUX%RCHR2EXCS
  DY_DT%TENS_1B = M_FLUX%RCHR2EXCS - M_FLUX%EVAP_1B - M_FLUX%TENS2FREE_1
  DY_DT%FREE_1  = M_FLUX%TENS2FREE_1 - M_FLUX%QPERC_12 - M_FLUX%QINTF_1 - M_FLUX%OFLOW_1
  !print *, M_FLUX%EFF_PPT, M_FLUX%QSURF, M_FLUX%EVAP_1A, M_FLUX%RCHR2EXCS
 CASE(iopt_tension1_1) ! upper layer broken up into tension and free storage
  DY_DT%TENS_1  = M_FLUX%EFF_PPT - M_FLUX%QSURF - M_FLUX%EVAP_1 - M_FLUX%TENS2FREE_1
  DY_DT%FREE_1  = M_FLUX%TENS2FREE_1 - M_FLUX%QPERC_12 - M_FLUX%QINTF_1 - M_FLUX%OFLOW_1
  !print *, 'in mstate_eqn, layer1 ', DY_DT%TENS_1, DY_DT%FREE_1, M_FLUX%EFF_PPT, M_FLUX%QSURF, M_FLUX%EVAP_1, &
  !                                   M_FLUX%TENS2FREE_1, M_FLUX%QPERC_12, M_FLUX%QINTF_1, M_FLUX%OFLOW_1
 CASE(iopt_onestate_1) ! upper layer defined by a single state variable
  DY_DT%WATR_1  = M_FLUX%EFF_PPT - M_FLUX%QSURF - M_FLUX%EVAP_1 - M_FLUX%QPERC_12 - M_FLUX%QINTF_1 &
                  - M_FLUX%OFLOW_1
  !print *, 'in mstate_eqn, layer1 ', DY_DT%WATR_1, M_FLUX%EFF_PPT, M_FLUX%QSURF, M_FLUX%EVAP_1, &
  !                                   M_FLUX%QPERC_12, M_FLUX%QINTF_1, M_FLUX%OFLOW_1
 CASE DEFAULT
  print *, "SMODL%iARCH1 must be iopt_tension2_1, iopt_tension1_1, or iopt_onestate_1"
  STOP
END SELECT  ! (upper layer architechure)
! ---------------------------------------------------------------------------------------
! (2) COMPUTE DERIVATIVES FOR STATES IN THE LOWER LAYER
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iARCH2)
 CASE(iopt_tens2pll_2) ! tension reservoir plus two parallel tanks
  DY_DT%TENS_2  = M_FLUX%QPERC_12*(1._SP-MPARAM%PERCFRAC) - M_FLUX%EVAP_2 - M_FLUX%TENS2FREE_2
  DY_DT%FREE_2A = M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP - M_FLUX%QBASE_2A &
                  - M_FLUX%OFLOW_2A
  DY_DT%FREE_2B = M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP - M_FLUX%QBASE_2B &
                  - M_FLUX%OFLOW_2B
  !print *, 'in mstate_eqn, layer2 ', M_FLUX%QPERC_12, M_FLUX%EVAP_2, M_FLUX%TENS2FREE_2, M_FLUX%QBASE_2A, M_FLUX%QBASE_2B
 CASE(iopt_unlimfrc_2,iopt_unlimpow_2,iopt_topmdexp_2,iopt_fixedsiz_2) ! single state
  ! (NOTE: M_FLUX%OFLOW_2=0 for 'unlimfrc_2','unlimpow_2','topmdexp_2') 
  DY_DT%WATR_2  = M_FLUX%QPERC_12 - M_FLUX%EVAP_2 - M_FLUX%QBASE_2 - M_FLUX%OFLOW_2
  !print *, 'in mstate_eqn, layer2 ', M_FLUX%EVAP_2, M_FLUX%QBASE_2, M_FLUX%OFLOW_2
 CASE DEFAULT
  print *, "SMODL%iARCH2 must be iopt_tens2pll_2, iopt_unlimfrc_2, iopt_unlimpow_2"
  print *, "  iopt_topmdexp_2, or iopt_fixedsiz_2"
  STOP
END SELECT
! ---------------------------------------------------------------------------------------
END SUBROUTINE MSTATE_EQN
