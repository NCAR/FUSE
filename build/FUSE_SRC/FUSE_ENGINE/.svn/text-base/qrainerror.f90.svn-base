SUBROUTINE QRAINERROR()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2008
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Computes the "effective" rainfall, following an error model
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multi_flux -- "effective" rainfall (eff_ppt) stored in  MODULE multi_flux
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition structure
USE model_defnames
USE multiforce                                        ! model forcing
USE multiparam                                        ! model parameters
USE multi_flux                                        ! model fluxes
IMPLICIT NONE
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iRFERR)
 CASE(iopt_additive_e) ! additive rainfall error
  M_FLUX%EFF_PPT = MAX(0.0_sp, MFORCE%PPT + MPARAM%RFERR_ADD)
 CASE(iopt_multiplc_e) ! multiplicative rainfall error
  M_FLUX%EFF_PPT = MFORCE%PPT * MPARAM%RFERR_MLT
 CASE DEFAULT       ! check for errors
  print *, "SMODL%iRFERR must be either iopt_additive_e or iopt_multiplc_e"
  STOP
END SELECT
! ---------------------------------------------------------------------------------------
END SUBROUTINE QRAINERROR
