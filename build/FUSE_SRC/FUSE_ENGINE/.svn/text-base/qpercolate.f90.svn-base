SUBROUTINE QPERCOLATE()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Computes the percolation from the upper soil layer to the lower soil layer
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multi_flux -- percolation stored in  MODULE multi_flux
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition structure
USE model_defnames
USE multiparam                                        ! model parameters
USE multistate                                        ! model states
USE multi_flux                                        ! model fluxes
IMPLICIT NONE
REAL(SP)                               :: LZ_PD       ! lower zone percolation demand
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iQPERC)
 CASE(iopt_perc_f2sat) ! water from (field cap to sat) avail for percolation
  M_FLUX%QPERC_12 = MPARAM%PERCRTE * (TSTATE%FREE_1/DPARAM%MAXFREE_1)**MPARAM%PERCEXP
 CASE(iopt_perc_w2sat) ! water from (wilt pt to sat) avail for percolation
  M_FLUX%QPERC_12 = MPARAM%PERCRTE * (TSTATE%WATR_1/MPARAM%MAXWATR_1)**MPARAM%PERCEXP
 CASE(iopt_perc_lower) ! perc defined by moisture content in lower layer (SAC)
  ! (compute lower-zone percolation demand -- multiplier on maximum percolation, then percolation)
  LZ_PD = 1._SP + MPARAM%SACPMLT*(1._SP - TSTATE%WATR_2/MPARAM%MAXWATR_2)**MPARAM%SACPEXP
  M_FLUX%QPERC_12 = DPARAM%QBSAT*LZ_PD * (TSTATE%FREE_1/DPARAM%MAXFREE_1)
  !print *, 'lz_pd = ', LZ_PD, MPARAM%SACPMLT, TSTATE%WATR_2/MPARAM%MAXWATR_2, MPARAM%SACPEXP
  !print *, 'qperc_12 = ', M_FLUX%QPERC_12, DPARAM%QBSAT, LZ_PD, TSTATE%FREE_1/DPARAM%MAXFREE_1
 CASE DEFAULT       ! check for errors
  print *, "SMODL%iQPERC must be iopt_perc_f2sat, iopt_perc_w2sat, or iopt_perc_lower"
  STOP
END SELECT
! ---------------------------------------------------------------------------------------
END SUBROUTINE QPERCOLATE
