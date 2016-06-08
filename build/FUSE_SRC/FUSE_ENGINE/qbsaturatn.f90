SUBROUTINE QBSATURATN()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Computes baseflow at saturation (used in the SAC percolation model)
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multiparam -- baseflow at saturation stored in MODULE multiparam
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition structures
USE model_defnames
USE multiparam                                        ! model parameters
IMPLICIT NONE
REAL(SP)                               :: TOPMDM      ! TOPMODEL "m" parameter
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iARCH2)
 ! --------------------------------------------------------------------------------------
 CASE(iopt_tens2pll_2) ! tension reservoir plus two parallel tanks
  DPARAM%QBSAT = MPARAM%QBRATE_2A*DPARAM%MAXFREE_2A + MPARAM%QBRATE_2B*DPARAM%MAXFREE_2B
 ! --------------------------------------------------------------------------------------
 CASE(iopt_unlimfrc_2) ! baseflow resvr of unlimited size
  DPARAM%QBSAT = MPARAM%QB_PRMS * MPARAM%MAXWATR_2
 ! --------------------------------------------------------------------------------------
 CASE(iopt_unlimpow_2) ! topmodel power-law transmissivity profile
  ! This is a bit tricky.  The capacity of the aquifer is m*n, where m is a scaling
  ! parameter.  We have the capacity, i.e., MPARAM%MAXWATR_2/1000., and need the
  ! TOPMODEL "m" parameter
  TOPMDM = (MPARAM%MAXWATR_2/1000._sp) / MPARAM%QB_POWR ! NOTE: mm --> m
  ! ...and, compute baseflow
  DPARAM%QBSAT = MPARAM%BASERTE * ( TOPMDM / (DPARAM%POWLAMB**MPARAM%QB_POWR) )
 ! --------------------------------------------------------------------------------------
 CASE(iopt_topmdexp_2) ! topmodel exponential transmissivity profile  (NOTE: mm --> m)
  ! for simplicity we use the CAPACITY as the TOPMODEL scaling parameter
  TOPMDM = MPARAM%MAXWATR_2/1000._sp                   ! NOTE: mm --> m
  ! ..., and compute baseflow
  DPARAM%QBSAT = MPARAM%BASERTE * TOPMDM * EXP(-MPARAM%LOGLAMB)
 ! --------------------------------------------------------------------------------------
 CASE(iopt_fixedsiz_2) ! baseflow reservoir of fixed size
  DPARAM%QBSAT = MPARAM%BASERTE
 ! --------------------------------------------------------------------------------------
 CASE DEFAULT
  print *, "SMODL%iARCH2 must be iopt_tens2pll_2, iopt_unlimfrc_2, iopt_unlimpow_2"
  print *, "  iopt_topmdexp_2, or iopt_fixedsiz_2"
  STOP
 ! --------------------------------------------------------------------------------------
END SELECT
! ---------------------------------------------------------------------------------------
END SUBROUTINE QBSATURATN
