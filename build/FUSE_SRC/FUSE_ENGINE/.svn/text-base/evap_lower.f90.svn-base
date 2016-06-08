SUBROUTINE EVAP_LOWER()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Computes evaporation from the lower soil layer
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multi_flux -- evaporation stored in  MODULE multi_flux
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition structure
USE model_defnames
USE multiparam                                        ! model parameters
USE multiforce                                        ! model forcing
USE multistate                                        ! model states
USE multi_flux                                        ! model fluxes
IMPLICIT NONE
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iARCH2)  ! lower layer architecture
 CASE(iopt_tens2pll_2,iopt_fixedsiz_2)
  ! -------------------------------------------------------------------------------------
  SELECT CASE(SMODL%iARCH1)
   ! ------------------------------------------------------------------------------------
   CASE(iopt_tension1_1,iopt_onestate_1) ! lower-layer evap is valid
   ! ------------------------------------------------------------------------------------
   ! use different evaporation schemes for the lower layer
   ! -----------------------------------------------------
   SELECT CASE(SMODL%iESOIL)
    CASE(iopt_sequential)
     M_FLUX%EVAP_2 = (MFORCE%PET-M_FLUX%EVAP_1) * (TSTATE%TENS_2/DPARAM%MAXTENS_2)
    CASE(iopt_rootweight)
     M_FLUX%EVAP_2 = MFORCE%PET * DPARAM%RTFRAC2 * (TSTATE%TENS_2/DPARAM%MAXTENS_2)
    CASE DEFAULT
     print *, "SMODL%iESOIL must be either iopt_sequential or iopt_rootweight"
   END SELECT  ! (evaporation schemes)
   ! ------------------------------------------------------------------------------------
   CASE(iopt_tension2_1)               ! lower-layer evap is zero
    M_FLUX%EVAP_2 = 0._sp
   ! ------------------------------------------------------------------------------------
   CASE DEFAULT
    print *, "SMODL%iARCH1 must be iopt_tension2_1, iopt_tension1_1, or iopt_onestate_1"
    STOP
   ! ------------------------------------------------------------------------------------
  END SELECT  ! (upper-layer architechure)
 ! --------------------------------------------------------------------------------------
 CASE(iopt_unlimfrc_2,iopt_unlimpow_2,iopt_topmdexp_2)
  M_FLUX%EVAP_2 = 0._sp
 ! --------------------------------------------------------------------------------------
 CASE DEFAULT
  print *, "SMODL%iARCH2 must be iopt_tens2pll_2, iopt_unlimfrc_2, iopt_unlimpow_2"
  print *, "  iopt_topmdexp_2, or iopt_fixedsiz_2"
  STOP
END SELECT
! ---------------------------------------------------------------------------------------
END SUBROUTINE EVAP_LOWER
