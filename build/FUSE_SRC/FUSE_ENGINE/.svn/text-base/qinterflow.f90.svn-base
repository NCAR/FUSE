SUBROUTINE QINTERFLOW()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Computes the interflow from free water in the upper soil layer
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multi_flux -- interflow stored in  MODULE multi_flux
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition structure
USE model_defnames
USE multiparam                                        ! model parameters
USE multistate                                        ! model states
USE multi_flux                                        ! model fluxes
IMPLICIT NONE
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iQINTF)
 CASE(iopt_intflwsome) ! interflow
  M_FLUX%QINTF_1 = MPARAM%IFLWRTE * (TSTATE%FREE_1/DPARAM%MAXFREE_1)
 CASE(iopt_intflwnone) ! no interflow
  M_FLUX%QINTF_1 = 0.
 CASE DEFAULT       ! check for errors
  print *, "SMODL%iQINTF must be either iopt_intflwsome or iopt_intflwnone"
  STOP
END SELECT
! ---------------------------------------------------------------------------------------
END SUBROUTINE QINTERFLOW
