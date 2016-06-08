SUBROUTINE ASSIGN_STT()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Build an array of strings that list model states used for the current model
!  configuration
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! Defines list and number of states in MODULE model_defn
! ---------------------------------------------------------------------------------------
USE model_defn                                        ! model definition
USE model_defnames
IMPLICIT NONE
! ---------------------------------------------------------------------------------------
NSTATE=0
!CSTATE(:)%SNAME(1:6) = 'NO_USE'
! ---------------------------------------------------------------------------------------
! (1) DEFINE STATE VARIABLES IN THE UPPER LAYER
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iARCH1)
 CASE(iopt_tension2_1)
  CSTATE(NSTATE+1)%iSNAME = iopt_TENS1A
  CSTATE(NSTATE+2)%iSNAME = iopt_TENS1B
  CSTATE(NSTATE+3)%iSNAME = iopt_FREE_1
  NSTATE = NSTATE+3
 CASE(iopt_tension1_1)
  CSTATE(NSTATE+1)%iSNAME = iopt_TENS_1
  CSTATE(NSTATE+2)%iSNAME = iopt_FREE_1
  NSTATE = NSTATE+2
 CASE(iopt_onestate_1)
  CSTATE(NSTATE+1)%iSNAME = iopt_WATR_1
  NSTATE = NSTATE+1
 CASE DEFAULT
  print *, "MDEFN(IMOD)%iARCH1 must be iopt_tension2_1, iopt_tension1_1, or iopt_onestate_1"
  STOP
END SELECT
! ---------------------------------------------------------------------------------------
! (2) DEFINE STATE VARIABLES IN THE LOWER LAYER
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iARCH2)
 CASE(iopt_tens2pll_2)
  CSTATE(NSTATE+1)%iSNAME = iopt_TENS_2
  CSTATE(NSTATE+2)%iSNAME = iopt_FREE2A
  CSTATE(NSTATE+3)%iSNAME = iopt_FREE2B
  NSTATE = NSTATE+3
 CASE(iopt_unlimfrc_2,iopt_unlimpow_2,iopt_topmdexp_2,iopt_fixedsiz_2)  
  CSTATE(NSTATE+1)%iSNAME = iopt_WATR_2
  NSTATE = NSTATE+1
 CASE DEFAULT
  print *, "SMODL%iARCH2 must be iopt_tens2pll_2, iopt_unlimfrc_2, iopt_unlimpow_2"
  print *, "  iopt_topmdexp_2, or iopt_fixedsiz_2"
  STOP
END SELECT
! ---------------------------------------------------------------------------------------
END SUBROUTINE ASSIGN_STT
