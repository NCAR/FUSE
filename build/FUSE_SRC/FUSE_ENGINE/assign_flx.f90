SUBROUTINE ASSIGN_FLX()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Build an array of strings that list model fluxes used for the current model
!  configuration
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! Defines list and number of states in MODULE model_defn
! ---------------------------------------------------------------------------------------
USE model_defn                                        ! model definition
USE model_defnames
IMPLICIT NONE
INTEGER(I4B) :: I_FLUX ! just used for testing
LOGICAL(LGT) :: L_TEST ! just used for testing
! ---------------------------------------------------------------------------------------
L_TEST=.FALSE.
N_FLUX=0
C_FLUX(:)%FNAME = '           '
! ---------------------------------------------------------------------------------------
! (1) DEFINE STATE VARIABLES IN THE UPPER LAYER
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iARCH1)
 CASE(iopt_tension2_1)
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'EFF_PPT    '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'EVAP_1A    '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'EVAP_1B    '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'RCHR2EXCS  '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'TENS2FREE_1'
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'QPERC_12   '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'QINTF_1    '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'OFLOW_1    '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'QSURF      '
 CASE(iopt_tension1_1)
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'EFF_PPT    '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'EVAP_1     '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'TENS2FREE_1'
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'QPERC_12   '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'QINTF_1    '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'OFLOW_1    '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'QSURF      '
 CASE(iopt_onestate_1)
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'EFF_PPT    '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'EVAP_1     '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'QPERC_12   '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'QINTF_1    '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'OFLOW_1    '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'QSURF      '
 CASE DEFAULT
  print *, "MDEFN(IMOD)%iARCH1 must be iopt_tension2_1, iopt_tension1_1, or iopt_onestate_1"
  STOP
END SELECT
! ---------------------------------------------------------------------------------------
! (2) DEFINE STATE VARIABLES IN THE LOWER LAYER
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iARCH2)
 CASE(iopt_tens2pll_2)
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'EVAP_2     '  
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'TENS2FREE_2'
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'QBASE_2A   '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'QBASE_2B   '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'QBASE_2    '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'OFLOW_2A   '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'OFLOW_2B   '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'OFLOW_2    '
 CASE(iopt_unlimfrc_2,iopt_unlimpow_2,iopt_topmdexp_2,iopt_fixedsiz_2)
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'EVAP_2     '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'QBASE_2    '
  N_FLUX=N_FLUX+1; C_FLUX(N_FLUX)%FNAME = 'OFLOW_2    '
 CASE DEFAULT
  print *, "SMODL%iARCH2 must be iopt_tens2pll_2, iopt_unlimfrc_2, iopt_unlimpow_2"
  print *, "  iopt_topmdexp_2, or iopt_fixedsiz_2"
  STOP
END SELECT
! ---------------------------------------------------------------------------------------
IF (L_TEST) THEN; DO I_FLUX=1,N_FLUX; WRITE(*,'(A20)') C_FLUX(I_FLUX)%FNAME; END DO; ENDIF
! ---------------------------------------------------------------------------------------
END SUBROUTINE ASSIGN_FLX
