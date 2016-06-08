SUBROUTINE QSATEXCESS()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Computes the saturated area and surface runoff
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multi_flux -- saturated area and surface runoff stored in  MODULE multi_flux
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE nr, ONLY : gammp                                  ! interface for the incomplete gamma function
USE model_defn                                        ! model definition structure
USE model_defnames
USE multiparam                                        ! model parameters
USE multiforce                                        ! model forcing
USE multistate                                        ! model states
USE multi_flux                                        ! model fluxes
IMPLICIT NONE
! internal variables
REAL(SP)                               :: TI_SAT      ! topographic index where saturated
REAL(SP)                               :: TI_LOG      ! critical value of topo index in log space
REAL(SP)                               :: TI_OFF      ! offset in the Gamma distribution
REAL(SP)                               :: TI_SHP      ! shape of the Gamma distribution
REAL(SP)                               :: TI_CHI      ! CHI, see Sivapalan et al., 1987
REAL(SP)                               :: TI_ARG      ! argument of the Gamma function
REAL(SP)                               :: NO_ZERO=1.E-8  ! avoid divide by zero
! ---------------------------------------------------------------------------------------
! saturated area method
SELECT CASE(SMODL%iQSURF)
 CASE(iopt_arno_x_vic) ! ARNO/Xzang/VIC parameterization (upper zone control)
  M_FLUX%SATAREA = 1._sp - ( 1._sp - MIN(TSTATE%WATR_1/MPARAM%MAXWATR_1, 1._sp) )**MPARAM%AXV_BEXP
 CASE(iopt_prms_varnt) ! PRMS variant (fraction of upper tension storage)
  M_FLUX%SATAREA = MIN(TSTATE%TENS_1/DPARAM%MAXTENS_1, 1._sp) * MPARAM%SAREAMAX
 CASE(iopt_tmdl_param) ! TOPMODEL parameterization (only valid for TOPMODEL qb)

  ! compute the minimum value of the topographic index where the basin is saturated
  ! (this is correct, as MPARAM%MAXWATR_2 is m*n -- units are meters**(1/n)
  TI_SAT = DPARAM%POWLAMB / (TSTATE%WATR_2/MPARAM%MAXWATR_2 + NO_ZERO)
  ! compute the saturated area
  IF (TI_SAT.GT.DPARAM%MAXPOW) THEN
   M_FLUX%SATAREA = 0.
  ELSE
   ! convert the topographic index to log space
   TI_LOG = LOG( TI_SAT**MPARAM%QB_POWR )   
   ! compute the saturated area (NOTE: critical value of the topographic index is in log space)
   TI_OFF = 3._sp           ! offset in the Gamma distribution (the "3rd" parameter)
   TI_SHP = MPARAM%TISHAPE  ! shape of the Gamma distribution (the "2nd" parameter)
   TI_CHI = (MPARAM%LOGLAMB - TI_OFF) / MPARAM%TISHAPE ! Chi -- loglamb is the first parameter (mean)
   TI_ARG = MAX(0._sp, TI_LOG - TI_OFF) / TI_CHI       ! argument to the incomplete Gamma function
   M_FLUX%SATAREA = 1._sp - GAMMP(TI_SHP, TI_ARG)      ! GAMMP is the incomplete Gamma function
  ENDIF

 ! check processed surface runoff selection 
 CASE DEFAULT
  print *, "SMODL%iQSURF must be iopt_arno_x_vic, iopt_prms_varnt, or iopt_tmdl_param"
  STOP

END SELECT  ! (different surface runoff options)

! ...and, compute surface runoff
! ------------------------------
M_FLUX%QSURF = M_FLUX%EFF_PPT * M_FLUX%SATAREA

END SUBROUTINE QSATEXCESS
