SUBROUTINE Q_MISSCELL()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007 (revised 2009 to include a residual method)
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Computes miscellaneous fluxes:
!   RCHR2EXCS   = flow from recharge to excess (mm day-1)
!   TENS2FREE_1 = flow from tension storage to free storage in the upper layer (mm day-1)
!   TENS2FREE_2 = flow from tension storage to free storage in the lower layer (mm day-1)
!   OFLOW_1     = overflow from the upper soil layer (mm day-1)
!   OFLOW_2     = overflow from the lower soil layer (mm day-1)
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multi_flux -- baseflow stored in  MODULE multi_flux
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition structure
USE model_defnames
USE multiparam, ONLY: MPARAM,DPARAM                   ! model parameters
USE multistate, ONLY: MSTATE,TSTATE                   ! model states
USE multi_flux, ONLY: M_FLUX,CURRENT_DT               ! model fluxes
USE model_numerix                                     ! access model numerix decisions
IMPLICIT NONE
REAL(SP)                               :: LOGISMOOTH  ! FUNCTION logistic smoothing
REAL(SP), PARAMETER                    :: PSMOOTH=0.01_SP ! smoothing parameter
REAL(SP)                               :: W_FUNC      ! result from LOGISMOOTH
REAL(SP)                               :: DT          ! current time step
INTEGER(I4B), PARAMETER                :: POP_CASE=9 ! just a temporary fix so the case statement is populated
! ---------------------------------------------------------------------------------------
SELECT CASE(SOLUTION_METHOD)
 CASE (EXPLICIT_EULER,IMPLICIT_EULER,EXPLICIT_HEUN,IMPLICIT_HEUN,SEMI_IMPLICIT)
 ! ---------------------------------------------------------------------------------------
 ! (1) OVERFLOW FLUXES AS A FRACTION OF INFLUXES
 ! ---------------------------------------------------------------------------------------
 SELECT CASE(SMODL%iARCH1)
  CASE(iopt_tension2_1) ! tension storage sub-divided into recharge and excess
   ! compute flow from recharge to excess (mm s-1)
   W_FUNC = LOGISMOOTH(TSTATE%TENS_1A,DPARAM%MAXTENS_1A,PSMOOTH)
   M_FLUX%RCHR2EXCS   = W_FUNC * (M_FLUX%EFF_PPT - M_FLUX%QSURF)
   ! compute flow from tension storage to free storage (mm s-1)
   W_FUNC = LOGISMOOTH(TSTATE%TENS_1B,DPARAM%MAXTENS_1B,PSMOOTH)
   M_FLUX%TENS2FREE_1 = W_FUNC * M_FLUX%RCHR2EXCS
   ! compute over-flow of free water
   W_FUNC = LOGISMOOTH(TSTATE%FREE_1,DPARAM%MAXFREE_1,PSMOOTH)
   M_FLUX%OFLOW_1     = W_FUNC * M_FLUX%TENS2FREE_1
  CASE(iopt_tension1_1) ! upper layer broken up into tension and free storage
   ! no separate recharge zone (flux should never be used)
   M_FLUX%RCHR2EXCS   = 0._SP
   ! compute flow from tension storage to free storage (mm s-1)
   W_FUNC = LOGISMOOTH(TSTATE%TENS_1,DPARAM%MAXTENS_1,PSMOOTH)
   M_FLUX%TENS2FREE_1 = W_FUNC * (M_FLUX%EFF_PPT - M_FLUX%QSURF)
   ! compute over-flow of free water
   W_FUNC = LOGISMOOTH(TSTATE%FREE_1,DPARAM%MAXFREE_1,PSMOOTH)
   M_FLUX%OFLOW_1     = W_FUNC * M_FLUX%TENS2FREE_1
  CASE(iopt_onestate_1) ! upper layer defined by a single state variable
   ! no tension stores
   M_FLUX%RCHR2EXCS   = 0._SP
   M_FLUX%TENS2FREE_1 = 0._SP
   ! compute over-flow of free water
   W_FUNC = LOGISMOOTH(TSTATE%WATR_1,MPARAM%MAXWATR_1,PSMOOTH)
   M_FLUX%OFLOW_1     = W_FUNC * (M_FLUX%EFF_PPT - M_FLUX%QSURF)
  CASE DEFAULT
   print *, "SMODL%iARCH1 must be iopt_tension2_1, iopt_tension1_1, or iopt_onestate_1"
   STOP
 END SELECT
 ! ---------------------------------------------------------------------------------------
 SELECT CASE(SMODL%iARCH2)
  CASE(iopt_tens2pll_2) ! tension reservoir plus two parallel tanks
   ! compute flow from tension storage to free storage (mm s-1)
   W_FUNC = LOGISMOOTH(TSTATE%TENS_2,DPARAM%MAXTENS_2,PSMOOTH)
   M_FLUX%TENS2FREE_2 = W_FUNC * M_FLUX%QPERC_12*(1._SP-MPARAM%PERCFRAC)
   ! compute over-flow of free water in the primary reservoir
   W_FUNC = LOGISMOOTH(TSTATE%FREE_2A,DPARAM%MAXFREE_2A,PSMOOTH)
   M_FLUX%OFLOW_2A    = W_FUNC * (M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP)
   ! compute over-flow of free water in the secondary reservoir
   W_FUNC = LOGISMOOTH(TSTATE%FREE_2B,DPARAM%MAXFREE_2B,PSMOOTH)
   M_FLUX%OFLOW_2B    = W_FUNC * (M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP)
   ! compute total overflow
   M_FLUX%OFLOW_2     = M_FLUX%OFLOW_2A + M_FLUX%OFLOW_2B
  CASE(iopt_fixedsiz_2)
   ! no tension store
   M_FLUX%TENS2FREE_2 = 0._SP
   M_FLUX%OFLOW_2A    = 0._SP
   M_FLUX%OFLOW_2B    = 0._SP
   ! compute over-flow of free water
   W_FUNC = LOGISMOOTH(TSTATE%WATR_2,MPARAM%MAXWATR_2,PSMOOTH)
   M_FLUX%OFLOW_2     = W_FUNC * M_FLUX%QPERC_12
  CASE(iopt_unlimfrc_2,iopt_unlimpow_2,iopt_topmdexp_2) ! unlimited size
   M_FLUX%TENS2FREE_2 = 0._SP
   M_FLUX%OFLOW_2     = 0._SP
   M_FLUX%OFLOW_2A    = 0._SP
   M_FLUX%OFLOW_2B    = 0._SP
  CASE DEFAULT
   print *, "SMODL%iARCH2 must be iopt_tens2pll_2, iopt_unlimfrc_2, iopt_unlimpow_2"
   print *, "  iopt_topmdexp_2, or iopt_fixedsiz_2"
   STOP
 END SELECT
 ! ---------------------------------------------------------------------------------------
 CASE (POP_CASE)
 ! ---------------------------------------------------------------------------------------
 ! (2) OVERFLOW FLUXES COMPUTED AS A RESIDUAL OF AVAILABLE STORAGE
 ! ---------------------------------------------------------------------------------------
 DT = CURRENT_DT
 ! ---------------------------------------------------------------------------------------
 SELECT CASE(SMODL%iARCH1)
  CASE(iopt_tension2_1) ! tension storage sub-divided into recharge and excess
   ! compute flow from recharge to excess (mm s-1)
   M_FLUX%RCHR2EXCS   = MAX(0._SP, (M_FLUX%EFF_PPT - M_FLUX%QSURF) - (DPARAM%MAXTENS_1A - MSTATE%TENS_1A)/DT)
   ! compute flow from tension storage to free storage (mm s-1)
   M_FLUX%TENS2FREE_1 = MAX(0._SP,  M_FLUX%RCHR2EXCS               - (DPARAM%MAXTENS_1B - MSTATE%TENS_1B)/DT)
   ! compute over-flow of free water
   M_FLUX%OFLOW_1     = MAX(0._SP,  M_FLUX%TENS2FREE_1             - (DPARAM%MAXFREE_1  - MSTATE%FREE_1) /DT)
  CASE(iopt_tension1_1) ! upper layer broken up into tension and free storage
   ! no separate recharge zone (flux should never be used)
   M_FLUX%RCHR2EXCS   = 0._SP
   ! compute flow from tension storage to free storage (mm s-1)
   M_FLUX%TENS2FREE_1 = MAX(0._SP, (M_FLUX%EFF_PPT - M_FLUX%QSURF) - (DPARAM%MAXTENS_1 - MSTATE%TENS_1)/DT)
   ! compute over-flow of free water
   M_FLUX%OFLOW_1     = MAX(0._SP,  M_FLUX%TENS2FREE_1             - (DPARAM%MAXFREE_1 - MSTATE%FREE_1)/DT)
  CASE(iopt_onestate_1) ! upper layer defined by a single state variable
   ! no tension stores
   M_FLUX%RCHR2EXCS   = 0._SP
   M_FLUX%TENS2FREE_1 = 0._SP
   ! compute over-flow of free water
   M_FLUX%OFLOW_1     = MAX(0._SP, (M_FLUX%EFF_PPT - M_FLUX%QSURF) - (MPARAM%MAXWATR_1 - MSTATE%WATR_1)/DT)
  CASE DEFAULT
   print *, "SMODL%iARCH1 must be iopt_tension2_1, iopt_tension1_1, or iopt_onestate_1"
   STOP
 END SELECT
 ! ---------------------------------------------------------------------------------------
 SELECT CASE(SMODL%iARCH2)
  CASE(iopt_tens2pll_2) ! tension reservoir plus two parallel tanks
   ! compute flow from tension storage to free storage (mm s-1)
   M_FLUX%TENS2FREE_2 = MAX(0._SP, M_FLUX%QPERC_12*(1._SP-MPARAM%PERCFRAC) - (DPARAM%MAXTENS_2  - MSTATE%TENS_2 )/DT)
   ! compute over-flow of free water in the primary reservoir
   M_FLUX%OFLOW_2A    = MAX(0._SP, (M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP) &
                                       - (DPARAM%MAXFREE_2A - MSTATE%FREE_2A)/DT)
   ! compute over-flow of free water in the secondary reservoir
   M_FLUX%OFLOW_2B    = MAX(0._SP, (M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP) & 
                                       - (DPARAM%MAXFREE_2B - MSTATE%FREE_2B)/DT)
   ! compute total overflow
   M_FLUX%OFLOW_2     = M_FLUX%OFLOW_2A + M_FLUX%OFLOW_2B
  CASE(iopt_fixedsiz_2)
   ! no tension store
   M_FLUX%TENS2FREE_2 = 0._SP
   M_FLUX%OFLOW_2A    = 0._SP
   M_FLUX%OFLOW_2B    = 0._SP
   ! compute over-flow of free water
   M_FLUX%OFLOW_2     = MAX(0._SP, M_FLUX%QPERC_12 - (MPARAM%MAXWATR_2 - MSTATE%WATR_2)/DT)
  CASE(iopt_unlimfrc_2,iopt_unlimpow_2,iopt_topmdexp_2) ! unlimited size
   M_FLUX%TENS2FREE_2 = 0._SP
   M_FLUX%OFLOW_2     = 0._SP
   M_FLUX%OFLOW_2A    = 0._SP
   M_FLUX%OFLOW_2B    = 0._SP
  CASE DEFAULT
   print *, "SMODL%iARCH2 must be iopt_tens2pll_2, iopt_unlimfrc_2, iopt_unlimpow_2"
   print *, "  iopt_topmdexp_2, or iopt_fixedsiz_2"
   STOP
 END SELECT
 ! ---------------------------------------------------------------------------------------
 CASE DEFAULT
  PRINT *, 'fatal error in q_misscell: unknown solution method; solution method must equal '//&
           '0 (explicit_euler), 1 (explicit heun), 2 (implicit_euler), 3 (implicit_heun), or '//&
           '4 (semi_implicit)'
  STOP
END SELECT
END SUBROUTINE Q_MISSCELL
