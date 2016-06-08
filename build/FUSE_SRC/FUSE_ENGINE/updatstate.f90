SUBROUTINE UPDATSTATE(DT)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Use model derivatives to update model states over time interval DT.
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable definitions, etc.
USE model_defn                                        ! model definition structures
USE model_defnames
USE multiparam                                        ! model parameters
USE multistate                                        ! model states
IMPLICIT NONE
! input
REAL(SP), INTENT(IN)                   :: DT          ! length of the time step
! internal
REAL(SP), PARAMETER                    :: XMIN=1.E-06 ! very small number
REAL(SP), PARAMETER                    :: missingValue=-9999._sp
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iARCH1)  ! (upper layer architecture)
 CASE(iopt_onestate_1) ! upper layer defined by a single state variable
  ! (update state)
  FSTATE%WATR_1  = MAX(XMIN*MPARAM%MAXWATR_1, FSTATE%WATR_1  + DY_DT%WATR_1 *DT)
  ! (derive state)
  FSTATE%TENS_1A = missingValue                              ! 1st tension store (undefined)
  FSTATE%TENS_1B = missingValue                              ! 2nd tension store (undefined)
  FSTATE%TENS_1  = MIN(FSTATE%WATR_1, DPARAM%MAXTENS_1)      ! tension storage
  FSTATE%FREE_1  = MAX(XMIN, FSTATE%WATR_1 - DPARAM%MAXTENS_1) ! free storage
 CASE(iopt_tension1_1) ! upper layer broken up into tension and free storage
  ! (update state)
  FSTATE%TENS_1  = MAX(XMIN*DPARAM%MAXTENS_1, FSTATE%TENS_1  + DY_DT%TENS_1 *DT)
  FSTATE%FREE_1  = MAX(XMIN*DPARAM%MAXFREE_1, FSTATE%FREE_1  + DY_DT%FREE_1 *DT)
  ! (derive state)
  FSTATE%TENS_1A = missingValue                              ! 1st tension store (undefined)
  FSTATE%TENS_1B = missingValue                              ! 2nd tension store (undefined)
  FSTATE%WATR_1  = FSTATE%TENS_1 + FSTATE%FREE_1             ! total storage
 CASE(iopt_tension2_1) ! tension storage sub-divided into recharge and excess
  ! (update state)
  FSTATE%TENS_1A = MAX(XMIN*DPARAM%MAXTENS_1A, FSTATE%TENS_1A + DY_DT%TENS_1A*DT)
  FSTATE%TENS_1B = MAX(XMIN*DPARAM%MAXTENS_1B, FSTATE%TENS_1B + DY_DT%TENS_1B*DT)
  FSTATE%FREE_1  = MAX(XMIN*DPARAM%MAXFREE_1,  FSTATE%FREE_1  + DY_DT%FREE_1 *DT)
  ! (derive state)
  FSTATE%TENS_1  = FSTATE%TENS_1A + FSTATE%TENS_1B           ! tension storage
  FSTATE%WATR_1  = FSTATE%TENS_1  + FSTATE%FREE_1            ! total storage
 CASE DEFAULT       ! (error check)
  print *, "MDEFN(IMOD)%iARCH1 must be iopt_tension2_1, iopt_tension1_1, or iopt_onestate_1"
  STOP
END SELECT  ! (upper layer architechure)
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iARCH2)  ! (lower layer architecture)
 CASE(iopt_tens2pll_2) ! tension reservoir plus two parallel tanks
  ! (update state)
  FSTATE%TENS_2  = MAX(XMIN*DPARAM%MAXTENS_2,  FSTATE%TENS_2  + DY_DT%TENS_2 *DT)
  FSTATE%FREE_2A = MAX(XMIN*DPARAM%MAXFREE_2A, FSTATE%FREE_2A + DY_DT%FREE_2A*DT)
  FSTATE%FREE_2B = MAX(XMIN*DPARAM%MAXFREE_2B, FSTATE%FREE_2B + DY_DT%FREE_2B*DT)
  ! (derive state)
  FSTATE%FREE_2  = FSTATE%FREE_2A + FSTATE%FREE_2B              ! free storage
  FSTATE%WATR_2  = FSTATE%TENS_2  + FSTATE%FREE_2               ! total storage
 CASE(iopt_unlimfrc_2,iopt_unlimpow_2,iopt_topmdexp_2,iopt_fixedsiz_2) ! single baseflow reservoir
  ! (update state)
  ! check in more detail -- unlimited states do not need truncation
  IF (SMODL%iARCH2.EQ.iopt_topmdexp_2) THEN
   FSTATE%WATR_2  = FSTATE%WATR_2  + DY_DT%WATR_2 *DT
  ELSE
   FSTATE%WATR_2  = MAX(XMIN*MPARAM%MAXWATR_2, FSTATE%WATR_2  + DY_DT%WATR_2 *DT)
  ENDIF
  ! (derive state)
  FSTATE%TENS_2  = MIN(FSTATE%WATR_2, DPARAM%MAXTENS_2)         ! tension storage
  FSTATE%FREE_2  = MAX(0._sp, FSTATE%WATR_2 - DPARAM%MAXTENS_2) ! free storage
  FSTATE%FREE_2A = missingValue                                 ! primary reservoir (undefined)
  FSTATE%FREE_2B = missingValue                                 ! secondary reservoir (undefined)
 CASE DEFAULT       ! (error check)
  print *, "MDEFN(IMOD)%iARCH2 must be iopt_tens2pll_2, iopt_unlimfrc_2, iopt_unlimpow_2"
  print *, "  iopt_topmdexp_2, or iopt_fixedsiz_2"
  STOP
END SELECT
! ---------------------------------------------------------------------------------------
END SUBROUTINE UPDATSTATE
