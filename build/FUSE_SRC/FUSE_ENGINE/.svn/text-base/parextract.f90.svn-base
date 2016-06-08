MODULE PAREXTRACT_MODULE
IMPLICIT NONE
CONTAINS
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
SUBROUTINE GET_PARSET(PARSET)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2008
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Extracts an entire parameter set from a data structure, based on the list of parameters
! in LPARAM
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE multiparam                                        ! model parameters
IMPLICIT NONE
! output
REAL(SP), INTENT(INOUT), DIMENSION(:)  :: PARSET      ! parameter set 
! local
INTEGER(I4B)                           :: IPAR        ! looping
! ---------------------------------------------------------------------------------------
DO IPAR=1,NUMPAR  ! NUMPAR is stored in module multiparam
 PARSET(IPAR) = PAREXTRACT(LPARAM(IPAR)%PARNAME)
END DO
! ---------------------------------------------------------------------------------------
END SUBROUTINE GET_PARSET
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
PURE FUNCTION PAREXTRACT(PARNAME)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Extracts parameter from data structures
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE multiparam                                        ! model parameters
USE model_numerix                                     ! model numerix parameters
IMPLICIT NONE
! input
CHARACTER(*), INTENT(IN)               :: PARNAME     ! parameter name
! internal
REAL(SP)                               :: XVAR        ! variable
! output
REAL(SP)                               :: PAREXTRACT  ! FUNCTION name
! ---------------------------------------------------------------------------------------
SELECT CASE (TRIM(PARNAME))
 ! model parameters
 CASE ('RFERR_ADD')  ; XVAR = MPARAM%RFERR_ADD
 CASE ('RFERR_MLT')  ; XVAR = MPARAM%RFERR_MLT
 CASE ('RFH1_MEAN')  ; XVAR = MPARAM%RFH1_MEAN
 CASE ('RFH2_SDEV')  ; XVAR = MPARAM%RFH2_SDEV
 CASE ('RH1P_MEAN')  ; XVAR = MPARAM%RH1P_MEAN
 CASE ('RH1P_SDEV')  ; XVAR = MPARAM%RH1P_SDEV
 CASE ('RH2P_MEAN')  ; XVAR = MPARAM%RH2P_MEAN
 CASE ('RH2P_SDEV')  ; XVAR = MPARAM%RH2P_SDEV
 CASE ('MAXWATR_1')  ; XVAR = MPARAM%MAXWATR_1
 CASE ('MAXWATR_2')  ; XVAR = MPARAM%MAXWATR_2
 CASE ('FRACTEN')    ; XVAR = MPARAM%FRACTEN
 CASE ('FRCHZNE')    ; XVAR = MPARAM%FRCHZNE
 CASE ('FPRIMQB')    ; XVAR = MPARAM%FPRIMQB
 CASE ('RTFRAC1')    ; XVAR = MPARAM%RTFRAC1
 CASE ('PERCRTE')    ; XVAR = MPARAM%PERCRTE
 CASE ('PERCEXP')    ; XVAR = MPARAM%PERCEXP
 CASE ('SACPMLT')    ; XVAR = MPARAM%SACPMLT
 CASE ('SACPEXP')    ; XVAR = MPARAM%SACPEXP
 CASE ('PERCFRAC')   ; XVAR = MPARAM%PERCFRAC
 CASE ('FRACLOWZ')   ; XVAR = MPARAM%FRACLOWZ
 CASE ('IFLWRTE')    ; XVAR = MPARAM%IFLWRTE
 CASE ('BASERTE')    ; XVAR = MPARAM%BASERTE
 CASE ('QB_POWR')    ; XVAR = MPARAM%QB_POWR
 CASE ('QB_PRMS')    ; XVAR = MPARAM%QB_PRMS
 CASE ('QBRATE_2A')  ; XVAR = MPARAM%QBRATE_2A
 CASE ('QBRATE_2B')  ; XVAR = MPARAM%QBRATE_2B
 CASE ('SAREAMAX')   ; XVAR = MPARAM%SAREAMAX
 CASE ('AXV_BEXP')   ; XVAR = MPARAM%AXV_BEXP
 CASE ('LOGLAMB')    ; XVAR = MPARAM%LOGLAMB
 CASE ('TISHAPE')    ; XVAR = MPARAM%TISHAPE
 CASE ('TIMEDELAY')  ; XVAR = MPARAM%TIMEDELAY
 ! derived parameters
 CASE ('MAXTENS_1')  ; XVAR = DPARAM%MAXTENS_1
 CASE ('MAXTENS_1A') ; XVAR = DPARAM%MAXTENS_1A 
 CASE ('MAXTENS_1B') ; XVAR = DPARAM%MAXTENS_1B
 CASE ('MAXFREE_1')  ; XVAR = DPARAM%MAXFREE_1
 CASE ('MAXTENS_2')  ; XVAR = DPARAM%MAXTENS_2
 CASE ('MAXFREE_2')  ; XVAR = DPARAM%MAXFREE_2
 CASE ('MAXFREE_2A') ; XVAR = DPARAM%MAXFREE_2A
 CASE ('MAXFREE_2B') ; XVAR = DPARAM%MAXFREE_2B
 CASE ('QBSAT')      ; XVAR = DPARAM%QBSAT
 CASE ('RTFRAC2')    ; XVAR = DPARAM%RTFRAC2
 CASE ('POWLAMB')    ; XVAR = DPARAM%POWLAMB
 CASE ('MAXPOW')     ; XVAR = DPARAM%MAXPOW
 ! numerical solution parameters
 CASE ('SOLUTION')   ; XVAR = REAL(SOLUTION_METHOD, KIND(SP))
 CASE ('TIMSTEP_TYP'); XVAR = REAL(TEMPORAL_ERROR_CONTROL, KIND(SP))
 CASE ('INITL_GUESS'); XVAR = REAL(INITIAL_NEWTON, KIND(SP))
 CASE ('JAC_RECOMPT'); XVAR = REAL(JAC_RECOMPUTE, KIND(SP))
 CASE ('CK_OVRSHOOT'); XVAR = REAL(CHECK_OVERSHOOT, KIND(SP)) 
 CASE ('SMALL_ESTEP'); XVAR = REAL(SMALL_ENDSTEP, KIND(SP)) 
 CASE ('ERRTRUNCABS'); XVAR = ERR_TRUNC_ABS
 CASE ('ERRTRUNCREL'); XVAR = ERR_TRUNC_REL
 CASE ('ERRITERFUNC'); XVAR = ERR_ITER_FUNC
 CASE ('ERR_ITER_DX'); XVAR = ERR_ITER_DX
 CASE ('THRESH_FRZE'); XVAR = THRESH_FRZE
 CASE ('FSTATE_MIN') ; XVAR = FRACSTATE_MIN
 CASE ('STEP_SAFETY'); XVAR = SAFETY
 CASE ('RMIN')       ; XVAR = RMIN
 CASE ('RMAX')       ; XVAR = RMAX
 CASE ('NITER_TOTAL'); XVAR = REAL(NITER_TOTAL, KIND(SP))
 CASE ('MIN_TSTEP')  ; XVAR = MIN_TSTEP
 CASE ('MAX_TSTEP')  ; XVAR = MAX_TSTEP
 ! Sobol identifier
 CASE ('SOBOL_INDX') ; XVAR = REAL(SOBOL_INDX, KIND(SP))
END SELECT
! and, save the output
PAREXTRACT = XVAR
! ---------------------------------------------------------------------------------------
END FUNCTION PAREXTRACT
END MODULE PAREXTRACT_MODULE
