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
! Modified by Brian Henn to include snow model, 6/2013
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
USE multibands                                        ! model basin band data
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
 CASE ('MBASE')      ; XVAR = MPARAM%MBASE
 CASE ('MFMAX')      ; XVAR = MPARAM%MFMAX
 CASE ('MFMIN')      ; XVAR = MPARAM%MFMIN
 CASE ('PXTEMP')     ; XVAR = MPARAM%PXTEMP
 CASE ('OPG')        ; XVAR = MPARAM%OPG
 CASE ('LAPSE')      ; XVAR = MPARAM%LAPSE
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
 ! basin band data
 CASE ('Z_MID01')   ; XVAR = MBANDS(1)%Z_MID
 CASE ('AF01')      ; XVAR = MBANDS(1)%AF
 CASE ('Z_MID02')   ; XVAR = MBANDS(2)%Z_MID
 CASE ('AF02')      ; XVAR = MBANDS(2)%AF
 CASE ('Z_MID03')   ; XVAR = MBANDS(3)%Z_MID
 CASE ('AF03')      ; XVAR = MBANDS(3)%AF
 CASE ('Z_MID04')   ; XVAR = MBANDS(4)%Z_MID
 CASE ('AF04')      ; XVAR = MBANDS(4)%AF
 CASE ('Z_MID05')   ; XVAR = MBANDS(5)%Z_MID
 CASE ('AF05')      ; XVAR = MBANDS(5)%AF
 CASE ('Z_MID06')   ; XVAR = MBANDS(6)%Z_MID
 CASE ('AF06')      ; XVAR = MBANDS(6)%AF
 CASE ('Z_MID07')   ; XVAR = MBANDS(7)%Z_MID
 CASE ('AF07')      ; XVAR = MBANDS(7)%AF
 CASE ('Z_MID08')   ; XVAR = MBANDS(8)%Z_MID
 CASE ('AF08')      ; XVAR = MBANDS(8)%AF
 CASE ('Z_MID09')   ; XVAR = MBANDS(9)%Z_MID
 CASE ('AF09')      ; XVAR = MBANDS(9)%AF
 CASE ('Z_MID10')   ; XVAR = MBANDS(10)%Z_MID
 CASE ('AF10')      ; XVAR = MBANDS(10)%AF
 CASE ('Z_MID11')   ; XVAR = MBANDS(11)%Z_MID
 CASE ('AF11')      ; XVAR = MBANDS(11)%AF
 CASE ('Z_MID12')   ; XVAR = MBANDS(12)%Z_MID
 CASE ('AF12')      ; XVAR = MBANDS(12)%AF
 CASE ('Z_MID13')   ; XVAR = MBANDS(13)%Z_MID
 CASE ('AF13')      ; XVAR = MBANDS(13)%AF
 CASE ('Z_MID14')   ; XVAR = MBANDS(14)%Z_MID
 CASE ('AF14')      ; XVAR = MBANDS(14)%AF
 CASE ('Z_MID15')   ; XVAR = MBANDS(15)%Z_MID
 CASE ('AF15')      ; XVAR = MBANDS(15)%AF
 CASE ('Z_MID16')   ; XVAR = MBANDS(16)%Z_MID
 CASE ('AF16')      ; XVAR = MBANDS(16)%AF
 CASE ('Z_MID17')   ; XVAR = MBANDS(17)%Z_MID
 CASE ('AF17')      ; XVAR = MBANDS(17)%AF
 CASE ('Z_MID18')   ; XVAR = MBANDS(18)%Z_MID
 CASE ('AF18')      ; XVAR = MBANDS(18)%AF
 CASE ('Z_MID19')   ; XVAR = MBANDS(19)%Z_MID
 CASE ('AF19')      ; XVAR = MBANDS(19)%AF
 CASE ('Z_MID20')   ; XVAR = MBANDS(20)%Z_MID
 CASE ('AF20')      ; XVAR = MBANDS(20)%AF
 CASE ('Z_MID21')   ; XVAR = MBANDS(21)%Z_MID
 CASE ('AF21')      ; XVAR = MBANDS(21)%AF
 CASE ('Z_MID22')   ; XVAR = MBANDS(22)%Z_MID
 CASE ('AF22')      ; XVAR = MBANDS(22)%AF
 CASE ('Z_MID23')   ; XVAR = MBANDS(23)%Z_MID
 CASE ('AF23')      ; XVAR = MBANDS(23)%AF
 CASE ('Z_MID24')   ; XVAR = MBANDS(24)%Z_MID
 CASE ('AF24')      ; XVAR = MBANDS(24)%AF
 CASE ('Z_MID25')   ; XVAR = MBANDS(25)%Z_MID
 CASE ('AF25')      ; XVAR = MBANDS(25)%AF
 CASE ('Z_MID26')   ; XVAR = MBANDS(26)%Z_MID
 CASE ('AF26')      ; XVAR = MBANDS(26)%AF
 CASE ('Z_MID27')   ; XVAR = MBANDS(27)%Z_MID
 CASE ('AF27')      ; XVAR = MBANDS(27)%AF
 CASE ('Z_MID28')   ; XVAR = MBANDS(28)%Z_MID
 CASE ('AF28')      ; XVAR = MBANDS(28)%AF
 CASE ('Z_MID29')   ; XVAR = MBANDS(29)%Z_MID
 CASE ('AF29')      ; XVAR = MBANDS(29)%AF
 CASE ('Z_MID30')   ; XVAR = MBANDS(30)%Z_MID
 CASE ('AF30')      ; XVAR = MBANDS(30)%AF
 CASE ('Z_MID31')   ; XVAR = MBANDS(31)%Z_MID
 CASE ('AF31')      ; XVAR = MBANDS(31)%AF
 CASE ('Z_MID32')   ; XVAR = MBANDS(32)%Z_MID
 CASE ('AF32')      ; XVAR = MBANDS(32)%AF
 CASE ('Z_MID33')   ; XVAR = MBANDS(33)%Z_MID
 CASE ('AF33')      ; XVAR = MBANDS(33)%AF
 CASE ('Z_MID34')   ; XVAR = MBANDS(34)%Z_MID
 CASE ('AF34')      ; XVAR = MBANDS(34)%AF
 CASE ('Z_MID35')   ; XVAR = MBANDS(35)%Z_MID
 CASE ('AF35')      ; XVAR = MBANDS(35)%AF
 CASE ('Z_MID36')   ; XVAR = MBANDS(36)%Z_MID
 CASE ('AF36')      ; XVAR = MBANDS(36)%AF
 CASE ('Z_MID37')   ; XVAR = MBANDS(37)%Z_MID
 CASE ('AF37')      ; XVAR = MBANDS(37)%AF
 CASE ('Z_MID38')   ; XVAR = MBANDS(38)%Z_MID
 CASE ('AF38')      ; XVAR = MBANDS(38)%AF
 CASE ('Z_MID39')   ; XVAR = MBANDS(39)%Z_MID
 CASE ('AF39')      ; XVAR = MBANDS(39)%AF
 CASE ('Z_MID40')   ; XVAR = MBANDS(40)%Z_MID
 CASE ('AF40')      ; XVAR = MBANDS(40)%AF
 CASE ('Z_MID41')   ; XVAR = MBANDS(41)%Z_MID
 CASE ('AF41')      ; XVAR = MBANDS(41)%AF
 CASE ('Z_MID42')   ; XVAR = MBANDS(42)%Z_MID
 CASE ('AF42')      ; XVAR = MBANDS(42)%AF
 CASE ('Z_MID43')   ; XVAR = MBANDS(43)%Z_MID
 CASE ('AF43')      ; XVAR = MBANDS(43)%AF
 CASE ('Z_MID44')   ; XVAR = MBANDS(44)%Z_MID
 CASE ('AF44')      ; XVAR = MBANDS(44)%AF
 CASE ('Z_MID45')   ; XVAR = MBANDS(45)%Z_MID
 CASE ('AF45')      ; XVAR = MBANDS(45)%AF
 CASE ('Z_MID46')   ; XVAR = MBANDS(46)%Z_MID
 CASE ('AF46')      ; XVAR = MBANDS(46)%AF
 CASE ('Z_MID47')   ; XVAR = MBANDS(47)%Z_MID
 CASE ('AF47')      ; XVAR = MBANDS(47)%AF
 CASE ('Z_MID48')   ; XVAR = MBANDS(48)%Z_MID
 CASE ('AF48')      ; XVAR = MBANDS(48)%AF
 CASE ('Z_MID49')   ; XVAR = MBANDS(49)%Z_MID
 CASE ('AF49')      ; XVAR = MBANDS(49)%AF
 CASE ('Z_MID50')   ; XVAR = MBANDS(50)%Z_MID
 CASE ('AF50')      ; XVAR = MBANDS(50)%AF
 CASE('N_BANDS')    ; XVAR = N_BANDS
 CASE('Z_FORCING')  ; XVAR = Z_FORCING
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
