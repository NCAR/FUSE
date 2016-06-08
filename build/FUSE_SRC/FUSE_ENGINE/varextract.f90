MODULE VAREXTRACT_MODULE
IMPLICIT NONE
CONTAINS
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
!PURE FUNCTION VAREXTRACT(VARNAME)
FUNCTION VAREXTRACT(VARNAME)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Extracts variable "VARNAME" from relevant data structures
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE metaoutput                                        ! metadata for all model variables
USE multiforce                                        ! model forcing data
USE multistate                                        ! model states
USE multi_flux                                        ! model fluxes
USE multiroute                                        ! routed runoff
USE model_numerix                                     ! model numerix parameters
IMPLICIT NONE
! input
CHARACTER(*), INTENT(IN)               :: VARNAME     ! variable name
! internal
REAL(SP)                               :: XVAR        ! variable
! output
REAL(SP)                               :: VAREXTRACT  ! FUNCTION name
! ---------------------------------------------------------------------------------------
! initialize XVAR
XVAR=-9999._sp
SELECT CASE (TRIM(VARNAME))
 ! extract forcing data
 CASE ('ppt')        ; XVAR = MFORCE%PPT
 CASE ('pet')        ; XVAR = MFORCE%PET
 CASE ('obsq')       ; XVAR = MFORCE%OBSQ
 ! extract model states
 CASE ('tens_1')     ; XVAR = FSTATE%TENS_1
 CASE ('tens_1a')    ; XVAR = FSTATE%TENS_1A
 CASE ('tens_1b')    ; XVAR = FSTATE%TENS_1B
 CASE ('free_1')     ; XVAR = FSTATE%FREE_1
 CASE ('watr_1')     ; XVAR = FSTATE%WATR_1
 CASE ('tens_2')     ; XVAR = FSTATE%TENS_2
 CASE ('free_2')     ; XVAR = FSTATE%FREE_2
 CASE ('free_2a')    ; XVAR = FSTATE%FREE_2A
 CASE ('free_2b')    ; XVAR = FSTATE%FREE_2B
 CASE ('watr_2')     ; XVAR = FSTATE%WATR_2
 ! extract model fluxes
 CASE ('eff_ppt')    ; XVAR = W_FLUX%EFF_PPT
 CASE ('satarea')    ; XVAR = W_FLUX%SATAREA
 CASE ('qsurf')      ; XVAR = W_FLUX%QSURF
 CASE ('evap_1a')    ; XVAR = W_FLUX%EVAP_1A
 CASE ('evap_1b')    ; XVAR = W_FLUX%EVAP_1B
 CASE ('evap_1')     ; XVAR = W_FLUX%EVAP_1
 CASE ('evap_2')     ; XVAR = W_FLUX%EVAP_2
 CASE ('rchr2excs')  ; XVAR = W_FLUX%RCHR2EXCS
 CASE ('tens2free_1'); XVAR = W_FLUX%TENS2FREE_1
 CASE ('oflow_1')    ; XVAR = W_FLUX%OFLOW_1
 CASE ('tens2free_2'); XVAR = W_FLUX%TENS2FREE_2
 CASE ('qintf_1')    ; XVAR = W_FLUX%QINTF_1
 CASE ('qperc_12')   ; XVAR = W_FLUX%QPERC_12
 CASE ('qbase_2')    ; XVAR = W_FLUX%QBASE_2
 CASE ('qbase_2a')   ; XVAR = W_FLUX%QBASE_2A
 CASE ('qbase_2b')   ; XVAR = W_FLUX%QBASE_2B
 CASE ('oflow_2')    ; XVAR = W_FLUX%OFLOW_2
 CASE ('oflow_2a')   ; XVAR = W_FLUX%OFLOW_2A
 CASE ('oflow_2b')   ; XVAR = W_FLUX%OFLOW_2B
 ! extract extrapolation errors
 CASE ('err_tens_1') ; XVAR = W_FLUX%ERR_TENS_1 
 CASE ('err_tens_1a'); XVAR = W_FLUX%ERR_TENS_1A
 CASE ('err_tens_1b'); XVAR = W_FLUX%ERR_TENS_1B
 CASE ('err_free_1') ; XVAR = W_FLUX%ERR_FREE_1
 CASE ('err_watr_1') ; XVAR = W_FLUX%ERR_WATR_1
 CASE ('err_tens_2') ; XVAR = W_FLUX%ERR_TENS_2
 CASE ('err_free_2') ; XVAR = W_FLUX%ERR_FREE_2
 CASE ('err_free_2a'); XVAR = W_FLUX%ERR_FREE_2A
 CASE ('err_free_2b'); XVAR = W_FLUX%ERR_FREE_2B
 CASE ('err_watr_2') ; XVAR = W_FLUX%ERR_WATR_2
 ! time check
 CASE ('chk_time')   ; XVAR = W_FLUX%CHK_TIME
 ! extract model runoff
 CASE ('q_instnt')   ; XVAR = MROUTE%Q_INSTNT
 CASE ('q_routed')   ; XVAR = MROUTE%Q_ROUTED
 ! extract information on numerical solution (shared in MODULE model_numerix)
 CASE ('num_funcs')  ; XVAR = NUM_FUNCS
 CASE ('numjacobian'); XVAR = NUM_JACOBIAN
 CASE ('sub_accept') ; XVAR = NUMSUB_ACCEPT
 CASE ('sub_reject') ; XVAR = NUMSUB_REJECT
 CASE ('sub_noconv') ; XVAR = NUMSUB_NOCONV
 CASE ('max_iterns') ; XVAR = MAXNUM_ITERNS
END SELECT
! and, save the output
VAREXTRACT = XVAR
! ---------------------------------------------------------------------------------------
END FUNCTION VAREXTRACT
END MODULE VAREXTRACT_MODULE
