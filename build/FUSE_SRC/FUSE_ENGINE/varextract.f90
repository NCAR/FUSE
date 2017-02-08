MODULE VAREXTRACT_MODULE
IMPLICIT NONE
CONTAINS
! ---------------------------------------------------------------------------------------
PURE FUNCTION VAREXTRACT(VARNAME)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! Modified by Brian Henn to include snow model, 6/2013
! Modified by Nans Addor to enable distributed modeling, 9/2016
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
USE multibands                                        ! model snow bands
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
 CASE ('temp')       ; XVAR = MFORCE%TEMP
 CASE ('pet')        ; XVAR = MFORCE%PET
 ! extract response data
 CASE ('obsq')       ; XVAR = valDat%OBSQ
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
 CASE ('swe_z01')    ; XVAR = MBANDS(1)%SWE
 CASE ('swe_z02')    ; XVAR = MBANDS(2)%SWE
 CASE ('swe_z03')    ; XVAR = MBANDS(3)%SWE
 CASE ('swe_z04')    ; XVAR = MBANDS(4)%SWE
 CASE ('swe_z05')    ; XVAR = MBANDS(5)%SWE
 CASE ('swe_z06')    ; XVAR = MBANDS(6)%SWE
 CASE ('swe_z07')    ; XVAR = MBANDS(7)%SWE
 CASE ('swe_z08')    ; XVAR = MBANDS(8)%SWE
 CASE ('swe_z09')    ; XVAR = MBANDS(9)%SWE
 CASE ('swe_z10')    ; XVAR = MBANDS(10)%SWE
 CASE ('swe_z11')    ; XVAR = MBANDS(11)%SWE
 CASE ('swe_z12')    ; XVAR = MBANDS(12)%SWE
 CASE ('swe_z13')    ; XVAR = MBANDS(13)%SWE
 CASE ('swe_z14')    ; XVAR = MBANDS(14)%SWE
 CASE ('swe_z15')    ; XVAR = MBANDS(15)%SWE
 CASE ('swe_z16')    ; XVAR = MBANDS(16)%SWE
 CASE ('swe_z17')    ; XVAR = MBANDS(17)%SWE
 CASE ('swe_z18')    ; XVAR = MBANDS(18)%SWE
 CASE ('swe_z19')    ; XVAR = MBANDS(19)%SWE
 CASE ('swe_z20')    ; XVAR = MBANDS(20)%SWE
 CASE ('swe_z21')    ; XVAR = MBANDS(21)%SWE
 CASE ('swe_z22')    ; XVAR = MBANDS(22)%SWE
 CASE ('swe_z23')    ; XVAR = MBANDS(23)%SWE
 CASE ('swe_z24')    ; XVAR = MBANDS(24)%SWE
 CASE ('swe_z25')    ; XVAR = MBANDS(25)%SWE
 CASE ('swe_z26')    ; XVAR = MBANDS(26)%SWE
 CASE ('swe_z27')    ; XVAR = MBANDS(27)%SWE
 CASE ('swe_z28')    ; XVAR = MBANDS(28)%SWE
 CASE ('swe_z29')    ; XVAR = MBANDS(29)%SWE
 CASE ('swe_z30')    ; XVAR = MBANDS(30)%SWE
 CASE ('swe_z31')    ; XVAR = MBANDS(31)%SWE
 CASE ('swe_z32')    ; XVAR = MBANDS(32)%SWE
 CASE ('swe_z33')    ; XVAR = MBANDS(33)%SWE
 CASE ('swe_z34')    ; XVAR = MBANDS(34)%SWE
 CASE ('swe_z35')    ; XVAR = MBANDS(35)%SWE
 CASE ('swe_z36')    ; XVAR = MBANDS(36)%SWE
 CASE ('swe_z37')    ; XVAR = MBANDS(37)%SWE
 CASE ('swe_z38')    ; XVAR = MBANDS(38)%SWE
 CASE ('swe_z39')    ; XVAR = MBANDS(39)%SWE
 CASE ('swe_z40')    ; XVAR = MBANDS(40)%SWE
 CASE ('swe_z41')    ; XVAR = MBANDS(41)%SWE
 CASE ('swe_z42')    ; XVAR = MBANDS(42)%SWE
 CASE ('swe_z43')    ; XVAR = MBANDS(43)%SWE
 CASE ('swe_z44')    ; XVAR = MBANDS(44)%SWE
 CASE ('swe_z45')    ; XVAR = MBANDS(45)%SWE
 CASE ('swe_z46')    ; XVAR = MBANDS(46)%SWE
 CASE ('swe_z47')    ; XVAR = MBANDS(47)%SWE
 CASE ('swe_z48')    ; XVAR = MBANDS(48)%SWE
 CASE ('swe_z49')    ; XVAR = MBANDS(49)%SWE
 CASE ('swe_z50')    ; XVAR = MBANDS(50)%SWE
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
 CASE ('snwacml_z01'); XVAR = MBANDS(1)%SNOWACCMLTN
 CASE ('snwacml_z02'); XVAR = MBANDS(2)%SNOWACCMLTN
 CASE ('snwacml_z03'); XVAR = MBANDS(3)%SNOWACCMLTN
 CASE ('snwacml_z04'); XVAR = MBANDS(4)%SNOWACCMLTN
 CASE ('snwacml_z05'); XVAR = MBANDS(5)%SNOWACCMLTN
 CASE ('snwacml_z06'); XVAR = MBANDS(6)%SNOWACCMLTN
 CASE ('snwacml_z07'); XVAR = MBANDS(7)%SNOWACCMLTN
 CASE ('snwacml_z08'); XVAR = MBANDS(8)%SNOWACCMLTN
 CASE ('snwacml_z09'); XVAR = MBANDS(9)%SNOWACCMLTN
 CASE ('snwacml_z10'); XVAR = MBANDS(10)%SNOWACCMLTN
 CASE ('snwacml_z11'); XVAR = MBANDS(11)%SNOWACCMLTN
 CASE ('snwacml_z12'); XVAR = MBANDS(12)%SNOWACCMLTN
 CASE ('snwacml_z13'); XVAR = MBANDS(13)%SNOWACCMLTN
 CASE ('snwacml_z14'); XVAR = MBANDS(14)%SNOWACCMLTN
 CASE ('snwacml_z15'); XVAR = MBANDS(15)%SNOWACCMLTN
 CASE ('snwacml_z16'); XVAR = MBANDS(16)%SNOWACCMLTN
 CASE ('snwacml_z17'); XVAR = MBANDS(17)%SNOWACCMLTN
 CASE ('snwacml_z18'); XVAR = MBANDS(18)%SNOWACCMLTN
 CASE ('snwacml_z19'); XVAR = MBANDS(19)%SNOWACCMLTN
 CASE ('snwacml_z20'); XVAR = MBANDS(20)%SNOWACCMLTN
 CASE ('snwacml_z21'); XVAR = MBANDS(21)%SNOWACCMLTN
 CASE ('snwacml_z22'); XVAR = MBANDS(22)%SNOWACCMLTN
 CASE ('snwacml_z23'); XVAR = MBANDS(23)%SNOWACCMLTN
 CASE ('snwacml_z24'); XVAR = MBANDS(24)%SNOWACCMLTN
 CASE ('snwacml_z25'); XVAR = MBANDS(25)%SNOWACCMLTN
 CASE ('snwacml_z26'); XVAR = MBANDS(26)%SNOWACCMLTN
 CASE ('snwacml_z27'); XVAR = MBANDS(27)%SNOWACCMLTN
 CASE ('snwacml_z28'); XVAR = MBANDS(28)%SNOWACCMLTN
 CASE ('snwacml_z29'); XVAR = MBANDS(29)%SNOWACCMLTN
 CASE ('snwacml_z30'); XVAR = MBANDS(30)%SNOWACCMLTN
 CASE ('snwacml_z31'); XVAR = MBANDS(31)%SNOWACCMLTN
 CASE ('snwacml_z32'); XVAR = MBANDS(32)%SNOWACCMLTN
 CASE ('snwacml_z33'); XVAR = MBANDS(33)%SNOWACCMLTN
 CASE ('snwacml_z34'); XVAR = MBANDS(34)%SNOWACCMLTN
 CASE ('snwacml_z35'); XVAR = MBANDS(35)%SNOWACCMLTN
 CASE ('snwacml_z36'); XVAR = MBANDS(36)%SNOWACCMLTN
 CASE ('snwacml_z37'); XVAR = MBANDS(37)%SNOWACCMLTN
 CASE ('snwacml_z38'); XVAR = MBANDS(38)%SNOWACCMLTN
 CASE ('snwacml_z39'); XVAR = MBANDS(39)%SNOWACCMLTN
 CASE ('snwacml_z40'); XVAR = MBANDS(40)%SNOWACCMLTN
 CASE ('snwacml_z41'); XVAR = MBANDS(41)%SNOWACCMLTN
 CASE ('snwacml_z42'); XVAR = MBANDS(42)%SNOWACCMLTN
 CASE ('snwacml_z43'); XVAR = MBANDS(43)%SNOWACCMLTN
 CASE ('snwacml_z44'); XVAR = MBANDS(44)%SNOWACCMLTN
 CASE ('snwacml_z45'); XVAR = MBANDS(45)%SNOWACCMLTN
 CASE ('snwacml_z46'); XVAR = MBANDS(46)%SNOWACCMLTN
 CASE ('snwacml_z47'); XVAR = MBANDS(47)%SNOWACCMLTN
 CASE ('snwacml_z48'); XVAR = MBANDS(48)%SNOWACCMLTN
 CASE ('snwacml_z49'); XVAR = MBANDS(49)%SNOWACCMLTN
 CASE ('snwacml_z50'); XVAR = MBANDS(50)%SNOWACCMLTN
 CASE ('snwmelt_z01'); XVAR = MBANDS(1)%SNOWMELT
 CASE ('snwmelt_z02'); XVAR = MBANDS(2)%SNOWMELT
 CASE ('snwmelt_z03'); XVAR = MBANDS(3)%SNOWMELT
 CASE ('snwmelt_z04'); XVAR = MBANDS(4)%SNOWMELT
 CASE ('snwmelt_z05'); XVAR = MBANDS(5)%SNOWMELT
 CASE ('snwmelt_z06'); XVAR = MBANDS(6)%SNOWMELT
 CASE ('snwmelt_z07'); XVAR = MBANDS(7)%SNOWMELT
 CASE ('snwmelt_z08'); XVAR = MBANDS(8)%SNOWMELT
 CASE ('snwmelt_z09'); XVAR = MBANDS(9)%SNOWMELT
 CASE ('snwmelt_z10'); XVAR = MBANDS(10)%SNOWMELT
 CASE ('snwmelt_z11'); XVAR = MBANDS(11)%SNOWMELT
 CASE ('snwmelt_z12'); XVAR = MBANDS(12)%SNOWMELT
 CASE ('snwmelt_z13'); XVAR = MBANDS(13)%SNOWMELT
 CASE ('snwmelt_z14'); XVAR = MBANDS(14)%SNOWMELT
 CASE ('snwmelt_z15'); XVAR = MBANDS(15)%SNOWMELT
 CASE ('snwmelt_z16'); XVAR = MBANDS(16)%SNOWMELT
 CASE ('snwmelt_z17'); XVAR = MBANDS(17)%SNOWMELT
 CASE ('snwmelt_z18'); XVAR = MBANDS(18)%SNOWMELT
 CASE ('snwmelt_z19'); XVAR = MBANDS(19)%SNOWMELT
 CASE ('snwmelt_z20'); XVAR = MBANDS(20)%SNOWMELT
 CASE ('snwmelt_z21'); XVAR = MBANDS(21)%SNOWMELT
 CASE ('snwmelt_z22'); XVAR = MBANDS(22)%SNOWMELT
 CASE ('snwmelt_z23'); XVAR = MBANDS(23)%SNOWMELT
 CASE ('snwmelt_z24'); XVAR = MBANDS(24)%SNOWMELT
 CASE ('snwmelt_z25'); XVAR = MBANDS(25)%SNOWMELT
 CASE ('snwmelt_z26'); XVAR = MBANDS(26)%SNOWMELT
 CASE ('snwmelt_z27'); XVAR = MBANDS(27)%SNOWMELT
 CASE ('snwmelt_z28'); XVAR = MBANDS(28)%SNOWMELT
 CASE ('snwmelt_z29'); XVAR = MBANDS(29)%SNOWMELT
 CASE ('snwmelt_z30'); XVAR = MBANDS(30)%SNOWMELT
 CASE ('snwmelt_z31'); XVAR = MBANDS(31)%SNOWMELT
 CASE ('snwmelt_z32'); XVAR = MBANDS(32)%SNOWMELT
 CASE ('snwmelt_z33'); XVAR = MBANDS(33)%SNOWMELT
 CASE ('snwmelt_z34'); XVAR = MBANDS(34)%SNOWMELT
 CASE ('snwmelt_z35'); XVAR = MBANDS(35)%SNOWMELT
 CASE ('snwmelt_z36'); XVAR = MBANDS(36)%SNOWMELT
 CASE ('snwmelt_z37'); XVAR = MBANDS(37)%SNOWMELT
 CASE ('snwmelt_z38'); XVAR = MBANDS(38)%SNOWMELT
 CASE ('snwmelt_z39'); XVAR = MBANDS(39)%SNOWMELT
 CASE ('snwmelt_z40'); XVAR = MBANDS(40)%SNOWMELT
 CASE ('snwmelt_z41'); XVAR = MBANDS(41)%SNOWMELT
 CASE ('snwmelt_z42'); XVAR = MBANDS(42)%SNOWMELT
 CASE ('snwmelt_z43'); XVAR = MBANDS(43)%SNOWMELT
 CASE ('snwmelt_z44'); XVAR = MBANDS(44)%SNOWMELT
 CASE ('snwmelt_z45'); XVAR = MBANDS(45)%SNOWMELT
 CASE ('snwmelt_z46'); XVAR = MBANDS(46)%SNOWMELT
 CASE ('snwmelt_z47'); XVAR = MBANDS(47)%SNOWMELT
 CASE ('snwmelt_z48'); XVAR = MBANDS(48)%SNOWMELT
 CASE ('snwmelt_z49'); XVAR = MBANDS(49)%SNOWMELT
 CASE ('snwmelt_z50'); XVAR = MBANDS(50)%SNOWMELT
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

! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
PURE FUNCTION VAREXTRACT_3d(VARNAME)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Nans Addor, based on Martyn Clark's 2007 VAREXTRACT
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
USE multibands                                        ! model snow bands
USE multiroute                                        ! routed runoff
USE model_numerix                                     ! model numerix parameters
IMPLICIT NONE
! input
CHARACTER(*), INTENT(IN)                :: VARNAME     ! variable name
! internal
real(sp),dimension(:,:,:),allocatable   :: XVAR_3d       ! variable
integer(i4b)                            :: ierr        ! error code
CHARACTER(LEN=1024)                     :: MESSAGE         ! error message
! output
real(sp),dimension(:,:,:),allocatable   :: VAREXTRACT_3d  ! FUNCTION name
! ---------------------------------------------------------------------------------------
! initialize XVAR_3d
allocate(XVAR_3d(nSpat1,nSpat2,NUMTIM), VAREXTRACT_3d(nSpat1,nSpat2,NUMTIM), stat=ierr)
if(ierr/=0)then; message=trim(message)//'problem allocating space for XVAR_3d'; return; endif

!XVAR=-9999._sp
SELECT CASE (TRIM(VARNAME))
 ! extract forcing data
 CASE ('ppt')        ; XVAR_3d = gForce_3d%PPT
 CASE ('temp')       ; XVAR_3d = gForce_3d%TEMP
 CASE ('pet')        ; XVAR_3d = gForce_3d%PET
 ! extract response data
 CASE ('obsq')       ; XVAR_3d = valDat%OBSQ
 ! extract model states
 CASE ('tens_1')     ; XVAR_3d = gState_3d%TENS_1
 CASE ('tens_1a')    ; XVAR_3d = gState_3d%TENS_1A
 CASE ('tens_1b')    ; XVAR_3d = gState_3d%TENS_1B
 CASE ('free_1')     ; XVAR_3d = gState_3d%FREE_1
 CASE ('watr_1')     ; XVAR_3d = gState_3d%WATR_1
 CASE ('tens_2')     ; XVAR_3d = gState_3d%TENS_2
 CASE ('free_2')     ; XVAR_3d = gState_3d%FREE_2
 CASE ('free_2a')    ; XVAR_3d = gState_3d%FREE_2A
 CASE ('free_2b')    ; XVAR_3d = gState_3d%FREE_2B
 CASE ('watr_2')     ; XVAR_3d = gState_3d%WATR_2
 CASE ('swe_z01')    ; XVAR_3d = MBANDS(1)%SWE
 CASE ('swe_z02')    ; XVAR_3d = MBANDS(2)%SWE
 CASE ('swe_z03')    ; XVAR_3d = MBANDS(3)%SWE
 CASE ('swe_z04')    ; XVAR_3d = MBANDS(4)%SWE
 CASE ('swe_z05')    ; XVAR_3d = MBANDS(5)%SWE
 CASE ('swe_z06')    ; XVAR_3d = MBANDS(6)%SWE
 CASE ('swe_z07')    ; XVAR_3d = MBANDS(7)%SWE
 CASE ('swe_z08')    ; XVAR_3d = MBANDS(8)%SWE
 CASE ('swe_z09')    ; XVAR_3d = MBANDS(9)%SWE
 CASE ('swe_z10')    ; XVAR_3d = MBANDS(10)%SWE
 CASE ('swe_z11')    ; XVAR_3d = MBANDS(11)%SWE
 CASE ('swe_z12')    ; XVAR_3d = MBANDS(12)%SWE
 CASE ('swe_z13')    ; XVAR_3d = MBANDS(13)%SWE
 CASE ('swe_z14')    ; XVAR_3d = MBANDS(14)%SWE
 CASE ('swe_z15')    ; XVAR_3d = MBANDS(15)%SWE
 CASE ('swe_z16')    ; XVAR_3d = MBANDS(16)%SWE
 CASE ('swe_z17')    ; XVAR_3d = MBANDS(17)%SWE
 CASE ('swe_z18')    ; XVAR_3d = MBANDS(18)%SWE
 CASE ('swe_z19')    ; XVAR_3d = MBANDS(19)%SWE
 CASE ('swe_z20')    ; XVAR_3d = MBANDS(20)%SWE
 CASE ('swe_z21')    ; XVAR_3d = MBANDS(21)%SWE
 CASE ('swe_z22')    ; XVAR_3d = MBANDS(22)%SWE
 CASE ('swe_z23')    ; XVAR_3d = MBANDS(23)%SWE
 CASE ('swe_z24')    ; XVAR_3d = MBANDS(24)%SWE
 CASE ('swe_z25')    ; XVAR_3d = MBANDS(25)%SWE
 CASE ('swe_z26')    ; XVAR_3d = MBANDS(26)%SWE
 CASE ('swe_z27')    ; XVAR_3d = MBANDS(27)%SWE
 CASE ('swe_z28')    ; XVAR_3d = MBANDS(28)%SWE
 CASE ('swe_z29')    ; XVAR_3d = MBANDS(29)%SWE
 CASE ('swe_z30')    ; XVAR_3d = MBANDS(30)%SWE
 CASE ('swe_z31')    ; XVAR_3d = MBANDS(31)%SWE
 CASE ('swe_z32')    ; XVAR_3d = MBANDS(32)%SWE
 CASE ('swe_z33')    ; XVAR_3d = MBANDS(33)%SWE
 CASE ('swe_z34')    ; XVAR_3d = MBANDS(34)%SWE
 CASE ('swe_z35')    ; XVAR_3d = MBANDS(35)%SWE
 CASE ('swe_z36')    ; XVAR_3d = MBANDS(36)%SWE
 CASE ('swe_z37')    ; XVAR_3d = MBANDS(37)%SWE
 CASE ('swe_z38')    ; XVAR_3d = MBANDS(38)%SWE
 CASE ('swe_z39')    ; XVAR_3d = MBANDS(39)%SWE
 CASE ('swe_z40')    ; XVAR_3d = MBANDS(40)%SWE
 CASE ('swe_z41')    ; XVAR_3d = MBANDS(41)%SWE
 CASE ('swe_z42')    ; XVAR_3d = MBANDS(42)%SWE
 CASE ('swe_z43')    ; XVAR_3d = MBANDS(43)%SWE
 CASE ('swe_z44')    ; XVAR_3d = MBANDS(44)%SWE
 CASE ('swe_z45')    ; XVAR_3d = MBANDS(45)%SWE
 CASE ('swe_z46')    ; XVAR_3d = MBANDS(46)%SWE
 CASE ('swe_z47')    ; XVAR_3d = MBANDS(47)%SWE
 CASE ('swe_z48')    ; XVAR_3d = MBANDS(48)%SWE
 CASE ('swe_z49')    ; XVAR_3d = MBANDS(49)%SWE
 CASE ('swe_z50')    ; XVAR_3d = MBANDS(50)%SWE
 ! extract model fluxes
 CASE ('eff_ppt')    ; XVAR_3d = W_FLUX_3d%EFF_PPT
 CASE ('satarea')    ; XVAR_3d = W_FLUX_3d%SATAREA
 CASE ('qsurf')      ; XVAR_3d = W_FLUX_3d%QSURF
 CASE ('evap_1a')    ; XVAR_3d = W_FLUX_3d%EVAP_1A
 CASE ('evap_1b')    ; XVAR_3d = W_FLUX_3d%EVAP_1B
 CASE ('evap_1')     ; XVAR_3d = W_FLUX_3d%EVAP_1
 CASE ('evap_2')     ; XVAR_3d = W_FLUX_3d%EVAP_2
 CASE ('rchr2excs')  ; XVAR_3d = W_FLUX_3d%RCHR2EXCS
 CASE ('tens2free_1'); XVAR_3d = W_FLUX_3d%TENS2FREE_1
 CASE ('oflow_1')    ; XVAR_3d = W_FLUX_3d%OFLOW_1
 CASE ('tens2free_2'); XVAR_3d = W_FLUX_3d%TENS2FREE_2
 CASE ('qintf_1')    ; XVAR_3d = W_FLUX_3d%QINTF_1
 CASE ('qperc_12')   ; XVAR_3d = W_FLUX_3d%QPERC_12
 CASE ('qbase_2')    ; XVAR_3d = W_FLUX_3d%QBASE_2
 CASE ('qbase_2a')   ; XVAR_3d = W_FLUX_3d%QBASE_2A
 CASE ('qbase_2b')   ; XVAR_3d = W_FLUX_3d%QBASE_2B
 CASE ('oflow_2')    ; XVAR_3d = W_FLUX_3d%OFLOW_2
 CASE ('oflow_2a')   ; XVAR_3d = W_FLUX_3d%OFLOW_2A
 CASE ('oflow_2b')   ; XVAR_3d = W_FLUX_3d%OFLOW_2B
CASE ('snwacml_z01'); XVAR_3d = MBANDS(1)%SNOWACCMLTN
 CASE ('snwacml_z02'); XVAR_3d = MBANDS(2)%SNOWACCMLTN
 CASE ('snwacml_z03'); XVAR_3d = MBANDS(3)%SNOWACCMLTN
 CASE ('snwacml_z04'); XVAR_3d = MBANDS(4)%SNOWACCMLTN
 CASE ('snwacml_z05'); XVAR_3d = MBANDS(5)%SNOWACCMLTN
 CASE ('snwacml_z06'); XVAR_3d = MBANDS(6)%SNOWACCMLTN
 CASE ('snwacml_z07'); XVAR_3d = MBANDS(7)%SNOWACCMLTN
 CASE ('snwacml_z08'); XVAR_3d = MBANDS(8)%SNOWACCMLTN
 CASE ('snwacml_z09'); XVAR_3d = MBANDS(9)%SNOWACCMLTN
 CASE ('snwacml_z10'); XVAR_3d = MBANDS(10)%SNOWACCMLTN
 CASE ('snwacml_z11'); XVAR_3d = MBANDS(11)%SNOWACCMLTN
 CASE ('snwacml_z12'); XVAR_3d = MBANDS(12)%SNOWACCMLTN
 CASE ('snwacml_z13'); XVAR_3d = MBANDS(13)%SNOWACCMLTN
 CASE ('snwacml_z14'); XVAR_3d = MBANDS(14)%SNOWACCMLTN
 CASE ('snwacml_z15'); XVAR_3d = MBANDS(15)%SNOWACCMLTN
 CASE ('snwacml_z16'); XVAR_3d = MBANDS(16)%SNOWACCMLTN
 CASE ('snwacml_z17'); XVAR_3d = MBANDS(17)%SNOWACCMLTN
 CASE ('snwacml_z18'); XVAR_3d = MBANDS(18)%SNOWACCMLTN
 CASE ('snwacml_z19'); XVAR_3d = MBANDS(19)%SNOWACCMLTN
 CASE ('snwacml_z20'); XVAR_3d = MBANDS(20)%SNOWACCMLTN
 CASE ('snwacml_z21'); XVAR_3d = MBANDS(21)%SNOWACCMLTN
 CASE ('snwacml_z22'); XVAR_3d = MBANDS(22)%SNOWACCMLTN
 CASE ('snwacml_z23'); XVAR_3d = MBANDS(23)%SNOWACCMLTN
 CASE ('snwacml_z24'); XVAR_3d = MBANDS(24)%SNOWACCMLTN
 CASE ('snwacml_z25'); XVAR_3d = MBANDS(25)%SNOWACCMLTN
 CASE ('snwacml_z26'); XVAR_3d = MBANDS(26)%SNOWACCMLTN
 CASE ('snwacml_z27'); XVAR_3d = MBANDS(27)%SNOWACCMLTN
 CASE ('snwacml_z28'); XVAR_3d = MBANDS(28)%SNOWACCMLTN
 CASE ('snwacml_z29'); XVAR_3d = MBANDS(29)%SNOWACCMLTN
 CASE ('snwacml_z30'); XVAR_3d = MBANDS(30)%SNOWACCMLTN
 CASE ('snwacml_z31'); XVAR_3d = MBANDS(31)%SNOWACCMLTN
 CASE ('snwacml_z32'); XVAR_3d = MBANDS(32)%SNOWACCMLTN
 CASE ('snwacml_z33'); XVAR_3d = MBANDS(33)%SNOWACCMLTN
 CASE ('snwacml_z34'); XVAR_3d = MBANDS(34)%SNOWACCMLTN
 CASE ('snwacml_z35'); XVAR_3d = MBANDS(35)%SNOWACCMLTN
 CASE ('snwacml_z36'); XVAR_3d = MBANDS(36)%SNOWACCMLTN
 CASE ('snwacml_z37'); XVAR_3d = MBANDS(37)%SNOWACCMLTN
 CASE ('snwacml_z38'); XVAR_3d = MBANDS(38)%SNOWACCMLTN
 CASE ('snwacml_z39'); XVAR_3d = MBANDS(39)%SNOWACCMLTN
 CASE ('snwacml_z40'); XVAR_3d = MBANDS(40)%SNOWACCMLTN
 CASE ('snwacml_z41'); XVAR_3d = MBANDS(41)%SNOWACCMLTN
 CASE ('snwacml_z42'); XVAR_3d = MBANDS(42)%SNOWACCMLTN
 CASE ('snwacml_z43'); XVAR_3d = MBANDS(43)%SNOWACCMLTN
 CASE ('snwacml_z44'); XVAR_3d = MBANDS(44)%SNOWACCMLTN
 CASE ('snwacml_z45'); XVAR_3d = MBANDS(45)%SNOWACCMLTN
 CASE ('snwacml_z46'); XVAR_3d = MBANDS(46)%SNOWACCMLTN
 CASE ('snwacml_z47'); XVAR_3d = MBANDS(47)%SNOWACCMLTN
 CASE ('snwacml_z48'); XVAR_3d = MBANDS(48)%SNOWACCMLTN
 CASE ('snwacml_z49'); XVAR_3d = MBANDS(49)%SNOWACCMLTN
 CASE ('snwacml_z50'); XVAR_3d = MBANDS(50)%SNOWACCMLTN
 CASE ('snwmelt_z01'); XVAR_3d = MBANDS(1)%SNOWMELT
 CASE ('snwmelt_z02'); XVAR_3d = MBANDS(2)%SNOWMELT
 CASE ('snwmelt_z03'); XVAR_3d = MBANDS(3)%SNOWMELT
 CASE ('snwmelt_z04'); XVAR_3d = MBANDS(4)%SNOWMELT
 CASE ('snwmelt_z05'); XVAR_3d = MBANDS(5)%SNOWMELT
 CASE ('snwmelt_z06'); XVAR_3d = MBANDS(6)%SNOWMELT
 CASE ('snwmelt_z07'); XVAR_3d = MBANDS(7)%SNOWMELT
 CASE ('snwmelt_z08'); XVAR_3d = MBANDS(8)%SNOWMELT
 CASE ('snwmelt_z09'); XVAR_3d = MBANDS(9)%SNOWMELT
 CASE ('snwmelt_z10'); XVAR_3d = MBANDS(10)%SNOWMELT
 CASE ('snwmelt_z11'); XVAR_3d = MBANDS(11)%SNOWMELT
 CASE ('snwmelt_z12'); XVAR_3d = MBANDS(12)%SNOWMELT
 CASE ('snwmelt_z13'); XVAR_3d = MBANDS(13)%SNOWMELT
 CASE ('snwmelt_z14'); XVAR_3d = MBANDS(14)%SNOWMELT
 CASE ('snwmelt_z15'); XVAR_3d = MBANDS(15)%SNOWMELT
 CASE ('snwmelt_z16'); XVAR_3d = MBANDS(16)%SNOWMELT
 CASE ('snwmelt_z17'); XVAR_3d = MBANDS(17)%SNOWMELT
 CASE ('snwmelt_z18'); XVAR_3d = MBANDS(18)%SNOWMELT
 CASE ('snwmelt_z19'); XVAR_3d = MBANDS(19)%SNOWMELT
 CASE ('snwmelt_z20'); XVAR_3d = MBANDS(20)%SNOWMELT
 CASE ('snwmelt_z21'); XVAR_3d = MBANDS(21)%SNOWMELT
 CASE ('snwmelt_z22'); XVAR_3d = MBANDS(22)%SNOWMELT
 CASE ('snwmelt_z23'); XVAR_3d = MBANDS(23)%SNOWMELT
 CASE ('snwmelt_z24'); XVAR_3d = MBANDS(24)%SNOWMELT
 CASE ('snwmelt_z25'); XVAR_3d = MBANDS(25)%SNOWMELT
 CASE ('snwmelt_z26'); XVAR_3d = MBANDS(26)%SNOWMELT
 CASE ('snwmelt_z27'); XVAR_3d = MBANDS(27)%SNOWMELT
 CASE ('snwmelt_z28'); XVAR_3d = MBANDS(28)%SNOWMELT
 CASE ('snwmelt_z29'); XVAR_3d = MBANDS(29)%SNOWMELT
 CASE ('snwmelt_z30'); XVAR_3d = MBANDS(30)%SNOWMELT
 CASE ('snwmelt_z31'); XVAR_3d = MBANDS(31)%SNOWMELT
 CASE ('snwmelt_z32'); XVAR_3d = MBANDS(32)%SNOWMELT
 CASE ('snwmelt_z33'); XVAR_3d = MBANDS(33)%SNOWMELT
 CASE ('snwmelt_z34'); XVAR_3d = MBANDS(34)%SNOWMELT
 CASE ('snwmelt_z35'); XVAR_3d = MBANDS(35)%SNOWMELT
 CASE ('snwmelt_z36'); XVAR_3d = MBANDS(36)%SNOWMELT
 CASE ('snwmelt_z37'); XVAR_3d = MBANDS(37)%SNOWMELT
 CASE ('snwmelt_z38'); XVAR_3d = MBANDS(38)%SNOWMELT
 CASE ('snwmelt_z39'); XVAR_3d = MBANDS(39)%SNOWMELT
 CASE ('snwmelt_z40'); XVAR_3d = MBANDS(40)%SNOWMELT
 CASE ('snwmelt_z41'); XVAR_3d = MBANDS(41)%SNOWMELT
 CASE ('snwmelt_z42'); XVAR_3d = MBANDS(42)%SNOWMELT
 CASE ('snwmelt_z43'); XVAR_3d = MBANDS(43)%SNOWMELT
 CASE ('snwmelt_z44'); XVAR_3d = MBANDS(44)%SNOWMELT
 CASE ('snwmelt_z45'); XVAR_3d = MBANDS(45)%SNOWMELT
 CASE ('snwmelt_z46'); XVAR_3d = MBANDS(46)%SNOWMELT
 CASE ('snwmelt_z47'); XVAR_3d = MBANDS(47)%SNOWMELT
 CASE ('snwmelt_z48'); XVAR_3d = MBANDS(48)%SNOWMELT
 CASE ('snwmelt_z49'); XVAR_3d = MBANDS(49)%SNOWMELT
 CASE ('snwmelt_z50'); XVAR_3d = MBANDS(50)%SNOWMELT
 ! extract extrapolation errors
 CASE ('err_tens_1') ; XVAR_3d = W_FLUX_3d%ERR_TENS_1
 CASE ('err_tens_1a'); XVAR_3d = W_FLUX_3d%ERR_TENS_1A
 CASE ('err_tens_1b'); XVAR_3d = W_FLUX_3d%ERR_TENS_1B
 CASE ('err_free_1') ; XVAR_3d = W_FLUX_3d%ERR_FREE_1
 CASE ('err_watr_1') ; XVAR_3d = W_FLUX_3d%ERR_WATR_1
 CASE ('err_tens_2') ; XVAR_3d = W_FLUX_3d%ERR_TENS_2
 CASE ('err_free_2') ; XVAR_3d = W_FLUX_3d%ERR_FREE_2
 CASE ('err_free_2a'); XVAR_3d = W_FLUX_3d%ERR_FREE_2A
 CASE ('err_free_2b'); XVAR_3d = W_FLUX_3d%ERR_FREE_2B
 CASE ('err_watr_2') ; XVAR_3d = W_FLUX_3d%ERR_WATR_2
 ! time check
 CASE ('chk_time')   ; XVAR_3d = W_FLUX_3d%CHK_TIME
 ! extract model runoff
 CASE ('q_instnt')   ; XVAR_3d = AROUTE_3d%Q_INSTNT
 CASE ('q_routed')   ; XVAR_3d = AROUTE_3d%Q_ROUTED
 ! extract information on numerical solution (shared in MODULE model_numerix)
 CASE ('num_funcs')  ; XVAR_3d = NUM_FUNCS
 CASE ('numjacobian'); XVAR_3d = NUM_JACOBIAN
 CASE ('sub_accept') ; XVAR_3d = NUMSUB_ACCEPT
 CASE ('sub_reject') ; XVAR_3d = NUMSUB_REJECT
 CASE ('sub_noconv') ; XVAR_3d = NUMSUB_NOCONV
 CASE ('max_iterns') ; XVAR_3d = MAXNUM_ITERNS
END SELECT
! and, save the output
VAREXTRACT_3d = XVAR_3d
! ---------------------------------------------------------------------------------------
END FUNCTION VAREXTRACT_3d

END MODULE VAREXTRACT_MODULE
