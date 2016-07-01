MODULE PUTPAR_STR_MODULE
IMPLICIT NONE
CONTAINS
SUBROUTINE PUTPAR_STR(METADAT,PARNAME)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! Modified by Brian Henn to include snow model, 6/2013
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Inserts parameter metadata into data structures
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE multiparam, ONLY : PARATT, PARMETA                ! derived type for parameter metadata
IMPLICIT NONE
! input
TYPE(PARATT), INTENT(IN)               :: METADAT     ! parameter metadata
CHARACTER(*), INTENT(IN)               :: PARNAME     ! parameter name
! ---------------------------------------------------------------------------------------
! model parameters
SELECTCASE(TRIM(PARNAME))
CASE('RFERR_ADD');  PARMETA%RFERR_ADD = METADAT 
CASE('RFERR_MLT');  PARMETA%RFERR_MLT = METADAT
CASE('RFH1_MEAN');  PARMETA%RFH1_MEAN = METADAT
CASE('RFH2_SDEV');  PARMETA%RFH2_SDEV = METADAT
CASE('RH1P_MEAN');  PARMETA%RH1P_MEAN = METADAT
CASE('RH1P_SDEV');  PARMETA%RH1P_SDEV = METADAT
CASE('RH2P_MEAN');  PARMETA%RH2P_MEAN = METADAT
CASE('RH2P_SDEV');  PARMETA%RH2P_SDEV = METADAT
CASE('MAXWATR_1');  PARMETA%MAXWATR_1 = METADAT
CASE('MAXWATR_2');  PARMETA%MAXWATR_2 = METADAT
CASE('FRACTEN');    PARMETA%FRACTEN   = METADAT
CASE('FRCHZNE');    PARMETA%FRCHZNE   = METADAT
CASE('FPRIMQB');    PARMETA%FPRIMQB   = METADAT
CASE('RTFRAC1');    PARMETA%RTFRAC1   = METADAT
CASE('PERCRTE');    PARMETA%PERCRTE   = METADAT
CASE('PERCEXP');    PARMETA%PERCEXP   = METADAT
CASE('SACPMLT');    PARMETA%SACPMLT   = METADAT
CASE('SACPEXP');    PARMETA%SACPEXP   = METADAT
CASE('PERCFRAC');   PARMETA%PERCFRAC  = METADAT
CASE('FRACLOWZ');   PARMETA%FRACLOWZ  = METADAT
CASE('IFLWRTE');    PARMETA%IFLWRTE   = METADAT
CASE('BASERTE');    PARMETA%BASERTE   = METADAT
CASE('QB_POWR');    PARMETA%QB_POWR   = METADAT
CASE('QB_PRMS');    PARMETA%QB_PRMS   = METADAT
CASE('QBRATE_2A');  PARMETA%QBRATE_2A = METADAT
CASE('QBRATE_2B');  PARMETA%QBRATE_2B = METADAT
CASE('SAREAMAX');   PARMETA%SAREAMAX  = METADAT
CASE('AXV_BEXP');   PARMETA%AXV_BEXP  = METADAT
CASE('LOGLAMB');    PARMETA%LOGLAMB   = METADAT
CASE('TISHAPE');    PARMETA%TISHAPE   = METADAT
CASE('TIMEDELAY');  PARMETA%TIMEDELAY = METADAT
CASE('MBASE');      PARMETA%MBASE     = METADAT
CASE('MFMAX');      PARMETA%MFMAX     = METADAT
CASE('MFMIN');      PARMETA%MFMIN     = METADAT
CASE('PXTEMP');     PARMETA%PXTEMP    = METADAT
CASE('OPG');        PARMETA%OPG       = METADAT
CASE('LAPSE');      PARMETA%LAPSE     = METADAT
CASE DEFAULT
 print *, 'parameter name (', TRIM(PARNAME), ') does not exist'; STOP
ENDSELECT
! ---------------------------------------------------------------------------------------
END SUBROUTINE PUTPAR_STR
END MODULE PUTPAR_STR_MODULE
