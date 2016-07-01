MODULE GETPAR_STR_MODULE
IMPLICIT NONE
CONTAINS
SUBROUTINE GETPAR_STR(PARNAME,METADAT)
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
CHARACTER(*), INTENT(IN)               :: PARNAME     ! parameter name
TYPE(PARATT), INTENT(OUT)              :: METADAT     ! parameter metadata
! ---------------------------------------------------------------------------------------
! model parameters
SELECTCASE(TRIM(PARNAME))
CASE('RFERR_ADD');  METADAT = PARMETA%RFERR_ADD
CASE('RFERR_MLT');  METADAT = PARMETA%RFERR_MLT
CASE('RFH1_MEAN');  METADAT = PARMETA%RFH1_MEAN
CASE('RFH2_SDEV');  METADAT = PARMETA%RFH2_SDEV
CASE('RH1P_MEAN');  METADAT = PARMETA%RH1P_MEAN
CASE('RH1P_SDEV');  METADAT = PARMETA%RH1P_SDEV
CASE('RH2P_MEAN');  METADAT = PARMETA%RH2P_MEAN
CASE('RH2P_SDEV');  METADAT = PARMETA%RH2P_SDEV
CASE('MAXWATR_1');  METADAT = PARMETA%MAXWATR_1
CASE('MAXWATR_2');  METADAT = PARMETA%MAXWATR_2
CASE('FRACTEN');    METADAT = PARMETA%FRACTEN  
CASE('FRCHZNE');    METADAT = PARMETA%FRCHZNE  
CASE('FPRIMQB');    METADAT = PARMETA%FPRIMQB  
CASE('RTFRAC1');    METADAT = PARMETA%RTFRAC1  
CASE('PERCRTE');    METADAT = PARMETA%PERCRTE  
CASE('PERCEXP');    METADAT = PARMETA%PERCEXP  
CASE('SACPMLT');    METADAT = PARMETA%SACPMLT  
CASE('SACPEXP');    METADAT = PARMETA%SACPEXP  
CASE('PERCFRAC');   METADAT = PARMETA%PERCFRAC 
CASE('FRACLOWZ');   METADAT = PARMETA%FRACLOWZ 
CASE('IFLWRTE');    METADAT = PARMETA%IFLWRTE  
CASE('BASERTE');    METADAT = PARMETA%BASERTE  
CASE('QB_POWR');    METADAT = PARMETA%QB_POWR  
CASE('QB_PRMS');    METADAT = PARMETA%QB_PRMS  
CASE('QBRATE_2A');  METADAT = PARMETA%QBRATE_2A
CASE('QBRATE_2B');  METADAT = PARMETA%QBRATE_2B
CASE('SAREAMAX');   METADAT = PARMETA%SAREAMAX 
CASE('AXV_BEXP');   METADAT = PARMETA%AXV_BEXP 
CASE('LOGLAMB');    METADAT = PARMETA%LOGLAMB 
CASE('TISHAPE');    METADAT = PARMETA%TISHAPE  
CASE('TIMEDELAY');  METADAT = PARMETA%TIMEDELAY
CASE('MBASE');      METADAT = PARMETA%MBASE
CASE('MFMAX');      METADAT = PARMETA%MFMAX
CASE('MFMIN');      METADAT = PARMETA%MFMIN
CASE('PXTEMP');     METADAT = PARMETA%PXTEMP
CASE('OPG');        METADAT = PARMETA%OPG
CASE('LAPSE');      METADAT = PARMETA%LAPSE
CASE DEFAULT
 print *, 'parameter name (', TRIM(PARNAME), ') does not exist '
 IF (TRIM(PARNAME).EQ.'NO_CHILD1' .OR. TRIM(PARNAME).EQ.'NO_CHILD2') &
  print *, ' * check the number of prior/hyper parameters specified '
 STOP
ENDSELECT
! ---------------------------------------------------------------------------------------
END SUBROUTINE GETPAR_STR
END MODULE GETPAR_STR_MODULE
