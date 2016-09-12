SUBROUTINE ASSIGN_PAR()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! Modified by Brian Henn to include snow model, 6/2013
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Gets a list of model parameters used for the unique model in the structure SMODL
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multi_flux -- list of model parameters is stored in MODULE multiparam
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition structure
USE model_defnames
USE multiparam, ONLY : lparam, paratt, numpar         ! model parameter structures
USE getpar_str_module                                 ! access to SUBROUTINE get_par_str
IMPLICIT NONE
INTEGER(I4B)                           :: MPAR        ! counter for number of parameters
TYPE(PARATT)                           :: PARAM_LEV1  ! parameter metadata (level 1)
TYPE(PARATT)                           :: PARAM_LEV2  ! parameter metadata (level 2)
! ---------------------------------------------------------------------------------------
MPAR = 0  ! initialize the number of model parameters
LPARAM(:)%PARNAME = 'PAR_NOUSE'
! ---------------------------------------------------------------------------------------
! (1) RAINFALL ERRORS
! ---------------------------------------------------------------------------------------

SELECT CASE(SMODL%iRFERR)
 CASE(iopt_additive_e) ! additive rainfall error
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'RFERR_ADD'  ! additive rainfall error (mm day-1)
 CASE(iopt_multiplc_e) ! multiplicative rainfall error
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'RFERR_MLT'  ! multiplicative rainfall error (-)
  ! check if RFERR_MLT has any prior/hyper-parameters, and, if so, save them
  CALL GETPAR_STR('RFERR_MLT',PARAM_LEV1)
  IF (PARAM_LEV1%NPRIOR.GT.0) THEN
   ! process 1st child
   MPAR=MPAR+1; LPARAM(MPAR)%PARNAME=PARAM_LEV1%CHILD1(1:9)  ! save 1st child
   CALL GETPAR_STR(PARAM_LEV1%CHILD1,PARAM_LEV2)             ! get metadata for 1st child
   IF (PARAM_LEV2%NPRIOR.GT.0) THEN  ! check if 1st child has prior/hyper-param
    MPAR=MPAR+1; LPARAM(MPAR)%PARNAME=PARAM_LEV2%CHILD1(1:9) ! save 1st grandchild (from 1st child)
    MPAR=MPAR+1; LPARAM(MPAR)%PARNAME=PARAM_LEV2%CHILD2(1:9) ! save 2nd grandchild (from 1st child)
   ENDIF
   ! process 2nd child
   MPAR=MPAR+1; LPARAM(MPAR)%PARNAME=PARAM_LEV1%CHILD2(1:9)  ! save 2nd child
   CALL GETPAR_STR(PARAM_LEV1%CHILD2,PARAM_LEV2)             ! get metadata for 1st child
   IF (PARAM_LEV2%NPRIOR.GT.0) THEN  ! check if 1st child has prior/hyper-param
    MPAR=MPAR+1; LPARAM(MPAR)%PARNAME=PARAM_LEV2%CHILD1(1:9) ! save 1st grandchild (from 2nd child)
    MPAR=MPAR+1; LPARAM(MPAR)%PARNAME=PARAM_LEV2%CHILD2(1:9) ! save 2nd grandchild (from 2nd child)
   ENDIF
  ENDIF
 CASE DEFAULT
  print *, "SMODL%RFERR must be 'additive_e' or 'multiplc_e'"
  STOP
END SELECT  ! (different upper-layer architecture)
! ---------------------------------------------------------------------------------------
! (2) UPPER-LAYER ARCHITECTURE
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iARCH1)
 CASE(iopt_tension2_1) ! tension storage sub-divided into recharge and excess 
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'FRCHZNE  '  ! PRMS: frac tension storage in recharge zone (-)
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'FRACTEN  '  ! frac total storage as tension storage (-)
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'MAXWATR_1'  ! maximum total storage in layer1 (mm)
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'FRACLOWZ '  ! fraction of soil excess to lower zone (-)
 CASE(iopt_tension1_1,iopt_onestate_1) ! (need to define tension and free storage -- even if one state)
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'FRACTEN  '  ! frac total storage as tension storage (-)
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'MAXWATR_1'  ! maximum total storage in layer1 (mm)
 CASE DEFAULT
  print *, "SMODL%iARCH1 must be iopt_tension2_1, iopt_tension1_1, or iopt_onestate_1"
  STOP
END SELECT  ! (different upper-layer architechure)
! ---------------------------------------------------------------------------------------
! (3) LOWER-LAYER ARCHITECTURE / BASEFLOW
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iARCH2)
 CASE(iopt_tens2pll_2) ! tension reservoir plus two parallel tanks
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'PERCFRAC '  ! fraction of percolation to tension storage (-)
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'FPRIMQB  '  ! SAC: fraction of baseflow in primary resvr (-)
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'MAXWATR_2'  ! maximum total storage in layer2 (mm)
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'QBRATE_2A ' ! baseflow depletion rate for primary resvr (day-1)
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'QBRATE_2B ' ! baseflow depletion rate for secondary resvr (day-1)
 CASE(iopt_unlimfrc_2) ! baseflow resvr of unlimited size (0-HUGE), frac rate
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'MAXWATR_2'  ! maximum total storage in layer2 (mm)
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'QB_PRMS  '  ! baseflow depletion rate (day-1)
 CASE(iopt_topmdexp_2,iopt_unlimpow_2) ! topmodel options
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'MAXWATR_2'  ! maximum total storage in layer2 (mm)
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'BASERTE  '  ! baseflow rate (mm day-1)
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'LOGLAMB  '  ! mean value of the log-transformed topographic index (m)
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'TISHAPE  '  ! shape parameter for the topo index Gamma distribution (-)
  ! (add extra paramater for the power-law transmissivity profile)
  IF (SMODL%iARCH2.EQ.iopt_unlimpow_2) THEN ! (power-law transmissivity profile)
   MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'QB_POWR  '  ! baseflow exponent (-)
  ENDIF
 CASE(iopt_fixedsiz_2)  ! power-law relation (no parameters needed for the topo index distribution)
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'MAXWATR_2'  ! maximum total storage in layer2 (mm)
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'BASERTE  '  ! baseflow rate (mm day-1)
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'QB_POWR  '  ! baseflow exponent (-)
 CASE DEFAULT
  print *, "SMODL%iARCH2 must be iopt_tens2pll_2, iopt_unlimfrc_2, iopt_unlimpow_2"
  print *, "  'topmdexp_2', or 'fixedsiz_2'"
  STOP
END SELECT  ! different lower-layer architecture / baseflow parameterizations)
! ---------------------------------------------------------------------------------------
! (4) EVAPORATION
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iESOIL)
 CASE(iopt_sequential)
  ! (no additional parameters for the sequential scheme)
 CASE(iopt_rootweight)
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'RTFRAC1  '  ! fraction of roots in the upper layer (-)
 CASE DEFAULT
  print *, "SMODL%iESOIL must be either iopt_sequential or iopt_rootweight'"
END SELECT  ! (different evaporation schemes)
! ---------------------------------------------------------------------------------------
! (5) PERCOLATION
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iQPERC)
 CASE(iopt_perc_f2sat,iopt_perc_w2sat) ! standard equation k(theta)**c
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'PERCRTE  '  ! percolation rate (mm day-1)
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'PERCEXP  '  ! percolation exponent (-)
 CASE(iopt_perc_lower) ! perc defined by moisture content in lower layer (SAC)
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'SACPMLT  '  ! multiplier in the SAC model for dry lower layer (-)
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'SACPEXP  '  ! exponent in the SAC model for dry lower layer (-)
 CASE DEFAULT       ! check for errors
  print *, "SMODL%iQPERC must be iopt_perc_f2sat, iopt_perc_w2sat, or iopt_perc_lower"
  STOP
END SELECT  ! (different percolation options)
! ---------------------------------------------------------------------------------------
! (6) INTERFLOW
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iQINTF)
 CASE(iopt_intflwsome) ! interflow
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'IFLWRTE  '  ! interflow rate (mm day-1)
 CASE(iopt_intflwnone) ! no interflow
  ! (no additional parameters for the case of no interflow)
 CASE DEFAULT       ! check for errors
  print *, "SMODL%iQINTF must be either iopt_intflwsome' or iopt_intflwnone'"
  STOP
END SELECT  ! (different interflow options)
! ---------------------------------------------------------------------------------------
! (7) SURFACE RUNOFF
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iQSURF)
 CASE(iopt_arno_x_vic) ! ARNO/Xzang/VIC parameterization (upper zone control)
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'AXV_BEXP '  ! ARNO/VIC "b" exponent
 CASE(iopt_prms_varnt) ! PRMS variant (fraction of upper tension storage)
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'SAREAMAX '  ! maximum saturated area
 CASE(iopt_tmdl_param) ! TOPMODEL parameterization
  ! need the topographic index if we don't have it for baseflow
  IF (SMODL%iARCH2.EQ.iopt_tens2pll_2 .OR. SMODL%iARCH2.EQ.iopt_unlimfrc_2 .OR. &
      SMODL%iARCH2.EQ.iopt_fixedsiz_2) THEN
   MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'LOGLAMB  '  ! mean value of the log-transformed topographic index (m)
   MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'TISHAPE  '  ! shape parameter for the topo index Gamma distribution (-)
  ENDIF
  ! need the topmodel power if we don't have it for baseflow
  IF (SMODL%iARCH2.EQ.iopt_tens2pll_2 .OR. SMODL%iARCH2.EQ.iopt_unlimfrc_2 .OR. &
      SMODL%iARCH2.EQ.iopt_topmdexp_2) THEN
   MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'QB_POWR  '  ! baseflow exponent (-), used to modify the topographic index
  ENDIF
 CASE DEFAULT
  print *, "SMODL%iQSURF must be iopt_arno_x_vic, iopt_prms_varnt, or iopt_tmdl_param"
  STOP
END SELECT  ! (different surface runoff options)
! ---------------------------------------------------------------------------------------
! (8) TIME DELAY IN RUNOFF
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iQ_TDH)
 CASE(iopt_rout_gamma) ! use a Gamma distribution with shape parameter = 2.5
  MPAR=MPAR+1; LPARAM(MPAR)%PARNAME = 'TIMEDELAY'  ! time delay in runoff
 CASE(iopt_no_routing) ! no routing
  ! (no additional parameters when there is no time delay in runoff)
 CASE DEFAULT       ! check for errors
  print *, "SMODL%iQ_TDH must be either iopt_rout_gamma or iopt_no_routing"
  STOP
END SELECT
! ---------------------------------------------------------------------------------------
! (9) SNOW MODEL 
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iSNOWM)
 CASE(iopt_temp_index) ! temperature index snow model
  MPAR = MPAR + 1; LPARAM(MPAR)%PARNAME = 'MBASE    ' ! snow base melting temperature
  MPAR = MPAR + 1; LPARAM(MPAR)%PARNAME = 'MFMAX    ' ! snow maximum melt factor
  MPAR = MPAR + 1; LPARAM(MPAR)%PARNAME = 'MFMIN    ' ! snow minimum melt factor
  MPAR = MPAR + 1; LPARAM(MPAR)%PARNAME = 'PXTEMP   ' ! rain snow partition temperature
  MPAR = MPAR + 1; LPARAM(MPAR)%PARNAME = 'OPG      ' ! precipitation gradient
  MPAR = MPAR + 1; LPARAM(MPAR)%PARNAME = 'LAPSE    ' ! temperature gradient
 CASE(iopt_no_snowmod) ! if no snow model, no additional parameters
 CASE DEFAULT
  print *, "SMODL%SNOWM must be either 'temp_index' or 'no_snowmod'"
  STOP
END SELECT
! ---------------------------------------------------------------------------------------
NUMPAR = MPAR  ! save the number of model parameters used in a given model SMODL
! ---------------------------------------------------------------------------------------
!DO MPAR=1,NUMPAR; WRITE(*,'(A11,1X)') LPARAM(MPAR)%PARNAME; END DO
! ---------------------------------------------------------------------------------------
END SUBROUTINE ASSIGN_PAR
