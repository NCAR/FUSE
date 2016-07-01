! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark
! Modified by Brian Henn to include snow model, 6/2013
! ---------------------------------------------------------------------------------------
MODULE multiparam
 USE nrtype
 USE model_defn,ONLY:NTDH_MAX
 ! --------------------------------------------------------------------------------------
 ! (1) PARAMETER METADATA
 ! --------------------------------------------------------------------------------------
 ! data structure to hold metadata for adjustable model parameters
 TYPE PARATT
  LOGICAL(LGT)                         :: PARFIT      ! flag to determine if parameter is fitted
  INTEGER(I4B)                         :: PARSTK      ! flag (0=deterministic, 1=stochastic)
  REAL(SP)                             :: PARDEF      ! default parameter set
  REAL(SP)                             :: PARLOW      ! lower limit of each parameter
  REAL(SP)                             :: PARUPP      ! upper limit of each parameter
  REAL(SP)                             :: FRSEED      ! fraction param space for "reasonable" bounds
  REAL(SP)                             :: PARSCL      ! typical scale of parameter
  INTEGER(I4B)                         :: PARVTN      ! method used for variable transformation
  INTEGER(I4B)                         :: PARDIS      ! parametric form of prob dist used for prior/hyper
  INTEGER(I4B)                         :: PARQTN      ! transformation applied before use of prob dist
  INTEGER(I4B)                         :: PARLAT      ! number of latent variables (0=onePerStep, -1=from data)
  INTEGER(I4B)                         :: PARMTH      ! imeth for all variables ???what is this???
  INTEGER(I4B)                         :: NPRIOR      ! number of prior/hyper-parameters
  CHARACTER(LEN=256)                   :: P_NAME      ! parameter name
  CHARACTER(LEN=256)                   :: CHILD1      ! name of 1st parameter child
  CHARACTER(LEN=256)                   :: CHILD2      ! name of 2nd parameter child
 END TYPE PARATT
 ! data structure to hold metadata for each parameter
 TYPE PARINFO
  ! rainfall error parameters (adjustable)
  TYPE(PARATT)                         :: RFERR_ADD   !  additive rainfall error (mm day-1)
  TYPE(PARATT)                         :: RFERR_MLT   ! multiplicative rainfall error (-)
  TYPE(PARATT)                         :: RFH1_MEAN   ! hyper parameter1: mean rainfall multiplier (-)
  TYPE(PARATT)                         :: RFH2_SDEV   ! hyper parameter2: sdev rainfall multiplier (-)
  TYPE(PARATT)                         :: RH1P_MEAN   ! prior param1 of hyper param1: prior mean of hypermean
  TYPE(PARATT)                         :: RH1P_SDEV   ! prior param2 of hyper param1: prior sdev of hypermean
  TYPE(PARATT)                         :: RH2P_MEAN   ! prior param1 of hyper param2: lower bound of hypersdev
  TYPE(PARATT)                         :: RH2P_SDEV   ! prior param2 of hyper param2: upper bound of hypersdev
  ! bucket sizes (adjustable)
  TYPE(PARATT)                         :: MAXWATR_1   ! maximum total storage in layer1 (mm)
  TYPE(PARATT)                         :: MAXWATR_2   ! maximum total storage in layer2 (mm)
  TYPE(PARATT)                         :: FRACTEN     ! frac total storage as tension storage (-)
  TYPE(PARATT)                         :: FRCHZNE     ! PRMS: frac tension storage in recharge zone (-)
  TYPE(PARATT)                         :: FPRIMQB     ! SAC: fraction of baseflow in primary resvr (-)
  ! evaporation (adjustable)
  TYPE(PARATT)                         :: RTFRAC1     ! fraction of roots in the upper layer (-)
  ! percolation (adjustable) 
  TYPE(PARATT)                         :: PERCRTE     ! percolation rate (mm day-1)
  TYPE(PARATT)                         :: PERCEXP     ! percolation exponent (-)
  TYPE(PARATT)                         :: SACPMLT     ! multiplier in the SAC model for dry lower layer (-)
  TYPE(PARATT)                         :: SACPEXP     ! exponent in the SAC model for dry lower layer (-)
  TYPE(PARATT)                         :: PERCFRAC    ! fraction of percolation to tension storage (-)
  TYPE(PARATT)                         :: FRACLOWZ    ! fraction of soil excess to lower zone (-)
  ! interflow (adjustable)
  TYPE(PARATT)                         :: IFLWRTE     ! interflow rate (mm day-1)
  ! baseflow (adjustable)
  TYPE(PARATT)                         :: BASERTE     ! baseflow rate (mm day-1)
  TYPE(PARATT)                         :: QB_POWR     ! baseflow exponent (-)
  TYPE(PARATT)                         :: QB_PRMS     ! baseflow depletion rate (day-1)
  TYPE(PARATT)                         :: QBRATE_2A   ! baseflow depletion rate for primary resvr (day-1)
  TYPE(PARATT)                         :: QBRATE_2B   ! baseflow depletion rate for secondary resvr (day-1)
  ! surface runoff (adjustable)
  TYPE(PARATT)                         :: SAREAMAX    ! maximum saturated area
  TYPE(PARATT)                         :: AXV_BEXP    ! ARNO/VIC "b" exponent
  TYPE(PARATT)                         :: LOGLAMB     ! mean value of the log-transformed topographic index (m)
  TYPE(PARATT)                         :: TISHAPE     ! shape parameter for the topo index Gamma distribution (-)
  ! time delay in runoff
  TYPE(PARATT)                         :: TIMEDELAY   ! time delay in runoff (days)
  ! snow model (adjustable)
  TYPE(PARATT)                         :: MBASE       ! base melt temperature (deg. C)
  TYPE(PARATT)                         :: MFMAX       ! maximum melt factor (mm melt deg C.-1 6hrs-1)
  TYPE(PARATT)                         :: MFMIN       ! minimum melt factor (mm melt deg C.-1 6hrs-1)
  TYPE(PARATT)                         :: PXTEMP      ! rain-snow partition temperature (deg. C)
  TYPE(PARATT)                         :: OPG         ! precipitation gradient (-) 
  TYPE(PARATT)                         :: LAPSE       ! temperature gradient (deg. C)
 ENDTYPE PARINFO
 ! --------------------------------------------------------------------------------------
 ! (2) ADJUSTABLE PARAMETERS
 ! --------------------------------------------------------------------------------------
 TYPE PARADJ
  ! rainfall error parameters (adjustable)
  REAL(SP)                             :: RFERR_ADD   ! additive rainfall error (mm day-1)
  REAL(SP)                             :: RFERR_MLT   ! multiplicative rainfall error (-)
  REAL(SP)                             :: RFH1_MEAN   ! hyper parameter1: mean rainfall multiplier (-)
  REAL(SP)                             :: RFH2_SDEV   ! hyper parameter2: sdev rainfall multiplier (-)
  REAL(SP)                             :: RH1P_MEAN   ! prior param1 of hyper param1: prior mean of hypermean
  REAL(SP)                             :: RH1P_SDEV   ! prior param2 of hyper param1: prior sdev of hypermean
  REAL(SP)                             :: RH2P_MEAN   ! prior param1 of hyper param2: lower bound of hypersdev
  REAL(SP)                             :: RH2P_SDEV   ! prior param2 of hyper param2: upper bound of hypersdev
  ! bucket sizes (adjustable)
  REAL(SP)                             :: MAXWATR_1   ! maximum total storage in layer1 (mm)
  REAL(SP)                             :: MAXWATR_2   ! maximum total storage in layer2 (mm)
  REAL(SP)                             :: FRACTEN     ! frac total storage as tension storage (-)
  REAL(SP)                             :: FRCHZNE     ! PRMS: frac tension storage in recharge zone (-)
  REAL(SP)                             :: FPRIMQB     ! SAC: fraction of baseflow in primary resvr (-)
  ! evaporation (adjustable)
  REAL(SP)                             :: RTFRAC1     ! fraction of roots in the upper layer (-)
  ! percolation (adjustable)
  REAL(SP)                             :: PERCRTE     ! percolation rate (mm day-1)
  REAL(SP)                             :: PERCEXP     ! percolation exponent (-)
  REAL(SP)                             :: SACPMLT     ! multiplier in the SAC model for dry lower layer (-)
  REAL(SP)                             :: SACPEXP     ! exponent in the SAC model for dry lower layer (-)
  REAL(SP)                             :: PERCFRAC    ! fraction of percolation to tension storage (-)
  REAL(SP)                             :: FRACLOWZ    ! fraction of soil excess to lower zone (-)
  ! interflow (adjustable)
  REAL(SP)                             :: IFLWRTE     ! interflow rate (mm day-1)
  ! baseflow (adjustable)
  REAL(SP)                             :: BASERTE     ! baseflow rate (mm day-1)
  REAL(SP)                             :: QB_POWR     ! baseflow exponent (-)
  REAL(SP)                             :: QB_PRMS     ! baseflow depletion rate (day-1)
  REAL(SP)                             :: QBRATE_2A   ! baseflow depletion rate for primary resvr (day-1)
  REAL(SP)                             :: QBRATE_2B   ! baseflow depletion rate for secondary resvr (day-1)
  ! surface runoff (adjustable)
  REAL(SP)                             :: SAREAMAX    ! maximum saturated area
  REAL(SP)                             :: AXV_BEXP    ! ARNO/VIC "b" exponent
  REAL(SP)                             :: LOGLAMB     ! mean value of the log-transformed topographic index (m)
  REAL(SP)                             :: TISHAPE     ! shape parameter for the topo index Gamma distribution (-)
  ! time delay in runoff
  REAL(SP)                             :: TIMEDELAY   ! time delay in runoff (days)
  ! snow model
  REAL(SP)                             :: MBASE       ! base melt temperature (deg. C)
  REAL(SP)                             :: MFMAX       ! maximum melt factor (mm melt deg C.-1 6hrs-1)
  REAL(SP)                             :: MFMIN       ! minimum melt factor (mm melt deg C.-1 6hrs-1)
  REAL(SP)                             :: PXTEMP      ! rain-snow partition temperature (deg. C)
  REAL(SP)                             :: OPG         ! precipitation gradient (-) 
  REAL(SP)                             :: LAPSE       ! temperature gradient (deg. C)
 END TYPE PARADJ
 ! --------------------------------------------------------------------------------------
 ! (3) DERIVED PARAMETERS
 ! --------------------------------------------------------------------------------------
 TYPE PARDVD
  ! bucket sizes (derived)
  REAL(SP)                             :: MAXTENS_1   ! maximum tension storage in layer1 (mm)
  REAL(SP)                             :: MAXTENS_2   ! maximum tension storage in layer2 (mm)
  REAL(SP)                             :: MAXFREE_1   ! maximum free storage in layer 1 (mm)
  REAL(SP)                             :: MAXFREE_2   ! maximum free storage in layer2 (mm)
  REAL(SP)                             :: MAXTENS_1A  ! maximum storage in the recharge zone (mm)
  REAL(SP)                             :: MAXTENS_1B  ! maximum storage in the lower zone (mm)
  REAL(SP)                             :: MAXFREE_2A  ! maximum storage in the primary resvr (mm)
  REAL(SP)                             :: MAXFREE_2B  ! maximum storage in the secondary resvr (mm)
  ! evaporation
  REAL(SP)                             :: RTFRAC2     ! fraction of roots in the lower layer (-)
  ! percolation/baseflow
  REAL(SP)                             :: QBSAT       ! baseflow at saturation
  ! surface runoff
  REAL(SP)                             :: POWLAMB     ! mean value of the power-transformed topographic index (m**(1/n))
  REAL(SP)                             :: MAXPOW      ! max value of the power-transformed topographic index (m**(1/n))
  ! routing
  REAL(SP), DIMENSION(NTDH_MAX)        :: FRAC_FUTURE ! fraction of runoff in future time steps
  INTEGER(I4B)                         :: NTDH_NEED   ! number of time-steps with non-zero routing contribution
 END TYPE PARDVD
 ! --------------------------------------------------------------------------------------
 ! (4) LIST OF PARAMETERS FOR A GIVEN MODEL
 ! --------------------------------------------------------------------------------------
 TYPE PAR_ID
  CHARACTER(LEN=9)                     :: PARNAME     ! list of parameter names
 ENDTYPE PAR_ID
 ! --------------------------------------------------------------------------------------
 ! (5) FINAL DATA STRUCTURES
 ! --------------------------------------------------------------------------------------
 INTEGER(I4B), PARAMETER               :: MAXPAR=50   ! maximum number of parameters for a single model
 TYPE(PARADJ), DIMENSION(:), POINTER   :: APARAM=>null()  ! all model parameter sets; DK/2008/10/21: explicit null
 TYPE(PARADJ)                          :: MPARAM      ! single model parameter set
 TYPE(PARDVD)                          :: DPARAM      ! derived model parameters
 TYPE(PARINFO)                         :: PARMETA     ! parameter metadata (all parameters)
 TYPE(PAR_ID), DIMENSION(MAXPAR)       :: LPARAM      ! list of model parameter names (need to modify to 16 for SCE)
 INTEGER(I4B)                          :: NUMPAR      ! number of model parameters for current model
 INTEGER(I4B)                          :: SOBOL_INDX  ! code to re-assemble Sobol parameters
 ! --------------------------------------------------------------------------------------
END MODULE multiparam
