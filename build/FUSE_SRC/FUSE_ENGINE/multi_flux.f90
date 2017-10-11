MODULE multi_flux
 USE nrtype
 TYPE FLUXES
  REAL(SP)                             :: EFF_PPT     ! effective precipitation (mm day-1)
  REAL(SP)                             :: SATAREA     ! saturated area (-)
  REAL(SP)                             :: QSURF       ! surface runoff (mm day-1)
  REAL(SP)                             :: EVAP_1A     ! evaporation from soil excess zone (mm day-1)
  REAL(SP)                             :: EVAP_1B     ! evaporation from soil recharge zone (mm day-1)
  REAL(SP)                             :: EVAP_1      ! evaporation from upper soil layer (mm day-1)
  REAL(SP)                             :: EVAP_2      ! evaporation from lower soil layer (mm day-1)
  REAL(SP)                             :: RCHR2EXCS   ! flow from recharge to excess (mm day-1)
  REAL(SP)                             :: TENS2FREE_1 ! flow from tension storage to free storage (mm day-1)
  REAL(SP)                             :: TENS2FREE_2 ! flow from tension storage to free storage (mm day-1)
  REAL(SP)                             :: QINTF_1     ! interflow from free water (mm day-1)
  REAL(SP)                             :: QPERC_12    ! percolation from upper to lower soil layers (mm day-1)
  REAL(SP)                             :: QBASE_2     ! baseflow (mm day-1)
  REAL(SP)                             :: QBASE_2A    ! baseflow from primary linear resvr (mm day-1)
  REAL(SP)                             :: QBASE_2B    ! baseflow from secondary linear resvr (mm day-1)
  REAL(SP)                             :: OFLOW_1     ! bucket overflow (mm day-1)
  REAL(SP)                             :: OFLOW_2     ! bucket overflow (mm day-1)
  REAL(SP)                             :: OFLOW_2A    ! bucket overflow (mm day-1)
  REAL(SP)                             :: OFLOW_2B    ! bucket overflow (mm day-1)
  REAL(SP)                             :: ERR_WATR_1  ! excessive extrapolation: total storage in layer1 (mm day-1)
  REAL(SP)                             :: ERR_TENS_1  ! excessive extrapolation: tension storage in layer1 (mm day-1)
  REAL(SP)                             :: ERR_FREE_1  ! excessive extrapolation: free storage in layer 1 (mm day-1)
  REAL(SP)                             :: ERR_TENS_1A ! excessive extrapolation: storage in the recharge zone (mm day-1)
  REAL(SP)                             :: ERR_TENS_1B ! excessive extrapolation: storage in the lower zone (mm day-1)
  REAL(SP)                             :: ERR_WATR_2  ! excessive extrapolation: total storage in layer2 (mm day-1)
  REAL(SP)                             :: ERR_TENS_2  ! excessive extrapolation: tension storage in layer2 (mm day-1)
  REAL(SP)                             :: ERR_FREE_2  ! excessive extrapolation: free storage in layer2 (mm day-1)
  REAL(SP)                             :: ERR_FREE_2A ! excessive extrapolation: storage in the primary resvr (mm day-1)
  REAL(SP)                             :: ERR_FREE_2B ! excessive extrapolation: storage in the secondary resvr (mm day-1)
  REAL(SP)                             :: CHK_TIME    ! time elapsed during time step (days)
 ENDTYPE FLUXES
 TYPE(FLUXES)                          :: M_FLUX      ! model fluxes
 TYPE(FLUXES)                          :: FLUX_0      ! model fluxes at start of step
 TYPE(FLUXES)                          :: FLUX_1      ! model fluxes at end of step
 TYPE(FLUXES), DIMENSION(:), POINTER   :: FDFLUX=>NULL() ! finite difference fluxes
 TYPE(FLUXES)                          :: W_FLUX      ! weighted sum of model fluxes over a time step
 TYPE(FLUXES), dimension(:,:,:), allocatable  :: W_FLUX_3d   ! weighted sum of model fluxes over a time step for several time steps
 REAL(SP)                              :: CURRENT_DT  ! current time step (days)
END MODULE multi_flux
