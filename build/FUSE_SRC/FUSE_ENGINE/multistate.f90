MODULE multistate
 USE nrtype
 ! --------------------------------------------------------------------------------------
 ! model state structure
 ! --------------------------------------------------------------------------------------
 TYPE STATEV
  ! snow layer
  REAL(SP)                             :: SWE_TOT    ! total storage as snow (mm)
  ! upper layer
  REAL(SP)                             :: WATR_1     ! total storage in layer1 (mm)
  REAL(SP)                             :: TENS_1     ! tension storage in layer1 (mm)
  REAL(SP)                             :: FREE_1     ! free storage in layer 1 (mm)
  REAL(SP)                             :: TENS_1A    ! storage in the recharge zone (mm)
  REAL(SP)                             :: TENS_1B    ! storage in the lower zone (mm)
  ! lower layer
  REAL(SP)                             :: WATR_2     ! total storage in layer2 (mm)
  REAL(SP)                             :: TENS_2     ! tension storage in layer2 (mm)
  REAL(SP)                             :: FREE_2     ! free storage in layer2 (mm)
  REAL(SP)                             :: FREE_2A    ! storage in the primary resvr (mm)
  REAL(SP)                             :: FREE_2B    ! storage in the secondary resvr (mm)
 END TYPE STATEV
 ! --------------------------------------------------------------------------------------
 ! model time structure
 ! --------------------------------------------------------------------------------------
 TYPE M_TIME
  REAL(SP)                             :: STEP       ! (time interval to advance model states)
 END TYPE M_TIME
 ! --------------------------------------------------------------------------------------
 ! variable definitions
 ! --------------------------------------------------------------------------------------
 type(statev),dimension(:,:),pointer   :: gState     ! (grid of model states)
 type(statev),dimension(:,:,:),pointer :: gState_3d  ! (grid of model states with a time dimension)
 TYPE(STATEV)                          :: ASTATE     ! (model states at the start of full timestep)
 TYPE(STATEV)                          :: FSTATE     ! (model states at start of sub-timestep)
 TYPE(STATEV)                          :: MSTATE     ! (model states at start/middle of sub-timestep)
 TYPE(STATEV)                          :: TSTATE     ! (temporary copy of model states)
 TYPE(STATEV)                          :: BSTATE     ! (temporary copy of model states)
 TYPE(STATEV)                          :: ESTATE     ! (temporary copy of model states)
 TYPE(STATEV)                          :: DSTATE     ! (default model states)
 TYPE(STATEV)                          :: DYDT_0     ! (derivative of model states at start of sub-step)
 TYPE(STATEV)                          :: DYDT_1     ! (derivative of model states at end of sub-step)
 TYPE(STATEV)                          :: DY_DT      ! (derivative of model states)
 TYPE(STATEV)                          :: DYDT_OLD   ! (derivative of model states for final solution)
 TYPE(M_TIME)                          :: HSTATE     ! (time interval to advance model states)
 ! --------------------------------------------------------------------------------------

 ! NetCDF
 integer(i4b)                          :: ncid_out=-1              ! NetCDF output file ID

 ! initial store fraction (initialization)
 real(sp),parameter::fracState0=0.25_sp

END MODULE multistate
