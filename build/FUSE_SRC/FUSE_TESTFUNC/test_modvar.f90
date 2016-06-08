MODULE test_modvar
 USE nrtype
 ! define model fluxes
 TYPE FLUXES
  REAL(SP)                             :: DRAINAGE    ! drainage rate (mm day-1)
  REAL(SP)                             :: CHECKTIM    ! time (day)
 END TYPE FLUXES
 ! define model states
 TYPE STATES
  REAL(SP)                             :: WATR_1      ! total storage in layer1
 END TYPE STATES
 ! define state names
 TYPE SNAMES
  CHARACTER(LEN=6)                     :: SNAME       ! state name
 END TYPE SNAMES
 ! define data structures
 TYPE(FLUXES)                          :: FLUX_0      ! model fluxes at the start of the sub-step
 TYPE(FLUXES)                          :: FLUX_1      ! model fluxes at the end of the sub-step
 TYPE(FLUXES)                          :: M_FLUX      ! model fluxes
 TYPE(FLUXES)                          :: W_FLUX      ! weighted fluxes
 TYPE(STATES)                          :: DSDT_0      ! model derivatives at the start of the sub-step
 TYPE(STATES)                          :: DSDT_1      ! model derivatives at the end of the sub-step
 TYPE(STATES)                          :: MDS_DT      ! model derivatives
 TYPE(STATES)                          :: MSTATE      ! model states
 TYPE(STATES)                          :: TSTATE      ! temporary model states (used to compute derivatives)
 TYPE(STATES)                          :: FSTATE      ! final model states (at the start/end of a full step)
 TYPE(STATES)                          :: MS_MIN      ! minimum values for model states
 TYPE(STATES)                          :: MS_MAX      ! maximum values for model states
 TYPE(SNAMES),DIMENSION(1)             :: CSTATE      ! state names
 REAL(SP)                              :: DT_SUB      ! length of sub-step
 REAL(SP)                              :: DT_FULL     ! length of full step
END MODULE test_modvar 
