MODULE multiroute
 USE nrtype
 USE model_defn,ONLY:NTDH_MAX
 TYPE RUNOFF
  REAL(SP)                             :: Q_INSTNT   ! instantaneous runoff
  REAL(SP)                             :: Q_ROUTED   ! routed runoff
  REAL(SP)                             :: Q_ACCURATE ! "accurate" runoff estimate (mm day-1)
 END TYPE RUNOFF
 REAL(SP), DIMENSION(NTDH_MAX)         :: FUTURE     ! runoff placed in future time steps
 TYPE(RUNOFF), DIMENSION(:), POINTER   :: AROUTE     ! runoff for all time steps
 TYPE(RUNOFF),dimension(:,:,:), allocatable   :: AROUTE_3d     ! runoff for all time steps on a grid
 TYPE(RUNOFF)                          :: MROUTE     ! runoff for one time step
END MODULE multiroute
