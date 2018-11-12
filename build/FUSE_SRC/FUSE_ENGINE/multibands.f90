! Created by Brian Henn to allow multi-band snow modeling, 6/2013
! Based on module MULTIFORCE by Martyn Clark
MODULE multibands
 USE nrtype
 TYPE BANDS ! for catchment scale modeling
  INTEGER(I4B)                         :: NUM             ! band number (-)
  REAL(SP)                             :: Z_MID           ! band mid-point elevation (m)
  REAL(SP)                             :: AF              ! fraction of basin area in band (-)
  REAL(SP)                             :: SWE             ! band snowpack water equivalent (mm)
  REAL(SP)                             :: SNOWACCMLTN     ! new snow accumulation in band (mm day-1)
  REAL(SP)                             :: SNOWMELT        ! snowmelt in band (mm day-1)
  REAL(SP)                             :: DSWE_DT         ! rate of change of band SWE (mm day-1)
 ENDTYPE BANDS

 ! for distributed modeling MBANDS is split between time-independent and time-dependent charactertistics

 TYPE BANDS_INFO ! invariant characteristics
  REAL(SP)                             :: Z_MID           ! band mid-point elevation (m)
  REAL(SP)                             :: AF              ! fraction of basin area in band (-)
 ENDTYPE BANDS_INFO

 TYPE BANDS_VAR ! time-dependent characteristics
  REAL(SP)                             :: SWE             ! band snowpack water equivalent (mm)
  REAL(SP)                             :: SNOWACCMLTN     ! new snow accumulation in band (mm day-1)
  REAL(SP)                             :: SNOWMELT        ! snowmelt in band (mm day-1)
  REAL(SP)                             :: DSWE_DT         ! rate of change of band SWE (mm day-1)
 ENDTYPE BANDS_VAR

 ! --------------------------------------------------------------------------------------
 TYPE(BANDS),DIMENSION(:),ALLOCATABLE  :: MBANDS          ! basin band information
 type(BANDS_INFO),dimension(:,:,:),ALLOCATABLE :: MBANDS_INFO_3d    ! basin band information in space
 type(BANDS_VAR),dimension(:,:,:,:),ALLOCATABLE :: MBANDS_VAR_4d    ! basin band information in space plus time

 INTEGER(I4B)                          :: N_BANDS=0       ! number of bands, initialize to zero
 REAL(SP)                              :: Z_FORCING       ! elevation of forcing data (m)
 REAL(SP),DIMENSION(:,:),ALLOCATABLE   :: Z_FORCING_grid  ! elevation of forcing data (m) for the 2D domain
 LOGICAL(LGT),DIMENSION(:,:),ALLOCATABLE   :: elev_mask   ! mask domain - TRUE means the cell must be masked, i.e. not run
 ! --------------------------------------------------------------------------------------
END MODULE multibands
