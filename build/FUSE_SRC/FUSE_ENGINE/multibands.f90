! Created by Brian Henn to allow multi-band snow modeling, 6/2013
! Based on module MULTIFORCE by Martyn Clark
MODULE multibands
 USE nrtype
 TYPE BANDS
  INTEGER(I4B)                         :: NUM             ! band number (-)
  REAL(SP)                             :: Z_MID           ! band mid-point elevation (m)
  REAL(SP)                             :: AF              ! fraction of basin area in band (-)
  REAL(SP)                             :: SWE             ! band snowpack water equivalent (mm)
  REAL(SP)                             :: SNOWACCMLTN     ! new snow accumulation in band (mm day-1)
  REAL(SP)                             :: SNOWMELT        ! snowmelt in band (mm day-1)
  REAL(SP)                             :: DSWE_DT         ! rate of change of band SWE (mm day-1) 
 ENDTYPE BANDS
 ! --------------------------------------------------------------------------------------
 TYPE(BANDS),DIMENSION(:),ALLOCATABLE  :: MBANDS          ! basin band information
 INTEGER(I4B)                          :: N_BANDS=0       ! number of bands, initialize to zero
 REAL(SP)                              :: Z_FORCING       ! elevation of forcing data (m)
 ! --------------------------------------------------------------------------------------
END MODULE multibands
