module set_all_module
USE nrtype
USE netcdf
implicit none

contains

  SUBROUTINE SET_STATE(VAL)

    ! ---------------------------------------------------------------------------------------
    ! Creator:
    ! --------
    ! Nans Addor based on Martyn Clark's INIT_STATE
    ! ---------------------------------------------------------------------------------------
    ! Purpose:
    ! --------
    ! Set model states to a given value - useful to set them to _FillValue
    ! ---------------------------------------------------------------------------------------

    USE multiparam                                        ! model parameters
    USE multistate                                        ! model states
    USE multibands                                        ! model snow bands
    IMPLICIT NONE
    REAL(SP), INTENT(IN)                  :: VAL          ! value
    INTEGER(I4B)                          :: ISNW         ! snow band index
    ! ---------------------------------------------------------------------------------------
    ! upper layer
    FSTATE%TENS_1A = VAL
    FSTATE%TENS_1B = VAL
    FSTATE%TENS_1  = VAL
    FSTATE%FREE_1  = VAL
    FSTATE%WATR_1  = VAL
    ! lower layer
    FSTATE%TENS_2  = VAL
    FSTATE%FREE_2  = VAL
    FSTATE%FREE_2A = VAL
    FSTATE%FREE_2B = VAL
    FSTATE%WATR_2  = VAL
    ! snow model

    DO ISNW=1,N_BANDS
     MBANDS(ISNW)%SWE = VAL
    END DO

    FSTATE%SWE_TOT = VAL

    ! ---------------------------------------------------------------------------------------
  END SUBROUTINE SET_STATE

  SUBROUTINE SET_FLUXES(VAL)

    ! ---------------------------------------------------------------------------------------
    ! Creator:
    ! --------
    ! Nans Addor based on Martyn Clark's INITFLUXES
    ! ---------------------------------------------------------------------------------------
    ! Purpose:
    ! --------
    ! Set model fluxes to a given value - useful to set them to _FillValue
    ! ---------------------------------------------------------------------------------------

    USE model_defn                                        ! model decision structures
    USE model_defnames                                    ! integer model definitions
    USE multi_flux                                        ! model fluxes
    USE multibands                                        ! model snow bands
    IMPLICIT NONE
    REAL(SP), INTENT(IN)                  :: VAL          ! value
    INTEGER(I4B)                          :: ISNW         ! index for looping though SWE
    ! ---------------------------------------------------------------------------------------
    M_FLUX%EFF_PPT     = VAL; W_FLUX%EFF_PPT     = VAL
    M_FLUX%SATAREA     = VAL; W_FLUX%SATAREA     = VAL
    M_FLUX%QSURF       = VAL; W_FLUX%QSURF       = VAL
    M_FLUX%EVAP_1A     = VAL; W_FLUX%EVAP_1A     = VAL
    M_FLUX%EVAP_1B     = VAL; W_FLUX%EVAP_1B     = VAL
    M_FLUX%EVAP_1      = VAL; W_FLUX%EVAP_1      = VAL
    M_FLUX%EVAP_2      = VAL; W_FLUX%EVAP_2      = VAL
    M_FLUX%RCHR2EXCS   = VAL; W_FLUX%RCHR2EXCS   = VAL
    M_FLUX%TENS2FREE_1 = VAL; W_FLUX%TENS2FREE_1 = VAL
    M_FLUX%TENS2FREE_2 = VAL; W_FLUX%TENS2FREE_2 = VAL
    M_FLUX%QINTF_1     = VAL; W_FLUX%QINTF_1     = VAL
    M_FLUX%QPERC_12    = VAL; W_FLUX%QPERC_12    = VAL
    M_FLUX%QBASE_2     = VAL; W_FLUX%QBASE_2     = VAL
    M_FLUX%QBASE_2A    = VAL; W_FLUX%QBASE_2A    = VAL
    M_FLUX%QBASE_2B    = VAL; W_FLUX%QBASE_2B    = VAL
    M_FLUX%OFLOW_1     = VAL; W_FLUX%OFLOW_1     = VAL
    M_FLUX%OFLOW_2     = VAL; W_FLUX%OFLOW_2     = VAL
    M_FLUX%OFLOW_2A    = VAL; W_FLUX%OFLOW_2A    = VAL
    M_FLUX%OFLOW_2B    = VAL; W_FLUX%OFLOW_2B    = VAL
    IF(SMODL%iSNOWM.EQ.iopt_temp_index) THEN !loop through snow model bands
     DO ISNW=1,N_BANDS
      MBANDS(ISNW)%SNOWACCMLTN  = VAL
      MBANDS(ISNW)%SNOWMELT     = VAL
     END DO
    ENDIF
    M_FLUX%ERR_WATR_1  = VAL; W_FLUX%ERR_WATR_1  = VAL
    M_FLUX%ERR_TENS_1  = VAL; W_FLUX%ERR_TENS_1  = VAL
    M_FLUX%ERR_FREE_1  = VAL; W_FLUX%ERR_FREE_1  = VAL
    M_FLUX%ERR_TENS_1A = VAL; W_FLUX%ERR_TENS_1A = VAL
    M_FLUX%ERR_TENS_1B = VAL; W_FLUX%ERR_TENS_1B = VAL
    M_FLUX%ERR_WATR_2  = VAL; W_FLUX%ERR_WATR_2  = VAL
    M_FLUX%ERR_TENS_2  = VAL; W_FLUX%ERR_TENS_2  = VAL
    M_FLUX%ERR_FREE_2  = VAL; W_FLUX%ERR_FREE_2  = VAL
    M_FLUX%ERR_FREE_2A = VAL; W_FLUX%ERR_FREE_2A = VAL
    M_FLUX%ERR_FREE_2B = VAL; W_FLUX%ERR_FREE_2B = VAL
    M_FLUX%CHK_TIME    = VAL; W_FLUX%CHK_TIME    = VAL
    ! ---------------------------------------------------------------------------------------

  END SUBROUTINE SET_FLUXES

  SUBROUTINE SET_ROUTE(VAL)

    ! ---------------------------------------------------------------------------------------
    ! Creator:
    ! --------
    ! Nans Addor based on Martyn Clark's INIT_STATE
    ! ---------------------------------------------------------------------------------------
    ! Purpose:
    ! --------
    ! Set runoff variables to a given value - useful to set them to _FillValue
    ! ---------------------------------------------------------------------------------------

    USE multiroute                                        ! routed runoff

    IMPLICIT NONE
    REAL(SP), INTENT(IN)                  :: VAL          ! value
    ! ---------------------------------------------------------------------------------------
    MROUTE%Q_INSTNT = VAL     ! instantaneous runoff
    MROUTE%Q_ROUTED = VAL     ! routed runoff
    MROUTE%Q_ACCURATE  = VAL  ! "accurate" runoff estimate (mm day-1)

    ! (routed runoff)
    ! FUTURE         = VAL
    ! ---------------------------------------------------------------------------------------
  END SUBROUTINE SET_ROUTE

  SUBROUTINE SET_SNOW(VAL)

    ! ---------------------------------------------------------------------------------------
    ! Creator:
    ! --------
    ! Nans Addor based on Martyn Clark's INIT_STATE
    ! ---------------------------------------------------------------------------------------
    ! Purpose:
    ! --------
    ! Set snow variables to a given value - useful to set them to _FillValue
    ! ---------------------------------------------------------------------------------------

    USE multibands                                        ! elevation bands for snow modeling

    IMPLICIT NONE
    REAL(SP), INTENT(IN)                  :: VAL          ! value
    INTEGER(I4B)                          :: IBANDS       ! snow band index

    ! ---------------------------------------------------------------------------------------
    DO IBANDS=1,N_BANDS
       MBANDS(IBANDS)%SWE=VAL           ! band snowpack water equivalent (mm)
       MBANDS(IBANDS)%SNOWACCMLTN=VAL   ! new snow accumulation in band (mm day-1)
       MBANDS(IBANDS)%SNOWMELT=VAL      ! snowmelt in band (mm day-1)
       MBANDS(IBANDS)%DSWE_DT=VAL       ! rate of change of band SWE (mm day-1)
    END DO

    ! ---------------------------------------------------------------------------------------
  END SUBROUTINE SET_SNOW

end module set_all_module
