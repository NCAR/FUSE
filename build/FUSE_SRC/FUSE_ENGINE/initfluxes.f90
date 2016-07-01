SUBROUTINE INITFLUXES()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! Modified by Brian Henn to include snow model, 6/2013
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Set all fluxes to zero at the start of each time step
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! Fluxes in MODULE multi_flux
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model decision structures
USE model_defnames                                    ! integer model definitions
USE multi_flux                                        ! model fluxes
USE multibands                                        ! model snow bands
IMPLICIT NONE
INTEGER(I4B)                          :: ISNW         ! index for looping though SWE 
! ---------------------------------------------------------------------------------------
M_FLUX%EFF_PPT     = 0._sp; W_FLUX%EFF_PPT     = 0._sp
M_FLUX%SATAREA     = 0._sp; W_FLUX%SATAREA     = 0._sp
M_FLUX%QSURF       = 0._sp; W_FLUX%QSURF       = 0._sp
M_FLUX%EVAP_1A     = 0._sp; W_FLUX%EVAP_1A     = 0._sp
M_FLUX%EVAP_1B     = 0._sp; W_FLUX%EVAP_1B     = 0._sp
M_FLUX%EVAP_1      = 0._sp; W_FLUX%EVAP_1      = 0._sp
M_FLUX%EVAP_2      = 0._sp; W_FLUX%EVAP_2      = 0._sp
M_FLUX%RCHR2EXCS   = 0._sp; W_FLUX%RCHR2EXCS   = 0._sp
M_FLUX%TENS2FREE_1 = 0._sp; W_FLUX%TENS2FREE_1 = 0._sp
M_FLUX%TENS2FREE_2 = 0._sp; W_FLUX%TENS2FREE_2 = 0._sp
M_FLUX%QINTF_1     = 0._sp; W_FLUX%QINTF_1     = 0._sp
M_FLUX%QPERC_12    = 0._sp; W_FLUX%QPERC_12    = 0._sp
M_FLUX%QBASE_2     = 0._sp; W_FLUX%QBASE_2     = 0._sp
M_FLUX%QBASE_2A    = 0._sp; W_FLUX%QBASE_2A    = 0._sp
M_FLUX%QBASE_2B    = 0._sp; W_FLUX%QBASE_2B    = 0._sp
M_FLUX%OFLOW_1     = 0._sp; W_FLUX%OFLOW_1     = 0._sp
M_FLUX%OFLOW_2     = 0._sp; W_FLUX%OFLOW_2     = 0._sp
M_FLUX%OFLOW_2A    = 0._sp; W_FLUX%OFLOW_2A    = 0._sp
M_FLUX%OFLOW_2B    = 0._sp; W_FLUX%OFLOW_2B    = 0._sp
IF(SMODL%iSNOWM.EQ.iopt_temp_index) THEN !loop through snow model bands
 DO ISNW=1,N_BANDS
  MBANDS(ISNW)%SNOWACCMLTN  = 0._sp 
  MBANDS(ISNW)%SNOWMELT     = 0._sp 
 END DO
ENDIF
M_FLUX%ERR_WATR_1  = 0._sp; W_FLUX%ERR_WATR_1  = 0._sp
M_FLUX%ERR_TENS_1  = 0._sp; W_FLUX%ERR_TENS_1  = 0._sp
M_FLUX%ERR_FREE_1  = 0._sp; W_FLUX%ERR_FREE_1  = 0._sp
M_FLUX%ERR_TENS_1A = 0._sp; W_FLUX%ERR_TENS_1A = 0._sp
M_FLUX%ERR_TENS_1B = 0._sp; W_FLUX%ERR_TENS_1B = 0._sp
M_FLUX%ERR_WATR_2  = 0._sp; W_FLUX%ERR_WATR_2  = 0._sp
M_FLUX%ERR_TENS_2  = 0._sp; W_FLUX%ERR_TENS_2  = 0._sp
M_FLUX%ERR_FREE_2  = 0._sp; W_FLUX%ERR_FREE_2  = 0._sp
M_FLUX%ERR_FREE_2A = 0._sp; W_FLUX%ERR_FREE_2A = 0._sp
M_FLUX%ERR_FREE_2B = 0._sp; W_FLUX%ERR_FREE_2B = 0._sp
M_FLUX%CHK_TIME    = 0._sp; W_FLUX%CHK_TIME    = 0._sp
! ---------------------------------------------------------------------------------------
END SUBROUTINE INITFLUXES
