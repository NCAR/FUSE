SUBROUTINE WGT_FLUXES()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Compute the contribution of the flux over the time interval HSTATE%STEP
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! Fluxes in MODULE multi_flux
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE multiforce, ONLY: DELTIM                          ! model forcing data
USE multi_flux                                        ! model fluxes
USE multistate, ONLY: HSTATE                          ! model states (use time step)
IMPLICIT NONE
REAL(SP)                                :: WEIGHT     ! weight (sub-step / full-step)
! ---------------------------------------------------------------------------------------
! determine the contribution of the sub-step to the full step
WEIGHT = HSTATE%STEP/DELTIM
! ---------------------------------------------------------------------------------------
! update data structure
W_FLUX%EFF_PPT     = W_FLUX%EFF_PPT     + M_FLUX%EFF_PPT     * WEIGHT
W_FLUX%SATAREA     = W_FLUX%SATAREA     + M_FLUX%SATAREA     * WEIGHT
W_FLUX%QSURF       = W_FLUX%QSURF       + M_FLUX%QSURF       * WEIGHT
W_FLUX%EVAP_1A     = W_FLUX%EVAP_1A     + M_FLUX%EVAP_1A     * WEIGHT
W_FLUX%EVAP_1B     = W_FLUX%EVAP_1B     + M_FLUX%EVAP_1B     * WEIGHT
W_FLUX%EVAP_1      = W_FLUX%EVAP_1      + M_FLUX%EVAP_1      * WEIGHT
W_FLUX%EVAP_2      = W_FLUX%EVAP_2      + M_FLUX%EVAP_2      * WEIGHT
W_FLUX%RCHR2EXCS   = W_FLUX%RCHR2EXCS   + M_FLUX%RCHR2EXCS   * WEIGHT
W_FLUX%TENS2FREE_1 = W_FLUX%TENS2FREE_1 + M_FLUX%TENS2FREE_1 * WEIGHT
W_FLUX%TENS2FREE_2 = W_FLUX%TENS2FREE_2 + M_FLUX%TENS2FREE_2 * WEIGHT
W_FLUX%QINTF_1     = W_FLUX%QINTF_1     + M_FLUX%QINTF_1     * WEIGHT
W_FLUX%QPERC_12    = W_FLUX%QPERC_12    + M_FLUX%QPERC_12    * WEIGHT
W_FLUX%QBASE_2     = W_FLUX%QBASE_2     + M_FLUX%QBASE_2     * WEIGHT
W_FLUX%QBASE_2A    = W_FLUX%QBASE_2A    + M_FLUX%QBASE_2A    * WEIGHT
W_FLUX%QBASE_2B    = W_FLUX%QBASE_2B    + M_FLUX%QBASE_2B    * WEIGHT
W_FLUX%OFLOW_1     = W_FLUX%OFLOW_1     + M_FLUX%OFLOW_1     * WEIGHT
W_FLUX%OFLOW_2     = W_FLUX%OFLOW_2     + M_FLUX%OFLOW_2     * WEIGHT
W_FLUX%OFLOW_2A    = W_FLUX%OFLOW_2A    + M_FLUX%OFLOW_2A    * WEIGHT
W_FLUX%OFLOW_2B    = W_FLUX%OFLOW_2B    + M_FLUX%OFLOW_2B    * WEIGHT
W_FLUX%ERR_WATR_1  = W_FLUX%ERR_WATR_1  + M_FLUX%ERR_WATR_1  * WEIGHT
W_FLUX%ERR_TENS_1  = W_FLUX%ERR_TENS_1  + M_FLUX%ERR_TENS_1  * WEIGHT
W_FLUX%ERR_FREE_1  = W_FLUX%ERR_FREE_1  + M_FLUX%ERR_FREE_1  * WEIGHT
W_FLUX%ERR_TENS_1A = W_FLUX%ERR_TENS_1A + M_FLUX%ERR_TENS_1A * WEIGHT
W_FLUX%ERR_TENS_1B = W_FLUX%ERR_TENS_1B + M_FLUX%ERR_TENS_1B * WEIGHT
W_FLUX%ERR_WATR_2  = W_FLUX%ERR_WATR_2  + M_FLUX%ERR_WATR_2  * WEIGHT
W_FLUX%ERR_TENS_2  = W_FLUX%ERR_TENS_2  + M_FLUX%ERR_TENS_2  * WEIGHT
W_FLUX%ERR_FREE_2  = W_FLUX%ERR_FREE_2  + M_FLUX%ERR_FREE_2  * WEIGHT
W_FLUX%ERR_FREE_2A = W_FLUX%ERR_FREE_2A + M_FLUX%ERR_FREE_2A * WEIGHT
W_FLUX%ERR_FREE_2B = W_FLUX%ERR_FREE_2B + M_FLUX%ERR_FREE_2B * WEIGHT
W_FLUX%CHK_TIME    = W_FLUX%CHK_TIME    + WEIGHT
! ---------------------------------------------------------------------------------------
END SUBROUTINE WGT_FLUXES
