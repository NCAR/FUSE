SUBROUTINE MEANFLUXES()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Compute 0.5*(FLUX_0 + FLUX_1)
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! Fluxes in MODULE multi_flux
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE multi_flux                                        ! model fluxes
USE multistate                                        ! model states (use time step)
IMPLICIT NONE
! ---------------------------------------------------------------------------------------
M_FLUX%EFF_PPT     = 0.5_sp * (FLUX_0%EFF_PPT     + FLUX_1%EFF_PPT    )
M_FLUX%SATAREA     = 0.5_sp * (FLUX_0%SATAREA     + FLUX_1%SATAREA    )
M_FLUX%QSURF       = 0.5_sp * (FLUX_0%QSURF       + FLUX_1%QSURF      )
M_FLUX%EVAP_1A     = 0.5_sp * (FLUX_0%EVAP_1A     + FLUX_1%EVAP_1A    )
M_FLUX%EVAP_1B     = 0.5_sp * (FLUX_0%EVAP_1B     + FLUX_1%EVAP_1B    )
M_FLUX%EVAP_1      = 0.5_sp * (FLUX_0%EVAP_1      + FLUX_1%EVAP_1     )
M_FLUX%EVAP_2      = 0.5_sp * (FLUX_0%EVAP_2      + FLUX_1%EVAP_2     )
M_FLUX%RCHR2EXCS   = 0.5_sp * (FLUX_0%RCHR2EXCS   + FLUX_1%RCHR2EXCS  )
M_FLUX%TENS2FREE_1 = 0.5_sp * (FLUX_0%TENS2FREE_1 + FLUX_1%TENS2FREE_1)
M_FLUX%TENS2FREE_2 = 0.5_sp * (FLUX_0%TENS2FREE_2 + FLUX_1%TENS2FREE_2)
M_FLUX%QINTF_1     = 0.5_sp * (FLUX_0%QINTF_1     + FLUX_1%QINTF_1    )
M_FLUX%QPERC_12    = 0.5_sp * (FLUX_0%QPERC_12    + FLUX_1%QPERC_12   )
M_FLUX%QBASE_2     = 0.5_sp * (FLUX_0%QBASE_2     + FLUX_1%QBASE_2    )
M_FLUX%QBASE_2A    = 0.5_sp * (FLUX_0%QBASE_2A    + FLUX_1%QBASE_2A   )
M_FLUX%QBASE_2B    = 0.5_sp * (FLUX_0%QBASE_2B    + FLUX_1%QBASE_2B   )
M_FLUX%OFLOW_1     = 0.5_sp * (FLUX_0%OFLOW_1     + FLUX_1%OFLOW_1    )
M_FLUX%OFLOW_2     = 0.5_sp * (FLUX_0%OFLOW_2     + FLUX_1%OFLOW_2    )
M_FLUX%OFLOW_2A    = 0.5_sp * (FLUX_0%OFLOW_2A    + FLUX_1%OFLOW_2A   )
M_FLUX%OFLOW_2B    = 0.5_sp * (FLUX_0%OFLOW_2B    + FLUX_1%OFLOW_2B   )
M_FLUX%ERR_WATR_1  = 0.5_sp * (FLUX_0%ERR_WATR_1  + FLUX_1%ERR_WATR_1 )
M_FLUX%ERR_TENS_1  = 0.5_sp * (FLUX_0%ERR_TENS_1  + FLUX_1%ERR_TENS_1 )
M_FLUX%ERR_FREE_1  = 0.5_sp * (FLUX_0%ERR_FREE_1  + FLUX_1%ERR_FREE_1 )
M_FLUX%ERR_TENS_1A = 0.5_sp * (FLUX_0%ERR_TENS_1A + FLUX_1%ERR_TENS_1A)
M_FLUX%ERR_TENS_1B = 0.5_sp * (FLUX_0%ERR_TENS_1B + FLUX_1%ERR_TENS_1B)
M_FLUX%ERR_WATR_2  = 0.5_sp * (FLUX_0%ERR_WATR_2  + FLUX_1%ERR_WATR_2 )
M_FLUX%ERR_TENS_2  = 0.5_sp * (FLUX_0%ERR_TENS_2  + FLUX_1%ERR_TENS_2 )
M_FLUX%ERR_FREE_2  = 0.5_sp * (FLUX_0%ERR_FREE_2  + FLUX_1%ERR_FREE_2 )
M_FLUX%ERR_FREE_2A = 0.5_sp * (FLUX_0%ERR_FREE_2A + FLUX_1%ERR_FREE_2A)
M_FLUX%ERR_FREE_2B = 0.5_sp * (FLUX_0%ERR_FREE_2B + FLUX_1%ERR_FREE_2B)
! ---------------------------------------------------------------------------------------
END SUBROUTINE MEANFLUXES
