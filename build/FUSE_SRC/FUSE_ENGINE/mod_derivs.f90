SUBROUTINE MOD_DERIVS()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! Modified by to include snow model by Brian Henn, 6/2013
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! compute the derivative (dydx) of all model states (y) at time (x)
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multi_flux -- populate structure M_FLUX%(*)
! MODULE multistate -- populate structure DY_DT%(*)
USE model_numerix, ONLY: NUM_FUNCS  ! (number of function evaluations)
! ---------------------------------------------------------------------------------------
! (1) COMPUTE FLUXES
! ---------------------------------------------------------------------------------------
CALL QSATEXCESS()  ! compute the saturated area and surface runoff
CALL EVAP_UPPER()  ! compute evaporation from the upper layer
CALL EVAP_LOWER()  ! compute evaporation from the lower layer
CALL QINTERFLOW()  ! compute interflow from free water in the upper layer
CALL QPERCOLATE()  ! compute percolation from the upper to lower soil layers
CALL Q_BASEFLOW()  ! compute baseflow from the lower soil layer
CALL Q_MISSCELL()  ! compute miscellaneous fluxes (NOTE: need sat area, evap, and perc)
! ---------------------------------------------------------------------------------------
! (2) COMPUTE DERIVATIVES FOR EACH OF THE MODEL STATES
! ---------------------------------------------------------------------------------------
CALL MSTATE_EQN()
! ---------------------------------------------------------------------------------------
! (3) KEEP TRACK OF THE NUMBER OF FUNCTION CALLS
! ---------------------------------------------------------------------------------------
NUM_FUNCS = NUM_FUNCS + 1  ! NUM_FUNCS is shared in module model_numerix
! ---------------------------------------------------------------------------------------
END SUBROUTINE MOD_DERIVS
