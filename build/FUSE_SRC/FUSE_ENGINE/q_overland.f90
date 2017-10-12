SUBROUTINE Q_OVERLAND()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! --------
! History
! 5 June 2013 AD: Modified by David McInerney to merge array loop operations
! 5 June 2013 AD: Modified by Dmitri Kavetski to avoid zero-element operations
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Computes the time delay in runoff in a basin (places runoff in future time steps)
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multiroute -- places runoff in array FUTURE(:)RUNOFF
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition structure
USE model_defnames
USE multiparam                                        ! model parameters
USE multi_flux                                        ! model fluxes
USE multiroute                                        ! routed runoff
USE multiforce, only:  MFORCE                         ! forcing - for debugging only

IMPLICIT NONE
INTEGER(I4B)                           :: NTDH        ! maximum number of future time steps
INTEGER(I4B)                           :: JTIM        ! (loop through future time steps)
REAL(SP), PARAMETER                    :: SNEG=-1.e-5 ! small negative number, used for checking
LOGICAL, PARAMETER                     :: USE_NTDH_NEED=.TRUE. ! flag to use NTDH_NEED to reduce array operations (loop length)
! ---------------------------------------------------------------------------------------
! compute total runoff (sum of surface runoff, overflow, interflow, and baseflow
MROUTE%Q_INSTNT = W_FLUX%QSURF + W_FLUX%OFLOW_1 + W_FLUX%QINTF_1 + W_FLUX%OFLOW_2 + W_FLUX%QBASE_2

if (W_FLUX%QSURF.lt.SNEG .or. W_FLUX%OFLOW_1.lt.SNEG .or. W_FLUX%QINTF_1.lt.SNEG .or. &
    W_FLUX%OFLOW_2.lt.SNEG .or. W_FLUX%QBASE_2.lt.SNEG) THEN

    !PRINT *, 'W_FLUX%QSURF = ', W_FLUX%QSURF
    !PRINT *, 'W_FLUX%OFLOW_1 = ', W_FLUX%OFLOW_1
    !PRINT *, 'W_FLUX%QINTF_1 = ', W_FLUX%QINTF_1
    !PRINT *, 'W_FLUX%OFLOW_2 = ', W_FLUX%OFLOW_2
    !PRINT *, 'W_FLUX%QBASE_2 = ', W_FLUX%QBASE_2
    !PRINT *, 'MROUTE%Q_INSTNT = ', MROUTE%Q_INSTNT

    !stop 'negative flux in q_overland'

  END IF
! ---------------------------------------------------------------------------------------
SELECT CASE(SMODL%iQ_TDH)
 CASE(iopt_rout_gamma)             ! use a Gamma distribution with shape parameter = 2.5
  NTDH = SIZE(DPARAM%FRAC_FUTURE)  ! maximum number of future time steps
  MROUTE%Q_ROUTED = FUTURE(1) + MROUTE%Q_INSTNT * DPARAM%FRAC_FUTURE(1)
  DO JTIM=2,MERGE(DPARAM%NTDH_NEED,NTDH,USE_NTDH_NEED) ! update and move array of states within the routing convolution
    FUTURE(JTIM-1) = FUTURE(JTIM) + MROUTE%Q_INSTNT * DPARAM%FRAC_FUTURE(JTIM)
  END DO
  FUTURE(JTIM-1) = 0._sp  ! last element (just in case) - the rest are never accessed (treated as 0)
 CASE(iopt_no_routing)             ! no routing
  MROUTE%Q_ROUTED   = MROUTE%Q_INSTNT
 CASE DEFAULT                      ! check for errors
  print *, "SMODL%iQ_TDH must be either iopt_rout_gamma or iopt_no_routing"
  STOP
END SELECT
! ---------------------------------------------------------------------------------------
END SUBROUTINE Q_OVERLAND
