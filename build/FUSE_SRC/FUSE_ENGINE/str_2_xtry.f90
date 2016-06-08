MODULE STR_2_XTRY_MODULE
IMPLICIT NONE
CONTAINS
SUBROUTINE STR_2_XTRY(TMPSTR,X_TRY)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Copy model states into vector X_TRY
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! Numerical Recipes data types
USE model_defn, ONLY: CSTATE,NSTATE                   ! model definitions
USE model_defnames
USE multistate, ONLY: STATEV                          ! model state structure
IMPLICIT NONE
! input
TYPE(STATEV), INTENT(IN)               :: TMPSTR      ! temporary state structure
! output
REAL(SP), DIMENSION(:), INTENT(OUT)    :: X_TRY       ! vector of model states
! internal
INTEGER(I4B)                           :: ISTT        ! loop through model states
! ---------------------------------------------------------------------------------------
DO ISTT=1,NSTATE
 SELECT CASE(CSTATE(ISTT)%iSNAME)
  CASE (iopt_TENS1A); X_TRY(ISTT) = TMPSTR%TENS_1A
  CASE (iopt_TENS1B); X_TRY(ISTT) = TMPSTR%TENS_1B
  CASE (iopt_TENS_1); X_TRY(ISTT) = TMPSTR%TENS_1
  CASE (iopt_FREE_1); X_TRY(ISTT) = TMPSTR%FREE_1
  CASE (iopt_WATR_1); X_TRY(ISTT) = TMPSTR%WATR_1
  CASE (iopt_TENS_2); X_TRY(ISTT) = TMPSTR%TENS_2
  CASE (iopt_FREE2A); X_TRY(ISTT) = TMPSTR%FREE_2A
  CASE (iopt_FREE2B); X_TRY(ISTT) = TMPSTR%FREE_2B
  CASE (iopt_WATR_2); X_TRY(ISTT) = TMPSTR%WATR_2
 END SELECT
END DO  ! istt
! ---------------------------------------------------------------------------------------
END SUBROUTINE STR_2_XTRY
END MODULE STR_2_XTRY_MODULE
