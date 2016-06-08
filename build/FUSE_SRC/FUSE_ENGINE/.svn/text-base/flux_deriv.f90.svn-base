SUBROUTINE FLUX_DERIV(J,DS)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Compute the flux derivatives, used in calculating time-step average fluxes in the
!  semi-implicit Euler scheme
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! Finite-difference fluxes in MODULE multi_flux
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE nrutil, ONLY : nrerror                            ! error control
USE model_defn, ONLY: NSTATE,N_FLUX,C_FLUX            ! number of state variables
USE multi_flux, ONLY: FLUX_0,M_FLUX,FDFLUX            ! model fluxes
IMPLICIT NONE
! input
INTEGER(I4B), INTENT(IN)               :: J           ! index of state variable
REAL(SP), INTENT(IN)                   :: DS          ! difference in state variable 
! internal
INTEGER(I4B)                           :: IFLUX       ! loop thru fluxes
INTEGER(I4B)                           :: IERR        ! error code for the allocate statement
! ---------------------------------------------------------------------------------------
! make sure that the finite-difference flux structure is allocated
IF (.NOT.ASSOCIATED(FDFLUX)) THEN
 ALLOCATE(FDFLUX(NSTATE), STAT=IERR ) ! NSTATE in structure model_defn
 IF (IERR.NE.0) CALL NRERROR('flux_deriv: problem allocating fdflux')
ENDIF
! make sure that there are some fluxes
IF (N_FLUX.EQ.0) CALL NRERROR('flux_deriv: number of fluxes is zero')
! ---------------------------------------------------------------------------------------
DO IFLUX=1,N_FLUX
 SELECT CASE(TRIM(C_FLUX(IFLUX)%FNAME))
  CASE('EFF_PPT')    ; FDFLUX(J)%EFF_PPT     = (M_FLUX%EFF_PPT     - FLUX_0%EFF_PPT)     / DS
  CASE('EVAP_1A')    ; FDFLUX(J)%EVAP_1A     = (M_FLUX%EVAP_1A     - FLUX_0%EVAP_1A)     / DS
  CASE('EVAP_1B')    ; FDFLUX(J)%EVAP_1B     = (M_FLUX%EVAP_1B     - FLUX_0%EVAP_1B)     / DS
  CASE('EVAP_1')     ; FDFLUX(J)%EVAP_1      = (M_FLUX%EVAP_1      - FLUX_0%EVAP_1)      / DS
  CASE('EVAP_2')     ; FDFLUX(J)%EVAP_2      = (M_FLUX%EVAP_2      - FLUX_0%EVAP_2)      / DS
  CASE('RCHR2EXCS')  ; FDFLUX(J)%RCHR2EXCS   = (M_FLUX%RCHR2EXCS   - FLUX_0%RCHR2EXCS)   / DS
  CASE('TENS2FREE_1'); FDFLUX(J)%TENS2FREE_1 = (M_FLUX%TENS2FREE_1 - FLUX_0%TENS2FREE_1) / DS
  CASE('TENS2FREE_2'); FDFLUX(J)%TENS2FREE_2 = (M_FLUX%TENS2FREE_2 - FLUX_0%TENS2FREE_2) / DS
  CASE('QSURF')      ; FDFLUX(J)%QSURF       = (M_FLUX%QSURF       - FLUX_0%QSURF)       / DS
  CASE('QPERC_12')   ; FDFLUX(J)%QPERC_12    = (M_FLUX%QPERC_12    - FLUX_0%QPERC_12)    / DS
  CASE('QINTF_1')    ; FDFLUX(J)%QINTF_1     = (M_FLUX%QINTF_1     - FLUX_0%QINTF_1)     / DS
  CASE('QBASE_2')    ; FDFLUX(J)%QBASE_2     = (M_FLUX%QBASE_2     - FLUX_0%QBASE_2)     / DS
  CASE('QBASE_2A')   ; FDFLUX(J)%QBASE_2A    = (M_FLUX%QBASE_2A    - FLUX_0%QBASE_2A)    / DS
  CASE('QBASE_2B')   ; FDFLUX(J)%QBASE_2B    = (M_FLUX%QBASE_2B    - FLUX_0%QBASE_2B)    / DS
  CASE('OFLOW_1')    ; FDFLUX(J)%OFLOW_1     = (M_FLUX%OFLOW_1     - FLUX_0%OFLOW_1)     / DS
  CASE('OFLOW_2')    ; FDFLUX(J)%OFLOW_2     = (M_FLUX%OFLOW_2     - FLUX_0%OFLOW_2)     / DS
  CASE('OFLOW_2A')   ; FDFLUX(J)%OFLOW_2A    = (M_FLUX%OFLOW_2A    - FLUX_0%OFLOW_2A)    / DS
  CASE('OFLOW_2B')   ; FDFLUX(J)%OFLOW_2B    = (M_FLUX%OFLOW_2B    - FLUX_0%OFLOW_2B)    / DS
  CASE DEFAULT       ; CALL NRERROR('flux_deriv: cannot find desired state')
 END SELECT
END DO ! (loop through fluxes)
! ---------------------------------------------------------------------------------------
END SUBROUTINE FLUX_DERIV
