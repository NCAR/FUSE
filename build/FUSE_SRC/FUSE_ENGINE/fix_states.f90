SUBROUTINE FIX_STATES(DT,ERROR_FLAG)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Ensure states are within bounds, and disaggregate fluxes if necessary
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multistate -- populates the MODULE multistate with derivatives DY_DT%(*)
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition structure
USE model_defnames
USE multiparam                                        ! model parameters
USE multiforce                                        ! model forcing data
USE multistate                                        ! model states
USE multi_flux                                        ! model fluxes
USE model_numerix                                     ! model numerix
IMPLICIT NONE
! input/output
REAL(SP), INTENT(IN)                   :: DT          ! time step
LOGICAL(LGT), INTENT(OUT)              :: ERROR_FLAG  ! .TRUE. if extrapolation error
! internal
REAL(SP)                               :: XMIN        ! very small number
INTEGER(I4B)                           :: ISTT        ! loop through model states
REAL(SP)                               :: ERROR_LOSS  ! error (L/T)
REAL(SP)                               :: TOTAL_LOSS  ! total loss (L/T)
! ---------------------------------------------------------------------------------------
ERROR_FLAG=.FALSE.   ! initialize with no extrapolation error
! ---------------------------------------------------------------------------------------
XMIN = FRACSTATE_MIN ! used to avoid zero derivatives
! ---------------------------------------------------------------------------------------
DO ISTT=1,NSTATE
 if (M_FLUX%QSURF.LT.0._sp) print *, 'start ', desc_int2str(cstate(istt)%isname), M_FLUX%QSURF
 ERROR_LOSS = 0._SP ! initialize state error
 SELECT CASE(CSTATE(ISTT)%iSNAME)
  ! ---------------------------------------------------------------------------------------
  ! (1) FIX STATES IN THE UPPER LAYER
  ! -------------------------------------------------------------------------------------
  CASE (iopt_TENS1A)
   IF (ESTATE%TENS_1A.LT.XMIN*DPARAM%MAXTENS_1A) THEN   ! too much drainage
    ERROR_LOSS         = (ESTATE%TENS_1A - XMIN*DPARAM%MAXTENS_1A)/DT  ! error (L/T)
    TOTAL_LOSS         = M_FLUX%QSURF + M_FLUX%EVAP_1A                 ! total loss (L/T)
    M_FLUX%QSURF       = M_FLUX%QSURF     + (M_FLUX%QSURF  /TOTAL_LOSS)*ERROR_LOSS
    M_FLUX%EVAP_1A     = M_FLUX%EVAP_1A   + (M_FLUX%EVAP_1A/TOTAL_LOSS)*ERROR_LOSS
    ESTATE%TENS_1A     = XMIN*DPARAM%MAXTENS_1A   ! (correct state)
    ERROR_FLAG         = .TRUE.
   ENDIF
   IF (ESTATE%TENS_1A.GT.DPARAM%MAXTENS_1A) THEN        ! too much input
    ERROR_LOSS         = (ESTATE%TENS_1A - DPARAM%MAXTENS_1A)/DT 
    M_FLUX%RCHR2EXCS   = M_FLUX%RCHR2EXCS + ERROR_LOSS
    ESTATE%TENS_1A     = DPARAM%MAXTENS_1A        ! (correct state)
    ESTATE%TENS_1B     = BSTATE%TENS_1B +   &     ! (correct subsequent states)
                         (M_FLUX%RCHR2EXCS - M_FLUX%EVAP_1B - M_FLUX%TENS2FREE_1)*DT
    ERROR_FLAG         = .TRUE.
   ENDIF
   M_FLUX%ERR_TENS_1A = ERROR_LOSS
  ! -------------------------------------------------------------------------------------
  CASE (iopt_TENS1B)
   IF (ESTATE%TENS_1B.LT.XMIN*DPARAM%MAXTENS_1B) THEN   ! too much drainage
    ERROR_LOSS         = (ESTATE%TENS_1B - XMIN*DPARAM%MAXTENS_1B)/DT
    M_FLUX%EVAP_1B     = M_FLUX%EVAP_1B + ERROR_LOSS
    ESTATE%TENS_1B     = XMIN*DPARAM%MAXTENS_1B   ! (correct state)
    ERROR_FLAG         = .TRUE.
   ENDIF
   IF (ESTATE%TENS_1B.GT.DPARAM%MAXTENS_1B) THEN        ! too much input 
    ERROR_LOSS         = (ESTATE%TENS_1B - DPARAM%MAXTENS_1B)/DT
    M_FLUX%TENS2FREE_1 = M_FLUX%TENS2FREE_1 + ERROR_LOSS
    ESTATE%TENS_1B     = DPARAM%MAXTENS_1B        ! (correct state)
    ESTATE%FREE_1      = BSTATE%FREE_1 + &        ! (correct subsequent states)
                          (M_FLUX%TENS2FREE_1 - M_FLUX%QPERC_12 - M_FLUX%QINTF_1 - M_FLUX%OFLOW_1)*DT
    ERROR_FLAG         = .TRUE.
   ENDIF
   M_FLUX%ERR_TENS_1B = ERROR_LOSS
  ! -------------------------------------------------------------------------------------
  CASE (iopt_TENS_1)
   IF (ESTATE%TENS_1.LT.XMIN*DPARAM%MAXTENS_1) THEN    ! too much drainage
    ERROR_LOSS         = (ESTATE%TENS_1 - XMIN*DPARAM%MAXTENS_1)/DT    ! error (L/T)
    TOTAL_LOSS         = M_FLUX%QSURF + M_FLUX%EVAP_1                  ! total loss (L/T)
    M_FLUX%QSURF       = M_FLUX%QSURF  + (M_FLUX%QSURF /TOTAL_LOSS)*ERROR_LOSS
    M_FLUX%EVAP_1      = M_FLUX%EVAP_1 + (M_FLUX%EVAP_1/TOTAL_LOSS)*ERROR_LOSS
    ESTATE%TENS_1      = XMIN*DPARAM%MAXTENS_1    ! (correct state)
    ERROR_FLAG         = .TRUE.
   ENDIF
   IF (ESTATE%TENS_1.GT.DPARAM%MAXTENS_1) THEN         ! too much input
    ERROR_LOSS         = (ESTATE%TENS_1 - DPARAM%MAXTENS_1)/DT
    M_FLUX%TENS2FREE_1 = M_FLUX%TENS2FREE_1 + (ESTATE%TENS_1 - DPARAM%MAXTENS_1)/DT
    ESTATE%TENS_1      = DPARAM%MAXTENS_1         ! (correct state)
    ESTATE%FREE_1      = BSTATE%FREE_1 + &        ! (correct subsequent states)
                          (M_FLUX%TENS2FREE_1 - M_FLUX%QPERC_12 - M_FLUX%QINTF_1 - M_FLUX%OFLOW_1)*DT
    ERROR_FLAG         = .TRUE.
   ENDIF
   M_FLUX%ERR_TENS_1 = ERROR_LOSS
  ! -------------------------------------------------------------------------------------
  CASE (iopt_FREE_1)
   IF (ESTATE%FREE_1.LT.XMIN*DPARAM%MAXFREE_1) THEN    ! too much drainage
    ERROR_LOSS         = (ESTATE%FREE_1 - XMIN*DPARAM%MAXFREE_1)/DT    ! error (L/T)
    TOTAL_LOSS         = M_FLUX%QPERC_12 + M_FLUX%QINTF_1              ! total loss (L/T)
    M_FLUX%QPERC_12    = M_FLUX%QPERC_12 + (M_FLUX%QPERC_12/TOTAL_LOSS)*ERROR_LOSS
    M_FLUX%QINTF_1     = M_FLUX%QINTF_1  + (M_FLUX%QINTF_1 /TOTAL_LOSS)*ERROR_LOSS
    ESTATE%FREE_1      = XMIN*DPARAM%MAXFREE_1    ! (correct state)
    ! correct subsequent states (deal appropriately with percolation)
    ! NOTE: do this here because only necessary to make corrections if M_FLUX%QPERC_12 changes
    SELECT CASE(SMODL%iARCH2)
     CASE(iopt_tens2pll_2) ! tension reservoir plus two parallel tanks
      ! fix overflow fluxes
      M_FLUX%TENS2FREE_2 = MAX(0._SP, M_FLUX%QPERC_12*(1._SP-MPARAM%PERCFRAC) - (DPARAM%MAXTENS_2  - BSTATE%TENS_2 )/DT)
      M_FLUX%OFLOW_2A    = MAX(0._SP, (M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP) &
                                          - (DPARAM%MAXFREE_2A - BSTATE%FREE_2A)/DT)  
      M_FLUX%OFLOW_2B    = MAX(0._SP, (M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP) &
                                          - (DPARAM%MAXFREE_2B - BSTATE%FREE_2B)/DT)
      M_FLUX%OFLOW_2     = M_FLUX%OFLOW_2A + M_FLUX%OFLOW_2B
      ! fix states
      ESTATE%TENS_2    = BSTATE%TENS_2 + &
                         (M_FLUX%QPERC_12*(1._SP-MPARAM%PERCFRAC) - M_FLUX%EVAP_2 - M_FLUX%TENS2FREE_2)*DT
      ESTATE%FREE_2A   = BSTATE%FREE_2A + &
                         (M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP - M_FLUX%QBASE_2A &
                             - M_FLUX%OFLOW_2A)*DT
      ESTATE%FREE_2B   = BSTATE%FREE_2B + &
                         (M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP - M_FLUX%QBASE_2B &
                             - M_FLUX%OFLOW_2B)*DT
     CASE(iopt_unlimfrc_2,iopt_unlimpow_2,iopt_fixedsiz_2) ! single state
      ! NOTE: M_FLUX%OFLOW_2 and M_FLUX%EVAP_2 only calculated for 'fixedsiz_2'
      ! fix overflow
      IF (SMODL%iARCH2.EQ.iopt_fixedsiz_2) &
       M_FLUX%OFLOW_2     = MAX(0._SP, M_FLUX%QPERC_12 - (MPARAM%MAXWATR_2 - BSTATE%WATR_2)/DT)
      ! fix states
      ESTATE%WATR_2    = BSTATE%WATR_2 + &
                         (M_FLUX%QPERC_12 - M_FLUX%EVAP_2 - M_FLUX%QBASE_2 - M_FLUX%OFLOW_2)*DT  
     CASE DEFAULT; stop ' SMODL%iARCH2 must be iopt_tens2pll_2, iopt_unlimfrc_2, iopt_unlimpow_2 or iopt_fixedsiz_2 '
    END SELECT  ! deal with modified percolation of water to the lower layer
    ERROR_FLAG         = .TRUE.
   ENDIF
   IF (ESTATE%FREE_1.GT.DPARAM%MAXFREE_1) THEN            ! too much input
    ERROR_LOSS         = (ESTATE%FREE_1 - DPARAM%MAXFREE_1)/DT
    M_FLUX%OFLOW_1     = M_FLUX%OFLOW_1 + ERROR_LOSS
    ESTATE%FREE_1      = DPARAM%MAXFREE_1         ! (correct state)
    ERROR_FLAG         = .TRUE.
   ENDIF
   M_FLUX%ERR_FREE_1 = ERROR_LOSS
  ! -------------------------------------------------------------------------------------
  CASE (iopt_WATR_1)
   IF (ESTATE%WATR_1.LT.XMIN*MPARAM%MAXWATR_1) THEN    ! too much drainage
    ERROR_LOSS         = (ESTATE%WATR_1 - XMIN*MPARAM%MAXWATR_1)/DT    ! error (L/T)
    TOTAL_LOSS         = M_FLUX%QSURF + M_FLUX%EVAP_1 + M_FLUX%QPERC_12 + M_FLUX%QINTF_1
    M_FLUX%QSURF       = M_FLUX%QSURF    + (M_FLUX%QSURF   /TOTAL_LOSS)*ERROR_LOSS
    M_FLUX%EVAP_1      = M_FLUX%EVAP_1   + (M_FLUX%EVAP_1  /TOTAL_LOSS)*ERROR_LOSS
    M_FLUX%QINTF_1     = M_FLUX%QINTF_1  + (M_FLUX%QINTF_1 /TOTAL_LOSS)*ERROR_LOSS
    M_FLUX%QPERC_12    = M_FLUX%QPERC_12 + (M_FLUX%QPERC_12/TOTAL_LOSS)*ERROR_LOSS
    ESTATE%WATR_1      = XMIN*MPARAM%MAXWATR_1    ! (correct state)
    ! correct subsequent states (deal appropriately with percolation)
    ! NOTE: do this here because only necessary to make corrections if M_FLUX%QPERC_12 changes
    SELECT CASE(SMODL%iARCH2)
     CASE(iopt_tens2pll_2) ! tension reservoir plus two parallel tanks
      ! fix overflow fluxes
      M_FLUX%TENS2FREE_2 = MAX(0._SP, M_FLUX%QPERC_12*(1._SP-MPARAM%PERCFRAC) - (DPARAM%MAXTENS_2  - BSTATE%TENS_2 )/DT)
      M_FLUX%OFLOW_2A    = MAX(0._SP, (M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP) &
                                          - (DPARAM%MAXFREE_2A - BSTATE%FREE_2A)/DT)  
      M_FLUX%OFLOW_2B    = MAX(0._SP, (M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP) &
                                          - (DPARAM%MAXFREE_2B - BSTATE%FREE_2B)/DT)
      M_FLUX%OFLOW_2     = M_FLUX%OFLOW_2A + M_FLUX%OFLOW_2B
      ! fix states
      ESTATE%TENS_2    = BSTATE%TENS_2 + &
                         (M_FLUX%QPERC_12*(1._SP-MPARAM%PERCFRAC) - M_FLUX%EVAP_2 - M_FLUX%TENS2FREE_2)*DT
      ESTATE%FREE_2A   = BSTATE%FREE_2A + &
                         (M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP - M_FLUX%QBASE_2A &
                             - M_FLUX%OFLOW_2A)*DT
      ESTATE%FREE_2B   = BSTATE%FREE_2B + &
                         (M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP - M_FLUX%QBASE_2B &
                             - M_FLUX%OFLOW_2B)*DT
     CASE(iopt_unlimfrc_2,iopt_unlimpow_2,iopt_fixedsiz_2) ! single state
      ! NOTE: M_FLUX%OFLOW_2 and M_FLUX%EVAP_2 only calculated for 'fixedsiz_2'
      ! fix overflow
      IF (SMODL%iARCH2.EQ.iopt_fixedsiz_2) &
       M_FLUX%OFLOW_2     = MAX(0._SP, M_FLUX%QPERC_12 - (MPARAM%MAXWATR_2 - BSTATE%WATR_2)/DT)
      ! fix states
      ESTATE%WATR_2    = BSTATE%WATR_2 + &
                         (M_FLUX%QPERC_12 - M_FLUX%EVAP_2 - M_FLUX%QBASE_2 - M_FLUX%OFLOW_2)*DT
     CASE DEFAULT; stop ' SMODL%iARCH2 must be iopt_tens2pll_2, iopt_unlimfrc_2, iopt_unlimpow_2 or iopt_fixedsiz_2 '
    END SELECT  ! deal with modified percolation of water to the lower layer
    ERROR_FLAG         = .TRUE.
   ENDIF
   IF (ESTATE%WATR_1.GT.MPARAM%MAXWATR_1) THEN            ! too much input
    ERROR_LOSS         = (ESTATE%WATR_1 - MPARAM%MAXWATR_1)/DT
    M_FLUX%OFLOW_1     = M_FLUX%OFLOW_1 + ERROR_LOSS
    ESTATE%WATR_1      = MPARAM%MAXWATR_1         ! (correct state)
    ERROR_FLAG         = .TRUE.
   ENDIF
   M_FLUX%ERR_WATR_1 = ERROR_LOSS
  ! -------------------------------------------------------------------------------------
  ! (2) FIX STATES IN THE LOWER LAYER
  ! -------------------------------------------------------------------------------------
  CASE (iopt_TENS_2)
   IF (ESTATE%TENS_2.LT.XMIN*DPARAM%MAXTENS_2) THEN    ! too much drainage
    ERROR_LOSS         = (ESTATE%TENS_2 - XMIN*DPARAM%MAXTENS_2)/DT
    M_FLUX%EVAP_2      = M_FLUX%EVAP_2 + ERROR_LOSS
    ESTATE%TENS_2      = XMIN*DPARAM%MAXTENS_2    ! (correct state)
    ERROR_FLAG         = .TRUE.
   ENDIF
   IF (ESTATE%TENS_2.GT.DPARAM%MAXTENS_2) THEN         ! too much input
    ERROR_LOSS         = (ESTATE%TENS_2 - DPARAM%MAXTENS_2)/DT
    M_FLUX%TENS2FREE_2 = M_FLUX%TENS2FREE_2 + ERROR_LOSS
    ESTATE%TENS_2      = DPARAM%MAXTENS_2         ! (correct state)
    ! ** correct subsequent states (NOTE: 2 parallel tanks always coupled with a tension store)
    ! fix overflow fluxes 
    M_FLUX%OFLOW_2A    = MAX(0._SP, (M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP) &
                                        - (DPARAM%MAXFREE_2A - BSTATE%FREE_2A)/DT)  
    M_FLUX%OFLOW_2B    = MAX(0._SP, (M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP) &
                                        - (DPARAM%MAXFREE_2B - BSTATE%FREE_2B)/DT)
    M_FLUX%OFLOW_2     = M_FLUX%OFLOW_2A + M_FLUX%OFLOW_2B
    ! fix states
    ESTATE%FREE_2A     = BSTATE%FREE_2A + &
                         (M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP &
                             - M_FLUX%QBASE_2A - M_FLUX%OFLOW_2A)*DT
    ESTATE%FREE_2B     = BSTATE%FREE_2B + &
                         (M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP &
                             - M_FLUX%QBASE_2B - M_FLUX%OFLOW_2B)*DT
    ERROR_FLAG         = .TRUE.
   ENDIF
   M_FLUX%ERR_TENS_2 = ERROR_LOSS
  ! -------------------------------------------------------------------------------------
  CASE (iopt_FREE2A)
   IF (ESTATE%FREE_2A.LT.XMIN*DPARAM%MAXFREE_2A) THEN    ! too much drainage
    ERROR_LOSS         = (ESTATE%FREE_2A - XMIN*DPARAM%MAXFREE_2A)/DT
    M_FLUX%QBASE_2A    = M_FLUX%QBASE_2A + ERROR_LOSS
    ESTATE%FREE_2A     = XMIN*DPARAM%MAXFREE_2A   ! (correct state)
    ERROR_FLAG         = .TRUE.
   ENDIF
   IF (ESTATE%FREE_2A.GT.DPARAM%MAXFREE_2A) THEN         ! too much input
    ERROR_LOSS         = (ESTATE%FREE_2A - DPARAM%MAXFREE_2A)/DT
    M_FLUX%OFLOW_2A    = M_FLUX%OFLOW_2A + ERROR_LOSS
    ESTATE%FREE_2A     = DPARAM%MAXFREE_2A        ! (correct state)
    ERROR_FLAG         = .TRUE.
   ENDIF
   M_FLUX%ERR_FREE_2A = ERROR_LOSS
  ! -------------------------------------------------------------------------------------
  CASE (iopt_FREE2B)
   IF (ESTATE%FREE_2B.LT.XMIN*DPARAM%MAXFREE_2B) THEN    ! too much drainage
    ERROR_LOSS         = (ESTATE%FREE_2B - XMIN*DPARAM%MAXFREE_2B)/DT
    M_FLUX%QBASE_2B    = M_FLUX%QBASE_2B + ERROR_LOSS
    ESTATE%FREE_2B     = XMIN*DPARAM%MAXFREE_2B  ! (correct state)
    ERROR_FLAG         = .TRUE.
   ENDIF
   IF (ESTATE%FREE_2B.GT.DPARAM%MAXFREE_2B) THEN         ! too much input
    ERROR_LOSS         = (ESTATE%FREE_2B - DPARAM%MAXFREE_2B)/DT
    M_FLUX%OFLOW_2B    = M_FLUX%OFLOW_2B + ERROR_LOSS
    ESTATE%FREE_2B     = DPARAM%MAXFREE_2B       ! (correct state)
    ERROR_FLAG         = .TRUE.
   ENDIF
   M_FLUX%ERR_FREE_2B = ERROR_LOSS
  ! -------------------------------------------------------------------------------------
  CASE (iopt_WATR_2)
   IF (ESTATE%WATR_2.LT.XMIN*MPARAM%MAXWATR_2) THEN    ! too much drainage
    ERROR_LOSS         = (ESTATE%WATR_2 - XMIN*MPARAM%MAXWATR_1)/DT    ! error (L/T)
    TOTAL_LOSS         = M_FLUX%EVAP_2 + M_FLUX%QBASE_2
    M_FLUX%EVAP_2      = M_FLUX%EVAP_2  + (M_FLUX%EVAP_2 /TOTAL_LOSS)*ERROR_LOSS
    M_FLUX%QBASE_2     = M_FLUX%QBASE_2 + (M_FLUX%QBASE_2/TOTAL_LOSS)*ERROR_LOSS
    ESTATE%WATR_2      = XMIN*MPARAM%MAXWATR_2   ! (correct state)
    ERROR_FLAG         = .TRUE.
   ENDIF
   IF (ESTATE%FREE_2B.GT.DPARAM%MAXFREE_2B) THEN
    ERROR_LOSS         = (ESTATE%WATR_2 - MPARAM%MAXWATR_2)/DT 
    M_FLUX%OFLOW_2     = M_FLUX%OFLOW_2 + ERROR_LOSS
    ESTATE%WATR_2      = MPARAM%MAXWATR_2        ! (correct state)
    ERROR_FLAG         = .TRUE.
   ENDIF
   M_FLUX%ERR_WATR_2 = ERROR_LOSS
  CASE DEFAULT; STOP ' cannot find state in fix_states() '
 END SELECT  ! select state variable for processing
 if (M_FLUX%QSURF.LT.0._sp) print *, 'end ', desc_int2str(cstate(istt)%isname), M_FLUX%QSURF
END DO     ! loop through state variables
! ---------------------------------------------------------------------------------------
! compute derived fluxes, if necessary
IF (SMODL%iARCH2.EQ.iopt_tens2pll_2) THEN ! tension reservoir plus two parallel tanks
 M_FLUX%QBASE_2 = M_FLUX%QBASE_2A + M_FLUX%QBASE_2B
 M_FLUX%OFLOW_2 = M_FLUX%OFLOW_2A + M_FLUX%OFLOW_2B
ENDIF
! ---------------------------------------------------------------------------------------
END SUBROUTINE FIX_STATES
