MODULE DISAGGFLUX_MODULE
IMPLICIT NONE
CONTAINS
SUBROUTINE DISAGGFLUX(DELS,EFLAG)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Disaggregate fluxes for the semi-implicit Euler method
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! Fluxes in MODULE multi_flux (M_FLUX)
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE nrutil, ONLY : nrerror                            ! error control
USE model_defn, ONLY: SMODL,&                         ! identify modelling decisions
                      C_FLUX,N_FLUX, &                ! loop through the fluxes
                      CSTATE,NSTATE                   ! loop through the states
USE model_defnames
USE multiforce, ONLY: MFORCE                          ! model forcing data
USE multi_flux, ONLY: FLUX_0,M_FLUX,FDFLUX            ! model fluxes
USE multiparam, ONLY: MPARAM                          ! model parameters
IMPLICIT NONE
! input/output
REAL(SP), DIMENSION(:), INTENT(IN)     :: DELS        ! difference in state vector
LOGICAL(LGT), INTENT(OUT)              :: EFLAG       ! error flag for unusual flux
! internal
INTEGER(I4B)                           :: IFLUX       ! loop thru fluxes
INTEGER(I4B)                           :: ISTT        ! loop through states
REAL(SP), PARAMETER                    :: ZERO=0._SP  ! zero
REAL(SP)                               :: IN_FLUX     ! influx to a given bucket
REAL(SP)                               :: TOTEVAP     ! total evaporation
! ---------------------------------------------------------------------------------------
! make sure that the finite-difference flux structure is allocated
IF (.NOT.ASSOCIATED(FDFLUX)) CALL NRERROR('disaggflux: fdflux is not allocated')
IF (N_FLUX.EQ.0) CALL NRERROR('disaggflux: no fluxes identified')
EFLAG=.FALSE.  ! (initialize error flag)
! ---------------------------------------------------------------------------------------
DO IFLUX=1,N_FLUX
 ! --------------------------------------------------------------------------------------
 ! (1) DISAGGREGATE FLUXES
 ! --------------------------------------------------------------------------------------
 SELECT CASE(TRIM(C_FLUX(IFLUX)%FNAME))
  CASE('EFF_PPT')    ; M_FLUX%EFF_PPT     = FLUX_0%EFF_PPT     + DOT_PRODUCT(FDFLUX(:)%EFF_PPT,    DELS(:))
  CASE('EVAP_1A')    ; M_FLUX%EVAP_1A     = FLUX_0%EVAP_1A     + DOT_PRODUCT(FDFLUX(:)%EVAP_1A,    DELS(:))
  CASE('EVAP_1B')    ; M_FLUX%EVAP_1B     = FLUX_0%EVAP_1B     + DOT_PRODUCT(FDFLUX(:)%EVAP_1B,    DELS(:))
  CASE('EVAP_1')     ; M_FLUX%EVAP_1      = FLUX_0%EVAP_1      + DOT_PRODUCT(FDFLUX(:)%EVAP_1,     DELS(:))
  CASE('EVAP_2')     ; M_FLUX%EVAP_2      = FLUX_0%EVAP_2      + DOT_PRODUCT(FDFLUX(:)%EVAP_2,     DELS(:))
  CASE('RCHR2EXCS')  ; M_FLUX%RCHR2EXCS   = FLUX_0%RCHR2EXCS   + DOT_PRODUCT(FDFLUX(:)%RCHR2EXCS,  DELS(:))
  CASE('TENS2FREE_1'); M_FLUX%TENS2FREE_1 = FLUX_0%TENS2FREE_1 + DOT_PRODUCT(FDFLUX(:)%TENS2FREE_1,DELS(:))
  CASE('TENS2FREE_2'); M_FLUX%TENS2FREE_2 = FLUX_0%TENS2FREE_2 + DOT_PRODUCT(FDFLUX(:)%TENS2FREE_2,DELS(:))
  CASE('QSURF')      ; M_FLUX%QSURF       = FLUX_0%QSURF       + DOT_PRODUCT(FDFLUX(:)%QSURF,      DELS(:))
  CASE('QPERC_12')   ; M_FLUX%QPERC_12    = FLUX_0%QPERC_12    + DOT_PRODUCT(FDFLUX(:)%QPERC_12,   DELS(:))
  CASE('QINTF_1')    ; M_FLUX%QINTF_1     = FLUX_0%QINTF_1     + DOT_PRODUCT(FDFLUX(:)%QINTF_1,    DELS(:))
  CASE('QBASE_2')    ; M_FLUX%QBASE_2     = FLUX_0%QBASE_2     + DOT_PRODUCT(FDFLUX(:)%QBASE_2,    DELS(:))
  CASE('QBASE_2A')   ; M_FLUX%QBASE_2A    = FLUX_0%QBASE_2A    + DOT_PRODUCT(FDFLUX(:)%QBASE_2A,   DELS(:))
  CASE('QBASE_2B')   ; M_FLUX%QBASE_2B    = FLUX_0%QBASE_2B    + DOT_PRODUCT(FDFLUX(:)%QBASE_2B,   DELS(:))
  CASE('OFLOW_1')    ; M_FLUX%OFLOW_1     = FLUX_0%OFLOW_1     + DOT_PRODUCT(FDFLUX(:)%OFLOW_1,    DELS(:))
  CASE('OFLOW_2')    ; M_FLUX%OFLOW_2     = FLUX_0%OFLOW_2     + DOT_PRODUCT(FDFLUX(:)%OFLOW_2,    DELS(:))
  CASE('OFLOW_2A')   ; M_FLUX%OFLOW_2A    = FLUX_0%OFLOW_2A    + DOT_PRODUCT(FDFLUX(:)%OFLOW_2A,   DELS(:))
  CASE('OFLOW_2B')   ; M_FLUX%OFLOW_2B    = FLUX_0%OFLOW_2B    + DOT_PRODUCT(FDFLUX(:)%OFLOW_2B,   DELS(:))
  CASE DEFAULT       ; CALL NRERROR('disaggflux: cannot find desired flux')
 END SELECT
 ! --------------------------------------------------------------------------------------
 ! (2) ENSURE THAT THE FLUXES ARE REALISTIC
 ! --------------------------------------------------------------------------------------
 SELECT CASE(TRIM(C_FLUX(IFLUX)%FNAME))
  CASE('EFF_PPT')    ; IF(M_FLUX%EFF_PPT     .LT.ZERO) THEN; M_FLUX%EFF_PPT     = ZERO; EFLAG=.TRUE.; ENDIF
  CASE('EVAP_1A')    ; IF(M_FLUX%EVAP_1A     .LT.ZERO) THEN; M_FLUX%EVAP_1A     = ZERO; EFLAG=.TRUE.; ENDIF
  CASE('EVAP_1B')    ; IF(M_FLUX%EVAP_1B     .LT.ZERO) THEN; M_FLUX%EVAP_1B     = ZERO; EFLAG=.TRUE.; ENDIF
  CASE('EVAP_1')     ; IF(M_FLUX%EVAP_1      .LT.ZERO) THEN; M_FLUX%EVAP_1      = ZERO; EFLAG=.TRUE.; ENDIF
  CASE('EVAP_2')     ; IF(M_FLUX%EVAP_2      .LT.ZERO) THEN; M_FLUX%EVAP_2      = ZERO; EFLAG=.TRUE.; ENDIF
  CASE('RCHR2EXCS')  ; IF(M_FLUX%RCHR2EXCS   .LT.ZERO) THEN; M_FLUX%RCHR2EXCS   = ZERO; EFLAG=.TRUE.; ENDIF
  CASE('TENS2FREE_1'); IF(M_FLUX%TENS2FREE_1 .LT.ZERO) THEN; M_FLUX%TENS2FREE_1 = ZERO; EFLAG=.TRUE.; ENDIF
  CASE('TENS2FREE_2'); IF(M_FLUX%TENS2FREE_2 .LT.ZERO) THEN; M_FLUX%TENS2FREE_2 = ZERO; EFLAG=.TRUE.; ENDIF
  CASE('QSURF')      ; IF(M_FLUX%QSURF       .LT.ZERO) THEN; M_FLUX%QSURF       = ZERO; EFLAG=.TRUE.; ENDIF
  CASE('QPERC_12')   ; IF(M_FLUX%QPERC_12    .LT.ZERO) THEN; M_FLUX%QPERC_12    = ZERO; EFLAG=.TRUE.; ENDIF
  CASE('QINTF_1')    ; IF(M_FLUX%QINTF_1     .LT.ZERO) THEN; M_FLUX%QINTF_1     = ZERO; EFLAG=.TRUE.; ENDIF
  CASE('QBASE_2')    ; IF(M_FLUX%QBASE_2     .LT.ZERO) THEN; M_FLUX%QBASE_2     = ZERO; EFLAG=.TRUE.; ENDIF
  CASE('QBASE_2A')   ; IF(M_FLUX%QBASE_2A    .LT.ZERO) THEN; M_FLUX%QBASE_2A    = ZERO; EFLAG=.TRUE.; ENDIF
  CASE('QBASE_2B')   ; IF(M_FLUX%QBASE_2B    .LT.ZERO) THEN; M_FLUX%QBASE_2B    = ZERO; EFLAG=.TRUE.; ENDIF
  CASE('OFLOW_1')    ; IF(M_FLUX%OFLOW_1     .LT.ZERO) THEN; M_FLUX%OFLOW_1     = ZERO; EFLAG=.TRUE.; ENDIF
  CASE('OFLOW_2')    ; IF(M_FLUX%OFLOW_2     .LT.ZERO) THEN; M_FLUX%OFLOW_2     = ZERO; EFLAG=.TRUE.; ENDIF
  CASE('OFLOW_2A')   ; IF(M_FLUX%OFLOW_2A    .LT.ZERO) THEN; M_FLUX%OFLOW_2A    = ZERO; EFLAG=.TRUE.; ENDIF
  CASE('OFLOW_2B')   ; IF(M_FLUX%OFLOW_2B    .LT.ZERO) THEN; M_FLUX%OFLOW_2B    = ZERO; EFLAG=.TRUE.; ENDIF
  CASE DEFAULT       ; CALL NRERROR('disaggflux: cannot find desired flux')
 END SELECT
END DO ! (loop through fluxes)
! deal with surface runoff
IF(M_FLUX%QSURF.GT.M_FLUX%EFF_PPT) THEN; M_FLUX%QSURF = M_FLUX%EFF_PPT; EFLAG=.TRUE.; ENDIF
! deal with evaporation
TOTEVAP = M_FLUX%EVAP_1+M_FLUX%EVAP_2
IF (TOTEVAP.GT.MFORCE%PET) THEN
 M_FLUX%EVAP_1 = (M_FLUX%EVAP_1/TOTEVAP) * MFORCE%PET
 M_FLUX%EVAP_2 = (M_FLUX%EVAP_2/TOTEVAP) * MFORCE%PET
 EFLAG=.TRUE.
ENDIF
! ---------------------------------------------------------------------------------------
! (2) ENSURE THAT THE bucket overflow fluxes are less than the bucket INFLUX
! ---------------------------------------------------------------------------------------
DO ISTT=1,NSTATE
 SELECT CASE(CSTATE(ISTT)%iSNAME)
  CASE (iopt_TENS1A); IN_FLUX = M_FLUX%EFF_PPT - M_FLUX%QSURF
                   IF (M_FLUX%RCHR2EXCS  .GT.IN_FLUX) THEN; M_FLUX%RCHR2EXCS  =IN_FLUX; EFLAG=.TRUE.; ENDIF
  CASE (iopt_TENS1B); IN_FLUX = M_FLUX%RCHR2EXCS
                   IF (M_FLUX%TENS2FREE_1.GT.IN_FLUX) THEN; M_FLUX%TENS2FREE_1=IN_FLUX; EFLAG=.TRUE.; ENDIF
  CASE (iopt_TENS_1); IN_FLUX = M_FLUX%EFF_PPT - M_FLUX%QSURF
                   IF (M_FLUX%TENS2FREE_1.GT.IN_FLUX) THEN; M_FLUX%TENS2FREE_1=IN_FLUX; EFLAG=.TRUE.; ENDIF
  CASE (iopt_FREE_1); IN_FLUX = M_FLUX%TENS2FREE_1
                   IF (M_FLUX%OFLOW_1    .GT.IN_FLUX) THEN; M_FLUX%OFLOW_1    =IN_FLUX; EFLAG=.TRUE.; ENDIF
  CASE (iopt_WATR_1); IN_FLUX = M_FLUX%EFF_PPT - M_FLUX%QSURF
                   IF (M_FLUX%OFLOW_1    .GT.IN_FLUX) THEN; M_FLUX%OFLOW_1    =IN_FLUX; EFLAG=.TRUE.; ENDIF
  CASE (iopt_TENS_2); IN_FLUX = M_FLUX%QPERC_12*(1._SP-MPARAM%PERCFRAC)
                   IF (M_FLUX%TENS2FREE_2.GT.IN_FLUX) THEN; M_FLUX%TENS2FREE_2=IN_FLUX; EFLAG=.TRUE.; ENDIF
  CASE (iopt_FREE2A); IN_FLUX = M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP)
                   IF (M_FLUX%OFLOW_2A   .GT.IN_FLUX) THEN; M_FLUX%OFLOW_2A   =IN_FLUX; EFLAG=.TRUE.; ENDIF
  CASE (iopt_FREE2B); IN_FLUX = M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP)
                   IF (M_FLUX%OFLOW_2B   .GT.IN_FLUX) THEN; M_FLUX%OFLOW_2B   =IN_FLUX; EFLAG=.TRUE.; ENDIF
  CASE (iopt_WATR_2); IN_FLUX = M_FLUX%QPERC_12
                   IF (M_FLUX%OFLOW_2    .GT.IN_FLUX) THEN; M_FLUX%OFLOW_2    =IN_FLUX; EFLAG=.TRUE.; ENDIF
 END SELECT
END DO
! ---------------------------------------------------------------------------------------
! compute total overflow from the lower zone
IF (SMODL%iARCH1.EQ.iopt_tension2_1) M_FLUX%OFLOW_2 = M_FLUX%OFLOW_2A + M_FLUX%OFLOW_2B
! ---------------------------------------------------------------------------------------
END SUBROUTINE DISAGGFLUX
END MODULE DISAGGFLUX_MODULE
