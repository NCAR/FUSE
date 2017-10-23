SUBROUTINE UPDATE_SWE(DT)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Brian Henn, as part of FUSE snow model implementation, 6/2013
! Based on subroutines QSATEXCESS and UPDATSTATE, by Martyn Clark
! Modified by Nans Addor to enable distributed modeling, 9/2016
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Computes the snow accumulation and melt from forcing data
! Then updates the SWE band states based on the fluxes
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multibands -- SWE bands stored in MODULE multibands
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition structure
USE model_defnames                                    ! integer model definitions
USE multiparam                                        ! model parameters
USE multibands                                        ! model basin band structure
USE multiforce                                        ! model forcing
USE multistate                                        ! model states
USE multi_flux                                        ! model fluxes
IMPLICIT NONE
! input
REAL(SP), INTENT(IN)                   :: DT             ! length of the time step
! internal variables
LOGICAL(LGT)                           :: LEAP           ! leap year flag
REAL(SP)                               :: JDAY           ! Julian day of year
REAL(SP), DIMENSION(12)                :: CUMD           ! cumulative days of year
REAL(SP)                               :: DZ             ! vert. distance from forcing
REAL(SP)                               :: TEMP_Z         ! band temperature at timestep
REAL(SP)                               :: PRECIP_Z       ! band precipitation at timestep
REAL(SP)                               :: MF             ! melt factor (mm/deg.C-6hr)
INTEGER(I4B)                           :: ISNW           ! loop through snow model bands
! ---------------------------------------------------------------------------------------
! snow accumulation and melt calculations for each band
! also calculates effective precipitation
! ---------------------------------------------------------------------------------------
! first calculate day of year for melt factor calculation
CUMD = real((/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /),sp)
IF (MOD(timDat%IY,4).EQ.0) THEN
 LEAP = .TRUE.
 CUMD(3:12) = CUMD(3:12) + 1._sp
ELSE
 LEAP = .FALSE.
ENDIF

JDAY = CUMD(timDat%IM) + timDat%ID

IF (LEAP) THEN   ! calculate melt factor from MFMAX, MFMIN and day of year
 MF = ((0.5_sp*SIN(((JDAY-81._sp)*2._sp*PI)/366._sp))+0.5_sp)*(MPARAM%MFMAX - MPARAM%MFMIN) + MPARAM%MFMIN
ELSE
 MF = ((0.5_sp*SIN(((JDAY-80._sp)*2._sp*PI)/365._sp))+0.5_sp)*(MPARAM%MFMAX - MPARAM%MFMIN) + MPARAM%MFMIN
ENDIF

! loop through model bands
DO ISNW=1,N_BANDS

 ! calculate forcing data for each band
 DZ = MBANDS(ISNW)%Z_MID - Z_FORCING
 TEMP_Z = MFORCE%TEMP + DZ*MPARAM%LAPSE/1000._sp ! adjust for elevation using lapse rate
 IF (DZ.GT.0._sp) THEN ! adjust for elevation using OPG
  PRECIP_Z = MFORCE%PPT * (1._sp + DZ*MPARAM%OPG/1000._sp)
 ELSE
  PRECIP_Z = MFORCE%PPT / (1._sp - DZ*MPARAM%OPG/1000._sp)
 ENDIF
 IF ((MBANDS(ISNW)%SWE.GT.0._sp).AND.(TEMP_Z.GT.MPARAM%MBASE)) THEN
  ! calculate the initial snowmelt rate from the melt factor and the temperature
  MBANDS(ISNW)%SNOWMELT = MF*(TEMP_Z - MPARAM%MBASE) ! MBANDS%SNOWMELT has units of mm day-1
 ELSE
  MBANDS(ISNW)%SNOWMELT = 0.0_sp
 ENDIF

 ! calculate the accumulation rate from the forcing data
 IF (TEMP_Z.LT.MPARAM%PXTEMP) THEN
  SELECT CASE(SMODL%iRFERR)
   CASE(iopt_additive_e) ! additive rainfall error
    MBANDS(ISNW)%SNOWACCMLTN = MAX(0.0_sp, PRECIP_Z + MPARAM%RFERR_ADD)
   CASE(iopt_multiplc_e) ! multiplicative rainfall error
    MBANDS(ISNW)%SNOWACCMLTN = PRECIP_Z * MPARAM%RFERR_MLT
   CASE DEFAULT       ! check for errors
    print *, "SMODL%iRFERR must be either iopt_additive_e or iopt_multiplc_e"
    STOP
  END SELECT
 ELSE
  MBANDS(ISNW)%SNOWACCMLTN = 0.0_sp
 ENDIF

 ! update SWE, and check to ensure non-negative values
 MBANDS(ISNW)%DSWE_DT = MBANDS(ISNW)%SNOWACCMLTN - MBANDS(ISNW)%SNOWMELT
 IF ((MBANDS(ISNW)%SWE + MBANDS(ISNW)%DSWE_DT*DT).GE.0._sp) THEN
  MBANDS(ISNW)%SWE = MBANDS(ISNW)%SWE + MBANDS(ISNW)%DSWE_DT*DT
 ELSE ! reduce melt rate in case of negative SWE
  MBANDS(ISNW)%SNOWMELT = MBANDS(ISNW)%SWE/DT + MBANDS(ISNW)%SNOWACCMLTN
  MBANDS(ISNW)%SWE = 0.0_sp
 ENDIF

 ! calculate rainfall plus snowmelt
 IF (TEMP_Z.GT.MPARAM%PXTEMP) THEN
  SELECT CASE(SMODL%iRFERR)
   CASE(iopt_additive_e) ! additive rainfall error
   M_FLUX%EFF_PPT = M_FLUX%EFF_PPT + MBANDS(ISNW)%AF * &
   (MAX(0.0_sp, PRECIP_Z + MPARAM%RFERR_ADD) + MBANDS(ISNW)%SNOWMELT)
   CASE(iopt_multiplc_e) ! multiplicative rainfall error
   M_FLUX%EFF_PPT = M_FLUX%EFF_PPT + MBANDS(ISNW)%AF * &
   (PRECIP_Z * MPARAM%RFERR_MLT +  MBANDS(ISNW)%SNOWMELT)
   CASE DEFAULT       ! check for errors
    print *, "SMODL%iRFERR must be either iopt_additive_e or iopt_multiplc_e"
    STOP
  END SELECT
 ELSE
  M_FLUX%EFF_PPT = M_FLUX%EFF_PPT + MBANDS(ISNW)%AF * MBANDS(ISNW)%SNOWMELT
 ENDIF
END DO
END SUBROUTINE UPDATE_SWE
