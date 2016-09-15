SUBROUTINE CLRSKY_RAD(MONTH,DAY,HOUR,DT,SLOPE,AZI,LAT,HRI,COSZEN)
! ----------------------------------------------------------------------------------------
! Used to get hourly radiation index
!
! Modification history
!  - comments added by David Rupp 2006.
!
! Procedure appears similar to Stull(1988) "An Introduction to Boundary Layer
! Meteorology " as seen in matlab routine obtained from Joe Kidston <joekidston@yahoo.co.uk>. 
! Note that equation of time is not used.  Also, solar declination is assumed to stay
! constant over the period of one day. What other assumptions are made?  Is this
! adequate for the level of precision we require?  Worth reading Stull(1988).
!
!  - Modified to integrate over time step up to, but not greater than, 24 hours (D. Rupp, July 2006)
!  
! ----------------------------------------------------------------------------------------
USE nrtype
IMPLICIT NONE
! Input variables
INTEGER(I4B), INTENT(IN)                  :: MONTH   ! month as mm integer
INTEGER(I4B), INTENT(IN)                  :: DAY     ! day of month as dd integer
REAL(SP), INTENT(IN)                      :: HOUR    ! hour of day as real
REAL(SP), INTENT(IN)                      :: DT      ! time step in units of hours
REAL(SP), INTENT(IN)                      :: SLOPE   ! slope of ground surface in degrees
REAL(SP), INTENT(IN)                      :: AZI     ! aspect (azimuth) of ground surface in degrees
REAL(SP), INTENT(IN)                      :: LAT     ! latitude in degrees (negative for southern hemisphere)
! Outputs
REAL(SP), INTENT(OUT)                     :: HRI     ! average radiation index over time step DT
REAL(SP), INTENT(OUT)                     :: COSZEN  ! average cosine of the zenith angle over time step DT
! Internal
REAL(SP)                                  :: CRAD    ! conversion from degrees to radians
REAL(SP)                                  :: YRAD    ! conversion from year to radians
REAL(SP)                                  :: T       ! time from noon in radians
REAL(SP)                                  :: DELT1   ! time step in radians
REAL(SP)                                  :: SLOPE1  ! slope of ground surface in radians
REAL(SP)                                  :: AZI1    ! aspect (azimuth) of ground surface in radians
REAL(SP)                                  :: LAT1    ! latitude in radians
REAL(SP)                                  :: FJULIAN ! julian date as real
REAL(SP)                                  :: D       ! solar declination
REAL(SP)                                  :: LP      ! latitude adjusted for non-level surface (= LAT1 for level surface)
REAL(SP)                                  :: TD      ! used to calculate sunrise/set
REAL(SP)                                  :: TPI     ! used to calculate sunrise/set
REAL(SP)                                  :: TP      ! used to calculate sunrise/set
REAL(SP)                                  :: DDT     ! used to calculate sunrise/set(= 0 for level surface)
REAL(SP)                                  :: T1      ! first time in time step or sunrise
REAL(SP)                                  :: T2      ! last time in time step or sunset
! ----------------------------------------------------------------------------------------
! CONVERSION FACTORS
!   degrees to radians
CRAD=PI/180._sp
!   days-of-year to radians
YRAD=2._sp*PI/365._sp
! CONVERT TIME TO RADIANS FROM NOON
T=(HOUR-12._sp)*PI/12._sp
! Convert time step to radians
DELT1=DT*PI/12._sp
! CONVERT ground slope, ground aspect, and latitude TO RADIANS
SLOPE1=SLOPE*CRAD  ! tilt angle
AZI1=AZI*CRAD ! surface-solar Azimuth ??
LAT1=LAT*CRAD ! latitude
! Calculate julian date
FJULIAN=real(JULIAN(MONTH,DAY), kind(sp))
! Calculate solar declination
D=CRAD*23.5_sp*SIN((FJULIAN-82._sp)*YRAD)
! Calculate latitude "adjustment" for ground slope, aspect and latitude (LP = LAT1 for level surface)
LP=ASIN(SIN(SLOPE1)*COS(AZI1)*COS(LAT1) + COS(SLOPE1)*SIN(LAT1)) ! angle between solar rays and surface (tilted) ??
! Calculate time of sunrise/sunset on level surface as radians from noon
TD=ACOS(-TAN(LAT1)*TAN(D))
! print *, 'Sunrise = ', TD
! Calculate time of sunrise/sunset adjusted for inclined ground surface as radians from noon???
TPI=-TAN(LP)*TAN(D)
IF(ABS(TPI).LT.1._sp) THEN
 TP=ACOS(TPI)
ELSE
 TP=0._sp
ENDIF
! Calculate time adjustment for ground slope, aspect and latitude (DDT = 0 for level surface)
DDT=ATAN(SIN(AZI1)*SIN(SLOPE1)/(COS(SLOPE1)*COS(LAT1)-COS(AZI1)*SIN(SLOPE1)*SIN(LAT1)))
! Set beginning time of time step (set to sunrise if before sunrise)
T1=MAX(T,-TP-DDT,-TD)
! Set end time of time step (adjust if after sunset)
T2=MIN(T+DELT1,TD,TP-DDT)
! print *, 'First t1 and t2 = ', t1, t2
IF(T2.LE.T1) THEN
 HRI=0._sp ! nighttime
ELSE
! Calculate integral of radiation index from T1 to T2 and divide by time step DELTA1
! NOTE: this assumes the declination does not change from T1 to T2
 HRI=(SIN(D)*SIN(LP)*(T2-T1)+COS(D)*COS(LP)*(SIN(T2+DDT) &
     -SIN(T1+DDT)))/(COS(SLOPE1)*DELT1)  ! radiation index 
ENDIF
! print *, hri
! ----------------- for time intervals that extend to following day ----------------------
! Check to see of timestep extends to following day
IF((T+DELT1).GT.PI) THEN
 ! Advance julian day by 1
 FJULIAN = FJULIAN + 1._sp
 ! Calculate solar declination
 D=CRAD*23.5_sp*SIN((FJULIAN-82._sp)*YRAD)
 ! Calculate time of sunrise/sunset on level surface as radians from noon
 TD=ACOS(-TAN(LAT1)*TAN(D))
 ! print *, 'Sunrise #2 = ', TD, DELT1
 ! Calculate time of sunrise/sunset adjusted for inclined ground surface as radians from noon???
 TPI=-TAN(LP)*TAN(D)
 IF(ABS(TPI).LT.1._sp) THEN
  TP=ACOS(TPI)
 ELSE
  TP=0._sp
 ENDIF
 ! Set beginning time to sunrise
 T1=MAX(-TP-DDT,-TD)
 ! Set end time of time step
 T2=MIN(T+DELT1-2*PI,TD,TP-DDT)
 ! print *, 'Second t1 and t2 = ', t1, t2
 IF(T2.LE.T1) THEN
  HRI=HRI ! still nighttime in day 2
 ELSE
  ! Calculate integral of radiation index from T1 to T2 and divide by time step DELTA1
  ! NOTE: this assumes the declination does not change from T1 to T2
  HRI=HRI+(SIN(D)*SIN(LP)*(T2-T1)+COS(D)*COS(LP)*(SIN(T2+DDT) &
          -SIN(T1+DDT)))/(COS(SLOPE1)*DELT1)  ! radiation index 
 ENDIF
 ! print *, hri
ENDIF
! ----------------------------------------------------------------------------------------
! Calculate cosine of solar zenith angle (= HRI for level surface)
COSZEN = HRI*COS(SLOPE1)
! this is assumed to be an appropriate representative value over the
! time step.  It is used for albedo calculations.
! ----------------------------------------------------------------------------------------
CONTAINS
 FUNCTION JULIAN(MONTH,DAY)
 USE nrtype
 IMPLICIT NONE
 ! input
 INTEGER(I4B)                             :: MONTH,DAY  ! month and day    
 ! output
 INTEGER(I4B)                             :: JULIAN     ! julian day
 ! internal
 INTEGER(I4B),DIMENSION(12)               :: MADD       ! julian day at start of each month
 ! specify the julian day at the start of each month (-1)
 MADD = (/0,31,59,90,120,151,181,212,243,273,304,334/)
 ! compute the julian day
 JULIAN=DAY+MADD(MONTH)
 RETURN
 END FUNCTION JULIAN
! ----------------------------------------------------------------------------------------
END SUBROUTINE CLRSKY_RAD
