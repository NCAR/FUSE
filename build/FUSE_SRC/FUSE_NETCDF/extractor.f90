!----------------------------------------------------------------------------------------
!   This is part of the code used as a replacement for the udunits
!   libraries.  It extracts the year, month and day from the reference
!   given in the netCDF data file
!
!   David Rupp -- 2006-March-07
!                    
!----------------------------------------------------------------------------------------
SUBROUTINE EXTRACTOR(REFDATE,YY,IM,DD,HH)
USE nrtype    
IMPLICIT NONE
CHARACTER(LEN=50)                    :: REFDATE  ! ref time and units string netCDF file
CHARACTER(LEN=50)                    :: REFD     ! temporary time and units string
CHARACTER(LEN=4)                     :: CYYYY    ! char year extracted from UNITSTR
CHARACTER(LEN=2)                     :: CMM, CDD, CHH ! char month and day and hour extracted from UNISTR
INTEGER(I4B)                         :: POSIT    ! used to extract date from UNITSTR
INTEGER(I4B)                         :: YY,IM,DD,HH   ! start time (year,month,day,hour)

! strip out time units, if they exist (seconds since , days since , hours since )
REFD = TRIM(REFDATE)
POSIT = INDEX(REFDATE, 'since')
IF (POSIT.GT.0) REFD = REFD(POSIT+6:50)  ! +6 because 'since' has 5 characters
! get the year
POSIT = INDEX(REFD, '-')     ! up to -
CYYYY = REFD(1:POSIT-1)
! get the month
REFD = REFD(POSIT+1:50)
POSIT = INDEX(REFD, '-')     ! up to -
CMM = REFD(1:POSIT-1)
! get the day
REFD = REFD(POSIT+1:50)      
POSIT = INDEX(REFD, ' ')     ! up to space
CDD = REFD(1:POSIT-1)
! get the hour
REFD = REFD(POSIT+1:50)
POSIT = INDEX(REFD, ':')    ! up to :
IF (POSIT.GT.0) THEN
 CHH = REFD(1:POSIT-1)
ELSE
 CHH = '00'
ENDIF
! convert to integers
READ(CYYYY,'(i4)') YY  
READ(CMM,'(i2)') IM
READ(CDD,'(i2)') DD  
READ(CHH,'(i2)') HH
END SUBROUTINE EXTRACTOR
