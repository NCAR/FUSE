SUBROUTINE GETMAHUDAT(NFORCE)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2008
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Read Mahurangi data from NetCDF files
!
! Data is stored in two files
!  Rain = rain_wra-mixed_1997010100_2002123100_02001770_hourly.nc
!  PET  = pet_wra-mixed_1997010100_2002123100_02001770_hourly.nc
!
! The rain file includes data from 13 stations, and the potential ET file includes PET
! estimates for the lowest elevation and highest elevation sub-basin in the Mahurangi.
!
! Simply average over the spatial dimension.
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multiforce -- populate structure AFORCE(*)%(*)
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! data types, etc.
USE ddirectory                                        ! defines data directory                               
USE multiforce                                        ! model forcing structures
USE multiroute                                        ! model routing structures
IMPLICIT NONE
! internal
INTEGER(I4B)                           :: I           ! looping
integer(i4b),parameter::lenPath=1024 !DK211008: allows longer file paths
INTEGER(I4B)                           :: IBEG,IEND   ! start/end indices of desired data
INTEGER(I4B)                           :: IVAR        ! loop through variables
CHARACTER(LEN=lenPath)                 :: FNAME_INPUT ! name of input file
CHARACTER(LEN=64)                      :: VARNAME     ! name of variable
INTEGER(I4B)                           :: IERR        ! error code
INTEGER(I4B)                           :: NCID        ! NetCDF file ID
INTEGER(I4B)                           :: IDIMID      ! NetCDF dimension ID
INTEGER(I4B)                           :: IVARID      ! NetCDF variable ID
INTEGER(I4B)                           :: NTIM        ! number of data intervals
INTEGER(I4B)                           :: NSTN        ! number of stations
REAL(DP),DIMENSION(:),ALLOCATABLE      :: ATIME       ! time vector
REAL(MSP),DIMENSION(:,:),ALLOCATABLE   :: TDATA       ! space-time data array
REAL(SP)                               :: TAVE        ! average of temporary data for one time interval
CHARACTER(LEN=256)                     :: TUNITS      ! time units
REAL(DP)                               :: REF_ZERO    ! ref date in sec since year dot
REAL(DP)                               :: JULDAYSS    ! FUNCTION NAME, used to compute REF_ZERO
REAL(DP)                               :: JUL_TIME    ! time stamp -- date in sec since year dot
INTEGER(I4B)                           :: ITIM        ! loop through time
INTEGER(I4B)                           :: JTIM        ! time index in output array
INTEGER(I4B)                           :: IY,IM,ID,IH ! reference time
INTEGER(I4B)                           :: JY,JM,JD,JH ! time for a given time step
INTEGER(I4B)                           :: JMIN        ! minute (NOT USED -- returned by caldatss.f)
REAL(DP)                               :: JSEC        ! second (NOT USED -- returned by caldatss.f)
INTEGER(I4B)                           :: ISTA        ! index of station desired
REAL(DP)                               :: AREA_K2     ! catchment area (km^2)
REAL(DP)                               :: AREA_M2     ! catchment area (m^2)
! output
INTEGER(I4B), INTENT(OUT)              :: NFORCE      ! number of time steps
include 'netcdf.inc'                                  ! use netCDF libraries
! ---------------------------------------------------------------------------------------
! define the start and end indices
IBEG =5088 ; IEND=40151; NFORCE=(IEND-IBEG)+1
! allocate space for the forcing structure (shared in module multiforce)
ALLOCATE(AFORCE(NFORCE),STAT=IERR); IF(IERR.NE.0) STOP ' problem allocating space for AFORCE '
! allocate space for the output structure (shared in module multiroute)
ALLOCATE(AROUTE(NFORCE),STAT=IERR); IF(IERR.NE.0) STOP ' problem allocating space for AROUTE '
! define catchment attributes
ISTA    = 1                     ! station #1 is Mahurangi at College
AREA_K2 = 46.650_dp             ! Mahurangi catchment area (km^2)
AREA_M2 = AREA_K2 * 1000000._dp ! Mahurangi catchment area (m^2)
! loop through variables (1=rain, 2=pet, 3=flow)
DO IVAR=1,3
 ! define variable names
 FORALL(I=1:LEN(VARNAME)) VARNAME(I:I) = ' '
 IF (IVAR.EQ.1) VARNAME='rain'
 IF (IVAR.EQ.2) VARNAME='pet'
 IF (IVAR.EQ.3) VARNAME='flow'
 ! ---------------------------------------------------------------------------------------
 ! (1) EXTRACT DATA FROM NETCDF FILES
 ! ---------------------------------------------------------------------------------------
 ! define filenames
 FORALL(I=1:LEN(FNAME_INPUT)) FNAME_INPUT(I:I) = ' '
 FNAME_INPUT = DATA_PATH(1:LEN_TRIM(DATA_PATH))//TRIM(VARNAME)//&
               '_wra-mixed_1997010100_2002123100_02001770_hourly.nc'
 ! open file
 IERR = NF_OPEN(TRIM(FNAME_INPUT),NF_NOWRITE,NCID); CALL HANDLE_ERR(IERR)
  ! get the number of time elements
  IERR = NF_INQ_DIMID(NCID,'time',IDIMID); CALL HANDLE_ERR(IERR)
  IERR = NF_INQ_DIMLEN(NCID,IDIMID,NTIM);  CALL HANDLE_ERR(IERR)
  ! get the number of "stations"
  IERR = NF_INQ_DIMID(NCID,'station',IDIMID); CALL HANDLE_ERR(IERR)
  IERR = NF_INQ_DIMLEN(NCID,IDIMID,NSTN);     CALL HANDLE_ERR(IERR)
  ! allocate space for temporary arrays
  ALLOCATE(ATIME(NTIM),STAT=IERR); IF(IERR.NE.0) STOP ' problem allocating space for ATIME '
  ALLOCATE(TDATA(NSTN,NTIM),STAT=IERR); IF(IERR.NE.0) STOP ' problem allocating space for TDATA '
  ! get the time data
  IERR = NF_INQ_VARID(NCID,'time',IVARID); CALL HANDLE_ERR(IERR)
  IERR = NF_GET_VARA_DOUBLE(NCID,IVARID,(/1/),(/NTIM/),ATIME); CALL HANDLE_ERR(IERR)
  IERR = NF_GET_ATT_TEXT(NCID,IVARID,'units',TUNITS); CALL HANDLE_ERR(IERR)
  ! get the data
  IERR = NF_INQ_VARID(NCID,TRIM(VARNAME),IVARID); CALL HANDLE_ERR(IERR)
  IERR = NF_GET_VARA_REAL(NCID,IVARID,(/1,1/),(/NSTN,NTIM/),TDATA); CALL HANDLE_ERR(IERR)
 ! close the NetCDF file
 IERR = NF_CLOSE(NCID); CALL HANDLE_ERR(IERR)
 ! ---------------------------------------------------------------------------------------
 ! (2) PUT DATA INTO DATA STRUCTURES
 ! ---------------------------------------------------------------------------------------
 ! convert the ref date in units of seconds since year dot
 CALL EXTRACTOR(TUNITS,IY,IM,ID,IH)  ! get year, month, day, hour, of reference date
 REF_ZERO = JULDAYSS(IY,IM,ID,IH)    ! get the ref date in units of seconds since year dot
 ! loop through time
 DO ITIM=MAX(1,IBEG),MIN(IEND,NTIM)
  ! define time index in output array
  JTIM = (ITIM-IBEG)+1
  ! put time in time arrays
  JUL_TIME = REF_ZERO+ATIME(ITIM)                ! get the julian time (double precision)
  IF (IVAR.EQ.1) THEN
   ! get the year/month/day/hour/minute/second (+0.1 sec to avoid min=59 sec=60)
   CALL CALDATSS(JUL_TIME+0.1_sp,JY,JM,JD,JH,JMIN,JSEC); JSEC = ANINT(JSEC)                                
   AFORCE(JTIM)%IY = JY; AFORCE(JTIM)%IM = JM; AFORCE(JTIM)%ID = JD; AFORCE(JTIM)%IH = JH
   AFORCE(JTIM)%IMIN = JMIN; AFORCE(JTIM)%DSEC = JSEC; AFORCE(JTIM)%DTIME = JUL_TIME
  ! check that the time matches
  ELSE
   IF (ABS(AFORCE(JTIM)%DTIME - JUL_TIME) .GT. 1.0D0) THEN  ! (one-second tolerance)
    WRITE(*,'(2(F20.1,1X))') AFORCE(JTIM)%DTIME, JUL_TIME
    STOP ' mis-match in time '
   ENDIF
  ENDIF 
  ! compute average from temporary data array (and convert mm/h --> mm/d)
  IF (TRIM(VARNAME).EQ.'rain' .OR. TRIM(VARNAME).EQ.'pet') THEN
   TAVE = (SUM(TDATA(:,ITIM))/NSTN)*24.            ! compute average
   IF (ANY(TDATA(:,ITIM) .LT. 0.)) STOP ' MISSING FORCING DATA IN DESIRED TIME RANGE '
  ENDIF
  ! select a station (and convert from m3/s to mm/h)
  IF (TRIM(VARNAME).EQ.'flow') THEN
   TAVE = (TDATA(ISTA,ITIM)/AREA_M2)*1000.*3600.   ! m3/s --> mm/h
   IF (TDATA(ISTA,ITIM) .LT. 0.) STOP ' MISSING VALIDATION DATA IN DESIRED TIME RANGE '
  ENDIF
  ! put data in the data structures
  IF (TRIM(VARNAME).EQ.'rain') AFORCE(JTIM)%PPT  = TAVE
  IF (TRIM(VARNAME).EQ.'pet')  AFORCE(JTIM)%PET  = TAVE
  IF (TRIM(VARNAME).EQ.'flow') AFORCE(JTIM)%OBSQ = TAVE
  !IF (IVAR.EQ.3) &
  ! WRITE(*,'(I10,1X,I4,1X,4(I2,1X),F9.3,1X,F15.1,1X,3(ES12.4,1X))') ITIM, AFORCE(JTIM)
 END DO  ! (looping through time)
 ! deallocate arrays
 DEALLOCATE(ATIME,TDATA, STAT=IERR); IF (IERR.NE.0) STOP ' problem deallocating ATIME/TDATA '
END DO  ! (looping through variables)
! flush buffer
CALL FLUSH(6)
! save the number of time steps
NUMTIM = NFORCE    ! (NUMTIM is stored in module multiforce)
! save the time step (DELTIM is stored in module multiforce) 
DELTIM = (AFORCE(2)%DTIME - AFORCE(1)%DTIME) / 86400._sp
!pause
! ---------------------------------------------------------------------------------------
END SUBROUTINE GETMAHUDAT
