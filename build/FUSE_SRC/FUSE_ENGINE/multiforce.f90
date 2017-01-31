! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark
! Modified by Brian Henn to include snow model, 6/2013
! Modified by Nans Addor to deal with missing values, 8/2016
! Modified by Nans Addor to enable distributed modeling, 9/2016
! ---------------------------------------------------------------------------------------
MODULE multiforce
 USE nrtype
 SAVE
 ! --------------------------------------------------------------------------------------
 ! the time data structure (will have no spatial dimension)
 TYPE TDATA
  INTEGER(I4B)                         :: IY         ! year
  INTEGER(I4B)                         :: IM         ! month
  INTEGER(I4B)                         :: ID         ! day
  INTEGER(I4B)                         :: IH         ! hour
  INTEGER(I4B)                         :: IMIN       ! minute
  REAL(SP)                             :: DSEC       ! second
  REAL(SP)                             :: DTIME      ! time in seconds since year dot
 ENDTYPE TDATA
 ! the response structure (will not have a spatial dimension)
 TYPE VDATA
  REAL(SP)                             :: OBSQ       ! observed runoff (mm day-1)
 END TYPE VDATA
 ! ancillary forcing variables used to compute ET (will have a spatial dimension)
 TYPE ADATA
  REAL(SP)                             :: AIRTEMP    ! air temperature (K)
  REAL(SP)                             :: SPECHUM    ! specific humidity (g/g)
  REAL(SP)                             :: AIRPRES    ! air pressure (Pa)
  REAL(SP)                             :: SWDOWN     ! downward sw radiation (W m-2)
  REAL(SP)                             :: NETRAD     ! net radiation (W m-2)
 END TYPE ADATA
 ! the forcing data structure (will have a spatial dimension)
 TYPE FDATA
  REAL(SP)                             :: PPT        ! water input: rain + melt (mm day-1)
  REAL(SP)                             :: TEMP       ! temperature for snow model (deg.C)
  REAL(SP)                             :: PET        ! energy input: potential ET (mm day-1)
 ENDTYPE FDATA
 ! --------------------------------------------------------------------------------------
 ! general
 INTEGER(I4B),PARAMETER                :: STRLEN=256 ! length of the character string
 ! time data structures
 type(tData)                           :: timDat     ! model time structure
 ! response data structures
 type(vData)                           :: valDat     ! validation structure
 type(vData), dimension(:), pointer    :: aValid     ! all model validation data
 ! forcing data structures
 TYPE(FDATA), DIMENSION(:), POINTER    :: CFORCE     ! COPY of model forcing data
 TYPE(FDATA), DIMENSION(:), POINTER    :: AFORCE     ! all model forcing data
 TYPE(FDATA)                           :: MFORCE     ! model forcing data for a single time step
 type(fData), dimension(:,:), pointer  :: gForce     ! model forcing data for a 2-d grid
 type(aData), dimension(:,:), pointer  :: ancilF     ! ancillary forcing data for the 2-d grid
 ! timing information
 real(sp)                              :: jdayRef                   ! reference time (days)
 real(sp)                              :: deltim=-1._dp             ! length of time step (days)
 integer(i4b)                          :: warmup_beg=-1             ! index for the start of the warm-up period
 integer(i4b)                          :: infern_beg=-1             ! index for the start of the inference period
 integer(i4b)                          :: infern_end=-1             ! index for the end of the inference period
 integer(i4b)                          :: istart=-1                 ! index for start of inference period (in reduced array)
 integer(i4b)                          :: numtim=-1                 ! number of time steps (in reduced array)
 character(len=strLen)                 :: timeUnits                 ! units string for time
 ! lat-lon
 real(sp)                              :: xlon                      ! longitude (degrees)
 real(sp)                              :: ylat                      ! latitude (degrees)
 ! dimension information
 integer(i4b)                          :: nSpat1=-1                 ! number of points in 1st spatial dimension
 integer(i4b)                          :: nSpat2=-1                 ! number of points in 2nd spatial dimension
 integer(i4b)                          :: nsteps=-1                 ! number of data steps
 ! filename
 character(len=StrLen)                 :: forcefile='undefined'     ! name of forcing file
 ! NetCDF
 integer(i4b)                          :: ncid_forc=-1              ! NetCDF forcing file ID
 ! name of time variables
 character(len=StrLen)                 :: vname_iy   ='undefined'   ! name of variable for year
 character(len=StrLen)                 :: vname_im   ='undefined'   ! name of variable for month
 character(len=StrLen)                 :: vname_id   ='undefined'   ! name of variable for day
 character(len=StrLen)                 :: vname_ih   ='undefined'   ! name of variable for hour
 character(len=StrLen)                 :: vname_imin ='undefined'   ! name of variable for minute
 character(len=StrLen)                 :: vname_dsec ='undefined'   ! name of variable for second
 character(len=StrLen)                 :: vname_dtime='undefined'   ! name of variable for time
 ! variable names
 character(len=StrLen)                 :: vname_aprecip='undefined' ! variable name: precipitation
 character(len=StrLen)                 :: vname_airtemp='undefined' ! variable name: temperature
 character(len=StrLen)                 :: vname_spechum='undefined' ! variable name: specific humidity
 character(len=StrLen)                 :: vname_airpres='undefined' ! variable name: surface pressure
 character(len=StrLen)                 :: vname_swdown ='undefined' ! variable name: downward shortwave radiation
 character(len=StrLen)                 :: vname_potevap='undefined' ! variable name: potential ET
 character(len=StrLen)                 :: vname_q      ='undefined' ! variable name: runoff
 ! indices for time data (only used in ASCII files)
 integer(i4b)                          :: ivarid_iy=-1              ! variable ID for year
 integer(i4b)                          :: ivarid_im=-1              ! variable ID for month
 integer(i4b)                          :: ivarid_id=-1              ! variable ID for day
 integer(i4b)                          :: ivarid_ih=-1              ! variable ID for hour
 integer(i4b)                          :: ivarid_imin=-1            ! variable ID for minute
 integer(i4b)                          :: ivarid_dsec=-1            ! variable ID for second
 ! indices for variables
 integer(i4b)                          :: ivarid_ppt=-1             ! variable ID for precipitation
 integer(i4b)                          :: ivarid_temp=-1            ! variable ID for temperature
 integer(i4b)                          :: ivarid_pet=-1             ! variable ID for potential ET
 integer(i4b)                          :: ivarid_q=-1               ! variable ID for runoff
 ! multipliers for variables to convert fluxes to mm/day
 real(sp)                              :: amult_ppt=-1._dp          ! convert precipitation to mm/day
 real(sp)                              :: amult_pet=-1._dp          ! convert potential ET to mm/day
 real(sp)                              :: amult_q=-1._dp            ! convert runoff to mm/day
 ! missing values
 INTEGER(I4B)                          :: NA_VALUE  		    ! integer designating missing values

 ! --------------------------------------------------------------------------------------
END MODULE multiforce
