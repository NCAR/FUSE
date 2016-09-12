! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark
! Modified by Brian Henn to include snow model, 6/2013
! Modified by Nans Addor to deal with missing values, 8/2016
! ---------------------------------------------------------------------------------------
MODULE multiforce
 USE nrtype
 TYPE FDATA
  INTEGER(I4B)                         :: IY         ! year
  INTEGER(I4B)                         :: IM         ! month
  INTEGER(I4B)                         :: ID         ! day
  INTEGER(I4B)                         :: IH         ! hour
  INTEGER(I4B)                         :: IMIN       ! minute
  REAL(SP)                             :: DSEC       ! second
  REAL(SP)                             :: DTIME      ! time in seconds since year dot
  REAL(SP)                             :: PPT        ! water input: rain + melt (mm day-1)
  REAL(SP)                             :: TEMP       ! temperature for snow model (deg.C)
  REAL(SP)                             :: PET        ! energy input: potential ET (mm day-1)
  REAL(SP)                             :: OBSQ       ! observed runoff (mm day-1)
 ENDTYPE FDATA
 ! --------------------------------------------------------------------------------------
 TYPE(FDATA), DIMENSION(:), POINTER    :: CFORCE     ! COPY of model forcing data
 TYPE(FDATA), DIMENSION(:), POINTER    :: AFORCE     ! all model forcing data
 TYPE(FDATA)                           :: MFORCE     ! model forcing data for a single time step 
 INTEGER(I4B)                          :: ISTART     ! index for start of the inference period
 INTEGER(I4B)                          :: NUMTIM     ! number of time steps
 INTEGER(I4B)                          :: NA_VALUE   ! integer designating missing values
 REAL(SP)                              :: DELTIM     ! length of time step (days)
 ! --------------------------------------------------------------------------------------
END MODULE multiforce
