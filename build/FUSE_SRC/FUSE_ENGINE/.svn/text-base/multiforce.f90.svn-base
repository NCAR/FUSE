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
  REAL(SP)                             :: PET        ! energy input: potential ET (mm day-1)
  REAL(SP)                             :: OBSQ       ! observed runoff (mm day-1)
 ENDTYPE FDATA
 ! --------------------------------------------------------------------------------------
 TYPE(FDATA), DIMENSION(:), POINTER    :: CFORCE     ! COPY of model forcing data
 TYPE(FDATA), DIMENSION(:), POINTER    :: AFORCE     ! all model forcing data
 TYPE(FDATA)                           :: MFORCE     ! model forcing data for a single time step 
 INTEGER(I4B)                          :: ISTART     ! index for start of the inference period
 INTEGER(I4B)                          :: NUMTIM     ! number of time steps
 REAL(SP)                              :: DELTIM     ! length of time step (days)
 ! --------------------------------------------------------------------------------------
END MODULE multiforce
