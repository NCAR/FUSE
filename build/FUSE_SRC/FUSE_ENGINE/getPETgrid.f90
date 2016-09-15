module getPETgrid_module
USE nrtype
implicit none
private
public::getPETgrid
contains

 SUBROUTINE getPETgrid(ierr,message)
 ! ---------------------------------------------------------------------------------------
 ! Creator:
 ! --------
 ! Martyn Clark, 2012
 ! ---------------------------------------------------------------------------------------
 ! Purpose:
 ! --------
 ! Compute potential ET on the grid
 ! ---------------------------------------------------------------------------------------
 ! Modules Modified:
 ! -----------------
 ! multiforce: gForce (gridded forcing data)
 ! ---------------------------------------------------------------------------------------
 USE multiforce,only:timDat             ! time
 USE multiforce,only:deltim             ! data interval (days)
 USE multiforce,only:xlon,ylat          ! longitude, latitude (degrees)
 USE multiforce,only:nspat1,nspat2      ! dimension lengths
 USE multiforce,only:gForce             ! gridded forcing data
 USE multiforce,only:ancilF             ! ancillary forcing data
 USE conv_funcs_module,only:sphm2relhm  ! convert specific humidity to relative humidity
 USE conv_funcs_module,only:rlhum2dewpt ! convert relative humidity to dewpoint
 USE conv_funcs_module,only:dewpt2vpair ! convert dewpoint to vapor pressure
 IMPLICIT NONE
 ! output
 integer(i4b), intent(out)              :: ierr                ! error code
 character(*), intent(out)              :: message             ! error message
 ! internal
 integer(i4b),parameter                 :: strLen=1024         ! length of character strings
 character(len=strLen)                  :: cmessage            ! error message of downwind routine
 integer(i4b)                           :: iSpat1,iSpat2       ! indices of spatial dimensions
 real(sp)                               :: csky_rad            ! clear-sky radiation (W m-2)
 real(sp)                               :: relhum              ! relative humidity (-)
 real(sp)                               :: dewpt               ! dewpoint temperature (K)
 real(sp)                               :: vpAir               ! vapor pressure of air (Pa)
 real(sp)                               :: vpSat               ! saturated vapor pressure (Pa)
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 ierr=0; message='getPETgrid/'
 ! ---------------------------------------------------------------------------------------
 ! compute theoretical shortwave radiation (constant over the grid)
 call csky_solar(timDat%im,timDat%id,timDat%ih,deltim,xlon,ylat,& ! input
                 csky_rad,ierr,cmessage)                          ! output
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! (loop through grid points)
 do iSpat1=1,nSpat1
  do iSpat2=1,nSpat2
   ! compute relative humidity (-)
   relhum = sphm2relhm(ancilF(iSpat1,iSpat2)%spechum, ancilF(iSpat1,iSpat2)%airpres, ancilF(iSpat1,iSpat2)%airtemp)
   ! compute dewpoint temperature (K)
   dewpt  = rlhum2dewpt(ancilF(iSpat1,iSpat2)%airtemp, relhum)
   ! compute vapor pressure of air (Pa)
   vpAir  = dewpt2vpair(dewpt)
   ! compute saturated vapor pressure (Pa)
   vpSat  = dewpt2vpair(ancilF(iSpat1,iSpat2)%airtemp) 
   ! calculate net radiation (empirical)
   call calcNetRad(ancilF(iSpat1,iSpat2)%swdown,csky_rad,ancilF(iSpat1,iSpat2)%airtemp,vpAir,&       ! input
                   ancilF(iSpat1,iSpat2)%netRad,ierr,cmessage)   ! output
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   ! calculate potential ET (Priestly-Taylor)
   call computePET(ancilF(iSpat1,iSpat2)%airtemp,ancilF(iSpat1,iSpat2)%airpres,ancilF(iSpat1,iSpat2)%netRad,&  ! input
                   gForce(iSpat1,iSpat2)%pet,ierr,message)    ! output
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end do  ! (looping through 2nd spatial dimension
 end do  ! (looping through 1st spatial dimension)
 end subroutine getPETgrid


 SUBROUTINE csky_solar(im,id,ih,dt,xlon,ylat,& ! input
                       crad,ierr,message)      ! output
 ! ---------------------------------------------------------------------------------------
 ! Creator:
 ! --------
 ! Martyn Clark, 2012
 ! ---------------------------------------------------------------------------------------
 ! Purpose:
 ! --------
 ! Compute theoretical shortwave radiation
 ! ---------------------------------------------------------------------------------------
 ! Modules Modified:
 ! -----------------
 ! None
 ! ---------------------------------------------------------------------------------------
 IMPLICIT NONE
 ! input
 integer(i4b),intent(in)                :: im                  ! month
 integer(i4b),intent(in)                :: id                  ! day
 integer(i4b),intent(in)                :: ih                  ! hour
 real(sp),intent(in)                    :: dt                  ! time interval (seconds)
 real(sp),intent(in)                    :: xlon                ! longitude (degrees)
 real(sp),intent(in)                    :: ylat                ! latitude (degrees)
 ! output
 real(sp),intent(out)                   :: crad                ! clear-sky radiation (W m-2)
 integer(i4b), intent(out)              :: ierr                ! error code
 character(*), intent(out)              :: message             ! error message
 ! internal: physical constants
 real(sp),parameter                     :: solarConst=1365._sp ! solar constant (W m-2)
 ! internal: general
 real(sp),parameter                     :: secprhour=3600._sp  ! number of seconds per hour
 real(sp),parameter                     :: hourprday=24._sp    ! number of hours per day
 integer(i4b),parameter                 :: strLen=1024         ! length of character strings
 character(len=strLen)                  :: cmessage            ! error message of downwind routine
 ! internal: compute time
 real(sp)                               :: tShift              ! time shift from Grenwich (in seconds)
 real(sp)                               :: lHour               ! local hour
 real(sp)                               :: dHour               ! time step (in hours)
 ! internal: compute theoretical clear-sky radiation
 real(sp),parameter                     :: slope=0._sp         ! slope of ground surface (degrees)
 real(sp),parameter                     :: azi=0._dp           ! aspect (azimuth) of ground surface in degrees
 real(sp)                               :: hri                 ! hourly radiation index
 real(sp)                               :: coszen              ! cosize of the solar zenith angle
 ! initialize error control
 ierr=0; message='csky_solar/'
 ! get time shift from Grenwich (in seconds)
 call utc_offset(xlon,tShift,ierr,cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! get local hour
 lHour = real(ih, kind(sp)) !+ tShift/secprhour
 if(lHour < 24._dp) lHour = lHour+24._sp
 if(lHour > 24._dp) lHour = lHour-24._sp
 ! get the time step in hours
 dhour = real(dt, kind(sp))*hourprday
 ! compute theoretical shortwave radiation
 call clrsky_rad(im,id,lHour,dHour,slope,azi,ylat,&  ! input
                 hri,coszen)                                       ! output
 crad = hri*solarConst
 end subroutine csky_solar


 SUBROUTINE calcNetRad(swdown,cskyRad,airtemp,vpAir,&  ! input
                       netRadi,ierr,message)           ! output
 ! ---------------------------------------------------------------------------------------
 ! Creator:
 ! --------
 ! Martyn Clark, 2012
 ! ---------------------------------------------------------------------------------------
 ! Purpose:
 ! --------
 ! Compute net radiation
 ! ---------------------------------------------------------------------------------------
 ! Modules Modified:
 ! -----------------
 ! None
 ! ---------------------------------------------------------------------------------------
 IMPLICIT NONE
 ! input
 real(sp),intent(in)                    :: swdown              ! downward sw radiation (W m-2)
 real(sp),intent(in)                    :: cskyRad             ! clear sky radiation (W m-2)
 real(sp),intent(in)                    :: airtemp             ! air temperature (K)
 real(sp),intent(in)                    :: vpAir               ! vapor pressure (Pa)
 ! output
 real(sp),intent(out)                   :: netRadi             ! net radiation (W m-2)
 integer(i4b), intent(out)              :: ierr                ! error code
 character(*), intent(out)              :: message             ! error message
 ! physical constants
 real(sp),parameter                     :: SBconst   = 5.6705e-8_sp  ! Stefan Boltzman W m-2 K-4
 ! empirical radiation parameters
 real(sp),parameter                     :: sAlbedo=0.3_sp      ! surface albedo
 real(sp),parameter                     :: maxatau=0.8_sp      ! Maximum atmospheric transmissivity
 real(sp),parameter                     :: acloudf=1.0_sp      ! Constant for cloud fraction -humid
 real(sp),parameter                     :: bcloudf=0.0_sp      ! Constant for cloud fraction -humid
 real(sp),parameter                     :: mult_ae=-0.14_sp    ! Multiplier in atmos. emmissivity eqn
 real(sp),parameter                     :: constae=0.34_sp     ! Constant in atmos. emmissivity eqn
 ! internal
 real(sp)                               :: atmTran             ! atmospheric transmissivity (-)
 real(sp)                               :: cldFrac             ! cloud fraction (-)
 real(sp)                               :: netEmis             ! net emissivity (-) 
 real(sp)                               :: netLwRd             ! net longwave radiation (W m-2)
 real(sp),parameter                     :: noRad=1._dp         ! threshold for no radiation
 ! initialize error control
 ierr=0; message='calcNetRad/'
 ! ----------------------------------------------------------------------------------------

 ! calculate atmospheric transmissivity
 if(cskyRad < noRad)then
  atmTran = 1._dp
 else
  atmTran = max(swdown/cskyRad, 1._sp)
 endif
 ! compute the cloud fraction
 cldFrac = acloudf*(atmTran/maxatau) + bcloudf
 ! compute net emissivity (dimensionless) - assumes pressure is in kPa, so divide by 1000
 netEmis = mult_ae*sqrt(vpAir/1000._sp) + constae
 ! compute net longwave radiation (W m-2)
 netLwRd = -cldFrac*netEmis*SBconst*airtemp**4._sp
 ! compute net radiation
 netRadi = (1._sp - sAlbedo)*swdown + netLwRd
 end SUBROUTINE calcNetRad 


 SUBROUTINE computePET(airtemp,airpres,netRadi,&  ! input
                       pet_force,ierr,message)    ! output
 ! ---------------------------------------------------------------------------------------
 ! Creator:
 ! --------
 ! Martyn Clark, 2012
 ! ---------------------------------------------------------------------------------------
 ! Purpose:
 ! --------
 ! Compute potential ET
 ! Note: This *DOES NOT* use the incoming longwave radiation. It is a completely different
 ! empirical computation, based on Shuttleworth (1993) 
 ! ---------------------------------------------------------------------------------------
 ! Modules Modified:
 ! -----------------
 ! None
 ! ---------------------------------------------------------------------------------------
 USE conv_funcs_module,only:dewpt2vpair ! convert dewpoint to vapor pressure
 IMPLICIT NONE
 ! input
 real(sp),intent(in)                    :: airtemp            ! air temperature (K)
 real(sp),intent(in)                    :: airpres            ! air pressure (Pa)
 real(sp),intent(in)                    :: netRadi            ! Net Radiation (W m-2)
 ! output
 real(sp),intent(out)                   :: pet_force          ! potential ET (mm/day)
 integer(i4b), intent(out)              :: ierr               ! error code
 character(*), intent(out)              :: message            ! error message
 ! internal: general parameters
 real(sp),parameter                     :: secprday=86400._sp ! number of seconds per day
 ! internal: physical constants
 real(sp),parameter                     :: TFREEZE   = 273.16_sp ! freezing point of pure water (K)
 ! internal: radiation parameters
 real(sp),parameter                     :: PT_alpha  = 1.26_sp   ! Priestly/Taylor "alpha" multiplier (humid)
 real(sp),parameter                     :: PT_beta   = 0.00_sp   ! Priestly/Taylor "beta" offset
 ! internal
 real(sp)                               :: SAT_VP             ! sat vapor pressure       (Pa)
 real(sp)                               :: TEMP_C             ! temperature              (oC)
 real(sp)                               :: LH_VAP             ! latent heat vaporizn   (J/kg)
 real(sp)                               :: DVP_DT             ! d vap press / d Temp (kPa/oC)
 real(sp)                               :: PSYCON             ! psychometric const   (kPa/oC)
 real(sp)                               :: R_MULT             ! radiation multiplier      (-)
 real(sp)                               :: pet_nrg            ! PET energy units (W m-2, or J m-2 s-1)
 real(sp)                               :: pet_liq            ! PET liquid units (kg m-2 s-1)
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 ierr=0; message='computePET/'
 ! ---------------------------------------------------------------------------------------
 ! compute saturation vapor pressure
 ! NOTE, when air temperature is used instead of dewpoint, returns sat vapour pressure
 SAT_VP = DEWPT2VPAIR(airtemp)                         ! sat vapor pressure      (Pa)
 ! ---------------------------------------------------------------------------------------- 
 ! *** The units in the following lines of code are in kPa, instead of Pa, which is
 !     inconsistent with the rest of the model.  However, the multiplier in the PET
 !     equation is dimensionless.  Still ugly, and needs attention.
 TEMP_C = airtemp-TFREEZE                                          ! air temperature          (oC)
 LH_VAP = (2501._sp - 2.361_sp*TEMP_C) * 1000._sp                  ! latent heat vaporizn   (J/kg)
 DVP_DT = 4098._sp * (SAT_VP/1000._sp) / (237.3_sp+TEMP_C)**2._sp  ! d vap press / d Temp (kPa/oC)
 PSYCON = 1.6286_sp * ((airpres/1000._sp)/(LH_VAP/1000._sp))       ! psychometric const   (kPa/oC)
 R_MULT = DVP_DT/(DVP_DT + PSYCON)                                 ! radiation multiplier      (-)
 ! *** (end of inconstencies and hard-coded constants) ***
 ! ---------------------------------------------------------------------------------------- 
 ! estimate potential ET
 pet_nrg   = PT_alpha*(R_MULT*netRadi) + PT_beta                   ! PET energy units (W m-2, or J m-2 s-1)
 pet_liq   = max(pet_nrg/LH_VAP, 0._sp)                            ! PET liquid units (kg m-2 s-1)
 ! put the PET in the structure
 pet_force = pet_liq*secprday                                      ! convert to mm/day
 end subroutine computePET


 subroutine utc_offset(xlon, tShift, ierr, message)
 ! ---------------------------------------------------------------------------------------
 ! Creator:
 ! --------
 ! Martyn Clark, 2012
 ! ---------------------------------------------------------------------------------------
 ! Purpose:
 ! --------
 ! calculates the offset in seconds from the UTC time based on longitude
 ! ---------------------------------------------------------------------------------------
 ! Modules Modified:
 ! -----------------
 ! None
 ! ---------------------------------------------------------------------------------------
 implicit none
 ! input
 real(sp),intent(in)           :: xlon                ! longitude (degrees)
 ! output
 real(sp),intent(out)          :: tShift              ! time shift (in seconds) from UTC
 integer(i4b), intent(out)     :: ierr                ! error code
 character(*), intent(out)     :: message             ! error message
 ! internal
 real(sp),parameter            :: secprday=86400._sp  ! number of seconds per day
 ! initialize error control
 ierr=0; message='utc_offset/'
 ! check if longitude is in range
 if(xlon > 180._sp .or. xlon < -180._sp)then
  message=trim(message)//'values for longitude only allowed to be in the interval [-180, 180]'
  ierr=20; return
 endif
 ! compute offset in seconds
 tShift = (xlon/360._sp)*secprday
 end subroutine utc_offset


end module getPETgrid_module 
