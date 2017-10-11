module get_gforce_module
USE nrtype
USE netcdf
implicit none
private
public::read_ginfo
public::get_dimIds
public::get_modtim
public::get_gforce
public::get_gforce_3d
public::get_varid

contains

 SUBROUTINE read_ginfo(ncid,ierr,message)
 ! ---------------------------------------------------------------------------------------
 ! Creator:
 ! --------
 ! Martyn Clark, 2012
 ! ---------------------------------------------------------------------------------------
 ! Purpose:
 ! --------
 ! Read grid info (spatial and temporal dimensions) from the NetCDF file

 ! ---------------------------------------------------------------------------------------
 ! Modules Modified:
 ! -----------------
 ! MODULE multiforce -- populate dimension lengths
 ! ---------------------------------------------------------------------------------------
 USE fuse_fileManager,only:SETNGS_PATH,FORCINGINFO,&   ! defines data directory
                           INPUT_PATH
 USE multiforce,only:forcefile,vname_aprecip           ! model forcing structures
 USE multiforce,only:nspat1,nspat2,numtim_in           ! dimension lengths
 USE multiforce,only:latitude,longitude,time_steps     ! dimension arrays
 USE multiforce,only:latUnits,lonUnits,timeUnits       ! units string for time
 USE multiforce,only:vname_dtime                       ! variable name: time sice reference time

 IMPLICIT NONE
 ! input
 integer(i4b),intent(in)                :: ncid        ! NetCDF file ID
! output
 integer(i4b), intent(out)              :: ierr        ! error code
 character(*), intent(out)              :: message     ! error message
 ! internal: general
 integer(i4b),parameter::lenPath=1024 ! DK211008: allows longer file paths
 INTEGER(I4B)                           :: I           ! looping
 CHARACTER(LEN=lenPath)                 :: cmessage    ! message of downwind routine
 ! internal: NetCDF read
 integer(i4b)                           :: ivarid      ! NetCDF variable ID
 integer(i4b),parameter                 :: ndims=3     ! number of dimensions for precipitation
 integer(i4b),dimension(ndims)          :: dimids_ppt  ! vector of dimension IDs for precipitation
 integer(i4b)                           :: iDimID      ! dimension ID
 integer(i4b)                           :: dimLen      ! dimension length
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 ierr=0; message='read_ginfo/'
 ! ---------------------------------------------------------------------------------------

 ! get the variable ID for precipitation
 ierr = nf90_inq_varid(ncid, vname_aprecip, ivarid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the dimension IDs for precipitation
 call get_dimIds(ncid, ivarid, ndims, dimids_ppt, ierr, cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! loop through dimensions
 do iDimID=1,ndims

  ! get the dimension lengths
  ierr = nf90_inquire_dimension(ncid,dimids_ppt(iDimID),len=dimLen)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! save the dimension lengths
  if(iDimID==1) nspat1 = dimLen  ! 1st spatial dimension
  if(iDimID==2) nspat2 = dimLen  ! 2nd spatial dimension
  if(iDimID==3) numtim_in = dimLen  ! record dimension (always last)

 end do

 allocate(longitude(nspat1),latitude(nspat2),time_steps(numtim_in))

 ! get longitude
 ierr = nf90_inq_varid(ncid, 'longitude', iVarID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'[variable=longitude]'; return; endif
 ierr = nf90_get_var(ncid, iVarID, longitude, start=(/1/), count=(/nSpat1/)); CALL HANDLE_ERR(IERR)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get latitude
 ierr = nf90_inq_varid(ncid, 'latitude', iVarID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'[variable=latitude]'; return; endif
 ierr = nf90_get_var(ncid, iVarID, latitude, start=(/1/), count=(/nSpat2/)); CALL HANDLE_ERR(IERR)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get time
 ierr = nf90_inq_varid(ncid, trim(vname_dtime), iVarID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'[variable='//trim(vname_dtime)//']'; return; endif
 ierr = nf90_get_var(ncid, iVarID, time_steps, start=(/1/), count=(/numtim_in/)); CALL HANDLE_ERR(IERR)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 ierr = nf90_get_att(ncid, iVarID, 'units', timeUnits)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'[variable='//trim(vname_dtime)//']'; return; endif

 end subroutine read_ginfo

  ! --------------------------------------------------------------------------------------
  subroutine get_dimIds(ncid, varid, nexpect, varDimIDs, ierr, message)
  ! used to get the vector of dimension ids for a given variable
  implicit none
  ! input
  integer(i4b),intent(in)   :: ncid     ! NetCDF file ID
  integer(i4b),intent(in)   :: varid    ! NetCDF variable ID
  integer(i4b),intent(in)   :: nexpect  ! number of dimensions expected
  ! output
  integer(i4b),intent(out)  :: varDimIDs(nexpect)  ! vector of dimension IDs
  integer(i4b),intent(out)  :: ierr     ! error code
  character(*), intent(out) :: message  ! error message
  ! internal variables
  integer(i4b)              :: nVarDims ! number of dimensions for given variable
  ! initialize error control
  ierr=0; message='get_dimIds/'
  ! get number of dimensions
  ierr = nf90_inquire_variable(ncid, varid, ndims=nVarDims)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
  ! check number of dimensions
  if(nVarDims/=nexpect)then; message=trim(message)//'unexpected number of dimensions for variable'; return; endif
  ! get vector of dimension IDs
  ierr = nf90_inquire_variable(ncid, varid, dimids=varDimIDs(:nVarDims))
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
  end subroutine get_dimIds
  ! --------------------------------------------------------------------------------------


  SUBROUTINE get_varID(ncid,ierr,message)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Nans Addor, 2017
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Get NetCDF ID for each variable of the forcing file
  ! ---------------------------------------------------------------------------------------
  ! Modules Modified:
  ! -----------------
  ! MODULE multiforce -- populate structure ncid_var%(*)
  ! ---------------------------------------------------------------------------------------
  USE multiforce, only: nForce                           ! number of forcing variables
  USE multiforce, only: ncid_var                         ! NetCDF forcing variable ID

  USE multiforce,only:forcefile                          ! name of forcing file
  USE multiforce,only:vname_aprecip                      ! variable name: precipitation
  USE multiforce,only:vname_airtemp                      ! variable name: temperature
  USE multiforce,only:vname_spechum                      ! variable name: specific humidity
  USE multiforce,only:vname_airpres                      ! variable name: surface pressure
  USE multiforce,only:vname_swdown                       ! variable name: downward shortwave radiation
  USE multiforce,only:vname_potevap                      ! variable name: potential ET
  USE multiforce,only:vname_q                            ! variable indice: observed discharge

  USE multiforce,only:ilook_aprecip                      ! variable indice: precipitation
  USE multiforce,only:ilook_airtemp                      ! variable indice: temperature
  USE multiforce,only:ilook_spechum                      ! variable indice: specific humidity
  USE multiforce,only:ilook_airpres                      ! variable indice: surface pressure
  USE multiforce,only:ilook_swdown                       ! variable indice: downward shortwave radiation
  USE multiforce,only:ilook_potevap                      ! variable indice: potential ET
  USE multiforce,only:ilook_q                            ! variable indice: observed discharge

  IMPLICIT NONE

  ! input
  integer(i4b), intent(in)                :: ncid        ! NetCDF file ID
  ! output
  integer(i4b), intent(out)               :: ierr        ! error code
  character(*), intent(out)               :: message     ! error message
  ! internal
  integer(i4b),parameter                  :: strLen=1024 ! length of character string
  type names
   character(len=strLen)                  :: vname       ! singlecharacter strings
  end type names
  type(names),dimension(nForce)           :: cVec        ! names of character strings
  integer(i4b)                            :: iVar        ! loop through forcing data

  ! ---------------------------------------------------------------------------------------
  ! initialize error control
  ierr=0; message='get_varID/'

  ! get the vector of variable names
  cVec(ilook_aprecip)%vname = trim(vname_aprecip)  ! variable name: precipitation
  cVec(ilook_potevap)%vname = trim(vname_potevap)  ! variable name: potential ET
  cVec(ilook_airtemp)%vname = trim(vname_airtemp)  ! variable name: temperature
  cVec(ilook_q)%vname = trim(vname_q)              ! variable name: observed discharge
  cVec(ilook_spechum)%vname = trim(vname_spechum)  ! variable name: specific humidity
  cVec(ilook_airpres)%vname = trim(vname_airpres)  ! variable name: surface pressure
  cVec(ilook_swdown)%vname  = trim(vname_swdown)   ! variable name: downward shortwave radiation

  !do ivar=1,nForce
  do ivar=1,4

    ! get the variable ID
    ierr = nf90_inq_varid(ncid, trim(cVec(iVar)%vname), ncid_var(ivar))
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'[variable='//trim(cVec(iVar)%vname)//']'; return; endif

  END DO

 END SUBROUTINE get_varID

 SUBROUTINE get_modtim(itim,ncid,ierr,message)
 ! ---------------------------------------------------------------------------------------
 ! Creator:
 ! --------
 ! Martyn Clark, 2012
 ! ---------------------------------------------------------------------------------------
 ! Purpose:
 ! --------
 ! Read NetCDF time variable for a given time step
 ! ---------------------------------------------------------------------------------------
 ! Modules Modified:
 ! -----------------
 ! MODULE multiforce -- populate structure timDat%(*)
 ! ---------------------------------------------------------------------------------------
 USE fuse_fileManager,only:INPUT_PATH                   ! defines data directory
 USE multiforce,only:forcefile                          ! name of forcing file
 USE multiforce,only:vname_dtime                        ! variable name: time since reference time
 USE multiforce,only:timDat                             ! time data strructure
 USE multiforce,only:jdayRef                            ! reference time (days)
 USE multiforce,only:latUnits,lonUnits,timeUnits        ! units string for time

 IMPLICIT NONE
 ! input
 integer(i4b), intent(in)                :: itim        ! index of model time step
 integer(i4b), intent(in)                :: ncid        ! NetCDF file ID
 ! output
 integer(i4b), intent(out)               :: ierr        ! error code
 character(*), intent(out)               :: message     ! error message
 ! internal
 integer(i4b),parameter                  :: strLen=1024 ! length of character string
 character(len=strLen)                   :: cmessage    ! error message of downwind routine
 integer(i4b)                            :: iVarID      ! NetCDF variable ID
 integer(i4b)                            :: iy,im,id    ! time of year
 integer(i4b)                            :: ih          ! time of day
 real(sp),dimension(1)                   :: atime       ! time array
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 ierr=0; message='get_modtim/'

 ! get variable ID for time
 ierr = nf90_inq_varid(ncid, trim(vname_dtime), iVarID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'[variable='//trim(vname_dtime)//']'; return; endif

 ! identify reference time
 call date_extractor(timeUnits,iy,im,id,ih)
 call juldayss(iy,im,id,ih,jdayRef,ierr,cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get the time
 ierr = nf90_get_var(ncid, iVarID, aTime, start=(/iTim/), count=(/1/))
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! put the time into the structure
 timDat%dtime = aTime(1)

 ! compute the year, month, day, hour, minute, second
 call caldatss(jdayRef+timDat%dtime,timDat%iy,timDat%im,timDat%id,timDat%ih,timDat%imin,timDat%dsec)
 !print*, timDat%iy,timDat%im,timDat%id!,timDat%ih,timDat%imin,timDat%dsec

 END SUBROUTINE get_modtim

 SUBROUTINE get_gforce(itim,ncid_forc,ierr,message)
 ! ---------------------------------------------------------------------------------------
 ! Creator:
 ! --------
 ! Martyn Clark, 2012
 ! ---------------------------------------------------------------------------------------
 ! Purpose:
 ! --------
 ! Read NetCDF gridded forcing data for a given time step
 ! ---------------------------------------------------------------------------------------
 ! Modules Modified:
 ! -----------------
 ! MODULE multiforce -- populate structure GFORCE(*,*)%(*)
 ! ---------------------------------------------------------------------------------------
 USE fuse_fileManager,only:INPUT_PATH                   ! defines data directory
 USE multiforce,only:forcefile                          ! name of forcing file
 USE multiforce,only:vname_aprecip                      ! variable name: precipitation
 USE multiforce,only:vname_airtemp                      ! variable name: temperature
 USE multiforce,only:vname_spechum                      ! variable name: specific humidity
 USE multiforce,only:vname_airpres                      ! variable name: surface pressure
 USE multiforce,only:vname_swdown                       ! variable name: downward shortwave radiation
 USE multiforce,only:vname_potevap                      ! variable name: potential ET

 USE multiforce,only:ilook_aprecip                      ! variable indice: precipitation
 USE multiforce,only:ilook_airtemp                      ! variable indice: temperature
 USE multiforce,only:ilook_spechum                      ! variable indice: specific humidity
 USE multiforce,only:ilook_airpres                      ! variable indice: surface pressure
 USE multiforce,only:ilook_swdown                       ! variable indice: downward shortwave radiation
 USE multiforce,only:ilook_potevap                      ! variable indice: potential ET

 USE multiforce,only:nspat1,nspat2                      ! dimension lengths
 USE multiforce,only:ncid_var                           ! NetCDF ID for forcing variables
 USE multiforce,only:amult_ppt,amult_pet                ! multipliers o convert to mm/day
 USE multiforce,only:gForce                             ! gridded forcing data
 USE multiforce,only:ancilF                             ! ancillary forcing data
 USE multiforce,only:nForce                             ! number of forcing variables

 IMPLICIT NONE
 ! input
 integer(i4b), intent(in)                :: itim        ! index of model time step
 integer(i4b), intent(in)                :: ncid_forc   ! NetCDF ID for the forcing file

 ! output
 integer(i4b), intent(out)               :: ierr        ! error code
 character(*), intent(out)               :: message     ! error message
 ! internal
 real(sp),parameter                      :: amiss=-9999._sp ! value for missing data
 integer(i4b),parameter                  :: strLen=1024 ! length of character string
 integer(i4b)                            :: iVar        ! loop through forcing data
 real(sp),dimension(:,:,:),allocatable   :: gTemp       ! temporary grid
 type names
  character(len=strLen)                  :: vname       ! singlecharacter strings
 end type names
 type(names),dimension(nForce)           :: cVec        ! names of character strings
 logical(lgt),dimension(nForce)          :: lCheck      ! check the existence of variables

 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 ierr=0; message='get_gforce/'
 ! ---------------------------------------------------------------------------------------

 ! initialize lCheck
 lCheck=.false.

 ! allocate space for the temporary grid
 allocate(gTemp(nSpat1,nSpat2,1), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating space for gTemp'; return; endif

 ! get the vector of variable names
 cVec(ilook_aprecip)%vname = trim(vname_aprecip)  ! variable name: precipitation
 cVec(ilook_potevap)%vname = trim(vname_potevap)  ! variable name: potential ET
 cVec(ilook_airtemp)%vname = trim(vname_airtemp)  ! variable name: temperature
 cVec(ilook_spechum)%vname = trim(vname_spechum)  ! variable name: specific humidity
 cVec(ilook_airpres)%vname = trim(vname_airpres)  ! variable name: surface pressure
 cVec(ilook_swdown)%vname  = trim(vname_swdown)   ! variable name: downward shortwave radiation

 ! get forcing grids
 ! do ivar=1,nForce
 do ivar=1,3

   ! get the data
   ierr = nf90_get_var(ncid_forc, ncid_var(ivar), gTemp, start=(/1,1,iTim/), count=(/nSpat1,nSpat2,1/)); CALL HANDLE_ERR(IERR)
   if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! save the data in the structure -- and convert fluxes to mm/day
  if(trim(cVec(iVar)%vname) == trim(vname_aprecip) )then; gForce(:,:)%ppt = gTemp(:,:,1)*amult_ppt; lCheck(ilook_aprecip) = .true.; endif
  if(trim(cVec(iVar)%vname) == trim(vname_potevap) )then; gForce(:,:)%pet = gTemp(:,:,1)*amult_pet; lCheck(ilook_potevap) = .true.; endif
  if(trim(cVec(iVar)%vname) == trim(vname_airtemp) )then; gForce(:,:)%temp = gTemp(:,:,1);       lCheck(ilook_airtemp) = .true.; endif

  ! save the other variables required to compute PET
  !if( trim(cVec(iVar)%vname) == trim(vname_airtemp) )then; ancilF(:,:)%airtemp = gTemp(:,:,1);       lCheck(ilook_airtemp) = .true.; endif
  !if( trim(cVec(iVar)%vname) == trim(vname_spechum) )then; ancilF(:,:)%spechum = gTemp(:,:,1);       lCheck(ilook_spechum) = .true.; endif
  !if( trim(cVec(iVar)%vname) == trim(vname_airpres) )then; ancilF(:,:)%airpres = gTemp(:,:,1);       lCheck(ilook_airpres) = .true.; endif
  !if( trim(cVec(iVar)%vname) == trim(vname_swdown)  )then; ancilF(:,:)%swdown  = gTemp(:,:,1);       lCheck(ilook_swdown)  = .true.; endif

 end do  ! (loop thru forcing variables)

 ! deallocate space for gTemp
 deallocate(gTemp, stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem deallocating space for gTemp'; return; endif

 end subroutine get_gforce

 SUBROUTINE get_gforce_3d(itim_start,numtim,ncid_forc,ierr,message)
 ! ---------------------------------------------------------------------------------------
 ! Creator:
 ! --------
 ! Nans Addor, based on Martyn Clark's get_gforce
 ! ---------------------------------------------------------------------------------------
 ! Purpose:
 ! --------
 ! Read NetCDF gridded forcing data for a range of time steps
 ! ---------------------------------------------------------------------------------------
 ! Modules Modified:
 ! -----------------
 ! MODULE multiforce -- populate structure GFORCE_3d(*,*)%(*)
 ! ---------------------------------------------------------------------------------------
 USE fuse_fileManager,only:INPUT_PATH                   ! defines data directory
 USE multiforce,only:forcefile                          ! name of forcing file
 USE multiforce,only:vname_aprecip                      ! variable name: precipitation
 USE multiforce,only:vname_airtemp                      ! variable name: temperature
 USE multiforce,only:vname_spechum                      ! variable name: specific humidity
 USE multiforce,only:vname_airpres                      ! variable name: surface pressure
 USE multiforce,only:vname_swdown                       ! variable name: downward shortwave radiation
 USE multiforce,only:vname_potevap                      ! variable name: potential ET
 USE multiforce,only:vname_q                            ! variable name: observed discharge

 USE multiforce,only:ilook_aprecip                      ! variable indice: precipitation
 USE multiforce,only:ilook_airtemp                      ! variable indice: temperature
 USE multiforce,only:ilook_spechum                      ! variable indice: specific humidity
 USE multiforce,only:ilook_airpres                      ! variable indice: surface pressure
 USE multiforce,only:ilook_swdown                       ! variable indice: downward shortwave radiation
 USE multiforce,only:ilook_potevap                      ! variable indice: potential ET
 USE multiforce,only:ilook_q                            ! variable indice: observed discharge

 USE multiforce,only:nspat1,nspat2                      ! dimension lengths
 USE multiforce,only:ncid_var                           ! NetCDF ID for forcing variables
 USE multiforce,only:amult_ppt,amult_pet                ! multipliers o convert to mm/day
 USE multiforce,only:gForce_3d                          ! gridded forcing data
 USE multiforce,only:ancilF_3d                          ! ancillary forcing data
 USE multiforce,only:nForce                             ! number of forcing variables
 USE multiforce,only:aValid                             ! time series of lumped forcing/response data

 IMPLICIT NONE
 ! input
 integer(i4b), intent(in)                :: itim_start  ! index of model time step - start of the period to extract
 integer(i4b), intent(in)                :: numtim      ! number of model time steps to extract
 integer(i4b), intent(in)                :: ncid_forc   ! NetCDF ID for the forcing file
 ! output
 integer(i4b), intent(out)               :: ierr        ! error code
 character(*), intent(out)               :: message     ! error message
 ! internal
 real(sp),parameter                      :: amiss=-9999._sp ! value for missing data
 integer(i4b),parameter                  :: strLen=1024 ! length of character string
 integer(i4b)                            :: iVar        ! loop through forcing data
 real(sp),dimension(:,:,:),allocatable   :: gTemp       ! temporary 3d grid
 type names
  character(len=strLen)                  :: vname       ! singlecharacter strings
 end type names
 type(names),dimension(nForce)           :: cVec        ! names of character strings
 logical(lgt),dimension(nForce)          :: lCheck      ! check the existence of variables

 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 ierr=0; message='get_gforce_3d/'
 ! ---------------------------------------------------------------------------------------

 ! initialize lCheck
 lCheck=.false.

 ! allocate space for the temporary grid
 allocate(gTemp(nSpat1,nSpat2,numtim), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating space for gTemp'; return; endif

 ! get the vector of variable names
 cVec(ilook_aprecip)%vname = trim(vname_aprecip)  ! variable name: precipitation
 cVec(ilook_potevap)%vname = trim(vname_potevap)  ! variable name: potential ET
 cVec(ilook_airtemp)%vname = trim(vname_airtemp)  ! variable name: temperature
 cVec(ilook_q)%vname = trim(vname_q)              ! variable name: observed discharge
 cVec(ilook_spechum)%vname = trim(vname_spechum)  ! variable name: specific humidity
 cVec(ilook_airpres)%vname = trim(vname_airpres)  ! variable name: surface pressure
 cVec(ilook_swdown)%vname  = trim(vname_swdown)   ! variable name: downward shortwave radiation

 ! get forcing grids
 ! do ivar=1,nForce

 do ivar=1,4

  ! get the data
  ierr = nf90_get_var(ncid_forc, ncid_var(ivar), gTemp, start=(/1,1,itim_start/), count=(/nSpat1,nSpat2,numtim/)); CALL HANDLE_ERR(IERR)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  !print *, trim(cVec(iVar)%vname)
  !print *, gTemp

  ! save the data in the structure -- and convert fluxes to mm/day
  if(trim(cVec(iVar)%vname) == trim(vname_aprecip) )then; gForce_3d(:,:,1:numtim)%ppt = gTemp(:,:,:)*amult_ppt; lCheck(ilook_aprecip) = .true.; endif
  if(trim(cVec(iVar)%vname) == trim(vname_potevap) )then; gForce_3d(:,:,1:numtim)%pet = gTemp(:,:,:)*amult_pet; lCheck(ilook_potevap) = .true.; endif
  if(trim(cVec(iVar)%vname) == trim(vname_airtemp) )then; gForce_3d(:,:,1:numtim)%temp = gTemp(:,:,:);       lCheck(ilook_airtemp) = .true.; endif
  if(trim(cVec(iVar)%vname) == trim(vname_q) )then;       aValid(:,:,1:numtim)%obsq = gTemp(:,:,:);       lCheck(ilook_q) = .true.; endif

  ! save the other variables required to compute PET
  !if( trim(cVec(iVar)%vname) == trim(vname_airtemp) )then; ancilF(:,:)%airtemp = gTemp(:,:,1);       lCheck(ilook_airtemp) = .true.; endif
  !if( trim(cVec(iVar)%vname) == trim(vname_spechum) )then; ancilF(:,:)%spechum = gTemp(:,:,1);       lCheck(ilook_spechum) = .true.; endif
  !if( trim(cVec(iVar)%vname) == trim(vname_airpres) )then; ancilF(:,:)%airpres = gTemp(:,:,1);       lCheck(ilook_airpres) = .true.; endif
  !if( trim(cVec(iVar)%vname) == trim(vname_swdown)  )then; ancilF(:,:)%swdown  = gTemp(:,:,1);       lCheck(ilook_swdown)  = .true.; endif

 end do  ! (loop thru forcing variables)

 ! deallocate space for gTemp
 deallocate(gTemp, stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem deallocating space for gTemp'; return; endif


  !print *,'aValid(1,1,1:10)',aValid(1,1,1:10)

 end subroutine get_gforce_3d

 subroutine date_extractor(refDate,iy,im,id,ih)
 ! used to extract the date from a units string
 ! (based on a routine written by David Rupp)
 IMPLICIT NONE
 ! input
 character(*),intent(in)             :: refDate     ! reference time string
 ! output
 integer(i4b),intent(out)            :: iy,im,id,ih ! reference time (year,month,day,hour)
 ! internal
 character(len=64)                   :: refd        ! temporary time and units string
 character(len=4)                    :: cyyyy       ! char year extracted from refDate
 character(len=2)                    :: cmm,cdd,chh ! char month and day and hour extracted from refDate
 integer(i4b)                        :: posit       ! position in refDate string

 ! strip out time units, if they exist (seconds since , days since , hours since )
 REFD = TRIM(REFDATE)
 POSIT = INDEX(REFDATE, 'since')
 IF (POSIT.GT.0) REFD = REFD(POSIT+6:len(refD))  ! +6 because 'since' has 5 characters
 ! get the year
 POSIT = INDEX(REFD, '-')     ! up to -
 CYYYY = REFD(1:POSIT-1)
 ! get the month
 REFD = REFD(POSIT+1:len(refD))
 POSIT = INDEX(REFD, '-')     ! up to -
 CMM = REFD(1:POSIT-1)
 ! get the day
 REFD = REFD(POSIT+1:len(refD))
 POSIT = INDEX(REFD, ' ')     ! up to space
 CDD = REFD(1:POSIT-1)
 ! get the hour
 REFD = REFD(POSIT+1:len(refD))
 POSIT = INDEX(REFD, ':')    ! up to :
 IF (POSIT.GT.0) THEN
  CHH = REFD(1:POSIT-1)
 ELSE
  CHH = '00'
 ENDIF
 ! convert to integers
 READ(CYYYY,'(i4)') IY
 READ(CMM,'(i2)') IM
 READ(CDD,'(i2)') ID
 READ(CHH,'(i2)') IH
 end subroutine date_extractor


 SUBROUTINE juldayss(iy,im,id,ih,            &  ! input
                     juldayFrac,ierr,message)   ! output
 !-----------------------------------------------------------------------
 ! compute the julian day
 !-----------------------------------------------------------------------
 !       Based on code from "Numerical Recipes in Fortran-77
 !
 !       Ref:    Press, W.H., S.A. Teukolsky, W.T. Vetterling, and B.P.
 !                 Flannery, 1992:  Numerical Recipes in Fortran 77:
 !                 The Art of Scientific Computing (2nd Ed.)  Cambridge
 !                 University Press, 933pp.
 !-----------------------------------------------------------------------
 !   Modified by David Rupp 2006-March-07 to account for hours with
 !   Output julian time in units of seconds from date
 !-----------------------------------------------------------------------
 IMPLICIT NONE
 ! input
 integer(i4b),intent(in)        :: iy,im,id,ih          ! year, month, day, hour
 ! output
 real(sp),intent(out)           :: juldayFrac           ! julian day (fraction of days)
 integer(i4b), intent(out)      :: ierr                 ! error code
 character(*), intent(out)      :: message              ! error message
 ! internal
 integer(i4b)                   :: julday               ! julian day (whole days)
 integer(i4b),parameter         :: IGREG=15+31*(10+12*1582) !IGREG = 588829
 integer(i4b)                   :: ja,jm,jy             ! temporary
 !-----------------------------------------------------------------------
 ! initialize error control
 !-----------------------------------------------------------------------
 ierr=0; message='juldayss/'
 !-----------------------------------------------------------------------
 jy=iy
 if (jy.eq.0)then; ierr=20; message=trim(message)//'there is no year zero'; return; endif
 if (jy.lt.0) jy=jy+1
 if (im.gt.2) then
  jm=im+1
 else
  jy=jy-1
  jm=im+13
 endif
 julday=int(365.25*jy)+int(30.6001*jm)+id+1720995
 if (id+31*(im+12*iy).ge.IGREG) then
  ja=int(0.01*jy)
  julday=julday+2-ja+int(0.25*ja)
 endif
 juldayFrac = real(julday, KIND(sp) ) &
            + real(ih,  KIND(sp) )/24._sp
 end SUBROUTINE juldayss

 SUBROUTINE caldatss(juliandd,iyyy,im,id,ih,imin,asec)
 !-----------------------------------------------------------------------
 !       Based on code from "Numerical Recipes in Fortran-77
 !
 !       Ref:    Press, W.H., S.A. Teukolsky, W.T. Vetterling, and B.P.
 !                 Flannery, 1992:  Numerical Recipes in Fortran 77:
 !                 The Art of Scientific Computing (2nd Ed.)  Cambridge
 !                 University Press, 933pp.
 !-----------------------------------------------------------------------
 !       Modified by David Rupp 2006-March-07 to account for hours
 !       Output is yyyy, mm, dd, hh
 !-----------------------------------------------------------------------
 IMPLICIT NONE
 ! input
 real(sp),intent(in)       :: juliandd   ! julian day in days since the beginning of time
 ! output
 integer(i4b),intent(out)  :: iyyy,im,id ! day of the year
 integer(i4b),intent(out)  :: ih,imin    ! time of day
 real(sp),intent(out)      :: asec       ! second
 ! internal
 integer(i4b)              :: julian     ! julian day
 integer(i4b),parameter    :: IGREG=2299161
 real(sp)                  :: hours, minutes
 integer(i4b)              :: ja,jalpha,jb,jc,jd,je
 ! get the julian day
 julian = int(juliandd)
 ! gets the hours, (remaining decimal)*24
 hours = (juliandd-julian)*24
 ih = int(hours)  ! convert to an integer
 ! get the minutes, (remaining decimal)*60
 minutes = (hours-ih)*60
 imin = int(minutes)
 ! get the seconds (keep as a decimal
 asec = (minutes-imin)*60
 ! uses the integer julian from above (below original num rec)
 if(julian.ge.IGREG)then
   jalpha=int(((julian-1867216)-0.25)/36524.25)
   ja=julian+1+jalpha-int(0.25*jalpha)
 else if(julian.lt.0)then
   ja=julian+36525*(1-julian/36525)
 else
   ja=julian
 endif
 jb=ja+1524
 jc=int(6680.+((jb-2439870)-122.1)/365.25)
 jd=365*jc+int(0.25*jc)
 je=int((jb-jd)/30.6001)
 id=jb-jd-int(30.6001*je)
 im=je-1
 if(im.gt.12)im=im-12
 iyyy=jc-4715
 if(im.gt.2)iyyy=iyyy-1
 if(iyyy.le.0)iyyy=iyyy-1
 if(julian.lt.0)iyyy=iyyy-100*(1-julian/36525)
 END SUBROUTINE caldatss

end module get_gforce_module
