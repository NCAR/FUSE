module get_gforce_module
USE nrtype
USE netcdf
USE time_io
implicit none
private
public::read_ginfo
public::get_dimIds
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
 USE multiforce,only:latitude,longitude                ! dimension arrays
 USE multiforce,only:time_steps,julian_time_steps      ! dimension arrays
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

 allocate(longitude(nspat1),latitude(nspat2),time_steps(numtim_in),julian_time_steps(numtim_in))

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

  ! save the data in the structure -- and convert fluxes to mm/day
  if(trim(cVec(iVar)%vname) == trim(vname_aprecip) )then
    if( ANY(gTemp(:,:,:).lt.0.0)) then; PRINT *, 'Negative precipitation in input file'; stop; endif
    gForce_3d(:,:,1:numtim)%ppt = gTemp(:,:,:)*amult_ppt; lCheck(ilook_aprecip) = .true.
  endif

  if(trim(cVec(iVar)%vname) == trim(vname_potevap) )then
    if( ANY(gTemp(:,:,:).lt.0.0)) then; PRINT *, 'Negative PET in input file'; stop; endif
    gForce_3d(:,:,1:numtim)%pet = gTemp(:,:,:)*amult_pet; lCheck(ilook_potevap) = .true.
  endif

  if(trim(cVec(iVar)%vname) == trim(vname_airtemp) )then
    gForce_3d(:,:,1:numtim)%temp = gTemp(:,:,:);       lCheck(ilook_airtemp) = .true.
  endif

  if(trim(cVec(iVar)%vname) == trim(vname_q) )then
    if( ANY(gTemp(:,:,:).lt.0.0)) then; PRINT *, 'Negative Q in input file'; stop; endif
    aValid(:,:,1:numtim)%obsq = gTemp(:,:,:);       lCheck(ilook_q) = .true.
  endif

  ! save the other variables required to compute PET
  !if( trim(cVec(iVar)%vname) == trim(vname_airtemp) )then; ancilF(:,:)%airtemp = gTemp(:,:,1);       lCheck(ilook_airtemp) = .true.; endif
  !if( trim(cVec(iVar)%vname) == trim(vname_spechum) )then; ancilF(:,:)%spechum = gTemp(:,:,1);       lCheck(ilook_spechum) = .true.; endif
  !if( trim(cVec(iVar)%vname) == trim(vname_airpres) )then; ancilF(:,:)%airpres = gTemp(:,:,1);       lCheck(ilook_airpres) = .true.; endif
  !if( trim(cVec(iVar)%vname) == trim(vname_swdown)  )then; ancilF(:,:)%swdown  = gTemp(:,:,1);       lCheck(ilook_swdown)  = .true.; endif

 end do  ! (loop thru forcing variables)

 ! deallocate space for gTemp
 deallocate(gTemp, stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem deallocating space for gTemp'; return; endif

 end subroutine get_gforce_3d


end module get_gforce_module
