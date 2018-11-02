module force_info_module
USE nrtype
USE netcdf
implicit none
private
public::force_info
contains

 SUBROUTINE force_info(fuse_mode,ierr,message)
 ! ---------------------------------------------------------------------------------------
 ! Creator:
 ! --------
 ! Martyn Clark, 2012
 ! Modified by Nans Addor to add numtim_sub for distributed modeling, 2017
 ! ---------------------------------------------------------------------------------------
 ! Purpose:
 ! --------
 ! Read information describing the forcing data file
 ! ---------------------------------------------------------------------------------------
 ! Modules Modified:
 ! -----------------
 ! MODULE multiforce -- populate variable names and time steps
 ! ---------------------------------------------------------------------------------------
 USE fuse_fileManager,only:SETNGS_PATH,FORCINGINFO,&        ! defines data directory
                           INPUT_PATH
 USE ascii_util_module,only:file_open                       ! open file (performs a few checks as well)
 USE ascii_util_module,only:get_vlines                      ! get a list of character strings from non-comment lines
 USE multiforce,only:vname_aprecip                          ! variable name: precipitation
 USE multiforce,only:vname_airtemp                          ! variable name: temperature
 USE multiforce,only:vname_spechum                          ! variable name: specific humidity
 USE multiforce,only:vname_airpres                          ! variable name: surface pressure
 USE multiforce,only:vname_swdown                           ! variable name: downward shortwave radiation
 USE multiforce,only:vname_potevap                          ! variable name: potential ET
 USE multiforce,only:vname_q                                ! variable name: runoff
 USE multiforce,only:vname_iy,vname_im,vname_id             ! names of time variables (day of year)
 USE multiforce,only:vname_ih,vname_imin,vname_dsec         ! names of time variables (time of day)
 USE multiforce,only:vname_dtime                            ! name of time variable (time since reference time)
 USE multiforce,only:deltim                                 ! model timestep (days)
 USE multiforce,only:istart,numtim_sim                      ! index for start of inference, and number steps in the reduced array
 USE multiforce,only:amult_ppt,amult_pet,amult_q            ! used to convert fluxes to mm/day
 USE multiforce,only:numtim_sub                             ! number of time steps of subperiod (will be kept in memory)

 IMPLICIT NONE
 ! input
 CHARACTER(LEN=10) , intent(in)         :: fuse_mode        ! fuse execution mode (run_def, run_best, calib_sce)
 ! output
 integer(i4b), intent(out)              :: ierr                 ! error code
 character(*), intent(out)              :: message              ! error message
 ! internal: general
 integer(i4b),parameter                 :: strLen=1024          ! length of character strings
 character(len=strLen)                  :: cmessage             ! message of downwind routine
 character(len=strLen),parameter        :: cVersion='FORCINGINFO.VERSION.2.2' ! version of forcinginfo file
 ! internal: read data from file
 integer(i4b)                           :: iunit                ! file unit
 character(len=strLen)                  :: cfile                ! name of control file
 character(len=strLen),allocatable      :: charlines(:)         ! vector of character strings
 ! internal: assign data
 integer(i4b)                           :: iLine                ! index of line in charlines
 integer(i4b)                           :: ibeg_name            ! start index of variable name in string charlines(iLine)
 integer(i4b)                           :: iend_name            ! end index of variable name in string charlines(iLine)
 integer(i4b)                           :: iend_data            ! end index of data in string charlines(iLine)
 character(len=strLen)                  :: cName,cData          ! name and data from charlines(iLine)
 ! internal: named variables
 integer(i4b),parameter                 :: maxinfo=22          ! maximum number of informational elements
 logical(lgt),dimension(maxinfo)        :: lCheck               ! vector to check that we have the infomation we need
 integer(i4b),parameter                 :: iVname_iy      =1    ! named variable for element of lCheck
 integer(i4b),parameter                 :: iVname_im      =2    ! named variable for element of lCheck
 integer(i4b),parameter                 :: iVname_id      =3    ! named variable for element of lCheck
 integer(i4b),parameter                 :: iVname_ih      =4    ! named variable for element of lCheck
 integer(i4b),parameter                 :: iVname_imin    =5    ! named variable for element of lCheck
 integer(i4b),parameter                 :: iVname_dsec    =6    ! named variable for element of lCheck
 integer(i4b),parameter                 :: iVname_dtime   =7    ! named variable for element of lCheck
 integer(i4b),parameter                 :: iVname_aprecip =8    ! named variable for element of lCheck
 integer(i4b),parameter                 :: iVname_airtemp =9   ! named variable for element of lCheck
 integer(i4b),parameter                 :: iVname_spechum =10   ! named variable for element of lCheck
 integer(i4b),parameter                 :: iVname_airpres =11   ! named variable for element of lCheck
 integer(i4b),parameter                 :: iVname_swdown  =12   ! named variable for element of lCheck
 integer(i4b),parameter                 :: iVname_potevap =13   ! named variable for element of lCheck
 integer(i4b),parameter                 :: iVname_q       =14   ! named variable for element of lCheck
 integer(i4b),parameter                 :: iUnits_aprecip =15   ! named variable for element of lCheck
 integer(i4b),parameter                 :: iUnits_airtemp =16   ! named variable for element of lCheck
 integer(i4b),parameter                 :: iUnits_spechum =17   ! named variable for element of lCheck
 integer(i4b),parameter                 :: iUnits_airpres =18   ! named variable for element of lCheck
 integer(i4b),parameter                 :: iUnits_swdown  =19   ! named variable for element of lCheck
 integer(i4b),parameter                 :: iUnits_potevap =20   ! named variable for element of lCheck
 integer(i4b),parameter                 :: iUnits_q       =21   ! named variable for element of lCheck
 integer(i4b),parameter                 :: iDeltim        =22   ! named variable for element of lCheck
 ! get units strings (used to define variable multipliers)
 character(len=strLen)                  :: units_aprecip='undefined' ! unit string for precipitation
 character(len=strLen)                  :: units_airtemp='undefined' ! unit string for air temperature
 character(len=strLen)                  :: units_spechum='undefined' ! unit string for specific humidity
 character(len=strLen)                  :: units_airpres='undefined' ! unit string for air pressure
 character(len=strLen)                  :: units_swdown ='undefined' ! unit string for downward sw radiation
 character(len=strLen)                  :: units_potevap='undefined' ! unit string for potential ET
 character(len=strLen)                  :: units_q      ='undefined' ! unit string for runoff
 integer(i4b)                           :: iVar                 ! loop through variables
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 ierr=0; message='force_info/'
 ! ---------------------------------------------------------------------------------------
 ! build filename
 cfile = trim(SETNGS_PATH)//trim(FORCINGINFO) ! uses paths and filenames from MODULE fuse_fileManager

 print *, 'Reading forcing info from:'
 print *, trim(cfile)
 ! open file (also returns un-used file unit used to open the file)
 call file_open(trim(cfile),iunit,ierr,cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! get a list of character strings from non-comment lines
 call get_vlines(iunit,charlines,ierr,cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! close the file unit
 close(iunit)
 ! ---------------------------------------------------------------------------------------
 ! initialize the check vector
 lCheck(:)=.false.
 ! loop through the non-comment lines in the input file
 do iLine=1,size(charlines)
  ! identify start and end of the name and the data
  ibeg_name = index(charlines(iLine),'<'); if(ibeg_name==0) ierr=20
  iend_name = index(charlines(iLine),'>'); if(iend_name==0) ierr=20
  iend_data = index(charlines(iLine),'!'); if(iend_data==0) ierr=20
  if(ierr/=0)then; message=trim(message)//'problem disentangling charlines(iLine) [string='//trim(charlines(iLine))//']'; return; endif
  ! extract name of the information, and the information itself
  cName = adjustl(charlines(iLine)(ibeg_name:iend_name))
  cData = adjustl(charlines(iLine)(iend_name+1:iend_data-1))
  ! put the information in its correct place
  select case(trim(cName))
   ! check version
   case('<version>')
    if(trim(cData)/=cVersion)then
     message=trim(message)//'version mis-match [version='//trim(cData)//'; expect "'//cVersion//'"]'
     ierr=20; return
    endif
   ! put character strings in their correct place
   case('<vname_iy>');      vname_iy      = trim(cData);              lCheck(iVname_iy)       = .true.
   case('<vname_im>');      vname_im      = trim(cData);              lCheck(iVname_im)       = .true.
   case('<vname_id>');      vname_id      = trim(cData);              lCheck(iVname_id)       = .true.
   case('<vname_ih>');      vname_ih      = trim(cData);              lCheck(iVname_ih)       = .true.
   case('<vname_imin>');    vname_imin    = trim(cData);              lCheck(iVname_imin)     = .true.
   case('<vname_dsec>');    vname_dsec    = trim(cData);              lCheck(iVname_dsec)     = .true.
   case('<vname_dtime>');   vname_dtime   = trim(cData);              lCheck(iVname_dtime)    = .true.
   case('<vname_aprecip>'); vname_aprecip = trim(cData);              lCheck(iVname_aprecip)  = .true.
   case('<vname_airtemp>'); vname_airtemp = trim(cData);              lCheck(iVname_airtemp)  = .true.
   case('<vname_spechum>'); vname_spechum = trim(cData);              lCheck(iVname_spechum)  = .true.
   case('<vname_airpres>'); vname_airpres = trim(cData);              lCheck(iVname_airpres)  = .true.
   case('<vname_swdown>');  vname_swdown  = trim(cData);              lCheck(iVname_swdown)   = .true.
   case('<vname_potevap>'); vname_potevap = trim(cData);              lCheck(iVname_potevap)  = .true.
   case('<vname_q>');       vname_q       = trim(cData);              lCheck(iVname_q)        = .true.
   case('<units_aprecip>'); units_aprecip = trim(cData);              lCheck(iUnits_aprecip)  = .true.
   case('<units_airtemp>'); units_airtemp = trim(cData);              lCheck(iUnits_airtemp)  = .true.
   case('<units_spechum>'); units_spechum = trim(cData);              lCheck(iUnits_spechum)  = .true.
   case('<units_airpres>'); units_airpres = trim(cData);              lCheck(iUnits_airpres)  = .true.
   case('<units_swdown>');  units_swdown  = trim(cData);              lCheck(iUnits_swdown)   = .true.
   case('<units_potevap>'); units_potevap = trim(cData);              lCheck(iUnits_potevap)  = .true.
   case('<units_q>');       units_q       = trim(cData);              lCheck(iUnits_q)        = .true.
   ! put real numbers and integers in their correct place
   case('<deltim>');     read(cData,*,iostat=ierr) deltim;            lCheck(iDeltim)         = .true.
   ! check for an unexpected string
   case default
    ierr=20; message=trim(message)//'do not have a case for string ['//trim(cName)//']'; return
  endselect
  ! check if there were any errors in the internal read statements
  if(ierr/=0)then
   message=trim(message)//'problem reading data for variable '//trim(cName)//'[data='//trim(cData)//']'
   ierr=50; return
  endif
 end do  ! (looping through non-comment lines in the file

 ! deallocate space for the variable line vector
 deallocate(charlines, stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem deallocating space for the variable line vector'; return; endif
 ! check that we got all desired variables
 if(any(lCheck .eqv. .false.))then
  ierr=20; message=trim(message)//'missing variable'
  write(*,'(a,1x,a,1x,L1)')    '<vname_iy>',      trim(vname_iy),      lCheck(iVname_iy)
  write(*,'(a,1x,a,1x,L1)')    '<vname_im>',      trim(vname_im),      lCheck(iVname_im)
  write(*,'(a,1x,a,1x,L1)')    '<vname_id>',      trim(vname_id),      lCheck(iVname_id)
  write(*,'(a,1x,a,1x,L1)')    '<vname_ih>',      trim(vname_ih),      lCheck(iVname_ih)
  write(*,'(a,1x,a,1x,L1)')    '<vname_imin>',    trim(vname_imin),    lCheck(iVname_im)
  write(*,'(a,1x,a,1x,L1)')    '<vname_dsec>',    trim(vname_dsec),    lCheck(iVname_dsec)
  write(*,'(a,1x,a,1x,L1)')    '<vname_dtime>',   trim(vname_dtime),   lCheck(iVname_dtime)
  write(*,'(a,1x,a,1x,L1)')    '<vname_aprecip>', trim(vname_aprecip), lCheck(iVname_aprecip)
  write(*,'(a,1x,a,1x,L1)')    '<vname_airtemp>', trim(vname_airtemp), lCheck(iVname_airtemp)
  write(*,'(a,1x,a,1x,L1)')    '<vname_spechum>', trim(vname_spechum), lCheck(iVname_spechum)
  write(*,'(a,1x,a,1x,L1)')    '<vname_airpres>', trim(vname_airpres), lCheck(iVname_airpres)
  write(*,'(a,1x,a,1x,L1)')    '<vname_swdown>',  trim(vname_swdown),  lCheck(iVname_swdown)
  write(*,'(a,1x,a,1x,L1)')    '<vname_potevap>', trim(vname_potevap), lCheck(iVname_potevap)
  write(*,'(a,1x,a,1x,L1)')    '<vname_q>',       trim(vname_q),       lCheck(iVname_q)
  write(*,'(a,1x,a,1x,L1)')    '<units_aprecip>', trim(units_aprecip), lCheck(iUnits_aprecip)
  write(*,'(a,1x,a,1x,L1)')    '<units_airtemp>', trim(units_airtemp), lCheck(iUnits_airtemp)
  write(*,'(a,1x,a,1x,L1)')    '<units_spechum>', trim(units_spechum), lCheck(iUnits_spechum)
  write(*,'(a,1x,a,1x,L1)')    '<units_airpres>', trim(units_airpres), lCheck(iUnits_airpres)
  write(*,'(a,1x,a,1x,L1)')    '<units_swdown>',  trim(units_swdown),  lCheck(iUnits_swdown)
  write(*,'(a,1x,a,1x,L1)')    '<units_potevap>', trim(units_potevap), lCheck(iUnits_potevap)
  write(*,'(a,1x,a,1x,L1)')    '<units_q>',       trim(units_q),       lCheck(iUnits_q)
  write(*,'(a,1x,f9.6,1x,L1)') '<deltim>',        deltim,              lCheck(iDeltim)
  print*, lCheck, size(lcheck)
  return
 endif  ! if we missed a variable

 ! get multipliers for each variable
 do ivar=1,3
  if(ivar==1) call get_multiplier(units_aprecip, amult_ppt, ierr, cmessage)
  if(ivar==2) call get_multiplier(units_potevap, amult_pet, ierr, cmessage)
  if(ivar==3) call get_multiplier(units_q,       amult_q,   ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do

 PRINT *, 'Done reading forcing info'

 end subroutine force_info

 ! ***** new subroutine: get multiplier for given flux variable (L/T)
 subroutine get_multiplier(cunits, amult, ierr, message)
 implicit none
 ! define input
 character(*),intent(in)                :: cunits    ! units
 ! define output
 real(sp),intent(out)                   :: amult     ! multiplier
 integer(i4b), intent(out)              :: ierr      ! error code
 character(*), intent(out)              :: message   ! error message
 ! define internal variables
 integer(i4b),parameter                 :: strLen=32 ! length of sub-strings
 integer(i4b)                           :: ipos      ! position of the "/" character
 character(strLen)                      :: cLength   ! length unit
 character(strLen)                      :: cTime     ! time unit
 real(sp),parameter                     :: secprday=86400._sp  ! number of seconds per day
 real(sp),parameter                     :: hrprday=24._sp      ! number of hours per day
 ! initialize error control
 ierr=0; message='get_multiplier/'
 ! if units are undefined, assume mm/day and have an early return
 if(trim(cunits)=='undefined')then; amult=1._sp; return; endif
 ! find the position of the "/" character
 ipos = index(trim(cunits),'/')
 ! check the "/" character exists
 if(ipos==0)then
  message=trim(message)//'expect the character "/" in the units string [units='//trim(cunits)//']'
  ierr=20; return
 endif
 ! get the length units
 cLength=cunits(1:ipos-1)
 if(cLength/='mm')then; ierr=20; message=trim(message)//'expect the length units to be "mm" [units='//trim(cLength)//']'; return; endif
 ! get the time units
 cTime=cunits(ipos+1:len_trim(cunits))
 ! get the multiplier
 select case(trim(cTime))
  case('d','day');    amult=1._sp
  case('h','hour');   amult=hrprday
  case('s','second'); amult=secprday
  case default
   ierr=20; message=trim(message)//'cannot identify the time units [time units = '//trim(cTime)//']'; return
 end select
 end subroutine get_multiplier

end module force_info_module
