module getf_ascii_module
USE nrtype
USE netcdf
implicit none
private
public:: prelim_asc
public:: close_file
public:: read_ascii
contains

 SUBROUTINE prelim_asc(ierr,message)
 ! ---------------------------------------------------------------------------------------
 ! Creator:
 ! --------
 ! Martyn Clark, 2012
 ! ---------------------------------------------------------------------------------------
 ! Purpose:
 ! --------
 ! Find columns for the input data
 ! ---------------------------------------------------------------------------------------
 ! Modules Modified:
 ! -----------------
 ! MODULE multiforce -- populates vairable indices
 ! ---------------------------------------------------------------------------------------
 USE fuse_fileManager,only:INPUT_PATH                       ! defines data directory
 USE ascii_util_module,only:file_open                       ! open file (performs a few checks as well)
 USE ascii_util_module,only:split_line                      ! get a vector of "words" from a line of character data
 USE multiforce,only:forcefile                              ! name of forcing file
 USE multiforce,only:vname_aprecip,vname_potevap,vname_q    ! names of variables in forcing file
 USE multiforce,only:vname_iy,vname_im,vname_id             ! names of time variables (day of year)
 USE multiforce,only:vname_ih,vname_imin,vname_dsec         ! names of time variables (time of day)
 USE multiforce,only:ivarid_iy,ivarid_im,ivarid_id          ! index of columns for day of year
 USE multiforce,only:ivarid_ih,ivarid_imin,ivarid_dsec      ! index of columns for time of day
 USE multiforce,only:ivarid_ppt,ivarid_pet,ivarid_q         ! index of columns for data
 USE multiforce,only:sim_beg                                ! index of the start of the simulation
 USE multiforce,only:deltim                                 ! time step of the data (data interval)
 IMPLICIT NONE
 ! output
 integer(i4b), intent(out)              :: ierr             ! error code
 character(*), intent(out)              :: message          ! error message
 ! internal: general
 integer(i4b),parameter                 :: strLen=1024      ! length of character strings
 character(len=strLen)                  :: cmessage         ! message of downwind routine
 ! internal: open file
 integer(i4b)                           :: iunit            ! file unit
 character(len=strLen)                  :: infile           ! name of input file
 ! internal: get column indices
 integer(i4b)                           :: ihead            ! loop through header lines
 integer(i4b)                           :: idata            ! loop through lines of data
 integer(i4b)                           :: iword            ! loop through words on a given line of data
 integer(i4b),parameter                 :: maxhead=100      ! maximum number of header lines
 integer(i4b),parameter                 :: maxdata=100000   ! maximum number of data lines
 character(len=strLen)                  :: cLine            ! a single line of data
 character(len=strLen),allocatable      :: words(:)         ! vector of "words" on a given line
 integer(i4b)                           :: nhead            ! number of header lines
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 ierr=0; message='prelim_asc/'
 ! ---------------------------------------------------------------------------------------
 ! build filename
 infile = trim(INPUT_PATH)//trim(forcefile) ! uses paths and filenames from MODULE fuse_fileManager
 ! open file (return file unit, and return errors if the file is already open, or if the file does not exist
 call file_open(trim(infile),iunit,ierr,cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! --------------------------------------------------------------------------------------
 ! (1) get column indices
 ! --------------------------------------------------------------------------------------
 do ihead=1,maxhead+1  ! loop through header lines
  ! read header line
  read(iunit,'(a)', iostat=ierr) cline
  ! check read statement
  if(ierr/=0)then
   write(message,'(a,i0,a)') trim(message)//'problem with data read after data line ',ihead,&
                             '; most likely reached the end of the file [file='//trim(infile)//']'
   return
  endif
  ! check that we found the header line within the allocated number of lines
  if(ihead==maxhead+1)then
   write(message,'(a,i0,a)') trim(message)//'cannot find header: looped through ',maxhead,' lines'
   ierr=20; return
  endif
  ! ignore comment lines (lines starting with "!")
  if(cline(1:1)=='!') cycle ! line of data is a comment
  ! ***** if get to here, found the header line
  ! split the line into words
  call split_line(cline,words,ierr,cmessage)
  if(ierr/=0)then; ierr=20; message=trim(message)//trim(cmessage); endif
  ! find the column indes for each desired variable
  do iword=1,size(words)
   ! get the time variables
   if(trim(words(iword))==trim(vname_iy))      ivarid_iy   = iword    ! year
   if(trim(words(iword))==trim(vname_im))      ivarid_im   = iword    ! month
   if(trim(words(iword))==trim(vname_id))      ivarid_id   = iword    ! day
   if(trim(words(iword))==trim(vname_ih))      ivarid_ih   = iword    ! hour
   if(trim(words(iword))==trim(vname_imin))    ivarid_imin = iword    ! minute
   if(trim(words(iword))==trim(vname_dsec))    ivarid_dsec = iword    ! second
   ! get the data variables
   if(trim(words(iword))==trim(vname_aprecip)) ivarid_ppt  = iword    ! precipitation
   if(trim(words(iword))==trim(vname_potevap)) ivarid_pet  = iword    ! potential ET
   if(trim(words(iword))==trim(vname_q))       ivarid_q    = iword    ! runoff
  end do
  ! deallocate the words vector
  deallocate(words, stat=ierr); if(ierr/=0)then; message=trim(message)//'problem deallocating words vector'; return; endif
  ! check found the columns for date variables
  if(ivarid_iy < 0 .or. ivarid_im < 0 .or. ivarid_id < 0)then
   message=trim(message)//'cannot find variables [vnames = ('//trim(vname_iy)//', '//trim(vname_im)//', '//trim(vname_id)//'); cLine = '//trim(cLine)
   ierr=20; return
  endif
  ! check found the columns for the time variables
  if(deltim < 1._sp        .and. ivarid_ih < 0)then; ierr=30; message=trim(message)//'cannot find variable "'//trim(vname_ih)//'"'; return; endif
  if(deltim < 1._sp/24._dp .and. ivarid_im < 0)then; ierr=40; message=trim(message)//'cannot find variable "'//trim(vname_im)//'"'; return; endif
  ! check found the columns for data variables
  if(ivarid_ppt < 0 .or. ivarid_pet < 0 .or. ivarid_q < 0)then
   message=trim(message)//'cannot find variables [vnames = ('//trim(vname_aprecip)//', '//trim(vname_potevap)//', '//trim(vname_q)//'); cLine = '//trim(cLine)
   ierr=20; return
  endif
  ! save the number of header lines
  nhead=ihead
  ! exit here (assume first non-comment line is the header)
  exit
 end do  ! (loop through header lines)
 ! --------------------------------------------------------------------------------------
 ! (2) read to the start of the warm-up period
 ! --------------------------------------------------------------------------------------
 if(sim_beg==1)return ! no need to read lines of data
 ! loop through data lines
 do idata=1,maxdata
  ! read a line of data and check that it is not a comment line
  read(iunit,*,iostat=ierr) cline              ! read a line of data
  ! check read statement
  if(ierr/=0)then
   write(message,'(a,i0,a)') trim(message)//'problem with data read after data line ',idata,&
                             '; most likely reached the end of the file [file='//trim(infile)//']'
   return
  endif
  ! check there are no comment lines in the data block
  if(cline(1:1)=='!')then; ierr=20; message=trim(message)//'do not allow comment lines after the header [file='//trim(infile)//']'; return; endif
  ! exit do loop if reached start of warm-up
  if(idata == warmup_beg-1) exit   ! exit if one line before the start of the warm-up period
  ! check that can find the warm-up period
  if(idata==maxdata)then
   write(message,'(a,i0,a)') trim(message)//'cannot find start of warm-up period: looped through ',maxdata,' lines'
   ierr=20; return
  endif
 end do ! (end looping through data lines)

 end subroutine prelim_asc


 SUBROUTINE close_file(ierr,message)
 ! ---------------------------------------------------------------------------------------
 ! Creator:
 ! --------
 ! Martyn Clark, 2012
 ! ---------------------------------------------------------------------------------------
 ! Purpose:
 ! --------
 ! Close file
 ! ---------------------------------------------------------------------------------------
 ! Modules Modified:
 ! -----------------
 ! None
 ! ---------------------------------------------------------------------------------------
 USE fuse_fileManager,only:INPUT_PATH                       ! defines data directory
 USE multiforce,only:forcefile                              ! name of forcing file
 IMPLICIT NONE
 ! output
 integer(i4b), intent(out)              :: ierr             ! error code
 character(*), intent(out)              :: message          ! error message
 ! internal
 integer(i4b),parameter                 :: strLen=1024      ! length of character strings
 integer(i4b)                           :: iunit            ! file unit
 character(len=strLen)                  :: infile           ! name of input file
 logical(lgt)                           :: xist             ! .true. if the file exists
 logical(lgt)                           :: xopn             ! .true. if the file is open
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 ierr=0; message='close_file/'
 ! ---------------------------------------------------------------------------------------
 ! build filename
 infile = trim(INPUT_PATH)//trim(forcefile) ! uses paths and filenames from MODULE fuse_fileManager
 ! check that the file is open, and get the file unit
 inquire(file=trim(infile),exist=xist,opened=xopn, number=iunit) ! Check if the file is open and exists
 if(.not.xist)then; ierr=10; message=trim(message)//'file "'//trim(infile)//'" does not exist!'; return; endif
 if(.not.xopn)then; ierr=20; message=trim(message)//'file "'//trim(infile)//'" is not open!'; return; endif
 ! close file
 close(iunit)
 end subroutine close_file


 SUBROUTINE read_ascii(ierr,message)
 ! ---------------------------------------------------------------------------------------
 ! Creator:
 ! --------
 ! Martyn Clark, 2012
 ! ---------------------------------------------------------------------------------------
 ! Purpose:
 ! --------
 ! Read a line of ASCII data and populate forcing structures
 ! ---------------------------------------------------------------------------------------
 ! Modules Modified:
 ! -----------------
 ! MODULE multiforce -- populates forcing structures
 ! ---------------------------------------------------------------------------------------
 USE fuse_fileManager,only:INPUT_PATH                       ! defines data directory
 USE ascii_util_module,only:split_line                      ! get a vector of "words" from a line of character data
 USE multiforce,only:forcefile                              ! name of forcing file
 USE multiforce,only:ivarid_iy,ivarid_im,ivarid_id          ! index of columns for day of year
 USE multiforce,only:ivarid_ih,ivarid_imin,ivarid_dsec      ! index of columns for time of day
 USE multiforce,only:ivarid_ppt,ivarid_pet,ivarid_q         ! index of columns for the data
 USE multiforce,only:gForce                                 ! gridded forcing data
 USE multiforce,only:valDat                                 ! response data
 USE multiforce,only:timDat                                 ! time structure
 IMPLICIT NONE
 ! output
 integer(i4b), intent(out)              :: ierr             ! error code
 character(*), intent(out)              :: message          ! error message
 ! internal: general
 integer(i4b),parameter                 :: strLen=1024      ! length of character strings
 character(len=strLen)                  :: cmessage         ! message of downwind routine
 ! internal: get the file unit
 logical(lgt)                           :: xist             ! .true. if the file exists
 logical(lgt)                           :: xopn             ! .true. if the file is open
 integer(i4b)                           :: iunit            ! file unit
 character(len=strLen)                  :: infile           ! name of input file
 ! internal: read a line of data
 character(len=strLen)                  :: cLine            ! a single line of data
 character(len=strLen),allocatable      :: words(:)         ! vector of "words" on a given line
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 ierr=0; message='read_ascii/'
 ! ---------------------------------------------------------------------------------------
 ! build filename
 infile = trim(INPUT_PATH)//trim(forcefile) ! uses paths and filenames from MODULE fuse_fileManager
 ! check that the file is open, and get the file unit
 inquire(file=trim(infile),exist=xist,opened=xopn, number=iunit) ! Check if the file is open and exists
 if(.not.xist)then; ierr=10; message=trim(message)//'file "'//trim(infile)//'" does not exist!'; return; endif
 if(.not.xopn)then; ierr=20; message=trim(message)//'file "'//trim(infile)//'" is not open!'; return; endif
 ! ---------------------------------------------------------------------------------------
 ! read a line of data
 read(iunit,'(a)',iostat=ierr) cline              ! read a line of data
 ! split the line into words
 call split_line(cline,words,ierr,cmessage)
 if(ierr/=0)then; ierr=20; message=trim(message)//trim(cmessage); endif
 ! put the data in the structure: day of year
 read(words(ivarid_iy),*) timDat%iy
 read(words(ivarid_im),*) timDat%im
 read(words(ivarid_id),*) timDat%id
 ! put the data in the structure: time of day
 if(ivarid_ih   > 0) read(words(ivarid_ih),*)   timDat%ih
 if(ivarid_imin > 0) read(words(ivarid_imin),*) timDat%imin
 if(ivarid_dsec > 0) read(words(ivarid_dsec),*) timDat%dsec
 ! put the data in the structure: forcing data
 read(words(ivarid_ppt),*) gForce(1,1)%ppt
 read(words(ivarid_pet),*) gForce(1,1)%pet
 ! put the data in the structure: response data
 read(words(ivarid_q),*)   valDat%obsq
 ! deallocate the words vector
 deallocate(words, stat=ierr); if(ierr/=0)then; message=trim(message)//'problem deallocating words vector'; return; endif
 end subroutine read_ascii


end module getf_ascii_module
