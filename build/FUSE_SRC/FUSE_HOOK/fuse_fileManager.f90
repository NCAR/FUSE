!******************************************************************
! (C) Copyright 2009-2010  ---  Dmitri Kavetski and Martyn Clark ---  All rights reserved
!******************************************************************
! Edited by Brian Henn to include snow model, 7/2013
! Edited by Nans Addor to set simulation and evaluation periods, 11/2017
MODULE fuse_filemanager
use kinds_dmsl_kit_FUSE,only:mik,mlk

implicit none
public
! FUSE-wide pathlength
integer(mik),parameter::fusePathLen=256
! defines the path for data files
CHARACTER(LEN=fusePathLen)  :: SETNGS_PATH
CHARACTER(LEN=fusePathLen)  :: INPUT_PATH
CHARACTER(LEN=fusePathLen)  :: OUTPUT_PATH
! content of input directory
CHARACTER(LEN=fusePathLen)  :: suffix_forcing    ! suffix for forcing file
CHARACTER(LEN=fusePathLen)  :: suffix_elev_bands ! suffix for elevation band file
! content of settings directory
CHARACTER(LEN=fusePathLen)  :: M_DECISIONS       ! definition of model decisions
CHARACTER(LEN=fusePathLen)  :: CONSTRAINTS       ! definition of parameter constraints
CHARACTER(LEN=fusePathLen)  :: MOD_NUMERIX       ! definition of numerical solution technique
CHARACTER(LEN=fusePathLen)  :: FORCINGINFO       ! info on forcing data files
CHARACTER(LEN=fusePathLen)  :: MBANDS_INFO       ! info on basin band data files ! not needed anymore
CHARACTER(LEN=fusePathLen)  :: MBANDS_NC         ! netcdf file defining the elevation bands
CHARACTER(LEN=fusePathLen)  :: BATEA_PARAM       ! definition of BATEA parameters ! remove this
! content of output directory
CHARACTER(LEN=64)           :: FMODEL_ID         ! string defining FUSE model
CHARACTER(LEN=64)           :: Q_ONLY_STR        ! TRUE = restrict attention to simulated runoff
LOGICAL                     :: Q_ONLY            ! .TRUE. = restrict attention to simulated runoff
! define simulation and evaluation periods
CHARACTER(len=20)           :: date_start_sim    ! date start simulation
CHARACTER(len=20)           :: date_end_sim      ! date end simulation
CHARACTER(len=20)           :: date_start_eval   ! date start evaluation period
CHARACTER(len=20)           :: date_end_eval     ! date end evaluation period
CHARACTER(len=20)           :: numtim_sub_str    ! number of time steps of subperiod (will be kept in memory)
! SCE parameters
CHARACTER(len=20)           :: KSTOP_str   ! number of shuffling loops the value must change by PCENTO
CHARACTER(len=20)           :: MAXN_str    ! maximum number of trials before optimization is terminated
CHARACTER(len=20)           :: PCENTO_str  ! the percentage

!----------------------------------------------------
contains
!----------------------------------------------------
subroutine fuse_SetDirsUndPhiles(fuseMusterDirektorIn,fuseFileManagerIn,err,message)
! Purpose: Sets direcotries and philenames for FUSE.
! ---
! Programmer: Dmitri Kavetski
! History:
!   Darby St,   18/10/2009 AD - leid out basik frammenverk
!   Sonnental,  17/06/2012 AD - more general path handling
! ---
! Usage
! fuseMusterDirektorIn  = master direktor file (path to filemanager)
! fuseFileManagerIn     = global names/path file
! ---
! Comments:
! 1. If present will try to use fuseMasterIn, otherwise default file.
!    if default not present in EXE path then uses default options
! ---
use utilities_dmsl_kit_FUSE,only:getSpareUnit
implicit none
! dummies
character(*),intent(in),optional::fuseMusterDirektorIn,fuseFileManagerIn
integer(mik),intent(out)::err
character(*),intent(out)::message
! registered settings
character(*),parameter::procnam="fuseSetDirsUndPhiles"
character(*),parameter::pathDelim="/\",defpathSymb="*",blank=" "
character(*),parameter::fuseMusterDirektorHeader="FUSE_MUSTERDIREKTOR_V1.0"
character(*),parameter::fuseFileManagerHeader="FUSE_FILEMANAGER_V1.5"
! locals
logical(mlk)::haveFMG,haveMUS
character(LEN=fusePathLen)::fuseMusterDirektor,fuseFileManager,defpath
character(LEN=100)::temp
integer(mik)::unt,i
! Start procedure here
err=0; message=procnam//"/ok"; defpath=blank
haveMUS=present(fuseMusterDirektorIn); haveFMG=present(fuseFileManagerIn)
if(haveMUS)haveMUS=len_trim(fuseMusterDirektorIn)>0
if(haveFMG)haveFMG=len_trim(fuseFileManagerIn)>0  ! check for zero-string
if(haveMUS.and.haveFMG)then
  message="f-"//procnam//"/mustSpecifyEither(notBoth)&
           &[fuseMusterDirektor.or.fuseFileManager]"
  err=10; return
elseif(haveFMG)then
  fuseFileManager=fuseFileManagerIn
  i=scan(fuseFileManager,pathDelim,back=.true.)
  if(i>0)defpath=fuseFileManager(:i-1)//pathDelim(1:1)
  print *, 'fuseFileManager:', TRIM(fuseFileManager)

elseif(haveMUS)then
  fuseMusterDirektor=fuseMusterDirektorIn
  i=scan(fuseMusterDirektor,pathDelim,back=.true.)
  if(i>0)defpath=fuseMusterDirektor(:i-1)//pathDelim(1:1)
  print *, 'fuseMusterDirektor:', TRIM(fuseMusterDirektor)

else
  message="f-"//procnam//"/mustSpecifyEither&
           &[fuseMusterDirektor.or.fuseFileManager]"
  err=20; return
endif
call getSpareUnit(unt,err,message) ! make sure 'unt' is actually available
if(err/=0)then
  message="f-"//procnam//"/weird/&"//message
  err=100; return
endif
if(.not.haveFMG)then  ! grab it from the muster-direktor

! 2. Open muster-direktor and read it
  open(unt,file=fuseMusterDirektor,status="old",action="read",iostat=err)
  if(err/=0)then
    message="f-"//procnam//"/musterDirektorFileOpenError['"//trim(fuseMusterDirektor)//"']"
    err=10; return
  endif
  read(unt,*)temp
  if(temp/=fuseMusterDirektorHeader)then
    message="f-"//procnam//"/unknownHeader&[file='"//trim(fuseMusterDirektor)//"']&&
      &[header='"//trim(temp)//"']"
    err=20; return
  endif
  read(unt,*)fuseFileManager
  close(unt)
endif
! open file manager file
open(unt,file=fuseFileManager,status="old",action="read",iostat=err)
if(err/=0)then
  message="f-"//procnam//"/fileManagerOpenError['"//trim(fuseFileManager)//"']"
  err=10; return
endif
read(unt,*)temp
if(temp/=fuseFileManagerHeader)then
  message="f-"//procnam//"/unknownHeader&[file='"//trim(fuseFileManager)//"']&&
    &[header="//trim(temp)//"]"

  message='This version of FUSE requires the file manager to follow the following format:  '//trim(fuseFileManagerHeader)//' not '//trim(temp)

  err=20; return
endif
read(unt,'(a)')temp
read(unt,*)SETNGS_PATH
read(unt,*)INPUT_PATH
read(unt,*)OUTPUT_PATH
read(unt,'(a)')temp
read(unt,*)suffix_forcing
read(unt,*)suffix_elev_bands
read(unt,'(a)')temp
read(unt,*)FORCINGINFO
read(unt,*)CONSTRAINTS
read(unt,*)MOD_NUMERIX
read(unt,*)M_DECISIONS
read(unt,'(a)')temp
read(unt,*)FMODEL_ID
read(unt,*)Q_ONLY_STR
read(unt,'(a)')temp
read(unt,*)date_start_sim
read(unt,*)date_end_sim
read(unt,*)date_start_eval
read(unt,*)date_end_eval
read(unt,*)numtim_sub_str
read(unt,'(a)')temp
read(unt,*)MAXN_STR
read(unt,*)KSTOP_STR
read(unt,*)PCENTO_STR
close(unt)

! Convert Q_ONLY to logical
if(Q_ONLY_STR=='TRUE')then
  Q_ONLY = .TRUE.
elseif(Q_ONLY_STR=='FALSE')then
  Q_ONLY = .FALSE.
else
  message="Q_ONLY must be either TRUE or FALSE"
  err=20; return
endif

PRINT*, 'Q_ONLY', Q_ONLY

! process paths a bit
if(SETNGS_PATH(1:1)==defpathSymb)SETNGS_PATH=trim(defpath)//SETNGS_PATH(2:)
if( INPUT_PATH(1:1)==defpathSymb) INPUT_PATH=trim(defpath)//INPUT_PATH (2:)
if(OUTPUT_PATH(1:1)==defpathSymb)OUTPUT_PATH=trim(defpath)//OUTPUT_PATH(2:)

PRINT *, 'Paths defined in file manager:'
PRINT *, 'SETNGS_PATH:', TRIM(SETNGS_PATH)
PRINT *, 'INPUT_PATH:', TRIM(INPUT_PATH)
PRINT *, 'OUTPUT_PATH:', TRIM(OUTPUT_PATH)

PRINT *, 'Dates defined in file manager:'
PRINT *, 'date_start_sim:', TRIM(date_start_sim)
PRINT *, 'date_end_sim:', TRIM(date_end_sim)
PRINT *, 'date_start_eval:', TRIM(date_start_eval)
PRINT *, 'date_end_eval:', TRIM(date_end_eval)
PRINT *, 'numtim_sub_str:', TRIM(numtim_sub_str)

! End procedure here
endsubroutine fuse_SetDirsUndPhiles
!----------------------------------------------------
END MODULE fuse_filemanager
