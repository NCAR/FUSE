!******************************************************************
! (C) Copyright 2009-2010  ---  Dmitri Kavetski and Martyn Clark ---  All rights reserved
!******************************************************************
! Edited by Brian Henn to include snow model, 7/2013
MODULE fuse_filemanager
use kinds_dmsl_kit_FUSE,only:mik,mlk
implicit none
public
! FUSE-wide pathlength
integer(mik),parameter::fusePathLen=256
! defines the path for data files (and default values)
CHARACTER(LEN=fusePathLen)  :: SETNGS_PATH='/home/naddor/fuse/settings/'
CHARACTER(LEN=fusePathLen)  :: INPUT_PATH ='/home/naddor/fuse/input/'
CHARACTER(LEN=fusePathLen)  :: OUTPUT_PATH='/home/naddor/fuse/output/'
! define name of control files    (and default values)
CHARACTER(LEN=fusePathLen)  :: M_DECISIONS='fuse_zDecisions.txt'    ! definition of model decisions
CHARACTER(LEN=fusePathLen)  :: CONSTRAINTS='fuse_zConstraints.txt'  ! definition of parameter constraints
CHARACTER(LEN=fusePathLen)  :: MOD_NUMERIX='fuse_zNumerix.txt'      ! definition of numerical solution technique
! additional control files (not needed by the FUSE engines)
CHARACTER(LEN=fusePathLen)  :: FORCINGINFO='forcing.info.txt'  ! info on forcing data files
CHARACTER(LEN=fusePathLen)  :: MBANDS_INFO='elevbands.info.txt'  ! info on basin band data files
CHARACTER(LEN=fusePathLen)  :: BATEA_PARAM='batea_param.txt'  ! definition of BATEA parameters
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
character(*),parameter::fuseFileManagerHeader="FUSE_FILEMANAGER_V1.1"
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
  err=20; return
endif
read(unt,'(a)')temp
read(unt,'(a)')temp
read(unt,*)SETNGS_PATH
read(unt,*)INPUT_PATH
read(unt,*)OUTPUT_PATH
read(unt,'(a)')temp
read(unt,*)FORCINGINFO
read(unt,*)MBANDS_INFO
read(unt,*)M_DECISIONS
read(unt,*)CONSTRAINTS
read(unt,*)MOD_NUMERIX
read(unt,*)BATEA_PARAM
close(unt)
! process paths a bit
if(SETNGS_PATH(1:1)==defpathSymb)SETNGS_PATH=trim(defpath)//SETNGS_PATH(2:)
if( INPUT_PATH(1:1)==defpathSymb) INPUT_PATH=trim(defpath)//INPUT_PATH (2:)
if(OUTPUT_PATH(1:1)==defpathSymb)OUTPUT_PATH=trim(defpath)//OUTPUT_PATH(2:)
! End procedure here
endsubroutine fuse_SetDirsUndPhiles
!----------------------------------------------------
END MODULE fuse_filemanager
