!******************************************************************
! (C) Copyright 2000-2020  ---  Dmitri Kavetski  ---  All rights reserved
! NB: CUSTOMIZED VERSION FOR FUSE SUITE OF MARTYN CLARK
!******************************************************************
module utilities_dmsl_kit_FUSE
! Purpose: Suite of DMSL-originated utilities used in FUSE
! Programmer: Dmitri Kavetski
! Last modified: 14 Sept 2014 Ad - Chiefleys
! Comments:
use kinds_dmsl_kit_FUSE
implicit none
!----------------------------------------------------
character(*),parameter::blankCH=" "
!----------------------------------------------------
contains
!----------------------------------------------------
subroutine getSpareUnit(unt,err,message)
! Purpose: Finds first-available unit for opening a file.
! Programmer: Dmitri Kavetski, circa 2000 AD - Newcastle
implicit none
! dummies
integer(mik),intent(inout)::unt
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
integer(mik)::i
logical(mlk)::opened,xist
integer(mik),parameter::minUnitsDef=7,maxUnits=2147483639
integer(mik)::minUnits
! Start procedure here
if(unt>minUnitsDef)then; minUnits=unt
else;                    minUnits=minUnitsDef; endif
do i=minUnits,maxUnits
  inquire(unit=i,opened=opened,exist=xist) ! check unit status
  if(.not.opened.and.xist)then ! un-opened existing unit found
    !print *, 'Free unit found ', i
    unt=i
    err=0
    !message="getSpareUnit/ok"
    exit
  endif
enddo
if(i>maxUnits)then  ! all units in use
  unt=-1; err=-10; message="getSpareUnit/allUnitsInUse&"//&
      "&(all 2.2billion-u've goda b jokin')"
endif
! End procedure here
endsubroutine getSpareUnit
!----------------------------------------------------
elemental subroutine stripTrailString(string,trailStart,trail)
! Purpose: trims string to remove a trailing substring
! that start with substring "trailStart"
! Typical use is to remove trailing comment from a character string
! Comments:
! 1. trailStart is searched for, not trim(trailStart), to allow blank parsing
! 2. resulting string padded with blanks
implicit none
! dummies
character(*),intent(inout)::string
character(*),intent(in)::trailStart
character(*),intent(out),optional::trail
! locals
integer(mik)::i
! Start procedure here
if(len(trailStart)==0)return ! do nothing if trailstart is nothing
i=index(string,trailStart)
selectcase(i)
case(0)   ! not encountered
  if(present(trail))trail=blankCH
case(1)   ! first thing encountered
  if(present(trail))trail=string
  string=blankCH
case(2:)  ! encountered later in string
  if(present(trail))trail=string(i:)
  string(i:)=blankCH
endselect
! End procedure here
endsubroutine stripTrailString
!----------------------------------------------------
endmodule utilities_dmsl_kit_FUSE
!******************************************************************
