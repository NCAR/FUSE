!-----------------------------------------------------------------------
!       Code from "Numerical Recipes in Fortran-77
!
!       Ref:    Press, W.H., S.A. Teukolsky, W.T. Vetterling, and B.P.
!                 Flannery, 1992:  Numerical Recipes in Fortran 77:
!                 The Art of Scientific Computing (2nd Ed.)  Cambridge
!                 University Press, 933pp.
!-----------------------------------------------------------------------
!       Modified by David Rupp 2006-March-07 to account for hours 
!       Output is yyyy, mm, dd, hh
!-----------------------------------------------------------------------
SUBROUTINE caldatss(julianss,iyyy,mm,id,ih,im,ss)
!SUBROUTINE caldat(julian,mm,id,iyyy)		
!INTEGER id,iyyy,julian,mm,IGREG			
INTEGER iyyy, mm, id, ih, im, julian, IGREG
DOUBLE PRECISION julianss, juliandd, hours, minutes, ss
PARAMETER (IGREG=2299161)
INTEGER ja,jalpha,jb,jc,jd,je

! gets the julian day in units of days since the beginning of time
juliandd = julianss / 86400
julian = int(juliandd)
! gets the hours, (remaining decimal)*24
hours = (juliandd-julian)*24
ih = int(hours)  ! convert to an integer
! get the minutes, (remaining decimal)*60  
minutes = (hours-ih)*60
im = int(minutes)
! get the seconds (keep as a decimal
ss = (minutes-im)*60

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
mm=je-1
if(mm.gt.12)mm=mm-12
iyyy=jc-4715
if(mm.gt.2)iyyy=iyyy-1
if(iyyy.le.0)iyyy=iyyy-1
if(julian.lt.0)iyyy=iyyy-100*(1-julian/36525)
ENDSUBROUTINE caldatss
