!-----------------------------------------------------------------------
!       Code from "Numerical Recipes in Fortran-77
!
!       Ref:    Press, W.H., S.A. Teukolsky, W.T. Vetterling, and B.P.
!                 Flannery, 1992:  Numerical Recipes in Fortran 77:
!                 The Art of Scientific Computing (2nd Ed.)  Cambridge
!                 University Press, 933pp.
!-----------------------------------------------------------------------
!	Modified by David Rupp 2006-March-07 to account for hours with
!	Output julian time in units of seconds from date
!-----------------------------------------------------------------------
FUNCTION juldayss(yyin,mmin,ddin,hhin)
INTEGER julday,iyyy,mm,id,ih,IGREG
INTEGER yyin,mmin,ddin,hhin
DOUBLE PRECISION juldayss
PARAMETER (IGREG=15+31*(10+12*1582))  !IGREG = 588829
INTEGER ja,jm,jy

iyyy = yyin
mm= mmin
id = ddin
ih = hhin

jy=iyyy
if (jy.eq.0) then
  write(*,*) 'julday: there is no year zero'
  stop
endif
if (jy.lt.0) jy=jy+1
if (mm.gt.2) then
  jm=mm+1
else
  jy=jy-1
  jm=mm+13
endif
julday=int(365.25*jy)+int(30.6001*jm)+id+1720995
if (id+31*(mm+12*iyyy).ge.IGREG) then
  ja=int(0.01*jy)
  julday=julday+2-ja+int(0.25*ja)
endif

juldayss = 86400.0D0*real(julday, KIND(JULDAYSS) ) &
                   + real(ih,  KIND(JULDAYSS) )*3600.0D0

return
END
