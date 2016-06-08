  FUNCTION pythag_sp(a,b)
  USE nrtype
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: a,b
  REAL(SP) :: pythag_sp
  REAL(SP) :: absa,absb
  absa=abs(a)
  absb=abs(b)
  if (absa > absb) then
    pythag_sp=absa*sqrt(1.0_sp+(absb/absa)**2)
  else
    if (absb == 0.0) then
      pythag_sp=0.0
    else
      pythag_sp=absb*sqrt(1.0_sp+(absa/absb)**2)
    end if
  end if
  END FUNCTION pythag_sp
