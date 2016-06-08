SUBROUTINE PAR_DERIVE(err,message)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Computes derived model parameters (bucket sizes, etc.)
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multiparam -- model parameters stored in MODULE multiparam
! ---------------------------------------------------------------------------------------
USE nrtype                                           ! define data types
USE model_defn, ONLY: SMODL                          ! model definition structures
USE model_defnames
USE multiparam, ONLY: MPARAM,DPARAM                  ! model parameter structures
IMPLICIT NONE
! dummies
integer(i4b),intent(out)::err
character(*),intent(out)::message
! ---------------------------------------------------------------------------------------
err=0
CALL BUCKETSIZE()        ! compute bucket size
CALL MEAN_TIPOW()        ! mean of the power-transformed topo index
CALL QBSATURATN()        ! compute baseflow at saturation (used in the SAC percolation model)
CALL QTIMEDELAY(err,message)        ! compute fraction of runoff in future time steps
if(err/=0)then
  err=10; message="f-PAR_DERIVE/&"//trim(message); return
endif
! ---------------------------------------------------------------------------------------
IF (SMODL%iESOIL.EQ.iopt_rootweight) DPARAM%RTFRAC2 = 1._SP - MPARAM%RTFRAC1
! ---------------------------------------------------------------------------------------
END SUBROUTINE PAR_DERIVE
