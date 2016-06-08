program make_batea_parfiles
! Martyn Clark, 2009
! used to make parameter files for BATEA
use nrtype                              ! variable types
use selectmodl_module                   ! access to SUBROUTINE selectmodl
implicit none
integer(i4b)                 :: nmod    ! number of possible models
integer(i4b)                 :: ierr    ! error code
integer(i1b)                 :: ipar    ! looping
character(len=256)           :: message ! error message
! ----------------------------------------------------------------------------------------
! get parameter metadata for all possible models
call getparmeta()
! identify the model used
call uniquemodl(nmod)         ! get nmod unique models
call selectmodl(ierr,message) ! identify single model (read control file m_decisions.txt)
if (ierr.ne.0) then; print *, trim(message); stop; endif
! identify the parameters used in the model selected
call assign_par()             ! parameters used are stored in module multiparam
! write parameter file for batea
call batea_file()
stop
end program make_batea_parfiles
