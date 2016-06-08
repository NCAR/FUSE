program batea_test
! BATEA modules
use ddirectory                      ! define directory that holds data
use kinds_dmsl_kit_FUSE             ! define data types
use fuse_stdDmdl_dmsl_mod           ! linking routines for FUSE
! FUSE data modules
use multiforce                      ! defines model forcing data
use multistate                      ! defines the states for the FUSE models
USE multiparam, only:paratt,lparam  ! parameter attribute structure
! FUSE informational modules
use str_2_xtry_module               ! gets state vector from structure in multistate
USE getpar_str_module               ! gets parameter metadata structure
USE par_insert_module               ! puts specific parameter into structure in multiparam
use parextract_module               ! gets specific parameter from structure in multiparam
implicit none
! general variables
integer(mik)                                :: modelID
integer(mik)                                :: err
character(len=256)                          :: message
! model info
character(len=256)                          :: modelName
integer(mik)                                :: ninputs,noutputs,nstate,npar,ninfo
character(len=256),dimension(:),allocatable :: inputName,outputName,stateName,parName,infoStateName
real(mrk),dimension(:),allocatable          :: parLo,parHi,stateLo,stateHi,&       ! ranges
                                               parScal,stateScal,inScal,outScal,&  ! scaling factors
                                               parDef,stateDef                     ! defaults 
integer(mik),dimension(:),allocatable       :: parTranDef                          ! param transform code
type(paratt)                                :: param_meta                          ! parameter metadata
! model control
integer(mik)                                :: iparset    ! case for type of parameter set
integer(mik), parameter                     :: irandom=0  ! random parameter set
integer(mik), parameter                     :: idefault=1 ! default parameter set
real(mrk),dimension(:),allocatable          :: parIn,parOut,stateIn,stateOut
logical(mlk)                                :: feas,setS0in,flexSin
real(mrk)                                   :: frac       ! used to provide an example state vector
! model run
character(len=8)                            :: cbasid          ! basin ID
integer(mik)                                :: itim,ntim       ! loop through time
real(mrk),dimension(:),allocatable          :: input,output,infoState
! local variables
integer(mik)                                :: i !,j,k           ! looping variables
integer(mik)                                :: ierr(10) !,icheck ! status codes for allocate statement
! real(mrk)                                   :: frac
real(mrk)                                   :: tA,tB
! ---------------------------------------------------------------------------------------
! (0) DEFINE PATH NAMES AND READ FORCING DATA
! ---------------------------------------------------------------------------------------
! Read data from the "BATEA-compliant" ASCII files
call getforcing(ntim)
! ---------------------------------------------------------------------------------------
! (1) SETUP AND MODEL INFO
! ---------------------------------------------------------------------------------------
modelID=-9999
! (a) get configuration and dimensions
call FUSE_getModelInfo(modelID,&
                       modelName,ninputs,noutputs,nstate,npar,ninfo,&
                       err=err,message=message)
! ---------------------------------------------------------------------------------------
! (b) allocate space
allocate(inputName(ninputs),outputName(noutputs),stateName(nstate),parName(npar),&
         infoStateName(ninfo), stat=ierr(1))
allocate(parLo(npar),parHi(npar), stat=ierr(2))
allocate(stateLo(nstate),stateHi(nstate), stat=ierr(3))
allocate(parScal(npar),stateScal(npar),inScal(npar),outScal(npar), stat=ierr(4))
allocate(parDef(npar),stateDef(nstate),parTranDef(npar), stat=ierr(5))
if (any(ierr(1:4).ne.0)) stop ' problem allocating space for model info '
write(*,*) len(modelName),len_trim(modelName),modelName(1:len_trim(modelName))
write(*,*) ninputs,noutputs,nstate,npar,ninfo
! (c) get model info
call FUSE_getModelInfo(modelID,&
                       modelName,ninputs,noutputs,nstate,npar,ninfo,&
                       inputName,outputName,stateName,parName,infoStateName,&
                       parLo,parHi,stateLo,stateHi,&
                       parScal,stateScal,inScal,outScal,&
                       parDef,stateDef,parTranDef,&
                       err=err,message=message)
write(*,*) 'after FUSE_getModelInfo'
!----------------------------------------------------------------------------------------
! (2) PRIME THE MODEL (TOPO DATA, ETC)
call FUSE_CebarModel(modelID,deltim,err=err,message=message)
if (err.ne.0) then
 write(*,*) trim(message)
 stop
endif
write(*,*) ' after FUSE_GetModelSetup '
! ---------------------------------------------------------------------------------------
! (3) GET MODEL CONTROL
! ---------------------------------------------------------------------------------------
! (a) allocate space
allocate(parIn(npar),parOut(npar),stateIn(nstate),stateOut(nstate), stat=ierr(1))
if (ierr(1).ne.0) stop ' problem allocating space for model control '
! (b) get an example model parameter set
! switch between random parameter set
!iparset  = 1  ! irandom = 0; idefault = 1
!select case(iparset)
! random parameter set
!case(irandom)
 !call get_params(1)       ! fill structure APARAM with just one parameter set
 !mparam=aparam(1)         ! set current parameter set to the parameter set just extracted
 !do i=1,npar; parIn(i) = parextract(parName(i)); end do  ! (extract parameters from mparam)
! default parameter set
!case(idefault)
 ! (use the default parameter values to set default states)
 do i=1,npar
  call getpar_str(lparam(i)%parname,param_meta)        ! extract full metadata structure
  call par_insert(param_meta%pardef,lparam(i)%parname) ! insert the default param to model param structure
  parIn(i) = param_meta%pardef
 end do
!case default
! write(*,*) 'case iparset must be either ', irandom, ' or ', idefault
! stop
!end select
! (c) get an example set of model states for that parameter set
call par_derive()         ! identify the derived parameters associated with mparam
frac = 0.5_mrk            ! define the fraction of capacity to initialize states
call init_state(frac)     ! initialize states at fraction (frac) of capacity
tstate=fstate             ! set current state to the first state
call str_2_xtry(stateIn)  ! extract a vector of states at the value tstate
! (d) define input flags
flexSin = .true.          ! (.true. = adjust states to be compatible w/ param values)
! setS0in = .false.        ! (.true. = states are re-initialized to default values)
setS0in = .true.          ! (.true. = states are re-initialized to default values)
! (e) call model control
call FUSE_controlModel(modelID,deltim,parIn,parOut,flexSin,setS0in,stateIn,stateOut,feas,&
                       err,message)
do i=1,ninputs
  write(*,*) i, trim(inputName(i))
end do
write(*,*) '----------'
do i=1,noutputs
  write(*,*) i, trim(outputName(i))
end do
write(*,*) '----------'
do i=1,nstate
  write(*,'(i2,1x,a9,1x,3(f9.3,1x))') i, stateName(i), stateDef(i), stateLo(i), stateHi(i)
end do
write(*,*) '----------'
do i=1,npar
  write(*,'(i2,1x,a9,1x,3(f9.3,1x))') i, parName(i), parIn(i), parLo(i), parHi(i)
end do
write(*,*) '----------'
do i=1,ninfo
  write(*,*) i, len(infoStateName(i)), len_trim(infoStateName(i)), trim(infoStateName(i))
end do
write(*,*) '----------'
pause
! ---------------------------------------------------------------------------------------
! (4) RUN MODEL
! ---------------------------------------------------------------------------------------
open(21,file=ModelName(1:8)//'.out',status='unknown')
! (a) allocate space for model inputs and outputs
allocate(input(ninputs),output(noutputs),infoState(ninfo), stat=ierr(1))
if (ierr(1).ne.0) stop ' problem allocating space for model control '
! (b) loop through time
! initialize sub-step length to the length of the time step
! hstate%step = deltim ! deltim is shared in module multiforce
call cpu_time(tA)
do itim=1,ntim
 ! (c) assign model forcing data
 do i=1,ninputs
  if (trim(inputName(i)).eq.'ppt') input(i) = aforce(itim)%ppt
  if (trim(inputName(i)).eq.'pet') input(i) = aforce(itim)%pet
 end do
 ! (d) run model
 call FUSE_runModel(modelID,deltim,input,output,infoState,err,message)
 ! (e) write output
 WRITE( *,'(I10,1X,I4,1X,4(I2,1X),F9.3,1X,F15.1,1X,4(ES12.4,1X))') ITIM, AFORCE(ITIM),OUTPUT
 WRITE(21,'(I10,1X,I4,1X,4(I2,1X),F9.3,1X,F15.1,1X,4(ES12.4,1X))') ITIM, AFORCE(ITIM),OUTPUT
end do  ! (looping through time)
call cpu_time(tB)
write(*,*)"CPU time, sec",tB-tA
close(21)
! ---------------------------------------------------------------------------------------
! deallocate space
deallocate(inputName,outputName,stateName,parName,infoStateName, stat=ierr(1))
deallocate(parLo,parHi,stateLo,stateHi, stat=ierr(2))
deallocate(parScal,stateScal,inScal,outScal, stat=ierr(3))
deallocate(parDef,stateDef, stat=ierr(4))
deallocate(parIn,parOut,stateIn,stateOut, stat=ierr(5))
deallocate(input,output,infoState, stat=ierr(6))
if (any(ierr(1:6).ne.0)) stop ' problem deallocating space '
stop
end program batea_test
