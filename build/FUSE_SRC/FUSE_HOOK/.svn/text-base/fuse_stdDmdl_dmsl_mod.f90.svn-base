!******************************************************************
! (C) Copyright 2009-2010  ---  Dmitri Kavetski and Martyn Clark ---  All rights reserved
!******************************************************************
module fuse_stdDmdl_dmsl_mod
! Purpose: Standard dynamic model template for FUSE.
use kinds_dmsl_kit_FUSE
use model_defn,only:FUSE_version,FUSE_enabled
implicit none
!----------------------------------------------------
private
public::FUSE_version,FUSE_enabled
public::FUSE_setModel,FUSE_getModelInfo,FUSE_cebarModel,FUSE_controlModel
public::FUSE_runModel,FUSE_runAllModel
!----------------------------------------------------
! * Basic properties: numbers of parameters and states
character(*),parameter::modelNameFUSE="FUSE_"
character(*),parameter::indxNameFUSE="time"
integer(mik),parameter::nInputFUSE=2,nOutputFUSE=1
integer(mik),parameter::parTranDefFUSE=0 ! default parameter transformations !DK: needs to be read from file
!----------------------------------------------------
contains
!-----------------------------------------------------------------------------------------
! ***** SET MODEL ******************************************************************
!-----------------------------------------------------------------------------------------
subroutine FUSE_setModel(modelID,setupCmd,chvarLibDef,err,message)
! Purpose: get setup information for the FUSE model
! At this stage, model parameters or even their number are not known by BATEA.
! This routine obtains the FUSE configuration from file.
USE model_defn,only:nstateFUSE=>nstate   ! defines the set of FUSE models
USE metaoutput,only:vardescribe          ! defines output for the FUSE models
! informational modules
use fuse_fileManager,only:fuse_SetDirsUndPhiles
USE selectmodl_module,only:selectmodl    ! identify the model using a control file
use model_numerix,only:JAC_RECOMPUTE,CONSTFULLSTEP,FJACCOPY,FJACDCMP,FJACINDX
! Purpose: get setup information for the FUSE model
implicit none
! dummies
integer(mik),intent(in)::modelID(:)
character(*),intent(in)::setupCmd
character(*),intent(in),optional::chvarLibDef(:,:)
integer(mik),intent(out)::err
character(*),intent(out)::message
! local variables
integer(mik)::nmod
character(5)::jchar
integer(mik)::lenJCH=len(jchar)
! Start procedure here
err=0; message="ok"
! check that the file exists
if(setupCmd/=" ")then
  jchar=setupCmd(:lenJCH) ! determine if musterfile or filemanager supplied
  selectcase(jchar)
  case("[fmf]","[FMF]") ! file manager file supplied
    call fuse_SetDirsUndPhiles(fuseFileManagerIn=setupCmd(lenJCH+1:),err=err,message=message)
  case("[mdf]","[MDF]") ! muster direktor file supplied
    call fuse_SetDirsUndPhiles(fuseMusterDirektorIn=setupCmd(lenJCH+1:),err=err,message=message)
  case default
    call fuse_SetDirsUndPhiles(fuseFileManagerIn=setupCmd,err=err,message=message)
  endselect
else
  call fuse_SetDirsUndPhiles(err=err,message=message)
endif
if(err>0)then ! somethign actually went wrong
  message="f-FUSE_setModel/&"//trim(message)
  err=100; return
else          ! just use default file (not a problem)
  err=0
endif
! Define model attributes (valid for all models)
call uniquemodl(nmod,err,message)                 ! get nmod unique models
if(err/=0)then
  message="f-FUSE_setModel/&"//trim(message)
  err=100; return
endif
call vardescribe()                                ! model variable descriptions (store in module metaoutput)
call getnumerix(err,message)                      ! decisions/parameters that define the numerical scheme
if(err/=0)then
  message="f-FUSE_setModel/&"//trim(message)
  err=100; return
endif
call getparmeta(err,message)                      ! read parameter metadata (parameter bounds, etc.) 
if(err/=0)then
  message="f-FUSE_setModel/&"//trim(message)
  err=100; return
endif
! Identify a single model (read control file)
call selectmodl(istatus=err,message=message)
if(err/=0)then
  message="f-FUSE_setModel/&"//trim(message)
  err=200; return
endif
!write(*,*) LEN_TRIM(SMODL%MNAME), ' - ', TRIM(SMODL%MNAME)
! determine number of states
call assign_stt()         ! state definitions stored in module model_defn [nstateFUSE]
! determine number of parameters
call assign_par()         ! parameter defintions stored in module multiparam [nparFUSE]
! Allocate Jacobian if necessary
IF (JAC_RECOMPUTE.EQ.CONSTFULLSTEP) THEN
 ALLOCATE(fjacCOPY(nstateFUSE,nstateFUSE),fjacDCMP(nstateFUSE,nstateFUSE),fjacINDX(nstateFUSE))
ENDIF
! End procedure here
endsubroutine FUSE_setModel
!-----------------------------------------------------------------------------------------
! ***** GET MODEL INFO *******************************************************************
!-----------------------------------------------------------------------------------------
subroutine FUSE_getModelInfo(modelID,infoCmd,&
  modelName,ninput,nstate,npar,&
  indxName,inputName,stateName,parName,&
  stateLo,stateHi,parLo,parHi,inScal,stateScal,parScal,&
  stateDef,parDef,parSD,parTranDef,parFitDef,&
  err,message)
! Purpose: Returns basic properties of the FUSE model
! data modules
USE model_defn,only:SMODL,CSTATE,nstateFUSE=>nstate   ! defines the set of FUSE models
USE model_defnames,only:desc_int2str
USE multiparam,only:paratt,lparam,numpar              ! parameter attribute structure
USE multistate,only:fstate,dstate,fracstate0          ! defines the states for the FUSE models
USE multiforce,only:DELTIM                            ! model time step (days)
USE metaoutput,only:VNAME,NOUTVAR         ! defines output for the FUSE models
! informational modules
USE str_2_xtry_module,only:str_2_xtry      ! gets state vector from structure in multistate
USE getpar_str_module,only:getpar_str      ! gets parameter metadata structure
USE par_insert_module,only:par_insert      ! puts specific parameter into structure in multiparam
USE parextract_module,only:get_parset      ! gets specific parameter from structure in multiparam
implicit none
! dummies
integer(mik),intent(in)::modelID(:)
character(*),intent(in)::infoCmd
character(*),intent(out),optional::modelName
integer(mik),intent(out),optional::ninput,nstate,npar
character(*),intent(out),optional::indxName ! this variable appeared in BATEAU v 502
character(*),intent(out),dimension(:),optional::inputName,stateName,parName
real(mrk),intent(out),dimension(:),optional::stateLo,stateHi,parLo,parHi,&
  inScal,stateScal,parScal,stateDef,parDef,parSD
integer(mik),intent(out),optional::parTranDef(:)
logical(mlk),intent(out),optional::parFitDef(:)
integer(mik),intent(out)::err
character(*),intent(out)::message
! local variables
! character(len=2) :: cnum               ! model number as text
! integer(mik)     :: nmod               ! number of unique models
integer(mik)     :: i !,j,k              ! looping variables
integer(mik)     :: istart             ! start index of variable list (to define output)
real(mrk)        :: frac               ! fraction of capacity to initialize states
type(paratt)     :: param_meta         ! parameter metadata
real(mrk)::dt !DK_BOTCH: hardcode 'dt'
! Start procedure here
err=0; message="ok"; dt=1._mrk
! Define model data step
DELTIM = dt
! Define model name
if(present(modelName)) modelName = smodl%mname ! smodl is in module model_defn  
! define model inputs (assume inputs are the ***first*** nInputs in varlist)
if(present(ninput))   nInput=nInputFUSE
if(present(indxName)) indxName=indxNameFUSE
if(present(inputName))forall(i=1:nInputFUSE) inputName(i) = vname(i)
! define model states
if(present(nstate))then
 nstate=nstateFUSE+nOutputFUSE ! +nOutputFUSE to include model outputs in "state" list
endif
! define model outputs (assume outputs are the ***last*** nOutputs in varlist)
if(present(stateName))then
  istart = (noutvar-nOutputFUSE)+1
  stateName(1:nOutputFUSE) = vname(istart:noutvar)
  stateName(nOutputFUSE+1:nOutputFUSE+nstateFUSE) = desc_int2str(cstate%isname)
endif
! define model parameters
if(present(npar))then
 npar=numpar               ! numpar from module multiparam
endif
if(present(parName))       forall(i=1:numpar) parName(i) = lparam(i)%parname
! define parameter ranges and default transformations
if(present(parLo) .and. present(parHi).and.present(parTranDef)) then
 do i=1,numpar
  call getpar_str(lparam(i)%parname,param_meta)
  parLo(i)   = param_meta%parlow
  parHi(i)   = param_meta%parupp
  parTranDef = param_meta%parvtn
 end do
endif
! define state ranges
if(present(stateLo) .and. present(stateHi)) then
 stateLo = 0._mrk                        ! set minimum states to zero
 ! (use the default parameter values to set bucket sizes)
 do i=1,numpar
  call getpar_str(lparam(i)%parname,param_meta)        ! extract full metadata structure
  call par_insert(param_meta%pardef,lparam(i)%parname) ! insert the default param to model param structure
 end do
 call par_derive(err,message)                       ! identify the derived parameters associated with mparam
 if(err/=0)then
  message="f-FUSE_getModelInfo/&"//trim(message); return
 endif
 frac = 1._mrk; call init_state(frac)    ! initialize states at fraction (frac) of capacity
 call str_2_xtry(fstate,stateHi)         ! extract a vector of states at the maximum value
endif
! define scaling factors
if(present(inScal))   inScal(1:nInputFUSE)    = 10._mrk
if(present(stateScal))stateScal(1:nstateFUSE) = 10._mrk
if(present(parScal))  parScal(1:numpar)       = 10._mrk
! define default parameter values
if(present(stateDef)) then
 ! (use the default parameter values to set default states)
 do i=1,numpar
  call getpar_str(lparam(i)%parname,param_meta)        ! extract full metadata structure
  call par_insert(param_meta%pardef,lparam(i)%parname) ! insert the default param to model param structure
 end do
 call par_derive(err,message)                       ! identify the derived parameters associated with mparam
 if(err/=0)then
  message="f-FUSE_getModelInfo/&"//trim(message); return
 endif
 call init_state(fracState0)             ! initialize states at fraction (frac) of capacity
 call str_2_xtry(fstate,stateDef)        ! extract a vector of states at the value tstate
 dstate=fstate                           ! save default states in module multistate
endif
if(present(parDef)) then
 do i=1,numpar
  call getpar_str(lparam(i)%parname,param_meta)  ! extract full metadata structure
  parDef(i)=param_meta%pardef                    ! set default param to value from structure
 end do
endif
if(present(parSD))parSD=undefRN
if(present(parFitDef))parFitDef=.true.
! End procedure here
endsubroutine FUSE_getModelInfo
!-----------------------------------------------------------------------------------------
! ***** PRIME MODEL WITH INPUT FILES, ETC  ******************************************************************
!-----------------------------------------------------------------------------------------
subroutine FUSE_cebarModel(modelID,cebarCmd,dataXY,dataProps,err,message)
! Purpose: This routine is used to prime the model for execution.
implicit none
! dummies
integer(mik),intent(in)::modelID(:)
character(*),intent(in)::cebarCmd
real(mrk),intent(in)::dataXY(:,:)
real(mrk),intent(in)::dataProps(:)
integer(mik),intent(out)::err
character(*),intent(out)::message
! local variables
! Start procedure here
err=0; message="ok"
! check that the file exists
if(cebarCmd/=" ")then
  message="w-FUSE_cebarModel/doesntUseFiles/&"//"[file'"//trim(cebarCmd)//"'notUsed]"
  err=100; return
endif
! End procedure here
endsubroutine FUSE_cebarModel
!-----------------------------------------------------------------------------------------
! ***** GET MODEL CONTROL ****************************************************************
!-----------------------------------------------------------------------------------------
subroutine FUSE_controlModel(modelID,inittCmd,dataXY,dataProps,parIn,dquanIn,&
  parOut,flexSin,setS0in,stateIn,stateOut,feas,err,message)
! Purpose: Sets/Gets model states and parameters.
! Usage:
!  - if(setS0in) then will  set all states to default values
!                    this is convenient when initialising the model without calibrating S0.
!  - if(flexSin) then will adjust states to be compatible with parameter values,
!                     eg, if state S exceeds its maximum value Smax, will reset S to Smax.
! data modules
USE multistate,only:fstate,mstate,fracstate0,hstate                      ! defines the states for the FUSE models
use multiforce,only:deltim
use multiparam,only:lparam
use multiroute,only:mroute
USE metaoutput,only:VNAME,NOUTVAR         ! defines output for the FUSE models
! informational modules
USE par_insert_module,only:par_insert,put_parset    ! puts specific parameter into structure in multiparam
USE parextract_module,only:get_parset               ! gets specific parameter from structure in multiparam
USE xtry_2_str_module,only:xtry_2_str               ! puts state vector into structure in multistate
USE str_2_xtry_module,only:str_2_xtry               ! gets state vector from structure in multistate
! DMSL
!use utilities_dmsl_kit,only:quickif
implicit none
! dummies
integer(mik),intent(in)::modelID(:)
character(*),intent(in),optional::inittCmd
real(mrk),intent(in),optional::dataXY(:,:),dataProps(:)
real(mrk),intent(in),optional::parIn(:),dquanIn(:),stateIn(:)
logical(mlk),intent(in),optional::flexSin,setS0in
real(mrk),intent(out),optional::parOut(:),stateOut(:)
logical(mlk),intent(out),optional::feas
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
logical(mlk)::haveSin,newRun !flexS,setS0,checkFeas
logical(mlk),parameter::flexSdef=.false.,setS0def=.false.
character(200)::parName1
integer(mik)::i,istart
real(mrk)::dt
! Start procedure here
! (a) Set FUSE parameters
err=0; message="FUSE_controlModel/ok"; if(present(feas))feas=.true.; haveSin=present(stateIn)
if(present(dataProps))then;  dt=dataProps(1)
else;                        dt=undefRN; endif
newRun=.false.;if(present(setS0in))newRun=setS0in
! newRun=quickif(setS0in,.false.) ! flag to avoid recomputing derived parameters:
!DK_NB: This is not a general fix, because prevents stochastic parameters other than rMult
! (b) Put/Get parameters into and out of the model structure
if (present(parIn)) then
  if(newRun)then  ! this happens at beginning of each new run
    call put_parset(parIn)  ! base parameters
    call par_derive(err,message)       ! corresponding derived parameters
    if(err/=0)then
      message="f-FUSE_controlModel/&"//trim(message); return
    endif
  else            ! this does stochastic parameters. currently only rain-error can be stochastic
    parName1=lparam(1)%parname
    select case (trim(parName1))
    case('RFERR_ADD','RFERR_MLT')
      call par_insert(parIn(1),parName1)
    case default
      err=100; message="f-FUSE_controlModel/unsupportedStochPar["//trim(parName1)//"]"
      return
    endselect
  endif
endif
if (present(parOut)) call get_parset(parOut)
! (c) Put/Get states into and out of the model structure
if (present(stateIn)) then
 call xtry_2_str(stateIn,fstate)  ! populates fstate
 mstate = fstate                  ! initialize the model state
endif
if (present(stateOut))then
 istart = (noutvar-nOutputFUSE)+1
 do i=istart,noutvar  ! noutvar is in module metaoutput
  if (trim(vname(i))=='q_routed')stateOut((i-istart)+1) = mroute%q_routed
 enddo
 call str_2_xtry(fstate,stateOut(nOutputFUSE+1:))
endif
! (d) Adjust states to be compatible w/ param values
if (present(flexSin)) then ! (needed for the case of stochastic parameters)
 if (flexSin) call adjust_stt()
endif
! (e) re-initialize states to default values
if (present(setS0in)) then ! (convenient when initialising w/o calibrating S0)
  if (setS0in)then
    call init_state(fracState0) ! initialize states at fraction (frac) of capacity
    !dstate=fstate ! save the initial state as the default state (not needed) MPC 2009/10/09
    mstate=fstate  ! save the initial state as the model state 2009/10/09
    hstate%step = dt ! deltim is shared in module multiforce.
!DK: NB: dt should NOT change between controlModel and runModel, must equal deltime.
! hstate%step is reduced/increased by error control. dont reset it between time steps.
  endif
endif
! (f) check parameter set obeys constraints
! no constrains currently
! End procedure here
endsubroutine FUSE_controlModel
!-----------------------------------------------------------------------------------------
! ***** RUN MODEL ************************************************************************
!-----------------------------------------------------------------------------------------
subroutine FUSE_runModel(modelID,runitCmd,iT,dataProps,input,state,feas,err,message)
! Purpose: Performs single step of FUSE model.
USE model_defn,only:nstate
USE metaoutput,only:vname,noutvar          ! defines output for the FUSE models
USE multiforce,only:mforce,deltim          ! forcing structure
USE multiroute,only:mroute          ! routing structure
USE multistate,only:hstate,fstate
USE interfaceb,only:ode_int,fuse_solve
USE str_2_xtry_module,only:STR_2_XTRY ! gets state vector from structure in multistate
implicit none
! dummies
integer(mik),intent(in)::modelID(:)
character(*),intent(in)::runitCmd
integer(mik),intent(in)::iT
real(mrk),intent(in)::dataProps(:)
real(mrk),intent(in)::input(:)
real(mrk),intent(out)::state(:)
logical(mlk),intent(out)::feas
integer(mik),intent(out)::err
character(*),intent(out)::message
! local
integer(mik)::i
real(mrk)::dt_sub,dt_full
real(mrk),dimension(nstate)::state0,state1
! Start procedure here
err=0; message="ok"; feas=.true.; state=undefRN
!DK_NB: hstate%step is reduced/increased by error control. dont reset it between time steps.
! get model inputs and put them in the structure
do i=1,nInputFUSE   ! (assume the first in the variable name list)
  if (trim(vname(i))=='ppt') mforce%ppt = input(i)
  if (trim(vname(i))=='pet') mforce%pet = input(i)
end do  ! (loop thru inputs)
DT_FULL = DELTIM
DT_SUB  = HSTATE%STEP
CALL STR_2_XTRY(FSTATE,STATE0)  ! get the vector of states from the FSTATE structure
CALL INITFLUXES()               ! set weighted sum of fluxes to zero
! temporally integrate the ordinary differential equations
CALL ODE_INT(FUSE_SOLVE,STATE0,STATE1,DT_SUB,DT_FULL,ERR,MESSAGE)
IF (ERR/=0) THEN
  message="f-FUSE_runModel/&"//TRIM(MESSAGE); return
ENDIF
HSTATE%STEP = DT_SUB
! perform overland flow routing
CALL Q_OVERLAND()
! get model outputs (assume the last in the variable name list)
call FUSE_controlModel(modelID,stateOut=state,err=err,message=message)
! compute summary statistics
CALL COMP_STATS()
! End procedure here
endsubroutine FUSE_runModel
!-----------------------------------------------------------------------------------------
! ***** RUNALL MODEL ************************************************************************
!-----------------------------------------------------------------------------------------
subroutine FUSE_runAllModel(modelID,runallCmd,dataProps,input,state,feas,err,message)
! Purpose: Performs all steps of FUSE model.
implicit none
! dummies
integer(mik),intent(in)::modelID(:)
character(*),intent(in)::runallCmd
real(mrk),intent(in)::dataProps(:)
real(mrk),intent(in)::input(:,:)
real(mrk),intent(out)::state(:,:)
logical(mlk),intent(out)::feas
integer(mik),intent(out)::err
character(*),intent(out)::message
! local
integer(mik)::nT,iT
! Start procedure here
nT=size(input,1)
do iT=1,nT
  call FUSE_runModel(modelID,runallCmd,iT,dataProps,&
    input(iT,:),state(iT,:),feas,err,message)
  if(err/=0)then
    write(message,'(a,i0,a)')"f-FUSE_runAllModel/[iT=",iT,"]/&"//trim(message)
    err=20; return
  endif
enddo
! End procedure here
endsubroutine FUSE_runAllModel
!----------------------------------------------------
endmodule fuse_stdDmdl_dmsl_mod
!******************************************************************
