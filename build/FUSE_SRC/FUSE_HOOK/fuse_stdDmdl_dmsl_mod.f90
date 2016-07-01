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
integer(mik),parameter::nInputFUSE_base=2,nInputFUSE_snow=4,nOutputFUSE=1
integer(mik),save::nInputFUSE=undefIN
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
USE model_defn,only:nstateFUSE=>nstate,SMODL   ! defines the set of FUSE models
USE model_defnames,only:iopt_temp_index        ! defines the integer model options
USE metaoutput,only:vardescribe                ! defines output for the FUSE models
! informational modules
use fuse_fileManager,only:fuse_SetDirsUndPhiles
USE selectmodl_module,only:selectmodl          ! identify the model using a control file
use model_numerix,only:JAC_RECOMPUTE,CONSTFULLSTEP,FJACCOPY,FJACDCMP,FJACINDX
! Purpose: get setup information for the FUSE model
implicit none
! dummies
integer(mik),intent(in)::modelID(:)
character(*),intent(in)::setupCmd
character(*),intent(in),optional::chvarLibDef(:,:)
integer(mik),intent(out)::err
character(*),intent(out)::message
! locals
character(*),parameter::procnam="FUSE_setModel"
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
  message="f-"//procnam//"/&"//message
  err=100; return
else          ! just use default file (not a problem)
  err=0
endif
! Define model attributes (valid for all models)
call uniquemodl(nmod,err,message)                 ! get nmod unique models
if(err/=0)then
  message="f-"//procnam//"/&"//message
  err=100; return
endif
call vardescribe()                                ! model variable descriptions (store in module metaoutput)
call getnumerix(err,message)                      ! decisions/parameters that define the numerical scheme
if(err/=0)then
  message="f-"//procnam//"/&"//message
  err=100; return
endif
call getparmeta(err,message)                      ! read parameter metadata (parameter bounds, etc.) 
if(err/=0)then
  message="f-"//procnam//"/&"//message
  err=100; return
endif
! Identify a single model (read control file)
call selectmodl(err=err,message=message)
if(err/=0)then
  message="f-"//procnam//"/&"//message
  err=200; return
endif
selectcase(SMODL%iSNOWM)
case(iopt_temp_index)
  CALL GET_MBANDS(err,message) ! BH: call of band data in snow model case
  if(err/=0)then
    message="f-"//procnam//"/&"//message
    return
  endif
endselect
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
USE model_defnames,only:desc_int2str,iopt_temp_index
USE multiparam,only:paratt,lparam,numpar              ! parameter attribute structure
USE multistate,only:fstate,dstate,fracstate0          ! defines the states for the FUSE models
USE multiforce,only:DELTIM                            ! model time step (days)
USE MULTIBANDS,only:N_BANDS               ! BH: access snow model bands 
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
! locals
character(*),parameter::procnam="FUSE_getModelInfo"
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
nInputFUSE=nInputFUSE_base
selectcase(SMODL%iSNOWM)
case(iopt_temp_index);nInputFUSE=nInputFUSE+nInputFUSE_snow
endselect
if(present(ninput))   nInput=nInputFUSE
if(present(indxName)) indxName=indxNameFUSE
if(present(inputName))then
  forall(i=1:nInputFUSE_base) inputName(i) = vname(i)
  selectcase(SMODL%iSNOWM)
  case(iopt_temp_index);inputName(nInputFUSE_base+1) = "temp"
  endselect
  if(nInputFUSE>=nInputFUSE_base+1+1)inputName(nInputFUSE_base+1+1)="iy"
  if(nInputFUSE>=nInputFUSE_base+1+2)inputName(nInputFUSE_base+1+2)="im"
  if(nInputFUSE>=nInputFUSE_base+1+3)inputName(nInputFUSE_base+1+3)="id"
endif
! define model states
if(present(nstate))then
 nstate=nstateFUSE+nOutputFUSE+2*N_BANDS ! +nOutputFUSE to include model outputs in "state" list
endif
! define model outputs (assume outputs are the ***last*** nOutputs in varlist)
if(present(stateName))then
  istart = (noutvar-nOutputFUSE)+1
  stateName(1:nOutputFUSE) = vname(istart:noutvar)
  DO i = 1,N_BANDS
    WRITE(stateName(nOutputFUSE+i),'(a,i2.2)') 'swe_z',I
  ENDDO
  DO i = 1,N_BANDS
    WRITE(stateName(nOutputFUSE+N_BANDS+i),'(a,i2.2)') 'temp_z',I
  ENDDO
  stateName((nOutputFUSE+2*N_BANDS+1):(nOutputFUSE+2*N_BANDS+nstateFUSE)) = desc_int2str(cstate%isname)
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
  enddo
  call par_derive(err,message)                       ! identify the derived parameters associated with mparam
  if(err/=0)then
    message="f-"//procnam//"/&"//message; return
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
    message="f-"//procnam//"/&"//message; return
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
! locals
character(*),parameter::procnam="FUSE_cebarModel"
! Start procedure here
err=0; message="ok"
! check that the file exists
if(cebarCmd/=" ")then
  message="w-"//procnam//"/doesntUseFiles/&"//"[cebarCmd='"//trim(cebarCmd)//"'notUsed]"
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
USE model_defn,ONLY:SMODL
USE model_defnames,ONLY:iopt_temp_index
USE multistate,only:fstate,mstate,fracstate0,hstate ! defines the states for the FUSE models
use multiforce,only:deltim,MFORCE
use multiparam,only:MPARAM
USE MULTIBANDS,only:N_BANDS,MBANDS,Z_FORCING        ! added to access snow model SWE, Brian Henn, Oct. 2013
use multiparam,only:lparam
use multiroute,only:mroute
USE metaoutput,only:VNAME,NOUTVAR                   ! defines output for the FUSE models
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
character(*),parameter::procnam="FUSE_controlModel"
logical(mlk)::haveSin,newRun !flexS,setS0,checkFeas
logical(mlk),parameter::flexSdef=.false.,setS0def=.false.
character(200)::parName1
integer(mik)::i,istart
real(mrk)::dt,DZ,TEMP_Z
! Start procedure here
! (a) Set FUSE parameters
err=0; if(present(feas))feas=.true.; haveSin=present(stateIn)
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
      message="f-"//procnam//"/&"//message; return
    endif
  else            ! this does stochastic parameters. currently only rain-error can be stochastic
    parName1=lparam(1)%parname
    selectcase(parName1)
    case('RFERR_ADD','RFERR_MLT')
      call par_insert(parIn(1),parName1)
    case default
      err=100; message="f-"//procnam//"/unsupportedStochPar["//trim(parName1)//"]"
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
    if (vname(i)=='q_routed')stateOut((i-istart)+1) = mroute%q_routed
  enddo
  IF(SMODL%iSNOWM.EQ.iopt_temp_index) then
    DO i = 1,N_BANDS
      stateOut(nOutputFUSE + i) = MBANDS(i)%SWE
    ENDDO
    DO i = 1,N_BANDS
      DZ = MBANDS(I)%Z_MID - Z_FORCING
      TEMP_Z = MFORCE%TEMP + DZ*MPARAM%LAPSE/1000                     
      stateOut(nOutputFUSE + N_BANDS + i) = TEMP_Z
    ENDDO
  ENDIF
  call str_2_xtry(fstate,stateOut((nOutputFUSE + 2*N_BANDS + 1):))
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
USE model_defn,only:nstate,SMODL
USE model_defnames,only:iopt_temp_index,iopt_no_snowmod
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
! locals
character(*),parameter::procnam="FUSE_runModel"
integer(mik)::i
real(mrk)::dt_sub,dt_full
real(mrk),dimension(nstate)::state0,state1
! Start procedure here
err=0; feas=.true.; state=undefRN
!DK_NB: hstate%step is reduced/increased by error control. DO NOT reset it between time steps.
! get model inputs and put them in the structure
do i=1,nInputFUSE   ! (assume the first in the variable name list)
  selectcase (vname(i))
  case('ppt');  mforce%ppt  = input(i)
  case('pet');  mforce%pet  = input(i)
  case('temp'); mforce%temp = input(i) ! temperature added to run snow model 
!   case default
!     write(message,'(a,i0,a)')"f-"//procnam//"/unknown[vname(i=",i,")='"//trim(vname(i))//"']"
!     err=100; return
  endselect
enddo  ! the next statements allow extra inputs to be read in
if(nInputFUSE>=nInputFUSE_base+1+1)mforce%iy=input(nInputFUSE_base+1+1)
if(nInputFUSE>=nInputFUSE_base+1+2)mforce%im=input(nInputFUSE_base+1+2)
if(nInputFUSE>=nInputFUSE_base+1+3)mforce%id=input(nInputFUSE_base+1+3)
DT_FULL = DELTIM
DT_SUB  = HSTATE%STEP
CALL STR_2_XTRY(FSTATE,STATE0)  ! get the vector of states from the FSTATE structure
CALL INITFLUXES()               ! set weighted sum of fluxes to zero
! if snow model, call UPDATE_SWE first to calculate snow fluxes and update snow bands 
! using explicit Euler approach; if not, call QRAINERROR
SELECTCASE(SMODL%iSNOWM)
CASE(iopt_temp_index)
  CALL UPDATE_SWE(DELTIM)
CASE(iopt_no_snowmod)
  CALL QRAINERROR()
CASE DEFAULT
  message="f-"//procnam//"/SMODL%iSNOWM must be either iopt_temp_index or iopt_no_snowmod"
  err=100; return
ENDSELECT
! temporally integrate the ordinary differential equations
CALL ODE_INT(FUSE_SOLVE,STATE0,STATE1,DT_SUB,DT_FULL,ERR,MESSAGE)
IF (ERR/=0) THEN
  message="f-"//procnam//"/&"//message; return
ENDIF
HSTATE%STEP = DT_SUB
! perform overland flow routing
CALL Q_OVERLAND()
! get model outputs (assume the last in the variable name list)
call FUSE_controlModel(modelID,stateOut=state,err=err,message=message)
IF(ERR/=0)THEN
  message="f-"//procnam//"/&"//message; return
ENDIF
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
! locals
character(*),parameter::procnam="FUSE_runAllModel"
integer(mik)::nT,iT
! Start procedure here
nT=size(input,1)
do iT=1,nT
  call FUSE_runModel(modelID,runallCmd,iT,dataProps,&
                     input(iT,:),state(iT,:),feas,err,message)
  if(err/=0)then
    write(message,'(a,i0,a)')"f-"//procnam//"/[iT=",iT,"]/&"//trim(message)
    err=20; return
  endif
enddo
! End procedure here
endsubroutine FUSE_runAllModel
!----------------------------------------------------
endmodule fuse_stdDmdl_dmsl_mod
!******************************************************************
