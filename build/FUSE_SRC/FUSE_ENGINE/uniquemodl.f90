SUBROUTINE UNIQUEMODL(NMOD)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007; modified in 2008 to include rainfall errors
! Modified by Brian Henn to include snow model, 6/2013
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Creates an array of character strings that define different model combinations
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE model_defn
!  LIST_*  = lists of options for * different model components
!  AMODL%* = structure that holds all (NMOD) unique combinations
! ---------------------------------------------------------------------------------------
USE nrtype                               
USE model_defn
USE model_defnames
IMPLICIT NONE
! Output
INTEGER(I4B)                           :: NMOD        ! number of model combinations
! Internal
INTEGER(I4B)                           :: ICOUNT      ! loop through unique models
INTEGER(I4B)                           :: ISW_RFERR   ! loop thru rainfall errors
INTEGER(I4B)                           :: ISW_ARCH1   ! loop thru upper layer architecture
INTEGER(I4B)                           :: ISW_ARCH2   ! loop thru lower layer architecture
INTEGER(I4B)                           :: ISW_QSURF   ! loop thru surface runoff
INTEGER(I4B)                           :: ISW_QPERC   ! loop thru percolation
INTEGER(I4B)                           :: ISW_ESOIL   ! loop thru evaporation
INTEGER(I4B)                           :: ISW_QINTF   ! loop thru interflow
INTEGER(I4B)                           :: ISW_Q_TDH   ! loop thru time delay options
INTEGER(I4B)                           :: ISW_SNOWM   ! loop thru snow model options
! Start procedure here
!err=0; message="UNIQUEMODL/ok"
! ---------------------------------------------------------------------------------------
! (1) POPULATE LISTS OF OPTIONS FOR THE DIFFERENT MODEL COMPONENTS
! ---------------------------------------------------------------------------------------
! rainfall error
LIST_RFERR(1)%MCOMPONENT = 'additive_e' ! additive rainfall error
LIST_RFERR(2)%MCOMPONENT = 'multiplc_e' ! multiplicative rainfall error
! upper-layer architecture
LIST_ARCH1(1)%MCOMPONENT = 'tension1_1' ! upper layer broken up into tension and free storage
LIST_ARCH1(2)%MCOMPONENT = 'tension2_1' ! tension storage sub-divided into recharge and excess
LIST_ARCH1(3)%MCOMPONENT = 'onestate_1' ! upper layer defined by a single state variable
! lower-layer architecture -- defines method for computing baseflow
LIST_ARCH2(1)%MCOMPONENT = 'tens2pll_2' ! tension reservoir plus two parallel tanks
LIST_ARCH2(2)%MCOMPONENT = 'unlimfrc_2' ! baseflow resvr of unlimited size (0-HUGE), frac rate
LIST_ARCH2(3)%MCOMPONENT = 'unlimpow_2' ! baseflow resvr of unlimited size (0-HUGE), power recession
LIST_ARCH2(4)%MCOMPONENT = 'fixedsiz_2' ! baseflow reservoir of fixed size
! surface runoff
LIST_QSURF(1)%MCOMPONENT = 'arno_x_vic' ! ARNO/Xzang/VIC parameterization (upper zone control)
LIST_QSURF(2)%MCOMPONENT = 'prms_varnt' ! PRMS variant (fraction of upper tension storage)
LIST_QSURF(3)%MCOMPONENT = 'tmdl_param' ! TOPMODEL parameterization (only valid for TOPMODEL qb)
! percolation
LIST_QPERC(1)%MCOMPONENT = 'perc_f2sat' ! water from (field cap to sat) avail for percolation
LIST_QPERC(2)%MCOMPONENT = 'perc_w2sat' ! water from (wilt pt to sat) avail for percolation
LIST_QPERC(3)%MCOMPONENT = 'perc_lower' ! perc defined by moisture content in lower layer (SAC)
! evaporation fluxes (lower layer evap = 0 for ['tension2_1','unlimfrc_2','unlimpow_2','topmdexp_2']
LIST_ESOIL(1)%MCOMPONENT = 'sequential' ! sequential evaporation model
LIST_ESOIL(2)%MCOMPONENT = 'rootweight' ! root weighting
! interflow
LIST_QINTF(1)%MCOMPONENT = 'intflwnone' ! no interflow
LIST_QINTF(2)%MCOMPONENT = 'intflwsome' ! interflow
! time delay in runoff
LIST_Q_TDH(1)%MCOMPONENT = 'rout_gamma' ! use a Gamma distribution with shape parameter = 2.5
LIST_Q_TDH(2)%MCOMPONENT = 'no_routing' ! no routing
! snow model switch
LIST_SNOWM(1)%MCOMPONENT = 'no_snowmod' ! no snow model
LIST_SNOWM(2)%MCOMPONENT = 'temp_index' ! temperature index snow model
! ---------------------------------------------------------------------------------------
! (2) LOOP THROUGH MODEL COMPONENTS AND DEFINE A SET OF UNIQUE MODELS
! ---------------------------------------------------------------------------------------
! sequence of model-building decisions
! a) define rainfall error
! b) define upper-layer architecture
! c) define lower-layer architecture
! d) define surface runoff method
! e) define percolation method
! f) define evaporation method
! g) define interflow method
! h) define time delay in runoff
ICOUNT = 0 ! initialize counter
! loop through snow model options
DO ISW_SNOWM=1,SIZE(LIST_SNOWM)
! (loop through time delay options)
DO ISW_Q_TDH=1,SIZE(LIST_Q_TDH)
 ! (loop through interflow options)
 DO ISW_QINTF=1,SIZE(LIST_QINTF)
  ! (loop through evaporation options)
  DO ISW_ESOIL=1,SIZE(LIST_ESOIL)
   ! (loop through percolation options)
   DO ISW_QPERC=1,SIZE(LIST_QPERC)
    ! (loop through surface runoff options)
    DO ISW_QSURF=1,SIZE(LIST_QSURF)
     ! (loop through lower-layer architecture options)
     DO ISW_ARCH2=1,SIZE(LIST_ARCH2)
      ! (loop through upper-layer architecture options)
      DO ISW_ARCH1=1,SIZE(LIST_ARCH1)
       ! (loop through rainfall error options)
       DO ISW_RFERR=1,SIZE(LIST_RFERR)
        ! don't allow a lower tension tank when there are two upper ones
        IF (LIST_ARCH1(ISW_ARCH1)%MCOMPONENT(1:10).EQ.'tension2_1'.AND. &
            LIST_ARCH2(ISW_ARCH2)%MCOMPONENT(1:10).EQ.'tens2pll_2') CYCLE
        ! don't allow percolation below field capacity if there are multiple upper tanks
        IF (LIST_ARCH1(ISW_ARCH1)%MCOMPONENT(1:10).NE.'onestate_1'.AND. &
            LIST_QPERC(ISW_QPERC)%MCOMPONENT(1:10).EQ.'perc_w2sat') CYCLE
        ICOUNT = ICOUNT + 1  ! (increment counter)
        IF (ICOUNT.LE.SIZE(AMODL)) THEN
         ! save unique model combinations
         AMODL(ICOUNT)%iRFERR = desc_str2int(LIST_RFERR(ISW_RFERR)%MCOMPONENT)
         AMODL(ICOUNT)%iARCH1 = desc_str2int(LIST_ARCH1(ISW_ARCH1)%MCOMPONENT)
         AMODL(ICOUNT)%iARCH2 = desc_str2int(LIST_ARCH2(ISW_ARCH2)%MCOMPONENT)
         AMODL(ICOUNT)%iQSURF = desc_str2int(LIST_QSURF(ISW_QSURF)%MCOMPONENT)
         AMODL(ICOUNT)%iQPERC = desc_str2int(LIST_QPERC(ISW_QPERC)%MCOMPONENT)
         AMODL(ICOUNT)%iESOIL = desc_str2int(LIST_ESOIL(ISW_ESOIL)%MCOMPONENT)
         AMODL(ICOUNT)%iQINTF = desc_str2int(LIST_QINTF(ISW_QINTF)%MCOMPONENT)
         AMODL(ICOUNT)%iQ_TDH = desc_str2int(LIST_Q_TDH(ISW_Q_TDH)%MCOMPONENT)
         AMODL(ICOUNT)%iSNOWM = desc_str2int(LIST_Q_TDH(ISW_SNOWM)%MCOMPONENT)
         !write(*,'(i3,1x,7(a10,1x))') icount, amodl(icount)
        ELSE
         ! need to allocate more space
         print *, 'insufficent space to hold model combinations'
         stop
        ENDIF
       END DO  ! RFERR
      END DO  ! ARCH1
     END DO  ! ARCH2
    END DO  ! QSURF
   END DO  ! QPERC
  END DO  ! ESOIL
 END DO  ! QINTF
END DO  ! Q_TDH
END DO ! SNOWM
! ---------------------------------------------------------------------------------------
NMOD = ICOUNT
!pause
END SUBROUTINE UNIQUEMODL
