MODULE SELECTMODL_MODULE
IMPLICIT NONE
CONTAINS
! ---------------------------------------------------------------------------------------
SUBROUTINE SELECTMODL(FUSE_ID,ERR,MESSAGE)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2008
! Modified by Brian Henn to include snow model, 6/2013
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Reads a control file and identifies a unique model index
! ---------------------------------------------------------------------------------------
USE nrtype,ONLY:I4B,LGT                               ! defines data types
USE utilities_dmsl_kit_FUSE,ONLY:getSpareUnit,stripTrailString
USE fuse_fileManager,only:SETNGS_PATH,M_DECISIONS     ! defines data directory
USE model_defn,ONLY:NDEC,SMODL,AMODL,&                ! defines model decisions
  LIST_RFERR,LIST_ARCH1,LIST_ARCH2,LIST_QSURF,LIST_QPERC,LIST_ESOIL,&
  LIST_QINTF,LIST_Q_TDH,LIST_SNOWM
USE model_defnames,ONLY:DESC_STR2INT
USE model_numerix,only:solution_method,temporal_error_control
IMPLICIT NONE
! Input
!INTEGER(I4B), INTENT(IN), OPTIONAL     :: FUSE_ID     ! identifier for FUSE model
CHARACTER(LEN=6), INTENT(IN), OPTIONAL :: FUSE_ID     ! identifier for FUSE model
! Output
INTEGER(I4B), INTENT(OUT)              :: ERR         ! error code
CHARACTER(LEN=*), INTENT(OUT)          :: MESSAGE     ! error message
! (1) read in the component choice and model component
LOGICAL(LGT)                           :: READ_FILE   ! .TRUE. if read model decisions from a file
integer(i4b),parameter::lenPath=1024 !DK/2008/10/21: allows longer file paths
INTEGER(I4B)                           :: IUNIT       ! file unit
CHARACTER(LEN=lenPath)                 :: CFILE       ! name of constraints file
LOGICAL(LGT)                           :: LEXIST      ! .TRUE. if file exists
CHARACTER(LEN=256)                     :: KEY         ! format code
INTEGER(I4B)                           :: IDEC        ! loop thru model decisions
CHARACTER(LEN=16)                      :: M_CHOICE    ! model choice
CHARACTER(LEN=8)                       :: DECISION    ! model decision
INTEGER(I4B)                           :: IERR        ! error code for read statement
! (2) loop through model components and identify an index for each decision
INTEGER(I4B)                           :: ICOUNT      ! model indx
INTEGER(I4B)                           :: ISW_MODEL   ! dummy loop
INTEGER(I4B)                           :: ISW_RFERR   ! loop thru rainfall error options
INTEGER(I4B)                           :: ISW_ARCH1   ! loop thru upper layer architecture
INTEGER(I4B)                           :: ISW_ARCH2   ! loop thru lower layer architecture
INTEGER(I4B)                           :: ISW_QSURF   ! loop thru surface runoff
INTEGER(I4B)                           :: ISW_QPERC   ! loop thru percolation
INTEGER(I4B)                           :: ISW_ESOIL   ! loop thru evaporation
INTEGER(I4B)                           :: ISW_QINTF   ! loop thru interflow
INTEGER(I4B)                           :: ISW_Q_TDH   ! loop thru time delay options
INTEGER(I4B)                           :: ISW_SNOWM   ! loop thru snow model options
INTEGER(I4B)                           :: IX_RFERR    ! index for rainfall error options
INTEGER(I4B)                           :: IX_ARCH1    ! index for upper layer architecture
INTEGER(I4B)                           :: IX_ARCH2    ! index for lower layer architecture
INTEGER(I4B)                           :: IX_QSURF    ! index for surface runoff
INTEGER(I4B)                           :: IX_QPERC    ! index for percolation
INTEGER(I4B)                           :: IX_ESOIL    ! index for evaporation
INTEGER(I4B)                           :: IX_QINTF    ! index for interflow
INTEGER(I4B)                           :: IX_Q_TDH    ! index for time delay options
INTEGER(I4B)                           :: IX_SNOWM    ! index for snow model options
INTEGER(I4B)                           :: IX_MODEL    ! model index
! (3) identify a unique model name
CHARACTER(LEN=8)                       :: CNUM        ! model index (converted to text)
INTEGER(I4B)                           :: NAME_FMT    ! format for the model name
! ---------------------------------------------------------------------------------------
! (0) INITIALIZE
! ---------------------------------------------------------------------------------------

NAME_FMT=1 ! format for the naming convention
ICOUNT  =0
IX_MODEL=0
IX_Q_TDH=0; IX_QINTF=0; IX_ESOIL=0; IX_QPERC=0; IX_QSURF=0; IX_ARCH2=0; IX_ARCH1=0; IX_RFERR=0
IX_SNOWM = 0
ERR =0
MESSAGE ='SELECTMODL/everything is fine'
! ---------------------------------------------------------------------------------------
! (1) READ IN THE COMPONENT CHOICE AND MODEL COMPONENT (AND SAVE IN A STRUCTURE)
! ---------------------------------------------------------------------------------------
! read in control file

!CFILE = TRIM(SETNGS_PATH)//M_DECISIONS      ! control file info shared in MODULE ddirectory
CFILE = TRIM(SETNGS_PATH)//'fuse_zDecisions_'//TRIM(FUSE_ID)//'.txt'      ! control file info shared in MODULE ddirectory

INQUIRE(FILE=CFILE,EXIST=LEXIST)  ! check that control file exists
IF (.not.LEXIST) THEN
  message="f-SELECTMODL/decisions file '"//trim(CFILE)//"' does not exist"
  err=100; return
ELSE
  print *, 'Now reading model decisions from:', trim(CFILE)
ENDIF

! open up model decisions file
CALL getSpareUnit(IUNIT,err,message) ! make sure IUNIT is actually available
IF (err/=0) THEN
message="f-SELECTMODL/weird/&"//message
err=100; return
ENDIF

OPEN(IUNIT,FILE=CFILE,STATUS='old')

! read format key (and strip out descriptive text)
READ(IUNIT,'(a)') KEY
CALL stripTrailString(string=KEY,trailStart='!')

! read model decisions
print *, 'Model decisions:'

DO IDEC=1,NDEC
  READ(IUNIT,KEY) M_CHOICE, DECISION
  !WRITE(*,KEY) M_CHOICE, DECISION
  PRINT *, DECISION, '-> ', M_CHOICE
  SELECT CASE (DECISION)
    CASE('RFERR'); SMODL%iRFERR = desc_str2int(M_CHOICE)
    CASE('ARCH1'); SMODL%iARCH1 = desc_str2int(M_CHOICE)
    CASE('ARCH2'); SMODL%iARCH2 = desc_str2int(M_CHOICE)
    CASE('QSURF'); SMODL%iQSURF = desc_str2int(M_CHOICE)
    CASE('QPERC'); SMODL%iQPERC = desc_str2int(M_CHOICE)
    CASE('ESOIL'); SMODL%iESOIL = desc_str2int(M_CHOICE)
    CASE('QINTF'); SMODL%iQINTF = desc_str2int(M_CHOICE)
    CASE('Q_TDH'); SMODL%iQ_TDH = desc_str2int(M_CHOICE)
    CASE('SNOWM'); SMODL%iSNOWM = desc_str2int(M_CHOICE)
    CASE DEFAULT
     message="f-SELECTMODL/UNRECOGNISED[DECISON='"//TRIM(DECISION)//"']&
                                    &[M_CHOICE='"//TRIM(M_CHOICE)//"']"
     err=100; RETURN
  ENDSELECT
END DO  ! (loop through model decisions)

! read the format for the naming convention
READ(IUNIT,*,IOSTAT=IERR) NAME_FMT; IF (IERR.NE.0) NAME_FMT=1
CLOSE(IUNIT)

!WRITE(*,'(7(A10,1X))') SMODL
! ---------------------------------------------------------------------------------------
! (2) LOOP THROUGH MODEL COMPONENTS AND IDENTIFY AN INDEX FOR EACH DECISION (AND THE MODEL)
! ---------------------------------------------------------------------------------------
! loop through model options
MODEL_OPTIONS: DO ISW_MODEL=1,1 ! (dummy loop to exit)
DO ISW_SNOWM=1,SIZE(LIST_SNOWM)  ! (snow model options)
 DO ISW_Q_TDH=1,SIZE(LIST_Q_TDH)  ! (time delay options)
  DO ISW_QINTF=1,SIZE(LIST_QINTF)  ! (interflow options)
   DO ISW_ESOIL=1,SIZE(LIST_ESOIL)  ! (evaporation options)
    DO ISW_QPERC=1,SIZE(LIST_QPERC)  ! (percolation options)
     DO ISW_QSURF=1,SIZE(LIST_QSURF)  ! (surface runoff options)
      DO ISW_ARCH2=1,SIZE(LIST_ARCH2)  ! (lower-layer architecture options)
       DO ISW_ARCH1=1,SIZE(LIST_ARCH1)  ! (upper-layer architecture options)
        DO ISW_RFERR=1,SIZE(LIST_RFERR)  ! (rainfall error options)
         ! don't allow a lower tension tank when there are two upper ones
         IF (LIST_ARCH1(ISW_ARCH1)%MCOMPONENT.EQ.'tension2_1'.AND. &
             LIST_ARCH2(ISW_ARCH2)%MCOMPONENT.EQ.'tens2pll_2') CYCLE
         ! don't allow percolation below field capacity if there are multiple upper tanks
         IF (LIST_ARCH1(ISW_ARCH1)%MCOMPONENT.NE.'onestate_1'.AND. &
             LIST_QPERC(ISW_QPERC)%MCOMPONENT.EQ.'perc_w2sat') CYCLE
         ICOUNT = ICOUNT + 1  ! (increment counter)
         ! identify unique model
         IF (SMODL%iRFERR.EQ.desc_str2int(LIST_RFERR(ISW_RFERR)%MCOMPONENT) .AND. &
             SMODL%iARCH1.EQ.desc_str2int(LIST_ARCH1(ISW_ARCH1)%MCOMPONENT) .AND. &
             SMODL%iARCH2.EQ.desc_str2int(LIST_ARCH2(ISW_ARCH2)%MCOMPONENT) .AND. &
             SMODL%iQSURF.EQ.desc_str2int(LIST_QSURF(ISW_QSURF)%MCOMPONENT) .AND. &
             SMODL%iQPERC.EQ.desc_str2int(LIST_QPERC(ISW_QPERC)%MCOMPONENT) .AND. &
             SMODL%iESOIL.EQ.desc_str2int(LIST_ESOIL(ISW_ESOIL)%MCOMPONENT) .AND. &
             SMODL%iQINTF.EQ.desc_str2int(LIST_QINTF(ISW_QINTF)%MCOMPONENT) .AND. &
             SMODL%iQ_TDH.EQ.desc_str2int(LIST_Q_TDH(ISW_Q_TDH)%MCOMPONENT) .AND. &
             SMODL%iSNOWM.EQ.desc_str2int(LIST_SNOWM(ISW_SNOWM)%MCOMPONENT)) THEN
          ! identify model components
          IX_RFERR = ISW_RFERR
          IX_ARCH1 = ISW_ARCH1
          IX_ARCH2 = ISW_ARCH2
          IX_QSURF = ISW_QSURF
          IX_QPERC = ISW_QPERC
          IX_ESOIL = ISW_ESOIL
          IX_QINTF = ISW_QINTF
          IX_Q_TDH = ISW_Q_TDH
          IX_SNOWM = ISW_SNOWM
          ! identify model
          IX_MODEL = ICOUNT
          ! exit main do loop
          EXIT MODEL_OPTIONS
         ENDIF
        END DO  ! RFERR
       END DO  ! ARCH1
      END DO  ! ARCH2
     END DO  ! QSURF
    END DO  ! QPERC
   END DO  ! ESOIL
  END DO  ! QINTF
 END DO  ! Q_TDH
END DO  ! SNOWM
END DO MODEL_OPTIONS

! check that a model was identified
IF (IX_MODEL.EQ.0) THEN
 MESSAGE='f-SELECTMODL/unable to find given model combination'
 ERR=1; RETURN
ENDIF
SMODL%MODIX = IX_MODEL  ! FUSE model index
! ---------------------------------------------------------------------------------------
! (3) IDENTIFY A UNIQUE MODEL NAME
! ---------------------------------------------------------------------------------------
! define an initial model name
SMODL%MNAME = 'FUSE_'
! convert model index to a character and append to string
WRITE(CNUM,'(I4.4)') IX_MODEL              ! convert model index to a character
SMODL%MNAME = TRIM(SMODL%MNAME)//TRIM(CNUM)//'_' ! (1 underscore)
! convert models to a character and append to the string
SELECT CASE (NAME_FMT)
 CASE(0)
  DECISION='RFERR'; CALL MAKE_MODEL(IX_RFERR,DECISION,SMODL%MNAME)
  DECISION='ARCH1'; CALL MAKE_MODEL(IX_ARCH1,DECISION,SMODL%MNAME)
  DECISION='ARCH2'; CALL MAKE_MODEL(IX_ARCH2,DECISION,SMODL%MNAME)
  DECISION='QSURF'; CALL MAKE_MODEL(IX_QSURF,DECISION,SMODL%MNAME)
  DECISION='QPERC'; CALL MAKE_MODEL(IX_QPERC,DECISION,SMODL%MNAME)
  DECISION='ESOIL'; CALL MAKE_MODEL(IX_ESOIL,DECISION,SMODL%MNAME)
  DECISION='QINTF'; CALL MAKE_MODEL(IX_QINTF,DECISION,SMODL%MNAME)
  DECISION='Q_TDH'; CALL MAKE_MODEL(IX_Q_TDH,DECISION,SMODL%MNAME)
  DECISION='SNOWM'; CALL MAKE_MODEL(IX_SNOWM,DECISION,SMODL%MNAME)
END SELECT
! add the numerix info
DECISION='NMETH'; CALL MAKE_MODEL(SOLUTION_METHOD,DECISION,SMODL%MNAME)
DECISION='SSTEP'; CALL MAKE_MODEL(TEMPORAL_ERROR_CONTROL,DECISION,SMODL%MNAME)
CONTAINS
! ---------------------------------------------------------------------------------------
SUBROUTINE MAKE_MODEL(IX,CNAME,MNAME)
! define variables
INTEGER(I4B),INTENT(IN)    :: IX     ! decision index
CHARACTER(*),INTENT(IN)    :: CNAME  ! decision name
CHARACTER(*),INTENT(INOUT) :: MNAME  ! model name
CHARACTER(LEN=8)           :: TXT_IX ! decision index as a character
WRITE(TXT_IX,'(I0)') IX                             ! convert decision index to a character
MNAME = TRIM(MNAME)//'_'//TRIM(CNAME)//'-'//TXT_IX  ! append to the string
END SUBROUTINE MAKE_MODEL
! ---------------------------------------------------------------------------------------
END SUBROUTINE SELECTMODL
END MODULE SELECTMODL_MODULE
