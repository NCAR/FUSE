MODULE SELECTMODL_MODULE
IMPLICIT NONE
CONTAINS
! ---------------------------------------------------------------------------------------
SUBROUTINE SELECTMODL(FUSE_ID,ISTATUS,MESSAGE)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2008
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Reads a control file and identifies a unique model index
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! defines data types
USE fuse_fileManager,only:SETNGS_PATH,M_DECISIONS                 ! defines data directory                               
USE model_defn                                        ! defines model decisions
USE model_defnames
USE model_numerix, only: solution_method,temporal_error_control
IMPLICIT NONE
! Input
INTEGER(I4B), INTENT(IN), OPTIONAL     :: FUSE_ID     ! identifier for FUSE model
! Output
INTEGER(I4B), INTENT(OUT)              :: ISTATUS     ! error code
CHARACTER(LEN=*), INTENT(OUT)          :: MESSAGE     ! error message
! (1) read in the component choice and model component
LOGICAL(LGT)                           :: READ_FILE   ! .TRUE. if read model decisions from a file
integer(i4b),parameter::lenPath=1024 !DK211008: allows longer file paths
INTEGER(I4B)                           :: IUNIT       ! file unit
CHARACTER(LEN=lenPath)                 :: CFILE       ! name of constraints file
LOGICAL(LGT)                           :: LEXIST      ! .TRUE. if file exists
CHARACTER(LEN=256)                     :: KEY         ! format code
INTEGER(I4B)                           :: IPOS,JPOS   ! used to manipulate strings
INTEGER(I4B)                           :: IDEC        ! loop thru model decisions
CHARACTER(LEN=10)                      :: M_CHOICE    ! model choice
CHARACTER(LEN=5)                       :: DECISION    ! model decision
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
INTEGER(I4B)                           :: IX_RFERR    ! index for rainfall error options
INTEGER(I4B)                           :: IX_ARCH1    ! index for upper layer architecture
INTEGER(I4B)                           :: IX_ARCH2    ! index for lower layer architecture
INTEGER(I4B)                           :: IX_QSURF    ! index for surface runoff
INTEGER(I4B)                           :: IX_QPERC    ! index for percolation
INTEGER(I4B)                           :: IX_ESOIL    ! index for evaporation
INTEGER(I4B)                           :: IX_QINTF    ! index for interflow
INTEGER(I4B)                           :: IX_Q_TDH    ! index for time delay options
INTEGER(I4B)                           :: IX_MODEL    ! model index
! (3) identify a unique model name
CHARACTER(LEN=3)                       :: CNUM        ! model index (converted to text)
INTEGER(I4B)                           :: ILEN        ! length of string
INTEGER(I4B)                           :: I !,J,K       ! looping variables
INTEGER(I4B)                           :: NAME_FMT    ! format for the model name
! ---------------------------------------------------------------------------------------
! (0) INITIALIZE
! ---------------------------------------------------------------------------------------
NAME_FMT=1 ! format for the naming convention
ICOUNT  =0
IX_MODEL=0
IX_Q_TDH=0; IX_QINTF=0; IX_ESOIL=0; IX_QPERC=0; IX_QSURF=0; IX_ARCH2=0; IX_ARCH1=0; IX_RFERR=0
ISTATUS =0
MESSAGE ='everything is fine'
! ---------------------------------------------------------------------------------------
! (1) READ IN THE COMPONENT CHOICE AND MODEL COMPONENT (AND SAVE IN A STRUCTURE)
! ---------------------------------------------------------------------------------------
! check if there is a need to read data from file
IF (.NOT.PRESENT(FUSE_ID)) THEN
 READ_FILE = .TRUE.
ELSE
 READ_FILE = (FUSE_ID < 1)
ENDIF
! read in control file
IF (READ_FILE) THEN
 IUNIT = 21  ! file unit
 CFILE = TRIM(SETNGS_PATH)//TRIM(M_DECISIONS)      ! control file info shared in MODULE ddirectory
 INQUIRE(FILE=CFILE,EXIST=LEXIST)  ! check that control file exists
 IF (LEXIST) THEN
  ! open up model decisions file
  OPEN(IUNIT,FILE=CFILE,STATUS='old')
   ! read format key (and strip out descriptive text)
   READ(IUNIT,'(a256)') KEY
   IPOS = INDEX(KEY,'!'); DO JPOS=IPOS,LEN(KEY); KEY(JPOS:JPOS)=' '; END DO
   !PRINT *, TRIM(KEY)
   ! read model decisions
   DO IDEC=1,8
    READ(IUNIT,KEY) M_CHOICE, DECISION
    !WRITE(*,TRIM(KEY)) M_CHOICE, DECISION 
    SELECT CASE (DECISION)
    CASE('RFERR'); SMODL%iRFERR = desc_str2int(M_CHOICE)
    CASE('ARCH1'); SMODL%iARCH1 = desc_str2int(M_CHOICE)
    CASE('ARCH2'); SMODL%iARCH2 = desc_str2int(M_CHOICE)
    CASE('QSURF'); SMODL%iQSURF = desc_str2int(M_CHOICE)
    CASE('QPERC'); SMODL%iQPERC = desc_str2int(M_CHOICE)
    CASE('ESOIL'); SMODL%iESOIL = desc_str2int(M_CHOICE)
    CASE('QINTF'); SMODL%iQINTF = desc_str2int(M_CHOICE)
    CASE('Q_TDH'); SMODL%iQ_TDH = desc_str2int(M_CHOICE)
    CASE DEFAULT
      WRITE(*,*)"ERROR: UNRECOGNISED[DECISON='"//TRIM(DECISION)//"']&
                                   &[M_CHOICE='"//TRIM(M_CHOICE)//"']"
    ENDSELECT
   END DO  ! (loop through model decisions)
   ! read the format for the naming convention
   READ(IUNIT,*,IOSTAT=IERR) NAME_FMT; IF (IERR.NE.0) NAME_FMT=1
  CLOSE(IUNIT)
 ELSE
  STOP ' model decision file does not exist '
 ENDIF
ELSE
 SMODL=AMODL(FUSE_ID)
ENDIF
!WRITE(*,'(7(A10,1X))') SMODL

! ---------------------------------------------------------------------------------------
! (2) LOOP THROUGH MODEL COMPONENTS AND IDENTIFY AN INDEX FOR EACH DECISION (AND THE MODEL)
! ---------------------------------------------------------------------------------------
! loop through model options
MODEL_OPTIONS: DO ISW_MODEL=1,1 ! (dummy loop to exit)
DO ISW_Q_TDH=1,SIZE(LIST_Q_TDH)  ! (time delay options)
 DO ISW_QINTF=1,SIZE(LIST_QINTF)  ! (interflow options)
  DO ISW_ESOIL=1,SIZE(LIST_ESOIL)  ! (evaporation options)
   DO ISW_QPERC=1,SIZE(LIST_QPERC)  ! (percolation options)
    DO ISW_QSURF=1,SIZE(LIST_QSURF)  ! (surface runoff options)
     DO ISW_ARCH2=1,SIZE(LIST_ARCH2)  ! (lower-layer architecture options)
      DO ISW_ARCH1=1,SIZE(LIST_ARCH1)  ! (upper-layer architecture options)
       DO ISW_RFERR=1,SIZE(LIST_RFERR)  ! (rainfall error options)
        ! don't allow a lower tension tank when there are two upper ones
        IF (LIST_ARCH1(ISW_ARCH1)%MCOMPONENT(1:10).EQ.'tension2_1'.AND. &
            LIST_ARCH2(ISW_ARCH2)%MCOMPONENT(1:10).EQ.'tens2pll_2') CYCLE
        ! don't allow percolation below field capacity if there are multiple upper tanks
        IF (LIST_ARCH1(ISW_ARCH1)%MCOMPONENT(1:10).NE.'onestate_1'.AND. &
            LIST_QPERC(ISW_QPERC)%MCOMPONENT(1:10).EQ.'perc_w2sat') CYCLE
        ICOUNT = ICOUNT + 1  ! (increment counter)
        ! identify unique model
        IF (SMODL%iRFERR.EQ.desc_str2int(LIST_RFERR(ISW_RFERR)%MCOMPONENT) .AND. &
            SMODL%iARCH1.EQ.desc_str2int(LIST_ARCH1(ISW_ARCH1)%MCOMPONENT) .AND. &
            SMODL%iARCH2.EQ.desc_str2int(LIST_ARCH2(ISW_ARCH2)%MCOMPONENT) .AND. &
            SMODL%iQSURF.EQ.desc_str2int(LIST_QSURF(ISW_QSURF)%MCOMPONENT) .AND. &
            SMODL%iQPERC.EQ.desc_str2int(LIST_QPERC(ISW_QPERC)%MCOMPONENT) .AND. &
            SMODL%iESOIL.EQ.desc_str2int(LIST_ESOIL(ISW_ESOIL)%MCOMPONENT) .AND. &
            SMODL%iQINTF.EQ.desc_str2int(LIST_QINTF(ISW_QINTF)%MCOMPONENT) .AND. &
            SMODL%iQ_TDH.EQ.desc_str2int(LIST_Q_TDH(ISW_Q_TDH)%MCOMPONENT)) THEN
         ! identify model components
         IX_RFERR = ISW_RFERR
         IX_ARCH1 = ISW_ARCH1
         IX_ARCH2 = ISW_ARCH2
         IX_QSURF = ISW_QSURF
         IX_QPERC = ISW_QPERC
         IX_ESOIL = ISW_ESOIL
         IX_QINTF = ISW_QINTF
         IX_Q_TDH = ISW_Q_TDH
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
END DO MODEL_OPTIONS
! check that a model was identified
IF (IX_MODEL.EQ.0) THEN
 ISTATUS=1
 MESSAGE='unable to find given model combination'
ENDIF
SMODL%MODIX = IX_MODEL  ! FUSE model index
! ---------------------------------------------------------------------------------------
! (3) IDENTIFY A UNIQUE MODEL NAME 
! ---------------------------------------------------------------------------------------
! define an initial model name
SMODL%MNAME = 'FUSE_'; DO I=6,LEN(SMODL%MNAME); SMODL%MNAME(I:I)=' '; END DO
IPOS=LEN_TRIM(SMODL%MNAME)+1  ! (get the position of the next character)
! convert model index to a character and append to string
WRITE(CNUM,'(I3)') ABS(IX_MODEL); CNUM=ADJUSTR(CNUM)  ! convert model index to a character
IF (ABS(IX_MODEL).LT.100) CNUM(1:1) = '0'             ! pad with zeroes
IF (ABS(IX_MODEL).LT. 10) CNUM(2:2) = '0'             ! pad with zeroes
ILEN=LEN(CNUM)+1; SMODL%MNAME(IPOS:IPOS+ILEN-1) = CNUM//'_' ! (1 underscore)
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
CHARACTER(LEN=1)           :: TXT_IX ! decision index as a character
! get the position of the next character
IPOS=LEN_TRIM(MNAME)+1
! convert decision index to a character
WRITE(TXT_IX,'(I1)') IX
! get the length of the string
ILEN = LEN(CNAME)+3
! append to the string
MNAME(IPOS:IPOS+ILEN-1) = '_'//CNAME//'-'//TXT_IX
END SUBROUTINE MAKE_MODEL
! ---------------------------------------------------------------------------------------
END SUBROUTINE SELECTMODL
END MODULE SELECTMODL_MODULE
