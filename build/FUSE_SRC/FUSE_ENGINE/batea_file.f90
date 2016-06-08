SUBROUTINE BATEA_FILE()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Used to write parameter files in BATEA format
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE fuse_fileManager,only:SETNGS_PATH,BATEA_PARAM     ! defines data directory                               
USE multiparam, ONLY: PARATT, LPARAM, NUMPAR          ! parameter attribute structure
USE getpar_str_module                                 ! provide access to SUBROUTINE getpar_str
IMPLICIT NONE
INTEGER(I4B)                           :: I           ! FORALL loop
CHARACTER(LEN=90)                      :: CNEW        ! new parameter delimiter
integer(i4b),parameter::lenPath=1024 !DK211008: allows longer file paths
CHARACTER(LEN=lenPath)                 :: CFILE       ! name of constraints file
INTEGER(I4B)                           :: IUNIT       ! file unit
INTEGER(I4B)                           :: IPAR        ! loop through parameters
INTEGER(I4B)                           :: IHYP        ! loop through hyper-parameters
INTEGER(I4B)                           :: IPRI        ! loop through prior-parameters
INTEGER(I4B)                           :: NPARFIT     ! number of fitted prior/hyper params
CHARACTER(LEN=256)                     :: PARNAME     ! parameter name
TYPE(PARATT)                           :: PARAM_MODEL ! parameter metadata (model parameters)
TYPE(PARATT)                           :: PARAM_HYPER ! parameter metadata (hyper-parameters)
TYPE(PARATT)                           :: PARAM_PRIOR ! parameter metadata (prior-parameters)
! ---------------------------------------------------------------------------------------
! initialize
CNEW(1:1)='!'
FORALL(I=2:LEN(CNEW)) CNEW(I:I)='*'                   ! define break
FORALL(I=1:LEN(PARNAME)) PARNAME(I:I)=' '
! open up batea output file
IUNIT = 21  ! file unit
CFILE = TRIM(SETNGS_PATH)//TRIM(BATEA_PARAM) ! file info shared in MODULE ddirectory
OPEN(IUNIT,FILE=CFILE,STATUS='unknown')
! write file header
WRITE(IUNIT,'(A)') 'DMDL_FARAMINEAUX_INFERN_FILE_V2'
WRITE(IUNIT,'(A)') '"Example of a faramineux infern file" ! file comment (not used)'
WRITE(IUNIT,'(A1)') ' '
WRITE(IUNIT,'(I1,1X,A19,1X,A)') 2, ' ', '! modelID (consult dynamicModelLibrary), 2=FUSE'
WRITE(IUNIT,'(A1)') ' '
WRITE(IUNIT,'(A90)') CNEW; WRITE(IUNIT,'(A1)') ' ' ! write delimiter plus blank line
! loop through parameters
DO IPAR=1,NUMPAR
 ! --------------------------------------------------------------------------------------
 ! get parameter metadata
 CALL GETPAR_STR(LPARAM(IPAR)%PARNAME,PARAM_MODEL)     ! get parameter metadata  
 ! write parameter title and parameter index
 WRITE(IUNIT,'(A,1X,I2,1X,A1,1X,A)') &
  'NEW_PAR_000 - Parameter', IPAR, '-', '"'//TRIM(LPARAM(IPAR)%PARNAME)//'"'
 WRITE(IUNIT,'(I2.2,1X,A18,1X,A)') IPAR, ' ', '! i, index of parameter'
 ! write parameter info
 CALL WRITE_PARINFO(PARAM_MODEL)                       ! write parameter info
 ! --------------------------------------------------------------------------------------
 ! check for hyper-parameter
 IF (PARAM_MODEL%NPRIOR.GT.0) THEN
  NPARFIT=0
  DO IHYP=1,PARAM_MODEL%NPRIOR
   ! identify name of child parameter
   IF (IHYP.EQ.1) PARNAME(1:LEN(PARAM_MODEL%CHILD1))=PARAM_MODEL%CHILD1(1:LEN(PARAM_MODEL%CHILD1))
   IF (IHYP.EQ.2) PARNAME(1:LEN(PARAM_MODEL%CHILD2))=PARAM_MODEL%CHILD2(1:LEN(PARAM_MODEL%CHILD2))
   IF (IHYP.GT.2) STOP ' only anticipate that there will ever by two hyper-parameters '
   ! get parameter metadata
   CALL GETPAR_STR(TRIM(PARNAME),PARAM_HYPER)          ! get parameter metadata
   ! keep track of the number of fitted parameters
   IF (PARAM_HYPER%PARFIT) NPARFIT = NPARFIT+1
   ! write parameter header
   WRITE(IUNIT,'(A1)') ' '
   FORALL(I=2:LEN(CNEW)) CNEW(I:I)='-'                 ! define new break
   WRITE(IUNIT,'(A90)') CNEW; WRITE(IUNIT,'(A1)') ' '  ! delimiter plus blank line
   ! write parameter title and parameter index
   WRITE(IUNIT,'(A,1X,I2,1X,A,1X,A1,1X,A)') 'NEW_PAR_000 - Hyper parameter', IHYP, &
    ' of "'//TRIM(LPARAM(IPAR)%PARNAME)//'"', '-', '"'//TRIM(PARAM_HYPER%P_NAME)//'"'
   WRITE(IUNIT,'(I2.2,1X,A18,1X,A)') IHYP, ' ', '! k, index of parameter'
   ! write parameter data
   CALL WRITE_PARINFO(PARAM_HYPER)
   ! --------------------------------------------------------------------------------------
   ! check for prior-parameter
   IF (PARAM_HYPER%NPRIOR.GT.0) THEN
    DO IPRI=1,PARAM_HYPER%NPRIOR
     ! identify name of child parameter
     IF (IPRI.EQ.1) PARNAME(1:LEN(PARAM_HYPER%CHILD1))=PARAM_HYPER%CHILD1(1:LEN(PARAM_HYPER%CHILD1))
     IF (IPRI.EQ.2) PARNAME(1:LEN(PARAM_HYPER%CHILD2))=PARAM_HYPER%CHILD2(1:LEN(PARAM_HYPER%CHILD2))
     IF (IPRI.GT.2) STOP ' only anticipate that there will ever by two prior-parameters '
     ! get parameter metadata
     CALL GETPAR_STR(TRIM(PARNAME),PARAM_PRIOR)         ! get parameter metadata
     ! write parameter header
     WRITE(IUNIT,'(A1)') ' '
     WRITE(IUNIT,'(A90)') CNEW; WRITE(IUNIT,'(A1)') ' ' ! delimiter plus blank line
     ! write parameter title and parameter index
     WRITE(IUNIT,'(A,1X,I2,1X,A,1X,I2,1X,A,1X,A1,1X,A)') 'NEW_PAR_000 - Prior parameter', IPRI, &
      ' of Hyper parameter', IHYP, ' of "'//TRIM(LPARAM(IPAR)%PARNAME)//'"', '-', &
      '"'//TRIM(PARAM_PRIOR%P_NAME)//'"'
     WRITE(IUNIT,'(I2.2,1X,A18,1X,A)') IHYP, ' ', '! k, index of parameter'
     ! write parameter data
     CALL WRITE_PARINFO(PARAM_PRIOR)
     ! write end text for prior parameter
     WRITE(IUNIT,'(A)') 'INF_LIST'
     WRITE(IUNIT,'(I1,1X,A19,1X,A)') 0, ' ', '! number of fitted prior/hyper-parameters'
     WRITE(IUNIT,'(I1,1X,A19,1X,A)') 0, ' ', '! list of fitted prior/hyper-parameters'
     WRITE(IUNIT,'(A,1X,I2,1X,A,1X,I2,1X,A,1X,A1,1X,A)') 'END_PAR_000 - Prior parameter', IPRI, &
      ' of Hyper parameter', IHYP, ' of "'//TRIM(LPARAM(IPAR)%PARNAME)//'"', '-', &
      '"'//TRIM(PARAM_PRIOR%P_NAME)//'"'
    END DO   ! (loop through prior parameters)
   ENDIF    ! (if there are prior parameters)
   ! write end text for hyper parameter
   IF (PARAM_HYPER%NPRIOR.GT.0) THEN
    WRITE(IUNIT,'(A1)') ' '
    WRITE(IUNIT,'(A90)') CNEW; WRITE(IUNIT,'(A1)') ' ' ! delimiter plus blank line
   ENDIF
   WRITE(IUNIT,'(A)') 'INF_LIST'
   WRITE(IUNIT,'(I1,1X,A19,1X,A)') 0, ' ', '! number of fitted prior/hyper-parameters'
   WRITE(IUNIT,'(I1,1X,A19,1X,A)') 0, ' ', '! list of fitted prior/hyper-parameters'
   WRITE(IUNIT,'(A,1X,I2,1X,A,1X,A1,1X,A)') 'END_PAR_000 - Hyper parameter', IHYP, &
    ' of "'//TRIM(LPARAM(IPAR)%PARNAME)//'"', '-', '"'//TRIM(PARAM_HYPER%P_NAME)//'"'
  END DO   ! (loop through hyper parameters)
  ! write end text for parameter
  WRITE(IUNIT,'(A1)') ' '
  WRITE(IUNIT,'(A90)') CNEW; WRITE(IUNIT,'(A1)') ' ' ! delimiter plus blank line
  WRITE(IUNIT,'(A)') 'INF_LIST'
  WRITE(IUNIT,'(I1,1X,A19,1X,A)') NPARFIT, ' ', '! number of fitted prior/hyper-parameters'
  IF (NPARFIT.EQ.0) WRITE(IUNIT,'(I1,1X,A19,1X,A)') 0, ' ', '! list of fitted prior/hyper-parameters'
  IF (NPARFIT.EQ.1) WRITE(IUNIT,'(I1,1X,A19,1X,A)') 1, ' ', '! list of fitted prior/hyper-parameters'
  IF (NPARFIT.EQ.2) WRITE(IUNIT,'(I1,A1,I1,A17,1X,A)') 1,',',2, ' ', '! list of fitted prior/hyper-parameters'
 ELSE
  ! write end text for parameter
  WRITE(IUNIT,'(A)') 'INF_LIST'
  WRITE(IUNIT,'(I1,1X,A19,1X,A)') 0, ' ', '! number of fitted prior/hyper-parameters'
  WRITE(IUNIT,'(I1,1X,A19,1X,A)') 0, ' ', '! list of fitted prior/hyper-parameters'
 ENDIF    ! (if there are hyper parameters)
 ! continue writing end text (same in both cases)
 WRITE(IUNIT,'(A,1X,I2,1X,A1,1X,A)') &
  'NEW_PAR_000 - Parameter', IPAR, '-', '"'//TRIM(LPARAM(IPAR)%PARNAME)//'"'
 WRITE(IUNIT,'(A1)') ' '
 FORALL(I=2:LEN(CNEW)) CNEW(I:I)='='                ! re-define delimiter
 WRITE(IUNIT,'(A90)') CNEW; WRITE(IUNIT,'(A1)') ' ' ! write delimiter plus blank line
END DO  ! loop through parameters
! write final delimiter
FORALL(I=2:LEN(CNEW)) CNEW(I:I)='*'                   ! define break
WRITE(IUNIT,'(A90)') CNEW; WRITE(IUNIT,'(A1)') ' ' ! delimiter plus blank line
CLOSE(IUNIT)
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
CONTAINS
 SUBROUTINE WRITE_PARINFO(PARAM_META)
 ! define parameter metadata structure
 TYPE(PARATT), INTENT(IN)             :: PARAM_META  ! parameter metadata
 REAL(SP)                             :: PAR_OFFSET  ! used to define "reasonable" parameter range
 ! write 1st block
 WRITE(IUNIT,'(A11, 1X, A9,1X,A)')  '"'//TRIM(PARAM_META%P_NAME)//'"     ', ' ', '! name of parameter'
 WRITE(IUNIT,'(L1,  1X,A19,1X,A)') PARAM_META%PARFIT, ' ', '! fit (T/F) [T=param fitted, F=param fixed at default]'
 WRITE(IUNIT,'(I1,  1X,A19,1X,A)') PARAM_META%PARSTK, ' ', '! flag (0=deterministic, 1=stochastic)'
 WRITE(IUNIT,'(A1)') ' '
 ! write 2nd block
 WRITE(IUNIT,'(E9.3,A1,E9.3,1X,A1,1X,A)') PARAM_META%PARLOW, ',', PARAM_META%PARUPP, ' ', &
  '! pLo and pHi:    Bounds on parameter'
 PAR_OFFSET = PARAM_META%FRSEED * (PARAM_META%PARUPP - PARAM_META%PARLOW)
 WRITE(IUNIT,'(E9.3,A1,E9.3,1X,A1,1X,A)')  PARAM_META%PARLOW+PAR_OFFSET, ',', PARAM_META%PARUPP-PAR_OFFSET, ' ', &
  '! pLoR and pHiR:  Reasonable bounds on parameter (seeding multi-sequences)'
 WRITE(IUNIT,'(E9.3,1X,A11,1X,A)') PARAM_META%PARSCL, ' ', '! typical scale of parameter'
 WRITE(IUNIT,'(E9.3,1X,A11,1X,A)') PARAM_META%PARDEF, ' ', '! initial value of parameter'
 WRITE(IUNIT,'(I1,  1X,A19,1X,A)') PARAM_META%PARVTN, ' ', '! ftran_v2z: fitting z-transform [see transformation library'
 WRITE(IUNIT,'(A1)') ' '
 ! write 3rd block
 WRITE(IUNIT,'(I1,  1X,A19,1X,A)') PARAM_META%PARDIS, ' ', '! prDistID - prior (det) or hyper (stok) [see distribution library]'
 WRITE(IUNIT,'(I1,  1X,A19,1X,A)') PARAM_META%PARQTN, ' ', '! ptran_v2q - probModel-transform [see transformation library]'
 WRITE(IUNIT,'(A1)') ' '
 ! write 4th block
 WRITE(IUNIT,'(I1,  1X,A19,1X,A)') PARAM_META%PARLAT, ' ', '! number of latents (ignored for det, stk: 0=onePerStep, -1=from data'
 WRITE(IUNIT,'(I1,  1X,A19,1X,A)') PARAM_META%PARMTH, ' ', '! imeth for all vars FXD_IMETH=0,EXP_IMETH=1,LIN_IMETH=2,FBF_IMETH=4'
 WRITE(IUNIT,'(A1)') ' '
 ! write 5th block
 WRITE(IUNIT,'(I1,  1X,A19,1X,A)') 0, ' ', '! number of auxiliaries needed'
 WRITE(IUNIT,'(I1,  1X,A19,1X,A)') 0, ' ', '! list of auxiliaries needed'
 WRITE(IUNIT,'(A1)') ' '
 ! write 6th block
 WRITE(IUNIT,'(I1,  1X,A19,1X,A)') PARAM_META%NPRIOR, ' ', '! number of prior/hyper-parameters'
 END SUBROUTINE WRITE_PARINFO
END SUBROUTINE BATEA_FILE
