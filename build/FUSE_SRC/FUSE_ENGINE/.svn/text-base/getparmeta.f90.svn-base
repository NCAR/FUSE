SUBROUTINE GETPARMETA(err,message)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Reads parameter metadata 
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multiparam -- model parameters stored in MODULE multiparam
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE fuse_fileManager,only:SETNGS_PATH,CONSTRAINTS     ! defines data directory                               
USE multiparam, ONLY: PARATT                          ! parameter attribute structure
USE putpar_str_module                                 ! provide access to SUBROUTINE putpar_str
USE par_insert_module                                 ! provide access to SUBROUTINE par_insert
IMPLICIT NONE
! dummies
integer(i4b),intent(out)::err
character(*),intent(out)::message
! locals
INTEGER(I4B)                           :: IUNIT       ! file unit
integer(i4b),parameter::lenPath=1024 !DK211008: allows longer file paths
INTEGER(I4B)                           :: IERR        ! error code for read statement\
CHARACTER(LEN=lenPath)                 :: CFILE       ! name of constraints file
LOGICAL(LGT)                           :: LEXIST      ! .TRUE. if file exists
CHARACTER(LEN=256)                     :: KEY         ! format code
TYPE(PARATT)                           :: PARAM_META  ! parameter metadata
INTEGER(I4B)                           :: IPOS,JPOS   ! indices of string
INTEGER(I4B)                           :: ICH         ! looping variable (do loop)
! ---------------------------------------------------------------------------------------
! read in control file
err=0
IUNIT = 21  ! file unit
CFILE = TRIM(SETNGS_PATH)//TRIM(CONSTRAINTS)      ! control file info shared in MODULE ddirectory
INQUIRE(FILE=CFILE,EXIST=LEXIST)  ! check that control file exists
IF (.not.LEXIST) THEN
 message="f-GETPARMETA/parameter constraints file '"//trim(CFILE)//"' does not exist "
 err=100; return
ENDIF 
! initialize parameter strings
DO ICH=1,LEN(PARAM_META%P_NAME); PARAM_META%P_NAME(ICH:ICH)=' '; END DO
DO ICH=1,LEN(PARAM_META%CHILD1); PARAM_META%CHILD1(ICH:ICH)=' '; END DO
DO ICH=1,LEN(PARAM_META%CHILD2); PARAM_META%CHILD2(ICH:ICH)=' '; END DO
! open up parameter metadata file
OPEN(IUNIT,FILE=CFILE,STATUS='old')
! read format key (and strip out descriptive text)
READ(IUNIT,'(a256)') KEY
IPOS = INDEX(KEY,'!'); DO JPOS=IPOS,LEN(KEY); KEY(JPOS:JPOS)=' '; END DO
!PRINT *, TRIM(KEY), len_trim(key)
DO
 ! read parameter constraints
 READ(IUNIT,TRIM(KEY), IOSTAT=IERR) &
  PARAM_META%PARFIT,  &  ! 'fit' (T/F) [T=parameter is fitted, F=parameter is fixed at the default value)
  PARAM_META%PARSTK,  &  ! flag (0=deterministic, 1=stochastic)
  PARAM_META%PARDEF,  &  ! default parameter set
  PARAM_META%PARLOW,  &  ! lower limit of each parameter
  PARAM_META%PARUPP,  &  ! upper limit of each parameter
  PARAM_META%FRSEED,  &  ! fraction param space used as offset for "reasonable" bounds
  PARAM_META%PARSCL,  &  ! typical scale of parameter
  PARAM_META%PARVTN,  &  ! method used for variable transformation
  PARAM_META%PARDIS,  &  ! parametric form of prob dist used for prior/hyper
  PARAM_META%PARQTN,  &  ! transformation applied before use of prob dist
  PARAM_META%PARLAT,  &  ! number of latent variables (0=onePerStep, -1=from data)
  PARAM_META%PARMTH,  &  ! imeth for all variables ???what is this???
  PARAM_META%NPRIOR,  &  ! number of prior/hyper-parameters
  PARAM_META%P_NAME,  &  ! parameter name
  PARAM_META%CHILD1,  &  ! name of 1st parameter child
  PARAM_META%CHILD2      ! name of 2nd parameter child
 IF (IERR.NE.0) EXIT
 !WRITE(*,TRIM(KEY)) PARAM_META
 ! put parameters in data structures
 CALL PUTPAR_STR(PARAM_META, PARAM_META%P_NAME)
 ! populate the model parameter structure with default values
 CALL PAR_INSERT(PARAM_META%PARDEF,PARAM_META%P_NAME)
END DO
CLOSE(IUNIT)
END SUBROUTINE GETPARMETA
