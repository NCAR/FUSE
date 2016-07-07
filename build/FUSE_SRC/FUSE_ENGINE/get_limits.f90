SUBROUTINE GET_LIMITS()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007; revised 2008 to make use of parameter names; 
!  revised 2009 to include extra information for BATEA
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Reads parameter constraints 
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multiparam -- model parameters stored in MODULE multiparam
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE fuse_fileManager,only:SETNGS_PATH,CONSTRAINTS     ! defines data directory                               
USE multiparam, ONLY: PARATT                          ! parameter attribute structure
USE putpar_str_module                                 ! provide access to SUBROUTINE putpar_str
IMPLICIT NONE
INTEGER(I4B)                           :: IUNIT       ! file unit
integer(i4b),parameter::lenPath=1024 !DK211008: allows longer file paths
INTEGER(I4B)                           :: IERR        ! error code for read statement\
REAL(SP)                               :: XVAR        ! argument for SUBROUTINE putpar_str
CHARACTER(LEN=lenPath)                 :: CFILE       ! name of constraints file
LOGICAL(LGT)                           :: LEXIST      ! .TRUE. if file exists
CHARACTER(LEN=256)                     :: KEY         ! format code
TYPE(PARATT)                           :: PARAM_META  ! parameter metadata
INTEGER(I4B)                           :: IPOS,JPOS   ! indices of string
INTEGER(I4B)                           :: ICH         ! looping variable (forall loop)
! ---------------------------------------------------------------------------------------
print *, 'in get_limits'
! read in control file
IUNIT = 21  ! file unit
CFILE = TRIM(SETNGS_PATH)//TRIM(CONSTRAINTS)      ! control file info shared in MODULE directory
INQUIRE(FILE=CFILE,EXIST=LEXIST)  ! check that control file exists
IF (LEXIST) THEN
 ! initialize parameter strings
 FORALL(ICH=1:LEN(PARAM_META%P_NAME)) PARAM_META%P_NAME(ICH:ICH)=' '
 FORALL(ICH=1:LEN(PARAM_META%CHILD1)) PARAM_META%CHILD1(ICH:ICH)=' '
 FORALL(ICH=1:LEN(PARAM_META%CHILD2)) PARAM_META%CHILD2(ICH:ICH)=' '
 ! open up model decisions file
 OPEN(IUNIT,FILE=CFILE,STATUS='old')
  ! read format key (and strip out descriptive text)
  READ(IUNIT,'(a256)') KEY
  IPOS = INDEX(KEY,'!'); FORALL(JPOS=IPOS:LEN(KEY)) KEY(JPOS:JPOS)=' '
  PRINT *, TRIM(KEY), len_trim(key)
  DO
   ! read parameter constraints
   READ(IUNIT,TRIM(KEY), IOSTAT=IERR) &
    PARAM_META%PARFIT,  &  ! 'fit' (T/F) [T=parameter is fitted, F=parameter is fixed at the default value)
    PARAM_META%PARSTK,  &  ! flag (0=deterministic, 1=stochastic)
    PARAM_META%PARDEF,  &  ! default parameter set
    PARAM_META%PARLOW,  &  ! lower limit of each parameter
    PARAM_META%PARUPP,  &  ! upper limit of each parameter
    PARAM_META%FRSEED,  &  ! fraction param space for "reasonable" bounds
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
   WRITE(*,TRIM(KEY)) PARAM_META
   ! put parameters in data structures
   CALL PUTPAR_STR(PARAM_META, PARAM_META%P_NAME)
  END DO
 CLOSE(IUNIT)
ELSE
 STOP ' parameter constraints file does not exist '
ENDIF 
END SUBROUTINE GET_LIMITS
