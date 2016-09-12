SUBROUTINE GET_MBANDS(err,message)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Created by Brian Henn, 7/2013
! Based on GETFORCING.f90 by Martyn Clark, 2009
! Updated by Dmitri Kavetski, 14 Sept 2014 AD - Chiefleys Newie
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Read ASCII basin band data in BATEA format
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multibands -- populate structure MBANDS(*)%(*)
! ---------------------------------------------------------------------------------------
use nrtype,only:I4B,LGT,SP
use utilities_dmsl_kit_FUSE,only:getSpareUnit,stripTrailString
USE fuse_fileManager,only:INPUT_PATH,SETNGS_PATH,MBANDS_INFO     ! defines data directory 
USE multibands,only:N_BANDS,MBANDS,Z_FORCING          ! model band structures
IMPLICIT NONE
! dummies
integer(I4B), intent(out)              :: err
character(*), intent(out)              :: message
! internal
integer(i4b),parameter::lenPath=1024 ! DK/2008/10/21: allows longer file paths
INTEGER(I4B),DIMENSION(2)              :: IERR        ! error codes
INTEGER(I4B)                           :: IUNIT       ! input file unit
CHARACTER(LEN=lenPath)                 :: CFILE       ! name of control file
CHARACTER(LEN=lenPath)                 :: BFILE       ! name of band file
LOGICAL(LGT)                           :: LEXIST      ! .TRUE. if control file exists
CHARACTER(LEN=lenPath)                 :: FNAME_INPUT ! name of band input file
INTEGER(I4B)                           :: NCOLB       ! number of band columns
INTEGER(I4B)                           :: IX_Z        ! column number for band elevation
INTEGER(I4B)                           :: IX_AF       ! column number for band area fraction
INTEGER(I4B)                           :: NHEADB      ! number of band header rows
INTEGER(I4B)                           :: BAND_START  ! index of start of band info
INTEGER(I4B)                           :: BAND_END    ! index of end of band info
INTEGER(I4B)                           :: IHEAD       ! header index
CHARACTER(LEN=lenPath)                 :: TMPTXT      ! descriptive text
INTEGER(I4B)                           :: IBANDS      ! band index (input data)
INTEGER(I4B)                           :: JBAND       ! band index (internal data structure)
REAL(SP),DIMENSION(:),ALLOCATABLE      :: TMPDAT      ! one line of data
! ---------------------------------------------------------------------------------------
! read in control file
err=0
CFILE = TRIM(SETNGS_PATH)//MBANDS_INFO      ! control file info shared in MODULE directory
print *, 'Elevation bands info file:',TRIM(CFILE)

INQUIRE(FILE=CFILE,EXIST=LEXIST)  ! check that control file exists
IF (.NOT.LEXIST) THEN
 message='f-GET_MBANDS/control file "'//TRIM(CFILE)//'" for band data does not exist ' 
 err=100; return
ENDIF
! read in parameters of the control files
CALL getSpareUnit(IUNIT,err,message) ! make sure IUNIT is actually available
IF (err/=0) THEN
 message="f-GET_MBANDS/weird/&"//message
 err=100; return
ENDIF
OPEN(IUNIT,FILE=CFILE,STATUS='old')           
READ(IUNIT,'(A)') FNAME_INPUT                         ! get input filename
! number of columns and column numbers
READ(IUNIT,*) NCOLB,IX_Z,IX_AF                           ! band data: number of columns, elevation, area fraction
READ(IUNIT,*) NHEADB,N_BANDS,BAND_START,BAND_END         ! number of headers, number of bands, first band line, last band line	
CLOSE(IUNIT)
! fill extra characters in filename with white space
CALL stripTrailString(string=FNAME_INPUT,trailStart='!')
IF (N_BANDS.NE.(BAND_END-BAND_START+1)) THEN
 message="f-GET_MBANDS/N_BANDS does not match the number of band lines in the band file"
 err=100; return
ENDIF
! ---------------------------------------------------------------------------------------
! read band data
ALLOCATE(MBANDS(N_BANDS),STAT=IERR(1))        ! (shared in module multibands)
ALLOCATE(TMPDAT(NCOLB),STAT=IERR(2))          ! (only used in this routine -- deallocate later)
IF (ANY(IERR.NE.0)) THEN
 message="f-GET_MBANDS/problem allocating data structures"
 err=100; return
ENDIF
JBAND = 0
BFILE = TRIM(INPUT_PATH)//FNAME_INPUT
INQUIRE(FILE=BFILE,EXIST=LEXIST)  ! check that control file exists
IF (.NOT.LEXIST) THEN
 print *, 'f-GET_MBANDS/band data file '//TRIM(BFILE)//' does not exist '
 err=100; return
ENDIF
CALL getSpareUnit(IUNIT,err,message) ! make sure IUNIT is actually available
IF (err/=0) THEN
 message="f-GET_MBANDS/weird/&"//message
 err=100; return
ENDIF
OPEN(IUNIT,FILE=BFILE,STATUS='old')
! read header
DO IHEAD=1,NHEADB
 IF (IHEAD.EQ.2) THEN
  READ(IUNIT,*) Z_FORCING ! elevation of the forcing data (shared in module multibands)
 ELSE
  READ(IUNIT,*) TMPTXT    ! descriptive text
 ENDIF
END DO

print *, 'Z_FORCING', Z_FORCING

! read data
DO IBANDS=1,N_BANDS
 READ(IUNIT,*) TMPDAT
 JBAND = JBAND+1
 MBANDS(JBAND)%NUM   = INT(TMPDAT(1))
 MBANDS(JBAND)%Z_MID = TMPDAT(IX_Z)
 MBANDS(JBAND)%AF    = TMPDAT(IX_AF)
END DO
CLOSE(IUNIT)
DEALLOCATE(TMPDAT, STAT=IERR(1))
IF (IERR(1).NE.0) THEN 
 message='f-GET_MBANDS/problem deallocating TMPDAT'
 err=100; return
END IF
! ---------------------------------------------------------------------------------------
END SUBROUTINE GET_MBANDS
