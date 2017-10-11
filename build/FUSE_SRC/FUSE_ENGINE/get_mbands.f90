module get_mbands_module
USE nrtype
USE netcdf
implicit none
private
public::GET_MBANDS
public::GET_MBANDS_INFO
contains

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
	print *, 'f-GET_MBANDS/control file ',TRIM(CFILE),' for band data does not exist '
	STOP
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


SUBROUTINE GET_MBANDS_INFO(ELEV_BANDS_NC,err,message)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Created by Nans Addor, 2/2017
! Based on Brian Henn's GET_MBANDS, 7/2013
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Read band data (elevation and area fraction) from a NetCDF grid
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multibands -- populate structure MBANDS_INFO_3d and Z_FORCING_grid
! ---------------------------------------------------------------------------------------
use nrtype,only:I4B,LGT,SP
use utilities_dmsl_kit_FUSE,only:getSpareUnit,stripTrailString
USE fuse_fileManager,only:INPUT_PATH,SETNGS_PATH,MBANDS_INFO     ! defines data directory
USE multibands,only:N_BANDS,MBANDS,MBANDS_INFO_3d,Z_FORCING,&
										Z_FORCING_grid,elev_mask          ! model band structures
USE multiforce,only:nspat1,nspat2                     ! dimension lengths

IMPLICIT NONE
! dummies
CHARACTER(LEN=1024),intent(in)      :: ELEV_BANDS_NC
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
INTEGER(I4B)                           :: BAND_START  ! index of start of band info
INTEGER(I4B)                           :: BAND_END    ! index of end of band info
INTEGER(I4B)                           :: IBANDS      ! band index (input data)
INTEGER(I4B)                           :: JBAND       ! band index (internal data structure)
INTEGER(I4B)                           :: NCID_EB     ! NetCDF ID for elevation bands file
INTEGER(I4B)                           :: iSpat1,iSpat2  ! loop through spatial dimensions
REAL(SP),dimension(:,:,:),allocatable  :: AF_TEMP, ME_TEMP ! Temporary data structures to store area_frac and mean_area

! internal: NetCDF read
integer(i4b)                           :: ivarid_af,ivarid_me  ! NetCDF variable ID for area_frac and mean_area
integer(i4b),parameter                 :: ndims=3     ! number of dimensions for frac_area
integer(i4b) 										       :: dimid_eb    ! ID elevation bands
integer(i4b)                           :: iDimID      ! dimension ID
integer(i4b)                           :: dimLen      ! dimension length

! ---------------------------------------------------------------------------------------
! read in NetCDF file defining the elevation bands - no info file needed for the gridded version
err=0
CFILE = TRIM(INPUT_PATH)//ELEV_BANDS_NC      ! control file info shared in MODULE directory
print *, 'Elevation bands info file:',TRIM(CFILE)

INQUIRE(FILE=CFILE,EXIST=LEXIST)  ! check that control file exists
IF (.NOT.LEXIST) THEN
	print *, 'f-GET_MBANDS_GRID/NetCDF file ',TRIM(CFILE),' for band data does not exist '
	STOP
ENDIF

!open netcdf file
err = nf90_open(CFILE, nf90_nowrite, NCID_EB)
if (err.ne.0) write(*,*) trim(message); if (err.gt.0) stop
PRINT *, 'NCID_EB is', NCID_EB

! get the dimension IDs for elevation_band
ierr = nf90_inq_dimid(NCID_EB, 'elevation_band', dimid_eb)
if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif

! get dimension length
ierr = nf90_inquire_dimension(ncid_eb,dimid_eb,len=dimLen)
if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif

! save the dimension lengths
N_BANDS = dimLen  ! number of elevation bands
print *, 'N_BANDS = ', N_BANDS

! get the variable ID for the fraction of the area contained in each elevation band
ierr = nf90_inq_varid(NCID_EB, 'area_frac', ivarid_af)
if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif
ierr = nf90_inq_varid(NCID_EB, 'mean_elev', ivarid_me)
if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif

! allocate 1 data stucture
ALLOCATE(MBANDS(N_BANDS),STAT=IERR(1))

! allocate data structures
ALLOCATE(Z_FORCING_grid(nspat1,nspat2),MBANDS_INFO_3d(nspat1,nspat2,n_bands),&
				 AF_TEMP(nspat1,nspat2,n_bands),ME_TEMP(nspat1,nspat2,n_bands),&
				 elev_mask(nspat1,nspat2),STAT=IERR(1))

IF (ANY(IERR.NE.0)) THEN
 message="f-GET_MBANDS/problem allocating elevation band data structures"
 err=100; return
ENDIF

! import data into temporary stuctures
ierr = nf90_get_var(NCID_EB, ivarid_af, AF_TEMP, start=(/1,1,1/), count=(/nSpat1,nSpat2,n_bands/)); CALL HANDLE_ERR(IERR)
if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif

! import data into temporary stuctures
ierr = nf90_get_var(NCID_EB, ivarid_me, me_TEMP, start=(/1,1,1/), count=(/nSpat1,nSpat2,n_bands/)); CALL HANDLE_ERR(IERR)
if(err/=0)then; message=trim(message)//trim(nf90_strerror(err)); return; endif

! populate MBANDS_INFO_3d, Z_FORCING_grid and elev_mask
DO iSpat2=1,nSpat2
	DO iSpat1=1,nSpat1

	 MBANDS_INFO_3d(iSpat1,iSpat2,:)%Z_MID = me_TEMP(iSpat1,iSpat2,:)
	 MBANDS_INFO_3d(iSpat1,iSpat2,:)%AF    = af_TEMP(iSpat1,iSpat2,:)
	 Z_FORCING_grid(iSpat1,iSpat2)    = sum(me_TEMP(iSpat1,iSpat2,:)*af_TEMP(iSpat1,iSpat2,:)) ! estimate mean elevation of forcing using weighted mean of EB elevation
	 elev_mask(iSpat1,iSpat2)=me_TEMP(iSpat1,iSpat2,1)

	 PRINT *, 'Z_FORCING_grid =', Z_FORCING_grid(iSpat1,iSpat2)
	 PRINT *, 'MBANDS_INFO_3d - ELEV =', MBANDS_INFO_3d(iSpat1,iSpat2,:)%Z_MID
	 PRINT *, 'MBANDS_INFO_3d - FRAC =', MBANDS_INFO_3d(iSpat1,iSpat2,:)%AF

	 if (abs(sum(MBANDS_INFO_3d(iSpat1,iSpat2,:)%AF)-1).GT.1E-6) then ! check that area fraction sum to 1

 	  print *, 'DIF EB = ', abs(sum(MBANDS_INFO_3d(iSpat1,iSpat2,:)%AF)-1)
	 	print *, "f-GET_MBANDS/area fraction of elevation bands do not sum to 1" ! TODO: use message instead?
		stop

	 end if

	END DO
END DO

err = nf90_close(ncid_eb)
if (err.ne.0) write(*,*) trim(message); if (err.gt.0) stop

DEALLOCATE(AF_TEMP, ME_TEMP)

print *, 'Done populating data structures for elevation bands'

END SUBROUTINE GET_MBANDS_INFO

end module get_mbands_module
