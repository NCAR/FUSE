SUBROUTINE GET_SMODEL(NETCDF_FILE,IMOD)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Read model decisions from a NetCDF output file
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE model_defn -- populate structure SMODL
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! data types, etc.
! USE fuse_fileManager,only                                        ! defines data directory
USE model_defn                                        ! model definition structures
IMPLICIT NONE
! input
CHARACTER(LEN=*), INTENT(IN)           :: NETCDF_FILE ! NetCDF file name
INTEGER(I4B), INTENT(IN)               :: IMOD        ! model index
! internal
LOGICAL(LGT)                           :: LEXIST      ! .TRUE. if NetCDF file exists
INTEGER(I4B)                           :: IERR        ! error code
INTEGER(I4B)                           :: NCID        ! NetCDF file ID
INTEGER(I4B)                           :: IDIMID      ! NetCDF dimension ID
INTEGER(I4B)                           :: IVARID      ! NetCDF variable ID
INTEGER(I4B)                           :: NDESC       ! number of model descriptors
INTEGER(I4B)                           :: NCHAR       ! length of model descriptors
INTEGER(I4B)                           :: IDESC       ! loop thru model descriptors
INTEGER(I4B), DIMENSION(3)             :: ISTART      ! start indices for data read
INTEGER(I4B), DIMENSION(3)             :: ICOUNT      ! number of elements read in each dimension
CHARACTER(LEN=50)                      :: TXTVEC      ! text vector
! ---------------------------------------------------------------------------------------
! CONTENT REMOVED FOR COPYRIGHT VIOLATION
! ---------------------------------------------------------------------------------------
END SUBROUTINE GET_SMODEL
