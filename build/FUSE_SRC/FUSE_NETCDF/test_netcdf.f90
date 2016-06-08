program test_netcdf
! test program for D
include 'netcdf.inc'                                  ! use netCDF libraries

! error code
integer  :: ierr

! NetCDF file id
integer  :: ncid

! NetCDF dimension ids
integer  :: ix_dim,iy_dim,it_dim

! NetCDF variable IDs
integer  :: ivar_id

! dimension lengths
integer  :: nx,ny

! character strings
integer  :: str_len=64
character(len=str_len) :: filename  ! filename
character(len=str_len) :: var_name  ! variable name
character(len=str_len) :: var_desc  ! variable description

! -------------- end of definitions ----------------------------------

! define some stuff

nx = 2  ! length of the x dimension
ny = 3  ! length of the y dimension

filename = 'D.nc'
var_name = 'stuff'
var_desc = 'just some random crap'

! create a file

! create file
ierr = nf_create(trim(filename),nf_clobber,ncid); call handle_err(ierr)

 ! define dimensions
 ierr = nf_def_dim(ncid,'x',nx,ix_dim); call handle_err(ierr)
 ierr = nf_def_dim(ncid,'y',ny,iy_dim); call handle_err(ierr)
 ierr = nf_def_dim(ncid,'t',nf_unlimited,it_dim); call handle_err(ierr)

 ! define a variable
 ierr = nf_def_var(ncid,'stuff',nf_real,3,(/ix_dim,iy_dim,it_dim/),ivar_id); call handle_err(ierr)
 ierr = nf_put_att_text(ncid,ivar_id,'long_name',len_trim(var_desc),trim(var_desc)); call handle_err(ierr)

! end definitions and close file
ierr = nf_enddef(ncid); call handle_err(ierr)
ierr = nf_close(ncid); call handle_err(ierr)

stop
end
