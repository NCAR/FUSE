program read_para_file

	character(length=100), intent(in)	:: fileName			!
	real(dp)							:: paramValues		! array for the parameters
	integer(kind=4)						:: fileId=93			! id for the input file


	open(fileId,file=trim(fileName))

	close(fileId)

	include

	subroutine 


end program read_para_file