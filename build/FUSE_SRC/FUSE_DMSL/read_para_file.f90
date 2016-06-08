program read_para_file

	implicit none

	character(length=100), intent(in)	:: fileName			! name of the input file
	integer(kind=4), intent(in)			:: numParam			! number of parameters
	real(dp), intent(out)				:: paramValues		! array for the parameters
	integer(kind=4)						:: fileId=93		! id for the input file

	include

	subroutine read_para_file(fileName,paramValues,paramValues)

		open(fileId,file=trim(fileName))

		read(*,fileId)

		close(fileId)

	end subroutine read_para_file

end program read_para_file