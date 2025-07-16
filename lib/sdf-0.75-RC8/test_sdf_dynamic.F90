      program test_dynamic_allocate
!
! - - Example program reading in an sdf variable dynamically in 
! - - a fortran 90 program.  Here, it is assumed you know nothing about
! - - the data you want to read in other than it is called ranxy and it is 
! - - a 2nd rank floating point array.
!
      implicit none
!
! - - match is an integer array, whose size is not yet determined
!
      integer, dimension(:), allocatable :: match
!
! - - rdata is a real array, to be read in from the sdf file, but whose size
! - - has yet to be determined.  Here, we'll read it in as a 1-d array.
!
      real, dimension(:,:), allocatable :: rdata
!
! - - dims is a statically allocated array that will exceed the expected 
! - - number of dimensions of any variable you might want to read in
!
      integer(kind=8) :: dims(20), arrsize
      integer :: ndat,order,ndim,nbpw,i,ii
      character(len=1) :: dtype
      character(len=100) :: label
!
! - - fname is the name of the sdf file that will be read
!
      character(len=100) :: fname
      fname="sdf_example_file.sdf"
!
! - - First call sdf_query to list the file dataset summary, and to find 
! - - the value of ndat, which is the number of datasets in the file:
!
      call sdf_query_f77(fname,ndat)
      write(6,*) 'ndat = ',ndat
!
! - - now, allocate the match array to be of size ndat
!
      allocate(match(ndat))
!
! - - Now find out which dataset order has the identifying label 'ranxy' :
!
      call sdf_labmatch_f77(fname,ndat,'ranxy',match)
      order=match(1)
      write(6,*) 'order = ',order
!
! - - OK, now we know the order of the dataset in the file.  Read in the 
! - - details about the ranxy variable, such as the number of dimensions ndim,
! - - and the actual value of the dimensions (contained in dims array):
!
      call sdf_details_f77(fname,order,label,dtype,nbpw,ndim,dims)
!
! - - ranxy will be read in as a 2-d array.  Could also be read in as 1d array.
!     Just for kicks, calculate the size of the 1-d array
! - - big enough to hold ranxy:
!
      arrsize=dims(1)
      if(ndim > 1) then
         do i=2,ndim
            arrsize=arrsize*dims(i)
         enddo
      endif
      write(6,*) 'arrsize = ',arrsize
! 
! - - Go ahead and allocate the array.
!
      allocate(rdata(dims(1),dims(2)))
!
! - - Now that the array is allocated, we can read in ranxy from the sdf
! - - file:

      call sdf_read_f77(fname,order,label,dtype,nbpw,ndim,dims,rdata)
!
! - - print out 1st row of elements of ranxy
!
      write(6,*) 'ranxy(1,1:dims(2)) = ',rdata(1,1:dims(2))
!
      deallocate(rdata)
      deallocate(match)
!
! - - we're done.  Note that, if you already know the dataset order, and
! - - you pre-dimension ranxy, the only thing you'd need to do to read in
! - - the data would be the single sdf_read_f77 call above!
!
      stop
      end program test_dynamic_allocate
