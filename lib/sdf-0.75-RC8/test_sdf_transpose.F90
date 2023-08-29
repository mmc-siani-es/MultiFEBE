      program test_sdf_transpose
!
! - - Example program reading in an sdf variable that is dynamically 
! - - allocated in 
! - - a fortran 90 program.  Here, it is assumed you know nothing about
! - - the data you want to read in other than its label is called bz in the
! - - sdf file and that it is 
! - - a 3rd rank floating point array.  In this example we read in the 3d array,
! - - then transpose all 3 dimensions, and write the transposed 
! - - array out to a new file.  Note that this transpose and write procedure
! - - would also work just fine if you read
! - - in the array into a long 1-d array of length = product of dimensions.
!
      implicit none
!
! - - match is an integer array, whose size is not yet determined
!
      integer, dimension(:), allocatable :: match
!
! - - rdata is a real array, to be read in from the sdf file, but whose size
! - - has yet to be determined.
!
! - - 1-d version commented out
!
!     real, dimension(:), allocatable :: rdata
!
      real, target, allocatable :: rdata(:,:,:)
!
! - - dims is a statically allocated array that will exceed the expected 
! - - number of dimensions of any variable you might want to read in
!
      integer(kind=8) :: arrsize, dims(20)
      integer :: indorder(20),directions(20)
      integer :: ndat,order,ndim,nbpw,i,ii
      character(len=1) :: dtype
      character(len=100) :: label
!
! - - fname is the name of the sdf file that will be read
!
      character(len=100) :: fnamein
!
! - - fnameout is the name of the output sdf file for the transposed array
!
      character(len=100) :: fnameout
!
      fnamein='bz.sdf'
      fnameout='bzt.sdf'
!
! - - first remove fnameout so it doesn't clog up with too much crap
!
      call sdf_rm_f77(fnameout)
!
! - - First call sdf_query to list the file dataset summary, and to find 
! - - the value of ndat (no. of datasets):
!
      call sdf_query_f77(fnamein,ndat)
      write(6,*) 'ndat = ',ndat
!
! - - now, allocate the match array to be of size ndat
!
      allocate(match(ndat))
!
! - - Now find out which dataset order has the identifying label 'bz' :
!
      call sdf_labmatch_f77(fnamein,ndat,'bz',match)
      order=match(1)
      write(6,*) 'order = ',order
!
! - - OK, now we know the order of the dataset in the file.  Read in the 
! - - details about the bz variable, such as the number of dimensions ndim,
! - - and the actual value of the dimensions (contained in dims array):
!
      call sdf_details_f77(fnamein,order,label,dtype,nbpw,ndim,dims)
!
! - - rdata could be read in as a 1-d array.  It could then be reshaped into its
! - - proper number of dimensions using the reshape function afterwards, if the
! - - user wanted to.  Here, will calculate the size of the 1-d array
! - - big enough to hold bz, but will actually go ahead and allocate it in 3d.
!
      arrsize=dims(1)
      if(ndim > 1) then
         do i=2,ndim
            arrsize=arrsize*dims(i)
         enddo
      endif
      write(6,*) '1-d array size would be ',arrsize
      write(6,*) '3-d input array dimensions = ',(dims(i),i=1,ndim)
! 
! - - (allocating the 1-d array commented out.)
!
!     allocate(rdata(arrsize))
!
      allocate(rdata(dims(1),dims(2),dims(3)))
!
! - - Now that the array is allocated, we can read in bz from the sdf
! - - file:

      call sdf_read_f77(fnamein,order,label,dtype,nbpw,ndim,dims,rdata)
!
! - - Here's where we try and transpose and mess around with the data, in place:
! - - This corresponds with making the z-index the 1st index instead of
! - - the last one, and transposes x and y, with the direction of y
! - - reversed.  This corresponds to a 90 degree clockwise rotation.
! - - indorder is the new desired index order, expressed in terms of the
! - - initial index directions of the input array.  Thus, for example, if
! - - our initial indices correponded to the x,y,z directions, the new
! - - indices in the output array correspond to the z,y,x directions, resp.
! - - The directions array, when its elements are negative corresponds to
! - - reversing the index in the *new* shuffled array.  To reverse the
! - - z-index in the transposed array, in this example, one would therefore set
! - - directions(1) = -1.
!
      indorder(1)=3
      indorder(2)=2
      indorder(3)=1
      directions(1)=1
      directions(2)=-1
      directions(3)=1
!
      call sdf_transpose_f77(indorder,directions,dtype,nbpw,ndim,dims,rdata)
!
! - - Even though the actual dimensions of rdata here in the calling program
! - - have not been changed, the contents of the array have been shuffled
! - - around to be consistent with the indorder and directions array.  
! - - The returned
! - - values in the dims array will reflect the new index order.
!
! - - Unless you are running right at the memory limit, it is probably
! - - faster and more convenient to use the reshape function in Fortran90/95, 
! - - to transpose or re-order an array,
! - - where the transpose can be done by specifying the new index order 
! - - in reshape.  In that case, however, the transpose can't be done
! - - in place, and you'll have to specify a new array to receive the
! - - transpose.
!
! - - Note that the in-place transpose done here
! - - is somewhat clunky to use:  even
! - - though the memory elements have been moved around to be consistent
! - - with the new index order and new values of the dims array, the calling
! - - program itself still thinks it is dimensioned the old way.  This will
! - - not affect what gets written out in the following sdf_write_f77 call,
! - - but you won't be able to use the new indices in this calling program
! - - on the transposed array.  To fix this, you could use an equivalence
! - - statement (like in the f77 example program) but that won't work for
! - - an allocatable array.  There might be a solution to this using fortran
! - - pointers, but all of the examples I've been able to do involve making
! - - a new copy of the array, which negates the advantage of having done
! - - the transpose in-place.  
!
! - - In any case, if you use the same code
! - - fragments as in this program to allocate and read in the transposed
! - - array from file fnameout, it should be read in properly with 
! - - the new order.
!
      call sdf_write_f77(fnameout,"bzt_yx",dtype,nbpw,ndim,dims,rdata)
! 
! - - Now the transposed array is written out to the sdf file.  
!
      deallocate(rdata)
      deallocate(match)
!
! - - done.
!
      stop
      end program test_sdf_transpose
