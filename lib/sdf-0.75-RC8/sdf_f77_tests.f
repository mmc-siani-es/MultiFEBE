      program testf77
c
c - - This is a simple fortran 77 program illustrating the use of sdf
c - - functions to read and write binary data to disk.  Also illustrates
c - - simple uses of in-place transposes using
c - - sdf_transpose_f77 on some 2d arrays.  The subroutine sdf_tranpose_f77
c - - can also perform in-place transposes for an arbitrary number of 
c - - dimensions.  Note the use of an equivalence statement so that the
c - - transposed array can be indexed correctly.
c
      implicit none
      integer ndim,fndim,i,j,ic,jc,mult,iorder
      real pi,sigma,xd,yd
      character*80 fname
      character*80 label
      character*80 flabel
      character*1 dt
      character*1 fdt
      integer*8 dims(10)
      integer*8 fdims(10)
      integer indorder(10),directions(10)
      integer nbpw
      integer fnbpw
      integer ndat
      real realdata(100,50)
      real fdata(100,50)
      real fdatat(50,500)
      complex cmpdata(100,50)
      integer*8 dsizes(100)
      integer dnbpw(100)
      character*40 chardata
      character*40 fchardata
      character*1 datatype
      character dtypes(100)*1
c
c - - fdata and fdatat are equivalenced to facilitate the different indexing
c - - of fdata that will occur after the in-place rotation/transpose:
c
      equivalence(fdata,fdatat)
c
      ndim=2
      nbpw=4
      dims(1)=100
      dims(2)=50
      pi=3.141592653
      dt='f'
c
c - - note that the blanks in label will be removed by the sdf software:
c
      label='test real 2d'
c
c - - fname is the name of the sdf file that will be used for the sdf i/o
c - - tests.
c
      fname="junkf77.sdf"
c
c - - Create some character data to output later.
c
      chardata='barf, barf, barf'
      ic=75
      jc=15
      sigma=10.
c
c - - compute a 2-d real array (an offset gaussian), and a 2-d complex array
c - - in which the real part looks like sin(x) and the imaginary part looks
c - - like cos(x):
c
      do i=1,100
         do j=1,50
           xd=float(i-ic)
           yd=float(j-jc)
           realdata(i,j)=exp(-(xd/sigma)**2-(yd/sigma)**2)
           cmpdata(i,j)=(1.,0.)*sin(pi*float(i)/float(100))
     1                 +(0.,1.)*cos(pi*float(i)/float(100))
         enddo
      enddo
c
c - - remove the sdf file fname, if it exists
c
      call sdf_rm_f77(fname)
c
c - - write out the 2-d real array to the sdf file
c
c - - here, dt = 'f' indicates a floating point array, nbpw=4 indicates
c - - no. of bytes per word, ndim = 2 indicates a 2-d array, and 
c - - (dims(1)=100 and dims(2)=50 indicate the dimension values for realdata.
c
c
      call sdf_write_f77(fname,label,dt,nbpw,ndim,dims,realdata)
c
c - - change dimensions to a 1-d array (of characters)
c
      dims(1)=40
      dims(2)=0
c
c - - write out the character data.  Note that the blank in the
c - - label field below will be removed by the sdf software.
c - - Here, use explicit values for dt, nbpw, ndim rather than variable names.
c
      call sdf_write_f77(fname,"char data","b",1,1,dims,chardata)
c
c - - now change the dimensions back to 2-d values for the complex array
c
      dims(1)=100
      dims(2)=50
c
c - - write out the complex 2-d array to the file.  Note the "c" value for dt.
c
      call sdf_write_f77(fname,"cmplxdat","c",4,2,dims,cmpdata)
c
c - - call sdf_query to print out the metadata of the datasets in the
c - - file
c
      call sdf_ndat_f77(fname,ndat)
      write(6,*) "ndat (from sdf_ndat_f77) = ",ndat
      call sdf_query_f77(fname,ndat)
      write(6,*) 'ndat (from sdf_query_f77) = ',ndat
      iorder=1
c
c - - with iorder set to 1, the next line reads in the 2nd dataset from the
c - - sdf file into the fchardata character string.  (datasets are numbered
c - - from 0, not 1):
c
      call sdf_read_f77(fname,iorder,flabel,fdt,fnbpw,fndim,fdims,
     1  fchardata)
c
c - - print out summaries of the properties of this dataset, along with
c - - what was read in
c
      write(6,*) "flabel = ", flabel
      write(6,*) "fdt = ", fdt
      write(6,*) "fnbpw = ", fnbpw
      write(6,*) "fndim = ",fndim
      write(6,*) "fdims = ",(fdims(i),i=1,fndim)
      write(6,*) "fdata = ",fchardata
c
      iorder=0
c 
c - - The next line reads in the 0th (ie 1st) dataset in the file into
c - - the array fdata.  This is the 2-D offset gaussian.
c
      call sdf_read_f77(fname,iorder,flabel,fdt,fnbpw,fndim,fdims,fdata)
c
c - - print out summaries of the metadata for this dataset, (and then print
c - - out all the data if you uncomment the last statement in this section)
c
      write(6,*) "flabel = ", flabel
      write(6,*) "fdt = ", fdt
      write(6,*) "fnbpw = ", fnbpw
      write(6,*) "fndim = ",fndim
      write(6,*) "fdims = ",(fdims(i),i=1,ndim)
c     write(6,*) "fdata = ",((fdata(i,j),j=1,fdims(2)),i=1,fdims(1))
c
c - - write out one row of fdata (for later comparison w/ transpose/rotation)
c
      write(6,*) "fdata row 15 = ",(fdata(i,15),i=1,100)

c 
c - - Now, write out one more dataset that is a transpose of the offset
c - - gaussian dataset.  To do this, use the in-place transpose
c - - function sdf_transpose_f77 - we will also flip the sign of new
c - - first index.  The end result of all of this should be a 90 degree
c - - rotation counter-clockwise.
c
c - - The next 2 statements interchange the roles of x and y in the tranpose:
c
      indorder(1)=2
      indorder(2)=1
c
c - - setting directions(1) = -1, in addition to the transpose,
c - - will yield a 90 degree counterclockwise
c - - rotation of the original array.  Note that the elements of directions
c - - refer to the order of dimensions *after* the transpose is done:
c
      directions(1)=-1
      directions(2)=1
c
c - - Now do the transpose in place.  Note that even though the dimensions
c - - of the array fdata in the calling program are not changed, the array 
c - - elements will have been shuffled around
c - - in place such that when written out to disk, it will have the 
c - - correct order for the transposed array.  The elements of the fdims
c - - array will have been changed on output to reflect the new dimensions,
c - - even though the dimensions of fdata are not changed in the main program.
c
      call sdf_transpose_f77(indorder,directions,fdt,fnbpw,fndim,
     $ fdims,fdata)
c
c - - print out 1 column of transpose/rot to compare with row of original data:
c - - (use equivalenced array fdatat for transposed indexing):
c
      write(6,*) " Column 36 of rotated floating array = ",
     $ (fdatat(36,i),i=1,100)
c
c - - Now the transpose has been done, write out the transposed 
c - - array to the sdf file.  Could use either fdata or fdatat on the write
c - - statement since they're equivalenced.
c
      call sdf_write_f77(fname,"testreal2d_yx",fdt,fnbpw,fndim,
     $ fdims,fdata)
c 
c - - At this point, if one wanted to, one could transpose fdata
c - - back to its original form by setting indorder(1)=2 and indorder(2)=1,
c - - and by setting directions(1)=1 and directions(2)=-1 , and calling
c - - sdf_transpose again, with the same arguments as above.
c
c
c - - Now transpose the complex data in place and write it out:
c
      indorder(1)=2
      indorder(2)=1
      dims(1)=100
      dims(2)=50
      call sdf_transpose_f77(indorder,directions,"c",4,2,dims,cmpdata)
c
c - - In-place transpose of complex data has been done.  Now write it out:
c
      call sdf_write_f77(fname,"cmplxdat_yx","c",4,2,dims,cmpdata)
c
c
c - - The next call determines the sizes of all the datasets in the file:
c
      call sdf_sizes_f77(fname,ndat,dtypes,dsizes,dnbpw)
c 
c - - print out all the size and type information for all the datasets
c
      write(6,*) "ndat = ",ndat
      write(6,*) "dtypes = ",(dtypes(i),i=1,ndat)
      write(6,*) "dsizes = ",(dsizes(i),i=1,ndat)
      write(6,*) "dnbpw = ",(dnbpw(i),i=1,ndat)

      do i=1,ndat
         mult=1
         if(dtypes(i) .eq. "c") mult=2
         write(6,*) "dataset:",i-1," type: ",dtypes(i), 
     1 " bytesize = ",dsizes(i)*dnbpw(i)*mult
      enddo
      stop
      end
