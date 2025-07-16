      program testwbrb
c
c - - Test the ability to create platform-independent large-endian binary
c - - data files, write to them, read from them.
c
c - - Note that for some reason, this test program runs *extremely* slowly
c - - on MS windows platforms.  From doing some experiments, this appears to
c - - be due to (1) very slow opening and closing of files on MS Windows,
c - - at least using ANSI-like fopen/fclose statements, and (2) writing out
c - - small amounts of data at a time.  From the other test programs, 
c - - in which MS windows performs similar to to linux, it appears that if
c - - the amount of data written in each write statement exceeds some 
c - - threshold, that the performance is reasonably good.  I think this is
c - - the only test program in which so many opens/closes are done and where
c - - the writes are so small.
c
      implicit none
c
c - - datapos, nelem must be 64-bit integers
c
      integer*8 datapos
      integer*8 nelemc,nelemf,nelemd,nelemi
c
      integer nbpwc,nbpwf,nbpwd,nbpwi,nrep,i
      character*80 fname
c
c - - variables to be written
c
      character*40 startlab
      real sp(100)
      double precision dp(100)
      integer isp(100)
c
c - - variables to be read
c 
      character*40 endlab
      real spr(100)
      double precision dpr(100)
      integer ispr(100)
c
c - - file descriptor:
c
      integer*8 fp
c
c - - Setting nrep = 3000000 will create a 4.8 GB file.  Adjust downward
c - - if you don't want a file that big!
c
      nrep=100000
      fname="testwbrb.dat"
      datapos=0
      nbpwc=1
      nbpwf=4
      nbpwd=8
      nbpwi=4
      nelemc=40
      nelemf=100
      nelemd=100
      nelemi=100
c
c - - create data to be written
c
      startlab="Initial string in file"
      do i=1,100 
         isp(i)=i
         sp(i)=sqrt(float(i))
         dp(i)=dsqrt(dble(i))
      enddo
c
c - - to test letting sdf_wb_f77 do its own opening and closing, set
c - - fp to zero.  To open and close yourself, use fp returned from 
c - - fopen_f77().
c
      fp=0
c
c - - uncomment next line to do your own open/close
c     call fopen_f77(fname,"wb",fp)
      write(6,*) "fp = ",fp

c
c - - write the file:
c
      call sdf_wb_f77(fname,datapos,nelemc,nbpwc,fp,startlab)
c
      do i=1,nrep
         call sdf_wb_f77(fname,datapos,nelemf,nbpwf,fp,sp)
         call sdf_wb_f77(fname,datapos,nelemi,nbpwi,fp,isp)
         call sdf_wb_f77(fname,datapos,nelemd,nbpwd,fp,dp)
      enddo
c
c - - uncomment next line to do your own open/close
c     call fclose_f77(fp)
c
c
c - - read the file:
c
      datapos=0
c
c
c - - uncomment next line to do your own open/close
c     call fopen_f77(fname,"rb",fp)
      write(6,*) "fp = ",fp
      call sdf_rb_f77(fname,datapos,nelemc,nbpwc,fp,endlab)
      do i=1,nrep
         call sdf_rb_f77(fname,datapos,nelemf,nbpwf,fp,spr)
         call sdf_rb_f77(fname,datapos,nelemi,nbpwi,fp,ispr)
         call sdf_rb_f77(fname,datapos,nelemd,nbpwd,fp,dpr)
      enddo
c
c - - uncomment next line to do your own open/close
c     call fclose_f77(fp)
c
c - - print out values that are read in from last read repetition:
c
      write(6,*) "endlab = ",endlab
      write(6,*) "spr = ", spr
      write(6,*) "ispr = ", ispr
      write(6,*) "dpr = ", dpr
      write(6,*) "final value of datapos = ", datapos
c
      stop
      end
