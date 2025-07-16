      program test_prim
c
c - - Test fortran primitive calls of C i/o ,byteswap, is_big_endian functions.
c - - ibe_f77 tests for big endian, byteswap_f77 does bytwswapping, 
c - - fopen_f77 provides an interface to the C fopen function, 
c - - fclose_f77 provides an interface to the C close function, 
c - - fwrite_f77 provides an interface to the C fwrite function,
c - - fread_f77 provides an interface to the C fread function, and
c - - fseeko_f77 provides an interface to the C fseeko function, and
c - - ftello_f77 provides an interface to the C ftello function, and
c - - sdf_ftc_f77 provides an interface to the SDF function sdf_file_truncate.
c
      implicit none
c
c - - fp is the file pointer, should be integer*8 to ensure it works on
c - - both 32-bit and 64-bit machines, as it must be big enough to
c - - hold an address.
c
      integer*8 fp
      real a
      double precision b
      integer c
      integer*8 d
c
      integer*8 lengthf
      integer*8 filesize
      integer*8 filepos
      integer nbpwf,nbpwd,nbpwi,nbpwi8
      integer ibe,ise
      character*80 fname
c
      fp=0
      fname="testprim.dat"
      lengthf=1
      nbpwf=4
      nbpwd=8
      nbpwi=4
      nbpwi8=8
      a=5.
      b=dble(11)
      c=4829636
      d=c
      d=d**2
c
c - - determine whether we are on big or small endian machine
c
      ise=1
      call ibe_f77(ibe)
      if(ibe .eq. 1) ise=0
      if(ise .eq. 1) write(6,*) "small endian: byteswapping needed"
      if(ise .eq. 0) write(6,*) "large endian: no byteswapping needed"
c
c - - byteswap before write if necessary
c
      if(ise .eq. 1) then
         call byteswap_f77(a,lengthf,nbpwf)
         call byteswap_f77(b,lengthf,nbpwd)
         call byteswap_f77(c,lengthf,nbpwi)
         call byteswap_f77(d,lengthf,nbpwi8)
      endif
c
c - - write out data
c
      call fopen_f77(fname,"wb+",fp)
c 
c - - Just for heck of it, print out file pointer (will be a big integer)
c
      write(6,*) "fp = ",fp
c
      call fwrite_f77(a,nbpwf,lengthf,fp)
      call fwrite_f77(b,nbpwd,lengthf,fp)
      call fwrite_f77(c,nbpwi,lengthf,fp)
      call fwrite_f77(d,nbpwi8,lengthf,fp)
c
c - - flush the stream fp.  This is not necessary.  It is just to 
c - - test that the function works ok.
c
      call fflush_f77(fp)
c
c - - truncate file and close
c
      call sdf_ftc_f77(fp,lengthf*(nbpwf+nbpwd+nbpwi+nbpwi8))
      call fclose_f77(fp)
c
c - - byteswap back if necessary, if calculations will then be done with
c - - a,b,c,d.  Otherwise, values of these variables will be completely
c - - messed up
c
      if(ise .eq. 1) then
         call byteswap_f77(a,lengthf,nbpwf)
         call byteswap_f77(b,lengthf,nbpwd)
         call byteswap_f77(c,lengthf,nbpwi)
         call byteswap_f77(d,lengthf,nbpwi8)
      endif
c
      a=0.
      b=dble(0.)
      c=0
      d=0
c
      call fopen_f77(fname,"rb",fp)
      call fread_f77(a,nbpwf,lengthf,fp)
      call fread_f77(b,nbpwd,lengthf,fp)
      call fread_f77(c,nbpwi,lengthf,fp)
      call fread_f77(d,nbpwi8,lengthf,fp)
c
c - - get file size:
c
      call file_sz_f77(fname,filesize)
      write(6,*) "filesize = ",filesize
c
      call fclose_f77(fp)
c
c - - byteswap if necessary, after reading in data
c
      if(ise .eq. 1) then
         call byteswap_f77(a,lengthf,nbpwf)
         call byteswap_f77(b,lengthf,nbpwd)
         call byteswap_f77(c,lengthf,nbpwi)
         call byteswap_f77(d,lengthf,nbpwi8)
      endif
      write(6,*) "a = ",a
      write(6,*) "b = ",b
      write(6,*) "c = ",c
      write(6,*) "d = ",d
c
c - - next call just flushes stdout.
c
      call ffstdout_f77()
c
c
c - - Go back and just read d, by seeking to the right place:
c
      call fopen_f77(fname,"rb",fp)
      call fseeko_f77(fp,filesize-8)
      call fread_f77(d,nbpwi8,lengthf,fp)
      if(ise .eq. 1) call byteswap_f77(d,lengthf,nbpwi8)
      call ftello_f77(fp,filepos)
      write(6,*) "filepos = ",filepos
      write(6,*) "d (again) = ",d
      call fclose_f77(fp)
      call ffstdout_f77()
      stop
      end
