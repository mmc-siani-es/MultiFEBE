      program testt
c
c - - To illustrate the usage of sdf_write_f77 and sdf_transpose_f77 to allow 
c - - one to
c - - order the output data such that all physical variables (rho,p,v)
c - - are written out for
c - - each gridpoint before the next gridpoint's data is written.  
c - - For illustration purposes, rho is set to 1, p to 2, vx to 3, vy to 4
c - - and vz to 5.  This will allow us to see if the variables get all
c - - scrambled up if the transpose is unsuccessful.

      implicit none

      real rho(100,100,100), p(100,100,100),v(100,100,100,3)
c
c - - pvart is dimensioned for the transpose of pvar
c
      real pvart(5,100,100,100), pvar(100,100,100,5)
c
c - - These next equivalence statements are the essential thing that lets us
c - - reorder the data in this program.  pvar and pvart occupy same memory
c - - space, but are indexed differently.  pvar is also equivalenced to the
c - - physical variables rho, p, and v.
c
      equivalence (pvar,pvart)
      equivalence (rho, pvar(1,1,1,1))
      equivalence (p,pvar(1,1,1,2))
      equivalence (v, pvar(1,1,1,3))
      integer nbpw, ndim, indorder(5),reverse(5)
      integer*8 dims(20)
      integer i,j,k
      character*1 dtype
      character*80 label
      character*80 fname
c
      do k=1,100
         do j=1,100
            do i=1,100
               rho(i,j,k)=1.0
               p(i,j,k)=2.0
               v(i,j,k,1)=3.0
               v(i,j,k,2)=4.0
               v(i,j,k,3)=5.0
            enddo
         enddo
      enddo
c
c - - Now, the values of the physical variables in their arrays have been
c - - set.  We want the pvar array to be ordered such that for each gridpoint,
c - - all the variables appear in order, with that order repeated when we
c - - move to the next gridpoint.  Even though pvart is dimensioned
c - - this way, this isn't the way the data is currently stored in memory:
c - - instead, rho is stored, then p, and then v.  To reorder the array, will
c - - use sdf_transpose_f77.  First, set up the indorder array to tell
c - - sdf_transpose_f77 what order the new indices should run in terms of
c - - the old indices:
c
      nbpw=4
      ndim=4
      dims(1)=100
      dims(2)=100
      dims(3)=100
      dims(4)=5
c
      indorder(1)=4
      indorder(2)=1
      indorder(3)=2
      indorder(4)=3
c
      reverse(1)=1
      reverse(2)=1
      reverse(3)=1
      reverse(4)=1
c
      label="pvar"
      dtype="f"
      fname="junkf77_2.sdf"
c
c - - OK, do the transpose of pvar:
c
c
      call sdf_transpose_f77(indorder,reverse,dtype,nbpw,ndim,dims,
     1 pvar)
c
c - - remove any pre-existing copies of file:
c
      call sdf_rm_f77(fname)
c
c - - Now the data should be arranged correctly. We're ready to write out 
c - - the 4d array to an sdf file:
c
      call sdf_write_f77(fname,label,dtype,nbpw,ndim,dims,pvar)
c
c - - Now, print out pvart(1-5,50,50,50) to show that it really is in
c - - the right order:
c
      write(6,*) "pvart(1,50,50,50) = ",pvart(1,50,50,50)
      write(6,*) "pvart(2,50,50,50) = ",pvart(2,50,50,50)
      write(6,*) "pvart(3,50,50,50) = ",pvart(3,50,50,50)
      write(6,*) "pvart(4,50,50,50) = ",pvart(4,50,50,50)
      write(6,*) "pvart(5,50,50,50) = ",pvart(5,50,50,50)
c
c - - Now transpose back to original order:
c
      indorder(1)=2
      indorder(2)=3
      indorder(3)=4
      indorder(4)=1
c
      call sdf_transpose_f77(indorder,reverse,dtype,nbpw,ndim,dims,
     1 pvar)
c
c - - Show that all the physical variables still retain correct
c - - values after being transposed back to orig. order:
c
      write(6,*) "rho(50,50,50) = ",rho(50,50,50)
      write(6,*) "p(50,50,50) = ", p(50,50,50)
      write(6,*) "v(50,50,50,1) = ", v(50,50,50,1)
      write(6,*) "v(50,50,50,2) = ", v(50,50,50,2)
      write(6,*) "v(50,50,50,3) = ", v(50,50,50,3)
c
      stop
      end
