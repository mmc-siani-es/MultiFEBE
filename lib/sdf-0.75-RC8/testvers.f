      program version
c
c - - test compiler version, and print out flag for 64-bit vs 32-bit
c - - string length hidden arguments
c
      implicit none
      character*100 strlen
      character*80 stringflag
      character*80 fflag_extra
c
# ifndef __GNUC__
# define __GNUC__ 0
# endif
      write(6,*) '__GNUC__ = ',__GNUC__
# ifndef __GFORTRAN__ 
# define __GFORTRAN__ 0
# endif
      write(6,*) '__GFORTRAN__ = ',__GFORTRAN__
# ifndef __INTEL_COMPILER
# define __INTEL_COMPILER  0
# endif
      write(6,*)  '__INTEL_COMPILER = ',__INTEL_COMPILER
      if((__GNUC__ .gt. 0) .and. (__GNUC__ .lt. 7)) then
        strlen='32-bit'
        stringflag='/* #define INTF_CLEN_64 1 */'
      else
        strlen='64-bit'
        stringflag='#define INTF_CLEN_64 1'
      endif
c
      if((__GNUC__ .gt. 0) .and. (__GNUC__ .ge. 10)) then
        fflag_extra = 'FFLAG_EXTRA = -fallow-argument-mismatch'
      else
        fflag_extra = 'FFLAG_EXTRA = '
      endif
c
      write(6,*) 'strlen = ',trim(strlen)
      write(6,*) 'stringflag = ',trim(stringflag)
      write(6,*) trim(fflag_extra)
c
      open(unit=1,file='strlenfmt.h',status='replace')
      write(unit=1,fmt='(a)') trim(stringflag)
      close(unit=1)
      write(6,*) 'file strlenfmt.h replaced'
c
      open(unit=1,file='fflag-extra.inc',status='replace')
      write(unit=1,fmt='(a)') trim(fflag_extra)
      close(unit=1)
      write(6,*) 'file fflag-extra.inc replaced'
     
c
      stop
      end program version
