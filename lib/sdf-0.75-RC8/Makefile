#
# To compile SDF library and sdf_browse executable, type "make".
# To compile SDF library, sdf_browse executable, and the test programs,
# type "make all".  The test programs all begin with the letter x.
# To install the libraries, sdf_browse, and the include
# file to the location denoted by PREFIX, type "make install".  To remove
# the installation, type "make uninstall".  Depending on permissions, you
# may need to be root to do "make install" and "make uninstall".  To
# remove all of the object files, executables, and output from test programs
# from the build directory, type "make clean".  To compile all of the test
# programs, you will need to define the Fortran compilers F77 and F95, but
# if you are only compiling the library and sdf_browse, you need only define
# the C-compiler COMPILER .  
# To compile the dynamic link libraries sdf.dll and sdf_if64.dll, 
# and definitions files sdf.def and sdf_if64.def,
# type "make sdf_dll.zip".  This will only work under msys/mingw 
# in MS Windows.  It
# will also create a zip file, sdf_dll.zip which will contain the .dll and .def
# files. 
#
# This Makefile is set to use gcc as the default
# C-compiler, gfortran as the default F77 compiler, and gfortran as the default
# F95 compiler.  The compiler flags FLAGS is set for optimization and no
# debugging.  The Fortran compiler flags F77FLAGS, F95FLAGS are set similarly.
# If you want to change the compilers or options for flags, delete the
# existing defaults and uncomment one of the other existing choices, depending
# on your platform and desired compilers and flags.  You can also add new
# choices if needed.
#
# Note:  If you change the default fortran integer size to 8 bytes,
# set the F77_IF64, F95_IF64 flags properly for your choices of fortran 
# compiler -- see examples later in this file.

# install directory: - can be changed to whatever you want
 PREFIX = /usr/local
# PREFIX = /usr/local/sdf-0.73

# object files
SDF_OBJECTS = sdf_subs.o

# C compiler (required) - uncomment the one you want
# COMPILER = gcc
COMPILER = gcc
# COMPILER = gcc-mp-10
# COMPILER = gcc4
# COMPILER = icc
# COMPILER = /opt/intel_cc_80/bin/icc
# COMPILER = cc
# COMPILER = pgcc
# COMPILER = xlC_r

# Fortran 77 compiler (optional if compiling only library and sdf_browse; 
# but required if compiling test programs) - uncomment the one you want

# F77 = g77
# F77 = g95
# F77 = gfortran
# F77 = gfortran-mp-10
F77 = gfortran
# F77 = ifort
# F77 = f77
# F77 = lf95
# F77 = pgf77
# F77 = pgf90
# F77 = xlf_r

# Fortran 90/95 compiler (optional if compiling only library and sdf_browse;
# but required for compiling the test programs) - uncomment the one you want

# F95 = g95
# F95 = gfortran
 F95 = gfortran
# F95 = ifort
# F95 = lf95
# F95 = pgf90
# F95 = f90
# F95 = xlf95_r -qsuffix=f=F90
#
# Now determine conditions for compiling SDF that depend on Fortran
# compiler details, by running a simple fortran code:

# ./xforttest will generate two small text files, fflag-extra.inc and 
# lenstringfmt.h which will facilitate the compilation of the SDF 
# library and the SDF test programs
#

# Flags for the C compiler:
# make sure flags are set for large file support
# UNCOMMENT ONE OF THESE SET OF COMPILER FLAGS OPTIONS

# debug version of flags, with profile support, and -Wall
#FLAGS = -g -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -I. -pg

# debug version of flags, with profile support
#FLAGS = -g -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -I. -pg

# debug version of flags, without  profile support
#FLAGS = -g -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -I.

# un-optimized, non-debug (default options) version of flags, with profiling
#FLAGS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -I. -pg

# un-optimized, non-debug (default options) version of flags, no profiling
#FLAGS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -I. 

# optimized version of flags, with profile support
#FLAGS = -O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -I. -pg

# optimized version of flags, without profile support
#FLAGS = -O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -I. 

# optimized version of flags, including -fomit-frame-pointer 
# also includes -fPIC (position independent code) - necessary for linking
# into C++ code
# FLAGS = -fPIC -fomit-frame-pointer -O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -I. 
# FLAGS = -fPIC -fomit-frame-pointer -O3 -INTF_CLEN_64 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -I. 
FLAGS = -fPIC -fomit-frame-pointer -O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -I. 

# optimized flags for SUN cc compiler, without profile support
#FLAGS = -O -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -I. 


# Flags for Fortran Compilers:  Pick one for each compiler and uncomment it.
#
# The following include file is designed for the GNU fortran compiler gfortran,
# and if __GNUC__ exceeds 10, will include the flag -fallow-argument-mismatch .
# include ./fflag-extra.inc
  INCLUDEFILE = fflag-extra.inc
  include $(INCLUDEFILE)

# Flags for optimization, no profiling
# F77FLAGS = -O3
# following is F77FLAGS for gfortran:
 F77FLAGS = -O3 $(FFLAG_EXTRA)
# F77FLAGS = -fallow-argument-mismatch -O3
# F95FLAGS = -O3
# following is F95FLAGS for gfortran:
 F95FLAGS = -O3 $(FFLAG_EXTRA)
#  F95FLAGS = -fallow-argument-mismatch -O3 

# Flags for optimization, and profiling
# F77FLAGS = -O3 -pg $(FFLAG_EXTRA)
# F95FLAGS = -O3 -pg $(FFLAG_EXTRA)

# Flags for debugging, no profiling
# F77FLAGS = -g $(FFLAG_EXTRA)
# F95FLAGS = -g $(FFLAG_EXTRA)

# Flags for debugging, and profiling
# F77FLAGS = -g -pg $(FFLAG_EXTRA)
# F95FLAGS = -g -pg $(FFLAG_EXTRA)

# Flags for default options, no profiling
# F77FLAGS = $(FFLAG_EXTRA) 
# F95FLAGS = $(FFLAG_EXTRA) 

# Flags for default options, and profiling
# F77FLAGS = -pg $(FFLAG_EXTRA)
# F95FLAGS = -pg $(FFLAG_EXTRA)

# Flags for optimization, no profiling, for IBM SPX 32-bit systems:
# F77FLAGS = -qarch=auto -O3  -bmaxdata:2000000000 -bmaxstack:25600000
# F95FLAGS = -qarch=auto -O3  -bmaxdata:2000000000 -bmaxstack:25600000

# Flags for optimization, no profiling, for IBM SPX 64-bit systems:
# F77FLAGS = -qarch=auto -O3 
# F95FLAGS = -qarch=auto -O3 

# NOTE - if you pick a version with profiling, you might need to add -lgprof
# to -lsdf (in MINGW, for example)

#
# F77, F95 Flags for setting default fortran integer size to 64-bit integers:
#

# This for ifort, g77, g95:
# F77_IF64 = -i8
# This for gfortran:
# F77_IF64 = -fdefault-integer-8
 F77_IF64 = -fdefault-integer-8

# This for ifort, g95:
# F95_IF64 = -i8
# This for gfortran:
# F95_IF64 = -fdefault-integer-8
 F95_IF64 = -fdefault-integer-8

sdf_browse : sdf_browse.c libsdf.a libsdf_if64.a
	echo "compiling sdf_browse" 
	$(COMPILER) $(FLAGS) sdf_browse.c -L. -lsdf -o sdf_browse

fflag-extra.inc : testvers.f
	echo "compiling xforttest"
	$(F95) -cpp testvers.f -o xforttest
	./xforttest

xforttest : testvers.f
	echo "compiling xforttest"
	$(F95) -cpp testvers.f -o xforttest
	./xforttest

libsdf.a : sdf_subs.o
	echo "making library archive" $@
	ar -rs libsdf.a sdf_subs.o
#	Intel compiler commands for creating the static library: (uncomment)
#	xiar -rsv libsdf.a sdf_subs.o

libsdf_if64.a : sdf_subs_if64.o
	echo "making library archive" $@
	ar -rs libsdf_if64.a sdf_subs_if64.o

sdf_subs.o : sdf_subs.c sdf_subs.h strlenfmt.h
#	echo "running xforttest"
#	./xforttest
#	echo "including fflag-extra.inc"
#	include $(INCLUDEFILE)
	echo "compiling sdf_subs.o" $@ 
	$(COMPILER) $(FLAGS) -c sdf_subs.c

sdf_subs_if64.o : sdf_subs.c sdf_subs.h
	echo "compiling sdf_subs.o" $@ 
	$(COMPILER) -DINTF_64 $(FLAGS) -c sdf_subs.c -o sdf_subs_if64.o

warnings : sdf_subs.c sdf_subs.h sdf_browse.c
	echo "generate all compiler warnings for sdf_subs.c:"
	$(COMPILER) -Wall $(FLAGS) -c sdf_subs.c
	$(COMPILER) -Wall $(FLAGS) -c sdf_browse.c

my_C_prog: libsdf.a my_C_prog.c sdf_subs.h
	echo "compiling my_C_prog" $@
	$(COMPILER) $(FLAGS) my_C_prog.c -I. -L. -lsdf -o my_C_prog

# make sdf.dll works using mingw/msys in MS Windows only:
sdf.dll: sdf_subs.c sdf_subs.h
	echo "compiling sdf.dll" $@ 
	$(COMPILER) $(FLAGS) -shared -o sdf.dll sdf_subs.c -Wl,--output-def,sdf.def,--out-implib,sdfsubs_dll_imp.a

# make sdf_if64.dll works using mingw/msys in MS Windows only:
sdf_if64.dll: sdf_subs.c sdf_subs.h
	echo "compiling sdf_if64.dll" $@ 
	$(COMPILER) -DINTF_64 $(FLAGS) -shared -o sdf_if64.dll sdf_subs.c -Wl,--output-def,sdf_if64.def,--out-implib,sdfsubs_if64_dll_imp.a

# make sdf_dll.zip (works using mingw/msys in MS Windows only):
sdf_dll.zip:  sdf.dll sdf_if64.dll 
	echo "creating zip file w/ sdf.dll, sdf_if64.dll sdf.def sdf_if64.def"
	zip sdf_dll sdf.def sdf.dll sdf_if64.def sdf_if64.dll

xtestmain: libsdf.a main_test_sdf.c sdf_subs.h
	echo "compiling xtestmain" $@
	$(COMPILER) $(FLAGS) main_test_sdf.c -I. -L. -lsdf  -o xtestmain

xtest_sdf_edit: libsdf.a sdf_test_edit.c sdf_subs.h
	echo "compiling xtest_sdf_edit" $@
	$(COMPILER) $(FLAGS) sdf_test_edit.c -I. -L. -lsdf -lm -o xtest_sdf_edit

xlarge: libsdf.a test_sdf_large.c sdf_subs.h
	echo "compiling xlarge" $@
	$(COMPILER) $(FLAGS) test_sdf_large.c -I. -L. -lsdf -lm -o xlarge

x3d: libsdf.a test_sdf_3d.c sdf_subs.h
	echo "compiling x3d" $@
	$(COMPILER) $(FLAGS) test_sdf_3d.c -I. -L. -lsdf -lm -o x3d

x4d: libsdf.a test_sdf_4d.c sdf_subs.h
	echo "compiling x4d" $@
	$(COMPILER) $(FLAGS) test_sdf_4d.c -I. -L. -lsdf -lm -o x4d

x5d: libsdf.a test_sdf_5d.c sdf_subs.h
	echo "compiling x4d" $@
	$(COMPILER) $(FLAGS) test_sdf_5d.c -I. -L. -lsdf -lm -o x5d

xtc89: libsdf.a test_complex_c89.c sdf_subs.h
	echo "compiling xtc89" $@
	$(COMPILER) $(FLAGS) test_complex_c89.c -I. -L. -lsdf -o xtc89

xtc99: libsdf.a test_complex_c99.c sdf_subs.h
	echo "compiling xtc99" $@
	$(COMPILER) $(FLAGS) test_complex_c99.c -I. -L. -lsdf -o xtc99

xtestf77: libsdf.a sdf_f77_tests.f
	echo "compiling xtestf77" $@
	$(F77) sdf_f77_tests.f $(F77FLAGS) -L. -lsdf -o xtestf77

xtestt: libsdf.a sdf_f77_testt.f
	echo "compiling xtestt" $@
	$(F77) sdf_f77_testt.f $(F77FLAGS) -L. -lsdf -o xtestt

xtestrbwb: libsdf.a test_rbwb.f
	echo "compiling xtestrbwb" $@
	$(F77) test_rbwb.f $(F77FLAGS) -L. -lsdf -o xtestrbwb

xtp: libsdf.a test_prim.f
	echo "compiling xtp" $@
	$(F77) test_prim.f $(F77FLAGS) -L. -lsdf -o xtp

xtest_dyn_f95: libsdf.a test_sdf_dynamic.F90
	echo "compiling xtest_dyn_f95" $@
	$(F95) test_sdf_dynamic.F90 $(F95FLAGS) -L. -lsdf -o xtest_dyn_f95

xtest_transp_f95: libsdf.a test_sdf_transpose.F90
	echo "compiling xtest_transp_f95" $@
	$(F95) test_sdf_transpose.F90 $(F95FLAGS) -L. -lsdf -o xtest_transp_f95

xtestf77_if64: libsdf_if64.a sdf_f77_tests.f
	echo "compiling xtestf77_if64" $@
	$(F77) sdf_f77_tests.f $(F77_IF64) $(F77FLAGS) -L. -lsdf_if64 -o xtestf77_if64

xtestt_if64: libsdf_if64.a sdf_f77_testt.f
	echo "compiling xtestt_if64" $@
	$(F77) sdf_f77_testt.f $(F77_IF64) $(F77FLAGS) -L. -lsdf_if64 -o xtestt_if64

xtestrbwb_if64: libsdf_if64.a test_rbwb.f
	echo "compiling xtestrbwb_if64" $@
	$(F77) test_rbwb.f $(F77_IF64) $(F77FLAGS) -L. -lsdf_if64 -o xtestrbwb_if64

xtp_if64: libsdf_if64.a test_prim.f
	echo "compiling xtp_if64" $@
	$(F77) test_prim.f $(F77_IF64) $(F77FLAGS) -L. -lsdf_if64 -o xtp_if64

xtest_dyn_f95_if64: libsdf_if64.a test_sdf_dynamic.F90
	echo "compiling xtest_dyn_f95_if64" $@
	$(F95) test_sdf_dynamic.F90 $(F95_If64) $(F95FLAGS) -L. -lsdf_if64 -o xtest_dyn_f95_if64

xtest_transp_f95_if64: libsdf_if64.a test_sdf_transpose.F90
	echo "compiling xtest_transp_f95_if64" $@
	$(F95) test_sdf_transpose.F90 $(F95_IF64) $(F95FLAGS) -L. -lsdf_if64 -o xtest_transp_f95_if64

all: sdf_browse libsdf_if64.a xtestmain xtestf77 xtest_dyn_f95 xtest_transp_f95  xtest_sdf_edit xlarge x3d x4d x5d xtc89 xtestt xtestrbwb xtp

all_if64: xtestf77_if64 xtest_dyn_f95_if64 xtest_transp_f95_if64 xtestt_if64 xtestrbwb_if64 xtp_if64
clean: 
	rm -f *.o  libsdf.a libsdf_if64.a sdf_browse sdf_browse.exe xtestmain xtestmain.exe sdf.def sdf_if64.def sdf.dll sdf_if64.dll sdf_imp_subs.lib xtestf77 xtestf77.exe xtestf77_if64 xtestf77_if64.exe xtest_transp_f95 xtest_transp_f95.exe xtest_transp_f95_if64 xtest_transp_f95_if64.exe xtest_dyn_f95 xtest_dyn_f95.exe xtest_dyn_f95_if64 xtest_dyn_f95_if64.exe xtest_sdf_edit xtest_sdf_edit.exe xlarge xlarge.exe x3d x3d.exe x4d x4d.exe x5d x5d.exe xtc89 xtc89.exe xtc99 xtc99.exe xtestt xtestt.exe xtestt_if64 xtestt_if64.exe xtestrbwb xtestrbwb.exe xtestrbwb_if64 xtestrbwb_if64.exe xtp xtp.exe xtp_if64 xtp_if64.exe xbe xbe.exe junkf77.sdf junk.sdf testoutlarge.sdf testout5d.sdf testout4d.sdf testout3d.sdf complex89.sdf complex99.sdf junkf77_2.sdf bzt.sdf testwbrb.dat testprim.dat gmon.out xforttest strlenfmt.h fflag-extra.inc
install:

	mkdir -p $(PREFIX)/lib
	mkdir -p $(PREFIX)/include
	mkdir -p $(PREFIX)/bin
	cp ./libsdf.a $(PREFIX)/lib
	cp ./libsdf_if64.a $(PREFIX)/lib
	cp ./sdf_subs.h $(PREFIX)/include
	cp ./sdf_browse $(PREFIX)/bin

uninstall:

	rm -f $(PREFIX)/lib/libsdf.a
	rm -f $(PREFIX)/lib/libsdf_if64.a
	rm -f $(PREFIX)/include/sdf_subs.h
	rm -f $(PREFIX)/bin/sdf_browse
