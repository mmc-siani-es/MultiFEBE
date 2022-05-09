Sept. 1, 2007 - GHF

Motivation:  ./docs/SDF_MANIFESTO.txt

General overview: ./docs/SDF-intro-spd-2006.pdf

To install:  ./docs/INSTALL.txt

Usage in Fortran, C, and IDL programs:  ./docs/SDF_USAGE_NOTES.txt

Recent Changes:  ./docs/CHANGELOG.txt

Examples:  See source code in sdf_browse.c, main_test_sdf.c, sdf_test_large.c,
test_sdf_3d.c, test_sdf_4d.c, test_sdf_5d.c test_sdf_complex_c89.c,
test_sdf_complex_c99.c, and sdf_test_edit.c (C examples) 
and sdf_f77_tests.f, sdf_f77_testt.f (Fortran 77 examples) and
test_sdf_dynamic.F90, test_sdf_transpose.F90
(using allocatable arrays, requires a Fortran 90/95 compiler).
For examples of primitive binary writes (outside the context of the SDF
file format) from Fortran, see the example files test_rbwb.f (using sdf_wb_f77
and sdf_rb_f77), and test_prim.f (using Fortran callable C disk I/O
subroutines that are now in the SDF library).

Try using the sdf file browser sdf_browse to examine the contents of the
sample file sdf_example_file.sdf, and/or using the IDL procedures to
read/write or edit this (or other files).

To easily create an sdf file with multiple datasets in IDL, create a number
of variables in an IDL session, and then enter the command
sdf_write_all,<fname> where <fname> is your test output file. (sdf_write_all
will require an IDL version of 6.1 or later).

Questions, contact George H. Fisher (fisher at ssl dot berkeley dot edu)
