# ------------------------------------------------------------------------------
# Main variables
# ==============

#
# Force compiler and compilation mode and compiler
#
ifndef MODE
# Compilation mode: debug or optimal
MODE = debug
# GNU Fortran Compiler
FC = gfortran
endif

# Documentation tool (f90doc, http://erikdemaine.org/software/f90doc/)
DOC = f90doc

# Library name
LIBNAME = libfbem

# Used suffixes
.SUFFIXES:
.SUFFIXES: .f90 .o .a

# Extensions of temporal files
TMPFILES= *.o *__genmod.f90 *.mod
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Paths
# =====

# Sources path
SRCPATH = src

# Path for f90doc documentation tool
DOCTOOLPATH = ../../utils/f90doc

# Path for f90doc generated documentation
DOCPATH = doc

# Standard libraries and includes
STDLIBPATH = /usr/lib
STDINCPATH = /usr/include
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Flags for compilation
# =====================

ifndef MODE

# Clean FFLAGS
FFLAGS =

#
# GENERAL FLAGS
#
# Option to Specify the standard to which the program is expected to conform,
# which may be one of `f95', `f2003', `f2008', `gnu', or `legacy'.
FFLAGS += -std=gnu -ffree-line-length-none
# Option to show diagnostic messages to be issued by the compiler
FFLAGS += -Wall -Wextra
# Profile the code (then $ gprof ./program gmon.out)
#FFLAGS += -p
# Option to see file source after preprocessor
#FFLAGS += -E
# Option to run preprocessor
FFLAGS += -cpp
# Option to control integer default size (by default compiler uses int32)
#FFLAGS += -fdefault-integer-8
# Option to control real default size (by default compiler uses real32)
#FFLAGS += -fdefault-real-8
# Option to specify that no implicit typing is allowed, unless overridden by
# explicit IMPLICIT statements. This is the equivalent of adding implicit
# none to the start of every procedure.
FFLAGS += -fimplicit-none
# Option to make an unified executable
FFLAGS += -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
#FFLAGS += -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
# GNU OpenMP
#FFLAGS += -fopenmp -fno-automatic
FFLAGS += -fopenmp
# SDF Library
#FFLAGS += -lsdf

#
# DEBUG MODE
#
ifeq ($(MODE),debug)
# Option to make executable for debugging
FFLAGS += -g
# Option to check for certain conditions at run time
FFLAGS += -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -fno-inline
# Option to control code optimization (-O, -O1, -O2, -O3, -Ofast, -Og (debug))
FFLAGS += -Og
endif

#
# OPTIMAL MODE
#
ifeq ($(MODE),optimal)
# Tuning for specific machine: native, core2, sandybridge, ... (Link: https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html)
#FFLAGS += -march=native
FFLAGS += -march=core2
# Option to control code optimization (-O, -O1, -O2, -O3, -Ofast, -Og (debug))
FFLAGS += -O3
# Option to run the standard link-time optimizer
FFLAGS += -flto
endif

endif

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Libraries and headers
# =====================

# Clean FFCOPT
FFCOPT =

# Standard libraries and headers
FFCOPT += -L $(STDLIBPATH) -I $(STDINCPATH)

# Resource files
FFCOPT += -I $(SRCPATH)/resources_quad_rules
FFCOPT += -I $(SRCPATH)/resources_shape_functions


# Check if is a 32-bit or 64-bit system
LBITS := $(shell getconf LONG_BIT)
ifeq ($(LBITS),32)
$(error "32-bit systems not supported")
endif

#
# Con librerias ATLAS (LAPACK MULTICORE)
#
FFCOPT += -m64 -llapack_atlas -latlas -llapack -lblas

#
# OpenBLAS
#
#FFCOPT += -m64 -lopenblas
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Files
# =====

# Source files
# ------------

SRCS =
SRCS += $(SRCPATH)/numerical.f90
SRCS += $(SRCPATH)/string_handling.f90
SRCS += $(SRCPATH)/quad_rules.f90
SRCS += $(SRCPATH)/shape_functions.f90
SRCS += $(SRCPATH)/data_structures.f90
SRCS += $(SRCPATH)/mesh.f90
SRCS += $(SRCPATH)/telles_transformation.f90
SRCS += $(SRCPATH)/polar_transformation.f90
SRCS += $(SRCPATH)/geometry.f90
SRCS += $(SRCPATH)/symmetry.f90
SRCS += $(SRCPATH)/quasisingular_integration.f90
SRCS += $(SRCPATH)/bem_general.f90
SRCS += $(SRCPATH)/bem_stapot2d.f90
SRCS += $(SRCPATH)/bem_stapot3d.f90
SRCS += $(SRCPATH)/bem_harpot2d.f90
SRCS += $(SRCPATH)/bem_harpot3d.f90
SRCS += $(SRCPATH)/bem_staela2d.f90
SRCS += $(SRCPATH)/bem_staela3d.f90
SRCS += $(SRCPATH)/bem_harela2d.f90
SRCS += $(SRCPATH)/bem_harela3d.f90
SRCS += $(SRCPATH)/bem_harpor2d.f90
SRCS += $(SRCPATH)/bem_harpor3d.f90
SRCS += $(SRCPATH)/fem_solids.f90
SRCS += $(SRCPATH)/fem_beams.f90
SRCS += $(SRCPATH)/fem_shells.f90
SRCS += $(SRCPATH)/harpot_incident_field.f90
SRCS += $(SRCPATH)/harela_incident_field.f90
SRCS += $(SRCPATH)/harpor_incident_field.f90
SRCS += $(SRCPATH)/gmsh.f90

# Object files
# ------------

OBJS=$(SRCS:.f90=.o)

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Rules
# =====

# Default target by convention is 'all'
all: $(LIBNAME)

# Documentation
documentation:
	@echo 'Generating documentation with $(DOC) ...'
	@perl $(DOCTOOLPATH)/$(DOC) $(SRCPATH)/*.f90
	@mkdir -p doc
	@mv *.html doc/.
	@echo 'Moving documentation to $(DOCPATH) ...'

# Build fbem library and copy *.o to ./bin and *.mod to ./include
$(LIBNAME): $(OBJS)
	@echo 'Copying object files to bin/ ...'
	@mkdir -p bin
	@cp $(SRCPATH)/*.o bin/.
	@echo 'Copying modules to include/ ...'
	@mkdir -p include
	@cp *.mod include/.

# General rule for ".f90" files of the program
%.o: %.f90
	@echo 'Compiling $< ...'
	@$(FC) $(FFLAGS) $(FFCOPT) -c $< -o $@

# Executed only when "$ make clean"
.PHONY: clean
clean:
	@echo 'Cleaning ...'
	@rm -f $(TMPFILES)
	@cd $(SRCPATH); rm -f $(TMPFILES)
	@cd bin; rm -f $(TMPFILES)
	@cd include; rm -f $(TMPFILES)
# ------------------------------------------------------------------------------
