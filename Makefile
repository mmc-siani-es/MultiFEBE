#-------------------------------------------------------------------------------
# Main variables
# ==============

# Compilation mode: debug or optimal
MODE = debug

# GNU Fortran Compiler
FC = gfortran

# Documentation tool (f90doc, http://erikdemaine.org/software/f90doc/)
DOC = f90doc

# Executable
EXE = multifebe

# Used suffixes
.SUFFIXES:
.SUFFIXES: .f90 .o

# Extensions of temporal files
TMPFILES= *.o *__genmod.f90 *.mod
#-------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Paths
# =====

# Sources path
SRCPATH = src

# Libraries
FBEMPATH = lib/fbem
SDFPATH  = lib/sdf-0.75-RC4

# Path for f90doc documentation tool
DOCTOOLPATH = utils/f90doc

# Path for f90doc generated documentation
DOCPATH = doc

# Standard libraries and includes
STDLIBPATH = /usr/lib
STDINCPATH = /usr/include
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Flags for compilation
# =====================

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

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Libraries and headers
# =====================

# Clean FFLAGS
FFCOPT =

# Standard libraries and headers
FFCOPT += -L $(STDLIBPATH) -I $(STDINCPATH)

# FGSL (Fortran interface for GNU Scientific Library)
#FFCOPT += -I `pkg-config --cflags fgsl` -L `pkg-config --libs fgsl`

# Libraries and resource files
FFCOPT += -L $(FBEMPATH)/bin
FFCOPT += -I $(FBEMPATH)/include
FFCOPT += -I $(FBEMPATH)/src/resources_quasisingular_integration
FFCOPT += -I $(FBEMPATH)/src/resources_shape_functions
FFCOPT += -I $(FBEMPATH)/src/resources_quad_rules
#FFCOPT += -L $(SDFPATH)

# Check if is a 32-bit or 64-bit system
LBITS := $(shell getconf LONG_BIT)
ifeq ($(LBITS),32)
$(error "32-bit systems not supported")
endif

#
# ATLAS LIBRARY
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
SRCS += $(SRCPATH)/problem_variables.f90
SRCS += $(SRCPATH)/version.f90
SRCS += $(SRCPATH)/getopt.f90
SRCS += $(SRCPATH)/help.f90
SRCS += $(SRCPATH)/process_command_line_options.f90
SRCS += $(SRCPATH)/read_problem.f90
SRCS += $(SRCPATH)/read_sensitivity.f90
SRCS += $(SRCPATH)/read_settings.f90
SRCS += $(SRCPATH)/read_frequencies.f90
SRCS += $(SRCPATH)/read_parts.f90
SRCS += $(SRCPATH)/read_nodes.f90
SRCS += $(SRCPATH)/read_elements.f90
SRCS += $(SRCPATH)/read_materials.f90
SRCS += $(SRCPATH)/read_cross_sections.f90
SRCS += $(SRCPATH)/read_regions.f90
SRCS += $(SRCPATH)/read_boundaries.f90
SRCS += $(SRCPATH)/read_commands.f90
SRCS += $(SRCPATH)/read_be_bodyloads.f90
SRCS += $(SRCPATH)/read_fe_subregions.f90
SRCS += $(SRCPATH)/read_symmetry_planes.f90
SRCS += $(SRCPATH)/read_groups.f90
SRCS += $(SRCPATH)/read_discontinuous_boundary_elements.f90
SRCS += $(SRCPATH)/read_special_boundary_elements.f90
SRCS += $(SRCPATH)/read_internal.f90
SRCS += $(SRCPATH)/read_internal_points.f90
SRCS += $(SRCPATH)/read_internal_points_from_mesh.f90
SRCS += $(SRCPATH)/read_internal_elements.f90
SRCS += $(SRCPATH)/transformation_element_to_quarter_point.f90
SRCS += $(SRCPATH)/assign_default_conditions_bem_boundaries_laplace.f90
SRCS += $(SRCPATH)/assign_default_conditions_bem_boundaries_mechanics_static.f90
SRCS += $(SRCPATH)/assign_default_conditions_bem_boundaries_mechanics_harmonic.f90
SRCS += $(SRCPATH)/assign_default_conditions_fem_nodes.f90
SRCS += $(SRCPATH)/assign_default_conditions_fem_elements.f90
SRCS += $(SRCPATH)/read_conditions_bem_boundaries_laplace.f90
SRCS += $(SRCPATH)/read_conditions_bem_boundaries_mechanics_static.f90
SRCS += $(SRCPATH)/read_conditions_bem_boundaries_mechanics_harmonic.f90
SRCS += $(SRCPATH)/transfer_conditions_bem_boundaries_laplace.f90
SRCS += $(SRCPATH)/transfer_conditions_bem_boundaries_mechanics_static.f90
SRCS += $(SRCPATH)/transfer_conditions_bem_boundaries_mechanics_harmonic.f90
SRCS += $(SRCPATH)/read_conditions_bem_nodes_laplace.f90
SRCS += $(SRCPATH)/read_conditions_bem_nodes_mechanics_static.f90
SRCS += $(SRCPATH)/read_conditions_bem_nodes_mechanics_harmonic.f90
SRCS += $(SRCPATH)/read_conditions_fem_nodes_mechanics_static.f90
SRCS += $(SRCPATH)/read_conditions_fem_nodes_mechanics_harmonic.f90
SRCS += $(SRCPATH)/read_conditions_fem_elements_mechanics.f90
SRCS += $(SRCPATH)/read_incident_mechanics_harmonic.f90
SRCS += $(SRCPATH)/assign_default_bem_formulation.f90
SRCS += $(SRCPATH)/read_bem_formulation_boundaries.f90
SRCS += $(SRCPATH)/read_bem_formulation_selected_nodes.f90
SRCS += $(SRCPATH)/read_fem_node_options.f90
SRCS += $(SRCPATH)/read_element_options.f90
SRCS += $(SRCPATH)/read_export.f90
SRCS += $(SRCPATH)/read_export_sif.f90
SRCS += $(SRCPATH)/build_connectivity_be_fe.f90
SRCS += $(SRCPATH)/build_connectivity_be_bodyload_element_fe.f90
SRCS += $(SRCPATH)/build_connectivities_be_bodyload_elements_and_fe.f90
SRCS += $(SRCPATH)/build_connectivities_be_and_fe.f90
SRCS += $(SRCPATH)/build_connectivities_rigid_regions.f90
SRCS += $(SRCPATH)/read_input_file.f90
SRCS += $(SRCPATH)/build_data.f90
SRCS += $(SRCPATH)/build_data_at_collocation_points.f90
SRCS += $(SRCPATH)/build_data_at_geometrical_nodes.f90
SRCS += $(SRCPATH)/build_data_at_functional_nodes.f90
SRCS += $(SRCPATH)/build_data_of_be_elements.f90
SRCS += $(SRCPATH)/build_auxiliary_variables_laplace.f90
SRCS += $(SRCPATH)/build_auxiliary_variables_mechanics_static.f90
SRCS += $(SRCPATH)/build_auxiliary_variables_mechanics_harmonic.f90
SRCS += $(SRCPATH)/assemble_bem_harpot_equation.f90
SRCS += $(SRCPATH)/assemble_bem_harpot_equation_sa.f90
SRCS += $(SRCPATH)/assemble_bem_staela_equation.f90
SRCS += $(SRCPATH)/assemble_bem_staela_equation_sa.f90
SRCS += $(SRCPATH)/assemble_bem_bl_staela_equation.f90
SRCS += $(SRCPATH)/assemble_bem_bl_harela_equation.f90
SRCS += $(SRCPATH)/assemble_fem_staela_stiffness_matrix.f90
SRCS += $(SRCPATH)/assemble_fem_staela_coupling.f90
SRCS += $(SRCPATH)/assemble_bem_harela_equation.f90
SRCS += $(SRCPATH)/assemble_bem_harela_equation_sa.f90
SRCS += $(SRCPATH)/assemble_fem_harela_stiffness_matrix.f90
SRCS += $(SRCPATH)/assemble_fem_harela_dKdx.f90
SRCS += $(SRCPATH)/assemble_fem_harela_coupling.f90
SRCS += $(SRCPATH)/assemble_fem_harela_dfdx.f90
SRCS += $(SRCPATH)/assemble_bem_harpor_equation.f90
SRCS += $(SRCPATH)/assemble_bem_laplace_equation.f90
SRCS += $(SRCPATH)/assemble_bem_laplace_equation_sa.f90
SRCS += $(SRCPATH)/calculate_incident_mechanics_harmonic.f90
SRCS += $(SRCPATH)/build_lse_laplace.f90
SRCS += $(SRCPATH)/build_lse_laplace_sa.f90
SRCS += $(SRCPATH)/build_lse_mechanics_static.f90
SRCS += $(SRCPATH)/build_lse_mechanics_bem_staela.f90
SRCS += $(SRCPATH)/build_lse_mechanics_fem_staela.f90
SRCS += $(SRCPATH)/build_fem_K_static.f90
SRCS += $(SRCPATH)/build_fem_f_distributed_static.f90
SRCS += $(SRCPATH)/build_lse_mechanics_static_sa.f90
SRCS += $(SRCPATH)/build_lse_mechanics_bem_harpot.f90
SRCS += $(SRCPATH)/build_lse_mechanics_bem_harpot_sa.f90
SRCS += $(SRCPATH)/build_lse_mechanics_bem_harela.f90
SRCS += $(SRCPATH)/build_lse_mechanics_bem_harela_sa.f90
SRCS += $(SRCPATH)/build_lse_mechanics_bem_harpor.f90
SRCS += $(SRCPATH)/build_lse_mechanics_fem_harela.f90
SRCS += $(SRCPATH)/build_fem_K_harmonic.f90
SRCS += $(SRCPATH)/build_fem_f_distributed_harmonic.f90
SRCS += $(SRCPATH)/build_lse_mechanics_harmonic.f90
SRCS += $(SRCPATH)/build_lse_mechanics_harmonic_sa.f90
SRCS += $(SRCPATH)/solve_lse_r.f90
SRCS += $(SRCPATH)/solve_lse_c.f90
SRCS += $(SRCPATH)/assign_solution_laplace.f90
SRCS += $(SRCPATH)/assign_solution_laplace_sa.f90
SRCS += $(SRCPATH)/assign_solution_mechanics_static.f90
SRCS += $(SRCPATH)/assign_solution_mechanics_static_sa.f90
SRCS += $(SRCPATH)/assign_solution_mechanics_harmonic.f90
SRCS += $(SRCPATH)/assign_solution_mechanics_harmonic_sa.f90
SRCS += $(SRCPATH)/calculate_stresses_mechanics_static.f90
SRCS += $(SRCPATH)/calculate_stresses_mechanics_harmonic.f90
SRCS += $(SRCPATH)/calculate_internal_points_laplace.f90
SRCS += $(SRCPATH)/calculate_internal_points_laplace_sa.f90
SRCS += $(SRCPATH)/calculate_internal_points_mechanics_static.f90
SRCS += $(SRCPATH)/calculate_internal_points_mechanics_harmonic.f90
SRCS += $(SRCPATH)/calculate_internal_points_mechanics_bem_harpot.f90
SRCS += $(SRCPATH)/calculate_internal_points_mechanics_bem_harela.f90
SRCS += $(SRCPATH)/calculate_internal_points_mechanics_bem_harpor.f90
SRCS += $(SRCPATH)/export_data_at_geometrical_nodes.f90
SRCS += $(SRCPATH)/export_data_at_functional_nodes.f90
SRCS += $(SRCPATH)/export_solution_laplace.f90
SRCS += $(SRCPATH)/export_solution_laplace_sa.f90
SRCS += $(SRCPATH)/export_solution_mechanics_static.f90
SRCS += $(SRCPATH)/export_solution_mechanics_static_sa.f90
SRCS += $(SRCPATH)/export_solution_mechanics_static_nso.f90
SRCS += $(SRCPATH)/export_solution_mechanics_static_eso.f90
SRCS += $(SRCPATH)/export_solution_mechanics_static_tot.f90
SRCS += $(SRCPATH)/export_solution_mechanics_static_gmsh.f90
SRCS += $(SRCPATH)/export_solution_mechanics_harmonic.f90
SRCS += $(SRCPATH)/export_solution_mechanics_harmonic_nso.f90
SRCS += $(SRCPATH)/export_solution_mechanics_harmonic_eso.f90
SRCS += $(SRCPATH)/export_solution_mechanics_harmonic_tot.f90
SRCS += $(SRCPATH)/export_solution_mechanics_harmonic_ier.f90
SRCS += $(SRCPATH)/export_solution_mechanics_harmonic_sif.f90
SRCS += $(SRCPATH)/export_solution_mechanics_harmonic_gmsh.f90
SRCS += $(SRCPATH)/export_solution_mechanics_harmonic_sa.f90
SRCS += $(SRCPATH)/print_frequency.f90
SRCS += $(SRCPATH)/$(EXE).f90

# Object files
# ------------

OBJS=$(SRCS:.f90=.o)

# Export compiler and flags (for other Makefiles)
# -----------------------------------------------

export FC FFLAGS MODE

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Rules
# =====

# Default target by convention is 'all'
all: $(EXE)

# Documentation
documentation:
	@echo 'Generating documentation with $(DOC) ...'
	@perl $(DOCTOOLPATH)/$(DOC) $(SRCPATH)/*.f90
	@mkdir -p doc
	@mv *.html doc/.
	@echo 'Moving documentation to $(DOCPATH) ...'

# Build exe and copy to ./bin
$(EXE): cleanhelp sdf fbem $(OBJS)
	@echo 'Building $(EXE) ...'
	@$(FC) $(OBJS) $(FFLAGS) $(FFCOPT) -o $(EXE) $(wildcard $(FBEMPATH)/bin/*.o) $(wildcard $(SDFPATH)/sdf_subs.o)
	@echo 'Copying $(EXE) to bin/$(EXE) ...'
	@mkdir -p bin
	@cp $(EXE) bin/$(EXE)

# Build fbem library
fbem:
	@echo 'Building FBEM library'
	@$(MAKE) -C $(FBEMPATH)

sdf:
	@echo 'Building SDF library'
	@$(MAKE) -C $(SDFPATH)

# General rule for ".f90" files of the program
%.o: %.f90
	@echo 'Compiling $< ...'
	@$(FC) $(FFLAGS) $(FFCOPT) -c $< -o $@

cleanhelp:
	@echo 'Cleaning  help.o ...'
	@cd $(SRCPATH); rm -f help.o

# Executed only when "$ make clean"
.PHONY: clean
clean:
	@echo 'Cleaning ...'
	@rm -f $(TMPFILES)
	@cd $(SRCPATH); rm -f $(TMPFILES)
	@$(MAKE) -C $(FBEMPATH) clean
	@$(MAKE) -C $(SDFPATH) clean
# ------------------------------------------------------------------------------
