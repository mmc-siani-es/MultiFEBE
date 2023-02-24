! ---------------------------------------------------------------------
! Copyright (C) 2014-2022 Universidad de Las Palmas de Gran Canaria:
!                         Jacob D.R. Bordon
!                         Guillermo M. Alamo
!                         Juan J. Aznarez
!                         Orlando Maeso.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
! ---------------------------------------------------------------------


!! @author Jacob David Rodriguez Bordon (jacobdavid.rodriguezbordon@ulpgc.es)
!!
!! @version 2.0
!!
!! <b> Problem variables module </b>
module problem_variables

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_mesh_module
  use fbem_quasisingular_integration

  ! No implicit variables are allowed
  implicit none

  ! The variables are arranged by topic.

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! Files configuration
  ! -------------------
  !
  character(len=fbem_path_max_length)     :: input_filename        !! Input file path
  character(len=fbem_path_max_length)     :: output_filename       !! Output files path (without extension)
  character(len=fbem_path_max_length)     :: input_filedir         !! Input file directory path
  character(len=8)                        :: timestamp_date_start  !! Timestamps
  character(len=10)                       :: timestamp_time_start  !! Timestamps
  character(len=8)                        :: timestamp_date_tick   !! Timestamps
  character(len=10)                       :: timestamp_time_tick   !! Timestamps
  character(len=8)                        :: timestamp_date_end    !! Timestamps
  character(len=10)                       :: timestamp_time_end    !! Timestamps
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! Problem definition
  ! ------------------
  !
  type(fbem_problem)                      :: problem
  type(fbem_mesh)                         :: design_mesh
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! Settings
  ! --------
  !
  ! Quasi-singular integration
  real(kind=real64)                       :: qsi_relative_error    !! Relative error of quasi-singular integration.
  type(fbem_qs_parameters)                :: qsi_parameters        !! Parameters for quasi-singular integration.
  integer                                 :: qsi_ns_max            !! Maximum number of subdivisions in quasi-singular integration (>=0).
  integer                                 :: n_precalsets          !! Number of precalculated datasets at integration points.
  integer, allocatable                    :: precalset_gln(:)      !! Number of Gauss-Legendre quadrature points of each set (<=30). The order is O=2N-1.
  integer                                 :: precalset_gln_min     !! Minimum of precalset_gln()
  integer                                 :: precalset_gln_max     !! Maximum of precalset_gln()
  ! Geometrical
  real(kind=real64)                       :: geometric_tolerance   !! General geometric tolerance.
  logical                                 :: collapse_nodal_pos    !! Collapse nodal positions of nodes.
  ! Linear system of equations
  logical                                 :: lse_scaling           !! Perform scaling of the LSE.
  logical                                 :: lse_condition         !! Estimate the condition number of the LSE.
  logical                                 :: lse_refine            !! Refine the solution of the LSE.
  ! Mesh file
  integer                                 :: mesh_file_mode        !! 0 (mesh is in the input file), 1 (mesh is in an auxiliary file), 2 (mesh is in an auxiliary file in gmsh format)
  character(len=fbem_path_max_length)     :: mesh_filename         !! Filename of the mesh
  ! Frequencies
  character(len=fbem_path_max_length)     :: frequencies_filename  !! Filename of frequencies
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! Linear system of equations (LSE)
  ! --------------------------------
  !
  integer                                 :: n_dof                 !! Number of active degrees of freedom.
  real(kind=real64), allocatable          :: A_r(:,:)              !! Matrix A (real) of the A·x=b linear system of equations.
  real(kind=real64), allocatable          :: b_r(:,:)              !! Vector b (real) of the A·x=b linear system of equations.
  real(kind=real64), allocatable          :: bsa_r(:,:)            !! Vector b (real) of the A·x=b linear system of equations for sensitivity analysis.
  complex(kind=real64), allocatable       :: A_c(:,:)              !! Matrix A (complex) of the A·x=b linear system of equations.
  complex(kind=real64), allocatable       :: b_c(:,:)              !! Vector b (complex) of the A·x=b linear system of equations.
  complex(kind=real64), allocatable       :: bsa_c(:,:)            !! Vector b (complex) of the A·x=b linear system of equations for sensitivity analysis.
  integer, allocatable                    :: fact_ipiv(:)          !! A reordering after factorization
  character(len=1)                        :: scal_equed            !! 'N':  No equilibration; 'R':  Row equilibration; 'C':  Column equilibration; 'B':  Both row and column equilibration.
  real(kind=real64), allocatable          :: scal_r(:)             !! The row scale factors for A.
  real(kind=real64), allocatable          :: scal_c(:)             !! The column scale factors for A.
  integer                                 :: Aodim                 !! Dimensions of Ao
  real(kind=real64), allocatable          :: Ao_r(:,:)             !! Matrix A (real) of the A·x=b linear system of equations (original - non-factorized A, required if condition or refine)
  complex(kind=real64), allocatable       :: Ao_c(:,:)             !! Matrix A (complex) of the A·x=b linear system of equations (original - non-factorized A, required if condition or refine)
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! Frequencies
  ! -----------
  !
  character                               :: frequency_units       !! Frequency units: 'f' for Hz, 'w' for rad/s.
  integer                                 :: n_frequencies         !! Number of frequencies.
  real(kind=real64), allocatable          :: frequency(:)          !! Frequencies, they must be in rad/s before any calculation.
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! Incident harmonic wave fields
  ! -----------------------------
  !
  integer                                 :: n_incidentfields      !! Number of incident wave fields
  type(fbem_incidentfield), allocatable   :: incidentfield(:)      !! Vector of incidentfield node structures
  integer                                 :: incidentfield_eid_min !! Minimum value of external identifiers.
  integer                                 :: incidentfield_eid_max !! Maximum  value of external identifiers.
  integer, allocatable                    :: incidentfield_iid(:)  !! Vector to convert external identifier to internal identifier.
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! Symmetry planes configuration
  ! -----------------------------
  !
  integer                                 :: n_symplanes           !! Number of symmetry planes: 0 ,1, 2 or 3.
  integer                                 :: n_symelements         !! Number of symmetrical elements
  integer                                 :: symplane_eid(3)       !! External identifiers of each symmetry plane: 1, 2, 3. The eid is 1 for symmetry of x coordinate, 2 for y, 3 for z.
  integer                                 :: symplane_iid(3)       !! Internal identifiers of each symmetry plane: 1, 2, 3.
  real(kind=real64)                       :: symplane_m(3,3)       !! Multiplier for each sym. plane for geometrical reflection: 1. or -1.
  real(kind=real64)                       :: symplane_s(3)         !! Multiplier for each sym. plane for scalar variables reflection: 1. or -1.
  real(kind=real64)                       :: symplane_t(3,3)       !! Multiplier for each sym. plane for translational vectorial variables reflection: 1. or -1.
  real(kind=real64)                       :: symplane_r(3,3)       !! Multiplier for each sym. plane for rotational vectorial variables reflection: 1. or -1.
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! Nodes
  ! -----
  !
  integer                                 :: n_nodes               !! Number of nodes.
  type(fbem_node), allocatable            :: node(:)               !! Vector of nodes data structures.
  integer                                 :: node_eid_min          !! Minimum value of external identifiers.
  integer                                 :: node_eid_max          !! Maximum  value of external identifiers.
  integer, allocatable                    :: node_iid(:)           !! Vector to convert external identifier to internal identifier.
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! Elements
  ! --------
  !
  integer                                 :: n_elements            !! Number of elements.
  type(fbem_element), allocatable         :: element(:)            !! Vector of elements data structures.
  integer                                 :: element_eid_min       !! Minimum value of external identifiers.
  integer                                 :: element_eid_max       !! Maximum value of external identifiers.
  integer, allocatable                    :: element_iid(:)        !! Vector to convert external identifier to internal identifier.
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! Sub-elements
  ! ------------
  !
  integer                                 :: n_subedges            !! Number of edge sub-elements.
  type(fbem_element), allocatable         :: subedge(:)            !! Vector of edge subelements data structures.
  integer                                 :: n_subfaces            !! Number of edge sub-elements.
  type(fbem_element), allocatable         :: subface(:)            !! Vector of edge subelements data structures.
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! Parts (mesh partitions)
  ! -----------------------
  !
  integer                                 :: n_parts               !! Number of parts.
  type(fbem_part), allocatable            :: part(:)               !! Vector of parts data structures.
  integer                                 :: part_eid_min          !! Minimum value of external identifiers.
  integer                                 :: part_eid_max          !! Maximum value of external identifiers.
  integer, allocatable                    :: part_iid(:)           !! Vector to convert external identifier to internal identifier.
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! BE body loads
  ! -------------
  !
  integer                                 :: n_be_bodyloads        !! Number of boundaries.
  type(fbem_be_bodyload), allocatable     :: be_bodyload(:)        !! Vector of elements data structures.
  integer                                 :: be_bodyload_eid_min   !! Minimum value of external identifiers.
  integer                                 :: be_bodyload_eid_max   !! Maximum value of external identifiers.
  integer, allocatable                    :: be_bodyload_iid(:)    !! Vector to convert external identifier to internal identifier.
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! BE boundaries
  ! -------------

  ! cambiar "boundary" por "BE_boundary" ?¿?¿?¿?¿?

  !
  integer                                 :: n_boundaries          !! Number of boundaries.
  type(fbem_boundary), allocatable        :: boundary(:)           !! Vector of elements data structures.
  integer                                 :: boundary_eid_min      !! Minimum value of external identifiers.
  integer                                 :: boundary_eid_max      !! Maximum value of external identifiers.
  integer, allocatable                    :: boundary_iid(:)       !! Vector to convert external identifier to internal identifier.
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! FE subregions
  ! -------------
  !
  integer                                 :: n_fe_subregions       !! Number of FE subregions.
  type(fbem_region), allocatable          :: fe_subregion(:)       !! Vector of elements data structures.
  integer                                 :: fe_subregion_eid_min  !! Minimum value of external identifiers.
  integer                                 :: fe_subregion_eid_max  !! Maximum value of external identifiers.
  integer, allocatable                    :: fe_subregion_iid(:)   !! Vector to convert external identifier to internal identifier.
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! Materials
  ! ---------
  !
  integer                                 :: n_materials           !! Number of regions.
  type(fbem_material), allocatable        :: material(:)           !! Vector of elements data structures.
  integer                                 :: material_eid_min      !! Minimum value of external identifiers.
  integer                                 :: material_eid_max      !! Maximum value of external identifiers.
  integer, allocatable                    :: material_iid(:)       !! Vector to convert external identifier to internal identifier.
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! Regions
  ! -------
  !
  integer                                 :: n_regions             !! Number of regions.
  integer                                 :: n_be_regions          !! Number of regions of class be.
  integer                                 :: n_fe_regions          !! Number of regions of class fe.
  type(fbem_region), allocatable          :: region(:)             !! Vector of elements data structures.
  integer                                 :: region_eid_min        !! Minimum value of external identifiers.
  integer                                 :: region_eid_max        !! Maximum value of external identifiers.
  integer, allocatable                    :: region_iid(:)         !! Vector to convert external identifier to internal identifier.
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! Internal points (BEM post-processing)
  ! -------------------------------------
  !
  integer                                 :: n_internalpoints      !! Number of internal points.
  type(fbem_internalpoint), allocatable   :: internalpoint(:)      !! Vector of internal points data structures.
  integer                                 :: internalpoint_eid_min !! Minimum value of external identifiers.
  integer                                 :: internalpoint_eid_max !! Maximum value of external identifiers.
  integer, allocatable                    :: internalpoint_iid(:)  !! Vector to convert external identifier to internal identifier.
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! Internal elements (BEM post-processing)
  ! ---------------------------------------
  !
  logical                                 :: internalelements            !! True if an internal elements mesh has been defined
  character(len=fbem_path_max_length)     :: internalelements_filename   !! Filename of the internal elements' mesh
  integer                                 :: internalelements_fileformat !! File format
  type(fbem_mesh)                         :: internalelements_mesh       !! Internal elements mesh data structure
  integer                                 :: internalelements_order      !! Order of the mesh for the stress resultants calculation ! to be deprecated
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! Groups
  ! ------
  !
  integer                                 :: n_groups              !! Number of groups.
  type(fbem_group), allocatable           :: group(:)              !! Vector of groups data structures.
  integer                                 :: group_eid_min         !! Minimum value of external identifiers.
  integer                                 :: group_eid_max         !! Maximum value of external identifiers.
  integer, allocatable                    :: group_iid(:)          !! Vector to convert external identifier to internal identifier.
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! Export settings
  ! ---------------
  !
  real(kind=real64)                       :: vectors_scale_factor  !! Scale factor for vectors
  character(len=32)                       :: fmt_integer           !! Format for integer export (Fortran specifier format)
  character(len=32)                       :: fmt_real              !! Format for real export (Fortran specifier format)
  integer                                 :: complex_notation      !! Notation for complex: 1: abs & arg, 2: real & imag.
  logical                                 :: export_overwrite      !! ¿Overwrite existent files (native file formats)?
  logical                                 :: export_geo_mesh_data  !! ¿Export geometrical e*.dat file?
  logical                                 :: export_fun_mesh_data  !! ¿Export functional e*.dat file?
  logical                                 :: export_nso            !! ¿Export .nso file?
  logical                                 :: export_eso            !! ¿Export .eso file?
  logical                                 :: export_pos            !! ¿Export .pos file (Gmsh)?
  character(len=32)                       :: fmt_real_pos          !! Format for real export (Fortran specifier format) - .pos file
  logical                                 :: export_sif            !! ¿Export .sif file?
  logical                                 :: export_wsp            !! ¿Export .wsp file?
  logical                                 :: export_tot            !! ¿Export .tot file?
  integer                                 :: tot_xm                !! Resultant moment in .tot file: 1 (origin), 2 (centroid of the boundary)
  logical                                 :: tot_apply_symmetry    !! Resultants (tot y ier) applying symmetry
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! Memory settings
  ! ---------------
  !
  integer(kind=int64)                     :: max_memory            !! Maximum memory in matrix A and vector b
  ! --------------------------------------------------------------------------------------------------------------------------------

end module problem_variables
