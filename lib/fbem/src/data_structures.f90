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
!! <b>This module implements data structures for problems implementation.</b>
!!
!! Note that this module is independent from calculation modules.
module fbem_data_structures

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_string_handling
  use fbem_shape_functions

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all are private
  private

  ! ================================================================================================================================
  ! Public data structures and parameters
  !
  ! PROBLEM
  public :: fbem_problem
  public :: fbem_mechanics
  public :: fbem_mechanics_plane_strain
  public :: fbem_mechanics_plane_stress
  public :: fbem_laplace
  public :: fbem_static
  public :: fbem_harmonic
  public :: fbem_modal
  public :: fbem_transient
  ! GENERAL PARAMETERS
  public :: fbem_be
  public :: fbem_fe
  public :: fbem_full_space
  public :: fbem_half_space
  public :: fbem_multilayered_half_space
  public :: fbem_point
  public :: fbem_line
  public :: fbem_plane
  ! NODE
  public :: fbem_node
  public :: fbem_sbie
  public :: fbem_sbie_mca
  public :: fbem_hbie
  public :: fbem_dual_burton_miller
  public :: fbem_dual_boundary
  ! ELEMENT
  public :: fbem_element
  ! PART
  public :: fbem_part
  public :: fbem_part_be_boundary
  public :: fbem_part_fe_subregion
  public :: fbem_part_be_bodyload
  ! BE BODY LOAD
  public :: fbem_be_bodyload
  public :: fbem_bl_coupling_beam_tip
  public :: fbem_bl_coupling_beam_line
  public :: fbem_bl_coupling_shell_edge
  public :: fbem_bl_coupling_shell_surface
  ! BOUNDARY
  public :: fbem_boundary
  public :: fbem_boundary_class_ordinary
  public :: fbem_boundary_class_cracklike
  public :: fbem_boundary_coupling_be
  public :: fbem_boundary_coupling_be_be
  public :: fbem_boundary_coupling_be_fe
  public :: fbem_boundary_coupling_be_fe_be
  ! MATERIAL
  public :: fbem_material
  ! REGION
  public :: fbem_region
  public :: fbem_potential
  public :: fbem_viscoelastic
  public :: fbem_poroelastic
  public :: fbem_rigid
  public :: fbem_elastic
  public :: fbem_bardet
  public :: fbem_brt_cp1
  public :: fbem_brt_cp2
  public :: fbem_brt_cpm
  public :: fbem_elastic_orthotropic
  ! INTERNAL POINT
  public :: fbem_internalpoint
  ! INCIDENT WAVE FIELD
  public :: fbem_incidentfield
  ! GROUP OF MESH/SUBMESH OBJECTS
  public :: fbem_group
  public :: fbem_group_type_nodes
  public :: fbem_group_type_elements
  public :: fbem_group_type_subedges
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! Public subroutines
  !
  ! CONNECTIVITY
  public :: fbem_node_nodes_connectivity
  public :: fbem_node_elements_connectivity
  public :: fbem_node_parts_connectivity
  public :: fbem_part_nodes_connectivity
  public :: fbem_node_symplanes_connectivity
  public :: fbem_element_symplanes_connectivity
  public :: fbem_check_nodes_symplanes_configuration
  public :: fbem_build_mesh_subelements
  public :: fbem_check_and_characterize_be_mesh
  public :: fbem_check_and_characterize_fe_mesh
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! PROBLEM
  !
  !! Problem data structure
  type fbem_problem
    character(len=fbem_string_max_length) :: description                !! Description of the problem
    integer                               :: type                       !! Type of problem
    integer                               :: subtype                    !! Subtype of problem
    integer                               :: analysis                   !! Type of analysis
    integer                               :: n                          !! Number of dimensions of the ambient space
    logical                               :: sensitivity                !! True if sensitivity analysis has to be performed
    integer                               :: n_designvariables          !! Number of design variables
  end type fbem_problem
  !
  ! Associated parameters
  !
  ! Internal identifiers for types of problem
  integer, parameter                  :: fbem_laplace  =1               !! Laplace problem
  integer, parameter                  :: fbem_mechanics=2               !! Mechanics problem
  ! Internal identifiers for sub-types of mechanics problems
  integer, parameter                  :: fbem_mechanics_plane_strain=1  !! 2D plane strain problem
  integer, parameter                  :: fbem_mechanics_plane_stress=2  !! 2D plane stress problem
  ! Internal identifiers for types of analysis
  integer, parameter                  :: fbem_static   =1               !! Static analysis
  integer, parameter                  :: fbem_harmonic =2               !! Time-harmonic analysis
  integer, parameter                  :: fbem_modal    =3               !! Modal analysis
  integer, parameter                  :: fbem_transient=4               !! Transient analysis
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! GENERAL PARAMETERS
  !
  ! Internal identifier for the class of entity
  integer, parameter                  :: fbem_be=1                      !! If a boundary element entity (node/element/region)
  integer, parameter                  :: fbem_fe=2                      !! If a finite element entity (node/element/region)
  ! Internal identifier for the type of space
  integer, parameter                  :: fbem_full_space=1              !! Full space
  integer, parameter                  :: fbem_half_space=2              !! Half-space or half-plane
  integer, parameter                  :: fbem_multilayered_half_space=3 !! Multilayered half-space or half-plane
  ! Internal identifier for the type of entity
  integer, parameter                  :: fbem_point=1                   !! Point
  integer, parameter                  :: fbem_line =2                   !! Line
  integer, parameter                  :: fbem_plane=3                   !! Plane
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! NODE
  !
  !! Node data structure
  type fbem_node
    !
    ! Basic data
    !
    integer                           :: id                             !! External identifier
    real(kind=real64), allocatable    :: x(:)                           !! Position vector (geometrical node)
    real(kind=real64), allocatable    :: dxda(:,:)                      !! Design variable velocity vector dx_k/da_j for sensitivity analysis.
    ! DOF management
    integer                           :: n_dof                          !! Number of DOF of the node (for FE nodes)
    real(kind=real64), allocatable    :: local_axes(:,:)                !! Local axes at the node (for FE nodes)
    integer, allocatable              :: row(:,:)                       !! Mapping between the DOF and the LSE, row of each equation: row(eq_index,group)
    integer, allocatable              :: col(:,:)                       !! Mapping between the DOF and the LSE, column of each degree of freedom: col(dof_index,group)
    !
    ! Master-slave relationship for rigid links
    !
    integer                           :: rigid_link                     !! 0 (normal node), 1 (master node of a rigid link), 2 (slave node of a rigid link)
    integer                           :: master                         !! Internal identifier of the master node (if it is a slave node)
    integer                           :: n_slaves                       !! Number of slave nodes (if it is a master node)
    integer, allocatable              :: slave(:)                       !! List of slave nodes
    ! Values of the node
    real   (kind=real64), allocatable :: value_r(:,:)                   !! Values of variables (real): value_r(dof_index,face)
    complex(kind=real64), allocatable :: value_c(:,:)                   !! Values of variables (complex): value_c(dof_index,face)
    complex(kind=real64), allocatable :: incident_c(:,:)                !! Values of incident field (complex): incident_c(dof_index,face)
    ! Sensitivity of the variables of the node
    real(kind=real64), allocatable    :: dvda_r(:,:,:)                  !! Derivative of variables with respect to design variables (real): dvda_r(dof_index,group,design_variable)
    complex(kind=real64), allocatable :: dvda_c(:,:,:)                  !! Derivative of variables with respect to design variables (complex): dvda_c(dof_index,group,design_variable)
    !
    ! Connectivity to the mesh
    !
    ! With other nodes sharing the same position (distance<=tolerance)
    integer                           :: n_nodes                        !! Number of nodes
    integer, allocatable              :: node(:)                        !! Nodes (internal identifiers)
    ! With a coupled node
    integer                           :: coupled_node                   !! Coupled node with a different type of part (internal identifier)
    ! With elements that contain the node
    integer                           :: n_elements                     !! Number of elements
    integer, allocatable              :: element(:)                     !! Elements (internal identifiers)
    integer, allocatable              :: element_node_iid(:)            !! Internal identifier of the node in each element
    integer, allocatable              :: element_node_loctype(:)        !! Type of node location in each element: 1: vertex, 2: edge interior, 3: face interior, 4: volume interior.
    !! Dimensional degree of the connection of the node to its elements:
    !! 0: Unconnected node.
    !! 1: Only connected to 1D elements.
    !! 2: Only connected to 2D elements.
    !! 3: Only connected to 3D elements.
    !! 4: Connected to 1D and 2D elements, i.e. the node is connected with >=2 elements which are 1D and 2D elements.
    !! 5: Connected to 1D and 3D elements, i.e. the node is connected with >=2 elements which are 1D and 3D elements.
    !! 6: Connected to 2D and 3D elements, i.e. the node is connected with >=2 elements which are 2D and 3D elements.
    !! 7: Connected to 1D, 2D and 3D elements, i.e. the node is connected with >=3 elements which are 1D, 2D and 3D elements.
    integer                           :: dimensional_degree
    ! With parts
    integer                           :: n_parts                        !! Number of parts
    integer, allocatable              :: part(:)                        !! Parts (internal identifier)
    ! With symmetry planes that pass through the node
    integer                           :: n_symplanes                    !! Number of symmetry planes
    integer, allocatable              :: symplane(:)                    !! Symmetry planes (internal identifiers)
    ! Mesh markers
    logical                           :: in_boundary                    !! True if the node is in the boundary of its mesh.
    logical                           :: is_singular                    !! True if the node is in a singular point of its mesh (point belonging to >2 elements for 1D elements, edge belonging to >2 2D elements).
    !
    ! Connectivity to the design mesh
    !
    integer                           :: dm_n_elements                  !! Number of elements
    integer, allocatable              :: dm_element(:)                  !! Elements (internal identifiers)
    integer, allocatable              :: dm_element_d(:)                !! Dimension of the elements (internal identifiers)
    real(kind=real64), allocatable    :: dm_element_xi(:,:)             !! Local coordinates at each element
    !
    ! Conditions at the node
    !
    integer, allocatable              :: ctype(:,:)                     !! Type of condition: ctype(dof_index,group)
    integer, allocatable              :: cvalue_i(:,:,:)                !! Condition values (integer): cvalue_i(dof_index,par_index,group)
    real(kind=real64), allocatable    :: cvalue_r(:,:,:)                !! Condition values (real   ): cvalue_r(dof_index,par_index,group)
    complex(kind=real64), allocatable :: cvalue_c(:,:,:)                !! Condition values (complex): cvalue_c(dof_index,par_index,group)
    !
    ! Vectors at functional node
    !
    real(kind=real64), allocatable    :: x_fn(:)                        !! Position vector (functional node)
    real(kind=real64), allocatable    :: n_fn(:)                        !! Unit normal
    real(kind=real64), allocatable    :: t1_fn(:)                       !! Unit tangent 1
    real(kind=real64), allocatable    :: t2_fn(:)                       !! Unit tangent 2
    real(kind=real64), allocatable    :: t3_fn(:)                       !! Unit tangent 3
    !
    ! BEM node formulation options
    !
    ! Possible BEM formulations of a node by combining SBIE, SBIE with MCA, and HBIE (which has always MCA):
    !  - Ordinary boundary:
    !    -SBIE
    !    -SBIE-MCA
    !    -HBIE
    !    -SBIE     + beta*HBIE (Dual Burton & Miller)
    !    -SBIE-MCA + beta*HBIE (Dual Burton & Miller)
    !  - Crack-like boundary (Dual Boundary Element Method):
    !    -SBIE / HBIE
    !    -SBIE-MCA / HBIE
    !
    ! Multiple Collocation Approach (MCA) consists on collocating the BIEs at collocation points inside the elements, and sum
    ! all the equations related with the associated node. It is applied only for nodes at edges/vertices. See Gallego and Dominguez,
    ! "Hypersingular BEM for transient elastodynamics". International Journal for Numerical Methods in Engineering 39, 1996.
    !
    integer                           :: sbie                           !! SBIE usage: 0: none, 1: SBIE, 2: SBIE-MCA
    integer                           :: hbie                           !! HBIE usage: 0: none, 1: HBIE
    integer                           :: dual                           !! Dual usage of BIEs: 0: none, 1: Dual B&M, 2: DBEM
    logical                           :: dual_is_common                 !! True if for a dual formulation, both the SBIE and the HBIE are collocated at the same points.
    real(kind=real64)                 :: alpha                          !! Parameter alpha for Dual B&M: if static analysis (beta=alpha), if harmonic analysis (beta=alpha*i/wavenumber).
    ! Special: for BE lineloads end nodes touching a boundary (sbie with 0.5 freeterm)
    logical                           :: sbie_lineload_end_boundary
    !
    ! Export settings
    !
    logical                           :: export                         !! Export flag
    logical                           :: export_sif                     !! Calculate and export SIF (Stress Intensity Factor)
  end type fbem_node
  !
  ! Associated parameters for the BEM
  !
  ! Internal identifier for the type of SBIE formulation
  integer, parameter                  :: fbem_sbie                  =1  !! SBIE collocated at nodes
  integer, parameter                  :: fbem_sbie_mca              =2  !! SBIE collocated away from edges/vertices (MCA)
  ! Internal identifier for the type of HBIE formulation
  integer, parameter                  :: fbem_hbie                  =1  !! HBIE collocated away from edges/vertices (MCA)
  ! Internal identifier for the type of DUAL formulation
  integer, parameter                  :: fbem_dual_burton_miller    =1  !! SBIE+alpha*HBIE
  integer, parameter                  :: fbem_dual_boundary         =2  !! SBIE and HBIE for crack-like boundaries
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! ELEMENT
  !
  !! Element data structure
  type fbem_element
    !
    ! Basic data
    !
    integer                           :: id                             !! External identifier
    integer                           :: type                           !! Base type of element
    integer                           :: n_dimension                    !! Number of dimensions of the element: 0 (point), 1 (line), 2 (surface) or 3 (volume).
    real(kind=real64)                 :: size                           !! Size of the element: none, length, area, or volume.
    !
    ! Connectivity to the mesh
    !
    ! With a part
    integer                           :: part                           !! Part (internal identifier)
    ! With nodes of the element (functional)
    integer                           :: n_nodes                        !! Number of nodes
    integer, allocatable              :: node(:)                        !! Nodes (internal identifiers)
    ! Connectivity to its sub-elements (if an element)
    integer                           :: n_subedges                     !! Number of edges
    integer, allocatable              :: subedge(:)                     !! Edges (internal identifiers)
    logical, allocatable              :: subedge_reversion(:)           !! Reversion of the edge for the element definition
    integer                           :: n_subfaces                     !! Number of faces
    integer, allocatable              :: subface(:)                     !! Faces (internal identifiers)
    logical, allocatable              :: subface_reversion(:)           !! Reversion of the face for the element definition
    ! Connectivity to its super-elements (if a sub-element)
    integer                           :: n_supelements                  !! Number of super-elements
    integer, allocatable              :: supelement(:)                  !! Super-element (internal identifiers)
    integer, allocatable              :: supelement_element_iid(:)      !! Internal identifier of the element in each super-element
    ! Connectivity to a FE with the same dimension or higher (if a BE)
    ! Connectivity to a BE with the same dimension (if a FE)
    integer                           :: element                        !! Element (internal identifier), if 0, it is uncoupled.
    integer, allocatable              :: element_node(:)                !! The coupled node of each corresponding node (element_node(k)=i means that the node k of the current element is connected to the node i)
    ! With symmetry planes that pass through the whole element
    integer                           :: n_symplanes                    !! Number of symmetry planes
    integer, allocatable              :: symplane(:)                    !! Symmetry planes (internal identifiers)
    !
    ! Connectivity to the design mesh
    !
    integer                           :: dm_mode                        !! 0 (default - isoparametric), 1 (parametization through design mesh)
    integer                           :: dm_n_elements                  !! Number of elements (0, 1 or 2)
    integer, allocatable              :: dm_element(:)                  !! Element (internal identifier)
    real(kind=real64), allocatable    :: dm_xi_gn(:,:,:)                !! Local coordinates of the nodes of the elements at each element
    logical                           :: sensitivity                    !! False if dxda is zero in it, true otherwise.
    !
    ! Connectivity to internal points (if an internal element)
    !
    integer                           :: n_internalpoints               !! Number of internal points
    integer, allocatable              :: internalpoint(:)               !! Internal identifier of the internal point
    !
    ! Geometrical data at the nodes of the element
    !
    real(kind=real64), allocatable    :: xi_gn(:,:)                     !! Reference coordinate: xi_gn(coordinate,node)
    real(kind=real64), allocatable    :: x_gn(:,:)                      !! Position vector: x_gn(coordinate,node)
    real(kind=real64), allocatable    :: n_gn(:,:)                      !! Unit normal: n_gn(coordinate,node)
    real(kind=real64), allocatable    :: t1_gn(:,:)                     !! Unit tangent 1 (xi1 curvilinear coordinate): t1_gn(coordinate,node)
    real(kind=real64), allocatable    :: t2_gn(:,:)                     !! Unit tangent 2 (xi2 curvilinear coordinate): t2_gn(coordinate,node)
    real(kind=real64), allocatable    :: t3_gn(:,:)                     !! Unit tangent 3 (xi3 curvilinear coordinate): t3_gn(coordinate,node)
    real(kind=real64), allocatable    :: tbp_gn(:,:)                    !! Unit tangent at element boundary with n+: tbp_gn(coordinate,node). For BE free-term calculation and other things ...
    real(kind=real64), allocatable    :: tbm_gn(:,:)                    !! Unit tangent at element boundary with n-: tbm_gn(coordinate,node). For BE free-term calculation and other things ...
    real(kind=real64), allocatable    :: jacobian_gn(:)                 !! Jacobian at each geometrical node: xi_gn(coordinate,node)
    logical                           :: plane                          !! True if the element is plane (null curvature)
    !
    ! Geometrical data of the element
    !
    real(kind=real64), allocatable    :: centroid(:)                    !! Centroid
    real(kind=real64), allocatable    :: bbox(:,:)                      !! Bounding box of the element
    real(kind=real64), allocatable    :: bball_centre(:)                !! Centre of a ball that contains the element
    real(kind=real64)                 :: bball_radius                   !! Radius of a ball that contains the element
    !
    ! For BE
    !
    ! Different functional and geometrical interpolation
    ! Nota: se supone que aunque se utilicen distintos tipos de interpolaciones para cada cosa, todos tienen el mismo numero de nodos,
    ! y la interpolacion funcional tiene la misma delta_f. Delta_g se asume 0, elementos geometricamente continuos (nodos)
    ! las rutinas de los modulos, estan generalizadas, por lo que hace falta es que el programa principal los tome en cuenta (lo
    ! cual no es tan sencillo)
    !
    integer                           :: type_g                         !! Geometrical interpolation
    integer                           :: type_f1                        !! Functional interpolation for primary variables
    integer                           :: type_f2                        !! Functional interpolation for secondary variables
    integer                           :: type_f1f                       !! Functional interpolation for tau (poroelastic media)
    integer                           :: type_f2f                       !! Functional interpolation for Un (poroelastic media)
    logical                           :: discontinuous                  !! True if the functional interpolation is discontinuous
    real(kind=real64)                 :: delta_f                        !! Shifting of functional nodes towards the centre of the element when discontinuous
    ! Variables needed for efficient integration
    real(kind=real64)                 :: csize                          !! Characteristic length (maximum length of all edges) --> change name by clength
    integer                           :: n_phi                          !! Number of quadrature points needed to integrate any phi
    ! Vectors at functional nodes
    real(kind=real64), allocatable    :: xi_fn(:,:)                     !! Reference coordinate: xi_fn(coordinate,node)
    real(kind=real64), allocatable    :: x_fn(:,:)                      !! Position vector: x_fn(coordinate,node)
    real(kind=real64), allocatable    :: n_fn(:,:)                      !! Unit normal: n_fn(coordinate,node)
    ! Vectors at SBIE collocation points
    real(kind=real64), allocatable    :: xi_i_sbie(:,:)                 !! Reference coordinate: xi_i(coordinate,node)
    real(kind=real64), allocatable    :: x_i_sbie(:,:)                  !! Position vector: x_i(coordinate,node)
    real(kind=real64), allocatable    :: n_i_sbie(:,:)                  !! Unit normal: n_i(coordinate,node)
    ! Vectors at SBIE MCA collocation points
    real(kind=real64), allocatable    :: delta_sbie_mca(:)              !! Colloc. point displacement towards the inside of the element for SBIE MCA
    real(kind=real64), allocatable    :: xi_i_sbie_mca(:,:)             !! Reference coordinate: xi_i(coordinate,node)
    real(kind=real64), allocatable    :: x_i_sbie_mca(:,:)              !! Position vector: x_i(coordinate,node)
    real(kind=real64), allocatable    :: n_i_sbie_mca(:,:)              !! Unit normal: n_i(coordinate,node)
    ! Vectors at HBIE collocation points
    real(kind=real64), allocatable    :: delta_hbie(:)                  !! Colloc. point displacement towards the inside of the element for HBIE.
    real(kind=real64), allocatable    :: xi_i_hbie(:,:)                 !! Reference coordinate: xi_i(coordinate,node)
    real(kind=real64), allocatable    :: x_i_hbie(:,:)                  !! Position vector: x_i(coordinate,node)
    real(kind=real64), allocatable    :: n_i_hbie(:,:)                  !! Unit normal: n_i(coordinate,node)
    !
    ! Related to sensitivity analysis
    !
    ! Local coordinate of the element of the design mesh at the collocation points
    real(kind=real64), allocatable    :: dm_xi_i_sbie(:,:)              !! Reference coordinate: xi_i(coordinate,node)
    real(kind=real64), allocatable    :: dm_xi_i_sbie_mca(:,:)          !! Reference coordinate: xi_i(coordinate,node)
    real(kind=real64), allocatable    :: dm_xi_i_hbie(:,:)              !! Reference coordinate: xi_i(coordinate,node)
    ! Design velocities at collocation points
    real(kind=real64), allocatable    :: dxda_i_sbie(:,:,:)             !! dxda at xi_i: dxda_i(coordinate,node,design_variable)
    real(kind=real64), allocatable    :: dxda_i_sbie_mca(:,:,:)         !! dxda at xi_i: dxda_i(coordinate,node,design_variable)
    real(kind=real64), allocatable    :: dxda_i_hbie(:,:,:)             !! dxda at xi_i: dxda_i(coordinate,node,design_variable)
    ! Gradient of design velocities of the unit normal at collocation point
    real(kind=real64), allocatable    :: d2xdadx_i_sbie(:,:,:,:)        !! d2xdadx at xi_i: d2xdadx_i(i,j,node,design_variable)
    real(kind=real64), allocatable    :: d2xdadx_i_sbie_mca(:,:,:,:)    !! d2xdadx at xi_i: d2xdadx_i(i,j,node,design_variable)
    real(kind=real64), allocatable    :: d2xdadx_i_hbie(:,:,:,:)        !! d2xdadx at xi_i: d2xdadx_i(i,j,node,design_variable)
    ! Design variation of the unit normal at collocation point
    real(kind=real64), allocatable    :: dnda_i_sbie(:,:,:)             !! dnda at xi_i: dnda_i(coordinate,node,design_variable)
    real(kind=real64), allocatable    :: dnda_i_sbie_mca(:,:,:)         !! dnda at xi_i: dnda_i(coordinate,node,design_variable)
    real(kind=real64), allocatable    :: dnda_i_hbie(:,:,:)             !! dnda at xi_i: dnda_i(coordinate,node,design_variable)
    ! Gradient of design velocities of the unit normal at collocation point
    real(kind=real64), allocatable    :: d3xdadxdx_i_sbie(:,:,:,:,:)    !! d3xdadxdx at xi_i: d3xdadxdx_i(i,j,k,node,design_variable)
    real(kind=real64), allocatable    :: d3xdadxdx_i_sbie_mca(:,:,:,:,:)!! d3xdadxdx at xi_i: d3xdadxdx_i(i,j,k,node,design_variable)
    real(kind=real64), allocatable    :: d3xdadxdx_i_hbie(:,:,:,:,:)    !! d3xdadxdx at xi_i: d3xdadxdx_i(i,j,k,node,design_variable)
    !
    ! For coupled BE body load elements
    !
    ! element()%type is the base interpolation in the length direction (isoparametric) (where original discretization nodes refer):
    !   -2D tip load: point1 (0D)
    !   -3D edge load: line2/line3/etc... (1D).
    !   -3D tip load: point1 (0D) (modelo Luis)
    ! element()%type_g is the geometrical interpolation in the thickness direction (additional nodes):
    !   -For now always line2 (1D)
    ! element()%type_f2 is the functional interpolation in the thickness direction (additional nodes):
    !   -point1/line2/line3/etc... (1D)
    ! [Coupled BE body load elements]: Vectors at SBIE collocation points , associated with type_f2 functional interpolation of the load
    integer, allocatable              :: cbl_n_cp(:)                    !! Number of collocation points per node: cbl_n(node)
    real(kind=real64), allocatable    :: cbl_xi_i(:,:,:)                !! Reference coordinate: cbl_xi_i(coordinate,collocation point,node)
    real(kind=real64), allocatable    :: cbl_x_i(:,:,:)                 !! Position vector: cbl_x_i(coordinate,collocation point,node)
    !
    ! For FE
    !
    ! DOF management
    integer                           :: n_dof                          !! Number of DOF of the element
    integer, allocatable              :: node_n_dof(:)                  !! Number of DOF of each element node
    !
    ! Masas puntuales a traves de un elemento 0D
    !
    real(kind=real64), allocatable    :: mass_matrix(:,:)               !! Mass matrix for point loads (defined as in the nodal axes)
    integer                           :: fe_type                        !! For 1D FE:
                                                                        !!   (0) Degenerated beam (degbeam)
                                                                        !!   (1) Euler-Bernoulli beam (strbeam_eb)
                                                                        !!   (2) Timoshenko beam (strbeam_t)
                                                                        !!   (3) Bar (bar)
                                                                        !!   (4) Discrete translational spring-dashpot (distra)
                                                                        !!   (5) Discrete rotational/translational spring-dashpot (disrotra)
                                                                        !!
                                                                        !! For 2D FE in 3D:
                                                                        !!   (0) Degenerated shell
                                                                        !!
                                                                        !! For 3D FE:
                                                                        !!   (0) Solid
    integer                           :: fe_options(1)                  !! Finite element options:
                                                                        !!   - strbeam:
                                                                        !!     fe_options(1): 0 (line3 with no mid-node rotation), 1 (line3 with mid-node rotation)
                                                                        !!   - distra and disrotra
                                                                        !!     fe_options(1): 1 (coupled electromechanical spring (moving coil) distra_em and disrotra_em)
    !
    ! Integration setup for finite elements
    !
    logical                           :: mitc                           !! True for using MITC elements (where it corresponds)
    integer                           :: K_intmode                      !! Integration mode (K matrix): 0 (full), (1) reduced, (2) selective, (3) user-defined
    integer                           :: K_intngp(3)                    !! Number of Gauss points for user-defined integration mode (K matrix).
    integer                           :: M_intmode                      !! Integration mode (M matrix): 0 (full), (1) reduced, (2) selective, (3) user-defined
    integer                           :: M_intngp(3)                    !! Number of Gauss points for user-defined integration mode (M matrix).
    integer                           :: Q_intmode                      !! Integration mode (Q matrix): 0 (full), (1) reduced, (2) selective, (3) user-defined
    integer                           :: Q_intngp(3)                    !! Number of Gauss points for user-defined integration mode (Q matrix).
    real(kind=real64), allocatable    :: ep(:,:)                        !! Element local orientation:
                                                                        !!   For 1D FE: Orientation of x', y' and z' system
                                                                        !!   For 2D FE in 3D: orientation of x' for calculating stress resultants
    logical                           :: flooded
    real(kind=real64)                 :: A_flooded(3)
    real(kind=real64)                 :: rho_flooded
    logical                           :: submerged
    real(kind=real64)                 :: A_submerged(3)
    real(kind=real64)                 :: rho_submerged
    integer                           :: cs_type                        !! Cross-section type: 0: generic, 1: circle, 2: hollow circle, 3: rectangular
    real(kind=real64)                 :: cs_param(2)                    !! Cross-section parameters
    real(kind=real64)                 :: A                              !! Area (bar, strbeam)
    real(kind=real64)                 :: I(6)                           !! Inertia tensor (strbeam): Ix'x',Iy'y',Iz'z',Ix'y',Iy'z',Ix'z'
    real(kind=real64)                 :: ksh(3)                         !! Shear correction factors:
                                                                        !!   -strbeam : kt  ky'  kz'
                                                                        !!   -degbeam : -   ky'  kz'
                                                                        !!   -degshell: -   -    kz'
    real(kind=real64)                 :: k_r(8)                         !! Stiffnesses for discrete elements (real):
                                                                        !!   -distra 2D  :  kx' ky' -    -      -    -    -      -
                                                                        !!   -distra 3D  :  kx' ky' kz'  -      -    -    -      -
                                                                        !!   -disrotra 2D:  kx' ky' krz' ky'rz' -    -    -      -
                                                                        !!   -disrotra 3D:  kx' ky' kz'  krx'   kry' krz' ky'rz' kz'ry'
    complex(kind=real64)              :: k_c(8)                         !! Stiffnesses for discrete elements (complex, for time-harmonic):
                                                                        !!   -distra 2D  :  kx' ky' -    -      -    -    -      -
                                                                        !!   -distra 3D  :  kx' ky' kz'  -      -    -    -      -
                                                                        !!   -disrotra 2D:  kx' ky' krz' ky'rz' -    -    -      -
                                                                        !!   -disrotra 3D:  kx' ky' kz'  krx'   kry' krz' ky'rz' kz'ry'
    real(kind=real64)                 :: c(8)                           !! Dashpot for discrete elements:
                                                                        !!   -distra 2D  :  kx' ky' -    -      -    -    -      -
                                                                        !!   -distra 3D  :  kx' ky' kz'  -      -    -    -      -
                                                                        !!   -disrotra 2D:  kx' ky' krz' ky'rz' -    -    -      -
                                                                        !!   -disrotra 3D:  kx' ky' kz'  krx'   kry' krz' ky'rz' kz'ry'
    !
    complex(kind=real64)              :: em_U                           !! For distra_em y disrotra_em: Usource, Bl, R, L (electromag. parameters)
    real(kind=real64)                 :: em_Bl, em_R, em_L
    ! Cross-section properties for degenerated beam and shell FE
    real(kind=real64), allocatable    :: v_midnode(:,:,:)               !! For degenerated elements: v_1, v_2 and v_3 director vectors which define local rotations and thicknesses direction.
    real(kind=real64), allocatable    :: tn_midnode(:,:)                !! For degenerated elements: thicknesses orthogonal to the mid-surface or mid-line.
    real(kind=real64), allocatable    :: tv_midnode(:,:)                !! For degenerated elements: thickness at each node along the corresponding director vector (tv=tn/(n·v_j))
    real(kind=real64)                 :: orthotropic_shell_fd1(3)       !! Fiber direction 1 of orthotropic shells
    ! BE line loads
    real(kind=real64)                 :: r_integration                  !! Radius of integration for coupled BE line loads
    ! Conditions over FE elements
    integer, allocatable              :: ctype(:,:)                     !! Type of condition: ctype(1,1) (0 or 1, load in local coordinates), ctype(2,1) (0 or 1, load in global global)
    integer, allocatable              :: cvalue_i(:,:,:)                !! Condition values (integer): cvalue_i(coordinate,node,type)
    real(kind=real64), allocatable    :: cvalue_r(:,:,:)                !! Condition values (real   ): cvalue_r(coordinate,node,type)
    complex(kind=real64), allocatable :: cvalue_c(:,:,:)                !! Condition values (complex): cvalue_c(coordinate,node,type)
    ! Values of variables
    real(kind=real64), allocatable    :: value_r(:,:,:)                 !! Condition values (real   ): value_r(coordinate,node,variable): variable==1 (stresses or stress resultants), variable==2 (equilibrating forces/moments)
    complex(kind=real64), allocatable :: value_c(:,:,:)                 !! Condition values (complex): value_c(coordinate,node,variable)
    !
    ! Incident field
    !
    complex(kind=real64), allocatable :: incident_c(:,:,:)              !! Incident field at the functional nodes of the element: incident_c(dof,node,face)
    !
    ! Resultants
    !
    real(kind=real64), allocatable    :: resultants_r(:,:,:)            !! Resultants values (real   ): resultants_r(component,variable,group)
    complex(kind=real64), allocatable :: resultants_c(:,:,:)            !! Resultants values (complex): resultants_c(component,variable,group)
    !
    ! Export settings
    !
    logical                           :: export                         !! Export flag
  end type fbem_element
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! PART (MESH PARTITION)
  !
  type fbem_part
    integer                           :: id                             !! External identifier
    character(len=fbem_stdcharlen)    :: name                           !! Name
    integer                           :: type                           !! Type of part / physical entity
    integer                           :: entity                         !! Physical entity indentifier (internal identifier)
    integer                           :: n_dimension                    !! Dimension of the elements it contains (1 (lines), 2 (surfaces), 3 (volumes))
    integer                           :: n_elements                     !! Number of elements
    integer, allocatable              :: element(:)                     !! Elements (internal identifiers)
    integer                           :: n_nodes                        !! Number of nodes
    integer, allocatable              :: node(:)                        !! Nodes (internal identifiers)
    logical                           :: export                         !! Export flag
  end type fbem_part
  ! Internal identifier for the type of part
  integer, parameter                  :: fbem_part_be_boundary =1       !! BE boundary
  integer, parameter                  :: fbem_part_fe_subregion=2       !! FE subregion
  integer, parameter                  :: fbem_part_be_bodyload =3       !! BE body load
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! BE BODY LOAD
  !
  !! Body load data structure
  type fbem_be_bodyload
    integer                           :: id                             !! External identifier
    character(len=fbem_stdcharlen)    :: name                           !! Name
    integer                           :: part                           !! Part (internal identifier)
    integer                           :: region                         !! Region (internal identifier)
    integer                           :: coupling                       !! Type of coupling: 0 (UNCOUPLED), 1 (WITH BEAM TIP), 3 (WITH SHELL EDGE), 2 (WITH BEAM LINE), 4 (WITH SHELL SURFACE)
    logical                           :: export                         !! Export flag
  end type fbem_be_bodyload
  integer, parameter                  :: fbem_bl_coupling_beam_tip     =1 !!
  integer, parameter                  :: fbem_bl_coupling_beam_line    =2 !!
  integer, parameter                  :: fbem_bl_coupling_shell_edge   =3 !!
  integer, parameter                  :: fbem_bl_coupling_shell_surface=4 !!
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! BOUNDARY
  !
  !! Boundary data structure
  type fbem_boundary
    !
    ! Basic data
    !
    integer                           :: id                             !! External identifier
    character(len=fbem_stdcharlen)    :: name                           !! Name
    integer                           :: class                          !! Class of boundary
    integer                           :: coupling                       !! Coupling type
    !
    ! Connectivity
    !
    ! With regions of the boundary
    integer                           :: n_regions                      !! Number of regions
    integer, allocatable              :: region(:)                      !! Regions (internal identifiers)
    integer, allocatable              :: region_boundary_idx(:)         !! Local index of the boundary in each region
    logical, allocatable              :: region_boundary_reversion(:)   !! Local index of the boundary in each region
    ! With part
    integer                           :: part                           !! Part (internal identifier)
    !
    ! Conditions over BE boundaries
    !
    integer, allocatable              :: ctype(:,:)                     !! Type of condition: ctype(coord,side)
    integer, allocatable              :: cvalue_i(:,:,:)                !! Condition values (integer): cvalue_i(coord,param,side)
    real(kind=real64), allocatable    :: cvalue_r(:,:,:)                !! Condition values (real): cvalue_r(coord,param,side)
    complex(kind=real64), allocatable :: cvalue_c(:,:,:)                !! Condition values (complex): cvalue_c(coord,param,side)
    !
    ! Resultants
    !
    real(kind=real64), allocatable    :: resultants_r(:,:,:)            !! Resultants values (real   ): resultants_r(component,variable,group)
    complex(kind=real64), allocatable :: resultants_c(:,:,:)            !! Resultants values (complex): resultants_c(component,variable,group)
    !
    ! Export settings
    !
    logical                           :: export                         !! Export flag
  end type fbem_boundary
  !
  ! Associated parameters
  !
  ! Internal identifier for the class of boundary
  integer, parameter                  :: fbem_boundary_class_ordinary =1   !! Ordinary boundary (to be changed by "type")
  integer, parameter                  :: fbem_boundary_class_cracklike=2   !! Crack-like boundary (to be changed by "type")
  ! Internal identifier for the type of coupling of a boundary
  integer, parameter                  :: fbem_boundary_coupling_be      =1 !! The BE boundary is uncoupled
  integer, parameter                  :: fbem_boundary_coupling_be_be   =2 !! The boundary couples two BE regions
  integer, parameter                  :: fbem_boundary_coupling_be_fe   =3 !! The boundary is coupled with FE
  integer, parameter                  :: fbem_boundary_coupling_be_fe_be=4 !! The boundary couples two BE regions with a FE
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! MATERIAL
  !
  !! Material data structure
  type fbem_material
    integer                           :: id                             !! External identifier
    character(len=fbem_stdcharlen)    :: type                           !! Type of material
    logical, allocatable              :: property_defined(:)            !! True if the property has been defined
    real(kind=real64), allocatable    :: property(:,:)                  !! Matrix with all the properties of the material
  end type fbem_material
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! REGION (OR SUBREGION)
  !
  !! Region data structure
  type fbem_region
    !
    ! Basic data
    !
    integer                           :: id                             !! External identifier
    character(len=fbem_stdcharlen)    :: name                           !! Name
    integer                           :: class                          !! Class of region: BE or FE
    ! TO BE DEPRECATED
    integer                           :: type                           !! Type of media     (to be deprecated)
    integer                           :: subtype                        !! Subtype of media  (to be deprecated)
    real(kind=real64), allocatable    :: property_r(:)                  !! Media properties (real) (to be deprecated)
    complex(kind=real64), allocatable :: property_c(:)                  !! Media properties (complex) (to be deprecated)
    ! NEW PARADIGM
    integer                           :: model
    integer                           :: n_materials
    integer, allocatable              :: material(:)
    !
    ! Connectivity to the mesh
    !
    ! FE REGION
    ! With its subregions
    integer                           :: n_fe_subregions                !! Number of subregions
    integer, allocatable              :: fe_subregion(:)                !! Subregions of the region (internal identifiers)
    ! BE region
    ! Type of fundamental solution
    integer                           :: space                          !! Space for BE regions: full-space, half-space
    integer                           :: halfspace_n                    !! For a multilayered or homogenous half-space: unit normal direction, e.g. -1 => n=(-1,0,0), -3 => n=(0,0,-1)
    real(kind=real64)                 :: halfspace_x                    !! For a half-space: Coordinate of the hyperplane
    integer                           :: halfspace_bc                   !! For a half-space: Boundary condition at the hyperplane
    real(kind=real64), allocatable    :: multilayered_x(:)              !! For a multilayered half-space: coordinate of the hyperplane (x top of the layer), including the depth of the half-space
    ! With BE boundaries
    integer                           :: n_boundaries                   !! Number of BE boundaries
    integer, allocatable              :: boundary(:)                    !! BE boundaries of the region (internal identifiers)
    logical, allocatable              :: boundary_reversion(:)          !! Reversion of normals of BE boundaries
    ! With BE body loads
    integer                           :: n_be_bodyloads                 !! Number of BE body loads
    integer, allocatable              :: be_bodyload(:)                 !! BE body loads of the region (internal identifiers)
    ! With internal points
    integer                           :: n_internalpoints               !! Number of internal points
    integer, allocatable              :: internalpoint(:)               !! Internal points of the region (internal identifiers)
    ! With incident wave fields
    integer                           :: n_incidentfields               !! Number of incident wave fields
    integer, allocatable              :: incidentfield(:)               !! Incident wave fields of the region (internal identifiers)
    ! FE SUBREGION
    ! With its part
    integer                           :: part                           !! Part (internal identifier)
    ! With its region
    integer                           :: region                         !! Region (internal identifier)
    !
    ! RIGID REGIONS
    !
    real(kind=real64), allocatable    :: master_node_x(:)               !! Coordinates of the master node
    integer                           :: master_node                    !! 0: means that it is not rigid, >0 indicates the master node internal identifier
    !
    ! Connectivity to the design mesh
    !
    logical                           :: sensitivity                    !! False if dxda is zero in it, true otherwise.
    !
    ! Export settings
    !
    logical                           :: export                         !! Export flag
  end type fbem_region
  !
  ! Associated parameters
  !
  ! Internal identifier for the type of region
  integer, parameter                  :: fbem_potential   =1            !! If a Laplace PDE (static) or Helmholtz PDE (inviscid fluid)
  integer, parameter                  :: fbem_viscoelastic=2            !! If a viscoelastic solid
  integer, parameter                  :: fbem_poroelastic =3            !! If a poroelastic medium
  integer, parameter                  :: fbem_rigid       =4            !! If a rigid solid
  ! Internal identifier for the subtype of a viscoelastic region
  integer, parameter                  :: fbem_elastic=1                 !! If a viscoelastic with null hysteretic damping
  integer, parameter                  :: fbem_bardet =2                 !! If a Bardet's model poroelastic->viscoelastic
  integer, parameter                  :: fbem_brt_cp1=3                 !! If a Bougacha-Roesset-Tassoulas poro->elastic with cp=cp1
  integer, parameter                  :: fbem_brt_cp2=4                 !! If a Bougacha-Roesset-Tassoulas poro->elastic with cp=cp2
  integer, parameter                  :: fbem_brt_cpm=5                 !! If a Bougacha-Roesset-Tassoulas poro->elastic with cp=(cp1+cp2)/2
  integer, parameter                  :: fbem_elastic_orthotropic=6     !! If a orthotropic elastic
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! INTERNAL POINT
  !
  !! Internal point data structure
  type fbem_internalpoint
    !
    ! Basic data
    !
    integer                           :: id                             !! External identifier
    real(kind=real64), allocatable    :: x(:)                           !! Position vector
    real(kind=real64), allocatable    :: n(:)                           !! Unit normal
    real(kind=real64), allocatable    :: value_r(:,:)                   !! Values of variables (real): value_r(coord,variable)
    complex(kind=real64), allocatable :: value_c(:,:)                   !! Values of variables (complex): value_c(coord,variable)
    complex(kind=real64), allocatable :: incident_c(:,:)                !! Values of incident field (complex): incident_c(coord,variable)
    ! Sensitivity of the variables of the node
    real(kind=real64), allocatable    :: dvda_r(:,:,:)                  !! Derivative of variables with respect to design variables (real): dvda_r(coord,variable,design_variable)
    complex(kind=real64), allocatable :: dvda_c(:,:,:)                  !! Derivative of variables with respect to design variables (complex): dvda_c(coord,variable,design_variable)
    !
    ! Connectivity to the mesh
    !
    ! With a region
    integer                           :: region                         !! Region (internal identifier)
    !
    ! Connectivity to the design mesh
    !
    integer                           :: dm_n_elements                  !! Number of elements
    integer, allocatable              :: dm_element(:)                  !! Elements (internal identifiers)
    integer, allocatable              :: dm_element_d(:)                !! Dimension of the elements (internal identifiers)
    real(kind=real64), allocatable    :: dm_element_xi(:,:)             !! Local coordinates at each element
    real(kind=real64), allocatable    :: dxda(:,:)                      !! dx/da     [delta x_{i}  ]: dxda(component,designvariable)
    real(kind=real64), allocatable    :: d2xdadx(:,:,:)                 !! d/dx_{j}(dx_{i}/da)      : d2xdadx(i,j,designvariable)
    real(kind=real64), allocatable    :: dnda(:,:,:)                    !! dn_{k}/da [delta n_{k}]  : dnda(component,axis,designvariable)
    logical                           :: sensitivity                    !! False only if dx/da is 0 for all design variables
    !
    ! Export settings
    !
    logical                           :: export                         !! Export flag
  end type fbem_internalpoint
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! INCIDENT FIELD
  !
  !! Incident wave field data structure
  type fbem_incidentfield
    !
    ! Basic data
    !
    integer                           :: id                             !! External identifier
    integer                           :: class                          !! Class of wave: point wave, line wave or plane wave.
    integer                           :: space                          !! Space of the source: full-space, half-space or multilayered half-space
    ! For half-spaces:
    integer                           :: np                             !! Direction of the plane unit normal.
    real(kind=real64)                 :: xp                             !! Coordinate of the plane.
    integer                           :: bc                             !! Boundary condition on the plane.
    ! For multilayered half-spaces:
    integer                           :: n_layers                       !! Number of layers
    real(kind=real64), allocatable    :: layer_ztop(:)                  !! Coordinate of the layer's top surface
    integer, allocatable              :: layer_material(:)              !! Material (internal identifier)
    complex(kind=real64), allocatable :: layer_properties(:,:)
    integer, allocatable              :: layer_ctype(:)
    complex(kind=real64), allocatable :: layer_cvalue(:)
    !
    complex(kind=real64), allocatable :: layer_lambda(:)                !!
    complex(kind=real64), allocatable :: layer_mu(:)                    !!
    real(kind=real64), allocatable    :: layer_rho(:)                   !!
    complex(kind=real64), allocatable :: layer_amplitudes(:,:)          !! Waves amplitudes at each layer
    !
    complex(kind=real64), allocatable :: layer_A(:,:,:,:)               !! Wave amplitudes for each variable mode, variable and layer
    complex(kind=real64), allocatable :: layer_k(:,:)                   !! Wavenumber of each layer
    ! For all:
    complex(kind=real64)              :: amplitude                      !! Amplitude of the incident wave
    integer                           :: variable                       !! 0 if primary variables, or 1 if secondary variables
    real(kind=real64), allocatable    :: x0(:)                          !! Origin of the incident wave (reference point coordinates)
    real(kind=real64)                 :: varphi                         !! Angle varphi
    real(kind=real64)                 :: theta                          !! Angle theta
    real(kind=real64), allocatable    :: xs(:)                          !! Origin of the symmetrical/antisymmetrical decomposition
    real(kind=real64), allocatable    :: symconf(:)                     !! Decomposition for each direction: -1 (anti-symmetry), 0 (total), 1 (symmetry)
    integer                           :: region_type                    !! Type of region: fluid, viscoelastic, poroelastic.
    integer                           :: wave_type                      !! Type of wave: eg a elastic plane wave could be a P, SH,..
    !
    ! Connectivity to regions
    !
    !integer                           :: n_regions
    !integer, allocatable              :: region()
  end type fbem_incidentfield
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! GROUP OF MESH/SUBMESH OBJECTS
  !
  !! Group of mesh/submesh objects
  type fbem_group
    integer                           :: id                             !! External identifier
    integer                           :: type                           !! Type of objects this group stores
    integer                           :: n_objects                      !! Number of objects
    integer, allocatable              :: object(:)                      !! Internal identifier of objects
  end type fbem_group
  integer, parameter                  :: fbem_group_type_nodes   =1     !! Group of nodes
  integer, parameter                  :: fbem_group_type_elements=2     !! Group of elements
  integer, parameter                  :: fbem_group_type_subedges=3     !! Group of subedges
  ! ================================================================================================================================

contains

  ! Build node->nodes connectivity
  subroutine fbem_node_nodes_connectivity(rn,n_nodes,node,tol)
    implicit none
    ! I/O
    integer              :: rn             !! Dimension of the space
    integer              :: n_nodes        !! Number of nodes
    type(fbem_node)      :: node(n_nodes)  !! List of nodes
    real(kind=real64)    :: tol            !! Tolerance (radius of the ball)
    ! Local
    integer              :: kni, knj       ! Counters
    integer              :: snj            ! Identifier
    real(kind=real64)    :: rv(rn)         ! Distance vector
    real(kind=real64)    :: r              ! Distance
    integer, allocatable :: node_nodes(:)  ! node -> nodes connectivity by distance
    !
    ! Initialise
    !
    allocate (node_nodes(n_nodes))
    do kni=1,n_nodes
      node(kni)%n_nodes=0
    end do
    !
    ! Build the connectivity vector of nodes for each node (distance symmetry property is fully exploited)
    !
    do kni=1,n_nodes
      if (node(kni)%n_nodes.eq.0) then
        node_nodes=0
        do knj=kni+1,n_nodes
          rv=node(kni)%x-node(knj)%x
          r=sqrt(dot_product(rv,rv))
          if (r.le.tol) then
            node(kni)%n_nodes=node(kni)%n_nodes+1
            node_nodes(node(kni)%n_nodes)=knj
          end if
        end do
        if (node(kni)%n_nodes.gt.0) then
          ! Allocate and copy the connectivity vector node i -> node j to the data structure
          allocate (node(kni)%node(node(kni)%n_nodes))
          node(kni)%node=node_nodes(1:node(kni)%n_nodes)
          ! Copy the connectivity node j -> node i
          do knj=1,node(kni)%n_nodes
            snj=node_nodes(knj)
            if (node(snj)%n_nodes.eq.0) then
              node(snj)%n_nodes=node(kni)%n_nodes
              allocate (node(snj)%node(node(snj)%n_nodes))
              node(snj)%node=node_nodes(1:node(snj)%n_nodes)
              ! Replace the appearance of node j in node_nodes by node i
              node(snj)%node(knj)=kni
            else
              call fbem_error_message(error_unit,0,'node',node(kni)%id,&
              'the set of nodes such that dist(x(node)-x)<=gtol must be the same for all nodes in the set: check gtol and/or mesh')
            end if
          end do
        end if
      end if
    end do
    !
    ! Finalise
    !
    deallocate (node_nodes)
  end subroutine fbem_node_nodes_connectivity

  ! Build node->elements connectivity
  subroutine fbem_node_elements_connectivity(n_nodes,node,n_elements,element)
    implicit none
    ! I/O
    integer            :: n_nodes
    type(fbem_node)    :: node(n_nodes)
    integer            :: n_elements
    type(fbem_element) :: element(n_elements)
    ! Local
    integer            :: kni, ke, se, knj
    integer            :: nd                  ! Number of dimensions of the element
    integer            :: n_nd_elements(3)    ! Number of 1d, 2d and 3d elements connected to a node
    ! Loop through nodes
    do kni=1,n_nodes
      !
      ! Calculate the number of elements that are connected with the node
      !
      ! Initialize
      node(kni)%n_elements=0
      ! Loop through elements
      do ke=1,n_elements
        ! Loop through the nodes of the element
        do knj=1,element(ke)%n_nodes
          if (element(ke)%node(knj).eq.kni) node(kni)%n_elements=node(kni)%n_elements+1
        end do
      end do
      !
      ! An unconnected node is not allowed
      !
      if (node(kni)%n_elements.eq.0) then
        call fbem_error_message(error_unit,0,'node',node(kni)%id,'this node is not connected any element')
      end if
      !
      ! Build the connectivity vectors
      !
      ! Allocate the connectivity vectors
      allocate (node(kni)%element(node(kni)%n_elements))
      allocate (node(kni)%element_node_iid(node(kni)%n_elements))
      allocate (node(kni)%element_node_loctype(node(kni)%n_elements))
      ! Initialize
      node(kni)%n_elements=0
      ! Loop through elements
      do ke=1,n_elements
        do knj=1,element(ke)%n_nodes
          if (element(ke)%node(knj).eq.kni) then
            node(kni)%n_elements=node(kni)%n_elements+1
            node(kni)%element(node(kni)%n_elements)=ke
            node(kni)%element_node_iid(node(kni)%n_elements)=knj
            if (element(ke)%type.ne.fbem_point1) then
              node(kni)%element_node_loctype(node(kni)%n_elements)=fbem_node_loctype(knj,element(ke)%type)
            else
              node(kni)%element_node_loctype(node(kni)%n_elements)=fbem_loctype_vertex
            end if
          end if
        end do
      end do
      !
      ! Determine what classes (line, surface, volume) of elements are connected to the node
      !
      ! Count the number of 1D (line), 2D (surface) and 3D (volume) elements connected to a node
      n_nd_elements=0
      do ke=1,node(kni)%n_elements
        se=node(kni)%element(ke)
        nd=element(se)%n_dimension
        if (element(se)%type.ne.fbem_point1) n_nd_elements(nd)=n_nd_elements(nd)+1
      end do
      ! Possible dimensional situations of the node connected with its elements:
      ! F F F = 0: Unconnected node.
      ! T F F = 1: Only connected to 1D elements.
      ! F T F = 2: Only connected to 2D elements.
      ! T T F = 3: Connected to 1D and 2D elements.
      ! F F T = 4: Only connected to 3D elements.
      ! T F T = 5: Connected to 1D and 3D elements.
      ! F T T = 6: Connected to 2D and 3D elements.
      ! T T T = 7: Connected to 1D, 2D and 3D elements.
      node(kni)%dimensional_degree=0
      do nd=1,3
        if (n_nd_elements(nd).gt.0) node(kni)%dimensional_degree=node(kni)%dimensional_degree+2**(nd-1)
      end do
      ! We interchange the code 3 and 4, to have more coherent codes.
      if (node(kni)%dimensional_degree.eq.3) then
        node(kni)%dimensional_degree=4
      else
        if (node(kni)%dimensional_degree.eq.4) node(kni)%dimensional_degree=3
      end if
      ! Definitive code for dimensional_degree:
      ! 0: Unconnected node
      ! 1: Only connected to 1D elements.
      ! 2: Only connected to 2D elements.
      ! 3: Only connected to 3D elements.
      ! 4: Connected to 1D and 2D elements, i.e. the node is connected with >=2 elements which are 1D and 2D elements.
      ! 5: Connected to 1D and 3D elements, i.e. the node is connected with >=2 elements which are 1D and 3D elements.
      ! 6: Connected to 2D and 3D elements, i.e. the node is connected with >=2 elements which are 2D and 3D elements.
      ! 7: Connected to 1D, 2D and 3D elements, i.e. the node is connected with >=3 elements which are 1D, 2D and 3D elements.
    end do
  end subroutine fbem_node_elements_connectivity

  !! Build connectivity node -> parts
  subroutine fbem_node_parts_connectivity(n_nodes,node,n_elements,element,n_parts,part)
    implicit none
    ! I/O
    integer              :: n_nodes
    type(fbem_node)      :: node(n_nodes)
    integer              :: n_elements
    type(fbem_element)   :: element(n_elements)
    integer              :: n_parts
    type(fbem_part)      :: part(n_parts)
    ! Local
    integer              :: kni, ke, se, tmp_n_parts, kp, kpj
    integer, allocatable :: tmp_part(:)
    !
    ! Find the maximum number of elements to which any node is connected and allocate
    !
    tmp_n_parts=0
    do kni=1,n_nodes
      if (node(kni)%n_elements.gt.tmp_n_parts) then
        tmp_n_parts=node(kni)%n_elements
      end if
    end do
    allocate (tmp_part(tmp_n_parts))
    !
    ! Build connectivity
    !
    do kni=1,n_nodes
      tmp_part=0
      kp=0
      do ke=1,node(kni)%n_elements
        se=node(kni)%element(ke)
        kp=kp+1
        tmp_part(kp)=element(se)%part
      end do
      do kp=1,tmp_n_parts
        if (tmp_part(kp).ne.0) then
          do kpj=kp+1,tmp_n_parts
            if (tmp_part(kp).eq.tmp_part(kpj)) tmp_part(kpj)=0
          end do
        end if
      end do
      node(kni)%n_parts=0
      do kp=1,tmp_n_parts
        if (tmp_part(kp).ne.0) then
          node(kni)%n_parts=node(kni)%n_parts+1
        end if
      end do
      allocate (node(kni)%part(node(kni)%n_parts))
      kpj=0
      do kp=1,tmp_n_parts
        if (tmp_part(kp).ne.0) then
          kpj=kpj+1
          node(kni)%part(kpj)=tmp_part(kp)
        end if
      end do
    end do
  end subroutine fbem_node_parts_connectivity

  !! Build connectivity part -> nodes
  subroutine fbem_part_nodes_connectivity(n_nodes,node,n_parts,part)
    implicit none
    ! I/O
    integer              :: n_nodes
    type(fbem_node)      :: node(n_nodes)
    integer              :: n_parts
    type(fbem_part)      :: part(n_parts)
    ! Local
    integer              :: kp, kn, kpj, knj
    do kp=1,n_parts
      if (part(kp)%entity.gt.0) then
        part(kp)%n_nodes=0
        do kn=1,n_nodes
          do kpj=1,node(kn)%n_parts
            if (node(kn)%part(kpj).eq.kp) then
              part(kp)%n_nodes=part(kp)%n_nodes+1
            end if
          end do
        end do
        allocate (part(kp)%node(part(kp)%n_nodes))
        knj=0
        do kn=1,n_nodes
          do kpj=1,node(kn)%n_parts
            if (node(kn)%part(kpj).eq.kp) then
              knj=knj+1
              part(kp)%node(knj)=kn
            end if
          end do
        end do
        !write(*,*) 'part', part(kp)%id, 'has', part(kp)%n_nodes, 'nodes'
      end if
    end do
  end subroutine fbem_part_nodes_connectivity

  ! Build the connectivity of nodes with symmetry planes.
  subroutine fbem_node_symplanes_connectivity(n_nodes,node,n_symplanes,symplane_eid,tol)
    implicit none
    ! I/O
    integer           :: n_nodes
    type(fbem_node)   :: node(n_nodes)
    integer           :: n_symplanes
    integer           :: symplane_eid(3)
    real(kind=real64) :: tol
    ! Local
    integer           :: i, j, k
    ! Build connectivity
    do i=1,n_nodes
      ! Initialize
      node(i)%n_symplanes=0
      do j=1,n_symplanes
        ! Check if the coordinate of the node in the direction of the normal of the symmetry plane is < tolerance
        if (abs(node(i)%x(symplane_eid(j))).le.tol) then
          ! Increment the counter
          node(i)%n_symplanes=node(i)%n_symplanes+1
        end if
      end do
      if (node(i)%n_symplanes.gt.0) then
        ! Allocate the vector of symmetry planes
        allocate (node(i)%symplane(node(i)%n_symplanes))
        ! Initialize
        k=0
        do j=1,n_symplanes
          ! Check if the coordinate of the node in the direction of the normal of the symmetry plane is < tolerance
          if (abs(node(i)%x(symplane_eid(j))).le.tol) then
            ! Increment the counter
            k=k+1
            ! Save the symmetry plane
            node(i)%symplane(k)=j
          end if
        end do
      end if
    end do
  end subroutine fbem_node_symplanes_connectivity

  ! Build the connectivity of elements with symmetry planes.
  subroutine fbem_element_symplanes_connectivity(n,n_nodes,node,n_elements,element)
    implicit none
    ! I/O
    integer            :: n
    integer            :: n_nodes
    type(fbem_node)    :: node(n_nodes)
    integer            :: n_elements
    type(fbem_element) :: element(n_elements)
    ! Local
    integer :: ke, kn, sn, ksp, kk, e_symplane(3), n_symplane(3)
    ! Build connectivity
    do ke=1,n_elements
      element(ke)%n_symplanes=0
      if (element(ke)%n_dimension.eq.n) cycle
      e_symplane=1
      do kn=1,element(ke)%n_nodes
        sn=element(ke)%node(kn)
        n_symplane=0
        do ksp=1,node(sn)%n_symplanes
          n_symplane(node(sn)%symplane(ksp))=1
        end do
        do ksp=1,3
          if ((e_symplane(ksp)+n_symplane(ksp)).ne.2) e_symplane(ksp)=0
        end do
        if (sum(e_symplane).eq.0) exit
      end do
      element(ke)%n_symplanes=sum(e_symplane)
      if (element(ke)%n_symplanes.gt.0) then
        allocate (element(ke)%symplane(element(ke)%n_symplanes))
        kk=1
        do ksp=1,3
          if (e_symplane(ksp).eq.1) then
            element(ke)%symplane(kk)=ksp
            kk=kk+1
          end if
        end do
      end if
    end do
  end subroutine fbem_element_symplanes_connectivity

  ! Check the position of all nodes with respect to the symmetry planes
  subroutine fbem_check_nodes_symplanes_configuration(n_nodes,node,n_symplanes,symplane_eid)
    implicit none
    ! I/O
    integer           :: n_nodes
    type(fbem_node)   :: node(n_nodes)
    integer           :: n_symplanes
    integer           :: symplane_eid(3)
    ! Local
    integer           :: ks, i, inisign
    ! Check that the connectivity is correct with respect to the existent symmetry configuration
    do ks=1,n_symplanes
      inisign=0
      do i=1,n_nodes
        if (node(i)%n_symplanes.eq.0) then
          if (inisign.eq.0) then
            if (node(i)%x(symplane_eid(ks)).gt.0.d0) then
              inisign=1
            else
              inisign=-1
            end if
          else
            if (node(i)%x(symplane_eid(ks)).gt.0.d0) then
              if (inisign.ne.1) then
                call fbem_error_message(error_unit,0,'node',node(i)%id,&
                                        'the node position is not compatible with the symmetry configuration')
              end if
            else
              if (inisign.ne.(-1)) then
                call fbem_error_message(error_unit,0,'node',node(i)%id,&
                                        'the node position is not compatible with the symmetry configuration')
              end if
            end if
          end if
        end if
      end do
    end do
  end subroutine fbem_check_nodes_symplanes_configuration

  !! Build the sub-elements of the mesh (edges of line/surface/volume elements, and faces of surface/volume elements),
  !! it has to be further developed
  subroutine fbem_build_mesh_subelements(n_nodes,node,n_elements,element,n_subedges,subedge,n_subfaces,subface)
    implicit none
    ! I/O
    integer                         :: n_nodes
    type(fbem_node)                 :: node(n_nodes)
    integer                         :: n_elements
    type(fbem_element)              :: element(n_elements)
    integer                         :: n_subedges
    type(fbem_element), allocatable :: subedge(:)
    integer                         :: n_subfaces
    type(fbem_element), allocatable :: subface(:)
    ! Local
    integer                         :: kei, sei, kd, kdi, sndi1, sndi2, kej, sej, kdj, sndj1, sndj2, kn
    !
    ! Build the faces. In the future when 3D solid FEM elements
    !
!    do kei=1,n_elements
!      if (element(kei)%n_dimension.lt.2) cycle
!      element(kei)%n_subfaces=fbem_n_faces(element(kei)%type)
!      allocate (element(kei)%subface(element(kei)%n_subfaces),element(kei)%subface_reversion(element(kei)%n_subfaces))
!      element(kei)%subface=0
!      element(kei)%subface_reversion=.false.
!    end do
    !
    ! Build the edges ¿?¿?¿?¿ This must be modified when faces 3D elements are implemented ?¿?¿?¿?¿ I don't think so.... a lo mejor incluso
    ! eliminar su calculo porque parece no ser necesario...
    !
    ! First calculate the number of edges and stores the element to edge information
    ! Allocate needed data structures
    do kei=1,n_elements
      if (element(kei)%n_dimension.lt.1) cycle
      element(kei)%n_subedges=fbem_n_edges(element(kei)%type)
      allocate (element(kei)%subedge(element(kei)%n_subedges),element(kei)%subedge_reversion(element(kei)%n_subedges))
      element(kei)%subedge=0
      element(kei)%subedge_reversion=.false.
    end do
    ! Start
    kd=0
    do kei=1,n_elements
      if (element(kei)%type.eq.fbem_point1) cycle
      do kdi=1,element(kei)%n_subedges
        ! Select an edge that has not been studied
        if (element(kei)%subedge(kdi).eq.0) then
          ! Increase counter and save the edge in the kei element
          kd=kd+1
          element(kei)%subedge(kdi)=kd
          element(kei)%subedge_reversion(kdi)=.false.
          ! Vertex nodes of edge kd
          sndi1=element(kei)%node(fbem_edge_node(1,kdi,element(kei)%type))
          sndi2=element(kei)%node(fbem_edge_node(2,kdi,element(kei)%type))
          ! Loop through the elements of the node sndi1
          do kej=1,node(sndi1)%n_elements
            sej=node(sndi1)%element(kej)
            ! Explore its edges
            do kdj=1,element(sej)%n_subedges
              ! If the edge has not been studied, study it
              if (element(sej)%subedge(kdj).eq.0) then
                ! Vertex nodes of edge kdj of the element sej
                sndj1=element(sej)%node(fbem_edge_node(1,kdj,element(sej)%type))
                sndj2=element(sej)%node(fbem_edge_node(2,kdj,element(sej)%type))
                ! Non-inverted edge
                if ((sndj1.eq.sndi1).and.(sndj2.eq.sndi2)) then
                  if (fbem_edge_type(kdj,element(sej)%type).eq.fbem_edge_type(kdi,element(kei)%type)) then
                    element(sej)%subedge(kdj)=kd
                    element(sej)%subedge_reversion(kdj)=.false.
                  else
                    call fbem_warning_message(error_unit,0,'element',sej,&
                    'an edge of this element shares vertices with other edge but they are not of the same type.')
                  end if
                  exit
                end if
                ! Inverted edge
                if ((sndj2.eq.sndi1).and.(sndj1.eq.sndi2)) then
                  if (fbem_edge_type(kdj,element(sej)%type).eq.fbem_edge_type(kdi,element(kei)%type)) then
                    element(sej)%subedge(kdj)=kd
                    element(sej)%subedge_reversion(kdj)=.true.
                  else
                    call fbem_warning_message(error_unit,0,'element',sej,&
                    'an edge of this element shares vertices with other edge but they are not of the same type.')
                  end if
                end if
              end if
            end do
          end do
        end if
      end do
    end do
    ! Now that we now the number of subedges, repeat it but now saving edge to element information
    n_subedges=kd
    allocate (subedge(n_subedges))
    ! Initialize
    do kd=1,n_subedges
      subedge(kd)%id=0
      subedge(kd)%type=0
      subedge(kd)%n_dimension=0
      subedge(kd)%size=0.d0
      subedge(kd)%part=0
      subedge(kd)%n_nodes=0
      subedge(kd)%element=0
      subedge(kd)%n_subedges=0
      subedge(kd)%n_subfaces=0
      subedge(kd)%n_supelements=0
      subedge(kd)%type_g=0
      subedge(kd)%type_f1=0
      subedge(kd)%type_f2=0
      subedge(kd)%type_f1f=0
      subedge(kd)%type_f2f=0
      subedge(kd)%discontinuous=.false.
      subedge(kd)%delta_f=0.d0
      subedge(kd)%csize=0.d0
      subedge(kd)%n_phi=0
      subedge(kd)%plane=.false.
      subedge(kd)%bball_radius=0.d0
      subedge(kd)%A=0.d0
      subedge(kd)%I=0.d0
      subedge(kd)%export=.true.
    end do
    ! Build
    do kd=1,n_subedges
      subedge(kd)%id=kd
      ! Find the edge in the element list
      subedge(kd)%n_supelements=0
      do kei=1,n_elements
        do kdi=1,element(kei)%n_subedges
          if (element(kei)%subedge(kdi).eq.kd) then
            subedge(kd)%n_supelements=subedge(kd)%n_supelements+1
            exit
          end if
        end do
      end do
      ! Allocate connectivity
      allocate (subedge(kd)%supelement(subedge(kd)%n_supelements),subedge(kd)%supelement_element_iid(subedge(kd)%n_supelements))
      subedge(kd)%n_supelements=0
      do kei=1,n_elements
        do kdi=1,element(kei)%n_subedges
          if (element(kei)%subedge(kdi).eq.kd) then
            subedge(kd)%n_supelements=subedge(kd)%n_supelements+1
            subedge(kd)%supelement(subedge(kd)%n_supelements)=kei
            subedge(kd)%supelement_element_iid(subedge(kd)%n_supelements)=kdi
            exit
          end if
        end do
      end do
      ! Allocate connectivity to nodes, using the edge in the first element
      sei=subedge(kd)%supelement(1)
      do kdi=1,element(sei)%n_subedges
        if (element(sei)%subedge(kdi).eq.kd) then
          subedge(kd)%type=fbem_edge_type(kdi,element(sei)%type)
          subedge(kd)%n_nodes=fbem_n_nodes(subedge(kd)%type)
          allocate (subedge(kd)%node(subedge(kd)%n_nodes))
          if (element(sei)%subedge_reversion(kdi).eqv.(.false.)) then
            do kn=1,subedge(kd)%n_nodes
              subedge(kd)%node(kn)=element(sei)%node(fbem_edge_node(kn,kdi,element(sei)%type))
            end do
          else
            do kn=1,subedge(kd)%n_nodes
              subedge(kd)%node(kn)=fbem_element_invert(element(sei)%node(fbem_edge_node(kn,kdi,element(sei)%type)),element(sei)%type)
            end do
          end if
          exit
        end if
      end do
    end do
    ! Write
!    write(*,*) 'Edge -> Nodes'
!    do kd=1,n_subedges
!      write(*,*) subedge(kd)%id, (subedge(kd)%node(kn),kn=1,subedge(kd)%n_nodes)
!    end do
!    write(*,*) 'Edge -> Elements'
!    do kd=1,n_subedges
!      write(*,*) subedge(kd)%id, (subedge(kd)%supelement(kei),kei=1,subedge(kd)%n_supelements), (subedge(kd)%supelement_element_iid(kei),kei=1,subedge(kd)%n_supelements)
!    end do
!    write(*,*) 'Elements -> Edges'
!    do kei=1,n_elements
!      write(*,*) element(kei)%id, (element(kei)%subedge(kd),kd=1,element(kei)%n_subedges), (element(kei)%subedge_reversion(kd),kd=1,element(kei)%n_subedges)
!    end do
!    stop
  end subroutine fbem_build_mesh_subelements

  !! Check and characterize the topology of the BE mesh
  subroutine fbem_check_and_characterize_be_mesh(rn,n_nodes,node,n_elements,element,n_subedges,subedge,n_parts,part)
    implicit none
    ! I/O
    integer             :: rn
    integer             :: n_nodes
    type(fbem_node)     :: node(n_nodes)
    integer             :: n_elements
    type(fbem_element)  :: element(n_elements)
    integer             :: n_subedges
    type(fbem_element)  :: subedge(n_subedges)
    integer             :: n_parts
    type(fbem_part)     :: part(n_parts)
    ! Local
    integer             :: kn, sn, kd, se1, se2, kp, ke, kpj, spj
    !
    ! Check that each BE BOUNDARY and BE BODY LOAD has unique nodes (except nodes belonging to parts with type 0 which are used for
    ! auxiliar purposes). Also, initialize "in_boundary" and "is_singular" flags.
    !
    do kp=1,n_parts
      if (part(kp)%entity.gt.0) then
        if ((part(kp)%type.eq.fbem_part_be_boundary).or.(part(kp)%type.eq.fbem_part_be_bodyload)) then
          do kn=1,part(kp)%n_nodes
            sn=part(kp)%node(kn)
            do kpj=1,node(sn)%n_parts
              spj=node(sn)%part(kpj)
              if (spj.eq.kp) cycle
              if (part(spj)%type.ne.0) then
                call fbem_error_message(error_unit,0,'node',node(sn)%id,&
                'this boundary element node must belong to only one "boundary" or "body load".')
              end if
            end do
          end do
          node(sn)%in_boundary=.false.
          node(sn)%is_singular=.false.
        end if
      end if
    end do
    !
    ! Check that the nodes of BE elements cannot have more than 1 symmetry plane in 2D, and 2 symmetry planes in  3D
    !
    !
    ! Check and characterize BE mesh
    !
    do kp=1,n_parts
      if ((part(kp)%type.eq.fbem_part_be_boundary).or.(part(kp)%type.eq.fbem_part_be_bodyload)) then
        select case (element(part(kp)%element(1))%n_dimension)
          ! 1D ELEMENTS
          case (1)
            do kn=1,part(kp)%n_nodes
              sn=part(kp)%node(kn)
              select case (node(sn)%n_elements)
                case (1)
                  if (node(sn)%element_node_loctype(1).eq.fbem_loctype_vertex) node(sn)%in_boundary=.true.
                case (2)
                  se1=node(sn)%element(1)
                  se2=node(sn)%element(2)
                  ! If connected to two elements, the node must be a vertex for both elements
                  if ((node(sn)%element_node_loctype(1).ne.fbem_loctype_vertex).or.&
                      (node(sn)%element_node_loctype(2).ne.fbem_loctype_vertex)) then
                    call fbem_error_message(error_unit,0,'node',node(sn)%id,'is connected wrongly with elements.')
                  end if
                  ! Check orientation
                  if (node(sn)%element_node_iid(1).eq.node(sn)%element_node_iid(2)) then
                    call fbem_warning_message(error_unit,0,'node',node(sn)%id,'connected to the following elements')
                    call fbem_warning_message(error_unit,0,'element',element(se1)%id,'wrong orientation')
                    call fbem_error_message  (error_unit,0,'element',element(se2)%id,'wrong orientation')
                  end if
                case default
                  call fbem_error_message(error_unit,0,'node',node(sn)%id,'is a singular node.')
              end select
            end do
          ! 2D ELEMENTS
          case (2)
            do kd=1,n_subedges
              if (element(subedge(kd)%supelement(1))%part.eq.kp) then
                select case (subedge(kd)%n_supelements)
                  case (1)
                    do kn=1,subedge(kd)%n_nodes
                      node(subedge(kd)%node(kn))%in_boundary=.true.
                    end do
                  case (2)
                    se1=subedge(kd)%supelement(1)
                    se2=subedge(kd)%supelement(2)
                    if (element(se1)%subedge_reversion(subedge(kd)%supelement_element_iid(1)).eqv.&
                        element(se2)%subedge_reversion(subedge(kd)%supelement_element_iid(2))) then
                      call fbem_warning_message(error_unit,0,'element',element(se1)%id,'wrong orientation')
                      call fbem_error_message  (error_unit,0,'element',element(se2)%id,'wrong orientation')
                    end if
                  case default
                    call fbem_error_message(error_unit,0,'element',element(subedge(kd)%supelement(1))%id,'has a singular edge')
                end select
              end if
            end do
        end select
      end if
    end do
  ! ?¿?¿??¿¿?¿?¿?
  ! Chequear que las BE body loads no traspasan los contornos, y que las body loads con la misma dimension que los BE boundary
  ! elements no coinciden....
  ! ?¿?¿??¿¿?¿?¿?
  ! ?¿?¿??¿¿?¿?¿?
  ! Check that the BE boundaries are connected together only through their own boundaries, this is related with geometric tolerance,
  ! and that this connection is conform (if required by SBIE collocation, otherwise SBIE MCA or discontinuous elements must be used)
  ! ?¿?¿??¿¿?¿?¿?
!    !
!    ! ?????????????????????????????????????
!    ! Check that the orientation of the boundaries in regions are compatible
!    ! ??????????????????????????????????????
!    !
!    ! Write nodes at the boundary of the mesh
!    do i=1,n_parts
!      if (part(i)%type.eq.fbem_part_be_boundary) then
!        do j=1,part(i)%n_nodes
!          if (node(part(i)%node(j))%in_boundary.eqv.(.true.)) then
!            write(10,'(i11,3e25.16)') node(part(i)%node(j))%id, node(part(i)%node(j))%x
!          else
!            write(11,'(i11,3e25.16)') node(part(i)%node(j))%id, node(part(i)%node(j))%x
!          end if
!        end do
!      end if
!    end do
!    stop
  end subroutine fbem_check_and_characterize_be_mesh

  !! Check and characterize the topology of the FE mesh
  subroutine fbem_check_and_characterize_fe_mesh(rn,n_nodes,node,n_elements,element,n_subedges,subedge,n_parts,part)
    implicit none
    ! I/O
    integer             :: rn
    integer             :: n_nodes
    type(fbem_node)     :: node(n_nodes)
    integer             :: n_elements
    type(fbem_element)  :: element(n_elements)
    integer             :: n_subedges
    type(fbem_element)  :: subedge(n_subedges)
    integer             :: n_parts
    type(fbem_part)     :: part(n_parts)
    ! Local
    integer             :: ke, se, kn, sn, kd, sd, se1, se2, kp

    ! falta....
    ! -Study symmetry configuration conditions for FE elements


    !
    ! By default all the nodes are not in the boundary of the mesh and are not singular
    !
    do kp=1,n_parts
      ! antes estaba fbem_part_be_boundary
      if (part(kp)%type.eq.fbem_part_fe_subregion) then
        do kn=1,part(kp)%n_nodes
          sn=part(kp)%node(kn)
          node(sn)%in_boundary=.false.
          node(sn)%is_singular=.false.
        end do
      end if
    end do
    !
    ! Assign values to node()%is_singular
    ! -----------------------------------
    !
    do kp=1,n_parts
      if (part(kp)%type.eq.fbem_part_fe_subregion) then
        do ke=1,part(kp)%n_elements
          se=part(kp)%element(ke)
          select case (element(se)%n_dimension)
            !
            ! Point elements elements
            !
            case (0)
              ! Nothing to do
            !
            ! Beam elements
            !
            case (1)
              do kn=1,element(se)%n_nodes
                sn=element(se)%node(kn)
                select case (node(sn)%dimensional_degree)
                  ! Connected to only 1D elements
                  case (1)
                    select case (node(sn)%n_elements)
                      case (1)
                        if (node(sn)%element_node_loctype(1).eq.fbem_loctype_vertex) node(sn)%in_boundary=.true.
                      case (2)
                        ! If connected to two elements, the node must be a vertex for both elements
                        if ((node(sn)%element_node_loctype(1).ne.fbem_loctype_vertex).or.&
                            (node(sn)%element_node_loctype(2).ne.fbem_loctype_vertex)) then
                          call fbem_error_message(error_unit,0,'node',node(sn)%id,'invalid mesh near this node.')
                        end if
                      case default
                        ! A joint between >2 elements (singular point)
                        node(sn)%is_singular=.true.
                    end select
                  ! Connected to 1D and 2D elements
                  case (4)
                    node(sn)%is_singular=.true.
                  ! Connected to 1D and 3D elements
                  case (5)
                    call fbem_error_message(error_unit,0,'node',node(sn)%id,'connection of 1D and 3D FEs not considered yet.')
                  ! Connected to 1D, 2D and 3D elements
                  case (7)
                    call fbem_error_message(error_unit,0,'node',node(sn)%id,'connection of 1D, 2D and 3D FEs not considered yet.')
                  case default
                    call fbem_error_message(error_unit,0,'node',node(sn)%id,'invalid mesh near this node.')
                end select
              end do
            !
            ! Solid elements in 2D case or Shell elements in 3D case
            !
            case (2)
              do kd=1,element(se)%n_subedges
                sd=element(se)%subedge(kd)
                do kn=1,subedge(sd)%n_nodes
                  sn=subedge(sd)%node(kn)
                  select case (node(sn)%dimensional_degree)
                    ! Connected to only 2D elements
                    case (2)
                      select case (subedge(sd)%n_supelements)
                        case (1)
                          node(sn)%in_boundary=.true.
                        case (2)
                          se1=subedge(sd)%supelement(1)
                          se2=subedge(sd)%supelement(2)
                          ! Check that the edge has opposite orientations for both elements
                          if (element(se1)%subedge_reversion(subedge(sd)%supelement_element_iid(1)).eqv.&
                              element(se2)%subedge_reversion(subedge(sd)%supelement_element_iid(2))) then
                            call fbem_error_message(error_unit,0,'element',element(se1)%id,&
                                                    'wrong orientation of this element.')
                          end if
                        case default
                          ! A joint between >2 elements (the node belong to a singular edge)
                          node(sn)%is_singular=.true.
                      end select
                    ! Connected to 1D and 2D elements
                    case (4)
                      node(sn)%is_singular=.true.
                    ! Connected to 2D and 3D elements
                    case (6)
                      call fbem_error_message(error_unit,0,'node',node(sn)%id,'connection of 2D and 3D FEs not considered yet.')
                    ! Connected to 1D, 2D and 3D elements
                    case (7)
                      call fbem_error_message(error_unit,0,'node',node(sn)%id,'connection of 1D, 2D and 3D FEs not considered yet.')
                    case default
                      call fbem_error_message(error_unit,0,'node',node(sn)%id,'invalid mesh near this node.')
                  end select
                end do
              end do
            case default
              call fbem_error_message(error_unit,0,'node',node(sn)%id,'other n_dimension for FE not considered yet.')
          end select
        end do
      end if
    end do
  end subroutine fbem_check_and_characterize_fe_mesh

end module fbem_data_structures
