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
!! <b> This subroutine reads the input file and processes the data. </b>
!!
!! The input file is divided into sections. Each section starts with a line "[section name]". Within each section, the data can be
!! arranged in two ways:
!!   - Using assignment statements: "variable name = values", or "variable name : values"; one per line.
!!   - Using a ordered sequence of numbers and strings.

subroutine read_input_file

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_geometry
  use fbem_string_handling
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! No implicit variables are allowed
  implicit none

  ! I/O variables
  ! None

  ! Local variables
  integer                               :: input_fileunit
  integer                               :: aux_fileunit
  integer                               :: iostat_var
  character(len=fbem_string_max_length) :: iomsg_var

  ! ===============
  ! OPEN INPUT FILE
  ! ===============

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'OPENING input file "'//trim(input_filename)//'"')
  input_fileunit=fbem_get_valid_unit()
  open(unit=input_fileunit,file=input_filename,action='read',status='old',recl=fbem_file_record_length,iostat=iostat_var,iomsg=iomsg_var)
  if (iostat_var.ne.0) call fbem_error_message(error_unit,0,'iostat',iostat_var,trim(iomsg_var))

  ! ===============
  ! READ INPUT FILE
  ! ===============

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'START reading input file')

  call read_commands(input_fileunit)

  call read_problem(input_fileunit)

  call read_settings(input_fileunit)

  call read_materials(input_fileunit)

  if (problem%analysis.eq.fbem_harmonic) then

    if (len_trim(frequencies_filename).eq.0) then
      call read_frequencies(input_fileunit)
    else
      if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'OPENING frequency file "'//trim(frequencies_filename)//'"')
      aux_fileunit=fbem_get_valid_unit()
      open(unit=aux_fileunit,file=frequencies_filename,action='read',status='old',recl=fbem_file_record_length,iostat=iostat_var,iomsg=iomsg_var)
      if (iostat_var.ne.0) call fbem_error_message(error_unit,0,'iostat',iostat_var,trim(iomsg_var))
      call read_frequencies(aux_fileunit)
      if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'CLOSING frequency file "'//trim(frequencies_filename)//'"')
      close(unit=aux_fileunit,iostat=iostat_var,iomsg=iomsg_var)
      if (iostat_var.ne.0) call fbem_error_message(error_unit,0,'iostat',iostat_var,trim(iomsg_var))
    end if

    call read_incident_mechanics_harmonic(input_fileunit)

  end if

  ! High level entities (physical entities)

  call read_regions(input_fileunit)

  call read_boundaries(input_fileunit)

  call read_fe_subregions(input_fileunit)

  call read_be_bodyloads(input_fileunit)

  ! Low level entities (mesh objects)

  select case (mesh_file_mode)
    case (0)
      call read_parts(input_fileunit,1)
      call read_elements(input_fileunit,1)
      call read_nodes(input_fileunit,1)
    case (1,2)
      if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'OPENING mesh file "'//trim(mesh_filename)//'"')
      aux_fileunit=fbem_get_valid_unit()
      open(unit=aux_fileunit,file=mesh_filename,action='read',status='old',recl=fbem_file_record_length,iostat=iostat_var,iomsg=iomsg_var)
      if (iostat_var.ne.0) call fbem_error_message(error_unit,0,'iostat',iostat_var,trim(iomsg_var))
      call read_parts(aux_fileunit,mesh_file_mode)
      call read_elements(aux_fileunit,mesh_file_mode)
      call read_nodes(aux_fileunit,mesh_file_mode)
      if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'CLOSING mesh file "'//trim(mesh_filename)//'"')
      close(unit=aux_fileunit,iostat=iostat_var,iomsg=iomsg_var)
      if (iostat_var.ne.0) call fbem_error_message(error_unit,0,'iostat',iostat_var,trim(iomsg_var))
    case  default
      call fbem_error_message(error_unit,1,'mesh_file_mode',0,'wrong type of mesh mode')
  end select

  call fbem_node_nodes_connectivity(problem%n,n_nodes,node,geometric_tolerance)

  call fbem_node_elements_connectivity(n_nodes,node,n_elements,element)

  call fbem_node_parts_connectivity(n_nodes,node,n_elements,element,n_parts,part)

  call fbem_part_nodes_connectivity(n_nodes,node,n_parts,part)

  call fbem_build_mesh_subelements(n_nodes,node,n_elements,element,n_subedges,subedge,n_subfaces,subface)

  ! Internal points and elements

  call read_internal(input_fileunit)

  ! Sensitivity analysis

  if (problem%sensitivity) call read_sensitivity(input_fileunit)

  ! Other data

  !call read_discontinuous_boundary_elements(input_fileunit)

  !call read_special_boundary_elements(input_fileunit)

  call read_cross_sections(input_fileunit)

  call read_symmetry_planes(input_fileunit)

  call fbem_node_symplanes_connectivity(n_nodes,node,n_symplanes,symplane_eid,geometric_tolerance)

  call fbem_element_symplanes_connectivity(problem%n,n_nodes,node,n_elements,element)

  if (collapse_nodal_pos) call fbem_transformation_collapse_nodal_positions(problem%n,n_nodes,node,n_symplanes,symplane_eid)

  call fbem_check_nodes_symplanes_configuration(n_nodes,node,n_symplanes,symplane_eid)

  if (n_be_regions.gt.0) then
    call fbem_check_and_characterize_be_mesh(problem%n,n_nodes,node,n_elements,element,n_subedges,subedge,n_parts,part)
  end if

  if (n_fe_regions.gt.0) then
    call fbem_check_and_characterize_fe_mesh(problem%n,n_nodes,node,n_elements,element,n_subedges,subedge,n_parts,part)
  end if

  call read_groups(input_fileunit)

  if (n_be_regions.gt.0) then

    call build_connectivities_be_bodyload_elements_and_fe

    call build_connectivities_be_and_fe

  end if

  call build_connectivities_rigid_regions

  ! ================================================================================================================================
  ! READ FEM NODE AND ELEMENT OPTIONS
  !
  call read_fem_node_options(input_fileunit)

  call read_element_options(input_fileunit)
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! PROCESS CONDITIONS AT BEM MESH
  !
  !   1. Assign default conditions over BEM boundaries.
  !   2. Read conditions over BEM boundaries.
  !   3. Transfer conditions over BEM boundaries to its nodes.
  !   4. Read conditions over BEM nodes.

  !   5. ?¿?¿??¿?Leer tambien las cargas de volumen..... haciendo lo pispo .....?¿?¿??¿?

  !
  if (n_be_regions.gt.0) then
    select case (problem%type)
      case (fbem_laplace)
        call assign_default_conditions_bem_boundaries_laplace
        call read_conditions_bem_boundaries_laplace(input_fileunit)
        call transfer_conditions_bem_boundaries_laplace
        call read_conditions_bem_nodes_laplace(input_fileunit)
      case (fbem_mechanics)
        select case (problem%analysis)
          case (fbem_static)

            call assign_default_conditions_bem_boundaries_mechanics_static
            call read_conditions_bem_boundaries_mechanics_static(input_fileunit)
            call transfer_conditions_bem_boundaries_mechanics_static

            call assign_default_conditions_bem_bodyloads_mechanics_static
            call read_conditions_bem_bodyloads_mechanics_static(input_fileunit)
            call transfer_conditions_bem_bodyloads_mechanics_static

            ! To-do: introduce bem bodyload nodes conditions using groups
            call read_conditions_bem_nodes_mechanics_static(input_fileunit)

          case (fbem_harmonic)

            call assign_default_conditions_bem_boundaries_mechanics_harmonic
            call read_conditions_bem_boundaries_mechanics_harmonic(input_fileunit)
            call transfer_conditions_bem_boundaries_mechanics_harmonic

            call assign_default_conditions_bem_bodyloads_mechanics_harmonic
            call read_conditions_bem_bodyloads_mechanics_harmonic(input_fileunit)
            call transfer_conditions_bem_bodyloads_mechanics_harmonic

            ! To-do: introduce bem bodyload nodes conditions
            call read_conditions_bem_nodes_mechanics_harmonic(input_fileunit)

        end select
    end select
  end if
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! PROCESS CONDITIONS AT FEM MESH
  !
  !   1. Assign default conditions over FEM elements.
  !   2. Read conditions over FEM elements.
  !   3. Assign default conditions over FEM nodes.
  !   4. Read conditions over FEM nodes.
  !
  if (n_fe_regions.gt.0) then
    ! FEM elements
    call assign_default_conditions_fem_elements
    call read_conditions_fem_elements_mechanics(input_fileunit)
    ! FEM nodes
    call assign_default_conditions_fem_nodes
    select case (problem%type)
      case (fbem_laplace)
        ! Not yet
      case (fbem_mechanics)
        select case (problem%analysis)
          case (fbem_static)
            call read_conditions_fem_nodes_mechanics_static(input_fileunit)
          case (fbem_harmonic)
            call read_conditions_fem_nodes_mechanics_harmonic(input_fileunit)
      end select
    end select
    call assign_default_conditions_fem_nodes_symmetry_planes
  end if
  ! --------------------------------------------------------------------------------------------------------------------------------

  call assign_default_bem_formulation
  ! anadir en la rutina de leer, poder leerlo tambien para las body loads
  call read_bem_formulation_boundaries(input_fileunit)

  !if (n_be_regions.gt.0) call read_export_sif(input_fileunit)

  call read_export(input_fileunit)

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'END reading input file')

  ! ================
  ! CLOSE INPUT FILE
  ! ================

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'CLOSING input file "'//trim(input_filename)//'"')
  close(unit=input_fileunit,iostat=iostat_var,iomsg=iomsg_var)
  if (iostat_var.ne.0) call fbem_error_message(error_unit,0,'iostat',iostat_var,trim(iomsg_var))

end subroutine read_input_file
