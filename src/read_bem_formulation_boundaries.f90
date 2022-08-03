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

subroutine read_bem_formulation_boundaries(input_fileunit)
  ! Fortran 2003 intrinsic module
  use iso_fortran_env
  ! fbem modules
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_data_structures
  ! Problem variables module
  use problem_variables
  ! No implicit variables are allowed
  implicit none
  ! I/O variables
  integer                               :: input_fileunit    !! Input file unit
  ! Local variables
  character(len=fbem_stdcharlen)        :: section_name      ! Name of the section
  logical                               :: found
  integer                               :: i
  character(len=fbem_stdcharlen)        :: keyword
  integer                               :: selected_n_nodes
  integer, allocatable                  :: selected_node(:)


  ! Nota: La dual B&M poroelastica no esta bien!!, la suma hay que ponderarla con otro numero de onda o hacerla de otra
  ! manera...


  ! Return if not needed
  if (n_be_regions.eq.0) return

  ! Locate the section
  section_name='bem formulation over boundaries'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(input_fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! Read formulation for BE nodes from boundaries
    do i=1,n_boundaries
      write(keyword,*) 'boundary ', boundary(i)%id
      call fbem_trim2b(keyword)
      call fbem_search_section(input_fileunit,section_name,found)
      call fbem_search_keyword(input_fileunit,keyword,':',found)
      if (found) then
        if (verbose_level.ge.2) write(*,*) 'boundary ', boundary(i)%id
        selected_n_nodes=part(boundary(i)%part)%n_nodes
        allocate (selected_node(selected_n_nodes))
        selected_node=part(boundary(i)%part)%node
        call read_bem_formulation_selected_nodes(input_fileunit,section_name,keyword,selected_n_nodes,selected_node)
        deallocate (selected_node)
      end if
    end do

    ! Read formulation for BE nodes from parts
    do i=1,n_parts
      write(keyword,*) 'part ', part(i)%id
      write(*,*) 'part ', part(i)%id
      call fbem_trim2b(keyword)
      call fbem_search_section(input_fileunit,section_name,found)
      call fbem_search_keyword(input_fileunit,keyword,':',found)
      if (found) then
        if (verbose_level.ge.2) write(*,*) 'part ', part(i)%id
        selected_n_nodes=part(i)%n_nodes
        allocate (selected_node(selected_n_nodes))
        selected_node=part(i)%node
        call read_bem_formulation_selected_nodes(input_fileunit,section_name,keyword,selected_n_nodes,selected_node)
        deallocate (selected_node)
      end if
    end do

    ! Read formulation for BE nodes from groups the integration point is not in the design macro-element
    do i=1,n_groups
      ! Process groups with BE nodes
      if (group(i)%type.ne.fbem_group_type_nodes) cycle
      if (part(node(group(i)%object(1))%part(1))%type.ne.fbem_part_be_boundary) cycle
      write(keyword,*) 'group ', group(i)%id
      call fbem_trim2b(keyword)
      call fbem_search_section(input_fileunit,section_name,found)
      call fbem_search_keyword(input_fileunit,keyword,':',found)
      if (found) then
        if (verbose_level.ge.3) write(*,*) 'group ', group(i)%id
        selected_n_nodes=group(i)%n_objects
        allocate (selected_node(selected_n_nodes))
        selected_node=group(i)%object
        call read_bem_formulation_selected_nodes(input_fileunit,section_name,keyword,selected_n_nodes,selected_node)
        deallocate (selected_node)
      end if
    end do

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  end if

end subroutine read_bem_formulation_boundaries
