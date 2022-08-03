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


subroutine read_conditions_fem_nodes_mechanics_static(input_fileunit)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_geometry
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! No implicit variables are allowed
  implicit none

  ! I/O variables
  integer                                 :: input_fileunit    ! Input file unit
  ! Local variables
  logical                                 :: found             ! Logical variable for sections and keywords
  integer                                 :: k
  integer                                 :: ks, sp
  integer                                 :: kn, sn
  character(len=fbem_stdcharlen)          :: keyword
  character(len=fbem_stdcharlen)          :: section_name
  integer                                 :: kg, kng
  real(kind=real64)                       :: center(3), axis(3), theta, x(3), urot(3)


  ! Locate the section
  section_name='conditions over nodes'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(input_fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! ==============================================================================================================================
    ! FIND GROUPS OF FEM NODES
    ! ==============================================================================================================================

    ! Loop through groups of FEM nodes
    do kg=1,n_groups
      ! Locate only group of nodes
      if (group(kg)%type.eq.fbem_group_type_nodes) then
        call fbem_search_section(input_fileunit,section_name,found)
        write(keyword,*) 'group ', group(kg)%id
        call fbem_trim2b(keyword)
        call fbem_search_keyword(input_fileunit,keyword,':',found)
        if (found) then
          ! Read only groups of nodes belonging to a FE subregion
          kn=group(kg)%object(1)
          sp=node(kn)%part(1)
          if (part(sp)%type.eq.fbem_part_fe_subregion) then

            ! ???? por ahora siempre se leen desplazamientos y rotaciones para todo.... incluso aunque el nodo no lo va a utilizar

            ! Read the condition for each degree of freedom
            do k=1,3*(problem%n-1)
              read(input_fileunit,*) node(kn)%ctype(k,1)
              if (k.eq.1) then
                call fbem_search_section(input_fileunit,section_name,found)
                call fbem_search_keyword(input_fileunit,keyword,':',found)
              else
                backspace(input_fileunit)
              end if
              ! Switch depending on the type of boundary condition
              select case (node(kn)%ctype(k,1))
                ! 0: Dirichlet (known displacement/rotation): value u_k
                ! 1: Neumann (known force/moment): value f_k
                case (0,1)
                  read(input_fileunit,*) node(kn)%ctype(k,1), node(kn)%cvalue_r(k,1,1)
                  ! Copy the condition to the rest of nodes
                  do kng=2,group(kg)%n_objects
                    sn=group(kg)%object(kng)
                    node(sn)%ctype(k,1)=node(kn)%ctype(k,1)
                    node(sn)%cvalue_r(k,1,1)=node(kn)%cvalue_r(k,1,1)
                  end do
                ! Infinitesimal rotation
                case (4)
                  ! Only for displacement DOF
                  if (k.le.problem%n) then
                    select case (problem%n)
                      case (2)
                        read(input_fileunit,*) node(kn)%ctype(k,1), center(1:2), theta
                        center(3)=0.d0
                        axis=[0.d0,0.d0,1.d0]
                      case (3)
                        read(input_fileunit,*) node(kn)%ctype(k,1), center, axis, theta
                    end select
                    axis=axis/sqrt(dot_product(axis,axis))
                    ! Copy the condition to all the nodes of the group
                    do kng=1,group(kg)%n_objects
                      sn=group(kg)%object(kng)
                      x=0.d0
                      x(1:problem%n)=node(sn)%x
                      urot=fbem_rotation_infinitesimal_displacement(center,axis,1.d0,x)
                      node(sn)%ctype(k,1)=0
                      node(sn)%cvalue_r(k,1,1)=urot(k)*theta
                    end do
                  else
                    call fbem_error_message(error_unit,0,'group',group(kg)%id,&
                                            'invalid type of boundary condition for this group of nodes')
                  end if

              end select
            end do

          end if
        end if
      end if
    end do

    ! ==============================================================================================================================
    ! FIND FEM NODES
    ! ==============================================================================================================================

    ! Loop through the nodes of each FE subregion
    do ks=1,n_fe_subregions
      sp=fe_subregion(ks)%part
      do kn=1,part(sp)%n_nodes
        sn=part(sp)%node(kn)
        call fbem_search_section(input_fileunit,section_name,found)
        write(keyword,*) 'node ', node(sn)%id
        call fbem_trim2b(keyword)
        call fbem_search_keyword(input_fileunit,keyword,':',found)
        if (found) then


          ! ???? por ahora siempre se leen desplazamientos y rotaciones para todo....

          ! Read the condition for each degree of freedom
          do k=1,3*(problem%n-1)
            read(input_fileunit,*) node(sn)%ctype(k,1)
            if (k.eq.1) then
              call fbem_search_section(input_fileunit,section_name,found)
              call fbem_search_keyword(input_fileunit,keyword,':',found)
            else
              backspace(input_fileunit)
            end if
            ! Switch depending on the type of boundary condition
            select case (node(sn)%ctype(k,1))
              ! Dirichlet (known displacement/rotation): value u_k
              case (0)
                read(input_fileunit,*) node(sn)%ctype(k,1), node(sn)%cvalue_r(k,1,1)
              ! Neumann (known force/moment): value f_k
              case (1)
                read(input_fileunit,*) node(sn)%ctype(k,1), node(sn)%cvalue_r(k,1,1)
            end select
          end do



        end if
      end do
    end do ! Loop through the nodes

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  end if

end subroutine read_conditions_fem_nodes_mechanics_static
