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


subroutine read_conditions_bem_nodes_laplace(input_fileunit)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_symmetry
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! No implicit variables are allowed
  implicit none

  ! I/O variables
  integer                                 :: input_fileunit    ! Input file unit
  ! Local variables
  logical                                 :: found             ! Logical variable for sections and keywords
  integer                                 :: kb, sp
  integer                                 :: kn, sn
  character(len=fbem_stdcharlen)          :: keyword
  character(len=fbem_stdcharlen)          :: section_name

  ! Locate the section
  section_name='conditions over nodes'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'READING section ['//trim(section_name)//']')
  call fbem_search_section(input_fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! ==============================================================================================================================
    ! FIND GROUPS OF BEM NODES
    ! ==============================================================================================================================

    ! ¿?¿?¿?¿? falta esto.....

    ! ==============================================================================================================================
    ! FIND BEM NODES
    ! ==============================================================================================================================

    ! Loop through the nodes of each BE boundary
    do kb=1,n_boundaries
      sp=boundary(kb)%part
      do kn=1,part(sp)%n_nodes
        sn=part(sp)%node(kn)
        ! Locate the node
        call fbem_search_section(input_fileunit,section_name,found)
        write(keyword,*) 'node ', node(sn)%id
        call fbem_trim2b(keyword)
        call fbem_search_keyword(input_fileunit,keyword,':',found)
        if (found) then
          select case (boundary(kb)%coupling)

            ! ===========
            ! BE BOUNDARY
            ! ===========

            case (fbem_boundary_coupling_be)
                select case (boundary(kb)%class)

                  ! -----------------
                  ! ORDINARY BOUNDARY
                  ! -----------------

                  case (fbem_boundary_class_ordinary)
                    read(input_fileunit,*) node(sn)%ctype(1,1)
                    call fbem_search_section(input_fileunit,section_name,found)
                    call fbem_search_keyword(input_fileunit,keyword,':',found)
                    select case (node(sn)%ctype(1,1))
                      ! 0: p=P
                      ! 1: j=J=k·dp/dn
                      case (0,1)
                        read(input_fileunit,*) node(sn)%ctype(1,1), node(sn)%cvalue_r(1,1,1)
                    end select

                  ! -------------------
                  ! CRACK-LIKE BOUNDARY
                  ! -------------------

                  case (fbem_boundary_class_cracklike)
                    ! Face +
                    read(input_fileunit,*) node(sn)%ctype(1,1)
                    call fbem_search_section(input_fileunit,section_name,found)
                    call fbem_search_keyword(input_fileunit,keyword,':',found)
                    select case (node(sn)%ctype(1,1))
                      ! 0: p=P
                      ! 1: j=J=k·dp/dn
                      case (0,1)
                        read(input_fileunit,*) node(sn)%ctype(1,1), node(sn)%cvalue_r(1,1,1)
                    end select
                    ! Face -
                    read(input_fileunit,*) node(sn)%ctype(1,2)
                    backspace(input_fileunit)
                    select case (node(sn)%ctype(1,2))
                      ! 0: p=P
                      ! 1: j=J=k·dp/dn
                      case (0,1)
                        read(input_fileunit,*) node(sn)%ctype(1,2), node(sn)%cvalue_r(1,1,2)
                    end select

                end select

            ! ==============
            ! BE-FE BOUNDARY
            ! ==============

            case (fbem_boundary_coupling_be_fe)
              stop 'not implemented yet'

            ! =================
            ! BE-FE-BE BOUNDARY
            ! =================

            case (fbem_boundary_coupling_be_fe_be)
              stop 'not implemented yet'

          end select
        end if
      end do
    end do

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  end if

end subroutine read_conditions_bem_nodes_laplace
