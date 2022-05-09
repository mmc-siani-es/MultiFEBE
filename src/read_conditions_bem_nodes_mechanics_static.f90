! ---------------------------------------------------------------------
! Copyright (C) 2014-2022 Universidad de Las Palmas de Gran Canaria:
!                         Jacob D.R. Bordon
!                         Guillermo M. Alamo
!                         Luis A. Padron
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


subroutine read_conditions_bem_nodes_mechanics_static(input_fileunit)

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
  character(len=fbem_stdcharlen)          :: section_name
  logical                                 :: found             ! Logical variable for sections and keywords
  integer                                 :: k
  integer                                 :: kb, sp, sb, kng, kg
  integer                                 :: kn, sn
  character(len=fbem_stdcharlen)          :: keyword

  ! Locate the section
  section_name='conditions over nodes'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(input_fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! ==============================================================================================================================
    ! FIND GROUPS OF BEM NODES
    ! ==============================================================================================================================

    ! Loop through groups of BEM nodes
    do kg=1,n_groups
      ! Locate only group of nodes
      if (group(kg)%type.eq.fbem_group_type_nodes) then
        call fbem_search_section(input_fileunit,section_name,found)
        write(keyword,*) 'group ', group(kg)%id
        call fbem_trim2b(keyword)
        call fbem_search_keyword(input_fileunit,keyword,':',found)
        if (found) then
          ! The first node of the group
          kn=group(kg)%object(1)
          sp=node(kn)%part(1)
          ! Read only groups of nodes belonging to a BE boundary
          if (part(sp)%type.eq.fbem_part_be_boundary) then
            ! BE boundary of part
            sb=part(sp)%entity
            select case (boundary(sb)%coupling)

              ! ===========
              ! BE BOUNDARY
              ! ===========

              case (fbem_boundary_coupling_be)

                    select case (boundary(sb)%class)

                      ! -----------------
                      ! ORDINARY BOUNDARY
                      ! -----------------

                      case (fbem_boundary_class_ordinary)
                        ! Read the type of boundary condition for k direction
                        do k=1,problem%n
                          read(input_fileunit,*) node(kn)%ctype(k,1)
                          if (k.eq.1) then
                            call fbem_search_section(input_fileunit,section_name,found)
                            call fbem_search_keyword(input_fileunit,keyword,':',found)
                          else
                            backspace(input_fileunit)
                          end if
                          ! Switch depending on the type of boundary condition
                          select case (node(kn)%ctype(k,1))
                            ! 0. Global axes                         : u_k = U
                            ! 1. Global axes                         : t_k = T
                            ! 2. Local axes (l_1=n, l_2=t_1, l_3=t_2): u·l = U
                            ! 3. Local axes (l_1=n, l_2=t_1, l_3=t_2): t·l = T
                            case (0,1,2,3)
                              read(input_fileunit,*) node(kn)%ctype(k,1), node(kn)%cvalue_r(k,1,1)
                          end select
                        end do
                        ! Check that the B.C. are for local or global axes, but not mixed.
                        if ((node(kn)%ctype(1,1).eq.0).or.(node(kn)%ctype(1,1).eq.1)) then
                          do k=2,problem%n
                            if ((node(kn)%ctype(k,1).ne.0).and.(node(kn)%ctype(k,1).ne.1)) then
                              call fbem_error_message(error_unit,0,'node',node(kn)%id,&
                                                      'local axes B.C. and global axes B.C. can not be mixed.')
                            end if
                          end do
                        end if
                        if ((node(kn)%ctype(1,1).eq.2).or.(node(kn)%ctype(1,1).eq.3)) then
                          do k=2,problem%n
                            if ((node(kn)%ctype(k,2).ne.2).and.(node(kn)%ctype(k,1).ne.3)) then
                              call fbem_error_message(error_unit,0,'node',node(kn)%id,&
                                                      'local axes B.C. and global axes B.C. can not be mixed.')
                            end if
                          end do
                        end if
                        ! Copy the condition to the rest of nodes
                        do kng=2,group(kg)%n_objects
                          sn=group(kg)%object(kng)
                          node(sn)%ctype=node(kn)%ctype
                          node(sn)%cvalue_r=node(kn)%cvalue_r
                        end do

                      ! -------------------
                      ! CRACK-LIKE BOUNDARY
                      ! -------------------

                      case (fbem_boundary_class_cracklike)
                        ! Face +
                        ! Read the type of boundary condition for k direction
                        do k=1,problem%n
                          read(input_fileunit,*) node(kn)%ctype(k,1)
                          if (k.eq.1) then
                            call fbem_search_section(input_fileunit,section_name,found)
                            call fbem_search_keyword(input_fileunit,keyword,':',found)
                          else
                            backspace(input_fileunit)
                          end if
                          ! Switch depending on the type of boundary condition
                          select case (node(kn)%ctype(k,1))
                            ! Dirichlet (known displacement): value u_k
                            case (0)
                              read(input_fileunit,*) node(kn)%ctype(k,1), node(kn)%cvalue_r(k,1,1)
                            ! Neumann (known pressure): value p_k
                            case (1)
                              read(input_fileunit,*) node(kn)%ctype(k,1), node(kn)%cvalue_r(k,1,1)
                          end select
                        end do
                        ! Face -
                        ! Read the type of boundary condition for k direction
                        do k=1,problem%n
                          read(input_fileunit,*) node(kn)%ctype(k,2)
                          backspace(input_fileunit)
                          ! Switch depending on the type of boundary condition
                          select case (node(kn)%ctype(k,2))
                            ! Dirichlet (known displacement): value u_k
                            case (0)
                              read(input_fileunit,*) node(kn)%ctype(k,2), node(kn)%cvalue_r(k,1,2)
                            ! Neumann (known pressure): value p_k
                            case (1)
                              read(input_fileunit,*) node(kn)%ctype(k,2), node(kn)%cvalue_r(k,1,2)
                          end select
                        end do
                        ! Copy the condition to the rest of nodes
                        do kng=2,group(kg)%n_objects
                          sn=group(kg)%object(kng)
                          node(sn)%ctype=node(kn)%ctype
                          node(sn)%cvalue_r=node(kn)%cvalue_r
                        end do

                    end select

              ! ==============
              ! BE-FE BOUNDARY
              ! ==============

              ! Coupled boundaries always have unknown pressures. The displacements of a node can be constrained, if so, the
              ! corresponding finite element node must have the same constraint.
              ! The same for ordinary boundaries or crack-like boundaries.
              case (fbem_boundary_coupling_be_fe)
                select case (boundary(sb)%class)

                  ! ================= !
                  ! ORDINARY BOUNDARY !
                  ! ================= !

                  case (fbem_boundary_class_ordinary)


                  ! =================== !
                  ! CRACK-LIKE BOUNDARY !
                  ! =================== !

                  case (fbem_boundary_class_cracklike)

                end select

              ! =================
              ! BE-FE-BE BOUNDARY
              ! =================

              ! Coupled boundaries always have unknown pressures. The displacements of a node can be constrained, if so, the
              ! corresponding finite element node must have the same constraint.
              case (fbem_boundary_coupling_be_fe_be)


            end select

          end if
        end if

      end if
    end do

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
          if (verbose_level.gt.1) write(*,*) keyword
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
                    ! Read the type of boundary condition for k direction
                    do k=1,problem%n
                      read(input_fileunit,*) node(sn)%ctype(k,1)
                      if (k.eq.1) then
                        call fbem_search_section(input_fileunit,section_name,found)
                        call fbem_search_keyword(input_fileunit,keyword,':',found)
                      else
                        backspace(input_fileunit)
                      end if
                      ! Switch depending on the type of boundary condition
                      select case (node(sn)%ctype(k,1))
                        ! Dirichlet (known displacement): value u_k
                        case (0)
                          read(input_fileunit,*) node(sn)%ctype(k,1), node(sn)%cvalue_r(k,1,1)
                        ! Neumann (known pressure): value p_k
                        case (1)
                          read(input_fileunit,*) node(sn)%ctype(k,1), node(sn)%cvalue_r(k,1,1)
                      end select
                    end do

                  ! -------------------
                  ! CRACK-LIKE BOUNDARY
                  ! -------------------

                  case (fbem_boundary_class_cracklike)
                    ! Face +
                    ! Read the type of boundary condition for k direction
                    do k=1,problem%n
                      read(input_fileunit,*) node(sn)%ctype(k,1)
                      if (k.eq.1) then
                        call fbem_search_section(input_fileunit,section_name,found)
                        call fbem_search_keyword(input_fileunit,keyword,':',found)
                      else
                        backspace(input_fileunit)
                      end if
                      ! Switch depending on the type of boundary condition
                      select case (node(sn)%ctype(k,1))
                        ! Dirichlet (known displacement): value u_k
                        case (0)
                          read(input_fileunit,*) node(sn)%ctype(k,1), node(sn)%cvalue_r(k,1,1)
                        ! Neumann (known pressure): value p_k
                        case (1)
                          read(input_fileunit,*) node(sn)%ctype(k,1), node(sn)%cvalue_r(k,1,1)
                      end select
                    end do
                    ! Face -
                    ! Read the type of boundary condition for k direction
                    do k=1,problem%n
                      read(input_fileunit,*) node(sn)%ctype(k,2)
                      backspace(input_fileunit)
                      ! Switch depending on the type of boundary condition
                      select case (node(sn)%ctype(k,2))
                        ! Dirichlet (known displacement): value u_k
                        case (0)
                          read(input_fileunit,*) node(sn)%ctype(k,2), node(sn)%cvalue_r(k,1,2)
                        ! Neumann (known pressure): value p_k
                        case (1)
                          read(input_fileunit,*) node(sn)%ctype(k,2), node(sn)%cvalue_r(k,1,2)
                      end select
                    end do

                end select

            ! ==============
            ! BE-FE BOUNDARY
            ! ==============
            ! Coupled boundaries always have unknown pressures. The displacements of a node can be constrained, if so, the
            ! corresponding finite element node must have the same constraint.
            ! The same for open boundaries or crack-like boundaries.

            case (fbem_boundary_coupling_be_fe)
              ! Read the type of boundary condition for k direction
              do k=1,problem%n
                read(input_fileunit,*) node(sn)%ctype(k,1)
                if (k.eq.1) then
                  call fbem_search_section(input_fileunit,section_name,found)
                  call fbem_search_keyword(input_fileunit,keyword,':',found)
                else
                  backspace(input_fileunit)
                end if
                ! If the displacement is constrained
                if (node(sn)%ctype(k,1).eq.0) then
                  read(input_fileunit,*) node(sn)%ctype(k,1), node(sn)%cvalue_r(k,1,1)
                end if
              end do

            ! =================
            ! BE-FE-BE BOUNDARY
            ! =================
            ! Coupled boundaries always have unknown pressures. The displacements of a node can be constrained, if so, the
            ! corresponding finite element node must have the same constraint.

            case (fbem_boundary_coupling_be_fe_be)
              ! Read the type of boundary condition for k direction
              do k=1,problem%n
                read(input_fileunit,*) node(sn)%ctype(k,1)
                if (k.eq.1) then
                  call fbem_search_section(input_fileunit,section_name,found)
                  call fbem_search_keyword(input_fileunit,keyword,':',found)
                else
                  backspace(input_fileunit)
                end if
                ! If the displacement is constrained
                if (node(sn)%ctype(k,1).eq.0) then
                  read(input_fileunit,*) node(sn)%ctype(k,1), node(sn)%cvalue_r(k,1,1)
                end if
              end do

          end select
        end if



      end do
    end do

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  end if


end subroutine read_conditions_bem_nodes_mechanics_static
