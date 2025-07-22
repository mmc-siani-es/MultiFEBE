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


subroutine read_conditions_bem_nodes_mechanics_harmonic(input_fileunit)

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
  integer                                 :: kb, sb, sp
  integer                                 :: kn, sn
  integer                                 :: kg, kng
  character(len=fbem_stdcharlen)          :: keyword
  character(len=fbem_stdcharlen)          :: section_name
  !real(kind=real64)                       :: c(3), a(3), theta, x(3), urot(3), tmp_theta

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

          ! --------------------------
          ! GROUP OF BE BOUNDARY NODES
          ! --------------------------

          ! Read only groups of nodes belonging to a BE boundary
          if (part(sp)%type.eq.fbem_part_be_boundary) then
            ! BE boundary of part
            sb=part(sp)%entity
            select case (boundary(sb)%coupling)

              ! ===========
              ! BE BOUNDARY
              ! ===========

              case (fbem_boundary_coupling_be)
                select case (region(boundary(sb)%region(1))%type)

                  ! --------------
                  ! INVISCID FLUID
                  ! --------------

                  case (fbem_potential)
                    select case (boundary(sb)%class)

                      ! -----------------
                      ! ORDINARY BOUNDARY
                      ! -----------------

                      case (fbem_boundary_class_ordinary)
                        ! Read the type of boundary condition
                        read(input_fileunit,*) node(kn)%ctype(1,1)
                        call fbem_search_section(input_fileunit,section_name,found)
                        call fbem_search_keyword(input_fileunit,keyword,':',found)
                        ! Switch depending on the type of boundary condition
                        select case (node(kn)%ctype(1,1))
                          ! Dirichlet (known pressure): value p
                          case (0)
                            read(input_fileunit,*) node(kn)%ctype(1,1), node(kn)%cvalue_c(1,1,1)
                          ! Neumann (known normal displacement): value Un
                          case (1)
                            read(input_fileunit,*) node(kn)%ctype(1,1), node(kn)%cvalue_c(1,1,1)
                        end select
                        ! Copy the condition to the rest of nodes
                        do kng=2,group(kg)%n_objects
                          sn=group(kg)%object(kng)
                          node(sn)%ctype=node(kn)%ctype
                          node(sn)%cvalue_c=node(kn)%cvalue_c
                        end do

                      ! -------------------
                      ! CRACK-LIKE BOUNDARY
                      ! -------------------

                      case (fbem_boundary_class_cracklike)
                        ! Face +
                        ! Read the type of boundary condition
                        read(input_fileunit,*) node(kn)%ctype(1,1)
                        call fbem_search_section(input_fileunit,section_name,found)
                        call fbem_search_keyword(input_fileunit,keyword,':',found)
                        ! Switch depending on the type of boundary condition
                        select case (node(kn)%ctype(1,1))
                          ! Dirichlet (known pressure): value p
                          case (0)
                            read(input_fileunit,*) node(kn)%ctype(1,1), node(kn)%cvalue_c(1,1,1)
                          ! Neumann (known normal displacement): value Un
                          case (1)
                            read(input_fileunit,*) node(kn)%ctype(1,1), node(kn)%cvalue_c(1,1,1)
                        end select
                        ! Face -
                        ! Read the type of boundary condition
                        read(input_fileunit,*) node(kn)%ctype(1,2)
                        backspace(input_fileunit)
                        ! Switch depending on the type of boundary condition
                        select case (node(kn)%ctype(1,2))
                          ! Dirichlet (known pressure): value p
                          case (0)
                            read(input_fileunit,*) node(kn)%ctype(1,2), node(kn)%cvalue_c(1,1,2)
                          ! Neumann (known normal displacement): value Un
                          case (1)
                            read(input_fileunit,*) node(kn)%ctype(1,2), node(kn)%cvalue_c(1,1,2)
                        end select
                        ! Copy the condition to the rest of nodes
                        do kng=2,group(kg)%n_objects
                          sn=group(kg)%object(kng)
                          node(sn)%ctype=node(kn)%ctype
                          node(sn)%cvalue_c=node(kn)%cvalue_c
                        end do

                    end select

                  ! ------------------
                  ! VISCOELASTIC SOLID
                  ! ------------------

                  case (fbem_viscoelastic)
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
                              read(input_fileunit,*) node(kn)%ctype(k,1), node(kn)%cvalue_c(k,1,1)
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
                          node(sn)%cvalue_c=node(kn)%cvalue_c
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
                              read(input_fileunit,*) node(kn)%ctype(k,1), node(kn)%cvalue_c(k,1,1)
                            ! Neumann (known pressure): value p_k
                            case (1)
                              read(input_fileunit,*) node(kn)%ctype(k,1), node(kn)%cvalue_c(k,1,1)
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
                              read(input_fileunit,*) node(kn)%ctype(k,2), node(kn)%cvalue_c(k,1,2)
                            ! Neumann (known pressure): value p_k
                            case (1)
                              read(input_fileunit,*) node(kn)%ctype(k,2), node(kn)%cvalue_c(k,1,2)
                          end select
                        end do
                        ! Copy the condition to the rest of nodes
                        do kng=2,group(kg)%n_objects
                          sn=group(kg)%object(kng)
                          node(sn)%ctype=node(kn)%ctype
                          node(sn)%cvalue_c=node(kn)%cvalue_c
                        end do

                    end select

                  ! -----------------
                  ! POROELASTIC MEDIA
                  ! -----------------

                  case (fbem_poroelastic)
                    stop 'not implemented yet'

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
                    select case (region(boundary(sb)%region(1))%type)

                      ! -------------- !
                      ! INVISCID FLUID !
                      ! -------------- !

                      case (fbem_potential)
                        ! It does not need B.C.

                      ! ------------------ !
                      ! VISCOELASTIC SOLID !
                      ! ------------------ !

                      case (fbem_viscoelastic)
                        ! It does not need B.C.

                      ! ----------------- !
                      ! POROELASTIC MEDIUM !
                      ! ----------------- !

                      case (fbem_poroelastic)
                        ! It does not need B.C.

                    end select

                  ! =================== !
                  ! CRACK-LIKE BOUNDARY !
                  ! =================== !

                  case (fbem_boundary_class_cracklike)
                    select case (region(boundary(sb)%region(1))%type)

                      ! -------------- !
                      ! INVISCID FLUID !
                      ! -------------- !

                      case (fbem_potential)
                        ! It does not need B.C.

                      ! ------------------ !
                      ! VISCOELASTIC SOLID !
                      ! ------------------ !

                      case (fbem_viscoelastic)
                        ! It does not need B.C.

                      ! ----------------- !
                      ! POROELASTIC MEDIUM !
                      ! ----------------- !

                      case (fbem_poroelastic)
                        ! It does not need B.C.

                    end select
                end select

              ! =================
              ! BE-FE-BE BOUNDARY
              ! =================

              ! Coupled boundaries always have unknown pressures. The displacements of a node can be constrained, if so, the
              ! corresponding finite element node must have the same constraint.
              case (fbem_boundary_coupling_be_fe_be)
                select case (region(boundary(sb)%region(1))%type)
                  case (fbem_potential)
                    select case (region(boundary(sb)%region(2))%type)

                      ! ------------------------------------------------- !
                      ! BE (inviscid fluid) - BE (inviscid fluid)         !
                      ! ------------------------------------------------- !

                      case (fbem_potential)
                        ! It does not need B.C.
                      ! ------------------------------------------------------- !
                      ! BE (1: inviscid fluid) - BE (2: viscoelastic solid)     !
                      ! ------------------------------------------------------- !

                      case (fbem_viscoelastic)
                        ! It does not need B.C.
                      ! -------------------------------------------------- !
                      ! BE (1: inviscid fluid) - BE (2: poroelastic medium) !
                      ! -------------------------------------------------- !

                      case (fbem_poroelastic)
                        ! It does not need B.C.
                    end select

                  case (fbem_viscoelastic)
                    select case (region(boundary(sb)%region(2))%type)

                      ! ---------------------------------------------------- !
                      ! BE (1: viscoelastic solid) - BE (2: inviscid fluid)  !
                      ! ---------------------------------------------------- !

                      case (fbem_potential)
                        ! It does not need B.C.
                      ! ------------------------------------------------- !
                      ! BE (viscoelastic solid) - BE (viscoelastic solid) !
                      ! ------------------------------------------------- !

                      case (fbem_viscoelastic)
                        ! It does not need B.C.

                      ! ------------------------------------------------------ !
                      ! BE (1: viscoelastic solid) - BE (2: poroelastic medium) !
                      ! ------------------------------------------------------ !

                      case (fbem_poroelastic)
                        ! It does not need B.C.
                    end select

                  case (fbem_poroelastic)
                    select case (region(boundary(sb)%region(2))%type)

                      ! -------------------------------------------------- !
                      ! BE (1: poroelastic medium) - BE (2: inviscid fluid) !
                      ! -------------------------------------------------- !

                      case (fbem_potential)
                        ! It does not need B.C.
                      ! ------------------------------------------------------ !
                      ! BE (1: poroelastic medium) - BE (2: viscoelastic solid) !
                      ! ------------------------------------------------------ !

                      case (fbem_viscoelastic)
                        ! It does not need B.C.
                      ! ----------------------------------------------------- !
                      ! BE (1: poroelastic medium) - BE (2: poroelastic medium) !
                      ! ----------------------------------------------------- !

                      case (fbem_poroelastic)
                        ! It does not need B.C.
                    end select

                end select

            end select

          end if

          ! ---------------------------
          ! GROUP OF BE BODY LOAD NODES
          ! ---------------------------

          if (part(sp)%type.eq.fbem_part_be_bodyload) then
            ! BE body load of part
            sb=part(sp)%entity

            ! to be implemented

          end if

        end if

      end if
    end do

    ! ==============================================================================================================================
    ! FIND BE BOUNDARY NODES
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
              select case (region(boundary(kb)%region(1))%type)

                ! --------------
                ! INVISCID FLUID
                ! --------------

                case (fbem_potential)
                  select case (boundary(kb)%class)

                    ! -----------------
                    ! ORDINARY BOUNDARY
                    ! -----------------

                    case (fbem_boundary_class_ordinary)
                      ! Read the type of boundary condition
                      read(input_fileunit,*) node(sn)%ctype(1,1)
                      call fbem_search_section(input_fileunit,section_name,found)
                      call fbem_search_keyword(input_fileunit,keyword,':',found)
                      ! Switch depending on the type of boundary condition
                      select case (node(sn)%ctype(1,1))
                        ! Dirichlet (known pressure): value p
                        case (0)
                          read(input_fileunit,*) node(sn)%ctype(1,1), node(sn)%cvalue_c(1,1,1)
                        ! Neumann (known normal displacement): value Un
                        case (1)
                          read(input_fileunit,*) node(sn)%ctype(1,1), node(sn)%cvalue_c(1,1,1)
                      end select

                    ! -------------------
                    ! CRACK-LIKE BOUNDARY
                    ! -------------------

                    case (fbem_boundary_class_cracklike)
                      ! Face +
                      ! Read the type of boundary condition
                      read(input_fileunit,*) node(sn)%ctype(1,1)
                      call fbem_search_section(input_fileunit,section_name,found)
                      call fbem_search_keyword(input_fileunit,keyword,':',found)
                      ! Switch depending on the type of boundary condition
                      select case (node(sn)%ctype(1,1))
                        ! Dirichlet (known pressure): value p
                        case (0)
                          read(input_fileunit,*) node(sn)%ctype(1,1), node(sn)%cvalue_c(1,1,1)
                        ! Neumann (known normal displacement): value Un
                        case (1)
                          read(input_fileunit,*) node(sn)%ctype(1,1), node(sn)%cvalue_c(1,1,1)
                      end select
                      ! Face -
                      ! Read the type of boundary condition
                      read(input_fileunit,*) node(sn)%ctype(1,2)
                      backspace(input_fileunit)
                      ! Switch depending on the type of boundary condition
                      select case (node(sn)%ctype(1,2))
                        ! Dirichlet (known pressure): value p
                        case (0)
                          read(input_fileunit,*) node(sn)%ctype(1,2), node(sn)%cvalue_c(1,1,2)
                        ! Neumann (known normal displacement): value Un
                        case (1)
                          read(input_fileunit,*) node(sn)%ctype(1,2), node(sn)%cvalue_c(1,1,2)
                      end select
                  end select

                ! ------------------
                ! VISCOELASTIC SOLID
                ! ------------------

                case (fbem_viscoelastic)
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
                          ! 0. Global axes                         : u_k = U
                          ! 1. Global axes                         : t_k = T
                          ! 2. Local axes (l_1=n, l_2=t_1, l_3=t_2): u·l = U
                          ! 3. Local axes (l_1=n, l_2=t_1, l_3=t_2): t·l = T
                          case (0,1,2,3)
                            read(input_fileunit,*) node(sn)%ctype(k,1), node(sn)%cvalue_c(k,1,1)
                        end select
                      end do
                      ! Check that the B.C. are for local or global axes, but not mixed.
                      if ((node(sn)%ctype(1,1).eq.0).or.(node(sn)%ctype(1,1).eq.1)) then
                        do k=2,problem%n
                          if ((node(sn)%ctype(k,1).ne.0).and.(node(sn)%ctype(k,1).ne.1)) then
                            call fbem_error_message(error_unit,0,'node',node(sn)%id,&
                                                    'local axes B.C. and global axes B.C. can not be mixed.')
                          end if
                        end do
                      end if
                      if ((node(sn)%ctype(1,1).eq.2).or.(node(sn)%ctype(1,1).eq.3)) then
                        do k=2,problem%n
                          if ((node(sn)%ctype(k,2).ne.2).and.(node(sn)%ctype(k,1).ne.3)) then
                            call fbem_error_message(error_unit,0,'node',node(sn)%id,&
                                                    'local axes B.C. and global axes B.C. can not be mixed.')
                          end if
                        end do
                      end if

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
                            read(input_fileunit,*) node(sn)%ctype(k,1), node(sn)%cvalue_c(k,1,1)
                          ! Neumann (known pressure): value p_k
                          case (1)
                            read(input_fileunit,*) node(sn)%ctype(k,1), node(sn)%cvalue_c(k,1,1)
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
                            read(input_fileunit,*) node(sn)%ctype(k,2), node(sn)%cvalue_c(k,1,2)
                          ! Neumann (known pressure): value p_k
                          case (1)
                            read(input_fileunit,*) node(sn)%ctype(k,2), node(sn)%cvalue_c(k,1,2)
                        end select
                      end do
                  end select

                ! -----------------
                ! POROELASTIC MEDIA
                ! -----------------

                case (fbem_poroelastic)
                  stop 'not implemented yet'

              end select

            ! ==============
            ! BE-FE BOUNDARY
            ! ==============

            ! Coupled boundaries always have unknown pressures. The displacements of a node can be constrained, if so, the
            ! corresponding finite element node must have the same constraint.
            ! The same for ordinary boundaries or crack-like boundaries.
            case (fbem_boundary_coupling_be_fe)
              select case (boundary(kb)%class)

                ! ================= !
                ! ORDINARY BOUNDARY !
                ! ================= !

                case (fbem_boundary_class_ordinary)
                  select case (region(boundary(kb)%region(1))%type)

                    ! -------------- !
                    ! INVISCID FLUID !
                    ! -------------- !

                    case (fbem_potential)
                      ! It does not need B.C.

                    ! ------------------ !
                    ! VISCOELASTIC SOLID !
                    ! ------------------ !

                    case (fbem_viscoelastic)
                      ! It does not need B.C.

                    ! ----------------- !
                    ! POROELASTIC MEDIUM !
                    ! ----------------- !

                    case (fbem_poroelastic)
                      ! It does not need B.C.

                  end select

                ! =================== !
                ! CRACK-LIKE BOUNDARY !
                ! =================== !

                case (fbem_boundary_class_cracklike)
                  select case (region(boundary(kb)%region(1))%type)

                    ! -------------- !
                    ! INVISCID FLUID !
                    ! -------------- !

                    case (fbem_potential)
                      ! It does not need B.C.

                    ! ------------------ !
                    ! VISCOELASTIC SOLID !
                    ! ------------------ !

                    case (fbem_viscoelastic)
                      ! It does not need B.C.

                    ! ----------------- !
                    ! POROELASTIC MEDIUM !
                    ! ----------------- !

                    case (fbem_poroelastic)
                      ! It does not need B.C.

                  end select
              end select

            ! =================
            ! BE-FE-BE BOUNDARY
            ! =================

            ! Coupled boundaries always have unknown pressures. The displacements of a node can be constrained, if so, the
            ! corresponding finite element node must have the same constraint.
            case (fbem_boundary_coupling_be_fe_be)
              select case (region(boundary(kb)%region(1))%type)
                case (fbem_potential)
                  select case (region(boundary(kb)%region(2))%type)

                    ! ------------------------------------------------- !
                    ! BE (inviscid fluid) - BE (inviscid fluid)         !
                    ! ------------------------------------------------- !

                    case (fbem_potential)
                      ! It does not need B.C.
                    ! ------------------------------------------------------- !
                    ! BE (1: inviscid fluid) - BE (2: viscoelastic solid)     !
                    ! ------------------------------------------------------- !

                    case (fbem_viscoelastic)
                      ! It does not need B.C.
                    ! -------------------------------------------------- !
                    ! BE (1: inviscid fluid) - BE (2: poroelastic medium) !
                    ! -------------------------------------------------- !

                    case (fbem_poroelastic)
                      ! It does not need B.C.
                  end select

                case (fbem_viscoelastic)
                  select case (region(boundary(kb)%region(2))%type)

                    ! ---------------------------------------------------- !
                    ! BE (1: viscoelastic solid) - BE (2: inviscid fluid)  !
                    ! ---------------------------------------------------- !

                    case (fbem_potential)
                      ! It does not need B.C.
                    ! ------------------------------------------------- !
                    ! BE (viscoelastic solid) - BE (viscoelastic solid) !
                    ! ------------------------------------------------- !

                    case (fbem_viscoelastic)
                      ! It does not need B.C.

                    ! ------------------------------------------------------ !
                    ! BE (1: viscoelastic solid) - BE (2: poroelastic medium) !
                    ! ------------------------------------------------------ !

                    case (fbem_poroelastic)
                      ! It does not need B.C.
                  end select

                case (fbem_poroelastic)
                  select case (region(boundary(kb)%region(2))%type)

                    ! -------------------------------------------------- !
                    ! BE (1: poroelastic medium) - BE (2: inviscid fluid) !
                    ! -------------------------------------------------- !

                    case (fbem_potential)
                      ! It does not need B.C.
                    ! ------------------------------------------------------ !
                    ! BE (1: poroelastic medium) - BE (2: viscoelastic solid) !
                    ! ------------------------------------------------------ !

                    case (fbem_viscoelastic)
                      ! It does not need B.C.
                    ! ----------------------------------------------------- !
                    ! BE (1: poroelastic medium) - BE (2: poroelastic medium) !
                    ! ----------------------------------------------------- !

                    case (fbem_poroelastic)
                      ! It does not need B.C.
                  end select

              end select

          end select
        end if

      end do
    end do

    ! ==============================================================================================================================
    ! FIND BE BODY LOAD NODES
    ! ==============================================================================================================================

    do kb=1,n_be_bodyloads
      sp=be_bodyload(kb)%part
      do kn=1,part(sp)%n_nodes
        sn=part(sp)%node(kn)
        ! Locate the node
        call fbem_search_section(input_fileunit,section_name,found)
        write(keyword,*) 'node ', node(sn)%id
        call fbem_trim2b(keyword)
        call fbem_search_keyword(input_fileunit,keyword,':',found)
        if (found) then
          if (verbose_level.gt.1) write(*,*) keyword
          select case (be_bodyload(kb)%coupling)
            !
            ! UNCOUPLED BE BODY LOAD
            !
            case (fbem_bl_uncoupled)

              select case (region(be_bodyload(kb)%region)%type)

                ! --------------
                ! INVISCID FLUID
                ! --------------

                case (fbem_potential)
                  read(input_fileunit,*) node(sn)%ctype(1,1)
                  call fbem_search_section(input_fileunit,section_name,found)
                  call fbem_search_keyword(input_fileunit,keyword,':',found)
                  ! Switch depending on the type of boundary condition
                  select case (node(sn)%ctype(1,1))
                    ! Constant amplitude
                    case (0)
                      read(input_fileunit,*) node(sn)%ctype(1,1), node(sn)%cvalue_c(1,1,1)
                    case default
                      call fbem_error_message(error_unit,0,trim(section_name),node(sn)%id,'invalid condition for this be body load node.')
                  end select

                ! ------------------
                ! VISCOELASTIC SOLID
                ! ------------------

                case (fbem_viscoelastic)
                  ! Read the type of condition for k direction
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
                      ! Force value
                      case (0)
                        read(input_fileunit,*) node(sn)%ctype(k,1), node(sn)%cvalue_c(k,1,1)
                      case default
                        call fbem_error_message(error_unit,0,trim(section_name),node(sn)%id,'invalid condition for this be body load node.')
                    end select
                  end do

                ! -----------------
                ! POROELASTIC MEDIUM
                ! -----------------

                case (fbem_poroelastic)
                  call fbem_error_message(error_unit,0,__FILE__,__LINE__,'body loads not available for poroelastic media')

              end select

            case (fbem_bl_coupling_beam_tip)
              ! N/A

            case (fbem_bl_coupling_beam_line)
              ! N/A

            case (fbem_bl_coupling_shell_edge)
              ! N/A

            case (fbem_bl_coupling_shell_surface)
              ! N/A

          end select

        end if
      end do
    end do

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  end if


end subroutine read_conditions_bem_nodes_mechanics_harmonic
