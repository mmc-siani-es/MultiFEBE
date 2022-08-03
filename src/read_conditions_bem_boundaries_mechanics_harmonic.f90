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


subroutine read_conditions_bem_boundaries_mechanics_harmonic(input_fileunit)

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
  integer                                 :: input_fileunit    !! Input file unit
  ! Local variables
  character(len=fbem_stdcharlen)          :: section_name      ! Name of the section
  logical                                 :: found
  integer                                 :: i, k
  character(len=fbem_stdcharlen)          :: keyword
  real(kind=real64)                       :: center(problem%n), axis(3)
  complex(kind=real64)                    :: theta

  ! Locate the section
  section_name='conditions over be boundaries'
if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(input_fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! Loop through BOUNDARIES
    do i=1,n_boundaries
      ! Find the boundary
      call fbem_search_section(input_fileunit,section_name,found)
      write(keyword,*) 'boundary ', boundary(i)%id
      call fbem_trim2b(keyword)
      call fbem_search_keyword(input_fileunit,keyword,':',found)
      if (found) then
        select case (boundary(i)%coupling)

          ! ========================================================================================================================
          ! BE BOUNDARY
          ! ========================================================================================================================

          case (fbem_boundary_coupling_be)
            select case (boundary(i)%class)

              ! =================
              ! ORDINARY BOUNDARY
              ! =================

              case (fbem_boundary_class_ordinary)
                select case (region(boundary(i)%region(1))%type)

                  ! --------------
                  ! INVISCID FLUID
                  ! --------------

                  case (fbem_potential)
                    ! Read the type of boundary condition
                    read(input_fileunit,*) boundary(i)%ctype(1,1)
                    call fbem_search_section(input_fileunit,section_name,found)
                    call fbem_search_keyword(input_fileunit,keyword,':',found)
                    ! Switch depending on the type of boundary condition
                    select case (boundary(i)%ctype(1,1))
                      ! 0: p=P
                      ! 1: Un=W
                      case (0,1)
                        read(input_fileunit,*) boundary(i)%ctype(1,1), boundary(i)%cvalue_c(1,1,1)
                      ! 2: Un+i/(rho*c*omega)p=0
                      case (2)
                        read(input_fileunit,*) boundary(i)%ctype(1,1)
                        boundary(i)%cvalue_c(1,1,1)=0.
                      ! 3: Un+(i/(rho*c*omega)-1/(2*R*rho*omega^2))p=0
                      case (3)
                        read(input_fileunit,*) boundary(i)%ctype(1,1), boundary(i)%cvalue_c(1,1,1)
                      case default
                        call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                    end select

                  ! ------------------
                  ! VISCOELASTIC SOLID
                  ! ------------------

                  case (fbem_viscoelastic)
                    ! Read the type of boundary condition for k direction
                    do k=1,problem%n
                      read(input_fileunit,*) boundary(i)%ctype(k,1)
                      if (k.eq.1) then
                        call fbem_search_section(input_fileunit,section_name,found)
                        call fbem_search_keyword(input_fileunit,keyword,':',found)
                      else
                        backspace(input_fileunit)
                      end if
                      ! Switch depending on the type of boundary condition
                      select case (boundary(i)%ctype(k,1))
                        ! 0: Global axes                         : u_k = U
                        ! 1: Global axes                         : t_k = T
                        ! 2: Local axes (l_1=n, l_2=t_1, l_3=t_2): u·l = U
                        ! 3: Local axes (l_1=n, l_2=t_1, l_3=t_2): t·l = T
                        case (0,1,2,3)
                          read(input_fileunit,*) boundary(i)%ctype(k,1), boundary(i)%cvalue_c(k,1,1)
                        ! 4: Global axes (inf. rotation field)   : u_k = theta (axis x (node_position-center))·e_k
                        case (4)
                          select case (problem%n)
                            case (2)
                              read(input_fileunit,*) boundary(i)%ctype(k,1), center, theta
                              axis=[0.d0,0.d0,1.d0]
                            case (3)
                              read(input_fileunit,*) boundary(i)%ctype(k,1), center, axis, theta
                          end select
                          boundary(i)%cvalue_c(k,:,1)=0.d0
                          boundary(i)%cvalue_c(k,1:problem%n,1)=center
                          boundary(i)%cvalue_c(k,4:6,1)=axis/sqrt(dot_product(axis,axis))
                          boundary(i)%cvalue_c(k,7,1)=theta
                        ! 10: normal pressure: t_k = P n_k
                        case (10)
                          read(input_fileunit,*) boundary(i)%ctype(1,1), boundary(i)%cvalue_c(1,1,1)
                          boundary(i)%ctype(:,1)=boundary(i)%ctype(1,1)
                          boundary(i)%cvalue_c(:,1,1)=boundary(i)%cvalue_c(1,1,1)
                          exit
                        case default
                          call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                      end select
                    end do
                    ! Check that the B.C. are for local or global axes, but not mixed.
                    if ((boundary(i)%ctype(1,1).eq.0).or.(boundary(i)%ctype(1,1).eq.1)) then
                      do k=2,problem%n
                        if ((boundary(i)%ctype(k,1).ne.0).and.(boundary(i)%ctype(k,1).ne.1)) then
                          call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'local axes B.C. and global axes B.C. can not be mixed.')
                        end if
                      end do
                    end if
                    if ((boundary(i)%ctype(1,1).eq.2).or.(boundary(i)%ctype(1,1).eq.3)) then
                      do k=2,problem%n
                        if ((boundary(i)%ctype(k,1).ne.2).and.(boundary(i)%ctype(k,1).ne.3)) then
                          call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'local axes B.C. and global axes B.C. can not be mixed.')
                        end if
                      end do
                    end if

                  ! -----------------
                  ! POROELASTIC MEDIUM
                  ! -----------------

                  case (fbem_poroelastic)
                    ! Read fluid phase B.C.
                    read(input_fileunit,*) boundary(i)%ctype(0,1)
                    call fbem_search_section(input_fileunit,section_name,found)
                    call fbem_search_keyword(input_fileunit,keyword,':',found)
                    ! Switch depending on the type of boundary condition
                    select case (boundary(i)%ctype(0,1))
                      ! open pore / permeable
                      ! 0: equivalent stress  : tau=T
                      ! 1: normal displacement: Un=U
                      case (0,1)
                        ! Fluid phase
                        read(input_fileunit,*) boundary(i)%ctype(0,1), boundary(i)%cvalue_c(0,1,1)
                        ! Solid skeleton
                        ! Read the type of boundary condition for k direction
                        do k=1,problem%n
                          read(input_fileunit,*) boundary(i)%ctype(k,1)
                          backspace(input_fileunit)
                          ! Switch depending on the type of boundary condition
                          select case (boundary(i)%ctype(k,1))
                            ! Global axes
                            ! 0: Displacement : u_k = U
                            ! 1: Traction     : t_k = T
                            ! Local axes
                            ! 2: Displacement : u·l = U
                            ! 3: Traction     : t·l = T
                            case (0,1,2,3)
                              read(input_fileunit,*) boundary(i)%ctype(k,1), boundary(i)%cvalue_c(k,1,1)
                            ! 8: Global axes (inf. rotation field)   : u_k = theta (axis x (node_position-center))·e_k
                            case (8)
                              select case (problem%n)
                                ! 2D
                                case (2)
                                  read(input_fileunit,*) boundary(i)%ctype(k,1), center, theta
                                  axis=[0.d0,0.d0,1.d0]
                                ! 3D
                                case (3)
                                  read(input_fileunit,*) boundary(i)%ctype(k,1), center, axis, theta
                              end select
                              boundary(i)%cvalue_c(k,:,1)=0.d0
                              boundary(i)%cvalue_c(k,1:problem%n,1)=center
                              boundary(i)%cvalue_c(k,4:6,1)=axis/sqrt(dot_product(axis,axis))
                              boundary(i)%cvalue_c(k,7,1)=theta
                            ! 10: normal pressure: t_k = P n_k
                            case (10)
                              read(input_fileunit,*) boundary(i)%ctype(1,1), boundary(i)%cvalue_c(1,1,1)
                              boundary(i)%ctype(1:problem%n,1)=boundary(i)%ctype(1,1)
                              boundary(i)%cvalue_c(1:problem%n,1,1)=boundary(i)%cvalue_c(1,1,1)
                              exit
                            case default
                              call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                          end select
                        end do
                        ! Check that the B.C. are for local or global axes, but not mixed.
                        if ((boundary(i)%ctype(1,1).eq.0).or.(boundary(i)%ctype(1,1).eq.1)) then
                          do k=2,problem%n
                            if ((boundary(i)%ctype(k,1).ne.0).and.(boundary(i)%ctype(k,1).ne.1)) then
                              call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'local axes B.C. and global axes B.C. can not be mixed.')
                            end if
                          end do
                        end if
                        if ((boundary(i)%ctype(1,1).eq.2).or.(boundary(i)%ctype(1,1).eq.3)) then
                          do k=2,problem%n
                            if ((boundary(i)%ctype(k,1).ne.2).and.(boundary(i)%ctype(k,1).ne.3)) then
                              call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'local axes B.C. and global axes B.C. can not be mixed.')
                            end if
                          end do
                        end if
                      ! close pore / impermeable
                      ! 2: impervious condition: Un=u_k·n_k
                      case (2)
                        read(input_fileunit,*) boundary(i)%ctype(0,1)
                        ! Solid skeleton
                        ! Read the type of boundary condition for k direction
                        do k=1,problem%n
                          read(input_fileunit,*) boundary(i)%ctype(k,1)
                          backspace(input_fileunit)
                          ! Switch depending on the type of boundary condition
                          select case (boundary(i)%ctype(k,1))
                            ! Global axes
                            ! 4: Displacement : u_k = U
                            ! 5: Traction     : t_k + tau n_k = T
                            ! Local axes
                            ! 6: Displacement : u·l = U
                            ! 7: Traction     : t·l + tau n·l = T
                            case (4,5,6,7)
                              read(input_fileunit,*) boundary(i)%ctype(k,1), boundary(i)%cvalue_c(k,1,1)
                            ! 40: Global axes (inf. rotation field)   : u_k = theta (axis x (node_position-center))·e_k
                            case (40)
                              select case (problem%n)
                                ! 2D
                                case (2)
                                  read(input_fileunit,*) boundary(i)%ctype(k,1), center, theta
                                  axis=[0.d0,0.d0,1.d0]
                                ! 3D
                                case (3)
                                  read(input_fileunit,*) boundary(i)%ctype(k,1), center, axis, theta
                              end select
                              boundary(i)%cvalue_c(k,:,1)=0.d0
                              boundary(i)%cvalue_c(k,1:problem%n,1)=center
                              boundary(i)%cvalue_c(k,4:6,1)=axis/sqrt(dot_product(axis,axis))
                              boundary(i)%cvalue_c(k,7,1)=theta
                            ! 50: total normal pressure: t_k + tau n_k = P n_k
                            case (50)
                              read(input_fileunit,*) boundary(i)%ctype(1,1), boundary(i)%cvalue_c(1,1,1)
                              boundary(i)%ctype(1:problem%n,1)=boundary(i)%ctype(1,1)
                              boundary(i)%cvalue_c(1:problem%n,1,1)=boundary(i)%cvalue_c(1,1,1)
                              exit
                            case default
                              call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                          end select
                        end do
                        ! Check that the B.C. are for local or global axes, but not mixed.
                        if ((boundary(i)%ctype(1,1).eq.4).or.(boundary(i)%ctype(1,1).eq.5)) then
                          do k=2,problem%n
                            if ((boundary(i)%ctype(k,1).ne.4).and.(boundary(i)%ctype(k,1).ne.5)) then
                              call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'local axes B.C. and global axes B.C. can not be mixed.')
                            end if
                          end do
                        end if
                        if ((boundary(i)%ctype(1,1).eq.6).or.(boundary(i)%ctype(1,1).eq.7)) then
                          do k=2,problem%n
                            if ((boundary(i)%ctype(k,1).ne.6).and.(boundary(i)%ctype(k,1).ne.7)) then
                              call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'local axes B.C. and global axes B.C. can not be mixed.')
                            end if
                          end do
                        end if
                      case default
                        call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                    end select

                end select

              ! ===================
              ! CRACK-LIKE BOUNDARY
              ! ===================

              case (fbem_boundary_class_cracklike)
                select case (region(boundary(i)%region(1))%type)

                  ! --------------
                  ! INVISCID FLUID
                  ! --------------

                  case (fbem_potential)
                    ! Face +
                    read(input_fileunit,*) boundary(i)%ctype(1,1)
                    call fbem_search_section(input_fileunit,section_name,found)
                    call fbem_search_keyword(input_fileunit,keyword,':',found)
                    ! Switch depending on the type of boundary condition
                    select case (boundary(i)%ctype(1,1))
                      ! 0: p^+=P^+
                      ! 1: Un^+=W^+
                      ! 20: Un^+ - Un^- = K
                      ! 21: p^+ - p^- = K
                      case (0,1,20,21)
                        read(input_fileunit,*) boundary(i)%ctype(1,1), boundary(i)%cvalue_c(1,1,1)
                      case default
                        call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                    end select
                    ! Face -
                    read(input_fileunit,*) boundary(i)%ctype(1,2)
                    backspace(input_fileunit)
                    ! Switch depending on the type of boundary condition
                    select case (boundary(i)%ctype(1,2))
                      ! 0: p^-=P^-
                      ! 1: Un^-=W^-
                      ! 20: Un^+ - Un^- = K
                      ! 21: p^+ - p^- = K
                      case (0,1,20,21)
                        read(input_fileunit,*) boundary(i)%ctype(1,2), boundary(i)%cvalue_c(1,1,2)
                      case default
                        call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                    end select
                    if ((boundary(i)%ctype(1,1).eq.20).and.(boundary(i)%ctype(1,2).eq.20)) then
                      call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'incompatible combination of boundary conditions')
                    end if
                    if ((boundary(i)%ctype(1,1).eq.21).and.(boundary(i)%ctype(1,2).eq.21)) then
                      call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'incompatible combination of boundary conditions')
                    end if

                  ! ------------------
                  ! VISCOELASTIC SOLID
                  ! ------------------

                  case (fbem_viscoelastic)
                    ! Face +
                    ! Read the type of boundary condition for k direction
                    do k=1,problem%n
                      read(input_fileunit,*) boundary(i)%ctype(k,1)
                      if (k.eq.1) then
                        call fbem_search_section(input_fileunit,section_name,found)
                        call fbem_search_keyword(input_fileunit,keyword,':',found)
                      else
                        backspace(input_fileunit)
                      end if
                      ! Switch depending on the type of boundary condition
                      select case (boundary(i)%ctype(k,1))
                        ! 0: u_k^+=U^+
                        ! 1: t_k^+=T^+
                        case (0,1)
                          read(input_fileunit,*) boundary(i)%ctype(k,1), boundary(i)%cvalue_c(k,1,1)
                        case default
                          call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                      end select
                    end do
                    ! Face -
                    ! Read the type of boundary condition for k direction
                    do k=1,problem%n
                      read(input_fileunit,*) boundary(i)%ctype(k,2)
                      backspace(input_fileunit)
                      ! Switch depending on the type of boundary condition
                      select case (boundary(i)%ctype(k,2))
                        ! 0: u_k^-=U^-
                        ! 1: t_k^-=T^-
                        case (0,1)
                          read(input_fileunit,*) boundary(i)%ctype(k,2), boundary(i)%cvalue_c(k,1,2)
                        case default
                          call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                      end select
                    end do

                  ! -----------------
                  ! POROELASTIC MEDIA
                  ! -----------------

                  case (fbem_poroelastic)
                    ! Face +
                    ! Read fluid phase B.C.
                    read(input_fileunit,*) boundary(i)%ctype(0,1)
                    call fbem_search_section(input_fileunit,section_name,found)
                    call fbem_search_keyword(input_fileunit,keyword,':',found)
                    ! Switch depending on the type of boundary condition
                    select case (boundary(i)%ctype(0,1))
                      ! open pore / permeable
                      ! 0: equivalent stress  : tau=T
                      ! 1: normal displacement: Un=U
                      case (0,1)
                        ! Fluid phase
                        read(input_fileunit,*) boundary(i)%ctype(0,1), boundary(i)%cvalue_c(0,1,1)
                        ! Solid skeleton
                        ! Read the type of boundary condition for k direction
                        do k=1,problem%n
                          read(input_fileunit,*) boundary(i)%ctype(k,1)
                          backspace(input_fileunit)
                          ! Switch depending on the type of boundary condition
                          select case (boundary(i)%ctype(k,1))
                            ! 0: Displacement : u_k = U
                            ! 1: Traction     : t_k = T
                            case (0,1)
                              read(input_fileunit,*) boundary(i)%ctype(k,1), boundary(i)%cvalue_c(k,1,1)
                            case default
                              call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                          end select
                        end do
                      ! close pore / impermeable
                      ! 2: impervious condition: Un=u_k·n_k
                      case (2)
                        read(input_fileunit,*) boundary(i)%ctype(0,1)
                        ! Solid skeleton
                        ! Read the type of boundary condition for k direction
                        do k=1,problem%n
                          read(input_fileunit,*) boundary(i)%ctype(k,1)
                          backspace(input_fileunit)
                          ! Switch depending on the type of boundary condition
                          select case (boundary(i)%ctype(k,1))
                            ! Global axes
                            ! 4. Displacement : u_k = U
                            ! 5. Traction     : t_k = T
                            case (2,3)
                              read(input_fileunit,*) boundary(i)%ctype(k,1), boundary(i)%cvalue_c(k,1,1)
                            case default
                              call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                          end select
                        end do
                      ! 4: Fluid(inviscid/incompressible)-filled permeable crack
                      case (4)
                        boundary(i)%ctype(0,2)=4
                        do k=1,problem%n
                          boundary(i)%ctype(k,1)=4
                          boundary(i)%ctype(k,2)=4
                        end do
                      ! 5: Fluid(inviscid/incompressible)-filled impermeable crack
                      case (5)
                        boundary(i)%ctype(0,2)=5
                        do k=1,problem%n
                          boundary(i)%ctype(k,1)=5
                          boundary(i)%ctype(k,2)=5
                        end do
                      case default
                        call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                    end select
                    ! Face -
                    if (boundary(i)%ctype(0,1).lt.4) then
                      ! Read fluid phase B.C.
                      read(input_fileunit,*) boundary(i)%ctype(0,2)
                      backspace(input_fileunit)
                      ! Switch depending on the type of boundary condition
                      select case (boundary(i)%ctype(0,2))
                        ! open pore / permeable
                        ! 0: equivalent stress  : tau=T
                        ! 1: normal displacement: Un=U
                        case (0,1)
                          ! Fluid phase
                          read(input_fileunit,*) boundary(i)%ctype(0,2), boundary(i)%cvalue_c(0,1,2)
                          ! Solid skeleton
                          ! Read the type of boundary condition for k direction
                          do k=1,problem%n
                            read(input_fileunit,*) boundary(i)%ctype(k,2)
                            backspace(input_fileunit)
                            ! Switch depending on the type of boundary condition
                            select case (boundary(i)%ctype(k,2))
                              ! 0: Displacement : u_k = U
                              ! 1: Traction     : t_k = T
                              case (0,1)
                                read(input_fileunit,*) boundary(i)%ctype(k,2), boundary(i)%cvalue_c(k,1,2)
                              case default
                                call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                            end select
                          end do
                        ! close pore / impermeable
                        ! 2. impervious condition: Un=u_k·n_k
                        case (2)
                          read(input_fileunit,*) boundary(i)%ctype(0,2)
                          ! Solid skeleton
                          ! Read the type of boundary condition for k direction
                          do k=1,problem%n
                            read(input_fileunit,*) boundary(i)%ctype(k,2)
                            backspace(input_fileunit)
                            ! Switch depending on the type of boundary condition
                            select case (boundary(i)%ctype(k,2))
                              ! Global axes
                              ! 4: Displacement : u_k = U
                              ! 5: Traction     : t_k = T
                              case (2,3)
                                read(input_fileunit,*) boundary(i)%ctype(k,2), boundary(i)%cvalue_c(k,1,2)
                              case default
                                call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                            end select
                          end do
                        case default
                          call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                      end select
                    end if

                end select

            end select

          ! ========================================================================================================================
          ! BE-BE BOUNDARY
          ! ========================================================================================================================

          case (fbem_boundary_coupling_be_be)
            select case (region(boundary(i)%region(1))%type)
              case (fbem_potential)
                select case (region(boundary(i)%region(2))%type)

                  ! ---------------------------------------
                  ! INVISCID FLUID (1) - INVISCID FLUID (2)
                  ! ---------------------------------------

                  case (fbem_potential)
                    ! It does not need B.C.

                  ! -------------------------------------------
                  ! INVISCID FLUID (1) - VISCOELASTIC SOLID (2)
                  ! -------------------------------------------

                  case (fbem_viscoelastic)
                    ! It does not need B.C.

                  ! ------------------------------------------
                  ! INVISCID FLUID (1) - POROELASTIC MEDIA (2)
                  ! ------------------------------------------

                  case (fbem_poroelastic)
                    read(input_fileunit,*) boundary(i)%ctype(1,1)
                    select case (boundary(i)%ctype(1,1))
                      ! 0: perfectly permeable
                      ! 1: perfectly impermeable
                      case (0,1)
                      case default
                        call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                    end select

                end select

              case (fbem_viscoelastic)
                select case (region(boundary(i)%region(2))%type)

                  ! -------------------------------------------
                  ! VISCOELASTIC SOLID (1) - INVISCID FLUID (2)
                  ! -------------------------------------------

                  case (fbem_potential)
                    ! It does not need B.C.

                  ! -----------------------------------------------
                  ! VISCOELASTIC SOLID (1) - VISCOELASTIC SOLID (2)
                  ! -----------------------------------------------

                  ! ------------------------------------------------- !
                  ! Global axes                                       !
                  ! 0 (def.) - perfect bonding   - no values needed   !
                  ! 1        - perfect debonding - no values needed   !
                  ! 2        - partial bonding   - K, C, M            !
                  ! Local axes                                        !
                  ! 3        - perfect bonding   - no values needed   !
                  ! 4        - perfect debonding - no values needed   !
                  ! 5        - partial bonding   - K, C, M            !
                  ! ------------------------------------------------- !

                  ! ----------------------------------------------
                  ! VISCOELASTIC SOLID (1) - POROELASTIC MEDIA (2)
                  ! ----------------------------------------------

                  ! ------------------------------------------------------ !
                  ! Global axes                                            !
                  ! 0 (def.) - perfect bonding   - no values needed        !
                  ! 1        - perfect debonding - no values needed        !
                  ! 2        - partial bonding   - K, C, M                 !
                  ! Local axes                                             !
                  ! 3        - perfect bonding   - no values needed        !
                  ! 4        - perfect debonding - no values needed        !
                  ! 5        - partial bonding   - K, C, M                 !
                  ! ------------------------------------------------------ !
                  ! Note: the viscoelastic solid is considered impervious  !
                  ! ------------------------------------------------------ !

                  case (fbem_viscoelastic,fbem_poroelastic)
                    ! Read the type of boundary condition for k direction
                    do k=1,problem%n
                      read(input_fileunit,*) boundary(i)%ctype(k,1)
                      if (k.eq.1) then
                        call fbem_search_section(input_fileunit,section_name,found)
                        call fbem_search_keyword(input_fileunit,keyword,':',found)
                      else
                        backspace(input_fileunit)
                      end if
                      ! Switch depending on the type of boundary condition
                      select case (boundary(i)%ctype(k,1))
                        ! Global axes
                        ! 0: perfect bonding
                        ! 1: perfect debonding
                        ! Local axes
                        ! 3: perfect bonding
                        ! 4: perfect debonding
                        case (0,1,3,4)
                          read(input_fileunit,*) boundary(i)%ctype(k,1)
                        ! Global axes
                        ! 2: partial bonding
                        ! Local axes
                        ! 5: partial bonding
                        case (2,5)
                          read(input_fileunit,*) boundary(i)%ctype(k,1), boundary(i)%cvalue_c(k,1,1), boundary(i)%cvalue_c(k,2,1), boundary(i)%cvalue_c(k,3,1)
                        case default
                          call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                      end select
                    end do
                    ! Check that the B.C. are for local or global axes, but not mixed.
                    if ((boundary(i)%ctype(1,1).eq.0).or.(boundary(i)%ctype(1,1).eq.1).or.(boundary(i)%ctype(1,1).eq.2)) then
                      do k=2,problem%n
                        if ((boundary(i)%ctype(k,1).ne.0).and.(boundary(i)%ctype(k,1).ne.1).and.(boundary(i)%ctype(k,1).ne.2)) then
                          call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'local axes B.C. and global axes B.C. can not be mixed.')
                        end if
                      end do
                    end if
                    if ((boundary(i)%ctype(1,1).eq.3).or.(boundary(i)%ctype(1,1).eq.4).or.(boundary(i)%ctype(1,1).eq.5)) then
                      do k=2,problem%n
                        if ((boundary(i)%ctype(k,1).ne.3).and.(boundary(i)%ctype(k,1).ne.4).and.(boundary(i)%ctype(k,1).ne.5)) then
                          call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'local axes B.C. and global axes B.C. can not be mixed.')
                        end if
                      end do
                    end if

                end select

              case (fbem_poroelastic)
                select case (region(boundary(i)%region(2))%type)

                  ! ------------------------------------------
                  ! POROELASTIC MEDIA (1) - INVISCID FLUID (2)
                  ! ------------------------------------------

                  case (fbem_potential)
                    read(input_fileunit,*) boundary(i)%ctype(1,1)
                    select case (boundary(i)%ctype(1,1))
                      ! 0. perfectly permeable
                      ! 1. perfectly impermeable
                      case (0,1)
                      case default
                        call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                    end select

                  ! ----------------------------------------------
                  ! POROELASTIC MEDIA (1) - VISCOELASTIC SOLID (2)
                  ! ----------------------------------------------

                  ! ------------------------------------------------------ !
                  ! Global axes                                            !
                  ! 0 (def.) - perfect bonding   - no values needed        !
                  ! 1        - perfect debonding - no values needed        !
                  ! 2        - partial bonding   - K, C, M                 !
                  ! Local axes                                             !
                  ! 3        - perfect bonding   - no values needed        !
                  ! 4        - perfect debonding - no values needed        !
                  ! 5        - partial bonding   - K, C, M                 !
                  ! ------------------------------------------------------ !
                  ! Note: the viscoelastic solid is considered impervios   !
                  ! ------------------------------------------------------ !

                  case (fbem_viscoelastic)
                    ! Read the type of boundary condition for k direction
                    do k=1,problem%n
                      read(input_fileunit,*) boundary(i)%ctype(k,1)
                      if (k.eq.1) then
                        call fbem_search_section(input_fileunit,section_name,found)
                        call fbem_search_keyword(input_fileunit,keyword,':',found)
                      else
                        backspace(input_fileunit)
                      end if
                      ! Switch depending on the type of boundary condition
                      select case (boundary(i)%ctype(k,1))
                        ! Global axes
                        ! 0: perfect bonding
                        ! 1: perfect debonding
                        ! Local axes
                        ! 3: perfect bonding
                        ! 4: perfect debonding
                        case (0,1,3,4)
                          read(input_fileunit,*) boundary(i)%ctype(k,1)
                        ! Global axes
                        ! 2: partial bonding
                        ! Local axes
                        ! 5: partial bonding
                        case (2,5)
                          read(input_fileunit,*) boundary(i)%ctype(k,1), boundary(i)%cvalue_c(k,1,1), boundary(i)%cvalue_c(k,2,1), boundary(i)%cvalue_c(k,3,1)
                        case default
                          call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                      end select
                    end do
                    ! Check that the B.C. are for local or global axes, but not mixed.
                    if ((boundary(i)%ctype(1,1).eq.0).or.(boundary(i)%ctype(1,1).eq.1).or.(boundary(i)%ctype(1,1).eq.2)) then
                      do k=2,problem%n
                        if ((boundary(i)%ctype(k,1).ne.0).and.(boundary(i)%ctype(k,1).ne.1).and.(boundary(i)%ctype(k,1).ne.2)) then
                          call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'local axes B.C. and global axes B.C. can not be mixed.')
                        end if
                      end do
                    end if
                    if ((boundary(i)%ctype(1,1).eq.3).or.(boundary(i)%ctype(1,1).eq.4).or.(boundary(i)%ctype(1,1).eq.5)) then
                      do k=2,problem%n
                        if ((boundary(i)%ctype(k,1).ne.3).and.(boundary(i)%ctype(k,1).ne.4).and.(boundary(i)%ctype(k,1).ne.5)) then
                          call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'local axes B.C. and global axes B.C. can not be mixed.')
                        end if
                      end do
                    end if

                  ! ---------------------------------------------
                  ! POROELASTIC MEDIA (1) - POROELASTIC MEDIA (2)
                  ! ---------------------------------------------

                  case (fbem_poroelastic)
                    read(input_fileunit,*) boundary(i)%ctype(1,1)
                    select case (boundary(i)%ctype(1,1))
                      ! 0. perfectly permeable
                      ! 1. perfectly impermeable
                      case (0,1)
                      ! 2. partially permeable - k (flow resistance)
                      case (2)
                        call fbem_search_section(input_fileunit,section_name,found)
                        call fbem_search_keyword(input_fileunit,keyword,':',found)
                        read(input_fileunit,*) boundary(i)%ctype(1,1), boundary(i)%cvalue_c(1,1,1)
                      case default
                        call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                    end select

                end select
            end select

          ! ========================================================================================================================
          ! BE-FE BOUNDARY
          ! ========================================================================================================================

          case (fbem_boundary_coupling_be_fe)
            select case (boundary(i)%class)

              ! ================= !
              ! ORDINARY BOUNDARY !
              ! ================= !

              case (fbem_boundary_class_ordinary)
                select case (region(boundary(i)%region(1))%type)

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
                select case (region(boundary(i)%region(1))%type)

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

          ! ========================================================================================================================
          ! BE-FE-BE BOUNDARY
          ! ========================================================================================================================

          case (fbem_boundary_coupling_be_fe_be)
            select case (region(boundary(i)%region(1))%type)
              case (fbem_potential)
                select case (region(boundary(i)%region(2))%type)

                  ! ---------------------------------------
                  ! INVISCID FLUID (1) - INVISCID FLUID (2)
                  ! ---------------------------------------

                  case (fbem_potential)
                    ! It does not need B.C.

                  ! -------------------------------------------
                  ! INVISCID FLUID (1) - VISCOELASTIC SOLID (2)
                  ! -------------------------------------------

                  case (fbem_viscoelastic)
                    ! It does not need B.C.

                  ! ------------------------------------------
                  ! INVISCID FLUID (1) - POROELASTIC MEDIA (2)
                  ! ------------------------------------------

                  case (fbem_poroelastic)
                    ! It does not need B.C.

                end select

              case (fbem_viscoelastic)
                select case (region(boundary(i)%region(2))%type)

                  ! -------------------------------------------
                  ! VISCOELASTIC SOLID (1) - INVISCID FLUID (2)
                  ! -------------------------------------------

                  case (fbem_potential)
                    ! It does not need B.C.

                  ! -----------------------------------------------
                  ! VISCOELASTIC SOLID (1) - VISCOELASTIC SOLID (2)
                  ! -----------------------------------------------

                  case (fbem_viscoelastic)
                    ! It does not need B.C.

                  ! ----------------------------------------------
                  ! VISCOELASTIC SOLID (1) - POROELASTIC MEDIA (2)
                  ! ----------------------------------------------

                  case (fbem_poroelastic)
                    ! It does not need B.C.

                end select

              case (fbem_poroelastic)
                select case (region(boundary(i)%region(2))%type)

                  ! ------------------------------------------
                  ! POROELASTIC MEDIA (1) - INVISCID FLUID (2)
                  ! ------------------------------------------

                  case (fbem_potential)
                    ! It does not need B.C.

                  ! ----------------------------------------------
                  ! POROELASTIC MEDIA (1) - VISCOELASTIC SOLID (2)
                  ! ----------------------------------------------

                  case (fbem_viscoelastic)
                    ! It does not need B.C.

                  ! ---------------------------------------------
                  ! POROELASTIC MEDIA (1) - POROELASTIC MEDIA (2)
                  ! ---------------------------------------------

                  case (fbem_poroelastic)
                    ! It does not need B.C.

                end select
            end select
        end select
      end if
    end do ! Loop through BOUNDARIES

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  end if

end subroutine read_conditions_bem_boundaries_mechanics_harmonic
