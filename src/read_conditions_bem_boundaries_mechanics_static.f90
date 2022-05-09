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


subroutine read_conditions_bem_boundaries_mechanics_static(input_fileunit)

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
  real(kind=real64)                       :: center(problem%n), axis(3), theta

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
                      ! 0: u_k=U_k
                      ! 1. t_k=T_k
                      case (0,1)
                        read(input_fileunit,*) boundary(i)%ctype(k,1), boundary(i)%cvalue_r(k,1,1)
                      ! 4: Global axes: Infinitesimal rotation field : u_k = theta (axis x (node_position-center))·e_k
                      case (4)
                        select case (problem%n)
                          case (2)
                            read(input_fileunit,*) boundary(i)%ctype(k,1), center, theta
                            axis=[0.d0,0.d0,1.d0]
                          case (3)
                            read(input_fileunit,*) boundary(i)%ctype(k,1), center, axis, theta
                        end select
                        boundary(i)%cvalue_r(k,:,1)=0.d0
                        boundary(i)%cvalue_r(k,1:problem%n,1)=center
                        boundary(i)%cvalue_r(k,4:6,1)=axis/sqrt(dot_product(axis,axis))
                        boundary(i)%cvalue_r(k,7,1)=theta
                      ! 10: normal pressure
                      case (10)
                        read(input_fileunit,*) boundary(i)%ctype(1,1), boundary(i)%cvalue_r(1,1,1)
                        boundary(i)%ctype(:,1)=boundary(i)%ctype(1,1)
                        boundary(i)%cvalue_r(:,1,1)=boundary(i)%cvalue_r(1,1,1)
                        exit
                      case default
                        call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                    end select
                  end do

                ! ===================
                ! CRACK-LIKE BOUNDARY
                ! ===================

                case (fbem_boundary_class_cracklike)
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
                      ! 0: u_k^+=U_k^+
                      ! 1. t_k^+=T_k^+
                      case (0,1)
                        read(input_fileunit,*) boundary(i)%ctype(k,1), boundary(i)%cvalue_r(k,1,1)
                      ! 4: Global axes: Infinitesimal rotation field : u_k^+ = theta (axis x (node_position-center))·e_k
                      case (4)
                        select case (problem%n)
                          case (2)
                            read(input_fileunit,*) boundary(i)%ctype(k,1), center, theta
                            axis=[0.d0,0.d0,1.d0]
                          case (3)
                            read(input_fileunit,*) boundary(i)%ctype(k,1), center, axis, theta
                        end select
                        boundary(i)%cvalue_r(k,:,1)=0.d0
                        boundary(i)%cvalue_r(k,1:problem%n,1)=center
                        boundary(i)%cvalue_r(k,4:6,1)=axis/sqrt(dot_product(axis,axis))
                        boundary(i)%cvalue_r(k,7,1)=theta
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
                      ! 0: u_k^-=U_k^-
                      ! 1. t_k^-=T_k^-
                      case (0,1)
                        read(input_fileunit,*) boundary(i)%ctype(k,2), boundary(i)%cvalue_r(k,1,2)
                      ! 4: Global axes: Infinitesimal rotation field : u_k^+ = theta (axis x (node_position-center))·e_k
                      case (4)
                        select case (problem%n)
                          case (2)
                            read(input_fileunit,*) boundary(i)%ctype(k,2), center, theta
                            axis=[0.d0,0.d0,1.d0]
                          case (3)
                            read(input_fileunit,*) boundary(i)%ctype(k,2), center, axis, theta
                        end select
                        boundary(i)%cvalue_r(k,:,2)=0.d0
                        boundary(i)%cvalue_r(k,1:problem%n,2)=center
                        boundary(i)%cvalue_r(k,4:6,2)=axis/sqrt(dot_product(axis,axis))
                        boundary(i)%cvalue_r(k,7,2)=theta
                      case default
                        call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                    end select
                  end do

              end select

          ! ========================================================================================================================
          ! BE-BE BOUNDARY
          ! ========================================================================================================================

          case (fbem_boundary_coupling_be_be)
            ! It does not need B.C.

          ! ========================================================================================================================
          ! BE-FE BOUNDARY
          ! ========================================================================================================================

          case (fbem_boundary_coupling_be_fe)
            ! It does not need B.C.

          ! ========================================================================================================================
          ! BE-BE BOUNDARY
          ! ========================================================================================================================

          case (fbem_boundary_coupling_be_fe_be)
            ! It does not need B.C.

        end select
      end if
    end do ! Loop through BOUNDARIES

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  end if

end subroutine read_conditions_bem_boundaries_mechanics_static
