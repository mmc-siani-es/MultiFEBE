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

program multifebe

  use iso_fortran_env
  use fbem_string_handling
  use fbem_data_structures
  use fbem_shape_functions
  use problem_variables

  implicit none
  integer :: kf

  call date_and_time(timestamp_date_start,timestamp_time_start)

  call process_command_line_options

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'START multifebe')

  call fbem_init_shape_functions_module

  call read_input_file

  call build_data

  select case (problem%type)

    ! ===============
    ! LAPLACE PROBLEM
    ! ===============

    case (fbem_laplace)

      ! TO-DO: check and make documentation

      call build_auxiliary_variables_laplace
      call build_lse_laplace
      call solve_lse_r(n_dof,A_r,fact_ipiv,Aodim,Ao_r,scal_equed,scal_r,scal_c,1,b_r,.true.,lse_scaling,lse_condition,lse_refine)
      call assign_solution_laplace
      call calculate_internal_points_laplace
      call export_solution_laplace
      if (problem%sensitivity) then
        call build_lse_laplace_sa
        call solve_lse_r(n_dof,A_r,fact_ipiv,Aodim,Ao_r,scal_equed,scal_r,scal_c,problem%n_designvariables,bsa_r,&
                        .false.,lse_scaling,lse_condition,lse_refine)
        call assign_solution_laplace_sa
        call calculate_internal_points_laplace_sa
        call export_solution_laplace_sa
      end if

    ! ======================================
    ! CONTINUUM/STRUCTURAL MECHANICS PROBLEM
    ! ======================================

    case (fbem_mechanics)

      select case (problem%analysis)

        ! ---------------
        ! STATIC ANALYSIS
        ! ---------------

        case (fbem_static)

          call build_auxiliary_variables_mechanics_static
          call build_lse_mechanics_static
          call solve_lse_r(n_dof,A_r,fact_ipiv,Aodim,Ao_r,scal_equed,scal_r,scal_c,1,b_r,&
                           .true.,lse_scaling,lse_condition,lse_refine)
          call assign_solution_mechanics_static
          call calculate_stresses_mechanics_static
          call calculate_internal_points_mechanics_static
          call export_solution_mechanics_static
          if (problem%sensitivity) then
            call build_lse_mechanics_static_sa
            call solve_lse_r(n_dof,A_r,fact_ipiv,Aodim,Ao_r,scal_equed,scal_r,scal_c,problem%n_designvariables,bsa_r,&
                             .false.,lse_scaling,lse_condition,lse_refine)
            call assign_solution_mechanics_static_sa
            call export_solution_mechanics_static_sa
          end if

        ! ----------------------
        ! TIME-HARMONIC ANALYSIS
        ! ----------------------

        case (fbem_harmonic)

          call build_auxiliary_variables_mechanics_harmonic
          do kf=1,n_frequencies
            call print_frequency(output_unit,1,kf)
            call calculate_incident_mechanics_harmonic(kf)
            call build_lse_mechanics_harmonic(kf)
            call solve_lse_c(n_dof,A_c,fact_ipiv,Aodim,Ao_c,scal_equed,scal_r,scal_c,1,b_c,&
                             .true.,lse_scaling,lse_condition,lse_refine)
            call assign_solution_mechanics_harmonic(kf)
            call calculate_stresses_mechanics_harmonic(kf)
            call calculate_internal_points_mechanics_harmonic(kf)
            call export_solution_mechanics_harmonic(kf)
            if (problem%sensitivity) then
              call build_lse_mechanics_harmonic_sa(kf)
              call solve_lse_c(n_dof,A_c,fact_ipiv,Aodim,Ao_c,scal_equed,scal_r,scal_c,problem%n_designvariables,bsa_c,&
                               .false.,lse_scaling,lse_condition,lse_refine)
              call assign_solution_mechanics_harmonic_sa(kf)
              call export_solution_mechanics_harmonic_sa(kf)
            end if
          end do

        case default

          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid problem analysis')

      end select

    case default

      call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid problem type')

  end select

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'END multifebe')

end program multifebe
