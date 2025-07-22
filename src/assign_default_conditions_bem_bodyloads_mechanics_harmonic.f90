! ---------------------------------------------------------------------
! Copyright (C) 2014-2024 Universidad de Las Palmas de Gran Canaria:
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

subroutine assign_default_conditions_bem_bodyloads_mechanics_harmonic

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! No implicit variables are allowed
  implicit none

  ! Local variables
  integer :: i

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START assigning default conditions (BE body loads)')

  ! Loop through BE BODY LOADS
  do i=1,n_be_bodyloads

    select case (be_bodyload(i)%coupling)

      case (fbem_bl_uncoupled)
        select case (region(be_bodyload(i)%region)%type)

          ! --------------
          ! INVISCID FLUID
          ! --------------

          case (fbem_potential)
            allocate (be_bodyload(i)%ctype(1,1))
            allocate (be_bodyload(i)%cvalue_c(1,1,1))
            !
            ! ctype=0: constant amplitude
            !
            be_bodyload(i)%ctype=0    ! constant body load
            be_bodyload(i)%cvalue_c=0 ! null body load

          ! ------------------
          ! VISCOELASTIC SOLID
          ! ------------------

          case (fbem_viscoelastic)
            allocate (be_bodyload(i)%ctype(problem%n,1))
            allocate (be_bodyload(i)%cvalue_c(problem%n,1,1))
            be_bodyload(i)%ctype=0    ! constant body load
            be_bodyload(i)%cvalue_c=0 ! null body load

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

  end do

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END assigning default conditions (BEM boundaries)')

end subroutine assign_default_conditions_bem_bodyloads_mechanics_harmonic
