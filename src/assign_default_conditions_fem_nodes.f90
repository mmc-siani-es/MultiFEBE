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

subroutine assign_default_conditions_fem_nodes

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_geometry
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_symmetry
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! No implicit variables are allowed
  implicit none

  ! Local variables
  integer :: i, k

  ! Starting message
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START assigning default conditions (FEM nodes)')

  !
  ! node(i)%ctype(k,1)=0 means not active DOF (with a given value)
  ! node(i)%ctype(k,1)=1 means active DOF (with a punctual load)    [by default]
  !
  ! FEM degrees of freedom indices:
  !
  ! maximum number of dof = 3*(problem%n-1)
  !
  ! Two types of nodes exist:
  !   1 without rotation: n_dof_node=problem%n
  !   2 with rotation:    n_dof_node=3*(problem%n-1) [all nodes are allocated this way]
  !
  ! 2D problems:
  !   1 u_1
  !   2 u_2
  !   3 theta (rotation 1-2)
  !
  ! 3D problems:
  !   1 u_1
  !   2 u_2
  !   3 u_3
  !   4 theta_1 (rotation 2-3) [if local theta_1' (rotation 2'-3')] : alpha rotation  (shells)
  !   5 theta_2 (rotation 3-1) [if local theta_2' (rotation 3'-1')] : beta rotation (shells)
  !   6 theta_3 (rotation 1-2) [if local theta_3' (rotation 1'-2')]
  !
  select case (problem%type)
    !
    ! Mechanics
    !
    case (fbem_mechanics)
      select case (problem%analysis)
        !
        ! Static
        !
        case (fbem_static)
          ! Loop through nodes
          do i=1,n_nodes
            if (part(node(i)%part(1))%type.eq.fbem_part_fe_subregion) then
              allocate (node(i)%ctype(3*(problem%n-1),1))
              allocate (node(i)%cvalue_r(3*(problem%n-1),1,1))
              do k=1,3*(problem%n-1)
                node(i)%ctype(k,1)=1
                node(i)%cvalue_r(k,1,1)=0
              end do
            end if
          end do
        !
        ! Harmonic
        !
        case (fbem_harmonic)
          ! Loop through nodes
          do i=1,n_nodes
            if (part(node(i)%part(1))%type.eq.fbem_part_fe_subregion) then
              allocate (node(i)%ctype(3*(problem%n-1),1))
              allocate (node(i)%cvalue_c(3*(problem%n-1),1,1))
              do k=1,3*(problem%n-1)
                node(i)%ctype(k,1)=1
                node(i)%cvalue_c(k,1,1)=0
              end do
            end if
          end do
      end select
  end select

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END assigning default conditions (FEM nodes)')

end subroutine assign_default_conditions_fem_nodes
