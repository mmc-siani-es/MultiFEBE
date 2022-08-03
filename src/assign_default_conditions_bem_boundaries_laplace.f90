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

subroutine assign_default_conditions_bem_boundaries_laplace

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

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START assigning default conditions (BEM boundaries)')

  do i=1,n_boundaries
    select case (boundary(i)%coupling)

      ! ===========
      ! BE BOUNDARY
      ! ===========

      case (fbem_boundary_coupling_be)
        select case (boundary(i)%class)

          ! =================
          ! ORDINARY BOUNDARY
          ! =================

          case (fbem_boundary_class_ordinary)
            ! Allocate structure members
            allocate (boundary(i)%ctype(1,1))
            allocate (boundary(i)%cvalue_r(1,1,1))
            ! j=k·dp/dn=0
            boundary(i)%ctype(1,1)=1
            boundary(i)%cvalue_r(1,1,1)=0.d0

          ! ===================
          ! CRACK-LIKE BOUNDARY
          ! ===================

          case (fbem_boundary_class_cracklike)
            ! Allocate structure members
            allocate (boundary(i)%ctype(1,2))
            allocate (boundary(i)%cvalue_r(1,1,2))
            ! Face + : j=k·dp/dn=0
            boundary(i)%ctype(1,1)=1
            boundary(i)%cvalue_r(1,1,1)=0.d0
            ! Face - : j=k·dp/dn=0
            boundary(i)%ctype(1,2)=1
            boundary(i)%cvalue_r(1,1,2)=0.d0
        end select

      ! ==============
      ! BE-BE BOUNDARY
      ! ==============

      case (fbem_boundary_coupling_be_be)
        ! It does not need B.C.

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

  end do

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END assigning default conditions (BEM boundaries)')

end subroutine assign_default_conditions_bem_boundaries_laplace
