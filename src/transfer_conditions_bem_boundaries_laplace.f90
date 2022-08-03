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

subroutine transfer_conditions_bem_boundaries_laplace

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_string_handling
  use fbem_geometry
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! No implicit variables are allowed
  implicit none

  ! Local variables
  integer           :: kb, sp
  integer           :: kn, sn

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START transferring conditions (BEM boundaries->BEM nodes)')

  ! Loop through BOUNDARIES
  do kb=1,n_boundaries
    sp=boundary(kb)%part
    select case (boundary(kb)%coupling)

      ! ============================================================================================================================
      ! BE BOUNDARY
      ! ============================================================================================================================

      case (fbem_boundary_coupling_be)
        select case (boundary(kb)%class)

          ! =================
          ! ORDINARY BOUNDARY
          ! =================

          case (fbem_boundary_class_ordinary)
            do kn=1,part(sp)%n_nodes
              sn=part(sp)%node(kn)
              allocate (node(sn)%ctype(1,1))
              allocate (node(sn)%cvalue_r(1,1,1))
              node(sn)%ctype=boundary(kb)%ctype
              node(sn)%cvalue_r=boundary(kb)%cvalue_r
            end do

          ! ===================
          ! CRACK-LIKE BOUNDARY
          ! ===================

          case (fbem_boundary_class_cracklike)
            do kn=1,part(sp)%n_nodes
              sn=part(sp)%node(kn)
              allocate (node(sn)%ctype(1,2))
              allocate (node(sn)%cvalue_r(1,1,2))
              node(sn)%ctype=boundary(kb)%ctype
              node(sn)%cvalue_r=boundary(kb)%cvalue_r
            end do

        end select

      ! ============================================================================================================================
      ! BE-BE BOUNDARY
      ! ============================================================================================================================

      case (fbem_boundary_coupling_be_be)
        ! It does not need B.C.

      ! ============================================================================================================================
      ! BE-FE BOUNDARY
      ! ============================================================================================================================

      case (fbem_boundary_coupling_be_fe)
        stop 'not implemented yet'

      ! ============================================================================================================================
      ! BE-FE-BE BOUNDARY
      ! ============================================================================================================================

      case (fbem_boundary_coupling_be_fe_be)
        stop 'not implemented yet'

    end select
  end do ! Loop through BOUNDARIES

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END transferring conditions (BEM boundaries->BEM nodes)')

end subroutine transfer_conditions_bem_boundaries_laplace
