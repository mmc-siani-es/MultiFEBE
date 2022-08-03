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

subroutine transfer_conditions_bem_boundaries_mechanics_static

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
  integer           :: k
  real(kind=real64) :: center(3), axis(3), x(3), urot(3)

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
              allocate (node(sn)%ctype(problem%n,1))
              allocate (node(sn)%cvalue_r(problem%n,7,1))
              node(sn)%ctype=boundary(kb)%ctype
              node(sn)%cvalue_r=boundary(kb)%cvalue_r
              ! Transformation of infinitesimal rotation to displacement at each node
              do k=1,problem%n
                if (node(sn)%ctype(k,1).eq.4) then
                  center=node(sn)%cvalue_r(k,1:3,1)
                  axis=node(sn)%cvalue_r(k,4:6,1)
                  x=node(sn)%x
                  urot=fbem_rotation_infinitesimal_displacement(center,axis,1.d0,x)
                  node(sn)%ctype(k,1)=0
                  node(sn)%cvalue_r(k,1,1)=urot(k)*node(sn)%cvalue_r(k,7,1)
                end if
              end do
            end do

          ! ===================
          ! CRACK-LIKE BOUNDARY
          ! ===================

          case (fbem_boundary_class_cracklike)
            do kn=1,part(sp)%n_nodes
              sn=part(sp)%node(kn)
              allocate (node(sn)%ctype(problem%n,2))
              allocate (node(sn)%cvalue_r(problem%n,7,2))
              node(sn)%ctype=boundary(kb)%ctype
              node(sn)%cvalue_r=boundary(kb)%cvalue_r
              ! Transformation of infinitesimal rotation to displacement at each node
              ! Face +
              do k=1,problem%n
                if (node(sn)%ctype(k,1).eq.4) then
                  center=node(sn)%cvalue_r(k,1:3,1)
                  axis=node(sn)%cvalue_r(k,4:6,1)
                  x=node(sn)%x
                  urot=fbem_rotation_infinitesimal_displacement(center,axis,1.d0,x)
                  node(sn)%ctype(k,1)=0
                  node(sn)%cvalue_r(k,1,1)=urot(k)*node(sn)%cvalue_r(k,7,1)
                end if
              end do
              ! Face -
              do k=1,problem%n
                if (node(sn)%ctype(k,2).eq.4) then
                  center=node(sn)%cvalue_r(k,1:3,2)
                  axis=node(sn)%cvalue_r(k,4:6,2)
                  x=node(sn)%x
                  urot=fbem_rotation_infinitesimal_displacement(center,axis,1.d0,x)
                  node(sn)%ctype(k,2)=0
                  node(sn)%cvalue_r(k,1,2)=urot(k)*node(sn)%cvalue_r(k,7,2)
                end if
              end do
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
        ! It does not need B.C.

      ! ============================================================================================================================
      ! BE-FE-BE BOUNDARY
      ! ============================================================================================================================

      case (fbem_boundary_coupling_be_fe_be)
        ! It does not need B.C.

    end select
  end do ! Loop through BOUNDARIES

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END transferring conditions (BEM boundaries->BEM nodes)')

end subroutine transfer_conditions_bem_boundaries_mechanics_static
