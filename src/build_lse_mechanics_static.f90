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

subroutine build_lse_mechanics_static

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_numerical
  use fbem_shape_functions
  use fbem_geometry

  ! Module of problem variables
  use problem_variables

  ! No implicit variables
  implicit none

  ! Local variables
  integer                        :: kr
  integer                        :: kn
  integer                        :: knm
  real(kind=real64), allocatable :: T(:,:)
  integer                        :: kc
  integer                        :: knc
  integer                        :: row

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'START building the linear system of equations')

  ! Initialization
  A_r=0.d0
  b_r=0.d0

  ! ==========================================================================================
  ! BUILD AND ASSEMBLE BEM INFLUENCE MATRICES AND FEM STIFFNESS MATRICES AND DISTRIBUTED LOADS
  ! ==========================================================================================

  do kr=1,n_regions
    select case (region(kr)%class)
      case (fbem_be)
        call build_lse_mechanics_bem_staela(kr)
      case (fbem_fe)
        call build_lse_mechanics_fem_staela(kr)
    end select
  end do

  ! ===================================
  ! ASSEMBLE FINITE ELEMENT NODAL LOADS
  ! ===================================

  do kn=1,n_nodes

    if (part(node(kn)%part(1))%type.eq.fbem_part_fe_subregion) then

      if (node(kn)%coupled_node.eq.0) then

        !
        ! Slave nodes
        !
        if (node(kn)%rigid_link.eq.2) then
          knm=node(kn)%master
          allocate (T(node(kn)%n_dof,node(knm)%n_dof))
          call fbem_rigid_solid_transformation_matrix(problem%n,node(kn)%x,node(knm)%x,node(kn)%n_dof,T)
          do kc=1,node(knm)%n_dof
            if (node(knm)%ctype(kc,1).eq.1) then
              row=node(knm)%row(kc,1)
              do knc=1,node(kn)%n_dof
                if (node(kn)%ctype(knc,1).eq.1) then
                  b_r(row,1)=b_r(row,1)+T(knc,kc)*node(kn)%cvalue_r(knc,1,1)
                else
                  call fbem_error_message(error_unit,0,'node',node(kn)%id,'slave nodes do not admit kinematic constraints')
                end if
              end do
            end if
          end do
          deallocate (T)
        !
        ! Normal and master nodes
        !
        else
          do kc=1,node(kn)%n_dof
            if (node(kn)%ctype(kc,1).eq.1) then
              row=node(kn)%row(kc,1)
              b_r(row,1)=b_r(row,1)+node(kn)%cvalue_r(kc,1,1)
            end if
          end do
        end if

      !
      ! Coupled nodal load
      !
      else

        stop 'not implemented yet'

      end if

    end if
  end do

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'END building the linear system of equations')

end subroutine build_lse_mechanics_static
