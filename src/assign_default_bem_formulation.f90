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

subroutine assign_default_bem_formulation

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! No implicit variables are allowed
  implicit none

  ! Local variables
  integer :: ke, se, kn, knj, sn, snj, i, k
  real(kind=real64), parameter :: default_mca_boundary_delta  = 0.05d0
  real(kind=real64), parameter :: default_mca_linear_delta    = 0.42264973d0
  real(kind=real64), parameter :: default_mca_quadratic_delta = 0.22540333d0
  real(kind=real64), parameter :: default_mca_cubic_delta     = 0.138863688d0

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START assigning default BEM formulations')

  !
  ! Each node may have one of the following BEM formulations (BIE & collocation strategy)
  !
  ! SBIE                : SBIE with nodal collocation (standard)
  ! SBIE MCA            : SBIE with Multiple Collocation Approach (MCA) (approach for avoiding corner problems at doubled nodes)
  ! HBIE                : HBIE with Multiple Collocation Approach (MCA) (approach for avoiding C1 requirement at collocation point)
  ! DUAL BURTON & MILLER: SBIE+alpha*HBIE formulation for avoiding spurious results for external problems with internal cavities
  ! DUAL BOUNDARY       : Simultaneous use of SBIE and HBIE for crack-like boundaries
  !

  ! =============
  ! BE BOUNDARIES
  ! =============

  ! Allocate delta vectors on data structures
  do ke=1,n_elements
    if (part(element(ke)%part)%type.eq.fbem_part_be_boundary) then
      allocate (element(ke)%delta_sbie_mca(element(ke)%n_nodes))
      allocate (element(ke)%delta_hbie(element(ke)%n_nodes))
    end if
  end do

  do i=1,n_boundaries
    select case (boundary(i)%class)
      !
      ! ORDINARY BOUNDARY
      !
      case (fbem_boundary_class_ordinary)
        do kn=1,part(boundary(i)%part)%n_nodes
          sn=part(boundary(i)%part)%node(kn)
          !
          ! SBIE with nodal collocation by default, except for nodes at the boundary of the boundary in continuous elements, which
          ! have SBIE with MCA with a small delta.
          !
          node(sn)%sbie=fbem_sbie
          node(sn)%hbie=0
          node(sn)%dual=0
          node(sn)%dual_is_common=.false.
          if (node(sn)%in_boundary.and.(.not.element(node(sn)%element(1))%discontinuous)) then
            node(sn)%sbie=fbem_sbie_mca
            do ke=1,node(sn)%n_elements
              se=node(sn)%element(ke)
              k=node(sn)%element_node_iid(ke)
              element(se)%delta_sbie_mca(k)=default_mca_boundary_delta
            end do
          end if
        end do
      !
      ! CRACK-LIKE BOUNDARY
      !
      case (fbem_boundary_class_cracklike)
        do kn=1,part(boundary(i)%part)%n_nodes
          ! Selected node
          sn=part(boundary(i)%part)%node(kn)
          !
          ! If the node belongs to continuous element
          !
          if (element(node(sn)%element(1))%discontinuous.eqv.(.false.)) then
            !
            ! SBIE & HBIE (Dual BEM) with MCA formulation. SBIE and HBIE have the same collocation point.
            !
            node(sn)%sbie=fbem_sbie_mca
            node(sn)%hbie=fbem_hbie
            node(sn)%dual=fbem_dual_boundary
            node(sn)%dual_is_common=.true.
            do ke=1,node(sn)%n_elements
              se=node(sn)%element(ke)
              k=node(sn)%element_node_iid(ke)
              select case (fbem_element_order(element(se)%type))
                ! Linear elements
                case (1)
                  element(se)%delta_sbie_mca(k)=default_mca_linear_delta
                  element(se)%delta_hbie(k)    =default_mca_linear_delta
                ! Quadratic elements
                case (2)
                  element(se)%delta_sbie_mca(k)=default_mca_quadratic_delta
                  element(se)%delta_hbie(k)    =default_mca_quadratic_delta
                ! Cubic elements
                case (3)
                  element(se)%delta_sbie_mca(k)=default_mca_cubic_delta
                  element(se)%delta_hbie(k)    =default_mca_cubic_delta
                case default
                  call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid element order')
              end select
            end do
          !
          ! If the node belongs to a discontinuous element
          !
          else
            ! The node has SBIE / HBIE (Dual BEM) with nodal collocation formulation.
            node(sn)%sbie=fbem_sbie
            node(sn)%hbie=fbem_hbie
            node(sn)%dual=fbem_dual_boundary
            node(sn)%dual_is_common=.true.
          end if
        end do
    end select
  end do

  ! =====================
  ! COUPLED BE BODY LOADS
  ! =====================

  ! Allocate delta vectors on data structures
  do ke=1,n_elements
    if (part(element(ke)%part)%type.eq.fbem_part_be_bodyload) then
      if (be_bodyload(part(element(ke)%part)%entity)%coupling.ne.0) then
        allocate (element(ke)%delta_sbie_mca(element(ke)%n_nodes))
        allocate (element(ke)%delta_hbie(element(ke)%n_nodes))
      end if
    end if
  end do

  ! Loop through COUPLED BE BODY LOADS
  do i=1,n_be_bodyloads
    select case (be_bodyload(i)%coupling)
      !
      ! BEM-FEM beam tip coupling
      !
      case (fbem_bl_coupling_beam_tip)
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'not implemented yet')
      !
      ! BEM-FEM shell edge coupling
      !
      case (fbem_bl_coupling_shell_edge)
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'not implemented yet')
      !
      ! BEM-FEM beam line or shell surface coupling
      !
      case (fbem_bl_coupling_beam_line,fbem_bl_coupling_shell_surface)
        do kn=1,part(be_bodyload(i)%part)%n_nodes
          sn=part(be_bodyload(i)%part)%node(kn)
          !
          ! SBIE with nodal collocation by default, except for nodes at the boundary of the element, which have SBIE with MCA with
          ! a small delta.
          !
          node(sn)%sbie=fbem_sbie
          node(sn)%hbie=0
          node(sn)%dual=0
          node(sn)%dual_is_common=.false.
          if (node(sn)%in_boundary) then
            node(sn)%sbie=fbem_sbie_mca
            do ke=1,node(sn)%n_elements
              se=node(sn)%element(ke)
              k=node(sn)%element_node_iid(ke)
              select case (fbem_element_order(element(se)%type))
                ! Linear elements
                case (1)
                  element(se)%delta_sbie_mca(k)=default_mca_linear_delta
                ! Quadratic elements
                case (2)
                  element(se)%delta_sbie_mca(k)=default_mca_quadratic_delta
                ! Cubic elements
                case (3)
                  element(se)%delta_sbie_mca(k)=default_mca_cubic_delta
              end select
            end do
            !
            ! Since it is not possible (as far as we know) to perform singular integrals along 3D line body loads (occur when end
            ! points coincides with a boundary element node), we use a MCA collocation for those boundary element nodes.
            !
            if (be_bodyload(i)%coupling.eq.fbem_bl_coupling_beam_line) then
              do knj=1,node(sn)%n_nodes
                snj=node(sn)%node(knj)
                !
                ! By selecting nodes with SBIE formulation it is assumed that nodal collocation was assigned to the node
                ! and it is necessary to change it to MCA.
                !
                if (node(snj)%sbie.eq.fbem_sbie) then
                  node(snj)%sbie=fbem_sbie_mca
                  do ke=1,node(snj)%n_elements
                    se=node(snj)%element(ke)
                    k=node(snj)%element_node_iid(ke)
                    select case (fbem_element_order(element(se)%type))
                      case (1)
                        element(se)%delta_sbie_mca(k)=default_mca_linear_delta
                      case (2)
                        element(se)%delta_sbie_mca(k)=default_mca_quadratic_delta
                      case (3)
                        element(se)%delta_sbie_mca(k)=default_mca_cubic_delta
                    end select
                  end do
                end if
              end do
            end if
            !
            ! Later it should be checked that node sn has a delta which produces collocation points inside the region
            !
          end if
        end do

    end select
  end do

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END assigning default BEM formulations')

end subroutine assign_default_bem_formulation
