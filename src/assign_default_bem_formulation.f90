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
  integer :: ke, se, kn, sn, i, k

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START assigning default BEM formulations')

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

  ! Loop through the BE BOUNDARIES
  do i=1,n_boundaries
    !
    ! Select between boundary class
    !
    select case (boundary(i)%class)
      !
      ! ORDINARY BOUNDARY
      !
      case (fbem_boundary_class_ordinary)
        ! Loop through the nodes of the boundary
        do kn=1,part(boundary(i)%part)%n_nodes
          ! Selected node
          sn=part(boundary(i)%part)%node(kn)
          ! If the node belongs to continuous elements
          if (element(node(sn)%element(1))%discontinuous.eqv.(.false.)) then
!            !
!            ! SBIE by default
!            !
!            node(sn)%sbie=fbem_sbie
!            node(sn)%hbie=0
!            node(sn)%dual=0
!            node(sn)%dual_is_common=.false.
            !
            ! SBIE by default, except nodes at the boundary of the boundary, which have SBIE MCA formulation.
            !
            ! The node has SBIE formulation
            node(sn)%sbie=fbem_sbie
            node(sn)%hbie=0
            node(sn)%dual=0
            node(sn)%dual_is_common=.false.
            ! If it is in the boundary of the boundary
            if (node(sn)%in_boundary) then
              ! The node has SBIE MCA formulation
              node(sn)%sbie=fbem_sbie_mca
              ! The displacement towards inside each element of the node for SBIE MCA formulation
              ! Loop through the elements of the node
              do ke=1,node(sn)%n_elements
                ! Selected element
                se=node(sn)%element(ke)
                ! Index of the node in the selected element
                k=node(sn)%element_node_iid(ke)
                ! It depends on the element type
                select case (element(se)%type)
                  ! Linear elements
                  case (fbem_line2,fbem_tri3,fbem_quad4)
                    element(se)%delta_sbie_mca(k)=0.42264973d0
                  ! Quadratic elements
                  case (fbem_line3,fbem_tri6,fbem_quad8,fbem_quad9)
                    element(se)%delta_sbie_mca(k)=0.22540333d0
                  ! Cubic elements
                  case (fbem_line4)
                    element(se)%delta_sbie_mca(k)=0.138863688d0
                end select
              end do
            end if
!            !
!            ! SBIE MCA by default
!            !
!            ! The node has SBIE MCA formulation
!            node(sn)%sbie=fbem_sbie_mca
!            node(sn)%hbie=0
!            node(sn)%dual=0
!            node(sn)%dual_is_common=.false.
!            ! The displacement towards inside each element of the node for SBIE MCA formulation
!            ! Loop through the elements of the node
!            do ke=1,node(sn)%n_elements
!              ! Selected element
!              se=node(sn)%element(ke)
!              ! Index of the node in the selected element
!              k=node(sn)%element_node_iid(ke)
!              ! It depends on the element type
!              select case (element(se)%type)
!                ! Linear elements
!                case (fbem_line2,fbem_tri3,fbem_quad4)
!                  element(se)%delta_sbie_mca(k)=0.42264973d0
!                ! Quadratic elements
!                case (fbem_line3,fbem_tri6,fbem_quad8,fbem_quad9)
!                  element(se)%delta_sbie_mca(k)=0.22540333d0
!                ! Cubic elements
!                case (fbem_line4)
!                  element(se)%delta_sbie_mca(k)=0.138863688d0
!              end select
!            end do
!            !
!            ! HBIE by default
!            !
!            ! The node has HBIE MCA formulation
!            node(sn)%sbie=0
!            node(sn)%hbie=fbem_hbie
!            node(sn)%dual=0
!            node(sn)%dual_is_common=.false.
!            ! The displacement towards inside each element of the node for HBIE MCA formulation
!            ! Loop through the elements of the node
!            do ke=1,node(sn)%n_elements
!              ! Selected element
!              se=node(sn)%element(ke)
!              ! Index of the node in the selected element
!              k=node(sn)%element_node_iid(ke)
!              ! It depends on the element type
!              select case (element(se)%type)
!                ! Linear elements
!                case (fbem_line2,fbem_tri3,fbem_quad4)
!                  element(se)%delta_hbie(k)=0.42264973d0
!                ! Quadratic elements
!                case (fbem_line3,fbem_tri6,fbem_quad8,fbem_quad9)
!                  element(se)%delta_hbie(k)=0.22540333d0
!                ! Cubic elements
!                case (fbem_line4)
!                  element(se)%delta_hbie(k)=0.138863688d0
!              end select
!            end do
!            !
!            ! SBIE MCA + HBIE MCA (Burton & Miller) by default
!            !
!            node(sn)%sbie=fbem_sbie_mca
!            node(sn)%hbie=fbem_hbie
!            node(sn)%dual=fbem_dual_burton_miller
!            node(sn)%dual_is_common=.true.
!            node(sn)%alpha=1.0d0
!            ! The displacement towards inside each element of the node for HBIE MCA formulation
!            ! Loop through the elements of the node
!            do ke=1,node(sn)%n_elements
!              ! Selected element
!              se=node(sn)%element(ke)
!              ! Index of the node in the selected element
!              k=node(sn)%element_node_iid(ke)
!              ! It depends on the element type
!              select case (element(se)%type)
!                ! Linear elements
!                case (fbem_line2,fbem_tri3,fbem_quad4)
!                  element(se)%delta_sbie_mca(k)=0.42264973d0
!                  element(se)%delta_hbie(k)=0.42264973d0
!                ! Quadratic elements
!                case (fbem_line3,fbem_tri6,fbem_quad8,fbem_quad9)
!                  element(se)%delta_sbie_mca(k)=0.22540333d0
!                  element(se)%delta_hbie(k)=0.22540333d0
!                ! Cubic elements
!                case (fbem_line4)
!                  element(se)%delta_sbie_mca(k)=0.138863688d0
!                  element(se)%delta_hbie(k)=0.138863688d0
!              end select
!            end do
          !
          ! If the node belongs to a discontinuous element
          !
          else
            !
            ! SBIE by default
            !
            node(sn)%sbie=fbem_sbie
            node(sn)%hbie=0
            node(sn)%dual=0
            node(sn)%dual_is_common=.false.
!            !
!            ! HBIE by default
!            !
!            node(sn)%sbie=0
!            node(sn)%hbie=fbem_hbie
!            node(sn)%dual=0
!            node(sn)%dual_is_common=.false.
!            !
!            ! SBIE + HBIE (Burton & Miller) by default
!            !
!            node(sn)%sbie=fbem_sbie
!            node(sn)%hbie=fbem_hbie
!            node(sn)%dual=fbem_dual_burton_miller
!            node(sn)%dual_is_common=.true.
!            node(sn)%alpha=1.0d0
          end if
        end do
      !
      ! CRACK-LIKE BOUNDARY
      !
      case (fbem_boundary_class_cracklike)
        ! Loop through the nodes of the boundary
        do kn=1,part(boundary(i)%part)%n_nodes
          ! Selected node
          sn=part(boundary(i)%part)%node(kn)
          ! If the node belongs to continuous elements
          if (element(node(sn)%element(1))%discontinuous.eqv.(.false.)) then
            !
            ! The node has SBIE MCA / HBIE (Dual BEM) formulation. SBIE and HBIE have the same collocation point.
            !
            node(sn)%sbie=fbem_sbie_mca
            node(sn)%hbie=fbem_hbie
            node(sn)%dual=fbem_dual_boundary
            node(sn)%dual_is_common=.true.
            ! The displacement towards inside each element of the node for HBIE formulation
            ! Loop through the elements of the node
            do ke=1,node(sn)%n_elements
              ! Selected element
              se=node(sn)%element(ke)
              ! Index of the node in the selected element
              k=node(sn)%element_node_iid(ke)
              ! It depends on the element type
              select case (element(se)%type)
                ! Linear elements
                case (fbem_line2,fbem_tri3,fbem_quad4)
                  element(se)%delta_sbie_mca(k)=0.42264973d0
                  element(se)%delta_hbie(k)    =0.42264973d0
                ! Quadratic elements
                case (fbem_line3,fbem_tri6,fbem_quad8,fbem_quad9)
                  element(se)%delta_sbie_mca(k)=0.22540333d0
                  element(se)%delta_hbie(k)    =0.22540333d0
                ! Cubic elements
                case (fbem_line4)
                  element(se)%delta_sbie_mca(k)=0.138863688d0
                  element(se)%delta_hbie(k)    =0.138863688d0
              end select
            end do
!            !
!            ! The node has SBIE / HBIE (Dual BEM) formulation
!            !
!            node(sn)%sbie=fbem_sbie
!            node(sn)%hbie=fbem_hbie
!            node(sn)%dual=fbem_dual_boundary
!            node(sn)%dual_is_common=.false.
!            ! The displacement towards inside each element of the node for HBIE formulation
!            ! Loop through the elements of the node
!            do ke=1,node(sn)%n_elements
!              ! Selected element
!              se=node(sn)%element(ke)
!              ! Index of the node in the selected element
!              k=node(sn)%element_node_iid(ke)
!              ! It depends on the element type
!              select case (element(se)%type)
!                ! Linear elements
!                case (fbem_line2,fbem_tri3,fbem_quad4)
!                  element(se)%delta_hbie(k)=0.42264973d0
!                ! Quadratic elements
!                case (fbem_line3,fbem_tri6,fbem_quad8,fbem_quad9)
!                  element(se)%delta_hbie(k)=0.22540333d0
!                ! Cubic elements
!                case (fbem_line4)
!                  element(se)%delta_hbie(k)=0.138863688d0
!              end select
!            end do
!            ! If it is in the boundary of the boundary
!            if (node(sn)%in_boundary) then
!              ! Instead of having SBIE formulation, the node has SBIE MCA formulation
!              node(sn)%sbie=fbem_sbie_mca
!              node(sn)%dual_is_common=.true.
!              ! The displacement towards inside each element of the node for SBIE MCA formulation
!              ! Loop through the elements of the node
!              do ke=1,node(sn)%n_elements
!                ! Selected element
!                se=node(sn)%element(ke)
!                ! Index of the node in the selected element
!                k=node(sn)%element_node_iid(ke)
!                ! It depends on the element type
!                select case (element(se)%type)
!                  ! Linear elements
!                  case (fbem_line2,fbem_tri3,fbem_quad4)
!                    element(se)%delta_sbie_mca(k)=0.42264973d0
!                  ! Quadratic elements
!                  case (fbem_line3,fbem_tri6,fbem_quad8,fbem_quad9)
!                    element(se)%delta_sbie_mca(k)=0.22540333d0
!                  ! Cubic elements
!                  case (fbem_line4)
!                    element(se)%delta_sbie_mca(k)=0.138863688d0
!                end select
!              end do
!            end if
          !
          ! If the node belongs to a discontinuous element
          !
          else
            ! The node has SBIE / HBIE (Dual BEM) formulation.
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
      case (fbem_bl_coupling_beam_tip)
        stop 'assign_default_bem_formulation:fbem_bl_coupling_beam_tip'
      case (fbem_bl_coupling_shell_edge)
        stop 'assign_default_bem_formulation:fbem_bl_coupling_shell_edge'
      case (fbem_bl_coupling_beam_line,fbem_bl_coupling_shell_surface)
        ! Loop through the nodes of the boundary
        do kn=1,part(be_bodyload(i)%part)%n_nodes
          ! Selected node
          sn=part(be_bodyload(i)%part)%node(kn)
          ! The node has SBIE formulation
          node(sn)%sbie=fbem_sbie
          ! HBIE & DUAL N/A here
          node(sn)%hbie=0
          node(sn)%dual=0
          node(sn)%dual_is_common=.false.
          ! If it is in the boundary of the boundary
          if (node(sn)%in_boundary) then
            ! The node has SBIE MCA formulation
            node(sn)%sbie=fbem_sbie_mca
            ! The displacement towards inside each element of the node for SBIE MCA formulation
            ! Loop through the elements of the node
            do ke=1,node(sn)%n_elements
              ! Selected element
              se=node(sn)%element(ke)
              ! Index of the node in the selected element
              k=node(sn)%element_node_iid(ke)
              ! It depends on the element type
              select case (element(se)%type)
                ! Linear elements
                case (fbem_line2,fbem_tri3,fbem_quad4)
                  element(se)%delta_sbie_mca(k)=0.42264973d0
                ! Quadratic elements
                case (fbem_line3,fbem_tri6,fbem_quad8,fbem_quad9)
                  element(se)%delta_sbie_mca(k)=0.22540333d0
                ! Cubic elements
                case (fbem_line4)
                  element(se)%delta_sbie_mca(k)=0.138863688d0
              end select
            end do
          end if
        end do

    end select
  end do

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END assigning default BEM formulations')

end subroutine assign_default_bem_formulation
