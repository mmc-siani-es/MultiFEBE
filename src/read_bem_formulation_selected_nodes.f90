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


subroutine read_bem_formulation_selected_nodes(input_fileunit,section_name,keyword,selected_n_nodes,selected_node)
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
  ! I/O variables
  integer                               :: input_fileunit    !! Input file unit
  character(len=fbem_stdcharlen)        :: section_name      ! Name of the section
  character(len=fbem_stdcharlen)        :: keyword
  integer                               :: selected_n_nodes
  integer                               :: selected_node(selected_n_nodes)
  ! Local variables
  logical                               :: found
  integer                               ::  k
  integer                               :: sb, ke, se, kn, sn
  character(len=fbem_stdcharlen)        :: tmp_type
  real(kind=real64)                     :: tmp_delta_sbie, tmp_delta_hbie
  integer                               :: tmp_type_value


  ! Nota: La dual B&M poroelastica no esta bien!!, la suma hay que ponderarla con otro numero de onda o hacerla de otra
  ! manera...


  ! Boundary of selected nodes
  sb=part(node(selected_node(1))%part(1))%entity
  select case (boundary(sb)%class)

    ! ------------------------------------------------------------------------------------------------------------------------
    ! ORDINARY BOUNDARY
    ! ------------------------------------------------------------------------------------------------------------------------

    case (fbem_boundary_class_ordinary)
      ! Read the type of formulation
      read(input_fileunit,*) tmp_type
      tmp_type_value=0
      if (trim(tmp_type).eq.'sbie'             ) tmp_type_value=1
      if (trim(tmp_type).eq.'sbie_boundary_mca') tmp_type_value=2
      if (trim(tmp_type).eq.'sbie_mca'         ) tmp_type_value=3
      if (trim(tmp_type).eq.'hbie'             ) tmp_type_value=4
      if (trim(tmp_type).eq.'dual_bm_same'     ) tmp_type_value=5
      if (trim(tmp_type).eq.'dual_bm_mixed'    ) tmp_type_value=6
      ! Check if valid
      if (tmp_type_value.eq.0) then
        call fbem_error_message(error_unit,0,'boundary',boundary(sb)%id,'the indicated BEM formulation does not exist.')
      end if
      ! Read parameters of the formulation if needed
      select case (tmp_type_value)
        ! Read delta_sbie
        case (2,3)
          call fbem_search_section(input_fileunit,section_name,found)
          call fbem_search_keyword(input_fileunit,keyword,':',found)
          read(input_fileunit,*) tmp_type, tmp_delta_sbie
        ! Read delta_hbie
        case (4)
          call fbem_search_section(input_fileunit,section_name,found)
          call fbem_search_keyword(input_fileunit,keyword,':',found)
          read(input_fileunit,*) tmp_type, tmp_delta_hbie
        ! Read delta_hbie (=delta_sbie)
        case (5)
          call fbem_search_section(input_fileunit,section_name,found)
          call fbem_search_keyword(input_fileunit,keyword,':',found)
          read(input_fileunit,*) tmp_type, tmp_delta_hbie
        ! Read delta_sbie y delta_hbie
        case (6)
          call fbem_search_section(input_fileunit,section_name,found)
          call fbem_search_keyword(input_fileunit,keyword,':',found)
          read(input_fileunit,*) tmp_type, tmp_delta_sbie, tmp_delta_hbie
      end select
      ! Loop through the nodes of the boundary
      do kn=1,selected_n_nodes
        ! Selected node
        sn=selected_node(kn)
        ! If the node belongs to continuous elements
        if (element(node(sn)%element(1))%discontinuous.eqv.(.false.)) then
          select case (tmp_type_value)
            !
            ! SBIE:
            ! All nodes with nodal collocation. This is unsafe since at doubled nodes, depending on the boundary conditions, it
            ! will lead to a singular linear system of equations.
            !
            case (1)
              node(sn)%sbie=fbem_sbie
              node(sn)%hbie=0
              node(sn)%dual=0
              node(sn)%dual_is_common=.false.
            !
            ! SBIE_BOUNDARY_MCA:
            ! All nodes with nodal collocation, except nodes at the boundary of the boundary, which have SBIE MCA.
            !
            case (2)
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
                  ! Depending on tmp_delta_sbie, use the automatic value
                  if (tmp_delta_sbie.le.0.0d0) then
                    ! It depends on the element type
                    select case (element(se)%type)
                      ! Linear elements
                      case (fbem_line2,fbem_tri3,fbem_quad4)
                        element(se)%delta_sbie_mca(k)=0.42264973d0
                      ! Quadratic elements
                      case (fbem_line3,fbem_tri6,fbem_quad8,fbem_quad9)
                        element(se)%delta_sbie_mca(k)=0.22540333d0
                    end select
                  else
                    element(se)%delta_sbie_mca(k)=tmp_delta_sbie
                  end if
                end do
              end if
            !
            ! SBIE_MCA:
            ! All nodes with SBIE MCA collocation.
            !
            case (3)
              ! The node has SBIE MCA formulation
              node(sn)%sbie=fbem_sbie_mca
              node(sn)%hbie=0
              node(sn)%dual=0
              node(sn)%dual_is_common=.false.
              ! The displacement towards inside each element of the node for SBIE MCA formulation
              ! Loop through the elements of the node
              do ke=1,node(sn)%n_elements
                ! Selected element
                se=node(sn)%element(ke)
                ! Index of the node in the selected element
                k=node(sn)%element_node_iid(ke)
                ! Depending on tmp_delta_sbie, use the automatic value
                if (tmp_delta_sbie.le.0.0d0) then
                  ! It depends on the element type
                  select case (element(se)%type)
                    ! Linear elements
                    case (fbem_line2,fbem_tri3,fbem_quad4)
                      element(se)%delta_sbie_mca(k)=0.42264973d0
                    ! Quadratic elements
                    case (fbem_line3,fbem_tri6,fbem_quad8,fbem_quad9)
                      element(se)%delta_sbie_mca(k)=0.22540333d0
                  end select
                else
                  element(se)%delta_sbie_mca(k)=tmp_delta_sbie
                end if
              end do
            !
            ! HBIE:
            ! All nodes with HBIE collocation.
            !
            case (4)
              ! The node has HBIE MCA formulation
              node(sn)%sbie=0
              node(sn)%hbie=fbem_hbie
              node(sn)%dual=0
              node(sn)%dual_is_common=.false.
              ! The displacement towards inside each element of the node for HBIE MCA formulation
              ! Loop through the elements of the node
              do ke=1,node(sn)%n_elements
                ! Selected element
                se=node(sn)%element(ke)
                ! Index of the node in the selected element
                k=node(sn)%element_node_iid(ke)
                ! Depending on tmp_delta_hbie, use the automatic value
                if (tmp_delta_hbie.le.0.0d0) then
                  ! It depends on the element type
                  select case (element(se)%type)
                    ! Linear elements
                    case (fbem_line2,fbem_tri3,fbem_quad4)
                      element(se)%delta_hbie(k)=0.42264973d0
                    ! Quadratic elements
                    case (fbem_line3,fbem_tri6,fbem_quad8,fbem_quad9)
                      element(se)%delta_hbie(k)=0.22540333d0
                  end select
                else
                  element(se)%delta_hbie(k)=tmp_delta_hbie
                end if
              end do
            !
            ! DUAL_BM_SAME:
            ! All nodes with SBIE MCA + HBIE MCA (Burton & Miller).
            !
            case (5)
              ! Assign values
              node(sn)%sbie=fbem_sbie_mca
              node(sn)%hbie=fbem_hbie
              node(sn)%dual=fbem_dual_burton_miller
              node(sn)%dual_is_common=.true.
              node(sn)%alpha=1.0d0
              ! The displacement towards inside each element of the node for HBIE MCA formulation
              ! Loop through the elements of the node
              do ke=1,node(sn)%n_elements
                ! Selected element
                se=node(sn)%element(ke)
                ! Index of the node in the selected element
                k=node(sn)%element_node_iid(ke)
                ! Depending on tmp_delta_hbie, use the automatic value
                if (tmp_delta_hbie.le.0.0d0) then
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
                  end select
                else
                  element(se)%delta_sbie_mca(k)=tmp_delta_hbie
                  element(se)%delta_hbie(k)    =tmp_delta_hbie
                end if
              end do
            !
            ! DUAL_BM_MIXED:
            ! All nodes with SBIE + HBIE MCA (Burton & Miller), except nodes at the boundary, which have MCA collocation.
            !
            case (6)
              ! Assign values
              node(sn)%sbie=fbem_sbie
              node(sn)%hbie=fbem_hbie
              node(sn)%dual=fbem_dual_burton_miller
              node(sn)%dual_is_common=.false.
              ! The displacement towards inside each element of the node for HBIE formulation
              ! Loop through the elements of the node
              do ke=1,node(sn)%n_elements
                ! Selected element
                se=node(sn)%element(ke)
                ! Index of the node in the selected element
                k=node(sn)%element_node_iid(ke)
                ! Depending on tmp_delta_hbie, use the automatic value
                if (tmp_delta_hbie.le.0.0d0) then
                  ! It depends on the element type
                  select case (element(se)%type)
                    ! Linear elements
                    case (fbem_line2,fbem_tri3,fbem_quad4)
                      element(se)%delta_hbie(k)=0.42264973d0
                    ! Quadratic elements
                    case (fbem_line3,fbem_tri6,fbem_quad8,fbem_quad9)
                      element(se)%delta_hbie(k)=0.22540333d0
                  end select
                else
                  element(se)%delta_hbie(k)=tmp_delta_hbie
                end if
              end do
              ! If it is in the boundary of the boundary
              if (node(sn)%in_boundary) then
                ! Instead of having SBIE formulation, the node has SBIE MCA formulation
                node(sn)%sbie=fbem_sbie_mca
                ! The displacement towards inside each element of the node for SBIE MCA formulation
                ! Loop through the elements of the node
                do ke=1,node(sn)%n_elements
                  ! Selected element
                  se=node(sn)%element(ke)
                  ! Index of the node in the selected element
                  k=node(sn)%element_node_iid(ke)
                  ! Depending on tmp_delta_sbie, use the automatic value
                  if (tmp_delta_sbie.le.0.0d0) then
                    ! It depends on the element type
                    select case (element(se)%type)
                      ! Linear elements
                      case (fbem_line2,fbem_tri3,fbem_quad4)
                        element(se)%delta_sbie_mca(k)=0.42264973d0
                      ! Quadratic elements
                      case (fbem_line3,fbem_tri6,fbem_quad8,fbem_quad9)
                        element(se)%delta_sbie_mca(k)=0.22540333d0
                    end select
                  else
                    element(se)%delta_sbie_mca(k)=tmp_delta_sbie
                  end if
                end do
              end if
          end select
        !
        ! If the node belongs to a discontinuous element
        !
        else
          select case (tmp_type_value)
            !
            ! SBIE
            !
            case (1,2,3)
              ! Assign values
              node(sn)%sbie=fbem_sbie
              node(sn)%hbie=0
              node(sn)%dual=0
              node(sn)%dual_is_common=.false.
            !
            ! HBIE
            !
            case (4)
              ! Assign values
              node(sn)%sbie=0
              node(sn)%hbie=fbem_hbie
              node(sn)%dual=0
              node(sn)%dual_is_common=.false.
            !
            ! SBIE + HBIE (Burton & Miller)
            !
            case (5,6)
              ! Assign values
              node(sn)%sbie=fbem_sbie
              node(sn)%hbie=fbem_hbie
              node(sn)%dual=fbem_dual_burton_miller
              node(sn)%dual_is_common=.true.
              node(sn)%alpha=1.0d0
          end select
        end if
      end do

    ! ------------------------------------------------------------------------------------------------------------------------
    ! CRACK-LIKE
    ! ------------------------------------------------------------------------------------------------------------------------

    case (fbem_boundary_class_cracklike)
      ! Read the type of formulation
      read(input_fileunit,*) tmp_type
      tmp_type_value=0
      if (trim(tmp_type).eq.'same' ) tmp_type_value=1
      if (trim(tmp_type).eq.'mixed') tmp_type_value=2
      ! Check if valid
      if (tmp_type_value.eq.0) then
        call fbem_error_message(error_unit,0,'boundary',boundary(sb)%id,'the indicated BEM formulation is not valid.')
      end if
      ! Read parameters of the formulation if needed
      select case (tmp_type_value)
        ! Read delta_hbie (=delta_sbie)
        case (1)
          call fbem_search_section(input_fileunit,section_name,found)
          call fbem_search_keyword(input_fileunit,keyword,':',found)
          read(input_fileunit,*) tmp_type, tmp_delta_hbie
        ! Read delta_sbie y delta_hbie
        case (2)
          call fbem_search_section(input_fileunit,section_name,found)
          call fbem_search_keyword(input_fileunit,keyword,':',found)
          read(input_fileunit,*) tmp_type, tmp_delta_sbie, tmp_delta_hbie
      end select
      ! Loop through the nodes of the boundary
      do kn=1,part(boundary(sb)%part)%n_nodes
        ! Selected node
        sn=part(boundary(sb)%part)%node(kn)
        ! If the node belongs to continuous elements
        if (element(node(sn)%element(1))%discontinuous.eqv.(.false.)) then
          select case (tmp_type_value)
            !
            ! SAME:
            ! The node has SBIE MCA / HBIE (Dual BEM) formulation. SBIE and HBIE have the same collocation point.
            !
            case (1)
              ! Assign values
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
                ! Depending on tmp_delta_hbie, use the automatic value
                if (tmp_delta_hbie.le.0.0d0) then
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
                  end select
                else
                  element(se)%delta_sbie_mca(k)=tmp_delta_hbie
                  element(se)%delta_hbie(k)    =tmp_delta_hbie
                end if
              end do
            !
            ! MIXED:
            ! The node has SBIE / HBIE (Dual BEM) formulation. The SBIE is collocated at nodal positions, except at the
            ! boundaries of the boundary.
            !
            case (2)
              ! Assign values
              node(sn)%sbie=fbem_sbie
              node(sn)%hbie=fbem_hbie
              node(sn)%dual=fbem_dual_boundary
              node(sn)%dual_is_common=.false.
              ! The displacement towards inside each element of the node for HBIE formulation
              ! Loop through the elements of the node
              do ke=1,node(sn)%n_elements
                ! Selected element
                se=node(sn)%element(ke)
                ! Index of the node in the selected element
                k=node(sn)%element_node_iid(ke)
                ! Depending on tmp_delta_hbie, use the automatic value
                if (tmp_delta_hbie.le.0.0d0) then
                  ! It depends on the element type
                  select case (element(se)%type)
                    ! Linear elements
                    case (fbem_line2,fbem_tri3,fbem_quad4)
                      element(se)%delta_hbie(k)=0.42264973d0
                    ! Quadratic elements
                    case (fbem_line3,fbem_tri6,fbem_quad8,fbem_quad9)
                      element(se)%delta_hbie(k)=0.22540333d0
                  end select
                else
                  element(se)%delta_hbie(k)=tmp_delta_hbie
                end if
              end do
              ! If it is in the boundary of the boundary
              if (node(sn)%in_boundary) then
                ! Instead of having SBIE formulation, the node has SBIE MCA formulation
                node(sn)%sbie=fbem_sbie_mca
                ! The displacement towards inside each element of the node for SBIE MCA formulation
                ! Loop through the elements of the node
                do ke=1,node(sn)%n_elements
                  ! Selected element
                  se=node(sn)%element(ke)
                  ! Index of the node in the selected element
                  k=node(sn)%element_node_iid(ke)
                  ! Depending on tmp_delta_sbie, use the automatic value
                  if (tmp_delta_sbie.le.0.0d0) then
                    ! It depends on the element type
                    select case (element(se)%type)
                      ! Linear elements
                      case (fbem_line2,fbem_tri3,fbem_quad4)
                        element(se)%delta_sbie_mca(k)=0.42264973d0
                      ! Quadratic elements
                      case (fbem_line3,fbem_tri6,fbem_quad8,fbem_quad9)
                        element(se)%delta_sbie_mca(k)=0.22540333d0
                    end select
                  else
                    element(se)%delta_sbie_mca(k)=tmp_delta_sbie
                  end if
                end do
              end if
          end select
        !
        ! If the node belongs to a discontinuous element
        !
        else
          ! The node has SBIE / HBIE (Dual BEM) formulation.
          ! Assign values
          node(sn)%sbie=fbem_sbie
          node(sn)%hbie=fbem_hbie
          node(sn)%dual=fbem_dual_boundary
          node(sn)%dual_is_common=.true.
        end if
      end do
  end select

end subroutine read_bem_formulation_selected_nodes
