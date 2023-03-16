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
  real(kind=real64), parameter          :: default_mca_boundary_delta  = 0.05d0
  real(kind=real64), parameter          :: default_mca_linear_delta    = 0.42264973d0
  real(kind=real64), parameter          :: default_mca_quadratic_delta = 0.22540333d0
  real(kind=real64), parameter          :: default_mca_cubic_delta     = 0.138863688d0
  logical                               :: found
  integer                               :: k
  integer                               :: sp, sb, ke, se, kn, sn
  character(len=fbem_stdcharlen)        :: tmp_type
  real(kind=real64)                     :: tmp_delta_sbie, tmp_delta_hbie
  integer                               :: tmp_type_value


  ! Nota: La dual B&M poroelastica no esta bien!!, la suma hay que ponderarla con otro numero de onda o hacerla de otra
  ! manera...

  ! Part and entity (BE boundary or BE bodyload) of the selected nodes
  sp=node(selected_node(1))%part(1)
  sb=part(sp)%entity

  select case (part(sp)%type)

    !
    ! BE BOUNDARY
    !
    case (fbem_part_be_boundary)

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
                  node(sn)%sbie=fbem_sbie
                  node(sn)%hbie=0
                  node(sn)%dual=0
                  node(sn)%dual_is_common=.false.
                  if (node(sn)%in_boundary) then
                    node(sn)%sbie=fbem_sbie_mca
                    do ke=1,node(sn)%n_elements
                      se=node(sn)%element(ke)
                      k=node(sn)%element_node_iid(ke)
                      if (tmp_delta_sbie.le.0.0d0) then
                        element(se)%delta_sbie_mca(k)=default_mca_boundary_delta
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
                  node(sn)%sbie=fbem_sbie_mca
                  node(sn)%hbie=0
                  node(sn)%dual=0
                  node(sn)%dual_is_common=.false.
                  do ke=1,node(sn)%n_elements
                    se=node(sn)%element(ke)
                    k=node(sn)%element_node_iid(ke)
                    if (tmp_delta_sbie.le.0.0d0) then
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
                        case default
                          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid element order')
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
                  node(sn)%sbie=0
                  node(sn)%hbie=fbem_hbie
                  node(sn)%dual=0
                  node(sn)%dual_is_common=.false.
                  do ke=1,node(sn)%n_elements
                    se=node(sn)%element(ke)
                    k=node(sn)%element_node_iid(ke)
                    if (tmp_delta_hbie.le.0.0d0) then
                      select case (fbem_element_order(element(se)%type))
                        ! Linear elements
                        case (1)
                          element(se)%delta_hbie(k)=default_mca_linear_delta
                        ! Quadratic elements
                        case (2)
                          element(se)%delta_hbie(k)=default_mca_quadratic_delta
                        ! Cubic elements
                        case (3)
                          element(se)%delta_hbie(k)=default_mca_cubic_delta
                        case default
                          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid element order')
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
                  node(sn)%sbie=fbem_sbie_mca
                  node(sn)%hbie=fbem_hbie
                  node(sn)%dual=fbem_dual_burton_miller
                  node(sn)%dual_is_common=.true.
                  node(sn)%alpha=1.0d0
                  do ke=1,node(sn)%n_elements
                    se=node(sn)%element(ke)
                    k=node(sn)%element_node_iid(ke)
                    if (tmp_delta_hbie.le.0.0d0) then
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
                  node(sn)%sbie=fbem_sbie
                  node(sn)%hbie=fbem_hbie
                  node(sn)%dual=fbem_dual_burton_miller
                  node(sn)%dual_is_common=.false.
                  do ke=1,node(sn)%n_elements
                    se=node(sn)%element(ke)
                    k=node(sn)%element_node_iid(ke)
                    if (tmp_delta_hbie.le.0.0d0) then
                      select case (fbem_element_order(element(se)%type))
                        ! Linear elements
                        case (1)
                          element(se)%delta_hbie(k)=default_mca_linear_delta
                        ! Quadratic elements
                        case (2)
                          element(se)%delta_hbie(k)=default_mca_quadratic_delta
                        ! Cubic elements
                        case (3)
                          element(se)%delta_hbie(k)=default_mca_cubic_delta
                        case default
                          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid element order')
                      end select
                    else
                      element(se)%delta_hbie(k)=tmp_delta_hbie
                    end if
                  end do
                  if (node(sn)%in_boundary) then
                    node(sn)%sbie=fbem_sbie_mca
                    do ke=1,node(sn)%n_elements
                      se=node(sn)%element(ke)
                      k=node(sn)%element_node_iid(ke)
                      if (tmp_delta_sbie.le.0.0d0) then
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
                          case default
                            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid element order')
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
                  node(sn)%sbie=fbem_sbie_mca
                  node(sn)%hbie=fbem_hbie
                  node(sn)%dual=fbem_dual_boundary
                  node(sn)%dual_is_common=.true.
                  do ke=1,node(sn)%n_elements
                    se=node(sn)%element(ke)
                    k=node(sn)%element_node_iid(ke)
                    if (tmp_delta_hbie.le.0.0d0) then
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
                  node(sn)%sbie=fbem_sbie
                  node(sn)%hbie=fbem_hbie
                  node(sn)%dual=fbem_dual_boundary
                  node(sn)%dual_is_common=.false.
                  do ke=1,node(sn)%n_elements
                    se=node(sn)%element(ke)
                    k=node(sn)%element_node_iid(ke)
                    if (tmp_delta_hbie.le.0.0d0) then
                      select case (fbem_element_order(element(se)%type))
                        ! Linear elements
                        case (1)
                          element(se)%delta_hbie(k)=default_mca_linear_delta
                        ! Quadratic elements
                        case (2)
                          element(se)%delta_hbie(k)=default_mca_quadratic_delta
                        ! Cubic elements
                        case (3)
                          element(se)%delta_hbie(k)=default_mca_cubic_delta
                        case default
                          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid element order')
                      end select
                    else
                      element(se)%delta_hbie(k)=tmp_delta_hbie
                    end if
                  end do
                  if (node(sn)%in_boundary) then
                    node(sn)%sbie=fbem_sbie_mca
                    do ke=1,node(sn)%n_elements
                      se=node(sn)%element(ke)
                      k=node(sn)%element_node_iid(ke)
                      if (tmp_delta_sbie.le.0.0d0) then
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
                          case default
                            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid element order')
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

    !
    ! BE BODYLOAD
    !
    case (fbem_part_be_bodyload)
      ! Read the type of formulation
      read(input_fileunit,*) tmp_type
      tmp_type_value=0
      if (trim(tmp_type).eq.'sbie'                      ) tmp_type_value=1
      if (trim(tmp_type).eq.'sbie_mca'                  ) tmp_type_value=2
      if (trim(tmp_type).eq.'sbie_lineload_end_boundary') tmp_type_value=3
      ! Check if valid
      if (tmp_type_value.eq.0) then
        call fbem_error_message(error_unit,0,'be_bodyload',be_bodyload(sb)%id,'the indicated BEM formulation does not exist.')
      end if
      ! Read parameters of the formulation if needed
      if (tmp_type_value.eq.2) then
        call fbem_search_section(input_fileunit,section_name,found)
        call fbem_search_keyword(input_fileunit,keyword,':',found)
        read(input_fileunit,*) tmp_type, tmp_delta_sbie
      end if
      ! Loop through the nodes of the boundary
      do kn=1,selected_n_nodes
        ! Selected node
        sn=selected_node(kn)
        select case (tmp_type_value)
          !
          ! SBIE:
          ! All nodes with nodal collocation. This is unsafe since at doubled nodes, depending on the boundary conditions, it
          ! will lead to a singular linear system of equations.
          !
          case (1,3)
            node(sn)%sbie=fbem_sbie
            node(sn)%hbie=0
            node(sn)%dual=0
            node(sn)%dual_is_common=.false.
            if (tmp_type_value.eq.3) node(sn)%sbie_lineload_end_boundary=.true.
          !
          ! SBIE_MCA:
          ! All nodes with SBIE MCA collocation.
          !
          case (2)
            node(sn)%sbie=fbem_sbie_mca
            node(sn)%hbie=0
            node(sn)%dual=0
            node(sn)%dual_is_common=.false.
            do ke=1,node(sn)%n_elements
              se=node(sn)%element(ke)
              k=node(sn)%element_node_iid(ke)
              if (tmp_delta_sbie.le.0.0d0) then
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
                  case default
                    call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid element order')
                end select
              else
                element(se)%delta_sbie_mca(k)=tmp_delta_sbie
              end if
            end do
        end select
      end do

  end select


end subroutine read_bem_formulation_selected_nodes
