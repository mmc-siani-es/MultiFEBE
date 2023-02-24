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


!! @author Jacob David Rodriguez Bordon (jacobdavid.rodriguezbordon@ulpgc.es)
!!
!! @version 2.0
!!
!! <b> Subroutine that reads the elements from a file. </b>
!!
!! Section format (similar to gmsh):
!!>
!! [elements]
!! <number of elements>
!! ...
!! <element id> <type> <number of tags (>0)> <tag 1 (part)> <node 1> <node 2> ... <node N>
!! ...
!!<

subroutine read_elements(fileunit,mode)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  implicit none
  ! I/O
  integer                        :: fileunit          !! Unit of the file to read from
  integer                        :: mode              !! Mode of reading: 1 (native format), 2 (gmsh format)
  ! Local
  character(len=fbem_stdcharlen) :: section_name      ! Name of the section
  logical                        :: found
  integer                        :: tmp_id            ! Temporary variable to read the id of the element
  integer                        :: tmp_n_tags        ! Temporary variable to read the number of tags of the element
  integer                        :: tmp_tags(8)       ! Temporary variable to read the number of tags of the element
  integer                        :: tmp_part          ! Temporary variable to read the part of the element
  character(len=fbem_stdcharlen) :: tmp_type          ! Temporary variable to read the type of the element
  integer                        :: tmp_n_elements    ! Temporary variable to read the number of elements in the file
  integer                        :: i, j, k, ke       ! Counters

  ! Locate the section
  select case (mode)
    case (1)
      section_name='elements'
      call fbem_search_section(fileunit,section_name,found)
    case (2)
      section_name='Elements'
      call fbem_search_section_gmsh(fileunit,section_name,found)
  end select
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! ===============================================================
    ! CALCULATE THE NUMBER OF ELEMENTS THAT BELONGS TO AN ACTIVE PART
    ! ===============================================================

    n_elements=0
    read(fileunit,*) tmp_n_elements
    do i=1,tmp_n_elements
      ! Read
      read(fileunit,*) tmp_id, tmp_type, tmp_n_tags, tmp_part
      ! Check number of tags
      if (tmp_n_tags.gt.8) then
        call fbem_error_message(error_unit,0,'element',tmp_id,'has >8 integer tags')
      end if
      ! Ignore if a part id is 0 for gmsh file format
      if ((tmp_part.eq.0).and.(mode.eq.2)) then
        ! Nothing to do
      else
        ! Check if the indicated part exists
        if ((tmp_part.ge.part_eid_min).and.(tmp_part.le.part_eid_max)) then
          if (part_iid(tmp_part).eq.0) then
            call fbem_warning_message(error_unit,0,'element',tmp_id,'wrong definition of its part')
            call fbem_error_message(error_unit,0,'part',tmp_part,'does not exist')
          end if
        else
          call fbem_warning_message(error_unit,0,'element',tmp_id,'wrong definition of its part')
          call fbem_error_message(error_unit,0,'part',tmp_part,'does not exist')
        end if

        !
        ! Include only physically connected elements
        !
        if (part(part_iid(tmp_part))%entity.gt.0) then
          n_elements=n_elements+1
        end if



!         !
!         ! Do not include physically unconnected elements (to be completely implemented in the future)
!         !
!         n_elements=n_elements+1



      end if
    end do

    ! Allocate and initialize
    if (n_elements.gt.0) then
      allocate (element(n_elements))
      do i=1,n_elements
        element(i)%id=0
        element(i)%type=0
        element(i)%n_dimension=0
        element(i)%size=0.d0
        element(i)%part=0
        element(i)%n_nodes=0
        element(i)%element=0
        element(i)%n_symplanes=0
        element(i)%dm_n_elements=0
        element(i)%dm_mode=0
        element(i)%plane=.false.
        element(i)%bball_radius=0.d0
        element(i)%n_subedges=0
        element(i)%n_subfaces=0
        element(i)%n_supelements=0
        element(i)%type_g=0
        element(i)%type_f1=0
        element(i)%type_f2=0
        element(i)%type_f1f=0
        element(i)%type_f2f=0
        element(i)%discontinuous=.false.
        element(i)%delta_f=0.d0
        element(i)%csize=0.d0
        element(i)%n_phi=0
        element(i)%n_dof=0
        element(i)%fe_type=0
        element(i)%fe_options=0
        element(i)%K_intmode=0
        element(i)%K_intngp=0
        element(i)%M_intmode=0
        element(i)%M_intngp=0
        element(i)%Q_intmode=0
        element(i)%Q_intngp=0
        element(i)%k_r=0
        element(i)%k_c=0
        element(i)%c=0
        element(i)%em_U=0
        element(i)%em_Bl=0
        element(i)%em_R=0
        element(i)%em_L=0
        element(i)%flooded=.false.
        element(i)%submerged=.false.
        element(i)%cs_type=0
        element(i)%cs_param=0
        element(i)%A=0
        element(i)%I=0
        element(i)%ksh=0
        element(i)%r_integration=0
        element(i)%orthotropic_shell_fd1=0
        element(i)%export=.true.
      end do
    else
      call fbem_error_message(error_unit,0,section_name,0,'the number of elements must be >0')
    end if

    ! =============================
    ! READ ELEMENTS IN ACTIVE PARTS
    ! =============================

    select case (mode)
      case (1)
        call fbem_search_section(fileunit,section_name,found)
      case (2)
        call fbem_search_section_gmsh(fileunit,section_name,found)
    end select
    ke=0
    read(fileunit,*) tmp_n_elements
    do i=1,tmp_n_elements
      read(fileunit,*) tmp_id, tmp_type, tmp_n_tags, (tmp_tags(j),j=1,tmp_n_tags)
      ! Ignore if a part id is 0 for gmsh file format
      if ((tmp_part.eq.0).and.(mode.eq.2)) then
        ! Nothing to do
        cycle
      else
        ! Check id
        if (tmp_id.le.0) then
          call fbem_error_message(error_unit,0,'element',tmp_id,'identifiers must be greater than 0')
        end if
        !
        ! Include only physically connected elements
        !
        if (part(part_iid(tmp_tags(1)))%entity.gt.0) then
          ke=ke+1
        else
          cycle
        end if
!       !
!       ! Do not include physically unconnected elements (to be completely implemented in the future)
!       !
!		ke=ke+1

      end if
      ! Save to the data structure the id and the part
      element(ke)%id=tmp_id
      element(ke)%part=part_iid(tmp_tags(1))
      ! Back to the line in order to read the nodes of the element
      backspace(fileunit)
      !
      ! Switch between the type of mesh partition (BE mesh or FE mesh)
      !
      select case (part(element(ke)%part)%type)

        ! ==========================================================================================================================
        ! BOUNDARY ELEMENT
        ! ==========================================================================================================================

        case(fbem_part_be_boundary)

          element(ke)%type=0
          select case (problem%n)
            case (2)
              if ((trim(tmp_type).eq.'line2').or.(trim(tmp_type).eq.'1' )) element(ke)%type=fbem_line2
              if ((trim(tmp_type).eq.'line3').or.(trim(tmp_type).eq.'8' )) element(ke)%type=fbem_line3
              if ((trim(tmp_type).eq.'line4').or.(trim(tmp_type).eq.'26')) element(ke)%type=fbem_line4
            case (3)
              if ((trim(tmp_type).eq.'tri3' ).or.(trim(tmp_type).eq.'2' )) element(ke)%type=fbem_tri3
              if ((trim(tmp_type).eq.'tri6' ).or.(trim(tmp_type).eq.'9' )) element(ke)%type=fbem_tri6
              if ((trim(tmp_type).eq.'quad4').or.(trim(tmp_type).eq.'3' )) element(ke)%type=fbem_quad4
              if ((trim(tmp_type).eq.'quad8').or.(trim(tmp_type).eq.'16')) element(ke)%type=fbem_quad8
              if ((trim(tmp_type).eq.'quad9').or.(trim(tmp_type).eq.'10')) element(ke)%type=fbem_quad9
          end select
          if (element(ke)%type.eq.0) then
            call fbem_error_message(error_unit,0,'element',element(ke)%id,'the indicated type of boundary element does not exist')
          end if
          element(ke)%type_g=element(ke)%type
          element(ke)%type_f1=element(ke)%type
          element(ke)%type_f2=element(ke)%type
          element(ke)%type_f1f=element(ke)%type
          element(ke)%type_f2f=element(ke)%type
          element(ke)%discontinuous=.false.
          element(ke)%delta_f=0

        ! ==========================================================================================================================
        ! FINITE ELEMENT
        ! ==========================================================================================================================

        case(fbem_part_fe_subregion,0)

          element(ke)%type=0
          select case (problem%n)
            case (2)
              if ((trim(tmp_type).eq.'point1').or.(trim(tmp_type).eq.'15')) element(ke)%type=fbem_point1
              if ((trim(tmp_type).eq.'line2' ).or.(trim(tmp_type).eq.'1' )) element(ke)%type=fbem_line2
              if ((trim(tmp_type).eq.'line3' ).or.(trim(tmp_type).eq.'8' )) element(ke)%type=fbem_line3
              if ((trim(tmp_type).eq.'line4' ).or.(trim(tmp_type).eq.'26')) element(ke)%type=fbem_line4
              if ((trim(tmp_type).eq.'tri3'  ).or.(trim(tmp_type).eq.'2' )) element(ke)%type=fbem_tri3
              if ((trim(tmp_type).eq.'tri6'  ).or.(trim(tmp_type).eq.'9' )) element(ke)%type=fbem_tri6
              if ((trim(tmp_type).eq.'quad4' ).or.(trim(tmp_type).eq.'3' )) element(ke)%type=fbem_quad4
              if ((trim(tmp_type).eq.'quad8' ).or.(trim(tmp_type).eq.'16')) element(ke)%type=fbem_quad8
              if ((trim(tmp_type).eq.'quad9' ).or.(trim(tmp_type).eq.'10')) element(ke)%type=fbem_quad9
            case (3)
              if ((trim(tmp_type).eq.'point1').or.(trim(tmp_type).eq.'15')) element(ke)%type=fbem_point1
              if ((trim(tmp_type).eq.'line2' ).or.(trim(tmp_type).eq.'1' )) element(ke)%type=fbem_line2
              if ((trim(tmp_type).eq.'line3' ).or.(trim(tmp_type).eq.'8' )) element(ke)%type=fbem_line3
              if ((trim(tmp_type).eq.'line4' ).or.(trim(tmp_type).eq.'26')) element(ke)%type=fbem_line4
              if ((trim(tmp_type).eq.'tri3'  ).or.(trim(tmp_type).eq.'2' )) element(ke)%type=fbem_tri3
              if ((trim(tmp_type).eq.'tri6'  ).or.(trim(tmp_type).eq.'9' )) element(ke)%type=fbem_tri6
              if ((trim(tmp_type).eq.'quad4' ).or.(trim(tmp_type).eq.'3' )) element(ke)%type=fbem_quad4
              if ((trim(tmp_type).eq.'quad8' ).or.(trim(tmp_type).eq.'16')) element(ke)%type=fbem_quad8
              if ((trim(tmp_type).eq.'quad9' ).or.(trim(tmp_type).eq.'10')) element(ke)%type=fbem_quad9
              !if ((trim(tmp_type).eq.'tet4' ).or.(trim(tmp_type).eq.'4' )) element(ke)%type=fbem_tet4
              !if ((trim(tmp_type).eq.'tet10').or.(trim(tmp_type).eq.'11')) element(ke)%type=fbem_tet10
              !if ((trim(tmp_type).eq.'hex8' ).or.(trim(tmp_type).eq.'5' )) element(ke)%type=fbem_hex8
              !if ((trim(tmp_type).eq.'hex20').or.(trim(tmp_type).eq.'17')) element(ke)%type=fbem_hex20
              !if ((trim(tmp_type).eq.'hex27').or.(trim(tmp_type).eq.'12')) element(ke)%type=fbem_hex27
          end select
          if (element(ke)%type.eq.0) then
            call fbem_error_message(error_unit,0,'element',element(ke)%id,'the indicated type of finite element does not exist')
          end if

        ! ==========================================================================================================================
        ! BEM BODY LOAD ELEMENT
        ! ==========================================================================================================================

        case(fbem_part_be_bodyload)

          element(ke)%type=0
          select case (problem%n)
            case (2)
              if ((trim(tmp_type).eq.'point1').or.(trim(tmp_type).eq.'15')) element(ke)%type=fbem_point1
              if ((trim(tmp_type).eq.'line2' ).or.(trim(tmp_type).eq.'1' )) element(ke)%type=fbem_line2
              if ((trim(tmp_type).eq.'line3' ).or.(trim(tmp_type).eq.'8' )) element(ke)%type=fbem_line3
              if ((trim(tmp_type).eq.'line4' ).or.(trim(tmp_type).eq.'26')) element(ke)%type=fbem_line4
              if ((trim(tmp_type).eq.'tri3'  ).or.(trim(tmp_type).eq.'2' )) element(ke)%type=fbem_tri3
              if ((trim(tmp_type).eq.'tri6'  ).or.(trim(tmp_type).eq.'9' )) element(ke)%type=fbem_tri6
              if ((trim(tmp_type).eq.'quad4' ).or.(trim(tmp_type).eq.'3' )) element(ke)%type=fbem_quad4
              if ((trim(tmp_type).eq.'quad8' ).or.(trim(tmp_type).eq.'16')) element(ke)%type=fbem_quad8
              if ((trim(tmp_type).eq.'quad9' ).or.(trim(tmp_type).eq.'10')) element(ke)%type=fbem_quad9
            case (3)
              if ((trim(tmp_type).eq.'point1').or.(trim(tmp_type).eq.'15')) element(ke)%type=fbem_point1
              if ((trim(tmp_type).eq.'line2' ).or.(trim(tmp_type).eq.'1' )) element(ke)%type=fbem_line2
              if ((trim(tmp_type).eq.'line3' ).or.(trim(tmp_type).eq.'8' )) element(ke)%type=fbem_line3
              if ((trim(tmp_type).eq.'line4' ).or.(trim(tmp_type).eq.'26')) element(ke)%type=fbem_line4
              if ((trim(tmp_type).eq.'tri3'  ).or.(trim(tmp_type).eq.'2' )) element(ke)%type=fbem_tri3
              if ((trim(tmp_type).eq.'tri6'  ).or.(trim(tmp_type).eq.'9' )) element(ke)%type=fbem_tri6
              if ((trim(tmp_type).eq.'quad4' ).or.(trim(tmp_type).eq.'3' )) element(ke)%type=fbem_quad4
              if ((trim(tmp_type).eq.'quad8' ).or.(trim(tmp_type).eq.'16')) element(ke)%type=fbem_quad8
              if ((trim(tmp_type).eq.'quad9' ).or.(trim(tmp_type).eq.'10')) element(ke)%type=fbem_quad9
              !if ((trim(tmp_type).eq.'tet4' ).or.(trim(tmp_type).eq.'4' )) element(ke)%type=fbem_tet4
              !if ((trim(tmp_type).eq.'tet10').or.(trim(tmp_type).eq.'11')) element(ke)%type=fbem_tet10
              !if ((trim(tmp_type).eq.'hex8' ).or.(trim(tmp_type).eq.'5' )) element(ke)%type=fbem_hex8
              !if ((trim(tmp_type).eq.'hex20').or.(trim(tmp_type).eq.'17')) element(ke)%type=fbem_hex20
              !if ((trim(tmp_type).eq.'hex27').or.(trim(tmp_type).eq.'12')) element(ke)%type=fbem_hex27
          end select
          if (element(ke)%type.eq.0) then
            call fbem_error_message(error_unit,0,'element',element(ke)%id,'the indicated BE body load element type does not exist')
          end if
          element(ke)%type_g=element(ke)%type
          element(ke)%type_f1=element(ke)%type
          element(ke)%type_f2=element(ke)%type
          element(ke)%type_f1f=element(ke)%type
          element(ke)%type_f2f=element(ke)%type
          element(ke)%discontinuous=.false.
          element(ke)%delta_f=0

      end select

      ! Comparando n_dimension con problem%n y sabiendo si BE or FE, se sabe si es un elemento estructural o no
      ! Si es estructural necesita la informacion acerca de la seccion (longitudes, areas, momentos de inercia, ...)

      ! Save the number of nodes of the element
      element(ke)%n_nodes=fbem_n_nodes(element(ke)%type)
      ! Save the number of dimensions of the element
      element(ke)%n_dimension=fbem_n_dimension(element(ke)%type)
      ! Allocate the vector of element nodes
      allocate(element(ke)%node(element(ke)%n_nodes))
      ! Read the nodes
      read(fileunit,*) tmp_id, tmp_type, tmp_n_tags, (tmp_tags(j),j=1,tmp_n_tags), (element(ke)%node(j),j=1,element(ke)%n_nodes)

    end do

    ! ===================================
    ! CHECK AND BUILD ELEMENT IDENTIFIERS
    ! ===================================

    element_eid_min=element(1)%id
    element_eid_max=element(1)%id
    do i=2,n_elements
      if (element(i)%id.lt.element_eid_min) element_eid_min=element(i)%id
      if (element(i)%id.gt.element_eid_max) element_eid_max=element(i)%id
    end do
    allocate (element_iid(element_eid_min:element_eid_max))
    element_iid=0
    do i=1,n_elements
      if (element_iid(element(i)%id).ne.0) then
        call fbem_error_message(error_unit,0,'element',element(i)%id,'is repeated')
      else
        element_iid(element(i)%id)=i
      end if
    end do

    ! ====================================
    ! BUILD CONNECTIVITY PART --> ELEMENTS
    ! ====================================

    do i=1,n_parts
      part(i)%n_elements=0
      do j=1,n_elements
        if (element(j)%part.eq.i) then
          part(i)%n_elements=part(i)%n_elements+1
        end if
      end do
      if (part(i)%n_elements.ne.0) then
        allocate (part(i)%element(part(i)%n_elements))
        k=0
        do j=1,n_elements
          if (element(j)%part.eq.i) then
            k=k+1
            part(i)%element(k)=j
          end if
        end do
      end if
    end do

    ! ============================================================
    ! CHECK THAT EACH PART CONTAINS ELEMENTS OF THE SAME DIMENSION
    ! ============================================================

    do i=1,n_parts
      do j=2,part(i)%n_elements
        if (element(part(i)%element(j))%n_dimension.ne.element(part(i)%element(1))%n_dimension) then
          call fbem_error_message(error_unit,0,'part',part(i)%id,'contain elements with diferent number of dimensions')
        end if
      end do
    end do

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  else
    call fbem_error_message(error_unit,0,trim(section_name),0,'this section is required')
  end if

end subroutine read_elements
