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


subroutine read_discontinuous_boundary_elements(input_fileunit)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_symmetry
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! No implicit variables are allowed
  implicit none

  ! I/O
  integer                                 :: input_fileunit              !! Input file unit
  ! Local
  character(len=fbem_stdcharlen)          :: section_name                ! Name of the section
  integer                                 :: n_discontinuous_boundaries
  integer                                 :: n_discontinuous_elements
  logical                                 :: found                       ! Logical variable for sections and keywords
  integer                                 :: kb, sb_eid, sb, sp
  integer                                 :: ke, se_eid, se
  integer                                 :: kn, sn, kej, sej, kek
  integer                                 :: newnode_eid
  integer                                 :: n_newnodes_max
  integer                                 :: n_newnodes
  integer, allocatable                    :: newnode_id(:)
  real(kind=real64), allocatable          :: newnode_x(:,:)
  integer, allocatable                    :: node_element_tmp(:)
  real(kind=real64)                       :: delta_f
  character(len=fbem_fmtstr)              :: fmt1, fmt2, fmt3            ! String used for write format string
  integer                                 :: output_fileunit
  character(len=fbem_filename_max_length) :: tmp_filename                ! Temporary file name
  logical                                 :: exist_discontinuous         ! Flag

  ! Return if not needed
  if (n_be_regions.eq.0) return





  ! Esto mismo se podria hacer para las body loads, pero por el momento solo continuos, habra que modificar cosas si lo metemos!!!



  ! Initialize
  exist_discontinuous=.false.

  ! ================================
  ! READ DISCONTINUOUS BE BOUNDARIES
  ! ================================

  ! Locate the section
  section_name='discontinuous be boundaries'
  call fbem_search_section(input_fileunit,section_name,found)
  if (found) then

    stop 'ERROR: discontinuous elements disabled (not completely implemented)'

    ! Starting message
    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')
    ! Read the number of boundaries where their elements are discontinuous
    read(input_fileunit,*) n_discontinuous_boundaries
    ! There exist discontinuous BE
    if (n_discontinuous_boundaries.gt.0) exist_discontinuous=.true.
    ! Read each boundary and its delta parameter
    do kb=1,n_discontinuous_boundaries
      read(input_fileunit,*) sb_eid, delta_f
      ! Check if the indicated boundary exists
      if ((sb_eid.ge.boundary_eid_min).and.(sb_eid.le.boundary_eid_max)) then
        if (boundary_iid(sb_eid).eq.0) then
          call fbem_warning_message(error_unit,0,section_name,0,'wrong definition of a boundary in this section')
          call fbem_error_message(error_unit,0,'boundary',sb_eid,'does not exist')
        end if
      else
        call fbem_warning_message(error_unit,0,section_name,0,'wrong definition of a boundary in this section')
        call fbem_error_message(error_unit,0,'boundary',sb_eid,'does not exist')
      end if
      sb=boundary_iid(sb_eid)
      sp=boundary(sb)%part
      ! If delta_f<=0. then set the default value
      if (delta_f.le.check_xi_tolerance) then
        do ke=1,part(sp)%n_elements
          element(part(sp)%element(ke))%discontinuous=.true.
          ! It depends on the element type
          select case (element(part(sp)%element(ke))%type)
            ! Linear elements
            case (fbem_line2,fbem_tri3,fbem_quad4)
              element(part(sp)%element(ke))%delta_f=0.42264973d0
            ! Quadratic elements
            case (fbem_line3,fbem_tri6,fbem_quad8,fbem_quad9)
              element(part(sp)%element(ke))%delta_f=0.22540333d0
          end select
        end do
      ! If not, set the indicated value
      else
        do ke=1,part(sp)%n_elements
          element(part(sp)%element(ke))%discontinuous=.true.
          ! It depends on the element type
          select case (element(part(sp)%element(ke))%type)
            ! Linear elements
            case (fbem_line2,fbem_tri3,fbem_quad4)
              element(part(sp)%element(ke))%delta_f=delta_f
            ! Quadratic elements
            case (fbem_line3,fbem_tri6,fbem_quad8,fbem_quad9)
              element(part(sp)%element(ke))%delta_f=delta_f
          end select
        end do
      end if
    end do
    ! Ending message
    if (verbose_level.ge.2) write(output_unit,'(1x,a)') 'done.'
  end if

  ! ==============================
  ! READ DISCONTINUOUS BE ELEMENTS
  ! ==============================

  ! Locate the section
  section_name='discontinuous be elements'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(input_fileunit,section_name,found)
  if (found) then

    stop 'ERROR: discontinuous elements disabled (not completely implemented)'

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! Read the number of discontinuous elements to be set
    read(input_fileunit,*) n_discontinuous_elements
    ! There exist discontinuous BE
    if (n_discontinuous_elements.gt.0) exist_discontinuous=.true.
    ! Read each element and its delta parameter
    do ke=1,n_discontinuous_elements
      read(input_fileunit,*) se_eid, delta_f
      se=element_iid(se_eid)
      if (part(element(se)%part)%type.eq.fbem_part_be_boundary) then
        element(se)%discontinuous=.true.
        ! If delta_f<=0. then set the default value
        if (delta_f.le.check_xi_tolerance) then
          ! It depends on the element type
          select case (element(se)%type)
            ! Linear elements
            case (fbem_line2,fbem_tri3,fbem_quad4)
              element(se)%delta_f=0.42264973d0
            ! Quadratic elements
            case (fbem_line3,fbem_tri6,fbem_quad8,fbem_quad9)
              element(se)%delta_f=0.22540333d0
          end select
        else
          ! It depends on the element type
          select case (element(se)%type)
            ! Linear elements
            case (fbem_line2,fbem_tri3,fbem_quad4)
              element(se)%delta_f=delta_f
            ! Quadratic elements
            case (fbem_line3,fbem_tri6,fbem_quad8,fbem_quad9)
              element(se)%delta_f=delta_f
          end select
        end if
      else
        call fbem_error_message(error_unit,0,'element',element(se)%id,'is not a element belonging to a BE boundary')
      end if
    end do
    ! Ending message
    if (verbose_level.ge.2) write(output_unit,'(1x,a)') 'done.'
  end if

  ! =======================================================
  ! CHECK THAT EACH DISCONTINUOUS ELEMENT HAVE UNIQUE NODES
  ! =======================================================

  if (exist_discontinuous) then
    ! Starting message
    if (verbose_level.ge.2) write(output_unit,'(1x,a)') 'Checking that discontinuous BE elements have unique nodes ...'
    !
    ! Calculate the maximum number of new nodes, allocate the auxiliary data structures, and build them.
    !
    n_newnodes_max=0
    do ke=1,n_elements
      n_newnodes_max=n_newnodes_max+element(ke)%n_nodes
    end do
    allocate(newnode_id(n_newnodes_max),newnode_x(problem%n,n_newnodes_max))
    do kn=1,n_nodes
      newnode_id(kn)=node(kn)%id
      newnode_x(:,kn)=node(kn)%x
    end do
    !
    ! Generate new mesh if needed
    !
    ! Existent maximum internal and external identifiers
    newnode_eid=node_eid_max
    n_newnodes=n_nodes
    ! Loop through elements
    do ke=1,n_elements
      ! If an element belonging to a BE boundary and discontinuous.
      if ((part(element(ke)%part)%type.eq.fbem_part_be_boundary).and.(element(ke)%discontinuous.eqv.(.true.))) then
        ! Loop through the nodes of the discontinuous element
        do kn=1,element(ke)%n_nodes
          ! Selected node
          sn=element(ke)%node(kn)
          ! If the node is shared by another element, create a new one for the present element, and change the connectivities.
          if (node(sn)%n_elements.gt.1) then
            ! Create new new node
            newnode_eid=newnode_eid+1
            n_newnodes=n_newnodes+1
            newnode_id(n_newnodes)=newnode_eid
            newnode_x(:,n_newnodes)=node(sn)%x
            ! Connect the new node to the element
            element(ke)%node(kn)=n_newnodes
            ! Disconnect the element from the existent node, and copy the new connectivity to a temporary variable
            allocate(node_element_tmp(node(sn)%n_elements-1))
            kek=1
            do kej=1,node(sn)%n_elements
              sej=node(sn)%element(kej)
              if (sej.ne.ke) then
                node_element_tmp(kek)=sej
                kek=kek+1
              end if
            end do
            ! Build the new connectivity vector
            node(sn)%n_elements=node(sn)%n_elements-1
            deallocate(node(sn)%element)
            allocate(node(sn)%element(node(sn)%n_elements))
            node(sn)%element=node_element_tmp
            deallocate(node_element_tmp)
          end if
        end do
      end if
    end do
    !
    ! If new nodes have been created, then export the new mesh
    !
    if (n_newnodes.gt.n_nodes) then
      !
      ! Open ".nmh" file
      !
      output_fileunit=fbem_get_valid_unit()
      write(tmp_filename,*) trim(output_filename), '.nmh'
      call fbem_trim2b(tmp_filename)
      open(unit=output_fileunit,file=trim(tmp_filename),action='write',recl=fbem_file_record_length)
      write(output_unit,'(a,a,a)') ' Opening "', trim(tmp_filename), '"'
      !
      ! Section [nodes]
      !
      write(output_fileunit,'(a)') '[nodes]'
      write(output_fileunit,'(i11)') n_newnodes
      select case (problem%n)
        case (2)
           do kn=1,n_newnodes
            write(output_fileunit,'(i11,2e25.16)') newnode_id(kn), newnode_x(:,kn)
          end do
        case (3)
          do kn=1,n_newnodes
            write(output_fileunit,'(i11,3e25.16)') newnode_id(kn), newnode_x(:,kn)
          end do
      end select
      !
      ! Section [elements]
      !
      write(output_fileunit,'(a)') '[elements]'
      write(output_fileunit,'(i11)') n_elements
      fmt1='(i11)'
      fmt2='(a8)'
      do ke=1,n_elements
        write(output_fileunit,fmt1,advance='no') element(ke)%id
        select case (element(ke)%type)
          case (fbem_line2)
            write(output_fileunit,fmt2,advance='no') 'line2'
          case (fbem_line3)
            write(output_fileunit,fmt2,advance='no') 'line3'
          case (fbem_tri3)
            write(output_fileunit,fmt2,advance='no') 'tri3'
          case (fbem_tri6)
            write(output_fileunit,fmt2,advance='no') 'tri6'
          case (fbem_quad4)
            write(output_fileunit,fmt2,advance='no') 'quad4'
          case (fbem_quad8)
            write(output_fileunit,fmt2,advance='no') 'quad8'
          case (fbem_quad9)
            write(output_fileunit,fmt2,advance='no') 'quad9'
        end select
        write(output_fileunit,'(i3,i11)',advance='no') 1, element(ke)%part
        write(fmt3,*) '(',element(ke)%n_nodes,'i11)'
        call fbem_trim2b(fmt3)
        write(output_fileunit,fmt3) (newnode_id(element(ke)%node(kn)),kn=1,element(ke)%n_nodes)
      end do
      !
      ! Close file
      !
      close(unit=output_fileunit)
      write(output_unit,'(a,a,a)') ' Closing "', trim(tmp_filename), '"'
      stop 'Note: replace [nodes] and [elements] sections by the new ones saved in the *.nmh file above'
    end if

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  end if

end subroutine read_discontinuous_boundary_elements
