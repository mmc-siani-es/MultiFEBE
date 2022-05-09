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

!! @author Jacob David Rodriguez Bordon (jacobdavid.rodriguezbordon@ulpgc.es)
!!
!! @version 2.0
!!
!! <b> Subroutine that reads the nodes from a file. </b>
!!
!!>
!! <number of nodes>
!! ...
!! <node i id> <node i x(1)> <node i x(2)> [<node i x(3)>]
!! ...
!!<
subroutine read_nodes(fileunit,mode)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! Local variables
  implicit none
  ! I/O
  integer                                :: fileunit             !! Unit of the file to read from
  integer                                :: mode                 !! Mode of reading: 1 (native format), 2 (gmsh format)
  ! Local
  character(len=fbem_stdcharlen)         :: section_name         ! Name of the section
  logical                                :: found
  integer                                :: i, j, kn             ! Counters
  character(len=fbem_file_record_length) :: linestr              ! Line
  character(len=fbem_file_record_length) :: tmp_str              ! Temporary string
  integer                                :: nwords
  integer                                :: tmp_n_nodes          ! Temporary variable to read the number of nodes in the file
  integer                                :: tmp_node_eid_max     ! Temporary
  integer                                :: tmp_node_eid_min     ! Temporary
  logical, allocatable                   :: tmp_active_nodes(:)  ! Temporary array
  integer                                :: tmp_id
  integer                                :: nc

  ! Locate the section
  select case (mode)
    case (1)
      section_name='nodes'
      call fbem_search_section(fileunit,section_name,found)
    case (2)
      section_name='Nodes'
      call fbem_search_section_gmsh(fileunit,section_name,found)
  end select
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  if (found) then

    ! Starting  message
    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! ============================================================
    ! CALCULATE THE NUMBER OF NODES THAT BELONGS TO AN ACTIVE PART
    ! ============================================================

    tmp_node_eid_min=1
    tmp_node_eid_max=1
    do i=1,n_elements
      do j=1,element(i)%n_nodes
        if (element(i)%node(j).lt.tmp_node_eid_min) tmp_node_eid_min=element(i)%node(j)
        if (element(i)%node(j).gt.tmp_node_eid_max) tmp_node_eid_max=element(i)%node(j)
      end do
    end do
    allocate (tmp_active_nodes(tmp_node_eid_min:tmp_node_eid_max))
    tmp_active_nodes=.false.
    do i=1,n_elements
      do j=1,element(i)%n_nodes
        tmp_active_nodes(element(i)%node(j))=.true.
      end do
    end do
    n_nodes=0
    do i=tmp_node_eid_min,tmp_node_eid_max
      if (tmp_active_nodes(i)) n_nodes=n_nodes+1
    end do

    ! ==============================
    ! READ THE NODES IN ACTIVE PARTS
    ! ==============================

    ! Allocate and initialize
    if (n_nodes.gt.1) then
      allocate (node(n_nodes))
      do i=1,n_nodes
        allocate(node(i)%x(problem%n),node(i)%x_fn(problem%n),node(i)%n_fn(problem%n),node(i)%t1_fn(problem%n))
        if (problem%n.eq.3) allocate(node(i)%t2_fn(problem%n))
        node(i)%x=0.d0
        node(i)%x_fn=0.d0
        node(i)%n_fn=0.d0
        node(i)%t1_fn=0.d0
        if (problem%n.eq.3) node(i)%t2_fn=0.d0
        node(i)%coupled_node=0
        node(i)%n_dof=0

        ! rigid feature
        node(i)%rigid_link=0
        node(i)%master=0
        node(i)%n_slaves=0

        node(i)%export=.true.
      end do
    else
      call fbem_error_message(error_unit,0,section_name,0,'the number of nodes must be >1')
    end if

    ! Read the number of nodes in the file
    read(fileunit,*) tmp_n_nodes
    if (tmp_n_nodes.lt.n_nodes) then
      call fbem_error_message(error_unit,0,'[nodes]',0,'the number of nodes must be greater')
    end if
    ! Node external identifier and coordinates
    kn=0
    do i=1,tmp_n_nodes
      read(fileunit,'(a)') linestr
      nwords=fbem_count_words(linestr)
      if (nwords.gt.4) then
        call fbem_error_message(error_unit,0,'[nodes]',i,'wrong number of arguments in this line')
      end if
      nwords=min0(nwords-1,problem%n)+1
      read(linestr,*) tmp_id
      if ((tmp_id.ge.tmp_node_eid_min).and.(tmp_id.le.tmp_node_eid_max)) then
        if (tmp_active_nodes(tmp_id)) then
          kn=kn+1
          node(kn)%id=tmp_id
          do j=2,nwords
            tmp_str=fbem_extract_word(linestr,j)
            call fbem_trimall(tmp_str)
            nc=scan(tmp_str,'eEdDqQ')
            if (nc.eq.0) then
              tmp_str(len_trim(tmp_str)+1:len_trim(tmp_str)+2)='D0'
            else
              nc=scan(tmp_str,'eE')
              if (nc.ne.0) tmp_str(nc:nc)='D'
            end if
            read(tmp_str,*) node(kn)%x(j-1)
          end do
          ! Check id
          if (node(kn)%id.le.0) then
            call fbem_error_message(error_unit,0,'node',node(kn)%id,'identifiers must be greater than 0')
          end if
        end if
      end if
    end do

    ! Ending message
    if (verbose_level.ge.2) write(output_unit,'(1x,a)') 'done.'

    ! =================================
    ! CHECK AND BUILD NODES IDENTIFIERS
    ! =================================

    node_eid_min=node(1)%id
    node_eid_max=node(1)%id
    do i=2,n_nodes
      if (node(i)%id.lt.node_eid_min) node_eid_min=node(i)%id
      if (node(i)%id.gt.node_eid_max) node_eid_max=node(i)%id
    end do
    allocate (node_iid(node_eid_min:node_eid_max))
    node_iid=0
    do i=1,n_nodes
      if (node_iid(node(i)%id).ne.0) then
        call fbem_error_message(error_unit,0,'node',node(i)%id,'is repeated')
      else
        node_iid(node(i)%id)=i
      end if
    end do

    ! ====================================================
    ! CHANGE EID TO IID IN THE ELEMENT->NODES CONNECTIVITY
    ! ====================================================

    do i=1,n_elements
      do j=1,element(i)%n_nodes
        if ((element(i)%node(j).ge.node_eid_min).and.(element(i)%node(j).le.node_eid_max)) then
          if (node_iid(element(i)%node(j)).eq.0) then
            call fbem_warning_message(error_unit,0,'node',element(i)%node(j),'does not exist')
            call fbem_error_message(error_unit,0,'element',element(i)%id,'contain the previous node')
          else
            element(i)%node(j)=node_iid(element(i)%node(j))
          end if
        else
          call fbem_warning_message(error_unit,0,'node',element(i)%node(j),'does not exist')
          call fbem_error_message(error_unit,0,'element',element(i)%id,'contain the previous node')
        end if
      end do
    end do

    ! Ending message
    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  else

    call fbem_error_message(error_unit,0,trim(section_name),0,'this section is required')

  end if

end subroutine read_nodes
