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
!! <b> Subroutine that reads the defined groups from a file. </b>

subroutine read_groups(fileunit)

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
  integer                        :: fileunit          ! Unit of the file to read from
  ! Local
  character(len=fbem_stdcharlen) :: section_name      ! Name of the section
  logical                        :: found             ! Logical variable for sections and keywords
  character(len=fbem_stdcharlen) :: tmp_type          ! Temporary variable to read the type of the group
  character(len=fbem_stdcharlen) :: tmp_built         ! Temporary variable to read the built mode of the group
  character(len=fbem_stdcharlen) :: tmp_top
  integer                        :: tmp_top_tag
  integer                        :: tmp_part
  integer                        :: i, j, k, l, m, kj, sn ! Counters
  logical                        :: tmp_check_built
  logical, allocatable           :: tmp_check_part(:)
  real(kind=real64)              :: tmp_radius, tmp_center(problem%n)
  real(kind=real64)              :: xlimits(2,3)
  integer                        :: tmp_n_nodes
  integer, allocatable           :: tmp_node(:)
  logical                        :: tmp_exists, tmp_complete

  section_name='groups'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! Read the number of groups
    read(fileunit,*) n_groups

    ! Allocate data structures
    if (n_groups.gt.0) then
      ! Allocate the vector of groups data structure
      allocate (group(n_groups))
    else
      call fbem_error_message(error_unit,0,section_name,0,'the number of groups must be >1')
    end if

    ! Read each group line
    do i=1,n_groups
      ! Read group id, type of group and the way it is built
      read(fileunit,*) group(i)%id, tmp_type, tmp_built
      ! Check id
      if (group(i)%id.le.0) then
        call fbem_error_message(error_unit,0,'group',group(i)%id,'identifiers must be greater than 0')
      end if
      ! Back to the line
      backspace(fileunit)

      ! ==============
      ! GROUP OF NODES
      ! ==============

      if (trim(tmp_type).eq.'nodes') then
        group(i)%type=fbem_group_type_nodes
        tmp_check_built=.false.
        !
        ! ALL
        !
        if (trim(tmp_built).eq.'all') then
          tmp_check_built=.true.
          group(i)%n_objects=0
          ! Just to read something and move to the next line
          read(fileunit,*) group(i)%id, tmp_type, tmp_built
          do j=1,n_nodes
            if (node(j)%n_parts.gt.0) then
              group(i)%n_objects=group(i)%n_objects+1
            end if
          end do
          allocate (group(i)%object(group(i)%n_objects))
          k=0
          do j=1,n_nodes
            if (node(j)%n_parts.gt.0) then
              k=k+1
              group(i)%object(k)=j
            end if
          end do
        end if
        !
        ! LIST: group of nodes defined by giving a list of nodes
        !
        ! Syntax:
        ! <id_group> nodes list <N> <node_1> <node_2> ... <node_N>
        !
        if (trim(tmp_built).eq.'list') then
          tmp_check_built=.true.
          ! Read and allocate the list
          read(fileunit,*) group(i)%id, tmp_type, tmp_built, group(i)%n_objects
          backspace(fileunit)
          allocate (group(i)%object(group(i)%n_objects))
          read(fileunit,*) group(i)%id, tmp_type, tmp_built, group(i)%n_objects, (group(i)%object(j),j=1,group(i)%n_objects)
          !
          ! CHECK IF THE REFERENCES EXIST
          !
          ! Transform the list from eid to iid, and check if they exist.
          do j=1,group(i)%n_objects
            if ((group(i)%object(j).ge.node_eid_min).and.(group(i)%object(j).le.node_eid_max)) then
              if (node_iid(group(i)%object(j)).eq.0) then
                call fbem_error_message(error_unit,0,'group',group(i)%id,'a node of this group does not exist')
              else
                group(i)%object(j)=node_iid(group(i)%object(j))
              end if
            else
              call fbem_error_message(error_unit,0,'group',group(i)%id,'a node of this group does not exist')
            end if
          end do
          !
          ! CHECK THAT ALL NODES BELONG TO THE SAME PARTS AS THE FIRST NODE
          !
          allocate (tmp_check_part(node(group(i)%object(1))%n_parts))
          do j=2,group(i)%n_objects
            if (node(group(i)%object(j))%n_parts.eq.node(group(i)%object(1))%n_parts) then
              tmp_check_part=.false.
              do k=1,node(group(i)%object(1))%n_parts
                do kj=1,node(group(i)%object(1))%n_parts
                  if (node(group(i)%object(j))%part(kj).eq.node(group(i)%object(1))%part(k)) tmp_check_part(k)=.true.
                end do
              end do
              if (all(tmp_check_part).eqv.(.false.)) then
                call fbem_error_message(error_unit,0,'group',group(i)%id,&
                                        'all nodes must have the same connectivity to parts as the first node')
              end if
            else
              call fbem_error_message(error_unit,0,'group',group(i)%id,&
                                      'all nodes must have the same connectivity to parts as the first node')
            end if
          end do
          deallocate (tmp_check_part)
        end if
        !
        ! PART: group of nodes defined by giving a part
        !
        ! Syntax:
        ! <id_group> nodes part <part_id>
        !
        if (trim(tmp_built).eq.'part') then
          tmp_check_built=.true.
          group(i)%n_objects=0
          ! Read
          read(fileunit,*) group(i)%id, tmp_type, tmp_built, tmp_part
          ! Check tmp_part
          if ((tmp_part.ge.part_eid_min).and.(tmp_part.le.part_eid_max)) then
            if (part_iid(tmp_part).eq.0) then
              call fbem_error_message(error_unit,0,'group',group(i)%id,'the part used for this group does not exist')
            else
              tmp_part=part_iid(tmp_part)
            end if
          else
            call fbem_error_message(error_unit,0,'group',group(i)%id,'the part used for this group does not exist')
          end if
          !write(*,*) 'part', part(tmp_part)%id, 'has', part(tmp_part)%n_nodes, 'nodes'
          group(i)%n_objects=part(tmp_part)%n_nodes
          allocate (group(i)%object(group(i)%n_objects))
          do j=1,part(tmp_part)%n_nodes
            group(i)%object(j)=part(tmp_part)%node(j)
            !write(*,*) 'part', part(tmp_part)%id, 'has node', node(part(tmp_part)%node(j))%id
          end do
        end if
        !
        ! BALL: group of nodes defined by node coordinates in an open region (ball shape)
        !
        ! Syntax:
        ! <id_group> 'nodes' 'ball' <'interior' or 'exterior'> <center> <radius> <id part>
        !
        if (trim(tmp_built).eq.'ball') then
          tmp_check_built=.true.
          group(i)%n_objects=0
          ! Read
          read(fileunit,*) group(i)%id, tmp_type, tmp_built, tmp_top, tmp_center, tmp_radius, tmp_part
          ! Check tmp_top
          tmp_top_tag=0
          if (trim(tmp_top).eq.'interior') tmp_top_tag=1
          if (trim(tmp_top).eq.'exterior') tmp_top_tag=2
          if (tmp_top_tag.eq.0) then
            call fbem_error_message(error_unit,0,'group',group(i)%id,'this parameter must be "interior" or "exterior"')
          end if
          ! Check tmp_part
          if ((tmp_part.ge.part_eid_min).and.(tmp_part.le.part_eid_max)) then
            if (part_iid(tmp_part).eq.0) then
              call fbem_error_message(error_unit,0,'group',group(i)%id,'the part used for this group does not exist')
            else
              tmp_part=part_iid(tmp_part)
            end if
          else
            call fbem_error_message(error_unit,0,'group',group(i)%id,'the part used for this group does not exist')
          end if
          ! Collect the nodes
          select case (tmp_top_tag)
            ! INTERIOR BALL
            case (1)
              do j=1,part(tmp_part)%n_nodes
                sn=part(tmp_part)%node(j)
                if (sqrt(dot_product(node(sn)%x-tmp_center,node(sn)%x-tmp_center)).lt.tmp_radius) then
                  group(i)%n_objects=group(i)%n_objects+1
                end if
              end do
              allocate (group(i)%object(group(i)%n_objects))
              k=0
              do j=1,part(tmp_part)%n_nodes
                sn=part(tmp_part)%node(j)
                if (sqrt(dot_product(node(sn)%x-tmp_center,node(sn)%x-tmp_center)).lt.tmp_radius) then
                  k=k+1
                  group(i)%object(k)=part(tmp_part)%node(j)
                end if
              end do
            ! EXTERIOR BALL
            case (2)
              do j=1,part(tmp_part)%n_nodes
                sn=part(tmp_part)%node(j)
                if (sqrt(dot_product(node(sn)%x-tmp_center,node(sn)%x-tmp_center)).gt.tmp_radius) then
                  group(i)%n_objects=group(i)%n_objects+1
                end if
              end do
              allocate (group(i)%object(group(i)%n_objects))
              k=0
              do j=1,part(tmp_part)%n_nodes
                sn=part(tmp_part)%node(j)
                if (sqrt(dot_product(node(sn)%x-tmp_center,node(sn)%x-tmp_center)).gt.tmp_radius) then
                  k=k+1
                  group(i)%object(k)=part(tmp_part)%node(j)
                end if
              end do
          end select
        end if
        !
        ! BOX : group of nodes defined by node coordinates in an open region (box shape)
        !
        ! Syntax:
        ! <id_group> 'nodes' 'box' <'interior' or 'exterior'> <xmin> <xmax> <ymin> <ymax> [<zmin> <zmax>] <id part>
        !
        if (trim(tmp_built).eq.'box') then
          tmp_check_built=.true.
          group(i)%n_objects=0
          ! Read

          ! falta chequear el numero de argumentos, utilizar la funcion de contar y extraer palabras
          read(fileunit,*) group(i)%id, tmp_type, tmp_built, tmp_top, xlimits(:,1:problem%n), tmp_part


          ! Check tmp_top
          tmp_top_tag=0
          if (trim(tmp_top).eq.'interior') tmp_top_tag=1
          if (trim(tmp_top).eq.'exterior') tmp_top_tag=2
          if (tmp_top_tag.eq.0) then
            call fbem_error_message(error_unit,0,'group',group(i)%id,'this parameter must be "interior" or "exterior"')
          end if
          ! Check tmp_part
          if ((tmp_part.ge.part_eid_min).and.(tmp_part.le.part_eid_max)) then
            if (part_iid(tmp_part).eq.0) then
              call fbem_error_message(error_unit,0,'group',group(i)%id,'the part used for this group does not exist')
            else
              tmp_part=part_iid(tmp_part)
            end if
          else
            call fbem_error_message(error_unit,0,'group',group(i)%id,'the part used for this group does not exist')
          end if
          ! Collect the nodes
          select case (tmp_top_tag)
            ! INTERIOR BOX
            case (1)
              select case (problem%n)
                case (2)
                  do j=1,part(tmp_part)%n_nodes
                    sn=part(tmp_part)%node(j)
                    if ((node(sn)%x(1).gt.xlimits(1,1)).and.(node(sn)%x(1).lt.xlimits(2,1)).and.&
                        (node(sn)%x(2).gt.xlimits(1,2)).and.(node(sn)%x(2).lt.xlimits(2,2))) then
                      group(i)%n_objects=group(i)%n_objects+1
                    end if
                  end do
                  allocate (group(i)%object(group(i)%n_objects))
                  k=0
                  do j=1,part(tmp_part)%n_nodes
                    sn=part(tmp_part)%node(j)
                    if ((node(sn)%x(1).gt.xlimits(1,1)).and.(node(sn)%x(1).lt.xlimits(2,1)).and.&
                        (node(sn)%x(2).gt.xlimits(1,2)).and.(node(sn)%x(2).lt.xlimits(2,2))) then
                      k=k+1
                      group(i)%object(k)=part(tmp_part)%node(j)
                    end if
                  end do
                case (3)
                  do j=1,part(tmp_part)%n_nodes
                    sn=part(tmp_part)%node(j)
                    if ((node(sn)%x(1).gt.xlimits(1,1)).and.(node(sn)%x(1).lt.xlimits(2,1)).and.&
                        (node(sn)%x(2).gt.xlimits(1,2)).and.(node(sn)%x(2).lt.xlimits(2,2)).and.&
                        (node(sn)%x(3).gt.xlimits(1,3)).and.(node(sn)%x(3).lt.xlimits(2,3))) then
                      group(i)%n_objects=group(i)%n_objects+1
                    end if
                  end do
                  allocate (group(i)%object(group(i)%n_objects))
                  k=0
                  do j=1,part(tmp_part)%n_nodes
                    sn=part(tmp_part)%node(j)
                    if ((node(sn)%x(1).gt.xlimits(1,1)).and.(node(sn)%x(1).lt.xlimits(2,1)).and.&
                        (node(sn)%x(2).gt.xlimits(1,2)).and.(node(sn)%x(2).lt.xlimits(2,2)).and.&
                        (node(sn)%x(3).gt.xlimits(1,3)).and.(node(sn)%x(3).lt.xlimits(2,3))) then
                      k=k+1
                      group(i)%object(k)=part(tmp_part)%node(j)
                    end if
                  end do
              end select
            ! EXTERIOR BOX
            case (2)
              select case (problem%n)
                case (2)
                  do j=1,part(tmp_part)%n_nodes
                    sn=part(tmp_part)%node(j)
                    if ((node(sn)%x(1).lt.xlimits(1,1)).or.(node(sn)%x(1).gt.xlimits(2,1)).or.&
                        (node(sn)%x(2).lt.xlimits(1,2)).or.(node(sn)%x(2).gt.xlimits(2,2))) then
                      group(i)%n_objects=group(i)%n_objects+1
                    end if
                  end do
                  allocate (group(i)%object(group(i)%n_objects))
                  k=0
                  do j=1,part(tmp_part)%n_nodes
                    sn=part(tmp_part)%node(j)
                    if ((node(sn)%x(1).lt.xlimits(1,1)).or.(node(sn)%x(1).gt.xlimits(2,1)).or.&
                        (node(sn)%x(2).lt.xlimits(1,2)).or.(node(sn)%x(2).gt.xlimits(2,2))) then
                      k=k+1
                      group(i)%object(k)=part(tmp_part)%node(j)
                    end if
                  end do
                case (3)
                  do j=1,part(tmp_part)%n_nodes
                    sn=part(tmp_part)%node(j)
                    if ((node(sn)%x(1).lt.xlimits(1,1)).or.(node(sn)%x(1).gt.xlimits(2,1)).or.&
                        (node(sn)%x(2).lt.xlimits(1,2)).or.(node(sn)%x(2).gt.xlimits(2,2)).or.&
                        (node(sn)%x(3).lt.xlimits(1,3)).or.(node(sn)%x(3).gt.xlimits(2,3))) then
                      group(i)%n_objects=group(i)%n_objects+1
                    end if
                  end do
                  allocate (group(i)%object(group(i)%n_objects))
                  k=0
                  do j=1,part(tmp_part)%n_nodes
                    sn=part(tmp_part)%node(j)
                    if ((node(sn)%x(1).lt.xlimits(1,1)).or.(node(sn)%x(1).gt.xlimits(2,1)).or.&
                        (node(sn)%x(2).lt.xlimits(1,2)).or.(node(sn)%x(2).gt.xlimits(2,2)).or.&
                        (node(sn)%x(3).lt.xlimits(1,3)).or.(node(sn)%x(3).gt.xlimits(2,3))) then
                      k=k+1
                      group(i)%object(k)=part(tmp_part)%node(j)
                    end if
                  end do
              end select
          end select
        end if
        !
        ! CASE DEFAULT
        !
        if (tmp_check_built.eqv.(.false.)) call fbem_error_message(error_unit,0,'section_name',group(i)%id,'unknown built mode')
      end if

      ! =================
      ! GROUP OF ELEMENTS
      ! =================

      if (trim(tmp_type).eq.'elements') then
        group(i)%type=fbem_group_type_elements
        tmp_check_built=.false.
        !
        ! LIST
        !
        if (trim(tmp_built).eq.'list') then
          tmp_check_built=.true.
          ! Read and allocate the list
          read(fileunit,*) group(i)%id, tmp_type, tmp_built, group(i)%n_objects
          backspace(fileunit)
          allocate (group(i)%object(group(i)%n_objects))
          read(fileunit,*) group(i)%id, tmp_type, tmp_built, group(i)%n_objects, (group(i)%object(j),j=1,group(i)%n_objects)
          !
          ! Check if the references exist
          !
          ! Transform the list from eid to iid, and check if they exist.
          do j=1,group(i)%n_objects
            if ((group(i)%object(j).ge.element_eid_min).and.(group(i)%object(j).le.element_eid_max)) then
              if (element_iid(group(i)%object(j)).eq.0) then
                call fbem_error_message(error_unit,0,'group',group(i)%id,'an element of this group does not exist')
              else
                group(i)%object(j)=element_iid(group(i)%object(j))
              end if
            else
              call fbem_error_message(error_unit,0,'group',group(i)%id,'an element of this group does not exist')
            end if
          end do
          !
          ! CHECK THAT ALL ELEMENTS BELONG TO THE SAME PART AS THE FIRST ELEMENT
          !
          do j=2,group(i)%n_objects
            if (element(group(i)%object(j))%part.ne.element(group(i)%object(1))%part) then
              call fbem_error_message(error_unit,0,'group',group(i)%id,&
                                      'all elements must belong to the same part as the first element')
            end if
          end do
        end if
        !
        ! ALL
        !
        if (trim(tmp_built).eq.'all') then
          tmp_check_built=.true.
          group(i)%n_objects=0

          ! Just to read something and move to the next line
          read(fileunit,*) group(i)%id, tmp_type, tmp_built

          do j=1,n_elements
            if (element(j)%part.gt.0) then
              group(i)%n_objects=group(i)%n_objects+1
            end if
          end do
          allocate (group(i)%object(group(i)%n_objects))
          k=0
          do j=1,n_elements
            if (element(j)%part.gt.0) then
              k=k+1
              group(i)%object(k)=j
            end if
          end do

        end if
        !
        ! CASE DEFAULT
        !
        if (tmp_check_built.eqv.(.false.)) call fbem_error_message(error_unit,0,'section_name',group(i)%id,'unknown built mode')
      end if

      ! =================
      ! GROUP OF SUBEDGES
      ! =================

      if (trim(tmp_type).eq.'subedges') then
        group(i)%type=fbem_group_type_subedges
        tmp_check_built=.false.
        !
        ! BALL (open region)
        !
        ! Syntax:
        ! <id_group> 'subedges' 'ball' <'interior' or 'exterior'> <center> <radius> <part where the subedges are selected>
        !
        if (trim(tmp_built).eq.'ball') then
          tmp_check_built=.true.
          ! Read
          read(fileunit,*) group(i)%id, tmp_type, tmp_built, tmp_top, tmp_center, tmp_radius, tmp_part
          ! Check tmp_top
          tmp_top_tag=0
          if (trim(tmp_top).eq.'interior') tmp_top_tag=1
          if (trim(tmp_top).eq.'exterior') tmp_top_tag=2
          if (tmp_top_tag.eq.0) then
            call fbem_error_message(error_unit,0,'group',group(i)%id,'this parameter must be "interior" or "exterior"')
          end if
          ! Check tmp_part
          if ((tmp_part.ge.part_eid_min).and.(tmp_part.le.part_eid_max)) then
            if (part_iid(tmp_part).eq.0) then
              call fbem_error_message(error_unit,0,'group',group(i)%id,'the part used for this group does not exist')
            else
              tmp_part=part_iid(tmp_part)
            end if
          else
            call fbem_error_message(error_unit,0,'group',group(i)%id,'the part used for this group does not exist')
          end if
          ! Collect the nodes
          select case (tmp_top_tag)
            ! INTERIOR BALL
            case (1)
              tmp_n_nodes=0
              do j=1,part(tmp_part)%n_nodes
                sn=part(tmp_part)%node(j)
                if (sqrt(dot_product(node(sn)%x-tmp_center,node(sn)%x-tmp_center)).lt.tmp_radius) then
                  tmp_n_nodes=tmp_n_nodes+1
                end if
              end do
              allocate (tmp_node(tmp_n_nodes))
              k=0
              do j=1,part(tmp_part)%n_nodes
                sn=part(tmp_part)%node(j)
                if (sqrt(dot_product(node(sn)%x-tmp_center,node(sn)%x-tmp_center)).lt.tmp_radius) then
                  k=k+1
                  tmp_node(k)=part(tmp_part)%node(j)
                end if
              end do
            ! EXTERIOR BALL
            case (2)
              tmp_n_nodes=0
              do j=1,part(tmp_part)%n_nodes
                sn=part(tmp_part)%node(j)
                if (sqrt(dot_product(node(sn)%x-tmp_center,node(sn)%x-tmp_center)).gt.tmp_radius) then
                  tmp_n_nodes=tmp_n_nodes+1
                end if
              end do
              allocate (tmp_node(tmp_n_nodes))
              k=0
              do j=1,part(tmp_part)%n_nodes
                sn=part(tmp_part)%node(j)
                if (sqrt(dot_product(node(sn)%x-tmp_center,node(sn)%x-tmp_center)).gt.tmp_radius) then
                  k=k+1
                  tmp_node(k)=part(tmp_part)%node(j)
                end if
              end do
          end select
        end if
        !
        ! BOX (open region)
        !
        ! Syntax:
        ! <id_group> 'subedges' 'box' <'interior' or 'exterior'> <xmin> <xmax> <ymin> <ymax> [<zmin> <zmax>]
        !
        if (trim(tmp_built).eq.'box') then
          tmp_check_built=.true.
          ! Read
          read(fileunit,*) group(i)%id, tmp_type, tmp_built, tmp_top, xlimits(:,1:problem%n), tmp_part
          ! Check tmp_top
          tmp_top_tag=0
          if (trim(tmp_top).eq.'interior') tmp_top_tag=1
          if (trim(tmp_top).eq.'exterior') tmp_top_tag=2
          if (tmp_top_tag.eq.0) then
            call fbem_error_message(error_unit,0,'group',group(i)%id,'this parameter must be "interior" or "exterior"')
          end if
          ! Check tmp_part
          if ((tmp_part.ge.part_eid_min).and.(tmp_part.le.part_eid_max)) then
            if (part_iid(tmp_part).eq.0) then
              call fbem_error_message(error_unit,0,'group',group(i)%id,'the part used for this group does not exist')
            else
              tmp_part=part_iid(tmp_part)
            end if
          else
            call fbem_error_message(error_unit,0,'group',group(i)%id,'the part used for this group does not exist')
          end if
          ! Collect the nodes
          select case (tmp_top_tag)
            ! INTERIOR BOX
            case (1)
              select case (problem%n)
                case (2)
                  tmp_n_nodes=0
                  do j=1,part(tmp_part)%n_nodes
                    sn=part(tmp_part)%node(j)
                    if ((node(sn)%x(1).gt.xlimits(1,1)).and.(node(sn)%x(1).lt.xlimits(2,1)).and.&
                        (node(sn)%x(2).gt.xlimits(1,2)).and.(node(sn)%x(2).lt.xlimits(2,2))) then
                      tmp_n_nodes=tmp_n_nodes+1
                    end if
                  end do
                  allocate (tmp_node(tmp_n_nodes))
                  k=0
                  do j=1,part(tmp_part)%n_nodes
                    sn=part(tmp_part)%node(j)
                    if ((node(sn)%x(1).gt.xlimits(1,1)).and.(node(sn)%x(1).lt.xlimits(2,1)).and.&
                        (node(sn)%x(2).gt.xlimits(1,2)).and.(node(sn)%x(2).lt.xlimits(2,2))) then
                      k=k+1
                      tmp_node(k)=part(tmp_part)%node(j)
                    end if
                  end do
                case (3)
                  tmp_n_nodes=0
                  do j=1,part(tmp_part)%n_nodes
                    sn=part(tmp_part)%node(j)
                    if ((node(sn)%x(1).gt.xlimits(1,1)).and.(node(sn)%x(1).lt.xlimits(2,1)).and.&
                        (node(sn)%x(2).gt.xlimits(1,2)).and.(node(sn)%x(2).lt.xlimits(2,2)).and.&
                        (node(sn)%x(3).gt.xlimits(1,3)).and.(node(sn)%x(3).lt.xlimits(2,3))) then
                      tmp_n_nodes=tmp_n_nodes+1
                    end if
                  end do
                  allocate (tmp_node(tmp_n_nodes))
                  k=0
                  do j=1,part(tmp_part)%n_nodes
                    sn=part(tmp_part)%node(j)
                    if ((node(sn)%x(1).gt.xlimits(1,1)).and.(node(sn)%x(1).lt.xlimits(2,1)).and.&
                        (node(sn)%x(2).gt.xlimits(1,2)).and.(node(sn)%x(2).lt.xlimits(2,2)).and.&
                        (node(sn)%x(3).gt.xlimits(1,3)).and.(node(sn)%x(3).lt.xlimits(2,3))) then
                      k=k+1
                      tmp_node(k)=part(tmp_part)%node(j)
                    end if
                  end do
              end select
            ! EXTERIOR BOX
            case (2)
              select case (problem%n)
                case (2)
                  tmp_n_nodes=0
                  do j=1,part(tmp_part)%n_nodes
                    sn=part(tmp_part)%node(j)
                    if ((node(sn)%x(1).lt.xlimits(1,1)).or.(node(sn)%x(1).gt.xlimits(2,1)).or.&
                        (node(sn)%x(2).lt.xlimits(1,2)).or.(node(sn)%x(2).gt.xlimits(2,2))) then
                      tmp_n_nodes=tmp_n_nodes+1
                    end if
                  end do
                  allocate (tmp_node(tmp_n_nodes))
                  k=0
                  do j=1,part(tmp_part)%n_nodes
                    sn=part(tmp_part)%node(j)
                    if ((node(sn)%x(1).lt.xlimits(1,1)).or.(node(sn)%x(1).gt.xlimits(2,1)).or.&
                        (node(sn)%x(2).lt.xlimits(1,2)).or.(node(sn)%x(2).gt.xlimits(2,2))) then
                      k=k+1
                      tmp_node(k)=part(tmp_part)%node(j)
                    end if
                  end do
                case (3)
                  tmp_n_nodes=0
                  do j=1,part(tmp_part)%n_nodes
                    sn=part(tmp_part)%node(j)
                    if ((node(sn)%x(1).lt.xlimits(1,1)).or.(node(sn)%x(1).gt.xlimits(2,1)).or.&
                        (node(sn)%x(2).lt.xlimits(1,2)).or.(node(sn)%x(2).gt.xlimits(2,2)).or.&
                        (node(sn)%x(3).lt.xlimits(1,3)).or.(node(sn)%x(3).gt.xlimits(2,3))) then
                      tmp_n_nodes=tmp_n_nodes+1
                    end if
                  end do
                  allocate (tmp_node(tmp_n_nodes))
                  k=0
                  do j=1,part(tmp_part)%n_nodes
                    sn=part(tmp_part)%node(j)
                    if ((node(sn)%x(1).lt.xlimits(1,1)).or.(node(sn)%x(1).gt.xlimits(2,1)).or.&
                        (node(sn)%x(2).lt.xlimits(1,2)).or.(node(sn)%x(2).gt.xlimits(2,2)).or.&
                        (node(sn)%x(3).lt.xlimits(1,3)).or.(node(sn)%x(3).gt.xlimits(2,3))) then
                      k=k+1
                      tmp_node(k)=part(tmp_part)%node(j)
                    end if
                  end do
              end select
          end select
        end if
        ! Collect those subedges belonging to the indicated part where all the nodes are in the previously built list
        group(i)%n_objects=0
        do j=1,n_subedges
          do kj=1,subedge(j)%n_supelements
            if (element(subedge(j)%supelement(kj))%part.eq.tmp_part) then
              tmp_complete=.true.
              do l=1,subedge(j)%n_nodes
                sn=subedge(j)%node(l)
                tmp_exists=.false.
                do m=1,tmp_n_nodes
                  if (tmp_node(m).eq.sn) then
                    tmp_exists=.true.
                  end if
                end do
                if (tmp_exists.eqv.(.false.)) then
                  tmp_complete=.false.
                  exit
                end if
              end do
              if (tmp_complete) group(i)%n_objects=group(i)%n_objects+1
              exit
            end if
          end do
        end do
        allocate (group(i)%object(group(i)%n_objects))
        k=0
        do j=1,n_subedges
          do kj=1,subedge(j)%n_supelements
            if (element(subedge(j)%supelement(kj))%part.eq.tmp_part) then
              tmp_complete=.true.
              do l=1,subedge(j)%n_nodes
                sn=subedge(j)%node(l)
                tmp_exists=.false.
                do m=1,tmp_n_nodes
                  if (tmp_node(m).eq.sn) then
                    tmp_exists=.true.
                  end if
                end do
                if (tmp_exists.eqv.(.false.)) then
                  tmp_complete=.false.
                  exit
                end if
              end do
              if (tmp_complete) then
                k=k+1
                group(i)%object(k)=j
              end if
              exit
            end if
          end do
        end do
        deallocate (tmp_node)
        !
        ! CASE DEFAULT
        !
        if (tmp_check_built.eqv.(.false.)) call fbem_error_message(error_unit,0,'section_name',group(i)%id,'unknown built mode')
      end if

      if (.not.tmp_check_built) then
        call fbem_error_message(error_unit,0,'group',group(i)%id,'unknown type of group')
      end if

    end do

    ! =================================
    ! CHECK AND BUILD GROUP IDENTIFIERS
    ! =================================

    group_eid_min=group(1)%id
    group_eid_max=group(1)%id
    do i=2,n_groups
      if (group(i)%id.lt.group_eid_min) group_eid_min=group(i)%id
      if (group(i)%id.gt.group_eid_max) group_eid_max=group(i)%id
    end do
    allocate (group_iid(group_eid_min:group_eid_max))
    group_iid=0
    do i=1,n_groups
      if (group_iid(group(i)%id).ne.0) then
        call fbem_error_message(error_unit,0,'group',group(i)%id,'is repeated')
      else
        group_iid(group(i)%id)=i
      end if
    end do

    ! ====================================
    ! CHECK THAT EACH GROUP HAS >0 OBJECTS
    ! ====================================

    do i=1,n_groups
      if (group(i)%n_objects.eq.0) then
        call fbem_error_message(error_unit,0,'group',group(i)%id,'has 0 objects')
      end if
    end do

    ! ========
    ! LOG FILE
    ! ========

    if (verbose_level.ge.4) then
      do i=1,n_groups
        select case (group(i)%type)
          case (fbem_group_type_nodes)
            write(*,*) 'Group of nodes', group(i)%id, '. N: ', group(i)%n_objects
            do j=1,group(i)%n_objects
              write(*,*) node(group(i)%object(j))%id
            end do
          case (fbem_group_type_elements)
            write(*,*) 'Group of elements', group(i)%id, '. N: ', group(i)%n_objects
            do j=1,group(i)%n_objects
              write(*,*) element(group(i)%object(j))%id
            end do
          case (fbem_group_type_subedges)
            write(*,*) 'Group of subedges', group(i)%id, '. N: ', group(i)%n_objects
            do j=1,group(i)%n_objects
              do k=1,subedge(group(i)%object(j))%n_nodes
                write(*,*) node(subedge(group(i)%object(j))%node(k))%x
              end do
              write(*,*)
            end do
        end select
      end do
    end if

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  else

    n_groups=0

  end if

end subroutine read_groups
