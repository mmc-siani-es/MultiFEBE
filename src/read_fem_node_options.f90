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

subroutine read_fem_node_options(input_fileunit)

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
  integer                                :: input_fileunit    !! Input file unit
  ! Local variables
  logical                                :: found
  integer                                :: kn, sn, kc
  integer                                :: nl, ios, nc, nw, eid, iid, tmp_int
  real(kind=real64)                      :: tmp_real, I2(3)
  character(len=fbem_string_max_length)  :: section_name
  character(len=fbem_file_record_length) :: line, word
  character(len=fbem_file_record_length) :: entity, option

  section_name='fem node options'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')
  call fbem_search_section(input_fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! Initialize line counter
    nl=1
    do
      ! Read the line and exit if the end of the file
      read(input_fileunit,'(a)',iostat=ios) line
      if (is_iostat_end(ios)) exit
      ! Remove beginning, ending and duplicated blanks of a string
      call fbem_trim2b(line)
      ! Count the number of characters
      nc=len_trim(line)
      ! If the line is not empty, continue processing the line
      if (nc.gt.0) then
        ! If the line starts with "[", then it is the end of the section
        if (line(1:1).eq.'[') exit
        ! If the line starts with "#", then it is a comment
        if (line(1:1).eq.'#') cycle
        ! Count the number of words
        nw=fbem_count_words(line)
        !
        ! First word is the type of entity
        !
        entity=fbem_extract_word(line,1)
        if ((trim(entity).eq.'node').or.(trim(entity).eq.'group')) then
          ! The second word must be a valid id
          word=fbem_extract_word(line,2)
          ! Read an integer from word
          read(word,*) eid
          ! Check if a valid node
          if (trim(entity).eq.'node') then
            if ((eid.ge.node_eid_min).and.(eid.le.node_eid_max)) then
              if (node_iid(eid).eq.0) then
                call fbem_error_message(error_unit,0,section_name,nl,'the node identifier given in this line does not exist.')
              else
                iid=node_iid(eid)
              end if
            else
              call fbem_error_message(error_unit,0,section_name,nl,'the node identifier given in this line does not exist.')
            end if
            if (part(node(iid)%part(1))%type.ne.fbem_part_fe_subregion) then
              call fbem_error_message(error_unit,0,section_name,nl,'this node is not a FEM node.')
            end if
          end if
          ! Check if valid group
          if (trim(entity).eq.'group') then
            if ((eid.ge.group_eid_min).and.(eid.le.group_eid_max)) then
              if (group_iid(eid).eq.0) then
                call fbem_error_message(error_unit,0,section_name,nl,'the group identifier given in this line does not exist.')
              else
                iid=group_iid(eid)
              end if
            else
              call fbem_error_message(error_unit,0,section_name,nl,'the group identifier given in this line does not exist.')
            end if
            if (group(iid)%type.ne.1) then
              call fbem_error_message(error_unit,0,section_name,nl,'this group is not a group of nodes.')
            end if
            if (part(node(group(iid)%object(1))%part(1))%type.ne.fbem_part_fe_subregion) then
              call fbem_error_message(error_unit,0,section_name,nl,'this group is not a FE subregion group of nodes.')
            end if
          end if
          ! The second word must be a valid option
          option=fbem_extract_word(line,3)
          ! Modify ndof for nodes or group of nodes
          !
          ! N_DOF
          !
          if (trim(option).eq.'n_dof') then
            word=fbem_extract_word(line,4)
            read(word,*) tmp_int
            if (trim(entity).eq.'node') then
              node(iid)%n_dof=tmp_int
            end if
            if (trim(entity).eq.'group') then
              do kn=1,group(iid)%n_objects
                sn=group(iid)%object(kn)
                node(sn)%n_dof=tmp_int
              end do
            end if
          !
          ! IS_SINGULAR
          !
          else if (trim(option).eq.'is_singular') then
            if (trim(entity).eq.'node') then
              node(iid)%is_singular=.true.
            end if
            if (trim(entity).eq.'group') then
              do kn=1,group(iid)%n_objects
                sn=group(iid)%object(kn)
                node(sn)%is_singular=.true.
              end do
            end if
          else
            call fbem_error_message(error_unit,0,section_name,nl,'invalid option')
          end if
        else
          call fbem_error_message(error_unit,0,section_name,nl,'invalid syntax in this line of the section.')
        end if
      end if
      ! Increment line counter
      nl=nl+1
    end do

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  end if

end subroutine read_fem_node_options
