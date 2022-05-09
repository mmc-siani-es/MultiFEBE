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

subroutine read_element_options(input_fileunit)

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
  integer                                :: kn, sn
  integer                                :: nl, ios, nc, nw, eid, iid, tmp_int
  real(kind=real64)                      :: x1ref(3)
  character(len=fbem_string_max_length)  :: section_name
  character(len=fbem_file_record_length) :: line, word
  character(len=fbem_file_record_length) :: entity, option

  ! Section name
  section_name='element options'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
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
        if ((trim(entity).eq.'element').or.(trim(entity).eq.'group')) then
          ! The second word must be a valid id
          word=fbem_extract_word(line,2)
          ! Read an integer from word
          read(word,*) eid

          ! Check if a valid element
          if (trim(entity).eq.'element') then
            if ((eid.ge.element_eid_min).and.(eid.le.element_eid_max)) then
              if (element_iid(eid).eq.0) then
                call fbem_error_message(error_unit,0,section_name,nl,'the element identifier given in this line does not exist.')
              else
                iid=element_iid(eid)
              end if
            else
              call fbem_error_message(error_unit,0,section_name,nl,'the element identifier given in this line does not exist.')
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
            if (group(iid)%type.ne.fbem_group_type_elements) then
              call fbem_error_message(error_unit,0,section_name,nl,'this group is not a group of elements.')
            end if
          end if

          ! The second word must be a valid option
          option=fbem_extract_word(line,3)

          !
          ! strbeam_line3_midnode_rotation (0: no rotation, 1: rotation)
          !
          if (trim(option).eq.'strbeam_line3_midnode_rotation') then

            word=fbem_extract_word(line,4)
            read(word,*) tmp_int

            if (.not.((tmp_int.eq.0).or.(tmp_int.eq.1))) then
              call fbem_error_message(error_unit,0,section_name,nl,'strbeam_line3_midnode_rotation value must be 0 or 1')
            end if

            if (trim(entity).eq.'element') then
              if (.not.((element(iid)%type.eq.fbem_line3).and.((element(iid)%fe_type.eq.1).or.(element(iid)%fe_type.eq.2)))) then
                call fbem_error_message(error_unit,0,section_name,nl,'this option is not applicable to this element')
              end if
              element(iid)%fe_options(1)=tmp_int
            end if

            if (trim(entity).eq.'group') then
              do kn=1,group(iid)%n_objects
                sn=group(iid)%object(kn)
                if (.not.((element(sn)%type.eq.fbem_line3).and.((element(sn)%fe_type.eq.1).or.(element(sn)%fe_type.eq.2)))) then
                  call fbem_error_message(error_unit,0,section_name,nl,'this option is not applicable to an element of this group')
                end if
                element(sn)%fe_options(1)=tmp_int
              end do
            end if
          !
          ! bar_mass_matrix (0: consistent, 1: lumped)
          !
          else if (trim(option).eq.'bar_mass_matrix') then

            word=fbem_extract_word(line,4)
            read(word,*) tmp_int

            if (.not.((tmp_int.eq.0).or.(tmp_int.eq.1))) then
              call fbem_error_message(error_unit,0,section_name,nl,'bar_mass_matrix value must be 0 (consistent) or 1 (lumped)')
            end if

            if (trim(entity).eq.'element') then
              if (.not.((element(iid)%type.eq.fbem_line2).and.(element(iid)%fe_type.eq.3))) then
                call fbem_error_message(error_unit,0,section_name,nl,'this option is not applicable to this element')
              end if
              element(iid)%fe_options(1)=tmp_int
            end if

            if (trim(entity).eq.'group') then
              do kn=1,group(iid)%n_objects
                sn=group(iid)%object(kn)
                if (.not.((element(sn)%type.eq.fbem_line2).and.(element(sn)%fe_type.eq.3))) then
                  call fbem_error_message(error_unit,0,section_name,nl,'this option is not applicable to an element of this group')
                end if
                element(sn)%fe_options(1)=tmp_int
              end do
            end if
          !
          ! degbeam_K_integration (full, reduced, selective, user-defined)
          !
          else if (trim(option).eq.'degbeam_K_integration') then

            word=fbem_extract_word(line,4)

            if (trim(entity).eq.'element') then

              if (.not.((element(iid)%n_dimension.eq.1).and.(element(iid)%fe_type.eq.0))) then
                call fbem_error_message(error_unit,0,section_name,nl,'this option is not applicable to this element')
              end if

              if (trim(word).eq.'full') then
                element(iid)%K_intmode=0
              else if (trim(word).eq.'reduced') then
                element(iid)%K_intmode=1
              else if (trim(word).eq.'selective') then
                element(iid)%K_intmode=2
              else if (trim(word).eq.'user_defined') then
                element(iid)%K_intmode=3
                ! Longitudinal
                word=fbem_extract_word(line,5)
                read(word,*) tmp_int
                if (tmp_int.le.0) then
                  call fbem_error_message(error_unit,0,section_name,nl,'the quadrature rule for degbeam_K_integration user-defined mode must be >0')
                end if
                element(iid)%K_intngp(1)=tmp_int
                ! Shear
                word=fbem_extract_word(line,6)
                read(word,*) tmp_int
                if (tmp_int.le.0) then
                  call fbem_error_message(error_unit,0,section_name,nl,'the quadrature rule for degbeam_K_integration user-defined mode must be >0')
                end if
                element(iid)%K_intngp(2)=tmp_int
                ! Thickness direction
                word=fbem_extract_word(line,7)
                read(word,*) tmp_int
                if (tmp_int.le.0) then
                  call fbem_error_message(error_unit,0,section_name,nl,'the quadrature rule for degbeam_K_integration user-defined mode must be >0')
                end if
                element(iid)%K_intngp(3)=tmp_int
              else
                call fbem_error_message(error_unit,0,section_name,nl,'degbeam_K_integration must be full, reduced, selective or user-defined')
              end if

            end if

            if (trim(entity).eq.'group') then
              do kn=1,group(iid)%n_objects
                sn=group(iid)%object(kn)

                if (.not.((element(sn)%n_dimension.eq.1).and.(element(sn)%fe_type.eq.0))) then
                  call fbem_error_message(error_unit,0,section_name,nl,'this option is not applicable to an element of this group')
                end if

                if (trim(word).eq.'full') then
                  element(sn)%K_intmode=0
                else if (trim(word).eq.'reduced') then
                  element(sn)%K_intmode=1
                else if (trim(word).eq.'selective') then
                  element(sn)%K_intmode=2
                else if (trim(word).eq.'user_defined') then
                  element(sn)%K_intmode=3
                  ! Longitudinal
                  word=fbem_extract_word(line,5)
                  read(word,*) tmp_int
                  if (tmp_int.le.0) then
                    call fbem_error_message(error_unit,0,section_name,nl,'the quadrature rule for degbeam_K_integration user-defined mode must be >0')
                  end if
                  element(sn)%K_intngp(1)=tmp_int
                  ! Shear
                  word=fbem_extract_word(line,6)
                  read(word,*) tmp_int
                  if (tmp_int.le.0) then
                    call fbem_error_message(error_unit,0,section_name,nl,'the quadrature rule for degbeam_K_integration user-defined mode must be >0')
                  end if
                  element(sn)%K_intngp(2)=tmp_int
                  ! Thickness direction
                  word=fbem_extract_word(line,7)
                  read(word,*) tmp_int
                  if (tmp_int.le.0) then
                    call fbem_error_message(error_unit,0,section_name,nl,'the quadrature rule for degbeam_K_integration user-defined mode must be >0')
                  end if
                  element(sn)%K_intngp(3)=tmp_int
                else
                  call fbem_error_message(error_unit,0,section_name,nl,'degbeam_K_integration must be full, reduced, selective or user-defined')
                end if

              end do

            end if
          !
          ! shell_orientation
          !
          else if (trim(option).eq.'shell_orientation') then
            word=fbem_extract_word(line,4)
            read(word,*) x1ref(1)
            word=fbem_extract_word(line,5)
            read(word,*) x1ref(2)
            word=fbem_extract_word(line,6)
            read(word,*) x1ref(3)
            if (trim(entity).eq.'element') then
              if (.not.allocated(element(iid)%ep)) then
                allocate(element(iid)%ep(3,3))
                element(iid)%ep=0
              end if
              element(iid)%ep(:,1)=x1ref
            end if
            if (trim(entity).eq.'group') then
              do kn=1,group(iid)%n_objects
                sn=group(iid)%object(kn)
                if (.not.allocated(element(sn)%ep)) then
                  allocate(element(sn)%ep(3,3))
                end if
                element(sn)%ep(:,1)=x1ref
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



end subroutine read_element_options
