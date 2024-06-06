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

subroutine read_commands(fileunit)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_data_structures
  use fbem_mesh_module

  ! Problem variables module
  use problem_variables

  implicit none
  ! I/O
  integer                                :: fileunit
  ! Local
  character(len=fbem_stdcharlen)         :: section_name
  logical                                :: found
  integer                                :: ios           ! Error flag
  integer                                :: nl
  character(len=fbem_file_record_length) :: line, word
  integer                                :: n_words
  integer                                :: ndim, np, kp
  integer, allocatable                   :: sp(:)
  logical                                :: exists

  section_name='commands'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! Read each line
    nl=0
    do while (.true.)
      ! Common processing of all lines
      read(fileunit,'(a)',iostat=ios) line
      if (is_iostat_end(ios)) then
        exit
      end if
      exists=.false.
      nl=nl+1
      call fbem_trim2b(line)
      if (line(1:1).eq.'[') exit
      n_words=fbem_count_words(line)
      if (n_words.eq.0) cycle
      ! Read the command keyword
      word=fbem_extract_word(line,1)
      !
      ! stop
      !
      if (trim(word).eq.'stop') then
        call fbem_timestamp_w_message(output_unit,2,'stop command in section [commands]')
        stop
      end if
      !
      ! exit
      !
      if (trim(word).eq.'exit') exit
      !
      ! convert_mesh_file_format <n> <Input file> <Input format> <Output file> <Output format>
      !
      if (trim(word).eq.'convert_mesh_file_format') then
        exists=.true.
        if (n_words.ne.6) then
          call fbem_error_message(error_unit,0,trim(section_name),nl,'invalid number of arguments in this line.')
        end if
        word=trim(fbem_extract_word(line,2))
        read(word,*) ndim
        if (.not.((ndim.eq.2).or.(ndim.eq.3))) then
          call fbem_error_message(error_unit,0,trim(section_name),nl,'invalid ndim in this line.')
        end if
        call fbem_convert_mesh_file_format(ndim,trim(fbem_extract_word(line,3)),trim(fbem_extract_word(line,4)),&
                                                trim(fbem_extract_word(line,5)),trim(fbem_extract_word(line,6)))
      end if
      !
      ! transform_mesh_parts_to_linear <n> <Input file> <Input format> <Output file> <Output format> <np: number of selected parts> <list selected parts>
      !
      if (trim(word).eq.'transform_mesh_parts_to_linear') then
        exists=.true.
        if (n_words.lt.8) then
          call fbem_error_message(error_unit,0,trim(section_name),nl,'invalid number of arguments in this line.')
        end if
        word=trim(fbem_extract_word(line,2))
        read(word,*) ndim
        if (.not.((ndim.eq.2).or.(ndim.eq.3))) then
          call fbem_error_message(error_unit,0,trim(section_name),nl,'invalid ndim in this line.')
        end if
        word=trim(fbem_extract_word(line,7))
        read(word,*) np
        if (np.le.0) then
          call fbem_error_message(error_unit,0,trim(section_name),nl,'invalid np in this line.')
        end if
        allocate(sp(np))
        do kp=1,np
          word=trim(fbem_extract_word(line,7+kp))
          read(word,*) sp(kp)
        end do
        call fbem_transform_mesh_parts_to_linear(ndim,1.d-6,trim(fbem_extract_word(line,3)),trim(fbem_extract_word(line,4)),&
                                                            trim(fbem_extract_word(line,5)),trim(fbem_extract_word(line,6)),np,sp)
        deallocate(sp)
      end if

      !
      ! Check if the command is not recognized
      !
      if (.not.exists) then
        call fbem_error_message(error_unit,0,trim(word),0,'unrecognized command in [commands].')
      end if
    end do

    ! Ending message
    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  end if

end subroutine read_commands
