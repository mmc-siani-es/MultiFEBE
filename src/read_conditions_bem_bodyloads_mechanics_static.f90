! ---------------------------------------------------------------------
! Copyright (C) 2014-2024 Universidad de Las Palmas de Gran Canaria:
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

subroutine read_conditions_bem_bodyloads_mechanics_static(input_fileunit)

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

  ! I/O variables
  integer                                :: input_fileunit    !! Input file unit
  ! Local variables
  character(len=fbem_stdcharlen)         :: section_name      ! Name of the section
  logical                                :: found
  integer                                :: ios           ! Error flag
  integer                                :: nl
  character(len=fbem_file_record_length) :: line, word
  integer                                :: n_words, kword
  integer                                :: entity_eid, i
  integer                                :: info
  logical                                :: exists

  ! Locate the section
  section_name='conditions over be bodyloads'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(input_fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! Read each line
    nl=0
    do while (.true.)
      ! Common processing of all lines
      read(input_fileunit,'(a)',iostat=ios) line
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
      ! be_bodyload
      !
      if (trim(word).eq.'be_bodyload') then
        exists=.true.
        if (n_words.lt.2) then
          call fbem_error_message(error_unit,0,trim(section_name),nl,'invalid number of arguments in this section line.')
        end if
        ! id
        word=trim(fbem_extract_word(line,2))
        read(word,*) entity_eid
        info = check_be_bodyload_eid(entity_eid)
        if (info.ne.0) call fbem_error_message(error_unit,0,trim(section_name),nl,'in this section line, the indicated be_bodyload does not exist.')
        i = be_bodyload_iid(entity_eid)
        kword = 3
        select case (be_bodyload(i)%coupling)
          !
          ! UNCOUPLED BE BODY LOAD
          !
          case (fbem_bl_uncoupled)
            do while (.true.)
              if (kword.gt.n_words) exit
              word=trim(fbem_extract_word(line,kword))
              if (trim(word).eq.'fx'.or.trim(word).eq.'f1') then
                kword=kword+1
                if (kword.gt.n_words) exit
                word=trim(fbem_extract_word(line,kword))
                be_bodyload(i)%ctype(1,1)=0
                read(word,*) be_bodyload(i)%cvalue_r(1,1,1)
                kword=kword+1
              else if (trim(word).eq.'fy'.or.trim(word).eq.'f2') then
                kword=kword+1
                if (kword.gt.n_words) exit
                word=trim(fbem_extract_word(line,kword))
                be_bodyload(i)%ctype(2,1)=0
                read(word,*) be_bodyload(i)%cvalue_r(2,1,1)
                kword=kword+1
              else if (trim(word).eq.'fz'.or.trim(word).eq.'f3') then
                if (problem%n.eq.2) call fbem_error_message(error_unit,0,trim(section_name),nl,trim(word)//' cannot be used in 2D.')
                kword=kword+1
                if (kword.gt.n_words) exit
                word=trim(fbem_extract_word(line,kword))
                be_bodyload(i)%ctype(3,1)=0
                read(word,*) be_bodyload(i)%cvalue_r(3,1,1)
                kword=kword+1
              else
                call fbem_error_message(error_unit,0,trim(section_name),nl,'syntax error in this section line.')
              end if
            end do

          case (fbem_bl_coupling_beam_tip)
            ! N/A

          case (fbem_bl_coupling_beam_line)
            ! N/A

          case (fbem_bl_coupling_shell_edge)
            ! N/A

          case (fbem_bl_coupling_shell_surface)
            ! N/A

        end select

      end if

      !
      ! Check if the command is not recognized
      !
      if (.not.exists) then
        call fbem_error_message(error_unit,0,trim(word),0,'unrecognized command section ['//trim(section_name)//']')
      end if
    end do

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  end if

end subroutine read_conditions_bem_bodyloads_mechanics_static
