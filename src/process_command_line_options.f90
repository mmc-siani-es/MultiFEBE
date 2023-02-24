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


subroutine process_command_line_options

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_string_handling

  ! Problem variables module
  use problem_variables

  ! A module created by Mark Gates to process program options
  use getopt_m

  ! No implicit variables are allowed
  implicit none

  integer                    :: i
  logical                    :: ok
  type(option_s)             :: opts(6)    ! Variable needed by getopt_m
  character(len=fbem_fmtstr) :: fmtstr     ! String used for write format string

  !
  ! Default values
  !
  input_filename=''
  output_filename=''
  max_memory=0
  verbose_level=1
  !
  ! Define options
  !
  opts(1)=option_s('input',.true.,'i')
  opts(2)=option_s('output',.true.,'o')
  opts(3)=option_s('maxmemory',.true.,'m')
  opts(4)=option_s('verbose',.true.,'b')
  opts(5)=option_s('version',.false.,'v')
  opts(6)=option_s('help',.false.,'h')
  !
  ! Process options
  !
  do
    select case (getopt("i:o:m:b:vh",opts))
      case (char(0))
        exit
      case ('i')
        input_filename=optarg
      case ('o')
        output_filename=optarg
      case ('m')
        read (optarg,*) max_memory
        max_memory=max_memory*1024_8**3
      case ('b')
        read (optarg,*) verbose_level
        if (verbose_level.lt.0) verbose_level=0
      case ('v')
        call version
        stop
      case ('h')
        call help
      case ('?')
        write(output_unit,'(a28)') 'Unknown command-line option.'
        write(output_unit,*)
        call help
    end select
  end do
  !
  ! Check options
  !
  ! CHECK INPUT FILE
  !
  if (len_trim(input_filename).eq.0) then
    write(output_unit,'(a24)') 'Input file is mandatory.'
    write(output_unit,*)
    call help
  else
    call fbem_trim(input_filename)
    if (.not.fbem_file_exists(input_filename)) then
      write(output_unit,'(a48)') 'Input file does not exist, check the given path.'
      write(output_unit,*)
      call help
    end if
    input_filedir = fbem_get_dirname(input_filename)
    call fbem_trim(input_filedir)
  end if
  !
  ! CHECK OUTPUT FILE
  !
  if (len_trim(output_filename).eq.0) then
    output_filename = input_filename
  else
    call fbem_trim(output_filename)
  end if
  !
  ! Print
  !
  if (verbose_level.ge.1) then
    ! Input file directory path
    if (len_trim(input_filedir).eq.0) then
      write(output_unit,'(a22)') 'Input file directory: '
    else
      write(fmtstr,*) '(a22,a',len_trim(input_filedir),')'
      call fbem_trimall(fmtstr)
      write(output_unit,fmtstr) 'Input file directory: ', trim(input_filedir)
    end if
    ! Input file path
    write(fmtstr,*) '(a22,a',len_trim(input_filename),')'
    call fbem_trimall(fmtstr)
    write(output_unit,fmtstr) 'Input file          : ', trim(input_filename)
    ! Output files base path
    write(fmtstr,*) '(a22,a',len_trim(output_filename),',a2)'
    call fbem_trimall(fmtstr)
    write(output_unit,fmtstr) 'Output files        : ', trim(output_filename), '.*'
    ! Verbose level
    write(fmtstr,*) '(a22,i',fbem_nchar_int(verbose_level),')'
    call fbem_trimall(fmtstr)
    write(output_unit,fmtstr) 'Verbose level       : ', verbose_level
  end if

end subroutine process_command_line_options
