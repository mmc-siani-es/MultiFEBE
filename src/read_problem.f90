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
!! <b> Subroutine that reads the problem definition section. </b>
!!
!!>
!!+-------------------------------+------------------------------+-----------------------------------------------------------------+
!!|                      VARIABLE |                       VALUES |                                                     DESCRIPTION |
!!+-------------------------------+------------------------------+-----------------------------------------------------------------+
!!|                          type |                      laplace |                                                 Laplace problem |
!!|                               |                    mechanics |                                               Mechanics problem |
!!|                             n |                           2D |                                         two-dimensional problem |
!!|                               |                           3D |                                       three-dimensional problem |
!!|                       subtype |                 plane strain |                            two-dimensional plane strain problem |
!!|                               |                 plane stress |                            two-dimensional plane stress problem |
!!|                      analysis |                       static |                                                 static analysis |
!!|                               |                     harmonic |                                               harmonic analysis |
!!|                   sensitivity |         <T or F> (default F) |                     True if a sensitivity analysis is performed |
!!+-------------------------------+------------------------------+-----------------------------------------------------------------+
!!<
subroutine read_problem(fileunit)

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
  integer                                 :: fileunit          !! Unit of the file to read from
  ! Local
  character(len=fbem_stdcharlen)          :: section_name      ! Name of the section
  logical                                 :: found             ! Logical variable for sections and keywords
  character(len=fbem_file_record_length)  :: line              ! Line

  section_name='problem'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! Read the variable "type"
    call fbem_search_section(fileunit,section_name,found)
    call fbem_search_keyword(fileunit,'type','=',found)
    if (found) then
      read(fileunit,'(a)') line
      call fbem_trim2b(line)
      problem%type=0
      if (trim(line).eq.'mechanics') problem%type=fbem_mechanics
      if (trim(line).eq.'laplace'  ) problem%type=fbem_laplace
      if (problem%type.eq.0) then
        call fbem_error_message(error_unit,0,'[problem]',0,'variable "type" has not a valid value')
      end if
      if (verbose_level.ge.3) write(output_unit,'(2x,a,a)') 'type = ', trim(line)
    else
      call fbem_error_message(error_unit,0,'[problem]',0,'variable "type" is required')
    end if

    ! Read the variable "analysis"
    if (problem%type.eq.fbem_mechanics) then
      call fbem_search_section(fileunit,section_name,found)
      call fbem_search_keyword(fileunit,'analysis','=',found)
      if (found) then
        problem%analysis=0
        read(fileunit,'(a)') line
        call fbem_trim2b(line)
        if (trim(line).eq.'static'  ) problem%analysis=fbem_static
        if (trim(line).eq.'harmonic') problem%analysis=fbem_harmonic
        if (problem%analysis.eq.0) then
          call fbem_error_message(error_unit,0,'[problem]',0,'variable "analysis" has not a valid value')
        end if
        if (verbose_level.ge.3) write(output_unit,'(2x,a,a)') 'analysis = ', trim(line)
      else
        call fbem_error_message(error_unit,0,'[problem]',0,'variable "analysis" is required')
      end if
    else
      problem%analysis=0
    end if

    ! Read the variable "n"

    ! Si se metiesen problemas 2D axisimetricos, o 2.5D, habría que cambiar un poco esto....

    call fbem_search_section(fileunit,section_name,found)
    call fbem_search_keyword(fileunit,'n','=',found)
    if (found) then
      problem%n=0
      read(fileunit,'(a)') line
      call fbem_trim2b(line)
      if (trim(line).eq.'2' ) problem%n=2
      if (trim(line).eq.'2D') problem%n=2
      if (trim(line).eq.'3' ) problem%n=3
      if (trim(line).eq.'3D') problem%n=3
      if (problem%n.eq.0) then
        call fbem_error_message(error_unit,0,'[problem]',0,'variable "n" has not a valid value')
      end if
      if (verbose_level.ge.3) write(output_unit,'(2x,a,i1)') 'n = ', problem%n
    else
      call fbem_error_message(error_unit,0,'[problem]',0,'variable "n" is required')
    end if

    ! Read subtype if it is required
    if ((problem%type.eq.fbem_mechanics).and.(problem%n.eq.2)) then
      call fbem_search_section(fileunit,section_name,found)
      call fbem_search_keyword(fileunit,'subtype','=',found)
      if (found) then
        problem%subtype=0
        read(fileunit,'(a)') line
        call fbem_trim2b(line)
        if ((trim(line).eq.'plane strain').or.(trim(line).eq.'plane_strain')) problem%subtype=fbem_mechanics_plane_strain
        if ((trim(line).eq.'plane stress').or.(trim(line).eq.'plane_stress')) problem%subtype=fbem_mechanics_plane_stress
        if (problem%subtype.eq.0) then
          call fbem_error_message(error_unit,0,'[problem]',0,'variable "subtype" has not a valid value')
        end if
        if (verbose_level.ge.3) write(output_unit,'(2x,a,a)') 'analysis = ', trim(line)
      else
        problem%subtype=fbem_mechanics_plane_strain
        call fbem_warning_message(error_unit,0,'[problem]',0,'by default subtype=plane_strain')
      end if
    else
      problem%subtype=0
    end if

    ! Read the variable "description"
    call fbem_search_section(fileunit,section_name,found)
    call fbem_search_keyword(fileunit,'description','=',found)
    if (found) then
      read(fileunit,'(a)') line
      call fbem_trim2b(line)
      problem%description=trim(line)
      if (verbose_level.ge.3) write(output_unit,'(2x,a,a)') 'description = ', trim(line)
    else
      problem%description=''
    end if

    ! Read the variable "sensitivity"
    call fbem_search_section(fileunit,section_name,found)
    call fbem_search_keyword(fileunit,'sensitivity','=',found)
    if (found) then
      read(fileunit,*) problem%sensitivity
      if (verbose_level.ge.2) write(output_unit,'(2x,a,l)') 'sensitivity = ', problem%sensitivity

      if (problem%n.ne.2) stop 'Sensitivity analysis only in 2D'

    else
      problem%sensitivity=.false.
    end if

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section [problem]')

  else

    call fbem_error_message(error_unit,0,'['//trim(section_name)//']',0,'this section is required')

  end if

end subroutine read_problem
