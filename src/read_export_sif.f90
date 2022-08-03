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


subroutine read_export_sif(input_fileunit)
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
  integer                                 :: input_fileunit    ! Input file unit
  ! Local
  character(len=fbem_stdcharlen)          :: section_name      ! Name of the section
  logical                                 :: found             ! Logical variable for sections and keywords
  integer                                 :: i                 ! Counters
  integer                                 :: tmp_int           ! Temporary integer
  integer                                 :: tmp_int2          ! Temporary integer

  section_name='export sif'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(input_fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    export_sif=.true.
    read(input_fileunit,*) tmp_int
    do i=1,tmp_int
      read(input_fileunit,*) tmp_int2





      node(node_iid(tmp_int2))%export_sif=.true.

      ! Falta chequear aqui si en el nodo indicado se puede hacer eso......
      ! Y que es solo por ahora para 2D
      if (problem%n.eq.3) stop 'not implemented yet'



    end do

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  else

    export_sif=.false.
    node%export_sif=.false.

  end if

end subroutine read_export_sif
