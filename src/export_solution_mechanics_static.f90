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

subroutine export_solution_mechanics_static

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_quad_rules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_geometry
  use fbem_symmetry

  ! Module of problem variables
  use problem_variables

  ! No implicit variables allowed
  implicit none

  ! Local variables
  integer                                 :: output_fileunit
  character(len=fbem_filename_max_length) :: tmp_filename

  if (verbose_level.ge.1)  call fbem_timestamp_w_message(output_unit,2,'START exporting solutions')

  ! ============================================================================================================================
  ! EXPORT NODAL SOLUTIONS
  ! ============================================================================================================================

  !
  ! Export nodal solutions in native format (*.nso)
  !
  if (export_nso) then
    output_fileunit=fbem_get_valid_unit()
    write(tmp_filename,*) trim(output_filename), '.nso'
    call fbem_trim2b(tmp_filename)
    if (export_overwrite) then
      open(unit=output_fileunit,file=trim(tmp_filename),action='write',recl=fbem_file_record_length)
    else
      open(unit=output_fileunit,file=trim(tmp_filename),access='append',recl=fbem_file_record_length)
    end if
    if (verbose_level.ge.2) then
      call fbem_timestamp_message(output_unit,2)
      write(output_unit,'(1x,a,a,a)') 'START exporting solution to ', trim(tmp_filename), '"'
    end if
    call export_solution_mechanics_static_nso(output_fileunit)
    close(unit=output_fileunit)
    if (verbose_level.ge.2) then
      call fbem_timestamp_message(output_unit,2)
      write(output_unit,'(1x,a,a,a)') 'END exporting solution to ', trim(tmp_filename), '"'
    end if
  end if

  ! ==========================================================================================================================
  ! EXPORT ELEMENT NODE SOLUTIONS
  ! ==========================================================================================================================

  !
  ! Export element node solutions in native format (*.eso)
  !
  if (export_eso) then
    output_fileunit=fbem_get_valid_unit()
    write(tmp_filename,*) trim(output_filename), '.eso'
    call fbem_trim2b(tmp_filename)
    if (export_overwrite) then
      open(unit=output_fileunit,file=trim(tmp_filename),action='write',recl=fbem_file_record_length)
    else
      open(unit=output_fileunit,file=trim(tmp_filename),access='append',recl=fbem_file_record_length)
    end if
    if (verbose_level.ge.2) then
      call fbem_timestamp_message(output_unit,2)
      write(output_unit,'(1x,a,a,a)') 'START exporting solution to ', trim(tmp_filename), '"'
    end if
    call export_solution_mechanics_static_eso(output_fileunit)
    close(unit=output_fileunit)
    if (verbose_level.ge.2) then
      call fbem_timestamp_message(output_unit,2)
      write(output_unit,'(1x,a,a,a)') 'END exporting solution to ', trim(tmp_filename), '"'
    end if
  end if

  ! ============================================================================================================================
  ! BOUNDARY RESULTANTS FILE
  ! ============================================================================================================================

  !
  ! Export boundary resultants in native format (*.tot)
  !
  if (export_eso) then
    output_fileunit=fbem_get_valid_unit()
    write(tmp_filename,*) trim(output_filename), '.tot'
    call fbem_trim2b(tmp_filename)
    if (export_overwrite) then
      open(unit=output_fileunit,file=trim(tmp_filename),action='write',recl=fbem_file_record_length)
    else
      open(unit=output_fileunit,file=trim(tmp_filename),access='append',recl=fbem_file_record_length)
    end if
    if (verbose_level.ge.2) then
      call fbem_timestamp_message(output_unit,2)
      write(output_unit,'(1x,a,a,a)') 'START exporting solution to ', trim(tmp_filename), '"'
    end if
    call export_solution_mechanics_static_tot(output_fileunit)
    close(unit=output_fileunit)
    if (verbose_level.ge.2) then
      call fbem_timestamp_message(output_unit,2)
      write(output_unit,'(1x,a,a,a)') 'END exporting solution to ', trim(tmp_filename), '"'
    end if
  end if

  ! ================================================================================================================================
  ! EXPORT SOLUTIONS IN GMSH FILE FORMAT
  ! ================================================================================================================================

  if (export_pos) then
    !
    ! Write the mesh in the *.pos file
    !
    ! Copy from input file
    select case (mesh_file_mode)
      ! From main input file (multifebe format)
      case (0)
        write(tmp_filename,*) trim(output_filename), '.pos'
        call fbem_trim2b(tmp_filename)
        call fbem_convert_mesh_file_format(problem%n,input_filename,'multifebe',tmp_filename,'gmsh')
      ! From auxiliary file (multifebe format)
      case (1)
        write(tmp_filename,*) trim(output_filename), '.pos'
        call fbem_trim2b(tmp_filename)
        call fbem_convert_mesh_file_format(problem%n,mesh_filename,'multifebe',tmp_filename,'gmsh')
      ! From auxiliary file (gmsh format)
      case (2)
        write(tmp_filename,*) trim(output_filename), '.pos'
        call fbem_trim2b(tmp_filename)
        call fbem_convert_mesh_file_format(problem%n,mesh_filename,'gmsh',tmp_filename,'gmsh')
    end select
    !
    ! Write the Node Data sections
    !
    output_fileunit=fbem_get_valid_unit()
    open(unit=output_fileunit,file=trim(tmp_filename),access='append',recl=fbem_file_record_length)
    if (verbose_level.ge.2) then
      call fbem_timestamp_message(output_unit,2)
      write(output_unit,'(1x,a,a,a)') 'START exporting solution to ', trim(tmp_filename), '"'
    end if
    call export_solution_mechanics_static_gmsh(output_fileunit)
    close(unit=output_fileunit)
    if (verbose_level.ge.2) then
      call fbem_timestamp_message(output_unit,2)
      write(output_unit,'(1x,a,a,a)') 'END exporting solution to ', trim(tmp_filename), '"'
    end if
  end if

  if (verbose_level.ge.1)  call fbem_timestamp_w_message(output_unit,2,'END exporting solutions')

end subroutine export_solution_mechanics_static
