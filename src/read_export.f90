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

subroutine read_export(input_fileunit)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_symmetry
  use fbem_data_structures
  use fbem_gmsh

  ! Problem variables module
  use problem_variables

  ! No implicit variables are allowed
  implicit none

  ! I/O
  integer                                 :: input_fileunit    ! Input file unit
  ! Local
  character(len=fbem_stdcharlen)          :: section_name      ! Name of the section
  logical                                 :: found             ! Logical variable for sections and keywords
  integer                                 :: i                 ! Counter
  character(len=fbem_filename_max_length) :: tmp_string        ! Temporary file name
  integer                                 :: output_fileunit
  character(len=fbem_filename_max_length) :: tmp_filename
  character(len=fbem_fmtstr)              :: fmtstr   ! String used for write format string


  ! Default settings
  ! ----------------

  ! Vectors scale factor
  vectors_scale_factor=1.0d0

  ! Integer export
  ! Automatic integer width
  write(fmt_integer,*) 'i', fbem_nchar_int(max(2,n_frequencies,region_eid_max,boundary_eid_max,element_eid_max,node_eid_max,internalpoint_eid_max))+1
  call fbem_trimall(fmt_integer)

  ! Real export
  ! Export with engineering single precision (all output files)
  fmt_real='en18.8e2'
  call fbem_export_gmsh_fmt_real(fmt_real)

  ! Complex notation
  complex_notation=2

  ! Default export files
  export_overwrite=.true.
  export_geo_mesh_data=.false.
  export_fun_mesh_data=.false.
  export_nso=.true.
  export_eso=.true.
  export_wsp=.false.
  export_tot=.false.
  tot_xm=0
  tot_apply_symmetry=.true.
  if (mesh_file_mode.eq.2) export_pos=.true.

  section_name='export'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  ! Find the "export" section in the file
  call fbem_search_section(input_fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    !
    ! Find "export_overwrite"
    !
    call fbem_search_section(input_fileunit,'export',found)
    call fbem_search_keyword(input_fileunit,'export_overwrite','=',found)
    if (found) read(input_fileunit,*) export_overwrite
    !
    ! Find "export_geo_mesh_data"
    !
    call fbem_search_section(input_fileunit,'export',found)
    call fbem_search_keyword(input_fileunit,'export_geo_mesh_data','=',found)
    if (found) read(input_fileunit,*) export_geo_mesh_data
    !
    ! Find "export_fun_mesh_data"
    !
    call fbem_search_section(input_fileunit,'export',found)
    call fbem_search_keyword(input_fileunit,'export_fun_mesh_data','=',found)
    if (found) read(input_fileunit,*) export_fun_mesh_data
    !
    ! Find "export_nso"
    !
    call fbem_search_section(input_fileunit,'export',found)
    call fbem_search_keyword(input_fileunit,'export_nso','=',found)
    if (found) read(input_fileunit,*) export_nso
    !
    ! Find "export_eso"
    !
    call fbem_search_section(input_fileunit,'export',found)
    call fbem_search_keyword(input_fileunit,'export_eso','=',found)
    if (found) read(input_fileunit,*) export_eso
    !
    ! Find "export_wsp"
    !
    call fbem_search_section(input_fileunit,'export',found)
    call fbem_search_keyword(input_fileunit,'export_wsp','=',found)
    if (found) read(input_fileunit,*) export_wsp
    !
    ! Find "export_tot"
    !
    call fbem_search_section(input_fileunit,'export',found)
    call fbem_search_keyword(input_fileunit,'export_tot','=',found)
    if (found) read(input_fileunit,*) export_tot
    !
    ! Find "export_pos"
    !
    call fbem_search_section(input_fileunit,'export',found)
    call fbem_search_keyword(input_fileunit,'export_pos','=',found)
    if (found) read(input_fileunit,*) export_pos
    !
    ! Find "tot_xm"
    !
    call fbem_search_section(input_fileunit,'export',found)
    call fbem_search_keyword(input_fileunit,'tot_xm','=',found)
    if (found) read(input_fileunit,*) tot_xm
    if ((tot_xm.ne.0).and.(tot_xm.ne.1)) then
      call fbem_error_message(error_unit,0,'[export]',tot_xm,'invalid value of tot_xm keyword')
    end if
    !
    ! Find "tot_apply_symmetry"
    !
    call fbem_search_section(input_fileunit,'export',found)
    call fbem_search_keyword(input_fileunit,'tot_apply_symmetry','=',found)
    if (found) read(input_fileunit,*) tot_apply_symmetry
    !
    ! Find "real_format"
    !
    call fbem_search_section(input_fileunit,'export',found)
    call fbem_search_keyword(input_fileunit,'real_format','=',found)
    ! If "real_format" is present, read from the indicated file
    if (found) then
      ! Read the temporary mode
      read(input_fileunit,*) fmt_real
      if (trim(fmt_real).eq.'sci_double') fmt_real='e25.16e3'
      if (trim(fmt_real).eq.'sci_simple') fmt_real='e16.8e2'
      if (trim(fmt_real).eq.'sci_less'  ) fmt_real='e11.3e2'
      if (trim(fmt_real).eq.'eng_double') fmt_real='en27.16e3'
      if (trim(fmt_real).eq.'eng_simple') fmt_real='en18.8e2'
      if (trim(fmt_real).eq.'eng_less'  ) fmt_real='en13.3e2'
      if (trim(fmt_real).eq.'auto'      ) fmt_real='en27.16e3'
      call fbem_trimall(fmt_real)
      ! Message
      if (verbose_level.ge.3) then
        write(output_unit,'(3x,a)') 'Real format testing:'
        write(output_unit,'('//fmt_real//')') -9999.d0
      end if
    end if
    !
    ! Find "real_format_pos"
    !
    call fbem_search_section(input_fileunit,'export',found)
    call fbem_search_keyword(input_fileunit,'real_format_pos','=',found)
    ! If "real_format" is present, read from the indicated file
    if (found) then
      ! Read the temporary mode
      read(input_fileunit,*) fmt_real_pos
      if (trim(fmt_real_pos).eq.'sci_double') fmt_real_pos='e25.16e3'
      if (trim(fmt_real_pos).eq.'sci_simple') fmt_real_pos='e16.8e2'
      if (trim(fmt_real_pos).eq.'sci_less'  ) fmt_real_pos='e11.3e2'
      if (trim(fmt_real_pos).eq.'eng_double') fmt_real_pos='en27.16e3'
      if (trim(fmt_real_pos).eq.'eng_simple') fmt_real_pos='en18.8e2'
      if (trim(fmt_real_pos).eq.'eng_less'  ) fmt_real_pos='en13.3e2'
      if (trim(fmt_real_pos).eq.'auto'      ) fmt_real_pos='en27.16e3'
      call fbem_trimall(fmt_real_pos)
      call fbem_export_gmsh_fmt_real(fmt_real_pos)
      ! Message
      if (verbose_level.ge.3) then
        write(output_unit,'(3x,a)') 'Real format testing:'
        write(output_unit,'('//fmt_real_pos//')') -9999.d0
      end if
    else
      call fbem_export_gmsh_fmt_real(fmt_real)
    end if
    !
    ! Find "integer_format"
    !
    call fbem_search_section(input_fileunit,'export',found)
    call fbem_search_keyword(input_fileunit,'integer_format','=',found)
    ! If "integer_format" is present, read from the indicated file
    if (found) then
      ! Read the temporary mode
      read(input_fileunit,*) fmt_integer
      if (trim(fmt_integer).eq.'max' ) fmt_integer='i11'
      if (trim(fmt_integer).eq.'auto') then
        write(fmt_integer,*) 'i', fbem_nchar_int(max(n_frequencies,region_eid_max,boundary_eid_max,element_eid_max,node_eid_max,internalpoint_eid_max))+1
      end if
      call fbem_trimall(fmt_integer)
      ! Message
      if (verbose_level.ge.3) then
        write(output_unit,'(3x,a)') 'Integer format testing:'
        write(output_unit,'('//fmt_integer//')') -9999
      end if
    end if
    !
    ! Find "complex_notation"
    !
    call fbem_search_section(input_fileunit,'export',found)
    call fbem_search_keyword(input_fileunit,'complex_notation','=',found)
    ! If "complex_notation" is present, read from the indicated file
    if (found) then
      ! Read the temporary the mode
      read(input_fileunit,*) tmp_string
      complex_notation=0
      if (trim(tmp_string).eq.'polar'     ) complex_notation=1
      if (trim(tmp_string).eq.'cartesian' ) complex_notation=2
      ! Check
      if (complex_notation.eq.0) then
        call fbem_error_message(error_unit,0,'[export]',0,'the indicated complex notation is not valid.')
      end if
    end if
    !
    ! Find "vectors_scale_factor"
    !
    call fbem_search_section(input_fileunit,'export',found)
    call fbem_search_keyword(input_fileunit,'vectors_scale_factor','=',found)
    ! If "vectors_scale_factor" is present, read from the indicated file
    if (found) then
      ! Read the temporary the mode
      read(input_fileunit,*) vectors_scale_factor
      ! Check
      if (vectors_scale_factor.le.0) call fbem_error_message(error_unit,0,'[export]',0,'the indicated scale factor is not valid.')
    end if
    !
    ! Find "nso_nodes"
    !
    call fbem_search_section(input_fileunit,'export',found)
    call fbem_search_keyword(input_fileunit,'nso_nodes','=',found)
    ! Read if found
    if (found) then
      read(input_fileunit,*) nso_nodes
      if (nso_nodes.le.0) then
        if (verbose_level.ge.3) then
          write(output_unit,'(3x,a)') 'Results will be printed for all nodes if export_nso is True'
        end if
      else if (nso_nodes.gt.0) then
        ! Allocate
        allocate(nso_nodes_export(nso_nodes))
        ! Set the cursor again at the correct position
        call fbem_search_section(input_fileunit,'export',found)
        call fbem_search_keyword(input_fileunit,'nso_nodes','=',found)
        read(input_fileunit,*) nso_nodes, (nso_nodes_export(i),i=1,nso_nodes)
        ! Sort the array in ascend order
        call fbem_quicksort(1,nso_nodes,nso_nodes,nso_nodes_export)
        ! Check if any value is repeated
        do i=2,nso_nodes
          if (nso_nodes_export(i).eq.nso_nodes_export(i-1)) then
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                'there are repeated nodes for export')
          end if
        end do
        ! Set for export only nodes in the list
        do i=1,n_nodes
          node(i)%export=.false.
        end do
          node(nso_nodes_export)%export=.true.
        ! Write
        if (verbose_level.ge.3) then
          write(fmtstr,*) '(2x,a,i2,',nso_nodes,'i3)'
          call fbem_trim2b(fmtstr)
          write(output_unit,fmtstr) 'nso_nodes = ', nso_nodes, (nso_nodes_export(i),i=1,nso_nodes)
        end if
      end if
    end if

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  end if

end subroutine read_export
