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
!! <b> Subroutine that reads the problem settings section. </b>
!!
!!>
!!+-------------------------------+------------------------------+-----------------------------------------------------------------+
!!|                      VARIABLE |                       VALUES |                                                     DESCRIPTION |
!!+-------------------------------+------------------------------+-----------------------------------------------------------------+
!!|            qsi_relative_error |  [1.d-15,1.d-3] (def. 1.d-6) |       relative error in quasi-singular integration when solving |
!!|       qsi_post_relative_error |  [1.d-15,1.d-3] (def. 1.d-6) |   rel. error in quasi-singular integration when post-processing |
!!|                    qsi_ns_max |            [0,inf] (def. 16) |                                  maximum number of subdivisions |
!!|                    precalsets |      n gln_1 gln_2 ... gln_n |    number of precalculated datasets (n) and number of points of |
!!|                               |                              |       the Gauss-Legendre 1D quadrature (gln_k) for each dataset |
!!|           geometric_tolerance |     >1.0d-12 (default 1.d-9) |                                             geometric tolerance |
!!|            collapse_nodal_pos |           T or F (default T) |         move each node to the mean position of its common nodes |
!!|                   lse_scaling |           T or F (default F) |               perform scaling of the linear system of equations |
!!|                 lse_condition |           T or F (default F) | estimate the condition number of the linear system of equations |
!!|                    lse_refine |           T or F (default F) |   perform solution refinement of the linear system of equations |
!!|                 mesh_file_mode|               <0,1,2> "file" | read mesh from auxiliary file if 1 (native format), if 2 (gmsh) |
!!|              frequencies_file |          "file" (default "") |                                                                 |
!!+-------------------------------+------------------------------+-----------------------------------------------------------------+
!!<
subroutine read_settings(fileunit)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  implicit none
  ! I/O
  integer                                 :: fileunit ! Unit of the file to read from
  ! Local
  logical                                 :: found    ! Logical variable for sections and keywords
  integer                                 :: i        ! Counter
  character(len=fbem_fmtstr)              :: fmtstr   ! String used for write format string
  logical                                 :: check
  integer                                 :: dbfileunit, dbfileiostat

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section [settings]')
  call fbem_search_section(fileunit,'settings',found)
  if (found.and.(verbose_level.ge.2)) call fbem_timestamp_w_message(output_unit,2,'START reading section [settings]')

  ! ------------------
  ! qsi_relative_error
  ! ------------------

  call fbem_search_section(fileunit,'settings',found)
  if (found) call fbem_search_keyword(fileunit,'qsi_relative_error','=',found)
  if (found) then
    read(fileunit,*) qsi_relative_error
    if ((qsi_relative_error.lt.1.0d-15).or.(qsi_relative_error.gt.1.0d-3)) then
      call fbem_error_message(error_unit,0,'qsi_relative_error',0,'this variable must be between 1.E-15 and 1.E-3')
    end if
    if (verbose_level.ge.3) write(output_unit,'(2x,a,e10.3)') 'qsi_relative_error =', qsi_relative_error
  else
    qsi_relative_error=1.0d-6
  end if
  call fbem_qs_calculate_parameters(qsi_relative_error,qsi_parameters)

  ! -----------------------
  ! qsi_post_relative_error
  ! -----------------------

  call fbem_search_section(fileunit,'settings',found)
  if (found) call fbem_search_keyword(fileunit,'qsi_post_relative_error','=',found)
  if (found) then
    read(fileunit,*) qsi_post_relative_error
    if ((qsi_post_relative_error.lt.1.0d-15).or.(qsi_post_relative_error.gt.1.0d-3)) then
      call fbem_error_message(error_unit,0,'qsi_post_relative_error',0,'this variable must be between 1.E-15 and 1.E-3')
    end if
    if (verbose_level.ge.3) write(output_unit,'(2x,a,e10.3)') 'qsi_post_relative_error =', qsi_post_relative_error
  else
    qsi_post_relative_error=1.0d-6
  end if
  call fbem_qs_calculate_parameters(qsi_post_relative_error,qsi_post_parameters)

  ! ----------
  ! qsi_ns_max
  ! ----------

  call fbem_search_section(fileunit,'settings',found)
  if (found) call fbem_search_keyword(fileunit,'qsi_ns_max','=',found)
  if (found) then
    read(fileunit,*) qsi_ns_max
    if (qsi_ns_max.lt.0) call fbem_error_message(error_unit,0,'qsi_ns_max',0,'this variable must be >=0')
    if (verbose_level.ge.3) write(output_unit,'(2x,a,i3)') 'qsi_ns_max =', qsi_ns_max
  else
    qsi_ns_max=16
  end if

  ! ----------
  ! precalsets
  ! ----------

  call fbem_search_section(fileunit,'settings',found)
  if (found) call fbem_search_keyword(fileunit,'precalsets','=',found)
  ! Read if found
  if (found) then
    read(fileunit,*) n_precalsets
    if ((n_precalsets.lt.2).or.(n_precalsets.gt.30)) then
      call fbem_error_message(error_unit,0,'precalsets',0,'the number of precalsets must be between 2 and 30')
    end if
    ! Allocate
    allocate(precalset_gln(n_precalsets))
    ! Set the cursor again at the correct position
    call fbem_search_section(fileunit,'settings',found)
    call fbem_search_keyword(fileunit,'precalsets','=',found)
    read(fileunit,*) n_precalsets, (precalset_gln(i),i=1,n_precalsets)
    if ((minval(precalset_gln).lt.2).or.(maxval(precalset_gln).gt.30)) then
      call fbem_error_message(error_unit,0,'precalsets',0,&
                              'the number of 1D Gauss-Legendre points in each precalset must be between 2 and 30')
    end if
    ! Sort the array in ascend order
    call fbem_quicksort(1,n_precalsets,n_precalsets,precalset_gln)
    ! Check if any value is repeated
    do i=2,n_precalsets
      if (precalset_gln(i).eq.precalset_gln(i-1)) then
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                'there are repeated 1D Gauss-Legendre points in the list of precalsets')
      end if
    end do
    ! Write
    if (verbose_level.ge.3) then
      write(fmtstr,*) '(2x,a,i2,',n_precalsets,'i3)'
      call fbem_trim2b(fmtstr)
      write(output_unit,fmtstr) 'precalsets = ', n_precalsets, (precalset_gln(i),i=1,n_precalsets)
    end if
  ! Default
  else
    n_precalsets=9
    allocate(precalset_gln(9))
    precalset_gln=(/2,3,4,5,6,7,8,9,10/)
  end if
  ! Save the maximum and minimum of precalset_gln()
  precalset_gln_min=precalset_gln(1)
  precalset_gln_max=precalset_gln(n_precalsets)

  ! -------------------
  ! geometric_tolerance
  ! -------------------

  call fbem_search_section(fileunit,'settings',found)
  if (found) call fbem_search_keyword(fileunit,'geometric_tolerance','=',found)
  if (found) then
    read(fileunit,*) geometric_tolerance
    if (geometric_tolerance.lt.0.0d0) then
      call fbem_error_message(error_unit,1,'geometric_tolerance',0,'the geometric tolerance must be positive')
    end if
    if (verbose_level.ge.3) write(output_unit,'(2x,a,e10.3)') 'geometric_tolerance =', geometric_tolerance
  else
    geometric_tolerance=1.d-6
  end if

  ! ------------------
  ! collapse_nodal_pos
  ! ------------------

  call fbem_search_section(fileunit,'settings',found)
  if (found) call fbem_search_keyword(fileunit,'collapse_nodal_pos','=',found)
  if (found) then
    read(fileunit,*) collapse_nodal_pos
    if (verbose_level.ge.3) write(output_unit,'(2x,a,l)') 'collapse_nodal_pos = ', collapse_nodal_pos
  else
    collapse_nodal_pos=.true.
  end if

  ! -----------
  ! lse_scaling
  ! -----------

  call fbem_search_section(fileunit,'settings',found)
  if (found) call fbem_search_keyword(fileunit,'lse_scaling','=',found)
  if (found) then
    read(fileunit,*) lse_scaling
    if (verbose_level.ge.3) write(output_unit,'(2x,a,l)') 'lse_scaling = ', lse_scaling
  else
    lse_scaling=.false.
  end if

  ! -------------
  ! lse_condition
  ! -------------

  call fbem_search_section(fileunit,'settings',found)
  if (found) call fbem_search_keyword(fileunit,'lse_condition','=',found)
  if (found) then
    read(fileunit,*) lse_condition
    if (verbose_level.ge.3) write(output_unit,'(2x,a,l)') 'lse_condition = ', lse_condition
  else
    lse_condition=.false.
  end if

  ! ----------
  ! lse_refine
  ! ----------

  call fbem_search_section(fileunit,'settings',found)
  if (found) call fbem_search_keyword(fileunit,'lse_refine','=',found)
  if (found) then
    read(fileunit,*) lse_refine
    if (verbose_level.ge.3) write(output_unit,'(2x,a,l)') 'lse_refine = ', lse_refine
  else
    lse_refine=.false.
  end if

  ! --------------
  ! mesh_file_mode
  ! --------------

  call fbem_search_section(fileunit,'settings',found)
  if (found) call fbem_search_keyword(fileunit,'mesh_file_mode','=',found)
  if (found) then
    read(fileunit,*) mesh_file_mode, mesh_filename
    call fbem_trimall(mesh_filename)
    ! If path to the same directory, remove it.
    if (mesh_filename(1:2).eq.'./') then
      mesh_filename=trim(mesh_filename(3:len_trim(mesh_filename)))
    end if
    ! Not yet allowed the full syntax to navigate between folders
    if (mesh_filename(1:2).eq.'..') then
      call fbem_error_message(error_unit,1,'mesh_filename',0,'wrong path to the file')
    end if
    ! Add path from the input file
    if (mesh_filename(1:1).ne.'/') then
      mesh_filename=trim(pwd)//trim(mesh_filename)
    end if
    call fbem_trimall(mesh_filename)
    ! Check
    if ((mesh_file_mode.lt.0).or.(mesh_file_mode.gt.2)) then
      call fbem_error_message(error_unit,1,'mesh_file_mode',0,'wrong type of mesh mode')
    end if
    if (verbose_level.ge.3) then
      write(output_unit,'(2x,a,i2)') 'mesh_file_mode = ', mesh_file_mode
      write(output_unit,'(2x,a,a)') 'mesh_filename = ', trim(mesh_filename)
    end if
  else
    mesh_file_mode=0
    mesh_filename=''
  end if

  ! ----------------
  ! frequencies_file
  ! ----------------

  call fbem_search_section(fileunit,'settings',found)
  if (found) call fbem_search_keyword(fileunit,'frequencies_file','=',found)
  if (found) then
    read(fileunit,*) frequencies_filename
    call fbem_trimall(frequencies_filename)
    ! If path to the same directory, remove it.
    if (frequencies_filename(1:2).eq.'./') then
      frequencies_filename=trim(frequencies_filename(3:len_trim(frequencies_filename)))
    end if
    ! Not yet allowed the full syntax to navigate between folders
    if (frequencies_filename(1:2).eq.'..') then
      call fbem_error_message(error_unit,1,'frequencies_filename',0,'wrong path to the file')
    end if
    ! Add path from the input file
    if (frequencies_filename(1:1).ne.'/') then
      frequencies_filename=trim(pwd)//trim(frequencies_filename)
    end if
    call fbem_trimall(frequencies_filename)
    if (verbose_level.ge.3) then
      write(output_unit,'(2x,a,a)') 'frequencies_filename = ', trim(frequencies_filename)
    end if
  else
    frequencies_filename=''
  end if

  call fbem_search_section(fileunit,'settings',found)
  if (found.and.(verbose_level.ge.2)) call fbem_timestamp_w_message(output_unit,2,'END reading section [settings]')

end subroutine read_settings
