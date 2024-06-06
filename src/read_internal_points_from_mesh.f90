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

!! @author Jacob David Rodriguez Bordon (jacobdavid.rodriguezbordon@ulpgc.es)
!!
!! @version 2.0
!!
!! <b> Read and create BEM internal points from a given mesh. </b>

subroutine read_internal_points_from_mesh(input_fileunit)
  ! Fortran 2003 intrinsic module
  use iso_fortran_env
  ! fbem modules
  use fbem_numerical
  use fbem_geometry
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_symmetry
  use fbem_data_structures
  use fbem_quad_rules
  ! Problem variables module
  use problem_variables
  ! No implicit variables are allowed
  implicit none
  ! I/O
  integer                                 :: input_fileunit    ! Input file unit
  ! Local
  character(len=fbem_stdcharlen)          :: section_name      ! Name of the section
  logical                                 :: found             ! Logical variable for sections and keywords
  character(len=fbem_filename_max_length) :: mesh_file_name    ! Mesh file name
  character(len=fbem_stdcharlen)          :: mesh_file_format  ! Mesh file format
  integer                                 :: i, j, j1, j2, k, sn
  integer                                 :: sp, sr, kp, kr
  type(fbem_mesh)                         :: internalpoints_mesh
  integer                                 :: n_additional_internalpoints
  type(fbem_internalpoint), allocatable   :: internalpoint_tmp(:)
  integer, allocatable                    :: internalpoint_iid_tmp(:)
  integer                                 :: iid, eid
  integer                                 :: etype
  real(kind=real64), allocatable          :: xi(:)
  real(kind=real64)                       :: centroid_part(problem%n), area_part


  ! Locate the section
  section_name='internal points from mesh'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(input_fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! -----------------------
    ! READ INTERNAL MESH FILE
    ! -----------------------

    !
    ! TO-DO: save mesh location in order to export the results to a *.pos there
    !

    call fbem_search_section(input_fileunit,section_name,found)
    if (found) call fbem_search_keyword(input_fileunit,'mesh_file','=',found)
    if (found) then
      read(input_fileunit,*) mesh_file_format, mesh_file_name
      call fbem_trim(mesh_file_name)
      if (.not.fbem_file_exists(mesh_file_name)) then
        write(output_unit,'(a89)') 'Mesh file does not exist ([internal points from mesh] : mesh_file), check the given path.'
        write(output_unit,*)
      end if
      if (fbem_path_is_relative(mesh_file_name)) then
        mesh_file_name=trim(input_filedir)//trim(mesh_file_name)
      end if
      if (verbose_level.ge.3) then
        write(output_unit,'(2x,a,a)') 'mesh_file_format = ', mesh_file_format
        write(output_unit,'(2x,a,a)') 'mesh_file_name = ', trim(mesh_file_name)
      end if
      ! Read the mesh from the file
      call internalpoints_mesh%read_from_file(problem%n,geometric_tolerance,mesh_file_name,mesh_file_format)
    end if

    ! ----------------------------------------------------------------
    ! READ THE CONNECTIVITY BETWEEN INTERNAL ELEMENTS PARTS TO REGIONS
    ! ----------------------------------------------------------------
    !
    ! region_association : <n_parts>
    ! <internal element mesh part id.> <region id.>
    ! ...

    call fbem_search_section(input_fileunit,section_name,found)
    if (found) call fbem_search_keyword(input_fileunit,'region_association',':',found)
    if (found) then
      ! For each part of the internal elements mesh, read the region where they are.
      read(input_fileunit,*) j
      if (j.ne.internalpoints_mesh%n_parts) then
        call fbem_error_message(error_unit,0,section_name,0,'the number of parts must be the same as in the internal elements mesh')
      end if
      do j=1,internalpoints_mesh%n_parts
        read(input_fileunit,*) sp, sr
        ! Check that sp is a valid part in the internal elements mesh, and then obtain its iid and save it to kp
        if ((sp.ge.internalpoints_mesh%part_eid_min).and.(sp.le.internalpoints_mesh%part_eid_max)) then
          if (internalpoints_mesh%part_iid(sp).eq.0) then
            call fbem_error_message(error_unit,0,section_name,sp,'the indicated part in the internal elements mesh does not exist')
          else
            kp=internalpoints_mesh%part_iid(sp)
          end if
        else
          call fbem_error_message(error_unit,0,section_name,sp,'the indicated part in the internal elements mesh does not exist')
        end if
        ! Check that sr is a valid region, obtain its iid, and save it to the internal element part member "entity"
        if ((sr.ge.region_eid_min).and.(sr.le.region_eid_max)) then
          if (region_iid(sr).eq.0) then
            call fbem_error_message(error_unit,0,section_name,sr,'the indicated region does not exist')
          else
            kr=region_iid(sr)
            internalpoints_mesh%part(kp)%entity=kr
          end if
        else
          call fbem_error_message(error_unit,0,section_name,sr,'the indicated region does not exist')
        end if
      end do
    else
      call fbem_error_message(error_unit,0,section_name,0,'the list "region_association" must be defined')
    end if

    ! -----------------------------------
    ! GENERATE ADDITIONAL INTERNAL POINTS
    ! -----------------------------------

    if (n_internalpoints.gt.0) then
      !
      ! Copy the internal points to temporary variables
      !
      allocate (internalpoint_tmp(n_internalpoints))
      allocate (internalpoint_iid_tmp(internalpoint_eid_min:internalpoint_eid_max))
      do i=1,n_internalpoints
        allocate (internalpoint_tmp(i)%x(problem%n))
        internalpoint_tmp(i)%id=internalpoint(i)%id
        internalpoint_tmp(i)%x=internalpoint(i)%x
        internalpoint_tmp(i)%region=internalpoint(i)%region
        internalpoint_tmp(i)%export=internalpoint(i)%export
      end do
      internalpoint_iid_tmp=internalpoint_iid
      !
      ! Deallocate the internal points
      !
      do i=1,n_internalpoints
        deallocate (internalpoint(i)%x)
      end do
      deallocate (internalpoint)
      deallocate (internalpoint_iid)
    end if

    !
    ! Count the number of additional internal points necessary
    !
    n_additional_internalpoints=internalpoints_mesh%n_nodes
    !
    ! Allocate the internal points list including the additional internal points, and save the internal points in it.
    !
    allocate (internalpoint(n_internalpoints+n_additional_internalpoints))
    do i=1,n_internalpoints
      allocate (internalpoint(i)%x(problem%n))
      internalpoint(i)%id=internalpoint_tmp(i)%id
      internalpoint(i)%x=internalpoint_tmp(i)%x
      internalpoint(i)%region=internalpoint_tmp(i)%region
      internalpoint(i)%export=internalpoint_tmp(i)%export
    end do
    !
    ! Generate additional internal points
    !
    do i=1,internalpoints_mesh%n_nodes
      allocate (internalpoint(n_internalpoints+i)%x(problem%n))
      internalpoint(n_internalpoints+i)%id=internalpoints_mesh%node(i)%id
      internalpoint(n_internalpoints+i)%x =internalpoints_mesh%node(i)%x
      j=internalpoints_mesh%node(i)%element(1)
      internalpoint(n_internalpoints+i)%region=internalpoints_mesh%part(internalpoints_mesh%element(j)%part)%entity
      internalpoint(n_internalpoints+i)%export=.true.
    end do
    n_internalpoints=n_internalpoints+n_additional_internalpoints

    ! =================================
    ! CHECK AND BUILD NODES IDENTIFIERS
    ! =================================

    internalpoint_eid_min=internalpoint(1)%id
    internalpoint_eid_max=internalpoint(1)%id
    do i=2,n_internalpoints
      if (internalpoint(i)%id.lt.internalpoint_eid_min) internalpoint_eid_min=internalpoint(i)%id
      if (internalpoint(i)%id.gt.internalpoint_eid_max) internalpoint_eid_max=internalpoint(i)%id
    end do
    if (internalpoint_eid_min.le.0) then
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                               'the internal points identifiers must be >0')
    end if
    if (allocated(internalpoint_iid)) deallocate(internalpoint_iid)
    allocate (internalpoint_iid(internalpoint_eid_min:internalpoint_eid_max))
    internalpoint_iid=0
    do i=1,n_internalpoints
      if (internalpoint_iid(internalpoint(i)%id).ne.0) then
        call fbem_error_message(error_unit,0,'internal point',internalpoint(i)%id,'is repeated')
      else
        internalpoint_iid(internalpoint(i)%id)=i
      end if
    end do

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  end if

end subroutine read_internal_points_from_mesh
