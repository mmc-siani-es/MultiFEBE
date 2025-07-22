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
!! <b> Subroutine that reads the boundaries from a file. </b>
subroutine read_boundaries(fileunit)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! Local variables
  implicit none
  ! I/O
  integer                        :: fileunit          ! Unit of the file to read from
  ! Local
  character(len=fbem_stdcharlen) :: section_name      ! Name of the section
  logical                        :: found
  character(len=fbem_stdcharlen) :: tmp_class         ! Temporary variable to read the class of boundary
  integer                        :: i, j, k           ! Counters
  integer                        :: n_region_be_boundaries
  integer                        :: kr, kb, sb
  integer                        :: tmp_int           ! Temporary integer
  logical                        :: tmp_bool          ! Temporary integer

  ! Default value
  n_boundaries=0

  ! Return if not needed
  if (n_be_regions.eq.0) return

  ! Detect if there are BE boundaries in BE regions
  n_region_be_boundaries=0
  do i=1,n_regions
    if (region(i)%class.eq.fbem_be) then
      n_region_be_boundaries=n_region_be_boundaries+region(i)%n_boundaries
    end if
  end do

  ! Locate the section
  section_name='boundaries'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! Read the number of regions
    read(fileunit,*) n_boundaries

    ! Allocate and initialize
    if (n_boundaries.gt.0) then
      allocate (boundary(n_boundaries))
      do i=1,n_boundaries
        boundary(i)%id=0
        boundary(i)%name=''
        boundary(i)%class=0
        boundary(i)%coupling=0
        boundary(i)%n_regions=0
        boundary(i)%part=0
        boundary(i)%export=.true.
      end do
    else
      call fbem_error_message(error_unit,0,section_name,0,'the number of boundaries must be >0')
    end if

    ! Read boundary id, part, and class.
    do i=1,n_boundaries
      read(fileunit,*) boundary(i)%id, boundary(i)%part, tmp_class
      ! Check id
      if (boundary(i)%id.le.0) then
        call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'identifiers must be greater than 0')
      end if
      ! Assign the value to the class structure member
      boundary(i)%class=0
      if (trim(tmp_class).eq.'ordinary'  ) boundary(i)%class=fbem_boundary_class_ordinary
      if (trim(tmp_class).eq.'crack-like') boundary(i)%class=fbem_boundary_class_cracklike
      if (boundary(i)%class.eq.0) then
        call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'the boundary must be "ordinary" or "crack-like"')
      end if
    end do

    ! ======================================
    ! CHECK AND BUILD BOUNDARIES IDENTIFIERS
    ! ======================================

    boundary_eid_min=boundary(1)%id
    boundary_eid_max=boundary(1)%id
    do i=2,n_boundaries
      if (boundary(i)%id.lt.boundary_eid_min) boundary_eid_min=boundary(i)%id
      if (boundary(i)%id.gt.boundary_eid_max) boundary_eid_max=boundary(i)%id
    end do
    allocate (boundary_iid(boundary_eid_min:boundary_eid_max))
    boundary_iid=0
    do i=1,n_boundaries
      if (boundary_iid(boundary(i)%id).ne.0) then
        call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'is repeated')
      else
        boundary_iid(boundary(i)%id)=i
      end if
    end do

    ! =============================================================
    ! CHECK AND BUILD CONNECTIVITIES BETWEEN BOUNDARIES AND REGIONS
    ! =============================================================

    ! Convert eid to iid in region -> boundaries connectivity, and check if the indicated boundaries actually exist.
    do kr=1,n_regions
      if (region(kr)%class.eq.fbem_be) then
        do kb=1,region(kr)%n_boundaries
          sb=region(kr)%boundary(kb)
          if ((sb.ge.boundary_eid_min).and.(sb.le.boundary_eid_max)) then
            if (boundary_iid(sb).eq.0) then
              call fbem_warning_message(error_unit,0,'region',region(kr)%id,'wrong definition of its boundaries')
              call fbem_error_message(error_unit,0,'boundary',sb,'does not exist')
            else
              region(kr)%boundary(kb)=boundary_iid(sb)
            end if
          else
            call fbem_warning_message(error_unit,0,'region',region(kr)%id,'wrong definition of its boundaries')
            call fbem_error_message(error_unit,0,'boundary',sb,'does not exist')
          end if
        end do
      end if
    end do

    ! Build BE boundary -> BE regions connectivities.
    ! Loop through boundaries
    do i=1,n_boundaries
      ! Find the number of BE regions of each BE boundary
      boundary(i)%n_regions=0
      do j=1,n_regions
        if (region(j)%class.eq.fbem_be) then
          do k=1,region(j)%n_boundaries
            if (region(j)%boundary(k).eq.i) then
              boundary(i)%n_regions=boundary(i)%n_regions+1
            end if
          end do
        end if
      end do
      ! Stops if boundary(i)%n_regions==0 or boundary(i)%n_regions>2
      if ((boundary(i)%n_regions.eq.0).or.(boundary(i)%n_regions.gt.2)) then
          call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'can only be connected to 1 or 2 BE regions')
      end if
      ! Allocate the connection vector
      allocate (boundary(i)%region(boundary(i)%n_regions))
      allocate (boundary(i)%region_boundary_idx(boundary(i)%n_regions))
      allocate (boundary(i)%region_boundary_reversion(boundary(i)%n_regions))
      ! Find the regions of each boundary
      boundary(i)%n_regions=0
      do j=1,n_regions
        if (region(j)%class.eq.fbem_be) then
          do k=1,region(j)%n_boundaries
            if (region(j)%boundary(k).eq.i) then
              boundary(i)%n_regions=boundary(i)%n_regions+1
              boundary(i)%region(boundary(i)%n_regions)=j
              boundary(i)%region_boundary_idx(boundary(i)%n_regions)=k
              boundary(i)%region_boundary_reversion(boundary(i)%n_regions)=region(j)%boundary_reversion(k)
            end if
          end do
        end if
      end do
      ! If the boundary is connected (limitar a contornos exteriores positivos)
      !if ((boundary(i)%n_regions.eq.1).and.(boundary(i)%region_boundary_reversion(1))) then
      !end if
      ! If the boundary is an interface between 2 BE regions
      if (boundary(i)%n_regions.eq.2) then
        ! Both regions must be different
        if (boundary(i)%region(1).eq.boundary(i)%region(2)) then
          call fbem_warning_message(error_unit,0,'region',region(boundary(i)%region(1))%id,'wrong definition of its boundaries')
          call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'is repeated')
        end if
        ! The boundary has inverted normal in one region, and not inverted in the other region.
        if (boundary(i)%region_boundary_reversion(1).eqv.boundary(i)%region_boundary_reversion(2)) then
          call fbem_warning_message(error_unit,0,'region',region(boundary(i)%region(1))%id,'wrong definition of its boundaries')
          call fbem_warning_message(error_unit,0,'region',region(boundary(i)%region(2))%id,'wrong definition of its boundaries')
          call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'wrong orientation in its regions')
        end if
        !
        ! The boundary(i)%region(1) is the region where the boundary has positive orientation
        !
        ! The region where the boundary has positive orientation goes to 1, and the other to 2. So swap regions if not.
        if (boundary(i)%region_boundary_reversion(1)) then
          ! Swap regions ids
          tmp_int=boundary(i)%region(1)
          boundary(i)%region(1)=boundary(i)%region(2)
          boundary(i)%region(2)=tmp_int
          ! Swap boundary index in the region
          tmp_int=boundary(i)%region_boundary_idx(1)
          boundary(i)%region_boundary_idx(1)=boundary(i)%region_boundary_idx(2)
          boundary(i)%region_boundary_idx(2)=tmp_int
          ! Swap orientations
          tmp_bool=boundary(i)%region_boundary_reversion(1)
          boundary(i)%region_boundary_reversion(1)=boundary(i)%region_boundary_reversion(2)
          boundary(i)%region_boundary_reversion(2)=tmp_bool
        end if
      end if
    end do

    ! ==================
    ! WRITE CONNECTIVITY
    ! ==================

    ! Write the results
    if (verbose_level.ge.4) then
      write(output_unit,'(3x,a,1x,a,1x,a,1x,a,1x,a,1x,a)') '___boundary','__n_regions','__region(1)','rev(1)','__region(2)','rev(2)'
      do i=1,n_boundaries
        write(output_unit,'(3x,i11,1x,i11)',advance='no') boundary(i)%id, boundary(i)%n_regions
        do j=1,boundary(i)%n_regions
          if (j.lt.boundary(i)%n_regions) then
            write(output_unit,'(1x,i11,l6)',advance='no') region(boundary(i)%region(j))%id, boundary(i)%region_boundary_reversion(j)
          else
            write(output_unit,'(1x,i11,l6)') region(boundary(i)%region(j))%id, boundary(i)%region_boundary_reversion(j)
          end if
        end do
      end do
    end if

    ! Ending message
    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  else

    if (n_region_be_boundaries.gt.0) then
      call fbem_error_message(error_unit,0,'['//trim(section_name)//']',0,'this section is required.')
    end if

  end if

end subroutine read_boundaries
