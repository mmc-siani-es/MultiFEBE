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
!! <b> Subroutine that reads the BE body loads from a file. </b>
subroutine read_be_bodyloads(fileunit)

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
  integer                        :: i, j, k
  integer                        :: n_region_be_bodyloads
  integer                        :: kr, kb, sb

  ! Default value
  n_be_bodyloads=0

  ! Detect if there are BE body loads in BE regions
  n_region_be_bodyloads=0
  do i=1,n_regions
    if (region(i)%class.eq.fbem_be) then
      n_region_be_bodyloads=n_region_be_bodyloads+region(i)%n_be_bodyloads
    end if
  end do

  ! Locate the section
  section_name='be body loads'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! Read the number of BE body loads
    read(fileunit,*) n_be_bodyloads

    ! Allocate and initialize
    if (n_be_bodyloads.gt.0) then
      allocate (be_bodyload(n_be_bodyloads))
      do i=1,n_be_bodyloads
        be_bodyload(i)%id=0
        be_bodyload(i)%name=''
        be_bodyload(i)%part=0
        be_bodyload(i)%region=0
        be_bodyload(i)%coupling=0
        be_bodyload(i)%export=.true.
      end do
    else
      call fbem_error_message(error_unit,0,section_name,0,'the number of BE body loads must be >0')
    end if

    ! Read BE body load id, part
    do i=1,n_be_bodyloads
      ! Read
      read(fileunit,*) be_bodyload(i)%id, be_bodyload(i)%part
      ! Check id
      if (be_bodyload(i)%id.le.0) then
        call fbem_error_message(error_unit,0,'BE body load',be_bodyload(i)%id,'identifiers must be greater than 0')
      end if
    end do

    ! =========================================
    ! CHECK AND BUILD BE BODY LOADS IDENTIFIERS
    ! =========================================

    be_bodyload_eid_min=be_bodyload(1)%id
    be_bodyload_eid_max=be_bodyload(1)%id
    do i=2,n_be_bodyloads
      if (be_bodyload(i)%id.lt.be_bodyload_eid_min) be_bodyload_eid_min=be_bodyload(i)%id
      if (be_bodyload(i)%id.gt.be_bodyload_eid_max) be_bodyload_eid_max=be_bodyload(i)%id
    end do
    allocate (be_bodyload_iid(be_bodyload_eid_min:be_bodyload_eid_max))
    be_bodyload_iid=0
    do i=1,n_be_bodyloads
      if (be_bodyload_iid(be_bodyload(i)%id).ne.0) then
        call fbem_error_message(error_unit,0,'BE body load',be_bodyload(i)%id,'is repeated')
      else
        be_bodyload_iid(be_bodyload(i)%id)=i
      end if
    end do

    ! ================================================================
    ! CHECK AND BUILD CONNECTIVITIES BETWEEN BE BODY LOADS AND REGIONS
    ! ================================================================

    ! Convert eid to iid in BE region -> BE body loads connectivity, and check if the indicated BE body loads actually exist.
    do kr=1,n_regions
      if (region(kr)%class.eq.fbem_be) then
        do kb=1,region(kr)%n_be_bodyloads
          sb=region(kr)%be_bodyload(kb)
          if ((sb.ge.be_bodyload_eid_min).and.(sb.le.be_bodyload_eid_max)) then
            if (be_bodyload_iid(sb).eq.0) then
              call fbem_warning_message(error_unit,0,'region',region(kr)%id,'wrong definition of its BE body loads')
              call fbem_error_message(error_unit,0,'BE body load',sb,'does not exist')
            else
              region(kr)%be_bodyload(kb)=be_bodyload_iid(sb)
            end if
          else
            call fbem_warning_message(error_unit,0,'region',region(kr)%id,'wrong definition of its BE body loads')
            call fbem_error_message(error_unit,0,'BE body load',sb,'does not exist')
          end if
        end do
      end if
    end do

    ! Build BE body loads -> BE region connectivities.
    do i=1,n_be_bodyloads
      do j=1,n_regions
        if (region(j)%class.eq.fbem_be) then
          do k=1,region(j)%n_be_bodyloads
            if (region(j)%be_bodyload(k).eq.i) then
              if (be_bodyload(i)%region.eq.0) then
                be_bodyload(i)%region=j
              else
                call fbem_warning_message(error_unit,0,'region',region(j)%id,'the following BE body load is present')
                call fbem_warning_message(error_unit,0,'region',region(be_bodyload(i)%region)%id,'the following BE body load is present')
                call fbem_error_message(error_unit,0,'BE body load',be_bodyload(i)%id,'is repeated')
              end if
            end if
          end do
        end if
      end do
    end do

    ! Check if each BE body load is connected to a region
    do i=1,n_be_bodyloads
      if (be_bodyload(i)%region.eq.0) then
        call fbem_warning_message(error_unit,0,'['//trim(section_name)//']',be_bodyload(i)%id,'this BE body load is not connected to any region.')
      end if
    end do

    ! Ending message
    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  else

    if (n_region_be_bodyloads.gt.0) then
      call fbem_error_message(error_unit,0,'['//trim(section_name)//']',0,'this section is required.')
    end if

  end if

end subroutine read_be_bodyloads
