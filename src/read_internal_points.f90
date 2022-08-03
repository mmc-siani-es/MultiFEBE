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


subroutine read_internal_points(input_fileunit)
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
  integer                        :: input_fileunit    ! Input file unit
  ! Local
  character(len=fbem_stdcharlen) :: section_name      ! Name of the section
  logical                        :: found             ! Logical variable for sections and keywords
  integer                        :: i                 ! Counters
  integer                        :: kr, kc

  ! Locate the section
  section_name='internal points'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(input_fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! Read the number of internal points
    read(input_fileunit,*) n_internalpoints
    ! Allocate the internal points data structure
    allocate (internalpoint(n_internalpoints))
    ! Loop through the internal points
    do i=1,n_internalpoints
      ! Allocate
      allocate (internalpoint(i)%x(problem%n))

      ! Read the node definition
      read(input_fileunit,*) internalpoint(i)%id, internalpoint(i)%region, (internalpoint(i)%x(kc),kc=1,problem%n)

      ! ====================================
      ! CHECK IF THE INDICATED REGION EXISTS
      ! ====================================

      if ((internalpoint(i)%region.ge.region_eid_min).and.(internalpoint(i)%region.le.region_eid_max)) then
        if (region_iid(internalpoint(i)%region).eq.0) then
          call fbem_error_message(error_unit,0,'internal point',internalpoint(i)%id,'the indicated region does not exist')
        else
          internalpoint(i)%region=region_iid(internalpoint(i)%region)
        end if
      else
        call fbem_error_message(error_unit,0,'internal point',internalpoint(i)%id,'the indicated region does not exist')
      end if

    end do

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
    allocate (internalpoint_iid(internalpoint_eid_min:internalpoint_eid_max))
    internalpoint_iid=0
    do i=1,n_internalpoints
      if (internalpoint_iid(internalpoint(i)%id).ne.0) then
        call fbem_error_message(error_unit,0,'internal point',internalpoint(i)%id,'is repeated')
      else
        internalpoint_iid(internalpoint(i)%id)=i
      end if
    end do

    ! Export them by default
    do i=1,n_internalpoints
      internalpoint(i)%export=.true.
    end do

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  else

    n_internalpoints=0
    do kr=1,n_regions
      region(kr)%n_internalpoints=0
    end do

  end if

end subroutine read_internal_points
