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

subroutine read_internal(input_fileunit)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_geometry
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_symmetry
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! No implicit variables are allowed
  implicit none

  integer                               :: input_fileunit    ! Input file unit
  integer                               :: kr, i

  n_internalpoints=0
  internalelements=.false.
  if (n_be_regions.gt.0) then
    call read_internal_points(input_fileunit)
    call read_internal_points_from_mesh(input_fileunit)
  end if
  call read_internal_elements(input_fileunit)

  if (n_be_regions.gt.0) then

    ! ==========================================
    ! BUILD REGION->INTERNAL POINTS CONNECTIVITY
    ! ==========================================

    do kr=1,n_regions
      if (region(kr)%class.eq.fbem_be) then
        ! Calculate the number of internal point in the region kr
        region(kr)%n_internalpoints=0
        do i=1,n_internalpoints
          if (internalpoint(i)%region.eq.kr) region(kr)%n_internalpoints=region(kr)%n_internalpoints+1
        end do
        ! Allocate
        allocate (region(kr)%internalpoint(region(kr)%n_internalpoints))
        ! Build connectivity
        if (region(kr)%n_internalpoints.gt.0) then
          region(kr)%n_internalpoints=0
          do i=1,n_internalpoints
            if (internalpoint(i)%region.eq.kr) then
              region(kr)%n_internalpoints=region(kr)%n_internalpoints+1
              region(kr)%internalpoint(region(kr)%n_internalpoints)=i
            end if
          end do
        end if
      end if
    end do

  end if

end subroutine read_internal
