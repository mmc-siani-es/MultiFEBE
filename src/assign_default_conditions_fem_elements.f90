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

subroutine assign_default_conditions_fem_elements

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

  ! Local variables
  integer :: ke

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START assigning default conditions (FEM elements)')

  ! ELEMENT LOAD (BODY LOAD FOR SOLID FE, SURFACE LOAD FOR SHELL AND BEAM (THIS SHOULD BE CHANGED))
  do ke=1,n_elements
    if (part(element(ke)%part)%type.eq.fbem_part_fe_subregion) then
      allocate (element(ke)%ctype(2,1))
      element(ke)%ctype(1,1)=0 ! type=1 (local ) ctype=0 (no load), ctype=1 (load).
      element(ke)%ctype(2,1)=0 ! type=2 (global) ctype=0 (no load), ctype=1 (load).
      allocate (element(ke)%cvalue_i(problem%n,element(ke)%n_nodes,2))
      allocate (element(ke)%cvalue_r(problem%n,element(ke)%n_nodes,2))
      allocate (element(ke)%cvalue_c(problem%n,element(ke)%n_nodes,2))
    end if
  end do

  ! SUBEDGE LOADS
  do ke=1,n_subedges
    if (part(element(subedge(ke)%supelement(1))%part)%type.eq.fbem_part_fe_subregion) then
      allocate (subedge(ke)%ctype(2,1))
      subedge(ke)%ctype(1,1)=0 ! type=1 (local ) ctype=0 (no load), ctype=1 (load).
      subedge(ke)%ctype(2,1)=0 ! type=2 (global) ctype=0 (no load), ctype=1 (load).
      allocate (subedge(ke)%cvalue_i(problem%n,subedge(ke)%n_nodes,2))
      allocate (subedge(ke)%cvalue_r(problem%n,subedge(ke)%n_nodes,2))
      allocate (subedge(ke)%cvalue_c(problem%n,subedge(ke)%n_nodes,2))
    end if
  end do

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END assigning default conditions (FEM elements)')

end subroutine assign_default_conditions_fem_elements
