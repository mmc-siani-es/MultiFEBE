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

subroutine build_data

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_string_handling
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! Build element data at geometrical nodes
  call build_data_at_geometrical_nodes
  if (export_geo_mesh_data) call export_data_at_geometrical_nodes

  ! Build element data at functional nodes
  call build_data_at_functional_nodes
  if (export_fun_mesh_data) call export_data_at_functional_nodes

  ! Deberia haber una para cada familia y tipo de elementos (poner aqui EB y T beams orientation)

  ! Build element data for boundary elements
  if (n_be_regions.gt.0) call build_data_of_be_elements

  ! Build element data at collocation points if there is any BE region
  if (n_be_regions.gt.0) call build_data_at_collocation_points

end subroutine build_data
