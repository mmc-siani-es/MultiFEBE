! ---------------------------------------------------------------------
! Copyright (C) 2014-2025 Universidad de Las Palmas de Gran Canaria:
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

subroutine assemble_bem_bl_harpot_equation(sb_int,se_int,se_int_n_snodes,sn_col,eq_index,gbl)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_shape_functions
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! No implicit variables
  implicit none

  ! I/O variables
  integer              :: sb_int
  integer              :: se_int
  integer              :: se_int_n_snodes
  integer              :: sn_col
  integer              :: eq_index
  complex(kind=real64) :: gbl(se_int_n_snodes)

  ! Local variables
  integer              :: il, ik
  integer              :: kn_int, sn_int
  integer              :: knt_int, kncbl, se_int_n_nodes
  integer              :: row, col

  ! ==========================
  ! BE BODY LOAD (integration)
  ! ==========================

  select case (be_bodyload(sb_int)%coupling)

    ! ---------
    ! UNCOUPLED
    ! ---------

    case (fbem_bl_uncoupled)
      row=node(sn_col)%row(1,eq_index)
      do kn_int=1,se_int_n_snodes
        sn_int=element(se_int)%node(kn_int)
        b_c(row,1)=b_c(row,1)+gbl(kn_int)*node(sn_int)%value_c(1,1)
      end do

    ! ---------------------------------
    ! BEAM TIP LOAD and SHELL EDGE LOAD
    ! ---------------------------------

    case (fbem_bl_coupling_beam_tip,fbem_bl_coupling_shell_edge)
      stop 'assemble_bem_bl_harpot_equation: not valid 1'

    ! ------------
    ! BEAM / SHELL
    !-------------

    case (fbem_bl_coupling_beam_line,fbem_bl_coupling_shell_surface)
      stop 'assemble_bem_bl_harpot_equation: not valid 2'

  end select

end subroutine assemble_bem_bl_harpot_equation
