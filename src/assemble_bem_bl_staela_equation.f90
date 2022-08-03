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

subroutine assemble_bem_bl_staela_equation(sb_int,se_int,se_int_n_snodes,sn_col,eq_index,gbl)

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
  integer           :: sb_int
  integer           :: se_int
  integer           :: se_int_n_snodes
  integer           :: sn_col
  integer           :: eq_index
  real(kind=real64) :: gbl(se_int_n_snodes,problem%n,problem%n)

  ! Local variables
  integer           :: il, ik
  integer           :: kn_int, sn_int
  integer           :: knt_int, kncbl, se_int_n_nodes
  integer           :: row, col

  ! ===================================
  ! BE BODY LOAD COUPLING (integration)
  ! ===================================

  select case (be_bodyload(sb_int)%coupling)

    ! ---------
    ! UNCOUPLED
    ! ---------

    case (0)
      do il=1,problem%n
        row=node(sn_col)%row(il,eq_index)
        do ik=1,problem%n
          do kn_int=1,se_int_n_snodes
            sn_int=element(se_int)%node(kn_int)
            b_r(row,1)=b_r(row,1)+gbl(kn_int,il,ik)*node(sn_int)%cvalue_r(problem%n+ik,1,1)
          end do
        end do
      end do

    ! ---------------------------------
    ! BEAM TIP LOAD and SHELL EDGE LOAD
    ! ---------------------------------

    case (fbem_bl_coupling_beam_tip,fbem_bl_coupling_shell_edge)
      se_int_n_nodes=element(se_int)%n_nodes
      do il=1,problem%n
        row=node(sn_col)%row(il,eq_index)
        do ik=1,problem%n
          do kn_int=1,se_int_n_nodes
            sn_int=element(se_int)%node(kn_int)
            do knt_int=1,se_int_n_snodes/se_int_n_nodes
              kncbl=(knt_int-1)*se_int_n_nodes+kn_int
              col=node(sn_int)%col(problem%n+ik,knt_int)
              A_r(row,col)=A_r(row,col)-gbl(kncbl,il,ik)
            end do
          end do
        end do
      end do

    ! ------------
    ! BEAM / SHELL
    ! ------------

    case (fbem_bl_coupling_beam_line,fbem_bl_coupling_shell_surface)
      do il=1,problem%n
        row=node(sn_col)%row(il,eq_index)
        do ik=1,problem%n
          do kn_int=1,se_int_n_snodes
            sn_int=element(se_int)%node(kn_int)
            col=node(sn_int)%col(problem%n+ik,1)
            A_r(row,col)=A_r(row,col)-gbl(kn_int,il,ik)
          end do
        end do
      end do

  end select

end subroutine assemble_bem_bl_staela_equation
