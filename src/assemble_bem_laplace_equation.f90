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

subroutine assemble_bem_laplace_equation(sb_int,sb_int_reversion,se_int,se_int_n_nodes,hp,gp,hm,gm,sn_col,eq_index)

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
  logical           :: sb_int_reversion
  integer           :: se_int
  integer           :: se_int_n_nodes
  real(kind=real64) :: hp(se_int_n_nodes)
  real(kind=real64) :: gp(se_int_n_nodes)
  real(kind=real64) :: hm(se_int_n_nodes)
  real(kind=real64) :: gm(se_int_n_nodes)
  integer           :: sn_col
  integer           :: eq_index

  ! Local variables
  integer           :: kn_int, sn_int
  integer           :: row, col

  select case (boundary(sb_int)%coupling)

    ! ===========
    ! BE BOUNDARY
    ! ===========

    case (fbem_boundary_coupling_be)
      select case (boundary(sb_int)%class)

        ! -----------------
        ! ORDINARY BOUNDARY
        ! -----------------

        case (fbem_boundary_class_ordinary)
          row=node(sn_col)%row(1,eq_index)
          do kn_int=1,se_int_n_nodes
            sn_int=element(se_int)%node(kn_int)
            select case (node(sn_int)%ctype(1,1))
              ! p known, j unknown
              case (0)
                col=node(sn_int)%col(2,1)
                A_r(row,col)=A_r(row,col)-gp(kn_int)
                b_r(row,1)=b_r(row,1)-hp(kn_int)*node(sn_int)%cvalue_r(1,1,1)
              ! j known, p unknown
              case (1)
                col=node(sn_int)%col(1,1)
                A_r(row,col)=A_r(row,col)+hp(kn_int)
                b_r(row,1)=b_r(row,1)+gp(kn_int)*node(sn_int)%cvalue_r(1,1,1)
            end select
          end do

        ! -------------------
        ! CRACK-LIKE BOUNDARY
        ! -------------------

        case (fbem_boundary_class_cracklike)
          row=node(sn_col)%row(1,eq_index)
          do kn_int=1,se_int_n_nodes
            sn_int=element(se_int)%node(kn_int)
            ! Face +
            select case (node(sn_int)%ctype(1,1))
              ! p known, j unknown
              case (0)
                col=node(sn_int)%col(2,1)
                A_r(row,col)=A_r(row,col)-gp(kn_int)
                b_r(row,1)=b_r(row,1)-hp(kn_int)*node(sn_int)%cvalue_r(1,1,1)
              ! j known, p unknown
              case (1)
                col=node(sn_int)%col(1,1)
                A_r(row,col)=A_r(row,col)+hp(kn_int)
                b_r(row,1)=b_r(row,1)+gp(kn_int)*node(sn_int)%cvalue_r(1,1,1)
            end select
            ! Face -
            select case (node(sn_int)%ctype(1,2))
              ! p known, j unknown
              case (0)
                col=node(sn_int)%col(2,2)
                A_r(row,col)=A_r(row,col)-gm(kn_int)
                b_r(row,1)=b_r(row,1)-hm(kn_int)*node(sn_int)%cvalue_r(1,1,2)
              ! j known, p unknown
              case (1)
                col=node(sn_int)%col(1,2)
                A_r(row,col)=A_r(row,col)+hm(kn_int)
                b_r(row,1)=b_r(row,1)+gm(kn_int)*node(sn_int)%cvalue_r(1,1,2)
            end select
          end do
      end select

    ! ==============
    ! BE-BE BOUNDARY
    ! ==============

    case (fbem_boundary_coupling_be_be)
      if (sb_int_reversion.eqv.(.false.)) then
        row=node(sn_col)%row(1,eq_index)
        do kn_int=1,se_int_n_nodes
          sn_int=element(se_int)%node(kn_int)
          ! p for region 1
          col=node(sn_int)%col(1,1)
          A_r(row,col)=A_r(row,col)+hp(kn_int)
          ! j for region 1
          col=node(sn_int)%col(2,1)
          A_r(row,col)=A_r(row,col)-gp(kn_int)
        end do
      else
        row=node(sn_col)%row(1,eq_index)
        do kn_int=1,se_int_n_nodes
          sn_int=element(se_int)%node(kn_int)
          ! p for region 2 (p_{region 1}=p_{region 2})
          col=node(sn_int)%col(1,1)
          A_r(row,col)=A_r(row,col)+hp(kn_int)
          ! j for region 2 (j_{region 1}=-j_{region 2})
          col=node(sn_int)%col(2,1)
          A_r(row,col)=A_r(row,col)+gp(kn_int)
        end do
      end if

    ! ==============
    ! BE-FE BOUNDARY
    ! ==============

    case (fbem_boundary_coupling_be_fe)
      stop 'not implemented yet'

    ! =================
    ! BE-FE-BE BOUNDARY
    ! =================

    case (fbem_boundary_coupling_be_fe_be)
      stop 'not implemented yet'

  end select

end subroutine assemble_bem_laplace_equation
