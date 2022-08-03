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

subroutine assemble_bem_staela_equation(kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,eq_index,hp,gp,hm,gm)

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
  integer           :: kr
  integer           :: sb_int
  logical           :: sb_int_reversion
  integer           :: se_int
  integer           :: se_int_n_nodes
  integer           :: sn_col
  integer           :: eq_index
  real(kind=real64) :: hp(se_int_n_nodes,problem%n,problem%n)
  real(kind=real64) :: gp(se_int_n_nodes,problem%n,problem%n)
  real(kind=real64) :: hm(se_int_n_nodes,problem%n,problem%n)
  real(kind=real64) :: gm(se_int_n_nodes,problem%n,problem%n)

  ! Local variables
  integer :: il, ik
  integer :: kn_int, sn_int
  integer :: row, col
  integer :: sn_int_fe

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
          do il=1,problem%n
            row=node(sn_col)%row(il,eq_index)
            do ik=1,problem%n
              do kn_int=1,se_int_n_nodes
                sn_int=element(se_int)%node(kn_int)
                select case (node(sn_int)%ctype(ik,1))
                  ! u_k known, p_k unknown
                  case (0)
                    col=node(sn_int)%col(problem%n+ik,1)
                    A_r(row,col)=A_r(row,col)-gp(kn_int,il,ik)
                    b_r(row,1)=b_r(row,1)-hp(kn_int,il,ik)*node(sn_int)%cvalue_r(ik,1,1)
                  ! p_k known, u_k unknown
                  case (1)
                    col=node(sn_int)%col(ik,1)
                    A_r(row,col)=A_r(row,col)+hp(kn_int,il,ik)
                    b_r(row,1)=b_r(row,1)+gp(kn_int,il,ik)*node(sn_int)%cvalue_r(ik,1,1)
                  ! p_k known (normal pressure), u_k unknown
                  case (10)
                    col=node(sn_int)%col(ik,1)
                    A_r(row,col)=A_r(row,col)+hp(kn_int,il,ik)
                    if (sb_int_reversion.eqv.(.false.)) then
                      b_r(row,1)=b_r(row,1)+gp(kn_int,il,ik)*node(sn_int)%cvalue_r(ik,1,1)*node(sn_int)%n_fn(ik)
                    else
                      b_r(row,1)=b_r(row,1)-gp(kn_int,il,ik)*node(sn_int)%cvalue_r(ik,1,1)*node(sn_int)%n_fn(ik)
                    end if
                end select
              end do
            end do
          end do

        ! -------------------
        ! CRACK-LIKE BOUNDARY
        ! -------------------

        case (fbem_boundary_class_cracklike)
          do il=1,problem%n
            row=node(sn_col)%row(il,eq_index)
            do ik=1,problem%n
              do kn_int=1,se_int_n_nodes
                sn_int=element(se_int)%node(kn_int)
                ! Face +
                select case (node(sn_int)%ctype(ik,1))
                  ! u_k known, p_k unknown
                  case (0)
                    col=node(sn_int)%col(problem%n+ik,1)
                    A_r(row,col)=A_r(row,col)-gp(kn_int,il,ik)
                    b_r(row,1)=b_r(row,1)-hp(kn_int,il,ik)*node(sn_int)%cvalue_r(ik,1,1)
                  ! p_k known, u_k unknown
                  case (1)
                    col=node(sn_int)%col(ik,1)
                    A_r(row,col)=A_r(row,col)+hp(kn_int,il,ik)
                    b_r(row,1)=b_r(row,1)+gp(kn_int,il,ik)*node(sn_int)%cvalue_r(ik,1,1)
                end select
                ! Face -
                select case (node(sn_int)%ctype(ik,2))
                  ! u_k known, p_k unknown
                  case (0)
                    col=node(sn_int)%col(problem%n+ik,2)
                    A_r(row,col)=A_r(row,col)-gm(kn_int,il,ik)
                    b_r(row,1)=b_r(row,1)-hm(kn_int,il,ik)*node(sn_int)%cvalue_r(ik,1,2)
                  ! p_k known, u_k unknown
                  case (1)
                    col=node(sn_int)%col(ik,2)
                    A_r(row,col)=A_r(row,col)+hm(kn_int,il,ik)
                    b_r(row,1)=b_r(row,1)+gm(kn_int,il,ik)*node(sn_int)%cvalue_r(ik,1,2)
                end select
              end do
            end do
          end do
      end select

    ! ==============
    ! BE-BE BOUNDARY
    ! ==============

    case (fbem_boundary_coupling_be_be)
      ! If the integration boundary has N+, then it is the region 1
      if (sb_int_reversion.eqv.(.false.)) then
        do il=1,problem%n
          row=node(sn_col)%row(il,eq_index)
          do ik=1,problem%n
            do kn_int=1,se_int_n_nodes
              sn_int=element(se_int)%node(kn_int)
              ! u_k for region 1
              col=node(sn_int)%col(          ik,1)
              A_r(row,col)=A_r(row,col)+hp(kn_int,il,ik)
              ! p_k for region 1
              col=node(sn_int)%col(problem%n+ik,1)
              A_r(row,col)=A_r(row,col)-gp(kn_int,il,ik)
            end do
          end do
        end do
      ! If the integration boundary has N-, then it is the region 2
      else
        do il=1,problem%n
          row=node(sn_col)%row(il,eq_index)
          do ik=1,problem%n
            do kn_int=1,se_int_n_nodes
              sn_int=element(se_int)%node(kn_int)
              ! u_k for region 2 (u_k_{region 1}=u_k_{region 2})
              col=node(sn_int)%col(          ik,1)
              A_r(row,col)=A_r(row,col)+hp(kn_int,il,ik)
              ! p_k for region 2 (p_k_{region 1}=-p_k_{region 2})
              col=node(sn_int)%col(problem%n+ik,1)
              A_r(row,col)=A_r(row,col)+gp(kn_int,il,ik)
            end do
          end do
        end do
      end if

    ! ==============
    ! BE-FE BOUNDARY
    ! ==============

    case (fbem_boundary_coupling_be_fe)
      select case (boundary(sb_int)%class)

        ! -----------------
        ! ORDINARY BOUNDARY
        ! -----------------

        case (fbem_boundary_class_ordinary)
         do il=1,problem%n
            row=node(sn_col)%row(il,eq_index)
            do ik=1,problem%n
              do kn_int=1,se_int_n_nodes
                sn_int=element(se_int)%node(kn_int)
                sn_int_fe=element(se_int)%element_node(kn_int)
                ! t_k is unknown
                col=node(sn_int)%col(problem%n+ik,1)
                A_r(row,col)=A_r(row,col)-gp(kn_int,il,ik)
                ! If u_k is unknown
                if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                  col=node(sn_int_fe)%col(ik,1)
                  A_r(row,col)=A_r(row,col)+hp(kn_int,il,ik)
                ! If u_k is known
                else
                  b_r(row,1)=b_r(row,1)-hp(kn_int,il,ik)*node(sn_int_fe)%cvalue_r(ik,1,1)
                end if
              end do
            end do
          end do

        ! -------------------
        ! CRACK-LIKE BOUNDARY
        ! -------------------


        case (fbem_boundary_class_cracklike)
!          !
!          ! SBIE (unknown Σt_k)
!          !
!          do il=1,problem%n
!            row=node(sn_col)%row(il,eq_index)
!            do ik=1,problem%n
!              do kn_int=1,se_int_n_nodes
!                sn_int=element(se_int)%node(kn_int)
!                sn_int_fe=element(se_int)%element_node(kn_int)
!                ! Σt_k = t_k^+ + t_k^- is unknown
!                col=node(sn_int)%col(problem%n+ik,1)
!                A_r(row,col)=A_r(row,col)-gp(kn_int,il,ik)
!                ! Note: h^+ + h^- = 0 always except when it contains the free-term, so this gives the required contribution
!                ! If u_k=u_k_{+}=u_k_{-} is unknown
!                if (node(sn_int_fe)%ctype(ik,1).eq.1) then
!                  col=node(sn_int_fe)%col(ik,1)
!                  A_r(row,col)=A_r(row,col)+hp(kn_int,il,ik)+hm(kn_int,il,ik)
!                ! If u_k=u_k_{+}=u_k_{-} is known
!                else
!                  b_r(row,1)=b_r(row,1)-(hp(kn_int,il,ik)+hm(kn_int,il,ik))*node(sn_int_fe)%cvalue_r(ik,1,1)
!                end if
!              end do
!            end do
!          end do
          !
          ! Caso más general con DBEM
          !
          do il=1,problem%n
            row=node(sn_col)%row(il,eq_index)
            do ik=1,problem%n
              do kn_int=1,se_int_n_nodes
                sn_int=element(se_int)%node(kn_int)
                sn_int_fe=element(se_int)%element_node(kn_int)
                ! t_k for face + is unknown
                col=node(sn_int)%col(problem%n+ik,1)
                A_r(row,col)=A_r(row,col)-gp(kn_int,il,ik)
                ! t_k for face - is unknown
                col=node(sn_int)%col(problem%n+ik,2)
                A_r(row,col)=A_r(row,col)-gm(kn_int,il,ik)
                ! If u_k=u_k_{+}=u_k_{-} is unknown
                if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                  col=node(sn_int_fe)%col(ik,1)
                  A_r(row,col)=A_r(row,col)+hp(kn_int,il,ik)+hm(kn_int,il,ik)
                ! If u_k=u_k_{+}=u_k_{-} is known
                else
                  b_r(row,1)=b_r(row,1)-(hp(kn_int,il,ik)+hm(kn_int,il,ik))*node(sn_int_fe)%cvalue_r(ik,1,1)
                end if
              end do
            end do
          end do

      end select

    ! =================
    ! BE-FE-BE BOUNDARY
    ! =================

    case (fbem_boundary_coupling_be_fe_be)
      ! If the integration boundary has N+, then it is the region 1
      if (sb_int_reversion.eqv.(.false.)) then

        do il=1,problem%n
          row=node(sn_col)%row(il,eq_index)
          do ik=1,problem%n
            do kn_int=1,se_int_n_nodes
              sn_int=element(se_int)%node(kn_int)
              sn_int_fe=element(se_int)%element_node(kn_int)
              ! t_k for region 1
              col=node(sn_int)%col(problem%n+ik,1)
              A_r(row,col)=A_r(row,col)-gp(kn_int,il,ik)
              ! If u_k is unknown
              if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                col=node(sn_int_fe)%col(ik,1)
                A_r(row,col)=A_r(row,col)+hp(kn_int,il,ik)
              ! If u_k is known
              else
                b_r(row,1)=b_r(row,1)-hp(kn_int,il,ik)*node(sn_int_fe)%cvalue_r(ik,1,1)
              end if
            end do
          end do
        end do

      ! If the integration boundary has N-, then it is the region 2
      else

        do il=1,problem%n
          row=node(sn_col)%row(il,eq_index)
          do ik=1,problem%n
            do kn_int=1,se_int_n_nodes
              sn_int=element(se_int)%node(kn_int)
              sn_int_fe=element(se_int)%element_node(kn_int)
              ! t_k for region 2
              col=node(sn_int)%col(problem%n+ik,2)
              A_r(row,col)=A_r(row,col)-gp(kn_int,il,ik)
              ! If u_k is unknown
              if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                col=node(sn_int_fe)%col(ik,1)
                A_r(row,col)=A_r(row,col)+hp(kn_int,il,ik)
              ! If u_k is known
              else
                b_r(row,1)=b_r(row,1)-hp(kn_int,il,ik)*node(sn_int_fe)%cvalue_r(ik,1,1)
              end if
            end do
          end do
        end do

      end if

  end select

end subroutine assemble_bem_staela_equation
