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

subroutine assemble_bem_harela_equation(omega,kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,eq_index,&
                                        hp,gp,hm,gm,u_inc,t_inc)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_shape_functions
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! No implicit variables
  implicit none

  ! I/O variables
  real(kind=real64)    :: omega
  integer              :: kr
  integer              :: sb_int
  logical              :: sb_int_reversion
  integer              :: se_int
  integer              :: se_int_n_nodes
  integer              :: sn_col
  integer              :: eq_index
  complex(kind=real64) :: hp (se_int_n_nodes,problem%n,problem%n)
  complex(kind=real64) :: gp (se_int_n_nodes,problem%n,problem%n)
  complex(kind=real64) :: hm(se_int_n_nodes,problem%n,problem%n)
  complex(kind=real64) :: gm(se_int_n_nodes,problem%n,problem%n)
  complex(kind=real64) :: u_inc(problem%n,se_int_n_nodes)
  complex(kind=real64) :: t_inc(problem%n,se_int_n_nodes)

  ! Local variables
  integer              :: il, ik
  integer              :: kn_int, sn_int, sn_int_fe
  integer              :: row, col
  complex(kind=real64) :: KT

  ! ========================================================================
  ! ASSEMBLE TOTAL (IF AN INCIDENT WAVE FIELD IS PRESENT) OR SCATTERED FIELD
  ! ========================================================================

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
                  ! u_k known, t_k unknown
                  case (0)
                    col=node(sn_int)%col(problem%n+ik,1)
                    A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
                    b_c(row,1)=b_c(row,1)-hp(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,1)
                  ! t_k known, u_k unknown
                  case (1)
                    col=node(sn_int)%col(          ik,1)
                    A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                    b_c(row,1)=b_c(row,1)+gp(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,1)
                  ! p_k known (normal pressure), u_k unknown
                  case (10)
                    col=node(sn_int)%col(ik,1)
                    A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                    if (sb_int_reversion.eqv.(.false.)) then
                      b_c(row,1)=b_c(row,1)+gp(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,1)*node(sn_int)%n_fn(ik)
                    else
                      b_c(row,1)=b_c(row,1)-gp(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,1)*node(sn_int)%n_fn(ik)
                    end if
                  ! u_k unknown, t_k unknown
                  case (2,3)
                    col=node(sn_int)%col(          ik,1)
                    A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                    col=node(sn_int)%col(problem%n+ik,1)
                    A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
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
                  ! u_k known, t_k unknown
                  case (0)
                    col=node(sn_int)%col(problem%n+ik,1)
                    A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
                    b_c(row,1)=b_c(row,1)-hp(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,1)
                  ! t_k known, u_k unknown
                  case (1)
                    col=node(sn_int)%col(ik,1)
                    A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                    b_c(row,1)=b_c(row,1)+gp(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,1)
                end select
                ! Face -
                select case (node(sn_int)%ctype(ik,2))
                  ! u_k known, t_k unknown
                  case (0)
                    col=node(sn_int)%col(problem%n+ik,2)
                    A_c(row,col)=A_c(row,col)-gm(kn_int,il,ik)
                    b_c(row,1)=b_c(row,1)-hm(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,2)
                  ! t_k known, u_k unknown
                  case (1)
                    col=node(sn_int)%col(ik,2)
                    A_c(row,col)=A_c(row,col)+hm(kn_int,il,ik)
                    b_c(row,1)=b_c(row,1)+gm(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,2)
                end select
              end do
            end do
          end do

      end select

    ! ==============
    ! BE-BE BOUNDARY
    ! ==============

    case (fbem_boundary_coupling_be_be)
      if (sb_int_reversion.eqv.(.false.)) then
        select case (region(boundary(sb_int)%region(2))%type)

          ! -------------------------------------------
          ! VISCOELASTIC SOLID (1) - INVISCID FLUID (2)
          ! -------------------------------------------

          case (fbem_potential)
            do il=1,problem%n
              row=node(sn_col)%row(il,eq_index)
              do ik=1,problem%n
                do kn_int=1,se_int_n_nodes
                  sn_int=element(se_int)%node(kn_int)
                  ! u_k for region 1
                  col=node(sn_int)%col(ik,1)
                  A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                  ! t_k for region 1 is substituted by t_k_{region 1}=-p_{region 2}·n_k_{region 1}
                  ! note that element(se_int)%n_fn(:,:) is the unit normal for region 1
                  col=node(sn_int)%col(1,2)
                  A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                end do
              end do
            end do

          ! -----------------------------------------------
          ! VISCOELASTIC SOLID (1) - VISCOELASTIC SOLID (2)
          ! -----------------------------------------------

          case (fbem_viscoelastic)

            do il=1,problem%n
              row=node(sn_col)%row(il,eq_index)
              do ik=1,problem%n
                do kn_int=1,se_int_n_nodes
                  sn_int=element(se_int)%node(kn_int)
                  select case (node(sn_int)%ctype(ik,1))
                    ! Global axes
                    ! 0 - perfect bonding - u_k^{1} and t_k^{1} are active since:
                    ! t_k^{(2)}=-t_k^{(1)}
                    ! u_k^{(2)}= u_k^{(1)}
                    case (0)
                      ! u_k for region 1
                      col=node(sn_int)%col(          ik,1)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                      ! t_k for region 1
                      col=node(sn_int)%col(problem%n+ik,1)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
                    ! Global axes
                    ! 1 - perfect debonding - u_k^{1} and u_k^{2} are active since:
                    ! t_k^{(1)}=0
                    ! t_k^{(2)}=0
                    case (1)
                      ! u_k for region 1
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                    ! Global axes
                    ! 2 - partial bonding - u_k^{1} and u_k^{2} are active since:
                    ! t_k^{(2)}+t_k^{(1)}=0
                    ! (u_k^{2}-u_k^{1})*KT=t_k^{1}
                    ! Then
                    ! t_k^{(1)}=(u_k^{2}-u_k^{1})*KT
                    ! t_k^{(2)}=-(u_k^{2}-u_k^{1})*KT
                    case (2)
                      ! KT = K + i*omega*C - omega**2*M
                      KT=node(sn_int)%cvalue_c(ik,1,1)+c_im*omega*node(sn_int)%cvalue_c(ik,2,1)&
                        -omega**2*node(sn_int)%cvalue_c(ik,3,1)
                      ! u_k for region 1
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                      ! t_k for region 1: t_k^{1}=u_k^{2}*KT-u_k^{1}*KT
                      ! u_k^{1}
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*KT
                      ! u_k^{2}
                      col=node(sn_int)%col(ik,2)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)*KT
                    ! Local axes
                    ! 3 - perfect bonding - u_k^{1}, u_k^{2} and t_k^{1} are active since:
                    ! t_k^{(2)}=-t_k^{(1)}
                    ! (u^{2}-u^{1})·l_k^(1)=0
                    ! 4 - perfect debonding - u_k^{1}, u_k^{2} and t_k^{1} are active since:
                    ! t_k^{(2)}=-t_k^{(1)}
                    ! t^(1)·l_k^(1)=0
                    ! 5 - partial bonding - u_k^{1}, u_k^{2} and t_k^{1} are active since:
                    ! t_k^{(2)}=-t_k^{(1)}
                    ! (u^{2}-u^{1})·l_k^(1)*KT=t^(1)·l_k^(1)
                    case (3,4,5)
                      ! u_k for region 1
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                      ! t_k for region 1
                      col=node(sn_int)%col(problem%n+ik,1)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
                  end select
                end do
              end do
            end do

          ! ----------------------------------------------
          ! VISCOELASTIC SOLID (1) - POROELASTIC MEDIA (2)
          ! ----------------------------------------------

          case (fbem_poroelastic)
            do il=1,problem%n
              row=node(sn_col)%row(il,eq_index)
              do ik=1,problem%n
                do kn_int=1,se_int_n_nodes
                  sn_int=element(se_int)%node(kn_int)
                  select case (node(sn_int)%ctype(ik,1))
                    ! Global axes
                    ! 0 - perfect bonding - u_k^{2} and t_k^{2} are active since:
                    ! t_k^{(1)}=-t_k^{(2)}+tau^(2)n_k^(1)
                    ! u_k^{(1)}= u_k^{(2)}
                    case (0)
                      ! u_k for region 1: u_k^(1)= u_k^(2)
                      col=node(sn_int)%col(ik,2)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                      ! t_k for region 1: t_k^(1)=-t_k^(2)+tau^(2)n_k^(1)
                      ! t_k^(2)
                      col=node(sn_int)%col(problem%n+1+ik,2)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)
                      ! tau^(2)
                      col=node(sn_int)%col(0,2)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                    ! Global axes
                    ! 1 - perfect debonding - u_k^{1} and u_k^{2} are active since:
                    ! t_k^{(1)}=0
                    ! t_k^{(2)}=tau^(2)n_k^(1)
                    case (1)
                      ! u_k for region 1
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                    ! Global axes
                    ! 2 - partial bonding - u_k^{1} and u_k^{2} are active since:
                    ! t_k^{(1)}=(u_k^{2}-u_k^{1})*KT
                    ! t_k^{(2)}=-(u_k^{2}-u_k^{1})*KT+tau^(2)n_k^(1)
                    case (2)
                      ! KT = K + i*omega*C - omega**2*M
                      KT=node(sn_int)%cvalue_c(ik,1,1)+c_im*omega*node(sn_int)%cvalue_c(ik,2,1)&
                        -omega**2*node(sn_int)%cvalue_c(ik,3,1)
                      ! u_k for region 1
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                      ! t_k for region 1: t_k^{1}=u_k^{2}*KT-u_k^{1}*KT
                      ! u_k^{1}
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*KT
                      ! u_k^{2}
                      col=node(sn_int)%col(ik,2)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)*KT
                    ! Local axes
                    ! In all cases, u_k^{1}, u_k^{2} and t_k^{1} are active, and t_k^{(2)}=-t_k^{(1)}+tau^(2)n_k^(1).
                    ! Each B.C. needs an additional equation:
                    ! 3 - perfect bonding  : (u^{2}-u^{1})·l_k^(1)=0
                    ! 4 - perfect debonding: t^(1)·l_k^(1)=0
                    ! 5 - partial bonding  : (u^{2}-u^{1})·l_k^(1)*KT=t^(1)·l_k^(1)
                    case (3,4,5)
                      ! u_k for region 1
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                      ! t_k for region 1
                      col=node(sn_int)%col(problem%n+ik,1)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
                  end select
                end do
              end do
            end do

        end select
      else
        select case (region(boundary(sb_int)%region(1))%type)

          ! -------------------------------------------
          ! INVISCID FLUID (1) - VISCOELASTIC SOLID (2)
          ! -------------------------------------------

          case (fbem_potential)
            do il=1,problem%n
              row=node(sn_col)%row(il,eq_index)
              do ik=1,problem%n
                do kn_int=1,se_int_n_nodes
                  sn_int=element(se_int)%node(kn_int)
                  ! u_k for region 2
                  col=node(sn_int)%col(ik,2)
                  A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                  ! t_k for region 2 is substituted by t_k_{region 2}=-p_{region 1}·n_k_{region 2}
                  ! note that element(se_int)%n_fn(:,:) is the unit normal for region 1
                  col=node(sn_int)%col(1,1)
                  A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                end do
              end do
            end do

          ! -----------------------------------------------
          ! VISCOELASTIC SOLID (1) - VISCOELASTIC SOLID (2)
          ! -----------------------------------------------

          case (fbem_viscoelastic)

            do il=1,problem%n
              row=node(sn_col)%row(il,eq_index)
              do ik=1,problem%n
                do kn_int=1,se_int_n_nodes
                  sn_int=element(se_int)%node(kn_int)
                  select case (node(sn_int)%ctype(ik,1))
                    ! Global axes
                    ! 0 - perfect bonding - u_k^{1} and t_k^{1} are active since:
                    ! t_k^{(2)}=-t_k^{(1)}
                    ! u_k^{(2)}= u_k^{(1)}
                    case (0)
                      ! u_k for region 2: u_k^{(2)}= u_k^{(1)}
                      col=node(sn_int)%col(          ik,1)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                      ! t_k for region 2: t_k^{(2)}=-t_k^{(1)}
                      col=node(sn_int)%col(problem%n+ik,1)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)
                    ! Global axes
                    ! 1 - perfect debonding - u_k^{1} and u_k^{2} are active since:
                    ! t_k^{(1)}=0
                    ! t_k^{(2)}=0
                    case (1)
                      ! u_k for region 2
                      col=node(sn_int)%col(ik,2)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                    ! Global axes
                    ! 2 - partial bonding - u_k^{1} and u_k^{2} are active since:
                    ! t_k^{(2)}+t_k^{(1)}=0
                    ! (u_k^{2}-u_k^{1})*KT=t_k^{1}
                    ! Then
                    ! t_k^{(1)}=(u_k^{2}-u_k^{1})*KT
                    ! t_k^{(2)}=-(u_k^{2}-u_k^{1})*KT
                    case (2)
                      ! KT = K + i*omega*C - omega**2*M
                      KT=node(sn_int)%cvalue_c(ik,1,1)+c_im*omega*node(sn_int)%cvalue_c(ik,2,1)&
                        -omega**2*node(sn_int)%cvalue_c(ik,3,1)
                      ! u_k for region 2
                      col=node(sn_int)%col(ik,2)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                      ! t_k for region 2: t_k^{(2)}=-u_k^{2}*KT+u_k^{1}*KT
                      ! u_k^{1}
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)*KT
                      ! u_k^{2}
                      col=node(sn_int)%col(ik,2)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*KT
                    ! Local axes
                    ! 3 - perfect bonding - u_k^{1}, u_k^{2} and t_k^{1} are active since:
                    ! t_k^{(2)}=-t_k^{(1)}
                    ! (u^{2}-u^{1})·l_k^(1)=0
                    ! 4 - perfect debonding - u_k^{1}, u_k^{2} and t_k^{1} are active since:
                    ! t_k^{(2)}=-t_k^{(1)}
                    ! t^(1)·l_k^(1)=0
                    ! 5 - partial bonding - u_k^{1}, u_k^{2} and t_k^{1} are active since:
                    ! t_k^{(2)}=-t_k^{(1)}
                    ! (u^{2}-u^{1})·l_k^(1)*KT=t^(1)·l_k^(1)
                    case (3,4,5)
                      ! u_k for region 2
                      col=node(sn_int)%col(          ik,2)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                      ! t_k for region 2: t_k^{(2)}=-t_k^{(1)}
                      col=node(sn_int)%col(problem%n+ik,1)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)
                  end select
                end do
              end do
            end do

          ! ----------------------------------------------
          ! POROELASTIC MEDIA (1) - VISCOELASTIC SOLID (2)
          ! ----------------------------------------------

          case (fbem_poroelastic)
            do il=1,problem%n
              row=node(sn_col)%row(il,eq_index)
              do ik=1,problem%n
                do kn_int=1,se_int_n_nodes
                  sn_int=element(se_int)%node(kn_int)
                  select case (node(sn_int)%ctype(ik,1))
                    ! Global axes
                    ! 0 - perfect bonding - u_k^{1} and t_k^{1} are active since:
                    ! t_k^{(2)}=-t_k^{(1)}-tau^(1)n_k^(1)
                    ! u_k^{(2)}= u_k^{(1)}
                    case (0)
                      ! u_k for region 2: u_k^(2)= u_k^(1)
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                      ! t_k for region 2: t_k^{(2)}=-t_k^{(1)}-tau^(1)n_k^(1)
                      ! t_k^(1)
                      col=node(sn_int)%col(problem%n+1+ik,1)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)
                      ! tau^(1)
                      col=node(sn_int)%col(0,1)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                    ! Global axes
                    ! 1 - perfect debonding - u_k^{1} and u_k^{2} are active since:
                    ! t_k^{(1)}=-tau^(1)n_k^(1)
                    ! t_k^{(2)}=0
                    case (1)
                      ! u_k for region 2
                      col=node(sn_int)%col(ik,2)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                    ! Global axes
                    ! 2 - partial bonding   - u_k^{1} and u_k^{2} are active since:
                    ! t_k^{(1)}=-(u_k^{1}-u_k^{2})*KT-tau^(1)n_k^(1)
                    ! t_k^{(2)}=(u_k^{1}-u_k^{2})*KT
                    case (2)
                      ! KT = K + i*omega*C - omega**2*M
                      KT=node(sn_int)%cvalue_c(ik,1,1)+c_im*omega*node(sn_int)%cvalue_c(ik,2,1)&
                        -omega**2*node(sn_int)%cvalue_c(ik,3,1)
                      ! u_k for region 2
                      col=node(sn_int)%col(ik,2)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                      ! t_k for region 2: t_k^{(2)}=(u_k^{1}-u_k^{2})*KT
                      ! u_k^{2}
                      col=node(sn_int)%col(ik,2)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*KT
                      ! u_k^{1}
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)*KT
                    ! Local axes
                    ! In all cases, u_k^{1}, u_k^{2} and t_k^{2} are active, and t_k^{(1)}=-t_k^{(2)}-tau^(1)n_k^(1).
                    case (3,4,5)
                      ! u_k for region 2
                      col=node(sn_int)%col(ik,2)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                      ! t_k for region 2
                      col=node(sn_int)%col(problem%n+ik,2)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
                  end select
                end do
              end do
            end do

        end select
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
                A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
                ! If u_k is unknown
                if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                  col=node(sn_int_fe)%col(ik,1)
                  A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                ! If u_k is known
                else
                  b_c(row,1)=b_c(row,1)-hp(kn_int,il,ik)*node(sn_int_fe)%cvalue_c(ik,1,1)
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
!                A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
!                ! Note: h^+ + h^- = 0 always except when it contains the free-term, so this gives the required contribution
!                ! If u_k=u_k_{+}=u_k_{-} is unknown
!                if (node(sn_int_fe)%ctype(ik,1).eq.1) then
!                  col=node(sn_int_fe)%col(ik,1)
!                  A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)+hm(kn_int,il,ik)
!                ! If u_k=u_k_{+}=u_k_{-} is known
!                else
!                  b_c(row,1)=b_c(row,1)-(hp(kn_int,il,ik)+hm(kn_int,il,ik))*node(sn_int_fe)%cvalue_c(ik,1,1)
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
                A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
                ! t_k for face - is unknown
                col=node(sn_int)%col(problem%n+ik,2)
                A_c(row,col)=A_c(row,col)-gm(kn_int,il,ik)
                ! If u_k=u_k_{+}=u_k_{-} is unknown
                if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                  col=node(sn_int_fe)%col(ik,1)
                  A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)+hm(kn_int,il,ik)
                ! If u_k=u_k_{+}=u_k_{-} is known
                else
                  b_c(row,1)=b_c(row,1)-(hp(kn_int,il,ik)+hm(kn_int,il,ik))*node(sn_int_fe)%cvalue_c(ik,1,1)
                end if
              end do
            end do
          end do

      end select

    ! =================
    ! BE-FE-BE BOUNDARY
    ! =================

    case (fbem_boundary_coupling_be_fe_be)
      if (sb_int_reversion.eqv.(.false.)) then

        do il=1,problem%n
          row=node(sn_col)%row(il,eq_index)
          do ik=1,problem%n
            do kn_int=1,se_int_n_nodes
              sn_int=element(se_int)%node(kn_int)
              sn_int_fe=element(se_int)%element_node(kn_int)
              ! t_k^(1) unknown
              col=node(sn_int)%col(problem%n+ik,1)
              A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
              ! u_k^(1) = u_k^(FE)
              ! u_k^(FE) is unknown
              if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                col=node(sn_int_fe)%col(ik,1)
                A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
              ! u_k^(FE) is prescribed
              else
                b_c(row,1)=b_c(row,1)-hp(kn_int,il,ik)*node(sn_int_fe)%cvalue_c(ik,1,1)
              end if
            end do
          end do
        end do

      else

        do il=1,problem%n
          row=node(sn_col)%row(il,eq_index)
          do ik=1,problem%n
            do kn_int=1,se_int_n_nodes
              sn_int=element(se_int)%node(kn_int)
              sn_int_fe=element(se_int)%element_node(kn_int)
              ! t_k^(2) is unknown
              col=node(sn_int)%col(problem%n+ik,2)
              A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
              ! u_k^(2) = u_k^(FE)
              ! u_k^(FE) is unknown
              if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                col=node(sn_int_fe)%col(ik,1)
                A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
              ! u_k^(FE) is prescribed
              else
                b_c(row,1)=b_c(row,1)-hp(kn_int,il,ik)*node(sn_int_fe)%cvalue_c(ik,1,1)
              end if
            end do
          end do
        end do

      end if

  end select

  ! ============================
  ! ASSEMBLE INCIDENT WAVE FIELD
  ! ============================

  if (region(kr)%n_incidentfields.gt.0) then
    select case (boundary(sb_int)%class)

      ! -----------------
      ! ORDINARY BOUNDARY
      ! -----------------

      case (fbem_boundary_class_ordinary)
        do il=1,problem%n
          row=node(sn_col)%row(il,eq_index)
          do ik=1,problem%n
            do kn_int=1,se_int_n_nodes
              b_c(row,1)=b_c(row,1)+hp(kn_int,il,ik)*u_inc(ik,kn_int)-gp(kn_int,il,ik)*t_inc(ik,kn_int)
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
              ! Face +
              b_c(row,1)=b_c(row,1)+hp(kn_int,il,ik)*u_inc(ik,kn_int)-gp(kn_int,il,ik)*t_inc(ik,kn_int)
              ! Face -
              b_c(row,1)=b_c(row,1)+hm(kn_int,il,ik)*u_inc(ik,kn_int)+gm(kn_int,il,ik)*t_inc(ik,kn_int)
            end do
          end do
        end do

    end select
  end if

end subroutine assemble_bem_harela_equation
