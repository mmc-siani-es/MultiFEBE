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

subroutine assemble_bem_harpor_equation(omega,kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,eq_index,&
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
  complex(kind=real64) :: hp(se_int_n_nodes,0:problem%n,0:problem%n)
  complex(kind=real64) :: gp(se_int_n_nodes,0:problem%n,0:problem%n)
  complex(kind=real64) :: hm(se_int_n_nodes,0:problem%n,0:problem%n)
  complex(kind=real64) :: gm(se_int_n_nodes,0:problem%n,0:problem%n)
  complex(kind=real64) :: u_inc(0:problem%n,se_int_n_nodes)
  complex(kind=real64) :: t_inc(0:problem%n,se_int_n_nodes)

  ! Local variables
  integer              :: il, ik, ikc
  integer              :: kn_int, sn_int, sn_int_fe
  integer              :: row, col
  real(kind=real64)    :: phi, phi1, phi2
  real(kind=real64)    :: d_phi1_phi2, d_phi2_phi1, p_1_phi1_phi2, p_1_phi2_phi1
  real(kind=real64)    :: ctephi
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
          do il=0,problem%n
            row=node(sn_col)%row(il,eq_index)
            !
            ! Fluid phase variables
            !
            ik=0
            do kn_int=1,se_int_n_nodes
              sn_int=element(se_int)%node(kn_int)
              select case (node(sn_int)%ctype(ik,1))
                !
                ! open pore / permeable
                !
                ! tau known, Un unknown
                case (0)
                  col=node(sn_int)%col(problem%n+1+ik,1)
                  A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
                  b_c(row,1)=b_c(row,1)-hp(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,1)
                ! Un known, tau unknown
                case (1)
                  col=node(sn_int)%col(ik,1)
                  A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                  b_c(row,1)=b_c(row,1)+gp(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,1)
                !
                ! close pore / impermeable
                !
                ! Un=u_k·n_k known, tau unknown
                case (2)
                  ! tau unknown
                  col=node(sn_int)%col(ik,1)
                  A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                  ! Un=u_k·n_k
                  if (sb_int_reversion.eqv.(.false.)) then
                    ! n_fn is positive
                    do ikc=1,problem%n
                      select case (node(sn_int)%ctype(ikc,1))
                        ! u_k known
                        case (4)
                          b_c(row,1)=b_c(row,1)+gp(kn_int,il,ik)*element(se_int)%n_fn(ikc,kn_int)*node(sn_int)%cvalue_c(ikc,1,1)
                        ! u_k unknown
                        case (5,50,6,7)
                          col=node(sn_int)%col(ikc,1)
                          A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)*element(se_int)%n_fn(ikc,kn_int)
                      end select
                    end do
                  else
                    ! n_fn is negative
                    do ikc=1,problem%n
                      select case (node(sn_int)%ctype(ikc,1))
                        ! u_k known
                        case (4)
                          b_c(row,1)=b_c(row,1)-gp(kn_int,il,ik)*element(se_int)%n_fn(ikc,kn_int)*node(sn_int)%cvalue_c(ikc,1,1)
                        ! u_k unknown
                        case (5,50,6,7)
                          col=node(sn_int)%col(ikc,1)
                          A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*element(se_int)%n_fn(ikc,kn_int)
                      end select
                    end do
                  end if
              end select
            end do
            !
            ! Solid skeleton variables
            !
            do ik=1,problem%n
              do kn_int=1,se_int_n_nodes
                sn_int=element(se_int)%node(kn_int)
                select case (node(sn_int)%ctype(ik,1))
                  !
                  ! open pore / permeable
                  !
                  ! u_k known, t^{hybrid}_k unknown
                  case (0)
                    col=node(sn_int)%col(problem%n+1+ik,1)
                    A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
                    b_c(row,1)=b_c(row,1)-hp(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,1)
                  ! t^{solid skeleton}_k known, u_k unknown
                  case (1)
                    ! Assemble u_k in A
                    col=node(sn_int)%col(ik,1)
                    A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                    ! Assemble known t^{solid skeleton}_k
                    b_c(row,1)=b_c(row,1)+gp(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,1)
                  ! t_k known (normal pressure), u_k unknown
                  case (10)
                    col=node(sn_int)%col(ik,1)
                    A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                    if (sb_int_reversion.eqv.(.false.)) then
                      b_c(row,1)=b_c(row,1)+gp(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,1)*node(sn_int)%n_fn(ik)
                    else
                      b_c(row,1)=b_c(row,1)-gp(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,1)*node(sn_int)%n_fn(ik)
                    end if
                  ! u_k unknown, t^{hybrid}_k unknown
                  case (2,3)
                    col=node(sn_int)%col(            ik,1)
                    A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                    col=node(sn_int)%col(problem%n+1+ik,1)
                    A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
                  !
                  ! close pore / impermeable
                  !
                  ! u_k known, t^{hybrid}_k unknown
                  case (4)
                    col=node(sn_int)%col(problem%n+1+ik,1)
                    A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
                    b_c(row,1)=b_c(row,1)-hp(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,1)
                  ! t^{total}_k known, u_k unknown
                  case (5)
                    ! u_k unknown
                    col=node(sn_int)%col(ik,1)
                    A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                    ! Assemble known t^{total}_k
                    b_c(row,1)=b_c(row,1)+gp(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,1)
                    ! Assemble tau for close pore / impermeable
                    ! tau unknown
                    col=node(sn_int)%col(0,1)
                    ! n_fn is positive
                    if (sb_int_reversion.eqv.(.false.)) then
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                    ! n_fn is negative
                    else
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                    end if
                  ! t_k+tau n_k=Pn_k known (normal pressure) tau unknown, u_k unknown
                  case (50)
                    ! Assemble u_k
                    col=node(sn_int)%col(ik,1)
                    A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                    ! Assemble t_k = Pn_k - tau n_k
                    ! P n_k
                    if (sb_int_reversion.eqv.(.false.)) then
                      b_c(row,1)=b_c(row,1)+gp(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,1)*node(sn_int)%n_fn(ik)
                    else
                      b_c(row,1)=b_c(row,1)-gp(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,1)*node(sn_int)%n_fn(ik)
                    end if
                    ! tau
                    col=node(sn_int)%col(0,1)
                    if (sb_int_reversion.eqv.(.false.)) then
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*node(sn_int)%n_fn(ik)
                    else
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)*node(sn_int)%n_fn(ik)
                    end if
                  ! u_k unknown, t^{hybrid}_k unknown
                  case (6,7)
                    col=node(sn_int)%col(            ik,1)
                    A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                    col=node(sn_int)%col(problem%n+1+ik,1)
                    A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
                end select
              end do
            end do
          end do

        ! -------------------
        ! CRACK-LIKE BOUNDARY
        ! -------------------

        case (fbem_boundary_class_cracklike)
          do il=0,problem%n
            row=node(sn_col)%row(il,eq_index)
            do ik=0,problem%n
              do kn_int=1,se_int_n_nodes
                sn_int=element(se_int)%node(kn_int)
                !
                ! Face +
                !
                select case (node(sn_int)%ctype(ik,1))
                  ! If k==0: tau^+ known, Un^+ unknown
                  ! If k>=1: u_k^+ known, t_k^+ unknown
                  case (0)
                    col=node(sn_int)%col(problem%n+1+ik,1)
                    A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
                    b_c(row,1)=b_c(row,1)-hp(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,1)
                  ! If k==0: Un^+ known, tau^+ unknown
                  ! If k>=1: t_k^+ known, u_k^+ unknown
                  case (1)
                    col=node(sn_int)%col(ik,1)
                    A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                    b_c(row,1)=b_c(row,1)+gp(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,1)
                  ! close pore / impermeable
                  ! If k==0: Un^+=u_k^+·n_k^+ known, tau^+ unknown
                  ! If k>=1: u_k^+ known, t_k^+ unknown
                  case (2)
                    ! Un^+=u_k^+·n_k^+ known, tau^+ unknown
                    if (ik.eq.0) then
                      ! tau^+ unknown
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                      ! Un^+=u_k^+·n_k^+
                      if (sb_int_reversion.eqv.(.false.)) then
                        ! n_fn is positive
                        do ikc=1,problem%n
                          select case (node(sn_int)%ctype(ikc,1))
                            ! u_k^+ known
                            case (2)
                              b_c(row,1)=b_c(row,1)+gp(kn_int,il,ik)*element(se_int)%n_fn(ikc,kn_int)*node(sn_int)%cvalue_c(ikc,1,1)
                            ! u_k^+ unknown
                            case (3)
                              col=node(sn_int)%col(ikc,1)
                              A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)*element(se_int)%n_fn(ikc,kn_int)
                          end select
                        end do
                      else
                        ! n_fn is negative
                        do ikc=1,problem%n
                          select case (node(sn_int)%ctype(ikc,1))
                            ! u_k^+ known
                            case (2)
                              b_c(row,1)=b_c(row,1)-gp(kn_int,il,ik)*element(se_int)%n_fn(ikc,kn_int)*node(sn_int)%cvalue_c(ikc,1,1)
                            ! u_k^+ unknown
                            case (3)
                              col=node(sn_int)%col(ikc,1)
                              A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*element(se_int)%n_fn(ikc,kn_int)
                          end select
                        end do
                      end if
                    ! u_k^+ known, t_k^+ unknown
                    else
                      col=node(sn_int)%col(problem%n+1+ik,1)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
                      b_c(row,1)=b_c(row,1)-hp(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,1)
                    end if
                  ! close pore / impermeable
                  ! Only for k>=1: t_k^++tau^+·n_k^+ known, u_k^+ unknown
                  case (3)
                    ! u_k^+ unknown
                    col=node(sn_int)%col(ik,1)
                    A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                    ! T^+ known
                    b_c(row,1)=b_c(row,1)+gp(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,1)
                    ! tau^+ unknown
                    col=node(sn_int)%col(0,1)
                    ! n_fn is positive
                    if (sb_int_reversion.eqv.(.false.)) then
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                    ! n_fn is negative
                    else
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                    end if
                  ! Fluid(inviscid/incompressible)-filled permeable crack
                  ! If k==0: tau^+ unknown, Un^+ unknown
                  ! If k>=1: u^+ unknown, t_k^+=(1/phi-1)tau^+·n_k^+ known
                  case (4)
                    if (ik.eq.0) then
                      ! tau^+ unknown
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                      ! Un^+ unknown
                      col=node(sn_int)%col(problem%n+1+ik,1)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
                    else
                      ! u_k^+ unknown
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                      ! t_k^+=(1/phi-1)tau^+·n_k^+ known
                      ! Calculate 1/phi-1
                      phi=region(boundary(sb_int)%region(1))%property_r(8)
                      ctephi=1.0d0/phi-1.0d0
                      ! tau^+
                      col=node(sn_int)%col(0,1)
                      ! n_fn is positive
                      if (sb_int_reversion.eqv.(.false.)) then
                        A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)*ctephi*element(se_int)%n_fn(ik,kn_int)
                      ! n_fn is negative
                      else
                        A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*ctephi*element(se_int)%n_fn(ik,kn_int)
                      end if
                    end if
                  ! Fluid(inviscid/incompressible)-filled impermeable crack
                  ! If k==0: tau^+ unknown, Un^+=u^+·n^+ known
                  ! If k>=1: u_k^+ unknown, t_k^+ unknown
                  case (5)
                    if (ik.eq.0) then
                      ! tau^+ unknown
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                      ! Un^+=u_k^+·n_k^+
                      ! Boundary orientation is positive
                      if (sb_int_reversion.eqv.(.false.)) then
                        do ikc=1,problem%n
                          col=node(sn_int)%col(ikc,1)
                          A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)*element(se_int)%n_fn(ikc,kn_int)
                        end do
                      ! Boundary orientation is negative
                      else
                        do ikc=1,problem%n
                          col=node(sn_int)%col(ikc,1)
                          A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*element(se_int)%n_fn(ikc,kn_int)
                        end do
                      end if
                    else
                      ! u_k^+ unknown
                      col=node(sn_int)%col(            ik,1)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                      ! t_k^+ unknown
                      col=node(sn_int)%col(problem%n+1+ik,1)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
                    end if
                end select
                !
                ! Face -
                !
                select case (node(sn_int)%ctype(ik,2))
                  ! If k==0: tau^- known, w^- unknown
                  ! If k>=1: u_k^- known, t_k^- unknown
                  case (0)
                    col=node(sn_int)%col(problem%n+1+ik,2)
                    A_c(row,col)=A_c(row,col)-gm(kn_int,il,ik)
                    b_c(row,1)=b_c(row,1)-hm(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,2)
                  ! If k==0: w^- known, tau^- unknown
                  ! If k>=1: t_k^- known, u_k^- unknown
                  case (1)
                    col=node(sn_int)%col(ik,2)
                    A_c(row,col)=A_c(row,col)+hm(kn_int,il,ik)
                    b_c(row,1)=b_c(row,1)+gm(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,2)
                  ! close pore / impermeable
                  ! If k==0: Un^-=u_k^-·n_k^- known, tau^- unknown
                  ! If k>=1: u_k^- known, t_k^- unknown
                  case (2)
                    ! Un^-=u_k^-·n_k^- known, tau^- unknown
                    if (ik.eq.0) then
                      ! tau^- unknown
                      col=node(sn_int)%col(ik,2)
                      A_c(row,col)=A_c(row,col)+hm(kn_int,il,ik)
                      ! Un^-=u_k^-·n_k^-
                      if (sb_int_reversion.eqv.(.false.)) then
                        ! n_fn is negative
                        do ikc=1,problem%n
                          select case (node(sn_int)%ctype(ikc,2))
                            ! u_k^- known
                            case (2)
                              b_c(row,1)=b_c(row,1)-gm(kn_int,il,ik)*element(se_int)%n_fn(ikc,kn_int)*node(sn_int)%cvalue_c(ikc,1,2)
                            ! u_k^- unknown
                            case (3)
                              col=node(sn_int)%col(ikc,2)
                              A_c(row,col)=A_c(row,col)+gm(kn_int,il,ik)*element(se_int)%n_fn(ikc,kn_int)
                          end select
                        end do
                      else
                        ! n_fn is positive
                        do ikc=1,problem%n
                          select case (node(sn_int)%ctype(ikc,2))
                            ! u_k^- known
                            case (2)
                              b_c(row,1)=b_c(row,1)+gm(kn_int,il,ik)*element(se_int)%n_fn(ikc,kn_int)*node(sn_int)%cvalue_c(ikc,1,2)
                            ! u_k^- unknown
                            case (3)
                              col=node(sn_int)%col(ikc,2)
                              A_c(row,col)=A_c(row,col)-gm(kn_int,il,ik)*element(se_int)%n_fn(ikc,kn_int)
                          end select
                        end do
                      end if
                    ! u_k^- known, t_k^- unknown
                    else
                      col=node(sn_int)%col(problem%n+1+ik,2)
                      A_c(row,col)=A_c(row,col)-gm(kn_int,il,ik)
                      b_c(row,1)=b_c(row,1)-hm(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,2)
                    end if
                  ! close pore / impermeable
                  ! Only for k>=1: T_k=t_k^-+tau^-·n_k^- known, u_k^- unknown
                  case (3)
                    ! u_k^- unknown
                    col=node(sn_int)%col(ik,2)
                    A_c(row,col)=A_c(row,col)+hm(kn_int,il,ik)
                    ! t^- known
                    b_c(row,1)=b_c(row,1)+gm(kn_int,il,ik)*node(sn_int)%cvalue_c(ik,1,2)
                    ! tau^- unknown
                    col=node(sn_int)%col(0,2)
                    ! n_fn is negative
                    if (sb_int_reversion.eqv.(.false.)) then
                      A_c(row,col)=A_c(row,col)-gm(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                    ! n_fn is positive
                    else
                      A_c(row,col)=A_c(row,col)+gm(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                    end if
                  ! Fluid(inviscid/incompressible)-filled permeable crack
                  ! If k==0: tau^- and Un^- known in function of other variables
                  ! If k>=1: u^- unknown, t_k^-=-(1/phi-1)tau^+·n_k^+ known
                  case (4)
                    if (ik.eq.0) then
                      ! tau^-=tau^+
                      ! tau^+ unknown
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)+hm(kn_int,il,ik)
                      ! Un^-= - Un^+ - (1/phi-1)u^+·n^+ + (1/phi-1)u^-·n^+
                      ! Calculate 1/phi-1
                      phi=region(boundary(sb_int)%region(1))%property_r(8)
                      ctephi=1.0d0/phi-1.0d0
                      ! Un^+ unknown
                      col=node(sn_int)%col(problem%n+1+ik,1)
                      A_c(row,col)=A_c(row,col)+gm(kn_int,il,ik)
                      ! u^+ and u^- unknowns
                      ! Boundary orientation is positive
                      if (sb_int_reversion.eqv.(.false.)) then
                        do ikc=1,problem%n
                          ! u_ikc^+
                          col=node(sn_int)%col(ikc,1)
                          A_c(row,col)=A_c(row,col)+gm(kn_int,il,ik)*ctephi*element(se_int)%n_fn(ikc,kn_int)
                          ! u_ikc^-
                          col=node(sn_int)%col(ikc,2)
                          A_c(row,col)=A_c(row,col)-gm(kn_int,il,ik)*ctephi*element(se_int)%n_fn(ikc,kn_int)
                        end do
                      else
                        do ikc=1,problem%n
                          ! u_ikc^+
                          col=node(sn_int)%col(ikc,1)
                          A_c(row,col)=A_c(row,col)-gm(kn_int,il,ik)*ctephi*element(se_int)%n_fn(ikc,kn_int)
                          ! u_ikc^-
                          col=node(sn_int)%col(ikc,2)
                          A_c(row,col)=A_c(row,col)+gm(kn_int,il,ik)*ctephi*element(se_int)%n_fn(ikc,kn_int)
                        end do
                      end if
                    else
                      ! u_k^- unknown
                      col=node(sn_int)%col(ik,2)
                      A_c(row,col)=A_c(row,col)+hm(kn_int,il,ik)
                      ! t_k^-=-(1/phi-1)tau^+·n_k^+ known
                      ! Calculate 1/phi-1
                      phi=region(boundary(sb_int)%region(1))%property_r(8)
                      ctephi=1.0d0/phi-1.0d0
                      ! tau^+
                      col=node(sn_int)%col(0,1)
                      ! n_fn is positive
                      if (sb_int_reversion.eqv.(.false.)) then
                        A_c(row,col)=A_c(row,col)+gm(kn_int,il,ik)*ctephi*element(se_int)%n_fn(ik,kn_int)
                      ! n_fn is negative
                      else
                        A_c(row,col)=A_c(row,col)-gm(kn_int,il,ik)*ctephi*element(se_int)%n_fn(ik,kn_int)
                      end if
                    end if
                  ! Fluid(inviscid/incompressible)-filled impermeable crack
                  ! If k==0: tau^- unknown, Un^-=u^-·n^- known
                  ! If k>=1: u^- unknown, t_k^-=(tau^--tau^+)n_k^+-t_k^+ known
                  case (5)
                    if (ik.eq.0) then
                      ! tau^- unknown
                      col=node(sn_int)%col(ik,2)
                      A_c(row,col)=A_c(row,col)+hm(kn_int,il,ik)
                      ! Un^-=-u_k^-·n_k^+
                      ! Boundary orientation is positive
                      if (sb_int_reversion.eqv.(.false.)) then
                        ! n_fn is positive
                        do ikc=1,problem%n
                          col=node(sn_int)%col(ikc,2)
                          A_c(row,col)=A_c(row,col)+gm(kn_int,il,ik)*element(se_int)%n_fn(ikc,kn_int)
                        end do
                      ! Boundary orientation is negative
                      else
                        do ikc=1,problem%n
                          col=node(sn_int)%col(ikc,2)
                          A_c(row,col)=A_c(row,col)-gm(kn_int,il,ik)*element(se_int)%n_fn(ikc,kn_int)
                        end do
                      end if
                    else
                      ! u_k^- unknown
                      col=node(sn_int)%col(ik,2)
                      A_c(row,col)=A_c(row,col)+hm(kn_int,il,ik)
                      ! t_k^-=tau^-*n_k^+ - tau^+*n_k^+-t_k^+
                      ! tau^+ and tau^- unknowns
                      ! Boundary orientation is positive
                      if (sb_int_reversion.eqv.(.false.)) then
                        ! tau^+ unknown
                        col=node(sn_int)%col(0,1)
                        A_c(row,col)=A_c(row,col)+gm(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                        ! tau^- unknown
                        col=node(sn_int)%col(0,2)
                        A_c(row,col)=A_c(row,col)-gm(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                      ! Boundary orientation is negative
                      else
                        ! tau^+ unknown
                        col=node(sn_int)%col(0,1)
                        A_c(row,col)=A_c(row,col)-gm(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                        ! tau^- unknown
                        col=node(sn_int)%col(0,2)
                        A_c(row,col)=A_c(row,col)+gm(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                      end if
                      ! t_k^+ unknown
                      col=node(sn_int)%col(problem%n+1+ik,1)
                      A_c(row,col)=A_c(row,col)+gm(kn_int,il,ik)
                    end if
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

          ! ------------------------------------------
          ! POROELASTIC MEDIA (1) - INVISCID FLUID (2)
          ! ------------------------------------------

          case (fbem_potential)
            phi=region(boundary(sb_int)%region(1))%property_r(8)
            do il=0,problem%n
              row=node(sn_col)%row(il,eq_index)
              do kn_int=1,se_int_n_nodes
                sn_int=element(se_int)%node(kn_int)
                ! Fluid phase variables
                ! tau^(1) unknown
                col=node(sn_int)%col(0,1)
                A_c(row,col)=A_c(row,col)+hp(kn_int,il,0)
                ! Un^(1) is assembled differently depending on the B.C.
                select case (node(sn_int)%ctype(1,1))
                  ! Perfectly permeable
                  case (0)
                    ! Un^(1) unknown
                    col=node(sn_int)%col(problem%n+1,1)
                    A_c(row,col)=A_c(row,col)-gp(kn_int,il,0)
                  ! Perfectly impermeable
                  case (1)
                    ! Un^(1)=u^(1)·n^(1)
                    do ik=1,problem%n
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,0)*element(se_int)%n_fn(ik,kn_int)
                    end do
                end select
              end do
              ! Solid skeleton variables
              do ik=1,problem%n
                do kn_int=1,se_int_n_nodes
                  sn_int=element(se_int)%node(kn_int)
                  ! u_k^(1) unknown
                  col=node(sn_int)%col(            ik,1)
                  A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                  ! t_k^(1) is assembled differently depending on the B.C.
                  select case (node(sn_int)%ctype(1,1))
                    ! Perfectly permeable
                    case (0)
                      ! t_k^(1)=(1-phi^(1))/phi^(1)*tau^(1)*n_k^(1)
                      col=node(sn_int)%col(0,1)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)*(1.0d0-phi)/phi*element(se_int)%n_fn(ik,kn_int)
                    ! Perfectly impermeable
                    case (1)
                      ! t_k^(1)=-p^(2)*n_k^(1)-tau^(1)*n_k^(1)
                      ! tau^(1)
                      col=node(sn_int)%col(1,2)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                      ! tau^(1)
                      col=node(sn_int)%col(0,1)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                  end select
                end do
              end do
            end do

          ! ----------------------------------------------
          ! POROELASTIC MEDIA (1) - VISCOELASTIC SOLID (2)
          ! ----------------------------------------------

          case (fbem_viscoelastic)
            do il=0,problem%n
              row=node(sn_col)%row(il,eq_index)
              do kn_int=1,se_int_n_nodes
                sn_int=element(se_int)%node(kn_int)
                ! Fluid phase variables
                ! tau^(1) unknown
                col=node(sn_int)%col(0,1)
                A_c(row,col)=A_c(row,col)+hp(kn_int,il,0)
                ! Un^(1)=u^(1)·n^(1)
                do ik=1,problem%n
                  ! u^(1) unknown
                  col=node(sn_int)%col(ik,1)
                  A_c(row,col)=A_c(row,col)-gp(kn_int,il,0)*element(se_int)%n_fn(ik,kn_int)
                end do
              end do
              ! Solid skeleton variables
              do ik=1,problem%n
                do kn_int=1,se_int_n_nodes
                  sn_int=element(se_int)%node(kn_int)
                  ! u_k^(1) unknown
                  col=node(sn_int)%col(ik,1)
                  A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                  ! t_k^(1) is assembled differently depending on the B.C.
                  select case (node(sn_int)%ctype(ik,1))
                    ! 0 - perfect bonding - t_k^(1) unknown
                    case (0)
                      col=node(sn_int)%col(problem%n+1+ik,1)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
                    ! 1 - perfect debonding - t_k^{(1)}=-tau^(1)n_k^(1)
                    case (1)
                      col=node(sn_int)%col(0,1)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                    ! 2 - partial bonding - t_k^{(1)}=-(u_k^{1}-u_k^{2})*KT-tau^(1)n_k^(1)
                    case (2)
                      ! KT = K + i*omega*C - omega**2*M
                      KT=node(sn_int)%cvalue_c(ik,1,1)+c_im*omega*node(sn_int)%cvalue_c(ik,2,1)&
                        -omega**2*node(sn_int)%cvalue_c(ik,3,1)
                      ! u_k^{1}
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*KT
                      ! u_k^{2}
                      col=node(sn_int)%col(ik,2)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)*KT
                      ! tau^(1)
                      col=node(sn_int)%col(0,1)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                    ! 3,4,5 - local axes - t_k^{(1)}=-t_k^{(2)}-tau^(1)n_k^(1)
                    case (3,4,5)
                      ! t_k^{2}
                      col=node(sn_int)%col(problem%n+ik,2)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)
                      ! tau^(1)
                      col=node(sn_int)%col(0,1)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                  end select
                end do
              end do
            end do

          ! ---------------------------------------------
          ! POROELASTIC MEDIA (1) - POROELASTIC MEDIA (2)
          ! ---------------------------------------------

          case (fbem_poroelastic)

            do il=0,problem%n
              row=node(sn_col)%row(il,eq_index)
              do kn_int=1,se_int_n_nodes
                sn_int=element(se_int)%node(kn_int)
                ! Fluid phase variables
                ! tau^(1) unknown
                col=node(sn_int)%col(0,1)
                A_c(row,col)=A_c(row,col)+hp(kn_int,il,0)
                ! Un^(1) is assembled differently depending on the B.C.
                select case (node(sn_int)%ctype(1,1))
                  ! Perfectly permeable or partially impermeable
                  case (0,2)
                    ! Un^(1) unknown
                    col=node(sn_int)%col(problem%n+1,1)
                    A_c(row,col)=A_c(row,col)-gp(kn_int,il,0)
                  ! Perfectly impermeable
                  case (1)
                    ! Un^(1)=u^(1)·n^(1)
                    do ik=1,problem%n
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,0)*element(se_int)%n_fn(ik,kn_int)
                    end do
                end select
              end do
              ! Solid skeleton variables
              do ik=1,problem%n
                do kn_int=1,se_int_n_nodes
                  sn_int=element(se_int)%node(kn_int)
                  ! u_k^(1) unknown
                  col=node(sn_int)%col(            ik,1)
                  A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                  ! t_k^(1) unknown
                  col=node(sn_int)%col(problem%n+1+ik,1)
                  A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
                end do
              end do
            end do

        end select
      else
        select case (region(boundary(sb_int)%region(1))%type)

          ! ------------------------------------------
          ! INVISCID FLUID (1) - POROELASTIC MEDIA (2)
          ! ------------------------------------------

          case (fbem_potential)
            phi=region(boundary(sb_int)%region(2))%property_r(8)
            do il=0,problem%n
              row=node(sn_col)%row(il,eq_index)
              do kn_int=1,se_int_n_nodes
                sn_int=element(se_int)%node(kn_int)
                ! Fluid phase variables
                ! tau^(2) unknown
                col=node(sn_int)%col(0,2)
                A_c(row,col)=A_c(row,col)+hp(kn_int,il,0)
                ! Un^(2) is assembled differently depending on the B.C.
                select case (node(sn_int)%ctype(1,1))
                  ! Perfectly permeable
                  case (0)
                    ! Un^(2) unknown
                    col=node(sn_int)%col(problem%n+1,2)
                    A_c(row,col)=A_c(row,col)-gp(kn_int,il,0)
                  ! Perfectly impermeable
                  case (1)
                    ! Un^(2)=-u^(2)·n^(1)
                    do ik=1,problem%n
                      col=node(sn_int)%col(ik,2)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,0)*element(se_int)%n_fn(ik,kn_int)
                    end do
                end select
              end do
              ! Solid skeleton variables
              do ik=1,problem%n
                do kn_int=1,se_int_n_nodes
                  sn_int=element(se_int)%node(kn_int)
                  ! u_k^(2) unknown
                  col=node(sn_int)%col(            ik,2)
                  A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                  ! t_k^(2) is assembled differently depending on the B.C.
                  select case (node(sn_int)%ctype(1,1))
                    ! Perfectly permeable
                    case (0)
                      ! t_k^(2)=-(1-phi^(2))/phi^(2)*tau^(2)*n_k^(1)
                      col=node(sn_int)%col(0,2)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*(1.0d0-phi)/phi*element(se_int)%n_fn(ik,kn_int)
                    ! Perfectly impermeable
                    case (1)
                      ! t_k^(2)=p^(1)*n_k^(1)+tau^(2)*n_k^(1)
                      ! tau^(1)
                      col=node(sn_int)%col(1,1)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                      ! tau^(1)
                      col=node(sn_int)%col(0,2)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                  end select
                end do
              end do
            end do

          ! ----------------------------------------------
          ! VISCOELASTIC SOLID (1) - POROELASTIC MEDIA (2)
          ! ----------------------------------------------

          case (fbem_viscoelastic)
            do il=0,problem%n
              row=node(sn_col)%row(il,eq_index)
              do kn_int=1,se_int_n_nodes
                sn_int=element(se_int)%node(kn_int)
                ! Fluid phase variables
                ! tau^(2) unknown
                col=node(sn_int)%col(0,2)
                A_c(row,col)=A_c(row,col)+hp(kn_int,il,0)
                ! Un^(2)=-u^(2)·n^(1)
                do ik=1,problem%n
                  ! u^(2) unknown
                  col=node(sn_int)%col(ik,2)
                  A_c(row,col)=A_c(row,col)+gp(kn_int,il,0)*element(se_int)%n_fn(ik,kn_int)
                end do
              end do
              ! Solid skeleton variables
              do ik=1,problem%n
                do kn_int=1,se_int_n_nodes
                  sn_int=element(se_int)%node(kn_int)
                  ! u_k^(2) unknown
                  col=node(sn_int)%col(ik,2)
                  A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                  ! t_k^(2) is assembled differently depending on the B.C.
                  select case (node(sn_int)%ctype(ik,1))
                    ! 0 - perfect bonding - t_k^(2) unknown
                    case (0)
                      col=node(sn_int)%col(problem%n+1+ik,2)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
                    ! 1 - perfect debonding - t_k^{(2)}=tau^(2)n_k^(1)
                    case (1)
                      col=node(sn_int)%col(0,2)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                    ! 2 - partial bonding - t_k^{(2)}=-(u_k^{2}-u_k^{1})*KT+tau^(2)n_k^(1)
                    case (2)
                      ! KT = K + i*omega*C - omega**2*M
                      KT=node(sn_int)%cvalue_c(ik,1,1)+c_im*omega*node(sn_int)%cvalue_c(ik,2,1)&
                        -omega**2*node(sn_int)%cvalue_c(ik,3,1)
                      ! u_k^{2}
                      col=node(sn_int)%col(ik,2)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*KT
                      ! u_k^{1}
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)*KT
                      ! tau^(2)
                      col=node(sn_int)%col(0,2)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                    ! 3,4,5 - local axes - t_k^{(2)}=-t_k^{(1)}+tau^(2)n_k^(1)
                    case (3,4,5)
                      ! t_k^{1}
                      col=node(sn_int)%col(problem%n+ik,1)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)
                      ! tau^(2)
                      col=node(sn_int)%col(0,2)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                  end select
                end do
              end do
            end do

          ! ---------------------------------------------
          ! POROELASTIC MEDIA (1) - POROELASTIC MEDIA (2)
          ! ---------------------------------------------

          case (fbem_poroelastic)

            ! Auxiliary constants
            phi1=region(boundary(sb_int)%region(1))%property_r(8)
            phi2=region(boundary(sb_int)%region(2))%property_r(8)
            d_phi1_phi2=phi1/phi2
            d_phi2_phi1=1.0d0/d_phi1_phi2
            p_1_phi2_phi1=1.0d0-d_phi2_phi1
            p_1_phi1_phi2=1.0d0-d_phi1_phi2
            do il=0,problem%n
              row=node(sn_col)%row(il,eq_index)
              do kn_int=1,se_int_n_nodes
                sn_int=element(se_int)%node(kn_int)
                ! Fluid phase variables
                ! Switch depending on the boundary condition of the node
                select case (node(sn_int)%ctype(1,1))
                  ! Perfectly permeable
                  case (0)
                    ! tau^(2): tau^(2)=phi^(2)/phi^(1)*tau^(1)
                    col=node(sn_int)%col(0,1)
                    A_c(row,col)=A_c(row,col)+hp(kn_int,il,0)*d_phi2_phi1
                    ! Un^(2): Un^(2)=-phi^(1)/phi^(2)*Un^(1)-(1-phi^(1)/phi^(2))*u_k(1)·n_k^(1)
                    ! Un^(1)
                    col=node(sn_int)%col(problem%n+1,1)
                    A_c(row,col)=A_c(row,col)+gp(kn_int,il,0)*d_phi1_phi2
                    ! u^(1)
                    do ik=1,problem%n
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,0)*p_1_phi1_phi2*element(se_int)%n_fn(ik,kn_int)
                    end do
                  ! Perfectly impermeable
                  case (1)
                    ! tau^(2) unknown
                    col=node(sn_int)%col(0,2)
                    A_c(row,col)=A_c(row,col)+hp(kn_int,il,0)
                    ! Un^(2): Un^(2)=-u^(1)·n^(1)
                    do ik=1,problem%n
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,0)*element(se_int)%n_fn(ik,kn_int)
                    end do
                  ! Partially permeable
                  case (2)
                    ! KT = k*phi^(1)*phi^(2)*i*omega
                    KT=node(sn_int)%cvalue_c(1,1,1)*phi1*phi2*c_im*omega
                    !
                    ! tau^(2): tau^(2)=phi^(2)/phi^(1)*tau^(1)+k*phi^(1)*phi^(2)*i*omega*(Un^(1)-u_k^(1)n_k^(1))
                    !
                    ! tau^(1)
                    col=node(sn_int)%col(0,1)
                    A_c(row,col)=A_c(row,col)+hp(kn_int,il,0)*d_phi2_phi1
                    ! Un^(1)
                    col=node(sn_int)%col(problem%n+1,1)
                    A_c(row,col)=A_c(row,col)+hp(kn_int,il,0)*KT
                    ! u^(1)
                    do ik=1,problem%n
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)-hp(kn_int,il,0)*KT*element(se_int)%n_fn(ik,kn_int)
                    end do
                    !
                    ! Un^(2): Un^(2)=-phi^(1)/phi^(2)*Un^(1)-(1-phi^(1)/phi^(2))*u_k(1)·n_k^(1)
                    !
                    ! Un^(1)
                    col=node(sn_int)%col(problem%n+1,1)
                    A_c(row,col)=A_c(row,col)+gp(kn_int,il,0)*d_phi1_phi2
                    ! u^(1)
                    do ik=1,problem%n
                      col=node(sn_int)%col(ik,1)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,0)*p_1_phi1_phi2*element(se_int)%n_fn(ik,kn_int)
                    end do
                end select
              end do
              ! Solid skeleton variables
              do ik=1,problem%n
                do kn_int=1,se_int_n_nodes
                  sn_int=element(se_int)%node(kn_int)
                  ! Switch depending on the boundary condition of the node
                  select case (node(sn_int)%ctype(1,1))
                    ! Perfectly permeable
                    case (0)
                      !
                      ! u_k^(2): u_k^(2)=u_k^(1)
                      !
                      col=node(sn_int)%col(            ik,1)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                      !
                      ! t_k^(2): t_k^(2)=-t_k^(1)-(1-phi^(2)/phi^(1))*n_k^(1)*tau^(1)
                      !
                      ! t_k^(1)
                      col=node(sn_int)%col(problem%n+1+ik,1)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)
                      ! tau^(1)
                      col=node(sn_int)%col(0,1)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*p_1_phi2_phi1*element(se_int)%n_fn(ik,kn_int)
                    ! Perfectly impermeable
                    case (1)
                      !
                      ! u_k^(2): u_k^(2)=u_k^(1)
                      !
                      col=node(sn_int)%col(            ik,1)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                      !
                      ! t_k^(2):t_k^(2)=-t_k^(1)-tau^(1)n_k^(1)+tau^(2)n_k^(1)
                      !
                      ! t_k^(1)
                      col=node(sn_int)%col(problem%n+1+ik,1)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)
                      ! tau^(1)
                      col=node(sn_int)%col(0,1)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                      ! tau^(2)
                      col=node(sn_int)%col(0,2)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)*element(se_int)%n_fn(ik,kn_int)
                    ! Partially permeable
                    case (2)
                      ! KT = k*phi^(1)*phi^(2)*i*omega
                      KT=node(sn_int)%cvalue_c(1,1,1)*phi1*phi2*c_im*omega
                      !
                      ! u_k^(2): u_k^(2)=u_k^(1)
                      !
                      col=node(sn_int)%col(            ik,1)
                      A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                      !
                      ! t_k^(2): t_k^(2)=-t_k^(1)+[k*phi^(1)*phi^(2)*i*omega*(Un^(1)-u_k^(1)n_k^(1))-(1-phi^(2)/phi^(1))*tau^(1)]*n_k^(1)
                      !
                      ! t_k^(1)
                      col=node(sn_int)%col(problem%n+1+ik,1)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)
                      ! Un^(1)
                      col=node(sn_int)%col(problem%n+1   ,1)
                      A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)*KT*element(se_int)%n_fn(ik,kn_int)
                      ! u^(1)
                      do ikc=1,problem%n
                        col=node(sn_int)%col(           ikc,1)
                        A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*KT*element(se_int)%n_fn(ikc,kn_int)*element(se_int)%n_fn(ik,kn_int)
                      end do
                      ! tau^(1)
                      col=node(sn_int)%col(             0,1)
                      A_c(row,col)=A_c(row,col)+gp(kn_int,il,ik)*p_1_phi2_phi1*element(se_int)%n_fn(ik,kn_int)
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
          do il=0,problem%n
            row=node(sn_col)%row(il,eq_index)
            ! Fluid phase
            do kn_int=1,se_int_n_nodes
              sn_int=element(se_int)%node(kn_int)
              sn_int_fe=element(se_int)%element_node(kn_int)
              ! tau is unknown
              col=node(sn_int)%col(0,1)
              A_c(row,col)=A_c(row,col)+hp(kn_int,il,0)
              ! Un=u·n
              if (sb_int_reversion.eqv.(.false.)) then
                do ik=1,problem%n
                  col=node(sn_int_fe)%col(ik,1)
                  if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                    A_c(row,col)=A_c(row,col)-gp(kn_int,il,0)*element(se_int)%n_fn(ik,kn_int)
                  else
                    b_c(row,1)=b_c(row,1)+gp(kn_int,il,0)*element(se_int)%n_fn(ik,kn_int)*node(sn_int_fe)%cvalue_c(ik,1,1)
                  end if
                end do
              ! Un=-u·n
              else
                do ik=1,problem%n
                  col=node(sn_int_fe)%col(ik,1)
                  if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                    A_c(row,col)=A_c(row,col)+gp(kn_int,il,0)*element(se_int)%n_fn(ik,kn_int)
                  else
                    b_c(row,1)=b_c(row,1)-gp(kn_int,il,0)*element(se_int)%n_fn(ik,kn_int)*node(sn_int_fe)%cvalue_c(ik,1,1)
                  end if
                end do
              end if
            end do
            ! Solid skeleton
            do ik=1,problem%n
              do kn_int=1,se_int_n_nodes
                sn_int=element(se_int)%node(kn_int)
                sn_int_fe=element(se_int)%element_node(kn_int)
                ! t_k is unknown
                col=node(sn_int)%col(problem%n+1+ik,1)
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
          do il=0,problem%n
            row=node(sn_col)%row(il,eq_index)
            ! Fluid phase
            do kn_int=1,se_int_n_nodes
              sn_int=element(se_int)%node(kn_int)
              sn_int_fe=element(se_int)%element_node(kn_int)
              ! Face +
              ! tau^+ is unknown
              col=node(sn_int)%col(0,1)
              A_c(row,col)=A_c(row,col)+hp(kn_int,il,0)
              ! Un^+=u^+·n^+
              if (sb_int_reversion.eqv.(.false.)) then
                do ik=1,problem%n
                  col=node(sn_int_fe)%col(ik,1)
                  if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                    A_c(row,col)=A_c(row,col)-gp(kn_int,il,0)*element(se_int)%n_fn(ik,kn_int)
                  else
                    b_c(row,1)=b_c(row,1)+gp(kn_int,il,0)*element(se_int)%n_fn(ik,kn_int)*node(sn_int_fe)%cvalue_c(ik,1,1)
                  end if
                end do
              ! Un^+=-u^+·n^+
              else
                do ik=1,problem%n
                  col=node(sn_int_fe)%col(ik,1)
                  if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                    A_c(row,col)=A_c(row,col)+gp(kn_int,il,0)*element(se_int)%n_fn(ik,kn_int)
                  else
                    b_c(row,1)=b_c(row,1)-gp(kn_int,il,0)*element(se_int)%n_fn(ik,kn_int)*node(sn_int_fe)%cvalue_c(ik,1,1)
                  end if
                end do
              end if
              ! Face -
              ! tau^- is unknown
              col=node(sn_int)%col(0,2)
              A_c(row,col)=A_c(row,col)+hm(kn_int,il,0)
              ! Un^-=-u^-·n^+
              if (sb_int_reversion.eqv.(.false.)) then
                do ik=1,problem%n
                  col=node(sn_int_fe)%col(ik,1)
                  if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                    A_c(row,col)=A_c(row,col)+gm(kn_int,il,0)*element(se_int)%n_fn(ik,kn_int)
                  else
                    b_c(row,1)=b_c(row,1)-gm(kn_int,il,0)*element(se_int)%n_fn(ik,kn_int)*node(sn_int_fe)%cvalue_c(ik,1,1)
                  end if
                end do
              ! Un^-=u^-·n^+
              else
                do ik=1,problem%n
                  col=node(sn_int_fe)%col(ik,1)
                  if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                    A_c(row,col)=A_c(row,col)-gm(kn_int,il,0)*element(se_int)%n_fn(ik,kn_int)
                  else
                    b_c(row,1)=b_c(row,1)+gm(kn_int,il,0)*element(se_int)%n_fn(ik,kn_int)*node(sn_int_fe)%cvalue_c(ik,1,1)
                  end if
                end do
              end if
            end do
            ! Solid skeleton
            do ik=1,problem%n
              do kn_int=1,se_int_n_nodes
                sn_int=element(se_int)%node(kn_int)
                sn_int_fe=element(se_int)%element_node(kn_int)
                ! Face +
                ! t_k^+ is unknown
                col=node(sn_int)%col(problem%n+1+ik,1)
                A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
                ! If u_k is unknown
                if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                  col=node(sn_int_fe)%col(ik,1)
                  A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
                ! If u_k is known
                else
                  b_c(row,1)=b_c(row,1)-hp(kn_int,il,ik)*node(sn_int_fe)%cvalue_c(ik,1,1)
                end if
                ! Face -
                ! t_k^- is unknown
                col=node(sn_int)%col(problem%n+1+ik,2)
                A_c(row,col)=A_c(row,col)-gm(kn_int,il,ik)
                ! If u_k is unknown
                if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                  col=node(sn_int_fe)%col(ik,1)
                  A_c(row,col)=A_c(row,col)+hm(kn_int,il,ik)
                ! If u_k is known
                else
                  b_c(row,1)=b_c(row,1)-hm(kn_int,il,ik)*node(sn_int_fe)%cvalue_c(ik,1,1)
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

        do il=0,problem%n
          row=node(sn_col)%row(il,eq_index)
          ! Fluid phase
          do kn_int=1,se_int_n_nodes
            sn_int=element(se_int)%node(kn_int)
            sn_int_fe=element(se_int)%element_node(kn_int)
            ! tau^(1) is unknown
            col=node(sn_int)%col(0,1)
            A_c(row,col)=A_c(row,col)+hp(kn_int,il,0)
            ! Un^(1)=u^(FE)·n^(1)
            do ik=1,problem%n
              col=node(sn_int_fe)%col(ik,1)
              ! u_k^(FE) is unknown
              if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                A_c(row,col)=A_c(row,col)-gp(kn_int,il,0)*element(se_int)%n_fn(ik,kn_int)
              ! u_k^(FE) is known
              else
                b_c(row,1)=b_c(row,1)+gp(kn_int,il,0)*element(se_int)%n_fn(ik,kn_int)*node(sn_int_fe)%cvalue_c(ik,1,1)
              end if
            end do
          end do
          ! Solid skeleton
          do ik=1,problem%n
            do kn_int=1,se_int_n_nodes
              sn_int=element(se_int)%node(kn_int)
              sn_int_fe=element(se_int)%element_node(kn_int)
              ! t_k^(1) is unknown
              col=node(sn_int)%col(problem%n+1+ik,1)
              A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
              ! u_k^(2) = u_k^(FE)
              ! u_k^(FE) is unknown
              if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                col=node(sn_int_fe)%col(ik,1)
                A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
              ! u_k^(FE) is known
              else
                b_c(row,1)=b_c(row,1)-hp(kn_int,il,ik)*node(sn_int_fe)%cvalue_c(ik,1,1)
              end if
            end do
          end do
        end do

      else

        do il=0,problem%n
          row=node(sn_col)%row(il,eq_index)
          ! Fluid phase
          do kn_int=1,se_int_n_nodes
            sn_int=element(se_int)%node(kn_int)
            sn_int_fe=element(se_int)%element_node(kn_int)
            ! tau^(2) is unknown
            col=node(sn_int)%col(0,2)
            A_c(row,col)=A_c(row,col)+hp(kn_int,il,0)
            ! Un^(2)=-u^(FE)·n^(1)
            do ik=1,problem%n
              col=node(sn_int_fe)%col(ik,1)
              ! u_k^(FE) is unknown
              if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                A_c(row,col)=A_c(row,col)+gp(kn_int,il,0)*element(se_int)%n_fn(ik,kn_int)
              ! u_k^(FE) is known
              else
                b_c(row,1)=b_c(row,1)-gp(kn_int,il,0)*element(se_int)%n_fn(ik,kn_int)*node(sn_int_fe)%cvalue_c(ik,1,1)
              end if
            end do
          end do
          ! Solid skeleton
          do ik=1,problem%n
            do kn_int=1,se_int_n_nodes
              sn_int=element(se_int)%node(kn_int)
              sn_int_fe=element(se_int)%element_node(kn_int)
              ! t_k^(2) is unknown
              col=node(sn_int)%col(problem%n+1+ik,2)
              A_c(row,col)=A_c(row,col)-gp(kn_int,il,ik)
              ! u_k^(2) = u_k^(FE)
              ! u_k^(FE) is unknown
              if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                col=node(sn_int_fe)%col(ik,1)
                A_c(row,col)=A_c(row,col)+hp(kn_int,il,ik)
              ! u_k^(FE) is known
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
        do il=0,problem%n
          row=node(sn_col)%row(il,eq_index)
          do ik=0,problem%n
            do kn_int=1,se_int_n_nodes
              b_c(row,1)=b_c(row,1)+hp(kn_int,il,ik)*u_inc(ik,kn_int)-gp(kn_int,il,ik)*t_inc(ik,kn_int)
            end do
          end do
        end do

      ! -------------------
      ! CRACK-LIKE BOUNDARY
      ! -------------------

      case (fbem_boundary_class_cracklike)
        do il=0,problem%n
          row=node(sn_col)%row(il,eq_index)
          do ik=0,problem%n
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

end subroutine assemble_bem_harpor_equation
