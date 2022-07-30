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

subroutine assemble_bem_harpot_equation(omega,kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,sn_col,eq_index,&
                                        hp,gp,hm,gm,p_inc,Un_inc)

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
  complex(kind=real64) :: hp (se_int_n_nodes)
  complex(kind=real64) :: gp (se_int_n_nodes)
  complex(kind=real64) :: hm(se_int_n_nodes)
  complex(kind=real64) :: gm(se_int_n_nodes)
  complex(kind=real64) :: p_inc  (se_int_n_nodes)
  complex(kind=real64) :: Un_inc (se_int_n_nodes)

  ! Local variables
  integer              :: kn_int, sn_int, ik, sn_int_fe
  integer              :: row, col
  real(kind=real64)    :: phi, rho
  complex(kind=real64) :: c

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

          row=node(sn_col)%row(1,eq_index)
          do kn_int=1,se_int_n_nodes
            sn_int=element(se_int)%node(kn_int)
            select case (node(sn_int)%ctype(1,1))
              ! p known, Un unknown
              case (0)
                col=node(sn_int)%col(2,1)
                A_c(row,col)=A_c(row,col)-gp(kn_int)
                b_c(row,1)=b_c(row,1)-hp(kn_int)*node(sn_int)%cvalue_c(1,1,1)
              ! Un known, p unknown
              case (1)
                col=node(sn_int)%col(1,1)
                A_c(row,col)=A_c(row,col)+hp(kn_int)
                b_c(row,1)=b_c(row,1)+gp(kn_int)*node(sn_int)%cvalue_c(1,1,1)
              ! p unknown, Un=-i/(rho*c*omega)p
              case (2)
                rho=region(kr)%property_r(1)
                c=region(kr)%property_c(4)
                col=node(sn_int)%col(1,1)
                A_c(row,col)=A_c(row,col)+hp(kn_int)+c_im/rho/c/omega*gp(kn_int)
              ! p unknown, Un=-(i/(rho*c*omega)+1/(2*R*rho*omega^2))p
              case (3)
                rho=region(kr)%property_r(1)
                c=region(kr)%property_c(4)
                col=node(sn_int)%col(1,1)
                A_c(row,col)=A_c(row,col)+hp(kn_int)+(c_im/rho/c/omega+1.d0/(2.d0*node(sn_int)%cvalue_c(1,1,1)*rho*omega**2))*gp(kn_int)
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
              ! p known, Un unknown
              case (0)
                col=node(sn_int)%col(2,1)
                A_c(row,col)=A_c(row,col)-gp(kn_int)
                b_c(row,1)=b_c(row,1)-hp(kn_int)*node(sn_int)%cvalue_c(1,1,1)
              ! Un known, p unknown
              case (1)
                col=node(sn_int)%col(1,1)
                A_c(row,col)=A_c(row,col)+hp(kn_int)
                b_c(row,1)=b_c(row,1)+gp(kn_int)*node(sn_int)%cvalue_c(1,1,1)
            end select
            ! Face -
            select case (node(sn_int)%ctype(1,2))
              ! p known, Un unknown
              case (0)
                col=node(sn_int)%col(2,2)
                A_c(row,col)=A_c(row,col)-gm(kn_int)
                b_c(row,1)=b_c(row,1)-hm(kn_int)*node(sn_int)%cvalue_c(1,1,2)
              ! Un known, p unknown
              case (1)
                col=node(sn_int)%col(1,2)
                A_c(row,col)=A_c(row,col)+hm(kn_int)
                b_c(row,1)=b_c(row,1)+gm(kn_int)*node(sn_int)%cvalue_c(1,1,2)
            end select
          end do

      end select

    ! ==============
    ! BE-BE BOUNDARY
    ! ==============

    case (fbem_boundary_coupling_be_be)
      if (sb_int_reversion.eqv.(.false.)) then
        select case (region(boundary(sb_int)%region(2))%type)

          ! ---------------------------------------
          ! INVISCID FLUID (1) - INVISCID FLUID (2)
          ! ---------------------------------------

          case (fbem_potential)

            row=node(sn_col)%row(1,eq_index)
            do kn_int=1,se_int_n_nodes
              sn_int=element(se_int)%node(kn_int)
              ! p for region 1
              col=node(sn_int)%col(1,1)
              A_c(row,col)=A_c(row,col)+hp(kn_int)
              ! Un for region 1
              col=node(sn_int)%col(2,1)
              A_c(row,col)=A_c(row,col)-gp(kn_int)
            end do

          ! -------------------------------------------
          ! INVISCID FLUID (1) - VISCOELASTIC SOLID (2)
          ! -------------------------------------------

          case (fbem_viscoelastic)
            row=node(sn_col)%row(1,eq_index)
            do kn_int=1,se_int_n_nodes
              sn_int=element(se_int)%node(kn_int)
              ! p for region 1
              col=node(sn_int)%col(1,1)
              A_c(row,col)=A_c(row,col)+hp(kn_int)
              ! Un for region 1 is substituted by Un_{region 1}=u_k_{region 2}·n_k_{region 1}
              ! note that element(se_int)%n_fn(:,:) is the unit normal for region 1
              do ik=1,problem%n
                col=node(sn_int)%col(ik,2)
                A_c(row,col)=A_c(row,col)-gp(kn_int)*element(se_int)%n_fn(ik,kn_int)
              end do
            end do

          ! ------------------------------------------
          ! INVISCID FLUID (1) - POROELASTIC MEDIA (2)
          ! ------------------------------------------

          case (fbem_poroelastic)
            ! Porosity of region 2
            phi=region(boundary(sb_int)%region(2))%property_r(8)
            ! BIE equation
            row=node(sn_col)%row(1,eq_index)
            do kn_int=1,se_int_n_nodes
              sn_int=element(se_int)%node(kn_int)
              ! Assemble depends on the B.C.
              select case (node(sn_int)%ctype(1,1))
                ! 0. Perfectly permeable
                case (0)
                  ! p^(1)=-tau^(2)/phi^(2)
                  col=node(sn_int)%col(0,2)
                  A_c(row,col)=A_c(row,col)-hp(kn_int)/phi
                  ! Un^(1)=-phi^(2)*Un^(2)+(1-phi^(2))*u^(2)·n^(1)
                  ! Un^(2)
                  col=node(sn_int)%col(problem%n+1,2)
                  A_c(row,col)=A_c(row,col)+gp(kn_int)*phi
                  ! u_k^(2)
                  do ik=1,problem%n
                    col=node(sn_int)%col(ik,2)
                    A_c(row,col)=A_c(row,col)-gp(kn_int)*(1.0d0-phi)*element(se_int)%n_fn(ik,kn_int)
                  end do
                ! 1. Perfectly impermeable
                case (1)
                  ! p^(1)
                  col=node(sn_int)%col(1,1)
                  A_c(row,col)=A_c(row,col)+hp(kn_int)
                  ! Un^(1)=u_k^(2)*n_k^(1)
                  do ik=1,problem%n
                    col=node(sn_int)%col(ik,2)
                    A_c(row,col)=A_c(row,col)-gp(kn_int)*element(se_int)%n_fn(ik,kn_int)
                  end do
              end select
            end do

        end select
      else
        select case (region(boundary(sb_int)%region(1))%type)

          ! ---------------------------------------
          ! INVISCID FLUID (1) - INVISCID FLUID (2)
          ! ---------------------------------------

          case (fbem_potential)

            row=node(sn_col)%row(1,eq_index)
            do kn_int=1,se_int_n_nodes
              sn_int=element(se_int)%node(kn_int)
              ! p for region 2 (p_{region 1}=p_{region 2})
              col=node(sn_int)%col(1,1)
              A_c(row,col)=A_c(row,col)+hp(kn_int)
              ! w for region 2 (w_{region 1}=-w_{region 2})
              col=node(sn_int)%col(2,1)
              A_c(row,col)=A_c(row,col)+gp(kn_int)
            end do

          ! -------------------------------------------
          ! VISCOELASTIC SOLID (1) - INVISCID FLUID (2)
          ! -------------------------------------------

          case (fbem_viscoelastic)
            row=node(sn_col)%row(1,eq_index)
            do kn_int=1,se_int_n_nodes
              sn_int=element(se_int)%node(kn_int)
              ! p for region 2
              col=node(sn_int)%col(1,2)
              A_c(row,col)=A_c(row,col)+hp(kn_int)
              ! Un for region 2 is substituted by Un_{region 2}=u_k_{region 1}·n_k_{region 2}=-u_k_{region 1}·n_k_{region 1}
              ! note that element(se_int)%n_fn(:,:) is the unit normal for region 1
              do ik=1,problem%n
                col=node(sn_int)%col(ik,1)
                A_c(row,col)=A_c(row,col)+gp(kn_int)*element(se_int)%n_fn(ik,kn_int)
              end do
            end do

          ! ------------------------------------------
          ! POROELASTIC MEDIA (1) - INVISCID FLUID (2)
          ! ------------------------------------------

          case (fbem_poroelastic)
            ! Porosity of region 1
            phi=region(boundary(sb_int)%region(1))%property_r(8)
            ! BIE equation
            row=node(sn_col)%row(1,eq_index)
            do kn_int=1,se_int_n_nodes
              sn_int=element(se_int)%node(kn_int)
              ! Assemble depends on the B.C.
              select case (node(sn_int)%ctype(1,1))
                ! 0. Perfectly permeable
                case (0)
                  ! p^(2)=-tau^(1)/phi^(1)
                  col=node(sn_int)%col(0,1)
                  A_c(row,col)=A_c(row,col)-hp(kn_int)/phi
                  ! Un^(2)=-phi^(1)*Un^(1)-(1-phi^(1))*u^(1)·n^(1)
                  ! Un^(1)
                  col=node(sn_int)%col(problem%n+1,1)
                  A_c(row,col)=A_c(row,col)+gp(kn_int)*phi
                  ! u_k^(1)
                  do ik=1,problem%n
                    col=node(sn_int)%col(ik,1)
                    A_c(row,col)=A_c(row,col)+gp(kn_int)*(1.0d0-phi)*element(se_int)%n_fn(ik,kn_int)
                  end do
                ! 1. Perfectly impermeable
                case (1)
                  ! p^(2)
                  col=node(sn_int)%col(1,2)
                  A_c(row,col)=A_c(row,col)+hp(kn_int)
                  ! Un^(2)=-u_k^(1)*n_k^(1)
                  do ik=1,problem%n
                    col=node(sn_int)%col(ik,1)
                    A_c(row,col)=A_c(row,col)+gp(kn_int)*element(se_int)%n_fn(ik,kn_int)
                  end do
              end select
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
          row=node(sn_col)%row(1,eq_index)
          do kn_int=1,se_int_n_nodes
            sn_int=element(se_int)%node(kn_int)
            sn_int_fe=element(se_int)%element_node(kn_int)
            ! p is unknown
            col=node(sn_int)%col(1,1)
            A_c(row,col)=A_c(row,col)+hp(kn_int)
            ! Un=u^(fe)·n
            if (sb_int_reversion.eqv.(.false.)) then
              do ik=1,problem%n
                if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                  col=node(sn_int_fe)%col(ik,1)
                  A_c(row,col)=A_c(row,col)-gp(kn_int)*element(se_int)%n_fn(ik,kn_int)
                else
                  b_c(row,1)=b_c(row,1)+gp(kn_int)*element(se_int)%n_fn(ik,kn_int)*node(sn_int_fe)%cvalue_c(ik,1,1)
                end if
              end do
            ! Un=-u^(fe)·n
            else
              do ik=1,problem%n
                if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                  col=node(sn_int_fe)%col(ik,1)
                  A_c(row,col)=A_c(row,col)+gp(kn_int)*element(se_int)%n_fn(ik,kn_int)
                else
                  b_c(row,1)=b_c(row,1)-gp(kn_int)*element(se_int)%n_fn(ik,kn_int)*node(sn_int_fe)%cvalue_c(ik,1,1)
                end if
              end do
            end if
          end do

        ! -------------------
        ! CRACK-LIKE BOUNDARY
        ! -------------------

        case (fbem_boundary_class_cracklike)
          row=node(sn_col)%row(1,eq_index)
          do kn_int=1,se_int_n_nodes
            sn_int=element(se_int)%node(kn_int)
            sn_int_fe=element(se_int)%element_node(kn_int)
            ! Face +
            ! p^+ is unknown
            col=node(sn_int)%col(1,1)
            A_c(row,col)=A_c(row,col)+hp(kn_int)
            ! Un^+=u^(fe)·n
            if (sb_int_reversion.eqv.(.false.)) then
              do ik=1,problem%n
                col=node(sn_int_fe)%col(ik,1)
                if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                  A_c(row,col)=A_c(row,col)-gp(kn_int)*element(se_int)%n_fn(ik,kn_int)
                else
                  b_c(row,1)=b_c(row,1)+gp(kn_int)*element(se_int)%n_fn(ik,kn_int)*node(sn_int_fe)%cvalue_c(ik,1,1)
                end if
              end do
            ! Un^+=-u^(fe)·n
            else
              do ik=1,problem%n
                col=node(sn_int_fe)%col(ik,1)
                if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                  A_c(row,col)=A_c(row,col)+gp(kn_int)*element(se_int)%n_fn(ik,kn_int)
                else
                  b_c(row,1)=b_c(row,1)-gp(kn_int)*element(se_int)%n_fn(ik,kn_int)*node(sn_int_fe)%cvalue_c(ik,1,1)
                end if
              end do
            end if
            ! Face -
            ! p^- is unknown
            col=node(sn_int)%col(1,2)
            A_c(row,col)=A_c(row,col)+hm(kn_int)
            ! Un^-=-u^(fe)·n
            if (sb_int_reversion.eqv.(.false.)) then
              do ik=1,problem%n
                col=node(sn_int_fe)%col(ik,1)
                if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                  A_c(row,col)=A_c(row,col)+gm(kn_int)*element(se_int)%n_fn(ik,kn_int)
                else
                  b_c(row,1)=b_c(row,1)-gm(kn_int)*element(se_int)%n_fn(ik,kn_int)*node(sn_int_fe)%cvalue_c(ik,1,1)
                end if
              end do
            ! Un^-=u^(fe)·n
            else
              do ik=1,problem%n
                col=node(sn_int_fe)%col(ik,1)
                if (node(sn_int_fe)%ctype(ik,1).eq.1) then
                  A_c(row,col)=A_c(row,col)-gm(kn_int)*element(se_int)%n_fn(ik,kn_int)
                else
                  b_c(row,1)=b_c(row,1)+gm(kn_int)*element(se_int)%n_fn(ik,kn_int)*node(sn_int_fe)%cvalue_c(ik,1,1)
                end if
              end do
            end if
          end do
      end select

    ! =================
    ! BE-FE-BE BOUNDARY
    ! =================

    case (fbem_boundary_coupling_be_fe_be)
      if (sb_int_reversion.eqv.(.false.)) then

        row=node(sn_col)%row(1,eq_index)
        do kn_int=1,se_int_n_nodes
          sn_int=element(se_int)%node(kn_int)
          sn_int_fe=element(se_int)%element_node(kn_int)
          ! p^(1)
          col=node(sn_int)%col(1,1)
          A_c(row,col)=A_c(row,col)+hp(kn_int)
          ! Un^(1)=u^(fe)·n^(1)
          do ik=1,problem%n
            if (node(sn_int_fe)%ctype(ik,1).eq.1) then
              col=node(sn_int_fe)%col(ik,1)
              A_c(row,col)=A_c(row,col)-gp(kn_int)*element(se_int)%n_fn(ik,kn_int)
            else
              b_c(row,1)=b_c(row,1)+gp(kn_int)*element(se_int)%n_fn(ik,kn_int)*node(sn_int_fe)%cvalue_c(ik,1,1)
            end if
          end do
        end do

      else

        row=node(sn_col)%row(1,eq_index)
        do kn_int=1,se_int_n_nodes
          sn_int=element(se_int)%node(kn_int)
          sn_int_fe=element(se_int)%element_node(kn_int)
          ! p^(2)
          col=node(sn_int)%col(1,2)
          A_c(row,col)=A_c(row,col)+hp(kn_int)
          ! Un^(2)=-u^(fe)·n^(1)
          do ik=1,problem%n
            if (node(sn_int_fe)%ctype(ik,1).eq.1) then
              col=node(sn_int_fe)%col(ik,1)
              A_c(row,col)=A_c(row,col)+gp(kn_int)*element(se_int)%n_fn(ik,kn_int)
            else
              b_c(row,1)=b_c(row,1)-gp(kn_int)*element(se_int)%n_fn(ik,kn_int)*node(sn_int_fe)%cvalue_c(ik,1,1)
            end if
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
        row=node(sn_col)%row(1,eq_index)
        do kn_int=1,se_int_n_nodes
          b_c(row,1)=b_c(row,1)+hp(kn_int)*p_inc(kn_int)-gp(kn_int)*Un_inc(kn_int)
        end do

      ! -------------------
      ! CRACK-LIKE BOUNDARY
      ! -------------------

      case (fbem_boundary_class_cracklike)
        row=node(sn_col)%row(1,eq_index)
        do kn_int=1,se_int_n_nodes
          ! Face +
          b_c(row,1)=b_c(row,1)+hp(kn_int)*p_inc(kn_int)-gp(kn_int)*Un_inc(kn_int)
          ! Face -
          b_c(row,1)=b_c(row,1)+hm(kn_int)*p_inc(kn_int)+gm(kn_int)*Un_inc(kn_int)
        end do

    end select
  end if

end subroutine assemble_bem_harpot_equation
