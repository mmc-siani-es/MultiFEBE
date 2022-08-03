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

subroutine assemble_bem_laplace_equation_sa(kr,sb_int,sb_int_reversion,&
                                           se_int,se_int_n_nodes,sdme_int_n_nodes,&
                                           h1p,h2p,h3p,g1p,g2p,g3p,h1m,h2m,h3m,g1m,g2m,g3m,dxda,&
                                           dxda_i,dnda_i,sn_col,eq_index)

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
  integer           :: kr
  integer           :: sb_int
  logical           :: sb_int_reversion
  integer           :: se_int
  integer           :: se_int_n_nodes
  integer           :: sdme_int_n_nodes
  real(kind=real64) :: h1p(sdme_int_n_nodes,se_int_n_nodes,problem%n)
  real(kind=real64) :: h2p(                 se_int_n_nodes,problem%n)
  real(kind=real64) :: h3p(                 se_int_n_nodes,problem%n)
  real(kind=real64) :: g1p(sdme_int_n_nodes,se_int_n_nodes,problem%n)
  real(kind=real64) :: g2p(                 se_int_n_nodes,problem%n)
  real(kind=real64) :: g3p(                 se_int_n_nodes,problem%n)
  real(kind=real64) :: h1m(sdme_int_n_nodes,se_int_n_nodes,problem%n)
  real(kind=real64) :: h2m(                 se_int_n_nodes,problem%n)
  real(kind=real64) :: h3m(                 se_int_n_nodes,problem%n)
  real(kind=real64) :: g1m(sdme_int_n_nodes,se_int_n_nodes,problem%n)
  real(kind=real64) :: g2m(                 se_int_n_nodes,problem%n)
  real(kind=real64) :: g3m(                 se_int_n_nodes,problem%n)
  real(kind=real64) :: dxda(problem%n,sdme_int_n_nodes,problem%n_designvariables)
  real(kind=real64) :: dxda_i(problem%n,problem%n_designvariables)
  real(kind=real64) :: dnda_i(problem%n,problem%n_designvariables)
  integer           :: sn_col
  integer           :: eq_index

  ! Local variables
  integer           :: ka
  integer           :: il, ik, im, ip, iq
  integer           :: kn_int, sn_int
  integer           :: tmp_sn_col
  integer           :: row, col
  real(kind=real64) :: hp(se_int_n_nodes)
  real(kind=real64) :: gp(se_int_n_nodes)
  real(kind=real64) :: hm(se_int_n_nodes)
  real(kind=real64) :: gm(se_int_n_nodes)

  ! Loop through the design variables
  do ka=1,problem%n_designvariables

    ! Assemble h and g taking into account the design variables velocities
    hp=0.d0
    gp=0.d0
    hm=0.d0
    gm=0.d0
    do im=1,problem%n
      do iq=1,sdme_int_n_nodes
        hp(:)=hp(:)+h1p(iq,:,im)*dxda(im,iq,ka)
        gp(:)=gp(:)+g1p(iq,:,im)*dxda(im,iq,ka)
        hm(:)=hm(:)+h1m(iq,:,im)*dxda(im,iq,ka)
        gm(:)=gm(:)+g1m(iq,:,im)*dxda(im,iq,ka)
      end do
      hp(:)=hp(:)-h2p(:,im)*dxda_i(im,ka)
      gp(:)=gp(:)-g2p(:,im)*dxda_i(im,ka)
      hm(:)=hm(:)-h2m(:,im)*dxda_i(im,ka)
      gm(:)=gm(:)-g2m(:,im)*dxda_i(im,ka)
      hp(:)=hp(:)+h3p(:,im)*dnda_i(im,ka)
      gp(:)=gp(:)+g3p(:,im)*dnda_i(im,ka)
      hm(:)=hm(:)+h3m(:,im)*dnda_i(im,ka)
      gm(:)=gm(:)+g3m(:,im)*dnda_i(im,ka)
    end do

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
              bsa_r(row,ka)=bsa_r(row,ka)-hp(kn_int)*node(sn_int)%value_r(1,1)
              bsa_r(row,ka)=bsa_r(row,ka)+gp(kn_int)*node(sn_int)%value_r(2,1)
            end do

          ! -------------------
          ! CRACK-LIKE BOUNDARY
          ! -------------------

          case (fbem_boundary_class_cracklike)
            row=node(sn_col)%row(1,eq_index)
            do kn_int=1,se_int_n_nodes
              sn_int=element(se_int)%node(kn_int)
              bsa_r(row,ka)=bsa_r(row,ka)-hp(kn_int)*node(sn_int)%value_r(1,1)
              bsa_r(row,ka)=bsa_r(row,ka)+gp(kn_int)*node(sn_int)%value_r(2,1)
              bsa_r(row,ka)=bsa_r(row,ka)-hm(kn_int)*node(sn_int)%value_r(1,2)
              bsa_r(row,ka)=bsa_r(row,ka)+gm(kn_int)*node(sn_int)%value_r(2,2)
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
            bsa_r(row,ka)=bsa_r(row,ka)-hp(kn_int)*node(sn_int)%value_r(1,1)
            bsa_r(row,ka)=bsa_r(row,ka)+gp(kn_int)*node(sn_int)%value_r(2,1)
          end do
        else
          row=node(sn_col)%row(1,eq_index)
          do kn_int=1,se_int_n_nodes
            sn_int=element(se_int)%node(kn_int)
            bsa_r(row,ka)=bsa_r(row,ka)-hp(kn_int)*node(sn_int)%value_r(1,2)
            bsa_r(row,ka)=bsa_r(row,ka)+gp(kn_int)*node(sn_int)%value_r(2,2)
          end do
        end if

      ! ==============
      ! BE-FE BOUNDARY
      ! ==============

      case (fbem_boundary_coupling_be_fe)
        stop 'not yet 15'

      ! =================
      ! BE-FE-BE BOUNDARY
      ! =================

      case (fbem_boundary_coupling_be_fe_be)
       stop 'not yet 16'

    end select

  end do

end subroutine assemble_bem_laplace_equation_sa
