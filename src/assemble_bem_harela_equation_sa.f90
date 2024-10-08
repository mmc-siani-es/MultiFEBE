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

subroutine assemble_bem_harela_equation_sa(omega,kr,sb_int,sb_int_reversion,&
                                           se_int,se_int_n_nodes,sdme_int_n_nodes,h1,h2,g1,g2,dxda,dxda_i,&
                                           sn_col,eq_index)

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
  integer              :: sdme_int_n_nodes
  complex(kind=real64) :: h1(sdme_int_n_nodes,se_int_n_nodes,problem%n,problem%n,problem%n)
  complex(kind=real64) :: h2(                 se_int_n_nodes,problem%n,problem%n,problem%n)
  complex(kind=real64) :: g1(sdme_int_n_nodes,se_int_n_nodes,problem%n,problem%n,problem%n)
  complex(kind=real64) :: g2(                 se_int_n_nodes,problem%n,problem%n,problem%n)
  real(kind=real64)    :: dxda(problem%n,sdme_int_n_nodes,problem%n_designvariables)
  real(kind=real64)    :: dxda_i(problem%n,problem%n_designvariables)
  integer              :: sn_col
  integer              :: eq_index

  ! Local variables
  integer              :: ka
  integer              :: il, ik, im, ip, iq
  integer              :: kn_int, sn_int
  integer              :: tmp_sn_col
  integer              :: row, col
  complex(kind=real64) :: h(se_int_n_nodes,problem%n,problem%n)
  complex(kind=real64) :: g(se_int_n_nodes,problem%n,problem%n)

  ! Loop through the design variables
  do ka=1,problem%n_designvariables

    ! Assemble h and g taking into account the design variables velocities
    h=(0.d0,0d0)
    g=(0.d0,0d0)
    do im=1,problem%n
      do iq=1,sdme_int_n_nodes
        h(:,:,:)=h(:,:,:)+h1(iq,:,:,:,im)*dxda(im,iq,ka)
        g(:,:,:)=g(:,:,:)+g1(iq,:,:,:,im)*dxda(im,iq,ka)
      end do
      h(:,:,:)=h(:,:,:)-h2(:,:,:,im)*dxda_i(im,ka)
      g(:,:,:)=g(:,:,:)-g2(:,:,:,im)*dxda_i(im,ka)
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
            do il=1,problem%n
              row=node(sn_col)%row(il,eq_index)
              do ik=1,problem%n
                do kn_int=1,se_int_n_nodes
                  sn_int=element(se_int)%node(kn_int)
                  bsa_c(row,ka)=bsa_c(row,ka)-h(kn_int,il,ik)*node(sn_int)%value_c(          ik,1)
                  bsa_c(row,ka)=bsa_c(row,ka)+g(kn_int,il,ik)*node(sn_int)%value_c(problem%n+ik,1)
                end do
              end do
            end do

          ! -------------------
          ! CRACK-LIKE BOUNDARY
          ! -------------------

          case (fbem_boundary_class_cracklike)
            stop 'not yet 1'

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
              stop 'not yet 2'

            ! -----------------------------------------------
            ! VISCOELASTIC SOLID (1) - VISCOELASTIC SOLID (2)
            ! -----------------------------------------------

            case (fbem_viscoelastic)
              do il=1,problem%n
                row=node(sn_col)%row(il,eq_index)
                do ik=1,problem%n
                  do kn_int=1,se_int_n_nodes
                    sn_int=element(se_int)%node(kn_int)
                    bsa_c(row,ka)=bsa_c(row,ka)-h(kn_int,il,ik)*node(sn_int)%value_c(          ik,1)
                    bsa_c(row,ka)=bsa_c(row,ka)+g(kn_int,il,ik)*node(sn_int)%value_c(problem%n+ik,1)
                  end do
                end do
              end do

            ! ----------------------------------------------
            ! VISCOELASTIC SOLID (1) - POROELASTIC MEDIA (2)
            ! ----------------------------------------------

            case (fbem_poroelastic)
              stop 'not yet 3'

          end select
        else
          select case (region(boundary(sb_int)%region(1))%type)

            ! -------------------------------------------
            ! INVISCID FLUID (1) - VISCOELASTIC SOLID (2)
            ! -------------------------------------------

            case (fbem_potential)
              stop 'not yet 4'

            ! -----------------------------------------------
            ! VISCOELASTIC SOLID (1) - VISCOELASTIC SOLID (2)
            ! -----------------------------------------------

            case (fbem_viscoelastic)
              do il=1,problem%n
                row=node(sn_col)%row(il,eq_index)
                do ik=1,problem%n
                  do kn_int=1,se_int_n_nodes
                    sn_int=element(se_int)%node(kn_int)
                    bsa_c(row,ka)=bsa_c(row,ka)-h(kn_int,il,ik)*node(sn_int)%value_c(          ik,2)
                    bsa_c(row,ka)=bsa_c(row,ka)+g(kn_int,il,ik)*node(sn_int)%value_c(problem%n+ik,2)
                  end do
                end do
              end do

            ! ----------------------------------------------
            ! POROELASTIC MEDIA (1) - VISCOELASTIC SOLID (2)
            ! ----------------------------------------------

            case (fbem_poroelastic)
              stop 'not yet 5'

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
                  bsa_c(row,ka)=bsa_c(row,ka)-h(kn_int,il,ik)*node(sn_int)%value_c(          ik,1)
                  bsa_c(row,ka)=bsa_c(row,ka)+g(kn_int,il,ik)*node(sn_int)%value_c(problem%n+ik,1)
                end do
              end do
            end do

          ! -------------------
          ! CRACK-LIKE BOUNDARY
          ! -------------------

          case (fbem_boundary_class_cracklike)
            stop 'not yet 6'

        end select

      ! =================
      ! BE-FE-BE BOUNDARY
      ! =================

      case (fbem_boundary_coupling_be_fe_be)
       stop 'not yet 7'

    end select

  end do

end subroutine assemble_bem_harela_equation_sa
