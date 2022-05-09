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

subroutine transfer_conditions_bem_boundaries_mechanics_harmonic

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_string_handling
  use fbem_geometry
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! No implicit variables are allowed
  implicit none

  ! Local variables
  integer           :: kb, sp
  integer           :: kn, sn
  integer           :: k
  real(kind=real64) :: center(3), axis(3), x(3), urot(3)

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START transferring conditions (BEM boundaries->BEM nodes)')

  ! Loop through BOUNDARIES
  do kb=1,n_boundaries
    sp=boundary(kb)%part
    select case (boundary(kb)%coupling)

      ! ============================================================================================================================
      ! BE BOUNDARY
      ! ============================================================================================================================

      case (fbem_boundary_coupling_be)
        select case (boundary(kb)%class)

          ! =================
          ! ORDINARY BOUNDARY
          ! =================

          case (fbem_boundary_class_ordinary)
            select case (region(boundary(kb)%region(1))%type)

              ! --------------
              ! INVISCID FLUID
              ! --------------

              case (fbem_potential)
                do kn=1,part(sp)%n_nodes
                  sn=part(sp)%node(kn)
                  allocate (node(sn)%ctype(1,1))
                  allocate (node(sn)%cvalue_c(1,1,1))
                  node(sn)%ctype=boundary(kb)%ctype
                  node(sn)%cvalue_c=boundary(kb)%cvalue_c
                end do

              ! ------------------
              ! VISCOELASTIC SOLID
              ! ------------------

              case (fbem_viscoelastic)
                do kn=1,part(sp)%n_nodes
                  sn=part(sp)%node(kn)
                  allocate (node(sn)%ctype(problem%n,1))
                  allocate (node(sn)%cvalue_c(problem%n,7,1))
                  node(sn)%ctype=boundary(kb)%ctype
                  node(sn)%cvalue_c=boundary(kb)%cvalue_c
                  ! Transformation of infinitesimal rotation to displacement at each node
                  do k=1,problem%n
                    if (node(sn)%ctype(k,1).eq.4) then
                      center=dreal(node(sn)%cvalue_c(k,1:3,1))
                      axis=dreal(node(sn)%cvalue_c(k,4:6,1))
                      x=node(sn)%x
                      urot=fbem_rotation_infinitesimal_displacement(center,axis,1.d0,x)
                      node(sn)%ctype(k,1)=0
                      node(sn)%cvalue_c(k,1,1)=urot(k)*node(sn)%cvalue_c(k,7,1)
                    end if
                  end do
                end do

              ! -----------------
              ! POROELASTIC MEDIUM
              ! -----------------

              case (fbem_poroelastic)
                do kn=1,part(sp)%n_nodes
                  sn=part(sp)%node(kn)
                  allocate (node(sn)%ctype(0:problem%n,1))
                  allocate (node(sn)%cvalue_c(0:problem%n,7,1))
                  node(sn)%ctype=boundary(kb)%ctype
                  node(sn)%cvalue_c=boundary(kb)%cvalue_c
                  ! Transformation of infinitesimal rotation to displacement at each node
                  do k=1,problem%n
                    if (node(sn)%ctype(k,1).eq.8) then
                      center=dreal(node(sn)%cvalue_c(k,1:3,1))
                      axis=dreal(node(sn)%cvalue_c(k,4:6,1))
                      x=node(sn)%x
                      urot=fbem_rotation_infinitesimal_displacement(center,axis,1.d0,x)
                      node(sn)%ctype(k,1)=0
                      node(sn)%cvalue_c(k,1,1)=urot(k)*node(sn)%cvalue_c(k,7,1)
                    end if
                    if (node(sn)%ctype(k,1).eq.40) then
                      center=dreal(node(sn)%cvalue_c(k,1:3,1))
                      axis=dreal(node(sn)%cvalue_c(k,4:6,1))
                      x=node(sn)%x
                      urot=fbem_rotation_infinitesimal_displacement(center,axis,1.d0,x)
                      node(sn)%ctype(k,1)=4
                      node(sn)%cvalue_c(k,1,1)=urot(k)*node(sn)%cvalue_c(k,7,1)
                    end if
                  end do

                end do

            end select

          ! ===================
          ! CRACK-LIKE BOUNDARY
          ! ===================

          case (fbem_boundary_class_cracklike)
            select case (region(boundary(kb)%region(1))%type)

              ! --------------
              ! INVISCID FLUID
              ! --------------

              case (fbem_potential)
                do kn=1,part(sp)%n_nodes
                  sn=part(sp)%node(kn)
                  allocate (node(sn)%ctype(1,2))
                  allocate (node(sn)%cvalue_c(1,1,2))
                  node(sn)%ctype=boundary(kb)%ctype
                  node(sn)%cvalue_c=boundary(kb)%cvalue_c
                end do

              ! ------------------
              ! VISCOELASTIC SOLID
              ! ------------------

              case (fbem_viscoelastic)
                do kn=1,part(sp)%n_nodes
                  sn=part(sp)%node(kn)
                  allocate (node(sn)%ctype(problem%n,2))
                  allocate (node(sn)%cvalue_c(problem%n,1,2))
                  node(sn)%ctype=boundary(kb)%ctype
                  node(sn)%cvalue_c=boundary(kb)%cvalue_c
                end do

              ! -----------------
              ! POROELASTIC MEDIA
              ! -----------------

              case (fbem_poroelastic)
                do kn=1,part(sp)%n_nodes
                  sn=part(sp)%node(kn)
                  allocate (node(sn)%ctype(0:problem%n,2))
                  allocate (node(sn)%cvalue_c(0:problem%n,1,2))
                  node(sn)%ctype=boundary(kb)%ctype
                  node(sn)%cvalue_c=boundary(kb)%cvalue_c
                end do

            end select

        end select

      ! ============================================================================================================================
      ! BE-BE BOUNDARY
      ! ============================================================================================================================

      case (fbem_boundary_coupling_be_be)
        select case (region(boundary(kb)%region(1))%type)
          case (fbem_potential)
            select case (region(boundary(kb)%region(2))%type)

              ! ---------------------------------------
              ! INVISCID FLUID (1) - INVISCID FLUID (2)
              ! ---------------------------------------

              case (fbem_potential)
                ! It does not need B.C.

              ! -------------------------------------------
              ! INVISCID FLUID (1) - VISCOELASTIC SOLID (2)
              ! -------------------------------------------

              case (fbem_viscoelastic)
                ! It does not need B.C.

              ! ------------------------------------------
              ! INVISCID FLUID (1) - POROELASTIC MEDIA (2)
              ! ------------------------------------------

              case (fbem_poroelastic)
                do kn=1,part(sp)%n_nodes
                  sn=part(sp)%node(kn)
                  allocate (node(sn)%ctype(1,1))
                  node(sn)%ctype=boundary(kb)%ctype
                end do

            end select

          case (fbem_viscoelastic)
            select case (region(boundary(kb)%region(2))%type)

              ! -------------------------------------------
              ! VISCOELASTIC SOLID (1) - INVISCID FLUID (2)
              ! -------------------------------------------

              case (fbem_potential)
                ! It does not need B.C.

              ! -----------------------------------------------
              ! VISCOELASTIC SOLID (1) - VISCOELASTIC SOLID (2)
              ! -----------------------------------------------

              case (fbem_viscoelastic)
                do kn=1,part(sp)%n_nodes
                  sn=part(sp)%node(kn)
                  allocate (node(sn)%ctype(problem%n,1))
                  allocate (node(sn)%cvalue_c(problem%n,3,1))
                  node(sn)%ctype=boundary(kb)%ctype
                  node(sn)%cvalue_c=boundary(kb)%cvalue_c
                end do

              ! ----------------------------------------------
              ! VISCOELASTIC SOLID (1) - POROELASTIC MEDIA (2)
              ! ----------------------------------------------

              case (fbem_poroelastic)
                do kn=1,part(sp)%n_nodes
                  sn=part(sp)%node(kn)
                  allocate (node(sn)%ctype(problem%n,1))
                  allocate (node(sn)%cvalue_c(problem%n,3,1))
                  node(sn)%ctype=boundary(kb)%ctype
                  node(sn)%cvalue_c=boundary(kb)%cvalue_c
                end do

            end select

          case (fbem_poroelastic)
            select case (region(boundary(kb)%region(2))%type)

              ! ------------------------------------------
              ! POROELASTIC MEDIA (1) - INVISCID FLUID (2)
              ! ------------------------------------------

              case (fbem_potential)
                do kn=1,part(sp)%n_nodes
                  sn=part(sp)%node(kn)
                  allocate (node(sn)%ctype(1,1))
                  node(sn)%ctype=boundary(kb)%ctype
                end do

              ! ----------------------------------------------
              ! POROELASTIC MEDIA (1) - VISCOELASTIC SOLID (2)
              ! ----------------------------------------------

              case (fbem_viscoelastic)
                do kn=1,part(sp)%n_nodes
                  sn=part(sp)%node(kn)
                  allocate (node(sn)%ctype(problem%n,1))
                  allocate (node(sn)%cvalue_c(problem%n,3,1))
                  node(sn)%ctype=boundary(kb)%ctype
                  node(sn)%cvalue_c=boundary(kb)%cvalue_c
                end do

              ! ---------------------------------------------
              ! POROELASTIC MEDIA (1) - POROELASTIC MEDIA (2)
              ! ---------------------------------------------

              case (fbem_poroelastic)
                do kn=1,part(sp)%n_nodes
                  sn=part(sp)%node(kn)
                  allocate (node(sn)%ctype(1,1))
                  allocate (node(sn)%cvalue_c(1,1,1))
                  node(sn)%ctype=boundary(kb)%ctype
                  node(sn)%cvalue_c=boundary(kb)%cvalue_c
                end do

            end select
        end select

      ! ============================================================================================================================
      ! BE-FE BOUNDARY
      ! ============================================================================================================================

      case (fbem_boundary_coupling_be_fe)
        select case (boundary(kb)%class)

          ! ================= !
          ! ORDINARY BOUNDARY !
          ! ================= !

          case (fbem_boundary_class_ordinary)
            select case (region(boundary(kb)%region(1))%type)

              ! -------------- !
              ! INVISCID FLUID !
              ! -------------- !

              case (fbem_potential)
                ! It does not need B.C.

              ! ------------------ !
              ! VISCOELASTIC SOLID !
              ! ------------------ !

              case (fbem_viscoelastic)
                ! It does not need B.C.

              ! ----------------- !
              ! POROELASTIC MEDIUM !
              ! ----------------- !

              case (fbem_poroelastic)
                ! It does not need B.C.

            end select

          ! =================== !
          ! CRACK-LIKE BOUNDARY !
          ! =================== !

          case (fbem_boundary_class_cracklike)
            select case (region(boundary(kb)%region(1))%type)

              ! -------------- !
              ! INVISCID FLUID !
              ! -------------- !

              case (fbem_potential)
                ! It does not need B.C.

              ! ------------------ !
              ! VISCOELASTIC SOLID !
              ! ------------------ !

              case (fbem_viscoelastic)
                ! It does not need B.C.

              ! ----------------- !
              ! POROELASTIC MEDIUM !
              ! ----------------- !

              case (fbem_poroelastic)
                ! It does not need B.C.

            end select
        end select

      ! ============================================================================================================================
      ! BE-FE-BE BOUNDARY
      ! ============================================================================================================================

      case (fbem_boundary_coupling_be_fe_be)
        select case (region(boundary(kb)%region(1))%type)
          case (fbem_potential)
            select case (region(boundary(kb)%region(2))%type)

              ! ---------------------------------------
              ! INVISCID FLUID (1) - INVISCID FLUID (2)
              ! ---------------------------------------

              case (fbem_potential)
                ! It does not need B.C.

              ! -------------------------------------------
              ! INVISCID FLUID (1) - VISCOELASTIC SOLID (2)
              ! -------------------------------------------

              case (fbem_viscoelastic)
                ! It does not need B.C.

              ! ------------------------------------------
              ! INVISCID FLUID (1) - POROELASTIC MEDIA (2)
              ! ------------------------------------------

              case (fbem_poroelastic)
                ! It does not need B.C.

            end select

          case (fbem_viscoelastic)
            select case (region(boundary(kb)%region(2))%type)

              ! -------------------------------------------
              ! VISCOELASTIC SOLID (1) - INVISCID FLUID (2)
              ! -------------------------------------------

              case (fbem_potential)
                ! It does not need B.C.

              ! -----------------------------------------------
              ! VISCOELASTIC SOLID (1) - VISCOELASTIC SOLID (2)
              ! -----------------------------------------------

              case (fbem_viscoelastic)
                ! It does not need B.C.

              ! ----------------------------------------------
              ! VISCOELASTIC SOLID (1) - POROELASTIC MEDIA (2)
              ! ----------------------------------------------

              case (fbem_poroelastic)
                ! It does not need B.C.

            end select

          case (fbem_poroelastic)
            select case (region(boundary(kb)%region(2))%type)

              ! ------------------------------------------
              ! POROELASTIC MEDIA (1) - INVISCID FLUID (2)
              ! ------------------------------------------

              case (fbem_potential)
                ! It does not need B.C.

              ! ----------------------------------------------
              ! POROELASTIC MEDIA (1) - VISCOELASTIC SOLID (2)
              ! ----------------------------------------------

              case (fbem_viscoelastic)
                ! It does not need B.C.

              ! ---------------------------------------------
              ! POROELASTIC MEDIA (1) - POROELASTIC MEDIA (2)
              ! ---------------------------------------------

              case (fbem_poroelastic)
                ! It does not need B.C.

            end select
        end select
    end select
  end do ! Loop through BOUNDARIES

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END transferring conditions (BEM boundaries->BEM nodes)')

end subroutine transfer_conditions_bem_boundaries_mechanics_harmonic
