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

subroutine assign_default_conditions_bem_boundaries_mechanics_harmonic

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! No implicit variables are allowed
  implicit none

  ! Local variables
  integer :: i

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START assigning default conditions (BEM boundaries)')

  ! Loop through BOUNDARIES
  do i=1,n_boundaries
    select case (boundary(i)%coupling)

      ! ===========
      ! BE BOUNDARY
      ! ===========

      case (fbem_boundary_coupling_be)
          select case (boundary(i)%class)

            ! =================
            ! ORDINARY BOUNDARY
            ! =================

            case (fbem_boundary_class_ordinary)
              select case (region(boundary(i)%region(1))%type)

                ! --------------
                ! INVISCID FLUID
                ! --------------

                case (fbem_potential)
                  ! Allocate structure members
                  allocate (boundary(i)%ctype(1,1))
                  allocate (boundary(i)%cvalue_c(1,1,1))
                  ! Un=0
                  boundary(i)%ctype(1,1)=1
                  boundary(i)%cvalue_c(1,1,1)=0.d0

                ! ------------------
                ! VISCOELASTIC SOLID
                ! ------------------

                case (fbem_viscoelastic)
                  ! Allocate structure members
                  allocate (boundary(i)%ctype(problem%n,1))
                  allocate (boundary(i)%cvalue_c(problem%n,7,1))
                  ! t_k=0
                  boundary(i)%ctype(:,1)=1
                  boundary(i)%cvalue_c(:,1,1)=0.d0

                ! -----------------
                ! POROELASTIC MEDIUM
                ! -----------------

                case (fbem_poroelastic)
                  ! Allocate structure members
                  allocate (boundary(i)%ctype(0:problem%n,1))
                  allocate (boundary(i)%cvalue_c(0:problem%n,7,1))
                  ! tau=0
                  boundary(i)%ctype(0,1)=0
                  boundary(i)%cvalue_c(0,1,1)=0.d0
                  ! t_k=0
                  boundary(i)%ctype(1:problem%n,1)=1
                  boundary(i)%cvalue_c(1:problem%n,1,1)=0.d0

              end select

            ! ===================
            ! CRACK-LIKE BOUNDARY
            ! ===================

            case (fbem_boundary_class_cracklike)
              select case (region(boundary(i)%region(1))%type)

                ! --------------
                ! INVISCID FLUID
                ! --------------

                case (fbem_potential)
                  ! Allocate structure members
                  allocate (boundary(i)%ctype(1,2))
                  allocate (boundary(i)%cvalue_c(1,1,2))
                  ! Face + : Un=0
                  boundary(i)%ctype(1,1)=1
                  boundary(i)%cvalue_c(1,1,1)=0.d0
                  ! Face - : Un=0
                  boundary(i)%ctype(1,2)=1
                  boundary(i)%cvalue_c(1,1,2)=0.d0

                ! ------------------
                ! VISCOELASTIC SOLID
                ! ------------------

                case (fbem_viscoelastic)
                  ! Allocate structure members
                  allocate (boundary(i)%ctype(problem%n,2))
                  allocate (boundary(i)%cvalue_c(problem%n,7,2))
                  ! Face + : t_k=0
                  boundary(i)%ctype(:,1)=1
                  boundary(i)%cvalue_c(:,1,1)=0.d0
                  ! Face - : t_k=0
                  boundary(i)%ctype(:,2)=1
                  boundary(i)%cvalue_c(:,1,2)=0.d0

                ! -----------------
                ! POROELASTIC MEDIUM
                ! -----------------

                case (fbem_poroelastic)
                  ! Allocate structure members
                  allocate (boundary(i)%ctype(0:problem%n,2))
                  allocate (boundary(i)%cvalue_c(0:problem%n,7,2))
                  ! Face + : tau=0, t_k=0
                  boundary(i)%ctype(0,1)=0
                  boundary(i)%cvalue_c(0,1,1)=0.d0
                  boundary(i)%ctype(1:problem%n,1)=1
                  boundary(i)%cvalue_c(1:problem%n,1,1)=0.d0
                  ! Face - : tau=0, t_k=0
                  boundary(i)%ctype(0,2)=0
                  boundary(i)%cvalue_c(0,1,2)=0.d0
                  boundary(i)%ctype(1:problem%n,2)=1
                  boundary(i)%cvalue_c(1:problem%n,1,2)=0.d0

              end select
            end select

      ! ==============
      ! BE-BE BOUNDARY
      ! ==============

      case (fbem_boundary_coupling_be_be)
        select case (region(boundary(i)%region(1))%type)
          case (fbem_potential)
            select case (region(boundary(i)%region(2))%type)

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
                ! Allocate structure members
                allocate (boundary(i)%ctype(1,1))
                ! p^(1)=-tau^(2)/phi^(2), -(1-phi^(2))p^(1)n_k^(1)+t_k^(2)=0, Un^(1)=-phi^(2)Un^(2)-(1-phi^(2))u_k^(2)n_k(2)
                boundary(i)%ctype(1,1)=0

            end select

          case (fbem_viscoelastic)
            select case (region(boundary(i)%region(2))%type)

              ! -------------------------------------------
              ! VISCOELASTIC SOLID (1) - INVISCID FLUID (2)
              ! -------------------------------------------

              case (fbem_potential)
                ! It does not need B.C.

              ! -----------------------------------------------
              ! VISCOELASTIC SOLID (1) - VISCOELASTIC SOLID (2)
              ! -----------------------------------------------

              case (fbem_viscoelastic)
                ! Allocate structure members
                allocate (boundary(i)%ctype(problem%n,1))
                allocate (boundary(i)%cvalue_c(problem%n,3,1))
                ! u_k^(1)=u_k^(2), t_k^(1)+t_k^(2)=0
                boundary(i)%ctype(:,1)=0
                boundary(i)%cvalue_c=(0.d0,0.d0)

              ! ----------------------------------------------
              ! VISCOELASTIC SOLID (1) - POROELASTIC MEDIA (2)
              ! ----------------------------------------------

              case (fbem_poroelastic)
                ! Allocate structure members
                allocate (boundary(i)%ctype(problem%n,1))
                allocate (boundary(i)%cvalue_c(problem%n,3,1))
                ! u_k^(1)=u_k^(2), t_k^(1)+tau^(2)n_k^(2)+t_k^(2)=0, Un^(2)=u_k^(2)n_k^(2)
                boundary(i)%ctype(:,1)=0
                boundary(i)%cvalue_c=(0.d0,0.d0)

            end select

          case (fbem_poroelastic)
            select case (region(boundary(i)%region(2))%type)

              ! ------------------------------------------
              ! POROELASTIC MEDIA (1) - INVISCID FLUID (2)
              ! ------------------------------------------

              case (fbem_potential)
                ! Allocate structure members
                allocate (boundary(i)%ctype(1,1))
                ! p^(2)=-tau^(1)/phi^(1), -(1-phi^(1))p^(2)n_k^(2)+t_k^(1)=0, Un^(2)=-phi^(1)Un^(1)-(1-phi^(1))u_k^(1)n_k(1)
                boundary(i)%ctype(1,1)=0

              ! ----------------------------------------------
              ! POROELASTIC MEDIA (1) - VISCOELASTIC SOLID (2)
              ! ----------------------------------------------

              case (fbem_viscoelastic)
                ! Allocate structure members
                allocate (boundary(i)%ctype(problem%n,1))
                allocate (boundary(i)%cvalue_c(problem%n,3,1))
                ! u_k^(1)=u_k^(2), t_k^(2)+tau^(1)n_k^(1)+t_k^(1)=0, Un^(1)=u_k^(1)n_k^(1)
                boundary(i)%ctype(:,1)=0
                boundary(i)%cvalue_c=(0.d0,0.d0)

              ! ---------------------------------------------
              ! POROELASTIC MEDIA (1) - POROELASTIC MEDIA (2)
              ! ---------------------------------------------

              case (fbem_poroelastic)
                ! Allocate structure members
                allocate (boundary(i)%ctype(1,1))
                allocate (boundary(i)%cvalue_c(1,1,1))
                ! perfectly permeable
                boundary(i)%ctype(1,1)=0
                boundary(i)%cvalue_c=0.

            end select
        end select

      ! ==============
      ! BE-FE BOUNDARY
      ! ==============

      case (fbem_boundary_coupling_be_fe)
        select case (boundary(i)%class)

          ! =================
          ! ORDINARY BOUNDARY
          ! =================

          case (fbem_boundary_class_ordinary)
            select case (region(boundary(i)%region(1))%type)

              ! --------------
              ! INVISCID FLUID
              ! --------------

              case (fbem_potential)
                ! It does not need B.C.

              ! ------------------
              ! VISCOELASTIC SOLID
              ! ------------------

              case (fbem_viscoelastic)
                ! It does not need B.C.

              ! -----------------
              ! POROELASTIC MEDIUM
              ! -----------------

              case (fbem_poroelastic)
                ! It does not need B.C.

            end select

          ! ===================
          ! CRACK-LIKE BOUNDARY
          ! ===================

          case (fbem_boundary_class_cracklike)
            select case (region(boundary(i)%region(1))%type)

              ! --------------
              ! INVISCID FLUID
              ! --------------

              case (fbem_potential)
                ! It does not need B.C.

              ! ------------------
              ! VISCOELASTIC SOLID
              ! ------------------

              case (fbem_viscoelastic)
                ! It does not need B.C.

              ! -----------------
              ! POROELASTIC MEDIUM
              ! -----------------

              case (fbem_poroelastic)
                ! It does not need B.C.

            end select
        end select

      ! =================
      ! BE-FE-BE BOUNDARY
      ! =================

      case (fbem_boundary_coupling_be_fe_be)
        select case (region(boundary(i)%region(1))%type)
          case (fbem_potential)
            select case (region(boundary(i)%region(2))%type)

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
            select case (region(boundary(i)%region(2))%type)

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
            select case (region(boundary(i)%region(2))%type)

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

  end do

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END assigning default conditions (BEM boundaries)')

end subroutine assign_default_conditions_bem_boundaries_mechanics_harmonic
