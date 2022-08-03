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
subroutine assign_solution_mechanics_harmonic_sa(kf)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling

  ! Module of problem variables
  use problem_variables

  ! No implicit variables are allowed
  implicit none

  ! I/O variables
  integer                    :: kf

  ! Local variables
  real(kind=real64)          :: omega
  integer                    :: kr, kb, ke, kn, ks, ss
  integer                    :: sb, se, sn
  logical                    :: sb_reversion
  integer                    :: k
  logical, allocatable       :: node_used(:)
  integer                    :: sn_fe
  real(kind=real64)          :: x_fn(problem%n), n_fn(problem%n)

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'START assigning the LSE solution to variables (sensitivity)')

  omega=frequency(kf)
  allocate (node_used(n_nodes))
  node_used=.false.

  do kr=1,n_regions
    select case (region(kr)%class)

      ! =========
      ! BE REGION
      ! =========

      case (fbem_be)
        do kb=1,region(kr)%n_boundaries
          sb=region(kr)%boundary(kb)
          sb_reversion=region(kr)%boundary_reversion(kb)
          do ke=1,part(boundary(sb)%part)%n_elements
            se=part(boundary(sb)%part)%element(ke)
            do kn=1,element(se)%n_nodes
              sn=element(se)%node(kn)
              if (.not.node_used(sn)) then
                node_used(sn)=.true.
                do k=1,problem%n
                  n_fn(k)=node(sn)%n_fn(k)
                  x_fn(k)=element(se)%x_fn(k,kn)
                end do
                !
                ! Boundary coupling
                !
                select case (boundary(sb)%coupling)

                  ! ===========
                  ! BE BOUNDARY
                  ! ===========

                  case (fbem_boundary_coupling_be)
                    select case (boundary(sb)%class)

                      ! =================
                      ! ORDINARY BOUNDARY
                      ! =================

                      case (fbem_boundary_class_ordinary)
                        select case (region(kr)%type)

                          ! --------------
                          ! INVISCID FLUID
                          ! --------------

                          case (fbem_potential)
                            ! Index of variables:
                            ! node(sn)%col(1,1): p
                            ! node(sn)%col(2,1): Un
                            select case (node(sn)%ctype(1,1))
                              case (0)
                                node(sn)%dvda_c(1,1,:)=(0.d0,0.d0)
                                node(sn)%dvda_c(2,1,:)=bsa_c(node(sn)%col(2,1),:)
                              case (1)
                                node(sn)%dvda_c(1,1,:)=bsa_c(node(sn)%col(1,1),1)
                                node(sn)%dvda_c(2,1,:)=(0.d0,0.d0)
                            end select

                          ! ------------------
                          ! VISCOELASTIC SOLID
                          ! ------------------

                          case (fbem_viscoelastic)
                            ! Index of variables for each coordinate k:
                            ! node(sn)%col(          k,1): u_k
                            ! node(sn)%col(problem%n+k,1): t_k
                            do k=1,problem%n
                              select case (node(sn)%ctype(k,1))
                                case (0)
                                  node(sn)%dvda_c(          k,1,:)=(0.d0,0.d0)
                                  node(sn)%dvda_c(problem%n+k,1,:)=bsa_c(node(sn)%col(problem%n+k,1),:)
                                case (1,10)
                                  node(sn)%dvda_c(          k,1,:)=bsa_c(node(sn)%col(k,1),:)
                                  node(sn)%dvda_c(problem%n+k,1,:)=(0.d0,0.d0)
                                case (2,3)
                                  node(sn)%dvda_c(          k,1,:)=bsa_c(node(sn)%col(k,1),:)
                                  node(sn)%dvda_c(problem%n+k,1,:)=bsa_c(node(sn)%col(problem%n+k,1),:)
                              end select
                            end do

                          ! ------------------
                          ! POROELASTIC MEDIUM
                          ! ------------------

                          case (fbem_poroelastic)
                            stop 'not yet 28'

                        end select

                      ! ===================
                      ! CRACK-LIKE BOUNDARY
                      ! ===================

                      case (fbem_boundary_class_cracklike)
                        select case (region(kr)%type)

                          ! --------------
                          ! INVISCID FLUID
                          ! --------------

                          case (fbem_potential)
                            ! Index of variables:
                            ! node(sn)%col(1,1): p for face +
                            ! node(sn)%col(2,1): Un for face +
                            ! node(sn)%col(1,2): p for face -
                            ! node(sn)%col(2,2): Un for face -
                            ! Face +
                            select case (node(sn)%ctype(1,1))
                              case (0)
                                node(sn)%dvda_c(1,1,:)=(0.d0,0.d0)
                                node(sn)%dvda_c(2,1,:)=bsa_c(node(sn)%col(2,1),:)
                              case (1)
                                node(sn)%dvda_c(1,1,:)=bsa_c(node(sn)%col(1,1),:)
                                node(sn)%dvda_c(2,1,:)=(0.d0,0.d0)
                            end select
                            ! Face -
                            select case (node(sn)%ctype(1,2))
                              case (0)
                                node(sn)%dvda_c(1,2,:)=(0.d0,0.d0)
                                node(sn)%dvda_c(2,2,:)=bsa_c(node(sn)%col(2,2),:)
                              case (1)
                                node(sn)%dvda_c(1,2,:)=bsa_c(node(sn)%col(1,2),:)
                                node(sn)%dvda_c(2,2,:)=(0.d0,0.d0)
                            end select

                          ! ------------------
                          ! VISCOELASTIC SOLID
                          ! ------------------

                          case (fbem_viscoelastic)
                            ! Index of variables for each coordinate k:
                            ! node(sn)%col(          k,1): u_k for face +
                            ! node(sn)%col(problem%n+k,1): t_k for face +
                            ! node(sn)%col(          k,2): u_k for face -
                            ! node(sn)%col(problem%n+k,2): t_k for face -
                            ! Face +
                            do k=1,problem%n
                              select case (node(sn)%ctype(k,1))
                                case (0)
                                  node(sn)%dvda_c(          k,1,:)=(0.d0,0.d0)
                                  node(sn)%dvda_c(problem%n+k,1,:)=bsa_c(node(sn)%col(problem%n+k,1),:)
                                case (1)
                                  node(sn)%dvda_c(          k,1,:)=bsa_c(node(sn)%col(k,1),:)
                                  node(sn)%dvda_c(problem%n+k,1,:)=(0.d0,0.d0)
                              end select
                            end do
                            ! Face -
                            do k=1,problem%n
                              select case (node(sn)%ctype(k,2))
                                case (0)
                                  node(sn)%dvda_c(          k,2,:)=(0.d0,0.d0)
                                  node(sn)%dvda_c(problem%n+k,2,:)=bsa_c(node(sn)%col(problem%n+k,2),:)
                                case (1)
                                  node(sn)%dvda_c(          k,2,:)=bsa_c(node(sn)%col(k,2),:)
                                  node(sn)%dvda_c(problem%n+k,2,:)=(0.d0,0.d0)
                              end select
                            end do

                          ! ------------------
                          ! POROELASTIC MEDIUM
                          ! ------------------

                          case (fbem_poroelastic)
                            stop 'not yet 29'

                        end select
                    end select

                  ! ==============
                  ! BE-BE BOUNDARY
                  ! ==============
                  ! The region 1 is the region where the boundary has N+
                  ! The region 2 is the region where the boundary has N-

                  case (fbem_boundary_coupling_be_be)
                    select case (region(boundary(sb)%region(1))%type)
                      case (fbem_potential)
                        select case (region(boundary(sb)%region(2))%type)

                          ! -----------------------------------------
                          ! BE (inviscid fluid) - BE (inviscid fluid)
                          ! -----------------------------------------

                          case (fbem_potential)
                            ! Index of variables:
                            ! node(sn)%col(1,1): p for region 1
                            ! node(sn)%col(2,1): Un for region 1
                            ! node(sn)%col(1,2): p for region 2 (not used since p_{region 1}=p_{region 2})
                            ! node(sn)%col(2,2): Un for region 2 (not used since Un_{region 1}=-Un_{region 2})
                            node(sn)%dvda_c(1,1,:)=bsa_c(node(sn)%col(1,1),:)
                            node(sn)%dvda_c(2,1,:)=bsa_c(node(sn)%col(2,1),:)
                            node(sn)%dvda_c(1,2,:)= node(sn)%dvda_c(1,1,:)
                            node(sn)%dvda_c(2,2,:)=-node(sn)%dvda_c(2,1,:)

                          ! ---------------------------------------------------
                          ! BE (1: inviscid fluid) - BE (2: viscoelastic solid)
                          ! ---------------------------------------------------

                          case (fbem_viscoelastic)
                            ! Index of variables:
                            ! node(sn)%col(          1,1): p for region 1
                            ! node(sn)%col(          2,1): Un for region 1 (not active since Un_{region 1}=u_k_{region 2}·n_k_{region 1})
                            ! node(sn)%col(          k,2): u_k for region 2
                            ! node(sn)%col(problem%n+k,2): t_k for region 2 (not active since t_k_{region 2}=-p_{region 1}·n_k_{region 2})
                            ! p of region 1
                            node(sn)%dvda_c(1,1,:)=bsa_c(node(sn)%col(1,1),:)
                            ! u_k of region 2
                            do k=1,problem%n
                              node(sn)%dvda_c(k,2,:)=bsa_c(node(sn)%col(k,2),:)
                            end do
                            ! Un of region 1
                            node(sn)%dvda_c(2,1,:)=0.0d0
                            do k=1,problem%n
                              node(sn)%dvda_c(2,1,:)=node(sn)%dvda_c(2,1,:)+node(sn)%dvda_c(k,2,:)*n_fn(k)
                            end do
                            ! t_k of region 2
                            do k=1,problem%n
                              node(sn)%dvda_c(problem%n+k,2,:)=node(sn)%dvda_c(1,1,:)*n_fn(k)
                            end do

                          ! ---------------------------------------------------
                          ! BE (1: inviscid fluid) - BE (2: poroelastic medium)
                          ! ---------------------------------------------------

                          case (fbem_poroelastic)
                            stop 'not yet 30'

                        end select

                      case (fbem_viscoelastic)
                        select case (region(boundary(sb)%region(2))%type)

                          ! ---------------------------------------------------
                          ! BE (1: viscoelastic solid) - BE (2: inviscid fluid)
                          ! ---------------------------------------------------

                          case (fbem_potential)
                            !
                            ! Index of variables:
                            ! node(sn)%col(          k,1): u_k for region 1
                            ! node(sn)%col(problem%n+k,1): t_k for region 1 (not active since t_k_{region 1}=-p_{region 2}·n_k_{region 1})
                            ! node(sn)%col(          1,2): p for region 2
                            ! node(sn)%col(          2,2): Un for region 2 (not active since Un_{region 2}=u_k_{region 1}·n_k_{region 2})
                            !
                            ! u_k of region 1
                            do k=1,problem%n
                              node(sn)%dvda_c(k,1,:)=bsa_c(node(sn)%col(k,1),:)
                            end do
                            ! p of region 2
                            node(sn)%dvda_c(1,2,:)=bsa_c(node(sn)%col(1,2),:)
                            ! t_k of region 1
                            do k=1,problem%n
                              node(sn)%dvda_c(problem%n+k,1,:)=-node(sn)%dvda_c(1,2,:)*n_fn(k)
                            end do
                            ! Un of region 2
                            node(sn)%dvda_c(2,2,:)=0.0d0
                            do k=1,problem%n
                              node(sn)%dvda_c(2,2,:)=node(sn)%dvda_c(2,2,:)-node(sn)%dvda_c(k,1,:)*n_fn(k)
                            end do

                          ! -------------------------------------------------
                          ! BE (viscoelastic solid) - BE (viscoelastic solid)
                          ! -------------------------------------------------

                          case (fbem_viscoelastic)
                            ! Index of variables for each coordinate k:
                            ! node(sn)%col(          k,1): u_k for region 1
                            ! node(sn)%col(problem%n+k,1): t_k for region 1
                            ! node(sn)%col(          k,2): u_k for region 2
                            ! node(sn)%col(problem%n+k,2): t_k for region 2
                            do k=1,problem%n
                              select case (node(sn)%ctype(k,1))
                                ! Global axes
                                ! 0 - perfect bonding - u_k^{1} and t_k^{1} are active since:
                                ! t_k^{(2)}=-t_k^{(1)}
                                ! u_k^{(2)}= u_k^{(1)}
                                case (0)
                                  node(sn)%dvda_c(          k,1,:)=bsa_c(node(sn)%col(k,1),:)
                                  node(sn)%dvda_c(problem%n+k,1,:)=bsa_c(node(sn)%col(problem%n+k,1),:)
                                  node(sn)%dvda_c(          k,2,:)= node(sn)%dvda_c(          k,1,:)
                                  node(sn)%dvda_c(problem%n+k,2,:)=-node(sn)%dvda_c(problem%n+k,1,:)
                                ! Global axes
                                ! 1 - perfect debonding - u_k^{1} and u_k^{2} are active since:
                                ! t_k^{(1)}=0
                                ! t_k^{(2)}=0
                                case (1)
                                  node(sn)%dvda_c(          k,1,:)=bsa_c(node(sn)%col(k,1),:)
                                  node(sn)%dvda_c(problem%n+k,1,:)=0.0d0
                                  node(sn)%dvda_c(          k,2,:)=bsa_c(node(sn)%col(k,2),:)
                                  node(sn)%dvda_c(problem%n+k,2,:)=0.0d0
                                case (2)
                                  stop 'not yet 31'
                                case (3,4,5)
                                  node(sn)%dvda_c(          k,1,:)=bsa_c(node(sn)%col(k,1),:)
                                  node(sn)%dvda_c(problem%n+k,1,:)=bsa_c(node(sn)%col(problem%n+k,1),:)
                                  node(sn)%dvda_c(          k,2,:)=bsa_c(node(sn)%col(k,2),:)
                                  node(sn)%dvda_c(problem%n+k,2,:)=-node(sn)%dvda_c(problem%n+k,1,:)
                              end select
                            end do

                          ! -------------------------------------------------------
                          ! BE (1: viscoelastic solid) - BE (2: poroelastic medium)
                          ! -------------------------------------------------------

                          case (fbem_poroelastic)
                            stop 'not yet 32'

                        end select

                      case (fbem_poroelastic)
                        stop 'not yet 33'

                    end select

                  ! ==============
                  ! BE-FE BOUNDARY
                  ! ==============

                  case (fbem_boundary_coupling_be_fe)
                    ! The fe node connected with the be node
                    sn_fe=element(se)%element_node(kn)
                    select case (boundary(sb)%class)

                      ! =================
                      ! ORDINARY BOUNDARY
                      ! =================

                      case (fbem_boundary_class_ordinary)
                        select case (region(kr)%type)

                          ! --------------
                          ! INVISCID FLUID
                          ! --------------

                          case (fbem_potential)
                            ! Index of variables:
                            ! node(sn)%col(1,1): p
                            ! node(sn)%col(2,1): Un (inactive since Un=u^(fe)·n)
                            !
                            ! p
                            node(sn)%dvda_c(1,1,:)=bsa_c(node(sn)%col(1,1),:)
                            ! Un=u^(fe)·n
                            ! Un initialization
                            node(sn)%dvda_c(2,1,:)=0.0d0
                            ! If the boundary is not reversed
                            if (sb_reversion.eqv.(.false.)) then
                              do k=1,problem%n
                                if (node(sn_fe)%ctype(k,1).eq.1) then
                                  node(sn)%dvda_c(2,1,:)=node(sn)%dvda_c(2,1,:)+bsa_c(node(sn_fe)%col(k,1),:)*n_fn(k)
                                end if
                              end do
                            ! If the boundary is reversed
                            else
                              do k=1,problem%n
                                if (node(sn_fe)%ctype(k,1).eq.1) then
                                  node(sn)%dvda_c(2,1,:)=node(sn)%dvda_c(2,1,:)-bsa_c(node(sn_fe)%col(k,1),:)*n_fn(k)
                                end if
                              end do
                            end if

                          ! ------------------
                          ! VISCOELASTIC SOLID
                          ! ------------------

                          case (fbem_viscoelastic)
                            ! Index of variables for each coordinate k:
                            ! node(sn)%col(          k,1): u_k
                            ! node(sn)%col(problem%n+k,1): t_k
                            do k=1,problem%n
                              if (node(sn_fe)%ctype(k,1).eq.1) then
                                node(sn)%dvda_c(k,1,:)=bsa_c(node(sn_fe)%col(k,1),:)
                              else
                                node(sn)%dvda_c(k,1,:)=(0.d0,0.d0)
                              end if
                              node(sn)%dvda_c(problem%n+k,1,:)=bsa_c(node(sn)%col(problem%n+k,1),:)
                            end do

                          ! ------------------
                          ! POROELASTIC MEDIUM
                          ! ------------------

                          case (fbem_poroelastic)
                            stop 'not yet 34'

                        end select

                      ! ===================
                      ! CRACK-LIKE BOUNDARY
                      ! ===================

                      case (fbem_boundary_class_cracklike)

                        select case (region(kr)%type)

                          ! --------------
                          ! INVISCID FLUID
                          ! --------------

                          case (fbem_potential)
                            ! Index of variables:
                            ! node(sn)%col(1,1): p^+
                            ! node(sn)%col(2,1): Un^+ (inactive since Un^+= u^(fe)·n)
                            ! node(sn)%col(1,2): p^-
                            ! node(sn)%col(2,2): Un^- (inactive since Un^-=-u^(fe)·n)
                            !
                            ! p^+
                            node(sn)%dvda_c(1,1,:)=bsa_c(node(sn)%col(1,1),:)
                            ! p^-
                            node(sn)%dvda_c(1,2,:)=bsa_c(node(sn)%col(1,2),:)
                            ! Un^+= u^(fe)·n
                            ! The fe node connected with the be node
                            sn_fe=element(se)%element_node(kn)
                            ! Un^+ initialization
                            node(sn)%dvda_c(2,1,:)=0.0d0
                            ! If the boundary is not reversed
                            if (sb_reversion.eqv.(.false.)) then
                              do k=1,problem%n
                                if (node(sn_fe)%ctype(k,1).eq.1) then
                                  node(sn)%dvda_c(2,1,:)=node(sn)%dvda_c(2,1,:)+bsa_c(node(sn_fe)%col(k,1),:)*n_fn(k)
                                end if
                              end do
                            ! If the boundary is reversed
                            else
                              do k=1,problem%n
                                if (node(sn_fe)%ctype(k,1).eq.1) then
                                  node(sn)%dvda_c(2,1,:)=node(sn)%dvda_c(2,1,:)-bsa_c(node(sn_fe)%col(k,:),1)*n_fn(k)
                                end if
                              end do
                            end if
                            ! Un^-=-Un^+
                            node(sn)%dvda_c(2,2,:)=-node(sn)%dvda_c(2,1,:)

                          ! ------------------
                          ! VISCOELASTIC SOLID
                          ! ------------------

                          case (fbem_viscoelastic)
                            ! Index of variables for each coordinate k:
                            ! node(sn)%col(          k,1): u_k for face +
                            ! node(sn)%col(problem%n+k,1): t_k for face +
                            ! node(sn)%col(          k,2): u_k for face - (not used since u_l_{+} = u_k_{-})
                            ! node(sn)%col(problem%n+k,2): t_k for face -
                            do k=1,problem%n
                              if (node(sn_fe)%ctype(k,1).eq.1) then
                                node(sn)%dvda_c(k,1,:)=bsa_c(node(sn_fe)%col(k,1),:)
                              else
                                node(sn)%dvda_c(k,1,:)=(0.d0,0.d0)
                              end if
                              node(sn)%dvda_c(problem%n+k,1,:)=bsa_c(node(sn)%col(problem%n+k,1),:)
                              node(sn)%dvda_c(          k,2,:)=node(sn)%dvda_c(k,1,:)
                              node(sn)%dvda_c(problem%n+k,2,:)=bsa_c(node(sn)%col(problem%n+k,2),:)
                            end do

                          ! ------------------
                          ! POROELASTIC MEDIUM
                          ! ------------------

                          case (fbem_poroelastic)
                            stop 'not yet 35'

                        end select

                    end select

                  ! =================
                  ! BE-FE-BE BOUNDARY
                  ! =================

                  case (fbem_boundary_coupling_be_fe_be)
                    stop 'not yet 36'

                end select
              end if
            end do
          end do
        end do

      ! =========
      ! FE REGION
      ! =========

      case (fbem_fe)
        do ks=1,region(kr)%n_fe_subregions
          ss=region(kr)%fe_subregion(ks)
          do ke=1,part(fe_subregion(ss)%part)%n_elements
            se=part(fe_subregion(ss)%part)%element(ke)
            do kn=1,element(se)%n_nodes
              sn=element(se)%node(kn)
              if (node_used(sn).eqv.(.false.)) then
                node_used(sn)=.true.
                do k=1,node(sn)%n_dof
                  if (node(sn)%ctype(k,1).eq.1) then
                    node(sn)%dvda_c(k,1,:)=bsa_c(node(sn)%col(k,1),:)
                  else
                    node(sn)%dvda_c(k,1,:)=(0.d0,0.d0)
                  end if
                end do
              end if
            end do
          end do
        end do

    end select

  end do

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'START assigning the LSE solution to variables (sensitivity)')

end subroutine assign_solution_mechanics_harmonic_sa
