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

subroutine assign_solution_mechanics_harmonic(kf)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_data_structures
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_geometry
  use fbem_symmetry
  use fbem_harela_incident_field

  ! Module of problem variables
  use problem_variables

  ! No implicit variables are allowed
  implicit none

  ! I/O variables
  integer                    :: kf

  ! Local variables
  real(kind=real64)          :: omega
  integer                    :: kr, kb, ke, kn, ks, ss
  integer                    :: sb, sp, se, sn, km, snm
  logical                    :: sb_reversion
  integer                    :: k, kc
  logical, allocatable       :: node_used(:)
  integer                    :: sn_fe
  real(kind=real64)          :: x_fn(problem%n), n_fn(problem%n)
  complex(kind=real64)       :: KT
  real(kind=real64)          :: phi, phi1, phi2, rho, d_phi1_phi2, d_phi2_phi1, p_1_phi1_phi2, p_1_phi2_phi1, ctephi
  complex(kind=real64)       :: c
  real(kind=real64), allocatable :: T(:,:)

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'START assigning the LSE solution to variables')

  omega=frequency(kf)
  allocate (node_used(n_nodes))
  node_used=.false.

  ! =============
  ! RIGID REGIONS
  ! =============

  do kr=1,n_regions
    if ((region(kr)%class.eq.fbem_fe).and.(region(kr)%type.eq.fbem_rigid)) then
      !
      ! Solution of the master node
      !
      sn=region(kr)%master_node
      node_used(sn)=.true.
      do k=1,node(sn)%n_dof
        if (node(sn)%ctype(k,1).eq.1) then
          node(sn)%value_c(k,1)=b_c(node(sn)%col(k,1),1)
        else
          node(sn)%value_c(k,1)=node(sn)%cvalue_c(k,1,1)
        end if
      end do
      !
      ! Solution of the slave nodes
      !
      snm=sn
      do kn=1,node(snm)%n_slaves
        sn=node(snm)%slave(kn)
        node_used(sn)=.true.
        allocate (T(node(sn)%n_dof,node(snm)%n_dof))
        call fbem_rigid_solid_transformation_matrix(problem%n,node(sn)%x,node(snm)%x,node(sn)%n_dof,T)
        do k=1,node(sn)%n_dof
          node(sn)%value_c(k,1)=0
          do km=1,node(snm)%n_dof
            node(sn)%value_c(k,1)=node(sn)%value_c(k,1)+T(k,km)*node(snm)%value_c(km,1)
          end do
        end do
        deallocate (T)
      end do
    end if
  end do

  ! ================
  ! FLEXIBLE REGIONS
  ! ================

  do kr=1,n_regions
    select case (region(kr)%class)

      ! =========
      ! BE REGION
      ! =========

      case (fbem_be)

        ! ==========
        ! BOUNDARIES
        ! ==========

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
                                node(sn)%value_c(1,1)=node(sn)%cvalue_c(1,1,1)
                                node(sn)%value_c(2,1)=b_c(node(sn)%col(2,1),1)
                              case (1)
                                node(sn)%value_c(1,1)=b_c(node(sn)%col(1,1),1)
                                node(sn)%value_c(2,1)=node(sn)%cvalue_c(1,1,1)
                              case (2)
                                rho=region(kr)%property_r(1)
                                c=region(kr)%property_c(4)
                                node(sn)%value_c(1,1)=b_c(node(sn)%col(1,1),1)
                                node(sn)%value_c(2,1)=c_im/rho/c/omega*node(sn)%value_c(1,1)
                              ! Un-(i/(rho*c*omega)-1/(2*R*rho*omega^2))p=0
                              case (3)
                                rho=region(kr)%property_r(1)
                                c=region(kr)%property_c(4)
                                node(sn)%value_c(1,1)=b_c(node(sn)%col(1,1),1)
                                node(sn)%value_c(2,1)=(c_im/rho/c/omega-1./(2.*node(sn)%cvalue_c(1,1,1)*rho*omega**2))*node(sn)%value_c(1,1)
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
                                  node(sn)%value_c(          k,1)=node(sn)%cvalue_c(k,1,1)
                                  node(sn)%value_c(problem%n+k,1)=b_c(node(sn)%col(problem%n+k,1),1)
                                case (1)
                                  node(sn)%value_c(          k,1)=b_c(node(sn)%col(k,1),1)
                                  node(sn)%value_c(problem%n+k,1)=node(sn)%cvalue_c(k,1,1)
                                case (10)
                                  node(sn)%value_c(          k,1)=b_c(node(sn)%col(k,1),1)
                                  if (.not.region(kr)%boundary_reversion(kb)) then
                                    node(sn)%value_c(problem%n+k,1)= node(sn)%cvalue_c(k,1,1)*node(sn)%n_fn(k)
                                  else
                                    node(sn)%value_c(problem%n+k,1)=-node(sn)%cvalue_c(k,1,1)*node(sn)%n_fn(k)
                                  end if
                                case (2,3)
                                  node(sn)%value_c(          k,1)=b_c(node(sn)%col(k,1),1)
                                  node(sn)%value_c(problem%n+k,1)=b_c(node(sn)%col(problem%n+k,1),1)
                              end select
                            end do

                          ! ------------------
                          ! POROELASTIC MEDIUM
                          ! ------------------

                          case (fbem_poroelastic)
                            ! Index of variables:
                            ! node(sn)%col(            0,1): tau
                            ! node(sn)%col(            k,1): u_k
                            ! node(sn)%col(problem%n+1+0,1): Un
                            ! node(sn)%col(problem%n+1+k,1): t_k
                            !
                            ! Fluid phase variables
                            !
                            k=0
                            select case (node(sn)%ctype(k,1))
                              !
                              ! open pore / permeable
                              !
                              ! tau known, Un unknown
                              case (0)
                                node(sn)%value_c(            k,1)=node(sn)%cvalue_c(k,1,1)
                                node(sn)%value_c(problem%n+1+k,1)=b_c(node(sn)%col(problem%n+1+k,1),1)
                              ! Un known, tau unknown
                              case (1)
                                node(sn)%value_c(            k,1)=b_c(node(sn)%col(k,1),1)
                                node(sn)%value_c(problem%n+1+k,1)=node(sn)%cvalue_c(k,1,1)
                              !
                              ! close pore / impermeable
                              !
                              ! Un=u_k·n_k known, tau unknown
                              case (2)
                                ! tau
                                node(sn)%value_c(            k,1)=b_c(node(sn)%col(k,1),1)
                                ! Un=u_k·n_k
                                node(sn)%value_c(problem%n+1+k,1)=0.0d0
                                ! n_fn is positive
                                if (.not.sb_reversion) then
                                  do kc=1,problem%n
                                    select case (node(sn)%ctype(kc,1))
                                      ! u_k known
                                      case (4)
                                        node(sn)%value_c(problem%n+1+k,1)=node(sn)%value_c(problem%n+1+k,1)&
                                                                         +node(sn)%cvalue_c(kc,1,1)*n_fn(kc)
                                      ! u_k unknown
                                      case (5,50,6,7)
                                        node(sn)%value_c(problem%n+1+k,1)=node(sn)%value_c(problem%n+1+k,1)&
                                                                         +b_c(node(sn)%col(kc,1),1)*n_fn(kc)
                                    end select
                                  end do
                                ! n_fn is negative
                                else
                                  do kc=1,problem%n
                                    select case (node(sn)%ctype(kc,1))
                                      ! u_k known
                                      case (4)
                                        node(sn)%value_c(problem%n+1+k,1)=node(sn)%value_c(problem%n+1+k,1)&
                                                                         -node(sn)%cvalue_c(kc,1,1)*n_fn(kc)
                                      ! u_k unknown
                                      case (5,50,6,7)
                                        node(sn)%value_c(problem%n+1+k,1)=node(sn)%value_c(problem%n+1+k,1)&
                                                                         -b_c(node(sn)%col(kc,1),1)*n_fn(kc)
                                    end select
                                  end do
                                end if
                            end select
                            !
                            ! Solid phase variables
                            !
                            do k=1,problem%n
                              select case (node(sn)%ctype(k,1))
                                !
                                ! open pore / permeable
                                !
                                ! u_k known, t_k unknown
                                case (0)
                                  node(sn)%value_c(            k,1)=node(sn)%cvalue_c(k,1,1)
                                  node(sn)%value_c(problem%n+1+k,1)=b_c(node(sn)%col(problem%n+1+k,1),1)
                                ! t_k known, u_k unknown
                                case (1)
                                  node(sn)%value_c(            k,1)=b_c(node(sn)%col(            k,1),1)
                                  node(sn)%value_c(problem%n+1+k,1)=node(sn)%cvalue_c(k,1,1)
                                ! t_k=Pn_k known, u_k unknown
                                case (10)
                                  node(sn)%value_c(          k,1)=b_c(node(sn)%col(k,1),1)
                                  if (.not.region(kr)%boundary_reversion(kb)) then
                                    node(sn)%value_c(problem%n+1+k,1)= node(sn)%cvalue_c(k,1,1)*node(sn)%n_fn(k)
                                  else
                                    node(sn)%value_c(problem%n+1+k,1)=-node(sn)%cvalue_c(k,1,1)*node(sn)%n_fn(k)
                                  end if
                                ! u_k unknown, t_k unknown
                                case (2,3)
                                  node(sn)%value_c(            k,1)=b_c(node(sn)%col(            k,1),1)
                                  node(sn)%value_c(problem%n+1+k,1)=b_c(node(sn)%col(problem%n+1+k,1),1)
                                !
                                ! close pore / impermeable
                                !
                                ! u_k known, t_k unknown
                                case (4)
                                  node(sn)%value_c(            k,1)=node(sn)%cvalue_c(k,1,1)
                                  node(sn)%value_c(problem%n+1+k,1)=b_c(node(sn)%col(problem%n+1+k,1),1)
                                ! t_k+tau·n_k=T known, u_k unknown
                                case (5)
                                  ! u_k unknown
                                  node(sn)%value_c(            k,1)=b_c(node(sn)%col(k,1),1)
                                  ! t_k=T-tau·n_k
                                  node(sn)%value_c(problem%n+1+k,1)=node(sn)%cvalue_c(k,1,1)
                                  ! n_fn is positive
                                  if (.not.sb_reversion) then
                                    node(sn)%value_c(problem%n+1+k,1)=node(sn)%value_c(problem%n+1+k,1)&
                                                                     -node(sn)%value_c(0,1)*n_fn(k)
                                  ! n_fn is negative
                                  else
                                    node(sn)%value_c(problem%n+1+k,1)=node(sn)%value_c(problem%n+1+k,1)&
                                                                     +node(sn)%value_c(0,1)*n_fn(k)
                                  end if
                                ! t_k + tau n_k = P n_k known, u_k unknown
                                case (50)
                                  node(sn)%value_c(          k,1)=b_c(node(sn)%col(k,1),1)
                                  if (.not.region(kr)%boundary_reversion(kb)) then
                                    node(sn)%value_c(problem%n+1+k,1)= (node(sn)%cvalue_c(k,1,1)-node(sn)%value_c(0,1))*node(sn)%n_fn(k)
                                  else
                                    node(sn)%value_c(problem%n+1+k,1)=-(node(sn)%cvalue_c(k,1,1)-node(sn)%value_c(0,1))*node(sn)%n_fn(k)
                                  end if
                                ! u_k unknown, t_k unknown
                                case (6,7)
                                  node(sn)%value_c(            k,1)=b_c(node(sn)%col(            k,1),1)
                                  node(sn)%value_c(problem%n+1+k,1)=b_c(node(sn)%col(problem%n+1+k,1),1)
                              end select
                            end do
                        end select

!                        ! Change of variables for traction quarter point elements
!                        if (element(se)%type.eq.fbem_line3) then
!                          select case (kn)
!                            ! Node 1
!                            case (1)
!                              select case (element(se)%type_f1f)
!                                case (fbem_line3_qp1t,fbem_line3_mqp1t)
!                                  node(sn)%value_c(0,1)=huge(1.d0)
!                              end select
!                              select case (element(se)%type_f1)
!                                case (fbem_line3_qp1t,fbem_line3_mqp1t)
!                                  do k=1,problem%n
!                                    node(sn)%value_c(k,1)=huge(1.d0)
!                                  end do
!                              end select
!                              select case (element(se)%type_f2f)
!                                case (fbem_line3_qp1t,fbem_line3_mqp1t)
!                                  node(sn)%value_c(problem%n+1,1)=huge(1.d0)
!                              end select
!                              select case (element(se)%type_f2)
!                                case (fbem_line3_qp1t,fbem_line3_mqp1t)
!                                  do k=1,problem%n
!                                    node(sn)%value_c(problem%n+1+k,1)=huge(1.d0)
!                                  end do
!                              end select
!                            ! Node 2
!                            case (2)
!                              select case (element(se)%type_f1f)
!                                case (fbem_line3_qp2t,fbem_line3_mqp2t)
!                                  node(sn)%value_c(0,1)=huge(1.d0)
!                              end select
!                              select case (element(se)%type_f1)
!                                case (fbem_line3_qp2t,fbem_line3_mqp2t)
!                                  do k=1,problem%n
!                                    node(sn)%value_c(k,1)=huge(1.d0)
!                                  end do
!                              end select
!                              select case (element(se)%type_f2f)
!                                case (fbem_line3_qp2t,fbem_line3_mqp2t)
!                                  node(sn)%value_c(problem%n+1,1)=huge(1.d0)
!                              end select
!                              select case (element(se)%type_f2)
!                                case (fbem_line3_qp2t,fbem_line3_mqp2t)
!                                  do k=1,problem%n
!                                    node(sn)%value_c(problem%n+1+k,1)=huge(1.d0)
!                                  end do
!                              end select
!                            ! Node 3
!                            case (3)
!                              select case (element(se)%type_f1f)
!                                case (fbem_line3_qp1t,fbem_line3_qp2t,fbem_line3_mqp1t,fbem_line3_mqp2t)
!                                  node(sn)%value_c(0,1)=2.d0*node(sn)%value_c(0,1)
!                              end select
!                              select case (element(se)%type_f1)
!                                case (fbem_line3_qp1t,fbem_line3_qp2t,fbem_line3_mqp1t,fbem_line3_mqp2t)
!                                  do k=1,problem%n
!                                    node(sn)%value_c(k,1)=2.d0*node(sn)%value_c(k,1)
!                                  end do
!                              end select
!                              select case (element(se)%type_f2f)
!                                case (fbem_line3_qp1t,fbem_line3_qp2t,fbem_line3_mqp1t,fbem_line3_mqp2t)
!                                  node(sn)%value_c(problem%n+1,1)=2.d0*node(sn)%value_c(problem%n+1,1)
!                              end select
!                              select case (element(se)%type_f2)
!                                case (fbem_line3_qp1t,fbem_line3_qp2t,fbem_line3_mqp1t,fbem_line3_mqp2t)
!                                  do k=1,problem%n
!                                    node(sn)%value_c(problem%n+1+k,1)=2.d0*node(sn)%value_c(problem%n+1+k,1)
!                                  end do
!                              end select
!                          end select
!                        end if

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
                                node(sn)%value_c(1,1)=node(sn)%cvalue_c(1,1,1)
                                node(sn)%value_c(2,1)=b_c(node(sn)%col(2,1),1)
                              case (1)
                                node(sn)%value_c(1,1)=b_c(node(sn)%col(1,1),1)
                                node(sn)%value_c(2,1)=node(sn)%cvalue_c(1,1,1)
                            end select
                            ! Face -
                            select case (node(sn)%ctype(1,2))
                              case (0)
                                node(sn)%value_c(1,2)=node(sn)%cvalue_c(1,1,2)
                                node(sn)%value_c(2,2)=b_c(node(sn)%col(2,2),1)
                              case (1)
                                node(sn)%value_c(1,2)=b_c(node(sn)%col(1,2),1)
                                node(sn)%value_c(2,2)=node(sn)%cvalue_c(1,1,2)
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
                                  node(sn)%value_c(          k,1)=node(sn)%cvalue_c(k,1,1)
                                  node(sn)%value_c(problem%n+k,1)=b_c(node(sn)%col(problem%n+k,1),1)
                                case (1)
                                  node(sn)%value_c(          k,1)=b_c(node(sn)%col(k,1),1)
                                  node(sn)%value_c(problem%n+k,1)=node(sn)%cvalue_c(k,1,1)
                              end select
                            end do
                            ! Face -
                            do k=1,problem%n
                              select case (node(sn)%ctype(k,2))
                                case (0)
                                  node(sn)%value_c(          k,2)=node(sn)%cvalue_c(k,1,2)
                                  node(sn)%value_c(problem%n+k,2)=b_c(node(sn)%col(problem%n+k,2),1)
                                case (1)
                                  node(sn)%value_c(          k,2)=b_c(node(sn)%col(k,2),1)
                                  node(sn)%value_c(problem%n+k,2)=node(sn)%cvalue_c(k,1,2)
                              end select
                            end do

                          ! ------------------
                          ! POROELASTIC MEDIUM
                          ! ------------------

                          case (fbem_poroelastic)
                            ! Index of variables:
                            ! node(sn)%col(            0,1): tau for face +
                            ! node(sn)%col(            k,1): u_k for face +
                            ! node(sn)%col(problem%n+1  ,1): Un  for face +
                            ! node(sn)%col(problem%n+1+k,1): t_k for face +
                            ! node(sn)%col(            0,2): tau for face -
                            ! node(sn)%col(            k,2): u_k for face -
                            ! node(sn)%col(problem%n+1  ,2): Un  for face -
                            ! node(sn)%col(problem%n+1+k,2): t_k for face -
                            !
                            ! Face +
                            !
                            do k=0,problem%n
                              select case (node(sn)%ctype(k,1))
                                case (0)
                                  node(sn)%value_c(            k,1)=node(sn)%cvalue_c(k,1,1)
                                  node(sn)%value_c(problem%n+1+k,1)=b_c(node(sn)%col(problem%n+1+k,1),1)
                                case (1)
                                  node(sn)%value_c(            k,1)=b_c(node(sn)%col(k,1),1)
                                  node(sn)%value_c(problem%n+1+k,1)=node(sn)%cvalue_c(k,1,1)
                                ! close pore / impermeable
                                ! If k==0: Un^+=u_k^+·n_k^+ known, tau^+ unknown
                                ! If k>=1: u_k^+ known, t_k^+ unknown
                                case (2)
                                  ! Fluid phase variables
                                  if (k.eq.0) then
                                    node(sn)%value_c(            k,1)=b_c(node(sn)%col(k,1),1)
                                    ! Un^+=u_k^+·n_k^+
                                    node(sn)%value_c(problem%n+1+k,1)=0.0d0
                                    ! n_fn is positive
                                    if (.not.sb_reversion) then
                                      do kc=1,problem%n
                                        select case (node(sn)%ctype(kc,1))
                                          ! u_k^+ known
                                          case (2)
                                            node(sn)%value_c(problem%n+1+k,1)=node(sn)%value_c(problem%n+1+k,1)&
                                                                             +node(sn)%cvalue_c(kc,1,1)*n_fn(kc)
                                          ! u_k^+ unknown
                                          case (3)
                                            node(sn)%value_c(problem%n+1+k,1)=node(sn)%value_c(problem%n+1+k,1)&
                                                                             +b_c(node(sn)%col(kc,1),1)*n_fn(kc)
                                        end select
                                      end do
                                    ! n_fn is negative
                                    else
                                      do kc=1,problem%n
                                        select case (node(sn)%ctype(kc,1))
                                          ! u_k^+ known
                                          case (2)
                                            node(sn)%value_c(problem%n+1+k,1)=node(sn)%value_c(problem%n+1+k,1)&
                                                                             -node(sn)%cvalue_c(kc,1,1)*n_fn(kc)
                                          ! u_k^+ unknown
                                          case (3)
                                            node(sn)%value_c(problem%n+1+k,1)=node(sn)%value_c(problem%n+1+k,1)&
                                                                             -b_c(node(sn)%col(kc,1),1)*n_fn(kc)
                                        end select
                                      end do
                                    end if
                                  ! Solid skeleton variables
                                  else
                                    node(sn)%value_c(            k,1)=node(sn)%cvalue_c(k,1,1)
                                    node(sn)%value_c(problem%n+1+k,1)=b_c(node(sn)%col(problem%n+1+k,1),1)
                                  end if
                                ! close pore / impermeable
                                ! If k>=1: t_k^++tau^+·n_k^+ known, u_k^+ unknown
                                case (3)
                                  ! u_k^+ unknown
                                  node(sn)%value_c(            k,1)=b_c(node(sn)%col(k,1),1)
                                  ! t_k^+=T^+-tau^+·n_k^+
                                  node(sn)%value_c(problem%n+1+k,1)=node(sn)%cvalue_c(k,1,1)
                                  ! n_fn is positive
                                  if (.not.sb_reversion) then
                                    node(sn)%value_c(problem%n+1+k,1)=node(sn)%value_c(problem%n+1+k,1)&
                                                                     -node(sn)%value_c(0,1)*n_fn(k)
                                  ! n_fn is negative
                                  else
                                    node(sn)%value_c(problem%n+1+k,1)=node(sn)%value_c(problem%n+1+k,1)&
                                                                     +node(sn)%value_c(0,1)*n_fn(k)
                                  end if
                                ! Fluid(inviscid/incompressible)-filled permeable crack
                                ! If k==0: tau^+ unknown, Un^+ unknown
                                ! If k>=1: u^+ unknown, t_k^+=(1/phi-1)tau^+·n_k^+ known
                                case (4)
                                  if (k.eq.0) then
                                    ! tau^+ unknown
                                    node(sn)%value_c(            k,1)=b_c(node(sn)%col(k,1),1)
                                    ! Un^+ unknown
                                    node(sn)%value_c(problem%n+1+k,1)=b_c(node(sn)%col(problem%n+1+k,1),1)
                                  else
                                    ! u_k^+ unknown
                                    node(sn)%value_c(            k,1)=b_c(node(sn)%col(k,1),1)
                                    ! Calculate 1/phi-1
                                    phi=region(boundary(sb)%region(1))%property_r(8)
                                    ctephi=1.0d0/phi-1.0d0
                                    ! t_k^+=(1/phi-1)tau^+·n_k^+
                                    if (.not.sb_reversion) then
                                      node(sn)%value_c(problem%n+1+k,1)= ctephi*node(sn)%value_c(0,1)*n_fn(k)
                                    else
                                      node(sn)%value_c(problem%n+1+k,1)=-ctephi*node(sn)%value_c(0,1)*n_fn(k)
                                    end if
                                  end if
                                ! Fluid(inviscid/incompressible)-filled impermeable crack
                                ! If k==0: tau^+ unknown, Un^+=u^+·n^+ known
                                ! If k>=1: u^+ unknown, t_k^+ unknown
                                case (5)
                                  if (k.eq.0) then
                                    ! tau^+ unknown
                                    node(sn)%value_c(            k,1)=b_c(node(sn)%col(k,1),1)
                                    ! Un^+=u_k^+·n_k^+
                                    node(sn)%value_c(problem%n+1+k,1)=0.0d0
                                    ! Boundary orientation is positive
                                    if (.not.sb_reversion) then
                                      do kc=1,problem%n
                                        node(sn)%value_c(problem%n+1+k,1)=node(sn)%value_c(problem%n+1+k,1)&
                                                                          +b_c(node(sn)%col(kc,1),1)*n_fn(kc)
                                      end do
                                    ! Boundary orientation is negative
                                    else
                                      do kc=1,problem%n
                                        node(sn)%value_c(problem%n+1+k,1)=node(sn)%value_c(problem%n+1+k,1)&
                                                                          -b_c(node(sn)%col(kc,1),1)*n_fn(kc)
                                      end do
                                    end if
                                  else
                                    ! u_k^+ unknown
                                    node(sn)%value_c(            k,1)=b_c(node(sn)%col(            k,1),1)
                                    ! t_k^+ unknown
                                    node(sn)%value_c(problem%n+1+k,1)=b_c(node(sn)%col(problem%n+1+k,1),1)
                                  end if
                              end select
                            end do
                            !
                            ! Face -
                            !
                            do k=0,problem%n
                              select case (node(sn)%ctype(k,2))
                                case (0)
                                  node(sn)%value_c(            k,2)=node(sn)%cvalue_c(k,1,2)
                                  node(sn)%value_c(problem%n+1+k,2)=b_c(node(sn)%col(problem%n+1+k,2),1)
                                case (1)
                                  node(sn)%value_c(            k,2)=b_c(node(sn)%col(k,2),1)
                                  node(sn)%value_c(problem%n+1+k,2)=node(sn)%cvalue_c(k,1,2)
                                ! close pore / impermeable
                                ! If k==0: Un^-=u_k^-·n_k^- known, tau^- unknown
                                ! If k>=1: u_k^- known, t_k^- unknown
                                case (2)
                                  ! Fluid phase variables
                                  if (k.eq.0) then
                                    node(sn)%value_c(            k,2)=b_c(node(sn)%col(k,2),1)
                                    ! Un^-=u_k^-·n_k^-
                                    node(sn)%value_c(problem%n+1+k,2)=0.0d0
                                    ! n_fn is negative
                                    if (.not.sb_reversion) then
                                      do kc=1,problem%n
                                        select case (node(sn)%ctype(kc,2))
                                          ! u_k^- known
                                          case (2)
                                            node(sn)%value_c(problem%n+1+k,2)=node(sn)%value_c(problem%n+1+k,2)&
                                                                             -node(sn)%cvalue_c(kc,1,2)*n_fn(kc)
                                          ! u_k^- unknown
                                          case (3)
                                            node(sn)%value_c(problem%n+1+k,2)=node(sn)%value_c(problem%n+1+k,2)&
                                                                             -b_c(node(sn)%col(kc,2),1)*n_fn(kc)
                                        end select
                                      end do
                                    ! n_fn is positive
                                    else
                                      do kc=1,problem%n
                                        select case (node(sn)%ctype(kc,2))
                                          ! u_k^- known
                                          case (2)
                                            node(sn)%value_c(problem%n+1+k,2)=node(sn)%value_c(problem%n+1+k,2)&
                                                                             +node(sn)%cvalue_c(kc,1,2)*n_fn(kc)
                                          ! u_k^- unknown
                                          case (3)
                                            node(sn)%value_c(problem%n+1+k,2)=node(sn)%value_c(problem%n+1+k,2)&
                                                                             +b_c(node(sn)%col(kc,2),1)*n_fn(kc)
                                        end select
                                      end do
                                    end if
                                  ! Solid skeleton variables
                                  else
                                    node(sn)%value_c(            k,2)=node(sn)%cvalue_c(k,1,2)
                                    node(sn)%value_c(problem%n+1+k,2)=b_c(node(sn)%col(problem%n+1+k,2),1)
                                  end if
                                ! close pore / impermeable
                                ! Only k>=1: t_k^-+tau^-·n_k^- known, u_k^- unknown
                                case (3)
                                  ! u_k^- unknown
                                  node(sn)%value_c(            k,2)=b_c(node(sn)%col(k,2),1)
                                  ! t_k^-=T^--tau^-·n_k^-
                                  node(sn)%value_c(problem%n+1+k,2)=node(sn)%cvalue_c(k,1,2)
                                  ! n_fn is negative
                                  if (.not.sb_reversion) then
                                    node(sn)%value_c(problem%n+1+k,2)=node(sn)%value_c(problem%n+1+k,2)&
                                                                     +node(sn)%value_c(0,2)*n_fn(k)
                                  ! n_fn is positive
                                  else
                                    node(sn)%value_c(problem%n+1+k,2)=node(sn)%value_c(problem%n+1+k,2)&
                                                                     -node(sn)%value_c(0,2)*n_fn(k)
                                  end if
                                ! Fluid(inviscid/incompressible)-filled permeable crack
                                ! If k==0: tau^- and Un^- known in function of other variables
                                ! If k>=1: u^- unknown, t_k^-=-(1/phi-1)tau^+·n_k^+ known
                                case (4)
                                  if (k.eq.0) then
                                    ! tau^-=tau^+
                                    node(sn)%value_c(            k,2)=node(sn)%value_c(            k,1)
                                    ! Un^-= - Un^+ - (1/phi+1)u^+·n^+ + (1/phi+1)u^-·n^+
                                    ! Calculate 1/phi-1
                                    phi=region(boundary(sb)%region(1))%property_r(8)
                                    ctephi=1.0d0/phi-1.0d0
                                    ! Add -Un^+
                                    node(sn)%value_c(problem%n+1+k,2)=-node(sn)%value_c(problem%n+1+k,1)
                                    ! Add u^+ and u^- components
                                    ! Boundary orientation is positive
                                    if (.not.sb_reversion) then
                                      do kc=1,problem%n
                                        ! Add -ctephi*u_kc^+*n_kc
                                        node(sn)%value_c(problem%n+1+k,2)=node(sn)%value_c(problem%n+1+k,2)&
                                                                          -ctephi*node(sn)%value_c(kc,1)*n_fn(kc)
                                        ! Add ctephi*u_kc^-*n_kc
                                        node(sn)%value_c(problem%n+1+k,2)=node(sn)%value_c(problem%n+1+k,2)&
                                                                          +ctephi*b_c(node(sn)%col(kc,2),1)*n_fn(kc)
                                      end do
                                    ! Boundary orientation is negative
                                    else
                                      do kc=1,problem%n
                                        ! Add -ctephi*u_kc^+*n_kc
                                        node(sn)%value_c(problem%n+1+k,2)=node(sn)%value_c(problem%n+1+k,2)&
                                                                          +ctephi*node(sn)%value_c(kc,1)*n_fn(kc)
                                        ! Add ctephi*u_kc^-*n_kc
                                        node(sn)%value_c(problem%n+1+k,2)=node(sn)%value_c(problem%n+1+k,2)&
                                                                          -ctephi*b_c(node(sn)%col(kc,2),1)*n_fn(kc)
                                      end do
                                    end if
                                  else
                                    ! u_k^- unknown
                                    node(sn)%value_c(            k,2)=b_c(node(sn)%col(k,2),1)
                                    ! Calculate 1/phi-1
                                    phi=region(boundary(sb)%region(1))%property_r(8)
                                    ctephi=1.0d0/phi-1.0d0
                                    ! t_k^-=-(1/phi-1)tau^+·n_k^+
                                    if (.not.sb_reversion) then
                                      node(sn)%value_c(problem%n+1+k,2)=-ctephi*node(sn)%value_c(0,1)*n_fn(k)
                                    else
                                      node(sn)%value_c(problem%n+1+k,2)=+ctephi*node(sn)%value_c(0,1)*n_fn(k)
                                    end if
                                  end if
                                ! Fluid(inviscid/incompressible)-filled impermeable crack
                                ! If k==0: tau^- unknown, Un^-=-u^-·n^+ known
                                ! If k>=1: u^- unknown, t_k^-=(tau^--tau^+)n_k^+-t_k^+ known
                                case (5)
                                  if (k.eq.0) then
                                    ! tau^- unknown
                                    node(sn)%value_c(            k,2)=b_c(node(sn)%col(k,2),1)
                                    ! Un^-=-Un^+
                                    node(sn)%value_c(problem%n+1+k,2)=-node(sn)%value_c(problem%n+1+k,1)
                                  else
                                    ! u_k^- unknown
                                    node(sn)%value_c(            k,2)=b_c(node(sn)%col(            k,2),1)
                                    ! t_k^-=(tau^--tau^+)n_k^+-t_k^+ known
                                    if (.not.sb_reversion) then
                                      node(sn)%value_c(problem%n+1+k,2)=(node(sn)%value_c(0,2)-node(sn)%value_c(0,1))*n_fn(k)
                                    else
                                      node(sn)%value_c(problem%n+1+k,2)=(node(sn)%value_c(0,1)-node(sn)%value_c(0,2))*n_fn(k)
                                    end if
                                    node(sn)%value_c(problem%n+1+k,2)=node(sn)%value_c(problem%n+1+k,2)&
                                                                     -node(sn)%value_c(problem%n+1+k,1)
                                  end if
                              end select
                            end do
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
                            node(sn)%value_c(1,1)=b_c(node(sn)%col(1,1),1)
                            node(sn)%value_c(2,1)=b_c(node(sn)%col(2,1),1)
                            node(sn)%value_c(1,2)= node(sn)%value_c(1,1)
                            node(sn)%value_c(2,2)=-node(sn)%value_c(2,1)

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
                            node(sn)%value_c(1,1)=b_c(node(sn)%col(1,1),1)
                            ! u_k of region 2
                            do k=1,problem%n
                              node(sn)%value_c(k,2)=b_c(node(sn)%col(k,2),1)
                            end do
                            ! Un of region 1
                            node(sn)%value_c(2,1)=0.0d0
                            do k=1,problem%n
                              node(sn)%value_c(2,1)=node(sn)%value_c(2,1)+node(sn)%value_c(k,2)*n_fn(k)
                            end do
                            ! t_k of region 2
                            do k=1,problem%n
                              node(sn)%value_c(problem%n+k,2)=node(sn)%value_c(1,1)*n_fn(k)
                            end do

                          ! ---------------------------------------------------
                          ! BE (1: inviscid fluid) - BE (2: poroelastic medium)
                          ! ---------------------------------------------------

                          case (fbem_poroelastic)
                            !
                            ! Switch depending on the B.C.
                            !
                            select case (node(sn)%ctype(1,1))
                              !
                              ! Perfectly permeable
                              !
                              case (0)
                                ! Index of variables:
                                ! node(sn)%col(            1,1): p   for region 1, (p^(1)=-tau^(2)/phi^(2))
                                ! node(sn)%col(            2,1): Un  for region 1, (Un^(1)=-phi^(2)*Un^(2)+(1-phi^(2))*u^(2)·n^(1))
                                ! node(sn)%col(            0,2): tau for region 2
                                ! node(sn)%col(            k,2): u_k for region 2
                                ! node(sn)%col(problem%n+1  ,2): Un  for region 2
                                ! node(sn)%col(problem%n+1+k,2): t_k for region 2
                                !
                                ! Region 2 variables (poroelastic medium)
                                phi=region(boundary(sb)%region(2))%property_r(8)
                                do k=0,problem%n
                                  node(sn)%value_c(            k,2)=b_c(node(sn)%col(            k,2),1)
                                  if (k.eq.0) then
                                    node(sn)%value_c(problem%n+1+k,2)=b_c(node(sn)%col(problem%n+1+k,2),1)
                                  else
                                    node(sn)%value_c(problem%n+1+k,2)=-(1.0d0-phi)/phi*node(sn)%value_c(0,2)*n_fn(k)
                                  end if
                                end do
                                ! Region 1 variables (fluid)
                                node(sn)%value_c(1,1)=-node(sn)%value_c(0,2)/phi
                                node(sn)%value_c(2,1)=-phi*node(sn)%value_c(problem%n+1,2)
                                do k=1,problem%n
                                  node(sn)%value_c(2,1)=node(sn)%value_c(2,1)+(1.0d0-phi)*node(sn)%value_c(k,2)*n_fn(k)
                                end do
                              !
                              ! Perfectly impermeable
                              !
                              case (1)
                                ! Index of variables:
                                ! node(sn)%col(            1,1): p   for region 1
                                ! node(sn)%col(            2,1): Un  for region 1, (Un^(1)=u^(2)·n^(1))
                                ! node(sn)%col(            0,2): tau for region 2
                                ! node(sn)%col(            k,2): u_k for region 2
                                ! node(sn)%col(problem%n+1  ,2): Un  for region 2, (Un^(2)=-u^(2)·n^(1))
                                ! node(sn)%col(problem%n+1+k,2): t_k for region 2, (t_k^(2)=p^(1)n_k^(1)+tau^(2)n_k^(1))
                                !
                                ! Region 2 variables (poroelastic medium)
                                ! Solid skeleton
                                do k=1,problem%n
                                  node(sn)%value_c(            k,2)=b_c(node(sn)%col(            k,2),1)
                                  node(sn)%value_c(problem%n+1+k,2)=(b_c(node(sn)%col(1,1),1)+b_c(node(sn)%col(0,2),1))*n_fn(k)
                                end do
                                ! Fluid phase
                                node(sn)%value_c(0,2)=b_c(node(sn)%col(0,2),1)
                                node(sn)%value_c(problem%n+1,2)=0.0d0
                                do k=1,problem%n
                                  node(sn)%value_c(problem%n+1,2)=node(sn)%value_c(problem%n+1,2)-node(sn)%value_c(k,2)*n_fn(k)
                                end do
                                ! Region 1 variables (fluid)
                                node(sn)%value_c(1,1)=b_c(node(sn)%col(1,1),1)
                                node(sn)%value_c(2,1)=0.0d0
                                do k=1,problem%n
                                  node(sn)%value_c(2,1)=node(sn)%value_c(2,1)+node(sn)%value_c(k,2)*n_fn(k)
                                end do
                            end select

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
                              node(sn)%value_c(k,1)=b_c(node(sn)%col(k,1),1)
                            end do
                            ! p of region 2
                            node(sn)%value_c(1,2)=b_c(node(sn)%col(1,2),1)
                            ! t_k of region 1
                            do k=1,problem%n
                              node(sn)%value_c(problem%n+k,1)=-node(sn)%value_c(1,2)*n_fn(k)
                            end do
                            ! Un of region 2
                            node(sn)%value_c(2,2)=0.0d0
                            do k=1,problem%n
                              node(sn)%value_c(2,2)=node(sn)%value_c(2,2)-node(sn)%value_c(k,1)*n_fn(k)
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
                                  node(sn)%value_c(          k,1)=b_c(node(sn)%col(          k,1),1)
                                  node(sn)%value_c(problem%n+k,1)=b_c(node(sn)%col(problem%n+k,1),1)
                                  node(sn)%value_c(          k,2)= node(sn)%value_c(          k,1)
                                  node(sn)%value_c(problem%n+k,2)=-node(sn)%value_c(problem%n+k,1)
                                ! Global axes
                                ! 1 - perfect debonding - u_k^{1} and u_k^{2} are active since:
                                ! t_k^{(1)}=0
                                ! t_k^{(2)}=0
                                case (1)
                                  node(sn)%value_c(          k,1)=b_c(node(sn)%col(k,1),1)
                                  node(sn)%value_c(problem%n+k,1)=0.0d0
                                  node(sn)%value_c(          k,2)=b_c(node(sn)%col(k,2),1)
                                  node(sn)%value_c(problem%n+k,2)=0.0d0
                                ! Global axes
                                ! 2 - partial bonding - u_k^{1} and u_k^{2} are active since:
                                ! t_k^{(2)}+t_k^{(1)}=0
                                ! (u_k^{2}-u_k^{1})*KT=t_k^{1}
                                ! Then
                                ! t_k^{(1)}=(u_k^{2}-u_k^{1})*KT
                                ! t_k^{(2)}=-(u_k^{2}-u_k^{1})*KT
                                case (2)
                                  ! KT = K + i*omega*C - omega**2*M
                                  KT=node(sn)%cvalue_c(k,1,1)+c_im*omega*node(sn)%cvalue_c(k,2,1)&
                                    -omega**2*node(sn)%cvalue_c(k,3,1)
                                  node(sn)%value_c(          k,1)=b_c(node(sn)%col(k,1),1)
                                  node(sn)%value_c(problem%n+k,1)=(b_c(node(sn)%col(k,2),1)-b_c(node(sn)%col(k,1),1))*KT
                                  node(sn)%value_c(          k,2)=b_c(node(sn)%col(k,2),1)
                                  node(sn)%value_c(problem%n+k,2)=-node(sn)%value_c(problem%n+k,1)
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
                                  node(sn)%value_c(          k,1)=b_c(node(sn)%col(k,1),1)
                                  node(sn)%value_c(problem%n+k,1)=b_c(node(sn)%col(problem%n+k,1),1)
                                  node(sn)%value_c(          k,2)=b_c(node(sn)%col(k,2),1)
                                  node(sn)%value_c(problem%n+k,2)=-node(sn)%value_c(problem%n+k,1)
                              end select
                            end do

                          ! -------------------------------------------------------
                          ! BE (1: viscoelastic solid) - BE (2: poroelastic medium)
                          ! -------------------------------------------------------

                          case (fbem_poroelastic)
                            ! node(sn)%col(            k,1): u_k for region 1
                            ! node(sn)%col(  problem%n+k,1): t_k for region 1
                            ! node(sn)%col(            0,2): tau for region 2
                            ! node(sn)%col(            k,2): u_k for region 2
                            ! node(sn)%col(problem%n+1  ,2): Un  for region 2
                            ! node(sn)%col(problem%n+1+k,2): t_k for region 2
                            !
                            ! Fluid phase variables of the poroelastic region 2
                            !
                            node(sn)%value_c(          0,2)=b_c(node(sn)%col(0,2),1)
                            ! Un^(2)=-u^(2)·n^(1)
                            node(sn)%value_c(problem%n+1,2)=0.0d0
                            do k=1,problem%n
                              node(sn)%value_c(problem%n+1,2)=node(sn)%value_c(problem%n+1,2)-b_c(node(sn)%col(k,2),1)*n_fn(k)
                            end do
                            !
                            ! Solid skeleton
                            !
                            do k=1,problem%n
                              select case (node(sn)%ctype(k,1))
                                ! Global axes
                                ! 0 - perfect bonding - u_k^{2} and t_k^{2} are active since:
                                ! t_k^{(1)}=-t_k^{(2)}+tau^(2)n_k^(1)
                                ! u_k^{(1)}= u_k^{(2)}
                                case (0)
                                  ! Region 1
                                  node(sn)%value_c(          k,1)= b_c(node(sn)%col(            k,2),1)
                                  node(sn)%value_c(problem%n+k,1)=-b_c(node(sn)%col(problem%n+1+k,2),1)+node(sn)%value_c(0,2)*n_fn(k)
                                  ! Region 2
                                  node(sn)%value_c(            k,2)=b_c(node(sn)%col(            k,2),1)
                                  node(sn)%value_c(problem%n+1+k,2)=b_c(node(sn)%col(problem%n+1+k,2),1)
                                ! Global axes
                                ! 1 - perfect debonding - u_k^{1} and u_k^{2} are active since:
                                ! t_k^{(1)}=0
                                ! t_k^{(2)}=tau^(2)n_k^(1)
                                case (1)
                                  ! Region 1
                                  node(sn)%value_c(          k,1)=b_c(node(sn)%col(k,1),1)
                                  node(sn)%value_c(problem%n+k,1)=0.0d0
                                  ! Region 2
                                  node(sn)%value_c(            k,2)=b_c(node(sn)%col(k,2),1)
                                  node(sn)%value_c(problem%n+1+k,2)=node(sn)%value_c(0,2)*n_fn(k)
                                ! Global axes
                                ! 2 - partial bonding   - u_k^{1} and u_k^{2} are active since:
                                ! t_k^{(1)}=(u_k^{2}-u_k^{1})*KT
                                ! t_k^{(2)}=-(u_k^{2}-u_k^{1})*KT+tau^(2)n_k^(1)
                                case (2)
                                  ! KT = K + i*omega*C - omega**2*M
                                  KT=node(sn)%cvalue_c(k,1,1)+c_im*omega*node(sn)%cvalue_c(k,2,1)&
                                    -omega**2*node(sn)%cvalue_c(k,3,1)
                                  ! Region 1
                                  node(sn)%value_c(          k,1)=b_c(node(sn)%col(k,1),1)
                                  node(sn)%value_c(problem%n+k,1)=(b_c(node(sn)%col(k,2),1)-b_c(node(sn)%col(k,1),1))*KT
                                  ! Region 2
                                  node(sn)%value_c(            k,2)=b_c(node(sn)%col(k,2),1)
                                  node(sn)%value_c(problem%n+1+k,2)=-node(sn)%value_c(problem%n+k,1)+node(sn)%value_c(0,2)*n_fn(k)
                              ! Local axes
                              ! In all cases, u_k^{1}, u_k^{2} and t_k^{1} are active, and t_k^{(2)}=-t_k^{(1)}+tau^(2)n_k^(1).
                                case (3,4,5)
                                  ! Region 1
                                  node(sn)%value_c(          k,1)=b_c(node(sn)%col(          k,1),1)
                                  node(sn)%value_c(problem%n+k,1)=b_c(node(sn)%col(problem%n+k,1),1)
                                  ! Region 2
                                  node(sn)%value_c(            k,2)=b_c(node(sn)%col(k,2),1)
                                  node(sn)%value_c(problem%n+1+k,2)=-node(sn)%value_c(problem%n+k,1)+node(sn)%value_c(0,2)*n_fn(k)
                              end select
                            end do

                        end select

                      case (fbem_poroelastic)
                        select case (region(boundary(sb)%region(2))%type)

                          ! ---------------------------------------------------
                          ! BE (1: poroelastic medium) - BE (2: inviscid fluid)
                          ! ---------------------------------------------------

                          case (fbem_potential)
                            !
                            ! Switch depending on the B.C.
                            !
                            select case (node(sn)%ctype(1,1))
                              !
                              ! Perfectly permeable
                              !
                              case (0)
                                ! Index of variables:
                                ! node(sn)%col(            0,1): tau for region 1
                                ! node(sn)%col(            k,1): u_k for region 1
                                ! node(sn)%col(problem%n+1  ,1): Un  for region 1
                                ! node(sn)%col(problem%n+1+k,1): t_k for region 1
                                ! node(sn)%col(            1,2): p   for region 2, (p^(2)=-tau^(1)/phi^(1))
                                ! node(sn)%col(            2,2): Un  for region 2, (Un^(2)=-phi^(1)*Un^(1)-(1-phi^(1))*u^(1)·n^(1))
                                !
                                ! Region 1 variables (poroelastic medium)
                                phi=region(boundary(sb)%region(1))%property_r(8)
                                do k=0,problem%n
                                  node(sn)%value_c(            k,1)=b_c(node(sn)%col(            k,1),1)
                                  if (k.eq.0) then
                                    node(sn)%value_c(problem%n+1+k,1)=b_c(node(sn)%col(problem%n+1+k,1),1)
                                  else
                                    node(sn)%value_c(problem%n+1+k,1)=(1.0d0-phi)/phi*node(sn)%value_c(0,1)*n_fn(k)
                                  end if
                                end do
                                ! Region 2 variables (fluid)
                                node(sn)%value_c(1,2)=-node(sn)%value_c(0,1)/phi
                                node(sn)%value_c(2,2)=-phi*node(sn)%value_c(problem%n+1,1)
                                do k=1,problem%n
                                  node(sn)%value_c(2,2)=node(sn)%value_c(2,2)-(1.0d0-phi)*node(sn)%value_c(k,1)*n_fn(k)
                                end do
                              !
                              ! Perfectly impermeable
                              !
                              case (1)
                                ! Index of variables:
                                ! node(sn)%col(            0,1): tau for region 1
                                ! node(sn)%col(            k,1): u_k for region 1
                                ! node(sn)%col(problem%n+1  ,1): Un  for region 1, (Un^(1)=u^(1)·n^(1))
                                ! node(sn)%col(problem%n+1+k,1): t_k for region 1
                                ! node(sn)%col(            1,2): p   for region 2
                                ! node(sn)%col(            2,2): Un  for region 2, (Un^(2)=-u^(1)·n^(1))
                                !
                                ! Region 1 variables (poroelastic medium)
                                ! Solid skeleton
                                do k=1,problem%n
                                  node(sn)%value_c(            k,1)=b_c(node(sn)%col(            k,1),1)
                                  node(sn)%value_c(problem%n+1+k,1)=-(b_c(node(sn)%col(1,2),1)+b_c(node(sn)%col(0,1),1))*n_fn(k)
                                end do
                                ! Fluid phase
                                node(sn)%value_c(0,1)=b_c(node(sn)%col(0,1),1)
                                node(sn)%value_c(problem%n+1,1)=0.0d0
                                do k=1,problem%n
                                  node(sn)%value_c(problem%n+1,1)=node(sn)%value_c(problem%n+1,1)+node(sn)%value_c(k,1)*n_fn(k)
                                end do
                                ! Region 2 variables (fluid)
                                node(sn)%value_c(1,2)=b_c(node(sn)%col(1,2),1)
                                node(sn)%value_c(2,2)=0.0d0
                                do k=1,problem%n
                                  node(sn)%value_c(2,2)=node(sn)%value_c(2,2)-node(sn)%value_c(k,1)*n_fn(k)
                                end do
                            end select

                          ! -------------------------------------------------------
                          ! BE (1: poroelastic medium) - BE (2: viscoelastic solid)
                          ! -------------------------------------------------------

                          case (fbem_viscoelastic)
                            ! node(sn)%col(            0,1): tau for region 1
                            ! node(sn)%col(            k,1): u_k for region 1
                            ! node(sn)%col(problem%n+1  ,1): Un  for region 1
                            ! node(sn)%col(problem%n+1+k,1): t_k for region 1
                            ! node(sn)%col(            k,2): u_k for region 2
                            ! node(sn)%col(  problem%n+k,2): t_k for region 2
                            !
                            ! Fluid phase variables of the poroelastic region 1
                            !
                            node(sn)%value_c(          0,1)=b_c(node(sn)%col(0,1),1)
                            ! Un^(1)=u^(1)·n^(1)
                            node(sn)%value_c(problem%n+1,1)=0.0d0
                            do k=1,problem%n
                              node(sn)%value_c(problem%n+1,1)=node(sn)%value_c(problem%n+1,1)+b_c(node(sn)%col(k,1),1)*n_fn(k)
                            end do
                            !
                            ! Solid skeleton
                            !
                            do k=1,problem%n
                              select case (node(sn)%ctype(k,1))
                                ! Global axes
                                ! 0 - perfect bonding - u_k^{1} and t_k^{1} are active since:
                                ! t_k^{(2)}=-t_k^{(1)}-tau^(1)n_k^(1)
                                ! u_k^{(2)}= u_k^{(1)}
                                case (0)
                                  ! Region 1
                                  node(sn)%value_c(            k,1)=b_c(node(sn)%col(            k,1),1)
                                  node(sn)%value_c(problem%n+1+k,1)=b_c(node(sn)%col(problem%n+1+k,1),1)
                                  ! Region 2
                                  node(sn)%value_c(          k,2)= b_c(node(sn)%col(            k,1),1)
                                  node(sn)%value_c(problem%n+k,2)=-b_c(node(sn)%col(problem%n+1+k,1),1)-node(sn)%value_c(0,1)*n_fn(k)
                                ! Global axes
                                ! 1 - perfect debonding - u_k^{1} and u_k^{2} are active since:
                                ! t_k^{(1)}=-tau^(1)n_k^(1)
                                ! t_k^{(2)}=0
                                case (1)
                                  ! Region 1
                                  node(sn)%value_c(            k,1)=b_c(node(sn)%col(k,1),1)
                                  node(sn)%value_c(problem%n+1+k,1)=-node(sn)%value_c(0,1)*n_fn(k)
                                  ! Region 2
                                  node(sn)%value_c(          k,2)=b_c(node(sn)%col(k,2),1)
                                  node(sn)%value_c(problem%n+k,2)=0.0d0
                                ! Global axes
                                ! 2 - partial bonding   - u_k^{1} and u_k^{2} are active since:
                                ! t_k^{(1)}=-(u_k^{1}-u_k^{2})*KT-tau^(1)n_k^(1)
                                ! t_k^{(2)}=(u_k^{1}-u_k^{2})*KT
                                case (2)
                                  ! KT = K + i*omega*C - omega**2*M
                                  KT=node(sn)%cvalue_c(k,1,1)+c_im*omega*node(sn)%cvalue_c(k,2,1)&
                                    -omega**2*node(sn)%cvalue_c(k,3,1)
                                  ! Region 2
                                  node(sn)%value_c(          k,2)=b_c(node(sn)%col(k,2),1)
                                  node(sn)%value_c(problem%n+k,2)=(b_c(node(sn)%col(k,1),1)-b_c(node(sn)%col(k,2),1))*KT
                                  ! Region 1
                                  node(sn)%value_c(            k,1)=b_c(node(sn)%col(k,1),1)
                                  node(sn)%value_c(problem%n+1+k,1)=-node(sn)%value_c(problem%n+k,2)-node(sn)%value_c(0,1)*n_fn(k)
                                ! Local axes
                                ! In all cases, u_k^{1}, u_k^{2} and t_k^{2} are active, and t_k^{(1)}=-t_k^{(2)}-tau^(1)n_k^(1).
                                case (3,4,5)
                                  ! Region 2
                                  node(sn)%value_c(          k,2)=b_c(node(sn)%col(          k,2),1)
                                  node(sn)%value_c(problem%n+k,2)=b_c(node(sn)%col(problem%n+k,2),1)
                                  ! Region 1
                                  node(sn)%value_c(            k,1)=b_c(node(sn)%col(k,1),1)
                                  node(sn)%value_c(problem%n+1+k,1)=-node(sn)%value_c(problem%n+k,2)-node(sn)%value_c(0,1)*n_fn(k)
                              end select
                            end do

                          ! -------------------------------------------------------
                          ! BE (1: poroelastic medium) - BE (2: poroelastic medium)
                          ! -------------------------------------------------------

                          case (fbem_poroelastic)
                            !
                            ! Switch depending on the B.C.
                            !
                            select case (node(sn)%ctype(1,1))
                              !
                              ! Perfectly permeable
                              !
                              case (0)
                                ! Index of variables:
                                ! node(sn)%col(            0,1): tau for region 1
                                ! node(sn)%col(            k,1): u_k for region 1
                                ! node(sn)%col(problem%n+1  ,1): Un  for region 1
                                ! node(sn)%col(problem%n+1+k,1): t_k for region 1
                                ! node(sn)%col(            0,2): tau for region 2 (not active since tau^(2)=phi^(2)/phi^(1)*tau^(1))
                                ! node(sn)%col(            k,2): u_k for region 2 (not active since u_k^(2)=u_k^(1))
                                ! node(sn)%col(problem%n+1  ,2): Un  for region 2 (not active since Un^(2)=-phi^(1)/phi^(2)*Un^(1)-(1-phi^(1)/phi^(2))*u^(1)·n^(1))
                                ! node(sn)%col(problem%n+1+k,2): t_k for region 2 (not active since t_k^(2)=-t_k^(1)-(1-phi^(2)/phi^(1))*n_k^(1)*tau^(1))
                                !
                                ! Region 1 variables
                                do k=0,problem%n
                                  node(sn)%value_c(            k,1)=b_c(node(sn)%col(            k,1),1)
                                  node(sn)%value_c(problem%n+1+k,1)=b_c(node(sn)%col(problem%n+1+k,1),1)
                                end do
                                ! Region 2 variables
                                ! Auxiliary constants
                                d_phi1_phi2=(region(boundary(sb)%region(1))%property_r(8))&
                                          /(region(boundary(sb)%region(2))%property_r(8))
                                d_phi2_phi1=1.0d0/d_phi1_phi2
                                p_1_phi2_phi1=1.0d0-d_phi2_phi1
                                p_1_phi1_phi2=1.0d0-d_phi1_phi2
                                ! Fluid phase
                                node(sn)%value_c(          0,2)= d_phi2_phi1*node(sn)%value_c(0,1)
                                node(sn)%value_c(problem%n+1,2)=-d_phi1_phi2*node(sn)%value_c(problem%n+1,1)
                                do k=1,problem%n
                                  node(sn)%value_c(problem%n+1,2)=node(sn)%value_c(problem%n+1,2)&
                                                                 -p_1_phi1_phi2*node(sn)%value_c(k,1)*n_fn(k)
                                end do
                                ! Solid skeleton
                                do k=1,problem%n
                                  node(sn)%value_c(            k,2)= node(sn)%value_c(            k,1)
                                  node(sn)%value_c(problem%n+1+k,2)=-node(sn)%value_c(problem%n+1+k,1)&
                                                                   -p_1_phi2_phi1*node(sn)%value_c(0,1)*n_fn(k)
                                end do
                              !
                              ! Perfectly impermeable
                              !
                              case (1)
                                ! Index of variables:
                                ! node(sn)%col(            0,1): tau for region 1 (active)
                                ! node(sn)%col(            k,1): u_k for region 1 (active)
                                ! node(sn)%col(problem%n+1  ,1): Un  for region 1 (not active since Un^(1)=u^(1)·n^(1))
                                ! node(sn)%col(problem%n+1+k,1): t_k for region 1 (active)
                                ! node(sn)%col(            0,2): tau for region 2 (active)
                                ! node(sn)%col(            k,2): u_k for region 2 (not active since u_k^(2)=u_k^(1))
                                ! node(sn)%col(problem%n+1  ,2): Un  for region 2 (not active since Un^(2)=-u^(1)·n^(1))
                                ! node(sn)%col(problem%n+1+k,2): t_k for region 2 (not active since t_k^(2)=-t_k^(1)-tau^(1)n_k^(1)+tau^(2)n_k^(1))
                                !
                                ! Region 1 variables
                                ! Solid skeleton
                                do k=1,problem%n
                                  node(sn)%value_c(            k,1)=b_c(node(sn)%col(            k,1),1)
                                  node(sn)%value_c(problem%n+1+k,1)=b_c(node(sn)%col(problem%n+1+k,1),1)
                                end do
                                ! Fluid phase
                                node(sn)%value_c(          0,1)=b_c(node(sn)%col(0,1),1)
                                node(sn)%value_c(problem%n+1,1)=0.0d0
                                do k=1,problem%n
                                  node(sn)%value_c(problem%n+1,1)=node(sn)%value_c(problem%n+1,1)+node(sn)%value_c(k,1)*n_fn(k)
                                end do
                                ! Region 2 variables
                                ! Fluid phase
                                node(sn)%value_c(          0,2)=b_c(node(sn)%col(0,2),1)
                                node(sn)%value_c(problem%n+1,2)=0.0d0
                                do k=1,problem%n
                                  node(sn)%value_c(problem%n+1,2)=node(sn)%value_c(problem%n+1,2)-node(sn)%value_c(k,1)*n_fn(k)
                                end do
                                ! Solid skeleton
                                do k=1,problem%n
                                  node(sn)%value_c(            k,2)= node(sn)%value_c(            k,1)
                                  node(sn)%value_c(problem%n+1+k,2)=-node(sn)%value_c(problem%n+1+k,1)&
                                                                   -node(sn)%value_c(0,1)*n_fn(k)+node(sn)%value_c(0,2)*n_fn(k)
                                end do
                              !
                              ! Partially permeable
                              !
                              case (2)
                                ! Index of variables:
                                ! node(sn)%col(            0,1): tau for region 1 (active)
                                ! node(sn)%col(            k,1): u_k for region 1 (active)
                                ! node(sn)%col(problem%n+1  ,1): Un  for region 1 (active)
                                ! node(sn)%col(problem%n+1+k,1): t_k for region 1 (active)
                                ! node(sn)%col(            0,2): tau for region 2 (not active since tau^(2)=k*phi1*phi2*i*omega*(Un^(1)-u_j^(1)*n_j^(1))+phi^(2)/phi^(1)*tau^(1))
                                ! node(sn)%col(            k,2): u_k for region 2 (not active since u_k^(2)=u_k^(1))
                                ! node(sn)%col(problem%n+1  ,2): Un  for region 2 (not active since Un^(2)=-phi^(1)/phi^(2)*Un^(1)-(1-phi^(1)/phi^(2))*u^(1)·n^(1))
                                ! node(sn)%col(problem%n+1+k,2): t_k for region 2 (not active since t_k^(2)=-t_k^(1)+[k*phi1*phi2*i*omega*(Un^(1)-u_j^(1)*n_j^(1))-(1-phi^(2)/phi^(1))*tau^(1)]*n_k^(1)
                                !
                                ! Region 1 variables
                                do k=0,problem%n
                                  node(sn)%value_c(            k,1)=b_c(node(sn)%col(            k,1),1)
                                  node(sn)%value_c(problem%n+1+k,1)=b_c(node(sn)%col(problem%n+1+k,1),1)
                                end do
                                ! Region 2 variables
                                ! Auxiliary constants
                                phi1=region(boundary(sb)%region(1))%property_r(8)
                                phi2=region(boundary(sb)%region(2))%property_r(8)
                                d_phi1_phi2=phi1/phi2
                                d_phi2_phi1=1.0d0/d_phi1_phi2
                                p_1_phi2_phi1=1.0d0-d_phi2_phi1
                                p_1_phi1_phi2=1.0d0-d_phi1_phi2
                                KT=node(sn)%cvalue_c(1,1,1)*phi1*phi2*c_im*omega
                                ! Fluid phase
                                ! tau^(2)
                                node(sn)%value_c(0,2)=d_phi2_phi1*node(sn)%value_c(0,1)
                                node(sn)%value_c(0,2)=node(sn)%value_c(0,2)+KT*node(sn)%value_c(problem%n+1,1)
                                do k=1,problem%n
                                  node(sn)%value_c(0,2)=node(sn)%value_c(0,2)-KT*node(sn)%value_c(k,1)*n_fn(k)
                                end do
                                ! Un^(2)
                                node(sn)%value_c(problem%n+1,2)=-d_phi1_phi2*node(sn)%value_c(problem%n+1,1)
                                do k=1,problem%n
                                  node(sn)%value_c(problem%n+1,2)=node(sn)%value_c(problem%n+1,2)-p_1_phi1_phi2*node(sn)%value_c(k,1)*n_fn(k)
                                end do
                                ! Solid skeleton
                                do k=1,problem%n
                                  node(sn)%value_c(            k,2)= node(sn)%value_c(            k,1)
                                  node(sn)%value_c(problem%n+1+k,2)=-node(sn)%value_c(problem%n+1+k,1)
                                  node(sn)%value_c(problem%n+1+k,2)= node(sn)%value_c(problem%n+1+k,2)-p_1_phi2_phi1*node(sn)%value_c(0,1)*n_fn(k)
                                  node(sn)%value_c(problem%n+1+k,2)= node(sn)%value_c(problem%n+1+k,2)+KT*node(sn)%value_c(problem%n+1,1)*n_fn(k)
                                  do kc=1,problem%n
                                    node(sn)%value_c(problem%n+1+k,2)= node(sn)%value_c(problem%n+1+k,2)-KT*node(sn)%value_c(kc,1)*n_fn(kc)*n_fn(k)
                                  end do
                                end do
                            end select

!                            ! Change of variables for traction quarter point elements
!                            if (element(se)%type.eq.fbem_line3) then
!                              select case (kn)
!                                ! Node 1
!                                case (1)
!                                  select case (element(se)%type_f1f)
!                                    case (fbem_line3_qp1t,fbem_line3_mqp1t)
!                                      node(sn)%value_c(0,1)=huge(1.d0)
!                                      node(sn)%value_c(0,2)=huge(1.d0)
!                                  end select
!                                  select case (element(se)%type_f1)
!                                    case (fbem_line3_qp1t,fbem_line3_mqp1t)
!                                      do k=1,problem%n
!                                        node(sn)%value_c(k,1)=huge(1.d0)
!                                        node(sn)%value_c(k,2)=huge(1.d0)
!                                      end do
!                                  end select
!                                  select case (element(se)%type_f2f)
!                                    case (fbem_line3_qp1t,fbem_line3_mqp1t)
!                                      node(sn)%value_c(problem%n+1,1)=huge(1.d0)
!                                      node(sn)%value_c(problem%n+1,2)=huge(1.d0)
!                                  end select
!                                  select case (element(se)%type_f2)
!                                    case (fbem_line3_qp1t,fbem_line3_mqp1t)
!                                      do k=1,problem%n
!                                        node(sn)%value_c(problem%n+1+k,1)=huge(1.d0)
!                                        node(sn)%value_c(problem%n+1+k,2)=huge(1.d0)
!                                      end do
!                                  end select
!                                ! Node 2
!                                case (2)
!                                  select case (element(se)%type_f1f)
!                                    case (fbem_line3_qp2t,fbem_line3_mqp2t)
!                                      node(sn)%value_c(0,1)=huge(1.d0)
!                                      node(sn)%value_c(0,2)=huge(1.d0)
!                                  end select
!                                  select case (element(se)%type_f1)
!                                    case (fbem_line3_qp2t,fbem_line3_mqp2t)
!                                      do k=1,problem%n
!                                        node(sn)%value_c(k,1)=huge(1.d0)
!                                        node(sn)%value_c(k,2)=huge(1.d0)
!                                      end do
!                                  end select
!                                  select case (element(se)%type_f2f)
!                                    case (fbem_line3_qp2t,fbem_line3_mqp2t)
!                                      node(sn)%value_c(problem%n+1,1)=huge(1.d0)
!                                      node(sn)%value_c(problem%n+1,2)=huge(1.d0)
!                                  end select
!                                  select case (element(se)%type_f2)
!                                    case (fbem_line3_qp2t,fbem_line3_mqp2t)
!                                      do k=1,problem%n
!                                        node(sn)%value_c(problem%n+1+k,1)=huge(1.d0)
!                                        node(sn)%value_c(problem%n+1+k,2)=huge(1.d0)
!                                      end do
!                                  end select
!                                ! Node 3
!                                case (3)
!                                  select case (element(se)%type_f1f)
!                                    case (fbem_line3_qp1t,fbem_line3_qp2t,fbem_line3_mqp1t,fbem_line3_mqp2t)
!                                      node(sn)%value_c(0,1)=2.d0*node(sn)%value_c(0,1)
!                                      node(sn)%value_c(0,2)=2.d0*node(sn)%value_c(0,2)
!                                  end select
!                                  select case (element(se)%type_f1)
!                                    case (fbem_line3_qp1t,fbem_line3_qp2t,fbem_line3_mqp1t,fbem_line3_mqp2t)
!                                      do k=1,problem%n
!                                        node(sn)%value_c(k,1)=2.d0*node(sn)%value_c(k,1)
!                                        node(sn)%value_c(k,2)=2.d0*node(sn)%value_c(k,2)
!                                      end do
!                                  end select
!                                  select case (element(se)%type_f2f)
!                                    case (fbem_line3_qp1t,fbem_line3_qp2t,fbem_line3_mqp1t,fbem_line3_mqp2t)
!                                      node(sn)%value_c(problem%n+1,1)=2.d0*node(sn)%value_c(problem%n+1,1)
!                                      node(sn)%value_c(problem%n+1,2)=2.d0*node(sn)%value_c(problem%n+1,2)
!                                  end select
!                                  select case (element(se)%type_f2)
!                                    case (fbem_line3_qp1t,fbem_line3_qp2t,fbem_line3_mqp1t,fbem_line3_mqp2t)
!                                      do k=1,problem%n
!                                        node(sn)%value_c(problem%n+1+k,1)=2.d0*node(sn)%value_c(problem%n+1+k,1)
!                                        node(sn)%value_c(problem%n+1+k,2)=2.d0*node(sn)%value_c(problem%n+1+k,2)
!                                      end do
!                                  end select
!                              end select
!                            end if

                        end select
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
                            node(sn)%value_c(1,1)=b_c(node(sn)%col(1,1),1)
                            ! Un=u^(fe)·n
                            ! Un initialization
                            node(sn)%value_c(2,1)=0.0d0
                            ! If the boundary is not reversed
                            if (.not.sb_reversion) then
                              do k=1,problem%n
                                if (node(sn_fe)%ctype(k,1).eq.1) then
                                  node(sn)%value_c(2,1)=node(sn)%value_c(2,1)+b_c(node(sn_fe)%col(k,1),1)*n_fn(k)
                                else
                                  node(sn)%value_c(2,1)=node(sn)%value_c(2,1)+node(sn_fe)%cvalue_c(k,1,1)*n_fn(k)
                                end if
                              end do
                            ! If the boundary is reversed
                            else
                              do k=1,problem%n
                                if (node(sn_fe)%ctype(k,1).eq.1) then
                                  node(sn)%value_c(2,1)=node(sn)%value_c(2,1)-b_c(node(sn_fe)%col(k,1),1)*n_fn(k)
                                else
                                  node(sn)%value_c(2,1)=node(sn)%value_c(2,1)-node(sn_fe)%cvalue_c(k,1,1)*n_fn(k)
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
                                node(sn)%value_c(k,1)=b_c(node(sn_fe)%col(k,1),1)
                              else
                                node(sn)%value_c(k,1)=node(sn_fe)%cvalue_c(k,1,1)
                              end if
                              node(sn)%value_c(problem%n+k,1)=b_c(node(sn)%col(problem%n+k,1),1)
                            end do

                          ! ------------------
                          ! POROELASTIC MEDIUM
                          ! ------------------

                          case (fbem_poroelastic)
                            ! Index of variables:
                            ! node(sn)%col(            0,1): tau
                            ! node(sn)%col(            k,1): u_k
                            ! node(sn)%col(problem%n+1  ,1): Un (inactive since Un=u·n, impermeable condition)
                            ! node(sn)%col(problem%n+1+k,1): t_k
                            ! tau
                            node(sn)%value_c(0,1)=b_c(node(sn)%col(0,1),1)
                            ! Un=u·n
                            ! Un initialization
                            node(sn)%value_c(problem%n+1,1)=0.0d0
                            ! If the boundary is not reversed
                            if (.not.sb_reversion) then
                              do k=1,problem%n
                                if (node(sn_fe)%ctype(k,1).eq.1) then
                                  node(sn)%value_c(problem%n+1,1)=node(sn)%value_c(problem%n+1,1)+b_c(node(sn_fe)%col(k,1),1)*n_fn(k)
                                else
                                  node(sn)%value_c(problem%n+1,1)=node(sn)%value_c(problem%n+1,1)+node(sn_fe)%cvalue_c(k,1,1)*n_fn(k)
                                end if
                              end do
                            ! If the boundary is reversed
                            else
                              do k=1,problem%n
                                if (node(sn_fe)%ctype(k,1).eq.1) then
                                  node(sn)%value_c(problem%n+1,1)=node(sn)%value_c(problem%n+1,1)-b_c(node(sn_fe)%col(k,1),1)*n_fn(k)
                                else
                                  node(sn)%value_c(problem%n+1,1)=node(sn)%value_c(problem%n+1,1)-node(sn_fe)%cvalue_c(k,1,1)*n_fn(k)
                                end if
                              end do
                            end if
                            ! Solid skeleton
                            do k=1,problem%n
                              if (node(sn_fe)%ctype(k,1).eq.1) then
                                node(sn)%value_c(k,1)=b_c(node(sn_fe)%col(k,1),1)
                              else
                                node(sn)%value_c(k,1)=node(sn_fe)%cvalue_c(k,1,1)
                              end if
                              node(sn)%value_c(problem%n+1+k,1)=b_c(node(sn)%col(problem%n+1+k,1),1)
                            end do

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
                            node(sn)%value_c(1,1)=b_c(node(sn)%col(1,1),1)
                            ! p^-
                            node(sn)%value_c(1,2)=b_c(node(sn)%col(1,2),1)
                            ! Un^+= u^(fe)·n
                            ! The fe node connected with the be node
                            sn_fe=element(se)%element_node(kn)
                            ! Un^+ initialization
                            node(sn)%value_c(2,1)=0.0d0
                            ! If the boundary is not reversed
                            if (.not.sb_reversion) then
                              do k=1,problem%n
                                if (node(sn_fe)%ctype(k,1).eq.1) then
                                  node(sn)%value_c(2,1)=node(sn)%value_c(2,1)+b_c(node(sn_fe)%col(k,1),1)*n_fn(k)
                                else
                                  node(sn)%value_c(2,1)=node(sn)%value_c(2,1)+node(sn_fe)%cvalue_c(k,1,1)*n_fn(k)
                                end if
                              end do
                            ! If the boundary is reversed
                            else
                              do k=1,problem%n
                                if (node(sn_fe)%ctype(k,1).eq.1) then
                                  node(sn)%value_c(2,1)=node(sn)%value_c(2,1)-b_c(node(sn_fe)%col(k,1),1)*n_fn(k)
                                else
                                  node(sn)%value_c(2,1)=node(sn)%value_c(2,1)-node(sn_fe)%cvalue_c(k,1,1)*n_fn(k)
                                end if
                              end do
                            end if
                            ! Un^-=-Un^+
                            node(sn)%value_c(2,2)=-node(sn)%value_c(2,1)

                          ! ------------------
                          ! VISCOELASTIC SOLID
                          ! ------------------

                          case (fbem_viscoelastic)
!                            !
!                            ! SBIE (solve for Σt = t^+ + t^-)
!                            !
!                            ! Index of variable values for each coordinate k:
!                            ! node(sn)%value_r(          k,1): u_k^+
!                            ! node(sn)%value_r(problem%n+k,1): t_k^+ (temporarily before calculating t_k^+ and t_k^- it will be Σt)
!                            ! node(sn)%value_r(          k,2): u_k^-
!                            ! node(sn)%value_r(problem%n+k,2): t_k^- (temporarily before calculating t_k^+ and t_k^- it will be 0)
!                            do k=1,problem%n
!                              if (node(sn_fe)%ctype(k,1).eq.1) then
!                                node(sn)%value_c(k,1)=b_c(node(sn_fe)%col(k,1),1)
!                              else
!                                node(sn)%value_c(k,1)=node(sn_fe)%cvalue_c(k,1,1)
!                              end if
!                              node(sn)%value_c(problem%n+k,1)=b_c(node(sn)%col(problem%n+k,1),1)
!                              node(sn)%value_c(          k,2)=node(sn)%value_c(k,1)
!                              node(sn)%value_c(problem%n+k,2)=0.d0
!                            end do
                            !
                            ! Usando Dual BEM (para CC mas generales)
                            !
                            ! Index of variables for each coordinate k:
                            ! node(sn)%col(          k,1): u_k for face +
                            ! node(sn)%col(problem%n+k,1): t_k for face +
                            ! node(sn)%col(          k,2): u_k for face - (not used since u_l_{+} = u_k_{-})
                            ! node(sn)%col(problem%n+k,2): t_k for face -
                            do k=1,problem%n
                              if (node(sn_fe)%ctype(k,1).eq.1) then
                                node(sn)%value_c(k,1)=b_c(node(sn_fe)%col(k,1),1)
                              else
                                node(sn)%value_c(k,1)=node(sn_fe)%cvalue_c(k,1,1)
                              end if
                              node(sn)%value_c(problem%n+k,1)=b_c(node(sn)%col(problem%n+k,1),1)
                              node(sn)%value_c(          k,2)=node(sn)%value_c(k,1)
                              node(sn)%value_c(problem%n+k,2)=b_c(node(sn)%col(problem%n+k,2),1)
                            end do

                          ! ------------------
                          ! POROELASTIC MEDIUM
                          ! ------------------

                          case (fbem_poroelastic)
                            ! Index of variables:
                            ! node(sn)%col(            0,1): tau^+
                            ! node(sn)%col(            k,1): u_k^+
                            ! node(sn)%col(problem%n+1  ,1): Un^+ (inactive since Un^+=u^+·n^+, impermeable condition)
                            ! node(sn)%col(problem%n+1+k,1): t_k^+
                            ! node(sn)%col(            0,2): tau^-
                            ! node(sn)%col(            k,2): u_k^-
                            ! node(sn)%col(problem%n+1  ,2): Un^- (inactive since Un^-=u^-·n^-, impermeable condition)
                            ! node(sn)%col(problem%n+1+k,2): t_k^-
                            !
                            ! Fluid phase
                            ! tau^+
                            node(sn)%value_c(0,1)=b_c(node(sn)%col(0,1),1)
                            ! tau^-
                            node(sn)%value_c(0,2)=b_c(node(sn)%col(0,2),1)
                            ! Un^+=u^+·n^+
                            ! Un^+ initialization
                            node(sn)%value_c(problem%n+1,1)=0.0d0
                            ! If the boundary is not reversed
                            if (.not.sb_reversion) then
                              do k=1,problem%n
                                if (node(sn_fe)%ctype(k,1).eq.1) then
                                  node(sn)%value_c(problem%n+1,1)=node(sn)%value_c(problem%n+1,1)+b_c(node(sn_fe)%col(k,1),1)*n_fn(k)
                                else
                                  node(sn)%value_c(problem%n+1,1)=node(sn)%value_c(problem%n+1,1)+node(sn_fe)%cvalue_c(k,1,1)*n_fn(k)
                                end if
                              end do
                            ! If the boundary is reversed
                            else
                              do k=1,problem%n
                                if (node(sn_fe)%ctype(k,1).eq.1) then
                                  node(sn)%value_c(problem%n+1,1)=node(sn)%value_c(problem%n+1,1)-b_c(node(sn_fe)%col(k,1),1)*n_fn(k)
                                else
                                  node(sn)%value_c(problem%n+1,1)=node(sn)%value_c(problem%n+1,1)-node(sn_fe)%cvalue_c(k,1,1)*n_fn(k)
                                end if
                              end do
                            end if
                            ! Un^-=-Un^+
                            node(sn)%value_c(problem%n+1,2)=-node(sn)%value_c(problem%n+1,1)
                            ! Solid skeleton
                            do k=1,problem%n
                              ! u_k^+ and u_k^-
                              if (node(sn_fe)%ctype(k,1).eq.1) then
                                node(sn)%value_c(k,1)=b_c(node(sn_fe)%col(k,1),1)
                                node(sn)%value_c(k,2)=b_c(node(sn_fe)%col(k,1),1)
                              else
                                node(sn)%value_c(k,1)=node(sn_fe)%cvalue_c(k,1,1)
                                node(sn)%value_c(k,2)=node(sn_fe)%cvalue_c(k,1,1)
                              end if
                              ! t_k^+ and t_k^-
                              node(sn)%value_c(problem%n+1+k,1)=b_c(node(sn)%col(problem%n+1+k,1),1)
                              node(sn)%value_c(problem%n+1+k,2)=b_c(node(sn)%col(problem%n+1+k,2),1)
                            end do

                        end select

                    end select

                  ! =================
                  ! BE-FE-BE BOUNDARY
                  ! =================

                  case (fbem_boundary_coupling_be_fe_be)
                    ! The fe node connected with the be node
                    sn_fe=element(se)%element_node(kn)

                    ! -------------------------------------------
                    ! BOUNDARY SEEN FROM THE REGION 1 (n^(1) = n)
                    ! -------------------------------------------

                    select case (region(boundary(sb)%region(1))%type)

                      ! --------------
                      ! INVISCID FLUID
                      ! --------------

                      case (fbem_potential)
                        ! p^(1)
                        node(sn)%value_c(1,1)=b_c(node(sn)%col(1,1),1)
                        ! Un^(1) = u_k^(fe)·n_k
                        node(sn)%value_c(2,1)=0.0d0
                        do k=1,problem%n
                          if (node(sn_fe)%ctype(k,1).eq.1) then
                            node(sn)%value_c(2,1)=node(sn)%value_c(2,1)+b_c(node(sn_fe)%col(k,1),1)*n_fn(k)
                          else
                            node(sn)%value_c(2,1)=node(sn)%value_c(2,1)+node(sn_fe)%cvalue_c(k,1,1)*n_fn(k)
                          end if
                        end do

                      ! ------------------
                      ! VISCOELASTIC SOLID
                      ! ------------------

                      case (fbem_viscoelastic)
                        do k=1,problem%n
                          ! u_k^(1) = u_k^(FE)
                          if (node(sn_fe)%ctype(k,1).eq.1) then
                            node(sn)%value_c(k,1)=b_c(node(sn_fe)%col(k,1),1)
                          else
                            node(sn)%value_c(k,1)=node(sn_fe)%cvalue_c(k,1,1)
                          end if
                          ! t_k^(1)
                          node(sn)%value_c(problem%n+k,1)=b_c(node(sn)%col(problem%n+k,1),1)
                        end do

                      ! ------------------
                      ! POROELASTIC MEDIUM
                      ! ------------------

                      case (fbem_poroelastic)
                        ! Fluid phase
                        ! tau^(1)
                        node(sn)%value_c(0,1)=b_c(node(sn)%col(0,1),1)
                        ! Un^(1) = u_k^(FE)·n_k
                        node(sn)%value_c(problem%n+1,1)=0.d0
                        do k=1,problem%n
                          if (node(sn_fe)%ctype(k,1).eq.1) then
                            node(sn)%value_c(problem%n+1,1)=node(sn)%value_c(problem%n+1,1)+b_c(node(sn_fe)%col(k,1),1)*n_fn(k)
                          else
                            node(sn)%value_c(problem%n+1,1)=node(sn)%value_c(problem%n+1,1)+node(sn_fe)%cvalue_c(k,1,1)*n_fn(k)
                          end if
                        end do
                        ! Solid phase
                        do k=1,problem%n
                          ! u_k^(1) = u_k^(FE)
                          if (node(sn_fe)%ctype(k,1).eq.1) then
                            node(sn)%value_c(k,1)=b_c(node(sn_fe)%col(k,1),1)
                          else
                            node(sn)%value_c(k,1)=node(sn_fe)%cvalue_c(k,1,1)
                          end if
                          ! t_k^(1)
                          node(sn)%value_c(problem%n+1+k,1)=b_c(node(sn)%col(problem%n+1+k,1),1)
                        end do

                    end select

                    ! --------------------------------------------
                    ! BOUNDARY SEEN FROM THE REGION 2 (n^(2) = -n)
                    ! --------------------------------------------

                    select case (region(boundary(sb)%region(2))%type)

                      ! --------------
                      ! INVISCID FLUID
                      ! --------------

                      case (fbem_potential)
                        ! p^(2)
                        node(sn)%value_c(1,2)=b_c(node(sn)%col(1,2),1)
                        ! Un^(2)= -u_k^(fe)·n_k
                        ! Un initialization
                        node(sn)%value_c(2,2)=0.0d0
                        do k=1,problem%n
                          if (node(sn_fe)%ctype(k,1).eq.1) then
                            node(sn)%value_c(2,2)=node(sn)%value_c(2,2)-b_c(node(sn_fe)%col(k,1),1)*n_fn(k)
                          else
                            node(sn)%value_c(2,2)=node(sn)%value_c(2,2)-node(sn_fe)%cvalue_c(k,1,1)*n_fn(k)
                          end if
                        end do

                      ! ------------------
                      ! VISCOELASTIC SOLID
                      ! ------------------

                      case (fbem_viscoelastic)
                        do k=1,problem%n
                          ! u_k^(2) = u_k^(FE)
                          if (node(sn_fe)%ctype(k,1).eq.1) then
                            node(sn)%value_c(k,2)=b_c(node(sn_fe)%col(k,1),1)
                          else
                            node(sn)%value_c(k,2)=node(sn_fe)%cvalue_c(k,1,1)
                          end if
                          ! t_k^(2)
                          node(sn)%value_c(problem%n+k,2)=b_c(node(sn)%col(problem%n+k,2),1)
                        end do

                      ! ------------------
                      ! POROELASTIC MEDIUM
                      ! ------------------

                      case (fbem_poroelastic)
                        ! Fluid phase
                        ! tau^(2)
                        node(sn)%value_c(0,2)=b_c(node(sn)%col(0,2),1)
                        ! Un^(2) = -u_k^(FE)·n_k
                        node(sn)%value_c(problem%n+1,2)=0.d0
                        do k=1,problem%n
                          if (node(sn_fe)%ctype(k,1).eq.1) then
                            node(sn)%value_c(problem%n+1,2)=node(sn)%value_c(problem%n+1,2)-b_c(node(sn_fe)%col(k,1),1)*n_fn(k)
                          else
                            node(sn)%value_c(problem%n+1,2)=node(sn)%value_c(problem%n+1,2)-node(sn_fe)%cvalue_c(k,1,1)*n_fn(k)
                          end if
                        end do
                        ! Solid phase
                        do k=1,problem%n
                          ! u_k^(2) = u_k^(FE)
                          if (node(sn_fe)%ctype(k,1).eq.1) then
                            node(sn)%value_c(k,2)=b_c(node(sn_fe)%col(k,1),1)
                          else
                            node(sn)%value_c(k,2)=node(sn_fe)%cvalue_c(k,1,1)
                          end if
                          ! t_k^(2)
                          node(sn)%value_c(problem%n+1+k,2)=b_c(node(sn)%col(problem%n+1+k,2),1)
                        end do

                    end select

                end select
              end if
            end do
          end do
        end do

        ! =====================
        ! COUPLED BE BODY LOADS
        ! =====================

        do kb=1,region(kr)%n_be_bodyloads
          sb=region(kr)%be_bodyload(kb)
          sp=be_bodyload(sb)%part
          select case (be_bodyload(sb)%coupling)

            ! ----------------------------------
            ! FE BEAM TIP - BE LINE/SURFACE LOAD
            ! ----------------------------------

            case (fbem_bl_coupling_beam_tip)
              stop 'not yet'

            ! -----------------------
            ! FE SHELL - BE EDGE LOAD
            ! -----------------------

            case (fbem_bl_coupling_shell_edge)
              stop 'not yet'

            ! -----------------------------------------------------
            ! FE BEAM - BE LINE LOAD AND FE SHELL - BE SURFACE LOAD
            ! -----------------------------------------------------

            case (fbem_bl_coupling_beam_line,fbem_bl_coupling_shell_surface)
              do ke=1,part(sp)%n_elements
                se=part(sp)%element(ke)
                do kn=1,element(se)%n_nodes
                  sn=element(se)%node(kn)
                  if (.not.node_used(sn)) then
                    node_used(sn)=.true.
                    sn_fe=element(se)%element_node(kn)
                    do k=1,problem%n
                      if (node(sn_fe)%ctype(k,1).eq.1) then
                        node(sn)%value_c(k,1)=b_c(node(sn_fe)%col(k,1),1)
                      else
                        node(sn)%value_c(k,1)=node(sn_fe)%cvalue_c(k,1,1)
                      end if
                      node(sn)%value_c(problem%n+k,1)=b_c(node(sn)%col(problem%n+k,1),1)
                    end do
                  end if
                end do
              end do

          end select

        end do

      ! ============================================================================================================================

      ! =========
      ! FE REGION
      ! =========

      case (fbem_fe)
        !
        ! FLEXIBLE REGIONS
        !
        if (region(kr)%type.ne.fbem_rigid) then
          do ks=1,region(kr)%n_fe_subregions
            ss=region(kr)%fe_subregion(ks)
            do ke=1,part(fe_subregion(ss)%part)%n_elements
              se=part(fe_subregion(ss)%part)%element(ke)
              do kn=1,element(se)%n_nodes
                sn=element(se)%node(kn)
                if (.not.node_used(sn)) then
                  node_used(sn)=.true.
                  do k=1,node(sn)%n_dof
                    if (node(sn)%ctype(k,1).eq.1) then
                      node(sn)%value_c(k,1)=b_c(node(sn)%col(k,1),1)
                    else
                      node(sn)%value_c(k,1)=node(sn)%cvalue_c(k,1,1)
                    end if
                  end do
                end if
              end do
            end do
          end do
        end if

    end select

  end do

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'END assigning the LSE solution to variables')

end subroutine assign_solution_mechanics_harmonic
