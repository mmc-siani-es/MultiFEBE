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
subroutine assign_solution_mechanics_static_sa

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_data_structures
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_geometry
  use fbem_symmetry

  ! Module of problem variables
  use problem_variables

  ! No implicit variables are allowed
  implicit none

  ! Local variables
  integer                    :: kr, kb, ke, kn, ks, ss
  integer                    :: sb, se, sn
  logical                    :: sb_reversion
  integer                    :: k
  logical, allocatable       :: node_used(:)
  integer                    :: sn_fe
  real(kind=real64)          :: x_fn(problem%n), n_fn(problem%n)

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'START assigning the LSE solution to variables (sensitivity)')

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
                        ! Index of variables for each coordinate k:
                        ! node(sn)%col(          k,1): u_k
                        ! node(sn)%col(problem%n+k,1): t_k
                        do k=1,problem%n
                          select case (node(sn)%ctype(k,1))
                            case (0)
                              node(sn)%dvda_r(          k,1,:)=0
                              node(sn)%dvda_r(problem%n+k,1,:)=bsa_r(node(sn)%col(problem%n+k,1),:)
                            case (1,10)
                              node(sn)%dvda_r(          k,1,:)=bsa_r(node(sn)%col(k,1),:)
                              node(sn)%dvda_r(problem%n+k,1,:)=0
                            case (2,3)
                              node(sn)%dvda_r(          k,1,:)=bsa_r(node(sn)%col(k,1),:)
                              node(sn)%dvda_r(problem%n+k,1,:)=bsa_r(node(sn)%col(problem%n+k,1),:)
                          end select
                        end do

                      ! ===================
                      ! CRACK-LIKE BOUNDARY
                      ! ===================

                      case (fbem_boundary_class_cracklike)
                        ! Index of variables for each coordinate k:
                        ! node(sn)%col(          k,1): u_k for face +
                        ! node(sn)%col(problem%n+k,1): t_k for face +
                        ! node(sn)%col(          k,2): u_k for face -
                        ! node(sn)%col(problem%n+k,2): t_k for face -
                        ! Face +
                        do k=1,problem%n
                          select case (node(sn)%ctype(k,1))
                            case (0)
                              node(sn)%dvda_r(          k,1,:)=0
                              node(sn)%dvda_r(problem%n+k,1,:)=bsa_r(node(sn)%col(problem%n+k,1),:)
                            case (1)
                              node(sn)%dvda_r(          k,1,:)=bsa_r(node(sn)%col(k,1),:)
                              node(sn)%dvda_r(problem%n+k,1,:)=0
                          end select
                        end do
                        ! Face -
                        do k=1,problem%n
                          select case (node(sn)%ctype(k,2))
                            case (0)
                              node(sn)%dvda_r(          k,2,:)=0
                              node(sn)%dvda_r(problem%n+k,2,:)=bsa_r(node(sn)%col(problem%n+k,2),:)
                            case (1)
                              node(sn)%dvda_r(          k,2,:)=bsa_r(node(sn)%col(k,2),:)
                              node(sn)%dvda_r(problem%n+k,2,:)=0
                          end select
                        end do

                    end select

                  ! ==============
                  ! BE-BE BOUNDARY
                  ! ==============
                  ! The region 1 is the region where the boundary has N+
                  ! The region 2 is the region where the boundary has N-

                  case (fbem_boundary_coupling_be_be)
                    ! Index of variables for each coordinate k:
                    ! node(sn)%col(          k,1): u_k for region 1
                    ! node(sn)%col(problem%n+k,1): t_k for region 1
                    ! node(sn)%col(          k,2): u_k for region 2
                    ! node(sn)%col(problem%n+k,2): t_k for region 2
                    ! perfect bonding
                    do k=1,problem%n
                      node(sn)%dvda_r(          k,1,:)=bsa_r(node(sn)%col(k,1),:)
                      node(sn)%dvda_r(problem%n+k,1,:)=bsa_r(node(sn)%col(problem%n+k,1),:)
                      node(sn)%dvda_r(          k,2,:)= node(sn)%dvda_r(          k,1,:)
                      node(sn)%dvda_r(problem%n+k,2,:)=-node(sn)%dvda_r(problem%n+k,1,:)
                    end do

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
                        ! Index of variables for each coordinate k:
                        ! node(sn)%col(          k,1): u_k
                        ! node(sn)%col(problem%n+k,1): t_k
                        do k=1,problem%n
                          if (node(sn_fe)%ctype(k,1).eq.1) then
                            node(sn)%dvda_r(k,1,:)=bsa_r(node(sn_fe)%col(k,1),:)
                          else
                            node(sn)%dvda_r(k,1,:)=0
                          end if
                          node(sn)%dvda_r(problem%n+k,1,:)=bsa_r(node(sn)%col(problem%n+k,1),:)
                        end do

                      ! ===================
                      ! CRACK-LIKE BOUNDARY
                      ! ===================

                      case (fbem_boundary_class_cracklike)
                        ! Index of variables for each coordinate k:
                        ! node(sn)%col(          k,1): u_k for face +
                        ! node(sn)%col(problem%n+k,1): t_k for face +
                        ! node(sn)%col(          k,2): u_k for face - (not used since u_l_{+} = u_k_{-})
                        ! node(sn)%col(problem%n+k,2): t_k for face -
                        do k=1,problem%n
                          if (node(sn_fe)%ctype(k,1).eq.1) then
                            node(sn)%dvda_r(k,1,:)=bsa_r(node(sn_fe)%col(k,1),:)
                          else
                            node(sn)%dvda_r(k,1,:)=0
                          end if
                          node(sn)%dvda_r(problem%n+k,1,:)=bsa_r(node(sn)%col(problem%n+k,1),:)
                          node(sn)%dvda_r(          k,2,:)=node(sn)%dvda_r(k,1,:)
                          node(sn)%dvda_r(problem%n+k,2,:)=bsa_r(node(sn)%col(problem%n+k,2),:)
                        end do

                    end select

                  ! =================
                  ! BE-FE-BE BOUNDARY
                  ! =================

                  case (fbem_boundary_coupling_be_fe_be)
                    stop 'not yet 38'

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
                    node(sn)%dvda_r(k,1,:)=bsa_r(node(sn)%col(k,1),:)
                  else
                    node(sn)%dvda_r(k,1,:)=0
                  end if
                end do
              end if
            end do
          end do
        end do

    end select

  end do

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'START assigning the LSE solution to variables (sensitivity)')

end subroutine assign_solution_mechanics_static_sa
