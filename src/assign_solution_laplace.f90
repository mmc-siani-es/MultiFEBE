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

subroutine assign_solution_laplace

  !
  ! !!!!!!!??¿?¿?¿?¿?¿
  ! para elementos discontinuos hay que transformar los valores en los nodos funcionales por los valores en los nodos geometricos
  ! hacerlo al final, despues de calcularlos para los nodos funcionales
  !
  ! poner también como variables de cada nodo, el campo incidente
  !

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

  ! Local variables
  integer              :: kr, kb, ke, kn
  integer              :: sb, se, sn
  logical              :: sb_reversion
  integer              :: k
  logical, allocatable :: node_used(:)
  real(kind=real64)    :: x_fn(problem%n), n_fn(problem%n)

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'START assigning the LSE solution to variables')

  allocate (node_used(n_nodes))
  node_used=.false.

  do kr=1,n_regions
    !
    ! Switch between region classes
    !
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

                      ! ================= !
                      ! ORDINARY BOUNDARY !
                      ! ================= !

                      case (fbem_boundary_class_ordinary)
                        ! Index of variables:
                        ! node(sn)%col(1,1): p
                        ! node(sn)%col(2,1): j
                        select case (node(sn)%ctype(1,1))
                          case (0)
                            node(sn)%value_r(1,1)=node(sn)%cvalue_r(1,1,1)
                            node(sn)%value_r(2,1)=b_r(node(sn)%col(2,1),1)
                          case (1)
                            node(sn)%value_r(1,1)=b_r(node(sn)%col(1,1),1)
                            node(sn)%value_r(2,1)=node(sn)%cvalue_r(1,1,1)
                        end select


                      ! =================== !
                      ! CRACK-LIKE BOUNDARY !
                      ! =================== !

                      case (fbem_boundary_class_cracklike)
                        ! Index of variables:
                        ! node(sn)%col(1,1): p for face +
                        ! node(sn)%col(2,1): j for face +
                        ! node(sn)%col(1,2): p for face -
                        ! node(sn)%col(2,2): j for face -
                        ! Face +
                        select case (node(sn)%ctype(1,1))
                          case (0)
                            node(sn)%value_r(1,1)=node(sn)%cvalue_r(1,1,1)
                            node(sn)%value_r(2,1)=b_r(node(sn)%col(2,1),1)
                          case (1)
                            node(sn)%value_r(1,1)=b_r(node(sn)%col(1,1),1)
                            node(sn)%value_r(2,1)=node(sn)%cvalue_r(1,1,1)
                        end select
                        ! Face -
                        select case (node(sn)%ctype(1,2))
                          case (0)
                            node(sn)%value_r(1,2)=node(sn)%cvalue_r(1,1,2)
                            node(sn)%value_r(2,2)=b_r(node(sn)%col(2,2),1)
                          case (1)
                            node(sn)%value_r(1,2)=b_r(node(sn)%col(1,2),1)
                            node(sn)%value_r(2,2)=node(sn)%cvalue_r(1,1,2)
                        end select

                    end select

                  ! ==============
                  ! BE-BE BOUNDARY
                  ! ==============
                  ! The region 1 is the region where the boundary has N+
                  ! The region 2 is the region where the boundary has N-

                  case (fbem_boundary_coupling_be_be)
                    ! Index of variables:
                    ! node(sn)%col(1,1): p for region 1
                    ! node(sn)%col(2,1): j for region 1
                    ! node(sn)%col(1,2): p for region 2 (not used since p_{region 1}=p_{region 2})
                    ! node(sn)%col(2,2): j for region 2 (not used since j_{region 1}=-j_{region 2})
                    node(sn)%value_r(1,1)=b_r(node(sn)%col(1,1),1)
                    node(sn)%value_r(2,1)=b_r(node(sn)%col(2,1),1)
                    node(sn)%value_r(1,2)= node(sn)%value_r(1,1)
                    node(sn)%value_r(2,2)=-node(sn)%value_r(2,1)


                  ! ==============
                  ! BE-FE BOUNDARY
                  ! ==============

                  case (fbem_boundary_coupling_be_fe)
                    stop 'not implemented yet'

                  ! =================
                  ! BE-FE-BE BOUNDARY
                  ! =================

                  case (fbem_boundary_coupling_be_fe_be)
                    stop 'not implemented yet'
                end select
              end if
            end do
          end do
        end do


      ! =========
      ! FE REGION
      ! =========

      case (fbem_fe)
        stop 'not implemented yet'

    end select
  end do

  ! Ending message
  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'END assigning the LSE solution to variables')

end subroutine assign_solution_laplace
