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

subroutine build_lse_mechanics_harmonic(kf)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_numerical
  use fbem_shape_functions
  use fbem_geometry
  use fbem_symmetry
  use fbem_quasisingular_integration
  use fbem_telles_transformation
  use fbem_harela_incident_field
  use fbem_bem_stapot2d
  use fbem_bem_harpot2d
  use fbem_bem_stapot3d
  use fbem_bem_staela2d
  use fbem_bem_harela2d
  use fbem_bem_staela3d
  use fbem_bem_harela3d
  use fbem_fem_beams
  use fbem_fem_shells
  use fbem_fem_solids

  ! Module of problem variables
  use problem_variables

  ! No implicit variables
  implicit none

  ! I/O variables
  integer                        :: kf
  ! Local variables
  real(kind=real64)              :: omega
  integer                        :: kr
  integer                        :: kn
  integer                        :: knm
  real(kind=real64), allocatable :: T(:,:)
  integer                        :: kc, kcj
  integer                        :: knc
  integer                        :: row, col
  ! Related to additional equations
  integer :: sb, sn
  integer :: kci
  complex(kind=real64) :: KT

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'START building the linear system of equations')

  ! Initialization
  A_c=0.
  b_c=0.
  omega=frequency(kf)

  ! ==========================================================================================
  ! BUILD AND ASSEMBLE BEM INFLUENCE MATRICES AND FEM STIFFNESS MATRICES AND DISTRIBUTED LOADS
  ! ==========================================================================================

  do kr=1,n_regions
    select case (region(kr)%class)
      case (fbem_be)
        select case (region(kr)%type)
          case (fbem_potential)
            call build_lse_mechanics_bem_harpot(kf,kr)
          case (fbem_viscoelastic)
            call build_lse_mechanics_bem_harela(kf,kr)
          case (fbem_poroelastic)
            call build_lse_mechanics_bem_harpor(kf,kr)
        end select
      case (fbem_fe)
        call build_lse_mechanics_fem_harela(kf,kr)
    end select
  end do

  ! ===================================
  ! ASSEMBLE FINITE ELEMENT NODAL LOADS
  ! ===================================
  !
  ! Nota: esto debe ser igual al caso estatico, excepto por cambiar _r por _c
  !
  do kn=1,n_nodes
    if (part(node(kn)%part(1))%type.eq.fbem_part_fe_subregion) then

      !
      ! Uncoupled nodal load
      !
      if (node(kn)%coupled_node.eq.0) then
        !
        ! Slave nodes
        !
        !
        ! ******
        ! chequear en la def de bc de nodos que solo cargas puntuales se pueden aplicar sobre nodos esclavos, no desplazamientos
        ! se hace aqui no obstante tambien
        ! ********
        !
        if (node(kn)%rigid_link.eq.2) then
          knm=node(kn)%master
          allocate (T(node(kn)%n_dof,node(knm)%n_dof))
          call fbem_rigid_solid_transformation_matrix(problem%n,node(kn)%x,node(knm)%x,node(kn)%n_dof,T)
          do kc=1,node(knm)%n_dof
            if (node(knm)%ctype(kc,1).eq.1) then
              row=node(knm)%row(kc,1)
              do knc=1,node(kn)%n_dof
                if (node(kn)%ctype(knc,1).eq.1) then
                  b_c(row,1)=b_c(row,1)+T(knc,kc)*node(kn)%cvalue_c(knc,1,1)
                else
                  call fbem_error_message(error_unit,0,'node',node(kn)%id,'slave nodes do not admit kinematic constraints')
                end if
              end do
            end if
          end do
          deallocate (T)
        !
        ! Normal and master nodes
        !
        else
          do kc=1,node(kn)%n_dof
            if (node(kn)%ctype(kc,1).eq.1) then
              row=node(kn)%row(kc,1)
              b_c(row,1)=b_c(row,1)+node(kn)%cvalue_c(kc,1,1)
            end if
          end do
        end if
      !
      ! Coupled nodal load
      !
      else
        stop 'coupled tip removed feature'
      end if
    end if
  end do

  ! ==========================
  ! BUILD ADDITIONAL EQUATIONS
  ! ==========================

  !
  ! Esto es para las condiciones de contacto deslizante y otras.... esto se deberia cambiar e introducir implicitamente
  ! haciendo un cambio de coordenadas nodales
  !

  ! Loop through NODES
  do kn=1,n_nodes
    sn=kn
    select case (part(node(sn)%part(1))%type)

      ! ============================================================================================================================
      ! BE BOUNDARY NODE
      !
      case (fbem_part_be_boundary)
        ! BE boundary of the node
        sb=part(node(sn)%part(1))%entity
        select case (boundary(sb)%class)

          ! =================
          ! ORDINARY BOUNDARY
          ! =================

          case (fbem_boundary_class_ordinary)

            select case (boundary(sb)%coupling)

              ! ----------- !
              ! BE BOUNDARY !
              ! ----------- !

              case (fbem_boundary_coupling_be)
                kr=boundary(sb)%region(1)
                select case (region(kr)%type)

                  ! -------------- !
                  ! INVISCID FLUID !
                  ! -------------- !

                  case (fbem_potential)

                  ! ------------------ !
                  ! VISCOELASTIC SOLID !
                  ! ------------------ !

                  case (fbem_viscoelastic)
                    ! Loop through each coordinate
                    do kc=1,problem%n
                      select case (node(sn)%ctype(kc,1))
                        ! Local axes: u·l = U
                        case (2)
                          row=node(sn)%row(kc,0)
                          select case (kc)
                            ! n direction
                            case (1)
                              do kci=1,problem%n
                                col=node(sn)%col(kci,1)
                                A_c(row,col)=node(sn)%n_fn(kci)
                              end do
                            ! t1 direction
                            case (2)
                              do kci=1,problem%n
                                col=node(sn)%col(kci,1)
                                A_c(row,col)=node(sn)%t1_fn(kci)
                              end do
                            ! t2 direction
                            case (3)
                              do kci=1,problem%n
                                col=node(sn)%col(kci,1)
                                A_c(row,col)=node(sn)%t2_fn(kci)
                              end do
                          end select
                          b_c(row,1)=node(sn)%cvalue_c(kc,1,1)
                        ! Local axes: t·l = T
                        case (3)
                          row=node(sn)%row(kc,0)
                          select case (kc)
                            ! n direction
                            case (1)
                              do kci=1,problem%n
                                col=node(sn)%col(problem%n+kci,1)
                                A_c(row,col)=node(sn)%n_fn(kci)
                              end do
                            ! t1 direction
                            case (2)
                              do kci=1,problem%n
                                col=node(sn)%col(problem%n+kci,1)
                                A_c(row,col)=node(sn)%t1_fn(kci)
                              end do
                            ! t2 direction
                            case (3)
                              do kci=1,problem%n
                                col=node(sn)%col(problem%n+kci,1)
                                A_c(row,col)=node(sn)%t2_fn(kci)
                              end do
                          end select
                          b_c(row,1)=node(sn)%cvalue_c(kc,1,1)
                      end select
                    end do

                  ! ----------------- !
                  ! POROELASTIC MEDIUM !
                  ! ----------------- !

                  case (fbem_poroelastic)
                    ! Loop through each coordinate k of the solid skeleton
                    do kc=1,problem%n
                      select case (node(sn)%ctype(kc,1))
                        ! Local axes: u·l = U
                        case (2,6)
                          row=node(sn)%row(kc,0)
                          select case (kc)
                            ! n direction
                            case (1)
                              do kci=1,problem%n
                                col=node(sn)%col(kci,1)
                                A_c(row,col)=node(sn)%n_fn(kci)
                              end do
                            ! t1 direction
                            case (2)
                              do kci=1,problem%n
                                col=node(sn)%col(kci,1)
                                A_c(row,col)=node(sn)%t1_fn(kci)
                              end do
                            ! t2 direction
                            case (3)
                              do kci=1,problem%n
                                col=node(sn)%col(kci,1)
                                A_c(row,col)=node(sn)%t2_fn(kci)
                              end do
                          end select
                          b_c(row,1)=node(sn)%cvalue_c(kc,1,1)
                        ! Local axes: t·l = T
                        case (3)
                          ! Row of the equation
                          row=node(sn)%row(kc,0)
                          select case (kc)
                            ! n direction
                            case (1)
                              do kci=1,problem%n
                                ! Assemble t_kci
                                col=node(sn)%col(problem%n+1+kci,1)
                                A_c(row,col)=node(sn)%n_fn(kci)
                              end do
                              ! Assemble tau depending on the boundary condition
                              select case (node(sn)%ctype(0,1))
                                ! Known tau
                                case (0)
                                  b_c(row,1)=b_c(row,1)+node(sn)%cvalue_c(0,1,1)
                                ! Unknown tau
                                case (1,2)
                                  col=node(sn)%col(0,1)
                                  A_c(row,col)=0
                              end select
                            ! t1 direction
                            case (2)
                              do kci=1,problem%n
                                col=node(sn)%col(problem%n+1+kci,1)
                                A_c(row,col)=node(sn)%t1_fn(kci)
                              end do
                            ! t2 direction
                            case (3)
                              do kci=1,problem%n
                                col=node(sn)%col(problem%n+1+kci,1)
                                A_c(row,col)=node(sn)%t2_fn(kci)
                              end do
                          end select
                          b_c(row,1)=b_c(row,1)+node(sn)%cvalue_c(kc,1,1)
                        ! Local axes: (t+tau*n)·l = T
                        case (7)
                          ! Row of the equation
                          row=node(sn)%row(kc,0)
                          select case (kc)
                            ! n direction
                            case (1)
                              do kci=1,problem%n
                                col=node(sn)%col(problem%n+1+kci,1)
                                A_c(row,col)=node(sn)%n_fn(kci)
                              end do
                              ! Assemble tau
                              col=node(sn)%col(0,1)
                              A_c(row,col)=1.d0
                            ! t1 direction
                            case (2)
                              do kci=1,problem%n
                                col=node(sn)%col(problem%n+1+kci,1)
                                A_c(row,col)=node(sn)%t1_fn(kci)
                              end do
                            ! t2 direction
                            case (3)
                              do kci=1,problem%n
                                col=node(sn)%col(problem%n+1+kci,1)
                                A_c(row,col)=node(sn)%t2_fn(kci)
                              end do
                          end select
                          b_c(row,1)=node(sn)%cvalue_c(kc,1,1)
                      end select
                    end do

                end select

              ! ================================================================================================================== !
              ! BE-BE BOUNDARY                                                                                                     !
              ! ================================================================================================================== !

              case (fbem_boundary_coupling_be_be)
                  select case (region(boundary(sb)%region(1))%type)
                    case (fbem_potential)
                      select case (region(boundary(sb)%region(2))%type)

                        ! ------------------------------------------------- !
                        ! BE (inviscid fluid) - BE (inviscid fluid)         !
                        ! ------------------------------------------------- !

                        case (fbem_potential)


                        ! ------------------------------------------------------- !
                        ! BE (1: inviscid fluid) - BE (2: viscoelastic solid)     !
                        ! ------------------------------------------------------- !

                        case (fbem_viscoelastic)

                        ! -------------------------------------------------- !
                        ! BE (1: inviscid fluid) - BE (2: poroelastic medium) !
                        ! -------------------------------------------------- !

                        case (fbem_poroelastic)


                      end select

                    case (fbem_viscoelastic)
                      select case (region(boundary(sb)%region(2))%type)

                        ! ---------------------------------------------------- !
                        ! BE (1: viscoelastic solid) - BE (2: inviscid fluid)  !
                        ! ---------------------------------------------------- !

                        case (fbem_potential)

                        ! ------------------------------------------------- !
                        ! BE (viscoelastic solid) - BE (viscoelastic solid) !
                        ! ------------------------------------------------- !

                        case (fbem_viscoelastic)
                          ! Loop through each coordinate
                          do kc=1,problem%n
                            select case (node(sn)%ctype(kc,1))
                              ! Local axes
                              ! 3 - perfect bonding - u_k^{1}, u_k^{2} and t_k^{1} are active since:
                              ! t_k^{(2)}=-t_k^{(1)}
                              ! (u^{2}-u^{1})·l_k^(1)=0
                              case (3)
                                row=node(sn)%row(kc,0)
                                select case (kc)
                                  ! n direction
                                  case (1)
                                    ! u^{1}
                                    do kci=1,problem%n
                                      col=node(sn)%col(kci,1)
                                      A_c(row,col)=-node(sn)%n_fn(kci)
                                    end do
                                    ! u^{2}
                                    do kci=1,problem%n
                                      col=node(sn)%col(kci,2)
                                      A_c(row,col)= node(sn)%n_fn(kci)
                                    end do
                                  ! t1 direction
                                  case (2)
                                    ! u^{1}
                                    do kci=1,problem%n
                                      col=node(sn)%col(kci,1)
                                      A_c(row,col)=-node(sn)%t1_fn(kci)
                                    end do
                                    ! u^{2}
                                    do kci=1,problem%n
                                      col=node(sn)%col(kci,2)
                                      A_c(row,col)= node(sn)%t1_fn(kci)
                                    end do
                                  ! t2 direction
                                  case (3)
                                    ! u^{1}
                                    do kci=1,problem%n
                                      col=node(sn)%col(kci,1)
                                      A_c(row,col)=-node(sn)%t2_fn(kci)
                                    end do
                                    ! u^{2}
                                    do kci=1,problem%n
                                      col=node(sn)%col(kci,2)
                                      A_c(row,col)= node(sn)%t2_fn(kci)
                                    end do
                                end select
                              ! Local axes
                              ! 4 - perfect debonding - u_k^{1}, u_k^{2} and t_k^{1} are active since:
                              ! t_k^{(2)}=-t_k^{(1)}
                              ! t^(1)·l_k^(1)=0
                              case (4)
                                row=node(sn)%row(kc,0)
                                select case (kc)
                                  ! n direction
                                  case (1)
                                    do kci=1,problem%n
                                      col=node(sn)%col(problem%n+kci,1)
                                      A_c(row,col)=node(sn)%n_fn(kci)
                                    end do
                                  ! t1 direction
                                  case (2)
                                    do kci=1,problem%n
                                      col=node(sn)%col(problem%n+kci,1)
                                      A_c(row,col)=node(sn)%t1_fn(kci)
                                    end do
                                  ! t2 direction
                                  case (3)
                                    do kci=1,problem%n
                                      col=node(sn)%col(problem%n+kci,1)
                                      A_c(row,col)=node(sn)%t2_fn(kci)
                                    end do
                                end select
                              ! Local axes
                              ! 5 - partial bonding - u_k^{1}, u_k^{2} and t_k^{1} are active since:
                              ! t_k^{(2)}=-t_k^{(1)}
                              ! (u^{2}-u^{1})·l_k^(1)*KT-t^(1)·l_k^(1)=0
                              case (5)
                                ! KT = K + i*omega*C - omega**2*M
                                KT=node(sn)%cvalue_c(kc,1,1)+c_im*omega*node(sn)%cvalue_c(kc,2,1)&
                                  -omega**2*node(sn)%cvalue_c(kc,3,1)
                                row=node(sn)%row(kc,0)
                                select case (kc)
                                  ! n direction
                                  case (1)
                                    ! u^{1}
                                    do kci=1,problem%n
                                      col=node(sn)%col(         kci,1)
                                      A_c(row,col)=-node(sn)%n_fn(kci)*KT
                                    end do
                                    ! u^{2}
                                    do kci=1,problem%n
                                      col=node(sn)%col(         kci,2)
                                      A_c(row,col)= node(sn)%n_fn(kci)*KT
                                    end do
                                    ! t^(1)
                                    do kci=1,problem%n
                                      col=node(sn)%col(problem%n+kci,1)
                                      A_c(row,col)=-node(sn)%n_fn(kci)
                                    end do
                                  ! t1 direction
                                  case (2)
                                    ! u^{1}
                                    do kci=1,problem%n
                                      col=node(sn)%col(          kci,1)
                                      A_c(row,col)=-node(sn)%t1_fn(kci)*KT
                                    end do
                                    ! u^{2}
                                    do kci=1,problem%n
                                      col=node(sn)%col(          kci,2)
                                      A_c(row,col)= node(sn)%t1_fn(kci)*KT
                                    end do
                                    ! t^(1)
                                    do kci=1,problem%n
                                      col=node(sn)%col(problem%n+kci,1)
                                      A_c(row,col)=-node(sn)%t1_fn(kci)
                                    end do
                                  ! t2 direction
                                  case (3)
                                    ! u^{1}
                                    do kci=1,problem%n
                                      col=node(sn)%col(          kci,1)
                                      A_c(row,col)=-node(sn)%t2_fn(kci)*KT
                                    end do
                                    ! u^{2}
                                    do kci=1,problem%n
                                      col=node(sn)%col(          kci,2)
                                      A_c(row,col)= node(sn)%t2_fn(kci)*KT
                                    end do
                                    ! t^(1)
                                    do kci=1,problem%n
                                      col=node(sn)%col(problem%n+kci,1)
                                      A_c(row,col)=-node(sn)%t2_fn(kci)
                                    end do
                                end select
                            end select
                          end do

                        ! ------------------------------------------------------ !
                        ! BE (1: viscoelastic solid) - BE (2: poroelastic medium) !
                        ! ------------------------------------------------------ !

                        case (fbem_poroelastic)
                          ! Loop through each coordinate
                          do kc=1,problem%n
                            select case (node(sn)%ctype(kc,1))
                              ! Local axes
                              ! In all cases, u_k^{1}, u_k^{2} and t_k^{1} are active, and t_k^{(2)}=-t_k^{(1)}+tau^(2)n_k^(1).
                              ! 3 - perfect bonding: (u^{2}-u^{1})·l_k^(1)=0
                              case (3)
                                row=node(sn)%row(kc,0)
                                select case (kc)
                                  ! n direction
                                  case (1)
                                    ! u^{1}
                                    do kci=1,problem%n
                                      col=node(sn)%col(kci,1)
                                      A_c(row,col)=-node(sn)%n_fn(kci)
                                    end do
                                    ! u^{2}
                                    do kci=1,problem%n
                                      col=node(sn)%col(kci,2)
                                      A_c(row,col)= node(sn)%n_fn(kci)
                                    end do
                                  ! t1 direction
                                  case (2)
                                    ! u^{1}
                                    do kci=1,problem%n
                                      col=node(sn)%col(kci,1)
                                      A_c(row,col)=-node(sn)%t1_fn(kci)
                                    end do
                                    ! u^{2}
                                    do kci=1,problem%n
                                      col=node(sn)%col(kci,2)
                                      A_c(row,col)= node(sn)%t1_fn(kci)
                                    end do
                                  ! t2 direction
                                  case (3)
                                    ! u^{1}
                                    do kci=1,problem%n
                                      col=node(sn)%col(kci,1)
                                      A_c(row,col)=-node(sn)%t2_fn(kci)
                                    end do
                                    ! u^{2}
                                    do kci=1,problem%n
                                      col=node(sn)%col(kci,2)
                                      A_c(row,col)= node(sn)%t2_fn(kci)
                                    end do
                                end select
                              ! Local axes
                              ! In all cases, u_k^{1}, u_k^{2} and t_k^{1} are active, and t_k^{(2)}=-t_k^{(1)}+tau^(2)n_k^(1).
                              ! 4 - perfect debonding: t^(1)·l_k^(1)=0
                              case (4)
                                row=node(sn)%row(kc,0)
                                select case (kc)
                                  ! n direction
                                  case (1)
                                    do kci=1,problem%n
                                      col=node(sn)%col(problem%n+kci,1)
                                      A_c(row,col)=node(sn)%n_fn(kci)
                                    end do
                                  ! t1 direction
                                  case (2)
                                    do kci=1,problem%n
                                      col=node(sn)%col(problem%n+kci,1)
                                      A_c(row,col)=node(sn)%t1_fn(kci)
                                    end do
                                  ! t2 direction
                                  case (3)
                                    do kci=1,problem%n
                                      col=node(sn)%col(problem%n+kci,1)
                                      A_c(row,col)=node(sn)%t2_fn(kci)
                                    end do
                                end select
                              ! Local axes
                              ! In all cases, u_k^{1}, u_k^{2} and t_k^{1} are active, and t_k^{(2)}=-t_k^{(1)}+tau^(2)n_k^(1).
                              ! 5 - partial bonding: (u^{2}-u^{1})·l_k^(1)*KT=t^(1)·l_k^(1)
                              case (5)
                                ! KT = K + i*omega*C - omega**2*M
                                KT=node(sn)%cvalue_c(kc,1,1)+c_im*omega*node(sn)%cvalue_c(kc,2,1)&
                                  -omega**2*node(sn)%cvalue_c(kc,3,1)
                                row=node(sn)%row(kc,0)
                                select case (kc)
                                  ! n direction
                                  case (1)
                                    ! u^{1}
                                    do kci=1,problem%n
                                      col=node(sn)%col(         kci,1)
                                      A_c(row,col)=-node(sn)%n_fn(kci)*KT
                                    end do
                                    ! u^{2}
                                    do kci=1,problem%n
                                      col=node(sn)%col(         kci,2)
                                      A_c(row,col)= node(sn)%n_fn(kci)*KT
                                    end do
                                    ! t^(1)
                                    do kci=1,problem%n
                                      col=node(sn)%col(problem%n+kci,1)
                                      A_c(row,col)=-node(sn)%n_fn(kci)
                                    end do
                                  ! t1 direction
                                  case (2)
                                    ! u^{1}
                                    do kci=1,problem%n
                                      col=node(sn)%col(          kci,1)
                                      A_c(row,col)=-node(sn)%t1_fn(kci)*KT
                                    end do
                                    ! u^{2}
                                    do kci=1,problem%n
                                      col=node(sn)%col(          kci,2)
                                      A_c(row,col)= node(sn)%t1_fn(kci)*KT
                                    end do
                                    ! t^(1)
                                    do kci=1,problem%n
                                      col=node(sn)%col(problem%n+kci,1)
                                      A_c(row,col)=-node(sn)%t1_fn(kci)
                                    end do
                                  ! t2 direction
                                  case (3)
                                    ! u^{1}
                                    do kci=1,problem%n
                                      col=node(sn)%col(          kci,1)
                                      A_c(row,col)=-node(sn)%t2_fn(kci)*KT
                                    end do
                                    ! u^{2}
                                    do kci=1,problem%n
                                      col=node(sn)%col(          kci,2)
                                      A_c(row,col)= node(sn)%t2_fn(kci)*KT
                                    end do
                                    ! t^(1)
                                    do kci=1,problem%n
                                      col=node(sn)%col(problem%n+kci,1)
                                      A_c(row,col)=-node(sn)%t2_fn(kci)
                                    end do
                                end select
                            end select
                          end do

                      end select

                    case (fbem_poroelastic)
                      select case (region(boundary(sb)%region(2))%type)

                        ! --------------------------------------------------- !
                        ! BE (1: poroelastic medium) - BE (2: inviscid fluid)  !
                        ! --------------------------------------------------- !

                        case (fbem_potential)

                        ! ------------------------------------------------------ !
                        ! BE (1: poroelastic medium) - BE (2: viscoelastic solid) !
                        ! ------------------------------------------------------ !

                        case (fbem_viscoelastic)
                          ! Loop through each coordinate
                          do kc=1,problem%n
                            select case (node(sn)%ctype(kc,1))
                              ! Local axes
                              ! In all cases, u_k^{1}, u_k^{2} and t_k^{2} are active, and t_k^{(1)}=-t_k^{(2)}-tau^(1)n_k^(1).
                              ! 3 - perfect bonding: (u^{1}-u^{2})·l_k^(1)=0
                              case (3)
                                row=node(sn)%row(kc,0)
                                select case (kc)
                                  ! n direction
                                  case (1)
                                    ! u^{2}
                                    do kci=1,problem%n
                                      col=node(sn)%col(kci,2)
                                      A_c(row,col)=-node(sn)%n_fn(kci)
                                    end do
                                    ! u^{1}
                                    do kci=1,problem%n
                                      col=node(sn)%col(kci,1)
                                      A_c(row,col)= node(sn)%n_fn(kci)
                                    end do
                                  ! t1 direction
                                  case (2)
                                    ! u^{2}
                                    do kci=1,problem%n
                                      col=node(sn)%col(kci,2)
                                      A_c(row,col)=-node(sn)%t1_fn(kci)
                                    end do
                                    ! u^{1}
                                    do kci=1,problem%n
                                      col=node(sn)%col(kci,1)
                                      A_c(row,col)= node(sn)%t1_fn(kci)
                                    end do
                                  ! t2 direction
                                  case (3)
                                    ! u^{2}
                                    do kci=1,problem%n
                                      col=node(sn)%col(kci,2)
                                      A_c(row,col)=-node(sn)%t2_fn(kci)
                                    end do
                                    ! u^{1}
                                    do kci=1,problem%n
                                      col=node(sn)%col(kci,1)
                                      A_c(row,col)= node(sn)%t2_fn(kci)
                                    end do
                                end select
                              ! Local axes
                              ! In all cases, u_k^{1}, u_k^{2} and t_k^{2} are active, and t_k^{(1)}=-t_k^{(2)}-tau^(1)n_k^(1).
                              ! 4 - perfect debonding: t^(2)·l_k^(1)=0
                              case (4)
                                row=node(sn)%row(kc,0)
                                select case (kc)
                                  ! n direction
                                  case (1)
                                    do kci=1,problem%n
                                      col=node(sn)%col(problem%n+kci,2)
                                      A_c(row,col)=node(sn)%n_fn(kci)
                                    end do
                                  ! t1 direction
                                  case (2)
                                    do kci=1,problem%n
                                      col=node(sn)%col(problem%n+kci,2)
                                      A_c(row,col)=node(sn)%t1_fn(kci)
                                    end do
                                  ! t2 direction
                                  case (3)
                                    do kci=1,problem%n
                                      col=node(sn)%col(problem%n+kci,2)
                                      A_c(row,col)=node(sn)%t2_fn(kci)
                                    end do
                                end select
                              ! Local axes
                              ! In all cases, u_k^{1}, u_k^{2} and t_k^{2} are active, and t_k^{(1)}=-t_k^{(2)}-tau^(1)n_k^(1).
                              ! 5 - partial bonding  : (u^{1}-u^{2})·l_k^(1)*KT=t^(2)·l_k^(1)
                              case (5)
                                ! KT = K + i*omega*C - omega**2*M
                                KT=node(sn)%cvalue_c(kc,1,1)+c_im*omega*node(sn)%cvalue_c(kc,2,1)&
                                  -omega**2*node(sn)%cvalue_c(kc,3,1)
                                row=node(sn)%row(kc,0)
                                select case (kc)
                                  ! n direction
                                  case (1)
                                    ! u^{2}
                                    do kci=1,problem%n
                                      col=node(sn)%col(         kci,2)
                                      A_c(row,col)=-node(sn)%n_fn(kci)*KT
                                    end do
                                    ! u^{1}
                                    do kci=1,problem%n
                                      col=node(sn)%col(         kci,1)
                                      A_c(row,col)= node(sn)%n_fn(kci)*KT
                                    end do
                                    ! t^(2)
                                    do kci=1,problem%n
                                      col=node(sn)%col(problem%n+kci,2)
                                      A_c(row,col)=-node(sn)%n_fn(kci)
                                    end do
                                  ! t1 direction
                                  case (2)
                                    ! u^{2}
                                    do kci=1,problem%n
                                      col=node(sn)%col(          kci,2)
                                      A_c(row,col)=-node(sn)%t1_fn(kci)*KT
                                    end do
                                    ! u^{1}
                                    do kci=1,problem%n
                                      col=node(sn)%col(          kci,1)
                                      A_c(row,col)= node(sn)%t1_fn(kci)*KT
                                    end do
                                    ! t^(2)
                                    do kci=1,problem%n
                                      col=node(sn)%col(problem%n+kci,2)
                                      A_c(row,col)=-node(sn)%t1_fn(kci)
                                    end do
                                  ! t2 direction
                                  case (3)
                                    ! u^{2}
                                    do kci=1,problem%n
                                      col=node(sn)%col(          kci,2)
                                      A_c(row,col)=-node(sn)%t2_fn(kci)*KT
                                    end do
                                    ! u^{1}
                                    do kci=1,problem%n
                                      col=node(sn)%col(          kci,1)
                                      A_c(row,col)= node(sn)%t2_fn(kci)*KT
                                    end do
                                    ! t^(2)
                                    do kci=1,problem%n
                                      col=node(sn)%col(problem%n+kci,2)
                                      A_c(row,col)=-node(sn)%t2_fn(kci)
                                    end do
                                end select
                            end select
                          end do

                        ! ----------------------------------------------- !
                        ! BE (poroelastic medium) - BE (poroelastic medium) !
                        ! ----------------------------------------------- !

                        case (fbem_poroelastic)

                      end select

                  end select

              case (fbem_boundary_coupling_be_fe)

              case (fbem_boundary_coupling_be_fe_be)

            end select

          ! =================== !
          ! CRACK-LIKE BOUNDARY !
          ! =================== !

          case (fbem_boundary_class_cracklike)

            select case (boundary(sb)%coupling)

              case (fbem_boundary_coupling_be)
                select case (region(boundary(sb)%region(1))%type)

                  ! -------------- !
                  ! INVISCID FLUID !
                  ! -------------- !

                  case (fbem_potential)

                  ! ------------------ !
                  ! VISCOELASTIC SOLID !
                  ! ------------------ !

                  case (fbem_viscoelastic)

                  ! ----------------- !
                  ! POROELASTIC MEDIUM !
                  ! ----------------- !

                  case (fbem_poroelastic)
                    ! Loop through each coordinate
                    do kc=1,problem%n
                      select case (node(sn)%ctype(kc,1))
                        ! Fluid(inviscid/incompressible)-filled impermeable crack
                        ! Local axes
                        ! kc==1: (u^--u^+)·n^+=0
                        ! kc==2: t^+·t_1^+=0
                        ! kc==3: t^+·t_2^+=0
                        case (5)
                          row=node(sn)%row(kc,0)
                          select case (kc)
                            ! n direction: (u^--u^+)·n^+=0
                            case (1)
                              ! u^+
                              do kci=1,problem%n
                                col=node(sn)%col(kci,1)
                                A_c(row,col)=-node(sn)%n_fn(kci)
                              end do
                              ! u^-
                              do kci=1,problem%n
                                col=node(sn)%col(kci,2)
                                A_c(row,col)= node(sn)%n_fn(kci)
                              end do
                            ! t1 direction: t^+·t_1^+=0
                            case (2)
                              do kci=1,problem%n
                                col=node(sn)%col(problem%n+1+kci,1)
                                A_c(row,col)=node(sn)%t1_fn(kci)
                              end do
                            ! t2 direction: t^+·t_2^+=0
                            case (3)
                              do kci=1,problem%n
                                col=node(sn)%col(problem%n+1+kci,1)
                                A_c(row,col)=node(sn)%t2_fn(kci)
                              end do
                          end select
                      end select
                    end do

                end select

              case (fbem_boundary_coupling_be_fe)

            end select

        end select
      ! ============================================================================================================================

      ! ============================================================================================================================
      ! FE SUBREGION NODE
      !
      case (fbem_part_fe_subregion)
      ! ============================================================================================================================
    end select

  end do ! Loop through NODES

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'END building the linear system of equations')

end subroutine build_lse_mechanics_harmonic
