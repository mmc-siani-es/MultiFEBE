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

subroutine build_auxiliary_variables_mechanics_harmonic

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
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
  integer                    :: kr, kp, kb, ke, kn, kn2
  integer                    :: sb, se, sn, ks, ss, sp, snm
  integer                    :: row, col
  integer                    :: k
  integer                    :: ndof_fe_node
  integer                    :: k_start, k_end, k_start2, k_end2, n_faces
  logical, allocatable       :: node_used(:)
  character(len=fbem_fmtstr) :: fmtstr
  integer(kind=int64)        :: memory

  if (verbose_level.ge.1)  call fbem_timestamp_w_message(output_unit,2,'START building auxiliary variables')

  ! Allocate auxiliary variable
  allocate (node_used(n_nodes))

  ! Initialize LSE row (equations) and columns (variables) counters
  row=1
  col=1
  ! Initialize the node usage indicator
  node_used=.false.

  ! ================================================================================================================================
  ! BE REGIONS
  ! ================================================================================================================================

  do kr=1,n_regions
    if (region(kr)%class.eq.fbem_be) then

      ! ====================
      ! BE BOUNDARY ELEMENTS
      ! ====================

      do kb=1,region(kr)%n_boundaries
        sb=region(kr)%boundary(kb)
        sp=boundary(sb)%part
        do ke=1,part(sp)%n_elements
          se=part(sp)%element(ke)
          do kn=1,element(se)%n_nodes
            sn=element(se)%node(kn)
            if (node_used(sn).eqv.(.false.)) then
              node_used(sn)=.true.
              select case (boundary(sb)%coupling)

                ! ================================================================================================================ !
                ! BE BOUNDARY                                                                                                      !
                ! ================================================================================================================ !

                case (fbem_boundary_coupling_be)
                  select case (boundary(sb)%class)

                    ! ================= !
                    ! ORDINARY BOUNDARY !
                    ! ================= !

                    case (fbem_boundary_class_ordinary)
                      select case (region(kr)%type)

                        ! -------------- !
                        ! INVISCID FLUID !
                        ! -------------- !

                        case (fbem_potential)

                          !
                          ! Index of equations:
                          ! node(sn)%row(1,1): SBIE, HBIE or SBIE+beta*HBIE
                          !
                          ! Index of variables:
                          ! node(sn)%col(1,1): p
                          ! node(sn)%col(2,1): Un
                          !
                          allocate (node(sn)%row(1,1))
                          allocate (node(sn)%col(2,1))
                          allocate (node(sn)%value_c(2,1))
                          allocate (node(sn)%dvda_c(2,1,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equation and column to the variables
                          ! Equation
                          node(sn)%row(1,1)=row
                          ! Increment counter
                          row=row+1
                          ! Assign values depending on the boundary condition
                          select case (node(sn)%ctype(1,1))
                            ! p known, Un unknown
                            case (0)
                              node(sn)%col(2,1)=col
                            ! Un known, p unknown
                            case (1)
                              node(sn)%col(1,1)=col
                            ! 2: p unknown, Un=i/(rho*c*omega)p
                            ! 3: p unknown, Un-(i/(rho*c*omega)-1/(2*R*rho*omega^2))p=0
                            case (2,3)
                              node(sn)%col(1,1)=col
                          end select
                          ! Increment counter
                          col=col+1

                        ! ------------------ !
                        ! VISCOELASTIC SOLID !
                        ! ------------------ !

                        case (fbem_viscoelastic)
                          !
                          ! Index of equations for each coordinate k:
                          ! node(sn)%row(k,1): SBIE, HBIE or SBIE+beta*HBIE
                          ! Additional equations for B.C. (only for local axis B.C.):
                          ! node(sn)%row(k,0): B.C.
                          !
                          ! Index of variables for each coordinate k:
                          ! node(sn)%col(          k,1): u_k
                          ! node(sn)%col(problem%n+k,1): t_k
                          !
                          allocate (node(sn)%row(problem%n,0:1))
                          allocate (node(sn)%col(2*problem%n,1))
                          allocate (node(sn)%value_c(2*problem%n,1))
                          allocate (node(sn)%dvda_c(2*problem%n,1,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equation and column to the variables
                          do k=1,problem%n
                            ! BIE equation
                            node(sn)%row(k,1)=row
                            ! Increment counter
                            row=row+1
                            ! Assign values depending on the boundary condition
                            select case (node(sn)%ctype(k,1))
                              ! u_k known, t_k unknown
                              case (0)
                                node(sn)%col(problem%n+k,1)=col
                                ! Increment counter
                                col=col+1
                              ! t_k known, u_k unknown
                              case (1,10)
                                node(sn)%col(          k,1)=col
                                ! Increment counter
                                col=col+1
                              ! u_k unknown, t_k unknown
                              case (2,3)
                                node(sn)%col(          k,1)=col
                                node(sn)%col(problem%n+k,1)=col+1
                                ! Increment counter
                                col=col+2
                                ! B.C. equation
                                node(sn)%row(k,0)=row
                                ! Increment counter
                                row=row+1
                            end select
                          end do

                        ! ----------------- !
                        ! POROELASTIC MEDIUM !
                        ! ----------------- !

                        case (fbem_poroelastic)
                          !
                          ! Index of equations:
                          ! node(sn)%row(0,1): SBIE, HBIE or SBIE+beta*HBIE (fluid)
                          ! node(sn)%row(k,1): SBIE, HBIE or SBIE+beta*HBIE (solid skeleton)
                          ! Additional equations for B.C. (only for local axis B.C.):
                          ! node(sn)%row(k,0): B.C.
                          !
                          ! Index of variables:
                          ! node(sn)%col(            0,1): tau
                          ! node(sn)%col(            k,1): u_k
                          ! node(sn)%col(problem%n+1  ,1): w
                          ! node(sn)%col(problem%n+1+k,1): t_k
                          !
                          allocate (node(sn)%row(0:problem%n,0:1))
                          allocate (node(sn)%col(0:(1+2*problem%n),1))
                          allocate (node(sn)%value_c(0:(1+2*problem%n),1))
                          allocate (node(sn)%dvda_c(0:(1+2*problem%n),1,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          !
                          ! Fluid phase
                          !
                          ! Assign row to the equation and column to the variables
                          ! BIE equation
                          node(sn)%row(0,1)=row
                          row=row+1
                          ! Unknowns depending on the boundary condition
                          select case (node(sn)%ctype(0,1))
                            !
                            ! open pore / permeable
                            !
                            ! tau known, Un unknown
                            case (0)
                              node(sn)%col(problem%n+1,1)=col
                            ! Un known, tau unknown
                            case (1)
                              node(sn)%col(          0,1)=col
                            !
                            ! close pore / impermeable
                            !
                            ! Un=u_k·n_k known, tau unknown
                            case (2)
                              node(sn)%col(          0,1)=col
                          end select
                          col=col+1
                          !
                          ! Solid skeleton
                          !
                          ! Assign row to the equation and column to the variables
                          do k=1,problem%n
                            ! BIE equation
                            node(sn)%row(k,1)=row
                            row=row+1
                            ! Assign values depending on the boundary condition
                            select case (node(sn)%ctype(k,1))
                              !
                              ! open pore / permeable
                              !
                              ! u_k known, t_k unknown
                              case (0)
                                node(sn)%col(problem%n+1+k,1)=col
                                ! Increment counter
                                col=col+1
                              ! t_k known, u_k unknown
                              case (1,10)
                                node(sn)%col(            k,1)=col
                                ! Increment counter
                                col=col+1
                              ! u_k unknown, t_k unknown
                              case (2,3)
                                node(sn)%col(            k,1)=col
                                node(sn)%col(problem%n+1+k,1)=col+1
                                ! Increment counter
                                col=col+2
                                ! B.C. equation
                                node(sn)%row(k,0)=row
                                ! Increment counter
                                row=row+1
                              !
                              ! close pore / impermeable
                              !
                              ! u_k known, t_k unknown
                              case (4)
                                node(sn)%col(problem%n+1+k,1)=col
                                ! Increment counter
                                col=col+1
                              ! t_k+tau·n_k known, u_k unknown
                              case (5,50)
                                node(sn)%col(            k,1)=col
                                ! Increment counter
                                col=col+1
                              ! u_k unknown, t_k unknown
                              case (6,7)
                                node(sn)%col(            k,1)=col
                                node(sn)%col(problem%n+1+k,1)=col+1
                                ! Increment counter
                                col=col+2
                                ! B.C. equation
                                node(sn)%row(k,0)=row
                                ! Increment counter
                                row=row+1
                            end select
                          end do

                      end select

                    ! =================== !
                    ! CRACK-LIKE BOUNDARY !
                    ! =================== !

                    case (fbem_boundary_class_cracklike)
                      select case (region(kr)%type)

                        ! -------------- !
                        ! INVISCID FLUID !
                        ! -------------- !

                        case (fbem_potential)
                          !
                          ! Index of equations (Dual BEM formulation):
                          ! node(sn)%row(1,1): SBIE
                          ! node(sn)%row(1,2): HBIE
                          !
                          ! Index of variables:
                          ! node(sn)%col(1,1): p for face +
                          ! node(sn)%col(2,1): Un for face +
                          ! node(sn)%col(1,2): p for face -
                          ! node(sn)%col(2,2): Un for face -
                          !
                          allocate (node(sn)%row(2,2))
                          allocate (node(sn)%col(2,2))
                          allocate (node(sn)%value_c(2,2))
                          allocate (node(sn)%dvda_c(2,2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equation and column to the variables
                          ! Equation for SBIE
                          node(sn)%row(1,1)=row
                          ! Equation for HBIE
                          node(sn)%row(1,2)=row+1
                          ! Increment counter
                          row=row+2
                          ! Face + / CC1
                          ! Assign values depending on the boundary condition
                          select case (node(sn)%ctype(1,1))
                            ! p known, Un unknown
                            case (0)
                              node(sn)%col(2,1)=col
                            ! Un known, p unknown
                            case (1)
                              node(sn)%col(1,1)=col
                          end select
                          ! Increment counter
                          col=col+1
                          ! Face - / CC2
                          ! Assign values depending on the boundary condition
                          select case (node(sn)%ctype(1,2))
                            ! p known, Un unknown
                            case (0)
                              node(sn)%col(2,2)=col
                            ! Un known, p unknown
                            case (1)
                              node(sn)%col(1,2)=col
                          end select
                          ! Increment counter
                          col=col+1

                        ! ------------------ !
                        ! VISCOELASTIC SOLID !
                        ! ------------------ !

                        case (fbem_viscoelastic)
                          !
                          ! Index of equations for each coordinate k (Dual BEM formulation):
                          ! node(sn)%row(k,1): SBIE
                          ! node(sn)%row(k,2): HBIE
                          !
                          ! Index of variables for each coordinate k:
                          ! node(sn)%col(          k,1): u_k for face +
                          ! node(sn)%col(problem%n+k,1): t_k for face +
                          ! node(sn)%col(          k,2): u_k for face -
                          ! node(sn)%col(problem%n+k,2): t_k for face -
                          !
                          allocate (node(sn)%row(problem%n,2))
                          allocate (node(sn)%col(2*problem%n,2))
                          allocate (node(sn)%value_c(2*problem%n,2))
                          allocate (node(sn)%dvda_c(2*problem%n,2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equation and column to the variables
                          do k=1,problem%n
                            ! Equation for SBIE
                            node(sn)%row(k,1)=row
                            ! Equation for HBIE
                            node(sn)%row(k,2)=row+1
                            ! Increment counter
                            row=row+2
                            ! Face +
                            ! Assign values depending on the boundary condition
                            select case (node(sn)%ctype(k,1))
                              ! u_k known, t_k unknown
                              case (0)
                                node(sn)%col(problem%n+k,1)=col
                              ! t_k known, u_k unknown
                              case (1)
                                node(sn)%col(          k,1)=col
                            end select
                            ! Increment counter
                            col=col+1
                            ! Face -
                            ! Assign values depending on the boundary condition
                            select case (node(sn)%ctype(k,2))
                              ! u_k known, t_k unknown
                              case (0)
                                node(sn)%col(problem%n+k,2)=col
                              ! t_k known, u_k unknown
                              case (1)
                                node(sn)%col(          k,2)=col
                            end select
                            ! Increment counter
                            col=col+1
                          end do

                        ! ------------------ !
                        ! POROELASTIC MEDIUM !
                        ! ------------------ !

                        case (fbem_poroelastic)
                          !
                          ! Index of equations (Dual BEM formulation):
                          ! node(sn)%row(0,1): SBIE (fluid)
                          ! node(sn)%row(k,1): SBIE (solid skeleton)
                          ! node(sn)%row(0,2): HBIE (fluid)
                          ! node(sn)%row(k,2): HBIE (solid skeleton)
                          ! Index of additional equations for each coordinate k:
                          ! node(sn)%row(k,0): compatibility or equilibrium equation
                          !
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
                          allocate (node(sn)%row(0:problem%n,0:2))
                          allocate (node(sn)%col(0:(1+2*problem%n),2))
                          allocate (node(sn)%value_c(0:(1+2*problem%n),2))
                          allocate (node(sn)%dvda_c(0:(1+2*problem%n),2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equation and column to the variables
                          do k=0,problem%n
                            ! Equation 1
                            node(sn)%row(k,1)=row
                            ! Equation 2
                            node(sn)%row(k,2)=row+1
                            ! Increment counter
                            row=row+2
                            ! Face +
                            ! Assign values depending on the boundary condition
                            select case (node(sn)%ctype(k,1))
                              ! open pore / permeable
                              ! If k==0: tau^+ known, w^+ unknown
                              ! If k>=1: u_k^+ known, t_k^+ unknown
                              case (0)
                                node(sn)%col(problem%n+1+k,1)=col
                                col=col+1
                              ! open pore / permeable
                              ! If k==0: w^+ known, tau^+ unknown
                              ! If k>=1: t_k^+ known, u_k^+ unknown
                              case (1)
                                node(sn)%col(            k,1)=col
                                col=col+1
                              ! close pore / impermeable
                              ! If k==0: Un^+=u_k^+·n_k^+ known, tau^+ unknown
                              ! If k>=1: u_k^+ known, t_k^+ unknown
                              case (2)
                                if (k.eq.0) then
                                  node(sn)%col(            k,1)=col
                                else
                                  node(sn)%col(problem%n+1+k,1)=col
                                end if
                                col=col+1
                              ! close pore / impermeable
                              ! Only k>=1: t_k^++tau^+·n_k^+ known, u_k^+ unknown
                              case (3)
                                node(sn)%col(            k,1)=col
                                col=col+1
                              ! Fluid(inviscid/incompressible)-filled permeable crack
                              ! If k==0: tau^+ unknown, Un^+ unknown
                              ! If k>=1: u^+ unknown, t_k^+=(1/phi-1)tau^+·n_k^+ known
                              case (4)
                                if (k.eq.0) then
                                  node(sn)%col(            k,1)=col
                                  node(sn)%col(problem%n+1+k,1)=col+1
                                  col=col+2
                                else
                                  node(sn)%col(            k,1)=col
                                  col=col+1
                                end if
                              ! Fluid(inviscid/incompressible)-filled impermeable crack
                              ! If k==0: tau^+ unknown, Un^+=u^+·n^+ known
                              ! If k>=1: u^+ unknown, t_k^+ unknown
                              case (5)
                                if (k.eq.0) then
                                  node(sn)%col(            k,1)=col
                                  col=col+1
                                else
                                  node(sn)%col(            k,1)=col
                                  node(sn)%col(problem%n+1+k,1)=col+1
                                  col=col+2
                                  ! Additional equations (k=1: compatibility in normal dir., k=2,3: tangential null stress)
                                  node(sn)%row(k,0)=row
                                  row=row+1
                                end if
                            end select
                            ! Face -
                            ! Assign values depending on the boundary condition
                            select case (node(sn)%ctype(k,2))
                              ! open pore / permeable
                              ! If k==0: tau^- known, w^- unknown
                              ! If k>=1: u_k^- known, t_k^- unknown
                              case (0)
                                node(sn)%col(problem%n+1+k,2)=col
                                col=col+1
                              ! open pore / permeable
                              ! If k==0: w^- known, tau^- unknown
                              ! If k>=1: t_k^- known, u_k^- unknown
                              case (1)
                                node(sn)%col(            k,2)=col
                                col=col+1
                              ! close pore / impermeable
                              ! If k==0: Un^-=u_k^-·n_k^- known, tau^- unknown
                              ! If k>=1: u_k^- known, t_k^- unknown
                              case (2)
                                if (k.eq.0) then
                                  node(sn)%col(            k,2)=col
                                else
                                  node(sn)%col(problem%n+1+k,2)=col
                                end if
                                col=col+1
                              ! close pore / impermeable
                              ! Only for k>=1: t_k^-+tau^-·n_k^- known, u_k^- unknown
                              case (3)
                                node(sn)%col(            k,2)=col
                                col=col+1
                              ! Fluid(inviscid/incompressible)-filled permeable crack
                              ! If k==0: tau^- and Un^- known in function of other variables
                              ! If k>=1: u^- unknown
                              case (4)
                                if (k.ge.1) then
                                  node(sn)%col(            k,2)=col
                                  col=col+1
                                end if
                              ! Fluid(inviscid/incompressible)-filled impermeable crack
                              ! If k==0: tau^- unknown, Un^-=u^-·n^- known
                              ! If k>=1: u^- unknown, t_k^-=(tau^--tau^+)n_k^+-t_k^+ known
                              case (5)
                                node(sn)%col(            k,2)=col
                                col=col+1
                            end select
                          end do

                      end select

                  end select

                ! ================================================================================================================ !
                ! BE-BE BOUNDARY                                                                                                   !
                ! ================================================================================================================ !
                !
                ! The region 1 is the region where the boundary has N+
                ! The region 2 is the region where the boundary has N-

                case (fbem_boundary_coupling_be_be)
                  select case (region(boundary(sb)%region(1))%type)
                    case (fbem_potential)
                      select case (region(boundary(sb)%region(2))%type)

                        ! ------------------------------------------------- !
                        ! BE (inviscid fluid) - BE (inviscid fluid)         !
                        ! ------------------------------------------------- !

                        case (fbem_potential)
                          !
                          ! Index of equations for each coordinate k:
                          ! node(sn)%row(1,1): BIE for region 1
                          ! node(sn)%row(1,2): BIE for region 2
                          !
                          ! Index of variables for each coordinate k:
                          ! node(sn)%col(1,1): p for region 1
                          ! node(sn)%col(2,1): Un for region 1
                          ! node(sn)%col(1,2): p for region 2 (not active since p_{region 1}=p_{region 2})
                          ! node(sn)%col(2,2): Un for region 2 (not active since Un_{region 1}=-Un_{region 2})
                          !
                          allocate (node(sn)%row(1,2))
                          allocate (node(sn)%col(2,2))
                          allocate (node(sn)%value_c(2,2))
                          allocate (node(sn)%dvda_c(2,2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Equation for region 1
                          node(sn)%row(1,1)=row
                          ! Equation for region 2
                          node(sn)%row(1,2)=row+1
                          ! Increment counter
                          row=row+2
                          ! Only the variables of region 1 are active
                          node(sn)%col(1,1)=col
                          node(sn)%col(2,1)=col+1
                          ! Increment counter
                          col=col+2

                        ! ------------------------------------------------------- !
                        ! BE (1: inviscid fluid) - BE (2: viscoelastic solid)     !
                        ! ------------------------------------------------------- !

                        case (fbem_viscoelastic)
                          !
                          ! Index of equations:
                          ! node(sn)%row(1,1): SBIE, HBIE or SBIE+beta*HBIE for region 1 (fluid)
                          ! node(sn)%row(k,2): SBIE, HBIE or SBIE+beta*HBIE for region 2 (viscoelastic)
                          !
                          ! Index of variables:
                          ! node(sn)%col(          1,1): p for region 1
                          ! node(sn)%col(          2,1): Un for region 1 (not active since Un_{region 1}=u_k_{region 2}·n_k_{region 1})
                          ! node(sn)%col(          k,2): u_k for region 2
                          ! node(sn)%col(problem%n+k,2): t_k for region 2 (not active since t_k_{region 2}=-p_{region 1}·n_k_{region 2})
                          !
                          allocate (node(sn)%row(problem%n,2))
                          allocate (node(sn)%col(2*problem%n,2))
                          allocate (node(sn)%value_c(2*problem%n,2))
                          allocate (node(sn)%dvda_c(2*problem%n,2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equations and column to the variables
                          ! Region 1
                          ! Equation for region 1
                          node(sn)%row(1,1)=row
                          ! Increment counter
                          row=row+1
                          ! The variable p of region 1 is active
                          node(sn)%col(1,1)=col
                          ! Increment counter
                          col=col+1
                          ! Region 2
                          do k=1,problem%n
                            ! Equation for region 2
                            node(sn)%row(k,2)=row
                            ! Increment counter
                            row=row+1
                            ! The variable u_k of region 2 is active
                            node(sn)%col(k,2)=col
                            ! Increment counter
                            col=col+1
                          end do

                        ! --------------------------------------------------- !
                        ! BE (1: inviscid fluid) - BE (2: poroelastic medium) !
                        ! ------------------------------------------.-------- !

                        case (fbem_poroelastic)
                          !
                          ! Index of equations:
                          ! node(sn)%row(1,1): SBIE, HBIE or SBIE+beta*HBIE for region 1 (fluid)
                          ! node(sn)%row(0,2): SBIE, HBIE or SBIE+beta*HBIE for region 2 (fluid phase)
                          ! node(sn)%row(k,2): SBIE, HBIE or SBIE+beta*HBIE for region 2 (solid skeleton)
                          !
                          ! Index of variables:
                          ! node(sn)%col(            1,1): p   for region 1
                          ! node(sn)%col(            2,1): Un  for region 1
                          ! node(sn)%col(            0,2): tau for region 2
                          ! node(sn)%col(            k,2): u_k for region 2
                          ! node(sn)%col(problem%n+1  ,2): Un  for region 2
                          ! node(sn)%col(problem%n+1+k,2): t_k for region 2
                          !
                          allocate (node(sn)%row(0:problem%n,2))
                          allocate (node(sn)%col(0:(1+2*problem%n),2))
                          allocate (node(sn)%value_c(0:(1+2*problem%n),2))
                          allocate (node(sn)%dvda_c(0:(1+2*problem%n),2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Region 1 (inviscid fluid)
                          ! BIE equation for region 1
                          node(sn)%row(1,1)=row
                          ! Increment counter
                          row=row+1
                          ! Switch depeding on the B.C.
                          select case (node(sn)%ctype(1,1))
                            ! 0. Perfectly permeable
                            case (0)
                              ! Both p and Un are inactive
                            ! 1. Perfectly impermeable
                            case (1)
                              ! The variable p is active
                              node(sn)%col(1,1)=col
                              ! Increment counter
                              col=col+1
                          end select
                          ! Region 2 (poroelastic medium)
                          do k=0,problem%n
                            ! BIE equation for region 2
                            node(sn)%row(k,2)=row
                            ! Increment counter
                            row=row+1
                            ! Switch depeding on the B.C.
                            select case (node(sn)%ctype(1,1))
                              ! 0. Perfectly permeable
                              case (0)
                                if (k.eq.0) then
                                  ! If k==0: tau and Un are active
                                  node(sn)%col(            k,2)=col
                                  node(sn)%col(problem%n+1+k,2)=col+1
                                  ! Increment counter
                                  col=col+2
                                else
                                  ! If k>=1: u_k is active
                                  node(sn)%col(            k,2)=col
                                  ! Increment counter
                                  col=col+1
                                end if
                              ! 1. Perfectly impermeable
                              case (1)
                                ! If k==0: tau is active
                                ! If k>=1: u_k is active
                                node(sn)%col(            k,2)=col
                                ! Increment counter
                                col=col+1
                            end select
                          end do

                      end select

                    case (fbem_viscoelastic)
                      select case (region(boundary(sb)%region(2))%type)

                        ! ---------------------------------------------------- !
                        ! BE (1: viscoelastic solid) - BE (2: inviscid fluid)  !
                        ! ---------------------------------------------------- !

                        case (fbem_potential)
                          !
                          ! Index of equations:
                          ! node(sn)%row(k,1): SBIE, HBIE or SBIE+beta*HBIE for region 1 (viscoelastic)
                          ! node(sn)%row(1,2): SBIE, HBIE or SBIE+beta*HBIE for region 2 (fluid)
                          !
                          ! Index of variables:
                          ! node(sn)%col(          k,1): u_k for region 1
                          ! node(sn)%col(problem%n+k,1): t_k for region 1 (not active since t_k_{region 1}=-p_{region 2}·n_k_{region 1})
                          ! node(sn)%col(          1,2): p for region 2
                          ! node(sn)%col(          2,2): Un for region 2 (not active since Un_{region 2}=u_k_{region 1}·n_k_{region 2})
                          !
                          allocate (node(sn)%row(problem%n,2))
                          allocate (node(sn)%col(2*problem%n,2))
                          allocate (node(sn)%value_c(2*problem%n,2))
                          allocate (node(sn)%dvda_c(2*problem%n,2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equations and column to the variables
                          ! Region 1
                          do k=1,problem%n
                            ! Equation for region 1
                            node(sn)%row(k,1)=row
                            ! Increment counter
                            row=row+1
                            ! The variable u_k of region 1 is active
                            node(sn)%col(k,1)=col
                            ! Increment counter
                            col=col+1
                          end do
                          ! Region 2
                          ! Equation for region 2
                          node(sn)%row(1,2)=row
                          ! Increment counter
                          row=row+1
                          ! The variable p of region 2 is active
                          node(sn)%col(1,2)=col
                          ! Increment counter
                          col=col+1

                        ! ------------------------------------------------- !
                        ! BE (viscoelastic solid) - BE (viscoelastic solid) !
                        ! ------------------------------------------------- !

                        case (fbem_viscoelastic)
                          !
                          ! Index of equations for each coordinate k:
                          ! node(sn)%row(k,1): SBIE, HBIE or SBIE+beta*HBIE for region 1
                          ! node(sn)%row(k,2): SBIE, HBIE or SBIE+beta*HBIE for region 2
                          ! Index of additional equations for each coordinate k:
                          ! node(sn)%row(k,0): compatibility or equilibrium equation
                          !
                          ! Index of variables for each coordinate k:
                          ! node(sn)%col(          k,1): u_k for region 1
                          ! node(sn)%col(problem%n+k,1): t_k for region 1
                          ! node(sn)%col(          k,2): u_k for region 2
                          ! node(sn)%col(problem%n+k,2): t_k for region 2
                          !
                          allocate (node(sn)%row(problem%n,0:2))
                          allocate (node(sn)%col(2*problem%n,2))
                          allocate (node(sn)%value_c(2*problem%n,2))
                          allocate (node(sn)%dvda_c(2*problem%n,2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equations and column to the variables
                          do k=1,problem%n
                            ! BIE for region 1
                            node(sn)%row(k,1)=row
                            ! BIE for region 2
                            node(sn)%row(k,2)=row+1
                            ! Increment counter
                            row=row+2
                            select case (node(sn)%ctype(k,1))
                              ! Global axes
                              ! 0 - perfect bonding - u_k^{1} and t_k^{1} are active since:
                              ! t_k^{(2)}=-t_k^{(1)}
                              ! u_k^{(2)}= u_k^{(1)}
                              case (0)
                                ! Active variables
                                node(sn)%col(          k,1)=col
                                node(sn)%col(problem%n+k,1)=col+1
                                col=col+2
                              ! Global axes
                              ! 1 - perfect debonding - u_k^{1} and u_k^{2} are active since:
                              ! t_k^{(1)}=0
                              ! t_k^{(2)}=0
                              ! 2 - partial bonding   - u_k^{1} and u_k^{2} are active since:
                              ! t_k^{(1)}=(u_k^{2}-u_k^{1})*KT
                              ! t_k^{(2)}=-(u_k^{2}-u_k^{1})*KT
                              case (1,2)
                                ! Active variables
                                node(sn)%col(          k,1)=col
                                node(sn)%col(          k,2)=col+1
                                col=col+2
                              ! Local axes
                              ! In the three cases, u_k^{1}, u_k^{2} and t_k^{1} are active, and t_k^{(2)}=-t_k^{(1)}.
                              ! Each B.C. needs an additional equation:
                              ! 3 - perfect bonding  : (u^{2}-u^{1})·l_k^(1)=0
                              ! 4 - perfect debonding: t^(1)·l_k^(1)=0
                              ! 5 - partial bonding  : (u^{2}-u^{1})·l_k^(1)*KT=t^(1)·l_k^(1)
                              case (3,4,5)
                                ! Active variables
                                node(sn)%col(          k,1)=col
                                node(sn)%col(problem%n+k,1)=col+1
                                node(sn)%col(          k,2)=col+2
                                col=col+3
                                ! Additional equation: (u^{2}-u^{1})·l_k^(1)=0
                                node(sn)%row(k,0)=row
                                row=row+1
                            end select
                          end do

                        ! ------------------------------------------------------- !
                        ! BE (1: viscoelastic solid) - BE (2: poroelastic medium) !
                        ! ------------------------------------------------------- !

                        case (fbem_poroelastic)
                          !
                          ! Index of equations for each coordinate k:
                          ! node(sn)%row(k,1): SBIE, HBIE or SBIE+beta*HBIE for region 1 (solid)
                          ! node(sn)%row(0,2): SBIE, HBIE or SBIE+beta*HBIE for region 2 (fluid phase)
                          ! node(sn)%row(k,2): SBIE, HBIE or SBIE+beta*HBIE for region 2 (solid skeleton)
                          ! Index of additional equations for each coordinate k:
                          ! node(sn)%row(k,0): compatibility or equilibrium equation
                          !
                          ! Index of variables for each coordinate k:
                          ! node(sn)%col(            k,1): u_k for region 1
                          ! node(sn)%col(  problem%n+k,1): t_k for region 1
                          ! node(sn)%col(            0,2): tau for region 2
                          ! node(sn)%col(            k,2): u_k for region 2
                          ! node(sn)%col(problem%n+1  ,2): Un  for region 2
                          ! node(sn)%col(problem%n+1+k,2): t_k for region 2
                          !
                          allocate (node(sn)%row(0:problem%n,0:2))
                          allocate (node(sn)%col(0:(1+2*problem%n),2))
                          allocate (node(sn)%value_c(0:(1+2*problem%n),2))
                          allocate (node(sn)%dvda_c(0:(1+2*problem%n),2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equations and column to the variables
                          !
                          ! Fluid phase of poroelastic region
                          !
                          ! BIE for fluid phase of region 2
                          node(sn)%row(0,2)=row
                          row=row+1
                          ! tau^(2) is always active
                          node(sn)%col(0,2)=col
                          col=col+1
                          ! Note: Un^(2)=-u^(2)·n^(1) always
                          !
                          ! Solid skeleton
                          !
                          do k=1,problem%n
                            ! BIE for region 1
                            node(sn)%row(k,1)=row
                            ! BIE for region 2
                            node(sn)%row(k,2)=row+1
                            ! Increment counter
                            row=row+2
                            ! Depending on the B.C. of each coordinate of each node
                            select case (node(sn)%ctype(k,1))
                              ! Global axes
                              ! 0 - perfect bonding - u_k^{2} and t_k^{2} are active since:
                              ! t_k^{(1)}=-t_k^{(2)}+tau^(2)n_k^(1)
                              ! u_k^{(1)}= u_k^{(2)}
                              case (0)
                                ! Active variables
                                node(sn)%col(            k,2)=col
                                node(sn)%col(problem%n+1+k,2)=col+1
                                col=col+2
                              ! Global axes
                              ! 1 - perfect debonding - u_k^{1} and u_k^{2} are active since:
                              ! t_k^{(1)}=0
                              ! t_k^{(2)}=tau^(2)n_k^(1)
                              ! 2 - partial bonding   - u_k^{1} and u_k^{2} are active since:
                              ! t_k^{(1)}=(u_k^{2}-u_k^{1})*KT
                              ! t_k^{(2)}=-(u_k^{2}-u_k^{1})*KT+tau^(2)n_k^(1)
                              case (1,2)
                                ! Active variables
                                node(sn)%col(          k,1)=col
                                node(sn)%col(          k,2)=col+1
                                col=col+2
                              ! Local axes
                              ! In all cases, u_k^{1}, u_k^{2} and t_k^{1} are active, and t_k^{(2)}=-t_k^{(1)}+tau^(2)n_k^(1).
                              ! Each B.C. needs an additional equation:
                              ! 3 - perfect bonding  : (u^{2}-u^{1})·l_k^(1)=0
                              ! 4 - perfect debonding: t^(1)·l_k^(1)=0
                              ! 5 - partial bonding  : (u^{2}-u^{1})·l_k^(1)*KT=t^(1)·l_k^(1)
                              case (3,4,5)
                                ! Active variables
                                node(sn)%col(          k,1)=col
                                node(sn)%col(problem%n+k,1)=col+1
                                node(sn)%col(          k,2)=col+2
                                col=col+3
                                ! Additional equation
                                node(sn)%row(k,0)=row
                                row=row+1
                            end select
                          end do

                      end select

                    case (fbem_poroelastic)
                      select case (region(boundary(sb)%region(2))%type)

                        ! --------------------------------------------------- !
                        ! BE (1: poroelastic medium) - BE (2: inviscid fluid) !
                        ! --------------------------------------------------- !

                        case (fbem_potential)
                          !
                          ! Index of equations:
                          ! node(sn)%row(0,1): SBIE, HBIE or SBIE+beta*HBIE for region 1 (fluid phase)
                          ! node(sn)%row(k,1): SBIE, HBIE or SBIE+beta*HBIE for region 1 (solid skeleton)
                          ! node(sn)%row(1,2): SBIE, HBIE or SBIE+beta*HBIE for region 2 (fluid)
                          !
                          ! Index of variables:
                          ! node(sn)%col(            0,1): tau for region 1
                          ! node(sn)%col(            k,1): u_k for region 1
                          ! node(sn)%col(problem%n+1  ,1): Un  for region 1
                          ! node(sn)%col(problem%n+1+k,1): t_k for region 1
                          ! node(sn)%col(            1,2): p   for region 2
                          ! node(sn)%col(            2,2): Un  for region 2
                          !
                          allocate (node(sn)%row(0:problem%n,2))
                          allocate (node(sn)%col(0:(1+2*problem%n),2))
                          allocate (node(sn)%value_c(0:(1+2*problem%n),2))
                          allocate (node(sn)%dvda_c(0:(1+2*problem%n),2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Region 2 (inviscid fluid)
                          ! BIE equation for region 2
                          node(sn)%row(1,2)=row
                          ! Increment counter
                          row=row+1
                          ! Switch depeding on the B.C.
                          select case (node(sn)%ctype(1,1))
                            ! 0. Perfectly permeable
                            case (0)
                              ! Both p and Un are inactive
                            ! 1. Perfectly impermeable
                            case (1)
                              ! The variable p is active
                              node(sn)%col(1,2)=col
                              ! Increment counter
                              col=col+1
                          end select
                          ! Region 1 (poroelastic medium)
                          do k=0,problem%n
                            ! BIE equation for region 1
                            node(sn)%row(k,1)=row
                            ! Increment counter
                            row=row+1
                            ! Switch depeding on the B.C.
                            select case (node(sn)%ctype(1,1))
                              ! 0. Perfectly permeable
                              case (0)
                                if (k.eq.0) then
                                  ! If k==0: tau and Un are active
                                  node(sn)%col(            k,1)=col
                                  node(sn)%col(problem%n+1+k,1)=col+1
                                  ! Increment counter
                                  col=col+2
                                else
                                  ! If k>=1: u_k is active
                                  node(sn)%col(            k,1)=col
                                  ! Increment counter
                                  col=col+1
                                end if
                              ! 1. Perfectly impermeable
                              case (1)
                                ! If k==0: tau is active
                                ! If k>=1: u_k is active
                                node(sn)%col(            k,1)=col
                                ! Increment counter
                                col=col+1
                            end select
                          end do

                        ! ------------------------------------------------------- !
                        ! BE (1: poroelastic medium) - BE (2: viscoelastic solid) !
                        ! ------------------------------------------------------- !

                        case (fbem_viscoelastic)
                          !
                          ! Index of equations for each coordinate k:
                          ! node(sn)%row(0,1): SBIE, HBIE or SBIE+beta*HBIE for region 1 (fluid phase)
                          ! node(sn)%row(k,1): SBIE, HBIE or SBIE+beta*HBIE for region 1 (solid skeleton)
                          ! node(sn)%row(k,2): SBIE, HBIE or SBIE+beta*HBIE for region 2 (solid)
                          ! Index of additional equations for each coordinate k:
                          ! node(sn)%row(k,0): compatibility or equilibrium equation
                          !
                          ! Index of variables for each coordinate k:
                          ! node(sn)%col(            0,1): tau for region 1
                          ! node(sn)%col(            k,1): u_k for region 1
                          ! node(sn)%col(problem%n+1  ,1): Un  for region 1
                          ! node(sn)%col(problem%n+1+k,1): t_k for region 1
                          ! node(sn)%col(            k,2): u_k for region 2
                          ! node(sn)%col(  problem%n+k,2): t_k for region 2
                          !
                          allocate (node(sn)%row(0:problem%n,0:2))
                          allocate (node(sn)%col(0:(1+2*problem%n),2))
                          allocate (node(sn)%value_c(0:(1+2*problem%n),2))
                          allocate (node(sn)%dvda_c(0:(1+2*problem%n),2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equations and column to the variables
                          !
                          ! Fluid phase of poroelastic region
                          !
                          ! BIE for fluid phase of region 1
                          node(sn)%row(0,1)=row
                          row=row+1
                          ! tau^(2) is always active
                          node(sn)%col(0,1)=col
                          col=col+1
                          ! Note: Un^(1)=u^(1)·n^(1) always
                          !
                          ! Solid skeleton
                          !
                          do k=1,problem%n
                            ! BIE for region 1
                            node(sn)%row(k,1)=row
                            ! BIE for region 2
                            node(sn)%row(k,2)=row+1
                            ! Increment counter
                            row=row+2
                            ! Depending on the B.C. of each coordinate of each node
                            select case (node(sn)%ctype(k,1))
                              ! Global axes
                              ! 0 - perfect bonding - u_k^{1} and t_k^{1} are active since:
                              ! t_k^{(2)}=-t_k^{(1)}-tau^(1)n_k^(1)
                              ! u_k^{(2)}= u_k^{(1)}
                              case (0)
                                ! Active variables
                                node(sn)%col(            k,1)=col
                                node(sn)%col(problem%n+1+k,1)=col+1
                                col=col+2
                              ! Global axes
                              ! 1 - perfect debonding - u_k^{1} and u_k^{2} are active since:
                              ! t_k^{(1)}=-tau^(1)n_k^(1)
                              ! t_k^{(2)}=0
                              ! 2 - partial bonding   - u_k^{1} and u_k^{2} are active since:
                              ! t_k^{(1)}=-(u_k^{1}-u_k^{2})*KT-tau^(1)n_k^(1)
                              ! t_k^{(2)}=(u_k^{1}-u_k^{2})*KT
                              case (1,2)
                                ! Active variables
                                node(sn)%col(          k,1)=col
                                node(sn)%col(          k,2)=col+1
                                col=col+2
                              ! Local axes
                              ! In all cases, u_k^{1}, u_k^{2} and t_k^{2} are active, and t_k^{(1)}=-t_k^{(2)}-tau^(1)n_k^(1).
                              ! Each B.C. needs an additional equation:
                              ! 3 - perfect bonding  : (u^{1}-u^{2})·l_k^(1)=0
                              ! 4 - perfect debonding: t^(2)·l_k^(1)=0
                              ! 5 - partial bonding  : (u^{1}-u^{2})·l_k^(1)*KT=t^(2)·l_k^(1)
                              case (3,4,5)
                                ! Active variables
                                node(sn)%col(          k,2)=col
                                node(sn)%col(problem%n+k,2)=col+1
                                node(sn)%col(          k,1)=col+2
                                col=col+3
                                ! Additional equation
                                node(sn)%row(k,0)=row
                                row=row+1
                            end select
                          end do

                        ! ------------------------------------------------------- !
                        ! BE (1: poroelastic medium) - BE (2: poroelastic medium) !
                        ! ------------------------------------------------------- !

                        case (fbem_poroelastic)
                          !
                          ! Index of equations for each coordinate k:
                          ! node(sn)%row(0,1): SBIE, HBIE or SBIE+beta*HBIE for region 1 (fluid)
                          ! node(sn)%row(k,1): SBIE, HBIE or SBIE+beta*HBIE for region 1 (solid skeleton)
                          ! node(sn)%row(0,2): SBIE, HBIE or SBIE+beta*HBIE for region 2 (fluid)
                          ! node(sn)%row(k,2): SBIE, HBIE or SBIE+beta*HBIE for region 2 (solid skeleton)
                          !
                          ! Index of variables:
                          !
                          ! B.C. 0 (perfectly permeable):
                          ! node(sn)%col(            0,1): tau for region 1 (active)
                          ! node(sn)%col(            k,1): u_k for region 1 (active)
                          ! node(sn)%col(problem%n+1  ,1): Un  for region 1 (active)
                          ! node(sn)%col(problem%n+1+k,1): t_k for region 1 (active)
                          ! node(sn)%col(            0,2): tau for region 2 (not active since tau^(2)=phi^(2)/phi^(1)*tau^(1))
                          ! node(sn)%col(            k,2): u_k for region 2 (not active since u_k^(2)=u_k^(1))
                          ! node(sn)%col(problem%n+1  ,2): Un  for region 2 (not active since Un^(2)=-phi^(1)/phi^(2)*Un^(1)-(1-phi^(1)/phi^(2))*u^(1)·n^(1))
                          ! node(sn)%col(problem%n+1+k,2): t_k for region 2 (not active since t_k^(2)=-t_k^(1)-(1-phi^(2)/phi^(1))*n_k^(1)*tau^(1))
                          !
                          ! B.C. 1 (perfectly impermeable):
                          ! node(sn)%col(            0,1): tau for region 1 (active)
                          ! node(sn)%col(            k,1): u_k for region 1 (active)
                          ! node(sn)%col(problem%n+1  ,1): Un  for region 1 (not active since Un^(1)=u^(1)·n^(1))
                          ! node(sn)%col(problem%n+1+k,1): t_k for region 1 (active)
                          ! node(sn)%col(            0,2): tau for region 2 (active)
                          ! node(sn)%col(            k,2): u_k for region 2 (not active since u_k^(2)=u_k^(1))
                          ! node(sn)%col(problem%n+1  ,2): Un  for region 2 (not active since Un^(2)=-u^(1)·n^(1))
                          ! node(sn)%col(problem%n+1+k,2): t_k for region 2 (not active since t_k^(2)=-t_k^(1)-tau^(1)n_k^(1)+tau^(2)n_k^(1))
                          !
                          ! B.C. 2 (partially permeable):
                          ! node(sn)%col(            0,1): tau for region 1 (active)
                          ! node(sn)%col(            k,1): u_k for region 1 (active)
                          ! node(sn)%col(problem%n+1  ,1): Un  for region 1 (active)
                          ! node(sn)%col(problem%n+1+k,1): t_k for region 1 (active)
                          ! node(sn)%col(            0,2): tau for region 2 (not active since tau^(2)=k*phi1*phi2*i*omega*(Un^(1)-u_j^(1)*n_j^(1))+phi^(2)/phi^(1)*tau^(1))
                          ! node(sn)%col(            k,2): u_k for region 2 (not active since u_k^(2)=u_k^(1))
                          ! node(sn)%col(problem%n+1  ,2): Un  for region 2 (not active since Un^(2)=-phi^(1)/phi^(2)*Un^(1)-(1-phi^(1)/phi^(2))*u^(1)·n^(1))
                          ! node(sn)%col(problem%n+1+k,2): t_k for region 2 (not active since t_k^(2)=-t_k^(1)+[k*phi1*phi2*i*omega*(Un^(1)-u_j^(1)*n_j^(1))-(1-phi^(2)/phi^(1))*tau^(1)]*n_k^(1)
                          !
                          allocate (node(sn)%row(0:problem%n,2))
                          allocate (node(sn)%col(0:(1+2*problem%n),2))
                          allocate (node(sn)%value_c(0:(1+2*problem%n),2))
                          allocate (node(sn)%dvda_c(0:(1+2*problem%n),2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equations and column to the variables
                          do k=0,problem%n
                            ! BIE for region 1
                            node(sn)%row(k,1)=row
                            ! BIE for region 2
                            node(sn)%row(k,2)=row+1
                            ! Increment counter
                            row=row+2
                            ! Active variables
                            select case (node(sn)%ctype(1,1))
                              case (0,2)
                                ! Only the variables of region 1 are active
                                node(sn)%col(            k,1)=col
                                node(sn)%col(problem%n+1+k,1)=col+1
                              case (1)
                                ! Fluid phase variables
                                if (k.eq.0) then
                                  ! tau for both regions are active
                                  node(sn)%col(            k,1)=col
                                  node(sn)%col(            k,2)=col+1
                                ! Solid skeleton variables
                                else
                                  ! u_k and t_k of region 1 variables are active
                                  node(sn)%col(            k,1)=col
                                  node(sn)%col(problem%n+1+k,1)=col+1
                                end if
                            end select
                            ! Increment counter
                            col=col+2
                          end do

                      end select

                  end select

                ! ================================================================================================================ !
                ! BE-FE BOUNDARY                                                                                                   !
                ! ================================================================================================================ !

                case (fbem_boundary_coupling_be_fe)
                  select case (boundary(sb)%class)

                    ! ================= !
                    ! ORDINARY BOUNDARY !
                    ! ================= !

                    case (fbem_boundary_class_ordinary)
                      select case (region(kr)%type)

                        ! -------------- !
                        ! INVISCID FLUID !
                        ! -------------- !

                        case (fbem_potential)
                          !
                          ! Index of equations:
                          ! node(sn)%row(1,1): SBIE, HBIE or SBIE+beta*HBIE
                          !
                          ! Index of variables:
                          ! node(sn)%col(1,1): p
                          ! node(sn)%col(2,1): Un (inactive since Un=u^(fe)·n, impermeable condition)
                          !
                          allocate (node(sn)%row(1,1))
                          allocate (node(sn)%col(2,1))
                          allocate (node(sn)%value_c(2,1))
                          allocate (node(sn)%dvda_c(2,1,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equation and column to the variables
                          ! Equation
                          node(sn)%row(1,1)=row
                          ! Increment counter
                          row=row+1
                          ! p unknown
                          node(sn)%col(1,1)=col
                          ! Increment counter
                          col=col+1

                        ! ------------------ !
                        ! VISCOELASTIC SOLID !
                        ! ------------------ !

                        case (fbem_viscoelastic)
                          !
                          ! Index of equations for each coordinate k:
                          ! node(sn)%row(k,1): SBIE, HBIE or SBIE+beta*HBIE
                          !
                          ! Index of variables for each coordinate k:
                          ! node(sn)%col(          k,1): u_k (inactive since u_k^(E)=u_k^(FE))
                          ! node(sn)%col(problem%n+k,1): t_k
                          !
                          allocate (node(sn)%row(problem%n,1))
                          allocate (node(sn)%col(2*problem%n,1))
                          allocate (node(sn)%value_c(2*problem%n,1))
                          allocate (node(sn)%dvda_c(2*problem%n,1,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equation and column to the variables
                          do k=1,problem%n
                            ! Equation
                            node(sn)%row(k,1)=row
                            ! Increment counter
                            row=row+1
                            ! t_k is always unknown
                            node(sn)%col(problem%n+k,1)=col
                            ! Increment counter
                            col=col+1
                            ! Perfect bonding: u_k^(E)=u_k^(FE)
                          end do

                        ! ------------------ !
                        ! POROELASTIC MEDIUM !
                        ! ------------------ !

                        case (fbem_poroelastic)
                          !
                          ! Index of equations:
                          ! node(sn)%row(0,1): SBIE, HBIE or SBIE+beta*HBIE (fluid phase)
                          ! node(sn)%row(k,1): SBIE, HBIE or SBIE+beta*HBIE (solid skeleton)
                          !
                          ! Index of variables:
                          ! node(sn)%col(            0,1): tau
                          ! node(sn)%col(            k,1): u_k
                          ! node(sn)%col(problem%n+1  ,1): Un (inactive since Un=u·n, impermeable condition)
                          ! node(sn)%col(problem%n+1+k,1): t_k
                          !
                          allocate (node(sn)%row    (0:problem%n      ,1))
                          allocate (node(sn)%col    (0:(1+2*problem%n),1))
                          allocate (node(sn)%value_c(0:(1+2*problem%n),1))
                          allocate (node(sn)%dvda_c(0:(1+2*problem%n),1,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equation and column to the variables
                          ! Fluid phase
                          ! BIE equation
                          node(sn)%row(0,1)=row
                          ! Increment counter
                          row=row+1
                          ! tau unknown
                          node(sn)%col(0,1)=col
                          ! Increment counter
                          col=col+1
                          ! Solid skeleton
                          do k=1,problem%n
                            ! BIE equation
                            node(sn)%row(k,1)=row
                            ! Increment counter
                            row=row+1
                            ! t_k is always unknown
                            node(sn)%col(problem%n+1+k,1)=col
                            ! Increment counter
                            col=col+1
                            ! Perfect bonding: u_k^(E)=u_k^(FE)
                          end do

                      end select

                    ! =================== !
                    ! CRACK-LIKE BOUNDARY !
                    ! =================== !

                    case (fbem_boundary_class_cracklike)
                      select case (region(kr)%type)

                        ! -------------- !
                        ! INVISCID FLUID !
                        ! -------------- !

                        case (fbem_potential)
                          !
                          ! Index of equations:
                          ! node(sn)%row(1,1): SBIE
                          ! node(sn)%row(1,2): HBIE
                          !
                          ! Index of variables:
                          ! node(sn)%col(1,1): p^+
                          ! node(sn)%col(2,1): Un^+ (inactive since Un^+= u^(fe)·n)
                          ! node(sn)%col(1,2): p^-
                          ! node(sn)%col(2,2): Un^- (inactive since Un^-=-u^(fe)·n)
                          !
                          allocate (node(sn)%row(1,2))
                          allocate (node(sn)%col(2,2))
                          allocate (node(sn)%value_c(2,2))
                          allocate (node(sn)%dvda_c(2,2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equation and column to the variables
                          ! SBIE equation
                          node(sn)%row(1,1)=row
                          ! HBIE equation
                          node(sn)%row(1,2)=row+1
                          ! Increment counter
                          row=row+2
                          ! p^+ unknown
                          node(sn)%col(1,1)=col
                          ! p^- unknown
                          node(sn)%col(1,2)=col+1
                          ! Increment counter
                          col=col+2

                        ! ------------------ !
                        ! VISCOELASTIC SOLID !
                        ! ------------------ !

                        case (fbem_viscoelastic)

!                          ! -------------------------------------------------------------------------------------- !
!                          ! SBIE (solve for Σt = t^+ + t^-) -> HBIE (calculate t^+ and t^- ) (for perfect bonding) !
!                          ! -------------------------------------------------------------------------------------- !

!                          !
!                          ! Index of equations for each coordinate k:
!                          ! node(sn)%row(k,1): SBIE
!                          !
!                          ! Index of variables for each coordinate k:
!                          ! node(sn)%col(          k,1): inactive
!                          ! node(sn)%col(problem%n+k,1): Σt_k = t_k^+ + t_k^-
!                          !
!                          ! Index of variable values for each coordinate k:
!                          !
!                          ! node(sn)%value_c(          k,1): u_k^+
!                          ! node(sn)%value_c(problem%n+k,1): t_k^+
!                          ! node(sn)%value_c(          k,2): u_k^-
!                          ! node(sn)%value_c(problem%n+k,2): t_k^-
!                          !
!                          ! Index of variable sensitivity values for each coordinate k:
!                          !
!                          ! node(sn)%dvda_c(          k,1,j): ∂u_k^+/∂a_j
!                          ! node(sn)%dvda_c(problem%n+k,1,j): ∂t_k^+/∂a_j
!                          ! node(sn)%dvda_c(          k,2,j): ∂u_k^-/∂a_j
!                          ! node(sn)%dvda_c(problem%n+k,2,j): ∂t_k^-/∂a_j
!                          !
!                          allocate (node(sn)%row(problem%n,1))
!                          allocate (node(sn)%col(2*problem%n,1))
!                          allocate (node(sn)%value_c(2*problem%n,2))
!                          allocate (node(sn)%dvda_c(2*problem%n,2,problem%n_designvariables))
!                          ! Initialize
!                          node(sn)%row=0
!                          node(sn)%col=0
!                          ! Assign row to the equation and column to the variables
!                          do k=1,problem%n
!                            ! SBIE
!                            node(sn)%row(k,1)=row
!                            row=row+1
!                            ! Σt_k
!                            node(sn)%col(problem%n+k,1)=col
!                            col=col+1
!                          end do

                          ! ------------------------------------------------------------- !
                          ! General model using Dual BEM (for general contact conditions) !
                          ! ------------------------------------------------------------- !

                          !
                          ! Index of equations for each coordinate k:
                          ! node(sn)%row(k,1): SBIE
                          ! node(sn)%row(k,2): HBIE
                          !
                          ! Index of variables for each coordinate k:
                          ! node(sn)%col(          k,1): u_k for face +
                          ! node(sn)%col(problem%n+k,1): t_k for face +
                          ! node(sn)%col(          k,2): u_k for face -
                          ! node(sn)%col(problem%n+k,2): t_k for face -
                          !
                          allocate (node(sn)%row(problem%n,2))
                          allocate (node(sn)%col(2*problem%n,2))
                          allocate (node(sn)%value_c(2*problem%n,2))
                          allocate (node(sn)%dvda_c(2*problem%n,2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equation and column to the variables
                          do k=1,problem%n
                            ! Equation for SBIE
                            node(sn)%row(k,1)=row
                            row=row+1
                            ! Equation for HBIE
                            node(sn)%row(k,2)=row
                            row=row+1
                            ! Assign values depending on the boundary condition.
                            ! t_k for face + is always unknown
                            node(sn)%col(problem%n+k,1)=col
                            col=col+1
                            ! t_k for face - is always unknown
                            node(sn)%col(problem%n+k,2)=col
                            col=col+1
                            ! Perfect bonding: u_k^+(E)=u_k^-(E)=u_k^(FE)
                          end do

                        ! ------------------ !
                        ! POROELASTIC MEDIUM !
                        ! ------------------ !

                        case (fbem_poroelastic)
                          !
                          ! Index of equations:
                          ! node(sn)%row(0,1): SBIE (fluid phase)
                          ! node(sn)%row(k,1): SBIE (solid skeleton)
                          ! node(sn)%row(0,2): HBIE (fluid phase)
                          ! node(sn)%row(k,2): HBIE (solid skeleton)
                          !
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
                          allocate (node(sn)%row    (0:problem%n      ,2))
                          allocate (node(sn)%col    (0:(1+2*problem%n),2))
                          allocate (node(sn)%value_c(0:(1+2*problem%n),2))
                          allocate (node(sn)%dvda_c(0:(1+2*problem%n),2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equation and column to the variables
                          ! Fluid phase
                          ! SBIE equation
                          node(sn)%row(0,1)=row
                          ! HBIE equation
                          node(sn)%row(0,2)=row+1
                          ! Increment counter
                          row=row+2
                          ! tau^+ unknown
                          node(sn)%col(0,1)=col
                          ! tau^- unknown
                          node(sn)%col(0,2)=col+1
                          ! Increment counter
                          col=col+2
                          ! Solid skeleton
                          do k=1,problem%n
                            ! SBIE equation
                            node(sn)%row(k,1)=row
                            ! HBIE equation
                            node(sn)%row(k,2)=row+1
                            ! Increment counter
                            row=row+2
                            ! t_k^+ is always unknown
                            node(sn)%col(problem%n+1+k,1)=col
                            ! t_k^- is always unknown
                            node(sn)%col(problem%n+1+k,2)=col+1
                            ! Increment counter
                            col=col+2
                            ! Perfect bonding: u_k^(E+)=u_k^(E-)=u_k^(FE)
                          end do

                      end select
                  end select

                ! ================================================================================================================ !
                ! BE-FE-BE BOUNDARY                                                                                                !
                ! ================================================================================================================ !
                !
                ! The region 1 is the region where the boundary has N+
                ! The region 2 is the region where the boundary has N-

                case (fbem_boundary_coupling_be_fe_be)
                  select case (region(boundary(sb)%region(1))%type)
                    case (fbem_potential)
                      select case (region(boundary(sb)%region(2))%type)

                        ! ------------------------------------------------- !
                        ! BE (inviscid fluid) - BE (inviscid fluid)         !
                        ! ------------------------------------------------- !

                        case (fbem_potential)
                          !
                          ! Index of equations:
                          ! node(sn)%row(1,1): BIE for region 1
                          ! node(sn)%row(1,2): BIE for region 2
                          !
                          ! Index of variables:
                          ! node(sn)%col(1,1): p for region 1
                          ! node(sn)%col(2,1): Un for region 1
                          ! node(sn)%col(1,2): p for region 2
                          ! node(sn)%col(2,2): Un for region 2
                          !
                          allocate (node(sn)%row(1,2))
                          allocate (node(sn)%col(2,2))
                          allocate (node(sn)%value_c(2,2))
                          allocate (node(sn)%dvda_c(2,2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equation and column to the variables
                          ! BIE equation for region 1
                          node(sn)%row(1,1)=row
                          ! BIE equation for region 2
                          node(sn)%row(1,2)=row+1
                          ! Increment counter
                          row=row+2
                          ! p for region 1 is always unknown
                          node(sn)%col(1,1)=col
                          ! p for region 2 is always unknown
                          node(sn)%col(1,2)=col+1
                          ! Increment counter
                          col=col+2
                          ! Perfect bonding: Un^(1)=-Un^(2)=u_k^(FE)n_k^(1)

                        ! ------------------------------------------------------- !
                        ! BE (1: inviscid fluid) - BE (2: viscoelastic solid)     !
                        ! ------------------------------------------------------- !

                        case (fbem_viscoelastic)
                          !
                          ! Index of equations:
                          ! node(sn)%row(1,1): BIE for region 1 (fluid)
                          ! node(sn)%row(k,2): BIE for region 2 (viscoelastic)
                          !
                          ! Index of variables:
                          ! node(sn)%col(          1,1): p   for region 1
                          ! node(sn)%col(          2,1): Un  for region 1
                          ! node(sn)%col(          k,2): u_k for region 2
                          ! node(sn)%col(problem%n+k,2): t_k for region 2
                          !
                          allocate (node(sn)%row(problem%n,2))
                          allocate (node(sn)%col(2*problem%n,2))
                          allocate (node(sn)%value_c(2*problem%n,2))
                          allocate (node(sn)%dvda_c(2*problem%n,2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equations and column to the variables
                          ! Region 1
                          ! Equation for region 1
                          node(sn)%row(1,1)=row
                          ! Increment counter
                          row=row+1
                          ! The variable p of region 1 is active
                          node(sn)%col(1,1)=col
                          ! Increment counter
                          col=col+1
                          ! Region 2
                          do k=1,problem%n
                            ! Equation for region 2
                            node(sn)%row(k,2)=row
                            ! Increment counter
                            row=row+1
                            ! The variable t_k of region 2
                            node(sn)%col(problem%n+k,2)=col
                            ! Increment counter
                            col=col+1
                          end do
                          ! Perfect bonding: Un^(1)=u_k^(FE)n_k^(1), u_k^(2)=u_k^(FE)

                        ! --------------------------------------------------- !
                        ! BE (1: inviscid fluid) - BE (2: poroelastic medium) !
                        ! --------------------------------------------------- !

                        case (fbem_poroelastic)
                          !
                          ! Index of equations for each coordinate k:
                          ! node(sn)%row(1,1): BIE for region 1 (inviscid fluid)
                          ! node(sn)%row(0,2): BIE region 2 (fluid)
                          ! node(sn)%row(k,2): BIE region 2 (solid skeleton)
                          !
                          ! Index of variables:
                          ! node(sn)%col(            1,1): p   for region 1 (active)
                          ! node(sn)%col(            2,1): Un  for region 1
                          ! node(sn)%col(            0,2): tau for region 2 (active)
                          ! node(sn)%col(            k,2): u_k for region 2
                          ! node(sn)%col(problem%n+1  ,2): Un  for region 2
                          ! node(sn)%col(problem%n+1+k,2): t_k for region 2 (active)
                          !
                          allocate (node(sn)%row    (0:problem%n      ,2))
                          allocate (node(sn)%col    (0:(1+2*problem%n),2))
                          allocate (node(sn)%value_c(0:(1+2*problem%n),2))
                          allocate (node(sn)%dvda_c (0:(1+2*problem%n),2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equation and column to the variables
                          ! Region 1
                          ! Equation for region 1
                          node(sn)%row(1,1)=row
                          ! Increment counter
                          row=row+1
                          ! The variable p of region 1 is active
                          node(sn)%col(1,1)=col
                          ! Increment counter
                          col=col+1
                          ! Region 2
                          ! Fluid phase
                          ! BIE equation for region 2
                          node(sn)%row(0,2)=row
                          ! Increment counter
                          row=row+1
                          ! tau^(2) unknown
                          node(sn)%col(0,2)=col
                          ! Increment counter
                          col=col+1
                          ! Solid skeleton
                          do k=1,problem%n
                            ! BIE equation for region 2
                            node(sn)%row(k,2)=row
                            ! Increment counter
                            row=row+1
                            ! t_k^(2) unknown
                            node(sn)%col(problem%n+1+k,2)=col
                            ! Increment counter
                            col=col+1
                          end do
                          ! Perfect bonding: Un^(1)=-Un^(2)=u_k^(FE)n_k^(1), u_k^(2)=u_k^(FE)

                      end select

                    case (fbem_viscoelastic)
                      select case (region(boundary(sb)%region(2))%type)

                        ! ---------------------------------------------------- !
                        ! BE (1: viscoelastic solid) - BE (2: inviscid fluid)  !
                        ! ---------------------------------------------------- !

                        case (fbem_potential)
                          !
                          ! Index of equations:
                          ! node(sn)%row(k,1): BIE for region 1 (viscoelastic)
                          ! node(sn)%row(1,2): BIE for region 2 (fluid)
                          !
                          ! Index of variables:
                          ! node(sn)%col(          k,1): u_k for region 1
                          ! node(sn)%col(problem%n+k,1): t_k for region 1
                          ! node(sn)%col(          1,2): p   for region 2
                          ! node(sn)%col(          2,2): Un  for region 2
                          !
                          allocate (node(sn)%row(problem%n,2))
                          allocate (node(sn)%col(2*problem%n,2))
                          allocate (node(sn)%value_c(2*problem%n,2))
                          allocate (node(sn)%dvda_c(2*problem%n,2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equations and column to the variables
                          ! Region 1
                          do k=1,problem%n
                            ! Equation for region 1
                            node(sn)%row(k,1)=row
                            ! Increment counter
                            row=row+1
                            ! The variable t_k of region 1
                            node(sn)%col(problem%n+k,1)=col
                            ! Increment counter
                            col=col+1
                          end do
                          ! Region 2
                          ! Equation for region 2
                          node(sn)%row(1,2)=row
                          ! Increment counter
                          row=row+1
                          ! The variable p of region 2 is active
                          node(sn)%col(1,2)=col
                          ! Increment counter
                          col=col+1
                          ! Perfect bonding: u_k^(1)=u_k^(FE), Un^(2)=-u_k^(FE)n_k^(1)

                        ! ------------------------------------------------- !
                        ! BE (viscoelastic solid) - BE (viscoelastic solid) !
                        ! ------------------------------------------------- !

                        case (fbem_viscoelastic)
                          !
                          ! Index of equations for each coordinate k:
                          ! node(sn)%row(k,1): BIE for region 1
                          ! node(sn)%row(k,2): BIE for region 2
                          !
                          ! Index of variables for each coordinate k:
                          ! node(sn)%col(k,1): u_k for region 1
                          ! node(sn)%col(k,1): t_k for region 1 (active)
                          ! node(sn)%col(k,2): u_k for region 2
                          ! node(sn)%col(k,2): t_k for region 2 (active)
                          !
                          allocate (node(sn)%row(problem%n,2))
                          allocate (node(sn)%col(2*problem%n,2))
                          allocate (node(sn)%value_c(2*problem%n,2))
                          allocate (node(sn)%dvda_c(2*problem%n,2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equation and column to the variables
                          do k=1,problem%n
                            ! BIE equation for region 1
                            node(sn)%row(k,1)=row
                            ! BIE equation for region 2
                            node(sn)%row(k,2)=row+1
                            ! Increment counter
                            row=row+2
                            ! t_k for region 1 is always unknown
                            node(sn)%col(problem%n+k,1)=col
                            ! t_k for region 2 is always unknown
                            node(sn)%col(problem%n+k,2)=col+1
                            ! Increment counter
                            col=col+2
                            ! Perfect bonding: u_k^(1)=u_k^(2)=u_k^(FE)
                          end do

                        ! ------------------------------------------------------- !
                        ! BE (1: viscoelastic solid) - BE (2: poroelastic medium) !
                        ! ------------------------------------------------------- !

                        case (fbem_poroelastic)
                          !
                          ! Index of equations for each coordinate k:
                          ! node(sn)%row(k,1): BIE for region 1 (viscoelastic)
                          ! node(sn)%row(0,2): BIE region 2 (fluid)
                          ! node(sn)%row(k,2): BIE region 2 (solid skeleton)
                          !
                          ! Index of variables:
                          ! node(sn)%col(            k,1): u_k for region 1
                          ! node(sn)%col(            k,1): t_k for region 1 (active)
                          ! node(sn)%col(            0,2): tau for region 2 (active)
                          ! node(sn)%col(            k,2): u_k for region 2
                          ! node(sn)%col(problem%n+1  ,2): Un  for region 2
                          ! node(sn)%col(problem%n+1+k,2): t_k for region 2 (active)
                          !
                          allocate (node(sn)%row    (0:problem%n      ,2))
                          allocate (node(sn)%col    (0:(1+2*problem%n),2))
                          allocate (node(sn)%value_c(0:(1+2*problem%n),2))
                          allocate (node(sn)%dvda_c (0:(1+2*problem%n),2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equation and column to the variables
                          ! Region 1
                          do k=1,problem%n
                            ! BIE equation for region 1
                            node(sn)%row(k,1)=row
                            ! Increment counter
                            row=row+1
                            ! t_k for region 2 is always unknown
                            node(sn)%col(problem%n+k,1)=col
                            ! Increment counter
                            col=col+1
                          end do
                          ! Region 2
                          ! Fluid phase
                          ! BIE equation for region 2
                          node(sn)%row(0,2)=row
                          ! Increment counter
                          row=row+1
                          ! tau^(2) unknown
                          node(sn)%col(0,2)=col
                          ! Increment counter
                          col=col+1
                          ! Solid skeleton
                          do k=1,problem%n
                            ! BIE equation for region 2
                            node(sn)%row(k,2)=row
                            ! Increment counter
                            row=row+1
                            ! t_k^(1) unknown
                            node(sn)%col(problem%n+1+k,2)=col
                            ! Increment counter
                            col=col+1
                          end do
                          ! Perfect bonding: Un^(2)=u_k^(FE)n_k^(2), u_k^(1)=u_k^(2)=u_k^(FE)

                      end select

                    case (fbem_poroelastic)
                      select case (region(boundary(sb)%region(2))%type)

                        ! --------------------------------------------------- !
                        ! BE (1: poroelastic medium) - BE (2: inviscid fluid) !
                        ! --------------------------------------------------- !

                        case (fbem_potential)
                          !
                          ! Index of equations for each coordinate k:
                          ! node(sn)%row(0,1): BIE region 1 (fluid)
                          ! node(sn)%row(k,1): BIE region 1 (solid skeleton)
                          ! node(sn)%row(1,2): BIE for region 2 (inviscid fluid)
                          !
                          ! Index of variables:
                          ! node(sn)%col(            0,1): tau for region 1 (active)
                          ! node(sn)%col(            k,1): u_k for region 1
                          ! node(sn)%col(problem%n+1  ,1): Un  for region 1
                          ! node(sn)%col(problem%n+1+k,1): t_k for region 1 (active)
                          ! node(sn)%col(            1,2): p for region 2 (active)
                          ! node(sn)%col(            2,2): Un for region 2
                          !
                          allocate (node(sn)%row    (0:problem%n      ,2))
                          allocate (node(sn)%col    (0:(1+2*problem%n),2))
                          allocate (node(sn)%value_c(0:(1+2*problem%n),2))
                          allocate (node(sn)%dvda_c (0:(1+2*problem%n),2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equation and column to the variables
                          ! Region 1
                          ! Fluid phase
                          ! BIE equation for region 1
                          node(sn)%row(0,1)=row
                          ! Increment counter
                          row=row+1
                          ! tau^(1) unknown
                          node(sn)%col(0,1)=col
                          ! Increment counter
                          col=col+1
                          ! Solid skeleton
                          do k=1,problem%n
                            ! BIE equation for region 1
                            node(sn)%row(k,1)=row
                            ! Increment counter
                            row=row+1
                            ! t_k^(1) unknown
                            node(sn)%col(problem%n+1+k,1)=col
                            ! Increment counter
                            col=col+1
                          end do
                          ! Region 2
                          ! Equation for region 2
                          node(sn)%row(1,2)=row
                          ! Increment counter
                          row=row+1
                          ! The variable p of region 2 is active
                          node(sn)%col(1,2)=col
                          ! Increment counter
                          col=col+1
                          ! Perfect bonding: Un^(1)=-Un^(2)=u_k^(FE)n_k^(1), u_k^(1)=u_k^(FE)

                        ! ------------------------------------------------------- !
                        ! BE (1: poroelastic medium) - BE (2: viscoelastic solid) !
                        ! ------------------------------------------------------- !

                        case (fbem_viscoelastic)
                          !
                          ! Index of equations for each coordinate k:
                          ! node(sn)%row(0,1): BIE region 1 (fluid)
                          ! node(sn)%row(k,1): BIE region 1 (solid skeleton)
                          ! node(sn)%row(k,2): BIE for region 2 (viscoelastic)
                          !
                          ! Index of variables:
                          ! node(sn)%col(            0,1): tau for region 1 (active)
                          ! node(sn)%col(            k,1): u_k for region 1
                          ! node(sn)%col(problem%n+1  ,1): Un  for region 1
                          ! node(sn)%col(problem%n+1+k,1): t_k for region 1 (active)
                          ! node(sn)%col(            k,2): u_k for region 2
                          ! node(sn)%col(            k,2): t_k for region 2 (active)
                          !
                          allocate (node(sn)%row    (0:problem%n      ,2))
                          allocate (node(sn)%col    (0:(1+2*problem%n),2))
                          allocate (node(sn)%value_c(0:(1+2*problem%n),2))
                          allocate (node(sn)%dvda_c (0:(1+2*problem%n),2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equation and column to the variables
                          ! Region 1
                          ! Fluid phase
                          ! BIE equation for region 1
                          node(sn)%row(0,1)=row
                          ! Increment counter
                          row=row+1
                          ! tau^(1) unknown
                          node(sn)%col(0,1)=col
                          ! Increment counter
                          col=col+1
                          ! Solid skeleton
                          do k=1,problem%n
                            ! BIE equation for region 1
                            node(sn)%row(k,1)=row
                            ! Increment counter
                            row=row+1
                            ! t_k^(1) unknown
                            node(sn)%col(problem%n+1+k,1)=col
                            ! Increment counter
                            col=col+1
                          end do
                          ! Region 2
                          do k=1,problem%n
                            ! BIE equation for region 2
                            node(sn)%row(k,2)=row
                            ! Increment counter
                            row=row+1
                            ! t_k for region 2 is always unknown
                            node(sn)%col(problem%n+k,2)=col
                            ! Increment counter
                            col=col+1
                          end do
                          ! Perfect bonding: Un^(1)=u_k^(FE)n_k^(1), u_k^(1)=u_k^(2)=u_k^(FE)

                        ! ------------------------------------------------------- !
                        ! BE (1: poroelastic medium) - BE (2: poroelastic medium) !
                        ! ------------------------------------------------------- !

                        case (fbem_poroelastic)
                          !
                          ! Index of equations for each coordinate k:
                          ! node(sn)%row(0,1): BIE region 1 (fluid)
                          ! node(sn)%row(k,1): BIE region 1 (solid skeleton)
                          ! node(sn)%row(0,2): BIE for region 2 (fluid)
                          ! node(sn)%row(k,2): BIE for region 2 (solid skeleton)
                          !
                          ! Index of variables:
                          ! node(sn)%col(            0,1): tau for region 1 (active)
                          ! node(sn)%col(            k,1): u_k for region 1
                          ! node(sn)%col(problem%n+1  ,1): Un  for region 1
                          ! node(sn)%col(problem%n+1+k,1): t_k for region 1 (active)
                          ! node(sn)%col(            0,2): tau for region 2 (active)
                          ! node(sn)%col(            k,2): u_k for region 2
                          ! node(sn)%col(problem%n+1  ,2): Un  for region 2
                          ! node(sn)%col(problem%n+1+k,2): t_k for region 2 (active)
                          !
                          allocate (node(sn)%row    (0:problem%n      ,2))
                          allocate (node(sn)%col    (0:(1+2*problem%n),2))
                          allocate (node(sn)%value_c(0:(1+2*problem%n),2))
                          allocate (node(sn)%dvda_c (0:(1+2*problem%n),2,problem%n_designvariables))
                          ! Initialize
                          node(sn)%row=0
                          node(sn)%col=0
                          ! Assign row to the equation and column to the variables
                          ! Fluid phase
                          ! BIE equation for region 1
                          node(sn)%row(0,1)=row
                          ! BIE equation for region 2
                          node(sn)%row(0,2)=row+1
                          ! Increment counter
                          row=row+2
                          ! tau^(1) unknown
                          node(sn)%col(0,1)=col
                          ! tau^(2) unknown
                          node(sn)%col(0,2)=col+1
                          ! Increment counter
                          col=col+2
                          ! Perfect bonding: Un^(1)=-Un^(2)=u_k^(FE)n_k^(1)
                          ! Solid skeleton
                          do k=1,problem%n
                            ! BIE equation for region 1
                            node(sn)%row(k,1)=row
                            ! BIE equation for region 2
                            node(sn)%row(k,2)=row+1
                            ! Increment counter
                            row=row+2
                            ! t_k^(1) unknown
                            node(sn)%col(problem%n+1+k,1)=col
                            ! t_k^(2) unknown
                            node(sn)%col(problem%n+1+k,2)=col+1
                            ! Increment counter
                            col=col+2
                            ! Perfect bonding: u_k^(1)=u_k^(2)=u_k^(FE)
                          end do

                      end select

                  end select
              end select
            end if
          end do
        end do
      end do

      ! =============================
      ! COUPLED BE BODY LOAD ELEMENTS
      ! =============================

      do kb=1,region(kr)%n_be_bodyloads
        sb=region(kr)%be_bodyload(kb)
        sp=be_bodyload(sb)%part
        select case (be_bodyload(sb)%coupling)

          ! ----------------------------------
          ! FE BEAM TIP - BE LINE/SURFACE LOAD
          ! ----------------------------------

          case (fbem_bl_coupling_beam_tip)
            stop 'not yet'

          ! ---------------------------------------------------------
          ! (FE BEAM - BE LINE LOAD) AND (FE SHELL - BE SURFACE LOAD)
          ! ---------------------------------------------------------

          case (fbem_bl_coupling_beam_line,fbem_bl_coupling_shell_surface)
                do ke=1,part(sp)%n_elements
                  se=part(sp)%element(ke)
                  do kn=1,element(se)%n_nodes
                    sn=element(se)%node(kn)
                    if (node_used(sn).eqv.(.false.)) then
                      node_used(sn)=.true.
                      !
                      ! Index of equations for each coordinate k:
                      ! node(sn)%row(k,1): SBIE
                      !
                      ! Index of variables for each coordinate k:
                      ! node(sn)%col(          k,1): u_k (not active since u_k == u_k^(FE))
                      ! node(sn)%col(problem%n+k,1): b_k
                      !
                      allocate (node(sn)%row(problem%n,1))
                      allocate (node(sn)%col(2*problem%n,1))
                      allocate (node(sn)%value_c(2*problem%n,1))
                      allocate (node(sn)%dvda_c(2*problem%n,1,problem%n_designvariables))
                      ! Initialize
                      node(sn)%row=0
                      node(sn)%col=0
                      ! Assign row to the equation and column to the variables
                      do k=1,problem%n
                        ! Equation
                        node(sn)%row(k,1)=row
                        row=row+1
                        ! b_k is unknown
                        node(sn)%col(problem%n+k,1)=col
                        col=col+1
                      end do

                    end if
                  end do
                end do

          ! -----------------------
          ! FE SHELL - BE EDGE LOAD
          ! -----------------------

          case (fbem_bl_coupling_shell_edge)
            stop 'not yet'

        end select

      end do


    end if
  end do

  ! ================================================================================================================================
  ! FE REGIONS
  ! ================================================================================================================================

  !
  ! Nota: igual que en el caso estatico, pero cambiar value_r -> value_c, y dvda_r -> dvda_c
  !

  !
  ! ELEMENT-WISE
  !
  do kr=1,n_regions
    if (region(kr)%class.eq.fbem_fe) then

      ! ----------------------------------------------------------------------------------------------------------------------------
      ! RIGID ELEMENTS
      ! ----------------------------------------------------------------------------------------------------------------------------

      if (region(kr)%type.eq.fbem_rigid) then

        ! Nothing to do (the displacement field is reduced to the master node)

      ! ----------------------------------------------------------------------------------------------------------------------------
      ! FLEXIBLE ELEMENTS
      ! ----------------------------------------------------------------------------------------------------------------------------

      else
        do ks=1,region(kr)%n_fe_subregions
          ss=region(kr)%fe_subregion(ks)
          sp=fe_subregion(ss)%part
          do ke=1,part(sp)%n_elements
            se=part(sp)%element(ke)

            select case (element(se)%n_dimension)

              ! ====================================================================================================================
              ! ZERO-DIMENSIONAL ELEMENTS
              ! ====================================================================================================================

              case (0)

                !---------------------------------------------------------------------------------------------------------------
                ! POINT MASS ELEMENT
                !
                !
                ! Values to store
                !
                ! done below
                !
                ! DOF handling
                !
                allocate(element(se)%node_n_dof(element(se)%n_nodes))
                element(se)%node_n_dof=0 ! The number of DOF is taken from the maximum from the other elements connected to the node

              ! ====================================================================================================================
              ! ONE-DIMENSIONAL ELEMENTS
              ! ====================================================================================================================

              case (1)

                select case (element(se)%fe_type)

                  case (0)

                    !---------------------------------------------------------------------------------------------------------------
                    ! DEGENERATED BEAM FINITE ELEMENT
                    !
                    !
                    ! Values to store
                    !
                    ! Stress resultants (local axes v_k): value_c(k,n,1): NX,VY,BZ or NX,VY,VZ,BX,BY,BZ
                    ! Equilibrating loads (global axes) : value_c(k,n,2): FX,FY,MZ or FX,FY,FZ,MX,MY,MZ
                    !
                    allocate (element(se)%value_c(3*(problem%n-1),element(se)%n_nodes,2))
                    element(se)%value_c=0
                    !
                    ! DOF handling
                    !
                    allocate(element(se)%node_n_dof(element(se)%n_nodes))
                    element(se)%node_n_dof=3*(problem%n-1)
                    !
                    !---------------------------------------------------------------------------------------------------------------

                  case (1,2)

                    !---------------------------------------------------------------------------------------------------------------
                    ! STRAIGHT EULER-BERNOULLI AND TIMOSHENKO BEAM FINITE ELEMENTS
                    !
                    !
                    ! Values to store
                    !
                    ! Stress resultants (local axes ep_k): value_c(k,n,1): NX,VY,BZ OR NX,VY,VZ,BX,BY,BZ
                    ! Equilibrating loads (global axes)  : value_c(k,n,2): FX,FY[,MZ] OR FX,FY,FZ[,MX,MY,MZ]
                    !
                    allocate (element(se)%value_c(3*(problem%n-1),element(se)%n_nodes,2))
                    element(se)%value_c=0
                    !
                    ! DOF handling
                    !
                    allocate(element(se)%node_n_dof(element(se)%n_nodes))
                    select case (element(se)%type)
                      case (fbem_line2)
                        element(se)%node_n_dof(1)=3*(problem%n-1)
                        element(se)%node_n_dof(2)=3*(problem%n-1)
                      case (fbem_line3)
                        if (element(se)%fe_options(1).eq.0) then
                          element(se)%node_n_dof(1)=3*(problem%n-1)
                          element(se)%node_n_dof(2)=3*(problem%n-1)
                          element(se)%node_n_dof(3)=problem%n
                        else
                          element(se)%node_n_dof(1)=3*(problem%n-1)
                          element(se)%node_n_dof(2)=3*(problem%n-1)
                          element(se)%node_n_dof(3)=3*(problem%n-1)
                        end if
                      case default
                        call fbem_error_message(error_unit,0,'element',element(se)%id,'strbeam element only available for line2 or line3 mesh element')
                    end select
                    !
                    !---------------------------------------------------------------------------------------------------------------

                  case (3)

                    !---------------------------------------------------------------------------------------------------------------
                    ! BAR FINITE ELEMENTS
                    !
                    !
                    ! Values to store
                    !
                    ! Stress resultants                  : value_c(k,n,1): NX
                    ! Equilibrating loads (global axes)  : value_c(k,n,2): FX,FY,FZ
                    !
                    allocate (element(se)%value_c(problem%n,element(se)%n_nodes,2))
                    element(se)%value_c=0
                    !
                    ! DOF handling
                    !
                    allocate(element(se)%node_n_dof(element(se)%n_nodes))
                    select case (element(se)%type)
                      case (fbem_line2)
                        element(se)%node_n_dof=problem%n
                      case default
                        call fbem_error_message(error_unit,0,'element',element(se)%id,'bar element only available for line2 mesh element')
                    end select
                    !
                    !---------------------------------------------------------------------------------------------------------------

                  case (4)

                    !---------------------------------------------------------------------------------------------------------------
                    ! DISCRETE TRANSLATIONAL SPRING FINITE ELEMENTS
                    !
                    !
                    ! Values to store
                    !
                    ! Spring forces                     : value_c(k,n,1): NX,NY,NZ
                    ! Equilibrating loads (nodal axes)  : value_c(k,n,2): FX,FY,FZ
                    !
                    allocate (element(se)%value_c(problem%n,2,2))
                    element(se)%value_c=0
                    !
                    ! DOF handling
                    !
                    allocate(element(se)%node_n_dof(2))
                    select case (element(se)%type)
                      case (fbem_line2)
                        element(se)%node_n_dof=problem%n
                      case default
                        call fbem_error_message(error_unit,0,'element',element(se)%id,'distra element only available for line2 mesh element')
                    end select
!                    allocate(element(se)%node_n_dof(element(se)%n_nodes))
!                    element(se)%node_n_dof=0
!                    element(se)%node_n_dof(1:2)=problem%n
                    !
                    !---------------------------------------------------------------------------------------------------------------

                  case (5)

                    !---------------------------------------------------------------------------------------------------------------
                    ! DISCRETE ROTATIONAL/TRANSLATIONAL SPRING FINITE ELEMENTS
                    !
                    !
                    ! Values to store
                    !
                    ! Spring forces/moments             : value_c(k,n,1): NX,NY,NZ,BX,BY,BZ
                    ! Equilibrating loads (nodal axes)  : value_c(k,n,2): FX,FY,FZ,MX,MY,MZ
                    !
                    allocate (element(se)%value_c(3*(problem%n-1),2,2))
                    element(se)%value_c=0
                    !
                    ! DOF handling
                    !
                    allocate(element(se)%node_n_dof(2))
                    select case (element(se)%type)
                      case (fbem_line2)
                        element(se)%node_n_dof=3*(problem%n-1)
                      case default
                        call fbem_error_message(error_unit,0,'element',element(se)%id,'disrotra element only available for line2 mesh element')
                    end select
!                    allocate(element(se)%node_n_dof(element(se)%n_nodes))
!                    element(se)%node_n_dof=0
!                    element(se)%node_n_dof(1:2)=3*(problem%n-1)
                    !
                    !---------------------------------------------------------------------------------------------------------------

                  case default

                    !---------------------------------------------------------------------------------------------------------------
                    ! OTHER TYPES
                    !
                    call fbem_error_message(error_unit,0,'element',element(se)%id,'invalid type of 1D element')
                    !
                    !---------------------------------------------------------------------------------------------------------------

                end select

              ! ====================================================================================================================
              ! TWO-DIMENSIONAL ELEMENTS
              ! ====================================================================================================================

              case (2)

                select case (problem%n)

                  case (2)

                    !---------------------------------------------------------------------------------------------------------------
                    ! SOLID / CONTINUUM ELEMENTS
                    !
                    !
                    ! Values to store
                    !
                    ! Stress tensor (global axes)      : value_c(k,n,1): SXX,SYY,SXY
                    ! Equilibrating loads (global axes): value_c(k,n,2): FX,FY
                    !
                    allocate (element(se)%value_c(3,element(se)%n_nodes,2))
                    element(se)%value_c=0
                    !
                    ! DOF handling
                    !
                    allocate(element(se)%node_n_dof(element(se)%n_nodes))
                    element(se)%node_n_dof=2
                    !
                    !---------------------------------------------------------------------------------------------------------------

                  case (3)

                    !---------------------------------------------------------------------------------------------------------------
                    ! DEGENERATED SHELL FINITE ELEMENT
                    !
                    !
                    ! Values to store
                    !
                    ! Stress resultants (local axes ep_k)               : value_c(k,n,1): NX,NY,NXY,MX,MY,MXY,VX,VY
                    ! Equilibrating loads (global u, global/local theta): value_c(k,n,2): FX,FY,FZ[,MA,MB][,MX,MY,MZ]
                    !
                    allocate (element(se)%value_c(8,element(se)%n_nodes,2))
                    element(se)%value_c=0
                    !
                    ! DOF handling
                    !
                    allocate(element(se)%node_n_dof(element(se)%n_nodes))
                    do kn=1,element(se)%n_nodes
                      sn=element(se)%node(kn)
                      if (node(sn)%n_dof.eq.0) then
                        element(se)%node_n_dof(kn)=5
                        if (node(sn)%is_singular.or.(node(sn)%rigid_link.ne.0).or.(node(sn)%n_symplanes.gt.0)) then
                          element(se)%node_n_dof(kn)=6
                        end if
                      else
                        element(se)%node_n_dof(kn)=node(sn)%n_dof
                        select case (node(sn)%n_dof)
                          case (5)
                            if (node(sn)%is_singular.or.(node(sn)%rigid_link.ne.0).or.(node(sn)%n_symplanes.gt.0)) then
                              call fbem_error_message(error_unit,0,'node',node(sn)%id,'this node can not be a 5 DOF node.')
                            end if
                          case (6)
                            ! nothing to do
                          case default
                            write(error_unit,*) node(sn)%id, node(sn)%n_dof
                            call fbem_error_message(error_unit,0,'node',node(sn)%id,'invalid type of shell node.')
                        end select
                      end if
                    end do
                    !
                    !---------------------------------------------------------------------------------------------------------------

                end select

              ! ====================================================================================================================
              ! THREE-DIMENSIONAL ELEMENTS
              ! ====================================================================================================================

              case (3)

                !-------------------------------------------------------------------------------------------------------------------
                ! SOLID / CONTINUUM ELEMENTS
                !
                !
                ! Values to store
                !
                ! Stress tensor (global axes)      : value_c(k,n,1): SXX,SYY,SZZ,SXY,SXZ,SYZ
                ! Equilibrating loads (global axes): value_c(k,n,2): FX,FY,FZ


                !
                ! NOTA: es mas logico poner las cargas equilibrantes en el indice 1
                !


                !
                allocate (element(se)%value_c(6,element(se)%n_nodes,2))
                element(se)%value_c=0
                !
                ! DOF handling
                !
                allocate(element(se)%node_n_dof(element(se)%n_nodes))
                element(se)%node_n_dof=3
                !
                call fbem_error_message(error_unit,0,'element',element(se)%id,'3D elements not available yet')
                !
                !-------------------------------------------------------------------------------------------------------------------

            end select

            element(se)%n_dof=sum(element(se)%node_n_dof)

          end do
        end do
      end if
    end if
  end do
  !
  ! NODE-WISE
  !
  !
  ! RIGID REGIONS
  !
  do kr=1,n_regions
    if (region(kr)%class.eq.fbem_fe) then
      if (region(kr)%type.eq.fbem_rigid) then
        !
        ! MASTER NODE
        !
        sn=region(kr)%master_node
        node_used(sn)=.true.
        !
        ! DOF handling
        !
        ! Index of equations:
        !   2D:
        ! node(sn)%row(1,1): equilibrium f_1
        ! node(sn)%row(2,1): equilibrium f_2
        ! node(sn)%row(3,1): equilibrium m_3
        !   3D:
        ! node(sn)%row(1,1): equilibrium f_1
        ! node(sn)%row(2,1): equilibrium f_2
        ! node(sn)%row(3,1): equilibrium f_3
        ! node(sn)%row(4,1): equilibrium m_1
        ! node(sn)%row(5,1): equilibrium m_2
        ! node(sn)%row(6,1): equilibrium m_3
        !
        ! Index of variables:
        !   2D:
        ! node(sn)%col(1,1): u_1
        ! node(sn)%col(2,1): u_2
        ! node(sn)%col(3,1): theta_3
        !   3D:
        ! node(sn)%col(1,1): u_1
        ! node(sn)%col(2,1): u_2
        ! node(sn)%col(3,1): u_3
        ! node(sn)%col(4,1): theta_1
        ! node(sn)%col(5,1): theta_2
        ! node(sn)%col(6,1): theta_3
        !
        node(sn)%n_dof=3*(problem%n-1)
        !
        ! Copy the n_dof to the point mass element if present
        !
        do ke=1,node(sn)%n_elements
          se=node(sn)%element(ke)
          if (element(se)%n_dimension.eq.0) then
            element(se)%node_n_dof=node(sn)%n_dof
            element(se)%n_dof=node(sn)%n_dof
            allocate (element(se)%value_c(element(se)%n_dof,element(se)%n_nodes,2))
            element(se)%value_c=0
          end if
        end do
        !
        !
        !
        allocate (node(sn)%row(node(sn)%n_dof,1))
        allocate (node(sn)%col(node(sn)%n_dof,1))
        node(sn)%row=0
        node(sn)%col=0
        ! Assign row and column to each DOF depending on the BC
        do k=1,node(sn)%n_dof
          select case (node(sn)%ctype(k,1))
            ! DOF is known
            case (0)
            ! DOF is unknown
            case (1)
              ! Equilibrium equation
              node(sn)%row(k,1)=row
              row=row+1
              ! DOF global number
              node(sn)%col(k,1)=col
              col=col+1
          end select
        end do
        !
        ! Values to store
        !
        ! Nodal gen. displacements: value_c(k,1): UX,UY,RZ OR UX,UY,UZ,RX,RY,RZ
        ! Nodal reactions/loads   : value_c(k,2): FX,FY,MZ OR FX,FY,FZ,MX,MY,MZ
        !
        allocate (node(sn)%value_c(node(sn)%n_dof,2))
        node(sn)%value_c=0
        allocate (node(sn)%dvda_c (node(sn)%n_dof,2,problem%n_designvariables))
        node(sn)%dvda_c=0
        !
        ! SLAVE NODES belonging to flexible regions are treated below
        !
      end if
    end if
  end do
  !
  ! FLEXIBLE REGIONS
  !
  do kr=1,n_regions
    if (region(kr)%class.eq.fbem_fe) then
      if (region(kr)%type.ne.fbem_rigid) then
        do ks=1,region(kr)%n_fe_subregions
          ss=region(kr)%fe_subregion(ks)
          sp=fe_subregion(ss)%part
          do kn=1,part(sp)%n_nodes
            sn=part(sp)%node(kn)
            if (.not.node_used(sn)) then
              ! Flag the node as already studied
              node_used(sn)=.true.
              ! If it is a slave node of a rigid region, it does not have its own active DOF
              if (node(sn)%master.ne.0) cycle
              ! Find the maximum number of DOFs between all elements that share the node (except 0D point mass elements)
              ndof_fe_node=0
              do ke=1,node(sn)%n_elements
                se=node(sn)%element(ke)
                kn2=node(sn)%element_node_iid(ke)
                if ((ndof_fe_node.lt.element(se)%node_n_dof(kn2)).and.(element(se)%n_dimension.gt.0)) then
                  ndof_fe_node=element(se)%node_n_dof(kn2)
                end if
              end do
              !
              ! Check compatibility between element nodes
              !
              select case (problem%n)
                case (2)
                  select case (ndof_fe_node)
                    case (0,2,3)
                      node(sn)%n_dof=ndof_fe_node
                    case default
                      call fbem_error_message(error_unit,0,'node',node(sn)%id,'invalid ndof_fe_node value.')
                  end select
                case (3)
                  select case (ndof_fe_node)
                    case (0,3,5)
                      node(sn)%n_dof=ndof_fe_node
                    case (6)
                      do ke=1,node(sn)%n_elements
                        se=node(sn)%element(ke)
                        kn2=node(sn)%element_node_iid(ke)
                        if (element(se)%node_n_dof(kn2).eq.2) then
                          call fbem_error_message(error_unit,0,'node',node(sn)%id,'invalid ndof_fe_node value.')
                        end if
                        if (element(se)%node_n_dof(kn2).eq.5) then
                          call fbem_warning_message(error_unit,0,'node',node(sn)%id,&
                          'this node is connected to element nodes of 6 DOF, it is automatically changed to 6 DOF node')
                          element(se)%node_n_dof(kn2)=6
                          element(se)%n_dof=sum(element(se)%node_n_dof)
                        end if
                      end do
                      node(sn)%n_dof=ndof_fe_node
                    case default
                      call fbem_error_message(error_unit,0,'node',node(sn)%id,'invalid ndof_fe_node value.')
                  end select
                case default
                  call fbem_error_message(error_unit,0,'problem%n',problem%n,'invalid value.')
              end select
              !
              ! Copy the n_dof to the point mass element if present
              !
              do ke=1,node(sn)%n_elements
                se=node(sn)%element(ke)
                if (element(se)%n_dimension.eq.0) then
                  element(se)%node_n_dof=node(sn)%n_dof
                  element(se)%n_dof=node(sn)%n_dof
                  allocate (element(se)%value_c(element(se)%n_dof,element(se)%n_nodes,2))
                  element(se)%value_c=0
                end if
              end do
              !
              ! Build node mapping to the LSE
              !
              allocate (node(sn)%row(node(sn)%n_dof,1))
              allocate (node(sn)%col(node(sn)%n_dof,1))
              node(sn)%row=0
              node(sn)%col=0
              ! Assign row and column to each degree of freedom depending on the B.C.
              do k=1,node(sn)%n_dof
                select case (node(sn)%ctype(k,1))
                  ! The displacement is already known
                  case (0)
                  ! The displacement is unknown
                  case (1)
                    ! FEM equation
                    node(sn)%row(k,1)=row
                    ! Increment counter
                    row=row+1
                    ! Degree of freedom k
                    node(sn)%col(k,1)=col
                    ! Increment counter
                    col=col+1
                end select
              end do
              !
              ! Values to store
              !
              !
              ! Nodal gen. displacements: value_c(k,1): UX,UY[,RZ] OR UX,UY,UZ[,RA,RB][,RX,RY,RZ]
              ! Nodal reactions/loads   : value_c(k,2): FX,FY[,MZ] OR FX,FY,FZ[,MA,MB][,MX,MY,MZ]
              !
              allocate (node(sn)%value_c(node(sn)%n_dof,2))
              node(sn)%value_c=0
              allocate (node(sn)%dvda_c (node(sn)%n_dof,2,problem%n_designvariables))
              node(sn)%dvda_c=0
            end if
          end do
        end do
      end if
    end if
  end do
  !
  ! RIGID REGIONS
  !
  do kr=1,n_regions
    if (region(kr)%class.eq.fbem_fe) then
      if (region(kr)%type.eq.fbem_rigid) then
        !
        ! SLAVE NODES not belonging to flexible regions
        !
        snm=region(kr)%master_node
        do kn=1,node(snm)%n_slaves
          sn=node(snm)%slave(kn)
          if (.not.node_used(sn)) then
            node_used(sn)=.true.
            select case (problem%n)
              case (2)
                !
                ! Index of equations:
                ! node(sn)%row(1,1): equilibrium f_1
                ! node(sn)%row(2,1): equilibrium f_2
                ! node(sn)%row(3,1): equilibrium m
                !
                ! Index of variables:
                ! node(sn)%col(1,1): u_1
                ! node(sn)%col(2,1): u_2
                ! node(sn)%col(3,1): theta
                !
                node(sn)%n_dof=3
              case (3)
                !
                ! Index of equations:
                ! node(sn)%row(1,1): equilibrium f_1
                ! node(sn)%row(2,1): equilibrium f_2
                ! node(sn)%row(3,1): equilibrium f_3
                ! node(sn)%row(4,1): equilibrium m_1
                ! node(sn)%row(5,1): equilibrium m_2
                ! node(sn)%row(6,1): equilibrium m_3
                !
                ! Index of variables:
                ! node(sn)%col(1,1): u_1
                ! node(sn)%col(2,1): u_2
                ! node(sn)%col(3,1): u_3
                ! node(sn)%col(4,1): theta_1
                ! node(sn)%col(5,1): theta_2
                ! node(sn)%col(6,1): theta_3
                !
                node(sn)%n_dof=6
              case default
                call fbem_error_message(error_unit,0,'',0,'invalid problem ambient space dimension')
            end select
            !
            ! Copy the n_dof to the point mass element if present
            !
            do ke=1,node(sn)%n_elements
              se=node(sn)%element(ke)
              if (element(se)%n_dimension.eq.0) then
                element(se)%node_n_dof=node(sn)%n_dof
                element(se)%n_dof=node(sn)%n_dof
                allocate (element(se)%value_c(element(se)%n_dof,element(se)%n_nodes,2))
                element(se)%value_c=0
              end if
            end do
            !
            !
            !
            allocate (node(sn)%row    (node(sn)%n_dof,1))
            allocate (node(sn)%col    (node(sn)%n_dof,1))
            ! Initialize
            node(sn)%row=0
            node(sn)%col=0
            allocate (node(sn)%value_c(node(sn)%n_dof,2))
            node(sn)%value_c=0
            allocate (node(sn)%dvda_c (node(sn)%n_dof,2,problem%n_designvariables))
            node(sn)%dvda_c=0
          end if
        end do
      end if
    end if
  end do

!~   !
!~   ! OLD ............................................
!~   !

!~   ! ================================================================================================================================
!~   ! FE REGIONS
!~   ! ================================================================================================================================

!~   !
!~   ! RIGID REGIONS
!~   !
!~   do kr=1,n_regions
!~     if (region(kr)%class.eq.fbem_fe) then
!~       if (region(kr)%type.eq.fbem_rigid) then
!~         sn=region(kr)%master_node
!~         node_used(sn)=.true.
!~         select case (problem%n)
!~           case (2)
!~             !
!~             ! Index of equations:
!~             ! node(sn)%row(1,1): equilibrium f_1
!~             ! node(sn)%row(2,1): equilibrium f_2
!~             ! node(sn)%row(3,1): equilibrium m
!~             !
!~             ! Index of variables:
!~             ! node(sn)%col(1,1): u_1
!~             ! node(sn)%col(2,1): u_2
!~             ! node(sn)%col(3,1): theta
!~             !
!~             node(sn)%n_dof=3
!~           case (3)
!~             !
!~             ! Index of equations:
!~             ! node(sn)%row(1,1): equilibrium f_1
!~             ! node(sn)%row(2,1): equilibrium f_2
!~             ! node(sn)%row(3,1): equilibrium f_3
!~             ! node(sn)%row(4,1): equilibrium m_1
!~             ! node(sn)%row(5,1): equilibrium m_2
!~             ! node(sn)%row(6,1): equilibrium m_3
!~             !
!~             ! Index of variables:
!~             ! node(sn)%col(1,1): u_1
!~             ! node(sn)%col(2,1): u_2
!~             ! node(sn)%col(3,1): u_3
!~             ! node(sn)%col(4,1): theta_1
!~             ! node(sn)%col(5,1): theta_2
!~             ! node(sn)%col(6,1): theta_3
!~             !
!~             node(sn)%n_dof=6
!~           case default
!~             call fbem_error_message(error_unit,0,'',0,'invalid problem ambient space dimension')
!~         end select
!~         allocate (node(sn)%row    (node(sn)%n_dof,1))
!~         allocate (node(sn)%col    (node(sn)%n_dof,1))
!~         allocate (node(sn)%value_c(node(sn)%n_dof,1))
!~         allocate (node(sn)%dvda_c (node(sn)%n_dof,1,problem%n_designvariables))
!~         ! Initialize
!~         node(sn)%row=0
!~         node(sn)%col=0
!~         ! Assign row and column to each DOF depending on the BC
!~         do k=1,node(sn)%n_dof
!~           select case (node(sn)%ctype(k,1))
!~             ! DOF is known
!~             case (0)
!~             ! DOF is unknown
!~             case (1)
!~               ! Equilibrium equation
!~               node(sn)%row(k,1)=row
!~               row=row+1
!~               ! DOF global number
!~               node(sn)%col(k,1)=col
!~               col=col+1
!~           end select
!~         end do
!~       end if
!~     end if
!~   end do
!~   !
!~   ! FLEXIBLE REGIONS
!~   !
!~   do kr=1,n_regions
!~     if (region(kr)%class.eq.fbem_fe) then

!~       if (region(kr)%type.ne.fbem_rigid) then
!~         do ks=1,region(kr)%n_fe_subregions
!~           ss=region(kr)%fe_subregion(ks)
!~           sp=fe_subregion(ss)%part
!~           do ke=1,part(sp)%n_elements
!~             se=part(sp)%element(ke)
!~             select case (problem%n)

!~               ! ====
!~               !  2D
!~               ! ====

!~               case (2)
!~                 select case (element(se)%n_dimension)

!~                   ! ARC/BEAM/ROD
!~                   case (1)
!~                     select case (element(se)%fe_type)

!~                       ! -------------------------------
!~                       ! DEGENERATED BEAM FINITE ELEMENT
!~                       ! -------------------------------

!~                       ! Base type can be any lineX, where all nodes has displacements and rotations.
!~                       case (0)
!~                         !
!~                         ! Index of equations:
!~                         ! node(sn)%row(1,1): u_1
!~                         ! node(sn)%row(2,1): u_2
!~                         ! node(sn)%row(3,1): theta
!~                         !
!~                         ! Index of variables:
!~                         ! node(sn)%col(1,1): u_1
!~                         ! node(sn)%col(2,1): u_2
!~                         ! node(sn)%col(3,1): theta
!~                         !
!~                         do kn=1,element(se)%n_nodes
!~                           sn=element(se)%node(kn)
!~                           if (.not.node_used(sn)) then
!~                             node_used(sn)=.true.
!~                             node(sn)%n_dof=3
!~                             allocate (node(sn)%row    (node(sn)%n_dof,1))
!~                             allocate (node(sn)%col    (node(sn)%n_dof,1))
!~                             allocate (node(sn)%value_c(node(sn)%n_dof,1))
!~                             allocate (node(sn)%dvda_c (node(sn)%n_dof,1,problem%n_designvariables))
!~                             ! Initialize
!~                             node(sn)%row=0
!~                             node(sn)%col=0
!~                             ! Assign row and column to each degree of freedom depending on the B.C.
!~                             do k=1,node(sn)%n_dof
!~                               select case (node(sn)%ctype(k,1))
!~                                 ! The displacement is already known
!~                                 case (0)
!~                                 ! The displacement is unknown
!~                                 case (1)
!~                                   ! FEM equation
!~                                   node(sn)%row(k,1)=row
!~                                   ! Increment counter
!~                                   row=row+1
!~                                   ! Degree of freedom k
!~                                   node(sn)%col(k,1)=col
!~                                   ! Increment counter
!~                                   col=col+1
!~                               end select
!~                             end do
!~                           end if
!~                         end do

!~                       ! ----------------------------
!~                       ! EULER-BERNOULLI BEAM ELEMENT
!~                       ! ----------------------------

!~                       ! Base type can only be line3, where vertices have displacements and rotations and the mid-node only displacements.
!~                       case (1,2)
!~                         do kn=1,element(se)%n_nodes
!~                           sn=element(se)%node(kn)
!~                           if (.not.node_used(sn)) then
!~                             node_used(sn)=.true.
!~                             ! Number of DOF of the node
!~                             select case (kn)
!~                               ! Without rotations
!~                               case (3)
!~                                 ndof_fe_node=problem%n
!~                               ! With rotations
!~                               case (1,2)
!~                                 ndof_fe_node=3*(problem%n-1)
!~                             end select
!~                             node(sn)%n_dof=ndof_fe_node
!~                             !
!~                             ! Se deben chequear los otros elementos que contienen al nodo y ver quien requiere mas grados de libertad....
!~                             ! o mejor aun, darle la vuelta a la tortilla y hacer el bucle en nodos "fe", y hacer esto mismo
!~                             ! ...para versiones futuras ....
!~                             !
!~                             ! Index of equations:
!~                             ! node(sn)%row(             1:problem%n,1): FEM equation for u_k
!~                             ! node(sn)%row(problem%n+1:ndof_fe_node,1): FEM equation for theta_k
!~                             !
!~                             ! Index of variables:
!~                             ! node(sn)%col(             1:problem%n,1): u_k
!~                             ! node(sn)%col(problem%n+1:ndof_fe_node,1): theta_k
!~                             allocate (node(sn)%row    (node(sn)%n_dof,1))
!~                             allocate (node(sn)%col    (node(sn)%n_dof,1))
!~                             allocate (node(sn)%value_c(node(sn)%n_dof,1))
!~                             allocate (node(sn)%dvda_c (node(sn)%n_dof,1,problem%n_designvariables))
!~                             ! Initialize
!~                             do k=1,ndof_fe_node
!~                               node(sn)%row(k,1)=0
!~                               node(sn)%col(k,1)=0
!~                             end do
!~                             ! Assign row and column to each degree of freedom depending on the B.C.
!~                             do k=1,ndof_fe_node
!~                               select case (node(sn)%ctype(k,1))
!~                                 ! The displacement is already known
!~                                 case (0)
!~                                 ! The displacement is unknown
!~                                 case (1)
!~                                   ! FEM equation
!~                                   node(sn)%row(k,1)=row
!~                                   ! Increment counter
!~                                   row=row+1
!~                                   ! Degree of freedom k
!~                                   node(sn)%col(k,1)=col
!~                                   ! Increment counter
!~                                   col=col+1
!~                               end select
!~                             end do
!~                           end if
!~                         end do
!~                       case default
!~                         stop 'Not valid 1D FEM element'
!~                     end select

!~                   ! SOLID
!~                   case (2)
!~                     !
!~                     ! Index of equations:
!~                     ! node(sn)%row(1:problem%n,1): FEM equation for u_k
!~                     !
!~                     ! Index of variables:
!~                     ! node(sn)%col(1:problem%n,1): u_k
!~                     !
!~                     do kn=1,element(se)%n_nodes
!~                       sn=element(se)%node(kn)
!~                       if (.not.node_used(sn)) then
!~                         node_used(sn)=.true.
!~                         node(sn)%n_dof=problem%n
!~                         allocate (node(sn)%row    (node(sn)%n_dof,1))
!~                         allocate (node(sn)%col    (node(sn)%n_dof,1))
!~                         allocate (node(sn)%value_c(node(sn)%n_dof,1))
!~                         allocate (node(sn)%dvda_c(node(sn)%n_dof,1,problem%n_designvariables))
!~                         ! Initialize
!~                         node(sn)%row=0
!~                         node(sn)%col=0
!~                         ! Assign row and column to each degree of freedom depending on the B.C.
!~                         do k=1,node(sn)%n_dof
!~                           select case (node(sn)%ctype(k,1))
!~                             ! The displacement is known
!~                             case (0)
!~                             ! The displacement is unknown
!~                             case (1)
!~                               ! FEM equation
!~                               node(sn)%row(k,1)=row
!~                               ! Increment counter
!~                               row=row+1
!~                               ! Degree of freedom k
!~                               node(sn)%col(k,1)=col
!~                               ! Increment counter
!~                               col=col+1
!~                           end select
!~                         end do
!~                       end if
!~                     end do

!~                 end select

!~               ! ====
!~               !  3D
!~               ! ====

!~               case (3)
!~                 select case (element(se)%n_dimension)

!~                   ! ARC/BEAM/ROD
!~                   case (1)
!~                     select case (element(se)%fe_type)
!~                       !
!~                       ! DEGENERATED BEAM FINITE ELEMENT
!~                       !
!~                       ! Base type can be any lineX, where all nodes has displacements and rotations.
!~                       case (0)
!~                         !
!~                         ! Index of equations:
!~                         ! node(sn)%row(1,1): u_1
!~                         ! node(sn)%row(2,1): u_2
!~                         ! node(sn)%row(3,1): u_3
!~                         ! node(sn)%row(4,1): theta_1
!~                         ! node(sn)%row(5,1): theta_2
!~                         ! node(sn)%row(6,1): theta_3
!~                         !
!~                         ! Index of variables:
!~                         ! node(sn)%col(1,1): u_1
!~                         ! node(sn)%col(2,1): u_2
!~                         ! node(sn)%col(3,1): u_3
!~                         ! node(sn)%col(4,1): theta_1
!~                         ! node(sn)%col(5,1): theta_2
!~                         ! node(sn)%col(6,1): theta_3
!~                         !
!~                         do kn=1,element(se)%n_nodes
!~                           sn=element(se)%node(kn)
!~                           if (.not.node_used(sn)) then
!~                             node_used(sn)=.true.
!~                             node(sn)%n_dof=6
!~                             allocate (node(sn)%row    (node(sn)%n_dof,1))
!~                             allocate (node(sn)%col    (node(sn)%n_dof,1))
!~                             allocate (node(sn)%value_c(node(sn)%n_dof,1))
!~                             allocate (node(sn)%dvda_c (node(sn)%n_dof,1,problem%n_designvariables))
!~                             ! Initialize
!~                             node(sn)%row=0
!~                             node(sn)%col=0
!~                             ! Assign row and column to each degree of freedom depending on the B.C.
!~                             do k=1,node(sn)%n_dof
!~                               select case (node(sn)%ctype(k,1))
!~                                 ! The displacement is already known
!~                                 case (0)
!~                                 ! The displacement is unknown
!~                                 case (1)
!~                                   ! FEM equation
!~                                   node(sn)%row(k,1)=row
!~                                   ! Increment counter
!~                                   row=row+1
!~                                   ! Degree of freedom k
!~                                   node(sn)%col(k,1)=col
!~                                   ! Increment counter
!~                                   col=col+1
!~                               end select
!~                             end do
!~                           end if
!~                         end do
!~                       !
!~                       ! EULER-BERNOULLI BEAM ELEMENT
!~                       !
!~                       ! Base type can only be line3, where vertices have displacements and rotations and the mid-node only displacements.
!~                       case (1,2)
!~                         do kn=1,element(se)%n_nodes
!~                           sn=element(se)%node(kn)
!~                           if (.not.node_used(sn)) then
!~                             node_used(sn)=.true.
!~                             ! Number of DOF of the node
!~                             select case (kn)
!~                               ! Without rotations
!~                               case (3)
!~                                 ndof_fe_node=problem%n
!~                               ! With rotations
!~                               case (1,2)
!~                                 ndof_fe_node=3*(problem%n-1)
!~                             end select
!~                             node(sn)%n_dof=ndof_fe_node
!~                             !
!~                             ! Se deben chequear los otros elementos que contienen al nodo y ver quien requiere mas grados de libertad....
!~                             ! o mejor aun, darle la vuelta a la tortilla y hacer el bucle en nodos "fe", y hacer esto mismo
!~                             ! ...para versiones futuras ....
!~                             !
!~                             ! Index of equations:
!~                             ! node(sn)%row(             1:problem%n,1): FEM equation for u_k
!~                             ! node(sn)%row(problem%n+1:ndof_fe_node,1): FEM equation for theta_k
!~                             !
!~                             ! Index of variables:
!~                             ! node(sn)%col(             1:problem%n,1): u_k
!~                             ! node(sn)%col(problem%n+1:ndof_fe_node,1): theta_k
!~                             allocate (node(sn)%row    (node(sn)%n_dof,1))
!~                             allocate (node(sn)%col    (node(sn)%n_dof,1))
!~                             allocate (node(sn)%value_c(node(sn)%n_dof,1))
!~                             allocate (node(sn)%dvda_c (node(sn)%n_dof,1,problem%n_designvariables))
!~                             ! Initialize
!~                             do k=1,ndof_fe_node
!~                               node(sn)%row(k,1)=0
!~                               node(sn)%col(k,1)=0
!~                             end do
!~                             ! Assign row and column to each degree of freedom depending on the B.C.
!~                             do k=1,ndof_fe_node
!~                               select case (node(sn)%ctype(k,1))
!~                                 ! The displacement is already known
!~                                 case (0)
!~                                 ! The displacement is unknown
!~                                 case (1)
!~                                   ! FEM equation
!~                                   node(sn)%row(k,1)=row
!~                                   ! Increment counter
!~                                   row=row+1
!~                                   ! Degree of freedom k
!~                                   node(sn)%col(k,1)=col
!~                                   ! Increment counter
!~                                   col=col+1
!~                               end select
!~                             end do
!~                           end if
!~                         end do
!~                      case default
!~                         stop 'Not valid 1D FEM element'
!~                     end select





!~                   ! SHELL/PLATE
!~                   case (2)
!~                   !
!~                   ! Shell obtained from degenerated 3D solid
!~                   !
!~                   ! 2 types of nodes depending on the connectivity of each shell
!~                   !   Type 5: 5 DOF, 3 global axes displacements and 2 local axes rotations
!~                   !   Type 6: 6 DOF, 3 global axes displacements and 3 global axes rotations
!~                   !
!~                   allocate (element(se)%value_c(8,element(se)%n_nodes,1))
!~                   do kn=1,element(se)%n_nodes
!~                     sn=element(se)%node(kn)
!~                     if (.not.node_used(sn)) then
!~                       node_used(sn)=.true.
!~                       ! The number of degrees of freedom of the node
!~                       ! If 0, apply default values
!~                       if (node(sn)%n_dof.eq.0) then
!~                         node(sn)%n_dof=5
!~                         if (node(sn)%is_singular) node(sn)%n_dof=6
!~                       else
!~                         select case (node(sn)%n_dof)
!~                           case (5)
!~                             if (node(sn)%is_singular) then
!~                               call fbem_error_message(error_unit,0,'node',node(sn)%id,'this node can not be a 5 DOF node.')
!~                             end if
!~                           case (6)
!~                             !write(*,*) 'node',node(sn)%id,'is a user defined 6 DOF node'
!~                           case default
!~                             call fbem_error_message(error_unit,0,'node',node(sn)%id,'invalid type of shell node.')
!~                         end select
!~                       end if
!~                       ndof_fe_node=node(sn)%n_dof
!~                       !
!~                       ! Index of equations:
!~                       ! node(sn)%row(1:problem%n,1): FEM equation for u_k
!~                       ! node(sn)%row(problem%n+1:ndof_fe_node,1): FEM equation for theta_k
!~                       !
!~                       ! Index of variables:
!~                       ! node(sn)%col(1:problem%n,1): u_k
!~                       ! node(sn)%col(problem%n+1:ndof_fe_node,1): FEM equation for theta_k
!~                       !
!~                       allocate (node(sn)%row(node(sn)%n_dof,1))
!~                       allocate (node(sn)%col(node(sn)%n_dof,1))
!~                       allocate (node(sn)%value_c(node(sn)%n_dof,1))
!~                       allocate (node(sn)%dvda_c(node(sn)%n_dof,1,problem%n_designvariables))
!~                       ! Initialize
!~                       node(sn)%row=0
!~                       node(sn)%col=0
!~                       ! Assign row and column to each degree of freedom
!~                       do k=1,ndof_fe_node
!~                         select case (node(sn)%ctype(k,1))
!~                           ! The displacement is already known
!~                           case (0)
!~                           ! The displacement is unknown
!~                           case (1)
!~                             ! FEM equation
!~                             node(sn)%row(k,1)=row
!~                             ! Increment counter
!~                             row=row+1
!~                             ! Degree of freedom k
!~                             node(sn)%col(k,1)=col
!~                             ! Increment counter
!~                             col=col+1
!~                         end select
!~                       end do
!~                     end if
!~                   end do

!~                 ! SOLID
!~                 case (3)
!~                   stop 'not yet 42' ! pero esta parte es igual que en 2D
!~                   !
!~                   ! Index of equations:
!~                   ! node(sn)%row(1:problem%n,1): FEM equation for u_k
!~                   !
!~                   ! Index of variables:
!~                   ! node(sn)%col(1:problem%n,1): u_k
!~                   !
!~                   do kn=1,element(se)%n_nodes
!~                     sn=element(se)%node(kn)
!~                     if (.not.node_used(sn)) then
!~                       node_used(sn)=.true.
!~                       node(sn)%n_dof=problem%n
!~                       allocate (node(sn)%row    (node(sn)%n_dof,1))
!~                       allocate (node(sn)%col    (node(sn)%n_dof,1))
!~                       allocate (node(sn)%value_c(node(sn)%n_dof,1))
!~                       allocate (node(sn)%dvda_c(node(sn)%n_dof,1,problem%n_designvariables))
!~                       ! Initialize
!~                       node(sn)%row=0
!~                       node(sn)%col=0
!~                       ! Assign row and column to each degree of freedom depending on the B.C.
!~                       do k=1,node(sn)%n_dof
!~                         select case (node(sn)%ctype(k,1))
!~                           ! The displacement is known
!~                           case (0)
!~                           ! The displacement is unknown
!~                           case (1)
!~                             ! FEM equation
!~                             node(sn)%row(k,1)=row
!~                             ! Increment counter
!~                             row=row+1
!~                             ! Degree of freedom k
!~                             node(sn)%col(k,1)=col
!~                             ! Increment counter
!~                             col=col+1
!~                         end select
!~                       end do
!~                     end if
!~                   end do

!~               end select
!~             end select
!~           end do
!~         end do
!~       end if
!~     end if
!~   end do

  ! Check if row==col
  if (row.ne.col) then
    write(output_unit,'(a,i8)') ' Rows   : ', row-1
    write(output_unit,'(a,i8)') ' Columns: ', col-1
    call fbem_error_message(error_unit,0,'fatal',0,'the mapping of the linear system of equations is wrong')
  end if

  ! Number of degrees of freedom
  n_dof=row-1
  if (n_dof.eq.0) stop 'n_dof=0'
  ! Print
  if (verbose_level.ge.1) then
    write(fmtstr,*) '(1x,a,i',fbem_nchar_int(n_dof),')'
    call fbem_trimall(fmtstr)
    write(output_unit,fmtstr) 'Number of degrees of freedom: ', n_dof
  end if

  ! Allocate and initialized variables for system of equations manipulations
  if (max_memory.ne.0) then
    memory=n_dof
    if (problem%sensitivity) then
      memory=16*((1+problem%n_designvariables)*memory+memory**2)
    else
      memory=16*(memory+memory**2)
    end if
    if (lse_condition.or.lse_refine) memory=2*memory
    if (memory.gt.max_memory) call fbem_error_message(error_unit,0,'memory',0,'required memory > memory limit')
  end if
  allocate (A_c(n_dof,n_dof),b_c(n_dof,1),fact_ipiv(n_dof),scal_r(n_dof),scal_c(n_dof))
  if (lse_condition.or.lse_refine) then
    allocate (Ao_c(n_dof,n_dof))
    Aodim=n_dof
  else
    allocate (Ao_c(1,1))
    Aodim=1
  end if
  if (problem%sensitivity) allocate (bsa_c(n_dof,problem%n_designvariables))

  ! =============================================================== !
  ! ALLOCATE INCIDENT FIELDS ON ELEMENTS, NODES AND INTERNAL POINTS !
  ! =============================================================== !

  do kr=1,n_regions
    if (region(kr)%class.eq.fbem_be) then
      ! BOUNDARIES
      do kb=1,region(kr)%n_boundaries
        sb=region(kr)%boundary(kb)
        sp=boundary(sb)%part
        select case (boundary(sb)%coupling)
          case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
            select case (region(kr)%type)
              case (fbem_potential)
                k_start=1
                k_end  =2
              case (fbem_viscoelastic)
                k_start=1
                k_end  =2*problem%n
              case (fbem_poroelastic)
                k_start=0
                k_end  =1+2*problem%n
            end select
            select case (boundary(sb)%class)
              case (fbem_boundary_class_ordinary)
                n_faces=1
              case (fbem_boundary_class_cracklike)
                n_faces=2
            end select
          case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
            select case (region(boundary(sb)%region(1))%type)
              case (fbem_potential)
                k_start=1
                k_end  =2
              case (fbem_viscoelastic)
                k_start=1
                k_end  =2*problem%n
              case (fbem_poroelastic)
                k_start=0
                k_end  =1+2*problem%n
            end select
            select case (region(boundary(sb)%region(2))%type)
              case (fbem_potential)
                k_start2=1
                k_end2  =2
              case (fbem_viscoelastic)
                k_start2=1
                k_end2  =2*problem%n
              case (fbem_poroelastic)
                k_start2=0
                k_end2  =1+2*problem%n
            end select
            k_start=min(k_start,k_start2)
            k_end=max(k_end,k_end2)
            n_faces=2
        end select
        do ke=1,part(sp)%n_elements
          se=part(sp)%element(ke)
          if (.not.allocated(element(se)%incident_c)) then
            allocate (element(se)%incident_c(k_start:k_end,element(se)%n_nodes,n_faces))
          end if
          do kn=1,element(se)%n_nodes
            sn=element(se)%node(kn)
            if (.not.allocated(node(sn)%incident_c)) then
              allocate (node(sn)%incident_c(k_start:k_end,n_faces))
            end if
          end do
        end do
      end do
      ! BODY LOADS
      do kb=1,region(kr)%n_be_bodyloads
        sb=region(kr)%be_bodyload(kb)
        sp=be_bodyload(sb)%part
        !
        ! The same code as BE element, except that n_faces=1
        !
        select case (region(kr)%type)
          case (fbem_potential)
            k_start=1
            k_end  =2
          case (fbem_viscoelastic)
            k_start=1
            k_end  =2*problem%n
          case (fbem_poroelastic)
            k_start=0
            k_end  =1+2*problem%n
        end select
        n_faces=1
        do ke=1,part(sp)%n_elements
          se=part(sp)%element(ke)
          if (.not.allocated(element(se)%incident_c)) then
            allocate (element(se)%incident_c(k_start:k_end,element(se)%n_nodes,n_faces))
          end if
          do kn=1,element(se)%n_nodes
            sn=element(se)%node(kn)
            if (.not.allocated(node(sn)%incident_c)) then
              allocate (node(sn)%incident_c(k_start:k_end,n_faces))
            end if
          end do
        end do
      end do
    end if
  end do

  !
  ! Allocate data for internal points
  !
  ! Description of value_c (and incident_c, and dvda):
  ! FLUID
  ! internalpoint(kip)%value_c(1,0): pressure                 : p
  ! internalpoint(kip)%value_c(1,i): normal displacement      : Un with normal n=e_i
  ! ELASTIC
  ! internalpoint(kip)%value_c(k,0): solid displacement       : u_k
  ! internalpoint(kip)%value_c(k,i): solid traction           : t_k with normal n=e_i
  ! POROELASTIC
  ! internalpoint(kip)%value_c(0,0): fluid equivalente stress : tau
  ! internalpoint(kip)%value_c(k,0): solid displacement       : u_k
  ! internalpoint(kip)%value_c(0,i): fluid normal displacement: Un with normal n=e_i
  ! internalpoint(kip)%value_c(k,i): solid traction           : t_k with normal n=e_i
  !
  do k=1,n_internalpoints
    select case (region(internalpoint(k)%region)%type)
      case (fbem_potential)
          k_start=1
          k_end  =1
      case (fbem_viscoelastic)
          k_start=1
          k_end  =problem%n
      case (fbem_poroelastic)
          k_start=0
          k_end  =problem%n
    end select
    allocate (internalpoint(k)%value_c   (k_start:k_end,0:problem%n))
    allocate (internalpoint(k)%incident_c(k_start:k_end,0:problem%n))
    allocate (internalpoint(k)%dvda_c    (k_start:k_end,0:problem%n,problem%n_designvariables))
  end do
  !
  ! Allocate data for internal elements
  !
  if (internalelements) then


    do kp=1,internalelements_mesh%n_parts
      kr=internalelements_mesh%part(kp)%entity
      if (kr.eq.0) cycle
      if (region(kr)%class.eq.fbem_be) then
        select case (region(kr)%type)
          case (fbem_potential)
              k_start=1
              k_end  =1
          case (fbem_viscoelastic)
              k_start=1
              k_end  =problem%n
          case (fbem_poroelastic)
              k_start=0
              k_end  =problem%n
        end select
        do ke=1,internalelements_mesh%part(kp)%n_elements
          se=internalelements_mesh%part(kp)%element(ke)

          allocate (internalelements_mesh%element(se)%value_c(k_start:k_end,internalelements_mesh%element(se)%n_nodes,0:problem%n))

          ! falta incident_c
          ! Falta allocatar en los nodos, para tambien sacar las soluciones por nodos (continuas)



        end do
      end if
    end do


  end if

  if (verbose_level.ge.1)  call fbem_timestamp_w_message(output_unit,2,'END building auxiliary variables')

end subroutine build_auxiliary_variables_mechanics_harmonic
