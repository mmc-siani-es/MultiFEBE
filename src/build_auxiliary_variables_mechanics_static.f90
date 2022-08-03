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

subroutine build_auxiliary_variables_mechanics_static

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
  integer                    :: kr, kp, kb, ke, kn, ks, ss, sp, kn2
  integer                    :: sb, se, sn, snm
  integer                    :: row, col
  integer                    :: k, kt
  integer                    :: ndof_fe_node
  logical, allocatable       :: node_used(:)
  character(len=fbem_fmtstr) :: fmtstr
  integer(kind=int64)        :: memory

  if (verbose_level.ge.1)  call fbem_timestamp_w_message(output_unit,2,'START building auxiliary variables')

  ! Allocate auxiliary variable
  allocate (node_used(n_nodes))

  ! Initialize LSE counters
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

                ! ------------------------------------------------------------------------------------------------------------------
                ! BE BOUNDARY
                ! ------------------------------------------------------------------------------------------------------------------

                case (fbem_boundary_coupling_be)
                  select case (boundary(sb)%class)

                    ! =================
                    ! ORDINARY BOUNDARY
                    ! =================

                    case (fbem_boundary_class_ordinary)
                      !
                      ! Index of equations for each coordinate k:
                      ! node(sn)%row(k,1): SBIE, HBIE or SBIE+beta*HBIE
                      !
                      ! Index of variables for each coordinate k (k=1,2 (2D), k=1,2,3 (3D)):
                      ! node(sn)%col(          k,1): u_k
                      ! node(sn)%col(problem%n+k,1): t_k
                      !
                      allocate (node(sn)%row(problem%n,1))
                      allocate (node(sn)%col(2*problem%n,1))
                      allocate (node(sn)%value_r(2*problem%n,1))
                      allocate (node(sn)%dvda_r(2*problem%n,1,problem%n_designvariables))
                      ! Initialize
                      do k=1,problem%n
                        node(sn)%row(k,1)=0
                        node(sn)%col(          k,1)=0
                        node(sn)%col(problem%n+k,1)=0
                      end do
                      ! Assign row to the equation and column to the variables
                      do k=1,problem%n
                        ! Equation
                        node(sn)%row(k,1)=row
                        ! Increment counter
                        row=row+1
                        ! Assign values depending on the boundary condition
                        select case (node(sn)%ctype(k,1))
                          ! u_k known, p_k unknown
                          case (0)
                            node(sn)%col(problem%n+k,1)=col
                          ! p_k known, u_k unknown
                          case (1,10)
                            node(sn)%col(          k,1)=col
                        end select
                        ! Increment counter
                        col=col+1
                      end do

                    ! ===================
                    ! CRACK-LIKE BOUNDARY
                    ! ===================

                    case (fbem_boundary_class_cracklike)
                      !
                      ! Index of equations for each coordinate k (Dual BEM):
                      ! node(sn)%row(k,1): SBIE
                      ! node(sn)%row(k,2): HBIE
                      !
                      ! Index of variables for each coordinate k:
                      ! node(sn)%col(          k,1): u_k for face +
                      ! node(sn)%col(problem%n+k,1): p_k for face +
                      ! node(sn)%col(          k,2): u_k for face -
                      ! node(sn)%col(problem%n+k,2): p_k for face -
                      !
                      allocate (node(sn)%row(problem%n,2))
                      allocate (node(sn)%col(2*problem%n,2))
                      allocate (node(sn)%value_r(2*problem%n,2))
                      allocate (node(sn)%dvda_r(2*problem%n,2,problem%n_designvariables))
                      ! Initialize
                      do k=1,problem%n
                        node(sn)%row(k,1)=0
                        node(sn)%row(k,2)=0
                        node(sn)%col(          k,1)=0
                        node(sn)%col(problem%n+k,1)=0
                        node(sn)%col(          k,2)=0
                        node(sn)%col(problem%n+k,2)=0
                      end do
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
                          ! u_k known, p_k unknown
                          case (0)
                            node(sn)%col(problem%n+k,1)=col
                          ! p_k known, u_k unknown
                          case (1)
                            node(sn)%col(          k,1)=col
                        end select
                        ! Increment counter
                        col=col+1
                        ! Face -
                        ! Assign values depending on the boundary condition
                        select case (node(sn)%ctype(k,2))
                          ! u_k known, p_k unknown
                          case (0)
                            node(sn)%col(problem%n+k,2)=col
                          ! p_k known, u_k unknown
                          case (1)
                            node(sn)%col(          k,2)=col
                        end select
                        ! Increment counter
                        col=col+1
                      end do
                  end select

                ! ------------------------------------------------------------------------------------------------------------------
                ! BE-BE BOUNDARY
                ! ------------------------------------------------------------------------------------------------------------------
                ! The region 1 is the region where the boundary has N+
                ! The region 2 is the region where the boundary has N-

                case (fbem_boundary_coupling_be_be)
                  !
                  ! Index of equations for each coordinate k:
                  ! node(sn)%row(k,1): SBIE, HBIE or SBIE+beta*HBIE for region 1
                  ! node(sn)%row(k,2): SBIE, HBIE or SBIE+beta*HBIE for region 2
                  !
                  ! Index of variables for each coordinate k:
                  ! node(sn)%col(          k,1): u_k for region 1
                  ! node(sn)%col(problem%n+k,1): p_k for region 1
                  ! node(sn)%col(          k,2): u_k for region 2 (not used since u_k_{region 1}=u_k_{region 2} is assumed)
                  ! node(sn)%col(problem%n+k,2): p_k for region 2 (not used since p_k_{region 1}=-p_k_{region 2} is assumed)
                  !
                  allocate (node(sn)%row(problem%n,2))
                  allocate (node(sn)%col(2*problem%n,2))
                  allocate (node(sn)%value_r(2*problem%n,2))
                  allocate (node(sn)%dvda_r(2*problem%n,2,problem%n_designvariables))
                  ! Initialize
                  do k=1,problem%n
                    node(sn)%row(k,1)=0
                    node(sn)%row(k,2)=0
                    node(sn)%col(          k,1)=0
                    node(sn)%col(problem%n+k,1)=0
                    node(sn)%col(          k,2)=0
                    node(sn)%col(problem%n+k,2)=0
                  end do
                  ! Assign row to the equations and column to the variables
                  do k=1,problem%n
                    ! Equation for region 1
                    node(sn)%row(k,1)=row
                    ! Equation for region 2
                    node(sn)%row(k,2)=row+1
                    ! Increment counter
                    row=row+2
                    ! Only the variables of region 1 are active
                    node(sn)%col(          k,1)=col
                    node(sn)%col(problem%n+k,1)=col+1
                    ! Increment counter
                    col=col+2
                  end do

                ! ------------------------------------------------------------------------------------------------------------------
                ! BE-FE BOUNDARY
                ! ------------------------------------------------------------------------------------------------------------------

                case (fbem_boundary_coupling_be_fe)
                  select case (boundary(sb)%class)

                    ! =================
                    ! ORDINARY BOUNDARY
                    ! =================

                    case (fbem_boundary_class_ordinary)
                      !
                      ! Index of equations for each coordinate k:
                      ! node(sn)%row(k,1): SBIE, HBIE or SBIE+beta*HBIE
                      !
                      ! Index of variables for each coordinate k:
                      ! node(sn)%col(          k,1): u_k
                      ! node(sn)%col(problem%n+k,1): t_k
                      !
                      allocate (node(sn)%row(problem%n,1))
                      allocate (node(sn)%col(2*problem%n,1))
                      allocate (node(sn)%value_r(2*problem%n,1))
                      allocate (node(sn)%dvda_r(2*problem%n,1,problem%n_designvariables))
                      ! Initialize
                      do k=1,problem%n
                        node(sn)%row(k,1)=0
                        node(sn)%col(          k,1)=0
                        node(sn)%col(problem%n+k,1)=0
                      end do
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

                    ! ===================
                    ! CRACK-LIKE BOUNDARY
                    ! ===================

                    case (fbem_boundary_class_cracklike)

!                      ! ------------------------------------------------------------------------------------- !
!                      ! SBIE (solve for Σt = t^+ + t^-) + HBIE (calculate t^+ and t^- ) (for perfect bonding) !
!                      ! ------------------------------------------------------------------------------------- !

!                      !
!                      ! Index of equations for each coordinate k:
!                      ! node(sn)%row(k,1): SBIE
!                      !
!                      ! Index of variables for each coordinate k:
!                      ! node(sn)%col(          k,1): inactive
!                      ! node(sn)%col(problem%n+k,1): Σt_k = t_k^+ + t_k^-
!                      !
!                      ! Index of variable values for each coordinate k:
!                      !
!                      ! node(sn)%value_r(          k,1): u_k^+
!                      ! node(sn)%value_r(problem%n+k,1): t_k^+
!                      ! node(sn)%value_r(          k,2): u_k^-
!                      ! node(sn)%value_r(problem%n+k,2): t_k^-
!                      !
!                      ! Index of variable sensitivity values for each coordinate k:
!                      !
!                      ! node(sn)%dvda_r(          k,1,j): ∂u_k^+/∂a_j
!                      ! node(sn)%dvda_r(problem%n+k,1,j): ∂t_k^+/∂a_j
!                      ! node(sn)%dvda_r(          k,2,j): ∂u_k^-/∂a_j
!                      ! node(sn)%dvda_r(problem%n+k,2,j): ∂t_k^-/∂a_j
!                      !
!                      allocate (node(sn)%row(problem%n,1))
!                      allocate (node(sn)%col(2*problem%n,1))
!                      allocate (node(sn)%value_r(2*problem%n,2))
!                      allocate (node(sn)%dvda_r(2*problem%n,2,problem%n_designvariables))
!                      ! Initialize
!                      node(sn)%row=0
!                      node(sn)%col=0
!                      ! Assign row to the equation and column to the variables
!                      do k=1,problem%n
!                        ! SBIE
!                        node(sn)%row(k,1)=row
!                        row=row+1
!                        ! Σt_k
!                        node(sn)%col(problem%n+k,1)=col
!                        col=col+1
!                      end do

                      ! ------------------------------------------------------------- !
                      ! General model using Dual BEM (for general contact conditions) !
                      ! ------------------------------------------------------------- !

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
                      allocate (node(sn)%value_r(2*problem%n,2))
                      allocate (node(sn)%dvda_r(2*problem%n,2,problem%n_designvariables))
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
                        ! Assign values depending on the boundary condition.
                        ! t_k for face + is always unknown
                        node(sn)%col(problem%n+k,1)=col
                        ! t_k for face - is always unknown
                        node(sn)%col(problem%n+k,2)=col+1
                        ! Increment counter
                        col=col+2
                        ! Perfect bonding: u_k^+(E)=u_k^-(E)=u_k^(FE)
                      end do

                  end select

                ! ------------------------------------------------------------------------------------------------------------------
                ! BE-FE-BE BOUNDARY
                ! ------------------------------------------------------------------------------------------------------------------

                case (fbem_boundary_coupling_be_fe_be)
                  !
                  ! Index of equations for each coordinate k:
                  ! node(sn)%row(k,1): SBIE, HBIE or SBIE+beta*HBIE for region 1
                  ! node(sn)%row(k,2): SBIE, HBIE or SBIE+beta*HBIE for region 2
                  !
                  ! Index of variables for each coordinate k:
                  ! node(sn)%col(          k,1): u_k for region 1
                  ! node(sn)%col(problem%n+k,1): t_k for region 1
                  ! node(sn)%col(          k,2): u_k for region 2 (not used since u_k_{region 1}=u_k_{region 2})
                  ! node(sn)%col(problem%n+k,2): t_k for region 2
                  !
                  allocate (node(sn)%row    (problem%n,2))
                  allocate (node(sn)%col    (2*problem%n,2))
                  allocate (node(sn)%value_c(2*problem%n,2))
                  ! Initialize
                  do k=1,problem%n
                    node(sn)%row(k,1)=0
                    node(sn)%row(k,2)=0
                    node(sn)%col(          k,1)=0
                    node(sn)%col(problem%n+k,1)=0
                    node(sn)%col(          k,2)=0
                    node(sn)%col(problem%n+k,2)=0
                  end do
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
            select case (problem%n)
              !
              ! 2D: LINE LOAD AT THE BEAM TIP
              !
              case (2)
                do ke=1,part(sp)%n_elements
                  se=part(sp)%element(ke)
                  do kn=1,element(se)%n_nodes
                    sn=element(se)%node(kn)
                    if (node_used(sn).eqv.(.false.)) then
                      node_used(sn)=.true.

                      !
                      ! Falta saber qué hacer con los en cada coordenada cuando la punta fem tiene condiciones de contorno
                      !

                      !
                      ! Index of equations:
                      ! - node(sn)%row(k,nf): SBIE for each coordinate k and each functional node (thickness) associated with node sn
                      !
                      ! Index of variables:
                      ! - node(sn)%col(          k,nf): u_k for each coordinate k and each functional node (thickness) associated with node sn (not active since equal to tip/edge displacement field)
                      ! - node(sn)%col(problem%n+k,nf): b_k for each coordinate k and each functional node (thickness) associated with node sn
                      !
                      allocate (node(sn)%row    (  problem%n,fbem_n_nodes(element(se)%type_f2)))
                      allocate (node(sn)%col    (2*problem%n,fbem_n_nodes(element(se)%type_f2)))
                      allocate (node(sn)%value_r(2*problem%n,fbem_n_nodes(element(se)%type_f2)))
                      ! Initialize
                      node(sn)%row=0
                      node(sn)%col=0
                      ! Assign row to the equation and column to the variables
                      do kt=1,fbem_n_nodes(element(se)%type_f2)/element(se)%n_nodes
                        do k=1,problem%n
                          ! Equation
                          node(sn)%row(k,kt)=row
                          ! Increment counter
                          row=row+1
                          ! u_k=u_k from tip/edge displacement field, b_k unknown
                          node(sn)%col(problem%n+k,kt)=col
                          ! Increment counter
                          col=col+1
                        end do
                      end do


                    end if
                  end do
                end do
              !
              ! 3D: SURFACE LOAD AT THE BEAM TIP
              !
              case (3)
                ! Modelo de Luis (en principio solo habria que crear un cuadrilatero que emule a una circunferencia, se puede hacer..)
                stop 'not yet 43'
            end select

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
                  allocate (node(sn)%row    (problem%n,1))
                  allocate (node(sn)%col    (2*problem%n,1))
                  allocate (node(sn)%value_r(2*problem%n,1))
                  allocate (node(sn)%dvda_r (2*problem%n,1,problem%n_designvariables))
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
            do ke=1,part(sp)%n_elements
              se=part(sp)%element(ke)
              do kn=1,element(se)%n_nodes
                sn=element(se)%node(kn)
                if (node_used(sn).eqv.(.false.)) then
                  node_used(sn)=.true.

                  !
                  ! Falta saber qué hacer con los en cada coordenada cuando la punta fem tiene condiciones de contorno
                  !
                  !
                  ! Index of equations:
                  ! - node(sn)%row(k,nf): SBIE for each coordinate k and each functional node in the thickness direction (nf) associated with node sn
                  !
                  ! Index of variables:
                  ! - node(sn)%col(          k,nf): u_k for each coordinate k and each functional node in the thickness direction (nf) associated with node sn (not active since equal to tip/edge displacement field)
                  ! - node(sn)%col(problem%n+k,nf): b_k for each coordinate k and each functional node in the thickness direction (nf) associated with node sn
                  !
                  allocate (node(sn)%row    (  problem%n,fbem_n_nodes(element(se)%type_f2)/element(se)%n_nodes))
                  allocate (node(sn)%col    (2*problem%n,fbem_n_nodes(element(se)%type_f2)/element(se)%n_nodes))
                  allocate (node(sn)%value_r(2*problem%n,fbem_n_nodes(element(se)%type_f2)/element(se)%n_nodes))
                  ! Initialize
                  node(sn)%row=0
                  node(sn)%col=0
                  ! Assign row to the equation and column to the variables
                  do kt=1,fbem_n_nodes(element(se)%type_f2)/element(se)%n_nodes
                    do k=1,problem%n
                      ! Equation
                      node(sn)%row(k,kt)=row
                      ! Increment counter
                      row=row+1
                      ! u_k^i=u_k^{Shell FE} from edge displacement field, b_k unknown
                      node(sn)%col(problem%n+k,kt)=col
                      ! Increment counter
                      col=col+1
                    end do
                  end do

                end if
              end do
            end do

        end select

      end do

    end if
  end do

  ! ================================================================================================================================
  ! FE REGIONS
  ! ================================================================================================================================

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
                    ! Stress resultants (local axes v_k): value_r(k,n,1): NX,VY,BZ or NX,VY,VZ,BX,BY,BZ
                    ! Equilibrating loads (global axes) : value_r(k,n,2): FX,FY,MZ or FX,FY,FZ,MX,MY,MZ
                    !
                    allocate (element(se)%value_r(3*(problem%n-1),element(se)%n_nodes,2))
                    element(se)%value_r=0
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
                    ! Stress resultants (local axes ep_k): value_r(k,n,1): NX,VY,BZ OR NX,VY,VZ,BX,BY,BZ
                    ! Equilibrating loads (global axes)  : value_r(k,n,2): FX,FY[,MZ] OR FX,FY,FZ[,MX,MY,MZ]
                    !
                    allocate (element(se)%value_r(3*(problem%n-1),element(se)%n_nodes,2))
                    element(se)%value_r=0
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
                    ! Stress resultants                  : value_r(k,n,1): NX
                    ! Equilibrating loads (global axes)  : value_r(k,n,2): FX,FY,FZ
                    !
                    allocate (element(se)%value_r(problem%n,element(se)%n_nodes,2))
                    element(se)%value_r=0
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
                    ! Spring forces                     : value_r(k,n,1): NX,NY,NZ
                    ! Equilibrating loads (nodal axes)  : value_r(k,n,2): FX,FY,FZ
                    !
                    allocate (element(se)%value_r(problem%n,2,2))
                    element(se)%value_r=0
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
                    ! Spring forces/moments             : value_r(k,n,1): NX,NY,NZ,BX,BY,BZ
                    ! Equilibrating loads (nodal axes)  : value_r(k,n,2): FX,FY,FZ,MX,MY,MZ
                    !
                    allocate (element(se)%value_r(3*(problem%n-1),2,2))
                    element(se)%value_r=0
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
                    ! Stress tensor (global axes)      : value_r(k,n,1): SXX,SYY,SXY
                    ! Equilibrating loads (global axes): value_r(k,n,2): FX,FY
                    !
                    allocate (element(se)%value_r(3,element(se)%n_nodes,2))
                    element(se)%value_r=0
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
                    ! Stress resultants (local axes ep_k)               : value_r(k,n,1): NX,NY,NXY,MX,MY,MXY,VX,VY
                    ! Equilibrating loads (global u, global/local theta): value_r(k,n,2): FX,FY,FZ[,MA,MB][,MX,MY,MZ]
                    !
                    allocate (element(se)%value_r(8,element(se)%n_nodes,2))
                    element(se)%value_r=0
                    !
                    ! DOF handling
                    !
                    allocate(element(se)%node_n_dof(element(se)%n_nodes))
                    do kn=1,element(se)%n_nodes
                      sn=element(se)%node(kn)
                      if (node(sn)%n_dof.eq.0) then
                        element(se)%node_n_dof(kn)=5
                        if (node(sn)%is_singular.or.(node(sn)%rigid_link.ne.0)) then
                          element(se)%node_n_dof(kn)=6
                        end if
                      else
                        element(se)%node_n_dof(kn)=node(sn)%n_dof
                        select case (node(sn)%n_dof)
                          case (5)
                            if (node(sn)%is_singular.or.(node(sn)%rigid_link.ne.0)) then
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
                ! Stress tensor (global axes)      : value_r(k,n,1): SXX,SYY,SZZ,SXY,SXZ,SYZ
                ! Equilibrating loads (global axes): value_r(k,n,2): FX,FY,FZ


                !
                ! NOTA: es mas logico poner las cargas equilibrantes en el indice 1
                !


                !
                allocate (element(se)%value_r(6,element(se)%n_nodes,2))
                element(se)%value_r=0
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
        ! Nodal gen. displacements: value_r(k,1): UX,UY,RZ OR UX,UY,UZ,RX,RY,RZ
        ! Nodal reactions/loads   : value_r(k,2): FX,FY,MZ OR FX,FY,FZ,MX,MY,MZ
        !
        allocate (node(sn)%value_r(node(sn)%n_dof,2))
        node(sn)%value_r=0
        allocate (node(sn)%dvda_r (node(sn)%n_dof,2,problem%n_designvariables))
        node(sn)%dvda_r=0
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
              ! Find the maximum number of DOFs between all elements that share the node
              ndof_fe_node=0
              do ke=1,node(sn)%n_elements
                se=node(sn)%element(ke)
                kn2=node(sn)%element_node_iid(ke)
                if (ndof_fe_node.lt.element(se)%node_n_dof(kn2)) then
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
              ! Nodal gen. displacements: value_r(k,1): UX,UY[,RZ] OR UX,UY,UZ[,RA,RB][,RX,RY,RZ]
              ! Nodal reactions/loads   : value_r(k,2): FX,FY[,MZ] OR FX,FY,FZ[,MA,MB][,MX,MY,MZ]
              !
              allocate (node(sn)%value_r(node(sn)%n_dof,2))
              node(sn)%value_r=0
              allocate (node(sn)%dvda_r (node(sn)%n_dof,2,problem%n_designvariables))
              node(sn)%dvda_r=0
            end if
          end do








          !
          ! Old code
          !

!~           do ke=1,part(sp)%n_elements

!~             se=part(sp)%element(ke)

!~             select case (element(se)%n_dimension)

!~               ! ====================================================================================================================
!~               ! ONE-DIMENSIONAL ELEMENTS
!~               ! ====================================================================================================================

!~               case (1)

!~                 select case (element(se)%fe_type)

!~                   case (0)

!~                     !---------------------------------------------------------------------------------------------------------------
!~                     ! DEGENERATED BEAM FINITE ELEMENT
!~                     !
!~                     ! Index of equations:
!~                     ! node(sn)%row(1,1): u_1
!~                     ! node(sn)%row(2,1): u_2
!~                     ! node(sn)%row(3,1): u_3
!~                     ! node(sn)%row(4,1): theta_1
!~                     ! node(sn)%row(5,1): theta_2
!~                     ! node(sn)%row(6,1): theta_3
!~                     !
!~                     ! Index of variables:
!~                     ! node(sn)%col(1,1): u_1
!~                     ! node(sn)%col(2,1): u_2
!~                     ! node(sn)%col(3,1): u_3
!~                     ! node(sn)%col(4,1): theta_1
!~                     ! node(sn)%col(5,1): theta_2
!~                     ! node(sn)%col(6,1): theta_3
!~                     !
!~                     do kn=1,element(se)%n_nodes
!~                       sn=element(se)%node(kn)
!~                       if (.not.node_used(sn)) then
!~                         node_used(sn)=.true.
!~                         !
!~                         ! DOF handling
!~                         !
!~                         node(sn)%n_dof=element(se)%node_n_dof(kn)
!~                         allocate (node(sn)%row(node(sn)%n_dof,1))
!~                         allocate (node(sn)%col(node(sn)%n_dof,1))
!~                         ! Initialize
!~                         node(sn)%row=0
!~                         node(sn)%col=0
!~                         ! If it is a slave node of a rigid region, it does not have its own active DOF
!~                         if (node(sn)%master.ne.0) cycle
!~                         ! Assign row and column to each degree of freedom depending on the B.C.
!~                         do k=1,node(sn)%n_dof
!~                           select case (node(sn)%ctype(k,1))
!~                             ! The displacement is already known
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
!~                         !
!~                         ! Values to store
!~                         !
!~                         !
!~                         ! Nodal gen. displacements: value_r(k,1): UX,UY,RZ OR UX,UY,UZ,RX,RY,RZ
!~                         ! Nodal reactions/loads   : value_r(k,2): FX,FY,MZ OR FX,FY,FZ,MX,MY,MZ
!~                         !
!~                         allocate (node(sn)%value_r(node(sn)%n_dof,2))
!~                         node(sn)%value_r=0
!~                         allocate (node(sn)%dvda_r (node(sn)%n_dof,2,problem%n_designvariables))
!~                         node(sn)%dvda_r=0
!~                         !
!~                       end if
!~                     end do
!~                     !
!~                     !---------------------------------------------------------------------------------------------------------------

!~                   case (1,2)

!~                     !---------------------------------------------------------------------------------------------------------------
!~                     ! STRAIGHT EULER-BERNOULLI AND TIMOSHENKO BEAM FINITE ELEMENTS
!~                     !
!~                     ! Index of equations:
!~                     ! node(sn)%row(             1:problem%n,1): FEM equation for u_k
!~                     ! node(sn)%row(problem%n+1:ndof_fe_node,1): FEM equation for theta_k
!~                     !
!~                     ! Index of variables:
!~                     ! node(sn)%col(             1:problem%n,1): u_k
!~                     ! node(sn)%col(problem%n+1:ndof_fe_node,1): theta_k
!~                     !
!~                     do kn=1,element(se)%n_nodes
!~                       sn=element(se)%node(kn)
!~                       if (.not.node_used(sn)) then
!~                         node_used(sn)=.true.
!~                         !
!~                         ! DOF handling
!~                         !
!~                         node(sn)%n_dof=element(se)%node_n_dof(kn)
!~                         allocate (node(sn)%row(node(sn)%n_dof,1))
!~                         allocate (node(sn)%col(node(sn)%n_dof,1))
!~                         ! Initialize
!~                         node(sn)%row=0
!~                         node(sn)%col=0
!~                         ! If it is a slave node of a rigid region, it does not have its own active DOF
!~                         if (node(sn)%master.ne.0) cycle
!~                         ! Assign row and column to each degree of freedom depending on the B.C.
!~                         do k=1,node(sn)%n_dof
!~                           select case (node(sn)%ctype(k,1))
!~                             ! The displacement is already known
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
!~                         !
!~                         ! Values to store
!~                         !
!~                         !
!~                         ! Nodal gen. displacements: value_r(k,1): UX,UY[,RZ] OR UX,UY,UZ[,RX,RY,RZ]
!~                         ! Nodal reactions/loads   : value_r(k,2): FX,FY[,MZ] OR FX,FY,FZ[,MX,MY,MZ]
!~                         !
!~                         allocate (node(sn)%value_r(node(sn)%n_dof,2))
!~                         node(sn)%value_r=0
!~                         allocate (node(sn)%dvda_r (node(sn)%n_dof,2,problem%n_designvariables))
!~                         node(sn)%dvda_r=0
!~                         !
!~                       end if
!~                     end do
!~                     !
!~                     !---------------------------------------------------------------------------------------------------------------


!~                   case default

!~                     !---------------------------------------------------------------------------------------------------------------
!~                     ! OTHER TYPES
!~                     !
!~                     !
!~                     ! AQUI FALTA METER LOS DISTRA, DISROTRA, Y BAR ELEMENTS
!~                     !
!~                     call fbem_error_message(error_unit,0,'element',element(se)%id,'invalid type of 1D element')
!~                     !
!~                     !---------------------------------------------------------------------------------------------------------------

!~                 end select

!~               ! ====================================================================================================================
!~               ! TWO-DIMENSIONAL ELEMENTS
!~               ! ====================================================================================================================

!~               case (2)

!~                 select case (problem%n)

!~                   case (2)

!~                     !---------------------------------------------------------------------------------------------------------------
!~                     ! SOLID / CONTINUUM ELEMENTS
!~                     !
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
!~                         !
!~                         ! DOF handling
!~                         !
!~                         node(sn)%n_dof=element(se)%node_n_dof(kn)
!~                         allocate (node(sn)%row(node(sn)%n_dof,1))
!~                         allocate (node(sn)%col(node(sn)%n_dof,1))
!~                         ! Initialize
!~                         node(sn)%row=0
!~                         node(sn)%col=0
!~                         ! If it is a slave node of a rigid region, it does not have its own active DOF
!~                         if (node(sn)%master.ne.0) cycle
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
!~                         !
!~                         ! Values to store
!~                         !
!~                         !
!~                         ! Nodal gen. displacements: value_r(k,1): UX,UY
!~                         ! Nodal reactions/loads   : value_r(k,2): FX,FY
!~                         !
!~                         allocate (node(sn)%value_r(node(sn)%n_dof,2))
!~                         node(sn)%value_r=0
!~                         allocate (node(sn)%dvda_r (node(sn)%n_dof,2,problem%n_designvariables))
!~                         node(sn)%dvda_r=0
!~                         !
!~                       end if
!~                     end do
!~                     !
!~                     !---------------------------------------------------------------------------------------------------------------

!~                   case (3)

!~                     !---------------------------------------------------------------------------------------------------------------
!~                     ! DEGENERATED SHELL FINITE ELEMENT
!~                     !
!~                     !
!~                     ! Shell obtained from degenerated 3D solid
!~                     !
!~                     ! 2 types of nodes depending on the connectivity of each shell
!~                     !   Type 5: 5 DOF, 3 global axes displacements and 2 local axes rotations
!~                     !   Type 6: 6 DOF, 3 global axes displacements and 3 global axes rotations
!~                     !
!~                     !
!~                     ! Index of equations:
!~                     ! node(sn)%row(1:problem%n,1): FEM equation for u_k
!~                     ! node(sn)%row(problem%n+1:ndof_fe_node,1): FEM equation for theta_k
!~                     !
!~                     ! Index of variables:
!~                     ! node(sn)%col(1:problem%n,1): u_k
!~                     ! node(sn)%col(problem%n+1:ndof_fe_node,1): FEM equation for theta_k
!~                     !
!~                     do kn=1,element(se)%n_nodes
!~                       sn=element(se)%node(kn)
!~                       if (.not.node_used(sn)) then
!~                         node_used(sn)=.true.
!~                         !
!~                         ! DOF handling
!~                         !
!~                         node(sn)%n_dof=element(se)%node_n_dof(kn)

!~                         allocate (node(sn)%row(node(sn)%n_dof,1))
!~                         allocate (node(sn)%col(node(sn)%n_dof,1))
!~                         ! Initialize
!~                         node(sn)%row=0
!~                         node(sn)%col=0
!~                         ! If it is a slave node of a rigid region, it does not have its own active DOF
!~                         if (node(sn)%master.ne.0) cycle
!~                         ! Assign row and column to each degree of freedom
!~                         do k=1,node(sn)%n_dof
!~                           select case (node(sn)%ctype(k,1))
!~                             ! The displacement is already known
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
!~                         !
!~                         ! Values to store
!~                         !
!~                         !
!~                         ! Nodal gen. displacements: value_r(k,1): UX,UY,UZ,RA,RB OR UX,UY,UZ,RX,RY,RZ
!~                         ! Nodal reactions/loads   : value_r(k,2): FX,FY,FZ,MA,MB OR FX,FY,FZ,MX,MY,MZ
!~                         !
!~                         allocate (node(sn)%value_r(node(sn)%n_dof,2))
!~                         allocate (node(sn)%dvda_r (node(sn)%n_dof,2,problem%n_designvariables))
!~                         !
!~                       end if
!~                     end do
!~                     !
!~                     !---------------------------------------------------------------------------------------------------------------

!~                 end select

!~               ! ====================================================================================================================
!~               ! THREE-DIMENSIONAL ELEMENTS
!~               ! ====================================================================================================================

!~               case (3)

!~                 !-------------------------------------------------------------------------------------------------------------------
!~                 ! SOLID / CONTINUUM ELEMENTS
!~                 !
!~                 !
!~                 ! Index of equations:
!~                 ! node(sn)%row(1:problem%n,1): FEM equation for u_k
!~                 !
!~                 ! Index of variables:
!~                 ! node(sn)%col(1:problem%n,1): u_k
!~                 !
!~                 do kn=1,element(se)%n_nodes
!~                   sn=element(se)%node(kn)
!~                   if (.not.node_used(sn)) then
!~                     node_used(sn)=.true.
!~                     !
!~                     ! DOF handling
!~                     !
!~                     node(sn)%n_dof=3
!~                     allocate (node(sn)%row(node(sn)%n_dof,1))
!~                     allocate (node(sn)%col(node(sn)%n_dof,1))
!~                     ! Initialize
!~                     node(sn)%row=0
!~                     node(sn)%col=0
!~                     ! If it is a slave node of a rigid region, it does not have its own active DOF
!~                     if (node(sn)%master.ne.0) cycle
!~                     ! Assign row and column to each degree of freedom depending on the B.C.
!~                     do k=1,node(sn)%n_dof
!~                       select case (node(sn)%ctype(k,1))
!~                         ! The displacement is known
!~                         case (0)
!~                         ! The displacement is unknown
!~                         case (1)
!~                           ! FEM equation
!~                           node(sn)%row(k,1)=row
!~                           ! Increment counter
!~                           row=row+1
!~                           ! Degree of freedom k
!~                           node(sn)%col(k,1)=col
!~                           ! Increment counter
!~                           col=col+1
!~                       end select
!~                     end do
!~                     !
!~                     ! Values to store
!~                     !
!~                     !
!~                     ! Nodal gen. displacements: value_r(k,1): UX,UY,UZ
!~                     ! Nodal reactions/loads   : value_r(k,2): FX,FY,FZ
!~                     !
!~                     allocate (node(sn)%value_r(node(sn)%n_dof,2))
!~                     node(sn)%value_r=0
!~                     allocate (node(sn)%dvda_r (node(sn)%n_dof,2,problem%n_designvariables))
!~                     node(sn)%dvda_r=0
!~                     !
!~                   end if
!~                 end do
!~                 !
!~                 call fbem_error_message(error_unit,0,'element',element(se)%id,'3D elements not available yet')
!~                 !
!~                 !-------------------------------------------------------------------------------------------------------------------

!~             end select
!~           end do





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
            allocate (node(sn)%row    (node(sn)%n_dof,1))
            allocate (node(sn)%col    (node(sn)%n_dof,1))
            ! Initialize
            node(sn)%row=0
            node(sn)%col=0
            allocate (node(sn)%value_r(node(sn)%n_dof,2))
            node(sn)%value_r=0
            allocate (node(sn)%dvda_r (node(sn)%n_dof,2,problem%n_designvariables))
            node(sn)%dvda_r=0
          end if
        end do
      end if
    end if
  end do

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

  ! Allocate and initialize variables for system of equations manipulations
  if (max_memory.ne.0) then
    memory=n_dof
    if (problem%sensitivity) then
      memory=8*((1+problem%n_designvariables)*memory+memory**2)
    else
      memory=8*(memory+memory**2)
    end if
    if (lse_condition.or.lse_refine) memory=2*memory
    if (memory.gt.max_memory) call fbem_error_message(error_unit,0,'memory',0,'required memory > memory limit')
  end if
  allocate (A_r(n_dof,n_dof),b_r(n_dof,1),fact_ipiv(n_dof),scal_r(n_dof),scal_c(n_dof))
  if (lse_condition.or.lse_refine) then
    allocate (Ao_r(n_dof,n_dof))
    Aodim=n_dof
  else
    allocate (Ao_r(1,1))
    Aodim=1
  end if
  if (problem%sensitivity) allocate (bsa_r(n_dof,problem%n_designvariables))

  ! Allocate data for internal points
  do k=1,n_internalpoints
    allocate (internalpoint(k)%value_r(problem%n,0:problem%n))
    allocate (internalpoint(k)%dvda_r(problem%n,0:problem%n,problem%n_designvariables))
  end do

  !
  ! Allocate data for internal elements
  !
  if (internalelements) then
    do kp=1,internalelements_mesh%n_parts
      kr=internalelements_mesh%part(kp)%entity
      if (kr.eq.0) cycle
      if (region(kr)%class.eq.fbem_be) then
        do ke=1,internalelements_mesh%part(kp)%n_elements
          se=internalelements_mesh%part(kp)%element(ke)
          allocate (internalelements_mesh%element(se)%value_r(1:problem%n,internalelements_mesh%element(se)%n_nodes,0:problem%n))
        end do
      end if
    end do
  end if

  if (verbose_level.ge.1)  call fbem_timestamp_w_message(output_unit,2,'END building auxiliary variables')

end subroutine build_auxiliary_variables_mechanics_static
