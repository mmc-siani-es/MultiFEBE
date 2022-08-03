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

subroutine build_auxiliary_variables_laplace

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
  integer                    :: kr, kb, ke, kn
  integer                    :: sb, se, sn
  integer                    :: row, col
  integer                    :: k
  logical, allocatable       :: node_used(:)
  character(len=fbem_fmtstr) :: fmtstr
  integer(kind=int64)        :: memory

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'START building auxiliary variables')

  ! Allocate auxiliary variable
  allocate (node_used(n_nodes))

  ! Initialize LSE row (equations) and columns (variables) counters
  row=1
  col=1
  ! Initialize the node usage indicator
  node_used=.false.
  !
  ! Loop through BE REGIONS
  !
  do kr=1,n_regions
    if (region(kr)%class.eq.fbem_be) then
      do kb=1,region(kr)%n_boundaries
        sb=region(kr)%boundary(kb)
        do ke=1,part(boundary(sb)%part)%n_elements
          se=part(boundary(sb)%part)%element(ke)
          do kn=1,element(se)%n_nodes
            sn=element(se)%node(kn)
            if (node_used(sn).eqv.(.false.)) then
              node_used(sn)=.true.
              select case (boundary(sb)%coupling)

                ! ==================================================================================================================
                ! BE BOUNDARY
                ! ==================================================================================================================

                case (fbem_boundary_coupling_be)
                  select case (boundary(sb)%class)

                    ! =================
                    ! ORDINARY BOUNDARY
                    ! =================

                    case (fbem_boundary_class_ordinary)
                      !
                      ! Index of equations:
                      ! node(sn)%row(1,1): SBIE, HBIE or SBIE+beta*HBIE
                      !
                      ! Index of variables:
                      ! node(sn)%col(1,1): p
                      ! node(sn)%col(2,1): j
                      !
                      allocate (node(sn)%row(1,1))
                      allocate (node(sn)%col(2,1))
                      allocate (node(sn)%value_r(2,1))
                      allocate (node(sn)%dvda_r(2,1,problem%n_designvariables))
                      ! Initialize
                      node(sn)%row(1,1)=0
                      node(sn)%col(1,1)=0
                      node(sn)%col(2,1)=0
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
                      end select
                      ! Increment counter
                      col=col+1

                    ! ===================
                    ! CRACK-LIKE BOUNDARY
                    ! ===================

                    case (fbem_boundary_class_cracklike)
                      !
                      ! Index of equations (Dual BEM formulation):
                      ! node(sn)%row(1,1): SBIE
                      ! node(sn)%row(1,2): HBIE
                      !
                      ! Index of variables:
                      ! node(sn)%col(1,1): p for face +
                      ! node(sn)%col(2,1): j for face +
                      ! node(sn)%col(1,2): p for face -
                      ! node(sn)%col(2,2): j for face -
                      !
                      allocate (node(sn)%row(2,2))
                      allocate (node(sn)%col(2,2))
                      allocate (node(sn)%value_r(2,2))
                      allocate (node(sn)%dvda_r(2,2,problem%n_designvariables))
                      ! Initialize
                      do k=1,problem%n
                        node(sn)%row(1,1)=0
                        node(sn)%row(1,2)=0
                        node(sn)%col(1,1)=0
                        node(sn)%col(2,1)=0
                        node(sn)%col(1,2)=0
                        node(sn)%col(2,2)=0
                      end do
                      ! Assign row to the equation and column to the variables
                      ! Equation for SBIE
                      node(sn)%row(1,1)=row
                      ! Equation for HBIE
                      node(sn)%row(1,2)=row+1
                      ! Increment counter
                      row=row+2
                      ! Face +
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
                      ! Face -
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

                  end select

                ! ==================================================================================================================
                ! BE-BE BOUNDARY
                ! ==================================================================================================================
                !
                ! The region 1 is the region where the boundary has N+
                ! The region 2 is the region where the boundary has N-

                case (fbem_boundary_coupling_be_be)
                  !
                  ! Index of equations for each coordinate k:
                  ! node(sn)%row(1,1): SBIE, HBIE or SBIE+beta*HBIE for region 1
                  ! node(sn)%row(1,2): SBIE, HBIE or SBIE+beta*HBIE for region 2
                  !
                  ! Index of variables for each coordinate k:
                  ! node(sn)%col(1,1): p for region 1
                  ! node(sn)%col(2,1): j for region 1
                  ! node(sn)%col(1,2): p for region 2 (not active since p_{region 1}=p_{region 2})
                  ! node(sn)%col(2,2): j for region 2 (not active since j_{region 1}=-j_{region 2})
                  !
                  allocate (node(sn)%row(1,2))
                  allocate (node(sn)%col(2,2))
                  allocate (node(sn)%value_r(2,2))
                  allocate (node(sn)%dvda_r(2,2,problem%n_designvariables))
                  ! Initialize
                  node(sn)%row(1,1)=0
                  node(sn)%row(1,2)=0
                  node(sn)%col(1,1)=0
                  node(sn)%col(2,1)=0
                  node(sn)%col(1,2)=0
                  node(sn)%col(2,2)=0
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

                ! ==================================================================================================================
                ! BE-FE BOUNDARY
                ! ==================================================================================================================

                case (fbem_boundary_coupling_be_fe)
                  stop 'not implemented yet'

                ! ==================================================================================================================
                ! BE-FE-BE BOUNDARY
                ! ==================================================================================================================
                !
                ! The region 1 is the region where the boundary has N+
                ! The region 2 is the region where the boundary has N-

                case (fbem_boundary_coupling_be_fe_be)
                  stop 'not implemented yet'

              end select
            end if
          end do
        end do
      end do
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
  ! Description of value_r:
  ! internalpoint(kip)%value_r(1,0): pressure                 : p
  ! internalpoint(kip)%value_r(1,i): flux                     : j with normal n=e_i
  do k=1,n_internalpoints
    allocate (internalpoint(k)%value_r(1,0:problem%n))
    allocate (internalpoint(k)%dvda_r(1,0:problem%n,problem%n_designvariables))
  end do

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'END building auxiliary variables')

end subroutine build_auxiliary_variables_laplace
