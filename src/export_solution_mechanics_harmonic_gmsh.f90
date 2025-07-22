! ---------------------------------------------------------------------
! Copyright (C) 2014-2024 Universidad de Las Palmas de Gran Canaria:
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

subroutine export_solution_mechanics_harmonic_gmsh(kf,output_fileunit)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_shape_functions
  use fbem_geometry
  use fbem_string_handling
  use fbem_numerical
  use fbem_fem_shells
  use fbem_gmsh

  ! Module of problem variables
  use problem_variables

  ! No implicit variables allowed
  implicit none

  ! I/O variables
  integer                                 :: kf
  integer                                 :: output_fileunit
  ! Local variables
  real(kind=real64)                       :: omega
  integer                                 :: k
  integer                                 :: kr, kb, ke, kn, kip, sip, kc, ki, kj, kk, kp
  integer                                 :: sb, se, sn, ks, ss
  logical                                 :: node_used(n_nodes)
  character(len=fbem_fmtstr)              :: fmt1, fmt2
  integer                                 :: ncint, ncreal, nc, ncmax
  character(len=fbem_string_max_length)   :: tmp_string
  integer                                 :: face, k_start, k_end

  integer                                 :: exp_n_nodes
  integer, allocatable                    :: exp_node_eid(:)
  real   (kind=real64), allocatable       :: exp_node_value_r(:,:)
  complex(kind=real64), allocatable       :: exp_node_value_c(:,:)

  integer                                 :: exp_n_elements
  integer, allocatable                    :: exp_element_eid(:)
  integer, allocatable                    :: exp_element_n_nodes(:)
  real   (kind=real64), allocatable       :: exp_element_value_r(:,:)
  complex(kind=real64), allocatable       :: exp_element_value_c(:,:)
  real   (kind=real64), allocatable       :: exp_element_node_value_r(:,:,:)
  complex(kind=real64), allocatable       :: exp_element_node_value_c(:,:,:)

  integer                                 :: kke, kkn, sse
  complex(kind=real64)                    :: dpdx_tangential(problem%n), dudx_tangential(problem%n,problem%n), dudx_tangential_dil
  complex(kind=real64)                    :: u_total(problem%n), sigma(problem%n,problem%n)
  complex(kind=real64)                    :: par_J, par_Z, par_rho2, par_rhoa, par_b, mu, lambda
  complex(kind=real64)                    :: cte
  real(kind=real64), allocatable          :: psi(:,:), xi(:)
  real(kind=real64)                       :: delta
  real(kind=real64)                       :: local_axes(problem%n,problem%n)
  real(kind=real64)                       :: ep1(3), ep2(3), ep3(3), N(3), T1(3), xi2d(2), xi1d

  !
  ! Nota: usar las siguientes funciones de la libreria fbem_gmsh
  !~   public :: fbem_export_gmsh_NodeData
  !~   public :: fbem_export_gmsh_ElementData
  !~   public :: fbem_export_gmsh_ElementNodeData
  !

  ! Frequency
  omega=frequency(kf)
  if (frequency_units.eq.'f') omega=omega*c_1_2pi

  ! ################################################################################################################################
  !
  ! SOLUTION MESH
  !
  ! ################################################################################################################################

  ! Allocate working variable
  allocate (exp_node_eid(n_nodes))
  allocate (exp_node_value_r(9,n_nodes))
  allocate (exp_node_value_c(9,n_nodes))
  allocate (exp_element_eid(n_elements))
  allocate (exp_element_n_nodes(n_elements))
  allocate (exp_element_value_r(9,n_elements))
  allocate (exp_element_value_c(9,n_elements))
  allocate (exp_element_node_value_r(9,maxval(fbem_n_nodes),n_elements))
  allocate (exp_element_node_value_c(9,maxval(fbem_n_nodes),n_elements))

  ! ============================================================================================================================
  ! Fluid pressure (positive faces)
  ! ============================================================================================================================

  exp_n_nodes=0
  exp_node_eid=0
  exp_node_value_c=0
  node_used=.false.
  ! Loop through REGIONS
  do kr=1,n_regions
    select case (region(kr)%class)

      ! ========================================================================================================================
      ! BE REGION
      !
      case (fbem_be)

        ! +----------------+
        ! | BOUNDARY NODES |
        ! +----------------+

        ! Select the start and end indices of variables
        select case (region(kr)%type)
          ! INVISCID FLUID
          case (fbem_potential)
            k_start=1
            cte = (1.d0,0.d0)
          ! VISCOELASTIC SOLID
          case (fbem_viscoelastic)
            cycle
          ! POROELASTIC MEDIA
          case (fbem_poroelastic)
            k_start=0
            cte = -(1.d0,0.d0)/region(kr)%property_r(8)
        end select
        ! Loop through boundaries
        do kb=1,region(kr)%n_boundaries
          sb=region(kr)%boundary(kb)
          ! Initialize
          node_used=.false.
          ! Loop through elements and nodes
          do ke=1,part(boundary(sb)%part)%n_elements
            se=part(boundary(sb)%part)%element(ke)
            do kn=1,element(se)%n_nodes
              sn=element(se)%node(kn)
              if (.not.node_used(sn)) then
                node_used(sn)=.true.
                select case (boundary(sb)%coupling)
                  !
                  ! Uncoupled boundary of BE-FE coupled boundary
                  !
                  case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
                    select case (boundary(sb)%class)
                      !
                      ! Ordinary boundary
                      !
                      case (fbem_boundary_class_ordinary)
                        exp_n_nodes=exp_n_nodes+1
                        exp_node_eid(exp_n_nodes)=node(sn)%id
                        exp_node_value_c(1,exp_n_nodes)=cte*node(sn)%value_c(k_start,1)
                      !
                      ! Crack-like boundaries
                      !
                      case (fbem_boundary_class_cracklike)
                        exp_n_nodes=exp_n_nodes+1
                        exp_node_eid(exp_n_nodes)=node(sn)%id
                        exp_node_value_c(1,exp_n_nodes)=cte*node(sn)%value_c(k_start,1)
                    end select
                  !
                  ! BE-BE coupled boundary of BE-FE-BE coupled boundary
                  !
                  case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                    ! The region is the region 1 of the boundary
                    if (.not.region(kr)%boundary_reversion(kb)) then
                      exp_n_nodes=exp_n_nodes+1
                      exp_node_eid(exp_n_nodes)=node(sn)%id
                      exp_node_value_c(1,exp_n_nodes)=cte*node(sn)%value_c(k_start,1)
                    end if
                end select
              end if
            end do
          end do
        end do


        ! +-----------------+
        ! | INTERNAL POINTS |
        ! +-----------------+

        !
        ! To do ...
        !

      ! ========================================================================================================================

      ! ========================================================================================================================
      ! FE region
      !
      case (fbem_fe)
        ! Nothing to do

      ! ========================================================================================================================

    end select
  end do ! Loop through REGIONS
  if (exp_n_nodes.gt.0) then
    exp_node_value_r=dreal(exp_node_value_c)
    call fbem_export_gmsh_NodeData(output_fileunit,'Fluid P (BE +n)',omega,2*kf-2,1,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
    exp_node_value_r=dimag(exp_node_value_c)
    call fbem_export_gmsh_NodeData(output_fileunit,'Fluid P (BE +n)',omega,2*kf-1,1,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
  end if

  ! ============================================================================================================================
  ! Fluid pressure (negative faces)
  ! ============================================================================================================================

  exp_n_nodes=0
  exp_node_eid=0
  exp_node_value_c=0
  node_used=.false.
  ! Loop through REGIONS
  do kr=1,n_regions
    select case (region(kr)%class)

      ! ========================================================================================================================
      ! BE REGION
      !
      case (fbem_be)

        ! +----------------+
        ! | BOUNDARY NODES |
        ! +----------------+

        ! Select the start and end indices of variables
        select case (region(kr)%type)
          ! INVISCID FLUID
          case (fbem_potential)
            k_start=1
            cte = (1.d0,0.d0)
          ! VISCOELASTIC SOLID
          case (fbem_viscoelastic)
            cycle
          ! POROELASTIC MEDIA
          case (fbem_poroelastic)
            k_start=0
            cte = -(1.d0,0.d0)/region(kr)%property_r(8)
        end select
        ! Loop through boundaries
        do kb=1,region(kr)%n_boundaries
          sb=region(kr)%boundary(kb)
          ! Initialize
          node_used=.false.
          ! Loop through elements and nodes
          do ke=1,part(boundary(sb)%part)%n_elements
            se=part(boundary(sb)%part)%element(ke)
            do kn=1,element(se)%n_nodes
              sn=element(se)%node(kn)
              if (.not.node_used(sn)) then
                node_used(sn)=.true.
                select case (boundary(sb)%coupling)
                  !
                  ! Uncoupled boundary of BE-FE coupled boundary
                  !
                  case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
                    select case (boundary(sb)%class)
                      !
                      ! Ordinary boundary
                      !
                      case (fbem_boundary_class_ordinary)

                      !
                      ! Crack-like boundaries
                      !
                      case (fbem_boundary_class_cracklike)
                        exp_n_nodes=exp_n_nodes+1
                        exp_node_eid(exp_n_nodes)=node(sn)%id
                        exp_node_value_c(1,exp_n_nodes)=cte*node(sn)%value_c(k_start,2)
                    end select
                  !
                  ! BE-BE coupled boundary of BE-FE-BE coupled boundary
                  !
                  case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                    ! The region is the region 1 of the boundary
                    if (region(kr)%boundary_reversion(kb)) then
                      exp_n_nodes=exp_n_nodes+1
                      exp_node_eid(exp_n_nodes)=node(sn)%id
                      exp_node_value_c(1,exp_n_nodes)=cte*node(sn)%value_c(k_start,2)
                    end if
                end select
              end if
            end do
          end do
        end do


        ! +-----------------+
        ! | INTERNAL POINTS |
        ! +-----------------+

        !
        ! To do ...
        !

      ! ========================================================================================================================

      ! ========================================================================================================================
      ! FE region
      !
      case (fbem_fe)
        ! Nothing to do

      ! ========================================================================================================================

    end select
  end do ! Loop through REGIONS
  if (exp_n_nodes.gt.0) then
    exp_node_value_r=dreal(exp_node_value_c)
    call fbem_export_gmsh_NodeData(output_fileunit,'Fluid P (BE -n)',omega,2*kf-2,1,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
    exp_node_value_r=dimag(exp_node_value_c)
    call fbem_export_gmsh_NodeData(output_fileunit,'Fluid P (BE -n)',omega,2*kf-1,1,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
  end if

  ! ============================================================================================================================
  ! Fluid total displacements (positive faces)
  ! ============================================================================================================================

  exp_n_elements=0
  exp_element_eid=0
  exp_element_n_nodes=0
  exp_element_node_value_c=0
  ! Loop through REGIONS
  do kr=1,n_regions
    select case (region(kr)%class)

      ! ========================================================================================================================
      ! BE REGION
      !
      case (fbem_be)

        ! +----------------+
        ! | BOUNDARY NODES |
        ! +----------------+

        ! Select the start and end indices of variables
        select case (region(kr)%type)
          ! INVISCID FLUID
          case (fbem_potential)
            k_start=1
            k_end  =problem%n
          ! VISCOELASTIC SOLID
          case (fbem_viscoelastic)
            k_start=0
            k_end  =0
          ! POROELASTIC MEDIA
          case (fbem_poroelastic)
            k_start=1
            k_end  =problem%n
        end select
        if (k_end.eq.0) cycle

        ! Loop through boundaries
        do kb=1,region(kr)%n_boundaries
          sb=region(kr)%boundary(kb)
          ! Loop through elements
          do ke=1,part(boundary(sb)%part)%n_elements
            se=part(boundary(sb)%part)%element(ke)

            face=1
            exp_n_elements=exp_n_elements+1
            exp_element_eid(exp_n_elements)=element(se)%id
            exp_element_n_nodes(exp_n_elements)=element(se)%n_nodes
            do kn=1,element(se)%n_nodes
              exp_element_node_value_c(k_start:k_end,kn,exp_n_elements)=element(se)%value_c(k_start:k_end,kn,face)
            end do

          end do
        end do

      ! ========================================================================================================================

      ! ========================================================================================================================
      ! FE region
      !
      case (fbem_fe)
        ! Nothing to do
      ! ========================================================================================================================

    end select
  end do ! Loop through REGIONS
  if (exp_n_elements.gt.0) then
    select case (problem%n)
      case (2)
        exp_element_node_value_r=dreal(exp_element_node_value_c)
        call fbem_export_gmsh_ElementNodeData(output_fileunit,'Fluid UX,UY (BE +n)',omega,2*kf-2,3,&
                                              n_elements,maxval(fbem_n_nodes),&
                                              exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value_r)
        exp_element_node_value_r=dimag(exp_element_node_value_c)
        call fbem_export_gmsh_ElementNodeData(output_fileunit,'Fluid UX,UY (BE +n)',omega,2*kf-1,3,&
                                              n_elements,maxval(fbem_n_nodes),&
                                              exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value_r)
      case (3)
        exp_element_node_value_r=dreal(exp_element_node_value_c)
        call fbem_export_gmsh_ElementNodeData(output_fileunit,'Fluid UX,UY,UZ (BE +n)',omega,2*kf-2,3,&
                                              n_elements,maxval(fbem_n_nodes),&
                                              exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value_r)
        exp_element_node_value_r=dimag(exp_element_node_value_c)
        call fbem_export_gmsh_ElementNodeData(output_fileunit,'Fluid UX,UY,UZ (BE +n)',omega,2*kf-1,3,&
                                              n_elements,maxval(fbem_n_nodes),&
                                              exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value_r)
      case default
        stop 'not valid problem%n'
    end select
  end if

  ! ============================================================================================================================
  ! Fluid total displacements (negative faces)
  ! ============================================================================================================================

  exp_n_elements=0
  exp_element_eid=0
  exp_element_n_nodes=0
  exp_element_node_value_c=0
  ! Loop through REGIONS
  do kr=1,n_regions
    select case (region(kr)%class)

      ! ========================================================================================================================
      ! BE REGION
      !
      case (fbem_be)

        ! +----------------+
        ! | BOUNDARY NODES |
        ! +----------------+

        ! Select the start and end indices of variables
        select case (region(kr)%type)
          ! INVISCID FLUID
          case (fbem_potential)
            k_start=1
            k_end  =problem%n
          ! VISCOELASTIC SOLID
          case (fbem_viscoelastic)
            k_start=0
            k_end  =0
          ! POROELASTIC MEDIA
          case (fbem_poroelastic)
            k_start=1
            k_end  =problem%n
        end select
        if (k_end.eq.0) cycle

        ! Loop through boundaries
        do kb=1,region(kr)%n_boundaries
          sb=region(kr)%boundary(kb)
          ! Loop through elements
          do ke=1,part(boundary(sb)%part)%n_elements
            se=part(boundary(sb)%part)%element(ke)

                select case (boundary(sb)%coupling)
                  !
                  ! Uncoupled boundary of BE-FE coupled boundary
                  !
                  case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
                    select case (boundary(sb)%class)
                      !
                      ! Ordinary boundary
                      !
                      case (fbem_boundary_class_ordinary)

                      !
                      ! Crack-like boundaries
                      !
                      case (fbem_boundary_class_cracklike)
                        face=2
                        exp_n_elements=exp_n_elements+1
                        exp_element_eid(exp_n_elements)=element(se)%id
                        exp_element_n_nodes(exp_n_elements)=element(se)%n_nodes
                        do kn=1,element(se)%n_nodes
                          exp_element_node_value_c(k_start:k_end,kn,exp_n_elements)=element(se)%value_c(k_start:k_end,kn,face)
                        end do

                    end select
                  !
                  ! BE-BE coupled boundary of BE-FE-BE coupled boundary
                  !
                  case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                    ! The region is the region 1 of the boundary
                    if (region(kr)%boundary_reversion(kb)) then
                      face=2
                      exp_n_elements=exp_n_elements+1
                      exp_element_eid(exp_n_elements)=element(se)%id
                      exp_element_n_nodes(exp_n_elements)=element(se)%n_nodes
                      do kn=1,element(se)%n_nodes
                        exp_element_node_value_c(k_start:k_end,kn,exp_n_elements)=element(se)%value_c(k_start:k_end,kn,face)
                      end do
                    end if
                end select

          end do
        end do

      ! ========================================================================================================================

      ! ========================================================================================================================
      ! FE region
      !
      case (fbem_fe)
        ! Nothing to do
      ! ========================================================================================================================

    end select
  end do ! Loop through REGIONS
  if (exp_n_elements.gt.0) then
    select case (problem%n)
      case (2)
        exp_element_node_value_r=dreal(exp_element_node_value_c)
        call fbem_export_gmsh_ElementNodeData(output_fileunit,'Fluid UX,UY (BE -n)',omega,2*kf-2,3,&
                                              n_elements,maxval(fbem_n_nodes),&
                                              exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value_r)
        exp_element_node_value_r=dimag(exp_element_node_value_c)
        call fbem_export_gmsh_ElementNodeData(output_fileunit,'Fluid UX,UY (BE -n)',omega,2*kf-1,3,&
                                              n_elements,maxval(fbem_n_nodes),&
                                              exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value_r)
      case (3)
        exp_element_node_value_r=dreal(exp_element_node_value_c)
        call fbem_export_gmsh_ElementNodeData(output_fileunit,'Fluid UX,UY,UZ (BE -n)',omega,2*kf-2,3,&
                                              n_elements,maxval(fbem_n_nodes),&
                                              exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value_r)
        exp_element_node_value_r=dimag(exp_element_node_value_c)
        call fbem_export_gmsh_ElementNodeData(output_fileunit,'Fluid UX,UY,UZ (BE -n)',omega,2*kf-1,3,&
                                              n_elements,maxval(fbem_n_nodes),&
                                              exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value_r)
      case default
        stop 'not valid problem%n'
    end select
  end if

  ! ============================================================================================================================
  ! Solid total displacements (positive faces)
  ! ============================================================================================================================

  exp_n_nodes=0
  exp_node_eid=0
  exp_node_value_c=0
  node_used=.false.
  ! Loop through REGIONS
  do kr=1,n_regions
    select case (region(kr)%class)

      ! ========================================================================================================================
      ! BE REGION
      !
      case (fbem_be)

        ! +----------------+
        ! | BOUNDARY NODES |
        ! +----------------+

        ! Select the start and end indices of variables
        select case (region(kr)%type)
          ! INVISCID FLUID
          case (fbem_potential)
            k_start=0
            k_end  =0
          ! VISCOELASTIC SOLID
          case (fbem_viscoelastic)
            k_start=1
            k_end  =problem%n
          ! POROELASTIC MEDIA
          case (fbem_poroelastic)
            k_start=1
            k_end  =problem%n
        end select
        if (k_end.eq.0) cycle
        ! Loop through boundaries
        do kb=1,region(kr)%n_boundaries
          sb=region(kr)%boundary(kb)
          ! Initialize
          !node_used=.false.
          ! Loop through elements and nodes
          do ke=1,part(boundary(sb)%part)%n_elements
            se=part(boundary(sb)%part)%element(ke)
            do kn=1,element(se)%n_nodes
              sn=element(se)%node(kn)
              if (.not.node_used(sn)) then
                node_used(sn)=.true.
                select case (boundary(sb)%coupling)
                  !
                  ! Uncoupled boundary of BE-FE coupled boundary
                  !
                  case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
                    select case (boundary(sb)%class)
                      !
                      ! Ordinary boundary
                      !
                      case (fbem_boundary_class_ordinary)
                        exp_n_nodes=exp_n_nodes+1
                        exp_node_eid(exp_n_nodes)=node(sn)%id
                        exp_node_value_c(1:problem%n,exp_n_nodes)=node(sn)%value_c(k_start:k_end,1)
                      !
                      ! Crack-like boundaries
                      !
                      case (fbem_boundary_class_cracklike)
                        exp_n_nodes=exp_n_nodes+1
                        exp_node_eid(exp_n_nodes)=node(sn)%id
                        exp_node_value_c(1:problem%n,exp_n_nodes)=node(sn)%value_c(k_start:k_end,1)
                    end select
                  !
                  ! BE-BE coupled boundary of BE-FE-BE coupled boundary
                  !
                  case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                    ! The region is the region 1 of the boundary
                    if (.not.region(kr)%boundary_reversion(kb)) then
                      exp_n_nodes=exp_n_nodes+1
                      exp_node_eid(exp_n_nodes)=node(sn)%id
                      exp_node_value_c(1:problem%n,exp_n_nodes)=node(sn)%value_c(k_start:k_end,1)
                    end if
                end select
              end if
            end do
          end do
        end do


        ! +-----------------+
        ! | INTERNAL POINTS |
        ! +-----------------+

        !
        ! To do ...
        !


      ! ========================================================================================================================

      ! ========================================================================================================================
      ! FE region
      !
      case (fbem_fe)
        !node_used=.false.
        do ks=1,region(kr)%n_fe_subregions
          ss=region(kr)%fe_subregion(ks)
          do ke=1,part(fe_subregion(ss)%part)%n_elements
            se=part(fe_subregion(ss)%part)%element(ke)
            do kn=1,element(se)%n_nodes
              sn=element(se)%node(kn)
              if ((.not.node_used(sn)).and.(node(sn)%n_dof.gt.0)) then
                node_used(sn)=.true.
                exp_n_nodes=exp_n_nodes+1
                exp_node_eid(exp_n_nodes)=node(sn)%id
                exp_node_value_c(1:problem%n,exp_n_nodes)=node(sn)%value_c(1:problem%n,1)
              end if
            end do
          end do
        end do
      ! ========================================================================================================================

    end select
  end do ! Loop through REGIONS
  if (exp_n_nodes.gt.0) then
    select case (problem%n)
      case (2)
        exp_node_value_r=dreal(exp_node_value_c)
        call fbem_export_gmsh_NodeData(output_fileunit,'Solid UX,UY (FE, BE +n)',omega,2*kf-2,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
        exp_node_value_r=dimag(exp_node_value_c)
        call fbem_export_gmsh_NodeData(output_fileunit,'Solid UX,UY (FE, BE +n)',omega,2*kf-1,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
      case (3)
        exp_node_value_r=dreal(exp_node_value_c)
        call fbem_export_gmsh_NodeData(output_fileunit,'Solid UX,UY,UZ (FE, BE +n)',omega,2*kf-2,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
        exp_node_value_r=dimag(exp_node_value_c)
        call fbem_export_gmsh_NodeData(output_fileunit,'Solid UX,UY,UZ (FE, BE +n)',omega,2*kf-1,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
      case default
        stop 'not valid problem%n'
    end select
  end if

  ! ============================================================================================================================
  ! Solid total displacements (negative faces)
  ! ============================================================================================================================

  exp_n_nodes=0
  exp_node_eid=0
  exp_node_value_c=0
  node_used=.false.
  ! Loop through REGIONS
  do kr=1,n_regions
    select case (region(kr)%class)

    ! ==========================================================================================================================
      ! BE REGION
      !
      case (fbem_be)

        ! +----------------+
        ! | BOUNDARY NODES |
        ! +----------------+

        ! Select the start and end indices of variables
        select case (region(kr)%type)
          ! INVISCID FLUID
          case (fbem_potential)
            k_start=0
            k_end  =0
          ! VISCOELASTIC SOLID
          case (fbem_viscoelastic)
            k_start=1
            k_end  =problem%n
          ! POROELASTIC MEDIA
          case (fbem_poroelastic)
            k_start=1
            k_end  =problem%n
        end select
        if (k_end.eq.0) cycle
        ! Loop through boundaries
        do kb=1,region(kr)%n_boundaries
          sb=region(kr)%boundary(kb)
          ! Initialize
          !node_used=.false.
          ! Loop through elements and nodes
          do ke=1,part(boundary(sb)%part)%n_elements
            se=part(boundary(sb)%part)%element(ke)
            do kn=1,element(se)%n_nodes
              sn=element(se)%node(kn)
              if (.not.node_used(sn)) then
                node_used(sn)=.true.
                select case (boundary(sb)%coupling)
                  !
                  ! Uncoupled boundary of BE-FE coupled boundary
                  !
                  case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
                    select case (boundary(sb)%class)
                      !
                      ! Ordinary boundary
                      !
                      case (fbem_boundary_class_ordinary)

                      !
                      ! Crack-like boundaries
                      !
                      case (fbem_boundary_class_cracklike)

                        exp_n_nodes=exp_n_nodes+1
                        exp_node_eid(exp_n_nodes)=node(sn)%id
                        exp_node_value_c(1:problem%n,exp_n_nodes)=node(sn)%value_c(k_start:k_end,2)

                    end select
                  !
                  ! BE-BE coupled boundary of BE-FE-BE coupled boundary
                  !
                  case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                    ! The region is the region 2 of the boundary
                    if (region(kr)%boundary_reversion(kb)) then
                      exp_n_nodes=exp_n_nodes+1
                      exp_node_eid(exp_n_nodes)=node(sn)%id
                      exp_node_value_c(1:problem%n,exp_n_nodes)=node(sn)%value_c(k_start:k_end,2)
                    end if
                end select
              end if
            end do
          end do
        end do


        ! +-----------------+
        ! | INTERNAL POINTS |
        ! +-----------------+

        !
        ! Nothing to do
        !

      ! ========================================================================================================================

      ! ========================================================================================================================
      ! FE region
      !
      case (fbem_fe)

        !
        ! Nothing to do
        !

      ! ========================================================================================================================

    end select
  end do ! Loop through REGIONS
  if (exp_n_nodes.gt.0) then
    select case (problem%n)
      case (2)
        exp_node_value_r=dreal(exp_node_value_c)
        call fbem_export_gmsh_NodeData(output_fileunit,'Solid UX,UY (BE -n)',omega,2*kf-2,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
        exp_node_value_r=dimag(exp_node_value_c)
        call fbem_export_gmsh_NodeData(output_fileunit,'Solid UX,UY (BE -n)',omega,2*kf-1,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
      case (3)
        exp_node_value_r=dreal(exp_node_value_c)
        call fbem_export_gmsh_NodeData(output_fileunit,'Solid UX,UY,UZ (BE -n)',omega,2*kf-2,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
        exp_node_value_r=dimag(exp_node_value_c)
        call fbem_export_gmsh_NodeData(output_fileunit,'Solid UX,UY,UZ (BE -n)',omega,2*kf-1,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
      case default
        stop 'not valid problem%n'
    end select
  end if

  ! ============================================================================================================================
  ! Solid total tractions (positive faces)
  ! ============================================================================================================================

  exp_n_nodes=0
  exp_node_eid=0
  exp_node_value_c=0
  ! Loop through REGIONS
  do kr=1,n_regions
    select case (region(kr)%class)

      ! ========================================================================================================================
      ! BE REGION
      !
      case (fbem_be)

        ! +----------------+
        ! | BOUNDARY NODES |
        ! +----------------+

        ! Select the start and end indices of variables
        select case (region(kr)%type)
          ! INVISCID FLUID
          case (fbem_potential)
            k_start=0
            k_end  =0
          ! VISCOELASTIC SOLID
          case (fbem_viscoelastic)
            k_start=problem%n+1
            k_end  =2*problem%n
          ! POROELASTIC MEDIA
          case (fbem_poroelastic)
            k_start=2+problem%n
            k_end  =1+2*problem%n
        end select
        if (k_end.eq.0) cycle
        ! Loop through boundaries
        do kb=1,region(kr)%n_boundaries
          sb=region(kr)%boundary(kb)
          ! Initialize
          node_used=.false.
          ! Loop through elements and nodes
          do ke=1,part(boundary(sb)%part)%n_elements
            se=part(boundary(sb)%part)%element(ke)
            do kn=1,element(se)%n_nodes
              sn=element(se)%node(kn)
              if (.not.node_used(sn)) then
                node_used(sn)=.true.
                select case (boundary(sb)%coupling)
                  !
                  ! Uncoupled boundary of BE-FE coupled boundary
                  !
                  case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
                    select case (boundary(sb)%class)
                      !
                      ! Ordinary boundary
                      !
                      case (fbem_boundary_class_ordinary)
                        exp_n_nodes=exp_n_nodes+1
                        exp_node_eid(exp_n_nodes)=node(sn)%id
                        exp_node_value_c(1:problem%n,exp_n_nodes)=node(sn)%value_c(k_start:k_end,1)
                      !
                      ! Crack-like boundaries
                      !
                      case (fbem_boundary_class_cracklike)
                        exp_n_nodes=exp_n_nodes+1
                        exp_node_eid(exp_n_nodes)=node(sn)%id
                        exp_node_value_c(1:problem%n,exp_n_nodes)=node(sn)%value_c(k_start:k_end,1)
                    end select
                  !
                  ! BE-BE coupled boundary of BE-FE-BE coupled boundary
                  !
                  case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                    ! The region is the region 1 of the boundary
                    if (.not.region(kr)%boundary_reversion(kb)) then
                      exp_n_nodes=exp_n_nodes+1
                      exp_node_eid(exp_n_nodes)=node(sn)%id
                      exp_node_value_c(1:problem%n,exp_n_nodes)=node(sn)%value_c(k_start:k_end,1)
                    end if
                end select
              end if
            end do
          end do
        end do


        ! +-----------------+
        ! | INTERNAL POINTS |
        ! +-----------------+

        !
        ! To do ...
        !


      ! ========================================================================================================================

      ! ========================================================================================================================
      ! FE region
      !
      case (fbem_fe)

        !
        ! Nothing to do...
        !

      ! ========================================================================================================================

    end select
  end do ! Loop through REGIONS
  if (exp_n_nodes.gt.0) then
    select case (problem%n)
      case (2)
        exp_node_value_r=dreal(exp_node_value_c)
        call fbem_export_gmsh_NodeData(output_fileunit,'Solid TX,TY (BE +n)',omega,2*kf-2,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
        exp_node_value_r=dimag(exp_node_value_c)
        call fbem_export_gmsh_NodeData(output_fileunit,'Solid TX,TY (BE +n)',omega,2*kf-1,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
      case (3)
        exp_node_value_r=dreal(exp_node_value_c)
        call fbem_export_gmsh_NodeData(output_fileunit,'Solid TX,TY,TZ (BE +n)',omega,2*kf-2,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
        exp_node_value_r=dimag(exp_node_value_c)
        call fbem_export_gmsh_NodeData(output_fileunit,'Solid TX,TY,TZ (BE +n)',omega,2*kf-1,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
      case default
        stop 'not valid problem%n'
    end select
  end if

  ! ============================================================================================================================
  ! Solid total tractions (negative faces)
  ! ============================================================================================================================

  exp_n_nodes=0
  exp_node_eid=0
  exp_node_value_c=0
  ! Loop through REGIONS
  do kr=1,n_regions
    select case (region(kr)%class)

      ! ========================================================================================================================
      ! BE REGION
      !
      case (fbem_be)

        ! +----------------+
        ! | BOUNDARY NODES |
        ! +----------------+

        ! Select the start and end indices of variables
        select case (region(kr)%type)
          ! INVISCID FLUID
          case (fbem_potential)
            k_start=0
            k_end  =0
          ! VISCOELASTIC SOLID
          case (fbem_viscoelastic)
            k_start=problem%n+1
            k_end  =2*problem%n
          ! POROELASTIC MEDIA
          case (fbem_poroelastic)
            k_start=2+problem%n
            k_end  =1+2*problem%n
        end select
        if (k_end.eq.0) cycle
        ! Loop through boundaries
        do kb=1,region(kr)%n_boundaries
          sb=region(kr)%boundary(kb)
          ! Initialize
          node_used=.false.
          ! Loop through elements and nodes
          do ke=1,part(boundary(sb)%part)%n_elements
            se=part(boundary(sb)%part)%element(ke)
            do kn=1,element(se)%n_nodes
              sn=element(se)%node(kn)
              if (.not.node_used(sn)) then
                node_used(sn)=.true.
                select case (boundary(sb)%coupling)
                  !
                  ! Uncoupled boundary of BE-FE coupled boundary
                  !
                  case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
                    select case (boundary(sb)%class)
                      !
                      ! Ordinary boundary
                      !
                      case (fbem_boundary_class_ordinary)

                        !
                        ! Nothing to do
                        !

                      !
                      ! Crack-like boundaries
                      !
                      case (fbem_boundary_class_cracklike)

                        exp_n_nodes=exp_n_nodes+1
                        exp_node_eid(exp_n_nodes)=node(sn)%id
                        exp_node_value_c(1:problem%n,exp_n_nodes)=node(sn)%value_c(k_start:k_end,2)

                    end select
                  !
                  ! BE-BE coupled boundary of BE-FE-BE coupled boundary
                  !
                  case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                    ! The region is the region 2 of the boundary
                    if (region(kr)%boundary_reversion(kb)) then
                      exp_n_nodes=exp_n_nodes+1
                      exp_node_eid(exp_n_nodes)=node(sn)%id
                      exp_node_value_c(1:problem%n,exp_n_nodes)=node(sn)%value_c(k_start:k_end,2)
                    end if
                end select
              end if
            end do
          end do
        end do


        ! +-----------------+
        ! | INTERNAL POINTS |
        ! +-----------------+

        !
        ! Nothing to do
        !

      ! ========================================================================================================================

      ! ========================================================================================================================
      ! FE region
      !
      case (fbem_fe)

        !
        ! Nothing to do
        !

      ! ========================================================================================================================

    end select
  end do ! Loop through REGIONS
  if (exp_n_nodes.gt.0) then
    select case (problem%n)
      case (2)
        exp_node_value_r=dreal(exp_node_value_c)
        call fbem_export_gmsh_NodeData(output_fileunit,'Solid TX,TY (BE -n)',omega,2*kf-2,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
        exp_node_value_r=dimag(exp_node_value_c)
        call fbem_export_gmsh_NodeData(output_fileunit,'Solid TX,TY (BE -n)',omega,2*kf-1,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
      case (3)
        exp_node_value_r=dreal(exp_node_value_c)
        call fbem_export_gmsh_NodeData(output_fileunit,'Solid TX,TY,TZ (BE -n)',omega,2*kf-2,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
        exp_node_value_r=dimag(exp_node_value_c)
        call fbem_export_gmsh_NodeData(output_fileunit,'Solid TX,TY,TZ (BE -n)',omega,2*kf-1,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
      case default
        stop 'not valid problem%n'
    end select
  end if

  ! ================================================================================================================================
  ! Stress tensor (positive faces)
  ! ================================================================================================================================

  exp_n_elements=0
  exp_element_eid=0
  exp_element_n_nodes=0
  exp_element_node_value_c=0
  ! Loop through REGIONS
  do kr=1,n_regions
    select case (region(kr)%class)

      ! ============================================================================================================================
      ! BE REGION
      !
      case (fbem_be)

        ! +-------------------+
        ! | BOUNDARY ELEMENTS |
        ! +-------------------+

        ! Select the start and end indices of variables
        select case (region(kr)%type)
          ! INVISCID FLUID
          case (fbem_potential)
            k_start=0
            k_end  =0
          ! VISCOELASTIC SOLID
          case (fbem_viscoelastic)
            k_start=1
            k_end  =9
          ! POROELASTIC MEDIA
          case (fbem_poroelastic)
            k_start=problem%n+1
            k_end  =problem%n+9
        end select
        if (k_end.eq.0) cycle
        do kb=1,region(kr)%n_boundaries
          sb=region(kr)%boundary(kb)
          do ke=1,part(boundary(sb)%part)%n_elements
            se=part(boundary(sb)%part)%element(ke)
            select case (boundary(sb)%coupling)

              case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)

                select case (boundary(sb)%class)

                  case (fbem_boundary_class_ordinary)

                    exp_n_elements=exp_n_elements+1
                    exp_element_eid(exp_n_elements)=element(se)%id
                    exp_element_n_nodes(exp_n_elements)=element(se)%n_nodes
                    face=1
                    do kn=1,element(se)%n_nodes
                      exp_element_node_value_c(k_start:k_end,kn,exp_n_elements)=element(se)%value_c(k_start:k_end,kn,face)
                    end do

                  case (fbem_boundary_class_cracklike)

                    exp_n_elements=exp_n_elements+1
                    exp_element_eid(exp_n_elements)=element(se)%id
                    exp_element_n_nodes(exp_n_elements)=element(se)%n_nodes
                    face=1
                    do kn=1,element(se)%n_nodes
                      exp_element_node_value_c(k_start:k_end,kn,exp_n_elements)=element(se)%value_c(k_start:k_end,kn,face)
                    end do

                end select
              case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                if (.not.region(kr)%boundary_reversion(kb)) then

                    exp_n_elements=exp_n_elements+1
                    exp_element_eid(exp_n_elements)=element(se)%id
                    exp_element_n_nodes(exp_n_elements)=element(se)%n_nodes
                    face=1
                    do kn=1,element(se)%n_nodes
                      exp_element_node_value_c(k_start:k_end,kn,exp_n_elements)=element(se)%value_c(k_start:k_end,kn,face)
                    end do

                end if
            end select
          end do
        end do

        ! +-----------------+
        ! | INTERNAL POINTS |
        ! +-----------------+

        ! to be done...

      ! ============================================================================================================================

      ! ============================================================================================================================
      ! FE region
      !
      case (fbem_fe)

        !
        ! Nothing to do
        !

      ! ============================================================================================================================

    end select
  end do ! Loop through REGIONS
  if (exp_n_elements.gt.0) then
    exp_element_node_value_r=dreal(exp_element_node_value_c)
    call fbem_export_gmsh_ElementNodeData(output_fileunit,'Solid Sij (BE +n)',omega,2*kf-2,9,&
                                      n_elements,maxval(fbem_n_nodes),&
                                      exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value_r)
    exp_element_node_value_r=dimag(exp_element_node_value_c)
    call fbem_export_gmsh_ElementNodeData(output_fileunit,'Solid Sij (BE +n)',omega,2*kf-1,9,&
                                      n_elements,maxval(fbem_n_nodes),&
                                      exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value_r)
  end if


  ! ================================================================================================================================
  ! Stress tensor (negative faces)
  ! ================================================================================================================================

  exp_n_elements=0
  exp_element_eid=0
  exp_element_n_nodes=0
  exp_element_node_value_c=0
  ! Loop through REGIONS
  do kr=1,n_regions
    select case (region(kr)%class)

      ! ============================================================================================================================
      ! BE REGION
      !
      case (fbem_be)

        ! +-------------------+
        ! | BOUNDARY ELEMENTS |
        ! +-------------------+

        ! Select the start and end indices of variables
        select case (region(kr)%type)
          ! INVISCID FLUID
          case (fbem_potential)
            k_start=0
            k_end  =0
          ! VISCOELASTIC SOLID
          case (fbem_viscoelastic)
            k_start=1
            k_end  =9
          ! POROELASTIC MEDIA
          case (fbem_poroelastic)
            k_start=problem%n+1
            k_end  =problem%n+9
        end select
        if (k_end.eq.0) cycle
        do kb=1,region(kr)%n_boundaries
          sb=region(kr)%boundary(kb)
          do ke=1,part(boundary(sb)%part)%n_elements
            se=part(boundary(sb)%part)%element(ke)
            select case (boundary(sb)%coupling)

              case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)

                select case (boundary(sb)%class)

                  case (fbem_boundary_class_ordinary)

                    !
                    ! Nothing to do
                    !

                  case (fbem_boundary_class_cracklike)

                    exp_n_elements=exp_n_elements+1
                    exp_element_eid(exp_n_elements)=element(se)%id
                    exp_element_n_nodes(exp_n_elements)=element(se)%n_nodes
                    face=2
                    do kn=1,element(se)%n_nodes
                      exp_element_node_value_c(k_start:k_end,kn,exp_n_elements)=element(se)%value_c(k_start:k_end,kn,face)
                    end do

                end select
              case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                if (region(kr)%boundary_reversion(kb)) then

                    exp_n_elements=exp_n_elements+1
                    exp_element_eid(exp_n_elements)=element(se)%id
                    exp_element_n_nodes(exp_n_elements)=element(se)%n_nodes
                    face=2
                    do kn=1,element(se)%n_nodes
                      exp_element_node_value_c(k_start:k_end,kn,exp_n_elements)=element(se)%value_c(k_start:k_end,kn,face)
                    end do

                end if
            end select
          end do
        end do

        ! +-----------------+
        ! | INTERNAL POINTS |
        ! +-----------------+

        ! to be done...

      ! ============================================================================================================================

      ! ============================================================================================================================
      ! FE region
      !
      case (fbem_fe)

        !
        ! Nothing to do
        !

      ! ============================================================================================================================

    end select
  end do ! Loop through REGIONS
  if (exp_n_elements.gt.0) then
    exp_element_node_value_r=dreal(exp_element_node_value_c)
    call fbem_export_gmsh_ElementNodeData(output_fileunit,'Solid Sij (BE -n)',omega,2*kf-2,9,&
                                      n_elements,maxval(fbem_n_nodes),&
                                      exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value_r)
    exp_element_node_value_r=dimag(exp_element_node_value_c)
    call fbem_export_gmsh_ElementNodeData(output_fileunit,'Solid Sij (BE -n)',omega,2*kf-1,9,&
                                      n_elements,maxval(fbem_n_nodes),&
                                      exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value_r)
  end if

  ! ================================================================================================================================
  ! FE nodal forces/reactions
  ! ================================================================================================================================

  exp_n_nodes=0
  exp_node_eid=0
  exp_node_value_c=0
  node_used=.false.
  ! Loop through REGIONS
  do kr=1,n_regions
    select case (region(kr)%class)

      ! ============================================================================================================================
      ! BE REGION
      !
      case (fbem_be)

        !
        ! Nothing to do...
        !

      ! ============================================================================================================================

      ! ============================================================================================================================
      ! FE region
      !
      case (fbem_fe)

        do ks=1,region(kr)%n_fe_subregions
          ss=region(kr)%fe_subregion(ks)
          do ke=1,part(fe_subregion(ss)%part)%n_elements
            se=part(fe_subregion(ss)%part)%element(ke)
            do kn=1,element(se)%n_nodes
              sn=element(se)%node(kn)
              if ((.not.node_used(sn)).and.(node(sn)%n_dof.gt.0)) then
                node_used(sn)=.true.
                if (node(sn)%n_dof.gt.problem%n) then
                  exp_n_nodes=exp_n_nodes+1
                  exp_node_eid(exp_n_nodes)=node(sn)%id
                  exp_node_value_c(1:problem%n,exp_n_nodes)=node(sn)%value_c(1:problem%n,2)
                end if
              end if
            end do
          end do
        end do

      ! ============================================================================================================================

    end select
  end do ! Loop through REGIONS
  if (exp_n_nodes.gt.0) then
    select case (problem%n)
      case (2)
        exp_node_value_r=dreal(exp_node_value_c)
        call fbem_export_gmsh_NodeData(output_fileunit,'Nodal forces FX,FY (FE)',omega,2*kf-2,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
        exp_node_value_r=dimag(exp_node_value_c)
        call fbem_export_gmsh_NodeData(output_fileunit,'Nodal forces FX,FY (FE)',omega,2*kf-1,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
      case (3)
        exp_node_value_r=dreal(exp_node_value_c)
        call fbem_export_gmsh_NodeData(output_fileunit,'Nodal forces FX,FY,FZ (FE)',omega,2*kf-2,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
        exp_node_value_r=dimag(exp_node_value_c)
        call fbem_export_gmsh_NodeData(output_fileunit,'Nodal forces FX,FY,FZ (FE)',omega,2*kf-1,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
      case default
        stop 'not valid problem%n'
    end select
  end if

!~   ! ================================================================================================================================
!~   ! FE nodal moment/reactions
!~   ! ================================================================================================================================

!~   exp_n_nodes=0
!~   exp_node_eid=0
!~   exp_node_value_c=0
!~   node_used=.false.
!~   ! Loop through REGIONS
!~   do kr=1,n_regions
!~     select case (region(kr)%class)

!~       ! ============================================================================================================================
!~       ! BE REGION
!~       !
!~       case (fbem_be)

!~         !
!~         ! Nothing to do...
!~         !

!~       ! ============================================================================================================================

!~       ! ============================================================================================================================
!~       ! FE region
!~       !
!~       case (fbem_fe)

!~         do ks=1,region(kr)%n_fe_subregions
!~           ss=region(kr)%fe_subregion(ks)
!~           do ke=1,part(fe_subregion(ss)%part)%n_elements
!~             se=part(fe_subregion(ss)%part)%element(ke)
!~             do kn=1,element(se)%n_nodes
!~               sn=element(se)%node(kn)
!~               if ((.not.node_used(sn)).and.(node(sn)%n_dof.gt.0)) then
!~                 node_used(sn)=.true.
!~                 if (node(sn)%n_dof.gt.problem%n) then
!~                   exp_n_nodes=exp_n_nodes+1
!~                   exp_node_eid(exp_n_nodes)=node(sn)%id
!~                   write(*,*)node(sn)%value_c(:,2)
!~                   exp_node_value_c(1:problem%n,exp_n_nodes)=node(sn)%value_c((problem%n+1):(3*(problem%n-1)),2)
!~                 end if
!~               end if
!~             end do
!~           end do
!~         end do

!~       ! ============================================================================================================================

!~     end select
!~   end do ! Loop through REGIONS
!~   if (exp_n_nodes.gt.0) then
!~     select case (problem%n)
!~       case (2)
!~         exp_node_value_r=dreal(exp_node_value_c)
!~         call fbem_export_gmsh_NodeData(output_fileunit,'Nodal moment MZ (FE)',omega,2*kf-2,1,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
!~         exp_node_value_r=dimag(exp_node_value_c)
!~         call fbem_export_gmsh_NodeData(output_fileunit,'Nodal moment MZ (FE)',omega,2*kf-1,1,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
!~       case (3)
!~         exp_node_value_r=dreal(exp_node_value_c)
!~         call fbem_export_gmsh_NodeData(output_fileunit,'Nodal moments MX,MY,MZ (FE)',omega,2*kf-2,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
!~         exp_node_value_r=dimag(exp_node_value_c)
!~         call fbem_export_gmsh_NodeData(output_fileunit,'Nodal moments MX,MY,MZ (FE)',omega,2*kf-1,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value_r)
!~       case default
!~         stop 'not valid problem%n'
!~     end select
!~   end if

  !
  ! Falta solo las cargas equilibrantes (no creo que sean muy necesarias de exportar)
  !

  ! ================================================================================================================================
  ! Beam stress resultants (as a tensor with 9 components: NX,VY,BZ,0,0,0,0,0,0 or NX,VY,VZ,BX,BY,BZ,0,0,0)
  ! ================================================================================================================================

  exp_n_elements=0
  exp_element_eid=0
  exp_element_n_nodes=0
  exp_element_node_value_c=0
  do kr=1,n_regions
    if (region(kr)%class.eq.fbem_fe) then
      do ks=1,region(kr)%n_fe_subregions
        ss=region(kr)%fe_subregion(ks)
        do ke=1,part(fe_subregion(ss)%part)%n_elements
          se=part(fe_subregion(ss)%part)%element(ke)
            if ((element(se)%n_dimension.eq.1).and.((element(se)%fe_type.eq.0).or.(element(se)%fe_type.eq.1).or.(element(se)%fe_type.eq.2))) then
              exp_n_elements=exp_n_elements+1
              exp_element_eid(exp_n_elements)=element(se)%id
              exp_element_n_nodes(exp_n_elements)=element(se)%n_nodes
              do kn=1,element(se)%n_nodes
                exp_element_node_value_c(1:3*(problem%n-1),kn,exp_n_elements)=element(se)%value_c(1:3*(problem%n-1),kn,1)
              end do
            end if
        end do
      end do
    end if
  end do
  if (exp_n_elements.gt.1) then
    select case (problem%n)
      case (2)
        exp_element_node_value_r=dreal(exp_element_node_value_c)
        call fbem_export_gmsh_ElementNodeData(output_fileunit,'Beam NX,VY,MZ (FE)',omega,2*kf-2,9,&
                                              n_elements,maxval(fbem_n_nodes),&
                                              exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value_r)
        exp_element_node_value_r=dimag(exp_element_node_value_c)
        call fbem_export_gmsh_ElementNodeData(output_fileunit,'Beam NX,VY,MZ (FE)',omega,2*kf-1,9,&
                                              n_elements,maxval(fbem_n_nodes),&
                                              exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value_r)
      case (3)
        exp_element_node_value_r=dreal(exp_element_node_value_c)
        call fbem_export_gmsh_ElementNodeData(output_fileunit,'Beam NX,VY,VZ,MX,MY,MZ (FE)',omega,2*kf-2,9,&
                                              n_elements,maxval(fbem_n_nodes),&
                                              exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value_r)
        exp_element_node_value_r=dimag(exp_element_node_value_c)
        call fbem_export_gmsh_ElementNodeData(output_fileunit,'Beam NX,VY,VZ,MX,MY,MZ (FE)',omega,2*kf-1,9,&
                                              n_elements,maxval(fbem_n_nodes),&
                                              exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value_r)
      case default
        stop 'not valid problem%n'
    end select

  end if


  if (kf.eq.1) then

    ! ================================================================================================================================
    ! Beam stress resultants local axis x
    ! ================================================================================================================================

    exp_n_elements=0
    exp_element_eid=0
    exp_element_value_r=0
    do kr=1,n_regions
      if (region(kr)%class.eq.fbem_fe) then
        do ks=1,region(kr)%n_fe_subregions
          ss=region(kr)%fe_subregion(ks)
          do ke=1,part(fe_subregion(ss)%part)%n_elements
            se=part(fe_subregion(ss)%part)%element(ke)
              if ((element(se)%n_dimension.eq.1).and.(element(se)%fe_type.eq.0)) then
                exp_n_elements=exp_n_elements+1
                exp_element_eid(exp_n_elements)=element(se)%id
                exp_element_value_r(1:3,exp_n_elements)=element(se)%v_midnode(:,1,1)
              end if
              if ((element(se)%n_dimension.eq.1).and.((element(se)%fe_type.eq.1).or.(element(se)%fe_type.eq.2))) then
                exp_n_elements=exp_n_elements+1
                exp_element_eid(exp_n_elements)=element(se)%id
                exp_element_value_r(1:3,exp_n_elements)=element(se)%ep(:,1)
              end if
          end do
        end do
      end if
    end do
    if (exp_n_elements.gt.1) then
      call fbem_export_gmsh_ElementData(output_fileunit,'Beam X axis',0.d0,0,3,&
                                        n_elements,exp_n_elements,exp_element_eid,exp_element_value_r)
    end if

    ! ================================================================================================================================
    ! Beam stress resultants local axis y
    ! ================================================================================================================================

    exp_n_elements=0
    exp_element_eid=0
    exp_element_value_r=0
    do kr=1,n_regions
      if (region(kr)%class.eq.fbem_fe) then
        do ks=1,region(kr)%n_fe_subregions
          ss=region(kr)%fe_subregion(ks)
          do ke=1,part(fe_subregion(ss)%part)%n_elements
            se=part(fe_subregion(ss)%part)%element(ke)
              if ((element(se)%n_dimension.eq.1).and.(element(se)%fe_type.eq.0)) then
                exp_n_elements=exp_n_elements+1
                exp_element_eid(exp_n_elements)=element(se)%id
                exp_element_value_r(1:3,exp_n_elements)=element(se)%v_midnode(:,2,1)
              end if
              if ((element(se)%n_dimension.eq.1).and.((element(se)%fe_type.eq.1).or.(element(se)%fe_type.eq.2))) then
                exp_n_elements=exp_n_elements+1
                exp_element_eid(exp_n_elements)=element(se)%id
                exp_element_value_r(1:3,exp_n_elements)=element(se)%ep(:,2)
              end if
          end do
        end do
      end if
    end do
    if (exp_n_elements.gt.1) then
      call fbem_export_gmsh_ElementData(output_fileunit,'Beam Y axis',0.d0,0,3,&
                                        n_elements,exp_n_elements,exp_element_eid,exp_element_value_r)
    end if

    ! ================================================================================================================================
    ! Beam stress resultants local axis z
    ! ================================================================================================================================

    if (problem%n.eq.3) then
      exp_n_elements=0
      exp_element_eid=0
      exp_element_value_r=0
      do kr=1,n_regions
        if (region(kr)%class.eq.fbem_fe) then
          do ks=1,region(kr)%n_fe_subregions
            ss=region(kr)%fe_subregion(ks)
            do ke=1,part(fe_subregion(ss)%part)%n_elements
              se=part(fe_subregion(ss)%part)%element(ke)
                if ((element(se)%n_dimension.eq.1).and.(element(se)%fe_type.eq.0)) then
                  exp_n_elements=exp_n_elements+1
                  exp_element_eid(exp_n_elements)=element(se)%id
                  exp_element_value_r(1:3,exp_n_elements)=element(se)%v_midnode(:,3,1)
                end if
                if ((element(se)%n_dimension.eq.1).and.((element(se)%fe_type.eq.1).or.(element(se)%fe_type.eq.2))) then
                  exp_n_elements=exp_n_elements+1
                  exp_element_eid(exp_n_elements)=element(se)%id
                  exp_element_value_r(1:3,exp_n_elements)=element(se)%ep(:,3)
                end if
            end do
          end do
        end if
      end do
      if (exp_n_elements.gt.1) then
        call fbem_export_gmsh_ElementData(output_fileunit,'Beam Z axis',0.d0,0,3,&
                                          n_elements,exp_n_elements,exp_element_eid,exp_element_value_r)
      end if
    end if

  end if

  ! ================================================================================================================================
  ! Bar stress resultants (as scalar)
  ! ================================================================================================================================

  exp_n_elements=0
  exp_element_eid=0
  exp_element_n_nodes=0
  exp_element_node_value_c=0
  do kr=1,n_regions
    if (region(kr)%class.eq.fbem_fe) then
      do ks=1,region(kr)%n_fe_subregions
        ss=region(kr)%fe_subregion(ks)
        do ke=1,part(fe_subregion(ss)%part)%n_elements
          se=part(fe_subregion(ss)%part)%element(ke)
            if ((element(se)%n_dimension.eq.1).and.(element(se)%fe_type.eq.3)) then
              exp_n_elements=exp_n_elements+1
              exp_element_eid(exp_n_elements)=element(se)%id
              exp_element_n_nodes(exp_n_elements)=element(se)%n_nodes
              do kn=1,element(se)%n_nodes
                exp_element_node_value_c(1,kn,exp_n_elements)=element(se)%value_c(1,kn,1)
              end do
            end if
        end do
      end do
    end if
  end do
  if (exp_n_elements.gt.1) then
    exp_element_node_value_r=dreal(exp_element_node_value_c)
    call fbem_export_gmsh_ElementNodeData(output_fileunit,'Bar NX',omega,2*kf-2,1,&
                                          n_elements,maxval(fbem_n_nodes),&
                                          exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value_r)
    exp_element_node_value_r=dimag(exp_element_node_value_c)
    call fbem_export_gmsh_ElementNodeData(output_fileunit,'Bar NX',omega,2*kf-1,1,&
                                          n_elements,maxval(fbem_n_nodes),&
                                          exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value_r)
  end if

  ! ================================================================================================================================
  ! Spring-dashpot stress resultants (as scalar)
  ! ================================================================================================================================

  exp_n_elements=0
  exp_element_eid=0
  exp_element_n_nodes=0
  exp_element_node_value_c=0
  do kr=1,n_regions
    if (region(kr)%class.eq.fbem_fe) then
      do ks=1,region(kr)%n_fe_subregions
        ss=region(kr)%fe_subregion(ks)
        do ke=1,part(fe_subregion(ss)%part)%n_elements
          se=part(fe_subregion(ss)%part)%element(ke)
            if ((element(se)%n_dimension.eq.1).and.(element(se)%fe_type.eq.6)) then
              exp_n_elements=exp_n_elements+1
              exp_element_eid(exp_n_elements)=element(se)%id
              exp_element_n_nodes(exp_n_elements)=element(se)%n_nodes
              do kn=1,element(se)%n_nodes
                exp_element_node_value_c(1,kn,exp_n_elements)=element(se)%value_c(1,kn,1)
              end do
            end if
        end do
      end do
    end if
  end do
  if (exp_n_elements.gt.1) then
    exp_element_node_value_r=dreal(exp_element_node_value_c)
    call fbem_export_gmsh_ElementNodeData(output_fileunit,'Spring/dashpot NX',omega,2*kf-2,1,&
                                          n_elements,maxval(fbem_n_nodes),&
                                          exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value_r)
    exp_element_node_value_r=dimag(exp_element_node_value_c)
    call fbem_export_gmsh_ElementNodeData(output_fileunit,'Spring/dashpot NX',omega,2*kf-1,1,&
                                          n_elements,maxval(fbem_n_nodes),&
                                          exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value_r)
  end if

  !
  ! Falta exportar esfuerzos y ejes en distra y disrotra
  !

  ! ================================================================================================================================
  ! Shell stress resultants (as a tensor with 9 components: Nx,Ny,Nxy,Mx,My,Mxy,Vx,Vy,0)
  ! ================================================================================================================================

  exp_n_elements=0
  exp_element_eid=0
  exp_element_n_nodes=0
  exp_element_node_value_c=0
  do kr=1,n_regions
    if (region(kr)%class.eq.fbem_fe) then
      do ks=1,region(kr)%n_fe_subregions
        ss=region(kr)%fe_subregion(ks)
        do ke=1,part(fe_subregion(ss)%part)%n_elements
          se=part(fe_subregion(ss)%part)%element(ke)
            if ((element(se)%n_dimension.eq.2).and.(problem%n.eq.3)) then
              exp_n_elements=exp_n_elements+1
              exp_element_eid(exp_n_elements)=element(se)%id
              exp_element_n_nodes(exp_n_elements)=element(se)%n_nodes
              do kn=1,element(se)%n_nodes
                exp_element_node_value_c(1:8,kn,exp_n_elements)=element(se)%value_c(1:8,kn,1)
              end do
            end if
        end do
      end do
    end if
  end do
  if (exp_n_elements.gt.1) then
    exp_element_node_value_r=dreal(exp_element_node_value_c)
    call fbem_export_gmsh_ElementNodeData(output_fileunit,'Shell Nx,Ny,Nxy,Mx,My,Mxy,Vx,Vy',omega,2*kf-2,9,&
                                          n_elements,maxval(fbem_n_nodes),&
                                          exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value_r)
    exp_element_node_value_r=dimag(exp_element_node_value_c)
    call fbem_export_gmsh_ElementNodeData(output_fileunit,'Shell Nx,Ny,Nxy,Mx,My,Mxy,Vx,Vy',omega,2*kf-1,9,&
                                          n_elements,maxval(fbem_n_nodes),&
                                          exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value_r)
  end if

  if (kf.eq.1) then

    ! ================================================================================================================================
    ! Shell stress resultants local axis z
    ! ================================================================================================================================

    exp_n_elements=0
    exp_element_eid=0
    exp_element_value_r=0
    do kr=1,n_regions
      if (region(kr)%class.eq.fbem_fe) then
        do ks=1,region(kr)%n_fe_subregions
          ss=region(kr)%fe_subregion(ks)
          do ke=1,part(fe_subregion(ss)%part)%n_elements
            se=part(fe_subregion(ss)%part)%element(ke)
              if ((element(se)%n_dimension.eq.2).and.(problem%n.eq.3)) then
                exp_n_elements=exp_n_elements+1
                exp_element_eid(exp_n_elements)=element(se)%id
                ! Local coordinate of the centroid in element reference coordinates
                select case (fbem_n_edges(element(se)%type))
                  case (3)
                    xi2d=1.d0/3.d0
                  case (4)
                    xi2d=0
                  case default
                    stop 'not valid n_edges'
                end select
                N=fbem_normal3d(element(se)%type,element(se)%x_gn,xi2d)
                T1=fbem_tangent_xi1(3,element(se)%type,element(se)%x_gn,xi2d)
                call fbem_fem_degshell_stress_resultants_ep(N,T1,element(se)%ep(:,1),ep1,ep2,ep3)
                exp_element_value_r(1:3,exp_n_elements)=ep1
              end if
          end do
        end do
      end if
    end do
    if (exp_n_elements.gt.1) then
      call fbem_export_gmsh_ElementData(output_fileunit,'Shell X axis',0.d0,0,3,&
                                        n_elements,exp_n_elements,exp_element_eid,exp_element_value_r)
    end if

    ! ================================================================================================================================
    ! Shell stress resultants local axis z
    ! ================================================================================================================================

    exp_n_elements=0
    exp_element_eid=0
    exp_element_value_r=0
    do kr=1,n_regions
      if (region(kr)%class.eq.fbem_fe) then
        do ks=1,region(kr)%n_fe_subregions
          ss=region(kr)%fe_subregion(ks)
          do ke=1,part(fe_subregion(ss)%part)%n_elements
            se=part(fe_subregion(ss)%part)%element(ke)
              if ((element(se)%n_dimension.eq.2).and.(problem%n.eq.3)) then
                exp_n_elements=exp_n_elements+1
                exp_element_eid(exp_n_elements)=element(se)%id
                ! Local coordinate of the centroid in element reference coordinates
                select case (fbem_n_edges(element(se)%type))
                  case (3)
                    xi2d=1.d0/3.d0
                  case (4)
                    xi2d=0
                  case default
                    stop 'not valid n_edges'
                end select
                N=fbem_normal3d(element(se)%type,element(se)%x_gn,xi2d)
                T1=fbem_tangent_xi1(3,element(se)%type,element(se)%x_gn,xi2d)
                call fbem_fem_degshell_stress_resultants_ep(N,T1,element(se)%ep(:,1),ep1,ep2,ep3)
                exp_element_value_r(1:3,exp_n_elements)=ep2
              end if
          end do
        end do
      end if
    end do
    if (exp_n_elements.gt.1) then
      call fbem_export_gmsh_ElementData(output_fileunit,'Shell Y axis',0.d0,0,3,&
                                        n_elements,exp_n_elements,exp_element_eid,exp_element_value_r)
    end if

    ! ================================================================================================================================
    ! Shell stress resultants local axis z
    ! ================================================================================================================================

    exp_n_elements=0
    exp_element_eid=0
    exp_element_value_r=0
    do kr=1,n_regions
      if (region(kr)%class.eq.fbem_fe) then
        do ks=1,region(kr)%n_fe_subregions
          ss=region(kr)%fe_subregion(ks)
          do ke=1,part(fe_subregion(ss)%part)%n_elements
            se=part(fe_subregion(ss)%part)%element(ke)
              if ((element(se)%n_dimension.eq.2).and.(problem%n.eq.3)) then
                exp_n_elements=exp_n_elements+1
                exp_element_eid(exp_n_elements)=element(se)%id
                ! Local coordinate of the centroid in element reference coordinates
                select case (fbem_n_edges(element(se)%type))
                  case (3)
                    xi2d=1.d0/3.d0
                  case (4)
                    xi2d=0
                  case default
                    stop 'not valid n_edges'
                end select
                N=fbem_normal3d(element(se)%type,element(se)%x_gn,xi2d)
                T1=fbem_tangent_xi1(3,element(se)%type,element(se)%x_gn,xi2d)
                call fbem_fem_degshell_stress_resultants_ep(N,T1,element(se)%ep(:,1),ep1,ep2,ep3)
                exp_element_value_r(1:3,exp_n_elements)=ep3
              end if
          end do
        end do
      end if
    end do
    if (exp_n_elements.gt.1) then
      call fbem_export_gmsh_ElementData(output_fileunit,'Shell Z axis',0.d0,0,3,&
                                        n_elements,exp_n_elements,exp_element_eid,exp_element_value_r)
    end if

  end if

  if (kf.eq.1) then

    ! ==============================================================================================================================
    ! Flooded finite elements
    ! ==============================================================================================================================

    exp_n_elements=0
    exp_element_eid=0
    exp_element_value_r=0
    do kr=1,n_regions
      if (region(kr)%class.eq.fbem_fe) then
        do ks=1,region(kr)%n_fe_subregions
          ss=region(kr)%fe_subregion(ks)
          do ke=1,part(fe_subregion(ss)%part)%n_elements
            se=part(fe_subregion(ss)%part)%element(ke)
            exp_n_elements=exp_n_elements+1
            exp_element_eid(exp_n_elements)=element(se)%id
            if (element(se)%flooded) then
              exp_element_value_r(1,exp_n_elements)=1
            end if
          end do
        end do
      end if
    end do
    if (exp_n_elements.gt.1) then
      call fbem_export_gmsh_ElementData(output_fileunit,'Flooded elements',0.d0,0,1,&
                                        n_elements,exp_n_elements,exp_element_eid,exp_element_value_r)
    end if

    ! ==============================================================================================================================
    ! Submerged finite elements
    ! ==============================================================================================================================

    exp_n_elements=0
    exp_element_eid=0
    exp_element_value_r=0
    do kr=1,n_regions
      if (region(kr)%class.eq.fbem_fe) then
        do ks=1,region(kr)%n_fe_subregions
          ss=region(kr)%fe_subregion(ks)
          do ke=1,part(fe_subregion(ss)%part)%n_elements
            se=part(fe_subregion(ss)%part)%element(ke)
            exp_n_elements=exp_n_elements+1
            exp_element_eid(exp_n_elements)=element(se)%id
            if (element(se)%submerged) then
              exp_element_value_r(1,exp_n_elements)=1
            end if
          end do
        end do
      end if
    end do
    if (exp_n_elements.gt.1) then
      call fbem_export_gmsh_ElementData(output_fileunit,'Submerged elements',0.d0,0,1,&
                                        n_elements,exp_n_elements,exp_element_eid,exp_element_value_r)
    end if

  end if

  ! Deallocate working variable
  deallocate (exp_node_eid)
  deallocate (exp_node_value_r)
  deallocate (exp_node_value_c)
  deallocate (exp_element_eid)
  deallocate (exp_element_n_nodes)
  deallocate (exp_element_value_r)
  deallocate (exp_element_value_c)
  deallocate (exp_element_node_value_r)
  deallocate (exp_element_node_value_c)

  ! ============================================================================================================================
  ! INTERNAL ELEMENTS
  ! ============================================================================================================================

  !
  ! To-do: update and use gmsh exporting subroutines
  !

  if (internalelements) then
    if (trim(mesh_filename).eq.trim(internalelements_filename)) then

      ! Allocate working variable
      allocate (exp_node_eid(internalelements_mesh%n_nodes))
      allocate (exp_node_value_c(9,internalelements_mesh%n_nodes))
      allocate (exp_element_eid(internalelements_mesh%n_elements))
      allocate (exp_element_n_nodes(internalelements_mesh%n_elements))
      allocate (exp_element_node_value_c(9,maxval(fbem_n_nodes),internalelements_mesh%n_elements))

      ! ==========================================================================================================================
      ! Fluid pressure (internal elements)
      ! ==========================================================================================================================

      exp_n_elements=0
      exp_element_eid=0
      exp_element_n_nodes=0
      exp_element_node_value_c=0
      do kp=1,internalelements_mesh%n_parts
        kr=internalelements_mesh%part(kp)%entity
        if (kr.eq.0) cycle
        if (region(kr)%class.eq.fbem_be) then
          select case (region(kr)%type)
            case (fbem_potential)
                k_start=1
                cte = (1.d0,0.d0)
            case (fbem_viscoelastic)
                cycle
            case (fbem_poroelastic)
                k_start=0
                cte = -(1.d0,0.d0)/region(kr)%property_r(8)
          end select
          do ke=1,internalelements_mesh%part(kp)%n_elements
            se=internalelements_mesh%part(kp)%element(ke)
            exp_n_elements=exp_n_elements+1
            exp_element_eid(exp_n_elements)=internalelements_mesh%element(se)%id
            exp_element_n_nodes(exp_n_elements)=internalelements_mesh%element(se)%n_nodes
            do kn=1,exp_element_n_nodes(exp_n_elements)
              exp_element_node_value_c(1,kn,exp_n_elements)=cte*internalelements_mesh%element(se)%value_c(k_start,kn,0)
            end do
          end do
        end if
      end do
      if (exp_n_elements.gt.0) then
        !
        ! Write to file the real part
        !
        write(output_fileunit,'(a16)' ) '$ElementNodeData'
        write(output_fileunit,'(a1)' ) '1'
        write(output_fileunit,'(a31)') '"p^{total} (internal elements)"'
        write(output_fileunit,'(a1)' ) '1'
        write(fmt1,*) '(',fmt_real,')'
        call fbem_trim2b(fmt1)
        write(output_fileunit,fmt1) omega
        write(output_fileunit,'(a1)') '3'
        write(fmt1,*) '(',fmt_integer,')'
        call fbem_trim2b(fmt1)
        write(output_fileunit,fmt1) 2*kf-2
        write(output_fileunit,'(a1)') '1'
        write(fmt1,*) '(',fmt_integer,')'
        call fbem_trim2b(fmt1)
        write(output_fileunit,fmt1) exp_n_elements
        do k=1,exp_n_elements
          write(fmt1,*) '(2',fmt_integer,',',exp_element_n_nodes(k),fmt_real,')'
          call fbem_trim2b(fmt1)
          write(output_fileunit,fmt1) exp_element_eid(k), exp_element_n_nodes(k), real(exp_element_node_value_c(1,1:exp_element_n_nodes(k),k))
        end do
        write(output_fileunit,'(a19)') '$EndElementNodeData'
        !
        ! Write to file the imaginary part
        !
        write(output_fileunit,'(a16)' ) '$ElementNodeData'
        write(output_fileunit,'(a1)' ) '1'
        write(output_fileunit,'(a31)') '"p^{total} (internal elements)"'
        write(output_fileunit,'(a1)' ) '1'
        write(fmt1,*) '(',fmt_real,')'
        call fbem_trim2b(fmt1)
        write(output_fileunit,fmt1) omega
        write(output_fileunit,'(a1)') '3'
        write(fmt1,*) '(',fmt_integer,')'
        call fbem_trim2b(fmt1)
        write(output_fileunit,fmt1) 2*kf-1
        write(output_fileunit,'(a1)') '1'
        write(fmt1,*) '(',fmt_integer,')'
        call fbem_trim2b(fmt1)
        write(output_fileunit,fmt1) exp_n_elements
        do k=1,exp_n_elements
          write(fmt1,*) '(2',fmt_integer,',',exp_element_n_nodes(k),fmt_real,')'
          call fbem_trim2b(fmt1)
          write(output_fileunit,fmt1) exp_element_eid(k), exp_element_n_nodes(k), imag(exp_element_node_value_c(1,1:exp_element_n_nodes(k),k))
        end do
        write(output_fileunit,'(a19)') '$EndElementNodeData'
      end if

      ! ==========================================================================================================================
      ! Solid total displacements (internal elements)
      ! ==========================================================================================================================

      exp_n_elements=0
      exp_element_eid=0
      exp_element_n_nodes=0
      exp_element_node_value_c=0
      do kp=1,internalelements_mesh%n_parts
        kr=internalelements_mesh%part(kp)%entity
        if (kr.eq.0) cycle
        if (region(kr)%class.eq.fbem_be) then
          select case (region(kr)%type)
            case (fbem_potential)
                k_start=0
                k_end  =0
            case (fbem_viscoelastic)
                k_start=1
                k_end  =problem%n
            case (fbem_poroelastic)
                k_start=1
                k_end  =problem%n
          end select
          if (k_end.eq.0) cycle
          do ke=1,internalelements_mesh%part(kp)%n_elements
            se=internalelements_mesh%part(kp)%element(ke)
            exp_n_elements=exp_n_elements+1
            exp_element_eid(exp_n_elements)=internalelements_mesh%element(se)%id
            exp_element_n_nodes(exp_n_elements)=internalelements_mesh%element(se)%n_nodes
            do kn=1,exp_element_n_nodes(exp_n_elements)
              exp_element_node_value_c(1:problem%n,kn,exp_n_elements)=internalelements_mesh%element(se)%value_c(1:problem%n,kn,0)
            end do
          end do
        end if
      end do
      if (exp_n_elements.gt.0) then
        !
        ! Write to file the real part
        !
        write(output_fileunit,'(a16)' ) '$ElementNodeData'
        write(output_fileunit,'(a1)' ) '1'
        write(output_fileunit,'(a31)') '"u^{total} (internal elements)"'
        write(output_fileunit,'(a1)' ) '1'
        write(fmt1,*) '(',fmt_real,')'
        call fbem_trim2b(fmt1)
        write(output_fileunit,fmt1) omega
        write(output_fileunit,'(a1)') '3'
        write(fmt1,*) '(',fmt_integer,')'
        call fbem_trim2b(fmt1)
        write(output_fileunit,fmt1) 2*kf-2
        write(output_fileunit,'(a1)') '3'
        write(fmt1,*) '(',fmt_integer,')'
        call fbem_trim2b(fmt1)
        write(output_fileunit,fmt1) exp_n_elements
        do k=1,exp_n_elements
          write(fmt1,*) '(2',fmt_integer,',',3*exp_element_n_nodes(k),fmt_real,')'
          call fbem_trim2b(fmt1)
          write(output_fileunit,fmt1) exp_element_eid(k), exp_element_n_nodes(k), real(exp_element_node_value_c(1:3,1:exp_element_n_nodes(k),k))
        end do
        write(output_fileunit,'(a19)') '$EndElementNodeData'
        !
        ! Write to file the imaginary part
        !
        write(output_fileunit,'(a16)' ) '$ElementNodeData'
        write(output_fileunit,'(a1)' ) '1'
        write(output_fileunit,'(a31)') '"u^{total} (internal elements)"'
        write(output_fileunit,'(a1)' ) '1'
        write(fmt1,*) '(',fmt_real,')'
        call fbem_trim2b(fmt1)
        write(output_fileunit,fmt1) omega
        write(output_fileunit,'(a1)') '3'
        write(fmt1,*) '(',fmt_integer,')'
        call fbem_trim2b(fmt1)
        write(output_fileunit,fmt1) 2*kf-1
        write(output_fileunit,'(a1)') '3'
        write(fmt1,*) '(',fmt_integer,')'
        call fbem_trim2b(fmt1)
        write(output_fileunit,fmt1) exp_n_elements
        do k=1,exp_n_elements
          write(fmt1,*) '(2',fmt_integer,',',3*exp_element_n_nodes(k),fmt_real,')'
          call fbem_trim2b(fmt1)
          write(output_fileunit,fmt1) exp_element_eid(k), exp_element_n_nodes(k), imag(exp_element_node_value_c(1:3,1:exp_element_n_nodes(k),k))
        end do
        write(output_fileunit,'(a19)') '$EndElementNodeData'
      end if

      ! ==========================================================================================================================
      ! Solid total stress tensor (internal elements)
      ! ==========================================================================================================================

      if (problem%n.ne.3) then
        close(unit=output_fileunit)
        stop 'internal elements only 3D'
      end if

      exp_n_elements=0
      exp_element_eid=0
      exp_element_n_nodes=0
      exp_element_node_value_c=0
      do kp=1,internalelements_mesh%n_parts
        kr=internalelements_mesh%part(kp)%entity
        if (kr.eq.0) cycle
        if (region(kr)%class.eq.fbem_be) then
          select case (region(kr)%type)
            case (fbem_potential)
                k_start=0
                k_end  =0
            case (fbem_viscoelastic)
                k_start=1
                k_end  =problem%n
            case (fbem_poroelastic)
                k_start=1
                k_end  =problem%n
          end select
          if (k_end.eq.0) cycle
          do ke=1,internalelements_mesh%part(kp)%n_elements
            se=internalelements_mesh%part(kp)%element(ke)
            exp_n_elements=exp_n_elements+1
            exp_element_eid(exp_n_elements)=internalelements_mesh%element(se)%id
            exp_element_n_nodes(exp_n_elements)=internalelements_mesh%element(se)%n_nodes
            do kn=1,exp_element_n_nodes(exp_n_elements)
              exp_element_node_value_c(1:3,kn,exp_n_elements)=internalelements_mesh%element(se)%value_c(1,kn,1:problem%n)
              exp_element_node_value_c(4:6,kn,exp_n_elements)=internalelements_mesh%element(se)%value_c(2,kn,1:problem%n)
              exp_element_node_value_c(7:9,kn,exp_n_elements)=internalelements_mesh%element(se)%value_c(3,kn,1:problem%n)
            end do
          end do
        end if
      end do
      if (exp_n_elements.gt.0) then
        !
        ! Write to file the real part
        !
        write(output_fileunit,'(a16)' ) '$ElementNodeData'
        write(output_fileunit,'(a1)' ) '1'
        write(output_fileunit,'(a35)') '"sigma^{total} (internal elements)"'
        write(output_fileunit,'(a1)' ) '1'
        write(fmt1,*) '(',fmt_real,')'
        call fbem_trim2b(fmt1)
        write(output_fileunit,fmt1) omega
        write(output_fileunit,'(a1)') '3'
        write(fmt1,*) '(',fmt_integer,')'
        call fbem_trim2b(fmt1)
        write(output_fileunit,fmt1) 2*kf-2
        write(output_fileunit,'(a1)') '9'
        write(fmt1,*) '(',fmt_integer,')'
        call fbem_trim2b(fmt1)
        write(output_fileunit,fmt1) exp_n_elements
        do k=1,exp_n_elements
          write(fmt1,*) '(2',fmt_integer,',',9*exp_element_n_nodes(k),fmt_real,')'
          call fbem_trim2b(fmt1)
          write(output_fileunit,fmt1) exp_element_eid(k), exp_element_n_nodes(k), real(exp_element_node_value_c(:,1:exp_element_n_nodes(k),k))
        end do
        write(output_fileunit,'(a19)') '$EndElementNodeData'
        !
        ! Write to file the imaginary part
        !
        write(output_fileunit,'(a16)' ) '$ElementNodeData'
        write(output_fileunit,'(a1)' ) '1'
        write(output_fileunit,'(a35)') '"sigma^{total} (internal elements)"'
        write(output_fileunit,'(a1)' ) '1'
        write(fmt1,*) '(',fmt_real,')'
        call fbem_trim2b(fmt1)
        write(output_fileunit,fmt1) omega
        write(output_fileunit,'(a1)') '3'
        write(fmt1,*) '(',fmt_integer,')'
        call fbem_trim2b(fmt1)
        write(output_fileunit,fmt1) 2*kf-1
        write(output_fileunit,'(a1)') '9'
        write(fmt1,*) '(',fmt_integer,')'
        call fbem_trim2b(fmt1)
        write(output_fileunit,fmt1) exp_n_elements
        do k=1,exp_n_elements
          write(fmt1,*) '(2',fmt_integer,',',9*exp_element_n_nodes(k),fmt_real,')'
          call fbem_trim2b(fmt1)
          write(output_fileunit,fmt1) exp_element_eid(k), exp_element_n_nodes(k), imag(exp_element_node_value_c(:,1:exp_element_n_nodes(k),k))
        end do
        write(output_fileunit,'(a19)') '$EndElementNodeData'
      end if

      deallocate (exp_node_eid)
      deallocate (exp_node_value_c)
      deallocate (exp_element_eid)
      deallocate (exp_element_n_nodes)
      deallocate (exp_element_node_value_c)

    else
      call fbem_warning_message(error_unit,0,'',0,'internal element *.pos file only with the same mesh file as the main mesh')
    end if

  end if


end subroutine export_solution_mechanics_harmonic_gmsh
