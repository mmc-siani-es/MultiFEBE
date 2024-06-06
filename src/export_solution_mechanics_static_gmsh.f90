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

subroutine export_solution_mechanics_static_gmsh(output_fileunit)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_shape_functions
  use fbem_geometry
  use fbem_string_handling
  use fbem_fem_shells
  use fbem_numerical
  use fbem_gmsh

  ! Module of problem variables
  use problem_variables

  ! No implicit variables allowed
  implicit none

  ! I/O variables
  integer                                 :: output_fileunit
  ! Local variables
  integer                                 :: k
  integer                                 :: kr, kb, ke, kn, kip, sip, kc, ki, kj, kk, kp
  integer                                 :: sb, se, sn, ks, ss
  logical, allocatable                    :: node_used(:)
  character(len=fbem_fmtstr)              :: fmt1, fmt2
  integer                                 :: ncint, ncreal, nc, ncmax
  character(len=fbem_string_max_length)   :: tmp_string
  integer                                 :: face, k_start, k_end

  integer                                 :: exp_n_nodes
  integer, allocatable                    :: exp_node_eid(:)
  real(kind=real64), allocatable          :: exp_node_value(:,:)

  integer                                 :: exp_n_elements
  integer, allocatable                    :: exp_element_eid(:)
  real(kind=real64), allocatable          :: exp_element_value(:,:)

  integer, allocatable                    :: exp_element_n_nodes(:)
  real(kind=real64), allocatable          :: exp_element_node_value(:,:,:), exp_element_node_value_s(:,:)

  integer                                 :: kke, kkn, sse
  real(kind=real64)                       :: dpdx_tangential(problem%n), dudx_tangential(problem%n,problem%n), dudx_tangential_dil
  real(kind=real64)                       :: u_total(problem%n), sigma(3,3), eps(3,3)
  real(kind=real64)                       :: eigval(3), eigvec(3,3), eigwork(8)

  integer :: eiginfo, eiglwork

  real(kind=real64)                       :: sigma_local(3,3), dudx_tangential_local(problem%n,problem%n)
  real(kind=real64)                       :: mu, lambda, nu
  real(kind=real64), allocatable          :: psi(:,:), xi(:)
  real(kind=real64)                       :: delta
  real(kind=real64)                       :: local_axes(problem%n,problem%n)
  real(kind=real64)                       :: ep1(3), ep2(3), ep3(3), N(3), T1(3), xi2d(2), xi1d


  ! ################################################################################################################################
  !
  ! SOLUTION MESH
  !
  ! ################################################################################################################################

  ! Allocate working variable
  allocate (node_used(n_nodes))
  allocate (exp_node_eid(n_nodes))
  allocate (exp_node_value(9,n_nodes))
  allocate (exp_element_eid(n_elements))
  allocate (exp_element_value(9,n_elements))
  allocate (exp_element_n_nodes(n_elements))
  allocate (exp_element_node_value(9,maxval(fbem_n_nodes),n_elements))
  allocate (exp_element_node_value_s(9,n_elements))

  ! ================================================================================================================================
  ! Displacements (positive faces)
  ! ================================================================================================================================

  exp_n_nodes=0
  exp_node_eid=0
  exp_node_value=0
  node_used=.false.
  ! Loop through REGIONS
  do kr=1,n_regions
    select case (region(kr)%class)

      ! ============================================================================================================================
      ! BE REGION
      !
      case (fbem_be)

        ! +----------------+
        ! | BOUNDARY NODES |
        ! +----------------+

        k_start=1
        k_end  =problem%n

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
                        exp_node_value(1:problem%n,exp_n_nodes)=node(sn)%value_r(k_start:k_end,1)
                      !
                      ! Crack-like boundaries
                      !
                      case (fbem_boundary_class_cracklike)
                        exp_n_nodes=exp_n_nodes+1
                        exp_node_eid(exp_n_nodes)=node(sn)%id
                        exp_node_value(1:problem%n,exp_n_nodes)=node(sn)%value_r(k_start:k_end,1)
                    end select
                  !
                  ! BE-BE coupled boundary of BE-FE-BE coupled boundary
                  !
                  case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                    ! The region is the region 1 of the boundary
                    if (.not.region(kr)%boundary_reversion(kb)) then
                      exp_n_nodes=exp_n_nodes+1
                      exp_node_eid(exp_n_nodes)=node(sn)%id
                      exp_node_value(1:problem%n,exp_n_nodes)=node(sn)%value_r(k_start:k_end,1)
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
                exp_n_nodes=exp_n_nodes+1
                exp_node_eid(exp_n_nodes)=node(sn)%id
                exp_node_value(1:problem%n,exp_n_nodes)=node(sn)%value_r(1:problem%n,1)
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
        call fbem_export_gmsh_NodeData(output_fileunit,'Nodal displacements UX,UY (FE, BE +n)',0.d0,0,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value)
      case (3)
        call fbem_export_gmsh_NodeData(output_fileunit,'Nodal displacements UX,UY,UZ (FE, BE +n)',0.d0,0,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value)
      case default
        stop 'not valid problem%n'
    end select
  end if

  ! ================================================================================================================================
  ! Displacements (negative faces)
  ! ================================================================================================================================

  exp_n_nodes=0
  exp_node_eid=0
  exp_node_value=0
  node_used=.false.
  ! Loop through REGIONS
  do kr=1,n_regions
    select case (region(kr)%class)

      ! ============================================================================================================================
      ! BE REGION
      !
      case (fbem_be)

        ! +----------------+
        ! | BOUNDARY NODES |
        ! +----------------+

        k_start=1
        k_end  =problem%n

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
                        ! Nothing to do
                        !

                      !
                      ! Crack-like boundaries
                      !
                      case (fbem_boundary_class_cracklike)

                        exp_n_nodes=exp_n_nodes+1
                        exp_node_eid(exp_n_nodes)=node(sn)%id
                        exp_node_value(1:problem%n,exp_n_nodes)=node(sn)%value_r(k_start:k_end,2)

                    end select
                  !
                  ! BE-BE coupled boundary of BE-FE-BE coupled boundary
                  !
                  case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                    ! The region is the region 2 of the boundary
                    if (region(kr)%boundary_reversion(kb)) then
                      exp_n_nodes=exp_n_nodes+1
                      exp_node_eid(exp_n_nodes)=node(sn)%id
                      exp_node_value(1:problem%n,exp_n_nodes)=node(sn)%value_r(k_start:k_end,2)
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
  if (exp_n_nodes.gt.0) then
    select case (problem%n)
      case (2)
        call fbem_export_gmsh_NodeData(output_fileunit,'Nodal displacements UX,UY (BE -n)',0.d0,0,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value)
      case (3)
        call fbem_export_gmsh_NodeData(output_fileunit,'Nodal displacements UX,UY,UZ (BE -n)',0.d0,0,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value)
      case default
        stop 'not valid problem%n'
    end select
  end if

  ! ================================================================================================================================
  ! FE nodal forces/reactions
  ! ================================================================================================================================

  exp_n_nodes=0
  exp_node_eid=0
  exp_node_value=0
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
                exp_n_nodes=exp_n_nodes+1
                exp_node_eid(exp_n_nodes)=node(sn)%id
                exp_node_value(1:problem%n,exp_n_nodes)=node(sn)%value_r(1:problem%n,2)
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
        call fbem_export_gmsh_NodeData(output_fileunit,'Nodal forces FX,FY (FE)',0.d0,0,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value)
      case (3)
        call fbem_export_gmsh_NodeData(output_fileunit,'Nodal forces FX,FY,FZ (FE)',0.d0,0,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value)
      case default
        stop 'not valid problem%n'
    end select
  end if

!~   ! ================================================================================================================================
!~   ! FE nodal moment/reactions
!~   ! ================================================================================================================================

!~   exp_n_nodes=0
!~   exp_node_eid=0
!~   exp_node_value=0
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
!~               if (.not.node_used(sn)) then
!~                 node_used(sn)=.true.
!~                 if (node(sn)%n_dof.gt.problem%n) then
!~                   exp_n_nodes=exp_n_nodes+1
!~                   exp_node_eid(exp_n_nodes)=node(sn)%id
!~                   exp_node_value(1:problem%n,exp_n_nodes)=node(sn)%value_r((problem%n+1):(3*(problem%n-1)),2)
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
!~         call fbem_export_gmsh_NodeData(output_fileunit,'Nodal moment MZ (FE)',0.d0,0,1,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value)
!~       case (3)
!~         call fbem_export_gmsh_NodeData(output_fileunit,'Nodal moments MX,MY,MZ (FE)',0.d0,0,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value)
!~       case default
!~         stop 'not valid problem%n'
!~     end select
!~   end if



  !
  ! Falta solo las cargas equilibrantes (no creo que sean muy necesarias)
  !




  ! ================================================================================================================================
  ! Beam stress resultants (as a tensor with 9 components: NX,VY,BZ,0,0,0,0,0,0 or NX,VY,VZ,BX,BY,BZ,0,0,0)
  ! ================================================================================================================================

  exp_n_elements=0
  exp_element_eid=0
  exp_element_n_nodes=0
  exp_element_node_value=0
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
                exp_element_node_value(1:3*(problem%n-1),kn,exp_n_elements)=element(se)%value_r(1:3*(problem%n-1),kn,1)
              end do
            end if
        end do
      end do
    end if
  end do
  if (exp_n_elements.gt.1) then
    select case (problem%n)
      case (2)
        call fbem_export_gmsh_ElementNodeData(output_fileunit,'Beam NX,VY,MZ (FE)',0.d0,0,9,&
                                              n_elements,maxval(fbem_n_nodes),&
                                              exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value)
      case (3)
        call fbem_export_gmsh_ElementNodeData(output_fileunit,'Beam NX,VY,VZ,MX,MY,MZ (FE)',0.d0,0,9,&
                                              n_elements,maxval(fbem_n_nodes),&
                                              exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value)
      case default
        stop 'not valid problem%n'
    end select

  end if

  ! ================================================================================================================================
  ! Beam stress resultants local axis x
  ! ================================================================================================================================

  exp_n_elements=0
  exp_element_eid=0
  exp_element_value=0
  do kr=1,n_regions
    if (region(kr)%class.eq.fbem_fe) then
      do ks=1,region(kr)%n_fe_subregions
        ss=region(kr)%fe_subregion(ks)
        do ke=1,part(fe_subregion(ss)%part)%n_elements
          se=part(fe_subregion(ss)%part)%element(ke)
            if ((element(se)%n_dimension.eq.1).and.(element(se)%fe_type.eq.0)) then
              exp_n_elements=exp_n_elements+1
              exp_element_eid(exp_n_elements)=element(se)%id
              exp_element_value(1:3,exp_n_elements)=element(se)%v_midnode(:,1,1)
            end if
            if ((element(se)%n_dimension.eq.1).and.((element(se)%fe_type.eq.1).or.(element(se)%fe_type.eq.2))) then
              exp_n_elements=exp_n_elements+1
              exp_element_eid(exp_n_elements)=element(se)%id
              exp_element_value(1:3,exp_n_elements)=element(se)%ep(:,1)
            end if
        end do
      end do
    end if
  end do
  if (exp_n_elements.gt.1) then
    call fbem_export_gmsh_ElementData(output_fileunit,'Beam X axis',0.d0,0,3,&
                                      n_elements,exp_n_elements,exp_element_eid,exp_element_value)
  end if

  ! ================================================================================================================================
  ! Beam stress resultants local axis y
  ! ================================================================================================================================

  exp_n_elements=0
  exp_element_eid=0
  exp_element_value=0
  do kr=1,n_regions
    if (region(kr)%class.eq.fbem_fe) then
      do ks=1,region(kr)%n_fe_subregions
        ss=region(kr)%fe_subregion(ks)
        do ke=1,part(fe_subregion(ss)%part)%n_elements
          se=part(fe_subregion(ss)%part)%element(ke)
            if ((element(se)%n_dimension.eq.1).and.(element(se)%fe_type.eq.0)) then
              exp_n_elements=exp_n_elements+1
              exp_element_eid(exp_n_elements)=element(se)%id
              exp_element_value(1:3,exp_n_elements)=element(se)%v_midnode(:,2,1)
            end if
            if ((element(se)%n_dimension.eq.1).and.((element(se)%fe_type.eq.1).or.(element(se)%fe_type.eq.2))) then
              exp_n_elements=exp_n_elements+1
              exp_element_eid(exp_n_elements)=element(se)%id
              exp_element_value(1:3,exp_n_elements)=element(se)%ep(:,2)
            end if
        end do
      end do
    end if
  end do
  if (exp_n_elements.gt.1) then
    call fbem_export_gmsh_ElementData(output_fileunit,'Beam Y axis',0.d0,0,3,&
                                      n_elements,exp_n_elements,exp_element_eid,exp_element_value)
  end if

  ! ================================================================================================================================
  ! Beam stress resultants local axis z
  ! ================================================================================================================================

  if (problem%n.eq.3) then
    exp_n_elements=0
    exp_element_eid=0
    exp_element_value=0
    do kr=1,n_regions
      if (region(kr)%class.eq.fbem_fe) then
        do ks=1,region(kr)%n_fe_subregions
          ss=region(kr)%fe_subregion(ks)
          do ke=1,part(fe_subregion(ss)%part)%n_elements
            se=part(fe_subregion(ss)%part)%element(ke)
              if ((element(se)%n_dimension.eq.1).and.(element(se)%fe_type.eq.0)) then
                exp_n_elements=exp_n_elements+1
                exp_element_eid(exp_n_elements)=element(se)%id
                exp_element_value(1:3,exp_n_elements)=element(se)%v_midnode(:,3,1)
              end if
              if ((element(se)%n_dimension.eq.1).and.((element(se)%fe_type.eq.1).or.(element(se)%fe_type.eq.2))) then
                exp_n_elements=exp_n_elements+1
                exp_element_eid(exp_n_elements)=element(se)%id
                exp_element_value(1:3,exp_n_elements)=element(se)%ep(:,3)
              end if
          end do
        end do
      end if
    end do
    if (exp_n_elements.gt.1) then
      call fbem_export_gmsh_ElementData(output_fileunit,'Beam Z axis',0.d0,0,3,&
                                        n_elements,exp_n_elements,exp_element_eid,exp_element_value)
    end if
  end if

  ! ================================================================================================================================
  ! Bar stress resultants (as scalar)
  ! ================================================================================================================================

  exp_n_elements=0
  exp_element_eid=0
  exp_element_n_nodes=0
  exp_element_node_value=0
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
                exp_element_node_value(1,kn,exp_n_elements)=element(se)%value_r(1,kn,1)
              end do
            end if
        end do
      end do
    end if
  end do
  if (exp_n_elements.gt.1) then
    call fbem_export_gmsh_ElementNodeData(output_fileunit,'Bar NX',0.d0,0,1,&
                                          n_elements,maxval(fbem_n_nodes),&
                                          exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value)
  end if

  ! ================================================================================================================================
  ! Discrete spring-dashpot stress resultants (as scalar)
  ! ================================================================================================================================

  exp_n_elements=0
  exp_element_eid=0
  exp_element_n_nodes=0
  exp_element_node_value=0
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
                exp_element_node_value(1,kn,exp_n_elements)=element(se)%value_r(1,kn,1)
              end do
            end if
        end do
      end do
    end if
  end do
  if (exp_n_elements.gt.1) then
    call fbem_export_gmsh_ElementNodeData(output_fileunit,'Spring-dashpot NX',0.d0,0,1,&
                                          n_elements,maxval(fbem_n_nodes),&
                                          exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value)
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
  exp_element_node_value=0
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
                exp_element_node_value(1:8,kn,exp_n_elements)=element(se)%value_r(1:8,kn,1)
              end do
            end if
        end do
      end do
    end if
  end do
  if (exp_n_elements.gt.1) then
    call fbem_export_gmsh_ElementNodeData(output_fileunit,'Shell Nx,Ny,Nxy,Mx,My,Mxy,Vx,Vy',0.d0,0,9,&
                                          n_elements,maxval(fbem_n_nodes),&
                                          exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value)
  end if

  ! ================================================================================================================================
  ! Shell stress resultants local axis z
  ! ================================================================================================================================

  exp_n_elements=0
  exp_element_eid=0
  exp_element_value=0
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
              exp_element_value(1:3,exp_n_elements)=ep1
            end if
        end do
      end do
    end if
  end do
  if (exp_n_elements.gt.1) then
    call fbem_export_gmsh_ElementData(output_fileunit,'Shell X axis',0.d0,0,3,&
                                      n_elements,exp_n_elements,exp_element_eid,exp_element_value)
  end if

  ! ================================================================================================================================
  ! Shell stress resultants local axis z
  ! ================================================================================================================================

  exp_n_elements=0
  exp_element_eid=0
  exp_element_value=0
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
              exp_element_value(1:3,exp_n_elements)=ep2
            end if
        end do
      end do
    end if
  end do
  if (exp_n_elements.gt.1) then
    call fbem_export_gmsh_ElementData(output_fileunit,'Shell Y axis',0.d0,0,3,&
                                      n_elements,exp_n_elements,exp_element_eid,exp_element_value)
  end if

  ! ================================================================================================================================
  ! Shell stress resultants local axis z
  ! ================================================================================================================================

  exp_n_elements=0
  exp_element_eid=0
  exp_element_value=0
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
              exp_element_value(1:3,exp_n_elements)=ep3
            end if
        end do
      end do
    end if
  end do
  if (exp_n_elements.gt.1) then
    call fbem_export_gmsh_ElementData(output_fileunit,'Shell Z axis',0.d0,0,3,&
                                      n_elements,exp_n_elements,exp_element_eid,exp_element_value)
  end if

  ! ================================================================================================================================
  ! Tractions (positive faces)
  ! ================================================================================================================================

  exp_n_nodes=0
  exp_node_eid=0
  exp_node_value=0
  ! Loop through REGIONS
  do kr=1,n_regions
    select case (region(kr)%class)

      ! ============================================================================================================================
      ! BE REGION
      !
      case (fbem_be)

        ! +----------------+
        ! | BOUNDARY NODES |
        ! +----------------+

        k_start=problem%n+1
        k_end  =2*problem%n

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
              if (node_used(sn).eqv.(.false.)) then
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
                        exp_node_value(1:problem%n,exp_n_nodes)=node(sn)%value_r(k_start:k_end,1)
                      !
                      ! Crack-like boundaries
                      !
                      case (fbem_boundary_class_cracklike)
                        exp_n_nodes=exp_n_nodes+1
                        exp_node_eid(exp_n_nodes)=node(sn)%id
                        exp_node_value(1:problem%n,exp_n_nodes)=node(sn)%value_r(k_start:k_end,1)
                    end select
                  !
                  ! BE-BE coupled boundary of BE-FE-BE coupled boundary
                  !
                  case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                    ! The region is the region 1 of the boundary
                    if (region(kr)%boundary_reversion(kb).eqv.(.false.)) then
                      exp_n_nodes=exp_n_nodes+1
                      exp_node_eid(exp_n_nodes)=node(sn)%id
                      exp_node_value(1:problem%n,exp_n_nodes)=node(sn)%value_r(k_start:k_end,1)
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


      ! ============================================================================================================================

      ! ============================================================================================================================
      ! FE region
      !
      case (fbem_fe)

        !
        ! Nothing to do...
        !

      ! ============================================================================================================================

    end select
  end do ! Loop through REGIONS
  if (exp_n_nodes.gt.0) then
    select case (problem%n)
      case (2)
        call fbem_export_gmsh_NodeData(output_fileunit,'Traction TX,TY (BE +n)',0.d0,0,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value)
      case (3)
        call fbem_export_gmsh_NodeData(output_fileunit,'Traction TX,TY,TZ (BE +n)',0.d0,0,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value)
      case default
        stop 'not valid problem%n'
    end select
  end if

  ! ================================================================================================================================
  ! Tractions (negative faces)
  ! ================================================================================================================================

  exp_n_nodes=0
  exp_node_eid=0
  exp_node_value=0
  ! Loop through REGIONS
  do kr=1,n_regions
    select case (region(kr)%class)

      ! ============================================================================================================================
      ! BE REGION
      !
      case (fbem_be)

        ! +----------------+
        ! | BOUNDARY NODES |
        ! +----------------+

        k_start=problem%n+1
        k_end  =2*problem%n

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
              if (node_used(sn).eqv.(.false.)) then
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
                        exp_node_value(1:problem%n,exp_n_nodes)=node(sn)%value_r(k_start:k_end,2)

                    end select
                  !
                  ! BE-BE coupled boundary of BE-FE-BE coupled boundary
                  !
                  case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                    ! The region is the region 2 of the boundary
                    if (region(kr)%boundary_reversion(kb)) then
                      exp_n_nodes=exp_n_nodes+1
                      exp_node_eid(exp_n_nodes)=node(sn)%id
                      exp_node_value(1:problem%n,exp_n_nodes)=node(sn)%value_r(k_start:k_end,2)
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
  if (exp_n_nodes.gt.0) then
    select case (problem%n)
      case (2)
        call fbem_export_gmsh_NodeData(output_fileunit,'Traction TX,TY (BE -n)',0.d0,0,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value)
      case (3)
        call fbem_export_gmsh_NodeData(output_fileunit,'Traction TX,TY,TZ (BE -n)',0.d0,0,3,n_nodes,exp_n_nodes,exp_node_eid,exp_node_value)
      case default
        stop 'not valid problem%n'
    end select
  end if

  ! ================================================================================================================================
  ! Stress tensor (positive faces)
  ! ================================================================================================================================

  ! Stress recovery from displacement gradient taken from pages 203-204 the reference book of Telles, Brebbia and Wrobel (1984).

  exp_n_elements=0
  exp_element_eid=0
  exp_element_n_nodes=0
  exp_element_node_value=0
  ! Loop through REGIONS
  do kr=1,n_regions
    select case (region(kr)%class)

      ! ============================================================================================================================
      ! BE REGION
      !
      case (fbem_be)

        ! +----------------+
        ! | BOUNDARY NODES |
        ! +----------------+

        k_start=problem%n+1
        k_end  =2*problem%n
        mu=region(kr)%property_r(2)
        lambda=region(kr)%property_r(6)
        nu=region(kr)%property_r(3)

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
                    face=1
                    allocate(psi(element(se)%n_nodes,problem%n),xi(element(se)%n_dimension))
                    exp_n_elements=exp_n_elements+1
                    exp_element_eid(exp_n_elements)=element(se)%id
                    exp_element_n_nodes(exp_n_elements)=element(se)%n_nodes
                    do kn=1,element(se)%n_nodes
                      xi=element(se)%xi_fn(:,kn)
                      psi=fbem_psi(problem%n,element(se)%type_g,element(se)%x_gn,element(se)%type_f1,element(se)%delta_f,xi)
                      dudx_tangential=0
                      do kkn=1,element(se)%n_nodes
                        do kc=1,problem%n
                          dudx_tangential(kc,:)=dudx_tangential(kc,:)+psi(kkn,:)*node(element(se)%node(kkn))%value_r(kc,face)
                        end do
                      end do
                      ! Local cartesian axes
                      local_axes=fbem_local_tangential_cartesian_axes(problem%n,element(se)%type_g,element(se)%x_gn,region(kr)%boundary_reversion(kb),xi)
                      select case (problem%n)
                        !
                        ! 2D (plane strain or plane stress)
                        !
                        case (2)
                          sigma_local=0
                          sigma=0
                          ! Local stress tensor due to traction at boundary
                          sigma_local(2,1)=dot_product(node(element(se)%node(kn))%value_r(3:4,face),local_axes(:,1))
                          sigma_local(2,2)=dot_product(node(element(se)%node(kn))%value_r(3:4,face),local_axes(:,2))
                          sigma_local(1,2)=sigma_local(2,1)
                          ! Tangential displacement gradient in local cartesian coordinates
                          dudx_tangential_local=matmul(matmul(transpose(local_axes),dudx_tangential),local_axes)
                          ! Local stress tensor due to traction at boundary and tangential displacement
                          select case (problem%subtype)
                            case (fbem_mechanics_plane_strain)
                              sigma_local(1,1)=1.d0/(1.d0-nu)*(2.d0*mu*dudx_tangential_local(1,1)+nu*sigma_local(2,2))
                            case (fbem_mechanics_plane_stress)
                              sigma_local(1,1)=2.d0*(1.d0+nu)*mu*dudx_tangential_local(1,1)+nu*sigma_local(2,2)
                            case default
                              call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid 2D subtype')
                          end select
                          ! Stress tensor in global coordinates
                          sigma(1:2,1:2)=matmul(matmul(local_axes,sigma_local(1:2,1:2)),transpose(local_axes))
                          if (problem%subtype.eq.fbem_mechanics_plane_strain) then
                            sigma(3,3)=nu*(sigma(1,1)+sigma(2,2))
                          end if
                        !
                        ! 3D
                        !
                        case (3)
                          ! Local stress tensor due to traction at boundary
                          sigma_local(3,1)=dot_product(node(element(se)%node(kn))%value_r(4:6,face),local_axes(:,1))
                          sigma_local(3,2)=dot_product(node(element(se)%node(kn))%value_r(4:6,face),local_axes(:,2))
                          sigma_local(3,3)=dot_product(node(element(se)%node(kn))%value_r(4:6,face),local_axes(:,3))
                          sigma_local(1,3)=sigma_local(3,1)
                          sigma_local(2,3)=sigma_local(3,2)
                          ! Tangential displacement gradient in local cartesian coordinates
                          dudx_tangential_local=matmul(matmul(transpose(local_axes),dudx_tangential),local_axes)
                          ! Local stress tensor due to traction at boundary and tangential displacement
                          sigma_local(1,1)=1.d0/(1.d0-nu)*(nu*sigma_local(3,3)+2.d0*mu*(dudx_tangential_local(1,1)+nu*dudx_tangential_local(2,2)))
                          sigma_local(2,2)=1.d0/(1.d0-nu)*(nu*sigma_local(3,3)+2.d0*mu*(dudx_tangential_local(2,2)+nu*dudx_tangential_local(1,1)))
                          sigma_local(1,2)=mu*(dudx_tangential_local(1,2)+dudx_tangential_local(2,1))
                          sigma_local(2,1)=sigma_local(1,2)
                          ! Stress tensor in global coordinates
                          sigma=matmul(matmul(local_axes,sigma_local),transpose(local_axes))
                        case default
                          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid dimension')
                      end select
                      ! Add
                      do ki=1,3
                        do kj=1,3
                          kk=(ki-1)*3+kj
                          exp_element_node_value(kk,kn,exp_n_elements)=sigma(ki,kj)
                        end do
                      end do
                    end do
                    deallocate(xi,psi)
                  !
                  ! Crack-like boundaries
                  !
                  case (fbem_boundary_class_cracklike)
                    face=1
                    allocate(psi(element(se)%n_nodes,problem%n),xi(element(se)%n_dimension))
                    exp_n_elements=exp_n_elements+1
                    exp_element_eid(exp_n_elements)=element(se)%id
                    exp_element_n_nodes(exp_n_elements)=element(se)%n_nodes
                    do kn=1,element(se)%n_nodes
                      xi=element(se)%xi_fn(:,kn)
                      psi=fbem_psi(problem%n,element(se)%type_g,element(se)%x_gn,element(se)%type_f1,element(se)%delta_f,xi)
                      dudx_tangential=0
                      do kkn=1,element(se)%n_nodes
                        do kc=1,problem%n
                          dudx_tangential(kc,:)=dudx_tangential(kc,:)+psi(kkn,:)*node(element(se)%node(kkn))%value_r(kc,face)
                        end do
                      end do
                      ! Local cartesian axes
                      local_axes=fbem_local_tangential_cartesian_axes(problem%n,element(se)%type_g,element(se)%x_gn,.false.,xi)
                      select case (problem%n)
                        !
                        ! 2D (plane strain or plane stress)
                        !
                        case (2)
                          sigma_local=0
                          sigma=0
                          ! Local stress tensor due to traction at boundary
                          sigma_local(2,1)=dot_product(node(element(se)%node(kn))%value_r(3:4,face),local_axes(:,1))
                          sigma_local(2,2)=dot_product(node(element(se)%node(kn))%value_r(3:4,face),local_axes(:,2))
                          sigma_local(1,2)=sigma_local(2,1)
                          ! Tangential displacement gradient in local cartesian coordinates
                          dudx_tangential_local=matmul(matmul(transpose(local_axes),dudx_tangential),local_axes)
                          ! Local stress tensor due to traction at boundary and tangential displacement
                          select case (problem%subtype)
                            case (fbem_mechanics_plane_strain)
                              sigma_local(1,1)=1.d0/(1.d0-nu)*(2.d0*mu*dudx_tangential_local(1,1)+nu*sigma_local(2,2))
                            case (fbem_mechanics_plane_stress)
                              sigma_local(1,1)=2.d0*(1.d0+nu)*mu*dudx_tangential_local(1,1)+nu*sigma_local(2,2)
                            case default
                              call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid 2D subtype')
                          end select
                          ! Stress tensor in global coordinates
                          sigma(1:2,1:2)=matmul(matmul(local_axes,sigma_local(1:2,1:2)),transpose(local_axes))
                          if (problem%subtype.eq.fbem_mechanics_plane_strain) then
                            sigma(3,3)=nu*(sigma(1,1)+sigma(2,2))
                          end if
                        !
                        ! 3D
                        !
                        case (3)
                          ! Local stress tensor due to traction at boundary
                          sigma_local(3,1)=dot_product(node(element(se)%node(kn))%value_r(4:6,face),local_axes(:,1))
                          sigma_local(3,2)=dot_product(node(element(se)%node(kn))%value_r(4:6,face),local_axes(:,2))
                          sigma_local(3,3)=dot_product(node(element(se)%node(kn))%value_r(4:6,face),local_axes(:,3))
                          sigma_local(1,3)=sigma_local(3,1)
                          sigma_local(2,3)=sigma_local(3,2)
                          ! Tangential displacement gradient in local cartesian coordinates
                          dudx_tangential_local=matmul(matmul(transpose(local_axes),dudx_tangential),local_axes)
                          ! Local stress tensor due to traction at boundary and tangential displacement
                          sigma_local(1,1)=1.d0/(1.d0-nu)*(nu*sigma_local(3,3)+2.d0*mu*(dudx_tangential_local(1,1)+nu*dudx_tangential_local(2,2)))
                          sigma_local(2,2)=1.d0/(1.d0-nu)*(nu*sigma_local(3,3)+2.d0*mu*(dudx_tangential_local(2,2)+nu*dudx_tangential_local(1,1)))
                          sigma_local(1,2)=mu*(dudx_tangential_local(1,2)+dudx_tangential_local(2,1))
                          sigma_local(2,1)=sigma_local(1,2)
                          ! Stress tensor in global coordinates
                          sigma=matmul(matmul(local_axes,sigma_local),transpose(local_axes))
                        case default
                          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid dimension')
                      end select
                      ! Add
                      do ki=1,3
                        do kj=1,3
                          kk=(ki-1)*3+kj
                          exp_element_node_value(kk,kn,exp_n_elements)=sigma(ki,kj)
                        end do
                      end do
                    end do
                    deallocate(xi,psi)

                end select
              !
              ! BE-BE coupled boundary of BE-FE-BE coupled boundary
              !
              case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                ! The region is the region 1 of the boundary
                if (region(kr)%boundary_reversion(kb).eqv.(.false.)) then
                  face=1
                  allocate(psi(element(se)%n_nodes,problem%n),xi(element(se)%n_dimension))
                  exp_n_elements=exp_n_elements+1
                  exp_element_eid(exp_n_elements)=element(se)%id
                  exp_element_n_nodes(exp_n_elements)=element(se)%n_nodes
                  do kn=1,element(se)%n_nodes
                    xi=element(se)%xi_fn(:,kn)
                    psi=fbem_psi(problem%n,element(se)%type_g,element(se)%x_gn,element(se)%type_f1,element(se)%delta_f,xi)
                    dudx_tangential=0
                    do kkn=1,element(se)%n_nodes
                      do kc=1,problem%n
                        dudx_tangential(kc,:)=dudx_tangential(kc,:)+psi(kkn,:)*node(element(se)%node(kkn))%value_r(kc,face)
                      end do
                    end do
                    ! Local cartesian axes
                    local_axes=fbem_local_tangential_cartesian_axes(problem%n,element(se)%type_g,element(se)%x_gn,.false.,xi)
                    select case (problem%n)
                      !
                      ! 2D (plane strain or plane stress)
                      !
                      case (2)
                        sigma_local=0
                        sigma=0
                        ! Local stress tensor due to traction at boundary
                        sigma_local(2,1)=dot_product(node(element(se)%node(kn))%value_r(3:4,face),local_axes(:,1))
                        sigma_local(2,2)=dot_product(node(element(se)%node(kn))%value_r(3:4,face),local_axes(:,2))
                        sigma_local(1,2)=sigma_local(2,1)
                        ! Tangential displacement gradient in local cartesian coordinates
                        dudx_tangential_local=matmul(matmul(transpose(local_axes),dudx_tangential),local_axes)
                        ! Local stress tensor due to traction at boundary and tangential displacement
                        select case (problem%subtype)
                          case (fbem_mechanics_plane_strain)
                            sigma_local(1,1)=1.d0/(1.d0-nu)*(2.d0*mu*dudx_tangential_local(1,1)+nu*sigma_local(2,2))
                          case (fbem_mechanics_plane_stress)
                            sigma_local(1,1)=2.d0*(1.d0+nu)*mu*dudx_tangential_local(1,1)+nu*sigma_local(2,2)
                          case default
                            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid 2D subtype')
                        end select
                        ! Stress tensor in global coordinates
                        sigma(1:2,1:2)=matmul(matmul(local_axes,sigma_local(1:2,1:2)),transpose(local_axes))
                        if (problem%subtype.eq.fbem_mechanics_plane_strain) then
                          sigma(3,3)=nu*(sigma(1,1)+sigma(2,2))
                        end if
                      !
                      ! 3D
                      !
                      case (3)
                        ! Local stress tensor due to traction at boundary
                        sigma_local(3,1)=dot_product(node(element(se)%node(kn))%value_r(4:6,face),local_axes(:,1))
                        sigma_local(3,2)=dot_product(node(element(se)%node(kn))%value_r(4:6,face),local_axes(:,2))
                        sigma_local(3,3)=dot_product(node(element(se)%node(kn))%value_r(4:6,face),local_axes(:,3))
                        sigma_local(1,3)=sigma_local(3,1)
                        sigma_local(2,3)=sigma_local(3,2)
                        ! Tangential displacement gradient in local cartesian coordinates
                        dudx_tangential_local=matmul(matmul(transpose(local_axes),dudx_tangential),local_axes)
                        ! Local stress tensor due to traction at boundary and tangential displacement
                        sigma_local(1,1)=1.d0/(1.d0-nu)*(nu*sigma_local(3,3)+2.d0*mu*(dudx_tangential_local(1,1)+nu*dudx_tangential_local(2,2)))
                        sigma_local(2,2)=1.d0/(1.d0-nu)*(nu*sigma_local(3,3)+2.d0*mu*(dudx_tangential_local(2,2)+nu*dudx_tangential_local(1,1)))
                        sigma_local(1,2)=mu*(dudx_tangential_local(1,2)+dudx_tangential_local(2,1))
                        sigma_local(2,1)=sigma_local(1,2)
                        ! Stress tensor in global coordinates
                        sigma=matmul(matmul(local_axes,sigma_local),transpose(local_axes))
                      case default
                        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid dimension')
                    end select
                    ! Add
                    do ki=1,3
                      do kj=1,3
                        kk=(ki-1)*3+kj
                        exp_element_node_value(kk,kn,exp_n_elements)=sigma(ki,kj)
                      end do
                    end do
                  end do
                  deallocate(xi,psi)

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
    call fbem_export_gmsh_ElementNodeData(output_fileunit,'Stress tensor SIJ (BE +n)',0.d0,0,9,&
                                          n_elements,maxval(fbem_n_nodes),&
                                          exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value)
  end if

  ! ================================================================================================================================
  ! Stress tensor (negative faces)
  ! ================================================================================================================================

  exp_n_elements=0
  exp_element_eid=0
  exp_element_n_nodes=0
  exp_element_node_value=0
  ! Loop through REGIONS
  do kr=1,n_regions
    select case (region(kr)%class)

      ! ============================================================================================================================
      ! BE REGION
      !
      case (fbem_be)

        ! +----------------+
        ! | BOUNDARY NODES |
        ! +----------------+

        k_start=problem%n+1
        k_end  =2*problem%n
        mu=region(kr)%property_r(2)
        lambda=region(kr)%property_r(6)

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
                    ! Nothing to do
                    !

                  !
                  ! Crack-like boundaries
                  !
                  case (fbem_boundary_class_cracklike)
                    face=2
                    allocate(psi(element(se)%n_nodes,problem%n),xi(element(se)%n_dimension))
                    exp_n_elements=exp_n_elements+1
                    exp_element_eid(exp_n_elements)=element(se)%id
                    exp_element_n_nodes(exp_n_elements)=element(se)%n_nodes
                    do kn=1,element(se)%n_nodes
                      xi=element(se)%xi_fn(:,kn)
                      psi=fbem_psi(problem%n,element(se)%type_g,element(se)%x_gn,element(se)%type_f1,element(se)%delta_f,xi)
                      dudx_tangential=0
                      do kkn=1,element(se)%n_nodes
                        do kc=1,problem%n
                          dudx_tangential(kc,:)=dudx_tangential(kc,:)+psi(kkn,:)*node(element(se)%node(kkn))%value_r(kc,face)
                        end do
                      end do
                      ! Local cartesian axes
                      local_axes=fbem_local_tangential_cartesian_axes(problem%n,element(se)%type_g,element(se)%x_gn,.true.,xi)
                      select case (problem%n)
                        !
                        ! 2D (plane strain or plane stress)
                        !
                        case (2)
                          sigma_local=0
                          sigma=0
                          ! Local stress tensor due to traction at boundary
                          sigma_local(2,1)=dot_product(node(element(se)%node(kn))%value_r(3:4,face),local_axes(:,1))
                          sigma_local(2,2)=dot_product(node(element(se)%node(kn))%value_r(3:4,face),local_axes(:,2))
                          sigma_local(1,2)=sigma_local(2,1)
                          ! Tangential displacement gradient in local cartesian coordinates
                          dudx_tangential_local=matmul(matmul(transpose(local_axes),dudx_tangential),local_axes)
                          ! Local stress tensor due to traction at boundary and tangential displacement
                          select case (problem%subtype)
                            case (fbem_mechanics_plane_strain)
                              sigma_local(1,1)=1.d0/(1.d0-nu)*(2.d0*mu*dudx_tangential_local(1,1)+nu*sigma_local(2,2))
                            case (fbem_mechanics_plane_stress)
                              sigma_local(1,1)=2.d0*(1.d0+nu)*mu*dudx_tangential_local(1,1)+nu*sigma_local(2,2)
                            case default
                              call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid 2D subtype')
                          end select
                          ! Stress tensor in global coordinates
                          sigma(1:2,1:2)=matmul(matmul(local_axes,sigma_local(1:2,1:2)),transpose(local_axes))
                          if (problem%subtype.eq.fbem_mechanics_plane_strain) then
                            sigma(3,3)=nu*(sigma(1,1)+sigma(2,2))
                          end if
                        !
                        ! 3D
                        !
                        case (3)
                          ! Local stress tensor due to traction at boundary
                          sigma_local(3,1)=dot_product(node(element(se)%node(kn))%value_r(4:6,face),local_axes(:,1))
                          sigma_local(3,2)=dot_product(node(element(se)%node(kn))%value_r(4:6,face),local_axes(:,2))
                          sigma_local(3,3)=dot_product(node(element(se)%node(kn))%value_r(4:6,face),local_axes(:,3))
                          sigma_local(1,3)=sigma_local(3,1)
                          sigma_local(2,3)=sigma_local(3,2)
                          ! Tangential displacement gradient in local cartesian coordinates
                          dudx_tangential_local=matmul(matmul(transpose(local_axes),dudx_tangential),local_axes)
                          ! Local stress tensor due to traction at boundary and tangential displacement
                          sigma_local(1,1)=1.d0/(1.d0-nu)*(nu*sigma_local(3,3)+2.d0*mu*(dudx_tangential_local(1,1)+nu*dudx_tangential_local(2,2)))
                          sigma_local(2,2)=1.d0/(1.d0-nu)*(nu*sigma_local(3,3)+2.d0*mu*(dudx_tangential_local(2,2)+nu*dudx_tangential_local(1,1)))
                          sigma_local(1,2)=mu*(dudx_tangential_local(1,2)+dudx_tangential_local(2,1))
                          sigma_local(2,1)=sigma_local(1,2)
                          ! Stress tensor in global coordinates
                          sigma=matmul(matmul(local_axes,sigma_local),transpose(local_axes))
                        case default
                          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid dimension')
                      end select
                      ! Add
                      do ki=1,3
                        do kj=1,3
                          kk=(ki-1)*3+kj
                          exp_element_node_value(kk,kn,exp_n_elements)=sigma(ki,kj)
                        end do
                      end do
                    end do
                    deallocate(xi,psi)

                end select
              !
              ! BE-BE coupled boundary of BE-FE-BE coupled boundary
              !
              case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                ! The region is the region 1 of the boundary
                if (region(kr)%boundary_reversion(kb)) then
                  face=2
                  allocate(psi(element(se)%n_nodes,problem%n),xi(element(se)%n_dimension))
                  exp_n_elements=exp_n_elements+1
                  exp_element_eid(exp_n_elements)=element(se)%id
                  exp_element_n_nodes(exp_n_elements)=element(se)%n_nodes
                  do kn=1,element(se)%n_nodes
                    xi=element(se)%xi_fn(:,kn)
                    psi=fbem_psi(problem%n,element(se)%type_g,element(se)%x_gn,element(se)%type_f1,element(se)%delta_f,xi)
                    dudx_tangential=0
                    do kkn=1,element(se)%n_nodes
                      do kc=1,problem%n
                        dudx_tangential(kc,:)=dudx_tangential(kc,:)+psi(kkn,:)*node(element(se)%node(kkn))%value_r(kc,face)
                      end do
                    end do
                    ! Local cartesian axes
                    local_axes=fbem_local_tangential_cartesian_axes(problem%n,element(se)%type_g,element(se)%x_gn,.true.,xi)
                    select case (problem%n)
                      !
                      ! 2D (plane strain or plane stress)
                      !
                      case (2)
                        sigma_local=0
                        sigma=0
                        ! Local stress tensor due to traction at boundary
                        sigma_local(2,1)=dot_product(node(element(se)%node(kn))%value_r(3:4,face),local_axes(:,1))
                        sigma_local(2,2)=dot_product(node(element(se)%node(kn))%value_r(3:4,face),local_axes(:,2))
                        sigma_local(1,2)=sigma_local(2,1)
                        ! Tangential displacement gradient in local cartesian coordinates
                        dudx_tangential_local=matmul(matmul(transpose(local_axes),dudx_tangential),local_axes)
                        ! Local stress tensor due to traction at boundary and tangential displacement
                        select case (problem%subtype)
                          case (fbem_mechanics_plane_strain)
                            sigma_local(1,1)=1.d0/(1.d0-nu)*(2.d0*mu*dudx_tangential_local(1,1)+nu*sigma_local(2,2))
                          case (fbem_mechanics_plane_stress)
                            sigma_local(1,1)=2.d0*(1.d0+nu)*mu*dudx_tangential_local(1,1)+nu*sigma_local(2,2)
                          case default
                            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid 2D subtype')
                        end select
                        ! Stress tensor in global coordinates
                        sigma(1:2,1:2)=matmul(matmul(local_axes,sigma_local(1:2,1:2)),transpose(local_axes))
                        if (problem%subtype.eq.fbem_mechanics_plane_strain) then
                          sigma(3,3)=nu*(sigma(1,1)+sigma(2,2))
                        end if
                      !
                      ! 3D
                      !
                      case (3)
                        ! Local stress tensor due to traction at boundary
                        sigma_local(3,1)=dot_product(node(element(se)%node(kn))%value_r(4:6,face),local_axes(:,1))
                        sigma_local(3,2)=dot_product(node(element(se)%node(kn))%value_r(4:6,face),local_axes(:,2))
                        sigma_local(3,3)=dot_product(node(element(se)%node(kn))%value_r(4:6,face),local_axes(:,3))
                        sigma_local(1,3)=sigma_local(3,1)
                        sigma_local(2,3)=sigma_local(3,2)
                        ! Tangential displacement gradient in local cartesian coordinates
                        dudx_tangential_local=matmul(matmul(transpose(local_axes),dudx_tangential),local_axes)
                        ! Local stress tensor due to traction at boundary and tangential displacement
                        sigma_local(1,1)=1.d0/(1.d0-nu)*(nu*sigma_local(3,3)+2.d0*mu*(dudx_tangential_local(1,1)+nu*dudx_tangential_local(2,2)))
                        sigma_local(2,2)=1.d0/(1.d0-nu)*(nu*sigma_local(3,3)+2.d0*mu*(dudx_tangential_local(2,2)+nu*dudx_tangential_local(1,1)))
                        sigma_local(1,2)=mu*(dudx_tangential_local(1,2)+dudx_tangential_local(2,1))
                        sigma_local(2,1)=sigma_local(1,2)
                        ! Stress tensor in global coordinates
                        sigma=matmul(matmul(local_axes,sigma_local),transpose(local_axes))
                      case default
                        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid dimension')
                    end select
                    ! Add
                    do ki=1,3
                      do kj=1,3
                        kk=(ki-1)*3+kj
                        exp_element_node_value(kk,kn,exp_n_elements)=sigma(ki,kj)
                      end do
                    end do
                  end do
                  deallocate(xi,psi)
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
  if (exp_n_elements.gt.1) then
    call fbem_export_gmsh_ElementNodeData(output_fileunit,'Stress tensor SIJ (BE -n)',0.d0,0,9,&
                                          n_elements,maxval(fbem_n_nodes),&
                                          exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value)
  end if

  ! Deallocate working variable
  deallocate (exp_element_node_value,exp_element_n_nodes,exp_element_eid,exp_element_value)
  deallocate (exp_node_value,exp_node_eid,exp_element_node_value_s)

  ! ################################################################################################################################
  !
  ! POST-PROCESSING MESH
  !
  ! ################################################################################################################################

  ! ================================================================================================================================
  ! INTERNAL ELEMENTS
  ! ================================================================================================================================

  if (internalelements) then

    if (trim(mesh_filename).ne.trim(internalelements_filename)) then
      call fbem_error_message(error_unit,0,'',0,'internal element *.pos file only with the same mesh file as the main mesh')
    end if

    ! Allocate working variable
    allocate (exp_node_eid(internalelements_mesh%n_nodes))
    allocate (exp_node_value(9,internalelements_mesh%n_nodes))
    allocate (exp_element_eid(internalelements_mesh%n_elements))
    allocate (exp_element_n_nodes(internalelements_mesh%n_elements))
    allocate (exp_element_node_value(9,maxval(fbem_n_nodes),internalelements_mesh%n_elements))
    allocate (exp_element_node_value_s(maxval(fbem_n_nodes),internalelements_mesh%n_elements))

    ! ==============================================================================================================================
    ! Displacements (internal elements)
    ! ==============================================================================================================================

    exp_n_elements=0
    exp_element_eid=0
    exp_element_n_nodes=0
    exp_element_node_value=0
    do kp=1,internalelements_mesh%n_parts
      kr=internalelements_mesh%part(kp)%entity
      if (kr.eq.0) cycle
      if (region(kr)%class.eq.fbem_be) then
        k_start=1
        k_end  =problem%n
        do ke=1,internalelements_mesh%part(kp)%n_elements
          se=internalelements_mesh%part(kp)%element(ke)
          exp_n_elements=exp_n_elements+1
          exp_element_eid(exp_n_elements)=internalelements_mesh%element(se)%id
          exp_element_n_nodes(exp_n_elements)=internalelements_mesh%element(se)%n_nodes
          do kn=1,exp_element_n_nodes(exp_n_elements)
            exp_element_node_value(1:problem%n,kn,exp_n_elements)=internalelements_mesh%element(se)%value_r(1:problem%n,kn,0)
          end do
        end do
      end if
    end do
    call fbem_export_gmsh_ElementNodeData(output_fileunit,'u (BE IE)',0.d0,0,3,&
                                          internalelements_mesh%n_elements,27,&
                                          exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value)

    ! ==============================================================================================================================
    ! Stress tensor (internal elements)
    ! ==============================================================================================================================

    if (problem%n.ne.3) stop 'internal elements only 3D'

    exp_n_elements=0
    exp_element_eid=0
    exp_element_n_nodes=0
    exp_element_node_value=0
    do kp=1,internalelements_mesh%n_parts
      kr=internalelements_mesh%part(kp)%entity
      if (kr.eq.0) cycle
      if (region(kr)%class.eq.fbem_be) then
        k_start=1
        k_end  =problem%n
        mu=region(kr)%property_r(2)
        lambda=region(kr)%property_r(6)
        nu=region(kr)%property_r(3)
        do ke=1,internalelements_mesh%part(kp)%n_elements
          se=internalelements_mesh%part(kp)%element(ke)
          exp_n_elements=exp_n_elements+1
          exp_element_eid(exp_n_elements)=internalelements_mesh%element(se)%id
          exp_element_n_nodes(exp_n_elements)=internalelements_mesh%element(se)%n_nodes
          do kn=1,exp_element_n_nodes(exp_n_elements)
            exp_element_node_value(1:3,kn,exp_n_elements)=internalelements_mesh%element(se)%value_r(1,kn,1:problem%n)
            exp_element_node_value(4:6,kn,exp_n_elements)=internalelements_mesh%element(se)%value_r(2,kn,1:problem%n)
            exp_element_node_value(7:9,kn,exp_n_elements)=internalelements_mesh%element(se)%value_r(3,kn,1:problem%n)

            sigma(1,1)=exp_element_node_value(1,kn,exp_n_elements)
            sigma(1,2)=exp_element_node_value(2,kn,exp_n_elements)
            sigma(1,3)=exp_element_node_value(3,kn,exp_n_elements)
            sigma(2,1)=exp_element_node_value(4,kn,exp_n_elements)
            sigma(2,2)=exp_element_node_value(5,kn,exp_n_elements)
            sigma(2,3)=exp_element_node_value(6,kn,exp_n_elements)
            sigma(3,1)=exp_element_node_value(7,kn,exp_n_elements)
            sigma(3,2)=exp_element_node_value(8,kn,exp_n_elements)
            sigma(3,3)=exp_element_node_value(9,kn,exp_n_elements)

            eps(1,1)=(sigma(1,1)-nu*(sigma(2,2)+sigma(3,3)))/(2.d0*mu*(1.d0+nu))
            eps(1,2)=sigma(1,2)/(2.d0*mu)
            eps(1,3)=sigma(1,3)/(2.d0*mu)
            eps(2,1)=sigma(2,1)/(2.d0*mu)
            eps(2,2)=(sigma(2,2)-nu*(sigma(1,1)+sigma(3,3)))/(2.d0*mu*(1.d0+nu))
            eps(2,3)=sigma(2,3)/(2.d0*mu)
            eps(3,1)=sigma(3,1)/(2.d0*mu)
            eps(3,2)=sigma(3,2)/(2.d0*mu)
            eps(3,3)=(sigma(3,3)-nu*(sigma(1,1)+sigma(2,2)))/(2.d0*mu*(1.d0+nu))

            call dsyevj3(eps,eigvec,eigval) ! substitute with lapack dsyev

            exp_element_node_value_s(kn,exp_n_elements)=maxval(eigval)-minval(eigval)

          end do
        end do
      end if
    end do
    call fbem_export_gmsh_ElementNodeData(output_fileunit,'Stress tensor SIJ (BE IE)',0.d0,0,9,&
                                          internalelements_mesh%n_elements,maxval(fbem_n_nodes),&
                                          exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value)

    exp_element_node_value(1,:,:)=exp_element_node_value_s
    call fbem_export_gmsh_ElementNodeData(output_fileunit,'Maximum shear strain (BE IE)',0.d0,0,1,&
                                          internalelements_mesh%n_elements,maxval(fbem_n_nodes),&
                                          exp_n_elements,exp_element_eid,exp_element_n_nodes,exp_element_node_value)

  end if

end subroutine export_solution_mechanics_static_gmsh
