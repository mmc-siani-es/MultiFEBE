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


subroutine build_lse_laplace

 ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_numerical
  use fbem_shape_functions
  use fbem_geometry
  use fbem_symmetry
  use fbem_bem_general
  use fbem_bem_stapot2d
  use fbem_bem_stapot3d

  ! Module of problem variables
  use problem_variables

  ! No implicit variables
  implicit none

  ! Local variables
  integer                        :: kr
  integer                        :: kb_int, sb_int
  logical                        :: sb_int_reversion
  integer                        :: ks
  integer                        :: ke_int, se_int
  type(fbem_bem_element)         :: se_int_data
  integer                        :: se_int_type_g, se_int_type_f1, se_int_type_f2, se_int_n_nodes
  logical                        :: se_int_reversion
  real(kind=real64)              :: se_int_delta_f
  real(kind=real64), allocatable :: x_gn_int(:,:), xi_gn_int(:,:)
  real(kind=real64), allocatable :: x_fn_int(:,:), n_fn_int(:,:)
  integer                        :: kb_col, sb_col
  logical                        :: sb_col_reversion
  integer                        :: ke_col, se_col
  integer                        :: se_col_type_g, se_col_type_f1, se_col_type_f2, se_col_n_nodes
  real(kind=real64)              :: se_col_delta_f
  integer                        :: kn_col, sn_col
  integer                        :: kb, sb
  integer                        :: ke, se
  integer                        :: kn, knc, sn
  integer                        :: ss1, ss2
  integer                        :: kc
  real(kind=real64), allocatable :: xi_i(:), barxi(:)
  real(kind=real64)              :: x_i(problem%n), n_i(problem%n)
  ! Flags
  logical, allocatable           :: node_freeterm_added(:)
  logical, allocatable           :: node_collocated(:)
  ! Region properties
  real(kind=real64)              :: conductivity
  ! Vectors for SBIE integration
  real(kind=real64), allocatable :: h(:), g(:)
  real(kind=real64), allocatable :: hp(:), gp(:)
  real(kind=real64), allocatable :: hm(:), gm(:)
  real(kind=real64)              :: c_plus, c_minus
  ! Vectors for HBIE integration
  real(kind=real64), allocatable :: m(:), l(:)
  real(kind=real64), allocatable :: mp(:), lp(:)
  real(kind=real64), allocatable :: mm(:), lm(:)
  ! Associated with free-term calculation
  real(kind=real64), allocatable :: pphi_i(:)
  real(kind=real64), allocatable :: sphi_i(:)
  integer                        :: n_c_elements
  real(kind=real64), allocatable :: n_set_at_gn(:,:), n_set_at_gn_reversed(:,:)
  real(kind=real64), allocatable :: t_set_at_gn(:,:), t_set_at_gn_reversed(:,:)
  ! Associated with symmetry
  real(kind=real64), allocatable :: symconf_m(:), symconf_t(:), symconf_r(:)
  real(kind=real64)              :: symconf_s
  logical                        :: reversed
  ! Assembling
  logical                        :: assemble
  ! Writing
  character(len=fbem_fmtstr)     :: fmtstr

  ! Starting message
  if (verbose_level.ge.1)  write(output_unit,'(a)') 'Building the LSE ...'

  ! Initialize A matrix and b vector.
  A_r=0.d0
  b_r=0.d0

  ! ===========================
  ! BUILD BEM AND FEM EQUATIONS
  ! ===========================

  ! Allocate auxiliary variables
  allocate (node_freeterm_added(n_nodes))
  allocate (node_collocated(n_nodes))
  allocate (symconf_m(problem%n))
  allocate (symconf_t(problem%n))
  allocate (symconf_r(problem%n))

  ! Loop through REGIONS
  do kr=1,n_regions
    select case (region(kr)%class)

      ! ============================================================================================================================
      ! BE REGION
      !
      case (fbem_be)

        ! Save the conductivity
        conductivity=region(kr)%property_r(1)

        ! Message
        if (verbose_level.ge.2) then
          write(fmtstr,*) '(1x,a,i',fbem_nchar_int(region(kr)%id),')'
          call fbem_trimall(fmtstr)
          write(output_unit,fmtstr) 'Region: ', region(kr)%id
        end if

        ! Initialize the free-term control variable (for nodes with SBIE collocation)
        node_freeterm_added=.false.

        ! Loop through the BOUNDARIES of the REGION
        do kb_int=1,region(kr)%n_boundaries
          sb_int=region(kr)%boundary(kb_int)
          sb_int_reversion=region(kr)%boundary_reversion(kb_int)
          ! Message
          if (verbose_level.ge.3) then
            write(fmtstr,*) '(2x,a,i',fbem_nchar_int(boundary(sb_int)%id),')'
            call fbem_trimall(fmtstr)
            write(output_unit,fmtstr) 'Boundary (integration): ', boundary(sb_int)%id
          end if
          ! Loop through the ELEMENTS of the BOUNDARY
          do ke_int=1,part(boundary(sb_int)%part)%n_elements
            ! INTEGRATION ELEMENT
            se_int=part(boundary(sb_int)%part)%element(ke_int)
            se_int_type_g=element(se_int)%type_g
            se_int_type_f1=element(se_int)%type_f1
            se_int_type_f2=element(se_int)%type_f2
            se_int_delta_f=element(se_int)%delta_f
            se_int_n_nodes=element(se_int)%n_nodes
            ! Initialize calculation element
            call se_int_data%init
            se_int_data%gtype=se_int_type_g
            se_int_data%d=element(se_int)%n_dimension
            se_int_data%n_gnodes=se_int_n_nodes
            se_int_data%n=problem%n
            allocate (se_int_data%x(problem%n,se_int_n_nodes))
            se_int_data%x=element(se_int)%x_gn
            se_int_data%ptype=se_int_type_f1
            se_int_data%ptype_delta=se_int_delta_f
            se_int_data%n_pnodes=se_int_n_nodes
            se_int_data%stype=se_int_type_f2
            se_int_data%stype_delta=se_int_delta_f
            se_int_data%n_snodes=se_int_n_nodes
            se_int_data%cl=element(se_int)%csize
            se_int_data%gln_far=element(se_int)%n_phi
            allocate (se_int_data%bball_centre(problem%n))
            se_int_data%bball_centre=element(se_int)%bball_centre
            se_int_data%bball_radius=element(se_int)%bball_radius
            ! Allocate element-wise variables
            allocate (x_gn_int(problem%n,se_int_n_nodes),xi_gn_int(element(se_int)%n_dimension,se_int_n_nodes))
            allocate (x_fn_int(problem%n,se_int_n_nodes),n_fn_int(problem%n,se_int_n_nodes))
            allocate (h(se_int_n_nodes),g(se_int_n_nodes))
            allocate (hp(se_int_n_nodes),gp(se_int_n_nodes))
            allocate (hm(se_int_n_nodes),gm(se_int_n_nodes))
            allocate (m(se_int_n_nodes),l(se_int_n_nodes))
            allocate (mp(se_int_n_nodes),lp(se_int_n_nodes))
            allocate (mm(se_int_n_nodes),lm(se_int_n_nodes))
            allocate (xi_i(element(se_int)%n_dimension),barxi(element(se_int)%n_dimension))
            allocate (pphi_i(se_int_n_nodes),sphi_i(se_int_n_nodes))
            ! Save to local variables
            xi_gn_int=element(se_int)%xi_gn
            x_fn_int=element(se_int)%x_fn
            n_fn_int=element(se_int)%n_fn
            if (sb_int_reversion) n_fn_int=-n_fn_int

            ! Loop through SYMMETRICAL ELEMENTS
            do ks=1,n_symelements
              ! SYMMETRY SETUP
              call fbem_symmetry_multipliers(ks,problem%n,n_symplanes,symplane_m,symplane_s,symplane_t,symplane_r,&
                                             symconf_m,symconf_s,symconf_t,symconf_r,reversed)
              ! Change of element orientation and coordinates due to symmetry
              se_int_reversion=sb_int_reversion.neqv.reversed
              do kn=1,se_int_n_nodes
                se_int_data%x(:,kn)=symconf_m*element(se_int)%x_gn(:,kn)
              end do
              ! Initialize precalculated datasets
              call se_int_data%init_precalculated_datasets(n_precalsets,precalset_gln)

              ! =====================================
              ! BE BOUNDARY ELEMENT NODES COLLOCATION
              ! =====================================

              ! Initialize the collocation control variable (SBIE nodal collocation)
              node_collocated=.false.

              ! Loop through the BOUNDARIES of the REGION
              do kb_col=1,region(kr)%n_boundaries
                sb_col=region(kr)%boundary(kb_col)
                sb_col_reversion=region(kr)%boundary_reversion(kb_col)
                ! Message
                if (verbose_level.ge.4) then
                  write(fmtstr,*) '(3x,a,i',fbem_nchar_int(boundary(sb_col)%id),')'
                  call fbem_trimall(fmtstr)
                  write(output_unit,fmtstr) 'Boundary (collocation): ', boundary(sb_col)%id
                end if
                ! Loop through the ELEMENTS of the BOUNDARY
                do ke_col=1,part(boundary(sb_col)%part)%n_elements
                  ! COLLOCATION ELEMENT
                  se_col=part(boundary(sb_col)%part)%element(ke_col)
                  se_col_type_g=element(se_col)%type_g
                  se_col_type_f1=element(se_col)%type_f1
                  se_col_type_f2=element(se_col)%type_f2
                  se_col_delta_f=element(se_col)%delta_f
                  se_col_n_nodes=element(se_col)%n_nodes
                  ! Loop through the NODES of the ELEMENT
                  do kn_col=1,se_col_n_nodes
                    ! COLLOCATION NODE
                    sn_col=element(se_col)%node(kn_col)
                    ! Initialize assemble flag
                    assemble=.false.

                    ! ------------------------------------------------
                    ! STEPS:
                    ! 1 - CALCULATE SBIE KERNELS (H,G) if required
                    ! 2 - CALCULATE HBIE KERNELS (M,L) if required
                    ! 3 - CALCULATE SBIE+beta*HBIE KERNELS if required
                    ! 4 - ASSEMBLE the KERNELS to A and b
                    ! ------------------------------------------------

                    ! ----------------------
                    ! SBIE COLLOCATION POINT
                    ! ----------------------

                    if ((node(sn_col)%sbie.eq.fbem_sbie).and.(.not.node_collocated(sn_col))) then
                      ! INITIALIZE
                      assemble=.true.
                      node_collocated(sn_col)=.true.
                      ! CALCULATE KERNELS
                      x_i=element(se_col)%x_i_sbie(:,kn_col)
                      select case (problem%n)
                        case (2)
                          call fbem_bem_stapot2d_sbie_auto(se_int_data,se_int_reversion,x_i,qsi_parameters,qsi_ns_max,h,g)
                        case (3)
                          call fbem_bem_stapot3d_sbie_auto(se_int_data,se_int_reversion,x_i,qsi_parameters,qsi_ns_max,h,g)
                      end select
                      ! BUILD KERNELS ACCORDING TO SYMMETRY
                      if (ks.gt.1) then
                        h=symconf_s*h
                        g=symconf_s*g
                      end if
                      ! BUILD KERNELS WITH N+ AND N-
                      hp=h
                      gp=g
                      if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
                        hm=-h
                        gm= g
                      end if
                      ! CALCULATE AND ADD THE FREE-TERM
                      ! It is added only when collocation element == integration element, and the root element is being integrated.
                      if ((ks.eq.1).and.(se_int.eq.se_col)) then
                        xi_i=element(se_col)%xi_i_sbie(:,kn_col)
                        ! If its associated free-term has not been added
                        if (node_freeterm_added(sn_col).eqv.(.false.)) then
                          ! Mark that the free-term of the node is going to be calculated and added
                          node_freeterm_added(sn_col)=.true.
                          !
                          ! If the collocation point is in an edge or a vertex, the free-term has to be calculated.
                          !
                          if (fbem_check_xi_on_element_boundary(se_col_type_g,xi_i).eqv.(.true.)) then
                            !
                            ! Count the number of elements connected to the node for the integration region
                            !
                            ! Elements connected directly
                            n_c_elements=node(sn_col)%n_elements
                            ! Elements connected through common nodes
                            do kn=1,node(sn_col)%n_nodes
                              ! Selected common node
                              sn=node(sn_col)%node(kn)
                              ! If a "be" node
                              if (part(node(sn)%part(1))%type.eq.fbem_part_be_boundary) then
                                ! If the boundary of the node is in the integration region, then add the number of elements
                                ! of the common node to the counter
                                sb=part(node(sn)%part(1))%entity
                                do kb=1,region(kr)%n_boundaries
                                  if (region(kr)%boundary(kb).eq.sb) then
                                    n_c_elements=n_c_elements+node(sn)%n_elements
                                  end if
                                end do
                              end if
                            end do
                            ! If the collocation node belongs to any symmetry plane
                            select case (node(sn_col)%n_symplanes)
                              ! If it belongs to 1 symmetry plane
                              case (1)
                                n_c_elements=2*n_c_elements
                              ! If it belongs to 2 symmetry planes
                              case (2)
                                n_c_elements=4*n_c_elements
                            end select
                            ! Check
                            if ((problem%n.eq.2).and.(n_c_elements.gt.2)) then
                              call fbem_error_message(error_unit,0,'node',node(sn_col)%id,&
                                                      'is connected to more than 2 "be" elements.')
                            end if
                            ! Allocate the normals and tangents temporary matrices
                            allocate (n_set_at_gn(problem%n,n_c_elements))
                            allocate (t_set_at_gn(problem%n,n_c_elements))
                            allocate (n_set_at_gn_reversed(problem%n,n_c_elements))
                            allocate (t_set_at_gn_reversed(problem%n,n_c_elements))
                            !
                            ! Loop through the elements connected directly to the collocation node
                            !
                            ! Initialize the counter
                            n_c_elements=0
                            do ke=1,node(sn_col)%n_elements
                              ! Selected element
                              se=node(sn_col)%element(ke)
                              ! The index of the node in the selected element
                              kn=node(sn_col)%element_node_iid(ke)
                              ! Increment the counter
                              n_c_elements=n_c_elements+1
                              ! Copy the normal and the tangent with the appropiate sign
                              if (sb_col_reversion.eqv.(.false.)) then
                                do kc=1,problem%n
                                  n_set_at_gn(kc,n_c_elements)=element(se)%n_gn(kc,kn)
                                  t_set_at_gn(kc,n_c_elements)=element(se)%tbp_gn(kc,kn)
                                  n_set_at_gn_reversed(kc,n_c_elements)=-element(se)%n_gn(kc,kn)
                                  t_set_at_gn_reversed(kc,n_c_elements)=element(se)%tbm_gn(kc,kn)
                                end do
                              else
                                do kc=1,problem%n
                                  n_set_at_gn(kc,n_c_elements)=-element(se)%n_gn(kc,kn)
                                  t_set_at_gn(kc,n_c_elements)=element(se)%tbm_gn(kc,kn)
                                  n_set_at_gn_reversed(kc,n_c_elements)=element(se)%n_gn(kc,kn)
                                  t_set_at_gn_reversed(kc,n_c_elements)=element(se)%tbp_gn(kc,kn)
                                end do
                              end if
                            end do
                            !
                            ! Loop through common nodes
                            !
                            do kn=1,node(sn_col)%n_nodes
                              ! Selected common node
                              sn=node(sn_col)%node(kn)
                              ! If a "be" node
                              if (part(node(sn)%part(1))%type.eq.fbem_part_be_boundary) then
                                ! Copy the boundary of the selected common node
                                sb=part(node(sn)%part(1))%entity
                                ! Loop through the boundaries of the integration region
                                do kb=1,region(kr)%n_boundaries
                                  ! If the boundary of the selected common node is in the integration region, then copy all the
                                  ! normals and tangents of the elements of the common node
                                  if (region(kr)%boundary(kb).eq.sb) then
                                    ! Loop through the elements of the common node
                                    do ke=1,node(sn)%n_elements
                                      ! Selected element
                                      se=node(sn)%element(ke)
                                      ! The index of the common node in the selected element
                                      knc=node(sn)%element_node_iid(ke)
                                      ! Increment the counter
                                      n_c_elements=n_c_elements+1
                                      ! Copy the normal and the tangent with the appropiate sign
                                      if (region(kr)%boundary_reversion(kb).eqv.(.false.)) then
                                        do kc=1,problem%n
                                          n_set_at_gn(kc,n_c_elements)=element(se)%n_gn(kc,knc)
                                          t_set_at_gn(kc,n_c_elements)=element(se)%tbp_gn(kc,knc)
                                          n_set_at_gn_reversed(kc,n_c_elements)=-element(se)%n_gn(kc,knc)
                                          t_set_at_gn_reversed(kc,n_c_elements)=element(se)%tbm_gn(kc,knc)
                                        end do
                                      else
                                        do kc=1,problem%n
                                          n_set_at_gn(kc,n_c_elements)=-element(se)%n_gn(kc,knc)
                                          t_set_at_gn(kc,n_c_elements)=element(se)%tbm_gn(kc,knc)
                                          n_set_at_gn_reversed(kc,n_c_elements)=element(se)%n_gn(kc,knc)
                                          t_set_at_gn_reversed(kc,n_c_elements)=element(se)%tbp_gn(kc,knc)
                                        end do
                                      end if
                                    end do
                                  end if
                                end do
                              end if
                            end do
                            !
                            ! If the collocation node belongs to any symmetry plane, it is necessary to build the normals and
                            ! tangents of symmetrical elements.
                            !
                            select case (node(sn_col)%n_symplanes)
                              !
                              ! If it belongs to 1 symmetry plane (it can happen in 2D and 3D)
                              !
                              case (1)
                                ! Symmetry plane of the collocation node
                                ss1=node(sn_col)%symplane(1)
                                ! Loop through the original elements
                                do ke=1,n_c_elements
                                  ! Reflect the normal and the tangent with reversed orientation of each root element.
                                  do kc=1,problem%n
                                    n_set_at_gn(kc,ke+n_c_elements)=symplane_m(kc,ss1)*n_set_at_gn(kc,ke)
                                    t_set_at_gn(kc,ke+n_c_elements)=symplane_m(kc,ss1)*t_set_at_gn_reversed(kc,ke)
                                    n_set_at_gn_reversed(kc,ke+n_c_elements)=symplane_m(kc,ss1)*n_set_at_gn_reversed(kc,ke)
                                    t_set_at_gn_reversed(kc,ke+n_c_elements)=symplane_m(kc,ss1)*t_set_at_gn(kc,ke)
                                  end do
                                end do
                                ! Save the total number of elements
                                n_c_elements=2*n_c_elements
                              !
                              ! If it belongs to 2 symmetry planes (it can happen only in 3D)
                              !
                              case (2)
                                ! Symmetry planes of the collocation node
                                ss1=node(sn_col)%symplane(1)
                                ss2=node(sn_col)%symplane(2)
                                ! Loop through the original elements
                                do ke=1,n_c_elements
                                  ! Reflect the normal and the tangent with reversed orientation of each root element with
                                  ! respect to the first symmetry plane.
                                  do kc=1,3
                                    n_set_at_gn(kc,ke+n_c_elements)=symplane_m(kc,ss1)*n_set_at_gn(kc,ke)
                                    t_set_at_gn(kc,ke+n_c_elements)=symplane_m(kc,ss1)*t_set_at_gn_reversed(kc,ke)
                                    n_set_at_gn_reversed(kc,ke+n_c_elements)=symplane_m(kc,ss1)*n_set_at_gn_reversed(kc,ke)
                                    t_set_at_gn_reversed(kc,ke+n_c_elements)=symplane_m(kc,ss1)*t_set_at_gn(kc,ke)
                                  end do
                                  ! Reflect the normal and the tangent of each root element with respect to the first and
                                  ! the second symmetry planes.
                                  do kc=1,3
                                    n_set_at_gn(kc,ke+2*n_c_elements)=symplane_m(kc,ss1)*symplane_m(kc,ss2)*n_set_at_gn(kc,ke)
                                    t_set_at_gn(kc,ke+2*n_c_elements)=symplane_m(kc,ss1)*symplane_m(kc,ss2)*t_set_at_gn(kc,ke)
                                    n_set_at_gn_reversed(kc,ke+2*n_c_elements)=symplane_m(kc,ss1)*symplane_m(kc,ss2)*&
                                                                               n_set_at_gn_reversed(kc,ke)
                                    t_set_at_gn_reversed(kc,ke+2*n_c_elements)=symplane_m(kc,ss1)*symplane_m(kc,ss2)*&
                                                                               t_set_at_gn_reversed(kc,ke)
                                  end do
                                  ! Reflect the normal and the tangent with reversed orientation of each root element with
                                  ! respect to the second symmetry plane.
                                  do kc=1,3
                                    n_set_at_gn(kc,ke+3*n_c_elements)=symplane_m(kc,ss2)*n_set_at_gn(kc,ke)
                                    t_set_at_gn(kc,ke+3*n_c_elements)=symplane_m(kc,ss2)*t_set_at_gn_reversed(kc,ke)
                                    n_set_at_gn_reversed(kc,ke+3*n_c_elements)=symplane_m(kc,ss2)*n_set_at_gn_reversed(kc,ke)
                                    t_set_at_gn_reversed(kc,ke+3*n_c_elements)=symplane_m(kc,ss2)*t_set_at_gn(kc,ke)
                                  end do
                                end do
                                ! Save the total number of elements
                                n_c_elements=4*n_c_elements
                            end select
                            !
                            ! Depending on the problem dimension, calculate the free-term of h+
                            !
                            select case (problem%n)
                              case (2)
                                call fbem_bem_pot2d_sbie_freeterm(n_c_elements,n_set_at_gn,t_set_at_gn,geometric_tolerance,c_plus)
                              case (3)
                                call fbem_bem_pot3d_sbie_freeterm(n_c_elements,n_set_at_gn,t_set_at_gn,geometric_tolerance,c_plus)
                            end select
                            !
                            ! Depending on the problem dimension, calculate the free-term of h-
                            !
                            select case (problem%n)
                              case (2)
                                call fbem_bem_pot2d_sbie_freeterm(n_c_elements,n_set_at_gn_reversed,t_set_at_gn_reversed,&
                                                                  geometric_tolerance,c_minus)
                              case (3)
                                call fbem_bem_pot3d_sbie_freeterm(n_c_elements,n_set_at_gn_reversed,t_set_at_gn_reversed,&
                                                                  geometric_tolerance,c_minus)
                            end select
                            ! Deallocate temporary variables
                            deallocate (n_set_at_gn,t_set_at_gn)
                            deallocate (n_set_at_gn_reversed,t_set_at_gn_reversed)
                          !
                          ! If the node is not in an edge or vertex of the element
                          !
                          else
                            ! The free-term is 1/2
                            c_plus=0.5d0
                            ! The free-term is 1/2
                            c_minus=0.5d0
                          end if ! If the node is on an edge or vertex
                          !
                          ! Add the free-term to h+ in the collocation node
                          !
                          hp(kn_col)=hp(kn_col)+c_plus
                          ! Add the free-term matrix of the inverted elements to h- if the integration boundary is a crack-like boundary
                          if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
                              hm(kn_col)=hm(kn_col)+c_minus
                          end if
                        end if ! If the free-term hasn't been calculated and added
                      end if ! If ks==1 and collocation element == integration element
                    end if

                    ! --------------------------
                    ! SBIE MCA COLLOCATION POINT
                    ! --------------------------

                    if (node(sn_col)%sbie.eq.fbem_sbie_mca) then
                      ! INITIALIZE
                      assemble=.true.
                      ! CALCULATE KERNELS
                      x_i=element(se_col)%x_i_sbie_mca(:,kn_col)
                      select case (problem%n)
                        case (2)
                          call fbem_bem_stapot2d_sbie_auto(se_int_data,se_int_reversion,x_i,qsi_parameters,qsi_ns_max,h,g)
                        case (3)
                          call fbem_bem_stapot3d_sbie_auto(se_int_data,se_int_reversion,x_i,qsi_parameters,qsi_ns_max,h,g)
                      end select
                      ! BUILD KERNELS ACCORDING TO SYMMETRY
                      if (ks.gt.1) then
                        h=symconf_s*h
                        g=symconf_s*g
                      end if
                      ! BUILD KERNELS WITH N+ AND N-
                      hp=h
                      gp=g
                      if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
                        hm=-h
                        gm= g
                      end if
                      ! CALCULATE AND ADD THE FREE-TERM
                      if ((ks.eq.1).and.(se_int.eq.se_col)) then
                        xi_i=element(se_col)%xi_i_sbie_mca(:,kn_col)
                        pphi_i=fbem_phi_hybrid(se_int_type_f1,se_int_delta_f,xi_i)
                        hp=hp+0.5d0*pphi_i
                        if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
                          hm=hm+0.5d0*pphi_i
                        end if
                      end if
                    end if

                    ! ----------------------
                    ! HBIE COLLOCATION POINT
                    ! ----------------------

                    if (node(sn_col)%hbie.eq.fbem_hbie) then
                      ! INITIALIZE
                      assemble=.true.
                      ! CALCULATE KERNELS
                      x_i=element(se_col)%x_i_hbie(:,kn_col)
                      n_i=element(se_col)%n_i_hbie(:,kn_col)
                      if (sb_col_reversion) n_i=-n_i
                      select case (problem%n)
                        case (2)
                          call fbem_bem_stapot2d_hbie_auto(se_int_data,se_int_reversion,x_i,n_i,qsi_parameters,qsi_ns_max,m,l)
                        case (3)
                          call fbem_bem_stapot3d_hbie_auto(se_int_data,se_int_reversion,x_i,n_i,qsi_parameters,qsi_ns_max,m,l)
                      end select
                      ! BUILD KERNELS ACCORDING TO SYMMETRY
                      if (ks.gt.1) then
                        m=symconf_s*m
                        l=symconf_s*l
                      end if
                      ! BUILD KERNELS WITH N+ AND N-
                      mp=m
                      lp=l
                      if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
                        mm=-m
                        lm= l
                      end if
                      ! CALCULATE AND ADD THE FREE-TERM
                      if ((ks.eq.1).and.(se_int.eq.se_col)) then
                        xi_i=element(se_col)%xi_i_hbie(:,kn_col)
                        sphi_i=fbem_phi_hybrid(se_int_type_f2,se_int_delta_f,xi_i)
                        lp=lp-0.5d0*sphi_i
                        if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
                          lm=lm+0.5d0*sphi_i
                        end if
                      end if
                      ! SBIE = HBIE or SBIE + beta*HBIE if required
                      ! If dual member is 0, and the HBIE has been calculated, is because only the HBIE is used. Copy
                      ! m and l kernels to h and g respectively, since these are assembled.
                      if (node(sn_col)%dual.eq.0) then
                        hp=mp
                        gp=lp
                        if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
                          hm=mm
                          gm=lm
                        end if
                      end if
                      ! If dual member is fbem_dual_burton_miller, add beta*HBIE to h and g, since these are assembled.
                      if (node(sn_col)%dual.eq.fbem_dual_burton_miller) then
                        stop 'Error: fbem_dual_burton_miller not available'
                      end if
                    end if

                    ! --------
                    ! ASSEMBLE
                    ! --------

                    if (assemble) then
                      select case (boundary(sb_col)%coupling)
                        case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
                          select case (boundary(sb_col)%class)
                            case (fbem_boundary_class_ordinary)
                              gp=gp/conductivity
                              gm=gm/conductivity
                              call assemble_bem_laplace_equation(sb_int,sb_int_reversion,se_int,se_int_n_nodes,hp,gp,hm,gm,sn_col,1)
                            case (fbem_boundary_class_cracklike)
                              gp=gp/conductivity
                              gm=gm/conductivity
                              lp=lp/conductivity
                              lm=lm/conductivity
                              call assemble_bem_laplace_equation(sb_int,sb_int_reversion,se_int,se_int_n_nodes,hp,gp,hm,gm,sn_col,1)
                              call assemble_bem_laplace_equation(sb_int,sb_int_reversion,se_int,se_int_n_nodes,mp,lp,mm,lm,sn_col,2)
                          end select
                        case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                          gp=gp/conductivity
                          gm=gm/conductivity
                          if (sb_col_reversion.eqv.(.false.)) then
                            call assemble_bem_laplace_equation(sb_int,sb_int_reversion,se_int,se_int_n_nodes,hp,gp,hm,gm,sn_col,1)
                          else
                            call assemble_bem_laplace_equation(sb_int,sb_int_reversion,se_int,se_int_n_nodes,hp,gp,hm,gm,sn_col,2)
                          end if
                      end select
                    end if ! If assembling is needed

                  end do ! Loop through the NODES of the ELEMENT for COLLOCATION

                end do ! Loop through the ELEMENTS of the BOUNDARY for COLLOCATION

              end do ! Loop through the BOUNDARIES of the REGION for COLLOCATION

            end do ! Loop through SYMMETRICAL ELEMENTS

            ! Deallocate element-wise data structures
            deallocate (x_gn_int,xi_gn_int,x_fn_int,n_fn_int,xi_i,barxi,pphi_i,sphi_i)
            deallocate (h,g,hp,gp,hm,gm)
            deallocate (m,l,mp,lp,mm,lm)

          end do ! Loop through the ELEMENTS of the BOUNDARY for COLLOCATION

        end do ! Loop through the BOUNDARIES of the REGION for INTEGRATION

      ! ============================================================================================================================

      ! ============================================================================================================================
      ! FE REGION
      !
      case (fbem_fe)
        stop 'laplace > FEM: not yet'

      ! ============================================================================================================================

    end select ! Switch between REGION CLASSES

  end do ! Loop through REGIONS

  ! Ending message
  if (verbose_level.ge.1) write(output_unit,'(a)') 'done.'

end subroutine build_lse_laplace
