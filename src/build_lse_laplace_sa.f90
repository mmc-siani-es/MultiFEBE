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


subroutine build_lse_laplace_sa

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
  use fbem_bem_general
  use fbem_bem_stapot2d
  use fbem_bem_stapot3d

  ! Module of problem variables
  use problem_variables

  ! No implicit variables
  implicit none

  ! Local
  integer                        :: kr
  integer                        :: i, il, ik, ij, im, iq, ip
  integer                        :: kb_int, sb_int
  logical                        :: sb_int_reversion
  integer                        :: ks
  type(fbem_bem_element)         :: se_int_data
  integer                        :: ke_int, se_int
  integer                        :: se_int_type_g, se_int_type_f1, se_int_type_f2
  integer                        :: se_int_n_nodes
  logical                        :: se_int_reversion
  real(kind=real64)              :: se_int_delta_f
  logical                        :: se_int_sensitivity
  real(kind=real64), allocatable :: x_gn_int(:,:)
  real(kind=real64), allocatable :: xi_gn_int(:,:)
  real(kind=real64), allocatable :: x_fn_int(:,:), n_fn_int(:,:)
  integer                        :: kb_col, sb_col
  logical                        :: sb_col_reversion
  integer                        :: ke_col, se_col
  integer                        :: se_col_type_g, se_col_type_f1, se_col_type_f2
  integer                        :: se_col_n_nodes
  real(kind=real64)              :: se_col_delta_f
  logical                        :: se_col_sensitivity
  integer                        :: kn_col, sn_col
  integer                        :: kb, sb
  integer                        :: ke, se
  integer                        :: kn, knc, sn
  integer                        :: ki, si
  integer                        :: ss1, ss2 ! Selected symmetry plane
  integer                        :: kc
  real(kind=real64), allocatable :: xi_i(:)
  real(kind=real64), allocatable :: x_i(:)
  real(kind=real64), allocatable :: n_i(:)
  ! Dataset at integration points
  logical, allocatable           :: node_freeterm_added(:)       ! Vector of flags to know if free-term has been calculated and added
  logical, allocatable           :: node_collocated(:)           ! Vector of flags to know if collocation has been done
  ! Region properties
  real(kind=real64)              :: conductivity
  ! Kernels for VSBIE integration
  real(kind=real64), allocatable :: h1(:,:,:), h2(:,:), h3(:,:), g1(:,:,:), g2(:,:), g3(:,:)
  real(kind=real64), allocatable :: h1p(:,:,:), h2p(:,:), h3p(:,:), g1p(:,:,:), g2p(:,:), g3p(:,:)
  real(kind=real64), allocatable :: h1m(:,:,:), h2m(:,:), h3m(:,:), g1m(:,:,:), g2m(:,:), g3m(:,:)
  real(kind=real64)              :: bp(problem%n,problem%n), bm(problem%n,problem%n)
  ! Kernels for VHBIE integration
  real(kind=real64), allocatable :: m1(:,:,:), m2(:,:), m3(:,:), l1(:,:,:), l2(:,:), l3(:,:)
  real(kind=real64), allocatable :: m1p(:,:,:), m2p(:,:), m3p(:,:), l1p(:,:,:), l2p(:,:), l3p(:,:)
  real(kind=real64), allocatable :: m1m(:,:,:), m2m(:,:), m3m(:,:), l1m(:,:,:), l2m(:,:), l3m(:,:)
  ! Associated with free-term calculation
  real(kind=real64), allocatable :: pphi_i(:)
  real(kind=real64), allocatable :: sphi_i(:)
  real(kind=real64), allocatable :: wqj_i(:,:)
  integer                        :: n_c_elements
  real(kind=real64), allocatable :: n_set_at_gn(:,:), n_set_at_gn_reversed(:,:)
  real(kind=real64), allocatable :: t_set_at_gn(:,:), t_set_at_gn_reversed(:,:)
  ! Associated with the DME
  integer                        :: sdme_int
  integer                        :: sdme_int_n_nodes
  integer                        :: sdme_col_n_nodes
  real(kind=real64), allocatable :: dme_phi_i(:)
  real(kind=real64), allocatable :: dme_wqj_i(:,:)
  real(kind=real64), allocatable :: dme_dxda(:,:,:)
  real(kind=real64)              :: dxda_i(problem%n,problem%n_designvariables)
  real(kind=real64)              :: dnda_i(problem%n,problem%n_designvariables)
  ! Associated with symmetry
  real(kind=real64), allocatable :: symconf_m(:), symconf_t(:), symconf_r(:)
  real(kind=real64)              :: symconf_s
  logical                        :: reversed
  ! Assembling control variable
  logical                        :: assemble
  ! Writing
  character(len=fbem_fmtstr)     :: fmtstr

  ! Message
  if (verbose_level.ge.1)  write(output_unit,'(a)') 'Building the b vector for sensitivity analysis ...'

  ! Initialize bsa vector
  bsa_r=0.d0

  ! Allocate auxiliary variables
  allocate (node_freeterm_added(n_nodes))
  allocate (node_collocated(n_nodes))
  allocate (symconf_m(problem%n))
  allocate (symconf_t(problem%n))
  allocate (symconf_r(problem%n))
  allocate (x_i(problem%n))
  allocate (n_i(problem%n))

  !
  ! Loop through REGIONS
  !
  do kr=1,n_regions
    if (region(kr)%sensitivity.eqv.(.false.)) cycle
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

        ! Initialize the free-term control variable
        node_freeterm_added=.false.

        !
        ! Loop through the BOUNDARIES of the REGION for INTEGRATION
        !
        do kb_int=1,region(kr)%n_boundaries
          sb_int=region(kr)%boundary(kb_int)
          sb_int_reversion=region(kr)%boundary_reversion(kb_int)

          ! Message
          if (verbose_level.ge.3) then
            write(fmtstr,*) '(2x,a,i',fbem_nchar_int(boundary(sb_int)%id),')'
            call fbem_trimall(fmtstr)
            write(output_unit,fmtstr,advance='no') 'Integration boundary: ', boundary(sb_int)%id
            if (sb_int_reversion.eqv.(.true.)) then
              write(output_unit,'(a)') ' (reversed)'
            else
              write(output_unit,*)
            end if
          end if

          !
          ! Loop through the ELEMENTS of the BOUNDARY for INTEGRATION
          !
          do ke_int=1,part(boundary(sb_int)%part)%n_elements
            ! INTEGRATION ELEMENT
            se_int=part(boundary(sb_int)%part)%element(ke_int)
            se_int_type_g=element(se_int)%type_g
            se_int_type_f1=element(se_int)%type_f1
            se_int_type_f2=element(se_int)%type_f2
            se_int_delta_f=element(se_int)%delta_f
            se_int_n_nodes=element(se_int)%n_nodes
            se_int_sensitivity=element(se_int)%sensitivity
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
            ! Design element
            ! Isoparametric
            if (element(se_int)%dm_mode.eq.0) then
              se_int_data%dmetype=se_int_data%gtype
              se_int_data%dme_d=se_int_data%d
              se_int_data%dme_n_gnodes=se_int_data%n_gnodes
              allocate (se_int_data%dme_x(problem%n,se_int_n_nodes))
              se_int_data%dme_x=element(se_int)%x_gn
              se_int_data%dme_cl=se_int_data%cl
            ! Macro Element
            else
              sdme_int=element(se_int)%dm_element(1)
              se_int_data%dmetype=design_mesh%element(sdme_int)%type
              se_int_data%dme_d=design_mesh%element(sdme_int)%n_dimension
              se_int_data%dme_n_gnodes=design_mesh%element(sdme_int)%n_nodes
              allocate (se_int_data%dme_x(problem%n,se_int_data%dme_n_gnodes))
              se_int_data%dme_x=design_mesh%element(sdme_int)%x_gn
              se_int_data%dme_cl=design_mesh%element(sdme_int)%csize
            end if
            sdme_int_n_nodes=se_int_data%dme_n_gnodes
            ! Allocate element-wise variables
            allocate (x_gn_int(problem%n,se_int_n_nodes),xi_gn_int(element(se_int)%n_dimension,se_int_n_nodes))
            allocate (x_fn_int(problem%n,se_int_n_nodes),n_fn_int(problem%n,se_int_n_nodes))
            allocate (h1(sdme_int_n_nodes,se_int_n_nodes,problem%n),h2(se_int_n_nodes,problem%n),h3(se_int_n_nodes,problem%n))
            allocate (g1(sdme_int_n_nodes,se_int_n_nodes,problem%n),g2(se_int_n_nodes,problem%n),g3(se_int_n_nodes,problem%n))
            allocate (h1p(sdme_int_n_nodes,se_int_n_nodes,problem%n),h2p(se_int_n_nodes,problem%n),h3p(se_int_n_nodes,problem%n))
            allocate (g1p(sdme_int_n_nodes,se_int_n_nodes,problem%n),g2p(se_int_n_nodes,problem%n),g3p(se_int_n_nodes,problem%n))
            allocate (h1m(sdme_int_n_nodes,se_int_n_nodes,problem%n),h2m(se_int_n_nodes,problem%n),h3m(se_int_n_nodes,problem%n))
            allocate (g1m(sdme_int_n_nodes,se_int_n_nodes,problem%n),g2m(se_int_n_nodes,problem%n),g3m(se_int_n_nodes,problem%n))
            allocate (m1(sdme_int_n_nodes,se_int_n_nodes,problem%n),m2(se_int_n_nodes,problem%n),m3(se_int_n_nodes,problem%n))
            allocate (l1(sdme_int_n_nodes,se_int_n_nodes,problem%n),l2(se_int_n_nodes,problem%n),l3(se_int_n_nodes,problem%n))
            allocate (m1p(sdme_int_n_nodes,se_int_n_nodes,problem%n),m2p(se_int_n_nodes,problem%n),m3p(se_int_n_nodes,problem%n))
            allocate (l1p(sdme_int_n_nodes,se_int_n_nodes,problem%n),l2p(se_int_n_nodes,problem%n),l3p(se_int_n_nodes,problem%n))
            allocate (m1m(sdme_int_n_nodes,se_int_n_nodes,problem%n),m2m(se_int_n_nodes,problem%n),m3m(se_int_n_nodes,problem%n))
            allocate (l1m(sdme_int_n_nodes,se_int_n_nodes,problem%n),l2m(se_int_n_nodes,problem%n),l3m(se_int_n_nodes,problem%n))
            allocate (dme_dxda(problem%n,sdme_int_n_nodes,problem%n_designvariables))
            allocate (xi_i(element(se_int)%n_dimension),pphi_i(se_int_n_nodes),sphi_i(se_int_n_nodes))
            allocate (wqj_i(se_int_n_nodes,problem%n))
            allocate (dme_wqj_i(sdme_int_n_nodes,problem%n))
            ! Save to local variables
            xi_gn_int=element(se_int)%xi_gn
            x_fn_int=element(se_int)%x_fn
            n_fn_int=element(se_int)%n_fn
            if (sb_int_reversion) n_fn_int=-n_fn_int
            ! Design element
            ! Isoparametric
            if (element(se_int)%dm_mode.eq.0) then
              do kn=1,se_int_n_nodes
                dme_dxda(:,kn,:)=node(element(se_int)%node(kn))%dxda
              end do
            ! Macro Element
            else
              do kn=1,sdme_int_n_nodes
                dme_dxda(:,kn,:)=design_mesh%node(design_mesh%element(sdme_int)%node(kn))%dxda
              end do
            end if

            !  Message
            if (verbose_level.ge.2) then
              write(fmtstr,*) '(3x,a,i',fbem_nchar_int(element(se_int)%id),')'
              call fbem_trimall(fmtstr)
              write(output_unit,fmtstr) 'Integration element: ', element(se_int)%id
            end if

            ! Loop through symmetrical elements
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

              ! Loop through the BOUNDARIES of the REGION for COLLOCATION
              do kb_col=1,region(kr)%n_boundaries
                sb_col=region(kr)%boundary(kb_col)
                sb_col_reversion=region(kr)%boundary_reversion(kb_col)
                ! Message
                if (verbose_level.ge.3) then
                  write(fmtstr,*) '(4x,a,i',fbem_nchar_int(boundary(sb_col)%id),')'
                  call fbem_trimall(fmtstr)
                  write(output_unit,fmtstr,advance='no') 'Collocation boundary: ', boundary(sb_col)%id
                  if (sb_col_reversion.eqv.(.true.)) then
                    write(output_unit,'(a)') ' (reversed)'
                  else
                    write(output_unit,*)
                  end if
                end if

                ! Loop through the ELEMENTS of the BOUNDARY for COLLOCATION
                do ke_col=1,part(boundary(sb_col)%part)%n_elements
                  ! COLLOCATION ELEMENT
                  se_col=part(boundary(sb_col)%part)%element(ke_col)
                  se_col_type_g=element(se_col)%type_g
                  se_col_type_f1=element(se_col)%type_f1
                  se_col_type_f2=element(se_col)%type_f2
                  se_col_delta_f=element(se_col)%delta_f
                  se_col_n_nodes=element(se_col)%n_nodes
                  se_col_sensitivity=element(se_col)%sensitivity
                  ! ONLY BUILD IF INTEGRATION AND COLLOCATION ELEMENTS HAVE DXDA!=0
                  if (.not.(se_int_sensitivity.or.se_col_sensitivity)) cycle
                  ! COLLOCATION DESIGN ELEMENT
                  if (element(se_col)%dm_mode.eq.0) then
                    sdme_col_n_nodes=se_col_n_nodes
                  else
                    sdme_col_n_nodes=design_mesh%element(element(se_col)%dm_element(1))%n_nodes
                  end if

                  ! Loop through the NODES of the ELEMENT for COLLOCATION
                  do kn_col=1,se_col_n_nodes
                    ! COLLOCATION NODE
                    sn_col=element(se_col)%node(kn_col)
                    ! Initialize assemble flag
                    assemble=.false.

                    ! True for dual formulations (Burton & Miller and DBEM) when the collocation point for SBIE and HBIE is the same.
                    if (node(sn_col)%dual_is_common) then
                      stop 'not yet 48'
                    else

                      ! ======
                      !  SBIE
                      ! ======

                      ! If the collocation node has SBIE formulation, then the SBIE kernels of the integration element have to be
                      ! calculated.
                      if ((node(sn_col)%sbie.eq.fbem_sbie).and.(.not.node_collocated(sn_col))) then
                        ! INITIALIZE
                        assemble=.true.
                        node_collocated(sn_col)=.true.
                        ! CALCULATE KERNELS
                        x_i=element(se_col)%x_i_sbie(:,kn_col)
                        select case (problem%n)
                          case (2)
                            call fbem_bem_stapot2d_vsbie_auto(se_int_data,se_int_reversion,x_i,qsi_parameters,qsi_ns_max,h1,h2,g1,g2)
                          case (3)
                            stop 'not yet 49'
                        end select
                        h3=0.d0
                        g3=0.d0
                        ! BUILD KERNELS ACCORDING TO SYMMETRY
                        !
                        ! Nota: no valido, hay que contar con la simetria de la velocidad de las variables de diseño
                        !
                        if (ks.gt.1) then
                          stop 'symetry not yet for sensibility'
                          h1=symconf_s*h1
                          h2=symconf_s*h2
                          h3=symconf_s*h3
                          g1=symconf_s*g1
                          g2=symconf_s*g2
                          g3=symconf_s*g3
                        end if
                        ! BUILD KERNELS WITH N+ AND N-
                        h1p=h1
                        h2p=h2
                        h3p=h3
                        g1p=g1
                        g2p=g2
                        g3p=g3
                        if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
                          h1m=-h1
                          h2m=-h2
                          h3m=-h3
                          g1m= g1
                          g2m= g2
                          g3m= g3
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
                            if (fbem_check_xi_on_element_boundary(se_col_type_g,xi_i)) then
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
                                  call fbem_bem_pot2d_vsbie_freeterm(n_c_elements,n_set_at_gn,t_set_at_gn,geometric_tolerance,bp)
                                case (3)
                                  stop 'not yet 50'
                              end select
                              !
                              ! Depending on the problem dimension, calculate the free-term of h-
                              !
                              select case (problem%n)
                                case (2)
                                  call fbem_bem_pot2d_vsbie_freeterm(n_c_elements,n_set_at_gn_reversed,t_set_at_gn_reversed,&
                                                                     geometric_tolerance,bm)
                                case (3)
                                  stop 'not yet 51'
                              end select
                              ! Deallocate temporary variables
                              deallocate (n_set_at_gn,t_set_at_gn)
                              deallocate (n_set_at_gn_reversed,t_set_at_gn_reversed)
                              ! Build
                              dme_wqj_i=fbem_dphidx(problem%n,se_int_data%dmetype,0.d0,se_int_data%dme_x,element(se_col)%dm_xi_i_sbie(:,kn_col))
                              ! Add
                              do im=1,problem%n
                                do ij=1,problem%n
                                  h1p(:,kn_col,im)=h1p(:,kn_col,im)+dme_wqj_i(:,ij)*bp(im,ij)
                                end do
                              end do
                              if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
                                do im=1,problem%n
                                  do ij=1,problem%n
                                    h1m(:,kn_col,im)=h1m(:,kn_col,im)+dme_wqj_i(:,ij)*bm(im,ij)
                                  end do
                                end do
                              end if
                            end if
                          end if ! If the free-term hasn't been calculated and added
                        end if ! If ks==1 and collocation element == integration element
                        ! DXDA AND DNDA AT THE COLLOCATION POINT
                        dxda_i=element(se_col)%dxda_i_sbie(:,kn_col,:)
                        dnda_i=0.d0
                      end if ! If the SBIE is required

                      ! ==========
                      !  SBIE MCA
                      ! ==========

                      ! If the collocation node has SBIE MCA formulation, then the SBIE MCA kernels of the integration element have to be
                      ! calculated.
                      if (node(sn_col)%sbie.eq.fbem_sbie_mca) then
                        ! INITIALIZE
                        assemble=.true.
                        ! CALCULATE KERNELS
                        x_i=element(se_col)%x_i_sbie_mca(:,kn_col)
                        select case (problem%n)
                          case (2)
                            call fbem_bem_stapot2d_vsbie_auto(se_int_data,se_int_reversion,x_i,qsi_parameters,qsi_ns_max,h1,h2,g1,g2)
                          case (3)
                            stop 'not yet 52'
                        end select
                        h3=0.d0
                        g3=0.d0
                        ! BUILD KERNELS ACCORDING TO SYMMETRY
                        !
                        ! Nota: no valido, hay que contar con la simetria de la velocidad de las variables de diseño
                        !
                        if (ks.gt.1) then
                          stop 'symetry not yet for sensibility'
                          h1=symconf_s*h1
                          h2=symconf_s*h2
                          h3=symconf_s*h3
                          g1=symconf_s*g1
                          g2=symconf_s*g2
                          g3=symconf_s*g3
                        end if
                        ! BUILD KERNELS WITH N+ AND N-
                        h1p=h1
                        h2p=h2
                        h3p=h3
                        g1p=g1
                        g2p=g2
                        g3p=g3
                        if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
                          h1m=-h1
                          h2m=-h2
                          h3m=-h3
                          g1m= g1
                          g2m= g2
                          g3m= g3
                        end if
                        ! DXDA AND DNDA AT THE COLLOCATION POINT
                        dxda_i=element(se_col)%dxda_i_sbie_mca(:,kn_col,:)
                        dnda_i=0.d0
                      end if ! If the SBIE MCA is required

                      ! ======
                      !  HBIE
                      ! ======

                      ! If the collocation node has HBIE (MCA) formulation, then the HBIE equation has to be integrated.
                      if (node(sn_col)%hbie.eq.fbem_hbie) then
                        ! INITIALIZE
                        assemble=.true.
                        ! CALCULATE KERNELS
                        x_i=element(se_col)%x_i_hbie(:,kn_col)
                        n_i=element(se_col)%n_i_hbie(:,kn_col)
                        if (sb_col_reversion) n_i=-n_i
                        select case (problem%n)
                          case (2)
                            call fbem_bem_stapot2d_vhbie_auto(se_int_data,se_int_reversion,x_i,n_i,qsi_parameters,qsi_ns_max,m1,m2,m3,l1,l2,l3)
                          case (3)
                            stop 'not yet 53'
                        end select
                        ! BUILD KERNELS ACCORDING TO SYMMETRY
                        !
                        ! Nota: no valido, hay que contar con la simetria de la velocidad de las variables de diseño
                        !
                        if (ks.gt.1) then
                          stop 'symetry not yet for sensibility'
                          m1=symconf_s*m1
                          m2=symconf_s*m2
                          m3=symconf_s*m3
                          l1=symconf_s*l1
                          l2=symconf_s*l2
                          l3=symconf_s*l3
                        end if
                        ! BUILD KERNELS WITH N+ AND N-
                        m1p=m1
                        m2p=m2
                        m3p=m3
                        l1p=l1
                        l2p=l2
                        l3p=l3
                        if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
                          m1m=-m1
                          m2m=-m2
                          m3m=-m3
                          l1m=l1
                          l2m=l2
                          l3m=l3
                        end if






                        ! CALCULATE AND ADD THE FREE-TERM
                        if ((ks.eq.1).and.(se_int.eq.se_col)) then




!                          ! Build free-term
!                          dme_wqj_i=fbem_dphidx(problem%n,se_int_data%dmetype,0.d0,se_int_data%dme_x,element(se_col)%dm_xi_i_hbie(:,kn_col))
!                          wqj_i=fbem_dphidx(problem%n,se_int_data%ptype,se_int_data%ptype_delta,se_int_data%x,element(se_col)%xi_i_hbie(:,kn_col))
!                          ! Add
!                          do iq=1,sdme_int_n_nodes
!                            do ip=1,se_int_n_nodes
!                              do im=1,problem%n
!                                m1p(iq,ip,im)=m1p(iq,ip,im)+0.25d0*wqj_i(ip,im)*dot_product(dme_wqj_i(iq,:),n_i(:))
!                              end do
!                            end do
!                          end do

!                          ! FALTA PARA CRACKS
!                          if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
!                            do im=1,problem%n
!                              do ij=1,problem%n
!                                m1m(:,kn_col,im)=m1m(:,kn_col,im)+dme_wqj_i(:,ij)*bm(im,ij)
!                              end do
!                            end do
!                          end if



                          ! Ignoramos los free-terms de segunda derivada del dme field.....



                        end if





                        ! VSBIE = VHBIE or VSBIE + beta*VHBIE if required
                        ! If dual member is 0, and the HBIE has been calculated, is because only the HBIE is used. Copy
                        ! m and l kernels to h and g respectively, since these are assembled.
                        if (node(sn_col)%dual.eq.0) then
                          h1p=m1p
                          h2p=m2p
                          h3p=m3p
                          g1p=l1p
                          g2p=l2p
                          g3p=l3p
                          if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
                            h1m=m1m
                            h2m=m2m
                            h3m=m3m
                            g1m=l1m
                            g2m=l2m
                            g3m=l3m
                          end if
                        end if
                        ! If dual member is fbem_dual_burton_miller, add beta*HBIE to h and g, since these are assembled.
                        if (node(sn_col)%dual.eq.fbem_dual_burton_miller) then
                          stop 'Error: fbem_dual_burton_miller not available'
                        end if
                        ! DXDA AND DNDA AT THE COLLOCATION POINT
                        dxda_i=element(se_col)%dxda_i_hbie(:,kn_col,:)
                        dnda_i=element(se_col)%dnda_i_hbie(:,kn_col,:)
                        if (sb_col_reversion) dnda_i=-dnda_i
                      end if ! If the HBIE is required

                    end if

                    ! ========
                    ! ASSEMBLE
                    ! ========

                    !
                    ! Cuidadito con la colocacion para esto, cuando se usen duales, solo valido para colocar en el mismo punto,
                    ! por aquello del dme_phi_i, y solo no nodal.....
                    !

                    if (assemble) then
                      select case (boundary(sb_col)%coupling)
                        case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
                          select case (boundary(sb_col)%class)
                            case (fbem_boundary_class_ordinary)
                              g1p=g1p/conductivity
                              g2p=g2p/conductivity
                              g3p=g3p/conductivity
                              call assemble_bem_laplace_equation_sa(kr,sb_int,sb_int_reversion,&
                                                                    se_int,se_int_n_nodes,sdme_int_n_nodes,&
                                                                    h1p,h2p,h3p,g1p,g2p,g3p,h1m,h2m,h3m,g1m,g2m,g3m,dme_dxda,&
                                                                    dxda_i,dnda_i,sn_col,1)
                            case (fbem_boundary_class_cracklike)
                              g1p=g1p/conductivity
                              g2p=g2p/conductivity
                              g3p=g3p/conductivity
                              l1p=l1p/conductivity
                              l2p=l2p/conductivity
                              l3p=l3p/conductivity
                              call assemble_bem_laplace_equation_sa(kr,sb_int,sb_int_reversion,&
                                                                    se_int,se_int_n_nodes,sdme_int_n_nodes,&
                                                                    h1p,h2p,h3p,g1p,g2p,g3p,h1m,h2m,h3m,g1m,g2m,g3m,dme_dxda,&
                                                                    dxda_i,dnda_i,sn_col,1)
                              call assemble_bem_laplace_equation_sa(kr,sb_int,sb_int_reversion,&
                                                                    se_int,se_int_n_nodes,sdme_int_n_nodes,&
                                                                    m1p,m2p,m3p,l1p,l2p,l3p,m1m,m2m,m3m,l1m,l2m,l3m,dme_dxda,&
                                                                    dxda_i,dnda_i,sn_col,2)
                          end select
                        case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                          g1p=g1p/conductivity
                          g2p=g2p/conductivity
                          g3p=g3p/conductivity
                          if (sb_col_reversion.eqv.(.false.)) then
                            call assemble_bem_laplace_equation_sa(kr,sb_int,sb_int_reversion,&
                                                                  se_int,se_int_n_nodes,sdme_int_n_nodes,&
                                                                  h1p,h2p,h3p,g1p,g2p,g3p,h1m,h2m,h3m,g1m,g2m,g3m,dme_dxda,&
                                                                  dxda_i,dnda_i,sn_col,1)
                          else
                            call assemble_bem_laplace_equation_sa(kr,sb_int,sb_int_reversion,&
                                                                  se_int,se_int_n_nodes,sdme_int_n_nodes,&
                                                                  h1p,h2p,h3p,g1p,g2p,g3p,h1m,h2m,h3m,g1m,g2m,g3m,dme_dxda,&
                                                                  dxda_i,dnda_i,sn_col,2)
                          end if
                      end select
                    end if

                  end do ! Loop through the NODES of the ELEMENT for COLLOCATION

                end do ! Loop through the ELEMENTS of the BOUNDARY for COLLOCATION

              end do ! Loop through the BOUNDARIES of the REGION for COLLOCATION

            end do ! Loop through SYMMETRICAL ELEMENTS

            ! Deallocate element-wise data structures
            deallocate (x_gn_int,xi_gn_int,x_fn_int,n_fn_int,xi_i,pphi_i,sphi_i,wqj_i,dme_wqj_i,dme_dxda)
            deallocate (h1,h2,h3,g1,g2,g3,h1p,h2p,h3p,g1p,g2p,g3p,h1m,h2m,h3m,g1m,g2m,g3m)
            deallocate (m1,m2,m3,l1,l2,l3,m1p,m2p,m3p,l1p,l2p,l3p,m1m,m2m,m3m,l1m,l2m,l3m)

          end do ! Loop through the ELEMENTS of the BOUNDARY for INTEGRATION

        end do ! Loop through the BOUNDARIES of the REGION for INTEGRATION

      ! ============================================================================================================================

      ! ============================================================================================================================
      ! FE REGION
      !
      case (fbem_fe)
        stop 'Error: Laplace > FEM > Shape Sensitivity Analysis: not yet'

      ! ============================================================================================================================

    end select ! Switch between REGION CLASSES

  end do ! Loop through REGIONS

  ! Ending message
  if (verbose_level.ge.1) write(output_unit,'(a)') 'done.'

end subroutine build_lse_laplace_sa
