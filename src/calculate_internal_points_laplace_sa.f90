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


subroutine calculate_internal_points_laplace_sa

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
  logical                        :: se_int_sensitivity
  real(kind=real64), allocatable :: x_gn_int(:,:), xi_gn_int(:,:)
  real(kind=real64), allocatable :: x_fn_int(:,:), n_fn_int(:,:)
  integer                        :: kn_int, sn_int
  integer                        :: kip, sip
  integer                        :: kn, kc, ka, im, iq
  real(kind=real64)              :: x_i(problem%n), n_i(problem%n)
  ! Region properties
  real(kind=real64)              :: conductivity
  ! Vectors for VSBIE integration
  real(kind=real64), allocatable :: h(:), g(:)
  real(kind=real64), allocatable :: hp(:), gp(:)
  real(kind=real64), allocatable :: hm(:), gm(:)
  real(kind=real64), allocatable :: h1(:,:,:), h2(:,:), g1(:,:,:), g2(:,:)
  real(kind=real64), allocatable :: h1p(:,:,:), h2p(:,:), g1p(:,:,:), g2p(:,:)
  real(kind=real64), allocatable :: h1m(:,:,:), h2m(:,:), g1m(:,:,:), g2m(:,:)
  real(kind=real64), allocatable :: hsap(:), gsap(:)
  real(kind=real64), allocatable :: hsam(:), gsam(:)
  ! Vectors for VHBIE integration
  real(kind=real64), allocatable :: m(:), l(:)
  real(kind=real64), allocatable :: mp(:), lp(:)
  real(kind=real64), allocatable :: mm(:), lm(:)
  real(kind=real64), allocatable :: m1(:,:,:), m2(:,:), m3(:,:), l1(:,:,:), l2(:,:), l3(:,:)
  real(kind=real64), allocatable :: m1p(:,:,:), m2p(:,:), m3p(:,:), l1p(:,:,:), l2p(:,:), l3p(:,:)
  real(kind=real64), allocatable :: m1m(:,:,:), m2m(:,:), m3m(:,:), l1m(:,:,:), l2m(:,:), l3m(:,:)
  real(kind=real64), allocatable :: msap(:), lsap(:)
  real(kind=real64), allocatable :: msam(:), lsam(:)
  ! Associated with the DME
  integer                        :: sdme_int
  integer                        :: sdme_int_n_nodes
  real(kind=real64), allocatable :: dme_dxda(:,:,:)
  real(kind=real64)              :: dxda_i(problem%n,problem%n_designvariables)
  real(kind=real64)              :: dnda_i(problem%n,problem%n_designvariables)
  ! Associated with symmetry
  real(kind=real64), allocatable :: symconf_m(:), symconf_t(:), symconf_r(:)
  real(kind=real64)              :: symconf_s
  logical                        :: reversed
  ! Writing
  character(len=fbem_fmtstr)     :: fmtstr

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'START calculating internal points sensitivity solutions at BE regions')

  ! Initialize internal points values
  do kip=1,n_internalpoints
    internalpoint(kip)%dvda_r=0.d0
  end do

  ! Allocate auxiliary variables
  allocate (symconf_m(problem%n))
  allocate (symconf_t(problem%n))
  allocate (symconf_r(problem%n))

  ! Loop through REGIONS
  do kr=1,n_regions

    ! Save the conductivity
    conductivity=region(kr)%property_r(1)

    ! Message
    if (verbose_level.ge.2) then
      write(fmtstr,*) '(1x,a,i',fbem_nchar_int(region(kr)%id),')'
      call fbem_trimall(fmtstr)
      write(output_unit,fmtstr) 'Region: ', region(kr)%id
    end if

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
        ! ???¿mirar que se puede ahorrar usando los marcadores de sensibilidad
        se_int_sensitivity=element(se_int)%sensitivity
        ! ???
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
        allocate (h(se_int_n_nodes),g(se_int_n_nodes))
        allocate (hp(se_int_n_nodes),gp(se_int_n_nodes))
        allocate (hm(se_int_n_nodes),gm(se_int_n_nodes))
        allocate (h1(sdme_int_n_nodes,se_int_n_nodes,problem%n),h2(se_int_n_nodes,problem%n))
        allocate (g1(sdme_int_n_nodes,se_int_n_nodes,problem%n),g2(se_int_n_nodes,problem%n))
        allocate (h1p(sdme_int_n_nodes,se_int_n_nodes,problem%n),h2p(se_int_n_nodes,problem%n))
        allocate (g1p(sdme_int_n_nodes,se_int_n_nodes,problem%n),g2p(se_int_n_nodes,problem%n))
        allocate (h1m(sdme_int_n_nodes,se_int_n_nodes,problem%n),h2m(se_int_n_nodes,problem%n))
        allocate (g1m(sdme_int_n_nodes,se_int_n_nodes,problem%n),g2m(se_int_n_nodes,problem%n))
        allocate (hsap(se_int_n_nodes),gsap(se_int_n_nodes))
        allocate (hsam(se_int_n_nodes),gsam(se_int_n_nodes))
        allocate (m(se_int_n_nodes),l(se_int_n_nodes))
        allocate (mp(se_int_n_nodes),lp(se_int_n_nodes))
        allocate (mm(se_int_n_nodes),lm(se_int_n_nodes))
        allocate (m1(sdme_int_n_nodes,se_int_n_nodes,problem%n),m2(se_int_n_nodes,problem%n),m3(se_int_n_nodes,problem%n))
        allocate (l1(sdme_int_n_nodes,se_int_n_nodes,problem%n),l2(se_int_n_nodes,problem%n),l3(se_int_n_nodes,problem%n))
        allocate (m1p(sdme_int_n_nodes,se_int_n_nodes,problem%n),m2p(se_int_n_nodes,problem%n),m3p(se_int_n_nodes,problem%n))
        allocate (l1p(sdme_int_n_nodes,se_int_n_nodes,problem%n),l2p(se_int_n_nodes,problem%n),l3p(se_int_n_nodes,problem%n))
        allocate (m1m(sdme_int_n_nodes,se_int_n_nodes,problem%n),m2m(se_int_n_nodes,problem%n),m3m(se_int_n_nodes,problem%n))
        allocate (l1m(sdme_int_n_nodes,se_int_n_nodes,problem%n),l2m(se_int_n_nodes,problem%n),l3m(se_int_n_nodes,problem%n))
        allocate (msap(se_int_n_nodes),lsap(se_int_n_nodes))
        allocate (msam(se_int_n_nodes),lsam(se_int_n_nodes))
        allocate (dme_dxda(problem%n,sdme_int_n_nodes,problem%n_designvariables))
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

          !
          ! Loop through INTERNAL POINTS
          !
          do kip=1,region(kr)%n_internalpoints
            ! Internal point iid
            sip=region(kr)%internalpoint(kip)
            ! Save x_i, dxda
            x_i=internalpoint(sip)%x
            dxda_i=internalpoint(sip)%dxda

            ! --------------
            ! CALCULATE DPDA
            ! --------------

            ! CALCULATE KERNELS
            select case (problem%n)
              case (2)
                call fbem_bem_stapot2d_sbie_auto(se_int_data,se_int_reversion,x_i,qsi_parameters,qsi_ns_max,h,g)
                call fbem_bem_stapot2d_vsbie_auto(se_int_data,se_int_reversion,x_i,qsi_parameters,qsi_ns_max,h1,h2,g1,g2)
              case (3)
                call fbem_bem_stapot3d_sbie_auto(se_int_data,se_int_reversion,x_i,qsi_parameters,qsi_ns_max,h,g)
                stop 'stapot3d_vsbie not yet'
            end select
            ! BUILD KERNELS ACCORDING TO SYMMETRY
            if (ks.gt.1) then
              h=symconf_s*h
              g=symconf_s*g
              stop 'symmetry not yet for sensibility'
              h1=symconf_s*h1
              h2=symconf_s*h2
              g1=symconf_s*g1
              g2=symconf_s*g2
            end if
            ! BUILD KERNELS WITH N+ AND N-
            hp =h
            gp =g
            h1p=h1
            h2p=h2
            g1p=g1
            g2p=g2
            if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
              hm =-h
              gm = g
              h1m=-h1
              h2m=-h2
              g1m= g1
              g2m= g2
            end if
            ! ASSEMBLE
            ! The flux variable is j=k·dp/dn
            gp=gp/conductivity
            gm=gm/conductivity
            g1p=g1p/conductivity
            g2p=g2p/conductivity
            g1m=g1m/conductivity
            g2m=g2m/conductivity
            ! Loop through DESIGN VARIABLES
            do ka=1,problem%n_designvariables
              ! Build h^{sa} and g^{sa} taken into account the design velocity field
              hsap=0.d0
              gsap=0.d0
              do im=1,problem%n
                do iq=1,sdme_int_n_nodes
                  hsap(:)=hsap(:)+h1p(iq,:,im)*dme_dxda(im,iq,ka)
                  gsap(:)=gsap(:)+g1p(iq,:,im)*dme_dxda(im,iq,ka)
                end do
                hsap(:)=hsap(:)-h2p(:,im)*dxda_i(im,ka)
                gsap(:)=gsap(:)-g2p(:,im)*dxda_i(im,ka)
              end do
              if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
                hsam=0.d0
                gsam=0.d0
                do im=1,problem%n
                  do iq=1,sdme_int_n_nodes
                    hsam(:)=hsam(:)+h1m(iq,:,im)*dme_dxda(im,iq,ka)
                    gsam(:)=gsam(:)+g1m(iq,:,im)*dme_dxda(im,iq,ka)
                  end do
                  hsam(:)=hsam(:)-h2m(:,im)*dxda_i(im,ka)
                  gsam(:)=gsam(:)-g2m(:,im)*dxda_i(im,ka)
                end do
              end if
              ! Assemble
              select case (boundary(sb_int)%coupling)
                ! BE OR BE-FE BOUNDARY
                case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
                  select case (boundary(sb_int)%class)
                    ! ORDINARY BOUNDARY
                    case (fbem_boundary_class_ordinary)
                      do kn_int=1,se_int_n_nodes
                        sn_int=element(se_int)%node(kn_int)
                        internalpoint(sip)%dvda_r(1,0,ka)=internalpoint(sip)%dvda_r(1,0,ka)&
                                                         +gp(kn_int)*node(sn_int)%dvda_r(2,1,ka)-hp(kn_int)*node(sn_int)%dvda_r(1,1,ka)&
                                                         +gsap(kn_int)*node(sn_int)%value_r(2,1)-hsap(kn_int)*node(sn_int)%value_r(1,1)
                      end do
                    ! CRACK-LIKE BOUNDARY
                    case (fbem_boundary_class_cracklike)
                      do kn_int=1,se_int_n_nodes
                        sn_int=element(se_int)%node(kn_int)
                        internalpoint(sip)%dvda_r(1,0,ka)=internalpoint(sip)%dvda_r(1,0,ka)&
                                                         +gp(kn_int)*node(sn_int)%dvda_r(2,1,ka)-hp(kn_int)*node(sn_int)%dvda_r(1,1,ka)&
                                                         +gsap(kn_int)*node(sn_int)%value_r(2,1)-hsap(kn_int)*node(sn_int)%value_r(1,1)&
                                                         +gm(kn_int)*node(sn_int)%dvda_r(2,2,ka)-hm(kn_int)*node(sn_int)%dvda_r(1,2,ka)&
                                                         +gsam(kn_int)*node(sn_int)%value_r(2,2)-hsam(kn_int)*node(sn_int)%value_r(1,2)
                      end do
                  end select
                ! BE-BE OR BE-FE-BE BOUNDARY
                case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                  if (sb_int_reversion.eqv.(.false.)) then
                    do kn_int=1,se_int_n_nodes
                      sn_int=element(se_int)%node(kn_int)
                      internalpoint(sip)%dvda_r(1,0,ka)=internalpoint(sip)%dvda_r(1,0,ka)&
                                                       +gp(kn_int)*node(sn_int)%dvda_r(2,1,ka)-hp(kn_int)*node(sn_int)%dvda_r(1,1,ka)&
                                                       +gsap(kn_int)*node(sn_int)%value_r(2,1)-hsap(kn_int)*node(sn_int)%value_r(1,1)
                    end do
                  else
                    do kn_int=1,se_int_n_nodes
                      sn_int=element(se_int)%node(kn_int)
                      internalpoint(sip)%dvda_r(1,0,ka)=internalpoint(sip)%dvda_r(1,0,ka)&
                                                       +gp(kn_int)*node(sn_int)%dvda_r(2,2,ka)-hp(kn_int)*node(sn_int)%dvda_r(1,2,ka)&
                                                       +gsap(kn_int)*node(sn_int)%value_r(2,2)-hsap(kn_int)*node(sn_int)%value_r(1,2)
                    end do
                  end if
              end select
            end do ! Loop through DESIGN VARIABLES

            ! ---------------------------
            ! CALCULATE DJDA WITH N=E_{K}
            ! ---------------------------

            ! Loop through AXES
            do kc=1,problem%n
              ! UNIT NORMAL
              n_i=0.d0
              n_i(kc)=1.d0
              ! DNDA
              dnda_i=internalpoint(sip)%dnda(:,kc,:)
              ! CALCULATE KERNELS
              select case (problem%n)
                case (2)
                  call fbem_bem_stapot2d_hbie_auto(se_int_data,se_int_reversion,x_i,n_i,qsi_parameters,qsi_ns_max,m,l)
                  call fbem_bem_stapot2d_vhbie_auto(se_int_data,se_int_reversion,x_i,n_i,qsi_parameters,qsi_ns_max,m1,m2,m3,l1,l2,l3)
                case (3)
                  call fbem_bem_stapot3d_hbie_auto(se_int_data,se_int_reversion,x_i,n_i,qsi_parameters,qsi_ns_max,m,l)
                  stop 'stapot3d_vhbie not yet'
              end select
              ! BUILD KERNELS ACCORDING TO SYMMETRY
              if (ks.gt.1) then
                m=symconf_s*m
                l=symconf_s*l
                stop 'symmetry not yet for sensibility'
                m1=symconf_s*m1
                m2=symconf_s*m2
                m3=symconf_s*m3
                l1=symconf_s*l1
                l2=symconf_s*l2
                l3=symconf_s*l3
              end if
              ! BUILD KERNELS WITH N+ AND N-
              mp =m
              lp =l
              m1p=m1
              m2p=m2
              m3p=m3
              l1p=l1
              l2p=l2
              l3p=l3
              if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
                mm =-m
                lm = l
                m1m=-m1
                m2m=-m2
                m3m=-m3
                l1m= l1
                l2m= l2
                l3m= l3
              end if
              ! ASSEMBLE
              ! The flux variable is j=k·dp/dn
              mp=mp*conductivity
              mm=mm*conductivity
              m1p=m1p*conductivity
              m2p=m2p*conductivity
              m3p=m3p*conductivity
              m1m=m1m*conductivity
              m2m=m2m*conductivity
              m3m=m3m*conductivity
              ! Loop through DESIGN VARIABLES
              do ka=1,problem%n_designvariables
                ! Build m^{sa} and l^{sa} taken into account the design velocity field
                msap=0.d0
                lsap=0.d0
                do im=1,problem%n
                  do iq=1,sdme_int_n_nodes
                    msap(:)=msap(:)+m1p(iq,:,im)*dme_dxda(im,iq,ka)
                    lsap(:)=lsap(:)+l1p(iq,:,im)*dme_dxda(im,iq,ka)
                  end do
                  msap(:)=msap(:)-m2p(:,im)*dxda_i(im,ka)
                  lsap(:)=lsap(:)-l2p(:,im)*dxda_i(im,ka)
                  msap(:)=msap(:)+m3p(:,im)*dnda_i(im,ka)
                  lsap(:)=lsap(:)+l3p(:,im)*dnda_i(im,ka)
                end do
                if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
                  msam=0.d0
                  lsam=0.d0
                  do im=1,problem%n
                    do iq=1,sdme_int_n_nodes
                      msam(:)=msam(:)+m1m(iq,:,im)*dme_dxda(im,iq,ka)
                      lsam(:)=lsam(:)+l1m(iq,:,im)*dme_dxda(im,iq,ka)
                    end do
                    msam(:)=msam(:)-m2m(:,im)*dxda_i(im,ka)
                    lsam(:)=lsam(:)-l2m(:,im)*dxda_i(im,ka)
                    msam(:)=msam(:)+m3m(:,im)*dnda_i(im,ka)
                    lsam(:)=lsam(:)+l3m(:,im)*dnda_i(im,ka)
                  end do
                end if
                ! Assemble
                select case (boundary(sb_int)%coupling)
                  ! BE OR BE-FE BOUNDARY
                  case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
                    select case (boundary(sb_int)%class)
                      ! ORDINARY BOUNDARY
                      case (fbem_boundary_class_ordinary)
                        do kn_int=1,se_int_n_nodes
                          sn_int=element(se_int)%node(kn_int)
                          internalpoint(sip)%dvda_r(1,kc,ka)=internalpoint(sip)%dvda_r(1,kc,ka)&
                                                            +lp(kn_int)*node(sn_int)%dvda_r(2,1,ka)-mp(kn_int)*node(sn_int)%dvda_r(1,1,ka)&
                                                            +lsap(kn_int)*node(sn_int)%value_r(2,1)-msap(kn_int)*node(sn_int)%value_r(1,1)
                        end do
                      ! CRACK-LIKE BOUNDARY
                      case (fbem_boundary_class_cracklike)
                        do kn_int=1,se_int_n_nodes
                          sn_int=element(se_int)%node(kn_int)
                          internalpoint(sip)%dvda_r(1,kc,ka)=internalpoint(sip)%dvda_r(1,kc,ka)&
                                                            +lp(kn_int)*node(sn_int)%dvda_r(2,1,ka)-mp(kn_int)*node(sn_int)%dvda_r(1,1,ka)&
                                                            +lsap(kn_int)*node(sn_int)%value_r(2,1)-msap(kn_int)*node(sn_int)%value_r(1,1)&
                                                            +lm(kn_int)*node(sn_int)%dvda_r(2,2,ka)-mm(kn_int)*node(sn_int)%dvda_r(1,2,ka)&
                                                            +lsam(kn_int)*node(sn_int)%value_r(2,2)-msam(kn_int)*node(sn_int)%value_r(1,2)
                        end do
                    end select
                  ! BE-BE OR BE-FE-BE BOUNDARY
                  case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                    if (sb_int_reversion.eqv.(.false.)) then
                      do kn_int=1,se_int_n_nodes
                        sn_int=element(se_int)%node(kn_int)
                        internalpoint(sip)%dvda_r(1,kc,ka)=internalpoint(sip)%dvda_r(1,kc,ka)&
                                                          +lp(kn_int)*node(sn_int)%dvda_r(2,1,ka)-mp(kn_int)*node(sn_int)%dvda_r(1,1,ka)&
                                                          +lsap(kn_int)*node(sn_int)%value_r(2,1)-msap(kn_int)*node(sn_int)%value_r(1,1)
                      end do
                    else
                      do kn_int=1,se_int_n_nodes
                        sn_int=element(se_int)%node(kn_int)
                        internalpoint(sip)%dvda_r(1,kc,ka)=internalpoint(sip)%dvda_r(1,kc,ka)&
                                                          +lp(kn_int)*node(sn_int)%dvda_r(2,2,ka)-mp(kn_int)*node(sn_int)%dvda_r(1,2,ka)&
                                                          +lsap(kn_int)*node(sn_int)%value_r(2,2)-msap(kn_int)*node(sn_int)%value_r(1,2)
                      end do
                    end if
                end select

              end do ! Loop through DESIGN VARIABLES

            end do ! Loop through AXES

          end do ! Loop through INTERNAL POINTS

        end do ! Loop through SYMMETRICAL ELEMENTS

        ! Deallocate element-wise data structures
        deallocate (x_gn_int,xi_gn_int,x_fn_int,n_fn_int,dme_dxda)
        deallocate (h,g,hp,gp,hm,gm)
        deallocate (m,l,mp,lp,mm,lm)
        deallocate (h1,h2,g1,g2,h1p,h2p,g1p,g2p,h1m,h2m,g1m,g2m)
        deallocate (m1,m2,m3,l1,l2,l3,m1p,m2p,m3p,l1p,l2p,l3p,m1m,m2m,m3m,l1m,l2m,l3m)
        deallocate (hsap,gsap,hsam,gsam)
        deallocate (msap,lsap,msam,lsam)

      end do ! Loop through the ELEMENTS of the BOUNDARY for COLLOCATION

    end do ! Loop through the BOUNDARIES of the REGION for INTEGRATION

  end do ! Loop through REGIONS

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'END calculating internal points sensitivity solutions at BE regions')

end subroutine calculate_internal_points_laplace_sa
