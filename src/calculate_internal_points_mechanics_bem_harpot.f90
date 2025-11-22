! ---------------------------------------------------------------------
! Copyright (C) 2014-2025 Universidad de Las Palmas de Gran Canaria:
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

subroutine calculate_internal_points_mechanics_bem_harpot(kf,kr)

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
  use fbem_bem_general
  use fbem_bem_harpot2d
  use fbem_bem_harpot3d

  ! Module of problem variables
  use problem_variables

  ! No implicit variables
  implicit none

  ! I/O variables
  integer                           :: kf
  integer                           :: kr
  ! Local variables
  real(kind=real64)                 :: omega
  integer                           :: kip, sip, kn, sn
  integer                           :: kb_int, sb_int
  logical                           :: sb_int_reversion
  integer                           :: sp_int
  integer                           :: ke_int, se_int
  integer                           :: se_int_n_nodes
  real(kind=real64)                 :: rd(problem%n), x_i(problem%n)
  ! Region properties
  real(kind=real64)                  :: rho
  complex(kind=real64)               :: c, k
  type(fbem_bem_harpot3d_parameters) :: p3d
  type(fbem_bem_harpot2d_parameters) :: p2d
  ! Writing
  character(len=fbem_fmtstr)        :: fmtstr

  ! Message
  if (verbose_level.ge.2) then
    write(fmtstr,*) '(1x,a6,i',fbem_nchar_int(region(kr)%id),',1x,a27)'
    call fbem_trimall(fmtstr)
    write(output_unit,fmtstr) 'Region', region(kr)%id, '(BE region, inviscid fluid)'
  end if

  ! Frequency
  omega=frequency(kf)

  ! TODO: put this into a function, it also is present in the calculate incident fields and other
  ! calculations (calculate stresses)

  ! Save the region properties to local variables
  rho=region(kr)%property_r(1)
  c=region(kr)%property_c(4)
  ! Wavenumber
  k=omega/c
  ! Calculate the region parameters
  select case (problem%n)
    case (2)
      call fbem_bem_harpot2d_calculate_parameters(rho,c,omega,p2d)
    case (3)
      call fbem_bem_harpot3d_calculate_parameters(rho,c,omega,p3d)
  end select

  ! ============================================= !
  ! OBTAIN DIFFRACTED FIELD INTERNAL POINT VALUES !
  ! ============================================= !

  ! --------------------------------------------
  ! CALCULATE AND ASSEMBLE ELEMENT BEM INTEGRALS
  ! --------------------------------------------

  ! REGION BOUNDARIES
  do kb_int=1,region(kr)%n_boundaries
    sb_int=region(kr)%boundary(kb_int)
    sb_int_reversion=region(kr)%boundary_reversion(kb_int)
    sp_int=boundary(sb_int)%part
    !$omp parallel do schedule (dynamic) default (shared) private (se_int,se_int_n_nodes,fmtstr)
    do ke_int=1,part(sp_int)%n_elements
      se_int=part(boundary(sb_int)%part)%element(ke_int)
      se_int_n_nodes=element(se_int)%n_nodes
      call calculate_internal_points_mechanics_bem_harpot_element(omega,kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,p2d,p3d)
    end do
    !$omp end parallel do
  end do
  ! REGION BODY LOADS
  do kb_int=1,region(kr)%n_be_bodyloads
    sb_int=region(kr)%be_bodyload(kb_int)
    sp_int=be_bodyload(sb_int)%part
    !$omp parallel do schedule (dynamic) default (shared) private (se_int,se_int_n_nodes)
    do ke_int=1,part(sp_int)%n_elements
      se_int=part(sp_int)%element(ke_int)
      se_int_n_nodes=element(se_int)%n_nodes
      call calculate_internal_points_mechanics_bem_harpot_bl(omega,kr,sb_int,se_int,se_int_n_nodes,p2d,p3d)
    end do
    !$omp end parallel do
  end do

  ! ===============================================
  ! CALCULATE THE TOTAL FIELD AT THE INTERNAL POINT
  ! ===============================================

  do kip=1,region(kr)%n_internalpoints
    sip=region(kr)%internalpoint(kip)
    internalpoint(sip)%value_c(:,:)=internalpoint(sip)%value_c(:,:)+internalpoint(sip)%incident_c(:,:)
  end do

  !
  ! Copia los valores nodales a los puntos internos en ellos (solo variable primaria)
  ! se podria hacer con la secundaria (cogiendo gradiente de desplazamientos y construyendo
  ! el tensor de tensiones en el contorno)
  !

  do kip=1,region(kr)%n_internalpoints
    sip=region(kr)%internalpoint(kip)
    x_i=internalpoint(sip)%x
    do kb_int=1,region(kr)%n_boundaries
      sb_int=region(kr)%boundary(kb_int)
      do kn=1,part(boundary(sb_int)%part)%n_nodes
        sn=part(boundary(sb_int)%part)%node(kn)
        rd=x_i-node(sn)%x
        if (sqrt(dot_product(rd,rd))/element(node(sn)%element(1))%csize.lt.1.d-12) then
          internalpoint(sip)%value_c=0
          internalpoint(sip)%value_c(1,0)=node(sn)%value_c(1,1)
        end if
      end do
    end do
  end do

end subroutine calculate_internal_points_mechanics_bem_harpot

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine calculate_internal_points_mechanics_bem_harpot_element(omega,kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,p2d,p3d)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_numerical
  use fbem_symmetry
  use fbem_bem_general
  use fbem_bem_harpot2d
  use fbem_bem_harpot3d

  ! Module of problem variables
  use problem_variables

  ! No implicit variables
  implicit none

  ! I/O variables
  real(kind=real64)                  :: omega
  integer                            :: kr
  integer                            :: sb_int
  logical                            :: sb_int_reversion
  integer                            :: se_int
  integer                            :: se_int_n_nodes
  type(fbem_bem_harpot3d_parameters) :: p3d
  type(fbem_bem_harpot2d_parameters) :: p2d
  ! Local variables
  integer                :: ks
  integer                :: kc
  integer                :: kn_int, sn_int
  type(fbem_bem_element) :: se_int_data
  logical                :: se_int_reversion
  integer                :: kip, sip
  integer                :: kn
  real(kind=real64)      :: x_i(problem%n), n_i(problem%n)
  ! Region properties
  real(kind=real64)      :: d1J
  ! Incident wave variables
  complex(kind=real64)   :: p_inc(se_int_n_nodes), Un_inc(se_int_n_nodes)
  ! Kernels for SBIE integration
  complex(kind=real64)   :: h (se_int_n_nodes), g (se_int_n_nodes)
  complex(kind=real64)   :: hp(se_int_n_nodes), gp(se_int_n_nodes)
  complex(kind=real64)   :: hm(se_int_n_nodes), gm(se_int_n_nodes)
  ! Kernels for HBIE integration
  complex(kind=real64)   :: m (se_int_n_nodes), l (se_int_n_nodes)
  complex(kind=real64)   :: mp(se_int_n_nodes), lp(se_int_n_nodes)
  complex(kind=real64)   :: mm(se_int_n_nodes), lm(se_int_n_nodes)
  ! Associated with symmetry
  real(kind=real64)      :: symconf_m(problem%n), symconf_t(problem%n), symconf_r(problem%n), symconf_s
  logical                :: reversed

  ! 1/J=rho*omega**2
  d1J=region(kr)%property_r(1)*omega**2

  ! Initialize calculation element
  call se_int_data%init
  se_int_data%gtype=element(se_int)%type_g
  se_int_data%d=element(se_int)%n_dimension
  se_int_data%n_gnodes=se_int_n_nodes
  se_int_data%n=problem%n
  allocate (se_int_data%x(problem%n,se_int_n_nodes))
  se_int_data%x=element(se_int)%x_gn
  se_int_data%ptype=element(se_int)%type_f1
  se_int_data%ptype_delta=element(se_int)%delta_f
  se_int_data%n_pnodes=se_int_n_nodes
  se_int_data%stype=element(se_int)%type_f2
  se_int_data%stype_delta=element(se_int)%delta_f
  se_int_data%n_snodes=se_int_n_nodes
  se_int_data%cl=element(se_int)%csize
  se_int_data%gln_far=element(se_int)%n_phi
  allocate (se_int_data%bball_centre(problem%n))
  se_int_data%bball_centre=element(se_int)%bball_centre
  se_int_data%bball_radius=element(se_int)%bball_radius
  ! Copy the incident field over the element
  select case (boundary(sb_int)%coupling)
    case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
      p_inc =element(se_int)%incident_c(1,:,1)
      Un_inc=element(se_int)%incident_c(2,:,1)
    case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
      if (sb_int_reversion) then
        p_inc =element(se_int)%incident_c(1,:,2)
        Un_inc=element(se_int)%incident_c(2,:,2)
      else
        p_inc =element(se_int)%incident_c(1,:,1)
        Un_inc=element(se_int)%incident_c(2,:,1)
      end if
  end select

  !
  ! Loop through symmetrical elements
  !
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
      sip=region(kr)%internalpoint(kip)
      x_i=internalpoint(sip)%x

      ! -------- !
      ! PRESSURE !
      ! -------- !

      ! CALCULATE KERNELS
      select case (problem%n)
        case (2)
          call fbem_bem_harpot2d_sbie_auto(se_int_data,se_int_reversion,x_i,p2d,qsi_parameters,qsi_ns_max,h,g)
        case (3)
          call fbem_bem_harpot3d_sbie_auto(se_int_data,se_int_reversion,x_i,p3d,qsi_parameters,qsi_ns_max,h,g)
      end select
      ! Additional kernels for half-space fundamental solution
      if (region(kr)%space.eq.fbem_half_space) then
        x_i(abs(region(kr)%halfspace_n))=2.d0*region(kr)%halfspace_x-x_i(abs(region(kr)%halfspace_n))
        select case (problem%n)
          case (2)
            call fbem_bem_harpot2d_sbie_auto(se_int_data,se_int_reversion,x_i,p2d,qsi_parameters,qsi_ns_max,hp,gp)
          case (3)
            call fbem_bem_harpot3d_sbie_auto(se_int_data,se_int_reversion,x_i,p3d,qsi_parameters,qsi_ns_max,hp,gp)
        end select
        select case (region(kr)%halfspace_bc)
          ! p=0
          case (0)
            h=h-hp
            g=g-gp
          ! Un=0
          case (1)
            h=h+hp
            g=g+gp
        end select
      end if
      ! BUILD KERNELS ACCORDING TO SYMMETRY
      if (ks.gt.1) then
        h(:)=symconf_s*h(:)
        g(:)=symconf_s*g(:)
      end if
      ! BUILD KERNELS WITH N+ AND N-
      hp=h
      gp=g
      ! If the integration boundary is a crack-like boundary, build N- kernels
      if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
        hm=-h
        gm= g
      end if
      ! ASSEMBLE
      ! Note: The flux variable is the normal displacement Un=1/(rho*omega**2)*dp/dn.
      gp=gp*d1J
      gm=gm*d1J
      !$omp critical
      select case (boundary(sb_int)%coupling)
        ! BE OR BE-FE BOUNDARY
        case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
          select case (boundary(sb_int)%class)
            ! ORDINARY BOUNDARY
            case (fbem_boundary_class_ordinary)
              do kn_int=1,se_int_n_nodes
                sn_int=element(se_int)%node(kn_int)
                internalpoint(sip)%value_c(1,0)=internalpoint(sip)%value_c(1,0)&
                                              +gp(kn_int)*(node(sn_int)%value_c(2,1)-Un_inc(kn_int))&
                                              -hp(kn_int)*(node(sn_int)%value_c(1,1)-p_inc(kn_int))
              end do
            ! CRACK-LIKE BOUNDARY
            case (fbem_boundary_class_cracklike)
              do kn_int=1,se_int_n_nodes
                sn_int=element(se_int)%node(kn_int)
                internalpoint(sip)%value_c(1,0)=internalpoint(sip)%value_c(1,0)&
                                               +gp(kn_int)*(node(sn_int)%value_c(2,1)-Un_inc(kn_int))&
                                               -hp(kn_int)*(node(sn_int)%value_c(1,1)-p_inc(kn_int))&
                                               +gm(kn_int)*(node(sn_int)%value_c(2,2)+Un_inc(kn_int))&
                                               -hm(kn_int)*(node(sn_int)%value_c(1,2)-p_inc(kn_int))
              end do
          end select
        ! BE-BE OR BE-FE-BE BOUNDARY
        case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
          if (sb_int_reversion.eqv.(.false.)) then
            do kn_int=1,se_int_n_nodes
              sn_int=element(se_int)%node(kn_int)
              internalpoint(sip)%value_c(1,0)=internalpoint(sip)%value_c(1,0)&
                                              +gp(kn_int)*(node(sn_int)%value_c(2,1)-Un_inc(kn_int))&
                                              -hp(kn_int)*(node(sn_int)%value_c(1,1)-p_inc(kn_int))
            end do
          else
            do kn_int=1,se_int_n_nodes
              sn_int=element(se_int)%node(kn_int)
              internalpoint(sip)%value_c(1,0)=internalpoint(sip)%value_c(1,0)&
                                             +gp(kn_int)*(node(sn_int)%value_c(2,2)-Un_inc(kn_int))&
                                             -hp(kn_int)*(node(sn_int)%value_c(1,2)-p_inc(kn_int))
            end do
          end if
      end select
      !$omp end critical

      ! =================== !
      ! NORMAL DISPLACEMENT !
      ! =================== !

      do kc=1,problem%n
        ! UNIT NORMAL
        n_i=0.0d0
        n_i(kc)=1.0d0
        ! CALCULATE KERNELS
        select case (problem%n)
          case (2)
            call fbem_bem_harpot2d_hbie_auto(se_int_data,se_int_reversion,x_i,n_i,p2d,qsi_parameters,qsi_ns_max,m,l)
          case (3)
            call fbem_bem_harpot3d_hbie_auto(se_int_data,se_int_reversion,x_i,n_i,p3d,qsi_parameters,qsi_ns_max,m,l)
        end select
        ! Additional kernels for half-space fundamental solution
        if (region(kr)%space.eq.fbem_half_space) then
          x_i(abs(region(kr)%halfspace_n))=2.d0*region(kr)%halfspace_x-x_i(abs(region(kr)%halfspace_n))
          n_i(abs(region(kr)%halfspace_n))=-n_i(abs(region(kr)%halfspace_n))
          select case (problem%n)
            case (2)
              call fbem_bem_harpot2d_hbie_auto(se_int_data,se_int_reversion,x_i,n_i,p2d,qsi_parameters,qsi_ns_max,mp,lp)
            case (3)
              call fbem_bem_harpot3d_hbie_auto(se_int_data,se_int_reversion,x_i,n_i,p3d,qsi_parameters,qsi_ns_max,mp,lp)
          end select
          select case (region(kr)%halfspace_bc)
            ! p=0
            case (0)
              m=m-mp
              l=l-lp
            ! Un=0
            case (1)
              m=m+mp
              l=l+lp
          end select
        end if
        ! BUILD KERNELS ACCORDING TO SYMMETRY
        if (ks.gt.1) then
          m(:)=symconf_s*m(:)
          l(:)=symconf_s*l(:)
        end if
        ! BUILD KERNELS WITH N+ AND N-
        mp=m
        lp=l
        ! If the integration boundary is a crack-like boundary, build N- kernels
        if (boundary(sb_int)%class.eq.fbem_boundary_class_cracklike) then
          mm=-m
          lm= l
        end if
        ! ASSEMBLE
        ! Note: The flux variable is the normal displacement Un=1/(rho*omega**2)*dp/dn.
        mp=mp/d1J
        mm=mm/d1J
        !$omp critical
        select case (boundary(sb_int)%coupling)
          ! BE OR BE-FE BOUNDARY
          case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
            select case (boundary(sb_int)%class)
              ! ORDINARY BOUNDARY
              case (fbem_boundary_class_ordinary)
                do kn_int=1,se_int_n_nodes
                  sn_int=element(se_int)%node(kn_int)
                  internalpoint(sip)%value_c(1,kc)=internalpoint(sip)%value_c(1,kc)&
                                                  +lp(kn_int)*(node(sn_int)%value_c(2,1)-Un_inc(kn_int))&
                                                  -mp(kn_int)*(node(sn_int)%value_c(1,1)-p_inc(kn_int))
                end do
              ! CRACK-LIKE BOUNDARY
              case (fbem_boundary_class_cracklike)
                do kn_int=1,se_int_n_nodes
                  sn_int=element(se_int)%node(kn_int)
                  internalpoint(sip)%value_c(1,kc)=internalpoint(sip)%value_c(1,kc)&
                                                  +lp(kn_int)*(node(sn_int)%value_c(2,1)-Un_inc(kn_int))&
                                                  -mp(kn_int)*(node(sn_int)%value_c(1,1)-p_inc(kn_int))&
                                                  +lm(kn_int)*(node(sn_int)%value_c(2,2)+Un_inc(kn_int))&
                                                  -mm(kn_int)*(node(sn_int)%value_c(1,2)-p_inc(kn_int))
                end do
            end select
          ! BE-BE OR BE-FE-BE BOUNDARY
          case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
            if (sb_int_reversion.eqv.(.false.)) then
              do kn_int=1,se_int_n_nodes
                sn_int=element(se_int)%node(kn_int)
                internalpoint(sip)%value_c(1,kc)=internalpoint(sip)%value_c(1,kc)&
                                                +lp(kn_int)*(node(sn_int)%value_c(2,1)-Un_inc(kn_int))&
                                                -mp(kn_int)*(node(sn_int)%value_c(1,1)-p_inc(kn_int))
              end do
            else
              do kn_int=1,se_int_n_nodes
                sn_int=element(se_int)%node(kn_int)
                internalpoint(sip)%value_c(1,kc)=internalpoint(sip)%value_c(1,kc)&
                                                +lp(kn_int)*(node(sn_int)%value_c(2,2)-Un_inc(kn_int))&
                                                -mp(kn_int)*(node(sn_int)%value_c(1,2)-p_inc(kn_int))
              end do
            end if
        end select
        !$omp end critical
      end do

    end do ! Loop through INTERNAL POINTS

  end do ! Loop through SYMMETRICAL ELEMENTS

end subroutine calculate_internal_points_mechanics_bem_harpot_element

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine calculate_internal_points_mechanics_bem_harpot_bl(omega,kr,sb_int,se_int,se_int_n_nodes,p2d,p3d)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_numerical
  use fbem_symmetry
  use fbem_bem_general
  use fbem_bem_harpot2d
  use fbem_bem_harpot3d

  ! Module of problem variables
  use problem_variables

  ! No implicit variables
  implicit none

  ! I/O variables
  real(kind=real64)                   :: omega
  integer                             :: kr
  integer                             :: sb_int
  integer                             :: se_int
  integer                             :: se_int_n_nodes
  type(fbem_bem_harpot3d_parameters)  :: p3d
  type(fbem_bem_harpot2d_parameters)  :: p2d
  ! Local variables
  integer                :: ks
  integer                :: il, ik, kc
  integer                :: kn_int, sn_int
  type(fbem_bem_element) :: se_int_data
  integer                :: kip, sip
  integer                :: kn
  real(kind=real64)      :: x_i(problem%n), n_i(problem%n)
  ! Kernels for SBIE integration
  complex(kind=real64)   :: g (se_int_n_nodes)
  ! Kernels for HBIE integration
  complex(kind=real64)   :: l (se_int_n_nodes)
  ! Symmetry plane configuration for the current element
  integer                :: se_n_symplanes
  integer                :: se_n_symelements
  real(kind=real64)      :: se_symplane_m(3,3)
  real(kind=real64)      :: se_symplane_s(3)
  real(kind=real64)      :: se_symplane_t(3,3)
  real(kind=real64)      :: se_symplane_r(3,3)
  ! Associated with symmetry
  real(kind=real64)      :: symconf_m(problem%n), symconf_t(problem%n), symconf_r(problem%n), symconf_s
  logical                :: reversed

  ! Initialize calculation element
  call se_int_data%init
  se_int_data%gtype=element(se_int)%type_g
  se_int_data%d=element(se_int)%n_dimension
  se_int_data%n_gnodes=se_int_n_nodes
  se_int_data%n=problem%n
  allocate (se_int_data%x(problem%n,se_int_n_nodes))
  se_int_data%x=element(se_int)%x_gn
  se_int_data%ptype=element(se_int)%type_f1
  se_int_data%ptype_delta=element(se_int)%delta_f
  se_int_data%n_pnodes=se_int_n_nodes
  se_int_data%stype=element(se_int)%type_f2
  se_int_data%stype_delta=element(se_int)%delta_f
  se_int_data%n_snodes=se_int_n_nodes
  se_int_data%cl=element(se_int)%csize
  se_int_data%gln_far=element(se_int)%n_phi
  allocate (se_int_data%bball_centre(problem%n))
  se_int_data%bball_centre=element(se_int)%bball_centre
  se_int_data%bball_radius=element(se_int)%bball_radius
  ! Copy the incident field over the element
  ! not needed..

  ! ACTIVE SYMMETRY PLANES FOR THE CURRENT ELEMENT
  call build_symplane_bodyload_elements(se_int,se_n_symplanes,se_n_symelements,se_symplane_m,se_symplane_s,se_symplane_t,se_symplane_r)

  !
  ! Loop through symmetrical elements
  !
  do ks=1,se_n_symelements
    ! SYMMETRY SETUP
    call fbem_symmetry_multipliers(ks,problem%n,se_n_symplanes,se_symplane_m,se_symplane_s,se_symplane_t,se_symplane_r,&
                                   symconf_m,symconf_s,symconf_t,symconf_r,reversed)
    do kn=1,se_int_n_nodes
      se_int_data%x(:,kn)=symconf_m*element(se_int)%x_gn(:,kn)
    end do
    ! Initialize precalculated datasets
    call se_int_data%init_precalculated_datasets(n_precalsets,precalset_gln)

    !
    ! Loop through INTERNAL POINTS
    !
    do kip=1,region(kr)%n_internalpoints
      sip=region(kr)%internalpoint(kip)
      x_i=internalpoint(sip)%x

      ! ------------ !
      ! DISPLACEMENT !
      ! ------------ !

      ! CALCULATE KERNELS
      select case (problem%n)
        case (2)
          call fbem_bem_harpot2d_sbie_bl_auto(se_int_data,x_i,p2d,qsi_parameters,qsi_ns_max,g)
        case (3)
          call fbem_bem_harpot3d_sbie_bl_auto(se_int_data,x_i,p3d,qsi_parameters,qsi_ns_max,g)
      end select
      ! BUILD KERNELS ACCORDING TO SYMMETRY
      if (ks.gt.1) then
        g(:)=symconf_s*g(:)
      end if
      ! ASSEMBLE
      !$omp critical
      do kn_int=1,se_int_n_nodes
        sn_int=element(se_int)%node(kn_int)
        internalpoint(sip)%value_c(1,0)=internalpoint(sip)%value_c(1,0)+g(kn_int)*node(sn_int)%value_c(1,1)
      end do
      !$omp end critical

      ! -------- !
      ! TRACTION !
      ! -------- !

      do kc=1,problem%n
        ! UNIT NORMAL
        n_i=0.0d0
        n_i(kc)=1.0d0
        ! CALCULATE KERNELS
        select case (problem%n)
          case (2)
            call fbem_bem_harpot2d_hbie_bl_auto(se_int_data,x_i,n_i,p2d,qsi_parameters,qsi_ns_max,l)
          case (3)
            call fbem_bem_harpot3d_hbie_bl_auto(se_int_data,x_i,n_i,p3d,qsi_parameters,qsi_ns_max,l)
        end select
        ! BUILD KERNELS ACCORDING TO SYMMETRY
        if (ks.gt.1) then
          l(:)=symconf_s*l(:)
        end if
        ! ASSEMBLE
        !$omp critical
        do kn_int=1,se_int_n_nodes
          sn_int=element(se_int)%node(kn_int)
          internalpoint(sip)%value_c(1,kc)=internalpoint(sip)%value_c(1,kc)+l(kn_int)*node(sn_int)%value_c(1,1)
        end do
        !$omp end critical
      end do

    end do ! Loop through INTERNAL POINTS

  end do ! Loop through SYMMETRICAL ELEMENTS

end subroutine calculate_internal_points_mechanics_bem_harpot_bl
