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

subroutine calculate_internal_points_mechanics_static

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_shape_functions

  ! Module of problem variables
  use problem_variables

  ! No implicit variables
  implicit none

  ! Local variables
  integer                    :: kr, kp
  integer                    :: kip
  integer                    :: kb_int, sb_int
  logical                    :: sb_int_reversion
  integer                    :: sp_int
  integer                    :: ke_int, se_int
  integer                    :: se_int_n_nodes
  integer                    :: etype, ke, se, i, j
  real(kind=real64)          :: delta_f, aux(10), xi1d
  real(kind=real64), allocatable ::xi(:), phi(:), value_r(:,:,:)

  ! Region properties
  real(kind=real64)          :: mu, nu
  ! Writing
  character(len=fbem_fmtstr) :: fmtstr            ! String used for write format string

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START calculating internal points solutions at BE regions')

  ! Initialize internal points
  do kip=1,n_internalpoints
    internalpoint(kip)%value_r=0
  end do

  ! REGIONS
  do kr=1,n_regions
    if ((region(kr)%class.eq.fbem_be).and.(region(kr)%n_internalpoints.gt.0)) then

      ! Copy material properties
      mu=region(kr)%property_r(2)
      nu=region(kr)%property_r(3)
      if ((problem%n.eq.2).and.(problem%subtype.eq.fbem_mechanics_plane_stress)) nu=nu/(1.d0+nu)

      ! REGION BOUNDARIES
      do kb_int=1,region(kr)%n_boundaries
        sb_int=region(kr)%boundary(kb_int)
        sb_int_reversion=region(kr)%boundary_reversion(kb_int)
        sp_int=boundary(sb_int)%part
        !$omp parallel do schedule (dynamic) default (shared) private (se_int,se_int_n_nodes)
        do ke_int=1,part(sp_int)%n_elements
          se_int=part(boundary(sb_int)%part)%element(ke_int)
          se_int_n_nodes=element(se_int)%n_nodes
          call calculate_internal_points_mechanics_bem_staela_element(kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,mu,nu)
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
          call calculate_internal_points_mechanics_bem_staela_bl(kr,sb_int,se_int,se_int_n_nodes,mu,nu)
        end do
        !$omp end parallel do
      end do

    end if
  end do

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END calculating internal points solutions at BE regions')

  ! ======================================================================
  ! TRANSFER INTERNAL POINT SOLUTIONS TO INTERNAL ELEMENTS INTERNAL POINTS
  ! ======================================================================

  if (internalelements) then

    if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'START mapping the internal points corresponding to internal elements')

    do kp=1,internalelements_mesh%n_parts
      kr=internalelements_mesh%part(kp)%entity
      if (kr.eq.0) cycle
      if (region(kr)%class.eq.fbem_be) then

        do ke=1,internalelements_mesh%part(kp)%n_elements
          se=internalelements_mesh%part(kp)%element(ke)
          etype=internalelements_mesh%element(se)%type
          delta_f=internalelements_mesh%element(se)%delta_f
          internalelements_mesh%element(se)%value_r=0
          allocate (xi(internalelements_mesh%element(se)%n_dimension))
          allocate (phi(internalelements_mesh%element(se)%n_nodes))
          allocate (value_r(1:problem%n,internalelements_mesh%element(se)%n_nodes,0:problem%n))
          ! Copy the corresponding internal point solutions to value_r
          value_r=0
          do j=1,internalelements_mesh%element(se)%n_nodes
            kip=internalelements_mesh%element(se)%internalpoint(j)
            value_r(:,j,:)=internalpoint(kip)%value_r
          end do
          ! Calculate values at element nodes
          select case (internalelements_mesh%element(se)%n_dimension)
            case (1)
              do j=1,internalelements_mesh%element(se)%n_nodes
#               define node j
#               define delta 0.d0
#               define xi xi1d
#               include <xi_1d_at_node.rc>
#               undef delta
#               undef node
#               define delta delta_f
#               include <phi_1d.rc>
#               undef xi
#               undef delta
                do i=1,internalelements_mesh%element(se)%n_nodes
                  internalelements_mesh%element(se)%value_r(:,j,:)=internalelements_mesh%element(se)%value_r(:,j,:)+phi(i)*value_r(:,j,:)
                end do
              end do
            case (2)
              do j=1,internalelements_mesh%element(se)%n_nodes
#               define node j
#               define delta 0.d0
#               include <xi_2d_at_node.rc>
#               undef delta
#               undef node
#               define delta delta_f
#               include <phi_2d.rc>
#               undef delta
                do i=1,internalelements_mesh%element(se)%n_nodes
                  internalelements_mesh%element(se)%value_r(:,j,:)=internalelements_mesh%element(se)%value_r(:,j,:)+phi(i)*value_r(:,j,:)
                end do
              end do
            case (3)
              do j=1,internalelements_mesh%element(se)%n_nodes
#               define node j
#               define delta 0.d0
#               include <xi_3d_at_node.rc>
#               undef delta
#               undef node
#               define delta delta_f
#               include <phi_3d.rc>
#               undef delta
                do i=1,internalelements_mesh%element(se)%n_nodes
                  internalelements_mesh%element(se)%value_r(:,j,:)=internalelements_mesh%element(se)%value_r(:,j,:)+phi(i)*value_r(:,j,:)
                end do
              end do
          end select
          deallocate(xi,phi,value_r)
        end do
      end if
    end do

    if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'END mapping the internal points corresponding to internal elements')

  end if

end subroutine calculate_internal_points_mechanics_static

subroutine calculate_internal_points_mechanics_bem_staela_element(kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,mu,nu)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_numerical
  use fbem_symmetry
  use fbem_bem_general
  use fbem_bem_staela2d
  use fbem_bem_staela3d

  ! Module of problem variables
  use problem_variables

  ! No implicit variables
  implicit none

  ! I/O variables
  integer           :: kr
  integer           :: sb_int
  logical           :: sb_int_reversion
  integer           :: se_int
  integer           :: se_int_n_nodes
  real(kind=real64) :: mu
  real(kind=real64) :: nu
  ! Local variables
  integer                :: ks
  integer                :: il, ik, kc
  integer                :: kn_int, sn_int
  type(fbem_bem_element) :: se_int_data
  logical                :: se_int_reversion
  integer                :: kip, sip
  integer                :: kn
  real(kind=real64)      :: x_i(problem%n), n_i(problem%n)
  ! Kernels for SBIE integration
  real(kind=real64)      :: h (se_int_n_nodes,problem%n,problem%n), g (se_int_n_nodes,problem%n,problem%n)
  real(kind=real64)      :: hp(se_int_n_nodes,problem%n,problem%n), gp(se_int_n_nodes,problem%n,problem%n)
  real(kind=real64)      :: hm(se_int_n_nodes,problem%n,problem%n), gm(se_int_n_nodes,problem%n,problem%n)
  ! Kernels for HBIE integration
  real(kind=real64)      :: m (se_int_n_nodes,problem%n,problem%n), l (se_int_n_nodes,problem%n,problem%n)
  real(kind=real64)      :: mp(se_int_n_nodes,problem%n,problem%n), lp(se_int_n_nodes,problem%n,problem%n)
  real(kind=real64)      :: mm(se_int_n_nodes,problem%n,problem%n), lm(se_int_n_nodes,problem%n,problem%n)
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

      ! ------------ !
      ! DISPLACEMENT !
      ! ------------ !

      ! CALCULATE KERNELS
      select case (problem%n)
        case (2)
          call fbem_bem_staela2d_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,qsi_parameters,qsi_ns_max,h,g)
        case (3)
          call fbem_bem_staela3d_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,qsi_parameters,qsi_ns_max,h,g)
      end select
      ! Additional kernels for half-space fundamental solution
      if (region(kr)%space.eq.fbem_half_space) then
        select case (problem%n)
          case (2)
            x_i(abs(region(kr)%halfspace_n))=2.d0*region(kr)%halfspace_x-x_i(abs(region(kr)%halfspace_n))
            n_i(abs(region(kr)%halfspace_n))=-n_i(abs(region(kr)%halfspace_n))
            call fbem_bem_staela2d_hfc_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,region(kr)%halfspace_n,region(kr)%halfspace_x,qsi_parameters,qsi_ns_max,hp,gp)
          case (3)
            call fbem_bem_staela3d_hsc_sbie_auto(se_int_data,se_int_reversion,x_i,mu,nu,region(kr)%halfspace_n,region(kr)%halfspace_x,qsi_parameters,qsi_ns_max,hp,gp)
        end select
        select case (region(kr)%halfspace_bc)
          ! u_k=0
          case (0)
            stop 'half-space with u_k=0 not available'
          ! t_k=0
          case (1)
            h=h+hp
            g=g+gp
        end select
      end if
      ! BUILD KERNELS ACCORDING TO SYMMETRY
      if (ks.gt.1) then
        do ik=1,problem%n
          h(:,:,ik)=symconf_t(ik)*h(:,:,ik)
          g(:,:,ik)=symconf_t(ik)*g(:,:,ik)
        end do
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
      !$omp critical
      select case (boundary(sb_int)%coupling)
        ! BE OR BE-FE BOUNDARY
        case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
          select case (boundary(sb_int)%class)
            ! ORDINARY BOUNDARY
            case (fbem_boundary_class_ordinary)
              do il=1,problem%n
                do ik=1,problem%n
                  do kn_int=1,se_int_n_nodes
                    sn_int=element(se_int)%node(kn_int)
                    internalpoint(sip)%value_r(il,0)=internalpoint(sip)%value_r(il,0)&
                                                    +gp(kn_int,il,ik)*node(sn_int)%value_r(problem%n+ik,1)&
                                                    -hp(kn_int,il,ik)*node(sn_int)%value_r(          ik,1)
                  end do
                end do
              end do
            ! CRACK-LIKE BOUNDARY
            case (fbem_boundary_class_cracklike)
              do il=1,problem%n
                do ik=1,problem%n
                  do kn_int=1,se_int_n_nodes
                    sn_int=element(se_int)%node(kn_int)
                    internalpoint(sip)%value_r(il,0)=internalpoint(sip)%value_r(il,0)&
                                                    +gp(kn_int,il,ik)*node(sn_int)%value_r(problem%n+ik,1)&
                                                    -hp(kn_int,il,ik)*node(sn_int)%value_r(          ik,1)&
                                                    +gm(kn_int,il,ik)*node(sn_int)%value_r(problem%n+ik,2)&
                                                    -hm(kn_int,il,ik)*node(sn_int)%value_r(          ik,2)
                  end do
                end do
              end do
          end select
        ! BE-BE OR BE-FE-BE BOUNDARY
        case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
          if (sb_int_reversion.eqv.(.false.)) then
            do il=1,problem%n
              do ik=1,problem%n
                do kn_int=1,se_int_n_nodes
                  sn_int=element(se_int)%node(kn_int)
                  internalpoint(sip)%value_r(il,0)=internalpoint(sip)%value_r(il,0)&
                                                  +gp(kn_int,il,ik)*node(sn_int)%value_r(problem%n+ik,1)&
                                                  -hp(kn_int,il,ik)*node(sn_int)%value_r(          ik,1)
                end do
              end do
            end do
          else
            do il=1,problem%n
              do ik=1,problem%n
                do kn_int=1,se_int_n_nodes
                  sn_int=element(se_int)%node(kn_int)
                  internalpoint(sip)%value_r(il,0)=internalpoint(sip)%value_r(il,0)&
                                                  +gp(kn_int,il,ik)*node(sn_int)%value_r(problem%n+ik,2)&
                                                  -hp(kn_int,il,ik)*node(sn_int)%value_r(          ik,2)
                end do
              end do
            end do
          end if
      end select
      !$omp end critical

      ! ======== !
      ! TRACTION !
      ! ======== !

      do kc=1,problem%n
        ! UNIT NORMAL
        n_i=0.0d0
        n_i(kc)=1.0d0
        ! CALCULATE KERNELS
        select case (problem%n)
          case (2)
            call fbem_bem_staela2d_hbie_auto(se_int_data,se_int_reversion,x_i,n_i,mu,nu,qsi_parameters,qsi_ns_max,m,l)
          case (3)
            call fbem_bem_staela3d_hbie_auto(se_int_data,se_int_reversion,x_i,n_i,mu,nu,qsi_parameters,qsi_ns_max,m,l)
        end select
        ! Additional kernels for half-space fundamental solution
        if (region(kr)%space.eq.fbem_half_space) then
          select case (problem%n)
            case (2)
              x_i(abs(region(kr)%halfspace_n))=2.d0*region(kr)%halfspace_x-x_i(abs(region(kr)%halfspace_n))
              n_i(abs(region(kr)%halfspace_n))=-n_i(abs(region(kr)%halfspace_n))
              stop 'elastostatics half-plane not yet'
            case (3)
              call fbem_bem_staela3d_hsc_hbie_auto(se_int_data,se_int_reversion,x_i,n_i,mu,nu,region(kr)%halfspace_n,region(kr)%halfspace_x,qsi_parameters,qsi_ns_max,mp,lp)
          end select
          select case (region(kr)%halfspace_bc)
            ! u_k=0
            case (0)
              stop 'half-space with u_k=0 not available'
            ! t_k=0
            case (1)
              m=m+mp
              l=l+lp
          end select
        end if
        ! BUILD KERNELS ACCORDING TO SYMMETRY
        if (ks.gt.1) then
          do ik=1,problem%n
            m(:,:,ik)=symconf_t(ik)*m(:,:,ik)
            l(:,:,ik)=symconf_t(ik)*l(:,:,ik)
          end do
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
        !$omp critical
        select case (boundary(sb_int)%coupling)
          ! BE OR BE-FE BOUNDARY
          case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
            select case (boundary(sb_int)%class)
              ! ORDINARY BOUNDARY
              case (fbem_boundary_class_ordinary)
                do il=1,problem%n
                  do ik=1,problem%n
                    do kn_int=1,se_int_n_nodes
                      sn_int=element(se_int)%node(kn_int)
                      internalpoint(sip)%value_r(il,kc)=internalpoint(sip)%value_r(il,kc)&
                                                       +lp(kn_int,il,ik)*node(sn_int)%value_r(problem%n+ik,1)&
                                                       -mp(kn_int,il,ik)*node(sn_int)%value_r(          ik,1)
                    end do
                  end do
                end do
              ! CRACK-LIKE BOUNDARY
              case (fbem_boundary_class_cracklike)
                do il=1,problem%n
                  do ik=1,problem%n
                    do kn_int=1,se_int_n_nodes
                      sn_int=element(se_int)%node(kn_int)
                      internalpoint(sip)%value_r(il,kc)=internalpoint(sip)%value_r(il,kc)&
                                                       +lp(kn_int,il,ik)*node(sn_int)%value_r(problem%n+ik,1)&
                                                       -mp(kn_int,il,ik)*node(sn_int)%value_r(          ik,1)&
                                                       +lm(kn_int,il,ik)*node(sn_int)%value_r(problem%n+ik,2)&
                                                       -mm(kn_int,il,ik)*node(sn_int)%value_r(          ik,2)
                    end do
                  end do
                end do
            end select
          ! BE-BE OR BE-FE-BE BOUNDARY
          case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
            if (sb_int_reversion.eqv.(.false.)) then
              do il=1,problem%n
                do ik=1,problem%n
                  do kn_int=1,se_int_n_nodes
                    sn_int=element(se_int)%node(kn_int)
                    internalpoint(sip)%value_r(il,kc)=internalpoint(sip)%value_r(il,kc)&
                                                     +lp(kn_int,il,ik)*node(sn_int)%value_r(problem%n+ik,1)&
                                                     -mp(kn_int,il,ik)*node(sn_int)%value_r(          ik,1)
                  end do
                end do
              end do
            else
              do il=1,problem%n
                do ik=1,problem%n
                  do kn_int=1,se_int_n_nodes
                    sn_int=element(se_int)%node(kn_int)
                    internalpoint(sip)%value_r(il,kc)=internalpoint(sip)%value_r(il,kc)&
                                                     +lp(kn_int,il,ik)*node(sn_int)%value_r(problem%n+ik,2)&
                                                     -mp(kn_int,il,ik)*node(sn_int)%value_r(          ik,2)
                  end do
                end do
              end do
            end if
        end select
        !$omp end critical
      end do

    end do ! Loop through INTERNAL POINTS

  end do ! Loop through SYMMETRICAL ELEMENTS

end subroutine calculate_internal_points_mechanics_bem_staela_element


subroutine calculate_internal_points_mechanics_bem_staela_bl(kr,sb_int,se_int,se_int_n_nodes,mu,nu)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_numerical
  use fbem_symmetry
  use fbem_bem_general
  use fbem_bem_staela2d
  use fbem_bem_staela3d

  ! Module of problem variables
  use problem_variables

  ! No implicit variables
  implicit none

  ! I/O variables
  integer           :: kr
  integer           :: sb_int
  integer           :: se_int
  integer           :: se_int_n_nodes
  real(kind=real64) :: mu
  real(kind=real64) :: nu
  ! Local variables
  integer                :: ks
  integer                :: il, ik, kc
  integer                :: kn_int, sn_int
  type(fbem_bem_element) :: se_int_data
  logical                :: se_int_reversion
  integer                :: kip, sip
  integer                :: kn
  real(kind=real64)      :: x_i(problem%n), n_i(problem%n)
  ! Kernels for SBIE integration
  real(kind=real64)      :: g (se_int_n_nodes,problem%n,problem%n)
  ! Kernels for HBIE integration
  real(kind=real64)      :: l (se_int_n_nodes,problem%n,problem%n)
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
          call fbem_bem_staela2d_sbie_bl_auto(se_int_data,x_i,mu,nu,qsi_parameters,qsi_ns_max,g)
        case (3)
          call fbem_bem_staela3d_sbie_bl_auto(se_int_data,x_i,mu,nu,qsi_parameters,qsi_ns_max,g)
      end select
      ! Additional kernels for half-space fundamental solution
      if (region(kr)%space.eq.fbem_half_space) then
        select case (problem%n)
          case (2)
            x_i(abs(region(kr)%halfspace_n))=2.d0*region(kr)%halfspace_x-x_i(abs(region(kr)%halfspace_n))
            n_i(abs(region(kr)%halfspace_n))=-n_i(abs(region(kr)%halfspace_n))
            stop 'Error: staela: 2D half-space HBIE for BE body loads not available'
          case (3)
            stop 'Error: staela: 3D half-space HBIE for BE body loads not available'
        end select
        select case (region(kr)%halfspace_bc)
          case (0)
            stop 'Error: staela: region(kr)%halfspace_bc=0 not yet'
          case (1)
        end select
      end if
      ! BUILD KERNELS ACCORDING TO SYMMETRY
      if (ks.gt.1) then
        do ik=1,problem%n
          g(:,:,ik)=symconf_t(ik)*g(:,:,ik)
        end do
      end if
      ! ASSEMBLE
      !$omp critical
      do il=1,problem%n
        do ik=1,problem%n
          do kn_int=1,se_int_n_nodes
            sn_int=element(se_int)%node(kn_int)
            internalpoint(sip)%value_r(il,0)=internalpoint(sip)%value_r(il,0)+g(kn_int,il,ik)*node(sn_int)%value_r(problem%n+ik,1)
          end do
        end do
      end do
      !$omp end critical

      ! ======== !
      ! TRACTION !
      ! ======== !

      do kc=1,problem%n
        ! UNIT NORMAL
        n_i=0.0d0
        n_i(kc)=1.0d0
        ! CALCULATE KERNELS
        select case (problem%n)
          case (2)
            call fbem_bem_staela2d_hbie_bl_auto(se_int_data,x_i,n_i,mu,nu,qsi_parameters,qsi_ns_max,l)
          case (3)
            call fbem_bem_staela3d_hbie_bl_auto(se_int_data,x_i,n_i,mu,nu,qsi_parameters,qsi_ns_max,l)
        end select
        ! Additional kernels for half-space fundamental solution
        if (region(kr)%space.eq.fbem_half_space) then
          select case (problem%n)
            case (2)
              x_i(abs(region(kr)%halfspace_n))=2.d0*region(kr)%halfspace_x-x_i(abs(region(kr)%halfspace_n))
              n_i(abs(region(kr)%halfspace_n))=-n_i(abs(region(kr)%halfspace_n))
              stop 'Error: staela: 2D half-space HBIE for BE body loads not available'
            case (3)
              stop 'Error: staela: 3D half-space HBIE for BE body loads not available'
          end select
          select case (region(kr)%halfspace_bc)
            case (0)
              stop 'Error: staela: region(kr)%halfspace_bc=0 not yet'
            case (1)
          end select
        end if
        ! BUILD KERNELS ACCORDING TO SYMMETRY
        if (ks.gt.1) then
          do ik=1,problem%n
            l(:,:,ik)=symconf_t(ik)*l(:,:,ik)
          end do
        end if
        ! ASSEMBLE
        !$omp critical
        do il=1,problem%n
          do ik=1,problem%n
            do kn_int=1,se_int_n_nodes
              sn_int=element(se_int)%node(kn_int)
              internalpoint(sip)%value_r(il,kc)=internalpoint(sip)%value_r(il,kc)+l(kn_int,il,ik)*node(sn_int)%value_r(problem%n+ik,1)
            end do
          end do
        end do
        !$omp end critical
      end do

    end do ! Loop through INTERNAL POINTS

  end do ! Loop through SYMMETRICAL ELEMENTS

end subroutine calculate_internal_points_mechanics_bem_staela_bl
