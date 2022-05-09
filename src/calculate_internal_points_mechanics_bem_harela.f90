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


subroutine calculate_internal_points_mechanics_bem_harela(kf,kr)

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
  use fbem_bem_harela2d
  use fbem_bem_harela3d
  use fbem_bem_harpor2d
  use fbem_bem_harpor3d

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
  integer                           :: ke_int, se_int
  integer                           :: se_int_n_nodes
  real(kind=real64)                 :: rd(problem%n), x_i(problem%n)
  ! Region properties
  real(kind=real64)                  :: rho, rho1, rho2, rhoa ! Densities (including poro->elastic)
  real(kind=real64)                  :: b                     ! Coupling parameters (poro->elastic)
  type(fbem_bem_harpor2d_parameters) :: p_harpor2d            ! Poroelastic parameters (poro->elastic)
  complex(kind=real64)               :: mu, lambda, R, Q, nu  ! Lame's elastic constants
  complex(kind=real64)               :: c1, c2                ! P and S Wave propagation speeds
  complex(kind=real64)               :: k1, k2                ! P and S wave number
  type(fbem_bem_harela3d_parameters) :: p3d
  type(fbem_bem_harela2d_parameters) :: p2d
  ! Writing
  character(len=fbem_fmtstr)        :: fmtstr            ! String used for write format string

  ! Message
  if (verbose_level.ge.2) then
    write(fmtstr,*) '(1x,a6,i',fbem_nchar_int(region(kr)%id),',1x,a26)'
    call fbem_trimall(fmtstr)
    write(output_unit,fmtstr) 'Region', region(kr)%id, '(BE region, elastic solid)'
  end if

  ! Frequency
  omega=frequency(kf)

  ! Save the region properties to local variables
  select case (region(kr)%subtype)
    !
    ! Elastic or viscoelastic
    !
    case (0,fbem_elastic)
      rho=region(kr)%property_r(1)
      nu=region(kr)%property_c(3)
      mu=region(kr)%property_c(2)
      lambda=region(kr)%property_c(6)
      c1=region(kr)%property_c(7)
      c2=region(kr)%property_c(8)
    !
    ! Bardet's viscoelasticity model of poroelasticity
    !
    case (fbem_bardet)
      rho=region(kr)%property_r(13)+region(kr)%property_r(14)
      c1=region(kr)%property_c(7)*zsqrt(1.0d0+c_im*omega*region(kr)%property_c(9 ))
      c2=region(kr)%property_c(8)*zsqrt(1.0d0+c_im*omega*region(kr)%property_c(10))
      mu=rho*c2**2
      lambda=rho*c1**2-2.0d0*mu
      nu=0.5d0*lambda/(lambda+mu)
    !
    ! Bougacha-Roesset-Tassoulas viscoelasticity model of poroelasticity
    !
    case (fbem_brt_cp1,fbem_brt_cp2,fbem_brt_cpm)
      ! Save the region properties to local variables
      lambda=region(kr)%property_c(3)
      mu=region(kr)%property_c(4)
      rho1=region(kr)%property_r(13)
      rho2=region(kr)%property_r(14)
      rhoa=region(kr)%property_r(9)
      R=region(kr)%property_r(10)
      Q=region(kr)%property_r(11)
      b=region(kr)%property_r(12)
      ! Obtain the wave velocities
      call fbem_bem_harpor2d_calculate_basic_parameters(lambda,mu,rho1,rho2,rhoa,R,Q,b,omega,p_harpor2d)
      ! Corresponding isomorphism
      rho=region(kr)%property_r(13)+region(kr)%property_r(14)
      select case (region(kr)%subtype)
        case (fbem_brt_cp1)
          c1=p_harpor2d%c1
        case (fbem_brt_cp2)
          c1=p_harpor2d%c2
        case (fbem_brt_cpm)
          c1=0.5d0*(p_harpor2d%c1+p_harpor2d%c2)
      end select
      c2=p_harpor2d%c3
      mu=rho*c2**2
      lambda=rho*c1**2-2.0d0*mu
      nu=0.5d0*lambda/(lambda+mu)
  end select
  select case (problem%n)
    case (2)
      call fbem_bem_harela2d_calculate_parameters(lambda,mu,rho,omega,p2d)
    case (3)
      call fbem_bem_harela3d_calculate_parameters(lambda,mu,rho,omega,p3d)
  end select

  ! ============================================= !
  ! OBTAIN DIFFRACTED FIELD INTERNAL POINT VALUES !
  ! ============================================= !

  !
  ! Loop through the BOUNDARIES of the REGION for INTEGRATION
  !
  do kb_int=1,region(kr)%n_boundaries
    ! INTEGRATION BOUNDARY
    sb_int=region(kr)%boundary(kb_int)
    sb_int_reversion=region(kr)%boundary_reversion(kb_int)

    ! Message
    if (verbose_level.ge.3) then
      write(fmtstr,*) '(2x,a20,1x,i',fbem_nchar_int(boundary(sb_int)%id),',1x,a8)'
      call fbem_trimall(fmtstr)
      if (sb_int_reversion) then
        write(output_unit,fmtstr) 'Integrating boundary', boundary(sb_int)%id, '(-n) ...'
      else
        write(output_unit,fmtstr) 'Integrating boundary', boundary(sb_int)%id, '(+n) ...'
      end if
    end if

    !
    ! Loop through the ELEMENTS of the BOUNDARY for INTEGRATION
    !
    !$omp parallel do schedule (dynamic) default (shared) private (se_int,se_int_n_nodes,fmtstr)
    do ke_int=1,part(boundary(sb_int)%part)%n_elements
      se_int=part(boundary(sb_int)%part)%element(ke_int)
      se_int_n_nodes=element(se_int)%n_nodes

      ! Message
      if (verbose_level.ge.4) then
        write(fmtstr,*) '(3x,a19,1x,i',fbem_nchar_int(element(se_int)%id),',1x,a3)'
        call fbem_trimall(fmtstr)
        if (sb_int_reversion) then
          write(output_unit,fmtstr) 'Integrating element', element(se_int)%id, '...'
        else
          write(output_unit,fmtstr) 'Integrating element', element(se_int)%id, '...'
        end if
      end if

      ! Build and assemble the element
      call calculate_internal_points_mechanics_bem_harela_element(omega,kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,p2d,p3d)

      ! Ending message
      if (verbose_level.ge.4) write(output_unit,'(3x,a)') 'done.'

    end do ! Loop through the ELEMENTS of the BOUNDARY for INTEGRATION
    !$omp end parallel do

  end do ! Loop through the BOUNDARIES of the REGION for INTEGRATION

  ! ===============================================
  ! OBTAIN THE TOTAL FIELD OF INTERNAL POINT VALUES
  ! ===============================================

  do kip=1,region(kr)%n_internalpoints
    sip=region(kr)%internalpoint(kip)
    internalpoint(sip)%value_c(:,:)=internalpoint(sip)%value_c(:,:)+internalpoint(sip)%incident_c(:,:)
  end do

      !
      ! Copia los valores nodales a los puntos internos en ellos (solo variable primaria)
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
          internalpoint(sip)%value_c=0.d0
          internalpoint(sip)%value_c(1:problem%n,0)=node(sn)%value_c(1:problem%n,1)
        end if
      end do
    end do
  end do




end subroutine calculate_internal_points_mechanics_bem_harela

subroutine calculate_internal_points_mechanics_bem_harela_element(omega,kr,sb_int,sb_int_reversion,se_int,se_int_n_nodes,p2d,p3d)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_numerical
  use fbem_symmetry
  use fbem_bem_general
  use fbem_bem_harela2d
  use fbem_bem_harela3d

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
  type(fbem_bem_harela3d_parameters) :: p3d
  type(fbem_bem_harela2d_parameters) :: p2d
  ! Local variables
  integer                :: ks
  integer                :: il, ik, kc
  integer                :: kn_int, sn_int
  type(fbem_bem_element) :: se_int_data
  logical                :: se_int_reversion
  integer                :: kip, sip
  integer                :: kn
  real(kind=real64)      :: x_i(problem%n), n_i(problem%n)
  ! Incident wave variables
  complex(kind=real64)   :: u_inc(problem%n,se_int_n_nodes), t_inc(problem%n,se_int_n_nodes)
  ! Kernels for SBIE integration
  complex(kind=real64)   :: h (se_int_n_nodes,problem%n,problem%n), g (se_int_n_nodes,problem%n,problem%n)
  complex(kind=real64)   :: hp(se_int_n_nodes,problem%n,problem%n), gp(se_int_n_nodes,problem%n,problem%n)
  complex(kind=real64)   :: hm(se_int_n_nodes,problem%n,problem%n), gm(se_int_n_nodes,problem%n,problem%n)
  ! Kernels for HBIE integration
  complex(kind=real64)   :: m (se_int_n_nodes,problem%n,problem%n), l (se_int_n_nodes,problem%n,problem%n)
  complex(kind=real64)   :: mp(se_int_n_nodes,problem%n,problem%n), lp(se_int_n_nodes,problem%n,problem%n)
  complex(kind=real64)   :: mm(se_int_n_nodes,problem%n,problem%n), lm(se_int_n_nodes,problem%n,problem%n)
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
  select case (boundary(sb_int)%coupling)
    case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
      u_inc=element(se_int)%incident_c(1:problem%n,:,1)
      t_inc=element(se_int)%incident_c((problem%n+1):(2*problem%n),:,1)
    case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
      if (sb_int_reversion) then
        u_inc=element(se_int)%incident_c(1:problem%n,:,2)
        t_inc=element(se_int)%incident_c((problem%n+1):(2*problem%n),:,2)
      else
        u_inc=element(se_int)%incident_c(1:problem%n,:,1)
        t_inc=element(se_int)%incident_c((problem%n+1):(2*problem%n),:,1)
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

      ! ------------ !
      ! DISPLACEMENT !
      ! ------------ !

      ! CALCULATE KERNELS
      select case (problem%n)
        case (2)
          call fbem_bem_harela2d_sbie_auto(1,se_int_data,se_int_reversion,x_i,p2d,qsi_parameters,qsi_ns_max,h,g)
        case (3)
          call fbem_bem_harela3d_sbie_auto(se_int_data,se_int_reversion,x_i,p3d,qsi_parameters,qsi_ns_max,h,g)
      end select
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
                    internalpoint(sip)%value_c(il,0)=internalpoint(sip)%value_c(il,0)&
                                                    +gp(kn_int,il,ik)*(node(sn_int)%value_c(problem%n+ik,1)-t_inc(ik,kn_int))&
                                                    -hp(kn_int,il,ik)*(node(sn_int)%value_c(          ik,1)-u_inc(ik,kn_int))
                  end do
                end do
              end do
            ! CRACK-LIKE BOUNDARY
            case (fbem_boundary_class_cracklike)
              do il=1,problem%n
                do ik=1,problem%n
                  do kn_int=1,se_int_n_nodes
                    sn_int=element(se_int)%node(kn_int)
                    internalpoint(sip)%value_c(il,0)=internalpoint(sip)%value_c(il,0)&
                                                    +gp(kn_int,il,ik)*(node(sn_int)%value_c(problem%n+ik,1)-t_inc(ik,kn_int))&
                                                    -hp(kn_int,il,ik)*(node(sn_int)%value_c(          ik,1)-u_inc(ik,kn_int))&
                                                    +gm(kn_int,il,ik)*(node(sn_int)%value_c(problem%n+ik,2)+t_inc(ik,kn_int))&
                                                    -hm(kn_int,il,ik)*(node(sn_int)%value_c(          ik,2)-u_inc(ik,kn_int))
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
                  internalpoint(sip)%value_c(il,0)=internalpoint(sip)%value_c(il,0)&
                                                  +gp(kn_int,il,ik)*(node(sn_int)%value_c(problem%n+ik,1)-t_inc(ik,kn_int))&
                                                  -hp(kn_int,il,ik)*(node(sn_int)%value_c(          ik,1)-u_inc(ik,kn_int))
                end do
              end do
            end do
          else
            do il=1,problem%n
              do ik=1,problem%n
                do kn_int=1,se_int_n_nodes
                  sn_int=element(se_int)%node(kn_int)
                  internalpoint(sip)%value_c(il,0)=internalpoint(sip)%value_c(il,0)&
                                                  +gp(kn_int,il,ik)*(node(sn_int)%value_c(problem%n+ik,2)-t_inc(ik,kn_int))&
                                                  -hp(kn_int,il,ik)*(node(sn_int)%value_c(          ik,2)-u_inc(ik,kn_int))
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
            call fbem_bem_harela2d_hbie_auto(1,se_int_data,se_int_reversion,x_i,n_i,p2d,qsi_parameters,qsi_ns_max,m,l)
          case (3)
            call fbem_bem_harela3d_hbie_auto(se_int_data,se_int_reversion,x_i,n_i,p3d,qsi_parameters,qsi_ns_max,m,l)
        end select
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
                      internalpoint(sip)%value_c(il,kc)=internalpoint(sip)%value_c(il,kc)&
                                                       +lp(kn_int,il,ik)*(node(sn_int)%value_c(problem%n+ik,1)-t_inc(ik,kn_int))&
                                                       -mp(kn_int,il,ik)*(node(sn_int)%value_c(          ik,1)-u_inc(ik,kn_int))
                    end do
                  end do
                end do
              ! CRACK-LIKE BOUNDARY
              case (fbem_boundary_class_cracklike)
                do il=1,problem%n
                  do ik=1,problem%n
                    do kn_int=1,se_int_n_nodes
                      sn_int=element(se_int)%node(kn_int)
                      internalpoint(sip)%value_c(il,kc)=internalpoint(sip)%value_c(il,kc)&
                                                       +lp(kn_int,il,ik)*(node(sn_int)%value_c(problem%n+ik,1)-t_inc(ik,kn_int))&
                                                       -mp(kn_int,il,ik)*(node(sn_int)%value_c(          ik,1)-u_inc(ik,kn_int))&
                                                       +lm(kn_int,il,ik)*(node(sn_int)%value_c(problem%n+ik,2)+t_inc(ik,kn_int))&
                                                       -mm(kn_int,il,ik)*(node(sn_int)%value_c(          ik,2)-u_inc(ik,kn_int))
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
                    internalpoint(sip)%value_c(il,kc)=internalpoint(sip)%value_c(il,kc)&
                                                     +lp(kn_int,il,ik)*(node(sn_int)%value_c(problem%n+ik,1)-t_inc(ik,kn_int))&
                                                     -mp(kn_int,il,ik)*(node(sn_int)%value_c(          ik,1)-u_inc(ik,kn_int))
                  end do
                end do
              end do
            else
              do il=1,problem%n
                do ik=1,problem%n
                  do kn_int=1,se_int_n_nodes
                    sn_int=element(se_int)%node(kn_int)
                    internalpoint(sip)%value_c(il,kc)=internalpoint(sip)%value_c(il,kc)&
                                                     +lp(kn_int,il,ik)*(node(sn_int)%value_c(problem%n+ik,2)-t_inc(ik,kn_int))&
                                                     -mp(kn_int,il,ik)*(node(sn_int)%value_c(          ik,2)-u_inc(ik,kn_int))
                  end do
                end do
              end do
            end if
        end select
        !$omp end critical
      end do

    end do ! Loop through INTERNAL POINTS

  end do ! Loop through SYMMETRICAL ELEMENTS

end subroutine calculate_internal_points_mechanics_bem_harela_element
