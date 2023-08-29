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

subroutine calculate_incident_mechanics_harmonic(kf)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_numerical
  use fbem_shape_functions
  use fbem_harpot_incident_field
  use fbem_harela_incident_field
  use fbem_harpor_incident_field
  use fbem_bem_harpot2d
  use fbem_bem_harpot3d
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
  ! Local variables
  real(kind=real64)                 :: omega
  integer                           :: kif, sif
  integer                           :: kr, kb, sb, sp, ke, se, kn, sn, knj, kc, kc2, kl
  integer                           :: kip, sip
  logical                           :: sb_reversion
  integer                           :: se_n_nodes
  integer                           :: k_start, k_end, face
  real(kind=real64), allocatable    :: x_fn(:,:), n_fn(:,:)
  complex(kind=real64), allocatable :: ui_ref(:), ti_ref(:) ! Campo incidente (ui,ti) en el punto de referencia
  complex(kind=real64)              :: a_norm               ! Amplitud usada para la normalizacion: (ui,ti)_normalizada=(ui,ti)*a_norm (a_norm=amplitud_definida/a_ref)
  complex(kind=real64), allocatable :: ui(:), ti(:)
  complex(kind=real64)              :: pi, Uni
  real(kind=real64)                 :: xn(3), nn(3)
  real(kind=real64)                 :: x_i(problem%n), n_i(problem%n)
  ! Material properties
  real(kind=real64)                 :: phi                           ! Porosity (poroelastic medium)
  real(kind=real64)                 :: rho                           ! Density (inviscid fluid and elastic solid)
  real(kind=real64)                 :: rhof, rhos, rho1, rho2, rhoa  ! Densities (porolastic medium)
  complex(kind=real64)              :: rhohat11, rhohat12, rhohat22  ! Effective densities (porelastic medium)
  complex(kind=real64)              :: mu, lambda, nu                ! Lame's elastic constants
  complex(kind=real64)              :: R, Q                          ! Biot's coupling parameters
  real(kind=real64)                 :: b                             ! Dissipation constant b
  type(fbem_bem_harpor2d_parameters):: p_harpor2d                    ! Poroelastic parameters (poro->elastic)
  complex(kind=real64)              :: c1, c2, c3                    ! Wave propagation speeds
  complex(kind=real64)              :: k1, k2, k3                    ! Wavenumbers
  complex(kind=real64)              :: c, k                          ! Wave propagation speed and wavenumber (inviscid fluid)



  ! Phurkhao incident field constants
  complex(kind=real64)              :: phi_s_1, phi_f_1
  complex(kind=real64)              :: phi_s_2, phi_f_2
  complex(kind=real64)              :: D1I

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'START calculating incident fields')

  ! Frequency
  omega=frequency(kf)

  ! Initialization of incident fields over nodes, elements, and internal points
  do kr=1,n_regions
    if (region(kr)%class.eq.fbem_be) then
      ! BOUNDARIES
      do kb=1,region(kr)%n_boundaries
        sb=region(kr)%boundary(kb)
        sp=boundary(sb)%part
        do ke=1,part(sp)%n_elements
          se=part(sp)%element(ke)
          element(se)%incident_c=0.
          do kn=1,element(se)%n_nodes
            sn=element(se)%node(kn)
            node(sn)%incident_c=0.
          end do
        end do
      end do
      ! BODY LOADS
      do kb=1,region(kr)%n_be_bodyloads
        sb=region(kr)%be_bodyload(kb)
        sp=be_bodyload(sb)%part
        do ke=1,part(sp)%n_elements
          se=part(sp)%element(ke)
          element(se)%incident_c=0.
          do kn=1,element(se)%n_nodes
            sn=element(se)%node(kn)
            node(sn)%incident_c=0.
          end do
        end do
      end do
      ! INTERNAL POINTS
      do kip=1,region(kr)%n_internalpoints
        sip=region(kr)%internalpoint(kip)
        internalpoint(sip)%incident_c=0.
      end do
    end if
  end do

  ! Initialization of incident fields
  do kif=1,n_incidentfields
    if ((incidentfield(kif)%class.eq.fbem_plane).and.(incidentfield(kif)%space.eq.fbem_multilayered_half_space)) then


      ! se asume campos incidentes con todos los materiales iguales
      if (material(incidentfield(kif)%layer_material(1))%type.eq.'biot_poroelastic_medium') then



        !
        ! Hay que meter la dependencia de normal positiva o negativa de la superficie libre como en el caso elastodinamico
        !




        if (.not.allocated(incidentfield(kif)%layer_A)) allocate(incidentfield(kif)%layer_A(3,2,8,incidentfield(kif)%n_layers))
        if (.not.allocated(incidentfield(kif)%layer_k)) allocate(incidentfield(kif)%layer_k(3,incidentfield(kif)%n_layers))
        call fbem_harpor_vertical_plane_wave_layered_halfspace_setup(incidentfield(kif)%n_layers,incidentfield(kif)%layer_ztop,&
        incidentfield(kif)%layer_properties,incidentfield(kif)%layer_cvalue,omega,&
        incidentfield(kif)%layer_A,incidentfield(kif)%layer_k)
      end if



    end if
  end do





!!! falta meter el caso estratificado!!!!!!!!!!!!!!!!!!!¿?¿?¿?¿?¿?









  ! Loop through REGIONS
  do kr=1,n_regions

    ! If BE REGIONS with >0 INCIDENT FIELDS
    if ((region(kr)%class.eq.fbem_be).and.(region(kr)%n_incidentfields.gt.0)) then

      ! =====================
      ! REGION INITIALIZATION
      ! =====================

      select case (region(kr)%type)

        ! --------------
        ! INVISCID FLUID
        ! --------------

        case (fbem_potential)
          ! DOF: p, Un
          k_start=1
          k_end  =2
          ! Save the region properties to local variables
          rho=region(kr)%property_r(1)
          c=region(kr)%property_c(4)
          ! Wavenumber
          k=omega/c

        ! -------------
        ! ELASTIC SOLID
        ! -------------

        case (fbem_viscoelastic)
          ! DOF: u_k, t_k
          k_start=1
          k_end  =2*problem%n


          if (region(kr)%space.eq.fbem_multilayered_half_space) then

          else
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
            ! Wavenumbers
            k1=omega/c1
            k2=omega/c2
          end if


          ! Initialize incident fields
          do kif=1,region(kr)%n_incidentfields
            sif=region(kr)%incidentfield(kif)


            ! Discriminar mejor los casos y tal.....
            if (incidentfield(sif)%space.eq.fbem_multilayered_half_space) then

              if (.not.allocated(incidentfield(sif)%layer_amplitudes)) then
               allocate(incidentfield(sif)%layer_amplitudes(incidentfield(sif)%n_layers,2))
              end if
              call fbem_harela_incident_plane_wave_multilayered_amplitudes(omega,incidentfield(sif)%np,&
              incidentfield(sif)%n_layers,incidentfield(sif)%layer_ztop,incidentfield(sif)%layer_lambda,incidentfield(sif)%layer_mu,&
              incidentfield(sif)%layer_rho,incidentfield(sif)%wave_type,incidentfield(sif)%layer_amplitudes)

            end if



          end do

        ! ------------------
        ! POROELASTIC MEDIUM
        ! ------------------

        case (fbem_poroelastic)
          ! DOF: tau, u_k, Un, t_k
          k_start=0
          k_end  =1+2*problem%n
          ! Save the region properties to local variables
          rhof=region(kr)%property_r(1)
          rhos=region(kr)%property_r(2)
          lambda=region(kr)%property_c(3)
          mu=region(kr)%property_c(4)
          phi=region(kr)%property_r(8)
          rhoa=region(kr)%property_r(9)
          rho1=region(kr)%property_r(13)
          rho2=region(kr)%property_r(14)
          R=region(kr)%property_c(10)
          Q=region(kr)%property_c(11)
          b=region(kr)%property_r(12)
          call fbem_bem_harpor2d_calculate_basic_parameters(lambda,mu,rho1,rho2,rhoa,R,Q,b,omega,p_harpor2d)
          k1=p_harpor2d%k1
          k2=p_harpor2d%k2
          k3=p_harpor2d%k3
          ! Densities rhohat
          rhohat11=p_harpor2d%rhohat11
          rhohat12=p_harpor2d%rhohat12
          rhohat22=p_harpor2d%rhohat22
          ! Phurkhao incident field constants
          D1I    =0.5d0
          phi_s_1=1.0d0
          phi_f_1=(k1**2*(lambda+2.0d0*mu+Q**2/R)-omega**2*rhohat11)/(omega**2*rhohat12-k1**2*Q)
          phi_s_2=1.0d0
          phi_f_2=(k2**2*(lambda+2.0d0*mu+Q**2/R)-omega**2*rhohat11)/(omega**2*rhohat12-k2**2*Q)
      end select

      ! ========== !
      ! BOUNDARIES !
      ! ========== !

      do kb=1,region(kr)%n_boundaries
        sb=region(kr)%boundary(kb)
        sb_reversion=region(kr)%boundary_reversion(kb)
        sp=boundary(sb)%part
        select case (boundary(sb)%coupling)
          case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
            face=1
          case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
            if (sb_reversion) then
              face=2
            else
              face=1
            end if
        end select

        ! =========================== !
        ! ELEMENT-WISE INCIDENT FIELD !
        ! =========================== !

        ! Loop through the ELEMENTS of the BOUNDARY
        do ke=1,part(sp)%n_elements
          se=part(sp)%element(ke)
          se_n_nodes=element(se)%n_nodes
          allocate (x_fn(problem%n,se_n_nodes),n_fn(problem%n,se_n_nodes))
          x_fn=element(se)%x_fn
          n_fn=element(se)%n_fn
          if (sb_reversion) n_fn=-n_fn
          select case (region(kr)%type)

            ! --------------
            ! INVISCID FLUID
            ! --------------

            case (fbem_potential)
              ! Loop through INCIDENT FIELDS of the REGION
              do kif=1,region(kr)%n_incidentfields
                sif=region(kr)%incidentfield(kif)
                ! Loop through the NODES of the ELEMENT
                do kn=1,se_n_nodes
                  sn=element(se)%node(kn)
                  select case (incidentfield(sif)%class)
                    case (fbem_point)
                      call fbem_harpot_pointwave(problem%n,omega,rho,c,incidentfield(sif)%amplitude,incidentfield(sif)%x0,&
                                                 incidentfield(sif)%space,incidentfield(sif)%np,incidentfield(sif)%xp,incidentfield(sif)%bc,&
                                                 incidentfield(sif)%symconf,incidentfield(sif)%xs,x_fn(:,kn),n_fn(:,kn),pi,Uni)
                    case (fbem_plane)
                      call fbem_harpot_planewave(problem%n,omega,rho,c,&
                           incidentfield(sif)%amplitude,incidentfield(sif)%x0,incidentfield(sif)%varphi,incidentfield(sif)%theta,&
                           incidentfield(sif)%space,incidentfield(sif)%np,incidentfield(sif)%xp,incidentfield(sif)%bc,&
                           incidentfield(sif)%symconf,incidentfield(sif)%xs,x_fn(:,kn),n_fn(:,kn),pi,Uni)
                  end select
                  ! Depending on the leading variable of the incident wave field
                  select case (incidentfield(sif)%variable)
                    ! Incident wave field in terms of pressure
                      case (0)
                        element(se)%incident_c(1,kn,face)=element(se)%incident_c(1,kn,face)+pi
                        element(se)%incident_c(2,kn,face)=element(se)%incident_c(2,kn,face)+Uni
                    ! Incident wave field in terms of normal displacements
                    case (1)
                      stop 'not implemented yet'
                  end select
                end do
              end do
              ! Treatment for traction quarter point elements
              select case (element(se)%type_f2)
                case (fbem_line2,fbem_line3,fbem_line4)
                case (fbem_line3_qp1t,fbem_line3_mqp1t)
                  element(se)%incident_c(2,1,face)=0.d0
                  element(se)%incident_c(2,3,face)=0.5d0*element(se)%incident_c(2,3,face)
                case (fbem_line3_qp2t,fbem_line3_mqp2t)
                  element(se)%incident_c(2,2,face)=0.d0
                  element(se)%incident_c(2,3,face)=0.5d0*element(se)%incident_c(2,3,face)
                case (fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9)
                case default
                  stop 'stop at calculate_incident_mechanics_harmoni: check element quarter point'
              end select
              ! Crack-like boundaries
              if (boundary(sb)%class.eq.fbem_boundary_class_cracklike) then
                do kn=1,se_n_nodes
                  do kc=1,problem%n
                    element(se)%incident_c(1,kn,2)= element(se)%incident_c(1,kn,1)
                    element(se)%incident_c(2,kn,2)=-element(se)%incident_c(2,kn,1)
                  end do
                end do
              end if

            ! -------------
            ! ELASTIC SOLID
            ! -------------

            case (fbem_viscoelastic)
              allocate (ui(3),ti(3))
              ! Loop through INCIDENT FIELDS of the REGION
              do kif=1,region(kr)%n_incidentfields
                sif=region(kr)%incidentfield(kif)
                ! Loop through the NODES of the ELEMENT
                do kn=1,se_n_nodes
                  sn=element(se)%node(kn)
                  ! Copy to arrays of dimension 3
                  do kc=1,problem%n
                    xn(kc)=x_fn(kc,kn)
                    nn(kc)=n_fn(kc,kn)
                  end do
                  ! Calculate incident field in terms of displacements
                  select case  (incidentfield(sif)%space)
                    case (fbem_full_space,fbem_half_space)

                    ! mirar esto para nu complejo

                      call fbem_harela_incident_plane_wave(problem%n,k1,k2,mu,lambda,dreal(nu),&
                           incidentfield(sif)%space,incidentfield(sif)%symconf(2),incidentfield(sif)%wave_type,&
                           incidentfield(sif)%varphi,incidentfield(sif)%theta,xn,incidentfield(sif)%xp,nn,ui,ti)
                      ui=0.5d0*ui
                      ti=0.5d0*ti
                    case (fbem_multilayered_half_space)
                      call fbem_harela_incident_plane_wave_vertical_multilayered(problem%n,omega,incidentfield(sif)%np,&
                      incidentfield(sif)%n_layers,incidentfield(sif)%layer_ztop,incidentfield(sif)%layer_lambda,&
                      incidentfield(sif)%layer_mu,incidentfield(sif)%layer_rho,incidentfield(sif)%symconf(2),&
                      incidentfield(sif)%wave_type,incidentfield(sif)%layer_amplitudes,incidentfield(sif)%varphi,xn,nn,ui,ti)
                  end select
                  ! Depending on the leading variable of the incident wave field
                  select case (incidentfield(sif)%variable)
                    ! Incident wave field in terms of displacements
                    case (0)
                      do kc=1,problem%n
                        element(se)%incident_c(          kc,kn,face)=element(se)%incident_c(          kc,kn,face)+ui(kc)
                        element(se)%incident_c(problem%n+kc,kn,face)=element(se)%incident_c(problem%n+kc,kn,face)+ti(kc)
                      end do
                    ! Incident wave field in terms of stresses
                    case (1)
                      if (incidentfield(sif)%space.eq.fbem_multilayered_half_space) stop 'This is not ready yet (multilayered with stress normal.)'
                      select case (incidentfield(sif)%wave_type)
                        ! P-waves
                        case (fbem_harela_p_wave)
                          do kc=1,problem%n
                            element(se)%incident_c(          kc,kn,face)=element(se)%incident_c(          kc,kn,face)+ui(kc)/(-c_im*k1*(lambda+2.d0*mu))
                            element(se)%incident_c(problem%n+kc,kn,face)=element(se)%incident_c(problem%n+kc,kn,face)+ti(kc)/(-c_im*k1*(lambda+2.d0*mu))
                          end do
                        ! S-waves
                        case (fbem_harela_sv_wave,fbem_harela_sh_wave)
                          do kc=1,problem%n
                            element(se)%incident_c(          kc,kn,face)=element(se)%incident_c(          kc,kn,face)+ui(kc)/(-c_im*k2*mu)
                            element(se)%incident_c(problem%n+kc,kn,face)=element(se)%incident_c(problem%n+kc,kn,face)+ti(kc)/(-c_im*k2*mu)
                          end do
                        ! Rayleigh-waves
                        case (fbem_harela_rayleigh_wave)
                          stop 'an incident rayleigh wave in terms of stresses is not implemented yet'
                      end select
                  end select
                end do
              end do
              ! Treatment for traction quarter point elements
              select case (element(se)%type_f2)
                case (fbem_line2,fbem_line3,fbem_line4)
                case (fbem_line3_qp1t,fbem_line3_mqp1t)
                  do kc=1,problem%n
                    element(se)%incident_c(problem%n+kc,1,face)=0.d0
                    element(se)%incident_c(problem%n+kc,3,face)=0.5d0*element(se)%incident_c(problem%n+kc,3,face)
                  end do
                case (fbem_line3_qp2t,fbem_line3_mqp2t)
                  do kc=1,problem%n
                    element(se)%incident_c(problem%n+kc,2,face)=0.d0
                    element(se)%incident_c(problem%n+kc,3,face)=0.5d0*element(se)%incident_c(problem%n+kc,3,face)
                  end do
                case (fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9)
                case default
                  stop 'stop at calculate_incident_mechanics_harmoni: check element quarter point'
              end select
              ! Crack-like boundaries
              if (boundary(sb)%class.eq.fbem_boundary_class_cracklike) then
                do kn=1,se_n_nodes
                  do kc=1,problem%n
                    element(se)%incident_c(          kc,kn,2)= element(se)%incident_c(          kc,kn,1)
                    element(se)%incident_c(problem%n+kc,kn,2)=-element(se)%incident_c(problem%n+kc,kn,1)
                  end do
                end do
              end if
              deallocate (ui,ti)

            ! -----------------
            ! POROELASTIC MEDIUM
            ! -----------------

            case (fbem_poroelastic)

              stop 'please, check poroelastic free-field'


              allocate (ui(0:3),ti(0:3))
              ! Loop through INCIDENT FIELDS of the REGION
              do kif=1,region(kr)%n_incidentfields
                sif=region(kr)%incidentfield(kif)
                ! Loop through the NODES of the ELEMENT
                do kn=1,se_n_nodes
                  sn=element(se)%node(kn)
                  do kc=1,problem%n
                    xn(kc)=x_fn(kc,kn)
                    nn(kc)=n_fn(kc,kn)
                  end do
                  ! Calculate incident field in terms of displacements
                  select case  (incidentfield(sif)%space)
                    case (fbem_full_space,fbem_half_space)
                      select case (incidentfield(sif)%wave_type)
                        case (fbem_harpor_p1_wave,fbem_harpor_p2_wave,fbem_harpor_shx_wave,fbem_harpor_shy_wave)
                          call fbem_harpor_incident_plane_wave(problem%n,omega,k1,k2,k3,lambda,mu,Q,R,rhohat11,rhohat22,rhohat12,&
                               incidentfield(sif)%space,incidentfield(sif)%wave_type,xn,incidentfield(sif)%xp,nn,ui,ti)
                        case (fbem_harpor_R_wave_permeable)
                          if (incidentfield(sif)%space.eq.fbem_full_space) then
                            call fbem_error_message(error_unit,0,'incident field',incidentfield(sif)%id,'Rayleigh incident field is present only in half-space')
                          end if
                          call fbem_harpor_incident_rayleigh_wave(problem%n,lambda,mu,Q,R,phi,rhos,rhof,rhoa,b,omega,xn,nn,ui,ti)
                      end select
                    case (fbem_multilayered_half_space)
                      ! Find kl
                      do kc=incidentfield(sif)%n_layers,1,-1
                        if (element(se)%centroid(problem%n).lt.(incidentfield(sif)%layer_ztop(kc)+1.d-6)) then
                          kl=kc
                          exit
                        end if
                      end do
                      call fbem_harpor_vertical_plane_wave_layered_halfspace(problem%n,incidentfield(sif)%n_layers,&
                      incidentfield(sif)%layer_ztop,incidentfield(sif)%layer_A,incidentfield(sif)%layer_k,incidentfield(sif)%wave_type,1.d-6,kl,xn,nn,ui,ti)
                  end select
                  ! Depending on the leading variable of the incident wave field
                  select case (incidentfield(sif)%variable)
                    ! Incident wave field in terms of displacements
                    case (0)
                      do kc=0,problem%n
                        element(se)%incident_c(            kc,kn,face)=element(se)%incident_c(            kc,kn,face)+ui(kc)
                        element(se)%incident_c(problem%n+1+kc,kn,face)=element(se)%incident_c(problem%n+1+kc,kn,face)+ti(kc)
                      end do
                    ! Incident wave field in terms of total stresses
                    case (1)
                      select case (incidentfield(sif)%wave_type)
                        ! P-waves
                        case (fbem_harpor_p1_wave)
                          do kc=0,problem%n
                            element(se)%incident_c(            kc,kn,face)=element(se)%incident_c(            kc,kn,face)+ui(kc)/(-c_im*k1*D1I*((lambda+2.d0*mu+Q*(1.d0+Q/R))*phi_s_1+(Q+R)*phi_f_1))
                            element(se)%incident_c(problem%n+1+kc,kn,face)=element(se)%incident_c(problem%n+1+kc,kn,face)+ti(kc)/(-c_im*k1*D1I*((lambda+2.d0*mu+Q*(1.d0+Q/R))*phi_s_1+(Q+R)*phi_f_1))
                          end do
                        ! S-waves
                        case (fbem_harpor_shx_wave,fbem_harpor_shy_wave)
                          do kc=0,problem%n
                            element(se)%incident_c(            kc,kn,face)=element(se)%incident_c(            kc,kn,face)+ui(kc)/(-c_im*mu*k3)
                            element(se)%incident_c(problem%n+1+kc,kn,face)=element(se)%incident_c(problem%n+1+kc,kn,face)+ti(kc)/(-c_im*mu*k3)
                          end do
                      end select
                    ! Incident wave field in terms of solid skeleton stresses
                    case (2)
                      select case (incidentfield(sif)%wave_type)
                        ! P-waves
                        case (fbem_harpor_p1_wave)
                          do kc=0,problem%n
                            element(se)%incident_c(            kc,kn,face)=element(se)%incident_c(            kc,kn,face)+ui(kc)/(-c_im*k1*D1I*(phi_s_1*(lambda+2.d0*mu+Q**2/R)+phi_f_1*Q))
                            element(se)%incident_c(problem%n+1+kc,kn,face)=element(se)%incident_c(problem%n+1+kc,kn,face)+ti(kc)/(-c_im*k1*D1I*(phi_s_1*(lambda+2.d0*mu+Q**2/R)+phi_f_1*Q))
                          end do
                        ! S-waves
                        case (fbem_harpor_shx_wave,fbem_harpor_shy_wave)
                          stop 'not implemented yet'
                      end select
                  end select
                end do
              end do
              ! Treatment for traction quarter point elements
              select case (element(se)%type_f2)
                case (fbem_line2,fbem_line3,fbem_line4)
                case (fbem_line3_qp1t,fbem_line3_mqp1t)
                  do kc=1,problem%n
                    element(se)%incident_c(problem%n+1+kc,1,face)=0.d0
                    element(se)%incident_c(problem%n+1+kc,3,face)=0.5d0*element(se)%incident_c(problem%n+1+kc,3,face)
                  end do
                case (fbem_line3_qp2t,fbem_line3_mqp2t)
                  do kc=1,problem%n
                    element(se)%incident_c(problem%n+1+kc,2,face)=0.d0
                    element(se)%incident_c(problem%n+1+kc,3,face)=0.5d0*element(se)%incident_c(problem%n+1+kc,3,face)
                  end do
                case (fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9)
                case default
                  stop 'stop at calculate_incident_mechanics_harmonic: check element quarter point'
              end select
              ! Crack-like boundaries
              if (boundary(sb)%class.eq.fbem_boundary_class_cracklike) then
                do kn=1,se_n_nodes
                  do kc=0,problem%n
                    element(se)%incident_c(            kc,kn,2)= element(se)%incident_c(            kc,kn,1)
                    element(se)%incident_c(problem%n+1+kc,kn,2)=-element(se)%incident_c(problem%n+1+kc,kn,1)
                  end do
                end do
              end if
              deallocate (ui,ti)
          end select

          ! Finalize
          deallocate (x_fn,n_fn)

        end do ! Loop through the ELEMENTS of the BOUNDARY

        ! ======================== !
        ! NODE-WISE INCIDENT FIELD !
        ! ======================== !

        ! Loop through the NODES of the BOUNDARY
        do kn=1,part(sp)%n_nodes
          sn=part(sp)%node(kn)
          do ke=1,node(sn)%n_elements
            se=node(sn)%element(ke)
            knj=node(sn)%element_node_iid(ke)
            node(sn)%incident_c(:,:)=node(sn)%incident_c(:,:)+element(se)%incident_c(:,knj,:)
          end do
          node(sn)%incident_c=node(sn)%incident_c/real(node(sn)%n_elements)
        end do

      end do ! Loop through the BOUNDARIES of the REGION

      ! ========== !
      ! BODY LOADS !
      ! ========== !

      do kb=1,region(kr)%n_be_bodyloads
        sb=region(kr)%be_bodyload(kb)
        sp=be_bodyload(sb)%part



        !
        ! The same code as BE element, except that n_faces=1, and that there is no quarter-point elements
        !

        face=1

        ! =========================== !
        ! ELEMENT-WISE INCIDENT FIELD !
        ! =========================== !

        ! Loop through the ELEMENTS of the BE BODY LOAD
        do ke=1,part(sp)%n_elements
          se=part(sp)%element(ke)
          se_n_nodes=element(se)%n_nodes
          allocate (x_fn(problem%n,se_n_nodes))
          x_fn=element(se)%x_fn
          n_fn=element(se)%n_fn
          select case (region(kr)%type)

            ! --------------
            ! INVISCID FLUID
            ! --------------

            case (fbem_potential)
              ! Loop through INCIDENT FIELDS of the REGION
              do kif=1,region(kr)%n_incidentfields
                sif=region(kr)%incidentfield(kif)
                ! Loop through the NODES of the ELEMENT
                do kn=1,se_n_nodes
                  sn=element(se)%node(kn)
                  select case (incidentfield(sif)%class)
                    case (fbem_point)
                      call fbem_harpot_pointwave(problem%n,omega,rho,c,incidentfield(sif)%amplitude,incidentfield(sif)%x0,&
                                                 incidentfield(sif)%space,incidentfield(sif)%np,incidentfield(sif)%xp,incidentfield(sif)%bc,&
                                                 incidentfield(sif)%symconf,incidentfield(sif)%xs,x_fn(:,kn),n_fn(:,kn),pi,Uni)
                    case (fbem_plane)
                      call fbem_harpot_planewave(problem%n,omega,rho,c,&
                           incidentfield(sif)%amplitude,incidentfield(sif)%x0,incidentfield(sif)%varphi,incidentfield(sif)%theta,&
                           incidentfield(sif)%space,incidentfield(sif)%np,incidentfield(sif)%xp,incidentfield(sif)%bc,&
                           incidentfield(sif)%symconf,incidentfield(sif)%xs,x_fn(:,kn),n_fn(:,kn),pi,Uni)
                  end select
                  ! Depending on the leading variable of the incident wave field
                  select case (incidentfield(sif)%variable)
                    ! Incident wave field in terms of pressure
                      case (0)
                        element(se)%incident_c(1,kn,face)=element(se)%incident_c(1,kn,face)+pi
                        element(se)%incident_c(2,kn,face)=element(se)%incident_c(2,kn,face)+Uni
                    ! Incident wave field in terms of normal displacements
                    case (1)
                      stop 'not implemented yet'
                  end select
                end do
              end do

            ! -------------
            ! ELASTIC SOLID
            ! -------------

            case (fbem_viscoelastic)
              allocate (ui(3),ti(3))
              ! Loop through INCIDENT FIELDS of the REGION
              do kif=1,region(kr)%n_incidentfields
                sif=region(kr)%incidentfield(kif)
                ! Loop through the NODES of the ELEMENT
                do kn=1,se_n_nodes
                  sn=element(se)%node(kn)
                  ! Copy to arrays of dimension 3
                  do kc=1,problem%n
                    xn(kc)=x_fn(kc,kn)
                    nn(kc)=n_fn(kc,kn)
                  end do
                  ! Calculate incident field in terms of displacements
                  select case  (incidentfield(sif)%space)
                    case (fbem_full_space,fbem_half_space)

                    ! mirar esto para nu complejo

                      call fbem_harela_incident_plane_wave(problem%n,k1,k2,mu,lambda,dreal(nu),&
                           incidentfield(sif)%space,incidentfield(sif)%symconf(2),incidentfield(sif)%wave_type,incidentfield(sif)%varphi,&
                           incidentfield(sif)%theta,xn,incidentfield(sif)%xp,nn,ui,ti)
                      ui=0.5d0*ui
                      ti=0.5d0*ti
                    case (fbem_multilayered_half_space)
                      call fbem_harela_incident_plane_wave_vertical_multilayered(problem%n,omega,incidentfield(sif)%np,&
                      incidentfield(sif)%n_layers,incidentfield(sif)%layer_ztop,incidentfield(sif)%layer_lambda,&
                      incidentfield(sif)%layer_mu,incidentfield(sif)%layer_rho,incidentfield(sif)%symconf(2),&
                      incidentfield(sif)%wave_type,incidentfield(sif)%layer_amplitudes,incidentfield(sif)%varphi,xn,nn,ui,ti)
                  end select
                  ! Depending on the leading variable of the incident wave field
                  select case (incidentfield(sif)%variable)
                    ! Incident wave field in terms of displacements
                    case (0)
                      do kc=1,problem%n
                        element(se)%incident_c(          kc,kn,face)=element(se)%incident_c(          kc,kn,face)+ui(kc)
                        element(se)%incident_c(problem%n+kc,kn,face)=element(se)%incident_c(problem%n+kc,kn,face)+ti(kc)
                      end do
                    ! Incident wave field in terms of stresses
                    case (1)
                      if (incidentfield(sif)%space.eq.fbem_multilayered_half_space) stop 'This is not ready yet (multilayered with stress normal.)'
                      select case (incidentfield(sif)%wave_type)
                        ! P-waves
                        case (fbem_harela_p_wave)
                          do kc=1,problem%n
                            element(se)%incident_c(          kc,kn,face)=element(se)%incident_c(          kc,kn,face)+ui(kc)/(-c_im*k1*(lambda+2.d0*mu))
                            element(se)%incident_c(problem%n+kc,kn,face)=element(se)%incident_c(problem%n+kc,kn,face)+ti(kc)/(-c_im*k1*(lambda+2.d0*mu))
                          end do
                        ! S-waves
                        case (fbem_harela_sv_wave,fbem_harela_sh_wave)
                          do kc=1,problem%n
                            element(se)%incident_c(          kc,kn,face)=element(se)%incident_c(          kc,kn,face)+ui(kc)/(-c_im*k2*mu)
                            element(se)%incident_c(problem%n+kc,kn,face)=element(se)%incident_c(problem%n+kc,kn,face)+ti(kc)/(-c_im*k2*mu)
                          end do
                        ! Rayleigh-waves
                        case (fbem_harela_rayleigh_wave)
                          stop 'an incident rayleigh wave in terms of stresses is not implemented yet'
                      end select
                  end select
                end do
              end do
              deallocate (ui,ti)

            ! -----------------
            ! POROELASTIC MEDIUM
            ! -----------------

            case (fbem_poroelastic)
              allocate (ui(0:3),ti(0:3))
              ! Loop through INCIDENT FIELDS of the REGION
              do kif=1,region(kr)%n_incidentfields
                sif=region(kr)%incidentfield(kif)
                ! Loop through the NODES of the ELEMENT
                do kn=1,se_n_nodes
                  sn=element(se)%node(kn)
                  do kc=1,problem%n
                    xn(kc)=x_fn(kc,kn)
                    nn(kc)=n_fn(kc,kn)
                  end do
                  ! Calculate incident field in terms of displacements
                  select case  (incidentfield(sif)%space)
                    case (fbem_full_space,fbem_half_space)
                      select case (incidentfield(sif)%wave_type)
                        case (fbem_harpor_p1_wave,fbem_harpor_p2_wave,fbem_harpor_shx_wave,fbem_harpor_shy_wave)
                          call fbem_harpor_incident_plane_wave(problem%n,omega,k1,k2,k3,lambda,mu,Q,R,rhohat11,rhohat22,rhohat12,&
                               incidentfield(sif)%space,incidentfield(sif)%wave_type,xn,incidentfield(sif)%xp,nn,ui,ti)
                        case (fbem_harpor_R_wave_permeable)
                          if (incidentfield(sif)%space.eq.fbem_full_space) then
                            call fbem_error_message(error_unit,0,'incident field',incidentfield(sif)%id,'Rayleigh incident field is present only in half-space')
                          end if
                          call fbem_harpor_incident_rayleigh_wave(problem%n,lambda,mu,Q,R,phi,rhos,rhof,rhoa,b,omega,xn,nn,ui,ti)
                      end select
                    case (fbem_multilayered_half_space)
                      ! Find kl
                      do kc=incidentfield(sif)%n_layers,1,-1
                        if (element(se)%centroid(problem%n).lt.(incidentfield(sif)%layer_ztop(kc)+1.d-6)) then
                          kl=kc
                          exit
                        end if
                      end do
                      call fbem_harpor_vertical_plane_wave_layered_halfspace(problem%n,incidentfield(sif)%n_layers,&
                      incidentfield(sif)%layer_ztop,incidentfield(sif)%layer_A,incidentfield(sif)%layer_k,incidentfield(sif)%wave_type,1.d-6,kl,xn,nn,ui,ti)
                  end select
                  ! Depending on the leading variable of the incident wave field
                  select case (incidentfield(sif)%variable)
                    ! Incident wave field in terms of displacements
                    case (0)
                      do kc=0,problem%n
                        element(se)%incident_c(            kc,kn,face)=element(se)%incident_c(            kc,kn,face)+ui(kc)
                        element(se)%incident_c(problem%n+1+kc,kn,face)=element(se)%incident_c(problem%n+1+kc,kn,face)+ti(kc)
                      end do
                    ! Incident wave field in terms of total stresses
                    case (1)
                      select case (incidentfield(sif)%wave_type)
                        ! P-waves
                        case (fbem_harpor_p1_wave)
                          do kc=0,problem%n
                            element(se)%incident_c(            kc,kn,face)=element(se)%incident_c(            kc,kn,face)+ui(kc)/(-c_im*k1*D1I*((lambda+2.d0*mu+Q*(1.d0+Q/R))*phi_s_1+(Q+R)*phi_f_1))
                            element(se)%incident_c(problem%n+1+kc,kn,face)=element(se)%incident_c(problem%n+1+kc,kn,face)+ti(kc)/(-c_im*k1*D1I*((lambda+2.d0*mu+Q*(1.d0+Q/R))*phi_s_1+(Q+R)*phi_f_1))
                          end do
                        ! S-waves
                        case (fbem_harpor_shx_wave,fbem_harpor_shy_wave)
                          do kc=0,problem%n
                            element(se)%incident_c(            kc,kn,face)=element(se)%incident_c(            kc,kn,face)+ui(kc)/(-c_im*mu*k3)
                            element(se)%incident_c(problem%n+1+kc,kn,face)=element(se)%incident_c(problem%n+1+kc,kn,face)+ti(kc)/(-c_im*mu*k3)
                          end do
                      end select
                    ! Incident wave field in terms of solid skeleton stresses
                    case (2)
                      select case (incidentfield(sif)%wave_type)
                        ! P-waves
                        case (fbem_harpor_p1_wave)
                          do kc=0,problem%n
                            element(se)%incident_c(            kc,kn,face)=element(se)%incident_c(            kc,kn,face)+ui(kc)/(-c_im*k1*D1I*(phi_s_1*(lambda+2.d0*mu+Q**2/R)+phi_f_1*Q))
                            element(se)%incident_c(problem%n+1+kc,kn,face)=element(se)%incident_c(problem%n+1+kc,kn,face)+ti(kc)/(-c_im*k1*D1I*(phi_s_1*(lambda+2.d0*mu+Q**2/R)+phi_f_1*Q))
                          end do
                        ! S-waves
                        case (fbem_harpor_shx_wave,fbem_harpor_shy_wave)
                          stop 'not implemented yet'
                      end select
                  end select
                end do
              end do
              deallocate (ui,ti)
          end select

          ! Finalize
          deallocate (x_fn,n_fn)

        end do ! Loop through the ELEMENTS of the BE BODY LOAD

        ! ======================== !
        ! NODE-WISE INCIDENT FIELD !
        ! ======================== !

        ! Loop through the NODES of the BE BODY LOAD
        do kn=1,part(sp)%n_nodes
          sn=part(sp)%node(kn)
          do ke=1,node(sn)%n_elements
            se=node(sn)%element(ke)
            knj=node(sn)%element_node_iid(ke)
            node(sn)%incident_c(:,:)=node(sn)%incident_c(:,:)+element(se)%incident_c(:,knj,:)
          end do
          node(sn)%incident_c=node(sn)%incident_c/real(node(sn)%n_elements)
        end do

      end do ! Loop through the BODY LOADS of the REGION

      ! ============================== !
      ! INTERNAL POINTS INCIDENT FIELD !
      ! ============================== !

      select case (region(kr)%type)

        ! --------------
        ! INVISCID FLUID
        ! --------------

        case (fbem_potential)
          do kip=1,region(kr)%n_internalpoints
            ! INTERNAL POINT
            sip=region(kr)%internalpoint(kip)
            x_i=internalpoint(sip)%x
            ! Calculate for each normal direction of the secondary variables
            do kc=1,problem%n
              n_i=0.d0
              n_i(kc)=1.d0
              ! Loop through the incident fields of the region
              do kif=1,region(kr)%n_incidentfields
                ! Selected incident field
                sif=region(kr)%incidentfield(kif)
                select case (incidentfield(sif)%class)
                  case (fbem_point)
                    call fbem_harpot_pointwave(problem%n,omega,rho,c,incidentfield(sif)%amplitude,incidentfield(sif)%x0,&
                                               incidentfield(sif)%space,incidentfield(sif)%np,incidentfield(sif)%xp,incidentfield(sif)%bc,&
                                               incidentfield(sif)%symconf,incidentfield(sif)%xs,x_i,n_i,pi,Uni)
                  case (fbem_plane)
                    call fbem_harpot_planewave(problem%n,omega,rho,c,&
                         incidentfield(sif)%amplitude,incidentfield(sif)%x0,incidentfield(sif)%varphi,incidentfield(sif)%theta,&
                         incidentfield(sif)%space,incidentfield(sif)%np,incidentfield(sif)%xp,incidentfield(sif)%bc,&
                         incidentfield(sif)%symconf,incidentfield(sif)%xs,x_i,n_i,pi,Uni)
                end select
                ! Depending on the leading variable of the incident wave field
                select case (incidentfield(sif)%variable)
                  ! Incident wave field in terms of pressure
                    case (0)
                      if (kc.eq.1) internalpoint(sip)%incident_c(1,0)=internalpoint(sip)%incident_c(1,0)+pi
                      internalpoint(sip)%incident_c(1,kc)=internalpoint(sip)%incident_c(1,kc)+Uni
                  ! Incident wave field in terms of normal displacements
                  case (1)
                    stop 'not implemented yet'
                end select
              end do
            end do
          end do

        ! -------------
        ! ELASTIC SOLID
        ! -------------

        case (fbem_viscoelastic)
          allocate (ui(3),ti(3))
          do kip=1,region(kr)%n_internalpoints
            ! INTERNAL POINT
            sip=region(kr)%internalpoint(kip)
            x_i=internalpoint(sip)%x
            xn=0.d0
            xn(1:problem%n)=x_i(1:problem%n)
            ! Calculate for each normal direction of the secondary variables
            do kc=1,problem%n
              nn=0.d0
              nn(kc)=1.d0
              ! Loop through the incident fields of the region
              do kif=1,region(kr)%n_incidentfields
                ! Selected incident field
                sif=region(kr)%incidentfield(kif)
                ! Calculate incident field in terms of displacements
                select case  (incidentfield(sif)%space)
                  case (fbem_full_space,fbem_half_space)

                  ! mirar esto para nu complejo!!!!

                    call fbem_harela_incident_plane_wave(problem%n,k1,k2,mu,lambda,dreal(nu),&
                         incidentfield(sif)%space,incidentfield(sif)%symconf(2),incidentfield(sif)%wave_type,incidentfield(sif)%varphi,&
                         incidentfield(sif)%theta,xn,incidentfield(sif)%xp,nn,ui,ti)
                  case (fbem_multilayered_half_space)
                    call fbem_harela_incident_plane_wave_vertical_multilayered(problem%n,omega,incidentfield(sif)%np,&
                    incidentfield(sif)%n_layers,incidentfield(sif)%layer_ztop,incidentfield(sif)%layer_lambda,&
                    incidentfield(sif)%layer_mu,incidentfield(sif)%layer_rho,incidentfield(sif)%symconf(2),&
                    incidentfield(sif)%wave_type,incidentfield(sif)%layer_amplitudes,incidentfield(sif)%varphi,xn,nn,ui,ti)
                end select
                ! Depending on the leading variable of the incident wave field
                select case (incidentfield(sif)%variable)
                  ! Incident wave field in terms of displacements
                  case (0)
                    if (kc.eq.1) internalpoint(sip)%incident_c(:,0)=internalpoint(sip)%incident_c(:,0)+ui
                    internalpoint(sip)%incident_c(:,kc)=internalpoint(sip)%incident_c(:,kc)+ti
                  ! Incident wave field in terms of stresses
                  case (1)
                    select case (incidentfield(sif)%wave_type)
                      ! P-waves
                      case (fbem_harela_p_wave)
                        if (kc.eq.1) internalpoint(sip)%incident_c(:,0)=internalpoint(sip)%incident_c(:,0)+ui/(-c_im*k1*(lambda+2.d0*mu))
                        internalpoint(sip)%incident_c(:,kc)=internalpoint(sip)%incident_c(:,kc)+ti/(-c_im*k1*(lambda+2.d0*mu))
                      ! S-waves
                      case (fbem_harela_sv_wave,fbem_harela_sh_wave)
                        if (kc.eq.1) internalpoint(sip)%incident_c(:,0)=internalpoint(sip)%incident_c(:,0)+ui/(-c_im*k2*mu)
                        internalpoint(sip)%incident_c(:,kc)=internalpoint(sip)%incident_c(:,kc)+ti/(-c_im*k2*mu)
                      ! Rayleigh-waves
                      case (fbem_harela_rayleigh_wave)
                        stop 'an incident rayleigh wave in terms of stresses is not implemented yet'
                    end select
                end select
              end do
            end do
          end do
          deallocate (ui,ti)

        ! ------------------
        ! POROELASTIC MEDIUM
        ! ------------------

        case (fbem_poroelastic)
          allocate (ui(0:3),ti(0:3))
          do kip=1,region(kr)%n_internalpoints
            ! INTERNAL POINT
            sip=region(kr)%internalpoint(kip)
            x_i=internalpoint(sip)%x
            xn=0.d0
            xn(1:problem%n)=x_i(1:problem%n)
            ! Calculate for each normal direction of the secondary variables
            do kc=1,problem%n
              nn=0.d0
              nn(kc)=1.d0
              ! Loop through INCIDENT FIELDS of the REGION
              do kif=1,region(kr)%n_incidentfields
                sif=region(kr)%incidentfield(kif)
                ! Calculate incident field in terms of displacements
                select case  (incidentfield(sif)%space)
                  case (fbem_full_space,fbem_half_space)
                    select case (incidentfield(sif)%wave_type)
                      case (fbem_harpor_p1_wave,fbem_harpor_p2_wave,fbem_harpor_shx_wave,fbem_harpor_shy_wave)
                        call fbem_harpor_incident_plane_wave(problem%n,omega,k1,k2,k3,lambda,mu,Q,R,rhohat11,rhohat22,rhohat12,&
                             incidentfield(sif)%space,incidentfield(sif)%wave_type,xn,incidentfield(sif)%xp,nn,ui,ti)
                      case (fbem_harpor_R_wave_permeable)
                        if (incidentfield(sif)%space.eq.fbem_full_space) then
                          call fbem_error_message(error_unit,0,'incident field',incidentfield(sif)%id,'Rayleigh incident field is present only in half-space')
                        end if
                        call fbem_harpor_incident_rayleigh_wave(problem%n,lambda,mu,Q,R,phi,rhos,rhof,rhoa,b,omega,xn,nn,ui,ti)
                    end select
                  case (fbem_multilayered_half_space)
                    ! Find kl
                    do kc2=incidentfield(sif)%n_layers,1,-1
                      if (element(se)%centroid(problem%n).lt.(incidentfield(sif)%layer_ztop(kc2)+1.d-6)) then
                        kl=kc2
                        exit
                      end if
                    end do
                    call fbem_harpor_vertical_plane_wave_layered_halfspace(problem%n,incidentfield(sif)%n_layers,&
                    incidentfield(sif)%layer_ztop,incidentfield(sif)%layer_A,incidentfield(sif)%layer_k,incidentfield(sif)%wave_type,1.d-6,kl,xn,nn,ui,ti)
                end select
                ! Depending on the leading variable of the incident wave field
                select case (incidentfield(sif)%variable)
                  ! Incident wave field in terms of displacements
                  case (0)
                    if (kc.eq.1) internalpoint(sip)%incident_c(:, 0)=internalpoint(sip)%incident_c(:, 0)+ui
                                 internalpoint(sip)%incident_c(:,kc)=internalpoint(sip)%incident_c(:,kc)+ti

                  ! Incident wave field in terms of total stresses
                  case (1)
                    select case (incidentfield(sif)%wave_type)
                      ! P-waves
                      case (fbem_harpor_p1_wave)
                        if (kc.eq.1) internalpoint(sip)%incident_c(:, 0)=internalpoint(sip)%incident_c(:, 0)+ui/(-c_im*k1*D1I*((lambda+2.d0*mu+Q*(1.d0+Q/R))*phi_s_1+(Q+R)*phi_f_1))
                                     internalpoint(sip)%incident_c(:,kc)=internalpoint(sip)%incident_c(:,kc)+ti/(-c_im*k1*D1I*((lambda+2.d0*mu+Q*(1.d0+Q/R))*phi_s_1+(Q+R)*phi_f_1))

                      ! S-waves
                      case (fbem_harpor_shx_wave,fbem_harpor_shy_wave)
                        do kc2=0,problem%n
                          if (kc.eq.1) internalpoint(sip)%incident_c(:, 0)=internalpoint(sip)%incident_c(:, 0)+ui/(-c_im*mu*k3)
                                       internalpoint(sip)%incident_c(:,kc)=internalpoint(sip)%incident_c(:,kc)+ti/(-c_im*mu*k3)
                        end do
                    end select
                  ! Incident wave field in terms of solid skeleton stresses
                  case (2)
                    select case (incidentfield(sif)%wave_type)
                      ! P-waves
                      case (fbem_harpor_p1_wave)
                        do kc2=0,problem%n
                          if (kc.eq.1) internalpoint(sip)%incident_c(:, 0)=internalpoint(sip)%incident_c(:, 0)+ui/(-c_im*k1*D1I*(phi_s_1*(lambda+2.d0*mu+Q**2/R)+phi_f_1*Q))
                                       internalpoint(sip)%incident_c(:,kc)=internalpoint(sip)%incident_c(:,kc)+ti/(-c_im*k1*D1I*(phi_s_1*(lambda+2.d0*mu+Q**2/R)+phi_f_1*Q))
                        end do
                      ! S-waves
                      case (fbem_harpor_shx_wave,fbem_harpor_shy_wave)
                        stop 'not implemented yet'
                    end select
                end select
              end do
            end do
          end do

          deallocate (ui,ti)

      end select

    end if

  end do

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'END calculating incident fields')

end subroutine calculate_incident_mechanics_harmonic
