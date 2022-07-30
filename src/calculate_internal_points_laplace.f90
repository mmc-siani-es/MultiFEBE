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


subroutine calculate_internal_points_laplace

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
  integer                        :: kn_int, sn_int
  integer                        :: kip, sip
  integer                        :: kn
  integer                        :: kc
  real(kind=real64)              :: x_i(problem%n), n_i(problem%n)
  ! Region properties
  real(kind=real64)              :: conductivity
  ! Vectors for SBIE integration
  real(kind=real64), allocatable :: h(:), g(:)
  real(kind=real64), allocatable :: hp(:), gp(:)
  real(kind=real64), allocatable :: hm(:), gm(:)
  ! Vectors for HBIE integration
  real(kind=real64), allocatable :: m(:), l(:)
  real(kind=real64), allocatable :: mp(:), lp(:)
  real(kind=real64), allocatable :: mm(:), lm(:)
  ! Associated with symmetry
  real(kind=real64), allocatable :: symconf_m(:), symconf_t(:), symconf_r(:)
  real(kind=real64)              :: symconf_s
  logical                        :: reversed
  ! Writing
  character(len=fbem_fmtstr)     :: fmtstr

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'START calculating internal points solutions at BE regions')

  ! Initialize internal points values
  do kip=1,n_internalpoints
    internalpoint(kip)%value_r=0.d0
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

          !
          ! Loop through INTERNAL POINTS
          !
          do kip=1,region(kr)%n_internalpoints
            ! Internal point iid
            sip=region(kr)%internalpoint(kip)
            ! Save x_i
            x_i=internalpoint(sip)%x

            ! -----------
            ! CALCULATE P
            ! -----------

            ! CALCULATE KERNELS
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
            ! ASSEMBLE
            ! The flux variable is j=k·dp/dn
            gp=gp/conductivity
            gm=gm/conductivity
            select case (boundary(sb_int)%coupling)
              ! BE OR BE-FE BOUNDARY
              case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
                select case (boundary(sb_int)%class)
                  ! ORDINARY BOUNDARY
                  case (fbem_boundary_class_ordinary)
                    do kn_int=1,se_int_n_nodes
                      sn_int=element(se_int)%node(kn_int)
                      internalpoint(sip)%value_r(1,0)=internalpoint(sip)%value_r(1,0)&
                                                     +gp(kn_int)*node(sn_int)%value_r(2,1)-hp(kn_int)*node(sn_int)%value_r(1,1)
                    end do
                  ! CRACK-LIKE BOUNDARY
                  case (fbem_boundary_class_cracklike)
                    do kn_int=1,se_int_n_nodes
                      sn_int=element(se_int)%node(kn_int)
                      internalpoint(sip)%value_r(1,0)=internalpoint(sip)%value_r(1,0)&
                                                     +gp(kn_int)*node(sn_int)%value_r(2,1)-hp(kn_int)*node(sn_int)%value_r(1,1)&
                                                     +gm(kn_int)*node(sn_int)%value_r(2,2)-hm(kn_int)*node(sn_int)%value_r(1,2)
                    end do
                end select
              ! BE-BE OR BE-FE-BE BOUNDARY
              case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                if (sb_int_reversion.eqv.(.false.)) then
                  do kn_int=1,se_int_n_nodes
                    sn_int=element(se_int)%node(kn_int)
                    internalpoint(sip)%value_r(1,0)=internalpoint(sip)%value_r(1,0)&
                                                   +gp(kn_int)*node(sn_int)%value_r(2,1)-hp(kn_int)*node(sn_int)%value_r(1,1)
                  end do
                else
                  do kn_int=1,se_int_n_nodes
                    sn_int=element(se_int)%node(kn_int)
                    internalpoint(sip)%value_r(1,0)=internalpoint(sip)%value_r(1,0)&
                                                   +gp(kn_int)*node(sn_int)%value_r(2,2)-hp(kn_int)*node(sn_int)%value_r(1,2)
                  end do
                end if
            end select

            ! ------------------------
            ! CALCULATE J WITH N=E_{K}
            ! ------------------------

            ! Loop through AXES
            do kc=1,problem%n
              ! UNIT NORMAL
              n_i=0.d0
              n_i(kc)=1.d0
              ! CALCULATE KERNELS
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
              ! ASSEMBLE
              ! The flux variable is j=k·dp/dn
              mp=mp*conductivity
              mm=mm*conductivity
              select case (boundary(sb_int)%coupling)
                ! BE OR BE-FE BOUNDARY
                case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
                  select case (boundary(sb_int)%class)
                    ! ORDINARY BOUNDARY
                    case (fbem_boundary_class_ordinary)
                      do kn_int=1,se_int_n_nodes
                        sn_int=element(se_int)%node(kn_int)
                        internalpoint(sip)%value_r(1,kc)=internalpoint(sip)%value_r(1,kc)&
                                                        +lp(kn_int)*node(sn_int)%value_r(2,1)-mp(kn_int)*node(sn_int)%value_r(1,1)
                      end do
                    ! CRACK-LIKE BOUNDARY
                    case (fbem_boundary_class_cracklike)
                      do kn_int=1,se_int_n_nodes
                        sn_int=element(se_int)%node(kn_int)
                        internalpoint(sip)%value_r(1,kc)=internalpoint(sip)%value_r(1,kc)&
                                                        +lp(kn_int)*node(sn_int)%value_r(2,1)-mp(kn_int)*node(sn_int)%value_r(1,1)&
                                                        +lm(kn_int)*node(sn_int)%value_r(2,2)-mm(kn_int)*node(sn_int)%value_r(1,2)
                      end do
                  end select
                ! BE-BE OR BE-FE-BE BOUNDARY
                case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                  if (sb_int_reversion.eqv.(.false.)) then
                    do kn_int=1,se_int_n_nodes
                      sn_int=element(se_int)%node(kn_int)
                      internalpoint(sip)%value_r(1,kc)=internalpoint(sip)%value_r(1,kc)&
                                                      +lp(kn_int)*node(sn_int)%value_r(2,1)-mp(kn_int)*node(sn_int)%value_r(1,1)
                    end do
                  else
                    do kn_int=1,se_int_n_nodes
                      sn_int=element(se_int)%node(kn_int)
                      internalpoint(sip)%value_r(1,kc)=internalpoint(sip)%value_r(1,kc)&
                                                      +lp(kn_int)*node(sn_int)%value_r(2,2)-mp(kn_int)*node(sn_int)%value_r(1,2)
                    end do
                  end if
              end select
            end do ! Loop through AXES

          end do ! Loop through INTERNAL POINTS

        end do ! Loop through SYMMETRICAL ELEMENTS

        ! Deallocate element-wise data structures
        deallocate (x_gn_int,xi_gn_int,x_fn_int,n_fn_int)
        deallocate (h,g,hp,gp,hm,gm)
        deallocate (m,l,mp,lp,mm,lm)

      end do ! Loop through the ELEMENTS of the BOUNDARY for COLLOCATION

    end do ! Loop through the BOUNDARIES of the REGION for INTEGRATION

  end do ! Loop through REGIONS

if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'END calculating internal points solutions at BE regions')

end subroutine calculate_internal_points_laplace
