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


subroutine build_lse_mechanics_harmonic_sa(kf)

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
  use fbem_bem_stapot2d
  use fbem_bem_harpot2d
  use fbem_bem_stapot3d
  use fbem_bem_staela2d
  use fbem_bem_harela2d
  use fbem_bem_staela3d
  use fbem_fem_beams
  use fbem_fem_shells
  use fbem_fem_solids

  ! Module of problem variables
  use problem_variables

  ! No implicit variables
  implicit none

  ! I/O variables
  integer                           :: kf
  ! Local variables
  integer                           :: ka
  integer                           :: sb
  integer                           :: kr
  integer                           :: ks_int, ss_int
  real(kind=real64)                 :: omega
  complex(kind=real64)              :: KT
  integer                           :: ke_int, se_int, se_int_type, se_int_n_nodes
  real(kind=real64), allocatable    :: x_gn_int(:,:)
  integer                           :: kn, sn, kni, knj, sni, snj, kdofei, kdofej
  integer                           :: kc, kci
  real(kind=real64)                 :: rho, nu ! Density and Poisson's ratio
  real(kind=real64)                 :: gamma
  ! Associated with finite elements
  complex(kind=real64)              :: E
  complex(kind=real64)              :: E11, E22, E33, G12, G13, G23       ! Properties of a orthotropic material
  real(kind=real64)                 :: nu12, nu21, nu13, nu31, nu23, nu32 ! Properties of a orthotropic material
  real(kind=real64)                 :: Area, Inertia(3), Volume, t1ref(3)
  real(kind=real64), allocatable    :: Lc(:,:), Lci(:,:), S(:,:), LSp(:,:), Ld(:,:), Qp(:,:), Q(:,:), LcQ(:,:)
  complex(kind=real64), allocatable :: Kp(:,:), K(:,:), Kph(:,:), D(:,:)
  real(kind=real64), allocatable    :: Mp(:,:), M(:,:)
  complex(kind=real64), allocatable :: t_e(:), f_e(:), tp_e(:)
  real(kind=real64), allocatable    :: local_axis(:,:), tv_midnode(:), v_midnode(:,:,:), t_nodes(:)
  complex(kind=real64), allocatable :: dKdx(:,:,:,:)
  real(kind=real64), allocatable    :: dMdx(:,:,:,:)
  integer                           :: se_int_n_dof
  integer, allocatable              :: se_int_n_dof_node(:)
  integer                           :: row, row_local, col, kdofe, gln_kip, gln_ksh, gln_a, gln_sh, kis, kie, kjs, kje
  character(len=fbem_fmtstr)        :: fmtstr            ! String used for write format string

  ! Message
  if (verbose_level.ge.1)  write(output_unit,'(a)') 'Building the b vector for sensitivity analysis ...'

  ! Initialize bsa vector
  bsa_c=(0.0d0,0.0d0)

  ! Frequency
  omega=frequency(kf)

  ! ===========================
  ! BUILD BEM AND FEM EQUATIONS
  ! ===========================

  ! FEM variables allocate
  allocate (local_axis(problem%n,problem%n))

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
        select case (region(kr)%type)
          ! INVISCID FLUID
          case (fbem_potential)
            call build_lse_mechanics_bem_harpot_sa(kf,kr)
          ! VISCOELASTIC SOLID
          case (fbem_viscoelastic)
            call build_lse_mechanics_bem_harela_sa(kf,kr)
          ! VISCOPOROELASTIC MEDIUM
          case (fbem_poroelastic)
            ! nothing to do yet
        end select
      ! ============================================================================================================================

      ! ============================================================================================================================
      ! FE REGION
      !
      case (fbem_fe)

        ! Use material of the FE REGION
        if (region(kr)%subtype.eq.fbem_elastic_orthotropic) then
          ! ORTHOTROPIC
          if (problem%n.ne.3) stop 'orthotropic material only 3D'
          E11 =region(kr)%property_c( 1)
          E22 =region(kr)%property_c( 2)
          E33 =region(kr)%property_c( 3)
          G12 =region(kr)%property_c( 4)
          G13 =region(kr)%property_c( 5)
          G23 =region(kr)%property_c( 6)
          nu12=region(kr)%property_r( 7)
          nu21=region(kr)%property_r( 8)
          nu13=region(kr)%property_r( 9)
          nu31=region(kr)%property_r(10)
          nu23=region(kr)%property_r(11)
          nu32=region(kr)%property_r(12)
          rho =region(kr)%property_r(13)
        else
          ! ISOTROPIC
          rho=region(kr)%property_r(1)
          E=region(kr)%property_c(5)
          nu=region(kr)%property_r(3)
        end if

        ! Loop through the FE SUBREGIONS of the FE REGION
        do ks_int=1,region(kr)%n_fe_subregions
          ss_int=region(kr)%fe_subregion(ks_int)

          ! Loop through the FINITE ELEMENTS of the FE SUBREGION
          do ke_int=1,part(fe_subregion(ss_int)%part)%n_elements

            ! --------------
            ! ELEMENT SET-UP
            ! --------------

            ! Integration element
            se_int=part(fe_subregion(ss_int)%part)%element(ke_int)
            ! Assemble only if dxda is not zero for some node, design variable and coordinate.
            if (.not.element(se_int)%sensitivity) cycle
            ! Integration element type
            se_int_type=element(se_int)%type
            ! Integration element number of nodes
            se_int_n_nodes=element(se_int)%n_nodes
            ! Coordinates of each (geometrical) node of the element
            allocate (x_gn_int(problem%n,se_int_n_nodes))
            x_gn_int=element(se_int)%x_gn
            ! DOF of the element
            allocate (se_int_n_dof_node(se_int_n_nodes))
            do kn=1,se_int_n_nodes
              sn=element(se_int)%node(kn)
              se_int_n_dof_node(kn)=node(sn)%n_dof
            end do
            se_int_n_dof=sum(se_int_n_dof_node)

            !  Message
            if (verbose_level.ge.2) then
              write(fmtstr,*) '(3x,a,i',fbem_nchar_int(element(se_int)%id),')'
              call fbem_trimall(fmtstr)
              write(output_unit,fmtstr) 'Integration element: ', element(se_int)%id
            end if

            select case (problem%n)

              ! -------------------------------------------------------------------------------------------------------------------- !
              ! 2D                                                                                                                   !
              ! -------------------------------------------------------------------------------------------------------------------- !

              case (2)
                select case (element(se_int)%n_dimension)

                  case (1)
                    ! Not yet

                  case (2)
                    ! Geometrical data
                    allocate (t_nodes(se_int_n_nodes))
                    select case (problem%subtype)
                      case (fbem_mechanics_plane_strain)
                        t_nodes=1.d0
                      case (fbem_mechanics_plane_stress)
                        t_nodes=element(se_int)%tn_midnode(3,:)
                    end select
                    ! Calculate D
                    allocate (D(3,3))
                    select case (problem%subtype)
                      case (fbem_mechanics_plane_strain)
                        call fbem_fem_harela2d_solid_D_isotropic(E,nu,1,D)
                      case (fbem_mechanics_plane_stress)
                        call fbem_fem_harela2d_solid_D_isotropic(E,nu,2,D)
                    end select
                    ! Integration setup
                    ! Full integration (FI)
                    select case (se_int_type)
                      case (fbem_tri3,fbem_quad4)
                        gln_a=1
                        gln_sh=1
                      case (fbem_tri6,fbem_quad8,fbem_quad9)
                        gln_a=2
                        gln_sh=2
                    end select
!                    ! Reduced integration (RI)
!                    select case (se_int_type)
!                      case (fbem_tri3,fbem_quad4)
!                        gln_a=1
!                        gln_sh=1
!                      case (fbem_tri6,fbem_quad8,fbem_quad9)
!                        gln_a=1
!                        gln_sh=1
!                    end select
!                    ! Selective integration
!                    select case (se_int_type)
!                      case (fbem_tri3,fbem_quad4)
!                        gln_a=1
!                        gln_sh=1
!                      case (fbem_tri6,fbem_quad8,fbem_quad9)
!                        gln_a=2
!                        gln_sh=1
!                    end select
                    ! Calculate dKdx
                    allocate (dKdx(2*se_int_n_nodes,2*se_int_n_nodes,2,se_int_n_nodes))
                    call fbem_fem_harela2d_solid_dKdx(se_int_type,x_gn_int,t_nodes,D,gln_a,gln_sh,dKdx)
                    ! Calculate dMdx
                    allocate (dMdx(2*se_int_n_nodes,2*se_int_n_nodes,2,se_int_n_nodes))
                    call fbem_fem_ela2d_solid_dMdx(se_int_type,x_gn_int,t_nodes,rho,gln_a,dMdx)
                    ! Assemble dKdx-omega**2*dMdx
                    dKdx=dKdx-omega**2*dMdx
                    call assemble_fem_harela_dKdx(se_int,se_int_n_dof,dKdx)
                    ! Assemble dfdx
                    call assemble_fem_harela_dfdx(se_int)
                    deallocate (t_nodes,D,dKdx,dMdx)

                end select

              ! -------------------------------------------------------------------------------------------------------------------- !
              ! 3D                                                                                                                   !
              ! -------------------------------------------------------------------------------------------------------------------- !

              case (3)
                select case (element(se_int)%n_dimension)

                  case (1)
                    stop 'not yet 76'

                  case (2)
                    ! Nothing

                  ! SOLID
                  case (3)
                    stop 'not yet 77'
                end select

            end select

            deallocate (x_gn_int,se_int_n_dof_node)

          end do ! Loop through the FINITE ELEMENTS of the FE SUBREGION

        end do ! Loop through the FE SUBREGIONS of the FE REGION

      ! ============================================================================================================================

    end select ! Switch between REGION CLASSES

  end do ! Loop through REGIONS

  ! Ending message
  if (verbose_level.ge.1) write(output_unit,'(a)') 'done.'

end subroutine build_lse_mechanics_harmonic_sa
