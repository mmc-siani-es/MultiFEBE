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

subroutine calculate_stresses_mechanics_harmonic(kf)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_numerical
  use fbem_shape_functions
  use fbem_geometry
  use fbem_symmetry
  use fbem_fem_beams
  use fbem_fem_shells
  use fbem_fem_solids

  ! Module of problem variables
  use problem_variables

  ! No implicit variables
  implicit none

  ! I/O
  integer                           :: kf
  ! Local variables
  real(kind=real64)                 :: omega
  integer                           :: kr
  integer                           :: sp
  integer                           :: ks_int, ss_int
  integer                           :: ke, ke_int, se, se_int
  integer                           :: kb, sb
  integer                           :: kn, sn, kn2, sn2
  integer                           :: ki, kj, kk
  integer                           :: face

  integer                           :: se_int_type
  integer                           :: ndof_mat
  integer                           :: kc, kdofe

  integer                           :: kdof
  real(kind=real64), allocatable    :: xi(:), phi(:), psi(:,:)
  real(kind=real64)                 :: aux(10)

  real(kind=real64)                 :: par_J

  complex(kind=real64)              :: E, lambda, mu, nu_c
  real(kind=real64)                 :: nu, rho
  complex(kind=real64)              :: E11, E22, E33, G12, G13, G23
  real(kind=real64)                 :: nu12, nu21, nu13, nu31, nu23, nu32

  complex(kind=real64)              :: R, Q
  real(kind=real64)                 :: rhos, rhof, phi_p, rho1, rho2, rhoa, b
  complex(kind=real64)              :: rhohat11, rhohat12, rhohat22
  complex(kind=real64)              :: pJ, pZ

  complex(kind=real64), allocatable :: a(:), ap(:), F(:,:), K(:,:), f_e(:), KLT(:,:)
  real(kind=real64), allocatable    :: Lc(:,:), x1ref(:)
  real(kind=real64), allocatable    :: Fsigma(:,:,:), nodal_axes(:,:,:)
  integer                           :: setype
  real(kind=real64)                 :: sedelta

  real(kind=real64)                 :: xi1d, xi2d(2), xi3d(3)

  complex(kind=real64)              :: dpdx_tangential(problem%n)
  complex(kind=real64)              :: dpdx_tangential_local(problem%n)
  complex(kind=real64)              :: U_normal(problem%n)
  complex(kind=real64)              :: U_local(problem%n), U(problem%n)
  complex(kind=real64)              :: dudx_tangential(problem%n,problem%n)
  complex(kind=real64)              :: dudx_tangential_local(problem%n,problem%n)
  complex(kind=real64)              :: sigma(3,3)
  complex(kind=real64)              :: sigma_local(3,3)
  real(kind=real64)                 :: delta
  real(kind=real64)                 :: local_axes(problem%n,problem%n)

  ! TODO: modify message

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'START calculating stresses at FE regions')

  omega=frequency(kf)

  ! ================================================================================================================================
  ! CALCULATE BEM SOLID STRESSES / FLUID DISPLACEMENTS ALONG BOUNDARIES
  ! CALCULATE FEM ELEMENT STRESSES / STRESS RESULTANTS / EQUILIBRATING LOADS AND NODAL REACTIONS
  ! ================================================================================================================================

  !
  ! Initialize
  !
  do kr=1,n_regions
    !
    ! Boundary elements
    !
    if (region(kr)%class.eq.fbem_be) then
      do kb=1,region(kr)%n_boundaries
        sb=region(kr)%boundary(kb)
        sp=boundary(sb)%part
        do ke=1,part(sp)%n_elements
          se=part(sp)%element(ke)
          element(se)%value_c=0
        end do
      end do
    end if
    !
    ! Finite elements
    !
    if (region(kr)%class.eq.fbem_fe) then
      do ks_int=1,region(kr)%n_fe_subregions
        ss_int=region(kr)%fe_subregion(ks_int)
        do ke_int=1,part(fe_subregion(ss_int)%part)%n_elements
          se_int=part(fe_subregion(ss_int)%part)%element(ke_int)
          element(se_int)%value_c(:,:,2)=0
          do kn=1,element(se_int)%n_nodes
            sn=element(se_int)%node(kn)
            node(sn)%value_c(:,2)=0
          end do
        end do
      end do
    end if
  end do

  ! ================================================================================================================================
  ! CALCULATE BEM SOLID STRESSES / FLUID DISPLACEMENTS ALONG BOUNDARIES
  ! ================================================================================================================================
  ! Stress recovery from displacement gradient taken from pages 203-204 the reference book of Telles, Brebbia and Wrobel (1984).

  ! Loop through REGIONS
  do kr=1,n_regions
    !
    ! BE region
    !
    if (region(kr)%class.eq.fbem_be) then

        ! Save material properties to local variables
        select case (region(kr)%type)

          case (fbem_potential)
            par_J=1.d0/(region(kr)%property_r(1)*omega**2)

          case (fbem_viscoelastic)
            if (region(kr)%space.eq.fbem_multilayered_half_space) then
              write(*,*) 'Warning: BEM boundary stress tensor, not multilayered halfspace for elastodynamics yet'
              cycle
            else
              select case (region(kr)%subtype)
                !
                ! Elastic or viscoelastic
                !
                case (0,fbem_elastic)
                  nu_c=region(kr)%property_c(3)
                  mu=region(kr)%property_c(2)
                  lambda=region(kr)%property_c(6)
                case default
                  write(*,*) 'Warning: BEM boundary stress tensor, only for elastic material (not for Bardet, etc...)'
                  cycle
              end select
            end if

          case (fbem_poroelastic)
            rhof  =region(kr)%property_r(1)
            rhos  =region(kr)%property_r(2)
            lambda=region(kr)%property_c(3)
            mu    =region(kr)%property_c(4)
            nu_c  =region(kr)%property_c(6)
            phi_p =region(kr)%property_r(8)
            rhoa  =region(kr)%property_r(9)
            rho1  =region(kr)%property_r(13)
            rho2  =region(kr)%property_r(14)
            R     =region(kr)%property_c(10)
            Q     =region(kr)%property_c(11)
            b     =region(kr)%property_r(12)
            rhohat11=rho1+rhoa-c_im*b/omega
            rhohat12=    -rhoa+c_im*b/omega
            rhohat22=rho2+rhoa-c_im*b/omega
            pJ=1.0d0/((omega**2)*rhohat22)
            pZ=rhohat12/rhohat22
        end select

        ! Loop through boundaries
        do kb=1,region(kr)%n_boundaries
          sb=region(kr)%boundary(kb)
          sp=boundary(sb)%part
          ! Loop through elements
          do ke=1,part(sp)%n_elements
            se=part(sp)%element(ke)
            allocate(psi(element(se)%n_nodes,problem%n))
            allocate(xi(element(se)%n_dimension))
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
                    select case (region(kr)%type)
                      !
                      ! INVISCID FLUID
                      !
                      case (fbem_potential)

                        face=1
                        do kn=1,element(se)%n_nodes
                          sn=element(se)%node(kn)
                          xi=element(se)%xi_fn(:,kn)
                          psi=fbem_psi(problem%n,element(se)%type_g,element(se)%x_gn,element(se)%type_f1,element(se)%delta_f,xi)
                          dpdx_tangential=0
                          do kn2=1,element(se)%n_nodes
                            sn2=element(se)%node(kn2)
                            dpdx_tangential = dpdx_tangential + psi(kn2,:)*node(sn2)%value_c(1,face)
                          end do
                          U_normal=node(sn)%value_c(2,face)*element(se)%n_gn(:,kn)
                          if (region(kr)%boundary_reversion(kb)) U_normal=-U_normal
                          element(se)%value_c(:,kn,face)=U_normal+dpdx_tangential*par_J
                        end do

                      !
                      ! VISCOELASTIC SOLID
                      !
                      case (fbem_viscoelastic)

                        face=1
                        do kn=1,element(se)%n_nodes
                          sn=element(se)%node(kn)
                          xi=element(se)%xi_fn(:,kn)
                          psi=fbem_psi(problem%n,element(se)%type_g,element(se)%x_gn,element(se)%type_f1,element(se)%delta_f,xi)
                          ! GRADIENT OF DISPLACEMENTS (TANGENTIAL PART)
                          dudx_tangential=0
                          do kn2=1,element(se)%n_nodes
                            sn2=element(se)%node(kn2)
                            do kc=1,problem%n
                              dudx_tangential(kc,:)=dudx_tangential(kc,:)+psi(kn2,:)*node(sn2)%value_c(kc,face)
                            end do
                          end do
                          ! LOCAL CARTESIAN AXES
                          local_axes=fbem_local_tangential_cartesian_axes(problem%n,element(se)%type_g,element(se)%x_gn,region(kr)%boundary_reversion(kb),xi)
                          ! CALCULATE STRESS TENSOR
                          select case (problem%n)
                            !
                            ! 2D (plane strain or plane stress)
                            !
                            case (2)
                              sigma_local=0
                              sigma=0
                              ! Local stress tensor due to traction at boundary
                              sigma_local(2,1)=dot_product(node(sn)%value_c(3:4,face),local_axes(:,1))
                              sigma_local(2,2)=dot_product(node(sn)%value_c(3:4,face),local_axes(:,2))
                              sigma_local(1,2)=sigma_local(2,1)
                              ! Tangential displacement gradient in local cartesian coordinates
                              dudx_tangential_local=matmul(matmul(transpose(local_axes),dudx_tangential),local_axes)
                              ! Local stress tensor due to traction at boundary and tangential displacement
                              select case (problem%subtype)
                                case (fbem_mechanics_plane_strain)
                                  sigma_local(1,1)=1.d0/(1.d0-nu_c)*(2.d0*mu*dudx_tangential_local(1,1)+nu_c*sigma_local(2,2))
                                case (fbem_mechanics_plane_stress)
                                  sigma_local(1,1)=2.d0*(1.d0+nu_c)*mu*dudx_tangential_local(1,1)+nu_c*sigma_local(2,2)
                                case default
                                  call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid 2D subtype')
                              end select
                              ! Stress tensor in global coordinates
                              sigma(1:2,1:2)=matmul(matmul(local_axes,sigma_local(1:2,1:2)),transpose(local_axes))
                              if (problem%subtype.eq.fbem_mechanics_plane_strain) then
                                sigma(3,3)=nu_c*(sigma(1,1)+sigma(2,2))
                              end if
                            !
                            ! 3D
                            !
                            case (3)
                              ! Local stress tensor due to traction at boundary
                              sigma_local(3,1)=dot_product(node(sn)%value_c(4:6,face),local_axes(:,1))
                              sigma_local(3,2)=dot_product(node(sn)%value_c(4:6,face),local_axes(:,2))
                              sigma_local(3,3)=dot_product(node(sn)%value_c(4:6,face),local_axes(:,3))
                              sigma_local(1,3)=sigma_local(3,1)
                              sigma_local(2,3)=sigma_local(3,2)
                              ! Tangential displacement gradient in local cartesian coordinates
                              dudx_tangential_local=matmul(matmul(transpose(local_axes),dudx_tangential),local_axes)
                              ! Local stress tensor due to traction at boundary and tangential displacement
                              sigma_local(1,1)=1.d0/(1.d0-nu_c)*(nu_c*sigma_local(3,3)+2.d0*mu*(dudx_tangential_local(1,1)+nu_c*dudx_tangential_local(2,2)))
                              sigma_local(2,2)=1.d0/(1.d0-nu_c)*(nu_c*sigma_local(3,3)+2.d0*mu*(dudx_tangential_local(2,2)+nu_c*dudx_tangential_local(1,1)))
                              sigma_local(1,2)=mu*(dudx_tangential_local(1,2)+dudx_tangential_local(2,1))
                              sigma_local(2,1)=sigma_local(1,2)
                              ! Stress tensor in global coordinates
                              sigma=matmul(matmul(local_axes,sigma_local),transpose(local_axes))
                            case default
                              call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid dimension')
                          end select
                          ! SAVE STRESS TENSOR
                          do ki=1,3
                            do kj=1,3
                              kk=(ki-1)*3+kj
                              element(se)%value_c(kk,kn,face)=sigma(ki,kj)
                            end do
                          end do
                        end do

                      ! POROELASTIC MEDIA
                      case (fbem_poroelastic)

                        face=1
                        do kn=1,element(se)%n_nodes
                          sn=element(se)%node(kn)
                          xi=element(se)%xi_fn(:,kn)
                          psi=fbem_psi(problem%n,element(se)%type_g,element(se)%x_gn,element(se)%type_f1,element(se)%delta_f,xi)
                          !
                          ! CALCULATE FLUID DISPLACEMENTS
                          !
                          ! GRADIENT OF FLUID EQUIVALENT PRESSURE (TAU)
                          dpdx_tangential=0
                          do kn2=1,element(se)%n_nodes
                            sn2=element(se)%node(kn2)
                            dpdx_tangential = dpdx_tangential + psi(kn2,:)*node(sn2)%value_c(0,face)
                          end do
                          ! ADD
                          U_normal=node(sn)%value_c(problem%n+1,face)*element(se)%n_gn(:,kn)
                          if (region(kr)%boundary_reversion(kb)) U_normal=-U_normal
                          ! SAVE
                          element(se)%value_c(1:problem%n,kn,face)=U_normal-dpdx_tangential*pJ-pZ*node(sn2)%value_c(1:problem%n,face)
                          !
                          ! CALCULATE SOLID STRESS TENSOR
                          !
                          ! GRADIENT OF SOLID DISPLACEMENTS (TANGENTIAL PART)
                          dudx_tangential=0
                          do kn2=1,element(se)%n_nodes
                            sn2=element(se)%node(kn2)
                            do kc=1,problem%n
                              dudx_tangential(kc,:)=dudx_tangential(kc,:)+psi(kn2,:)*node(sn2)%value_c(kc,face)
                            end do
                          end do
                          ! LOCAL CARTESIAN AXES
                          local_axes=fbem_local_tangential_cartesian_axes(problem%n,element(se)%type_g,element(se)%x_gn,region(kr)%boundary_reversion(kb),xi)
                          ! CALCULATE SOLID STRESS TENSOR
                          select case (problem%n)
                            !
                            ! 2D (plane strain)
                            !
                            case (2)
                              sigma_local=0
                              sigma=0
                              ! Local stress tensor due to traction at boundary
                              sigma_local(2,1)=dot_product(node(sn)%value_c(4:5,face),local_axes(:,1))
                              sigma_local(2,2)=dot_product(node(sn)%value_c(4:5,face),local_axes(:,2))
                              sigma_local(1,2)=sigma_local(2,1)
                              ! Tangential displacement gradient in local cartesian coordinates
                              dudx_tangential_local=matmul(matmul(transpose(local_axes),dudx_tangential),local_axes)
                              ! Local stress tensor due to traction at boundary and tangential displacement
                              select case (problem%subtype)
                                case (fbem_mechanics_plane_strain)
                                  ! Ojo: falta chequear si la nu usada es la correcta para def. plana.
                                  sigma_local(1,1)=1.d0/(1.d0-nu_c)*(2.d0*mu*dudx_tangential_local(1,1)+nu_c*sigma_local(2,2))
                                case default
                                  call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid 2D subtype (plane stress), not available for poroelastodynamics.')
                              end select
                              ! Stress tensor in global coordinates
                              sigma(1:2,1:2)=matmul(matmul(local_axes,sigma_local(1:2,1:2)),transpose(local_axes))
                              ! Ojo: falta chequear esto para def. plana.
                              !if (problem%subtype.eq.fbem_mechanics_plane_strain) then
                              !  sigma(3,3)=nu_c*(sigma(1,1)+sigma(2,2))
                              !end if
                              sigma(1,1)=sigma(1,1)+Q/R*node(sn)%value_c(0,face)
                              sigma(2,2)=sigma(2,2)+Q/R*node(sn)%value_c(0,face)
                            !
                            ! 3D
                            !
                            case (3)
                              ! Local stress tensor due to traction at boundary
                              sigma_local(3,1)=dot_product(node(sn)%value_c(5:7,face),local_axes(:,1))
                              sigma_local(3,2)=dot_product(node(sn)%value_c(5:7,face),local_axes(:,2))
                              sigma_local(3,3)=dot_product(node(sn)%value_c(5:7,face),local_axes(:,3))
                              sigma_local(1,3)=sigma_local(3,1)
                              sigma_local(2,3)=sigma_local(3,2)
                              ! Tangential displacement gradient in local cartesian coordinates
                              dudx_tangential_local=matmul(matmul(transpose(local_axes),dudx_tangential),local_axes)
                              ! Local stress tensor due to traction at boundary and tangential displacement
                              sigma_local(1,1)=1.d0/(1.d0-nu_c)*(nu_c*sigma_local(3,3)+2.d0*mu*(dudx_tangential_local(1,1)+nu_c*dudx_tangential_local(2,2)))
                              sigma_local(2,2)=1.d0/(1.d0-nu_c)*(nu_c*sigma_local(3,3)+2.d0*mu*(dudx_tangential_local(2,2)+nu_c*dudx_tangential_local(1,1)))
                              sigma_local(1,2)=mu*(dudx_tangential_local(1,2)+dudx_tangential_local(2,1))
                              sigma_local(2,1)=sigma_local(1,2)
                              ! Stress tensor in global coordinates
                              sigma=matmul(matmul(local_axes,sigma_local),transpose(local_axes))
                              sigma(1,1)=sigma(1,1)+Q/R*node(sn)%value_c(0,face)
                              sigma(2,2)=sigma(2,2)+Q/R*node(sn)%value_c(0,face)
                              sigma(3,3)=sigma(3,3)+Q/R*node(sn)%value_c(0,face)
                            case default
                              call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid dimension')
                          end select
                          ! SAVE
                          do ki=1,3
                            do kj=1,3
                              kk=problem%n+(ki-1)*3+kj
                              element(se)%value_c(kk,kn,face)=sigma(ki,kj)
                            end do
                          end do
                        end do

                    end select

                  !
                  ! Crack-like boundaries
                  !
                  case (fbem_boundary_class_cracklike)

                    ! Similar to BE boundary or BE-FE coupled boundary, except that face=2 needs normal reversion with respecto to the boundary.

                    select case (region(kr)%type)
                      !
                      ! INVISCID FLUID
                      !
                      case (fbem_potential)

                        do face=1,2
                          do kn=1,element(se)%n_nodes
                            sn=element(se)%node(kn)
                            xi=element(se)%xi_fn(:,kn)
                            psi=fbem_psi(problem%n,element(se)%type_g,element(se)%x_gn,element(se)%type_f1,element(se)%delta_f,xi)
                            dpdx_tangential=0
                            do kn2=1,element(se)%n_nodes
                              sn2=element(se)%node(kn2)
                              dpdx_tangential = dpdx_tangential + psi(kn2,:)*node(sn2)%value_c(1,face)
                            end do
                            U_normal=node(sn)%value_c(2,face)*element(se)%n_gn(:,kn)
                            if (face.eq.1) then
                              if (region(kr)%boundary_reversion(kb)) U_normal=-U_normal
                            else
                              if (.not.region(kr)%boundary_reversion(kb)) U_normal=-U_normal
                            end if
                            element(se)%value_c(:,kn,face)=U_normal+dpdx_tangential*par_J
                          end do
                        end do

                      !
                      ! VISCOELASTIC SOLID
                      !
                      case (fbem_viscoelastic)

                        do face=1,2
                          do kn=1,element(se)%n_nodes
                            sn=element(se)%node(kn)
                            xi=element(se)%xi_fn(:,kn)
                            psi=fbem_psi(problem%n,element(se)%type_g,element(se)%x_gn,element(se)%type_f1,element(se)%delta_f,xi)
                            ! GRADIENT OF DISPLACEMENTS (TANGENTIAL PART)
                            dudx_tangential=0
                            do kn2=1,element(se)%n_nodes
                              sn2=element(se)%node(kn2)
                              do kc=1,problem%n
                                dudx_tangential(kc,:)=dudx_tangential(kc,:)+psi(kn2,:)*node(sn2)%value_c(kc,face)
                              end do
                            end do
                            ! LOCAL CARTESIAN AXES
                            if (face.eq.1) then
                              local_axes=fbem_local_tangential_cartesian_axes(problem%n,element(se)%type_g,element(se)%x_gn,region(kr)%boundary_reversion(kb),xi)
                            else
                              local_axes=fbem_local_tangential_cartesian_axes(problem%n,element(se)%type_g,element(se)%x_gn,.not.region(kr)%boundary_reversion(kb),xi)
                            end if
                            ! CALCULATE STRESS TENSOR
                            select case (problem%n)
                              !
                              ! 2D (plane strain or plane stress)
                              !
                              case (2)
                                sigma_local=0
                                sigma=0
                                ! Local stress tensor due to traction at boundary
                                sigma_local(2,1)=dot_product(node(sn)%value_c(3:4,face),local_axes(:,1))
                                sigma_local(2,2)=dot_product(node(sn)%value_c(3:4,face),local_axes(:,2))
                                sigma_local(1,2)=sigma_local(2,1)
                                ! Tangential displacement gradient in local cartesian coordinates
                                dudx_tangential_local=matmul(matmul(transpose(local_axes),dudx_tangential),local_axes)
                                ! Local stress tensor due to traction at boundary and tangential displacement
                                select case (problem%subtype)
                                  case (fbem_mechanics_plane_strain)
                                    sigma_local(1,1)=1.d0/(1.d0-nu_c)*(2.d0*mu*dudx_tangential_local(1,1)+nu_c*sigma_local(2,2))
                                  case (fbem_mechanics_plane_stress)
                                    sigma_local(1,1)=2.d0*(1.d0+nu_c)*mu*dudx_tangential_local(1,1)+nu_c*sigma_local(2,2)
                                  case default
                                    call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid 2D subtype')
                                end select
                                ! Stress tensor in global coordinates
                                sigma(1:2,1:2)=matmul(matmul(local_axes,sigma_local(1:2,1:2)),transpose(local_axes))
                                if (problem%subtype.eq.fbem_mechanics_plane_strain) then
                                  sigma(3,3)=nu_c*(sigma(1,1)+sigma(2,2))
                                end if
                              !
                              ! 3D
                              !
                              case (3)
                                ! Local stress tensor due to traction at boundary
                                sigma_local(3,1)=dot_product(node(sn)%value_c(4:6,face),local_axes(:,1))
                                sigma_local(3,2)=dot_product(node(sn)%value_c(4:6,face),local_axes(:,2))
                                sigma_local(3,3)=dot_product(node(sn)%value_c(4:6,face),local_axes(:,3))
                                sigma_local(1,3)=sigma_local(3,1)
                                sigma_local(2,3)=sigma_local(3,2)
                                ! Tangential displacement gradient in local cartesian coordinates
                                dudx_tangential_local=matmul(matmul(transpose(local_axes),dudx_tangential),local_axes)
                                ! Local stress tensor due to traction at boundary and tangential displacement
                                sigma_local(1,1)=1.d0/(1.d0-nu_c)*(nu_c*sigma_local(3,3)+2.d0*mu*(dudx_tangential_local(1,1)+nu_c*dudx_tangential_local(2,2)))
                                sigma_local(2,2)=1.d0/(1.d0-nu_c)*(nu_c*sigma_local(3,3)+2.d0*mu*(dudx_tangential_local(2,2)+nu_c*dudx_tangential_local(1,1)))
                                sigma_local(1,2)=mu*(dudx_tangential_local(1,2)+dudx_tangential_local(2,1))
                                sigma_local(2,1)=sigma_local(1,2)
                                ! Stress tensor in global coordinates
                                sigma=matmul(matmul(local_axes,sigma_local),transpose(local_axes))
                              case default
                                call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid dimension')
                            end select
                            ! SAVE STRESS TENSOR
                            do ki=1,3
                              do kj=1,3
                                kk=(ki-1)*3+kj
                                element(se)%value_c(kk,kn,face)=sigma(ki,kj)
                              end do
                            end do
                          end do
                        end do

                      !
                      ! POROELASTIC MEDIA
                      !
                      case (fbem_poroelastic)

                        do face=1,2
                          do kn=1,element(se)%n_nodes
                            sn=element(se)%node(kn)
                            xi=element(se)%xi_fn(:,kn)
                            psi=fbem_psi(problem%n,element(se)%type_g,element(se)%x_gn,element(se)%type_f1,element(se)%delta_f,xi)
                            !
                            ! CALCULATE FLUID DISPLACEMENTS
                            !
                            ! GRADIENT OF FLUID EQUIVALENT PRESSURE (TAU)
                            dpdx_tangential=0
                            do kn2=1,element(se)%n_nodes
                              sn2=element(se)%node(kn2)
                              dpdx_tangential = dpdx_tangential + psi(kn2,:)*node(sn2)%value_c(0,face)
                            end do
                            ! ADD
                            U_normal=node(sn)%value_c(problem%n+1,face)*element(se)%n_gn(:,kn)
                            if (face.eq.1) then
                              if (region(kr)%boundary_reversion(kb)) U_normal=-U_normal
                            else
                              if (.not.region(kr)%boundary_reversion(kb)) U_normal=-U_normal
                            end if
                            ! SAVE
                            element(se)%value_c(1:problem%n,kn,face)=U_normal-dpdx_tangential*pJ-pZ*node(sn2)%value_c(1:problem%n,face)
                            !
                            ! CALCULATE SOLID STRESS TENSOR
                            !
                            ! GRADIENT OF SOLID DISPLACEMENTS (TANGENTIAL PART)
                            dudx_tangential=0
                            do kn2=1,element(se)%n_nodes
                              sn2=element(se)%node(kn2)
                              do kc=1,problem%n
                                dudx_tangential(kc,:)=dudx_tangential(kc,:)+psi(kn2,:)*node(sn2)%value_c(kc,face)
                              end do
                            end do
                            ! LOCAL CARTESIAN AXES
                            if (face.eq.1) then
                              local_axes=fbem_local_tangential_cartesian_axes(problem%n,element(se)%type_g,element(se)%x_gn,region(kr)%boundary_reversion(kb),xi)
                            else
                              local_axes=fbem_local_tangential_cartesian_axes(problem%n,element(se)%type_g,element(se)%x_gn,.not.region(kr)%boundary_reversion(kb),xi)
                            end if
                            ! CALCULATE SOLID STRESS TENSOR
                            select case (problem%n)
                              !
                              ! 2D (plane strain)
                              !
                              case (2)
                                sigma_local=0
                                sigma=0
                                ! Local stress tensor due to traction at boundary
                                sigma_local(2,1)=dot_product(node(sn)%value_c(4:5,face),local_axes(:,1))
                                sigma_local(2,2)=dot_product(node(sn)%value_c(4:5,face),local_axes(:,2))
                                sigma_local(1,2)=sigma_local(2,1)
                                ! Tangential displacement gradient in local cartesian coordinates
                                dudx_tangential_local=matmul(matmul(transpose(local_axes),dudx_tangential),local_axes)
                                ! Local stress tensor due to traction at boundary and tangential displacement
                                select case (problem%subtype)
                                  case (fbem_mechanics_plane_strain)
                                    ! Ojo: falta chequear si la nu usada es la correcta para def. plana.
                                    sigma_local(1,1)=1.d0/(1.d0-nu_c)*(2.d0*mu*dudx_tangential_local(1,1)+nu_c*sigma_local(2,2))
                                  case default
                                    call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid 2D subtype (plane stress), not available for poroelastodynamics.')
                                end select
                                ! Stress tensor in global coordinates
                                sigma(1:2,1:2)=matmul(matmul(local_axes,sigma_local(1:2,1:2)),transpose(local_axes))
                                ! Ojo: falta chequear esto para def. plana.
                                !if (problem%subtype.eq.fbem_mechanics_plane_strain) then
                                !  sigma(3,3)=nu_c*(sigma(1,1)+sigma(2,2))
                                !end if
                                sigma(1,1)=sigma(1,1)+Q/R*node(sn)%value_c(0,face)
                                sigma(2,2)=sigma(2,2)+Q/R*node(sn)%value_c(0,face)
                              !
                              ! 3D
                              !
                              case (3)
                                ! Local stress tensor due to traction at boundary
                                sigma_local(3,1)=dot_product(node(sn)%value_c(5:7,face),local_axes(:,1))
                                sigma_local(3,2)=dot_product(node(sn)%value_c(5:7,face),local_axes(:,2))
                                sigma_local(3,3)=dot_product(node(sn)%value_c(5:7,face),local_axes(:,3))
                                sigma_local(1,3)=sigma_local(3,1)
                                sigma_local(2,3)=sigma_local(3,2)
                                ! Tangential displacement gradient in local cartesian coordinates
                                dudx_tangential_local=matmul(matmul(transpose(local_axes),dudx_tangential),local_axes)
                                ! Local stress tensor due to traction at boundary and tangential displacement
                                sigma_local(1,1)=1.d0/(1.d0-nu_c)*(nu_c*sigma_local(3,3)+2.d0*mu*(dudx_tangential_local(1,1)+nu_c*dudx_tangential_local(2,2)))
                                sigma_local(2,2)=1.d0/(1.d0-nu_c)*(nu_c*sigma_local(3,3)+2.d0*mu*(dudx_tangential_local(2,2)+nu_c*dudx_tangential_local(1,1)))
                                sigma_local(1,2)=mu*(dudx_tangential_local(1,2)+dudx_tangential_local(2,1))
                                sigma_local(2,1)=sigma_local(1,2)
                                ! Stress tensor in global coordinates
                                sigma=matmul(matmul(local_axes,sigma_local),transpose(local_axes))
                                sigma(1,1)=sigma(1,1)+Q/R*node(sn)%value_c(0,face)
                                sigma(2,2)=sigma(2,2)+Q/R*node(sn)%value_c(0,face)
                                sigma(3,3)=sigma(3,3)+Q/R*node(sn)%value_c(0,face)
                              case default
                                call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid dimension')
                            end select
                            ! SAVE
                            do ki=1,3
                              do kj=1,3
                                kk=problem%n+(ki-1)*3+kj
                                element(se)%value_c(kk,kn,face)=sigma(ki,kj)
                              end do
                            end do
                          end do
                        end do

                    end select

                end select

              !
              ! BE-BE coupled boundary of BE-FE-BE coupled boundary
              !
              case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)

                ! The region is the region 1 of the boundary
                if (.not.region(kr)%boundary_reversion(kb)) then
                  face=1
                ! The region is the region 2 of the boundary
                else
                  face=2
                end if

                select case (region(kr)%type)
                  !
                  ! INVISCID FLUID
                  !
                  case (fbem_potential)

                    do kn=1,element(se)%n_nodes
                      sn=element(se)%node(kn)
                      xi=element(se)%xi_fn(:,kn)
                      psi=fbem_psi(problem%n,element(se)%type_g,element(se)%x_gn,element(se)%type_f1,element(se)%delta_f,xi)
                      dpdx_tangential=0
                      do kn2=1,element(se)%n_nodes
                        sn2=element(se)%node(kn2)
                        dpdx_tangential = dpdx_tangential + psi(kn2,:)*node(sn2)%value_c(1,face)
                      end do
                      U_normal=node(sn)%value_c(2,face)*element(se)%n_gn(:,kn)
                      if (region(kr)%boundary_reversion(kb)) U_normal=-U_normal
                      element(se)%value_c(:,kn,face)=U_normal+dpdx_tangential*par_J
                    end do

                  !
                  ! VISCOELASTIC SOLID
                  !
                  case (fbem_viscoelastic)

                    do kn=1,element(se)%n_nodes
                      sn=element(se)%node(kn)
                      xi=element(se)%xi_fn(:,kn)
                      psi=fbem_psi(problem%n,element(se)%type_g,element(se)%x_gn,element(se)%type_f1,element(se)%delta_f,xi)
                      ! GRADIENT OF DISPLACEMENTS (TANGENTIAL PART)
                      dudx_tangential=0
                      do kn2=1,element(se)%n_nodes
                        sn2=element(se)%node(kn2)
                        do kc=1,problem%n
                          dudx_tangential(kc,:)=dudx_tangential(kc,:)+psi(kn2,:)*node(sn2)%value_c(kc,face)
                        end do
                      end do
                      ! LOCAL CARTESIAN AXES
                      local_axes=fbem_local_tangential_cartesian_axes(problem%n,element(se)%type_g,element(se)%x_gn,region(kr)%boundary_reversion(kb),xi)
                      ! CALCULATE STRESS TENSOR
                      select case (problem%n)
                        !
                        ! 2D (plane strain or plane stress)
                        !
                        case (2)
                          sigma_local=0
                          sigma=0
                          ! Local stress tensor due to traction at boundary
                          sigma_local(2,1)=dot_product(node(sn)%value_c(3:4,face),local_axes(:,1))
                          sigma_local(2,2)=dot_product(node(sn)%value_c(3:4,face),local_axes(:,2))
                          sigma_local(1,2)=sigma_local(2,1)
                          ! Tangential displacement gradient in local cartesian coordinates
                          dudx_tangential_local=matmul(matmul(transpose(local_axes),dudx_tangential),local_axes)
                          ! Local stress tensor due to traction at boundary and tangential displacement
                          select case (problem%subtype)
                            case (fbem_mechanics_plane_strain)
                              sigma_local(1,1)=1.d0/(1.d0-nu_c)*(2.d0*mu*dudx_tangential_local(1,1)+nu_c*sigma_local(2,2))
                            case (fbem_mechanics_plane_stress)
                              sigma_local(1,1)=2.d0*(1.d0+nu_c)*mu*dudx_tangential_local(1,1)+nu_c*sigma_local(2,2)
                            case default
                              call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid 2D subtype')
                          end select
                          ! Stress tensor in global coordinates
                          sigma(1:2,1:2)=matmul(matmul(local_axes,sigma_local(1:2,1:2)),transpose(local_axes))
                          if (problem%subtype.eq.fbem_mechanics_plane_strain) then
                            sigma(3,3)=nu_c*(sigma(1,1)+sigma(2,2))
                          end if
                        !
                        ! 3D
                        !
                        case (3)
                          ! Local stress tensor due to traction at boundary
                          sigma_local(3,1)=dot_product(node(sn)%value_c(4:6,face),local_axes(:,1))
                          sigma_local(3,2)=dot_product(node(sn)%value_c(4:6,face),local_axes(:,2))
                          sigma_local(3,3)=dot_product(node(sn)%value_c(4:6,face),local_axes(:,3))
                          sigma_local(1,3)=sigma_local(3,1)
                          sigma_local(2,3)=sigma_local(3,2)
                          ! Tangential displacement gradient in local cartesian coordinates
                          dudx_tangential_local=matmul(matmul(transpose(local_axes),dudx_tangential),local_axes)
                          ! Local stress tensor due to traction at boundary and tangential displacement
                          sigma_local(1,1)=1.d0/(1.d0-nu_c)*(nu_c*sigma_local(3,3)+2.d0*mu*(dudx_tangential_local(1,1)+nu_c*dudx_tangential_local(2,2)))
                          sigma_local(2,2)=1.d0/(1.d0-nu_c)*(nu_c*sigma_local(3,3)+2.d0*mu*(dudx_tangential_local(2,2)+nu_c*dudx_tangential_local(1,1)))
                          sigma_local(1,2)=mu*(dudx_tangential_local(1,2)+dudx_tangential_local(2,1))
                          sigma_local(2,1)=sigma_local(1,2)
                          ! Stress tensor in global coordinates
                          sigma=matmul(matmul(local_axes,sigma_local),transpose(local_axes))
                        case default
                          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid dimension')
                      end select
                      ! SAVE STRESS TENSOR
                      do ki=1,3
                        do kj=1,3
                          kk=(ki-1)*3+kj
                          element(se)%value_c(kk,kn,face)=sigma(ki,kj)
                        end do
                      end do
                    end do

                  ! POROELASTIC MEDIA
                  case (fbem_poroelastic)

                    do kn=1,element(se)%n_nodes
                      sn=element(se)%node(kn)
                      xi=element(se)%xi_fn(:,kn)
                      psi=fbem_psi(problem%n,element(se)%type_g,element(se)%x_gn,element(se)%type_f1,element(se)%delta_f,xi)
                      !
                      ! CALCULATE FLUID DISPLACEMENTS
                      !
                      ! GRADIENT OF FLUID EQUIVALENT PRESSURE (TAU)
                      dpdx_tangential=0
                      do kn2=1,element(se)%n_nodes
                        sn2=element(se)%node(kn2)
                        dpdx_tangential = dpdx_tangential + psi(kn2,:)*node(sn2)%value_c(0,face)
                      end do
                      ! ADD
                      U_normal=node(sn)%value_c(problem%n+1,face)*element(se)%n_gn(:,kn)
                      if (region(kr)%boundary_reversion(kb)) U_normal=-U_normal
                      ! SAVE
                      element(se)%value_c(1:problem%n,kn,face)=U_normal-dpdx_tangential*pJ-pZ*node(sn2)%value_c(1:problem%n,face)
                      !
                      ! CALCULATE SOLID STRESS TENSOR
                      !
                      ! GRADIENT OF SOLID DISPLACEMENTS (TANGENTIAL PART)
                      dudx_tangential=0
                      do kn2=1,element(se)%n_nodes
                        sn2=element(se)%node(kn2)
                        do kc=1,problem%n
                          dudx_tangential(kc,:)=dudx_tangential(kc,:)+psi(kn2,:)*node(sn2)%value_c(kc,face)
                        end do
                      end do
                      ! LOCAL CARTESIAN AXES
                      local_axes=fbem_local_tangential_cartesian_axes(problem%n,element(se)%type_g,element(se)%x_gn,region(kr)%boundary_reversion(kb),xi)
                      ! CALCULATE SOLID STRESS TENSOR
                      select case (problem%n)
                        !
                        ! 2D (plane strain)
                        !
                        case (2)
                          sigma_local=0
                          sigma=0
                          ! Local stress tensor due to traction at boundary
                          sigma_local(2,1)=dot_product(node(sn)%value_c(4:5,face),local_axes(:,1))
                          sigma_local(2,2)=dot_product(node(sn)%value_c(4:5,face),local_axes(:,2))
                          sigma_local(1,2)=sigma_local(2,1)
                          ! Tangential displacement gradient in local cartesian coordinates
                          dudx_tangential_local=matmul(matmul(transpose(local_axes),dudx_tangential),local_axes)
                          ! Local stress tensor due to traction at boundary and tangential displacement
                          select case (problem%subtype)
                            case (fbem_mechanics_plane_strain)
                              ! Ojo: falta chequear si la nu usada es la correcta para def. plana.
                              sigma_local(1,1)=1.d0/(1.d0-nu_c)*(2.d0*mu*dudx_tangential_local(1,1)+nu_c*sigma_local(2,2))
                            case default
                              call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid 2D subtype (plane stress), not available for poroelastodynamics.')
                          end select
                          ! Stress tensor in global coordinates
                          sigma(1:2,1:2)=matmul(matmul(local_axes,sigma_local(1:2,1:2)),transpose(local_axes))
                          ! Ojo: falta chequear esto para def. plana.
                          !if (problem%subtype.eq.fbem_mechanics_plane_strain) then
                          !  sigma(3,3)=nu_c*(sigma(1,1)+sigma(2,2))
                          !end if
                          sigma(1,1)=sigma(1,1)+Q/R*node(sn)%value_c(0,face)
                          sigma(2,2)=sigma(2,2)+Q/R*node(sn)%value_c(0,face)
                        !
                        ! 3D
                        !
                        case (3)
                          ! Local stress tensor due to traction at boundary
                          sigma_local(3,1)=dot_product(node(sn)%value_c(5:7,face),local_axes(:,1))
                          sigma_local(3,2)=dot_product(node(sn)%value_c(5:7,face),local_axes(:,2))
                          sigma_local(3,3)=dot_product(node(sn)%value_c(5:7,face),local_axes(:,3))
                          sigma_local(1,3)=sigma_local(3,1)
                          sigma_local(2,3)=sigma_local(3,2)
                          ! Tangential displacement gradient in local cartesian coordinates
                          dudx_tangential_local=matmul(matmul(transpose(local_axes),dudx_tangential),local_axes)
                          ! Local stress tensor due to traction at boundary and tangential displacement
                          sigma_local(1,1)=1.d0/(1.d0-nu_c)*(nu_c*sigma_local(3,3)+2.d0*mu*(dudx_tangential_local(1,1)+nu_c*dudx_tangential_local(2,2)))
                          sigma_local(2,2)=1.d0/(1.d0-nu_c)*(nu_c*sigma_local(3,3)+2.d0*mu*(dudx_tangential_local(2,2)+nu_c*dudx_tangential_local(1,1)))
                          sigma_local(1,2)=mu*(dudx_tangential_local(1,2)+dudx_tangential_local(2,1))
                          sigma_local(2,1)=sigma_local(1,2)
                          ! Stress tensor in global coordinates
                          sigma=matmul(matmul(local_axes,sigma_local),transpose(local_axes))
                          sigma(1,1)=sigma(1,1)+Q/R*node(sn)%value_c(0,face)
                          sigma(2,2)=sigma(2,2)+Q/R*node(sn)%value_c(0,face)
                          sigma(3,3)=sigma(3,3)+Q/R*node(sn)%value_c(0,face)
                        case default
                          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid dimension')
                      end select
                      ! SAVE
                      do ki=1,3
                        do kj=1,3
                          kk=problem%n+(ki-1)*3+kj
                          element(se)%value_c(kk,kn,face)=sigma(ki,kj)
                        end do
                      end do
                    end do

                end select

            end select

            deallocate(xi,psi)

          end do

        end do

    end if

  end do ! Loop through REGIONS

  ! ================================================================================================================================
  ! CALCULATE FEM ELEMENT STRESSES / STRESS RESULTANTS / EQUILIBRATING LOADS AND NODAL REACTIONS
  ! ================================================================================================================================

  ! Loop through the REGIONS
  do kr=1,n_regions
    !
    ! Finite elements
    !
    if (region(kr)%class.eq.fbem_fe) then

      !
      ! Igual habria que hacer algo....
      !
      if (region(kr)%type.eq.fbem_rigid) cycle

      !
      ! TODO: Usar el paradigma de region()%material(1)
      ! Use material of the FE REGION
      !
      ! Use material of the FE REGION
      if (region(kr)%subtype.eq.fbem_elastic_orthotropic) then
        ! ORTHOTROPIC
        stop 'orthotropic material only 3D'
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
      else
        ! ISOTROPIC
        E  =region(kr)%property_c(5)
        nu =region(kr)%property_r(3)
        rho=region(kr)%property_r(1)
      end if

      ! Loop through the FE SUBREGIONS of the FE REGION
      do ks_int=1,region(kr)%n_fe_subregions
        ss_int=region(kr)%fe_subregion(ks_int)

        ! Loop through the FINITE ELEMENTS of the FE SUBREGION
        do ke_int=1,part(fe_subregion(ss_int)%part)%n_elements
          se_int=part(fe_subregion(ss_int)%part)%element(ke_int)

          ! ==================================================================
          ! STRESSES / STRESS RESULTANTS: element(*)%value_c(component,node,1)
          ! ==================================================================

          select case (element(se_int)%n_dimension)

            ! ======================================================================================================================
            ! ONE-DIMENSIONAL ELEMENTS
            ! ======================================================================================================================

            case (0)

              !---------------------------------------------------------------------------------------------------------------------
              ! POINT MASS ELEMENT
              !
              ! Nothing to do
              !---------------------------------------------------------------------------------------------------------------------

            ! ======================================================================================================================
            ! ONE-DIMENSIONAL ELEMENTS
            ! ======================================================================================================================

            case (1)

              select case (element(se_int)%fe_type)

                case (0)

                  !-----------------------------------------------------------------------------------------------------------------
                  ! DEGENERATED BEAM FINITE ELEMENT
                  !
                  ! Build the element DOF vector
                  allocate (a(element(se_int)%n_dof))
                  kdof=1
                  do kn=1,element(se_int)%n_nodes
                    sn=element(se_int)%node(kn)
                    a((kdof):(kdof-1+element(se_int)%node_n_dof(kn)))=node(sn)%value_c(1:element(se_int)%node_n_dof(kn),1)
                    kdof=kdof+element(se_int)%node_n_dof(kn)
                  end do
                  ! Calculate the STRESS RESULTANTS INTERPOLATION SCHEME
                  select case (problem%n)
                    case (2)
                      select case (problem%subtype)
                        case (fbem_mechanics_plane_strain)
                          call fbem_fem_degbeam_stress_resultants(2,1,&
                                                                  element(se_int)%type,&
                                                                  element(se_int)%x_gn,&
                                                                  element(se_int)%v_midnode,&
                                                                  element(se_int)%tv_midnode,&
                                                                  nu,element(se_int)%ksh,setype,sedelta,Fsigma)
                        case (fbem_mechanics_plane_stress)
                          call fbem_fem_degbeam_stress_resultants(2,2,&
                                                                  element(se_int)%type,&
                                                                  element(se_int)%x_gn,&
                                                                  element(se_int)%v_midnode,&
                                                                  element(se_int)%tv_midnode,&
                                                                  nu,element(se_int)%ksh,setype,sedelta,Fsigma)
                        case default
                          call fbem_error_message(error_unit,0,'element',element(se_int)%id,'invalid type of 2D analysis')
                      end select
                    case (3)
                      call fbem_fem_degbeam_stress_resultants(3,0,&
                                                              element(se_int)%type,&
                                                              element(se_int)%x_gn,&
                                                              element(se_int)%v_midnode,&
                                                              element(se_int)%tv_midnode,&
                                                              nu,element(se_int)%ksh,setype,sedelta,Fsigma)
                    case default
                      call fbem_error_message(error_unit,0,'element',element(se_int)%id,'invalid problem%n')
                  end select


                  allocate (F(3*(problem%n-1),fbem_n_nodes(setype)))
                  ! Calculate at Gauss points
                  do kn=1,fbem_n_nodes(setype)
                    F(:,kn)=E*matmul(FSigma(:,:,kn),a)
                  end do
                  ! Extrapolate to nodes
                  se_int_type=element(se_int)%type
                  allocate (phi(fbem_n_nodes(setype)))
                  element(se_int)%value_c=0.d0
                  do kn=1,element(se_int)%n_nodes
                    ! Local coordinates at node kn
#                   define xi xi1d
#                   define etype se_int_type
#                   define node kn
#                   define delta 0.d0
#                   include <xi_1d_at_node.rc>
#                   undef etype
#                   undef node
#                   undef delta
                    ! Shape functions of stress resultants at node kn
#                   define etype setype
#                   define delta sedelta
#                   include <phi_1d.rc>
#                   undef etype
#                   undef delta
#                   undef xi
                    ! Calculate
                    do kn2=1,fbem_n_nodes(setype)
                      element(se_int)%value_c(1:3*(problem%n-1),kn,1)=element(se_int)%value_c(1:3*(problem%n-1),kn,1)+phi(kn2)*F(:,kn2)
                    end do
                  end do
                  deallocate (a,F,Fsigma,phi)
                  !
                  !-----------------------------------------------------------------------------------------------------------------

                case (1,2)

                  !-----------------------------------------------------------------------------------------------------------------
                  ! STRAIGHT EULER-BERNOULLI AND TIMOSHENKO BEAM FINITE ELEMENTS
                  !
                  ! Build the element DOF vector
                  allocate (a(3*(problem%n-1)*element(se_int)%n_nodes))
                  kdof=1
                  do kn=1,element(se_int)%n_nodes
                    sn=element(se_int)%node(kn)
                    a((kdof):(kdof-1+element(se_int)%node_n_dof(kn)))=node(sn)%value_c(1:element(se_int)%node_n_dof(kn),1)
                    kdof=kdof+element(se_int)%node_n_dof(kn)
                  end do
                  !
                  ! TO-DO: Los ejes nodales habría que guardarlos para cada nodo y construirlos para cada elemento
                  !
                  allocate(nodal_axes(problem%n,problem%n,element(se_int)%n_nodes))
                  nodal_axes=0
                  do kn=1,element(se_int)%n_nodes
                    do kc=1,problem%n
                      nodal_axes(kc,kc,kn)=1
                    end do
                  end do
                  !
                  ! Calculate the STRESS RESULTANTS INTERPOLATION SCHEME
                  !
                  ! NOTA: viga plane strain (A debe introducirse ya con espesor unitario, y para plane stress con el espesor que toque)
                  call fbem_fem_strbeam_stress_resultants(problem%n,&
                                                          element(se_int)%type,element(se_int)%fe_options(1),element(se_int)%fe_type,&
                                                          element(se_int)%x_gn,element(se_int)%ep,&
                                                          element(se_int)%A,element(se_int)%I(1:3),element(se_int)%ksh,&
                                                          nu,nodal_axes,setype,sedelta,Fsigma)
                  allocate (F(3*(problem%n-1),fbem_n_nodes(setype)))
                  ! Calculate at Gauss points
                  do kn=1,fbem_n_nodes(setype)
                    F(:,kn)=E*matmul(FSigma(:,:,kn),a)
                  end do
                  ! Extrapolate to nodes
                  se_int_type=element(se_int)%type
                  allocate (phi(fbem_n_nodes(setype)))
                  element(se_int)%value_c=0.d0
                  do kn=1,element(se_int)%n_nodes
                    ! Local coordinates at node kn
#                   define xi xi1d
#                   define etype se_int_type
#                   define node kn
#                   define delta 0.d0
#                   include <xi_1d_at_node.rc>
#                   undef etype
#                   undef node
#                   undef delta
                    ! Shape functions of stress resultants at node kn
#                   define etype setype
#                   define delta sedelta
#                   include <phi_1d.rc>
#                   undef etype
#                   undef delta
#                   undef xi
                    ! Calculate
                    do kn2=1,fbem_n_nodes(setype)
                      element(se_int)%value_c(1:3*(problem%n-1),kn,1)=element(se_int)%value_c(1:3*(problem%n-1),kn,1)+phi(kn2)*F(:,kn2)
                    end do
                  end do
                  deallocate (a,nodal_axes,F,Fsigma,phi)
                  !
                  !-----------------------------------------------------------------------------------------------------------------

                case (3)

                  !-----------------------------------------------------------------------------------------------------------------
                  ! BAR FINITE ELEMENTS
                  !
                  ! Build the element DOF vector
                  allocate (a(problem%n*2))
                  kdof=1
                  do kn=1,2
                    sn=element(se_int)%node(kn)
                    a((kdof):(kdof-1+element(se_int)%node_n_dof(kn)))=node(sn)%value_c(1:element(se_int)%node_n_dof(kn),1)
                    kdof=kdof+element(se_int)%node_n_dof(kn)
                  end do
                  !
                  ! TO-DO: Los ejes nodales habría que guardarlos para cada nodo y construirlos para cada elemento
                  !
                  allocate(nodal_axes(problem%n,problem%n,2))
                  nodal_axes=0
                  do kn=1,2
                    do kc=1,problem%n
                      nodal_axes(kc,kc,kn)=1
                    end do
                  end do
                  !
                  ! Calculate the SPRINGS FORCES
                  !
                  allocate (f_e(2*problem%n),KLT(2,2*problem%n))
                  call fbem_fem_bar_KLT_harmonic(problem%n,omega,element(se_int)%x_gn,element(se_int)%A,nodal_axes,E,KLT)
                  f_e=matmul(KLT,a)
                  element(se_int)%value_c(1,1,1)=-f_e(1)
                  element(se_int)%value_c(1,2,1)= f_e(2)
                  deallocate (a,nodal_axes,f_e,KLT)
                  !
                  !-----------------------------------------------------------------------------------------------------------------

                case (4)

                  !-----------------------------------------------------------------------------------------------------------------
                  ! DISCRETE TRANSLATIONAL SPRING FINITE ELEMENTS
                  !
                  ! Build the element DOF vector
                  allocate (a(2*problem%n))
                  kdof=1
                  do kn=1,2
                    sn=element(se_int)%node(kn)
                    a((kdof):(kdof-1+element(se_int)%node_n_dof(kn)))=node(sn)%value_c(1:element(se_int)%node_n_dof(kn),1)
                    kdof=kdof+element(se_int)%node_n_dof(kn)
                  end do
                  !
                  ! TO-DO: Los ejes nodales habría que guardarlos para cada nodo y construirlos para cada elemento
                  !
                  allocate(nodal_axes(problem%n,problem%n,2))
                  nodal_axes=0
                  do kn=1,2
                    do kc=1,problem%n
                      nodal_axes(kc,kc,kn)=1
                    end do
                  end do
                  !
                  ! Calculate the SPRINGS FORCES
                  !
                  allocate (f_e(2*problem%n),KLT(2*problem%n,2*problem%n))
                  call fbem_fem_distra_KLT_harmonic(problem%n,element(se_int)%ep,nodal_axes,&
                  element(se_int)%k_c+c_im*omega*element(se_int)%c,KLT)
                  f_e=matmul(KLT,a)
                  element(se_int)%value_c(1:problem%n,1,1)=-f_e(1:problem%n)
                  element(se_int)%value_c(1:problem%n,2,1)= f_e((problem%n+1):(2*problem%n))
                  deallocate (a,nodal_axes,f_e,KLT)
                  !
                  !-----------------------------------------------------------------------------------------------------------------

                case (5)

                  !-----------------------------------------------------------------------------------------------------------------
                  ! DISCRETE ROTATIONAL/TRANSLATIONAL SPRING FINITE ELEMENTS
                  !
                  ! Build the element DOF vector
                  allocate (a(2*3*(problem%n-1)))
                  kdof=1
                  do kn=1,2
                    sn=element(se_int)%node(kn)
                    a((kdof):(kdof-1+element(se_int)%node_n_dof(kn)))=node(sn)%value_c(1:element(se_int)%node_n_dof(kn),1)
                    kdof=kdof+element(se_int)%node_n_dof(kn)
                  end do
                  !
                  ! TO-DO: Los ejes nodales habría que guardarlos para cada nodo y construirlos para cada elemento
                  !
                  allocate(nodal_axes(problem%n,problem%n,2))
                  nodal_axes=0
                  do kn=1,2
                    do kc=1,problem%n
                      nodal_axes(kc,kc,kn)=1
                    end do
                  end do
                  !
                  ! Calculate the SPRINGS FORCES
                  !
                  allocate (f_e(2*3*(problem%n-1)),KLT(2*3*(problem%n-1),2*3*(problem%n-1)))
                  call fbem_fem_disrotra_KLT_harmonic(problem%n,element(se_int)%ep,nodal_axes,&
                  element(se_int)%k_c+c_im*omega*element(se_int)%c,KLT)
                  f_e=matmul(KLT,a)
                  element(se_int)%value_c(1:3*(problem%n-1),1,1)=-f_e(1:3*(problem%n-1))
                  element(se_int)%value_c(1:3*(problem%n-1),2,1)= f_e((3*(problem%n-1)+1):(2*3*(problem%n-1)))
                  deallocate (a,nodal_axes,f_e,KLT)
                  !
                  !-----------------------------------------------------------------------------------------------------------------

                case (6)

                  !-----------------------------------------------------------------------------------------------------------------
                  ! DISCRETE SPRING-DASHPOT
                  !
                  ! Build the element DOF vector
                  allocate (a(problem%n*2))
                  kdof=1
                  do kn=1,2
                    sn=element(se_int)%node(kn)
                    a((kdof):(kdof-1+element(se_int)%node_n_dof(kn)))=node(sn)%value_c(1:element(se_int)%node_n_dof(kn),1)
                    kdof=kdof+element(se_int)%node_n_dof(kn)
                  end do
                  !
                  ! TO-DO: Los ejes nodales habría que guardarlos para cada nodo y construirlos para cada elemento
                  !
                  allocate(nodal_axes(problem%n,problem%n,2))
                  nodal_axes=0
                  do kn=1,2
                    do kc=1,problem%n
                      nodal_axes(kc,kc,kn)=1
                    end do
                  end do
                  !
                  ! Calculate the SPRINGS FORCES
                  !
                  allocate (f_e(2*problem%n),KLT(2,2*problem%n))

                  call fbem_fem_springdashpot_KLT_harmonic(problem%n,omega,element(se_int)%x_gn,&
                                                           element(se_int)%k_c(1),element(se_int)%c(1),nodal_axes,KLT)
                  f_e=matmul(KLT,a)
                  element(se_int)%value_c(1,1,1)=-f_e(1)
                  element(se_int)%value_c(1,2,1)= f_e(2)
                  deallocate (a,nodal_axes,f_e,KLT)
                  !
                  !-----------------------------------------------------------------------------------------------------------------


                case default

                  !-----------------------------------------------------------------------------------------------------------------
                  ! OTHER TYPES
                  !
                  call fbem_error_message(error_unit,0,'element',element(se_int)%id,'invalid type of 1D element')
                  !
                  !-----------------------------------------------------------------------------------------------------------------

              end select

            ! ======================================================================================================================
            ! TWO-DIMENSIONAL ELEMENTS
            ! ======================================================================================================================

            case (2)

              select case (problem%n)

                case (2)

                  !-----------------------------------------------------------------------------------------------------------------
                  ! SOLID / CONTINUUM ELEMENTS
                  !

                  !
                  ! TO BE DONE...
                  !

                  !
                  !
                  !-----------------------------------------------------------------------------------------------------------------

                case (3)

                  !-----------------------------------------------------------------------------------------------------------------
                  ! DEGENERATED SHELL FINITE ELEMENT
                  !
                  ! Shell geometrical data
                  allocate (x1ref(3))
                  x1ref=element(se_int)%ep(:,1)
                  ! Calculate COORDINATE TRANSFORMATION MATRIX (Lc)
                  allocate (Lc(element(se_int)%n_dof,5*element(se_int)%n_nodes))
                  call fbem_fem_degshell_Lc(element(se_int)%type,element(se_int)%n_dof,element(se_int)%node_n_dof,&
                                            element(se_int)%v_midnode,Lc)
                  ! Build the element DOF vector with displacements and local rotations
                  allocate (a(element(se_int)%n_dof),ap(5*element(se_int)%n_nodes))
                  kdof=1
                  do kn=1,element(se_int)%n_nodes
                    sn=element(se_int)%node(kn)
                    a((kdof):(kdof-1+element(se_int)%node_n_dof(kn)))=node(sn)%value_c(1:element(se_int)%node_n_dof(kn),1)
                    kdof=kdof+element(se_int)%node_n_dof(kn)
                  end do
                  ap=matmul(transpose(Lc),a)
                  ! Calculate the STRESS RESULTANTS INTERPOLATION SCHEME
                  call fbem_fem_degshell_stress_resultants(element(se_int)%type,&
                                                           element(se_int)%x_gn,element(se_int)%v_midnode,&
                                                           element(se_int)%tv_midnode,x1ref,nu,element(se_int)%ksh,&
                                                           setype,sedelta,Fsigma)
                  allocate (F(8,fbem_n_nodes(setype)))
                  ! Calculate at Gauss points
                  do kn=1,fbem_n_nodes(setype)
                    F(:,kn)=E*matmul(FSigma(:,:,kn),ap)
                  end do
                  ! Extrapolate to nodes
                  se_int_type=element(se_int)%type
                  allocate (phi(fbem_n_nodes(setype)))
                  element(se_int)%value_c=0.d0
                  do kn=1,element(se_int)%n_nodes
                    ! Local coordinates at node kn
#                   define xi xi2d
#                   define etype se_int_type
#                   define node kn
#                   define delta 0.d0
#                   include <xi_2d_at_node.rc>
#                   undef etype
#                   undef node
#                   undef delta
                    ! Shape functions of stress resultants at node kn
#                   define etype setype
#                   define delta sedelta
#                   include <phi_2d.rc>
#                   undef etype
#                   undef delta
#                   undef xi
                    ! Calculate
                    do kn2=1,fbem_n_nodes(setype)
                      element(se_int)%value_c(1:8,kn,1)=element(se_int)%value_c(1:8,kn,1)+phi(kn2)*F(:,kn2)
                    end do
                  end do
                  deallocate (a,ap,F,Lc,Fsigma,x1ref,phi)
                  !
                  !-----------------------------------------------------------------------------------------------------------------

              end select

            ! ======================================================================================================================
            ! THREE-DIMENSIONAL ELEMENTS
            ! ======================================================================================================================

            case (3)

              !---------------------------------------------------------------------------------------------------------------------
              ! SOLID / CONTINUUM ELEMENTS
              !
              call fbem_error_message(error_unit,0,'element',element(se_int)%id,'3D elements not available yet')
              !
              !---------------------------------------------------------------------------------------------------------------------

          end select

          ! =============================================================================
          ! EQUILIBRATING LOADS (q=K·a-f_bem-f_ext): element(*)%value_c(component,node,2)
          ! =============================================================================

          !
          ! K·a
          !
          ! Build the stiffness matrix (K)
          ndof_mat=3*(problem%n-1)*element(se_int)%n_nodes
          allocate (K(ndof_mat,ndof_mat))
          call build_fem_K_harmonic(kf,se_int,ndof_mat,K)
          ! Build the element DOF vector (a)
          allocate (a(ndof_mat))
          a=0
          kdofe=1
          do kn=1,element(se_int)%n_nodes
            sn=element(se_int)%node(kn)
            do kc=1,element(se_int)%node_n_dof(kn)
              a(kdofe)=node(sn)%value_c(kc,1)
              kdofe=kdofe+1
            end do
          end do
          ! Multiply K·a
          a=matmul(K,a)
          ! Save to q
          kdofe=1
          do kn=1,element(se_int)%n_nodes
            sn=element(se_int)%node(kn)
            do kc=1,element(se_int)%node_n_dof(kn)
              element(se_int)%value_c(kc,kn,2)=a(kdofe)
              kdofe=kdofe+1
            end do
          end do
          deallocate(a,K)
          !
          ! f_bem
          !
          !
          ! TODO: Falta validar
          !
          call assemble_fem_harela_coupling(se_int,1)
          !
          ! f_ext
          !
          !
          ! Build the load vector (f_ext)
          ndof_mat=3*(problem%n-1)*element(se_int)%n_nodes
          allocate (f_e(ndof_mat))
          call build_fem_f_distributed_harmonic(kf,se_int,ndof_mat,f_e)
          !
          ! TODO: creo que en este paso no influye regiones rigidas
          !
          ! Save to q
          kdofe=1
          do kn=1,element(se_int)%n_nodes
            sn=element(se_int)%node(kn)
            do kc=1,element(se_int)%node_n_dof(kn)
              element(se_int)%value_c(kc,kn,2)=element(se_int)%value_c(kc,kn,2)-f_e(kdofe)
              kdofe=kdofe+1
            end do
          end do
          deallocate(f_e)

          ! =============================================
          ! NODAL REACTIONS: node(*)%value_c(component,2)
          ! =============================================

          do kn=1,element(se_int)%n_nodes
            sn=element(se_int)%node(kn)
            do kc=1,element(se_int)%node_n_dof(kn)
              node(sn)%value_c(kc,2)=node(sn)%value_c(kc,2)+element(se_int)%value_c(kc,kn,2)
            end do
          end do

!~           !
!~           ! Write
!~           !
!~           do kn=1,element(se_int)%n_nodes
!~             sn=element(se_int)%node(kn)
!~             write(63,*) element(se_int)%id,node(sn)%id,element(se_int)%value_c(1:element(se_int)%node_n_dof(kn),kn,2)
!~           end do

        end do ! Loop through the FINITE ELEMENTS of the FE SUBREGION

      end do ! Loop through the FE SUBREGIONS of the FE REGION

    end if
  end do ! Loop through the REGIONS

  ! Ending message
  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'END calculating stresses at FE regions')

end subroutine calculate_stresses_mechanics_harmonic
