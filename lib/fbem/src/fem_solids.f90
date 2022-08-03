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

!! @author Jacob David Rodriguez Bordon (jacobdavid.rodriguezbordon@ulpgc.es)
!!
!! @version 2.0
!!
!! <b> Solid linear elastic finite elements </b>
module fbem_fem_solids

  ! Fortran 2003 standard
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_quad_rules
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_geometry

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all are private
  private

  !
  public  :: fbem_fem_solid2d_K_static
  public  :: fbem_fem_solid2d_K_harmonic
  !
  public :: fbem_fem_staela2d_solid_D_isotropic
  public :: fbem_fem_staela2d_solid_K
  !
  public :: fbem_fem_harela2d_solid_D_isotropic
  public :: fbem_fem_harela2d_solid_K
  !
  public :: fbem_fem_ela2d_solid_Q_body
  public :: fbem_fem_ela2d_solid_Q_edge
  public :: fbem_fem_ela2d_solid_M
  !
  public :: fbem_fem_staela2d_solid_dKdx
  public :: fbem_fem_staela2d_solid_dKdx_FD
  public :: fbem_fem_harela2d_solid_dKdx
  public :: fbem_fem_harela2d_solid_dKdx_FD
  !
  public :: fbem_fem_ela2d_solid_dMdx
  public :: fbem_fem_ela2d_solid_dMdx_FD
  public :: fbem_fem_ela2d_solid_dQdx_body
  public :: fbem_fem_ela2d_solid_dQdx_body_FD
  public :: fbem_fem_ela2d_solid_dQdx_edge
  public :: fbem_fem_ela2d_solid_dQdx_edge_FD

contains

  ! Calculate the stiffness matrix K for 2D static analysis of isotropic solids.
  subroutine fbem_fem_solid2d_K_static(etype,plane,x_nodes,t_nodes,E,nu,K_intmode,K_intngp,K)
    implicit none
    ! I/O
    integer           :: etype                                          !! Type of element: tri3, tri6, quad4, quad8, quad9.
    integer           :: plane                                          !! 1 (plane strain), 2 (plane_stress)
    real(kind=real64) :: x_nodes(2,fbem_n_nodes(etype))                 !! Coordinates of the nodes.
    real(kind=real64) :: t_nodes(fbem_n_nodes(etype))                   !! Thickness at each node (t=1 for 2D plane strain problem)
    real(kind=real64) :: E                                              !! Young' modulus
    real(kind=real64) :: nu                                             !! Poisson's ratio
    integer           :: K_intmode                                      !! Integration mode (K): 0 (full), (1) reduced, (2) selective, (3) user-defined
    integer           :: K_intngp(3)                                    !! Number of integration points for user-defined mode (K): (1) membrane, (2) shear
    real(kind=real64) :: K(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype)) !! Stiffness matrix
    ! Local
    real(kind=real64) :: D(3,3)                                         !! Strain-stress constitutive matrix in global axes.
    integer           :: gln_a                                          !! Number of Gauss-Legendre integration points for axial contribution
    integer           :: gln_sh                                         !! Number of Gauss-Legendre integration points for shear contribution.
    ! Calculate D
    select case (plane)
      case (1)
        call fbem_fem_staela2d_solid_D_isotropic(E,nu,1,D)
      case (2)
        call fbem_fem_staela2d_solid_D_isotropic(E,nu,2,D)
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of plane')
    end select
    ! Integration mode
    select case (K_intmode)
      ! Full integration
      case (0)
        select case (etype)
          case (fbem_tri3,fbem_quad4)
            gln_a=2
            gln_sh=2
          case (fbem_tri6,fbem_quad8,fbem_quad9)
            gln_a=3
            gln_sh=3
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of element for full integration')
        end select
      ! Reduced integration
      case (1)
        select case (etype)
          case (fbem_tri3,fbem_quad4)
            gln_a=1
            gln_sh=1
          case (fbem_tri6,fbem_quad8,fbem_quad9)
            gln_a=2
            gln_sh=2
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of element for reduced integration')
        end select
      ! Selective integration
      case (2)
        select case (etype)
          case (fbem_tri3,fbem_quad4)
            gln_a=2
            gln_sh=1
          case (fbem_tri6,fbem_quad8,fbem_quad9)
            gln_a=3
            gln_sh=2
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of element for selective integration')
        end select
      ! User-defined integration
      case (3)
        gln_a=K_intngp(1)
        gln_sh=K_intngp(2)
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid integration mode')
    end select
    call fbem_fem_staela2d_solid_K(etype,x_nodes,t_nodes,D,gln_a,gln_sh,K)
  end subroutine fbem_fem_solid2d_K_static

  ! Calculate the stiffness matrix K for 2D time harmonic analysis of isotropic solids.
  subroutine fbem_fem_solid2d_K_harmonic(omega,etype,plane,x_nodes,t_nodes,E,nu,rho,K_intmode,K_intngp,M_intmode,M_intngp,K)
    implicit none
    ! I/O
    real(kind=real64)    :: omega                                          !! Circular frequency
    integer              :: etype                                          !! Type of element: tri3, tri6, quad4, quad8, quad9.
    integer              :: plane                                          !! 1 (plane strain), 2 (plane_stress)
    real(kind=real64)    :: x_nodes(2,fbem_n_nodes(etype))                 !! Coordinates of the nodes.
    real(kind=real64)    :: t_nodes(fbem_n_nodes(etype))                   !! Thickness at each node (t=1 for 2D plane strain problem)
    complex(kind=real64) :: E                                              !! Young' modulus
    real(kind=real64)    :: nu                                             !! Poisson's ratio
    real(kind=real64)    :: rho                                            !! Poisson's ratio
    integer              :: K_intmode                                      !! Integration mode (K): 0 (full), (1) reduced, (2) selective, (3) user-defined
    integer              :: K_intngp(3)                                    !! Number of integration points for user-defined mode (K): (1) membrane, (2) shear
    integer              :: M_intmode                                      !! Integration mode (M): 0 (full), (1) reduced, (2) selective, (3) user-defined
    integer              :: M_intngp(3)                                    !! Number of integration points for user-defined mode (M): (1) used
    complex(kind=real64) :: K(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype)) !! Time harmonic stiffness matrix
    ! Local
    real(kind=real64)    :: M(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype)) !! Mass matrix
    complex(kind=real64) :: D(3,3)                                         !! Strain-stress constitutive matrix in global axes.
    integer              :: gln_a                                          !! Number of Gauss-Legendre integration points for axial contribution
    integer              :: gln_sh                                         !! Number of Gauss-Legendre integration points for shear contribution.
    ! Calculate D
    select case (plane)
      case (1)
        call fbem_fem_harela2d_solid_D_isotropic(E,nu,1,D)
      case (2)
        call fbem_fem_harela2d_solid_D_isotropic(E,nu,2,D)
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of plane')
    end select
    ! Integration mode (K)
    select case (K_intmode)
      ! Full integration
      case (0)
        select case (etype)
          case (fbem_tri3,fbem_quad4)
            gln_a=2
            gln_sh=2
          case (fbem_tri6,fbem_quad8,fbem_quad9)
            gln_a=3
            gln_sh=3
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of element for full integration')
        end select
      ! Reduced integration
      case (1)
        select case (etype)
          case (fbem_tri3,fbem_quad4)
            gln_a=1
            gln_sh=1
          case (fbem_tri6,fbem_quad8,fbem_quad9)
            gln_a=2
            gln_sh=2
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of element for reduced integration')
        end select
      ! Selective integration
      case (2)
        select case (etype)
          case (fbem_tri3,fbem_quad4)
            gln_a=2
            gln_sh=1
          case (fbem_tri6,fbem_quad8,fbem_quad9)
            gln_a=3
            gln_sh=2
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of element for selective integration')
        end select
      ! User-defined integration
      case (3)
        gln_a=K_intngp(1)
        gln_sh=K_intngp(1)
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid integration mode')
    end select
    call fbem_fem_harela2d_solid_K(etype,x_nodes,t_nodes,D,gln_a,gln_sh,K)
    ! Integration mode (M)
    select case (K_intmode)
      ! Full integration
      case (0)
        select case (etype)
          case (fbem_tri3,fbem_quad4)
            gln_a=2
          case (fbem_tri6,fbem_quad8,fbem_quad9)
            gln_a=3
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of element for full integration')
        end select
      ! Reduced integration
      case (1,2)
        select case (etype)
          case (fbem_tri3,fbem_quad4)
            gln_a=1
          case (fbem_tri6,fbem_quad8,fbem_quad9)
            gln_a=2
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid type of element for reduced integration')
        end select
      ! User-defined integration
      case (3)
        gln_a=M_intngp(1)
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid integration mode')
    end select
    call fbem_fem_ela2d_solid_M(etype,x_nodes,t_nodes,rho,gln_a,M)
    K=K-omega**2*M
  end subroutine fbem_fem_solid2d_K_harmonic

  !! Constitutive law matrix for 2D static analysis of isotropic solids.
  subroutine fbem_fem_staela2d_solid_D_isotropic(E,nu,plane,D)
    implicit none
    real(kind=real64) :: E      ! Young's modulus
    real(kind=real64) :: nu     ! Poisson's ratio
    integer           :: plane  ! 1 (plane strain), 2 (plane_stress)
    real(kind=real64) :: D(3,3) ! Constitutive matrix
    D=0.d0
    select case (plane)
      case (1)
        D(1,1)=E*(1.d0-nu)/(1.d0+nu)/(1.d0-2.d0*nu)
        D(2,2)=D(1,1)
        D(1,2)=nu/(1.d0-nu)*D(1,1)
        D(2,1)=D(1,2)
        D(3,3)=0.5d0*E/(1.d0+nu)
      case (2)
        D(1,1)=E/(1.d0-nu**2)
        D(2,2)=D(1,1)
        D(1,2)=nu*D(1,1)
        D(2,1)=D(1,2)
        D(3,3)=0.5d0*E/(1.d0+nu)
    end select
  end subroutine fbem_fem_staela2d_solid_D_isotropic

  !! Calculate the stiffness matrix K for 2D static analysis of isotropic solids.
  subroutine fbem_fem_staela2d_solid_K(etype,x_nodes,t_nodes,D,gln_a,gln_sh,K)
    implicit none
    ! I/O
    integer           :: etype                                          !! Type of element: tri3, tri6, quad4, quad8, quad9.
    real(kind=real64) :: x_nodes(2,fbem_n_nodes(etype))                 !! Coordinates of the nodes.
    real(kind=real64) :: t_nodes(fbem_n_nodes(etype))                   !! Thickness at each node (t=1. for 2D plane strain problem)
    real(kind=real64) :: D(3,3)                                         !! Strain-stress constitutive matrix in global axes.
    integer           :: gln_a                                          !! Number of Gauss-Legendre integration points for axial contribution
    integer           :: gln_sh                                         !! Number of Gauss-Legendre integration points for shear contribution.
    real(kind=real64) :: K(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype)) !! Stiffness matrix
    ! Local
    integer           :: n_nodes                                        ! Number of nodes
    integer           :: i                                              ! Counter
    integer           :: contribution                                   ! Contribution part
    integer           :: gln                                            ! Number of Gauss points
    integer           :: rule                                           ! Rule of Wandzura quadrature
    integer           :: kxi1, kxi2, kxit                               ! Integration points counters
    real(kind=real64) :: xi(2), w(2), wt                                ! Curvilinear coordinates and quadrature weights
    real(kind=real64) :: aux(10)                                        ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64) :: phi(fbem_n_nodes(etype))                       ! Shape functions
    real(kind=real64) :: dphidxi1(fbem_n_nodes(etype))                  ! Shape functions derivatives with respect to xi_1
    real(kind=real64) :: dphidxi2(fbem_n_nodes(etype))                  ! Shape functions derivatives with respect to xi_2
    real(kind=real64) :: J(2,2), invJ(2,2), detJ                        ! Jacobian matrix, its inverse and the jacobian determinant
    real(kind=real64) :: jw                                             ! det(J) * weights
    real(kind=real64) :: t                                              ! Thickness at the integration point
    real(kind=real64) :: Dcopy(3,3)                                     ! Strain-stress constitutive matrix in global axes (local copy).
    real(kind=real64) :: dNdx(fbem_n_nodes(etype),2)                    ! dN/dx
    real(kind=real64) :: B(3,2,fbem_n_nodes(etype))                     ! B matrix
    integer           :: ki, kis, kie, kj, kjs, kje                     ! Counters and nodal DOF limits
    ! Initialization
    K=0.d0
    n_nodes=fbem_n_nodes(etype)
    ! Axial/shear contributions
    do contribution=1,2
      ! Decide how the element is going to be integrated
      Dcopy=0.d0
      if (gln_a.ne.gln_sh) then
        select case (contribution)
          case (1)
            gln=gln_a
            Dcopy(1:2,1:2)=D(1:2,1:2)
          case (2)
            gln=gln_sh
            Dcopy(1:2,3)=D(1:2,3)
            Dcopy(  3,:)=D(  3,:)
        end select
      else
        gln=gln_a
        Dcopy=D
        if (contribution.eq.2) exit
      end if
      ! Switch between types of elements
      select case (fbem_n_edges(etype))
        ! TRIANGULAR ELEMENTS
        case (3)
          ! Loop through integration points
          ! Transform to order for Wandzura rules
          rule=2*gln-1
          do kxit=1,wantri_n(rule)
            ! xi_1, xi_2 and w_t
            xi(1)=wantri_xi1(kxit,rule)
            xi(2)=wantri_xi2(kxit,rule)
            wt=wantri_w(kxit,rule)
            ! Shape functions and their first derivatives
#           define delta 0.0d0
#           include <phi_and_dphidxik_2d.rc>
#           undef delta
            ! Calculate jacobian matrix, its inverse and determinant
            J=0.d0
            do i=1,n_nodes
              J(1,:)=J(1,:)+dphidxi1(i)*x_nodes(:,i)
              J(2,:)=J(2,:)+dphidxi2(i)*x_nodes(:,i)
            end do
            call fbem_invert_2x2_matrix(J,invJ,detJ)
            ! dN/dx
            do i=1,n_nodes
              dNdx(i,1)=invJ(1,1)*dphidxi1(i)+invJ(1,2)*dphidxi2(i)
              dNdx(i,2)=invJ(2,1)*dphidxi1(i)+invJ(2,2)*dphidxi2(i)
            end do
            ! Build matrix B for each node
            B=0.d0
            do i=1,n_nodes
              B(1,1,i)=dNdx(i,1)
              B(2,2,i)=dNdx(i,2)
              B(3,1,i)=dNdx(i,2)
              B(3,2,i)=dNdx(i,1)
            end do
            ! Thickness at the integration point
            t=dot_product(phi,t_nodes)
            ! det(J) * weights
            jw=detJ*t*wt
            ! Build the element stiffness matrix
            do ki=1,n_nodes
              kis=(ki-1)*2+1
              kie=kis+1
              do kj=ki,n_nodes
                kjs=(kj-1)*2+1
                kje=kjs+1
                K(kis:kie,kjs:kje)=K(kis:kie,kjs:kje)+matmul(transpose(B(:,:,ki)),matmul(Dcopy,B(:,:,kj)))*jw
              end do
            end do
          end do ! Loop through integration points
        ! QUADRILATERAL ELEMENTS
        case (4)
          ! Loop through integration points
          ! Transform to order for Wandzura rules
          do kxi1=1,gln
            ! xi_1 and w_1
            xi(1)=gl11_xi(kxi1,gln)
            w(1)=gl11_w(kxi1,gln)
            do kxi2=1,gln
              ! xi_2 and w_2
              xi(2)=gl11_xi(kxi2,gln)
              w(2)=gl11_w(kxi2,gln)
              ! Shape functions and their first derivatives
#             define delta 0.0d0
#             include <phi_and_dphidxik_2d.rc>
#             undef delta
              ! Calculate jacobian matrix, its inverse and determinant
              J=0.d0
              do i=1,n_nodes
                J(1,:)=J(1,:)+dphidxi1(i)*x_nodes(:,i)
                J(2,:)=J(2,:)+dphidxi2(i)*x_nodes(:,i)
              end do
              call fbem_invert_2x2_matrix(J,invJ,detJ)
              ! dN/dx
              do i=1,n_nodes
                dNdx(i,1)=invJ(1,1)*dphidxi1(i)+invJ(1,2)*dphidxi2(i)
                dNdx(i,2)=invJ(2,1)*dphidxi1(i)+invJ(2,2)*dphidxi2(i)
              end do
              ! Build matrix B for each node
              B=0.d0
              do i=1,n_nodes
                B(1,1,i)=dNdx(i,1)
                B(2,2,i)=dNdx(i,2)
                B(3,1,i)=dNdx(i,2)
                B(3,2,i)=dNdx(i,1)
              end do
              ! Thickness at the integration point
              t=dot_product(phi,t_nodes)
              ! det(J) * weights
              jw=detJ*t*w(1)*w(2)
              ! Build the element stiffness matrix
              do ki=1,n_nodes
                kis=(ki-1)*2+1
                kie=kis+1
                do kj=ki,n_nodes
                  kjs=(kj-1)*2+1
                  kje=kjs+1
                  K(kis:kie,kjs:kje)=K(kis:kie,kjs:kje)+matmul(transpose(B(:,:,ki)),matmul(Dcopy,B(:,:,kj)))*jw
                end do
              end do
            end do ! Loop through xi_2 coordinate
          end do ! Loop through xi_1 coordinate
        ! Unknown
        case default
          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={tri3,tri6,quad4,quad8,quad9}')
      end select
    end do ! Axial/shear contributions
    ! Copy symmetric part
    do ki=1,2*n_nodes
      do kj=ki,2*n_nodes
        K(kj,ki)=K(ki,kj)
      end do
    end do
  end subroutine fbem_fem_staela2d_solid_K

 !! Constitutive law matrix for 2D isotropic solids.
  subroutine fbem_fem_harela2d_solid_D_isotropic(E,nu,plane,D)
    implicit none
    complex(kind=real64) :: E      ! Young's modulus
    real(kind=real64)    :: nu     ! Poisson's ratio
    integer              :: plane  ! 1 (plane strain), 2 (plane_stress)
    complex(kind=real64) :: D(3,3) ! Constitutive matrix
    D=(0.d0,0.d0)
    select case (plane)
      case (1)
        D(1,1)=E*(1.d0-nu)/(1.d0+nu)/(1.d0-2.d0*nu)
        D(2,2)=D(1,1)
        D(1,2)=nu/(1.d0-nu)*D(1,1)
        D(2,1)=D(1,2)
        D(3,3)=0.5d0*E/(1.d0+nu)
      case (2)
        D(1,1)=E/(1.d0-nu**2)
        D(2,2)=D(1,1)
        D(1,2)=nu*D(1,1)
        D(2,1)=D(1,2)
        D(3,3)=0.5d0*E/(1.d0+nu)
    end select
  end subroutine fbem_fem_harela2d_solid_D_isotropic

  !! Calculate the stiffness matrix K for 2D time hormonic / modal analysis of isotropic solids.
  subroutine fbem_fem_harela2d_solid_K(etype,x_nodes,t_nodes,D,gln_a,gln_sh,K)
    implicit none
    ! I/O
    integer              :: etype                                          !! Type of element: tri3, tri6, quad4, quad8, quad9.
    real(kind=real64)    :: x_nodes(2,fbem_n_nodes(etype))                 !! Coordinates of the nodes.
    real(kind=real64)    :: t_nodes(fbem_n_nodes(etype))                   !! Thickness at each node (t=1. for 2D plane strain problem)
    complex(kind=real64) :: D(3,3)                                         !! Strain-stress constitutive matrix in global axes.
    integer              :: gln_a                                          !! Number of Gauss-Legendre integration points for axial contribution
    integer              :: gln_sh                                         !! Number of Gauss-Legendre integration points for shear contribution.
    complex(kind=real64) :: K(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype)) !! Stiffness matrix
    ! Local
    integer              :: n_nodes                                        ! Number of nodes
    integer              :: i                                              ! Counter
    integer              :: contribution                                   ! Contribution part
    integer              :: gln                                            ! Number of Gauss points
    integer              :: rule                                           ! Rule of Wandzura quadrature
    integer              :: kxi1, kxi2, kxit                               ! Integration points counters
    real(kind=real64)    :: xi(2), w(2), wt                                ! Curvilinear coordinates and quadrature weights
    real(kind=real64)    :: aux(10)                                        ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)    :: phi(fbem_n_nodes(etype))                       ! Shape functions
    real(kind=real64)    :: dphidxi1(fbem_n_nodes(etype))                  ! Shape functions derivatives with respect to xi_1
    real(kind=real64)    :: dphidxi2(fbem_n_nodes(etype))                  ! Shape functions derivatives with respect to xi_2
    real(kind=real64)    :: J(2,2), invJ(2,2), detJ                        ! Jacobian matrix, its inverse and the jacobian determinant
    real(kind=real64)    :: jw                                             ! det(J) * weights
    real(kind=real64)    :: t                                              ! Thickness at the integration point
    complex(kind=real64) :: Dcopy(3,3)                                     ! Strain-stress constitutive matrix in global axes (local copy).
    real(kind=real64)    :: dNdx(fbem_n_nodes(etype),2)                    ! dN/dx
    real(kind=real64)    :: B(3,2,fbem_n_nodes(etype))                     ! B matrix
    integer              :: ki, kis, kie, kj, kjs, kje                     ! Counters and nodal DOF limits
    ! Initialization
    K=(0.d0,0.d0)
    n_nodes=fbem_n_nodes(etype)
    ! Axial/shear contributions
    do contribution=1,2
      ! Decide how the element is going to be integrated
      Dcopy=0.d0
      if (gln_a.ne.gln_sh) then
        select case (contribution)
          case (1)
            gln=gln_a
            Dcopy(1:2,1:2)=D(1:2,1:2)
          case (2)
            gln=gln_sh
            Dcopy(1:2,3)=D(1:2,3)
            Dcopy(  3,:)=D(  3,:)
        end select
      else
        gln=gln_a
        Dcopy=D
        if (contribution.eq.2) exit
      end if
      ! Switch between types of elements
      select case (fbem_n_edges(etype))
        ! TRIANGULAR ELEMENTS
        case (3)
          ! Loop through integration points
          ! Transform to order for Wandzura rules
          rule=2*gln-1
          do kxit=1,wantri_n(rule)
            ! xi_1, xi_2 and w_t
            xi(1)=wantri_xi1(kxit,rule)
            xi(2)=wantri_xi2(kxit,rule)
            wt=wantri_w(kxit,rule)
            ! Shape functions and their first derivatives
#           define delta 0.0d0
#           include <phi_and_dphidxik_2d.rc>
#           undef delta
            ! Calculate jacobian matrix, its inverse and determinant
            J=0.d0
            do i=1,n_nodes
              J(1,:)=J(1,:)+dphidxi1(i)*x_nodes(:,i)
              J(2,:)=J(2,:)+dphidxi2(i)*x_nodes(:,i)
            end do
            call fbem_invert_2x2_matrix(J,invJ,detJ)
            ! dN/dx
            do i=1,n_nodes
              dNdx(i,1)=invJ(1,1)*dphidxi1(i)+invJ(1,2)*dphidxi2(i)
              dNdx(i,2)=invJ(2,1)*dphidxi1(i)+invJ(2,2)*dphidxi2(i)
            end do
            ! Build matrix B for each node
            B=0.d0
            do i=1,n_nodes
              B(1,1,i)=dNdx(i,1)
              B(2,2,i)=dNdx(i,2)
              B(3,1,i)=dNdx(i,2)
              B(3,2,i)=dNdx(i,1)
            end do
            ! Thickness at the integration point
            t=dot_product(phi,t_nodes)
            ! det(J) * weights
            jw=detJ*t*wt
            ! Build the element stiffness matrix
            do ki=1,n_nodes
              kis=(ki-1)*2+1
              kie=kis+1
              do kj=ki,n_nodes
                kjs=(kj-1)*2+1
                kje=kjs+1
                K(kis:kie,kjs:kje)=K(kis:kie,kjs:kje)+matmul(transpose(B(:,:,ki)),matmul(Dcopy,B(:,:,kj)))*jw
              end do
            end do
          end do ! Loop through integration points
        ! QUADRILATERAL ELEMENTS
        case (4)
          ! Loop through integration points
          ! Transform to order for Wandzura rules
          do kxi1=1,gln
            ! xi_1 and w_1
            xi(1)=gl11_xi(kxi1,gln)
            w(1)=gl11_w(kxi1,gln)
            do kxi2=1,gln
              ! xi_2 and w_2
              xi(2)=gl11_xi(kxi2,gln)
              w(2)=gl11_w(kxi2,gln)
              ! Shape functions and their first derivatives
#             define delta 0.0d0
#             include <phi_and_dphidxik_2d.rc>
#             undef delta
              ! Calculate jacobian matrix, its inverse and determinant
              J=0.d0
              do i=1,n_nodes
                J(1,:)=J(1,:)+dphidxi1(i)*x_nodes(:,i)
                J(2,:)=J(2,:)+dphidxi2(i)*x_nodes(:,i)
              end do
              call fbem_invert_2x2_matrix(J,invJ,detJ)
              ! dN/dx
              do i=1,n_nodes
                dNdx(i,1)=invJ(1,1)*dphidxi1(i)+invJ(1,2)*dphidxi2(i)
                dNdx(i,2)=invJ(2,1)*dphidxi1(i)+invJ(2,2)*dphidxi2(i)
              end do
              ! Build matrix B for each node
              B=0.d0
              do i=1,n_nodes
                B(1,1,i)=dNdx(i,1)
                B(2,2,i)=dNdx(i,2)
                B(3,1,i)=dNdx(i,2)
                B(3,2,i)=dNdx(i,1)
              end do
              ! Thickness at the integration point
              t=dot_product(phi,t_nodes)
              ! thickness * det(J) * weights
              jw=detJ*t*w(1)*w(2)
              ! Build the element stiffness matrix
              do ki=1,n_nodes
                kis=(ki-1)*2+1
                kie=kis+1
                do kj=ki,n_nodes
                  kjs=(kj-1)*2+1
                  kje=kjs+1
                  K(kis:kie,kjs:kje)=K(kis:kie,kjs:kje)+matmul(transpose(B(:,:,ki)),matmul(Dcopy,B(:,:,kj)))*jw
                end do
              end do
            end do ! Loop through xi_2 coordinate
          end do ! Loop through xi_1 coordinate
        ! Unknown
        case default
          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={tri3,tri6,quad4,quad8,quad9}')
      end select
    end do ! Axial/shear contributions
    ! Copy symmetric part
    do ki=1,2*n_nodes
      do kj=ki,2*n_nodes
        K(kj,ki)=K(ki,kj)
      end do
    end do
  end subroutine fbem_fem_harela2d_solid_K

  !! Calculate the mass matrix M for 2D solids.
  subroutine fbem_fem_ela2d_solid_M(etype,x_nodes,t_nodes,rho,gln,M)
    implicit none
    ! I/O
    integer           :: etype                                          !! Type of element: tri3, tri6, quad4, quad8, quad9.
    real(kind=real64) :: x_nodes(2,fbem_n_nodes(etype))                 !! Coordinates of the nodes.
    real(kind=real64) :: t_nodes(fbem_n_nodes(etype))                   !! Thickness at each node (t=1. for 2D plane strain problem)
    real(kind=real64) :: rho                                            !! Density
    integer           :: gln                                            !! Number of Gauss-Legendre integration points
    real(kind=real64) :: M(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype)) !! Mass matrix
    ! Local
    integer           :: n_nodes                                        ! Number of nodes
    integer           :: i                                              ! Counter
    integer           :: rule                                           ! Rule of Wandzura quadrature
    integer           :: kxi1, kxi2, kxit                               ! Integration points counters
    real(kind=real64) :: xi(2), w(2), wt                                ! Curvilinear coordinates and quadrature weights
    real(kind=real64) :: aux(10)                                        ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64) :: phi(fbem_n_nodes(etype))                       ! Shape functions
    real(kind=real64) :: dphidxi1(fbem_n_nodes(etype))                  ! Shape functions derivatives with respect to xi_1
    real(kind=real64) :: dphidxi2(fbem_n_nodes(etype))                  ! Shape functions derivatives with respect to xi_2
    real(kind=real64) :: J(2,2), detJ                                   ! Jacobian matrix and its determinant
    real(kind=real64) :: jw                                             ! det(J) * weights
    real(kind=real64) :: t                                              ! Thickness at the integration point
    real(kind=real64) :: N(2,2,fbem_n_nodes(etype))                     ! N matrix (matrix of shape functions)
    integer           :: ki, kis, kie, kj, kjs, kje                     ! Counters and nodal DOF limits
    ! Initialization
    M=0.d0
    n_nodes=fbem_n_nodes(etype)
    ! Switch between types of elements
    select case (fbem_n_edges(etype))
      ! TRIANGULAR ELEMENTS
      case (3)
        ! Loop through integration points
        ! Transform to order for Wandzura rules
        rule=2*gln-1
        do kxit=1,wantri_n(rule)
          ! xi_1, xi_2 and w_t
          xi(1)=wantri_xi1(kxit,rule)
          xi(2)=wantri_xi2(kxit,rule)
          wt=wantri_w(kxit,rule)
          ! Shape functions and their first derivatives
#         define delta 0.0d0
#         include <phi_and_dphidxik_2d.rc>
#         undef delta
          ! Calculate jacobian matrix and determinant
          J=0.d0
          do i=1,n_nodes
            J(1,:)=J(1,:)+dphidxi1(i)*x_nodes(:,i)
            J(2,:)=J(2,:)+dphidxi2(i)*x_nodes(:,i)
          end do
          detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1)
          ! Build matrix N of shape functions
          N=0.d0
          do i=1,n_nodes
            N(1,1,i)=phi(i)
            N(2,2,i)=phi(i)
          end do
          ! Thickness at the integration point
          t=dot_product(phi,t_nodes)
          ! thickness * det(J) * weights
          jw=detJ*t*wt
          ! Build the element mass matrix
          do ki=1,n_nodes
            kis=(ki-1)*2+1
            kie=kis+1
            do kj=ki,n_nodes
              kjs=(kj-1)*2+1
              kje=kjs+1
              M(kis:kie,kjs:kje)=M(kis:kie,kjs:kje)+matmul(transpose(N(:,:,ki)),N(:,:,kj))*jw
            end do
          end do
        end do ! Loop through integration points
      ! QUADRILATERAL ELEMENTS
      case (4)
        ! Loop through integration points
        ! Transform to order for Wandzura rules
        do kxi1=1,gln
          ! xi_1 and w_1
          xi(1)=gl11_xi(kxi1,gln)
          w(1)=gl11_w(kxi1,gln)
          do kxi2=1,gln
            ! xi_2 and w_2
            xi(2)=gl11_xi(kxi2,gln)
            w(2)=gl11_w(kxi2,gln)
            ! Shape functions and their first derivatives
#           define delta 0.0d0
#           include <phi_and_dphidxik_2d.rc>
#           undef delta
            ! Calculate jacobian matrix and determinant
            J=0.d0
            do i=1,n_nodes
              J(1,:)=J(1,:)+dphidxi1(i)*x_nodes(:,i)
              J(2,:)=J(2,:)+dphidxi2(i)*x_nodes(:,i)
            end do
            detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1)
            ! Build matrix N of shape functions
            N=0.d0
            do i=1,n_nodes
              N(1,1,i)=phi(i)
              N(2,2,i)=phi(i)
            end do
            ! Thickness at the integration point
            t=dot_product(phi,t_nodes)
            ! thickness * det(J) * weights
            jw=detJ*t*w(1)*w(2)
            ! Build the element mass matrix
            do ki=1,n_nodes
              kis=(ki-1)*2+1
              kie=kis+1
              do kj=ki,n_nodes
                kjs=(kj-1)*2+1
                kje=kjs+1
                M(kis:kie,kjs:kje)=M(kis:kie,kjs:kje)+matmul(transpose(N(:,:,ki)),N(:,:,kj))*jw
              end do
            end do
          end do ! Loop through xi_2 coordinate
        end do ! Loop through xi_1 coordinate
      ! Unknown
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={tri3,tri6,quad4,quad8,quad9}')
    end select
    ! Copy symmetric part
    do ki=1,2*n_nodes
      do kj=ki,2*n_nodes
        M(kj,ki)=M(ki,kj)
      end do
    end do
    M=rho*M
  end subroutine fbem_fem_ela2d_solid_M

  !! Calculate the body load matrix Q for 2D solids.
  subroutine fbem_fem_ela2d_solid_Q_body(etype,x_nodes,t_nodes,gln,Q)
    implicit none
    ! I/O
    integer           :: etype                                          !! Type of element: tri3, tri6, quad4, quad8, quad9.
    real(kind=real64) :: x_nodes(2,fbem_n_nodes(etype))                 !! Coordinates of the nodes.
    real(kind=real64) :: t_nodes(fbem_n_nodes(etype))                   !! Thickness at each node (t=1. for 2D plane strain problem)
    integer           :: gln                                            !! Number of Gauss-Legendre integration points
    real(kind=real64) :: Q(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype)) !! Body load matrix
    ! Local
    integer           :: n_nodes                                        ! Number of nodes
    integer           :: i                                              ! Counter
    integer           :: rule                                           ! Rule of Wandzura quadrature
    integer           :: kxi1, kxi2, kxit                               ! Integration points counters
    real(kind=real64) :: xi(2), w(2), wt                                ! Curvilinear coordinates and quadrature weights
    real(kind=real64) :: aux(10)                                        ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64) :: phi(fbem_n_nodes(etype))                       ! Shape functions
    real(kind=real64) :: dphidxi1(fbem_n_nodes(etype))                  ! Shape functions derivatives with respect to xi_1
    real(kind=real64) :: dphidxi2(fbem_n_nodes(etype))                  ! Shape functions derivatives with respect to xi_2
    real(kind=real64) :: J(2,2), detJ                                   ! Jacobian matrix and its determinant
    real(kind=real64) :: jw                                             ! det(J) * weights
    real(kind=real64) :: t                                              ! Thickness at the integration point
    real(kind=real64) :: N(2,2,fbem_n_nodes(etype))                     ! N matrix (matrix of shape functions)
    integer           :: ki, kis, kie, kj, kjs, kje                     ! Counters and nodal DOF limits
    ! Initialization
    Q=0.d0
    n_nodes=fbem_n_nodes(etype)
    ! Switch between types of elements
    select case (fbem_n_edges(etype))
      ! TRIANGULAR ELEMENTS
      case (3)
        ! Loop through integration points
        ! Transform to order for Wandzura rules
        rule=2*gln-1
        do kxit=1,wantri_n(rule)
          ! xi_1, xi_2 and w_t
          xi(1)=wantri_xi1(kxit,rule)
          xi(2)=wantri_xi2(kxit,rule)
          wt=wantri_w(kxit,rule)
          ! Shape functions and their first derivatives
#         define delta 0.0d0
#         include <phi_and_dphidxik_2d.rc>
#         undef delta
          ! Calculate jacobian matrix and determinant
          J=0.d0
          do i=1,n_nodes
            J(1,:)=J(1,:)+dphidxi1(i)*x_nodes(:,i)
            J(2,:)=J(2,:)+dphidxi2(i)*x_nodes(:,i)
          end do
          detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1)
          ! Build matrix N of shape functions
          N=0.d0
          do i=1,n_nodes
            N(1,1,i)=phi(i)
            N(2,2,i)=phi(i)
          end do
          ! Thickness at the integration point
          t=dot_product(phi,t_nodes)
          ! thickness * det(J) * weights
          jw=detJ*t*wt
          ! Build the element mass matrix
          do ki=1,n_nodes
            kis=(ki-1)*2+1
            kie=kis+1
            do kj=ki,n_nodes
              kjs=(kj-1)*2+1
              kje=kjs+1
              Q(kis:kie,kjs:kje)=Q(kis:kie,kjs:kje)+matmul(transpose(N(:,:,ki)),N(:,:,kj))*jw
            end do
          end do
        end do ! Loop through integration points
      ! QUADRILATERAL ELEMENTS
      case (4)
        ! Loop through integration points
        ! Transform to order for Wandzura rules
        do kxi1=1,gln
          ! xi_1 and w_1
          xi(1)=gl11_xi(kxi1,gln)
          w(1)=gl11_w(kxi1,gln)
          do kxi2=1,gln
            ! xi_2 and w_2
            xi(2)=gl11_xi(kxi2,gln)
            w(2)=gl11_w(kxi2,gln)
            ! Shape functions and their first derivatives
#           define delta 0.0d0
#           include <phi_and_dphidxik_2d.rc>
#           undef delta
            ! Calculate jacobian matrix and determinant
            J=0.d0
            do i=1,n_nodes
              J(1,:)=J(1,:)+dphidxi1(i)*x_nodes(:,i)
              J(2,:)=J(2,:)+dphidxi2(i)*x_nodes(:,i)
            end do
            detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1)
            ! Build matrix N of shape functions
            N=0.d0
            do i=1,n_nodes
              N(1,1,i)=phi(i)
              N(2,2,i)=phi(i)
            end do
            ! Thickness at the integration point
            t=dot_product(phi,t_nodes)
            ! thickness * det(J) * weights
            jw=detJ*t*w(1)*w(2)
            ! Build the element mass matrix
            do ki=1,n_nodes
              kis=(ki-1)*2+1
              kie=kis+1
              do kj=ki,n_nodes
                kjs=(kj-1)*2+1
                kje=kjs+1
                Q(kis:kie,kjs:kje)=Q(kis:kie,kjs:kje)+matmul(transpose(N(:,:,ki)),N(:,:,kj))*jw
              end do
            end do
          end do ! Loop through xi_2 coordinate
        end do ! Loop through xi_1 coordinate
      ! Unknown
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={tri3,tri6,quad4,quad8,quad9}')
    end select
    ! Copy symmetric part
    do ki=1,2*n_nodes
      do kj=ki,2*n_nodes
        Q(kj,ki)=Q(ki,kj)
      end do
    end do
  end subroutine fbem_fem_ela2d_solid_Q_body

  !! Load matrix Q for an edge distributed load
  subroutine fbem_fem_ela2d_solid_Q_edge(etype,x_nodes,t_nodes,kedge,dtype,ngp,Q)
    implicit none
    ! I/O
    integer           :: etype                                          !! Type of element: tri3, tri6, quad4, quad8, quad9.
    real(kind=real64) :: x_nodes(2,fbem_n_nodes(etype))                 !! Position vectors of the mid-plane nodes.
    real(kind=real64) :: t_nodes(fbem_n_nodes(etype))                   !! Thickness of each mid-node in the v_3 direction.
    integer           :: kedge                                          !! Selected edge
    integer           :: dtype                                          !! Type of edge load element: point1, line2, line3.
    integer           :: ngp                                            !! Number of Gauss points for edge coordinate integration
    real(kind=real64) :: Q(2*fbem_n_nodes(etype),2*fbem_n_nodes(dtype)) !! Load matrix Q for a global coordinates distributed load over an edge
    ! Local
    integer           :: n_nodes                              ! Number of nodes of the midplane element
    integer           :: n_dnodes                             ! Number of nodes of the load element
    integer           :: k_node                               ! Counter of mid-nodes
    integer           :: kp                                   ! Integration points counter
    integer           :: ki                                   ! Counter
    real(kind=real64) :: xip                                  ! xip coordinate
    real(kind=real64) :: w                                    ! Weight
    real(kind=real64) :: xi_e(2,2)                            ! xi coordinates of the edge
    real(kind=real64) :: xi_nodes(2,fbem_n_nodes(etype))      ! xi1,xi2 coordinates of the nodes of the element
    real(kind=real64) :: dxidxip(2)                           ! xip->xi1,xi2 tangent
    real(kind=real64) :: xi(2)                                ! xi coordinates
    real(kind=real64) :: aux(10)                              ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64) :: phi(fbem_n_nodes(etype))             ! phi shape functions
    real(kind=real64) :: dphidxi1(fbem_n_nodes(etype))        ! phi shape functions derivatives with respect to xi_1
    real(kind=real64) :: dphidxi2(fbem_n_nodes(etype))        ! phi shape functions derivatives with respect to xi_2
    real(kind=real64) :: dxdxi1(2), dxdxi2(2)                 ! dx/dxik
    real(kind=real64) :: dxdxip(2)                            ! dx/dxip
    real(kind=real64) :: t                                    ! Thickness at the integration point
    real(kind=real64) :: jg                                   ! Area jacobian
    real(kind=real64) :: jw                                   ! jacobian * weight
    real(kind=real64) :: Ne(2,2*fbem_n_nodes(etype))          ! Element shape functions matrix
    real(kind=real64) :: Nedge(2,2*fbem_n_nodes(dtype))       ! Edge load shape functions matrix
    ! Check
    if (kedge.gt.fbem_n_edges(etype)) then
      call fbem_error_message(error_unit,0,'fbem_fem_ela2d_solid_Q_edge',kedge,'invalid edge index')
    end if
    ! Initialize load matrix
    Q=0.d0
    ! Number of nodes of the midplane element
    n_nodes=fbem_n_nodes(etype)
    ! Number of nodes of the load element
    n_dnodes=fbem_n_nodes(dtype)
    ! xi1,xi2 coordinates of kedge
    xi_nodes=fbem_xi1xi2(etype,0.d0)
    xi_e(:,1)=xi_nodes(:,fbem_edge_node(1,kedge,etype))
    xi_e(:,2)=xi_nodes(:,fbem_edge_node(2,kedge,etype))
    ! Loop through xip
    do kp=1,gl11_n(ngp)
      ! XIP COORDINATE AND WEIGHT
      xip=gl11_xi(kp,ngp)
      w=gl11_w(kp,ngp)
      ! XIP->XI(1),XI(2) TRANSFORMATION
      xi=0.5d0*(1.d0-xip)*xi_e(:,1)+0.5d0*(1.d0+xip)*xi_e(:,2)
      dxidxip=0.5d0*(xi_e(:,2)-xi_e(:,1))
      ! Solid element shape functions and its derivatives
#     define delta 0.0d0
#     include <phi_and_dphidxik_2d.rc>
#     undef delta
      ! XIP -> XI -> X TRANSFORMATIONS
      dxdxi1=0.d0
      dxdxi2=0.d0
      do k_node=1,n_nodes
        dxdxi1=dxdxi1+dphidxi1(k_node)*x_nodes(:,k_node)
        dxdxi2=dxdxi2+dphidxi2(k_node)*x_nodes(:,k_node)
      end do
      dxdxip=dxdxi1*dxidxip(1)+dxdxi2*dxidxip(2)
      jg=sqrt(dot_product(dxdxip,dxdxip))
      ! Thickness at the integration point
      t=dot_product(phi,t_nodes)
      ! Build element shape function matrix Ne
      Ne=0.d0
      do k_node=1,n_nodes
        ki=2*(k_node-1)
        Ne(1,ki+1)=phi(k_node)
        Ne(2,ki+2)=phi(k_node)
      end do
      ! Edge load shape functions
      phi=0.d0
#     define etype dtype
#     define delta 0.0d0
#     define xi xip
#     include <phi_1d.rc>
#     undef etype
#     undef delta
#     undef xi
      ! Build load shape function matrix Nedge
      Nedge=0.d0
      do k_node=1,n_dnodes
        ki=2*(k_node-1)
        Nedge(1,ki+1)=phi(k_node)
        Nedge(2,ki+2)=phi(k_node)
      end do
      ! hhickness * jacobian * weight
      jw=t*jg*w
      ! Add
      Q=Q+matmul(transpose(Ne),Nedge)*jw
    end do
  end subroutine fbem_fem_ela2d_solid_Q_edge

  subroutine fbem_fem_staela2d_solid_dKdx(etype,x_nodes,t_nodes,D,gln_a,gln_sh,dKdx)
    implicit none
    ! I/O
    integer           :: etype                                                                   !! Type of element: tri3, tri6, quad4, quad8, quad9.
    real(kind=real64) :: x_nodes(2,fbem_n_nodes(etype))                                          !! Coordinates of the nodes.
    real(kind=real64) :: t_nodes(fbem_n_nodes(etype))                                            !! Thickness at each node (t=1. for 2D plane strain problem)
    real(kind=real64) :: D(3,3)                                                                  !! Strain-stress constitutive matrix in global axes.
    integer           :: gln_a                                                                   !! Number of Gauss-Legendre integration points for axial contribution
    integer           :: gln_sh                                                                  !! Number of Gauss-Legendre integration points for shear contribution.
    real(kind=real64) :: dKdx(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype),2,fbem_n_nodes(etype)) !! Stiffness matrix sensitivity dK/dx_k^i
    ! Local
    integer           :: n_nodes                                        ! Number of nodes
    integer           :: i                                              ! Counter
    integer           :: contribution                                   ! Contribution part
    integer           :: gln                                            ! Number of Gauss points
    integer           :: rule                                           ! Rule of Wandzura quadrature
    integer           :: kxi1, kxi2, kxit                               ! Integration points counters
    real(kind=real64) :: xi(2), w(2), wt                                ! Curvilinear coordinates and quadrature weights
    real(kind=real64) :: aux(10)                                        ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64) :: phi(fbem_n_nodes(etype))                       ! Shape functions
    real(kind=real64) :: dphidxi1(fbem_n_nodes(etype))                  ! Shape functions derivatives with respect to xi_1
    real(kind=real64) :: dphidxi2(fbem_n_nodes(etype))                  ! Shape functions derivatives with respect to xi_2
    real(kind=real64) :: J(2,2), invJ(2,2), detJ                        ! Jacobian matrix, its inverse and the jacobian determinant
    real(kind=real64) :: jw                                             ! det(J) * weights
    real(kind=real64) :: t                                              ! Thickness at the integration point
    real(kind=real64) :: Dcopy(3,3)                                     ! Strain-stress constitutive matrix in global axes (local copy).
    real(kind=real64) :: dNdx(fbem_n_nodes(etype),2)                    ! dN/dx
    real(kind=real64) :: d1JdJdxjk                                      ! 1/|J| d|J|/dx_j^k
    real(kind=real64) :: B(3,2,fbem_n_nodes(etype))                     ! B matrix
    real(kind=real64) :: dBdxjk(3,2,fbem_n_nodes(etype))                ! dB/dx_j^k
    real(kind=real64) :: hatB(3,2,fbem_n_nodes(etype))                  ! har{B} matrix
    integer           :: ki, kis, kie, kj, kjs, kje, kdj, kdk           ! Counters and nodal DOF limits
    ! Initialization
    dKdx=0.d0
    n_nodes=fbem_n_nodes(etype)
    ! Axial/shear contributions
    do contribution=1,2
      ! Decide how the element is going to be integrated
      Dcopy=0.d0
      if (gln_a.ne.gln_sh) then
        select case (contribution)
          case (1)
            gln=gln_a
            Dcopy(1:2,1:2)=D(1:2,1:2)
          case (2)
            gln=gln_sh
            Dcopy(1:2,3)=D(1:2,3)
            Dcopy(  3,:)=D(  3,:)
        end select
      else
        gln=gln_a
        Dcopy=D
        if (contribution.eq.2) exit
      end if
      ! Switch between types of elements
      select case (fbem_n_edges(etype))
        ! TRIANGULAR ELEMENTS
        case (3)
          ! Loop through integration points
          ! Transform to order for Wandzura rules
          rule=2*gln-1
          do kxit=1,wantri_n(rule)
            ! xi_1, xi_2 and w_t
            xi(1)=wantri_xi1(kxit,rule)
            xi(2)=wantri_xi2(kxit,rule)
            wt=wantri_w(kxit,rule)
            ! Shape functions and their first derivatives
#           define delta 0.0d0
#           include <phi_and_dphidxik_2d.rc>
#           undef delta
            ! Calculate jacobian matrix, its inverse and determinant
            J=0.d0
            do i=1,n_nodes
              J(1,:)=J(1,:)+dphidxi1(i)*x_nodes(:,i)
              J(2,:)=J(2,:)+dphidxi2(i)*x_nodes(:,i)
            end do
            call fbem_invert_2x2_matrix(J,invJ,detJ)
            ! dN/dx
            do i=1,n_nodes
              dNdx(i,1)=invJ(1,1)*dphidxi1(i)+invJ(1,2)*dphidxi2(i)
              dNdx(i,2)=invJ(2,1)*dphidxi1(i)+invJ(2,2)*dphidxi2(i)
            end do
            ! Build matrix B for each node
            B=0.d0
            do i=1,n_nodes
              B(1,1,i)=dNdx(i,1)
              B(2,2,i)=dNdx(i,2)
              B(3,1,i)=dNdx(i,2)
              B(3,2,i)=dNdx(i,1)
            end do
            ! Thickness at the integration point
            t=dot_product(phi,t_nodes)
            ! det(J) * weights
            jw=detJ*t*wt
            ! Integrate dK/dx_j^k
            do kdj=1,2
              do kdk=1,n_nodes
                ! 1/|J|*d|J|/dx_j^k
                d1JdJdxjk=dNdx(kdk,kdj)
                ! Build matrix dB/dx_j^k for each node
                dBdxjk=0.d0
                do i=1,n_nodes
                  dBdxjk(1,1,i)=-dNdx(i,kdj)*dNdx(kdk,1)
                  dBdxjk(2,2,i)=-dNdx(i,kdj)*dNdx(kdk,2)
                  dBdxjk(3,1,i)=-dNdx(i,kdj)*dNdx(kdk,2)
                  dBdxjk(3,2,i)=-dNdx(i,kdj)*dNdx(kdk,1)
                end do
                ! hat{B}
                hatB=2.d0*dBdxjk+d1JdJdxjk*B
                ! Build the element stiffness matrix
                do ki=1,n_nodes
                  kis=(ki-1)*2+1
                  kie=kis+1
                  do kj=1,n_nodes
                    kjs=(kj-1)*2+1
                    kje=kjs+1
                    dKdx(kis:kie,kjs:kje,kdj,kdk)=dKdx(kis:kie,kjs:kje,kdj,kdk)+matmul(transpose(B(:,:,ki)),matmul(Dcopy,hatB(:,:,kj)))*jw
                  end do
                end do
              end do
            end do
          end do ! Loop through integration points
        ! QUADRILATERAL ELEMENTS
        case (4)
          ! Loop through integration points
          ! Transform to order for Wandzura rules
          do kxi1=1,gln
            ! xi_1 and w_1
            xi(1)=gl11_xi(kxi1,gln)
            w(1)=gl11_w(kxi1,gln)
            do kxi2=1,gln
              ! xi_2 and w_2
              xi(2)=gl11_xi(kxi2,gln)
              w(2)=gl11_w(kxi2,gln)
              ! Shape functions and their first derivatives
#             define delta 0.0d0
#             include <phi_and_dphidxik_2d.rc>
#             undef delta
              ! Calculate jacobian matrix, its inverse and determinant
              J=0.d0
              do i=1,n_nodes
                J(1,:)=J(1,:)+dphidxi1(i)*x_nodes(:,i)
                J(2,:)=J(2,:)+dphidxi2(i)*x_nodes(:,i)
              end do
              call fbem_invert_2x2_matrix(J,invJ,detJ)
              ! dN/dx
              do i=1,n_nodes
                dNdx(i,1)=invJ(1,1)*dphidxi1(i)+invJ(1,2)*dphidxi2(i)
                dNdx(i,2)=invJ(2,1)*dphidxi1(i)+invJ(2,2)*dphidxi2(i)
              end do
              ! Build matrix B for each node
              B=0.d0
              do i=1,n_nodes
                B(1,1,i)=dNdx(i,1)
                B(2,2,i)=dNdx(i,2)
                B(3,1,i)=dNdx(i,2)
                B(3,2,i)=dNdx(i,1)
              end do
              ! Thickness at the integration point
              t=dot_product(phi,t_nodes)
              ! det(J) * weights
              jw=detJ*t*w(1)*w(2)
              ! Integrate dK/dx_j^k
              do kdj=1,2
                do kdk=1,n_nodes
                  ! 1/|J|*d|J|/dx_j^k
                  d1JdJdxjk=dNdx(kdk,kdj)
                  ! Build matrix dB/dx_j^k for each node
                  dBdxjk=0.d0
                  do i=1,n_nodes
                    dBdxjk(1,1,i)=-dNdx(i,kdj)*dNdx(kdk,1)
                    dBdxjk(2,2,i)=-dNdx(i,kdj)*dNdx(kdk,2)
                    dBdxjk(3,1,i)=-dNdx(i,kdj)*dNdx(kdk,2)
                    dBdxjk(3,2,i)=-dNdx(i,kdj)*dNdx(kdk,1)
                  end do
                  ! hat{B}
                  hatB=2.d0*dBdxjk+d1JdJdxjk*B
                  ! Build the element stiffness matrix
                  do ki=1,n_nodes
                    kis=(ki-1)*2+1
                    kie=kis+1
                    do kj=1,n_nodes
                      kjs=(kj-1)*2+1
                      kje=kjs+1
                      dKdx(kis:kie,kjs:kje,kdj,kdk)=dKdx(kis:kie,kjs:kje,kdj,kdk)+matmul(transpose(B(:,:,ki)),matmul(Dcopy,hatB(:,:,kj)))*jw
                    end do
                  end do
                end do
              end do
            end do ! Loop through xi_2 coordinate
          end do ! Loop through xi_1 coordinate
        ! Unknown
        case default
          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={tri3,tri6,quad4,quad8,quad9}')
      end select
    end do ! Axial/shear contributions
    ! Built K = 1/2*[(B^T·D·hatB)+(B^T·D·hatB)^T]
    do kdj=1,2
      do kdk=1,n_nodes
        dKdx(:,:,kdj,kdk)=dKdx(:,:,kdj,kdk)+transpose(dKdx(:,:,kdj,kdk))
      end do
    end do
    dKdx=0.5d0*dKdx
  end subroutine fbem_fem_staela2d_solid_dKdx

  !! Calculate the stiffness matrix sensitivity of nodal coordinates dK/dx_k^i for 2D static analysis of isotropic solids. It is
  !! calculated using Finite Differences (4-points) with a deltax=1.e-4.
  subroutine fbem_fem_staela2d_solid_dKdx_FD(etype,x_nodes,t_nodes,D,gln_a,gln_sh,dKdx)
    implicit none
    ! I/O
    integer           :: etype                                                                   !! Type of element: tri3, tri6, quad4, quad8, quad9.
    real(kind=real64) :: x_nodes(2,fbem_n_nodes(etype))                                          !! Coordinates of the nodes.
    real(kind=real64) :: t_nodes(fbem_n_nodes(etype))                                            !! Thickness at each node (t=1. for 2D plane strain problem)
    real(kind=real64) :: D(3,3)                                                                  !! Strain-stress constitutive matrix in global axes.
    integer           :: gln_a                                                                   !! Number of Gauss-Legendre integration points for axial contribution
    integer           :: gln_sh                                                                  !! Number of Gauss-Legendre integration points for shear contribution.
    real(kind=real64) :: dKdx(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype),2,fbem_n_nodes(etype)) !! Stiffness matrix sensitivity dK/dx_k^i
    ! Local
    real(kind=real64) :: cl
    real(kind=real64) :: deltax
    real(kind=real64) :: x_nodes_FD(2,fbem_n_nodes(etype))
    real(kind=real64) :: K1(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    real(kind=real64) :: K2(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    real(kind=real64) :: K3(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    real(kind=real64) :: K4(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    integer           :: i, j
    ! Characteristic length
    cl=fbem_characteristic_length(2,etype,x_nodes,1.d-6)
    deltax=1.d-4*cl
    do i=1,fbem_n_nodes(etype)
      do j=1,2
        x_nodes_FD=x_nodes
        x_nodes_FD(j,i)=x_nodes_FD(j,i)+2.d0*deltax
        call fbem_fem_staela2d_solid_K(etype,x_nodes_FD,t_nodes,D,gln_a,gln_sh,K1)
        x_nodes_FD=x_nodes
        x_nodes_FD(j,i)=x_nodes_FD(j,i)+deltax
        call fbem_fem_staela2d_solid_K(etype,x_nodes_FD,t_nodes,D,gln_a,gln_sh,K2)
        x_nodes_FD=x_nodes
        x_nodes_FD(j,i)=x_nodes_FD(j,i)-deltax
        call fbem_fem_staela2d_solid_K(etype,x_nodes_FD,t_nodes,D,gln_a,gln_sh,K3)
        x_nodes_FD=x_nodes
        x_nodes_FD(j,i)=x_nodes_FD(j,i)-2.d0*deltax
        call fbem_fem_staela2d_solid_K(etype,x_nodes_FD,t_nodes,D,gln_a,gln_sh,K4)
        dKdx(:,:,j,i)=(K4-K1+8.d0*(K2-K3))/(12.d0*deltax)
      end do
    end do
  end subroutine fbem_fem_staela2d_solid_dKdx_FD

  subroutine fbem_fem_harela2d_solid_dKdx(etype,x_nodes,t_nodes,D,gln_a,gln_sh,dKdx)
    implicit none
    ! I/O
    integer              :: etype                                                                   !! Type of element: tri3, tri6, quad4, quad8, quad9.
    real(kind=real64)    :: x_nodes(2,fbem_n_nodes(etype))                                          !! Coordinates of the nodes.
    real(kind=real64)    :: t_nodes(fbem_n_nodes(etype))                                            !! Thickness at each node (t=1. for 2D plane strain problem)
    complex(kind=real64) :: D(3,3)                                                                  !! Strain-stress constitutive matrix in global axes.
    integer              :: gln_a                                                                   !! Number of Gauss-Legendre integration points for axial contribution
    integer              :: gln_sh                                                                  !! Number of Gauss-Legendre integration points for shear contribution.
    complex(kind=real64) :: dKdx(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype),2,fbem_n_nodes(etype)) !! Stiffness matrix sensitivity dK/dx_k^i
    ! Local
    integer              :: n_nodes                                        ! Number of nodes
    integer              :: i                                              ! Counter
    integer              :: contribution                                   ! Contribution part
    integer              :: gln                                            ! Number of Gauss points
    integer              :: rule                                           ! Rule of Wandzura quadrature
    integer              :: kxi1, kxi2, kxit                               ! Integration points counters
    real(kind=real64)    :: xi(2), w(2), wt                                ! Curvilinear coordinates and quadrature weights
    real(kind=real64)    :: aux(10)                                        ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)    :: phi(fbem_n_nodes(etype))                       ! Shape functions
    real(kind=real64)    :: dphidxi1(fbem_n_nodes(etype))                  ! Shape functions derivatives with respect to xi_1
    real(kind=real64)    :: dphidxi2(fbem_n_nodes(etype))                  ! Shape functions derivatives with respect to xi_2
    real(kind=real64)    :: J(2,2), invJ(2,2), detJ                        ! Jacobian matrix, its inverse and the jacobian determinant
    real(kind=real64)    :: jw                                             ! det(J) * weights
    real(kind=real64)    :: t                                              ! Thickness at the integration point
    complex(kind=real64) :: Dcopy(3,3)                                     ! Strain-stress constitutive matrix in global axes (local copy).
    real(kind=real64)    :: dNdx(fbem_n_nodes(etype),2)                    ! dN/dx
    real(kind=real64)    :: d1JdJdxjk                                      ! 1/|J| d|J|/dx_j^k
    real(kind=real64)    :: B(3,2,fbem_n_nodes(etype))                     ! B matrix
    real(kind=real64)    :: dBdxjk(3,2,fbem_n_nodes(etype))                ! dB/dx_j^k
    real(kind=real64)    :: hatB(3,2,fbem_n_nodes(etype))                  ! har{B} matrix
    integer              :: ki, kis, kie, kj, kjs, kje, kdj, kdk           ! Counters and nodal DOF limits
    ! Initialization
    dKdx=(0.d0,0.d0)
    n_nodes=fbem_n_nodes(etype)
    ! Axial/shear contributions
    do contribution=1,2
      ! Decide how the element is going to be integrated
      Dcopy=0.d0
      if (gln_a.ne.gln_sh) then
        select case (contribution)
          case (1)
            gln=gln_a
            Dcopy(1:2,1:2)=D(1:2,1:2)
          case (2)
            gln=gln_sh
            Dcopy(1:2,3)=D(1:2,3)
            Dcopy(  3,:)=D(  3,:)
        end select
      else
        gln=gln_a
        Dcopy=D
        if (contribution.eq.2) exit
      end if
      ! Switch between types of elements
      select case (fbem_n_edges(etype))
        ! TRIANGULAR ELEMENTS
        case (3)
          ! Loop through integration points
          ! Transform to order for Wandzura rules
          rule=2*gln-1
          do kxit=1,wantri_n(rule)
            ! xi_1, xi_2 and w_t
            xi(1)=wantri_xi1(kxit,rule)
            xi(2)=wantri_xi2(kxit,rule)
            wt=wantri_w(kxit,rule)
            ! Shape functions and their first derivatives
#           define delta 0.0d0
#           include <phi_and_dphidxik_2d.rc>
#           undef delta
            ! Calculate jacobian matrix, its inverse and determinant
            J=0.d0
            do i=1,n_nodes
              J(1,:)=J(1,:)+dphidxi1(i)*x_nodes(:,i)
              J(2,:)=J(2,:)+dphidxi2(i)*x_nodes(:,i)
            end do
            call fbem_invert_2x2_matrix(J,invJ,detJ)
            ! dN/dx
            do i=1,n_nodes
              dNdx(i,1)=invJ(1,1)*dphidxi1(i)+invJ(1,2)*dphidxi2(i)
              dNdx(i,2)=invJ(2,1)*dphidxi1(i)+invJ(2,2)*dphidxi2(i)
            end do
            ! Build matrix B for each node
            B=0.d0
            do i=1,n_nodes
              B(1,1,i)=dNdx(i,1)
              B(2,2,i)=dNdx(i,2)
              B(3,1,i)=dNdx(i,2)
              B(3,2,i)=dNdx(i,1)
            end do
            ! Thickness at the integration point
            t=dot_product(phi,t_nodes)
            ! det(J) * weights
            jw=detJ*t*wt
            ! Integrate dK/dx_j^k
            do kdj=1,2
              do kdk=1,n_nodes
                ! 1/|J|*d|J|/dx_j^k
                d1JdJdxjk=dNdx(kdk,kdj)
                ! Build matrix dB/dx_j^k for each node
                dBdxjk=0.d0
                do i=1,n_nodes
                  dBdxjk(1,1,i)=-dNdx(i,kdj)*dNdx(kdk,1)
                  dBdxjk(2,2,i)=-dNdx(i,kdj)*dNdx(kdk,2)
                  dBdxjk(3,1,i)=-dNdx(i,kdj)*dNdx(kdk,2)
                  dBdxjk(3,2,i)=-dNdx(i,kdj)*dNdx(kdk,1)
                end do
                ! hat{B}
                hatB=2.d0*dBdxjk+d1JdJdxjk*B
                ! Build the element stiffness matrix
                do ki=1,n_nodes
                  kis=(ki-1)*2+1
                  kie=kis+1
                  do kj=1,n_nodes
                    kjs=(kj-1)*2+1
                    kje=kjs+1
                    dKdx(kis:kie,kjs:kje,kdj,kdk)=dKdx(kis:kie,kjs:kje,kdj,kdk)+matmul(transpose(B(:,:,ki)),matmul(Dcopy,hatB(:,:,kj)))*jw
                  end do
                end do
              end do
            end do
          end do ! Loop through integration points
        ! QUADRILATERAL ELEMENTS
        case (4)
          ! Loop through integration points
          ! Transform to order for Wandzura rules
          do kxi1=1,gln
            ! xi_1 and w_1
            xi(1)=gl11_xi(kxi1,gln)
            w(1)=gl11_w(kxi1,gln)
            do kxi2=1,gln
              ! xi_2 and w_2
              xi(2)=gl11_xi(kxi2,gln)
              w(2)=gl11_w(kxi2,gln)
              ! Shape functions and their first derivatives
#             define delta 0.0d0
#             include <phi_and_dphidxik_2d.rc>
#             undef delta
              ! Calculate jacobian matrix, its inverse and determinant
              J=0.d0
              do i=1,n_nodes
                J(1,:)=J(1,:)+dphidxi1(i)*x_nodes(:,i)
                J(2,:)=J(2,:)+dphidxi2(i)*x_nodes(:,i)
              end do
              call fbem_invert_2x2_matrix(J,invJ,detJ)
              ! dN/dx
              do i=1,n_nodes
                dNdx(i,1)=invJ(1,1)*dphidxi1(i)+invJ(1,2)*dphidxi2(i)
                dNdx(i,2)=invJ(2,1)*dphidxi1(i)+invJ(2,2)*dphidxi2(i)
              end do
              ! Build matrix B for each node
              B=0.d0
              do i=1,n_nodes
                B(1,1,i)=dNdx(i,1)
                B(2,2,i)=dNdx(i,2)
                B(3,1,i)=dNdx(i,2)
                B(3,2,i)=dNdx(i,1)
              end do
              ! Thickness at the integration point
              t=dot_product(phi,t_nodes)
              ! det(J) * weights
              jw=detJ*t*w(1)*w(2)
              ! Integrate dK/dx_j^k
              do kdj=1,2
                do kdk=1,n_nodes
                  ! 1/|J|*d|J|/dx_j^k
                  d1JdJdxjk=dNdx(kdk,kdj)
                  ! Build matrix dB/dx_j^k for each node
                  dBdxjk=0.d0
                  do i=1,n_nodes
                    dBdxjk(1,1,i)=-dNdx(i,kdj)*dNdx(kdk,1)
                    dBdxjk(2,2,i)=-dNdx(i,kdj)*dNdx(kdk,2)
                    dBdxjk(3,1,i)=-dNdx(i,kdj)*dNdx(kdk,2)
                    dBdxjk(3,2,i)=-dNdx(i,kdj)*dNdx(kdk,1)
                  end do
                  ! hat{B}
                  hatB=2.d0*dBdxjk+d1JdJdxjk*B
                  ! Build the element stiffness matrix
                  do ki=1,n_nodes
                    kis=(ki-1)*2+1
                    kie=kis+1
                    do kj=1,n_nodes
                      kjs=(kj-1)*2+1
                      kje=kjs+1
                      dKdx(kis:kie,kjs:kje,kdj,kdk)=dKdx(kis:kie,kjs:kje,kdj,kdk)+matmul(transpose(B(:,:,ki)),matmul(Dcopy,hatB(:,:,kj)))*jw
                    end do
                  end do
                end do
              end do
            end do ! Loop through xi_2 coordinate
          end do ! Loop through xi_1 coordinate
        ! Unknown
        case default
          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={tri3,tri6,quad4,quad8,quad9}')
      end select
    end do ! Axial/shear contributions
    ! Built K = 1/2*[(B^T·D·hatB)+(B^T·D·hatB)^T]
    do kdj=1,2
      do kdk=1,n_nodes
        dKdx(:,:,kdj,kdk)=dKdx(:,:,kdj,kdk)+transpose(dKdx(:,:,kdj,kdk))
      end do
    end do
    dKdx=0.5d0*dKdx
  end subroutine fbem_fem_harela2d_solid_dKdx

  !! Calculate the stiffness matrix sensitivity of nodal coordinates dK/dx_k^i for 2D dynamic analysis of isotropic solids. It is
  !! calculated using Finite Differences (4-points) with a deltax=1.e-4.
  subroutine fbem_fem_harela2d_solid_dKdx_FD(etype,x_nodes,t_nodes,D,gln_a,gln_sh,dKdx)
    implicit none
    ! I/O
    integer              :: etype                                                                   !! Type of element: tri3, tri6, quad4, quad8, quad9.
    real(kind=real64)    :: x_nodes(2,fbem_n_nodes(etype))                                          !! Coordinates of the nodes.
    real(kind=real64)    :: t_nodes(fbem_n_nodes(etype))                                            !! Thickness at each node (t=1. for 2D plane strain problem)
    complex(kind=real64) :: D(3,3)                                                                  !! Strain-stress constitutive matrix in global axes.
    integer              :: gln_a                                                                   !! Number of Gauss-Legendre integration points for axial contribution
    integer              :: gln_sh                                                                  !! Number of Gauss-Legendre integration points for shear contribution.
    complex(kind=real64) :: dKdx(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype),2,fbem_n_nodes(etype)) !! Stiffness matrix sensitivity dK/dx_k^i
    ! Local
    real(kind=real64)    :: cl
    real(kind=real64)    :: deltax
    real(kind=real64)    :: x_nodes_FD(2,fbem_n_nodes(etype))
    complex(kind=real64) :: K1(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    complex(kind=real64) :: K2(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    complex(kind=real64) :: K3(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    complex(kind=real64) :: K4(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    integer              :: i, j
!    if (etype.eq.fbem_tri3) then
!      call fbem_fem_harela2d_solid_dKdx_tri3(x_nodes,t_nodes,D,dKdx)
!      return
!    end if
    ! Characteristic length
    cl=fbem_characteristic_length(2,etype,x_nodes,1.d-6)
    deltax=1.d-4*cl
    do i=1,fbem_n_nodes(etype)
      do j=1,2
        x_nodes_FD=x_nodes
        x_nodes_FD(j,i)=x_nodes_FD(j,i)+2.d0*deltax
        call fbem_fem_harela2d_solid_K(etype,x_nodes_FD,t_nodes,D,gln_a,gln_sh,K1)
        x_nodes_FD=x_nodes
        x_nodes_FD(j,i)=x_nodes_FD(j,i)+deltax
        call fbem_fem_harela2d_solid_K(etype,x_nodes_FD,t_nodes,D,gln_a,gln_sh,K2)
        x_nodes_FD=x_nodes
        x_nodes_FD(j,i)=x_nodes_FD(j,i)-deltax
        call fbem_fem_harela2d_solid_K(etype,x_nodes_FD,t_nodes,D,gln_a,gln_sh,K3)
        x_nodes_FD=x_nodes
        x_nodes_FD(j,i)=x_nodes_FD(j,i)-2.d0*deltax
        call fbem_fem_harela2d_solid_K(etype,x_nodes_FD,t_nodes,D,gln_a,gln_sh,K4)
        dKdx(:,:,j,i)=(K4-K1+8.d0*(K2-K3))/(12.d0*deltax)
      end do
    end do
  end subroutine fbem_fem_harela2d_solid_dKdx_FD

  subroutine fbem_fem_ela2d_solid_dMdx(etype,x_nodes,t_nodes,rho,gln,dMdx)
    implicit none
    ! I/O
    integer           :: etype                                                                   !! Type of element: tri3, tri6, quad4, quad8, quad9.
    real(kind=real64) :: x_nodes(2,fbem_n_nodes(etype))                                          !! Coordinates of the nodes.
    real(kind=real64) :: t_nodes(fbem_n_nodes(etype))                                            !! Thickness at each node (t=1. for 2D plane strain problem)
    real(kind=real64) :: rho                                                                     !! Density
    integer           :: gln                                                                     !! Number of Gauss-Legendre integration points
    real(kind=real64) :: dMdx(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype),2,fbem_n_nodes(etype)) !! Mass matrix sensitivity dM/dx_k^i
    ! Local
    integer           :: n_nodes                                        ! Number of nodes
    integer           :: i                                              ! Counter
    integer           :: rule                                           ! Rule of Wandzura quadrature
    integer           :: kxi1, kxi2, kxit                               ! Integration points counters
    real(kind=real64) :: xi(2), w(2), wt                                ! Curvilinear coordinates and quadrature weights
    real(kind=real64) :: aux(10)                                        ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64) :: phi(fbem_n_nodes(etype))                       ! Shape functions
    real(kind=real64) :: dphidxi1(fbem_n_nodes(etype))                  ! Shape functions derivatives with respect to xi_1
    real(kind=real64) :: dphidxi2(fbem_n_nodes(etype))                  ! Shape functions derivatives with respect to xi_2
    real(kind=real64) :: J(2,2), invJ(2,2), detJ                        ! Jacobian matrix and its determinant
    real(kind=real64) :: jw                                             ! det(J) * weights
    real(kind=real64) :: t                                              ! Thickness at the integration point
    real(kind=real64) :: N(2,2,fbem_n_nodes(etype))                     ! N matrix (matrix of shape functions)
    real(kind=real64) :: dNdx(fbem_n_nodes(etype),2)                    ! dN/dx
    integer           :: ki, kis, kie, kj, kjs, kje, kdj, kdk           ! Counters and nodal DOF limits
    ! Initialization
    dMdx=0.d0
    n_nodes=fbem_n_nodes(etype)
    ! Switch between types of elements
    select case (fbem_n_edges(etype))
      ! TRIANGULAR ELEMENTS
      case (3)
        ! Loop through integration points
        ! Transform to order for Wandzura rules
        rule=2*gln-1
        do kxit=1,wantri_n(rule)
          ! xi_1, xi_2 and w_t
          xi(1)=wantri_xi1(kxit,rule)
          xi(2)=wantri_xi2(kxit,rule)
          wt=wantri_w(kxit,rule)
          ! Shape functions and their first derivatives
#         define delta 0.0d0
#         include <phi_and_dphidxik_2d.rc>
#         undef delta
          ! Calculate jacobian matrix and determinant
          J=0.d0
          do i=1,n_nodes
            J(1,:)=J(1,:)+dphidxi1(i)*x_nodes(:,i)
            J(2,:)=J(2,:)+dphidxi2(i)*x_nodes(:,i)
          end do
          call fbem_invert_2x2_matrix(J,invJ,detJ)
          ! dN/dx
          do i=1,n_nodes
            dNdx(i,1)=invJ(1,1)*dphidxi1(i)+invJ(1,2)*dphidxi2(i)
            dNdx(i,2)=invJ(2,1)*dphidxi1(i)+invJ(2,2)*dphidxi2(i)
          end do
          ! Build matrix N of shape functions
          N=0.d0
          do i=1,n_nodes
            N(1,1,i)=phi(i)
            N(2,2,i)=phi(i)
          end do
          ! Thickness at the integration point
          t=dot_product(phi,t_nodes)
          ! thickness * det(J) * weights
          jw=detJ*t*wt
          ! Integrate dM/dx_j^k
          do kdj=1,2
            do kdk=1,n_nodes
              ! Build the element mass matrix
              do ki=1,n_nodes
                kis=(ki-1)*2+1
                kie=kis+1
                do kj=ki,n_nodes
                  kjs=(kj-1)*2+1
                  kje=kjs+1
                  dMdx(kis:kie,kjs:kje,kdj,kdk)=dMdx(kis:kie,kjs:kje,kdj,kdk)+matmul(transpose(N(:,:,ki)),N(:,:,kj))*dNdx(kdk,kdj)*jw
                end do
              end do
            end do
          end do
        end do ! Loop through integration points
      ! QUADRILATERAL ELEMENTS
      case (4)
        ! Loop through integration points
        ! Transform to order for Wandzura rules
        do kxi1=1,gln
          ! xi_1 and w_1
          xi(1)=gl11_xi(kxi1,gln)
          w(1)=gl11_w(kxi1,gln)
          do kxi2=1,gln
            ! xi_2 and w_2
            xi(2)=gl11_xi(kxi2,gln)
            w(2)=gl11_w(kxi2,gln)
            ! Shape functions and their first derivatives
#           define delta 0.0d0
#           include <phi_and_dphidxik_2d.rc>
#           undef delta
            ! Calculate jacobian matrix and determinant
            J=0.d0
            do i=1,n_nodes
              J(1,:)=J(1,:)+dphidxi1(i)*x_nodes(:,i)
              J(2,:)=J(2,:)+dphidxi2(i)*x_nodes(:,i)
            end do
            call fbem_invert_2x2_matrix(J,invJ,detJ)
            ! dN/dx
            do i=1,n_nodes
              dNdx(i,1)=invJ(1,1)*dphidxi1(i)+invJ(1,2)*dphidxi2(i)
              dNdx(i,2)=invJ(2,1)*dphidxi1(i)+invJ(2,2)*dphidxi2(i)
            end do
            ! Build matrix N of shape functions
            N=0.d0
            do i=1,n_nodes
              N(1,1,i)=phi(i)
              N(2,2,i)=phi(i)
            end do
            ! Thickness at the integration point
            t=dot_product(phi,t_nodes)
            ! thickness * det(J) * weights
            jw=detJ*t*w(1)*w(2)
            ! Integrate dM/dx_j^k
            do kdj=1,2
              do kdk=1,n_nodes
                ! Build the element mass matrix
                do ki=1,n_nodes
                  kis=(ki-1)*2+1
                  kie=kis+1
                  do kj=ki,n_nodes
                    kjs=(kj-1)*2+1
                    kje=kjs+1
                    dMdx(kis:kie,kjs:kje,kdj,kdk)=dMdx(kis:kie,kjs:kje,kdj,kdk)+matmul(transpose(N(:,:,ki)),N(:,:,kj))*dNdx(kdk,kdj)*jw
                  end do
                end do
              end do
            end do
          end do ! Loop through xi_2 coordinate
        end do ! Loop through xi_1 coordinate
      ! Unknown
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={tri3,tri6,quad4,quad8,quad9}')
    end select
    ! Copy symmetric part
    do kdj=1,2
      do kdk=1,n_nodes
        do ki=1,2*n_nodes
          do kj=ki,2*n_nodes
            dMdx(kj,ki,kdj,kdk)=dMdx(ki,kj,kdj,kdk)
          end do
        end do
      end do
    end do
    dMdx=rho*dMdx
  end subroutine fbem_fem_ela2d_solid_dMdx

  !! Calculate the mass matrix sensitivity of nodal coordinates dM/dx_k^i for 2D static analysis of isotropic solids. It is
  !! calculated using Finite Differences (4-points) with a deltax=1.e-4.
  subroutine fbem_fem_ela2d_solid_dMdx_FD(etype,x_nodes,t_nodes,rho,gln,dMdx)
    implicit none
    ! I/O
    integer           :: etype                                                                   !! Type of element: tri3, tri6, quad4, quad8, quad9.
    real(kind=real64) :: x_nodes(2,fbem_n_nodes(etype))                                          !! Coordinates of the nodes.
    real(kind=real64) :: t_nodes(fbem_n_nodes(etype))                                            !! Thickness at each node (t=1. for 2D plane strain problem)
    real(kind=real64) :: rho                                                                     !! Density
    integer           :: gln                                                                     !! Number of Gauss-Legendre integration points
    real(kind=real64) :: dMdx(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype),2,fbem_n_nodes(etype)) !! Mass matrix sensitivity dM/dx_k^i
    ! Local
    real(kind=real64) :: cl
    real(kind=real64) :: deltax
    real(kind=real64) :: x_nodes_FD(2,fbem_n_nodes(etype))
    real(kind=real64) :: M1(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    real(kind=real64) :: M2(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    real(kind=real64) :: M3(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    real(kind=real64) :: M4(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    integer           :: i, j
    ! Characteristic length
    cl=fbem_characteristic_length(2,etype,x_nodes,1.d-6)
    deltax=1.d-4*cl
    do i=1,fbem_n_nodes(etype)
      do j=1,2
        x_nodes_FD=x_nodes
        x_nodes_FD(j,i)=x_nodes_FD(j,i)+2.d0*deltax
        call fbem_fem_ela2d_solid_M(etype,x_nodes_FD,t_nodes,rho,gln,M1)
        x_nodes_FD=x_nodes
        x_nodes_FD(j,i)=x_nodes_FD(j,i)+deltax
        call fbem_fem_ela2d_solid_M(etype,x_nodes_FD,t_nodes,rho,gln,M2)
        x_nodes_FD=x_nodes
        x_nodes_FD(j,i)=x_nodes_FD(j,i)-deltax
        call fbem_fem_ela2d_solid_M(etype,x_nodes_FD,t_nodes,rho,gln,M3)
        x_nodes_FD=x_nodes
        x_nodes_FD(j,i)=x_nodes_FD(j,i)-2.d0*deltax
        call fbem_fem_ela2d_solid_M(etype,x_nodes_FD,t_nodes,rho,gln,M4)
        dMdx(:,:,j,i)=(M4-M1+8.d0*(M2-M3))/(12.d0*deltax)
      end do
    end do
  end subroutine fbem_fem_ela2d_solid_dMdx_FD

  subroutine fbem_fem_ela2d_solid_dQdx_body(etype,x_nodes,t_nodes,gln,dQdx)
    implicit none
    ! I/O
    integer           :: etype                                                                   !! Type of element: tri3, tri6, quad4, quad8, quad9.
    real(kind=real64) :: x_nodes(2,fbem_n_nodes(etype))                                          !! Coordinates of the nodes.
    real(kind=real64) :: t_nodes(fbem_n_nodes(etype))                                            !! Thickness at each node (t=1. for 2D plane strain problem)
    integer           :: gln                                                                     !! Number of Gauss-Legendre integration points
    real(kind=real64) :: dQdx(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype),2,fbem_n_nodes(etype)) !! Mass matrix sensitivity dQ^{body}/dx_k^i
    ! Local
    integer           :: n_nodes                                        ! Number of nodes
    integer           :: i                                              ! Counter
    integer           :: rule                                           ! Rule of Wandzura quadrature
    integer           :: kxi1, kxi2, kxit                               ! Integration points counters
    real(kind=real64) :: xi(2), w(2), wt                                ! Curvilinear coordinates and quadrature weights
    real(kind=real64) :: aux(10)                                        ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64) :: phi(fbem_n_nodes(etype))                       ! Shape functions
    real(kind=real64) :: dphidxi1(fbem_n_nodes(etype))                  ! Shape functions derivatives with respect to xi_1
    real(kind=real64) :: dphidxi2(fbem_n_nodes(etype))                  ! Shape functions derivatives with respect to xi_2
    real(kind=real64) :: J(2,2), invJ(2,2), detJ                        ! Jacobian matrix and its determinant
    real(kind=real64) :: jw                                             ! det(J) * weights
    real(kind=real64) :: t                                              ! Thickness at the integration point
    real(kind=real64) :: N(2,2,fbem_n_nodes(etype))                     ! N matrix (matrix of shape functions)
    real(kind=real64) :: dNdx(fbem_n_nodes(etype),2)                    ! dN/dx
    integer           :: ki, kis, kie, kj, kjs, kje, kdj, kdk           ! Counters and nodal DOF limits
    ! Initialization
    dQdx=0.d0
    n_nodes=fbem_n_nodes(etype)
    ! Switch between types of elements
    select case (fbem_n_edges(etype))
      ! TRIANGULAR ELEMENTS
      case (3)
        ! Loop through integration points
        ! Transform to order for Wandzura rules
        rule=2*gln-1
        do kxit=1,wantri_n(rule)
          ! xi_1, xi_2 and w_t
          xi(1)=wantri_xi1(kxit,rule)
          xi(2)=wantri_xi2(kxit,rule)
          wt=wantri_w(kxit,rule)
          ! Shape functions and their first derivatives
#         define delta 0.0d0
#         include <phi_and_dphidxik_2d.rc>
#         undef delta
          ! Calculate jacobian matrix and determinant
          J=0.d0
          do i=1,n_nodes
            J(1,:)=J(1,:)+dphidxi1(i)*x_nodes(:,i)
            J(2,:)=J(2,:)+dphidxi2(i)*x_nodes(:,i)
          end do
          call fbem_invert_2x2_matrix(J,invJ,detJ)
          ! dN/dx
          do i=1,n_nodes
            dNdx(i,1)=invJ(1,1)*dphidxi1(i)+invJ(1,2)*dphidxi2(i)
            dNdx(i,2)=invJ(2,1)*dphidxi1(i)+invJ(2,2)*dphidxi2(i)
          end do
          ! Build matrix N of shape functions
          N=0.d0
          do i=1,n_nodes
            N(1,1,i)=phi(i)
            N(2,2,i)=phi(i)
          end do
          ! Thickness at the integration point
          t=dot_product(phi,t_nodes)
          ! thickness * det(J) * weights
          jw=detJ*t*wt
          ! Integrate dM/dx_j^k
          do kdj=1,2
            do kdk=1,n_nodes
              ! Build the element mass matrix
              do ki=1,n_nodes
                kis=(ki-1)*2+1
                kie=kis+1
                do kj=ki,n_nodes
                  kjs=(kj-1)*2+1
                  kje=kjs+1
                  dQdx(kis:kie,kjs:kje,kdj,kdk)=dQdx(kis:kie,kjs:kje,kdj,kdk)+matmul(transpose(N(:,:,ki)),N(:,:,kj))*dNdx(kdk,kdj)*jw
                end do
              end do
            end do
          end do
        end do ! Loop through integration points
      ! QUADRILATERAL ELEMENTS
      case (4)
        ! Loop through integration points
        ! Transform to order for Wandzura rules
        do kxi1=1,gln
          ! xi_1 and w_1
          xi(1)=gl11_xi(kxi1,gln)
          w(1)=gl11_w(kxi1,gln)
          do kxi2=1,gln
            ! xi_2 and w_2
            xi(2)=gl11_xi(kxi2,gln)
            w(2)=gl11_w(kxi2,gln)
            ! Shape functions and their first derivatives
#           define delta 0.0d0
#           include <phi_and_dphidxik_2d.rc>
#           undef delta
            ! Calculate jacobian matrix and determinant
            J=0.d0
            do i=1,n_nodes
              J(1,:)=J(1,:)+dphidxi1(i)*x_nodes(:,i)
              J(2,:)=J(2,:)+dphidxi2(i)*x_nodes(:,i)
            end do
            call fbem_invert_2x2_matrix(J,invJ,detJ)
            ! dN/dx
            do i=1,n_nodes
              dNdx(i,1)=invJ(1,1)*dphidxi1(i)+invJ(1,2)*dphidxi2(i)
              dNdx(i,2)=invJ(2,1)*dphidxi1(i)+invJ(2,2)*dphidxi2(i)
            end do
            ! Build matrix N of shape functions
            N=0.d0
            do i=1,n_nodes
              N(1,1,i)=phi(i)
              N(2,2,i)=phi(i)
            end do
            ! Thickness at the integration point
            t=dot_product(phi,t_nodes)
            ! thickness * det(J) * weights
            jw=detJ*t*w(1)*w(2)
            ! Integrate dM/dx_j^k
            do kdj=1,2
              do kdk=1,n_nodes
                ! Build the element mass matrix
                do ki=1,n_nodes
                  kis=(ki-1)*2+1
                  kie=kis+1
                  do kj=ki,n_nodes
                    kjs=(kj-1)*2+1
                    kje=kjs+1
                    dQdx(kis:kie,kjs:kje,kdj,kdk)=dQdx(kis:kie,kjs:kje,kdj,kdk)+matmul(transpose(N(:,:,ki)),N(:,:,kj))*dNdx(kdk,kdj)*jw
                  end do
                end do
              end do
            end do
          end do ! Loop through xi_2 coordinate
        end do ! Loop through xi_1 coordinate
      ! Unknown
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={tri3,tri6,quad4,quad8,quad9}')
    end select
    ! Copy symmetric part
    do kdj=1,2
      do kdk=1,n_nodes
        do ki=1,2*n_nodes
          do kj=ki,2*n_nodes
            dQdx(kj,ki,kdj,kdk)=dQdx(ki,kj,kdj,kdk)
          end do
        end do
      end do
    end do
  end subroutine fbem_fem_ela2d_solid_dQdx_body

  !! Calculate the body load matrix dQdx for 2D solids.
  subroutine fbem_fem_ela2d_solid_dQdx_body_FD(etype,x_nodes,t_nodes,gln,dQdx)
    implicit none
    ! I/O
    integer           :: etype                                                                   !! Type of element: tri3, tri6, quad4, quad8, quad9.
    real(kind=real64) :: x_nodes(2,fbem_n_nodes(etype))                                          !! Coordinates of the nodes.
    real(kind=real64) :: t_nodes(fbem_n_nodes(etype))                                            !! Thickness at each node (t=1. for 2D plane strain problem)
    integer           :: gln                                                                     !! Number of Gauss-Legendre integration points
    real(kind=real64) :: dQdx(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype),2,fbem_n_nodes(etype)) !! Body load matrix
    ! Local
    real(kind=real64) :: cl
    real(kind=real64) :: deltax
    real(kind=real64) :: x_nodes_FD(2,fbem_n_nodes(etype))
    real(kind=real64) :: Q1(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    real(kind=real64) :: Q2(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    real(kind=real64) :: Q3(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    real(kind=real64) :: Q4(2*fbem_n_nodes(etype),2*fbem_n_nodes(etype))
    integer           :: i, j
    ! Characteristic length
    cl=fbem_characteristic_length(2,etype,x_nodes,1.d-6)
    deltax=1.d-4*cl
    do i=1,fbem_n_nodes(etype)
      do j=1,2
        x_nodes_FD=x_nodes
        x_nodes_FD(j,i)=x_nodes_FD(j,i)+2.d0*deltax
        call fbem_fem_ela2d_solid_Q_body(etype,x_nodes_FD,t_nodes,gln,Q1)
        x_nodes_FD=x_nodes
        x_nodes_FD(j,i)=x_nodes_FD(j,i)+deltax
        call fbem_fem_ela2d_solid_Q_body(etype,x_nodes_FD,t_nodes,gln,Q2)
        x_nodes_FD=x_nodes
        x_nodes_FD(j,i)=x_nodes_FD(j,i)-deltax
        call fbem_fem_ela2d_solid_Q_body(etype,x_nodes_FD,t_nodes,gln,Q3)
        x_nodes_FD=x_nodes
        x_nodes_FD(j,i)=x_nodes_FD(j,i)-2.d0*deltax
        call fbem_fem_ela2d_solid_Q_body(etype,x_nodes_FD,t_nodes,gln,Q4)
        dQdx(:,:,j,i)=(Q4-Q1+8.d0*(Q2-Q3))/(12.d0*deltax)
      end do
    end do
  end subroutine fbem_fem_ela2d_solid_dQdx_body_FD

  subroutine fbem_fem_ela2d_solid_dQdx_edge(etype,x_nodes,t_nodes,kedge,dtype,ngp,dQdx)
    implicit none
    ! I/O
    integer           :: etype                                                                   !! Type of element: tri3, tri6, quad4, quad8, quad9.
    real(kind=real64) :: x_nodes(2,fbem_n_nodes(etype))                                          !! Position vectors of the mid-plane nodes.
    real(kind=real64) :: t_nodes(fbem_n_nodes(etype))                                            !! Thickness of each mid-node in the v_3 direction.
    integer           :: kedge                                                                   !! Selected edge
    integer           :: dtype                                                                   !! Type of edge load element: point1, line2, line3.
    integer           :: ngp                                                                     !! Number of Gauss points for edge coordinate integration
    real(kind=real64) :: dQdx(2*fbem_n_nodes(etype),2*fbem_n_nodes(dtype),2,fbem_n_nodes(etype)) !! Load matrix Q for a global coordinates distributed load over an edge
    ! Local
    integer           :: n_nodes                              ! Number of nodes of the midplane element
    integer           :: n_dnodes                             ! Number of nodes of the load element
    integer           :: k_node                               ! Counter of mid-nodes
    integer           :: kp                                   ! Integration points counter
    integer           :: ki, kdj, kdk                         ! Counter
    real(kind=real64) :: xip                                  ! xip coordinate
    real(kind=real64) :: w                                    ! Weight
    real(kind=real64) :: xi_e(2,2)                            ! xi coordinates of the edge
    real(kind=real64) :: xi_nodes(2,fbem_n_nodes(etype))      ! xi1,xi2 coordinates of the nodes of the element
    real(kind=real64) :: dxidxip(2)                           ! xip->xi1,xi2 tangent
    real(kind=real64) :: xi(2)                                ! xi coordinates
    real(kind=real64) :: aux(10)                              ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64) :: phi(fbem_n_nodes(etype))             ! phi shape functions
    real(kind=real64) :: dphidxi1(fbem_n_nodes(etype))        ! phi shape functions derivatives with respect to xi_1
    real(kind=real64) :: dphidxi2(fbem_n_nodes(etype))        ! phi shape functions derivatives with respect to xi_2
    real(kind=real64) :: dxdxi1(2), dxdxi2(2)                 ! dx/dxik
    real(kind=real64) :: dxdxip(2)                            ! dx/dxip
    real(kind=real64) :: t                                    ! Thickness at the integration point
    real(kind=real64) :: jg                                   ! Area jacobian
    real(kind=real64) :: jw                                   ! d |J^{be}|/dx_j^k * thickness * weight
    real(kind=real64) :: dJbedxjk                             ! d |J^{be}|/dx_j^k
    real(kind=real64) :: Ne(2,2*fbem_n_nodes(etype))          ! Element shape functions matrix
    real(kind=real64) :: Nedge(2,2*fbem_n_nodes(dtype))       ! Edge load shape functions matrix
    ! Check
    if (kedge.gt.fbem_n_edges(etype)) then
      call fbem_error_message(error_unit,0,'fbem_fem_ela2d_solid_Q_edge',kedge,'invalid edge index')
    end if
    ! Initialize load matrix
    dQdx=0.d0
    ! Number of nodes of the midplane element
    n_nodes=fbem_n_nodes(etype)
    ! Number of nodes of the load element
    n_dnodes=fbem_n_nodes(dtype)
    ! xi1,xi2 coordinates of kedge
    xi_nodes=fbem_xi1xi2(etype,0.d0)
    xi_e(:,1)=xi_nodes(:,fbem_edge_node(1,kedge,etype))
    xi_e(:,2)=xi_nodes(:,fbem_edge_node(2,kedge,etype))
    ! Loop through xip
    do kp=1,gl11_n(ngp)
      ! XIP COORDINATE AND WEIGHT
      xip=gl11_xi(kp,ngp)
      w=gl11_w(kp,ngp)
      ! XIP->XI(1),XI(2) TRANSFORMATION
      xi=0.5d0*(1.d0-xip)*xi_e(:,1)+0.5d0*(1.d0+xip)*xi_e(:,2)
      dxidxip=0.5d0*(xi_e(:,2)-xi_e(:,1))
      ! Solid element shape functions and its derivatives
#     define delta 0.0d0
#     include <phi_and_dphidxik_2d.rc>
#     undef delta
      ! XIP -> XI -> X TRANSFORMATIONS
      dxdxi1=0.d0
      dxdxi2=0.d0
      do k_node=1,n_nodes
        dxdxi1=dxdxi1+dphidxi1(k_node)*x_nodes(:,k_node)
        dxdxi2=dxdxi2+dphidxi2(k_node)*x_nodes(:,k_node)
      end do
      dxdxip=dxdxi1*dxidxip(1)+dxdxi2*dxidxip(2)
      jg=sqrt(dot_product(dxdxip,dxdxip))
      ! Thickness at the integration point
      t=dot_product(phi,t_nodes)
      ! Build element shape function matrix Ne
      Ne=0.d0
      do k_node=1,n_nodes
        ki=2*(k_node-1)
        Ne(1,ki+1)=phi(k_node)
        Ne(2,ki+2)=phi(k_node)
      end do
      ! Edge load shape functions
      phi=0.d0
#     define etype dtype
#     define delta 0.0d0
#     define xi xip
#     include <phi_1d.rc>
#     undef etype
#     undef delta
#     undef xi
      ! Build load shape function matrix Nedge
      Nedge=0.d0
      do k_node=1,n_dnodes
        ki=2*(k_node-1)
        Nedge(1,ki+1)=phi(k_node)
        Nedge(2,ki+2)=phi(k_node)
      end do
      ! Add
      do kdj=1,2
        do kdk=1,n_nodes
          dJbedxjk=dxdxip(kdj)/jg*(dphidxi1(kdk)*dxidxip(1)+dphidxi2(kdk)*dxidxip(2))
          jw=dJbedxjk*t*w
          dQdx(:,:,kdj,kdk)=dQdx(:,:,kdj,kdk)+matmul(transpose(Ne),Nedge)*jw
        end do
      end do

    end do
  end subroutine fbem_fem_ela2d_solid_dQdx_edge

  !! Calculate the edge load matrix sensitivity of nodal coordinates dQ^{edge}/dx_k^i for 2D solid FE. It is
  !! calculated using Finite Differences (4-points) with a deltax=1.e-4.
  subroutine fbem_fem_ela2d_solid_dQdx_edge_FD(etype,x_nodes,t_nodes,kedge,dtype,ngp,dQdx)
    implicit none
    ! I/O
    integer           :: etype                                                                   !! Type of element: tri3, tri6, quad4, quad8, quad9.
    real(kind=real64) :: x_nodes(2,fbem_n_nodes(etype))                                          !! Position vectors of the mid-plane nodes.
    real(kind=real64) :: t_nodes(fbem_n_nodes(etype))                                            !! Thickness of each mid-node in the v_3 direction.
    integer           :: kedge                                                                   !! Selected edge
    integer           :: dtype                                                                   !! Type of edge load element: point1, line2, line3.
    integer           :: ngp                                                                     !! Number of Gauss points for edge coordinate integration
    real(kind=real64) :: dQdx(2*fbem_n_nodes(etype),2*fbem_n_nodes(dtype),2,fbem_n_nodes(etype)) !! Load matrix Q for a global coordinates distributed load over an edge
    ! Local
    real(kind=real64) :: cl
    real(kind=real64) :: deltax
    real(kind=real64) :: x_nodes_FD(2,fbem_n_nodes(etype))
    real(kind=real64) :: Q1(2*fbem_n_nodes(etype),2*fbem_n_nodes(dtype))
    real(kind=real64) :: Q2(2*fbem_n_nodes(etype),2*fbem_n_nodes(dtype))
    real(kind=real64) :: Q3(2*fbem_n_nodes(etype),2*fbem_n_nodes(dtype))
    real(kind=real64) :: Q4(2*fbem_n_nodes(etype),2*fbem_n_nodes(dtype))
    integer           :: i, j
    ! Characteristic length
    cl=fbem_characteristic_length(2,etype,x_nodes,1.d-6)
    deltax=1.d-4*cl
    do i=1,fbem_n_nodes(etype)
      do j=1,2
        x_nodes_FD=x_nodes
        x_nodes_FD(j,i)=x_nodes_FD(j,i)+2.d0*deltax
        call fbem_fem_ela2d_solid_Q_edge(etype,x_nodes_FD,t_nodes,kedge,dtype,ngp,Q1)
        x_nodes_FD=x_nodes
        x_nodes_FD(j,i)=x_nodes_FD(j,i)+deltax
        call fbem_fem_ela2d_solid_Q_edge(etype,x_nodes_FD,t_nodes,kedge,dtype,ngp,Q2)
        x_nodes_FD=x_nodes
        x_nodes_FD(j,i)=x_nodes_FD(j,i)-deltax
        call fbem_fem_ela2d_solid_Q_edge(etype,x_nodes_FD,t_nodes,kedge,dtype,ngp,Q3)
        x_nodes_FD=x_nodes
        x_nodes_FD(j,i)=x_nodes_FD(j,i)-2.d0*deltax
        call fbem_fem_ela2d_solid_Q_edge(etype,x_nodes_FD,t_nodes,kedge,dtype,ngp,Q4)
        dQdx(:,:,j,i)=(Q4-Q1+8.d0*(Q2-Q3))/(12.d0*deltax)
      end do
    end do
  end subroutine fbem_fem_ela2d_solid_dQdx_edge_FD

end module fbem_fem_solids
