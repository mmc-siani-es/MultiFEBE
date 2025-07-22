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
!! <b> This module implements the calculation of element-wise integrals of Boundary Integral Equations of the 3D Helmholtz
!! problem.</b>
module fbem_bem_harpot3d

  ! Fortran 2003 standard
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_quad_rules
  use fbem_telles_transformation
  use fbem_polar_transformation
  use fbem_quasisingular_integration
  use fbem_geometry
  use fbem_bem_general
  use fbem_bem_stapot3d

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all are private
  private

  ! ================================================================================================================================
  ! INITIAL SETUP
  public :: fbem_bem_harpot3d_parameters
  public :: fbem_bem_harpot3d_calculate_parameters
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE)
  ! Fundamental solution
  public :: fbem_bem_harpot3d_sbie_p
  public :: fbem_bem_harpot3d_sbie_q
  ! BOUNDARY ELEMENTS
  ! Exterior integration
  public :: fbem_bem_harpot3d_sbie_ext_pre
  public :: fbem_bem_harpot3d_sbie_ext_st
  public :: fbem_bem_harpot3d_sbie_ext_adp
  ! Interior integration
  public :: fbem_bem_harpot3d_sbie_int
  ! Automatic integration
  public :: fbem_bem_harpot3d_sbie_auto
  ! BODY LOADS
  public :: fbem_bem_harpot3d_sbie_bl_auto
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (HBIE)
  public :: fbem_bem_harpot3d_hbie_d
  public :: fbem_bem_harpot3d_hbie_s
  ! BOUNDARY ELEMENTS
  ! Exterior integration
  public :: fbem_bem_harpot3d_hbie_ext_pre
  public :: fbem_bem_harpot3d_hbie_ext_st
  public :: fbem_bem_harpot3d_hbie_ext_adp
  ! Interior integration
  public :: fbem_bem_harpot3d_hbie_int
  ! Automatic integration
  public :: fbem_bem_harpot3d_hbie_auto
  ! BODY LOADS
  public :: fbem_bem_harpot3d_hbie_bl_auto
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! BOTH SBIE AND HBIE SIMULTANEOUSLY
  ! BOUNDARY ELEMENTS
  ! Exterior integration
  public :: fbem_bem_harpot3d_shbie_ext_pre
  ! Automatic integration
  public :: fbem_bem_harpot3d_shbie_auto
  ! ================================================================================================================================

  ! ================================================================================================================================
  !! Region parameters (material and frequency dependant)
  type fbem_bem_harpot3d_parameters
    ! Region parameters
    real(kind=real64)    :: rho
    complex(kind=real64) :: c
    ! Frequency
    real(kind=real64)    :: omega
    ! Wavenumber
    complex(kind=real64) :: k
    ! Coefficients of the components of the fundamental solution
    complex(kind=real64) :: P(1)
    complex(kind=real64) :: Q(2)
    complex(kind=real64) :: S1(5), S2(3)
  end type fbem_bem_harpot3d_parameters
  ! ================================================================================================================================

contains

  ! ================================================================================================================================
  ! INITIAL SETUP

  subroutine fbem_bem_harpot3d_calculate_parameters(rho,c,omega,pars)
    implicit none
    ! I/O
    real(kind=real64)                  :: rho
    complex(kind=real64)               :: c
    real(kind=real64)                  :: omega
    type(fbem_bem_harpot3d_parameters) :: pars
    ! Local
    complex(kind=real64)               :: k
    !
    ! Frequency-dependant parameters
    !
    k=omega/c
    !
    ! Save to structure
    !
    pars%omega=omega
    pars%rho=rho
    pars%c=c
    pars%k=k
    !
    ! Coefficients of equations of fundamental solution's components
    !
    ! P
    pars%P(1)=-c_im*k
    ! Q
    pars%Q(1)= 0.5d0*k**2
    pars%Q(2)= c_im*k
    ! S1
    pars%S1(1)= 3.d0
    pars%S1(2)= 0.5d0*k**2
    pars%S1(3)=-k**2
    pars%S1(4)= 3.d0*c_im*k
    pars%S1(5)= 3.d0
    ! S2
    pars%S2(1)= 0.5d0*k**2
    pars%S2(2)=-0.333333333333333333d0*c_im*k**3
    pars%S2(3)= c_im*k

    ! Incluir c_1_4pi en las constantes y quitarlas de las rutinas


  end subroutine fbem_bem_harpot3d_calculate_parameters

!  subroutine fbem_bem_harpot3d_test_krlim()
!    implicit none
!    ! I/O
!    real(kind=real64)                  :: rho
!    real(kind=real64)                  :: c
!    real(kind=real64)                  :: omega
!    type(fbem_bem_harpot3d_parameters) :: p
!    integer                            :: ko, no, kr, nr
!    real(kind=real64)                  :: o_max, o_min, r_max, r_min
!    real(kind=real64)                  :: deltao, deltar
!    real(kind=real64)                  :: r, d1r, d1r2, d1r3, k
!    complex(kind=real64)               :: z, expz, E0, E1, E2, E3, E4, FSP, FSQ, FSS1, FSS2
!    ! Properties
!    rho=1000.
!    c  =1400.
!    ! omega setup
!    o_min=1.d-12
!    o_max=1.d12
!    no=100
!    ! r setup
!    r_min=1.d-12
!    r_max=1.d12
!    nr=100
!    ! Delta of omega and r
!    deltao=(log10(o_max)-log10(o_min))/dble(no-1)
!    deltar=(log10(r_max)-log10(r_min))/dble(nr-1)
!    ! Loop through omega
!    do ko=1,no
!      omega=10**(dlog10(o_min)+deltao*dble(ko-1))
!      call fbem_bem_harpot3d_calculate_parameters(rho,c,omega,p)
!      k=p%k
!      ! Loop through r
!      do kr=1,nr
!        r=10**(dlog10(r_min)+deltar*dble(kr-1))
!        d1r=1.d0/r
!        d1r2=1.d0/r**2
!        d1r3=1.d0/r**3
!        z=-c_im*k*r
!        expz=exp(z)
!        call fbem_decomposed_zexp(z,E0,E1,E2,E3,E4)
!        ! Full way
!        FSP=d1r*expz
!        FSQ=(d1r2+p%Q(2)*d1r)*expz
!        FSS1=(p%S1(3)*d1r+p%S1(4)*d1r2+p%S1(5)*d1r3)*expz
!        FSS2=(p%S2(3)*d1r2+d1r3)*expz
!        write(10,'(i11,e25.16,i11,10e25.16)') ko, omega, kr, r, omega*r, dreal(FSP), dreal(FSQ), dreal(FSS1), dreal(FSS2), dimag(FSP), dimag(FSQ), dimag(FSS1), dimag(FSS2)
!        ! Decomposed way
!        FSP=d1r+p%P(1)+d1r*E2
!        FSQ=d1r2+p%Q(1)+p%Q(2)*d1r*E2+d1r2*E3
!        FSS1=p%S1(1)*d1r3+p%S1(2)*d1r+p%S1(3)*d1r*E2+p%S1(4)*d1r2*E3+p%S1(5)*d1r3*E4
!        FSS2=d1r3+p%S2(1)*d1r+p%S2(2)+p%S2(3)*d1r2*E3+d1r3*E4
!        write(11,'(i11,e25.16,i11,10e25.16)') ko, omega, kr, r, omega*r, dreal(FSP), dreal(FSQ), dreal(FSS1), dreal(FSS2), dimag(FSP), dimag(FSQ), dimag(FSS1), dimag(FSS2)
!      end do
!    end do
!  end subroutine

  ! ================================================================================================================================

  ! ================================================================================================================================
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE)

  !! Fundamental solution p*
  subroutine fbem_bem_harpot3d_sbie_p(x,x_i,p,po)
    implicit none
    ! I/O
    real(kind=real64)                  :: x(3)    !! Observation point
    real(kind=real64)                  :: x_i(3)  !! Collocation point
    type(fbem_bem_harpot3d_parameters) :: p       !! Parameters of the region
    complex(kind=real64)               :: po      !! p*
    ! Local
    real(kind=real64)    :: rv(3), r, d1r1, d1r2
    complex(kind=real64) :: z(1), EnR(0:6,1)
    rv=x-x_i
    r=sqrt(dot_product(rv,rv))
    d1r1=1.d0/r
    d1r2=d1r1**2
    z(1)=-c_im*p%k*r
    call fbem_zexp_decomposed(1,z,EnR)
    EnR(2,:)=EnR(2,:)*d1r1
    EnR(3,:)=EnR(3,:)*d1r2
    po=c_1_4pi*(d1r1+p%P(1)+EnR(2,1))
  end subroutine fbem_bem_harpot3d_sbie_p

  !! Fundamental solution q*
  subroutine fbem_bem_harpot3d_sbie_q(x,n,x_i,p,qo)
    implicit none
    ! I/O
    real(kind=real64)                  :: x(3)    !! Observation point
    real(kind=real64)                  :: n(3)    !! Observation point normal
    real(kind=real64)                  :: x_i(3)  !! Collocation point
    type(fbem_bem_harpot3d_parameters) :: p       !! Parameters of the region
    complex(kind=real64)               :: qo      !! q*
    ! Local
    real(kind=real64)    :: rv(3), r, d1r1, d1r2, drdx(3), drdn
    complex(kind=real64) :: z(1), EnR(0:6,1)
    rv=x-x_i
    r=sqrt(dot_product(rv,rv))
    d1r1=1.d0/r
    d1r2=d1r1**2
    drdx=rv*d1r1
    drdn=dot_product(drdx,n)
    z(1)=-c_im*p%k*r
    call fbem_zexp_decomposed(1,z,EnR)
    EnR(2,:)=EnR(2,:)*d1r1
    EnR(3,:)=EnR(3,:)*d1r2
    qo=-c_1_4pi*(d1r2+p%Q(1)+p%Q(2)*EnR(2,1)+EnR(3,1))*drdn
  end subroutine fbem_bem_harpot3d_sbie_q


  subroutine fbem_bem_harpot3d_sbie_ext_pre(ps,e,reverse,x_i,p,h,g)
    implicit none
    ! I/O
    integer                            :: ps            !! Selected precalculated dataset
    type(fbem_bem_element)             :: e             !! Element
    logical                            :: reverse       !! Reverse normal vector
    real(kind=real64)                  :: x_i(3)        !! Collocation point position vector
    type(fbem_bem_harpot3d_parameters) :: p             !! Parameters of the region
    complex(kind=real64)               :: h(e%n_pnodes) !! h integration kernels matrix
    complex(kind=real64)               :: g(e%n_snodes) !! g integration kernels matrix
    ! Local
    integer              :: kip                ! Counter variable for integration points loop
    real(kind=real64)    :: x(3)               ! Position vector at integration point
    real(kind=real64)    :: n(3)               ! Unit normal vector at integration point
    real(kind=real64)    :: pphijw(e%n_pnodes) ! phi^p * jacobian * weight at integration point
    real(kind=real64)    :: sphijw(e%n_snodes) ! phi^s * jacobian * weight at integration point
    real(kind=real64)    :: rv(3)              ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, d1r1, d1r2      ! Distance vector module and its inverse
    real(kind=real64)    :: drdx(3)            ! Distance vector derivatives with respect to x_k
    real(kind=real64)    :: drdn               ! Partial derivative of r respect to unit normal
    complex(kind=real64) :: z(1)               ! Arguments z=-ikr
    complex(kind=real64) :: EnR(0:6,1)         ! Exponential function decomposition for each wavenumber
    complex(kind=real64) :: fs_P, fs_Q         ! Components of the fundamental solutions
    ! Initialize
    h=(0.d0,0.d0)
    g=(0.d0,0.d0)
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      x=e%ps_x(:,kip,ps)
      n=e%ps_n(:,kip,ps)
      pphijw=e%ps_pphijw(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      d1r1=1.d0/r
      d1r2=d1r1**2
      drdx=rv*d1r1
      drdn=dot_product(drdx,n)
      z(1)=-c_im*p%k*r
      call fbem_zexp_decomposed(1,z,EnR)
      EnR(2,:)=EnR(2,:)*d1r1
      EnR(3,:)=EnR(3,:)*d1r2
      fs_P=d1r1+p%P(1)+EnR(2,1)
      fs_Q=d1r2+p%Q(1)+p%Q(2)*EnR(2,1)+EnR(3,1)
      h=h+fs_Q*drdn*pphijw
      g=g+fs_P*sphijw
    end do
    ! Multiply by constants
    h=-h*c_1_4pi
    g= g*c_1_4pi
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_harpot3d_sbie_ext_pre

  subroutine fbem_bem_harpot3d_sbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,p,gln,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                                !! Integration element
    logical                            :: reverse                          !! Reverse normal vector
    real(kind=real64)                  :: xi_s(2,fbem_n_vertices(e%gtype)) !! Coordinates of the subdivision of the parent element
    real(kind=real64)                  :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)                  :: barxip(2)                        !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)                  :: barr                             !! Telles jacobian at barxip
    type(fbem_bem_harpot3d_parameters) :: p                                !! Parameters of the region
    integer                            :: gln                              !! 1D Gauss-Legendre number of integration points (<=32)
    complex(kind=real64)               :: h(e%n_pnodes)                    !! h kernel vector
    complex(kind=real64)               :: g(e%n_snodes)                    !! g kernel vector
    ! Local
    integer                      :: kphi                      ! Counter variable for shape functions loops
    integer                      :: k1                        ! Counter variable for reference coordinate xi_1
    integer                      :: k2                        ! Counter variable for reference coordinate xi_2
    real(kind=real64)            :: aux(10)                   ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)          ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi1(e%n_gnodes)     ! Geometrical shape functions derivatives values
    real(kind=real64)            :: dgphidxi2(e%n_gnodes)     ! Geometrical shape functions derivatives values
    real(kind=real64)            :: pphi(e%n_pnodes)          ! Functional shape functions values
    real(kind=real64)            :: sphi(e%n_snodes)          ! Functional shape functions values
    real(kind=real64)            :: gamma(2)                  ! Vector of gamma_1,gamma_2 coordinates
    real(kind=real64)            :: w(2)                      ! Weights of the integration rule
    real(kind=real64)            :: xip(2)                    ! Vector of xip_1,xip_2 coordinates
    real(kind=real64)            :: dxidxi1p(2), dxidxi2p(2)  ! xi derivatives with respect to xip
    real(kind=real64)            :: js                        ! Subdivision jacobian
    real(kind=real64)            :: xi(2)                     ! Vector of xi_1,xi_2 coordinates
    real(kind=real64)            :: xipp(2)                   ! Coordinate xipp used for quadrilateral-triangle transformation
    real(kind=real64)            :: barxipp(2)                ! Coordinate xipp of collocation point
    real(kind=real64)            :: jqt                       ! Jacobian of the quadrilateral-triangle transformation
    type(fbem_telles_parameters) :: telles_parameters(2)      ! Telles parameters for each coordinate
    real(kind=real64)            :: jt(2)                     ! Telles jacobian for each coordinate: xi_1->gamma_1 and xi_2->gamma_2
    real(kind=real64)            :: x(3)                      ! Position vector at xi_1,xi_2
    real(kind=real64)            :: T1(3), T2(3)              ! Tangent vectors at xi_1,xi_2
    real(kind=real64)            :: N(3)                      ! Normal vector at xi_1,xi_2
    real(kind=real64)            :: rv(3)                     ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)            :: r, d1r1, d1r2             ! Distance vector module and its inverse
    real(kind=real64)            :: drdx(3)                   ! Distance vector derivatives with respect to x_k
    real(kind=real64)            :: jg                        ! Geometric jacobian
    real(kind=real64)            :: jw                        ! Jacobians * weights
    real(kind=real64)            :: drdn                      ! Partial derivative of r respect to unit normal
    real(kind=real64)            :: pphijw(e%n_pnodes)        ! Functional shape functions * jw
    real(kind=real64)            :: sphijw(e%n_snodes)        ! Functional shape functions * jw
    complex(kind=real64)         :: z(1)                      ! Arguments z=-ikr
    complex(kind=real64)         :: EnR(0:6,1)                ! Exponential function decomposition for each wavenumber
    complex(kind=real64)         :: fs_P, fs_Q                ! Components of the fundamental solutions
    ! Initialization
    h=(0.d0,0.d0)
    g=(0.d0,0.d0)

    ! Numerical integration
    select case (fbem_n_vertices(e%gtype))
      ! QUADRILATERAL ELEMENTS
      case (4)
        ! Calculate Telles parameters for each direction
        telles_parameters(1)=fbem_telles11_calculate_parameters(barxip(1),barr)
        telles_parameters(2)=fbem_telles11_calculate_parameters(barxip(2),barr)
        ! Loop through INTEGRATION POINTS
        do k1=1,gl11_n(gln)
          ! GAMMA1 COORDINATE
          gamma(1)=gl11_xi(k1,gln)
          w(1)=gl11_w(k1,gln)
          ! GAMMA1->XIP1 TRANSFORMATION
          call fbem_telles_xi_and_jacobian(telles_parameters(1),gamma(1),xip(1),jt(1))
          do k2=1,gl11_n(gln)
            ! GAMMA2 COORDINATE
            gamma(2)=gl11_xi(k2,gln)
            w(2)=gl11_w(k2,gln)
            ! GAMMA2->XIP2 TRANSFORMATION
            call fbem_telles_xi_and_jacobian(telles_parameters(2),gamma(2),xip(2),jt(2))
            ! XIP->XI TRANSFORMATION
            ! Shape functions and its derivatives
#           define delta 0.d0
#           define xi xip
#           define phi gphi
#           define dphidxi1 dgphidxi1
#           define dphidxi2 dgphidxi2
#           include <phi_quad4.rc>
#           include <dphidxi1_quad4.rc>
#           include <dphidxi2_quad4.rc>
#           undef delta
#           undef xi
#           undef phi
#           undef dphidxi1
#           undef dphidxi2
            ! xi coordinates, and xi derivatives
            xi=0.d0
            dxidxi1p=0.d0
            dxidxi2p=0.d0
            do kphi=1,4
              xi=xi+gphi(kphi)*xi_s(:,kphi)
              dxidxi1p=dxidxi1p+dgphidxi1(kphi)*xi_s(:,kphi)
              dxidxi2p=dxidxi2p+dgphidxi2(kphi)*xi_s(:,kphi)
            end do
            ! xip->xi jacobian
            js=dxidxi1p(1)*dxidxi2p(2)-dxidxi1p(2)*dxidxi2p(1)
            ! XI->X TRANSFORMATION
            ! Geometrical shape functions and first derivatives at xi
#           define etype e%gtype
#           define delta 0.d0
#           define phi gphi
#           define dphidxi1 dgphidxi1
#           define dphidxi2 dgphidxi2
#           include <phi_and_dphidxik_2d.rc>
#           undef etype
#           undef delta
#           undef phi
#           undef dphidxi1
#           undef dphidxi2
            ! Components calculation of x, T1 and T2 at xi
            x=0.d0
            T1=0.d0
            T2=0.d0
            do kphi=1,e%n_gnodes
              x=x+gphi(kphi)*e%x(:,kphi)
              T1=T1+dgphidxi1(kphi)*e%x(:,kphi)
              T2=T2+dgphidxi2(kphi)*e%x(:,kphi)
            end do
            ! Normal vector as T1 x T2 at xi
            N(1)=T1(2)*T2(3)-T1(3)*T2(2)
            N(2)=T1(3)*T2(1)-T1(1)*T2(3)
            N(3)=T1(1)*T2(2)-T1(2)*T2(1)
            ! Geometric jacobian
            jg=sqrt(dot_product(N,N))
            ! Unit normal
            n=N/jg
            ! Distance vector
            rv=x-x_i
            ! Distance vector norm
            r=sqrt(dot_product(rv,rv))
            d1r1=1.d0/r
            d1r2=d1r1**2
            drdx=rv*d1r1
            drdn=dot_product(drdx,n)
            ! Jacobians * weights
            jw=jg*js*jt(1)*jt(2)*w(1)*w(2)
            ! FUNCTIONAL SHAPE FUNCTIONS
            ! Functional shape functions (primary variables) at xi
#           define etype e%ptype
#           define delta e%ptype_delta
#           define phi pphi
#           include <phi_2d.rc>
#           undef etype
#           undef delta
#           undef phi
            ! Functional shape functions (secondary variables) at xi
#           define etype e%stype
#           define delta e%stype_delta
#           define phi sphi
#           include <phi_2d.rc>
#           undef etype
#           undef delta
#           undef phi
            ! Functional shape functions * jacobians * weights
            pphijw=pphi*jw
            sphijw=sphi*jw
            ! COMPONENTS OF THE FUNDAMENTAL SOLUTION
            z(1)=-c_im*p%k*r
            call fbem_zexp_decomposed(1,z,EnR)
            EnR(2,:)=EnR(2,:)*d1r1
            EnR(3,:)=EnR(3,:)*d1r2
            fs_P=d1r1+p%P(1)+EnR(2,1)
            fs_Q=d1r2+p%Q(1)+p%Q(2)*EnR(2,1)+EnR(3,1)
            ! Add
            h=h+fs_Q*drdn*pphijw
            g=g+fs_P*sphijw
          end do
        end do
      ! TRIANGULAR ELEMENTS
      case (3)
        ! Telles transformation is applied to Gauss-Legendre*Gauss-Legendre quadrature before the quadrilateral->triangle
        ! transformation. Because barxi for triangles are given in area triangle coordinates, for Telles transformation
        ! they must be transformed to quadrilateral coordinates. A special treatment is needed when barxi_2 is near 1, because
        ! transformation diverges.
        if (barxip(2).gt.0.995d0) then
          barxipp(1)=0.5d0
          barxipp(2)=1.0d0
        else
          barxipp(1)=barxip(1)/(1.0d0-barxip(2))
          barxipp(2)=barxip(2)
        end if
        ! Calculate Telles parameters
        telles_parameters(1)=fbem_telles01_calculate_parameters(barxipp(1),barr)
        telles_parameters(2)=fbem_telles01_calculate_parameters(barxipp(2),barr)
        ! Loop through INTEGRATION POINTS
        do k1=1,gl01_n(gln)
          ! GAMMA1 COORDINATE
          gamma(1)=gl01_xi(k1,gln)
          w(1)=gl01_w(k1,gln)
          ! GAMMA1->XIPP1 TRANSFORMATION
          ! xipp_1 coordinate and jacobian from Telles transformation
          call fbem_telles_xi_and_jacobian(telles_parameters(1),gamma(1),xipp(1),jt(1))
          do k2=1,gl01_n(gln)
            ! GAMMA2 COORDINATE
            gamma(2)=gl01_xi(k2,gln)
            w(2)=gl01_w(k2,gln)
            ! GAMMA2->XIPP2 TRANSFORMATION
            call fbem_telles_xi_and_jacobian(telles_parameters(2),gamma(2),xipp(2),jt(2))
            ! XIPP->XIP TRANSFORMATION (QUAD->TRI)
            xip(1)=(1.d0-xipp(2))*xipp(1)
            xip(2)=xipp(2)
            jqt=1.d0-xipp(2)
            ! XIP->XI TRANSFORMATION
            ! Shape functions and its derivatives
#           define delta 0.d0
#           define xi xip
#           define phi gphi
#           define dphidxi1 dgphidxi1
#           define dphidxi2 dgphidxi2
#           include <phi_tri3.rc>
#           include <dphidxi1_tri3.rc>
#           include <dphidxi2_tri3.rc>
#           undef delta
#           undef xi
#           undef phi
#           undef dphidxi1
#           undef dphidxi2
            ! xi coordinates, and xi derivatives
            xi=0.d0
            dxidxi1p=0.d0
            dxidxi2p=0.d0
            do kphi=1,3
              xi=xi+gphi(kphi)*xi_s(:,kphi)
              dxidxi1p=dxidxi1p+dgphidxi1(kphi)*xi_s(:,kphi)
              dxidxi2p=dxidxi2p+dgphidxi2(kphi)*xi_s(:,kphi)
            end do
            ! xip->xi jacobian
            js=dxidxi1p(1)*dxidxi2p(2)-dxidxi1p(2)*dxidxi2p(1)
            ! XI->X transformation
            ! Geometrical shape functions and first derivatives at xi
#           define etype e%gtype
#           define delta 0.d0
#           define phi gphi
#           define dphidxi1 dgphidxi1
#           define dphidxi2 dgphidxi2
#           include <phi_and_dphidxik_2d.rc>
#           undef etype
#           undef delta
#           undef phi
#           undef dphidxi1
#           undef dphidxi2
            ! Components calculation of x, T1 and T2 at xi
            x=0.d0
            T1=0.d0
            T2=0.d0
            do kphi=1,e%n_gnodes
              x=x+gphi(kphi)*e%x(:,kphi)
              T1=T1+dgphidxi1(kphi)*e%x(:,kphi)
              T2=T2+dgphidxi2(kphi)*e%x(:,kphi)
            end do
            ! Normal vector as T1 x T2 at xi
            N(1)=T1(2)*T2(3)-T1(3)*T2(2)
            N(2)=T1(3)*T2(1)-T1(1)*T2(3)
            N(3)=T1(1)*T2(2)-T1(2)*T2(1)
            ! Geometric jacobian
            jg=sqrt(dot_product(N,N))
            ! Unit normal
            n=N/jg
            ! Distance vector
            rv=x-x_i
            ! Distance vector norm
            r=sqrt(dot_product(rv,rv))
            d1r1=1.d0/r
            d1r2=d1r1**2
            drdx=rv*d1r1
            drdn=dot_product(drdx,n)
            ! Jacobians * weights
            jw=jg*js*jqt*jt(1)*jt(2)*w(1)*w(2)
            ! FUNCTIONAL SHAPE FUNCTIONS
            ! Functional shape functions (primary variables) at xi
#           define etype e%ptype
#           define delta e%ptype_delta
#           define phi pphi
#           include <phi_2d.rc>
#           undef etype
#           undef delta
#           undef phi
            ! Functional shape functions (secondary variables) at xi
#           define etype e%stype
#           define delta e%stype_delta
#           define phi sphi
#           include <phi_2d.rc>
#           undef etype
#           undef delta
#           undef phi
            ! Functional shape functions * jacobians * weights
            pphijw=pphi*jw
            sphijw=sphi*jw
            ! COMPONENTS OF THE FUNDAMENTAL SOLUTION
            z(1)=-c_im*p%k*r
            call fbem_zexp_decomposed(1,z,EnR)
            EnR(2,:)=EnR(2,:)*d1r1
            EnR(3,:)=EnR(3,:)*d1r2
            fs_P=d1r1+p%P(1)+EnR(2,1)
            fs_Q=d1r2+p%Q(1)+p%Q(2)*EnR(2,1)+EnR(3,1)
            ! Add
            h=h+fs_Q*drdn*pphijw
            g=g+fs_P*sphijw
          end do
        end do
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'n_edges not valid')
    end select
    ! Multiply by constants
    h=-h*c_1_4pi
    g= g*c_1_4pi
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_harpot3d_sbie_ext_st

  recursive subroutine fbem_bem_harpot3d_sbie_ext_adp(e,reverse,xi_s,x_i,p,qsp,ks,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                                !! Element
    logical                            :: reverse                          !! Reverse orientation
    real(kind=real64)                  :: xi_s(2,fbem_n_vertices(e%gtype)) !! Subdivision of the parent element
    real(kind=real64)                  :: x_i(3)                           !! Collocation point position vector
    type(fbem_bem_harpot3d_parameters) :: p                                !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp                              !! Quasi-singular integration parameters
    integer                            :: ks                               !! Current level of subdivisions
    integer                            :: ns                               !! Maximum level of subdivision
    complex(kind=real64)               :: h(e%n_pnodes)                    !! h integration kernels matrix
    complex(kind=real64)               :: g(e%n_snodes)                    !! g integration kernels matrix
    ! Local
    integer              :: gln_near                             ! 1D Gauss-Legendre integ. points required to integrate only  the quasi-singular integrand
    integer              :: gln                                  ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical              :: subdivide                            ! True if subdivision has to be performed
    real(kind=real64)    :: barxi(2)                             ! Nearest element coordinates with respect to collocation point
    real(kind=real64)    :: barxip(2)                            ! Nearest element subdivision local coordinates with respect to collocation point
    real(kind=real64)    :: rmin                                 ! Minimum distance between collocation point and barxi on the element
    real(kind=real64)    :: barr                                 ! Telles jacobian at barxi
    real(kind=real64)    :: cl                                   ! Characteristic length
    real(kind=real64)    :: d                                    ! Normalized distance between collocation point and element subdivision
    integer              :: method                               ! Method used in nearest point algorithm
    real(kind=real64)    :: tmp_xi_s(2,fbem_n_vertices(e%gtype)) ! Subdivision
    real(kind=real64)    :: x_s(3,e%n_gnodes)                    ! Coordinates of the element subdivision
    complex(kind=real64) :: h_tmp(e%n_pnodes)                    ! h integration kernels matrix (temporary)
    complex(kind=real64) :: g_tmp(e%n_snodes)                    ! g integration kernels matrix (temporary)
    ! Initialize
    if (ks.eq.1) then
      h=(0.d0,0.d0)
      g=(0.d0,0.d0)
      select case (fbem_n_vertices(e%gtype))
        case (3)
          xi_s(:,1)=[1.d0,0.d0]
          xi_s(:,2)=[0.d0,1.d0]
          xi_s(:,3)=[0.d0,0.d0]
        case (4)
          xi_s(:,1)=[-1.d0,-1.d0]
          xi_s(:,2)=[ 1.d0,-1.d0]
          xi_s(:,3)=[ 1.d0, 1.d0]
          xi_s(:,4)=[-1.d0, 1.d0]
      end select
      call fbem_nearest_element_point_bem(3,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
      barxip=barxi
    else
      call fbem_obtain_element_subdivision_coordinates(3,e%gtype,e%x,xi_s,x_s)
      cl=fbem_characteristic_length(3,e%gtype,x_s,1.d-12)
      call fbem_nearest_element_point_bem(3,e%gtype,x_s,cl,x_i,barxip,rmin,d,method)
    end if
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,3,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_harpot3d_sbie_ext_adp',ns,'maximum number of subdivisions reached')
        gln_near=30
      end if
    else
      if (gln_near.eq.0) subdivide=.true.
    end if
    ! Subdivide by 1/4
    if (subdivide) then
      select case (fbem_n_vertices(e%gtype))
        ! TRIANGULAR ELEMENTS
        case (3)
          ! SUBTRI 1
          tmp_xi_s(:,1)=xi_s(:,1)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          call fbem_bem_harpot3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
          ! SUBTRI 2
          tmp_xi_s(:,1)=xi_s(:,2)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          call fbem_bem_harpot3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
          ! SUBTRI 3
          tmp_xi_s(:,1)=xi_s(:,3)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          call fbem_bem_harpot3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
          ! SUBTRI 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          call fbem_bem_harpot3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
        ! QUADRILATERALS
        case (4)
          ! SUBQUAD 1
          tmp_xi_s(:,1)=xi_s(:,1)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,3)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          call fbem_bem_harpot3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
          ! SUBQUAD 2
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=xi_s(:,2)
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,4)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          call fbem_bem_harpot3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
          ! SUBQUAD 3
          tmp_xi_s(:,1)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=xi_s(:,3)
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          call fbem_bem_harpot3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
          ! SUBQUAD 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          tmp_xi_s(:,2)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=xi_s(:,4)
          call fbem_bem_harpot3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
      end select
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_harpot3d_sbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,p,gln,h_tmp,g_tmp)
      h=h+h_tmp
      g=g+g_tmp
    end if
  end subroutine fbem_bem_harpot3d_sbie_ext_adp

  !! This subroutine calculates kernels vectors h and g of singular formulation for interior integration.
  subroutine fbem_bem_harpot3d_sbie_int(type_g,type_f1,type_f2,delta_f,x_nodes,reverse,xi_i,p,h,g)
    implicit none
    ! I/O
    integer                            :: type_g                          !! Geometrial interpolation
    integer                            :: type_f1                         !! Functional interpolation (primary variables)
    integer                            :: type_f2                         !! Functional interpolation (secondary variables)
    real(kind=real64)                  :: delta_f                         !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64)                  :: x_nodes(3,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
    logical                            :: reverse                         !! Normal vector inversion
    real(kind=real64)                  :: xi_i(2)                         !! Reference coordinates of the singular point.
    type(fbem_bem_harpot3d_parameters) :: p                               !! Parameters of the region
    complex(kind=real64)               :: h(fbem_n_nodes(type_f1))        !! h kernel vector
    complex(kind=real64)               :: g(fbem_n_nodes(type_f2))        !! g kernel vector
    ! Local
    integer              :: ksubtri                          ! Counter variable for subtriangles loop
    integer              :: nsubtri                          ! Number of subtriangles
    integer              :: subtriangle(8)                   ! Vector that contains what subtriangles need to be integrated
    real(kind=real64)    :: theta_subtri(2,8)                ! Matrix that contains the angles theta of the subtriangles to be integrated
    real(kind=real64)    :: thetap_subtri(2,8)               ! Matrix that contains the angles thetap of the subtriangles to be integrated
    integer              :: ktheta                           ! Counter variable for theta coordinate loop
    integer              :: krho                             ! Counter variable for rho coordinate loop
    integer              :: kphi                             ! Counter coordinates for shape functions loops
    integer              :: nnodes_g                         ! Number of nodes of geometrical interpolation
    integer              :: ngp_theta                        ! Number of Gauss points for theta coordinate
    integer              :: ngp_rho                          ! Number of Gauss points for rho coordinate
    real(kind=real64)    :: thetai, thetaf, thetapi, thetapf ! Initial and final angles for subtriangle integration
    real(kind=real64)    :: w_angular                        ! Weight of the angular coordinate
    real(kind=real64)    :: w_radial                         ! Weight of the radial coordinate
    real(kind=real64)    :: theta                            ! Angle coordinate theta
    real(kind=real64)    :: thetap                           ! Angle coordinate thetap
    real(kind=real64)    :: thetapp                          ! Angle coordinate thetap on [0,1] domain
    real(kind=real64)    :: jthetap                          ! thetap->thetapp jacobian
    real(kind=real64)    :: rhoij                            ! Maximum rho (radial) value for each edge
    real(kind=real64)    :: rho                              ! Radial coordinate rho
    real(kind=real64)    :: rhop                             ! Radial coordinate rho on [0,1] domain
    real(kind=real64)    :: aux(10)                          ! Auxiliary variable for shape functions resources
    real(kind=real64)    :: xi(2)                            ! Reference xi_1,xi_2 coordinates
    real(kind=real64)    :: phi_f1(fbem_n_nodes(type_f1))    ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)    :: phi_f2(fbem_n_nodes(type_f2))    ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)    :: phi_g(fbem_n_nodes(type_g))      ! Geometrical shape functions values at xi_1,xi_2
    real(kind=real64)    :: dphidxi1_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64)    :: dphidxi2_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64)    :: x_i(3)                           ! Collocation point position vector
    real(kind=real64)    :: x(3)                             ! Position vector at xi_1,xi_2
    real(kind=real64)    :: T1(3), T2(3)                     ! Tangent vectors at xi_1,xi_2
    real(kind=real64)    :: N(3)                             ! Normal vector at xi_1,xi_2
    real(kind=real64)    :: rv(3)                            ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, d1r, d1r2                     ! Distance vector module and its inverse
    real(kind=real64)    :: jg                               ! Geometric jacobian
    real(kind=real64)    :: jw                               ! Jacobian * weights
    real(kind=real64)    :: drdn                             ! Partial derivative of r respect to unit normal
    real(kind=real64)    :: phif1jw(fbem_n_nodes(type_f1))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)    :: phif2jw(fbem_n_nodes(type_f2))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)    :: costheta, sintheta               ! cos(theta), sin(theta)
    complex(kind=real64) :: z                                ! z=-i·k·r
    complex(kind=real64) :: E0, E1, E2, E3, E4               ! exp(z) expansion terms
    complex(kind=real64) :: fs_P, fs_Q                       ! Components of the fundamental solutions
    !
    ! Initialization
    !
    ! Kernel vectors
    h=0.d0
    g=0.d0
    ! Number of nodes of geometrical interpolation
    nnodes_g=fbem_n_nodes(type_g)
    ! Calculate real coordinates of collocation point
    ! Geometrical shape functions at xi_i
#   define etype type_g
#   define delta 0.0d0
#   define xi xi_i
#   define phi phi_g
#   include <phi_2d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
    ! Calculate x_i
    x_i=0.d0
    do kphi=1,nnodes_g
      x_i=x_i+phi_g(kphi)*x_nodes(:,kphi)
    end do
    ! Setup of the polar transformation
    call fbem_polar_transformation_setup(type_g,xi_i,nsubtri,subtriangle,theta_subtri,thetap_subtri)
    !
    ! Numerical integration
    !
    ! Loop through triangles
    do ksubtri=1,nsubtri
      ! Initial and final angles of the subtriangle in the theta and thetap space
      thetai=theta_subtri(1,ksubtri)
      thetaf=theta_subtri(2,ksubtri)
      thetapi=thetap_subtri(1,ksubtri)
      thetapf=thetap_subtri(2,ksubtri)
      ! Select the number of Gauss points (this has been determined by studying relative errors over a triangle rectangle)
      ngp_rho=15
      ngp_theta=5+nint(25.0d0*(thetaf-thetai)/c_pi_2)
      ! Loop through angular coordinate
      do ktheta=1,gl01_n(ngp_theta)
        thetapp=gl01_xi(ktheta,ngp_theta)
        w_angular=gl01_w(ktheta,ngp_theta)
        jthetap=(thetapf-thetapi)
        thetap=jthetap*thetapp+thetapi
        ! Angular transformation
        call fbem_polar_transformation_angular(type_g,xi_i,subtriangle(ksubtri),thetap,theta,rhoij)
        ! Save cos(theta) and sin(theta)
        costheta=cos(theta)
        sintheta=sin(theta)
        ! Loop through radial coordinate
        do krho=1,gl01_n(ngp_rho)
          ! Coordinates rhop and rho, and jacobian
          rhop=gl01_xi(krho,ngp_rho)
          w_radial=gl01_w(krho,ngp_rho)
          rho=rhoij*rhop
          ! Coordinates xi1 and xi2, and jacobian
          xi(1)=xi_i(1)+rho*costheta
          xi(2)=xi_i(2)+rho*sintheta
          ! XI->X transformation
          ! Geometrical shape functions and first derivatives at xi
#         define etype type_g
#         define delta 0.0d0
#         define phi phi_g
#         define dphidxi1 dphidxi1_g
#         define dphidxi2 dphidxi2_g
#         include <phi_and_dphidxik_2d.rc>
#         undef etype
#         undef delta
#         undef phi
#         undef dphidxi1
#         undef dphidxi2
          ! Components calculation of x, T1 and T2 at xi
          x=0.d0
          T1=0.d0
          T2=0.d0
          do kphi=1,nnodes_g
            x=x+phi_g(kphi)*x_nodes(:,kphi)
            T1=T1+dphidxi1_g(kphi)*x_nodes(:,kphi)
            T2=T2+dphidxi2_g(kphi)*x_nodes(:,kphi)
          end do
          ! Normal vector as T1 x T2 at xi
          N(1)=T1(2)*T2(3)-T1(3)*T2(2)
          N(2)=T1(3)*T2(1)-T1(1)*T2(3)
          N(3)=T1(1)*T2(2)-T1(2)*T2(1)
          ! Geometric jacobian
          jg=dot_product(N,N)
          jg=sqrt(jg)
          ! Distance vector
          rv=x-x_i
          ! Distance vector norm
          r=dot_product(rv,rv)
          r=sqrt(r)
          d1r=1.d0/r
          d1r2=d1r**2
          ! dr/dn
          drdn=dot_product(rv,N)*d1r/jg
          ! Jacobians * weights
          jw=jg*rho*jthetap*w_angular*w_radial
          ! Components of the fundamental solutions
          z=-c_im*p%k*r
          call fbem_decomposed_zexp(z,E0,E1,E2,E3,E4)
          fs_P=d1r+p%P(1)+d1r*E2
          fs_Q=d1r2+p%Q(1)+p%Q(2)*d1r*E2+d1r2*E3
          ! FUNCTIONAL SHAPE FUNCTIONS
          ! Functional shape functions (primary variables) at xi
#         define etype type_f1
#         define delta delta_f
#         define phi phi_f1
#         include <phi_2d.rc>
#         undef etype
#         undef delta
#         undef phi
          ! Functional shape functions (secondary variables) at xi
#         define etype type_f2
#         define delta delta_f
#         define phi phi_f2
#         include <phi_2d.rc>
#         undef etype
#         undef delta
#         undef phi
          ! Functional shape functions * jacobians * weights
          phif1jw=phi_f1*jw
          phif2jw=phi_f2*jw
          ! Add to kernels
          h=h+fs_Q*drdn*phif1jw
          g=g+fs_P*phif2jw
        end do ! Loop through rho coordinate
      end do ! Loop through theta coordinate
    end do ! Loop through rectangle triangles
    ! Multiply h by -1/(4*pi) and g by 1/(4*pi)
    h=-h*c_1_4pi
    g= g*c_1_4pi
    ! If the normal has to be reversed, then h=-h
    if (reverse) h=-h
  end subroutine fbem_bem_harpot3d_sbie_int

  subroutine fbem_bem_harpot3d_sbie_auto(e,reverse,x_i,p,qsp,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e             !! Integration element
    logical                            :: reverse       !! Reverse orientation
    real(kind=real64)                  :: x_i(3)        !! Collocation point
    type(fbem_bem_harpot3d_parameters) :: p             !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp           !! Quasi-singular integration parameters
    integer                            :: ns            !! Maximum level of subdivisions
    complex(kind=real64)               :: h(e%n_pnodes) !! h integration kernel
    complex(kind=real64)               :: g(e%n_snodes) !! g integration kernel
    ! Local
    real(kind=real64) :: r(3)                             ! Distance vector
    real(kind=real64) :: rmin                             ! Minimum distance between element and x_i
    real(kind=real64) :: barxi(2)                         ! Local coordinates of the nearest element point with respect to x_i
    real(kind=real64) :: d                                ! Dimensionless distance
    integer           :: delta                            ! Control variable
    real(kind=real64) :: xi_s(2,fbem_n_vertices(e%gtype)) ! Local coordinates of the element subdivision
    integer           :: method                           ! Method used when calculating the nearest element point
    integer           :: gln_near                         ! 1D Gauss-Legendre integration points required by the quasi-singular function
    integer           :: gln                              ! 1D Gauss-Legendre integration points used in the integration
    integer           :: ps                               ! Selected precalculated dataset
    integer           :: i                                ! Counter
    ! Determine if interior or exterior integration
    !   - Interior integration (delta=1) requires: xi_i
    !   - Exterior integration (delta=0) requires: x_i, barxi, rmin and d
    ! Use the element ball
    r=e%bball_centre-x_i
    rmin=sqrt(dot_product(r,r))-e%bball_radius
    if (rmin.gt.(4.d0*e%bball_radius)) then
      delta=0
      barxi=0.d0
      d=rmin/e%cl
    else
      ! Use an adaptative algorithm that combines sampling and minimization algorithms
      call fbem_nearest_element_point_bem(3,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
      if (d.le.1.d-12) then
        delta=1
      else
        delta=0
      end if
    end if
    ! Integrate
    select case (delta)
      case (1)
        call fbem_bem_harpot3d_sbie_int(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi,p,h,g)
      case (0)
        ! Estimate the required integration rule
        gln_near=fbem_qs_n_estimation_standard(e%n,e%gtype,3,qsp,d,barxi)
        gln=max(e%gln_far,gln_near)
        ! Integrate using a conservative precalculated dataset
        if ((gln.le.e%ps_gln_max).and.(gln_near.gt.0)) then
          do i=1,e%n_ps
            if (e%ps_gln(i).ge.gln) then
              ps=i
              exit
            end if
          end do
          call fbem_bem_harpot3d_sbie_ext_pre(ps,e,reverse,x_i,p,h,g)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_harpot3d_sbie_ext_adp(e,reverse,xi_s,x_i,p,qsp,1,ns,h,g)
        end if
    end select
  end subroutine fbem_bem_harpot3d_sbie_auto

  subroutine fbem_bem_harpot3d_sbie_bl_auto(e,x_i,p,qsp,ns,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e             !! Integration element
    real(kind=real64)                  :: x_i(3)        !! Collocation point
    type(fbem_bem_harpot3d_parameters) :: p             !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp           !! Quasi-singular integration parameters
    integer                            :: ns            !! Maximum level of subdivisions
    complex(kind=real64)               :: g(e%n_snodes) !! g integration kernel
    ! Local
    real(kind=real64) :: r(3)                             ! Distance vector
    real(kind=real64) :: rmin                             ! Minimum distance between element and x_i
    real(kind=real64) :: barxi(2)                         ! Local coordinates of the nearest element point with respect to x_i
    real(kind=real64) :: d                                ! Dimensionless distance
    integer           :: delta                            ! Control variable
    real(kind=real64) :: xi_s(2,fbem_n_vertices(e%gtype)) ! Local coordinates of the element subdivision
    integer           :: method                           ! Method used when calculating the nearest element point
    integer           :: gln_near                         ! 1D Gauss-Legendre integration points required by the quasi-singular function
    integer           :: gln                              ! 1D Gauss-Legendre integration points used in the integration
    integer           :: ps                               ! Selected precalculated dataset
    integer           :: i                                ! Counter
    ! POINT BODY LOAD
    if (e%d.eq.0) then
      r=e%x(:,1)-x_i
      rmin=sqrt(dot_product(r,r))
      if (rmin.eq.0.d0) then
        call fbem_error_message(output_unit,0,'fbem_bem_harpot3d_sbie_bl_auto',0,'it is not possible to collocate at a point load')
      else
        call fbem_bem_harpot3d_sbie_p(e%x(:,1),x_i,p,g(1))
        return
      end if
    ! LINE, SURFACE OR VOLUME BODY LOAD
    ! Determine if interior or exterior integration
    !   - Interior integration (delta=1) requires: xi_i
    !   - Exterior integration (delta=0) requires: x_i, barxi, rmin and d
    ! Use the element ball
    else
      !
      ! Only point loads
      !
      call fbem_error_message(output_unit,0,'fbem_bem_harpot3d_sbie_bl_auto',0,'only point loads are available')
      r=e%bball_centre-x_i
      rmin=sqrt(dot_product(r,r))-e%bball_radius
      if (rmin.gt.(4.d0*e%bball_radius)) then
        delta=0
        barxi=0.d0
        d=rmin/e%cl
      else
        ! Use an adaptative algorithm that combines sampling and minimization algorithms
        call fbem_nearest_element_point_bem(3,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
        if (d.le.1.d-12) then
          delta=1
        else
          delta=0
        end if
      end if
      ! Integrate
      select case (delta)
        case (1)
          !call fbem_bem_harpot3d_sbie_bl_int(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,barxi,p,g)
        case (0)
          ! Estimate the required integration rule
          gln_near=fbem_qs_n_estimation_standard(e%n,e%gtype,3,qsp,d,barxi)
          gln=max(e%gln_far,gln_near)
          ! Integrate using a conservative precalculated dataset
          if ((gln.le.e%ps_gln_max).and.(gln_near.gt.0)) then
            do i=1,e%n_ps
              if (e%ps_gln(i).ge.gln) then
                ps=i
                exit
              end if
            end do
            !call fbem_bem_harpot3d_sbie_bl_ext_pre(ps,e,x_i,p,g)
          ! Integrate using an adaptative algorithm
          else
            !call fbem_bem_harpot3d_sbie_bl_ext_adp(e,xi_s,x_i,p,qsp,1,ns,g)
          end if
      end select
    end if
  end subroutine fbem_bem_harpot3d_sbie_bl_auto


  ! ================================================================================================================================

  ! ================================================================================================================================
  ! HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (HBIE)





  !! Fundamental solution d*
  subroutine fbem_bem_harpot3d_hbie_d(x,x_i,n_i,p,do)
    implicit none
    ! I/O
    real(kind=real64)                  :: x(3)    !! Observation point
    real(kind=real64)                  :: x_i(3)  !! Collocation point
    real(kind=real64)                  :: n_i(3)  !! Collocation point unit normal
    type(fbem_bem_harpot3d_parameters) :: p       !! Parameters of the region
    complex(kind=real64)               :: do      !! d
    ! Local
    real(kind=real64)    :: rv(3), r, d1r1, d1r2, d1r3, drdx(3), drdni
    complex(kind=real64) :: z(1), EnR(0:6,1)
    rv=x-x_i
    r=sqrt(dot_product(rv,rv))
    d1r1=1.d0/r
    d1r2=d1r1**2
    d1r3=d1r2*d1r1
    drdx=rv*d1r1
    drdni=-dot_product(drdx,n_i)
    z(1)=-c_im*p%k*r
    call fbem_zexp_decomposed(1,z,EnR)
    EnR(2,:)=EnR(2,:)*d1r1
    EnR(3,:)=EnR(3,:)*d1r2
    EnR(4,:)=EnR(4,:)*d1r3
    do=-c_1_4pi*(d1r2+p%Q(1)+p%Q(2)*EnR(2,1)+EnR(3,1))*drdni
  end subroutine fbem_bem_harpot3d_hbie_d

  !! Fundamental solution s*
  subroutine fbem_bem_harpot3d_hbie_s(x,n,x_i,n_i,p,so)
    implicit none
    ! I/O
    real(kind=real64)                  :: x(3)    !! Observation point
    real(kind=real64)                  :: n(3)    !! Observation point unit normal
    real(kind=real64)                  :: x_i(3)  !! Collocation point
    real(kind=real64)                  :: n_i(3)  !! Collocation point unit normal
    type(fbem_bem_harpot3d_parameters) :: p       !! Parameters of the region
    complex(kind=real64)               :: so      !! s*
    ! Local
    real(kind=real64)    :: rv(3), r, d1r1, d1r2, d1r3, drdx(3), drdn, drdni, n_dot_ni
    complex(kind=real64) :: z(1), EnR(0:6,1)
    complex(kind=real64) :: S1, S2
    rv=x-x_i
    r=sqrt(dot_product(rv,rv))
    d1r1=1.d0/r
    d1r2=d1r1**2
    d1r3=d1r2*d1r1
    drdx=rv*d1r1
    drdn=dot_product(drdx,n)
    drdni=-dot_product(drdx,n_i)
    n_dot_ni=dot_product(n,n_i)
    z(1)=-c_im*p%k*r
    call fbem_zexp_decomposed(1,z,EnR)
    EnR(2,:)=EnR(2,:)*d1r1
    EnR(3,:)=EnR(3,:)*d1r2
    EnR(4,:)=EnR(4,:)*d1r3
    S1=p%S1(1)*d1r3+p%S1(2)*d1r1+p%S1(3)*EnR(2,1)+p%S1(4)*EnR(3,1)+p%S1(5)*EnR(4,1)
    S2=d1r3+p%S2(1)*d1r1+p%S2(2)+p%S2(3)*EnR(3,1)+EnR(4,1)
    so=c_1_4pi*(S1*drdn*drdni+S2*n_dot_ni)
  end subroutine fbem_bem_harpot3d_hbie_s

  subroutine fbem_bem_harpot3d_hbie_ext_pre(ps,e,reverse,x_i,n_i,p,m,l)
    implicit none
    ! I/O
    integer                            :: ps            !! Selected precalculated dataset
    type(fbem_bem_element)             :: e             !! Element
    logical                            :: reverse       !! Reverse normal vector
    real(kind=real64)                  :: x_i(3)        !! Collocation point position vector
    real(kind=real64)                  :: n_i(3)        !! Collocation point unit normal vector
    type(fbem_bem_harpot3d_parameters) :: p             !! Parameters of the region
    complex(kind=real64)               :: m(e%n_pnodes) !! m integration kernels vector
    complex(kind=real64)               :: l(e%n_snodes) !! l integration kernels vector
    ! Local
    integer              :: kip                 ! Counter variable for integration points loop
    real(kind=real64)    :: x(3)                ! Position vector at integration point
    real(kind=real64)    :: n(3)                ! Unit normal vector at integration point
    real(kind=real64)    :: pphijw(e%n_pnodes)  ! phi^p * jacobian * weight at integration point
    real(kind=real64)    :: sphijw(e%n_snodes)  ! phi^s * jacobian * weight at integration point
    real(kind=real64)    :: rv(3)               ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, d1r1, d1r2, d1r3 ! Distance vector module and its inverse
    real(kind=real64)    :: drdx(3)             ! Distance vector derivatives with respect to x_k
    real(kind=real64)    :: drdn                ! Partial derivative of r respect to unit normal
    real(kind=real64)    :: drdni               ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)    :: n_dot_ni            ! Dot product of n and n_i
    complex(kind=real64) :: z(1)                ! Arguments z=-ikr
    complex(kind=real64) :: EnR(0:6,1)          ! Exponential function decomposition for each wavenumber
    complex(kind=real64) :: Q, S1, S2           ! Components of the fundamental solutions
    complex(kind=real64) :: fs_d, fs_s          ! Fundamental solutions
    ! Initialize
    m=(0.d0,0.d0)
    l=(0.d0,0.d0)
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      x=e%ps_x(:,kip,ps)
      n=e%ps_n(:,kip,ps)
      pphijw=e%ps_pphijw(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      d1r1=1.d0/r
      d1r2=d1r1**2
      d1r3=d1r2*d1r1
      drdx=rv*d1r1
      drdn=dot_product(drdx,n)
      drdni=-dot_product(drdx,n_i)
      n_dot_ni=dot_product(n,n_i)
      z(1)=-c_im*p%k*r
      call fbem_zexp_decomposed(1,z,EnR)
      EnR(2,:)=EnR(2,:)*d1r1
      EnR(3,:)=EnR(3,:)*d1r2
      EnR(4,:)=EnR(4,:)*d1r3
      Q=d1r2+p%Q(1)+p%Q(2)*EnR(2,1)+EnR(3,1)
      S1=p%S1(1)*d1r3+p%S1(2)*d1r1+p%S1(3)*EnR(2,1)+p%S1(4)*EnR(3,1)+p%S1(5)*EnR(4,1)
      S2=d1r3+p%S2(1)*d1r1+p%S2(2)+p%S2(3)*EnR(3,1)+EnR(4,1)
      fs_d=Q*drdni
      fs_s=S1*drdn*drdni+S2*n_dot_ni
      m=m+fs_s*pphijw
      l=l+fs_d*sphijw
    end do
    ! Multiply by constants
    m= m*c_1_4pi
    l=-l*c_1_4pi
    ! Reverse if needed
    if (reverse) m=-m
  end subroutine fbem_bem_harpot3d_hbie_ext_pre

  subroutine fbem_bem_harpot3d_hbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,p,gln,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                                !! Integration element
    logical                            :: reverse                          !! Reverse normal vector
    real(kind=real64)                  :: xi_s(2,fbem_n_vertices(e%gtype)) !! Coordinates of the subdivision of the parent element
    real(kind=real64)                  :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)                  :: n_i(3)                           !! Unit normal at the collocation point
    real(kind=real64)                  :: barxip(2)                        !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)                  :: barr                             !! Telles jacobian at barxip
    type(fbem_bem_harpot3d_parameters) :: p                                !! Parameters of the region
    integer                            :: gln                              !! 1D Gauss-Legendre number of integration points (<=32)
    complex(kind=real64)               :: m(e%n_pnodes)                    !! m kernel vector
    complex(kind=real64)               :: l(e%n_snodes)                    !! l kernel vector
    ! Local
    integer                      :: kphi                            ! Counter variable for shape functions loops
    integer                      :: k1                              ! Counter variable for reference coordinate xi_1
    integer                      :: k2                              ! Counter variable for reference coordinate xi_2
    real(kind=real64)            :: aux(10)                         ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)                ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi1(e%n_gnodes)           ! Geometrical shape functions derivatives values
    real(kind=real64)            :: dgphidxi2(e%n_gnodes)           ! Geometrical shape functions derivatives values
    real(kind=real64)            :: pphi(e%n_pnodes)                ! Functional shape functions values
    real(kind=real64)            :: sphi(e%n_snodes)                ! Functional shape functions values
    real(kind=real64)            :: gamma(2)                        ! Vector of gamma_1,gamma_2 coordinates
    real(kind=real64)            :: w(2)                            ! Weights of the integration rule
    real(kind=real64)            :: xip(2)                          ! Vector of xip_1,xip_2 coordinates
    real(kind=real64)            :: dxidxi1p(2), dxidxi2p(2)        ! xi derivatives with respect to xip
    real(kind=real64)            :: js                              ! Subdivision jacobian
    real(kind=real64)            :: xi(2)                           ! Vector of xi_1,xi_2 coordinates
    real(kind=real64)            :: xipp(2)                         ! Coordinate xipp used for quadrilateral-triangle transformation
    real(kind=real64)            :: barxipp(2)                      ! Coordinate xipp of collocation point
    real(kind=real64)            :: jqt                             ! Jacobian of the quadrilateral-triangle transformation
    type(fbem_telles_parameters) :: telles_parameters(2)            ! Telles parameters for each coordinate
    real(kind=real64)            :: jt(2)                           ! Telles jacobian for each coordinate: xi_1->gamma_1 and xi_2->gamma_2
    real(kind=real64)            :: x(3)                            ! Position vector at xi_1,xi_2
    real(kind=real64)            :: T1(3), T2(3)                    ! Tangent vectors at xi_1,xi_2
    real(kind=real64)            :: N(3)                            ! Normal vector at xi_1,xi_2
    real(kind=real64)            :: rv(3)                           ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)            :: r, d1r1, d1r2, d1r3             ! Distance vector module and its inverse
    real(kind=real64)            :: drdx(3)                         ! Distance vector derivatives with respect to x_k
    real(kind=real64)            :: jg                              ! Geometric jacobian
    real(kind=real64)            :: jw                              ! Jacobians * weights
    real(kind=real64)            :: drdn                            ! Partial derivative of r respect to unit normal
    real(kind=real64)            :: drdni                           ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)            :: n_dot_ni                        ! Dot product of n and n_i
    real(kind=real64)            :: pphijw(e%n_pnodes)              ! Functional shape functions * jw
    real(kind=real64)            :: sphijw(e%n_snodes)              ! Functional shape functions * jw
    complex(kind=real64)         :: z(1)                            ! Arguments z=-ikr
    complex(kind=real64)         :: EnR(0:6,1)                      ! Exponential function decomposition for each wavenumber
    complex(kind=real64)         :: Q, S1, S2                             ! Components of the fundamental solutions
    complex(kind=real64)         :: fs_d, fs_s                            ! Fundamental solutions
    ! Initialization
    m=(0.d0,0.d0)
    l=(0.d0,0.d0)
    ! Numerical integration
    select case (fbem_n_vertices(e%gtype))
      ! QUADRILATERAL ELEMENTS
      case (4)
        ! Calculate Telles parameters for each direction
        telles_parameters(1)=fbem_telles11_calculate_parameters(barxip(1),barr)
        telles_parameters(2)=fbem_telles11_calculate_parameters(barxip(2),barr)
        ! Loop through INTEGRATION POINTS
        do k1=1,gl11_n(gln)
          ! GAMMA1 COORDINATE
          gamma(1)=gl11_xi(k1,gln)
          w(1)=gl11_w(k1,gln)
          ! GAMMA1->XIP1 TRANSFORMATION
          call fbem_telles_xi_and_jacobian(telles_parameters(1),gamma(1),xip(1),jt(1))
          do k2=1,gl11_n(gln)
            ! GAMMA2 COORDINATE
            gamma(2)=gl11_xi(k2,gln)
            w(2)=gl11_w(k2,gln)
            ! GAMMA2->XIP2 TRANSFORMATION
            call fbem_telles_xi_and_jacobian(telles_parameters(2),gamma(2),xip(2),jt(2))
            ! XIP->XI TRANSFORMATION
            ! Shape functions and its derivatives
#           define delta 0.d0
#           define xi xip
#           define phi gphi
#           define dphidxi1 dgphidxi1
#           define dphidxi2 dgphidxi2
#           include <phi_quad4.rc>
#           include <dphidxi1_quad4.rc>
#           include <dphidxi2_quad4.rc>
#           undef delta
#           undef xi
#           undef phi
#           undef dphidxi1
#           undef dphidxi2
            ! xi coordinates, and xi derivatives
            xi=0.d0
            dxidxi1p=0.d0
            dxidxi2p=0.d0
            do kphi=1,4
              xi=xi+gphi(kphi)*xi_s(:,kphi)
              dxidxi1p=dxidxi1p+dgphidxi1(kphi)*xi_s(:,kphi)
              dxidxi2p=dxidxi2p+dgphidxi2(kphi)*xi_s(:,kphi)
            end do
            ! xip->xi jacobian
            js=dxidxi1p(1)*dxidxi2p(2)-dxidxi1p(2)*dxidxi2p(1)
            ! XI->X TRANSFORMATION
            ! Geometrical shape functions and first derivatives at xi
#           define etype e%gtype
#           define delta 0.d0
#           define phi gphi
#           define dphidxi1 dgphidxi1
#           define dphidxi2 dgphidxi2
#           include <phi_and_dphidxik_2d.rc>
#           undef etype
#           undef delta
#           undef phi
#           undef dphidxi1
#           undef dphidxi2
            ! Components calculation of x, T1 and T2 at xi
            x=0.d0
            T1=0.d0
            T2=0.d0
            do kphi=1,e%n_gnodes
              x=x+gphi(kphi)*e%x(:,kphi)
              T1=T1+dgphidxi1(kphi)*e%x(:,kphi)
              T2=T2+dgphidxi2(kphi)*e%x(:,kphi)
            end do
            ! Normal vector as T1 x T2 at xi
            N(1)=T1(2)*T2(3)-T1(3)*T2(2)
            N(2)=T1(3)*T2(1)-T1(1)*T2(3)
            N(3)=T1(1)*T2(2)-T1(2)*T2(1)
            ! Geometric jacobian
            jg=sqrt(dot_product(N,N))
            ! Unit normal
            n=N/jg
            ! Distance vector
            rv=x-x_i
            ! Distance vector norm
            r=sqrt(dot_product(rv,rv))
            d1r1=1.d0/r
            d1r2=d1r1**2
            d1r3=d1r2*d1r1
            drdx=rv*d1r1
            drdn=dot_product(drdx,n)
            drdni=-dot_product(drdx,n_i)
            n_dot_ni=dot_product(n,n_i)
            ! Jacobians * weights
            jw=jg*js*jt(1)*jt(2)*w(1)*w(2)
            ! FUNCTIONAL SHAPE FUNCTIONS
            ! Functional shape functions (primary variables) at xi
#           define etype e%ptype
#           define delta e%ptype_delta
#           define phi pphi
#           include <phi_2d.rc>
#           undef etype
#           undef delta
#           undef phi
            ! Functional shape functions (secondary variables) at xi
#           define etype e%stype
#           define delta e%stype_delta
#           define phi sphi
#           include <phi_2d.rc>
#           undef etype
#           undef delta
#           undef phi
            ! Functional shape functions * jacobians * weights
            pphijw=pphi*jw
            sphijw=sphi*jw
            ! COMPONENTS OF THE FUNDAMENTAL SOLUTION
            z(1)=-c_im*p%k*r
            call fbem_zexp_decomposed(1,z,EnR)
            EnR(2,:)=EnR(2,:)*d1r1
            EnR(3,:)=EnR(3,:)*d1r2
            EnR(4,:)=EnR(4,:)*d1r3
            Q=d1r2+p%Q(1)+p%Q(2)*EnR(2,1)+EnR(3,1)
            S1=p%S1(1)*d1r3+p%S1(2)*d1r1+p%S1(3)*EnR(2,1)+p%S1(4)*EnR(3,1)+p%S1(5)*EnR(4,1)
            S2=d1r3+p%S2(1)*d1r1+p%S2(2)+p%S2(3)*EnR(3,1)+EnR(4,1)
            fs_d=Q*drdni
            fs_s=S1*drdn*drdni+S2*n_dot_ni
            ! Add
            m=m+fs_s*pphijw
            l=l+fs_d*sphijw
          end do
        end do
       ! TRIANGULAR ELEMENTS
      case (3)
        ! Telles transformation is applied to Gauss-Legendre*Gauss-Legendre quadrature before the quadrilateral->triangle
        ! transformation. Because barxi for triangles are given in area triangle coordinates, for Telles transformation
        ! they must be transformed to quadrilateral coordinates. A special treatment is needed when barxi_2 is near 1, because
        ! transformation diverges.
        if (barxip(2).gt.0.995d0) then
          barxipp(1)=0.5d0
          barxipp(2)=1.0d0
        else
          barxipp(1)=barxip(1)/(1.0d0-barxip(2))
          barxipp(2)=barxip(2)
        end if
        ! Calculate Telles parameters
        telles_parameters(1)=fbem_telles01_calculate_parameters(barxipp(1),barr)
        telles_parameters(2)=fbem_telles01_calculate_parameters(barxipp(2),barr)
        ! Loop through INTEGRATION POINTS
        do k1=1,gl01_n(gln)
          ! GAMMA1 COORDINATE
          gamma(1)=gl01_xi(k1,gln)
          w(1)=gl01_w(k1,gln)
          ! GAMMA1->XIPP1 TRANSFORMATION
          ! xipp_1 coordinate and jacobian from Telles transformation
          call fbem_telles_xi_and_jacobian(telles_parameters(1),gamma(1),xipp(1),jt(1))
          do k2=1,gl01_n(gln)
            ! GAMMA2 COORDINATE
            gamma(2)=gl01_xi(k2,gln)
            w(2)=gl01_w(k2,gln)
            ! GAMMA2->XIPP2 TRANSFORMATION
            call fbem_telles_xi_and_jacobian(telles_parameters(2),gamma(2),xipp(2),jt(2))
            ! XIPP->XIP TRANSFORMATION (QUAD->TRI)
            xip(1)=(1.d0-xipp(2))*xipp(1)
            xip(2)=xipp(2)
            jqt=1.d0-xipp(2)
            ! XIP->XI TRANSFORMATION
            ! Shape functions and its derivatives
#           define delta 0.d0
#           define xi xip
#           define phi gphi
#           define dphidxi1 dgphidxi1
#           define dphidxi2 dgphidxi2
#           include <phi_tri3.rc>
#           include <dphidxi1_tri3.rc>
#           include <dphidxi2_tri3.rc>
#           undef delta
#           undef xi
#           undef phi
#           undef dphidxi1
#           undef dphidxi2
            ! xi coordinates, and xi derivatives
            xi=0.d0
            dxidxi1p=0.d0
            dxidxi2p=0.d0
            do kphi=1,3
              xi=xi+gphi(kphi)*xi_s(:,kphi)
              dxidxi1p=dxidxi1p+dgphidxi1(kphi)*xi_s(:,kphi)
              dxidxi2p=dxidxi2p+dgphidxi2(kphi)*xi_s(:,kphi)
            end do
            ! xip->xi jacobian
            js=dxidxi1p(1)*dxidxi2p(2)-dxidxi1p(2)*dxidxi2p(1)
            ! XI->X transformation
            ! Geometrical shape functions and first derivatives at xi
#           define etype e%gtype
#           define delta 0.d0
#           define phi gphi
#           define dphidxi1 dgphidxi1
#           define dphidxi2 dgphidxi2
#           include <phi_and_dphidxik_2d.rc>
#           undef etype
#           undef delta
#           undef phi
#           undef dphidxi1
#           undef dphidxi2
            ! Components calculation of x, T1 and T2 at xi
            x=0.d0
            T1=0.d0
            T2=0.d0
            do kphi=1,e%n_gnodes
              x=x+gphi(kphi)*e%x(:,kphi)
              T1=T1+dgphidxi1(kphi)*e%x(:,kphi)
              T2=T2+dgphidxi2(kphi)*e%x(:,kphi)
            end do
            ! Normal vector as T1 x T2 at xi
            N(1)=T1(2)*T2(3)-T1(3)*T2(2)
            N(2)=T1(3)*T2(1)-T1(1)*T2(3)
            N(3)=T1(1)*T2(2)-T1(2)*T2(1)
            ! Geometric jacobian
            jg=sqrt(dot_product(N,N))
            ! Unit normal
            n=N/jg
            ! Distance vector
            rv=x-x_i
            ! Distance vector norm
            r=sqrt(dot_product(rv,rv))
            d1r1=1.d0/r
            d1r2=d1r1**2
            d1r3=d1r2*d1r1
            drdx=rv*d1r1
            drdn=dot_product(drdx,n)
            drdni=-dot_product(drdx,n_i)
            n_dot_ni=dot_product(n,n_i)
            ! Jacobians * weights
            jw=jg*js*jqt*jt(1)*jt(2)*w(1)*w(2)
            ! FUNCTIONAL SHAPE FUNCTIONS
            ! Functional shape functions (primary variables) at xi
#           define etype e%ptype
#           define delta e%ptype_delta
#           define phi pphi
#           include <phi_2d.rc>
#           undef etype
#           undef delta
#           undef phi
            ! Functional shape functions (secondary variables) at xi
#           define etype e%stype
#           define delta e%stype_delta
#           define phi sphi
#           include <phi_2d.rc>
#           undef etype
#           undef delta
#           undef phi
            ! Functional shape functions * jacobians * weights
            pphijw=pphi*jw
            sphijw=sphi*jw
            ! COMPONENTS OF THE FUNDAMENTAL SOLUTION
            z(1)=-c_im*p%k*r
            call fbem_zexp_decomposed(1,z,EnR)
            EnR(2,:)=EnR(2,:)*d1r1
            EnR(3,:)=EnR(3,:)*d1r2
            EnR(4,:)=EnR(4,:)*d1r3
            Q=d1r2+p%Q(1)+p%Q(2)*EnR(2,1)+EnR(3,1)
            S1=p%S1(1)*d1r3+p%S1(2)*d1r1+p%S1(3)*EnR(2,1)+p%S1(4)*EnR(3,1)+p%S1(5)*EnR(4,1)
            S2=d1r3+p%S2(1)*d1r1+p%S2(2)+p%S2(3)*EnR(3,1)+EnR(4,1)
            fs_d=Q*drdni
            fs_s=S1*drdn*drdni+S2*n_dot_ni
            ! Add
            m=m+fs_s*pphijw
            l=l+fs_d*sphijw
          end do
        end do
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'n_edges not valid')
    end select
    ! Multiply by constants
    m= m*c_1_4pi
    l=-l*c_1_4pi
    ! Reverse if needed
    if (reverse) m=-m
  end subroutine fbem_bem_harpot3d_hbie_ext_st

  recursive subroutine fbem_bem_harpot3d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,p,qsp,ks,ns,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                                !! Element
    logical                            :: reverse                          !! Reverse orientation
    real(kind=real64)                  :: xi_s(2,fbem_n_vertices(e%gtype)) !! Subdivision of the parent element
    real(kind=real64)                  :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)                  :: n_i(3)                           !! Unit normal at the collocation point
    type(fbem_bem_harpot3d_parameters) :: p                                !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp                              !! Quasi-singular integration parameters
    integer                            :: ks                               !! Current level of subdivisions
    integer                            :: ns                               !! Maximum level of subdivision
    complex(kind=real64)               :: m(e%n_pnodes)                    !! m integration kernels matrix
    complex(kind=real64)               :: l(e%n_snodes)                    !! l integration kernels matrix
    ! Local
    integer              :: gln_near                             ! 1D Gauss-Legendre integ. points required to integrate only  the quasi-singular integrand
    integer              :: gln                                  ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical              :: subdivide                            ! True if subdivision has to be performed
    real(kind=real64)    :: barxi(2)                             ! Nearest element coordinates with respect to collocation point
    real(kind=real64)    :: barxip(2)                            ! Nearest element subdivision local coordinates with respect to collocation point
    real(kind=real64)    :: rmin                                 ! Minimum distance between collocation point and barxi on the element
    real(kind=real64)    :: barr                                 ! Telles jacobian at barxi
    real(kind=real64)    :: cl                                   ! Characteristic length
    real(kind=real64)    :: d                                    ! Normalized distance between collocation point and element subdivision
    integer              :: method                               ! Method used in nearest point algorithm
    real(kind=real64)    :: tmp_xi_s(2,fbem_n_vertices(e%gtype)) ! Subdivision
    real(kind=real64)    :: x_s(3,e%n_gnodes)                    ! Coordinates of the element subdivision
    complex(kind=real64) :: m_tmp(e%n_pnodes)                    ! m integration kernels matrix (temporary)
    complex(kind=real64) :: l_tmp(e%n_snodes)                    ! l integration kernels matrix (temporary)
    ! Initialize
    if (ks.eq.1) then
      m=(0.d0,0.d0)
      l=(0.d0,0.d0)
      select case (fbem_n_vertices(e%gtype))
        case (3)
          xi_s(:,1)=[1.d0,0.d0]
          xi_s(:,2)=[0.d0,1.d0]
          xi_s(:,3)=[0.d0,0.d0]
        case (4)
          xi_s(:,1)=[-1.d0,-1.d0]
          xi_s(:,2)=[ 1.d0,-1.d0]
          xi_s(:,3)=[ 1.d0, 1.d0]
          xi_s(:,4)=[-1.d0, 1.d0]
      end select
      call fbem_nearest_element_point_bem(3,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
      barxip=barxi
    else
      call fbem_obtain_element_subdivision_coordinates(3,e%gtype,e%x,xi_s,x_s)
      cl=fbem_characteristic_length(3,e%gtype,x_s,1.d-12)
      call fbem_nearest_element_point_bem(3,e%gtype,x_s,cl,x_i,barxip,rmin,d,method)
    end if
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,5,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_harpot3d_hbie_ext_adp',ns,'maximum number of subdivisions reached')
        gln_near=30
      end if
    else
      if (gln_near.eq.0) subdivide=.true.
    end if
    ! Subdivide by 1/4
    if (subdivide) then
      select case (fbem_n_vertices(e%gtype))
        ! TRIANGULAR ELEMENTS
        case (3)
          ! SUBTRI 1
          tmp_xi_s(:,1)=xi_s(:,1)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          call fbem_bem_harpot3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
          ! SUBTRI 2
          tmp_xi_s(:,1)=xi_s(:,2)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          call fbem_bem_harpot3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
          ! SUBTRI 3
          tmp_xi_s(:,1)=xi_s(:,3)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          call fbem_bem_harpot3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
          ! SUBTRI 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          call fbem_bem_harpot3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
        ! QUADRILATERALS
        case (4)
          ! SUBQUAD 1
          tmp_xi_s(:,1)=xi_s(:,1)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,3)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          call fbem_bem_harpot3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
          ! SUBQUAD 2
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=xi_s(:,2)
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,4)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          call fbem_bem_harpot3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
          ! SUBQUAD 3
          tmp_xi_s(:,1)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=xi_s(:,3)
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          call fbem_bem_harpot3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
          ! SUBQUAD 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          tmp_xi_s(:,2)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=xi_s(:,4)
          call fbem_bem_harpot3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
      end select
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_harpot3d_hbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,p,gln,m_tmp,l_tmp)
      m=m+m_tmp
      l=l+l_tmp
    end if
  end subroutine fbem_bem_harpot3d_hbie_ext_adp

  subroutine fbem_bem_harpot3d_hbie_int(type_g,type_f1,type_f2,delta_f,x_nodes,reverse,xi_i,p,m,l)
    implicit none
    ! I/O
    integer                            :: type_g                          !! Geometrial interpolation
    integer                            :: type_f1                         !! Functional interpolation (primary variables)
    integer                            :: type_f2                         !! Functional interpolation (secondary variables)
    real(kind=real64)                  :: delta_f                         !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64)                  :: x_nodes(3,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
    logical                            :: reverse                         !! Normal vector inversion
    real(kind=real64)                  :: xi_i(2)                         !! Reference coordinates of the singular point.
    type(fbem_bem_harpot3d_parameters) :: p                               !! Parameters of the region
    complex(kind=real64)               :: m(fbem_n_nodes(type_f1))        !! m kernel vector
    complex(kind=real64)               :: l(fbem_n_nodes(type_f2))        !! l kernel vector
    ! Local
    integer              :: ksubtri                          ! Counter variable for subtriangles loop
    integer              :: nsubtri                          ! Number of subtriangles
    integer              :: subtriangle(8)                   ! Vector that contains subtriangles that need to be integrated
    real(kind=real64)    :: theta_subtri(2,8)
    real(kind=real64)    :: thetap_subtri(2,8)
    integer              :: ktheta                           ! Counter variable for theta coordinate loop
    integer              :: krho                             ! Counter variable for rho coordinate loop
    integer              :: kphi                             ! Counter coordinates for shape functions loops
    integer              :: nnodes_g                         ! Number of nodes of geometrical interpolation
    integer              :: nnodes_f1                        ! Number of nodes of functional (primary variables) interpolation
    integer              :: ngp_theta                        ! Number of Gauss points for theta coordinate
    integer              :: ngp_rho                          ! Number of Gauss points for rho coordinate
    real(kind=real64)    :: thetai, thetaf, thetapi, thetapf ! Initial and ending angles for subtriangle integration
    real(kind=real64)    :: w_angular                        ! Weight of the angular coordinate
    real(kind=real64)    :: w_radial                         ! Weight of the radial coordinate
    real(kind=real64)    :: theta                            ! Angle coordinate theta
    real(kind=real64)    :: thetap                           ! Angle coordinate thetap
    real(kind=real64)    :: thetapp                          ! Angle coordinate thetap on [0,1] domain
    real(kind=real64)    :: jthetap                          ! thetap->thetapp jacobian
    real(kind=real64)    :: rhoij                            ! Maximum rho (radial) value for each edge
    real(kind=real64)    :: rho                              ! Radial coordinate rho
    real(kind=real64)    :: rhop                             ! Radial coordinate rho on [0,1] domain
    real(kind=real64)    :: aux(10)                          ! Auxiliary variable for shape functions resources
    real(kind=real64)    :: xi(2)                            ! Reference xi_1,xi_2 coordinates
    real(kind=real64)    :: phi_f1(fbem_n_nodes(type_f1))    ! Functional shape functions values (primary variables) at xi_1,xi_2
    real(kind=real64)    :: phi_f2(fbem_n_nodes(type_f2))    ! Functional shape functions values (secondary variables) at xi_1,xi_2
    real(kind=real64)    :: phi_f1_i(fbem_n_nodes(type_f1))  ! Functional shape functions values (primary variables) at xi_i1,xi_i2
    real(kind=real64)    :: psi_i(3,fbem_n_nodes(type_f1))   ! psi vectors at xi_i1,xi_i2
    real(kind=real64)    :: phi_g(fbem_n_nodes(type_g))      ! Geometrical shape functions values at xi_1,xi_2
    real(kind=real64)    :: dphidxi1_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64)    :: dphidxi2_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64)    :: jxi1, jxi2                       ! xi1->x, xi2->x jacobians
    real(kind=real64)    :: x_i(3)                           ! Collocation point position vector
    real(kind=real64)    :: n_i(3)                           ! Unit normal at collocation point
    real(kind=real64)    :: x(3)                             ! Position vector at xi_1,xi_2
    real(kind=real64)    :: T1(3), T2(3)                     ! Tangent vectors at xi_1,xi_2
    real(kind=real64)    :: v1(3), v2(3)                     ! Orthogonal tangent vectors
    real(kind=real64)    :: Mvt(2,2)                         ! Orthogonal coordinates transformation matrix
    real(kind=real64)    :: dxidv(2,2)                       ! xi coordinates derivatives with respect to v orthogonal coordinates
    real(kind=real64)    :: detMvt                           ! Mvt(1,1)*Mvt(2,2)-Mvt(1,2)*Mvt(2,1)
    real(kind=real64)    :: n(3)                             ! Unit normal vector at xi_1,xi_2
    real(kind=real64)    :: rv(3)                            ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, d1r, d1r2, d1r3               ! Distance vector module and its inverse
    real(kind=real64)    :: jg                               ! Geometric jacobian
    real(kind=real64)    :: drdn                             ! Partial derivative of r respect to unit normal
    real(kind=real64)    :: drdni                            ! Partial derivative of r respect to unit normal on collocation point
    real(kind=real64)    :: n_dot_ni                         ! Scalar product between integration point unit normal and collocation point unit normal
    real(kind=real64)    :: jw                               ! Jacobians * weights
    real(kind=real64)    :: phif1jw(fbem_n_nodes(type_f1))   ! Functional shape functions values (primary variables) at xi_1,xi_2
    real(kind=real64)    :: phif2jw(fbem_n_nodes(type_f2))   ! Functional shape functions values (secondary variables) at xi_1,xi_2
    real(kind=real64)    :: costheta, sintheta               ! cos(theta), sin(theta)
    complex(kind=real64) :: z                                ! z=-i·k·r
    complex(kind=real64) :: E0, E1, E2, E3, E4               ! exp(z) expansion terms
    complex(kind=real64) :: Q, S1, S2                        ! Components of the fundamental solutions
    complex(kind=real64) :: fs_d, fs_s                       ! Fundamental solutions
    ! Local variables associated with line integrals
    real(kind=real64)              :: re                                   ! Quasi-singular integration relative error
    integer                        :: gln_min                              ! Minimum Number of 1D Gauss-Legendre number of points
    integer                        :: ns_max                               ! Maximum number of subdivisions
    type(fbem_qs_parameters)       :: qsp                                  ! Quasi-singular integration parameters
    real(kind=real64)              :: xi_s(1,2)                            ! Edge subdivision
    integer                        :: kedge                                ! Counter variable for edges
    integer                        :: nedges                               ! Number of edges
    integer                        :: type_edge                            ! Line element type for line integrals
    integer                        :: nnodes_edge                          ! Number of nodes of the edge
    real(kind=real64), allocatable :: x_nodes_edge(:,:)                    ! Coordinates of the edge elements
    real(kind=real64)              :: mli1, mli2(3)                        ! Line integrals values
    !
    ! Initialization
    !
    ! Kernel vectors
    m=0.d0
    l=0.d0
    ! Number of nodes of geometrical interpolation
    nnodes_g=fbem_n_nodes(type_g)
    ! Number of nodes of functional (primary variables) interpolation
    nnodes_f1=fbem_n_nodes(type_f1)
    ! Calculate real coordinates of collocation point, and tangent vectors
    ! Geometrical shape functions and first derivatives at xi_i
#   define etype type_g
#   define delta 0.0d0
#   define xi xi_i
#   define phi phi_g
#   define dphidxi1 dphidxi1_g
#   define dphidxi2 dphidxi2_g
#   include <phi_and_dphidxik_2d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
#   undef dphidxi1
#   undef dphidxi2
    ! Calculate x_i, T1 and T2
    x_i=0.d0
    T1=0.d0
    T2=0.d0
    do kphi=1,nnodes_g
      x_i=x_i+phi_g(kphi)*x_nodes(:,kphi)
      T1=T1+dphidxi1_g(kphi)*x_nodes(:,kphi)
      T2=T2+dphidxi2_g(kphi)*x_nodes(:,kphi)
    end do
    ! Normal vector as T1 x T2 at xi_i
    N(1)=T1(2)*T2(3)-T1(3)*T2(2)
    N(2)=T1(3)*T2(1)-T1(1)*T2(3)
    N(3)=T1(1)*T2(2)-T1(2)*T2(1)
    ! Geometric jacobian
    jg=dot_product(N,N)
    jg=sqrt(jg)
    ! Unit normal
    n_i=N/jg
    ! In jt1 and jt2 it is saved the tangent vectors norm
    jxi1=dot_product(T1,T1)
    jxi1=sqrt(jxi1)
    t1=T1/jxi1
    jxi2=dot_product(T2,T2)
    jxi2=sqrt(jxi2)
    t2=T2/jxi2
    ! Calculate psi shape functions vector (discretized tangential gradient of the primary variable)
    ! Functional shape functions and its first derivatives at xi_i
#   define etype type_f1
#   define delta delta_f
#   define xi xi_i
#   define phi phi_f1_i
#   define dphidxi1 dphidxi1_g
#   define dphidxi2 dphidxi2_g
#   include <phi_and_dphidxik_2d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
#   undef dphidxi1
#   undef dphidxi2
    ! Components calculation of psi at xi_i
    ! The orthogonal system at xi_i
    ! v1=t1
    v1=t1
    ! v2 = ni x t1
    v2(1)=n_i(2)*t1(3)-n_i(3)*t1(2)
    v2(2)=n_i(3)*t1(1)-n_i(1)*t1(3)
    v2(3)=n_i(1)*t1(2)-n_i(2)*t1(1)
    ! Orthogonal coordinates transformation matrix
    Mvt(1,1)=1.d0               ! v1·t1
    Mvt(1,2)=dot_product(v1,t2) ! v1·t2
    Mvt(2,1)=0.d0               ! v2·t1
    Mvt(2,2)=dot_product(v2,t2) ! v2·t2
    ! xi derivatives with respect to surface orthogonal coordinates
    detMvt=Mvt(1,1)*Mvt(2,2)-Mvt(1,2)*Mvt(2,1)
    dxidv(1,1)= Mvt(2,2)/(jxi1*detMvt)
    dxidv(1,2)=-Mvt(1,2)/(jxi1*detMvt)
    dxidv(2,1)=-Mvt(2,1)/(jxi2*detMvt)
    dxidv(2,2)= Mvt(1,1)/(jxi2*detMvt)
    ! Build psi_i
    do kphi=1,nnodes_f1
      psi_i(:,kphi)=(dphidxi1_g(kphi)*dxidv(1,1)+dphidxi2_g(kphi)*dxidv(2,1))*v1+&
                    (dphidxi1_g(kphi)*dxidv(1,2)+dphidxi2_g(kphi)*dxidv(2,2))*v2
    end do
    ! Check if xi_i is at vertex
    if (fbem_check_xi1xi2_edge(type_g,xi_i).eqv.(.true.)) then
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                              'the HBIE cannot be collocated at a vertex')
    end if
    ! Setup of the polar transformation
    call fbem_polar_transformation_setup(type_g,xi_i,nsubtri,subtriangle,theta_subtri,thetap_subtri)
    !
    ! Numerical integration
    !
    !
    ! SURFACE INTEGRALS (all regular or weakly singular)
    !
    ! Loop through triangles
    do ksubtri=1,nsubtri
      ! Initial and final angles of the subtriangle in the theta and thetap space
      thetai=theta_subtri(1,ksubtri)
      thetaf=theta_subtri(2,ksubtri)
      thetapi=thetap_subtri(1,ksubtri)
      thetapf=thetap_subtri(2,ksubtri)
      ! Select the number of Gauss points
      ngp_rho=15
      ngp_theta=5+nint(25.0d0*(thetaf-thetai)/c_pi_2)
      ! Loop through theta coordinate
      do ktheta=1,gl01_n(ngp_theta)
        thetapp=gl01_xi(ktheta,ngp_theta)
        w_angular=gl01_w(ktheta,ngp_theta)
        jthetap=(thetapf-thetapi)
        thetap=jthetap*thetapp+thetapi
        ! Angular transformation
        call fbem_polar_transformation_angular(type_g,xi_i,subtriangle(ksubtri),thetap,theta,rhoij)
        ! Save cos(theta) and sin(theta)
        costheta=cos(theta)
        sintheta=sin(theta)
        ! Loop through rho coordinate
        do krho=1,gl01_n(ngp_rho)
         ! Coordinates rhop and rho, and jacobian
          rhop=gl01_xi(krho,ngp_rho)
          w_radial=gl01_w(krho,ngp_rho)
          rho=rhoij*rhop
          ! Coordinates xi1 and xi2, and jacobian
          xi(1)=xi_i(1)+rho*costheta
          xi(2)=xi_i(2)+rho*sintheta
          ! XI->X transformation
          ! Geometrical shape functions and first derivatives at xi
#         define etype type_g
#         define delta 0.0d0
#         define phi phi_g
#         define dphidxi1 dphidxi1_g
#         define dphidxi2 dphidxi2_g
#         include <phi_and_dphidxik_2d.rc>
#         undef etype
#         undef delta
#         undef phi
#         undef dphidxi1
#         undef dphidxi2
          ! Components calculation of x, T1 and T2 at xi
          x=0.d0
          T1=0.d0
          T2=0.d0
          do kphi=1,nnodes_g
            x=x+phi_g(kphi)*x_nodes(:,kphi)
            T1=T1+dphidxi1_g(kphi)*x_nodes(:,kphi)
            T2=T2+dphidxi2_g(kphi)*x_nodes(:,kphi)
          end do
          ! Normal vector as T1 x T2 at xi
          N(1)=T1(2)*T2(3)-T1(3)*T2(2)
          N(2)=T1(3)*T2(1)-T1(1)*T2(3)
          N(3)=T1(1)*T2(2)-T1(2)*T2(1)
          ! Geometric jacobian
          jg=dot_product(N,N)
          jg=sqrt(jg)
          ! Unit normal
          n=N/jg
          ! Distance vector
          rv=x-x_i
          ! Distance vector norm
          r=dot_product(rv,rv)
          r=sqrt(r)
          d1r=1.d0/r
          d1r2=d1r**2
          d1r3=d1r2*d1r
          ! dr/dn
          drdn=dot_product(rv,n)*d1r
          ! dr/dni
          drdni=-dot_product(rv,n_i)*d1r
          ! n·n_i
          n_dot_ni=dot_product(n,n_i)
          ! Jacobians * weights
          jw=jg*rho*jthetap*w_angular*w_radial
          ! Components of the fundamental solutions
          z=-c_im*p%k*r
          call fbem_decomposed_zexp(z,E0,E1,E2,E3,E4)
          Q=d1r2+p%Q(1)+p%Q(2)*d1r*E2+d1r2*E3
          S1=p%S1(2)*d1r+p%S1(3)*d1r*E2+p%S1(4)*d1r2*E3+p%S1(5)*d1r3*E4
          S2=p%S2(1)*d1r+p%S2(2)+p%S2(3)*d1r2*E3+d1r3*E4
          ! Fundamental solution
          fs_d=Q*drdni
          fs_s=S1*drdn*drdni+S2*n_dot_ni
          ! FUNCTIONAL SHAPE FUNCTIONS
          ! Functional shape functions (primary variables) at xi
#         define etype type_f1
#         define delta delta_f
#         define phi phi_f1
#         include <phi_2d.rc>
#         undef etype
#         undef delta
#         undef phi
          ! Functional shape functions (secondary variables) at xi
#         define etype type_f2
#         define delta delta_f
#         define phi phi_f2
#         include <phi_2d.rc>
#         undef etype
#         undef delta
#         undef phi
          ! Functional shape functions * jacobians * weights
          phif1jw=phi_f1*jw
          phif2jw=phi_f2*jw
          ! Add to kernels
          m=m+3.0d0*d1r3*drdn*drdni*(phi_f1-phi_f1_i)*jw
          m=m+d1r3*n_dot_ni*(phi_f1-phi_f1_i-psi_i(1,:)*rv(1)-psi_i(2,:)*rv(2)-psi_i(3,:)*rv(3))*jw
          m=m-(psi_i(1,:)*n(1)+psi_i(2,:)*n(2)+psi_i(3,:)*n(3))*d1r2*drdni*jw
          m=m+fs_s*phif2jw
          l=l+fs_d*drdni*phif2jw
        end do ! Loop through rho coordinate
      end do ! Loop through theta coordinate
    end do ! Loop through triangles
    !
    ! LINE INTEGRALS (all quasi-singular integrals, all integrations performed with an adaptive Telles+subdivision)
    !
    ! Initialize
    mli1=0.d0
    mli2=0.d0
    re=1.d-15
    gln_min=5
    ns_max=16
    call fbem_qs_calculate_parameters(re,qsp)
    ! Number of edges of the element
    nedges=fbem_n_edges(type_g)
    ! Loop through edges
    do kedge=1,nedges
      ! Type of edge
      type_edge=fbem_edge_type(kedge,type_g)
      ! Number of nodes
      nnodes_edge=fbem_n_nodes(type_edge)
      ! Allocate
      allocate (x_nodes_edge(3,nnodes_edge))
      ! Copy its coordinates
      do kphi=1,nnodes_edge
        x_nodes_edge(:,kphi)=x_nodes(:,fbem_edge_node(kphi,kedge,type_g))
      end do
      ! Calculate the line integral (r x ni)·t/r**3 and (e_k x ni)·t/r
      call fbem_bem_stapot3d_hbie_int_li(type_edge,x_nodes_edge,xi_s,x_i,n_i,gln_min,qsp,1,ns_max,mli1,mli2)
      ! Deallocate the edge coordinates
      deallocate (x_nodes_edge)
    end do ! Loop through edges
    !
    ! Add line integrals
    !
    m=m+phi_f1_i*mli1+psi_i(1,:)*mli2(1)+psi_i(2,:)*mli2(2)+psi_i(3,:)*mli2(3)
    ! Multiply m by 1/(4*pi) and l by -1/(4*pi)
    m= m*c_1_4pi
    l=-l*c_1_4pi
    ! If the normal has to be reversed, then l=-l
    if (reverse) l=-l
  end subroutine fbem_bem_harpot3d_hbie_int

  subroutine fbem_bem_harpot3d_hbie_auto(e,reverse,x_i,n_i,p,qsp,ns,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e             !! Integration element
    logical                            :: reverse       !! Reverse orientation
    real(kind=real64)                  :: x_i(3)        !! Collocation point
    real(kind=real64)                  :: n_i(3)        !! Unit normal at the collocation point
    type(fbem_bem_harpot3d_parameters) :: p             !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp           !! Quasi-singular integration parameters
    integer                            :: ns            !! Maximum level of subdivisions
    complex(kind=real64)               :: m(e%n_pnodes) !! m integration kernel
    complex(kind=real64)               :: l(e%n_snodes) !! l integration kernel
    ! Local
    real(kind=real64) :: r(3)                             ! Distance vector
    real(kind=real64) :: rmin                             ! Minimum distance between element and x_i
    real(kind=real64) :: barxi(2)                         ! Local coordinates of the nearest element point with respect to x_i
    real(kind=real64) :: d                                ! Dimensionless distance
    integer           :: delta                            ! Control variable
    real(kind=real64) :: xi_s(2,fbem_n_vertices(e%gtype)) ! Local coordinates of the element subdivision
    integer           :: method                           ! Method used when calculating the nearest element point
    integer           :: gln_near                         ! 1D Gauss-Legendre integration points required by the quasi-singular function
    integer           :: gln                              ! 1D Gauss-Legendre integration points used in the integration
    integer           :: ps                               ! Selected precalculated dataset
    integer           :: i                                ! Counter
    ! Determine if interior or exterior integration
    !   - Interior integration (delta=1) requires: xi_i
    !   - Exterior integration (delta=0) requires: x_i, barxi, rmin and d
    ! Use the element ball
    r=e%bball_centre-x_i
    rmin=sqrt(dot_product(r,r))-e%bball_radius
    if (rmin.gt.(4.d0*e%bball_radius)) then
      delta=0
      barxi=0.d0
      d=rmin/e%cl
    else
      ! Use an adaptative algorithm that combines sampling and minimization algorithms
      call fbem_nearest_element_point_bem(3,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
      if (d.le.1.d-12) then
        delta=1
      else
        delta=0
      end if
    end if
    ! Integrate
    select case (delta)
      case (1)
        call fbem_bem_harpot3d_hbie_int(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi,p,m,l)
      case (0)
        ! Estimate the required integration rule
        gln_near=fbem_qs_n_estimation_standard(e%n,e%gtype,5,qsp,d,barxi)
        gln=max(e%gln_far,gln_near)
        ! Integrate using a conservative precalculated dataset
        if ((gln.le.e%ps_gln_max).and.(gln_near.gt.0)) then
          do i=1,e%n_ps
            if (e%ps_gln(i).ge.gln) then
              ps=i
              exit
            end if
          end do
          call fbem_bem_harpot3d_hbie_ext_pre(ps,e,reverse,x_i,n_i,p,m,l)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_harpot3d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,p,qsp,1,ns,m,l)
        end if
    end select
  end subroutine fbem_bem_harpot3d_hbie_auto

  subroutine fbem_bem_harpot3d_hbie_bl_auto(e,x_i,n_i,p,qsp,ns,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e             !! Integration element
    real(kind=real64)                  :: x_i(3)        !! Collocation point
    real(kind=real64)                  :: n_i(3)        !! Unit normal at the collocation point
    type(fbem_bem_harpot3d_parameters) :: p             !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp           !! Quasi-singular integration parameters
    integer                            :: ns            !! Maximum level of subdivisions
    complex(kind=real64)               :: l(e%n_snodes) !! l integration kernel
    ! Local
    real(kind=real64) :: r(3)                             ! Distance vector
    real(kind=real64) :: rmin                             ! Minimum distance between element and x_i
    real(kind=real64) :: barxi(2)                         ! Local coordinates of the nearest element point with respect to x_i
    real(kind=real64) :: d                                ! Dimensionless distance
    integer           :: delta                            ! Control variable
    real(kind=real64) :: xi_s(2,fbem_n_vertices(e%gtype)) ! Local coordinates of the element subdivision
    integer           :: method                           ! Method used when calculating the nearest element point
    integer           :: gln_near                         ! 1D Gauss-Legendre integration points required by the quasi-singular function
    integer           :: gln                              ! 1D Gauss-Legendre integration points used in the integration
    integer           :: ps                               ! Selected precalculated dataset
    integer           :: i                                ! Counter
    ! POINT BODY LOAD
    if (e%d.eq.0) then
      r=e%x(:,1)-x_i
      rmin=sqrt(dot_product(r,r))
      if (rmin.eq.0.d0) then
        call fbem_error_message(output_unit,0,'fbem_bem_harpot3d_hbie_bl_auto',0,'it is not possible to collocate at a point load')
      else
        call fbem_bem_harpot3d_hbie_d(e%x(:,1),x_i,n_i,p,l(1))
        return
      end if
    ! LINE, SURFACE OR VOLUME BODY LOAD
    ! Determine if interior or exterior integration
    !   - Interior integration (delta=1) requires: xi_i
    !   - Exterior integration (delta=0) requires: x_i, barxi, rmin and d
    ! Use the element ball
    else
      !
      ! Only point loads
      !
      call fbem_error_message(output_unit,0,'fbem_bem_harpot3d_hbie_bl_auto',0,'only point loads are available')
      r=e%bball_centre-x_i
      rmin=sqrt(dot_product(r,r))-e%bball_radius
      if (rmin.gt.(4.d0*e%bball_radius)) then
        delta=0
        barxi=0.d0
        d=rmin/e%cl
      else
        ! Use an adaptative algorithm that combines sampling and minimization algorithms
        call fbem_nearest_element_point_bem(3,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
        if (d.le.1.d-12) then
          delta=1
        else
          delta=0
        end if
      end if
      ! Integrate
      select case (delta)
        case (1)
          !call fbem_bem_harpot3d_hbie_bl_int(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,barxi,p,l)
        case (0)
          ! Estimate the required integration rule
          gln_near=fbem_qs_n_estimation_standard(e%n,e%gtype,5,qsp,d,barxi)
          gln=max(e%gln_far,gln_near)
          ! Integrate using a conservative precalculated dataset
          if ((gln.le.e%ps_gln_max).and.(gln_near.gt.0)) then
            do i=1,e%n_ps
              if (e%ps_gln(i).ge.gln) then
                ps=i
                exit
              end if
            end do
            !call fbem_bem_harpot3d_hbie_bl_ext_pre(ps,e,x_i,n_i,p,l)
          ! Integrate using an adaptative algorithm
          else
            !call fbem_bem_harpot3d_hbie_bl_ext_adp(e,xi_s,x_i,n_i,p,qsp,1,ns,l)
          end if
      end select
    end if
  end subroutine fbem_bem_harpot3d_hbie_bl_auto

  ! ================================================================================================================================

  ! ================================================================================================================================
  ! BOTH SBIE AND HBIE SIMULTANEOUSLY

  subroutine fbem_bem_harpot3d_shbie_ext_pre(ps,e,reverse,x_i,n_i,p,h,g,m,l)
    implicit none
    ! I/O
    integer                            :: ps            !! Selected precalculated dataset
    type(fbem_bem_element)             :: e             !! Element
    logical                            :: reverse       !! Reverse normal vector
    real(kind=real64)                  :: x_i(3)        !! Collocation point position vector
    real(kind=real64)                  :: n_i(3)        !! Collocation point unit normal vector
    type(fbem_bem_harpot3d_parameters) :: p             !! Parameters of the region
    complex(kind=real64)               :: h(e%n_pnodes) !! h integration kernels vector
    complex(kind=real64)               :: g(e%n_snodes) !! g integration kernels vector
    complex(kind=real64)               :: m(e%n_pnodes) !! m integration kernels vector
    complex(kind=real64)               :: l(e%n_snodes) !! l integration kernels vector
    ! Local
    integer              :: kip                 ! Counter variable for integration points loop
    real(kind=real64)    :: x(3)                ! Position vector at integration point
    real(kind=real64)    :: n(3)                ! Unit normal vector at integration point
    real(kind=real64)    :: pphijw(e%n_pnodes)  ! phi^p * jacobian * weight at integration point
    real(kind=real64)    :: sphijw(e%n_snodes)  ! phi^s * jacobian * weight at integration point
    real(kind=real64)    :: rv(3)               ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, d1r1, d1r2, d1r3 ! Distance vector module and its inverse
    real(kind=real64)    :: drdx(3)             ! Distance vector derivatives with respect to x_k
    real(kind=real64)    :: drdn                ! Partial derivative of r respect to unit normal
    real(kind=real64)    :: drdni               ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)    :: n_dot_ni            ! Dot product of n and n_i
    complex(kind=real64) :: z(1)                ! Arguments z=-ikr
    complex(kind=real64) :: EnR(0:6,1)          ! Exponential function decomposition for each wavenumber
    complex(kind=real64) :: PP, Q, S1, S2       ! Components of the fundamental solutions
    complex(kind=real64) :: fs_p, fs_q          ! Fundamental solutions
    complex(kind=real64) :: fs_d, fs_s         ! Fundamental solutions
    ! Initialize
    h=(0.d0,0.d0)
    g=(0.d0,0.d0)
    m=(0.d0,0.d0)
    l=(0.d0,0.d0)
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      x=e%ps_x(:,kip,ps)
      n=e%ps_n(:,kip,ps)
      pphijw=e%ps_pphijw(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      d1r1=1.d0/r
      d1r2=d1r1**2
      d1r3=d1r2*d1r1
      drdx=rv*d1r1
      drdn=dot_product(drdx,n)
      drdni=-dot_product(drdx,n_i)
      n_dot_ni=dot_product(n,n_i)
      z(1)=-c_im*p%k*r
      call fbem_zexp_decomposed(1,z,EnR)
      PP=d1r1+p%P(1)+d1r1*EnR(2,1)
      Q=d1r2+p%Q(1)+p%Q(2)*d1r1*EnR(2,1)+d1r2*EnR(3,1)
      S1=p%S1(1)*d1r3+p%S1(2)*d1r1+p%S1(3)*d1r1*EnR(2,1)+p%S1(4)*d1r2*EnR(3,1)+p%S1(5)*d1r3*EnR(4,1)
      S2=d1r3+p%S2(1)*d1r1+p%S2(2)+p%S2(3)*d1r2*EnR(3,1)+d1r3*EnR(4,1)
      fs_p=PP
      fs_q=Q*drdn
      fs_d=Q*drdni
      fs_s=S1*drdn*drdni+S2*n_dot_ni
      h=h+fs_q*pphijw
      g=g+fs_p*sphijw
      m=m+fs_s*pphijw
      l=l+fs_d*sphijw
    end do
    ! Multiply by constants
    h=-h*c_1_4pi
    g= g*c_1_4pi
    m= m*c_1_4pi
    l=-l*c_1_4pi
    ! Reverse if needed
    if (reverse) then
      h=-h
      m=-m
    end if
  end subroutine fbem_bem_harpot3d_shbie_ext_pre

  subroutine fbem_bem_harpot3d_shbie_auto(e,reverse,x_i,n_i,p,qsp,ns,h,g,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e             !! Integration element
    logical                            :: reverse       !! Reverse orientation
    real(kind=real64)                  :: x_i(3)        !! Collocation point
    real(kind=real64)                  :: n_i(3)        !! Unit normal at the collocation point
    type(fbem_bem_harpot3d_parameters) :: p             !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp           !! Quasi-singular integration parameters
    integer                            :: ns            !! Maximum level of subdivisions
    complex(kind=real64)               :: h(e%n_pnodes) !! h integration kernel
    complex(kind=real64)               :: g(e%n_snodes) !! g integration kernel
    complex(kind=real64)               :: m(e%n_pnodes) !! m integration kernel
    complex(kind=real64)               :: l(e%n_snodes) !! l integration kernel
    ! Local
    real(kind=real64) :: r(3)                             ! Distance vector
    real(kind=real64) :: rmin                             ! Minimum distance between element and x_i
    real(kind=real64) :: barxi(2)                         ! Local coordinates of the nearest element point with respect to x_i
    real(kind=real64) :: d                                ! Dimensionless distance
    integer           :: delta                            ! Control variable
    real(kind=real64) :: xi_s(2,fbem_n_vertices(e%gtype)) ! Local coordinates of the element subdivision
    integer           :: method                           ! Method used when calculating the nearest element point
    integer           :: gln_near                         ! 1D Gauss-Legendre integration points required by the quasi-singular function
    integer           :: gln                              ! 1D Gauss-Legendre integration points used in the integration
    integer           :: ps                               ! Selected precalculated dataset
    integer           :: i                                ! Counter
    ! Determine if interior or exterior integration
    !   - Interior integration (delta=1) requires: xi_i
    !   - Exterior integration (delta=0) requires: x_i, barxi, rmin and d
    ! Use the element ball
    r=e%bball_centre-x_i
    rmin=sqrt(dot_product(r,r))-e%bball_radius
    if (rmin.gt.(4.d0*e%bball_radius)) then
      delta=0
      barxi=0.d0
      d=rmin/e%cl
    else
      ! Use an adaptative algorithm that combines sampling and minimization algorithms
      call fbem_nearest_element_point_bem(3,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
      if (d.le.1.d-12) then
        delta=1
      else
        delta=0
      end if
    end if
    ! Integrate
    select case (delta)
      case (1)
        call fbem_bem_harpot3d_sbie_int(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi,p,h,g)
        call fbem_bem_harpot3d_hbie_int(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi,p,m,l)
      case (0)
        ! Estimate the required integration rule
        gln_near=fbem_qs_n_estimation_standard(e%n,e%gtype,5,qsp,d,barxi)
        gln=max(e%gln_far,gln_near)
        ! Integrate using a conservative precalculated dataset
        if ((gln.le.e%ps_gln_max).and.(gln_near.gt.0)) then
          do i=1,e%n_ps
            if (e%ps_gln(i).ge.gln) then
              ps=i
              exit
            end if
          end do
          call fbem_bem_harpot3d_shbie_ext_pre(ps,e,reverse,x_i,n_i,p,h,g,m,l)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_harpot3d_sbie_ext_adp(e,reverse,xi_s,x_i,p,qsp,1,ns,h,g)
          call fbem_bem_harpot3d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,p,qsp,1,ns,m,l)
        end if
    end select
  end subroutine fbem_bem_harpot3d_shbie_auto

  ! ================================================================================================================================

end module fbem_bem_harpot3d
