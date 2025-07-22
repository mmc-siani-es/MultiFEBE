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
!! <b> This module implements the calculation of element-wise integrals of Boundary Integral Equations of the 2D Helmholtz
!! problem. </b>
module fbem_bem_harpot2d

  ! Fortran 2003 standard
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_quad_rules
  use fbem_telles_transformation
  use fbem_quasisingular_integration
  use fbem_geometry
  use fbem_bem_general

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all are private
  private

  ! ================================================================================================================================
  ! INITIAL SETUP
  public :: fbem_bem_harpot2d_parameters
  public :: fbem_bem_harpot2d_calculate_parameters
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE)
  ! Fundamental solution
  public :: fbem_bem_harpot2d_sbie_p
  public :: fbem_bem_harpot2d_sbie_q
  ! BOUNDARY ELEMENTS
  ! Exterior integration
  public :: fbem_bem_harpot2d_sbie_ext_pre
  public :: fbem_bem_harpot2d_sbie_ext_st
  public :: fbem_bem_harpot2d_sbie_ext_adp
  ! Interior integration
  public :: fbem_bem_harpot2d_sbie_int
  ! Automatic integration
  public :: fbem_bem_harpot2d_sbie_auto
  ! BODY LOADS
  public :: fbem_bem_harpot2d_sbie_bl_auto
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (HBIE)
  ! Fundamental solution
  public :: fbem_bem_harpot2d_hbie_d
  public :: fbem_bem_harpot2d_hbie_s
  ! BOUNDARY ELEMENTS
  ! Exterior integration
  public :: fbem_bem_harpot2d_hbie_ext_pre
  public :: fbem_bem_harpot2d_hbie_ext_st
  public :: fbem_bem_harpot2d_hbie_ext_adp
  ! Interior integration
  public :: fbem_bem_harpot2d_hbie_int
  ! Automatic integration
  public :: fbem_bem_harpot2d_hbie_auto
  ! BODY LOADS
  public :: fbem_bem_harpot2d_hbie_bl_auto
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! BOTH SBIE AND HBIE SIMULTANEOUSLY
  ! BOUNDARY ELEMENTS
  ! Exterior integration
  public :: fbem_bem_harpot2d_shbie_ext_pre
  ! Automatic integration
  public :: fbem_bem_harpot2d_shbie_auto
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! VARIATION SINGULAR BOUNDARY INTEGRAL EQUATION (VSBIE)
  ! BOUNDARY ELEMENTS
  ! Automatic integration
  public :: fbem_bem_harpot2d_vsbie_auto
  ! Exterior integration
  public :: fbem_bem_harpot2d_vsbie_ext_pre
  public :: fbem_bem_harpot2d_vsbie_ext_st
  public :: fbem_bem_harpot2d_vsbie_ext_adp
  ! Interior integration
  public :: fbem_bem_harpot2d_vsbie_int
  ! ================================================================================================================================

  ! ================================================================================================================================
  !! Region parameters (material and frequency dependant)
  type fbem_bem_harpot2d_parameters
    ! Region parameters
    real(kind=real64)    :: rho
    complex(kind=real64) :: c
    real(kind=real64)    :: J
    ! Frequency
    real(kind=real64)    :: omega
    ! Wavenumber
    complex(kind=real64) :: k
    ! Coefficients of the components of the fundamental solution
    complex(kind=real64) :: P(3)
    complex(kind=real64) :: Q(4)
    complex(kind=real64) :: S1(5), S2(4)
  end type fbem_bem_harpot2d_parameters
  ! ================================================================================================================================

contains

  ! ================================================================================================================================
  ! INITIAL SETUP

  subroutine fbem_bem_harpot2d_calculate_parameters(rho,c,omega,pars)
    implicit none
    ! I/O
    real(kind=real64)                  :: rho
    complex(kind=real64)               :: c
    real(kind=real64)                  :: omega
    type(fbem_bem_harpot2d_parameters) :: pars
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
    pars%J=1.0d0/(rho*(omega**2))
    !
    ! Coefficients of equations of fundamental solution's components
    !
    ! P
    pars%P(1)=-c_1_2pi
    pars%P(2)=-c_1_2pi*(log(0.5d0*c_im*k)+c_gamma)
    pars%P(3)= c_1_2pi
    ! Q
    pars%Q(1)=-c_1_2pi
    pars%Q(2)= c_1_2pi*0.5d0*k**2
    pars%Q(3)= c_1_2pi*0.5d0*k**2*(log(0.5d0*c_im*k)+c_gamma-0.5d0)
    pars%Q(4)=-c_1_2pi*c_im*k
    ! S1
    pars%S1(1)= c_1_2pi*2.0d0
    pars%S1(2)= c_1_2pi*0.5d0*k**2
    pars%S1(3)=-c_1_2pi*0.125d0*k**4
    pars%S1(4)=-c_1_2pi*0.125d0*k**4*(log(0.5d0*c_im*k)+c_gamma-0.75d0)
    pars%S1(5)=-c_1_2pi*k**2
    ! S2
    pars%S2(1)= c_1_2pi
    pars%S2(2)=-c_1_2pi*0.5d0*k**2
    pars%S2(3)=-c_1_2pi*0.5d0*k**2*(log(0.5d0*c_im*k)+c_gamma-0.5d0)
    pars%S2(4)= c_1_2pi*c_im*k
  end subroutine fbem_bem_harpot2d_calculate_parameters

  ! ================================================================================================================================

  ! ================================================================================================================================
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE)

  !! Fundamental solution p*
  subroutine fbem_bem_harpot2d_sbie_p(x,x_i,pars,po)
    implicit none
    ! I/O
    real(kind=real64)                  :: x(2)    !! Observation point
    real(kind=real64)                  :: x_i(2)  !! Collocation point
    type(fbem_bem_harpot2d_parameters) :: pars    !! Parameters of the region
    complex(kind=real64)               :: po      !! p*
    ! Local
    real(kind=real64)    :: rv(2), r, d1r, logr, drdx(2)
    complex(kind=real64) :: z(1), KnR(0:2,1)
    rv=x-x_i
    r=sqrt(dot_product(rv,rv))
    d1r=1.0d0/r
    logr=log(r)
    drdx=rv*d1r
    z(1)=c_im*pars%k*r
    call fbem_BesselKnR_decomposed(1,z,KnR)
    po=pars%P(1)*logr+pars%P(2)+pars%P(3)*KnR(0,1)
  end subroutine fbem_bem_harpot2d_sbie_p

  !! Fundamental solution q*
  subroutine fbem_bem_harpot2d_sbie_q(x,n,x_i,pars,qo)
    implicit none
    ! I/O
    real(kind=real64)                  :: x(2)    !! Observation point
    real(kind=real64)                  :: n(2)    !! Observation point normal
    real(kind=real64)                  :: x_i(2)  !! Collocation point
    type(fbem_bem_harpot2d_parameters) :: pars    !! Parameters of the region
    complex(kind=real64)               :: qo      !! q*
    ! Local
    real(kind=real64)    :: rv(2), r, d1r, logr, drdx(2), drdn
    complex(kind=real64) :: z(1), KnR(0:2,1)
    rv=x-x_i
    r=sqrt(dot_product(rv,rv))
    d1r=1.0d0/r
    logr=log(r)
    drdx=rv*d1r
    drdn=dot_product(drdx,n)
    z(1)=c_im*pars%k*r
    call fbem_BesselKnR_decomposed(1,z,KnR)
    qo=(pars%Q(1)*d1r+(pars%Q(2)*logr+pars%Q(3))*r+pars%Q(4)*KnR(1,1))*drdn
  end subroutine fbem_bem_harpot2d_sbie_q

  subroutine fbem_bem_harpot2d_sbie_ext_pre(ps,e,reverse,x_i,pars,h,g)
    implicit none
    ! I/O
    integer                            :: ps            !! Selected precalculated dataset
    type(fbem_bem_element)             :: e             !! Element
    logical                            :: reverse       !! Reverse normal vector
    real(kind=real64)                  :: x_i(2)        !! Collocation point position vector
    type(fbem_bem_harpot2d_parameters) :: pars          !! Parameters of the region
    complex(kind=real64)               :: h(e%n_pnodes) !! h integration kernels matrix
    complex(kind=real64)               :: g(e%n_snodes) !! g integration kernels matrix
    ! Local
    integer              :: kip                  ! Counter variable for integration points loop
    real(kind=real64)    :: x(2)                 ! Position vector at integration point
    real(kind=real64)    :: n(2)                 ! Unit normal vector at integration point
    real(kind=real64)    :: pphijw(e%n_pnodes)   ! phi^p * jacobian * weight at integration point
    real(kind=real64)    :: sphijw(e%n_snodes)   ! phi^s * jacobian * weight at integration point
    real(kind=real64)    :: rv(2)                ! Radiovector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, d1r, logr         ! Radiovector module, squared, inverses and log(r)
    real(kind=real64)    :: drdx(2)              ! Radiovector derivatives with respect to x_k
    real(kind=real64)    :: drdn                 ! Partial derivative of r respect to unit normal
    complex(kind=real64) :: z(1)                 ! Bessel functions arguments: z=ikr
    complex(kind=real64) :: KnR(0:2,1)           ! Bessel functions decomposition
    complex(kind=real64) :: P, Q                 ! Components of the fundamental solutions
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
      d1r=1.0d0/r
      logr=log(r)
      drdx=rv*d1r
      drdn=dot_product(drdx,n)
      z(1)=c_im*pars%k*r
      call fbem_BesselKnR_decomposed(1,z,KnR)
      P=pars%P(1)*logr+pars%P(2)+pars%P(3)*KnR(0,1)
      Q=pars%Q(1)*d1r+(pars%Q(2)*logr+pars%Q(3))*r+pars%Q(4)*KnR(1,1)
      h=h+Q*drdn*pphijw
      g=g+P*sphijw
    end do
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_harpot2d_sbie_ext_pre

  subroutine fbem_bem_harpot2d_sbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,pars,gln,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e             !! Integration element
    logical                            :: reverse       !! Reverse normal vector
    real(kind=real64)                  :: xi_s(1,2)     !! Coordinates of the subdivision of the parent element (xi space [-1,1])
    real(kind=real64)                  :: x_i(2)        !! Collocation point position vector
    real(kind=real64)                  :: barxip(1)     !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)                  :: barr          !! Telles jacobian at barxip
    type(fbem_bem_harpot2d_parameters) :: pars          !! Parameters of the region
    integer                            :: gln           !! 1D Gauss-Legendre number of integration points (<=32)
    complex(kind=real64)               :: h(e%n_pnodes) !! h kernel vector
    complex(kind=real64)               :: g(e%n_snodes) !! g kernel vector
    ! Local
    integer                      :: kphi                 ! Counter variable for shape functions loops
    integer                      :: kip                  ! Counter variable of integration points
    real(kind=real64)            :: gamma                ! Coordinate gamma (Telles transformation space [0,1])
    real(kind=real64)            :: w                    ! Weights of an integration point
    type(fbem_telles_parameters) :: telles_parameters    ! Telles parameters
    real(kind=real64)            :: jt                   ! Telles jacobian: xip->gamma
    real(kind=real64)            :: xip                  ! Coordinate xip (subdivision space [0,1])
    real(kind=real64)            :: js                   ! Subdivision jacobian: xi->xip
    real(kind=real64)            :: xi                   ! Coordinate xi [xi1,xi2]
    real(kind=real64)            :: aux(10)              ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)     ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi(e%n_gnodes) ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: pphi(e%n_pnodes)     ! Functional shape functions values
    real(kind=real64)            :: sphi(e%n_snodes)     ! Functional shape functions values
    real(kind=real64)            :: x(2)                 ! Position vector at xi
    real(kind=real64)            :: T(2)                 ! Tangent vector at xi
    real(kind=real64)            :: N(2)                 ! Normal vector at xi
    real(kind=real64)            :: rv(2)                ! Radiovector between collocation point and integration point (x-x_i)
    real(kind=real64)            :: r, d1r, logr         ! Radiovector module, squared, its inverses and log(r)
    real(kind=real64)            :: drdx(2)              ! Radiovector derivatives with respect to x_k
    real(kind=real64)            :: jg                   ! Geometric jacobian
    real(kind=real64)            :: drdn                 ! Partial derivative of r respect to unit normal
    real(kind=real64)            :: jw                   ! Jacobian * weight
    real(kind=real64)            :: pphijw(e%n_pnodes)   ! Auxiliary variables for integrand evaluation
    real(kind=real64)            :: sphijw(e%n_snodes)   ! Auxiliary variables for integrand evaluation
    complex(kind=real64)         :: z(1)                 ! Bessel functions arguments: z=ikr
    complex(kind=real64)         :: KnR(0:2,1)           ! Bessel functions decomposition
    complex(kind=real64)         :: P, Q                 ! Components of the fundamental solutions
    ! Initialization
    h=(0.d0,0.d0)
    g=(0.d0,0.d0)
    telles_parameters=fbem_telles11_calculate_parameters(barxip(1),barr)
    ! Numerical integration
    do kip=1,gl11_n(gln)
      ! GAMMA COORDINATE
      gamma=gl11_xi(kip,gln)
      w=gl11_w(kip,gln)
      ! GAMMA->XIP TRANSFORMATION
      call fbem_telles_xi_and_jacobian(telles_parameters,gamma,xip,jt)
      ! XIP->XI TRANSFORMATION
      xi=0.5d0*(1.d0-xip)*xi_s(1,1)+0.5d0*(1.d0+xip)*xi_s(1,2)
      js=0.5d0*(xi_s(1,2)-xi_s(1,1))
      ! XI->X TRANSFORMATION
#     define etype e%gtype
#     define delta 0.d0
#     define phi gphi
#     define dphidxi dgphidxi
#     include <phi_and_dphidxi_1d.rc>
#     undef etype
#     undef delta
#     undef phi
#     undef dphidxi
      x=0.d0
      T=0.d0
      do kphi=1,e%n_gnodes
        x=x+gphi(kphi)*e%x(:,kphi)
        T=T+dgphidxi(kphi)*e%x(:,kphi)
      end do
      ! Normal vector
      N(1)=T(2)
      N(2)=-T(1)
      ! Geometric jacobian
      jg=sqrt(dot_product(T,T))
      ! Unit normal vector
      n=N/jg
      ! Distance vector
      rv=x-x_i
      ! Radiovector norm and its inverse
      r=sqrt(dot_product(rv,rv))
      d1r=1.0d0/r
      logr=log(r)
      drdx=rv*d1r
      drdn=dot_product(drdx,n)
      ! Jacobian * weight
      jw=jg*js*jt*w
      ! FUNCTIONAL SHAPE FUNCTIONS
      ! Functional shape functions (primary variables) at xi
#     define etype e%ptype
#     define delta e%ptype_delta
#     define phi pphi
#     include <phi_1d.rc>
#     undef etype
#     undef delta
#     undef phi
      ! Functional shape functions (secondary variables) at xi
#     define etype e%stype
#     define delta e%stype_delta
#     define phi sphi
#     include <phi_1d.rc>
#     undef etype
#     undef delta
#     undef phi
      ! Functional shape functions * jacobians* weights
      pphijw=pphi*jw
      sphijw=sphi*jw
      ! COMPONENTS OF THE FUNDAMENTAL SOLUTION
      z(1)=c_im*pars%k*r
      call fbem_BesselKnR_decomposed(1,z,KnR)
      ! Components of the fundamental solutions
      P=pars%P(1)*logr+pars%P(2)+pars%P(3)*KnR(0,1)
      Q=pars%Q(1)*d1r+(pars%Q(2)*logr+pars%Q(3))*r+pars%Q(4)*KnR(1,1)
      ! Add
      h=h+Q*drdn*pphijw
      g=g+P*sphijw
    end do
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_harpot2d_sbie_ext_st

  recursive subroutine fbem_bem_harpot2d_sbie_ext_adp(e,reverse,xi_s,x_i,p,qsp,ks,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e             !! Element
    logical                            :: reverse       !! Reverse orientation
    real(kind=real64)                  :: xi_s(1,2)     !! Subdivision of the parent element
    real(kind=real64)                  :: x_i(2)        !! Collocation point position vector
    type(fbem_bem_harpot2d_parameters) :: p             !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp           !! Quasi-singular integration parameters
    integer                            :: ks            !! Current level of subdivisions
    integer                            :: ns            !! Maximum level of subdivision
    complex(kind=real64)               :: h(e%n_pnodes) !! h integration kernels matrix
    complex(kind=real64)               :: g(e%n_snodes) !! g integration kernels matrix
    ! Local
    integer              :: gln_near          ! 1D Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer              :: gln               ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical              :: subdivide         ! True if subdivision has to be performed
    real(kind=real64)    :: barxi(1)          ! Nearest element coordinate with respect to collocation point
    real(kind=real64)    :: barxip(1)         ! Nearest element subdivision local coordinate with respect to collocation point
    real(kind=real64)    :: rmin              ! Minimum distance between collocation point and barxi on the element
    real(kind=real64)    :: barr              ! Telles jacobian at barxi
    real(kind=real64)    :: cl                ! Characteristic length
    real(kind=real64)    :: d                 ! Normalized distance between collocation point and element subdivision
    integer              :: method            ! Method used in nearest point algorithm
    real(kind=real64)    :: tmp_xi_s(1,2)     ! Subdivision
    real(kind=real64)    :: x_s(2,e%n_gnodes) ! Coordinates of the element subdivision
    complex(kind=real64) :: h_tmp(e%n_pnodes) ! h integration kernels matrix (temporary)
    complex(kind=real64) :: g_tmp(e%n_snodes) ! g integration kernels matrix (temporary)
    ! Initialize
    if (ks.eq.1) then
      h=(0.d0,0.d0)
      g=(0.d0,0.d0)
      xi_s(1,1)=-1.d0
      xi_s(1,2)= 1.d0
      call fbem_nearest_element_point_bem(2,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
      barxip=barxi
    else
      call fbem_obtain_element_subdivision_coordinates(2,e%gtype,e%x,xi_s,x_s)
      cl=fbem_characteristic_length(2,e%gtype,x_s,1.d-12)
      call fbem_nearest_element_point_bem(2,e%gtype,x_s,cl,x_i,barxip,rmin,d,method)
    end if
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,2,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_harpot2d_sbie_ext_adp',ns,'maximum number of subdivisions reached')
        gln_near=30
      end if
    else
      if (gln_near.eq.0) subdivide=.true.
    end if
    ! Subdivide by 1/2
    if (subdivide) then
      ! SUBLINE 1
      tmp_xi_s(1,1)=xi_s(1,1)
      tmp_xi_s(1,2)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      call fbem_bem_harpot2d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_harpot2d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_harpot2d_sbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,p,gln,h_tmp,g_tmp)
      h=h+h_tmp
      g=g+g_tmp
    end if
  end subroutine fbem_bem_harpot2d_sbie_ext_adp

  subroutine fbem_bem_harpot2d_sbie_int(ngp,type_g,type_f1,type_f2,delta_f,x_nodes,reverse,xi_i,pars,h,g)
    implicit none
    ! I/O
    integer                            :: ngp                             !! Number of Gauss point to be used
    integer                            :: type_g                          !! Geometrial interpolation
    integer                            :: type_f1                         !! Functional interpolation (primary variables)
    integer                            :: type_f2                         !! Functional interpolation (secondary variables)
    real(kind=real64)                  :: delta_f                         !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64)                  :: x_nodes(2,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
    logical                            :: reverse                         !! Normal vector inversion
    real(kind=real64)                  :: xi_i                            !! Reference coordinate of the singular point.
    type(fbem_bem_harpot2d_parameters) :: pars                            !! Parameters of the region
    complex(kind=real64)               :: h(fbem_n_nodes(type_f1))        !! h kernel vector
    complex(kind=real64)               :: g(fbem_n_nodes(type_f2))        !! g kernel vector
    ! Local
    integer                            :: nnodes_g                        ! Number of nodes of the element
    integer                            :: kphi                            ! Counter variable for shape functions loops
    integer                            :: kip                             ! Counter of gamma coordinates
    real(kind=real64)                  :: x_i(2)                          ! Real coordinates of collocation point
    real(kind=real64)                  :: phi_f1(fbem_n_nodes(type_f1))   ! Functional shape functions values at xi
    real(kind=real64)                  :: phi_f2(fbem_n_nodes(type_f2))   ! Functional shape functions values at xi
    real(kind=real64)                  :: phi_g(fbem_n_nodes(type_g))     ! Geometrical shape functions values at xi
    real(kind=real64)                  :: dphidxi_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi
    integer                            :: nsub                            ! Number of subdivision of the element
    integer                            :: ksub                            ! Counter of subdivision
    real(kind=real64)                  :: gamma                           ! Reference coordinate
    real(kind=real64)                  :: w                               ! Weights of each integration point
    type(fbem_telles_parameters)       :: telles_parameters               ! Telles parameters
    real(kind=real64)                  :: jt                              ! Telles jacobian
    real(kind=real64)                  :: xip                             ! Reference coordinate in subdivided element [0,1]
    real(kind=real64)                  :: js                              ! Jacobian of the xi [xisub(1,:),xisub[2,:]] -> xip [0,1] transformation
    real(kind=real64)                  :: xi                              ! Reference coordinate
    real(kind=real64)                  :: xisub(2,2)                      ! Coordinates of subdivisions in xi space
    real(kind=real64)                  :: xip_i(2)                        ! Singular point in xi_sub space
    real(kind=real64)                  :: aux(10)                         ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)                  :: x(2)                            ! Position vector at xi
    real(kind=real64)                  :: T(2)                            ! Tangent vector at xi
    real(kind=real64)                  :: N(2)                            ! Normal vector at xi
    real(kind=real64)                  :: rv(2)                           ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)                  :: r, d1r, logr                    ! Distance vector module
    real(kind=real64)                  :: jg                              ! Geometric jacobian
    real(kind=real64)                  :: drdn                            ! Partial derivative of r respect to unit normal
    real(kind=real64)                  :: jw                              ! Jacobians * weight
    real(kind=real64)                  :: phif1jw(fbem_n_nodes(type_f1))  ! Auxiliary variables for integrand evaluation
    real(kind=real64)                  :: phif2jw(fbem_n_nodes(type_f2))  ! Auxiliary variables for integrand evaluation
    complex(kind=real64)               :: z                               ! Argument of Modified Bessel functions
    complex(kind=real64)               :: K0r, K1r, K2r                   ! Modified Bessel functions K values
    complex(kind=real64)               :: P, Q                            ! Components of the fundamental solutions
    !
    ! Initialize
    !
    ! Initialize kernel vectors
    h=(0.d0,0.d0)
    g=(0.d0,0.d0)
    ! Number of nodes of the element
    nnodes_g=fbem_n_nodes(type_g)
    ! Calculate real coordinates of collocation point
    ! Geometrical shape functions at xi_i
#   define etype type_g
#   define delta 0.0d0
#   define xi xi_i
#   define phi phi_g
#   include <phi_1d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
    x_i=0.d0
    do kphi=1,nnodes_g
      x_i=x_i+phi_g(kphi)*x_nodes(:,kphi)
    end do
    ! If xi_i belongs to one of the vertices, subdivision is not needed
    if (fbem_check_xi_vertex(xi_i).eqv.(.true.)) then
      nsub=1
    ! If xi_i is inside the element, 2 subdivisions are needed
    else
      nsub=2
    end if
    ! Setup the subdivisions transformation to [0,1] domain
    select case (nsub)
      ! If just 1 subdivision
      case (1)
       ! The subdivision is the whole element
        xisub(1,1)=-1.0d0
        xisub(2,1)=1.0d0
        ! If xi_i is at xi=-1
        if (xi_i.lt.0.0d0) xip_i(1)=0.0d0
        ! If xi_i is at xi=1
        if (xi_i.gt.0.0d0) xip_i(1)=1.0d0
      ! If 2 subdivisions
      case (2)
        ! Subdivision 1
        xisub(1,1)=-1.0d0
        xisub(2,1)=xi_i
        ! The xip coordinate of the collocation point
        xip_i(1)=1.0d0
        ! Subdivision 2
        xisub(1,2)=xi_i
        xisub(2,2)=1.0d0
        ! The xip coordinate of the collocation point
        xip_i(2)=0.0d0
    end select
    !
    ! Integrate weakly singular integrals
    !
    ! Loop through subdivisions
    do ksub=1,nsub
      ! Subdivision jacobian
      js=xisub(2,ksub)-xisub(1,ksub)
      ! Telles transformation parameters
      telles_parameters=fbem_telles01_calculate_parameters(xip_i(ksub),0.0d0)
      ! Loop through integration points
      do kip=1,gl01_n(ngp)
        ! GAMMA->XIP->XI COORDINATE TRANSFORMATION
        ! Gamma coordinate and weight
        gamma=gl01_xi(kip,ngp)
        w=gl01_w(kip,ngp)
        ! xip coordinate, weight and jacobian from Telles transformation
        call fbem_telles_xi_and_jacobian(telles_parameters,gamma,xip,jt)
        ! xi coordinate
        xi=js*xip+xisub(1,ksub)
        ! COMPONENTS OF THE FUNDAMENTAL SOLUTIONS
        ! Geometrical shape functions and first derivatives at xi
#       define etype type_g
#       define delta 0.0d0
#       define phi phi_g
#       define dphidxi dphidxi_g
#       include <phi_and_dphidxi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
#       undef dphidxi
        ! Components calculation of x and T at xi
        x=0.d0
        T=0.d0
        do kphi=1,nnodes_g
          x=x+phi_g(kphi)*x_nodes(:,kphi)
          T=T+dphidxi_g(kphi)*x_nodes(:,kphi)
        end do
        ! Normal vector
        N(1)=T(2)
        N(2)=-T(1)
        ! Geometric jacobian
        jg=dot_product(T,T)
        jg=sqrt(jg)
        ! Distance vector
        rv=x-x_i
        ! Distance vector norm
        r=dot_product(rv,rv)
        r=sqrt(r)
        d1r=1.0d0/r
        logr=log(r)
        ! dr/dn
        drdn=dot_product(rv,N)*d1r/jg
        ! Components of the fundamental solution
        P=pars%P(1)*logr
        ! FUNCTIONAL SHAPE FUNCTIONS
        ! Functional shape functions (secondary variables) at xi
#       define etype type_f2
#       define delta delta_f
#       define phi phi_f2
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Jacobians * weight
        jw=jg*js*jt*w
        ! Functional shape functions * jacobians* weights
        phif2jw=phi_f2*jw
        ! Add to kernels
        g=g+P*phif2jw
      end do ! Loop through integration points
    end do ! Loop through subdivisions
    !
    ! Integrate regular integrals
    !
    ! Loop through subdivisions
    do ksub=1,nsub
      ! Subdivision jacobian
      js=xisub(2,ksub)-xisub(1,ksub)
      ! Loop through integration points
      do kip=1,gl01_n(ngp)
        ! XIP->XI COORDINATE TRANSFORMATION
        ! Coordinate and weight in gamma [0,1]
        xip=gl01_xi(kip,ngp)
        w=gl01_w(kip,ngp)
        ! Calculate xi
        xi=js*xip+xisub(1,ksub)
        ! COMPONENTS OF THE FUNDAMENTAL SOLUTIONS
        ! Geometrical shape functions and first derivatives at xi
#       define etype type_g
#       define delta 0.0d0
#       define phi phi_g
#       define dphidxi dphidxi_g
#       include <phi_and_dphidxi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
#       undef dphidxi
        ! Components calculation of x and T at xi
        x=0.d0
        T=0.d0
        do kphi=1,nnodes_g
          x=x+phi_g(kphi)*x_nodes(:,kphi)
          T=T+dphidxi_g(kphi)*x_nodes(:,kphi)
        end do
        ! Normal vector
        N(1)=T(2)
        N(2)=-T(1)
        ! Geometric jacobian
        jg=dot_product(T,T)
        jg=sqrt(jg)
        ! Distance vector
        rv=x-x_i
        ! Distance vector norm
        r=dot_product(rv,rv)
        r=sqrt(r)
        d1r=1.0d0/r
        logr=log(r)
        ! dr/dn
        drdn=dot_product(rv,N)*d1r/jg
        ! Modified Bessel functions argument, z=ikr
        z=c_im*pars%k*r
        call fbem_modified_bessel_K0r_K1r_K2r_2(z,K0r,K1r,K2r)
        ! Components of the fundamental solutions
        P=pars%P(2)+pars%P(3)*K0r
        Q=pars%Q(1)*d1r+(pars%Q(2)*logr+pars%Q(3))*r+pars%Q(4)*K1r
        ! FUNCTIONAL SHAPE FUNCTIONS
        ! Functional shape functions (primary variables) at xi
#       define etype type_f1
#       define delta delta_f
#       define phi phi_f1
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Functional shape functions (secondary variables) at xi
#       define etype type_f2
#       define delta delta_f
#       define phi phi_f2
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Jacobians * weight
        jw=jg*js*w
        ! Functional shape functions * jacobians* weights
        phif1jw=phi_f1*jw
        phif2jw=phi_f2*jw
        ! Add to kernels
        h=h+Q*drdn*phif1jw
        g=g+P*phif2jw
      end do ! Loop through integration points
    end do ! Loop through subdivisions
    ! If the normal is inverted, then h=-h
    if (reverse) h=-h
  end subroutine fbem_bem_harpot2d_sbie_int

  subroutine fbem_bem_harpot2d_sbie_auto(e,reverse,x_i,p,qsp,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e             !! Integration element
    logical                            :: reverse       !! Reverse orientation
    real(kind=real64)                  :: x_i(2)        !! Collocation point
    type(fbem_bem_harpot2d_parameters) :: p             !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp           !! Quasi-singular integration parameters
    integer                            :: ns            !! Maximum level of subdivisions
    complex(kind=real64)               :: h(e%n_pnodes) !! h integration kernel
    complex(kind=real64)               :: g(e%n_snodes) !! g integration kernel
    ! Local
    real(kind=real64) :: r(2)      ! Distance vector
    real(kind=real64) :: rmin      ! Minimum distance between element and x_i
    real(kind=real64) :: barxi(1)  ! Local coordinates of the nearest element point with respect to x_i
    real(kind=real64) :: d         ! Dimensionless distance
    integer           :: delta     ! Control variable
    real(kind=real64) :: xi_s(1,2) ! Local coordinates of the element subdivision
    integer           :: method    ! Method used when calculating the nearest element point
    integer           :: gln_near  ! 1D Gauss-Legendre integration points required by the quasi-singular function
    integer           :: gln       ! 1D Gauss-Legendre integration points used in the integration
    integer           :: ps        ! Selected precalculated dataset
    integer           :: i         ! Counter
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
      call fbem_nearest_element_point_bem(2,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
      if (d.le.1.d-12) then
        delta=1
      else
        delta=0
      end if
    end if
    ! Integrate
    select case (delta)
      case (1)
        call fbem_bem_harpot2d_sbie_int(30,e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi(1),p,h,g)
      case (0)
        ! Estimate the required integration rule
        gln_near=fbem_qs_n_estimation_standard(e%n,e%gtype,2,qsp,d,barxi)
        gln=max(e%gln_far,gln_near)
        ! Integrate using a conservative precalculated dataset
        if ((gln.le.e%ps_gln_max).and.(gln_near.gt.0)) then
          do i=1,e%n_ps
            if (e%ps_gln(i).ge.gln) then
              ps=i
              exit
            end if
          end do
          call fbem_bem_harpot2d_sbie_ext_pre(ps,e,reverse,x_i,p,h,g)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_harpot2d_sbie_ext_adp(e,reverse,xi_s,x_i,p,qsp,1,ns,h,g)
        end if
    end select
  end subroutine fbem_bem_harpot2d_sbie_auto

  !
  ! BODY LOADS
  !





  subroutine fbem_bem_harpot2d_sbie_bl_auto(e,x_i,p,qsp,ns,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e             !! Integration element
    real(kind=real64)                  :: x_i(2)        !! Collocation point
    type(fbem_bem_harpot2d_parameters) :: p             !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp           !! Quasi-singular integration parameters
    integer                            :: ns            !! Maximum level of subdivisions
    complex(kind=real64)               :: g(e%n_snodes) !! g integration kernel
    ! Local
    real(kind=real64) :: r(2)      ! Distance vector
    real(kind=real64) :: rmin      ! Minimum distance between element and x_i
    real(kind=real64) :: barxi(1)  ! Local coordinates of the nearest element point with respect to x_i
    real(kind=real64) :: d         ! Dimensionless distance
    integer           :: delta     ! Control variable
    real(kind=real64) :: xi_s(1,2) ! Local coordinates of the element subdivision
    integer           :: method    ! Method used when calculating the nearest element point
    integer           :: gln_near  ! 1D Gauss-Legendre integration points required by the quasi-singular function
    integer           :: gln       ! 1D Gauss-Legendre integration points used in the integration
    integer           :: ps        ! Selected precalculated dataset
    integer           :: i         ! Counter
    ! POINT BODY LOAD
    if (e%d.eq.0) then
      r=e%x(:,1)-x_i
      rmin=sqrt(dot_product(r,r))
      if (rmin.eq.0.d0) then
        call fbem_error_message(output_unit,0,'fbem_bem_harpot2d_sbie_bl_auto',0,'it is not possible to collocate at a point load')
      else
        call fbem_bem_harpot2d_sbie_p(e%x(:,1),x_i,p,g(1))
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
      call fbem_error_message(output_unit,0,'fbem_bem_harpot2d_sbie_bl_auto',0,'only point loads are available')
      r=e%bball_centre-x_i
      rmin=sqrt(dot_product(r,r))-e%bball_radius
      if (rmin.gt.(4.d0*e%bball_radius)) then
        delta=0
        barxi=0.d0
        d=rmin/e%cl
      else
        ! Use an adaptative algorithm that combines sampling and minimization algorithms
        call fbem_nearest_element_point_bem(2,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
        if (d.le.1.d-12) then
          delta=1
        else
          delta=0
        end if
      end if
      ! Integrate
      select case (delta)
        case (1)
          !call fbem_bem_harpot2d_sbie_bl_int(30,e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,barxi(1),p,g)
        case (0)
          ! Estimate the required integration rule
          gln_near=fbem_qs_n_estimation_standard(e%n,e%gtype,2,qsp,d,barxi)
          gln=max(e%gln_far,gln_near)
          ! Integrate using a conservative precalculated dataset
          if ((gln.le.e%ps_gln_max).and.(gln_near.gt.0)) then
            do i=1,e%n_ps
              if (e%ps_gln(i).ge.gln) then
                ps=i
                exit
              end if
            end do
            !call fbem_bem_harpot2d_sbie_bl_ext_pre(ps,e,x_i,p,g)
          ! Integrate using an adaptative algorithm
          else
            !call fbem_bem_harpot2d_sbie_bl_ext_adp(e,xi_s,x_i,p,qsp,1,ns,g)
          end if
      end select
    end if
  end subroutine fbem_bem_harpot2d_sbie_bl_auto






  ! ================================================================================================================================

  ! ================================================================================================================================
  ! HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (HBIE)

  !! Fundamental solution d*
  subroutine fbem_bem_harpot2d_hbie_d(x,x_i,n_i,pars,do)
    implicit none
    ! I/O
    real(kind=real64)                  :: x(2)    !! Observation point
    real(kind=real64)                  :: x_i(2)  !! Collocation point
    real(kind=real64)                  :: n_i(2)  !! Collocation point unit normal
    type(fbem_bem_harpot2d_parameters) :: pars    !! Parameters of the region
    complex(kind=real64)               :: do      !! d*
    ! Local
    real(kind=real64)    :: rv(2), r, d1r1, d1r2, logr, drdx(2), drdni
    complex(kind=real64) :: z(1), KnR(0:2,1)
    rv=x-x_i
    r=sqrt(dot_product(rv,rv))
    d1r1=1.0d0/r
    d1r2=d1r1*d1r1
    logr=log(r)
    drdx=rv*d1r1
    drdni=-dot_product(drdx,n_i)
    z(1)=c_im*pars%k*r
    call fbem_BesselKnR_decomposed(1,z,KnR)
    do=(pars%Q(1)*d1r1+(pars%Q(2)*logr+pars%Q(3))*r+pars%Q(4)*KnR(1,1))*drdni
  end subroutine fbem_bem_harpot2d_hbie_d

  !! Fundamental solution s*
  subroutine fbem_bem_harpot2d_hbie_s(x,n,x_i,n_i,pars,so)
    implicit none
    ! I/O
    real(kind=real64)                  :: x(2)    !! Observation point
    real(kind=real64)                  :: n(2)    !! Observation point unit normal
    real(kind=real64)                  :: x_i(2)  !! Collocation point
    real(kind=real64)                  :: n_i(2)  !! Collocation point unit normal
    type(fbem_bem_harpot2d_parameters) :: pars    !! Parameters of the region
    complex(kind=real64)               :: so      !! s*
    ! Local
    real(kind=real64)    :: rv(2), r, d1r1, d1r2, logr, drdx(2), drdn, drdni, n_dot_ni
    complex(kind=real64) :: z(1), KnR(0:2,1)
    complex(kind=real64) :: S1, S2
    rv=x-x_i
    r=sqrt(dot_product(rv,rv))
    d1r1=1.0d0/r
    d1r2=d1r1*d1r1
    logr=log(r)
    drdx=rv*d1r1
    drdn=dot_product(drdx,n)
    drdni=-dot_product(drdx,n_i)
    n_dot_ni=dot_product(n,n_i)
    z(1)=c_im*pars%k*r
    call fbem_BesselKnR_decomposed(1,z,KnR)
    S1=pars%S1(1)*d1r2+pars%S1(2)+(pars%S1(3)*logr+pars%S1(4))*r**2+pars%S1(5)*KnR(2,1)
    S2=pars%S2(1)*d1r2+pars%S2(2)*logr+pars%S2(3)+pars%S2(4)*d1r1*KnR(1,1)
    so=S1*drdn*drdni+S2*n_dot_ni
  end subroutine fbem_bem_harpot2d_hbie_s

  subroutine fbem_bem_harpot2d_hbie_ext_pre(ps,e,reverse,x_i,n_i,pars,m,l)
    implicit none
    ! I/O
    integer                            :: ps            !! Selected precalculated dataset
    type(fbem_bem_element)             :: e             !! Element
    logical                            :: reverse       !! Reverse normal vector
    real(kind=real64)                  :: x_i(2)        !! Collocation point position vector
    real(kind=real64)                  :: n_i(2)        !! Collocation point unit normal vector
    type(fbem_bem_harpot2d_parameters) :: pars          !! Parameters of the region
    complex(kind=real64)               :: m(e%n_pnodes) !! m kernels matrix
    complex(kind=real64)               :: l(e%n_snodes) !! l kernels matrix
    ! Local
    integer              :: kip                            ! Counter variable for integration points loop
    real(kind=real64)    :: x(2)                           ! Position vector at integration point
    real(kind=real64)    :: n(2)                           ! Unit normal vector at integration point
    real(kind=real64)    :: pphijw(e%n_pnodes)             ! phi^p * jacobian * weight at integration point
    real(kind=real64)    :: sphijw(e%n_snodes)             ! phi^s * jacobian * weight at integration point
    real(kind=real64)    :: rv(2)                          ! Radiovector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, d1r1, d1r2, logr            ! Radiovector module, squared, inverses and log(r)
    real(kind=real64)    :: drdx(2)                        ! Radiovector derivatives with respect to x_k
    real(kind=real64)    :: drdn                           ! Partial derivative of r respect to unit normal
    real(kind=real64)    :: drdni                          ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)    :: n_dot_ni                       ! Dot product of unit normals
    complex(kind=real64) :: z(1)                           ! Bessel functions arguments: z=ikr
    complex(kind=real64) :: KnR(0:2,1)                     ! Bessel functions decomposition
    complex(kind=real64) :: Q, S1, S2                      ! Components of the fundamental solution
    complex(kind=real64) :: fs_d, fs_s                     ! Fundamental solution values
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
      d1r1=1.0d0/r
      d1r2=d1r1*d1r1
      logr=log(r)
      drdx=rv*d1r1
      drdn=dot_product(drdx,n)
      drdni=-dot_product(drdx,n_i)
      n_dot_ni=dot_product(n,n_i)
      z(1)=c_im*pars%k*r
      call fbem_BesselKnR_decomposed(1,z,KnR)
      Q=pars%Q(1)*d1r1+(pars%Q(2)*logr+pars%Q(3))*r+pars%Q(4)*KnR(1,1)
      S1=pars%S1(1)*d1r2+pars%S1(2)+(pars%S1(3)*logr+pars%S1(4))*r**2+pars%S1(5)*KnR(2,1)
      S2=pars%S2(1)*d1r2+pars%S2(2)*logr+pars%S2(3)+pars%S2(4)*d1r1*KnR(1,1)
      fs_d=Q*drdni
      fs_s=S1*drdn*drdni+S2*n_dot_ni
      m=m+fs_s*pphijw
      l=l+fs_d*sphijw
    end do
    ! Reverse if needed
    if (reverse) m=-m
  end subroutine fbem_bem_harpot2d_hbie_ext_pre

  subroutine fbem_bem_harpot2d_hbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,pars,gln,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e             !! Integration element
    logical                            :: reverse       !! Reverse normal vector
    real(kind=real64)                  :: xi_s(1,2)     !! Coordinates of the subdivision of the parent element (xi space [-1,1])
    real(kind=real64)                  :: x_i(2)        !! Collocation point position vector
    real(kind=real64)                  :: n_i(2)        !! Collocation point unit normal vector
    real(kind=real64)                  :: barxip(1)     !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)                  :: barr          !! Telles jacobian at barxip
    type(fbem_bem_harpot2d_parameters) :: pars          !! Parameters of the region
    integer                            :: gln           !! 1D Gauss-Legendre number of integration points (<=32)
    complex(kind=real64)               :: m(e%n_pnodes) !! m kernels matrix
    complex(kind=real64)               :: l(e%n_snodes) !! l kernels matrix
    ! Local
    integer                      :: kphi                    ! Counter variable for shape functions loops
    integer                      :: kip                     ! Counter variable of integration points
    real(kind=real64)            :: gamma                   ! Coordinate gamma (Telles transformation space [0,1])
    real(kind=real64)            :: w                       ! Weights of an integration point
    type(fbem_telles_parameters) :: telles_parameters       ! Telles parameters
    real(kind=real64)            :: jt                      ! Telles jacobian: xip->gamma
    real(kind=real64)            :: xip                     ! Coordinate xip (subdivision space [0,1])
    real(kind=real64)            :: js                      ! Subdivision jacobian: xi->xip
    real(kind=real64)            :: xi                      ! Coordinate xi [xi1,xi2]
    real(kind=real64)            :: aux(10)                 ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)        ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi(e%n_gnodes)    ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: pphi(e%n_pnodes)        ! Functional shape functions values
    real(kind=real64)            :: sphi(e%n_snodes)        ! Functional shape functions values
    real(kind=real64)            :: x(2)                    ! Position vector at xi
    real(kind=real64)            :: T(2)                    ! Tangent vector at xi
    real(kind=real64)            :: N(2)                    ! Normal vector at xi
    real(kind=real64)            :: rv(2)                   ! Radiovector between collocation point and integration point (x-x_i)
    real(kind=real64)            :: r, d1r1, d1r2, logr     ! Radiovector module, squared, inverses and log(r)
    real(kind=real64)            :: drdx(2)                 ! Radiovector derivatives with respect to x_k
    real(kind=real64)            :: jg                      ! Geometric jacobian
    real(kind=real64)            :: drdn                    ! Partial derivative of r respect to unit normal
    real(kind=real64)            :: drdni                   ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)            :: n_dot_ni                ! Dot product of unit normals
    real(kind=real64)            :: jw                      ! Jacobians * weights
    real(kind=real64)            :: pphijw(e%n_pnodes)      ! Auxiliary variables for integrand evaluation
    real(kind=real64)            :: sphijw(e%n_snodes)      ! Auxiliary variables for integrand evaluation
    complex(kind=real64)         :: z(1)                    ! Bessel functions arguments: z=ikr
    complex(kind=real64)         :: KnR(0:2,1)              ! Bessel functions decomposition
    complex(kind=real64)         :: Q, S1, S2               ! Components of the fundamental solution
    complex(kind=real64)         :: fs_d, fs_s              ! Fundamental solution values
    ! Initialize
    m=(0.0d0,0.0d0)
    l=(0.0d0,0.0d0)
    telles_parameters=fbem_telles11_calculate_parameters(barxip(1),barr)
    ! Numerical integration
    do kip=1,gl11_n(gln)
      ! GAMMA COORDINATE
      gamma=gl11_xi(kip,gln)
      w=gl11_w(kip,gln)
      ! GAMMA->XIP TRANSFORMATION
      call fbem_telles_xi_and_jacobian(telles_parameters,gamma,xip,jt)
      ! XIP->XI TRANSFORMATION
      xi=0.5d0*(1.d0-xip)*xi_s(1,1)+0.5d0*(1.d0+xip)*xi_s(1,2)
      js=0.5d0*(xi_s(1,2)-xi_s(1,1))
      ! XI->X TRANSFORMATION
#     define etype e%gtype
#     define delta 0.d0
#     define phi gphi
#     define dphidxi dgphidxi
#     include <phi_and_dphidxi_1d.rc>
#     undef etype
#     undef delta
#     undef phi
#     undef dphidxi
      x=0.d0
      T=0.d0
      do kphi=1,e%n_gnodes
        x=x+gphi(kphi)*e%x(:,kphi)
        T=T+dgphidxi(kphi)*e%x(:,kphi)
      end do
      ! Normal vector
      N(1)=T(2)
      N(2)=-T(1)
      ! Geometric jacobian
      jg=sqrt(dot_product(T,T))
      ! Unit normal vector
      n=N/jg
      ! Radiovector
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      d1r1=1.0d0/r
      d1r2=d1r1*d1r1
      logr=log(r)
      drdx=rv*d1r1
      drdn=dot_product(drdx,n)
      drdni=-dot_product(drdx,n_i)
      n_dot_ni=dot_product(n,n_i)
      ! Jacobian * weight
      jw=jg*js*jt*w
      ! FUNCTIONAL SHAPE FUNCTIONS
      ! Functional shape functions (primary variables) at xi
#     define etype e%ptype
#     define delta e%ptype_delta
#     define phi pphi
#     include <phi_1d.rc>
#     undef etype
#     undef delta
#     undef phi
      ! Functional shape functions (secondary variables) at xi
#     define etype e%stype
#     define delta e%stype_delta
#     define phi sphi
#     include <phi_1d.rc>
#     undef etype
#     undef delta
#     undef phi
      ! Functional shape functions * jacobians* weights
      pphijw=pphi*jw
      sphijw=sphi*jw
      ! COMPONENTS OF THE FUNDAMENTAL SOLUTION
      z(1)=c_im*pars%k*r
      call fbem_BesselKnR_decomposed(1,z,KnR)
      Q=pars%Q(1)*d1r1+(pars%Q(2)*logr+pars%Q(3))*r+pars%Q(4)*KnR(1,1)
      S1=pars%S1(1)*d1r2+pars%S1(2)+(pars%S1(3)*logr+pars%S1(4))*r**2+pars%S1(5)*KnR(2,1)
      S2=pars%S2(1)*d1r2+pars%S2(2)*logr+pars%S2(3)+pars%S2(4)*d1r1*KnR(1,1)
      fs_d=Q*drdni
      fs_s=S1*drdn*drdni+S2*n_dot_ni
      ! Add
      m=m+fs_s*pphijw
      l=l+fs_d*sphijw
    end do
    ! Reverse if needed
    if (reverse) m=-m
  end subroutine fbem_bem_harpot2d_hbie_ext_st

  recursive subroutine fbem_bem_harpot2d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,p,qsp,ks,ns,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e             !! Element
    logical                            :: reverse       !! Reverse orientation
    real(kind=real64)                  :: xi_s(1,2)     !! Subdivision of the parent element
    real(kind=real64)                  :: x_i(2)        !! Collocation point position vector
    real(kind=real64)                  :: n_i(2)        !! Collocation point unit normal vector
    type(fbem_bem_harpot2d_parameters) :: p             !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp           !! Quasi-singular integration parameters
    integer                            :: ks            !! Current level of subdivisions
    integer                            :: ns            !! Maximum level of subdivision
    complex(kind=real64)               :: m(e%n_pnodes) !! m kernels matrix
    complex(kind=real64)               :: l(e%n_snodes) !! l kernels matrix
    ! Local
    integer              :: gln_near          ! 1D Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer              :: gln               ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical              :: subdivide         ! True if subdivision has to be performed
    real(kind=real64)    :: barxi(1)          ! Nearest element coordinate with respect to collocation point
    real(kind=real64)    :: barxip(1)         ! Nearest element subdivision local coordinate with respect to collocation point
    real(kind=real64)    :: rmin              ! Minimum distance between collocation point and barxi on the element
    real(kind=real64)    :: barr              ! Telles jacobian at barxi
    real(kind=real64)    :: cl                ! Characteristic length
    real(kind=real64)    :: d                 ! Normalized distance between collocation point and element subdivision
    integer              :: method            ! Method used in nearest point algorithm
    real(kind=real64)    :: tmp_xi_s(1,2)     ! Subdivision
    real(kind=real64)    :: x_s(2,e%n_gnodes) ! Coordinates of the element subdivision
    complex(kind=real64) :: m_tmp(e%n_pnodes) ! m kernels matrix (temporary)
    complex(kind=real64) :: l_tmp(e%n_snodes) ! l kernels matrix (temporary)
    ! Initialize
    if (ks.eq.1) then
      m=(0.d0,0.d0)
      l=(0.d0,0.d0)
      xi_s(1,1)=-1.d0
      xi_s(1,2)= 1.d0
      call fbem_nearest_element_point_bem(2,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
      barxip=barxi
    else
      call fbem_obtain_element_subdivision_coordinates(2,e%gtype,e%x,xi_s,x_s)
      cl=fbem_characteristic_length(2,e%gtype,x_s,1.d-12)
      call fbem_nearest_element_point_bem(2,e%gtype,x_s,cl,x_i,barxip,rmin,d,method)
    end if
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,4,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_harpot2d_hbie_ext_adp',ns,'maximum number of subdivisions reached')
        gln_near=30
      end if
    else
      if (gln_near.eq.0) subdivide=.true.
    end if
    ! Subdivide by 1/2
    if (subdivide) then
      ! SUBLINE 1
      tmp_xi_s(1,1)=xi_s(1,1)
      tmp_xi_s(1,2)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      call fbem_bem_harpot2d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_harpot2d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_harpot2d_hbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,p,gln,m_tmp,l_tmp)
      m=m+m_tmp
      l=l+l_tmp
    end if
  end subroutine fbem_bem_harpot2d_hbie_ext_adp

  subroutine fbem_bem_harpot2d_hbie_int(ngp,type_g,type_f1,type_f2,delta_f,x_nodes,reverse,xi_i,p,m,l)
    implicit none
    ! I/O
    integer                            :: ngp                                 !! Number of Gauss point to be used
    integer                            :: type_g                              !! Geometrial interpolation
    integer                            :: type_f1                             !! Functional interpolation (primary variables)
    integer                            :: type_f2                             !! Functional interpolation (secondary variables)
    real(kind=real64)                  :: delta_f                             !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64)                  :: x_nodes(2,fbem_n_nodes(type_g))     !! Position vectors of geometrical nodes
    logical                            :: reverse                             !! Normal vector inversion
    real(kind=real64)                  :: xi_i                                !! Reference coordinate of the singular point.
    type(fbem_bem_harpot2d_parameters) :: p                                   !! Parameters of the region
    complex(kind=real64)               :: m(fbem_n_nodes(type_f1))            !! m kernel vector
    complex(kind=real64)               :: l(fbem_n_nodes(type_f2))            !! l kernel vector
    ! Local
    complex(kind=real64)               :: mr(fbem_n_nodes(type_f1))           ! m kernel vector (regular part)
    complex(kind=real64)               :: lr(fbem_n_nodes(type_f2))           ! l kernel vector (regular part)
    integer                            :: nnodes_g                            ! Number of nodes of the element
    integer                            :: kphi                                ! Counter variable for shape functions loops
    real(kind=real64)                  :: x_i(2)                              ! Real coordinates of collocation point
    real(kind=real64)                  :: n_i(2)                              ! Collocation point unit normal vector
    real(kind=real64)                  :: t_i(2)                              ! Unit tangent vector at xi_i
    real(kind=real64)                  :: phi_g(fbem_n_nodes(type_g))         ! Geometrical shape functions values at xi
    real(kind=real64)                  :: dphidxi_g(fbem_n_nodes(type_g))     ! Geometrical shape functions derivatives values at xi
    real(kind=real64)                  :: phi_f1(fbem_n_nodes(type_f1))       ! Functional shape functions values at xi
    real(kind=real64)                  :: phi_f2(fbem_n_nodes(type_f2))       ! Functional shape functions values at xi
    real(kind=real64)                  :: phi_f1_i(fbem_n_nodes(type_f1))     ! Functional shape functions values at xi_i
    real(kind=real64)                  :: dphidxi_f1_i(fbem_n_nodes(type_f1)) ! Derivatives of functional shape functions values at xi_i
    integer                            :: nsub                                ! Number of subdivisions
    integer                            :: ksub                                ! Counter of subdivision
    integer                            :: kip                                 ! Counter of integration points
    real(kind=real64)                  :: gamma                               ! Coordinate gamma
    real(kind=real64)                  :: xip                                 ! Reference coordinate in subdivided element [0,1]
    real(kind=real64)                  :: js                                  ! Jacobian of the xi [xisub(1,:),xisub[2,:]] -> xip [0,1] transformation
    real(kind=real64)                  :: xip_i(2)                            ! Singular point in xip space
    real(kind=real64)                  :: w                                   ! Weights of each integration point
    real(kind=real64)                  :: xisub(2,2)                          ! Coordinates of subdivisions in xi space
    real(kind=real64)                  :: xi                                  ! Reference coordinate
    real(kind=real64)                  :: aux(10)                             ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)                  :: x(2)                                ! Position vector at xi
    real(kind=real64)                  :: T(2)                                ! Tangent vector at xi
    real(kind=real64)                  :: N(2)                                ! Normal vector at xi
    real(kind=real64)                  :: rv(2)                               ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)                  :: ra, rb                              ! Distance vector from collocation point to element vertices
    real(kind=real64)                  :: r, d1r, d1r2, logr                  ! Distance vector module and its inverse
    real(kind=real64)                  :: jg                                  ! Geometric jacobian
    real(kind=real64)                  :: jgi                                 ! Geometric jacobian at collocation point
    real(kind=real64)                  :: drdn                                ! Partial derivative of r respect to unit normal
    real(kind=real64)                  :: drdni                               ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)                  :: drdt                                ! Partial derivative of r respect to unit tangent
    real(kind=real64)                  :: n_dot_ni                            ! Scalar product between integration point unit normal and collocation point unit normal
    type(fbem_telles_parameters)       :: telles_parameters                   ! Telles parameters
    real(kind=real64)                  :: jt                                  ! Telles jacobian
    real(kind=real64)                  :: jw                                  ! Jacobians * weight
    real(kind=real64)                  :: phif1jw(fbem_n_nodes(type_f1))      ! Auxiliary variables for integrand evaluation
    real(kind=real64)                  :: phif2jw(fbem_n_nodes(type_f2))      ! Auxiliary variables for integrand evaluation
    complex(kind=real64)               :: z                                   ! Argument of Modified Bessel functions
    complex(kind=real64)               :: K0r, K1r, K2r                       ! Modified Bessel functions K values
    complex(kind=real64)               :: Q, S1, S2                           ! Components of the fundamental solution
    complex(kind=real64)               :: fs_d, fs_s                          ! Fundamental solution values
    !
    ! Initialize
    !
    ! Initialize kernel vectors
    m=(0.d0,0.d0)
    l=(0.d0,0.d0)
    mr=(0.d0,0.d0)
    lr=(0.d0,0.d0)
    ! Number of nodes of the geometrical element
    nnodes_g=fbem_n_nodes(type_g)
    ! Calculate real coordinates of collocation point
    ! Geometrical shape functions at xi_i
#   define etype type_g
#   define delta 0.0d0
#   define xi xi_i
#   define phi phi_g
#   define dphidxi dphidxi_g
#   include <phi_and_dphidxi_1d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
#   undef dphidxi
    x_i=0.d0
    T=0.d0
    do kphi=1,nnodes_g
      x_i=x_i+phi_g(kphi)*x_nodes(:,kphi)
      T=T+dphidxi_g(kphi)*x_nodes(:,kphi)
    end do
    ! Jacobian
    jgi=dot_product(T,T)
    jgi=sqrt(jgi)
    ! Calculate unit tangent
    t_i=T/jgi
    ! Unit normal
    n_i(1)=t_i(2)
    n_i(2)=-t_i(1)
    ! Functional shape functions and its first derivatives at xi_i
#   define etype type_f1
#   define delta delta_f
#   define xi xi_i
#   define phi phi_f1_i
#   define dphidxi dphidxi_f1_i
#   include <phi_and_dphidxi_1d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
#   undef dphidxi
    ! The HBIE can not be collocated at vertices
    if (fbem_check_xi_vertex(xi_i).eqv.(.true.)) then
      call fbem_error_message(output_unit,0,__FILE__,__LINE__,&
                              'the HBIE cannot be collocated at a vertex')
    ! If xi_i is inside the element, 2 subdivisions are needed
    else
      nsub=2
      ! Subdivision 1
      xisub(1,1)=-1.0d0
      xisub(2,1)=xi_i
      ! The xip coordinate of the collocation point
      xip_i(1)=1.0d0
      ! Subdivision 2
      xisub(1,2)=xi_i
      xisub(2,2)=1.0d0
      ! The xip coordinate of the collocation point
      xip_i(2)=0.0d0
    end if
    !
    ! Integrate weakly singular integrals
    !
    ! Loop through subdivisions
    do ksub=1,nsub
      ! Subdivision jacobian
      js=xisub(2,ksub)-xisub(1,ksub)
      ! Telles transformation parameters
      telles_parameters=fbem_telles01_calculate_parameters(xip_i(ksub),0.0d0)
      ! Loop through integration points
      do kip=1,gl01_n(ngp)
        ! GAMMA->XIP->XI COORDINATE TRANSFORMATION
        ! Gamma coordinate and weight
        gamma=gl01_xi(kip,ngp)
        w=gl01_w(kip,ngp)
        ! xip coordinate, weight and jacobian from Telles transformation
        call fbem_telles_xi_and_jacobian(telles_parameters,gamma,xip,jt)
        ! xi coordinate
        xi=js*xip+xisub(1,ksub)
        ! COMPONENTS OF THE FUNDAMENTAL SOLUTIONS
        ! Geometrical shape functions and first derivatives at xi
#       define etype type_g
#       define delta 0.0d0
#       define phi phi_g
#       define dphidxi dphidxi_g
#       include <phi_and_dphidxi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
#       undef dphidxi
        ! Components calculation of x and T at xi
        x=0.d0
        T=0.d0
        do kphi=1,nnodes_g
          x=x+phi_g(kphi)*x_nodes(:,kphi)
          T=T+dphidxi_g(kphi)*x_nodes(:,kphi)
        end do
        ! Normal vector
        N(1)=T(2)
        N(2)=-T(1)
        ! Geometric jacobian
        jg=dot_product(T,T)
        jg=sqrt(jg)
        ! Distance vector
        rv=x-x_i
        ! Distance vector norm
        r=dot_product(rv,rv)
        r=sqrt(r)
        logr=log(r)
        ! n dot n_i
        n_dot_ni=dot_product(N,n_i)/jg
        ! Fundamental solution value
        fs_s=p%S2(2)*logr*n_dot_ni
        ! FUNCTIONAL SHAPE FUNCTIONS
        ! Functional shape functions (primary variables) at xi
#       define etype type_f1
#       define delta delta_f
#       define phi phi_f1
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Jacobians * weight
        jw=jg*js*jt*w
        ! Functional shape functions * jacobians * weights
        phif1jw=phi_f1*jw
        ! Add to kernels
        m=m+fs_s*phif1jw
      end do ! Loop through integration points
    end do ! Loop through subdivisions
    !
    ! Integrate regular integrals of singular parts
    !
    ! Loop through subdivisions
    do ksub=1,nsub
      ! Subdivision jacobian (constant)
      js=xisub(2,ksub)-xisub(1,ksub)
      ! Loop through integration points
      do kip=1,gl01_n(ngp)
        ! XIP->XI COORDINATE TRANSFORMATION
        ! Coordinate and weight in gamma [0,1]
        xip=gl01_xi(kip,ngp)
        w=gl01_w(kip,ngp)
        ! Calculate xi
        xi=js*xip+xisub(1,ksub)
        ! COMPONENTS OF THE FUNDAMENTAL SOLUTIONS
        ! Geometrical shape functions and first derivatives at xi
#       define etype type_g
#       define delta 0.0d0
#       define phi phi_g
#       define dphidxi dphidxi_g
#       include <phi_and_dphidxi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
#       undef dphidxi
        ! Components calculation of x and T at xi
        x=0.d0
        T=0.d0
        do kphi=1,nnodes_g
          x=x+phi_g(kphi)*x_nodes(:,kphi)
          T=T+dphidxi_g(kphi)*x_nodes(:,kphi)
        end do
        ! Normal vector
        N(1)=T(2)
        N(2)=-T(1)
        ! Geometric jacobian
        jg=dot_product(T,T)
        jg=sqrt(jg)
        ! Distance vector
        rv=x-x_i
        ! Distance vector norm
        r=dot_product(rv,rv)
        r=sqrt(r)
        d1r=1.0d0/r
        d1r2=d1r**2
        logr=log(r)
        ! dr/dn
        drdn=dot_product(rv,N)*d1r/jg
        ! dr/dn_i
        drdni=-dot_product(rv,n_i)*d1r
        ! dr/dt
        drdt=dot_product(rv,T)*d1r/jg
        ! n dot n_i
        n_dot_ni=dot_product(N,n_i)/jg
        ! Components of the fundamental solutions (singular parts that lead to regular integrand)
        Q=p%Q(1)*d1r
        S1=p%S1(1)*d1r2
        ! Fundamental solution values
        fs_d=Q*drdni
        fs_s=S1*drdn*drdni
        ! FUNCTIONAL SHAPE FUNCTIONS
        ! Functional shape functions (primary variables) at xi
#       define etype type_f1
#       define delta delta_f
#       define phi phi_f1
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Functional shape functions (secondary variables) at xi
#       define etype type_f2
#       define delta delta_f
#       define phi phi_f2
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Jacobians * weight
        jw=jg*js*w
        ! Functional shape functions * jacobians* weights
        phif1jw=phi_f1*jw
        phif2jw=phi_f2*jw
        ! Add to kernels
        m=m+fs_s*phif1jw
        l=l+fs_d*phif2jw
        ! Regular parts of the HFP integral
        ! First part
        fs_s=p%S2(1)*d1r2*(n_dot_ni-abs(drdt))
        ! Kernels calculation
        m=m+fs_s*phif1jw
        ! Second part
        fs_s=p%S2(1)*d1r2*abs(drdt)
        if (ksub.eq.1) then
          m=m+fs_s*(phi_f1-phi_f1_i+dphidxi_f1_i/jgi*r)*jw
        else
          m=m+fs_s*(phi_f1-phi_f1_i-dphidxi_f1_i/jgi*r)*jw
        end if
      end do ! Loop through integration points
    end do ! Loop through subdivisions
    !
    ! Add analytical terms
    !
    ! Calculate ra and rb
    ra=sqrt((x_nodes(1,1)-x_i(1))**2+(x_nodes(2,1)-x_i(2))**2)
    rb=sqrt((x_nodes(1,2)-x_i(1))**2+(x_nodes(2,2)-x_i(2))**2)
    ! Loop through functional shape functions
    m=m+p%S2(1)*(-phi_f1_i*(1.d0/ra+1.d0/rb)+dphidxi_f1_i/jgi*(log(rb)-log(ra)))
    !
    ! Integrate regular integrals of regular parts
    !
    ! Loop through subdivisions
    do ksub=1,nsub
      ! Subdivision jacobian (constant)
      js=xisub(2,ksub)-xisub(1,ksub)
      ! Loop through integration points
      do kip=1,gl01_n(ngp)
        ! XIP->XI COORDINATE TRANSFORMATION
        ! Coordinate and weight in gamma [0,1]
        xip=gl01_xi(kip,ngp)
        w=gl01_w(kip,ngp)
        ! Calculate xi
        xi=js*xip+xisub(1,ksub)
        ! COMPONENTS OF THE FUNDAMENTAL SOLUTIONS
        ! Geometrical shape functions and first derivatives at xi
#       define etype type_g
#       define delta 0.0d0
#       define phi phi_g
#       define dphidxi dphidxi_g
#       include <phi_and_dphidxi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
#       undef dphidxi
        ! Components calculation of x and T at xi
        x=0.d0
        T=0.d0
        do kphi=1,nnodes_g
          x=x+phi_g(kphi)*x_nodes(:,kphi)
          T=T+dphidxi_g(kphi)*x_nodes(:,kphi)
        end do
        ! Normal vector
        N(1)=T(2)
        N(2)=-T(1)
        ! Geometric jacobian
        jg=dot_product(T,T)
        jg=sqrt(jg)
        ! Distance vector
        rv=x-x_i
        ! Distance vector norm
        r=dot_product(rv,rv)
        r=sqrt(r)
        d1r=1.0d0/r
        d1r2=d1r**2
        logr=log(r)
        ! dr/dn
        drdn=dot_product(rv,N)*d1r/jg
        ! dr/dn_i
        drdni=-dot_product(rv,n_i)*d1r
        ! n dot n_i
        n_dot_ni=dot_product(N,n_i)/jg
        ! Modified Bessel functions argument, z=ikr
        z=c_im*p%k*r
        call fbem_modified_bessel_K0r_K1r_K2r_2(z,K0r,K1r,K2r)
        ! Components of the fundamental solutions (regular parts)
        Q=(p%Q(2)*logr+p%Q(3))*r+p%Q(4)*K1r
        S1=p%S1(2)+(p%S1(3)*logr+p%S1(4))*r**2+p%S1(5)*K2r
        S2=p%S2(3)+p%S2(4)*d1r*K1r
        ! Fundamental solution values
        fs_d=Q*drdni
        fs_s=S1*drdn*drdni+S2*n_dot_ni
        ! FUNCTIONAL SHAPE FUNCTIONS
        ! Functional shape functions (primary variables) at xi
#       define etype type_f1
#       define delta delta_f
#       define phi phi_f1
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Functional shape functions (secondary variables) at xi
#       define etype type_f2
#       define delta delta_f
#       define phi phi_f2
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Jacobians * weight
        jw=jg*js*w
        ! Functional shape functions * jacobians* weights
        phif1jw=phi_f1*jw
        phif2jw=phi_f2*jw
        ! Add to kernels
        mr=mr+fs_s*phif1jw
        lr=lr+fs_d*phif2jw
      end do ! Loop through integration points
    end do ! Loop through subdivisions
    !
    ! Add singular and regular parts
    !
    m=m+mr
    l=l+lr
    !
    ! If the normal is inverted
    !
    if (reverse) l=-l
  end subroutine fbem_bem_harpot2d_hbie_int

  subroutine fbem_bem_harpot2d_hbie_auto(e,reverse,x_i,n_i,p,qsp,ns,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e             !! Integration element
    logical                            :: reverse       !! Reverse orientation
    real(kind=real64)                  :: x_i(2)        !! Collocation point
    real(kind=real64)                  :: n_i(2)        !! Unit normal at the collocation point
    type(fbem_bem_harpot2d_parameters) :: p             !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp           !! Quasi-singular integration parameters
    integer                            :: ns            !! Maximum level of subdivisions
    complex(kind=real64)               :: m(e%n_pnodes) !! m integration kernel
    complex(kind=real64)               :: l(e%n_snodes) !! l integration kernel
    ! Local
    real(kind=real64) :: r(2)      ! Distance vector
    real(kind=real64) :: rmin      ! Minimum distance between element and x_i
    real(kind=real64) :: barxi(1)  ! Local coordinates of the nearest element point with respect to x_i
    real(kind=real64) :: d         ! Dimensionless distance
    integer           :: delta     ! Control variable
    real(kind=real64) :: xi_s(1,2) ! Local coordinates of the element subdivision
    integer           :: method    ! Method used when calculating the nearest element point
    integer           :: gln_near  ! 1D Gauss-Legendre integration points required by the quasi-singular function
    integer           :: gln       ! 1D Gauss-Legendre integration points used in the integration
    integer           :: ps        ! Selected precalculated dataset
    integer           :: i         ! Counter
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
      call fbem_nearest_element_point_bem(2,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
      if (d.le.1.d-12) then
        delta=1
      else
        delta=0
      end if
    end if
    ! Integrate
    select case (delta)
      case (1)
        call fbem_bem_harpot2d_hbie_int(30,e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi(1),p,m,l)
      case (0)
        ! Estimate the required integration rule
        gln_near=fbem_qs_n_estimation_standard(e%n,e%gtype,4,qsp,d,barxi)
        gln=max(e%gln_far,gln_near)
        ! Integrate using a conservative precalculated dataset
        if ((gln.le.e%ps_gln_max).and.(gln_near.gt.0)) then
          do i=1,e%n_ps
            if (e%ps_gln(i).ge.gln) then
              ps=i
              exit
            end if
          end do
          call fbem_bem_harpot2d_hbie_ext_pre(ps,e,reverse,x_i,n_i,p,m,l)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_harpot2d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,p,qsp,1,ns,m,l)
        end if
    end select
  end subroutine fbem_bem_harpot2d_hbie_auto


  subroutine fbem_bem_harpot2d_hbie_bl_auto(e,x_i,n_i,p,qsp,ns,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e             !! Integration element
    real(kind=real64)                  :: x_i(2)        !! Collocation point
    real(kind=real64)                  :: n_i(2)        !! Unit normal at the collocation point
    type(fbem_bem_harpot2d_parameters) :: p             !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp           !! Quasi-singular integration parameters
    integer                            :: ns            !! Maximum level of subdivisions
    complex(kind=real64)               :: l(e%n_snodes) !! l integration kernel
    ! Local
    real(kind=real64) :: r(2)      ! Distance vector
    real(kind=real64) :: rmin      ! Minimum distance between element and x_i
    real(kind=real64) :: barxi(1)  ! Local coordinates of the nearest element point with respect to x_i
    real(kind=real64) :: d         ! Dimensionless distance
    integer           :: delta     ! Control variable
    real(kind=real64) :: xi_s(1,2) ! Local coordinates of the element subdivision
    integer           :: method    ! Method used when calculating the nearest element point
    integer           :: gln_near  ! 1D Gauss-Legendre integration points required by the quasi-singular function
    integer           :: gln       ! 1D Gauss-Legendre integration points used in the integration
    integer           :: ps        ! Selected precalculated dataset
    integer           :: i         ! Counter
    ! POINT BODY LOAD
    if (e%d.eq.0) then
      r=e%x(:,1)-x_i
      rmin=sqrt(dot_product(r,r))
      if (rmin.eq.0.d0) then
        call fbem_error_message(output_unit,0,'fbem_bem_harpot2d_hbie_bl_auto',0,'it is not possible to collocate at a point load')
      else
        call fbem_bem_harpot2d_hbie_d(e%x(:,1),x_i,n_i,p,l(1))
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
      call fbem_error_message(output_unit,0,'fbem_bem_harpot2d_hbie_bl_auto',0,'only point loads are available')
      r=e%bball_centre-x_i
      rmin=sqrt(dot_product(r,r))-e%bball_radius
      if (rmin.gt.(4.d0*e%bball_radius)) then
        delta=0
        barxi=0.d0
        d=rmin/e%cl
      else
        ! Use an adaptative algorithm that combines sampling and minimization algorithms
        call fbem_nearest_element_point_bem(2,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
        if (d.le.1.d-12) then
          delta=1
        else
          delta=0
        end if
      end if
      ! Integrate
      select case (delta)
        case (1)
          !call fbem_bem_harpot2d_hbie_bl_int(30,e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,barxi(1),p,l)
        case (0)
          ! Estimate the required integration rule
          gln_near=fbem_qs_n_estimation_standard(e%n,e%gtype,4,qsp,d,barxi)
          gln=max(e%gln_far,gln_near)
          ! Integrate using a conservative precalculated dataset
          if ((gln.le.e%ps_gln_max).and.(gln_near.gt.0)) then
            do i=1,e%n_ps
              if (e%ps_gln(i).ge.gln) then
                ps=i
                exit
              end if
            end do
            !call fbem_bem_harpot2d_hbie_bl_ext_pre(ps,e,x_i,n_i,p,l)
          ! Integrate using an adaptative algorithm
          else
            !call fbem_bem_harpot2d_hbie_bl_ext_adp(e,xi_s,x_i,n_i,p,qsp,1,ns,l)
          end if
      end select
    end if
  end subroutine fbem_bem_harpot2d_hbie_bl_auto

  ! ================================================================================================================================





  ! ================================================================================================================================
  ! BOTH SBIE AND HBIE SIMULTANEOUSLY

  subroutine fbem_bem_harpot2d_shbie_ext_pre(ps,e,reverse,x_i,n_i,pars,h,g,m,l)
    implicit none
    ! I/O
    integer                            :: ps            !! Selected precalculated dataset
    type(fbem_bem_element)             :: e             !! Element
    logical                            :: reverse       !! Reverse normal vector
    real(kind=real64)                  :: x_i(2)        !! Collocation point position vector
    real(kind=real64)                  :: n_i(2)        !! Collocation point unit normal vector
    type(fbem_bem_harpot2d_parameters) :: pars          !! Parameters of the region
    complex(kind=real64)               :: h(e%n_pnodes) !! h integration kernels matrix
    complex(kind=real64)               :: g(e%n_snodes) !! g integration kernels matrix
    complex(kind=real64)               :: m(e%n_pnodes) !! m kernels matrix
    complex(kind=real64)               :: l(e%n_snodes) !! l kernels matrix
    ! Local
    integer              :: kip                            ! Counter variable for integration points loop
    real(kind=real64)    :: x(2)                           ! Position vector at integration point
    real(kind=real64)    :: n(2)                           ! Unit normal vector at integration point
    real(kind=real64)    :: pphijw(e%n_pnodes)             ! phi^p * jacobian * weight at integration point
    real(kind=real64)    :: sphijw(e%n_snodes)             ! phi^s * jacobian * weight at integration point
    real(kind=real64)    :: rv(2)                          ! Radiovector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, d1r1, d1r2, logr            ! Radiovector module, squared, inverses and log(r)
    real(kind=real64)    :: drdx(2)                        ! Radiovector derivatives with respect to x_k
    real(kind=real64)    :: drdn                           ! Partial derivative of r respect to unit normal
    real(kind=real64)    :: drdni                          ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)    :: n_dot_ni                       ! Dot product of unit normals
    complex(kind=real64) :: z(1)                           ! Bessel functions arguments: z=ikr
    complex(kind=real64) :: KnR(0:2,1)                     ! Bessel functions decomposition
    complex(kind=real64) :: P, Q, S1, S2                   ! Components of the fundamental solution
    complex(kind=real64) :: fs_d, fs_s                     ! Fundamental solution values
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
      d1r1=1.0d0/r
      d1r2=d1r1*d1r1
      logr=log(r)
      drdx=rv*d1r1
      drdn=dot_product(drdx,n)
      drdni=-dot_product(drdx,n_i)
      n_dot_ni=dot_product(n,n_i)
      z(1)=c_im*pars%k*r
      call fbem_BesselKnR_decomposed(1,z,KnR)
      P=pars%P(1)*logr+pars%P(2)+pars%P(3)*KnR(0,1)
      Q=pars%Q(1)*d1r1+(pars%Q(2)*logr+pars%Q(3))*r+pars%Q(4)*KnR(1,1)
      S1=pars%S1(1)*d1r2+pars%S1(2)+(pars%S1(3)*logr+pars%S1(4))*r**2+pars%S1(5)*KnR(2,1)
      S2=pars%S2(1)*d1r2+pars%S2(2)*logr+pars%S2(3)+pars%S2(4)*d1r1*KnR(1,1)
      fs_d=Q*drdni
      fs_s=S1*drdn*drdni+S2*n_dot_ni
      h=h+Q*drdn*pphijw
      g=g+P*sphijw
      m=m+fs_s*pphijw
      l=l+fs_d*sphijw
    end do
    ! Reverse if needed
    if (reverse) then
      h=-h
      m=-m
    end if
  end subroutine fbem_bem_harpot2d_shbie_ext_pre

  subroutine fbem_bem_harpot2d_shbie_auto(e,reverse,x_i,n_i,p,qsp,ns,h,g,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e             !! Integration element
    logical                            :: reverse       !! Reverse orientation
    real(kind=real64)                  :: x_i(2)        !! Collocation point
    real(kind=real64)                  :: n_i(2)        !! Unit normal at the collocation point
    type(fbem_bem_harpot2d_parameters) :: p             !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp           !! Quasi-singular integration parameters
    integer                            :: ns            !! Maximum level of subdivisions
    complex(kind=real64)               :: h(e%n_pnodes) !! h integration kernel
    complex(kind=real64)               :: g(e%n_snodes) !! g integration kernel
    complex(kind=real64)               :: m(e%n_pnodes) !! m integration kernel
    complex(kind=real64)               :: l(e%n_snodes) !! l integration kernel
    ! Local
    real(kind=real64) :: r(2)      ! Distance vector
    real(kind=real64) :: rmin      ! Minimum distance between element and x_i
    real(kind=real64) :: barxi(1)  ! Local coordinates of the nearest element point with respect to x_i
    real(kind=real64) :: d         ! Dimensionless distance
    integer           :: delta     ! Control variable
    real(kind=real64) :: xi_s(1,2) ! Local coordinates of the element subdivision
    integer           :: method    ! Method used when calculating the nearest element point
    integer           :: gln_near  ! 1D Gauss-Legendre integration points required by the quasi-singular function
    integer           :: gln       ! 1D Gauss-Legendre integration points used in the integration
    integer           :: ps        ! Selected precalculated dataset
    integer           :: i         ! Counter
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
      call fbem_nearest_element_point_bem(2,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
      if (d.le.1.d-12) then
        delta=1
      else
        delta=0
      end if
    end if
    ! Integrate
    select case (delta)
      case (1)
        call fbem_bem_harpot2d_sbie_int(30,e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi(1),p,h,g)
        call fbem_bem_harpot2d_hbie_int(30,e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi(1),p,m,l)
      case (0)
        ! Estimate the required integration rule
        gln_near=fbem_qs_n_estimation_standard(e%n,e%gtype,4,qsp,d,barxi)
        gln=max(e%gln_far,gln_near)
        ! Integrate using a conservative precalculated dataset
        if ((gln.le.e%ps_gln_max).and.(gln_near.gt.0)) then
          do i=1,e%n_ps
            if (e%ps_gln(i).ge.gln) then
              ps=i
              exit
            end if
          end do
          call fbem_bem_harpot2d_shbie_ext_pre(ps,e,reverse,x_i,n_i,p,h,g,m,l)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_harpot2d_sbie_ext_adp(e,reverse,xi_s,x_i,p,qsp,1,ns,h,g)
          call fbem_bem_harpot2d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,p,qsp,1,ns,m,l)
        end if
    end select
  end subroutine fbem_bem_harpot2d_shbie_auto

  ! ================================================================================================================================







  ! ================================================================================================================================
  ! VARIATION SINGULAR BOUNDARY INTEGRAL EQUATION (VSBIE)

  subroutine fbem_bem_harpot2d_vsbie_ext_pre(ps,e,reverse,x_i,p,h1,h2,g1,g2)
    implicit none
    ! I/O
    integer                            :: ps                          !! Selected precalculated dataset
    type(fbem_bem_element)             :: e                           !! Element
    logical                            :: reverse                     !! Reverse normal vector
    real(kind=real64)                  :: x_i(2)                      !! Collocation point position vector
    type(fbem_bem_harpot2d_parameters) :: p                           !! Parameters of the region
    complex(kind=real64)               :: h1(e%dme_n_gnodes,e%n_pnodes,2) !! h1 matrix
    complex(kind=real64)               :: h2(               e%n_pnodes,2) !! h2 matrix
    complex(kind=real64)               :: g1(e%dme_n_gnodes,e%n_snodes,2) !! g1 matrix
    complex(kind=real64)               :: g2(               e%n_snodes,2) !! g2 matrix
    ! Local
    integer              :: im, iq                    ! Counter variables
    integer              :: kip                       ! Counter variable for integration points loop
    real(kind=real64)    :: x(2)                      ! Position vector at integration point
    real(kind=real64)    :: t(2)                      ! Unit tangent vector at integration point
    real(kind=real64)    :: n(2)                      ! Unit normal vector at integration point
    real(kind=real64)    :: pphijw(e%n_pnodes)        ! phi^p * jacobian * weight at integration point
    real(kind=real64)    :: sphijw(e%n_snodes)        ! phi^s at integration point
    real(kind=real64)    :: rv(2)                     ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, dr1, dr2, logr         ! Distance vector module, inverse and log(1/r)
    real(kind=real64)    :: drdx(2)                   ! Distance vector derivatives with respect to x_k
    real(kind=real64)    :: drdn                      ! Partial derivative of r respect to unit normal
    real(kind=real64)    :: drdt                      ! Partial derivative of r respect to unit tangent
    real(kind=real64)    :: dme_gphi(e%dme_n_gnodes)  ! phi^{dme-g} at integration point
    real(kind=real64)    :: dme_wqj(e%dme_n_gnodes,2) ! dphi{dme-g}/dx_j at integration point
    real(kind=real64)    :: dme_vq(e%dme_n_gnodes)    ! vq=w_{j,q}·t_j at integration point
    complex(kind=real64) :: z(1)                      ! Wavenumbers
    complex(kind=real64) :: KnR(0:2,1)                ! Bessel functions decomposition
    complex(kind=real64) :: PP, Q, S1, S2             ! Fundamental solution components
    complex(kind=real64) :: fs_p, fs_q                ! Fundamental solutions values
    complex(kind=real64) :: fs_qt, fs_dpdm, fs_dqdm   ! Fundamental solutions values
    ! Initialize matrices
    h1=(0.d0,0.d0)
    h2=(0.d0,0.d0)
    g1=(0.d0,0.d0)
    g2=(0.d0,0.d0)
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      ! Components of the integrand
      x=e%ps_x(:,kip,ps)
      t=e%ps_t(:,1,kip,ps)
      n=e%ps_n(:,kip,ps)
      pphijw=e%ps_pphijw(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      dme_gphi=e%ps_dme_gphi(:,kip,ps)
      dme_wqj=e%ps_dme_dgphidx(:,:,kip,ps)
      dme_vq(:)=dme_wqj(:,1)*t(1)+dme_wqj(:,2)*t(2)
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      dr1=1.d0/r
      dr2=1.d0/r**2
      logr=log(r)
      drdx=rv*dr1
      drdn=dot_product(drdx,n)
      drdt=dot_product(drdx,t)
      ! Calculation of the components of the fundamental solution
      z(1)=c_im*p%k*r
      call fbem_BesselKnR_decomposed(1,z,KnR)
      PP=p%P(1)*logr+p%P(2)+p%P(3)*KnR(0,1)
      Q=p%Q(1)*dr1+(p%Q(2)*logr+p%Q(3))*r+p%Q(4)*KnR(1,1)
      S1=p%S1(1)*dr2+p%S1(2)+(p%S1(3)*logr+p%S1(4))*r**2+p%S1(5)*KnR(2,1)
      S2=p%S2(1)*dr2+p%S2(2)*logr+p%S2(3)+p%S2(4)*dr1*KnR(1,1)
      ! Build kernels and add integrand*weight
      fs_p=PP
      fs_q=Q*drdn
      fs_qt=Q*drdt
      do im=1,2
        fs_dpdm=Q*drdx(im)
        fs_dqdm=S1*drdx(im)*drdn-S2*n(im)
        do iq=1,e%dme_n_gnodes
          h1(iq,:,im)=h1(iq,:,im)+fs_q*dme_vq(iq)*t(im)*pphijw(:)
          h1(iq,:,im)=h1(iq,:,im)+fs_dqdm*dme_gphi(iq)*pphijw(:)
          h1(iq,:,im)=h1(iq,:,im)-fs_qt*dme_vq(iq)*n(im)*pphijw(:)
          g1(iq,:,im)=g1(iq,:,im)+fs_p*dme_vq(iq)*t(im)*sphijw(:)
          g1(iq,:,im)=g1(iq,:,im)+fs_dpdm*dme_gphi(iq)*sphijw(:)
        end do
        h2(:,im)=h2(:,im)+fs_dqdm*pphijw(:)
        g2(:,im)=g2(:,im)+fs_dpdm*sphijw(:)
      end do
    end do
    ! Reverse element orientation
    if (reverse) then
      h1=-h1
      h2=-h2
    end if
  end subroutine fbem_bem_harpot2d_vsbie_ext_pre

  subroutine fbem_bem_harpot2d_vsbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,p,gln,h1,h2,g1,g2)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                           !! Integration element
    logical                            :: reverse                     !! Reverse normal vector
    real(kind=real64)                  :: xi_s(1,2)                   !! Coordinates of the subdivision of the parent element (xi space [-1,1])
    real(kind=real64)                  :: x_i(2)                      !! Collocation point position vector
    real(kind=real64)                  :: barxip(1)                   !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)                  :: barr                        !! Telles jacobian at barxip
    type(fbem_bem_harpot2d_parameters) :: p                           !! Parameters of the region
    integer                            :: gln                         !! 1D Gauss-Legendre number of integration points (<=32)
    complex(kind=real64)               :: h1(e%dme_n_gnodes,e%n_pnodes,2) !! h1 matrix
    complex(kind=real64)               :: h2(               e%n_pnodes,2) !! h2 matrix
    complex(kind=real64)               :: g1(e%dme_n_gnodes,e%n_snodes,2) !! g1 matrix
    complex(kind=real64)               :: g2(               e%n_snodes,2) !! g2 matrix
    ! Local
    integer                      :: im, iq                    ! Counter variables
    integer                      :: kphi                      ! Counter variable for shape functions loops
    integer                      :: kip                       ! Counter variable of integration points
    real(kind=real64)            :: gamma                     ! Coordinate gamma (Telles transformation space [-1,1])
    real(kind=real64)            :: w                         ! Weights of an integration point
    type(fbem_telles_parameters) :: telles_parameters         ! Telles parameters
    real(kind=real64)            :: jt                        ! Telles jacobian: xip->gamma
    real(kind=real64)            :: xip                       ! Coordinate xip (subdivision space [-1,1])
    real(kind=real64)            :: js                        ! Subdivision jacobian: xi->xip
    real(kind=real64)            :: xi                        ! Coordinate xi [xi_s(1,1),xi_s(1,2)]
    real(kind=real64)            :: aux(10)                   ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)          ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi(e%n_gnodes)      ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: pphi(e%n_pnodes)          ! Functional shape functions values
    real(kind=real64)            :: sphi(e%n_snodes)          ! Functional shape functions values
    real(kind=real64)            :: x(2)                      ! Position vector at xi
    real(kind=real64)            :: T(2)                      ! Tangent vector at xi
    real(kind=real64)            :: N(2)                      ! Normal vector at xi
    real(kind=real64)            :: rv(2)                     ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)            :: r, dr1, dr2, logr         ! Distance vector module, its inverse and log(1/r)
    real(kind=real64)            :: drdx(2)                   ! Distance vector derivatives with respect to x_k
    real(kind=real64)            :: jg                        ! Geometric jacobian
    real(kind=real64)            :: drdn                      ! Partial derivative of r respect to unit normal
    real(kind=real64)            :: drdt                      ! Partial derivative of r respect to unit tangent
    real(kind=real64)            :: jw                        ! Jacobians (except geometric jacobian) * weight
    real(kind=real64)            :: pphijw(e%n_pnodes)        ! Functional shape functions values * jw
    real(kind=real64)            :: sphijw(e%n_snodes)        ! Functional shape functions values * jw
    real(kind=real64)            :: dme_gphi(e%dme_n_gnodes)  ! phi^{dme-g} at integration point
    real(kind=real64)            :: dme_wqj(e%dme_n_gnodes,2) ! w_{q,j}=dphi^{dme}_q/dx_j at integration point
    real(kind=real64)            :: dme_vq(e%dme_n_gnodes)    ! v_{q}=w_{q,j}·t_j at integration point
    real(kind=real64)            :: dme_d                     ! Dimensionless distance between DME and a point
    real(kind=real64)            :: dme_xi(e%dme_d)           ! Local coordinate of the DME
    complex(kind=real64)         :: z(1)                      ! Wavenumbers
    complex(kind=real64)         :: KnR(0:2,1)                ! Bessel functions decomposition
    complex(kind=real64)         :: PP, Q, S1, S2             ! Fundamental solution components
    complex(kind=real64)         :: fs_p, fs_q                ! Fundamental solutions values
    complex(kind=real64)         :: fs_qt, fs_dpdm, fs_dqdm   ! Fundamental solutions values
    ! Initialize matrices
    h1=(0.d0,0.d0)
    h2=(0.d0,0.d0)
    g1=(0.d0,0.d0)
    g2=(0.d0,0.d0)
    ! Initialize dme_xi
    dme_xi=0.d0
    ! Calculate Telles parameters
    telles_parameters=fbem_telles11_calculate_parameters(barxip(1),barr)
    ! Numerical integration
    do kip=1,gl11_n(gln)
      ! GAMMA COORDINATE
      gamma=gl11_xi(kip,gln)
      w=gl11_w(kip,gln)
      ! GAMMA->XIP TRANSFORMATION
      call fbem_telles_xi_and_jacobian(telles_parameters,gamma,xip,jt)
      ! XIP->XI TRANSFORMATION
      xi=0.5d0*(1.d0-xip)*xi_s(1,1)+0.5d0*(1.d0+xip)*xi_s(1,2)
      js=0.5d0*(xi_s(1,2)-xi_s(1,1))
      ! XI->X TRANSFORMATION
#     define etype e%gtype
#     define delta 0.d0
#     define phi gphi
#     define dphidxi dgphidxi
#     include <phi_and_dphidxi_1d.rc>
#     undef etype
#     undef delta
#     undef phi
#     undef dphidxi
      x=0.d0
      T=0.d0
      do kphi=1,e%n_gnodes
        x=x+gphi(kphi)*e%x(:,kphi)
        T=T+dgphidxi(kphi)*e%x(:,kphi)
      end do
      ! Geometric jacobian
      jg=sqrt(dot_product(T,T))
      ! Unit tangent
      t=T/jg
      ! Normal vector
      n(1)=t(2)
      n(2)=-t(1)
      ! Distance vector and relatives
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      dr1=1.d0/r
      dr2=1.d0/r**2
      logr=log(r)
      drdx=rv*dr1
      drdn=dot_product(drdx,n)
      drdt=dot_product(drdx,t)
      ! jacobians * weight
      jw=jg*js*jt*w
      ! DESIGN MACRO-ELEMENT VALUES AT INTEGRATION POINT
      ! Calculate the local coordinates of x_i in the design macro-element
      call fbem_local_coordinates(2,e%dmetype,e%dme_x,e%dme_cl,x,dme_xi,dme_d)
      if (dme_d.gt.1.d-12) then
        write(error_unit,*) 'd=',dme_d
        write(error_unit,*) 'dme_x=',e%dme_x
        write(error_unit,*) 'x=',x
        write(error_unit,*) 'xi=',dme_xi
        stop 'the integration point is not in the design macro-element'
      end if
      ! Shape functions at dme_xi
      dme_gphi=fbem_phi_hybrid(e%dmetype,0.d0,dme_xi)
      ! Derivatives of the shape functions with respect to x_j
      dme_wqj=fbem_dphidx(2,e%dmetype,0.d0,e%dme_x,dme_xi)
      ! v_{q}=w_{q,j}·t_j
      dme_vq(:)=dme_wqj(:,1)*t(1)+dme_wqj(:,2)*t(2)
      ! FUNCTIONAL SHAPE FUNCTIONS
      ! Functional shape functions (primary variables) at xi
#     define etype e%ptype
#     define delta e%ptype_delta
#     define phi pphi
#     include <phi_1d.rc>
#     undef etype
#     undef delta
#     undef phi
      ! Functional shape functions (secondary variables) at xi
#     define etype e%stype
#     define delta e%stype_delta
#     define phi sphi
#     include <phi_1d.rc>
#     undef etype
#     undef delta
#     undef phi
      ! functional shape functions * jw
      pphijw=pphi*jw
      sphijw=sphi*jw
      ! Calculation of the components of the fundamental solution
      z(1)=c_im*p%k*r
      call fbem_BesselKnR_decomposed(1,z,KnR)
      PP=p%P(1)*logr+p%P(2)+p%P(3)*KnR(0,1)
      Q=p%Q(1)*dr1+(p%Q(2)*logr+p%Q(3))*r+p%Q(4)*KnR(1,1)
      S1=p%S1(1)*dr2+p%S1(2)+(p%S1(3)*logr+p%S1(4))*r**2+p%S1(5)*KnR(2,1)
      S2=p%S2(1)*dr2+p%S2(2)*logr+p%S2(3)+p%S2(4)*dr1*KnR(1,1)
      ! Build kernels and add integrand*weight
      fs_p=PP
      fs_q=Q*drdn
      fs_qt=Q*drdt
      do im=1,2
        fs_dpdm=Q*drdx(im)
        fs_dqdm=S1*drdx(im)*drdn-S2*n(im)
        do iq=1,e%dme_n_gnodes
          h1(iq,:,im)=h1(iq,:,im)+fs_q*dme_vq(iq)*t(im)*pphijw(:)
          h1(iq,:,im)=h1(iq,:,im)+fs_dqdm*dme_gphi(iq)*pphijw(:)
          h1(iq,:,im)=h1(iq,:,im)-fs_qt*dme_vq(iq)*n(im)*pphijw(:)
          g1(iq,:,im)=g1(iq,:,im)+fs_p*dme_vq(iq)*t(im)*sphijw(:)
          g1(iq,:,im)=g1(iq,:,im)+fs_dpdm*dme_gphi(iq)*sphijw(:)
        end do
        h2(:,im)=h2(:,im)+fs_dqdm*pphijw(:)
        g2(:,im)=g2(:,im)+fs_dpdm*sphijw(:)
      end do
    end do
    ! Reverse element orientation
    if (reverse) then
      h1=-h1
      h2=-h2
    end if
  end subroutine fbem_bem_harpot2d_vsbie_ext_st

  recursive subroutine fbem_bem_harpot2d_vsbie_ext_adp(e,reverse,xi_s,x_i,p,qsp,ks,ns,h1,h2,g1,g2)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                           !! Element
    logical                            :: reverse                     !! Reverse orientation
    real(kind=real64)                  :: xi_s(1,2)                   !! Subdivision of the parent element
    real(kind=real64)                  :: x_i(2)                      !! Collocation point position vector
    type(fbem_bem_harpot2d_parameters) :: p                           !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp                         !! Quasi-singular integration parameters
    integer                            :: ks                          !! Current level of subdivisions
    integer                            :: ns                          !! Maximum level of subdivision
    complex(kind=real64)               :: h1(e%dme_n_gnodes,e%n_pnodes,2) !! h1 matrix
    complex(kind=real64)               :: h2(               e%n_pnodes,2) !! h2 matrix
    complex(kind=real64)               :: g1(e%dme_n_gnodes,e%n_snodes,2) !! g1 matrix
    complex(kind=real64)               :: g2(               e%n_snodes,2) !! g2 matrix
    ! Local
    integer              :: gln_near                        ! 1D Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer              :: gln                             ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical              :: subdivide                       ! True if subdivision has to be performed
    real(kind=real64)    :: barxi(1)                        ! Nearest element coordinate with respect to collocation point
    real(kind=real64)    :: barxip(1)                       ! Nearest element subdivision local coordinate with respect to collocation point
    real(kind=real64)    :: rmin                            ! Minimum distance between collocation point and barxi on the element
    real(kind=real64)    :: barr                            ! Telles jacobian at barxi
    real(kind=real64)    :: cl                              ! Characteristic length
    real(kind=real64)    :: d                               ! Normalized distance between collocation point and element subdivision
    integer              :: method                          ! Method used in nearest point algorithm
    real(kind=real64)    :: tmp_xi_s(1,2)                   ! Subdivision
    real(kind=real64)    :: x_s(2,e%n_gnodes)               ! Coordinates of the element subdivision
    complex(kind=real64) :: h1_tmp(e%dme_n_gnodes,e%n_pnodes,2) ! h1 matrix (temporary)
    complex(kind=real64) :: h2_tmp(               e%n_pnodes,2) ! h2 matrix (temporary)
    complex(kind=real64) :: g1_tmp(e%dme_n_gnodes,e%n_snodes,2) ! g1 matrix (temporary)
    complex(kind=real64) :: g2_tmp(               e%n_snodes,2) ! g2 matrix (temporary)
    ! Initialize
    if (ks.eq.1) then
      h1=(0.d0,0.d0)
      h2=(0.d0,0.d0)
      g1=(0.d0,0.d0)
      g2=(0.d0,0.d0)
      xi_s(1,1)=-1.d0
      xi_s(1,2)= 1.d0
      call fbem_nearest_element_point_bem(2,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
      barxip=barxi
    else
      call fbem_obtain_element_subdivision_coordinates(2,e%gtype,e%x,xi_s,x_s)
      cl=fbem_characteristic_length(2,e%gtype,x_s,1.d-12)
      call fbem_nearest_element_point_bem(2,e%gtype,x_s,cl,x_i,barxip,rmin,d,method)
    end if
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,4,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_harpot2d_vsbie_ext_adp',ns,'maximum number of subdivisions reached')
        gln_near=30
      end if
    else
      if (gln_near.eq.0) subdivide=.true.
    end if
    ! Subdivide by 1/2
    if (subdivide) then
      ! SUBLINE 1
      tmp_xi_s(1,1)=xi_s(1,1)
      tmp_xi_s(1,2)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      call fbem_bem_harpot2d_vsbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h1,h2,g1,g2)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_harpot2d_vsbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h1,h2,g1,g2)
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_harpot2d_vsbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,p,gln,h1_tmp,h2_tmp,g1_tmp,g2_tmp)
      h1=h1+h1_tmp
      h2=h2+h2_tmp
      g1=g1+g1_tmp
      g2=g2+g2_tmp
    end if
  end subroutine fbem_bem_harpot2d_vsbie_ext_adp

  subroutine fbem_bem_harpot2d_vsbie_int(e,reverse,xi_i,p,h1,g1)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                           !! Integration element
    logical                            :: reverse                     !! Reverse normal vector
    real(kind=real64)                  :: xi_i(1)                     !! Local coordinate of the singular point.
    type(fbem_bem_harpot2d_parameters) :: p                           !! Parameters of the region
    complex(kind=real64)               :: h1(e%dme_n_gnodes,e%n_pnodes,2) !! h1 matrix
    complex(kind=real64)               :: g1(e%dme_n_gnodes,e%n_snodes,2) !! g1 matrix
    ! Local
    integer              :: gln                              ! 1D Gauss-Legendre number of integration points (<=32)
    integer              :: im, iq                           ! Counter variables
    integer              :: kphi                             ! Counter variable for shape functions loops
    integer              :: kip                              ! Counter of integration points
    real(kind=real64)    :: x_i(2)                           ! Real coordinates of collocation point
    real(kind=real64)    :: gphi(e%n_gnodes)                 ! Geometrical shape functions values at xi
    real(kind=real64)    :: dgphidxi(e%n_gnodes)             ! Geometrical shape functions derivatives values at xi
    real(kind=real64)    :: gphi_i(e%n_gnodes)               ! Geometrical shape functions values at xi_i
    real(kind=real64)    :: dgphidxi_i(e%n_gnodes)           ! Geometrical shape functions derivatives values at xi_i
    real(kind=real64)    :: pphi(e%n_pnodes)                 ! Functional shape functions values at xi
    real(kind=real64)    :: sphi(e%n_snodes)                 ! Functional shape functions values at xi
    real(kind=real64)    :: pphi_i(e%n_pnodes)               ! Functional shape functions values at xi_i
    integer              :: nsub                             ! Number of subdivision of the element
    integer              :: ksub                             ! Counter of subdivision
    real(kind=real64)    :: w                                ! Weights of each integration point
    real(kind=real64)    :: xip                              ! Coordinate  xip of subdivided element [0,1]
    real(kind=real64)    :: js                               ! Jacobian of the xi [xisub(1,:),xisub[2,:]] -> xip [0,1] transformation
    real(kind=real64)    :: xip_i(2)                         ! Singular point in xip space
    real(kind=real64)    :: xi                               ! Coordinate xi
    real(kind=real64)    :: xisub(2,2)                       ! Coordinates of element subdivisions
    real(kind=real64)    :: aux(10)                          ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)    :: x(2)                             ! Position vector at xi
    real(kind=real64)    :: T(2)                             ! Tangent vector at xi
    real(kind=real64)    :: N(2)                             ! Normal vector at xi
    real(kind=real64)    :: t_i(2), n_i(2)                   ! Unit tangent and normal at xi_i
    real(kind=real64)    :: rv(2)                            ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, dr1, dr2, logr                ! Distance vector module, its inverse and log(1/r)
    real(kind=real64)    :: drdx(2)                          ! Distance vector derivatives with respect to x_k
    real(kind=real64)    :: drdti                            ! (dr/dGamma)^i=+-1, -1 antes de i, 1 despues de i
    real(kind=real64)    :: jg                               ! Geometric jacobian
    real(kind=real64)    :: jgi                              ! Geometric jacobian at xi_i
    real(kind=real64)    :: drdn                             ! Partial derivative of r respect to unit normal
    real(kind=real64)    :: drdt                             ! Partial derivative of r respect to unit tangent
    real(kind=real64)    :: jw                               ! Jacobians * weight
    real(kind=real64)    :: dme_gphi(e%dme_n_gnodes)         ! phi^{dme-g} at integration point
    real(kind=real64)    :: dme_gphi_i(e%dme_n_gnodes)       ! phi^{dme-g} at collocation point
    real(kind=real64)    :: dme_wqj(e%dme_n_gnodes,2)        ! w_{q,j}=dphi^{dme}_q/dx_j at integration point
    real(kind=real64)    :: dme_wqj_i(e%dme_n_gnodes,2)      ! w_{q,j}=dphi^{dme}_q/dx_j at collocation point
    real(kind=real64)    :: dme_vq(e%dme_n_gnodes)           ! v_{q}=w_{q,j}·t_j at integration point
    real(kind=real64)    :: dme_vq_i(e%dme_n_gnodes)         ! v_{q}=w_{q,j}·t_j at collocation point
    real(kind=real64)    :: dme_d                            ! Dimensionless distance between DME and a point
    real(kind=real64)    :: dme_xi(e%dme_d)                  ! Local coordinate of the DME
    complex(kind=real64) :: z(1)                             ! Wavenumbers
    complex(kind=real64) :: KnR(0:2,1)                       ! Bessel functions decomposition
    complex(kind=real64) :: PP, Q, S1, S2                    ! Fundamental solution components
    complex(kind=real64) :: fs_p, fs_q                       ! Fundamental solutions values
    complex(kind=real64) :: fs_qt, fs_dpdm, fs_dqdm          ! Fundamental solutions values
    complex(kind=real64) :: h(e%n_pnodes)                    ! h matrix
    complex(kind=real64) :: g(e%n_snodes)                    ! g matrix
    !
    ! Initialization
    !
    ! Integration points to be used
    gln=32
    ! Initialize matrices
    h1=(0.d0,0.d0)
    g1=(0.d0,0.d0)
    ! Initialize dme_xi
    dme_xi=0.d0
    ! If xi_i belongs to one of the vertices, subdivision is not needed
    if (fbem_check_xi_vertex(xi_i(1))) then
      nsub=1
      xisub(1,1)=-1.0d0
      xisub(2,1)=1.0d0
      if (xi_i(1).lt.0.0d0) xip_i(1)=0.0d0
      if (xi_i(1).gt.0.0d0) xip_i(1)=1.0d0
    ! If xi_i is inside the element, 2 subdivisions are needed
    else
      nsub=2
      xisub(1,1)=-1.0d0
      xisub(2,1)=xi_i(1)
      xip_i(1)=1.0d0
      xisub(1,2)=xi_i(1)
      xisub(2,2)=1.0d0
      xip_i(2)=0.0d0
    end if
    ! INTEGRATION ELEMENT VALUES AT COLLOCATION POINT
    ! Calculate spatial coordinates at xi_i
#   define etype e%gtype
#   define delta 0.d0
#   define xi xi_i(1)
#   define phi gphi_i
#   define dphidxi dgphidxi_i
#   include <phi_and_dphidxi_1d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
#   undef dphidxi
    x_i=0.d0
    T_i=0.d0
    do kphi=1,e%n_gnodes
      x_i=x_i+gphi_i(kphi)*e%x(:,kphi)
      T_i=T_i+dgphidxi_i(kphi)*e%x(:,kphi)
    end do
    ! Jacobian
    jgi=sqrt(dot_product(T_i,T_i))
    ! Calculate unit tangent
    t_i=T_i/jgi
    ! Unit normal
    n_i(1)=t_i(2)
    n_i(2)=-t_i(1)
    ! Functional shape functions (primary variable) at xi_i
#   define etype e%ptype
#   define delta e%ptype_delta
#   define xi xi_i(1)
#   define phi pphi_i
#   include <phi_1d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
#   undef dphidxi
    ! DESIGN MACRO-ELEMENT VALUES AT COLLOCATION POINT
    ! Calculate the local coordinates of x_i in the design macro-element
    call fbem_local_coordinates(2,e%dmetype,e%dme_x,e%dme_cl,x_i,dme_xi,dme_d)
    if (dme_d.gt.1.d-12) then
      write(error_unit,*) 'd=',dme_d
      write(error_unit,*) 'dme_x=',e%dme_x
      write(error_unit,*) 'x_i=',x_i
      write(error_unit,*) 'xi=',dme_xi
      stop 'the collocation point is not in the design macro-element'
    end if
    ! Shape functions at dme_xi
    dme_gphi_i=fbem_phi_hybrid(e%dmetype,0.d0,dme_xi)
    ! Derivatives of the shape functions with respect to x_j
    dme_wqj_i=fbem_dphidx(2,e%dmetype,0.d0,e%dme_x,dme_xi)
    ! v_{q}=w_{q,j}·t_j
    dme_vq_i(:)=dme_wqj_i(:,1)*t_i(1)+dme_wqj_i(:,2)*t_i(2)
    !
    ! Numerical integration of regular integrals
    !
    do ksub=1,nsub
      do kip=1,gl01_n(gln)
        ! XIP COORDINATE
        xip=gl01_xi(kip,gln)
        w=gl01_w(kip,gln)
        ! XIP->XI TRANSFORMATION
        js=xisub(2,ksub)-xisub(1,ksub)
        xi=js*xip+xisub(1,ksub)
        ! dr/dGamma
        if (xi.lt.xi_i(1)) then
          drdti=-1.d0
        else
          drdti= 1.d0
        end if
        ! XI->X TRANSFORMATION
#       define etype e%gtype
#       define delta 0.d0
#       define phi gphi
#       define dphidxi dgphidxi
#       include <phi_and_dphidxi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
#       undef dphidxi
        x=0.d0
        T=0.d0
        do kphi=1,e%n_gnodes
          x=x+gphi(kphi)*e%x(:,kphi)
          T=T+dgphidxi(kphi)*e%x(:,kphi)
        end do
        ! xi->x jacobian
        jg=sqrt(dot_product(T,T))
        ! Unit tangent vector
        t=T/jg
        ! Unit normal vector
        n(1)=t(2)
        n(2)=-t(1)
        ! Distance vector and relatives
        rv=x-x_i
        r=sqrt(dot_product(rv,rv))
        dr1=1.d0/r
        dr2=1.d0/r**2
        logr=log(r)
        drdx=rv*dr1
        drdn=dot_product(drdx,n)
        drdt=dot_product(drdx,t)
        ! Jacobian * weight
        jw=jg*js*w
        ! DESIGN MACRO-ELEMENT VALUES AT INTEGRATION POINT
        ! Calculate the local coordinates of x_i in the design macro-element
        call fbem_local_coordinates(2,e%dmetype,e%dme_x,e%dme_cl,x,dme_xi,dme_d)
        if (dme_d.gt.1.d-12) then
          write(error_unit,*) 'd=',dme_d
          write(error_unit,*) 'dme_x=',e%dme_x
          write(error_unit,*) 'x=',x
          write(error_unit,*) 'xi=',dme_xi
          stop 'the integration point is not in the design macro-element'
        end if
        ! Shape functions at dme_xi
        dme_gphi=fbem_phi_hybrid(e%dmetype,0.d0,dme_xi)
        ! Derivatives of the shape functions with respect to x_j
        dme_wqj=fbem_dphidx(2,e%dmetype,0.d0,e%dme_x,dme_xi)
        ! v_{q}=w_{q,j}·t_j
        dme_vq(:)=dme_wqj(:,1)*t(1)+dme_wqj(:,2)*t(2)
        ! FUNCTIONAL SHAPE FUNCTIONS
        ! Functional shape functions (primary variables) at xi
#       define etype e%ptype
#       define delta e%ptype_delta
#       define phi pphi
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Functional shape functions (secondary variables) at xi
#       define etype e%stype
#       define delta e%stype_delta
#       define phi sphi
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! COMPONENTS OF THE FUNDAMENTAL SOLUTION
        ! Calculation of the components of the fundamental solution (only regular part)
        z(1)=c_im*p%k*r
        call fbem_BesselKnR_decomposed(1,z,KnR)
        PP=p%P(1)*logr+p%P(2)+p%P(3)*KnR(0,1)
        Q=(p%Q(2)*logr+p%Q(3))*r+p%Q(4)*KnR(1,1)
        S1=p%S1(2)+(p%S1(3)*logr+p%S1(4))*r**2+p%S1(5)*KnR(2,1)
        S2=p%S2(2)*logr+p%S2(3)+p%S2(4)*dr1*KnR(1,1)
        ! Build kernels and add integrand*weight
        fs_p=PP
        fs_q=(p%Q(1)*dr1+Q)*drdn
        fs_qt=Q*drdt
        do im=1,2
          fs_dpdm=(p%Q(1)*dr1+Q)*drdx(im)
          fs_dqdm=S1*drdx(im)*drdn-S2*n(im)
          do iq=1,e%dme_n_gnodes
            ! G^s
            g1(iq,:,im)=g1(iq,:,im)+fs_p*(dme_vq(iq)*t(im)-dme_vq_i(iq)*t_i(im))*sphi(:)*jw
            ! H^s
            h1(iq,:,im)=h1(iq,:,im)+fs_q*(dme_vq(iq)*t(im)-dme_vq_i(iq)*t_i(im))*pphi(:)*jw
            ! G^r
            g1(iq,:,im)=g1(iq,:,im)+fs_dpdm*(dme_gphi(iq)-dme_gphi_i(iq))*sphi(:)*jw
            ! H^r
            h1(iq,:,im)=h1(iq,:,im)+fs_dqdm*(dme_gphi(iq)-dme_gphi_i(iq))*pphi(:)*jw
            h1(iq,:,im)=h1(iq,:,im)+c_1_2pi*2.d0*drdx(im)*dr2*drdn*(dme_gphi(iq)-dme_gphi_i(iq))*pphi(:)*jw
            h1(iq,:,im)=h1(iq,:,im)-c_1_2pi*n(im)*dr2*(dme_gphi(iq)-dme_gphi_i(iq)-dot_product(dme_wqj_i(iq,:),rv))*pphi(:)*jw
            h1(iq,:,im)=h1(iq,:,im)-c_1_2pi*dr1*(dot_product(dme_wqj_i(iq,:),drdx)*n(im)*pphi(:)-drdti*dme_vq_i(iq)*n_i(im)*pphi_i(:))*jw
            h1(iq,:,im)=h1(iq,:,im)-c_1_2pi*dme_vq_i(iq)*n_i(im)*pphi_i(:)*dr1*(drdti-drdt)*jw
            ! H^n
            h1(iq,:,im)=h1(iq,:,im)-fs_qt*dme_vq(iq)*n(im)*pphi(:)*jw
            h1(iq,:,im)=h1(iq,:,im)+c_1_2pi*dr1*drdt*(dme_vq(iq)*n(im)*pphi(:)-dme_vq_i(iq)*n_i(im)*pphi_i(:))*jw
          end do
        end do
      end do
    end do
    ! Assemble matrices from SBIE
    call fbem_bem_harpot2d_sbie_int(gln,e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,.false.,xi_i(1),p,h,g)
    do im=1,2
      do iq=1,e%dme_n_gnodes
        h1(iq,:,im)=h1(iq,:,im)+dme_vq_i(iq)*t_i(im)*h(:)
        g1(iq,:,im)=g1(iq,:,im)+dme_vq_i(iq)*t_i(im)*g(:)
      end do
    end do
    ! Reverse element orientation
    if (reverse) h1=-h1
  end subroutine fbem_bem_harpot2d_vsbie_int

  subroutine fbem_bem_harpot2d_vsbie_auto(e,reverse,x_i,p,qsp,ns,h1,h2,g1,g2)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                           !! Integration element
    logical                            :: reverse                     !! Reverse orientation
    real(kind=real64)                  :: x_i(2)                      !! Collocation point
    type(fbem_bem_harpot2d_parameters) :: p                           !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp                         !! Quasi-singular integration parameters
    integer                            :: ns                          !! Maximum level of subdivisions
    complex(kind=real64)               :: h1(e%dme_n_gnodes,e%n_pnodes,2) !! h1 matrix
    complex(kind=real64)               :: h2(               e%n_pnodes,2) !! h2 matrix
    complex(kind=real64)               :: g1(e%dme_n_gnodes,e%n_snodes,2) !! g1 matrix
    complex(kind=real64)               :: g2(               e%n_snodes,2) !! g2 matrix
    ! Local
    real(kind=real64)        :: r(2)              ! Distance vector
    real(kind=real64)        :: rmin              ! Minimum distance between element and x_i
    real(kind=real64)        :: barxi(1)          ! Local coordinates of the nearest element point with respect to x_i
    real(kind=real64)        :: d                 ! Dimensionless distance
    integer                  :: delta             ! Control variable
    real(kind=real64)        :: xi_s(1,2)         ! Local coordinates of the element subdivision
    integer                  :: method            ! Method used when calculating the nearest element point
    integer                  :: gln_near          ! 1D Gauss-Legendre integration points required by the quasi-singular function
    integer                  :: gln               ! 1D Gauss-Legendre integration points used in the integration
    integer                  :: ps                ! Selected precalculated dataset
    integer                  :: i                 ! Counter
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
      call fbem_nearest_element_point_bem(2,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
      if (d.le.1.d-12) then
        delta=1
      else
        delta=0
      end if
    end if
    ! Integrate
    select case (delta)
      case (1)
        h2=(0.d0,0.d0)
        g2=(0.d0,0.d0)
        call fbem_bem_harpot2d_vsbie_int(e,reverse,barxi,p,h1,g1)
      case (0)
        ! Estimate the required integration rule
        gln_near=fbem_qs_n_estimation_standard(e%n,e%gtype,4,qsp,d,barxi)
        gln=max(e%gln_far,gln_near)
        ! Integrate using a conservative precalculated dataset
        if ((gln.le.e%ps_gln_max).and.(gln_near.gt.0)) then
          do i=1,e%n_ps
            if (e%ps_gln(i).ge.gln) then
              ps=i
              exit
            end if
          end do
          call fbem_bem_harpot2d_vsbie_ext_pre(ps,e,reverse,x_i,p,h1,h2,g1,g2)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_harpot2d_vsbie_ext_adp(e,reverse,xi_s,x_i,p,qsp,1,ns,h1,h2,g1,g2)
        end if
    end select
  end subroutine fbem_bem_harpot2d_vsbie_auto

  ! ================================================================================================================================

end module fbem_bem_harpot2d
