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

!! @author Jacob David Rodriguez Bordon (jacobdavid.rodriguezbordon@ulpgc.es)
!!
!! @version 2.0
!!
!! <b> This module implements the calculation of BEM integrals for 2D Laplace problems.</b>
module fbem_bem_stapot2d

  ! Fortran 2003 standard
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_quad_rules
  use fbem_geometry
  use fbem_telles_transformation
  use fbem_quasisingular_integration
  use fbem_bem_general

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all are private
  private

  ! ================================================================================================================================
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE)
  ! Free-term c
  public :: fbem_bem_pot2d_sbie_freeterm
  ! Exterior integration
  public :: fbem_bem_stapot2d_sbie_ext_pre
  public :: fbem_bem_stapot2d_sbie_ext_st
  public :: fbem_bem_stapot2d_sbie_ext_adp
  ! Interior integration
  public :: fbem_bem_stapot2d_sbie_int
  ! Automatic integration
  public :: fbem_bem_stapot2d_sbie_auto
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (HBIE)
  ! Exterior integration
  public :: fbem_bem_stapot2d_hbie_ext_pre
  public :: fbem_bem_stapot2d_hbie_ext_st
  public :: fbem_bem_stapot2d_hbie_ext_adp
  ! Interior integration
  public :: fbem_bem_stapot2d_hbie_int
  ! Automatic integration
  public :: fbem_bem_stapot2d_hbie_auto
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! VARIATION SINGULAR BOUNDARY INTEGRAL EQUATION (VSBIE)
  ! Free-term b
  public :: fbem_bem_pot2d_vsbie_freeterm
  ! Exterior integration
  public :: fbem_bem_stapot2d_vsbie_ext_pre
  public :: fbem_bem_stapot2d_vsbie_ext_st
  public :: fbem_bem_stapot2d_vsbie_ext_adp
  ! Interior integration
  public :: fbem_bem_stapot2d_vsbie_int
  ! Automatic integration
  public :: fbem_bem_stapot2d_vsbie_auto
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! VARIATION HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (VHBIE)
  ! Free-term b
!  public :: fbem_bem_pot2d_vhbie_freeterm
  ! Exterior integration
  public :: fbem_bem_stapot2d_vhbie_ext_pre
  public :: fbem_bem_stapot2d_vhbie_ext_st
  public :: fbem_bem_stapot2d_vhbie_ext_adp
  ! Interior integration
  public :: fbem_bem_stapot2d_vhbie_int
  ! Automatic integration
  public :: fbem_bem_stapot2d_vhbie_auto
  ! ================================================================================================================================

contains

  ! ================================================================================================================================
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE)

  subroutine fbem_bem_pot2d_sbie_freeterm(n_elements,n,tc,tol,c)
    implicit none
    ! I/O
    integer           :: n_elements  !! Number of elements (1 or 2)
    real(kind=real64) :: n(2,2)      !! Unit normals of elements at collocation point
    real(kind=real64) :: tc(2,2)     !! Unit tangents of elements at collocation point towards inside the element
    real(kind=real64) :: tol         !! Geometric tolerance [1.0e-12,1.0e-3] (default 1.0e-6)
    real(kind=real64) :: c           !! Free-term
    ! Local
    real(kind=real64) :: nsum(2)     ! Unit sum of normals
    real(kind=real64) :: local_tol   ! Local copy of geometric tolerance
    real(kind=real64) :: vectornorm  ! Norm of a vector
    real(kind=real64) :: alpha, beta ! Angles between tangents and unit sum of normals
    integer           :: i           ! Counter
    ! Check n_elements
    if ((n_elements.lt.1).or.(n_elements.gt.2)) then
      call fbem_error_message(output_unit,0,__FILE__,__LINE__,&
                              'n_elements must be 1 or 2 for free-term calculation.')
    end if
    ! Check geometric tolerance
    if ((tol.lt.1.0d-12).or.(tol.gt.1.0d-3)) then
      local_tol=1.0d-6
    else
      local_tol=tol
    end if
    ! Check that n and tc are orthogonal (using scalar product)
    do i=1,n_elements
      if (dabs(tc(1,i)*n(1,i)+tc(2,i)*n(2,i)).gt.local_tol) then
        call fbem_error_message(output_unit,0,__FILE__,__LINE__,&
                                'a pair of vectors n and tc given to the free-term routine are not orthogonal.')
      end if
    end do
    ! The calculation is different for 1 or 2 elements
    select case (n_elements)
      ! For 1 element, it is supposed that the boundary is smooth
      case (1)
        c=0.5d0
      ! For 2 elements, the freeterm is calculated as a vertex
      case (2)
        ! Add up normals and normalize it
        nsum(1)=n(1,1)+n(1,2)
        nsum(2)=n(2,1)+n(2,2)
        vectornorm=dsqrt(nsum(1)**2+nsum(2)**2)
        ! If the norm is too low, is because elements are almost anti-parallel.
        if (vectornorm.lt.local_tol) then
          call fbem_error_message(output_unit,0,__FILE__,__LINE__,&
                                  'when calculating free-term, it was found that two connected elements are almost parallel.')
        end if
        nsum(1)=nsum(1)/vectornorm
        nsum(2)=nsum(2)/vectornorm
        ! Angles between each tangent and unit normal sum
        alpha=dacos(tc(1,1)*nsum(1)+tc(2,1)*nsum(2))
        beta =dacos(tc(1,2)*nsum(1)+tc(2,2)*nsum(2))
        ! Free-term
        c=1.0d0-(alpha+beta)/(2.0d0*c_pi)
    end select
  end subroutine fbem_bem_pot2d_sbie_freeterm

  subroutine fbem_bem_stapot2d_sbie_ext_pre(ps,e,reverse,x_i,h,g)
    implicit none
    ! I/O
    integer                :: ps            !! Selected precalculated dataset
    type(fbem_bem_element) :: e             !! Element
    logical                :: reverse       !! Reverse element orientation
    real(kind=real64)      :: x_i(2)        !! Position vector of the collocation point
    real(kind=real64)      :: h(e%n_pnodes) !! h vector
    real(kind=real64)      :: g(e%n_snodes) !! g vector
    ! Local
    integer           :: kip                 ! Counter
    real(kind=real64) :: x(2), n(2)          ! Position and unit normal vectors at the integration point
    real(kind=real64) :: rv(2), r, dr1, logr ! Distance vector and its modulus, inverse and natural logarithm
    real(kind=real64) :: drdn                ! dr/dn
    real(kind=real64) :: pphijw(e%n_pnodes)  ! shape functions (primary variable) * jacobian * weight
    real(kind=real64) :: sphijw(e%n_snodes)  ! shape functions (secondary variable) * jacobian * weight
    ! Initialization
    h=0.d0
    g=0.d0
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      ! Components of the integrand
      x=e%ps_x(:,kip,ps)
      n=e%ps_n(:,kip,ps)
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      dr1=1.d0/r
      logr=log(r)
      drdn=dot_product(rv,n)*dr1
      pphijw=e%ps_pphijw(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      ! Add the integration point
      h=h+dr1*drdn*pphijw
      g=g-logr*sphijw
    end do
    ! Multiply by constants
    h=-h*c_1_2pi
    g= g*c_1_2pi
    ! Reverse element orientation
    if (reverse) h=-h
  end subroutine fbem_bem_stapot2d_sbie_ext_pre

  subroutine fbem_bem_stapot2d_sbie_ext_pre_exp(ps,e,reverse,x_i,hm,gm,he,ge)
    implicit none
    ! I/O
    integer                :: ps               !! Selected precalculated dataset
    type(fbem_bem_element) :: e                !! Element
    logical                :: reverse          !! Reverse element orientation
    real(kind=real64)      :: x_i(2)           !! Position vector of the collocation point
    real(kind=real64)      :: gm(e%n_snodes)   !! g vector
    real(kind=real64)      :: hm(e%n_snodes)   !! h vector
    real(kind=real64)      :: ge(e%n_snodes,2) !! expansion of g matrix
    real(kind=real64)      :: he(e%n_snodes,2) !! expansion of h matrix
    ! Local
    integer           :: kip                      ! Counter
    real(kind=real64) :: x(2), n(2)               ! Position and unit normal vectors at the integration point
    real(kind=real64) :: rv(2), r, dr1, dr2, logr ! Distance vector and its modulus, inverse and natural logarithm
    real(kind=real64) :: drdx(2)                  ! dr/dx
    real(kind=real64) :: drdn                     ! dr/dn
    real(kind=real64) :: pphijw(e%n_pnodes)       ! shape functions (primary variable) * jacobian * weight
    real(kind=real64) :: sphijw(e%n_snodes)       ! shape functions (secondary variable) * jacobian * weight
    ! Initialization
    gm=0.d0
    hm=0.d0
    ge=0.d0
    he=0.d0
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      ! Components of the integrand
      x=e%ps_x(:,kip,ps)
      n=e%ps_n(:,kip,ps)
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      dr1=1.d0/r
      dr2=dr1**2
      logr=log(r)
      drdx=rv*dr1
      drdn=dot_product(rv,n)*dr1
      pphijw=e%ps_pphijw(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      ! Add the integration point
      hm=hm+dr1*drdn*pphijw
      gm=gm-logr*sphijw
      he(:,1)=he(:,1)-dr2*(n(1)-2.d0*drdn*drdx(1))*pphijw
      he(:,2)=he(:,2)-dr2*(n(2)-2.d0*drdn*drdx(2))*pphijw
      ge(:,1)=ge(:,1)+dr1*drdx(1)*sphijw
      ge(:,2)=ge(:,2)+dr1*drdx(2)*sphijw
    end do
    ! Multiply by constants
    hm=-hm*c_1_2pi
    he=-he*c_1_2pi
    gm= gm*c_1_2pi
    ge= ge*c_1_2pi
    ! Reverse element orientation
    if (reverse) then
      hm=-hm
      he=-he
    end if
  end subroutine fbem_bem_stapot2d_sbie_ext_pre_exp

  subroutine fbem_bem_stapot2d_sbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,gln,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element) :: e             !! Element
    logical                :: reverse       !! Reverse element orientation
    real(kind=real64)      :: xi_s(1,2)     !! Subdivision [xi_s(1,1),xi_s(1,2)] of the element (xi space [-1,1])
    real(kind=real64)      :: x_i(2)        !! Position vector of the collocation point
    real(kind=real64)      :: barxip(1)     !! Nearest coordinate of the subdivision with respect to the collocation point (xip space [-1,1])
    real(kind=real64)      :: barr          !! Telles jacobian at barxip
    integer                :: gln           !! Number of Gauss-Legendre integration points (<=32)
    real(kind=real64)      :: h(e%n_pnodes) !! h vector
    real(kind=real64)      :: g(e%n_snodes) !! g vector
    ! Local
    integer                      :: kip, kphi            ! Counters
    real(kind=real64)            :: gamma                ! gamma coordinate [-1,1] (Telles space)
    real(kind=real64)            :: w                    ! Integration point weight
    type(fbem_telles_parameters) :: telles_parameters    ! Telles parameters
    real(kind=real64)            :: xip                  ! xip coordinate [-1,1] (subdivision space)
    real(kind=real64)            :: jt                   ! Telles jacobian: gamma -> xip
    real(kind=real64)            :: xi                   ! xi coordinate [-1,1] (element space)
    real(kind=real64)            :: js                   ! Subdivision jacobian: xip -> xi
    real(kind=real64)            :: x(2), T(2), N(2)     ! Position, tangent and normal vectors
    real(kind=real64)            :: jg                   ! Geometric jacobian: xi -> x
    real(kind=real64)            :: rv(2), r, dr1, logr  ! Distance vector and its modulus, inverse and natural logarithm
    real(kind=real64)            :: drdn                 ! dr/dn
    real(kind=real64)            :: aux(10)              ! Auxiliary variable
    real(kind=real64)            :: gphi(e%n_gnodes)     ! Geometrical shape functions
    real(kind=real64)            :: dgphidxi(e%n_gnodes) ! Geometrical shape functions derivatives
    real(kind=real64)            :: pphi(e%n_pnodes)     ! Functional (primary variable) shape functions
    real(kind=real64)            :: sphi(e%n_snodes)     ! Functional (secondary variable) shape functions
    real(kind=real64)            :: jw                   ! Jacobians * weight
    real(kind=real64)            :: pphijw(e%n_pnodes)   ! shape functions (primary variable) * jw
    real(kind=real64)            :: sphijw(e%n_snodes)   ! shape functions (secondary variable) * jw
    ! Initialization
    h=0.d0
    g=0.d0
    telles_parameters=fbem_telles11_calculate_parameters(barxip(1),barr)
    ! Numerical integration
    do kip=1,gl11_n(gln)
      ! gamma [-1,1]
      gamma=gl11_xi(kip,gln)
      w=gl11_w(kip,gln)
      ! gamma [-1,1] -> xip [-1,1]
      call fbem_telles_xi_and_jacobian(telles_parameters,gamma,xip,jt)
      ! xip [-1,1] -> xi [-1,1]
      xi=0.5d0*(1.d0-xip)*xi_s(1,1)+0.5d0*(1.d0+xip)*xi_s(1,2)
      js=0.5d0*(xi_s(1,2)-xi_s(1,1))
      ! xi [-1,1] -> x, T
      ! Geometrical shape functions and their derivatives
#     define etype e%gtype
#     define delta 0.0d0
#     define phi gphi
#     define dphidxi dgphidxi
#     include <phi_and_dphidxi_1d.rc>
#     undef etype
#     undef delta
#     undef phi
#     undef dphidxi
      ! Position vector x and tangent vector T
      x=0.d0
      T=0.d0
      do kphi=1,e%n_gnodes
        x=x+gphi(kphi)*e%x(:,kphi)
        T=T+dgphidxi(kphi)*e%x(:,kphi)
      end do
      ! Jacobian: xi [-1,1] -> x
      jg=sqrt(dot_product(T,T))
      ! Unit tangent and unit normal
      t=T/jg
      n(1)=t(2)
      n(2)=-t(1)
      ! Distance vector and its modulus, inverse and natural logarithm
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      dr1=1.d0/r
      logr=log(r)
      ! dr/dn
      drdn=dot_product(rv,n)*dr1
      ! Functional shape functions (primary variable)
#     define etype e%ptype
#     define delta e%ptype_delta
#     define phi pphi
#     include <phi_1d.rc>
#     undef etype
#     undef delta
#     undef phi
      ! Functional shape functions (secondary variable)
#     define etype e%stype
#     define delta e%stype_delta
#     define phi sphi
#     include <phi_1d.rc>
#     undef etype
#     undef delta
#     undef phi
      ! jacobians * weight
      jw=jg*js*jt*w
      ! Functional shape functions * jacobians * weights
      pphijw=pphi*jw
      sphijw=sphi*jw
      ! Add the integration point
      h=h+dr1*drdn*pphijw
      g=g-logr*sphijw
    end do
    ! Multiply by constants
    h=-h*c_1_2pi
    g= g*c_1_2pi
    ! Reverse element orientation
    if (reverse) h=-h
  end subroutine fbem_bem_stapot2d_sbie_ext_st

  recursive subroutine fbem_bem_stapot2d_sbie_ext_adp(e,reverse,xi_s,x_i,qsp,ks,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e             !! Element
    logical                  :: reverse       !! Reverse element orientation
    real(kind=real64)        :: xi_s(1,2)     !! Subdivision [xi_s(1,1),xi_s(1,2)] of the element (xi space [-1,1])
    real(kind=real64)        :: x_i(2)        !! Position vector of the collocation point
    type(fbem_qs_parameters) :: qsp           !! Quasi-singular integration parameters
    integer                  :: ks            !! Current level of subdivision
    integer                  :: ns            !! Maximum level of subdivision
    real(kind=real64)        :: h(e%n_pnodes) !! h vector
    real(kind=real64)        :: g(e%n_snodes) !! g vector
    ! Local
    integer           :: gln_near          ! Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer           :: gln               ! Gauss-Legendre integ. points required to integrate the whole integrand
    logical           :: subdivide         ! True if subdivision has to be performed
    real(kind=real64) :: barxi(1)          ! Nearest element coordinate with respect to collocation point
    real(kind=real64) :: barxip(1)         ! Nearest element subdivision local coordinate with respect to collocation point
    real(kind=real64) :: rmin              ! Minimum distance between collocation point and barxi on the element
    real(kind=real64) :: barr              ! Telles jacobian at barxi
    real(kind=real64) :: cl                ! Characteristic length
    real(kind=real64) :: d                 ! Normalized distance between collocation point and element subdivision
    integer           :: method            ! Method used in nearest point algorithm
    real(kind=real64) :: tmp_xi_s(1,2)     ! Subdivision
    real(kind=real64) :: x_s(2,e%n_gnodes) ! Coordinates of the element subdivision
    real(kind=real64) :: h_tmp(e%n_pnodes) ! h vector (temporary)
    real(kind=real64) :: g_tmp(e%n_snodes) ! g vector (temporary)
    ! Initialize
    if (ks.eq.1) then
      h=0.d0
      g=0.d0
      xi_s(1,1)=-1.d0
      xi_s(1,2)= 1.d0
      call fbem_nearest_element_point_bem(2,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
      barxip=barxi
    else
      call fbem_obtain_element_subdivision_coordinates(2,e%gtype,e%x,xi_s,x_s)
      cl=fbem_characteristic_length(2,e%gtype,x_s,1.d-6)
      call fbem_nearest_element_point_bem(2,e%gtype,x_s,cl,x_i,barxip,rmin,d,method)
    end if
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,2,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_stapot2d_sbie_ext_adp',ns,'maximum number of subdivisions reached')
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
      call fbem_bem_stapot2d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,qsp,ks+1,ns,h,g)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_stapot2d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,qsp,ks+1,ns,h,g)
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_stapot2d_sbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,gln,h_tmp,g_tmp)
      h=h+h_tmp
      g=g+g_tmp
    end if
  end subroutine fbem_bem_stapot2d_sbie_ext_adp

  subroutine fbem_bem_stapot2d_sbie_int(e,reverse,xi_i,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element) :: e             !! Element
    logical                :: reverse       !! Reverse element orientation
    real(kind=real64)      :: xi_i(1)       !! Position vector of the collocation point
    real(kind=real64)      :: h(e%n_pnodes) !! h vector
    real(kind=real64)      :: g(e%n_snodes) !! g vector
    ! Local
    integer           :: gln                  ! Number of Gauss-Legendre integration points (<=32)
    integer           :: kphi, kip            ! Counters
    real(kind=real64) :: x_i(2)               ! Real coordinates of collocation point
    integer           :: nsub                 ! Number of subdivision of the element
    integer           :: ksub                 ! Counter of subdivision
    real(kind=real64) :: xip                  ! Coordinate xip of subdivided element [0,1]
    real(kind=real64) :: w                    ! Weights of each integration point
    real(kind=real64) :: js                   ! Jacobian of the xi [xisub(1,:),xisub[2,:]] -> xip [0,1] transformation
    real(kind=real64) :: xi                   ! Coordinate xi
    real(kind=real64) :: xisub(2,2)           ! Coordinates of subdivisions in xi space
    real(kind=real64) :: x(2), T(2), N(2)     ! Position, tangent and normal vectors
    real(kind=real64) :: jg                   ! Geometric jacobian: xi -> x
    real(kind=real64) :: rv(2), r, dr1, logr  ! Distance vector and its modulus, inverse and natural logarithm
    real(kind=real64) :: drdn                 ! dr/dn
    real(kind=real64) :: aux(10)              ! Auxiliary variable
    real(kind=real64) :: gphi(e%n_gnodes)     ! Geometrical shape functions
    real(kind=real64) :: dgphidxi(e%n_gnodes) ! Geometrical shape functions derivatives
    real(kind=real64) :: pphi(e%n_pnodes)     ! Functional (primary variable) shape functions
    real(kind=real64) :: sphi(e%n_snodes)     ! Functional (secondary variable) shape functions
    real(kind=real64) :: sphi_i(e%n_snodes)   ! Functional (secondary variable) shape functions at the collocation point
    real(kind=real64) :: jw                   ! Jacobians * weight
    real(kind=real64) :: Ilogr                ! Int log(r) dGamma
    ! Initialization
    gln=32
    h=0.d0
    g=0.d0
    ! Setup the subdivisions for transformation xip -> xi: [0,1] -> [xisub_1,xisub_2]
    if (fbem_check_xi_vertex(xi_i(1))) then
      nsub=1
      xisub(1,1)=-1.0d0
      xisub(2,1)= 1.0d0
    else
      nsub=2
      xisub(1,1)=-1.0d0
      xisub(2,1)=xi_i(1)
      xisub(1,2)=xi_i(1)
      xisub(2,2)=1.0d0
    end if
    ! Calculate x_i
#   define etype e%gtype
#   define delta 0.0d0
#   define xi xi_i(1)
#   define phi gphi
#   include <phi_1d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
    x_i=0.d0
    do kphi=1,e%gtype
      x_i=x_i+gphi(kphi)*e%x(:,kphi)
    end do
    ! Functional shape functions (secondary variable) at the collocation point
#   define etype e%stype
#   define delta e%stype_delta
#   define xi xi_i(1)
#   define phi sphi_i
#   include <phi_1d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
    ! Numerical integration
    do ksub=1,nsub
      do kip=1,gl01_n(gln)
        ! xip [0,1]
        xip=gl01_xi(kip,gln)
        w=gl01_w(kip,gln)
        ! xip [0,1] -> xi [-1,1]
        xi=(xisub(2,ksub)-xisub(1,ksub))*xip+xisub(1,ksub)
        js=xisub(2,ksub)-xisub(1,ksub)
        ! xi [-1,1] -> x, T
        ! Geometrical shape functions and their derivatives
#       define etype e%gtype
#       define delta 0.0d0
#       define phi gphi
#       define dphidxi dgphidxi
#       include <phi_and_dphidxi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
#       undef dphidxi
        ! Position vector x and tangent vector T
        x=0.d0
        T=0.d0
        do kphi=1,e%n_gnodes
          x=x+gphi(kphi)*e%x(:,kphi)
          T=T+dgphidxi(kphi)*e%x(:,kphi)
        end do
        ! Jacobian: xi [-1,1] -> x
        jg=sqrt(dot_product(T,T))
        ! Unit tangent and unit normal
        t=T/jg
        n(1)=t(2)
        n(2)=-t(1)
        ! Distance vector and its modulus, inverse and natural logarithm
        rv=x-x_i
        r=sqrt(dot_product(rv,rv))
        dr1=1.d0/r
        logr=log(r)
        ! dr/dn
        drdn=dot_product(rv,n)*dr1
        ! Functional shape functions (primary variable) at xi
#       define etype e%ptype
#       define delta e%ptype_delta
#       define phi pphi
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Functional shape functions (secondary variable) at xi
#       define etype e%stype
#       define delta e%stype_delta
#       define phi sphi
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Jacobians * weight
        jw=jg*js*w
        ! Add the integration point
        h=h+dr1*drdn*pphi*jw
        g=g-logr*(sphi-sphi_i)*jw
      end do
    end do
    ! Add log(r) integrals
    Ilogr=fbem_bem_int_logr(2,e%gtype,e%x,xi_i)
    g=g-sphi_i*Ilogr
    ! Multiply by constants
    h=-h*c_1_2pi
    g= g*c_1_2pi
    ! Reverse element orientation
    if (reverse) h=-h
  end subroutine fbem_bem_stapot2d_sbie_int

  subroutine fbem_bem_stapot2d_sbie_auto(e,reverse,x_i,qsp,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e             !! Element
    logical                  :: reverse       !! Reverse element orientation
    real(kind=real64)        :: x_i(2)        !! Position vector of the collocation point
    type(fbem_qs_parameters) :: qsp           !! Quasi-singular integration parameters
    integer                  :: ns            !! Maximum level of subdivisions
    real(kind=real64)        :: h(e%n_pnodes) !! h vector
    real(kind=real64)        :: g(e%n_snodes) !! g vector
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
        call fbem_bem_stapot2d_sbie_int(e,reverse,barxi,h,g)
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
          call fbem_bem_stapot2d_sbie_ext_pre(ps,e,reverse,x_i,h,g)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_stapot2d_sbie_ext_adp(e,reverse,xi_s,x_i,qsp,1,ns,h,g)
        end if
    end select
  end subroutine fbem_bem_stapot2d_sbie_auto

  ! ================================================================================================================================

  ! ================================================================================================================================
  ! HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (HBIE)

  subroutine fbem_bem_stapot2d_hbie_ext_pre(ps,e,reverse,x_i,n_i,m,l)
    implicit none
    ! I/O
    integer                :: ps            !! Selected precalculated dataset
    type(fbem_bem_element) :: e             !! Element
    logical                :: reverse       !! Reverse element orientation
    real(kind=real64)      :: x_i(2)        !! Position vector of the collocation point
    real(kind=real64)      :: n_i(2)        !! Unit normal vector of the collocation point
    real(kind=real64)      :: m(e%n_pnodes) !! m vector
    real(kind=real64)      :: l(e%n_snodes) !! l vector
    ! Local
    integer           :: kip                ! Counter
    real(kind=real64) :: x(2), n(2)         ! Position and unit normal vectors at the integration point
    real(kind=real64) :: rv(2), r, dr1, dr2 ! Distance vector and its modulus, 1/r and 1/r^2
    real(kind=real64) :: drdn               ! dr/dn
    real(kind=real64) :: drdni              ! dr/dn^i
    real(kind=real64) :: n_dot_ni           ! n·n^i
    real(kind=real64) :: pphijw(e%n_pnodes) ! shape functions (primary variable) * jacobian * weight
    real(kind=real64) :: sphijw(e%n_snodes) ! shape functions (secondary variable) * jacobian * weight
    ! Initialization
    m=0.d0
    l=0.d0
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      ! Components of the integrand
      x=e%ps_x(:,kip,ps)
      n=e%ps_n(:,kip,ps)
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      dr1=1.d0/r
      dr2=dr1**2
      drdn=dot_product(rv,n)*dr1
      drdni=-dot_product(rv,n_i)*dr1
      n_dot_ni=dot_product(n,n_i)
      pphijw=e%ps_pphijw(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      ! Add the integration point
      m=m+dr2*(2.d0*drdn*drdni+n_dot_ni)*pphijw
      l=l+dr1*drdni*sphijw
    end do
    ! Multiply by constants
    m= m*c_1_2pi
    l=-l*c_1_2pi
    ! Reverse element orientation
    if (reverse) m=-m
  end subroutine fbem_bem_stapot2d_hbie_ext_pre

  subroutine fbem_bem_stapot2d_hbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,gln,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element) :: e             !! Element
    logical                :: reverse       !! Reverse element orientation
    real(kind=real64)      :: xi_s(1,2)     !! Subdivision [xi_s(1,1),xi_s(1,2)] of the element (xi space [-1,1])
    real(kind=real64)      :: x_i(2)        !! Position vector of the collocation point
    real(kind=real64)      :: n_i(2)        !! Unit normal vector of the collocation point
    real(kind=real64)      :: barxip(1)     !! Nearest coordinate of the subdivision with respect to the collocation point (xip space [-1,1])
    real(kind=real64)      :: barr          !! Telles jacobian at barxip
    integer                :: gln           !! Number of Gauss-Legendre integration points (<=32)
    real(kind=real64)      :: m(e%n_pnodes) !! m vector
    real(kind=real64)      :: l(e%n_snodes) !! l vector
    ! Local
    integer                      :: kip, kphi            ! Counters
    real(kind=real64)            :: gamma                ! gamma coordinate [-1,1] (Telles space)
    real(kind=real64)            :: w                    ! Integration point weight
    type(fbem_telles_parameters) :: telles_parameters    ! Telles parameters
    real(kind=real64)            :: xip                  ! xip coordinate [-1,1] (subdivision space)
    real(kind=real64)            :: jt                   ! Telles jacobian: gamma -> xip
    real(kind=real64)            :: xi                   ! xi coordinate [-1,1] (element space)
    real(kind=real64)            :: js                   ! Subdivision jacobian: xip -> xi
    real(kind=real64)            :: x(2), T(2), N(2)     ! Position, tangent and normal vectors
    real(kind=real64)            :: jg                   ! Geometric jacobian: xi -> x
    real(kind=real64)            :: rv(2), r, dr1, dr2   ! Distance vector and its modulus, 1/r and 1/r^2
    real(kind=real64)            :: drdn                 ! dr/dn
    real(kind=real64)            :: drdni                ! dr/dn^i
    real(kind=real64)            :: n_dot_ni             ! n·n^i
    real(kind=real64)            :: aux(10)              ! Auxiliary variable
    real(kind=real64)            :: gphi(e%n_gnodes)     ! Geometrical shape functions
    real(kind=real64)            :: dgphidxi(e%n_gnodes) ! Geometrical shape functions derivatives
    real(kind=real64)            :: pphi(e%n_pnodes)     ! Functional (primary variable) shape functions
    real(kind=real64)            :: sphi(e%n_snodes)     ! Functional (secondary variable) shape functions
    real(kind=real64)            :: jw                   ! Jacobians * weight
    real(kind=real64)            :: pphijw(e%n_pnodes)   ! shape functions (primary variable) * jw
    real(kind=real64)            :: sphijw(e%n_snodes)   ! shape functions (secondary variable) * jw
    ! Initialization
    m=0.d0
    l=0.d0
    telles_parameters=fbem_telles11_calculate_parameters(barxip(1),barr)
    ! Loop through gamma coordinate
    do kip=1,gl11_n(gln)
      ! gamma [-1,1]
      gamma=gl11_xi(kip,gln)
      w=gl11_w(kip,gln)
      ! gamma [-1,1] -> xip [-1,1]
      call fbem_telles_xi_and_jacobian(telles_parameters,gamma,xip,jt)
      ! xip [-1,1] -> xi [-1,1]
      xi=0.5d0*(1.d0-xip)*xi_s(1,1)+0.5d0*(1.d0+xip)*xi_s(1,2)
      js=0.5d0*(xi_s(1,2)-xi_s(1,1))
      ! xi [-1,1] -> x, T
      ! Geometrical shape functions and their derivatives
#     define etype e%gtype
#     define delta 0.0d0
#     define phi gphi
#     define dphidxi dgphidxi
#     include <phi_and_dphidxi_1d.rc>
#     undef etype
#     undef delta
#     undef phi
#     undef dphidxi
      ! Position vector x and tangent vector T
      x=0.d0
      T=0.d0
      do kphi=1,e%n_gnodes
        x=x+gphi(kphi)*e%x(:,kphi)
        T=T+dgphidxi(kphi)*e%x(:,kphi)
      end do
      ! Jacobian: xi [-1,1] -> x
      jg=sqrt(dot_product(T,T))
      ! Unit tangent and unit normal
      t=T/jg
      n(1)=t(2)
      n(2)=-t(1)
      ! Distance vector and its modulus, 1/r and 1/r^2
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      dr1=1.d0/r
      dr2=dr1**2
      ! dr/dn
      drdn=dot_product(rv,n)*dr1
      ! dr/dn_i
      drdni=-dot_product(rv,n_i)*dr1
      ! n·n_i
      n_dot_ni=dot_product(n,n_i)
      ! Functional shape functions (primary variable)
#     define etype e%ptype
#     define delta e%ptype_delta
#     define phi pphi
#     include <phi_1d.rc>
#     undef etype
#     undef delta
#     undef phi
      ! Functional shape functions (secondary variable)
#     define etype e%stype
#     define delta e%stype_delta
#     define phi sphi
#     include <phi_1d.rc>
#     undef etype
#     undef delta
#     undef phi
      ! jacobians * weight
      jw=jg*js*jt*w
      ! Functional shape functions * jacobians * weights
      pphijw=pphi*jw
      sphijw=sphi*jw
      ! Add the integration point
      m=m+dr2*(2.d0*drdn*drdni+n_dot_ni)*pphijw
      l=l+dr1*drdni*sphijw
    end do ! Loop through gamma coordinate
    ! Multiply m constants
    m= m*c_1_2pi
    l=-l*c_1_2pi
    ! Reverse element orientation
    if (reverse) m=-m
  end subroutine fbem_bem_stapot2d_hbie_ext_st

  recursive subroutine fbem_bem_stapot2d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,qsp,ks,ns,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e             !! Element
    logical                  :: reverse       !! Reverse element orientation
    real(kind=real64)        :: xi_s(1,2)     !! Subdivision [xi_s(1,1),xi_s(1,2)] of the element (xi space [-1,1])
    real(kind=real64)        :: x_i(2)        !! Position vector of the collocation point
    real(kind=real64)        :: n_i(2)        !! Unit normal vector of the collocation point
    type(fbem_qs_parameters) :: qsp           !! Quasi-singular integration parameters
    integer                  :: ks            !! Current level of subdivision
    integer                  :: ns            !! Maximum level of subdivision
    real(kind=real64)        :: m(e%n_pnodes) !! m vector
    real(kind=real64)        :: l(e%n_snodes) !! l vector
    ! Local
    integer           :: gln_near          ! Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer           :: gln               ! Gauss-Legendre integ. points required to integrate the whole integrand
    logical           :: subdivide         ! True if subdivision has to be performed
    real(kind=real64) :: barxi(1)          ! Nearest element coordinate with respect to collocation point
    real(kind=real64) :: barxip(1)         ! Nearest element subdivision local coordinate with respect to collocation point
    real(kind=real64) :: rmin              ! Minimum distance between collocation point and barxi on the element
    real(kind=real64) :: barr              ! Telles jacobian at barxi
    real(kind=real64) :: cl                ! Characteristic length
    real(kind=real64) :: d                 ! Normalized distance between collocation point and element subdivision
    integer           :: method            ! Method used in nearest point algorithm
    real(kind=real64) :: tmp_xi_s(1,2)     ! Subdivision
    real(kind=real64) :: x_s(2,e%n_gnodes) ! Coordinates of the element subdivision
    real(kind=real64) :: m_tmp(e%n_pnodes) ! m vector (temporary)
    real(kind=real64) :: l_tmp(e%n_snodes) ! l vector (temporary)
    ! Initialize
    if (ks.eq.1) then
      m=0.d0
      l=0.d0
      xi_s(1,1)=-1.d0
      xi_s(1,2)= 1.d0
      call fbem_nearest_element_point_bem(2,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
      barxip=barxi
    else
      call fbem_obtain_element_subdivision_coordinates(2,e%gtype,e%x,xi_s,x_s)
      cl=fbem_characteristic_length(2,e%gtype,x_s,1.d-6)
      call fbem_nearest_element_point_bem(2,e%gtype,x_s,cl,x_i,barxip,rmin,d,method)
    end if
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,4,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_stapot2d_hbie_ext_adp',ns,'maximum number of subdivisions reached')
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
      call fbem_bem_stapot2d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,qsp,ks+1,ns,m,l)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_stapot2d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,qsp,ks+1,ns,m,l)
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_stapot2d_hbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,gln,m_tmp,l_tmp)
      m=m+m_tmp
      l=l+l_tmp
    end if
  end subroutine fbem_bem_stapot2d_hbie_ext_adp

  subroutine fbem_bem_stapot2d_hbie_int(e,reverse,xi_i,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element) :: e             !! Element
    logical                :: reverse       !! Reverse element orientation
    real(kind=real64)      :: xi_i(1)       !! Position vector of the collocation point
    real(kind=real64)      :: m(e%n_pnodes) !! m vector
    real(kind=real64)      :: l(e%n_snodes) !! l vector
    ! Local
    integer           :: gln                    ! Number of Gauss-Legendre integration points (<=32)
    integer           :: kphi, kip              ! Counters
    integer           :: nsub                   ! Number of subdivision of the element
    integer           :: ksub                   ! Counter of subdivision
    real(kind=real64) :: xip                    ! Coordinate xip of subdivided element [0,1]
    real(kind=real64) :: w                      ! Weights of each integration point
    real(kind=real64) :: js                     ! Jacobian of the xi [xisub(1,:),xisub[2,:]] -> xip [0,1] transformation
    real(kind=real64) :: xi                     ! Coordinate xi
    real(kind=real64) :: xisub(2,2)             ! Coordinates of subdivisions in xi space
    real(kind=real64) :: x(2), T(2), N(2)       ! Position, tangent and normal vectors
    real(kind=real64) :: x_i(2), T_i(2), N_i(2) ! Position, tangent and normal vectors at the collocation point
    real(kind=real64) :: jg                     ! Geometric jacobian: xi -> x
    real(kind=real64) :: jgi                    ! Geometric jacobian: xi -> x at the collocation point
    real(kind=real64) :: rv(2), r, dr1, dr2     ! Distance vector and its modulus, 1/r and 1/r^2
    real(kind=real64) :: drdn                   ! dr/dn
    real(kind=real64) :: drdni                  ! dr/dn^i
    real(kind=real64) :: drdt                   ! dr/dt
    real(kind=real64) :: drdti                  ! dr/dt^i
    real(kind=real64) :: n_dot_ni               ! n·n^i
    real(kind=real64) :: aux(10)                ! Auxiliary variable
    real(kind=real64) :: gphi(e%n_gnodes)       ! Geometrical shape functions
    real(kind=real64) :: dgphidxi(e%n_gnodes)   ! Geometrical shape functions derivatives
    real(kind=real64) :: pphi(e%n_pnodes)       ! Functional (primary variable) shape functions
    real(kind=real64) :: pphi_i(e%n_pnodes)     ! Functional (primary variable) shape functions at the collocation point
    real(kind=real64) :: dpphidxi_i(e%n_pnodes) ! Functional (primary variable) shape functions derivatives at the collocation point
    real(kind=real64) :: sphi(e%n_snodes)       ! Functional (secondary variable) shape functions
    real(kind=real64) :: jw                     ! Jacobians * weight
    real(kind=real64) :: ra, rb                 ! Distance from collocation point to element vertices
    ! Initialization
    gln=32
    m=0.d0
    l=0.d0
    ! Setup the subdivisions for transformation xip -> xi: [0,1] -> [xisub_1,xisub_2]
    if (fbem_check_xi_vertex(xi_i(1))) then
      call fbem_error_message(output_unit,0,'fbem_bem_stapot2d_hbie_int',0,'the HBIE does not allow vertex collocation.')
    end if
    nsub=2
    xisub(1,1)=-1.0d0
    xisub(2,1)=xi_i(1)
    xisub(1,2)=xi_i(1)
    xisub(2,2)=1.0d0
    ! Calculate x_i, t_i, n_i and jgi
#   define etype e%gtype
#   define delta 0.0d0
#   define xi xi_i(1)
#   define phi gphi
#   define dphidxi dgphidxi
#   include <phi_and_dphidxi_1d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
#   undef dphidxi
    x_i=0.d0
    T_i=0.d0
    do kphi=1,e%n_gnodes
      x_i=x_i+gphi(kphi)*e%x(:,kphi)
      T_i=T_i+dgphidxi(kphi)*e%x(:,kphi)
    end do
    ! Jacobian: xi [-1,1] -> x at the collocation point
    jgi=sqrt(dot_product(T_i,T_i))
    ! Unit tangent and unit normal at the collocation point
    t_i=T_i/jgi
    n_i(1)=t_i(2)
    n_i(2)=-t_i(1)
    ! Functional (primary variable) shape functions and its first derivatives at the collocation point
#   define etype e%ptype
#   define delta e%ptype_delta
#   define xi xi_i(1)
#   define phi pphi_i
#   define dphidxi dpphidxi_i
#   include <phi_and_dphidxi_1d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
#   undef dphidxi
    ! Numerical integration
    do ksub=1,nsub
      if (ksub.eq.1) then
        drdti=-1.d0
      else
        drdti= 1.d0
      end if
      do kip=1,gl01_n(gln)
        ! xip [0,1]
        xip=gl01_xi(kip,gln)
        w=gl01_w(kip,gln)
        ! xip [0,1] -> xi [-1,1]
        xi=(xisub(2,ksub)-xisub(1,ksub))*xip+xisub(1,ksub)
        js=xisub(2,ksub)-xisub(1,ksub)
        ! xi [-1,1] -> x, T
        ! Geometrical shape functions and their derivatives
#       define etype e%gtype
#       define delta 0.0d0
#       define phi gphi
#       define dphidxi dgphidxi
#       include <phi_and_dphidxi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
#       undef dphidxi
        ! Position vector x and tangent vector T
        x=0.d0
        T=0.d0
        do kphi=1,e%n_gnodes
          x=x+gphi(kphi)*e%x(:,kphi)
          T=T+dgphidxi(kphi)*e%x(:,kphi)
        end do
        ! Jacobian: xi [-1,1] -> x
        jg=sqrt(dot_product(T,T))
        ! Unit tangent and unit normal
        t=T/jg
        n(1)=t(2)
        n(2)=-t(1)
        ! Distance vector and its modulus, inverse and natural logarithm
        rv=x-x_i
        r=sqrt(dot_product(rv,rv))
        dr1=1.d0/r
        dr2=dr1**2
        ! dr/dn
        drdn=dot_product(rv,n)*dr1
        ! dr/dn_i
        drdni=-dot_product(rv,n_i)*dr1
        ! dr/dt
        drdt=dot_product(rv,t)*dr1
        ! n dot n_i
        n_dot_ni=dot_product(n,n_i)
        ! Functional shape functions (primary variable) at xi
#       define etype e%ptype
#       define delta e%ptype_delta
#       define phi pphi
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Functional shape functions (secondary variable) at xi
#       define etype e%stype
#       define delta e%stype_delta
#       define phi sphi
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Jacobians * weight
        jw=jg*js*w
        ! Add integration points
        m=m+2.d0*dr2*drdn*drdni*pphi*jw
        m=m+dr2*(n_dot_ni-abs(drdt))*pphi*jw
        m=m+dr2*abs(drdt)*(pphi-pphi_i-drdti*dpphidxi_i/jgi*r)*jw
        l=l+dr1*drdni*sphi*jw
      end do
    end do
    ! Add analytical integrals
    ra=sqrt((e%x(1,1)-x_i(1))**2+(e%x(2,1)-x_i(2))**2)
    rb=sqrt((e%x(1,2)-x_i(1))**2+(e%x(2,2)-x_i(2))**2)
    m=m-pphi_i*(1.0d0/ra+1.0d0/rb)
    m=m+dpphidxi_i/jgi*(log(rb)-log(ra))
    ! Multiply by constants
    m= m*c_1_2pi
    l=-l*c_1_2pi
    ! Reverse element orientation
    if (reverse) l=-l
  end subroutine fbem_bem_stapot2d_hbie_int

  subroutine fbem_bem_stapot2d_hbie_auto(e,reverse,x_i,n_i,qsp,ns,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e             !! Element
    logical                  :: reverse       !! Reverse element orientation
    real(kind=real64)        :: x_i(2)        !! Position vector of the collocation point
    real(kind=real64)        :: n_i(2)        !! Unit normal vector of the collocation point
    type(fbem_qs_parameters) :: qsp           !! Quasi-singular integration parameters
    integer                  :: ns            !! Maximum level of subdivisions
    real(kind=real64)        :: m(e%n_pnodes) !! m vector
    real(kind=real64)        :: l(e%n_snodes) !! l vector
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
        call fbem_bem_stapot2d_hbie_int(e,reverse,barxi,m,l)
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
          call fbem_bem_stapot2d_hbie_ext_pre(ps,e,reverse,x_i,n_i,m,l)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_stapot2d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,qsp,1,ns,m,l)
        end if
    end select
  end subroutine fbem_bem_stapot2d_hbie_auto

  ! ================================================================================================================================

  ! ================================================================================================================================
  ! VARIATION SINGULAR BOUNDARY INTEGRAL EQUATION (VSBIE)

  subroutine fbem_bem_pot2d_vsbie_freeterm(n_elements,n,tc,tol,b)
    implicit none
    ! I/O
    integer           :: n_elements  !! Number of elements (1 or 2)
    real(kind=real64) :: n(2,2)      !! Unit normals of elements at collocation point
    real(kind=real64) :: tc(2,2)     !! Unit tangents of elements at collocation point towards inside the element
    real(kind=real64) :: tol         !! Geometric tolerance [1.0e-12,1.0e-3] (default 1.0e-6)
    real(kind=real64) :: b(2,2)      !! Free-term
    ! Local
    real(kind=real64) :: local_n(2,2)  ! Local copy of n
    real(kind=real64) :: local_tc(2,2) ! Local copy of tc
    real(kind=real64) :: local_tol     ! Local copy of geometric tolerance
    real(kind=real64) :: nsum(2)       ! Unit sum of normals
    real(kind=real64) :: vectornorm    ! Norm of a vector
    real(kind=real64) :: theta(2)      ! Angles of unit tangents
    integer           :: i             ! Counter
    ! Check n_elements
    if ((n_elements.lt.1).or.(n_elements.gt.2)) then
      call fbem_error_message(output_unit,0,__FILE__,__LINE__,&
                              'n_elements must be 1 or 2 for free-term calculation.')
    end if
    ! Check geometric tolerance
    if ((tol.lt.1.0d-12).or.(tol.gt.1.0d-3)) then
      local_tol=1.0d-6
    else
      local_tol=tol
    end if
    ! Check that n and tc are orthogonal (using scalar product)
    do i=1,n_elements
      if (abs(tc(1,i)*n(1,i)+tc(2,i)*n(2,i)).gt.local_tol) then
        call fbem_error_message(output_unit,0,__FILE__,__LINE__,&
                                'a pair of vectors n and tc given to the free-term routine are not orthogonal.')
      end if
    end do
    ! The calculation is different for 1 or 2 elements
    select case (n_elements)
      ! For 1 element, it is supposed that the boundary is smooth
      case (1)
        b=0.d0
      ! For 2 elements, the freeterm is calculated as a vertex
      case (2)
        !
        ! Know who is the element 1 and 2.
        !
        ! Add up normals and normalize it
        nsum(1)=n(1,1)+n(1,2)
        nsum(2)=n(2,1)+n(2,2)
        vectornorm=sqrt(nsum(1)**2+nsum(2)**2)
        ! If the norm is too low, is because elements are almost anti-parallel.
        if (vectornorm.lt.local_tol) then
          call fbem_error_message(output_unit,0,__FILE__,__LINE__,&
                                  'when calculating free-term, it was found that two connected elements are almost parallel.')
        end if
        ! The element 1 has the cross product: t x n > 0. This fact can be used to
        if ((tc(1,1)*n(2,1)-tc(2,1)*n(1,1)).gt.0) then
          local_n(1,1) =n(1,1)
          local_n(2,1) =n(2,1)
          local_tc(1,1)=tc(1,1)
          local_tc(2,1)=tc(2,1)
          local_n(1,2) =n(1,2)
          local_n(2,2) =n(2,2)
          local_tc(1,2)=tc(1,2)
          local_tc(2,2)=tc(2,2)
        else
          local_n(1,1) =n(1,2)
          local_n(2,1) =n(2,2)
          local_tc(1,1)=tc(1,2)
          local_tc(2,1)=tc(2,2)
          local_n(1,2) =n(1,1)
          local_n(2,2) =n(2,1)
          local_tc(1,2)=tc(1,1)
          local_tc(2,2)=tc(2,1)
        end if
        ! Note: if t x nsum = 0, no matter who is the first and the second.
        !
        ! Calculate angles
        !
        ! Angles theta in (-pi,pi]
        theta(1)=atan2(local_tc(2,1),local_tc(1,1))
        theta(2)=atan2(local_tc(2,2),local_tc(1,2))
        ! Convert to (0,2*pi)
        if (theta(1).lt.0.0d0) theta(1)=theta(1)+c_2pi
        if (theta(2).lt.0.0d0) theta(2)=theta(2)+c_2pi
        ! Use theta2>theta1
        if (theta(1).gt.theta(2)) theta(2)=theta(2)+c_2pi
        ! Calculations
        b(1,1)=sin(2.d0*theta(2))-sin(2.d0*theta(1))
        b(1,2)=cos(2.d0*theta(1))-cos(2.d0*theta(2))
        b(2,1)=b(1,2)
        b(2,2)=-b(1,1)
        b=c_1_4pi*b
    end select
  end subroutine fbem_bem_pot2d_vsbie_freeterm

  subroutine fbem_bem_stapot2d_vsbie_ext_pre(ps,e,reverse,x_i,h1,h2,g1,g2)
    implicit none
    ! I/O
    integer                :: ps                              !! Selected precalculated dataset
    type(fbem_bem_element) :: e                               !! Element
    logical                :: reverse                         !! Reverse normal vector
    real(kind=real64)      :: x_i(2)                          !! Collocation point position vector
    real(kind=real64)      :: h1(e%dme_n_gnodes,e%n_pnodes,2) !! h1 matrix
    real(kind=real64)      :: h2(               e%n_pnodes,2) !! h2 matrix
    real(kind=real64)      :: g1(e%dme_n_gnodes,e%n_snodes,2) !! g1 matrix
    real(kind=real64)      :: g2(               e%n_snodes,2) !! g2 matrix
    ! Local
    integer           :: im, iq                        ! Counter variables
    integer           :: kip                           ! Counter variable for integration points loop
    real(kind=real64) :: x(2)                          ! Position vector at integration point
    real(kind=real64) :: t(2)                          ! Unit tangent vector at integration point
    real(kind=real64) :: n(2)                          ! Unit normal vector at integration point
    real(kind=real64) :: pphijw(e%n_pnodes)            ! phi^p * jacobian * weight at integration point
    real(kind=real64) :: sphijw(e%n_snodes)            ! phi^s at integration point
    real(kind=real64) :: rv(2), r, dr1, dr2, logr      ! Distance vector module, 1/r^n and log(1/r)
    real(kind=real64) :: drdx(2)                       ! Distance vector derivatives with respect to x_k
    real(kind=real64) :: drdn                          ! Partial derivative of r respect to unit normal
    real(kind=real64) :: drdt                          ! Partial derivative of r respect to unit tangent
    real(kind=real64) :: dme_gphi(e%dme_n_gnodes)      ! phi^{dme-g} at integration point
    real(kind=real64) :: dme_wqj(e%dme_n_gnodes,2)     ! dphi{dme-g}/dx_j at integration point
    real(kind=real64) :: dme_vq(e%dme_n_gnodes)        ! vq=w_{j,q}·t_j at integration point
    real(kind=real64) :: fs_q, fs_qt, fs_dpdm, fs_dqdm ! Fundamental solutions values
    ! Initialization
    h1=0.d0
    h2=0.d0
    g1=0.d0
    g2=0.d0
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
      dr2=dr1**2
      logr=log(r)
      drdx=rv*dr1
      drdn=dot_product(drdx,n)
      drdt=dot_product(drdx,t)
      ! Add the integration point
      fs_q=dr1*drdn
      fs_qt=dr1*drdt
      do im=1,2
        fs_dpdm=dr1*drdx(im)
        fs_dqdm=dr2*(-2.d0*drdx(im)*drdn+n(im))
        do iq=1,e%dme_n_gnodes
          h1(iq,:,im)=h1(iq,:,im)+fs_q*dme_vq(iq)*t(im)*pphijw(:)
          h1(iq,:,im)=h1(iq,:,im)+fs_dqdm*dme_gphi(iq)*pphijw(:)
          h1(iq,:,im)=h1(iq,:,im)-fs_qt*dme_vq(iq)*n(im)*pphijw(:)
          g1(iq,:,im)=g1(iq,:,im)+logr*dme_vq(iq)*t(im)*sphijw(:)
          g1(iq,:,im)=g1(iq,:,im)+fs_dpdm*dme_gphi(iq)*sphijw(:)
        end do
        h2(:,im)=h2(:,im)+fs_dqdm*pphijw(:)
        g2(:,im)=g2(:,im)+fs_dpdm*sphijw(:)
      end do
    end do
    ! Multiply by constants
    h1=-c_1_2pi*h1
    h2=-c_1_2pi*h2
    g1=-c_1_2pi*g1
    g2=-c_1_2pi*g2
    ! Reverse element orientation
    if (reverse) then
      h1=-h1
      h2=-h2
    end if
  end subroutine fbem_bem_stapot2d_vsbie_ext_pre

  subroutine fbem_bem_stapot2d_vsbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,gln,h1,h2,g1,g2)
    implicit none
    ! I/O
    type(fbem_bem_element) :: e                               !! Integration element
    logical                :: reverse                         !! Reverse normal vector
    real(kind=real64)      :: xi_s(1,2)                       !! Coordinates of the subdivision of the parent element (xi space [-1,1])
    real(kind=real64)      :: x_i(2)                          !! Collocation point position vector
    real(kind=real64)      :: barxip(1)                       !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)      :: barr                            !! Telles jacobian at barxip
    integer                :: gln                             !! 1D Gauss-Legendre number of integration points (<=32)
    real(kind=real64)      :: h1(e%dme_n_gnodes,e%n_pnodes,2) !! h1 matrix
    real(kind=real64)      :: h2(               e%n_pnodes,2) !! h2 matrix
    real(kind=real64)      :: g1(e%dme_n_gnodes,e%n_snodes,2) !! g1 matrix
    real(kind=real64)      :: g2(               e%n_snodes,2) !! g2 matrix
    ! Local
    integer                      :: im, iq                        ! Counter variables
    integer                      :: kphi                          ! Counter variable for shape functions loops
    integer                      :: kip                           ! Counter variable of integration points
    real(kind=real64)            :: gamma                         ! Coordinate gamma (Telles transformation space [-1,1])
    real(kind=real64)            :: w                             ! Weights of an integration point
    type(fbem_telles_parameters) :: telles_parameters             ! Telles parameters
    real(kind=real64)            :: jt                            ! Telles jacobian: xip->gamma
    real(kind=real64)            :: xip                           ! Coordinate xip (subdivision space [-1,1])
    real(kind=real64)            :: js                            ! Subdivision jacobian: xi->xip
    real(kind=real64)            :: xi                            ! Coordinate xi [xi_s(1,1),xi_s(1,2)]
    real(kind=real64)            :: aux(10)                       ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)              ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi(e%n_gnodes)          ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: pphi(e%n_pnodes)              ! Functional shape functions values
    real(kind=real64)            :: sphi(e%n_snodes)              ! Functional shape functions values
    real(kind=real64)            :: x(2), T(2), N(2)              ! Normal vector at xi
    real(kind=real64)            :: rv(2), r, dr1, dr2, logr      ! Distance vector and 1/r^n
    real(kind=real64)            :: drdx(2)                       ! Distance vector derivatives with respect to x_k
    real(kind=real64)            :: jg                            ! Geometric jacobian
    real(kind=real64)            :: drdn                          ! Partial derivative of r respect to unit normal
    real(kind=real64)            :: drdt                          ! Partial derivative of r respect to unit tangent
    real(kind=real64)            :: jw                            ! Jacobians (except geometric jacobian) * weight
    real(kind=real64)            :: pphijw(e%n_pnodes)            ! Functional shape functions values * jw
    real(kind=real64)            :: sphijw(e%n_snodes)            ! Functional shape functions values * jw
    real(kind=real64)            :: dme_gphi(e%dme_n_gnodes)      ! phi^{dme-g} at integration point
    real(kind=real64)            :: dme_wqj(e%dme_n_gnodes,2)     ! w_{q,j}=dphi^{dme}_q/dx_j at integration point
    real(kind=real64)            :: dme_vq(e%dme_n_gnodes)        ! v_{q}=w_{q,j}·t_j at integration point
    real(kind=real64)            :: dme_d                         ! Dimensionless distance between DME and a point
    real(kind=real64)            :: dme_xi(e%dme_d)               ! Local coordinate of the DME
    real(kind=real64)            :: fs_q, fs_qt, fs_dpdm, fs_dqdm ! Fundamental solutions values
    ! Initialization
    h1=0.d0
    h2=0.d0
    g1=0.d0
    g2=0.d0
    ! Initialize dme_xi
    dme_xi=0.d0
    telles_parameters=fbem_telles11_calculate_parameters(barxip(1),barr)
    ! Numerical integration
    do kip=1,gl11_n(gln)
      ! gamma [-1,1]
      gamma=gl11_xi(kip,gln)
      w=gl11_w(kip,gln)
      ! gamma [-1,1] -> xip [-1,1]
      call fbem_telles_xi_and_jacobian(telles_parameters,gamma,xip,jt)
      ! xip [-1,1] -> xi [-1,1]
      xi=0.5d0*(1.d0-xip)*xi_s(1,1)+0.5d0*(1.d0+xip)*xi_s(1,2)
      js=0.5d0*(xi_s(1,2)-xi_s(1,1))
      ! xi [-1,1] -> x, T
      ! Geometrical shape functions and their derivatives
#     define etype e%gtype
#     define delta 0.d0
#     define phi gphi
#     define dphidxi dgphidxi
#     include <phi_and_dphidxi_1d.rc>
#     undef etype
#     undef delta
#     undef phi
#     undef dphidxi
      ! Position vector x and tangent vector T
      x=0.d0
      T=0.d0
      do kphi=1,e%n_gnodes
        x=x+gphi(kphi)*e%x(:,kphi)
        T=T+dgphidxi(kphi)*e%x(:,kphi)
      end do
      ! Jacobian: xi [-1,1] -> x
      jg=sqrt(dot_product(T,T))
      ! Unit tangent and unit normal
      t=T/jg
      n(1)=t(2)
      n(2)=-t(1)
      ! Distance vector and its modulus, inverse and natural logarithm
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      dr1=1.d0/r
      dr2=1.d0/r**2
      logr=log(r)
      ! dr/dx
      drdx=rv*dr1
      ! dr/dn
      drdn=dot_product(drdx,n)
      ! dr/dt
      drdt=dot_product(drdx,t)
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
      ! Functional shape functions (primary variables)
#     define etype e%ptype
#     define delta e%ptype_delta
#     define phi pphi
#     include <phi_1d.rc>
#     undef etype
#     undef delta
#     undef phi
      ! Functional shape functions (secondary variables)
#     define etype e%stype
#     define delta e%stype_delta
#     define phi sphi
#     include <phi_1d.rc>
#     undef etype
#     undef delta
#     undef phi
      ! jacobians * weight
      jw=jg*js*jt*w
      ! functional shape functions * jw
      pphijw=pphi*jw
      sphijw=sphi*jw
      ! Add the integration point
      fs_q=dr1*drdn
      fs_qt=dr1*drdt
      do im=1,2
        fs_dpdm=dr1*drdx(im)
        fs_dqdm=dr2*(-2.d0*drdx(im)*drdn+n(im))
        do iq=1,e%dme_n_gnodes
          h1(iq,:,im)=h1(iq,:,im)+fs_q*dme_vq(iq)*t(im)*pphijw(:)
          h1(iq,:,im)=h1(iq,:,im)+fs_dqdm*dme_gphi(iq)*pphijw(:)
          h1(iq,:,im)=h1(iq,:,im)-fs_qt*dme_vq(iq)*n(im)*pphijw(:)
          g1(iq,:,im)=g1(iq,:,im)+logr*dme_vq(iq)*t(im)*sphijw(:)
          g1(iq,:,im)=g1(iq,:,im)+fs_dpdm*dme_gphi(iq)*sphijw(:)
        end do
        h2(:,im)=h2(:,im)+fs_dqdm*pphijw(:)
        g2(:,im)=g2(:,im)+fs_dpdm*sphijw(:)
      end do
    end do
    ! Multiply by constants
    h1=-c_1_2pi*h1
    h2=-c_1_2pi*h2
    g1=-c_1_2pi*g1
    g2=-c_1_2pi*g2
    ! Reverse element orientation
    if (reverse) then
      h1=-h1
      h2=-h2
    end if
  end subroutine fbem_bem_stapot2d_vsbie_ext_st

  recursive subroutine fbem_bem_stapot2d_vsbie_ext_adp(e,reverse,xi_s,x_i,qsp,ks,ns,h1,h2,g1,g2)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                               !! Element
    logical                  :: reverse                         !! Reverse orientation
    real(kind=real64)        :: xi_s(1,2)                       !! Subdivision of the parent element
    real(kind=real64)        :: x_i(2)                          !! Collocation point position vector
    type(fbem_qs_parameters) :: qsp                             !! Quasi-singular integration parameters
    integer                  :: ks                              !! Current level of subdivisions
    integer                  :: ns                              !! Maximum level of subdivision
    real(kind=real64)        :: h1(e%dme_n_gnodes,e%n_pnodes,2) !! h1 matrix
    real(kind=real64)        :: h2(               e%n_pnodes,2) !! h2 matrix
    real(kind=real64)        :: g1(e%dme_n_gnodes,e%n_snodes,2) !! g1 matrix
    real(kind=real64)        :: g2(               e%n_snodes,2) !! g2 matrix
    ! Local
    integer           :: gln_near                            ! Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer           :: gln                                 ! Gauss-Legendre integ. points required to integrate the whole integrand
    logical           :: subdivide                           ! True if subdivision has to be performed
    real(kind=real64) :: barxi(1)                            ! Nearest element coordinate with respect to collocation point
    real(kind=real64) :: barxip(1)                           ! Nearest element subdivision local coordinate with respect to collocation point
    real(kind=real64) :: rmin                                ! Minimum distance between collocation point and barxi on the element
    real(kind=real64) :: barr                                ! Telles jacobian at barxi
    real(kind=real64) :: cl                                  ! Characteristic length
    real(kind=real64) :: d                                   ! Normalized distance between collocation point and element subdivision
    integer           :: method                              ! Method used in nearest point algorithm
    real(kind=real64) :: tmp_xi_s(1,2)                       ! Subdivision
    real(kind=real64) :: x_s(2,e%n_gnodes)                   ! Coordinates of the element subdivision
    real(kind=real64) :: h1_tmp(e%dme_n_gnodes,e%n_pnodes,2) ! h1 matrix (temporary)
    real(kind=real64) :: h2_tmp(               e%n_pnodes,2) ! h2 matrix (temporary)
    real(kind=real64) :: g1_tmp(e%dme_n_gnodes,e%n_snodes,2) ! g1 matrix (temporary)
    real(kind=real64) :: g2_tmp(               e%n_snodes,2) ! g2 matrix (temporary)
    ! Initialize
    if (ks.eq.1) then
      h1=0.d0
      h2=0.d0
      g1=0.d0
      g2=0.d0
      xi_s(1,1)=-1.d0
      xi_s(1,2)= 1.d0
      call fbem_nearest_element_point_bem(2,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
      barxip=barxi
    else
      call fbem_obtain_element_subdivision_coordinates(2,e%gtype,e%x,xi_s,x_s)
      cl=fbem_characteristic_length(2,e%gtype,x_s,1.d-6)
      call fbem_nearest_element_point_bem(2,e%gtype,x_s,cl,x_i,barxip,rmin,d,method)
    end if
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,4,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_stapot2d_vsbie_ext_adp',ns,'maximum number of subdivisions reached')
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
      call fbem_bem_stapot2d_vsbie_ext_adp(e,reverse,tmp_xi_s,x_i,qsp,ks+1,ns,h1,h2,g1,g2)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_stapot2d_vsbie_ext_adp(e,reverse,tmp_xi_s,x_i,qsp,ks+1,ns,h1,h2,g1,g2)
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_stapot2d_vsbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,gln,h1_tmp,h2_tmp,g1_tmp,g2_tmp)
      h1=h1+h1_tmp
      h2=h2+h2_tmp
      g1=g1+g1_tmp
      g2=g2+g2_tmp
    end if
  end subroutine fbem_bem_stapot2d_vsbie_ext_adp

  subroutine fbem_bem_stapot2d_vsbie_int(e,reverse,xi_i,h1,g1)
    implicit none
    ! I/O
    type(fbem_bem_element) :: e                               !! Integration element
    logical                :: reverse                         !! Reverse normal vector
    real(kind=real64)      :: xi_i(1)                         !! Local coordinate of the singular point.
    real(kind=real64)      :: h1(e%dme_n_gnodes,e%n_pnodes,2) !! h1 matrix
    real(kind=real64)      :: g1(e%dme_n_gnodes,e%n_snodes,2) !! g1 matrix
    ! Local
    integer           :: gln                         ! Gauss-Legendre number of integration points (<=32)
    integer           :: im, iq                      ! Counter variables
    integer           :: kphi                        ! Counter variable for shape functions loops
    integer           :: kip                         ! Counter of integration points
    real(kind=real64) :: gphi(e%n_gnodes)            ! Geometrical shape functions values at xi
    real(kind=real64) :: dgphidxi(e%n_gnodes)        ! Geometrical shape functions derivatives values at xi
    real(kind=real64) :: gphi_i(e%n_gnodes)          ! Geometrical shape functions values at xi_i
    real(kind=real64) :: dgphidxi_i(e%n_gnodes)      ! Geometrical shape functions derivatives values at xi_i
    real(kind=real64) :: pphi(e%n_pnodes)            ! Functional shape functions values at xi
    real(kind=real64) :: pphi_i(e%n_pnodes)          ! Functional shape functions values at xi_i
    real(kind=real64) :: sphi(e%n_snodes)            ! Functional shape functions values at xi
    integer           :: nsub                        ! Number of subdivision of the element
    integer           :: ksub                        ! Counter of subdivision
    real(kind=real64) :: w                           ! Weights of each integration point
    real(kind=real64) :: xip                         ! Coordinate  xip of subdivided element [0,1]
    real(kind=real64) :: js                          ! Jacobian of the xi [xisub(1,:),xisub[2,:]] -> xip [0,1] transformation
    real(kind=real64) :: xip_i(2)                    ! Singular point in xip space
    real(kind=real64) :: xi                          ! Coordinate xi
    real(kind=real64) :: xisub(2,2)                  ! Coordinates of element subdivisions
    real(kind=real64) :: aux(10)                     ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64) :: x(2), T(2), N(2)            ! Position, tangent and normal vectors
    real(kind=real64) :: x_i(2), t_i(2), n_i(2)      ! Position, tangent and normal vectors at collocation point
    real(kind=real64) :: rv(2), r, dr1, dr2, logr    ! Distance vector
    real(kind=real64) :: drdx(2)                     ! Distance vector derivatives with respect to x_k
    real(kind=real64) :: drdti                       ! (dr/dGamma)^i=+-1, -1 antes de i, 1 despues de i
    real(kind=real64) :: jg                          ! Geometric jacobian
    real(kind=real64) :: jgi                         ! Geometric jacobian at xi_i
    real(kind=real64) :: drdn                        ! Partial derivative of r respect to unit normal
    real(kind=real64) :: drdt                        ! Partial derivative of r respect to unit tangent
    real(kind=real64) :: jw                          ! Jacobians * weight
    real(kind=real64) :: dme_gphi(e%dme_n_gnodes)    ! phi^{dme-g} at integration point
    real(kind=real64) :: dme_gphi_i(e%dme_n_gnodes)  ! phi^{dme-g} at collocation point
    real(kind=real64) :: dme_wqj(e%dme_n_gnodes,2)   ! w_{q,j}=dphi^{dme}_q/dx_j at integration point
    real(kind=real64) :: dme_wqj_i(e%dme_n_gnodes,2) ! w_{q,j}=dphi^{dme}_q/dx_j at collocation point
    real(kind=real64) :: dme_vq(e%dme_n_gnodes)      ! v_{q}=w_{q,j}·t_j at integration point
    real(kind=real64) :: dme_vq_i(e%dme_n_gnodes)    ! v_{q}=w_{q,j}·t_j at collocation point
    real(kind=real64) :: dme_d                       ! Dimensionless distance between DME and a point
    real(kind=real64) :: dme_xi(e%dme_d)             ! Local coordinate of the DME
    real(kind=real64) :: h(e%n_pnodes)               ! h matrix
    real(kind=real64) :: g(e%n_snodes)               ! g matrix
    ! Initialization
    gln=32
    h1=0.d0
    g1=0.d0
    ! Initialize dme_xi
    dme_xi=0.d0
    ! Setup the subdivisions for transformation xip -> xi: [0,1] -> [xisub_1,xisub_2]
    if (fbem_check_xi_vertex(xi_i(1))) then
      nsub=1
      xisub(1,1)=-1.0d0
      xisub(2,1)=1.0d0
      if (xi_i(1).lt.0.0d0) xip_i=0.0d0
      if (xi_i(1).gt.0.0d0) xip_i=1.0d0
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
    ! Unit tangent and normal vectors
    t_i=T_i/jgi
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
    ! Numerical integration
    do ksub=1,nsub
      do kip=1,gl01_n(gln)
        ! xip [0,1]
        xip=gl01_xi(kip,gln)
        w=gl01_w(kip,gln)
        ! xip [0,1] -> xi [-1,1]
        js=xisub(2,ksub)-xisub(1,ksub)
        xi=js*xip+xisub(1,ksub)
        ! dr/dGamma
        if (xi.lt.xi_i(1)) then
          drdti=-1.d0
        else
          drdti= 1.d0
        end if
        ! xi [-1,1] -> x, T
        ! Geometrical shape functions and their derivatives
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
        ! Unit tangent and normal vectors
        t=T/jg
        n(1)=t(2)
        n(2)=-t(1)
        ! Distance vector
        rv=x-x_i
        r=sqrt(dot_product(rv,rv))
        dr1=1.d0/r
        dr2=dr1**2
        logr=log(r)
        drdx=rv*dr1
        drdn=dot_product(drdx,n)
        drdt=dot_product(drdx,t)
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
        ! Functional shape functions (primary variables)
#       define etype e%ptype
#       define delta e%ptype_delta
#       define phi pphi
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Functional shape functions (secondary variables)
#       define etype e%stype
#       define delta e%stype_delta
#       define phi sphi
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Jacobian * weight
        jw=jg*js*w
        ! Add the integration point
        do im=1,2
          do iq=1,e%dme_n_gnodes
            ! G^s
            g1(iq,:,im)=g1(iq,:,im)+logr*(dme_vq(iq)*t(im)-dme_vq_i(iq)*t_i(im))*sphi(:)*jw
            ! H^s
            h1(iq,:,im)=h1(iq,:,im)+dr1*drdn*(dme_vq(iq)*t(im)-dme_vq_i(iq)*t_i(im))*pphi(:)*jw
            ! G^r
            g1(iq,:,im)=g1(iq,:,im)+dr1*drdx(im)*(dme_gphi(iq)-dme_gphi_i(iq))*sphi(:)*jw
            ! H^r
            h1(iq,:,im)=h1(iq,:,im)-2.d0*drdx(im)*dr2*drdn*(dme_gphi(iq)-dme_gphi_i(iq))*pphi(:)*jw
            h1(iq,:,im)=h1(iq,:,im)+n(im)*dr2*(dme_gphi(iq)-dme_gphi_i(iq)-dot_product(dme_wqj_i(iq,:),rv))*pphi(:)*jw
            h1(iq,:,im)=h1(iq,:,im)+dr1*(dot_product(dme_wqj_i(iq,:),drdx)*n(im)*pphi(:)-drdti*dme_vq_i(iq)*n_i(im)*pphi_i(:))*jw
            h1(iq,:,im)=h1(iq,:,im)+dme_vq_i(iq)*n_i(im)*pphi_i(:)*dr1*(drdti-drdt)*jw
            ! H^n
            h1(iq,:,im)=h1(iq,:,im)-dr1*drdt*(dme_vq(iq)*n(im)*pphi(:)-dme_vq_i(iq)*n_i(im)*pphi_i(:))*jw
          end do
        end do
      end do
    end do
    ! Multiply by constants
    h1=-c_1_2pi*h1
    g1=-c_1_2pi*g1
    ! Reverse element orientation
    if (reverse) h1=-h1
    ! Assemble matrices from SBIE
    call fbem_bem_stapot2d_sbie_int(e,reverse,xi_i,h,g)
    do im=1,2
      do iq=1,e%dme_n_gnodes
        h1(iq,:,im)=h1(iq,:,im)+dme_vq_i(iq)*t_i(im)*h(:)
        g1(iq,:,im)=g1(iq,:,im)+dme_vq_i(iq)*t_i(im)*g(:)
      end do
    end do
  end subroutine fbem_bem_stapot2d_vsbie_int

  subroutine fbem_bem_stapot2d_vsbie_auto(e,reverse,x_i,qsp,ns,h1,h2,g1,g2)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                               !! Integration element
    logical                  :: reverse                         !! Reverse orientation
    real(kind=real64)        :: x_i(2)                          !! Collocation point
    type(fbem_qs_parameters) :: qsp                             !! Quasi-singular integration parameters
    integer                  :: ns                              !! Maximum level of subdivisions
    real(kind=real64)        :: h1(e%dme_n_gnodes,e%n_pnodes,2) !! h1 matrix
    real(kind=real64)        :: h2(               e%n_pnodes,2) !! h2 matrix
    real(kind=real64)        :: g1(e%dme_n_gnodes,e%n_snodes,2) !! g1 matrix
    real(kind=real64)        :: g2(               e%n_snodes,2) !! g2 matrix
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
        h2=0.d0
        g2=0.d0
        call fbem_bem_stapot2d_vsbie_int(e,reverse,barxi,h1,g1)
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
          call fbem_bem_stapot2d_vsbie_ext_pre(ps,e,reverse,x_i,h1,h2,g1,g2)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_stapot2d_vsbie_ext_adp(e,reverse,xi_s,x_i,qsp,1,ns,h1,h2,g1,g2)
        end if
    end select
  end subroutine fbem_bem_stapot2d_vsbie_auto

  ! ================================================================================================================================

  ! ================================================================================================================================
  ! VARIATION HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (VHBIE)

  subroutine fbem_bem_stapot2d_vhbie_ext_pre(ps,e,reverse,x_i,n_i,m1,m2,m3,l1,l2,l3)
    implicit none
    ! I/O
    integer                :: ps                              !! Selected precalculated dataset
    type(fbem_bem_element) :: e                               !! Element
    logical                :: reverse                         !! Reverse normal vector
    real(kind=real64)      :: x_i(2)                          !! Position vector of the collocation point
    real(kind=real64)      :: n_i(2)                          !! Unit normal vector at the collocation point
    real(kind=real64)      :: m1(e%dme_n_gnodes,e%n_pnodes,2) !! m1 matrix
    real(kind=real64)      :: m2(               e%n_pnodes,2) !! m2 matrix
    real(kind=real64)      :: m3(               e%n_pnodes,2) !! m3 matrix
    real(kind=real64)      :: l1(e%dme_n_gnodes,e%n_snodes,2) !! l1 matrix
    real(kind=real64)      :: l2(               e%n_snodes,2) !! l2 matrix
    real(kind=real64)      :: l3(               e%n_snodes,2) !! l3 matrix
    ! Local
    integer           :: im, iq                        ! Counter variables
    integer           :: kip                           ! Counter variable for integration points loop
    real(kind=real64) :: x(2)                          ! Position vector at integration point
    real(kind=real64) :: t(2)                          ! Unit tangent vector at integration point
    real(kind=real64) :: n(2)                          ! Unit normal vector at integration point
    real(kind=real64) :: pphijw(e%n_pnodes)            ! phi^p * jacobian * weight at integration point
    real(kind=real64) :: sphijw(e%n_snodes)            ! phi^s at integration point
    real(kind=real64) :: rv(2), r, dr1, dr2, dr3       ! Distance vector module and 1/r^n
    real(kind=real64) :: drdx(2)                       ! Distance vector derivatives with respect to x_k
    real(kind=real64) :: drdn, drdni, n_dot_ni         ! dr/dn, dr/dn^i, n·n_i
    real(kind=real64) :: drdt, ni_dot_t                ! dr/dt and n^i·t
    real(kind=real64) :: dme_gphi(e%dme_n_gnodes)      ! phi^{dme-g} at integration point
    real(kind=real64) :: dme_wqj(e%dme_n_gnodes,2)     ! dphi{dme-g}/dx_j at integration point
    real(kind=real64) :: dme_vq(e%dme_n_gnodes)        ! vq=w_{j,q}·t_j at integration point
    real(kind=real64) :: fs_s, fs_d                    ! Fundamental solutions values
    real(kind=real64) :: fs_dr, fs_dni                 ! Fundamental solutions values
    real(kind=real64) :: fs_sr, fs_sni, fs_sn          ! Fundamental solutions values
    ! Initialization
    m1=0.d0
    m2=0.d0
    m3=0.d0
    l1=0.d0
    l2=0.d0
    l3=0.d0
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
      dr2=dr1**2
      dr3=dr2*dr1
      drdx=rv*dr1
      drdn=dot_product(drdx,n)
      drdni=-dot_product(drdx,n_i)
      drdt=dot_product(drdx,t)
      n_dot_ni=dot_product(n,n_i)
      ni_dot_t=dot_product(n_i,t)
      ! Add the integration point
      fs_d=-dr1*drdni
      fs_s=dr2*(2.d0*drdn*drdni+n_dot_ni)
      fs_sn=-dr2*(ni_dot_t+2.d0*drdt*drdni)
      do im=1,2
        fs_dr=dr2*(n_i(im)+2.d0*drdx(im)*drdni)
        fs_dni=dr1*drdx(im)
        fs_sr=-2.d0*dr3*((4.d0*drdn*drdni+n_dot_ni)*drdx(im)+n_i(im)*drdn-n(im)*drdni)
        fs_sni=dr2*(n(im)-2.d0*drdx(im)*drdn)
        do iq=1,e%dme_n_gnodes
          m1(iq,:,im)=m1(iq,:,im)+fs_s*dme_vq(iq)*t(im)*pphijw(:)
          m1(iq,:,im)=m1(iq,:,im)+fs_sr*dme_gphi(iq)*pphijw(:)
          m1(iq,:,im)=m1(iq,:,im)+fs_sn*dme_vq(iq)*n(im)*pphijw(:)
          l1(iq,:,im)=l1(iq,:,im)+fs_d*dme_vq(iq)*t(im)*sphijw(:)
          l1(iq,:,im)=l1(iq,:,im)+fs_dr*dme_gphi(iq)*sphijw(:)
        end do
        m2(:,im)=m2(:,im)+fs_sr*pphijw(:)
        l2(:,im)=l2(:,im)+fs_dr*sphijw(:)
        m3(:,im)=m3(:,im)+fs_sni*pphijw(:)
        l3(:,im)=l3(:,im)+fs_dni*sphijw(:)
      end do
    end do
    ! Multiply by constants
    m1=c_1_2pi*m1
    m2=c_1_2pi*m2
    m3=c_1_2pi*m3
    l1=c_1_2pi*l1
    l2=c_1_2pi*l2
    l3=c_1_2pi*l3
    ! Reverse element orientation
    if (reverse) then
      m1=-m1
      m2=-m2
      m3=-m3
    end if
  end subroutine fbem_bem_stapot2d_vhbie_ext_pre

  subroutine fbem_bem_stapot2d_vhbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,gln,m1,m2,m3,l1,l2,l3)
    implicit none
    ! I/O
    type(fbem_bem_element) :: e                               !! Integration element
    logical                :: reverse                         !! Reverse normal vector
    real(kind=real64)      :: xi_s(1,2)                       !! Coordinates of the subdivision of the parent element (xi space [-1,1])
    real(kind=real64)      :: x_i(2)                          !! Position vector of he collocation point
    real(kind=real64)      :: n_i(2)                          !! Unit normal vector at the collocation point
    real(kind=real64)      :: barxip(1)                       !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)      :: barr                            !! Telles jacobian at barxip
    integer                :: gln                             !! Gauss-Legendre number of integration points (<=32)
    real(kind=real64)      :: m1(e%dme_n_gnodes,e%n_pnodes,2) !! m1 matrix
    real(kind=real64)      :: m2(               e%n_pnodes,2) !! m2 matrix
    real(kind=real64)      :: m3(               e%n_pnodes,2) !! m3 matrix
    real(kind=real64)      :: l1(e%dme_n_gnodes,e%n_snodes,2) !! l1 matrix
    real(kind=real64)      :: l2(               e%n_snodes,2) !! l2 matrix
    real(kind=real64)      :: l3(               e%n_snodes,2) !! l3 matrix
    ! Local
    integer                      :: im, iq                        ! Counter variables
    integer                      :: kphi                          ! Counter variable for shape functions loops
    integer                      :: kip                           ! Counter variable of integration points
    real(kind=real64)            :: gamma                         ! Coordinate gamma (Telles transformation space [-1,1])
    real(kind=real64)            :: w                             ! Weights of an integration point
    type(fbem_telles_parameters) :: telles_parameters             ! Telles parameters
    real(kind=real64)            :: jt                            ! Telles jacobian: xip->gamma
    real(kind=real64)            :: xip                           ! Coordinate xip (subdivision space [-1,1])
    real(kind=real64)            :: js                            ! Subdivision jacobian: xi->xip
    real(kind=real64)            :: xi                            ! Coordinate xi [xi_s(1,1),xi_s(1,2)]
    real(kind=real64)            :: aux(10)                       ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)              ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi(e%n_gnodes)          ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: pphi(e%n_pnodes)              ! Functional shape functions values
    real(kind=real64)            :: sphi(e%n_snodes)              ! Functional shape functions values
    real(kind=real64)            :: x(2), T(2), N(2)              ! Position, tangent and normal vectors
    real(kind=real64)            :: rv(2), r, dr1, dr2, dr3       ! Distance vector and 1/r^n
    real(kind=real64)            :: drdx(2)                       ! Distance vector derivatives with respect to x_k
    real(kind=real64)            :: jg                            ! Geometric jacobian
    real(kind=real64)            :: drdn, drdni, n_dot_ni         ! dr/dn, dr/dn^i, n·n_i
    real(kind=real64)            :: drdt, ni_dot_t                ! dr/dt and n^i·t
    real(kind=real64)            :: jw                            ! Jacobians (except geometric jacobian) * weight
    real(kind=real64)            :: pphijw(e%n_pnodes)            ! Functional shape functions values * jw
    real(kind=real64)            :: sphijw(e%n_snodes)            ! Functional shape functions values * jw
    real(kind=real64)            :: dme_gphi(e%dme_n_gnodes)      ! phi^{dme-g} at integration point
    real(kind=real64)            :: dme_wqj(e%dme_n_gnodes,2)     ! w_{q,j}=dphi^{dme}_q/dx_j at integration point
    real(kind=real64)            :: dme_vq(e%dme_n_gnodes)        ! v_{q}=w_{q,j}·t_j at integration point
    real(kind=real64)            :: dme_d                         ! Dimensionless distance between DME and a point
    real(kind=real64)            :: dme_xi(e%dme_d)               ! Local coordinate of the DME
    real(kind=real64)            :: fs_s, fs_d                    ! Fundamental solutions values
    real(kind=real64)            :: fs_dr, fs_dni                 ! Fundamental solutions values
    real(kind=real64)            :: fs_sr, fs_sni, fs_sn          ! Fundamental solutions values
    ! Initialization
    m1=0.d0
    m2=0.d0
    m3=0.d0
    l1=0.d0
    l2=0.d0
    l3=0.d0
    ! Initialize dme_xi
    dme_xi=0.d0
    telles_parameters=fbem_telles11_calculate_parameters(barxip(1),barr)
    ! Numerical integration
    do kip=1,gl11_n(gln)
      ! gamma [-1,1]
      gamma=gl11_xi(kip,gln)
      w=gl11_w(kip,gln)
      ! gamma [-1,1] -> xip [-1,1]
      call fbem_telles_xi_and_jacobian(telles_parameters,gamma,xip,jt)
      ! xip [-1,1] -> xi [-1,1]
      xi=0.5d0*(1.d0-xip)*xi_s(1,1)+0.5d0*(1.d0+xip)*xi_s(1,2)
      js=0.5d0*(xi_s(1,2)-xi_s(1,1))
      ! xi [-1,1] -> x, T
      ! Geometrical shape functions and their derivatives
#     define etype e%gtype
#     define delta 0.d0
#     define phi gphi
#     define dphidxi dgphidxi
#     include <phi_and_dphidxi_1d.rc>
#     undef etype
#     undef delta
#     undef phi
#     undef dphidxi
      ! Position vector x and tangent vector T
      x=0.d0
      T=0.d0
      do kphi=1,e%n_gnodes
        x=x+gphi(kphi)*e%x(:,kphi)
        T=T+dgphidxi(kphi)*e%x(:,kphi)
      end do
      ! Jacobian: xi [-1,1] -> x
      jg=sqrt(dot_product(T,T))
      ! Unit tangent and unit normal
      t=T/jg
      n(1)=t(2)
      n(2)=-t(1)
      ! Distance vector and relatives
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      dr1=1.d0/r
      dr2=dr1**2
      dr3=dr2*dr1
      drdx=rv*dr1
      drdn=dot_product(drdx,n)
      drdni=-dot_product(drdx,n_i)
      drdt=dot_product(drdx,t)
      n_dot_ni=dot_product(n,n_i)
      ni_dot_t=dot_product(n_i,t)
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
      ! Functional shape functions (primary variables)
#     define etype e%ptype
#     define delta e%ptype_delta
#     define phi pphi
#     include <phi_1d.rc>
#     undef etype
#     undef delta
#     undef phi
      ! Functional shape functions (secondary variables)
#     define etype e%stype
#     define delta e%stype_delta
#     define phi sphi
#     include <phi_1d.rc>
#     undef etype
#     undef delta
#     undef phi
      ! jacobians * weight
      jw=jg*js*jt*w
      ! functional shape functions * jw
      pphijw=pphi*jw
      sphijw=sphi*jw
      ! Add the integration point
      fs_d=-dr1*drdni
      fs_s=dr2*(2.d0*drdn*drdni+n_dot_ni)
      fs_sn=-dr2*(ni_dot_t+2.d0*drdt*drdni)
      do im=1,2
        fs_dr=dr2*(n_i(im)+2.d0*drdx(im)*drdni)
        fs_dni=dr1*drdx(im)
        fs_sr=-2.d0*dr3*((4.d0*drdn*drdni+n_dot_ni)*drdx(im)+n_i(im)*drdn-n(im)*drdni)
        fs_sni=dr2*(n(im)-2.d0*drdx(im)*drdn)
        do iq=1,e%dme_n_gnodes
          m1(iq,:,im)=m1(iq,:,im)+fs_s*dme_vq(iq)*t(im)*pphijw(:)
          m1(iq,:,im)=m1(iq,:,im)+fs_sr*dme_gphi(iq)*pphijw(:)
          m1(iq,:,im)=m1(iq,:,im)+fs_sn*dme_vq(iq)*n(im)*pphijw(:)
          l1(iq,:,im)=l1(iq,:,im)+fs_d*dme_vq(iq)*t(im)*sphijw(:)
          l1(iq,:,im)=l1(iq,:,im)+fs_dr*dme_gphi(iq)*sphijw(:)
        end do
        m2(:,im)=m2(:,im)+fs_sr*pphijw(:)
        l2(:,im)=l2(:,im)+fs_dr*sphijw(:)
        m3(:,im)=m3(:,im)+fs_sni*pphijw(:)
        l3(:,im)=l3(:,im)+fs_dni*sphijw(:)
      end do
    end do
    ! Multiply by constants
    m1=c_1_2pi*m1
    m2=c_1_2pi*m2
    m3=c_1_2pi*m3
    l1=c_1_2pi*l1
    l2=c_1_2pi*l2
    l3=c_1_2pi*l3
    ! Reverse element orientation
    if (reverse) then
      m1=-m1
      m2=-m2
      m3=-m3
    end if
  end subroutine fbem_bem_stapot2d_vhbie_ext_st

  recursive subroutine fbem_bem_stapot2d_vhbie_ext_adp(e,reverse,xi_s,x_i,n_i,qsp,ks,ns,m1,m2,m3,l1,l2,l3)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                               !! Element
    logical                  :: reverse                         !! Reverse orientation
    real(kind=real64)        :: xi_s(1,2)                       !! Subdivision of the parent element
    real(kind=real64)        :: x_i(2)                          !! Position vector of he collocation point
    real(kind=real64)        :: n_i(2)                          !! Unit normal vector at the collocation point
    type(fbem_qs_parameters) :: qsp                             !! Quasi-singular integration parameters
    integer                  :: ks                              !! Current level of subdivisions
    integer                  :: ns                              !! Maximum level of subdivision
    real(kind=real64)        :: m1(e%dme_n_gnodes,e%n_pnodes,2) !! m1 matrix
    real(kind=real64)        :: m2(               e%n_pnodes,2) !! m2 matrix
    real(kind=real64)        :: m3(               e%n_pnodes,2) !! m3 matrix
    real(kind=real64)        :: l1(e%dme_n_gnodes,e%n_snodes,2) !! l1 matrix
    real(kind=real64)        :: l2(               e%n_snodes,2) !! l2 matrix
    real(kind=real64)        :: l3(               e%n_snodes,2) !! l3 matrix
    ! Local
    integer           :: gln_near                            ! Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer           :: gln                                 ! Gauss-Legendre integ. points required to integrate the whole integrand
    logical           :: subdivide                           ! True if subdivision has to be performed
    real(kind=real64) :: barxi(1)                            ! Nearest element coordinate with respect to collocation point
    real(kind=real64) :: barxip(1)                           ! Nearest element subdivision local coordinate with respect to collocation point
    real(kind=real64) :: rmin                                ! Minimum distance between collocation point and barxi on the element
    real(kind=real64) :: barr                                ! Telles jacobian at barxi
    real(kind=real64) :: cl                                  ! Characteristic length
    real(kind=real64) :: d                                   ! Normalized distance between collocation point and element subdivision
    integer           :: method                              ! Method used in nearest point algorithm
    real(kind=real64) :: tmp_xi_s(1,2)                       ! Subdivision
    real(kind=real64) :: x_s(2,e%n_gnodes)                   ! Coordinates of the element subdivision
    real(kind=real64) :: m1_tmp(e%dme_n_gnodes,e%n_pnodes,2) ! m1 matrix (temporary)
    real(kind=real64) :: m2_tmp(               e%n_pnodes,2) ! m2 matrix (temporary)
    real(kind=real64) :: m3_tmp(               e%n_pnodes,2) ! m3 matrix (temporary)
    real(kind=real64) :: l1_tmp(e%dme_n_gnodes,e%n_snodes,2) ! l1 matrix (temporary)
    real(kind=real64) :: l2_tmp(               e%n_snodes,2) ! l2 matrix (temporary)
    real(kind=real64) :: l3_tmp(               e%n_snodes,2) ! l3 matrix (temporary)
    ! Initialize
    if (ks.eq.1) then
      m1=0.d0
      m2=0.d0
      m3=0.d0
      l1=0.d0
      l2=0.d0
      l3=0.d0
      xi_s(1,1)=-1.d0
      xi_s(1,2)= 1.d0
      call fbem_nearest_element_point_bem(2,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
      barxip=barxi
    else
      call fbem_obtain_element_subdivision_coordinates(2,e%gtype,e%x,xi_s,x_s)
      cl=fbem_characteristic_length(2,e%gtype,x_s,1.d-6)
      call fbem_nearest_element_point_bem(2,e%gtype,x_s,cl,x_i,barxip,rmin,d,method)
    end if
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,6,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_stapot2d_vhbie_ext_adp',ns,'maximum number of subdivisions reached')
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
      call fbem_bem_stapot2d_vhbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,qsp,ks+1,ns,m1,m2,m3,l1,l2,l3)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_stapot2d_vhbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,qsp,ks+1,ns,m1,m2,m3,l1,l2,l3)
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_stapot2d_vhbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,gln,m1_tmp,m2_tmp,m3_tmp,l1_tmp,l2_tmp,l3_tmp)
      m1=m1+m1_tmp
      m2=m2+m2_tmp
      m3=m3+m3_tmp
      l1=l1+l1_tmp
      l2=l2+l2_tmp
      l3=l3+l3_tmp
    end if
  end subroutine fbem_bem_stapot2d_vhbie_ext_adp

  subroutine fbem_bem_stapot2d_vhbie_auto(e,reverse,x_i,n_i,qsp,ns,m1,m2,m3,l1,l2,l3)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                               !! Integration element
    logical                  :: reverse                         !! Reverse orientation
    real(kind=real64)        :: x_i(2)                          !! Position vector of the collocation point
    real(kind=real64)        :: n_i(2)                          !! Unit normal at the collocation point
    type(fbem_qs_parameters) :: qsp                             !! Quasi-singular integration parameters
    integer                  :: ns                              !! Maximum level of subdivisions
    real(kind=real64)        :: m1(e%dme_n_gnodes,e%n_pnodes,2) !! m1 matrix
    real(kind=real64)        :: m2(               e%n_pnodes,2) !! m2 matrix
    real(kind=real64)        :: m3(               e%n_pnodes,2) !! m3 matrix
    real(kind=real64)        :: l1(e%dme_n_gnodes,e%n_snodes,2) !! l1 matrix
    real(kind=real64)        :: l2(               e%n_snodes,2) !! l2 matrix
    real(kind=real64)        :: l3(               e%n_snodes,2) !! l3 matrix
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
        m2=0.d0
        m3=0.d0
        l2=0.d0
        l3=0.d0
        call fbem_bem_stapot2d_vhbie_int(e,reverse,barxi,m1,l1)
      case (0)
        ! Estimate the required integration rule
        gln_near=fbem_qs_n_estimation_standard(e%n,e%gtype,6,qsp,d,barxi)
        gln=max(e%gln_far,gln_near)
        ! Integrate using a conservative precalculated dataset
        if ((gln.le.e%ps_gln_max).and.(gln_near.gt.0)) then
          do i=1,e%n_ps
            if (e%ps_gln(i).ge.gln) then
              ps=i
              exit
            end if
          end do
          call fbem_bem_stapot2d_vhbie_ext_pre(ps,e,reverse,x_i,n_i,m1,m2,m3,l1,l2,l3)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_stapot2d_vhbie_ext_adp(e,reverse,xi_s,x_i,n_i,qsp,1,ns,m1,m2,m3,l1,l2,l3)
        end if
    end select
  end subroutine fbem_bem_stapot2d_vhbie_auto

  subroutine fbem_bem_stapot2d_vhbie_int(e,reverse,xi_i,m1,l1)
    implicit none
    ! I/O
    type(fbem_bem_element) :: e                               !! Integration element
    logical                :: reverse                         !! Reverse normal vector
    real(kind=real64)      :: xi_i(1)                         !! Local coordinate of the singular point.
    real(kind=real64)      :: m1(e%dme_n_gnodes,e%n_pnodes,2) !! m1 matrix
    real(kind=real64)      :: l1(e%dme_n_gnodes,e%n_snodes,2) !! l1 matrix
    ! Local
    integer           :: gln                             ! Gauss-Legendre number of integration points (<=32)
    integer           :: im, iq                          ! Counter variables
    integer           :: kphi                            ! Counter variable for shape functions loops
    integer           :: kip                             ! Counter of integration points
    real(kind=real64) :: gphi(e%n_gnodes)                ! Geometrical shape functions values at xi
    real(kind=real64) :: dgphidxi(e%n_gnodes)            ! Geometrical shape functions derivatives values at xi
    real(kind=real64) :: gphi_i(e%n_gnodes)              ! Geometrical shape functions values at xi_i
    real(kind=real64) :: dgphidxi_i(e%n_gnodes)          ! Geometrical shape functions derivatives values at xi_i
    real(kind=real64) :: d2gphidxi2_i(e%n_gnodes)        ! Geometrical shape functions second derivatives values at xi
    real(kind=real64) :: pphi(e%n_pnodes)                ! Functional shape functions values at xi
    real(kind=real64) :: sphi(e%n_snodes)                ! Functional shape functions values at xi
    real(kind=real64) :: pphi_i(e%n_pnodes)              ! Functional shape functions values at xi_i
    real(kind=real64) :: sphi_i(e%n_snodes)              ! Functional shape functions values at xi_i
    real(kind=real64) :: dpphidx_i(e%n_pnodes,2)         ! d pphi / dx at xi_i
    integer           :: nsub                            ! Number of subdivision of the element
    integer           :: ksub                            ! Counter of subdivision
    real(kind=real64) :: w                               ! Weights of each integration point
    real(kind=real64) :: xip                             ! Coordinate  xip of subdivided element [0,1]
    real(kind=real64) :: js                              ! Jacobian of the xi [xisub(1,:),xisub[2,:]] -> xip [0,1] transformation
    real(kind=real64) :: xip_i(2)                        ! Singular point in xip space
    real(kind=real64) :: xi                              ! Coordinate xi
    real(kind=real64) :: xisub(2,2)                      ! Coordinates of element subdivisions
    real(kind=real64) :: aux(10)                         ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64) :: x(2), T(2), N(2)                ! Position, tangent and normal vectors
    real(kind=real64) :: x_i(2), t_i(2), n_i(2)          ! Position, tangent and normal vectors at xi_i
    real(kind=real64) :: dtdr_i(2), dTdxi(2)             ! Derivative of the unit tangent vector
    real(kind=real64) :: rv(2), r, dr1, dr2, dr3         ! Distance vector
    real(kind=real64) :: drdx(2)                         ! Distance vector derivatives with respect to x_k
    real(kind=real64) :: jg                              ! Geometric jacobian
    real(kind=real64) :: jgi                             ! Geometric jacobian at xi_i
    real(kind=real64) :: drdn                            ! Partial derivative of r respect to unit normal
    real(kind=real64) :: drdt                            ! Partial derivative of r respect to unit tangent
    real(kind=real64) :: drdni                           ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64) :: n_dot_ni                        ! n·n_i
    real(kind=real64) :: drdti                           ! (dr/dGamma)^i=+-1, -1 antes de i, 1 despues de i
    real(kind=real64) :: drdx_dot_t, drdx_dot_ti         ! Auxiliary dot products
    real(kind=real64) :: n_dot_ti, t_dot_ni              ! Auxiliary dot products
    real(kind=real64) :: jw                              ! Jacobians * weight
    real(kind=real64) :: dme_gphi(e%dme_n_gnodes)        ! phi^{dme-g} at integration point
    real(kind=real64) :: dme_gphi_i(e%dme_n_gnodes)      ! phi^{dme-g} at collocation point
    real(kind=real64) :: dme_wqj(e%dme_n_gnodes,2)       ! w_{q,j}=dphi^{dme}_q/dx_j at integration point
    real(kind=real64) :: dme_wqj_i(e%dme_n_gnodes,2)     ! w_{q,j}=dphi^{dme}_q/dx_j at collocation point
    real(kind=real64) :: dme_vq(e%dme_n_gnodes)          ! v_{q}=w_{q,j}·t_j at integration point
    real(kind=real64) :: dme_vq_i(e%dme_n_gnodes)        ! v_{q}=w_{q,j}·t_j at collocation point
    real(kind=real64) :: dme_d                           ! Dimensionless distance between DME and a point
    real(kind=real64) :: dme_xi(e%dme_d)                 ! Local coordinate of the DME
    real(kind=real64) :: Idr1dr, Idr2sgndr               ! Int 1/r dr, Int 1/r^2 sgn(drdt) dr
    real(kind=real64) :: ra, rb                          ! Distance vector from collocation point to element vertices
    ! Initialization
    gln=32
    m1=0.d0
    l1=0.d0
    ! Initialize dme_xi
    dme_xi=0.d0
    ! Setup the subdivisions for transformation xip -> xi: [0,1] -> [xisub_1,xisub_2]
    if (fbem_check_xi_vertex(xi_i(1))) then
      call fbem_error_message(output_unit,0,'fbem_bem_stapot2d_vhbie_int',0,'the VHBIE does not allow vertex collocation.')
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
    ! Calculate geometrical vectors
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
    ! Unit tangent and normal vectors
    t_i=T_i/jgi
    n_i(1)=t_i(2)
    n_i(2)=-t_i(1)
    ! Derivative of the unit tangent vector
#   define etype e%gtype
#   define delta 0.d0
#   define xi xi_i(1)
#   define d2phidxi2 d2gphidxi2_i
#   include <d2phidxi2_1d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef d2phidxi2
    dTdxi=0.d0
    do kphi=1,e%n_gnodes
      dTdxi=dTdxi+d2gphidxi2_i(kphi)*e%x(:,kphi)
    end do
    dtdr_i=(dTdxi-t_i*dot_product(dTdxi,t_i))/jgi**2
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
    dpphidx_i=fbem_dphidx(2,e%ptype,e%ptype_delta,e%x,xi_i)
    ! Functional shape functions (secondary variable) at xi_i
#   define etype e%stype
#   define delta e%stype_delta
#   define xi xi_i(1)
#   define phi sphi_i
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
    ! Numerical integration
    do ksub=1,nsub
      do kip=1,gl01_n(gln)
        ! xip [0,1]
        xip=gl01_xi(kip,gln)
        w=gl01_w(kip,gln)
        ! xip [0,1] -> xi [-1,1]
        js=xisub(2,ksub)-xisub(1,ksub)
        xi=js*xip+xisub(1,ksub)
        ! dr/dGamma
        if (xi.lt.xi_i(1)) then
          drdti=-1.d0
        else
          drdti= 1.d0
        end if
        ! xi [-1,1] -> x, T
        ! Geometrical shape functions and their derivatives
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
        ! Unit tangent and normal vectors
        t=T/jg
        n(1)=t(2)
        n(2)=-t(1)
        ! Distance vector
        rv=x-x_i
        r=sqrt(dot_product(rv,rv))
        dr1=1.d0/r
        dr2=1.d0/r**2
        dr3=1.d0/r**3
        drdx=rv*dr1
        drdn=dot_product(drdx,n)
        drdni=-dot_product(drdx,n_i)
        drdt=dot_product(drdx,t)
        n_dot_ni=dot_product(n,n_i)
        drdx_dot_t=dot_product(drdx,t)
        drdx_dot_ti=dot_product(drdx,t_i)
        n_dot_ti=dot_product(n,t_i)
        t_dot_ni=dot_product(t,n_i)
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
        ! Functional shape functions (primary variables)
#       define etype e%ptype
#       define delta e%ptype_delta
#       define phi pphi
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Functional shape functions (secondary variables)
#       define etype e%stype
#       define delta e%stype_delta
#       define phi sphi
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Jacobian * weight
        jw=jg*js*w
        ! Add the integration point
        do im=1,2
          do iq=1,e%dme_n_gnodes
            ! L^s
            l1(iq,:,im)=l1(iq,:,im)-dr1*drdni*dme_vq(iq)*t(im)*sphi(:)*jw
            ! L^r
            l1(iq,:,im)=l1(iq,:,im)+2.d0*dr2*drdni*(dme_gphi(iq)-dme_gphi_i(iq))*drdx(im)*sphi(:)*jw
            l1(iq,:,im)=l1(iq,:,im)+n_i(im)*dr2*(dme_gphi(iq)-dme_gphi_i(iq)-dot_product(dme_wqj_i(iq,:),rv))*sphi(:)*jw
            l1(iq,:,im)=l1(iq,:,im)+n_i(im)*dr1*(dot_product(dme_wqj_i(iq,:),drdx)*sphi(:)-drdti*dme_vq_i(iq)*sphi_i(:))*jw
            ! L^ni
            l1(iq,:,im)=l1(iq,:,im)-dr1*(drdx_dot_ti*sphi(:)-drdti*sphi_i(:))*n_i(im)*dme_vq_i(iq)*jw

            ! M^s
            m1(iq,:,im)=m1(iq,:,im)+2.d0*dr2*drdn*drdni*dme_vq(iq)*t(im)*pphi(:)*jw
            m1(iq,:,im)=m1(iq,:,im)+dr2*n_dot_ni*dme_vq(iq)*t(im)*(pphi(:)-pphi_i(:)-dpphidx_i(:,1)*rv(1)-dpphidx_i(:,2)*rv(2))*jw

                                                                                          ! here  ,,,, it should be the second derivative dme field
            m1(iq,:,im)=m1(iq,:,im)+pphi_i(:)*dr2*n_dot_ni*(dme_vq(iq)*t(im)-dme_vq_i(iq)*t_i(im)-0.d0*r)*jw
            m1(iq,:,im)=m1(iq,:,im)+pphi_i(:)*dme_vq_i(:)*t_i(im)*dr2*(n_dot_ni-abs(drdt))*jw

            ! it gives zero because null second derivative of dme field...
            !m1(iq,:,im)=m1(iq,:,im)+pphi_i(:)*dr1*(n_dot_ni*dme_d2gphidxdx(iq,)...)
            !m1(iq,:,im)=m1(iq,:,im)+pphi_i(:)*dme_d2gphi

            m1(iq,:,im)=m1(iq,:,im)+dr1*(n_dot_ni*dme_vq(iq)*t(im)*(dpphidx_i(:,1)*drdx(1)+dpphidx_i(:,2)*drdx(2))-drdti*dme_vq_i(iq)*t_i(im)*(dpphidx_i(:,1)*t_i(1)+dpphidx_i(:,2)*t_i(2)))*jw
            m1(iq,:,im)=m1(iq,:,im)+dme_vq_i(iq)*t_i(im)*(dpphidx_i(:,1)*t_i(1)+dpphidx_i(:,2)*t_i(2))*dr1*(drdti-drdt)*jw

            ! M^r
            m1(iq,:,im)=m1(iq,:,im)+dr3*(-8.d0*drdn*drdni*drdx(im)+2.d0*(n(im)*drdni-n_i(im)*drdn))*(dme_gphi(iq)-dme_gphi_i(iq))*pphi(:)*jw
                                                                                                                            ! here  ,,,, it should be the second derivative dme field
            m1(iq,:,im)=m1(iq,:,im)-2.d0*dr3*n_dot_ni*drdx(im)*pphi(:)*(dme_gphi(iq)-dme_gphi_i(iq)-dot_product(dme_wqj_i(iq,:),rv)-0.d0*r)*jw
            m1(iq,:,im)=m1(iq,:,im)-2.d0*n_dot_ni*drdx(im)*dot_product(dme_wqj_i(iq,:),drdx)*dr2*(pphi(:)-pphi_i(:)-dpphidx_i(:,1)*rv(1)-dpphidx_i(:,2)*rv(2))*jw

            m1(iq,:,im)=m1(iq,:,im)-pphi_i(:)*2.d0*dr2*n_dot_ni*dot_product(dme_wqj_i(iq,:),drdx)*(drdx(im)-drdti*t_i(im)-0.5d0*dtdr_i(im)*r)*jw
            m1(iq,:,im)=m1(iq,:,im)-pphi_i(:)*t_i(im)*2.d0*dr2*n_dot_ni*drdti*(dot_product(dme_wqj_i(iq,:),drdx)-drdti*dme_vq_i(iq)-0.5d0*dot_product(dme_wqj_i(iq,:),dtdr_i)*r)*jw

            m1(iq,:,im)=m1(iq,:,im)-pphi_i(:)*t_i(im)*2.d0*dme_vq_i(iq)*dr2*(n_dot_ni-abs(drdt))*jw
            m1(iq,:,im)=m1(iq,:,im)-pphi_i(:)*t_i(im)*dot_product(dme_wqj_i(iq,:),dtdr_i)*dr1*(n_dot_ni-1.d0)*drdti*jw
            m1(iq,:,im)=m1(iq,:,im)-pphi_i(:)*t_i(im)*dot_product(dme_wqj_i(iq,:),dtdr_i)*dr1*(drdti-drdt)*jw
            m1(iq,:,im)=m1(iq,:,im)-pphi_i(:)*dtdr_i(im)*dr1*(n_dot_ni*dot_product(dme_wqj_i(iq,:),drdx)-drdti*dme_vq_i(iq))*jw
            m1(iq,:,im)=m1(iq,:,im)-pphi_i(:)*dtdr_i(im)*dme_vq_i(iq)*dr1*(drdti-drdt)*jw

            m1(iq,:,im)=m1(iq,:,im)-2.d0*dr1*(n_dot_ni*drdx(im)*(dpphidx_i(:,1)*drdx(1)+dpphidx_i(:,2)*drdx(2))*dot_product(dme_wqj_i(iq,:),drdx)-drdti*t_i(im)*(dpphidx_i(:,1)*t_i(1)+dpphidx_i(:,2)*t_i(2))*dme_vq_i(iq))*jw
            m1(iq,:,im)=m1(iq,:,im)-2.d0*t_i(im)*(dpphidx_i(:,1)*t_i(1)+dpphidx_i(:,2)*t_i(2))*dme_vq_i(iq)*dr1*(drdti-drdt)*jw

            ! falta integrales segunda derivada dme field
            !m1(iq,:,im)=m1(iq,:,im)

            ! M^n
            m1(iq,:,im)=m1(iq,:,im)-dr2*((n_dot_ti+t_dot_ni-2.d0*(drdx_dot_ti*drdn-drdx_dot_t*drdni))*n_i(im)*dme_vq_i(iq)&
                                        +(t_dot_ni+2.d0*drdt*drdni)*(n(im)*dme_vq(iq)-n_i(im)*dme_vq_i(iq)))*pphi(:)*jw

          end do
        end do
      end do
    end do
    write(50,*)
    !
    ! Analytical integrals
    !
    ra=sqrt((e%x(1,1)-x_i(1))**2+(e%x(2,1)-x_i(2))**2)
    rb=sqrt((e%x(1,2)-x_i(1))**2+(e%x(2,2)-x_i(2))**2)
    Idr1dr=log(rb)-log(ra)
    Idr2sgndr=-(1.d0/ra+1.d0/rb)
    do im=1,2
      do iq=1,e%dme_n_gnodes
        m1(iq,:,im)=m1(iq,:,im)-dme_vq_i(iq)*t_i(im)*pphi_i(:)*Idr2sgndr
        ! Falta termino segunda derivada de dme funcion
        !m1(iq,:,im)=m1(iq,:,im)+ ... * Idr1dr
        m1(iq,:,im)=m1(iq,:,im)-t_i(im)*dme_vq_i(iq)*(dpphidx_i(:,1)*t_i(1)+dpphidx_i(:,2)*t_i(2))*Idr1dr
        m1(iq,:,im)=m1(iq,:,im)-pphi_i(:)*dtdr_i(im)*dme_vq_i(iq)*Idr1dr
        m1(iq,:,im)=m1(iq,:,im)-pphi_i(:)*t_i(im)*dot_product(dme_wqj_i(iq,:),dtdr_i)*Idr1dr
      end do
    end do
    ! Multiply by constants
    m1=c_1_2pi*m1
    l1=c_1_2pi*l1
    ! Reverse element orientation
    if (reverse) l1=-l1
  end subroutine fbem_bem_stapot2d_vhbie_int

  ! ================================================================================================================================

end module fbem_bem_stapot2d

