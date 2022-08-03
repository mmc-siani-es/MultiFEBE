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
!! <b> This module implements the calculation of BEM integrals for 3D Laplace problems.</b>
module fbem_bem_stapot3d

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

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all are private
  private

  ! ================================================================================================================================
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE)
  ! Free-term c
  public :: fbem_bem_pot3d_sbie_freeterm
  ! Exterior integration
  public :: fbem_bem_stapot3d_sbie_ext_pre
  public :: fbem_bem_stapot3d_sbie_ext_st
  public :: fbem_bem_stapot3d_sbie_ext_adp
  ! Interior integration
  public :: fbem_bem_stapot3d_sbie_int
  ! Automatic integration
  public :: fbem_bem_stapot3d_sbie_auto
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (HBIE)
  ! Exterior integration
  public :: fbem_bem_stapot3d_hbie_ext_pre
  public :: fbem_bem_stapot3d_hbie_ext_st
  public :: fbem_bem_stapot3d_hbie_ext_adp
  ! Interior integration
  public :: fbem_bem_stapot3d_hbie_int
  public :: fbem_bem_stapot3d_hbie_int_li
  public :: fbem_bem_stapot3d_hbie_int_li1
  public :: fbem_bem_stapot3d_hbie_int_li2
  ! Automatic integration
  public :: fbem_bem_stapot3d_hbie_auto
  ! ================================================================================================================================

contains

  ! ================================================================================================================================
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE)

  !! Calculation of the BEM elastic free-term based on this amazing work:
  !!
  !! V. Mantic, "A new formula for the C-matrix in the Somigliana identity", Journal of Elasticity 33, 1993.
  !!
  !! In the implementation, "element boundary tangent" t is the unit tangent of the boundary of the element at a given point, where
  !! the unit tangent is taken always forward.
  !!
  !!>
  !! GENERAL
  !! Note: (+) is the point where t is calculated
  !!
  !! Edge point:
  !!
  !! |         t     |
  !! |       ---->   |
  !! +------(+)------+
  !!
  !! Corner point:
  !!
  !!        ^ |
  !!      t | |
  !!        | |
  !! --------(+)
  !!
  !! For our elements:
  !!
  !! QUADRILATERAL ELEMENTS (in the ilustration a QUAD9 is used)
  !!
  !!    xi2
  !!     ^
  !!     |
  !! 4---7---3
  !! |   |   |
  !! 8   9---6---> xi1
  !! |       |
  !! 1---5---2
  !!
  !! Node 1: t= (dx/dxi1)/|dx/dxi1|
  !! Node 5: t= (dx/dxi1)/|dx/dxi1|
  !! Node 2: t= (dx/dxi2)/|dx/dxi2|
  !! Node 6: t= (dx/dxi2)/|dx/dxi2|
  !! Node 3: t=-(dx/dxi1)/|dx/dxi1|
  !! Node 7: t=-(dx/dxi1)/|dx/dxi1|
  !! Node 4: t=-(dx/dxi2)/|dx/dxi2|
  !! Node 8: t=-(dx/dxi2)/|dx/dxi2|
  !!
  !! TRIANGULAR (in the ilustration a TRI6 is used)
  !!
  !! xi2
  !!  ^
  !!  |
  !!  2
  !!  |\
  !!  | \
  !!  |  \
  !!  5   4
  !!  |    \
  !!  |     \
  !!  |      \
  !!  3---6---1---> xi1
  !!           \
  !!            *
  !!           xi3
  !!
  !! Node 1: t=-(dx/dxi3)/|dx/dxi3|
  !! Node 4: t=-(dx/dxi3)/|dx/dxi3|
  !! Node 2: t=-(dx/dxi2)/|dx/dxi2|
  !! Node 5: t=-(dx/dxi2)/|dx/dxi2|
  !! Node 3: t= (dx/dxi1)/|dx/dxi1|
  !! Node 6: t= (dx/dxi1)/|dx/dxi1|
  !!
  !!<
  subroutine fbem_bem_pot3d_sbie_freeterm(n_elements,n,t,tol,c)
    implicit none
    ! Subroutine I/O variables
    integer            :: n_elements            !! Number of elements
    real(kind=real64)  :: n(3,n_elements)       !! Unit normal of each element
    real(kind=real64)  :: t(3,n_elements)       !! Unit element boundary tangent of each element
    real(kind=real64)  :: tol                   !! Geometric tolerance [1.0e-12,1.0e-3] (default 1.0e-6)
    real(kind=real64)  :: c                     !! Free-term
    ! Local variables
    integer            :: kc, ki, kj            ! Counters
    real(kind=real64)  :: ltol                  ! Geometric tolerance (local copy)
    real(kind=real64)  :: ln(3,0:n_elements+1)  ! Unit normal of each element (local copy)
    real(kind=real64)  :: lt(3,0:n_elements+1)  ! Unit element boundary tangent of each element (local copy)
    real(kind=real64)  :: e1(3), e2(3), e3(3)   ! Directional cosines
    real(kind=real64)  :: lti(3,n_elements)     ! lt in the system of coordinates (t_i, n_i x t_i, n_i)
    integer            :: n_t_common            ! Number of tangents in the same plane as the element i
    integer            :: t_common(n_elements)  ! The index of each common tangent
    real(kind=real64)  :: theta(n_elements)     ! The angle of each common tangent
    real(kind=real64)  :: mintheta              ! The minimum angle of all common tangent
    real(kind=real64)  :: tmp                   ! Temporary variable
    integer            :: minkj                 ! The index of the element with the minimum angle
    real(kind=real64)  :: nxn(3), nxndr, ndn
    real(kind=real64)  :: sum_a
    !
    ! Check geometric tolerance
    !
    if ((tol.lt.1.0d-12).or.(tol.gt.1.0d-3)) then
      ltol=1.0d-6
    else
      ltol=tol
    end if
    !
    ! Copy n and t to local copies ln and lt, respectively
    !
    do ki=1,n_elements
      do kc=1,3
        ln(kc,ki)=n(kc,ki)
        lt(kc,ki)=t(kc,ki)
      end do
    end do
    !
    ! Sort all the elements
    !
    ! Loop through elements (element i)
    do ki=1,n_elements-1
      !
      ! Express the tangents in the system of coordinates formed by (t_i, n_i x t_i, n_i)
      !
      ! Build directional cosines
      do kc=1,3
        e1(kc)=lt(kc,ki)
        e3(kc)=ln(kc,ki)
      end do
      e2(1)=e3(2)*e1(3)-e3(3)*e1(2)
      e2(2)=e3(3)*e1(1)-e3(1)*e1(3)
      e2(3)=e3(1)*e1(2)-e3(2)*e1(1)
      ! Transform lt to lti
      do kj=1,n_elements
        lti(1,kj)=e1(1)*lt(1,kj)+e1(2)*lt(2,kj)+e1(3)*lt(3,kj)
        lti(2,kj)=e2(1)*lt(1,kj)+e2(2)*lt(2,kj)+e2(3)*lt(3,kj)
        lti(3,kj)=e3(1)*lt(1,kj)+e3(2)*lt(2,kj)+e3(3)*lt(3,kj)
      end do
      !
      ! Find the tangents of other elements that are in the same plane as the element i
      !
      n_t_common=0
      do kj=ki+1,n_elements
        if (dabs(lti(3,kj)).le.ltol) then
          n_t_common=n_t_common+1
          t_common(n_t_common)=kj
        end if
      end do
      if (n_t_common.eq.0) then
        do kj=1,n_elements
          write(error_unit,*) 'Local element:', kj
          write(error_unit,*) 'n = ', n(:,kj)
          write(error_unit,*) 't = ', t(:,kj)
        end do
        stop 'the normals/tangents configuration is not valid'
      end if
      !
      ! Find the tangent that is next to the tangent i, which is in e1 in the new system of coordinates
      !
      ! Calculate the angle [0,2*pi] of each tangent j in the same plane as the element i
      do kj=1,n_t_common
        theta(t_common(kj))=datan2(lti(2,t_common(kj)),lti(1,t_common(kj)))
        if (theta(t_common(kj)).lt.0.0d0) theta(t_common(kj))=theta(t_common(kj))+2.0d0*c_pi
      end do
      ! Find the element with the minimum theta
      mintheta=theta(t_common(1))
      minkj=t_common(1)
      do kj=2,n_t_common
        if (theta(t_common(kj)).lt.mintheta) then
          mintheta=theta(t_common(kj))
          minkj=t_common(kj)
        end if
      end do
      !
      ! Swap the found element with the element i+1
      !
      do kc=1,3
        tmp=lt(kc,ki+1)
        lt(kc,ki+1)=lt(kc,minkj)
        lt(kc,minkj)=tmp
        tmp=ln(kc,ki+1)
        ln(kc,ki+1)=ln(kc,minkj)
        ln(kc,minkj)=tmp
      end do
    end do
    ! Save the element 0 as the element n_elements
    do kc=1,3
      ln(kc,0)=ln(kc,n_elements)
      lt(kc,0)=lt(kc,n_elements)
    end do
    ! Save the element n_elements+1 as the element 1
    do kc=1,3
      ln(kc,n_elements+1)=ln(kc,1)
      lt(kc,n_elements+1)=lt(kc,1)
    end do
    !
    ! Calculation
    !
    sum_a=0.0d0
    do ki=1,n_elements
      ! nxn = n_{i-1} x n_{i}
      nxn(1)=ln(2,ki-1)*ln(3,ki)-ln(3,ki-1)*ln(2,ki)
      nxn(2)=ln(3,ki-1)*ln(1,ki)-ln(1,ki-1)*ln(3,ki)
      nxn(3)=ln(1,ki-1)*ln(2,ki)-ln(2,ki-1)*ln(1,ki)
      ! nxn · t_{i}
      nxndr=nxn(1)*lt(1,ki)+nxn(2)*lt(2,ki)+nxn(3)*lt(3,ki)
      ! sgn(nxndr)
      if (nxndr.lt.0.0d0) nxndr=-1.0d0
      if (nxndr.gt.0.0d0) nxndr= 1.0d0
      ! ndn = n_{i-1} · n_{i}
      ndn=ln(1,ki-1)*ln(1,ki)+ln(2,ki-1)*ln(2,ki)+ln(3,ki-1)*ln(3,ki)
      ! To prevent errors with acos
      if (ndn.gt.1.d0) ndn=1.d0
      ! sgn(nxndr)*acos(ndn)
      sum_a=sum_a+nxndr*dacos(ndn)
    end do
    c=1.0d0/(4.0d0*c_pi)*(2.0d0*c_pi+sum_a)
  end subroutine fbem_bem_pot3d_sbie_freeterm

  subroutine fbem_bem_stapot3d_sbie_ext_pre(ps,e,reverse,x_i,h,g)
    implicit none
    ! I/O
    integer                :: ps            !! Selected precalculated dataset
    type(fbem_bem_element) :: e             !! Element
    logical                :: reverse       !! Reverse normal vector
    real(kind=real64)      :: x_i(3)        !! Position vector of the collocation point
    real(kind=real64)      :: h(e%n_pnodes) !! h vector
    real(kind=real64)      :: g(e%n_snodes) !! g vector
    ! Local
    integer           :: kip                 ! Counter
    real(kind=real64) :: x(3), n(3)          ! Position and unit normal vectors at the integration point
    real(kind=real64) :: rv(3), r, dr1, dr2  ! Distance vector and its modulus, 1/r and 1/r^2
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
      dr2=dr1**2
      drdn=dot_product(rv,n)*dr1
      pphijw=e%ps_pphijw(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      ! Add the integration point
      h=h+dr2*drdn*pphijw
      g=g+dr1*sphijw
    end do
    ! Multiply by constants
    h=-h*c_1_4pi
    g= g*c_1_4pi
    ! Reverse element orientation
    if (reverse) h=-h
  end subroutine fbem_bem_stapot3d_sbie_ext_pre

  subroutine fbem_bem_stapot3d_sbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,gln,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element) :: e                                !! Integration element
    logical                :: reverse                          !! Reverse normal vector
    real(kind=real64)      :: xi_s(2,fbem_n_vertices(e%gtype)) !! Coordinates of the subdivision of the parent element
    real(kind=real64)      :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)      :: barxip(2)                        !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)      :: barr                             !! Telles jacobian at barxip
    integer                :: gln                              !! Gauss-Legendre number of integration points (<=32)
    real(kind=real64)      :: h(e%n_pnodes)                    !! h vector
    real(kind=real64)      :: g(e%n_snodes)                    !! g vector
    ! Local
    integer                      :: kphi, k1, k2             ! Counters
    real(kind=real64)            :: aux(10)                  ! Auxiliary variable
    real(kind=real64)            :: gphi(e%n_gnodes)         ! Geometrical shape functions
    real(kind=real64)            :: dgphidxi1(e%n_gnodes)    ! Geometrical shape functions derivatives
    real(kind=real64)            :: dgphidxi2(e%n_gnodes)    ! Geometrical shape functions derivatives
    real(kind=real64)            :: pphi(e%n_pnodes)         ! Functional shape functions values
    real(kind=real64)            :: sphi(e%n_snodes)         ! Functional shape functions values
    real(kind=real64)            :: gamma(2)                 ! Vector of gamma_1,gamma_2 coordinates
    real(kind=real64)            :: w(2)                     ! Weights of the integration rule
    real(kind=real64)            :: xip(2)                   ! Vector of xip_1,xip_2 coordinates
    real(kind=real64)            :: phi_subquad(4)           ! Quadrilateral subdivision shape functions
    real(kind=real64)            :: dphidxi1_subquad(4)      ! Quadrilateral subdivision shape functions derivatives
    real(kind=real64)            :: dphidxi2_subquad(4)      ! Quadrilateral subdivision shape functions derivatives
    real(kind=real64)            :: phi_subtri(3)            ! Triangular subdivision shape functions
    real(kind=real64)            :: dphidxi1_subtri(3)       ! Triangular subdivision shape functions derivatives
    real(kind=real64)            :: dphidxi2_subtri(3)       ! Triangular subdivision shape functions derivatives
    real(kind=real64)            :: dxidxi1p(2)              ! xi derivatives with respect to xi1p
    real(kind=real64)            :: dxidxi2p(2)              ! xi derivatives with respect to xi2p
    real(kind=real64)            :: js                       ! Subdivision jacobian
    real(kind=real64)            :: xi(2)                    ! Vector of xi_1,xi_2 coordinates
    real(kind=real64)            :: xipp(2)                  ! Coordinate xipp used for quadrilateral-triangle transformation
    real(kind=real64)            :: barxipp(2)               ! Coordinate xipp of collocation point
    real(kind=real64)            :: jqt                      ! Jacobian of the quadrilateral-triangle transformation
    type(fbem_telles_parameters) :: telles_parameters(2)     ! Telles parameters for each coordinate
    real(kind=real64)            :: jt(2)                    ! Telles jacobian for each coordinate: xi_1->gamma_1 and xi_2->gamma_2
    real(kind=real64)            :: x(3), T1(3), T2(3), N(3) ! Position, tangents and normal vectors
    real(kind=real64)            :: jg                       ! Geometric jacobian
    real(kind=real64)            :: rv(3), r, dr1, dr2       ! Distance vector
    real(kind=real64)            :: drdn                     ! Partial derivative of r respect to unit normal
    real(kind=real64)            :: jw                       ! Jacobians * weights
    real(kind=real64)            :: pphijw(e%n_pnodes)       ! shape functions (primary variable) * jacobian * weight
    real(kind=real64)            :: sphijw(e%n_snodes)       ! shape functions (secondary variable) * jacobian * weight
    ! Initialization
    h=0.d0
    g=0.d0
    ! Numerical integration
    select case (fbem_n_vertices(e%gtype))
      ! QUADRILATERAL ELEMENTS
      case (4)
        ! Calculate Telles parameters
        telles_parameters(1)=fbem_telles11_calculate_parameters(barxip(1),barr)
        telles_parameters(2)=fbem_telles11_calculate_parameters(barxip(2),barr)
        ! Loop through integration points (gamma1)
        do k1=1,gl11_n(gln)
          ! gamma1
          gamma(1)=gl11_xi(k1,gln)
          w(1)=gl11_w(k1,gln)
          ! gamma1 -> xi1p
          call fbem_telles_xi_and_jacobian(telles_parameters(1),gamma(1),xip(1),jt(1))
          ! Loop through integration points (gamma2)
          do k2=1,gl11_n(gln)
            ! gamma2
            gamma(2)=gl11_xi(k2,gln)
            w(2)=gl11_w(k2,gln)
            ! gamma2 -> xi2p
            call fbem_telles_xi_and_jacobian(telles_parameters(2),gamma(2),xip(2),jt(2))
            ! xip -> xi
#           define delta 0.0d0
#           define xi xip
#           define phi phi_subquad
#           define dphidxi1 dphidxi1_subquad
#           define dphidxi2 dphidxi2_subquad
#           include <phi_quad4.rc>
#           include <dphidxi1_quad4.rc>
#           include <dphidxi2_quad4.rc>
#           undef delta
#           undef xi
#           undef phi
#           undef dphidxi1
#           undef dphidxi2
            ! Calculation of xi, dxidxi1p and dxidxi2p
            xi=0.d0
            dxidxi1p=0.d0
            dxidxi2p=0.d0
            do kphi=1,4
              xi=xi+phi_subquad(kphi)*xi_s(:,kphi)
              dxidxi1p=dxidxi1p+dphidxi1_subquad(kphi)*xi_s(:,kphi)
              dxidxi2p=dxidxi2p+dphidxi2_subquad(kphi)*xi_s(:,kphi)
            end do
            ! xip -> xi jacobian
            js=dxidxi1p(1)*dxidxi2p(2)-dxidxi1p(2)*dxidxi2p(1)
            ! xi -> x
#           define etype e%gtype
#           define delta 0.0d0
#           define phi gphi
#           define dphidxi1 dgphidxi1
#           define dphidxi2 dgphidxi2
#           include <phi_and_dphidxik_2d.rc>
#           undef etype
#           undef delta
#           undef phi
#           undef dphidxi1
#           undef dphidxi2
            ! Calculation of x, T1 and T2
            x=0.d0
            T1=0.d0
            T2=0.d0
            do kphi=1,e%n_gnodes
              x=x+gphi(kphi)*e%x(:,kphi)
              T1=T1+dgphidxi1(kphi)*e%x(:,kphi)
              T2=T2+dgphidxi2(kphi)*e%x(:,kphi)
            end do
            ! Normal vector is N = T1 x T2
            N(1)=T1(2)*T2(3)-T1(3)*T2(2)
            N(2)=T1(3)*T2(1)-T1(1)*T2(3)
            N(3)=T1(1)*T2(2)-T1(2)*T2(1)
            ! xi -> x jacobian
            jg=sqrt(dot_product(N,N))
            ! Unit normal
            n=N/jg
            ! Distance vector
            rv=x-x_i
            r=sqrt(dot_product(rv,rv))
            dr1=1.d0/r
            dr2=dr1**2
            ! dr/dn
            drdn=dot_product(rv,n)*dr1
            ! Functional shape functions (primary variables)
#           define etype e%ptype
#           define delta e%ptype_delta
#           define phi pphi
#           include <phi_2d.rc>
#           undef etype
#           undef delta
#           undef phi
            ! Functional shape functions (secondary variables)
#           define etype e%stype
#           define delta e%stype_delta
#           define phi sphi
#           include <phi_2d.rc>
#           undef etype
#           undef delta
#           undef phi
            ! Jacobians * weights
            jw=jg*js*jt(1)*jt(2)*w(1)*w(2)
            ! Functional shape functions * jacobians * weights
            pphijw=pphi*jw
            sphijw=sphi*jw
            ! Add the integration points
            h=h+dr2*drdn*pphijw
            g=g+dr1*sphijw
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
        ! Loop through integration points (gamma1)
        do k1=1,gl01_n(gln)
          ! gamma1
          gamma(1)=gl01_xi(k1,gln)
          w(1)=gl01_w(k1,gln)
          ! gamma1 -> xi1pp
          call fbem_telles_xi_and_jacobian(telles_parameters(1),gamma(1),xipp(1),jt(1))
          ! Loop through integration points (gamma2)
          do k2=1,gl01_n(gln)
            ! gamma2
            gamma(2)=gl01_xi(k2,gln)
            w(2)=gl01_w(k2,gln)
            ! gamma2 -> xi2pp
            call fbem_telles_xi_and_jacobian(telles_parameters(2),gamma(2),xipp(2),jt(2))
            ! xipp -> xip (QUAD -> TRI)
            xip(1)=(1.0d0-xipp(2))*xipp(1)
            xip(2)=xipp(2)
            jqt=1.d0-xipp(2)
            ! xip -> xi
#           define delta 0.0d0
#           define xi xip
#           define phi phi_subtri
#           define dphidxi1 dphidxi1_subtri
#           define dphidxi2 dphidxi2_subtri
#           include <phi_tri3.rc>
#           include <dphidxi1_tri3.rc>
#           include <dphidxi2_tri3.rc>
#           undef delta
#           undef xi
#           undef phi
#           undef dphidxi1
#           undef dphidxi2
            ! Calculation of xi, dxidxi1p and dxidxi2p
            xi=0.d0
            dxidxi1p=0.d0
            dxidxi2p=0.d0
            do kphi=1,3
              xi=xi+phi_subtri(kphi)*xi_s(:,kphi)
              dxidxi1p=dxidxi1p+dphidxi1_subtri(kphi)*xi_s(:,kphi)
              dxidxi2p=dxidxi2p+dphidxi2_subtri(kphi)*xi_s(:,kphi)
            end do
            ! xip -> xi jacobian
            js=dxidxi1p(1)*dxidxi2p(2)-dxidxi1p(2)*dxidxi2p(1)
            ! xi -> x
#           define etype e%gtype
#           define delta 0.0d0
#           define phi gphi
#           define dphidxi1 dgphidxi1
#           define dphidxi2 dgphidxi2
#           include <phi_and_dphidxik_2d.rc>
#           undef etype
#           undef delta
#           undef phi
#           undef dphidxi1
#           undef dphidxi2
            ! Calculation of x, T1 and T2
            x=0.d0
            T1=0.d0
            T2=0.d0
            do kphi=1,e%n_gnodes
              x=x+gphi(kphi)*e%x(:,kphi)
              T1=T1+dgphidxi1(kphi)*e%x(:,kphi)
              T2=T2+dgphidxi2(kphi)*e%x(:,kphi)
            end do
            ! Normal vector is N = T1 x T2
            N(1)=T1(2)*T2(3)-T1(3)*T2(2)
            N(2)=T1(3)*T2(1)-T1(1)*T2(3)
            N(3)=T1(1)*T2(2)-T1(2)*T2(1)
            ! xi -> x jacobian
            jg=sqrt(dot_product(N,N))
            ! Unit normal
            n=N/jg
            ! Distance vector
            rv=x-x_i
            r=sqrt(dot_product(rv,rv))
            dr1=1.d0/r
            dr2=dr1**2
            ! dr/dn
            drdn=dot_product(rv,n)*dr1
            ! Functional shape functions (primary variables)
#           define etype e%ptype
#           define delta e%ptype_delta
#           define phi pphi
#           include <phi_2d.rc>
#           undef etype
#           undef delta
#           undef phi
            ! Functional shape functions (secondary variables)
#           define etype e%stype
#           define delta e%stype_delta
#           define phi sphi
#           include <phi_2d.rc>
#           undef etype
#           undef delta
#           undef phi
            ! Jacobians * weights
            jw=jg*js*jqt*jt(1)*jt(2)*w(1)*w(2)
            ! Functional shape functions * jacobians * weights
            pphijw=pphi*jw
            sphijw=sphi*jw
            ! Add the integration points
            h=h+dr2*drdn*pphijw
            g=g+dr1*sphijw
          end do
        end do
      case default
        call fbem_error_message(error_unit,0,'fbem_bem_stapot3d_sbie_ext_st',0,'n_edges not valid')
    end select
    ! Multiply by constants
    h=-h*c_1_4pi
    g= g*c_1_4pi
    ! Reverse element orientation
    if (reverse) h=-h
  end subroutine fbem_bem_stapot3d_sbie_ext_st

  recursive subroutine fbem_bem_stapot3d_sbie_ext_adp(e,reverse,xi_s,x_i,qsp,ks,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                !! Element
    logical                  :: reverse                          !! Reverse orientation
    real(kind=real64)        :: xi_s(2,fbem_n_vertices(e%gtype)) !! Subdivision of the parent element
    real(kind=real64)        :: x_i(3)                           !! Collocation point position vector
    type(fbem_qs_parameters) :: qsp                              !! Quasi-singular integration parameters
    integer                  :: ks                               !! Current level of subdivisions
    integer                  :: ns                               !! Maximum level of subdivision
    real(kind=real64)        :: h(e%n_pnodes)                    !! h vector
    real(kind=real64)        :: g(e%n_snodes)                    !! g vector
    ! Local
    integer           :: gln_near                             ! Gauss-Legendre integ. points required to integrate only  the quasi-singular integrand
    integer           :: gln                                  ! Gauss-Legendre integ. points required to integrate the whole integrand
    logical           :: subdivide                            ! True if subdivision has to be performed
    real(kind=real64) :: barxi(2)                             ! Nearest element coordinates with respect to collocation point
    real(kind=real64) :: barxip(2)                            ! Nearest element subdivision local coordinates with respect to collocation point
    real(kind=real64) :: rmin                                 ! Minimum distance between collocation point and barxi on the element
    real(kind=real64) :: barr                                 ! Telles jacobian at barxi
    real(kind=real64) :: cl                                   ! Characteristic length
    real(kind=real64) :: d                                    ! Normalized distance between collocation point and element subdivision
    integer           :: method                               ! Method used in nearest point algorithm
    real(kind=real64) :: tmp_xi_s(2,fbem_n_vertices(e%gtype)) ! Subdivision
    real(kind=real64) :: x_s(3,e%n_gnodes)                    ! Coordinates of the element subdivision
    real(kind=real64) :: h_tmp(e%n_pnodes)                    ! h vector (temporary)
    real(kind=real64) :: g_tmp(e%n_snodes)                    ! g vector (temporary)
    ! Initialize
    if (ks.eq.1) then
      h=0.d0
      g=0.d0
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
      cl=fbem_characteristic_length(3,e%gtype,x_s,1.d-6)
      call fbem_nearest_element_point_bem(3,e%gtype,x_s,cl,x_i,barxip,rmin,d,method)
    end if
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,3,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_stapot3d_sbie_ext_adp',ns,'maximum number of subdivisions reached')
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
          call fbem_bem_stapot3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,qsp,ks+1,ns,h,g)
          ! SUBTRI 2
          tmp_xi_s(:,1)=xi_s(:,2)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          call fbem_bem_stapot3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,qsp,ks+1,ns,h,g)
          ! SUBTRI 3
          tmp_xi_s(:,1)=xi_s(:,3)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          call fbem_bem_stapot3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,qsp,ks+1,ns,h,g)
          ! SUBTRI 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          call fbem_bem_stapot3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,qsp,ks+1,ns,h,g)
        ! QUADRILATERALS
        case (4)
          ! SUBQUAD 1
          tmp_xi_s(:,1)=xi_s(:,1)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,3)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          call fbem_bem_stapot3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,qsp,ks+1,ns,h,g)
          ! SUBQUAD 2
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=xi_s(:,2)
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,4)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          call fbem_bem_stapot3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,qsp,ks+1,ns,h,g)
          ! SUBQUAD 3
          tmp_xi_s(:,1)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=xi_s(:,3)
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          call fbem_bem_stapot3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,qsp,ks+1,ns,h,g)
          ! SUBQUAD 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          tmp_xi_s(:,2)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=xi_s(:,4)
          call fbem_bem_stapot3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,qsp,ks+1,ns,h,g)
      end select
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_stapot3d_sbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,gln,h_tmp,g_tmp)
      h=h+h_tmp
      g=g+g_tmp
    end if
  end subroutine fbem_bem_stapot3d_sbie_ext_adp

  subroutine fbem_bem_stapot3d_sbie_int(e,reverse,xi_i,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element) :: e             !! Element
    logical                :: reverse       !! Normal vector inversion
    real(kind=real64)      :: xi_i(2)       !! Reference coordinates of the singular point.
    real(kind=real64)      :: h(e%n_pnodes) !! h vector
    real(kind=real64)      :: g(e%n_snodes) !! g vector
    ! Local
    integer           :: ksubtri                          ! Counter variable for subtriangles loop
    integer           :: nsubtri                          ! Number of subtriangles
    integer           :: subtriangle(8)                   ! Vector that contains what subtriangles need to be integrated
    real(kind=real64) :: theta_subtri(2,8)                ! Matrix that contains the angles theta of the subtriangles to be integrated
    real(kind=real64) :: thetap_subtri(2,8)               ! Matrix that contains the angles thetap of the subtriangles to be integrated
    integer           :: ktheta                           ! Counter variable for theta coordinate loop
    integer           :: krho                             ! Counter variable for rho coordinate loop
    integer           :: kphi                             ! Counter coordinates for shape functions loops
    integer           :: ngp_theta                        ! Number of Gauss points for theta coordinate
    integer           :: ngp_rho                          ! Number of Gauss points for rho coordinate
    real(kind=real64) :: thetai, thetaf, thetapi, thetapf ! Initial and final angles for subtriangle integration
    real(kind=real64) :: w_angular                        ! Weight of the angular coordinate
    real(kind=real64) :: w_radial                         ! Weight of the radial coordinate
    real(kind=real64) :: theta                            ! Angle coordinate theta
    real(kind=real64) :: thetap                           ! Angle coordinate thetap
    real(kind=real64) :: thetapp                          ! Angle coordinate thetap on [0,1] domain
    real(kind=real64) :: jthetap                          ! thetap->thetapp jacobian
    real(kind=real64) :: rhoij                            ! Maximum rho (radial) value for each edge
    real(kind=real64) :: rho                              ! Radial coordinate rho
    real(kind=real64) :: rhop                             ! Radial coordinate rho on [0,1] domain
    real(kind=real64) :: aux(10)                          ! Auxiliary variable for shape functions resources
    real(kind=real64) :: xi(2)                            ! Reference xi_1,xi_2 coordinates
    real(kind=real64) :: gphi(e%n_gnodes)                 ! Geometrical shape functions values at xi_1,xi_2
    real(kind=real64) :: dgphidxi1(e%n_gnodes)            ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64) :: dgphidxi2(e%n_gnodes)            ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64) :: pphi(e%n_pnodes)                 ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: sphi(e%n_snodes)                 ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: x_i(3)                           ! Position vector at the collocation point
    real(kind=real64) :: x(3), T1(3), T2(3), N(3)         ! Position, tangents and normal vectors
    real(kind=real64) :: jg                               ! Geometric jacobian
    real(kind=real64) :: rv(3), r, dr1, dr2               ! Distance vector
    real(kind=real64) :: jw                               ! Jacobians * weights
    real(kind=real64) :: drdn                             ! dr/dn
    real(kind=real64) :: pphijw(e%n_pnodes)               ! Functional shape functions * jw
    real(kind=real64) :: sphijw(e%n_snodes)               ! Functional shape functions * jw
    real(kind=real64) :: costheta, sintheta               ! cos(theta), sin(theta)
    ! Initialization
    h=0.d0
    g=0.d0
    ! Calculate x_i
#   define etype e%gtype
#   define delta 0.0d0
#   define xi xi_i
#   define phi gphi
#   include <phi_2d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
    x_i=0.d0
    do kphi=1,e%n_gnodes
      x_i=x_i+gphi(kphi)*e%x(:,kphi)
    end do
    ! Setup of the polar transformation
    call fbem_polar_transformation_setup(e%gtype,xi_i,nsubtri,subtriangle,theta_subtri,thetap_subtri)
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
      ! Select the number of Gauss points
      ngp_rho=15
      ngp_theta=5+nint(25.0d0*(thetaf-thetai)/c_pi_2)
      ! Loop through angular coordinate
      do ktheta=1,gl01_n(ngp_theta)
        thetapp=gl01_xi(ktheta,ngp_theta)
        w_angular=gl01_w(ktheta,ngp_theta)
        jthetap=(thetapf-thetapi)
        thetap=jthetap*thetapp+thetapi
        ! Angular transformation
        call fbem_polar_transformation_angular(e%gtype,xi_i,subtriangle(ksubtri),thetap,theta,rhoij)
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
          ! xi -> x
#         define etype e%gtype
#         define delta 0.0d0
#         define phi gphi
#         define dphidxi1 dgphidxi1
#         define dphidxi2 dgphidxi2
#         include <phi_and_dphidxik_2d.rc>
#         undef etype
#         undef delta
#         undef phi
#         undef dphidxi1
#         undef dphidxi2
          ! Calculation of x, T1 and T2
          x=0.d0
          T1=0.d0
          T2=0.d0
          do kphi=1,e%n_gnodes
            x=x+gphi(kphi)*e%x(:,kphi)
            T1=T1+dgphidxi1(kphi)*e%x(:,kphi)
            T2=T2+dgphidxi2(kphi)*e%x(:,kphi)
          end do
          ! Normal vector is N = T1 x T2
          N(1)=T1(2)*T2(3)-T1(3)*T2(2)
          N(2)=T1(3)*T2(1)-T1(1)*T2(3)
          N(3)=T1(1)*T2(2)-T1(2)*T2(1)
          ! Geometric jacobian
          jg=sqrt(dot_product(N,N))
          ! Unit normal
          n=N/jg
          ! Distance vector
          rv=x-x_i
          r=sqrt(dot_product(rv,rv))
          dr1=1.d0/r
          dr2=dr1**2
          ! dr/dn
          drdn=dot_product(rv,n)*dr1
          ! Functional shape functions (primary variable)
#         define etype e%ptype
#         define delta e%ptype_delta
#         define phi pphi
#         include <phi_2d.rc>
#         undef etype
#         undef delta
#         undef phi
          ! Functional shape functions (secondary variable)
#         define etype e%stype
#         define delta e%stype_delta
#         define phi sphi
#         include <phi_2d.rc>
#         undef etype
#         undef delta
#         undef phi
          ! Jacobians * weights
          jw=jg*rho*jthetap*w_angular*w_radial
          ! Functional shape functions * jacobians * weights
          pphijw=pphi*jw
          sphijw=sphi*jw
          ! Add to kernels
          h=h+dr2*drdn*pphijw
          g=g+dr1*sphijw
        end do ! Loop through rho coordinate
      end do ! Loop through theta coordinate
    end do ! Loop through rectangle triangles
    ! Multiply by constants
    h=-h*c_1_4pi
    g= g*c_1_4pi
    ! Reverse element orientation
    if (reverse) h=-h
  end subroutine fbem_bem_stapot3d_sbie_int

  subroutine fbem_bem_stapot3d_sbie_auto(e,reverse,x_i,qsp,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e             !! Element
    logical                  :: reverse       !! Reverse orientation
    real(kind=real64)        :: x_i(3)        !! Collocation point
    type(fbem_qs_parameters) :: qsp           !! Quasi-singular integration parameters
    integer                  :: ns            !! Maximum level of subdivisions
    real(kind=real64)        :: h(e%n_pnodes) !! h vector
    real(kind=real64)        :: g(e%n_snodes) !! g vector
    ! Local
    real(kind=real64) :: r(3)                             ! Distance vector
    real(kind=real64) :: rmin                             ! Minimum distance between element and x_i
    real(kind=real64) :: barxi(2)                         ! Local coordinates of the nearest element point with respect to x_i
    real(kind=real64) :: d                                ! Dimensionless distance
    integer           :: delta                            ! Control variable
    real(kind=real64) :: xi_s(2,fbem_n_vertices(e%gtype)) ! Local coordinates of the element subdivision
    integer           :: method                           ! Method used when calculating the nearest element point
    integer           :: gln_near                         ! Gauss-Legendre integration points required by the quasi-singular function
    integer           :: gln                              ! Gauss-Legendre integration points used in the integration
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
        call fbem_bem_stapot3d_sbie_int(e,reverse,barxi,h,g)
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
          call fbem_bem_stapot3d_sbie_ext_pre(ps,e,reverse,x_i,h,g)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_stapot3d_sbie_ext_adp(e,reverse,xi_s,x_i,qsp,1,ns,h,g)
        end if
    end select
  end subroutine fbem_bem_stapot3d_sbie_auto

  ! ================================================================================================================================

  ! ================================================================================================================================
  ! HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (HBIE)

  subroutine fbem_bem_stapot3d_hbie_ext_pre(ps,e,reverse,x_i,n_i,m,l)
    implicit none
    ! I/O
    integer                :: ps            !! Selected precalculated dataset
    type(fbem_bem_element) :: e             !! Element
    logical                :: reverse       !! Reverse normal vector
    real(kind=real64)      :: x_i(3)        !! Position vector of the collocation point
    real(kind=real64)      :: n_i(3)        !! Unit normal vector of the collocation point
    real(kind=real64)      :: m(e%n_pnodes) !! m vector
    real(kind=real64)      :: l(e%n_snodes) !! l vector
    ! Local
    integer           :: kip                     ! Counter
    real(kind=real64) :: x(3), n(3)              ! Position and unit normal vectors at the integration point
    real(kind=real64) :: rv(3), r, dr1, dr2, dr3 ! Distance vector and its modulus, 1/r, 1/r^2 and 1/r^3
    real(kind=real64) :: drdn                    ! dr/dn
    real(kind=real64) :: drdni                   ! dr/dn^i
    real(kind=real64) :: n_dot_ni                ! n·n^i
    real(kind=real64) :: pphijw(e%n_pnodes)      ! shape functions (primary variable) * jacobian * weight
    real(kind=real64) :: sphijw(e%n_snodes)      ! shape functions (secondary variable) * jacobian * weight
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
      dr3=dr2*dr1
      drdn=dot_product(rv,n)*dr1
      drdni=-dot_product(rv,n_i)*dr1
      n_dot_ni=dot_product(n,n_i)
      pphijw=e%ps_pphijw(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      ! Add the integration point
      m=m+dr3*(3.d0*drdn*drdni+n_dot_ni)*pphijw
      l=l+dr2*drdni*sphijw
    end do
    ! Multiply by constants
    m= m*c_1_4pi
    l=-l*c_1_4pi
    ! Reverse element orientation
    if (reverse) m=-m
  end subroutine fbem_bem_stapot3d_hbie_ext_pre

  subroutine fbem_bem_stapot3d_hbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,gln,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element) :: e                                !! Integration element
    logical                :: reverse                          !! Reverse normal vector
    real(kind=real64)      :: xi_s(2,fbem_n_vertices(e%gtype)) !! Coordinates of the subdivision of the parent element
    real(kind=real64)      :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)      :: n_i(3)                           !! Unit normal vector of the collocation point
    real(kind=real64)      :: barxip(2)                        !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)      :: barr                             !! Telles jacobian at barxip
    integer                :: gln                              !! Gauss-Legendre number of integration points (<=32)
    real(kind=real64)      :: m(e%n_pnodes)                    !! m vector
    real(kind=real64)      :: l(e%n_snodes)                    !! l vector
    ! Local
    integer                      :: kphi, k1, k2             ! Counters
    real(kind=real64)            :: aux(10)                  ! Auxiliary variable
    real(kind=real64)            :: gphi(e%n_gnodes)         ! Geometrical shape functions
    real(kind=real64)            :: dgphidxi1(e%n_gnodes)    ! Geometrical shape functions derivatives
    real(kind=real64)            :: dgphidxi2(e%n_gnodes)    ! Geometrical shape functions derivatives
    real(kind=real64)            :: pphi(e%n_pnodes)         ! Functional shape functions values
    real(kind=real64)            :: sphi(e%n_snodes)         ! Functional shape functions values
    real(kind=real64)            :: gamma(2)                 ! Vector of gamma_1,gamma_2 coordinates
    real(kind=real64)            :: w(2)                     ! Weights of the integration rule
    real(kind=real64)            :: xip(2)                   ! Vector of xip_1,xip_2 coordinates
    real(kind=real64)            :: phi_subquad(4)           ! Quadrilateral subdivision shape functions
    real(kind=real64)            :: dphidxi1_subquad(4)      ! Quadrilateral subdivision shape functions derivatives
    real(kind=real64)            :: dphidxi2_subquad(4)      ! Quadrilateral subdivision shape functions derivatives
    real(kind=real64)            :: phi_subtri(3)            ! Triangular subdivision shape functions
    real(kind=real64)            :: dphidxi1_subtri(3)       ! Triangular subdivision shape functions derivatives
    real(kind=real64)            :: dphidxi2_subtri(3)       ! Triangular subdivision shape functions derivatives
    real(kind=real64)            :: dxidxi1p(2)              ! xi derivatives with respect to xi1p
    real(kind=real64)            :: dxidxi2p(2)              ! xi derivatives with respect to xi2p
    real(kind=real64)            :: js                       ! Subdivision jacobian
    real(kind=real64)            :: xi(2)                    ! Vector of xi_1,xi_2 coordinates
    real(kind=real64)            :: xipp(2)                  ! Coordinate xipp used for quadrilateral-triangle transformation
    real(kind=real64)            :: barxipp(2)               ! Coordinate xipp of collocation point
    real(kind=real64)            :: jqt                      ! Jacobian of the quadrilateral-triangle transformation
    type(fbem_telles_parameters) :: telles_parameters(2)     ! Telles parameters for each coordinate
    real(kind=real64)            :: jt(2)                    ! Telles jacobian for each coordinate: xi_1->gamma_1 and xi_2->gamma_2
    real(kind=real64)            :: x(3), T1(3), T2(3), N(3) ! Position, tangents and normal vectors
    real(kind=real64)            :: jg                       ! Geometric jacobian
    real(kind=real64)            :: rv(3), r, dr1, dr2,  dr3 ! Distance vector
    real(kind=real64)            :: drdn                     ! dr/dn
    real(kind=real64)            :: drdni                    ! dr/dn^i
    real(kind=real64)            :: n_dot_ni                 ! n·n^i
    real(kind=real64)            :: jw                       ! Jacobians * weights
    real(kind=real64)            :: pphijw(e%n_pnodes)       ! shape functions (primary variable) * jacobian * weight
    real(kind=real64)            :: sphijw(e%n_snodes)       ! shape functions (secondary variable) * jacobian * weight
    ! Initialization
    m=0.d0
    l=0.d0
    ! Numerical integration
    select case (fbem_n_vertices(e%gtype))
      ! QUADRILATERAL ELEMENTS
      case (4)
        ! Calculate Telles parameters
        telles_parameters(1)=fbem_telles11_calculate_parameters(barxip(1),barr)
        telles_parameters(2)=fbem_telles11_calculate_parameters(barxip(2),barr)
        ! Loop through integration points (gamma1)
        do k1=1,gl11_n(gln)
          ! gamma1
          gamma(1)=gl11_xi(k1,gln)
          w(1)=gl11_w(k1,gln)
          ! gamma1 -> xi1p
          call fbem_telles_xi_and_jacobian(telles_parameters(1),gamma(1),xip(1),jt(1))
          ! Loop through integration points (gamma2)
          do k2=1,gl11_n(gln)
            ! gamma2
            gamma(2)=gl11_xi(k2,gln)
            w(2)=gl11_w(k2,gln)
            ! gamma2 -> xi2p
            call fbem_telles_xi_and_jacobian(telles_parameters(2),gamma(2),xip(2),jt(2))
            ! xip -> xi
#           define delta 0.0d0
#           define xi xip
#           define phi phi_subquad
#           define dphidxi1 dphidxi1_subquad
#           define dphidxi2 dphidxi2_subquad
#           include <phi_quad4.rc>
#           include <dphidxi1_quad4.rc>
#           include <dphidxi2_quad4.rc>
#           undef delta
#           undef xi
#           undef phi
#           undef dphidxi1
#           undef dphidxi2
            ! Calculation of xi, dxidxi1p and dxidxi2p
            xi=0.d0
            dxidxi1p=0.d0
            dxidxi2p=0.d0
            do kphi=1,4
              xi=xi+phi_subquad(kphi)*xi_s(:,kphi)
              dxidxi1p=dxidxi1p+dphidxi1_subquad(kphi)*xi_s(:,kphi)
              dxidxi2p=dxidxi2p+dphidxi2_subquad(kphi)*xi_s(:,kphi)
            end do
            ! xip -> xi jacobian
            js=dxidxi1p(1)*dxidxi2p(2)-dxidxi1p(2)*dxidxi2p(1)
            ! xi -> x
#           define etype e%gtype
#           define delta 0.0d0
#           define phi gphi
#           define dphidxi1 dgphidxi1
#           define dphidxi2 dgphidxi2
#           include <phi_and_dphidxik_2d.rc>
#           undef etype
#           undef delta
#           undef phi
#           undef dphidxi1
#           undef dphidxi2
            ! Calculation of x, T1 and T2
            x=0.d0
            T1=0.d0
            T2=0.d0
            do kphi=1,e%n_gnodes
              x=x+gphi(kphi)*e%x(:,kphi)
              T1=T1+dgphidxi1(kphi)*e%x(:,kphi)
              T2=T2+dgphidxi2(kphi)*e%x(:,kphi)
            end do
            ! Normal vector is N = T1 x T2
            N(1)=T1(2)*T2(3)-T1(3)*T2(2)
            N(2)=T1(3)*T2(1)-T1(1)*T2(3)
            N(3)=T1(1)*T2(2)-T1(2)*T2(1)
            ! xi -> x jacobian
            jg=sqrt(dot_product(N,N))
            ! Unit normal
            n=N/jg
            ! Distance vector
            rv=x-x_i
            r=sqrt(dot_product(rv,rv))
            dr1=1.d0/r
            dr2=dr1**2
            dr3=dr2*dr1
            ! dr/dn
            drdn=dot_product(rv,n)*dr1
            ! dr/dni
            drdni=-dot_product(rv,n_i)*dr1
            ! n·n_i
            n_dot_ni=dot_product(n,n_i)
            ! Functional shape functions (primary variables)
#           define etype e%ptype
#           define delta e%ptype_delta
#           define phi pphi
#           include <phi_2d.rc>
#           undef etype
#           undef delta
#           undef phi
            ! Functional shape functions (secondary variables)
#           define etype e%stype
#           define delta e%stype_delta
#           define phi sphi
#           include <phi_2d.rc>
#           undef etype
#           undef delta
#           undef phi
            ! Jacobians * weights
            jw=jg*js*jt(1)*jt(2)*w(1)*w(2)
            ! Functional shape functions * jacobians * weights
            pphijw=pphi*jw
            sphijw=sphi*jw
            ! Add the integration points
            m=m+dr3*(3.d0*drdn*drdni+n_dot_ni)*pphijw
            l=l+dr2*drdni*sphijw
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
        ! Loop through integration points (gamma1)
        do k1=1,gl01_n(gln)
          ! gamma1
          gamma(1)=gl01_xi(k1,gln)
          w(1)=gl01_w(k1,gln)
          ! gamma1 -> xi1pp
          call fbem_telles_xi_and_jacobian(telles_parameters(1),gamma(1),xipp(1),jt(1))
          ! Loop through integration points (gamma2)
          do k2=1,gl01_n(gln)
            ! gamma2
            gamma(2)=gl01_xi(k2,gln)
            w(2)=gl01_w(k2,gln)
            ! gamma2 -> xi2pp
            call fbem_telles_xi_and_jacobian(telles_parameters(2),gamma(2),xipp(2),jt(2))
            ! xipp -> xip (QUAD -> TRI)
            xip(1)=(1.0d0-xipp(2))*xipp(1)
            xip(2)=xipp(2)
            jqt=1.d0-xipp(2)
            ! xip -> xi
#           define delta 0.0d0
#           define xi xip
#           define phi phi_subtri
#           define dphidxi1 dphidxi1_subtri
#           define dphidxi2 dphidxi2_subtri
#           include <phi_tri3.rc>
#           include <dphidxi1_tri3.rc>
#           include <dphidxi2_tri3.rc>
#           undef delta
#           undef xi
#           undef phi
#           undef dphidxi1
#           undef dphidxi2
            ! Calculation of xi, dxidxi1p and dxidxi2p
            xi=0.d0
            dxidxi1p=0.d0
            dxidxi2p=0.d0
            do kphi=1,3
              xi=xi+phi_subtri(kphi)*xi_s(:,kphi)
              dxidxi1p=dxidxi1p+dphidxi1_subtri(kphi)*xi_s(:,kphi)
              dxidxi2p=dxidxi2p+dphidxi2_subtri(kphi)*xi_s(:,kphi)
            end do
            ! xip -> xi jacobian
            js=dxidxi1p(1)*dxidxi2p(2)-dxidxi1p(2)*dxidxi2p(1)
            ! xi -> x
#           define etype e%gtype
#           define delta 0.0d0
#           define phi gphi
#           define dphidxi1 dgphidxi1
#           define dphidxi2 dgphidxi2
#           include <phi_and_dphidxik_2d.rc>
#           undef etype
#           undef delta
#           undef phi
#           undef dphidxi1
#           undef dphidxi2
            ! Calculation of x, T1 and T2
            x=0.d0
            T1=0.d0
            T2=0.d0
            do kphi=1,e%n_gnodes
              x=x+gphi(kphi)*e%x(:,kphi)
              T1=T1+dgphidxi1(kphi)*e%x(:,kphi)
              T2=T2+dgphidxi2(kphi)*e%x(:,kphi)
            end do
            ! Normal vector is N = T1 x T2
            N(1)=T1(2)*T2(3)-T1(3)*T2(2)
            N(2)=T1(3)*T2(1)-T1(1)*T2(3)
            N(3)=T1(1)*T2(2)-T1(2)*T2(1)
            ! xi -> x jacobian
            jg=sqrt(dot_product(N,N))
            ! Unit normal
            n=N/jg
            ! Distance vector
            rv=x-x_i
            r=sqrt(dot_product(rv,rv))
            dr1=1.d0/r
            dr2=dr1**2
            dr3=dr2*dr1
            ! dr/dn
            drdn=dot_product(rv,n)*dr1
            ! dr/dni
            drdni=-dot_product(rv,n_i)*dr1
            ! n·n_i
            n_dot_ni=dot_product(n,n_i)
            ! Functional shape functions (primary variables)
#           define etype e%ptype
#           define delta e%ptype_delta
#           define phi pphi
#           include <phi_2d.rc>
#           undef etype
#           undef delta
#           undef phi
            ! Functional shape functions (secondary variables)
#           define etype e%stype
#           define delta e%stype_delta
#           define phi sphi
#           include <phi_2d.rc>
#           undef etype
#           undef delta
#           undef phi
            ! Jacobians * weights
            jw=jg*js*jqt*jt(1)*jt(2)*w(1)*w(2)
            ! Functional shape functions * jacobians * weights
            pphijw=pphi*jw
            sphijw=sphi*jw
            ! Add the integration points
            m=m+dr3*(3.d0*drdn*drdni+n_dot_ni)*pphijw
            l=l+dr2*drdni*sphijw
          end do
        end do
      case default
        call fbem_error_message(error_unit,0,'fbem_bem_stapot3d_hbie_ext_st',0,'n_edges not valid')
    end select
    ! Multiply by constants
    m= m*c_1_4pi
    l=-l*c_1_4pi
    ! Reverse element orientation
    if (reverse) m=-m
  end subroutine fbem_bem_stapot3d_hbie_ext_st

  recursive subroutine fbem_bem_stapot3d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,qsp,ks,ns,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                !! Element
    logical                  :: reverse                          !! Reverse orientation
    real(kind=real64)        :: xi_s(2,fbem_n_vertices(e%gtype)) !! Subdivision of the parent element
    real(kind=real64)        :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)        :: n_i(3)                           !! Unit normal vector of the collocation point
    type(fbem_qs_parameters) :: qsp                              !! Quasi-singular integration parameters
    integer                  :: ks                               !! Current level of subdivisions
    integer                  :: ns                               !! Maximum level of subdivision
    real(kind=real64)        :: m(e%n_pnodes)                    !! m vector
    real(kind=real64)        :: l(e%n_snodes)                    !! l vector
    ! Local
    integer           :: gln_near                             ! Gauss-Legendre integ. points required to integrate only  the quasi-singular integrand
    integer           :: gln                                  ! Gauss-Legendre integ. points required to integrate the whole integrand
    logical           :: subdivide                            ! True if subdivision has to be performed
    real(kind=real64) :: barxi(2)                             ! Nearest element coordinates with respect to collocation point
    real(kind=real64) :: barxip(2)                            ! Nearest element subdivision local coordinates with respect to collocation point
    real(kind=real64) :: rmin                                 ! Minimum distance between collocation point and barxi on the element
    real(kind=real64) :: barr                                 ! Telles jacobian at barxi
    real(kind=real64) :: cl                                   ! Characteristic length
    real(kind=real64) :: d                                    ! Normalized distance between collocation point and element subdivision
    integer           :: method                               ! Method used in nearest point algorithm
    real(kind=real64) :: tmp_xi_s(2,fbem_n_vertices(e%gtype)) ! Subdivision
    real(kind=real64) :: x_s(3,e%n_gnodes)                    ! Coordinates of the element subdivision
    real(kind=real64) :: m_tmp(e%n_pnodes)                    ! m vector (temporary)
    real(kind=real64) :: l_tmp(e%n_snodes)                    ! l vector (temporary)
    ! Initialize
    if (ks.eq.1) then
      m=0.d0
      l=0.d0
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
      cl=fbem_characteristic_length(3,e%gtype,x_s,1.d-6)
      call fbem_nearest_element_point_bem(3,e%gtype,x_s,cl,x_i,barxip,rmin,d,method)
    end if
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,5,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_stapot3d_hbie_ext_adp',ns,'maximum number of subdivisions reached')
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
          call fbem_bem_stapot3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,qsp,ks+1,ns,m,l)
          ! SUBTRI 2
          tmp_xi_s(:,1)=xi_s(:,2)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          call fbem_bem_stapot3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,qsp,ks+1,ns,m,l)
          ! SUBTRI 3
          tmp_xi_s(:,1)=xi_s(:,3)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          call fbem_bem_stapot3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,qsp,ks+1,ns,m,l)
          ! SUBTRI 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          call fbem_bem_stapot3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,qsp,ks+1,ns,m,l)
        ! QUADRILATERALS
        case (4)
          ! SUBQUAD 1
          tmp_xi_s(:,1)=xi_s(:,1)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,3)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          call fbem_bem_stapot3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,qsp,ks+1,ns,m,l)
          ! SUBQUAD 2
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=xi_s(:,2)
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,4)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          call fbem_bem_stapot3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,qsp,ks+1,ns,m,l)
          ! SUBQUAD 3
          tmp_xi_s(:,1)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=xi_s(:,3)
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          call fbem_bem_stapot3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,qsp,ks+1,ns,m,l)
          ! SUBQUAD 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          tmp_xi_s(:,2)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=xi_s(:,4)
          call fbem_bem_stapot3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,qsp,ks+1,ns,m,l)
      end select
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_stapot3d_hbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,gln,m_tmp,l_tmp)
      m=m+m_tmp
      l=l+l_tmp
    end if
  end subroutine fbem_bem_stapot3d_hbie_ext_adp

  subroutine fbem_bem_stapot3d_hbie_int(type_g,type_f1,type_f2,delta_f,x_nodes,reverse,xi_i,m,l)
    implicit none
    ! I/O
    integer           :: type_g                          !! Geometrical interpolation
    integer           :: type_f1                         !! Functional interpolation (primary variables)
    integer           :: type_f2                         !! Functional interpolation (secondary variables)
    real(kind=real64) :: delta_f                         !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64) :: x_nodes(3,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
    logical           :: reverse                         !! Normal vector inversion
    real(kind=real64) :: xi_i(2)                         !! Reference coordinates of the singular point.
    real(kind=real64) :: m(fbem_n_nodes(type_f1))        !! m kernel vector
    real(kind=real64) :: l(fbem_n_nodes(type_f2))        !! l kernel vector
    ! Local
    integer           :: ksubtri                          ! Counter variable for subtriangles loop
    integer           :: nsubtri                          ! Number of subtriangles
    integer           :: subtriangle(8)                   ! Vector that contains subtriangles that need to be integrated
    real(kind=real64) :: theta_subtri(2,8)                ! Angles theta of each subtriangle
    real(kind=real64) :: thetap_subtri(2,8)               ! Angles theta' of each subtriangle
    integer           :: ktheta                           ! Counter variable for theta coordinate loop
    integer           :: krho                             ! Counter variable for rho coordinate loop
    integer           :: kphi                             ! Counter coordinates for shape functions loops
    integer           :: nnodes_g                         ! Number of nodes of geometrical interpolation
    integer           :: nnodes_f1                        ! Number of nodes of functional (primary variables) interpolation
    integer           :: ngp_theta                        ! Number of Gauss points for theta coordinate
    integer           :: ngp_rho                          ! Number of Gauss points for rho coordinate
    real(kind=real64) :: thetai, thetaf, thetapi, thetapf ! Initial and ending angles for subtriangle integration
    real(kind=real64) :: w_angular                        ! Weight of the angular coordinate
    real(kind=real64) :: w_radial                         ! Weight of the radial coordinate
    real(kind=real64) :: theta                            ! Angle coordinate theta
    real(kind=real64) :: thetap                           ! Angle coordinate thetap
    real(kind=real64) :: thetapp                          ! Angle coordinate thetap on [0,1] domain
    real(kind=real64) :: jthetap                          ! thetap->thetapp jacobian
    real(kind=real64) :: rhoij                            ! Maximum rho (radial) value for each edge
    real(kind=real64) :: rho                              ! Radial coordinate rho
    real(kind=real64) :: rhop                             ! Radial coordinate rho on [0,1] domain
    real(kind=real64) :: aux(10)                          ! Auxiliary variable for shape functions resources
    real(kind=real64) :: xi(2)                            ! Reference xi_1,xi_2 coordinates
    real(kind=real64) :: phi_f1(fbem_n_nodes(type_f1))    ! Functional shape functions values (primary variables) at xi_1,xi_2
    real(kind=real64) :: phi_f2(fbem_n_nodes(type_f2))    ! Functional shape functions values (secondary variables) at xi_1,xi_2
    real(kind=real64) :: phi_f1_i(fbem_n_nodes(type_f1))  ! Functional shape functions values (primary variables) at xi_i1,xi_i2
    real(kind=real64) :: psi_i(3,fbem_n_nodes(type_f1))   ! psi vectors at xi_i1,xi_i2
    real(kind=real64) :: phi_g(fbem_n_nodes(type_g))      ! Geometrical shape functions values at xi_1,xi_2
    real(kind=real64) :: dphidxi1_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64) :: dphidxi2_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64) :: jxi1, jxi2                       ! xi1->x, xi2->x jacobians
    real(kind=real64) :: x_i(3)                           ! Collocation point position vector
    real(kind=real64) :: n_i(3)                           ! Unit normal at collocation point
    real(kind=real64) :: x(3)                             ! Position vector at xi_1,xi_2
    real(kind=real64) :: T1(3), T2(3)                     ! Tangent vectors at xi_1,xi_2
    real(kind=real64) :: v1(3), v2(3)                     ! Orthogonal tangent vectors
    real(kind=real64) :: Mvt(2,2)                         ! Orthogonal coordinates transformation matrix
    real(kind=real64) :: dxidv(2,2)                       ! xi coordinates derivatives with respect to v orthogonal coordinates
    real(kind=real64) :: detMvt                           ! Mvt(1,1)*Mvt(2,2)-Mvt(1,2)*Mvt(2,1)
    real(kind=real64) :: n(3)                             ! Unit normal vector at xi_1,xi_2
    real(kind=real64) :: rv(3)                            ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64) :: r, dr1, dr2, dr3               ! Distance vector module and its inverse
    real(kind=real64) :: jg                               ! Geometric jacobian
    real(kind=real64) :: drdn                             ! Partial derivative of r respect to unit normal
    real(kind=real64) :: drdni                            ! Partial derivative of r respect to unit normal on collocation point
    real(kind=real64) :: n_dot_ni                         ! Scalar product between integration point unit normal and collocation point unit normal
    real(kind=real64) :: jw                               ! Jacobians * weights
    real(kind=real64) :: phif1jw(fbem_n_nodes(type_f1))   ! Functional shape functions values (primary variables) at xi_1,xi_2
    real(kind=real64) :: phif2jw(fbem_n_nodes(type_f2))   ! Functional shape functions values (secondary variables) at xi_1,xi_2
    real(kind=real64) :: costheta, sintheta               ! cos(theta), sin(theta)
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
    ! Initialization
    m=0.d0
    l=0.d0
    ! Number of nodes of geometrical interpolation
    nnodes_g=fbem_n_nodes(type_g)
    ! Number of nodes of functional (primary variables) interpolation
    nnodes_f1=fbem_n_nodes(type_f1)
    ! Calculate real coordinates, tangent vectors, and unit normal of the collocation point
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
    ! v1 = t1
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
    if (fbem_check_xi1xi2_edge(type_g,xi_i)) then
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
          dr1=1.d0/r
          dr2=dr1**2
          dr3=dr2*dr1
          ! dr/dn
          drdn=dot_product(rv,n)*dr1
          ! dr/dni
          drdni=-dot_product(rv,n_i)*dr1
          ! n·n_i
          n_dot_ni=dot_product(n,n_i)
          ! Jacobians * weights
          jw=jg*rho*jthetap*w_angular*w_radial
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
          m=m+3.0d0*dr3*drdn*drdni*(phi_f1-phi_f1_i)*jw
          m=m+dr3*n_dot_ni*(phi_f1-phi_f1_i-psi_i(1,:)*rv(1)-psi_i(2,:)*rv(2)-psi_i(3,:)*rv(3))*jw
          m=m-(psi_i(1,:)*n(1)+psi_i(2,:)*n(2)+psi_i(3,:)*n(3))*dr2*drdni*jw
          l=l+dr2*drdni*phif2jw
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
  end subroutine fbem_bem_stapot3d_hbie_int

  !! This subroutine calculates adaptatively the line integral (r x ni)·t/r**3 for the HBIE interior integration using Telles
  !! transformation and subdivision. Before calling this subroutine, mli1 must be set to zero.
  recursive subroutine fbem_bem_stapot3d_hbie_int_li1(gtype,x_nodes,xi_s,x_i,n_i,ngp_min,qsp,ks,ns,mli1)
    implicit none
    ! I/O
    integer                      :: gtype                           !! Geometrial interpolation
    real(kind=real64)            :: x_nodes(3,fbem_n_nodes(gtype))  !! Position vectors of geometrical nodes
    real(kind=real64)            :: xi_s(1,2)                       !! Subdivision of the parent element
    real(kind=real64)            :: x_i(3)                          !! Collocation point position vector
    real(kind=real64)            :: n_i(3)                          !! Collocation point normal vector
    integer                      :: ngp_min                         !! Minimum number of Gauss points
    type(fbem_qs_parameters)     :: qsp                             !! Quasi-singular integration parameters
    integer                      :: ks                              !! Current level of subdivisions
    integer                      :: ns                              !! Maximum level of subdivision
    real(kind=real64)            :: mli1                            !! Line integral (r x ni)·t/r**3
    ! Local
    integer                      :: gln_near                        ! 1D Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer                      :: gln                             ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical                      :: subdivide                       ! True if subdivision has to be performed
    real(kind=real64)            :: barxip(1)                       ! Nearest element subdivision local coordinate with respect to collocation point
    real(kind=real64)            :: rmin                            ! Minimum distance between collocation point and barxi on the element
    real(kind=real64)            :: barr                            ! Telles jacobian at barxi
    real(kind=real64)            :: cl                              ! Characteristic length
    real(kind=real64)            :: d                               ! Normalized distance between collocation point and element subdivision
    integer                      :: method                          ! Method used in nearest point algorithm
    real(kind=real64)            :: tmp_xi_s(1,2)                   ! Subdivision
    real(kind=real64)            :: x_s(3,fbem_n_nodes(gtype))      ! Coordinates of the element subdivision
    real(kind=real64)            :: mli1_tmp                        ! Line integral (r x ni)·t/r**3 (temporary)
    ! Local calculation
    integer                      :: kphi                            ! Counter variable for shape functions loops
    integer                      :: n_gnodes                        ! Number of nodes of the geometrical element
    integer                      :: kip                             ! Counter variable for integration points
    real(kind=real64)            :: gamma                           ! Reference coordinate
    type(fbem_telles_parameters) :: telles_parameters               ! Telles parameters
    real(kind=real64)            :: jt                              ! Telles jacobian for each coordinate: xi->gamma
    real(kind=real64)            :: xip                             ! Reference coordinate of the subdivision
    real(kind=real64)            :: js                              ! Jacobian of xip->xi transformation
    real(kind=real64)            :: xi                              ! Reference coordinate
    real(kind=real64)            :: w                               ! Weights of each integration point
    real(kind=real64)            :: aux(10)                         ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(fbem_n_nodes(gtype))       ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi(fbem_n_nodes(gtype))   ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: x(3), T(3), rv(3), r, dr3       ! Position, tangent, distance vector, distance and 1/r^3 at xi
    real(kind=real64)            :: rxni(3)                         ! r x ni
    real(kind=real64)            :: jg                              ! Geometric jacobian
    real(kind=real64)            :: jw                              ! Jacobians * weight
    ! Initialize
    if (ks.eq.1) then
      xi_s(1,1)=-1.d0
      xi_s(1,2)= 1.d0
      x_s=x_nodes
    else
      call fbem_obtain_element_subdivision_coordinates(3,gtype,x_nodes,xi_s,x_s)
    end if
    cl=fbem_characteristic_length(3,gtype,x_s,1.d-12)
    call fbem_nearest_element_point_bem(3,gtype,x_s,cl,x_i,barxip,rmin,d,method)
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(3,gtype,2,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_stapot3d_hbie_int_li1',ns,'maximum number of subdivisions reached')
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
      call fbem_bem_stapot3d_hbie_int_li1(gtype,x_nodes,tmp_xi_s,x_i,n_i,ngp_min,qsp,ks+1,ns,mli1)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_stapot3d_hbie_int_li1(gtype,x_nodes,tmp_xi_s,x_i,n_i,ngp_min,qsp,ks+1,ns,mli1)
    ! Calculate the subdivided element using Telles transformation
    else
      ! Initialization
      n_gnodes=fbem_n_nodes(gtype)
      mli1_tmp=0.d0
      gln=max(gln_near,ngp_min)
      barr=fbem_telles_barr(d,fbem_f_any)
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
#       define etype gtype
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
        do kphi=1,n_gnodes
          x=x+gphi(kphi)*x_nodes(:,kphi)
          T=T+dgphidxi(kphi)*x_nodes(:,kphi)
        end do
        ! Geometric jacobian
        jg=sqrt(dot_product(T,T))
        ! Unit tangent
        t=T/jg
        ! Distance vector
        rv=x-x_i
        ! Distance vector norm
        r=sqrt(dot_product(rv,rv))
        dr3=1.0d0/r**3
        ! r x ni
        rxni(1)=rv(2)*n_i(3)-rv(3)*n_i(2)
        rxni(2)=rv(3)*n_i(1)-rv(1)*n_i(3)
        rxni(3)=rv(1)*n_i(2)-rv(2)*n_i(1)
        ! Jacobians * weight
        jw=jg*js*jt*w
        ! Add
        mli1_tmp=mli1_tmp+dr3*dot_product(rxni,t)*jw
      end do
      mli1=mli1+mli1_tmp
    end if
  end subroutine fbem_bem_stapot3d_hbie_int_li1

  !! This subroutine calculates adaptatively the line integral (e_k x ni)·t/r for the HBIE interior integration using Telles
  !! transformation and subdivision. Before calling this subroutine, mli1 must be set to zero.
  recursive subroutine fbem_bem_stapot3d_hbie_int_li2(gtype,x_nodes,xi_s,x_i,n_i,ngp_min,qsp,ks,ns,mli2)
    implicit none
    ! I/O
    integer                      :: gtype                           !! Geometrial interpolation
    real(kind=real64)            :: x_nodes(3,fbem_n_nodes(gtype))  !! Position vectors of geometrical nodes
    real(kind=real64)            :: xi_s(1,2)                       !! Subdivision of the parent element
    real(kind=real64)            :: x_i(3)                          !! Collocation point position vector
    real(kind=real64)            :: n_i(3)                          !! Collocation point normal vector
    integer                      :: ngp_min                         !! Minimum number of Gauss points
    type(fbem_qs_parameters)     :: qsp                             !! Quasi-singular integration parameters
    integer                      :: ks                              !! Current level of subdivisions
    integer                      :: ns                              !! Maximum level of subdivision
    real(kind=real64)            :: mli2(3)                         !! Line integral (e_k x ni)·t/r
    ! Local
    integer                      :: gln_near                        ! 1D Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer                      :: gln                             ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical                      :: subdivide                       ! True if subdivision has to be performed
    real(kind=real64)            :: barxip(1)                       ! Nearest element subdivision local coordinate with respect to collocation point
    real(kind=real64)            :: rmin                            ! Minimum distance between collocation point and barxi on the element
    real(kind=real64)            :: barr                            ! Telles jacobian at barxi
    real(kind=real64)            :: cl                              ! Characteristic length
    real(kind=real64)            :: d                               ! Normalized distance between collocation point and element subdivision
    integer                      :: method                          ! Method used in nearest point algorithm
    real(kind=real64)            :: tmp_xi_s(1,2)                   ! Subdivision
    real(kind=real64)            :: x_s(3,fbem_n_nodes(gtype))      ! Coordinates of the element subdivision
    real(kind=real64)            :: mli2_tmp(3)                     ! Line integral (e_k x ni)·t/r (temporary)
    ! Local calculation
    integer                      :: kphi                            ! Counter variable for shape functions loops
    integer                      :: n_gnodes                        ! Number of nodes of the geometrical element
    integer                      :: kip                             ! Counter variable for integration points
    integer                      :: k                               ! Coordinate counter
    real(kind=real64)            :: gamma                           ! Reference coordinate
    type(fbem_telles_parameters) :: telles_parameters               ! Telles parameters
    real(kind=real64)            :: jt                              ! Telles jacobian for each coordinate: xi->gamma
    real(kind=real64)            :: xip                             ! Reference coordinate of the subdivision
    real(kind=real64)            :: js                              ! Jacobian of xip->xi transformation
    real(kind=real64)            :: xi                              ! Reference coordinate
    real(kind=real64)            :: w                               ! Weights of each integration point
    real(kind=real64)            :: aux(10)                         ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(fbem_n_nodes(gtype))       ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi(fbem_n_nodes(gtype))   ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: x(3), T(3), rv(3), r, dr1       ! Position, tangent, distance vector, distance and 1/r at xi
    real(kind=real64)            :: e_k(3)                          ! e_k
    real(kind=real64)            :: ekni(3)                         ! e_k x ni
    real(kind=real64)            :: jg                              ! Geometric jacobian
    real(kind=real64)            :: jw                              ! Jacobians * weight
    ! Initialize
    if (ks.eq.1) then
      xi_s(1,1)=-1.d0
      xi_s(1,2)= 1.d0
      x_s=x_nodes
    else
      call fbem_obtain_element_subdivision_coordinates(3,gtype,x_nodes,xi_s,x_s)
    end if
    cl=fbem_characteristic_length(3,gtype,x_s,1.d-12)
    call fbem_nearest_element_point_bem(3,gtype,x_s,cl,x_i,barxip,rmin,d,method)
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(3,gtype,1,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_stapot3d_hbie_int_li2',ns,'maximum number of subdivisions reached')
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
      call fbem_bem_stapot3d_hbie_int_li2(gtype,x_nodes,tmp_xi_s,x_i,n_i,ngp_min,qsp,ks+1,ns,mli2)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_stapot3d_hbie_int_li2(gtype,x_nodes,tmp_xi_s,x_i,n_i,ngp_min,qsp,ks+1,ns,mli2)
    ! Calculate the subdivided element using Telles transformation
    else
      ! Initialization
      n_gnodes=fbem_n_nodes(gtype)
      mli2_tmp=0.d0
      gln=max(gln_near,ngp_min)
      barr=fbem_telles_barr(d,fbem_f_any)
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
#       define etype gtype
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
        do kphi=1,n_gnodes
          x=x+gphi(kphi)*x_nodes(:,kphi)
          T=T+dgphidxi(kphi)*x_nodes(:,kphi)
        end do
        ! Geometric jacobian
        jg=sqrt(dot_product(T,T))
        ! Unit tangent
        t=T/jg
        ! Distance vector
        rv=x-x_i
        ! Distance vector norm
        r=sqrt(dot_product(rv,rv))
        dr1=1.d0/r
        ! Jacobians * weight
        jw=jg*js*jt*w
        ! Loop through coordinates
        do k=1,3
          ! e_k
          e_k=0.d0
          e_k(k)=1.d0
          ! e_k x ni
          ekni(1)=e_k(2)*n_i(3)-e_k(3)*n_i(2)
          ekni(2)=e_k(3)*n_i(1)-e_k(1)*n_i(3)
          ekni(3)=e_k(1)*n_i(2)-e_k(2)*n_i(1)
          ! Add kernels
          mli2_tmp(k)=mli2_tmp(k)+dr1*dot_product(ekni,t)*jw
        end do
      end do
      mli2=mli2+mli2_tmp
    end if
  end subroutine fbem_bem_stapot3d_hbie_int_li2

  !! This subroutine calculates adaptatively the line integrals (r x ni)·t/r**3 and (e_k x ni)·t/r for the HBIE interior integration
  !! using Telles transformation and subdivision. Before calling this subroutine, mli1 must be set to zero.
  recursive subroutine fbem_bem_stapot3d_hbie_int_li(gtype,x_nodes,xi_s,x_i,n_i,ngp_min,qsp,ks,ns,mli1,mli2)
    implicit none
    ! I/O
    integer                      :: gtype                           !! Geometrial interpolation
    real(kind=real64)            :: x_nodes(3,fbem_n_nodes(gtype))  !! Position vectors of geometrical nodes
    real(kind=real64)            :: xi_s(1,2)                       !! Subdivision of the parent element
    real(kind=real64)            :: x_i(3)                          !! Collocation point position vector
    real(kind=real64)            :: n_i(3)                          !! Collocation point normal vector
    integer                      :: ngp_min                         !! Minimum number of Gauss points
    type(fbem_qs_parameters)     :: qsp                             !! Quasi-singular integration parameters
    integer                      :: ks                              !! Current level of subdivisions
    integer                      :: ns                              !! Maximum level of subdivision
    real(kind=real64)            :: mli1                            !! Line integral (r x ni)·t/r**3
    real(kind=real64)            :: mli2(3)                         !! Line integral (e_k x ni)·t/r
    ! Local
    integer                      :: gln_near                        ! 1D Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer                      :: gln                             ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical                      :: subdivide                       ! True if subdivision has to be performed
    real(kind=real64)            :: barxip(1)                       ! Nearest element subdivision local coordinate with respect to collocation point
    real(kind=real64)            :: rmin                            ! Minimum distance between collocation point and barxi on the element
    real(kind=real64)            :: barr                            ! Telles jacobian at barxi
    real(kind=real64)            :: cl                              ! Characteristic length
    real(kind=real64)            :: d                               ! Normalized distance between collocation point and element subdivision
    integer                      :: method                          ! Method used in nearest point algorithm
    real(kind=real64)            :: tmp_xi_s(1,2)                   ! Subdivision
    real(kind=real64)            :: x_s(3,fbem_n_nodes(gtype))      ! Coordinates of the element subdivision
    real(kind=real64)            :: mli1_tmp                        ! Line integral (r x ni)·t/r**3 (temporary)
    real(kind=real64)            :: mli2_tmp(3)                     ! Line integral (e_k x ni)·t/r (temporary)
    ! Local calculation
    integer                      :: kphi                            ! Counter variable for shape functions loops
    integer                      :: n_gnodes                        ! Number of nodes of the geometrical element
    integer                      :: kip                             ! Counter variable for integration points
    integer                      :: k                               ! Coordinate counter
    real(kind=real64)            :: gamma                           ! Reference coordinate
    type(fbem_telles_parameters) :: telles_parameters               ! Telles parameters
    real(kind=real64)            :: jt                              ! Telles jacobian for each coordinate: xi->gamma
    real(kind=real64)            :: xip                             ! Reference coordinate of the subdivision
    real(kind=real64)            :: js                              ! Jacobian of xip->xi transformation
    real(kind=real64)            :: xi                              ! Reference coordinate
    real(kind=real64)            :: w                               ! Weights of each integration point
    real(kind=real64)            :: aux(10)                         ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(fbem_n_nodes(gtype))       ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi(fbem_n_nodes(gtype))   ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: x(3), T(3), rv(3), r, dr1, dr3  ! Position, tangent, distance vector, distance and 1/r^n at xi
    real(kind=real64)            :: rxni(3)                         ! r x ni
    real(kind=real64)            :: e_k(3)                          ! e_k
    real(kind=real64)            :: ekni(3)                         ! e_k x ni
    real(kind=real64)            :: jg                              ! Geometric jacobian
    real(kind=real64)            :: jw                              ! Jacobians * weight
    ! Initialize
    if (ks.eq.1) then
      xi_s(1,1)=-1.d0
      xi_s(1,2)= 1.d0
      x_s=x_nodes
    else
      call fbem_obtain_element_subdivision_coordinates(3,gtype,x_nodes,xi_s,x_s)
    end if
    cl=fbem_characteristic_length(3,gtype,x_s,1.d-12)
    call fbem_nearest_element_point_bem(3,gtype,x_s,cl,x_i,barxip,rmin,d,method)
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(3,gtype,2,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_stapot3d_hbie_int_li',ns,'maximum number of subdivisions reached')
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
      call fbem_bem_stapot3d_hbie_int_li(gtype,x_nodes,tmp_xi_s,x_i,n_i,ngp_min,qsp,ks+1,ns,mli1,mli2)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_stapot3d_hbie_int_li(gtype,x_nodes,tmp_xi_s,x_i,n_i,ngp_min,qsp,ks+1,ns,mli1,mli2)
    ! Calculate the subdivided element using Telles transformation
    else
      ! Initialization
      n_gnodes=fbem_n_nodes(gtype)
      mli1_tmp=0.d0
      mli2_tmp=0.d0
      gln=max(gln_near,ngp_min)
      barr=fbem_telles_barr(d,fbem_f_any)
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
#       define etype gtype
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
        do kphi=1,n_gnodes
          x=x+gphi(kphi)*x_nodes(:,kphi)
          T=T+dgphidxi(kphi)*x_nodes(:,kphi)
        end do
        ! Geometric jacobian
        jg=sqrt(dot_product(T,T))
        ! Jacobians * weight
        jw=jg*js*jt*w
        ! Unit tangent
        t=T/jg
        ! Distance vector
        rv=x-x_i
        ! Distance vector norm
        r=sqrt(dot_product(rv,rv))
        dr1=1.d0/r
        dr3=1.0d0/r**3
        ! r x ni
        rxni(1)=rv(2)*n_i(3)-rv(3)*n_i(2)
        rxni(2)=rv(3)*n_i(1)-rv(1)*n_i(3)
        rxni(3)=rv(1)*n_i(2)-rv(2)*n_i(1)
        ! Add
        mli1_tmp=mli1_tmp+dr3*dot_product(rxni,t)*jw
        ! Loop through coordinates
        do k=1,3
          ! e_k
          e_k=0.d0
          e_k(k)=1.d0
          ! e_k x ni
          ekni(1)=e_k(2)*n_i(3)-e_k(3)*n_i(2)
          ekni(2)=e_k(3)*n_i(1)-e_k(1)*n_i(3)
          ekni(3)=e_k(1)*n_i(2)-e_k(2)*n_i(1)
          ! Add
          mli2_tmp(k)=mli2_tmp(k)+dr1*dot_product(ekni,t)*jw
        end do
      end do
      mli1=mli1+mli1_tmp
      mli2=mli2+mli2_tmp
    end if
  end subroutine fbem_bem_stapot3d_hbie_int_li

  subroutine fbem_bem_stapot3d_hbie_auto(e,reverse,x_i,n_i,qsp,ns,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e             !! Element
    logical                  :: reverse       !! Reverse orientation
    real(kind=real64)        :: x_i(3)        !! Collocation point
    real(kind=real64)        :: n_i(3)        !! Unit normal vector of the collocation point
    type(fbem_qs_parameters) :: qsp           !! Quasi-singular integration parameters
    integer                  :: ns            !! Maximum level of subdivisions
    real(kind=real64)        :: m(e%n_pnodes) !! m vector
    real(kind=real64)        :: l(e%n_snodes) !! l vector
    ! Local
    real(kind=real64) :: r(3)                             ! Distance vector
    real(kind=real64) :: rmin                             ! Minimum distance between element and x_i
    real(kind=real64) :: barxi(2)                         ! Local coordinates of the nearest element point with respect to x_i
    real(kind=real64) :: d                                ! Dimensionless distance
    integer           :: delta                            ! Control variable
    real(kind=real64) :: xi_s(2,fbem_n_vertices(e%gtype)) ! Local coordinates of the element subdivision
    integer           :: method                           ! Method used when calculating the nearest element point
    integer           :: gln_near                         ! Gauss-Legendre integration points required by the quasi-singular function
    integer           :: gln                              ! Gauss-Legendre integration points used in the integration
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
        call fbem_bem_stapot3d_hbie_int(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi,m,l)
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
          call fbem_bem_stapot3d_hbie_ext_pre(ps,e,reverse,x_i,n_i,m,l)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_stapot3d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,qsp,1,ns,m,l)
        end if
    end select
  end subroutine fbem_bem_stapot3d_hbie_auto

  ! ================================================================================================================================

end module fbem_bem_stapot3d
