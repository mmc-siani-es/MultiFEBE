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
!! <b> This module implements the calculation of element-wise integrals of Boundary Integral Equations of 3D elastodynamics.</b>
module fbem_bem_harela3d

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
  use fbem_bem_staela3d

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all are private
  private

  ! ================================================================================================================================
  ! INITIAL SETUP
  public :: fbem_bem_harela3d_parameters
  public :: fbem_bem_harela3d_calculate_parameters
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE)
  ! Free-term
  public :: fbem_bem_harela3d_sbie_freeterm
  ! BOUNDARY ELEMENTS
  ! Exterior integration
  public :: fbem_bem_harela3d_sbie_ext_pre
  public :: fbem_bem_harela3d_sbie_ext_st
  public :: fbem_bem_harela3d_sbie_ext_adp
  ! Interior integration
  public :: fbem_bem_harela3d_sbie_int
  ! Automatic integration
  public :: fbem_bem_harela3d_sbie_auto
  ! BODY LOADS
  ! Exterior integration
  public :: fbem_bem_harela3d_sbie_bl_ext_pre
  public :: fbem_bem_harela3d_sbie_bl_ext_st
  public :: fbem_bem_harela3d_sbie_bl_ext_adp
  ! Interior integration
  public :: fbem_bem_harela3d_sbie_bl_int
  ! Automatic integration
  public :: fbem_bem_harela3d_sbie_bl_auto
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (HBIE)
  ! BOUNDARY ELEMENTS
  ! Exterior integration
  public :: fbem_bem_harela3d_hbie_ext_pre
  public :: fbem_bem_harela3d_hbie_ext_st
  public :: fbem_bem_harela3d_hbie_ext_adp
  ! Interior integration
  public :: fbem_bem_harela3d_hbie_int
  ! Automatic integration
  public :: fbem_bem_harela3d_hbie_auto
  ! BODY LOADS
  ! Exterior integration
  public :: fbem_bem_harela3d_hbie_bl_ext_pre
  public :: fbem_bem_harela3d_hbie_bl_ext_st
  public :: fbem_bem_harela3d_hbie_bl_ext_adp
  ! Interior integration
  public :: fbem_bem_harela3d_hbie_bl_int
  ! Automatic integration
  public :: fbem_bem_harela3d_hbie_bl_auto
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! BOTH SBIE AND HBIE
  ! Exterior integration
  public :: fbem_bem_harela3d_shbie_ext_pre
  ! Automatic integration
  public :: fbem_bem_harela3d_shbie_auto
  ! ================================================================================================================================

  ! ================================================================================================================================
  !! Data structure for region parameters (material properties and frequency dependant parameters)
  type fbem_bem_harela3d_parameters
    ! Region parameters
    complex(kind=real64) :: lambda
    complex(kind=real64) :: mu
    real(kind=real64)    :: rho
    ! Frequency
    real(kind=real64)    :: omega
    ! Wave propagation speeds
    complex(kind=real64) :: c1, c2
    ! Wavenumbers
    complex(kind=real64) :: k1, k2
    ! Coefficients of the components of the fundamental solutions u*, t*, d* and s*
    complex(kind=real64) :: psi(6), chi(6)
    complex(kind=real64) :: T1(10), T2(9), T3(9)
    complex(kind=real64) :: S1(11), S2(11), S3(12), S4(10), S5(11)
    complex(kind=real64) :: cte_u, cte_t, cte_d, cte_s
  end type fbem_bem_harela3d_parameters
  ! ================================================================================================================================


contains

  ! ================================================================================================================================
  ! INITIAL SETUP

  !! Subroutine that calculate all parameters for given properties and frequency
  subroutine fbem_bem_harela3d_calculate_parameters(lambda,mu,rho,omega,p)
    implicit none
    ! I/O
    complex(kind=real64)               :: lambda
    complex(kind=real64)               :: mu
    real(kind=real64)                  :: rho
    real(kind=real64)                  :: omega
    type(fbem_bem_harela3d_parameters) :: p
    ! Local
    complex(kind=real64)               :: c1, c2, k1, k2
    !
    ! Local parameters
    !
    c1=sqrt((lambda+2.d0*mu)/rho)
    c2=sqrt(mu/rho)
    k1=omega/c1
    k2=omega/c2
    !
    ! Save to structure the main parameters
    !
    p%omega=omega
    p%lambda=lambda
    p%mu=mu
    p%rho=rho
    p%c1=c1
    p%c2=c2
    p%k1=k1
    p%k2=k2
    !
    ! Coefficients of the components of the fundamental solutions u*, t*, d* and s*
    !
    ! psi
    p%psi(1)=0.5d0*(1.d0+c2**2/c1**2)
    p%psi(2)=-1.d0/3.d0*(2.d0/c2+c2**2/c1**3)*c_im*omega
    p%psi(3)=c_im*k1/k2**2
    p%psi(4)=1.d0/(c_im*k2)
    p%psi(5)=1.d0/k2**2
    p%psi(6)=-1.d0/k2**2
    ! chi
    p%chi(1)=-0.5d0*(1.d0-c2**2/c1**2)
    p%chi(2)=-c2**2/c1**2
    p%chi(3)=-3.d0*c2**2/c1**2/(c_im*k1)
    p%chi(4)=3.d0/(c_im*k2)
    p%chi(5)=-3.d0*c2**2/c1**2/(c_im*k1)**2
    p%chi(6)=3.d0/(c_im*k2)**2
    ! T1
    p%T1( 1)=3.d0*(c2**2/c1**2-1.d0)
    p%T1( 2)=-0.25d0*(1.d0/c2**2-c2**2/c1**4)*omega**2
    p%T1( 3)=-2.d0*c_im*k1*c2**2/c1**2
    p%T1( 4)=2.d0*c_im*k2
    p%T1( 5)=-12.d0*c2**2/c1**2
    p%T1( 6)=12.d0
    p%T1( 7)=-c2**2/c1**2*30.d0/(c_im*k1)
    p%T1( 8)=30.d0/(c_im*k2)
    p%T1( 9)=-c2**2/c1**2*30.d0/(c_im*k1)**2
    p%T1(10)=30.d0/(c_im*k2)**2
    ! T2
    p%T2( 1)=-c2**2/c1**2
    p%T2( 2)=-0.25d0*(1.d0/c2**2+c2**2/c1**4)*omega**2
    p%T2( 3)=-c_im*k2
    p%T2( 4)=2.d0*c2**2/c1**2
    p%T2( 5)=-3.d0
    p%T2( 6)=c2**2/c1**2*6.d0/(c_im*k1)
    p%T2( 7)=-6.d0/(c_im*k2)
    p%T2( 8)=c2**2/c1**2*6.d0/(c_im*k1)**2
    p%T2( 9)=-6.d0/(c_im*k2)**2
    ! T3
    p%T3( 1)=c2**2/c1**2
    p%T3( 2)=0.25d0*(3.d0*c2**2/c1**4-2.d0/c1**2+1.d0/c2**2)*omega**2
    p%T3( 3)=(2.d0*c2**2/c1**2-1.d0)*c_im*k1
    p%T3( 4)=4.d0*c2**2/c1**2-1.d0
    p%T3( 5)=-2.d0
    p%T3( 6)=c2**2/c1**2*6.d0/(c_im*k1)
    p%T3( 7)=-6.d0/(c_im*k2)
    p%T3( 8)=c2**2/c1**2*6.d0/(c_im*k1)**2
    p%T3( 9)=-6.d0/(c_im*k2)**2
    ! S1
    p%S1( 1)=3.d0*(1.d0-2.d0*c2**2/c1**2)
    p%S1( 2)=-0.5d0*c2**2/c1**4*omega**2
    p%S1( 3)=k2**2
    p%S1( 4)=4.d0*c_im*k1*k1**2/k2**2
    p%S1( 5)=-7.d0*c_im*k2
    p%S1( 6)=24.d0*k1**2/k2**2
    p%S1( 7)=-27.d0
    p%S1( 8)=60.d0*k1**2/k2**2/(c_im*k1)
    p%S1( 9)=-60.d0/(c_im*k2)
    p%S1(10)=60.d0*k1**2/k2**2/(c_im*k1)**2
    p%S1(11)=-60.d0/(c_im*k2)**2
    ! S2
    p%S2( 1)=6.d0*c2**2/c1**2
    p%S2( 2)=(0.5d0/c2**2+1.5d0*c2**2/c1**4-1.d0/c1**2)*omega**2
    p%S2( 3)=2.d0*c2**2/c1**4*(c1**2/c2**2-2.d0)*omega**2
    p%S2( 4)=2.d0*c_im*c2**2/c1**3*(8.d0-3.d0*c1**2/c2**2)*omega
    p%S2( 5)=-4.d0*c_im*k2
    p%S2( 6)=6.d0*c2**2/c1**2*(6.d0-c1**2/c2**2)
    p%S2( 7)=-24.d0
    p%S2( 8)=c2**2/c1**2*60.d0/(c_im*k1)
    p%S2( 9)=-60d0/(c_im*k2)
    p%S2(10)=60.d0/(c_im*k2)**2
    p%S2(11)=-60.d0/(c_im*k2)**2
    ! S3
    p%S3( 1)=30.d0*(1.d0-c2**2/c1**2)
    p%S3( 2)=1.5d0*(1.d0/c2**2-c2**2/c1**4)*omega**2
    p%S3( 3)=-4.d0*c2**2/c1**2*k1**2
    p%S3( 4)=4.d0*k2**2
    p%S3( 5)=40.d0*c2**2/c1**2*c_im*k1
    p%S3( 6)=-40.d0*c_im*k2
    p%S3( 7)=180.d0*c2**2/c1**2
    p%S3( 8)=-180.d0
    p%S3( 9)=420.d0*c2**2/c1**2/(c_im*k1)
    p%S3(10)=-420.d0/(c_im*k2)
    p%S3(11)=420.d0/(c_im*k2)**2
    p%S3(12)=-420.d0/(c_im*k2)**2
    ! S4
    p%S4( 1)=2.d0*c2**2/c1**2
    p%S4( 2)=0.5d0*(c2**2/c1**4+1.d0/c2**2)*omega**2
    p%S4( 3)=-2.d0/5.d0*(1.d0/c2**3+2.d0/3.d0*c2**2/c1**5)*c_im*omega**3
    p%S4( 4)=2.d0*c_im*k2
    p%S4( 5)=-4.d0*c2**2/c1**2
    p%S4( 6)=6.d0
    p%S4( 7)=-12.d0*c2**2/c1**2/(c_im*k1)
    p%S4( 8)=12.d0/(c_im*k2)
    p%S4( 9)=-12.d0/(c_im*k2)**2
    p%S4(10)=12.d0/(c_im*k2)**2
    ! S5
    p%S5( 1)=2.d0*(1.d0-3.d0*c2**2/c1**2)
    p%S5( 2)=(-2.d0/c1**2+0.5d0/c2**2+0.5d0*c2**2/c1**4)*omega**2
    p%S5( 3)=(8.d0/3.d0/c1**3+4.d0/15.d0/c2**3-1.d0/c1/c2**2-24.d0/15.d0*c2**2/c1**5)*c_im*omega**3
    p%S5( 4)=(-4.d0/c1**2+1.d0/c2**2+4.d0*c2**2/c1**4)*omega**2
    p%S5( 5)=4.d0*c_im*(1.d0/c1-2.d0*c2**2/c1**3)*omega
    p%S5( 6)=4.d0*(1.d0-3.d0*c2**2/c1**2)
    p%S5( 7)=4.d0
    p%S5( 8)=12.d0*c_im*k1/k2**2
    p%S5( 9)=-12.d0*c_im/k2
    p%S5(10)=12.d0/k2**2
    p%S5(11)=-12.d0/k2**2
    ! Constant of u*: 1/(4*pi*mu)
    p%cte_u=c_1_4pi/mu
    ! Constant of t*: 1/(4*pi)
    p%cte_t=c_1_4pi
    ! Constant of d*: 1/(4*pi)
    p%cte_d=c_1_4pi
    ! Constant of s*: mu/(4*pi)
    p%cte_s=c_1_4pi*mu
  end subroutine fbem_bem_harela3d_calculate_parameters
  ! ================================================================================================================================

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
  subroutine fbem_bem_harela3d_sbie_freeterm(n_elements,n,t,tol,nu,c)
    implicit none
    ! I/O
    integer               :: n_elements            !! Number of elements
    real(kind=real64)     :: n(3,n_elements)       !! Unit normal of each element
    real(kind=real64)     :: t(3,n_elements)       !! Unit element boundary tangent of each element
    real(kind=real64)     :: tol                   !! Geometric tolerance [1.0e-12,1.0e-3] (default 1.0e-6)
    complex(kind=real64)  :: nu                    !! Poisson's ratio
    complex(kind=real64)  :: c(3,3)                !! Free-term matrix
    ! Local
    integer            :: kc, kc1, kc2, ki, kj  ! Counters
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
    real(kind=real64)  :: nxn(3), nxndr, ndn, sum_a, cp
    real(kind=real64)  :: rmr(3), rmrxn(3), vtv(3,3), sum_b(3,3)
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
    ! Potential part
    sum_a=0.0d0
    do ki=1,n_elements
      ! nxn = n_{i-1} x n_{i}
      nxn(1)=ln(2,ki-1)*ln(3,ki)-ln(3,ki-1)*ln(2,ki)
      nxn(2)=ln(3,ki-1)*ln(1,ki)-ln(1,ki-1)*ln(3,ki)
      nxn(3)=ln(1,ki-1)*ln(2,ki)-ln(2,ki-1)*ln(1,ki)
      ! nxndr = nxn · t_{i}
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
    cp=1.0d0/(4.0d0*c_pi)*(2.0d0*c_pi+sum_a)
    ! Additional part
    do kc1=1,3
      do kc2=1,3
        sum_b(kc1,kc2)=0.0d0
      end do
    end do
    do ki=1,n_elements
      ! rmr = t_{i+1}-t_{i}
      do kc=1,3
        rmr(kc)=lt(kc,ki+1)-lt(kc,ki)
      end do
      ! rmrxn = rmr x n_{i}
      rmrxn(1)=rmr(2)*ln(3,ki)-rmr(3)*ln(2,ki)
      rmrxn(2)=rmr(3)*ln(1,ki)-rmr(1)*ln(3,ki)
      rmrxn(3)=rmr(1)*ln(2,ki)-rmr(2)*ln(1,ki)
      ! Tensor product of rmrxn and n_{i}
      do kc1=1,3
        do kc2=1,3
          vtv(kc1,kc2)=rmrxn(kc1)*ln(kc2,ki)
        end do
      end do
      ! vtv sum
      do kc1=1,3
        do kc2=1,3
          sum_b(kc1,kc2)=sum_b(kc1,kc2)+vtv(kc1,kc2)
        end do
      end do
    end do
    do kc1=1,3
      do kc2=1,3
        c(kc1,kc2)=-1.0d0/(8.0d0*c_pi*(1.0d0-nu))*sum_b(kc1,kc2)
        if (kc1.eq.kc2) c(kc1,kc2)=c(kc1,kc2)+cp
      end do
    end do
  end subroutine fbem_bem_harela3d_sbie_freeterm

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! BOUNDARY ELEMENTS

  subroutine fbem_bem_harela3d_sbie_ext_pre(ps,e,reverse,x_i,p,h,g)
    implicit none
    ! I/O
    integer                            :: ps                !! Selected precalculated dataset
    type(fbem_bem_element)             :: e                 !! Element
    logical                            :: reverse           !! Reverse normal vector
    real(kind=real64)                  :: x_i(3)            !! Collocation point position vector
    type(fbem_bem_harela3d_parameters) :: p                 !! Parameters of the region
    complex(kind=real64)               :: h(e%n_pnodes,3,3) !! h integration kernels matrix
    complex(kind=real64)               :: g(e%n_snodes,3,3) !! g integration kernels matrix
    ! Local
    integer              :: kip                        ! Counter variable for integration points loop
    real(kind=real64)    :: x(3)                       ! Position vector at integration point
    real(kind=real64)    :: n(3)                       ! Unit normal vector at integration point
    real(kind=real64)    :: pphijw(e%n_pnodes)         ! phi^p * jacobian * weight at integration point
    real(kind=real64)    :: sphijw(e%n_snodes)         ! phi^s * jacobian * weight at integration point
    real(kind=real64)    :: rv(3)                      ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, d1r1, d1r2, d1r3, d1r4  ! Distance vector module and its inverse
    real(kind=real64)    :: drdx(3)                    ! Distance vector derivatives with respect to x_k
    real(kind=real64)    :: drdn                       ! Partial derivative of r respect to unit normal
    integer              :: il, ik                     ! Counter for load / observation components
    complex(kind=real64) :: z(2)                       ! Arguments z=-ikr
    complex(kind=real64) :: EnR(0:6,2)                 ! Exponential function decomposition for each wavenumber
    complex(kind=real64) :: psi, chi                   ! Components of the fundamental solutions
    complex(kind=real64) :: TT1, TT2, TT3              ! Components of the fundamental solutions
    complex(kind=real64) :: fs_u, fs_t                 ! Fundamental solutions values
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
      d1r3=d1r2*d1r1
      d1r4=d1r3*d1r1
      drdx=rv*d1r1
      drdn=dot_product(drdx,n)
      z(1)=-c_im*p%k1*r
      z(2)=-c_im*p%k2*r
      call fbem_zexp_decomposed(2,z,EnR)
      EnR(2,:)=EnR(2,:)*d1r1
      EnR(3,:)=EnR(3,:)*d1r2
      EnR(4,:)=EnR(4,:)*d1r3
      EnR(5,:)=EnR(5,:)*d1r4
      psi=p%psi(1)*d1r1+p%psi(2)+EnR(2,2)+p%psi(3)*EnR(3,1)+p%psi(4)*EnR(3,2)+p%psi(5)*EnR(4,1)+p%psi(6)*EnR(4,2)
      chi=p%chi(1)*d1r1+p%chi(2)*EnR(2,1)+EnR(2,2)+p%chi(3)*EnR(3,1)+p%chi(4)*EnR(3,2)+p%chi(5)*EnR(4,1)+p%chi(6)*EnR(4,2)
      TT1=p%T1(1)*d1r2+p%T1(2)+p%T1(3)*EnR(2,1)+p%T1(4)*EnR(2,2)+p%T1(5)*EnR(3,1)+p%T1(6)*EnR(3,2)+p%T1(7)*EnR(4,1)&
         +p%T1(8)*EnR(4,2)+p%T1(9)*EnR(5,1)+p%T1(10)*EnR(5,2)
      TT2=p%T2(1)*d1r2+p%T2(2)+p%T2(3)*EnR(2,2)+p%T2(4)*EnR(3,1)+p%T2(5)*EnR(3,2)+p%T2(6)*EnR(4,1)+p%T2(7)*EnR(4,2)&
         +p%T2(8)*EnR(5,1)+p%T2(9)*EnR(5,2)
      TT3=p%T3(1)*d1r2+p%T3(2)+p%T3(3)*EnR(2,1)+p%T3(4)*EnR(3,1)+p%T3(5)*EnR(3,2)+p%T3(6)*EnR(4,1)+p%T3(7)*EnR(4,2)&
         +p%T3(8)*EnR(5,1)+p%T3(9)*EnR(5,2)
      do il=1,3
        do ik=1,3
          fs_u=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
          fs_t=TT1*drdx(il)*drdx(ik)*drdn+TT2*(drdn*c_dkr(il,ik)+drdx(ik)*n(il))+TT3*drdx(il)*n(ik)
          h(:,il,ik)=h(:,il,ik)+fs_t*pphijw
          g(:,il,ik)=g(:,il,ik)+fs_u*sphijw
        end do
      end do
    end do
    ! Multiply by constants
    h=p%cte_t*h
    g=p%cte_u*g
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_harela3d_sbie_ext_pre

  subroutine fbem_bem_harela3d_sbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,p,gln,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                                !! Integration element
    logical                            :: reverse                          !! Reverse normal vector
    real(kind=real64)                  :: xi_s(2,fbem_n_vertices(e%gtype)) !! Coordinates of the subdivision of the parent element
    real(kind=real64)                  :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)                  :: barxip(2)                        !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)                  :: barr                             !! Telles jacobian at barxip
    type(fbem_bem_harela3d_parameters) :: p                                !! Parameters of the region
    integer                            :: gln                              !! 1D Gauss-Legendre number of integration points (<=32)
    complex(kind=real64)               :: h(e%n_pnodes,3,3)                !! h kernel vector
    complex(kind=real64)               :: g(e%n_snodes,3,3)                !! g kernel vector
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
    real(kind=real64)            :: r, d1r1, d1r2, d1r3, d1r4 ! Distance vector module and its inverse
    real(kind=real64)            :: drdx(3)                   ! Distance vector derivatives with respect to x_k
    real(kind=real64)            :: jg                        ! Geometric jacobian
    real(kind=real64)            :: jw                        ! Jacobians * weights
    real(kind=real64)            :: drdn                      ! Partial derivative of r respect to unit normal
    real(kind=real64)            :: pphijw(e%n_pnodes)        ! Functional shape functions * jw
    real(kind=real64)            :: sphijw(e%n_snodes)        ! Functional shape functions * jw
    integer                      :: il, ik                    ! Counter for load / observation components
    complex(kind=real64)         :: z(2)                      ! Arguments z=-ikr
    complex(kind=real64)         :: EnR(0:6,2)                ! Exponential function decomposition for each wavenumber
    complex(kind=real64)         :: psi, chi                  ! Components of the fundamental solutions
    complex(kind=real64)         :: TT1, TT2, TT3             ! Components of the fundamental solutions
    complex(kind=real64)         :: fs_u, fs_t                ! Fundamental solutions values
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
            d1r3=d1r2*d1r1
            d1r4=d1r3*d1r1
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
            z(1)=-c_im*p%k1*r
            z(2)=-c_im*p%k2*r
            call fbem_zexp_decomposed(2,z,EnR)
            EnR(2,:)=EnR(2,:)*d1r1
            EnR(3,:)=EnR(3,:)*d1r2
            EnR(4,:)=EnR(4,:)*d1r3
            EnR(5,:)=EnR(5,:)*d1r4
            psi=p%psi(1)*d1r1+p%psi(2)+EnR(2,2)+p%psi(3)*EnR(3,1)+p%psi(4)*EnR(3,2)+p%psi(5)*EnR(4,1)+p%psi(6)*EnR(4,2)
            chi=p%chi(1)*d1r1+p%chi(2)*EnR(2,1)+EnR(2,2)+p%chi(3)*EnR(3,1)+p%chi(4)*EnR(3,2)+p%chi(5)*EnR(4,1)+p%chi(6)*EnR(4,2)
            TT1=p%T1(1)*d1r2+p%T1(2)+p%T1(3)*EnR(2,1)+p%T1(4)*EnR(2,2)+p%T1(5)*EnR(3,1)+p%T1(6)*EnR(3,2)+p%T1(7)*EnR(4,1)&
               +p%T1(8)*EnR(4,2)+p%T1(9)*EnR(5,1)+p%T1(10)*EnR(5,2)
            TT2=p%T2(1)*d1r2+p%T2(2)+p%T2(3)*EnR(2,2)+p%T2(4)*EnR(3,1)+p%T2(5)*EnR(3,2)+p%T2(6)*EnR(4,1)+p%T2(7)*EnR(4,2)&
               +p%T2(8)*EnR(5,1)+p%T2(9)*EnR(5,2)
            TT3=p%T3(1)*d1r2+p%T3(2)+p%T3(3)*EnR(2,1)+p%T3(4)*EnR(3,1)+p%T3(5)*EnR(3,2)+p%T3(6)*EnR(4,1)+p%T3(7)*EnR(4,2)&
               +p%T3(8)*EnR(5,1)+p%T3(9)*EnR(5,2)
            ! Add
            do il=1,3
              do ik=1,3
                fs_u=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
                fs_t=TT1*drdx(il)*drdx(ik)*drdn+TT2*(drdn*c_dkr(il,ik)+drdx(ik)*n(il))+TT3*drdx(il)*n(ik)
                h(:,il,ik)=h(:,il,ik)+fs_t*pphijw
                g(:,il,ik)=g(:,il,ik)+fs_u*sphijw
              end do
            end do
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
            d1r4=d1r3*d1r1
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
            z(1)=-c_im*p%k1*r
            z(2)=-c_im*p%k2*r
            call fbem_zexp_decomposed(2,z,EnR)
            EnR(2,:)=EnR(2,:)*d1r1
            EnR(3,:)=EnR(3,:)*d1r2
            EnR(4,:)=EnR(4,:)*d1r3
            EnR(5,:)=EnR(5,:)*d1r4
            psi=p%psi(1)*d1r1+p%psi(2)+EnR(2,2)+p%psi(3)*EnR(3,1)+p%psi(4)*EnR(3,2)+p%psi(5)*EnR(4,1)+p%psi(6)*EnR(4,2)
            chi=p%chi(1)*d1r1+p%chi(2)*EnR(2,1)+EnR(2,2)+p%chi(3)*EnR(3,1)+p%chi(4)*EnR(3,2)+p%chi(5)*EnR(4,1)+p%chi(6)*EnR(4,2)
            TT1=p%T1(1)*d1r2+p%T1(2)+p%T1(3)*EnR(2,1)+p%T1(4)*EnR(2,2)+p%T1(5)*EnR(3,1)+p%T1(6)*EnR(3,2)+p%T1(7)*EnR(4,1)&
               +p%T1(8)*EnR(4,2)+p%T1(9)*EnR(5,1)+p%T1(10)*EnR(5,2)
            TT2=p%T2(1)*d1r2+p%T2(2)+p%T2(3)*EnR(2,2)+p%T2(4)*EnR(3,1)+p%T2(5)*EnR(3,2)+p%T2(6)*EnR(4,1)+p%T2(7)*EnR(4,2)&
               +p%T2(8)*EnR(5,1)+p%T2(9)*EnR(5,2)
            TT3=p%T3(1)*d1r2+p%T3(2)+p%T3(3)*EnR(2,1)+p%T3(4)*EnR(3,1)+p%T3(5)*EnR(3,2)+p%T3(6)*EnR(4,1)+p%T3(7)*EnR(4,2)&
               +p%T3(8)*EnR(5,1)+p%T3(9)*EnR(5,2)
            ! Add
            do il=1,3
              do ik=1,3
                fs_u=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
                fs_t=TT1*drdx(il)*drdx(ik)*drdn+TT2*(drdn*c_dkr(il,ik)+drdx(ik)*n(il))+TT3*drdx(il)*n(ik)
                h(:,il,ik)=h(:,il,ik)+fs_t*pphijw
                g(:,il,ik)=g(:,il,ik)+fs_u*sphijw
              end do
            end do
          end do
        end do
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'n_edges not valid')
    end select
    ! Multiply by constants
    h=p%cte_t*h
    g=p%cte_u*g
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_harela3d_sbie_ext_st

  recursive subroutine fbem_bem_harela3d_sbie_ext_adp(e,reverse,xi_s,x_i,p,qsp,ks,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                                !! Element
    logical                            :: reverse                          !! Reverse orientation
    real(kind=real64)                  :: xi_s(2,fbem_n_vertices(e%gtype)) !! Subdivision of the parent element
    real(kind=real64)                  :: x_i(3)                           !! Collocation point position vector
    type(fbem_bem_harela3d_parameters) :: p                                !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp                              !! Quasi-singular integration parameters
    integer                            :: ks                               !! Current level of subdivisions
    integer                            :: ns                               !! Maximum level of subdivision
    complex(kind=real64)               :: h(e%n_pnodes,3,3)                !! h integration kernels matrix
    complex(kind=real64)               :: g(e%n_snodes,3,3)                !! g integration kernels matrix
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
    complex(kind=real64) :: h_tmp(e%n_pnodes,3,3)                ! h integration kernels matrix (temporary)
    complex(kind=real64) :: g_tmp(e%n_snodes,3,3)                ! g integration kernels matrix (temporary)
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
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,5,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_harela3d_sbie_ext_adp',ns,'maximum number of subdivisions reached')
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
          call fbem_bem_harela3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
          ! SUBTRI 2
          tmp_xi_s(:,1)=xi_s(:,2)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          call fbem_bem_harela3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
          ! SUBTRI 3
          tmp_xi_s(:,1)=xi_s(:,3)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          call fbem_bem_harela3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
          ! SUBTRI 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          call fbem_bem_harela3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
        ! QUADRILATERALS
        case (4)
          ! SUBQUAD 1
          tmp_xi_s(:,1)=xi_s(:,1)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,3)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          call fbem_bem_harela3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
          ! SUBQUAD 2
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=xi_s(:,2)
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,4)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          call fbem_bem_harela3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
          ! SUBQUAD 3
          tmp_xi_s(:,1)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=xi_s(:,3)
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          call fbem_bem_harela3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
          ! SUBQUAD 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          tmp_xi_s(:,2)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=xi_s(:,4)
          call fbem_bem_harela3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
      end select
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_harela3d_sbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,p,gln,h_tmp,g_tmp)
      h=h+h_tmp
      g=g+g_tmp
    end if
  end subroutine fbem_bem_harela3d_sbie_ext_adp

  subroutine fbem_bem_harela3d_sbie_int(type_g,type_f1,type_f2,delta_f,x_nodes,reverse,xi_i,p,h,g)
    implicit none
    ! I/O
    integer                            :: type_g                           !! Geometrial interpolation
    integer                            :: type_f1                          !! Functional interpolation (primary variables)
    integer                            :: type_f2                          !! Functional interpolation (secondary variables)
    real(kind=real64)                  :: delta_f                          !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64)                  :: x_nodes(3,fbem_n_nodes(type_g))  !! Position vectors of geometrical nodes
    logical                            :: reverse                          !! Normal vector inversion
    real(kind=real64)                  :: xi_i(2)                          !! Reference coordinates of the singular point.
    type(fbem_bem_harela3d_parameters) :: p                                !! Parameters of the region
    complex(kind=real64)               :: h(fbem_n_nodes(type_f1),3,3)     !! h kernel vector
    complex(kind=real64)               :: g(fbem_n_nodes(type_f2),3,3)     !! g kernel vector
    ! Local
    integer                            :: ksubtri                          ! Counter variable for subtriangles loop
    integer                            :: nsubtri                          ! Number of subtriangles
    integer                            :: subtriangle(8)                   ! Vector that contains what subtriangles need to be integrated
    real(kind=real64)                  :: theta_subtri(2,8)                ! Matrix that contains the angles theta of the subtriangles to be integrated
    real(kind=real64)                  :: thetap_subtri(2,8)               ! Matrix that contains the angles thetap of the subtriangles to be integrated
    integer                            :: ktheta                           ! Counter variable for theta coordinate loop
    integer                            :: krho                             ! Counter variable for rho coordinate loop
    integer                            :: kphi                             ! Counter coordinates for shape functions loops
    integer                            :: nnodes_g                         ! Number of nodes of geometrical interpolation
    integer                            :: ngp_theta                        ! Number of Gauss points for theta coordinate
    integer                            :: ngp_rho                          ! Number of Gauss points for rho coordinate
    real(kind=real64)                  :: thetai, thetaf, thetapi, thetapf ! Initial and final angles for subtriangle integration
    real(kind=real64)                  :: w_angular                        ! Weight of the angular coordinate
    real(kind=real64)                  :: w_radial                         ! Weight of the radial coordinate
    real(kind=real64)                  :: theta                            ! Angle coordinate theta
    real(kind=real64)                  :: thetap                           ! Angle coordinate thetap
    real(kind=real64)                  :: thetapp                          ! Angle coordinate thetap on [0,1] domain
    real(kind=real64)                  :: jthetap                          ! thetap->thetapp jacobian
    real(kind=real64)                  :: rhoij                            ! Maximum rho (radial) value for each edge
    real(kind=real64)                  :: rho                              ! Radial coordinate rho
    real(kind=real64)                  :: rhop                             ! Radial coordinate rho on [0,1] domain
    real(kind=real64)                  :: aux(10)                          ! Auxiliary variable for shape functions resources
    real(kind=real64)                  :: xi(2)                            ! Reference xi_1,xi_2 coordinates
    real(kind=real64)                  :: phi_f1(fbem_n_nodes(type_f1))    ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)                  :: phi_f2(fbem_n_nodes(type_f2))    ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)                  :: phi_f1_i(fbem_n_nodes(type_f1))  ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)                  :: phi_g(fbem_n_nodes(type_g))      ! Geometrical shape functions values at xi_1,xi_2
    real(kind=real64)                  :: dphidxi1_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64)                  :: dphidxi2_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64)                  :: x_i(3)                           ! Collocation point position vector
    real(kind=real64)                  :: x(3)                             ! Position vector at xi_1,xi_2
    real(kind=real64)                  :: T1(3), T2(3)                     ! Tangent vectors at xi_1,xi_2
    real(kind=real64)                  :: N(3)                             ! Normal vector at xi_1,xi_2
    real(kind=real64)                  :: rv(3)                            ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)                  :: r, dr1, dr2, dr3, dr4            ! Distance vector module and its inverse
    real(kind=real64)                  :: drdx(3)                          ! Distance vector derivatives with respect to x_k
    real(kind=real64)                  :: jg                               ! Geometric jacobian
    real(kind=real64)                  :: jw                               ! Jacobian * weights
    real(kind=real64)                  :: drdn                             ! Partial derivative of r respect to unit normal
    real(kind=real64)                  :: phif1jw(fbem_n_nodes(type_f1))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)                  :: phif2jw(fbem_n_nodes(type_f2))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)                  :: costheta, sintheta               ! cos(theta), sin(theta)
    integer                            :: il, ik                           ! Counter for load / observation components
    complex(kind=real64)               :: z(2)                             ! Wavenumbers
    complex(kind=real64)               :: E(0:6,2)                         ! Exponential function decomposition for each wavenumber
    complex(kind=real64)               :: psi, chi, TT1, TT2, TT3          ! Components of the fundamental solutions
    complex(kind=real64)               :: fs_u, fs_t                       ! Fundamental solutions values
    ! Local variables associated with line integrals
    real(kind=real64)                  :: re                               ! Quasi-singular integration relative error
    integer                            :: gln_min                          ! Minimum Number of 1D Gauss-Legendre number of points
    integer                            :: ns_max                           ! Maximum number of subdivisions
    type(fbem_qs_parameters)           :: qsp                              ! Quasi-singular integration parameters
    real(kind=real64)                  :: xi_s(1,2)                        ! Edge subdivision
    integer                            :: kedge                            ! Counter variable for edges
    integer                            :: nedges                           ! Number of edges
    logical                            :: integrate                        ! Control variable
    integer                            :: type_edge                        ! Line element type for line integrals
    integer                            :: nnodes_edge                      ! Number of nodes of the edge
    real(kind=real64), allocatable     :: x_nodes_edge(:,:)                ! Coordinates of the edge elements
    real(kind=real64)                  :: hli(3,3)                         ! Line integrals values
    !
    ! Initialization
    !
    ! Kernel vectors
    h=(0.d0,0.d0)
    g=(0.d0,0.d0)
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
    ! Functional shape functions at xi_i
#   define etype type_f1
#   define delta delta_f
#   define xi xi_i
#   define phi phi_f1_i
#   include <phi_2d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
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
          ! Unit normal vector
          n=N/jg
          ! Distance vector
          rv=x-x_i
          ! Distance vector norm
          r=dot_product(rv,rv)
          r=sqrt(r)
          dr1=1.d0/r
          dr2=dr1**2
          dr3=dr2*dr1
          dr4=dr3*dr1
          ! r_{,k}
          drdx=rv*dr1
          ! dr/dn
          drdn=dot_product(drdx,n)
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
          ! COMPONENTS OF THE FUNDAMENTAL SOLUTION (ONLY REGULAR PARTS)
          z(1)=-c_im*p%k1*r
          z(2)=-c_im*p%k2*r
          call fbem_zexp_decomposed(2,z,E)
          E(2,:)=E(2,:)*dr1
          E(3,:)=E(3,:)*dr2
          E(4,:)=E(4,:)*dr3
          E(5,:)=E(5,:)*dr4
          psi=p%psi(1)*dr1+p%psi(2)+E(2,2)+p%psi(3)*E(3,1)+p%psi(4)*E(3,2)+p%psi(5)*E(4,1)+p%psi(6)*E(4,2)
          chi=p%chi(1)*dr1+p%chi(2)*E(2,1)+E(2,2)+p%chi(3)*E(3,1)+p%chi(4)*E(3,2)+p%chi(5)*E(4,1)+p%chi(6)*E(4,2)
          TT1=p%T1(2)+p%T1(3)*E(2,1)+p%T1(4)*E(2,2)+p%T1(5)*E(3,1)+p%T1(6)*E(3,2)+p%T1(7)*E(4,1)+p%T1(8)*E(4,2)&
             +p%T1(9)*E(5,1)+p%T1(10)*E(5,2)
          TT2=p%T2(2)+p%T2(3)*E(2,2)+p%T2(4)*E(3,1)+p%T2(5)*E(3,2)+p%T2(6)*E(4,1)+p%T2(7)*E(4,2)+p%T2(8)*E(5,1)&
             +p%T2(9)*E(5,2)
          TT3=p%T3(2)+p%T3(3)*E(2,1)+p%T3(4)*E(3,1)+p%T3(5)*E(3,2)+p%T3(6)*E(4,1)+p%T3(7)*E(4,2)+p%T3(8)*E(5,1)&
             +p%T3(9)*E(5,2)
          do il=1,3
            do ik=1,3
              ! Regular parts of fundamental solutions values without constants
              fs_u=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
              fs_t=(TT1+p%T1(1)*dr2)*drdx(il)*drdx(ik)*drdn+(TT2+p%T2(1)*dr2)*drdn*c_dkr(il,ik)+TT2*drdx(ik)*n(il)&
                  +TT3*drdx(il)*n(ik)
              ! Add to kernels
              h(:,il,ik)=h(:,il,ik)+fs_t*phif1jw
              g(:,il,ik)=g(:,il,ik)+fs_u*phif2jw
              ! The regular part of the CPV integral of h
              fs_t=p%T2(1)*dr2*(n(il)*drdx(ik)-n(ik)*drdx(il))
              h(:,il,ik)=h(:,il,ik)+fs_t*(phi_f1-phi_f1_i)*jw
            end do
          end do
        end do ! Loop through rho coordinate
      end do ! Loop through theta coordinate
    end do ! Loop through rectangle triangles
    !
    ! LINE INTEGRALS (all quasi-singular integrals, all integrations performed with an adaptive Telles+subdivision)
    !
    ! Initialize
    hli=0.d0
    re=1.d-15
    gln_min=5
    ns_max=16
    call fbem_qs_calculate_parameters(re,qsp)
    ! Number of edges of the element
    nedges=fbem_n_edges(type_g)
    ! Loop through edges
    do kedge=1,nedges
      ! Check if it is necessary to calculate it
      integrate=.false.
      do ksubtri=1,nsubtri
        if ((subtriangle(ksubtri).eq.(2*kedge-1)).or.(subtriangle(ksubtri).eq.(2*kedge))) then
          integrate=.true.
          exit
        end if
      end do
      ! Integrate if needed
      if (integrate) then
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
        ! Calculate the line integral epsilon_{ijk}*e_k·t/r
        call fbem_bem_staela3d_sbie_int_li(type_edge,x_nodes_edge,xi_s,x_i,gln_min,qsp,1,ns_max,hli)
        ! Deallocate the edge coordinates
        deallocate (x_nodes_edge)
      end if
    end do ! Loop through edges
    !
    ! Add line integrals
    !
    ! Loop through load direction and observation direction
    do il=1,3
      do ik=1,3
        h(:,il,ik)=h(:,il,ik)+phi_f1_i*p%T2(1)*hli(il,ik)
      end do
    end do
    ! Multiply h and g by constants of t* and u* respectively
    h=p%cte_t*h
    g=p%cte_u*g
    ! If the normal has to be reversed, then h=-h
    if (reverse) h=-h
  end subroutine fbem_bem_harela3d_sbie_int

  subroutine fbem_bem_harela3d_sbie_auto(e,reverse,x_i,p,qsp,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                 !! Integration element
    logical                            :: reverse           !! Reverse orientation
    real(kind=real64)                  :: x_i(3)            !! Collocation point
    type(fbem_bem_harela3d_parameters) :: p                 !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp               !! Quasi-singular integration parameters
    integer                            :: ns                !! Maximum level of subdivisions
    complex(kind=real64)               :: h(e%n_pnodes,3,3) !! h integration kernel
    complex(kind=real64)               :: g(e%n_snodes,3,3) !! g integration kernel
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
        call fbem_bem_harela3d_sbie_int(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi,p,h,g)
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
          call fbem_bem_harela3d_sbie_ext_pre(ps,e,reverse,x_i,p,h,g)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_harela3d_sbie_ext_adp(e,reverse,xi_s,x_i,p,qsp,1,ns,h,g)
        end if
    end select
  end subroutine fbem_bem_harela3d_sbie_auto
  ! --------------------------------------------------------------------------------------------------------------------------------

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! BODY LOAD ELEMENTS

  subroutine fbem_bem_harela3d_sbie_bl_ext_pre(ps,e,x_i,p,g)
    implicit none
    ! I/O
    integer                            :: ps                !! Selected precalculated dataset
    type(fbem_bem_element)             :: e                 !! Element
    real(kind=real64)                  :: x_i(3)            !! Collocation point position vector
    type(fbem_bem_harela3d_parameters) :: p                 !! Parameters of the region
    complex(kind=real64)               :: g(e%n_snodes,3,3) !! g integration kernels matrix
    ! Local
    integer              :: kip                        ! Counter variable for integration points loop
    real(kind=real64)    :: x(3)                       ! Position vector at integration point
    real(kind=real64)    :: n(3)                       ! Unit normal vector at integration point
    real(kind=real64)    :: sphijw(e%n_snodes)         ! phi^s * jacobian * weight at integration point
    real(kind=real64)    :: rv(3)                      ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, d1r1, d1r2, d1r3, d1r4  ! Distance vector module and its inverse
    real(kind=real64)    :: drdx(3)                    ! Distance vector derivatives with respect to x_k
    real(kind=real64)    :: drdn                       ! Partial derivative of r respect to unit normal
    integer              :: il, ik                     ! Counter for load / observation components
    complex(kind=real64) :: z(2)                       ! Arguments z=-ikr
    complex(kind=real64) :: EnR(0:6,2)                 ! Exponential function decomposition for each wavenumber
    complex(kind=real64) :: psi, chi                   ! Components of the fundamental solutions
    complex(kind=real64) :: fs_u                       ! Fundamental solutions values
    ! Initialize
    g=(0.d0,0.d0)
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      x=e%ps_x(:,kip,ps)
      n=e%ps_n(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      d1r1=1.d0/r
      d1r2=d1r1**2
      d1r3=d1r2*d1r1
      d1r4=d1r3*d1r1
      drdx=rv*d1r1
      drdn=dot_product(drdx,n)
      z(1)=-c_im*p%k1*r
      z(2)=-c_im*p%k2*r
      call fbem_zexp_decomposed(2,z,EnR)
      EnR(2,:)=EnR(2,:)*d1r1
      EnR(3,:)=EnR(3,:)*d1r2
      EnR(4,:)=EnR(4,:)*d1r3
      EnR(5,:)=EnR(5,:)*d1r4
      psi=p%psi(1)*d1r1+p%psi(2)+EnR(2,2)+p%psi(3)*EnR(3,1)+p%psi(4)*EnR(3,2)+p%psi(5)*EnR(4,1)+p%psi(6)*EnR(4,2)
      chi=p%chi(1)*d1r1+p%chi(2)*EnR(2,1)+EnR(2,2)+p%chi(3)*EnR(3,1)+p%chi(4)*EnR(3,2)+p%chi(5)*EnR(4,1)+p%chi(6)*EnR(4,2)
      do il=1,3
        do ik=1,3
          fs_u=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
          g(:,il,ik)=g(:,il,ik)+fs_u*sphijw
        end do
      end do
    end do
    g=p%cte_u*g
  end subroutine fbem_bem_harela3d_sbie_bl_ext_pre

  subroutine fbem_bem_harela3d_sbie_bl_ext_st(e,xi_s,x_i,barxip,barr,p,gln,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                                !! Integration element
    real(kind=real64)                  :: xi_s(2,fbem_n_vertices(e%gtype)) !! Coordinates of the subdivision of the parent element
    real(kind=real64)                  :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)                  :: barxip(2)                        !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)                  :: barr                             !! Telles jacobian at barxip
    type(fbem_bem_harela3d_parameters) :: p                                !! Parameters of the region
    integer                            :: gln                              !! 1D Gauss-Legendre number of integration points (<=32)
    complex(kind=real64)               :: g(e%n_snodes,3,3)                !! g kernel vector
    ! Local
    integer                      :: kphi                      ! Counter variable for shape functions loops
    integer                      :: k1                        ! Counter variable for reference coordinate xi_1
    integer                      :: k2                        ! Counter variable for reference coordinate xi_2
    real(kind=real64)            :: aux(10)                   ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)          ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi(e%n_gnodes)      ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: dgphidxi1(e%n_gnodes)     ! Geometrical shape functions derivatives values
    real(kind=real64)            :: dgphidxi2(e%n_gnodes)     ! Geometrical shape functions derivatives values
    real(kind=real64)            :: sphi(e%n_snodes)          ! Functional shape functions values
    real(kind=real64)            :: gamma(e%d)                ! Vector of gamma_1,gamma_2 coordinates
    real(kind=real64)            :: w(e%d)                    ! Weights of the integration rule
    real(kind=real64)            :: xip(e%d)                  ! Vector of xip_1,xip_2 coordinates
    real(kind=real64)            :: dxidxi1p(2), dxidxi2p(2)  ! xi derivatives with respect to xip
    real(kind=real64)            :: js                        ! Subdivision jacobian
    real(kind=real64)            :: xin(e%d)                  ! Coordinate xi [xi_s(1,1),xi_s(1,2)]
    real(kind=real64)            :: xipp(2)                   ! Coordinate xipp used for quadrilateral-triangle transformation
    real(kind=real64)            :: barxipp(2)                ! Coordinate xipp of collocation point
    real(kind=real64)            :: jqt                       ! Jacobian of the quadrilateral-triangle transformation
    type(fbem_telles_parameters) :: telles_parameters(e%d)    ! Telles parameters for each coordinate
    real(kind=real64)            :: jt(e%d)                   ! Telles jacobian for each coordinate: xi_1->gamma_1 and xi_2->gamma_2
    real(kind=real64)            :: x(3)                      ! Position vector at xi_1,xi_2
    real(kind=real64)            :: T(3), T1(3), T2(3)              ! Tangent vectors at xi_1,xi_2
    real(kind=real64)            :: N(3)                      ! Normal vector at xi_1,xi_2
    real(kind=real64)            :: rv(3)                     ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)            :: r, d1r1, d1r2, d1r3, d1r4 ! Distance vector module and its inverse
    real(kind=real64)            :: drdx(3)                   ! Distance vector derivatives with respect to x_k
    real(kind=real64)            :: jg                        ! Geometric jacobian
    real(kind=real64)            :: jw                        ! Jacobians * weights
    real(kind=real64)            :: drdn                      ! Partial derivative of r respect to unit normal
    real(kind=real64)            :: sphijw(e%n_snodes)        ! Functional shape functions * jw
    integer                      :: il, ik                    ! Counter for load / observation components
    complex(kind=real64)         :: z(2)                      ! Arguments z=-ikr
    complex(kind=real64)         :: EnR(0:6,2)                ! Exponential function decomposition for each wavenumber
    complex(kind=real64)         :: psi, chi                  ! Components of the fundamental solutions
    complex(kind=real64)         :: fs_u                      ! Fundamental solutions values
    ! Initialization
    g=(0.d0,0.d0)
    ! Calculate
    select case (e%d)
      !
      ! LINE LOAD
      !
      case (1)
        ! Calculate Telles parameters
        telles_parameters(1)=fbem_telles11_calculate_parameters(barxip(1),barr)
        ! Numerical integration
        do k1=1,gl11_n(gln)
          ! GAMMA COORDINATE
          gamma(1)=gl11_xi(k1,gln)
          w(1)=gl11_w(k1,gln)
          ! GAMMA->XIP TRANSFORMATION
          call fbem_telles_xi_and_jacobian(telles_parameters(1),gamma(1),xip(1),jt(1))
          ! XIP->XI TRANSFORMATION
          xin(1)=0.5d0*(1.d0-xip(1))*xi_s(1,1)+0.5d0*(1.d0+xip(1))*xi_s(1,2)
          js=0.5d0*(xi_s(1,2)-xi_s(1,1))
          ! XI->X TRANSFORMATION
#         define etype e%gtype
#         define delta 0.d0
#         define xi xin(1)
#         define phi gphi
#         define dphidxi dgphidxi
#         include <phi_and_dphidxi_1d.rc>
#         undef etype
#         undef delta
#         undef xi
#         undef phi
#         undef dphidxi
          x=0.d0
          T=0.d0
          do kphi=1,e%n_gnodes
            x=x+gphi(kphi)*e%x(:,kphi)
            T=T+dgphidxi(kphi)*e%x(:,kphi)
          end do
          ! Geometric jacobian
          jg=sqrt(dot_product(T,T))
          ! Distance and functions of distance and unit normal
          rv=x-x_i
          r=sqrt(dot_product(rv,rv))
          d1r1=1.d0/r
          drdx=rv*d1r1
          ! Jacobians * weight
          jw=jg*js*jt(1)*w(1)
          ! FUNCTIONAL SHAPE FUNCTIONS
          ! Functional shape functions (body load variables) at xi
#         define etype e%stype
#         define delta e%stype_delta
#         define xi xin(1)
#         define phi sphi
#         include <phi_1d.rc>
#         undef etype
#         undef delta
#         undef xi
#         undef phi
          ! Functional shape functions * jacobians* weights
          sphijw=sphi*jw
          ! COMPONENTS OF THE FUNDAMENTAL SOLUTION
          z(1)=-c_im*p%k1*r
          z(2)=-c_im*p%k2*r
          call fbem_zexp_decomposed(2,z,EnR)
          EnR(2,:)=EnR(2,:)*d1r1
          EnR(3,:)=EnR(3,:)*d1r2
          EnR(4,:)=EnR(4,:)*d1r3
          EnR(5,:)=EnR(5,:)*d1r4
          psi=p%psi(1)*d1r1+p%psi(2)+EnR(2,2)+p%psi(3)*EnR(3,1)+p%psi(4)*EnR(3,2)+p%psi(5)*EnR(4,1)+p%psi(6)*EnR(4,2)
          chi=p%chi(1)*d1r1+p%chi(2)*EnR(2,1)+EnR(2,2)+p%chi(3)*EnR(3,1)+p%chi(4)*EnR(3,2)+p%chi(5)*EnR(4,1)+p%chi(6)*EnR(4,2)
          ! Add
          do il=1,3
            do ik=1,3
              fs_u=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
              g(:,il,ik)=g(:,il,ik)+fs_u*sphijw
            end do
          end do
        end do
      !
      ! SURFACE LOAD
      !
      case (2)
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
#               define delta 0.d0
#               define xi xip
#               define phi gphi
#               define dphidxi1 dgphidxi1
#               define dphidxi2 dgphidxi2
#               include <phi_quad4.rc>
#               include <dphidxi1_quad4.rc>
#               include <dphidxi2_quad4.rc>
#               undef delta
#               undef xi
#               undef phi
#               undef dphidxi1
#               undef dphidxi2
                ! xi coordinates, and xi derivatives
                xin=0.d0
                dxidxi1p=0.d0
                dxidxi2p=0.d0
                do kphi=1,4
                  xin=xin+gphi(kphi)*xi_s(:,kphi)
                  dxidxi1p=dxidxi1p+dgphidxi1(kphi)*xi_s(:,kphi)
                  dxidxi2p=dxidxi2p+dgphidxi2(kphi)*xi_s(:,kphi)
                end do
                ! xip->xi jacobian
                js=dxidxi1p(1)*dxidxi2p(2)-dxidxi1p(2)*dxidxi2p(1)
                ! XI->X TRANSFORMATION
                ! Geometrical shape functions and first derivatives at xi
#               define etype e%gtype
#               define xi xin
#               define delta 0.d0
#               define phi gphi
#               define dphidxi1 dgphidxi1
#               define dphidxi2 dgphidxi2
#               include <phi_and_dphidxik_2d.rc>
#               undef etype
#               undef xi
#               undef delta
#               undef phi
#               undef dphidxi1
#               undef dphidxi2
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
                d1r4=d1r3*d1r1
                drdx=rv*d1r1
                drdn=dot_product(drdx,n)
                ! Jacobians * weights
                jw=jg*js*jt(1)*jt(2)*w(1)*w(2)
                ! FUNCTIONAL SHAPE FUNCTIONS
                ! Functional shape functions (secondary variables) at xi
#               define etype e%stype
#               define xi xin
#               define delta e%stype_delta
#               define phi sphi
#               include <phi_2d.rc>
#               undef etype
#               undef xi
#               undef delta
#               undef phi
                ! Functional shape functions * jacobians * weights
                sphijw=sphi*jw
                ! COMPONENTS OF THE FUNDAMENTAL SOLUTION
                z(1)=-c_im*p%k1*r
                z(2)=-c_im*p%k2*r
                call fbem_zexp_decomposed(2,z,EnR)
                EnR(2,:)=EnR(2,:)*d1r1
                EnR(3,:)=EnR(3,:)*d1r2
                EnR(4,:)=EnR(4,:)*d1r3
                EnR(5,:)=EnR(5,:)*d1r4
                psi=p%psi(1)*d1r1+p%psi(2)+EnR(2,2)+p%psi(3)*EnR(3,1)+p%psi(4)*EnR(3,2)+p%psi(5)*EnR(4,1)+p%psi(6)*EnR(4,2)
                chi=p%chi(1)*d1r1+p%chi(2)*EnR(2,1)+EnR(2,2)+p%chi(3)*EnR(3,1)+p%chi(4)*EnR(3,2)+p%chi(5)*EnR(4,1)+p%chi(6)*EnR(4,2)
                ! Add
                do il=1,3
                  do ik=1,3
                    fs_u=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
                    g(:,il,ik)=g(:,il,ik)+fs_u*sphijw
                  end do
                end do
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
#               define delta 0.d0
#               define xi xip
#               define phi gphi
#               define dphidxi1 dgphidxi1
#               define dphidxi2 dgphidxi2
#               include <phi_tri3.rc>
#               include <dphidxi1_tri3.rc>
#               include <dphidxi2_tri3.rc>
#               undef delta
#               undef xi
#               undef phi
#               undef dphidxi1
#               undef dphidxi2
                ! xi coordinates, and xi derivatives
                xin=0.d0
                dxidxi1p=0.d0
                dxidxi2p=0.d0
                do kphi=1,3
                  xin=xin+gphi(kphi)*xi_s(:,kphi)
                  dxidxi1p=dxidxi1p+dgphidxi1(kphi)*xi_s(:,kphi)
                  dxidxi2p=dxidxi2p+dgphidxi2(kphi)*xi_s(:,kphi)
                end do
                ! xip->xi jacobian
                js=dxidxi1p(1)*dxidxi2p(2)-dxidxi1p(2)*dxidxi2p(1)
                ! XI->X transformation
                ! Geometrical shape functions and first derivatives at xi
#               define etype e%gtype
#               define xi xin
#               define delta 0.d0
#               define phi gphi
#               define dphidxi1 dgphidxi1
#               define dphidxi2 dgphidxi2
#               include <phi_and_dphidxik_2d.rc>
#               undef etype
#               undef xi
#               undef delta
#               undef phi
#               undef dphidxi1
#               undef dphidxi2
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
                d1r4=d1r3*d1r1
                drdx=rv*d1r1
                drdn=dot_product(drdx,n)
                ! Jacobians * weights
                jw=jg*js*jqt*jt(1)*jt(2)*w(1)*w(2)
                ! FUNCTIONAL SHAPE FUNCTIONS
                ! Functional shape functions (secondary variables) at xi
#               define etype e%stype
#               define xi xin
#               define delta e%stype_delta
#               define phi sphi
#               include <phi_2d.rc>
#               undef etype
#               undef xi
#               undef delta
#               undef phi
                ! Functional shape functions * jacobians * weights
                sphijw=sphi*jw
                ! COMPONENTS OF THE FUNDAMENTAL SOLUTION
                z(1)=-c_im*p%k1*r
                z(2)=-c_im*p%k2*r
                call fbem_zexp_decomposed(2,z,EnR)
                EnR(2,:)=EnR(2,:)*d1r1
                EnR(3,:)=EnR(3,:)*d1r2
                EnR(4,:)=EnR(4,:)*d1r3
                EnR(5,:)=EnR(5,:)*d1r4
                psi=p%psi(1)*d1r1+p%psi(2)+EnR(2,2)+p%psi(3)*EnR(3,1)+p%psi(4)*EnR(3,2)+p%psi(5)*EnR(4,1)+p%psi(6)*EnR(4,2)
                chi=p%chi(1)*d1r1+p%chi(2)*EnR(2,1)+EnR(2,2)+p%chi(3)*EnR(3,1)+p%chi(4)*EnR(3,2)+p%chi(5)*EnR(4,1)+p%chi(6)*EnR(4,2)
                ! Add
                do il=1,3
                  do ik=1,3
                    fs_u=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
                    g(:,il,ik)=g(:,il,ik)+fs_u*sphijw
                  end do
                end do
              end do
            end do
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'n_edges not valid')
        end select
      !
      ! VOLUME LOADS
      !
      case (3)
        stop 'not yet : fbem_bem_harela3d_sbie_bl_ext_st'
      !
      ! OTHERS
      !
      case default
        call fbem_error_message(output_unit,0,'fbem_bem_harela3d_sbie_bl_ext_st',0,&
                                'it is only possible to integrate line, surface or volume loads')
    end select
    g=p%cte_u*g
  end subroutine fbem_bem_harela3d_sbie_bl_ext_st

  recursive subroutine fbem_bem_harela3d_sbie_bl_ext_adp(e,xi_s,x_i,p,qsp,ks,ns,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                                !! Element
    real(kind=real64)                  :: xi_s(2,fbem_n_vertices(e%gtype)) !! Subdivision of the parent element
    real(kind=real64)                  :: x_i(3)                           !! Collocation point position vector
    type(fbem_bem_harela3d_parameters) :: p                                !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp                              !! Quasi-singular integration parameters
    integer                            :: ks                               !! Current level of subdivisions
    integer                            :: ns                               !! Maximum level of subdivision
    complex(kind=real64)               :: g(e%n_snodes,3,3)                !! g integration kernels matrix
    ! Local
    integer              :: gln_near                               ! 1D Gauss-Legendre integ. points required to integrate only  the quasi-singular integrand
    integer              :: gln                                    ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical              :: subdivide                              ! True if subdivision has to be performed
    real(kind=real64)    :: barxi(e%d)                             ! Nearest element coordinates with respect to collocation point
    real(kind=real64)    :: barxip(e%d)                            ! Nearest element subdivision local coordinates with respect to collocation point
    real(kind=real64)    :: rmin                                   ! Minimum distance between collocation point and barxi on the element
    real(kind=real64)    :: barr                                   ! Telles jacobian at barxi
    real(kind=real64)    :: cl                                     ! Characteristic length
    real(kind=real64)    :: d                                      ! Normalized distance between collocation point and element subdivision
    integer              :: method                                 ! Method used in nearest point algorithm
    real(kind=real64)    :: tmp_xi_s(e%d,fbem_n_vertices(e%gtype)) ! Subdivision
    real(kind=real64)    :: x_s(3,e%n_gnodes)                      ! Coordinates of the element subdivision
    complex(kind=real64) :: g_tmp(e%n_snodes,3,3)                  ! g integration kernels matrix (temporary)
    ! Initialize
    if (ks.eq.1) then
      g=(0.d0,0.d0)
      select case (e%d)
        case (1)
          xi_s(1,1)=-1.d0
          xi_s(1,2)= 1.d0
        case (2)
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
        case (3)
          stop 'not yet : fbem_bem_staela3d_sbie_bl_ext_adp'
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
        call fbem_warning_message(error_unit,0,'fbem_bem_harela3d_sbie_bl_ext_adp',ns,'maximum number of subdivisions reached')
        gln_near=30
      end if
    else
      if (gln_near.eq.0) subdivide=.true.
    end if
    ! Subdivide by 1/4
    if (subdivide) then
      select case (e%d)
        !
        ! LINE LOAD
        !
        case (1)
          ! 1D
          ! SUBLINE 1
          tmp_xi_s(1,1)=xi_s(1,1)
          tmp_xi_s(1,2)=0.5d0*(xi_s(1,1)+xi_s(1,2))
          call fbem_bem_harela3d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,p,qsp,ks+1,ns,g)
          ! SUBLINE 2
          tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
          tmp_xi_s(1,2)=xi_s(1,2)
          call fbem_bem_harela3d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,p,qsp,ks+1,ns,g)
        !
        ! SURFACE LOAD
        !
        case (2)
          select case (fbem_n_vertices(e%gtype))
            ! TRIANGULAR ELEMENTS
            case (3)
              ! SUBTRI 1
              tmp_xi_s(:,1)=xi_s(:,1)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
              call fbem_bem_harela3d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,p,qsp,ks+1,ns,g)
              ! SUBTRI 2
              tmp_xi_s(:,1)=xi_s(:,2)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              call fbem_bem_harela3d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,p,qsp,ks+1,ns,g)
              ! SUBTRI 3
              tmp_xi_s(:,1)=xi_s(:,3)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              call fbem_bem_harela3d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,p,qsp,ks+1,ns,g)
              ! SUBTRI 4
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
              call fbem_bem_harela3d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,p,qsp,ks+1,ns,g)
            ! QUADRILATERALS
            case (4)
              ! SUBQUAD 1
              tmp_xi_s(:,1)=xi_s(:,1)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,3)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,4)=0.50d0*(xi_s(:,1)+xi_s(:,4))
              call fbem_bem_harela3d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,p,qsp,ks+1,ns,g)
              ! SUBQUAD 2
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,2)=xi_s(:,2)
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,4)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              call fbem_bem_harela3d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,p,qsp,ks+1,ns,g)
              ! SUBQUAD 3
              tmp_xi_s(:,1)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=xi_s(:,3)
              tmp_xi_s(:,4)=0.50d0*(xi_s(:,3)+xi_s(:,4))
              call fbem_bem_harela3d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,p,qsp,ks+1,ns,g)
              ! SUBQUAD 4
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,4))
              tmp_xi_s(:,2)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,4)=xi_s(:,4)
              call fbem_bem_harela3d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,p,qsp,ks+1,ns,g)
          end select
        !
        ! VOLUME LOAD
        !
        case (3)
          stop 'not yet : fbem_bem_staela3d_sbie_bl_ext_adp'
      end select
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_harela3d_sbie_bl_ext_st(e,xi_s,x_i,barxip,barr,p,gln,g_tmp)
      g=g+g_tmp
    end if
  end subroutine fbem_bem_harela3d_sbie_bl_ext_adp

  subroutine fbem_bem_harela3d_sbie_bl_int(type_g,type_f1,type_f2,delta_f,x_nodes,xi_i,p,g)
    implicit none
    ! I/O
    integer                            :: type_g                           !! Geometrial interpolation
    integer                            :: type_f1                          !! Functional interpolation (primary variables)
    integer                            :: type_f2                          !! Functional interpolation (secondary variables)
    real(kind=real64)                  :: delta_f                          !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64)                  :: x_nodes(3,fbem_n_nodes(type_g))  !! Position vectors of geometrical nodes
    real(kind=real64)                  :: xi_i(2)                          !! Reference coordinates of the singular point.
    type(fbem_bem_harela3d_parameters) :: p                                !! Parameters of the region
    complex(kind=real64)               :: g(fbem_n_nodes(type_f2),3,3)     !! g kernel vector
    ! Local
    integer                            :: ksubtri                          ! Counter variable for subtriangles loop
    integer                            :: nsubtri                          ! Number of subtriangles
    integer                            :: subtriangle(8)                   ! Vector that contains what subtriangles need to be integrated
    real(kind=real64)                  :: theta_subtri(2,8)                ! Matrix that contains the angles theta of the subtriangles to be integrated
    real(kind=real64)                  :: thetap_subtri(2,8)               ! Matrix that contains the angles thetap of the subtriangles to be integrated
    integer                            :: ktheta                           ! Counter variable for theta coordinate loop
    integer                            :: krho                             ! Counter variable for rho coordinate loop
    integer                            :: kphi                             ! Counter coordinates for shape functions loops
    integer                            :: nnodes_g                         ! Number of nodes of geometrical interpolation
    integer                            :: ngp_theta                        ! Number of Gauss points for theta coordinate
    integer                            :: ngp_rho                          ! Number of Gauss points for rho coordinate
    real(kind=real64)                  :: thetai, thetaf, thetapi, thetapf ! Initial and final angles for subtriangle integration
    real(kind=real64)                  :: w_angular                        ! Weight of the angular coordinate
    real(kind=real64)                  :: w_radial                         ! Weight of the radial coordinate
    real(kind=real64)                  :: theta                            ! Angle coordinate theta
    real(kind=real64)                  :: thetap                           ! Angle coordinate thetap
    real(kind=real64)                  :: thetapp                          ! Angle coordinate thetap on [0,1] domain
    real(kind=real64)                  :: jthetap                          ! thetap->thetapp jacobian
    real(kind=real64)                  :: rhoij                            ! Maximum rho (radial) value for each edge
    real(kind=real64)                  :: rho                              ! Radial coordinate rho
    real(kind=real64)                  :: rhop                             ! Radial coordinate rho on [0,1] domain
    real(kind=real64)                  :: aux(10)                          ! Auxiliary variable for shape functions resources
    real(kind=real64)                  :: xi(2)                            ! Reference xi_1,xi_2 coordinates
    real(kind=real64)                  :: phi_f2(fbem_n_nodes(type_f2))    ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)                  :: phi_f1_i(fbem_n_nodes(type_f1))  ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)                  :: phi_g(fbem_n_nodes(type_g))      ! Geometrical shape functions values at xi_1,xi_2
    real(kind=real64)                  :: dphidxi1_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64)                  :: dphidxi2_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64)                  :: x_i(3)                           ! Collocation point position vector
    real(kind=real64)                  :: x(3)                             ! Position vector at xi_1,xi_2
    real(kind=real64)                  :: T1(3), T2(3)                     ! Tangent vectors at xi_1,xi_2
    real(kind=real64)                  :: N(3)                             ! Normal vector at xi_1,xi_2
    real(kind=real64)                  :: rv(3)                            ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)                  :: r, dr1, dr2, dr3, dr4            ! Distance vector module and its inverse
    real(kind=real64)                  :: drdx(3)                          ! Distance vector derivatives with respect to x_k
    real(kind=real64)                  :: jg                               ! Geometric jacobian
    real(kind=real64)                  :: jw                               ! Jacobian * weights
    real(kind=real64)                  :: drdn                             ! Partial derivative of r respect to unit normal
    real(kind=real64)                  :: phif2jw(fbem_n_nodes(type_f2))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)                  :: costheta, sintheta               ! cos(theta), sin(theta)
    integer                            :: il, ik                           ! Counter for load / observation components
    complex(kind=real64)               :: z(2)                             ! Wavenumbers
    complex(kind=real64)               :: E(0:6,2)                         ! Exponential function decomposition for each wavenumber
    complex(kind=real64)               :: psi, chi                         ! Components of the fundamental solutions
    complex(kind=real64)               :: fs_u                             ! Fundamental solutions values
    !
    ! Initialization
    !
    ! Kernel vectors
    g=(0.d0,0.d0)
    ! Calculate
    select case (fbem_n_dimension(type_g))
      !
      ! SURFACE LOAD
      !
      case (2)
        ! Number of nodes of geometrical interpolation
        nnodes_g=fbem_n_nodes(type_g)
        ! Calculate real coordinates of collocation point
        ! Geometrical shape functions at xi_i
#       define etype type_g
#       define delta 0.0d0
#       define xi xi_i
#       define phi phi_g
#       include <phi_2d.rc>
#       undef etype
#       undef delta
#       undef xi
#       undef phi
        ! Calculate x_i
        x_i=0.d0
        do kphi=1,nnodes_g
          x_i=x_i+phi_g(kphi)*x_nodes(:,kphi)
        end do
        ! Functional shape functions at xi_i
#       define etype type_f1
#       define delta delta_f
#       define xi xi_i
#       define phi phi_f1_i
#       include <phi_2d.rc>
#       undef etype
#       undef delta
#       undef xi
#       undef phi
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
#             define etype type_g
#             define delta 0.0d0
#             define phi phi_g
#             define dphidxi1 dphidxi1_g
#             define dphidxi2 dphidxi2_g
#             include <phi_and_dphidxik_2d.rc>
#             undef etype
#             undef delta
#             undef phi
#             undef dphidxi1
#             undef dphidxi2
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
              ! Unit normal vector
              n=N/jg
              ! Distance vector
              rv=x-x_i
              ! Distance vector norm
              r=dot_product(rv,rv)
              r=sqrt(r)
              dr1=1.d0/r
              dr2=dr1**2
              dr3=dr2*dr1
              dr4=dr3*dr1
              ! r_{,k}
              drdx=rv*dr1
              ! dr/dn
              drdn=dot_product(drdx,n)
              ! Jacobians * weights
              jw=jg*rho*jthetap*w_angular*w_radial
              ! FUNCTIONAL SHAPE FUNCTIONS
              ! Functional shape functions (secondary variables) at xi
#             define etype type_f2
#             define delta delta_f
#             define phi phi_f2
#             include <phi_2d.rc>
#             undef etype
#             undef delta
#             undef phi
              ! Functional shape functions * jacobians * weights
              phif2jw=phi_f2*jw
              ! COMPONENTS OF THE FUNDAMENTAL SOLUTION (ONLY REGULAR PARTS)
              z(1)=-c_im*p%k1*r
              z(2)=-c_im*p%k2*r
              call fbem_zexp_decomposed(2,z,E)
              E(2,:)=E(2,:)*dr1
              E(3,:)=E(3,:)*dr2
              E(4,:)=E(4,:)*dr3
              E(5,:)=E(5,:)*dr4
              psi=p%psi(1)*dr1+p%psi(2)+E(2,2)+p%psi(3)*E(3,1)+p%psi(4)*E(3,2)+p%psi(5)*E(4,1)+p%psi(6)*E(4,2)
              chi=p%chi(1)*dr1+p%chi(2)*E(2,1)+E(2,2)+p%chi(3)*E(3,1)+p%chi(4)*E(3,2)+p%chi(5)*E(4,1)+p%chi(6)*E(4,2)
              do il=1,3
                do ik=1,3
                  ! Regular parts of fundamental solutions values without constants
                  fs_u=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
                  g(:,il,ik)=g(:,il,ik)+fs_u*phif2jw
                end do
              end do
            end do ! Loop through rho coordinate
          end do ! Loop through theta coordinate
        end do ! Loop through rectangle triangles
        g=p%cte_u*g
      !
      ! VOLUME LOAD
      !
      case (3)
        stop 'not yet : fbem_bem_harela3d_sbie_bl_int'
      !
      ! OTHERS
      !
      case default

        write(*,*) xi_i
        write(*,*) x_nodes



        call fbem_error_message(output_unit,0,'fbem_bem_harela3d_sbie_bl_int',0,'it is only possible to integrate surface or volume loads')



    end select
  end subroutine fbem_bem_harela3d_sbie_bl_int

  subroutine fbem_bem_harela3d_sbie_bl_auto(e,x_i,p,qsp,ns,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                 !! Integration element
    real(kind=real64)                  :: x_i(3)            !! Collocation point
    type(fbem_bem_harela3d_parameters) :: p                 !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp               !! Quasi-singular integration parameters
    integer                            :: ns                !! Maximum level of subdivisions
    complex(kind=real64)               :: g(e%n_snodes,3,3) !! g integration kernel
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
        call fbem_bem_harela3d_sbie_bl_int(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,barxi,p,g)
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
          call fbem_bem_harela3d_sbie_bl_ext_pre(ps,e,x_i,p,g)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_harela3d_sbie_bl_ext_adp(e,xi_s,x_i,p,qsp,1,ns,g)
        end if
    end select
  end subroutine fbem_bem_harela3d_sbie_bl_auto
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (HBIE)

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! BOUNDARY ELEMENTS

  subroutine fbem_bem_harela3d_hbie_ext_pre(ps,e,reverse,x_i,n_i,p,m,l)
    implicit none
    ! I/O
    integer                            :: ps                !! Selected precalculated dataset
    type(fbem_bem_element)             :: e                 !! Element
    logical                            :: reverse           !! Reverse normal vector
    real(kind=real64)                  :: x_i(3)            !! Collocation point position vector
    real(kind=real64)                  :: n_i(3)            !! Collocation point unit normal vector
    type(fbem_bem_harela3d_parameters) :: p                 !! Parameters of the region
    complex(kind=real64)               :: m(e%n_pnodes,3,3) !! m integration kernels vector
    complex(kind=real64)               :: l(e%n_snodes,3,3) !! l integration kernels vector
    ! Local
    integer              :: kip                             ! Counter variable for integration points loop
    real(kind=real64)    :: x(3)                            ! Position vector at integration point
    real(kind=real64)    :: n(3)                            ! Unit normal vector at integration point
    real(kind=real64)    :: pphijw(e%n_pnodes)              ! phi^p * jacobian * weight at integration point
    real(kind=real64)    :: sphijw(e%n_snodes)              ! phi^s * jacobian * weight at integration point
    real(kind=real64)    :: rv(3)                           ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, d1r1, d1r2, d1r3, d1r4, d1r5 ! Distance vector module and its inverse
    real(kind=real64)    :: drdx(3)                         ! Distance vector derivatives with respect to x_k
    real(kind=real64)    :: drdn                            ! Partial derivative of r respect to unit normal
    real(kind=real64)    :: drdni                           ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)    :: n_dot_ni                        ! Dot product of n and n_i
    integer              :: il, ik                          ! Counter for load / observation components
    complex(kind=real64) :: z(2)                            ! Arguments z=ikr
    complex(kind=real64) :: EnR(0:6,2)                      ! Exponential function decomposition for each wavenumber
    complex(kind=real64) :: TT1, TT2, TT3                   ! Components of the fundamental solutions
    complex(kind=real64) :: S1, S2, S3, S4, S5              ! Components of the fundamental solutions
    complex(kind=real64) :: fs_d, fs_s                      ! Fundamental solutions values
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
      d1r4=d1r3*d1r1
      d1r5=d1r4*d1r1
      drdx=rv*d1r1
      drdn=dot_product(drdx,n)
      drdni=-dot_product(drdx,n_i)
      n_dot_ni=dot_product(n,n_i)
      z(1)=-c_im*p%k1*r
      z(2)=-c_im*p%k2*r
      call fbem_zexp_decomposed(2,z,EnR)
      EnR(2,:)=EnR(2,:)*d1r1
      EnR(3,:)=EnR(3,:)*d1r2
      EnR(4,:)=EnR(4,:)*d1r3
      EnR(5,:)=EnR(5,:)*d1r4
      EnR(6,:)=EnR(6,:)*d1r5
      TT1=p%T1(1)*d1r2+p%T1(2)+p%T1(3)*EnR(2,1)+p%T1(4)*EnR(2,2)+p%T1(5)*EnR(3,1)+p%T1(6)*EnR(3,2)+p%T1(7)*EnR(4,1)&
         +p%T1(8)*EnR(4,2)+p%T1(9)*EnR(5,1)+p%T1(10)*EnR(5,2)
      TT2=p%T2(1)*d1r2+p%T2(2)+p%T2(3)*EnR(2,2)+p%T2(4)*EnR(3,1)+p%T2(5)*EnR(3,2)+p%T2(6)*EnR(4,1)+p%T2(7)*EnR(4,2)&
         +p%T2(8)*EnR(5,1)+p%T2(9)*EnR(5,2)
      TT3=p%T3(1)*d1r2+p%T3(2)+p%T3(3)*EnR(2,1)+p%T3(4)*EnR(3,1)+p%T3(5)*EnR(3,2)+p%T3(6)*EnR(4,1)+p%T3(7)*EnR(4,2)&
         +p%T3(8)*EnR(5,1)+p%T3(9)*EnR(5,2)
      S1=p%S1(1)*d1r3+p%S1(2)*d1r1+p%S1(3)*EnR(2,2)+p%S1(4)*EnR(3,1)+p%S1(5)*EnR(3,2)+p%S1(6)*EnR(4,1)+p%S1(7)*EnR(4,2)&
        +p%S1(8)*EnR(5,1)+p%S1(9)*EnR(5,2)+p%S1(10)*EnR(6,1)+p%S1(11)*EnR(6,2)
      S2=p%S2(1)*d1r3+p%S2(2)*d1r1+p%S2(3)*EnR(2,1)+p%S2(4)*EnR(3,1)+p%S2(5)*EnR(3,2)+p%S2(6)*EnR(4,1)+p%S2(7)*EnR(4,2)&
        +p%S2(8)*EnR(5,1)+p%S2(9)*EnR(5,2)+p%S2(10)*EnR(6,1)+p%S2(11)*EnR(6,2)
      S3=p%S3(1)*d1r3+p%S3(2)*d1r1+p%S3(3)*EnR(2,1)+p%S3(4)*EnR(2,2)+p%S3(5)*EnR(3,1)+p%S3(6)*EnR(3,2)+p%S3(7)*EnR(4,1)&
        +p%S3(8)*EnR(4,2)+p%S3(9)*EnR(5,1)+p%S3(10)*EnR(5,2)+p%S3(11)*EnR(6,1)+p%S3(12)*EnR(6,2)
      S4=p%S4(1)*d1r3+p%S4(2)*d1r1+p%S4(3)+p%S4(4)*EnR(3,2)+p%S4(5)*EnR(4,1)+p%S4(6)*EnR(4,2)+p%S4(7)*EnR(5,1)+p%S4(8)*EnR(5,2)&
        +p%S4(9)*EnR(6,1)+p%S4(10)*EnR(6,2)
      S5=p%S5(1)*d1r3+p%S5(2)*d1r1+p%S5(3)+p%S5(4)*EnR(2,1)+p%S5(5)*EnR(3,1)+p%S5(6)*EnR(4,1)+p%S5(7)*EnR(4,2)+p%S5(8)*EnR(5,1)&
        +p%S5(9)*EnR(5,2)+p%S5(10)*EnR(6,1)+p%S5(11)*EnR(6,2)
      do il=1,3
        do ik=1,3
          fs_d=TT1*drdx(il)*drdx(ik)*drdni-TT2*(-drdni*c_dkr(il,ik)+drdx(il)*n_i(ik))-TT3*drdx(ik)*n_i(il)
          fs_s=S1*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni-c_dkr(il,ik)*drdn*drdni+drdx(il)*drdx(ik)*n_dot_ni)&
              +S2*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni)+S3*drdx(il)*drdx(ik)*drdn*drdni&
              +S4*(c_dkr(il,ik)*n_dot_ni+n_i(ik)*n(il))+S5*n(ik)*n_i(il)
          m(:,il,ik)=m(:,il,ik)+fs_s*pphijw
          l(:,il,ik)=l(:,il,ik)+fs_d*sphijw
        end do
      end do
    end do
    ! Multiply by constants
    m=p%cte_s*m
    l=p%cte_d*l
    ! Reverse if needed
    if (reverse) m=-m
  end subroutine fbem_bem_harela3d_hbie_ext_pre

  subroutine fbem_bem_harela3d_hbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,p,gln,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                                !! Integration element
    logical                            :: reverse                          !! Reverse normal vector
    real(kind=real64)                  :: xi_s(2,fbem_n_vertices(e%gtype)) !! Coordinates of the subdivision of the parent element
    real(kind=real64)                  :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)                  :: n_i(3)                           !! Unit normal at the collocation point
    real(kind=real64)                  :: barxip(2)                        !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)                  :: barr                             !! Telles jacobian at barxip
    type(fbem_bem_harela3d_parameters) :: p                                !! Parameters of the region
    integer                            :: gln                              !! 1D Gauss-Legendre number of integration points (<=32)
    complex(kind=real64)               :: m(e%n_pnodes,3,3)                !! m kernel vector
    complex(kind=real64)               :: l(e%n_snodes,3,3)                !! l kernel vector
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
    real(kind=real64)            :: r, d1r1, d1r2, d1r3, d1r4, d1r5 ! Distance vector module and its inverse
    real(kind=real64)            :: drdx(3)                         ! Distance vector derivatives with respect to x_k
    real(kind=real64)            :: jg                              ! Geometric jacobian
    real(kind=real64)            :: jw                              ! Jacobians * weights
    real(kind=real64)            :: drdn                            ! Partial derivative of r respect to unit normal
    real(kind=real64)            :: drdni                           ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)            :: n_dot_ni                        ! Dot product of n and n_i
    real(kind=real64)            :: pphijw(e%n_pnodes)              ! Functional shape functions * jw
    real(kind=real64)            :: sphijw(e%n_snodes)              ! Functional shape functions * jw
    integer                      :: il, ik                          ! Counter for load / observation components
    complex(kind=real64)         :: z(2)                            ! Arguments z=-ikr
    complex(kind=real64)         :: EnR(0:6,2)                      ! Exponential function decomposition for each wavenumber
    complex(kind=real64)         :: TT1, TT2, TT3                   ! Components of the fundamental solutions
    complex(kind=real64)         :: S1, S2, S3, S4, S5              ! Components of the fundamental solutions
    complex(kind=real64)         :: fs_d, fs_s                      ! Fundamental solutions values
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
            d1r4=d1r3*d1r1
            d1r5=d1r4*d1r1
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
            z(1)=-c_im*p%k1*r
            z(2)=-c_im*p%k2*r
            call fbem_zexp_decomposed(2,z,EnR)
            EnR(2,:)=EnR(2,:)*d1r1
            EnR(3,:)=EnR(3,:)*d1r2
            EnR(4,:)=EnR(4,:)*d1r3
            EnR(5,:)=EnR(5,:)*d1r4
            EnR(6,:)=EnR(6,:)*d1r5
            TT1=p%T1(1)*d1r2+p%T1(2)+p%T1(3)*EnR(2,1)+p%T1(4)*EnR(2,2)+p%T1(5)*EnR(3,1)+p%T1(6)*EnR(3,2)+p%T1(7)*EnR(4,1)&
               +p%T1(8)*EnR(4,2)+p%T1(9)*EnR(5,1)+p%T1(10)*EnR(5,2)
            TT2=p%T2(1)*d1r2+p%T2(2)+p%T2(3)*EnR(2,2)+p%T2(4)*EnR(3,1)+p%T2(5)*EnR(3,2)+p%T2(6)*EnR(4,1)+p%T2(7)*EnR(4,2)&
               +p%T2(8)*EnR(5,1)+p%T2(9)*EnR(5,2)
            TT3=p%T3(1)*d1r2+p%T3(2)+p%T3(3)*EnR(2,1)+p%T3(4)*EnR(3,1)+p%T3(5)*EnR(3,2)+p%T3(6)*EnR(4,1)+p%T3(7)*EnR(4,2)&
               +p%T3(8)*EnR(5,1)+p%T3(9)*EnR(5,2)
            S1=p%S1(1)*d1r3+p%S1(2)*d1r1+p%S1(3)*EnR(2,2)+p%S1(4)*EnR(3,1)+p%S1(5)*EnR(3,2)+p%S1(6)*EnR(4,1)+p%S1(7)*EnR(4,2)&
              +p%S1(8)*EnR(5,1)+p%S1(9)*EnR(5,2)+p%S1(10)*EnR(6,1)+p%S1(11)*EnR(6,2)
            S2=p%S2(1)*d1r3+p%S2(2)*d1r1+p%S2(3)*EnR(2,1)+p%S2(4)*EnR(3,1)+p%S2(5)*EnR(3,2)+p%S2(6)*EnR(4,1)+p%S2(7)*EnR(4,2)&
              +p%S2(8)*EnR(5,1)+p%S2(9)*EnR(5,2)+p%S2(10)*EnR(6,1)+p%S2(11)*EnR(6,2)
            S3=p%S3(1)*d1r3+p%S3(2)*d1r1+p%S3(3)*EnR(2,1)+p%S3(4)*EnR(2,2)+p%S3(5)*EnR(3,1)+p%S3(6)*EnR(3,2)+p%S3(7)*EnR(4,1)&
              +p%S3(8)*EnR(4,2)+p%S3(9)*EnR(5,1)+p%S3(10)*EnR(5,2)+p%S3(11)*EnR(6,1)+p%S3(12)*EnR(6,2)
            S4=p%S4(1)*d1r3+p%S4(2)*d1r1+p%S4(3)+p%S4(4)*EnR(3,2)+p%S4(5)*EnR(4,1)+p%S4(6)*EnR(4,2)+p%S4(7)*EnR(5,1)&
              +p%S4(8)*EnR(5,2)+p%S4(9)*EnR(6,1)+p%S4(10)*EnR(6,2)
            S5=p%S5(1)*d1r3+p%S5(2)*d1r1+p%S5(3)+p%S5(4)*EnR(2,1)+p%S5(5)*EnR(3,1)+p%S5(6)*EnR(4,1)+p%S5(7)*EnR(4,2)&
              +p%S5(8)*EnR(5,1)+p%S5(9)*EnR(5,2)+p%S5(10)*EnR(6,1)+p%S5(11)*EnR(6,2)
            ! Add
            do il=1,3
              do ik=1,3
                fs_d=TT1*drdx(il)*drdx(ik)*drdni-TT2*(drdx(il)*n_i(ik)-c_dkr(il,ik)*drdni)-TT3*drdx(ik)*n_i(il)
                fs_s=S1*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni-c_dkr(il,ik)*drdn*drdni+drdx(il)*drdx(ik)*n_dot_ni)&
                    +S2*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni)+S3*drdx(il)*drdx(ik)*drdn*drdni&
                    +S4*(c_dkr(il,ik)*n_dot_ni+n_i(ik)*n(il))+S5*n(ik)*n_i(il)
                m(:,il,ik)=m(:,il,ik)+fs_s*pphijw
                l(:,il,ik)=l(:,il,ik)+fs_d*sphijw
              end do
            end do
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
            d1r4=d1r3*d1r1
            d1r5=d1r4*d1r1
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
            z(1)=-c_im*p%k1*r
            z(2)=-c_im*p%k2*r
            call fbem_zexp_decomposed(2,z,EnR)
            EnR(2,:)=EnR(2,:)*d1r1
            EnR(3,:)=EnR(3,:)*d1r2
            EnR(4,:)=EnR(4,:)*d1r3
            EnR(5,:)=EnR(5,:)*d1r4
            EnR(6,:)=EnR(6,:)*d1r5
            TT1=p%T1(1)*d1r2+p%T1(2)+p%T1(3)*EnR(2,1)+p%T1(4)*EnR(2,2)+p%T1(5)*EnR(3,1)+p%T1(6)*EnR(3,2)+p%T1(7)*EnR(4,1)&
               +p%T1(8)*EnR(4,2)+p%T1(9)*EnR(5,1)+p%T1(10)*EnR(5,2)
            TT2=p%T2(1)*d1r2+p%T2(2)+p%T2(3)*EnR(2,2)+p%T2(4)*EnR(3,1)+p%T2(5)*EnR(3,2)+p%T2(6)*EnR(4,1)+p%T2(7)*EnR(4,2)&
               +p%T2(8)*EnR(5,1)+p%T2(9)*EnR(5,2)
            TT3=p%T3(1)*d1r2+p%T3(2)+p%T3(3)*EnR(2,1)+p%T3(4)*EnR(3,1)+p%T3(5)*EnR(3,2)+p%T3(6)*EnR(4,1)+p%T3(7)*EnR(4,2)&
               +p%T3(8)*EnR(5,1)+p%T3(9)*EnR(5,2)
            S1=p%S1(1)*d1r3+p%S1(2)*d1r1+p%S1(3)*EnR(2,2)+p%S1(4)*EnR(3,1)+p%S1(5)*EnR(3,2)+p%S1(6)*EnR(4,1)+p%S1(7)*EnR(4,2)&
              +p%S1(8)*EnR(5,1)+p%S1(9)*EnR(5,2)+p%S1(10)*EnR(6,1)+p%S1(11)*EnR(6,2)
            S2=p%S2(1)*d1r3+p%S2(2)*d1r1+p%S2(3)*EnR(2,1)+p%S2(4)*EnR(3,1)+p%S2(5)*EnR(3,2)+p%S2(6)*EnR(4,1)+p%S2(7)*EnR(4,2)&
              +p%S2(8)*EnR(5,1)+p%S2(9)*EnR(5,2)+p%S2(10)*EnR(6,1)+p%S2(11)*EnR(6,2)
            S3=p%S3(1)*d1r3+p%S3(2)*d1r1+p%S3(3)*EnR(2,1)+p%S3(4)*EnR(2,2)+p%S3(5)*EnR(3,1)+p%S3(6)*EnR(3,2)+p%S3(7)*EnR(4,1)&
              +p%S3(8)*EnR(4,2)+p%S3(9)*EnR(5,1)+p%S3(10)*EnR(5,2)+p%S3(11)*EnR(6,1)+p%S3(12)*EnR(6,2)
            S4=p%S4(1)*d1r3+p%S4(2)*d1r1+p%S4(3)+p%S4(4)*EnR(3,2)+p%S4(5)*EnR(4,1)+p%S4(6)*EnR(4,2)+p%S4(7)*EnR(5,1)&
              +p%S4(8)*EnR(5,2)+p%S4(9)*EnR(6,1)+p%S4(10)*EnR(6,2)
            S5=p%S5(1)*d1r3+p%S5(2)*d1r1+p%S5(3)+p%S5(4)*EnR(2,1)+p%S5(5)*EnR(3,1)+p%S5(6)*EnR(4,1)+p%S5(7)*EnR(4,2)&
              +p%S5(8)*EnR(5,1)+p%S5(9)*EnR(5,2)+p%S5(10)*EnR(6,1)+p%S5(11)*EnR(6,2)
            ! Add
            do il=1,3
              do ik=1,3
                fs_d=TT1*drdx(il)*drdx(ik)*drdni-TT2*(drdx(il)*n_i(ik)-c_dkr(il,ik)*drdni)-TT3*drdx(ik)*n_i(il)
                fs_s=S1*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni-c_dkr(il,ik)*drdn*drdni+drdx(il)*drdx(ik)*n_dot_ni)&
                    +S2*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni)+S3*drdx(il)*drdx(ik)*drdn*drdni&
                    +S4*(c_dkr(il,ik)*n_dot_ni+n_i(ik)*n(il))+S5*n(ik)*n_i(il)
                m(:,il,ik)=m(:,il,ik)+fs_s*pphijw
                l(:,il,ik)=l(:,il,ik)+fs_d*sphijw
              end do
            end do
          end do
        end do
        case default
          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'n_edges not valid')
    end select
    ! Multiply by constants
    m=p%cte_s*m
    l=p%cte_d*l
    ! Reverse if needed
    if (reverse) m=-m
  end subroutine fbem_bem_harela3d_hbie_ext_st

  recursive subroutine fbem_bem_harela3d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,p,qsp,ks,ns,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                                !! Element
    logical                            :: reverse                          !! Reverse orientation
    real(kind=real64)                  :: xi_s(2,fbem_n_vertices(e%gtype)) !! Subdivision of the parent element
    real(kind=real64)                  :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)                  :: n_i(3)                           !! Unit normal at the collocation point
    type(fbem_bem_harela3d_parameters) :: p                                !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp                              !! Quasi-singular integration parameters
    integer                            :: ks                               !! Current level of subdivisions
    integer                            :: ns                               !! Maximum level of subdivision
    complex(kind=real64)               :: m(e%n_pnodes,3,3)                !! m integration kernels matrix
    complex(kind=real64)               :: l(e%n_snodes,3,3)                !! l integration kernels matrix
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
    complex(kind=real64) :: m_tmp(e%n_pnodes,3,3)                ! m integration kernels matrix (temporary)
    complex(kind=real64) :: l_tmp(e%n_snodes,3,3)                ! l integration kernels matrix (temporary)
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
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,7,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_harela3d_hbie_ext_adp',ns,'maximum number of subdivisions reached')
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
          call fbem_bem_harela3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
          ! SUBTRI 2
          tmp_xi_s(:,1)=xi_s(:,2)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          call fbem_bem_harela3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
          ! SUBTRI 3
          tmp_xi_s(:,1)=xi_s(:,3)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          call fbem_bem_harela3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
          ! SUBTRI 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          call fbem_bem_harela3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
        ! QUADRILATERALS
        case (4)
          ! SUBQUAD 1
          tmp_xi_s(:,1)=xi_s(:,1)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,3)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          call fbem_bem_harela3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
          ! SUBQUAD 2
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=xi_s(:,2)
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,4)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          call fbem_bem_harela3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
          ! SUBQUAD 3
          tmp_xi_s(:,1)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=xi_s(:,3)
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          call fbem_bem_harela3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
          ! SUBQUAD 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          tmp_xi_s(:,2)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=xi_s(:,4)
          call fbem_bem_harela3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
      end select
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_harela3d_hbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,p,gln,m_tmp,l_tmp)
      m=m+m_tmp
      l=l+l_tmp
    end if
  end subroutine fbem_bem_harela3d_hbie_ext_adp

  subroutine fbem_bem_harela3d_hbie_int(type_g,type_f1,type_f2,delta_f,x_nodes,reverse,xi_i,p,m,l)
    implicit none
    ! I/O
    integer                            :: type_g                            !! Geometrial interpolation
    integer                            :: type_f1                           !! Functional interpolation (primary variables)
    integer                            :: type_f2                           !! Functional interpolation (secondary variables)
    real(kind=real64)                  :: delta_f                           !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64)                  :: x_nodes(3,fbem_n_nodes(type_g))   !! Position vectors of geometrical nodes
    logical                            :: reverse                           !! Normal vector inversion
    real(kind=real64)                  :: xi_i(2)                           !! Reference coordinates of the singular point.
    type(fbem_bem_harela3d_parameters) :: p                                 !! Parameters of the region
    complex(kind=real64)               :: m(fbem_n_nodes(type_f1),3,3)      !! m kernel vector
    complex(kind=real64)               :: l(fbem_n_nodes(type_f2),3,3)      !! l kernel vector
    ! Local
    integer                            :: ksubtri                           ! Counter variable for subtriangles loop
    integer                            :: nsubtri                           ! Number of subtriangles
    integer                            :: subtriangle(8)                    ! Vector that contains what subtriangles need to be integrated
    real(kind=real64)                  :: theta_subtri(2,8)                 ! Matrix that contains the angles theta of the subtriangles to be integrated
    real(kind=real64)                  :: thetap_subtri(2,8)                ! Matrix that contains the angles thetap of the subtriangles to be integrated
    integer                            :: ktheta                            ! Counter variable for theta coordinate loop
    integer                            :: krho                              ! Counter variable for rho coordinate loop
    integer                            :: kphi                              ! Counter coordinates for shape functions loops
    integer                            :: nnodes_g                          ! Number of nodes of geometrical interpolation
    integer                            :: nnodes_f1                         ! Number of nodes of primary variables interpolation
    integer                            :: ngp_theta                         ! Number of Gauss points for theta coordinate
    integer                            :: ngp_rho                           ! Number of Gauss points for rho coordinate
    real(kind=real64)                  :: thetai, thetaf, thetapi, thetapf  ! Initial and final angles for subtriangle integration
    real(kind=real64)                  :: w_angular                         ! Weight of the angular coordinate
    real(kind=real64)                  :: w_radial                          ! Weight of the radial coordinate
    real(kind=real64)                  :: theta                             ! Angle coordinate theta
    real(kind=real64)                  :: thetap                            ! Angle coordinate thetap
    real(kind=real64)                  :: thetapp                           ! Angle coordinate thetap on [0,1] domain
    real(kind=real64)                  :: jthetap                           ! thetap->thetapp jacobian
    real(kind=real64)                  :: rhoij                             ! Maximum rho (radial) value for each edge
    real(kind=real64)                  :: rho                               ! Radial coordinate rho
    real(kind=real64)                  :: rhop                              ! Radial coordinate rho on [0,1] domain
    real(kind=real64)                  :: aux(10)                           ! Auxiliary variable for shape functions resources
    real(kind=real64)                  :: xi(2)                             ! Reference xi_1,xi_2 coordinates
    real(kind=real64)                  :: phi_f1(fbem_n_nodes(type_f1))     ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)                  :: phi_f2(fbem_n_nodes(type_f2))     ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)                  :: phi_f1_i(fbem_n_nodes(type_f1))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)                  :: psi_i(3,fbem_n_nodes(type_f1))    ! psi vectors at xi_i1,xi_i2
    real(kind=real64)                  :: phi_f2_i(fbem_n_nodes(type_f2))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)                  :: phi_g(fbem_n_nodes(type_g))       ! Geometrical shape functions values at xi_1,xi_2
    real(kind=real64)                  :: dphidxi1_g(fbem_n_nodes(type_g))  ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64)                  :: dphidxi2_g(fbem_n_nodes(type_g))  ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64)                  :: jxi1, jxi2                        ! xi1->x, xi2->x jacobians
    real(kind=real64)                  :: x_i(3)                            ! Collocation point position vector
    real(kind=real64)                  :: n_i(3)                            ! Collocation point unit normal vector
    real(kind=real64)                  :: x(3)                              ! Position vector at xi_1,xi_2
    real(kind=real64)                  :: T1(3), T2(3)                      ! Tangent vectors at xi_1,xi_2
    real(kind=real64)                  :: v1(3), v2(3)                      ! Orthogonal tangent vectors
    real(kind=real64)                  :: Mvt(2,2)                          ! Orthogonal coordinates transformation matrix
    real(kind=real64)                  :: dxidv(2,2)                        ! xi coordinates derivatives with respect to v orthogonal coordinates
    real(kind=real64)                  :: detMvt                            ! Mvt(1,1)*Mvt(2,2)-Mvt(1,2)*Mvt(2,1)
    real(kind=real64)                  :: N(3)                              ! Normal vector at xi_1,xi_2
    real(kind=real64)                  :: rv(3)                             ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)                  :: r, dr1, dr2, dr3, dr4, dr5        ! Distance vector module and its inverse
    real(kind=real64)                  :: drdx(3)                           ! Distance vector derivatives with respect to x_k
    real(kind=real64)                  :: jg                                ! Geometric jacobian
    real(kind=real64)                  :: jw                                ! Jacobian * weights
    real(kind=real64)                  :: drdn                              ! Partial derivative of r respect to unit normal
    real(kind=real64)                  :: drdni                             ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)                  :: n_dot_ni                          ! Dot product n · n_i
    real(kind=real64)                  :: phif1jw(fbem_n_nodes(type_f1))    ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)                  :: phif2jw(fbem_n_nodes(type_f2))    ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)                  :: costheta, sintheta                ! cos(theta), sin(theta)
    integer                            :: il, ik, ij                        ! Counter for load / observation components
    complex(kind=real64)               :: z(2)                              ! Wavenumbers
    complex(kind=real64)               :: E(0:6,2)                          ! Exponential function decomposition for each wavenumber
    complex(kind=real64)               :: TT1, TT2, TT3, S1, S2, S3, S4, S5 ! Components of the fundamental solutions
    complex(kind=real64)               :: fs_d, fs_s                        ! Fundamental solutions values
    !
    real(kind=real64)                  :: re                                ! Quasi-singular integration relative error
    integer                            :: gln_min                           ! Minimum Number of 1D Gauss-Legendre number of points
    integer                            :: ns_max                            ! Maximum number of subdivisions
    type(fbem_qs_parameters)           :: qsp                               ! Quasi-singular integration parameters
    real(kind=real64)                  :: xi_s(1,2)                         ! Edge subdivision
    integer                            :: kedge                             ! Counter variable for edges
    integer                            :: nedges                            ! Number of edges
    integer                            :: type_edge                         ! Line element type for line integrals
    integer                            :: nnodes_edge                       ! Number of nodes of the edge
    real(kind=real64), allocatable     :: x_nodes_edge(:,:)                 ! Coordinates of the edge elements
    real(kind=real64)                  :: lli(3,3)                          ! Line integral epsilon_{ijk}·e_k·t/r
    real(kind=real64)                  :: mli1                              ! Line integral (r x ni)·t/r**3
    real(kind=real64)                  :: mli2(3)                           ! Line integral (e_k x ni)·t/r
    real(kind=real64)                  :: mli3(3,3)                         ! Line integral (r x ni)·t*r_{,l}*r_{,k}/r**3
    real(kind=real64)                  :: mli4(3)                           ! Line integral (r x e_k)·t/r**3
    real(kind=real64)                  :: mli5(3,3,3)                       ! Line integral (e_k x ni)·t*r_{,l}*r_{,j}/r
    !
    ! Initialization
    !
    ! Kernels
    m=(0.d0,0.d0)
    l=(0.d0,0.d0)
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
    ! Functional shape functions (secondary variables) at xi_i
#   define etype type_f2
#   define delta delta_f
#   define xi xi_i
#   define phi phi_f2_i
#   include <phi_2d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
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
          ! Unit normal vector
          n=N/jg
          ! Distance vector
          rv=x-x_i
          ! Distance vector norm
          r=dot_product(rv,rv)
          r=sqrt(r)
          dr1=1.d0/r
          dr2=dr1**2
          dr3=dr2*dr1
          dr4=dr3*dr1
          dr5=dr4*dr1
          ! r_{,k}
          drdx=rv*dr1
          ! dr/dn
          drdn=dot_product(drdx,n)
          ! dr/dni
          drdni=-dot_product(drdx,n_i)
          ! n·ni
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
          ! COMPONENTS OF THE FUNDAMENTAL SOLUTION (ONLY REGULAR PARTS)
          z(1)=-c_im*p%k1*r
          z(2)=-c_im*p%k2*r
          call fbem_zexp_decomposed(2,z,E)
          E(2,:)=E(2,:)*dr1
          E(3,:)=E(3,:)*dr2
          E(4,:)=E(4,:)*dr3
          E(5,:)=E(5,:)*dr4
          E(6,:)=E(6,:)*dr5
          TT1=p%T1(2)+p%T1(3)*E(2,1)+p%T1(4)*E(2,2)+p%T1(5)*E(3,1)+p%T1(6)*E(3,2)+p%T1(7)*E(4,1)+p%T1(8)*E(4,2)&
             +p%T1(9)*E(5,1)+p%T1(10)*E(5,2)
          TT2=p%T2(2)+p%T2(3)*E(2,2)+p%T2(4)*E(3,1)+p%T2(5)*E(3,2)+p%T2(6)*E(4,1)+p%T2(7)*E(4,2)+p%T2(8)*E(5,1)&
             +p%T2(9)*E(5,2)
          TT3=p%T3(2)+p%T3(3)*E(2,1)+p%T3(4)*E(3,1)+p%T3(5)*E(3,2)+p%T3(6)*E(4,1)+p%T3(7)*E(4,2)+p%T3(8)*E(5,1)&
             +p%T3(9)*E(5,2)
          S1=p%S1(2)*dr1+p%S1(3)*E(2,2)+p%S1(4)*E(3,1)+p%S1(5)*E(3,2)+p%S1(6)*E(4,1)+p%S1(7)*E(4,2)&
            +p%S1(8)*E(5,1)+p%S1(9)*E(5,2)+p%S1(10)*E(6,1)+p%S1(11)*E(6,2)
          S2=p%S2(2)*dr1+p%S2(3)*E(2,1)+p%S2(4)*E(3,1)+p%S2(5)*E(3,2)+p%S2(6)*E(4,1)+p%S2(7)*E(4,2)&
            +p%S2(8)*E(5,1)+p%S2(9)*E(5,2)+p%S2(10)*E(6,1)+p%S2(11)*E(6,2)
          S3=p%S3(2)*dr1+p%S3(3)*E(2,1)+p%S3(4)*E(2,2)+p%S3(5)*E(3,1)+p%S3(6)*E(3,2)+p%S3(7)*E(4,1)&
            +p%S3(8)*E(4,2)+p%S3(9)*E(5,1)+p%S3(10)*E(5,2)+p%S3(11)*E(6,1)+p%S3(12)*E(6,2)
          S4=p%S4(2)*dr1+p%S4(3)+p%S4(4)*E(3,2)+p%S4(5)*E(4,1)+p%S4(6)*E(4,2)+p%S4(7)*E(5,1)+p%S4(8)*E(5,2)&
            +p%S4(9)*E(6,1)+p%S4(10)*E(6,2)
          S5=p%S5(2)*dr1+p%S5(3)+p%S5(4)*E(2,1)+p%S5(5)*E(3,1)+p%S5(6)*E(4,1)+p%S5(7)*E(4,2)+p%S5(8)*E(5,1)&
            +p%S5(9)*E(5,2)+p%S5(10)*E(6,1)+p%S5(11)*E(6,2)
          ! Loop through load direction and observation direction
          do il=1,3
            do ik=1,3
              !
              ! Regular parts
              !
              ! Regular parts of fundamental solutions values without constants
              fs_d=(TT1+p%T1(1)*dr2)*drdx(il)*drdx(ik)*drdni+(TT2+p%T2(1)*dr2)*c_dkr(il,ik)*drdni-TT2*drdx(il)*n_i(ik)&
                  -TT3*drdx(ik)*n_i(il)+p%T2(1)*dr2*((n_i(il)-n(il))*drdx(ik)-(n_i(ik)-n(ik))*drdx(il))
              fs_s=S1*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni-c_dkr(il,ik)*drdn*drdni+drdx(il)*drdx(ik)*n_dot_ni)&
                  +S2*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni)+S3*drdx(il)*drdx(ik)*drdn*drdni&
                  +S4*(c_dkr(il,ik)*n_dot_ni+n_i(ik)*n(il))+S5*n(ik)*n_i(il)-c_dkr(il,ik)*p%S1(1)*dr3*drdn*drdni&
                  +p%S3(1)*dr3*drdx(il)*drdx(ik)*drdn*drdni
              ! Add to kernels
              m(:,il,ik)=m(:,il,ik)+fs_s*phif1jw
              l(:,il,ik)=l(:,il,ik)+fs_d*phif2jw
              !
              ! Regular integrals after regularization (L)
              !
              ! L21
              if (il.ne.ik) then
                fs_d=p%T2(1)*dr2*(n(il)*drdx(ik)-n(ik)*drdx(il))
                l(:,il,ik)=l(:,il,ik)+fs_d*(phi_f2-phi_f2_i)*jw
              end if
              !
              ! Regular integrals after regularization (M)
              !
              ! Regular parts of M2
              if (il.eq.ik) then
                ! M21
                fs_s=p%S4(1)*dr3*n_dot_ni
                m(:,il,ik)=m(:,il,ik)+fs_s*(phi_f1-phi_f1_i-psi_i(1,:)*rv(1)-psi_i(2,:)*rv(2)-psi_i(3,:)*rv(3))*jw
                ! M221
                fs_s=-p%S4(1)*3.d0*dr3*drdn*drdni
                m(:,il,ik)=m(:,il,ik)+phi_f1_i*fs_s*jw
                ! M231
                fs_s=-p%S4(1)*dr2*drdni
                m(:,il,ik)=m(:,il,ik)+fs_s*(psi_i(1,:)*n(1)+psi_i(2,:)*n(2)+psi_i(3,:)*n(3))*jw
              end if
              ! Regular parts of M3
              ! M31
              fs_s=dr3*(p%S2(1)*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni)&
                        +p%S1(1)*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni)&
                        +p%S1(1)*drdx(il)*drdx(ik)*n_dot_ni&
                        +p%S4(1)*n_i(ik)*n(il)+p%S5(1)*n_i(il)*n(ik))
              m(:,il,ik)=m(:,il,ik)+fs_s*(phi_f1-phi_f1_i-psi_i(1,:)*rv(1)-psi_i(2,:)*rv(2)-psi_i(3,:)*rv(3))*jw
              ! M321
              fs_s=dr3*(-5.d0*p%S1(1)*drdx(il)*drdx(ik)*drdn*drdni+p%S1(1)*drdx(ik)*(n_i(il)*drdn-n(il)*drdni)&
                         +p%S2(1)*drdx(il)*(n_i(ik)*drdn-n(ik)*drdni))
              m(:,il,ik)=m(:,il,ik)+phi_f1_i*fs_s*jw
              ! M331
              do ij=1,3
                fs_s=dr2*(drdx(ij)*(p%S2(1)*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni)&
                                    +p%S1(1)*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni)&
                                    -p%S1(1)*drdx(il)*n(ik)*drdni)&
                          -((p%S1(1)/3.d0-p%S5(1))*n_i(il)*drdx(ij)+p%S1(1)/3.d0*n_i(ij)*drdx(il))*(n(ik)-n_i(ik)*n_dot_ni)&
                          +p%S4(1)*n_i(ik)*drdx(ij)*(n(il)-n_i(il)*n_dot_ni)&
                          -p%S1(1)/3.d0*(c_dkr(il,ik)+n_i(il)*n_i(ik))*n(ij)*drdni&
                          -p%S1(1)/3.d0*(c_dkr(ij,ik)-n_i(ij)*n_i(ik))*n(il)*drdni)
                m(:,il,ik)=m(:,il,ik)+psi_i(ij,:)*fs_s*jw
              end do
            end do
          end do
        end do ! Loop through rho coordinate
      end do ! Loop through theta coordinate
    end do ! Loop through rectangle triangles
    !
    ! LINE INTEGRALS (all quasi-singular integrals, all integrations performed with an adaptive Telles+subdivision)
    !
    ! Initialize
    lli=0.d0
    mli1=0.d0
    mli2=0.d0
    mli3=0.d0
    mli4=0.d0
    mli5=0.d0
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
      ! Calculate all line integrals
      call fbem_bem_staela3d_hbie_int_li(type_edge,x_nodes_edge,xi_s,x_i,n_i,gln_min,qsp,1,ns_max,lli,mli1,mli2,mli3,mli4,mli5)
      ! Deallocate the edge coordinates
      deallocate (x_nodes_edge)
    end do ! Loop through edges
    !
    ! Add line integrals
    !
    ! Loop through load direction and observation direction
    do il=1,3
      do ik=1,3
        ! L22
        l(:,il,ik)=l(:,il,ik)+phi_f2_i*p%T2(1)*lli(il,ik)
        if (il.eq.ik) then
          ! M222
          m(:,il,ik)=m(:,il,ik)+phi_f1_i*p%S4(1)*mli1
          ! M232
          do ij=1,3
            m(:,il,ik)=m(:,il,ik)+psi_i(ij,:)*p%S4(1)*mli2(ij)
          end do
        end if
        ! M322
        m(:,il,ik)=m(:,il,ik)+phi_f1_i*p%S1(1)*mli3(il,ik)
        ! M323
        m(:,il,ik)=m(:,il,ik)+phi_f1_i*p%S4(1)*n_i(ik)*mli4(il)
        ! M324
        m(:,il,ik)=m(:,il,ik)+phi_f1_i*p%S5(1)*n_i(il)*mli4(ik)
        ! M33X
        do ij=1,3
          ! M332
          m(:,il,ik)=m(:,il,ik)+psi_i(ij,:)*p%S1(1)/3.d0*mli5(ij,il,ik)
          ! M333
          m(:,il,ik)=m(:,il,ik)+psi_i(ij,:)*p%S1(1)/3.d0*(c_dkr(il,ik)+n_i(il)*n_i(ik))*mli2(ij)
          ! M334
          m(:,il,ik)=m(:,il,ik)+psi_i(ij,:)*p%S1(1)/3.d0*(c_dkr(ij,ik)-n_i(ij)*n_i(ik))*mli2(il)
        end do
      end do
    end do
    ! Multiply m and l by constants of s* and d* respectively
    m=p%cte_s*m
    l=p%cte_d*l
    ! If the normal has to be reversed, then l=-l
    if (reverse) l=-l
  end subroutine fbem_bem_harela3d_hbie_int

  subroutine fbem_bem_harela3d_hbie_auto(e,reverse,x_i,n_i,p,qsp,ns,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                 !! Integration element
    logical                            :: reverse           !! Reverse orientation
    real(kind=real64)                  :: x_i(3)            !! Collocation point
    real(kind=real64)                  :: n_i(3)            !! Unit normal at the collocation point
    type(fbem_bem_harela3d_parameters) :: p                 !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp               !! Quasi-singular integration parameters
    integer                            :: ns                !! Maximum level of subdivisions
    complex(kind=real64)               :: m(e%n_pnodes,3,3) !! m integration kernel
    complex(kind=real64)               :: l(e%n_snodes,3,3) !! l integration kernel
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
        call fbem_bem_harela3d_hbie_int(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi,p,m,l)
      case (0)
        ! Estimate the required integration rule
        gln_near=fbem_qs_n_estimation_standard(e%n,e%gtype,7,qsp,d,barxi)
        gln=max(e%gln_far,gln_near)
        ! Integrate using a conservative precalculated dataset
        if ((gln.le.e%ps_gln_max).and.(gln_near.gt.0)) then
          do i=1,e%n_ps
            if (e%ps_gln(i).ge.gln) then
              ps=i
              exit
            end if
          end do
          call fbem_bem_harela3d_hbie_ext_pre(ps,e,reverse,x_i,n_i,p,m,l)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_harela3d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,p,qsp,1,ns,m,l)
        end if
    end select
  end subroutine fbem_bem_harela3d_hbie_auto
  ! --------------------------------------------------------------------------------------------------------------------------------

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! BODY LOAD ELEMENTS

  subroutine fbem_bem_harela3d_hbie_bl_ext_pre(ps,e,x_i,n_i,p,l)
    implicit none
    ! I/O
    integer                            :: ps                !! Selected precalculated dataset
    type(fbem_bem_element)             :: e                 !! Element
    real(kind=real64)                  :: x_i(3)            !! Collocation point position vector
    real(kind=real64)                  :: n_i(3)            !! Collocation point unit normal vector
    type(fbem_bem_harela3d_parameters) :: p                 !! Parameters of the region
    complex(kind=real64)               :: l(e%n_snodes,3,3) !! l integration kernels vector
    ! Local
    integer              :: kip                             ! Counter variable for integration points loop
    real(kind=real64)    :: x(3)                            ! Position vector at integration point
    real(kind=real64)    :: n(3)                            ! Unit normal vector at integration point
    real(kind=real64)    :: pphijw(e%n_pnodes)              ! phi^p * jacobian * weight at integration point
    real(kind=real64)    :: sphijw(e%n_snodes)              ! phi^s * jacobian * weight at integration point
    real(kind=real64)    :: rv(3)                           ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, d1r1, d1r2, d1r3, d1r4, d1r5 ! Distance vector module and its inverse
    real(kind=real64)    :: drdx(3)                         ! Distance vector derivatives with respect to x_k
    real(kind=real64)    :: drdn                            ! Partial derivative of r respect to unit normal
    real(kind=real64)    :: drdni                           ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)    :: n_dot_ni                        ! Dot product of n and n_i
    integer              :: il, ik                          ! Counter for load / observation components
    complex(kind=real64) :: z(2)                            ! Arguments z=ikr
    complex(kind=real64) :: EnR(0:6,2)                      ! Exponential function decomposition for each wavenumber
    complex(kind=real64) :: TT1, TT2, TT3                   ! Components of the fundamental solutions
    complex(kind=real64) :: fs_d                            ! Fundamental solutions values
    ! Initialize
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
      d1r4=d1r3*d1r1
      d1r5=d1r4*d1r1
      drdx=rv*d1r1
      drdn=dot_product(drdx,n)
      drdni=-dot_product(drdx,n_i)
      n_dot_ni=dot_product(n,n_i)
      z(1)=-c_im*p%k1*r
      z(2)=-c_im*p%k2*r
      call fbem_zexp_decomposed(2,z,EnR)
      EnR(2,:)=EnR(2,:)*d1r1
      EnR(3,:)=EnR(3,:)*d1r2
      EnR(4,:)=EnR(4,:)*d1r3
      EnR(5,:)=EnR(5,:)*d1r4
      EnR(6,:)=EnR(6,:)*d1r5
      TT1=p%T1(1)*d1r2+p%T1(2)+p%T1(3)*EnR(2,1)+p%T1(4)*EnR(2,2)+p%T1(5)*EnR(3,1)+p%T1(6)*EnR(3,2)+p%T1(7)*EnR(4,1)&
         +p%T1(8)*EnR(4,2)+p%T1(9)*EnR(5,1)+p%T1(10)*EnR(5,2)
      TT2=p%T2(1)*d1r2+p%T2(2)+p%T2(3)*EnR(2,2)+p%T2(4)*EnR(3,1)+p%T2(5)*EnR(3,2)+p%T2(6)*EnR(4,1)+p%T2(7)*EnR(4,2)&
         +p%T2(8)*EnR(5,1)+p%T2(9)*EnR(5,2)
      TT3=p%T3(1)*d1r2+p%T3(2)+p%T3(3)*EnR(2,1)+p%T3(4)*EnR(3,1)+p%T3(5)*EnR(3,2)+p%T3(6)*EnR(4,1)+p%T3(7)*EnR(4,2)&
         +p%T3(8)*EnR(5,1)+p%T3(9)*EnR(5,2)
      do il=1,3
        do ik=1,3
          fs_d=TT1*drdx(il)*drdx(ik)*drdni-TT2*(-drdni*c_dkr(il,ik)+drdx(il)*n_i(ik))-TT3*drdx(ik)*n_i(il)
          l(:,il,ik)=l(:,il,ik)+fs_d*sphijw
        end do
      end do
    end do
    l=p%cte_d*l
  end subroutine fbem_bem_harela3d_hbie_bl_ext_pre

  subroutine fbem_bem_harela3d_hbie_bl_ext_st(e,xi_s,x_i,n_i,barxip,barr,p,gln,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                                !! Integration element
    real(kind=real64)                  :: xi_s(2,fbem_n_vertices(e%gtype)) !! Coordinates of the subdivision of the parent element
    real(kind=real64)                  :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)                  :: n_i(3)                           !! Unit normal at the collocation point
    real(kind=real64)                  :: barxip(2)                        !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)                  :: barr                             !! Telles jacobian at barxip
    type(fbem_bem_harela3d_parameters) :: p                                !! Parameters of the region
    integer                            :: gln                              !! 1D Gauss-Legendre number of integration points (<=32)
    complex(kind=real64)               :: l(e%n_snodes,3,3)                !! l kernel vector
    ! Local
    integer                      :: kphi                            ! Counter variable for shape functions loops
    integer                      :: k1                              ! Counter variable for reference coordinate xi_1
    integer                      :: k2                              ! Counter variable for reference coordinate xi_2
    real(kind=real64)            :: aux(10)                         ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)                ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi1(e%n_gnodes)           ! Geometrical shape functions derivatives values
    real(kind=real64)            :: dgphidxi2(e%n_gnodes)           ! Geometrical shape functions derivatives values
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
    real(kind=real64)            :: r, d1r1, d1r2, d1r3, d1r4, d1r5 ! Distance vector module and its inverse
    real(kind=real64)            :: drdx(3)                         ! Distance vector derivatives with respect to x_k
    real(kind=real64)            :: jg                              ! Geometric jacobian
    real(kind=real64)            :: jw                              ! Jacobians * weights
    real(kind=real64)            :: drdn                            ! Partial derivative of r respect to unit normal
    real(kind=real64)            :: drdni                           ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)            :: n_dot_ni                        ! Dot product of n and n_i
    real(kind=real64)            :: sphijw(e%n_snodes)              ! Functional shape functions * jw
    integer                      :: il, ik                          ! Counter for load / observation components
    complex(kind=real64)         :: z(2)                            ! Arguments z=-ikr
    complex(kind=real64)         :: EnR(0:6,2)                      ! Exponential function decomposition for each wavenumber
    complex(kind=real64)         :: TT1, TT2, TT3                   ! Components of the fundamental solutions
    complex(kind=real64)         :: fs_d                            ! Fundamental solutions values
    ! Initialization
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
            d1r4=d1r3*d1r1
            d1r5=d1r4*d1r1
            drdx=rv*d1r1
            drdn=dot_product(drdx,n)
            drdni=-dot_product(drdx,n_i)
            n_dot_ni=dot_product(n,n_i)
            ! Jacobians * weights
            jw=jg*js*jt(1)*jt(2)*w(1)*w(2)
            ! FUNCTIONAL SHAPE FUNCTIONS
            ! Functional shape functions (secondary variables) at xi
#           define etype e%stype
#           define delta e%stype_delta
#           define phi sphi
#           include <phi_2d.rc>
#           undef etype
#           undef delta
#           undef phi
            ! Functional shape functions * jacobians * weights
            sphijw=sphi*jw
            ! COMPONENTS OF THE FUNDAMENTAL SOLUTION
            z(1)=-c_im*p%k1*r
            z(2)=-c_im*p%k2*r
            call fbem_zexp_decomposed(2,z,EnR)
            EnR(2,:)=EnR(2,:)*d1r1
            EnR(3,:)=EnR(3,:)*d1r2
            EnR(4,:)=EnR(4,:)*d1r3
            EnR(5,:)=EnR(5,:)*d1r4
            EnR(6,:)=EnR(6,:)*d1r5
            TT1=p%T1(1)*d1r2+p%T1(2)+p%T1(3)*EnR(2,1)+p%T1(4)*EnR(2,2)+p%T1(5)*EnR(3,1)+p%T1(6)*EnR(3,2)+p%T1(7)*EnR(4,1)&
               +p%T1(8)*EnR(4,2)+p%T1(9)*EnR(5,1)+p%T1(10)*EnR(5,2)
            TT2=p%T2(1)*d1r2+p%T2(2)+p%T2(3)*EnR(2,2)+p%T2(4)*EnR(3,1)+p%T2(5)*EnR(3,2)+p%T2(6)*EnR(4,1)+p%T2(7)*EnR(4,2)&
               +p%T2(8)*EnR(5,1)+p%T2(9)*EnR(5,2)
            TT3=p%T3(1)*d1r2+p%T3(2)+p%T3(3)*EnR(2,1)+p%T3(4)*EnR(3,1)+p%T3(5)*EnR(3,2)+p%T3(6)*EnR(4,1)+p%T3(7)*EnR(4,2)&
               +p%T3(8)*EnR(5,1)+p%T3(9)*EnR(5,2)
            ! Add
            do il=1,3
              do ik=1,3
                fs_d=TT1*drdx(il)*drdx(ik)*drdni-TT2*(drdx(il)*n_i(ik)-c_dkr(il,ik)*drdni)-TT3*drdx(ik)*n_i(il)
                l(:,il,ik)=l(:,il,ik)+fs_d*sphijw
              end do
            end do
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
            d1r4=d1r3*d1r1
            d1r5=d1r4*d1r1
            drdx=rv*d1r1
            drdn=dot_product(drdx,n)
            drdni=-dot_product(drdx,n_i)
            n_dot_ni=dot_product(n,n_i)
            ! Jacobians * weights
            jw=jg*js*jqt*jt(1)*jt(2)*w(1)*w(2)
            ! FUNCTIONAL SHAPE FUNCTIONS
            ! Functional shape functions (secondary variables) at xi
#           define etype e%stype
#           define delta e%stype_delta
#           define phi sphi
#           include <phi_2d.rc>
#           undef etype
#           undef delta
#           undef phi
            ! Functional shape functions * jacobians * weights
            sphijw=sphi*jw
            ! COMPONENTS OF THE FUNDAMENTAL SOLUTION
            z(1)=-c_im*p%k1*r
            z(2)=-c_im*p%k2*r
            call fbem_zexp_decomposed(2,z,EnR)
            EnR(2,:)=EnR(2,:)*d1r1
            EnR(3,:)=EnR(3,:)*d1r2
            EnR(4,:)=EnR(4,:)*d1r3
            EnR(5,:)=EnR(5,:)*d1r4
            EnR(6,:)=EnR(6,:)*d1r5
            TT1=p%T1(1)*d1r2+p%T1(2)+p%T1(3)*EnR(2,1)+p%T1(4)*EnR(2,2)+p%T1(5)*EnR(3,1)+p%T1(6)*EnR(3,2)+p%T1(7)*EnR(4,1)&
               +p%T1(8)*EnR(4,2)+p%T1(9)*EnR(5,1)+p%T1(10)*EnR(5,2)
            TT2=p%T2(1)*d1r2+p%T2(2)+p%T2(3)*EnR(2,2)+p%T2(4)*EnR(3,1)+p%T2(5)*EnR(3,2)+p%T2(6)*EnR(4,1)+p%T2(7)*EnR(4,2)&
               +p%T2(8)*EnR(5,1)+p%T2(9)*EnR(5,2)
            TT3=p%T3(1)*d1r2+p%T3(2)+p%T3(3)*EnR(2,1)+p%T3(4)*EnR(3,1)+p%T3(5)*EnR(3,2)+p%T3(6)*EnR(4,1)+p%T3(7)*EnR(4,2)&
               +p%T3(8)*EnR(5,1)+p%T3(9)*EnR(5,2)
            ! Add
            do il=1,3
              do ik=1,3
                fs_d=TT1*drdx(il)*drdx(ik)*drdni-TT2*(drdx(il)*n_i(ik)-c_dkr(il,ik)*drdni)-TT3*drdx(ik)*n_i(il)
                l(:,il,ik)=l(:,il,ik)+fs_d*sphijw
              end do
            end do
          end do
        end do
        case default
          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'n_edges not valid')
    end select
    l=p%cte_d*l
  end subroutine fbem_bem_harela3d_hbie_bl_ext_st

  recursive subroutine fbem_bem_harela3d_hbie_bl_ext_adp(e,xi_s,x_i,n_i,p,qsp,ks,ns,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                                !! Element
    real(kind=real64)                  :: xi_s(2,fbem_n_vertices(e%gtype)) !! Subdivision of the parent element
    real(kind=real64)                  :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)                  :: n_i(3)                           !! Unit normal at the collocation point
    type(fbem_bem_harela3d_parameters) :: p                                !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp                              !! Quasi-singular integration parameters
    integer                            :: ks                               !! Current level of subdivisions
    integer                            :: ns                               !! Maximum level of subdivision
    complex(kind=real64)               :: l(e%n_snodes,3,3)                !! l integration kernels matrix
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
    complex(kind=real64) :: l_tmp(e%n_snodes,3,3)                ! l integration kernels matrix (temporary)
    ! Initialize
    if (ks.eq.1) then
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
        call fbem_warning_message(error_unit,0,'fbem_bem_harela3d_hbie_bl_ext_adp',ns,'maximum number of subdivisions reached')
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
          call fbem_bem_harela3d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,l)
          ! SUBTRI 2
          tmp_xi_s(:,1)=xi_s(:,2)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          call fbem_bem_harela3d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,l)
          ! SUBTRI 3
          tmp_xi_s(:,1)=xi_s(:,3)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          call fbem_bem_harela3d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,l)
          ! SUBTRI 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          call fbem_bem_harela3d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,l)
        ! QUADRILATERALS
        case (4)
          ! SUBQUAD 1
          tmp_xi_s(:,1)=xi_s(:,1)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,3)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          call fbem_bem_harela3d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,l)
          ! SUBQUAD 2
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=xi_s(:,2)
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,4)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          call fbem_bem_harela3d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,l)
          ! SUBQUAD 3
          tmp_xi_s(:,1)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=xi_s(:,3)
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          call fbem_bem_harela3d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,l)
          ! SUBQUAD 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          tmp_xi_s(:,2)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=xi_s(:,4)
          call fbem_bem_harela3d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,l)
      end select
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_harela3d_hbie_bl_ext_st(e,xi_s,x_i,n_i,barxip,barr,p,gln,l_tmp)
      l=l+l_tmp
    end if
  end subroutine fbem_bem_harela3d_hbie_bl_ext_adp

  subroutine fbem_bem_harela3d_hbie_bl_int(type_g,type_f1,type_f2,delta_f,x_nodes,xi_i,p,l)
    implicit none
    ! I/O
    integer                            :: type_g                            !! Geometrial interpolation
    integer                            :: type_f1                           !! Functional interpolation (primary variables)
    integer                            :: type_f2                           !! Functional interpolation (secondary variables)
    real(kind=real64)                  :: delta_f                           !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64)                  :: x_nodes(3,fbem_n_nodes(type_g))   !! Position vectors of geometrical nodes
    real(kind=real64)                  :: xi_i(2)                           !! Reference coordinates of the singular point.
    type(fbem_bem_harela3d_parameters) :: p                                 !! Parameters of the region
    complex(kind=real64)               :: l(fbem_n_nodes(type_f2),3,3)      !! l kernel vector
    ! Local
    integer                            :: ksubtri                           ! Counter variable for subtriangles loop
    integer                            :: nsubtri                           ! Number of subtriangles
    integer                            :: subtriangle(8)                    ! Vector that contains what subtriangles need to be integrated
    real(kind=real64)                  :: theta_subtri(2,8)                 ! Matrix that contains the angles theta of the subtriangles to be integrated
    real(kind=real64)                  :: thetap_subtri(2,8)                ! Matrix that contains the angles thetap of the subtriangles to be integrated
    integer                            :: ktheta                            ! Counter variable for theta coordinate loop
    integer                            :: krho                              ! Counter variable for rho coordinate loop
    integer                            :: kphi                              ! Counter coordinates for shape functions loops
    integer                            :: nnodes_g                          ! Number of nodes of geometrical interpolation
    integer                            :: nnodes_f1                         ! Number of nodes of primary variables interpolation
    integer                            :: ngp_theta                         ! Number of Gauss points for theta coordinate
    integer                            :: ngp_rho                           ! Number of Gauss points for rho coordinate
    real(kind=real64)                  :: thetai, thetaf, thetapi, thetapf  ! Initial and final angles for subtriangle integration
    real(kind=real64)                  :: w_angular                         ! Weight of the angular coordinate
    real(kind=real64)                  :: w_radial                          ! Weight of the radial coordinate
    real(kind=real64)                  :: theta                             ! Angle coordinate theta
    real(kind=real64)                  :: thetap                            ! Angle coordinate thetap
    real(kind=real64)                  :: thetapp                           ! Angle coordinate thetap on [0,1] domain
    real(kind=real64)                  :: jthetap                           ! thetap->thetapp jacobian
    real(kind=real64)                  :: rhoij                             ! Maximum rho (radial) value for each edge
    real(kind=real64)                  :: rho                               ! Radial coordinate rho
    real(kind=real64)                  :: rhop                              ! Radial coordinate rho on [0,1] domain
    real(kind=real64)                  :: aux(10)                           ! Auxiliary variable for shape functions resources
    real(kind=real64)                  :: xi(2)                             ! Reference xi_1,xi_2 coordinates
    real(kind=real64)                  :: phi_f2(fbem_n_nodes(type_f2))     ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)                  :: phi_f1_i(fbem_n_nodes(type_f1))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)                  :: psi_i(3,fbem_n_nodes(type_f1))    ! psi vectors at xi_i1,xi_i2
    real(kind=real64)                  :: phi_f2_i(fbem_n_nodes(type_f2))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)                  :: phi_g(fbem_n_nodes(type_g))       ! Geometrical shape functions values at xi_1,xi_2
    real(kind=real64)                  :: dphidxi1_g(fbem_n_nodes(type_g))  ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64)                  :: dphidxi2_g(fbem_n_nodes(type_g))  ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64)                  :: jxi1, jxi2                        ! xi1->x, xi2->x jacobians
    real(kind=real64)                  :: x_i(3)                            ! Collocation point position vector
    real(kind=real64)                  :: n_i(3)                            ! Collocation point unit normal vector
    real(kind=real64)                  :: x(3)                              ! Position vector at xi_1,xi_2
    real(kind=real64)                  :: T1(3), T2(3)                      ! Tangent vectors at xi_1,xi_2
    real(kind=real64)                  :: v1(3), v2(3)                      ! Orthogonal tangent vectors
    real(kind=real64)                  :: Mvt(2,2)                          ! Orthogonal coordinates transformation matrix
    real(kind=real64)                  :: dxidv(2,2)                        ! xi coordinates derivatives with respect to v orthogonal coordinates
    real(kind=real64)                  :: detMvt                            ! Mvt(1,1)*Mvt(2,2)-Mvt(1,2)*Mvt(2,1)
    real(kind=real64)                  :: N(3)                              ! Normal vector at xi_1,xi_2
    real(kind=real64)                  :: rv(3)                             ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)                  :: r, dr1, dr2, dr3, dr4, dr5        ! Distance vector module and its inverse
    real(kind=real64)                  :: drdx(3)                           ! Distance vector derivatives with respect to x_k
    real(kind=real64)                  :: jg                                ! Geometric jacobian
    real(kind=real64)                  :: jw                                ! Jacobian * weights
    real(kind=real64)                  :: drdn                              ! Partial derivative of r respect to unit normal
    real(kind=real64)                  :: drdni                             ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)                  :: n_dot_ni                          ! Dot product n · n_i
    real(kind=real64)                  :: phif2jw(fbem_n_nodes(type_f2))    ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)                  :: costheta, sintheta                ! cos(theta), sin(theta)
    integer                            :: il, ik                            ! Counter for load / observation components
    complex(kind=real64)               :: z(2)                              ! Wavenumbers
    complex(kind=real64)               :: E(0:6,2)                          ! Exponential function decomposition for each wavenumber
    complex(kind=real64)               :: TT1, TT2, TT3                     ! Components of the fundamental solutions
    complex(kind=real64)               :: fs_d                              ! Fundamental solutions values
    !
    real(kind=real64)                  :: re                                ! Quasi-singular integration relative error
    integer                            :: gln_min                           ! Minimum Number of 1D Gauss-Legendre number of points
    integer                            :: ns_max                            ! Maximum number of subdivisions
    type(fbem_qs_parameters)           :: qsp                               ! Quasi-singular integration parameters
    real(kind=real64)                  :: xi_s(1,2)                         ! Edge subdivision
    integer                            :: kedge                             ! Counter variable for edges
    integer                            :: nedges                            ! Number of edges
    integer                            :: type_edge                         ! Line element type for line integrals
    integer                            :: nnodes_edge                       ! Number of nodes of the edge
    real(kind=real64), allocatable     :: x_nodes_edge(:,:)                 ! Coordinates of the edge elements
    real(kind=real64)                  :: lli(3,3)                          ! Line integral epsilon_{ijk}·e_k·t/r
    real(kind=real64)                  :: mli1                              ! Line integral (r x ni)·t/r**3
    real(kind=real64)                  :: mli2(3)                           ! Line integral (e_k x ni)·t/r
    real(kind=real64)                  :: mli3(3,3)                         ! Line integral (r x ni)·t*r_{,l}*r_{,k}/r**3
    real(kind=real64)                  :: mli4(3)                           ! Line integral (r x e_k)·t/r**3
    real(kind=real64)                  :: mli5(3,3,3)                       ! Line integral (e_k x ni)·t*r_{,l}*r_{,j}/r
    !
    ! Initialization
    !
    ! Kernels
    l=(0.d0,0.d0)
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
    ! Functional shape functions (secondary variables) at xi_i
#   define etype type_f2
#   define delta delta_f
#   define xi xi_i
#   define phi phi_f2_i
#   include <phi_2d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
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
          ! Unit normal vector
          n=N/jg
          ! Distance vector
          rv=x-x_i
          ! Distance vector norm
          r=dot_product(rv,rv)
          r=sqrt(r)
          dr1=1.d0/r
          dr2=dr1**2
          dr3=dr2*dr1
          dr4=dr3*dr1
          dr5=dr4*dr1
          ! r_{,k}
          drdx=rv*dr1
          ! dr/dn
          drdn=dot_product(drdx,n)
          ! dr/dni
          drdni=-dot_product(drdx,n_i)
          ! n·ni
          n_dot_ni=dot_product(n,n_i)
          ! Jacobians * weights
          jw=jg*rho*jthetap*w_angular*w_radial
          ! FUNCTIONAL SHAPE FUNCTIONS
          ! Functional shape functions (secondary variables) at xi
#         define etype type_f2
#         define delta delta_f
#         define phi phi_f2
#         include <phi_2d.rc>
#         undef etype
#         undef delta
#         undef phi
          ! Functional shape functions * jacobians * weights
          phif2jw=phi_f2*jw
          ! COMPONENTS OF THE FUNDAMENTAL SOLUTION (ONLY REGULAR PARTS)
          z(1)=-c_im*p%k1*r
          z(2)=-c_im*p%k2*r
          call fbem_zexp_decomposed(2,z,E)
          E(2,:)=E(2,:)*dr1
          E(3,:)=E(3,:)*dr2
          E(4,:)=E(4,:)*dr3
          E(5,:)=E(5,:)*dr4
          E(6,:)=E(6,:)*dr5
          TT1=p%T1(2)+p%T1(3)*E(2,1)+p%T1(4)*E(2,2)+p%T1(5)*E(3,1)+p%T1(6)*E(3,2)+p%T1(7)*E(4,1)+p%T1(8)*E(4,2)&
             +p%T1(9)*E(5,1)+p%T1(10)*E(5,2)
          TT2=p%T2(2)+p%T2(3)*E(2,2)+p%T2(4)*E(3,1)+p%T2(5)*E(3,2)+p%T2(6)*E(4,1)+p%T2(7)*E(4,2)+p%T2(8)*E(5,1)&
             +p%T2(9)*E(5,2)
          TT3=p%T3(2)+p%T3(3)*E(2,1)+p%T3(4)*E(3,1)+p%T3(5)*E(3,2)+p%T3(6)*E(4,1)+p%T3(7)*E(4,2)+p%T3(8)*E(5,1)&
             +p%T3(9)*E(5,2)
          ! Loop through load direction and observation direction
          do il=1,3
            do ik=1,3
              !
              ! Regular parts
              !
              ! Regular parts of fundamental solutions values without constants
              fs_d=(TT1+p%T1(1)*dr2)*drdx(il)*drdx(ik)*drdni+(TT2+p%T2(1)*dr2)*c_dkr(il,ik)*drdni-TT2*drdx(il)*n_i(ik)&
                  -TT3*drdx(ik)*n_i(il)+p%T2(1)*dr2*((n_i(il)-n(il))*drdx(ik)-(n_i(ik)-n(ik))*drdx(il))
              l(:,il,ik)=l(:,il,ik)+fs_d*phif2jw
              !
              ! Regular integrals after regularization (L)
              !
              ! L21
              if (il.ne.ik) then
                fs_d=p%T2(1)*dr2*(n(il)*drdx(ik)-n(ik)*drdx(il))
                l(:,il,ik)=l(:,il,ik)+fs_d*(phi_f2-phi_f2_i)*jw
              end if
            end do
          end do
        end do ! Loop through rho coordinate
      end do ! Loop through theta coordinate
    end do ! Loop through rectangle triangles
    !
    ! LINE INTEGRALS (all quasi-singular integrals, all integrations performed with an adaptive Telles+subdivision)
    !
    ! Initialize
    lli=0.d0
    mli1=0.d0
    mli2=0.d0
    mli3=0.d0
    mli4=0.d0
    mli5=0.d0
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
      ! Calculate all line integrals
      call fbem_bem_staela3d_hbie_int_li(type_edge,x_nodes_edge,xi_s,x_i,n_i,gln_min,qsp,1,ns_max,lli,mli1,mli2,mli3,mli4,mli5)
      ! Deallocate the edge coordinates
      deallocate (x_nodes_edge)
    end do ! Loop through edges
    !
    ! Add line integrals
    !
    ! Loop through load direction and observation direction
    do il=1,3
      do ik=1,3
        ! L22
        l(:,il,ik)=l(:,il,ik)+phi_f2_i*p%T2(1)*lli(il,ik)
      end do
    end do

    l=p%cte_d*l
  end subroutine fbem_bem_harela3d_hbie_bl_int

  subroutine fbem_bem_harela3d_hbie_bl_auto(e,x_i,n_i,p,qsp,ns,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                 !! Integration element
    real(kind=real64)                  :: x_i(3)            !! Collocation point
    real(kind=real64)                  :: n_i(3)            !! Unit normal at the collocation point
    type(fbem_bem_harela3d_parameters) :: p                 !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp               !! Quasi-singular integration parameters
    integer                            :: ns                !! Maximum level of subdivisions
    complex(kind=real64)               :: l(e%n_snodes,3,3) !! l integration kernel
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
        call fbem_bem_harela3d_hbie_bl_int(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,barxi,p,l)
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
          call fbem_bem_harela3d_hbie_bl_ext_pre(ps,e,x_i,n_i,p,l)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_harela3d_hbie_bl_ext_adp(e,xi_s,x_i,n_i,p,qsp,1,ns,l)
        end if
    end select
  end subroutine fbem_bem_harela3d_hbie_bl_auto
  ! --------------------------------------------------------------------------------------------------------------------------------

  ! ================================================================================================================================

  ! ================================================================================================================================
  ! BOTH SBIE AND HBIE

  ! Note:
  ! The following set of subroutines calculates SBIE and HBIE kernels for the same collocation point.
  ! If each integral equation has different collocation points, tipically it would be SBIE has nodal collocation and the
  ! HBIE has MCA collocation, then individual routines for SBIE and HBIE must be used. A comparison between both methodologies
  ! should be done particularly on the performance (computing effort) side..... It seems that, in general, using the same
  ! collocation point saves time since there are a lot of common calculations.

  subroutine fbem_bem_harela3d_shbie_ext_pre(ps,e,reverse,x_i,n_i,p,h,g,m,l)
    implicit none
    ! I/O
    integer                            :: ps                !! Selected precalculated dataset
    type(fbem_bem_element)             :: e                 !! Element
    logical                            :: reverse           !! Reverse normal vector
    real(kind=real64)                  :: x_i(3)            !! Collocation point position vector
    real(kind=real64)                  :: n_i(3)            !! Collocation point unit normal vector
    type(fbem_bem_harela3d_parameters) :: p                 !! Parameters of the region
    complex(kind=real64)               :: h(e%n_pnodes,3,3) !! h integration kernels vector
    complex(kind=real64)               :: g(e%n_snodes,3,3) !! g integration kernels vector
    complex(kind=real64)               :: m(e%n_pnodes,3,3) !! m integration kernels vector
    complex(kind=real64)               :: l(e%n_snodes,3,3) !! l integration kernels vector
    ! Local
    integer              :: kip                                            ! Counter variable for integration points loop
    real(kind=real64)    :: x(3)                            ! Position vector at integration point
    real(kind=real64)    :: n(3)                            ! Unit normal vector at integration point
    real(kind=real64)    :: pphijw(e%n_pnodes)              ! phi^p * jacobian * weight at integration point
    real(kind=real64)    :: sphijw(e%n_snodes)              ! phi^s * jacobian * weight at integration point
    real(kind=real64)    :: rv(3)                                          ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, d1r1, d1r2, d1r3, d1r4, d1r5                     ! Distance vector module and its inverse
    real(kind=real64)    :: drdx(3)                                        ! Distance vector derivatives with respect to x_k
    real(kind=real64)    :: drdn                                           ! Partial derivative of r respect to unit normal
    real(kind=real64)    :: drdni                                          ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)    :: n_dot_ni                                       ! Dot product n · n_i
    integer              :: il, ik                                         ! Counter for load / observation components
    complex(kind=real64) :: z(2)                                           ! Wavenumbers
    complex(kind=real64) :: EnR(0:6,2)                                       ! Exponential function decomposition for each wavenumber
    complex(kind=real64) :: psi, chi
    complex(kind=real64) :: TT1, TT2, TT3
    complex(kind=real64) :: S1, S2, S3, S4, S5        ! Components of the fundamental solutions
    complex(kind=real64) :: fs_u, fs_t, fs_d, fs_s                         ! Fundamental solutions values
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
      d1r4=d1r3*d1r1
      d1r5=d1r4*d1r1
      drdx=rv*d1r1
      drdn=dot_product(drdx,n)
      drdni=-dot_product(drdx,n_i)
      n_dot_ni=dot_product(n,n_i)
      z(1)=-c_im*p%k1*r
      z(2)=-c_im*p%k2*r
      call fbem_zexp_decomposed(2,z,EnR)
      EnR(2,:)=EnR(2,:)*d1r1
      EnR(3,:)=EnR(3,:)*d1r2
      EnR(4,:)=EnR(4,:)*d1r3
      EnR(5,:)=EnR(5,:)*d1r4
      EnR(6,:)=EnR(6,:)*d1r5
      psi=p%psi(1)*d1r1+p%psi(2)+EnR(2,2)+p%psi(3)*EnR(3,1)+p%psi(4)*EnR(3,2)+p%psi(5)*EnR(4,1)+p%psi(6)*EnR(4,2)
      chi=p%chi(1)*d1r1+p%chi(2)*EnR(2,1)+EnR(2,2)+p%chi(3)*EnR(3,1)+p%chi(4)*EnR(3,2)+p%chi(5)*EnR(4,1)+p%chi(6)*EnR(4,2)
      TT1=p%T1(1)*d1r2+p%T1(2)+p%T1(3)*EnR(2,1)+p%T1(4)*EnR(2,2)+p%T1(5)*EnR(3,1)+p%T1(6)*EnR(3,2)+p%T1(7)*EnR(4,1)&
         +p%T1(8)*EnR(4,2)+p%T1(9)*EnR(5,1)+p%T1(10)*EnR(5,2)
      TT2=p%T2(1)*d1r2+p%T2(2)+p%T2(3)*EnR(2,2)+p%T2(4)*EnR(3,1)+p%T2(5)*EnR(3,2)+p%T2(6)*EnR(4,1)+p%T2(7)*EnR(4,2)&
         +p%T2(8)*EnR(5,1)+p%T2(9)*EnR(5,2)
      TT3=p%T3(1)*d1r2+p%T3(2)+p%T3(3)*EnR(2,1)+p%T3(4)*EnR(3,1)+p%T3(5)*EnR(3,2)+p%T3(6)*EnR(4,1)+p%T3(7)*EnR(4,2)&
         +p%T3(8)*EnR(5,1)+p%T3(9)*EnR(5,2)
      S1=p%S1(1)*d1r3+p%S1(2)*d1r1+p%S1(3)*EnR(2,2)+p%S1(4)*EnR(3,1)+p%S1(5)*EnR(3,2)+p%S1(6)*EnR(4,1)+p%S1(7)*EnR(4,2)&
        +p%S1(8)*EnR(5,1)+p%S1(9)*EnR(5,2)+p%S1(10)*EnR(6,1)+p%S1(11)*EnR(6,2)
      S2=p%S2(1)*d1r3+p%S2(2)*d1r1+p%S2(3)*EnR(2,1)+p%S2(4)*EnR(3,1)+p%S2(5)*EnR(3,2)+p%S2(6)*EnR(4,1)+p%S2(7)*EnR(4,2)&
        +p%S2(8)*EnR(5,1)+p%S2(9)*EnR(5,2)+p%S2(10)*EnR(6,1)+p%S2(11)*EnR(6,2)
      S3=p%S3(1)*d1r3+p%S3(2)*d1r1+p%S3(3)*EnR(2,1)+p%S3(4)*EnR(2,2)+p%S3(5)*EnR(3,1)+p%S3(6)*EnR(3,2)+p%S3(7)*EnR(4,1)&
        +p%S3(8)*EnR(4,2)+p%S3(9)*EnR(5,1)+p%S3(10)*EnR(5,2)+p%S3(11)*EnR(6,1)+p%S3(12)*EnR(6,2)
      S4=p%S4(1)*d1r3+p%S4(2)*d1r1+p%S4(3)+p%S4(4)*EnR(3,2)+p%S4(5)*EnR(4,1)+p%S4(6)*EnR(4,2)+p%S4(7)*EnR(5,1)+p%S4(8)*EnR(5,2)&
        +p%S4(9)*EnR(6,1)+p%S4(10)*EnR(6,2)
      S5=p%S5(1)*d1r3+p%S5(2)*d1r1+p%S5(3)+p%S5(4)*EnR(2,1)+p%S5(5)*EnR(3,1)+p%S5(6)*EnR(4,1)+p%S5(7)*EnR(4,2)+p%S5(8)*EnR(5,1)&
        +p%S5(9)*EnR(5,2)+p%S5(10)*EnR(6,1)+p%S5(11)*EnR(6,2)
      do il=1,3
        do ik=1,3
          fs_u=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
          fs_t=TT1*drdx(il)*drdx(ik)*drdn+TT2*(drdn*c_dkr(il,ik)+drdx(ik)*n(il))+TT3*drdx(il)*n(ik)
          fs_d=TT1*drdx(il)*drdx(ik)*drdni-TT2*(-drdni*c_dkr(il,ik)+drdx(il)*n_i(ik))-TT3*drdx(ik)*n_i(il)
          fs_s=S1*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni-c_dkr(il,ik)*drdn*drdni+drdx(il)*drdx(ik)*n_dot_ni)&
              +S2*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni)+S3*drdx(il)*drdx(ik)*drdn*drdni&
              +S4*(c_dkr(il,ik)*n_dot_ni+n_i(ik)*n(il))+S5*n(ik)*n_i(il)
          h(:,il,ik)=h(:,il,ik)+fs_t*pphijw
          g(:,il,ik)=g(:,il,ik)+fs_u*sphijw
          m(:,il,ik)=m(:,il,ik)+fs_s*pphijw
          l(:,il,ik)=l(:,il,ik)+fs_d*sphijw
        end do
      end do
    end do
    ! Multiply by constants
    h=p%cte_t*h
    g=p%cte_u*g
    m=p%cte_s*m
    l=p%cte_d*l
    ! Reverse if needed
    if (reverse) then
      h=-h
      m=-m
    end if
  end subroutine fbem_bem_harela3d_shbie_ext_pre

  subroutine fbem_bem_harela3d_shbie_auto(e,reverse,x_i,n_i,p,qsp,ns,h,g,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                 !! Integration element
    logical                            :: reverse           !! Reverse orientation
    real(kind=real64)                  :: x_i(3)            !! Collocation point
    real(kind=real64)                  :: n_i(3)            !! Unit normal at the collocation point
    type(fbem_bem_harela3d_parameters) :: p                 !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp               !! Quasi-singular integration parameters
    integer                            :: ns                !! Maximum level of subdivisions
    complex(kind=real64)               :: h(e%n_pnodes,3,3) !! h integration kernel
    complex(kind=real64)               :: g(e%n_snodes,3,3) !! g integration kernel
    complex(kind=real64)               :: m(e%n_pnodes,3,3) !! m integration kernel
    complex(kind=real64)               :: l(e%n_snodes,3,3) !! l integration kernel
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
        call fbem_bem_harela3d_sbie_int(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi,p,h,g)
        call fbem_bem_harela3d_hbie_int(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi,p,m,l)
      case (0)
        ! Estimate the required integration rule
        gln_near=fbem_qs_n_estimation_standard(e%n,e%gtype,7,qsp,d,barxi)
        gln=max(e%gln_far,gln_near)
        ! Integrate using a conservative precalculated dataset
        if ((gln.le.e%ps_gln_max).and.(gln_near.gt.0)) then
          do i=1,e%n_ps
            if (e%ps_gln(i).ge.gln) then
              ps=i
              exit
            end if
          end do
          call fbem_bem_harela3d_shbie_ext_pre(ps,e,reverse,x_i,n_i,p,h,g,m,l)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_harela3d_sbie_ext_adp(e,reverse,xi_s,x_i,p,qsp,1,ns,h,g)
          call fbem_bem_harela3d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,p,qsp,1,ns,m,l)
        end if
    end select
  end subroutine fbem_bem_harela3d_shbie_auto

  ! ================================================================================================================================

end module fbem_bem_harela3d
