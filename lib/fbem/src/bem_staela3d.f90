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
!! <b> This module implements the calculation of BEM integration kernels vectors for 3D elastostatic problems.</b>
module fbem_bem_staela3d

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
  !use fbem_bem_stapot3d

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all are private
  private

  ! ================================================================================================================================
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE)
  ! Free-term
  public :: fbem_bem_staela3d_sbie_freeterm
  ! Fundamental solution
  public :: fbem_bem_staela3d_sbie_u, fbem_bem_staela3d_sbie_u_cmplx
  public :: fbem_bem_staela3d_sbie_t, fbem_bem_staela3d_sbie_t_cmplx
  ! BOUNDARY ELEMENTS
  ! Automatic integration
  public :: fbem_bem_staela3d_sbie_auto, fbem_bem_staela3d_sbie_auto_cmplx
  ! Exterior integration
  public :: fbem_bem_staela3d_sbie_ext_pre, fbem_bem_staela3d_sbie_ext_pre_cmplx
  public :: fbem_bem_staela3d_sbie_ext_st, fbem_bem_staela3d_sbie_ext_st_cmplx
  public :: fbem_bem_staela3d_sbie_ext_adp, fbem_bem_staela3d_sbie_ext_adp_cmplx
  ! Interior integration
  public :: fbem_bem_staela3d_sbie_int, fbem_bem_staela3d_sbie_int_cmplx
  public :: fbem_bem_staela3d_sbie_int_li
  ! BODY LOAD ELEMENTS
  ! Automatic integration
  public :: fbem_bem_staela3d_sbie_bl_auto, fbem_bem_staela3d_sbie_bl_auto_cmplx
  ! Exterior integration
  public :: fbem_bem_staela3d_sbie_bl_ext_pre, fbem_bem_staela3d_sbie_bl_ext_pre_cmplx
  public :: fbem_bem_staela3d_sbie_bl_ext_st, fbem_bem_staela3d_sbie_bl_ext_st_cmplx
  public :: fbem_bem_staela3d_sbie_bl_ext_adp, fbem_bem_staela3d_sbie_bl_ext_adp_cmplx
  ! Interior integration
  public :: fbem_bem_staela3d_sbie_bl_int, fbem_bem_staela3d_sbie_bl_int_cmplx
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (HBIE)
  ! Fundamental solution
  public :: fbem_bem_staela3d_hbie_d, fbem_bem_staela3d_hbie_d_cmplx
  public :: fbem_bem_staela3d_hbie_s, fbem_bem_staela3d_hbie_s_cmplx
  ! BOUNDARY ELEMENTS
  ! Automatic integration
  public :: fbem_bem_staela3d_hbie_auto
  ! Exterior integration
  public :: fbem_bem_staela3d_hbie_ext_pre
  public :: fbem_bem_staela3d_hbie_ext_st
  public :: fbem_bem_staela3d_hbie_ext_adp
  ! Interior integration
  public :: fbem_bem_staela3d_hbie_int
  public :: fbem_bem_staela3d_hbie_int_li
  public :: fbem_bem_staela3d_hbie_int_mli3
  public :: fbem_bem_staela3d_hbie_int_mli4
  public :: fbem_bem_staela3d_hbie_int_mli5
  ! BODY LOAD ELEMENTS
  ! Automatic integration
  public :: fbem_bem_staela3d_hbie_bl_auto
  ! Exterior integration
  public :: fbem_bem_staela3d_hbie_bl_ext_pre
  public :: fbem_bem_staela3d_hbie_bl_ext_st
  public :: fbem_bem_staela3d_hbie_bl_ext_adp
  ! Interior integration
  !public :: fbem_bem_staela3d_hbie_bl_int
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE) - TRACTION-FREE HALF-SPACE FUNDAMENTAL SOLUTION
  public :: fbem_bem_staela3d_hs_sbie, fbem_bem_staela3d_hs_sbie_cmplx
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE) - TRACTION-FREE HALF-SPACE FUNDAMENTAL SOLUTION (ONLY COMPLEMENTARY PART)
  public :: fbem_bem_staela3d_hsc_sbie, fbem_bem_staela3d_hsc_sbie_cmplx
  ! BOUNDARY ELEMENTS
  ! Automatic integration
  public :: fbem_bem_staela3d_hsc_sbie_auto, fbem_bem_staela3d_hsc_sbie_auto_cmplx
  ! Exterior integration
  public :: fbem_bem_staela3d_hsc_sbie_ext_pre, fbem_bem_staela3d_hsc_sbie_ext_pre_cmplx
  public :: fbem_bem_staela3d_hsc_sbie_ext_st, fbem_bem_staela3d_hsc_sbie_ext_st_cmplx
  public :: fbem_bem_staela3d_hsc_sbie_ext_adp, fbem_bem_staela3d_hsc_sbie_ext_adp_cmplx
  ! Interior integration
  public :: fbem_bem_staela3d_hsc_sbie_int, fbem_bem_staela3d_hsc_sbie_int_cmplx
  ! BODY LOAD ELEMENTS
  ! Automatic integration
  public :: fbem_bem_staela3d_hsc_sbie_bl_auto, fbem_bem_staela3d_hsc_sbie_bl_auto_cmplx
  ! Exterior integration
  public :: fbem_bem_staela3d_hsc_sbie_bl_ext_pre, fbem_bem_staela3d_hsc_sbie_bl_ext_pre_cmplx
  public :: fbem_bem_staela3d_hsc_sbie_bl_ext_st, fbem_bem_staela3d_hsc_sbie_bl_ext_st_cmplx
  public :: fbem_bem_staela3d_hsc_sbie_bl_ext_adp, fbem_bem_staela3d_hsc_sbie_bl_ext_adp_cmplx
  ! Interior integration
  public :: fbem_bem_staela3d_hsc_sbie_bl_int, fbem_bem_staela3d_hsc_sbie_bl_int_cmplx
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (HBIE) - TRACTION-FREE HALF-SPACE FUNDAMENTAL SOLUTION (COMPLEMENTARY PART)
  public :: fbem_bem_staela3d_hsc_hbie
  ! BOUNDARY ELEMENTS
  ! Automatic integration
  public :: fbem_bem_staela3d_hsc_hbie_auto
  ! Exterior integration
  public :: fbem_bem_staela3d_hsc_hbie_ext_pre
  public :: fbem_bem_staela3d_hsc_hbie_ext_st
  public :: fbem_bem_staela3d_hsc_hbie_ext_adp
  ! ================================================================================================================================




contains

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE)
  ! --------------------------------------------------------------------------------------------------------------------------------

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
  subroutine fbem_bem_staela3d_sbie_freeterm(n_elements,n,t,tol,nu,c)
    implicit none
    ! I/O
    integer            :: n_elements            !! Number of elements
    real(kind=real64)  :: n(3,n_elements)       !! Unit normal of each element
    real(kind=real64)  :: t(3,n_elements)       !! Unit element boundary tangent of each element
    real(kind=real64)  :: tol                   !! Geometric tolerance [1.0e-12,1.0e-3] (default 1.0e-6)
    real(kind=real64)  :: nu                    !! Poisson's ratio
    real(kind=real64)  :: c(3,3)                !! Free-term matrix
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
  end subroutine fbem_bem_staela3d_sbie_freeterm

  !! Fundamental solution u*
  subroutine fbem_bem_staela3d_sbie_u(x,x_i,mu,nu,uo)
    implicit none
    ! I/O
    real(kind=real64) :: x(3)    !! Observation point
    real(kind=real64) :: x_i(3)  !! Collocation point
    real(kind=real64) :: mu      !! Shear modulus
    real(kind=real64) :: nu      !! Poisson's ratio
    real(kind=real64) :: uo(3,3) !! u*_{lk}
    ! Local
    integer           :: l, k
    real(kind=real64) :: cteu1, cteu2
    real(kind=real64) :: rv(3), r, d1r, drdx(3)
    cteu1=1.d0/(16.d0*c_pi*mu*(1.d0-nu))
    cteu2=3.d0-4.d0*nu
    rv=x-x_i
    r=sqrt(dot_product(rv,rv))
    d1r=1.d0/r
    drdx=rv*d1r
    do l=1,3
      do k=1,3
        uo(l,k)=d1r*(cteu2*c_dkr(l,k)+drdx(l)*drdx(k))
      end do
    end do
    uo=cteu1*uo
  end subroutine fbem_bem_staela3d_sbie_u

  !! Fundamental solution t*
  subroutine fbem_bem_staela3d_sbie_t(x,n,x_i,nu,to)
    implicit none
    ! I/O
    real(kind=real64) :: x(3)    !! Observation point
    real(kind=real64) :: n(3)    !! Observation point normal
    real(kind=real64) :: x_i(3)  !! Collocation point
    real(kind=real64) :: nu      !! Poisson's ratio
    real(kind=real64) :: to(3,3) !! t*_{lk}
    ! Local
    integer           :: l, k
    real(kind=real64) :: rv(3), r, d1r, d1r2, drdx(3), drdn
    real(kind=real64) :: ctet1, ctet2
    ctet1=-1.d0/(8.d0*c_pi*(1.d0-nu))
    ctet2=1.d0-2.d0*nu
    rv=x-x_i
    r=sqrt(dot_product(rv,rv))
    d1r=1.d0/r
    d1r2=d1r**2
    drdx=rv*d1r
    drdn=dot_product(drdx,n)
    do l=1,3
      do k=1,3
        to(l,k)=d1r2*((ctet2*c_dkr(l,k)+3.d0*drdx(l)*drdx(k))*drdn+ctet2*(n(l)*drdx(k)-n(k)*drdx(l)))
      end do
    end do
    to=ctet1*to
  end subroutine fbem_bem_staela3d_sbie_t

  !! Fundamental solution u*
  subroutine fbem_bem_staela3d_sbie_u_cmplx(x,x_i,mu,nu,uo)
    implicit none
    ! I/O
    real(kind=real64)    :: x(3)    !! Observation point
    real(kind=real64)    :: x_i(3)  !! Collocation point
    complex(kind=real64) :: mu      !! Shear modulus
    complex(kind=real64) :: nu      !! Poisson's ratio
    complex(kind=real64) :: uo(3,3) !! u*_{lk}
    ! Local
    integer              :: l, k
    complex(kind=real64) :: cteu1, cteu2
    real(kind=real64)    :: rv(3), r, d1r, drdx(3)
    cteu1=1.d0/(16.d0*c_pi*mu*(1.d0-nu))
    cteu2=3.d0-4.d0*nu
    rv=x-x_i
    r=sqrt(dot_product(rv,rv))
    d1r=1.d0/r
    drdx=rv*d1r
    do l=1,3
      do k=1,3
        uo(l,k)=d1r*(cteu2*c_dkr(l,k)+drdx(l)*drdx(k))
      end do
    end do
    uo=cteu1*uo
  end subroutine fbem_bem_staela3d_sbie_u_cmplx

  !! Fundamental solution t*
  subroutine fbem_bem_staela3d_sbie_t_cmplx(x,n,x_i,nu,to)
    implicit none
    ! I/O
    real(kind=real64)    :: x(3)    !! Observation point
    real(kind=real64)    :: n(3)    !! Observation point normal
    real(kind=real64)    :: x_i(3)  !! Collocation point
    complex(kind=real64) :: nu      !! Poisson's ratio
    complex(kind=real64) :: to(3,3) !! t*_{lk}
    ! Local
    integer              :: l, k
    real(kind=real64)    :: rv(3), r, d1r, d1r2, drdx(3), drdn
    complex(kind=real64) :: ctet1, ctet2
    ctet1=-1.d0/(8.d0*c_pi*(1.d0-nu))
    ctet2=1.d0-2.d0*nu
    rv=x-x_i
    r=sqrt(dot_product(rv,rv))
    d1r=1.d0/r
    d1r2=d1r**2
    drdx=rv*d1r
    drdn=dot_product(drdx,n)
    do l=1,3
      do k=1,3
        to(l,k)=d1r2*((ctet2*c_dkr(l,k)+3.d0*drdx(l)*drdx(k))*drdn+ctet2*(n(l)*drdx(k)-n(k)*drdx(l)))
      end do
    end do
    to=ctet1*to
  end subroutine fbem_bem_staela3d_sbie_t_cmplx

  ! ==================== !
  ! BE BOUNDARY ELEMENTS !
  ! ==================== !

  !! Efficient automatic integration
  subroutine fbem_bem_staela3d_sbie_auto(e,reverse,x_i,mu,nu,qsp,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                !! Integration element
    logical                  :: reverse                          !! Reverse orientation
    real(kind=real64)        :: x_i(3)                           !! Collocation point
    real(kind=real64)        :: mu                               !! Shear modulus
    real(kind=real64)        :: nu                               !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp                              !! Quasi-singular integration parameters
    integer                  :: ns                               !! Maximum level of subdivisions
    real(kind=real64)        :: h(e%n_pnodes,3,3)                !! h integration kernel
    real(kind=real64)        :: g(e%n_snodes,3,3)                !! g integration kernel
    ! Local
    real(kind=real64)        :: r(3)                             ! Distance vector
    real(kind=real64)        :: rmin                             ! Minimum distance between element and x_i
    real(kind=real64)        :: barxi(2)                         ! Local coordinates of the nearest element point with respect to x_i
    real(kind=real64)        :: d                                ! Dimensionless distance
    integer                  :: delta                            ! Control variable
    real(kind=real64)        :: xi_s(2,fbem_n_vertices(e%gtype)) ! Local coordinates of the element subdivision
    integer                  :: method                           ! Method used when calculating the nearest element point
    integer                  :: gln_near                         ! 1D Gauss-Legendre integration points required by the quasi-singular function
    integer                  :: gln                              ! 1D Gauss-Legendre integration points used in the integration
    integer                  :: ps                               ! Selected precalculated dataset
    integer                  :: i                                ! Counter
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
        call fbem_bem_staela3d_sbie_int(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi,mu,nu,h,g)
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
          call fbem_bem_staela3d_sbie_ext_pre(ps,e,reverse,x_i,mu,nu,h,g)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_staela3d_sbie_ext_adp(e,reverse,xi_s,x_i,mu,nu,qsp,1,ns,h,g)
        end if
    end select
  end subroutine fbem_bem_staela3d_sbie_auto

  !! This subroutine calculates the kernels for SBIE exterior integration needing precalculated data at integration points.
  subroutine fbem_bem_staela3d_sbie_ext_pre(ps,e,reverse,x_i,mu,nu,h,g)
    implicit none
    ! I/O
    integer                :: ps                !! Selected precalculated dataset
    type(fbem_bem_element) :: e                 !! Element
    logical                :: reverse           !! Reverse normal vector
    real(kind=real64)      :: x_i(3)            !! Collocation point position vector
    real(kind=real64)      :: mu                !! Shear modulus
    real(kind=real64)      :: nu                !! Poisson's ratio
    real(kind=real64)      :: h(e%n_pnodes,3,3) !! h integration kernels matrix
    real(kind=real64)      :: g(e%n_snodes,3,3) !! g integration kernels matrix
    ! Local
    integer           :: il, ik                     ! Counter for load / observation components
    integer           :: kip                        ! Counter variable for integration points loop
    real(kind=real64) :: x(3)                       ! Position vector at integration point
    real(kind=real64) :: n(3)                       ! Unit normal vector at integration point
    real(kind=real64) :: pphijw(e%n_pnodes)         ! phi^p * jacobian * weight at integration point
    real(kind=real64) :: sphijw(e%n_snodes)         ! phi^s * jacobian * weight at integration point
    real(kind=real64) :: rv(3)                      ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64) :: r, d1r, d1r2               ! Distance vector module and its inverse
    real(kind=real64) :: drdx(3)                    ! Distance vector derivatives with respect to x_k
    real(kind=real64) :: drdn                       ! Partial derivative of r respect to unit normal
    real(kind=real64) :: cteu1, cteu2, ctet1, ctet2 ! Auxiliary constants
    real(kind=real64) :: fs_u, fs_t                 ! Fundamental solutions values
    ! Initialize kernels
    h=0.d0
    g=0.d0
    ! Initialize auxiliary constants for fundamental solutions calculation
    cteu1=1.d0/(16.d0*c_pi*mu*(1.d0-nu))
    cteu2=3.d0-4.d0*nu
    ctet1=-1.d0/(8.d0*c_pi*(1.d0-nu))
    ctet2=1.d0-2.d0*nu
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      x=e%ps_x(:,kip,ps)
      n=e%ps_n(:,kip,ps)
      pphijw=e%ps_pphijw(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      d1r=1.d0/r
      d1r2=d1r**2
      drdx=rv*d1r
      drdn=dot_product(drdx,n)
      do il=1,3
        do ik=1,3
          fs_u=d1r*(cteu2*c_dkr(il,ik)+drdx(il)*drdx(ik))
          fs_t=d1r2*((ctet2*c_dkr(il,ik)+3.d0*drdx(il)*drdx(ik))*drdn+ctet2*(n(il)*drdx(ik)-n(ik)*drdx(il)))
          h(:,il,ik)=h(:,il,ik)+fs_t*pphijw
          g(:,il,ik)=g(:,il,ik)+fs_u*sphijw
        end do
      end do
    end do ! Loop through integrations points
    ! Multiply by constants
    h=ctet1*h
    g=cteu1*g
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_staela3d_sbie_ext_pre

  !! This subroutine calculates the kernels for SBIE exterior integration (near collocation points) using Telles transformation
  !! within a subdivision of the element, needing only needs basic data.
  subroutine fbem_bem_staela3d_sbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,mu,nu,gln,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)       :: e                                !! Integration element
    logical                      :: reverse                          !! Reverse normal vector
    real(kind=real64)            :: xi_s(2,fbem_n_vertices(e%gtype)) !! Coordinates of the subdivision of the parent element
    real(kind=real64)            :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)            :: barxip(2)                        !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)            :: barr                             !! Telles jacobian at barxip
    real(kind=real64)            :: mu                               !! Shear modulus
    real(kind=real64)            :: nu                               !! Poisson's ratio
    integer                      :: gln                              !! 1D Gauss-Legendre number of integration points (<=32)
    real(kind=real64)            :: h(e%n_pnodes,3,3)                !! h kernel vector
    real(kind=real64)            :: g(e%n_snodes,3,3)                !! g kernel vector
    ! Local
    integer                      :: kphi                       ! Counter variable for shape functions loops
    integer                      :: k1                         ! Counter variable for reference coordinate xi_1
    integer                      :: k2                         ! Counter variable for reference coordinate xi_2
    real(kind=real64)            :: aux(10)                    ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)           ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi1(e%n_gnodes)      ! Geometrical shape functions derivatives values
    real(kind=real64)            :: dgphidxi2(e%n_gnodes)      ! Geometrical shape functions derivatives values
    real(kind=real64)            :: pphi(e%n_pnodes)           ! Functional shape functions values
    real(kind=real64)            :: sphi(e%n_snodes)           ! Functional shape functions values
    real(kind=real64)            :: gamma(2)                   ! Vector of gamma_1,gamma_2 coordinates
    real(kind=real64)            :: w(2)                       ! Weights of the integration rule
    real(kind=real64)            :: xip(2)                     ! Vector of xip_1,xip_2 coordinates
    real(kind=real64)            :: dxidxi1p(2)                ! xi derivatives with respect to xi1p
    real(kind=real64)            :: dxidxi2p(2)                ! xi derivatives with respect to xi2p
    real(kind=real64)            :: js                         ! Subdivision jacobian
    real(kind=real64)            :: xi(2)                      ! Vector of xi_1,xi_2 coordinates
    real(kind=real64)            :: xipp(2)                    ! Coordinate xipp used for quadrilateral-triangle transformation
    real(kind=real64)            :: barxipp(2)                 ! Coordinate xipp of collocation point
    real(kind=real64)            :: jqt                        ! Jacobian of the quadrilateral-triangle transformation
    type(fbem_telles_parameters) :: telles_parameters(2)       ! Telles parameters for each coordinate
    real(kind=real64)            :: jt(2)                      ! Telles jacobian for each coordinate: xi_1->gamma_1 and xi_2->gamma_2
    real(kind=real64)            :: x(3)                       ! Position vector at xi_1,xi_2
    real(kind=real64)            :: T1(3), T2(3)               ! Tangent vectors at xi_1,xi_2
    real(kind=real64)            :: N(3)                       ! Normal vector at xi_1,xi_2
    real(kind=real64)            :: rv(3)                      ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)            :: r, d1r, d1r2               ! Distance vector module and its inverse
    real(kind=real64)            :: drdx(3)                    ! Distance vector derivatives with respect to x_k
    real(kind=real64)            :: jg                         ! Geometric jacobian
    real(kind=real64)            :: jw                         ! Jacobians * weights
    real(kind=real64)            :: drdn                       ! Partial derivative of r respect to unit normal
    real(kind=real64)            :: pphijw(e%n_pnodes)         ! Functional shape functions * jw
    real(kind=real64)            :: sphijw(e%n_snodes)         ! Functional shape functions * jw
    integer                      :: il, ik                     ! Counter for load / observation components
    real(kind=real64)            :: cteu1, cteu2, ctet1, ctet2 ! Auxiliary constants
    real(kind=real64)            :: fs_u, fs_t                 ! Fundamental solutions values
    ! Initialize kernels
    h=0.d0
    g=0.d0
    ! Initialize auxiliary constants for fundamental solutions calculation
    cteu1=1.d0/(16.d0*c_pi*mu*(1.d0-nu))
    cteu2=3.d0-4.d0*nu
    ctet1=-1.d0/(8.d0*c_pi*(1.d0-nu))
    ctet2=1.d0-2.d0*nu
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
            d1r=1.d0/r
            d1r2=d1r**2
            drdx=rv*d1r
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
            ! Add integration points
            do il=1,3
              do ik=1,3
                fs_u=d1r*(cteu2*c_dkr(il,ik)+drdx(il)*drdx(ik))
                fs_t=d1r2*((ctet2*c_dkr(il,ik)+3.d0*drdx(il)*drdx(ik))*drdn+ctet2*(n(il)*drdx(ik)-n(ik)*drdx(il)))
                h(:,il,ik)=h(:,il,ik)+fs_t*pphijw
                g(:,il,ik)=g(:,il,ik)+fs_u*sphijw
              end do
            end do
          end do
        end do
      ! TRIANGULAR ELEMENTS
      case (3)
        ! Telles transformation is applied to Gauss-Legendre*Gauss-Legendre quadrature before the quadrilateral->triangle
        ! transformation. Because barxip for triangles are given in area triangle coordinates, for Telles transformation
        ! they must be transformed to quadrilateral coordinates. A special treatment is needed when barxi_2 is near 1, because
        ! transformation diverges.
        ! BARXIP->BARXIPP TRANSFORMATION
        if (barxip(2).gt.0.995d0) then
          barxipp(1)=0.5d0
          barxipp(2)=1.0d0
        else
          barxipp(1)=barxip(1)/(1.d0-barxip(2))
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
            d1r=1.d0/r
            d1r2=d1r**2
            drdx=rv*d1r
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
            ! Add integration points
            do il=1,3
              do ik=1,3
                fs_u=d1r*(cteu2*c_dkr(il,ik)+drdx(il)*drdx(ik))
                fs_t=d1r2*((ctet2*c_dkr(il,ik)+3.d0*drdx(il)*drdx(ik))*drdn+ctet2*(n(il)*drdx(ik)-n(ik)*drdx(il)))
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
    h=ctet1*h
    g=cteu1*g
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_staela3d_sbie_ext_st

  !! This subroutine calculates adaptatively the kernels for SBIE exterior integration using Telles transformation and subdivision
  !! if needed, using only basic data.
  recursive subroutine fbem_bem_staela3d_sbie_ext_adp(e,reverse,xi_s,x_i,mu,nu,qsp,ks,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                !! Element
    logical                  :: reverse                          !! Reverse orientation
    real(kind=real64)        :: xi_s(2,fbem_n_vertices(e%gtype)) !! Subdivision of the parent element
    real(kind=real64)        :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)        :: mu                               !! Shear modulus
    real(kind=real64)        :: nu                               !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp                              !! Quasi-singular integration parameters
    integer                  :: ks                               !! Current level of subdivisions
    integer                  :: ns                               !! Maximum level of subdivision
    real(kind=real64)        :: h(e%n_pnodes,3,3)                !! h integration kernels matrix
    real(kind=real64)        :: g(e%n_snodes,3,3)                !! g integration kernels matrix
    ! Local
    integer           :: gln_near                             ! 1D Gauss-Legendre integ. points required to integrate only  the quasi-singular integrand
    integer           :: gln                                  ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
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
    real(kind=real64) :: h_tmp(e%n_pnodes,3,3)                ! h integration kernels matrix (temporary)
    real(kind=real64) :: g_tmp(e%n_snodes,3,3)                ! g integration kernels matrix (temporary)
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
      cl=fbem_characteristic_length(3,e%gtype,x_s,1.d-12)
      call fbem_nearest_element_point_bem(3,e%gtype,x_s,cl,x_i,barxip,rmin,d,method)
    end if
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,5,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_staela3d_sbie_ext_adp',ns,'maximum number of subdivisions reached')
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
          call fbem_bem_staela3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,h,g)
          ! SUBTRI 2
          tmp_xi_s(:,1)=xi_s(:,2)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          call fbem_bem_staela3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,h,g)
          ! SUBTRI 3
          tmp_xi_s(:,1)=xi_s(:,3)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          call fbem_bem_staela3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,h,g)
          ! SUBTRI 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          call fbem_bem_staela3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,h,g)
        ! QUADRILATERALS
        case (4)
          ! SUBQUAD 1
          tmp_xi_s(:,1)=xi_s(:,1)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,3)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          call fbem_bem_staela3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,h,g)
          ! SUBQUAD 2
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=xi_s(:,2)
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,4)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          call fbem_bem_staela3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,h,g)
          ! SUBQUAD 3
          tmp_xi_s(:,1)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=xi_s(:,3)
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          call fbem_bem_staela3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,h,g)
          ! SUBQUAD 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          tmp_xi_s(:,2)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=xi_s(:,4)
          call fbem_bem_staela3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,h,g)
      end select
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_staela3d_sbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,mu,nu,gln,h_tmp,g_tmp)
      h=h+h_tmp
      g=g+g_tmp
    end if
  end subroutine fbem_bem_staela3d_sbie_ext_adp

  !! This subroutine calculates static elastic integration kernels vectors h and g of singular formulation for interior
  !! integration.
  subroutine fbem_bem_staela3d_sbie_int(type_g,type_f1,type_f2,delta_f,x_nodes,reverse,xi_i,mu,nu,h,g)
    implicit none
    ! I/O
    integer           :: type_g                          !! Geometrial interpolation
    integer           :: type_f1                         !! Functional interpolation (primary variables)
    integer           :: type_f2                         !! Functional interpolation (secondary variables)
    real(kind=real64) :: delta_f                         !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64) :: x_nodes(3,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
    logical           :: reverse                         !! Normal vector inversion
    real(kind=real64) :: xi_i(2)                         !! Reference coordinates of the singular point.
    real(kind=real64) :: mu                              !! Shear modulus
    real(kind=real64) :: nu                              !! Poisson's ratio
    real(kind=real64) :: h(fbem_n_nodes(type_f1),3,3)    !! h kernel vector
    real(kind=real64) :: g(fbem_n_nodes(type_f2),3,3)    !! g kernel vector
    ! Local
    integer           :: ksubtri                          ! Counter variable for subtriangles loop
    integer           :: nsubtri                          ! Number of subtriangles
    integer           :: subtriangle(8)                   ! Vector that contains what subtriangles need to be integrated
    real(kind=real64) :: theta_subtri(2,8)                ! Matrix that contains the angles theta of the subtriangles to be integrated
    real(kind=real64) :: thetap_subtri(2,8)               ! Matrix that contains the angles thetap of the subtriangles to be integrated
    integer           :: ktheta                           ! Counter variable for theta coordinate loop
    integer           :: krho                             ! Counter variable for rho coordinate loop
    integer           :: kphi                             ! Counter coordinates for shape functions loops
    integer           :: nnodes_g                         ! Number of nodes of geometrical interpolation
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
    real(kind=real64) :: phi_f1(fbem_n_nodes(type_f1))    ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phi_f2(fbem_n_nodes(type_f2))    ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phi_f1_i(fbem_n_nodes(type_f1))  ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phi_g(fbem_n_nodes(type_g))      ! Geometrical shape functions values at xi_1,xi_2
    real(kind=real64) :: dphidxi1_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64) :: dphidxi2_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64) :: x_i(3)                           ! Collocation point position vector
    real(kind=real64) :: x(3)                             ! Position vector at xi_1,xi_2
    real(kind=real64) :: T1(3), T2(3)                     ! Tangent vectors at xi_1,xi_2
    real(kind=real64) :: N(3)                             ! Normal vector at xi_1,xi_2
    real(kind=real64) :: rv(3)                            ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64) :: r, d1r, d1r2                     ! Distance vector module and its inverse
    real(kind=real64) :: drdx(3)                          ! Distance vector derivatives with respect to x_k
    real(kind=real64) :: jg                               ! Geometric jacobian
    real(kind=real64) :: jw                               ! Jacobian * weights
    real(kind=real64) :: drdn                             ! Partial derivative of r respect to unit normal
    real(kind=real64) :: phif1jw(fbem_n_nodes(type_f1))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phif2jw(fbem_n_nodes(type_f2))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: costheta, sintheta               ! cos(theta), sin(theta)
    integer           :: il, ik                           ! Counter for load / observation components
    real(kind=real64) :: cteu1, cteu2, ctet1, ctet2       ! Auxiliary constants
    real(kind=real64) :: fs_u, fs_t                       ! Fundamental solutions values
    ! Local variables associated with line integrals
    real(kind=real64)              :: re                                   ! Quasi-singular integration relative error
    integer                        :: gln_min                              ! Minimum Number of 1D Gauss-Legendre number of points
    integer                        :: ns_max                               ! Maximum number of subdivisions
    type(fbem_qs_parameters)       :: qsp                                  ! Quasi-singular integration parameters
    real(kind=real64)              :: xi_s(1,2)                            ! Edge subdivision
    integer                        :: kedge                                ! Counter variable for edges
    integer                        :: nedges                               ! Number of edges
    logical                        :: integrate                            ! Control variable
    integer                        :: type_edge                            ! Line element type for line integrals
    integer                        :: nnodes_edge                          ! Number of nodes of the edge
    real(kind=real64), allocatable :: x_nodes_edge(:,:)                    ! Coordinates of the edge elements
    real(kind=real64)              :: hli(3,3)                             ! Line integrals values
    !
    ! Initialization
    !
    ! Kernel vectors
    h=0.d0
    g=0.d0
    ! Number of nodes of geometrical interpolation
    nnodes_g=fbem_n_nodes(type_g)
    ! Initialize auxiliary constants for fundamental solutions calculation
    cteu1=1.d0/(16.d0*c_pi*mu*(1.d0-nu))
    cteu2=3.d0-4.d0*nu
    ctet1=-1.d0/(8.d0*c_pi*(1.d0-nu))
    ctet2=1.d0-2.d0*nu
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
          jg=sqrt(dot_product(N,N))
          ! Unit normal vector
          n=N/jg
          ! Distance vector
          rv=x-x_i
          ! Distance vector norm
          r=dot_product(rv,rv)
          r=sqrt(r)
          d1r=1.d0/r
          d1r2=d1r**2
          ! r_{,k}
          drdx=rv*d1r
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
          ! Loop through load direction and observation direction
          do il=1,3
            do ik=1,3
              ! Fundamental solutions values without cteu1 and ctep1, respectively
              fs_u=d1r*(cteu2*c_dkr(il,ik)+drdx(il)*drdx(ik))
              ! Regular part of the fundamental solution
              fs_t=d1r2*drdn*(ctet2*c_dkr(il,ik)+3.d0*drdx(il)*drdx(ik))
              ! Add to kernels
              h(:,il,ik)=h(:,il,ik)+fs_t*phif1jw
              g(:,il,ik)=g(:,il,ik)+fs_u*phif2jw
              ! The regular part of the CPV integral of h
              fs_t=d1r2*ctet2*(n(il)*drdx(ik)-n(ik)*drdx(il))
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
        h(:,il,ik)=h(:,il,ik)+phi_f1_i*ctet2*hli(il,ik)
      end do
    end do
    ! Multiply h by ctet1 and g by cteu1
    h=ctet1*h
    g=cteu1*g
    ! If the normal has to be reversed, then h=-h
    if (reverse) h=-h
  end subroutine fbem_bem_staela3d_sbie_int

  !! Efficient automatic integration
  subroutine fbem_bem_staela3d_sbie_auto_cmplx(e,reverse,x_i,mu,nu,qsp,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                !! Integration element
    logical                  :: reverse                          !! Reverse orientation
    real(kind=real64)        :: x_i(3)                           !! Collocation point
    complex(kind=real64)        :: mu                               !! Shear modulus
    complex(kind=real64)        :: nu                               !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp                              !! Quasi-singular integration parameters
    integer                  :: ns                               !! Maximum level of subdivisions
    complex(kind=real64)        :: h(e%n_pnodes,3,3)                !! h integration kernel
    complex(kind=real64)        :: g(e%n_snodes,3,3)                !! g integration kernel
    ! Local
    real(kind=real64)        :: r(3)                             ! Distance vector
    real(kind=real64)        :: rmin                             ! Minimum distance between element and x_i
    real(kind=real64)        :: barxi(2)                         ! Local coordinates of the nearest element point with respect to x_i
    real(kind=real64)        :: d                                ! Dimensionless distance
    integer                  :: delta                            ! Control variable
    real(kind=real64)        :: xi_s(2,fbem_n_vertices(e%gtype)) ! Local coordinates of the element subdivision
    integer                  :: method                           ! Method used when calculating the nearest element point
    integer                  :: gln_near                         ! 1D Gauss-Legendre integration points required by the quasi-singular function
    integer                  :: gln                              ! 1D Gauss-Legendre integration points used in the integration
    integer                  :: ps                               ! Selected precalculated dataset
    integer                  :: i                                ! Counter
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
        call fbem_bem_staela3d_sbie_int_cmplx(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi,mu,nu,h,g)
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
          call fbem_bem_staela3d_sbie_ext_pre_cmplx(ps,e,reverse,x_i,mu,nu,h,g)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_staela3d_sbie_ext_adp_cmplx(e,reverse,xi_s,x_i,mu,nu,qsp,1,ns,h,g)
        end if
    end select
  end subroutine fbem_bem_staela3d_sbie_auto_cmplx

  !! This subroutine calculates the kernels for SBIE exterior integration needing precalculated data at integration points.
  subroutine fbem_bem_staela3d_sbie_ext_pre_cmplx(ps,e,reverse,x_i,mu,nu,h,g)
    implicit none
    ! I/O
    integer                :: ps                !! Selected precalculated dataset
    type(fbem_bem_element) :: e                 !! Element
    logical                :: reverse           !! Reverse normal vector
    real(kind=real64)      :: x_i(3)            !! Collocation point position vector
    complex(kind=real64)      :: mu                !! Shear modulus
    complex(kind=real64)      :: nu                !! Poisson's ratio
    complex(kind=real64)      :: h(e%n_pnodes,3,3) !! h integration kernels matrix
    complex(kind=real64)      :: g(e%n_snodes,3,3) !! g integration kernels matrix
    ! Local
    integer           :: il, ik                     ! Counter for load / observation components
    integer           :: kip                        ! Counter variable for integration points loop
    real(kind=real64) :: x(3)                       ! Position vector at integration point
    real(kind=real64) :: n(3)                       ! Unit normal vector at integration point
    real(kind=real64) :: pphijw(e%n_pnodes)         ! phi^p * jacobian * weight at integration point
    real(kind=real64) :: sphijw(e%n_snodes)         ! phi^s * jacobian * weight at integration point
    real(kind=real64) :: rv(3)                      ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64) :: r, d1r, d1r2               ! Distance vector module and its inverse
    real(kind=real64) :: drdx(3)                    ! Distance vector derivatives with respect to x_k
    real(kind=real64) :: drdn                       ! Partial derivative of r respect to unit normal
    complex(kind=real64) :: cteu1, cteu2, ctet1, ctet2 ! Auxiliary constants
    complex(kind=real64) :: fs_u, fs_t                 ! Fundamental solutions values
    ! Initialize kernels
    h=0.d0
    g=0.d0
    ! Initialize auxiliary constants for fundamental solutions calculation
    cteu1=1.d0/(16.d0*c_pi*mu*(1.d0-nu))
    cteu2=3.d0-4.d0*nu
    ctet1=-1.d0/(8.d0*c_pi*(1.d0-nu))
    ctet2=1.d0-2.d0*nu
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      x=e%ps_x(:,kip,ps)
      n=e%ps_n(:,kip,ps)
      pphijw=e%ps_pphijw(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      d1r=1.d0/r
      d1r2=d1r**2
      drdx=rv*d1r
      drdn=dot_product(drdx,n)
      do il=1,3
        do ik=1,3
          fs_u=d1r*(cteu2*c_dkr(il,ik)+drdx(il)*drdx(ik))
          fs_t=d1r2*((ctet2*c_dkr(il,ik)+3.d0*drdx(il)*drdx(ik))*drdn+ctet2*(n(il)*drdx(ik)-n(ik)*drdx(il)))
          h(:,il,ik)=h(:,il,ik)+fs_t*pphijw
          g(:,il,ik)=g(:,il,ik)+fs_u*sphijw
        end do
      end do
    end do ! Loop through integrations points
    ! Multiply by constants
    h=ctet1*h
    g=cteu1*g
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_staela3d_sbie_ext_pre_cmplx

  !! This subroutine calculates the kernels for SBIE exterior integration (near collocation points) using Telles transformation
  !! within a subdivision of the element, needing only needs basic data.
  subroutine fbem_bem_staela3d_sbie_ext_st_cmplx(e,reverse,xi_s,x_i,barxip,barr,mu,nu,gln,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)       :: e                                !! Integration element
    logical                      :: reverse                          !! Reverse normal vector
    real(kind=real64)            :: xi_s(2,fbem_n_vertices(e%gtype)) !! Coordinates of the subdivision of the parent element
    real(kind=real64)            :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)            :: barxip(2)                        !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)            :: barr                             !! Telles jacobian at barxip
    complex(kind=real64)            :: mu                               !! Shear modulus
    complex(kind=real64)            :: nu                               !! Poisson's ratio
    integer                      :: gln                              !! 1D Gauss-Legendre number of integration points (<=32)
    complex(kind=real64)            :: h(e%n_pnodes,3,3)                !! h kernel vector
    complex(kind=real64)            :: g(e%n_snodes,3,3)                !! g kernel vector
    ! Local
    integer                      :: kphi                       ! Counter variable for shape functions loops
    integer                      :: k1                         ! Counter variable for reference coordinate xi_1
    integer                      :: k2                         ! Counter variable for reference coordinate xi_2
    real(kind=real64)            :: aux(10)                    ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)           ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi1(e%n_gnodes)      ! Geometrical shape functions derivatives values
    real(kind=real64)            :: dgphidxi2(e%n_gnodes)      ! Geometrical shape functions derivatives values
    real(kind=real64)            :: pphi(e%n_pnodes)           ! Functional shape functions values
    real(kind=real64)            :: sphi(e%n_snodes)           ! Functional shape functions values
    real(kind=real64)            :: gamma(2)                   ! Vector of gamma_1,gamma_2 coordinates
    real(kind=real64)            :: w(2)                       ! Weights of the integration rule
    real(kind=real64)            :: xip(2)                     ! Vector of xip_1,xip_2 coordinates
    real(kind=real64)            :: dxidxi1p(2)                ! xi derivatives with respect to xi1p
    real(kind=real64)            :: dxidxi2p(2)                ! xi derivatives with respect to xi2p
    real(kind=real64)            :: js                         ! Subdivision jacobian
    real(kind=real64)            :: xi(2)                      ! Vector of xi_1,xi_2 coordinates
    real(kind=real64)            :: xipp(2)                    ! Coordinate xipp used for quadrilateral-triangle transformation
    real(kind=real64)            :: barxipp(2)                 ! Coordinate xipp of collocation point
    real(kind=real64)            :: jqt                        ! Jacobian of the quadrilateral-triangle transformation
    type(fbem_telles_parameters) :: telles_parameters(2)       ! Telles parameters for each coordinate
    real(kind=real64)            :: jt(2)                      ! Telles jacobian for each coordinate: xi_1->gamma_1 and xi_2->gamma_2
    real(kind=real64)            :: x(3)                       ! Position vector at xi_1,xi_2
    real(kind=real64)            :: T1(3), T2(3)               ! Tangent vectors at xi_1,xi_2
    real(kind=real64)            :: N(3)                       ! Normal vector at xi_1,xi_2
    real(kind=real64)            :: rv(3)                      ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)            :: r, d1r, d1r2               ! Distance vector module and its inverse
    real(kind=real64)            :: drdx(3)                    ! Distance vector derivatives with respect to x_k
    real(kind=real64)            :: jg                         ! Geometric jacobian
    real(kind=real64)            :: jw                         ! Jacobians * weights
    real(kind=real64)            :: drdn                       ! Partial derivative of r respect to unit normal
    real(kind=real64)            :: pphijw(e%n_pnodes)         ! Functional shape functions * jw
    real(kind=real64)            :: sphijw(e%n_snodes)         ! Functional shape functions * jw
    integer                      :: il, ik                     ! Counter for load / observation components
    complex(kind=real64)            :: cteu1, cteu2, ctet1, ctet2 ! Auxiliary constants
    complex(kind=real64)            :: fs_u, fs_t                 ! Fundamental solutions values
    ! Initialize kernels
    h=0.d0
    g=0.d0
    ! Initialize auxiliary constants for fundamental solutions calculation
    cteu1=1.d0/(16.d0*c_pi*mu*(1.d0-nu))
    cteu2=3.d0-4.d0*nu
    ctet1=-1.d0/(8.d0*c_pi*(1.d0-nu))
    ctet2=1.d0-2.d0*nu
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
            d1r=1.d0/r
            d1r2=d1r**2
            drdx=rv*d1r
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
            ! Add integration points
            do il=1,3
              do ik=1,3
                fs_u=d1r*(cteu2*c_dkr(il,ik)+drdx(il)*drdx(ik))
                fs_t=d1r2*((ctet2*c_dkr(il,ik)+3.d0*drdx(il)*drdx(ik))*drdn+ctet2*(n(il)*drdx(ik)-n(ik)*drdx(il)))
                h(:,il,ik)=h(:,il,ik)+fs_t*pphijw
                g(:,il,ik)=g(:,il,ik)+fs_u*sphijw
              end do
            end do
          end do
        end do
      ! TRIANGULAR ELEMENTS
      case (3)
        ! Telles transformation is applied to Gauss-Legendre*Gauss-Legendre quadrature before the quadrilateral->triangle
        ! transformation. Because barxip for triangles are given in area triangle coordinates, for Telles transformation
        ! they must be transformed to quadrilateral coordinates. A special treatment is needed when barxi_2 is near 1, because
        ! transformation diverges.
        ! BARXIP->BARXIPP TRANSFORMATION
        if (barxip(2).gt.0.995d0) then
          barxipp(1)=0.5d0
          barxipp(2)=1.0d0
        else
          barxipp(1)=barxip(1)/(1.d0-barxip(2))
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
            d1r=1.d0/r
            d1r2=d1r**2
            drdx=rv*d1r
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
            ! Add integration points
            do il=1,3
              do ik=1,3
                fs_u=d1r*(cteu2*c_dkr(il,ik)+drdx(il)*drdx(ik))
                fs_t=d1r2*((ctet2*c_dkr(il,ik)+3.d0*drdx(il)*drdx(ik))*drdn+ctet2*(n(il)*drdx(ik)-n(ik)*drdx(il)))
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
    h=ctet1*h
    g=cteu1*g
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_staela3d_sbie_ext_st_cmplx

  !! This subroutine calculates adaptatively the kernels for SBIE exterior integration using Telles transformation and subdivision
  !! if needed, using only basic data.
  recursive subroutine fbem_bem_staela3d_sbie_ext_adp_cmplx(e,reverse,xi_s,x_i,mu,nu,qsp,ks,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                !! Element
    logical                  :: reverse                          !! Reverse orientation
    real(kind=real64)        :: xi_s(2,fbem_n_vertices(e%gtype)) !! Subdivision of the parent element
    real(kind=real64)        :: x_i(3)                           !! Collocation point position vector
    complex(kind=real64)        :: mu                               !! Shear modulus
    complex(kind=real64)        :: nu                               !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp                              !! Quasi-singular integration parameters
    integer                  :: ks                               !! Current level of subdivisions
    integer                  :: ns                               !! Maximum level of subdivision
    complex(kind=real64)        :: h(e%n_pnodes,3,3)                !! h integration kernels matrix
    complex(kind=real64)        :: g(e%n_snodes,3,3)                !! g integration kernels matrix
    ! Local
    integer           :: gln_near                             ! 1D Gauss-Legendre integ. points required to integrate only  the quasi-singular integrand
    integer           :: gln                                  ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
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
    complex(kind=real64) :: h_tmp(e%n_pnodes,3,3)                ! h integration kernels matrix (temporary)
    complex(kind=real64) :: g_tmp(e%n_snodes,3,3)                ! g integration kernels matrix (temporary)
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
      cl=fbem_characteristic_length(3,e%gtype,x_s,1.d-12)
      call fbem_nearest_element_point_bem(3,e%gtype,x_s,cl,x_i,barxip,rmin,d,method)
    end if
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,5,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_staela3d_sbie_ext_adp',ns,'maximum number of subdivisions reached')
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
          call fbem_bem_staela3d_sbie_ext_adp_cmplx(e,reverse,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,h,g)
          ! SUBTRI 2
          tmp_xi_s(:,1)=xi_s(:,2)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          call fbem_bem_staela3d_sbie_ext_adp_cmplx(e,reverse,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,h,g)
          ! SUBTRI 3
          tmp_xi_s(:,1)=xi_s(:,3)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          call fbem_bem_staela3d_sbie_ext_adp_cmplx(e,reverse,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,h,g)
          ! SUBTRI 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          call fbem_bem_staela3d_sbie_ext_adp_cmplx(e,reverse,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,h,g)
        ! QUADRILATERALS
        case (4)
          ! SUBQUAD 1
          tmp_xi_s(:,1)=xi_s(:,1)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,3)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          call fbem_bem_staela3d_sbie_ext_adp_cmplx(e,reverse,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,h,g)
          ! SUBQUAD 2
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=xi_s(:,2)
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,4)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          call fbem_bem_staela3d_sbie_ext_adp_cmplx(e,reverse,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,h,g)
          ! SUBQUAD 3
          tmp_xi_s(:,1)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=xi_s(:,3)
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          call fbem_bem_staela3d_sbie_ext_adp_cmplx(e,reverse,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,h,g)
          ! SUBQUAD 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          tmp_xi_s(:,2)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=xi_s(:,4)
          call fbem_bem_staela3d_sbie_ext_adp_cmplx(e,reverse,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,h,g)
      end select
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_staela3d_sbie_ext_st_cmplx(e,reverse,xi_s,x_i,barxip,barr,mu,nu,gln,h_tmp,g_tmp)
      h=h+h_tmp
      g=g+g_tmp
    end if
  end subroutine fbem_bem_staela3d_sbie_ext_adp_cmplx

  !! This subroutine calculates static elastic integration kernels vectors h and g of singular formulation for interior
  !! integration.
  subroutine fbem_bem_staela3d_sbie_int_cmplx(type_g,type_f1,type_f2,delta_f,x_nodes,reverse,xi_i,mu,nu,h,g)
    implicit none
    ! I/O
    integer           :: type_g                          !! Geometrial interpolation
    integer           :: type_f1                         !! Functional interpolation (primary variables)
    integer           :: type_f2                         !! Functional interpolation (secondary variables)
    real(kind=real64) :: delta_f                         !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64) :: x_nodes(3,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
    logical           :: reverse                         !! Normal vector inversion
    real(kind=real64) :: xi_i(2)                         !! Reference coordinates of the singular point.
    complex(kind=real64) :: mu                              !! Shear modulus
    complex(kind=real64) :: nu                              !! Poisson's ratio
    complex(kind=real64) :: h(fbem_n_nodes(type_f1),3,3)    !! h kernel vector
    complex(kind=real64) :: g(fbem_n_nodes(type_f2),3,3)    !! g kernel vector
    ! Local
    integer           :: ksubtri                          ! Counter variable for subtriangles loop
    integer           :: nsubtri                          ! Number of subtriangles
    integer           :: subtriangle(8)                   ! Vector that contains what subtriangles need to be integrated
    real(kind=real64) :: theta_subtri(2,8)                ! Matrix that contains the angles theta of the subtriangles to be integrated
    real(kind=real64) :: thetap_subtri(2,8)               ! Matrix that contains the angles thetap of the subtriangles to be integrated
    integer           :: ktheta                           ! Counter variable for theta coordinate loop
    integer           :: krho                             ! Counter variable for rho coordinate loop
    integer           :: kphi                             ! Counter coordinates for shape functions loops
    integer           :: nnodes_g                         ! Number of nodes of geometrical interpolation
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
    real(kind=real64) :: phi_f1(fbem_n_nodes(type_f1))    ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phi_f2(fbem_n_nodes(type_f2))    ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phi_f1_i(fbem_n_nodes(type_f1))  ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phi_g(fbem_n_nodes(type_g))      ! Geometrical shape functions values at xi_1,xi_2
    real(kind=real64) :: dphidxi1_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64) :: dphidxi2_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64) :: x_i(3)                           ! Collocation point position vector
    real(kind=real64) :: x(3)                             ! Position vector at xi_1,xi_2
    real(kind=real64) :: T1(3), T2(3)                     ! Tangent vectors at xi_1,xi_2
    real(kind=real64) :: N(3)                             ! Normal vector at xi_1,xi_2
    real(kind=real64) :: rv(3)                            ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64) :: r, d1r, d1r2                     ! Distance vector module and its inverse
    real(kind=real64) :: drdx(3)                          ! Distance vector derivatives with respect to x_k
    real(kind=real64) :: jg                               ! Geometric jacobian
    real(kind=real64) :: jw                               ! Jacobian * weights
    real(kind=real64) :: drdn                             ! Partial derivative of r respect to unit normal
    real(kind=real64) :: phif1jw(fbem_n_nodes(type_f1))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phif2jw(fbem_n_nodes(type_f2))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: costheta, sintheta               ! cos(theta), sin(theta)
    integer           :: il, ik                           ! Counter for load / observation components
    complex(kind=real64) :: cteu1, cteu2, ctet1, ctet2       ! Auxiliary constants
    complex(kind=real64) :: fs_u, fs_t                       ! Fundamental solutions values
    ! Local variables associated with line integrals
    real(kind=real64)              :: re                                   ! Quasi-singular integration relative error
    integer                        :: gln_min                              ! Minimum Number of 1D Gauss-Legendre number of points
    integer                        :: ns_max                               ! Maximum number of subdivisions
    type(fbem_qs_parameters)       :: qsp                                  ! Quasi-singular integration parameters
    real(kind=real64)              :: xi_s(1,2)                            ! Edge subdivision
    integer                        :: kedge                                ! Counter variable for edges
    integer                        :: nedges                               ! Number of edges
    logical                        :: integrate                            ! Control variable
    integer                        :: type_edge                            ! Line element type for line integrals
    integer                        :: nnodes_edge                          ! Number of nodes of the edge
    real(kind=real64), allocatable :: x_nodes_edge(:,:)                    ! Coordinates of the edge elements
    real(kind=real64)              :: hli(3,3)                             ! Line integrals values
    !
    ! Initialization
    !
    ! Kernel vectors
    h=0.d0
    g=0.d0
    ! Number of nodes of geometrical interpolation
    nnodes_g=fbem_n_nodes(type_g)
    ! Initialize auxiliary constants for fundamental solutions calculation
    cteu1=1.d0/(16.d0*c_pi*mu*(1.d0-nu))
    cteu2=3.d0-4.d0*nu
    ctet1=-1.d0/(8.d0*c_pi*(1.d0-nu))
    ctet2=1.d0-2.d0*nu
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
          jg=sqrt(dot_product(N,N))
          ! Unit normal vector
          n=N/jg
          ! Distance vector
          rv=x-x_i
          ! Distance vector norm
          r=dot_product(rv,rv)
          r=sqrt(r)
          d1r=1.d0/r
          d1r2=d1r**2
          ! r_{,k}
          drdx=rv*d1r
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
          ! Loop through load direction and observation direction
          do il=1,3
            do ik=1,3
              ! Fundamental solutions values without cteu1 and ctep1, respectively
              fs_u=d1r*(cteu2*c_dkr(il,ik)+drdx(il)*drdx(ik))
              ! Regular part of the fundamental solution
              fs_t=d1r2*drdn*(ctet2*c_dkr(il,ik)+3.d0*drdx(il)*drdx(ik))
              ! Add to kernels
              h(:,il,ik)=h(:,il,ik)+fs_t*phif1jw
              g(:,il,ik)=g(:,il,ik)+fs_u*phif2jw
              ! The regular part of the CPV integral of h
              fs_t=d1r2*ctet2*(n(il)*drdx(ik)-n(ik)*drdx(il))
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
        h(:,il,ik)=h(:,il,ik)+phi_f1_i*ctet2*hli(il,ik)
      end do
    end do
    ! Multiply h by ctet1 and g by cteu1
    h=ctet1*h
    g=cteu1*g
    ! If the normal has to be reversed, then h=-h
    if (reverse) h=-h
  end subroutine fbem_bem_staela3d_sbie_int_cmplx

  !! This subroutine calculates adaptatively the line integral epsilon_{ijk}*e_k·t/r for the SBIE interior integration using Telles
  !! transformation and subdivision. Before calling this subroutine, mli1 must be set to zero.
  recursive subroutine fbem_bem_staela3d_sbie_int_li(gtype,x_nodes,xi_s,x_i,ngp_min,qsp,ks,ns,hli)
    implicit none
    ! I/O
    integer                      :: gtype                           !! Geometrial interpolation
    real(kind=real64)            :: x_nodes(3,fbem_n_nodes(gtype))  !! Position vectors of geometrical nodes
    real(kind=real64)            :: xi_s(1,2)                       !! Subdivision of the parent element
    real(kind=real64)            :: x_i(3)                          !! Collocation point position vector
    integer                      :: ngp_min                         !! Minimum number of Gauss points
    type(fbem_qs_parameters)     :: qsp                             !! Quasi-singular integration parameters
    integer                      :: ks                              !! Current level of subdivisions
    integer                      :: ns                              !! Maximum level of subdivision
    real(kind=real64)            :: hli(3,3)                        !! Line integral epsilon_{ijk}*e_k·t/r
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
    real(kind=real64)            :: hli_tmp(3,3)                    ! Line integral epsilon_{ijk}*e_k·t/r (temporary)
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
    real(kind=real64)            :: x(3), T(3), rv(3), r, dr1       ! Position, tangent, distance vector, distance and 1/r at xi
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
        call fbem_warning_message(error_unit,0,'fbem_bem_staela3d_sbie_int_li',ns,'maximum number of subdivisions reached')
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
      call fbem_bem_staela3d_sbie_int_li(gtype,x_nodes,tmp_xi_s,x_i,ngp_min,qsp,ks+1,ns,hli)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_staela3d_sbie_int_li(gtype,x_nodes,tmp_xi_s,x_i,ngp_min,qsp,ks+1,ns,hli)
    ! Calculate the subdivided element using Telles transformation
    else
      ! Initialization
      n_gnodes=fbem_n_nodes(gtype)
      hli_tmp=0.d0
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
        dr1=1.0d0/r
        ! Jacobians * weight
        jw=jg*js*jt*w
        ! Add integration points
        hli_tmp(1,2)=hli_tmp(1,2)-dr1*t(3)*jw
        hli_tmp(1,3)=hli_tmp(1,3)+dr1*t(2)*jw
        hli_tmp(2,3)=hli_tmp(2,3)-dr1*t(1)*jw
      end do
      ! Copy the antisymmetrical parts
      hli_tmp(2,1)=-hli_tmp(1,2)
      hli_tmp(3,1)=-hli_tmp(1,3)
      hli_tmp(3,2)=-hli_tmp(2,3)
      ! Add hli
      hli=hli+hli_tmp
    end if
  end subroutine fbem_bem_staela3d_sbie_int_li

!  !! This subroutine calculates adaptatively the line integrals epsilon_{ijk}*e_k·t/r for the SBIE interior integration using Telles
!  !! transformation and subdivision. Before calling this subroutine, mli1 must be set to zero, and xi1=-1.0d0 and xi2=1.0d0
!  recursive subroutine fbem_bem_staela3d_sbie_int_li(type_g,x_nodes,xi1,xi2,x_i,ngp_min,re,hli)
!    implicit none
!    ! I/O
!    integer           :: type_g                          !! Geometrial interpolation
!    real(kind=real64) :: x_nodes(3,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
!    real(kind=real64) :: xi1                             !! Minimum coordinate xi of the subdivision
!    real(kind=real64) :: xi2                             !! Maximum coordinate xi of the subdivision
!    real(kind=real64) :: x_i(3)                          !! Collocation point position vector
!    integer           :: ngp_min                         !! Minimum number of Gauss points
!    real(kind=real64) :: re                              !! Integration relative error
!    real(kind=real64) :: hli(3,3)                        !! Line integrals e_k·t/r
!    ! Local
!    real(kind=real64) :: xise(2,4)                       ! xi coordinates of the subdivision
!    integer           :: ngp                             ! Number of Gauss points
!    logical           :: achieve                         ! True if the calculated ngp achieve the relative error
!    real(kind=real64) :: barxi                           ! Nearest element coordinate with respect to collocation point
!    real(kind=real64) :: r_min                           ! Minimum distance between collocation point and barxi on the element
!    real(kind=real64) :: barr                            ! Telles jacobian at barxi
!    real(kind=real64) :: csize                           ! Characteristic size
!    real(kind=real64) :: d                               ! Normalized distance between collocation point and element subdivision
!    real(kind=real64) :: barxisub                        ! Nearest subdivision coordinate with respect to the collocation point
!    ! Local calculation
!    real(kind=real64)            :: hli_tmp(3,3)                    ! Line integrals e_k·t/r (temporary)
!    integer                      :: kphi                            ! Counter variable for shape functions loops
!    integer                      :: nnodes_g                        ! Number of nodes of the geometrical element
!    integer                      :: kip                             ! Counter variable for integration points
!    real(kind=real64)            :: gamma                           ! Reference coordinate
!    type(fbem_telles_parameters) :: telles_parameters               ! Telles parameters
!    real(kind=real64)            :: jt                              ! Telles jacobian for each coordinate: xi->gamma
!    real(kind=real64)            :: xip                             ! Reference coordinate of the subdivision
!    real(kind=real64)            :: barxip                          ! Collocation point in xip space
!    real(kind=real64)            :: js                              ! Jacobian of xip->xi transformation
!    real(kind=real64)            :: xi                              ! Reference coordinate
!    real(kind=real64)            :: w                               ! Weights of each integration point
!    real(kind=real64)            :: aux(10)                         ! Auxiliary variable needed for shape_functions module resources
!    real(kind=real64)            :: phi_g(fbem_n_nodes(type_g))     ! Geometrical shape functions values
!    real(kind=real64)            :: dphidxi_g(fbem_n_nodes(type_g)) ! Geometrical shape functions first derivatives values
!    real(kind=real64)            :: x(3)                            ! Position vector at xi
!    real(kind=real64)            :: T(3)                            ! Tangent vector at xi
!    real(kind=real64)            :: rv(3)                           ! Distance vector between collocation point and integration point (x-x_i)
!    real(kind=real64)            :: r, d1r                          ! Distance vector module and its inverse
!    real(kind=real64)            :: jw                              ! Jacobians * weight
!    real(kind=real64)            :: jg                              ! Geometric jacobian
!    !
!    ! Integration setup
!    !
!    ! Find the minimum distance between the element subdivision and the collocation point, and the nearest element coordinate in
!    ! element reference space [-1,1]
!    call fbem_nearest_xi_subdivision(3,type_g,x_nodes,xi1,xi2,x_i,25,14,2,barxi,r_min)
!    ! Obtain the characteristic size of the element subdivision
!    xise(1,1)=xi1
!    xise(1,2)=xi2
!    csize=fbem_characteristic_length_subdivision(3,type_g,x_nodes,xise,1.d-6)
!    ! Calculate the normalized minimum distance
!    d=r_min/csize
!    ! Obtain the barxi but in the element subdivision space (barxisub) [-1,1]
!    barxisub=2.0d0*(barxi-xi1)/(xi2-xi1)-1.0d0
!    ! Obtain the estimated number of Gaussian points
!    call fbem_qs_f_ngp_telles_1d(type_g,fbem_f_1r1,barxisub,d,re,ngp,achieve)
!    !
!    ! If not achieve the indicated relative error, then subdivide by 1/2
!    !
!    if (achieve.eqv.(.false.)) then
!      call fbem_bem_staela3d_sbie_int_li(type_g,x_nodes,xi1,0.5d0*(xi1+xi2),x_i,ngp_min,re,hli)
!      call fbem_bem_staela3d_sbie_int_li(type_g,x_nodes,0.5d0*(xi1+xi2),xi2,x_i,ngp_min,re,hli)
!    !
!    ! If achieve, then integrate the present subdivision
!    !
!    else
!      ! Obtain the optimum Telles jacobian
!      barr=fbem_telles_barr(d,fbem_f_any)
!      ! Use ngp_min if ngp is less than ngp_min
!      if (ngp.lt.ngp_min) ngp=ngp_min
!      ! Integrals initialization
!      hli_tmp=0.d0
!      ! Number of nodes of geometrical interpolation
!      nnodes_g=fbem_n_nodes(type_g)
!      ! Subdivision jacobian (constant)
!      js=xi2-xi1
!      ! Barxi in subdivision reference space
!      barxip=(barxi-xi1)/js
!      ! Calculate Telles parameters
!      telles_parameters=fbem_telles01_calculate_parameters(barxip,barr)
!      ! Loop through gamma coordinate
!      do kip=1,gl01_n(ngp)
!        ! GAMMA->XIP->XI COORDINATE TRANSFORMATION
!        ! Gamma coordinate and weight
!        gamma=gl01_xi(kip,ngp)
!        w=gl01_w(kip,ngp)
!        ! xip coordinate, weight and jacobian from Telles transformation
!        call fbem_telles_xi_and_jacobian(telles_parameters,gamma,xip,jt)
!        ! xi coordinate
!        xi=js*xip+xi1
!        ! COMPONENTS OF THE FUNDAMENTAL SOLUTIONS
!        ! Geometrical shape functions and first derivatives at xi
!#       define etype type_g
!#       define delta 0.0d0
!#       define phi phi_g
!#       define dphidxi dphidxi_g
!#       include <phi_and_dphidxi_1d.rc>
!#       undef etype
!#       undef delta
!#       undef phi
!#       undef dphidxi
!        ! Components calculation of x and T at xi
!        x=0.d0
!        T=0.d0
!        do kphi=1,nnodes_g
!          x=x+phi_g(kphi)*x_nodes(:,kphi)
!          T=T+dphidxi_g(kphi)*x_nodes(:,kphi)
!        end do
!        ! Geometric jacobian
!        jg=dot_product(T,T)
!        jg=sqrt(jg)
!        ! Unit tangent
!        t=T/jg
!        ! Distance vector
!        rv=x-x_i
!        ! Distance vector norm
!        r=dot_product(rv,rv)
!        r=sqrt(r)
!        d1r=1.d0/r
!        ! Jacobians * weight
!        jw=jg*js*jt*w
!        ! Add integration points
!        hli_tmp(1,2)=hli_tmp(1,2)-d1r*t(3)*jw
!        hli_tmp(1,3)=hli_tmp(1,3)+d1r*t(2)*jw
!        hli_tmp(2,3)=hli_tmp(2,3)-d1r*t(1)*jw
!      end do ! Loop through gamma coordinate
!      ! Copy the antisymmetrical parts
!      hli_tmp(2,1)=-hli_tmp(1,2)
!      hli_tmp(3,1)=-hli_tmp(1,3)
!      hli_tmp(3,2)=-hli_tmp(2,3)
!      ! Add hli
!      hli=hli+hli_tmp
!    end if
!  end subroutine fbem_bem_staela3d_sbie_int_li

  ! =====================
  ! BE BODY LOAD ELEMENTS
  ! =====================

  !! Efficient automatic integration of BE body load elements
  subroutine fbem_bem_staela3d_sbie_bl_auto(e,x_i,mu,nu,qsp,ns,gbl)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                   !! Integration element
    real(kind=real64)        :: x_i(3)              !! Collocation point
    real(kind=real64)        :: mu                  !! Shear modulus
    real(kind=real64)        :: nu                  !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp                 !! Quasi-singular integration parameters
    integer                  :: ns                  !! Maximum level of subdivisions
    real(kind=real64)        :: gbl(e%n_snodes,3,3) !! gbl integration kernel
    ! Local
    real(kind=real64)        :: r(3)                               ! Distance vector
    real(kind=real64)        :: rmin                               ! Minimum distance between element and x_i
    real(kind=real64)        :: barxi(e%d)                         ! Local coordinates of the nearest element point with respect to x_i
    real(kind=real64)        :: d                                  ! Dimensionless distance
    integer                  :: delta                              ! Control variable
    real(kind=real64)        :: xi_s(e%d,fbem_n_vertices(e%gtype)) ! Local coordinates of the element subdivision
    integer                  :: method                             ! Method used when calculating the nearest element point
    integer                  :: gln_near                           ! 1D Gauss-Legendre integration points required by the quasi-singular function
    integer                  :: gln                                ! 1D Gauss-Legendre integration points used in the integration
    integer                  :: ps                                 ! Selected precalculated dataset
    integer                  :: i                                  ! Counter
    ! POINT BODY LOAD
    if (e%d.eq.0) then
      r=e%x(:,1)-x_i
      rmin=sqrt(dot_product(r,r))
      if (rmin.eq.0.d0) then
        call fbem_error_message(output_unit,0,'fbem_bem_staela3d_sbie_bl_auto',0,'it is not possible to collocate at a point load')
      else
        call fbem_bem_staela3d_sbie_u(e%x(:,1),x_i,mu,nu,gbl(1,:,:))
        return
      end if
    ! LINE, SURFACE OR VOLUME BODY LOAD
    ! Determine if interior or exterior integration
    !   - Interior integration (delta=1) requires: xi_i
    !   - Exterior integration (delta=0) requires: x_i, barxi, rmin and d
    ! Use the element ball
    else
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
    end if
    ! Integrate
    select case (delta)
      case (1)
        call fbem_bem_staela3d_sbie_bl_int(e,barxi,mu,nu,qsp,ns,gbl)
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
          call fbem_bem_staela3d_sbie_bl_ext_pre(ps,e,x_i,mu,nu,gbl)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_staela3d_sbie_bl_ext_adp(e,xi_s,x_i,mu,nu,qsp,1,ns,gbl)
        end if
    end select
  end subroutine fbem_bem_staela3d_sbie_bl_auto

  !! This subroutine calculates the kernels for SBIE exterior integration needing precalculated data at integration points.
  subroutine fbem_bem_staela3d_sbie_bl_ext_pre(ps,e,x_i,mu,nu,gbl)
    implicit none
    ! I/O
    integer                :: ps                  !! Selected precalculated dataset
    type(fbem_bem_element) :: e                   !! Element
    real(kind=real64)      :: x_i(3)              !! Collocation point position vector
    real(kind=real64)      :: mu                  !! Shear modulus
    real(kind=real64)      :: nu                  !! Poisson's ratio
    real(kind=real64)      :: gbl(e%n_snodes,3,3) !! gbl kernels
    ! Local
    integer                :: il, ik              ! Counter for load / observation components
    integer                :: kip                 ! Counter variable for integration points loop
    real(kind=real64)      :: x(3)                ! Position vector at integration point
    real(kind=real64)      :: sphijw(e%n_snodes)  ! phi^s * jacobian * weight at integration point
    real(kind=real64)      :: rv(3)               ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)      :: r, d1r              ! Distance vector module and its inverse
    real(kind=real64)      :: drdx(3)             ! Distance vector derivatives with respect to x_k
    real(kind=real64)      :: cteu1, cteu2        ! Auxiliary constants
    real(kind=real64)      :: fs_u                ! Fundamental solutions values
    ! Initialize kernels
    gbl=0.d0
    ! Initialize auxiliary constants for fundamental solutions calculation
    cteu1=1.d0/(16.d0*c_pi*mu*(1.d0-nu))
    cteu2=3.d0-4.d0*nu
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      x=e%ps_x(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      d1r=1.d0/r
      drdx=rv*d1r
      do il=1,3
        do ik=1,3
          fs_u=d1r*(cteu2*c_dkr(il,ik)+drdx(il)*drdx(ik))
          gbl(:,il,ik)=gbl(:,il,ik)+fs_u*sphijw
        end do
      end do
    end do ! Loop through integrations points
    ! Multiply by constants
    gbl=cteu1*gbl
  end subroutine fbem_bem_staela3d_sbie_bl_ext_pre

  !! This subroutine calculates adaptatively the body load kernels for SBIE exterior integration using Telles transformation and
  !! subdivision if needed.
  recursive subroutine fbem_bem_staela3d_sbie_bl_ext_adp(e,xi_s,x_i,mu,nu,qsp,ks,ns,gbl)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                      !! Element
    real(kind=real64)        :: xi_s(e%d,fbem_n_vertices(e%gtype))     !! Subdivision of the parent element
    real(kind=real64)        :: x_i(3)                                 !! Collocation point position vector
    real(kind=real64)        :: mu                                     !! Shear modulus
    real(kind=real64)        :: nu                                     !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp                                    !! Quasi-singular integration parameters
    integer                  :: ks                                     !! Current level of subdivisions
    integer                  :: ns                                     !! Maximum level of subdivision
    real(kind=real64)        :: gbl(e%n_snodes,3,3)                    !! gbl kernels
    ! Local
    integer                  :: gln_near                               ! 1D Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer                  :: gln                                    ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical                  :: subdivide                              ! True if subdivision has to be performed
    real(kind=real64)        :: barxi(e%d)                             ! Nearest element coordinate with respect to collocation point
    real(kind=real64)        :: barxip(e%d)                            ! Nearest element subdivision local coordinate with respect to collocation point
    real(kind=real64)        :: rmin                                   ! Minimum distance between collocation point and barxi on the element
    real(kind=real64)        :: barr                                   ! Telles jacobian at barxi
    real(kind=real64)        :: cl                                     ! Characteristic length
    real(kind=real64)        :: d                                      ! Normalized distance between collocation point and element subdivision
    integer                  :: method                                 ! Method used in nearest point algorithm
    real(kind=real64)        :: tmp_xi_s(e%d,fbem_n_vertices(e%gtype)) ! Subdivision
    real(kind=real64)        :: x_s(3,e%n_gnodes)                      ! Coordinates of the element subdivision
    real(kind=real64)        :: gbl_tmp(e%n_snodes,3,3)                ! gbl kernels (temporary)
    ! Initialize
    if (ks.eq.1) then
      gbl=0.d0
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
        call fbem_warning_message(error_unit,0,'fbem_bem_staela3d_sbie_bl_ext_adp',ns,'maximum number of subdivisions reached')
        gln_near=30
      end if
    else
      if (gln_near.eq.0) subdivide=.true.
    end if
    ! Subdivide
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
          call fbem_bem_staela3d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
          ! SUBLINE 2
          tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
          tmp_xi_s(1,2)=xi_s(1,2)
          call fbem_bem_staela3d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
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
              call fbem_bem_staela3d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
              ! SUBTRI 2
              tmp_xi_s(:,1)=xi_s(:,2)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              call fbem_bem_staela3d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
              ! SUBTRI 3
              tmp_xi_s(:,1)=xi_s(:,3)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              call fbem_bem_staela3d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
              ! SUBTRI 4
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
              call fbem_bem_staela3d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
            ! QUADRILATERALS
            case (4)
              ! SUBQUAD 1
              tmp_xi_s(:,1)=xi_s(:,1)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,3)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,4)=0.50d0*(xi_s(:,1)+xi_s(:,4))
              call fbem_bem_staela3d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
              ! SUBQUAD 2
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,2)=xi_s(:,2)
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,4)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              call fbem_bem_staela3d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
              ! SUBQUAD 3
              tmp_xi_s(:,1)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=xi_s(:,3)
              tmp_xi_s(:,4)=0.50d0*(xi_s(:,3)+xi_s(:,4))
              call fbem_bem_staela3d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
              ! SUBQUAD 4
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,4))
              tmp_xi_s(:,2)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,4)=xi_s(:,4)
              call fbem_bem_staela3d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
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
      call fbem_bem_staela3d_sbie_bl_ext_st(e,xi_s,x_i,barxip,barr,mu,nu,gln,gbl_tmp)
      gbl=gbl+gbl_tmp
    end if
  end subroutine fbem_bem_staela3d_sbie_bl_ext_adp

  !! This subroutine calculates the body load kernels for SBIE exterior integration (near collocation points) using Telles
  !! transformation within a subdivision of the element, needing only basic data.
  subroutine fbem_bem_staela3d_sbie_bl_ext_st(e,xi_s,x_i,barxip,barr,mu,nu,gln,gbl)
    implicit none
    ! I/O
    type(fbem_bem_element)       :: e                                  !! Integration element
    real(kind=real64)            :: xi_s(e%d,fbem_n_vertices(e%gtype)) !! Coordinates of the subdivision of the parent element (xi space [-1,1])
    real(kind=real64)            :: x_i(3)                             !! Collocation point position vector
    real(kind=real64)            :: barxip(e%d)                        !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)            :: barr                               !! Telles jacobian at barxip
    real(kind=real64)            :: mu                                 !! Shear modulus
    real(kind=real64)            :: nu                                 !! Poisson's ratio
    integer                      :: gln                                !! 1D Gauss-Legendre number of integration points (<=32)
    real(kind=real64)            :: gbl(e%n_snodes,3,3)                !! gbl kernels
    ! Local
    integer                      :: il, ik                             ! Counter variables for load direction and observation direction
    integer                      :: kphi                               ! Counter variable for shape functions loops
    integer                      :: kip                                ! Counter variable of integration points
    integer                      :: k1                                 ! Counter variable for reference coordinate xi_1
    integer                      :: k2                                 ! Counter variable for reference coordinate xi_2
    real(kind=real64)            :: gamma(e%d)                         ! Coordinate gamma (Telles transformation space [-1,1])
    real(kind=real64)            :: w(e%d)                             ! Weights of an integration point
    type(fbem_telles_parameters) :: telles_parameters(e%d)             ! Telles parameters
    real(kind=real64)            :: jt(e%d)                            ! Telles jacobian: xip->gamma
    real(kind=real64)            :: xip(e%d)                           ! Coordinate xip (subdivision space [-1,1])
    real(kind=real64)            :: js                                 ! Subdivision jacobian: xi->xip
    real(kind=real64)            :: xin(e%d)                           ! Coordinate xi [xi_s(1,1),xi_s(1,2)]
    real(kind=real64)            :: aux(10)                            ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)                   ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi(e%n_gnodes)               ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: dgphidxi1(e%n_gnodes)              ! Geometrical shape functions derivatives values
    real(kind=real64)            :: dgphidxi2(e%n_gnodes)              ! Geometrical shape functions derivatives values
    real(kind=real64)            :: dxidxi1p(2)                        ! xi derivatives with respect to xi1p
    real(kind=real64)            :: dxidxi2p(2)                        ! xi derivatives with respect to xi2p
    real(kind=real64)            :: xipp(2)                            ! Coordinate xipp used for quadrilateral-triangle transformation
    real(kind=real64)            :: barxipp(2)                         ! Coordinate xipp of collocation point
    real(kind=real64)            :: jqt                                ! Jacobian of the quadrilateral-triangle transformation
    real(kind=real64)            :: sphi(e%n_snodes)                   ! Functional shape functions values
    real(kind=real64)            :: x(3)                               ! Position vector at xi
    real(kind=real64)            :: N(3), T(3), T1(3), T2(3)           ! Normal and tangent vectors at xi
    real(kind=real64)            :: rv(3)                              ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)            :: r, d1r                             ! Distance vector module, its inverse and log(1/r)
    real(kind=real64)            :: drdx(3)                            ! Distance vector derivatives with respect to x_k
    real(kind=real64)            :: jg                                 ! Geometric jacobian
    real(kind=real64)            :: jw                                 ! Jacobians * weight
    real(kind=real64)            :: sphijw(e%n_snodes)                 ! secondary shape functions * jw
    real(kind=real64)            :: cteu1, cteu2                       ! Auxiliary constants
    real(kind=real64)            :: fs_u                               ! Fundamental solutions values
    ! Initialize kernels
    gbl=0.d0
    ! Initialize auxiliary constants for fundamental solutions calculation
    cteu1=1.d0/(16.d0*c_pi*mu*(1.d0-nu))
    cteu2=3.d0-4.d0*nu
    ! Calculate
    select case (e%d)
      !
      ! LINE LOAD
      !
      case (1)
        ! Calculate Telles parameters
        telles_parameters(1)=fbem_telles11_calculate_parameters(barxip(1),barr)
        ! Numerical integration
        do kip=1,gl11_n(gln)
          ! GAMMA COORDINATE
          gamma(1)=gl11_xi(kip,gln)
          w(1)=gl11_w(kip,gln)
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
          d1r=1.d0/r
          drdx=rv*d1r
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
          ! Add integration points
          do il=1,3
            do ik=1,3
              fs_u=d1r*(cteu2*c_dkr(il,ik)+drdx(il)*drdx(ik))
              gbl(:,il,ik)=gbl(:,il,ik)+fs_u*sphijw
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
#               define etype e%gtype
#               define delta 0.d0
#               define xi xin
#               define phi gphi
#               define dphidxi1 dgphidxi1
#               define dphidxi2 dgphidxi2
#               include <phi_and_dphidxik_2d.rc>
#               undef etype
#               undef delta
#               undef xi
#               undef phi
#               undef dphidxi1
#               undef dphidxi2
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
                ! xi->x jacobian
                jg=sqrt(dot_product(N,N))
                ! Distance vector and derived terms
                rv=x-x_i
                r=sqrt(dot_product(rv,rv))
                d1r=1.d0/r
                drdx=rv*d1r
                ! Jacobians * weights
                jw=jg*js*jt(1)*jt(2)*w(1)*w(2)
                ! FUNCTIONAL SHAPE FUNCTIONS
                ! Functional shape functions (body load variables) at xi
#               define etype e%stype
#               define delta e%stype_delta
#               define xi xin
#               define phi sphi
#               include <phi_2d.rc>
#               undef etype
#               undef delta
#               undef xi
#               undef phi
                ! Functional shape functions * jacobians * weights
                sphijw=sphi*jw
                ! Add integration points
                do il=1,3
                  do ik=1,3
                    fs_u=d1r*(cteu2*c_dkr(il,ik)+drdx(il)*drdx(ik))
                    gbl(:,il,ik)=gbl(:,il,ik)+fs_u*sphijw
                  end do
                end do
              end do
            end do
          ! TRIANGULAR ELEMENTS
          case (3)
            ! Telles transformation is applied to Gauss-Legendre*Gauss-Legendre quadrature before the quadrilateral->triangle
            ! transformation. Because barxip for triangles are given in area triangle coordinates, for Telles transformation
            ! they must be transformed to quadrilateral coordinates. A special treatment is needed when barxi_2 is near 1, because
            ! transformation diverges
            ! BARXIP->BARXIPP TRANSFORMATION
            if (barxip(2).gt.0.995d0) then
              barxipp(1)=0.5d0
              barxipp(2)=1.0d0
            else
              barxipp(1)=barxip(1)/(1.d0-barxip(2))
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
#               define delta 0.d0
#               define xi xin
#               define phi gphi
#               define dphidxi1 dgphidxi1
#               define dphidxi2 dgphidxi2
#               include <phi_and_dphidxik_2d.rc>
#               undef etype
#               undef delta
#               undef xi
#               undef phi
#               undef dphidxi1
#               undef dphidxi2
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
                ! xi->x jacobian
                jg=sqrt(dot_product(N,N))
                ! Distance vector and derived terms
                rv=x-x_i
                r=sqrt(dot_product(rv,rv))
                d1r=1.d0/r
                drdx=rv*d1r
                ! Jacobians * weights
                jw=jg*js*jqt*jt(1)*jt(2)*w(1)*w(2)
                ! FUNCTIONAL SHAPE FUNCTIONS
                ! Functional shape functions (body load variables) at xi
#               define etype e%stype
#               define delta e%stype_delta
#               define xi xin
#               define phi sphi
#               include <phi_2d.rc>
#               undef etype
#               undef delta
#               undef xi
#               undef phi
                ! Functional shape functions * jacobians * weights
                sphijw=sphi*jw
                ! Add integration points
                do il=1,3
                  do ik=1,3
                    fs_u=d1r*(cteu2*c_dkr(il,ik)+drdx(il)*drdx(ik))
                    gbl(:,il,ik)=gbl(:,il,ik)+fs_u*sphijw
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
        stop 'not yet : fbem_bem_staela3d_sbie_bl_ext_st'
      !
      ! OTHERS
      !
      case default
        call fbem_error_message(output_unit,0,'fbem_bem_staela3d_sbie_bl_int',0,&
                                'it is only possible to integrate line, surface or volume loads')
    end select
    ! Multiply by constants
    gbl=cteu1*gbl
  end subroutine fbem_bem_staela3d_sbie_bl_ext_st

  !! This subroutine calculates the body load kernels for SBIE interior integration, needing only basic data.
  subroutine fbem_bem_staela3d_sbie_bl_int(e,xi_i,mu,nu,qsp,ns,gbl)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                !! Integration element
    real(kind=real64)        :: xi_i(e%d)                        !! Local coordinate of the singular point.
    real(kind=real64)        :: mu                               !! Shear modulus
    real(kind=real64)        :: nu                               !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp                              !! Quasi-singular integration parameters
    integer                  :: ns                               !! Maximum level of subdivisions
    real(kind=real64)        :: gbl(e%n_snodes,3,3)              !! gbl kernels
    ! Local
    integer                :: il, ik                           ! Counter variables for load direction and observation direction
    integer                :: kphi                             ! Counter variable for shape functions loops
    real(kind=real64)      :: x_i(3), x_ic(3)                  ! Real coordinates of collocation point
    real(kind=real64)      :: ep1(3), ep2(3), ep3(3)           ! Local system at the collocation point
    real(kind=real64)      :: gphi(e%n_gnodes)                 ! Geometrical shape functions values at xi
    real(kind=real64)      :: dgphidxi1(e%n_gnodes)            ! Geometrical shape functions derivatives values at xi_1
    real(kind=real64)      :: dgphidxi2(e%n_gnodes)            ! Geometrical shape functions derivatives values at xi_2
    real(kind=real64)      :: sphi(e%n_snodes)                 ! Functional shape functions values at xi
    integer                :: ksubtri                          ! Counter variable for subtriangles loop
    integer                :: nsubtri                          ! Number of subtriangles
    integer                :: subtriangle(8)                   ! Vector that contains what subtriangles need to be integrated
    real(kind=real64)      :: theta_subtri(2,8)                ! Matrix that contains the angles theta of the subtriangles to be integrated
    real(kind=real64)      :: thetap_subtri(2,8)               ! Matrix that contains the angles thetap of the subtriangles to be integrated
    integer                :: ktheta                           ! Counter variable for theta coordinate loop
    integer                :: krho                             ! Counter variable for rho coordinate loop
    integer                :: ngp_theta                        ! Number of Gauss points for theta coordinate
    integer                :: ngp_rho                          ! Number of Gauss points for rho coordinate
    real(kind=real64)      :: thetai, thetaf, thetapi, thetapf ! Initial and final angles for subtriangle integration
    real(kind=real64)      :: w_angular                        ! Weight of the angular coordinate
    real(kind=real64)      :: w_radial                         ! Weight of the radial coordinate
    real(kind=real64)      :: theta                            ! Angle coordinate theta
    real(kind=real64)      :: thetap                           ! Angle coordinate thetap
    real(kind=real64)      :: thetapp                          ! Angle coordinate thetap on [0,1] domain
    real(kind=real64)      :: jthetap                          ! thetap->thetapp jacobian
    real(kind=real64)      :: costheta, sintheta               ! cos(theta), sin(theta)
    real(kind=real64)      :: rhoij                            ! Maximum rho (radial) value for each edge
    real(kind=real64)      :: rho                              ! Radial coordinate rho
    real(kind=real64)      :: rhop                             ! Radial coordinate rho on [0,1] domain
    real(kind=real64)      :: xi(e%d)                          ! Coordinate xi
    real(kind=real64)      :: aux(10)                          ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)      :: x(3)                             ! Position vector at xi
    real(kind=real64)      :: N(3), T1(3), T2(3)               ! Normal and tangent vectors at xi
    real(kind=real64)      :: rv(3)                            ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)      :: r, d1r                           ! Distance vector module, its inverse and log(1/r)
    real(kind=real64)      :: drdx(3)                          ! Distance vector derivatives with respect to x_k
    real(kind=real64)      :: jg                               ! Geometric jacobian
    real(kind=real64)      :: jw                               ! Jacobians * weight
    real(kind=real64)      :: sphijw(e%n_snodes)               ! phi^s * jw
    real(kind=real64)      :: cteu1, cteu2                     ! Auxiliary constants
    real(kind=real64)      :: fs_u                             ! Fundamental solutions values and auxiliary parts
    real(kind=real64)      :: gbl_tmp(e%n_snodes,3,3)          ! gbl kernels
    ! Initialize kernels
    gbl=0.d0
    ! Calculate
    select case (e%d)
      !
      ! SURFACE LOAD
      !
      case (2)
        ! Initialize auxiliary constants for fundamental solutions calculation
        cteu1=1.d0/(16.d0*c_pi*mu*(1.d0-nu))
        cteu2=3.d0-4.d0*nu
        ! Calculate spatial coordinates at xi_i
#       define etype e%gtype
#       define delta 0.d0
#       define xi xi_i
#       define phi gphi
#       include <phi_2d.rc>
#       undef etype
#       undef delta
#       undef xi
#       undef phi
        ! Calculate x_i
        x_i=0.d0
        do kphi=1,e%n_gnodes
          x_i=x_i+gphi(kphi)*e%x(:,kphi)
        end do
        ! Setup of the polar transformation
        call fbem_polar_transformation_setup(e%gtype,xi_i,nsubtri,subtriangle,theta_subtri,thetap_subtri)
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
            ! THETAPP COORDINATE
            thetapp=gl01_xi(ktheta,ngp_theta)
            w_angular=gl01_w(ktheta,ngp_theta)
            ! THETAPP -> THETAP TRANSFORMATION
            jthetap=(thetapf-thetapi)
            thetap=jthetap*thetapp+thetapi
            ! THETAP -> THETA TRANSFORMATION
            call fbem_polar_transformation_angular(e%gtype,xi_i,subtriangle(ksubtri),thetap,theta,rhoij)
            ! Save cos(theta) and sin(theta)
            costheta=cos(theta)
            sintheta=sin(theta)
            ! Loop through radial coordinate
            do krho=1,gl01_n(ngp_rho)
              ! RHOP -> RHO TRANSFORMATION
              rhop=gl01_xi(krho,ngp_rho)
              w_radial=gl01_w(krho,ngp_rho)
              rho=rhoij*rhop
              ! RHO,THETA -> XI TRANSFORMATION
              xi(1)=xi_i(1)+rho*costheta
              xi(2)=xi_i(2)+rho*sintheta
              ! XI->X TRANSFORMATION
#             define etype e%gtype
#             define delta 0.d0
#             define phi gphi
#             define dphidxi1 dgphidxi1
#             define dphidxi2 dgphidxi2
#             include <phi_and_dphidxik_2d.rc>
#             undef etype
#             undef delta
#             undef phi
#             undef dphidxi1
#             undef dphidxi2
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
              ! xi->x jacobian
              jg=sqrt(dot_product(N,N))
              ! Distance vector and other derived terms
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              d1r=1.d0/r
              drdx=rv*d1r
              ! Jacobians * weights
              jw=jg*rho*jthetap*w_angular*w_radial
              ! FUNCTIONAL SHAPE FUNCTIONS
              ! Functional shape functions (body load variables) at xi
#             define etype e%stype
#             define delta e%stype_delta
#             define phi sphi
#             include <phi_2d.rc>
#             undef etype
#             undef delta
#             undef phi
              ! Functional shape functions * jacobians * weights
              sphijw=sphi*jw
              do il=1,3
                do ik=1,3
                  fs_u=d1r*(cteu2*c_dkr(il,ik)+drdx(il)*drdx(ik))
                  gbl(:,il,ik)=gbl(:,il,ik)+fs_u*sphijw
                end do
              end do
            end do
          end do
        end do
        ! Multiply by constants
        gbl=cteu1*gbl
      !
      ! VOLUME LOAD
      !
      case (3)
        stop 'not yet : fbem_bem_staela3d_sbie_bl_int'
      !
      ! OTHERS
      !
      case default
        call fbem_error_message(output_unit,0,'fbem_bem_staela3d_sbie_bl_int',0,'it is only possible to integrate surface or volume loads')
    end select
  end subroutine fbem_bem_staela3d_sbie_bl_int

  !! Efficient automatic integration of BE body load elements
  subroutine fbem_bem_staela3d_sbie_bl_auto_cmplx(e,x_i,mu,nu,qsp,ns,gbl)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                   !! Integration element
    real(kind=real64)        :: x_i(3)              !! Collocation point
    complex(kind=real64)     :: mu                  !! Shear modulus
    complex(kind=real64)     :: nu                  !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp                 !! Quasi-singular integration parameters
    integer                  :: ns                  !! Maximum level of subdivisions
    complex(kind=real64)     :: gbl(e%n_snodes,3,3) !! gbl integration kernel
    ! Local
    real(kind=real64)        :: r(3)                               ! Distance vector
    real(kind=real64)        :: rmin                               ! Minimum distance between element and x_i
    real(kind=real64)        :: barxi(e%d)                         ! Local coordinates of the nearest element point with respect to x_i
    real(kind=real64)        :: d                                  ! Dimensionless distance
    integer                  :: delta                              ! Control variable
    real(kind=real64)        :: xi_s(e%d,fbem_n_vertices(e%gtype)) ! Local coordinates of the element subdivision
    integer                  :: method                             ! Method used when calculating the nearest element point
    integer                  :: gln_near                           ! 1D Gauss-Legendre integration points required by the quasi-singular function
    integer                  :: gln                                ! 1D Gauss-Legendre integration points used in the integration
    integer                  :: ps                                 ! Selected precalculated dataset
    integer                  :: i                                  ! Counter
    ! POINT BODY LOAD
    if (e%d.eq.0) then
      r=e%x(:,1)-x_i
      rmin=sqrt(dot_product(r,r))
      if (rmin.eq.0.d0) then
        call fbem_error_message(output_unit,0,'fbem_bem_staela3d_sbie_bl_auto',0,'it is not possible to collocate at a point load')
      else
        call fbem_bem_staela3d_sbie_u_cmplx(e%x(:,1),x_i,mu,nu,gbl(1,:,:))
        return
      end if
    ! LINE, SURFACE OR VOLUME BODY LOAD
    ! Determine if interior or exterior integration
    !   - Interior integration (delta=1) requires: xi_i
    !   - Exterior integration (delta=0) requires: x_i, barxi, rmin and d
    ! Use the element ball
    else
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
    end if
    ! Integrate
    select case (delta)
      case (1)
        call fbem_bem_staela3d_sbie_bl_int_cmplx(e,barxi,mu,nu,qsp,ns,gbl)
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
          call fbem_bem_staela3d_sbie_bl_ext_pre_cmplx(ps,e,x_i,mu,nu,gbl)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_staela3d_sbie_bl_ext_adp_cmplx(e,xi_s,x_i,mu,nu,qsp,1,ns,gbl)
        end if
    end select
  end subroutine fbem_bem_staela3d_sbie_bl_auto_cmplx

  !! This subroutine calculates the kernels for SBIE exterior integration needing precalculated data at integration points.
  subroutine fbem_bem_staela3d_sbie_bl_ext_pre_cmplx(ps,e,x_i,mu,nu,gbl)
    implicit none
    ! I/O
    integer                :: ps                  !! Selected precalculated dataset
    type(fbem_bem_element) :: e                   !! Element
    real(kind=real64)      :: x_i(3)              !! Collocation point position vector
    complex(kind=real64)      :: mu                  !! Shear modulus
    complex(kind=real64)      :: nu                  !! Poisson's ratio
    complex(kind=real64)      :: gbl(e%n_snodes,3,3) !! gbl kernels
    ! Local
    integer                :: il, ik              ! Counter for load / observation components
    integer                :: kip                 ! Counter variable for integration points loop
    real(kind=real64)      :: x(3)                ! Position vector at integration point
    real(kind=real64)      :: sphijw(e%n_snodes)  ! phi^s * jacobian * weight at integration point
    real(kind=real64)      :: rv(3)               ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)      :: r, d1r              ! Distance vector module and its inverse
    real(kind=real64)      :: drdx(3)             ! Distance vector derivatives with respect to x_k
    complex(kind=real64)      :: cteu1, cteu2        ! Auxiliary constants
    complex(kind=real64)      :: fs_u                ! Fundamental solutions values
    ! Initialize kernels
    gbl=0.d0
    ! Initialize auxiliary constants for fundamental solutions calculation
    cteu1=1.d0/(16.d0*c_pi*mu*(1.d0-nu))
    cteu2=3.d0-4.d0*nu
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      x=e%ps_x(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      d1r=1.d0/r
      drdx=rv*d1r
      do il=1,3
        do ik=1,3
          fs_u=d1r*(cteu2*c_dkr(il,ik)+drdx(il)*drdx(ik))
          gbl(:,il,ik)=gbl(:,il,ik)+fs_u*sphijw
        end do
      end do
    end do ! Loop through integrations points
    ! Multiply by constants
    gbl=cteu1*gbl
  end subroutine fbem_bem_staela3d_sbie_bl_ext_pre_cmplx

  !! This subroutine calculates adaptatively the body load kernels for SBIE exterior integration using Telles transformation and
  !! subdivision if needed.
  recursive subroutine fbem_bem_staela3d_sbie_bl_ext_adp_cmplx(e,xi_s,x_i,mu,nu,qsp,ks,ns,gbl)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                      !! Element
    real(kind=real64)        :: xi_s(e%d,fbem_n_vertices(e%gtype))     !! Subdivision of the parent element
    real(kind=real64)        :: x_i(3)                                 !! Collocation point position vector
    complex(kind=real64)        :: mu                                     !! Shear modulus
    complex(kind=real64)        :: nu                                     !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp                                    !! Quasi-singular integration parameters
    integer                  :: ks                                     !! Current level of subdivisions
    integer                  :: ns                                     !! Maximum level of subdivision
    complex(kind=real64)        :: gbl(e%n_snodes,3,3)                    !! gbl kernels
    ! Local
    integer                  :: gln_near                               ! 1D Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer                  :: gln                                    ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical                  :: subdivide                              ! True if subdivision has to be performed
    real(kind=real64)        :: barxi(e%d)                             ! Nearest element coordinate with respect to collocation point
    real(kind=real64)        :: barxip(e%d)                            ! Nearest element subdivision local coordinate with respect to collocation point
    real(kind=real64)        :: rmin                                   ! Minimum distance between collocation point and barxi on the element
    real(kind=real64)        :: barr                                   ! Telles jacobian at barxi
    real(kind=real64)        :: cl                                     ! Characteristic length
    real(kind=real64)        :: d                                      ! Normalized distance between collocation point and element subdivision
    integer                  :: method                                 ! Method used in nearest point algorithm
    real(kind=real64)        :: tmp_xi_s(e%d,fbem_n_vertices(e%gtype)) ! Subdivision
    real(kind=real64)        :: x_s(3,e%n_gnodes)                      ! Coordinates of the element subdivision
    complex(kind=real64)        :: gbl_tmp(e%n_snodes,3,3)                ! gbl kernels (temporary)
    ! Initialize
    if (ks.eq.1) then
      gbl=0.d0
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
        call fbem_warning_message(error_unit,0,'fbem_bem_staela3d_sbie_bl_ext_adp_cmplx',ns,'maximum number of subdivisions reached')
        gln_near=30
      end if
    else
      if (gln_near.eq.0) subdivide=.true.
    end if
    ! Subdivide
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
          call fbem_bem_staela3d_sbie_bl_ext_adp_cmplx(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
          ! SUBLINE 2
          tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
          tmp_xi_s(1,2)=xi_s(1,2)
          call fbem_bem_staela3d_sbie_bl_ext_adp_cmplx(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
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
              call fbem_bem_staela3d_sbie_bl_ext_adp_cmplx(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
              ! SUBTRI 2
              tmp_xi_s(:,1)=xi_s(:,2)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              call fbem_bem_staela3d_sbie_bl_ext_adp_cmplx(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
              ! SUBTRI 3
              tmp_xi_s(:,1)=xi_s(:,3)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              call fbem_bem_staela3d_sbie_bl_ext_adp_cmplx(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
              ! SUBTRI 4
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
              call fbem_bem_staela3d_sbie_bl_ext_adp_cmplx(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
            ! QUADRILATERALS
            case (4)
              ! SUBQUAD 1
              tmp_xi_s(:,1)=xi_s(:,1)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,3)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,4)=0.50d0*(xi_s(:,1)+xi_s(:,4))
              call fbem_bem_staela3d_sbie_bl_ext_adp_cmplx(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
              ! SUBQUAD 2
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,2)=xi_s(:,2)
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,4)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              call fbem_bem_staela3d_sbie_bl_ext_adp_cmplx(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
              ! SUBQUAD 3
              tmp_xi_s(:,1)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=xi_s(:,3)
              tmp_xi_s(:,4)=0.50d0*(xi_s(:,3)+xi_s(:,4))
              call fbem_bem_staela3d_sbie_bl_ext_adp_cmplx(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
              ! SUBQUAD 4
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,4))
              tmp_xi_s(:,2)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,4)=xi_s(:,4)
              call fbem_bem_staela3d_sbie_bl_ext_adp_cmplx(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
          end select
        !
        ! VOLUME LOAD
        !
        case (3)
          stop 'not yet : fbem_bem_staela3d_sbie_bl_ext_adp_cmplx'
      end select
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_staela3d_sbie_bl_ext_st_cmplx(e,xi_s,x_i,barxip,barr,mu,nu,gln,gbl_tmp)
      gbl=gbl+gbl_tmp
    end if
  end subroutine fbem_bem_staela3d_sbie_bl_ext_adp_cmplx

  !! This subroutine calculates the body load kernels for SBIE exterior integration (near collocation points) using Telles
  !! transformation within a subdivision of the element, needing only basic data.
  subroutine fbem_bem_staela3d_sbie_bl_ext_st_cmplx(e,xi_s,x_i,barxip,barr,mu,nu,gln,gbl)
    implicit none
    ! I/O
    type(fbem_bem_element)       :: e                                  !! Integration element
    real(kind=real64)            :: xi_s(e%d,fbem_n_vertices(e%gtype)) !! Coordinates of the subdivision of the parent element (xi space [-1,1])
    real(kind=real64)            :: x_i(3)                             !! Collocation point position vector
    real(kind=real64)            :: barxip(e%d)                        !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)            :: barr                               !! Telles jacobian at barxip
    complex(kind=real64)            :: mu                                 !! Shear modulus
    complex(kind=real64)            :: nu                                 !! Poisson's ratio
    integer                      :: gln                                !! 1D Gauss-Legendre number of integration points (<=32)
    complex(kind=real64)            :: gbl(e%n_snodes,3,3)                !! gbl kernels
    ! Local
    integer                      :: il, ik                             ! Counter variables for load direction and observation direction
    integer                      :: kphi                               ! Counter variable for shape functions loops
    integer                      :: kip                                ! Counter variable of integration points
    integer                      :: k1                                 ! Counter variable for reference coordinate xi_1
    integer                      :: k2                                 ! Counter variable for reference coordinate xi_2
    real(kind=real64)            :: gamma(e%d)                         ! Coordinate gamma (Telles transformation space [-1,1])
    real(kind=real64)            :: w(e%d)                             ! Weights of an integration point
    type(fbem_telles_parameters) :: telles_parameters(e%d)             ! Telles parameters
    real(kind=real64)            :: jt(e%d)                            ! Telles jacobian: xip->gamma
    real(kind=real64)            :: xip(e%d)                           ! Coordinate xip (subdivision space [-1,1])
    real(kind=real64)            :: js                                 ! Subdivision jacobian: xi->xip
    real(kind=real64)            :: xin(e%d)                           ! Coordinate xi [xi_s(1,1),xi_s(1,2)]
    real(kind=real64)            :: aux(10)                            ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)                   ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi(e%n_gnodes)               ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: dgphidxi1(e%n_gnodes)              ! Geometrical shape functions derivatives values
    real(kind=real64)            :: dgphidxi2(e%n_gnodes)              ! Geometrical shape functions derivatives values
    real(kind=real64)            :: dxidxi1p(2)                        ! xi derivatives with respect to xi1p
    real(kind=real64)            :: dxidxi2p(2)                        ! xi derivatives with respect to xi2p
    real(kind=real64)            :: xipp(2)                            ! Coordinate xipp used for quadrilateral-triangle transformation
    real(kind=real64)            :: barxipp(2)                         ! Coordinate xipp of collocation point
    real(kind=real64)            :: jqt                                ! Jacobian of the quadrilateral-triangle transformation
    real(kind=real64)            :: sphi(e%n_snodes)                   ! Functional shape functions values
    real(kind=real64)            :: x(3)                               ! Position vector at xi
    real(kind=real64)            :: N(3), T(3), T1(3), T2(3)           ! Normal and tangent vectors at xi
    real(kind=real64)            :: rv(3)                              ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)            :: r, d1r                             ! Distance vector module, its inverse and log(1/r)
    real(kind=real64)            :: drdx(3)                            ! Distance vector derivatives with respect to x_k
    real(kind=real64)            :: jg                                 ! Geometric jacobian
    real(kind=real64)            :: jw                                 ! Jacobians * weight
    real(kind=real64)            :: sphijw(e%n_snodes)                 ! secondary shape functions * jw
    complex(kind=real64)            :: cteu1, cteu2                       ! Auxiliary constants
    complex(kind=real64)            :: fs_u                               ! Fundamental solutions values
    ! Initialize kernels
    gbl=0.d0
    ! Initialize auxiliary constants for fundamental solutions calculation
    cteu1=1.d0/(16.d0*c_pi*mu*(1.d0-nu))
    cteu2=3.d0-4.d0*nu
    ! Calculate
    select case (e%d)
      !
      ! LINE LOAD
      !
      case (1)
        ! Calculate Telles parameters
        telles_parameters(1)=fbem_telles11_calculate_parameters(barxip(1),barr)
        ! Numerical integration
        do kip=1,gl11_n(gln)
          ! GAMMA COORDINATE
          gamma(1)=gl11_xi(kip,gln)
          w(1)=gl11_w(kip,gln)
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
          d1r=1.d0/r
          drdx=rv*d1r
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
          ! Add integration points
          do il=1,3
            do ik=1,3
              fs_u=d1r*(cteu2*c_dkr(il,ik)+drdx(il)*drdx(ik))
              gbl(:,il,ik)=gbl(:,il,ik)+fs_u*sphijw
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
#               define etype e%gtype
#               define delta 0.d0
#               define xi xin
#               define phi gphi
#               define dphidxi1 dgphidxi1
#               define dphidxi2 dgphidxi2
#               include <phi_and_dphidxik_2d.rc>
#               undef etype
#               undef delta
#               undef xi
#               undef phi
#               undef dphidxi1
#               undef dphidxi2
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
                ! xi->x jacobian
                jg=sqrt(dot_product(N,N))
                ! Distance vector and derived terms
                rv=x-x_i
                r=sqrt(dot_product(rv,rv))
                d1r=1.d0/r
                drdx=rv*d1r
                ! Jacobians * weights
                jw=jg*js*jt(1)*jt(2)*w(1)*w(2)
                ! FUNCTIONAL SHAPE FUNCTIONS
                ! Functional shape functions (body load variables) at xi
#               define etype e%stype
#               define delta e%stype_delta
#               define xi xin
#               define phi sphi
#               include <phi_2d.rc>
#               undef etype
#               undef delta
#               undef xi
#               undef phi
                ! Functional shape functions * jacobians * weights
                sphijw=sphi*jw
                ! Add integration points
                do il=1,3
                  do ik=1,3
                    fs_u=d1r*(cteu2*c_dkr(il,ik)+drdx(il)*drdx(ik))
                    gbl(:,il,ik)=gbl(:,il,ik)+fs_u*sphijw
                  end do
                end do
              end do
            end do
          ! TRIANGULAR ELEMENTS
          case (3)
            ! Telles transformation is applied to Gauss-Legendre*Gauss-Legendre quadrature before the quadrilateral->triangle
            ! transformation. Because barxip for triangles are given in area triangle coordinates, for Telles transformation
            ! they must be transformed to quadrilateral coordinates. A special treatment is needed when barxi_2 is near 1, because
            ! transformation diverges
            ! BARXIP->BARXIPP TRANSFORMATION
            if (barxip(2).gt.0.995d0) then
              barxipp(1)=0.5d0
              barxipp(2)=1.0d0
            else
              barxipp(1)=barxip(1)/(1.d0-barxip(2))
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
#               define delta 0.d0
#               define xi xin
#               define phi gphi
#               define dphidxi1 dgphidxi1
#               define dphidxi2 dgphidxi2
#               include <phi_and_dphidxik_2d.rc>
#               undef etype
#               undef delta
#               undef xi
#               undef phi
#               undef dphidxi1
#               undef dphidxi2
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
                ! xi->x jacobian
                jg=sqrt(dot_product(N,N))
                ! Distance vector and derived terms
                rv=x-x_i
                r=sqrt(dot_product(rv,rv))
                d1r=1.d0/r
                drdx=rv*d1r
                ! Jacobians * weights
                jw=jg*js*jqt*jt(1)*jt(2)*w(1)*w(2)
                ! FUNCTIONAL SHAPE FUNCTIONS
                ! Functional shape functions (body load variables) at xi
#               define etype e%stype
#               define delta e%stype_delta
#               define xi xin
#               define phi sphi
#               include <phi_2d.rc>
#               undef etype
#               undef delta
#               undef xi
#               undef phi
                ! Functional shape functions * jacobians * weights
                sphijw=sphi*jw
                ! Add integration points
                do il=1,3
                  do ik=1,3
                    fs_u=d1r*(cteu2*c_dkr(il,ik)+drdx(il)*drdx(ik))
                    gbl(:,il,ik)=gbl(:,il,ik)+fs_u*sphijw
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
        stop 'not yet : fbem_bem_staela3d_sbie_bl_ext_st_cmplx'
      !
      ! OTHERS
      !
      case default
        call fbem_error_message(output_unit,0,'fbem_bem_staela3d_sbie_bl_int_cmplx',0,&
                                'it is only possible to integrate line, surface or volume loads')
    end select
    ! Multiply by constants
    gbl=cteu1*gbl
  end subroutine fbem_bem_staela3d_sbie_bl_ext_st_cmplx

  !! This subroutine calculates the body load kernels for SBIE interior integration, needing only basic data.
  subroutine fbem_bem_staela3d_sbie_bl_int_cmplx(e,xi_i,mu,nu,qsp,ns,gbl)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                !! Integration element
    real(kind=real64)        :: xi_i(e%d)                        !! Local coordinate of the singular point.
    complex(kind=real64)        :: mu                               !! Shear modulus
    complex(kind=real64)        :: nu                               !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp                              !! Quasi-singular integration parameters
    integer                  :: ns                               !! Maximum level of subdivisions
    complex(kind=real64)        :: gbl(e%n_snodes,3,3)              !! gbl kernels
    ! Local
    integer                :: il, ik                           ! Counter variables for load direction and observation direction
    integer                :: kphi                             ! Counter variable for shape functions loops
    real(kind=real64)      :: x_i(3), x_ic(3)                  ! Real coordinates of collocation point
    real(kind=real64)      :: ep1(3), ep2(3), ep3(3)           ! Local system at the collocation point
    real(kind=real64)      :: gphi(e%n_gnodes)                 ! Geometrical shape functions values at xi
    real(kind=real64)      :: dgphidxi1(e%n_gnodes)            ! Geometrical shape functions derivatives values at xi_1
    real(kind=real64)      :: dgphidxi2(e%n_gnodes)            ! Geometrical shape functions derivatives values at xi_2
    real(kind=real64)      :: sphi(e%n_snodes)                 ! Functional shape functions values at xi
    integer                :: ksubtri                          ! Counter variable for subtriangles loop
    integer                :: nsubtri                          ! Number of subtriangles
    integer                :: subtriangle(8)                   ! Vector that contains what subtriangles need to be integrated
    real(kind=real64)      :: theta_subtri(2,8)                ! Matrix that contains the angles theta of the subtriangles to be integrated
    real(kind=real64)      :: thetap_subtri(2,8)               ! Matrix that contains the angles thetap of the subtriangles to be integrated
    integer                :: ktheta                           ! Counter variable for theta coordinate loop
    integer                :: krho                             ! Counter variable for rho coordinate loop
    integer                :: ngp_theta                        ! Number of Gauss points for theta coordinate
    integer                :: ngp_rho                          ! Number of Gauss points for rho coordinate
    real(kind=real64)      :: thetai, thetaf, thetapi, thetapf ! Initial and final angles for subtriangle integration
    real(kind=real64)      :: w_angular                        ! Weight of the angular coordinate
    real(kind=real64)      :: w_radial                         ! Weight of the radial coordinate
    real(kind=real64)      :: theta                            ! Angle coordinate theta
    real(kind=real64)      :: thetap                           ! Angle coordinate thetap
    real(kind=real64)      :: thetapp                          ! Angle coordinate thetap on [0,1] domain
    real(kind=real64)      :: jthetap                          ! thetap->thetapp jacobian
    real(kind=real64)      :: costheta, sintheta               ! cos(theta), sin(theta)
    real(kind=real64)      :: rhoij                            ! Maximum rho (radial) value for each edge
    real(kind=real64)      :: rho                              ! Radial coordinate rho
    real(kind=real64)      :: rhop                             ! Radial coordinate rho on [0,1] domain
    real(kind=real64)      :: xi(e%d)                          ! Coordinate xi
    real(kind=real64)      :: aux(10)                          ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)      :: x(3)                             ! Position vector at xi
    real(kind=real64)      :: N(3), T1(3), T2(3)               ! Normal and tangent vectors at xi
    real(kind=real64)      :: rv(3)                            ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)      :: r, d1r                           ! Distance vector module, its inverse and log(1/r)
    real(kind=real64)      :: drdx(3)                          ! Distance vector derivatives with respect to x_k
    real(kind=real64)      :: jg                               ! Geometric jacobian
    real(kind=real64)      :: jw                               ! Jacobians * weight
    real(kind=real64)      :: sphijw(e%n_snodes)               ! phi^s * jw
    complex(kind=real64)      :: cteu1, cteu2                     ! Auxiliary constants
    complex(kind=real64)      :: fs_u                             ! Fundamental solutions values and auxiliary parts
    complex(kind=real64)      :: gbl_tmp(e%n_snodes,3,3)          ! gbl kernels
    ! Initialize kernels
    gbl=0.d0
    ! Calculate
    select case (e%d)
      !
      ! SURFACE LOAD
      !
      case (2)
        ! Initialize auxiliary constants for fundamental solutions calculation
        cteu1=1.d0/(16.d0*c_pi*mu*(1.d0-nu))
        cteu2=3.d0-4.d0*nu
        ! Calculate spatial coordinates at xi_i
#       define etype e%gtype
#       define delta 0.d0
#       define xi xi_i
#       define phi gphi
#       include <phi_2d.rc>
#       undef etype
#       undef delta
#       undef xi
#       undef phi
        ! Calculate x_i
        x_i=0.d0
        do kphi=1,e%n_gnodes
          x_i=x_i+gphi(kphi)*e%x(:,kphi)
        end do
        ! Setup of the polar transformation
        call fbem_polar_transformation_setup(e%gtype,xi_i,nsubtri,subtriangle,theta_subtri,thetap_subtri)
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
            ! THETAPP COORDINATE
            thetapp=gl01_xi(ktheta,ngp_theta)
            w_angular=gl01_w(ktheta,ngp_theta)
            ! THETAPP -> THETAP TRANSFORMATION
            jthetap=(thetapf-thetapi)
            thetap=jthetap*thetapp+thetapi
            ! THETAP -> THETA TRANSFORMATION
            call fbem_polar_transformation_angular(e%gtype,xi_i,subtriangle(ksubtri),thetap,theta,rhoij)
            ! Save cos(theta) and sin(theta)
            costheta=cos(theta)
            sintheta=sin(theta)
            ! Loop through radial coordinate
            do krho=1,gl01_n(ngp_rho)
              ! RHOP -> RHO TRANSFORMATION
              rhop=gl01_xi(krho,ngp_rho)
              w_radial=gl01_w(krho,ngp_rho)
              rho=rhoij*rhop
              ! RHO,THETA -> XI TRANSFORMATION
              xi(1)=xi_i(1)+rho*costheta
              xi(2)=xi_i(2)+rho*sintheta
              ! XI->X TRANSFORMATION
#             define etype e%gtype
#             define delta 0.d0
#             define phi gphi
#             define dphidxi1 dgphidxi1
#             define dphidxi2 dgphidxi2
#             include <phi_and_dphidxik_2d.rc>
#             undef etype
#             undef delta
#             undef phi
#             undef dphidxi1
#             undef dphidxi2
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
              ! xi->x jacobian
              jg=sqrt(dot_product(N,N))
              ! Distance vector and other derived terms
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              d1r=1.d0/r
              drdx=rv*d1r
              ! Jacobians * weights
              jw=jg*rho*jthetap*w_angular*w_radial
              ! FUNCTIONAL SHAPE FUNCTIONS
              ! Functional shape functions (body load variables) at xi
#             define etype e%stype
#             define delta e%stype_delta
#             define phi sphi
#             include <phi_2d.rc>
#             undef etype
#             undef delta
#             undef phi
              ! Functional shape functions * jacobians * weights
              sphijw=sphi*jw
              do il=1,3
                do ik=1,3
                  fs_u=d1r*(cteu2*c_dkr(il,ik)+drdx(il)*drdx(ik))
                  gbl(:,il,ik)=gbl(:,il,ik)+fs_u*sphijw
                end do
              end do
            end do
          end do
        end do
        ! Multiply by constants
        gbl=cteu1*gbl
      !
      ! VOLUME LOAD
      !
      case (3)
        stop 'not yet : fbem_bem_staela3d_sbie_bl_int_cmplx'
      !
      ! OTHERS
      !
      case default
        call fbem_error_message(output_unit,0,'fbem_bem_staela3d_sbie_bl_int_cmplx',0,'it is only possible to integrate surface or volume loads')
    end select
  end subroutine fbem_bem_staela3d_sbie_bl_int_cmplx

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (HBIE)
  !

  !! Fundamental solution d*
  subroutine fbem_bem_staela3d_hbie_d(x,x_i,n_i,nu,do)
    implicit none
    ! I/O
    real(kind=real64) :: x(3)    !! Observation point
    real(kind=real64) :: x_i(3)  !! Collocation point
    real(kind=real64) :: n_i(3)  !! Collocation point unit normal
    real(kind=real64) :: nu      !! Poisson's ratio
    real(kind=real64) :: do(3,3) !! d*_{lk}
    ! Local
    integer           :: il, ik
    real(kind=real64) :: rv(3), r, d1r, d1r2, drdx(3), drdni
    real(kind=real64) :: cted1, cted2
    cted1=-1.d0/(8.d0*c_pi*(1.d0-nu))
    cted2=1.d0-2.d0*nu
    rv=x-x_i
    r=sqrt(dot_product(rv,rv))
    d1r =1.d0/r
    d1r2=d1r**2
    drdx=rv*d1r
    drdni=-dot_product(drdx,n_i)
    do il=1,3
      do ik=1,3
        do(il,ik)=d1r2*((cted2*c_dkr(il,ik)+3.d0*drdx(il)*drdx(ik))*drdni+cted2*(n_i(il)*drdx(ik)-n_i(ik)*drdx(il)))
      end do
    end do
    do=cted1*do
  end subroutine fbem_bem_staela3d_hbie_d

  !! Fundamental solution s*
  subroutine fbem_bem_staela3d_hbie_s(x,n,x_i,n_i,mu,nu,so)
    implicit none
    ! I/O
    real(kind=real64) :: x(3)    !! Observation point
    real(kind=real64) :: n(3)    !! Observation point unit normal
    real(kind=real64) :: x_i(3)  !! Collocation point
    real(kind=real64) :: n_i(3)  !! Collocation point unit normal
    real(kind=real64) :: mu      !! Shear modulus
    real(kind=real64) :: nu      !! Poisson's ratio
    real(kind=real64) :: so(3,3) !! s*_{lk}
    ! Local
    integer           :: il, ik
    real(kind=real64) :: rv(3), r, d1r, d1r2, d1r3, drdx(3), drdn, drdni, n_dot_ni
    real(kind=real64) :: ctes1, ctes2, ctes3
    ctes1=mu/(4.d0*c_pi*(1.d0-nu))
    ctes2=1.d0-2.d0*nu
    ctes3=1.d0-4.d0*nu
    rv=x-x_i
    r=sqrt(dot_product(rv,rv))
    d1r =1.d0/r
    d1r2=d1r**2
    d1r3=d1r2*d1r
    drdx=rv*d1r
    drdn=dot_product(drdx,n)
    drdni=-dot_product(drdx,n_i)
    n_dot_ni=dot_product(n,n_i)
    do il=1,3
      do ik=1,3
        so(il,ik)=d1r3*(3.d0*(5.d0*drdx(il)*drdx(ik)-nu*c_dkr(il,ik))*drdn*drdni &
                       +3.d0*ctes2*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni) &
                       +3.d0*nu*   (drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni) &
                       +(3.d0*nu*drdx(il)*drdx(ik)+ctes2*c_dkr(il,ik))*n_dot_ni &
                       +ctes2*n(il)*n_i(ik)-ctes3*n(ik)*n_i(il))
      end do
    end do
    so=ctes1*so
  end subroutine fbem_bem_staela3d_hbie_s

  !! Fundamental solution d*
  subroutine fbem_bem_staela3d_hbie_d_cmplx(x,x_i,n_i,nu,do)
    implicit none
    ! I/O
    real(kind=real64)    :: x(3)    !! Observation point
    real(kind=real64)    :: x_i(3)  !! Collocation point
    real(kind=real64)    :: n_i(3)  !! Collocation point unit normal
    complex(kind=real64) :: nu      !! Poisson's ratio
    complex(kind=real64) :: do(3,3) !! d*_{lk}
    ! Local
    integer              :: il, ik
    real(kind=real64)    :: rv(3), r, d1r, d1r2, drdx(3), drdni
    complex(kind=real64) :: cted1, cted2
    cted1=-1.d0/(8.d0*c_pi*(1.d0-nu))
    cted2=1.d0-2.d0*nu
    rv=x-x_i
    r=sqrt(dot_product(rv,rv))
    d1r =1.d0/r
    d1r2=d1r**2
    drdx=rv*d1r
    drdni=-dot_product(drdx,n_i)
    do il=1,3
      do ik=1,3
        do(il,ik)=d1r2*((cted2*c_dkr(il,ik)+3.d0*drdx(il)*drdx(ik))*drdni+cted2*(n_i(il)*drdx(ik)-n_i(ik)*drdx(il)))
      end do
    end do
    do=cted1*do
  end subroutine fbem_bem_staela3d_hbie_d_cmplx

  !! Fundamental solution s*
  subroutine fbem_bem_staela3d_hbie_s_cmplx(x,n,x_i,n_i,mu,nu,so)
    implicit none
    ! I/O
    real(kind=real64)    :: x(3)    !! Observation point
    real(kind=real64)    :: n(3)    !! Observation point unit normal
    real(kind=real64)    :: x_i(3)  !! Collocation point
    real(kind=real64)    :: n_i(3)  !! Collocation point unit normal
    complex(kind=real64) :: mu      !! Shear modulus
    complex(kind=real64) :: nu      !! Poisson's ratio
    complex(kind=real64) :: so(3,3) !! s*_{lk}
    ! Local
    integer           :: il, ik
    real(kind=real64) :: rv(3), r, d1r, d1r2, d1r3, drdx(3), drdn, drdni, n_dot_ni
    real(kind=real64) :: ctes1, ctes2, ctes3
    ctes1=mu/(4.d0*c_pi*(1.d0-nu))
    ctes2=1.d0-2.d0*nu
    ctes3=1.d0-4.d0*nu
    rv=x-x_i
    r=sqrt(dot_product(rv,rv))
    d1r =1.d0/r
    d1r2=d1r**2
    d1r3=d1r2*d1r
    drdx=rv*d1r
    drdn=dot_product(drdx,n)
    drdni=-dot_product(drdx,n_i)
    n_dot_ni=dot_product(n,n_i)
    do il=1,3
      do ik=1,3
        so(il,ik)=d1r3*(3.d0*(5.d0*drdx(il)*drdx(ik)-nu*c_dkr(il,ik))*drdn*drdni &
                       +3.d0*ctes2*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni) &
                       +3.d0*nu*   (drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni) &
                       +(3.d0*nu*drdx(il)*drdx(ik)+ctes2*c_dkr(il,ik))*n_dot_ni &
                       +ctes2*n(il)*n_i(ik)-ctes3*n(ik)*n_i(il))
      end do
    end do
    so=ctes1*so
  end subroutine fbem_bem_staela3d_hbie_s_cmplx

  !! Efficient automatic integration
  subroutine fbem_bem_staela3d_hbie_auto(e,reverse,x_i,n_i,mu,nu,qsp,ns,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                             !! Integration element
    logical                  :: reverse                       !! Reverse orientation
    real(kind=real64)        :: x_i(3)                        !! Collocation point
    real(kind=real64)        :: n_i(3)                        !! Unit normal at the collocation point
    real(kind=real64)        :: mu                            !! Shear modulus
    real(kind=real64)        :: nu                            !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp                           !! Quasi-singular integration parameters
    integer                  :: ns                            !! Maximum level of subdivisions
    real(kind=real64)        :: m(e%n_pnodes,3,3)             !! m integration kernel
    real(kind=real64)        :: l(e%n_snodes,3,3)             !! l integration kernel
    ! Local
    real(kind=real64)        :: r(3)                             ! Distance vector
    real(kind=real64)        :: rmin                             ! Minimum distance between element and x_i
    real(kind=real64)        :: barxi(2)                         ! Local coordinates of the nearest element point with respect to x_i
    real(kind=real64)        :: d                                ! Dimensionless distance
    integer                  :: delta                            ! Control variable
    real(kind=real64)        :: xi_s(2,fbem_n_vertices(e%gtype)) ! Local coordinates of the element subdivision
    integer                  :: method                           ! Method used when calculating the nearest element point
    integer                  :: gln_near                         ! 1D Gauss-Legendre integration points required by the quasi-singular function
    integer                  :: gln                              ! 1D Gauss-Legendre integration points used in the integration
    integer                  :: ps                               ! Selected precalculated dataset
    integer                  :: i                                ! Counter
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
        call fbem_bem_staela3d_hbie_int(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi,mu,nu,m,l)
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
          call fbem_bem_staela3d_hbie_ext_pre(ps,e,reverse,x_i,n_i,mu,nu,m,l)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_staela3d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,mu,nu,qsp,1,ns,m,l)
        end if
    end select
  end subroutine fbem_bem_staela3d_hbie_auto

  !! This subroutine calculates the kernels for HBIE exterior integration, needing precalculated data at integration points.
  subroutine fbem_bem_staela3d_hbie_ext_pre(ps,e,reverse,x_i,n_i,mu,nu,m,l)
    implicit none
    ! I/O
    integer                :: ps                !! Selected precalculated dataset
    type(fbem_bem_element) :: e                 !! Element
    logical                :: reverse           !! Reverse normal vector
    real(kind=real64)      :: x_i(3)            !! Collocation point position vector
    real(kind=real64)      :: n_i(3)            !! Unit normal at the collocation point
    real(kind=real64)      :: mu                !! Shear modulus
    real(kind=real64)      :: nu                !! Poisson's ratio
    real(kind=real64)      :: m(e%n_pnodes,3,3) !! m integration kernels matrix
    real(kind=real64)      :: l(e%n_snodes,3,3) !! l integration kernels matrix
    ! Local
    integer           :: il, ik              ! Counter for load / observation components
    integer           :: kip                 ! Counter variable for integration points loop
    real(kind=real64) :: x(3)                ! Position vector at integration point
    real(kind=real64) :: n(3)                ! Unit normal vector at integration point
    real(kind=real64) :: pphijw(e%n_pnodes)  ! phi^p * jacobian * weight at integration point
    real(kind=real64) :: sphijw(e%n_snodes)  ! phi^s * jacobian * weight at integration point
    real(kind=real64) :: rv(3)               ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64) :: r, d1r, d1r2, d1r3  ! Distance vector module and its inverse
    real(kind=real64) :: drdx(3)             ! Distance vector derivatives with respect to x_k
    real(kind=real64) :: drdn                ! Partial derivative of r respect to unit normal
    real(kind=real64) :: drdni               ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64) :: n_dot_ni            ! Dot product of n and n_i
    real(kind=real64) :: cted1, cted2        ! Auxiliary constants
    real(kind=real64) :: ctes1, ctes2, ctes3 ! Auxiliary constants
    real(kind=real64) :: fs_d, fs_s          ! Fundamental solutions values
    ! Initialize kernels
    m=0.d0
    l=0.d0
    ! Initialize auxiliary constants for fundamental solutions calculation
    cted1=-1.d0/(8.d0*c_pi*(1.d0-nu))
    cted2=1.d0-2.d0*nu
    ctes1=mu/(4.d0*c_pi*(1.d0-nu))
    ctes2=1.d0-2.d0*nu
    ctes3=1.d0-4.d0*nu
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      x=e%ps_x(:,kip,ps)
      n=e%ps_n(:,kip,ps)
      pphijw=e%ps_pphijw(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      d1r =1.d0/r
      d1r2=d1r**2
      d1r3=d1r2*d1r
      drdx=rv*d1r
      drdn=dot_product(drdx,n)
      drdni=-dot_product(drdx,n_i)
      n_dot_ni=dot_product(n,n_i)
      do il=1,3
        do ik=1,3
          fs_d=d1r2*((cted2*c_dkr(il,ik)+3.d0*drdx(il)*drdx(ik))*drdni+cted2*(n_i(il)*drdx(ik)-n_i(ik)*drdx(il)))
          fs_s=d1r3*(3.d0*(5.d0*drdx(il)*drdx(ik)-nu*c_dkr(il,ik))*drdn*drdni &
                    +3.d0*ctes2*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni) &
                    +3.d0*nu*   (drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni) &
                    +(3.d0*nu*drdx(il)*drdx(ik)+ctes2*c_dkr(il,ik))*n_dot_ni &
                    +ctes2*n(il)*n_i(ik)-ctes3*n(ik)*n_i(il))
          m(:,il,ik)=m(:,il,ik)+fs_s*pphijw
          l(:,il,ik)=l(:,il,ik)+fs_d*sphijw
        end do
      end do
    end do
    ! Multiply by constants
    m=ctes1*m
    l=cted1*l
    ! Reverse if needed
    if (reverse) m=-m
  end subroutine fbem_bem_staela3d_hbie_ext_pre

  !! This subroutine calculates the kernels for HBIE exterior integration (near collocation points) using Telles transformation
  !! within a subdivision of the element, needing only needs basic data.
  subroutine fbem_bem_staela3d_hbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,mu,nu,gln,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)       :: e                                !! Integration element
    logical                      :: reverse                          !! Reverse normal vector
    real(kind=real64)            :: xi_s(2,fbem_n_vertices(e%gtype)) !! Coordinates of the subdivision of the parent element
    real(kind=real64)            :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)            :: n_i(3)                           !! Unit normal at the collocation point
    real(kind=real64)            :: barxip(2)                        !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)            :: barr                             !! Telles jacobian at barxip
    real(kind=real64)            :: mu                               !! Shear modulus
    real(kind=real64)            :: nu                               !! Poisson's ratio
    integer                      :: gln                              !! 1D Gauss-Legendre number of integration points (<=32)
    real(kind=real64)            :: m(e%n_pnodes,3,3)                !! m kernel vector
    real(kind=real64)            :: l(e%n_snodes,3,3)                !! l kernel vector
    ! Local
    integer                      :: kphi                  ! Counter variable for shape functions loops
    integer                      :: k1                    ! Counter variable for reference coordinate xi_1
    integer                      :: k2                    ! Counter variable for reference coordinate xi_2
    real(kind=real64)            :: aux(10)               ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)      ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi1(e%n_gnodes) ! Geometrical shape functions derivatives values
    real(kind=real64)            :: dgphidxi2(e%n_gnodes) ! Geometrical shape functions derivatives values
    real(kind=real64)            :: pphi(e%n_pnodes)      ! Functional shape functions values
    real(kind=real64)            :: sphi(e%n_snodes)      ! Functional shape functions values
    real(kind=real64)            :: gamma(2)              ! Vector of gamma_1,gamma_2 coordinates
    real(kind=real64)            :: w(2)                  ! Weights of the integration rule
    real(kind=real64)            :: xip(2)                ! Vector of xip_1,xip_2 coordinates
    real(kind=real64)            :: dxidxi1p(2)           ! xi derivatives with respect to xi1p
    real(kind=real64)            :: dxidxi2p(2)           ! xi derivatives with respect to xi2p
    real(kind=real64)            :: js                    ! Subdivision jacobian
    real(kind=real64)            :: xi(2)                 ! Vector of xi_1,xi_2 coordinates
    real(kind=real64)            :: xipp(2)               ! Coordinate xipp used for quadrilateral-triangle transformation
    real(kind=real64)            :: barxipp(2)            ! Coordinate xipp of collocation point
    real(kind=real64)            :: jqt                   ! Jacobian of the quadrilateral-triangle transformation
    type(fbem_telles_parameters) :: telles_parameters(2)  ! Telles parameters for each coordinate
    real(kind=real64)            :: jt(2)                 ! Telles jacobian for each coordinate: xi_1->gamma_1 and xi_2->gamma_2
    real(kind=real64)            :: x(3)                  ! Position vector at xi_1,xi_2
    real(kind=real64)            :: T1(3), T2(3)          ! Tangent vectors at xi_1,xi_2
    real(kind=real64)            :: N(3)                  ! Normal vector at xi_1,xi_2
    real(kind=real64)            :: rv(3)                 ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)            :: r, d1r, d1r2, d1r3    ! Distance vector module and its inverse
    real(kind=real64)            :: drdx(3)               ! Distance vector derivatives with respect to x_k
    real(kind=real64)            :: jg                    ! Geometric jacobian
    real(kind=real64)            :: jw                    ! Jacobians * weights
    real(kind=real64)            :: drdn                  ! Partial derivative of r respect to unit normal
    real(kind=real64)            :: drdni                 ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)            :: n_dot_ni              ! Dot product of n and n_i
    real(kind=real64)            :: pphijw(e%n_pnodes)    ! Functional shape functions * jw
    real(kind=real64)            :: sphijw(e%n_snodes)    ! Functional shape functions * jw
    integer                      :: il, ik                ! Counter for load / observation components
    real(kind=real64)            :: cted1, cted2          ! Auxiliary constants
    real(kind=real64)            :: ctes1, ctes2, ctes3   ! Auxiliary constants
    real(kind=real64)            :: fs_d, fs_s            ! Fundamental solutions values
    ! Initialize kernels
    m=0.d0
    l=0.d0
    ! Auxiliary constants for fundamental solutions calculation
    cted1=-1.d0/(8.d0*c_pi*(1.d0-nu))
    cted2=1.d0-2.d0*nu
    ctes1=mu/(4.d0*c_pi*(1.d0-nu))
    ctes2=1.d0-2.d0*nu
    ctes3=1.d0-4.d0*nu
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
            r=dot_product(rv,rv)
            r=sqrt(r)
            d1r =1.d0/r
            d1r2=d1r**2
            d1r3=d1r2*d1r
            drdx=rv*d1r
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
            ! Add integration points
            do il=1,3
              do ik=1,3
                fs_d=d1r2*((cted2*c_dkr(il,ik)+3.d0*drdx(il)*drdx(ik))*drdni+cted2*(n_i(il)*drdx(ik)-n_i(ik)*drdx(il)))
                fs_s=d1r3*(3.d0*(5.d0*drdx(il)*drdx(ik)-nu*c_dkr(il,ik))*drdn*drdni &
                          +3.d0*ctes2*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni) &
                          +3.d0*nu*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni) &
                          +(3.d0*nu*drdx(il)*drdx(ik)+ctes2*c_dkr(il,ik))*n_dot_ni &
                          +ctes2*n(il)*n_i(ik)-ctes3*n(ik)*n_i(il))
                m(:,il,ik)=m(:,il,ik)+fs_s*pphijw
                l(:,il,ik)=l(:,il,ik)+fs_d*sphijw
              end do
            end do
          end do
        end do
      ! TRIANGULAR ELEMENTS
      case (3)
        ! Telles transformation is applied to Gauss-Legendre*Gauss-Legendre quadrature before the quadrilateral->triangle
        ! transformation. Because barxip for triangles are given in area triangle coordinates, for Telles transformation
        ! they must be transformed to quadrilateral coordinates. A special treatment is needed when barxi_2 is near 1, because
        ! transformation diverges
        ! BARXIP->BARXIPP TRANSFORMATION
        if (barxip(2).gt.0.995d0) then
          barxipp(1)=0.5d0
          barxipp(2)=1.0d0
        else
          barxipp(1)=barxip(1)/(1.d0-barxip(2))
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
            r=dot_product(rv,rv)
            r=sqrt(r)
            d1r =1.d0/r
            d1r2=d1r**2
            d1r3=d1r2*d1r
            drdx=rv*d1r
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
            ! Add integration points
            do il=1,3
              do ik=1,3
                fs_d=d1r2*((cted2*c_dkr(il,ik)+3.d0*drdx(il)*drdx(ik))*drdni+cted2*(n_i(il)*drdx(ik)-n_i(ik)*drdx(il)))
                fs_s=d1r3*(3.d0*(5.d0*drdx(il)*drdx(ik)-nu*c_dkr(il,ik))*drdn*drdni &
                          +3.d0*ctes2*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni) &
                          +3.d0*nu*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni) &
                          +(3.d0*nu*drdx(il)*drdx(ik)+ctes2*c_dkr(il,ik))*n_dot_ni &
                          +ctes2*n(il)*n_i(ik)-ctes3*n(ik)*n_i(il))
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
    m=ctes1*m
    l=cted1*l
    ! Reverse if needed
    if (reverse) m=-m
  end subroutine fbem_bem_staela3d_hbie_ext_st

  !! This subroutine calculates adaptatively the kernels for HBIE exterior integration using Telles transformation and subdivision
  !! if needed, using only basic data.
  recursive subroutine fbem_bem_staela3d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,mu,nu,qsp,ks,ns,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                !! Element
    logical                  :: reverse                          !! Reverse orientation
    real(kind=real64)        :: xi_s(2,fbem_n_vertices(e%gtype)) !! Subdivision of the parent element
    real(kind=real64)        :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)        :: n_i(3)                           !! Unit normal at the collocation point
    real(kind=real64)        :: mu                               !! Shear modulus
    real(kind=real64)        :: nu                               !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp                              !! Quasi-singular integration parameters
    integer                  :: ks                               !! Current level of subdivisions
    integer                  :: ns                               !! Maximum level of subdivision
    real(kind=real64)        :: m(e%n_pnodes,3,3)                !! m integration kernels matrix
    real(kind=real64)        :: l(e%n_snodes,3,3)                !! l integration kernels matrix
    ! Local
    integer           :: gln_near                             ! 1D Gauss-Legendre integ. points required to integrate only  the quasi-singular integrand
    integer           :: gln                                  ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
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
    real(kind=real64) :: m_tmp(e%n_pnodes,3,3)                ! m integration kernels matrix (temporary)
    real(kind=real64) :: l_tmp(e%n_snodes,3,3)                ! l integration kernels matrix (temporary)
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
      cl=fbem_characteristic_length(3,e%gtype,x_s,1.d-12)
      call fbem_nearest_element_point_bem(3,e%gtype,x_s,cl,x_i,barxip,rmin,d,method)
    end if
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,7,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_staela3d_hbie_ext_adp',ns,'maximum number of subdivisions reached')
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
          call fbem_bem_staela3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,m,l)
          ! SUBTRI 2
          tmp_xi_s(:,1)=xi_s(:,2)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          call fbem_bem_staela3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,m,l)
          ! SUBTRI 3
          tmp_xi_s(:,1)=xi_s(:,3)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          call fbem_bem_staela3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,m,l)
          ! SUBTRI 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          call fbem_bem_staela3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,m,l)
        ! QUADRILATERALS
        case (4)
          ! SUBQUAD 1
          tmp_xi_s(:,1)=xi_s(:,1)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,3)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          call fbem_bem_staela3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,m,l)
          ! SUBQUAD 2
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=xi_s(:,2)
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,4)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          call fbem_bem_staela3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,m,l)
          ! SUBQUAD 3
          tmp_xi_s(:,1)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=xi_s(:,3)
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          call fbem_bem_staela3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,m,l)
          ! SUBQUAD 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          tmp_xi_s(:,2)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=xi_s(:,4)
          call fbem_bem_staela3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,m,l)
      end select
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_staela3d_hbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,mu,nu,gln,m_tmp,l_tmp)
      m=m+m_tmp
      l=l+l_tmp
    end if
  end subroutine fbem_bem_staela3d_hbie_ext_adp

  !! This subroutine calculates static elastic HBIE integration kernels vectors m and l of singular formulation for interior
  !! integration.
  subroutine fbem_bem_staela3d_hbie_int(type_g,type_f1,type_f2,delta_f,x_nodes,reverse,xi_i,mu,nu,m,l)
    implicit none
    ! I/O
    integer           :: type_g                          !! Geometrial interpolation
    integer           :: type_f1                         !! Functional interpolation (primary variables)
    integer           :: type_f2                         !! Functional interpolation (secondary variables)
    real(kind=real64) :: delta_f                         !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64) :: x_nodes(3,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
    logical           :: reverse                         !! Normal vector inversion
    real(kind=real64) :: xi_i(2)                         !! Reference coordinates of the singular point.
    real(kind=real64) :: mu                              !! Shear modulus
    real(kind=real64) :: nu                              !! Poisson's ratio
    real(kind=real64) :: m(fbem_n_nodes(type_f1),3,3)    !! m kernel vector
    real(kind=real64) :: l(fbem_n_nodes(type_f2),3,3)    !! l kernel vector
    ! Local
    integer           :: ksubtri                          ! Counter variable for subtriangles loop
    integer           :: nsubtri                          ! Number of subtriangles
    integer           :: subtriangle(8)                   ! Vector that contains what subtriangles need to be integrated
    real(kind=real64) :: theta_subtri(2,8)                ! Matrix that contains the angles theta of the subtriangles to be integrated
    real(kind=real64) :: thetap_subtri(2,8)               ! Matrix that contains the angles thetap of the subtriangles to be integrated
    integer           :: ktheta                           ! Counter variable for theta coordinate loop
    integer           :: krho                             ! Counter variable for rho coordinate loop
    integer           :: kphi                             ! Counter coordinates for shape functions loops
    integer           :: nnodes_g                         ! Number of nodes of geometrical interpolation
    integer           :: nnodes_f1                        ! Number of nodes of primary variables interpolation
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
    real(kind=real64) :: phi_f1(fbem_n_nodes(type_f1))    ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phi_f2(fbem_n_nodes(type_f2))    ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phi_f1_i(fbem_n_nodes(type_f1))  ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: psi_i(3,fbem_n_nodes(type_f1))   ! psi vectors at xi_i1,xi_i2
    real(kind=real64) :: phi_f2_i(fbem_n_nodes(type_f2))  ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phi_g(fbem_n_nodes(type_g))      ! Geometrical shape functions values at xi_1,xi_2
    real(kind=real64) :: dphidxi1_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64) :: dphidxi2_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64) :: jxi1, jxi2                       ! xi1->x, xi2->x jacobians
    real(kind=real64) :: x_i(3)                           ! Collocation point position vector
    real(kind=real64) :: n_i(3)                           ! Collocation point unit normal vector
    real(kind=real64) :: x(3)                             ! Position vector at xi_1,xi_2
    real(kind=real64) :: T1(3), T2(3)                     ! Tangent vectors at xi_1,xi_2
    real(kind=real64) :: v1(3), v2(3)                     ! Orthogonal tangent vectors
    real(kind=real64) :: Mvt(2,2)                         ! Orthogonal coordinates transformation matrix
    real(kind=real64) :: dxidv(2,2)                       ! xi coordinates derivatives with respect to v orthogonal coordinates
    real(kind=real64) :: detMvt                           ! Mvt(1,1)*Mvt(2,2)-Mvt(1,2)*Mvt(2,1)
    real(kind=real64) :: N(3)                             ! Normal vector at xi_1,xi_2
    real(kind=real64) :: rv(3)                            ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64) :: r, d1r, d1r2, d1r3               ! Distance vector module and its inverse
    real(kind=real64) :: drdx(3)                          ! Distance vector derivatives with respect to x_k
    real(kind=real64) :: jg                               ! Geometric jacobian
    real(kind=real64) :: jw                               ! Jacobian * weights
    real(kind=real64) :: drdn                             ! Partial derivative of r respect to unit normal
    real(kind=real64) :: drdni                            ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64) :: n_dot_ni                         ! Dot product n · n_i
    real(kind=real64) :: phif1jw(fbem_n_nodes(type_f1))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phif2jw(fbem_n_nodes(type_f2))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: costheta, sintheta               ! cos(theta), sin(theta)
    integer           :: il, ik, ij                       ! Counter for load / observation components
    real(kind=real64) :: cted1, cted2                     ! Auxiliary constants
    real(kind=real64) :: ctes1, ctes2, ctes3              ! Auxiliary constants
    real(kind=real64) :: fs_d, fs_s                       ! Fundamental solutions values
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
    real(kind=real64)              :: lli(3,3)                             ! Line integral epsilon_{ijk}·e_k·t/r
    real(kind=real64)              :: mli1                                 ! Line integral (r x ni)·t/r**3
    real(kind=real64)              :: mli2(3)                              ! Line integral (e_k x ni)·t/r
    real(kind=real64)              :: mli3(3,3)                            ! Line integral (r x ni)·t*r_{,l}*r_{,k}/r**3
    real(kind=real64)              :: mli4(3)                              ! Line integral (r x e_k)·t/r**3
    real(kind=real64)              :: mli5(3,3,3)                          ! Line integral (e_k x ni)·t*r_{,l}*r_{,j}/r
    !
    ! Initialization
    !
    ! Kernels
    m=0.d0
    l=0.d0
    ! Auxiliary constants for fundamental solutions calculation
    cted1=-1./(8.*c_pi*(1.-nu))
    cted2=1.-2.*nu
    ctes1=mu/(4.*c_pi*(1.-nu))
    ctes2=1.-2.*nu
    ctes3=1.-4.*nu
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
          d1r=1.d0/r
          d1r2=d1r**2
          d1r3=d1r2*d1r
          ! r_{,k}
          drdx=rv*d1r
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
          ! Loop through load direction and observation direction
          do il=1,3
            do ik=1,3
              !
              ! L integrals
              !
              ! L1
              fs_d=d1r2*((cted2*c_dkr(il,ik)+3.*drdx(il)*drdx(ik))*drdni+cted2*((n_i(il)-n(il))*drdx(ik)-(n_i(ik)-n(ik))*drdx(il)))
              l(:,il,ik)=l(:,il,ik)+fs_d*phif2jw
              ! The regular parts of the CPV integral
              ! L21
              if (il.ne.ik) then
                fs_d=d1r2*cted2*(n(il)*drdx(ik)-n(ik)*drdx(il))
                l(:,il,ik)=l(:,il,ik)+fs_d*(phi_f2-phi_f2_i)*jw
              end if
              !
              ! M integrals
              !
              ! M1
              fs_s=3.d0*d1r3*drdn*drdni*(5.d0*drdx(il)*drdx(ik)-nu*c_dkr(il,ik))
              m(:,il,ik)=m(:,il,ik)+fs_s*phif1jw
              ! Regular parts of M2
              if (il.eq.ik) then
                ! M21
                fs_s=ctes2*d1r3*n_dot_ni
                m(:,il,ik)=m(:,il,ik)+fs_s*(phi_f1-phi_f1_i-psi_i(1,:)*rv(1)-psi_i(2,:)*rv(2)-psi_i(3,:)*rv(3))*jw
                ! M221
                fs_s=-ctes2*3.d0*d1r3*drdn*drdni
                m(:,il,ik)=m(:,il,ik)+phi_f1_i*fs_s*jw
                ! M231
                fs_s=-ctes2*d1r2*drdni
                m(:,il,ik)=m(:,il,ik)+fs_s*(psi_i(1,:)*n(1)+psi_i(2,:)*n(2)+psi_i(3,:)*n(3))*jw
              end if
              ! Regular parts of M3
              ! M31
              fs_s=d1r3*(3.d0*ctes2*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni)&
                        +3.d0*nu*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni)&
                        +3.d0*nu*drdx(il)*drdx(ik)*n_dot_ni&
                        +ctes2*n_i(ik)*n(il)-ctes3*n_i(il)*n(ik))
              m(:,il,ik)=m(:,il,ik)+fs_s*(phi_f1-phi_f1_i-psi_i(1,:)*rv(1)-psi_i(2,:)*rv(2)-psi_i(3,:)*rv(3))*jw
              ! M321
              fs_s=d1r3*(-15.d0*nu*drdx(il)*drdx(ik)*drdn*drdni&
                         +3.d0*nu*drdx(ik)*(n_i(il)*drdn-n(il)*drdni)&
                         +3.d0*ctes2*drdx(il)*(n_i(ik)*drdn-n(ik)*drdni))
              m(:,il,ik)=m(:,il,ik)+phi_f1_i*fs_s*jw
              ! M331
              do ij=1,3
                fs_s=d1r2*(drdx(ij)*(3.d0*ctes2*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni)&
                                    +3.d0*nu*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni)&
                                    -3.d0*nu*drdx(il)*n(ik)*drdni)&
                          -((1.d0-3.d0*nu)*n_i(il)*drdx(ij)+nu*n_i(ij)*drdx(il))*(n(ik)-n_i(ik)*n_dot_ni)&
                          +ctes2*n_i(ik)*drdx(ij)*(n(il)-n_i(il)*n_dot_ni)&
                          -nu*(c_dkr(il,ik)+n_i(il)*n_i(ik))*n(ij)*drdni&
                          -nu*(c_dkr(ij,ik)-n_i(ij)*n_i(ik))*n(il)*drdni)
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
        l(:,il,ik)=l(:,il,ik)+phi_f2_i*cted2*lli(il,ik)
        if (il.eq.ik) then
          ! M222
          m(:,il,ik)=m(:,il,ik)+phi_f1_i*ctes2*mli1
          ! M232
          do ij=1,3
            m(:,il,ik)=m(:,il,ik)+psi_i(ij,:)*ctes2*mli2(ij)
          end do
        end if
        ! M322
        m(:,il,ik)=m(:,il,ik)+phi_f1_i*3.d0*nu*mli3(il,ik)
        ! M323
        m(:,il,ik)=m(:,il,ik)+phi_f1_i*ctes2*n_i(ik)*mli4(il)
        ! M324
        m(:,il,ik)=m(:,il,ik)-phi_f1_i*ctes3*n_i(il)*mli4(ik)
        ! M33X
        do ij=1,3
          ! M332
          m(:,il,ik)=m(:,il,ik)+psi_i(ij,:)*nu*mli5(ij,il,ik)
          ! M333
          m(:,il,ik)=m(:,il,ik)+psi_i(ij,:)*nu*(c_dkr(il,ik)+n_i(il)*n_i(ik))*mli2(ij)
          ! M334
          m(:,il,ik)=m(:,il,ik)+psi_i(ij,:)*nu*(c_dkr(ij,ik)-n_i(ij)*n_i(ik))*mli2(il)
        end do
      end do
    end do
    ! Multiply m by ctes1 and l by cted1
    m=ctes1*m
    l=cted1*l
    ! If the normal has to be reversed, then l=-l
    if (reverse) l=-l
  end subroutine fbem_bem_staela3d_hbie_int


  !! This subroutine calculates adaptatively all line integrals for the HBIE interior integration
  !! using Telles transformation and subdivision. Before calling this subroutine, mli3 must be set to zero.
  recursive subroutine fbem_bem_staela3d_hbie_int_li(gtype,x_nodes,xi_s,x_i,n_i,ngp_min,qsp,ks,ns,lli,mli1,mli2,mli3,mli4,mli5)
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
    real(kind=real64)            :: lli(3,3)                        !! Line integral epsilon_{ijk}*e_k·t/r
    real(kind=real64)            :: mli1                            !! Line integral (r x ni)·t/r**3
    real(kind=real64)            :: mli2(3)                         !! Line integral (e_k x ni)·t/r
    real(kind=real64)            :: mli3(3,3)                       !! Line integrals (r x ni)·t*r_{,l}*r_{,k}/r**3
    real(kind=real64)            :: mli4(3)                         !! Line integrals (r x e_k)·t/r**3
    real(kind=real64)            :: mli5(3,3,3)                     !! Line integrals (e_k x ni)·t*r_{,l}*r_{,j}/r
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
    real(kind=real64)            :: lli_tmp(3,3)                    ! Line integral epsilon_{ijk}*e_k·t/r
    real(kind=real64)            :: mli1_tmp                        ! Line integral (r x ni)·t/r**3
    real(kind=real64)            :: mli2_tmp(3)                     ! Line integral (e_k x ni)·t/r
    real(kind=real64)            :: mli3_tmp(3,3)                   ! Line integrals (r x ni)·t*r_{,l}*r_{,k}/r**3
    real(kind=real64)            :: mli4_tmp(3)                     ! Line integrals (r x e_k)·t/r**3
    real(kind=real64)            :: mli5_tmp(3,3,3)                 ! Line integrals (e_k x ni)·t*r_{,l}*r_{,j}/r
    ! Local calculation
    integer                      :: kphi                            ! Counter variable for shape functions loops
    integer                      :: n_gnodes                        ! Number of nodes of the geometrical element
    integer                      :: kip                             ! Counter variable for integration points
    integer                      :: ik                              ! Coordinate counter
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
    real(kind=real64)            :: x(3), T(3), rv(3), r, dr1, dr3  ! Position, tangent, distance vector, distance and 1/r at xi
    real(kind=real64)            :: drdx(3)                         ! dr/dx
    real(kind=real64)            :: ek(3)                           ! e_k vector
    real(kind=real64)            :: rxni(3)                         ! r x ni
    real(kind=real64)            :: rxnit                           ! (r x n_i)·t
    real(kind=real64)            :: rxek(3)                         ! r x e_k
    real(kind=real64)            :: rxekt                           ! (r x e_k)·t
    real(kind=real64)            :: ekxni(3)                        ! e_k x ni
    real(kind=real64)            :: ekxnit                          ! (e_k x ni)·t
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
        call fbem_warning_message(error_unit,0,'fbem_bem_staela3d_hbie_int_mli',ns,'maximum number of subdivisions reached')
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
      call fbem_bem_staela3d_hbie_int_li(gtype,x_nodes,tmp_xi_s,x_i,n_i,ngp_min,qsp,ks+1,ns,lli,mli1,mli2,mli3,mli4,mli5)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_staela3d_hbie_int_li(gtype,x_nodes,tmp_xi_s,x_i,n_i,ngp_min,qsp,ks+1,ns,lli,mli1,mli2,mli3,mli4,mli5)
    ! Calculate the subdivided element using Telles transformation
    else
      ! Initialization
      n_gnodes=fbem_n_nodes(gtype)
      lli_tmp=0.d0
      mli1_tmp=0.d0
      mli2_tmp=0.d0
      mli3_tmp=0.d0
      mli4_tmp=0.d0
      mli5_tmp=0.d0
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
        dr3=1.d0/r**3
        ! r_{,k}
        drdx=rv*dr1
        ! Jacobians * weight
        jw=jg*js*jt*w
        !
        ! lli
        !
        ! Add integration points
        lli_tmp(1,2)=lli_tmp(1,2)-dr1*t(3)*jw
        lli_tmp(1,3)=lli_tmp(1,3)+dr1*t(2)*jw
        lli_tmp(2,3)=lli_tmp(2,3)-dr1*t(1)*jw
        !
        ! mli1 and mli3
        !
        ! r x n_i
        rxni(1)=rv(2)*n_i(3)-rv(3)*n_i(2)
        rxni(2)=rv(3)*n_i(1)-rv(1)*n_i(3)
        rxni(3)=rv(1)*n_i(2)-rv(2)*n_i(1)
        ! (r x n_i)·t
        rxnit=dot_product(rxni,t)
        ! Add
        mli1_tmp=mli1_tmp+dr3*rxnit*jw
        mli3_tmp(1,1)=mli3_tmp(1,1)+dr3*rxnit*drdx(1)*drdx(1)*jw
        mli3_tmp(1,2)=mli3_tmp(1,2)+dr3*rxnit*drdx(1)*drdx(2)*jw
        mli3_tmp(1,3)=mli3_tmp(1,3)+dr3*rxnit*drdx(1)*drdx(3)*jw
        mli3_tmp(2,2)=mli3_tmp(2,2)+dr3*rxnit*drdx(2)*drdx(2)*jw
        mli3_tmp(2,3)=mli3_tmp(2,3)+dr3*rxnit*drdx(2)*drdx(3)*jw
        mli3_tmp(3,3)=mli3_tmp(3,3)+dr3*rxnit*drdx(3)*drdx(3)*jw
        !
        ! mli2, mli4 and mli5
        !
        ! Loop through coordinates
        do ik=1,3
          ! ek
          ek=0.d0
          ek(ik)=1.d0
          ! e_k x ni
          ekxni(1)=ek(2)*n_i(3)-ek(3)*n_i(2)
          ekxni(2)=ek(3)*n_i(1)-ek(1)*n_i(3)
          ekxni(3)=ek(1)*n_i(2)-ek(2)*n_i(1)
          ! (e_k x ni)·t
          ekxnit=dot_product(ekxni,t)
          ! r x e_k
          rxek(1)=rv(2)*ek(3)-rv(3)*ek(2)
          rxek(2)=rv(3)*ek(1)-rv(1)*ek(3)
          rxek(3)=rv(1)*ek(2)-rv(2)*ek(1)
          ! (r x e_k)·t
          rxekt=dot_product(rxek,t)
          ! Add
          mli2_tmp(ik)=mli2_tmp(ik)+dr1*ekxnit*jw
          mli4_tmp(ik)=mli4_tmp(ik)+dr3*rxekt*jw
          mli5_tmp(1,1,ik)=mli5_tmp(1,1,ik)+dr1*ekxnit*drdx(1)*drdx(1)*jw
          mli5_tmp(1,2,ik)=mli5_tmp(1,2,ik)+dr1*ekxnit*drdx(1)*drdx(2)*jw
          mli5_tmp(1,3,ik)=mli5_tmp(1,3,ik)+dr1*ekxnit*drdx(1)*drdx(3)*jw
          mli5_tmp(2,2,ik)=mli5_tmp(2,2,ik)+dr1*ekxnit*drdx(2)*drdx(2)*jw
          mli5_tmp(2,3,ik)=mli5_tmp(2,3,ik)+dr1*ekxnit*drdx(2)*drdx(3)*jw
          mli5_tmp(3,3,ik)=mli5_tmp(3,3,ik)+dr1*ekxnit*drdx(3)*drdx(3)*jw
        end do
      end do
      ! Copy antisymmetrical parts
      lli_tmp(2,1)=-lli_tmp(1,2)
      lli_tmp(3,1)=-lli_tmp(1,3)
      lli_tmp(3,2)=-lli_tmp(2,3)
      ! Copy symmetrical parts
      mli3_tmp(2,1)=mli3_tmp(1,2)
      mli3_tmp(3,1)=mli3_tmp(1,3)
      mli3_tmp(3,2)=mli3_tmp(2,3)
      mli5_tmp(2,1,:)=mli5_tmp(1,2,:)
      mli5_tmp(3,1,:)=mli5_tmp(1,3,:)
      mli5_tmp(3,2,:)=mli5_tmp(2,3,:)
      ! Add
      lli=lli+lli_tmp
      mli1=mli1+mli1_tmp
      mli2=mli2+mli2_tmp
      mli3=mli3+mli3_tmp
      mli4=mli4+mli4_tmp
      mli5=mli5+mli5_tmp
    end if
  end subroutine fbem_bem_staela3d_hbie_int_li

  !! This subroutine calculates adaptatively the line integrals (r x ni)·t*r_{,l}*r_{,k}/r**3 for the HBIE interior integration
  !! using Telles transformation and subdivision. Before calling this subroutine, mli3 must be set to zero.
  recursive subroutine fbem_bem_staela3d_hbie_int_mli3(gtype,x_nodes,xi_s,x_i,n_i,ngp_min,qsp,ks,ns,mli3)
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
    real(kind=real64)            :: mli3(3,3)                       !! Line integrals (r x ni)·t*r_{,l}*r_{,k}/r**3
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
    real(kind=real64)            :: mli3_tmp(3,3)                   ! Line integrals (r x ni)·t*r_{,l}*r_{,k}/r**3 (temporary)
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
    real(kind=real64)            :: x(3), T(3), rv(3), r, dr1, dr3  ! Position, tangent, distance vector, distance and 1/r at xi
    real(kind=real64)            :: drdx(3)                         ! dr/dx
    real(kind=real64)            :: rxni(3)                         ! r x n_i
    real(kind=real64)            :: rxnit                           ! (r x n_i)·t
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
        call fbem_warning_message(error_unit,0,'fbem_bem_staela3d_hbie_int_mli3',ns,'maximum number of subdivisions reached')
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
      call fbem_bem_staela3d_hbie_int_mli3(gtype,x_nodes,tmp_xi_s,x_i,n_i,ngp_min,qsp,ks+1,ns,mli3)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_staela3d_hbie_int_mli3(gtype,x_nodes,tmp_xi_s,x_i,n_i,ngp_min,qsp,ks+1,ns,mli3)
    ! Calculate the subdivided element using Telles transformation
    else
      ! Initialization
      n_gnodes=fbem_n_nodes(gtype)
      mli3_tmp=0.d0
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
        dr3=1.d0/r**3
        ! r_{,k}
        drdx=rv*dr1
        ! r x n_i
        rxni(1)=rv(2)*n_i(3)-rv(3)*n_i(2)
        rxni(2)=rv(3)*n_i(1)-rv(1)*n_i(3)
        rxni(3)=rv(1)*n_i(2)-rv(2)*n_i(1)
        ! (r x n_i)·t
        rxnit=dot_product(rxni,t)
        ! Jacobians * weight
        jw=jg*js*jt*w
        ! Add integration points
        mli3_tmp(1,1)=mli3_tmp(1,1)+dr3*rxnit*drdx(1)*drdx(1)*jw
        mli3_tmp(1,2)=mli3_tmp(1,2)+dr3*rxnit*drdx(1)*drdx(2)*jw
        mli3_tmp(1,3)=mli3_tmp(1,3)+dr3*rxnit*drdx(1)*drdx(3)*jw
        mli3_tmp(2,2)=mli3_tmp(2,2)+dr3*rxnit*drdx(2)*drdx(2)*jw
        mli3_tmp(2,3)=mli3_tmp(2,3)+dr3*rxnit*drdx(2)*drdx(3)*jw
        mli3_tmp(3,3)=mli3_tmp(3,3)+dr3*rxnit*drdx(3)*drdx(3)*jw
      end do ! Loop through gamma coordinate
      ! Copy symmetrical part
      mli3_tmp(2,1)=mli3_tmp(1,2)
      mli3_tmp(3,1)=mli3_tmp(1,3)
      mli3_tmp(3,2)=mli3_tmp(2,3)
      ! Add mli3
      mli3=mli3+mli3_tmp
    end if
  end subroutine fbem_bem_staela3d_hbie_int_mli3

!  !! This subroutine calculates adaptatively the line integrals (r x ni)·t*r_{,l}*r_{,k}/r**3 for the HBIE interior integration
!  !! using Telles transformation and subdivision. Before calling this subroutine, mli1 must be set to zero, and xi1=-1.0d0 and
!  !! xi2=1.0d0
!  recursive subroutine fbem_bem_staela3d_hbie_int_mli3(type_g,x_nodes,xi1,xi2,x_i,n_i,ngp_min,re,mli3)
!    implicit none
!    ! I/O
!    integer           :: type_g                          !! Geometrial interpolation
!    real(kind=real64) :: x_nodes(3,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
!    real(kind=real64) :: xi1                             !! Minimum coordinate xi of the subdivision
!    real(kind=real64) :: xi2                             !! Maximum coordinate xi of the subdivision
!    real(kind=real64) :: x_i(3)                          !! Collocation point position vector
!    real(kind=real64) :: n_i(3)                          !! Collocation point unit normal vector
!    integer           :: ngp_min                         !! Minimum number of Gauss points
!    real(kind=real64) :: re                              !! Integration relative error
!    real(kind=real64) :: mli3(3,3)                       !! Line integrals (r x ni)·t*r_{,l}*r_{,k}/r**3
!    ! Local
!    real(kind=real64) :: xise(2,4)                       ! xi coordinates of the subdivision
!    integer           :: ngp                             ! Number of Gauss points
!    logical           :: achieve                         ! True if the calculated ngp achieve the relative error
!    real(kind=real64) :: barxi                           ! Nearest element coordinate with respect to collocation point
!    real(kind=real64) :: r_min                           ! Minimum distance between collocation point and barxi on the element
!    real(kind=real64) :: barr                            ! Telles jacobian at barxi
!    real(kind=real64) :: csize                           ! Characteristic size
!    real(kind=real64) :: d                               ! Normalized distance between collocation point and element subdivision
!    real(kind=real64) :: barxisub                        ! Nearest subdivision coordinate with respect to the collocation point
!    ! Local calculation
!    real(kind=real64)            :: mli3_tmp(3,3)                   ! Line integrals (r x ni)·t*r_{,l}*r_{,k}/r**3 (temporary)
!    integer                      :: kphi                            ! Counter variable for shape functions loops
!    integer                      :: nnodes_g                        ! Number of nodes of the geometrical element
!    integer                      :: kip                             ! Counter variable for integration points
!    integer                      :: il, ik                          ! Counters
!    real(kind=real64)            :: gamma                           ! Reference coordinate
!    type(fbem_telles_parameters) :: telles_parameters               ! Telles parameters
!    real(kind=real64)            :: jt                              ! Telles jacobian for each coordinate: xi->gamma
!    real(kind=real64)            :: xip                             ! Reference coordinate of the subdivision
!    real(kind=real64)            :: barxip                          ! Collocation point in xip space
!    real(kind=real64)            :: js                              ! Jacobian of xip->xi transformation
!    real(kind=real64)            :: xi                              ! Reference coordinate
!    real(kind=real64)            :: w                               ! Weights of each integration point
!    real(kind=real64)            :: aux(10)                         ! Auxiliary variable needed for shape_functions module resources
!    real(kind=real64)            :: phi_g(fbem_n_nodes(type_g))     ! Geometrical shape functions values
!    real(kind=real64)            :: dphidxi_g(fbem_n_nodes(type_g)) ! Geometrical shape functions first derivatives values
!    real(kind=real64)            :: x(3)                            ! Position vector at xi
!    real(kind=real64)            :: T(3)                            ! Tangent vector at xi
!    real(kind=real64)            :: rv(3)                           ! Distance vector between collocation point and integration point (x-x_i)
!    real(kind=real64)            :: drdx(3)                         ! dr/dx
!    real(kind=real64)            :: r, d1r, d1r3                    ! Distance vector module and its inverse
!    real(kind=real64)            :: rxni(3)                         ! r x n_i
!    real(kind=real64)            :: rxnit                           ! (r x n_i)·t
!    real(kind=real64)            :: jw                              ! Jacobians * weight
!    real(kind=real64)            :: jg                              ! Geometric jacobian
!    !
!    ! Integration setup
!    !
!    ! Find the minimum distance between the element subdivision and the collocation point, and the nearest element coordinate in
!    ! element reference space [-1,1]
!    call fbem_nearest_xi_subdivision(3,type_g,x_nodes,xi1,xi2,x_i,25,14,2,barxi,r_min)
!    ! Obtain the characteristic size of the element subdivision
!    xise(1,1)=xi1
!    xise(1,2)=xi2
!    csize=fbem_characteristic_length_subdivision(3,type_g,x_nodes,xise,1.d-6)
!    ! Calculate the normalized minimum distance
!    d=r_min/csize
!    ! Obtain the barxi but in the element subdivision space (barxisub) [-1,1]
!    barxisub=2.0d0*(barxi-xi1)/(xi2-xi1)-1.0d0
!    ! Obtain the estimated number of Gaussian points
!    call fbem_qs_f_ngp_telles_1d(type_g,fbem_f_1r2,barxisub,d,re,ngp,achieve)
!    !
!    ! If not achieve the indicated relative error, then subdivide by 1/2
!    !
!    if (achieve.eqv.(.false.)) then
!      call fbem_bem_staela3d_hbie_int_mli3(type_g,x_nodes,xi1,0.5d0*(xi1+xi2),x_i,n_i,ngp_min,re,mli3)
!      call fbem_bem_staela3d_hbie_int_mli3(type_g,x_nodes,0.5d0*(xi1+xi2),xi2,x_i,n_i,ngp_min,re,mli3)
!    !
!    ! If achieve, then integrate the present subdivision
!    !
!    else
!      ! Obtain the optimum Telles jacobian
!      barr=fbem_telles_barr(d,fbem_f_any)
!      ! Use ngp_min if ngp is less than ngp_min
!      if (ngp.lt.ngp_min) ngp=ngp_min
!      ! Integrals initialization
!      mli3_tmp=0.d0
!      ! Number of nodes of geometrical interpolation
!      nnodes_g=fbem_n_nodes(type_g)
!      ! Subdivision jacobian (constant)
!      js=xi2-xi1
!      ! Barxi in subdivision reference space
!      barxip=(barxi-xi1)/js
!      ! Calculate Telles parameters
!      telles_parameters=fbem_telles01_calculate_parameters(barxip,barr)
!      ! Loop through gamma coordinate
!      do kip=1,gl01_n(ngp)
!        ! GAMMA->XIP->XI COORDINATE TRANSFORMATION
!        ! Gamma coordinate and weight
!        gamma=gl01_xi(kip,ngp)
!        w=gl01_w(kip,ngp)
!        ! xip coordinate, weight and jacobian from Telles transformation
!        call fbem_telles_xi_and_jacobian(telles_parameters,gamma,xip,jt)
!        ! xi coordinate
!        xi=js*xip+xi1
!        ! COMPONENTS OF THE FUNDAMENTAL SOLUTIONS
!        ! Geometrical shape functions and first derivatives at xi
!#       define etype type_g
!#       define delta 0.0d0
!#       define phi phi_g
!#       define dphidxi dphidxi_g
!#       include <phi_and_dphidxi_1d.rc>
!#       undef etype
!#       undef delta
!#       undef phi
!#       undef dphidxi
!        ! Components calculation of x and T at xi
!        x=0.d0
!        T=0.d0
!        do kphi=1,nnodes_g
!          x=x+phi_g(kphi)*x_nodes(:,kphi)
!          T=T+dphidxi_g(kphi)*x_nodes(:,kphi)
!        end do
!        ! Geometric jacobian
!        jg=dot_product(T,T)
!        jg=sqrt(jg)
!        ! Unit tangent
!        t=T/jg
!        ! Distance vector
!        rv=x-x_i
!        ! Distance vector norm
!        r=dot_product(rv,rv)
!        r=sqrt(r)
!        d1r=1.d0/r
!        d1r3=d1r**3
!        ! r_{,k}
!        drdx=rv*d1r
!        ! Jacobians * weight
!        jw=jg*js*jt*w
!        ! r x n_i
!        rxni(1)=rv(2)*n_i(3)-rv(3)*n_i(2)
!        rxni(2)=rv(3)*n_i(1)-rv(1)*n_i(3)
!        rxni(3)=rv(1)*n_i(2)-rv(2)*n_i(1)
!        ! (r x n_i)·t
!        rxnit=dot_product(rxni,t)
!        ! Add integration points
!        do il=1,3
!          do ik=1,3
!            mli3_tmp(il,ik)=mli3_tmp(il,ik)+d1r3*rxnit*drdx(il)*drdx(ik)*jw
!          end do
!        end do
!      end do ! Loop through gamma coordinate
!      ! Add mli3
!      mli3=mli3+mli3_tmp
!    end if
!  end subroutine fbem_bem_staela3d_hbie_int_mli3

  !! This subroutine calculates adaptatively the line integrals (r x e_k)·t/r**3 for the HBIE interior integration
  !! using Telles transformation and subdivision. Before calling this subroutine, mli4 must be set to zero.
  recursive subroutine fbem_bem_staela3d_hbie_int_mli4(gtype,x_nodes,xi_s,x_i,ngp_min,qsp,ks,ns,mli4)
    implicit none
    ! I/O
    integer                      :: gtype                           !! Geometrial interpolation
    real(kind=real64)            :: x_nodes(3,fbem_n_nodes(gtype))  !! Position vectors of geometrical nodes
    real(kind=real64)            :: xi_s(1,2)                       !! Subdivision of the parent element
    real(kind=real64)            :: x_i(3)                          !! Collocation point position vector
    integer                      :: ngp_min                         !! Minimum number of Gauss points
    type(fbem_qs_parameters)     :: qsp                             !! Quasi-singular integration parameters
    integer                      :: ks                              !! Current level of subdivisions
    integer                      :: ns                              !! Maximum level of subdivision
    real(kind=real64)            :: mli4(3)                         !! Line integrals (r x e_k)·t/r**3
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
    real(kind=real64)            :: mli4_tmp(3)                     ! Line integrals (r x e_k)·t/r**3 (temporary)
    ! Local calculation
    integer                      :: kphi                            ! Counter variable for shape functions loops
    integer                      :: n_gnodes                        ! Number of nodes of the geometrical element
    integer                      :: kip                             ! Counter variable for integration points
    integer                      :: ik
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
    real(kind=real64)            :: x(3), T(3), rv(3), r, dr3       ! Position, tangent, distance vector, distance and 1/r at xi
    real(kind=real64)            :: ek(3)                           ! e_k vector
    real(kind=real64)            :: rxek(3)                         ! r x e_k
    real(kind=real64)            :: rxekt                           ! (r x e_k)·t
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
        call fbem_warning_message(error_unit,0,'fbem_bem_staela3d_hbie_int_mli4',ns,'maximum number of subdivisions reached')
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
      call fbem_bem_staela3d_hbie_int_mli4(gtype,x_nodes,tmp_xi_s,x_i,ngp_min,qsp,ks+1,ns,mli4)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_staela3d_hbie_int_mli4(gtype,x_nodes,tmp_xi_s,x_i,ngp_min,qsp,ks+1,ns,mli4)
    ! Calculate the subdivided element using Telles transformation
    else
      ! Initialization
      n_gnodes=fbem_n_nodes(gtype)
      mli4_tmp=0.d0
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
        dr3=1.d0/r**3
        ! Jacobians * weight
        jw=jg*js*jt*w
        ! Loop through components
        do ik=1,3
          ! e_k
          ek=0.d0
          ek(ik)=1.d0
          ! r x e_k
          rxek(1)=rv(2)*ek(3)-rv(3)*ek(2)
          rxek(2)=rv(3)*ek(1)-rv(1)*ek(3)
          rxek(3)=rv(1)*ek(2)-rv(2)*ek(1)
          ! (r x e_k)·t
          rxekt=dot_product(rxek,t)
          ! Add integration points
          mli4_tmp(ik)=mli4_tmp(ik)+dr3*rxekt*jw
        end do
      end do
      ! Add mli4
      mli4=mli4+mli4_tmp
    end if
  end subroutine fbem_bem_staela3d_hbie_int_mli4

!  !! This subroutine calculates adaptatively the line integrals (r x e_k)·t/r**3 for the HBIE interior integration
!  !! using Telles transformation and subdivision. Before calling this subroutine, mli1 must be set to zero, and xi1=-1.0d0 and
!  !! xi2=1.0d0
!  recursive subroutine fbem_bem_staela3d_hbie_int_mli4(type_g,x_nodes,xi1,xi2,x_i,ngp_min,re,mli4)
!    implicit none
!    ! I/O
!    integer           :: type_g                          !! Geometrial interpolation
!    real(kind=real64) :: x_nodes(3,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
!    real(kind=real64) :: xi1                             !! Minimum coordinate xi of the subdivision
!    real(kind=real64) :: xi2                             !! Maximum coordinate xi of the subdivision
!    real(kind=real64) :: x_i(3)                          !! Collocation point position vector
!    integer           :: ngp_min                         !! Minimum number of Gauss points
!    real(kind=real64) :: re                              !! Integration relative error
!    real(kind=real64) :: mli4(3)                         !! Line integrals (r x e_k)·t/r**3
!    ! Local
!    real(kind=real64) :: xise(2,4)                       ! xi coordinates of the subdivision
!    integer           :: ngp                             ! Number of Gauss points
!    logical           :: achieve                         ! True if the calculated ngp achieve the relative error
!    real(kind=real64) :: barxi                           ! Nearest element coordinate with respect to collocation point
!    real(kind=real64) :: r_min                           ! Minimum distance between collocation point and barxi on the element
!    real(kind=real64) :: barr                            ! Telles jacobian at barxi
!    real(kind=real64) :: csize                           ! Characteristic size
!    real(kind=real64) :: d                               ! Normalized distance between collocation point and element subdivision
!    real(kind=real64) :: barxisub                        ! Nearest subdivision coordinate with respect to the collocation point
!    ! Local calculation
!    real(kind=real64)            :: mli4_tmp(3)                     ! Line integrals (r x e_k)·t/r**3 (temporary)
!    integer                      :: kphi                            ! Counter variable for shape functions loops
!    integer                      :: nnodes_g                        ! Number of nodes of the geometrical element
!    integer                      :: kip                             ! Counter variable for integration points
!    integer                      :: ik                              ! Counter
!    real(kind=real64)            :: gamma                           ! Reference coordinate
!    type(fbem_telles_parameters) :: telles_parameters               ! Telles parameters
!    real(kind=real64)            :: jt                              ! Telles jacobian for each coordinate: xi->gamma
!    real(kind=real64)            :: xip                             ! Reference coordinate of the subdivision
!    real(kind=real64)            :: barxip                          ! Collocation point in xip space
!    real(kind=real64)            :: js                              ! Jacobian of xip->xi transformation
!    real(kind=real64)            :: xi                              ! Reference coordinate
!    real(kind=real64)            :: w                               ! Weights of each integration point
!    real(kind=real64)            :: aux(10)                         ! Auxiliary variable needed for shape_functions module resources
!    real(kind=real64)            :: phi_g(fbem_n_nodes(type_g))     ! Geometrical shape functions values
!    real(kind=real64)            :: dphidxi_g(fbem_n_nodes(type_g)) ! Geometrical shape functions first derivatives values
!    real(kind=real64)            :: x(3)                            ! Position vector at xi
!    real(kind=real64)            :: T(3)                            ! Tangent vector at xi
!    real(kind=real64)            :: rv(3)                           ! Distance vector between collocation point and integration point (x-x_i)
!    real(kind=real64)            :: r, d1r, d1r3                    ! Distance vector module and its inverse
!    real(kind=real64)            :: ek(3)                           ! e_k vector
!    real(kind=real64)            :: rxek(3)                         ! r x e_k
!    real(kind=real64)            :: rxekt                           ! (r x e_k)·t
!    real(kind=real64)            :: jw                              ! Jacobians * weight
!    real(kind=real64)            :: jg                              ! Geometric jacobian
!    !
!    ! Integration setup
!    !
!    ! Find the minimum distance between the element subdivision and the collocation point, and the nearest element coordinate in
!    ! element reference space [-1,1]
!    call fbem_nearest_xi_subdivision(3,type_g,x_nodes,xi1,xi2,x_i,25,14,2,barxi,r_min)
!    ! Obtain the characteristic size of the element subdivision
!    xise(1,1)=xi1
!    xise(1,2)=xi2
!    csize=fbem_characteristic_length_subdivision(3,type_g,x_nodes,xise,1.d-6)
!    ! Calculate the normalized minimum distance
!    d=r_min/csize
!    ! Obtain the barxi but in the element subdivision space (barxisub) [-1,1]
!    barxisub=2.0d0*(barxi-xi1)/(xi2-xi1)-1.0d0
!    ! Obtain the estimated number of Gaussian points
!    call fbem_qs_f_ngp_telles_1d(type_g,fbem_f_1r2,barxisub,d,re,ngp,achieve)
!    !
!    ! If not achieve the indicated relative error, then subdivide by 1/2
!    !
!    if (achieve.eqv.(.false.)) then
!      call fbem_bem_staela3d_hbie_int_mli4(type_g,x_nodes,xi1,0.5d0*(xi1+xi2),x_i,ngp_min,re,mli4)
!      call fbem_bem_staela3d_hbie_int_mli4(type_g,x_nodes,0.5d0*(xi1+xi2),xi2,x_i,ngp_min,re,mli4)
!    !
!    ! If achieve, then integrate the present subdivision
!    !
!    else
!      ! Obtain the optimum Telles jacobian
!      barr=fbem_telles_barr(d,fbem_f_any)
!      ! Use ngp_min if ngp is less than ngp_min
!      if (ngp.lt.ngp_min) ngp=ngp_min
!      ! Integrals initialization
!      mli4_tmp=0.d0
!      ! Number of nodes of geometrical interpolation
!      nnodes_g=fbem_n_nodes(type_g)
!      ! Subdivision jacobian (constant)
!      js=xi2-xi1
!      ! Barxi in subdivision reference space
!      barxip=(barxi-xi1)/js
!      ! Calculate Telles parameters
!      telles_parameters=fbem_telles01_calculate_parameters(barxip,barr)
!      ! Loop through gamma coordinate
!      do kip=1,gl01_n(ngp)
!        ! GAMMA->XIP->XI COORDINATE TRANSFORMATION
!        ! Gamma coordinate and weight
!        gamma=gl01_xi(kip,ngp)
!        w=gl01_w(kip,ngp)
!        ! xip coordinate, weight and jacobian from Telles transformation
!        call fbem_telles_xi_and_jacobian(telles_parameters,gamma,xip,jt)
!        ! xi coordinate
!        xi=js*xip+xi1
!        ! COMPONENTS OF THE FUNDAMENTAL SOLUTIONS
!        ! Geometrical shape functions and first derivatives at xi
!#       define etype type_g
!#       define delta 0.0d0
!#       define phi phi_g
!#       define dphidxi dphidxi_g
!#       include <phi_and_dphidxi_1d.rc>
!#       undef etype
!#       undef delta
!#       undef phi
!#       undef dphidxi
!        ! Components calculation of x and T at xi
!        x=0.d0
!        T=0.d0
!        do kphi=1,nnodes_g
!          x=x+phi_g(kphi)*x_nodes(:,kphi)
!          T=T+dphidxi_g(kphi)*x_nodes(:,kphi)
!        end do
!        ! Geometric jacobian
!        jg=dot_product(T,T)
!        jg=sqrt(jg)
!        ! Unit tangent
!        t=T/jg
!        ! Distance vector
!        rv=x-x_i
!        ! Distance vector norm
!        r=dot_product(rv,rv)
!        r=sqrt(r)
!        d1r=1.d0/r
!        d1r3=d1r**3
!        ! Jacobians * weight
!        jw=jg*js*jt*w
!        ! Loop through components
!        do ik=1,3
!          ! e_k
!          ek=0.d0
!          ek(ik)=1.d0
!          ! r x e_k
!          rxek(1)=rv(2)*ek(3)-rv(3)*ek(2)
!          rxek(2)=rv(3)*ek(1)-rv(1)*ek(3)
!          rxek(3)=rv(1)*ek(2)-rv(2)*ek(1)
!          ! (r x e_k)·t
!          rxekt=dot_product(rxek,t)
!          ! Add integration points
!          mli4_tmp(ik)=mli4_tmp(ik)+d1r3*rxekt*jw
!        end do
!      end do ! Loop through gamma coordinate
!      ! Add mli4
!      mli4=mli4+mli4_tmp
!    end if
!  end subroutine fbem_bem_staela3d_hbie_int_mli4

  !! This subroutine calculates adaptatively the line integrals (e_k x ni)·t*r_{,l}*r_{,j}/r for the HBIE interior integration
  !! using Telles transformation and subdivision. Before calling this subroutine, mli5 must be set to zero.
  recursive subroutine fbem_bem_staela3d_hbie_int_mli5(gtype,x_nodes,xi_s,x_i,n_i,ngp_min,qsp,ks,ns,mli5)
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
    real(kind=real64)            :: mli5(3,3,3)                     !! Line integrals (e_k x ni)·t*r_{,l}*r_{,j}/r
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
    real(kind=real64)            :: mli5_tmp(3,3,3)                 ! Line integrals (e_k x ni)·t*r_{,l}*r_{,j}/r (temporary)
    ! Local calculation
    integer                      :: kphi                            ! Counter variable for shape functions loops
    integer                      :: n_gnodes                        ! Number of nodes of the geometrical element
    integer                      :: kip                             ! Counter variable for integration points
    integer                      :: il, ik, ij                      ! Counters
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
    real(kind=real64)            :: drdx(3)                         ! dr/dx
    real(kind=real64)            :: ek(3)                           ! e_k
    real(kind=real64)            :: ekxni(3)                        ! e_k x ni
    real(kind=real64)            :: ekxnit                          ! (e_k x ni)·t
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
        call fbem_warning_message(error_unit,0,'fbem_bem_staela3d_hbie_int_mli5',ns,'maximum number of subdivisions reached')
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
      call fbem_bem_staela3d_hbie_int_mli5(gtype,x_nodes,tmp_xi_s,x_i,n_i,ngp_min,qsp,ks+1,ns,mli5)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_staela3d_hbie_int_mli5(gtype,x_nodes,tmp_xi_s,x_i,n_i,ngp_min,qsp,ks+1,ns,mli5)
    ! Calculate the subdivided element using Telles transformation
    else
      ! Initialization
      n_gnodes=fbem_n_nodes(gtype)
      mli5_tmp=0.d0
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
        ! r_{,k}
        drdx=rv*dr1
        ! Jacobians * weight
        jw=jg*js*jt*w
        ! Loop through ik components
        do ik=1,3
          ! e_k
          ek=0.d0
          ek(ik)=1.d0
          ! e_k x ni
          ekxni(1)=ek(2)*n_i(3)-ek(3)*n_i(2)
          ekxni(2)=ek(3)*n_i(1)-ek(1)*n_i(3)
          ekxni(3)=ek(1)*n_i(2)-ek(2)*n_i(1)
          ! (e_k x ni)·t
          ekxnit=dot_product(ekxni,t)
          ! Loop through il and ij components
          do il=1,3
            do ij=1,3
              ! Add integration points
              mli5_tmp(ij,il,ik)=mli5_tmp(ij,il,ik)+dr1*ekxnit*drdx(il)*drdx(ij)*jw
            end do
          end do
        end do
      end do ! Loop through gamma coordinate
      ! Add mli5
      mli5=mli5+mli5_tmp
    end if
  end subroutine fbem_bem_staela3d_hbie_int_mli5

!  !! This subroutine calculates adaptatively the line integrals (e_k x ni)·t*r_{,l}*r_{,j}/r for the HBIE interior integration
!  !! using Telles transformation and subdivision. Before calling this subroutine, mli1 must be set to zero, and xi1=-1.0d0 and
!  !! xi2=1.0d0
!  recursive subroutine fbem_bem_staela3d_hbie_int_mli5(type_g,x_nodes,xi1,xi2,x_i,n_i,ngp_min,re,mli5)
!    implicit none
!    ! I/O
!    integer           :: type_g                          !! Geometrial interpolation
!    real(kind=real64) :: x_nodes(3,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
!    real(kind=real64) :: xi1                             !! Minimum coordinate xi of the subdivision
!    real(kind=real64) :: xi2                             !! Maximum coordinate xi of the subdivision
!    real(kind=real64) :: x_i(3)                          !! Collocation point position vector
!    real(kind=real64) :: n_i(3)                          !! Collocation point unit normal vector
!    integer           :: ngp_min                         !! Minimum number of Gauss points
!    real(kind=real64) :: re                              !! Integration relative error
!    real(kind=real64) :: mli5(3,3,3)                     !! Line integrals (e_k x ni)·t*r_{,l}*r_{,j}/r
!    ! Local
!    real(kind=real64) :: xise(2,4)                       ! xi coordinates of the subdivision
!    integer           :: ngp                             ! Number of Gauss points
!    logical           :: achieve                         ! True if the calculated ngp achieve the relative error
!    real(kind=real64) :: barxi                           ! Nearest element coordinate with respect to collocation point
!    real(kind=real64) :: r_min                           ! Minimum distance between collocation point and barxi on the element
!    real(kind=real64) :: barr                            ! Telles jacobian at barxi
!    real(kind=real64) :: csize                           ! Characteristic size
!    real(kind=real64) :: d                               ! Normalized distance between collocation point and element subdivision
!    real(kind=real64) :: barxisub                        ! Nearest subdivision coordinate with respect to the collocation point
!    ! Local calculation
!    real(kind=real64)            :: mli5_tmp(3,3,3)                 ! Line integrals (e_k x ni)·t*r_{,l}*r_{,j}/r (temporary)
!    integer                      :: kphi                            ! Counter variable for shape functions loops
!    integer                      :: nnodes_g                        ! Number of nodes of the geometrical element
!    integer                      :: kip                             ! Counter variable for integration points
!    integer                      :: il, ik, ij                      ! Counters
!    real(kind=real64)            :: gamma                           ! Reference coordinate
!    type(fbem_telles_parameters) :: telles_parameters               ! Telles parameters
!    real(kind=real64)            :: jt                              ! Telles jacobian for each coordinate: xi->gamma
!    real(kind=real64)            :: xip                             ! Reference coordinate of the subdivision
!    real(kind=real64)            :: barxip                          ! Collocation point in xip space
!    real(kind=real64)            :: js                              ! Jacobian of xip->xi transformation
!    real(kind=real64)            :: xi                              ! Reference coordinate
!    real(kind=real64)            :: w                               ! Weights of each integration point
!    real(kind=real64)            :: aux(10)                         ! Auxiliary variable needed for shape_functions module resources
!    real(kind=real64)            :: phi_g(fbem_n_nodes(type_g))     ! Geometrical shape functions values
!    real(kind=real64)            :: dphidxi_g(fbem_n_nodes(type_g)) ! Geometrical shape functions first derivatives values
!    real(kind=real64)            :: x(3)                            ! Position vector at xi
!    real(kind=real64)            :: T(3)                            ! Tangent vector at xi
!    real(kind=real64)            :: rv(3)                           ! Distance vector between collocation point and integration point (x-x_i)
!    real(kind=real64)            :: drdx(3)                         ! dr/dx
!    real(kind=real64)            :: r, d1r                          ! Distance vector module and its inverse
!    real(kind=real64)            :: ek(3)                           ! e_k
!    real(kind=real64)            :: ekxni(3)                        ! e_k x ni
!    real(kind=real64)            :: ekxnit                          ! (e_k x ni)·t
!    real(kind=real64)            :: jw                              ! Jacobians * weight
!    real(kind=real64)            :: jg                              ! Geometric jacobian
!    !
!    ! Integration setup
!    !
!    ! Find the minimum distance between the element subdivision and the collocation point, and the nearest element coordinate in
!    ! element reference space [-1,1]
!    call fbem_nearest_xi_subdivision(3,type_g,x_nodes,xi1,xi2,x_i,25,14,2,barxi,r_min)
!    ! Obtain the characteristic size of the element subdivision
!    xise(1,1)=xi1
!    xise(1,2)=xi2
!    csize=fbem_characteristic_length_subdivision(3,type_g,x_nodes,xise,1.d-6)
!    ! Calculate the normalized minimum distance
!    d=r_min/csize
!    ! Obtain the barxi but in the element subdivision space (barxisub) [-1,1]
!    barxisub=2.0d0*(barxi-xi1)/(xi2-xi1)-1.0d0
!    ! Obtain the estimated number of Gaussian points
!    call fbem_qs_f_ngp_telles_1d(type_g,fbem_f_1r1,barxisub,d,re,ngp,achieve)
!    !
!    ! If not achieve the indicated relative error, then subdivide by 1/2
!    !
!    if (achieve.eqv.(.false.)) then
!      call fbem_bem_staela3d_hbie_int_mli5(type_g,x_nodes,xi1,0.5d0*(xi1+xi2),x_i,n_i,ngp_min,re,mli5)
!      call fbem_bem_staela3d_hbie_int_mli5(type_g,x_nodes,0.5d0*(xi1+xi2),xi2,x_i,n_i,ngp_min,re,mli5)
!    !
!    ! If achieve, then integrate the present subdivision
!    !
!    else
!      ! Obtain the optimum Telles jacobian
!      barr=fbem_telles_barr(d,fbem_f_any)
!      ! Use ngp_min if ngp is less than ngp_min
!      if (ngp.lt.ngp_min) ngp=ngp_min
!      ! Integrals initialization
!      mli5_tmp=0.d0
!      ! Number of nodes of geometrical interpolation
!      nnodes_g=fbem_n_nodes(type_g)
!      ! Subdivision jacobian (constant)
!      js=xi2-xi1
!      ! Barxi in subdivision reference space
!      barxip=(barxi-xi1)/js
!      ! Calculate Telles parameters
!      telles_parameters=fbem_telles01_calculate_parameters(barxip,barr)
!      ! Loop through gamma coordinate
!      do kip=1,gl01_n(ngp)
!        ! GAMMA->XIP->XI COORDINATE TRANSFORMATION
!        ! Gamma coordinate and weight
!        gamma=gl01_xi(kip,ngp)
!        w=gl01_w(kip,ngp)
!        ! xip coordinate, weight and jacobian from Telles transformation
!        call fbem_telles_xi_and_jacobian(telles_parameters,gamma,xip,jt)
!        ! xi coordinate
!        xi=js*xip+xi1
!        ! COMPONENTS OF THE FUNDAMENTAL SOLUTIONS
!        ! Geometrical shape functions and first derivatives at xi
!#       define etype type_g
!#       define delta 0.0d0
!#       define phi phi_g
!#       define dphidxi dphidxi_g
!#       include <phi_and_dphidxi_1d.rc>
!#       undef etype
!#       undef delta
!#       undef phi
!#       undef dphidxi
!        ! Components calculation of x and T at xi
!        x=0.d0
!        T=0.d0
!        do kphi=1,nnodes_g
!          x=x+phi_g(kphi)*x_nodes(:,kphi)
!          T=T+dphidxi_g(kphi)*x_nodes(:,kphi)
!        end do
!        ! Geometric jacobian
!        jg=dot_product(T,T)
!        jg=sqrt(jg)
!        ! Unit tangent
!        t=T/jg
!        ! Distance vector
!        rv=x-x_i
!        ! Distance vector norm
!        r=dot_product(rv,rv)
!        r=sqrt(r)
!        d1r=1.d0/r
!        ! r_{,k}
!        drdx=rv*d1r
!        ! Jacobians * weight
!        jw=jg*js*jt*w
!        ! Loop through ik components
!        do ik=1,3
!          ! e_k
!          ek=0.d0
!          ek(ik)=1.d0
!          ! e_k x ni
!          ekxni(1)=ek(2)*n_i(3)-ek(3)*n_i(2)
!          ekxni(2)=ek(3)*n_i(1)-ek(1)*n_i(3)
!          ekxni(3)=ek(1)*n_i(2)-ek(2)*n_i(1)
!          ! (e_k x ni)·t
!          ekxnit=dot_product(ekxni,t)
!          ! Loop through il and ij components
!          do il=1,3
!            do ij=1,3
!              ! Add integration points
!              mli5_tmp(ij,il,ik)=mli5_tmp(ij,il,ik)+d1r*ekxnit*drdx(il)*drdx(ij)*jw
!            end do
!          end do
!        end do
!      end do ! Loop through gamma coordinate
!      ! Add mli5
!      mli5=mli5+mli5_tmp
!    end if
!  end subroutine fbem_bem_staela3d_hbie_int_mli5

  ! =====================
  ! BE BODY LOAD ELEMENTS
  ! =====================

  !! Efficient automatic integration of BE body load elements
  subroutine fbem_bem_staela3d_hbie_bl_auto(e,x_i,n_i,mu,nu,qsp,ns,lbl)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                   !! Integration element
    real(kind=real64)        :: x_i(3)              !! Collocation point
    real(kind=real64)        :: n_i(3)              !! Collocation point unit normal vector
    real(kind=real64)        :: mu                  !! Shear modulus
    real(kind=real64)        :: nu                  !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp                 !! Quasi-singular integration parameters
    integer                  :: ns                  !! Maximum level of subdivisions
    real(kind=real64)        :: lbl(e%n_snodes,3,3) !! lbl integration kernel
    ! Local
    real(kind=real64)        :: r(3)                               ! Distance vector
    real(kind=real64)        :: rmin                               ! Minimum distance between element and x_i
    real(kind=real64)        :: barxi(e%d)                         ! Local coordinates of the nearest element point with respect to x_i
    real(kind=real64)        :: d                                  ! Dimensionless distance
    integer                  :: delta                              ! Control variable
    real(kind=real64)        :: xi_s(e%d,fbem_n_vertices(e%gtype)) ! Local coordinates of the element subdivision
    integer                  :: method                             ! Method used when calculating the nearest element point
    integer                  :: gln_near                           ! 1D Gauss-Legendre integration points required by the quasi-singular function
    integer                  :: gln                                ! 1D Gauss-Legendre integration points used in the integration
    integer                  :: ps                                 ! Selected precalculated dataset
    integer                  :: i                                  ! Counter
    ! POINT BODY LOAD
    if (e%d.eq.0) then
      r=e%x(:,1)-x_i
      rmin=sqrt(dot_product(r,r))
      if (rmin.eq.0.d0) then
        call fbem_error_message(output_unit,0,'fbem_bem_staela3d_hbie_bl_auto',0,'it is not possible to collocate at a point load')
      else
        call fbem_bem_staela3d_hbie_d(e%x(:,1),x_i,n_i,nu,lbl(1,:,:))
        return
      end if
    ! LINE, SURFACE OR VOLUME BODY LOAD
    ! Determine if interior or exterior integration
    !   - Interior integration (delta=1) requires: xi_i
    !   - Exterior integration (delta=0) requires: x_i, barxi, rmin and d
    ! Use the element ball
    else
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
    end if
    ! Integrate
    select case (delta)
      case (1)
        stop 'not implemented yet fbem_bem_staela3d_hbie_bl_int'
        !call fbem_bem_staela3d_hbie_bl_int(e,barxi,n_i,mu,nu,lbl)
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
          call fbem_bem_staela3d_hbie_bl_ext_pre(ps,e,x_i,n_i,mu,nu,lbl)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_staela3d_hbie_bl_ext_adp(e,xi_s,x_i,n_i,mu,nu,qsp,1,ns,lbl)
        end if
    end select
  end subroutine fbem_bem_staela3d_hbie_bl_auto

  !! This subroutine calculates the kernels for HBIE exterior integration, needing precalculated data at integration points.
  subroutine fbem_bem_staela3d_hbie_bl_ext_pre(ps,e,x_i,n_i,mu,nu,lbl)
    implicit none
    ! I/O
    integer                :: ps                  !! Selected precalculated dataset
    type(fbem_bem_element) :: e                   !! Element
    real(kind=real64)      :: x_i(3)              !! Collocation point position vector
    real(kind=real64)      :: n_i(3)              !! Unit normal at the collocation point
    real(kind=real64)      :: mu                  !! Shear modulus
    real(kind=real64)      :: nu                  !! Poisson's ratio
    real(kind=real64)      :: lbl(e%n_snodes,3,3) !! l kernels
    ! Local
    integer                :: il, ik              ! Counter for load / observation components
    integer                :: kip                 ! Counter variable for integration points loop
    real(kind=real64)      :: x(3)                ! Position vector at integration point
    real(kind=real64)      :: sphijw(e%n_snodes)  ! phi^s * jacobian * weight at integration point
    real(kind=real64)      :: rv(3)               ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)      :: r, d1r, d1r2        ! Distance vector module and its inverse
    real(kind=real64)      :: drdx(3)             ! Distance vector derivatives with respect to x_k
    real(kind=real64)      :: drdni               ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)      :: cted1, cted2        ! Auxiliary constants
    real(kind=real64)      :: fs_d                ! Fundamental solutions values
    ! Initialize kernels
    lbl=0.d0
    ! Initialize auxiliary constants for fundamental solutions calculation
    cted1=-1.d0/(8.d0*c_pi*(1.d0-nu))
    cted2=1.d0-2.d0*nu
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      x=e%ps_x(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      d1r =1.d0/r
      d1r2=d1r**2
      drdx=rv*d1r
      drdni=-dot_product(drdx,n_i)
      do il=1,3
        do ik=1,3
          fs_d=d1r2*((cted2*c_dkr(il,ik)+3.d0*drdx(il)*drdx(ik))*drdni+cted2*(n_i(il)*drdx(ik)-n_i(ik)*drdx(il)))
          lbl(:,il,ik)=lbl(:,il,ik)+fs_d*sphijw
        end do
      end do
    end do
    ! Multiply by constants
    lbl=cted1*lbl
  end subroutine fbem_bem_staela3d_hbie_bl_ext_pre

  !! This subroutine calculates adaptatively the body load kernels for HBIE exterior integration using Telles transformation and
  !! subdivision if needed.
  recursive subroutine fbem_bem_staela3d_hbie_bl_ext_adp(e,xi_s,x_i,n_i,mu,nu,qsp,ks,ns,lbl)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                      !! Element
    real(kind=real64)        :: xi_s(e%d,fbem_n_vertices(e%gtype))     !! Subdivision of the parent element
    real(kind=real64)        :: x_i(3)                                 !! Collocation point position vector
    real(kind=real64)        :: n_i(3)                                 !! Collocation point unit normal vector
    real(kind=real64)        :: mu                                     !! Shear modulus
    real(kind=real64)        :: nu                                     !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp                                    !! Quasi-singular integration parameters
    integer                  :: ks                                     !! Current level of subdivisions
    integer                  :: ns                                     !! Maximum level of subdivision
    real(kind=real64)        :: lbl(e%n_snodes,3,3)                    !! lbl kernels
    ! Local
    integer                  :: gln_near                               ! 1D Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer                  :: gln                                    ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical                  :: subdivide                              ! True if subdivision has to be performed
    real(kind=real64)        :: barxi(e%d)                             ! Nearest element coordinate with respect to collocation point
    real(kind=real64)        :: barxip(e%d)                            ! Nearest element subdivision local coordinate with respect to collocation point
    real(kind=real64)        :: rmin                                   ! Minimum distance between collocation point and barxi on the element
    real(kind=real64)        :: barr                                   ! Telles jacobian at barxi
    real(kind=real64)        :: cl                                     ! Characteristic length
    real(kind=real64)        :: d                                      ! Normalized distance between collocation point and element subdivision
    integer                  :: method                                 ! Method used in nearest point algorithm
    real(kind=real64)        :: tmp_xi_s(e%d,fbem_n_vertices(e%gtype)) ! Subdivision
    real(kind=real64)        :: x_s(3,e%n_gnodes)                      ! Coordinates of the element subdivision
    real(kind=real64)        :: lbl_tmp(e%n_snodes,3,3)                ! lbl kernels (temporary)
    ! Initialize
    if (ks.eq.1) then
      lbl=0.d0
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
          stop 'not yet : fbem_bem_staela3d_hbie_bl_ext_adp'
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
        call fbem_warning_message(error_unit,0,'fbem_bem_staela3d_hbie_bl_ext_adp',ns,'maximum number of subdivisions reached')
        gln_near=30
      end if
    else
      if (gln_near.eq.0) subdivide=.true.
    end if
    ! Subdivide
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
          call fbem_bem_staela3d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,lbl)
          ! SUBLINE 2
          tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
          tmp_xi_s(1,2)=xi_s(1,2)
          call fbem_bem_staela3d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,lbl)
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
              call fbem_bem_staela3d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,lbl)
              ! SUBTRI 2
              tmp_xi_s(:,1)=xi_s(:,2)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              call fbem_bem_staela3d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,lbl)
              ! SUBTRI 3
              tmp_xi_s(:,1)=xi_s(:,3)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              call fbem_bem_staela3d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,lbl)
              ! SUBTRI 4
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
              call fbem_bem_staela3d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,lbl)
            ! QUADRILATERALS
            case (4)
              ! SUBQUAD 1
              tmp_xi_s(:,1)=xi_s(:,1)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,3)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,4)=0.50d0*(xi_s(:,1)+xi_s(:,4))
              call fbem_bem_staela3d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,lbl)
              ! SUBQUAD 2
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,2)=xi_s(:,2)
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,4)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              call fbem_bem_staela3d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,lbl)
              ! SUBQUAD 3
              tmp_xi_s(:,1)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=xi_s(:,3)
              tmp_xi_s(:,4)=0.50d0*(xi_s(:,3)+xi_s(:,4))
              call fbem_bem_staela3d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,lbl)
              ! SUBQUAD 4
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,4))
              tmp_xi_s(:,2)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,4)=xi_s(:,4)
              call fbem_bem_staela3d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,lbl)
          end select
        !
        ! VOLUME LOAD
        !
        case (3)
          stop 'not yet : fbem_bem_staela3d_hbie_bl_ext_adp'
      end select
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_staela3d_hbie_bl_ext_st(e,xi_s,x_i,n_i,barxip,barr,mu,nu,gln,lbl_tmp)
      lbl=lbl+lbl_tmp
    end if
  end subroutine fbem_bem_staela3d_hbie_bl_ext_adp

  !! This subroutine calculates the body load kernels for HBIE exterior integration (near collocation points) using Telles
  !! transformation within a subdivision of the element, needing only basic data.
  subroutine fbem_bem_staela3d_hbie_bl_ext_st(e,xi_s,x_i,n_i,barxip,barr,mu,nu,gln,lbl)
    implicit none
    ! I/O
    type(fbem_bem_element)       :: e                                  !! Integration element
    real(kind=real64)            :: xi_s(e%d,fbem_n_vertices(e%gtype)) !! Coordinates of the subdivision of the parent element (xi space [-1,1])
    real(kind=real64)            :: x_i(3)                             !! Collocation point position vector
    real(kind=real64)            :: n_i(3)                             !! Collocation point unit normal vector
    real(kind=real64)            :: barxip(e%d)                        !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)            :: barr                               !! Telles jacobian at barxip
    real(kind=real64)            :: mu                                 !! Shear modulus
    real(kind=real64)            :: nu                                 !! Poisson's ratio
    integer                      :: gln                                !! 1D Gauss-Legendre number of integration points (<=32)
    real(kind=real64)            :: lbl(e%n_snodes,3,3)                !! lbl kernels
    ! Local
    integer                      :: il, ik                             ! Counter variables for load direction and observation direction
    integer                      :: kphi                               ! Counter variable for shape functions loops
    integer                      :: kip                                ! Counter variable of integration points
    integer                      :: k1                                 ! Counter variable for reference coordinate xi_1
    integer                      :: k2                                 ! Counter variable for reference coordinate xi_2
    real(kind=real64)            :: gamma(e%d)                         ! Coordinate gamma (Telles transformation space [-1,1])
    real(kind=real64)            :: w(e%d)                             ! Weights of an integration point
    type(fbem_telles_parameters) :: telles_parameters(e%d)             ! Telles parameters
    real(kind=real64)            :: jt(e%d)                            ! Telles jacobian: xip->gamma
    real(kind=real64)            :: xip(e%d)                           ! Coordinate xip (subdivision space [-1,1])
    real(kind=real64)            :: js                                 ! Subdivision jacobian: xi->xip
    real(kind=real64)            :: xin(e%d)                           ! Coordinate xi [xi_s(1,1),xi_s(1,2)]
    real(kind=real64)            :: aux(10)                            ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)                   ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi(e%n_gnodes)               ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: dgphidxi1(e%n_gnodes)              ! Geometrical shape functions derivatives values
    real(kind=real64)            :: dgphidxi2(e%n_gnodes)              ! Geometrical shape functions derivatives values
    real(kind=real64)            :: dxidxi1p(2)                        ! xi derivatives with respect to xi1p
    real(kind=real64)            :: dxidxi2p(2)                        ! xi derivatives with respect to xi2p
    real(kind=real64)            :: xipp(2)                            ! Coordinate xipp used for quadrilateral-triangle transformation
    real(kind=real64)            :: barxipp(2)                         ! Coordinate xipp of collocation point
    real(kind=real64)            :: jqt                                ! Jacobian of the quadrilateral-triangle transformation
    real(kind=real64)            :: sphi(e%n_snodes)                   ! Functional shape functions values
    real(kind=real64)            :: x(3)                               ! Position vector at xi
    real(kind=real64)            :: N(3), T(3), T1(3), T2(3)           ! Tangent vector at xi
    real(kind=real64)            :: rv(3)                              ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)            :: r, d1r, d1r2                       ! Distance vector module, its inverse and log(1/r)
    real(kind=real64)            :: drdx(3)                            ! Distance vector derivatives with respect to x_k
    real(kind=real64)            :: drdni                              ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)            :: jg                                 ! Geometric jacobian
    real(kind=real64)            :: jw                                 ! Jacobians * weight
    real(kind=real64)            :: sphijw(e%n_snodes)                 ! secondary shape functions * jw
    real(kind=real64)            :: cted1, cted2                       ! Auxiliary constants
    real(kind=real64)            :: fs_d                               ! Fundamental solutions values
    ! Initialize kernels
    lbl=0.d0
    ! Initialize auxiliary constants for fundamental solutions calculation
    cted1=-1.d0/(8.d0*c_pi*(1.d0-nu))
    cted2=1.d0-2.d0*nu
    ! Calculate
    select case (e%d)
      !
      ! LINE LOAD
      !
      case (1)
        ! Calculate Telles parameters
        telles_parameters(1)=fbem_telles11_calculate_parameters(barxip(1),barr)
        ! Numerical integration
        do kip=1,gl11_n(gln)
          ! GAMMA COORDINATE
          gamma(1)=gl11_xi(kip,gln)
          w(1)=gl11_w(kip,gln)
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
          d1r=1.d0/r
          d1r2=d1r**2
          drdx=rv*d1r
          drdni=-dot_product(drdx,n_i)
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
          ! Add integration points
          do il=1,3
            do ik=1,3
              fs_d=d1r2*((cted2*c_dkr(il,ik)+3.d0*drdx(il)*drdx(ik))*drdni+cted2*(n_i(il)*drdx(ik)-n_i(ik)*drdx(il)))
              lbl(:,il,ik)=lbl(:,il,ik)+fs_d*sphijw
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
#               define etype e%gtype
#               define delta 0.d0
#               define xi xin
#               define phi gphi
#               define dphidxi1 dgphidxi1
#               define dphidxi2 dgphidxi2
#               include <phi_and_dphidxik_2d.rc>
#               undef etype
#               undef delta
#               undef xi
#               undef phi
#               undef dphidxi1
#               undef dphidxi2
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
                ! xi->x jacobian
                jg=sqrt(dot_product(N,N))
                ! Distance vector and derived terms
                rv=x-x_i
                r=sqrt(dot_product(rv,rv))
                d1r=1.d0/r
                d1r2=d1r**2
                drdx=rv*d1r
                drdni=-dot_product(drdx,n_i)
                ! Jacobians * weights
                jw=jg*js*jt(1)*jt(2)*w(1)*w(2)
                ! FUNCTIONAL SHAPE FUNCTIONS
                ! Functional shape functions (body load variables) at xi
#               define etype e%stype
#               define delta e%stype_delta
#               define xi xin
#               define phi sphi
#               include <phi_2d.rc>
#               undef etype
#               undef delta
#               undef xi
#               undef phi
                ! Functional shape functions * jacobians * weights
                sphijw=sphi*jw
                ! Add integration points
                do il=1,3
                  do ik=1,3
                    fs_d=d1r2*((cted2*c_dkr(il,ik)+3.d0*drdx(il)*drdx(ik))*drdni+cted2*(n_i(il)*drdx(ik)-n_i(ik)*drdx(il)))
                    lbl(:,il,ik)=lbl(:,il,ik)+fs_d*sphijw
                  end do
                end do
              end do
            end do
          ! TRIANGULAR ELEMENTS
          case (3)
            ! Telles transformation is applied to Gauss-Legendre*Gauss-Legendre quadrature before the quadrilateral->triangle
            ! transformation. Because barxip for triangles are given in area triangle coordinates, for Telles transformation
            ! they must be transformed to quadrilateral coordinates. A special treatment is needed when barxi_2 is near 1, because
            ! transformation diverges
            ! BARXIP->BARXIPP TRANSFORMATION
            if (barxip(2).gt.0.995d0) then
              barxipp(1)=0.5d0
              barxipp(2)=1.0d0
            else
              barxipp(1)=barxip(1)/(1.d0-barxip(2))
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
#               define delta 0.d0
#               define xi xin
#               define phi gphi
#               define dphidxi1 dgphidxi1
#               define dphidxi2 dgphidxi2
#               include <phi_and_dphidxik_2d.rc>
#               undef etype
#               undef delta
#               undef xi
#               undef phi
#               undef dphidxi1
#               undef dphidxi2
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
                ! xi->x jacobian
                jg=sqrt(dot_product(N,N))
                ! Distance vector and derived terms
                rv=x-x_i
                r=sqrt(dot_product(rv,rv))
                d1r=1.d0/r
                d1r2=d1r**2
                drdx=rv*d1r
                drdni=-dot_product(drdx,n_i)
                ! Jacobians * weights
                jw=jg*js*jqt*jt(1)*jt(2)*w(1)*w(2)
                ! FUNCTIONAL SHAPE FUNCTIONS
                ! Functional shape functions (body load variables) at xi
#               define etype e%stype
#               define delta e%stype_delta
#               define xi xin
#               define phi sphi
#               include <phi_2d.rc>
#               undef etype
#               undef delta
#               undef xi
#               undef phi
                ! Functional shape functions * jacobians * weights
                sphijw=sphi*jw
                ! Add integration points
                do il=1,3
                  do ik=1,3
                    fs_d=d1r2*((cted2*c_dkr(il,ik)+3.d0*drdx(il)*drdx(ik))*drdni+cted2*(n_i(il)*drdx(ik)-n_i(ik)*drdx(il)))
                    lbl(:,il,ik)=lbl(:,il,ik)+fs_d*sphijw
                  end do
                end do
              end do
            end do
            case default
              call fbem_error_message(error_unit,0,__FILE__,__LINE__,'n_edges not valid')
        end select
      !
      ! VOLUME LOAD
      !
      case (3)
        stop 'not yet : fbem_bem_staela3d_hbie_bl_ext_st'
      !
      ! OTHERS
      !
      case default
        call fbem_error_message(output_unit,0,'fbem_bem_staela3d_hbie_bl_int',0,&
                                'it is only possible to integrate at line, surface or volume loads')
    end select
    ! Multiply by constants
    lbl=cted1*lbl
  end subroutine fbem_bem_staela3d_hbie_bl_ext_st



  ! --------------------------------------------------------------------------------------------------------------------------------
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE) - TRACTION-FREE HALF-SPACE FUNDAMENTAL SOLUTION
  ! --------------------------------------------------------------------------------------------------------------------------------

  !! Three-dimensional traction-free half-space fundamental solution (elastostatics)
  subroutine fbem_bem_staela3d_hs_sbie(mu,nu,np,xp,x_i,x,n,u,t)
    implicit none
    real(kind=real64), intent( in) :: mu      !! Shear modulus
    real(kind=real64), intent( in) :: nu      !! Poisson's ratio
    integer          , intent( in) :: np      !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64), intent( in) :: xp      !! Coordinate of the plane: x(np)
    real(kind=real64), intent( in) :: x_i(3)  !! Collocation point (original, not the image)
    real(kind=real64), intent( in) :: x(3)
    real(kind=real64), intent( in) :: n(3)
    real(kind=real64), intent(out) :: u(3,3)
    real(kind=real64), intent(out) :: t(3,3)
    ! Local
    integer           :: k1, k2, k3
    real(kind=real64) :: x_i_1(3)
    real(kind=real64) :: x_1(3)
    real(kind=real64) :: n_1(3)
    real(kind=real64) :: L(3,3), LT(3,3), theta, un(3,3),sn(3,3)
    integer           :: i, j, k, o, p, q
    real(kind=real64) :: sigma(3,3,3)
    real(kind=real64) :: xp_i(3)
    real(kind=real64) :: rv(3), r
    real(kind=real64) :: rpv(3)
    real(kind=real64) :: rp, rp2, rp3, rp5, rp7
    real(kind=real64) :: c
    real(kind=real64) :: barx, z
    real(kind=real64) :: c1, c2, c3, c4, c5
    ! Checkings
    if (np.ne.-3) then
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,'only the case of a half-space with unit normal n=(0,0,-1) is allowed')
    end if
    if (x_i(3).lt.xp) then
      write(error_unit,*) 'Error: x_i must belong to the half-space'
      stop
    end if
    if (x(3).lt.xp) then
      write(error_unit,*) 'Error: x must belong to the half-space (1)'
      stop
    end if
    ! Initialization
    c1=1.-nu
    c2=1.-2.*nu
    c3=3.-2.*nu
    c4=3.-4.*nu
    c5=5.-4.*nu
    c=x_i(3)-xp
    z=x(3)-xp
    rv=x-x_i
    r=sqrt(rv(1)**2+rv(2)**2)
    xp_i=x_i
    xp_i(3)=-xp_i(3)
    rpv=x-xp_i
    rp=norm2(rpv)
    rp2=rp**2
    rp3=rp2*rp
    rp5=rp3*rp2
    rp7=rp5*rp2
    ! Normal component
    ! In local cylindrical coordinates
    un=0.
    un(3,1)=r*(c4*(z-c)/rp3-4.*c1*c2/rp/(rp+z+c)+6.*c*z*(z+c)/rp5)
    un(3,3)=(8.*c1**2-c4)/rp+(c4*(z+c)**2-2.*c*z)/rp3+6.*c*z*(z+c)**2/rp5
    sn=0.
    sn(1,1)=-c2*(z+7.*c)/rp3+4.*c1*c2/rp/(rp+z+c)+(6.*c*c2*(z+c)**2-6.*c**2*(z+c)-3.*c4*r**2*(z-c))/rp5-30.*c*r**2*z*(z+c)/rp7
    sn(1,3)=r*(c2/rp3-(3.*c4*z*(z+c)-3.*c*(3.*z+c))/rp5-30.*c*z*(z+c)**2/rp7)
    sn(3,1)=sn(1,3)
    sn(2,2)=c2*(c4*(z+c)-6.*c)/rp3-4.*c1*c2/rp/(rp+z+c)+6.*c2*c*(z+c)**2/rp5-6.*c**2*(z+c)/rp5
    sn(3,3)=c2*(z-c)/rp3-(3.*c4*z*(z+c)**2-3.*c*(z+c)*(5.*z-c))/rp5-30.*c*z*(z+c)**3/rp7
    ! In global coordinates
    ! n=(0,0,-1)
    theta=atan2(x(2)-x_i(2),x(1)-x_i(1))
    L=0.
    L(1,1)=cos(theta)
    L(1,2)=-sin(theta)
    L(2,1)=sin(theta)
    L(2,2)=cos(theta)
    L(3,3)=1.
    LT=transpose(L)
    un=matmul(L,matmul(un,LT))
    sn=matmul(L,matmul(sn,LT))
    ! Tangential components (in global coordinates)
    u=0.
    u(1,1)=1./rp+c4*rv(1)**2/rp3+2.*c*z/rp3*(1.-3.*rv(1)**2/rp2)+4.*c1*c2/(rp+z+c)*(1.-rv(1)**2/rp/(rp+z+c))
    u(1,2)=rv(1)*rv(2)*(c4/rp3-6.*c*z/rp5-4.*c1*c2/rp/(rp+z+c)**2)
    u(1,3)=rv(1)*(c4*(z-c)/rp3-6.*c*z*(z+c)/rp5+4.*c1*c2/rp/(rp+z+c))
    u(2,1)=rv(1)*rv(2)*(c4/rp3-6.*c*z/rp5-4.*c1*c2/rp/(rp+z+c)**2)
    u(2,2)=1./rp+c4*rv(2)**2/rp3+2.*c*z/rp3*(1.-3.*rv(2)**2/rp2)+4.*c1*c2/(rp+z+c)*(1.-rv(2)**2/rp/(rp+z+c))
    u(2,3)=rv(2)*(c4*(z-c)/rp3-6.*c*z*(z+c)/rp5+4.*c1*c2/rp/(rp+z+c))
    sigma=0.
    sigma(1,1,1)=rv(1)*(c2*c5/rp3-3.*c4*rv(1)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(3.-rv(1)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                +6.*c/rp5*(3.*c-c3*(z+c)+5.*rv(1)**2*z/rp2))
    sigma(2,2,1)=rv(1)*(c2*c4/rp3-3.*c4*rv(2)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(1.-rv(2)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                +6.*c/rp5*(c-c2*(z+c)+5.*rv(2)**2*z/rp2))
    sigma(3,3,1)=rv(1)*(-c2/rp3-3.*c4*(z+c)**2/rp5+6.*c/rp5*(c+c2*(z+c)+5.*z*(z+c)**2/rp2))
    sigma(2,3,1)=rv(1)*rv(2)*(-3.*c4*(z+c)/rp5+6.*c/rp5*(c2+5.*z*(z+c)/rp2))
    sigma(3,1,1)=c2*(z-c)/rp3-3.*c4*rv(1)**2*(z+c)/rp5-6.*c/rp5*(z*(z+c)-c2*rv(1)**2-5.*rv(1)**2*z*(z+c)/rp2)
    sigma(1,2,1)=rv(2)*(c2/rp3-3.*c4*rv(1)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(1.-rv(1)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                -6.*c*z/rp5*(1.-5.*rv(1)**2/rp2))
    sigma(3,2,1)=sigma(2,3,1)
    sigma(1,3,1)=sigma(3,1,1)
    sigma(2,1,1)=sigma(1,2,1)

    sigma(1,1,2)=rv(2)*(c2*c4/rp3-3.*c4*rv(1)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(1.-rv(1)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                +6.*c/rp5*(c-c2*(z+c)+5.*rv(1)**2*z/rp2))
    sigma(2,2,2)=rv(2)*(c2*c5/rp3-3.*c4*rv(2)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(3.-rv(2)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                +6.*c/rp5*(3.*c-c3*(z+c)+5.*rv(2)**2*z/rp2))
    sigma(3,3,2)=rv(2)*(-c2/rp3-3.*c4*(z+c)**2/rp5+6.*c/rp5*(c+c2*(z+c)+5.*z*(z+c)**2/rp2))
    sigma(1,3,2)=rv(1)*rv(2)*(-3.*c4*(z+c)/rp5+6.*c/rp5*(c2+5.*z*(z+c)/rp2))
    sigma(3,2,2)=c2*(z-c)/rp3-3.*c4*rv(2)**2*(z+c)/rp5-6.*c/rp5*(z*(z+c)-c2*rv(2)**2-5.*rv(2)**2*z*(z+c)/rp2)
    sigma(2,1,2)=rv(1)*(c2/rp3-3.*c4*rv(2)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(1.-rv(2)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                -6.*c*z/rp5*(1.-5.*rv(2)**2/rp2))
    sigma(3,1,2)=sigma(1,3,2)
    sigma(2,3,2)=sigma(3,2,2)
    sigma(1,2,2)=sigma(2,1,2)
    ! Add normal and tangential components
    sigma(:,:,3)=sn
    u=u+un
    u=u/(16.*c_pi*c1*mu)
    sigma=sigma/(8.*c_pi*c1)
    ! Calculate traction
    do i=1,3
      do j=1,3
        t(i,j)=dot_product(sigma(j,:,i),n)
      end do
    end do
    ! Add Kelvin's fundamental solution
    do i=1,3
      do j=1,3
        u(i,j)=u(i,j)+(1.d0/norm2(rv))*(c4*c_dkr(i,j)+rv(i)/norm2(rv)*rv(j)/norm2(rv))/(16.*c_pi*c1*mu)
        t(i,j)=t(i,j)-(1.d0/norm2(rv)**2)*((c2*c_dkr(i,j)+3.d0*rv(i)/norm2(rv)*rv(j)/norm2(rv))*dot_product(rv/norm2(rv),n)&
                            +c2*(n(i)*rv(j)/norm2(rv)-n(j)*rv(i)/norm2(rv)))/(8.*c_pi*c1)
      end do
    end do
  end subroutine fbem_bem_staela3d_hs_sbie

    !! Three-dimensional traction-free half-space fundamental solution (elastostatics) - with complex modulus
  subroutine fbem_bem_staela3d_hs_sbie_cmplx(mu,nu,np,xp,x_i,x,n,u,t)
    implicit none
    complex(kind=real64), intent( in) :: mu      !! Shear modulus
    complex(kind=real64), intent( in) :: nu      !! Poisson's ratio
    integer             , intent( in) :: np      !! Unit normal of the plane: -1,1,-2,2,-3,3
    real   (kind=real64), intent( in) :: xp      !! Coordinate of the plane: x(np)
    real   (kind=real64), intent( in) :: x_i(3)  !! Collocation point (original, not the image)
    real   (kind=real64), intent( in) :: x(3)
    real   (kind=real64), intent( in) :: n(3)
    complex(kind=real64), intent(out) :: u(3,3)
    complex(kind=real64), intent(out) :: t(3,3)
    ! Local
    integer           :: k1, k2, k3
    real(kind=real64) :: x_i_1(3)
    real(kind=real64) :: x_1(3)
    real(kind=real64) :: n_1(3)
    real(kind=real64) :: L(3,3), LT(3,3), theta
    complex(kind=real64) :: un(3,3),sn(3,3)
    integer           :: i, j, k, o, p, q
    complex(kind=real64) :: sigma(3,3,3)
    real(kind=real64) :: xp_i(3)
    real(kind=real64) :: rv(3), r
    real(kind=real64) :: rpv(3)
    real(kind=real64) :: rp, rp2, rp3, rp5, rp7
    real(kind=real64) :: c
    real(kind=real64) :: barx, z
    complex(kind=real64) :: c1, c2, c3, c4, c5
    ! Checkings
    if (np.ne.-3) then
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,'only the case of a half-space with unit normal n=(0,0,-1) is allowed')
    end if
    if (x_i(3).lt.xp) then
      write(error_unit,*) 'Error: x_i must belong to the half-space'
      stop
    end if
    if (x(3).lt.xp) then
      write(error_unit,*) 'Error: x must belong to the half-space (1)'
      stop
    end if
    ! Initialization
    c1=1.-nu
    c2=1.-2.*nu
    c3=3.-2.*nu
    c4=3.-4.*nu
    c5=5.-4.*nu
    c=x_i(3)-xp
    z=x(3)-xp
    rv=x-x_i
    r=sqrt(rv(1)**2+rv(2)**2)
    xp_i=x_i
    xp_i(3)=-xp_i(3)
    rpv=x-xp_i
    rp=norm2(rpv)
    rp2=rp**2
    rp3=rp2*rp
    rp5=rp3*rp2
    rp7=rp5*rp2
    ! Normal component
    ! In local cylindrical coordinates
    un=0.
    un(3,1)=r*(c4*(z-c)/rp3-4.*c1*c2/rp/(rp+z+c)+6.*c*z*(z+c)/rp5)
    un(3,3)=(8.*c1**2-c4)/rp+(c4*(z+c)**2-2.*c*z)/rp3+6.*c*z*(z+c)**2/rp5
    sn=0.
    sn(1,1)=-c2*(z+7.*c)/rp3+4.*c1*c2/rp/(rp+z+c)+(6.*c*c2*(z+c)**2-6.*c**2*(z+c)-3.*c4*r**2*(z-c))/rp5-30.*c*r**2*z*(z+c)/rp7
    sn(1,3)=r*(c2/rp3-(3.*c4*z*(z+c)-3.*c*(3.*z+c))/rp5-30.*c*z*(z+c)**2/rp7)
    sn(3,1)=sn(1,3)
    sn(2,2)=c2*(c4*(z+c)-6.*c)/rp3-4.*c1*c2/rp/(rp+z+c)+6.*c2*c*(z+c)**2/rp5-6.*c**2*(z+c)/rp5
    sn(3,3)=c2*(z-c)/rp3-(3.*c4*z*(z+c)**2-3.*c*(z+c)*(5.*z-c))/rp5-30.*c*z*(z+c)**3/rp7
    ! In global coordinates
    ! n=(0,0,-1)
    theta=atan2(x(2)-x_i(2),x(1)-x_i(1))
    L=0.
    L(1,1)=cos(theta)
    L(1,2)=-sin(theta)
    L(2,1)=sin(theta)
    L(2,2)=cos(theta)
    L(3,3)=1.
    LT=transpose(L)
    un=matmul(L,matmul(un,LT))
    sn=matmul(L,matmul(sn,LT))
    ! Tangential components (in global coordinates)
    u=0.
    u(1,1)=1./rp+c4*rv(1)**2/rp3+2.*c*z/rp3*(1.-3.*rv(1)**2/rp2)+4.*c1*c2/(rp+z+c)*(1.-rv(1)**2/rp/(rp+z+c))
    u(1,2)=rv(1)*rv(2)*(c4/rp3-6.*c*z/rp5-4.*c1*c2/rp/(rp+z+c)**2)
    u(1,3)=rv(1)*(c4*(z-c)/rp3-6.*c*z*(z+c)/rp5+4.*c1*c2/rp/(rp+z+c))
    u(2,1)=rv(1)*rv(2)*(c4/rp3-6.*c*z/rp5-4.*c1*c2/rp/(rp+z+c)**2)
    u(2,2)=1./rp+c4*rv(2)**2/rp3+2.*c*z/rp3*(1.-3.*rv(2)**2/rp2)+4.*c1*c2/(rp+z+c)*(1.-rv(2)**2/rp/(rp+z+c))
    u(2,3)=rv(2)*(c4*(z-c)/rp3-6.*c*z*(z+c)/rp5+4.*c1*c2/rp/(rp+z+c))
    sigma=0.
    sigma(1,1,1)=rv(1)*(c2*c5/rp3-3.*c4*rv(1)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(3.-rv(1)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                +6.*c/rp5*(3.*c-c3*(z+c)+5.*rv(1)**2*z/rp2))
    sigma(2,2,1)=rv(1)*(c2*c4/rp3-3.*c4*rv(2)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(1.-rv(2)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                +6.*c/rp5*(c-c2*(z+c)+5.*rv(2)**2*z/rp2))
    sigma(3,3,1)=rv(1)*(-c2/rp3-3.*c4*(z+c)**2/rp5+6.*c/rp5*(c+c2*(z+c)+5.*z*(z+c)**2/rp2))
    sigma(2,3,1)=rv(1)*rv(2)*(-3.*c4*(z+c)/rp5+6.*c/rp5*(c2+5.*z*(z+c)/rp2))
    sigma(3,1,1)=c2*(z-c)/rp3-3.*c4*rv(1)**2*(z+c)/rp5-6.*c/rp5*(z*(z+c)-c2*rv(1)**2-5.*rv(1)**2*z*(z+c)/rp2)
    sigma(1,2,1)=rv(2)*(c2/rp3-3.*c4*rv(1)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(1.-rv(1)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                -6.*c*z/rp5*(1.-5.*rv(1)**2/rp2))
    sigma(3,2,1)=sigma(2,3,1)
    sigma(1,3,1)=sigma(3,1,1)
    sigma(2,1,1)=sigma(1,2,1)

    sigma(1,1,2)=rv(2)*(c2*c4/rp3-3.*c4*rv(1)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(1.-rv(1)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                +6.*c/rp5*(c-c2*(z+c)+5.*rv(1)**2*z/rp2))
    sigma(2,2,2)=rv(2)*(c2*c5/rp3-3.*c4*rv(2)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(3.-rv(2)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                +6.*c/rp5*(3.*c-c3*(z+c)+5.*rv(2)**2*z/rp2))
    sigma(3,3,2)=rv(2)*(-c2/rp3-3.*c4*(z+c)**2/rp5+6.*c/rp5*(c+c2*(z+c)+5.*z*(z+c)**2/rp2))
    sigma(1,3,2)=rv(1)*rv(2)*(-3.*c4*(z+c)/rp5+6.*c/rp5*(c2+5.*z*(z+c)/rp2))
    sigma(3,2,2)=c2*(z-c)/rp3-3.*c4*rv(2)**2*(z+c)/rp5-6.*c/rp5*(z*(z+c)-c2*rv(2)**2-5.*rv(2)**2*z*(z+c)/rp2)
    sigma(2,1,2)=rv(1)*(c2/rp3-3.*c4*rv(2)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(1.-rv(2)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                -6.*c*z/rp5*(1.-5.*rv(2)**2/rp2))
    sigma(3,1,2)=sigma(1,3,2)
    sigma(2,3,2)=sigma(3,2,2)
    sigma(1,2,2)=sigma(2,1,2)
    ! Add normal and tangential components
    sigma(:,:,3)=sn
    u=u+un
    u=u/(16.*c_pi*c1*mu)
    sigma=sigma/(8.*c_pi*c1)
    ! Calculate traction
    do i=1,3
      do j=1,3
        t(i,j)=dot_product(sigma(j,:,i),n)
      end do
    end do
    ! Add Kelvin's fundamental solution
    do i=1,3
      do j=1,3
        u(i,j)=u(i,j)+(1.d0/norm2(rv))*(c4*c_dkr(i,j)+rv(i)/norm2(rv)*rv(j)/norm2(rv))/(16.d0*c_pi*c1*mu)
        t(i,j)=t(i,j)-(1.d0/norm2(rv)**2)*((c2*c_dkr(i,j)+3.d0*rv(i)/norm2(rv)*rv(j)/norm2(rv))*dot_product(rv/norm2(rv),n)&
                            +c2*(n(i)*rv(j)/norm2(rv)-n(j)*rv(i)/norm2(rv)))/(8.d0*c_pi*c1)
      end do
    end do
  end subroutine fbem_bem_staela3d_hs_sbie_cmplx

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE) - TRACTION-FREE HALF-SPACE FUNDAMENTAL SOLUTION (COMPLEMENTARY PART)
  ! --------------------------------------------------------------------------------------------------------------------------------

  !! Complementary part of the three-dimensional half-space fundamental solution of elastostatics
  subroutine fbem_bem_staela3d_hsc_sbie(mu,nu,np,xp,x_i,x,n,u,t)
    implicit none
    real(kind=real64), intent(in) :: mu      !! Shear modulus
    real(kind=real64), intent(in) :: nu      !! Poisson's ratio
    integer          , intent(in) :: np      !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64), intent(in) :: xp      !! Coordinate of the plane: x(np)
    real(kind=real64), intent(in) :: x_i(3)  !! Collocation point (original, not the image)
    real(kind=real64), intent(in) :: x(3)
    real(kind=real64), intent(in) :: n(3)
    real(kind=real64), intent(out) :: u(3,3)
    real(kind=real64), intent(out) :: t(3,3)
    ! Local
    integer           :: k1, k2, k3
    real(kind=real64) :: x_i_1(3)
    real(kind=real64) :: x_1(3)
    real(kind=real64) :: n_1(3)
    real(kind=real64) :: L(3,3), LT(3,3), theta, un(3,3),sn(3,3)
    integer           :: i, j, k, o, p, q
    real(kind=real64) :: sigma(3,3,3)
    real(kind=real64) :: xp_i(3)
    real(kind=real64) :: rv(3), r
    real(kind=real64) :: rpv(3)
    real(kind=real64) :: rp, rp2, rp3, rp5, rp7
    real(kind=real64) :: c
    real(kind=real64) :: barx, z
    real(kind=real64) :: c1, c2, c3, c4, c5
    real(kind=real64) :: ub(3,3),tb(3,3)
    ! Checkings
    if (np.ne.-3) then
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,'only the case of a half-space with unit normal n=(0,0,-1) is allowed')
    end if
    if (x_i(3).lt.xp) then
      write(error_unit,*) 'Error: x_i must belong to the half-space'
      stop
    end if
    if (x(3).lt.xp) then
      write(error_unit,*) 'Error: x must belong to the half-space (2)'
      stop
    end if
    ! Initialization
    c1=1.-nu
    c2=1.-2.*nu
    c3=3.-2.*nu
    c4=3.-4.*nu
    c5=5.-4.*nu
    c=x_i(3)-xp
    z=x(3)-xp
    rv=x-x_i
    r=sqrt(rv(1)**2+rv(2)**2)
    xp_i=x_i
    xp_i(3)=-xp_i(3)
    rpv=x-xp_i
    rp=norm2(rpv)
    rp2=rp**2
    rp3=rp2*rp
    rp5=rp3*rp2
    rp7=rp5*rp2
    ! Normal component
    ! In local cylindrical coordinates
    un=0.
    un(3,1)=r*(c4*(z-c)/rp3-4.*c1*c2/rp/(rp+z+c)+6.*c*z*(z+c)/rp5)
    un(3,3)=(8.*c1**2-c4)/rp+(c4*(z+c)**2-2.*c*z)/rp3+6.*c*z*(z+c)**2/rp5
    sn=0.
    sn(1,1)=-c2*(z+7.*c)/rp3+4.*c1*c2/rp/(rp+z+c)+(6.*c*c2*(z+c)**2-6.*c**2*(z+c)-3.*c4*r**2*(z-c))/rp5-30.*c*r**2*z*(z+c)/rp7
    sn(1,3)=r*(c2/rp3-(3.*c4*z*(z+c)-3.*c*(3.*z+c))/rp5-30.*c*z*(z+c)**2/rp7)
    sn(3,1)=sn(1,3)
    sn(2,2)=c2*(c4*(z+c)-6.*c)/rp3-4.*c1*c2/rp/(rp+z+c)+6.*c2*c*(z+c)**2/rp5-6.*c**2*(z+c)/rp5
    sn(3,3)=c2*(z-c)/rp3-(3.*c4*z*(z+c)**2-3.*c*(z+c)*(5.*z-c))/rp5-30.*c*z*(z+c)**3/rp7
    ! In global coordinates
    ! n=(0,0,-1)
    theta=atan2(x(2)-x_i(2),x(1)-x_i(1))
    L=0.
    L(1,1)=cos(theta)
    L(1,2)=-sin(theta)
    L(2,1)=sin(theta)
    L(2,2)=cos(theta)
    L(3,3)=1.
    LT=transpose(L)
    un=matmul(L,matmul(un,LT))
    sn=matmul(L,matmul(sn,LT))
    ! Tangential components (in global coordinates)
    u=0.
    u(1,1)=1./rp+c4*rv(1)**2/rp3+2.*c*z/rp3*(1.-3.*rv(1)**2/rp2)+4.*c1*c2/(rp+z+c)*(1.-rv(1)**2/rp/(rp+z+c))
    u(1,2)=rv(1)*rv(2)*(c4/rp3-6.*c*z/rp5-4.*c1*c2/rp/(rp+z+c)**2)
    u(1,3)=rv(1)*(c4*(z-c)/rp3-6.*c*z*(z+c)/rp5+4.*c1*c2/rp/(rp+z+c))
    u(2,1)=rv(1)*rv(2)*(c4/rp3-6.*c*z/rp5-4.*c1*c2/rp/(rp+z+c)**2)
    u(2,2)=1./rp+c4*rv(2)**2/rp3+2.*c*z/rp3*(1.-3.*rv(2)**2/rp2)+4.*c1*c2/(rp+z+c)*(1.-rv(2)**2/rp/(rp+z+c))
    u(2,3)=rv(2)*(c4*(z-c)/rp3-6.*c*z*(z+c)/rp5+4.*c1*c2/rp/(rp+z+c))
    sigma=0.
    sigma(1,1,1)=rv(1)*(c2*c5/rp3-3.*c4*rv(1)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(3.-rv(1)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                +6.*c/rp5*(3.*c-c3*(z+c)+5.*rv(1)**2*z/rp2))
    sigma(2,2,1)=rv(1)*(c2*c4/rp3-3.*c4*rv(2)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(1.-rv(2)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                +6.*c/rp5*(c-c2*(z+c)+5.*rv(2)**2*z/rp2))
    sigma(3,3,1)=rv(1)*(-c2/rp3-3.*c4*(z+c)**2/rp5+6.*c/rp5*(c+c2*(z+c)+5.*z*(z+c)**2/rp2))
    sigma(2,3,1)=rv(1)*rv(2)*(-3.*c4*(z+c)/rp5+6.*c/rp5*(c2+5.*z*(z+c)/rp2))
    sigma(3,1,1)=c2*(z-c)/rp3-3.*c4*rv(1)**2*(z+c)/rp5-6.*c/rp5*(z*(z+c)-c2*rv(1)**2-5.*rv(1)**2*z*(z+c)/rp2)
    sigma(1,2,1)=rv(2)*(c2/rp3-3.*c4*rv(1)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(1.-rv(1)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                -6.*c*z/rp5*(1.-5.*rv(1)**2/rp2))
    sigma(3,2,1)=sigma(2,3,1)
    sigma(1,3,1)=sigma(3,1,1)
    sigma(2,1,1)=sigma(1,2,1)

    sigma(1,1,2)=rv(2)*(c2*c4/rp3-3.*c4*rv(1)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(1.-rv(1)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                +6.*c/rp5*(c-c2*(z+c)+5.*rv(1)**2*z/rp2))
    sigma(2,2,2)=rv(2)*(c2*c5/rp3-3.*c4*rv(2)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(3.-rv(2)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                +6.*c/rp5*(3.*c-c3*(z+c)+5.*rv(2)**2*z/rp2))
    sigma(3,3,2)=rv(2)*(-c2/rp3-3.*c4*(z+c)**2/rp5+6.*c/rp5*(c+c2*(z+c)+5.*z*(z+c)**2/rp2))
    sigma(1,3,2)=rv(1)*rv(2)*(-3.*c4*(z+c)/rp5+6.*c/rp5*(c2+5.*z*(z+c)/rp2))
    sigma(3,2,2)=c2*(z-c)/rp3-3.*c4*rv(2)**2*(z+c)/rp5-6.*c/rp5*(z*(z+c)-c2*rv(2)**2-5.*rv(2)**2*z*(z+c)/rp2)
    sigma(2,1,2)=rv(1)*(c2/rp3-3.*c4*rv(2)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(1.-rv(2)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                -6.*c*z/rp5*(1.-5.*rv(2)**2/rp2))
    sigma(3,1,2)=sigma(1,3,2)
    sigma(2,3,2)=sigma(3,2,2)
    sigma(1,2,2)=sigma(2,1,2)
    ! Add normal and tangential components
    sigma(:,:,3)=sn
    u=u+un
    u=u/(16.*c_pi*c1*mu)
    sigma=sigma/(8.*c_pi*c1)
    ! Calculate traction
    do i=1,3
      do j=1,3
        t(i,j)=dot_product(sigma(j,:,i),n)
      end do
    end do
  end subroutine fbem_bem_staela3d_hsc_sbie

  !! Complementary part of the three-dimensional half-space fundamental solution of elastostatics
  subroutine fbem_bem_staela3d_hsc_sbie_cmplx(mu,nu,np,xp,x_i,x,n,u,t)
    implicit none
    complex(kind=real64), intent(in) :: mu      !! Shear modulus
    complex(kind=real64), intent(in) :: nu      !! Poisson's ratio
    integer          , intent(in) :: np      !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64), intent(in) :: xp      !! Coordinate of the plane: x(np)
    real(kind=real64), intent(in) :: x_i(3)  !! Collocation point (original, not the image)
    real(kind=real64), intent(in) :: x(3)
    real(kind=real64), intent(in) :: n(3)
    complex(kind=real64), intent(out) :: u(3,3)
    complex(kind=real64), intent(out) :: t(3,3)
    ! Local
    integer           :: k1, k2, k3
    real(kind=real64) :: x_i_1(3)
    real(kind=real64) :: x_1(3)
    real(kind=real64) :: n_1(3)
    real(kind=real64) :: L(3,3), LT(3,3), theta
    complex(kind=real64) :: un(3,3),sn(3,3)
    integer           :: i, j, k, o, p, q
    complex(kind=real64) :: sigma(3,3,3)
    real(kind=real64) :: xp_i(3)
    real(kind=real64) :: rv(3), r
    real(kind=real64) :: rpv(3)
    real(kind=real64) :: rp, rp2, rp3, rp5, rp7
    real(kind=real64) :: c
    real(kind=real64) :: barx, z
    complex(kind=real64) :: c1, c2, c3, c4, c5
    complex(kind=real64) :: ub(3,3),tb(3,3)
    ! Checkings
    if (np.ne.-3) then
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,'only the case of a half-space with unit normal n=(0,0,-1) is allowed')
    end if
    if (x_i(3).lt.xp) then
      write(error_unit,*) 'Error: x_i must belong to the half-space'
      stop
    end if
    if (x(3).lt.xp) then
      write(error_unit,*) 'Error: x must belong to the half-space (2)'
      stop
    end if
    ! Initialization
    c1=1.-nu
    c2=1.-2.*nu
    c3=3.-2.*nu
    c4=3.-4.*nu
    c5=5.-4.*nu
    c=x_i(3)-xp
    z=x(3)-xp
    rv=x-x_i
    r=sqrt(rv(1)**2+rv(2)**2)
    xp_i=x_i
    xp_i(3)=-xp_i(3)
    rpv=x-xp_i
    rp=norm2(rpv)
    rp2=rp**2
    rp3=rp2*rp
    rp5=rp3*rp2
    rp7=rp5*rp2
    ! Normal component
    ! In local cylindrical coordinates
    un=0.
    un(3,1)=r*(c4*(z-c)/rp3-4.*c1*c2/rp/(rp+z+c)+6.*c*z*(z+c)/rp5)
    un(3,3)=(8.*c1**2-c4)/rp+(c4*(z+c)**2-2.*c*z)/rp3+6.*c*z*(z+c)**2/rp5
    sn=0.
    sn(1,1)=-c2*(z+7.*c)/rp3+4.*c1*c2/rp/(rp+z+c)+(6.*c*c2*(z+c)**2-6.*c**2*(z+c)-3.*c4*r**2*(z-c))/rp5-30.*c*r**2*z*(z+c)/rp7
    sn(1,3)=r*(c2/rp3-(3.*c4*z*(z+c)-3.*c*(3.*z+c))/rp5-30.*c*z*(z+c)**2/rp7)
    sn(3,1)=sn(1,3)
    sn(2,2)=c2*(c4*(z+c)-6.*c)/rp3-4.*c1*c2/rp/(rp+z+c)+6.*c2*c*(z+c)**2/rp5-6.*c**2*(z+c)/rp5
    sn(3,3)=c2*(z-c)/rp3-(3.*c4*z*(z+c)**2-3.*c*(z+c)*(5.*z-c))/rp5-30.*c*z*(z+c)**3/rp7
    ! In global coordinates
    ! n=(0,0,-1)
    theta=atan2(x(2)-x_i(2),x(1)-x_i(1))
    L=0.
    L(1,1)=cos(theta)
    L(1,2)=-sin(theta)
    L(2,1)=sin(theta)
    L(2,2)=cos(theta)
    L(3,3)=1.
    LT=transpose(L)
    un=matmul(L,matmul(un,LT))
    sn=matmul(L,matmul(sn,LT))
    ! Tangential components (in global coordinates)
    u=0.
    u(1,1)=1./rp+c4*rv(1)**2/rp3+2.*c*z/rp3*(1.-3.*rv(1)**2/rp2)+4.*c1*c2/(rp+z+c)*(1.-rv(1)**2/rp/(rp+z+c))
    u(1,2)=rv(1)*rv(2)*(c4/rp3-6.*c*z/rp5-4.*c1*c2/rp/(rp+z+c)**2)
    u(1,3)=rv(1)*(c4*(z-c)/rp3-6.*c*z*(z+c)/rp5+4.*c1*c2/rp/(rp+z+c))
    u(2,1)=rv(1)*rv(2)*(c4/rp3-6.*c*z/rp5-4.*c1*c2/rp/(rp+z+c)**2)
    u(2,2)=1./rp+c4*rv(2)**2/rp3+2.*c*z/rp3*(1.-3.*rv(2)**2/rp2)+4.*c1*c2/(rp+z+c)*(1.-rv(2)**2/rp/(rp+z+c))
    u(2,3)=rv(2)*(c4*(z-c)/rp3-6.*c*z*(z+c)/rp5+4.*c1*c2/rp/(rp+z+c))
    sigma=0.
    sigma(1,1,1)=rv(1)*(c2*c5/rp3-3.*c4*rv(1)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(3.-rv(1)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                +6.*c/rp5*(3.*c-c3*(z+c)+5.*rv(1)**2*z/rp2))
    sigma(2,2,1)=rv(1)*(c2*c4/rp3-3.*c4*rv(2)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(1.-rv(2)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                +6.*c/rp5*(c-c2*(z+c)+5.*rv(2)**2*z/rp2))
    sigma(3,3,1)=rv(1)*(-c2/rp3-3.*c4*(z+c)**2/rp5+6.*c/rp5*(c+c2*(z+c)+5.*z*(z+c)**2/rp2))
    sigma(2,3,1)=rv(1)*rv(2)*(-3.*c4*(z+c)/rp5+6.*c/rp5*(c2+5.*z*(z+c)/rp2))
    sigma(3,1,1)=c2*(z-c)/rp3-3.*c4*rv(1)**2*(z+c)/rp5-6.*c/rp5*(z*(z+c)-c2*rv(1)**2-5.*rv(1)**2*z*(z+c)/rp2)
    sigma(1,2,1)=rv(2)*(c2/rp3-3.*c4*rv(1)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(1.-rv(1)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                -6.*c*z/rp5*(1.-5.*rv(1)**2/rp2))
    sigma(3,2,1)=sigma(2,3,1)
    sigma(1,3,1)=sigma(3,1,1)
    sigma(2,1,1)=sigma(1,2,1)

    sigma(1,1,2)=rv(2)*(c2*c4/rp3-3.*c4*rv(1)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(1.-rv(1)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                +6.*c/rp5*(c-c2*(z+c)+5.*rv(1)**2*z/rp2))
    sigma(2,2,2)=rv(2)*(c2*c5/rp3-3.*c4*rv(2)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(3.-rv(2)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                +6.*c/rp5*(3.*c-c3*(z+c)+5.*rv(2)**2*z/rp2))
    sigma(3,3,2)=rv(2)*(-c2/rp3-3.*c4*(z+c)**2/rp5+6.*c/rp5*(c+c2*(z+c)+5.*z*(z+c)**2/rp2))
    sigma(1,3,2)=rv(1)*rv(2)*(-3.*c4*(z+c)/rp5+6.*c/rp5*(c2+5.*z*(z+c)/rp2))
    sigma(3,2,2)=c2*(z-c)/rp3-3.*c4*rv(2)**2*(z+c)/rp5-6.*c/rp5*(z*(z+c)-c2*rv(2)**2-5.*rv(2)**2*z*(z+c)/rp2)
    sigma(2,1,2)=rv(1)*(c2/rp3-3.*c4*rv(2)**2/rp5-4.*c1*c2/rp/(rp+z+c)**2*(1.-rv(2)**2*(3.*rp+z+c)/rp2/(rp+z+c))&
                -6.*c*z/rp5*(1.-5.*rv(2)**2/rp2))
    sigma(3,1,2)=sigma(1,3,2)
    sigma(2,3,2)=sigma(3,2,2)
    sigma(1,2,2)=sigma(2,1,2)
    ! Add normal and tangential components
    sigma(:,:,3)=sn
    u=u+un
    u=u/(16.*c_pi*c1*mu)
    sigma=sigma/(8.*c_pi*c1)
    ! Calculate traction
    do i=1,3
      do j=1,3
        t(i,j)=dot_product(sigma(j,:,i),n)
      end do
    end do
  end subroutine fbem_bem_staela3d_hsc_sbie_cmplx

  subroutine fbem_bem_staela3d_hsc_sbie_ext_pre(ps,e,reverse,x_i,mu,nu,np,xp,h,g)
    implicit none
    ! I/O
    integer                :: ps                !! Selected precalculated dataset
    type(fbem_bem_element) :: e                 !! Element
    logical                :: reverse           !! Reverse normal vector
    real(kind=real64)      :: x_i(3)            !! Collocation point position vector
    real(kind=real64)      :: mu                !! Shear modulus
    real(kind=real64)      :: nu                !! Poisson's ratio
    integer                :: np                !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64)      :: xp                !! Coordinate of the plane: x(np)
    real(kind=real64)      :: h(e%n_pnodes,3,3) !! h integration kernels matrix
    real(kind=real64)      :: g(e%n_snodes,3,3) !! g integration kernels matrix
    ! Local
    integer           :: il, ik               ! Counter for load / observation components
    integer           :: kip                  ! Counter variable for integration points loop
    real(kind=real64) :: x(3)                 ! Position vector at integration point
    real(kind=real64) :: n(3)                 ! Unit normal vector at integration point
    real(kind=real64) :: pphijw(e%n_pnodes)   ! phi^p * jacobian * weight at integration point
    real(kind=real64) :: sphijw(e%n_snodes)   ! phi^s * jacobian * weight at integration point
    real(kind=real64) :: fs_u(3,3), fs_t(3,3) ! Fundamental solution
    ! Initialize
    h=0.
    g=0.
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      x=e%ps_x(:,kip,ps)
      n=e%ps_n(:,kip,ps)
      pphijw=e%ps_pphijw(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      call fbem_bem_staela3d_hsc_sbie(mu,nu,np,xp,x_i,x,n,fs_u,fs_t)
      do il=1,3
        do ik=1,3
          h(:,il,ik)=h(:,il,ik)+fs_t(il,ik)*pphijw
          g(:,il,ik)=g(:,il,ik)+fs_u(il,ik)*sphijw
        end do
      end do
    end do
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_staela3d_hsc_sbie_ext_pre

  subroutine fbem_bem_staela3d_hsc_sbie_ext_pre_cmplx(ps,e,reverse,x_i,mu,nu,np,xp,h,g)
    implicit none
    ! I/O
    integer                :: ps                !! Selected precalculated dataset
    type(fbem_bem_element) :: e                 !! Element
    logical                :: reverse           !! Reverse normal vector
    real(kind=real64)      :: x_i(3)            !! Collocation point position vector
    complex(kind=real64)      :: mu                !! Shear modulus
    complex(kind=real64)      :: nu                !! Poisson's ratio
    integer                :: np                !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64)      :: xp                !! Coordinate of the plane: x(np)
    complex(kind=real64)      :: h(e%n_pnodes,3,3) !! h integration kernels matrix
    complex(kind=real64)      :: g(e%n_snodes,3,3) !! g integration kernels matrix
    ! Local
    integer           :: il, ik               ! Counter for load / observation components
    integer           :: kip                  ! Counter variable for integration points loop
    real(kind=real64) :: x(3)                 ! Position vector at integration point
    real(kind=real64) :: n(3)                 ! Unit normal vector at integration point
    real(kind=real64) :: pphijw(e%n_pnodes)   ! phi^p * jacobian * weight at integration point
    real(kind=real64) :: sphijw(e%n_snodes)   ! phi^s * jacobian * weight at integration point
    complex(kind=real64) :: fs_u(3,3), fs_t(3,3) ! Fundamental solution
    ! Initialize
    h=0.
    g=0.
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      x=e%ps_x(:,kip,ps)
      n=e%ps_n(:,kip,ps)
      pphijw=e%ps_pphijw(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      call fbem_bem_staela3d_hsc_sbie_cmplx(mu,nu,np,xp,x_i,x,n,fs_u,fs_t)
      do il=1,3
        do ik=1,3
          h(:,il,ik)=h(:,il,ik)+fs_t(il,ik)*pphijw
          g(:,il,ik)=g(:,il,ik)+fs_u(il,ik)*sphijw
        end do
      end do
    end do
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_staela3d_hsc_sbie_ext_pre_cmplx

  !! This subroutine calculates the kernels for SBIE exterior integration (near collocation points) using Telles transformation
  !! within a subdivision of the element, needing only needs basic data.
  subroutine fbem_bem_staela3d_hsc_sbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,mu,nu,np,xp,gln,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)       :: e                                !! Integration element
    logical                      :: reverse                          !! Reverse normal vector
    real(kind=real64)            :: xi_s(2,fbem_n_vertices(e%gtype)) !! Coordinates of the subdivision of the parent element
    real(kind=real64)            :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)            :: barxip(2)                        !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)            :: barr                             !! Telles jacobian at barxip
    real(kind=real64)            :: mu                               !! Shear modulus
    real(kind=real64)            :: nu                               !! Poisson's ratio
    integer                      :: np                               !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64)            :: xp                               !! Coordinate of the plane: x(np)
    integer                      :: gln                              !! 1D Gauss-Legendre number of integration points (<=32)
    real(kind=real64)            :: h(e%n_pnodes,3,3)                !! h kernel vector
    real(kind=real64)            :: g(e%n_snodes,3,3)                !! g kernel vector
    ! Local
    integer                      :: kphi                       ! Counter variable for shape functions loops
    integer                      :: k1                         ! Counter variable for reference coordinate xi_1
    integer                      :: k2                         ! Counter variable for reference coordinate xi_2
    real(kind=real64)            :: aux(10)                    ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)           ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi1(e%n_gnodes)      ! Geometrical shape functions derivatives values
    real(kind=real64)            :: dgphidxi2(e%n_gnodes)      ! Geometrical shape functions derivatives values
    real(kind=real64)            :: pphi(e%n_pnodes)           ! Functional shape functions values
    real(kind=real64)            :: sphi(e%n_snodes)           ! Functional shape functions values
    real(kind=real64)            :: gamma(2)                   ! Vector of gamma_1,gamma_2 coordinates
    real(kind=real64)            :: w(2)                       ! Weights of the integration rule
    real(kind=real64)            :: xip(2)                     ! Vector of xip_1,xip_2 coordinates
    real(kind=real64)            :: dxidxi1p(2)                ! xi derivatives with respect to xi1p
    real(kind=real64)            :: dxidxi2p(2)                ! xi derivatives with respect to xi2p
    real(kind=real64)            :: js                         ! Subdivision jacobian
    real(kind=real64)            :: xi(2)                      ! Vector of xi_1,xi_2 coordinates
    real(kind=real64)            :: xipp(2)                    ! Coordinate xipp used for quadrilateral-triangle transformation
    real(kind=real64)            :: barxipp(2)                 ! Coordinate xipp of collocation point
    real(kind=real64)            :: jqt                        ! Jacobian of the quadrilateral-triangle transformation
    type(fbem_telles_parameters) :: telles_parameters(2)       ! Telles parameters for each coordinate
    real(kind=real64)            :: jt(2)                      ! Telles jacobian for each coordinate: xi_1->gamma_1 and xi_2->gamma_2
    real(kind=real64)            :: x(3)                       ! Position vector at xi_1,xi_2
    real(kind=real64)            :: T1(3), T2(3)               ! Tangent vectors at xi_1,xi_2
    real(kind=real64)            :: N(3)                       ! Normal vector at xi_1,xi_2
    real(kind=real64)            :: jg                         ! Geometric jacobian
    real(kind=real64)            :: jw                         ! Jacobians * weights
    real(kind=real64)            :: pphijw(e%n_pnodes)         ! Functional shape functions * jw
    real(kind=real64)            :: sphijw(e%n_snodes)         ! Functional shape functions * jw
    integer                      :: il, ik                     ! Counter for load / observation components
    real(kind=real64)            :: fs_u(3,3), fs_t(3,3)       ! Fundamental solutions values
    ! Initialize
    h=0.
    g=0.
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
            ! Fundamental solution
            call fbem_bem_staela3d_hsc_sbie(mu,nu,np,xp,x_i,x,n,fs_u,fs_t)
            ! Add integration points
            do il=1,3
              do ik=1,3
                h(:,il,ik)=h(:,il,ik)+fs_t(il,ik)*pphijw
                g(:,il,ik)=g(:,il,ik)+fs_u(il,ik)*sphijw
              end do
            end do
          end do
        end do
      ! TRIANGULAR ELEMENTS
      case (3)
        ! Telles transformation is applied to Gauss-Legendre*Gauss-Legendre quadrature before the quadrilateral->triangle
        ! transformation. Because barxip for triangles are given in area triangle coordinates, for Telles transformation
        ! they must be transformed to quadrilateral coordinates. A special treatment is needed when barxi_2 is near 1, because
        ! transformation diverges.
        ! BARXIP->BARXIPP TRANSFORMATION
        if (barxip(2).gt.0.995d0) then
          barxipp(1)=0.5d0
          barxipp(2)=1.0d0
        else
          barxipp(1)=barxip(1)/(1.d0-barxip(2))
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
            ! Fundamental solution
            call fbem_bem_staela3d_hsc_sbie(mu,nu,np,xp,x_i,x,n,fs_u,fs_t)
            ! Add integration points
            do il=1,3
              do ik=1,3
                h(:,il,ik)=h(:,il,ik)+fs_t(il,ik)*pphijw
                g(:,il,ik)=g(:,il,ik)+fs_u(il,ik)*sphijw
              end do
            end do
          end do
        end do
        case default
          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'n_edges not valid')
    end select
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_staela3d_hsc_sbie_ext_st

  !! This subroutine calculates the kernels for SBIE exterior integration (near collocation points) using Telles transformation
  !! within a subdivision of the element, needing only needs basic data.
  subroutine fbem_bem_staela3d_hsc_sbie_ext_st_cmplx(e,reverse,xi_s,x_i,barxip,barr,mu,nu,np,xp,gln,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)       :: e                                !! Integration element
    logical                      :: reverse                          !! Reverse normal vector
    real(kind=real64)            :: xi_s(2,fbem_n_vertices(e%gtype)) !! Coordinates of the subdivision of the parent element
    real(kind=real64)            :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)            :: barxip(2)                        !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)            :: barr                             !! Telles jacobian at barxip
    complex(kind=real64)            :: mu                               !! Shear modulus
    complex(kind=real64)            :: nu                               !! Poisson's ratio
    integer                      :: np                               !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64)            :: xp                               !! Coordinate of the plane: x(np)
    integer                      :: gln                              !! 1D Gauss-Legendre number of integration points (<=32)
    complex(kind=real64)            :: h(e%n_pnodes,3,3)                !! h kernel vector
    complex(kind=real64)            :: g(e%n_snodes,3,3)                !! g kernel vector
    ! Local
    integer                      :: kphi                       ! Counter variable for shape functions loops
    integer                      :: k1                         ! Counter variable for reference coordinate xi_1
    integer                      :: k2                         ! Counter variable for reference coordinate xi_2
    real(kind=real64)            :: aux(10)                    ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)           ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi1(e%n_gnodes)      ! Geometrical shape functions derivatives values
    real(kind=real64)            :: dgphidxi2(e%n_gnodes)      ! Geometrical shape functions derivatives values
    real(kind=real64)            :: pphi(e%n_pnodes)           ! Functional shape functions values
    real(kind=real64)            :: sphi(e%n_snodes)           ! Functional shape functions values
    real(kind=real64)            :: gamma(2)                   ! Vector of gamma_1,gamma_2 coordinates
    real(kind=real64)            :: w(2)                       ! Weights of the integration rule
    real(kind=real64)            :: xip(2)                     ! Vector of xip_1,xip_2 coordinates
    real(kind=real64)            :: dxidxi1p(2)                ! xi derivatives with respect to xi1p
    real(kind=real64)            :: dxidxi2p(2)                ! xi derivatives with respect to xi2p
    real(kind=real64)            :: js                         ! Subdivision jacobian
    real(kind=real64)            :: xi(2)                      ! Vector of xi_1,xi_2 coordinates
    real(kind=real64)            :: xipp(2)                    ! Coordinate xipp used for quadrilateral-triangle transformation
    real(kind=real64)            :: barxipp(2)                 ! Coordinate xipp of collocation point
    real(kind=real64)            :: jqt                        ! Jacobian of the quadrilateral-triangle transformation
    type(fbem_telles_parameters) :: telles_parameters(2)       ! Telles parameters for each coordinate
    real(kind=real64)            :: jt(2)                      ! Telles jacobian for each coordinate: xi_1->gamma_1 and xi_2->gamma_2
    real(kind=real64)            :: x(3)                       ! Position vector at xi_1,xi_2
    real(kind=real64)            :: T1(3), T2(3)               ! Tangent vectors at xi_1,xi_2
    real(kind=real64)            :: N(3)                       ! Normal vector at xi_1,xi_2
    real(kind=real64)            :: jg                         ! Geometric jacobian
    real(kind=real64)            :: jw                         ! Jacobians * weights
    real(kind=real64)            :: pphijw(e%n_pnodes)         ! Functional shape functions * jw
    real(kind=real64)            :: sphijw(e%n_snodes)         ! Functional shape functions * jw
    integer                      :: il, ik                     ! Counter for load / observation components
    complex(kind=real64)            :: fs_u(3,3), fs_t(3,3)       ! Fundamental solutions values
    ! Initialize
    h=0.
    g=0.
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
            ! Fundamental solution
            call fbem_bem_staela3d_hsc_sbie_cmplx(mu,nu,np,xp,x_i,x,n,fs_u,fs_t)
            ! Add integration points
            do il=1,3
              do ik=1,3
                h(:,il,ik)=h(:,il,ik)+fs_t(il,ik)*pphijw
                g(:,il,ik)=g(:,il,ik)+fs_u(il,ik)*sphijw
              end do
            end do
          end do
        end do
      ! TRIANGULAR ELEMENTS
      case (3)
        ! Telles transformation is applied to Gauss-Legendre*Gauss-Legendre quadrature before the quadrilateral->triangle
        ! transformation. Because barxip for triangles are given in area triangle coordinates, for Telles transformation
        ! they must be transformed to quadrilateral coordinates. A special treatment is needed when barxi_2 is near 1, because
        ! transformation diverges.
        ! BARXIP->BARXIPP TRANSFORMATION
        if (barxip(2).gt.0.995d0) then
          barxipp(1)=0.5d0
          barxipp(2)=1.0d0
        else
          barxipp(1)=barxip(1)/(1.d0-barxip(2))
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
            ! Fundamental solution
            call fbem_bem_staela3d_hsc_sbie_cmplx(mu,nu,np,xp,x_i,x,n,fs_u,fs_t)
            ! Add integration points
            do il=1,3
              do ik=1,3
                h(:,il,ik)=h(:,il,ik)+fs_t(il,ik)*pphijw
                g(:,il,ik)=g(:,il,ik)+fs_u(il,ik)*sphijw
              end do
            end do
          end do
        end do
        case default
          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'n_edges not valid')
    end select
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_staela3d_hsc_sbie_ext_st_cmplx

  recursive subroutine fbem_bem_staela3d_hsc_sbie_ext_adp(e,reverse,xi_s,x_i,mu,nu,np,xp,qsp,ks,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                !! Element
    logical                  :: reverse                          !! Reverse orientation
    real(kind=real64)        :: xi_s(2,fbem_n_vertices(e%gtype)) !! Subdivision of the parent element
    real(kind=real64)        :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)        :: mu                               !! Shear modulus
    real(kind=real64)        :: nu                               !! Poisson's ratio
    integer                  :: np                               !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64)        :: xp                               !! Coordinate of the plane: x(np)
    type(fbem_qs_parameters) :: qsp                              !! Quasi-singular integration parameters
    integer                  :: ks                               !! Current level of subdivisions
    integer                  :: ns                               !! Maximum level of subdivision
    real(kind=real64)        :: h(e%n_pnodes,3,3)                !! h integration kernels matrix
    real(kind=real64)        :: g(e%n_snodes,3,3)                !! g integration kernels matrix
    ! Local
    real(kind=real64) :: xp_i(3)                              ! Image collocation point position vector
    integer           :: gln_near                             ! 1D Gauss-Legendre integ. points required to integrate only  the quasi-singular integrand
    integer           :: gln                                  ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
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
    real(kind=real64) :: h_tmp(e%n_pnodes,3,3)                ! h integration kernels matrix (temporary)
    real(kind=real64) :: g_tmp(e%n_snodes,3,3)                ! g integration kernels matrix (temporary)
    ! Collocation point image in global axes
    xp_i=x_i
    xp_i(abs(np))=2.*xp-xp_i(abs(np))
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
      call fbem_nearest_element_point_bem(3,e%gtype,e%x,e%cl,xp_i,barxi,rmin,d,method)
      barxip=barxi
    else
      call fbem_obtain_element_subdivision_coordinates(3,e%gtype,e%x,xi_s,x_s)
      cl=fbem_characteristic_length(3,e%gtype,x_s,1.d-12)
      call fbem_nearest_element_point_bem(3,e%gtype,x_s,cl,xp_i,barxip,rmin,d,method)
    end if
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,5,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_staela3d_hsc_sbie_ext_adp',ns,'maximum number of subdivisions reached')
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
          call fbem_bem_staela3d_hsc_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,h,g)
          ! SUBTRI 2
          tmp_xi_s(:,1)=xi_s(:,2)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          call fbem_bem_staela3d_hsc_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,h,g)
          ! SUBTRI 3
          tmp_xi_s(:,1)=xi_s(:,3)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          call fbem_bem_staela3d_hsc_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,h,g)
          ! SUBTRI 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          call fbem_bem_staela3d_hsc_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,h,g)
        ! QUADRILATERALS
        case (4)
          ! SUBQUAD 1
          tmp_xi_s(:,1)=xi_s(:,1)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,3)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          call fbem_bem_staela3d_hsc_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,h,g)
          ! SUBQUAD 2
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=xi_s(:,2)
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,4)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          call fbem_bem_staela3d_hsc_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,h,g)
          ! SUBQUAD 3
          tmp_xi_s(:,1)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=xi_s(:,3)
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          call fbem_bem_staela3d_hsc_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,h,g)
          ! SUBQUAD 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          tmp_xi_s(:,2)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=xi_s(:,4)
          call fbem_bem_staela3d_hsc_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,h,g)
      end select
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_staela3d_hsc_sbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,mu,nu,np,xp,gln,h_tmp,g_tmp)
      h=h+h_tmp
      g=g+g_tmp
    end if
  end subroutine fbem_bem_staela3d_hsc_sbie_ext_adp

  recursive subroutine fbem_bem_staela3d_hsc_sbie_ext_adp_cmplx(e,reverse,xi_s,x_i,mu,nu,np,xp,qsp,ks,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                !! Element
    logical                  :: reverse                          !! Reverse orientation
    real(kind=real64)        :: xi_s(2,fbem_n_vertices(e%gtype)) !! Subdivision of the parent element
    real(kind=real64)        :: x_i(3)                           !! Collocation point position vector
    complex(kind=real64)        :: mu                               !! Shear modulus
    complex(kind=real64)        :: nu                               !! Poisson's ratio
    integer                  :: np                               !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64)        :: xp                               !! Coordinate of the plane: x(np)
    type(fbem_qs_parameters) :: qsp                              !! Quasi-singular integration parameters
    integer                  :: ks                               !! Current level of subdivisions
    integer                  :: ns                               !! Maximum level of subdivision
    complex(kind=real64)        :: h(e%n_pnodes,3,3)                !! h integration kernels matrix
    complex(kind=real64)        :: g(e%n_snodes,3,3)                !! g integration kernels matrix
    ! Local
    real(kind=real64) :: xp_i(3)                              ! Image collocation point position vector
    integer           :: gln_near                             ! 1D Gauss-Legendre integ. points required to integrate only  the quasi-singular integrand
    integer           :: gln                                  ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
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
    complex(kind=real64) :: h_tmp(e%n_pnodes,3,3)                ! h integration kernels matrix (temporary)
    complex(kind=real64) :: g_tmp(e%n_snodes,3,3)                ! g integration kernels matrix (temporary)
    ! Collocation point image in global axes
    xp_i=x_i
    xp_i(abs(np))=2.*xp-xp_i(abs(np))
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
      call fbem_nearest_element_point_bem(3,e%gtype,e%x,e%cl,xp_i,barxi,rmin,d,method)
      barxip=barxi
    else
      call fbem_obtain_element_subdivision_coordinates(3,e%gtype,e%x,xi_s,x_s)
      cl=fbem_characteristic_length(3,e%gtype,x_s,1.d-12)
      call fbem_nearest_element_point_bem(3,e%gtype,x_s,cl,xp_i,barxip,rmin,d,method)
    end if
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,5,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_staela3d_hsc_sbie_ext_adp',ns,'maximum number of subdivisions reached')
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
          call fbem_bem_staela3d_hsc_sbie_ext_adp_cmplx(e,reverse,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,h,g)
          ! SUBTRI 2
          tmp_xi_s(:,1)=xi_s(:,2)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          call fbem_bem_staela3d_hsc_sbie_ext_adp_cmplx(e,reverse,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,h,g)
          ! SUBTRI 3
          tmp_xi_s(:,1)=xi_s(:,3)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          call fbem_bem_staela3d_hsc_sbie_ext_adp_cmplx(e,reverse,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,h,g)
          ! SUBTRI 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          call fbem_bem_staela3d_hsc_sbie_ext_adp_cmplx(e,reverse,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,h,g)
        ! QUADRILATERALS
        case (4)
          ! SUBQUAD 1
          tmp_xi_s(:,1)=xi_s(:,1)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,3)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          call fbem_bem_staela3d_hsc_sbie_ext_adp_cmplx(e,reverse,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,h,g)
          ! SUBQUAD 2
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=xi_s(:,2)
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,4)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          call fbem_bem_staela3d_hsc_sbie_ext_adp_cmplx(e,reverse,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,h,g)
          ! SUBQUAD 3
          tmp_xi_s(:,1)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=xi_s(:,3)
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          call fbem_bem_staela3d_hsc_sbie_ext_adp_cmplx(e,reverse,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,h,g)
          ! SUBQUAD 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          tmp_xi_s(:,2)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=xi_s(:,4)
          call fbem_bem_staela3d_hsc_sbie_ext_adp_cmplx(e,reverse,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,h,g)
      end select
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_staela3d_hsc_sbie_ext_st_cmplx(e,reverse,xi_s,x_i,barxip,barr,mu,nu,np,xp,gln,h_tmp,g_tmp)
      h=h+h_tmp
      g=g+g_tmp
    end if
  end subroutine fbem_bem_staela3d_hsc_sbie_ext_adp_cmplx

  subroutine fbem_bem_staela3d_hsc_sbie_int(type_g,type_f1,type_f2,delta_f,x_nodes,reverse,xi_i,mu,nu,np,xp,h,g)
    implicit none
    ! I/O
    integer           :: type_g                          !! Geometrial interpolation
    integer           :: type_f1                         !! Functional interpolation (primary variables)
    integer           :: type_f2                         !! Functional interpolation (secondary variables)
    real(kind=real64) :: delta_f                         !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64) :: x_nodes(3,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
    logical           :: reverse                         !! Normal vector inversion
    real(kind=real64) :: xi_i(2)                         !! Reference coordinates of the singular point.
    real(kind=real64) :: mu                              !! Shear modulus
    real(kind=real64) :: nu                              !! Poisson's ratio
    integer           :: np                              !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64) :: xp                              !! Coordinate of the plane: x(np)
    real(kind=real64) :: h(fbem_n_nodes(type_f1),3,3)    !! h kernel vector
    real(kind=real64) :: g(fbem_n_nodes(type_f2),3,3)    !! g kernel vector
    ! Local
    real(kind=real64) :: hk(fbem_n_nodes(type_f1),3,3)    ! h kernel vector (full-space)
    real(kind=real64) :: gk(fbem_n_nodes(type_f2),3,3)    ! g kernel vector (full-space)
    integer           :: ksubtri                          ! Counter variable for subtriangles loop
    integer           :: nsubtri                          ! Number of subtriangles
    integer           :: subtriangle(8)                   ! Vector that contains what subtriangles need to be integrated
    real(kind=real64) :: theta_subtri(2,8)                ! Matrix that contains the angles theta of the subtriangles to be integrated
    real(kind=real64) :: thetap_subtri(2,8)               ! Matrix that contains the angles thetap of the subtriangles to be integrated
    integer           :: ktheta                           ! Counter variable for theta coordinate loop
    integer           :: krho                             ! Counter variable for rho coordinate loop
    integer           :: kphi                             ! Counter coordinates for shape functions loops
    integer           :: nnodes_g                         ! Number of nodes of geometrical interpolation
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
    real(kind=real64) :: phi_f1(fbem_n_nodes(type_f1))    ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phi_f2(fbem_n_nodes(type_f2))    ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phi_f1_i(fbem_n_nodes(type_f1))  ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phi_g(fbem_n_nodes(type_g))      ! Geometrical shape functions values at xi_1,xi_2
    real(kind=real64) :: dphidxi1_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64) :: dphidxi2_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64) :: x_i(3)                           ! Collocation point position vector
    real(kind=real64) :: x(3)                             ! Position vector at xi_1,xi_2
    real(kind=real64) :: T1(3), T2(3)                     ! Tangent vectors at xi_1,xi_2
    real(kind=real64) :: N(3)                             ! Normal vector at xi_1,xi_2
    real(kind=real64) :: rv(3)                            ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64) :: r, d1r, d1r2                     ! Distance vector module and its inverse
    real(kind=real64) :: drdx(3)                          ! Distance vector derivatives with respect to x_k
    real(kind=real64) :: jg                               ! Geometric jacobian
    real(kind=real64) :: jw                               ! Jacobian * weights
    real(kind=real64) :: drdn                             ! Partial derivative of r respect to unit normal
    real(kind=real64) :: phif1jw(fbem_n_nodes(type_f1))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phif2jw(fbem_n_nodes(type_f2))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: costheta, sintheta               ! cos(theta), sin(theta)
    integer           :: il, ik                           ! Counter for load / observation components
    !
    ! Initialization
    !
    ! Kernel vectors
    h=0.d0
    g=0.d0
    ! Number of nodes of geometrical interpolation
    nnodes_g=fbem_n_nodes(type_g)
    !
    ! Checkings
    do kphi=1,nnodes_g
      if (abs(x_nodes(3,kphi)-xp).gt.1.d-12) then
        write(error_unit,*) 'Error: x must belong to the half-space (3)'
        stop
      end if
    end do
    ! Calculate full-space
    call fbem_bem_staela3d_sbie_int(type_g,type_f1,type_f2,delta_f,x_nodes,reverse,xi_i,mu,nu,hk,gk)
    ! At the free-surface, the only not null part is equal to the full-space part with opposite sign (which makes null the
    ! complete half-space fundamental solution in terms of traction at the free-surface)
    h=-hk
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
          jg=sqrt(dot_product(N,N))
          ! Unit normal vector
          n=N/jg
          ! Distance vector
          rv=x-x_i
          ! Distance vector norm
          r=dot_product(rv,rv)
          r=sqrt(r)
          d1r=1.d0/r
          ! r_{,k}
          drdx=rv*d1r
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
          ! Loop through load direction and observation direction
          g(:,1,1)=g(:,1,1)+(d1r+(3.-4.*nu)*drdx(1)**2*d1r+4.*(1.-nu)*(1.-2.*nu)/(r+x(3)-xp)*(1.-drdx(1)**2*r/(r+x(3)-xp)))*phif2jw
          g(:,1,2)=g(:,1,2)+((3.-4.*nu)*drdx(1)*drdx(2)*d1r-4.*(1.-nu)*(1.-2.*nu)*drdx(1)*drdx(2)*r/(r+x(3)-xp)**2)*phif2jw
          g(:,1,3)=g(:,1,3)+((3.-4.*nu)*drdx(1)*drdx(3)*d1r+4.*(1.-nu)*(1.-2.*nu)*drdx(1)/(r+x(3)-xp))*phif2jw
          g(:,2,1)=g(:,2,1)+((3.-4.*nu)*drdx(1)*drdx(2)*d1r-4.*(1.-nu)*(1.-2.*nu)*drdx(1)*drdx(2)*r/(r+x(3)-xp)**2)*phif2jw
          g(:,2,2)=g(:,2,2)+(d1r+(3.-4.*nu)*drdx(2)**2*d1r+4.*(1.-nu)*(1.-2.*nu)/(r+x(3)-xp)*(1.-drdx(2)**2*r/(r+x(3)-xp)))*phif2jw
          g(:,2,3)=g(:,2,3)+((3.-4.*nu)*drdx(2)*drdx(3)*d1r+4.*(1.-nu)*(1.-2.*nu)*drdx(2)/(r+x(3)-xp))*phif2jw
          g(:,3,1)=g(:,3,1)+((3.-4.*nu)*drdx(1)*drdx(3)*d1r-4.*(1.-nu)*(1.-2.*nu)*drdx(1)/(r+x(3)-xp))*phif2jw
          g(:,3,2)=g(:,3,2)+((3.-4.*nu)*drdx(2)*drdx(3)*d1r-4.*(1.-nu)*(1.-2.*nu)*drdx(2)/(r+x(3)-xp))*phif2jw
          g(:,3,3)=g(:,3,3)+((8.*(1.-nu)**2-(3.-4.*nu))*d1r+(3.-4.*nu)*drdx(3)**2*d1r)*phif2jw
        end do ! Loop through rho coordinate
      end do ! Loop through theta coordinate
    end do ! Loop through rectangle triangles
    g=g/(16.*c_pi*mu*(1.-nu))
  end subroutine fbem_bem_staela3d_hsc_sbie_int

  subroutine fbem_bem_staela3d_hsc_sbie_int_cmplx(type_g,type_f1,type_f2,delta_f,x_nodes,reverse,xi_i,mu,nu,np,xp,h,g)
    implicit none
    ! I/O
    integer           :: type_g                          !! Geometrial interpolation
    integer           :: type_f1                         !! Functional interpolation (primary variables)
    integer           :: type_f2                         !! Functional interpolation (secondary variables)
    real(kind=real64) :: delta_f                         !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64) :: x_nodes(3,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
    logical           :: reverse                         !! Normal vector inversion
    real(kind=real64) :: xi_i(2)                         !! Reference coordinates of the singular point.
    complex(kind=real64) :: mu                              !! Shear modulus
    complex(kind=real64) :: nu                              !! Poisson's ratio
    integer           :: np                              !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64) :: xp                              !! Coordinate of the plane: x(np)
    complex(kind=real64) :: h(fbem_n_nodes(type_f1),3,3)    !! h kernel vector
    complex(kind=real64) :: g(fbem_n_nodes(type_f2),3,3)    !! g kernel vector
    ! Local
    complex(kind=real64) :: hk(fbem_n_nodes(type_f1),3,3)    ! h kernel vector (full-space)
    complex(kind=real64) :: gk(fbem_n_nodes(type_f2),3,3)    ! g kernel vector (full-space)
    integer           :: ksubtri                          ! Counter variable for subtriangles loop
    integer           :: nsubtri                          ! Number of subtriangles
    integer           :: subtriangle(8)                   ! Vector that contains what subtriangles need to be integrated
    real(kind=real64) :: theta_subtri(2,8)                ! Matrix that contains the angles theta of the subtriangles to be integrated
    real(kind=real64) :: thetap_subtri(2,8)               ! Matrix that contains the angles thetap of the subtriangles to be integrated
    integer           :: ktheta                           ! Counter variable for theta coordinate loop
    integer           :: krho                             ! Counter variable for rho coordinate loop
    integer           :: kphi                             ! Counter coordinates for shape functions loops
    integer           :: nnodes_g                         ! Number of nodes of geometrical interpolation
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
    real(kind=real64) :: phi_f1(fbem_n_nodes(type_f1))    ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phi_f2(fbem_n_nodes(type_f2))    ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phi_f1_i(fbem_n_nodes(type_f1))  ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phi_g(fbem_n_nodes(type_g))      ! Geometrical shape functions values at xi_1,xi_2
    real(kind=real64) :: dphidxi1_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64) :: dphidxi2_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64) :: x_i(3)                           ! Collocation point position vector
    real(kind=real64) :: x(3)                             ! Position vector at xi_1,xi_2
    real(kind=real64) :: T1(3), T2(3)                     ! Tangent vectors at xi_1,xi_2
    real(kind=real64) :: N(3)                             ! Normal vector at xi_1,xi_2
    real(kind=real64) :: rv(3)                            ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64) :: r, d1r, d1r2                     ! Distance vector module and its inverse
    real(kind=real64) :: drdx(3)                          ! Distance vector derivatives with respect to x_k
    real(kind=real64) :: jg                               ! Geometric jacobian
    real(kind=real64) :: jw                               ! Jacobian * weights
    real(kind=real64) :: drdn                             ! Partial derivative of r respect to unit normal
    real(kind=real64) :: phif1jw(fbem_n_nodes(type_f1))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phif2jw(fbem_n_nodes(type_f2))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: costheta, sintheta               ! cos(theta), sin(theta)
    integer           :: il, ik                           ! Counter for load / observation components
    !
    ! Initialization
    !
    ! Kernel vectors
    h=0.d0
    g=0.d0
    ! Number of nodes of geometrical interpolation
    nnodes_g=fbem_n_nodes(type_g)
    !
    ! Checkings
    do kphi=1,nnodes_g
      if (abs(x_nodes(3,kphi)-xp).gt.1.d-12) then
        write(error_unit,*) 'Error: x must belong to the half-space (3)'
        stop
      end if
    end do
    ! Calculate full-space
    call fbem_bem_staela3d_sbie_int_cmplx(type_g,type_f1,type_f2,delta_f,x_nodes,reverse,xi_i,mu,nu,hk,gk)
    ! At the free-surface, the only not null part is equal to the full-space part with opposite sign (which makes null the
    ! complete half-space fundamental solution in terms of traction at the free-surface)
    h=-hk
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
          jg=sqrt(dot_product(N,N))
          ! Unit normal vector
          n=N/jg
          ! Distance vector
          rv=x-x_i
          ! Distance vector norm
          r=dot_product(rv,rv)
          r=sqrt(r)
          d1r=1.d0/r
          ! r_{,k}
          drdx=rv*d1r
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
          ! Loop through load direction and observation direction
          g(:,1,1)=g(:,1,1)+(d1r+(3.-4.*nu)*drdx(1)**2*d1r+4.*(1.-nu)*(1.-2.*nu)/(r+x(3)-xp)*(1.-drdx(1)**2*r/(r+x(3)-xp)))*phif2jw
          g(:,1,2)=g(:,1,2)+((3.-4.*nu)*drdx(1)*drdx(2)*d1r-4.*(1.-nu)*(1.-2.*nu)*drdx(1)*drdx(2)*r/(r+x(3)-xp)**2)*phif2jw
          g(:,1,3)=g(:,1,3)+((3.-4.*nu)*drdx(1)*drdx(3)*d1r+4.*(1.-nu)*(1.-2.*nu)*drdx(1)/(r+x(3)-xp))*phif2jw
          g(:,2,1)=g(:,2,1)+((3.-4.*nu)*drdx(1)*drdx(2)*d1r-4.*(1.-nu)*(1.-2.*nu)*drdx(1)*drdx(2)*r/(r+x(3)-xp)**2)*phif2jw
          g(:,2,2)=g(:,2,2)+(d1r+(3.-4.*nu)*drdx(2)**2*d1r+4.*(1.-nu)*(1.-2.*nu)/(r+x(3)-xp)*(1.-drdx(2)**2*r/(r+x(3)-xp)))*phif2jw
          g(:,2,3)=g(:,2,3)+((3.-4.*nu)*drdx(2)*drdx(3)*d1r+4.*(1.-nu)*(1.-2.*nu)*drdx(2)/(r+x(3)-xp))*phif2jw
          g(:,3,1)=g(:,3,1)+((3.-4.*nu)*drdx(1)*drdx(3)*d1r-4.*(1.-nu)*(1.-2.*nu)*drdx(1)/(r+x(3)-xp))*phif2jw
          g(:,3,2)=g(:,3,2)+((3.-4.*nu)*drdx(2)*drdx(3)*d1r-4.*(1.-nu)*(1.-2.*nu)*drdx(2)/(r+x(3)-xp))*phif2jw
          g(:,3,3)=g(:,3,3)+((8.*(1.-nu)**2-(3.-4.*nu))*d1r+(3.-4.*nu)*drdx(3)**2*d1r)*phif2jw
        end do ! Loop through rho coordinate
      end do ! Loop through theta coordinate
    end do ! Loop through rectangle triangles
    g=g/(16.*c_pi*mu*(1.-nu))
  end subroutine fbem_bem_staela3d_hsc_sbie_int_cmplx

  !! Efficient automatic integration
  subroutine fbem_bem_staela3d_hsc_sbie_auto(e,reverse,x_i,mu,nu,np,xp,qsp,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                !! Integration element
    logical                  :: reverse                          !! Reverse orientation
    real(kind=real64)        :: x_i(3)                           !! Collocation point
    real(kind=real64)        :: mu                               !! Shear modulus
    real(kind=real64)        :: nu                               !! Poisson's ratio
    integer                  :: np                               !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64)        :: xp                               !! Coordinate of the plane: x(np)
    type(fbem_qs_parameters) :: qsp                              !! Quasi-singular integration parameters
    integer                  :: ns                               !! Maximum level of subdivisions
    real(kind=real64)        :: h(e%n_pnodes,3,3)                !! h integration kernel
    real(kind=real64)        :: g(e%n_snodes,3,3)                !! g integration kernel
    ! Local
    real(kind=real64) :: xp_i(3)                          ! Collocation point image in global axes
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
    ! Collocation point image in global axes
    xp_i=x_i
    xp_i(abs(np))=2.*xp-xp_i(abs(np))
    ! Determine if interior or exterior integration
    !   - Interior integration (delta=1) requires: xi_i
    !   - Exterior integration (delta=0) requires: xp_i, barxi, rmin and d
    ! Use the element ball
    r=e%bball_centre-xp_i
    rmin=sqrt(dot_product(r,r))-e%bball_radius
    if (rmin.gt.(4.d0*e%bball_radius)) then
      delta=0
      barxi=0.d0
      d=rmin/e%cl
    else
      ! Use an adaptative algorithm that combines sampling and minimization algorithms
      call fbem_nearest_element_point_bem(3,e%gtype,e%x,e%cl,xp_i,barxi,rmin,d,method)
      if (d.le.1.d-12) then
        delta=1
      else
        delta=0
      end if
    end if
    ! Integrate
    select case (delta)
      case (1)
        call fbem_bem_staela3d_hsc_sbie_int(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi,mu,nu,np,xp,h,g)
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
          call fbem_bem_staela3d_hsc_sbie_ext_pre(ps,e,reverse,x_i,mu,nu,np,xp,h,g)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_staela3d_hsc_sbie_ext_adp(e,reverse,xi_s,x_i,mu,nu,np,xp,qsp,1,ns,h,g)
        end if
    end select
  end subroutine fbem_bem_staela3d_hsc_sbie_auto

    !! Efficient automatic integration
  subroutine fbem_bem_staela3d_hsc_sbie_auto_cmplx(e,reverse,x_i,mu,nu,np,xp,qsp,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                !! Integration element
    logical                  :: reverse                          !! Reverse orientation
    real(kind=real64)        :: x_i(3)                           !! Collocation point
    complex(kind=real64)        :: mu                               !! Shear modulus
    complex(kind=real64)        :: nu                               !! Poisson's ratio
    integer                  :: np                               !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64)        :: xp                               !! Coordinate of the plane: x(np)
    type(fbem_qs_parameters) :: qsp                              !! Quasi-singular integration parameters
    integer                  :: ns                               !! Maximum level of subdivisions
    complex(kind=real64)        :: h(e%n_pnodes,3,3)                !! h integration kernel
    complex(kind=real64)        :: g(e%n_snodes,3,3)                !! g integration kernel
    ! Local
    real(kind=real64) :: xp_i(3)                          ! Collocation point image in global axes
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
    ! Collocation point image in global axes
    xp_i=x_i
    xp_i(abs(np))=2.*xp-xp_i(abs(np))
    ! Determine if interior or exterior integration
    !   - Interior integration (delta=1) requires: xi_i
    !   - Exterior integration (delta=0) requires: xp_i, barxi, rmin and d
    ! Use the element ball
    r=e%bball_centre-xp_i
    rmin=sqrt(dot_product(r,r))-e%bball_radius
    if (rmin.gt.(4.d0*e%bball_radius)) then
      delta=0
      barxi=0.d0
      d=rmin/e%cl
    else
      ! Use an adaptative algorithm that combines sampling and minimization algorithms
      call fbem_nearest_element_point_bem(3,e%gtype,e%x,e%cl,xp_i,barxi,rmin,d,method)
      if (d.le.1.d-12) then
        delta=1
      else
        delta=0
      end if
    end if
    ! Integrate
    select case (delta)
      case (1)
        call fbem_bem_staela3d_hsc_sbie_int_cmplx(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi,mu,nu,np,xp,h,g)
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
          call fbem_bem_staela3d_hsc_sbie_ext_pre_cmplx(ps,e,reverse,x_i,mu,nu,np,xp,h,g)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_staela3d_hsc_sbie_ext_adp_cmplx(e,reverse,xi_s,x_i,mu,nu,np,xp,qsp,1,ns,h,g)
        end if
    end select
  end subroutine fbem_bem_staela3d_hsc_sbie_auto_cmplx

  ! --------------------- !
  ! BE BODY LOAD ELEMENTS !
  ! --------------------- !

  subroutine fbem_bem_staela3d_hsc_sbie_bl_ext_pre(ps,e,x_i,mu,nu,np,xp,g)
    implicit none
    ! I/O
    integer                :: ps                !! Selected precalculated dataset
    type(fbem_bem_element) :: e                 !! Element
    real(kind=real64)      :: x_i(3)            !! Collocation point position vector
    real(kind=real64)      :: mu                !! Shear modulus
    real(kind=real64)      :: nu                !! Poisson's ratio
    integer                :: np                !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64)      :: xp                !! Coordinate of the plane: x(np)
    real(kind=real64)      :: g(e%n_snodes,3,3) !! g integration kernels matrix
    ! Local
    integer           :: il, ik               ! Counter for load / observation components
    integer           :: kip                  ! Counter variable for integration points loop
    real(kind=real64) :: x(3)                 ! Position vector at integration point
    real(kind=real64) :: n(3)                 ! Unit normal vector at integration point
    real(kind=real64) :: sphijw(e%n_snodes)   ! phi^s * jacobian * weight at integration point
    real(kind=real64) :: fs_u(3,3), fs_t(3,3) ! Fundamental solution
    ! Initialize
    g=0.
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      x=e%ps_x(:,kip,ps)
      n=e%ps_n(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      call fbem_bem_staela3d_hsc_sbie(mu,nu,np,xp,x_i,x,n,fs_u,fs_t)
      do il=1,3
        do ik=1,3
          g(:,il,ik)=g(:,il,ik)+fs_u(il,ik)*sphijw
        end do
      end do
    end do
  end subroutine fbem_bem_staela3d_hsc_sbie_bl_ext_pre

  subroutine fbem_bem_staela3d_hsc_sbie_bl_ext_pre_cmplx(ps,e,x_i,mu,nu,np,xp,g)
    implicit none
    ! I/O
    integer                :: ps                !! Selected precalculated dataset
    type(fbem_bem_element) :: e                 !! Element
    real(kind=real64)      :: x_i(3)            !! Collocation point position vector
    complex(kind=real64)      :: mu                !! Shear modulus
    complex(kind=real64)      :: nu                !! Poisson's ratio
    integer                :: np                !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64)      :: xp                !! Coordinate of the plane: x(np)
    complex(kind=real64)      :: g(e%n_snodes,3,3) !! g integration kernels matrix
    ! Local
    integer           :: il, ik               ! Counter for load / observation components
    integer           :: kip                  ! Counter variable for integration points loop
    real(kind=real64) :: x(3)                 ! Position vector at integration point
    real(kind=real64) :: n(3)                 ! Unit normal vector at integration point
    real(kind=real64) :: sphijw(e%n_snodes)   ! phi^s * jacobian * weight at integration point
    complex(kind=real64) :: fs_u(3,3), fs_t(3,3) ! Fundamental solution
    ! Initialize
    g=0.
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      x=e%ps_x(:,kip,ps)
      n=e%ps_n(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      call fbem_bem_staela3d_hsc_sbie_cmplx(mu,nu,np,xp,x_i,x,n,fs_u,fs_t)
      do il=1,3
        do ik=1,3
          g(:,il,ik)=g(:,il,ik)+fs_u(il,ik)*sphijw
        end do
      end do
    end do
  end subroutine fbem_bem_staela3d_hsc_sbie_bl_ext_pre_cmplx

  !! This subroutine calculates the kernels for SBIE exterior integration (near collocation points) using Telles transformation
  !! within a subdivision of the element, needing only needs basic data.
  subroutine fbem_bem_staela3d_hsc_sbie_bl_ext_st(e,xi_s,x_i,barxip,barr,mu,nu,np,xp,gln,g)
    implicit none
    ! I/O
    type(fbem_bem_element)       :: e                                !! Integration element
    real(kind=real64)            :: xi_s(2,fbem_n_vertices(e%gtype)) !! Coordinates of the subdivision of the parent element
    real(kind=real64)            :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)            :: barxip(2)                        !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)            :: barr                             !! Telles jacobian at barxip
    real(kind=real64)            :: mu                               !! Shear modulus
    real(kind=real64)            :: nu                               !! Poisson's ratio
    integer                      :: np                               !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64)            :: xp                               !! Coordinate of the plane: x(np)
    integer                      :: gln                              !! 1D Gauss-Legendre number of integration points (<=32)
    real(kind=real64)            :: g(e%n_snodes,3,3)                !! g kernel vector
    ! Local
    integer                      :: kphi                       ! Counter variable for shape functions loops
    integer                      :: k1                         ! Counter variable for reference coordinate xi_1
    integer                      :: k2                         ! Counter variable for reference coordinate xi_2
    real(kind=real64)            :: aux(10)                    ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)           ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi (e%n_gnodes)      ! Geometrical shape functions derivatives values
    real(kind=real64)            :: dgphidxi1(e%n_gnodes)      ! Geometrical shape functions derivatives values
    real(kind=real64)            :: dgphidxi2(e%n_gnodes)      ! Geometrical shape functions derivatives values
    real(kind=real64)            :: sphi(e%n_snodes)           ! Functional shape functions values
    real(kind=real64)            :: gamma(e%d)                 ! Vector of gamma_1,gamma_2 coordinates
    real(kind=real64)            :: w(e%d)                     ! Weights of the integration rule
    real(kind=real64)            :: xip(e%d)                   ! Vector of xip_1,xip_2 coordinates
    real(kind=real64)            :: dxidxi1p(2)                ! xi derivatives with respect to xi1p
    real(kind=real64)            :: dxidxi2p(2)                ! xi derivatives with respect to xi2p
    real(kind=real64)            :: js                         ! Subdivision jacobian
    real(kind=real64)            :: xi(e%d)                    ! Vector of xi_1,xi_2 coordinates
    real(kind=real64)            :: xipp(2)                    ! Coordinate xipp used for quadrilateral-triangle transformation
    real(kind=real64)            :: barxipp(2)                 ! Coordinate xipp of collocation point
    real(kind=real64)            :: jqt                        ! Jacobian of the quadrilateral-triangle transformation
    type(fbem_telles_parameters) :: telles_parameters(e%d)     ! Telles parameters for each coordinate
    real(kind=real64)            :: jt(e%d)                    ! Telles jacobian for each coordinate: xi_1->gamma_1 and xi_2->gamma_2
    real(kind=real64)            :: x(3)                       ! Position vector at xi_1,xi_2
    real(kind=real64)            :: T(3), T1(3), T2(3), N(3)   ! Tangent vectors
    real(kind=real64)            :: jg                         ! Geometric jacobian
    real(kind=real64)            :: jw                         ! Jacobians * weights
    real(kind=real64)            :: sphijw(e%n_snodes)         ! Functional shape functions * jw
    integer                      :: il, ik                     ! Counter for load / observation components
    real(kind=real64)            :: fs_u(3,3), fs_t(3,3)       ! Fundamental solutions values
    ! Initialize
    g=0
    n=0
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
          xipp(1)=0.5d0*(1.d0-xip(1))*xi_s(1,1)+0.5d0*(1.d0+xip(1))*xi_s(1,2)
          js=0.5d0*(xi_s(1,2)-xi_s(1,1))
          ! XI->X TRANSFORMATION
#         define etype e%gtype
#         define delta 0.d0
#         define xi xipp(1)
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
          ! Jacobians * weight
          jw=jg*js*jt(1)*w(1)
          ! FUNCTIONAL SHAPE FUNCTIONS
          ! Functional shape functions (body load variables) at xi
#         define etype e%stype
#         define delta e%stype_delta
#         define xi xipp(1)
#         define phi sphi
#         include <phi_1d.rc>
#         undef etype
#         undef delta
#         undef xi
#         undef phi
          ! Functional shape functions * jacobians* weights
          sphijw=sphi*jw
          ! Fundamental solution
          call fbem_bem_staela3d_hsc_sbie(mu,nu,np,xp,x_i,x,n,fs_u,fs_t)
          ! Add integration points
          do il=1,3
            do ik=1,3
              g(:,il,ik)=g(:,il,ik)+fs_u(il,ik)*sphijw
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
#               define etype e%gtype
#               define delta 0.d0
#               define phi gphi
#               define dphidxi1 dgphidxi1
#               define dphidxi2 dgphidxi2
#               include <phi_and_dphidxik_2d.rc>
#               undef etype
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
                ! Jacobians * weights
                jw=jg*js*jt(1)*jt(2)*w(1)*w(2)
                ! FUNCTIONAL SHAPE FUNCTIONS
                ! Functional shape functions (secondary variables) at xi
#               define etype e%stype
#               define delta e%stype_delta
#               define phi sphi
#               include <phi_2d.rc>
#               undef etype
#               undef delta
#               undef phi
                ! Functional shape functions * jacobians * weights
                sphijw=sphi*jw
                ! Fundamental solution
                call fbem_bem_staela3d_hsc_sbie(mu,nu,np,xp,x_i,x,n,fs_u,fs_t)
                ! Add integration points
                do il=1,3
                  do ik=1,3
                    g(:,il,ik)=g(:,il,ik)+fs_u(il,ik)*sphijw
                  end do
                end do
              end do
            end do
          ! TRIANGULAR ELEMENTS
          case (3)
            ! Telles transformation is applied to Gauss-Legendre*Gauss-Legendre quadrature before the quadrilateral->triangle
            ! transformation. Because barxip for triangles are given in area triangle coordinates, for Telles transformation
            ! they must be transformed to quadrilateral coordinates. A special treatment is needed when barxi_2 is near 1, because
            ! transformation diverges.
            ! BARXIP->BARXIPP TRANSFORMATION
            if (barxip(2).gt.0.995d0) then
              barxipp(1)=0.5d0
              barxipp(2)=1.0d0
            else
              barxipp(1)=barxip(1)/(1.d0-barxip(2))
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
#               define etype e%gtype
#               define delta 0.d0
#               define phi gphi
#               define dphidxi1 dgphidxi1
#               define dphidxi2 dgphidxi2
#               include <phi_and_dphidxik_2d.rc>
#               undef etype
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
                ! Jacobians * weights
                jw=jg*js*jqt*jt(1)*jt(2)*w(1)*w(2)
                ! FUNCTIONAL SHAPE FUNCTIONS
                ! Functional shape functions (secondary variables) at xi
#               define etype e%stype
#               define delta e%stype_delta
#               define phi sphi
#               include <phi_2d.rc>
#               undef etype
#               undef delta
#               undef phi
                ! Functional shape functions * jacobians * weights
                sphijw=sphi*jw
                ! Fundamental solution
                call fbem_bem_staela3d_hsc_sbie(mu,nu,np,xp,x_i,x,n,fs_u,fs_t)
                ! Add integration points
                do il=1,3
                  do ik=1,3
                    g(:,il,ik)=g(:,il,ik)+fs_u(il,ik)*sphijw
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
        stop 'not yet : fbem_bem_staela3d_sbie_bl_ext_st'
      !
      ! OTHERS
      !
      case default
        call fbem_error_message(output_unit,0,'fbem_bem_staela3d_hsc_sbie_bl_ext_st',0,&
                                'it is only possible to integrate line, surface or volume loads')
    end select
  end subroutine fbem_bem_staela3d_hsc_sbie_bl_ext_st

  !! This subroutine calculates the kernels for SBIE exterior integration (near collocation points) using Telles transformation
  !! within a subdivision of the element, needing only needs basic data.
  subroutine fbem_bem_staela3d_hsc_sbie_bl_ext_st_cmplx(e,xi_s,x_i,barxip,barr,mu,nu,np,xp,gln,g)
    implicit none
    ! I/O
    type(fbem_bem_element)       :: e                                !! Integration element
    real(kind=real64)            :: xi_s(2,fbem_n_vertices(e%gtype)) !! Coordinates of the subdivision of the parent element
    real(kind=real64)            :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)            :: barxip(2)                        !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)            :: barr                             !! Telles jacobian at barxip
    complex(kind=real64)            :: mu                               !! Shear modulus
    complex(kind=real64)            :: nu                               !! Poisson's ratio
    integer                      :: np                               !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64)            :: xp                               !! Coordinate of the plane: x(np)
    integer                      :: gln                              !! 1D Gauss-Legendre number of integration points (<=32)
    complex(kind=real64)            :: g(e%n_snodes,3,3)                !! g kernel vector
    ! Local
    integer                      :: kphi                       ! Counter variable for shape functions loops
    integer                      :: k1                         ! Counter variable for reference coordinate xi_1
    integer                      :: k2                         ! Counter variable for reference coordinate xi_2
    real(kind=real64)            :: aux(10)                    ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)           ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi (e%n_gnodes)      ! Geometrical shape functions derivatives values
    real(kind=real64)            :: dgphidxi1(e%n_gnodes)      ! Geometrical shape functions derivatives values
    real(kind=real64)            :: dgphidxi2(e%n_gnodes)      ! Geometrical shape functions derivatives values
    real(kind=real64)            :: sphi(e%n_snodes)           ! Functional shape functions values
    real(kind=real64)            :: gamma(e%d)                 ! Vector of gamma_1,gamma_2 coordinates
    real(kind=real64)            :: w(e%d)                     ! Weights of the integration rule
    real(kind=real64)            :: xip(e%d)                   ! Vector of xip_1,xip_2 coordinates
    real(kind=real64)            :: dxidxi1p(2)                ! xi derivatives with respect to xi1p
    real(kind=real64)            :: dxidxi2p(2)                ! xi derivatives with respect to xi2p
    real(kind=real64)            :: js                         ! Subdivision jacobian
    real(kind=real64)            :: xi(e%d)                    ! Vector of xi_1,xi_2 coordinates
    real(kind=real64)            :: xipp(2)                    ! Coordinate xipp used for quadrilateral-triangle transformation
    real(kind=real64)            :: barxipp(2)                 ! Coordinate xipp of collocation point
    real(kind=real64)            :: jqt                        ! Jacobian of the quadrilateral-triangle transformation
    type(fbem_telles_parameters) :: telles_parameters(e%d)     ! Telles parameters for each coordinate
    real(kind=real64)            :: jt(e%d)                    ! Telles jacobian for each coordinate: xi_1->gamma_1 and xi_2->gamma_2
    real(kind=real64)            :: x(3)                       ! Position vector at xi_1,xi_2
    real(kind=real64)            :: T(3), T1(3), T2(3), N(3)   ! Tangent vectors
    real(kind=real64)            :: jg                         ! Geometric jacobian
    real(kind=real64)            :: jw                         ! Jacobians * weights
    real(kind=real64)            :: sphijw(e%n_snodes)         ! Functional shape functions * jw
    integer                      :: il, ik                     ! Counter for load / observation components
    complex(kind=real64)            :: fs_u(3,3), fs_t(3,3)       ! Fundamental solutions values
    ! Initialize
    g=0
    n=0
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
          xipp(1)=0.5d0*(1.d0-xip(1))*xi_s(1,1)+0.5d0*(1.d0+xip(1))*xi_s(1,2)
          js=0.5d0*(xi_s(1,2)-xi_s(1,1))
          ! XI->X TRANSFORMATION
#         define etype e%gtype
#         define delta 0.d0
#         define xi xipp(1)
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
          ! Jacobians * weight
          jw=jg*js*jt(1)*w(1)
          ! FUNCTIONAL SHAPE FUNCTIONS
          ! Functional shape functions (body load variables) at xi
#         define etype e%stype
#         define delta e%stype_delta
#         define xi xipp(1)
#         define phi sphi
#         include <phi_1d.rc>
#         undef etype
#         undef delta
#         undef xi
#         undef phi
          ! Functional shape functions * jacobians* weights
          sphijw=sphi*jw
          ! Fundamental solution
          call fbem_bem_staela3d_hsc_sbie_cmplx(mu,nu,np,xp,x_i,x,n,fs_u,fs_t)
          ! Add integration points
          do il=1,3
            do ik=1,3
              g(:,il,ik)=g(:,il,ik)+fs_u(il,ik)*sphijw
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
#               define etype e%gtype
#               define delta 0.d0
#               define phi gphi
#               define dphidxi1 dgphidxi1
#               define dphidxi2 dgphidxi2
#               include <phi_and_dphidxik_2d.rc>
#               undef etype
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
                ! Jacobians * weights
                jw=jg*js*jt(1)*jt(2)*w(1)*w(2)
                ! FUNCTIONAL SHAPE FUNCTIONS
                ! Functional shape functions (secondary variables) at xi
#               define etype e%stype
#               define delta e%stype_delta
#               define phi sphi
#               include <phi_2d.rc>
#               undef etype
#               undef delta
#               undef phi
                ! Functional shape functions * jacobians * weights
                sphijw=sphi*jw
                ! Fundamental solution
                call fbem_bem_staela3d_hsc_sbie_cmplx(mu,nu,np,xp,x_i,x,n,fs_u,fs_t)
                ! Add integration points
                do il=1,3
                  do ik=1,3
                    g(:,il,ik)=g(:,il,ik)+fs_u(il,ik)*sphijw
                  end do
                end do
              end do
            end do
          ! TRIANGULAR ELEMENTS
          case (3)
            ! Telles transformation is applied to Gauss-Legendre*Gauss-Legendre quadrature before the quadrilateral->triangle
            ! transformation. Because barxip for triangles are given in area triangle coordinates, for Telles transformation
            ! they must be transformed to quadrilateral coordinates. A special treatment is needed when barxi_2 is near 1, because
            ! transformation diverges.
            ! BARXIP->BARXIPP TRANSFORMATION
            if (barxip(2).gt.0.995d0) then
              barxipp(1)=0.5d0
              barxipp(2)=1.0d0
            else
              barxipp(1)=barxip(1)/(1.d0-barxip(2))
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
#               define etype e%gtype
#               define delta 0.d0
#               define phi gphi
#               define dphidxi1 dgphidxi1
#               define dphidxi2 dgphidxi2
#               include <phi_and_dphidxik_2d.rc>
#               undef etype
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
                ! Jacobians * weights
                jw=jg*js*jqt*jt(1)*jt(2)*w(1)*w(2)
                ! FUNCTIONAL SHAPE FUNCTIONS
                ! Functional shape functions (secondary variables) at xi
#               define etype e%stype
#               define delta e%stype_delta
#               define phi sphi
#               include <phi_2d.rc>
#               undef etype
#               undef delta
#               undef phi
                ! Functional shape functions * jacobians * weights
                sphijw=sphi*jw
                ! Fundamental solution
                call fbem_bem_staela3d_hsc_sbie_cmplx(mu,nu,np,xp,x_i,x,n,fs_u,fs_t)
                ! Add integration points
                do il=1,3
                  do ik=1,3
                    g(:,il,ik)=g(:,il,ik)+fs_u(il,ik)*sphijw
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
        stop 'not yet : fbem_bem_staela3d_sbie_bl_ext_st_cmplx'
      !
      ! OTHERS
      !
      case default
        call fbem_error_message(output_unit,0,'fbem_bem_staela3d_hsc_sbie_bl_ext_st_cmplx',0,&
                                'it is only possible to integrate line, surface or volume loads')
    end select
  end subroutine fbem_bem_staela3d_hsc_sbie_bl_ext_st_cmplx

  recursive subroutine fbem_bem_staela3d_hsc_sbie_bl_ext_adp(e,xi_s,x_i,mu,nu,np,xp,qsp,ks,ns,g)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                !! Element
    real(kind=real64)        :: xi_s(2,fbem_n_vertices(e%gtype)) !! Subdivision of the parent element
    real(kind=real64)        :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)        :: mu                               !! Shear modulus
    real(kind=real64)        :: nu                               !! Poisson's ratio
    integer                  :: np                               !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64)        :: xp                               !! Coordinate of the plane: x(np)
    type(fbem_qs_parameters) :: qsp                              !! Quasi-singular integration parameters
    integer                  :: ks                               !! Current level of subdivisions
    integer                  :: ns                               !! Maximum level of subdivision
    real(kind=real64)        :: g(e%n_snodes,3,3)                !! g integration kernels matrix
    ! Local
    real(kind=real64) :: xp_i(3)                              ! Image collocation point position vector
    integer           :: gln_near                             ! 1D Gauss-Legendre integ. points required to integrate only  the quasi-singular integrand
    integer           :: gln                                  ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical           :: subdivide                            ! True if subdivision has to be performed
    real(kind=real64) :: barxi(e%d)                           ! Nearest element coordinates with respect to collocation point
    real(kind=real64) :: barxip(e%d)                          ! Nearest element subdivision local coordinates with respect to collocation point
    real(kind=real64) :: rmin                                 ! Minimum distance between collocation point and barxi on the element
    real(kind=real64) :: barr                                 ! Telles jacobian at barxi
    real(kind=real64) :: cl                                   ! Characteristic length
    real(kind=real64) :: d                                    ! Normalized distance between collocation point and element subdivision
    integer           :: method                               ! Method used in nearest point algorithm
    real(kind=real64) :: tmp_xi_s(2,fbem_n_vertices(e%gtype)) ! Subdivision
    real(kind=real64) :: x_s(3,e%n_gnodes)                    ! Coordinates of the element subdivision
    real(kind=real64) :: g_tmp(e%n_snodes,3,3)                ! g integration kernels matrix (temporary)
    ! Collocation point image in global axes
    xp_i=x_i
    xp_i(abs(np))=2.*xp-xp_i(abs(np))
    ! Initialize
    if (ks.eq.1) then
      g=0.d0
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
          stop 'not yet : fbem_bem_staela3d_hsc_sbie_bl_ext_adp'
      end select
      call fbem_nearest_element_point_bem(3,e%gtype,e%x,e%cl,xp_i,barxi,rmin,d,method)
      barxip=barxi
    else
      call fbem_obtain_element_subdivision_coordinates(3,e%gtype,e%x,xi_s,x_s)
      cl=fbem_characteristic_length(3,e%gtype,x_s,1.d-12)
      call fbem_nearest_element_point_bem(3,e%gtype,x_s,cl,xp_i,barxip,rmin,d,method)
    end if
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,5,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_staela3d_hsc_sbie_bl_ext_adp',ns,'maximum number of subdivisions reached')
        gln_near=30
      end if
    else
      if (gln_near.eq.0) subdivide=.true.
    end if
    ! Subdivide
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
          call fbem_bem_staela3d_hsc_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,g)
          ! SUBLINE 2
          tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
          tmp_xi_s(1,2)=xi_s(1,2)
          call fbem_bem_staela3d_hsc_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,g)
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
              call fbem_bem_staela3d_hsc_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,g)
              ! SUBTRI 2
              tmp_xi_s(:,1)=xi_s(:,2)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              call fbem_bem_staela3d_hsc_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,g)
              ! SUBTRI 3
              tmp_xi_s(:,1)=xi_s(:,3)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              call fbem_bem_staela3d_hsc_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,g)
              ! SUBTRI 4
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
              call fbem_bem_staela3d_hsc_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,g)
            ! QUADRILATERALS
            case (4)
              ! SUBQUAD 1
              tmp_xi_s(:,1)=xi_s(:,1)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,3)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,4)=0.50d0*(xi_s(:,1)+xi_s(:,4))
              call fbem_bem_staela3d_hsc_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,g)
              ! SUBQUAD 2
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,2)=xi_s(:,2)
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,4)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              call fbem_bem_staela3d_hsc_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,g)
              ! SUBQUAD 3
              tmp_xi_s(:,1)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=xi_s(:,3)
              tmp_xi_s(:,4)=0.50d0*(xi_s(:,3)+xi_s(:,4))
              call fbem_bem_staela3d_hsc_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,g)
              ! SUBQUAD 4
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,4))
              tmp_xi_s(:,2)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,4)=xi_s(:,4)
              call fbem_bem_staela3d_hsc_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,g)
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
      call fbem_bem_staela3d_hsc_sbie_bl_ext_st(e,xi_s,x_i,barxip,barr,mu,nu,np,xp,gln,g_tmp)
      g=g+g_tmp
    end if
  end subroutine fbem_bem_staela3d_hsc_sbie_bl_ext_adp

  recursive subroutine fbem_bem_staela3d_hsc_sbie_bl_ext_adp_cmplx(e,xi_s,x_i,mu,nu,np,xp,qsp,ks,ns,g)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                !! Element
    real(kind=real64)        :: xi_s(2,fbem_n_vertices(e%gtype)) !! Subdivision of the parent element
    real(kind=real64)        :: x_i(3)                           !! Collocation point position vector
    complex(kind=real64)        :: mu                               !! Shear modulus
    complex(kind=real64)        :: nu                               !! Poisson's ratio
    integer                  :: np                               !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64)        :: xp                               !! Coordinate of the plane: x(np)
    type(fbem_qs_parameters) :: qsp                              !! Quasi-singular integration parameters
    integer                  :: ks                               !! Current level of subdivisions
    integer                  :: ns                               !! Maximum level of subdivision
    complex(kind=real64)        :: g(e%n_snodes,3,3)                !! g integration kernels matrix
    ! Local
    real(kind=real64) :: xp_i(3)                              ! Image collocation point position vector
    integer           :: gln_near                             ! 1D Gauss-Legendre integ. points required to integrate only  the quasi-singular integrand
    integer           :: gln                                  ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical           :: subdivide                            ! True if subdivision has to be performed
    real(kind=real64) :: barxi(e%d)                           ! Nearest element coordinates with respect to collocation point
    real(kind=real64) :: barxip(e%d)                          ! Nearest element subdivision local coordinates with respect to collocation point
    real(kind=real64) :: rmin                                 ! Minimum distance between collocation point and barxi on the element
    real(kind=real64) :: barr                                 ! Telles jacobian at barxi
    real(kind=real64) :: cl                                   ! Characteristic length
    real(kind=real64) :: d                                    ! Normalized distance between collocation point and element subdivision
    integer           :: method                               ! Method used in nearest point algorithm
    real(kind=real64) :: tmp_xi_s(2,fbem_n_vertices(e%gtype)) ! Subdivision
    real(kind=real64) :: x_s(3,e%n_gnodes)                    ! Coordinates of the element subdivision
    complex(kind=real64) :: g_tmp(e%n_snodes,3,3)                ! g integration kernels matrix (temporary)
    ! Collocation point image in global axes
    xp_i=x_i
    xp_i(abs(np))=2.*xp-xp_i(abs(np))
    ! Initialize
    if (ks.eq.1) then
      g=0.d0
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
          stop 'not yet : fbem_bem_staela3d_hsc_sbie_bl_ext_adp'
      end select
      call fbem_nearest_element_point_bem(3,e%gtype,e%x,e%cl,xp_i,barxi,rmin,d,method)
      barxip=barxi
    else
      call fbem_obtain_element_subdivision_coordinates(3,e%gtype,e%x,xi_s,x_s)
      cl=fbem_characteristic_length(3,e%gtype,x_s,1.d-12)
      call fbem_nearest_element_point_bem(3,e%gtype,x_s,cl,xp_i,barxip,rmin,d,method)
    end if
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,5,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_staela3d_hsc_sbie_bl_ext_adp_cmplx',ns,'maximum number of subdivisions reached')
        gln_near=30
      end if
    else
      if (gln_near.eq.0) subdivide=.true.
    end if
    ! Subdivide
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
          call fbem_bem_staela3d_hsc_sbie_bl_ext_adp_cmplx(e,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,g)
          ! SUBLINE 2
          tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
          tmp_xi_s(1,2)=xi_s(1,2)
          call fbem_bem_staela3d_hsc_sbie_bl_ext_adp_cmplx(e,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,g)
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
              call fbem_bem_staela3d_hsc_sbie_bl_ext_adp_cmplx(e,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,g)
              ! SUBTRI 2
              tmp_xi_s(:,1)=xi_s(:,2)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              call fbem_bem_staela3d_hsc_sbie_bl_ext_adp_cmplx(e,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,g)
              ! SUBTRI 3
              tmp_xi_s(:,1)=xi_s(:,3)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              call fbem_bem_staela3d_hsc_sbie_bl_ext_adp_cmplx(e,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,g)
              ! SUBTRI 4
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
              call fbem_bem_staela3d_hsc_sbie_bl_ext_adp_cmplx(e,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,g)
            ! QUADRILATERALS
            case (4)
              ! SUBQUAD 1
              tmp_xi_s(:,1)=xi_s(:,1)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,3)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,4)=0.50d0*(xi_s(:,1)+xi_s(:,4))
              call fbem_bem_staela3d_hsc_sbie_bl_ext_adp_cmplx(e,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,g)
              ! SUBQUAD 2
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,2)=xi_s(:,2)
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,4)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              call fbem_bem_staela3d_hsc_sbie_bl_ext_adp_cmplx(e,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,g)
              ! SUBQUAD 3
              tmp_xi_s(:,1)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=xi_s(:,3)
              tmp_xi_s(:,4)=0.50d0*(xi_s(:,3)+xi_s(:,4))
              call fbem_bem_staela3d_hsc_sbie_bl_ext_adp_cmplx(e,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,g)
              ! SUBQUAD 4
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,4))
              tmp_xi_s(:,2)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,4)=xi_s(:,4)
              call fbem_bem_staela3d_hsc_sbie_bl_ext_adp_cmplx(e,tmp_xi_s,x_i,mu,nu,np,xp,qsp,ks+1,ns,g)
          end select
        !
        ! VOLUME LOAD
        !
        case (3)
          stop 'not yet : fbem_bem_staela3d_sbie_bl_ext_adp_cmplx'
      end select
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_staela3d_hsc_sbie_bl_ext_st_cmplx(e,xi_s,x_i,barxip,barr,mu,nu,np,xp,gln,g_tmp)
      g=g+g_tmp
    end if
  end subroutine fbem_bem_staela3d_hsc_sbie_bl_ext_adp_cmplx

  subroutine fbem_bem_staela3d_hsc_sbie_bl_int(type_g,type_f1,type_f2,delta_f,x_nodes,xi_i,mu,nu,np,xp,g)
    implicit none
    ! I/O
    integer           :: type_g                          !! Geometrial interpolation
    integer           :: type_f1                         !! Functional interpolation (primary variables)
    integer           :: type_f2                         !! Functional interpolation (secondary variables)
    real(kind=real64) :: delta_f                         !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64) :: x_nodes(3,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
    real(kind=real64) :: xi_i(2)                         !! Reference coordinates of the singular point.
    real(kind=real64) :: mu                              !! Shear modulus
    real(kind=real64) :: nu                              !! Poisson's ratio
    integer           :: np                              !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64) :: xp                              !! Coordinate of the plane: x(np)
    real(kind=real64) :: g(fbem_n_nodes(type_f2),3,3)    !! g kernel vector
    ! Local
    real(kind=real64) :: hk(fbem_n_nodes(type_f1),3,3)    ! h kernel vector (full-space)
    real(kind=real64) :: gk(fbem_n_nodes(type_f2),3,3)    ! g kernel vector (full-space)
    integer           :: ksubtri                          ! Counter variable for subtriangles loop
    integer           :: nsubtri                          ! Number of subtriangles
    integer           :: subtriangle(8)                   ! Vector that contains what subtriangles need to be integrated
    real(kind=real64) :: theta_subtri(2,8)                ! Matrix that contains the angles theta of the subtriangles to be integrated
    real(kind=real64) :: thetap_subtri(2,8)               ! Matrix that contains the angles thetap of the subtriangles to be integrated
    integer           :: ktheta                           ! Counter variable for theta coordinate loop
    integer           :: krho                             ! Counter variable for rho coordinate loop
    integer           :: kphi                             ! Counter coordinates for shape functions loops
    integer           :: nnodes_g                         ! Number of nodes of geometrical interpolation
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
    real(kind=real64) :: phi_f1(fbem_n_nodes(type_f1))    ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phi_f2(fbem_n_nodes(type_f2))    ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phi_f1_i(fbem_n_nodes(type_f1))  ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phi_g(fbem_n_nodes(type_g))      ! Geometrical shape functions values at xi_1,xi_2
    real(kind=real64) :: dphidxi1_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64) :: dphidxi2_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64) :: x_i(3)                           ! Collocation point position vector
    real(kind=real64) :: x(3)                             ! Position vector at xi_1,xi_2
    real(kind=real64) :: T1(3), T2(3)                     ! Tangent vectors at xi_1,xi_2
    real(kind=real64) :: N(3)                             ! Normal vector at xi_1,xi_2
    real(kind=real64) :: rv(3)                            ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64) :: r, d1r, d1r2                     ! Distance vector module and its inverse
    real(kind=real64) :: drdx(3)                          ! Distance vector derivatives with respect to x_k
    real(kind=real64) :: jg                               ! Geometric jacobian
    real(kind=real64) :: jw                               ! Jacobian * weights
    real(kind=real64) :: drdn                             ! Partial derivative of r respect to unit normal
    real(kind=real64) :: phif1jw(fbem_n_nodes(type_f1))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phif2jw(fbem_n_nodes(type_f2))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: costheta, sintheta               ! cos(theta), sin(theta)
    integer           :: il, ik                           ! Counter for load / observation components
    ! Initialization
    g=0.d0
    ! Number of nodes of geometrical interpolation
    nnodes_g=fbem_n_nodes(type_g)
    ! Checkings
    do kphi=1,nnodes_g
      if (abs(x_nodes(3,kphi)-xp).gt.1.d-12) then
        write(error_unit,*) 'Error: x must belong to the half-space (4)'
        write(*,*) x_nodes
        write(*,*) xi_i
        stop
      end if
    end do
    ! Calculate
    select case (fbem_n_dimension(type_g))
      !
      ! SURFACE LOAD
      !
      case (2)
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
              jg=sqrt(dot_product(N,N))
              ! Unit normal vector
              n=N/jg
              ! Distance vector
              rv=x-x_i
              ! Distance vector norm
              r=dot_product(rv,rv)
              r=sqrt(r)
              d1r=1.d0/r
              ! r_{,k}
              drdx=rv*d1r
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
              ! Loop through load direction and observation direction
              g(:,1,1)=g(:,1,1)+(d1r+(3.-4.*nu)*drdx(1)**2*d1r+4.*(1.-nu)*(1.-2.*nu)/(r+x(3)-xp)*(1.-drdx(1)**2*r/(r+x(3)-xp)))*phif2jw
              g(:,1,2)=g(:,1,2)+((3.-4.*nu)*drdx(1)*drdx(2)*d1r-4.*(1.-nu)*(1.-2.*nu)*drdx(1)*drdx(2)*r/(r+x(3)-xp)**2)*phif2jw
              g(:,1,3)=g(:,1,3)+((3.-4.*nu)*drdx(1)*drdx(3)*d1r+4.*(1.-nu)*(1.-2.*nu)*drdx(1)/(r+x(3)-xp))*phif2jw
              g(:,2,1)=g(:,2,1)+((3.-4.*nu)*drdx(1)*drdx(2)*d1r-4.*(1.-nu)*(1.-2.*nu)*drdx(1)*drdx(2)*r/(r+x(3)-xp)**2)*phif2jw
              g(:,2,2)=g(:,2,2)+(d1r+(3.-4.*nu)*drdx(2)**2*d1r+4.*(1.-nu)*(1.-2.*nu)/(r+x(3)-xp)*(1.-drdx(2)**2*r/(r+x(3)-xp)))*phif2jw
              g(:,2,3)=g(:,2,3)+((3.-4.*nu)*drdx(2)*drdx(3)*d1r+4.*(1.-nu)*(1.-2.*nu)*drdx(2)/(r+x(3)-xp))*phif2jw
              g(:,3,1)=g(:,3,1)+((3.-4.*nu)*drdx(1)*drdx(3)*d1r-4.*(1.-nu)*(1.-2.*nu)*drdx(1)/(r+x(3)-xp))*phif2jw
              g(:,3,2)=g(:,3,2)+((3.-4.*nu)*drdx(2)*drdx(3)*d1r-4.*(1.-nu)*(1.-2.*nu)*drdx(2)/(r+x(3)-xp))*phif2jw
              g(:,3,3)=g(:,3,3)+((8.*(1.-nu)**2-(3.-4.*nu))*d1r+(3.-4.*nu)*drdx(3)**2*d1r)*phif2jw
            end do ! Loop through rho coordinate
          end do ! Loop through theta coordinate
        end do ! Loop through rectangle triangles
        g=g/(16.*c_pi*mu*(1.-nu))
      !
      ! VOLUME LOAD
      !
      case (3)
        stop 'not yet : fbem_bem_staela3d_sbie_bl_int'
      !
      ! OTHERS
      !
      case default
        call fbem_error_message(output_unit,0,'fbem_bem_staela3d_sbie_bl_int',0,'it is only possible to integrate surface or volume loads')
    end select
  end subroutine fbem_bem_staela3d_hsc_sbie_bl_int

  subroutine fbem_bem_staela3d_hsc_sbie_bl_int_cmplx(type_g,type_f1,type_f2,delta_f,x_nodes,xi_i,mu,nu,np,xp,g)
    implicit none
    ! I/O
    integer           :: type_g                          !! Geometrial interpolation
    integer           :: type_f1                         !! Functional interpolation (primary variables)
    integer           :: type_f2                         !! Functional interpolation (secondary variables)
    real(kind=real64) :: delta_f                         !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64) :: x_nodes(3,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
    real(kind=real64) :: xi_i(2)                         !! Reference coordinates of the singular point.
    complex(kind=real64) :: mu                              !! Shear modulus
    complex(kind=real64) :: nu                              !! Poisson's ratio
    integer           :: np                              !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64) :: xp                              !! Coordinate of the plane: x(np)
    complex(kind=real64) :: g(fbem_n_nodes(type_f2),3,3)    !! g kernel vector
    ! Local
    complex(kind=real64) :: hk(fbem_n_nodes(type_f1),3,3)    ! h kernel vector (full-space)
    complex(kind=real64) :: gk(fbem_n_nodes(type_f2),3,3)    ! g kernel vector (full-space)
    integer           :: ksubtri                          ! Counter variable for subtriangles loop
    integer           :: nsubtri                          ! Number of subtriangles
    integer           :: subtriangle(8)                   ! Vector that contains what subtriangles need to be integrated
    real(kind=real64) :: theta_subtri(2,8)                ! Matrix that contains the angles theta of the subtriangles to be integrated
    real(kind=real64) :: thetap_subtri(2,8)               ! Matrix that contains the angles thetap of the subtriangles to be integrated
    integer           :: ktheta                           ! Counter variable for theta coordinate loop
    integer           :: krho                             ! Counter variable for rho coordinate loop
    integer           :: kphi                             ! Counter coordinates for shape functions loops
    integer           :: nnodes_g                         ! Number of nodes of geometrical interpolation
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
    real(kind=real64) :: phi_f1(fbem_n_nodes(type_f1))    ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phi_f2(fbem_n_nodes(type_f2))    ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phi_f1_i(fbem_n_nodes(type_f1))  ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phi_g(fbem_n_nodes(type_g))      ! Geometrical shape functions values at xi_1,xi_2
    real(kind=real64) :: dphidxi1_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64) :: dphidxi2_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64) :: x_i(3)                           ! Collocation point position vector
    real(kind=real64) :: x(3)                             ! Position vector at xi_1,xi_2
    real(kind=real64) :: T1(3), T2(3)                     ! Tangent vectors at xi_1,xi_2
    real(kind=real64) :: N(3)                             ! Normal vector at xi_1,xi_2
    real(kind=real64) :: rv(3)                            ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64) :: r, d1r, d1r2                     ! Distance vector module and its inverse
    real(kind=real64) :: drdx(3)                          ! Distance vector derivatives with respect to x_k
    real(kind=real64) :: jg                               ! Geometric jacobian
    real(kind=real64) :: jw                               ! Jacobian * weights
    real(kind=real64) :: drdn                             ! Partial derivative of r respect to unit normal
    real(kind=real64) :: phif1jw(fbem_n_nodes(type_f1))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: phif2jw(fbem_n_nodes(type_f2))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64) :: costheta, sintheta               ! cos(theta), sin(theta)
    integer           :: il, ik                           ! Counter for load / observation components
    ! Initialization
    g=0.d0
    ! Number of nodes of geometrical interpolation
    nnodes_g=fbem_n_nodes(type_g)
    ! Checkings
    do kphi=1,nnodes_g
      if (abs(x_nodes(3,kphi)-xp).gt.1.d-12) then
        write(error_unit,*) 'Error: x must belong to the half-space (4)'
        write(*,*) x_nodes
        write(*,*) xi_i
        stop
      end if
    end do
    ! Calculate
    select case (fbem_n_dimension(type_g))
      !
      ! SURFACE LOAD
      !
      case (2)
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
              jg=sqrt(dot_product(N,N))
              ! Unit normal vector
              n=N/jg
              ! Distance vector
              rv=x-x_i
              ! Distance vector norm
              r=dot_product(rv,rv)
              r=sqrt(r)
              d1r=1.d0/r
              ! r_{,k}
              drdx=rv*d1r
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
              ! Loop through load direction and observation direction
              g(:,1,1)=g(:,1,1)+(d1r+(3.-4.*nu)*drdx(1)**2*d1r+4.*(1.-nu)*(1.-2.*nu)/(r+x(3)-xp)*(1.-drdx(1)**2*r/(r+x(3)-xp)))*phif2jw
              g(:,1,2)=g(:,1,2)+((3.-4.*nu)*drdx(1)*drdx(2)*d1r-4.*(1.-nu)*(1.-2.*nu)*drdx(1)*drdx(2)*r/(r+x(3)-xp)**2)*phif2jw
              g(:,1,3)=g(:,1,3)+((3.-4.*nu)*drdx(1)*drdx(3)*d1r+4.*(1.-nu)*(1.-2.*nu)*drdx(1)/(r+x(3)-xp))*phif2jw
              g(:,2,1)=g(:,2,1)+((3.-4.*nu)*drdx(1)*drdx(2)*d1r-4.*(1.-nu)*(1.-2.*nu)*drdx(1)*drdx(2)*r/(r+x(3)-xp)**2)*phif2jw
              g(:,2,2)=g(:,2,2)+(d1r+(3.-4.*nu)*drdx(2)**2*d1r+4.*(1.-nu)*(1.-2.*nu)/(r+x(3)-xp)*(1.-drdx(2)**2*r/(r+x(3)-xp)))*phif2jw
              g(:,2,3)=g(:,2,3)+((3.-4.*nu)*drdx(2)*drdx(3)*d1r+4.*(1.-nu)*(1.-2.*nu)*drdx(2)/(r+x(3)-xp))*phif2jw
              g(:,3,1)=g(:,3,1)+((3.-4.*nu)*drdx(1)*drdx(3)*d1r-4.*(1.-nu)*(1.-2.*nu)*drdx(1)/(r+x(3)-xp))*phif2jw
              g(:,3,2)=g(:,3,2)+((3.-4.*nu)*drdx(2)*drdx(3)*d1r-4.*(1.-nu)*(1.-2.*nu)*drdx(2)/(r+x(3)-xp))*phif2jw
              g(:,3,3)=g(:,3,3)+((8.*(1.-nu)**2-(3.-4.*nu))*d1r+(3.-4.*nu)*drdx(3)**2*d1r)*phif2jw
            end do ! Loop through rho coordinate
          end do ! Loop through theta coordinate
        end do ! Loop through rectangle triangles
        g=g/(16.*c_pi*mu*(1.-nu))
      !
      ! VOLUME LOAD
      !
      case (3)
        stop 'not yet : fbem_bem_staela3d_sbie_bl_int'
      !
      ! OTHERS
      !
      case default
        call fbem_error_message(output_unit,0,'fbem_bem_staela3d_sbie_bl_int_cmplx',0,'it is only possible to integrate surface or volume loads')
    end select
  end subroutine fbem_bem_staela3d_hsc_sbie_bl_int_cmplx

  !! Efficient automatic integration
  subroutine fbem_bem_staela3d_hsc_sbie_bl_auto(e,x_i,mu,nu,np,xp,qsp,ns,g)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                !! Integration element
    real(kind=real64)        :: x_i(3)                           !! Collocation point
    real(kind=real64)        :: mu                               !! Shear modulus
    real(kind=real64)        :: nu                               !! Poisson's ratio
    integer                  :: np                               !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64)        :: xp                               !! Coordinate of the plane: x(np)
    type(fbem_qs_parameters) :: qsp                              !! Quasi-singular integration parameters
    integer                  :: ns                               !! Maximum level of subdivisions
    real(kind=real64)        :: g(e%n_snodes,3,3)                !! g integration kernel
    ! Local
    real(kind=real64) :: xp_i(3)                          ! Collocation point image in global axes
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
    ! Collocation point image in global axes
    xp_i=x_i
    xp_i(abs(np))=2.*xp-xp_i(abs(np))
    ! POINT BODY LOAD
    if (e%d.eq.0) then
      r=e%x(:,1)-xp_i
      rmin=sqrt(dot_product(r,r))
      if (rmin.eq.0.d0) then
        call fbem_error_message(output_unit,0,'fbem_bem_staela3d_hsc_sbie_bl_auto',0,'it is not possible to collocate at a point load')
      else
        stop 'not yet : fbem_bem_staela3d_hsc_sbie_bl_auto'
        !call fbem_bem_staela3d_sbie_u(e%x(:,1),x_i,mu,nu,gbl(1,:,:))
        return
      end if
    ! LINE, SURFACE OR VOLUME BODY LOAD
    ! Determine if interior or exterior integration
    !   - Interior integration (delta=1) requires: xi_i
    !   - Exterior integration (delta=0) requires: x_i, barxi, rmin and d
    ! Use the element ball
    else
      r=e%bball_centre-xp_i
      rmin=sqrt(dot_product(r,r))-e%bball_radius
      if (rmin.gt.(4.d0*e%bball_radius)) then
        delta=0
        barxi=0.d0
        d=rmin/e%cl
      else
        ! Use an adaptative algorithm that combines sampling and minimization algorithms
        call fbem_nearest_element_point_bem(3,e%gtype,e%x,e%cl,xp_i,barxi,rmin,d,method)
        if (d.le.1.d-12) then
          delta=1
        else
          delta=0
        end if
      end if
    end if
    ! Integrate
    select case (delta)
      case (1)
        call fbem_bem_staela3d_hsc_sbie_bl_int(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,barxi,mu,nu,np,xp,g)
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
          call fbem_bem_staela3d_hsc_sbie_bl_ext_pre(ps,e,x_i,mu,nu,np,xp,g)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_staela3d_hsc_sbie_bl_ext_adp(e,xi_s,x_i,mu,nu,np,xp,qsp,1,ns,g)
        end if
    end select
  end subroutine fbem_bem_staela3d_hsc_sbie_bl_auto

    !! Efficient automatic integration
  subroutine fbem_bem_staela3d_hsc_sbie_bl_auto_cmplx(e,x_i,mu,nu,np,xp,qsp,ns,g)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                !! Integration element
    real(kind=real64)        :: x_i(3)                           !! Collocation point
    complex(kind=real64)        :: mu                               !! Shear modulus
    complex(kind=real64)        :: nu                               !! Poisson's ratio
    integer                  :: np                               !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64)        :: xp                               !! Coordinate of the plane: x(np)
    type(fbem_qs_parameters) :: qsp                              !! Quasi-singular integration parameters
    integer                  :: ns                               !! Maximum level of subdivisions
    complex(kind=real64)        :: g(e%n_snodes,3,3)                !! g integration kernel
    ! Local
    real(kind=real64) :: xp_i(3)                          ! Collocation point image in global axes
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
    ! Collocation point image in global axes
    xp_i=x_i
    xp_i(abs(np))=2.*xp-xp_i(abs(np))
    ! POINT BODY LOAD
    if (e%d.eq.0) then
      r=e%x(:,1)-xp_i
      rmin=sqrt(dot_product(r,r))
      if (rmin.eq.0.d0) then
        call fbem_error_message(output_unit,0,'fbem_bem_staela3d_hsc_sbie_bl_auto_cmplx',0,'it is not possible to collocate at a point load')
      else
        stop 'not yet : fbem_bem_staela3d_hsc_sbie_bl_auto_cmplx'
        !call fbem_bem_staela3d_sbie_u(e%x(:,1),x_i,mu,nu,gbl(1,:,:))
        return
      end if
    ! LINE, SURFACE OR VOLUME BODY LOAD
    ! Determine if interior or exterior integration
    !   - Interior integration (delta=1) requires: xi_i
    !   - Exterior integration (delta=0) requires: x_i, barxi, rmin and d
    ! Use the element ball
    else
      r=e%bball_centre-xp_i
      rmin=sqrt(dot_product(r,r))-e%bball_radius
      if (rmin.gt.(4.d0*e%bball_radius)) then
        delta=0
        barxi=0.d0
        d=rmin/e%cl
      else
        ! Use an adaptative algorithm that combines sampling and minimization algorithms
        call fbem_nearest_element_point_bem(3,e%gtype,e%x,e%cl,xp_i,barxi,rmin,d,method)
        if (d.le.1.d-12) then
          delta=1
        else
          delta=0
        end if
      end if
    end if
    ! Integrate
    select case (delta)
      case (1)
        call fbem_bem_staela3d_hsc_sbie_bl_int_cmplx(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,barxi,mu,nu,np,xp,g)
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
          call fbem_bem_staela3d_hsc_sbie_bl_ext_pre_cmplx(ps,e,x_i,mu,nu,np,xp,g)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_staela3d_hsc_sbie_bl_ext_adp_cmplx(e,xi_s,x_i,mu,nu,np,xp,qsp,1,ns,g)
        end if
    end select
  end subroutine fbem_bem_staela3d_hsc_sbie_bl_auto_cmplx

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (HBIE) - TRACTION-FREE HALF-SPACE FUNDAMENTAL SOLUTION (COMPLEMENTARY PART)
  ! --------------------------------------------------------------------------------------------------------------------------------

  subroutine fbem_bem_staela3d_hsc_hbie(mu,nu,np,xp,x_i,n_i,x,n,d,s)
    implicit none
    ! I/O
    real(kind=real64), intent( in) :: mu      !! Shear modulus
    real(kind=real64), intent( in) :: nu      !! Poisson's ratio
    integer          , intent( in) :: np      !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64), intent( in) :: xp      !! Coordinate of the plane: x(np)
    real(kind=real64), intent( in) :: x_i(3)  !! Collocation point (original, not the image)
    real(kind=real64), intent( in) :: n_i(3)  !! Collocation point unit normal
    real(kind=real64), intent( in) :: x(3)
    real(kind=real64), intent( in) :: n(3)
    real(kind=real64), intent(out) :: d(3,3)
    real(kind=real64), intent(out) :: s(3,3)
    ! Local
    real(kind=real64)              :: c1, c2, c3, c4, c5
    real(kind=real64)              :: rx, ry, z, zi, R2p1, R2p2, R2p3, R2p4, R2p5, R2p6, R2p7, R2p9
    real(kind=real64)              :: Ku, lambda
    integer                        :: il, ik, ij, ip
    real(kind=real64)              :: dudxi(3,3,3)     ! d u*_{lk} / d x^i_p
    real(kind=real64)              :: dudxdxi(3,3,3,3) ! d u*_{lk} / d x_j / d x^i_p
    real(kind=real64)              :: dsdxi(3,3,3,3)   ! d sigma*_{lkj} / d x^i_p
    real(kind=real64)              :: dtdxi(3,3,3)     ! d t*_{lk} / d x^i_p
    real(kind=real64)              :: dlkp(3,3,3)      ! d*_{lkp}
    real(kind=real64)              :: slkp(3,3,3)      ! s*_{lkp}
    ! Checkings
    if (np.ne.-3) then
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,'only the case of a half-space with unit normal n=(0,0,-1) is allowed')
    end if
    if (x_i(3).lt.xp) then
      write(error_unit,*) 'Error: x_i must belong to the half-space'
      stop
    end if
    if (x(3).lt.xp) then
      write(error_unit,*) 'Error: x must belong to the half-space (5)'
      stop
    end if
    !
    ! u*_{lk} derivatives required for hypersingular formulation. Obtained using Mathematica. Without 1/(16*pi*mu*(1-nu)) factor
    !
    ! Parameters
    c1=1.-nu
    c2=1.-2.*nu
    c3=3.-2.*nu
    c4=3.-4.*nu
    c5=5.-4.*nu
    rx=x(1)-x_i(1)
    ry=x(2)-x_i(2)
    z=x(3)-xp
    zi=x_i(3)-xp
    R2p1=sqrt(rx**2+ry**2+(z+zi)**2)
    R2p2=R2p1**2
    R2p3=R2p2*R2p1
    R2p4=R2p3*R2p1
    R2p5=R2p4*R2p1
    R2p6=R2p5*R2p1
    R2p7=R2p6*R2p1
    R2p9=R2p7*R2p2
    ! d u*_{lk} / d x^i_p
    dudxi(1,1,1)=(rx*(R2p2*(R2p2-2*c4*R2p2+3*c4*rx**2)+6*(3*R2p2-5*rx**2)*z*zi-(8*c1*c2*R2p5*rx**2)/(R2p1+z+zi)**3+(4*c1*c2*R2p4*(3*R2p2-rx**2))/(R2p1+z+zi)**2))/R2p7
    dudxi(2,1,1)=(ry*(-(c4*R2p2*(R2p2-3*rx**2))+6*(R2p2-5*rx**2)*z*zi-(8*c1*c2*R2p5*rx**2)/(R2p1+z+zi)**3+(4*c1*c2*R2p4*(R2p2-rx**2))/(R2p1+z+zi)**2))/R2p7
    dudxi(3,1,1)=((-(c4*R2p2*(R2p2-3*rx**2)*(z-zi))-6*(R2p2-5*rx**2)*z**2*zi-6*(R2p2-5*rx**2)*z*zi**2-(4*c1*c2*R2p5*rx**2)/(R2p1+z+zi)**2+(4*c1*c2*R2p4*(R2p2-rx**2))/(R2p1+z+zi)))/R2p7
    dudxi(1,2,1)=(ry*(-(c4*R2p2*(R2p2-3*rx**2))+6*(R2p2-5*rx**2)*z*zi-(8*c1*c2*R2p5*rx**2)/(R2p1+z+zi)**3+(4*c1*c2*R2p4*(R2p2-rx**2))/(R2p1+z+zi)**2))/R2p7
    dudxi(2,2,1)=(rx*(R2p2*(R2p2+3*c4*ry**2)+6*(R2p2-5*ry**2)*z*zi-(8*c1*c2*R2p5*ry**2)/(R2p1+z+zi)**3+(4*c1*c2*R2p4*(R2p2-ry**2))/(R2p1+z+zi)**2))/R2p7
    dudxi(3,2,1)=(rx*ry*(3*c4*R2p2*(z-zi)+30*z**2*zi+30*z*zi**2-(4*c1*c2*R2p4*(2*R2p1+z+zi))/(R2p1+z+zi)**2))/R2p7
    dudxi(1,3,1)=((-(c4*R2p2*(R2p2-3*rx**2)*(z-zi))+6*(R2p2-5*rx**2)*z**2*zi+6*(R2p2-5*rx**2)*z*zi**2+(4*c1*c2*R2p5*rx**2)/(R2p1+z+zi)**2+(4*c1*c2*R2p4*(-R2p2+rx**2))/(R2p1+z+zi)))/R2p7
    dudxi(2,3,1)=(rx*ry*(3*c4*R2p2*(z-zi)-30*z**2*zi-30*z*zi**2+(4*c1*c2*R2p4*(2*R2p1+z+zi))/(R2p1+z+zi)**2))/R2p7
    dudxi(3,3,1)=(rx*(8*c1**2*R2p4+c4*R2p2*(-R2p2+3*(z+zi)**2)+6*z*zi*(-R2p2+5*(z+zi)**2)))/R2p7
    dudxi(1,1,2)=(ry*(R2p2*(R2p2+3*c4*rx**2)+6*(R2p2-5*rx**2)*z*zi-(8*c1*c2*R2p5*rx**2)/(R2p1+z+zi)**3+(4*c1*c2*R2p4*(R2p2-rx**2))/(R2p1+z+zi)**2))/R2p7
    dudxi(2,1,2)=(rx*(-(c4*R2p2*(R2p2-3*ry**2))+6*(R2p2-5*ry**2)*z*zi-(8*c1*c2*R2p5*ry**2)/(R2p1+z+zi)**3+(4*c1*c2*R2p4*(R2p2-ry**2))/(R2p1+z+zi)**2))/R2p7
    dudxi(3,1,2)=(rx*ry*(3*c4*R2p2*(z-zi)+30*z**2*zi+30*z*zi**2-(4*c1*c2*R2p4*(2*R2p1+z+zi))/(R2p1+z+zi)**2))/R2p7
    dudxi(1,2,2)=(rx*(-(c4*R2p2*(R2p2-3*ry**2))+6*(R2p2-5*ry**2)*z*zi-(8*c1*c2*R2p5*ry**2)/(R2p1+z+zi)**3+(4*c1*c2*R2p4*(R2p2-ry**2))/(R2p1+z+zi)**2))/R2p7
    dudxi(2,2,2)=(ry*(R2p2*(R2p2-2*c4*R2p2+3*c4*ry**2)+6*(3*R2p2-5*ry**2)*z*zi-(8*c1*c2*R2p5*ry**2)/(R2p1+z+zi)**3+(4*c1*c2*R2p4*(3*R2p2-ry**2))/(R2p1+z+zi)**2))/R2p7
    dudxi(3,2,2)=((-(c4*R2p2*(R2p2-3*ry**2)*(z-zi))-6*(R2p2-5*ry**2)*z**2*zi-6*(R2p2-5*ry**2)*z*zi**2-(4*c1*c2*R2p5*ry**2)/(R2p1+z+zi)**2+(4*c1*c2*R2p4*(R2p2-ry**2))/(R2p1+z+zi)))/R2p7
    dudxi(1,3,2)=(rx*ry*(3*c4*R2p2*(z-zi)-30*z**2*zi-30*z*zi**2+(4*c1*c2*R2p4*(2*R2p1+z+zi))/(R2p1+z+zi)**2))/R2p7
    dudxi(2,3,2)=((-(c4*R2p2*(R2p2-3*ry**2)*(z-zi))+6*(R2p2-5*ry**2)*z**2*zi+6*(R2p2-5*ry**2)*z*zi**2+(4*c1*c2*R2p5*ry**2)/(R2p1+z+zi)**2+(4*c1*c2*R2p4*(-R2p2+ry**2))/(R2p1+z+zi)))/R2p7
    dudxi(3,3,2)=(ry*(8*c1**2*R2p4+c4*R2p2*(-R2p2+3*(z+zi)**2)+6*z*zi*(-R2p2+5*(z+zi)**2)))/R2p7
    dudxi(1,1,3)=((30*rx**2*z*zi*(z+zi)+(4*c1*c2*R2p5*rx**2)/(R2p1+z+zi)**2-(4*c1*c2*R2p6)/(R2p1+z+zi)+R2p4*(z-zi+(4*c1*c2*rx**2)/(R2p1+z+zi))-3*R2p2*(2*z*zi*(z+zi)+rx**2*((2+c4)*z+c4*zi))))/R2p7
    dudxi(2,1,3)=(rx*ry*(30*z*zi*(z+zi)+(4*c1*c2*R2p5)/(R2p1+z+zi)**2+(4*c1*c2*R2p4)/(R2p1+z+zi)-3*R2p2*((2+c4)*z+c4*zi)))/R2p7
    dudxi(3,1,3)=(rx*(4*c1*c2*R2p4-c4*R2p2*(R2p2+3*(z-zi)*(z+zi))+6*z*(-5*zi*(z+zi)**2+R2p2*(z+2*zi))))/R2p7
    dudxi(1,2,3)=(rx*ry*(30*z*zi*(z+zi)+(4*c1*c2*R2p5)/(R2p1+z+zi)**2+(4*c1*c2*R2p4)/(R2p1+z+zi)-3*R2p2*((2+c4)*z+c4*zi)))/R2p7
    dudxi(2,2,3)=((30*ry**2*z*zi*(z+zi)+(4*c1*c2*R2p5*ry**2)/(R2p1+z+zi)**2-(4*c1*c2*R2p6)/(R2p1+z+zi)+R2p4*(z-zi+(4*c1*c2*ry**2)/(R2p1+z+zi))-3*R2p2*(2*z*zi*(z+zi)+ry**2*((2+c4)*z+c4*zi))))/R2p7
    dudxi(3,2,3)=(ry*(4*c1*c2*R2p4-c4*R2p2*(R2p2+3*(z-zi)*(z+zi))+6*z*(-5*zi*(z+zi)**2+R2p2*(z+2*zi))))/R2p7
    dudxi(1,3,3)=-((rx*(4*c1*c2*R2p4+c4*R2p2*(R2p2+3*(z-zi)*(z+zi))+6*z*(-5*zi*(z+zi)**2+R2p2*(z+2*zi))))/R2p7)
    dudxi(2,3,3)=-((ry*(4*c1*c2*R2p4+c4*R2p2*(R2p2+3*(z-zi)*(z+zi))+6*z*(-5*zi*(z+zi)**2+R2p2*(z+2*zi))))/R2p7)
    dudxi(3,3,3)=-(((30*z*zi*(z+zi)**3+R2p4*((2+8*c1**2-3*c4)*z+(8*c1**2-3*c4)*zi)+3*R2p2*(z+zi)*((-2+c4)*z**2+2*(-4+c4)*z*zi+c4*zi**2)))/R2p7)
    ! d u*_{lk} / d x_j / d x^i_1
    dudxdxi(1,1,1,1)=((R2p2*((1-2*c4)*R2p4+3*(-1+5*c4)*R2p2*rx**2-15*c4*rx**4)+6*(3*R2p4-30*R2p2*rx**2+35*rx**4)*z*zi+(24*c1*c2*R2p6*rx**4)/(R2p1+z+zi)**4+(24*c1*c2*R2p5*rx**2*(-2*R2p2+rx**2))/(R2p1+z+zi)**3+(12*c1*c2*R2p4*(R2p2-rx**2)**2)/(R2p1+z+zi)**2))/R2p9
    dudxdxi(2,1,1,1)=(3*rx*ry*(c4*R2p2*(3*R2p2-5*rx**2)+10*(-3*R2p2+7*rx**2)*z*zi+(8*c1*c2*R2p6*rx**2)/(R2p1+z+zi)**4+(8*c1*c2*R2p5*(-R2p2+rx**2))/(R2p1+z+zi)**3+(4*c1*c2*R2p4*(-R2p2+rx**2))/(R2p1+z+zi)**2))/R2p9
    dudxdxi(3,1,1,1)=(rx*(3*c4*R2p2*(3*R2p2-5*rx**2)*(z-zi)+30*(3*R2p2-7*rx**2)*z**2*zi+30*(3*R2p2-7*rx**2)*z*zi**2+(8*c1*c2*R2p6*rx**2)/(R2p1+z+zi)**3+(12*c1*c2*R2p5*(-R2p2+rx**2))/(R2p1+z+zi)**2+(12*c1*c2*R2p4*(-R2p2+rx**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(1,2,1,1)=(3*rx*ry*(c4*R2p2*(3*R2p2-5*rx**2)+10*(-3*R2p2+7*rx**2)*z*zi+(8*c1*c2*R2p6*rx**2)/(R2p1+z+zi)**4+(8*c1*c2*R2p5*(-R2p2+rx**2))/(R2p1+z+zi)**3+(4*c1*c2*R2p4*(-R2p2+rx**2))/(R2p1+z+zi)**2))/R2p9
    dudxdxi(2,2,1,1)=((R2p2*(R2p2*(R2p2-3*rx**2)+3*c4*(R2p2-5*rx**2)*ry**2)+6*(R2p4+35*rx**2*ry**2-5*R2p2*(rx**2+ry**2))*z*zi+(24*c1*c2*R2p6*rx**2*ry**2)/(R2p1+z+zi)**4-(8*c1*c2*R2p5*(R2p2*rx**2+(R2p2-3*rx**2)*ry**2))/(R2p1+z+zi)**3+(4*c1*c2*R2p4*(R2p4+3*rx**2*ry**2-R2p2*(rx**2+ry**2)))/(R2p1+z+zi)**2))/R2p9
    dudxdxi(3,2,1,1)=(ry*(3*c4*R2p2*(R2p2-5*rx**2)*(z-zi)+30*(R2p2-7*rx**2)*z**2*zi+30*(R2p2-7*rx**2)*z*zi**2+(8*c1*c2*R2p6*rx**2)/(R2p1+z+zi)**3-(4*c1*c2*R2p5*(R2p2-3*rx**2))/(R2p1+z+zi)**2-(4*c1*c2*R2p4*(R2p2-3*rx**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(1,3,1,1)=(rx*(3*c4*R2p2*(3*R2p2-5*rx**2)*(z-zi)+30*(-3*R2p2+7*rx**2)*z**2*zi+30*(-3*R2p2+7*rx**2)*z*zi**2-(8*c1*c2*R2p6*rx**2)/(R2p1+z+zi)**3+(12*c1*c2*R2p5*(R2p2-rx**2))/(R2p1+z+zi)**2+(12*c1*c2*R2p4*(R2p2-rx**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(2,3,1,1)=(ry*(3*c4*R2p2*(R2p2-5*rx**2)*(z-zi)-30*(R2p2-7*rx**2)*z**2*zi-30*(R2p2-7*rx**2)*z*zi**2-(8*c1*c2*R2p6*rx**2)/(R2p1+z+zi)**3+(4*c1*c2*R2p5*(R2p2-3*rx**2))/(R2p1+z+zi)**2+(4*c1*c2*R2p4*(R2p2-3*rx**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(3,3,1,1)=((8*c1**2*R2p4*(R2p2-3*rx**2)-6*z*zi*(R2p4+35*rx**2*(z+zi)**2-5*R2p2*(rx**2+(z+zi)**2))-c4*R2p2*(R2p4+15*rx**2*(z+zi)**2-3*R2p2*(rx**2+(z+zi)**2))))/R2p9
    dudxdxi(1,1,2,1)=(3*rx*ry*(R2p2*((-1+2*c4)*R2p2-5*c4*rx**2)+10*(-3*R2p2+7*rx**2)*z*zi+(8*c1*c2*R2p6*rx**2)/(R2p1+z+zi)**4+(8*c1*c2*R2p5*(-R2p2+rx**2))/(R2p1+z+zi)**3+(4*c1*c2*R2p4*(-R2p2+rx**2))/(R2p1+z+zi)**2))/R2p9
    dudxdxi(2,1,2,1)=((c4*R2p2*(-(R2p2*(R2p2-3*rx**2))+3*(R2p2-5*rx**2)*ry**2)+6*(R2p4+35*rx**2*ry**2-5*R2p2*(rx**2+ry**2))*z*zi+(24*c1*c2*R2p6*rx**2*ry**2)/(R2p1+z+zi)**4-(8*c1*c2*R2p5*(R2p2*rx**2+(R2p2-3*rx**2)*ry**2))/(R2p1+z+zi)**3+(4*c1*c2*R2p4*(R2p4+3*rx**2*ry**2-R2p2*(rx**2+ry**2)))/(R2p1+z+zi)**2))/R2p9
    dudxdxi(3,1,2,1)=(ry*(3*c4*R2p2*(R2p2-5*rx**2)*(z-zi)+30*(R2p2-7*rx**2)*z**2*zi+30*(R2p2-7*rx**2)*z*zi**2+(8*c1*c2*R2p6*rx**2)/(R2p1+z+zi)**3-(4*c1*c2*R2p5*(R2p2-3*rx**2))/(R2p1+z+zi)**2-(4*c1*c2*R2p4*(R2p2-3*rx**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(1,2,2,1)=((c4*R2p2*(-(R2p2*(R2p2-3*rx**2))+3*(R2p2-5*rx**2)*ry**2)+6*(R2p4+35*rx**2*ry**2-5*R2p2*(rx**2+ry**2))*z*zi+(24*c1*c2*R2p6*rx**2*ry**2)/(R2p1+z+zi)**4-(8*c1*c2*R2p5*(R2p2*rx**2+(R2p2-3*rx**2)*ry**2))/(R2p1+z+zi)**3+(4*c1*c2*R2p4*(R2p4+3*rx**2*ry**2-R2p2*(rx**2+ry**2)))/(R2p1+z+zi)**2))/R2p9
    dudxdxi(2,2,2,1)=(3*rx*ry*(R2p2*((-1+2*c4)*R2p2-5*c4*ry**2)+10*(-3*R2p2+7*ry**2)*z*zi+(8*c1*c2*R2p6*ry**2)/(R2p1+z+zi)**4+(8*c1*c2*R2p5*(-R2p2+ry**2))/(R2p1+z+zi)**3+(4*c1*c2*R2p4*(-R2p2+ry**2))/(R2p1+z+zi)**2))/R2p9
    dudxdxi(3,2,2,1)=(rx*(3*c4*R2p2*(R2p2-5*ry**2)*(z-zi)+30*(R2p2-7*ry**2)*z**2*zi+30*(R2p2-7*ry**2)*z*zi**2+(8*c1*c2*R2p6*ry**2)/(R2p1+z+zi)**3-(4*c1*c2*R2p5*(R2p2-3*ry**2))/(R2p1+z+zi)**2-(4*c1*c2*R2p4*(R2p2-3*ry**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(1,3,2,1)=(ry*(3*c4*R2p2*(R2p2-5*rx**2)*(z-zi)-30*(R2p2-7*rx**2)*z**2*zi-30*(R2p2-7*rx**2)*z*zi**2-(8*c1*c2*R2p6*rx**2)/(R2p1+z+zi)**3+(4*c1*c2*R2p5*(R2p2-3*rx**2))/(R2p1+z+zi)**2+(4*c1*c2*R2p4*(R2p2-3*rx**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(2,3,2,1)=(rx*(3*c4*R2p2*(R2p2-5*ry**2)*(z-zi)-30*(R2p2-7*ry**2)*z**2*zi-30*(R2p2-7*ry**2)*z*zi**2-(8*c1*c2*R2p6*ry**2)/(R2p1+z+zi)**3+(4*c1*c2*R2p5*(R2p2-3*ry**2))/(R2p1+z+zi)**2+(4*c1*c2*R2p4*(R2p2-3*ry**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(3,3,2,1)=(-3*rx*ry*(8*c1**2*R2p4+c4*R2p2*(-R2p2+5*(z+zi)**2)+10*z*zi*(-R2p2+7*(z+zi)**2)))/R2p9
    dudxdxi(1,1,3,1)=-((rx*(4*c1*c2*R2p4*(6*R2p4+9*R2p3*(z+zi)-9*R2p1*rx**2*(z+zi)-3*rx**2*(z+zi)**2+R2p2*(-8*rx**2+3*(z+zi)**2))-3*(R2p1+z+zi)**3*(70*rx**2*z*zi*(z+zi)+R2p4*((-1+2*c4)*z+(5+2*c4)*zi)-5*R2p2*(c4*rx**2*(z+zi)+2*zi*(rx**2+3*z*(z+zi))))))/(R2p9*(R2p1+z+zi)**3))
    dudxdxi(2,1,3,1)=(ry*(6*(R2p4+35*rx**2*z**2-5*R2p2*(rx**2+z**2))*zi-30*(R2p2-7*rx**2)*z*zi**2+3*c4*R2p2*(R2p2-5*rx**2)*(z+zi)+(8*c1*c2*R2p6*rx**2)/(R2p1+z+zi)**3-(4*c1*c2*R2p5*(R2p2-3*rx**2))/(R2p1+z+zi)**2-(4*c1*c2*R2p4*(R2p2-3*rx**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(3,1,3,1)=-(((4*c1*c2*R2p4*(R2p2-3*rx**2)+c4*R2p2*(R2p4+15*rx**2*(z-zi)*(z+zi)-3*R2p2*(rx**2+z**2-zi**2))+6*zi*(35*rx**2*z*(z+zi)**2+R2p4*(2*z+zi)-5*R2p2*(z*(z+zi)**2+rx**2*(2*z+zi)))))/R2p9)
    dudxdxi(1,2,3,1)=(ry*(6*(R2p4+35*rx**2*z**2-5*R2p2*(rx**2+z**2))*zi-30*(R2p2-7*rx**2)*z*zi**2+3*c4*R2p2*(R2p2-5*rx**2)*(z+zi)+(8*c1*c2*R2p6*rx**2)/(R2p1+z+zi)**3-(4*c1*c2*R2p5*(R2p2-3*rx**2))/(R2p1+z+zi)**2-(4*c1*c2*R2p4*(R2p2-3*rx**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(2,2,3,1)=-((rx*(4*c1*c2*R2p4*(2*R2p4+3*R2p3*(z+zi)-9*R2p1*ry**2*(z+zi)-3*ry**2*(z+zi)**2+R2p2*(-8*ry**2+(z+zi)**2))+3*(R2p1+z+zi)**3*(R2p4*(z-zi)-70*ry**2*z*zi*(z+zi)+5*R2p2*(c4*ry**2*(z+zi)+2*zi*(ry**2+z*(z+zi))))))/(R2p9*(R2p1+z+zi)**3))
    dudxdxi(3,2,3,1)=(3*rx*ry*(4*c1*c2*R2p4+c4*R2p2*(R2p2-5*z**2+5*zi**2)+10*zi*(-7*z*(z+zi)**2+R2p2*(2*z+zi))))/R2p9
    dudxdxi(1,3,3,1)=((4*c1*c2*R2p4*(R2p2-3*rx**2)-c4*R2p2*(R2p4+15*rx**2*(z-zi)*(z+zi)-3*R2p2*(rx**2+z**2-zi**2))+6*zi*(35*rx**2*z*(z+zi)**2+R2p4*(2*z+zi)-5*R2p2*(z*(z+zi)**2+rx**2*(2*z+zi)))))/R2p9
    dudxdxi(2,3,3,1)=(-3*rx*ry*(4*c1*c2*R2p4-c4*R2p2*(R2p2-5*z**2+5*zi**2)+10*zi*(-7*z*(z+zi)**2+R2p2*(2*z+zi))))/R2p9
    dudxdxi(3,3,3,1)=(3*rx*(-8*c1**2*R2p4*(z+zi)+c4*R2p2*(z+zi)*(3*R2p2-5*(z+zi)**2)-2*zi*(R2p4+35*z*(z+zi)**3-5*R2p2*(z+zi)*(4*z+zi))))/R2p9
    ! d u*_{lk} / d x_j / d x^i_2
    dudxdxi(1,1,1,2)=(3*rx*ry*(R2p2*((-1+2*c4)*R2p2-5*c4*rx**2)+10*(-3*R2p2+7*rx**2)*z*zi+(8*c1*c2*R2p6*rx**2)/(R2p1+z+zi)**4+(8*c1*c2*R2p5*(-R2p2+rx**2))/(R2p1+z+zi)**3+(4*c1*c2*R2p4*(-R2p2+rx**2))/(R2p1+z+zi)**2))/R2p9
    dudxdxi(2,1,1,2)=((c4*R2p2*(-(R2p2*(R2p2-3*rx**2))+3*(R2p2-5*rx**2)*ry**2)+6*(R2p4+35*rx**2*ry**2-5*R2p2*(rx**2+ry**2))*z*zi+(24*c1*c2*R2p6*rx**2*ry**2)/(R2p1+z+zi)**4-(8*c1*c2*R2p5*(R2p2*rx**2+(R2p2-3*rx**2)*ry**2))/(R2p1+z+zi)**3+(4*c1*c2*R2p4*(R2p4+3*rx**2*ry**2-R2p2*(rx**2+ry**2)))/(R2p1+z+zi)**2))/R2p9
    dudxdxi(3,1,1,2)=(ry*(3*c4*R2p2*(R2p2-5*rx**2)*(z-zi)+30*(R2p2-7*rx**2)*z**2*zi+30*(R2p2-7*rx**2)*z*zi**2+(8*c1*c2*R2p6*rx**2)/(R2p1+z+zi)**3-(4*c1*c2*R2p5*(R2p2-3*rx**2))/(R2p1+z+zi)**2-(4*c1*c2*R2p4*(R2p2-3*rx**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(1,2,1,2)=((c4*R2p2*(-(R2p2*(R2p2-3*rx**2))+3*(R2p2-5*rx**2)*ry**2)+6*(R2p4+35*rx**2*ry**2-5*R2p2*(rx**2+ry**2))*z*zi+(24*c1*c2*R2p6*rx**2*ry**2)/(R2p1+z+zi)**4-(8*c1*c2*R2p5*(R2p2*rx**2+(R2p2-3*rx**2)*ry**2))/(R2p1+z+zi)**3+(4*c1*c2*R2p4*(R2p4+3*rx**2*ry**2-R2p2*(rx**2+ry**2)))/(R2p1+z+zi)**2))/R2p9
    dudxdxi(2,2,1,2)=(3*rx*ry*(R2p2*((-1+2*c4)*R2p2-5*c4*ry**2)+10*(-3*R2p2+7*ry**2)*z*zi+(8*c1*c2*R2p6*ry**2)/(R2p1+z+zi)**4+(8*c1*c2*R2p5*(-R2p2+ry**2))/(R2p1+z+zi)**3+(4*c1*c2*R2p4*(-R2p2+ry**2))/(R2p1+z+zi)**2))/R2p9
    dudxdxi(3,2,1,2)=(rx*(3*c4*R2p2*(R2p2-5*ry**2)*(z-zi)+30*(R2p2-7*ry**2)*z**2*zi+30*(R2p2-7*ry**2)*z*zi**2+(8*c1*c2*R2p6*ry**2)/(R2p1+z+zi)**3-(4*c1*c2*R2p5*(R2p2-3*ry**2))/(R2p1+z+zi)**2-(4*c1*c2*R2p4*(R2p2-3*ry**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(1,3,1,2)=(ry*(3*c4*R2p2*(R2p2-5*rx**2)*(z-zi)-30*(R2p2-7*rx**2)*z**2*zi-30*(R2p2-7*rx**2)*z*zi**2-(8*c1*c2*R2p6*rx**2)/(R2p1+z+zi)**3+(4*c1*c2*R2p5*(R2p2-3*rx**2))/(R2p1+z+zi)**2+(4*c1*c2*R2p4*(R2p2-3*rx**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(2,3,1,2)=(rx*(3*c4*R2p2*(R2p2-5*ry**2)*(z-zi)-30*(R2p2-7*ry**2)*z**2*zi-30*(R2p2-7*ry**2)*z*zi**2-(8*c1*c2*R2p6*ry**2)/(R2p1+z+zi)**3+(4*c1*c2*R2p5*(R2p2-3*ry**2))/(R2p1+z+zi)**2+(4*c1*c2*R2p4*(R2p2-3*ry**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(3,3,1,2)=(-3*rx*ry*(8*c1**2*R2p4+c4*R2p2*(-R2p2+5*(z+zi)**2)+10*z*zi*(-R2p2+7*(z+zi)**2)))/R2p9
    dudxdxi(1,1,2,2)=((R2p2*(R2p2*(R2p2+3*c4*rx**2)-3*(R2p2+5*c4*rx**2)*ry**2)+6*(R2p4+35*rx**2*ry**2-5*R2p2*(rx**2+ry**2))*z*zi+(24*c1*c2*R2p6*rx**2*ry**2)/(R2p1+z+zi)**4-(8*c1*c2*R2p5*(R2p2*rx**2+(R2p2-3*rx**2)*ry**2))/(R2p1+z+zi)**3+(4*c1*c2*R2p4*(R2p4+3*rx**2*ry**2-R2p2*(rx**2+ry**2)))/(R2p1+z+zi)**2))/R2p9
    dudxdxi(2,1,2,2)=(3*rx*ry*(c4*R2p2*(3*R2p2-5*ry**2)+10*(-3*R2p2+7*ry**2)*z*zi+(8*c1*c2*R2p6*ry**2)/(R2p1+z+zi)**4+(8*c1*c2*R2p5*(-R2p2+ry**2))/(R2p1+z+zi)**3+(4*c1*c2*R2p4*(-R2p2+ry**2))/(R2p1+z+zi)**2))/R2p9
    dudxdxi(3,1,2,2)=(rx*(3*c4*R2p2*(R2p2-5*ry**2)*(z-zi)+30*(R2p2-7*ry**2)*z**2*zi+30*(R2p2-7*ry**2)*z*zi**2+(8*c1*c2*R2p6*ry**2)/(R2p1+z+zi)**3-(4*c1*c2*R2p5*(R2p2-3*ry**2))/(R2p1+z+zi)**2-(4*c1*c2*R2p4*(R2p2-3*ry**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(1,2,2,2)=(3*rx*ry*(c4*R2p2*(3*R2p2-5*ry**2)+10*(-3*R2p2+7*ry**2)*z*zi+(8*c1*c2*R2p6*ry**2)/(R2p1+z+zi)**4+(8*c1*c2*R2p5*(-R2p2+ry**2))/(R2p1+z+zi)**3+(4*c1*c2*R2p4*(-R2p2+ry**2))/(R2p1+z+zi)**2))/R2p9
    dudxdxi(2,2,2,2)=((R2p2*((1-2*c4)*R2p4+3*(-1+5*c4)*R2p2*ry**2-15*c4*ry**4)+6*(3*R2p4-30*R2p2*ry**2+35*ry**4)*z*zi+(24*c1*c2*R2p6*ry**4)/(R2p1+z+zi)**4+(24*c1*c2*R2p5*ry**2*(-2*R2p2+ry**2))/(R2p1+z+zi)**3+(12*c1*c2*R2p4*(R2p2-ry**2)**2)/(R2p1+z+zi)**2))/R2p9
    dudxdxi(3,2,2,2)=(ry*(3*c4*R2p2*(3*R2p2-5*ry**2)*(z-zi)+30*(3*R2p2-7*ry**2)*z**2*zi+30*(3*R2p2-7*ry**2)*z*zi**2+(8*c1*c2*R2p6*ry**2)/(R2p1+z+zi)**3+(12*c1*c2*R2p5*(-R2p2+ry**2))/(R2p1+z+zi)**2+(12*c1*c2*R2p4*(-R2p2+ry**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(1,3,2,2)=(rx*(3*c4*R2p2*(R2p2-5*ry**2)*(z-zi)-30*(R2p2-7*ry**2)*z**2*zi-30*(R2p2-7*ry**2)*z*zi**2-(8*c1*c2*R2p6*ry**2)/(R2p1+z+zi)**3+(4*c1*c2*R2p5*(R2p2-3*ry**2))/(R2p1+z+zi)**2+(4*c1*c2*R2p4*(R2p2-3*ry**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(2,3,2,2)=(ry*(3*c4*R2p2*(3*R2p2-5*ry**2)*(z-zi)+30*(-3*R2p2+7*ry**2)*z**2*zi+30*(-3*R2p2+7*ry**2)*z*zi**2-(8*c1*c2*R2p6*ry**2)/(R2p1+z+zi)**3+(12*c1*c2*R2p5*(R2p2-ry**2))/(R2p1+z+zi)**2+(12*c1*c2*R2p4*(R2p2-ry**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(3,3,2,2)=((8*c1**2*R2p4*(R2p2-3*ry**2)-6*z*zi*(R2p4+35*ry**2*(z+zi)**2-5*R2p2*(ry**2+(z+zi)**2))-c4*R2p2*(R2p4+15*ry**2*(z+zi)**2-3*R2p2*(ry**2+(z+zi)**2))))/R2p9
    dudxdxi(1,1,3,2)=-((ry*(4*c1*c2*R2p4*(2*R2p4+3*R2p3*(z+zi)-9*R2p1*rx**2*(z+zi)-3*rx**2*(z+zi)**2+R2p2*(-8*rx**2+(z+zi)**2))+3*(R2p1+z+zi)**3*(R2p4*(z-zi)-70*rx**2*z*zi*(z+zi)+5*R2p2*(c4*rx**2*(z+zi)+2*zi*(rx**2+z*(z+zi))))))/(R2p9*(R2p1+z+zi)**3))
    dudxdxi(2,1,3,2)=(rx*(6*(R2p4+35*ry**2*z**2-5*R2p2*(ry**2+z**2))*zi-30*(R2p2-7*ry**2)*z*zi**2+3*c4*R2p2*(R2p2-5*ry**2)*(z+zi)+(8*c1*c2*R2p6*ry**2)/(R2p1+z+zi)**3-(4*c1*c2*R2p5*(R2p2-3*ry**2))/(R2p1+z+zi)**2-(4*c1*c2*R2p4*(R2p2-3*ry**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(3,1,3,2)=(3*rx*ry*(4*c1*c2*R2p4+c4*R2p2*(R2p2-5*z**2+5*zi**2)+10*zi*(-7*z*(z+zi)**2+R2p2*(2*z+zi))))/R2p9
    dudxdxi(1,2,3,2)=(rx*(6*(R2p4+35*ry**2*z**2-5*R2p2*(ry**2+z**2))*zi-30*(R2p2-7*ry**2)*z*zi**2+3*c4*R2p2*(R2p2-5*ry**2)*(z+zi)+(8*c1*c2*R2p6*ry**2)/(R2p1+z+zi)**3-(4*c1*c2*R2p5*(R2p2-3*ry**2))/(R2p1+z+zi)**2-(4*c1*c2*R2p4*(R2p2-3*ry**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(2,2,3,2)=-((ry*(4*c1*c2*R2p4*(6*R2p4+9*R2p3*(z+zi)-9*R2p1*ry**2*(z+zi)-3*ry**2*(z+zi)**2+R2p2*(-8*ry**2+3*(z+zi)**2))-3*(R2p1+z+zi)**3*(70*ry**2*z*zi*(z+zi)+R2p4*((-1+2*c4)*z+(5+2*c4)*zi)-5*R2p2*(c4*ry**2*(z+zi)+2*zi*(ry**2+3*z*(z+zi))))))/(R2p9*(R2p1+z+zi)**3))
    dudxdxi(3,2,3,2)=-(((4*c1*c2*R2p4*(R2p2-3*ry**2)+c4*R2p2*(R2p4+15*ry**2*(z-zi)*(z+zi)-3*R2p2*(ry**2+z**2-zi**2))+6*zi*(35*ry**2*z*(z+zi)**2+R2p4*(2*z+zi)-5*R2p2*(z*(z+zi)**2+ry**2*(2*z+zi)))))/R2p9)
    dudxdxi(1,3,3,2)=(-3*rx*ry*(4*c1*c2*R2p4-c4*R2p2*(R2p2-5*z**2+5*zi**2)+10*zi*(-7*z*(z+zi)**2+R2p2*(2*z+zi))))/R2p9
    dudxdxi(2,3,3,2)=((4*c1*c2*R2p4*(R2p2-3*ry**2)-c4*R2p2*(R2p4+15*ry**2*(z-zi)*(z+zi)-3*R2p2*(ry**2+z**2-zi**2))+6*zi*(35*ry**2*z*(z+zi)**2+R2p4*(2*z+zi)-5*R2p2*(z*(z+zi)**2+ry**2*(2*z+zi)))))/R2p9
    dudxdxi(3,3,3,2)=(3*ry*(-8*c1**2*R2p4*(z+zi)+c4*R2p2*(z+zi)*(3*R2p2-5*(z+zi)**2)-2*zi*(R2p4+35*z*(z+zi)**3-5*R2p2*(z+zi)*(4*z+zi))))/R2p9
    ! d u*_{lk} / d x_j / d x^i_3
    dudxdxi(1,1,1,3)=(rx*(4*c1*c2*R2p4*(6*R2p4+9*R2p3*(z+zi)-9*R2p1*rx**2*(z+zi)-3*rx**2*(z+zi)**2+R2p2*(-8*rx**2+3*(z+zi)**2))-3*(R2p1+z+zi)**3*(70*rx**2*z*zi*(z+zi)+R2p4*((5+2*c4)*z+(-1+2*c4)*zi)-5*R2p2*(6*z*zi*(z+zi)+rx**2*((2+c4)*z+c4*zi)))))/(R2p9*(R2p1+z+zi)**3)
    dudxdxi(2,1,1,3)=(ry*(-3*(2+c4)*R2p2*(R2p2-5*rx**2)*z-3*(c4*R2p2*(R2p2-5*rx**2)-10*(R2p2-7*rx**2)*z**2)*zi+30*(R2p2-7*rx**2)*z*zi**2-(8*c1*c2*R2p6*rx**2)/(R2p1+z+zi)**3+(4*c1*c2*R2p5*(R2p2-3*rx**2))/(R2p1+z+zi)**2+(4*c1*c2*R2p4*(R2p2-3*rx**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(3,1,1,3)=((4*c1*c2*R2p4*(R2p2-3*rx**2)-c4*R2p2*(R2p4+15*rx**2*(-z**2+zi**2)-3*R2p2*(rx**2-z**2+zi**2))+6*z*(35*rx**2*zi*(z+zi)**2+R2p4*(z+2*zi)-5*R2p2*(zi*(z+zi)**2+rx**2*(z+2*zi)))))/R2p9
    dudxdxi(1,2,1,3)=(ry*(-3*(2+c4)*R2p2*(R2p2-5*rx**2)*z-3*(c4*R2p2*(R2p2-5*rx**2)-10*(R2p2-7*rx**2)*z**2)*zi+30*(R2p2-7*rx**2)*z*zi**2-(8*c1*c2*R2p6*rx**2)/(R2p1+z+zi)**3+(4*c1*c2*R2p5*(R2p2-3*rx**2))/(R2p1+z+zi)**2+(4*c1*c2*R2p4*(R2p2-3*rx**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(2,2,1,3)=(rx*(4*c1*c2*R2p4*(2*R2p4+3*R2p3*(z+zi)-9*R2p1*ry**2*(z+zi)-3*ry**2*(z+zi)**2+R2p2*(-8*ry**2+(z+zi)**2))-3*(R2p1+z+zi)**3*(R2p4*(z-zi)+70*ry**2*z*zi*(z+zi)-5*R2p2*(2*z*zi*(z+zi)+ry**2*((2+c4)*z+c4*zi)))))/(R2p9*(R2p1+z+zi)**3)
    dudxdxi(3,2,1,3)=(-3*rx*ry*(4*c1*c2*R2p4-c4*R2p2*(R2p2+5*(z-zi)*(z+zi))+10*z*(-7*zi*(z+zi)**2+R2p2*(z+2*zi))))/R2p9
    dudxdxi(1,3,1,3)=-(((4*c1*c2*R2p4*(R2p2-3*rx**2)+c4*R2p2*(R2p4+15*rx**2*(-z**2+zi**2)-3*R2p2*(rx**2-z**2+zi**2))+6*z*(35*rx**2*zi*(z+zi)**2+R2p4*(z+2*zi)-5*R2p2*(zi*(z+zi)**2+rx**2*(z+2*zi)))))/R2p9)
    dudxdxi(2,3,1,3)=(3*rx*ry*(4*c1*c2*R2p4+c4*R2p2*(R2p2+5*(z-zi)*(z+zi))+10*z*(-7*zi*(z+zi)**2+R2p2*(z+2*zi))))/R2p9
    dudxdxi(3,3,1,3)=(3*rx*(70*z*zi*(z+zi)**3+R2p4*((2+8*c1**2-3*c4)*z+(8*c1**2-3*c4)*zi)+5*R2p2*(z+zi)*((-2+c4)*z**2+2*(-4+c4)*z*zi+c4*zi**2)))/R2p9
    dudxdxi(1,1,2,3)=(ry*(4*c1*c2*R2p4*(2*R2p4+3*R2p3*(z+zi)-9*R2p1*rx**2*(z+zi)-3*rx**2*(z+zi)**2+R2p2*(-8*rx**2+(z+zi)**2))-3*(R2p1+z+zi)**3*(R2p4*(z-zi)+70*rx**2*z*zi*(z+zi)-5*R2p2*(2*z*zi*(z+zi)+rx**2*((2+c4)*z+c4*zi)))))/(R2p9*(R2p1+z+zi)**3)
    dudxdxi(2,1,2,3)=(rx*(-3*(2+c4)*R2p2*(R2p2-5*ry**2)*z-3*(c4*R2p2*(R2p2-5*ry**2)-10*(R2p2-7*ry**2)*z**2)*zi+30*(R2p2-7*ry**2)*z*zi**2-(8*c1*c2*R2p6*ry**2)/(R2p1+z+zi)**3+(4*c1*c2*R2p5*(R2p2-3*ry**2))/(R2p1+z+zi)**2+(4*c1*c2*R2p4*(R2p2-3*ry**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(3,1,2,3)=(-3*rx*ry*(4*c1*c2*R2p4-c4*R2p2*(R2p2+5*(z-zi)*(z+zi))+10*z*(-7*zi*(z+zi)**2+R2p2*(z+2*zi))))/R2p9
    dudxdxi(1,2,2,3)=(rx*(-3*(2+c4)*R2p2*(R2p2-5*ry**2)*z-3*(c4*R2p2*(R2p2-5*ry**2)-10*(R2p2-7*ry**2)*z**2)*zi+30*(R2p2-7*ry**2)*z*zi**2-(8*c1*c2*R2p6*ry**2)/(R2p1+z+zi)**3+(4*c1*c2*R2p5*(R2p2-3*ry**2))/(R2p1+z+zi)**2+(4*c1*c2*R2p4*(R2p2-3*ry**2))/(R2p1+z+zi)))/R2p9
    dudxdxi(2,2,2,3)=(ry*(4*c1*c2*R2p4*(6*R2p4+9*R2p3*(z+zi)-9*R2p1*ry**2*(z+zi)-3*ry**2*(z+zi)**2+R2p2*(-8*ry**2+3*(z+zi)**2))-3*(R2p1+z+zi)**3*(70*ry**2*z*zi*(z+zi)+R2p4*((5+2*c4)*z+(-1+2*c4)*zi)-5*R2p2*(6*z*zi*(z+zi)+ry**2*((2+c4)*z+c4*zi)))))/(R2p9*(R2p1+z+zi)**3)
    dudxdxi(3,2,2,3)=((4*c1*c2*R2p4*(R2p2-3*ry**2)-c4*R2p2*(R2p4+15*ry**2*(-z**2+zi**2)-3*R2p2*(ry**2-z**2+zi**2))+6*z*(35*ry**2*zi*(z+zi)**2+R2p4*(z+2*zi)-5*R2p2*(zi*(z+zi)**2+ry**2*(z+2*zi)))))/R2p9
    dudxdxi(1,3,2,3)=(3*rx*ry*(4*c1*c2*R2p4+c4*R2p2*(R2p2+5*(z-zi)*(z+zi))+10*z*(-7*zi*(z+zi)**2+R2p2*(z+2*zi))))/R2p9
    dudxdxi(2,3,2,3)=-(((4*c1*c2*R2p4*(R2p2-3*ry**2)+c4*R2p2*(R2p4+15*ry**2*(-z**2+zi**2)-3*R2p2*(ry**2-z**2+zi**2))+6*z*(35*ry**2*zi*(z+zi)**2+R2p4*(z+2*zi)-5*R2p2*(zi*(z+zi)**2+ry**2*(z+2*zi)))))/R2p9)
    dudxdxi(3,3,2,3)=(3*ry*(70*z*zi*(z+zi)**3+R2p4*((2+8*c1**2-3*c4)*z+(8*c1**2-3*c4)*zi)+5*R2p2*(z+zi)*((-2+c4)*z**2+2*(-4+c4)*z*zi+c4*zi**2)))/R2p9
    dudxdxi(1,1,3,3)=(((1+4*c1*c2)*R2p6-210*rx**2*z*zi*(z+zi)**2-3*R2p4*((2+4*c1*c2+c4)*rx**2+z**2+4*z*zi+zi**2)+15*R2p2*(2*z*zi*(z+zi)**2+rx**2*((2+c4)*z**2+2*(3+c4)*z*zi+(2+c4)*zi**2))))/R2p9
    dudxdxi(2,1,3,3)=(-3*rx*ry*((2+4*c1*c2+c4)*R2p4+70*z*zi*(z+zi)**2-5*R2p2*((2+c4)*z**2+2*(3+c4)*z*zi+(2+c4)*zi**2)))/R2p9
    dudxdxi(3,1,3,3)=(3*rx*(70*z*zi*(z+zi)**3+R2p4*(-((-4+4*c1*c2+c4)*z)+(4-4*c1*c2+c4)*zi)+5*R2p2*(z+zi)*((-2+c4)*z**2-10*z*zi-(2+c4)*zi**2)))/R2p9
    dudxdxi(1,2,3,3)=(-3*rx*ry*((2+4*c1*c2+c4)*R2p4+70*z*zi*(z+zi)**2-5*R2p2*((2+c4)*z**2+2*(3+c4)*z*zi+(2+c4)*zi**2)))/R2p9
    dudxdxi(2,2,3,3)=(((1+4*c1*c2)*R2p6-210*ry**2*z*zi*(z+zi)**2-3*R2p4*((2+4*c1*c2+c4)*ry**2+z**2+4*z*zi+zi**2)+15*R2p2*(2*z*zi*(z+zi)**2+ry**2*((2+c4)*z**2+2*(3+c4)*z*zi+(2+c4)*zi**2))))/R2p9
    dudxdxi(3,2,3,3)=(3*ry*(70*z*zi*(z+zi)**3+R2p4*(-((-4+4*c1*c2+c4)*z)+(4-4*c1*c2+c4)*zi)+5*R2p2*(z+zi)*((-2+c4)*z**2-10*z*zi-(2+c4)*zi**2)))/R2p9
    dudxdxi(1,3,3,3)=(3*rx*(-70*z*zi*(z+zi)**3+R2p4*((-4+4*c1*c2-c4)*z+(-4+4*c1*c2+c4)*zi)+5*R2p2*(z+zi)*((2+c4)*z**2+10*z*zi-(-2+c4)*zi**2)))/R2p9
    dudxdxi(2,3,3,3)=(3*ry*(-70*z*zi*(z+zi)**3+R2p4*((-4+4*c1*c2-c4)*z+(-4+4*c1*c2+c4)*zi)+5*R2p2*(z+zi)*((2+c4)*z**2+10*z*zi-(-2+c4)*zi**2)))/R2p9
    dudxdxi(3,3,3,3)=-((((2+8*c1**2-3*c4)*R2p6-210*z*zi*(z+zi)**4-6*R2p4*((4+4*c1**2-3*c4)*z**2+(11+8*c1**2-6*c4)*z*zi+(4+4*c1**2-3*c4)*zi**2)-15*R2p2*(z+zi)**2*((-2+c4)*z**2+2*(-8+c4)*z*zi+(-2+c4)*zi**2)))/R2p9)
    ! Constants
    lambda=2.d0*mu*nu/(1.d0-2.d0*nu)
    Ku=1.d0/(16.d0*c_pi*mu*(1.d0-nu))
    dudxi=Ku*dudxi
    dudxdxi=Ku*dudxdxi
    ! d sigma*_{lkj} / d x^i_p
    do ip=1,3
      do ik=1,3
        do il=1,3
          do ij=1,3
            dsdxi(il,ik,ij,ip)=lambda*c_dkr(ik,ij)*(dudxdxi(il,1,1,ip)+dudxdxi(il,2,2,ip)+dudxdxi(il,3,3,ip))+mu*(dudxdxi(il,ik,ij,ip)+dudxdxi(il,ij,ik,ip))
          end do
        end do
      end do
    end do
    ! d t*_{lk} / d x^i_p
    do ip=1,3
      do ik=1,3
        do il=1,3
          dtdxi(il,ik,ip)=dot_product(dsdxi(il,ik,:,ip),n(:))
        end do
      end do
    end do
    ! d*_{lkp}
    do ip=1,3
      do ik=1,3
        do il=1,3
          dlkp(il,ik,ip)=lambda*c_dkr(il,ip)*(dudxi(1,ik,1)+dudxi(2,ik,2)+dudxi(3,ik,3))+mu*(dudxi(il,ik,ip)+dudxi(ip,ik,il))
        end do
      end do
    end do
    ! d*_{lk}
    do ik=1,3
      do il=1,3
        d(il,ik)=dot_product(dlkp(il,ik,:),n_i)
      end do
    end do
    ! s*_{lkp}
    do ip=1,3
      do ik=1,3
        do il=1,3
          slkp(il,ik,ip)=lambda*c_dkr(il,ip)*(dtdxi(1,ik,1)+dtdxi(2,ik,2)+dtdxi(3,ik,3))+mu*(dtdxi(il,ik,ip)+dtdxi(ip,ik,il))
        end do
      end do
    end do
    ! s*_{lk}
    do ik=1,3
      do il=1,3
        s(il,ik)=dot_product(slkp(il,ik,:),n_i)
      end do
    end do
  end subroutine

  subroutine fbem_bem_staela3d_hsc_hbie_ext_pre(ps,e,reverse,x_i,n_i,mu,nu,np,xp,m,l)
    implicit none
    ! I/O
    integer                :: ps                !! Selected precalculated dataset
    type(fbem_bem_element) :: e                 !! Element
    logical                :: reverse           !! Reverse normal vector
    real(kind=real64)      :: x_i(3)            !! Collocation point position vector
    real(kind=real64)      :: n_i(3)            !! Collocation point unit normal vector
    real(kind=real64)      :: mu                !! Shear modulus
    real(kind=real64)      :: nu                !! Poisson's ratio
    integer                :: np                !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64)      :: xp                !! Coordinate of the plane: x(np)
    real(kind=real64)      :: m(e%n_pnodes,3,3) !! m integration kernels matrix
    real(kind=real64)      :: l(e%n_snodes,3,3) !! l integration kernels matrix
    ! Local
    integer           :: il, ik               ! Counter for load / observation components
    integer           :: kip                  ! Counter variable for integration points loop
    real(kind=real64) :: x(3)                 ! Position vector at integration point
    real(kind=real64) :: n(3)                 ! Unit normal vector at integration point
    real(kind=real64) :: pphijw(e%n_pnodes)   ! phi^p * jacobian * weight at integration point
    real(kind=real64) :: sphijw(e%n_snodes)   ! phi^s * jacobian * weight at integration point
    real(kind=real64) :: fs_d(3,3), fs_s(3,3) ! Fundamental solution
    ! Initialize
    m=0.
    l=0.
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      x=e%ps_x(:,kip,ps)
      n=e%ps_n(:,kip,ps)
      pphijw=e%ps_pphijw(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      call fbem_bem_staela3d_hsc_hbie(mu,nu,np,xp,x_i,n_i,x,n,fs_d,fs_s)
      do il=1,3
        do ik=1,3
          m(:,il,ik)=m(:,il,ik)+fs_s(il,ik)*pphijw
          l(:,il,ik)=l(:,il,ik)+fs_d(il,ik)*sphijw
        end do
      end do
    end do
    ! Reverse if needed
    if (reverse) m=-m
  end subroutine fbem_bem_staela3d_hsc_hbie_ext_pre

  !! This subroutine calculates the kernels for SBIE exterior integration (near collocation points) using Telles transformation
  !! within a subdivision of the element, needing only needs basic data.
  subroutine fbem_bem_staela3d_hsc_hbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,mu,nu,np,xp,gln,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)       :: e                                !! Integration element
    logical                      :: reverse                          !! Reverse normal vector
    real(kind=real64)            :: xi_s(2,fbem_n_vertices(e%gtype)) !! Coordinates of the subdivision of the parent element
    real(kind=real64)            :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)            :: n_i(3)                           !! Collocation point unit normal
    real(kind=real64)            :: barxip(2)                        !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)            :: barr                             !! Telles jacobian at barxip
    real(kind=real64)            :: mu                               !! Shear modulus
    real(kind=real64)            :: nu                               !! Poisson's ratio
    integer                      :: np                               !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64)            :: xp                               !! Coordinate of the plane: x(np)
    integer                      :: gln                              !! 1D Gauss-Legendre number of integration points (<=32)
    real(kind=real64)            :: m(e%n_pnodes,3,3)                !! m kernel vector
    real(kind=real64)            :: l(e%n_snodes,3,3)                !! l kernel vector
    ! Local
    integer                      :: kphi                       ! Counter variable for shape functions loops
    integer                      :: k1                         ! Counter variable for reference coordinate xi_1
    integer                      :: k2                         ! Counter variable for reference coordinate xi_2
    real(kind=real64)            :: aux(10)                    ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)           ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi1(e%n_gnodes)      ! Geometrical shape functions derivatives values
    real(kind=real64)            :: dgphidxi2(e%n_gnodes)      ! Geometrical shape functions derivatives values
    real(kind=real64)            :: pphi(e%n_pnodes)           ! Functional shape functions values
    real(kind=real64)            :: sphi(e%n_snodes)           ! Functional shape functions values
    real(kind=real64)            :: gamma(2)                   ! Vector of gamma_1,gamma_2 coordinates
    real(kind=real64)            :: w(2)                       ! Weights of the integration rule
    real(kind=real64)            :: xip(2)                     ! Vector of xip_1,xip_2 coordinates
    real(kind=real64)            :: dxidxi1p(2)                ! xi derivatives with respect to xi1p
    real(kind=real64)            :: dxidxi2p(2)                ! xi derivatives with respect to xi2p
    real(kind=real64)            :: js                         ! Subdivision jacobian
    real(kind=real64)            :: xi(2)                      ! Vector of xi_1,xi_2 coordinates
    real(kind=real64)            :: xipp(2)                    ! Coordinate xipp used for quadrilateral-triangle transformation
    real(kind=real64)            :: barxipp(2)                 ! Coordinate xipp of collocation point
    real(kind=real64)            :: jqt                        ! Jacobian of the quadrilateral-triangle transformation
    type(fbem_telles_parameters) :: telles_parameters(2)       ! Telles parameters for each coordinate
    real(kind=real64)            :: jt(2)                      ! Telles jacobian for each coordinate: xi_1->gamma_1 and xi_2->gamma_2
    real(kind=real64)            :: x(3)                       ! Position vector at xi_1,xi_2
    real(kind=real64)            :: T1(3), T2(3)               ! Tangent vectors at xi_1,xi_2
    real(kind=real64)            :: N(3)                       ! Normal vector at xi_1,xi_2
    real(kind=real64)            :: jg                         ! Geometric jacobian
    real(kind=real64)            :: jw                         ! Jacobians * weights
    real(kind=real64)            :: pphijw(e%n_pnodes)         ! Functional shape functions * jw
    real(kind=real64)            :: sphijw(e%n_snodes)         ! Functional shape functions * jw
    integer                      :: il, ik                     ! Counter for load / observation components
    real(kind=real64)            :: fs_d(3,3), fs_s(3,3)       ! Fundamental solutions values
    ! Initialize
    m=0.
    l=0.
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
            ! Fundamental solution
            call fbem_bem_staela3d_hsc_hbie(mu,nu,np,xp,x_i,n_i,x,n,fs_d,fs_s)
            ! Add integration points
            do il=1,3
              do ik=1,3
                m(:,il,ik)=m(:,il,ik)+fs_s(il,ik)*pphijw
                l(:,il,ik)=l(:,il,ik)+fs_d(il,ik)*sphijw
              end do
            end do
          end do
        end do
      ! TRIANGULAR ELEMENTS
      case (3)
        ! Telles transformation is applied to Gauss-Legendre*Gauss-Legendre quadrature before the quadrilateral->triangle
        ! transformation. Because barxip for triangles are given in area triangle coordinates, for Telles transformation
        ! they must be transformed to quadrilateral coordinates. A special treatment is needed when barxi_2 is near 1, because
        ! transformation diverges.
        ! BARXIP->BARXIPP TRANSFORMATION
        if (barxip(2).gt.0.995d0) then
          barxipp(1)=0.5d0
          barxipp(2)=1.0d0
        else
          barxipp(1)=barxip(1)/(1.d0-barxip(2))
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
            ! Fundamental solution
            call fbem_bem_staela3d_hsc_hbie(mu,nu,np,xp,x_i,n_i,x,n,fs_d,fs_s)
            ! Add integration points
            do il=1,3
              do ik=1,3
                m(:,il,ik)=m(:,il,ik)+fs_s(il,ik)*pphijw
                l(:,il,ik)=l(:,il,ik)+fs_d(il,ik)*sphijw
              end do
            end do
          end do
        end do
        case default
          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'n_edges not valid')
    end select
    ! Reverse if needed
    if (reverse) m=-m
  end subroutine fbem_bem_staela3d_hsc_hbie_ext_st

  recursive subroutine fbem_bem_staela3d_hsc_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,mu,nu,np,xp,qsp,ks,ns,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                !! Element
    logical                  :: reverse                          !! Reverse orientation
    real(kind=real64)        :: xi_s(2,fbem_n_vertices(e%gtype)) !! Subdivision of the parent element
    real(kind=real64)        :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)        :: n_i(3)                           !! Collocation point unit normal
    real(kind=real64)        :: mu                               !! Shear modulus
    real(kind=real64)        :: nu                               !! Poisson's ratio
    integer                  :: np                               !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64)        :: xp                               !! Coordinate of the plane: x(np)
    type(fbem_qs_parameters) :: qsp                              !! Quasi-singular integration parameters
    integer                  :: ks                               !! Current level of subdivisions
    integer                  :: ns                               !! Maximum level of subdivision
    real(kind=real64)        :: m(e%n_pnodes,3,3)                !! m integration kernels matrix
    real(kind=real64)        :: l(e%n_snodes,3,3)                !! l integration kernels matrix
    ! Local
    real(kind=real64) :: xp_i(3)                              ! Image collocation point position vector
    integer           :: gln_near                             ! 1D Gauss-Legendre integ. points required to integrate only  the quasi-singular integrand
    integer           :: gln                                  ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
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
    real(kind=real64) :: m_tmp(e%n_pnodes,3,3)                ! m integration kernels matrix (temporary)
    real(kind=real64) :: l_tmp(e%n_snodes,3,3)                ! l integration kernels matrix (temporary)
    ! Collocation point image in global axes
    xp_i=x_i
    xp_i(abs(np))=2.*xp-xp_i(abs(np))
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
      call fbem_nearest_element_point_bem(3,e%gtype,e%x,e%cl,xp_i,barxi,rmin,d,method)
      barxip=barxi
    else
      call fbem_obtain_element_subdivision_coordinates(3,e%gtype,e%x,xi_s,x_s)
      cl=fbem_characteristic_length(3,e%gtype,x_s,1.d-12)
      call fbem_nearest_element_point_bem(3,e%gtype,x_s,cl,xp_i,barxip,rmin,d,method)
    end if
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,7,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_staela3d_hsc_hbie_ext_adp',ns,'maximum number of subdivisions reached')
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
          call fbem_bem_staela3d_hsc_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,mu,nu,np,xp,qsp,ks+1,ns,m,l)
          ! SUBTRI 2
          tmp_xi_s(:,1)=xi_s(:,2)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          call fbem_bem_staela3d_hsc_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,mu,nu,np,xp,qsp,ks+1,ns,m,l)
          ! SUBTRI 3
          tmp_xi_s(:,1)=xi_s(:,3)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          call fbem_bem_staela3d_hsc_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,mu,nu,np,xp,qsp,ks+1,ns,m,l)
          ! SUBTRI 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          call fbem_bem_staela3d_hsc_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,mu,nu,np,xp,qsp,ks+1,ns,m,l)
        ! QUADRILATERALS
        case (4)
          ! SUBQUAD 1
          tmp_xi_s(:,1)=xi_s(:,1)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,3)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          call fbem_bem_staela3d_hsc_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,mu,nu,np,xp,qsp,ks+1,ns,m,l)
          ! SUBQUAD 2
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=xi_s(:,2)
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,4)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          call fbem_bem_staela3d_hsc_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,mu,nu,np,xp,qsp,ks+1,ns,m,l)
          ! SUBQUAD 3
          tmp_xi_s(:,1)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=xi_s(:,3)
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          call fbem_bem_staela3d_hsc_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,mu,nu,np,xp,qsp,ks+1,ns,m,l)
          ! SUBQUAD 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          tmp_xi_s(:,2)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=xi_s(:,4)
          call fbem_bem_staela3d_hsc_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,mu,nu,np,xp,qsp,ks+1,ns,m,l)
      end select
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_staela3d_hsc_hbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,mu,nu,np,xp,gln,m_tmp,l_tmp)
      m=m+m_tmp
      l=l+l_tmp
    end if
  end subroutine fbem_bem_staela3d_hsc_hbie_ext_adp

  !! Efficient automatic integration
  subroutine fbem_bem_staela3d_hsc_hbie_auto(e,reverse,x_i,n_i,mu,nu,np,xp,qsp,ns,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                !! Integration element
    logical                  :: reverse                          !! Reverse orientation
    real(kind=real64)        :: x_i(3)                           !! Collocation point
    real(kind=real64)        :: n_i(3)                           !! Collocation point unit normal
    real(kind=real64)        :: mu                               !! Shear modulus
    real(kind=real64)        :: nu                               !! Poisson's ratio
    integer                  :: np                               !! Unit normal of the plane: -1,1,-2,2,-3,3
    real(kind=real64)        :: xp                               !! Coordinate of the plane: x(np)
    type(fbem_qs_parameters) :: qsp                              !! Quasi-singular integration parameters
    integer                  :: ns                               !! Maximum level of subdivisions
    real(kind=real64)        :: m(e%n_pnodes,3,3)                !! m integration kernel
    real(kind=real64)        :: l(e%n_snodes,3,3)                !! l integration kernel
    ! Local
    real(kind=real64) :: xp_i(3)                          ! Collocation point image in global axes
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
    ! Collocation point image in global axes
    xp_i=x_i
    xp_i(abs(np))=2.*xp-xp_i(abs(np))
    ! Determine if interior or exterior integration
    !   - Interior integration (delta=1) requires: xi_i
    !   - Exterior integration (delta=0) requires: xp_i, barxi, rmin and d
    ! Use the element ball
    r=e%bball_centre-xp_i
    rmin=sqrt(dot_product(r,r))-e%bball_radius
    if (rmin.gt.(4.d0*e%bball_radius)) then
      delta=0
      barxi=0.d0
      d=rmin/e%cl
    else
      ! Use an adaptative algorithm that combines sampling and minimization algorithms
      call fbem_nearest_element_point_bem(3,e%gtype,e%x,e%cl,xp_i,barxi,rmin,d,method)
      if (d.le.1.d-12) then
        delta=1
      else
        delta=0
      end if
    end if
    ! Integrate
    select case (delta)
      case (1)
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'HBIE of complementary Mindlin cannot be integrated for interior collocation points yet')
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
          call fbem_bem_staela3d_hsc_hbie_ext_pre(ps,e,reverse,x_i,n_i,mu,nu,np,xp,m,l)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_staela3d_hsc_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,mu,nu,np,xp,qsp,1,ns,m,l)
        end if
    end select
  end subroutine fbem_bem_staela3d_hsc_hbie_auto

end module fbem_bem_staela3d
