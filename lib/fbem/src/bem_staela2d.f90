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
!! <b> This module implements the calculation of BEM integration kernels matrixs for 2D elastostatic problems.</b>
module fbem_bem_staela2d

  ! Fortran 2003 standard
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_quad_rules
  use fbem_telles_transformation
  use fbem_polar_transformation
  use fbem_geometry
  use fbem_quasisingular_integration
  use fbem_bem_general

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all are private
  private

  ! ================================================================================================================================
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE) - FULL-SPACE FUNDAMENTAL SOLUTION
  ! Free-term c
  public :: fbem_bem_staela2d_sbie_freeterm
  ! Fundamental solution
  public :: fbem_bem_staela2d_sbie_u
  public :: fbem_bem_staela2d_sbie_t
  ! BOUNDARY ELEMENTS
  ! Automatic integration
  public :: fbem_bem_staela2d_sbie_auto
  ! Exterior integration
  public :: fbem_bem_staela2d_sbie_ext_pre
  public :: fbem_bem_staela2d_sbie_ext_st
  public :: fbem_bem_staela2d_sbie_ext_adp
  ! Interior integration
  public :: fbem_bem_staela2d_sbie_int
  ! BODY LOAD ELEMENTS
  ! Automatic integration
  public :: fbem_bem_staela2d_sbie_bl_auto
  ! Exterior integration
  public :: fbem_bem_staela2d_sbie_bl_ext_pre
  public :: fbem_bem_staela2d_sbie_bl_ext_st
  public :: fbem_bem_staela2d_sbie_bl_ext_adp
  ! Interior integration
  public :: fbem_bem_staela2d_sbie_bl_int
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE) - HALF-SPACE FUNDAMENTAL SOLUTION (COMPLEMENTARY PART)
  ! Fundamental solution
  public :: fbem_bem_staela2d_hfc_sbie
  ! Automatic integration
  public :: fbem_bem_staela2d_hfc_sbie_auto
  ! Exterior integration
  public :: fbem_bem_staela2d_hfc_sbie_ext_pre
  public :: fbem_bem_staela2d_hfc_sbie_ext_st
  public :: fbem_bem_staela2d_hfc_sbie_ext_adp
  ! Interior integration
  !public :: fbem_bem_staela2d_hfc_sbie_int
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (HBIE) - FULL-SPACE FUNDAMENTAL SOLUTION
  ! Fundamental solution
  public :: fbem_bem_staela2d_hbie_d
  public :: fbem_bem_staela2d_hbie_s
  ! BOUNDARY ELEMENTS
  ! Automatic integration
  public :: fbem_bem_staela2d_hbie_auto
  ! Exterior integration
  public :: fbem_bem_staela2d_hbie_ext_pre
  public :: fbem_bem_staela2d_hbie_ext_st
  public :: fbem_bem_staela2d_hbie_ext_adp
  ! Interior integration
  public :: fbem_bem_staela2d_hbie_int
  ! BODY LOAD ELEMENTS
  ! Automatic integration
  public :: fbem_bem_staela2d_hbie_bl_auto
  ! Exterior integration
  public :: fbem_bem_staela2d_hbie_bl_ext_pre
  public :: fbem_bem_staela2d_hbie_bl_ext_st
  public :: fbem_bem_staela2d_hbie_bl_ext_adp
  ! Interior integration
  public :: fbem_bem_staela2d_hbie_bl_int
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! VARIATION SINGULAR BOUNDARY INTEGRAL EQUATION (VSBIE)
  ! Free-term b
  public :: fbem_bem_staela2d_vsbie_freeterm
  ! Automatic integration
  public :: fbem_bem_staela2d_vsbie_auto
  ! Exterior integration
  public :: fbem_bem_staela2d_vsbie_ext_pre
  public :: fbem_bem_staela2d_vsbie_ext_st
  public :: fbem_bem_staela2d_vsbie_ext_adp
  ! Interior integration
  public :: fbem_bem_staela2d_vsbie_int
  ! ================================================================================================================================

contains

  ! ================================================================================================================================
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE) - FULL-SPACE FUNDAMENTAL SOLUTION

  !! Subroutine that calculates free-term of 2D potential problem
  subroutine fbem_bem_staela2d_sbie_freeterm(n_elements,n,tc,tol,nu,c)
    implicit none
    ! I/O
    real(kind=real64) :: c(2,2)        !! Free-term
    integer           :: n_elements    !! Number of elements (1 or 2)
    real(kind=real64) :: n(2,2)        !! Normals of elements at collocation point
    real(kind=real64) :: tc(2,2)       !! Tangents of elements at collocation point towards inside the element
    real(kind=real64) :: tol           !! Geometric tolerance [1.0e-12,1.0e-3] (default 1.0e-6)
    real(kind=real64) :: nu            !! Poisson's ratio
    ! Local
    real(kind=real64) :: local_n(2,2)  ! Local copy of n
    real(kind=real64) :: local_tc(2,2) ! Local copy of tc
    real(kind=real64) :: local_tol     ! Local copy of geometric tolerance
    real(kind=real64) :: nsum(2)       ! Unit sum of normals
    real(kind=real64) :: vectornorm    ! Norm of a vector
    real(kind=real64) :: theta(2)      ! Angles of unit tangents
    real(kind=real64) :: theta_ext     ! Angle theta exterior
    real(kind=real64) :: cte, ctesin   ! Auxiliary variable
    integer           :: i             ! Counter
    !
    ! Check n_elements
    !
    if ((n_elements.lt.1).or.(n_elements.gt.2)) then
      call fbem_error_message(output_unit,0,__FILE__,__LINE__,&
                                'n_elements must be 1 or 2 for free-term calculation.')
    end if
    !
    ! Check geometric tolerance
    !
    if ((tol.lt.1.0d-12).or.(tol.gt.1.0d-3)) then
      local_tol=1.0d-6
    else
      local_tol=tol
    end if
    !
    ! Check that n and tc are orthogonal (using scalar product)
    !
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
        c(1,1)=0.5d0
        c(1,2)=0.0d0
        c(2,1)=0.0d0
        c(2,2)=0.5d0
      ! For 2 elements, the freeterm is calculated as a vertex
      case (2)
        !
        ! Know who is the element 1 and 2.
        !
        ! Add up normals and normalize it
        nsum(1)=n(1,1)+n(1,2)
        nsum(2)=n(2,1)+n(2,2)
        vectornorm=dsqrt(nsum(1)**2+nsum(2)**2)
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
        theta(1)=datan2(local_tc(2,1),local_tc(1,1))
        theta(2)=datan2(local_tc(2,2),local_tc(1,2))
        ! Convert to (0,2*pi)
        if (theta(1).lt.0.0d0) theta(1)=theta(1)+c_2pi
        if (theta(2).lt.0.0d0) theta(2)=theta(2)+c_2pi
        ! Use theta2>theta1
        if (theta(1).gt.theta(2)) theta(2)=theta(2)+c_2pi
        ! Calculations
        theta_ext=theta(2)-theta(1)
        cte=1.0d0/(8.0d0*c_pi*(1.0d0-nu))
        theta(1)=2.0d0*theta(1)
        theta(2)=2.0d0*theta(2)
        ctesin=dsin(theta(2))-dsin(theta(1))
        ! Assign
        c(1,1)=1.0d0-theta_ext*c_1_2pi-cte*ctesin
        c(1,2)=cte*(dcos(theta(2))-dcos(theta(1)))
        c(2,1)=c(1,2)
        c(2,2)=1.0d0-theta_ext*c_1_2pi+cte*ctesin
    end select
  end subroutine fbem_bem_staela2d_sbie_freeterm

  !! Fundamental solution u*
  subroutine fbem_bem_staela2d_sbie_u(x,x_i,mu,nu,uo)
    implicit none
    ! I/O
    real(kind=real64) :: x(2)    !! Observation point
    real(kind=real64) :: x_i(2)  !! Collocation point
    real(kind=real64) :: mu      !! Shear modulus
    real(kind=real64) :: nu      !! Poisson's ratio
    real(kind=real64) :: uo(2,2) !! u*_{lk}
    ! Local
    integer           :: l, k
    real(kind=real64) :: cteu1, cteu2
    real(kind=real64) :: rv(2), r, d1r, logd1r, drdx(2)
    cteu1=1.0d0/(8.0d0*c_pi*mu*(1.0d0-nu))
    cteu2=3.0d0-4.0d0*nu
    rv=x-x_i
    r=dot_product(rv,rv)
    r=sqrt(r)
    d1r=1.d0/r
    logd1r=dlog(d1r)
    drdx=rv*d1r
    do l=1,2
      do k=1,2
        uo(l,k)=cteu2*logd1r*c_dkr(l,k)+drdx(l)*drdx(k)
      end do
    end do
    uo=cteu1*uo
  end subroutine fbem_bem_staela2d_sbie_u

  !! Fundamental solution t*
  subroutine fbem_bem_staela2d_sbie_t(x,n,x_i,nu,to)
    implicit none
    ! I/O
    real(kind=real64) :: x(2)    !! Observation point
    real(kind=real64) :: n(2)    !! Observation point normal
    real(kind=real64) :: x_i(2)  !! Collocation point
    real(kind=real64) :: nu      !! Poisson's ratio
    real(kind=real64) :: to(2,2) !! t*_{lk}
    ! Local
    integer           :: l, k
    real(kind=real64) :: ctep1, ctep2
    real(kind=real64) :: rv(2), r, d1r, logd1r, drdx(2), drdn
    ctep1=-1.0d0/(4.0d0*c_pi*(1.0d0-nu))
    ctep2=1.0d0-2.0d0*nu
    rv=x-x_i
    r=dot_product(rv,rv)
    r=sqrt(r)
    d1r=1.d0/r
    logd1r=dlog(d1r)
    drdx=rv*d1r
    drdn=dot_product(drdx,n)
    do l=1,2
      do k=1,2
        to(l,k)=d1r*(drdn*(ctep2*c_dkr(l,k)+2.d0*drdx(l)*drdx(k))+ctep2*(n(l)*drdx(k)-n(k)*drdx(l)))
      end do
    end do
    to=ctep1*to
  end subroutine fbem_bem_staela2d_sbie_t

  ! ====================
  ! BE BOUNDARY ELEMENTS
  ! ====================

  !! Efficient automatic integration of boundary elements
  subroutine fbem_bem_staela2d_sbie_auto(e,reverse,x_i,mu,nu,qsp,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                 !! Integration element
    logical                  :: reverse           !! Reverse orientation
    real(kind=real64)        :: x_i(2)            !! Collocation point
    real(kind=real64)        :: mu                !! Shear modulus
    real(kind=real64)        :: nu                !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp               !! Quasi-singular integration parameters
    integer                  :: ns                !! Maximum level of subdivisions
    real(kind=real64)        :: h(e%n_pnodes,2,2) !! h integration kernel
    real(kind=real64)        :: g(e%n_snodes,2,2) !! g integration kernel
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
        call fbem_bem_staela2d_sbie_int(e,reverse,barxi,mu,nu,h,g)
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
          call fbem_bem_staela2d_sbie_ext_pre(ps,e,reverse,x_i,mu,nu,h,g)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_staela2d_sbie_ext_adp(e,reverse,xi_s,x_i,mu,nu,qsp,1,ns,h,g)
        end if
    end select
  end subroutine fbem_bem_staela2d_sbie_auto

!  !! This subroutine calculates the kernels for SBIE exterior integration (far collocation points), needing precalculated data at
!  !! integration points.
!  !!
!  !! It is useful because for a given element, the vast majority of kernels are associated with far collocation points, which
!  !! need the same number of Gauss points. To generate integrand components, the subroutine <tt>fbem_gauss_points_values2d</tt>
!  !! can be called (from geometry module).
!  subroutine fbem_bem_staela2d_sbie_ext_pre(rule,max_rule,ngp,max_ngp,type_f1,type_f2,phi_f1,phi_f2,x,n,jw,reverse,x_i,mu,nu,h,g)
!    implicit none
!    ! I/O
!    integer           :: rule                                           !! Selected rule
!    integer           :: max_rule                                       !! Total number of rules in the rules set-up
!    integer           :: ngp                                            !! Number of Gaussian point of the selected rule
!    integer           :: max_ngp                                        !! Maximum number of Gaussian points in the rules set-up
!    integer           :: type_f1                                        !! Functional interpolation (primary variables)
!    integer           :: type_f2                                        !! Geometrical interpolation (secondary variables)
!    real(kind=real64) :: phi_f1(fbem_n_nodes(type_f1),max_ngp,max_rule) !! Functional shape functions at integration points
!    real(kind=real64) :: phi_f2(fbem_n_nodes(type_f2),max_ngp,max_rule) !! Functional shape functions at integration points
!    real(kind=real64) :: x(2,max_ngp,max_rule)                          !! Position vectors at integration points
!    real(kind=real64) :: n(2,max_ngp,max_rule)                          !! Unit normal vector at integration points
!    real(kind=real64) :: jw(max_ngp,max_rule)                           !! Geometric jacobian multiplied by weight at integration points
!    logical           :: reverse                                        !! Reverse normal vector
!    real(kind=real64) :: x_i(2)                                         !! Collocation point position vector
!    real(kind=real64) :: mu                                             !! Shear modulus
!    real(kind=real64) :: nu                                             !! Poisson's ratio
!    real(kind=real64) :: h(fbem_n_nodes(type_f1),2,2)                   !! h integration kernels matrix
!    real(kind=real64) :: g(fbem_n_nodes(type_f2),2,2)                   !! g integration kernels matrix
!    ! Local
!    integer           :: l, k                                        ! Counter variables for load direction and observation direction
!    integer           :: kip                                         ! Counter variable for integration points loop
!    real(kind=real64) :: rv(2)                                       ! Distance vector between collocation point and integration point (x-x_i)
!    real(kind=real64) :: r, d1r, logd1r                              ! Distance vector module, inverse and log(1/r)
!    real(kind=real64) :: drdx(2)                                     ! Distance vector derivatives with respect to x_k
!    real(kind=real64) :: drdn                                        ! Partial derivative of r respect to unit normal
!    real(kind=real64) :: phif1jw(fbem_n_nodes(type_f1))              ! Auxiliary variables for integrand evaluation
!    real(kind=real64) :: phif2jw(fbem_n_nodes(type_f2))              ! Auxiliary variables for integrand evaluation
!    real(kind=real64) :: cteu1, cteu2, ctep1, ctep2                  ! Auxiliary constants
!    real(kind=real64) :: fs_u, fs_p                                  ! Fundamental solutions values
!    !
!    ! Initialize
!    !
!    ! Kernel matrices
!    h=0.0d0
!    g=0.0d0
!    ! Initialize auxiliary constants for fundamental solutions calculation
!    cteu1=1.0d0/(8.0d0*c_pi*mu*(1.0d0-nu))
!    cteu2=3.0d0-4.0d0*nu
!    ctep1=-1.0d0/(4.0d0*c_pi*(1.0d0-nu))
!    ctep2=1.0d0-2.0d0*nu
!    ! Numerical integration
!    ! Loop through integration points
!    do kip=1,ngp
!      ! Distance vector
!      rv=x(:,kip,rule)-x_i
!      ! Distance vector norm and its inverse
!      r=dot_product(rv,rv)
!      r=sqrt(r)
!      d1r=1.0d0/r
!      logd1r=dlog(d1r)
!      ! r_{,k}
!      drdx=rv*d1r
!      ! dr/dn
!      drdn=dot_product(drdx,n(:,kip,rule))
!      ! Shape functions * jw
!      phif1jw=phi_f1(:,kip,rule)*jw(kip,rule)
!      phif2jw=phi_f2(:,kip,rule)*jw(kip,rule)
!      ! Loop through load direction and observation direction
!      do l=1,2
!        do k=1,2
!          ! Fundamental solutions values without cteu1 and ctep1, respectively
!          fs_u=cteu2*logd1r*c_dkr(l,k)+drdx(l)*drdx(k)
!          fs_p=d1r*(drdn*(ctep2*c_dkr(l,k)+2.d0*drdx(l)*drdx(k))+ctep2*(n(l,kip,rule)*drdx(k)-n(k,kip,rule)*drdx(l)))
!          ! Add to kernels
!          h(:,l,k)=h(:,l,k)+fs_p*phif1jw
!          g(:,l,k)=g(:,l,k)+fs_u*phif2jw
!        end do
!      end do
!    end do ! Loop through integrations points
!    ! Multiply h by ctep1 and g by cteu1
!    h=ctep1*h
!    g=cteu1*g
!    ! If the normal has to be reversed, then h=-h
!    if (reverse) h=-h
!  end subroutine fbem_bem_staela2d_sbie_ext_pre

  !! This subroutine calculates the kernels for SBIE exterior integration using precalculated data at integration points.
  subroutine fbem_bem_staela2d_sbie_ext_pre(ps,e,reverse,x_i,mu,nu,h,g)
    implicit none
    ! I/O
    integer                :: ps                !! Selected precalculated dataset
    type(fbem_bem_element) :: e                 !! Element
    logical                :: reverse           !! Reverse normal vector
    real(kind=real64)      :: x_i(2)            !! Collocation point position vector
    real(kind=real64)      :: mu                !! Shear modulus
    real(kind=real64)      :: nu                !! Poisson's ratio
    real(kind=real64)      :: h(e%n_pnodes,2,2) !! h integration kernels matrix
    real(kind=real64)      :: g(e%n_snodes,2,2) !! g integration kernels matrix
    ! Local
    integer           :: l, k                       ! Counter variables for load direction and observation direction
    integer           :: kip                        ! Counter variable for integration points loop
    real(kind=real64) :: x(2)                       ! Position vector at integration point
    real(kind=real64) :: n(2)                       ! Unit normal vector at integration point
    real(kind=real64) :: pphijw(e%n_pnodes)         ! phi^p * jacobian * weight at integration point
    real(kind=real64) :: sphijw(e%n_snodes)         ! phi^s * jacobian * weight at integration point
    real(kind=real64) :: rv(2)                      ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64) :: r, d1r, logd1r             ! Distance vector module, inverse and log(1/r)
    real(kind=real64) :: drdx(2)                    ! Distance vector derivatives with respect to x_k
    real(kind=real64) :: drdn                       ! Partial derivative of r respect to unit normal
    real(kind=real64) :: cteu1, cteu2, ctet1, ctet2 ! Auxiliary constants
    real(kind=real64) :: fs_u, fs_t                 ! Fundamental solutions values
    ! Initialize kernels
    h=0.d0
    g=0.d0
    ! Initialize auxiliary constants for fundamental solutions calculation
    cteu1=1.0d0/(8.0d0*c_pi*mu*(1.0d0-nu))
    cteu2=3.0d0-4.0d0*nu
    ctet1=-1.0d0/(4.0d0*c_pi*(1.0d0-nu))
    ctet2=1.0d0-2.0d0*nu
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      x=e%ps_x(:,kip,ps)
      n=e%ps_n(:,kip,ps)
      pphijw=e%ps_pphijw(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      d1r=1.d0/r
      logd1r=log(d1r)
      drdx=rv*d1r
      drdn=dot_product(drdx,n)
      do l=1,2
        do k=1,2
          fs_u=cteu2*logd1r*c_dkr(l,k)+drdx(l)*drdx(k)
          fs_t=d1r*(drdn*(ctet2*c_dkr(l,k)+2.d0*drdx(l)*drdx(k))+ctet2*(n(l)*drdx(k)-n(k)*drdx(l)))
          h(:,l,k)=h(:,l,k)+fs_t*pphijw
          g(:,l,k)=g(:,l,k)+fs_u*sphijw
        end do
      end do
    end do
    ! Multiply by constants
    h=ctet1*h
    g=cteu1*g
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_staela2d_sbie_ext_pre

!  !! This subroutine calculates the kernels for SBIE exterior integration (near collocation points) using Telles transformation
!  !! within a subdivision of the element, needing only needs basic data.
!  subroutine fbem_bem_staela2d_sbie_ext_st(ngp,xi1,xi2,barxi,barr,type_g,type_f1,type_f2,delta_f,x_nodes,reverse,x_i,mu,nu,h,g)
!    implicit none
!    ! I/O
!    integer                      :: ngp                             !! Number of Gauss points (<=32)
!    real(kind=real64)            :: barxi                           !! Nearest element coordinate with respect to collocation point
!    real(kind=real64)            :: barr                            !! Telles jacobian at barxi (xi1<=barxi<=xi2)
!    real(kind=real64)            :: xi1                             !! Minimum xi coordinate of the subdivision (xi1)
!    real(kind=real64)            :: xi2                             !! Maximum xi coordinate of the subdivision (xi2)
!    integer                      :: type_g                          !! Geometrial interpolation
!    integer                      :: type_f1                         !! Functional interpolation (primary variables)
!    integer                      :: type_f2                         !! Functional interpolation (secondary variables)
!    real(kind=real64)            :: delta_f                         !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
!    real(kind=real64)            :: x_nodes(2,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
!    logical                      :: reverse                         !! Reverse normal vector
!    real(kind=real64)            :: x_i(2)                          !! Collocation point position vector
!    real(kind=real64)            :: mu                              !! Shear modulus
!    real(kind=real64)            :: nu                              !! Poisson's ratio
!    real(kind=real64)            :: h(fbem_n_nodes(type_f1),2,2)    !! h kernel vector
!    real(kind=real64)            :: g(fbem_n_nodes(type_f2),2,2)    !! g kernel vector
!    ! Local
!    integer                      :: l, k                            ! Counter variables for load direction and observation direction
!    integer                      :: kphi                            ! Counter variable for shape functions loops
!    integer                      :: nnodes_g                        ! Number of nodes of the element
!    integer                      :: kip                             ! Counter variable of integration points
!    real(kind=real64)            :: gamma                           ! Coordinate gamma (Telles transformation space [0,1])
!    real(kind=real64)            :: w                               ! Weights of an integration point
!    type(fbem_telles_parameters) :: telles_parameters               ! Telles parameters
!    real(kind=real64)            :: jt                              ! Telles jacobian: xip->gamma
!    real(kind=real64)            :: xip                             ! Coordinate xip (subdivision space [0,1])
!    real(kind=real64)            :: barxip                          ! Barxi in subdivision space
!    real(kind=real64)            :: js                              ! Subdivision jacobian: xi->xip
!    real(kind=real64)            :: xi                              ! Coordinate xi [xi1,xi2]
!    real(kind=real64)            :: aux(10)                         ! Auxiliary variable needed for shape_functions module resources
!    real(kind=real64)            :: phi_g(fbem_n_nodes(type_g))     ! Geometrical shape functions values
!    real(kind=real64)            :: dphidxi_g(fbem_n_nodes(type_g)) ! Geometrical shape functions first derivatives values
!    real(kind=real64)            :: phi_f1(fbem_n_nodes(type_f1))   ! Functional shape functions values
!    real(kind=real64)            :: phi_f2(fbem_n_nodes(type_f2))   ! Functional shape functions values
!    real(kind=real64)            :: x(2)                            ! Position vector at xi
!    real(kind=real64)            :: T(2)                            ! Tangent vector at xi
!    real(kind=real64)            :: N(2)                            ! Normal vector at xi
!    real(kind=real64)            :: rv(2)                           ! Distance vector between collocation point and integration point (x-x_i)
!    real(kind=real64)            :: r, d1r, logd1r                  ! Distance vector module, its inverse and log(1/r)
!    real(kind=real64)            :: drdx(2)                         ! Distance vector derivatives with respect to x_k
!    real(kind=real64)            :: jg                              ! Geometric jacobian
!    real(kind=real64)            :: drdn                            ! Partial derivative of r respect to unit normal
!    real(kind=real64)            :: jw                              ! Jacobians * weight
!    real(kind=real64)            :: phif1jw(fbem_n_nodes(type_f1))  ! shape functions * jw
!    real(kind=real64)            :: phif2jw(fbem_n_nodes(type_f2))  ! shape functions * jw
!    real(kind=real64)            :: cteu1, cteu2, ctep1, ctep2      ! Auxiliary constants
!    real(kind=real64)            :: fs_u, fs_p                      ! Fundamental solutions values
!    !
!    ! Initialization
!    !
!    ! Number of nodes of the geometrical interpolation element
!    nnodes_g=fbem_n_nodes(type_g)
!    ! Initialize kernel matrices
!    h=0.d0
!    g=0.d0
!    ! Initialize auxiliary constants for fundamental solutions calculation
!    cteu1=1.0d0/(8.0d0*c_pi*mu*(1.0d0-nu))
!    cteu2=3.0d0-4.0d0*nu
!    ctep1=-1.0d0/(4.0d0*c_pi*(1.0d0-nu))
!    ctep2=1.0d0-2.0d0*nu
!    ! Subdivision jacobian (is constant)
!    js=xi2-xi1
!    ! Convert barxi (which is in xi space [xi1,xi2]) to the subdivision space xip [0,1]
!    barxip=(barxi-xi1)/js
!    ! Calculate Telles parameters
!    telles_parameters=fbem_telles01_calculate_parameters(barxip,barr)
!    ! Loop through integrations points
!    do kip=1,gl01_n(ngp)
!      ! GAMMA->XIP->XI COORDINATE TRANSFORMATION
!      ! Gamma coordinate and weight
!      gamma=gl01_xi(kip,ngp)
!      w=gl01_w(kip,ngp)
!      ! xip coordinate, weight and jacobian from Telles transformation
!      call fbem_telles_xi_and_jacobian(telles_parameters,gamma,xip,jt)
!      ! xi coordinate
!      xi=js*xip+xi1
!      ! COMPONENTS OF THE FUNDAMENTAL SOLUTIONS
!      ! Geometrical shape functions and first derivatives at xi
!#     define etype type_g
!#     define delta 0.0d0
!#     define phi phi_g
!#     define dphidxi dphidxi_g
!#     include <phi_and_dphidxi_1d.rc>
!#     undef etype
!#     undef delta
!#     undef phi
!#     undef dphidxi
!      ! Components calculation of x and T at xi
!      x=0.d0
!      T=0.d0
!      do kphi=1,nnodes_g
!        x=x+phi_g(kphi)*x_nodes(:,kphi)
!        T=T+dphidxi_g(kphi)*x_nodes(:,kphi)
!      end do
!      ! Normal vector
!      N(1)=T(2)
!      N(2)=-T(1)
!      ! Geometric jacobian
!      jg=dot_product(T,T)
!      jg=sqrt(jg)
!      ! Unit normal vector
!      n=N/jg
!      ! Distance vector
!      rv=x-x_i
!      ! Distance vector norm and its inverse
!      r=dot_product(rv,rv)
!      r=sqrt(r)
!      d1r=1.d0/r
!      logd1r=dlog(d1r)
!      ! r_{,k}
!      drdx=rv*d1r
!      ! dr/dn
!      drdn=dot_product(drdx,n)
!      ! Jacobian * weight
!      jw=jg*js*jt*w
!      ! FUNCTIONAL SHAPE FUNCTIONS
!      ! Functional shape functions (primary variables) at xi
!#     define etype type_f1
!#     define delta delta_f
!#     define phi phi_f1
!#     include <phi_1d.rc>
!#     undef etype
!#     undef delta
!#     undef phi
!      ! Functional shape functions (secondary variables) at xi
!#     define etype type_f2
!#     define delta delta_f
!#     define phi phi_f2
!#     include <phi_1d.rc>
!#     undef etype
!#     undef delta
!#     undef phi
!      ! Functional shape functions * jacobians* weights
!      phif1jw=phi_f1*jw
!      phif2jw=phi_f2*jw
!      ! Loop through load and observation directions
!      do l=1,2
!        do k=1,2
!          ! Fundamental solutions values without cteu1 and ctep1, respectively
!          fs_u=cteu2*logd1r*c_dkr(l,k)+drdx(l)*drdx(k)
!          fs_p=d1r*(drdn*(ctep2*c_dkr(l,k)+2.d0*drdx(l)*drdx(k))+ctep2*(n(l)*drdx(k)-n(k)*drdx(l)))
!          ! Add to kernels
!          h(:,l,k)=h(:,l,k)+fs_p*phif1jw
!          g(:,l,k)=g(:,l,k)+fs_u*phif2jw
!        end do
!      end do
!    end do ! Loop through integrations points
!    ! Multiply h by ctep1 and g by cteu1
!    h=ctep1*h
!    g=cteu1*g
!    ! If the normal has to be reversed, then h=-h
!    if (reverse) h=-h
!  end subroutine fbem_bem_staela2d_sbie_ext_st

  !! This subroutine calculates the kernels for SBIE exterior integration (near collocation points) using Telles transformation
  !! within a subdivision of the element, needing only needs basic data.
  subroutine fbem_bem_staela2d_sbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,mu,nu,gln,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)       :: e                 !! Integration element
    logical                      :: reverse           !! Reverse normal vector
    real(kind=real64)            :: xi_s(1,2)         !! Coordinates of the subdivision of the parent element (xi space [-1,1])
    real(kind=real64)            :: x_i(2)            !! Collocation point position vector
    real(kind=real64)            :: barxip(1)         !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)            :: barr              !! Telles jacobian at barxip
    real(kind=real64)            :: mu                !! Shear modulus
    real(kind=real64)            :: nu                !! Poisson's ratio
    integer                      :: gln               !! 1D Gauss-Legendre number of integration points (<=32)
    real(kind=real64)            :: h(e%n_pnodes,2,2) !! h kernel vector
    real(kind=real64)            :: g(e%n_snodes,2,2) !! g kernel vector
    ! Local
    integer                      :: l, k                       ! Counter variables for load direction and observation direction
    integer                      :: kphi                       ! Counter variable for shape functions loops
    integer                      :: kip                        ! Counter variable of integration points
    real(kind=real64)            :: gamma                      ! Coordinate gamma (Telles transformation space [-1,1])
    real(kind=real64)            :: w                          ! Weights of an integration point
    type(fbem_telles_parameters) :: telles_parameters          ! Telles parameters
    real(kind=real64)            :: jt                         ! Telles jacobian: xip->gamma
    real(kind=real64)            :: xip                        ! Coordinate xip (subdivision space [-1,1])
    real(kind=real64)            :: js                         ! Subdivision jacobian: xi->xip
    real(kind=real64)            :: xi                         ! Coordinate xi [xi_s(1,1),xi_s(1,2)]
    real(kind=real64)            :: aux(10)                    ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)           ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi(e%n_gnodes)       ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: pphi(e%n_pnodes)           ! Functional shape functions values
    real(kind=real64)            :: sphi(e%n_snodes)           ! Functional shape functions values
    real(kind=real64)            :: x(2)                       ! Position vector at xi
    real(kind=real64)            :: T(2)                       ! Tangent vector at xi
    real(kind=real64)            :: N(2)                       ! Normal vector at xi
    real(kind=real64)            :: rv(2)                      ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)            :: r, d1r, logd1r             ! Distance vector module, its inverse and log(1/r)
    real(kind=real64)            :: drdx(2)                    ! Distance vector derivatives with respect to x_k
    real(kind=real64)            :: jg                         ! Geometric jacobian
    real(kind=real64)            :: drdn                       ! Partial derivative of r respect to unit normal
    real(kind=real64)            :: jw                         ! Jacobians * weight
    real(kind=real64)            :: pphijw(e%n_pnodes)         ! primary shape functions * jw
    real(kind=real64)            :: sphijw(e%n_snodes)         ! secondary shape functions * jw
    real(kind=real64)            :: cteu1, cteu2, ctet1, ctet2 ! Auxiliary constants
    real(kind=real64)            :: fs_u, fs_t                 ! Fundamental solutions values
    ! Initialization
    h=0.d0
    g=0.d0
    cteu1=1.0d0/(8.0d0*c_pi*mu*(1.0d0-nu))
    cteu2=3.0d0-4.0d0*nu
    ctet1=-1.0d0/(4.0d0*c_pi*(1.0d0-nu))
    ctet2=1.0d0-2.0d0*nu
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
      ! Distance and functions of distance and unit normal
      r=sqrt(dot_product(rv,rv))
      d1r=1.d0/r
      logd1r=dlog(d1r)
      drdx=rv*d1r
      drdn=dot_product(drdx,n)
      ! Jacobians * weight
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
      ! Add integration points
      do l=1,2
        do k=1,2
          fs_u=cteu2*logd1r*c_dkr(l,k)+drdx(l)*drdx(k)
          fs_t=d1r*(drdn*(ctet2*c_dkr(l,k)+2.d0*drdx(l)*drdx(k))+ctet2*(n(l)*drdx(k)-n(k)*drdx(l)))
          h(:,l,k)=h(:,l,k)+fs_t*pphijw
          g(:,l,k)=g(:,l,k)+fs_u*sphijw
        end do
      end do
    end do
    ! Multiply by constants
    h=ctet1*h
    g=cteu1*g
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_staela2d_sbie_ext_st

  !! This subroutine calculates adaptatively the kernels for SBIE exterior integration using Telles transformation and subdivision
  !! if needed.
  recursive subroutine fbem_bem_staela2d_sbie_ext_adp(e,reverse,xi_s,x_i,mu,nu,qsp,ks,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                 !! Element
    logical                  :: reverse           !! Reverse orientation
    real(kind=real64)        :: xi_s(1,2)         !! Subdivision of the parent element
    real(kind=real64)        :: x_i(2)            !! Collocation point position vector
    real(kind=real64)        :: mu                !! Shear modulus
    real(kind=real64)        :: nu                !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp               !! Quasi-singular integration parameters
    integer                  :: ks                !! Current level of subdivisions
    integer                  :: ns                !! Maximum level of subdivision
    real(kind=real64)        :: h(e%n_pnodes,2,2) !! h integration kernels matrix
    real(kind=real64)        :: g(e%n_snodes,2,2) !! g integration kernels matrix
    ! Local
    integer           :: gln_near                 ! 1D Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer           :: gln                      ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical           :: subdivide                ! True if subdivision has to be performed
    real(kind=real64) :: barxi(1)                 ! Nearest element coordinate with respect to collocation point
    real(kind=real64) :: barxip(1)                ! Nearest element subdivision local coordinate with respect to collocation point
    real(kind=real64) :: rmin                     ! Minimum distance between collocation point and barxi on the element
    real(kind=real64) :: barr                     ! Telles jacobian at barxi
    real(kind=real64) :: cl                       ! Characteristic length
    real(kind=real64) :: d                        ! Normalized distance between collocation point and element subdivision
    integer           :: method                   ! Method used in nearest point algorithm
    real(kind=real64) :: tmp_xi_s(1,2)            ! Subdivision
    real(kind=real64) :: x_s(2,e%n_gnodes)        ! Coordinates of the element subdivision
    real(kind=real64) :: h_tmp(e%n_pnodes,2,2)    ! h integration kernels matrix (temporary)
    real(kind=real64) :: g_tmp(e%n_snodes,2,2)    ! g integration kernels matrix (temporary)
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
      cl=fbem_characteristic_length(2,e%gtype,x_s,1.d-12)
      call fbem_nearest_element_point_bem(2,e%gtype,x_s,cl,x_i,barxip,rmin,d,method)
    end if
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,4,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_staela2d_sbie_ext_adp',ns,'maximum number of subdivisions reached')
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
      call fbem_bem_staela2d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,h,g)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_staela2d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,h,g)
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_staela2d_sbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,mu,nu,gln,h_tmp,g_tmp)
      h=h+h_tmp
      g=g+g_tmp
    end if
  end subroutine fbem_bem_staela2d_sbie_ext_adp

!  !! This subroutine calculates the kernels for SBIE interior integration, needing only basic data.
!  subroutine fbem_bem_staela2d_sbie_int(ngp,type_g,type_f1,type_f2,delta_f,x_nodes,reverse,xi_i,mu,nu,h,g)
!    implicit none
!    ! I/O
!    integer                      :: ngp                             !! Number of Gauss point to be used (<=32)
!    integer                      :: type_g                          !! Geometrial interpolation
!    integer                      :: type_f1                         !! Functional interpolation (primary variables)
!    integer                      :: type_f2                         !! Functional interpolation (secondary variables)
!    real(kind=real64)            :: delta_f                         !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
!    real(kind=real64)            :: x_nodes(2,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
!    logical                      :: reverse                         !! Reverse normal vector
!    real(kind=real64)            :: xi_i                            !! Reference coordinate of the singular point.
!    real(kind=real64)            :: mu                              !! Shear modulus
!    real(kind=real64)            :: nu                              !! Poisson's ratio
!    real(kind=real64)            :: h(fbem_n_nodes(type_f1),2,2)    !! h kernel vector
!    real(kind=real64)            :: g(fbem_n_nodes(type_f2),2,2)    !! g kernel vector
!    ! Local
!    integer                      :: l, k                            ! Counter variables for load direction and observation direction
!    integer                      :: nnodes_g                        ! Number of nodes of the element interpolation
!    integer                      :: kphi                            ! Counter variable for shape functions loops
!    integer                      :: kip                             ! Counter of integration points
!    real(kind=real64)            :: x_i(2)                          ! Real coordinates of collocation point
!    real(kind=real64)            :: phi_f1(fbem_n_nodes(type_f1))   ! Functional shape functions values at xi
!    real(kind=real64)            :: phi_f2(fbem_n_nodes(type_f2))   ! Functional shape functions values at xi
!    real(kind=real64)            :: phi_f1_i(fbem_n_nodes(type_f1)) ! Functional shape functions values at xi_i
!    real(kind=real64)            :: phi_g(fbem_n_nodes(type_g))     ! Geometrical shape functions values at xi
!    real(kind=real64)            :: dphidxi_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi
!    integer                      :: nsub                            ! Number of subdivision of the element
!    integer                      :: ksub                            ! Counter of subdivision
!    real(kind=real64)            :: gamma                           ! Coordinate gamma
!    real(kind=real64)            :: w                               ! Weights of each integration point
!    real(kind=real64)            :: xip                             ! Coordinate  xip of subdivided element [0,1]
!    real(kind=real64)            :: js                              ! Jacobian of the xi [xisub(1,:),xisub[2,:]] -> xip [0,1] transformation
!    real(kind=real64)            :: xip_i(2)                        ! Singular point in xip space
!    real(kind=real64)            :: xi                              ! Coordinate xi
!    real(kind=real64)            :: xisub(2,2)                      ! Coordinates of subdivisions in xi space
!    real(kind=real64)            :: aux(10)                         ! Auxiliary variable needed for shape_functions module resources
!    real(kind=real64)            :: x(2)                            ! Position vector at xi
!    real(kind=real64)            :: T(2)                            ! Tangent vector at xi
!    real(kind=real64)            :: N(2)                            ! Normal vector at xi
!    real(kind=real64)            :: rv(2)                           ! Distance vector between collocation point and integration point (x-x_i)
!    real(kind=real64)            :: r, d1r, logd1r                  ! Distance vector module, its inverse and log(1/r)
!    real(kind=real64)            :: ra, rb                          ! Distance vector from collocation point to element vertices
!    real(kind=real64)            :: drdx(2)                         ! Distance vector derivatives with respect to x_k
!    real(kind=real64)            :: jg                              ! Geometric jacobian
!    real(kind=real64)            :: drdn                            ! Partial derivative of r respect to unit normal
!    real(kind=real64)            :: drdt                            ! Partial derivative of r respect to unit tangent
!    real(kind=real64)            :: jw                              ! Jacobians * weight
!    real(kind=real64)            :: phif1jw(fbem_n_nodes(type_f1))  ! shape functions * jw
!    real(kind=real64)            :: phif2jw(fbem_n_nodes(type_f2))  ! shape functions * jw
!    type(fbem_telles_parameters) :: telles_parameters               ! Telles parameters
!    real(kind=real64)            :: jt                              ! Telles jacobian
!    real(kind=real64)            :: cteu1, cteu2, ctep1, ctep2      ! Auxiliary constants
!    real(kind=real64)            :: fs_u                            ! Fundamental solutions values and auxiliary parts
!    real(kind=real64)            :: fs_p1(fbem_n_nodes(type_f1))    ! Fundamental solutions values and auxiliary parts
!    real(kind=real64)            :: fs_p2(fbem_n_nodes(type_f1))    ! Fundamental solutions values and auxiliary parts
!    !
!    ! Initialize
!    !
!    ! Number of nodes of the element
!    nnodes_g=fbem_n_nodes(type_g)
!    ! Initialize kernel matrices
!    h=0.d0
!    g=0.d0
!    ! Initialize auxiliary constants for fundamental solutions calculation
!    cteu1=1.0d0/(8.0d0*c_pi*mu*(1.0d0-nu))
!    cteu2=3.0d0-4.0d0*nu
!    ctep1=-1.0d0/(4.0d0*c_pi*(1.0d0-nu))
!    ctep2=1.0d0-2.0d0*nu
!    ! If xi_i belongs to one of the vertices, subdivision is not needed
!    if (fbem_check_xi_vertex(xi_i).eqv.(.true.)) then
!      nsub=1
!    ! If xi_i is inside the element, 2 subdivisions are needed
!    else
!      nsub=2
!    end if
!    ! Setup the subdivisions transformation from xi to xip: [xisub1,xisub2] -> [0,1]
!    select case (nsub)
!      ! If just 1 subdivision
!      case (1)
!        ! Coodinates xi of the subdivision
!        xisub(1,1)=-1.0d0
!        xisub(2,1)= 1.0d0
!        ! Coordinate xip of the collocation point
!        ! If xi_i is at xi=-1
!        if (xi_i.lt.0.0d0) xip_i(1)=0.0d0
!        ! If xi_i is at xi=1
!        if (xi_i.gt.0.0d0) xip_i(1)=1.0d0
!      ! If 2 subdivisions
!      case (2)
!        ! Coordinates xi of the subdivision 1
!        xisub(1,1)=-1.0d0
!        xisub(2,1)=xi_i
!        ! Coordinate xip of the collocation point
!        xip_i(1)=1.0d0
!        ! Coordinates xi of the subdivision 2
!        xisub(1,2)=xi_i
!        xisub(2,2)=1.0d0
!        ! Coordinate xip of the collocation point
!        xip_i(2)=0.0d0
!    end select
!    ! Calculate real coordinates of collocation point
!    ! Geometrical shape functions at xi_i
!#   define etype type_g
!#   define delta 0.0d0
!#   define xi xi_i
!#   define phi phi_g
!#   include <phi_1d.rc>
!#   undef etype
!#   undef delta
!#   undef xi
!#   undef phi
!    x_i=0.d0
!    do kphi=1,nnodes_g
!      x_i=x_i+phi_g(kphi)*x_nodes(:,kphi)
!    end do
!    ! Functional shape functions at xi_i
!#   define etype type_f1
!#   define delta delta_f
!#   define xi xi_i
!#   define phi phi_f1_i
!#   include <phi_1d.rc>
!#   undef etype
!#   undef delta
!#   undef xi
!#   undef phi
!    !
!    ! Weakly singular part
!    !
!    ! Loop through subdivisions
!    do ksub=1,nsub
!      ! Telles transformation parameters (barr=0 at the collocation point)
!      telles_parameters=fbem_telles01_calculate_parameters(xip_i(ksub),0.0d0)
!      ! Jacobian of xip->xi transformation (is constant)
!      js=xisub(2,ksub)-xisub(1,ksub)
!      ! Loop through integration points
!      do kip=1,gl01_n(ngp)
!        ! GAMMA->XIP->XI COORDINATE TRANSFORMATION
!        ! Gamma coordinate and weight
!        gamma=gl01_xi(kip,ngp)
!        w=gl01_w(kip,ngp)
!        ! xip coordinate, weight and jacobian from Telles transformation
!        call fbem_telles_xi_and_jacobian(telles_parameters,gamma,xip,jt)
!        ! xi coordinate
!        xi=js*xip+xisub(1,ksub)
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
!        ! Normal vector
!        N(1)=T(2)
!        N(2)=-T(1)
!        ! Geometric jacobian
!        jg=dot_product(T,T)
!        jg=sqrt(jg)
!        ! Unit normal vector
!        n=N/jg
!        ! Unit tangent vector
!        t=T/jg
!        ! Distance vector
!        rv=x-x_i
!        ! Distance vector norm and its inverse
!        r=dot_product(rv,rv)
!        r=sqrt(r)
!        d1r=1.d0/r
!        logd1r=dlog(d1r)
!        ! r_{,k}
!        drdx=rv*d1r
!        ! dr/dn
!        drdn=dot_product(drdx,n)
!        ! dr/dGamma
!        drdt=dot_product(drdx,t)
!        ! Jacobian * weight
!        jw=jg*js*jt*w
!        ! FUNCTIONAL SHAPE FUNCTIONS
!        ! Functional shape functions (primary variables) at xi
!#       define etype type_f1
!#       define delta delta_f
!#       define phi phi_f1
!#       include <phi_1d.rc>
!#       undef etype
!#       undef delta
!#       undef phi
!        ! Functional shape functions (secondary variables) at xi
!#       define etype type_f2
!#       define delta delta_f
!#       define phi phi_f2
!#       include <phi_1d.rc>
!#       undef etype
!#       undef delta
!#       undef phi
!        ! Functional shape functions * jacobians* weights
!        phif1jw=phi_f1*jw
!        phif2jw=phi_f2*jw
!        ! For g integral kernel
!        ! Loop through load and observation directions
!        do l=1,2
!          do k=1,2
!            ! Fundamental solutions values without cteu1 and ctep1, respectively
!            fs_u=cteu2*logd1r*c_dkr(l,k)
!            ! Add to kernels
!            g(:,l,k)=g(:,l,k)+fs_u*phif2jw
!          end do
!        end do
!      end do ! Loop through integration points
!    end do ! Loop through subdivisions
!    !
!    ! Regular part
!    !
!    ! Loop through subdivisions
!    do ksub=1,nsub
!      ! Jacobian of xip->xi transformation (is constant)
!      js=xisub(2,ksub)-xisub(1,ksub)
!      ! Loop through integration points
!      do kip=1,gl01_n(ngp)
!        ! XIP coordinate and weight
!        xip=gl01_xi(kip,ngp)
!        w=gl01_w(kip,ngp)
!        ! XIP->XI COORDINATE TRANSFORMATION
!        ! xi coordinate
!        xi=js*xip+xisub(1,ksub)
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
!        ! Normal vector
!        N(1)=T(2)
!        N(2)=-T(1)
!        ! Geometric jacobian
!        jg=dot_product(T,T)
!        jg=sqrt(jg)
!        ! Unit normal vector
!        n=N/jg
!        ! Unit tangent vector
!        t=T/jg
!        ! Distance vector
!        rv=x-x_i
!        ! Distance vector norm and its inverse
!        r=dot_product(rv,rv)
!        r=sqrt(r)
!        d1r=1.d0/r
!        logd1r=dlog(d1r)
!        ! r_{,k}
!        drdx=rv*d1r
!        ! dr/dn
!        drdn=dot_product(drdx,n)
!        ! dr/dGamma
!        drdt=dot_product(drdx,t)
!        ! Jacobian * weight
!        jw=jg*js*w
!        ! FUNCTIONAL SHAPE FUNCTIONS
!        ! Functional shape functions (primary variables) at xi
!#       define etype type_f1
!#       define delta delta_f
!#       define phi phi_f1
!#       include <phi_1d.rc>
!#       undef etype
!#       undef delta
!#       undef phi
!        ! Functional shape functions (secondary variables) at xi
!#       define etype type_f2
!#       define delta delta_f
!#       define phi phi_f2
!#       include <phi_1d.rc>
!#       undef etype
!#       undef delta
!#       undef phi
!        ! Functional shape functions * jacobians* weights
!        phif1jw=phi_f1*jw
!        phif2jw=phi_f2*jw
!        ! For g integral kernel
!        ! Loop through load and observation directions
!        do l=1,2
!          do k=1,2
!            ! Fundamental solutions values without cteu1 and ctep1, respectively
!            fs_u=drdx(l)*drdx(k)
!            ! Add to kernels
!            g(:,l,k)=g(:,l,k)+fs_u*phif2jw
!          end do
!        end do
!        ! For h integral kernel
!        ! Common part of H11 and H22
!        fs_p1=d1r*drdn*phif1jw
!        ! H11
!        h(:,1,1)=h(:,1,1)+fs_p1*(ctep2+2.d0*drdx(1)*drdx(1))
!        ! H22
!        h(:,2,2)=h(:,2,2)+fs_p1*(ctep2+2.d0*drdx(2)*drdx(2))
!        ! Common terms of H12 and H21
!        fs_p1=2.d0*d1r*drdn*drdx(1)*drdx(2)*phif1jw
!        fs_p2=ctep2*d1r*drdt*(phi_f1-phi_f1_i)*jw
!        ! H12
!        ! A part
!        h(:,1,2)=h(:,1,2)+fs_p1
!        ! B numerical part
!        h(:,1,2)=h(:,1,2)+fs_p2
!        ! H21
!        ! A part
!        h(:,2,1)=h(:,2,1)+fs_p1
!        ! B numerical part
!        h(:,2,1)=h(:,2,1)-fs_p2
!      end do ! Loop through integration points
!    end do ! Loop through subdivisions
!    ! Add analytical terms to H12 and H21
!    ! Calculate ra and rb
!    ra=sqrt((x_nodes(1,1)-x_i(1))**2+(x_nodes(2,1)-x_i(2))**2)
!    rb=sqrt((x_nodes(1,2)-x_i(1))**2+(x_nodes(2,2)-x_i(2))**2)
!    ! The addition depends on the position of xi_i
!    select case (nsub)
!      ! If just 1 subdivision, xi_i=-1 or xi_i=1
!      case (1)
!        ! If xi_i is at xi=-1
!        if (xi_i.lt.0.0d0) then
!          h(:,1,2)=h(:,1,2)+ctep2*phi_f1_i*log(rb)
!          h(:,2,1)=h(:,2,1)-ctep2*phi_f1_i*log(rb)
!        end if
!        ! If xi_i is at xi=1
!        if (xi_i.gt.0.0d0) then
!          h(:,1,2)=h(:,1,2)-ctep2*phi_f1_i*log(ra)
!          h(:,2,1)=h(:,2,1)+ctep2*phi_f1_i*log(ra)
!        end if
!      ! If 2 subdivisions, xi_i is inside de element -1<xi_i<1
!      case (2)
!        h(:,1,2)=h(:,1,2)+ctep2*phi_f1_i*(log(rb)-log(ra))
!        h(:,2,1)=h(:,2,1)-ctep2*phi_f1_i*(log(rb)-log(ra))
!    end select
!    ! Multiply h by ctep1 and g by cteu1
!    h=ctep1*h
!    g=cteu1*g
!    ! If the normal has to be reversed, then h=-h
!    if (reverse) h=-h
!  end subroutine fbem_bem_staela2d_sbie_int

  !! This subroutine calculates the kernels for SBIE interior integration, needing only basic data.
  subroutine fbem_bem_staela2d_sbie_int(e,reverse,xi_i,mu,nu,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element) :: e                          !! Integration element
    logical                :: reverse                    !! Reverse normal vector
    real(kind=real64)      :: xi_i(1)                    !! Local coordinate of the singular point.
    real(kind=real64)      :: mu                         !! Shear modulus
    real(kind=real64)      :: nu                         !! Poisson's ratio
    real(kind=real64)      :: h(e%n_pnodes,2,2)          !! h kernel vector
    real(kind=real64)      :: g(e%n_snodes,2,2)          !! g kernel vector
    ! Local
    integer                :: gln                        ! Gauss-Legendre number of integration points (<=32)
    integer                :: l, k                       ! Counter variables for load direction and observation direction
    integer                :: kphi                       ! Counter variable for shape functions loops
    integer                :: kip                        ! Counter of integration points
    real(kind=real64)      :: x_i(2)                     ! Real coordinates of collocation point
    real(kind=real64)      :: gphi(e%n_gnodes)           ! Geometrical shape functions values at xi
    real(kind=real64)      :: dgphidxi(e%n_gnodes)       ! Geometrical shape functions derivatives values at xi
    real(kind=real64)      :: pphi(e%n_pnodes)           ! Functional shape functions values at xi
    real(kind=real64)      :: sphi(e%n_snodes)           ! Functional shape functions values at xi
    real(kind=real64)      :: pphi_i(e%n_pnodes)         ! Functional shape functions values at xi_i
    integer                :: nsub                       ! Number of subdivision of the element
    integer                :: ksub                       ! Counter of subdivision
    real(kind=real64)      :: w                          ! Weights of each integration point
    real(kind=real64)      :: xip                        ! Coordinate  xip of subdivided element [0,1]
    real(kind=real64)      :: js                         ! Jacobian of the xi [xisub(1,:),xisub[2,:]] -> xip [0,1] transformation
    real(kind=real64)      :: xip_i(2)                   ! Singular point in xip space
    real(kind=real64)      :: xi                         ! Coordinate xi
    real(kind=real64)      :: xisub(2,2)                 ! Coordinates of element subdivisions
    real(kind=real64)      :: aux(10)                    ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)      :: x(2)                       ! Position vector at xi
    real(kind=real64)      :: T(2)                       ! Tangent vector at xi
    real(kind=real64)      :: N(2)                       ! Normal vector at xi
    real(kind=real64)      :: rv(2)                      ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)      :: r, d1r, logd1r             ! Distance vector module, its inverse and log(1/r)
    real(kind=real64)      :: ra, rb                     ! Distance vector from collocation point to element vertices
    real(kind=real64)      :: drdx(2)                    ! Distance vector derivatives with respect to x_k
    real(kind=real64)      :: jg                         ! Geometric jacobian
    real(kind=real64)      :: drdn                       ! Partial derivative of r respect to unit normal
    real(kind=real64)      :: drdt                       ! Partial derivative of r respect to unit tangent
    real(kind=real64)      :: jw                         ! Jacobians * weight
    real(kind=real64)      :: pphijw(e%n_pnodes)         ! phi^p * jw
    real(kind=real64)      :: sphijw(e%n_snodes)         ! phi^s * jw
    real(kind=real64)      :: lsphi(e%n_snodes)          ! log(1/r)*phi^s integrals
    real(kind=real64)      :: cteu1, cteu2, ctet1, ctet2 ! Auxiliary constants
    real(kind=real64)      :: fs_u                       ! Fundamental solutions values and auxiliary parts
    real(kind=real64)      :: fs_t                       ! Fundamental solutions values and auxiliary parts
    ! Integration points to be used
    gln=32
    ! Initialize kernels
    h=0.d0
    g=0.d0
    ! Initialize auxiliary constants for fundamental solutions calculation
    cteu1=1.0d0/(8.0d0*c_pi*mu*(1.0d0-nu))
    cteu2=3.0d0-4.0d0*nu
    ctet1=-1.0d0/(4.0d0*c_pi*(1.0d0-nu))
    ctet2=1.0d0-2.0d0*nu
    ! Setup the subdivisions for transformation xip -> xi: [0,1] -> [xisub_1,xisub_2]
    if (fbem_check_xi_vertex(xi_i(1))) then
      nsub=1
      xisub(1,1)=-1.0d0
      xisub(2,1)= 1.0d0
      if (xi_i(1).lt.0.0d0) xip_i(1)=0.0d0
      if (xi_i(1).gt.0.0d0) xip_i(1)=1.0d0
    else
      nsub=2
      xisub(1,1)=-1.0d0
      xisub(2,1)=xi_i(1)
      xip_i(1)=1.0d0
      xisub(1,2)=xi_i(1)
      xisub(2,2)=1.0d0
      xip_i(2)=0.0d0
    end if
    ! Calculate spatial coordinates at xi_i
#   define etype e%gtype
#   define delta 0.d0
#   define xi xi_i(1)
#   define phi gphi
#   include <phi_1d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
    x_i=0.d0
    do kphi=1,e%n_gnodes
      x_i=x_i+gphi(kphi)*e%x(:,kphi)
    end do
    ! Save functional shape functions at xi_i
#   define etype e%ptype
#   define delta e%ptype_delta
#   define xi xi_i(1)
#   define phi pphi_i
#   include <phi_1d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
    !
    ! Weakly singular integrals
    !
    call fbem_bem_logd1rphi_1d_int(2,e%gtype,e%x,e%stype,e%stype_delta,xi_i,lsphi)
    g(:,1,1)=cteu2*lsphi
    g(:,2,2)=cteu2*lsphi
    !
    ! Regular integrals
    !
    do ksub=1,nsub
      ! Numerical integration
      do kip=1,gl01_n(gln)
        ! XIP COORDINATE
        xip=gl01_xi(kip,gln)
        w=gl01_w(kip,gln)
        ! XIP->XI TRANSFORMATION
        js=xisub(2,ksub)-xisub(1,ksub)
        xi=js*xip+xisub(1,ksub)
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
        ! Normal vector
        N(1)=T(2)
        N(2)=-T(1)
        ! xi->x jacobian
        jg=sqrt(dot_product(T,T))
        ! Unit normal vector
        n=N/jg
        ! Unit tangent vector
        t=T/jg
        ! Distance vector and other derived terms
        rv=x-x_i
        r=dot_product(rv,rv)
        r=sqrt(r)
        d1r=1.d0/r
        logd1r=log(d1r)
        drdx=rv*d1r
        drdn=dot_product(drdx,n)
        drdt=dot_product(drdx,t)
        ! Jacobian * weight
        jw=jg*js*w
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
        ! Functional shape functions * jacobians* weights
        pphijw=pphi*jw
        sphijw=sphi*jw
        ! Add integration points (regular part)
        do l=1,2
          do k=1,2
            fs_u=drdx(l)*drdx(k)
            fs_t=d1r*(drdn*(ctet2*c_dkr(l,k)+2.d0*drdx(l)*drdx(k)))
            h(:,l,k)=h(:,l,k)+fs_t*pphijw
            g(:,l,k)=g(:,l,k)+fs_u*sphijw
          end do
        end do
        ! Regularized part of the CPV integrals
        h(:,1,2)=h(:,1,2)+ctet2*d1r*drdt*(pphi-pphi_i)*jw
        h(:,2,1)=h(:,2,1)-ctet2*d1r*drdt*(pphi-pphi_i)*jw
      end do
    end do
    !
    ! Analytical integrals
    !
    ra=sqrt((e%x(1,1)-x_i(1))**2+(e%x(2,1)-x_i(2))**2)
    rb=sqrt((e%x(1,2)-x_i(1))**2+(e%x(2,2)-x_i(2))**2)
    select case (nsub)
      case (1)
        if (xi_i(1).lt.0.0d0) then
          h(:,1,2)=h(:,1,2)+ctet2*pphi_i*log(rb)
          h(:,2,1)=h(:,2,1)-ctet2*pphi_i*log(rb)
        end if
        if (xi_i(1).gt.0.0d0) then
          h(:,1,2)=h(:,1,2)-ctet2*pphi_i*log(ra)
          h(:,2,1)=h(:,2,1)+ctet2*pphi_i*log(ra)
        end if
      case (2)
        h(:,1,2)=h(:,1,2)+ctet2*pphi_i*(log(rb)-log(ra))
        h(:,2,1)=h(:,2,1)-ctet2*pphi_i*(log(rb)-log(ra))
    end select
    ! Multiply by constants
    h=ctet1*h
    g=cteu1*g
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_staela2d_sbie_int

  ! =====================
  ! BE BODY LOAD ELEMENTS
  ! =====================

  !! Efficient automatic integration of BE body load elements
  subroutine fbem_bem_staela2d_sbie_bl_auto(e,x_i,mu,nu,qsp,ns,gbl)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                   !! Integration element
    real(kind=real64)        :: x_i(2)              !! Collocation point
    real(kind=real64)        :: mu                  !! Shear modulus
    real(kind=real64)        :: nu                  !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp                 !! Quasi-singular integration parameters
    integer                  :: ns                  !! Maximum level of subdivisions
    real(kind=real64)        :: gbl(e%n_snodes,2,2) !! gbl integration kernel
    ! Local
    real(kind=real64)        :: r(2)                               ! Distance vector
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
        call fbem_error_message(output_unit,0,'fbem_bem_staela2d_sbie_bl_auto',0,'it is not possible to collocate at a point load')
      else
        call fbem_bem_staela2d_sbie_u(e%x(:,1),x_i,mu,nu,gbl(1,:,:))
        return
      end if
    ! LINE OR SURFACE BODY LOAD
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
        call fbem_nearest_element_point_bem(2,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
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
        call fbem_bem_staela2d_sbie_bl_int(e,barxi,mu,nu,gbl)
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
          call fbem_bem_staela2d_sbie_bl_ext_pre(ps,e,x_i,mu,nu,gbl)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_staela2d_sbie_bl_ext_adp(e,xi_s,x_i,mu,nu,qsp,1,ns,gbl)
        end if
    end select
  end subroutine fbem_bem_staela2d_sbie_bl_auto

  !! This subroutine calculates the body load kernels for SBIE interior integration, needing only basic data.
  subroutine fbem_bem_staela2d_sbie_bl_int(e,xi_i,mu,nu,gbl)
    implicit none
    ! I/O
    type(fbem_bem_element) :: e                                !! Integration element
    real(kind=real64)      :: xi_i(e%d)                        !! Local coordinate of the singular point.
    real(kind=real64)      :: mu                               !! Shear modulus
    real(kind=real64)      :: nu                               !! Poisson's ratio
    real(kind=real64)      :: gbl(e%n_snodes,2,2)              !! gbl kernel vector
    ! Local
    integer                :: gln                              ! 1D Gauss-Legendre number of integration points (<=32)
    integer                :: l, k                             ! Counter variables for load direction and observation direction
    integer                :: kphi                             ! Counter variable for shape functions loops
    integer                :: kip                              ! Counter of integration points
    real(kind=real64)      :: x_i(2)                           ! Real coordinates of collocation point
    real(kind=real64)      :: gphi(e%n_gnodes)                 ! Geometrical shape functions values at xi
    real(kind=real64)      :: dgphidxi(e%n_gnodes)             ! Geometrical shape functions derivatives values at xi
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
    integer                :: nsub                             ! Number of subdivision of the element
    integer                :: ksub                             ! Counter of subdivision
    real(kind=real64)      :: w                                ! Weights of each integration point
    real(kind=real64)      :: xip                              ! Coordinate  xip of subdivided element [0,1]
    real(kind=real64)      :: js                               ! Jacobian of the xi [xisub(1,:),xisub[2,:]] -> xip [0,1] transformation
    real(kind=real64)      :: xip_i(2)                         ! Singular point in xip space
    real(kind=real64)      :: xin(e%d)                         ! Coordinate xi
    real(kind=real64)      :: xisub(2,2)                       ! Coordinates of element subdivisions
    real(kind=real64)      :: aux(10)                          ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)      :: x(2)                             ! Position vector at xi
    real(kind=real64)      :: T(2), T1(2), T2(2)               ! Tangent vectors at xi
    real(kind=real64)      :: rv(2)                            ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)      :: r, d1r, logd1r                   ! Distance vector module, its inverse and log(1/r)
    real(kind=real64)      :: drdx(2)                          ! Distance vector derivatives with respect to x_k
    real(kind=real64)      :: jg                               ! Geometric jacobian
    real(kind=real64)      :: jw                               ! Jacobians * weight
    real(kind=real64)      :: sphijw(e%n_snodes)               ! phi^s * jw
    real(kind=real64)      :: lsphi(e%n_snodes)                ! log(1/r)*phi^s integrals
    real(kind=real64)      :: cteu1, cteu2                     ! Auxiliary constants
    real(kind=real64)      :: fs_u                             ! Fundamental solutions values and auxiliary parts
    ! Initialize kernels
    gbl=0.d0
    ! Initialize auxiliary constants for fundamental solution calculation
    cteu1=1.0d0/(8.0d0*c_pi*mu*(1.0d0-nu))
    cteu2=3.0d0-4.0d0*nu
    ! Calculate
    select case (e%d)
      !
      ! LINE LOAD
      !
      case (1)
        ! Integration points to be used
        gln=30
        ! Setup the subdivisions for transformation xip -> xi: [0,1] -> [xisub_1,xisub_2]
        if (fbem_check_xi_vertex(xi_i(1))) then
          nsub=1
          xisub(1,1)=-1.d0
          xisub(2,1)= 1.d0
          if (xi_i(1).lt.0.d0) xip_i(1)=0.d0
          if (xi_i(1).gt.0.d0) xip_i(1)=1.d0
        else
          nsub=2
          xisub(1,1)=-1.d0
          xisub(2,1)=xi_i(1)
          xip_i(1)=1.d0
          xisub(1,2)=xi_i(1)
          xisub(2,2)=1.d0
          xip_i(2)=0.d0
        end if
        ! Calculate spatial coordinates at xi_i
#       define etype e%gtype
#       define delta 0.d0
#       define xi xi_i(1)
#       define phi gphi
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef xi
#       undef phi
        x_i=0.d0
        do kphi=1,e%n_gnodes
          x_i=x_i+gphi(kphi)*e%x(:,kphi)
        end do
        ! Weakly singular integrals
        call fbem_bem_logd1rphi_1d_int(2,e%gtype,e%x,e%stype,e%stype_delta,xi_i,lsphi)
        gbl(:,1,1)=cteu2*lsphi
        gbl(:,2,2)=cteu2*lsphi
        ! Regular integrals
        do ksub=1,nsub
          ! Numerical integration
          do kip=1,gl01_n(gln)
            ! XIP COORDINATE
            xip=gl01_xi(kip,gln)
            w=gl01_w(kip,gln)
            ! XIP->XI TRANSFORMATION
            js=xisub(2,ksub)-xisub(1,ksub)
            xin(1)=js*xip+xisub(1,ksub)
            ! XI->X TRANSFORMATION
#           define etype e%gtype
#           define delta 0.d0
#           define xi xin(1)
#           define phi gphi
#           define dphidxi dgphidxi
#           include <phi_and_dphidxi_1d.rc>
#           undef etype
#           undef delta
#           undef xi
#           undef phi
#           undef dphidxi
            x=0.d0
            T=0.d0
            do kphi=1,e%n_gnodes
              x=x+gphi(kphi)*e%x(:,kphi)
              T=T+dgphidxi(kphi)*e%x(:,kphi)
            end do
            ! xi->x jacobian
            jg=sqrt(dot_product(T,T))
            ! Distance vector and other derived terms
            rv=x-x_i
            r=sqrt(dot_product(rv,rv))
            d1r=1.d0/r
            logd1r=log(d1r)
            drdx=rv*d1r
            ! Jacobian * weight
            jw=jg*js*w
            ! FUNCTIONAL SHAPE FUNCTIONS
            ! Functional shape functions (body load variables) at xi
#           define etype e%stype
#           define delta e%stype_delta
#           define xi xin(1)
#           define phi sphi
#           include <phi_1d.rc>
#           undef etype
#           undef delta
#           undef xi
#           undef phi
            ! Functional shape functions * jacobians* weights
            sphijw=sphi*jw
            ! Add integration points
            do l=1,2
              do k=1,2
                gbl(:,l,k)=gbl(:,l,k)+drdx(l)*drdx(k)*sphijw
              end do
            end do
          end do
        end do
      !
      ! SURFACE LOAD
      !
      case (2)
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
              xin(1)=xi_i(1)+rho*costheta
              xin(2)=xi_i(2)+rho*sintheta
              ! XI->X TRANSFORMATION
#             define etype e%gtype
#             define delta 0.d0
#             define xi xin
#             define phi gphi
#             define dphidxi1 dgphidxi1
#             define dphidxi2 dgphidxi2
#             include <phi_and_dphidxik_2d.rc>
#             undef etype
#             undef delta
#             undef xi
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
              ! xi->x jacobian
              jg=T1(1)*T2(2)-T1(2)*T2(1)
              ! Distance vector and other derived terms
              rv=x-x_i
              r=dot_product(rv,rv)
              r=sqrt(r)
              d1r=1.d0/r
              logd1r=dlog(d1r)
              drdx=rv*d1r
              ! Jacobians * weights
              jw=jg*rho*jthetap*w_angular*w_radial
              ! FUNCTIONAL SHAPE FUNCTIONS
              ! Functional shape functions (body load variables) at xi
#             define etype e%stype
#             define delta e%stype_delta
#             define xi xin
#             define phi sphi
#             include <phi_2d.rc>
#             undef etype
#             undef delta
#             undef xi
#             undef phi
              ! Functional shape functions * jacobians * weights
              sphijw=sphi*jw
              ! Add integration points
              do l=1,2
                do k=1,2
                  fs_u=cteu2*logd1r*c_dkr(l,k)+drdx(l)*drdx(k)
                  gbl(:,l,k)=gbl(:,l,k)+fs_u*sphijw
                end do
              end do
            end do
          end do
        end do
      !
      ! OTHERS
      !
      case default
        call fbem_error_message(output_unit,0,'fbem_bem_staela2d_sbie_bl_int',0,&
                                'it is only possible to integrate line or surface loads')
    end select
    ! Multiply by constants
    gbl=cteu1*gbl
  end subroutine fbem_bem_staela2d_sbie_bl_int

  !! This subroutine calculates the body load kernels for SBIE exterior integration using precalculated data at integration points.
  subroutine fbem_bem_staela2d_sbie_bl_ext_pre(ps,e,x_i,mu,nu,gbl)
    implicit none
    ! I/O
    integer                :: ps                  !! Selected precalculated dataset
    type(fbem_bem_element) :: e                   !! Element
    real(kind=real64)      :: x_i(2)              !! Collocation point position vector
    real(kind=real64)      :: mu                  !! Shear modulus
    real(kind=real64)      :: nu                  !! Poisson's ratio
    real(kind=real64)      :: gbl(e%n_snodes,2,2) !! gbl kernels
    ! Local
    integer                :: l, k                ! Counter variables for load direction and observation direction
    integer                :: kip                 ! Counter variable for integration points loop
    real(kind=real64)      :: x(2)                ! Position vector at integration point
    real(kind=real64)      :: sphijw(e%n_snodes)  ! phi^s * jacobian * weight at integration point
    real(kind=real64)      :: rv(2)               ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)      :: r, d1r, logd1r      ! Distance vector module, inverse and log(1/r)
    real(kind=real64)      :: drdx(2)             ! Distance vector derivatives with respect to x_k
    real(kind=real64)      :: cteu1, cteu2        ! Auxiliary constants
    real(kind=real64)      :: fs_u                ! Fundamental solutions values
    ! Initialize kernels
    gbl=0.d0
    ! Initialize auxiliary constants for fundamental solutions calculation
    cteu1=1.0d0/(8.0d0*c_pi*mu*(1.0d0-nu))
    cteu2=3.0d0-4.0d0*nu
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      x=e%ps_x(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      d1r=1.d0/r
      logd1r=log(d1r)
      drdx=rv*d1r
      do l=1,2
        do k=1,2
          fs_u=cteu2*logd1r*c_dkr(l,k)+drdx(l)*drdx(k)
          gbl(:,l,k)=gbl(:,l,k)+fs_u*sphijw
        end do
      end do
    end do
    ! Multiply by constants
    gbl=cteu1*gbl
  end subroutine fbem_bem_staela2d_sbie_bl_ext_pre

  !! This subroutine calculates the body load kernels for SBIE exterior integration (near collocation points) using Telles
  !! transformation within a subdivision of the element, needing only basic data.
  subroutine fbem_bem_staela2d_sbie_bl_ext_st(e,xi_s,x_i,barxip,barr,mu,nu,gln,gbl)
    implicit none
    ! I/O
    type(fbem_bem_element)       :: e                                  !! Integration element
    real(kind=real64)            :: xi_s(e%d,fbem_n_vertices(e%gtype)) !! Coordinates of the subdivision of the parent element (xi space [-1,1])
    real(kind=real64)            :: x_i(2)                             !! Collocation point position vector
    real(kind=real64)            :: barxip(e%d)                        !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)            :: barr                               !! Telles jacobian at barxip
    real(kind=real64)            :: mu                                 !! Shear modulus
    real(kind=real64)            :: nu                                 !! Poisson's ratio
    integer                      :: gln                                !! 1D Gauss-Legendre number of integration points (<=32)
    real(kind=real64)            :: gbl(e%n_snodes,2,2)                !! gbl kernels
    ! Local
    integer                      :: l, k                               ! Counter variables for load direction and observation direction
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
    real(kind=real64)            :: x(2)                               ! Position vector at xi
    real(kind=real64)            :: T(2), T1(2), T2(2)                 ! Tangent vector at xi
    real(kind=real64)            :: rv(2)                              ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)            :: r, d1r, logd1r                     ! Distance vector module, its inverse and log(1/r)
    real(kind=real64)            :: drdx(2)                            ! Distance vector derivatives with respect to x_k
    real(kind=real64)            :: jg                                 ! Geometric jacobian
    real(kind=real64)            :: jw                                 ! Jacobians * weight
    real(kind=real64)            :: sphijw(e%n_snodes)                 ! secondary shape functions * jw
    real(kind=real64)            :: cteu1, cteu2                       ! Auxiliary constants
    real(kind=real64)            :: fs_u                               ! Fundamental solutions values
    ! Initialize kernels
    gbl=0.d0
    ! Initialize auxiliary constants for fundamental solutions calculation
    cteu1=1.0d0/(8.0d0*c_pi*mu*(1.0d0-nu))
    cteu2=3.0d0-4.0d0*nu
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
          logd1r=dlog(d1r)
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
          do l=1,2
            do k=1,2
              fs_u=cteu2*logd1r*c_dkr(l,k)+drdx(l)*drdx(k)
              gbl(:,l,k)=gbl(:,l,k)+fs_u*sphijw
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
                ! Geometric jacobian
                jg=T1(1)*T2(2)-T1(2)*T2(1)
                ! Distance vector and derived terms
                rv=x-x_i
                r=sqrt(dot_product(rv,rv))
                d1r=1.d0/r
                logd1r=dlog(d1r)
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
                do l=1,2
                  do k=1,2
                    fs_u=cteu2*logd1r*c_dkr(l,k)+drdx(l)*drdx(k)
                    gbl(:,l,k)=gbl(:,l,k)+fs_u*sphijw
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
                jg=T1(1)*T2(2)-T1(2)*T2(1)
                ! Distance vector and derived terms
                rv=x-x_i
                r=sqrt(dot_product(rv,rv))
                d1r=1.d0/r
                logd1r=dlog(d1r)
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
                do l=1,2
                  do k=1,2
                    fs_u=cteu2*logd1r*c_dkr(l,k)+drdx(l)*drdx(k)
                    gbl(:,l,k)=gbl(:,l,k)+fs_u*sphijw
                  end do
                end do
              end do
            end do
            case default
              call fbem_error_message(error_unit,0,__FILE__,__LINE__,'n_edges not valid')
        end select
      !
      ! OTHERS
      !
      case default
        call fbem_error_message(output_unit,0,'fbem_bem_staela2d_sbie_bl_int',0,&
                                'it is only possible to integrate line or surface loads')
    end select
    ! Multiply by constants
    gbl=cteu1*gbl
  end subroutine fbem_bem_staela2d_sbie_bl_ext_st

  !! This subroutine calculates adaptatively the body load kernels for SBIE exterior integration using Telles transformation and
  !! subdivision if needed.
  recursive subroutine fbem_bem_staela2d_sbie_bl_ext_adp(e,xi_s,x_i,mu,nu,qsp,ks,ns,gbl)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                      !! Element
    real(kind=real64)        :: xi_s(e%d,fbem_n_vertices(e%gtype))     !! Subdivision of the parent element
    real(kind=real64)        :: x_i(2)                                 !! Collocation point position vector
    real(kind=real64)        :: mu                                     !! Shear modulus
    real(kind=real64)        :: nu                                     !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp                                    !! Quasi-singular integration parameters
    integer                  :: ks                                     !! Current level of subdivisions
    integer                  :: ns                                     !! Maximum level of subdivision
    real(kind=real64)        :: gbl(e%n_snodes,2,2)                    !! gbl kernels
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
    real(kind=real64)        :: x_s(2,e%n_gnodes)                      ! Coordinates of the element subdivision
    real(kind=real64)        :: gbl_tmp(e%n_snodes,2,2)                ! gbl kernels (temporary)
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
      end select
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
        call fbem_warning_message(error_unit,0,'fbem_bem_staela2d_sbie_bl_ext_adp',ns,'maximum number of subdivisions reached')
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
          call fbem_bem_staela2d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
          ! SUBLINE 2
          tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
          tmp_xi_s(1,2)=xi_s(1,2)
          call fbem_bem_staela2d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
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
              call fbem_bem_staela2d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
              ! SUBTRI 2
              tmp_xi_s(:,1)=xi_s(:,2)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              call fbem_bem_staela2d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
              ! SUBTRI 3
              tmp_xi_s(:,1)=xi_s(:,3)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              call fbem_bem_staela2d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
              ! SUBTRI 4
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
              call fbem_bem_staela2d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
            ! QUADRILATERALS
            case (4)
              ! SUBQUAD 1
              tmp_xi_s(:,1)=xi_s(:,1)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,3)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,4)=0.50d0*(xi_s(:,1)+xi_s(:,4))
              call fbem_bem_staela2d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
              ! SUBQUAD 2
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,2)=xi_s(:,2)
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,4)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              call fbem_bem_staela2d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
              ! SUBQUAD 3
              tmp_xi_s(:,1)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=xi_s(:,3)
              tmp_xi_s(:,4)=0.50d0*(xi_s(:,3)+xi_s(:,4))
              call fbem_bem_staela2d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
              ! SUBQUAD 4
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,4))
              tmp_xi_s(:,2)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,4)=xi_s(:,4)
              call fbem_bem_staela2d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,gbl)
          end select
      end select
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_staela2d_sbie_bl_ext_st(e,xi_s,x_i,barxip,barr,mu,nu,gln,gbl_tmp)
      gbl=gbl+gbl_tmp
    end if
  end subroutine fbem_bem_staela2d_sbie_bl_ext_adp
  ! --------------------------------------------------------------------------------------------------------------------------------

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (HBIE)

  !! Fundamental solution d*
  subroutine fbem_bem_staela2d_hbie_d(x,x_i,n_i,nu,do)
    implicit none
    ! I/O
    real(kind=real64) :: x(2)    !! Observation point
    real(kind=real64) :: x_i(2)  !! Collocation point
    real(kind=real64) :: n_i(2)  !! Collocation point unit normal
    real(kind=real64) :: nu      !! Poisson's ratio
    real(kind=real64) :: do(2,2) !! d*_{lk}
    ! Local
    integer           :: il, ik
    real(kind=real64) :: cte1, cted
    real(kind=real64) :: rv(2), r, d1r, drdx(2), drdni
    cte1=1.0d0-2.0d0*nu
    cted=-1.0d0/(4.0d0*c_pi*(1.0d0-nu))
    rv=x-x_i
    r=dot_product(rv,rv)
    r=sqrt(r)
    d1r=1.d0/r
    drdx=rv*d1r
    drdni=-dot_product(drdx,n_i)
    do il=1,2
      do ik=1,2
        do(il,ik)=d1r*(drdni*(cte1*c_dkr(il,ik)+2.d0*drdx(il)*drdx(ik))+cte1*(drdx(ik)*n_i(il)-drdx(il)*n_i(ik)))
      end do
    end do
    do=cted*do
  end subroutine fbem_bem_staela2d_hbie_d

  !! Fundamental solution s*
  subroutine fbem_bem_staela2d_hbie_s(x,n,x_i,n_i,mu,nu,so)
    implicit none
    ! I/O
    real(kind=real64) :: x(2)    !! Observation point
    real(kind=real64) :: n(2)    !! Observation point unit normal
    real(kind=real64) :: x_i(2)  !! Collocation point
    real(kind=real64) :: n_i(2)  !! Collocation point unit normal
    real(kind=real64) :: mu      !! Shear modulus
    real(kind=real64) :: nu      !! Poisson's ratio
    real(kind=real64) :: so(2,2) !! s*_{lk}
    ! Local
    integer           :: il, ik
    real(kind=real64) :: cte1, cte2, ctes
    real(kind=real64) :: rv(2), r, d1r, d1r2, drdx(2), drdn, drdni, n_dot_ni
    cte1=1.0d0-2.0d0*nu
    cte2=1.0d0-4.0d0*nu
    ctes=mu/(2.0d0*c_pi*(1.0d0-nu))
    rv=x-x_i
    r=dot_product(rv,rv)
    r=sqrt(r)
    d1r=1.d0/r
    d1r2=d1r**2
    drdx=rv*d1r
    drdn=dot_product(drdx,n)
    drdni=-dot_product(drdx,n_i)
    n_dot_ni=dot_product(n,n_i)
    do il=1,2
      do ik=1,2
        so(il,ik)=d1r2*(2.d0*drdn*(drdni*(4.d0*drdx(il)*drdx(ik)-nu*c_dkr(il,ik))+nu*drdx(il)*n_i(ik)+cte1*drdx(ik)*n_i(il))&
                       -2.d0*drdni*(nu*drdx(ik)*n(il)+cte1*drdx(il)*n(ik))+n_dot_ni*(2.d0*nu*drdx(il)*drdx(ik)+cte1*c_dkr(il,ik))&
                       +cte1*n(il)*n_i(ik)-cte2*n(ik)*n_i(il))
      end do
    end do
    so=ctes*so
  end subroutine fbem_bem_staela2d_hbie_s

  !! Efficient automatic integration
  subroutine fbem_bem_staela2d_hbie_auto(e,reverse,x_i,n_i,mu,nu,qsp,ns,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                 !! Integration element
    logical                  :: reverse           !! Reverse orientation
    real(kind=real64)        :: x_i(2)            !! Collocation point
    real(kind=real64)        :: n_i(2)            !! Unit normal at the collocation point
    real(kind=real64)        :: mu                !! Shear modulus
    real(kind=real64)        :: nu                !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp               !! Quasi-singular integration parameters
    integer                  :: ns                !! Maximum level of subdivisions
    real(kind=real64)        :: m(e%n_pnodes,2,2) !! h integration kernel
    real(kind=real64)        :: l(e%n_snodes,2,2) !! g integration kernel
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
        call fbem_bem_staela2d_hbie_int(e,reverse,barxi,mu,nu,m,l)
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
          call fbem_bem_staela2d_hbie_ext_pre(ps,e,reverse,x_i,n_i,mu,nu,m,l)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_staela2d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,mu,nu,qsp,1,ns,m,l)
        end if
    end select
  end subroutine fbem_bem_staela2d_hbie_auto

!  !! This subroutine calculates the kernels for HBIE exterior integration (far collocation points), needing precalculated data at
!  !! integration points.
!  !!
!  !! It is useful because for a given element, the vast majority of kernels are associated with far collocation points, which
!  !! need the same number of Gauss points. To generate integrand components, the subroutine <tt>fbem_gauss_points_values2d</tt>
!  !! can be called (from geometry module).
!  subroutine fbem_bem_staela2d_hbie_ext_pre(rule,max_rule,ngp,max_ngp,type_f1,type_f2,phi_f1,phi_f2,x,n,jw,reverse,x_i,n_i,mu,nu,m,l)
!    implicit none
!    ! I/O
!    integer           :: rule                                           !! Selected rule
!    integer           :: max_rule                                       !! Total number of rules in the rules set-up
!    integer           :: ngp                                            !! Number of Gaussian point of the selected rule
!    integer           :: max_ngp                                        !! Maximum number of Gaussian points in the rules set-up
!    integer           :: type_f1                                        !! Functional interpolation (primary variables)
!    integer           :: type_f2                                        !! Geometrical interpolation (secondary variables)
!    real(kind=real64) :: phi_f1(fbem_n_nodes(type_f1),max_ngp,max_rule) !! Functional shape functions at integration points
!    real(kind=real64) :: phi_f2(fbem_n_nodes(type_f2),max_ngp,max_rule) !! Functional shape functions at integration points
!    real(kind=real64) :: x(2,max_ngp,max_rule)                          !! Position vectors at integration points
!    real(kind=real64) :: n(2,max_ngp,max_rule)                          !! Unit normal vector at integration points
!    real(kind=real64) :: jw(max_ngp,max_rule)                           !! Geometric jacobian multiplied by weight at integration points
!    logical           :: reverse                                        !! Reverse normal vector
!    real(kind=real64) :: x_i(2)                                         !! Collocation point position vector
!    real(kind=real64) :: n_i(2)                                         !! Collocation point unit normal vector
!    real(kind=real64) :: mu                                             !! Shear modulus
!    real(kind=real64) :: nu                                             !! Poisson's ratio
!    real(kind=real64) :: m(fbem_n_nodes(type_f1),2,2)                   !! m integration kernels matrix
!    real(kind=real64) :: l(fbem_n_nodes(type_f2),2,2)                   !! l integration kernels matrix
!    ! Local
!    integer           :: il, ik                                      ! Counter variables for load direction and observation direction
!    integer           :: kip                                         ! Counter variable for integration points loop
!    real(kind=real64) :: rv(2)                                       ! Distance vector between collocation point and integration point (x-x_i)
!    real(kind=real64) :: r, d1r, d1r2                                ! Distance vector module, inverse
!    real(kind=real64) :: drdx(2)                                     ! Distance vector derivatives with respect to x_k
!    real(kind=real64) :: drdn                                        ! Partial derivative of r respect to unit normal
!    real(kind=real64) :: drdni                                       ! Partial derivative of r respect to unit normal at collocation point
!    real(kind=real64) :: n_dot_ni                                    ! Dot product of unit normals
!    real(kind=real64) :: phif1jw(fbem_n_nodes(type_f1))              ! Auxiliary variables for integrand evaluation
!    real(kind=real64) :: phif2jw(fbem_n_nodes(type_f2))              ! Auxiliary variables for integrand evaluation
!    real(kind=real64) :: cte1, cte2, ctes, cted                      ! Auxiliary constants
!    real(kind=real64) :: fs_s, fs_d                                  ! Fundamental solutions values
!    !
!    ! Initialize
!    !
!    ! Initialize kernel matrices
!    m=0.d0
!    l=0.d0
!    ! Initialize auxiliary constants for fundamental solutions calculation
!    cte1=1.0d0-2.0d0*nu
!    cte2=1.0d0-4.0d0*nu
!    ctes=mu/(2.0d0*c_pi*(1.0d0-nu))
!    cted=-1.0d0/(4.0d0*c_pi*(1.0d0-nu))
!    !
!    ! Numerical integration
!    !
!    ! Loop through integration points
!    do kip=1,ngp
!      ! Distance vector
!      rv=x(:,kip,rule)-x_i
!      ! Distance vector norm and its inverse
!      r=dot_product(rv,rv)
!      r=sqrt(r)
!      d1r=1.d0/r
!      d1r2=d1r**2
!      ! r_{,k}
!      drdx=rv*d1r
!      ! dr/dn
!      drdn=dot_product(drdx,n(:,kip,rule))
!      ! dr/dni
!      drdni=-dot_product(drdx,n_i)
!      ! n · n_i
!      n_dot_ni=dot_product(n(:,kip,rule),n_i)
!      ! Shape functions * jw
!      phif1jw=phi_f1(:,kip,rule)*jw(kip,rule)
!      phif2jw=phi_f2(:,kip,rule)*jw(kip,rule)
!      ! Loop through load and observation direction
!      do il=1,2
!        do ik=1,2
!          ! Fundamental solutions values without ctes and cted, respectively
!          fs_s=d1r2*(2.d0*drdn*(drdni*(4.d0*drdx(il)*drdx(ik)-nu*c_dkr(il,ik))+nu*drdx(il)*n_i(ik)+cte1*drdx(ik)*n_i(il))&
!                    -2.d0*drdni*(nu*drdx(ik)*n(il,kip,rule)+cte1*drdx(il)*n(ik,kip,rule))&
!                    +n_dot_ni*(2.d0*nu*drdx(il)*drdx(ik)+cte1*c_dkr(il,ik))+cte1*n(il,kip,rule)*n_i(ik)-cte2*n(ik,kip,rule)*n_i(il))
!          fs_d=d1r*(drdni*(cte1*c_dkr(il,ik)+2.d0*drdx(il)*drdx(ik))+cte1*(drdx(ik)*n_i(il)-drdx(il)*n_i(ik)))
!          ! Add to kernels
!          m(:,il,ik)=m(:,il,ik)+fs_s*phif1jw
!          l(:,il,ik)=l(:,il,ik)+fs_d*phif2jw
!        end do
!      end do
!    end do ! Loop through integrations points
!    ! Multiply m by ctes and l by cted
!    m=ctes*m
!    l=cted*l
!    ! If the normal has to be reversed, then m=-m
!    if (reverse) m=-m
!  end subroutine fbem_bem_staela2d_hbie_ext_pre

  !! This subroutine calculates the kernels for HBIE exterior integration, needing precalculated data at integration points.
  subroutine fbem_bem_staela2d_hbie_ext_pre(ps,e,reverse,x_i,n_i,mu,nu,m,l)
    implicit none
    ! I/O
    integer                :: ps                !! Selected precalculated dataset
    type(fbem_bem_element) :: e                 !! Element
    logical                :: reverse           !! Reverse normal vector
    real(kind=real64)      :: x_i(2)            !! Collocation point position vector
    real(kind=real64)      :: n_i(2)            !! Collocation point unit normal vector
    real(kind=real64)      :: mu                !! Shear modulus
    real(kind=real64)      :: nu                !! Poisson's ratio
    real(kind=real64)      :: m(e%n_pnodes,2,2) !! m integration kernels matrix
    real(kind=real64)      :: l(e%n_snodes,2,2) !! l integration kernels matrix
    ! Local
    integer           :: il, ik                 ! Counter variables for load direction and observation direction
    integer           :: kip                    ! Counter variable for integration points loop
    real(kind=real64) :: x(2)                   ! Position vector at integration point
    real(kind=real64) :: n(2)                   ! Unit normal vector at integration point
    real(kind=real64) :: pphijw(e%n_pnodes)     ! phi^p * jacobian * weight at integration point
    real(kind=real64) :: sphijw(e%n_snodes)     ! phi^s * jacobian * weight at integration point
    real(kind=real64) :: rv(2)                  ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64) :: r, d1r, d1r2           ! Distance vector module, inverse
    real(kind=real64) :: drdx(2)                ! Distance vector derivatives with respect to x_k
    real(kind=real64) :: drdn                   ! Partial derivative of r respect to unit normal
    real(kind=real64) :: drdni                  ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64) :: n_dot_ni               ! Dot product of unit normals
    real(kind=real64) :: cte1, cte2, ctes, cted ! Auxiliary constants
    real(kind=real64) :: fs_s, fs_d             ! Fundamental solutions values
    ! Initialize kernels
    m=0.d0
    l=0.d0
    ! Initialize auxiliary constants for fundamental solutions calculation
    cte1=1.0d0-2.0d0*nu
    cte2=1.0d0-4.0d0*nu
    ctes=mu/(2.0d0*c_pi*(1.0d0-nu))
    cted=-1.0d0/(4.0d0*c_pi*(1.0d0-nu))
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
      drdni=-dot_product(drdx,n_i)
      n_dot_ni=dot_product(n,n_i)
      do il=1,2
        do ik=1,2
          fs_s=d1r2*(2.d0*drdn*(drdni*(4.d0*drdx(il)*drdx(ik)-nu*c_dkr(il,ik))+nu*drdx(il)*n_i(ik)+cte1*drdx(ik)*n_i(il))&
                    -2.d0*drdni*(nu*drdx(ik)*n(il)+cte1*drdx(il)*n(ik))+n_dot_ni*(2.d0*nu*drdx(il)*drdx(ik)+cte1*c_dkr(il,ik))&
                    +cte1*n(il)*n_i(ik)-cte2*n(ik)*n_i(il))
          fs_d=d1r*(drdni*(cte1*c_dkr(il,ik)+2.d0*drdx(il)*drdx(ik))+cte1*(drdx(ik)*n_i(il)-drdx(il)*n_i(ik)))
          m(:,il,ik)=m(:,il,ik)+fs_s*pphijw
          l(:,il,ik)=l(:,il,ik)+fs_d*sphijw
        end do
      end do
    end do
    ! Multiply by constants
    m=ctes*m
    l=cted*l
    ! Reverse if needed
    if (reverse) m=-m
  end subroutine fbem_bem_staela2d_hbie_ext_pre

!  !! This subroutine calculates the kernels for HBIE exterior integration (near collocation points) using Telles transformation
!  !! within a subdivision of the element, needing only needs basic data.
!  subroutine fbem_bem_staela2d_hbie_ext_st(ngp,xi1,xi2,barxi,barr,type_g,type_f1,type_f2,delta_f,x_nodes,reverse,x_i,n_i,mu,nu,m,l)
!    implicit none
!    ! I/O
!    integer                      :: ngp                             !! Number of Gauss points (<=32)
!    real(kind=real64)            :: barxi                           !! Nearest element coordinate with respect to collocation point
!    real(kind=real64)            :: barr                            !! Telles jacobian at barxi (xi1<=barxi<=xi2)
!    real(kind=real64)            :: xi1                             !! Minimum xi coordinate of the subdivision (xi1)
!    real(kind=real64)            :: xi2                             !! Maximum xi coordinate of the subdivision (xi2)
!    integer                      :: type_g                          !! Geometrial interpolation
!    integer                      :: type_f1                         !! Functional interpolation (primary variables)
!    integer                      :: type_f2                         !! Functional interpolation (secondary variables)
!    real(kind=real64)            :: delta_f                         !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
!    real(kind=real64)            :: x_nodes(2,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
!    logical                      :: reverse                         !! Reverse normal vector
!    real(kind=real64)            :: x_i(2)                          !! Collocation point position vector
!    real(kind=real64)            :: n_i(2)                          !! Collocation point unit normal vector
!    real(kind=real64)            :: mu                              !! Shear modulus
!    real(kind=real64)            :: nu                              !! Poisson's ratio
!    real(kind=real64)            :: m(fbem_n_nodes(type_f1),2,2)    !! h kernel vector
!    real(kind=real64)            :: l(fbem_n_nodes(type_f2),2,2)    !! g kernel vector
!    ! Local
!    integer                      :: il, ik                          ! Counter variables for load direction and observation direction
!    integer                      :: kphi                            ! Counter variable for shape functions loops
!    integer                      :: nnodes_g                        ! Number of nodes of the element
!    integer                      :: kip                             ! Counter variable of integration points
!    real(kind=real64)            :: gamma                           ! Coordinate gamma (Telles transformation space [0,1])
!    real(kind=real64)            :: w                               ! Weights of an integration point
!    type(fbem_telles_parameters) :: telles_parameters               ! Telles parameters
!    real(kind=real64)            :: jt                              ! Telles jacobian: xip->gamma
!    real(kind=real64)            :: xip                             ! Coordinate xip (subdivision space [0,1])
!    real(kind=real64)            :: barxip                          ! Barxi in subdivision space
!    real(kind=real64)            :: js                              ! Subdivision jacobian: xi->xip
!    real(kind=real64)            :: xi                              ! Coordinate xi [xi1,xi2]
!    real(kind=real64)            :: aux(10)                         ! Auxiliary variable needed for shape_functions module resources
!    real(kind=real64)            :: phi_g(fbem_n_nodes(type_g))     ! Geometrical shape functions values
!    real(kind=real64)            :: dphidxi_g(fbem_n_nodes(type_g)) ! Geometrical shape functions first derivatives values
!    real(kind=real64)            :: phi_f1(fbem_n_nodes(type_f1))   ! Functional shape functions values
!    real(kind=real64)            :: phi_f2(fbem_n_nodes(type_f2))   ! Functional shape functions values
!    real(kind=real64)            :: x(2)                            ! Position vector at xi
!    real(kind=real64)            :: T(2)                            ! Tangent vector at xi
!    real(kind=real64)            :: N(2)                            ! Normal vector at xi
!    real(kind=real64)            :: rv(2)                           ! Distance vector between collocation point and integration point (x-x_i)
!    real(kind=real64)            :: r, d1r, d1r2                    ! Distance vector module, its inverse
!    real(kind=real64)            :: drdx(2)                         ! Distance vector derivatives with respect to x_k
!    real(kind=real64)            :: jg                              ! Geometric jacobian
!    real(kind=real64)            :: drdn                            ! Partial derivative of r respect to unit normal
!    real(kind=real64)            :: drdni                           ! Partial derivative of r respect to unit normal at collocation point
!    real(kind=real64)            :: n_dot_ni                        ! Dot product of unit normals
!    real(kind=real64)            :: jw                              ! Jacobians * weight
!    real(kind=real64)            :: phif1jw(fbem_n_nodes(type_f1))  ! Shape functions * jw
!    real(kind=real64)            :: phif2jw(fbem_n_nodes(type_f2))  ! Shape functions * jw
!    real(kind=real64)            :: cte1, cte2, ctes, cted          ! Auxiliary constants
!    real(kind=real64)            :: fs_s, fs_d                      ! Fundamental solutions values
!    !
!    ! Initialize
!    !
!    ! Number of nodes of the element
!    nnodes_g=fbem_n_nodes(type_g)
!    ! Initialize kernel matrices
!    m=0.d0
!    l=0.d0
!    ! Initialize auxiliary constants for fundamental solutions calculation
!    cte1=1.0d0-2.0d0*nu
!    cte2=1.0d0-4.0d0*nu
!    ctes=mu/(2.0d0*c_pi*(1.0d0-nu))
!    cted=-1.0d0/(4.0d0*c_pi*(1.0d0-nu))
!    ! Subdivision jacobian (is constant)
!    js=xi2-xi1
!    ! Convert barxi (which is in xi space [xi1,xi2]) to the subdivision space xip [0,1]
!    barxip=(barxi-xi1)/js
!    ! Calculate Telles parameters
!    telles_parameters=fbem_telles01_calculate_parameters(barxip,barr)
!    ! Loop through integrations points
!    do kip=1,gl01_n(ngp)
!      ! GAMMA->XIP->XI COORDINATE TRANSFORMATION
!      ! Gamma coordinate and weight
!      gamma=gl01_xi(kip,ngp)
!      w=gl01_w(kip,ngp)
!      ! xip coordinate, weight and jacobian from Telles transformation
!      call fbem_telles_xi_and_jacobian(telles_parameters,gamma,xip,jt)
!      ! xi coordinate
!      xi=js*xip+xi1
!      ! COMPONENTS OF THE FUNDAMENTAL SOLUTIONS
!      ! Geometrical shape functions and first derivatives at xi
!#     define etype type_g
!#     define delta 0.0d0
!#     define phi phi_g
!#     define dphidxi dphidxi_g
!#     include <phi_and_dphidxi_1d.rc>
!#     undef etype
!#     undef delta
!#     undef phi
!#     undef dphidxi
!      ! Components calculation of x and T at xi
!      x=0.d0
!      T=0.d0
!      do kphi=1,nnodes_g
!        x=x+phi_g(kphi)*x_nodes(:,kphi)
!        T=T+dphidxi_g(kphi)*x_nodes(:,kphi)
!      end do
!      ! Normal vector
!      N(1)=T(2)
!      N(2)=-T(1)
!      ! Geometric jacobian
!      jg=dot_product(T,T)
!      jg=sqrt(jg)
!      ! Unit normal vector
!      n=N/jg
!      ! Distance vector
!      rv=x-x_i
!      ! Distance vector norm and its inverse
!      r=dot_product(rv,rv)
!      r=sqrt(r)
!      d1r=1.d0/r
!      d1r2=d1r**2
!      ! r_{,k}
!      drdx=rv*d1r
!      ! dr/dn
!      drdn=dot_product(drdx,n)
!      ! dr/dn_i
!      drdni=-dot_product(drdx,n_i)
!      ! n · n_i
!      n_dot_ni=dot_product(n,n_i)
!      ! Jacobian * weight
!      jw=jg*js*jt*w
!      ! FUNCTIONAL SHAPE FUNCTIONS
!      ! Functional shape functions (primary variables) at xi
!#     define etype type_f1
!#     define delta delta_f
!#     define phi phi_f1
!#     include <phi_1d.rc>
!#     undef etype
!#     undef delta
!#     undef phi
!      ! Functional shape functions (secondary variables) at xi
!#     define etype type_f2
!#     define delta delta_f
!#     define phi phi_f2
!#     include <phi_1d.rc>
!#     undef etype
!#     undef delta
!#     undef phi
!      ! Functional shape functions * jacobians* weights
!      phif1jw=phi_f1*jw
!      phif2jw=phi_f2*jw
!      ! Loop through load and observation directions
!      do il=1,2
!        do ik=1,2
!          ! Fundamental solutions values without ctes and cted, respectively
!          fs_s=d1r2*( 2.d0*drdn*(drdni*(4.d0*drdx(il)*drdx(ik)-nu*c_dkr(il,ik))+nu*drdx(il)*n_i(ik)+cte1*drdx(ik)*n_i(il))&
!                     -2.d0*drdni*(nu*drdx(ik)*n(il)+cte1*drdx(il)*n(ik))+n_dot_ni*(2.d0*nu*drdx(il)*drdx(ik)+cte1*c_dkr(il,ik))&
!                     +cte1*n(il)*n_i(ik)-cte2*n(ik)*n_i(il))
!          fs_d=d1r*(drdni*(cte1*c_dkr(il,ik)+2.0d0*drdx(il)*drdx(ik))+cte1*(drdx(ik)*n_i(il)-drdx(il)*n_i(ik)))
!          ! Add to kernels
!          m(:,il,ik)=m(:,il,ik)+fs_s*phif1jw
!          l(:,il,ik)=l(:,il,ik)+fs_d*phif2jw
!        end do
!      end do
!    end do ! Loop through integrations points
!    ! Multiply m by ctes and l by cted
!    m=ctes*m
!    l=cted*l
!    ! If the normal has to be reversed, then m=-m
!    if (reverse) m=-m
!  end subroutine fbem_bem_staela2d_hbie_ext_st

  !! This subroutine calculates the kernels for HBIE exterior integration (near collocation points) using Telles transformation
  !! within a subdivision of the element, needing only needs basic data.
  subroutine fbem_bem_staela2d_hbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,mu,nu,gln,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)       :: e                 !! Integration element
    logical                      :: reverse           !! Reverse normal vector
    real(kind=real64)            :: xi_s(1,2)         !! Coordinates of the subdivision of the parent element (xi space [-1,1])
    real(kind=real64)            :: x_i(2)            !! Collocation point position vector
    real(kind=real64)            :: n_i(2)            !! Collocation point unit normal vector
    real(kind=real64)            :: barxip(1)         !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)            :: barr              !! Telles jacobian at barxip
    real(kind=real64)            :: mu                !! Shear modulus
    real(kind=real64)            :: nu                !! Poisson's ratio
    integer                      :: gln               !! 1D Gauss-Legendre number of integration points (<=32)
    real(kind=real64)            :: m(e%n_pnodes,2,2) !! m kernel vector
    real(kind=real64)            :: l(e%n_snodes,2,2) !! l kernel vector
    ! Local
    integer                      :: il, ik                 ! Counter variables for load direction and observation direction
    integer                      :: kphi                   ! Counter variable for shape functions loops
    integer                      :: kip                    ! Counter variable of integration points
    real(kind=real64)            :: gamma                  ! Coordinate gamma (Telles transformation space [-1,1])
    real(kind=real64)            :: w                      ! Weights of an integration point
    type(fbem_telles_parameters) :: telles_parameters      ! Telles parameters
    real(kind=real64)            :: jt                     ! Telles jacobian: xip->gamma
    real(kind=real64)            :: xip                    ! Coordinate xip (subdivision space [-1,1])
    real(kind=real64)            :: js                     ! Subdivision jacobian: xi->xip
    real(kind=real64)            :: xi                     ! Coordinate xi [xi_s(1,1),xi_s(1,2)]
    real(kind=real64)            :: aux(10)                ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)       ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi(e%n_gnodes)   ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: pphi(e%n_pnodes)       ! Functional shape functions values
    real(kind=real64)            :: sphi(e%n_snodes)       ! Functional shape functions values
    real(kind=real64)            :: x(2)                   ! Position vector at xi
    real(kind=real64)            :: T(2)                   ! Tangent vector at xi
    real(kind=real64)            :: N(2)                   ! Normal vector at xi
    real(kind=real64)            :: rv(2)                  ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)            :: r, d1r, d1r2           ! Distance vector module, its inverse
    real(kind=real64)            :: drdx(2)                ! Distance vector derivatives with respect to x_k
    real(kind=real64)            :: jg                     ! Geometric jacobian
    real(kind=real64)            :: drdn                   ! Partial derivative of r respect to unit normal
    real(kind=real64)            :: drdni                  ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)            :: n_dot_ni               ! Dot product of unit normals
    real(kind=real64)            :: jw                     ! Jacobians * weight
    real(kind=real64)            :: pphijw(e%n_pnodes)     ! Shape functions * jw
    real(kind=real64)            :: sphijw(e%n_snodes)     ! Shape functions * jw
    real(kind=real64)            :: cte1, cte2, ctes, cted ! Auxiliary constants
    real(kind=real64)            :: fs_s, fs_d             ! Fundamental solutions values
    ! Initialize kernel
    m=0.d0
    l=0.d0
    ! Initialize auxiliary constants for fundamental solutions calculation
    cte1=1.0d0-2.0d0*nu
    cte2=1.0d0-4.0d0*nu
    ctes=mu/(2.0d0*c_pi*(1.0d0-nu))
    cted=-1.0d0/(4.0d0*c_pi*(1.0d0-nu))
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
      ! Normal vector
      N(1)=T(2)
      N(2)=-T(1)
      ! Geometric jacobian
      jg=sqrt(dot_product(T,T))
      ! Unit normal vector
      n=N/jg
      ! Distance vector
      rv=x-x_i
      ! Distance vector norm and its inverse
      r=sqrt(dot_product(rv,rv))
      d1r=1.d0/r
      d1r2=d1r**2
      drdx=rv*d1r
      drdn=dot_product(drdx,n)
      drdni=-dot_product(drdx,n_i)
      n_dot_ni=dot_product(n,n_i)
      ! Jacobians * weight
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
      ! Add integration points
      do il=1,2
        do ik=1,2
          fs_s=d1r2*( 2.d0*drdn*(drdni*(4.d0*drdx(il)*drdx(ik)-nu*c_dkr(il,ik))+nu*drdx(il)*n_i(ik)+cte1*drdx(ik)*n_i(il))&
                     -2.d0*drdni*(nu*drdx(ik)*n(il)+cte1*drdx(il)*n(ik))+n_dot_ni*(2.d0*nu*drdx(il)*drdx(ik)+cte1*c_dkr(il,ik))&
                     +cte1*n(il)*n_i(ik)-cte2*n(ik)*n_i(il))
          fs_d=d1r*(drdni*(cte1*c_dkr(il,ik)+2.0d0*drdx(il)*drdx(ik))+cte1*(drdx(ik)*n_i(il)-drdx(il)*n_i(ik)))
          m(:,il,ik)=m(:,il,ik)+fs_s*pphijw
          l(:,il,ik)=l(:,il,ik)+fs_d*sphijw
        end do
      end do
    end do ! Loop through integrations points
    ! Multiply by constants
    m=ctes*m
    l=cted*l
    ! Reverse if needed
    if (reverse) m=-m
  end subroutine fbem_bem_staela2d_hbie_ext_st

  !! This subroutine calculates adaptatively the kernels for HBIE exterior integration using Telles transformation and subdivision
  !! if needed.
  recursive subroutine fbem_bem_staela2d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,mu,nu,qsp,ks,ns,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                 !! Element
    logical                  :: reverse           !! Reverse orientation
    real(kind=real64)        :: xi_s(1,2)         !! Subdivision of the parent element
    real(kind=real64)        :: x_i(2)            !! Collocation point position vector
    real(kind=real64)        :: n_i(2)            !! Collocation point unit normal vector
    real(kind=real64)        :: mu                !! Shear modulus
    real(kind=real64)        :: nu                !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp               !! Quasi-singular integration parameters
    integer                  :: ks                !! Current level of subdivisions
    integer                  :: ns                !! Maximum level of subdivision
    real(kind=real64)        :: m(e%n_pnodes,2,2) !! m integration kernels matrix
    real(kind=real64)        :: l(e%n_snodes,2,2) !! l integration kernels matrix
    ! Local
    integer           :: gln_near                 ! 1D Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer           :: gln                      ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical           :: subdivide                ! True if subdivision has to be performed
    real(kind=real64) :: barxi(1)                 ! Nearest element coordinate with respect to collocation point
    real(kind=real64) :: barxip(1)                ! Nearest element subdivision local coordinate with respect to collocation point
    real(kind=real64) :: rmin                     ! Minimum distance between collocation point and barxi on the element
    real(kind=real64) :: barr                     ! Telles jacobian at barxi
    real(kind=real64) :: cl                       ! Characteristic length
    real(kind=real64) :: d                        ! Normalized distance between collocation point and element subdivision
    integer           :: method                   ! Method used in nearest point algorithm
    real(kind=real64) :: tmp_xi_s(1,2)            ! Subdivision
    real(kind=real64) :: x_s(2,e%n_gnodes)        ! Coordinates of the element subdivision
    real(kind=real64) :: m_tmp(e%n_pnodes,2,2)    ! m integration kernels matrix (temporary)
    real(kind=real64) :: l_tmp(e%n_snodes,2,2)    ! l integration kernels matrix (temporary)
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
      cl=fbem_characteristic_length(2,e%gtype,x_s,1.d-12)
      call fbem_nearest_element_point_bem(2,e%gtype,x_s,cl,x_i,barxip,rmin,d,method)
    end if
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,6,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_staela2d_hbie_ext_adp',ns,'maximum number of subdivisions reached')
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
      call fbem_bem_staela2d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,m,l)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_staela2d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,m,l)
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_staela2d_hbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,mu,nu,gln,m_tmp,l_tmp)
      m=m+m_tmp
      l=l+l_tmp
    end if
  end subroutine fbem_bem_staela2d_hbie_ext_adp

!  !! This subroutine calculates the kernels for HBIE interior integration, needing only basic data.
!  subroutine fbem_bem_staela2d_hbie_int(ngp,type_g,type_f1,type_f2,delta_f,x_nodes,reverse,xi_i,mu,nu,m,l)
!    implicit none
!    ! I/O
!    integer                      :: ngp                             !! Number of Gauss point to be used (<=32)
!    integer                      :: type_g                          !! Geometrial interpolation
!    integer                      :: type_f1                         !! Functional interpolation (primary variables)
!    integer                      :: type_f2                         !! Functional interpolation (secondary variables)
!    real(kind=real64)            :: delta_f                         !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
!    real(kind=real64)            :: x_nodes(2,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
!    logical                      :: reverse                         !! Reverse normal vector
!    real(kind=real64)            :: xi_i                            !! Reference coordinate of the singular point.
!    real(kind=real64)            :: mu                              !! Shear modulus
!    real(kind=real64)            :: nu                              !! Poisson's ratio
!    real(kind=real64)            :: m(fbem_n_nodes(type_f1),2,2)    !! m kernel vector
!    real(kind=real64)            :: l(fbem_n_nodes(type_f2),2,2)    !! l kernel vector
!    ! Local
!    integer                      :: nnodes_g                            ! Number of nodes of the element interpolation
!    integer                      :: kphi                                ! Counter variable for shape functions loops
!    integer                      :: kip                                 ! Counter of integration points
!    real(kind=real64)            :: x_i(2)                              ! Real coordinates of collocation point
!    real(kind=real64)            :: t_i(2)                              ! Unit tangent at collocation point
!    real(kind=real64)            :: n_i(2)                              ! Unit normal at collocation point
!    real(kind=real64)            :: phi_f1(fbem_n_nodes(type_f1))       ! Functional shape functions values at xi
!    real(kind=real64)            :: phi_f2(fbem_n_nodes(type_f2))       ! Functional shape functions values at xi
!    real(kind=real64)            :: phi_f1_i(fbem_n_nodes(type_f1))     ! Functional shape functions values at xi_i
!    real(kind=real64)            :: dphidxi_f1_i(fbem_n_nodes(type_f1)) ! Functional shape functions derivatives values at xi_i
!    real(kind=real64)            :: phi_f2_i(fbem_n_nodes(type_f2))     ! Functional shape functions values at xi_i
!    real(kind=real64)            :: phi_g(fbem_n_nodes(type_g))         ! Geometrical shape functions values at xi
!    real(kind=real64)            :: dphidxi_g(fbem_n_nodes(type_g))     ! Geometrical shape functions derivatives values at xi
!    integer                      :: nsub                                ! Number of subdivision of the element
!    integer                      :: ksub                                ! Counter of subdivision
!    real(kind=real64)            :: w                                   ! Weights of each integration point
!    real(kind=real64)            :: xip                                 ! Coordinate  xip of subdivided element [0,1]
!    real(kind=real64)            :: js                                  ! Jacobian of the xi [xisub(1,:),xisub[2,:]] -> xip [0,1] transformation
!    real(kind=real64)            :: xip_i(2)                            ! Singular point in xip space
!    real(kind=real64)            :: xi                                  ! Coordinate xi
!    real(kind=real64)            :: xisub(2,2)                          ! Coordinates of subdivisions in xi space
!    real(kind=real64)            :: aux(10)                             ! Auxiliary variable needed for shape_functions module resources
!    real(kind=real64)            :: x(2)                                ! Position vector at xi
!    real(kind=real64)            :: T(2)                                ! Tangent vector at xi
!    real(kind=real64)            :: N(2)                                ! Normal vector at xi
!    real(kind=real64)            :: rv(2)                               ! Distance vector between collocation point and integration point (x-x_i)
!    real(kind=real64)            :: r, d1r, d1r2                        ! Distance vector module, its inverse
!    real(kind=real64)            :: ra, rb                              ! Distance vector from collocation point to element vertices
!    real(kind=real64)            :: drdx(2)                             ! Distance vector derivatives with respect to x_k
!    real(kind=real64)            :: jg                                  ! Geometric jacobian
!    real(kind=real64)            :: jgi                                 ! Geometric jacobian at the collocation point
!    real(kind=real64)            :: drdn                                ! Partial derivative of r respect to unit normal
!    real(kind=real64)            :: drdni                               ! Partial derivative of r respect to unit normal at collocation point
!    real(kind=real64)            :: n_dot_ni                            ! Dot product of unit normals
!    real(kind=real64)            :: drdt                                ! Partial derivative of r respect to unit tangent
!    real(kind=real64)            :: jw                                  ! Jacobians * weight
!    real(kind=real64)            :: phif1jw(fbem_n_nodes(type_f1))      ! Shape functions * jw
!    real(kind=real64)            :: phif2jw(fbem_n_nodes(type_f2))      ! Shape functions * jw
!    real(kind=real64)            :: cte1, cte2, ctes, cted              ! Auxiliary constants
!    real(kind=real64)            :: fs_d1(fbem_n_nodes(type_f2))        ! Fundamental solutions auxiliary variables
!    real(kind=real64)            :: fs_d2(fbem_n_nodes(type_f2))        ! Fundamental solutions auxiliary variables
!    real(kind=real64)            :: fs_s(fbem_n_nodes(type_f1))         ! Fundamental solutions auxiliary variables
!    !
!    ! Initialization
!    !
!    ! Number of nodes of the element
!    nnodes_g=fbem_n_nodes(type_g)
!    ! Initialize kernel matrices
!    m=0.d0
!    l=0.d0
!    ! Initialize auxiliary constants for fundamental solutions calculation
!    cte1=1.0d0-2.0d0*nu
!    cte2=1.0d0-4.0d0*nu
!    ctes=mu/(2.0d0*c_pi*(1.0d0-nu))
!    cted=-1.0d0/(4.0d0*c_pi*(1.0d0-nu))
!    ! Calculate real coordinates of collocation point
!    ! Geometrical shape functions at xi_i
!#   define etype type_g
!#   define delta 0.0d0
!#   define xi xi_i
!#   define phi phi_g
!#   define dphidxi dphidxi_g
!#   include <phi_and_dphidxi_1d.rc>
!#   undef etype
!#   undef delta
!#   undef xi
!#   undef phi
!#   undef dphidxi
!    x_i=0.d0
!    T=0.d0
!    do kphi=1,nnodes_g
!      x_i=x_i+phi_g(kphi)*x_nodes(:,kphi)
!      T=T+dphidxi_g(kphi)*x_nodes(:,kphi)
!    end do
!    ! Jacobian
!    jgi=dot_product(T,T)
!    jgi=sqrt(jgi)
!    ! Calculate unit tangent
!    t_i=T/jgi
!    ! Unit normal
!    n_i(1)=t_i(2)
!    n_i(2)=-t_i(1)
!    ! Functional shape functions and its first derivatives at xi_i
!#   define etype type_f1
!#   define delta delta_f
!#   define xi xi_i
!#   define phi phi_f1_i
!#   define dphidxi dphidxi_f1_i
!#   include <phi_and_dphidxi_1d.rc>
!#   undef etype
!#   undef delta
!#   undef xi
!#   undef phi
!#   undef dphidxi
!    ! Functional shape functions (secondary variable) at xi_i
!#   define etype type_f2
!#   define delta delta_f
!#   define xi xi_i
!#   define phi phi_f2_i
!#   include <phi_1d.rc>
!#   undef etype
!#   undef delta
!#   undef xi
!#   undef phi
!    ! If xi_i belongs to one of the vertices, subdivision is not needed
!    if (fbem_check_xi_vertex(xi_i).eqv.(.true.)) then
!      call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
!                              'the HBIE cannot be collocated at vertex')
!    ! If xi_i is inside the element, 2 subdivisions are needed
!    else
!      nsub=2
!      ! Coordinates xi of the subdivision 1
!      xisub(1,1)=-1.0d0
!      xisub(2,1)=xi_i
!      ! Coordinate xip of the collocation point
!      xip_i(1)=1.0d0
!      ! Coordinates xi of the subdivision 2
!      xisub(1,2)=xi_i
!      xisub(2,2)=1.0d0
!      ! Coordinate xip of the collocation point
!      xip_i(2)=0.0d0
!    end if
!    !
!    ! Numerical integration
!    !
!    ! Loop through subdivisions
!    do ksub=1,nsub
!      ! Jacobian of xip->xi transformation (is constant)
!      js=xisub(2,ksub)-xisub(1,ksub)
!      ! Loop through integration points
!      do kip=1,gl01_n(ngp)
!        ! XIP->XI COORDINATE TRANSFORMATION
!        ! Coordinate and weight in xip [0,1]
!        xip=gl01_xi(kip,ngp)
!        w=gl01_w(kip,ngp)
!        ! Calculate xi and xi_sub->xi jacobian
!        xi=js*xip+xisub(1,ksub)
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
!        ! Normal vector
!        N(1)=T(2)
!        N(2)=-T(1)
!        ! Geometric jacobian
!        jg=dot_product(T,T)
!        jg=sqrt(jg)
!        ! Unit normal vector
!        n=N/jg
!        ! Unit tangent vector
!        t=T/jg
!        ! Distance vector
!        rv=x-x_i
!        ! Distance vector norm and its inverse
!        r=dot_product(rv,rv)
!        r=dsqrt(r)
!        d1r=1.d0/r
!        d1r2=d1r**2
!        ! r_{,k}
!        drdx=rv*d1r
!        ! dr/dn
!        drdn=dot_product(drdx,n)
!        ! dr/dn_i
!        drdni=-dot_product(drdx,n_i)
!        ! n · n_i
!        n_dot_ni=dot_product(n,n_i)
!        ! dr/dGamma
!        drdt=dot_product(drdx,t)
!        ! Jacobian * weight
!        jw=jg*js*w
!        ! FUNCTIONAL SHAPE FUNCTIONS
!        ! Functional shape functions (primary variables) at xi
!#       define etype type_f1
!#       define delta delta_f
!#       define phi phi_f1
!#       include <phi_1d.rc>
!#       undef etype
!#       undef delta
!#       undef phi
!        ! Functional shape functions (secondary variables) at xi
!#       define etype type_f2
!#       define delta delta_f
!#       define phi phi_f2
!#       include <phi_1d.rc>
!#       undef etype
!#       undef delta
!#       undef phi
!        ! Functional shape functions * jacobians* weights
!        phif1jw=phi_f1*jw
!        phif2jw=phi_f2*jw
!        ! For l integral kernel
!        ! Common part of L11 and L22
!        fs_d1=d1r*drdni*phif2jw
!        ! L11
!        l(:,1,1)=l(:,1,1)+fs_d1*(cte1+2.d0*drdx(1)*drdx(1))
!        ! L22
!        l(:,2,2)=l(:,2,2)+fs_d1*(cte1+2.d0*drdx(2)*drdx(2))
!        ! Common terms of L12 and L21
!        fs_d1=2.d0*d1r*drdni*drdx(1)*drdx(2)*phif2jw
!        fs_d2=cte1*d1r*((n_i(1)*drdx(2)-n_i(2)*drdx(1)-drdt)*phi_f2+drdt*(phi_f2-phi_f2_i))*jw
!        ! L12
!        ! A part
!        l(:,1,2)=l(:,1,2)+fs_d1
!        ! B numerical part
!        l(:,1,2)=l(:,1,2)+fs_d2
!        ! L21
!        ! A part
!        l(:,2,1)=l(:,2,1)+fs_d1
!        ! B numerical part
!        l(:,2,1)=l(:,2,1)-fs_d2
!        ! For m integral kernel
!        ! M11
!        m(:,1,1)=m(:,1,1)+d1r2*(2.d0*(4.d0*(drdx(1)**2)*drdn*drdni+drdx(1)*(n_i(1)*drdn-n(1)*drdni))+n_dot_ni-dabs(drdt))*phif1jw
!        ! M22
!        m(:,2,2)=m(:,2,2)+d1r2*(2.d0*(4.d0*(drdx(2)**2)*drdn*drdni+drdx(2)*(n_i(2)*drdn-n(2)*drdni))+n_dot_ni-dabs(drdt))*phif1jw
!        ! Common part of M11 and M22
!        if (ksub.eq.1) then
!          fs_s=d1r2*dabs(drdt)*(phi_f1-phi_f1_i+dphidxi_f1_i/jgi*r)*jw
!        else
!          fs_s=d1r2*dabs(drdt)*(phi_f1-phi_f1_i-dphidxi_f1_i/jgi*r)*jw
!        end if
!        ! M11
!        m(:,1,1)=m(:,1,1)+fs_s
!        ! M22
!        m(:,2,2)=m(:,2,2)+fs_s
!        ! Common terms of M12 and M21
!        fs_s=d1r2*(8.d0*drdx(1)*drdx(2)*drdn*drdni+2.0d0*drdx(1)*drdx(2)*n_dot_ni+n(1)*n_i(2)+n(2)*n_i(1))
!        m(:,1,2)=m(:,1,2)+fs_s*phif1jw
!        m(:,2,1)=m(:,2,1)+fs_s*phif1jw
!      end do ! Loop through integration points
!    end do ! Loop through subdivisions
!    ! Add analytical terms
!    ! Calculate ra and rb
!    ra=sqrt((x_nodes(1,1)-x_i(1))**2+(x_nodes(2,1)-x_i(2))**2)
!    rb=sqrt((x_nodes(1,2)-x_i(1))**2+(x_nodes(2,2)-x_i(2))**2)
!    ! To L12 and L21
!    fs_d1=cte1*phi_f2_i*(log(rb)-log(ra))
!    l(:,1,2)=l(:,1,2)+fs_d1
!    l(:,2,1)=l(:,2,1)-fs_d1
!    ! To M11 and M22
!    fs_s=-phi_f1_i*(1.d0/ra+1.d0/rb)+dphidxi_f1_i/jgi*(log(rb)-log(ra))
!    m(:,1,1)=m(:,1,1)+fs_s
!    m(:,2,2)=m(:,2,2)+fs_s
!    ! Multiply m by ctes and l by cted
!    m=ctes*m
!    l=cted*l
!    ! If the normal has to be reversed, then l=-l
!    if (reverse) l=-l
!  end subroutine fbem_bem_staela2d_hbie_int

  !! This subroutine calculates the kernels for HBIE interior integration, needing only basic data.
  subroutine fbem_bem_staela2d_hbie_int(e,reverse,xi_i,mu,nu,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element) :: e                       !! Integration element
    logical                :: reverse                 !! Reverse normal vector
    real(kind=real64)      :: xi_i(1)                 !! Local coordinate of the singular point.
    real(kind=real64)      :: mu                      !! Shear modulus
    real(kind=real64)      :: nu                      !! Poisson's ratio
    real(kind=real64)      :: m(e%n_pnodes,2,2)       !! m kernels
    real(kind=real64)      :: l(e%n_snodes,2,2)       !! l kernels
    ! Local
    integer                :: gln                     ! 1D Gauss-Legendre number of integration points (<=32)
    integer                :: il, ik                  ! Counters for load and observation loops
    integer                :: kphi                    ! Counter variable for shape functions loops
    integer                :: kip                     ! Counter of integration points
    real(kind=real64)      :: x_i(2)                  ! Real coordinates of collocation point
    real(kind=real64)      :: t_i(2)                  ! Unit tangent at collocation point
    real(kind=real64)      :: n_i(2)                  ! Unit normal at collocation point
    real(kind=real64)      :: gphi(e%n_gnodes)        ! Geometrical shape functions values at xi
    real(kind=real64)      :: dgphidxi(e%n_gnodes)    ! Geometrical shape functions derivatives values at xi
    real(kind=real64)      :: pphi(e%n_pnodes)        ! Functional shape functions values at xi
    real(kind=real64)      :: sphi(e%n_snodes)        ! Functional shape functions values at xi
    real(kind=real64)      :: pphi_i(e%n_pnodes)      ! Functional shape functions values at xi_i
    real(kind=real64)      :: dpphidxi_i(e%n_pnodes)  ! Functional shape functions derivatives values at xi_i
    real(kind=real64)      :: sphi_i(e%n_snodes)      ! Functional shape functions values at xi_i
    integer                :: nsub                    ! Number of subdivision of the element
    integer                :: ksub                    ! Counter of subdivision
    real(kind=real64)      :: w                       ! Weights of each integration point
    real(kind=real64)      :: xip                     ! Coordinate  xip of subdivided element [0,1]
    real(kind=real64)      :: js                      ! Jacobian of the xi [xisub(1,:),xisub[2,:]] -> xip [0,1] transformation
    real(kind=real64)      :: xip_i(2)                ! Singular point in xip space
    real(kind=real64)      :: xi                      ! Coordinate xi
    real(kind=real64)      :: xisub(2,2)              ! Coordinates of subdivisions in xi space
    real(kind=real64)      :: aux(10)                 ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)      :: x(2)                    ! Position vector at xi
    real(kind=real64)      :: T(2)                    ! Tangent vector at xi
    real(kind=real64)      :: N(2)                    ! Normal vector at xi
    real(kind=real64)      :: rv(2)                   ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)      :: r, d1r, d1r2            ! Distance vector module, its inverse
    real(kind=real64)      :: ra, rb                  ! Distance vector from collocation point to element vertices
    real(kind=real64)      :: drdx(2)                 ! Distance vector derivatives with respect to x_k
    real(kind=real64)      :: jg                      ! Geometric jacobian
    real(kind=real64)      :: jgi                     ! Geometric jacobian at the collocation point
    real(kind=real64)      :: drdn                    ! Partial derivative of r respect to unit normal
    real(kind=real64)      :: drdni                   ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)      :: n_dot_ni                ! Dot product of unit normals
    real(kind=real64)      :: drdt                    ! Partial derivative of r respect to unit tangent
    real(kind=real64)      :: jw                      ! Jacobians * weight
    real(kind=real64)      :: pphijw(e%n_pnodes)      ! phi^p * jw
    real(kind=real64)      :: sphijw(e%n_snodes)      ! phi^s * jw
    real(kind=real64)      :: cte1, cte2, ctes, cted  ! Auxiliary constants
    real(kind=real64)      :: fs_d                    ! Fundamental solutions auxiliary variables
    real(kind=real64)      :: fs_s                    ! Fundamental solutions auxiliary variables
    ! Integration points to be used
    gln=30
    ! Initialize kernels
    m=0.d0
    l=0.d0
    ! Initialize auxiliary constants for fundamental solutions calculation
    cte1=1.0d0-2.0d0*nu
    cte2=1.0d0-4.0d0*nu
    ctes=mu/(2.0d0*c_pi*(1.0d0-nu))
    cted=-1.0d0/(4.0d0*c_pi*(1.0d0-nu))
    ! Calculate spatial coordinates at xi_i
#   define etype e%gtype
#   define delta 0.d0
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
    T=0.d0
    do kphi=1,e%n_gnodes
      x_i=x_i+gphi(kphi)*e%x(:,kphi)
      T=T+dgphidxi(kphi)*e%x(:,kphi)
    end do
    ! Jacobian
    jgi=sqrt(dot_product(T,T))
    ! Calculate unit tangent
    t_i=T/jgi
    ! Unit normal
    n_i(1)=t_i(2)
    n_i(2)=-t_i(1)
    ! Functional shape functions and its first derivatives at xi_i
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
    ! Setup the subdivisions for transformation xip -> xi: [0,1] -> [xisub_1,xisub_2]
    if (fbem_check_xi_vertex(xi_i(1))) then
      call fbem_error_message(error_unit,0,'fbem_bem_staela2d_hbie_int',0,'the HBIE cannot be collocated at vertices')
    else
      nsub=2
      xisub(1,1)=-1.0d0
      xisub(2,1)=xi_i(1)
      xip_i(1)=1.0d0
      xisub(1,2)=xi_i(1)
      xisub(2,2)=1.0d0
      xip_i(2)=0.0d0
    end if
    ! Numerical integration
    do ksub=1,nsub
      do kip=1,gl01_n(gln)
        ! XIP COORDINATE
        xip=gl01_xi(kip,gln)
        w=gl01_w(kip,gln)
        ! XIP->XI TRANSFORMATION
        js=xisub(2,ksub)-xisub(1,ksub)
        xi=js*xip+xisub(1,ksub)
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
        ! Normal vector
        N(1)=T(2)
        N(2)=-T(1)
        ! Geometric jacobian
        jg=sqrt(dot_product(T,T))
        ! Unit normal vector
        n=N/jg
        ! Unit tangent vector
        t=T/jg
        ! Distance vector and derived terms
        rv=x-x_i
        r=sqrt(dot_product(rv,rv))
        d1r=1.d0/r
        d1r2=d1r**2
        drdx=rv*d1r
        drdn=dot_product(drdx,n)
        drdni=-dot_product(drdx,n_i)
        n_dot_ni=dot_product(n,n_i)
        drdt=dot_product(drdx,t)
        ! Jacobian * weight
        jw=jg*js*w
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
        ! Functional shape functions * jacobians* weights
        pphijw=pphi*jw
        sphijw=sphi*jw
        ! Add integration points (regular parts of l)
        do il=1,2
          do ik=1,2
            fs_d=d1r*(drdni*(cte1*c_dkr(il,ik)+2.d0*drdx(il)*drdx(ik))+cte1*((n_i(il)-n(il))*drdx(ik)-(n_i(ik)-n(ik))*drdx(il)))
            l(:,il,ik)=l(:,il,ik)+fs_d*sphijw
          end do
        end do
        ! Regularized part of the CPV integrals of l
        l(:,1,2)=l(:,1,2)+cte1*d1r*drdt*(sphi-sphi_i)*jw
        l(:,2,1)=l(:,2,1)-cte1*d1r*drdt*(sphi-sphi_i)*jw
        ! Regularized part of HFP integrals of m
        do il=1,2
          fs_s=d1r2*(2.d0*(4.d0*(drdx(il)**2)*drdn*drdni+drdx(il)*(n_i(il)*drdn-n(il)*drdni))+n_dot_ni-abs(drdt))
          m(:,il,il)=m(:,il,il)+fs_s*pphijw
        end do
        do il=1,2
          if (ksub.eq.1) then
            m(:,il,il)=m(:,il,il)+d1r2*abs(drdt)*(pphi-pphi_i+dpphidxi_i/jgi*r)*jw
          else
            m(:,il,il)=m(:,il,il)+d1r2*abs(drdt)*(pphi-pphi_i-dpphidxi_i/jgi*r)*jw
          end if
        end do
        fs_s=d1r2*(8.d0*drdx(1)*drdx(2)*drdn*drdni+2.0d0*drdx(1)*drdx(2)*n_dot_ni+n(1)*n_i(2)+n(2)*n_i(1))
        m(:,1,2)=m(:,1,2)+fs_s*pphijw
        m(:,2,1)=m(:,2,1)+fs_s*pphijw
      end do ! Loop through integration points
    end do ! Loop through subdivisions
    !
    ! Analytical integrals
    !
    ra=sqrt((e%x(1,1)-x_i(1))**2+(e%x(2,1)-x_i(2))**2)
    rb=sqrt((e%x(1,2)-x_i(1))**2+(e%x(2,2)-x_i(2))**2)
    ! To L12 and L21
    l(:,1,2)=l(:,1,2)+cte1*sphi_i*(log(rb)-log(ra))
    l(:,2,1)=l(:,2,1)-cte1*sphi_i*(log(rb)-log(ra))
    ! To M11 and M22
    m(:,1,1)=m(:,1,1)-pphi_i*(1.d0/ra+1.d0/rb)+dpphidxi_i/jgi*(log(rb)-log(ra))
    m(:,2,2)=m(:,2,2)-pphi_i*(1.d0/ra+1.d0/rb)+dpphidxi_i/jgi*(log(rb)-log(ra))
    ! Multiply by constants
    m=ctes*m
    l=cted*l
    ! Reverse if needed
    if (reverse) l=-l
  end subroutine fbem_bem_staela2d_hbie_int

  ! =====================
  ! BE BODY LOAD ELEMENTS
  ! =====================

  !! Efficient automatic integration of BE body load elements
  subroutine fbem_bem_staela2d_hbie_bl_auto(e,x_i,n_i,mu,nu,qsp,ns,lbl)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                   !! Integration element
    real(kind=real64)        :: x_i(2)              !! Collocation point
    real(kind=real64)        :: n_i(2)              !! Collocation point unit normal vector
    real(kind=real64)        :: mu                  !! Shear modulus
    real(kind=real64)        :: nu                  !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp                 !! Quasi-singular integration parameters
    integer                  :: ns                  !! Maximum level of subdivisions
    real(kind=real64)        :: lbl(e%n_snodes,2,2) !! lbl integration kernel
    ! Local
    real(kind=real64)        :: r(2)                               ! Distance vector
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
        call fbem_error_message(output_unit,0,'fbem_bem_staela2d_hbie_bl_auto',0,'it is not possible to collocate at a point load')
      else
        call fbem_bem_staela2d_hbie_d(e%x(:,1),x_i,n_i,nu,lbl(1,:,:))
        return
      end if
    ! LINE OR SURFACE BODY LOAD
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
        call fbem_nearest_element_point_bem(2,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
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
        stop 'not implemented yet fbem_bem_staela2d_hbie_bl_int'
        !call fbem_bem_staela2d_hbie_bl_int(e,barxi,n_i,mu,nu,lbl)
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
          call fbem_bem_staela2d_hbie_bl_ext_pre(ps,e,x_i,n_i,mu,nu,lbl)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_staela2d_hbie_bl_ext_adp(e,xi_s,x_i,n_i,mu,nu,qsp,1,ns,lbl)
        end if
    end select
  end subroutine fbem_bem_staela2d_hbie_bl_auto

  !! This subroutine calculates the kernels for HBIE exterior integration, needing precalculated data at integration points.
  subroutine fbem_bem_staela2d_hbie_bl_ext_pre(ps,e,x_i,n_i,mu,nu,lbl)
    implicit none
    ! I/O
    integer                :: ps                  !! Selected precalculated dataset
    type(fbem_bem_element) :: e                   !! Element
    real(kind=real64)      :: x_i(2)              !! Collocation point position vector
    real(kind=real64)      :: n_i(2)              !! Collocation point unit normal vector
    real(kind=real64)      :: mu                  !! Shear modulus
    real(kind=real64)      :: nu                  !! Poisson's ratio
    real(kind=real64)      :: lbl(e%n_snodes,2,2) !! lbl integration kernels matrix
    ! Local
    integer                :: il, ik              ! Counter variables for load direction and observation direction
    integer                :: kip                 ! Counter variable for integration points loop
    real(kind=real64)      :: x(2)                ! Position vector at integration point
    real(kind=real64)      :: sphijw(e%n_snodes)  ! phi^s * jacobian * weight at integration point
    real(kind=real64)      :: rv(2)               ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)      :: r, d1r              ! Distance vector module, inverse
    real(kind=real64)      :: drdx(2)             ! Distance vector derivatives with respect to x_k
    real(kind=real64)      :: drdni               ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)      :: cte1, cted          ! Auxiliary constants
    real(kind=real64)      :: fs_d                ! Fundamental solutions values
    ! Initialize kernels
    lbl=0.d0
    ! Initialize auxiliary constants for fundamental solutions calculation
    cte1=1.0d0-2.0d0*nu
    cted=-1.0d0/(4.0d0*c_pi*(1.0d0-nu))
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      x=e%ps_x(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      d1r=1.d0/r
      drdx=rv*d1r
      drdni=-dot_product(drdx,n_i)
      do il=1,2
        do ik=1,2
          fs_d=d1r*(drdni*(cte1*c_dkr(il,ik)+2.d0*drdx(il)*drdx(ik))+cte1*(drdx(ik)*n_i(il)-drdx(il)*n_i(ik)))
          lbl(:,il,ik)=lbl(:,il,ik)+fs_d*sphijw
        end do
      end do
    end do
    ! Multiply by constants
    lbl=cted*lbl
  end subroutine fbem_bem_staela2d_hbie_bl_ext_pre

  !! This subroutine calculates the body load kernels for HBIE exterior integration (near collocation points) using Telles
  !! transformation within a subdivision of the element, needing only basic data.
  subroutine fbem_bem_staela2d_hbie_bl_ext_st(e,xi_s,x_i,n_i,barxip,barr,mu,nu,gln,lbl)
    implicit none
    ! I/O
    type(fbem_bem_element)       :: e                                  !! Integration element
    real(kind=real64)            :: xi_s(e%d,fbem_n_vertices(e%gtype)) !! Coordinates of the subdivision of the parent element (xi space [-1,1])
    real(kind=real64)            :: x_i(2)                             !! Collocation point position vector
    real(kind=real64)            :: n_i(2)                             !! Collocation point unit normal vector
    real(kind=real64)            :: barxip(e%d)                        !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)            :: barr                               !! Telles jacobian at barxip
    real(kind=real64)            :: mu                                 !! Shear modulus
    real(kind=real64)            :: nu                                 !! Poisson's ratio
    integer                      :: gln                                !! 1D Gauss-Legendre number of integration points (<=32)
    real(kind=real64)            :: lbl(e%n_snodes,2,2)                !! lbl kernels
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
    real(kind=real64)            :: x(2)                               ! Position vector at xi
    real(kind=real64)            :: T(2), T1(2), T2(2)                 ! Tangent vector at xi
    real(kind=real64)            :: rv(2)                              ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)            :: r, d1r                             ! Distance vector module, its inverse and log(1/r)
    real(kind=real64)            :: drdx(2)                            ! Distance vector derivatives with respect to x_k
    real(kind=real64)            :: drdni                              ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)            :: jg                                 ! Geometric jacobian
    real(kind=real64)            :: jw                                 ! Jacobians * weight
    real(kind=real64)            :: sphijw(e%n_snodes)                 ! secondary shape functions * jw
    real(kind=real64)            :: cte1, cted                         ! Auxiliary constants
    real(kind=real64)            :: fs_d                               ! Fundamental solutions values
    ! Initialize kernels
    lbl=0.d0
    ! Initialize auxiliary constants for fundamental solutions calculation
    cte1=1.0d0-2.0d0*nu
    cted=-1.0d0/(4.0d0*c_pi*(1.0d0-nu))
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
          do il=1,2
            do ik=1,2
              fs_d=d1r*(drdni*(cte1*c_dkr(il,ik)+2.d0*drdx(il)*drdx(ik))+cte1*(drdx(ik)*n_i(il)-drdx(il)*n_i(ik)))
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
                ! Geometric jacobian
                jg=T1(1)*T2(2)-T1(2)*T2(1)
                ! Distance vector and derived terms
                rv=x-x_i
                r=sqrt(dot_product(rv,rv))
                d1r=1.d0/r
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
                do il=1,2
                  do ik=1,2
                    fs_d=d1r*(drdni*(cte1*c_dkr(il,ik)+2.d0*drdx(il)*drdx(ik))+cte1*(drdx(ik)*n_i(il)-drdx(il)*n_i(ik)))
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
                jg=T1(1)*T2(2)-T1(2)*T2(1)
                ! Distance vector and derived terms
                rv=x-x_i
                r=sqrt(dot_product(rv,rv))
                d1r=1.d0/r
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
                do il=1,2
                  do ik=1,2
                    fs_d=d1r*(drdni*(cte1*c_dkr(il,ik)+2.d0*drdx(il)*drdx(ik))+cte1*(drdx(ik)*n_i(il)-drdx(il)*n_i(ik)))
                    lbl(:,il,ik)=lbl(:,il,ik)+fs_d*sphijw
                  end do
                end do
              end do
            end do
            case default
              call fbem_error_message(error_unit,0,__FILE__,__LINE__,'n_edges not valid')
        end select
      !
      ! OTHERS
      !
      case default
        call fbem_error_message(output_unit,0,'fbem_bem_staela2d_hbie_bl_int',0,&
                                'it is only possible to collocate at line or surface loads')
    end select
    ! Multiply by constants
    lbl=cted*lbl
  end subroutine fbem_bem_staela2d_hbie_bl_ext_st

  !! This subroutine calculates adaptatively the body load kernels for HBIE exterior integration using Telles transformation and
  !! subdivision if needed.
  recursive subroutine fbem_bem_staela2d_hbie_bl_ext_adp(e,xi_s,x_i,n_i,mu,nu,qsp,ks,ns,lbl)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                                      !! Element
    real(kind=real64)        :: xi_s(e%d,fbem_n_vertices(e%gtype))     !! Subdivision of the parent element
    real(kind=real64)        :: x_i(2)                                 !! Collocation point position vector
    real(kind=real64)        :: n_i(2)                                 !! Collocation point unit normal vector
    real(kind=real64)        :: mu                                     !! Shear modulus
    real(kind=real64)        :: nu                                     !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp                                    !! Quasi-singular integration parameters
    integer                  :: ks                                     !! Current level of subdivisions
    integer                  :: ns                                     !! Maximum level of subdivision
    real(kind=real64)        :: lbl(e%n_snodes,2,2)                    !! lbl kernels
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
    real(kind=real64)        :: x_s(2,e%n_gnodes)                      ! Coordinates of the element subdivision
    real(kind=real64)        :: lbl_tmp(e%n_snodes,2,2)                ! lbl kernels (temporary)
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
      end select
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
        call fbem_warning_message(error_unit,0,'fbem_bem_staela2d_hbie_bl_ext_adp',ns,'maximum number of subdivisions reached')
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
          call fbem_bem_staela2d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,lbl)
          ! SUBLINE 2
          tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
          tmp_xi_s(1,2)=xi_s(1,2)
          call fbem_bem_staela2d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,lbl)
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
              call fbem_bem_staela2d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,lbl)
              ! SUBTRI 2
              tmp_xi_s(:,1)=xi_s(:,2)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              call fbem_bem_staela2d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,lbl)
              ! SUBTRI 3
              tmp_xi_s(:,1)=xi_s(:,3)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              call fbem_bem_staela2d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,lbl)
              ! SUBTRI 4
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
              call fbem_bem_staela2d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,lbl)
            ! QUADRILATERALS
            case (4)
              ! SUBQUAD 1
              tmp_xi_s(:,1)=xi_s(:,1)
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,3)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,4)=0.50d0*(xi_s(:,1)+xi_s(:,4))
              call fbem_bem_staela2d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,lbl)
              ! SUBQUAD 2
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
              tmp_xi_s(:,2)=xi_s(:,2)
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,4)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              call fbem_bem_staela2d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,lbl)
              ! SUBQUAD 3
              tmp_xi_s(:,1)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
              tmp_xi_s(:,3)=xi_s(:,3)
              tmp_xi_s(:,4)=0.50d0*(xi_s(:,3)+xi_s(:,4))
              call fbem_bem_staela2d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,lbl)
              ! SUBQUAD 4
              tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,4))
              tmp_xi_s(:,2)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,3)=0.50d0*(xi_s(:,3)+xi_s(:,4))
              tmp_xi_s(:,4)=xi_s(:,4)
              call fbem_bem_staela2d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,mu,nu,qsp,ks+1,ns,lbl)
          end select
      end select
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_staela2d_hbie_bl_ext_st(e,xi_s,x_i,n_i,barxip,barr,mu,nu,gln,lbl_tmp)
      lbl=lbl+lbl_tmp
    end if
  end subroutine fbem_bem_staela2d_hbie_bl_ext_adp

  !! This subroutine calculates the body load kernels for HBIE interior integration, needing only basic data.
  subroutine fbem_bem_staela2d_hbie_bl_int(e,xi_i,n_i,mu,nu,lbl)
    implicit none
    ! I/O
    type(fbem_bem_element) :: e                   !! Integration element
    real(kind=real64)      :: xi_i(e%d)           !! Local coordinate of the singular point.
    real(kind=real64)      :: n_i(2)              !! Unit normal at the collocation point
    real(kind=real64)      :: mu                  !! Shear modulus
    real(kind=real64)      :: nu                  !! Poisson's ratio
    real(kind=real64)      :: lbl(e%n_snodes,2,2) !! lbl kernel vector
    ! Local
    integer                :: il, ik                           ! Counter variables for load direction and observation direction
    integer                :: kphi                             ! Counter variable for shape functions loops
    real(kind=real64)      :: x_i(2)                           ! Spatial coordinates at the collocation point
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
    real(kind=real64)      :: x(2)                             ! Position vector at xi
    real(kind=real64)      :: T1(2), T2(2)                     ! Tangent vectors at xi
    real(kind=real64)      :: rv(2)                            ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)      :: r, d1r                           ! Distance vector module, its inverse and log(1/r)
    real(kind=real64)      :: drdx(2)                          ! Distance vector derivatives with respect to x_k
    real(kind=real64)      :: drdni                            ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)      :: jg                               ! Geometric jacobian
    real(kind=real64)      :: jw                               ! Jacobians * weight
    real(kind=real64)      :: sphijw(e%n_snodes)               ! phi^s * jw
    real(kind=real64)      :: cte1, cted                       ! Auxiliary constants
    real(kind=real64)      :: fs_d                             ! Fundamental solutions values and auxiliary parts
    ! Initialize kernels
    lbl=0.d0
    ! Initialize auxiliary constants for fundamental solutions calculation
    cte1=1.0d0-2.0d0*nu
    cted=-1.0d0/(4.0d0*c_pi*(1.0d0-nu))
    ! Calculate
    select case (e%d)
      !
      ! SURFACE LOAD
      !
      case (2)
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
              ! xi->x jacobian
              jg=T1(1)*T2(2)-T1(2)*T2(1)
              ! Distance vector and other derived terms
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              d1r=1.d0/r
              drdx=rv*d1r
              drdni=-dot_product(drdx,n_i)
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
              ! Add integration points
              do il=1,2
                do ik=1,2
                  fs_d=d1r*(drdni*(cte1*c_dkr(il,ik)+2.d0*drdx(il)*drdx(ik))+cte1*(drdx(ik)*n_i(il)-drdx(il)*n_i(ik)))
                  lbl(:,il,ik)=lbl(:,il,ik)+fs_d*sphijw
                end do
              end do
            end do
          end do
        end do
      !
      ! OTHERS
      !
      case default
        call fbem_error_message(output_unit,0,'fbem_bem_staela2d_hbie_bl_int',0,&
                                'it is only possible to integrate surface loads')
    end select
    ! Multiply by constants
    lbl=cted*lbl
  end subroutine fbem_bem_staela2d_hbie_bl_int

  ! ================================================================================================================================
  ! VARIATION SINGULAR BOUNDARY INTEGRAL EQUATION (VSBIE)
  ! ================================================================================================================================

  !! Function that calculates one of the free-terms of the VSBIE of 2D elastic problem
  subroutine fbem_bem_staela2d_vsbie_freeterm(n_elements,n,tc,tol,nu,b)
    implicit none
    ! I/O
    integer           :: n_elements !! Number of elements (1 or 2)
    real(kind=real64) :: n(2,2)     !! Normals of elements at collocation point
    real(kind=real64) :: tc(2,2)    !! Tangents of elements at collocation point towards inside the element
    real(kind=real64) :: tol        !! Geometric tolerance [1.0e-12,1.0e-3] (default 1.0e-6)
    real(kind=real64) :: nu         !! Poisson's ratio
    real(kind=real64) :: b(2,2,2,2) !! Free-term b_{l,k,j,m}
    ! Local
    real(kind=real64) :: local_n(2,2)     ! Local copy of n
    real(kind=real64) :: local_tc(2,2)    ! Local copy of tc
    real(kind=real64) :: local_tol        ! Local copy of geometric tolerance
    real(kind=real64) :: nsum(2)          ! Unit sum of normals
    real(kind=real64) :: vectornorm       ! Norm of a vector
    real(kind=real64) :: theta(2)         ! Angles of unit tangents
    real(kind=real64) :: t1, t2           ! Auxiliary variables
    integer           :: i                ! Counter

    !
    ! Check n_elements
    !
    if ((n_elements.lt.1).or.(n_elements.gt.2)) then
      call fbem_error_message(output_unit,0,__FILE__,__LINE__,&
                                'n_elements must be 1 or 2 for free-term calculation.')
    end if
    !
    ! Check geometric tolerance
    !
    if ((tol.lt.1.0d-12).or.(tol.gt.1.0d-3)) then
      local_tol=1.0d-6
    else
      local_tol=tol
    end if
    !
    ! Check that n and tc are orthogonal (using scalar product)
    !
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
        b=0.d0
      ! For 2 elements, the freeterm is calculated as a vertex
      case (2)
        !
        ! Know who is the element 1 and 2.
        !
        ! Add up normals and normalize it
        nsum(1)=n(1,1)+n(1,2)
        nsum(2)=n(2,1)+n(2,2)
        vectornorm=dsqrt(nsum(1)**2+nsum(2)**2)
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
        ! Calculation
        t1=theta(1)
        t2=theta(2)
        b(1,1,1,1)=((-4.*(-1.+nu)*cos(t1+t2)+cos(3.*t1+t2)+cos(t1+3.*t2))*sin(t1-t2))/2.
        b(1,1,1,2)=-((cos(2.*t1)-cos(2.*t2))*(1.-2.*nu+cos(2.*t1)+cos(2.*t2)))/2.
        b(1,1,2,1)=-((cos(2.*t1)-cos(2.*t2))*(3.-2.*nu+cos(2.*t1)+cos(2.*t2)))/2.
        b(1,1,2,2)=-((-4.*(-1.+nu)*cos(t1+t2)+cos(3.*t1+t2)+cos(t1+3.*t2))*sin(t1-t2))/2.
        b(1,2,1,1)=-((cos(2.*t1)-cos(2.*t2))*(1.-2.*nu+cos(2.*t1)+cos(2.*t2)))/2.
        b(1,2,1,2)=(4.*nu*sin(2.*t1)-sin(4.*t1)-4.*nu*sin(2.*t2)+sin(4.*t2))/4.
        b(1,2,2,1)=-((-4.*(-1.+nu)*cos(t1+t2)+cos(3.*t1+t2)+cos(t1+3.*t2))*sin(t1-t2))/2.
        b(1,2,2,2)=((cos(2.*t1)-cos(2.*t2))*(1.-2.*nu+cos(2.*t1)+cos(2.*t2)))/2.
        b(2,1,1,1)=-((cos(2.*t1)-cos(2.*t2))*(-1.+2.*nu+cos(2.*t1)+cos(2.*t2)))/2.
        b(2,1,1,2)=-((4.*(-1.+nu)*cos(t1+t2)+cos(3.*t1+t2)+cos(t1+3.*t2))*sin(t1-t2))/2.
        b(2,1,2,1)=(-4.*nu*sin(2.*t1)-sin(4.*t1)+4.*nu*sin(2.*t2)+sin(4.*t2))/4.
        b(2,1,2,2)=((cos(2.*t1)-cos(2.*t2))*(-1.+2.*nu+cos(2.*t1)+cos(2.*t2)))/2.
        b(2,2,1,1)=-((4.*(-1.+nu)*cos(t1+t2)+cos(3.*t1+t2)+cos(t1+3.*t2))*sin(t1-t2))/2.
        b(2,2,1,2)=((cos(2.*t1)-cos(2.*t2))*(-3.+2.*nu+cos(2.*t1)+cos(2.*t2)))/2.
        b(2,2,2,1)=((cos(2.*t1)-cos(2.*t2))*(-1.+2.*nu+cos(2.*t1)+cos(2.*t2)))/2.
        b(2,2,2,2)=((4.*(-1.+nu)*cos(t1+t2)+cos(3.*t1+t2)+cos(t1+3.*t2))*sin(t1-t2))/2.
        b=-1.d0/(4.d0*c_pi*(1.d0-nu))*b
    end select
  end subroutine fbem_bem_staela2d_vsbie_freeterm

  !! Exterior integration using precalculated data at integration points.
  subroutine fbem_bem_staela2d_vsbie_ext_pre(ps,e,reverse,x_i,mu,nu,h1,h2,g1,g2)
    implicit none
    ! I/O
    integer                :: ps                              !! Selected precalculated dataset
    type(fbem_bem_element) :: e                               !! Element
    logical                :: reverse                         !! Reverse normal vector
    real(kind=real64)      :: x_i(2)                          !! Collocation point position vector
    real(kind=real64)      :: mu                              !! Shear modulus
    real(kind=real64)      :: nu                              !! Poisson's ratio
    real(kind=real64)      :: h1(e%dme_n_gnodes,e%n_pnodes,2,2,2) !! h1 matrix
    real(kind=real64)      :: h2(               e%n_pnodes,2,2,2) !! h2 matrix
    real(kind=real64)      :: g1(e%dme_n_gnodes,e%n_snodes,2,2,2) !! g1 matrix
    real(kind=real64)      :: g2(               e%n_snodes,2,2,2) !! g2 matrix
    ! Local
    integer           :: il, ik, im, iq                   ! Counter variables
    integer           :: kip                              ! Counter variable for integration points loop.
    real(kind=real64) :: x(2)                             ! Position vector at integration point
    real(kind=real64) :: t(2)                             ! Unit tangent vector at integration point
    real(kind=real64) :: n(2)                             ! Unit normal vector at integration point
    real(kind=real64) :: pphijw(e%n_pnodes)               ! phi^p*jacobian*weight at integration point
    real(kind=real64) :: sphijw(e%n_snodes)               ! phi^s*jacobian*weight at integration point
    real(kind=real64) :: rv(2)                            ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64) :: r, dr1, dr2, logr                ! Distance vector module, inverse and log(1/r)
    real(kind=real64) :: drdx(2)                          ! Distance vector derivatives with respect to x_k
    real(kind=real64) :: drdn                             ! Partial derivative of r respect to unit normal
    real(kind=real64) :: drdt                             ! Partial derivative of r respect to unit tangent
    real(kind=real64) :: dme_gphi(e%dme_n_gnodes)         ! phi^{dme-g} at integration point
    real(kind=real64) :: dme_wqj(e%dme_n_gnodes,2)        ! dphi{dme-g}/dx_j at integration point
    real(kind=real64) :: dme_vq(e%dme_n_gnodes)           ! vq=w_{j,q}·t_j at integration point
    real(kind=real64) :: fs_u, fs_t                       ! Fundamental solutions values
    real(kind=real64) :: fs_dudm, fs_tt, fs_dsigdm        ! Fundamental solutions values
    ! Initialize matrices
    h1=0.d0
    h2=0.d0
    g1=0.d0
    g2=0.d0
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
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
      ! Add integration points
      do il=1,2
        do ik=1,2
          ! u*_{lk}
          fs_u=-(3.d0-4.d0*nu)*c_dkr(il,ik)*logr+drdx(il)*drdx(ik)
          ! t*_{lk}
          fs_t=dr1*(2.d0*drdx(il)*drdx(ik)*drdn+(1.d0-2.d0*nu)*(drdn*c_dkr(il,ik)+drdx(ik)*n(il)-drdx(il)*n(ik)))
          ! sigma*_{lkj}·t_j
          fs_tt=dr1*(2.*drdx(il)*drdx(ik)*drdt+(1.d0-2.d0*nu)*(drdt*c_dkr(il,ik)+drdx(ik)*t(il)-drdx(il)*t(ik)))
          do im=1,2
            ! u*_{lk,m}
            fs_dudm=dr1*(-(3.d0-4.d0*nu)*drdx(im)*c_dkr(il,ik)-2.d0*drdx(il)*drdx(ik)*drdx(im)&
                        +c_dkr(il,im)*drdx(ik)+c_dkr(ik,im)*drdx(il))
            ! sigma*_{lkj,m}·n_{j}
            fs_dsigdm=dr2*(-8.d0*drdx(il)*drdx(ik)*drdx(im)*drdn&
                          -2.d0*(1.d0-2.d0*nu)*(c_dkr(il,ik)*drdn*drdx(im)+n(il)*drdx(ik)*drdx(im)-n(ik)*drdx(il)*drdx(im))&
                          +2.d0*(c_dkr(il,im)*drdx(ik)*drdn+c_dkr(ik,im)*drdx(il)*drdn+n(im)*drdx(il)*drdx(ik))&
                          +(1.d0-2.d0*nu)*(c_dkr(il,ik)*n(im)+c_dkr(ik,im)*n(il)-c_dkr(il,im)*n(ik)))
            do iq=1,e%dme_n_gnodes
              h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)+fs_t*dme_vq(iq)*t(im)*pphijw(:)
              h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)+fs_dsigdm*dme_gphi(iq)*pphijw(:)
              h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)-fs_tt*dme_vq(iq)*n(im)*pphijw(:)
              g1(iq,:,il,ik,im)=g1(iq,:,il,ik,im)+fs_u*dme_vq(iq)*t(im)*sphijw(:)
              g1(iq,:,il,ik,im)=g1(iq,:,il,ik,im)+fs_dudm*dme_gphi(iq)*sphijw(:)
            end do
            h2(:,il,ik,im)=h2(:,il,ik,im)+fs_dsigdm*pphijw(:)
            g2(:,il,ik,im)=g2(:,il,ik,im)+fs_dudm*sphijw(:)
          end do
        end do
      end do
    end do
    ! Multiply by constants
    h1=(-1.d0/(4.d0*c_pi*(1.d0-nu)))*h1
    h2=(-1.d0/(4.d0*c_pi*(1.d0-nu)))*h2
    g1=(1.d0/(8.d0*c_pi*mu*(1.d0-nu)))*g1
    g2=(1.d0/(8.d0*c_pi*mu*(1.d0-nu)))*g2
    ! Assemble depending on reversion or not
    if (reverse) then
      h1=-h1
      h2=-h2
    end if
  end subroutine fbem_bem_staela2d_vsbie_ext_pre

  !! Exterior integration (near collocation points) using Telles transformatioN within a subdivision of the element, needing only
  !! needs basic data.
  subroutine fbem_bem_staela2d_vsbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,mu,nu,gln,h1,h2,g1,g2)
    implicit none
    ! I/O
    type(fbem_bem_element) :: e                               !! Integration element
    logical                :: reverse                         !! Reverse normal vector
    real(kind=real64)      :: xi_s(1,2)                       !! Coordinates of the subdivision of the parent element (xi space [-1,1])
    real(kind=real64)      :: x_i(2)                          !! Collocation point position vector
    real(kind=real64)      :: barxip(1)                       !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)      :: barr                            !! Telles jacobian at barxip
    real(kind=real64)      :: mu                              !! Shear modulus
    real(kind=real64)      :: nu                              !! Poisson's ratio
    integer                :: gln                             !! 1D Gauss-Legendre number of integration points (<=32)
    real(kind=real64)      :: h1(e%dme_n_gnodes,e%n_pnodes,2,2,2) !! h1 matrix
    real(kind=real64)      :: h2(               e%n_pnodes,2,2,2) !! h2 matrix
    real(kind=real64)      :: g1(e%dme_n_gnodes,e%n_snodes,2,2,2) !! g1 matrix
    real(kind=real64)      :: g2(               e%n_snodes,2,2,2) !! g2 matrix
    ! Local
    integer                      :: il, ik, im, iq                   ! Counter variables
    integer                      :: kphi                             ! Counter variable for shape functions loops
    integer                      :: kip                              ! Counter variable of integration points
    real(kind=real64)            :: gamma                            ! Coordinate gamma (Telles transformation space [-1,1])
    real(kind=real64)            :: w                                ! Weights of an integration point
    type(fbem_telles_parameters) :: telles_parameters                ! Telles parameters
    real(kind=real64)            :: jt                               ! Telles jacobian: xip->gamma
    real(kind=real64)            :: xip                              ! Coordinate xip (subdivision space [-1,1])
    real(kind=real64)            :: js                               ! Subdivision jacobian: xi->xip
    real(kind=real64)            :: xi                               ! Coordinate xi [xi_s(1,1),xi_s(1,2)]
    real(kind=real64)            :: aux(10)                          ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)                 ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi(e%n_gnodes)             ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: pphi(e%n_pnodes)                 ! Functional shape functions values
    real(kind=real64)            :: sphi(e%n_snodes)                 ! Functional shape functions values
    real(kind=real64)            :: x(2)                             ! Position vector at xi
    real(kind=real64)            :: T(2)                             ! Tangent vector at xi
    real(kind=real64)            :: N(2)                             ! Normal vector at xi
    real(kind=real64)            :: rv(2)                            ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)            :: r, dr1, dr2, logr                ! Distance vector module, its inverse and log(1/r)
    real(kind=real64)            :: drdx(2)                          ! Distance vector derivatives with respect to x_k
    real(kind=real64)            :: jg                               ! Geometric jacobian
    real(kind=real64)            :: drdn                             ! Partial derivative of r respect to unit normal
    real(kind=real64)            :: drdt                             ! Partial derivative of r respect to unit tangent
    real(kind=real64)            :: jw                               ! Jacobians * weight
    real(kind=real64)            :: pphijw(e%n_pnodes)               ! Functional shape functions values * jw
    real(kind=real64)            :: sphijw(e%n_snodes)               ! Functional shape functions values * jw
    real(kind=real64)            :: dme_gphi(e%dme_n_gnodes)         ! phi^{dme-g} at integration point
    real(kind=real64)            :: dme_wqj(e%dme_n_gnodes,2)        ! w_{q,j}=dphi^{dme}_q/dx_j at integration point
    real(kind=real64)            :: dme_vq(e%dme_n_gnodes)           ! v_{q}=w_{q,j}·t_j at integration point
    real(kind=real64)            :: dme_d                            ! Dimensionless distance between DME and a point
    real(kind=real64)            :: dme_xi(e%dme_d)                  ! Local coordinate of the DME
    real(kind=real64)            :: fs_u, fs_t                       ! Fundamental solutions values
    real(kind=real64)            :: fs_dudm, fs_tt, fs_dsigdm        ! Fundamental solutions values
    ! Initialize matrices
    h1=0.d0
    h2=0.d0
    g1=0.d0
    g2=0.d0
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
      ! Add integration points
      do il=1,2
        do ik=1,2
          ! u*_{lk}
          fs_u=-(3.d0-4.d0*nu)*c_dkr(il,ik)*logr+drdx(il)*drdx(ik)
          ! t*_{lk}
          fs_t=dr1*(2.d0*drdx(il)*drdx(ik)*drdn+(1.d0-2.d0*nu)*(drdn*c_dkr(il,ik)+drdx(ik)*n(il)-drdx(il)*n(ik)))
          ! sigma*_{lkj}·t_j
          fs_tt=dr1*(2.*drdx(il)*drdx(ik)*drdt+(1.d0-2.d0*nu)*(drdt*c_dkr(il,ik)+drdx(ik)*t(il)-drdx(il)*t(ik)))
          do im=1,2
            ! u*_{lk,m}
            fs_dudm=dr1*(-(3.d0-4.d0*nu)*drdx(im)*c_dkr(il,ik)-2.d0*drdx(il)*drdx(ik)*drdx(im)&
                        +c_dkr(il,im)*drdx(ik)+c_dkr(ik,im)*drdx(il))
            ! sigma*_{lkj,m}·n_{j}
            fs_dsigdm=dr2*(-8.d0*drdx(il)*drdx(ik)*drdx(im)*drdn&
                          -2.d0*(1.d0-2.d0*nu)*(c_dkr(il,ik)*drdn*drdx(im)+n(il)*drdx(ik)*drdx(im)-n(ik)*drdx(il)*drdx(im))&
                          +2.d0*(c_dkr(il,im)*drdx(ik)*drdn+c_dkr(ik,im)*drdx(il)*drdn+n(im)*drdx(il)*drdx(ik))&
                          +(1.d0-2.d0*nu)*(c_dkr(il,ik)*n(im)+c_dkr(ik,im)*n(il)-c_dkr(il,im)*n(ik)))
            do iq=1,e%dme_n_gnodes
              h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)+fs_t*dme_vq(iq)*t(im)*pphijw(:)
              h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)+fs_dsigdm*dme_gphi(iq)*pphijw(:)
              h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)-fs_tt*dme_vq(iq)*n(im)*pphijw(:)
              g1(iq,:,il,ik,im)=g1(iq,:,il,ik,im)+fs_u*dme_vq(iq)*t(im)*sphijw(:)
              g1(iq,:,il,ik,im)=g1(iq,:,il,ik,im)+fs_dudm*dme_gphi(iq)*sphijw(:)
            end do
            h2(:,il,ik,im)=h2(:,il,ik,im)+fs_dsigdm*pphijw(:)
            g2(:,il,ik,im)=g2(:,il,ik,im)+fs_dudm*sphijw(:)
          end do
        end do
      end do
    end do
    ! Multiply by constants
    h1=(-1.d0/(4.d0*c_pi*(1.d0-nu)))*h1
    h2=(-1.d0/(4.d0*c_pi*(1.d0-nu)))*h2
    g1=(1.d0/(8.d0*c_pi*mu*(1.d0-nu)))*g1
    g2=(1.d0/(8.d0*c_pi*mu*(1.d0-nu)))*g2
    ! Assemble depending on reversion or not
    if (reverse) then
      h1=-h1
      h2=-h2
    end if
  end subroutine fbem_bem_staela2d_vsbie_ext_st

  !! Exterior integration using Telles transformation and subdivision if needed.
  recursive subroutine fbem_bem_staela2d_vsbie_ext_adp(e,reverse,xi_s,x_i,mu,nu,qsp,ks,ns,h1,h2,g1,g2)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                               !! Element
    logical                  :: reverse                         !! Reverse orientation
    real(kind=real64)        :: xi_s(1,2)                       !! Subdivision of the parent element
    real(kind=real64)        :: x_i(2)                          !! Collocation point position vector
    real(kind=real64)        :: mu                              !! Shear modulus
    real(kind=real64)        :: nu                              !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp                             !! Quasi-singular integration parameters
    integer                  :: ks                              !! Current level of subdivisions
    integer                  :: ns                              !! Maximum level of subdivision
    real(kind=real64)        :: h1(e%dme_n_gnodes,e%n_pnodes,2,2,2) !! h1 matrix
    real(kind=real64)        :: h2(               e%n_pnodes,2,2,2) !! h2 matrix
    real(kind=real64)        :: g1(e%dme_n_gnodes,e%n_snodes,2,2,2) !! g1 matrix
    real(kind=real64)        :: g2(               e%n_snodes,2,2,2) !! g2 matrix
    ! Local
    integer           :: gln_near                            ! 1D Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer           :: gln                                 ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
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
    real(kind=real64) :: h1_tmp(e%dme_n_gnodes,e%n_pnodes,2,2,2) ! h1 matrix (temporary)
    real(kind=real64) :: h2_tmp(               e%n_pnodes,2,2,2) ! h2 matrix (temporary)
    real(kind=real64) :: g1_tmp(e%dme_n_gnodes,e%n_snodes,2,2,2) ! g1 matrix (temporary)
    real(kind=real64) :: g2_tmp(               e%n_snodes,2,2,2) ! g2 matrix (temporary)
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
      cl=fbem_characteristic_length(2,e%gtype,x_s,1.d-12)
      call fbem_nearest_element_point_bem(2,e%gtype,x_s,cl,x_i,barxip,rmin,d,method)
    end if
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,6,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_staela2d_vsbie_ext_adp',ns,'maximum number of subdivisions reached')
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
      call fbem_bem_staela2d_vsbie_ext_adp(e,reverse,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,h1,h2,g1,g2)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_staela2d_vsbie_ext_adp(e,reverse,tmp_xi_s,x_i,mu,nu,qsp,ks+1,ns,h1,h2,g1,g2)
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_staela2d_vsbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,mu,nu,gln,h1_tmp,h2_tmp,g1_tmp,g2_tmp)
      h1=h1+h1_tmp
      h2=h2+h2_tmp
      g1=g1+g1_tmp
      g2=g2+g2_tmp
    end if
  end subroutine fbem_bem_staela2d_vsbie_ext_adp

  !! Interior integration
  subroutine fbem_bem_staela2d_vsbie_int(e,reverse,xi_i,mu,nu,h1,g1)
    implicit none
    ! I/O
    type(fbem_bem_element) :: e                               !! Integration element
    logical                :: reverse                         !! Reverse normal vector
    real(kind=real64)      :: xi_i(1)                         !! Local coordinate of the singular point.
    real(kind=real64)      :: mu                              !! Shear modulus
    real(kind=real64)      :: nu                              !! Poisson's ratio
    real(kind=real64)      :: h1(e%dme_n_gnodes,e%n_pnodes,2,2,2) !! h1 matrix
    real(kind=real64)      :: g1(e%dme_n_gnodes,e%n_snodes,2,2,2) !! g1 matrix
    ! Local
    integer           :: gln                              ! 1D Gauss-Legendre number of integration points (<=32)
    integer           :: il, ik, im, iq                   ! Counter variables
    integer           :: kphi                             ! Counter variable for shape functions loops
    integer           :: kip                              ! Counter of integration points
    real(kind=real64) :: x_i(2)                           ! Real coordinates of collocation point
    real(kind=real64) :: gphi(e%n_gnodes)                 ! Geometrical shape functions values at xi
    real(kind=real64) :: dgphidxi(e%n_gnodes)             ! Geometrical shape functions derivatives values at xi
    real(kind=real64) :: gphi_i(e%n_gnodes)               ! Geometrical shape functions values at xi_i
    real(kind=real64) :: dgphidxi_i(e%n_gnodes)           ! Geometrical shape functions derivatives values at xi_i
    real(kind=real64) :: pphi(e%n_pnodes)                 ! Functional shape functions values at xi
    real(kind=real64) :: sphi(e%n_snodes)                 ! Functional shape functions values at xi
    real(kind=real64) :: pphi_i(e%n_pnodes)               ! Functional shape functions values at xi_i
    integer           :: nsub                             ! Number of subdivision of the element
    integer           :: ksub                             ! Counter of subdivision
    real(kind=real64) :: w                                ! Weights of each integration point
    real(kind=real64) :: xip                              ! Coordinate  xip of subdivided element [0,1]
    real(kind=real64) :: js                               ! Jacobian of the xi [xisub(1,:),xisub[2,:]] -> xip [0,1] transformation
    real(kind=real64) :: xip_i(2)                         ! Singular point in xip space
    real(kind=real64) :: xi                               ! Coordinate xi
    real(kind=real64) :: xisub(2,2)                       ! Coordinates of element subdivisions
    real(kind=real64) :: aux(10)                          ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64) :: x(2)                             ! Position vector at xi
    real(kind=real64) :: T(2)                             ! Tangent vector at xi
    real(kind=real64) :: N(2)                             ! Normal vector at xi
    real(kind=real64) :: t_i(2), n_i(2)                   ! Unit tangent and normal at xi_i
    real(kind=real64) :: rv(2)                            ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64) :: r, dr1, dr2, logr                ! Distance vector module, its inverse and log(1/r)
    real(kind=real64) :: drdx(2)                          ! Distance vector derivatives with respect to x_k
    real(kind=real64) :: drdti                            ! (dr/dGamma)^i=+-1, -1 antes de i, 1 despues de i
    real(kind=real64) :: jg                               ! Geometric jacobian
    real(kind=real64) :: jgi                              ! Geometric jacobian at xi_i
    real(kind=real64) :: drdn                             ! Partial derivative of r respect to unit normal
    real(kind=real64) :: drdt                             ! Partial derivative of r respect to unit tangent
    real(kind=real64) :: jw                               ! Jacobians * weight
    real(kind=real64) :: dme_gphi(e%dme_n_gnodes)         ! phi^{dme-g} at integration point
    real(kind=real64) :: dme_gphi_i(e%dme_n_gnodes)       ! phi^{dme-g} at collocation point
    real(kind=real64) :: dme_wqj(e%dme_n_gnodes,2)        ! w_{q,j}=dphi^{dme}_q/dx_j at integration point
    real(kind=real64) :: dme_wqj_i(e%dme_n_gnodes,2)      ! w_{q,j}=dphi^{dme}_q/dx_j at collocation point
    real(kind=real64) :: dme_vq(e%dme_n_gnodes)           ! v_{q}=w_{q,j}·t_j at integration point
    real(kind=real64) :: dme_vq_i(e%dme_n_gnodes)         ! v_{q}=w_{q,j}·t_j at collocation point
    real(kind=real64) :: dme_d                            ! Dimensionless distance between DME and a point
    real(kind=real64) :: dme_xi(e%dme_d)                  ! Local coordinate of the DME
    real(kind=real64) :: fs_u, fs_t                       ! Fundamental solutions values
    real(kind=real64) :: fs_dudm, fs_tt                   ! Fundamental solutions values
    real(kind=real64) :: Idr1dr                           ! Int 1/r dr
    real(kind=real64) :: ra, rb                           ! Distance vector from collocation point to element vertices
    real(kind=real64) :: h(e%n_pnodes,2,2)                ! h matrix
    real(kind=real64) :: g(e%n_snodes,2,2)                ! g matrix
    !
    ! Initialization
    !
    ! Integration points to be used
    gln=32
    ! Initialize matrices
    h1=0.d0
    g1=0.d0
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
        dr2=dr1**2
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
        ! Add integration points
        do il=1,2
          do ik=1,2
          ! u*_{lk}
          fs_u=-(3.d0-4.d0*nu)*c_dkr(il,ik)*logr+drdx(il)*drdx(ik)
          ! t*_{lk}
          fs_t =dr1*(2.d0*drdx(il)*drdx(ik)*drdn+(1.d0-2.d0*nu)*(drdn*c_dkr(il,ik)+drdx(ik)*n(il)-drdx(il)*n(ik)))
          ! sigma*_{lkj}·t_j
          fs_tt=dr1*(2.d0*drdx(il)*drdx(ik)*drdt+(1.d0-2.d0*nu)*(drdt*c_dkr(il,ik)+drdx(ik)*t(il)-drdx(il)*t(ik)))
          do im=1,2
            ! u*_{lk,m}
            fs_dudm=dr1*(-(3.d0-4.d0*nu)*drdx(im)*c_dkr(il,ik)-2.d0*drdx(il)*drdx(ik)*drdx(im)&
                   +c_dkr(il,im)*drdx(ik)+c_dkr(ik,im)*drdx(il))
              do iq=1,e%dme_n_gnodes
                ! H^r
                h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)+dr2*drdn*(dme_gphi(iq)-dme_gphi_i(iq))*(-8.d0*drdx(il)*drdx(ik)*drdx(im)-2.d0*(1.d0-2.d0*nu)*c_dkr(il,ik)*drdx(im)+2.d0*(c_dkr(il,im)*drdx(ik)+c_dkr(ik,im)*drdx(il)))*pphi(:)*jw
                h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)+dr2*(dme_gphi(iq)-dme_gphi_i(iq)-dot_product(dme_wqj_i(iq,:),rv)           )*(-2.d0*(1.d0-2.d0*nu)*drdx(im)*(n(il)*drdx(ik)-n(ik)*drdx(il))+2.d0*n(im)*drdx(il)*drdx(ik)+(1.d0-2.d0*nu)*(c_dkr(il,ik)*n(im)+c_dkr(ik,im)*n(il)-c_dkr(il,im)*n(ik)))*pphi(:)*jw
                h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)+dr1*(dot_product(dme_wqj_i(iq,:),drdx)*pphi(:)-drdti*dme_vq_i(iq)*pphi_i(:))*(-2.d0*(1.d0-2.d0*nu)*drdx(im)*(n(il)*drdx(ik)-n(ik)*drdx(il))+2.d0*n(im)*drdx(il)*drdx(ik)+(1.d0-2.d0*nu)*(c_dkr(il,ik)*n(im)+c_dkr(ik,im)*n(il)-c_dkr(il,im)*n(ik)))*jw
                h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)+dme_vq_i(iq)*pphi_i(:)*dr1*drdti*(-2.d0*(1.d0-2.d0*nu)*(drdx(im)-drdti*t_i(im))*(n(il)*drdx(ik)-n(ik)*drdx(il))+2.d0*(n(im)*drdx(il)*drdx(ik)-n_i(im)*t_i(il)*t_i(ik))+(1.d0-2.d0*nu)*(c_dkr(il,ik)*(n(im)-n_i(im))+c_dkr(ik,im)*(n(il)-n_i(il))-c_dkr(il,im)*(n(ik)-n_i(ik))))*jw
                h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)+dme_vq_i(iq)*pphi_i(:)*(2.d0*n_i(im)*t_i(il)*t_i(ik)+(1.d0-2.d0*nu)*(c_dkr(il,ik)*n_i(im)+c_dkr(ik,im)*n_i(il)-c_dkr(il,im)*n_i(ik)))*dr1*(drdti-drdt)*jw
                ! H^n
                h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)-fs_tt*(dme_vq(iq)*n(im)*pphi(:)-dme_vq_i(iq)*n_i(im)*pphi_i(:))*jw
                h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)-dme_vq_i(iq)*n_i(im)*pphi_i(:)*dr1*(2.d0*drdt*(drdx(il)*drdx(ik)-t_i(il)*t_i(ik))+(1.d0-2.d0*nu)*(t(il)*drdx(ik)-t(ik)*drdx(il)))*jw
                ! H^s
                h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)+fs_t*(dme_vq(iq)*t(im)-dme_vq_i(iq)*t_i(im))*pphi(:)*jw
                ! G^r
                g1(iq,:,il,ik,im)=g1(iq,:,il,ik,im)+fs_dudm*(dme_gphi(iq)-dme_gphi_i(iq))*sphi(:)*jw
                ! G^s
                g1(iq,:,il,ik,im)=g1(iq,:,il,ik,im)+fs_u*(dme_vq(iq)*t(im)-dme_vq_i(iq)*t_i(im))*sphi(:)*jw
              end do
            end do
          end do
        end do
      end do
    end do
    !
    ! Analytical integrals
    !
    ra=sqrt((e%x(1,1)-x_i(1))**2+(e%x(2,1)-x_i(2))**2)
    rb=sqrt((e%x(1,2)-x_i(1))**2+(e%x(2,2)-x_i(2))**2)
    ! The integral of 1/r·dr depends on the position of xi_i
    select case (nsub)
      case (1)
        if (xi_i(1).lt.0.d0) Idr1dr= log(rb)
        if (xi_i(1).gt.0.d0) Idr1dr=-log(ra)
      case (2)
        Idr1dr=log(rb)-log(ra)
    end select
    do il=1,2
    do ik=1,2
    do im=1,2
    do iq=1,e%dme_n_gnodes
      h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)-2.d0*(1.d0-2.d0*nu)*dme_vq_i(iq)*pphi_i(:)*t_i(im)*(1.d0-c_dkr(il,ik))*(c_dkr(il,1)-c_dkr(il,2))*Idr1dr
      h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)+dme_vq_i(iq)*pphi_i(:)*(1.d0-2.d0*nu)*(c_dkr(ik,im)*n_i(il)-c_dkr(il,im)*n_i(ik))*Idr1dr
    end do
    end do
    end do
    end do
    ! Multiply by constants
    h1=(-1.d0/(4.d0*c_pi*(1.d0-nu)))*h1
    g1=(1.d0/(8.d0*c_pi*mu*(1.d0-nu)))*g1
    ! Assemble depending on reversion or not
    if (reverse) h1=-h1
    ! Assemble matrices from SBIE
    call fbem_bem_staela2d_sbie_int(e,reverse,xi_i,mu,nu,h,g)
    do im=1,2
      do iq=1,e%dme_n_gnodes
        h1(iq,:,:,:,im)=h1(iq,:,:,:,im)+dme_vq_i(iq)*t_i(im)*h(:,:,:)
        g1(iq,:,:,:,im)=g1(iq,:,:,:,im)+dme_vq_i(iq)*t_i(im)*g(:,:,:)
      end do
    end do
  end subroutine fbem_bem_staela2d_vsbie_int

  !! Efficient automatic integration of boundary elements
  subroutine fbem_bem_staela2d_vsbie_auto(e,reverse,x_i,mu,nu,qsp,ns,h1,h2,g1,g2)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                               !! Integration element
    logical                  :: reverse                         !! Reverse orientation
    real(kind=real64)        :: x_i(2)                          !! Collocation point
    real(kind=real64)        :: mu                              !! Shear modulus
    real(kind=real64)        :: nu                              !! Poisson's ratio
    type(fbem_qs_parameters) :: qsp                             !! Quasi-singular integration parameters
    integer                  :: ns                              !! Maximum level of subdivisions
    real(kind=real64)        :: h1(e%dme_n_gnodes,e%n_pnodes,2,2,2) !! h1 matrix
    real(kind=real64)        :: h2(               e%n_pnodes,2,2,2) !! h2 matrix
    real(kind=real64)        :: g1(e%dme_n_gnodes,e%n_snodes,2,2,2) !! g1 matrix
    real(kind=real64)        :: g2(               e%n_snodes,2,2,2) !! g2 matrix
    ! Local
    real(kind=real64) :: r(2)              ! Distance vector
    real(kind=real64) :: rmin              ! Minimum distance between element and x_i
    real(kind=real64) :: barxi(1)          ! Local coordinates of the nearest element point with respect to x_i
    real(kind=real64) :: d                 ! Dimensionless distance
    integer           :: delta             ! Control variable
    real(kind=real64) :: xi_s(1,2)         ! Local coordinates of the element subdivision
    integer           :: method            ! Method used when calculating the nearest element point
    integer           :: gln_near          ! 1D Gauss-Legendre integration points required by the quasi-singular function
    integer           :: gln               ! 1D Gauss-Legendre integration points used in the integration
    integer           :: ps                ! Selected precalculated dataset
    integer           :: i                 ! Counter
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
        call fbem_bem_staela2d_vsbie_int(e,reverse,barxi,mu,nu,h1,g1)
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
          call fbem_bem_staela2d_vsbie_ext_pre(ps,e,reverse,x_i,mu,nu,h1,h2,g1,g2)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_staela2d_vsbie_ext_adp(e,reverse,xi_s,x_i,mu,nu,qsp,1,ns,h1,h2,g1,g2)
        end if
    end select
  end subroutine fbem_bem_staela2d_vsbie_auto

  ! ================================================================================================================================

  ! ================================================================================================================================
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE) - HALF-SPACE FUNDAMENTAL SOLUTION (COMPLEMENTARY PART)

  ! Half-plane complementary fundamental solution for elastostatics, Telles (1981).
  subroutine fbem_bem_staela2d_hfc_sbie(mu,nu,hsn,hsx,x_inp,n_inp,x_ip_inp,cte,reg,u_out,t_out)
    implicit none
    ! I/O variables
    real(kind=real64) :: mu
    real(kind=real64) :: nu
    integer           :: hsn
    real(kind=real64) :: hsx
    real(kind=real64) :: x_inp(2)
    real(kind=real64) :: n_inp(2)
    real(kind=real64) :: x_ip_inp(2)
    real(kind=real64) :: u_out(2,2)
    real(kind=real64) :: dudx_out(2,2,2)
    logical           :: cte ! True if common constants included
    logical           :: reg ! True if only the regular part
    real(kind=real64) :: sigma_out(2,2,2)
    real(kind=real64) :: t_out(2,2)
    ! Local variables
    real(kind=real64) :: x(2)
    real(kind=real64) :: n(2)
    real(kind=real64) :: x_i(2)
    real(kind=real64) :: u(2,2)
    real(kind=real64) :: dudx(2,2,2)
    real(kind=real64) :: sigma(2,2,2)
    real(kind=real64) :: t(2,2)
    integer           :: il, ik, ij
    real(kind=real64) :: lambda
    real(kind=real64) :: p_cte_u, p_cte_t, p_u(3)
    real(kind=real64) :: x_ip(2)
    !real(kind=real64) :: rv(2), r, drdx(2), drdn
    real(kind=real64) :: sv(2), s, dsdx(2), dsdn
    real(kind=real64) :: logs, d1s1, d1s2, theta
    !
    ! Transform the real half-space configuration (input) to the reference configuration (x1>=0)
    !
    select case (hsn)
      case (-1)
        x  (1)  = x_inp(1)-hsx
        x  (2)  = x_inp(2)
        n  (1)  = n_inp(1)
        n  (2)  = n_inp(2)
        x_ip(1) = x_ip_inp(1)-hsx
        x_ip(2) = x_ip_inp(2)
      case ( 1)
        x  (1)  =-(x_inp(1)-hsx)
        x  (2)  =-x_inp(2)
        n  (1)  =-n_inp(1)
        n  (2)  =-n_inp(2)
        x_ip(1) =-(x_ip_inp(1)-hsx)
        x_ip(2) =-x_ip_inp(2)
      case (-2)
        x  (1)  = x_inp(2)-hsx
        x  (2)  =-x_inp(1)
        n  (1)  = n_inp(2)
        n  (2)  =-n_inp(1)
        x_ip(1) = x_ip_inp(2)-hsx
        x_ip(2) =-x_ip_inp(1)
      case ( 2)
        x  (1)  =-(x_inp(2)-hsx)
        x  (2)  = x_inp(1)
        n  (1)  =-n_inp(2)
        n  (2)  = n_inp(1)
        x_ip(1) =-(x_ip_inp(2)-hsx)
        x_ip(2) = x_ip_inp(1)
      case default
        call fbem_error_message(error_unit,0,'hsn',hsn,'invalid value (hfn={-1,1,-2,2})')
    end select
    !
    ! Half-plane complementary fundamental solution for x1>=0
    !
    ! Source collocation point
    x_i(1)=-x_ip(1)
    x_i(2)= x_ip(2)
    ! Check collocation and observation point position
    if ((x(1).lt.0.d0).or.(x_i(1).lt.0.d0)) then
      call fbem_error_message(error_unit,0,'fbem_bem_staela2d_hfc_sbie',0,'col. and obs. points must be inside half-space.')
    end if
    ! Parameters
    lambda=2.*mu*nu/(1.d0-2.*nu)
    p_cte_u=1./(8.*c_pi*mu*(1.-nu))
    p_cte_t=-1./(4.*c_pi*(1.-nu))
    p_u(1)=-8.*(1.-nu)**2+3.-4.*nu
    p_u(2)=3.-4.*nu
    p_u(3)=4.*(1.-nu)*(1.-2.*nu)
    ! Components
    sv=x-x_ip
    s=sqrt(dot_product(sv,sv))
    logs=log(s)
    d1s1=1.d0/s
    d1s2=d1s1*d1s1
    dsdx=sv/s
    theta=atan2(sv(2),sv(1))
    dsdn=dot_product(dsdx,n)
    if (reg) then
      do il=1,2
        do ik=1,2
          u(il,ik)=p_u(2)*dsdx(il)*dsdx(ik)-(1.-c_dkr(il,ik))*2.*p_u(2)*x_i(1)*d1s1*dsdx(2)&
                  +2.*x(1)*x_i(1)*d1s2*(2.*dsdx(ik)*(c_dkr(il,1)*dsdx(1)-c_dkr(il,2)*dsdx(2))&
                  -c_dkr(il,1)*c_dkr(ik,1)+c_dkr(il,2)*c_dkr(ik,2))-c_lc2(il,ik)*p_u(3)*theta
          t(il,ik)=(p_u(2)*2.*dsdx(il)*dsdx(ik)-(1.-2.*nu)*c_dkr(il,ik))*d1s1*dsdn&
                  -2.*x_i(1)*d1s2*(-nu*p_u(2)/(1.-2.*nu)*n(ik)*(c_dkr(il,1)+2.*dsdx(2)*(dsdx(il)-dsdx(1)-dsdx(2)))&
                  -0.5*p_u(2)*((1.-c_dkr(il,ik))*(n(2)-2.*dsdx(2)*dsdn)+(n(1)+n(2)-n(il))*(c_dkr(2,ik)-2.*dsdx(2)*dsdx(ik)))&
                  +(c_dkr(il,1)*dsdx(1)-c_dkr(il,2)*dsdx(2))*(2.*nu/(1.-2.*nu)*n(ik)*dsdx(1)+n(1)*dsdx(ik)+c_dkr(1,ik)*dsdn)&
                  -nu/(1.-2.*nu)*c_dkr(il,1)*n(ik)-c_dkr(il,1)*c_dkr(ik,1)*n(1)+0.5*c_dkr(il,2)*(c_dkr(ik,2)*n(1)+c_dkr(ik,1)*n(2))&
                  +2.*x(1)*d1s1*((n(ik)-4.*dsdx(ik)*dsdn)*(c_dkr(il,1)*dsdx(1)-c_dkr(il,2)*dsdx(2))&
                  +c_dkr(il,1)*(c_dkr(ik,1)*dsdn+n(1)*dsdx(ik))-c_dkr(il,2)*(c_dkr(ik,2)*dsdn+n(2)*dsdx(ik))))
        end do
      end do
    else
      do il=1,2
        do ik=1,2
          u(il,ik)=c_dkr(il,ik)*p_u(1)*logs+p_u(2)*dsdx(il)*dsdx(ik)-(1.-c_dkr(il,ik))*2.*p_u(2)*x_i(1)*d1s1*dsdx(2)&
                  +2.*x(1)*x_i(1)*d1s2*(2.*dsdx(ik)*(c_dkr(il,1)*dsdx(1)-c_dkr(il,2)*dsdx(2))&
                  -c_dkr(il,1)*c_dkr(ik,1)+c_dkr(il,2)*c_dkr(ik,2))-c_lc2(il,ik)*p_u(3)*theta
          t(il,ik)=(p_u(2)*2.*dsdx(il)*dsdx(ik)-(1.-2.*nu)*c_dkr(il,ik))*d1s1*dsdn-(1.-2.*nu)*d1s1*(n(il)*dsdx(ik)-n(ik)*dsdx(il))&
                  -2.*x_i(1)*d1s2*(-nu*p_u(2)/(1.-2.*nu)*n(ik)*(c_dkr(il,1)+2.*dsdx(2)*(dsdx(il)-dsdx(1)-dsdx(2)))&
                  -0.5*p_u(2)*((1.-c_dkr(il,ik))*(n(2)-2.*dsdx(2)*dsdn)+(n(1)+n(2)-n(il))*(c_dkr(2,ik)-2.*dsdx(2)*dsdx(ik)))&
                  +(c_dkr(il,1)*dsdx(1)-c_dkr(il,2)*dsdx(2))*(2.*nu/(1.-2.*nu)*n(ik)*dsdx(1)+n(1)*dsdx(ik)+c_dkr(1,ik)*dsdn)&
                  -nu/(1.-2.*nu)*c_dkr(il,1)*n(ik)-c_dkr(il,1)*c_dkr(ik,1)*n(1)+0.5*c_dkr(il,2)*(c_dkr(ik,2)*n(1)+c_dkr(ik,1)*n(2))&
                  +2.*x(1)*d1s1*((n(ik)-4.*dsdx(ik)*dsdn)*(c_dkr(il,1)*dsdx(1)-c_dkr(il,2)*dsdx(2))&
                  +c_dkr(il,1)*(c_dkr(ik,1)*dsdn+n(1)*dsdx(ik))-c_dkr(il,2)*(c_dkr(ik,2)*dsdn+n(2)*dsdx(ik))))
        end do
      end do
    end if
    if (cte) then
      u=p_cte_u*u
      t=p_cte_t*t
    end if
    !
    ! Transform the u and t from the reference configuration (x1>=0) to the real half-space configuration
    !
    if (abs(hsn).eq.2) then
      u_out(1,1) = u(2,2)
      u_out(1,2) =-u(2,1)
      u_out(2,1) =-u(1,2)
      u_out(2,2) = u(1,1)
      t_out(1,1) = t(2,2)
      t_out(1,2) =-t(2,1)
      t_out(2,1) =-t(1,2)
      t_out(2,2) = t(1,1)
    else
      u_out=u
      t_out=t
    end if
  end subroutine fbem_bem_staela2d_hfc_sbie

  subroutine fbem_bem_staela2d_hfc_sbie_ext_pre(ps,e,reverse,x_ip,mu,nu,hsn,hsx,h,g)
    implicit none
    ! I/O
    integer                :: ps                   !! Selected precalculated dataset
    type(fbem_bem_element) :: e                    !! Element
    logical                :: reverse              !! Reverse normal vector
    real(kind=real64)      :: x_ip(2)              !! Image collocation point position vector
    real(kind=real64)      :: mu                   !! Shear modulus
    real(kind=real64)      :: nu                   !! Poisson's ratio
    integer                :: hsn                  !! Unit normal of the half-plane: -1,1,-2,2
    real(kind=real64)      :: hsx                  !! Coordinate of the half-plane
    real(kind=real64)      :: h(e%n_pnodes,2,2)    !! h integration kernels matrix
    real(kind=real64)      :: g(e%n_snodes,2,2)    !! g integration kernels matrix
    ! Local
    integer                :: l, k                 ! Counter variables for load direction and observation direction
    integer                :: kip                  ! Counter variable for integration points loop
    real(kind=real64)      :: x(2)                 ! Position vector at integration point
    real(kind=real64)      :: n(2)                 ! Unit normal vector at integration point
    real(kind=real64)      :: pphijw(e%n_pnodes)   ! phi^p * jacobian * weight at integration point
    real(kind=real64)      :: sphijw(e%n_snodes)   ! phi^s * jacobian * weight at integration point
    real(kind=real64)      :: fs_u(2,2)            ! Complementary fundamental solution u*
    real(kind=real64)      :: fs_t(2,2)            ! Complementary fundamental solution t*
    ! Initialize kernels
    h=0.d0
    g=0.d0
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      x=e%ps_x(:,kip,ps)
      n=e%ps_n(:,kip,ps)
      pphijw=e%ps_pphijw(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      call fbem_bem_staela2d_hfc_sbie(mu,nu,hsn,hsx,x,n,x_ip,.true.,.false.,fs_u,fs_t)
      do l=1,2
        do k=1,2
          h(:,l,k)=h(:,l,k)+fs_t(l,k)*pphijw
          g(:,l,k)=g(:,l,k)+fs_u(l,k)*sphijw
        end do
      end do
    end do
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_staela2d_hfc_sbie_ext_pre

  subroutine fbem_bem_staela2d_hfc_sbie_ext_st(e,reverse,xi_s,x_ip,barxip,barr,mu,nu,hsn,hsx,gln,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)       :: e                 !! Integration element
    logical                      :: reverse           !! Reverse normal vector
    real(kind=real64)            :: xi_s(1,2)         !! Coordinates of the subdivision of the parent element (xi space [-1,1])
    real(kind=real64)            :: x_ip(2)           !! Image collocation point position vector
    real(kind=real64)            :: barxip(1)         !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)            :: barr              !! Telles jacobian at barxip
    real(kind=real64)            :: mu                !! Shear modulus
    real(kind=real64)            :: nu                !! Poisson's ratio
    integer                      :: hsn               !! Unit normal of the half-plane: -1,1,-2,2
    real(kind=real64)            :: hsx               !! Coordinate of the half-plane
    integer                      :: gln               !! 1D Gauss-Legendre number of integration points (<=32)
    real(kind=real64)            :: h(e%n_pnodes,2,2) !! h kernel vector
    real(kind=real64)            :: g(e%n_snodes,2,2) !! g kernel vector
    ! Local
    integer                      :: l, k                       ! Counter variables for load direction and observation direction
    integer                      :: kphi                       ! Counter variable for shape functions loops
    integer                      :: kip                        ! Counter variable of integration points
    real(kind=real64)            :: gamma                      ! Coordinate gamma (Telles transformation space [-1,1])
    real(kind=real64)            :: w                          ! Weights of an integration point
    type(fbem_telles_parameters) :: telles_parameters          ! Telles parameters
    real(kind=real64)            :: jt                         ! Telles jacobian: xip->gamma
    real(kind=real64)            :: xip                        ! Coordinate xip (subdivision space [-1,1])
    real(kind=real64)            :: js                         ! Subdivision jacobian: xi->xip
    real(kind=real64)            :: xi                         ! Coordinate xi [xi_s(1,1),xi_s(1,2)]
    real(kind=real64)            :: aux(10)                    ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)           ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi(e%n_gnodes)       ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: pphi(e%n_pnodes)           ! Functional shape functions values
    real(kind=real64)            :: sphi(e%n_snodes)           ! Functional shape functions values
    real(kind=real64)            :: x(2)                       ! Position vector at xi
    real(kind=real64)            :: T(2)                       ! Tangent vector at xi
    real(kind=real64)            :: N(2)                       ! Normal vector at xi
    real(kind=real64)            :: jg                         ! Geometric jacobian
    real(kind=real64)            :: jw                         ! Jacobians * weight
    real(kind=real64)            :: pphijw(e%n_pnodes)         ! primary shape functions * jw
    real(kind=real64)            :: sphijw(e%n_snodes)         ! secondary shape functions * jw
    real(kind=real64)            :: fs_u(2,2), fs_t(2,2)                 ! Fundamental solutions values
    ! Initialization
    h=0.d0
    g=0.d0
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
      ! Jacobians * weight
      jw=jg*js*jt*w
      ! Complementary fundamental solutions
      call fbem_bem_staela2d_hfc_sbie(mu,nu,hsn,hsx,x,n,x_ip,.true.,.false.,fs_u,fs_t)
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
      ! Add integration points
      do l=1,2
        do k=1,2
          h(:,l,k)=h(:,l,k)+fs_t(l,k)*pphijw
          g(:,l,k)=g(:,l,k)+fs_u(l,k)*sphijw
        end do
      end do
    end do
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_staela2d_hfc_sbie_ext_st

  recursive subroutine fbem_bem_staela2d_hfc_sbie_ext_adp(e,reverse,xi_s,x_ip,mu,nu,hsn,hsx,qsp,ks,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                 !! Element
    logical                  :: reverse           !! Reverse orientation
    real(kind=real64)        :: xi_s(1,2)         !! Subdivision of the parent element
    real(kind=real64)        :: x_ip(2)           !! Image collocation point position vector
    real(kind=real64)        :: mu                !! Shear modulus
    real(kind=real64)        :: nu                !! Poisson's ratio
    integer                  :: hsn               !! Unit normal of the half-plane: -1,1,-2,2
    real(kind=real64)        :: hsx               !! Coordinate of the half-plane
    type(fbem_qs_parameters) :: qsp               !! Quasi-singular integration parameters
    integer                  :: ks                !! Current level of subdivisions
    integer                  :: ns                !! Maximum level of subdivision
    real(kind=real64)        :: h(e%n_pnodes,2,2) !! h integration kernels matrix
    real(kind=real64)        :: g(e%n_snodes,2,2) !! g integration kernels matrix
    ! Local
    integer           :: gln_near                 ! 1D Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer           :: gln                      ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical           :: subdivide                ! True if subdivision has to be performed
    real(kind=real64) :: barxi(1)                 ! Nearest element coordinate with respect to collocation point
    real(kind=real64) :: barxip(1)                ! Nearest element subdivision local coordinate with respect to collocation point
    real(kind=real64) :: rmin                     ! Minimum distance between collocation point and barxi on the element
    real(kind=real64) :: barr                     ! Telles jacobian at barxi
    real(kind=real64) :: cl                       ! Characteristic length
    real(kind=real64) :: d                        ! Normalized distance between collocation point and element subdivision
    integer           :: method                   ! Method used in nearest point algorithm
    real(kind=real64) :: tmp_xi_s(1,2)            ! Subdivision
    real(kind=real64) :: x_s(2,e%n_gnodes)        ! Coordinates of the element subdivision
    real(kind=real64) :: h_tmp(e%n_pnodes,2,2)    ! h integration kernels matrix (temporary)
    real(kind=real64) :: g_tmp(e%n_snodes,2,2)    ! g integration kernels matrix (temporary)
    ! Initialize
    if (ks.eq.1) then
      h=0.d0
      g=0.d0
      xi_s(1,1)=-1.d0
      xi_s(1,2)= 1.d0
      call fbem_nearest_element_point_bem(2,e%gtype,e%x,e%cl,x_ip,barxi,rmin,d,method)
      barxip=barxi
    else
      call fbem_obtain_element_subdivision_coordinates(2,e%gtype,e%x,xi_s,x_s)
      cl=fbem_characteristic_length(2,e%gtype,x_s,1.d-12)
      call fbem_nearest_element_point_bem(2,e%gtype,x_s,cl,x_ip,barxip,rmin,d,method)
    end if
    ! Obtain an estimation of the number of Gaussian points
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,4,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_staela2d_hfc_sbie_ext_adp',ns,'maximum number of subdivisions reached')
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
      call fbem_bem_staela2d_hfc_sbie_ext_adp(e,reverse,tmp_xi_s,x_ip,mu,nu,hsn,hsx,qsp,ks+1,ns,h,g)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_staela2d_hfc_sbie_ext_adp(e,reverse,tmp_xi_s,x_ip,mu,nu,hsn,hsx,qsp,ks+1,ns,h,g)
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_staela2d_hfc_sbie_ext_st(e,reverse,xi_s,x_ip,barxip,barr,mu,nu,hsn,hsx,gln,h_tmp,g_tmp)
      h=h+h_tmp
      g=g+g_tmp
    end if
  end subroutine fbem_bem_staela2d_hfc_sbie_ext_adp

  subroutine fbem_bem_staela2d_hfc_sbie_int(e,reverse,xi_i,mu,nu,hsn,hsx,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element) :: e                          !! Integration element
    logical                :: reverse                    !! Reverse normal vector
    real(kind=real64)      :: xi_i(1)                    !! Local coordinate of the singular point.
    real(kind=real64)      :: mu                         !! Shear modulus
    real(kind=real64)      :: nu                         !! Poisson's ratio
    integer                :: hsn                        !! Unit normal of the half-plane: -1,1,-2,2
    real(kind=real64)      :: hsx                        !! Coordinate of the half-plane
    real(kind=real64)      :: h(e%n_pnodes,2,2)          !! h kernel vector
    real(kind=real64)      :: g(e%n_snodes,2,2)          !! g kernel vector
    ! Local
    integer                :: gln                        ! Gauss-Legendre number of integration points (<=32)
    integer                :: l, k                       ! Counter variables for load direction and observation direction
    integer                :: kphi                       ! Counter variable for shape functions loops
    integer                :: kip                        ! Counter of integration points
    real(kind=real64)      :: x_i(2)                     ! Real coordinates of collocation point
    real(kind=real64)      :: gphi(e%n_gnodes)           ! Geometrical shape functions values at xi
    real(kind=real64)      :: dgphidxi(e%n_gnodes)       ! Geometrical shape functions derivatives values at xi
    real(kind=real64)      :: pphi(e%n_pnodes)           ! Functional shape functions values at xi
    real(kind=real64)      :: sphi(e%n_snodes)           ! Functional shape functions values at xi
    real(kind=real64)      :: pphi_i(e%n_pnodes)         ! Functional shape functions values at xi_i
    integer                :: nsub                       ! Number of subdivision of the element
    integer                :: ksub                       ! Counter of subdivision
    real(kind=real64)      :: w                          ! Weights of each integration point
    real(kind=real64)      :: xip                        ! Coordinate  xip of subdivided element [0,1]
    real(kind=real64)      :: js                         ! Jacobian of the xi [xisub(1,:),xisub[2,:]] -> xip [0,1] transformation
    real(kind=real64)      :: xip_i(2)                   ! Singular point in xip space
    real(kind=real64)      :: xi                         ! Coordinate xi
    real(kind=real64)      :: xisub(2,2)                 ! Coordinates of element subdivisions
    real(kind=real64)      :: aux(10)                    ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)      :: x(2)                       ! Position vector at xi
    real(kind=real64)      :: T(2)                       ! Tangent vector at xi
    real(kind=real64)      :: N(2)                       ! Normal vector at xi
    real(kind=real64)      :: rv(2)                      ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)      :: r, d1r, logd1r             ! Distance vector module, its inverse and log(1/r)
    real(kind=real64)      :: ra, rb                     ! Distance vector from collocation point to element vertices
    real(kind=real64)      :: drdx(2)                    ! Distance vector derivatives with respect to x_k
    real(kind=real64)      :: jg                         ! Geometric jacobian
    real(kind=real64)      :: drdn                       ! Partial derivative of r respect to unit normal
    real(kind=real64)      :: drdt                       ! Partial derivative of r respect to unit tangent
    real(kind=real64)      :: jw                         ! Jacobians * weight
    real(kind=real64)      :: pphijw(e%n_pnodes)         ! phi^p * jw
    real(kind=real64)      :: sphijw(e%n_snodes)         ! phi^s * jw
    real(kind=real64)      :: lsphi(e%n_snodes)          ! log(1/r)*phi^s integrals
    real(kind=real64)      :: cteu1, cteu2, cteu3        ! Auxiliary constants
    real(kind=real64)      :: ctet1, ctet2               ! Auxiliary constants
    real(kind=real64)      :: fs_u(2,2)                  ! Fundamental solutions values
    real(kind=real64)      :: fs_t(2,2)                  ! Fundamental solutions values
    ! Integration points to be used
    gln=32
    ! Initialize kernels
    h=0.d0
    g=0.d0
    ! Initialize auxiliary constants for fundamental solutions calculation
    cteu1=1.0d0/(8.0d0*c_pi*mu*(1.0d0-nu))
    cteu2=3.0d0-4.0d0*nu
    cteu3=4.*(1.-nu)*(1.-2.*nu)
    ctet1=-1.0d0/(4.0d0*c_pi*(1.0d0-nu))
    ctet2=1.0d0-2.0d0*nu
    ! Setup the subdivisions for transformation xip -> xi: [0,1] -> [xisub_1,xisub_2]
    if (fbem_check_xi_vertex(xi_i(1))) then
      nsub=1
      xisub(1,1)=-1.0d0
      xisub(2,1)= 1.0d0
      if (xi_i(1).lt.0.0d0) xip_i(1)=0.0d0
      if (xi_i(1).gt.0.0d0) xip_i(1)=1.0d0
    else
      nsub=2
      xisub(1,1)=-1.0d0
      xisub(2,1)=xi_i(1)
      xip_i(1)=1.0d0
      xisub(1,2)=xi_i(1)
      xisub(2,2)=1.0d0
      xip_i(2)=0.0d0
    end if
    ! Calculate spatial coordinates at xi_i
#   define etype e%gtype
#   define delta 0.d0
#   define xi xi_i(1)
#   define phi gphi
#   include <phi_1d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
    x_i=0.d0
    do kphi=1,e%n_gnodes
      x_i=x_i+gphi(kphi)*e%x(:,kphi)
    end do
    ! Save functional shape functions at xi_i
#   define etype e%ptype
#   define delta e%ptype_delta
#   define xi xi_i(1)
#   define phi pphi_i
#   include <phi_1d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
    !
    ! Weakly singular integrals
    !
    call fbem_bem_logd1rphi_1d_int(2,e%gtype,e%x,e%stype,e%stype_delta,xi_i,lsphi)
    g(:,1,1)=-(-8.*(1.-nu)**2+3.-4.*nu)*lsphi
    g(:,2,2)=-(-8.*(1.-nu)**2+3.-4.*nu)*lsphi
    !
    ! Regular integrals
    !
    do ksub=1,nsub
      ! Numerical integration
      do kip=1,gl01_n(gln)
        ! XIP COORDINATE
        xip=gl01_xi(kip,gln)
        w=gl01_w(kip,gln)
        ! XIP->XI TRANSFORMATION
        js=xisub(2,ksub)-xisub(1,ksub)
        xi=js*xip+xisub(1,ksub)
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
        ! Normal vector
        N(1)=T(2)
        N(2)=-T(1)
        ! xi->x jacobian
        jg=sqrt(dot_product(T,T))
        ! Unit normal vector
        n=N/jg
        ! Unit tangent vector
        t=T/jg
        ! Distance vector and other derived terms
        rv=x-x_i
        r=dot_product(rv,rv)
        r=sqrt(r)
        d1r=1.d0/r
        logd1r=log(d1r)
        drdx=rv*d1r
        drdn=dot_product(drdx,n)
        drdt=dot_product(drdx,t)
        ! Jacobian * weight
        jw=jg*js*w
        ! Complementary fundamental solutions
        call fbem_bem_staela2d_hfc_sbie(mu,nu,hsn,hsx,x,n,x_i,.false.,.true.,fs_u,fs_t)
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
        ! Functional shape functions * jacobians* weights
        pphijw=pphi*jw
        sphijw=sphi*jw
        ! Add integration points (regular part)
        do l=1,2
          do k=1,2
            h(:,l,k)=h(:,l,k)+fs_t(l,k)*pphijw
            g(:,l,k)=g(:,l,k)+fs_u(l,k)*sphijw
          end do
        end do
        ! Regularized part of the CPV integrals
        h(:,1,2)=h(:,1,2)-ctet2*d1r*drdt*(pphi-pphi_i)*jw
        h(:,2,1)=h(:,2,1)+ctet2*d1r*drdt*(pphi-pphi_i)*jw
      end do
    end do
    !
    ! Analytical integrals
    !
    ra=sqrt((e%x(1,1)-x_i(1))**2+(e%x(2,1)-x_i(2))**2)
    rb=sqrt((e%x(1,2)-x_i(1))**2+(e%x(2,2)-x_i(2))**2)
    select case (nsub)
      case (1)
        if (xi_i(1).lt.0.0d0) then
          h(:,1,2)=h(:,1,2)-ctet2*pphi_i*log(rb)
          h(:,2,1)=h(:,2,1)+ctet2*pphi_i*log(rb)
        end if
        if (xi_i(1).gt.0.0d0) then
          h(:,1,2)=h(:,1,2)+ctet2*pphi_i*log(ra)
          h(:,2,1)=h(:,2,1)-ctet2*pphi_i*log(ra)
        end if
      case (2)
        h(:,1,2)=h(:,1,2)-ctet2*pphi_i*(log(rb)-log(ra))
        h(:,2,1)=h(:,2,1)+ctet2*pphi_i*(log(rb)-log(ra))
    end select
    ! Multiply by constants
    h=ctet1*h
    g=cteu1*g
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_staela2d_hfc_sbie_int

  subroutine fbem_bem_staela2d_hfc_sbie_auto(e,reverse,x_i,mu,nu,hsn,hsx,qsp,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)   :: e                 !! Integration element
    logical                  :: reverse           !! Reverse orientation
    real(kind=real64)        :: x_i(2)            !! Image collocation point
    real(kind=real64)        :: mu                !! Shear modulus
    real(kind=real64)        :: nu                !! Poisson's ratio
    integer                  :: hsn               !! Unit normal of the half-plane: -1,1,-2,2
    real(kind=real64)        :: hsx               !! Coordinate of the half-plane
    type(fbem_qs_parameters) :: qsp               !! Quasi-singular integration parameters
    integer                  :: ns                !! Maximum level of subdivisions
    real(kind=real64)        :: h(e%n_pnodes,2,2) !! h integration kernel
    real(kind=real64)        :: g(e%n_snodes,2,2) !! g integration kernel
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
        call fbem_bem_staela2d_hfc_sbie_int(e,reverse,barxi,mu,nu,hsn,hsx,h,g)
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
          call fbem_bem_staela2d_hfc_sbie_ext_pre(ps,e,reverse,x_i,mu,nu,hsn,hsx,h,g)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_staela2d_hfc_sbie_ext_adp(e,reverse,xi_s,x_i,mu,nu,hsn,hsx,qsp,1,ns,h,g)
        end if
    end select
  end subroutine fbem_bem_staela2d_hfc_sbie_auto

  ! ================================================================================================================================

end module fbem_bem_staela2d
