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
!! <b>This module implements the calculation of element-wise integrals of Boundary Integral Equations of the 2D elastodynamics.</b>
module fbem_bem_harela2d

  ! Fortran 2003 standard
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_quad_rules
  use fbem_telles_transformation
  use fbem_geometry
  use fbem_quasisingular_integration
  use fbem_bem_general

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all are private
  private

  ! ================================================================================================================================
  !! Data structure for region parameters (material properties and frequency dependant parameters)
  type fbem_bem_harela2d_parameters
    ! Region parameters
    complex(kind=real64) :: lambda
    complex(kind=real64) :: mu
    real(kind=real64)    :: rho
    real(kind=real64)    :: nu
    ! Frequency
    real(kind=real64)    :: omega
    ! Wave propagation speeds
    complex(kind=real64) :: c1, c2
    ! Wavenumbers
    complex(kind=real64) :: k1, k2
    ! Coefficients of the components of the fundamental solutions u*, t*, d* and s*
    complex(kind=real64) :: psi(4), chi(4)
    complex(kind=real64) :: T1(6), T2(6), T3(6)
    complex(kind=real64) :: S1(7), S2(9), S3(8), S4(6), S5(7)
    complex(kind=real64) :: cte_u, cte_t, cte_d, cte_s
    ! Coefficients of the components of the fundamental solutions u*_{lk,j}, sigma*_{lkm,j}
    complex(kind=real64) :: U1(6), U2(6), U3(5)
    complex(kind=real64) :: R1(6), R2(7), R3(8), R4(6), R5(7), R6(6)
  end type fbem_bem_harela2d_parameters
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! INITIAL SETUP
  !public :: test_harela2d_fs
  public :: fbem_bem_harela2d_parameters
  public :: fbem_bem_harela2d_calculate_parameters
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE)
  ! Free-term
  public :: fbem_bem_harela2d_sbie_freeterm
  ! Fundamental solution
  public :: fbem_bem_harela2d_sbie_u
  public :: fbem_bem_harela2d_sbie_t
  ! BOUNDARY ELEMENTS
  ! Exterior integration
  public :: fbem_bem_harela2d_sbie_ext_pre
  public :: fbem_bem_harela2d_sbie_ext_st
  public :: fbem_bem_harela2d_sbie_ext_adp
  ! Interior integration
  public :: fbem_bem_harela2d_sbie_int
  ! Automatic integration
  public :: fbem_bem_harela2d_sbie_auto
  ! BODY LOADS
  ! Exterior integration
  public :: fbem_bem_harela2d_sbie_bl_ext_pre
  public :: fbem_bem_harela2d_sbie_bl_ext_st
  public :: fbem_bem_harela2d_sbie_bl_ext_adp
  ! Interior integration
  public :: fbem_bem_harela2d_sbie_bl_int
  ! Automatic integration
  public :: fbem_bem_harela2d_sbie_bl_auto
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (HBIE)
  ! Fundamental solution
  public :: fbem_bem_harela2d_hbie_d
  public :: fbem_bem_harela2d_hbie_s
  ! Exterior integration
  public :: fbem_bem_harela2d_hbie_ext_pre
  public :: fbem_bem_harela2d_hbie_ext_st
  public :: fbem_bem_harela2d_hbie_ext_adp
  ! Interior integration
  public :: fbem_bem_harela2d_hbie_int
  ! Automatic integration
  public :: fbem_bem_harela2d_hbie_auto
  ! BODY LOADS
  ! Exterior integration
  public :: fbem_bem_harela2d_hbie_bl_ext_pre
  public :: fbem_bem_harela2d_hbie_bl_ext_st
  public :: fbem_bem_harela2d_hbie_bl_ext_adp
  ! Interior integration
  public :: fbem_bem_harela2d_hbie_bl_int
  ! Automatic integration
  public :: fbem_bem_harela2d_hbie_bl_auto
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! BOTH SBIE AND HBIE SIMULTANEOUSLY
  ! Exterior integration
  public :: fbem_bem_harela2d_shbie_ext_pre
  ! Automatic integration
  public :: fbem_bem_harela2d_shbie_auto
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! VARIATION SINGULAR BOUNDARY INTEGRAL EQUATION (VSBIE)
  public :: fbem_bem_harela2d_vsbie_freeterm
  ! Exterior integration
  public :: fbem_bem_harela2d_vsbie_ext_pre
  public :: fbem_bem_harela2d_vsbie_ext_st
  public :: fbem_bem_harela2d_vsbie_ext_adp
  ! Interior integration
  public :: fbem_bem_harela2d_vsbie_int
  ! Automatic integration
  public :: fbem_bem_harela2d_vsbie_auto
  ! ================================================================================================================================

contains

  ! ================================================================================================================================
  ! INITIAL SETUP

  !! Subroutine that calculate all parameters for given properties and frequency
  subroutine fbem_bem_harela2d_calculate_parameters(lambda,mu,rho,omega,p)
    implicit none
    ! I/O
    complex(kind=real64)               :: lambda
    complex(kind=real64)               :: mu
    real(kind=real64)                  :: rho
    real(kind=real64)                  :: omega
    type(fbem_bem_harela2d_parameters) :: p
    ! Local
    complex(kind=real64)               :: c1, c2, k1, k2, logk1, logk2
    !
    ! Local parameters
    !
    c1=sqrt((lambda+2.d0*mu)/rho)
    c2=sqrt(mu/rho)
    k1=omega/c1
    k2=omega/c2
    logk1=log(0.5d0*c_im*k1)
    logk2=log(0.5d0*c_im*k2)
    !
    ! Save to structure the main parameters
    !
    p%omega=omega
    p%lambda=lambda
    p%mu=mu
    p%rho=rho
    p%nu=dble(0.5d0*lambda/(lambda+mu))
    p%c1=c1
    p%c2=c2
    p%k1=k1
    p%k2=k2
    !
    ! Coefficients of the components of the fundamental solutions u*, t*, d* and s*
    !
    ! psi
    p%psi(1)=-0.5d0*(1.d0+c2**2/c1**2)
    p%psi(2)=-0.5d0*(c_gamma+logk2+0.5d0+k1**2/k2**2*(c_gamma+logk1-0.5d0))
    p%psi(3)=c_im*k1/k2**2
    p%psi(4)=-c_im/k2
    ! chi
    p%chi(1)=0.5d0*(c2**2/c1**2-1.d0)
    p%chi(2)=0.125d0*(k2**2-k1**4/k2**2)
    p%chi(3)=0.125d0*(k2**2*(c_gamma+logk2-0.75d0)-k1**4/k2**2*(c_gamma+logk1-0.75d0))
    p%chi(4)=-k1**2/k2**2
    ! T1
    p%T1(1)=2.d0*(c2**2/c1**2-1.d0)
    p%T1(2)=0.25d0*(k1**4/k2**2-k2**2)
    p%T1(3)=-2.d0*c_im*k1**3/k2**2
    p%T1(4)=2.d0*c_im*k2
    p%T1(5)=-8.d0*k1**2/k2**2
    p%T1(6)=8.d0
    ! T2
    p%T2(1)=-c2**2/c1**2
    p%T2(2)=0.25d0*(k1**4/k2**2+k2**2)
    p%T2(3)=0.25d0*(k1**4/k2**2*(c_gamma+logk1-0.75d0)+k2**2*(c_gamma+logk2-0.25d0))
    p%T2(4)=-c_im*k2
    p%T2(5)=2.d0*k1**2/k2**2
    p%T2(6)=-2.d0
    ! T3
    p%T3(1)=c2**2/c1**2
    p%T3(2)=0.25d0*(2.d0*k1**2-k2**2-3.d0*k1**4/k2**2)
    p%T3(3)=0.25d0*(2.d0*k1**2*(c_gamma+logk1-0.5d0)-k2**2*(c_gamma+logk2-0.75d0)-3.d0*k1**4/k2**2*(c_gamma+logk1-5.d0/12.d0))
    p%T3(4)=c_im*k1*(2.d0*k1**2/k2**2-1.d0)
    p%T3(5)=2.d0*k1**2/k2**2
    p%T3(6)=-2.d0
    ! S1
    p%S1(1)=2.d0*(1.d0-2.d0*c2**2/c1**2)
    p%S1(2)=-0.5d0*k1**4/k2**2
    p%S1(3)=k2**2
    p%S1(4)=4.d0*c_im*k1**3/k2**2
    p%S1(5)=-6.d0*c_im*k2
    p%S1(6)=16.d0*k1**2/k2**2
    p%S1(7)=-16.d0
    ! S2
    p%S2(1)=4.d0*c2**2/c1**2
    p%S2(2)=-0.5d0*(2.d0*k1**2-k2**2-3.d0*k1**4/k2**2)
    p%S2(3)=0.5d0*k1**4*(0.5d0-k1**2/k2**2)
    p%S2(4)=0.5d0*k1**4*(0.5d0-k1**2/k2**2)*(c_gamma+logk1-0.75d0)
    p%S2(5)=4.d0*c_im*k1**3/k2**2
    p%S2(6)=-4.d0*c_im*k2
    p%S2(7)=16.d0*k1**2/k2**2
    p%S2(8)=2.d0*k1**2*(1.d0-2.d0*k1**2/k2**2)
    p%S2(9)=-16.d0
    ! S3
    p%S3(1)=16.d0*(1.d0-c2**2/c1**2)
    p%S3(2)=k2**2-k1**4/k2**2
    p%S3(3)=-4.d0*k1**4/k2**2
    p%S3(4)=4.d0*k2**2
    p%S3(5)=32.d0*c_im*k1**3/k2**2
    p%S3(6)=-32.d0*c_im*k2
    p%S3(7)=96.d0*k1**2/k2**2
    p%S3(8)=-96.d0
    ! S4
    p%S4(1)=2.d0*c2**2/c1**2
    p%S4(2)=-0.5d0*(k2**2+k1**4/k2**2)
    p%S4(3)=-0.5d0*(k2**2*(c_gamma+logk2-0.25d0)+k1**4/k2**2*(c_gamma+logk1-0.75d0))
    p%S4(4)=2.d0*c_im*k2
    p%S4(5)=-4.d0*k1**2/k2**2
    p%S4(6)=4.d0
    ! S5
    p%S5(1)=2.d0*(1.d0-3.d0*c2**2/c1**2)
    p%S5(2)=0.5d0*(4.d0*k1**2-k2**2-k1**4/k2**2)
    p%S5(3)=0.5d0*(4.d0*k1**2*(c_gamma+logk1+0.5d0)-k2**2*(c_gamma+2.d0*logk1-logk2+0.75d0)-k1**4/k2**2*(c_gamma+logk1+13.d0/4.d0))
    p%S5(4)=-4.d0*k1**2+k2**2+4.d0*k1**4/k2**2
    p%S5(5)=4.d0*c_im*k1*(1.d0-2.d0*k1**2/k2**2)
    p%S5(6)=-4.d0*k1**2/k2**2
    p%S5(7)=4.d0
    ! Constant of u*: 1/(2*pi*mu)
    p%cte_u=c_1_2pi/mu
    ! Constant of t*: 1/(2*pi)
    p%cte_t=c_1_2pi
    ! Constant of d*: 1/(2*pi)
    p%cte_d=c_1_2pi
    ! Constant of s*: mu/(2*pi)
    p%cte_s=c_1_2pi*mu
    !
    ! Coefficients of the components of u*_{lk,j} and sigma_{lkm,j}
    !
    ! U1
    p%U1(1)=-0.5d0*(c2**2/c1**2+1.d0)
    p%U1(2)=0.125d0*(k1**4/k2**2+3.d0*k2**2)
    p%U1(3)=0.125d0*(k1**4/k2**2*(logk1+c_gamma-0.75d0)-k2**2*(logk2+c_gamma-0.75d0))+0.5d0*k2**2*(logk2+c_gamma-0.5d0)
    p%U1(4)=-c_im*k2
    p%U1(5)=k1**2/k2**2
    p%U1(6)=-1.d0
    ! U2
    p%U2(1)=c2**2/c1**2-1.d0
    p%U2(2)=0.125d0*(k1**4/k2**2-k2**2)
    p%U2(3)=-k1**2/k2**2*c_im*k1
    p%U2(4)=c_im*k2
    p%U2(5)=-4.d0*k1**2/k2**2
    p%U2(6)=4.d0
    ! U3
    p%U3(1)=0.5d0*(1.d0-c2**2/c1**2)
    p%U3(2)=0.125d0*(k1**4/k2**2-k2**2)
    p%U3(3)=0.125d0*(k1**4/k2**2*logk1-k2**2*logk2+(k1**4/k2**2-k2**2)*(c_gamma-0.75d0))
    p%U3(4)=k1**2/k2**2
    p%U3(5)=-1.d0
    ! R1
    p%R1(1)=2.d0*(c2**2/c1**2-1.d0)
    p%R1(2)=0.25d0*(k1**4/k2**2-k2**2)
    p%R1(3)=-2.d0*k1**2/k2**2*c_im*k1
    p%R1(4)=2.d0*c_im*k2
    p%R1(5)=-8.d0*k1**2/k2**2
    p%R1(6)=8.d0
    ! R2
    p%R2(1)=2.d0*c2**2/c1**2
    p%R2(2)=0.25d0*(k1**4/k2**2+k2**2)
    p%R2(3)=-k2**2
    p%R2(4)=-2.d0*k1**2/k2**2*c_im*k1
    p%R2(5)=4.d0*c_im*k2
    p%R2(6)=-8.d0*k1**2/k2**2
    p%R2(7)=8.d0
    ! R3
    p%R3(1)=8.d0*(1.d0-c2**2/c1**2)
    p%R3(2)=0.5d0*(k2**2-k1**4/k2**2)
    p%R3(3)=-2.d0*k1**4/k2**2
    p%R3(4)=2.d0*k2**2
    p%R3(5)=16.d0*k1**2/k2**2*c_im*k1
    p%R3(6)=-16.d0*c_im*k2
    p%R3(7)=48.d0*k1**2/k2**2
    p%R3(8)=-48.d0
    ! R4
    p%R4(1)=c2**2/c1**2
    p%R4(2)=0.25d0*(2.d0*k1**2-k2**2-3.d0*k1**4/k2**2)
    p%R4(3)=0.25d0*(2.d0*k1**2*(c_gamma+logk1-0.5d0)-k2**2*(c_gamma+logk2-0.75d0)-3.d0*k1**4/k2**2*(c_gamma+logk1-5.d0/12.d0))
    p%R4(4)=c_im*k1*(2.d0*k1**2/k2**2-1.d0)
    p%R4(5)=2.d0*k1**2/k2**2
    p%R4(6)=-2.d0
    ! R5
    p%R5(1)=-2.d0*c2**2/c1**2
    p%R5(2)=0.25d0*(-3.d0*k1**4/k2**2+2.d0*k1**2-k2**2)
    p%R5(3)=2.d0*k1**4/k2**2-k1**2
    p%R5(4)=2.d0*c_im*k1*(1.d0-3.d0*k1**2/k2**2)
    p%R5(5)=2.d0*c_im*k2
    p%R5(6)=-8.d0*k1**2/k2**2
    p%R5(7)=8.d0
    ! R6
    p%R6(1)=-c2**2/c1**2
    p%R6(2)=0.25d0*(k1**4/k2**2+k2**2)
    p%R6(3)=0.25d0*(k1**4/k2**2*(c_gamma+logk1-0.75d0)+k2**2*(c_gamma+logk2-0.25d0))
    p%R6(4)=-c_im*k2
    p%R6(5)=2.d0*k1**2/k2**2
    p%R6(6)=-2.d0
  end subroutine fbem_bem_harela2d_calculate_parameters

!  !! Test of fundamental solutions
!  subroutine test_harela2d_fs
!    implicit none
!    complex(kind=real64)               :: lambda
!    complex(kind=real64)               :: mu
!    complex(kind=real64)               :: k1, k2, c1, c2
!    real(kind=real64)                  :: rho
!    type(fbem_bem_harela2d_parameters) :: p
!    real(kind=real64)                  :: omega, o_min, o_max, deltao
!    integer                            :: ko, no
!    real(kind=real64)                  :: r, r_min, r_max, deltar
!    integer                            :: kr, nr
!    real(kind=real64)                  :: logr, r2, dr1, dr2, dr3, dr4, dr5
!    complex(kind=real64)               :: z(2), Knr(0:2,2), K0z1, K1z1, K2z1, K0z2, K1z2, K2z2
!    complex(kind=real64)               :: psi, chi, dpsidr, dchidr, d2psid2r, d2chid2r
!    complex(kind=real64)               :: T1, T2, T3, S1, S2, S3, S4, S5
!    complex(kind=real64)               :: U1, U2, U3, RR1, RR2, RR3, RR4, RR5, RR6
!    ! Files
!    open(unit=10,file='highkr.dat')
!    open(unit=11,file='lowkr.dat')
!    ! Parameters
!    lambda=100000.d0
!    mu=100000.d0
!    rho=1.d0
!    ! omega setup
!    o_min=1.d-12
!    o_max=1.d12
!    no=25
!    ! r setup
!    r_min=1.d-12
!    r_max=1.d12
!    nr=1000
!    ! Delta of omega and r
!    deltao=(log10(o_max)-log10(o_min))/dble(no-1)
!    deltar=(log10(r_max)-log10(r_min))/dble(nr-1)
!    ! Loop through omega
!    do ko=1,no
!      omega=10**(dlog10(o_min)+deltao*dble(ko-1))
!      call fbem_bem_harela2d_calculate_parameters(lambda,mu,rho,omega,p)
!      k1=p%k1
!      k2=p%k2
!      c1=p%c1
!      c2=p%c2
!      ! Loop through r
!      do kr=1,nr
!        r=10**(dlog10(r_min)+deltar*dble(kr-1))
!        logr=log(r)
!        r2=r**2
!        dr1=1.d0/r
!        dr2=1.d0/r**2
!        dr3=1.d0/r**3
!        dr4=1.d0/r**4
!        dr5=1.d0/r**5
!        z(1)=c_im*k1*r
!        z(2)=c_im*k2*r
!        ! Compact form of the components. Almost all of them break for k·r << 1
!        call fbem_modified_bessel_K0_and_K1(z(1), K0z1, K1z1)
!        K2z1=K0z1+2.d0/z(1)*K1z1
!        call fbem_modified_bessel_K0_and_K1(z(2), K0z2, K1z2)
!        K2z2=K0z2+2.d0/z(2)*K1z2
!        psi=K0z2+1.d0/z(2)*(K1z2-k1/k2*K1z1)
!        chi=K2z2-k1**2/k2**2*K2z1
!        dpsidr=-c_im*k2*K1z2-dr1*chi
!        dchidr=c_im*k1*k1**2/k2**2*K1z1-c_im*k2*K1z2-2.d0*dr1*chi
!        d2psid2r=-k2**2*K0z2+c_im*k2*dr1*K1z2+dr1*(-dchidr+dr1*chi)
!        d2chid2r=k1**4/k2**2*(K0z1+1.d0/z(1)*K1z1)-k2**2*(K0z2+1.d0/z(2)*K1z2)+2.d0*dr1*(dr1*chi-dchidr)
!        T1=-2.d0*(dchidr-2.d0*dr1*chi)
!        T2=dpsidr-dr1*chi
!        T3=-2.d0*dr1*chi+(c1**2/c2**2-2.d0)*(dpsidr-dchidr-dr1*chi)
!        S1=-d2psid2r+dr1*dpsidr+3.d0*dr1*dchidr-6.d0*dr2*chi
!        S2=-2.d0*dr1*T1+2.d0*(c1**2/c2**2-2.d0)*(-d2psid2r+d2chid2r+dr1*(dpsidr-2.d0*dr1*chi))
!        S3=-4.d0*(d2chid2r-dr1*(5.d0*dchidr-8.d0*dr1*chi))
!        S4=-2.d0*dr1*T2
!        S5=4.d0*dr2*chi+(c1**2/c2**2-2.d0)*4.d0*dr1*(-dpsidr+dchidr+dr1*chi)&
!          +(c1**2/c2**2-2.d0)**2*(-d2psid2r+d2chid2r-dr1*dpsidr+2.d0*dr1*dchidr)
!        U1=dpsidr
!        U2=2.*dr1*chi-dchidr
!        U3=-dr1*chi
!        RR1=2.*dr1*(-dchidr+2.*dr1*chi)
!        RR2=d2psid2r-dr1*(dpsidr+dchidr-2.*dr1*chi)
!        RR3=2.*(-d2chid2r+dr1*(5.*dchidr-8.*dr1*chi))
!        RR4=dr1*(lambda/mu*(dpsidr-dchidr-dr1*chi)-2.*dr1*chi)
!        RR5=lambda/mu*(d2psid2r-d2chid2r-dr1*(dpsidr-2.*dr1*chi))-2.*dr1*(dchidr-2.*dr1*chi)
!        RR6=dr1*(dpsidr-dr1*chi)
!        write(10,'(i11,e25.16,i11,41e25.16)') ko,omega,kr,r,dreal(k1),dreal(k2),dreal(psi),dimag(psi),dreal(chi),dimag(chi),&
!                                                                                dreal( T1),dimag( T1),dreal( T2),dimag( T2),&
!                                                                                dreal( T3),dimag( T3),dreal( S1),dimag( S1),&
!                                                                                dreal( S2),dimag( S2),dreal( S3),dimag( S3),&
!                                                                                dreal( S4),dimag( S4),dreal( S5),dimag( S5),&
!                                                                                dreal( U1),dimag( U1),dreal( U2),dimag( U2),&
!                                                                                dreal( U3),dimag( U3),dreal(RR1),dimag(RR1),&
!                                                                                dreal(RR2),dimag(RR2),dreal(RR3),dimag(RR3),&
!                                                                                dreal(RR4),dimag(RR4),dreal(RR5),dimag(RR5),&
!                                                                                dreal(RR6),dimag(RR6)
!        ! Expanded form of the components. Mandatory for k·r << 1, some of them break for k·r >> 1 (for very very large k·r).
!        call fbem_BesselKnR_decomposed(2,z,KnR)
!        psi=p%psi(1)*logr+p%psi(2)+KnR(0,2)+dr1*(p%psi(3)*KnR(1,1)+p%psi(4)*KnR(1,2))
!        chi=p%chi(1)+(p%chi(2)*logr+p%chi(3))*r2+p%chi(4)*KnR(2,1)+KnR(2,2)
!        T1=p%T1(1)*dr1+p%T1(2)*r+p%T1(3)*KnR(1,1)+p%T1(4)*KnR(1,2)+p%T1(5)*dr1*KnR(2,1)+p%T1(6)*dr1*KnR(2,2)
!        T2=p%T2(1)*dr1+(p%T2(2)*logr+p%T2(3))*r+p%T2(4)*KnR(1,2)+p%T2(5)*dr1*KnR(2,1)+p%T2(6)*dr1*KnR(2,2)
!        T3=p%T3(1)*dr1+(p%T3(2)*logr+p%T3(3))*r+p%T3(4)*KnR(1,1)+p%T3(5)*dr1*KnR(2,1)+p%T3(6)*dr1*KnR(2,2)
!        S1=p%S1(1)*dr2+p%S1(2)+p%S1(3)*KnR(0,2)+p%S1(4)*dr1*KnR(1,1)+p%S1(5)*dr1*KnR(1,2)+p%S1(6)*dr2*KnR(2,1)+p%S1(7)*dr2*KnR(2,2)
!        S2=p%S2(1)*dr2+p%S2(2)+(p%S2(3)*logr+p%S2(4))*r2+p%S2(5)*dr1*KnR(1,1)+p%S2(6)*dr1*KnR(1,2)+p%S2(7)*dr2*KnR(2,1)+p%S2(8)*KnR(2,1)+p%S2(9)*dr2*KnR(2,2)
!        S3=p%S3(1)*dr2+p%S3(2)+p%S3(3)*KnR(0,1)+p%S3(4)*KnR(0,2)+p%S3(5)*dr1*KnR(1,1)+p%S3(6)*dr1*KnR(1,2)+p%S3(7)*dr2*KnR(2,1)+p%S3(8)*dr2*KnR(2,2)
!        S4=p%S4(1)*dr2+p%S4(2)*logr+p%S4(3)+p%S4(4)*dr1*KnR(1,2)+p%S4(5)*dr2*KnR(2,1)+p%S4(6)*dr2*KnR(2,2)
!        S5=p%S5(1)*dr2+p%S5(2)*logr+p%S5(3)+p%S5(4)*KnR(0,1)+p%S5(5)*dr1*KnR(1,1)+p%S5(6)*dr2*KnR(2,1)+p%S5(7)*dr2*KnR(2,2)
!        U1=p%U1(1)*dr1+(p%U1(2)*logr+p%U1(3))*r+p%U1(4)*KnR(1,2)+p%U1(5)*dr1*KnR(2,1)+p%U1(6)*dr1*KnR(2,2)
!        U2=p%U2(1)*dr1+p%U2(2)*r+p%U2(3)*KnR(1,1)+p%U2(4)*KnR(1,2)+p%U2(5)*dr1*KnR(2,1)+p%U2(6)*dr1*KnR(2,2)
!        U3=p%U3(1)*dr1+(p%U3(2)*logr+p%U3(3))*r+p%U3(4)*dr1*KnR(2,1)+p%U3(5)*dr1*KnR(2,2)
!        RR1=p%R1(1)*dr2+p%R1(2)+p%R1(3)*dr1*KnR(1,1)+p%R1(4)*dr1*KnR(1,2)+p%R1(5)*dr2*KnR(2,1)+p%R1(6)*dr2*KnR(2,2)
!        RR2=p%R2(1)*dr2+p%R2(2)+p%R2(3)*KnR(0,2)+p%R2(4)*dr1*KnR(1,1)+p%R2(5)*dr1*KnR(1,2)+p%R2(6)*dr2*KnR(2,1)+p%R2(7)*dr2*KnR(2,2)
!        RR3=p%R3(1)*dr2+p%R3(2)+p%R3(3)*KnR(0,1)+p%R3(4)*KnR(0,2)+p%R3(5)*dr1*KnR(1,1)+p%R3(6)*dr1*KnR(1,2)+p%R3(7)*dr2*KnR(2,1)+p%R3(8)*dr2*KnR(2,2)
!        RR4=p%R4(1)*dr2+p%R4(2)*logr+p%R4(3)+p%R4(4)*dr1*KnR(1,1)+p%R4(5)*dr2*KnR(2,1)+p%R4(6)*dr2*KnR(2,2)
!        RR5=p%R5(1)*dr2+p%R5(2)+p%R5(3)*KnR(0,1)+p%R5(4)*dr1*KnR(1,1)+p%R5(5)*dr1*KnR(1,2)+p%R5(6)*dr2*KnR(2,1)+p%R5(7)*dr2*KnR(2,2)
!        RR6=p%R6(1)*dr2+p%R6(2)*logr+p%R6(3)+p%R6(4)*dr1*KnR(1,2)+p%R6(5)*dr2*KnR(2,1)+p%R6(6)*dr2*KnR(2,2)
!        write(11,'(i11,e25.16,i11,41e25.16)') ko,omega,kr,r,dreal(k1),dreal(k2),dreal(psi),dimag(psi),dreal(chi),dimag(chi),&
!                                                                                dreal( T1),dimag( T1),dreal( T2),dimag( T2),&
!                                                                                dreal( T3),dimag( T3),dreal( S1),dimag( S1),&
!                                                                                dreal( S2),dimag( S2),dreal( S3),dimag( S3),&
!                                                                                dreal( S4),dimag( S4),dreal( S5),dimag( S5),&
!                                                                                dreal( U1),dimag( U1),dreal( U2),dimag( U2),&
!                                                                                dreal( U3),dimag( U3),dreal(RR1),dimag(RR1),&
!                                                                                dreal(RR2),dimag(RR2),dreal(RR3),dimag(RR3),&
!                                                                                dreal(RR4),dimag(RR4),dreal(RR5),dimag(RR5),&
!                                                                                dreal(RR6),dimag(RR6)
!      end do
!    end do
!    !
!    ! GNUPLOT FILE SCRIPT
!!    !
!!set terminal pdfcairo enhanced size 29.7cm,21.0cm
!!set log xy
!!set output "test_harela2d_fs.pdf"
!!set xlabel "k_2*r"
!!set title "Re(psi)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0):(abs($7)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==1?$6*$4:1/0):(abs($7)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($7)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==12?$6*$4:1/0):(abs($7)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($7)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==25?$6*$4:1/0):(abs($7)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Im(psi)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0):(abs($8)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==1?$6*$4:1/0):(abs($8)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($8)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==12?$6*$4:1/0):(abs($8)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($8)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==25?$6*$4:1/0):(abs($8)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Re(chi)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0):(abs($9)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==1?$6*$4:1/0):(abs($9)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($9)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==12?$6*$4:1/0):(abs($9)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($9)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==25?$6*$4:1/0):(abs($9)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Im(chi)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0):(abs($10)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==1?$6*$4:1/0):(abs($10)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($10)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==12?$6*$4:1/0):(abs($10)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($10)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==25?$6*$4:1/0):(abs($10)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Re(T1)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0):(abs($11)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==1?$6*$4:1/0):(abs($11)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($11)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==12?$6*$4:1/0):(abs($11)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($11)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==25?$6*$4:1/0):(abs($11)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Im(T1)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0):(abs($12)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==1?$6*$4:1/0):(abs($12)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($12)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==12?$6*$4:1/0):(abs($12)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($12)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==25?$6*$4:1/0):(abs($12)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Re(T2)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0):(abs($13)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==1?$6*$4:1/0):(abs($13)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($13)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==12?$6*$4:1/0):(abs($13)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($13)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==25?$6*$4:1/0):(abs($13)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Im(T2)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0):(abs($14)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==1?$6*$4:1/0):(abs($14)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($14)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==12?$6*$4:1/0):(abs($14)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($14)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==25?$6*$4:1/0):(abs($14)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Re(T3)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0):(abs($15)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==1?$6*$4:1/0):(abs($15)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($15)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==12?$6*$4:1/0):(abs($15)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($15)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==25?$6*$4:1/0):(abs($15)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Im(T3)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0):(abs($16)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==1?$6*$4:1/0):(abs($16)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($16)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==12?$6*$4:1/0):(abs($16)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($16)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==25?$6*$4:1/0):(abs($16)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Re(S1)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0):(abs($17)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==1?$6*$4:1/0):(abs($17)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($17)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==12?$6*$4:1/0):(abs($17)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($17)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==25?$6*$4:1/0):(abs($17)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Im(S1)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0):(abs($18)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==1?$6*$4:1/0):(abs($18)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($18)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==12?$6*$4:1/0):(abs($18)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($18)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==25?$6*$4:1/0):(abs($18)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Re(S2)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0):(abs($19)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==1?$6*$4:1/0):(abs($19)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($19)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==12?$6*$4:1/0):(abs($19)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($19)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==25?$6*$4:1/0):(abs($19)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Im(S2)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0):(abs($20)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==1?$6*$4:1/0):(abs($20)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($20)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==12?$6*$4:1/0):(abs($20)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($20)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==25?$6*$4:1/0):(abs($20)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Re(S3)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0):(abs($21)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==1?$6*$4:1/0):(abs($21)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($21)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==12?$6*$4:1/0):(abs($21)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($21)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==25?$6*$4:1/0):(abs($21)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Im(S3)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0):(abs($22)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==1?$6*$4:1/0):(abs($22)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($22)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==12?$6*$4:1/0):(abs($22)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($22)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==25?$6*$4:1/0):(abs($22)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Re(S4)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0):(abs($23)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==1?$6*$4:1/0):(abs($23)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($23)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==12?$6*$4:1/0):(abs($23)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($23)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==25?$6*$4:1/0):(abs($23)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Im(S4)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0):(abs($24)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==1?$6*$4:1/0):(abs($24)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($24)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==12?$6*$4:1/0):(abs($24)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($24)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==25?$6*$4:1/0):(abs($24)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Re(S5)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0):(abs($25)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==1?$6*$4:1/0):(abs($25)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($25)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==12?$6*$4:1/0):(abs($25)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($25)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==25?$6*$4:1/0):(abs($25)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Im(S5)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0):(abs($26)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==1?$6*$4:1/0):(abs($26)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($26)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==12?$6*$4:1/0):(abs($26)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($26)) w l t "k*r>>1", \
!!     "./lowkr.dat" u ($1==25?$6*$4:1/0):(abs($26)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Re(U1)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0 ):(abs($27)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==1?$6*$4:1/0 ):(abs($27)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($27)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==12?$6*$4:1/0):(abs($27)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($27)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==25?$6*$4:1/0):(abs($27)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Im(U1)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0 ):(abs($28)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==1?$6*$4:1/0 ):(abs($28)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($28)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==12?$6*$4:1/0):(abs($28)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($28)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==25?$6*$4:1/0):(abs($28)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Re(U2)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0 ):(abs($29)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==1?$6*$4:1/0 ):(abs($29)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($29)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==12?$6*$4:1/0):(abs($29)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($29)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==25?$6*$4:1/0):(abs($29)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Im(U2)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0 ):(abs($30)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==1?$6*$4:1/0 ):(abs($30)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($30)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==12?$6*$4:1/0):(abs($30)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($30)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==25?$6*$4:1/0):(abs($30)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Re(U3)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0 ):(abs($31)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==1?$6*$4:1/0 ):(abs($31)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($31)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==12?$6*$4:1/0):(abs($31)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($31)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==25?$6*$4:1/0):(abs($31)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Im(U3)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0 ):(abs($32)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==1?$6*$4:1/0 ):(abs($32)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($32)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==12?$6*$4:1/0):(abs($32)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($32)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==25?$6*$4:1/0):(abs($32)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Re(RR1)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0 ):(abs($33)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==1?$6*$4:1/0 ):(abs($33)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($33)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==12?$6*$4:1/0):(abs($33)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($33)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==25?$6*$4:1/0):(abs($33)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Im(RR1)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0 ):(abs($34)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==1?$6*$4:1/0 ):(abs($34)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($34)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==12?$6*$4:1/0):(abs($34)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($34)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==25?$6*$4:1/0):(abs($34)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Re(RR2)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0 ):(abs($35)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==1?$6*$4:1/0 ):(abs($35)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($35)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==12?$6*$4:1/0):(abs($35)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($35)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==25?$6*$4:1/0):(abs($35)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Im(RR2)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0 ):(abs($36)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==1?$6*$4:1/0 ):(abs($36)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($36)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==12?$6*$4:1/0):(abs($36)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($36)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==25?$6*$4:1/0):(abs($36)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Re(RR3)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0 ):(abs($37)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==1?$6*$4:1/0 ):(abs($37)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($37)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==12?$6*$4:1/0):(abs($37)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($37)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==25?$6*$4:1/0):(abs($37)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Im(RR3)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0 ):(abs($38)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==1?$6*$4:1/0 ):(abs($38)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($38)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==12?$6*$4:1/0):(abs($38)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($38)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==25?$6*$4:1/0):(abs($38)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Re(RR4)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0 ):(abs($39)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==1?$6*$4:1/0 ):(abs($39)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($39)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==12?$6*$4:1/0):(abs($39)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($39)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==25?$6*$4:1/0):(abs($39)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Im(RR4)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0 ):(abs($40)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==1?$6*$4:1/0 ):(abs($40)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($40)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==12?$6*$4:1/0):(abs($40)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($40)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==25?$6*$4:1/0):(abs($40)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Re(RR5)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0 ):(abs($41)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==1?$6*$4:1/0 ):(abs($41)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($41)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==12?$6*$4:1/0):(abs($41)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($41)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==25?$6*$4:1/0):(abs($41)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Im(RR5)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0 ):(abs($42)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==1?$6*$4:1/0 ):(abs($42)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($42)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==12?$6*$4:1/0):(abs($42)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($42)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==25?$6*$4:1/0):(abs($42)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Re(RR6)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0 ):(abs($43)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==1?$6*$4:1/0 ):(abs($43)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($43)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==12?$6*$4:1/0):(abs($43)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($43)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==25?$6*$4:1/0):(abs($43)) w l t "k*r<<1 (omega=10^{12})"
!!set title "Im(RR6)"
!!plot "./highkr.dat" u ($1==1?$6*$4:1/0 ):(abs($44)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==1?$6*$4:1/0 ):(abs($44)) w l t "k*r<<1 (omega=10^{-12})", \
!!     "./highkr.dat" u ($1==12?$6*$4:1/0):(abs($44)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==12?$6*$4:1/0):(abs($44)) w l t "k*r<<1 (omega=1)", \
!!     "./highkr.dat" u ($1==25?$6*$4:1/0):(abs($44)) w l t "k*r>>1", \
!!     "./lowkr.dat"  u ($1==25?$6*$4:1/0):(abs($44)) w l t "k*r<<1 (omega=10^{12})"
!!set output
!  end subroutine test_harela2d_fs
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE)

  !! Subroutine that calculates free-term of 2D potential problem
  subroutine fbem_bem_harela2d_sbie_freeterm(n_elements,n,tc,tol,nu,c)
    implicit none
    ! I/O
    integer              :: n_elements    !! Number of elements (1 or 2)
    real(kind=real64)    :: n(2,2)        !! Normals of elements at collocation point
    real(kind=real64)    :: tc(2,2)       !! Tangents of elements at collocation point towards inside the element
    real(kind=real64)    :: tol           !! Geometric tolerance [1.0e-12,1.0e-3] (default 1.0e-6)
    complex(kind=real64) :: nu            !! Poisson's ratio
    complex(kind=real64) :: c(2,2)        !! Free-term
    ! Local
    real(kind=real64)    :: local_n(2,2)  ! Local copy of n
    real(kind=real64)    :: local_tc(2,2) ! Local copy of tc
    real(kind=real64)    :: local_tol     ! Local copy of geometric tolerance
    real(kind=real64)    :: nsum(2)       ! Unit sum of normals
    real(kind=real64)    :: vectornorm    ! Norm of a vector
    real(kind=real64)    :: theta(2)      ! Angles of unit tangents
    real(kind=real64)    :: theta_ext     ! Angle theta exterior
    complex(kind=real64) :: cte
    real(kind=real64)    :: ctesin   ! Auxiliary variable
    integer              :: i             ! Counter
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
  end subroutine fbem_bem_harela2d_sbie_freeterm

  !! Fundamental solution u*
  subroutine fbem_bem_harela2d_sbie_u(x,x_i,p,uo)
    implicit none
    ! I/O
    real(kind=real64)                  :: x(2)    !! Observation point
    real(kind=real64)                  :: x_i(2)  !! Collocation point
    type(fbem_bem_harela2d_parameters) :: p       !! Parameters of the region
    complex(kind=real64)               :: uo(2,2) !! u*_{lk}
    ! Local
    integer              :: il, ik
    real(kind=real64)    :: rv(2), r, r2, d1r1, logr, drdx(2)
    complex(kind=real64) :: z(2), KnR(0:2,2)
    complex(kind=real64) :: psi, chi
    rv=x-x_i
    r=sqrt(dot_product(rv,rv))
    r2=r**2
    d1r1=1.0d0/r
    logr=log(r)
    drdx=rv*d1r1
    z(1)=c_im*p%k1*r
    z(2)=c_im*p%k2*r
    call fbem_BesselKnR_decomposed(2,z,KnR)
    psi=p%psi(1)*logr+p%psi(2)+KnR(0,2)+d1r1*p%psi(3)*KnR(1,1)+d1r1*p%psi(4)*KnR(1,2)
    chi=p%chi(1)+(p%chi(2)*logr+p%chi(3))*r2+p%chi(4)*KnR(2,1)+KnR(2,2)
    do il=1,2
      do ik=1,2
        uo(il,ik)=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
      end do
    end do
    uo=p%cte_u*uo
  end subroutine fbem_bem_harela2d_sbie_u

  !! Fundamental solution t*
  subroutine fbem_bem_harela2d_sbie_t(x,n,x_i,p,to)
    implicit none
    ! I/O
    real(kind=real64)                  :: x(2)    !! Observation point
    real(kind=real64)                  :: n(2)    !! Observation point normal
    real(kind=real64)                  :: x_i(2)  !! Collocation point
    type(fbem_bem_harela2d_parameters) :: p       !! Parameters of the region
    complex(kind=real64)               :: to(2,2) !! t*_{lk}
    ! Local
    integer              :: il, ik
    real(kind=real64)    :: rv(2), r, r2, d1r1, logr, drdx(2), drdn
    complex(kind=real64) :: z(2), KnR(0:2,2)
    complex(kind=real64) :: T1, T2, T3
    rv=x-x_i
    r=sqrt(dot_product(rv,rv))
    r2=r**2
    d1r1=1.0d0/r
    logr=log(r)
    drdx=rv*d1r1
    drdn=dot_product(drdx,n)
    z(1)=c_im*p%k1*r
    z(2)=c_im*p%k2*r
    call fbem_BesselKnR_decomposed(2,z,KnR)
    T1=p%T1(1)*d1r1+p%T1(2)*r+p%T1(3)*KnR(1,1)+p%T1(4)*KnR(1,2)+p%T1(5)*d1r1*KnR(2,1)+p%T1(6)*d1r1*KnR(2,2)
    T2=p%T2(1)*d1r1+(p%T2(2)*logr+p%T2(3))*r+p%T2(4)*KnR(1,2)+p%T2(5)*d1r1*KnR(2,1)+p%T2(6)*d1r1*KnR(2,2)
    T3=p%T3(1)*d1r1+(p%T3(2)*logr+p%T3(3))*r+p%T3(4)*KnR(1,1)+p%T3(5)*d1r1*KnR(2,1)+p%T3(6)*d1r1*KnR(2,2)
    do il=1,2
      do ik=1,2
        to(il,ik)=T1*drdx(il)*drdx(ik)*drdn+T2*(drdn*c_dkr(il,ik)+drdx(ik)*n(il))+T3*drdx(il)*n(ik)
      end do
    end do
    to=p%cte_t*to
  end subroutine fbem_bem_harela2d_sbie_t

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! BOUNDARY ELEMENTS

  subroutine fbem_bem_harela2d_sbie_ext_pre(ps,e,reverse,x_i,p,h,g)
    implicit none
    ! I/O
    integer                            :: ps                !! Selected precalculated dataset
    type(fbem_bem_element)             :: e                 !! Element
    logical                            :: reverse           !! Reverse normal vector
    real(kind=real64)                  :: x_i(2)            !! Collocation point position vector
    type(fbem_bem_harela2d_parameters) :: p                 !! Parameters of the region
    complex(kind=real64)               :: h(e%n_pnodes,2,2) !! h integration kernels matrix
    complex(kind=real64)               :: g(e%n_snodes,2,2) !! g integration kernels matrix
    ! Local
    integer              :: kip                  ! Counter variable for integration points loop
    real(kind=real64)    :: x(2)                 ! Position vector at integration point
    real(kind=real64)    :: n(2)                 ! Unit normal vector at integration point
    real(kind=real64)    :: pphijw(e%n_pnodes)   ! phi^p * jacobian * weight at integration point
    real(kind=real64)    :: sphijw(e%n_snodes)   ! phi^s * jacobian * weight at integration point
    real(kind=real64)    :: rv(2)                ! Radiovector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, r2, d1r1, logr    ! Radiovector module, squared, inverses and log(r)
    real(kind=real64)    :: drdx(2)              ! Radiovector derivatives with respect to x_k
    real(kind=real64)    :: drdn                 ! Partial derivative of r respect to unit normal
    integer              :: il, ik               ! Counter variables for load direction and observation direction
    complex(kind=real64) :: z(2)                 ! Bessel functions arguments: z=ikr
    complex(kind=real64) :: KnR(0:2,2)           ! Bessel functions decomposition
    complex(kind=real64) :: psi, chi, T1, T2, T3 ! Fundamental solution components
    complex(kind=real64) :: fs_u, fs_t           ! Fundamental solutions values
    ! Initialization
    h=(0.d0,0.d0)
    g=(0.d0,0.d0)
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      ! Functional shape functions * jacobian * weight
      pphijw=e%ps_pphijw(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      ! Components of fundamental solutions
      x=e%ps_x(:,kip,ps)
      n=e%ps_n(:,kip,ps)
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      r2=r**2
      d1r1=1.0d0/r
      logr=log(r)
      drdx=rv*d1r1
      drdn=dot_product(drdx,n)
      z(1)=c_im*p%k1*r
      z(2)=c_im*p%k2*r
      call fbem_BesselKnR_decomposed(2,z,KnR)
      psi=p%psi(1)*logr+p%psi(2)+KnR(0,2)+d1r1*p%psi(3)*KnR(1,1)+d1r1*p%psi(4)*KnR(1,2)
      chi=p%chi(1)+(p%chi(2)*logr+p%chi(3))*r2+p%chi(4)*KnR(2,1)+KnR(2,2)
      T1=p%T1(1)*d1r1+p%T1(2)*r+p%T1(3)*KnR(1,1)+p%T1(4)*KnR(1,2)+p%T1(5)*d1r1*KnR(2,1)+p%T1(6)*d1r1*KnR(2,2)
      T2=p%T2(1)*d1r1+(p%T2(2)*logr+p%T2(3))*r+p%T2(4)*KnR(1,2)+p%T2(5)*d1r1*KnR(2,1)+p%T2(6)*d1r1*KnR(2,2)
      T3=p%T3(1)*d1r1+(p%T3(2)*logr+p%T3(3))*r+p%T3(4)*KnR(1,1)+p%T3(5)*d1r1*KnR(2,1)+p%T3(6)*d1r1*KnR(2,2)
      ! Add the term of the numerical integration summation (fundamental solutions * functional shape functions * jacobian * weight)
      do il=1,2
        do ik=1,2
          ! Fundamental solutions
          fs_u=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
          fs_t=T1*drdx(il)*drdx(ik)*drdn+T2*(drdn*c_dkr(il,ik)+drdx(ik)*n(il))+T3*drdx(il)*n(ik)
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
  end subroutine fbem_bem_harela2d_sbie_ext_pre

  subroutine fbem_bem_harela2d_sbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,p,gln,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                 !! Integration element
    logical                            :: reverse           !! Reverse normal vector
    real(kind=real64)                  :: xi_s(1,2)         !! Coordinates of the subdivision of the parent element (xi space [-1,1])
    real(kind=real64)                  :: x_i(2)            !! Collocation point position vector
    real(kind=real64)                  :: barxip(1)         !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)                  :: barr              !! Telles jacobian at barxip
    type(fbem_bem_harela2d_parameters) :: p                 !! Parameters of the region
    integer                            :: gln               !! 1D Gauss-Legendre number of integration points (<=32)
    complex(kind=real64)               :: h(e%n_pnodes,2,2) !! h kernel vector
    complex(kind=real64)               :: g(e%n_snodes,2,2) !! g kernel vector
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
    real(kind=real64)            :: r, r2, d1r, logr     ! Radiovector module, squared, its inverses and log(r)
    real(kind=real64)            :: drdx(2)              ! Radiovector derivatives with respect to x_k
    real(kind=real64)            :: jg                   ! Geometric jacobian
    real(kind=real64)            :: drdn                 ! Partial derivative of r respect to unit normal
    real(kind=real64)            :: jw                   ! Jacobian * weight
    real(kind=real64)            :: pphijw(e%n_pnodes)   ! Auxiliary variables for integrand evaluation
    real(kind=real64)            :: sphijw(e%n_snodes)   ! Auxiliary variables for integrand evaluation
    integer                      :: il, ik               ! Counter variables for load direction and observation direction
    complex(kind=real64)         :: z(2)                 ! Bessel functions arguments: z=ikr
    complex(kind=real64)         :: KnR(0:2,2)           ! Bessel functions decomposition
    complex(kind=real64)         :: psi, chi             ! Fundamental solution components
    complex(kind=real64)         :: T1, T2, T3           ! Fundamental solution components
    complex(kind=real64)         :: fs_u, fs_t           ! Fundamental solutions values
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
      r2=r**2
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
      z(1)=c_im*p%k1*r
      z(2)=c_im*p%k2*r
      call fbem_BesselKnR_decomposed(2,z,KnR)
      psi=p%psi(1)*logr+p%psi(2)+KnR(0,2)+d1r*p%psi(3)*KnR(1,1)+d1r*p%psi(4)*KnR(1,2)
      chi=p%chi(1)+(p%chi(2)*logr+p%chi(3))*r2+p%chi(4)*KnR(2,1)+KnR(2,2)
      T1=p%T1(1)*d1r+p%T1(2)*r+p%T1(3)*KnR(1,1)+p%T1(4)*KnR(1,2)+p%T1(5)*d1r*KnR(2,1)+p%T1(6)*d1r*KnR(2,2)
      T2=p%T2(1)*d1r+(p%T2(2)*logr+p%T2(3))*r+p%T2(4)*KnR(1,2)+p%T2(5)*d1r*KnR(2,1)+p%T2(6)*d1r*KnR(2,2)
      T3=p%T3(1)*d1r+(p%T3(2)*logr+p%T3(3))*r+p%T3(4)*KnR(1,1)+p%T3(5)*d1r*KnR(2,1)+p%T3(6)*d1r*KnR(2,2)
      ! Add
      do il=1,2
        do ik=1,2
          fs_u=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
          fs_t=T1*drdx(il)*drdx(ik)*drdn+T2*(drdn*c_dkr(il,ik)+drdx(ik)*n(il))+T3*drdx(il)*n(ik)
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
  end subroutine fbem_bem_harela2d_sbie_ext_st

  recursive subroutine fbem_bem_harela2d_sbie_ext_adp(e,reverse,xi_s,x_i,p,qsp,ks,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                 !! Element
    logical                            :: reverse           !! Reverse orientation
    real(kind=real64)                  :: xi_s(1,2)         !! Subdivision of the parent element
    real(kind=real64)                  :: x_i(2)            !! Collocation point position vector
    type(fbem_bem_harela2d_parameters) :: p                 !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp               !! Quasi-singular integration parameters
    integer                            :: ks                !! Current level of subdivisions
    integer                            :: ns                !! Maximum level of subdivision
    complex(kind=real64)               :: h(e%n_pnodes,2,2) !! h integration kernels matrix
    complex(kind=real64)               :: g(e%n_snodes,2,2) !! g integration kernels matrix
    ! Local
    integer              :: gln_near              ! 1D Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer              :: gln                   ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical              :: subdivide             ! True if subdivision has to be performed
    real(kind=real64)    :: barxi(1)              ! Nearest element coordinate with respect to collocation point
    real(kind=real64)    :: barxip(1)             ! Nearest element subdivision local coordinate with respect to collocation point
    real(kind=real64)    :: rmin                  ! Minimum distance between collocation point and barxi on the element
    real(kind=real64)    :: barr                  ! Telles jacobian at barxi
    real(kind=real64)    :: cl                    ! Characteristic length
    real(kind=real64)    :: d                     ! Normalized distance between collocation point and element subdivision
    integer              :: method                ! Method used in nearest point algorithm
    real(kind=real64)    :: tmp_xi_s(1,2)         ! Subdivision
    real(kind=real64)    :: x_s(2,e%n_gnodes)     ! Coordinates of the element subdivision
    complex(kind=real64) :: h_tmp(e%n_pnodes,2,2) ! h integration kernels matrix (temporary)
    complex(kind=real64) :: g_tmp(e%n_snodes,2,2) ! g integration kernels matrix (temporary)
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
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,4,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_harela2d_sbie_ext_adp',ns,'maximum number of subdivisions reached')
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
      call fbem_bem_harela2d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_harela2d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_harela2d_sbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,p,gln,h_tmp,g_tmp)
      h=h+h_tmp
      g=g+g_tmp
    end if
  end subroutine fbem_bem_harela2d_sbie_ext_adp

  !! Calculation of h and g when the collocation point is on the integration element.
  subroutine fbem_bem_harela2d_sbie_int(ngp,type_g,type_f1,type_f2,delta_f,x_nodes,reverse,xi_i,p,h,g)
    implicit none
    ! I/O
    integer                            :: ngp                             !! Number of Gauss point to be used (<=32)
    integer                            :: type_g                          !! Geometrical interpolation
    integer                            :: type_f1                         !! Functional interpolation (primary variables)
    integer                            :: type_f2                         !! Geometrical interpolation (secondary variables)
    real(kind=real64)                  :: delta_f                         !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64)                  :: x_nodes(2,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
    logical                            :: reverse                         !! Reverse normal vector
    real(kind=real64)                  :: xi_i                            !! Reference coordinate of the singular point.
    type(fbem_bem_harela2d_parameters) :: p                               !! Parameters of the region
    complex(kind=real64)               :: h(fbem_n_nodes(type_f1),2,2)    !! h integration kernels matrix
    complex(kind=real64)               :: g(fbem_n_nodes(type_f2),2,2)    !! g integration kernels matrix
    ! Local
    integer                            :: nnodes_g                        ! Number of nodes of the element interpolation
    integer                            :: kphi                            ! Counter variable for shape functions loops
    integer                            :: kip                             ! Counter of integration points
    real(kind=real64)                  :: x_i(2)                          ! Real coordinates of collocation point
    real(kind=real64)                  :: phi_f1(fbem_n_nodes(type_f1))   ! Functional shape functions values at xi
    real(kind=real64)                  :: phi_f1_i(fbem_n_nodes(type_f1)) ! Functional shape functions values at xi_i
    real(kind=real64)                  :: phi_f2(fbem_n_nodes(type_f2))   ! Functional shape functions values at xi
    real(kind=real64)                  :: phi_g(fbem_n_nodes(type_g))     ! Geometrical shape functions values at xi
    real(kind=real64)                  :: dphidxi_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi
    integer                            :: nsub                            ! Number of subdivision of the element
    integer                            :: ksub                            ! Counter of subdivision
    real(kind=real64)                  :: gamma                           ! Coordinate gamma
    real(kind=real64)                  :: w                               ! Weights of each integration point
    real(kind=real64)                  :: xip                             ! Coordinate  xip of subdivided element [0,1]
    real(kind=real64)                  :: js                              ! Jacobian of the xi [xisub(1,:),xisub[2,:]] -> xip [0,1] transformation
    real(kind=real64)                  :: xip_i(2)                        ! Singular point in xip space
    real(kind=real64)                  :: xi                              ! Coordinate xi
    real(kind=real64)                  :: xisub(2,2)                      ! Coordinates of subdivisions in xi space
    real(kind=real64)                  :: aux(10)                         ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)                  :: x(2)                            ! Position vector at xi
    real(kind=real64)                  :: T(2)                            ! Tangent vector at xi
    real(kind=real64)                  :: N(2)                            ! Normal vector at xi
    real(kind=real64)                  :: rv(2)                           ! Radiovector between collocation point and integration point (x-x_i)
    real(kind=real64)                  :: r, r2, dr1, logr                ! Radiovector module, its inverse and log(1/r)
    real(kind=real64)                  :: ra, rb                          ! Radiovector from collocation point to element vertices
    real(kind=real64)                  :: drdx(2)                         ! Radiovector derivatives with respect to x_k
    real(kind=real64)                  :: jg                              ! Geometric jacobian
    real(kind=real64)                  :: drdn                            ! Partial derivative of r respect to unit normal
    real(kind=real64)                  :: drdt                            ! Partial derivative of r respect to unit tangent
    real(kind=real64)                  :: jw                              ! Jacobians * weight
    real(kind=real64)                  :: phif1jw(fbem_n_nodes(type_f1))  ! Auxiliary variables for integrand evaluation
    real(kind=real64)                  :: phif2jw(fbem_n_nodes(type_f2))  ! Auxiliary variables for integrand evaluation
    type(fbem_telles_parameters)       :: telles_parameters               ! Telles parameters
    real(kind=real64)                  :: jt                              ! Telles jacobian
    integer                            :: il, ik                          ! Counter variables for load direction and observation direction
    complex(kind=real64)               :: z(2)                            ! Wavenumbers
    complex(kind=real64)               :: KnR(0:2,2)                      ! Bessel functions decomposition
    complex(kind=real64)               :: psi, chi, T1, T2, T3            ! Fundamental solution components
    complex(kind=real64)               :: fs_u, fs_t                      ! Fundamental solutions values
    !
    ! Initialize
    !
    ! Initialize kernel matrices
    h=(0.d0,0.d0)
    g=(0.d0,0.d0)
    ! Number of nodes of the element
    nnodes_g=fbem_n_nodes(type_g)
    ! If xi_i belongs to one of the vertices, subdivision is not needed
    if (fbem_check_xi_vertex(xi_i).eqv.(.true.)) then
      nsub=1
    ! If xi_i is inside the element, 2 subdivisions are needed
    else
      nsub=2
    end if
    ! Setup the subdivisions transformation from xi to xip: [xisub1,xisub2] -> [0,1]
    select case (nsub)
      ! If just 1 subdivision
      case (1)
        ! Coodinates xi of the subdivision
        xisub(1,1)=-1.0d0
        xisub(2,1)= 1.0d0
        ! Coordinate xip of the collocation point
        ! If xi_i is at xi=-1
        if (xi_i.lt.0.0d0) xip_i(1)=0.0d0
        ! If xi_i is at xi=1
        if (xi_i.gt.0.0d0) xip_i(1)=1.0d0
      ! If 2 subdivisions
      case (2)
        ! Coordinates xi of the subdivision 1
        xisub(1,1)=-1.0d0
        xisub(2,1)=xi_i
        ! Coordinate xip of the collocation point
        xip_i(1)=1.0d0
        ! Coordinates xi of the subdivision 2
        xisub(1,2)=xi_i
        xisub(2,2)=1.0d0
        ! Coordinate xip of the collocation point
        xip_i(2)=0.0d0
    end select
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
    ! x_i
    x_i=0.0d0
    do kphi=1,nnodes_g
      x_i=x_i+phi_g(kphi)*x_nodes(:,kphi)
    end do
    ! Functional shape functions (primary variables) at xi_i
#   define etype type_f1
#   define delta delta_f
#   define xi xi_i
#   define phi phi_f1_i
#   include <phi_1d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
    !
    ! INTEGRATE WEAKLY SINGULAR INTEGRALS (log(r) integrals)
    !
    ! Loop through subdivisions
    do ksub=1,nsub
      ! Telles transformation parameters (barr=0 at the collocation point)
      telles_parameters=fbem_telles01_calculate_parameters(xip_i(ksub),0.0d0)
      ! Jacobian of xip->xi transformation (is constant)
      js=xisub(2,ksub)-xisub(1,ksub)
      ! Loop through integration points
      do kip=1,gl01_n(ngp)
        ! Coordinate and weight in gamma [0,1]
        gamma=gl01_xi(kip,ngp)
        w=gl01_w(kip,ngp)
        ! Calculate xip and gamma->xip jacobian
        call fbem_telles_xi_and_jacobian(telles_parameters,gamma,xip,jt)
        ! Calculate xi
        xi=js*xip+xisub(1,ksub)
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
        ! Radiovector
        rv=x-x_i
        ! Radiovector norm and its inverse
        r=dot_product(rv,rv)
        r=sqrt(r)
        logr=log(r)
        ! Jacobians * weight
        jw=jg*js*jt*w
        ! Functional shape functions (secondary variables) at xi
#       define etype type_f2
#       define delta delta_f
#       define phi phi_f2
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Functional shape functions * jacobians* weights
        phif2jw=phi_f2*jw
        ! Loop through load direction
        do il=1,2
          g(:,il,il)=g(:,il,il)+p%psi(1)*logr*phif2jw
        end do
      end do ! Loop through integration points
    end do ! Loop through subdivisions
    !
    ! INTEGRATE REGULAR INTEGRALS
    !
    ! Loop through subdivisions
    do ksub=1,nsub
      ! Jacobian of xip->xi transformation (is constant)
      js=xisub(2,ksub)-xisub(1,ksub)
      ! Loop through integration points
      do kip=1,gl01_n(ngp)
        ! Coordinate and weight in xip [0,1]
        xip=gl01_xi(kip,ngp)
        w=gl01_w(kip,ngp)
        ! Calculate xi
        xi=js*xip+xisub(1,ksub)
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
        ! Unit normal vector
        n=N/jg
        ! Unit tangent vector
        t=T/jg
        ! Radiovector
        rv=x-x_i
        ! Radiovector norm and its inverse
        r=dot_product(rv,rv)
        r=sqrt(r)
        r2=r**2
        dr1=1.0d0/r
        logr=log(r)
        ! r_{,k}
        drdx=rv*dr1
        ! dr/dn
        drdn=dot_product(drdx,n)
        ! dr/dGamma
        drdt=dot_product(drdx,t)
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
        ! COMPONENTS OF THE FUNDAMENTAL SOLUTION
        z(1)=c_im*p%k1*r
        z(2)=c_im*p%k2*r
        call fbem_BesselKnR_decomposed(2,z,KnR)
        psi=p%psi(2)+KnR(0,2)+dr1*p%psi(3)*KnR(1,1)+dr1*p%psi(4)*KnR(1,2)
        chi=p%chi(1)+(p%chi(2)*logr+p%chi(3))*r2+p%chi(4)*KnR(2,1)+KnR(2,2)
        T1=p%T1(2)*r+p%T1(3)*KnR(1,1)+p%T1(4)*KnR(1,2)+p%T1(5)*dr1*KnR(2,1)+p%T1(6)*dr1*KnR(2,2)
        T2=(p%T2(2)*logr+p%T2(3))*r+p%T2(4)*KnR(1,2)+p%T2(5)*dr1*KnR(2,1)+p%T2(6)*dr1*KnR(2,2)
        T3=(p%T3(2)*logr+p%T3(3))*r+p%T3(4)*KnR(1,1)+p%T3(5)*dr1*KnR(2,1)+p%T3(6)*dr1*KnR(2,2)
        ! Loop through load and observation direction
        do il=1,2
          do ik=1,2
            ! Fundamental solution values without cteu
            fs_u=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
            fs_t=T1*drdx(il)*drdx(ik)*drdn+T2*(drdn*c_dkr(il,ik)+drdx(ik)*n(il))+T3*drdx(il)*n(ik)&
                +p%T1(1)*dr1*drdx(il)*drdx(ik)*drdn+p%T2(1)*dr1*drdn*c_dkr(il,ik)
            ! Loop through nodes
            h(:,il,ik)=h(:,il,ik)+fs_t*phif1jw
            g(:,il,ik)=g(:,il,ik)+fs_u*phif2jw
          end do
        end do
        ! Regular and numerically evaluable part of the CPV integral
        h(:,1,2)=h(:,1,2)-p%T3(1)*dr1*drdt*(phi_f1-phi_f1_i)*jw
        h(:,2,1)=h(:,2,1)+p%T3(1)*dr1*drdt*(phi_f1-phi_f1_i)*jw
      end do ! Loop through integration points
    end do ! Loop through subdivisions
    ! Add CPV analytical terms to H12 and H21
    ! Calculate ra and rb
    ra=sqrt((x_nodes(1,1)-x_i(1))**2+(x_nodes(2,1)-x_i(2))**2)
    rb=sqrt((x_nodes(1,2)-x_i(1))**2+(x_nodes(2,2)-x_i(2))**2)
    ! The addition depends on the position of xi_i
    select case (nsub)
      ! If just 1 subdivision, xi_i=-1 or xi_i=1
      case (1)
        ! If xi_i is at xi=-1
        if (xi_i.lt.0.0d0) then
          h(:,1,2)=h(:,1,2)-p%T3(1)*phi_f1_i*log(rb)
          h(:,2,1)=h(:,2,1)+p%T3(1)*phi_f1_i*log(rb)
        end if
        ! If xi_i is at xi=1
        if (xi_i.gt.0.0d0) then
          h(:,1,2)=h(:,1,2)+p%T3(1)*phi_f1_i*log(ra)
          h(:,2,1)=h(:,2,1)-p%T3(1)*phi_f1_i*log(ra)
        end if
      ! If 2 subdivisions, xi_i is inside de element -1<xi_i<1
      case (2)
        h(:,1,2)=h(:,1,2)-p%T3(1)*phi_f1_i*(log(rb)-log(ra))
        h(:,2,1)=h(:,2,1)+p%T3(1)*phi_f1_i*(log(rb)-log(ra))
    end select
    ! Multiply h and g by constants of t* and u* respectively
    h=p%cte_t*h
    g=p%cte_u*g
    ! If the normal has to be reversed, then h=-h
    if (reverse) h=-h
  end subroutine fbem_bem_harela2d_sbie_int

  subroutine fbem_bem_harela2d_sbie_auto(mode,e,reverse,x_i,p,qsp,ns,h,g)
    implicit none
    ! I/O
    integer                            :: mode              !! 0: Solve, 1: Internal points
    type(fbem_bem_element)             :: e                 !! Integration element
    logical                            :: reverse           !! Reverse orientation
    real(kind=real64)                  :: x_i(2)            !! Collocation point
    type(fbem_bem_harela2d_parameters) :: p                 !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp               !! Quasi-singular integration parameters
    integer                            :: ns                !! Maximum level of subdivisions
    complex(kind=real64)               :: h(e%n_pnodes,2,2) !! h integration kernel
    complex(kind=real64)               :: g(e%n_snodes,2,2) !! g integration kernel
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
        if (mode.eq.1) then
          h=0.d0
          g=0.d0
          return
        end if
      else
        delta=0
      end if
    end if
    ! Integrate
    select case (delta)
      case (1)
        call fbem_bem_harela2d_sbie_int(30,e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi(1),p,h,g)
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
          call fbem_bem_harela2d_sbie_ext_pre(ps,e,reverse,x_i,p,h,g)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_harela2d_sbie_ext_adp(e,reverse,xi_s,x_i,p,qsp,1,ns,h,g)
        end if
    end select
  end subroutine fbem_bem_harela2d_sbie_auto

  ! --------------------------------------------------------------------------------------------------------------------------------

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! BODY LOADS

  subroutine fbem_bem_harela2d_sbie_bl_ext_pre(ps,e,x_i,p,g)
    implicit none
    ! I/O
    integer                            :: ps                !! Selected precalculated dataset
    type(fbem_bem_element)             :: e                 !! Element
    real(kind=real64)                  :: x_i(2)            !! Collocation point position vector
    type(fbem_bem_harela2d_parameters) :: p                 !! Parameters of the region
    complex(kind=real64)               :: g(e%n_snodes,2,2) !! g integration kernels matrix
    ! Local
    integer              :: kip                  ! Counter variable for integration points loop
    real(kind=real64)    :: x(2)                 ! Position vector at integration point
    real(kind=real64)    :: sphijw(e%n_snodes)   ! phi^s * jacobian * weight at integration point
    real(kind=real64)    :: rv(2)                ! Radiovector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, r2, d1r1, logr    ! Radiovector module, squared, inverses and log(r)
    real(kind=real64)    :: drdx(2)              ! Radiovector derivatives with respect to x_k
    integer              :: il, ik               ! Counter variables for load direction and observation direction
    complex(kind=real64) :: z(2)                 ! Bessel functions arguments: z=ikr
    complex(kind=real64) :: KnR(0:2,2)           ! Bessel functions decomposition
    complex(kind=real64) :: psi, chi             ! Fundamental solution components
    complex(kind=real64) :: fs_u                 ! Fundamental solutions values
    g=0
    do kip=1,e%ps_ngp(ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      x=e%ps_x(:,kip,ps)
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      r2=r**2
      d1r1=1.0d0/r
      logr=log(r)
      drdx=rv*d1r1
      z(1)=c_im*p%k1*r
      z(2)=c_im*p%k2*r
      call fbem_BesselKnR_decomposed(2,z,KnR)
      psi=p%psi(1)*logr+p%psi(2)+KnR(0,2)+d1r1*p%psi(3)*KnR(1,1)+d1r1*p%psi(4)*KnR(1,2)
      chi=p%chi(1)+(p%chi(2)*logr+p%chi(3))*r2+p%chi(4)*KnR(2,1)+KnR(2,2)
      do il=1,2
        do ik=1,2
          fs_u=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
          g(:,il,ik)=g(:,il,ik)+fs_u*sphijw
        end do
      end do
    end do
    g=p%cte_u*g
  end subroutine fbem_bem_harela2d_sbie_bl_ext_pre

  subroutine fbem_bem_harela2d_sbie_bl_ext_st(e,xi_s,x_i,barxip,barr,p,gln,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                 !! Integration element
    real(kind=real64)                  :: xi_s(1,2)         !! Coordinates of the subdivision of the parent element (xi space [-1,1])
    real(kind=real64)                  :: x_i(2)            !! Collocation point position vector
    real(kind=real64)                  :: barxip(1)         !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)                  :: barr              !! Telles jacobian at barxip
    type(fbem_bem_harela2d_parameters) :: p                 !! Parameters of the region
    integer                            :: gln               !! 1D Gauss-Legendre number of integration points (<=32)
    complex(kind=real64)               :: g(e%n_snodes,2,2) !! g kernel vector
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
    real(kind=real64)            :: sphi(e%n_snodes)     ! Functional shape functions values
    real(kind=real64)            :: x(2)                 ! Position vector at xi
    real(kind=real64)            :: T(2)                 ! Tangent vector at xi
    real(kind=real64)            :: N(2)                 ! Normal vector at xi
    real(kind=real64)            :: rv(2)                ! Radiovector between collocation point and integration point (x-x_i)
    real(kind=real64)            :: r, r2, d1r, logr     ! Radiovector module, squared, its inverses and log(r)
    real(kind=real64)            :: drdx(2)              ! Radiovector derivatives with respect to x_k
    real(kind=real64)            :: jg                   ! Geometric jacobian
    real(kind=real64)            :: jw                   ! Jacobian * weight
    real(kind=real64)            :: sphijw(e%n_snodes)   ! Auxiliary variables for integrand evaluation
    integer                      :: il, ik               ! Counter variables for load direction and observation direction
    complex(kind=real64)         :: z(2)                 ! Bessel functions arguments: z=ikr
    complex(kind=real64)         :: KnR(0:2,2)           ! Bessel functions decomposition
    complex(kind=real64)         :: psi, chi             ! Fundamental solution components
    complex(kind=real64)         :: fs_u                 ! Fundamental solutions values
    ! Initialization
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
      r2=r**2
      d1r=1.0d0/r
      logr=log(r)
      drdx=rv*d1r
      ! Jacobian * weight
      jw=jg*js*jt*w
      ! FUNCTIONAL SHAPE FUNCTIONS
      ! Functional shape functions (secondary variables) at xi
#     define etype e%stype
#     define delta e%stype_delta
#     define phi sphi
#     include <phi_1d.rc>
#     undef etype
#     undef delta
#     undef phi
      ! Functional shape functions * jacobians* weights
      sphijw=sphi*jw
      ! COMPONENTS OF THE FUNDAMENTAL SOLUTION
      z(1)=c_im*p%k1*r
      z(2)=c_im*p%k2*r
      call fbem_BesselKnR_decomposed(2,z,KnR)
      psi=p%psi(1)*logr+p%psi(2)+KnR(0,2)+d1r*p%psi(3)*KnR(1,1)+d1r*p%psi(4)*KnR(1,2)
      chi=p%chi(1)+(p%chi(2)*logr+p%chi(3))*r2+p%chi(4)*KnR(2,1)+KnR(2,2)
      ! Add
      do il=1,2
        do ik=1,2
          fs_u=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
          g(:,il,ik)=g(:,il,ik)+fs_u*sphijw
        end do
      end do
    end do
    g=p%cte_u*g
  end subroutine fbem_bem_harela2d_sbie_bl_ext_st

  recursive subroutine fbem_bem_harela2d_sbie_bl_ext_adp(e,xi_s,x_i,p,qsp,ks,ns,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                 !! Element
    real(kind=real64)                  :: xi_s(1,2)         !! Subdivision of the parent element
    real(kind=real64)                  :: x_i(2)            !! Collocation point position vector
    type(fbem_bem_harela2d_parameters) :: p                 !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp               !! Quasi-singular integration parameters
    integer                            :: ks                !! Current level of subdivisions
    integer                            :: ns                !! Maximum level of subdivision
    complex(kind=real64)               :: g(e%n_snodes,2,2) !! g integration kernels matrix
    ! Local
    integer              :: gln_near              ! 1D Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer              :: gln                   ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical              :: subdivide             ! True if subdivision has to be performed
    real(kind=real64)    :: barxi(1)              ! Nearest element coordinate with respect to collocation point
    real(kind=real64)    :: barxip(1)             ! Nearest element subdivision local coordinate with respect to collocation point
    real(kind=real64)    :: rmin                  ! Minimum distance between collocation point and barxi on the element
    real(kind=real64)    :: barr                  ! Telles jacobian at barxi
    real(kind=real64)    :: cl                    ! Characteristic length
    real(kind=real64)    :: d                     ! Normalized distance between collocation point and element subdivision
    integer              :: method                ! Method used in nearest point algorithm
    real(kind=real64)    :: tmp_xi_s(1,2)         ! Subdivision
    real(kind=real64)    :: x_s(2,e%n_gnodes)     ! Coordinates of the element subdivision
    complex(kind=real64) :: g_tmp(e%n_snodes,2,2) ! g integration kernels matrix (temporary)
    ! Initialize
    if (ks.eq.1) then
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
        call fbem_warning_message(error_unit,0,'fbem_bem_harela2d_sbie_ext_adp',ns,'maximum number of subdivisions reached')
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
      call fbem_bem_harela2d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,p,qsp,ks+1,ns,g)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_harela2d_sbie_bl_ext_adp(e,tmp_xi_s,x_i,p,qsp,ks+1,ns,g)
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_harela2d_sbie_bl_ext_st(e,xi_s,x_i,barxip,barr,p,gln,g_tmp)
      g=g+g_tmp
    end if
  end subroutine fbem_bem_harela2d_sbie_bl_ext_adp

  !! Calculation of g when the collocation point is on the integration element.
  subroutine fbem_bem_harela2d_sbie_bl_int(ngp,type_g,type_f1,type_f2,delta_f,x_nodes,xi_i,p,g)
    implicit none
    ! I/O
    integer                            :: ngp                             !! Number of Gauss point to be used (<=32)
    integer                            :: type_g                          !! Geometrical interpolation
    integer                            :: type_f1                         !! Functional interpolation (primary variables)
    integer                            :: type_f2                         !! Geometrical interpolation (secondary variables)
    real(kind=real64)                  :: delta_f                         !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64)                  :: x_nodes(2,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
    real(kind=real64)                  :: xi_i                            !! Reference coordinate of the singular point.
    type(fbem_bem_harela2d_parameters) :: p                               !! Parameters of the region
    complex(kind=real64)               :: g(fbem_n_nodes(type_f2),2,2)    !! g integration kernels matrix
    ! Local
    integer                            :: nnodes_g                        ! Number of nodes of the element interpolation
    integer                            :: kphi                            ! Counter variable for shape functions loops
    integer                            :: kip                             ! Counter of integration points
    real(kind=real64)                  :: x_i(2)                          ! Real coordinates of collocation point
    real(kind=real64)                  :: phi_f1(fbem_n_nodes(type_f1))   ! Functional shape functions values at xi
    real(kind=real64)                  :: phi_f1_i(fbem_n_nodes(type_f1)) ! Functional shape functions values at xi_i
    real(kind=real64)                  :: phi_f2(fbem_n_nodes(type_f2))   ! Functional shape functions values at xi
    real(kind=real64)                  :: phi_g(fbem_n_nodes(type_g))     ! Geometrical shape functions values at xi
    real(kind=real64)                  :: dphidxi_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi
    integer                            :: nsub                            ! Number of subdivision of the element
    integer                            :: ksub                            ! Counter of subdivision
    real(kind=real64)                  :: gamma                           ! Coordinate gamma
    real(kind=real64)                  :: w                               ! Weights of each integration point
    real(kind=real64)                  :: xip                             ! Coordinate  xip of subdivided element [0,1]
    real(kind=real64)                  :: js                              ! Jacobian of the xi [xisub(1,:),xisub[2,:]] -> xip [0,1] transformation
    real(kind=real64)                  :: xip_i(2)                        ! Singular point in xip space
    real(kind=real64)                  :: xi                              ! Coordinate xi
    real(kind=real64)                  :: xisub(2,2)                      ! Coordinates of subdivisions in xi space
    real(kind=real64)                  :: aux(10)                         ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)                  :: x(2)                            ! Position vector at xi
    real(kind=real64)                  :: T(2)                            ! Tangent vector at xi
    real(kind=real64)                  :: N(2)                            ! Normal vector at xi
    real(kind=real64)                  :: rv(2)                           ! Radiovector between collocation point and integration point (x-x_i)
    real(kind=real64)                  :: r, r2, dr1, logr                ! Radiovector module, its inverse and log(1/r)
    real(kind=real64)                  :: ra, rb                          ! Radiovector from collocation point to element vertices
    real(kind=real64)                  :: drdx(2)                         ! Radiovector derivatives with respect to x_k
    real(kind=real64)                  :: jg                              ! Geometric jacobian
    real(kind=real64)                  :: jw                              ! Jacobians * weight
    real(kind=real64)                  :: phif2jw(fbem_n_nodes(type_f2))  ! Auxiliary variables for integrand evaluation
    type(fbem_telles_parameters)       :: telles_parameters               ! Telles parameters
    real(kind=real64)                  :: jt                              ! Telles jacobian
    integer                            :: il, ik                          ! Counter variables for load direction and observation direction
    complex(kind=real64)               :: z(2)                            ! Wavenumbers
    complex(kind=real64)               :: KnR(0:2,2)                      ! Bessel functions decomposition
    complex(kind=real64)               :: psi, chi                        ! Fundamental solution components
    complex(kind=real64)               :: fs_u                            ! Fundamental solutions values
    !
    ! Initialize
    !
    ! Initialize kernel matrices
    g=(0.d0,0.d0)
    ! Number of nodes of the element
    nnodes_g=fbem_n_nodes(type_g)
    ! If xi_i belongs to one of the vertices, subdivision is not needed
    if (fbem_check_xi_vertex(xi_i).eqv.(.true.)) then
      nsub=1
    ! If xi_i is inside the element, 2 subdivisions are needed
    else
      nsub=2
    end if
    ! Setup the subdivisions transformation from xi to xip: [xisub1,xisub2] -> [0,1]
    select case (nsub)
      ! If just 1 subdivision
      case (1)
        ! Coodinates xi of the subdivision
        xisub(1,1)=-1.0d0
        xisub(2,1)= 1.0d0
        ! Coordinate xip of the collocation point
        ! If xi_i is at xi=-1
        if (xi_i.lt.0.0d0) xip_i(1)=0.0d0
        ! If xi_i is at xi=1
        if (xi_i.gt.0.0d0) xip_i(1)=1.0d0
      ! If 2 subdivisions
      case (2)
        ! Coordinates xi of the subdivision 1
        xisub(1,1)=-1.0d0
        xisub(2,1)=xi_i
        ! Coordinate xip of the collocation point
        xip_i(1)=1.0d0
        ! Coordinates xi of the subdivision 2
        xisub(1,2)=xi_i
        xisub(2,2)=1.0d0
        ! Coordinate xip of the collocation point
        xip_i(2)=0.0d0
    end select
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
    ! x_i
    x_i=0.0d0
    do kphi=1,nnodes_g
      x_i=x_i+phi_g(kphi)*x_nodes(:,kphi)
    end do
    ! Functional shape functions (primary variables) at xi_i
#   define etype type_f1
#   define delta delta_f
#   define xi xi_i
#   define phi phi_f1_i
#   include <phi_1d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
    !
    ! INTEGRATE WEAKLY SINGULAR INTEGRALS (log(r) integrals)
    !
    ! Loop through subdivisions
    do ksub=1,nsub
      ! Telles transformation parameters (barr=0 at the collocation point)
      telles_parameters=fbem_telles01_calculate_parameters(xip_i(ksub),0.0d0)
      ! Jacobian of xip->xi transformation (is constant)
      js=xisub(2,ksub)-xisub(1,ksub)
      ! Loop through integration points
      do kip=1,gl01_n(ngp)
        ! Coordinate and weight in gamma [0,1]
        gamma=gl01_xi(kip,ngp)
        w=gl01_w(kip,ngp)
        ! Calculate xip and gamma->xip jacobian
        call fbem_telles_xi_and_jacobian(telles_parameters,gamma,xip,jt)
        ! Calculate xi
        xi=js*xip+xisub(1,ksub)
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
        ! Radiovector
        rv=x-x_i
        ! Radiovector norm and its inverse
        r=dot_product(rv,rv)
        r=sqrt(r)
        logr=log(r)
        ! Jacobians * weight
        jw=jg*js*jt*w
        ! Functional shape functions (secondary variables) at xi
#       define etype type_f2
#       define delta delta_f
#       define phi phi_f2
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Functional shape functions * jacobians* weights
        phif2jw=phi_f2*jw
        ! Loop through load direction
        do il=1,2
          g(:,il,il)=g(:,il,il)+p%psi(1)*logr*phif2jw
        end do
      end do ! Loop through integration points
    end do ! Loop through subdivisions
    !
    ! INTEGRATE REGULAR INTEGRALS
    !
    ! Loop through subdivisions
    do ksub=1,nsub
      ! Jacobian of xip->xi transformation (is constant)
      js=xisub(2,ksub)-xisub(1,ksub)
      ! Loop through integration points
      do kip=1,gl01_n(ngp)
        ! Coordinate and weight in xip [0,1]
        xip=gl01_xi(kip,ngp)
        w=gl01_w(kip,ngp)
        ! Calculate xi
        xi=js*xip+xisub(1,ksub)
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
        ! Unit normal vector
        n=N/jg
        ! Unit tangent vector
        t=T/jg
        ! Radiovector
        rv=x-x_i
        ! Radiovector norm and its inverse
        r=dot_product(rv,rv)
        r=sqrt(r)
        r2=r**2
        dr1=1.0d0/r
        logr=log(r)
        ! r_{,k}
        drdx=rv*dr1
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
        phif2jw=phi_f2*jw
        ! COMPONENTS OF THE FUNDAMENTAL SOLUTION
        z(1)=c_im*p%k1*r
        z(2)=c_im*p%k2*r
        call fbem_BesselKnR_decomposed(2,z,KnR)
        psi=p%psi(2)+KnR(0,2)+dr1*p%psi(3)*KnR(1,1)+dr1*p%psi(4)*KnR(1,2)
        chi=p%chi(1)+(p%chi(2)*logr+p%chi(3))*r2+p%chi(4)*KnR(2,1)+KnR(2,2)
        ! Loop through load and observation direction
        do il=1,2
          do ik=1,2
            ! Fundamental solution values without cteu
            fs_u=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
            g(:,il,ik)=g(:,il,ik)+fs_u*phif2jw
          end do
        end do
      end do ! Loop through integration points
    end do ! Loop through subdivisions
    g=p%cte_u*g
  end subroutine fbem_bem_harela2d_sbie_bl_int

  subroutine fbem_bem_harela2d_sbie_bl_auto(mode,e,x_i,p,qsp,ns,g)
    implicit none
    ! I/O
    integer                            :: mode              !! 0: Solve, 1: Internal points
    type(fbem_bem_element)             :: e                 !! Integration element
    real(kind=real64)                  :: x_i(2)            !! Collocation point
    type(fbem_bem_harela2d_parameters) :: p                 !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp               !! Quasi-singular integration parameters
    integer                            :: ns                !! Maximum level of subdivisions
    complex(kind=real64)               :: g(e%n_snodes,2,2) !! g integration kernel
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
        call fbem_error_message(output_unit,0,'fbem_bem_harela2d_sbie_bl_auto',0,'it is not possible to collocate at a point load')
      else
        call fbem_bem_harela2d_sbie_u(e%x(:,1),x_i,p,g(1,:,:))
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
        call fbem_nearest_element_point_bem(2,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
        if (d.le.1.d-12) then
          delta=1
          if (mode.eq.1) then
            g=0.d0
            return
          end if
        else
          delta=0
        end if
      end if
      ! Integrate
      select case (delta)
        case (1)
          call fbem_bem_harela2d_sbie_bl_int(30,e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,barxi(1),p,g)
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
            call fbem_bem_harela2d_sbie_bl_ext_pre(ps,e,x_i,p,g)
          ! Integrate using an adaptative algorithm
          else
            call fbem_bem_harela2d_sbie_bl_ext_adp(e,xi_s,x_i,p,qsp,1,ns,g)
          end if
      end select
    end if
  end subroutine fbem_bem_harela2d_sbie_bl_auto

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (HBIE)

  !! Fundamental solution d*
  subroutine fbem_bem_harela2d_hbie_d(x,x_i,n_i,p,do)
    implicit none
    ! I/O
    real(kind=real64)                  :: x(2)    !! Observation point
    real(kind=real64)                  :: x_i(2)  !! Collocation point
    real(kind=real64)                  :: n_i(2)  !! Collocation point unit normal
    type(fbem_bem_harela2d_parameters) :: p       !! Parameters of the region
    complex(kind=real64)               :: do(2,2) !! d*_{lk}
    ! Local
    integer              :: il, ik
    real(kind=real64)    :: rv(2), r, r2, d1r1, d1r2, logr, drdx(2), drdni                          ! Partial derivative of r respect to unit normal at collocation point
    complex(kind=real64) :: z(2), KnR(0:2,2)
    complex(kind=real64) :: T1, T2, T3
    rv=x-x_i
    r=sqrt(dot_product(rv,rv))
    r2=r**2
    d1r1=1.0d0/r
    d1r2=d1r1*d1r1
    logr=log(r)
    drdx=rv*d1r1
    drdni=-dot_product(drdx,n_i)
    z(1)=c_im*p%k1*r
    z(2)=c_im*p%k2*r
    call fbem_BesselKnR_decomposed(2,z,KnR)
    T1=p%T1(1)*d1r1+p%T1(2)*r+p%T1(3)*KnR(1,1)+p%T1(4)*KnR(1,2)+p%T1(5)*d1r1*KnR(2,1)+p%T1(6)*d1r1*KnR(2,2)
    T2=p%T2(1)*d1r1+(p%T2(2)*logr+p%T2(3))*r+p%T2(4)*KnR(1,2)+p%T2(5)*d1r1*KnR(2,1)+p%T2(6)*d1r1*KnR(2,2)
    T3=p%T3(1)*d1r1+(p%T3(2)*logr+p%T3(3))*r+p%T3(4)*KnR(1,1)+p%T3(5)*d1r1*KnR(2,1)+p%T3(6)*d1r1*KnR(2,2)
    do il=1,2
      do ik=1,2
        do(il,ik)=T1*drdx(il)*drdx(ik)*drdni-T2*(-drdni*c_dkr(il,ik)+drdx(il)*n_i(ik))-T3*drdx(ik)*n_i(il)
      end do
    end do
    do=p%cte_d*do
  end subroutine fbem_bem_harela2d_hbie_d

  !! Fundamental solution s*
  subroutine fbem_bem_harela2d_hbie_s(x,n,x_i,n_i,p,so)
    implicit none
    ! I/O
    real(kind=real64)                  :: x(2)    !! Observation point
    real(kind=real64)                  :: n(2)    !! Observation point unit normal
    real(kind=real64)                  :: x_i(2)  !! Collocation point
    real(kind=real64)                  :: n_i(2)  !! Collocation point unit normal
    type(fbem_bem_harela2d_parameters) :: p       !! Parameters of the region
    complex(kind=real64)               :: so(2,2) !! s*_{lk}
    ! Local
    integer              :: il, ik
    real(kind=real64)    :: rv(2), r, r2, d1r1, d1r2, logr, drdx(2), drdn, drdni, n_dot_ni
    complex(kind=real64) :: z(2), KnR(0:2,2)
    complex(kind=real64) :: S1, S2, S3, S4, S5
    rv=x-x_i
    r=sqrt(dot_product(rv,rv))
    r2=r**2
    d1r1=1.0d0/r
    d1r2=d1r1*d1r1
    logr=log(r)
    drdx=rv*d1r1
    drdn=dot_product(drdx,n)
    drdni=-dot_product(drdx,n_i)
    n_dot_ni=dot_product(n,n_i)
    z(1)=c_im*p%k1*r
    z(2)=c_im*p%k2*r
    call fbem_BesselKnR_decomposed(2,z,KnR)
    S1=p%S1(1)*d1r2+p%S1(2)+p%S1(3)*KnR(0,2)+p%S1(4)*d1r1*KnR(1,1)+p%S1(5)*d1r1*KnR(1,2)+p%S1(6)*d1r2*KnR(2,1)&
      +p%S1(7)*d1r2*KnR(2,2)
    S2=p%S2(1)*d1r2+p%S2(2)+(p%S2(3)*logr+p%S2(4))*r2+p%S2(5)*d1r1*KnR(1,1)+p%S2(6)*d1r1*KnR(1,2)+p%S2(7)*d1r2*KnR(2,1)&
      +p%S2(8)*KnR(2,1)+p%S2(9)*d1r2*KnR(2,2)
    S3=p%S3(1)*d1r2+p%S3(2)+p%S3(3)*KnR(0,1)+p%S3(4)*KnR(0,2)+p%S3(5)*d1r1*KnR(1,1)+p%S3(6)*d1r1*KnR(1,2)&
      +p%S3(7)*d1r2*KnR(2,1)+p%S3(8)*d1r2*KnR(2,2)
    S4=p%S4(1)*d1r2+p%S4(2)*logr+p%S4(3)+p%S4(4)*d1r1*KnR(1,2)+p%S4(5)*d1r2*KnR(2,1)+p%S4(6)*d1r2*KnR(2,2)
    S5=p%S5(1)*d1r2+p%S5(2)*logr+p%S5(3)+p%S5(4)*KnR(0,1)+p%S5(5)*d1r1*KnR(1,1)+p%S5(6)*d1r2*KnR(2,1)+p%S5(7)*d1r2*KnR(2,2)
    do il=1,2
      do ik=1,2
        so(il,ik)=S1*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni-c_dkr(il,ik)*drdn*drdni+drdx(il)*drdx(ik)*n_dot_ni)&
                 +S2*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni)+S3*drdx(il)*drdx(ik)*drdn*drdni&
                 +S4*(c_dkr(il,ik)*n_dot_ni+n_i(ik)*n(il))+S5*n(ik)*n_i(il)
      end do
    end do
    so=p%cte_s*so
  end subroutine fbem_bem_harela2d_hbie_s

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! BOUNDARY ELEMENTS

  !! Calculation of m and l when the collocation point is outside the integration element, using a precalculated dataset of data.
  subroutine fbem_bem_harela2d_hbie_ext_pre(ps,e,reverse,x_i,n_i,p,m,l)
    implicit none
    ! I/O
    integer                            :: ps                !! Selected precalculated dataset
    type(fbem_bem_element)             :: e                 !! Element
    logical                            :: reverse           !! Reverse normal vector
    real(kind=real64)                  :: x_i(2)            !! Collocation point position vector
    real(kind=real64)                  :: n_i(2)            !! Collocation point unit normal vector
    type(fbem_bem_harela2d_parameters) :: p                 !! Parameters of the region
    complex(kind=real64)               :: m(e%n_pnodes,2,2) !! m kernels matrix
    complex(kind=real64)               :: l(e%n_snodes,2,2) !! l kernels matrix
    ! Local
    integer              :: kip                            ! Counter variable for integration points loop
    real(kind=real64)    :: x(2)                           ! Position vector at integration point
    real(kind=real64)    :: n(2)                           ! Unit normal vector at integration point
    real(kind=real64)    :: pphijw(e%n_pnodes)             ! phi^p * jacobian * weight at integration point
    real(kind=real64)    :: sphijw(e%n_snodes)             ! phi^s * jacobian * weight at integration point
    real(kind=real64)    :: rv(2)                          ! Radiovector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, r2, d1r1, d1r2, logr        ! Radiovector module, squared, inverses and log(r)
    real(kind=real64)    :: drdx(2)                        ! Radiovector derivatives with respect to x_k
    real(kind=real64)    :: drdn                           ! Partial derivative of r respect to unit normal
    real(kind=real64)    :: drdni                          ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)    :: n_dot_ni                       ! Dot product of unit normals
    integer              :: il, ik                         ! Counter variables for load direction and observation direction
    complex(kind=real64) :: z(2)                           ! Wavenumbers
    complex(kind=real64) :: KnR(0:2,2)                     ! Bessel functions decomposition
    complex(kind=real64) :: T1, T2, T3, S1, S2, S3, S4, S5 ! Fundamental solution components
    complex(kind=real64) :: fs_s, fs_d                     ! Fundamental solutions values
    ! Initialize
    m=(0.0d0,0.0d0)
    l=(0.0d0,0.0d0)
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      x=e%ps_x(:,kip,ps)
      n=e%ps_n(:,kip,ps)
      pphijw=e%ps_pphijw(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      r2=r**2
      d1r1=1.0d0/r
      d1r2=d1r1*d1r1
      logr=log(r)
      drdx=rv*d1r1
      drdn=dot_product(drdx,n)
      drdni=-dot_product(drdx,n_i)
      n_dot_ni=dot_product(n,n_i)
      z(1)=c_im*p%k1*r
      z(2)=c_im*p%k2*r
      call fbem_BesselKnR_decomposed(2,z,KnR)
      T1=p%T1(1)*d1r1+p%T1(2)*r+p%T1(3)*KnR(1,1)+p%T1(4)*KnR(1,2)+p%T1(5)*d1r1*KnR(2,1)+p%T1(6)*d1r1*KnR(2,2)
      T2=p%T2(1)*d1r1+(p%T2(2)*logr+p%T2(3))*r+p%T2(4)*KnR(1,2)+p%T2(5)*d1r1*KnR(2,1)+p%T2(6)*d1r1*KnR(2,2)
      T3=p%T3(1)*d1r1+(p%T3(2)*logr+p%T3(3))*r+p%T3(4)*KnR(1,1)+p%T3(5)*d1r1*KnR(2,1)+p%T3(6)*d1r1*KnR(2,2)
      S1=p%S1(1)*d1r2+p%S1(2)+p%S1(3)*KnR(0,2)+p%S1(4)*d1r1*KnR(1,1)+p%S1(5)*d1r1*KnR(1,2)+p%S1(6)*d1r2*KnR(2,1)&
        +p%S1(7)*d1r2*KnR(2,2)
      S2=p%S2(1)*d1r2+p%S2(2)+(p%S2(3)*logr+p%S2(4))*r2+p%S2(5)*d1r1*KnR(1,1)+p%S2(6)*d1r1*KnR(1,2)+p%S2(7)*d1r2*KnR(2,1)&
        +p%S2(8)*KnR(2,1)+p%S2(9)*d1r2*KnR(2,2)
      S3=p%S3(1)*d1r2+p%S3(2)+p%S3(3)*KnR(0,1)+p%S3(4)*KnR(0,2)+p%S3(5)*d1r1*KnR(1,1)+p%S3(6)*d1r1*KnR(1,2)&
        +p%S3(7)*d1r2*KnR(2,1)+p%S3(8)*d1r2*KnR(2,2)
      S4=p%S4(1)*d1r2+p%S4(2)*logr+p%S4(3)+p%S4(4)*d1r1*KnR(1,2)+p%S4(5)*d1r2*KnR(2,1)+p%S4(6)*d1r2*KnR(2,2)
      S5=p%S5(1)*d1r2+p%S5(2)*logr+p%S5(3)+p%S5(4)*KnR(0,1)+p%S5(5)*d1r1*KnR(1,1)+p%S5(6)*d1r2*KnR(2,1)+p%S5(7)*d1r2*KnR(2,2)
      do il=1,2
        do ik=1,2
          fs_d=T1*drdx(il)*drdx(ik)*drdni-T2*(-drdni*c_dkr(il,ik)+drdx(il)*n_i(ik))-T3*drdx(ik)*n_i(il)
          fs_s=S1*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni-c_dkr(il,ik)*drdn*drdni+drdx(il)*drdx(ik)*n_dot_ni)&
              +S2*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni)+S3*drdx(il)*drdx(ik)*drdn*drdni&
              +S4*(c_dkr(il,ik)*n_dot_ni+n_i(ik)*n(il))+S5*n(ik)*n_i(il)
          m(:,il,ik)=m(:,il,ik)+fs_s*pphijw
          l(:,il,ik)=l(:,il,ik)+fs_d*sphijw
        end do
      end do
    end do
    ! Multiply y constants
    m=p%cte_s*m
    l=p%cte_d*l
    ! Reverse if needed
    if (reverse) m=-m
  end subroutine fbem_bem_harela2d_hbie_ext_pre

  subroutine fbem_bem_harela2d_hbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,p,gln,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                 !! Integration element
    logical                            :: reverse           !! Reverse normal vector
    real(kind=real64)                  :: xi_s(1,2)         !! Coordinates of the subdivision of the parent element (xi space [-1,1])
    real(kind=real64)                  :: x_i(2)            !! Collocation point position vector
    real(kind=real64)                  :: n_i(2)            !! Collocation point unit normal vector
    real(kind=real64)                  :: barxip(1)         !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)                  :: barr              !! Telles jacobian at barxip
    type(fbem_bem_harela2d_parameters) :: p                 !! Parameters of the region
    integer                            :: gln               !! 1D Gauss-Legendre number of integration points (<=32)
    complex(kind=real64)               :: m(e%n_pnodes,2,2) !! m kernels matrix
    complex(kind=real64)               :: l(e%n_snodes,2,2) !! l kernels matrix
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
    real(kind=real64)            :: r, r2, d1r1, d1r2, logr ! Radiovector module, squared, inverses and log(r)
    real(kind=real64)            :: drdx(2)                 ! Radiovector derivatives with respect to x_k
    real(kind=real64)            :: jg                      ! Geometric jacobian
    real(kind=real64)            :: drdn                    ! Partial derivative of r respect to unit normal
    real(kind=real64)            :: drdni                   ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)            :: n_dot_ni                ! Dot product of unit normals
    real(kind=real64)            :: jw                      ! Jacobians * weights
    real(kind=real64)            :: pphijw(e%n_pnodes)      ! Auxiliary variables for integrand evaluation
    real(kind=real64)            :: sphijw(e%n_snodes)      ! Auxiliary variables for integrand evaluation
    integer                      :: il, ik                  ! Counter variables for load direction and observation direction
    complex(kind=real64)         :: z(2)                    ! Bessel functions arguments: z=ikr
    complex(kind=real64)         :: KnR(0:2,2)              ! Bessel functions decomposition
    complex(kind=real64)         :: T1, T2, T3              ! Fundamental solution components
    complex(kind=real64)         :: S1, S2, S3, S4, S5      ! Fundamental solution components
    complex(kind=real64)         :: fs_s, fs_d              ! Fundamental solutions values
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
      r2=r**2
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
      z(1)=c_im*p%k1*r
      z(2)=c_im*p%k2*r
      call fbem_BesselKnR_decomposed(2,z,KnR)
      T1=p%T1(1)*d1r1+p%T1(2)*r+p%T1(3)*KnR(1,1)+p%T1(4)*KnR(1,2)+p%T1(5)*d1r1*KnR(2,1)+p%T1(6)*d1r1*KnR(2,2)
      T2=p%T2(1)*d1r1+(p%T2(2)*logr+p%T2(3))*r+p%T2(4)*KnR(1,2)+p%T2(5)*d1r1*KnR(2,1)+p%T2(6)*d1r1*KnR(2,2)
      T3=p%T3(1)*d1r1+(p%T3(2)*logr+p%T3(3))*r+p%T3(4)*KnR(1,1)+p%T3(5)*d1r1*KnR(2,1)+p%T3(6)*d1r1*KnR(2,2)
      S1=p%S1(1)*d1r2+p%S1(2)+p%S1(3)*KnR(0,2)+p%S1(4)*d1r1*KnR(1,1)+p%S1(5)*d1r1*KnR(1,2)+p%S1(6)*d1r2*KnR(2,1)&
        +p%S1(7)*d1r2*KnR(2,2)
      S2=p%S2(1)*d1r2+p%S2(2)+(p%S2(3)*logr+p%S2(4))*r2+p%S2(5)*d1r1*KnR(1,1)+p%S2(6)*d1r1*KnR(1,2)+p%S2(7)*d1r2*KnR(2,1)&
        +p%S2(8)*KnR(2,1)+p%S2(9)*d1r2*KnR(2,2)
      S3=p%S3(1)*d1r2+p%S3(2)+p%S3(3)*KnR(0,1)+p%S3(4)*KnR(0,2)+p%S3(5)*d1r1*KnR(1,1)+p%S3(6)*d1r1*KnR(1,2)&
        +p%S3(7)*d1r2*KnR(2,1)+p%S3(8)*d1r2*KnR(2,2)
      S4=p%S4(1)*d1r2+p%S4(2)*logr+p%S4(3)+p%S4(4)*d1r1*KnR(1,2)+p%S4(5)*d1r2*KnR(2,1)+p%S4(6)*d1r2*KnR(2,2)
      S5=p%S5(1)*d1r2+p%S5(2)*logr+p%S5(3)+p%S5(4)*KnR(0,1)+p%S5(5)*d1r1*KnR(1,1)+p%S5(6)*d1r2*KnR(2,1)+p%S5(7)*d1r2*KnR(2,2)
      ! Add
      do il=1,2
        do ik=1,2
          fs_d=T1*drdx(il)*drdx(ik)*drdni-T2*(drdx(il)*n_i(ik)-c_dkr(il,ik)*drdni)-T3*drdx(ik)*n_i(il)
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
  end subroutine fbem_bem_harela2d_hbie_ext_st

  recursive subroutine fbem_bem_harela2d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,p,qsp,ks,ns,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                 !! Element
    logical                            :: reverse           !! Reverse orientation
    real(kind=real64)                  :: xi_s(1,2)         !! Subdivision of the parent element
    real(kind=real64)                  :: x_i(2)            !! Collocation point position vector
    real(kind=real64)                  :: n_i(2)            !! Collocation point unit normal vector
    type(fbem_bem_harela2d_parameters) :: p                 !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp               !! Quasi-singular integration parameters
    integer                            :: ks                !! Current level of subdivisions
    integer                            :: ns                !! Maximum level of subdivision
    complex(kind=real64)               :: m(e%n_pnodes,2,2) !! m integration kernels matrix
    complex(kind=real64)               :: l(e%n_snodes,2,2) !! l integration kernels matrix
    ! Local
    integer              :: gln_near              ! 1D Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer              :: gln                   ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical              :: subdivide             ! True if subdivision has to be performed
    real(kind=real64)    :: barxi(1)              ! Nearest element coordinate with respect to collocation point
    real(kind=real64)    :: barxip(1)             ! Nearest element subdivision local coordinate with respect to collocation point
    real(kind=real64)    :: rmin                  ! Minimum distance between collocation point and barxi on the element
    real(kind=real64)    :: barr                  ! Telles jacobian at barxi
    real(kind=real64)    :: cl                    ! Characteristic length
    real(kind=real64)    :: d                     ! Normalized distance between collocation point and element subdivision
    integer              :: method                ! Method used in nearest point algorithm
    real(kind=real64)    :: tmp_xi_s(1,2)         ! Subdivision
    real(kind=real64)    :: x_s(2,e%n_gnodes)     ! Coordinates of the element subdivision
    complex(kind=real64) :: m_tmp(e%n_pnodes,2,2) ! m kernels matrix (temporary)
    complex(kind=real64) :: l_tmp(e%n_snodes,2,2) ! l kernels matrix (temporary)
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
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,6,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_harela2d_hbie_ext_adp',ns,'maximum number of subdivisions reached')
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
      call fbem_bem_harela2d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_harela2d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_harela2d_hbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,p,gln,m_tmp,l_tmp)
      m=m+m_tmp
      l=l+l_tmp
    end if
  end subroutine fbem_bem_harela2d_hbie_ext_adp

  !! Calculation of m and l when the collocation point is on the integration element.
  subroutine fbem_bem_harela2d_hbie_int(ngp,type_g,type_f1,type_f2,delta_f,x_nodes,reverse,xi_i,p,m,l)
    implicit none
    ! I/O
    integer                            :: ngp                                 !! Number of Gauss point to be used (<=32)
    integer                            :: type_g                              !! Geometrical interpolation
    integer                            :: type_f1                             !! Functional interpolation (primary variables)
    integer                            :: type_f2                             !! Geometrical interpolation (secondary variables)
    real(kind=real64)                  :: delta_f                             !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64)                  :: x_nodes(2,fbem_n_nodes(type_g))     !! Position vectors of geometrical nodes
    logical                            :: reverse                             !! Reverse normal vector
    real(kind=real64)                  :: xi_i                                !! Reference coordinate of the singular point.
    type(fbem_bem_harela2d_parameters) :: p                               !! Parameters of the region
    complex(kind=real64)               :: m(fbem_n_nodes(type_f1),2,2)        !! m integration kernels matrix
    complex(kind=real64)               :: l(fbem_n_nodes(type_f2),2,2)        !! l integration kernels matrix
    ! Local
    integer                            :: nnodes_g                            ! Number of nodes of the element interpolation
    integer                            :: nnodes_f1                           ! Number of nodes of the element interpolation
    integer                            :: nnodes_f2                           ! Number of nodes of the element interpolation
    integer                            :: kphi                                ! Counter variable for shape functions loops
    integer                            :: kip                                 ! Counter of integration points
    real(kind=real64)                  :: x_i(2)                              ! Real coordinates of collocation point
    real(kind=real64)                  :: t_i(2)                              ! Unit tangent at collocation point
    real(kind=real64)                  :: n_i(2)                              ! Unit normal at collocation point
    real(kind=real64)                  :: phi_f1(fbem_n_nodes(type_f1))       ! Functional shape functions values at xi
    real(kind=real64)                  :: phi_f1_i(fbem_n_nodes(type_f1))     ! Functional shape functions values at xi_i
    real(kind=real64)                  :: dphidxi_f1_i(fbem_n_nodes(type_f1)) ! Functional shape functions derivatives values at xi_i
    real(kind=real64)                  :: phi_f2(fbem_n_nodes(type_f2))       ! Functional shape functions values at xi
    real(kind=real64)                  :: phi_f2_i(fbem_n_nodes(type_f2))     ! Functional shape functions values at xi_i
    real(kind=real64)                  :: phi_g(fbem_n_nodes(type_g))         ! Geometrical shape functions values at xi
    real(kind=real64)                  :: dphidxi_g(fbem_n_nodes(type_g))     ! Geometrical shape functions derivatives values at xi
    integer                            :: nsub                                ! Number of subdivision of the element
    integer                            :: ksub                                ! Counter of subdivision
    real(kind=real64)                  :: gamma                               ! Coordinate gamma
    real(kind=real64)                  :: w                                   ! Weights of each integration point
    real(kind=real64)                  :: xip                                 ! Coordinate  xip of subdivided element [0,1]
    real(kind=real64)                  :: js                                  ! Jacobian of the xi [xisub(1,:),xisub[2,:]] -> xip [0,1] transformation
    real(kind=real64)                  :: xip_i(2)                            ! Singular point in xip space
    real(kind=real64)                  :: xi                                  ! Coordinate xi
    real(kind=real64)                  :: xisub(2,2)                          ! Coordinates of subdivisions in xi space
    real(kind=real64)                  :: aux(10)                             ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)                  :: x(2)                                ! Position vector at xi
    real(kind=real64)                  :: T(2)                                ! Tangent vector at xi
    real(kind=real64)                  :: N(2)                                ! Normal vector at xi
    real(kind=real64)                  :: rv(2)                               ! Radiovector between collocation point and integration point (x-x_i)
    real(kind=real64)                  :: r, r2, dr1, dr2, logr               ! Radiovector module, squared, inverses and log(r)
    real(kind=real64)                  :: ra, rb                              ! Radiovector from collocation point to element vertices
    real(kind=real64)                  :: drdx(2)                             ! Radiovector derivatives with respect to x_k
    real(kind=real64)                  :: jg                                  ! Geometric jacobian
    real(kind=real64)                  :: jgi                                 ! Geometric jacobian at the collocation point
    real(kind=real64)                  :: drdn                                ! Partial derivative of r respect to unit normal
    real(kind=real64)                  :: drdni                               ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)                  :: n_dot_ni                            ! Dot product of unit normals
    real(kind=real64)                  :: drdt                                ! Partial derivative of r respect to unit tangent
    real(kind=real64)                  :: jw                                  ! Jacobians*weight
    real(kind=real64)                  :: phif1jw(fbem_n_nodes(type_f1))      ! Auxiliary variables for integrand evaluation
    real(kind=real64)                  :: phif2jw(fbem_n_nodes(type_f2))      ! Auxiliary variables for integrand evaluation
    type(fbem_telles_parameters)       :: telles_parameters                   ! Telles parameters
    real(kind=real64)                  :: jt                                  ! Telles jacobian
    integer                            :: il, ik                              ! Counter variables for load direction and observation direction
    complex(kind=real64)               :: z(2)                                ! Wavenumbers
    complex(kind=real64)               :: KnR(0:2,2)                          ! Bessel functions decomposition
    complex(kind=real64)               :: T1, T2, T3, S1, S2, S3, S4, S5      ! Fundamental solution components
    complex(kind=real64)               :: fs_s, fs_d                          ! Fundamental solutions values
    !
    ! Initialize
    !
    ! Initialize kernel matrices
    m=(0.0d0,0.0d0)
    l=(0.0d0,0.0d0)
    ! Number of nodes of the element
    nnodes_g=fbem_n_nodes(type_g)
    nnodes_f1=fbem_n_nodes(type_f1)
    nnodes_f2=fbem_n_nodes(type_f2)
    ! If xi_i belongs to one of the vertices, subdivision is not needed
    if (fbem_check_xi_vertex(xi_i).eqv.(.true.)) then
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                              'the HBIE cannot be collocated at a vertex')
    ! If xi_i is inside the element, 2 subdivisions are needed
    else
      nsub=2
      ! Coordinates xi of the subdivision 1
      xisub(1,1)=-1.0d0
      xisub(2,1)=xi_i
      ! Coordinate xip of the collocation point
      xip_i(1)=1.0d0
      ! Coordinates xi of the subdivision 2
      xisub(1,2)=xi_i
      xisub(2,2)=1.0d0
      ! Coordinate xip of the collocation point
      xip_i(2)=0.0d0
    end if
    !
    ! Calculate real coordinates, geometric jacobian, unit tangent and unit normal at collocation point
    !
    ! Geometrical shape functions and shape functions derivatives at xi_i
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
    ! x_i and T_i
    x_i=0.0d0
    T_i=0.0d0
    do kphi=1,nnodes_g
      x_i=x_i+phi_g(kphi)*x_nodes(:,kphi)
      T_i=T_i+dphidxi_g(kphi)*x_nodes(:,kphi)
    end do
    ! Normal vector
    N_i(1)=T_i(2)
    N_i(2)=-T_i(1)
    ! Geometric jacobian
    jgi=sqrt(T_i(1)**2+T_i(2)**2)
    ! Unit normal vector at collocation point
    n_i=N_i/jgi
    ! Unit tangent vector at collocation point
    t_i=T_i/jgi
    !
    ! Functional shape functions and its first derivatives (primary variable) at xi_i
    !
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
    !
    ! Functional shape functions (secondary variable) at xi_i
    !
#   define etype type_f2
#   define delta delta_f
#   define xi xi_i
#   define phi phi_f2_i
#   include <phi_1d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
    !
    ! INTEGRATE WEAKLY SINGULAR INTEGRALS (log(r) integrals)
    !
    ! Loop through subdivisions
    do ksub=1,nsub
      ! Telles transformation parameters (barr=0 at the collocation point)
      telles_parameters=fbem_telles01_calculate_parameters(xip_i(ksub),0.0d0)
      ! Jacobian of xip->xi transformation (is constant)
      js=xisub(2,ksub)-xisub(1,ksub)
      ! Loop through integration points
      do kip=1,gl01_n(ngp)
        ! Coordinate and weight in gamma [0,1]
        gamma=gl01_xi(kip,ngp)
        w=gl01_w(kip,ngp)
        ! Calculate xip and gamma->xip jacobian
        call fbem_telles_xi_and_jacobian(telles_parameters,gamma,xip,jt)
        ! Calculate xi
        xi=js*xip+xisub(1,ksub)
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
        ! Unit normal vector
        n=N/jg
        ! Radiovector
        rv=x-x_i
        ! Radiovector norm and its inverse
        r=dot_product(rv,rv)
        r=sqrt(r)
        logr=log(r)
        ! n · n_i
        n_dot_ni=dot_product(n,n_i)
        ! Functional shape functions (primary variable) at xi
#       define etype type_f1
#       define delta delta_f
#       define phi phi_f1
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Jacobian * weight
        jw=jg*js*jt*w
        ! Functional shape functions * jacobians* weight
        phif1jw=phi_f1*jw
        ! Loop through load direction
        do il=1,2
          ! Loop through observation direction
          do ik=1,2
            fs_s=p%S4(2)*logr*(c_dkr(il,ik)*n_dot_ni+n_i(ik)*n(il))+p%S5(2)*logr*n(ik)*n_i(il)
            ! Add to kernels
            m(:,il,ik)=m(:,il,ik)+fs_s*phif1jw
          end do
        end do
      end do ! Loop through integration points
    end do ! Loop through subdivisions
    !
    ! INTEGRATE REGULAR INTEGRALS
    !
    ! Loop through subdivisions
    do ksub=1,nsub
      ! Jacobian of xip->xi transformation (is constant)
      js=xisub(2,ksub)-xisub(1,ksub)
      ! Loop through integration points
      do kip=1,gl01_n(ngp)
        ! Coordinate and weight in xip [0,1]
        xip=gl01_xi(kip,ngp)
        w=gl01_w(kip,ngp)
        ! Calculate xi
        xi=js*xip+xisub(1,ksub)
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
        ! Unit normal vector
        n=N/jg
        ! Unit tangent vector
        t=T/jg
        ! Radiovector
        rv=x-x_i
        ! Radiovector norm and its inverse
        r=dot_product(rv,rv)
        r=sqrt(r)
        r2=r**2
        dr1=1.d0/r
        dr2=dr1**2
        logr=log(r)
        ! r_{,k}
        drdx=rv*dr1
        ! dr/dn
        drdn=dot_product(drdx,n)
        ! dr/dn_i
        drdni=-dot_product(drdx,n_i)
        ! n · n_i
        n_dot_ni=dot_product(n,n_i)
        ! dr/dGamma
        drdt=dot_product(drdx,t)
        ! Functional shape functions (primary variable) at xi
#       define etype type_f1
#       define delta delta_f
#       define phi phi_f1
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Functional shape functions (secondary variable) at xi
#       define etype type_f2
#       define delta delta_f
#       define phi phi_f2
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Jacobian * weight
        jw=jg*js*w
        ! Functional shape functions * jacobians* weight
        phif1jw=phi_f1*jw
        phif2jw=phi_f2*jw
        ! COMPONENTS OF THE FUNDAMENTAL SOLUTION
        z(1)=c_im*p%k1*r
        z(2)=c_im*p%k2*r
        call fbem_BesselKnR_decomposed(2,z,KnR)
        T1=p%T1(2)*r+p%T1(3)*KnR(1,1)+p%T1(4)*KnR(1,2)+p%T1(5)*dr1*KnR(2,1)+p%T1(6)*dr1*KnR(2,2)
        T2=(p%T2(2)*logr+p%T2(3))*r+p%T2(4)*KnR(1,2)+p%T2(5)*dr1*KnR(2,1)+p%T2(6)*dr1*KnR(2,2)
        T3=(p%T3(2)*logr+p%T3(3))*r+p%T3(4)*KnR(1,1)+p%T3(5)*dr1*KnR(2,1)+p%T3(6)*dr1*KnR(2,2)
        S1=p%S1(2)+p%S1(3)*KnR(0,2)+p%S1(4)*dr1*KnR(1,1)+p%S1(5)*dr1*KnR(1,2)+p%S1(6)*dr2*KnR(2,1)&
          +p%S1(7)*dr2*KnR(2,2)
        S2=p%S2(2)+(p%S2(3)*logr+p%S2(4))*r2+p%S2(5)*dr1*KnR(1,1)+p%S2(6)*dr1*KnR(1,2)+p%S2(7)*dr2*KnR(2,1)&
          +p%S2(8)*KnR(2,1)+p%S2(9)*dr2*KnR(2,2)
        S3=p%S3(2)+p%S3(3)*KnR(0,1)+p%S3(4)*KnR(0,2)+p%S3(5)*dr1*KnR(1,1)+p%S3(6)*dr1*KnR(1,2)&
          +p%S3(7)*dr2*KnR(2,1)+p%S3(8)*dr2*KnR(2,2)
        S4=p%S4(3)+p%S4(4)*dr1*KnR(1,2)+p%S4(5)*dr2*KnR(2,1)+p%S4(6)*dr2*KnR(2,2)
        S5=p%S5(3)+p%S5(4)*KnR(0,1)+p%S5(5)*dr1*KnR(1,1)+p%S5(6)*dr2*KnR(2,1)+p%S5(7)*dr2*KnR(2,2)
        ! Loop through load direction
        do il=1,2
          ! Loop through observation direction
          do ik=1,2
            ! Fundamental solutions values without constants
            fs_d=T1*drdx(il)*drdx(ik)*drdni-T2*(drdx(il)*n_i(ik)-c_dkr(il,ik)*drdni)-T3*drdx(ik)*n_i(il)&
                +p%T1(1)*dr1*drdx(il)*drdx(ik)*drdni+p%T2(1)*dr1*c_dkr(il,ik)*drdni&
                -p%T2(1)*dr1*drdx(il)*(n_i(ik)-n(ik))-p%T3(1)*dr1*drdx(ik)*(n_i(il)-n(il))
            fs_s=S1*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni-c_dkr(il,ik)*drdn*drdni+drdx(il)*drdx(ik)*n_dot_ni)&
                +S2*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni)+S3*drdx(il)*drdx(ik)*drdn*drdni&
                +S4*(c_dkr(il,ik)*n_dot_ni+n_i(ik)*n(il))+S5*n(ik)*n_i(il)&
                +p%S3(1)*dr2*drdx(il)*drdx(ik)*drdn*drdni&
                +(1.d0-c_dkr(il,ik))*0.125d0*p%S3(1)*dr2*(2.d0*drdx(il)*drdx(ik)*n_dot_ni+n(il)*n_i(ik)+n(ik)*n_i(il))&
                +c_dkr(il,ik)*0.25d0*p%S3(1)*dr2*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni)&
                +c_dkr(il,ik)*0.125d0*p%S3(1)*dr2*(n_dot_ni-dabs(drdt))
            ! Loop through nodes
            m(:,il,ik)=m(:,il,ik)+fs_s*phif1jw
            l(:,il,ik)=l(:,il,ik)+fs_d*phif2jw
          end do
        end do
        ! Regular and numerically evaluable parts of the CPV integrals of L
        l(:,1,2)=l(:,1,2)-p%T3(1)*dr1*drdt*(phi_f2-phi_f2_i)*jw
        l(:,2,1)=l(:,2,1)+p%T3(1)*dr1*drdt*(phi_f2-phi_f2_i)*jw
        ! Regular and numerically evaluable parts of the HFP integrals of M
        if (ksub.eq.1) then
          m(:,1,1)=m(:,1,1)+0.125d0*p%S3(1)*dr2*dabs(drdt)*(phi_f1-phi_f1_i+dphidxi_f1_i/jgi*r)*jw
          m(:,2,2)=m(:,2,2)+0.125d0*p%S3(1)*dr2*dabs(drdt)*(phi_f1-phi_f1_i+dphidxi_f1_i/jgi*r)*jw
        else
          m(:,1,1)=m(:,1,1)+0.125d0*p%S3(1)*dr2*dabs(drdt)*(phi_f1-phi_f1_i-dphidxi_f1_i/jgi*r)*jw
          m(:,2,2)=m(:,2,2)+0.125d0*p%S3(1)*dr2*dabs(drdt)*(phi_f1-phi_f1_i-dphidxi_f1_i/jgi*r)*jw
        end if
      end do ! Loop through integration points
    end do ! Loop through subdivisions
    ! Add CPV and HFP analytical terms
    ! Calculate ra and rb
    ra=sqrt((x_nodes(1,1)-x_i(1))**2+(x_nodes(2,1)-x_i(2))**2)
    rb=sqrt((x_nodes(1,2)-x_i(1))**2+(x_nodes(2,2)-x_i(2))**2)
    ! To L12 and L21
    l(:,1,2)=l(:,1,2)-p%T3(1)*phi_f2_i*(log(rb)-log(ra))
    l(:,2,1)=l(:,2,1)+p%T3(1)*phi_f2_i*(log(rb)-log(ra))
    ! To M11 and M22
    m(:,1,1)=m(:,1,1)+0.125d0*p%S3(1)*(-phi_f1_i*(1.0d0/ra+1.0d0/rb)+dphidxi_f1_i/jgi*(log(rb)-log(ra)))
    m(:,2,2)=m(:,2,2)+0.125d0*p%S3(1)*(-phi_f1_i*(1.0d0/ra+1.0d0/rb)+dphidxi_f1_i/jgi*(log(rb)-log(ra)))
    ! Multiply m and l by constants of s* and d* respectively
    m=p%cte_s*m
    l=p%cte_d*l
    ! If the normal has to be reversed, then l=-l
    if (reverse) l=-l
  end subroutine fbem_bem_harela2d_hbie_int

  subroutine fbem_bem_harela2d_hbie_auto(mode,e,reverse,x_i,n_i,p,qsp,ns,m,l)
    implicit none
    ! I/O
    integer                            :: mode              !! 0: Solve, 1: Internal points
    type(fbem_bem_element)             :: e                 !! Integration element
    logical                            :: reverse           !! Reverse orientation
    real(kind=real64)                  :: x_i(2)            !! Collocation point
    real(kind=real64)                  :: n_i(2)            !! Unit normal at the collocation point
    type(fbem_bem_harela2d_parameters) :: p                 !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp               !! Quasi-singular integration parameters
    integer                            :: ns                !! Maximum level of subdivisions
    complex(kind=real64)               :: m(e%n_pnodes,2,2) !! m integration kernel
    complex(kind=real64)               :: l(e%n_snodes,2,2) !! l integration kernel
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
        if (mode.eq.1) then
          m=0.d0
          l=0.d0
          return
        end if
      else
        delta=0
      end if
    end if
    ! Integrate
    select case (delta)
      case (1)
        call fbem_bem_harela2d_hbie_int(30,e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi(1),p,m,l)
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
          call fbem_bem_harela2d_hbie_ext_pre(ps,e,reverse,x_i,n_i,p,m,l)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_harela2d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,p,qsp,1,ns,m,l)
        end if
    end select
  end subroutine fbem_bem_harela2d_hbie_auto
  ! --------------------------------------------------------------------------------------------------------------------------------

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! BODY LOADS

  subroutine fbem_bem_harela2d_hbie_bl_ext_pre(ps,e,x_i,n_i,p,l)
    implicit none
    ! I/O
    integer                            :: ps                !! Selected precalculated dataset
    type(fbem_bem_element)             :: e                 !! Element
    real(kind=real64)                  :: x_i(2)            !! Collocation point position vector
    real(kind=real64)                  :: n_i(2)            !! Collocation point unit normal vector
    type(fbem_bem_harela2d_parameters) :: p                 !! Parameters of the region
    complex(kind=real64)               :: l(e%n_snodes,2,2) !! l kernels matrix
    ! Local
    integer              :: kip                            ! Counter variable for integration points loop
    real(kind=real64)    :: x(2)                           ! Position vector at integration point
    real(kind=real64)    :: n(2)                           ! Unit normal vector at integration point
    real(kind=real64)    :: pphijw(e%n_pnodes)             ! phi^p * jacobian * weight at integration point
    real(kind=real64)    :: sphijw(e%n_snodes)             ! phi^s * jacobian * weight at integration point
    real(kind=real64)    :: rv(2)                          ! Radiovector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, r2, d1r1, d1r2, logr        ! Radiovector module, squared, inverses and log(r)
    real(kind=real64)    :: drdx(2)                        ! Radiovector derivatives with respect to x_k
    real(kind=real64)    :: drdn                           ! Partial derivative of r respect to unit normal
    real(kind=real64)    :: drdni                          ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)    :: n_dot_ni                       ! Dot product of unit normals
    integer              :: il, ik                         ! Counter variables for load direction and observation direction
    complex(kind=real64) :: z(2)                           ! Wavenumbers
    complex(kind=real64) :: KnR(0:2,2)                     ! Bessel functions decomposition
    complex(kind=real64) :: T1, T2, T3                     ! Fundamental solution components
    complex(kind=real64) :: fs_d                           ! Fundamental solutions values
    ! Initialize
    l=(0.0d0,0.0d0)
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      x=e%ps_x(:,kip,ps)
      n=e%ps_n(:,kip,ps)
      pphijw=e%ps_pphijw(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      r2=r**2
      d1r1=1.0d0/r
      d1r2=d1r1*d1r1
      logr=log(r)
      drdx=rv*d1r1
      drdn=dot_product(drdx,n)
      drdni=-dot_product(drdx,n_i)
      n_dot_ni=dot_product(n,n_i)
      z(1)=c_im*p%k1*r
      z(2)=c_im*p%k2*r
      call fbem_BesselKnR_decomposed(2,z,KnR)
      T1=p%T1(1)*d1r1+p%T1(2)*r+p%T1(3)*KnR(1,1)+p%T1(4)*KnR(1,2)+p%T1(5)*d1r1*KnR(2,1)+p%T1(6)*d1r1*KnR(2,2)
      T2=p%T2(1)*d1r1+(p%T2(2)*logr+p%T2(3))*r+p%T2(4)*KnR(1,2)+p%T2(5)*d1r1*KnR(2,1)+p%T2(6)*d1r1*KnR(2,2)
      T3=p%T3(1)*d1r1+(p%T3(2)*logr+p%T3(3))*r+p%T3(4)*KnR(1,1)+p%T3(5)*d1r1*KnR(2,1)+p%T3(6)*d1r1*KnR(2,2)
      do il=1,2
        do ik=1,2
          fs_d=T1*drdx(il)*drdx(ik)*drdni-T2*(-drdni*c_dkr(il,ik)+drdx(il)*n_i(ik))-T3*drdx(ik)*n_i(il)
          l(:,il,ik)=l(:,il,ik)+fs_d*sphijw
        end do
      end do
    end do
    ! Multiply y constants
    l=p%cte_d*l
  end subroutine fbem_bem_harela2d_hbie_bl_ext_pre

  subroutine fbem_bem_harela2d_hbie_bl_ext_st(e,xi_s,x_i,n_i,barxip,barr,p,gln,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                 !! Integration element
    real(kind=real64)                  :: xi_s(1,2)         !! Coordinates of the subdivision of the parent element (xi space [-1,1])
    real(kind=real64)                  :: x_i(2)            !! Collocation point position vector
    real(kind=real64)                  :: n_i(2)            !! Collocation point unit normal vector
    real(kind=real64)                  :: barxip(1)         !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)                  :: barr              !! Telles jacobian at barxip
    type(fbem_bem_harela2d_parameters) :: p                 !! Parameters of the region
    integer                            :: gln               !! 1D Gauss-Legendre number of integration points (<=32)
    complex(kind=real64)               :: l(e%n_snodes,2,2) !! l kernels matrix
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
    real(kind=real64)            :: r, r2, d1r1, d1r2, logr ! Radiovector module, squared, inverses and log(r)
    real(kind=real64)            :: drdx(2)                 ! Radiovector derivatives with respect to x_k
    real(kind=real64)            :: jg                      ! Geometric jacobian
    real(kind=real64)            :: drdn                    ! Partial derivative of r respect to unit normal
    real(kind=real64)            :: drdni                   ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)            :: n_dot_ni                ! Dot product of unit normals
    real(kind=real64)            :: jw                      ! Jacobians * weights
    real(kind=real64)            :: pphijw(e%n_pnodes)      ! Auxiliary variables for integrand evaluation
    real(kind=real64)            :: sphijw(e%n_snodes)      ! Auxiliary variables for integrand evaluation
    integer                      :: il, ik                  ! Counter variables for load direction and observation direction
    complex(kind=real64)         :: z(2)                    ! Bessel functions arguments: z=ikr
    complex(kind=real64)         :: KnR(0:2,2)              ! Bessel functions decomposition
    complex(kind=real64)         :: T1, T2, T3              ! Fundamental solution components
    complex(kind=real64)         :: fs_d                    ! Fundamental solutions values
    ! Initialize
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
      r2=r**2
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
      z(1)=c_im*p%k1*r
      z(2)=c_im*p%k2*r
      call fbem_BesselKnR_decomposed(2,z,KnR)
      T1=p%T1(1)*d1r1+p%T1(2)*r+p%T1(3)*KnR(1,1)+p%T1(4)*KnR(1,2)+p%T1(5)*d1r1*KnR(2,1)+p%T1(6)*d1r1*KnR(2,2)
      T2=p%T2(1)*d1r1+(p%T2(2)*logr+p%T2(3))*r+p%T2(4)*KnR(1,2)+p%T2(5)*d1r1*KnR(2,1)+p%T2(6)*d1r1*KnR(2,2)
      T3=p%T3(1)*d1r1+(p%T3(2)*logr+p%T3(3))*r+p%T3(4)*KnR(1,1)+p%T3(5)*d1r1*KnR(2,1)+p%T3(6)*d1r1*KnR(2,2)
      ! Add
      do il=1,2
        do ik=1,2
          fs_d=T1*drdx(il)*drdx(ik)*drdni-T2*(drdx(il)*n_i(ik)-c_dkr(il,ik)*drdni)-T3*drdx(ik)*n_i(il)
          l(:,il,ik)=l(:,il,ik)+fs_d*sphijw
        end do
      end do
    end do
    l=p%cte_d*l
  end subroutine fbem_bem_harela2d_hbie_bl_ext_st

  recursive subroutine fbem_bem_harela2d_hbie_bl_ext_adp(e,xi_s,x_i,n_i,p,qsp,ks,ns,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                 !! Element
    real(kind=real64)                  :: xi_s(1,2)         !! Subdivision of the parent element
    real(kind=real64)                  :: x_i(2)            !! Collocation point position vector
    real(kind=real64)                  :: n_i(2)            !! Collocation point unit normal vector
    type(fbem_bem_harela2d_parameters) :: p                 !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp               !! Quasi-singular integration parameters
    integer                            :: ks                !! Current level of subdivisions
    integer                            :: ns                !! Maximum level of subdivision
    complex(kind=real64)               :: l(e%n_snodes,2,2) !! l integration kernels matrix
    ! Local
    integer              :: gln_near              ! 1D Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer              :: gln                   ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical              :: subdivide             ! True if subdivision has to be performed
    real(kind=real64)    :: barxi(1)              ! Nearest element coordinate with respect to collocation point
    real(kind=real64)    :: barxip(1)             ! Nearest element subdivision local coordinate with respect to collocation point
    real(kind=real64)    :: rmin                  ! Minimum distance between collocation point and barxi on the element
    real(kind=real64)    :: barr                  ! Telles jacobian at barxi
    real(kind=real64)    :: cl                    ! Characteristic length
    real(kind=real64)    :: d                     ! Normalized distance between collocation point and element subdivision
    integer              :: method                ! Method used in nearest point algorithm
    real(kind=real64)    :: tmp_xi_s(1,2)         ! Subdivision
    real(kind=real64)    :: x_s(2,e%n_gnodes)     ! Coordinates of the element subdivision
    complex(kind=real64) :: l_tmp(e%n_snodes,2,2) ! l kernels matrix (temporary)
    ! Initialize
    if (ks.eq.1) then
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
        call fbem_warning_message(error_unit,0,'fbem_bem_harela2d_hbie_bl_ext_adp',ns,'maximum number of subdivisions reached')
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
      call fbem_bem_harela2d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,l)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_harela2d_hbie_bl_ext_adp(e,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,l)
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_harela2d_hbie_bl_ext_st(e,xi_s,x_i,n_i,barxip,barr,p,gln,l_tmp)
      l=l+l_tmp
    end if
  end subroutine fbem_bem_harela2d_hbie_bl_ext_adp

  !! Calculation of m and l when the collocation point is on the integration element.
  subroutine fbem_bem_harela2d_hbie_bl_int(ngp,type_g,type_f1,type_f2,delta_f,x_nodes,xi_i,p,l)
    implicit none
    ! I/O
    integer                            :: ngp                                 !! Number of Gauss point to be used (<=32)
    integer                            :: type_g                              !! Geometrical interpolation
    integer                            :: type_f1                             !! Functional interpolation (primary variables)
    integer                            :: type_f2                             !! Geometrical interpolation (secondary variables)
    real(kind=real64)                  :: delta_f                             !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64)                  :: x_nodes(2,fbem_n_nodes(type_g))     !! Position vectors of geometrical nodes
    real(kind=real64)                  :: xi_i                                !! Reference coordinate of the singular point.
    type(fbem_bem_harela2d_parameters) :: p                               !! Parameters of the region
    complex(kind=real64)               :: l(fbem_n_nodes(type_f2),2,2)        !! l integration kernels matrix
    ! Local
    integer                            :: nnodes_g                            ! Number of nodes of the element interpolation
    integer                            :: nnodes_f1                           ! Number of nodes of the element interpolation
    integer                            :: nnodes_f2                           ! Number of nodes of the element interpolation
    integer                            :: kphi                                ! Counter variable for shape functions loops
    integer                            :: kip                                 ! Counter of integration points
    real(kind=real64)                  :: x_i(2)                              ! Real coordinates of collocation point
    real(kind=real64)                  :: t_i(2)                              ! Unit tangent at collocation point
    real(kind=real64)                  :: n_i(2)                              ! Unit normal at collocation point
    real(kind=real64)                  :: phi_f1(fbem_n_nodes(type_f1))       ! Functional shape functions values at xi
    real(kind=real64)                  :: phi_f1_i(fbem_n_nodes(type_f1))     ! Functional shape functions values at xi_i
    real(kind=real64)                  :: dphidxi_f1_i(fbem_n_nodes(type_f1)) ! Functional shape functions derivatives values at xi_i
    real(kind=real64)                  :: phi_f2(fbem_n_nodes(type_f2))       ! Functional shape functions values at xi
    real(kind=real64)                  :: phi_f2_i(fbem_n_nodes(type_f2))     ! Functional shape functions values at xi_i
    real(kind=real64)                  :: phi_g(fbem_n_nodes(type_g))         ! Geometrical shape functions values at xi
    real(kind=real64)                  :: dphidxi_g(fbem_n_nodes(type_g))     ! Geometrical shape functions derivatives values at xi
    integer                            :: nsub                                ! Number of subdivision of the element
    integer                            :: ksub                                ! Counter of subdivision
    real(kind=real64)                  :: gamma                               ! Coordinate gamma
    real(kind=real64)                  :: w                                   ! Weights of each integration point
    real(kind=real64)                  :: xip                                 ! Coordinate  xip of subdivided element [0,1]
    real(kind=real64)                  :: js                                  ! Jacobian of the xi [xisub(1,:),xisub[2,:]] -> xip [0,1] transformation
    real(kind=real64)                  :: xip_i(2)                            ! Singular point in xip space
    real(kind=real64)                  :: xi                                  ! Coordinate xi
    real(kind=real64)                  :: xisub(2,2)                          ! Coordinates of subdivisions in xi space
    real(kind=real64)                  :: aux(10)                             ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)                  :: x(2)                                ! Position vector at xi
    real(kind=real64)                  :: T(2)                                ! Tangent vector at xi
    real(kind=real64)                  :: N(2)                                ! Normal vector at xi
    real(kind=real64)                  :: rv(2)                               ! Radiovector between collocation point and integration point (x-x_i)
    real(kind=real64)                  :: r, r2, dr1, dr2, logr               ! Radiovector module, squared, inverses and log(r)
    real(kind=real64)                  :: ra, rb                              ! Radiovector from collocation point to element vertices
    real(kind=real64)                  :: drdx(2)                             ! Radiovector derivatives with respect to x_k
    real(kind=real64)                  :: jg                                  ! Geometric jacobian
    real(kind=real64)                  :: jgi                                 ! Geometric jacobian at the collocation point
    real(kind=real64)                  :: drdn                                ! Partial derivative of r respect to unit normal
    real(kind=real64)                  :: drdni                               ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)                  :: n_dot_ni                            ! Dot product of unit normals
    real(kind=real64)                  :: drdt                                ! Partial derivative of r respect to unit tangent
    real(kind=real64)                  :: jw                                  ! Jacobians*weight
    real(kind=real64)                  :: phif1jw(fbem_n_nodes(type_f1))      ! Auxiliary variables for integrand evaluation
    real(kind=real64)                  :: phif2jw(fbem_n_nodes(type_f2))      ! Auxiliary variables for integrand evaluation
    type(fbem_telles_parameters)       :: telles_parameters                   ! Telles parameters
    real(kind=real64)                  :: jt                                  ! Telles jacobian
    integer                            :: il, ik                              ! Counter variables for load direction and observation direction
    complex(kind=real64)               :: z(2)                                ! Wavenumbers
    complex(kind=real64)               :: KnR(0:2,2)                          ! Bessel functions decomposition
    complex(kind=real64)               :: T1, T2, T3                          ! Fundamental solution components
    complex(kind=real64)               :: fs_d                                ! Fundamental solutions values
    !
    ! Initialize
    !
    ! Initialize kernel matrices
    l=(0.0d0,0.0d0)
    ! Number of nodes of the element
    nnodes_g=fbem_n_nodes(type_g)
    nnodes_f1=fbem_n_nodes(type_f1)
    nnodes_f2=fbem_n_nodes(type_f2)
    ! If xi_i belongs to one of the vertices, subdivision is not needed
    if (fbem_check_xi_vertex(xi_i).eqv.(.true.)) then
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                              'the HBIE cannot be collocated at a vertex')
    ! If xi_i is inside the element, 2 subdivisions are needed
    else
      nsub=2
      ! Coordinates xi of the subdivision 1
      xisub(1,1)=-1.0d0
      xisub(2,1)=xi_i
      ! Coordinate xip of the collocation point
      xip_i(1)=1.0d0
      ! Coordinates xi of the subdivision 2
      xisub(1,2)=xi_i
      xisub(2,2)=1.0d0
      ! Coordinate xip of the collocation point
      xip_i(2)=0.0d0
    end if
    !
    ! Calculate real coordinates, geometric jacobian, unit tangent and unit normal at collocation point
    !
    ! Geometrical shape functions and shape functions derivatives at xi_i
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
    ! x_i and T_i
    x_i=0.0d0
    T_i=0.0d0
    do kphi=1,nnodes_g
      x_i=x_i+phi_g(kphi)*x_nodes(:,kphi)
      T_i=T_i+dphidxi_g(kphi)*x_nodes(:,kphi)
    end do
    ! Normal vector
    N_i(1)=T_i(2)
    N_i(2)=-T_i(1)
    ! Geometric jacobian
    jgi=sqrt(T_i(1)**2+T_i(2)**2)
    ! Unit normal vector at collocation point
    n_i=N_i/jgi
    ! Unit tangent vector at collocation point
    t_i=T_i/jgi
    !
    ! Functional shape functions and its first derivatives (primary variable) at xi_i
    !
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
    !
    ! Functional shape functions (secondary variable) at xi_i
    !
#   define etype type_f2
#   define delta delta_f
#   define xi xi_i
#   define phi phi_f2_i
#   include <phi_1d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
    !
    ! INTEGRATE REGULAR INTEGRALS
    !
    ! Loop through subdivisions
    do ksub=1,nsub
      ! Jacobian of xip->xi transformation (is constant)
      js=xisub(2,ksub)-xisub(1,ksub)
      ! Loop through integration points
      do kip=1,gl01_n(ngp)
        ! Coordinate and weight in xip [0,1]
        xip=gl01_xi(kip,ngp)
        w=gl01_w(kip,ngp)
        ! Calculate xi
        xi=js*xip+xisub(1,ksub)
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
        ! Unit normal vector
        n=N/jg
        ! Unit tangent vector
        t=T/jg
        ! Radiovector
        rv=x-x_i
        ! Radiovector norm and its inverse
        r=dot_product(rv,rv)
        r=sqrt(r)
        r2=r**2
        dr1=1.d0/r
        dr2=dr1**2
        logr=log(r)
        ! r_{,k}
        drdx=rv*dr1
        ! dr/dn
        drdn=dot_product(drdx,n)
        ! dr/dn_i
        drdni=-dot_product(drdx,n_i)
        ! n · n_i
        n_dot_ni=dot_product(n,n_i)
        ! dr/dGamma
        drdt=dot_product(drdx,t)
        ! Functional shape functions (secondary variable) at xi
#       define etype type_f2
#       define delta delta_f
#       define phi phi_f2
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Jacobian * weight
        jw=jg*js*w
        ! Functional shape functions * jacobians* weight
        phif1jw=phi_f1*jw
        phif2jw=phi_f2*jw
        ! COMPONENTS OF THE FUNDAMENTAL SOLUTION
        z(1)=c_im*p%k1*r
        z(2)=c_im*p%k2*r
        call fbem_BesselKnR_decomposed(2,z,KnR)
        T1=p%T1(2)*r+p%T1(3)*KnR(1,1)+p%T1(4)*KnR(1,2)+p%T1(5)*dr1*KnR(2,1)+p%T1(6)*dr1*KnR(2,2)
        T2=(p%T2(2)*logr+p%T2(3))*r+p%T2(4)*KnR(1,2)+p%T2(5)*dr1*KnR(2,1)+p%T2(6)*dr1*KnR(2,2)
        T3=(p%T3(2)*logr+p%T3(3))*r+p%T3(4)*KnR(1,1)+p%T3(5)*dr1*KnR(2,1)+p%T3(6)*dr1*KnR(2,2)
        ! Loop through load direction
        do il=1,2
          ! Loop through observation direction
          do ik=1,2
            ! Fundamental solutions values without constants
            fs_d=T1*drdx(il)*drdx(ik)*drdni-T2*(drdx(il)*n_i(ik)-c_dkr(il,ik)*drdni)-T3*drdx(ik)*n_i(il)&
                +p%T1(1)*dr1*drdx(il)*drdx(ik)*drdni+p%T2(1)*dr1*c_dkr(il,ik)*drdni&
                -p%T2(1)*dr1*drdx(il)*(n_i(ik)-n(ik))-p%T3(1)*dr1*drdx(ik)*(n_i(il)-n(il))
            ! Loop through nodes
            l(:,il,ik)=l(:,il,ik)+fs_d*phif2jw
          end do
        end do
        ! Regular and numerically evaluable parts of the CPV integrals of L
        l(:,1,2)=l(:,1,2)-p%T3(1)*dr1*drdt*(phi_f2-phi_f2_i)*jw
        l(:,2,1)=l(:,2,1)+p%T3(1)*dr1*drdt*(phi_f2-phi_f2_i)*jw
      end do ! Loop through integration points
    end do ! Loop through subdivisions
    ! Add CPV and HFP analytical terms
    ! Calculate ra and rb
    ra=sqrt((x_nodes(1,1)-x_i(1))**2+(x_nodes(2,1)-x_i(2))**2)
    rb=sqrt((x_nodes(1,2)-x_i(1))**2+(x_nodes(2,2)-x_i(2))**2)
    ! To L12 and L21
    l(:,1,2)=l(:,1,2)-p%T3(1)*phi_f2_i*(log(rb)-log(ra))
    l(:,2,1)=l(:,2,1)+p%T3(1)*phi_f2_i*(log(rb)-log(ra))
    ! Multiply m and l by constants of s* and d* respectively
    l=p%cte_d*l
  end subroutine fbem_bem_harela2d_hbie_bl_int

  subroutine fbem_bem_harela2d_hbie_bl_auto(mode,e,x_i,n_i,p,qsp,ns,l)
    implicit none
    ! I/O
    integer                            :: mode              !! 0: Solve, 1: Internal points
    type(fbem_bem_element)             :: e                 !! Integration element
    real(kind=real64)                  :: x_i(2)            !! Collocation point
    real(kind=real64)                  :: n_i(2)            !! Unit normal at the collocation point
    type(fbem_bem_harela2d_parameters) :: p                 !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp               !! Quasi-singular integration parameters
    integer                            :: ns                !! Maximum level of subdivisions
    complex(kind=real64)               :: l(e%n_snodes,2,2) !! l integration kernel
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
        call fbem_error_message(output_unit,0,'fbem_bem_harela2d_hbie_bl_auto',0,'it is not possible to collocate at a point load')
      else
        call fbem_bem_harela2d_hbie_d(e%x(:,1),x_i,n_i,p,l(1,:,:))
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
        call fbem_nearest_element_point_bem(2,e%gtype,e%x,e%cl,x_i,barxi,rmin,d,method)
        if (d.le.1.d-12) then
          delta=1
          if (mode.eq.1) then
            l=0.d0
            return
          end if
        else
          delta=0
        end if
      end if
      ! Integrate
      select case (delta)
        case (1)
          call fbem_bem_harela2d_hbie_bl_int(30,e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,barxi(1),p,l)
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
            call fbem_bem_harela2d_hbie_bl_ext_pre(ps,e,x_i,n_i,p,l)
          ! Integrate using an adaptative algorithm
          else
            call fbem_bem_harela2d_hbie_bl_ext_adp(e,xi_s,x_i,n_i,p,qsp,1,ns,l)
          end if
      end select
    end if
  end subroutine fbem_bem_harela2d_hbie_bl_auto
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

  !! Calculation of h, g, m and l when the collocation point is outside the integration element, using a precalculated dataset of
  !! data.
  subroutine fbem_bem_harela2d_shbie_ext_pre(ps,e,reverse,x_i,n_i,p,h,g,m,l)
    implicit none
    ! I/O
    integer                            :: ps                !! Selected precalculated dataset
    type(fbem_bem_element)             :: e                 !! Element
    logical                            :: reverse           !! Reverse normal vector
    real(kind=real64)                  :: x_i(2)            !! Collocation point position vector
    real(kind=real64)                  :: n_i(2)            !! Collocation point unit normal vector
    type(fbem_bem_harela2d_parameters) :: p                 !! Parameters of the region
    complex(kind=real64)               :: h(e%n_pnodes,2,2) !! h integration kernels matrix
    complex(kind=real64)               :: g(e%n_snodes,2,2) !! g integration kernels matrix
    complex(kind=real64)               :: m(e%n_pnodes,2,2) !! m integration kernels matrix
    complex(kind=real64)               :: l(e%n_snodes,2,2) !! l integration kernels matrix
    ! Local
    integer              :: kip                     ! Counter variable for integration points loop
    real(kind=real64)    :: x(2)                    ! Position vector at integration point
    real(kind=real64)    :: n(2)                    ! Unit normal vector at integration point
    real(kind=real64)    :: pphijw(e%n_pnodes)      ! phi^p * jacobian * weight at integration point
    real(kind=real64)    :: sphijw(e%n_snodes)      ! phi^s * jacobian * weight at integration point
    real(kind=real64)    :: rv(2)                   ! Radiovector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, r2, d1r1, d1r2, logr ! Radiovector module, squared, inverses and log(r)
    real(kind=real64)    :: drdx(2)                 ! Radiovector derivatives with respect to x_k
    real(kind=real64)    :: drdn                    ! Partial derivative of r respect to unit normal
    real(kind=real64)    :: drdni                   ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)    :: n_dot_ni                ! Dot product of unit normals
    integer              :: il, ik                  ! Counter variables for load direction and observation direction
    complex(kind=real64) :: z(2)                    ! Wavenumbers
    complex(kind=real64) :: KnR(0:2,2)              ! Bessel functions decomposition
    complex(kind=real64) :: psi, chi                ! Fundamental solution components
    complex(kind=real64) :: T1, T2, T3              ! Fundamental solution components
    complex(kind=real64) :: S1, S2, S3, S4, S5      ! Fundamental solution components
    complex(kind=real64) :: fs_u, fs_t, fs_d, fs_s  ! Fundamental solutions values
    ! Initialize
    h=(0.0d0,0.0d0)
    g=(0.0d0,0.0d0)
    m=(0.0d0,0.0d0)
    l=(0.0d0,0.0d0)
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      ! Functional shape functions * jacobian * weight
      pphijw=e%ps_pphijw(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      ! Components of fundamental solutions
      x=e%ps_x(:,kip,ps)
      n=e%ps_n(:,kip,ps)
      rv=x-x_i
      r=sqrt(dot_product(rv,rv))
      r2=r**2
      d1r1=1.0d0/r
      d1r2=d1r1*d1r1
      logr=log(r)
      drdx=rv*d1r1
      drdn=dot_product(drdx,n)
      drdni=-dot_product(drdx,n_i)
      n_dot_ni=dot_product(n,n_i)
      z(1)=c_im*p%k1*r
      z(2)=c_im*p%k2*r
      call fbem_BesselKnR_decomposed(2,z,KnR)
      psi=p%psi(1)*logr+p%psi(2)+KnR(0,2)+d1r1*p%psi(3)*KnR(1,1)+d1r1*p%psi(4)*KnR(1,2)
      chi=p%chi(1)+(p%chi(2)*logr+p%chi(3))*r2+p%chi(4)*KnR(2,1)+KnR(2,2)
      T1=p%T1(1)*d1r1+p%T1(2)*r+p%T1(3)*KnR(1,1)+p%T1(4)*KnR(1,2)+p%T1(5)*d1r1*KnR(2,1)+p%T1(6)*d1r1*KnR(2,2)
      T2=p%T2(1)*d1r1+(p%T2(2)*logr+p%T2(3))*r+p%T2(4)*KnR(1,2)+p%T2(5)*d1r1*KnR(2,1)+p%T2(6)*d1r1*KnR(2,2)
      T3=p%T3(1)*d1r1+(p%T3(2)*logr+p%T3(3))*r+p%T3(4)*KnR(1,1)+p%T3(5)*d1r1*KnR(2,1)+p%T3(6)*d1r1*KnR(2,2)
      S1=p%S1(1)*d1r2+p%S1(2)+p%S1(3)*KnR(0,2)+p%S1(4)*d1r1*KnR(1,1)+p%S1(5)*d1r1*KnR(1,2)+p%S1(6)*d1r2*KnR(2,1)&
        +p%S1(7)*d1r2*KnR(2,2)
      S2=p%S2(1)*d1r2+p%S2(2)+(p%S2(3)*logr+p%S2(4))*r2+p%S2(5)*d1r1*KnR(1,1)+p%S2(6)*d1r1*KnR(1,2)+p%S2(7)*d1r2*KnR(2,1)&
        +p%S2(8)*KnR(2,1)+p%S2(9)*d1r2*KnR(2,2)
      S3=p%S3(1)*d1r2+p%S3(2)+p%S3(3)*KnR(0,1)+p%S3(4)*KnR(0,2)+p%S3(5)*d1r1*KnR(1,1)+p%S3(6)*d1r1*KnR(1,2)&
        +p%S3(7)*d1r2*KnR(2,1)+p%S3(8)*d1r2*KnR(2,2)
      S4=p%S4(1)*d1r2+p%S4(2)*logr+p%S4(3)+p%S4(4)*d1r1*KnR(1,2)+p%S4(5)*d1r2*KnR(2,1)+p%S4(6)*d1r2*KnR(2,2)
      S5=p%S5(1)*d1r2+p%S5(2)*logr+p%S5(3)+p%S5(4)*KnR(0,1)+p%S5(5)*d1r1*KnR(1,1)+p%S5(6)*d1r2*KnR(2,1)+p%S5(7)*d1r2*KnR(2,2)
      ! Add the term of the numerical integration summation (fundamental solutions * functional shape functions * jacobian * weight)
      do il=1,2
        do ik=1,2
          ! Fundamental solutions
          fs_u=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
          fs_t=T1*drdx(il)*drdx(ik)*drdn+T2*(drdn*c_dkr(il,ik)+drdx(ik)*n(il))+T3*drdx(il)*n(ik)
          fs_d=T1*drdx(il)*drdx(ik)*drdni-T2*(-drdni*c_dkr(il,ik)+drdx(il)*n_i(ik))-T3*drdx(ik)*n_i(il)
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
  end subroutine fbem_bem_harela2d_shbie_ext_pre

  !! Efficient automatic calculation of h, g, m and l.
  subroutine fbem_bem_harela2d_shbie_auto(e,reverse,x_i,n_i,p,qsp,ns,h,g,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                 !! Integration element
    logical                            :: reverse           !! Reverse orientation
    real(kind=real64)                  :: x_i(2)            !! Collocation point
    real(kind=real64)                  :: n_i(2)            !! Unit normal at the collocation point
    type(fbem_bem_harela2d_parameters) :: p                 !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp               !! Quasi-singular integration parameters
    integer                            :: ns                !! Maximum level of subdivisions
    complex(kind=real64)               :: h(e%n_pnodes,2,2) !! h integration kernel
    complex(kind=real64)               :: g(e%n_snodes,2,2) !! g integration kernel
    complex(kind=real64)               :: m(e%n_pnodes,2,2) !! m integration kernel
    complex(kind=real64)               :: l(e%n_snodes,2,2) !! l integration kernel
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
        call fbem_bem_harela2d_sbie_int(30,e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi(1),p,h,g)
        call fbem_bem_harela2d_hbie_int(30,e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi(1),p,m,l)
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
          call fbem_bem_harela2d_shbie_ext_pre(ps,e,reverse,x_i,n_i,p,h,g,m,l)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_harela2d_sbie_ext_adp(e,reverse,xi_s,x_i,p,qsp,1,ns,h,g)
          call fbem_bem_harela2d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,p,qsp,1,ns,m,l)
        end if
    end select
  end subroutine fbem_bem_harela2d_shbie_auto

  ! ================================================================================================================================

  ! ================================================================================================================================
  ! VARIATION SINGULAR BOUNDARY INTEGRAL EQUATION (VSBIE)

  !! Function that calculates one of the free-terms of the VSBIE of 2D elastic problem
  subroutine fbem_bem_harela2d_vsbie_freeterm(n_elements,n,tc,tol,nu,b)
    implicit none
    ! I/O
    integer              :: n_elements !! Number of elements (1 or 2)
    real(kind=real64)    :: n(2,2)     !! Normals of elements at collocation point
    real(kind=real64)    :: tc(2,2)    !! Tangents of elements at collocation point towards inside the element
    real(kind=real64)    :: tol        !! Geometric tolerance [1.0e-12,1.0e-3] (default 1.0e-6)
    complex(kind=real64) :: nu         !! Poisson's ratio
    complex(kind=real64) :: b(2,2,2,2) !! Free-term b_{l,k,j,m}
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
  end subroutine fbem_bem_harela2d_vsbie_freeterm

  !! Calculation of h1, h2, g1 and g2 when the collocation point is outside the integration element, using a precalculated dataset
  !! of data.
  subroutine fbem_bem_harela2d_vsbie_ext_pre(ps,e,reverse,x_i,p,h1,h2,g1,g2)
    implicit none
    ! I/O
    integer                            :: ps                                  !! Selected precalculated dataset
    type(fbem_bem_element)             :: e                                   !! Element
    logical                            :: reverse                             !! Reverse normal vector
    real(kind=real64)                  :: x_i(2)                              !! Collocation point position vector
    type(fbem_bem_harela2d_parameters) :: p                                   !! Parameters of the region
    complex(kind=real64)               :: h1(e%dme_n_gnodes,e%n_pnodes,2,2,2) !! h1 matrix
    complex(kind=real64)               :: h2(               e%n_pnodes,2,2,2) !! h2 matrix
    complex(kind=real64)               :: g1(e%dme_n_gnodes,e%n_snodes,2,2,2) !! g1 matrix
    complex(kind=real64)               :: g2(               e%n_snodes,2,2,2) !! g2 matrix
    ! Local
    integer              :: il, ik, im, iq                   ! Counter variables
    integer              :: kip                              ! Counter variable for integration points loop
    real(kind=real64)    :: x(2), t(2), n(2)                 ! Position, unit tangent and unit normal vectors
    real(kind=real64)    :: pphijw(e%n_pnodes)               ! phi^p*jacobian*weight at integration point
    real(kind=real64)    :: sphijw(e%n_snodes)               ! phi^s*jacobian*weight at integration point
    real(kind=real64)    :: rv(2), r, dr1, dr2, logr         ! Distance vector
    real(kind=real64)    :: drdx(2)                          ! Distance vector derivatives with respect to x_k
    real(kind=real64)    :: drdn                             ! Partial derivative of r respect to unit normal
    real(kind=real64)    :: drdt                             ! Partial derivative of r respect to unit tangent
    real(kind=real64)    :: dme_gphi(e%dme_n_gnodes)         ! phi^{dme-g} at integration point
    real(kind=real64)    :: dme_wqj(e%dme_n_gnodes,2)        ! dphi{dme-g}/dx_j at integration point
    real(kind=real64)    :: dme_vq(e%dme_n_gnodes)           ! vq=w_{j,q}·t_j at integration point
    complex(kind=real64) :: z(2)                             ! Wavenumbers
    complex(kind=real64) :: KnR(0:2,2)                       ! Bessel functions decomposition
    complex(kind=real64) :: psi, chi, T1, T2, T3             ! Fundamental solution components
    complex(kind=real64) :: U1, U2, U3                       ! Fundamental solution components
    complex(kind=real64) :: R1, R2, R3, R4, R5, R6           ! Fundamental solution components
    complex(kind=real64) :: fs_u, fs_t                       ! Fundamental solutions values
    complex(kind=real64) :: fs_dudm, fs_tt, fs_dsigdm        ! Fundamental solutions values
    ! Initialize matrices
    h1=(0.d0,0.d0)
    h2=(0.d0,0.d0)
    g1=(0.d0,0.d0)
    g2=(0.d0,0.d0)
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
      dr2=dr1**2
      logr=log(r)
      drdx=rv*dr1
      drdn=dot_product(drdx,n)
      drdt=dot_product(drdx,t)
      ! Calculation of the components of the fundamental solution
      z(1)=c_im*p%k1*r
      z(2)=c_im*p%k2*r
      call fbem_BesselKnR_decomposed(2,z,KnR)
      psi=p%psi(1)*logr+p%psi(2)+KnR(0,2)+dr1*p%psi(3)*KnR(1,1)+dr1*p%psi(4)*KnR(1,2)
      chi=p%chi(1)+(p%chi(2)*logr+p%chi(3))*r**2+p%chi(4)*KnR(2,1)+KnR(2,2)
      T1=p%T1(1)*dr1+p%T1(2)*r+p%T1(3)*KnR(1,1)+p%T1(4)*KnR(1,2)+p%T1(5)*dr1*KnR(2,1)+p%T1(6)*dr1*KnR(2,2)
      T2=p%T2(1)*dr1+(p%T2(2)*logr+p%T2(3))*r+p%T2(4)*KnR(1,2)+p%T2(5)*dr1*KnR(2,1)+p%T2(6)*dr1*KnR(2,2)
      T3=p%T3(1)*dr1+(p%T3(2)*logr+p%T3(3))*r+p%T3(4)*KnR(1,1)+p%T3(5)*dr1*KnR(2,1)+p%T3(6)*dr1*KnR(2,2)
      U1=p%U1(1)*dr1+(p%U1(2)*logr+p%U1(3))*r+p%U1(4)*KnR(1,2)+p%U1(5)*dr1*KnR(2,1)+p%U1(6)*dr1*KnR(2,2)
      U2=p%U2(1)*dr1+p%U2(2)*r+p%U2(3)*KnR(1,1)+p%U2(4)*KnR(1,2)+p%U2(5)*dr1*KnR(2,1)+p%U2(6)*dr1*KnR(2,2)
      U3=p%U3(1)*dr1+(p%U3(2)*logr+p%U3(3))*r+p%U3(4)*dr1*KnR(2,1)+p%U3(5)*dr1*KnR(2,2)
      R1=p%R1(1)*dr2+p%R1(2)+p%R1(3)*dr1*KnR(1,1)+p%R1(4)*dr1*KnR(1,2)+p%R1(5)*dr2*KnR(2,1)+p%R1(6)*dr2*KnR(2,2)
      R2=p%R2(1)*dr2+p%R2(2)+p%R2(3)*KnR(0,2)+p%R2(4)*dr1*KnR(1,1)+p%R2(5)*dr1*KnR(1,2)+p%R2(6)*dr2*KnR(2,1)+p%R2(7)*dr2*KnR(2,2)
      R3=p%R3(1)*dr2+p%R3(2)+p%R3(3)*KnR(0,1)+p%R3(4)*KnR(0,2)+p%R3(5)*dr1*KnR(1,1)+p%R3(6)*dr1*KnR(1,2)+p%R3(7)*dr2*KnR(2,1)&
        +p%R3(8)*dr2*KnR(2,2)
      R4=p%R4(1)*dr2+p%R4(2)*logr+p%R4(3)+p%R4(4)*dr1*KnR(1,1)+p%R4(5)*dr2*KnR(2,1)+p%R4(6)*dr2*KnR(2,2)
      R5=p%R5(1)*dr2+p%R5(2)+p%R5(3)*KnR(0,1)+p%R5(4)*dr1*KnR(1,1)+p%R5(5)*dr1*KnR(1,2)+p%R5(6)*dr2*KnR(2,1)+p%R5(7)*dr2*KnR(2,2)
      R6=p%R6(1)*dr2+p%R6(2)*logr+p%R6(3)+p%R6(4)*dr1*KnR(1,2)+p%R6(5)*dr2*KnR(2,1)+p%R6(6)*dr2*KnR(2,2)
      ! Add integration points
      do il=1,2
        do ik=1,2
          ! u*_{lk}
          fs_u=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
          ! t*_{lk}
          fs_t=T1*drdx(il)*drdx(ik)*drdn+T2*(drdn*c_dkr(il,ik)+drdx(ik)*n(il))+T3*drdx(il)*n(ik)
          ! sigma*_{lkj}·t_j
          fs_tt=T1*drdx(il)*drdx(ik)*drdt+T2*(drdt*c_dkr(il,ik)+drdx(ik)*t(il))+T3*drdx(il)*t(ik)
          do im=1,2
            ! u*_{lk,m}
            fs_dudm=U1*drdx(im)*c_dkr(il,ik)+U2*drdx(il)*drdx(ik)*drdx(im)+U3*(c_dkr(il,im)*drdx(ik)+c_dkr(ik,im)*drdx(il))
            ! sigma*_{lkj,m}·n_{j}
            fs_dsigdm=R1*(c_dkr(il,im)*drdx(ik)*drdn+c_dkr(ik,im)*drdx(il)*drdn+drdx(il)*drdx(ik)*n(im))&
                     +R2*(c_dkr(il,ik)*drdx(im)*drdn+n(il)*drdx(ik)*drdx(im))+R3*drdx(il)*drdx(ik)*drdx(im)*drdn&
                     +R4*c_dkr(il,im)*n(ik)+R5*n(ik)*drdx(il)*drdx(im)+R6*(c_dkr(il,ik)*n(im)+c_dkr(ik,im)*n(il))
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
    h1=p%cte_t*h1
    h2=p%cte_t*h2
    g1=p%cte_u*g1
    g2=p%cte_u*g2
    ! Assemble depending on reversion or not
    if (reverse) then
      h1=-h1
      h2=-h2
    end if
  end subroutine fbem_bem_harela2d_vsbie_ext_pre

  !! Calculation of h1, h2, g1 and g2 over a subdivision of the integration element when the collocation point is outside but near
  !! the integration element. It uses Telles transformation with a given Telles jacobian.
  subroutine fbem_bem_harela2d_vsbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,p,gln,h1,h2,g1,g2)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                                   !! Integration element
    logical                            :: reverse                             !! Reverse normal vector
    real(kind=real64)                  :: xi_s(1,2)                           !! Coordinates of the subdivision of the parent element (xi space [-1,1])
    real(kind=real64)                  :: x_i(2)                              !! Collocation point position vector
    real(kind=real64)                  :: barxip(1)                           !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)                  :: barr                                !! Telles jacobian at barxip
    type(fbem_bem_harela2d_parameters) :: p                                   !! Parameters of the region
    integer                            :: gln                                 !! 1D Gauss-Legendre number of integration points (<=32)
    complex(kind=real64)               :: h1(e%dme_n_gnodes,e%n_pnodes,2,2,2) !! h1 matrix
    complex(kind=real64)               :: h2(               e%n_pnodes,2,2,2) !! h2 matrix
    complex(kind=real64)               :: g1(e%dme_n_gnodes,e%n_snodes,2,2,2) !! g1 matrix
    complex(kind=real64)               :: g2(               e%n_snodes,2,2,2) !! g2 matrix
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
    complex(kind=real64)         :: z(2)                             ! Wavenumbers
    complex(kind=real64)         :: KnR(0:2,2)                       ! Bessel functions decomposition
    complex(kind=real64)         :: psi, chi, T1, T2, T3             ! Fundamental solution components
    complex(kind=real64)         :: U1, U2, U3                       ! Fundamental solution components
    complex(kind=real64)         :: R1, R2, R3, R4, R5, R6           ! Fundamental solution components
    complex(kind=real64)         :: fs_u, fs_t                       ! Fundamental solutions values
    complex(kind=real64)         :: fs_dudm, fs_tt, fs_dsigdm        ! Fundamental solutions values
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
      z(1)=c_im*p%k1*r
      z(2)=c_im*p%k2*r
      call fbem_BesselKnR_decomposed(2,z,KnR)
      psi=p%psi(1)*logr+p%psi(2)+KnR(0,2)+dr1*p%psi(3)*KnR(1,1)+dr1*p%psi(4)*KnR(1,2)
      chi=p%chi(1)+(p%chi(2)*logr+p%chi(3))*r**2+p%chi(4)*KnR(2,1)+KnR(2,2)
      T1=p%T1(1)*dr1+p%T1(2)*r+p%T1(3)*KnR(1,1)+p%T1(4)*KnR(1,2)+p%T1(5)*dr1*KnR(2,1)+p%T1(6)*dr1*KnR(2,2)
      T2=p%T2(1)*dr1+(p%T2(2)*logr+p%T2(3))*r+p%T2(4)*KnR(1,2)+p%T2(5)*dr1*KnR(2,1)+p%T2(6)*dr1*KnR(2,2)
      T3=p%T3(1)*dr1+(p%T3(2)*logr+p%T3(3))*r+p%T3(4)*KnR(1,1)+p%T3(5)*dr1*KnR(2,1)+p%T3(6)*dr1*KnR(2,2)
      U1=p%U1(1)*dr1+(p%U1(2)*logr+p%U1(3))*r+p%U1(4)*KnR(1,2)+p%U1(5)*dr1*KnR(2,1)+p%U1(6)*dr1*KnR(2,2)
      U2=p%U2(1)*dr1+p%U2(2)*r+p%U2(3)*KnR(1,1)+p%U2(4)*KnR(1,2)+p%U2(5)*dr1*KnR(2,1)+p%U2(6)*dr1*KnR(2,2)
      U3=p%U3(1)*dr1+(p%U3(2)*logr+p%U3(3))*r+p%U3(4)*dr1*KnR(2,1)+p%U3(5)*dr1*KnR(2,2)
      R1=p%R1(1)*dr2+p%R1(2)+p%R1(3)*dr1*KnR(1,1)+p%R1(4)*dr1*KnR(1,2)+p%R1(5)*dr2*KnR(2,1)+p%R1(6)*dr2*KnR(2,2)
      R2=p%R2(1)*dr2+p%R2(2)+p%R2(3)*KnR(0,2)+p%R2(4)*dr1*KnR(1,1)+p%R2(5)*dr1*KnR(1,2)+p%R2(6)*dr2*KnR(2,1)+p%R2(7)*dr2*KnR(2,2)
      R3=p%R3(1)*dr2+p%R3(2)+p%R3(3)*KnR(0,1)+p%R3(4)*KnR(0,2)+p%R3(5)*dr1*KnR(1,1)+p%R3(6)*dr1*KnR(1,2)+p%R3(7)*dr2*KnR(2,1)&
        +p%R3(8)*dr2*KnR(2,2)
      R4=p%R4(1)*dr2+p%R4(2)*logr+p%R4(3)+p%R4(4)*dr1*KnR(1,1)+p%R4(5)*dr2*KnR(2,1)+p%R4(6)*dr2*KnR(2,2)
      R5=p%R5(1)*dr2+p%R5(2)+p%R5(3)*KnR(0,1)+p%R5(4)*dr1*KnR(1,1)+p%R5(5)*dr1*KnR(1,2)+p%R5(6)*dr2*KnR(2,1)+p%R5(7)*dr2*KnR(2,2)
      R6=p%R6(1)*dr2+p%R6(2)*logr+p%R6(3)+p%R6(4)*dr1*KnR(1,2)+p%R6(5)*dr2*KnR(2,1)+p%R6(6)*dr2*KnR(2,2)
      ! Add integration points
      do il=1,2
        do ik=1,2
          ! u*_{lk}
          fs_u=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
          ! t*_{lk}
          fs_t=T1*drdx(il)*drdx(ik)*drdn+T2*(drdn*c_dkr(il,ik)+drdx(ik)*n(il))+T3*drdx(il)*n(ik)
          ! sigma*_{lkj}·t_j
          fs_tt=T1*drdx(il)*drdx(ik)*drdt+T2*(drdt*c_dkr(il,ik)+drdx(ik)*t(il))+T3*drdx(il)*t(ik)
          do im=1,2
            ! u*_{lk,m}
            fs_dudm=U1*drdx(im)*c_dkr(il,ik)+U2*drdx(il)*drdx(ik)*drdx(im)+U3*(c_dkr(il,im)*drdx(ik)+c_dkr(ik,im)*drdx(il))
            ! sigma*_{lkj,m}·n_{j}
            fs_dsigdm=R1*(c_dkr(il,im)*drdx(ik)*drdn+c_dkr(ik,im)*drdx(il)*drdn+drdx(il)*drdx(ik)*n(im))&
                     +R2*(c_dkr(il,ik)*drdx(im)*drdn+n(il)*drdx(ik)*drdx(im))+R3*drdx(il)*drdx(ik)*drdx(im)*drdn&
                     +R4*c_dkr(il,im)*n(ik)+R5*n(ik)*drdx(il)*drdx(im)+R6*(c_dkr(il,ik)*n(im)+c_dkr(ik,im)*n(il))
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
    h1=p%cte_t*h1
    h2=p%cte_t*h2
    g1=p%cte_u*g1
    g2=p%cte_u*g2
    ! Assemble depending on reversion or not
    if (reverse) then
      h1=-h1
      h2=-h2
    end if
  end subroutine fbem_bem_harela2d_vsbie_ext_st

  !! Adaptive calculation of h1, h2, g1 and g2 when the collocation point is outside the integration element.
  recursive subroutine fbem_bem_harela2d_vsbie_ext_adp(e,reverse,xi_s,x_i,p,qsp,ks,ns,h1,h2,g1,g2)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                               !! Element
    logical                            :: reverse                         !! Reverse orientation
    real(kind=real64)                  :: xi_s(1,2)                       !! Subdivision of the parent element
    real(kind=real64)                  :: x_i(2)                          !! Collocation point position vector
    type(fbem_bem_harela2d_parameters) :: p                               !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp                             !! Quasi-singular integration parameters
    integer                            :: ks                              !! Current level of subdivisions
    integer                            :: ns                              !! Maximum level of subdivision
    complex(kind=real64)               :: h1(e%dme_n_gnodes,e%n_pnodes,2,2,2) !! h1 matrix
    complex(kind=real64)               :: h2(               e%n_pnodes,2,2,2) !! h2 matrix
    complex(kind=real64)               :: g1(e%dme_n_gnodes,e%n_snodes,2,2,2) !! g1 matrix
    complex(kind=real64)               :: g2(               e%n_snodes,2,2,2) !! g2 matrix
    ! Local
    integer              :: gln_near                            ! 1D Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer              :: gln                                 ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical              :: subdivide                           ! True if subdivision has to be performed
    real(kind=real64)    :: barxi(1)                            ! Nearest element coordinate with respect to collocation point
    real(kind=real64)    :: barxip(1)                           ! Nearest element subdivision local coordinate with respect to collocation point
    real(kind=real64)    :: rmin                                ! Minimum distance between collocation point and barxi on the element
    real(kind=real64)    :: barr                                ! Telles jacobian at barxi
    real(kind=real64)    :: cl                                  ! Characteristic length
    real(kind=real64)    :: d                                   ! Normalized distance between collocation point and element subdivision
    integer              :: method                              ! Method used in nearest point algorithm
    real(kind=real64)    :: tmp_xi_s(1,2)                       ! Subdivision
    real(kind=real64)    :: x_s(2,e%n_gnodes)                   ! Coordinates of the element subdivision
    complex(kind=real64) :: h1_tmp(e%dme_n_gnodes,e%n_pnodes,2,2,2) ! h1 matrix (temporary)
    complex(kind=real64) :: h2_tmp(               e%n_pnodes,2,2,2) ! h2 matrix (temporary)
    complex(kind=real64) :: g1_tmp(e%dme_n_gnodes,e%n_snodes,2,2,2) ! g1 matrix (temporary)
    complex(kind=real64) :: g2_tmp(               e%n_snodes,2,2,2) ! g2 matrix (temporary)
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
    gln_near=fbem_qs_n_estimation_telles(e%n,e%gtype,6,qsp,d,barxip)
    ! Decide if subdivide or calculate the subdivision
    subdivide=.false.
    if (ks.eq.ns) then
      if (gln_near.eq.0) then
        call fbem_warning_message(error_unit,0,'fbem_bem_harela2d_vsbie_ext_adp',ns,'maximum number of subdivisions reached')
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
      call fbem_bem_harela2d_vsbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h1,h2,g1,g2)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_harela2d_vsbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h1,h2,g1,g2)
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_harela2d_vsbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,p,gln,h1_tmp,h2_tmp,g1_tmp,g2_tmp)
      h1=h1+h1_tmp
      h2=h2+h2_tmp
      g1=g1+g1_tmp
      g2=g2+g2_tmp
    end if
  end subroutine fbem_bem_harela2d_vsbie_ext_adp

  !! Calculation of h1 and g1 when the collocation point is on the integration element.
  subroutine fbem_bem_harela2d_vsbie_int(e,reverse,xi_i,p,h1,g1)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                                   !! Integration element
    logical                            :: reverse                             !! Reverse normal vector
    real(kind=real64)                  :: xi_i(1)                             !! Local coordinate of the singular point.
    type(fbem_bem_harela2d_parameters) :: p                                   !! Parameters of the region
    complex(kind=real64)               :: h1(e%dme_n_gnodes,e%n_pnodes,2,2,2) !! h1 matrix
    complex(kind=real64)               :: g1(e%dme_n_gnodes,e%n_snodes,2,2,2) !! g1 matrix
    ! Local
    integer              :: gln                              ! 1D Gauss-Legendre number of integration points (<=32)
    integer              :: il, ik, im, iq                   ! Counter variables
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
    real(kind=real64)    :: xip(1)                           ! Coordinate  xip of subdivided element [0,1]
    real(kind=real64)    :: js                               ! Jacobian of the xi [xisub(1,:),xisub[2,:]] -> xip [0,1] transformation
    real(kind=real64)    :: xip_i(2)                         ! Singular point in xip space
    real(kind=real64)    :: xin(1)                           ! Coordinate xi
    real(kind=real64)    :: xisub(2,2)                       ! Coordinates of element subdivisions
    real(kind=real64)    :: aux(10)                          ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)    :: x(2)                             ! Position vector at xi
    real(kind=real64)    :: T(2)                             ! Tangent vector at xi
    real(kind=real64)    :: N(2)                             ! Normal vector at xi
    real(kind=real64)    :: t_i(2), n_i(2)                   ! Unit tangent and normal at xi_i
    real(kind=real64)    :: rv(2)                            ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, dr1, dr2, logr                ! Distance vector module, its inverse and log(1/r)
    real(kind=real64)    :: ra, rb                           ! Distance vector from collocation point to element vertices
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
    complex(kind=real64) :: z(2)                             ! Wavenumbers
    complex(kind=real64) :: KnR(0:2,2)                       ! Bessel functions decomposition
    complex(kind=real64) :: psi, chi, T1, T2, T3             ! Fundamental solution components
    complex(kind=real64) :: U1, U2, U3                       ! Fundamental solution components
    complex(kind=real64) :: R1, R2, R3, R4, R5, R6           ! Fundamental solution components
    complex(kind=real64) :: fs_u, fs_t                       ! Fundamental solutions values
    complex(kind=real64) :: fs_dudm, fs_tt, fs_dsigdm        ! Fundamental solutions values
    complex(kind=real64) :: fs_tt_s                          ! Fundamental solutions values
    real(kind=real64)    :: Idr1dr                           ! Int 1/r dr
    complex(kind=real64) :: h(e%n_pnodes,2,2)                ! h matrix
    complex(kind=real64) :: g(e%n_snodes,2,2)                ! g matrix
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
        xin=js*xip+xisub(1,ksub)
        ! dr/dGamma
        if (xin(1).lt.xi_i(1)) then
          drdti=-1.d0
        else
          drdti= 1.d0
        end if
        ! XI->X TRANSFORMATION
#       define etype e%gtype
#       define delta 0.d0
#       define xi xin(1)
#       define phi gphi
#       define dphidxi dgphidxi
#       include <phi_and_dphidxi_1d.rc>
#       undef etype
#       undef delta
#       undef xi
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
#       define xi xin(1)
#       define phi pphi
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef xi
#       undef phi
        ! Functional shape functions (secondary variables) at xi
#       define etype e%stype
#       define delta e%stype_delta
#       define xi xin(1)
#       define phi sphi
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef xi
#       undef phi
        ! Calculation of the components of the fundamental solution
        z(1)=c_im*p%k1*r
        z(2)=c_im*p%k2*r
        call fbem_BesselKnR_decomposed(2,z,KnR)
        psi=p%psi(1)*logr+p%psi(2)+KnR(0,2)+dr1*p%psi(3)*KnR(1,1)+dr1*p%psi(4)*KnR(1,2)
        chi=p%chi(1)+(p%chi(2)*logr+p%chi(3))*r**2+p%chi(4)*KnR(2,1)+KnR(2,2)
        T1=p%T1(2)*r+p%T1(3)*KnR(1,1)+p%T1(4)*KnR(1,2)+p%T1(5)*dr1*KnR(2,1)+p%T1(6)*dr1*KnR(2,2)
        T2=(p%T2(2)*logr+p%T2(3))*r+p%T2(4)*KnR(1,2)+p%T2(5)*dr1*KnR(2,1)+p%T2(6)*dr1*KnR(2,2)
        T3=(p%T3(2)*logr+p%T3(3))*r+p%T3(4)*KnR(1,1)+p%T3(5)*dr1*KnR(2,1)+p%T3(6)*dr1*KnR(2,2)
        U1=p%U1(1)*dr1+(p%U1(2)*logr+p%U1(3))*r+p%U1(4)*KnR(1,2)+p%U1(5)*dr1*KnR(2,1)+p%U1(6)*dr1*KnR(2,2)
        U2=p%U2(1)*dr1+p%U2(2)*r+p%U2(3)*KnR(1,1)+p%U2(4)*KnR(1,2)+p%U2(5)*dr1*KnR(2,1)+p%U2(6)*dr1*KnR(2,2)
        U3=p%U3(1)*dr1+(p%U3(2)*logr+p%U3(3))*r+p%U3(4)*dr1*KnR(2,1)+p%U3(5)*dr1*KnR(2,2)
        R1=p%R1(2)+p%R1(3)*dr1*KnR(1,1)+p%R1(4)*dr1*KnR(1,2)+p%R1(5)*dr2*KnR(2,1)+p%R1(6)*dr2*KnR(2,2)
        R2=p%R2(2)+p%R2(3)*KnR(0,2)+p%R2(4)*dr1*KnR(1,1)+p%R2(5)*dr1*KnR(1,2)+p%R2(6)*dr2*KnR(2,1)+p%R2(7)*dr2*KnR(2,2)
        R3=p%R3(2)+p%R3(3)*KnR(0,1)+p%R3(4)*KnR(0,2)+p%R3(5)*dr1*KnR(1,1)+p%R3(6)*dr1*KnR(1,2)+p%R3(7)*dr2*KnR(2,1)&
          +p%R3(8)*dr2*KnR(2,2)
        R4=p%R4(2)*logr+p%R4(3)+p%R4(4)*dr1*KnR(1,1)+p%R4(5)*dr2*KnR(2,1)+p%R4(6)*dr2*KnR(2,2)
        R5=p%R5(2)+p%R5(3)*KnR(0,1)+p%R5(4)*dr1*KnR(1,1)+p%R5(5)*dr1*KnR(1,2)+p%R5(6)*dr2*KnR(2,1)+p%R5(7)*dr2*KnR(2,2)
        R6=p%R6(2)*logr+p%R6(3)+p%R6(4)*dr1*KnR(1,2)+p%R6(5)*dr2*KnR(2,1)+p%R6(6)*dr2*KnR(2,2)
        ! Add integration points
        do il=1,2
          do ik=1,2
            ! u*_{lk} (all parts)
            fs_u=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
            ! t*_{lk} (all parts)
            fs_t=(p%T1(1)*dr1+T1)*drdx(il)*drdx(ik)*drdn+(p%T2(1)*dr1+T2)*(drdn*c_dkr(il,ik)+drdx(ik)*n(il))&
                +(p%T3(1)*dr1+T3)*drdx(il)*n(ik)
            ! sigma*_{lkj}·t_j (regular part)
            fs_tt=T1*drdx(il)*drdx(ik)*drdt+T2*(drdt*c_dkr(il,ik)+drdx(ik)*t(il))+T3*drdx(il)*t(ik)
            ! sigma*_{lkj}·t_j (singular part)
            fs_tt_s=dr1*(p%T1(1)*drdx(il)*drdx(ik)*drdt+p%T2(1)*(drdt*c_dkr(il,ik)+drdx(ik)*t(il))+p%T3(1)*drdx(il)*t(ik))
            do im=1,2
              ! u*_{lk,m} (all parts)
              fs_dudm=U1*drdx(im)*c_dkr(il,ik)+U2*drdx(il)*drdx(ik)*drdx(im)+U3*(c_dkr(il,im)*drdx(ik)+c_dkr(ik,im)*drdx(il))
              ! sigma*_{lkj,m}·n_{j} (regular part)
              fs_dsigdm=R1*(c_dkr(il,im)*drdx(ik)*drdn+c_dkr(ik,im)*drdx(il)*drdn+drdx(il)*drdx(ik)*n(im))&
                       +R2*(c_dkr(il,ik)*drdx(im)*drdn+n(il)*drdx(ik)*drdx(im))+R3*drdx(il)*drdx(ik)*drdx(im)*drdn&
                       +R4*c_dkr(il,im)*n(ik)+R5*n(ik)*drdx(il)*drdx(im)+R6*(c_dkr(il,ik)*n(im)+c_dkr(ik,im)*n(il))
              do iq=1,e%dme_n_gnodes
                ! H^r
                h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)+fs_dsigdm*(dme_gphi(iq)-dme_gphi_i(iq))*pphi(:)*jw
                h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)+dr2*drdn*(dme_gphi(iq)-dme_gphi_i(iq))*(p%R3(1)*drdx(il)*drdx(ik)*drdx(im)&
                                                             +p%R2(1)*c_dkr(il,ik)*drdx(im)&
                                                             +p%R1(1)*(c_dkr(il,im)*drdx(ik)+c_dkr(ik,im)*drdx(il)))*pphi(:)*jw
                h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)+dr2*(dme_gphi(iq)-dme_gphi_i(iq)-dot_product(dme_wqj_i(iq,:),rv))*(&
                                                    p%R2(1)*drdx(im)*(n(il)*drdx(ik)-n(ik)*drdx(il))&
                                                   +p%R1(1)*n(im)*drdx(il)*drdx(ik)&
                                                   +p%R6(1)*(c_dkr(il,ik)*n(im)+c_dkr(ik,im)*n(il)-c_dkr(il,im)*n(ik)))*pphi(:)*jw
                h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)+dr1*(dot_product(dme_wqj_i(iq,:),drdx)*pphi(:)-drdti*dme_vq_i(iq)*pphi_i(:))*(&
                                                    p%R2(1)*drdx(im)*(n(il)*drdx(ik)-n(ik)*drdx(il))&
                                                   +p%R1(1)*n(im)*drdx(il)*drdx(ik)&
                                                   +p%R6(1)*(c_dkr(il,ik)*n(im)+c_dkr(ik,im)*n(il)-c_dkr(il,im)*n(ik)))*jw
                h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)+dme_vq_i(iq)*pphi_i(:)*dr1*drdti*(&
                                                    p%R2(1)*(drdx(im)-drdti*t_i(im))*(n(il)*drdx(ik)-n(ik)*drdx(il))&
                                                   +p%R1(1)*(n(im)*drdx(il)*drdx(ik)-n_i(im)*t_i(il)*t_i(ik))&
                                                   +p%R6(1)*(c_dkr(il,ik)*(n(im)-n_i(im))+c_dkr(ik,im)*(n(il)-n_i(il))&
                                                   -c_dkr(il,im)*(n(ik)-n_i(ik))))*jw
                h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)+dme_vq_i(iq)*pphi_i(:)*(p%R1(1)*n_i(im)*t_i(il)*t_i(ik)+p%R6(1)*(&
                                                    c_dkr(il,ik)*n_i(im)+c_dkr(ik,im)*n_i(il)&
                                                   -c_dkr(il,im)*n_i(ik)))*dr1*(drdti-drdt)*jw
                ! H^n
                h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)-fs_tt*dme_vq(iq)*n(im)*pphi(:)*jw
                h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)-fs_tt_s*(dme_vq(iq)*n(im)*pphi(:)-dme_vq_i(iq)*n_i(im)*pphi_i(:))*jw
                h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)-dme_vq_i(iq)*n_i(im)*pphi_i(:)*dr1*(&
                                                    drdt*p%T1(1)*(drdx(il)*drdx(ik)-t_i(il)*t_i(ik))&
                                                   +p%T2(1)*(t(il)*drdx(ik)-t(ik)*drdx(il)))*jw
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
    h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)+p%R2(1)*dme_vq_i(iq)*pphi_i(:)*t_i(im)*(1.d0-c_dkr(il,ik))*(c_dkr(il,1)-c_dkr(il,2))*Idr1dr
    h1(iq,:,il,ik,im)=h1(iq,:,il,ik,im)+dme_vq_i(iq)*pphi_i(:)*p%R6(1)*(c_dkr(ik,im)*n_i(il)-c_dkr(il,im)*n_i(ik))*Idr1dr
    end do
    end do
    end do
    end do
    ! Multiply by constants
    h1=p%cte_t*h1
    g1=p%cte_u*g1
    ! Assemble depending on reversion or not
    if (reverse) h1=-h1
    ! Assemble matrices from SBIE
    call fbem_bem_harela2d_sbie_int(gln,e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,xi_i(1),p,h,g)
    do im=1,2
      do iq=1,e%dme_n_gnodes
        h1(iq,:,:,:,im)=h1(iq,:,:,:,im)+dme_vq_i(iq)*t_i(im)*h(:,:,:)
        g1(iq,:,:,:,im)=g1(iq,:,:,:,im)+dme_vq_i(iq)*t_i(im)*g(:,:,:)
      end do
    end do
  end subroutine fbem_bem_harela2d_vsbie_int

  !! Efficient automatic calculation of h, g, m and l.
  subroutine fbem_bem_harela2d_vsbie_auto(e,reverse,x_i,p,qsp,ns,h1,h2,g1,g2)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                                   !! Integration element
    logical                            :: reverse                             !! Reverse orientation
    real(kind=real64)                  :: x_i(2)                              !! Collocation point
    type(fbem_bem_harela2d_parameters) :: p                                   !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp                                 !! Quasi-singular integration parameters
    integer                            :: ns                                  !! Maximum level of subdivisions
    complex(kind=real64)               :: h1(e%dme_n_gnodes,e%n_pnodes,2,2,2) !! h1 matrix
    complex(kind=real64)               :: h2(               e%n_pnodes,2,2,2) !! h2 matrix
    complex(kind=real64)               :: g1(e%dme_n_gnodes,e%n_snodes,2,2,2) !! g1 matrix
    complex(kind=real64)               :: g2(               e%n_snodes,2,2,2) !! g2 matrix
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
        h2=(0.d0,0.d0)
        g2=(0.d0,0.d0)
        call fbem_bem_harela2d_vsbie_int(e,reverse,barxi,p,h1,g1)
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
          call fbem_bem_harela2d_vsbie_ext_pre(ps,e,reverse,x_i,p,h1,h2,g1,g2)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_harela2d_vsbie_ext_adp(e,reverse,xi_s,x_i,p,qsp,1,ns,h1,h2,g1,g2)
        end if
    end select
  end subroutine fbem_bem_harela2d_vsbie_auto

end module fbem_bem_harela2d
