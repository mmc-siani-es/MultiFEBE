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
!! <b> This module implements the calculation of element-wise integrals of Boundary Integral Equations of the 3D Biot's
!! poroelastodynamics.</b>
module fbem_bem_harpor3d

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
  use fbem_bem_staela3d

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all are private
  private

  ! ================================================================================================================================
  ! INITIAL SETUP
  public :: fbem_bem_harpor3d_parameters
  public :: fbem_bem_harpor3d_calculate_parameters
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE)
  ! Exterior integration
  public :: fbem_bem_harpor3d_sbie_ext_pre
  public :: fbem_bem_harpor3d_sbie_ext_st
  public :: fbem_bem_harpor3d_sbie_ext_adp
  ! Interior integration
  public :: fbem_bem_harpor3d_sbie_int
  ! Automatic integration
  public :: fbem_bem_harpor3d_sbie_auto
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (HBIE)
  ! Exterior integration
  public :: fbem_bem_harpor3d_hbie_ext_pre
  public :: fbem_bem_harpor3d_hbie_ext_st
  public :: fbem_bem_harpor3d_hbie_ext_adp
  ! Interior integration
  public :: fbem_bem_harpor3d_hbie_int
  ! Automatic integration
  public :: fbem_bem_harpor3d_hbie_auto
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! BOTH SBIE AND HBIE
  ! Exterior integration
  public :: fbem_bem_harpor3d_shbie_ext_pre
  ! Automatic integration
  public :: fbem_bem_harpor3d_shbie_auto
  ! ================================================================================================================================

  ! ================================================================================================================================
  !! Data structure for region parameters (material properties and frequency dependant parameters)
  type fbem_bem_harpor3d_parameters
    ! Region parameters
    complex(kind=real64) :: lambda
    complex(kind=real64) :: mu
    complex(kind=real64) :: nu
    real(kind=real64)    :: rho1
    real(kind=real64)    :: rho2
    real(kind=real64)    :: rhoa
    complex(kind=real64) :: R
    complex(kind=real64) :: Q
    real(kind=real64)    :: b
    complex(kind=real64) :: rhohat11
    complex(kind=real64) :: rhohat22
    complex(kind=real64) :: rhohat12
    complex(kind=real64) :: Z
    complex(kind=real64) :: J
    ! Frequency
    real(kind=real64)    :: omega
    ! Wavenumbers
    complex(kind=real64) :: k1, k2, k3
    ! Wave propagation speeds
    complex(kind=real64) :: c1, c2, c3
    ! Coefficients of the components of the fundamental solution
    complex(kind=real64) :: eta(3)
    complex(kind=real64) :: vartheta(5)
    complex(kind=real64) :: psi(8), chi(9)
    complex(kind=real64) :: W0(6)
    complex(kind=real64) :: T01(7), T02(8)
    complex(kind=real64) :: W1(10), W2(9)
    complex(kind=real64) :: T1(14), T2(12), T3(13)
    complex(kind=real64) :: Q1(11), Q2(10)
    complex(kind=real64) :: S01(14), S02(13), S03(12)
    complex(kind=real64) :: S1(15), S2(16), S3(17), S4(13), S5(16)
    complex(kind=real64) :: cte_u(0:3,0:3), cte_t(0:3,0:3), cte_d(0:3,0:3), cte_s(0:3,0:3)
  end type fbem_bem_harpor3d_parameters
  ! ================================================================================================================================


contains

  ! ============================================================================================================================== !
  ! INITIAL SETUP                                                                                                                  !
  ! ============================================================================================================================== !

  !! Subroutine that calculate basic parameters for a given region properties and a given frequency
  subroutine fbem_bem_harpor3d_calculate_basic_parameters(lambda,mu,rho1,rho2,rhoa,R,Q,b,omega,pars)
    implicit none
    ! I/O
    complex(kind=real64)               :: lambda
    complex(kind=real64)               :: mu
    real(kind=real64)                  :: rho1
    real(kind=real64)                  :: rho2
    real(kind=real64)                  :: rhoa
    complex(kind=real64)               :: R
    complex(kind=real64)               :: Q
    real(kind=real64)                  :: b
    real(kind=real64)                  :: omega
    type(fbem_bem_harpor3d_parameters) :: pars
    ! Local
    complex(kind=real64)               :: rhohat11, rhohat22, rhohat12
    complex(kind=real64)               :: Z, J
    complex(kind=real64)               :: ca, cb, cc
    complex(kind=real64)               :: k1, k2, k3, tmpk
    !
    ! Frequency-dependant parameters
    !
    rhohat11=rho1+rhoa-c_im*b/omega
    rhohat12=    -rhoa+c_im*b/omega
    rhohat22=rho2+rhoa-c_im*b/omega
    J=1.0d0/((omega**2)*rhohat22)
    Z=rhohat12/rhohat22
    k3=zsqrt(1.0d0/(mu*J)*(rhohat11/rhohat22-Z**2))
    if (real(k3).lt.0.d0) k3=-k3
    ca=lambda+2.0d0*mu
    cb=(lambda+2.0d0*mu)/(J*R)+mu*k3**2+1.0d0/J*(Q/R-Z)**2
    cc=mu/(J*R)*k3**2
    k1=zsqrt(0.5d0*(cb-zsqrt(cb*cb-4.0d0*ca*cc))/ca)
    k2=zsqrt(0.5d0*(cb+zsqrt(cb*cb-4.0d0*ca*cc))/ca)
    if (real(k1).lt.0.d0) k1=-k1
    if (real(k2).lt.0.d0) k2=-k2
    if (dreal(k1).gt.dreal(k2)) then
      tmpk=k1
      k1=k2
      k2=tmpk
    end if
    !
    ! Save to structure
    !
    pars%omega=omega
    pars%lambda=lambda
    pars%mu=mu
    pars%nu=0.5d0*lambda/(lambda+mu)
    pars%rho1=rho1
    pars%rho2=rho2
    pars%rhoa=rhoa
    pars%R=R
    pars%Q=Q
    pars%b=b
    pars%rhohat11=rhohat11
    pars%rhohat22=rhohat22
    pars%rhohat12=rhohat12
    pars%Z=Z
    pars%J=J
    pars%k1=k1
    pars%k2=k2
    pars%k3=k3
    pars%c1=omega/k1
    pars%c2=omega/k2
    pars%c3=omega/k3
  end subroutine fbem_bem_harpor3d_calculate_basic_parameters

  !! Subroutine that calculate all parameters for a given region properties and a given frequency
  subroutine fbem_bem_harpor3d_calculate_parameters(lambda,mu,rho1,rho2,rhoa,R,Q,b,omega,pars)
    implicit none
    ! I/O
    complex(kind=real64)               :: lambda
    complex(kind=real64)               :: mu
    real(kind=real64)                  :: rho1
    real(kind=real64)                  :: rho2
    real(kind=real64)                  :: rhoa
    complex(kind=real64)               :: R
    complex(kind=real64)               :: Q
    real(kind=real64)                  :: b
    real(kind=real64)                  :: omega
    type(fbem_bem_harpor3d_parameters) :: pars
    ! Local
    complex(kind=real64)               :: rhohat11, rhohat22, rhohat12
    complex(kind=real64)               :: Z, J
    complex(kind=real64)               :: ca, cb, cc
    complex(kind=real64)               :: k1, k2, k3, tmpk
    complex(kind=real64)               :: alpha1, alpha2, beta1, beta2, vartheta_c, k1mk2
    !
    ! Frequency-dependant parameters
    !
    rhohat11=rho1+rhoa-c_im*b/omega
    rhohat12=    -rhoa+c_im*b/omega
    rhohat22=rho2+rhoa-c_im*b/omega
    J=1.0d0/(rhohat22*omega**2)
    Z=rhohat12/rhohat22
    k3=zsqrt(1.0d0/(mu*J)*(rhohat11/rhohat22-Z**2))
    if (real(k3).lt.0.d0) k3=-k3
    ca=lambda+2.0d0*mu
    cb=(lambda+2.0d0*mu)/(J*R)+mu*k3**2+1.0d0/J*(Q/R-Z)**2
    cc=mu/(J*R)*k3**2
    k1=zsqrt(0.5d0*(cb-zsqrt(cb*cb-4.0d0*ca*cc))/ca)
    k2=zsqrt(0.5d0*(cb+zsqrt(cb*cb-4.0d0*ca*cc))/ca)
    if (real(k1).lt.0.d0) k1=-k1
    if (real(k2).lt.0.d0) k2=-k2
    if (dreal(k1).gt.dreal(k2)) then
      tmpk=k1
      k1=k2
      k2=tmpk
    end if
    !
    ! Save to structure
    !
    pars%omega=omega
    pars%lambda=lambda
    pars%mu=mu
    pars%nu=0.5d0*lambda/(lambda+mu)
    pars%rho1=rho1
    pars%rho2=rho2
    pars%rhoa=rhoa
    pars%R=R
    pars%Q=Q
    pars%b=b
    pars%rhohat11=rhohat11
    pars%rhohat22=rhohat22
    pars%rhohat12=rhohat12
    pars%Z=Z
    pars%J=J
    pars%k1=k1
    pars%k2=k2
    pars%k3=k3
    pars%c1=omega/k1
    pars%c2=omega/k2
    pars%c3=omega/k3
    !
    ! Coefficients of equations of fundamental solution's components
    !
    ! Constants
    alpha1=k1**2-mu/(lambda+2.0*mu)*k3**2
    alpha2=k2**2-mu/(lambda+2.0*mu)*k3**2
    beta1=mu/(lambda+2.0*mu)*k1**2-k1**2*k2**2/k3**2
    beta2=mu/(lambda+2.0*mu)*k2**2-k1**2*k2**2/k3**2
    vartheta_c=(Q/R-Z)/(lambda+2.0d0*mu)
    k1mk2=k1**2-k2**2
    ! eta
    pars%eta(1)=-(c_im*k1*alpha1-c_im*k2*alpha2)/k1mk2
    pars%eta(2)= alpha1/k1mk2
    pars%eta(3)=-alpha2/k1mk2
    ! vartheta
    pars%vartheta(1)= 0.5d0*vartheta_c
    pars%vartheta(2)= vartheta_c*c_im*k1/k1mk2
    pars%vartheta(3)=-vartheta_c*c_im*k2/k1mk2
    pars%vartheta(4)= vartheta_c/k1mk2
    pars%vartheta(5)=-vartheta_c/k1mk2
    ! psi
    pars%psi(1)= 0.5d0*(lambda+3.d0*mu)/(lambda+2.d0*mu)
    pars%psi(2)=-1.d0/3.d0*((c_im*k1*beta1-c_im*k2*beta2)/k1mk2+2.d0*c_im*k3)
    pars%psi(3)=-beta1/(c_im*k1)/k1mk2
    pars%psi(4)= beta2/(c_im*k2)/k1mk2
    pars%psi(5)= 1.d0/(c_im*k3)
    pars%psi(6)= beta1/k1**2/k1mk2
    pars%psi(7)=-beta2/k2**2/k1mk2
    pars%psi(8)=-1.d0/k3**2
    ! chi
    pars%chi(1)=-0.5d0*(lambda+mu)/(lambda+2.d0*mu)
    pars%chi(2)=-beta1/k1mk2
    pars%chi(3)= beta2/k1mk2
    pars%chi(4)=-3.d0*beta1/(c_im*k1)/k1mk2
    pars%chi(5)= 3.d0*beta2/(c_im*k2)/k1mk2
    pars%chi(6)= 3.d0/(c_im*k3)
    pars%chi(7)= 3.d0*beta1/k1**2/k1mk2
    pars%chi(8)=-3.d0*beta2/k2**2/k1mk2
    pars%chi(9)=-3.d0/k3**2
    ! W0
    pars%W0(1)= J
    pars%W0(2)= 0.5d0*(Z*vartheta_c+J*(k1**2*alpha1-k2**2*alpha2)/k1mk2)
    pars%W0(3)= (Z*vartheta_c+J*alpha1)*c_im*k1/k1mk2
    pars%W0(4)=-(Z*vartheta_c+J*alpha2)*c_im*k2/k1mk2
    pars%W0(5)= (Z*vartheta_c+J*alpha1)/k1mk2
    pars%W0(6)=-(Z*vartheta_c+J*alpha2)/k1mk2
    ! T01
    pars%T01(1)= mu*vartheta_c
    pars%T01(2)=-2.d0*mu*vartheta_c*k1**2/k1mk2
    pars%T01(3)= 2.d0*mu*vartheta_c*k2**2/k1mk2
    pars%T01(4)= 6.d0*mu*vartheta_c*c_im*k1/k1mk2
    pars%T01(5)=-6.d0*mu*vartheta_c*c_im*k2/k1mk2
    pars%T01(6)= 6.d0*mu*vartheta_c/k1mk2
    pars%T01(7)=-6.d0*mu*vartheta_c/k1mk2
    ! T02
    pars%T02(1)= mu*vartheta_c+Z
    pars%T02(2)= (lambda+2.d0/3.d0*mu)*vartheta_c*c_im*(k1**3-k2**3)/k1mk2-Q/R*(c_im*k1*alpha1-c_im*k2*alpha2)/k1mk2
    pars%T02(3)= (Q/R*alpha1-lambda*vartheta_c*k1**2)/k1mk2
    pars%T02(4)=-(Q/R*alpha2-lambda*vartheta_c*k2**2)/k1mk2
    pars%T02(5)=-2.d0*mu*vartheta_c*c_im*k1/k1mk2
    pars%T02(6)= 2.d0*mu*vartheta_c*c_im*k2/k1mk2
    pars%T02(7)=-2.d0*mu*vartheta_c/k1mk2
    pars%T02(8)= 2.d0*mu*vartheta_c/k1mk2
    ! W1
    pars%W1( 1)= 0.5d0*(Q/R*mu/(lambda+2.d0*mu)-Z)
    pars%W1( 2)=-(mu*vartheta_c*k1**2+Z*beta1)/k1mk2
    pars%W1( 3)= (mu*vartheta_c*k2**2+Z*beta2)/k1mk2
    pars%W1( 4)= Z
    pars%W1( 5)= 3.d0*(mu*vartheta_c*c_im*k1-Z*beta1/(c_im*k1))/k1mk2
    pars%W1( 6)=-3.d0*(mu*vartheta_c*c_im*k2-Z*beta2/(c_im*k2))/k1mk2
    pars%W1( 7)= 3.d0*Z/(c_im*k3)
    pars%W1( 8)= 3.d0*(mu*vartheta_c+Z*beta1/k1**2)/k1mk2
    pars%W1( 9)=-3.d0*(mu*vartheta_c+Z*beta2/k2**2)/k1mk2
    pars%W1(10)=-3.d0*Z/k3**2
    ! W2
    pars%W2(1)=-0.5d0*(Q/R*mu/(lambda+2.d0*mu)+Z)
    pars%W2(2)= 1.d0/3.d0*(mu*vartheta_c*c_im*(k1**3-k2**3)/k1mk2+Z*((c_im*k1*beta1-c_im*k2*beta2)/k1mk2+2.d0*c_im*k3))
    pars%W2(3)=-Z
    pars%W2(4)=-(mu*vartheta_c*c_im*k1-Z*beta1/(c_im*k1))/k1mk2
    pars%W2(5)= (mu*vartheta_c*c_im*k2-Z*beta2/(c_im*k2))/k1mk2
    pars%W2(6)=-Z/(c_im*k3)
    pars%W2(7)=-(mu*vartheta_c+Z*beta1/k1**2)/k1mk2
    pars%W2(8)= (mu*vartheta_c+Z*beta2/k2**2)/k1mk2
    pars%W2(9)= Z/k3**2
    ! T1
    pars%T1( 1)=-3.d0*(lambda+mu)/(lambda+2.d0*mu)
    pars%T1( 2)= 0.25d0*((k1**2*beta1-k2**2*beta2)/k1mk2-k3**2)
    pars%T1( 3)=-2.d0*c_im*k1*beta1/k1mk2
    pars%T1( 4)= 2.d0*c_im*k2*beta2/k1mk2
    pars%T1( 5)= 2.d0*c_im*k3
    pars%T1( 6)=-12.d0*beta1/k1mk2
    pars%T1( 7)= 12.d0*beta2/k1mk2
    pars%T1( 8)= 12.d0
    pars%T1( 9)=-30.d0*beta1/(c_im*k1)/k1mk2
    pars%T1(10)= 30.d0*beta2/(c_im*k2)/k1mk2
    pars%T1(11)= 30.d0/(c_im*k3)
    pars%T1(12)= 30.d0*beta1/k1**2/k1mk2
    pars%T1(13)=-30.d0*beta2/k2**2/k1mk2
    pars%T1(14)=-30.d0/k3**2
    ! T2
    pars%T2( 1)=-mu/(lambda+2.d0*mu)
    pars%T2( 2)=-0.25d0*((k1**2*beta1-k2**2*beta2)/k1mk2+k3**2)
    pars%T2( 3)=-c_im*k3
    pars%T2( 4)= 2.d0*beta1/k1mk2
    pars%T2( 5)=-2.d0*beta2/k1mk2
    pars%T2( 6)=-3.d0
    pars%T2( 7)= 6.d0*beta1/(c_im*k1)/k1mk2
    pars%T2( 8)=-6.d0*beta2/(c_im*k2)/k1mk2
    pars%T2( 9)=-6.d0/(c_im*k3)
    pars%T2(10)=-6.d0*beta1/k1**2/k1mk2
    pars%T2(11)= 6.d0*beta2/k2**2/k1mk2
    pars%T2(12)= 6.d0/k3**2
    ! T3
    pars%T3( 1)= mu/(lambda+2.d0*mu)
    pars%T3( 2)= 0.25d0*(2.d0*Q/R*vartheta_c/J-(2.d0*lambda/mu+1.d0)*(k1**2*beta1-k2**2*beta2)/k1mk2+k3**2)
    pars%T3( 3)= (Q/R*vartheta_c/J-lambda/mu*beta1)*c_im*k1/k1mk2
    pars%T3( 4)=-(Q/R*vartheta_c/J-lambda/mu*beta2)*c_im*k2/k1mk2
    pars%T3( 5)= (Q/R*vartheta_c/J-(lambda/mu-2.d0)*beta1)/k1mk2
    pars%T3( 6)=-(Q/R*vartheta_c/J-(lambda/mu-2.d0)*beta2)/k1mk2
    pars%T3( 7)=-2.d0
    pars%T3( 8)= 6.d0*beta1/(c_im*k1)/k1mk2
    pars%T3( 9)=-6.d0*beta2/(c_im*k2)/k1mk2
    pars%T3(10)=-6.d0/(c_im*k3)
    pars%T3(11)=-6.d0*beta1/k1**2/k1mk2
    pars%T3(12)= 6.d0*beta2/k2**2/k1mk2
    pars%T3(13)= 6.d0/k3**2
    ! Q1
    pars%Q1( 1)= 3.d0*J
    pars%Q1( 2)= Z*vartheta_c+0.5d0*J*(k1**2*alpha1-k2**2*alpha2)/k1mk2-0.5d0*Z**2/mu*(lambda+mu)/(lambda+2.d0*mu)
    pars%Q1( 3)=-((2.d0*Z*vartheta_c+J*alpha1)*k1**2+Z**2/mu*beta1)/k1mk2
    pars%Q1( 4)= ((2.d0*Z*vartheta_c+J*alpha2)*k2**2+Z**2/mu*beta2)/k1mk2
    pars%Q1( 5)= Z**2/mu
    pars%Q1( 6)= 3.d0*((2.d0*Z*vartheta_c+J*alpha1)*c_im*k1-Z**2/mu*beta1/(c_im*k1))/k1mk2
    pars%Q1( 7)=-3.d0*((2.d0*Z*vartheta_c+J*alpha2)*c_im*k2-Z**2/mu*beta2/(c_im*k2))/k1mk2
    pars%Q1( 8)= 3.d0*Z**2/mu/(c_im*k3)
    pars%Q1( 9)= 3.d0*((2.d0*Z*vartheta_c+J*alpha1)+Z**2/mu*beta1/k1**2)/k1mk2
    pars%Q1(10)=-3.d0*((2.d0*Z*vartheta_c+J*alpha2)+Z**2/mu*beta2/k2**2)/k1mk2
    pars%Q1(11)=-3.d0*Z**2/mu/k3**2
    ! Q2
    pars%Q2( 1)= J
    pars%Q2( 2)= Z*vartheta_c+0.5d0*J*(k1**2*alpha1-k2**2*alpha2)/k1mk2+0.5d0*Z**2/mu*(lambda+3.d0*mu)/(lambda+2.d0*mu)
    pars%Q2( 3)=-1.d0/3.d0*(2.d0*Z*vartheta_c*c_im*(k1**3-k2**3)/k1mk2+J*c_im*(k1**3*alpha1-k2**3*alpha2)/k1mk2&
                            +Z**2/mu*((c_im*k1*beta1-c_im*k2*beta2)/k1mk2+2.d0*c_im*k3))
    pars%Q2( 4)= Z**2/mu
    pars%Q2( 5)= ((2.d0*Z*vartheta_c+J*alpha1)*c_im*k1-Z**2/mu*beta1/(c_im*k1))/k1mk2
    pars%Q2( 6)=-((2.d0*Z*vartheta_c+J*alpha2)*c_im*k2-Z**2/mu*beta2/(c_im*k2))/k1mk2
    pars%Q2( 7)= Z**2/mu/(c_im*k3)
    pars%Q2( 8)= (2.d0*Z*vartheta_c+J*alpha1+Z**2/mu*beta1/k1**2)/k1mk2
    pars%Q2( 9)=-(2.d0*Z*vartheta_c+J*alpha2+Z**2/mu*beta2/k2**2)/k1mk2
    pars%Q2(10)=-Z**2/mu/k3**2
    ! S01
    pars%S01( 1)= 3.d0*(Q/R*mu/(lambda+2.d0*mu)-Z)
    pars%S01( 2)= 0.25d0*(mu*vartheta_c*(k1**4-k2**4)/k1mk2+Z*((k1**2*beta1-k2**2*beta2)/k1mk2-k3**2))
    pars%S01( 3)=-2.d0*(mu*vartheta_c*k1**2+Z*beta1)*c_im*k1/k1mk2
    pars%S01( 4)= 2.d0*(mu*vartheta_c*k2**2+Z*beta2)*c_im*k2/k1mk2
    pars%S01( 5)= 2.d0*Z*c_im*k3
    pars%S01( 6)=-12.d0*(mu*vartheta_c*k1**2+Z*beta1)/k1mk2
    pars%S01( 7)= 12.d0*(mu*vartheta_c*k2**2+Z*beta2)/k1mk2
    pars%S01( 8)= 12.d0*Z
    pars%S01( 9)= 30.d0*(mu*vartheta_c*c_im*k1-Z*beta1/(c_im*k1))/k1mk2
    pars%S01(10)=-30.d0*(mu*vartheta_c*c_im*k2-Z*beta2/(c_im*k2))/k1mk2
    pars%S01(11)= 30.d0*Z/(c_im*k3)
    pars%S01(12)= 30.d0*(mu*vartheta_c+Z*beta1/k1**2)/k1mk2
    pars%S01(13)=-30.d0*(mu*vartheta_c+Z*beta2/k2**2)/k1mk2
    pars%S01(14)=-30.d0*Z/k3**2
    ! S02
    pars%S02( 1)= Q/R*mu/(lambda+2.d0*mu)+Z
    pars%S02( 2)= 0.25d0*(2.d0*Q/R*Z/J*vartheta_c-(2.d0*lambda+mu)*vartheta_c*(k1**4-k2**4)/k1mk2&
                          +2.d0*Q/R*(k1**2*alpha1-k2**2*alpha2)/k1mk2&
                          -Z*((2.d0*lambda/mu+1.d0)*(k1**2*beta1-k2**2*beta2)/k1mk2-k3**2))
    pars%S02( 3)= (Q/R*Z/J*vartheta_c-lambda*vartheta_c*k1**2+Q/R*alpha1-Z*lambda/mu*beta1)*c_im*k1/k1mk2
    pars%S02( 4)=-(Q/R*Z/J*vartheta_c-lambda*vartheta_c*k2**2+Q/R*alpha2-Z*lambda/mu*beta2)*c_im*k2/k1mk2
    pars%S02( 5)= (Q/R*Z/J*vartheta_c-(lambda-2.d0*mu)*vartheta_c*k1**2+Q/R*alpha1+Z*(2.d0-lambda/mu)*beta1)/k1mk2
    pars%S02( 6)=-(Q/R*Z/J*vartheta_c-(lambda-2.d0*mu)*vartheta_c*k2**2+Q/R*alpha2+Z*(2.d0-lambda/mu)*beta2)/k1mk2
    pars%S02( 7)=-2.d0*Z
    pars%S02( 8)=-6.d0*(mu*vartheta_c+Z*beta1/k1**2)*c_im*k1/k1mk2
    pars%S02( 9)= 6.d0*(mu*vartheta_c+Z*beta2/k2**2)*c_im*k2/k1mk2
    pars%S02(10)=-6.d0*Z/(c_im*k3)
    pars%S02(11)=-6.d0*(mu*vartheta_c+Z*beta1/k1**2)/k1mk2
    pars%S02(12)= 6.d0*(mu*vartheta_c+Z*beta2/k2**2)/k1mk2
    pars%S02(13)= 6.d0*Z/k3**2
    ! S03
    pars%S03( 1)= Q/R*mu/(lambda+2.d0*mu)
    pars%S03( 2)= 0.25d0*(mu*vartheta_c*(k1**4-k2**4)/k1mk2+Z*((k1**2*beta1-k2**2*beta2)/k1mk2+k3**2))
    pars%S03( 3)= Z*c_im*k3
    pars%S03( 4)=-2.d0*(mu*vartheta_c*k1**2+Z*beta1)/k1mk2
    pars%S03( 5)= 2.d0*(mu*vartheta_c*k2**2+Z*beta2)/k1mk2
    pars%S03( 6)= 3.d0*Z
    pars%S03( 7)= 6.d0*(mu*vartheta_c+Z*beta1/k1**2)*c_im*k1/k1mk2
    pars%S03( 8)=-6.d0*(mu*vartheta_c+Z*beta2/k2**2)*c_im*k2/k1mk2
    pars%S03( 9)= 6.d0*Z/(c_im*k3)
    pars%S03(10)= 6.d0*(mu*vartheta_c+Z*beta1/k1**2)/k1mk2
    pars%S03(11)=-6.d0*(mu*vartheta_c+Z*beta2/k2**2)/k1mk2
    pars%S03(12)=-6.d0*Z/k3**2
    ! S1
    pars%S1( 1)= 3.d0*lambda/(lambda+2.d0*mu)
    pars%S1( 2)=-0.5d0*(k1**2*beta1-k2**2*beta2)/k1mk2
    pars%S1( 3)= k3**2
    pars%S1( 4)= 4.d0*c_im*k1*beta1/k1mk2
    pars%S1( 5)=-4.d0*c_im*k2*beta2/k1mk2
    pars%S1( 6)=-7.d0*c_im*k3
    pars%S1( 7)= 24.d0*beta1/k1mk2
    pars%S1( 8)=-24.d0*beta2/k1mk2
    pars%S1( 9)=-27.d0
    pars%S1(10)= 60.d0*beta1/(c_im*k1)/k1mk2
    pars%S1(11)=-60.d0*beta2/(c_im*k2)/k1mk2
    pars%S1(12)=-60.d0/(c_im*k3)
    pars%S1(13)=-60.d0*beta1/k1**2/k1mk2
    pars%S1(14)= 60.d0*beta2/k2**2/k1mk2
    pars%S1(15)= 60.d0/k3**2
    ! S2
    pars%S2( 1)= 6.d0*mu/(lambda+2.d0*mu)
    pars%S2( 2)= Q/R*vartheta_c/J-(lambda/mu+0.5d0)*(k1**2*beta1-k2**2*beta2)/k1mk2+0.5d0*k3**2
    pars%S2( 3)=-2.d0*(Q/R*vartheta_c/J-lambda/mu*beta1)*k1**2/k1mk2
    pars%S2( 4)= 2.d0*(Q/R*vartheta_c/J-lambda/mu*beta2)*k2**2/k1mk2
    pars%S2( 5)= 2.d0*(3.d0*Q/R*vartheta_c/J-(3.d0*lambda/mu-2.d0)*beta1)*c_im*k1/k1mk2
    pars%S2( 6)=-2.d0*(3.d0*Q/R*vartheta_c/J-(3.d0*lambda/mu-2.d0)*beta2)*c_im*k2/k1mk2
    pars%S2( 7)=-4.d0*c_im*k3
    pars%S2( 8)= 6.d0*(Q/R*vartheta_c/J-(lambda/mu-4.d0)*beta1)/k1mk2
    pars%S2( 9)=-6.d0*(Q/R*vartheta_c/J-(lambda/mu-4.d0)*beta2)/k1mk2
    pars%S2(10)=-24.d0
    pars%S2(11)= 60.d0*beta1/(c_im*k1)/k1mk2
    pars%S2(12)=-60.d0*beta2/(c_im*k2)/k1mk2
    pars%S2(13)=-60.d0/(c_im*k3)
    pars%S2(14)=-60.d0*beta1/k1**2/k1mk2
    pars%S2(15)= 60.d0*beta2/k2**2/k1mk2
    pars%S2(16)= 60.d0/k3**2
    ! S3
    pars%S3( 1)= 30.d0*(lambda+mu)/(lambda+2.d0*mu)
    pars%S3( 2)=-1.5d0*((k1**2*beta1-k2**2*beta2)/k1mk2-k3**2)
    pars%S3( 3)=-4.d0*k1**2*beta1/k1mk2
    pars%S3( 4)= 4.d0*k2**2*beta2/k1mk2
    pars%S3( 5)= 4.d0*k3**2
    pars%S3( 6)= 40.d0*c_im*k1*beta1/k1mk2
    pars%S3( 7)=-40.d0*c_im*k2*beta2/k1mk2
    pars%S3( 8)=-40.d0*c_im*k3
    pars%S3( 9)= 180.d0*beta1/k1mk2
    pars%S3(10)=-180.d0*beta2/k1mk2
    pars%S3(11)=-180.d0
    pars%S3(12)= 420.d0*beta1/(c_im*k1)/k1mk2
    pars%S3(13)=-420.d0*beta2/(c_im*k2)/k1mk2
    pars%S3(14)=-420.d0/(c_im*k3)
    pars%S3(15)=-420.d0*beta1/k1**2/k1mk2
    pars%S3(16)= 420.d0*beta2/k2**2/k1mk2
    pars%S3(17)= 420.d0/k3**2
    ! S4
    pars%S4( 1)= 2.0d0*mu/(lambda+2.d0*mu)
    pars%S4( 2)= 0.5d0*((k1**2*beta1-k2**2*beta2)/k1mk2+k3**2)
    pars%S4( 3)=-2.d0/5.d0*c_im*(2.d0/3.d0*(k1**3*beta1-k2**3*beta2)/k1mk2+k3**3)
    pars%S4( 4)= 2.d0*c_im*k3
    pars%S4( 5)=-4.d0*beta1/k1mk2
    pars%S4( 6)= 4.d0*beta2/k1mk2
    pars%S4( 7)= 6.d0
    pars%S4( 8)=-12.d0*beta1/(c_im*k1)/k1mk2
    pars%S4( 9)= 12.d0*beta2/(c_im*k2)/k1mk2
    pars%S4(10)= 12.d0/(c_im*k3)
    pars%S4(11)= 12.d0*beta1/k1**2/k1mk2
    pars%S4(12)=-12.d0*beta2/k2**2/k1mk2
    pars%S4(13)=-12.d0/k3**2
    ! S5
    pars%S5( 1)= 2.d0*(lambda-mu)/(lambda+2.d0*mu)
    pars%S5( 2)= Q**2/R**2/mu/J-2.d0*(lambda/mu+1.d0)*Q/R*vartheta_c/J&
               +(2.d0*lambda/mu+lambda**2/mu**2+0.5d0)*(k1**2*beta1-k2**2*beta2)/k1mk2-0.5d0*k3**2
    pars%S5( 3)=-c_im*(Q**2/R**2/mu/J*(k1*alpha1-k2*alpha2)/k1mk2-2.d0*(lambda/mu+2.d0/3.d0)*Q/R*vartheta_c/J*(k1**3-k2**3)/k1mk2&
                       +(4.d0/3.d0*lambda/mu+lambda**2/mu**2+4.d0/15.d0)*(k1**3*beta1-k2**3*beta2)/k1mk2-4.0d0/15.d0*k3**3)
    pars%S5( 4)= (Q**2/R**2/mu/J*alpha1-2.d0*Q/R*lambda/mu*vartheta_c/J*k1**2+lambda**2/mu**2*k1**2*beta1)/k1mk2
    pars%S5( 5)=-(Q**2/R**2/mu/J*alpha2-2.d0*Q/R*lambda/mu*vartheta_c/J*k2**2+lambda**2/mu**2*k2**2*beta2)/k1mk2
    pars%S5( 6)=-4.d0*(Q/R*vartheta_c/J-lambda/mu*beta1)*c_im*k1/k1mk2
    pars%S5( 7)= 4.d0*(Q/R*vartheta_c/J-lambda/mu*beta2)*c_im*k2/k1mk2
    pars%S5( 8)=-4.d0*(Q/R*vartheta_c/J-(lambda/mu-1.d0)*beta1)/k1mk2
    pars%S5( 9)= 4.d0*(Q/R*vartheta_c/J-(lambda/mu-1.d0)*beta2)/k1mk2
    pars%S5(10)= 4.d0
    pars%S5(11)=-12.d0*beta1/(c_im*k1)/k1mk2
    pars%S5(12)= 12.d0*beta2/(c_im*k2)/k1mk2
    pars%S5(13)= 12.d0/(c_im*k3)
    pars%S5(14)= 12.d0*beta1/k1**2/k1mk2
    pars%S5(15)=-12.d0*beta2/k2**2/k1mk2
    pars%S5(16)=-12.d0/k3**2
    ! Constants (common factors) of fundamental solutions
    ! u*
    pars%cte_u(0,0)=-c_1_4pi
    pars%cte_u(0,1:3)=-c_1_4pi
    pars%cte_u(1:3,0)=-c_1_4pi/pars%J
    pars%cte_u(1:3,1:3)=c_1_4pi/pars%mu
    ! t*
    pars%cte_t(0,0)=-c_1_4pi
    pars%cte_t(0,1:3)= c_1_4pi
    pars%cte_t(1:3,0)=-c_1_4pi/pars%mu
    pars%cte_t(1:3,1:3)=c_1_4pi
    ! d*
    pars%cte_d(0,0)=-c_1_4pi/pars%J
    pars%cte_d(0,1:3)=c_1_4pi/pars%mu
    pars%cte_d(1:3,0)=-c_1_4pi/pars%J
    pars%cte_d(1:3,1:3)=c_1_4pi
    ! s*
    pars%cte_s(0,0)=-c_1_4pi
    pars%cte_s(0,1:3)=c_1_4pi
    pars%cte_s(1:3,0)=-c_1_4pi
    pars%cte_s(1:3,1:3)=c_1_4pi*pars%mu
  end subroutine fbem_bem_harpor3d_calculate_parameters

!  ! Fundamental solutions tests
!  subroutine test_harpor3d_fs
!    implicit none
!    complex(kind=real64)               :: lambda, mu, RR, Q
!    real(kind=real64)                  :: rho1, rho2, rhoa, b
!    complex(kind=real64)               :: k1, k2, k3, c1, c2, c3
!    type(fbem_bem_harpor3d_parameters) :: p
!    real(kind=real64)                  :: omega, o_min, o_max, deltao
!    integer                            :: ko, no
!    real(kind=real64)                  :: r, r_min, r_max, deltar
!    integer                            :: kr, nr
!    real(kind=real64)                  :: d1r1, d1r2, d1r3, d1r4, d1r5
!    complex(kind=real64)               :: z1, z2, z3
!    complex(kind=real64)               :: E0z1, E1z1, E2z1, E3z1, E4z1, E5z1, E6z1
!    complex(kind=real64)               :: E0z2, E1z2, E2z2, E3z2, E4z2, E5z2, E6z2
!    complex(kind=real64)               :: E0z3, E1z3, E2z3, E3z3, E4z3, E5z3, E6z3
!    complex(kind=real64)               :: eta, vartheta, psi, chi
!    complex(kind=real64)               :: W0, T01, T02, W1, W2, T1, T2, T3
!    complex(kind=real64)               :: Q1, Q2, S01, S02, S03, S1, S2, S3, S4, S5
!    ! Parameters
!    lambda=3.5d7
!    mu=3.5d7
!    RR=7.7d8
!    Q=1.43d9
!    rho1=1723.d0
!    rho2=350.d0
!    rhoa=100.d0
!    b=1.201725d6
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
!      call fbem_bem_harpor3d_calculate_parameters(lambda,mu,rho1,rho2,rhoa,RR,Q,b,omega,p)
!      k1=p%k1
!      k2=p%k2
!      k3=p%k3
!      c1=p%c1
!      c2=p%c2
!      c3=p%c3
!      ! Loop through r
!      do kr=1,nr
!        r=10**(dlog10(r_min)+deltar*dble(kr-1))
!        d1r1=1.d0/r
!        d1r2=1.d0/r**2
!        d1r3=1.d0/r**3
!        d1r4=1.d0/r**4
!        d1r5=1.d0/r**5
!        z1=-c_im*k1*r
!        z2=-c_im*k2*r
!        z3=-c_im*k3*r
!        call fbem_decomposed_zexp_6(z1,E0z1,E1z1,E2z1,E3z1,E4z1,E5z1,E6z1)
!        call fbem_decomposed_zexp_6(z2,E0z2,E1z2,E2z2,E3z2,E4z2,E5z2,E6z2)
!        call fbem_decomposed_zexp_6(z3,E0z3,E1z3,E2z3,E3z3,E4z3,E5z3,E6z3)
!        !
!        ! Compact form of the components. Almost all of them break for k·r << 1
!        !
!        eta=p%eta(2)*d1r1*E0z1+p%eta(3)*d1r1*E0z2
!        vartheta=p%vartheta(2)*d1r1*E0z1+p%vartheta(3)*d1r1*E0z2+p%vartheta(4)*d1r2*E0z1+p%vartheta(5)*d1r2*E0z2
!        psi=d1r1*E0z3+p%psi(3)*d1r2*E0z1+p%psi(4)*d1r2*E0z2+p%psi(5)*d1r2*E0z3+p%psi(6)*d1r3*E0z1+p%psi(7)*d1r3*E0z2+p%psi(8)*d1r3*E0z3
!        chi=p%chi(2)*d1r1*E0z1+p%chi(3)*d1r1*E0z2+d1r1*E0z3+p%chi(4)*d1r2*E0z1+p%chi(5)*d1r2*E0z2+p%chi(6)*d1r2*E0z3+p%chi(7)*d1r3*E0z1+p%chi(8)*d1r3*E0z2+p%chi(9)*d1r3*E0z3
!        write(10,'(i11,e25.16e3,i11,12e25.16e3)') ko,omega,kr,r,zabs(k1),zabs(k2),zabs(k3),dreal(eta),dimag(eta),dreal(vartheta),dimag(vartheta),dreal(psi),dimag(psi),dreal(chi),dimag(chi)
!        !
!        ! Expanded form of the components. Mandatory for k·r << 1, some of them break for k·r >> 1 (for very large k·r).
!        !
!        ! Components of U*
!        eta=d1r1&
!           +p%eta(1)&
!           +p%eta(2)*d1r1*E2z1&
!           +p%eta(3)*d1r1*E2z2
!        vartheta=p%vartheta(1)&
!                +p%vartheta(2)*d1r1*E2z1&
!                +p%vartheta(3)*d1r1*E2z2&
!                +p%vartheta(4)*d1r2*E3z1&
!                +p%vartheta(5)*d1r2*E3z2
!        psi=p%psi(1)*d1r1&
!           +p%psi(2)+d1r1*E2z3&
!           +p%psi(3)*d1r2*E3z1&
!           +p%psi(4)*d1r2*E3z2&
!           +p%psi(5)*d1r2*E3z3&
!           +p%psi(6)*d1r3*E4z1&
!           +p%psi(7)*d1r3*E4z2&
!           +p%psi(8)*d1r3*E4z3
!        chi=p%chi(1)*d1r1&
!           +p%chi(2)*d1r1*E2z1&
!           +p%chi(3)*d1r1*E2z2&
!           +         d1r1*E2z3&
!           +p%chi(4)*d1r2*E3z1&
!           +p%chi(5)*d1r2*E3z2&
!           +p%chi(6)*d1r2*E3z3&
!           +p%chi(7)*d1r3*E4z1&
!           +p%chi(8)*d1r3*E4z2&
!           +p%chi(9)*d1r3*E4z3
!        ! Components of T* and D*
!        W0=p%W0(1)*d1r2&
!          +p%W0(2)&
!          +p%W0(3)*d1r1*E2z1&
!          +p%W0(4)*d1r1*E2z2&
!          +p%W0(5)*d1r2*E3z1&
!          +p%W0(6)*d1r2*E3z2
!        T01=p%T01(1)*d1r1&
!           +p%T01(2)*d1r1*E2z1&
!           +p%T01(3)*d1r1*E2z2&
!           +p%T01(4)*d1r2*E3z1&
!           +p%T01(5)*d1r2*E3z2&
!           +p%T01(6)*d1r3*E4z1&
!           +p%T01(7)*d1r3*E4z2
!        T02=p%T02(1)*d1r1&
!           +p%T02(2)&
!           +p%T02(3)*d1r1*E2z1&
!           +p%T02(4)*d1r1*E2z2&
!           +p%T02(5)*d1r2*E3z1&
!           +p%T02(6)*d1r2*E3z2&
!           +p%T02(7)*d1r3*E4z1&
!           +p%T02(8)*d1r3*E4z2
!        W1=p%W1( 1)*d1r1&
!          +p%W1( 2)*d1r1*E2z1&
!          +p%W1( 3)*d1r1*E2z2&
!          +p%W1( 4)*d1r1*E2z3&
!          +p%W1( 5)*d1r2*E3z1&
!          +p%W1( 6)*d1r2*E3z2&
!          +p%W1( 7)*d1r2*E3z3&
!          +p%W1( 8)*d1r3*E4z1&
!          +p%W1( 9)*d1r3*E4z2&
!          +p%W1(10)*d1r3*E4z3
!        W2=p%W2(1)*d1r1&
!          +p%W2(2)&
!          +p%W2(3)*d1r1*E2z3&
!          +p%W2(4)*d1r2*E3z1&
!          +p%W2(5)*d1r2*E3z2&
!          +p%W2(6)*d1r2*E3z3&
!          +p%W2(7)*d1r3*E4z1&
!          +p%W2(8)*d1r3*E4z2&
!          +p%W2(9)*d1r3*E4z3
!        T1=p%T1( 1)*d1r2&
!          +p%T1( 2)&
!          +p%T1( 3)*d1r1*E2z1&
!          +p%T1( 4)*d1r1*E2z2&
!          +p%T1( 5)*d1r1*E2z3&
!          +p%T1( 6)*d1r2*E3z1&
!          +p%T1( 7)*d1r2*E3z2&
!          +p%T1( 8)*d1r2*E3z3&
!          +p%T1( 9)*d1r3*E4z1&
!          +p%T1(10)*d1r3*E4z2&
!          +p%T1(11)*d1r3*E4z3&
!          +p%T1(12)*d1r4*E5z1&
!          +p%T1(13)*d1r4*E5z2&
!          +p%T1(14)*d1r4*E5z3
!        T2=p%T2( 1)*d1r2&
!          +p%T2( 2)&
!          +p%T2( 3)*d1r1*E2z3&
!          +p%T2( 4)*d1r2*E3z1&
!          +p%T2( 5)*d1r2*E3z2&
!          +p%T2( 6)*d1r2*E3z3&
!          +p%T2( 7)*d1r3*E4z1&
!          +p%T2( 8)*d1r3*E4z2&
!          +p%T2( 9)*d1r3*E4z3&
!          +p%T2(10)*d1r4*E5z1&
!          +p%T2(11)*d1r4*E5z2&
!          +p%T2(12)*d1r4*E5z3
!        T3=p%T3( 1)*d1r2&
!          +p%T3( 2)&
!          +p%T3( 3)*d1r1*E2z1&
!          +p%T3( 4)*d1r1*E2z2&
!          +p%T3( 5)*d1r2*E3z1&
!          +p%T3( 6)*d1r2*E3z2&
!          +p%T3( 7)*d1r2*E3z3&
!          +p%T3( 8)*d1r3*E4z1&
!          +p%T3( 9)*d1r3*E4z2&
!          +p%T3(10)*d1r3*E4z3&
!          +p%T3(11)*d1r4*E5z1&
!          +p%T3(12)*d1r4*E5z2&
!          +p%T3(13)*d1r4*E5z3
!        ! Components of S*
!        Q1=p%Q1( 1)*d1r3&
!          +p%Q1( 2)*d1r1&
!          +p%Q1( 3)*d1r1*E2z1&
!          +p%Q1( 4)*d1r1*E2z2&
!          +p%Q1( 5)*d1r1*E2z3&
!          +p%Q1( 6)*d1r2*E3z1&
!          +p%Q1( 7)*d1r2*E3z2&
!          +p%Q1( 8)*d1r2*E3z3&
!          +p%Q1( 9)*d1r3*E4z1&
!          +p%Q1(10)*d1r3*E4z2&
!          +p%Q1(11)*d1r3*E4z3
!        Q2=p%Q2( 1)*d1r3&
!          +p%Q2( 2)*d1r1&
!          +p%Q2( 3)&
!          +p%Q2( 4)*d1r1*E2z3&
!          +p%Q2( 5)*d1r2*E3z1&
!          +p%Q2( 6)*d1r2*E3z2&
!          +p%Q2( 7)*d1r2*E3z3&
!          +p%Q2( 8)*d1r3*E4z1&
!          +p%Q2( 9)*d1r3*E4z2&
!          +p%Q2(10)*d1r3*E4z3
!        S01=p%S01( 1)*d1r2&
!           +p%S01( 2)&
!           +p%S01( 3)*d1r1*E2z1&
!           +p%S01( 4)*d1r1*E2z2&
!           +p%S01( 5)*d1r1*E2z3&
!           +p%S01( 6)*d1r2*E3z1&
!           +p%S01( 7)*d1r2*E3z2&
!           +p%S01( 8)*d1r2*E3z3&
!           +p%S01( 9)*d1r3*E4z1&
!           +p%S01(10)*d1r3*E4z2&
!           +p%S01(11)*d1r3*E4z3&
!           +p%S01(12)*d1r4*E5z1&
!           +p%S01(13)*d1r4*E5z2&
!           +p%S01(14)*d1r4*E5z3
!        S02=p%S02( 1)*d1r2&
!           +p%S02( 2)&
!           +p%S02( 3)*d1r1*E2z1&
!           +p%S02( 4)*d1r1*E2z2&
!           +p%S02( 5)*d1r2*E3z1&
!           +p%S02( 6)*d1r2*E3z2&
!           +p%S02( 7)*d1r2*E3z3&
!           +p%S02( 8)*d1r3*E4z1&
!           +p%S02( 9)*d1r3*E4z2&
!           +p%S02(10)*d1r3*E4z3&
!           +p%S02(11)*d1r4*E5z1&
!           +p%S02(12)*d1r4*E5z2&
!           +p%S02(13)*d1r4*E5z3
!        S03=p%S03( 1)*d1r2&
!           +p%S03( 2)&
!           +p%S03( 3)*d1r1*E2z3&
!           +p%S03( 4)*d1r2*E3z1&
!           +p%S03( 5)*d1r2*E3z2&
!           +p%S03( 6)*d1r2*E3z3&
!           +p%S03( 7)*d1r3*E4z1&
!           +p%S03( 8)*d1r3*E4z2&
!           +p%S03( 9)*d1r3*E4z3&
!           +p%S03(10)*d1r4*E5z1&
!           +p%S03(11)*d1r4*E5z2&
!           +p%S03(12)*d1r4*E5z3
!        S1=p%S1( 1)*d1r3&
!          +p%S1( 2)*d1r1&
!          +p%S1( 3)*d1r1*E2z3&
!          +p%S1( 4)*d1r2*E3z1&
!          +p%S1( 5)*d1r2*E3z2&
!          +p%S1( 6)*d1r2*E3z3&
!          +p%S1( 7)*d1r3*E4z1&
!          +p%S1( 8)*d1r3*E4z2&
!          +p%S1( 9)*d1r3*E4z3&
!          +p%S1(10)*d1r4*E5z1&
!          +p%S1(11)*d1r4*E5z2&
!          +p%S1(12)*d1r4*E5z3&
!          +p%S1(13)*d1r5*E6z1&
!          +p%S1(14)*d1r5*E6z2&
!          +p%S1(15)*d1r5*E6z3
!        S2=p%S2( 1)*d1r3&
!          +p%S2( 2)*d1r1&
!          +p%S2( 3)*d1r1*E2z1&
!          +p%S2( 4)*d1r1*E2z2&
!          +p%S2( 5)*d1r2*E3z1&
!          +p%S2( 6)*d1r2*E3z2&
!          +p%S2( 7)*d1r2*E3z3&
!          +p%S2( 8)*d1r3*E4z1&
!          +p%S2( 9)*d1r3*E4z2&
!          +p%S2(10)*d1r3*E4z3&
!          +p%S2(11)*d1r4*E5z1&
!          +p%S2(12)*d1r4*E5z2&
!          +p%S2(13)*d1r4*E5z3&
!          +p%S2(14)*d1r5*E6z1&
!          +p%S2(15)*d1r5*E6z2&
!          +p%S2(16)*d1r5*E6z3
!        S3=p%S3( 1)*d1r3&
!          +p%S3( 2)*d1r1&
!          +p%S3( 3)*d1r1*E2z1&
!          +p%S3( 4)*d1r1*E2z2&
!          +p%S3( 5)*d1r1*E2z3&
!          +p%S3( 6)*d1r2*E3z1&
!          +p%S3( 7)*d1r2*E3z2&
!          +p%S3( 8)*d1r2*E3z3&
!          +p%S3( 9)*d1r3*E4z1&
!          +p%S3(10)*d1r3*E4z2&
!          +p%S3(11)*d1r3*E4z3&
!          +p%S3(12)*d1r4*E5z1&
!          +p%S3(13)*d1r4*E5z2&
!          +p%S3(14)*d1r4*E5z3&
!          +p%S3(15)*d1r5*E6z1&
!          +p%S3(16)*d1r5*E6z2&
!          +p%S3(17)*d1r5*E6z3
!        S4=p%S4( 1)*d1r3&
!          +p%S4( 2)*d1r1&
!          +p%S4( 3)&
!          +p%S4( 4)*d1r2*E3z3&
!          +p%S4( 5)*d1r3*E4z1&
!          +p%S4( 6)*d1r3*E4z2&
!          +p%S4( 7)*d1r3*E4z3&
!          +p%S4( 8)*d1r4*E5z1&
!          +p%S4( 9)*d1r4*E5z2&
!          +p%S4(10)*d1r4*E5z3&
!          +p%S4(11)*d1r5*E6z1&
!          +p%S4(12)*d1r5*E6z2&
!          +p%S4(13)*d1r5*E6z3
!        S5=p%S5( 1)*d1r3&
!          +p%S5( 2)*d1r1&
!          +p%S5( 3)&
!          +p%S5( 4)*d1r1*E2z1&
!          +p%S5( 5)*d1r1*E2z2&
!          +p%S5( 6)*d1r2*E3z1&
!          +p%S5( 7)*d1r2*E3z2&
!          +p%S5( 8)*d1r3*E4z1&
!          +p%S5( 9)*d1r3*E4z2&
!          +p%S5(10)*d1r3*E4z3&
!          +p%S5(11)*d1r4*E5z1&
!          +p%S5(12)*d1r4*E5z2&
!          +p%S5(13)*d1r4*E5z3&
!          +p%S5(14)*d1r5*E6z1&
!          +p%S5(15)*d1r5*E6z2&
!          +p%S5(16)*d1r5*E6z3
!        write(11,'(i11,e25.16e3,i11,12e25.16e3)') ko,omega,kr,r,zabs(k1),zabs(k2),zabs(k3),dreal(eta),dimag(eta),&
!                                                                                           dreal(vartheta),dimag(vartheta),&
!                                                                                           dreal(psi),dimag(psi),&
!                                                                                           dreal(chi),dimag(chi)
!      end do
!      write(10,*)
!      write(10,*)
!      write(11,*)
!      write(11,*)
!    end do
!  end subroutine test_harpor3d_fs

  ! ============================================================================================================================== !
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE)                                                                                     !
  ! ============================================================================================================================== !

  subroutine fbem_bem_harpor3d_sbie_ext_pre(ps,e,reverse,x_i,p,h,g)
    implicit none
    ! I/O
    integer                            :: ps                    !! Selected precalculated dataset (quasi-singular part)
    type(fbem_bem_element)             :: e                     !! Element
    logical                            :: reverse               !! Reverse normal vector
    real(kind=real64)                  :: x_i(3)                !! Collocation point position vector
    type(fbem_bem_harpor3d_parameters) :: p                     !! Parameters of the region
    complex(kind=real64)               :: h(e%n_pnodes,0:3,0:3) !! h integration kernels matrix
    complex(kind=real64)               :: g(e%n_snodes,0:3,0:3) !! g integration kernels matrix
    ! Local
    integer              :: kip, kphi                           ! Counter variable for integration points loop
    real(kind=real64)    :: x(3)                                ! Position vector at integration point
    real(kind=real64)    :: n(3)                                ! Unit normal vector at integration point
    real(kind=real64)    :: pphijw(e%n_pnodes)                  ! phi^p * jacobian * weight at integration point
    real(kind=real64)    :: sphijw(e%n_snodes)                  ! phi^s * jacobian * weight at integration point
    real(kind=real64)    :: rv(3)                               ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, d1r1, d1r2, d1r3, d1r4           ! Distance vector module and its inverse
    real(kind=real64)    :: drdx(3)
    real(kind=real64)    :: drdn                                ! Partial derivative of r respect to unit normal
    integer              :: il, ik                              ! Counter for load / observation components
    complex(kind=real64) :: z(3)                                ! Arguments z=-ikr
    complex(kind=real64) :: EnR(0:6,3)                          ! Exponential function decomposition for each wavenumber
    complex(kind=real64) :: eta, vartheta, psi, chi             ! Components of the fundamental solution U*
    complex(kind=real64) :: W0, T01, T02, W1, W2, TT1, TT2, TT3 ! Components of the fundamental solution T*
    complex(kind=real64) :: fs_u(0:3,0:3), fs_t(0:3,0:3)        ! Fundamental solutions matrices
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
      d1r1=1.d0/r
      d1r2=d1r1**2
      d1r3=d1r2*d1r1
      d1r4=d1r3*d1r1
      drdx=rv*d1r1
      drdn=dot_product(drdx,n)
      z(1)=-c_im*p%k1*r
      z(2)=-c_im*p%k2*r
      z(3)=-c_im*p%k3*r
      call fbem_zexp_decomposed(3,z,EnR)
      EnR(2,:)=EnR(2,:)*d1r1
      EnR(3,:)=EnR(3,:)*d1r2
      EnR(4,:)=EnR(4,:)*d1r3
      EnR(5,:)=EnR(5,:)*d1r4
      eta=d1r1+p%eta(1)+p%eta(2)*EnR(2,1)+p%eta(3)*EnR(2,2)
      vartheta=p%vartheta(1)+p%vartheta(2)*EnR(2,1)+p%vartheta(3)*EnR(2,2)+p%vartheta(4)*EnR(3,1)+p%vartheta(5)*EnR(3,2)
      psi=p%psi(1)*d1r1+p%psi(2)+EnR(2,3)+p%psi(3)*EnR(3,1)+p%psi(4)*EnR(3,2)+p%psi(5)*EnR(3,3)+p%psi(6)*EnR(4,1)&
         +p%psi(7)*EnR(4,2)+p%psi(8)*EnR(4,3)
      chi=p%chi(1)*d1r1+p%chi(2)*EnR(2,1)+p%chi(3)*EnR(2,2)+EnR(2,3)+p%chi(4)*EnR(3,1)+p%chi(5)*EnR(3,2)+p%chi(6)*EnR(3,3)&
         +p%chi(7)*EnR(4,1)+p%chi(8)*EnR(4,2)+p%chi(9)*EnR(4,3)
      W0 =p%W0(1)*d1r2+p%W0(2)+p%W0(3)*EnR(2,1)+p%W0(4)*EnR(2,2)+p%W0(5)*EnR(3,1)+p%W0(6)*EnR(3,2)
      T01=p%T01(1)*d1r1+p%T01(2)*EnR(2,1)+p%T01(3)*EnR(2,2)+p%T01(4)*EnR(3,1)+p%T01(5)*EnR(3,2)+p%T01(6)*EnR(4,1)&
         +p%T01(7)*EnR(4,2)
      T02=p%T02(1)*d1r1+p%T02(2)+p%T02(3)*EnR(2,1)+p%T02(4)*EnR(2,2)+p%T02(5)*EnR(3,1)+p%T02(6)*EnR(3,2)+p%T02(7)*EnR(4,1)&
         +p%T02(8)*EnR(4,2)
      W1 =p%W1(1)*d1r1+p%W1(2)*EnR(2,1)+p%W1(3)*EnR(2,2)+p%W1(4)*EnR(2,3)+p%W1(5)*EnR(3,1)+p%W1(6)*EnR(3,2)+p%W1(7)*EnR(3,3)&
         +p%W1(8)*EnR(4,1)+p%W1(9)*EnR(4,2)+p%W1(10)*EnR(4,3)
      W2 =p%W2(1)*d1r1+p%W2(2)+p%W2(3)*EnR(2,3)+p%W2(4)*EnR(3,1)+p%W2(5)*EnR(3,2)+p%W2(6)*EnR(3,3)+p%W2(7)*EnR(4,1)&
         +p%W2(8)*EnR(4,2)+p%W2(9)*EnR(4,3)
      TT1=p%T1(1)*d1r2+p%T1(2)+p%T1(3)*EnR(2,1)+p%T1(4)*EnR(2,2)+p%T1(5)*EnR(2,3)+p%T1(6)*EnR(3,1)+p%T1(7)*EnR(3,2)&
         +p%T1(8)*EnR(3,3)+p%T1(9)*EnR(4,1)+p%T1(10)*EnR(4,2)+p%T1(11)*EnR(4,3)+p%T1(12)*EnR(5,1)+p%T1(13)*EnR(5,2)&
         +p%T1(14)*EnR(5,3)
      TT2=p%T2(1)*d1r2+p%T2(2)+p%T2(3)*EnR(2,3)+p%T2(4)*EnR(3,1)+p%T2(5)*EnR(3,2)+p%T2(6)*EnR(3,3)+p%T2(7)*EnR(4,1)&
         +p%T2(8)*EnR(4,2)+p%T2(9)*EnR(4,3)+p%T2(10)*EnR(5,1)+p%T2(11)*EnR(5,2)+p%T2(12)*EnR(5,3)
      TT3=p%T3(1)*d1r2+p%T3(2)+p%T3(3)*EnR(2,1)+p%T3(4)*EnR(2,2)+p%T3(5)*EnR(3,1)+p%T3(6)*EnR(3,2)+p%T3(7)*EnR(3,3)&
         +p%T3(8)*EnR(4,1)+p%T3(9)*EnR(4,2)+p%T3(10)*EnR(4,3)+p%T3(11)*EnR(5,1)+p%T3(12)*EnR(5,2)+p%T3(13)*EnR(5,3)
      ! Fundamental solutions
      fs_u(0,0)=eta
      fs_u(0,1:3)=vartheta*drdx
      fs_u(1:3,0)=vartheta*drdx
      fs_t(0,0)=W0*drdn
      fs_t(0,1:3)=T01*drdx*drdn+T02*n
      fs_t(1:3,0)=W1*drdx*drdn+W2*n
      do ik=1,3
        do il=1,3
          fs_u(il,ik)=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
          fs_t(il,ik)=TT1*drdx(il)*drdx(ik)*drdn+TT2*(drdn*c_dkr(il,ik)+drdx(ik)*n(il))+TT3*drdx(il)*n(ik)
        end do
      end do
      ! Add the term of the numerical integration summation (fundamental solutions * functional shape functions * jacobian * weight)
      do ik=0,3
        do il=0,3
          h(:,il,ik)=h(:,il,ik)+fs_t(il,ik)*pphijw
          g(:,il,ik)=g(:,il,ik)+fs_u(il,ik)*sphijw
        end do
      end do
    end do
    ! Multiply by constants
    do ik=0,3
      do il=0,3
        h(:,il,ik)=p%cte_t(il,ik)*h(:,il,ik)
        g(:,il,ik)=p%cte_u(il,ik)*g(:,il,ik)
      end do
    end do
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_harpor3d_sbie_ext_pre

  subroutine fbem_bem_harpor3d_sbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,p,gln,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                                !! Integration element
    logical                            :: reverse                          !! Reverse normal vector
    real(kind=real64)                  :: xi_s(2,fbem_n_vertices(e%gtype)) !! Coordinates of the subdivision of the parent element
    real(kind=real64)                  :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)                  :: barxip(2)                        !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)                  :: barr                             !! Telles jacobian at barxip
    type(fbem_bem_harpor3d_parameters) :: p                                !! Parameters of the region
    integer                            :: gln                              !! 1D Gauss-Legendre number of integration points (<=32)
    complex(kind=real64)               :: h(e%n_pnodes,0:3,0:3)            !! h kernel vector
    complex(kind=real64)               :: g(e%n_snodes,0:3,0:3)            !! g kernel vector
    ! Local
    integer                      :: kphi                                ! Counter variable for shape functions loops
    integer                      :: k1                                  ! Counter variable for reference coordinate xi_1
    integer                      :: k2                                  ! Counter variable for reference coordinate xi_2
    real(kind=real64)            :: aux(10)                             ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)                    ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi1(e%n_gnodes)               ! Geometrical shape functions derivatives values
    real(kind=real64)            :: dgphidxi2(e%n_gnodes)               ! Geometrical shape functions derivatives values
    real(kind=real64)            :: pphi(e%n_pnodes)                    ! Functional shape functions values
    real(kind=real64)            :: sphi(e%n_snodes)                    ! Functional shape functions values
    real(kind=real64)            :: gamma(2)                            ! Vector of gamma_1,gamma_2 coordinates
    real(kind=real64)            :: w(2)                                ! Weights of the integration rule
    real(kind=real64)            :: xip(2)                              ! Vector of xip_1,xip_2 coordinates
    real(kind=real64)            :: dxidxi1p(2), dxidxi2p(2)            ! xi derivatives with respect to xip
    real(kind=real64)            :: js                                  ! Subdivision jacobian
    real(kind=real64)            :: xi(2)                               ! Vector of xi_1,xi_2 coordinates
    real(kind=real64)            :: xipp(2)                             ! Coordinate xipp used for quadrilateral-triangle transformation
    real(kind=real64)            :: barxipp(2)                          ! Coordinate xipp of collocation point
    real(kind=real64)            :: jqt                                 ! Jacobian of the quadrilateral-triangle transformation
    type(fbem_telles_parameters) :: telles_parameters(2)                ! Telles parameters for each coordinate
    real(kind=real64)            :: jt(2)                               ! Telles jacobian for each coordinate: xi_1->gamma_1 and xi_2->gamma_2
    real(kind=real64)            :: x(3)                                ! Position vector at xi_1,xi_2
    real(kind=real64)            :: T1(3), T2(3)                        ! Tangent vectors at xi_1,xi_2
    real(kind=real64)            :: N(3)                                ! Normal vector at xi_1,xi_2
    real(kind=real64)            :: rv(3)                               ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)            :: r, d1r1, d1r2, d1r3, d1r4           ! Distance vector module and its inverse
    real(kind=real64)            :: drdx(3)                             ! Distance vector derivatives with respect to x_k
    real(kind=real64)            :: jg                                  ! Geometric jacobian
    real(kind=real64)            :: jw                                  ! Jacobians * weights
    real(kind=real64)            :: drdn                                ! Partial derivative of r respect to unit normal
    real(kind=real64)            :: pphijw(e%n_pnodes)                  ! Functional shape functions * jw
    real(kind=real64)            :: sphijw(e%n_snodes)                  ! Functional shape functions * jw
    integer                      :: il, ik                              ! Counter for load / observation components
    complex(kind=real64)         :: z(3)                                ! Arguments z=-ikr
    complex(kind=real64)         :: EnR(0:6,3)                          ! Exponential function decomposition for each wavenumber
    complex(kind=real64)         :: eta, vartheta, psi, chi             ! Components of the fundamental solution U*
    complex(kind=real64)         :: W0, T01, T02, W1, W2, TT1, TT2, TT3 ! Components of the fundamental solution T*
    complex(kind=real64)         :: fs_u(0:3,0:3), fs_t(0:3,0:3)        ! Fundamental solutions matrices
    complex(kind=real64)         :: cte_u(0:3,0:3), cte_t(0:3,0:3)      ! Constants of fundamental solutions matrices
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
            z(3)=-c_im*p%k3*r
            call fbem_zexp_decomposed(3,z,EnR)
            EnR(2,:)=EnR(2,:)*d1r1
            EnR(3,:)=EnR(3,:)*d1r2
            EnR(4,:)=EnR(4,:)*d1r3
            EnR(5,:)=EnR(5,:)*d1r4
            eta=d1r1+p%eta(1)+p%eta(2)*EnR(2,1)+p%eta(3)*EnR(2,2)
            vartheta=p%vartheta(1)+p%vartheta(2)*EnR(2,1)+p%vartheta(3)*EnR(2,2)+p%vartheta(4)*EnR(3,1)+p%vartheta(5)*EnR(3,2)
            psi=p%psi(1)*d1r1+p%psi(2)+EnR(2,3)+p%psi(3)*EnR(3,1)+p%psi(4)*EnR(3,2)+p%psi(5)*EnR(3,3)+p%psi(6)*EnR(4,1)&
               +p%psi(7)*EnR(4,2)+p%psi(8)*EnR(4,3)
            chi=p%chi(1)*d1r1+p%chi(2)*EnR(2,1)+p%chi(3)*EnR(2,2)+EnR(2,3)+p%chi(4)*EnR(3,1)+p%chi(5)*EnR(3,2)+p%chi(6)*EnR(3,3)&
               +p%chi(7)*EnR(4,1)+p%chi(8)*EnR(4,2)+p%chi(9)*EnR(4,3)
            W0 =p%W0(1)*d1r2+p%W0(2)+p%W0(3)*EnR(2,1)+p%W0(4)*EnR(2,2)+p%W0(5)*EnR(3,1)+p%W0(6)*EnR(3,2)
            T01=p%T01(1)*d1r1+p%T01(2)*EnR(2,1)+p%T01(3)*EnR(2,2)+p%T01(4)*EnR(3,1)+p%T01(5)*EnR(3,2)+p%T01(6)*EnR(4,1)&
               +p%T01(7)*EnR(4,2)
            T02=p%T02(1)*d1r1+p%T02(2)+p%T02(3)*EnR(2,1)+p%T02(4)*EnR(2,2)+p%T02(5)*EnR(3,1)+p%T02(6)*EnR(3,2)+p%T02(7)*EnR(4,1)&
               +p%T02(8)*EnR(4,2)
            W1 =p%W1(1)*d1r1+p%W1(2)*EnR(2,1)+p%W1(3)*EnR(2,2)+p%W1(4)*EnR(2,3)+p%W1(5)*EnR(3,1)+p%W1(6)*EnR(3,2)+p%W1(7)*EnR(3,3)&
               +p%W1(8)*EnR(4,1)+p%W1(9)*EnR(4,2)+p%W1(10)*EnR(4,3)
            W2 =p%W2(1)*d1r1+p%W2(2)+p%W2(3)*EnR(2,3)+p%W2(4)*EnR(3,1)+p%W2(5)*EnR(3,2)+p%W2(6)*EnR(3,3)+p%W2(7)*EnR(4,1)&
               +p%W2(8)*EnR(4,2)+p%W2(9)*EnR(4,3)
            TT1=p%T1(1)*d1r2+p%T1(2)+p%T1(3)*EnR(2,1)+p%T1(4)*EnR(2,2)+p%T1(5)*EnR(2,3)+p%T1(6)*EnR(3,1)+p%T1(7)*EnR(3,2)&
               +p%T1(8)*EnR(3,3)+p%T1(9)*EnR(4,1)+p%T1(10)*EnR(4,2)+p%T1(11)*EnR(4,3)+p%T1(12)*EnR(5,1)+p%T1(13)*EnR(5,2)&
               +p%T1(14)*EnR(5,3)
            TT2=p%T2(1)*d1r2+p%T2(2)+p%T2(3)*EnR(2,3)+p%T2(4)*EnR(3,1)+p%T2(5)*EnR(3,2)+p%T2(6)*EnR(3,3)+p%T2(7)*EnR(4,1)&
               +p%T2(8)*EnR(4,2)+p%T2(9)*EnR(4,3)+p%T2(10)*EnR(5,1)+p%T2(11)*EnR(5,2)+p%T2(12)*EnR(5,3)
            TT3=p%T3(1)*d1r2+p%T3(2)+p%T3(3)*EnR(2,1)+p%T3(4)*EnR(2,2)+p%T3(5)*EnR(3,1)+p%T3(6)*EnR(3,2)+p%T3(7)*EnR(3,3)&
               +p%T3(8)*EnR(4,1)+p%T3(9)*EnR(4,2)+p%T3(10)*EnR(4,3)+p%T3(11)*EnR(5,1)+p%T3(12)*EnR(5,2)+p%T3(13)*EnR(5,3)
            ! Fundamental solutions
            fs_u(0,0)=eta
            fs_t(0,0)=W0*drdn
            fs_u(0,1:3)=vartheta*drdx
            fs_t(0,1:3)=T01*drdx*drdn+T02*n
            fs_u(1:3,0)=vartheta*drdx
            fs_t(1:3,0)=W1*drdx*drdn+W2*n
            do il=1,3
              do ik=1,3
                fs_u(il,ik)=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
                fs_t(il,ik)=TT1*drdx(il)*drdx(ik)*drdn+TT2*(drdn*c_dkr(il,ik)+drdx(ik)*n(il))+TT3*drdx(il)*n(ik)
              end do
            end do
            ! Add the term of the numerical integration summation (fundamental solutions * functional shape functions * jacobian * weight)
            do il=0,3
              do ik=0,3
                h(:,il,ik)=h(:,il,ik)+fs_t(il,ik)*pphijw
                g(:,il,ik)=g(:,il,ik)+fs_u(il,ik)*sphijw
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
            z(3)=-c_im*p%k3*r
            call fbem_zexp_decomposed(3,z,EnR)
            EnR(2,:)=EnR(2,:)*d1r1
            EnR(3,:)=EnR(3,:)*d1r2
            EnR(4,:)=EnR(4,:)*d1r3
            EnR(5,:)=EnR(5,:)*d1r4
            eta=d1r1+p%eta(1)+p%eta(2)*EnR(2,1)+p%eta(3)*EnR(2,2)
            vartheta=p%vartheta(1)+p%vartheta(2)*EnR(2,1)+p%vartheta(3)*EnR(2,2)+p%vartheta(4)*EnR(3,1)+p%vartheta(5)*EnR(3,2)
            psi=p%psi(1)*d1r1+p%psi(2)+EnR(2,3)+p%psi(3)*EnR(3,1)+p%psi(4)*EnR(3,2)+p%psi(5)*EnR(3,3)+p%psi(6)*EnR(4,1)&
               +p%psi(7)*EnR(4,2)+p%psi(8)*EnR(4,3)
            chi=p%chi(1)*d1r1+p%chi(2)*EnR(2,1)+p%chi(3)*EnR(2,2)+EnR(2,3)+p%chi(4)*EnR(3,1)+p%chi(5)*EnR(3,2)+p%chi(6)*EnR(3,3)&
               +p%chi(7)*EnR(4,1)+p%chi(8)*EnR(4,2)+p%chi(9)*EnR(4,3)
            W0 =p%W0(1)*d1r2+p%W0(2)+p%W0(3)*EnR(2,1)+p%W0(4)*EnR(2,2)+p%W0(5)*EnR(3,1)+p%W0(6)*EnR(3,2)
            T01=p%T01(1)*d1r1+p%T01(2)*EnR(2,1)+p%T01(3)*EnR(2,2)+p%T01(4)*EnR(3,1)+p%T01(5)*EnR(3,2)+p%T01(6)*EnR(4,1)&
               +p%T01(7)*EnR(4,2)
            T02=p%T02(1)*d1r1+p%T02(2)+p%T02(3)*EnR(2,1)+p%T02(4)*EnR(2,2)+p%T02(5)*EnR(3,1)+p%T02(6)*EnR(3,2)+p%T02(7)*EnR(4,1)&
               +p%T02(8)*EnR(4,2)
            W1 =p%W1(1)*d1r1+p%W1(2)*EnR(2,1)+p%W1(3)*EnR(2,2)+p%W1(4)*EnR(2,3)+p%W1(5)*EnR(3,1)+p%W1(6)*EnR(3,2)+p%W1(7)*EnR(3,3)&
               +p%W1(8)*EnR(4,1)+p%W1(9)*EnR(4,2)+p%W1(10)*EnR(4,3)
            W2 =p%W2(1)*d1r1+p%W2(2)+p%W2(3)*EnR(2,3)+p%W2(4)*EnR(3,1)+p%W2(5)*EnR(3,2)+p%W2(6)*EnR(3,3)+p%W2(7)*EnR(4,1)&
               +p%W2(8)*EnR(4,2)+p%W2(9)*EnR(4,3)
            TT1=p%T1(1)*d1r2+p%T1(2)+p%T1(3)*EnR(2,1)+p%T1(4)*EnR(2,2)+p%T1(5)*EnR(2,3)+p%T1(6)*EnR(3,1)+p%T1(7)*EnR(3,2)&
               +p%T1(8)*EnR(3,3)+p%T1(9)*EnR(4,1)+p%T1(10)*EnR(4,2)+p%T1(11)*EnR(4,3)+p%T1(12)*EnR(5,1)+p%T1(13)*EnR(5,2)&
               +p%T1(14)*EnR(5,3)
            TT2=p%T2(1)*d1r2+p%T2(2)+p%T2(3)*EnR(2,3)+p%T2(4)*EnR(3,1)+p%T2(5)*EnR(3,2)+p%T2(6)*EnR(3,3)+p%T2(7)*EnR(4,1)&
               +p%T2(8)*EnR(4,2)+p%T2(9)*EnR(4,3)+p%T2(10)*EnR(5,1)+p%T2(11)*EnR(5,2)+p%T2(12)*EnR(5,3)
            TT3=p%T3(1)*d1r2+p%T3(2)+p%T3(3)*EnR(2,1)+p%T3(4)*EnR(2,2)+p%T3(5)*EnR(3,1)+p%T3(6)*EnR(3,2)+p%T3(7)*EnR(3,3)&
               +p%T3(8)*EnR(4,1)+p%T3(9)*EnR(4,2)+p%T3(10)*EnR(4,3)+p%T3(11)*EnR(5,1)+p%T3(12)*EnR(5,2)+p%T3(13)*EnR(5,3)
            ! Fundamental solutions
            fs_u(0,0)=eta
            fs_t(0,0)=W0*drdn
            fs_u(0,1:3)=vartheta*drdx
            fs_t(0,1:3)=T01*drdx*drdn+T02*n
            fs_u(1:3,0)=vartheta*drdx
            fs_t(1:3,0)=W1*drdx*drdn+W2*n
            do il=1,3
              do ik=1,3
                fs_u(il,ik)=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
                fs_t(il,ik)=TT1*drdx(il)*drdx(ik)*drdn+TT2*(drdn*c_dkr(il,ik)+drdx(ik)*n(il))+TT3*drdx(il)*n(ik)
              end do
            end do
            ! Add the term of the numerical integration summation (fundamental solutions * functional shape functions * jacobian * weight)
            do il=0,3
              do ik=0,3
                h(:,il,ik)=h(:,il,ik)+fs_t(il,ik)*pphijw
                g(:,il,ik)=g(:,il,ik)+fs_u(il,ik)*sphijw
              end do
            end do
          end do
        end do
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'n_edges not valid')
    end select
    ! Multiply by constants
    do ik=0,3
      do il=0,3
        h(:,il,ik)=p%cte_t(il,ik)*h(:,il,ik)
        g(:,il,ik)=p%cte_u(il,ik)*g(:,il,ik)
      end do
    end do
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_harpor3d_sbie_ext_st

  recursive subroutine fbem_bem_harpor3d_sbie_ext_adp(e,reverse,xi_s,x_i,p,qsp,ks,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                                !! Element
    logical                            :: reverse                          !! Reverse orientation
    real(kind=real64)                  :: xi_s(2,fbem_n_vertices(e%gtype)) !! Subdivision of the parent element
    real(kind=real64)                  :: x_i(3)                           !! Collocation point position vector
    type(fbem_bem_harpor3d_parameters) :: p                                !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp                              !! Quasi-singular integration parameters
    integer                            :: ks                               !! Current level of subdivisions
    integer                            :: ns                               !! Maximum level of subdivision
    complex(kind=real64)               :: h(e%n_pnodes,0:3,0:3)            !! h integration kernels matrix
    complex(kind=real64)               :: g(e%n_snodes,0:3,0:3)            !! g integration kernels matrix
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
    complex(kind=real64) :: h_tmp(e%n_pnodes,0:3,0:3)            ! h integration kernels matrix (temporary)
    complex(kind=real64) :: g_tmp(e%n_snodes,0:3,0:3)            ! g integration kernels matrix (temporary)
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
        call fbem_warning_message(error_unit,0,'fbem_bem_harpor3d_sbie_ext_adp',ns,'maximum number of subdivisions reached')
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
          call fbem_bem_harpor3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
          ! SUBTRI 2
          tmp_xi_s(:,1)=xi_s(:,2)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          call fbem_bem_harpor3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
          ! SUBTRI 3
          tmp_xi_s(:,1)=xi_s(:,3)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          call fbem_bem_harpor3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
          ! SUBTRI 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          call fbem_bem_harpor3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
        ! QUADRILATERALS
        case (4)
          ! SUBQUAD 1
          tmp_xi_s(:,1)=xi_s(:,1)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,3)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          call fbem_bem_harpor3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
          ! SUBQUAD 2
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=xi_s(:,2)
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,4)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          call fbem_bem_harpor3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
          ! SUBQUAD 3
          tmp_xi_s(:,1)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=xi_s(:,3)
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          call fbem_bem_harpor3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
          ! SUBQUAD 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          tmp_xi_s(:,2)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=xi_s(:,4)
          call fbem_bem_harpor3d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
      end select
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_harpor3d_sbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,p,gln,h_tmp,g_tmp)
      h=h+h_tmp
      g=g+g_tmp
    end if
  end subroutine fbem_bem_harpor3d_sbie_ext_adp

  subroutine fbem_bem_harpor3d_sbie_int(type_g,type_f1,type_f2,delta_f,x_nodes,reverse,xi_i,p,h,g)
    implicit none
    ! I/O
    integer                            :: type_g                           !! Geometrial interpolation
    integer                            :: type_f1                          !! Functional interpolation (primary variables)
    integer                            :: type_f2                          !! Functional interpolation (secondary variables)
    real(kind=real64)                  :: delta_f                          !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64)                  :: x_nodes(3,fbem_n_nodes(type_g))  !! Position vectors of geometrical nodes
    logical                            :: reverse                          !! Normal vector inversion
    real(kind=real64)                  :: xi_i(2)                          !! Reference coordinates of the singular point.
    type(fbem_bem_harpor3d_parameters) :: p                                !! Parameters of the region
    complex(kind=real64)               :: h(fbem_n_nodes(type_f1),0:3,0:3) !! h kernel vector
    complex(kind=real64)               :: g(fbem_n_nodes(type_f2),0:3,0:3) !! g kernel vector
    ! Local
    integer              :: ksubtri                                   ! Counter variable for subtriangles loop
    integer              :: nsubtri                                   ! Number of subtriangles
    integer              :: subtriangle(8)                            ! Vector that contains what subtriangles need to be integrated
    real(kind=real64)    :: theta_subtri(2,8)                         ! Matrix that contains the angles theta of the subtriangles to be integrated
    real(kind=real64)    :: thetap_subtri(2,8)                        ! Matrix that contains the angles thetap of the subtriangles to be integrated
    integer              :: ktheta                                    ! Counter variable for theta coordinate loop
    integer              :: krho                                      ! Counter variable for rho coordinate loop
    integer              :: kphi                                      ! Counter coordinates for shape functions loops
    integer              :: nnodes_g, n_pnodes, n_snodes              ! Number of nodes of geometrical interpolation
    integer              :: ngp_theta                                 ! Number of Gauss points for theta coordinate
    integer              :: ngp_rho                                   ! Number of Gauss points for rho coordinate
    real(kind=real64)    :: thetai, thetaf, thetapi, thetapf          ! Initial and final angles for subtriangle integration
    real(kind=real64)    :: w_angular                                 ! Weight of the angular coordinate
    real(kind=real64)    :: w_radial                                  ! Weight of the radial coordinate
    real(kind=real64)    :: theta                                     ! Angle coordinate theta
    real(kind=real64)    :: thetap                                    ! Angle coordinate thetap
    real(kind=real64)    :: thetapp                                   ! Angle coordinate thetap on [0,1] domain
    real(kind=real64)    :: jthetap                                   ! thetap->thetapp jacobian
    real(kind=real64)    :: rhoij                                     ! Maximum rho (radial) value for each edge
    real(kind=real64)    :: rho                                       ! Radial coordinate rho
    real(kind=real64)    :: rhop                                      ! Radial coordinate rho on [0,1] domain
    real(kind=real64)    :: aux(10)                                   ! Auxiliary variable for shape functions resources
    real(kind=real64)    :: xi(2)                                     ! Reference xi_1,xi_2 coordinates
    real(kind=real64)    :: phi_f1(fbem_n_nodes(type_f1))             ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)    :: phi_f2(fbem_n_nodes(type_f2))             ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)    :: phi_f1_i(fbem_n_nodes(type_f1))           ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)    :: phi_g(fbem_n_nodes(type_g))               ! Geometrical shape functions values at xi_1,xi_2
    real(kind=real64)    :: dphidxi1_g(fbem_n_nodes(type_g))          ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64)    :: dphidxi2_g(fbem_n_nodes(type_g))          ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64)    :: x_i(3)                                    ! Collocation point position vector
    real(kind=real64)    :: x(3)                                      ! Position vector at xi_1,xi_2
    real(kind=real64)    :: T1(3), T2(3)                              ! Tangent vectors at xi_1,xi_2
    real(kind=real64)    :: N(3)                                      ! Normal vector at xi_1,xi_2
    real(kind=real64)    :: rv(3)                                     ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, d1r1, d1r2, d1r3, d1r4                 ! Distance vector module and its inverse
    real(kind=real64)    :: drdx(3)                                   ! Distance vector derivatives with respect to x_k
    real(kind=real64)    :: jg                                        ! Geometric jacobian
    real(kind=real64)    :: jw                                        ! Jacobian * weights
    real(kind=real64)    :: drdn                                      ! Partial derivative of r respect to unit normal
    real(kind=real64)    :: phif1jw(fbem_n_nodes(type_f1))            ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)    :: phif2jw(fbem_n_nodes(type_f2))            ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)    :: costheta, sintheta                        ! cos(theta), sin(theta)
    integer              :: il, ik                                    ! Counter for load / observation components
    complex(kind=real64) :: z(3)                                ! Arguments z=-ikr
    complex(kind=real64) :: EnR(0:6,3)                          ! Exponential function decomposition for each wavenumber
    complex(kind=real64) :: eta, vartheta, psi, chi             ! Components of the fundamental solution U*
    complex(kind=real64) :: W0, T01, T02, W1, W2, TT1, TT2, TT3 ! Components of the fundamental solution T*
    complex(kind=real64) :: fs_u(0:3,0:3), fs_t(0:3,0:3)              ! Fundamental solutions matrices
    complex(kind=real64) :: cte_u(0:3,0:3), cte_t(0:3,0:3)            ! Constants of fundamental solutions matrices
    ! Local variables associated with line integrals
    real(kind=real64)              :: re                              ! Quasi-singular integration relative error
    integer                        :: gln_min                         ! Minimum Number of 1D Gauss-Legendre number of points
    integer                        :: ns_max                          ! Maximum number of subdivisions
    type(fbem_qs_parameters)       :: qsp                             ! Quasi-singular integration parameters
    real(kind=real64)              :: xi_s(1,2)                       ! Edge subdivision
    integer                        :: kedge                           ! Counter variable for edges
    integer                        :: nedges                          ! Number of edges
    logical                        :: integrate                       ! Control variable
    integer                        :: type_edge                       ! Line element type for line integrals
    integer                        :: nnodes_edge                     ! Number of nodes of the edge
    real(kind=real64), allocatable :: x_nodes_edge(:,:)               ! Coordinates of the edge elements
    real(kind=real64)              :: hli(3,3)                        ! Line integrals values
    !
    ! Initialization
    !
    ! Kernel vectors
    h=(0.d0,0.d0)
    g=(0.d0,0.d0)
    ! Number of nodes of geometrical interpolation
    nnodes_g=fbem_n_nodes(type_g)
    n_pnodes=fbem_n_nodes(type_f1)
    n_snodes=fbem_n_nodes(type_f2)
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
          jg=sqrt(dot_product(N,N))
          ! Unit normal vector
          n=N/jg
          ! Distance vector
          rv=x-x_i
          ! Distance vector norm
          r=sqrt(dot_product(rv,rv))
          d1r1=1.d0/r
          d1r2=d1r1**2
          d1r3=d1r2*d1r1
          d1r4=d1r3*d1r1
          ! r_{,k}
          drdx=rv*d1r1
          ! dr/dn
          drdn=dot_product(drdx,n)
          ! Jacobians * weights
          jw=jg*rho*jthetap*w_angular*w_radial
          ! FUNCTIONAL SHAPE FUNCTIONS
          ! Functional shape functions (primary variables) at xi
!#         define etype type_f1
!#         define delta delta_f
!#         define phi phi_f1
!#         include <phi_2d.rc>
!#         undef etype
!#         undef delta
!#         undef phi
!          ! Functional shape functions (secondary variables) at xi
!#         define etype type_f2
!#         define delta delta_f
!#         define phi phi_f2
!#         include <phi_2d.rc>
!#         undef etype
!#         undef delta
!#         undef phi
          phi_f1=phi_g
          phi_f2=phi_g
          ! Functional shape functions * jacobians * weights
          phif1jw=phi_f1*jw
          phif2jw=phi_f2*jw
          !
          ! ONLY AT MOST WEAKLY SINGULAR PARTS
          !
          ! COMPONENTS OF THE FUNDAMENTAL SOLUTION
          z(1)=-c_im*p%k1*r
          z(2)=-c_im*p%k2*r
          z(3)=-c_im*p%k3*r
          call fbem_zexp_decomposed(3,z,EnR)
          EnR(2,:)=EnR(2,:)*d1r1
          EnR(3,:)=EnR(3,:)*d1r2
          EnR(4,:)=EnR(4,:)*d1r3
          EnR(5,:)=EnR(5,:)*d1r4
          eta=d1r1+p%eta(1)+p%eta(2)*EnR(2,1)+p%eta(3)*EnR(2,2)
          vartheta=p%vartheta(1)+p%vartheta(2)*EnR(2,1)+p%vartheta(3)*EnR(2,2)+p%vartheta(4)*EnR(3,1)+p%vartheta(5)*EnR(3,2)
          psi=p%psi(1)*d1r1+p%psi(2)+EnR(2,3)+p%psi(3)*EnR(3,1)+p%psi(4)*EnR(3,2)+p%psi(5)*EnR(3,3)+p%psi(6)*EnR(4,1)&
             +p%psi(7)*EnR(4,2)+p%psi(8)*EnR(4,3)
          chi=p%chi(1)*d1r1+p%chi(2)*EnR(2,1)+p%chi(3)*EnR(2,2)+EnR(2,3)+p%chi(4)*EnR(3,1)+p%chi(5)*EnR(3,2)+p%chi(6)*EnR(3,3)&
             +p%chi(7)*EnR(4,1)+p%chi(8)*EnR(4,2)+p%chi(9)*EnR(4,3)
          W0 =p%W0(2)+p%W0(3)*EnR(2,1)+p%W0(4)*EnR(2,2)+p%W0(5)*EnR(3,1)+p%W0(6)*EnR(3,2)
          T01=p%T01(1)*d1r1+p%T01(2)*EnR(2,1)+p%T01(3)*EnR(2,2)+p%T01(4)*EnR(3,1)+p%T01(5)*EnR(3,2)+p%T01(6)*EnR(4,1)&
             +p%T01(7)*EnR(4,2)
          T02=p%T02(1)*d1r1+p%T02(2)+p%T02(3)*EnR(2,1)+p%T02(4)*EnR(2,2)+p%T02(5)*EnR(3,1)+p%T02(6)*EnR(3,2)+p%T02(7)*EnR(4,1)&
             +p%T02(8)*EnR(4,2)
          W1 =p%W1(1)*d1r1+p%W1(2)*EnR(2,1)+p%W1(3)*EnR(2,2)+p%W1(4)*EnR(2,3)+p%W1(5)*EnR(3,1)+p%W1(6)*EnR(3,2)+p%W1(7)*EnR(3,3)&
             +p%W1(8)*EnR(4,1)+p%W1(9)*EnR(4,2)+p%W1(10)*EnR(4,3)
          W2 =p%W2(1)*d1r1+p%W2(2)+p%W2(3)*EnR(2,3)+p%W2(4)*EnR(3,1)+p%W2(5)*EnR(3,2)+p%W2(6)*EnR(3,3)+p%W2(7)*EnR(4,1)&
             +p%W2(8)*EnR(4,2)+p%W2(9)*EnR(4,3)
          TT1=p%T1(2)+p%T1(3)*EnR(2,1)+p%T1(4)*EnR(2,2)+p%T1(5)*EnR(2,3)+p%T1(6)*EnR(3,1)+p%T1(7)*EnR(3,2)&
             +p%T1(8)*EnR(3,3)+p%T1(9)*EnR(4,1)+p%T1(10)*EnR(4,2)+p%T1(11)*EnR(4,3)+p%T1(12)*EnR(5,1)+p%T1(13)*EnR(5,2)&
             +p%T1(14)*EnR(5,3)
          TT2=p%T2(2)+p%T2(3)*EnR(2,3)+p%T2(4)*EnR(3,1)+p%T2(5)*EnR(3,2)+p%T2(6)*EnR(3,3)+p%T2(7)*EnR(4,1)&
             +p%T2(8)*EnR(4,2)+p%T2(9)*EnR(4,3)+p%T2(10)*EnR(5,1)+p%T2(11)*EnR(5,2)+p%T2(12)*EnR(5,3)
          TT3=p%T3(2)+p%T3(3)*EnR(2,1)+p%T3(4)*EnR(2,2)+p%T3(5)*EnR(3,1)+p%T3(6)*EnR(3,2)+p%T3(7)*EnR(3,3)&
             +p%T3(8)*EnR(4,1)+p%T3(9)*EnR(4,2)+p%T3(10)*EnR(4,3)+p%T3(11)*EnR(5,1)+p%T3(12)*EnR(5,2)+p%T3(13)*EnR(5,3)
          ! FUNDAMENTAL SOLUTIONS
          ! Fluid load / Fluid response
          fs_u(0,0)=eta
          fs_t(0,0)=(p%W0(1)*d1r2+W0)*drdn
          ! Fluid load / Solid response
          fs_u(0,1:3)=vartheta*drdx
          fs_t(0,1:3)=T01*drdx*drdn+T02*n
          ! Solid load / Fluid response
          fs_u(1:3,0)=vartheta*drdx
          fs_t(1:3,0)=W1*drdx*drdn+W2*n
          ! Solid load / Solid response
          do ik=1,3
            do il=1,3
              fs_u(il,ik)=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
              fs_t(il,ik)=(p%T1(1)*d1r2+TT1)*drdx(il)*drdx(ik)*drdn+(p%T2(1)*d1r2+TT2)*drdn*c_dkr(il,ik)+TT2*drdx(ik)*n(il)+TT3*drdx(il)*n(ik)
            end do
          end do
          ! ADD INTEGRAND EVALUATION
          do ik=0,3
            do il=0,3
              h(:,il,ik)=h(:,il,ik)+fs_t(il,ik)*phif1jw
              g(:,il,ik)=g(:,il,ik)+fs_u(il,ik)*phif2jw
            end do
          end do

!      do kphi=1,n_pnodes
!        h(kphi,:,:)=h(kphi,:,:)+fs_t*phif1jw(kphi)
!      end do
!      do kphi=1,n_snodes
!        g(kphi,:,:)=g(kphi,:,:)+fs_u*phif2jw(kphi)
!      end do

          ! ADD THE REGULAR PART OF THE CPV INTEGRAL PRESENT IN THE SOLID SKELETON PART
          do ik=1,3
            do il=1,3
              h(:,il,ik)=h(:,il,ik)+p%T2(1)*d1r2*(n(il)*drdx(ik)-n(ik)*drdx(il))*(phi_f1-phi_f1_i)*jw
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
    do ik=1,3
      do il=1,3
        h(:,il,ik)=h(:,il,ik)+phi_f1_i*p%T2(1)*hli(il,ik)
      end do
    end do
    ! Multiply by constants
    do ik=0,3
      do il=0,3
        h(:,il,ik)=p%cte_t(il,ik)*h(:,il,ik)
        g(:,il,ik)=p%cte_u(il,ik)*g(:,il,ik)
      end do
    end do
    ! If the normal has to be reversed, then h=-h
    if (reverse) h=-h
  end subroutine fbem_bem_harpor3d_sbie_int









  subroutine fbem_bem_harpor3d_sbie_auto(e,reverse,x_i,p,qsp,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                     !! Integration element
    logical                            :: reverse               !! Reverse orientation
    real(kind=real64)                  :: x_i(3)                !! Collocation point
    type(fbem_bem_harpor3d_parameters) :: p                     !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp                   !! Quasi-singular integration parameters
    integer                            :: ns                    !! Maximum level of subdivisions
    complex(kind=real64)               :: h(e%n_pnodes,0:3,0:3) !! h integration kernel
    complex(kind=real64)               :: g(e%n_snodes,0:3,0:3) !! g integration kernel
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
        call fbem_bem_harpor3d_sbie_int(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi,p,h,g)
        !call fbem_bem_harpor3d_sbie_int_s(e%gtype,e%n_gnodes,e%x,reverse,barxi,p,h,g)
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
          call fbem_bem_harpor3d_sbie_ext_pre(ps,e,reverse,x_i,p,h,g)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_harpor3d_sbie_ext_adp(e,reverse,xi_s,x_i,p,qsp,1,ns,h,g)
        end if
    end select
  end subroutine fbem_bem_harpor3d_sbie_auto

  ! ================================================================================================================================

  ! ==============================================================================================================================
  ! HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (HBIE)

  subroutine fbem_bem_harpor3d_hbie_ext_pre(ps,e,reverse,x_i,n_i,p,m,l)
    implicit none
    ! I/O
    integer                            :: ps                    !! Selected precalculated dataset
    type(fbem_bem_element)             :: e                     !! Element
    logical                            :: reverse               !! Reverse normal vector
    real(kind=real64)                  :: x_i(3)                !! Collocation point position vector
    real(kind=real64)                  :: n_i(3)                !! Collocation point unit normal vector
    type(fbem_bem_harpor3d_parameters) :: p                     !! Parameters of the region
    complex(kind=real64)               :: m(e%n_pnodes,0:3,0:3) !! m integration kernels vector
    complex(kind=real64)               :: l(e%n_snodes,0:3,0:3) !! l integration kernels vector              !! l integration kernels vector
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
    complex(kind=real64) :: z(3)                            ! Arguments z=ikr
    complex(kind=real64) :: EnR(0:6,3)                      ! Exponential function decomposition for each wavenumber
    complex(kind=real64) :: W0, T01, T02, W1, W2, TT1, TT2, TT3       ! Components of the fundamental solution D*
    complex(kind=real64) :: Q1, Q2, S01, S02, S03, S1, S2, S3, S4, S5 ! Components of the fundamental solution S*
    complex(kind=real64) :: fs_d(0:3,0:3), fs_s(0:3,0:3)              ! Fundamental solutions matrices
    complex(kind=real64) :: cte_d(0:3,0:3), cte_s(0:3,0:3)            ! Constants of fundamental solutions matrices
    ! Initialize
    m=(0.d0,0.d0)
    l=(0.d0,0.d0)
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      ! Functional shape functions * jacobian * weight
      pphijw=e%ps_pphijw(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      ! Calculation of the components of the fundamental solution
      x=e%ps_x(:,kip,ps)
      n=e%ps_n(:,kip,ps)
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
      z(3)=-c_im*p%k3*r
      call fbem_zexp_decomposed(3,z,EnR)
      EnR(2,:)=EnR(2,:)*d1r1
      EnR(3,:)=EnR(3,:)*d1r2
      EnR(4,:)=EnR(4,:)*d1r3
      EnR(5,:)=EnR(5,:)*d1r4
      EnR(6,:)=EnR(6,:)*d1r5
      W0=p%W0(1)*d1r2+p%W0(2)+p%W0(3)*EnR(2,1)+p%W0(4)*EnR(2,2)+p%W0(5)*EnR(3,1)+p%W0(6)*EnR(3,2)
      T01=p%T01(1)*d1r1+p%T01(2)*EnR(2,1)+p%T01(3)*EnR(2,2)+p%T01(4)*EnR(3,1)+p%T01(5)*EnR(3,2)+p%T01(6)*EnR(4,1)&
         +p%T01(7)*EnR(4,2)
      T02=p%T02(1)*d1r1+p%T02(2)+p%T02(3)*EnR(2,1)+p%T02(4)*EnR(2,2)+p%T02(5)*EnR(3,1)+p%T02(6)*EnR(3,2)&
         +p%T02(7)*EnR(4,1)+p%T02(8)*EnR(4,2)
      W1=p%W1(1)*d1r1+p%W1(2)*EnR(2,1)+p%W1(3)*EnR(2,2)+p%W1(4)*EnR(2,3)+p%W1(5)*EnR(3,1)+p%W1(6)*EnR(3,2)&
        +p%W1(7)*EnR(3,3)+p%W1(8)*EnR(4,1)+p%W1(9)*EnR(4,2)+p%W1(10)*EnR(4,3)
      W2=p%W2(1)*d1r1+p%W2(2)+p%W2(3)*EnR(2,3)+p%W2(4)*EnR(3,1)+p%W2(5)*EnR(3,2)+p%W2(6)*EnR(3,3)+p%W2(7)*EnR(4,1)&
        +p%W2(8)*EnR(4,2)+p%W2(9)*EnR(4,3)
      TT1=p%T1(1)*d1r2+p%T1(2)+p%T1(3)*EnR(2,1)+p%T1(4)*EnR(2,2)+p%T1(5)*EnR(2,3)+p%T1(6)*EnR(3,1)+p%T1(7)*EnR(3,2)&
         +p%T1(8)*EnR(3,3)+p%T1(9)*EnR(4,1)+p%T1(10)*EnR(4,2)+p%T1(11)*EnR(4,3)+p%T1(12)*EnR(5,1)+p%T1(13)*EnR(5,2)&
         +p%T1(14)*EnR(5,3)
      TT2=p%T2(1)*d1r2+p%T2(2)+p%T2(3)*EnR(2,3)+p%T2(4)*EnR(3,1)+p%T2(5)*EnR(3,2)+p%T2(6)*EnR(3,3)+p%T2(7)*EnR(4,1)&
         +p%T2(8)*EnR(4,2)+p%T2(9)*EnR(4,3)+p%T2(10)*EnR(5,1)+p%T2(11)*EnR(5,2)+p%T2(12)*EnR(5,3)
      TT3=p%T3(1)*d1r2+p%T3(2)+p%T3(3)*EnR(2,1)+p%T3(4)*EnR(2,2)+p%T3(5)*EnR(3,1)+p%T3(6)*EnR(3,2)+p%T3(7)*EnR(3,3)&
         +p%T3(8)*EnR(4,1)+p%T3(9)*EnR(4,2)+p%T3(10)*EnR(4,3)+p%T3(11)*EnR(5,1)+p%T3(12)*EnR(5,2)+p%T3(13)*EnR(5,3)
      Q1=p%Q1(1)*d1r3+p%Q1(2)*d1r1+p%Q1(3)*EnR(2,1)+p%Q1(4)*EnR(2,2)+p%Q1(5)*EnR(2,3)+p%Q1(6)*EnR(3,1)&
        +p%Q1(7)*EnR(3,2)+p%Q1(8)*EnR(3,3)+p%Q1(9)*EnR(4,1)+p%Q1(10)*EnR(4,2)+p%Q1(11)*EnR(4,3)
      Q2=p%Q2(1)*d1r3+p%Q2(2)*d1r1+p%Q2(3)+p%Q2(4)*EnR(2,3)+p%Q2(5)*EnR(3,1)+p%Q2(6)*EnR(3,2)+p%Q2(7)*EnR(3,3)&
        +p%Q2(8)*EnR(4,1)+p%Q2(9)*EnR(4,2)+p%Q2(10)*EnR(4,3)
      S01=p%S01(1)*d1r2+p%S01(2)+p%S01(3)*EnR(2,1)+p%S01(4)*EnR(2,2)+p%S01(5)*EnR(2,3)+p%S01(6)*EnR(3,1)&
         +p%S01(7)*EnR(3,2)+p%S01(8)*EnR(3,3)+p%S01(9)*EnR(4,1)+p%S01(10)*EnR(4,2)+p%S01(11)*EnR(4,3)&
         +p%S01(12)*EnR(5,1)+p%S01(13)*EnR(5,2)+p%S01(14)*EnR(5,3)
      S02=p%S02(1)*d1r2+p%S02(2)+p%S02(3)*EnR(2,1)+p%S02(4)*EnR(2,2)+p%S02(5)*EnR(3,1)+p%S02(6)*EnR(3,2)&
         +p%S02(7)*EnR(3,3)+p%S02(8)*EnR(4,1)+p%S02(9)*EnR(4,2)+p%S02(10)*EnR(4,3)+p%S02(11)*EnR(5,1)&
         +p%S02(12)*EnR(5,2)+p%S02(13)*EnR(5,3)
      S03=p%S03(1)*d1r2+p%S03(2)+p%S03(3)*EnR(2,3)+p%S03(4)*EnR(3,1)+p%S03(5)*EnR(3,2)+p%S03(6)*EnR(3,3)&
         +p%S03(7)*EnR(4,1)+p%S03(8)*EnR(4,2)+p%S03(9)*EnR(4,3)+p%S03(10)*EnR(5,1)+p%S03(11)*EnR(5,2)&
         +p%S03(12)*EnR(5,3)
      S1=p%S1(1)*d1r3+p%S1(2)*d1r1+p%S1(3)*EnR(2,3)+p%S1(4)*EnR(3,1)+p%S1(5)*EnR(3,2)+p%S1(6)*EnR(3,3)&
        +p%S1(7)*EnR(4,1)+p%S1(8)*EnR(4,2)+p%S1(9)*EnR(4,3)+p%S1(10)*EnR(5,1)+p%S1(11)*EnR(5,2)+p%S1(12)*EnR(5,3)&
        +p%S1(13)*EnR(6,1)+p%S1(14)*EnR(6,2)+p%S1(15)*EnR(6,3)
      S2=p%S2(1)*d1r3+p%S2(2)*d1r1+p%S2(3)*EnR(2,1)+p%S2(4)*EnR(2,2)+p%S2(5)*EnR(3,1)+p%S2(6)*EnR(3,2)&
        +p%S2(7)*EnR(3,3)+p%S2(8)*EnR(4,1)+p%S2(9)*EnR(4,2)+p%S2(10)*EnR(4,3)+p%S2(11)*EnR(5,1)+p%S2(12)*EnR(5,2)&
        +p%S2(13)*EnR(5,3)+p%S2(14)*EnR(6,1)+p%S2(15)*EnR(6,2)+p%S2(16)*EnR(6,3)
      S3=p%S3(1)*d1r3+p%S3(2)*d1r1+p%S3(3)*EnR(2,1)+p%S3(4)*EnR(2,2)+p%S3(5)*EnR(2,3)+p%S3(6)*EnR(3,1)&
        +p%S3(7)*EnR(3,2)+p%S3(8)*EnR(3,3)+p%S3(9)*EnR(4,1)+p%S3(10)*EnR(4,2)+p%S3(11)*EnR(4,3)+p%S3(12)*EnR(5,1)&
        +p%S3(13)*EnR(5,2)+p%S3(14)*EnR(5,3)+p%S3(15)*EnR(6,1)+p%S3(16)*EnR(6,2)+p%S3(17)*EnR(6,3)
      S4=p%S4(1)*d1r3+p%S4(2)*d1r1+p%S4(3)+p%S4(4)*EnR(3,3)+p%S4(5)*EnR(4,1)+p%S4(6)*EnR(4,2)+p%S4(7)*EnR(4,3)&
        +p%S4(8)*EnR(5,1)+p%S4(9)*EnR(5,2)+p%S4(10)*EnR(5,3)+p%S4(11)*EnR(6,1)+p%S4(12)*EnR(6,2)+p%S4(13)*EnR(6,3)
      S5=p%S5(1)*d1r3+p%S5(2)*d1r1+p%S5(3)+p%S5(4)*EnR(2,1)+p%S5(5)*EnR(2,2)+p%S5(6)*EnR(3,1)+p%S5(7)*EnR(3,2)&
        +p%S5(8)*EnR(4,1)+p%S5(9)*EnR(4,2)+p%S5(10)*EnR(4,3)+p%S5(11)*EnR(5,1)+p%S5(12)*EnR(5,2)+p%S5(13)*EnR(5,3)&
        +p%S5(14)*EnR(6,1)+p%S5(15)*EnR(6,2)+p%S5(16)*EnR(6,3)
      ! Fundamental solutions
      fs_d(0,0)=W0*drdni
      fs_s(0,0)=Q1*drdn*drdni+Q2*n_dot_ni
      fs_d(0,1:3)=-W1*drdx*drdni+W2*n_i
      fs_s(0,1:3)=S01*drdx*drdn*drdni+S02*n*drdni+S03*(n_i*drdn+drdx*n_dot_ni)
      fs_d(1:3,0)=-T01*drdx*drdni+T02*n_i
      fs_s(1:3,0)=-S01*drdx*drdn*drdni+S02*n_i*drdn-S03*(-n*drdni+drdx*n_dot_ni)
      do il=1,3
        do ik=1,3
          fs_d(il,ik)=TT1*drdx(il)*drdx(ik)*drdni-TT2*(drdx(il)*n_i(ik)-c_dkr(il,ik)*drdni)-TT3*drdx(ik)*n_i(il)
          fs_s(il,ik)=S1*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni-c_dkr(il,ik)*drdn*drdni+drdx(il)*drdx(ik)*n_dot_ni)&
                     +S2*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni)+S3*drdx(il)*drdx(ik)*drdn*drdni&
                     +S4*(c_dkr(il,ik)*n_dot_ni+n_i(ik)*n(il))+S5*n(ik)*n_i(il)
        end do
      end do
      ! Add the term of the numerical integration summation (fundamental solutions * functional shape functions * jacobian * weight)
      do il=0,3
        do ik=0,3
          m(:,il,ik)=m(:,il,ik)+fs_s(il,ik)*pphijw
          l(:,il,ik)=l(:,il,ik)+fs_d(il,ik)*sphijw
        end do
      end do
    end do
    ! Multiply by constants
    cte_d(0,0)=-c_1_4pi/p%J
    cte_s(0,0)=-c_1_4pi
    cte_d(0,1:3)=c_1_4pi/p%mu
    cte_s(0,1:3)=c_1_4pi
    cte_d(1:3,0)=-c_1_4pi/p%J
    cte_s(1:3,0)=-c_1_4pi
    cte_d(1:3,1:3)=c_1_4pi
    cte_s(1:3,1:3)=c_1_4pi*p%mu
    do il=0,3
      do ik=0,3
        m(:,il,ik)=cte_s(il,ik)*m(:,il,ik)
        l(:,il,ik)=cte_d(il,ik)*l(:,il,ik)
      end do
    end do
    ! Reverse if needed
    if (reverse) m=-m
  end subroutine fbem_bem_harpor3d_hbie_ext_pre

  subroutine fbem_bem_harpor3d_hbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,p,gln,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                                !! Integration element
    logical                            :: reverse                          !! Reverse normal vector
    real(kind=real64)                  :: xi_s(2,fbem_n_vertices(e%gtype)) !! Coordinates of the subdivision of the parent element
    real(kind=real64)                  :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)                  :: n_i(3)                           !! Unit normal at the collocation point
    real(kind=real64)                  :: barxip(2)                        !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)                  :: barr                             !! Telles jacobian at barxip
    type(fbem_bem_harpor3d_parameters) :: p                                !! Parameters of the region
    integer                            :: gln                              !! 1D Gauss-Legendre number of integration points (<=32)
    complex(kind=real64)               :: m(e%n_pnodes,0:3,0:3)            !! m kernel vector
    complex(kind=real64)               :: l(e%n_snodes,0:3,0:3)            !! l kernel vector
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
    complex(kind=real64)         :: z(3)                            ! Arguments z=-ikr
    complex(kind=real64)         :: EnR(0:6,3)                      ! Exponential function decomposition for each wavenumber
    complex(kind=real64)         :: W0, T01, T02, W1, W2, TT1, TT2, TT3       ! Components of the fundamental solution D*
    complex(kind=real64)         :: Q1, Q2, S01, S02, S03, S1, S2, S3, S4, S5 ! Components of the fundamental solution S*
    complex(kind=real64)         :: fs_d(0:3,0:3), fs_s(0:3,0:3)              ! Fundamental solutions matrices
    complex(kind=real64)         :: cte_d(0:3,0:3), cte_s(0:3,0:3)            ! Constants of fundamental solutions matrices
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
            z(3)=-c_im*p%k3*r
            call fbem_zexp_decomposed(3,z,EnR)
            EnR(2,:)=EnR(2,:)*d1r1
            EnR(3,:)=EnR(3,:)*d1r2
            EnR(4,:)=EnR(4,:)*d1r3
            EnR(5,:)=EnR(5,:)*d1r4
            EnR(6,:)=EnR(6,:)*d1r5
            W0=p%W0(1)*d1r2+p%W0(2)+p%W0(3)*EnR(2,1)+p%W0(4)*EnR(2,2)+p%W0(5)*EnR(3,1)+p%W0(6)*EnR(3,2)
            T01=p%T01(1)*d1r1+p%T01(2)*EnR(2,1)+p%T01(3)*EnR(2,2)+p%T01(4)*EnR(3,1)+p%T01(5)*EnR(3,2)+p%T01(6)*EnR(4,1)&
               +p%T01(7)*EnR(4,2)
            T02=p%T02(1)*d1r1+p%T02(2)+p%T02(3)*EnR(2,1)+p%T02(4)*EnR(2,2)+p%T02(5)*EnR(3,1)+p%T02(6)*EnR(3,2)&
               +p%T02(7)*EnR(4,1)+p%T02(8)*EnR(4,2)
            W1=p%W1(1)*d1r1+p%W1(2)*EnR(2,1)+p%W1(3)*EnR(2,2)+p%W1(4)*EnR(2,3)+p%W1(5)*EnR(3,1)+p%W1(6)*EnR(3,2)&
              +p%W1(7)*EnR(3,3)+p%W1(8)*EnR(4,1)+p%W1(9)*EnR(4,2)+p%W1(10)*EnR(4,3)
            W2=p%W2(1)*d1r1+p%W2(2)+p%W2(3)*EnR(2,3)+p%W2(4)*EnR(3,1)+p%W2(5)*EnR(3,2)+p%W2(6)*EnR(3,3)+p%W2(7)*EnR(4,1)&
              +p%W2(8)*EnR(4,2)+p%W2(9)*EnR(4,3)
            TT1=p%T1(1)*d1r2+p%T1(2)+p%T1(3)*EnR(2,1)+p%T1(4)*EnR(2,2)+p%T1(5)*EnR(2,3)+p%T1(6)*EnR(3,1)+p%T1(7)*EnR(3,2)&
               +p%T1(8)*EnR(3,3)+p%T1(9)*EnR(4,1)+p%T1(10)*EnR(4,2)+p%T1(11)*EnR(4,3)+p%T1(12)*EnR(5,1)+p%T1(13)*EnR(5,2)&
               +p%T1(14)*EnR(5,3)
            TT2=p%T2(1)*d1r2+p%T2(2)+p%T2(3)*EnR(2,3)+p%T2(4)*EnR(3,1)+p%T2(5)*EnR(3,2)+p%T2(6)*EnR(3,3)+p%T2(7)*EnR(4,1)&
               +p%T2(8)*EnR(4,2)+p%T2(9)*EnR(4,3)+p%T2(10)*EnR(5,1)+p%T2(11)*EnR(5,2)+p%T2(12)*EnR(5,3)
            TT3=p%T3(1)*d1r2+p%T3(2)+p%T3(3)*EnR(2,1)+p%T3(4)*EnR(2,2)+p%T3(5)*EnR(3,1)+p%T3(6)*EnR(3,2)+p%T3(7)*EnR(3,3)&
               +p%T3(8)*EnR(4,1)+p%T3(9)*EnR(4,2)+p%T3(10)*EnR(4,3)+p%T3(11)*EnR(5,1)+p%T3(12)*EnR(5,2)+p%T3(13)*EnR(5,3)
            Q1=p%Q1(1)*d1r3+p%Q1(2)*d1r1+p%Q1(3)*EnR(2,1)+p%Q1(4)*EnR(2,2)+p%Q1(5)*EnR(2,3)+p%Q1(6)*EnR(3,1)&
              +p%Q1(7)*EnR(3,2)+p%Q1(8)*EnR(3,3)+p%Q1(9)*EnR(4,1)+p%Q1(10)*EnR(4,2)+p%Q1(11)*EnR(4,3)
            Q2=p%Q2(1)*d1r3+p%Q2(2)*d1r1+p%Q2(3)+p%Q2(4)*EnR(2,3)+p%Q2(5)*EnR(3,1)+p%Q2(6)*EnR(3,2)+p%Q2(7)*EnR(3,3)&
              +p%Q2(8)*EnR(4,1)+p%Q2(9)*EnR(4,2)+p%Q2(10)*EnR(4,3)
            S01=p%S01(1)*d1r2+p%S01(2)+p%S01(3)*EnR(2,1)+p%S01(4)*EnR(2,2)+p%S01(5)*EnR(2,3)+p%S01(6)*EnR(3,1)&
               +p%S01(7)*EnR(3,2)+p%S01(8)*EnR(3,3)+p%S01(9)*EnR(4,1)+p%S01(10)*EnR(4,2)+p%S01(11)*EnR(4,3)&
               +p%S01(12)*EnR(5,1)+p%S01(13)*EnR(5,2)+p%S01(14)*EnR(5,3)
            S02=p%S02(1)*d1r2+p%S02(2)+p%S02(3)*EnR(2,1)+p%S02(4)*EnR(2,2)+p%S02(5)*EnR(3,1)+p%S02(6)*EnR(3,2)&
               +p%S02(7)*EnR(3,3)+p%S02(8)*EnR(4,1)+p%S02(9)*EnR(4,2)+p%S02(10)*EnR(4,3)+p%S02(11)*EnR(5,1)&
               +p%S02(12)*EnR(5,2)+p%S02(13)*EnR(5,3)
            S03=p%S03(1)*d1r2+p%S03(2)+p%S03(3)*EnR(2,3)+p%S03(4)*EnR(3,1)+p%S03(5)*EnR(3,2)+p%S03(6)*EnR(3,3)&
               +p%S03(7)*EnR(4,1)+p%S03(8)*EnR(4,2)+p%S03(9)*EnR(4,3)+p%S03(10)*EnR(5,1)+p%S03(11)*EnR(5,2)&
               +p%S03(12)*EnR(5,3)
            S1=p%S1(1)*d1r3+p%S1(2)*d1r1+p%S1(3)*EnR(2,3)+p%S1(4)*EnR(3,1)+p%S1(5)*EnR(3,2)+p%S1(6)*EnR(3,3)&
              +p%S1(7)*EnR(4,1)+p%S1(8)*EnR(4,2)+p%S1(9)*EnR(4,3)+p%S1(10)*EnR(5,1)+p%S1(11)*EnR(5,2)+p%S1(12)*EnR(5,3)&
              +p%S1(13)*EnR(6,1)+p%S1(14)*EnR(6,2)+p%S1(15)*EnR(6,3)
            S2=p%S2(1)*d1r3+p%S2(2)*d1r1+p%S2(3)*EnR(2,1)+p%S2(4)*EnR(2,2)+p%S2(5)*EnR(3,1)+p%S2(6)*EnR(3,2)&
              +p%S2(7)*EnR(3,3)+p%S2(8)*EnR(4,1)+p%S2(9)*EnR(4,2)+p%S2(10)*EnR(4,3)+p%S2(11)*EnR(5,1)+p%S2(12)*EnR(5,2)&
              +p%S2(13)*EnR(5,3)+p%S2(14)*EnR(6,1)+p%S2(15)*EnR(6,2)+p%S2(16)*EnR(6,3)
            S3=p%S3(1)*d1r3+p%S3(2)*d1r1+p%S3(3)*EnR(2,1)+p%S3(4)*EnR(2,2)+p%S3(5)*EnR(2,3)+p%S3(6)*EnR(3,1)&
              +p%S3(7)*EnR(3,2)+p%S3(8)*EnR(3,3)+p%S3(9)*EnR(4,1)+p%S3(10)*EnR(4,2)+p%S3(11)*EnR(4,3)+p%S3(12)*EnR(5,1)&
              +p%S3(13)*EnR(5,2)+p%S3(14)*EnR(5,3)+p%S3(15)*EnR(6,1)+p%S3(16)*EnR(6,2)+p%S3(17)*EnR(6,3)
            S4=p%S4(1)*d1r3+p%S4(2)*d1r1+p%S4(3)+p%S4(4)*EnR(3,3)+p%S4(5)*EnR(4,1)+p%S4(6)*EnR(4,2)+p%S4(7)*EnR(4,3)&
              +p%S4(8)*EnR(5,1)+p%S4(9)*EnR(5,2)+p%S4(10)*EnR(5,3)+p%S4(11)*EnR(6,1)+p%S4(12)*EnR(6,2)+p%S4(13)*EnR(6,3)
            S5=p%S5(1)*d1r3+p%S5(2)*d1r1+p%S5(3)+p%S5(4)*EnR(2,1)+p%S5(5)*EnR(2,2)+p%S5(6)*EnR(3,1)+p%S5(7)*EnR(3,2)&
              +p%S5(8)*EnR(4,1)+p%S5(9)*EnR(4,2)+p%S5(10)*EnR(4,3)+p%S5(11)*EnR(5,1)+p%S5(12)*EnR(5,2)+p%S5(13)*EnR(5,3)&
              +p%S5(14)*EnR(6,1)+p%S5(15)*EnR(6,2)+p%S5(16)*EnR(6,3)
            ! Fundamental solutions
            fs_d(0,0)=W0*drdni
            fs_s(0,0)=Q1*drdn*drdni+Q2*n_dot_ni
            fs_d(0,1:3)=-W1*drdx*drdni+W2*n_i
            fs_s(0,1:3)=S01*drdx*drdn*drdni+S02*n*drdni+S03*(n_i*drdn+drdx*n_dot_ni)
            fs_d(1:3,0)=-T01*drdx*drdni+T02*n_i
            fs_s(1:3,0)=-S01*drdx*drdn*drdni+S02*n_i*drdn-S03*(-n*drdni+drdx*n_dot_ni)
            do il=1,3
              do ik=1,3
                fs_d(il,ik)=TT1*drdx(il)*drdx(ik)*drdni-TT2*(drdx(il)*n_i(ik)-c_dkr(il,ik)*drdni)-TT3*drdx(ik)*n_i(il)
                fs_s(il,ik)=S1*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni-c_dkr(il,ik)*drdn*drdni+drdx(il)*drdx(ik)*n_dot_ni)&
                           +S2*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni)+S3*drdx(il)*drdx(ik)*drdn*drdni&
                           +S4*(c_dkr(il,ik)*n_dot_ni+n_i(ik)*n(il))+S5*n(ik)*n_i(il)
              end do
            end do
            ! Add the term of the numerical integration summation (fundamental solutions * functional shape functions * jacobian * weight)
            do il=0,3
              do ik=0,3
                m(:,il,ik)=m(:,il,ik)+fs_s(il,ik)*pphijw
                l(:,il,ik)=l(:,il,ik)+fs_d(il,ik)*sphijw
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
            z(3)=-c_im*p%k3*r
            call fbem_zexp_decomposed(3,z,EnR)
            EnR(2,:)=EnR(2,:)*d1r1
            EnR(3,:)=EnR(3,:)*d1r2
            EnR(4,:)=EnR(4,:)*d1r3
            EnR(5,:)=EnR(5,:)*d1r4
            EnR(6,:)=EnR(6,:)*d1r5
            W0=p%W0(1)*d1r2+p%W0(2)+p%W0(3)*EnR(2,1)+p%W0(4)*EnR(2,2)+p%W0(5)*EnR(3,1)+p%W0(6)*EnR(3,2)
            T01=p%T01(1)*d1r1+p%T01(2)*EnR(2,1)+p%T01(3)*EnR(2,2)+p%T01(4)*EnR(3,1)+p%T01(5)*EnR(3,2)+p%T01(6)*EnR(4,1)&
               +p%T01(7)*EnR(4,2)
            T02=p%T02(1)*d1r1+p%T02(2)+p%T02(3)*EnR(2,1)+p%T02(4)*EnR(2,2)+p%T02(5)*EnR(3,1)+p%T02(6)*EnR(3,2)&
               +p%T02(7)*EnR(4,1)+p%T02(8)*EnR(4,2)
            W1=p%W1(1)*d1r1+p%W1(2)*EnR(2,1)+p%W1(3)*EnR(2,2)+p%W1(4)*EnR(2,3)+p%W1(5)*EnR(3,1)+p%W1(6)*EnR(3,2)&
              +p%W1(7)*EnR(3,3)+p%W1(8)*EnR(4,1)+p%W1(9)*EnR(4,2)+p%W1(10)*EnR(4,3)
            W2=p%W2(1)*d1r1+p%W2(2)+p%W2(3)*EnR(2,3)+p%W2(4)*EnR(3,1)+p%W2(5)*EnR(3,2)+p%W2(6)*EnR(3,3)+p%W2(7)*EnR(4,1)&
              +p%W2(8)*EnR(4,2)+p%W2(9)*EnR(4,3)
            TT1=p%T1(1)*d1r2+p%T1(2)+p%T1(3)*EnR(2,1)+p%T1(4)*EnR(2,2)+p%T1(5)*EnR(2,3)+p%T1(6)*EnR(3,1)+p%T1(7)*EnR(3,2)&
               +p%T1(8)*EnR(3,3)+p%T1(9)*EnR(4,1)+p%T1(10)*EnR(4,2)+p%T1(11)*EnR(4,3)+p%T1(12)*EnR(5,1)+p%T1(13)*EnR(5,2)&
               +p%T1(14)*EnR(5,3)
            TT2=p%T2(1)*d1r2+p%T2(2)+p%T2(3)*EnR(2,3)+p%T2(4)*EnR(3,1)+p%T2(5)*EnR(3,2)+p%T2(6)*EnR(3,3)+p%T2(7)*EnR(4,1)&
               +p%T2(8)*EnR(4,2)+p%T2(9)*EnR(4,3)+p%T2(10)*EnR(5,1)+p%T2(11)*EnR(5,2)+p%T2(12)*EnR(5,3)
            TT3=p%T3(1)*d1r2+p%T3(2)+p%T3(3)*EnR(2,1)+p%T3(4)*EnR(2,2)+p%T3(5)*EnR(3,1)+p%T3(6)*EnR(3,2)+p%T3(7)*EnR(3,3)&
               +p%T3(8)*EnR(4,1)+p%T3(9)*EnR(4,2)+p%T3(10)*EnR(4,3)+p%T3(11)*EnR(5,1)+p%T3(12)*EnR(5,2)+p%T3(13)*EnR(5,3)
            Q1=p%Q1(1)*d1r3+p%Q1(2)*d1r1+p%Q1(3)*EnR(2,1)+p%Q1(4)*EnR(2,2)+p%Q1(5)*EnR(2,3)+p%Q1(6)*EnR(3,1)&
              +p%Q1(7)*EnR(3,2)+p%Q1(8)*EnR(3,3)+p%Q1(9)*EnR(4,1)+p%Q1(10)*EnR(4,2)+p%Q1(11)*EnR(4,3)
            Q2=p%Q2(1)*d1r3+p%Q2(2)*d1r1+p%Q2(3)+p%Q2(4)*EnR(2,3)+p%Q2(5)*EnR(3,1)+p%Q2(6)*EnR(3,2)+p%Q2(7)*EnR(3,3)&
              +p%Q2(8)*EnR(4,1)+p%Q2(9)*EnR(4,2)+p%Q2(10)*EnR(4,3)
            S01=p%S01(1)*d1r2+p%S01(2)+p%S01(3)*EnR(2,1)+p%S01(4)*EnR(2,2)+p%S01(5)*EnR(2,3)+p%S01(6)*EnR(3,1)&
               +p%S01(7)*EnR(3,2)+p%S01(8)*EnR(3,3)+p%S01(9)*EnR(4,1)+p%S01(10)*EnR(4,2)+p%S01(11)*EnR(4,3)&
               +p%S01(12)*EnR(5,1)+p%S01(13)*EnR(5,2)+p%S01(14)*EnR(5,3)
            S02=p%S02(1)*d1r2+p%S02(2)+p%S02(3)*EnR(2,1)+p%S02(4)*EnR(2,2)+p%S02(5)*EnR(3,1)+p%S02(6)*EnR(3,2)&
               +p%S02(7)*EnR(3,3)+p%S02(8)*EnR(4,1)+p%S02(9)*EnR(4,2)+p%S02(10)*EnR(4,3)+p%S02(11)*EnR(5,1)&
               +p%S02(12)*EnR(5,2)+p%S02(13)*EnR(5,3)
            S03=p%S03(1)*d1r2+p%S03(2)+p%S03(3)*EnR(2,3)+p%S03(4)*EnR(3,1)+p%S03(5)*EnR(3,2)+p%S03(6)*EnR(3,3)&
               +p%S03(7)*EnR(4,1)+p%S03(8)*EnR(4,2)+p%S03(9)*EnR(4,3)+p%S03(10)*EnR(5,1)+p%S03(11)*EnR(5,2)&
               +p%S03(12)*EnR(5,3)
            S1=p%S1(1)*d1r3+p%S1(2)*d1r1+p%S1(3)*EnR(2,3)+p%S1(4)*EnR(3,1)+p%S1(5)*EnR(3,2)+p%S1(6)*EnR(3,3)&
              +p%S1(7)*EnR(4,1)+p%S1(8)*EnR(4,2)+p%S1(9)*EnR(4,3)+p%S1(10)*EnR(5,1)+p%S1(11)*EnR(5,2)+p%S1(12)*EnR(5,3)&
              +p%S1(13)*EnR(6,1)+p%S1(14)*EnR(6,2)+p%S1(15)*EnR(6,3)
            S2=p%S2(1)*d1r3+p%S2(2)*d1r1+p%S2(3)*EnR(2,1)+p%S2(4)*EnR(2,2)+p%S2(5)*EnR(3,1)+p%S2(6)*EnR(3,2)&
              +p%S2(7)*EnR(3,3)+p%S2(8)*EnR(4,1)+p%S2(9)*EnR(4,2)+p%S2(10)*EnR(4,3)+p%S2(11)*EnR(5,1)+p%S2(12)*EnR(5,2)&
              +p%S2(13)*EnR(5,3)+p%S2(14)*EnR(6,1)+p%S2(15)*EnR(6,2)+p%S2(16)*EnR(6,3)
            S3=p%S3(1)*d1r3+p%S3(2)*d1r1+p%S3(3)*EnR(2,1)+p%S3(4)*EnR(2,2)+p%S3(5)*EnR(2,3)+p%S3(6)*EnR(3,1)&
              +p%S3(7)*EnR(3,2)+p%S3(8)*EnR(3,3)+p%S3(9)*EnR(4,1)+p%S3(10)*EnR(4,2)+p%S3(11)*EnR(4,3)+p%S3(12)*EnR(5,1)&
              +p%S3(13)*EnR(5,2)+p%S3(14)*EnR(5,3)+p%S3(15)*EnR(6,1)+p%S3(16)*EnR(6,2)+p%S3(17)*EnR(6,3)
            S4=p%S4(1)*d1r3+p%S4(2)*d1r1+p%S4(3)+p%S4(4)*EnR(3,3)+p%S4(5)*EnR(4,1)+p%S4(6)*EnR(4,2)+p%S4(7)*EnR(4,3)&
              +p%S4(8)*EnR(5,1)+p%S4(9)*EnR(5,2)+p%S4(10)*EnR(5,3)+p%S4(11)*EnR(6,1)+p%S4(12)*EnR(6,2)+p%S4(13)*EnR(6,3)
            S5=p%S5(1)*d1r3+p%S5(2)*d1r1+p%S5(3)+p%S5(4)*EnR(2,1)+p%S5(5)*EnR(2,2)+p%S5(6)*EnR(3,1)+p%S5(7)*EnR(3,2)&
              +p%S5(8)*EnR(4,1)+p%S5(9)*EnR(4,2)+p%S5(10)*EnR(4,3)+p%S5(11)*EnR(5,1)+p%S5(12)*EnR(5,2)+p%S5(13)*EnR(5,3)&
              +p%S5(14)*EnR(6,1)+p%S5(15)*EnR(6,2)+p%S5(16)*EnR(6,3)
            ! Fundamental solutions
            fs_d(0,0)=W0*drdni
            fs_s(0,0)=Q1*drdn*drdni+Q2*n_dot_ni
            fs_d(0,1:3)=-W1*drdx*drdni+W2*n_i
            fs_s(0,1:3)=S01*drdx*drdn*drdni+S02*n*drdni+S03*(n_i*drdn+drdx*n_dot_ni)
            fs_d(1:3,0)=-T01*drdx*drdni+T02*n_i
            fs_s(1:3,0)=-S01*drdx*drdn*drdni+S02*n_i*drdn-S03*(-n*drdni+drdx*n_dot_ni)
            do il=1,3
              do ik=1,3
                fs_d(il,ik)=TT1*drdx(il)*drdx(ik)*drdni-TT2*(drdx(il)*n_i(ik)-c_dkr(il,ik)*drdni)-TT3*drdx(ik)*n_i(il)
                fs_s(il,ik)=S1*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni-c_dkr(il,ik)*drdn*drdni+drdx(il)*drdx(ik)*n_dot_ni)&
                           +S2*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni)+S3*drdx(il)*drdx(ik)*drdn*drdni&
                           +S4*(c_dkr(il,ik)*n_dot_ni+n_i(ik)*n(il))+S5*n(ik)*n_i(il)
              end do
            end do
            ! Add the term of the numerical integration summation (fundamental solutions * functional shape functions * jacobian * weight)
            do il=0,3
              do ik=0,3
                m(:,il,ik)=m(:,il,ik)+fs_s(il,ik)*pphijw
                l(:,il,ik)=l(:,il,ik)+fs_d(il,ik)*sphijw
              end do
            end do
          end do
        end do
        case default
          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'n_edges not valid')
    end select
    ! Multiply by constants
    cte_d(0,0)=-c_1_4pi/p%J
    cte_s(0,0)=-c_1_4pi
    cte_d(0,1:3)=c_1_4pi/p%mu
    cte_s(0,1:3)=c_1_4pi
    cte_d(1:3,0)=-c_1_4pi/p%J
    cte_s(1:3,0)=-c_1_4pi
    cte_d(1:3,1:3)=c_1_4pi
    cte_s(1:3,1:3)=c_1_4pi*p%mu
    do il=0,3
      do ik=0,3
        m(:,il,ik)=cte_s(il,ik)*m(:,il,ik)
        l(:,il,ik)=cte_d(il,ik)*l(:,il,ik)
      end do
    end do
    ! Reverse if needed
    if (reverse) m=-m
  end subroutine fbem_bem_harpor3d_hbie_ext_st

  recursive subroutine fbem_bem_harpor3d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,p,qsp,ks,ns,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                                !! Element
    logical                            :: reverse                          !! Reverse orientation
    real(kind=real64)                  :: xi_s(2,fbem_n_vertices(e%gtype)) !! Subdivision of the parent element
    real(kind=real64)                  :: x_i(3)                           !! Collocation point position vector
    real(kind=real64)                  :: n_i(3)                           !! Unit normal at the collocation point
    type(fbem_bem_harpor3d_parameters) :: p                                !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp                              !! Quasi-singular integration parameters
    integer                            :: ks                               !! Current level of subdivisions
    integer                            :: ns                               !! Maximum level of subdivision
    complex(kind=real64)               :: m(e%n_pnodes,0:3,0:3)            !! m integration kernels matrix
    complex(kind=real64)               :: l(e%n_snodes,0:3,0:3)            !! l integration kernels matrix
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
    complex(kind=real64) :: m_tmp(e%n_pnodes,0:3,0:3)            ! m integration kernels matrix (temporary)
    complex(kind=real64) :: l_tmp(e%n_snodes,0:3,0:3)            ! l integration kernels matrix (temporary)
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
        call fbem_warning_message(error_unit,0,'fbem_bem_harpor3d_hbie_ext_adp',ns,'maximum number of subdivisions reached')
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
          call fbem_bem_harpor3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
          ! SUBTRI 2
          tmp_xi_s(:,1)=xi_s(:,2)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          call fbem_bem_harpor3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
          ! SUBTRI 3
          tmp_xi_s(:,1)=xi_s(:,3)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          call fbem_bem_harpor3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
          ! SUBTRI 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,1)+xi_s(:,3))
          call fbem_bem_harpor3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
        ! QUADRILATERALS
        case (4)
          ! SUBQUAD 1
          tmp_xi_s(:,1)=xi_s(:,1)
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,3)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          call fbem_bem_harpor3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
          ! SUBQUAD 2
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,2))
          tmp_xi_s(:,2)=xi_s(:,2)
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,4)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          call fbem_bem_harpor3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
          ! SUBQUAD 3
          tmp_xi_s(:,1)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,2)=0.50d0*(xi_s(:,2)+xi_s(:,3))
          tmp_xi_s(:,3)=xi_s(:,3)
          tmp_xi_s(:,4)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          call fbem_bem_harpor3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
          ! SUBQUAD 4
          tmp_xi_s(:,1)=0.50d0*(xi_s(:,1)+xi_s(:,4))
          tmp_xi_s(:,2)=0.25d0*(xi_s(:,1)+xi_s(:,2)+xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,3)=0.50d0*(xi_s(:,3)+xi_s(:,4))
          tmp_xi_s(:,4)=xi_s(:,4)
          call fbem_bem_harpor3d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
      end select
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_harpor3d_hbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,p,gln,m_tmp,l_tmp)
      m=m+m_tmp
      l=l+l_tmp
    end if
  end subroutine fbem_bem_harpor3d_hbie_ext_adp

  subroutine fbem_bem_harpor3d_hbie_int(type_g,type_f1,type_f2,delta_f,x_nodes,reverse,xi_i,p,m,l)
    implicit none
    ! I/O
    integer                            :: type_g                           !! Geometrial interpolation
    integer                            :: type_f1                          !! Functional interpolation (primary variables)
    integer                            :: type_f2                          !! Functional interpolation (secondary variables)
    real(kind=real64)                  :: delta_f                          !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64)                  :: x_nodes(3,fbem_n_nodes(type_g))  !! Position vectors of geometrical nodes
    logical                            :: reverse                          !! Normal vector inversion
    real(kind=real64)                  :: xi_i(2)                          !! Reference coordinates of the singular point.
    type(fbem_bem_harpor3d_parameters) :: p                                !! Parameters of the region
    complex(kind=real64)               :: m(fbem_n_nodes(type_f1),0:3,0:3) !! m kernel vector
    complex(kind=real64)               :: l(fbem_n_nodes(type_f2),0:3,0:3) !! l kernel vector
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
    integer              :: nnodes_f1                        ! Number of nodes of primary variables interpolation
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
    real(kind=real64)    :: phi_f1_i(fbem_n_nodes(type_f1))  ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)    :: psi_i(3,fbem_n_nodes(type_f1))   ! psi vectors at xi_i1,xi_i2
    real(kind=real64)    :: phi_f2_i(fbem_n_nodes(type_f2))  ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)    :: phi_g(fbem_n_nodes(type_g))      ! Geometrical shape functions values at xi_1,xi_2
    real(kind=real64)    :: dphidxi1_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64)    :: dphidxi2_g(fbem_n_nodes(type_g)) ! Geometrical shape functions derivatives values at xi_1,xi_2
    real(kind=real64)    :: jxi1, jxi2                       ! xi1->x, xi2->x jacobians
    real(kind=real64)    :: x_i(3)                           ! Collocation point position vector
    real(kind=real64)    :: n_i(3)                           ! Collocation point unit normal vector
    real(kind=real64)    :: x(3)                             ! Position vector at xi_1,xi_2
    real(kind=real64)    :: T1(3), T2(3)                     ! Tangent vectors at xi_1,xi_2
    real(kind=real64)    :: v1(3), v2(3)                     ! Orthogonal tangent vectors
    real(kind=real64)    :: Mvt(2,2)                         ! Orthogonal coordinates transformation matrix
    real(kind=real64)    :: dxidv(2,2)                       ! xi coordinates derivatives with respect to v orthogonal coordinates
    real(kind=real64)    :: detMvt                           ! Mvt(1,1)*Mvt(2,2)-Mvt(1,2)*Mvt(2,1)
    real(kind=real64)    :: N(3)                             ! Normal vector at xi_1,xi_2
    real(kind=real64)    :: rv(3)                            ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, d1r1, d1r2, d1r3, d1r4, d1r5  ! Distance vector module and its inverse
    real(kind=real64)    :: drdx(3)                          ! Distance vector derivatives with respect to x_k
    real(kind=real64)    :: jg                               ! Geometric jacobian
    real(kind=real64)    :: jw                               ! Jacobian * weights
    real(kind=real64)    :: drdn                             ! Partial derivative of r respect to unit normal
    real(kind=real64)    :: drdni                            ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)    :: n_dot_ni                         ! Dot product n · n_i
    real(kind=real64)    :: phif1jw(fbem_n_nodes(type_f1))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)    :: phif2jw(fbem_n_nodes(type_f2))   ! Functional shape functions values at xi_1,xi_2
    real(kind=real64)    :: costheta, sintheta               ! cos(theta), sin(theta)
    integer              :: il, ik, ij                       ! Counter for load / observation components
    complex(kind=real64) :: z1, z2, z3                                ! zj=-i*kj*r
    complex(kind=real64) :: E0z1, E1z1, E2z1, E3z1, E4z1, E5z1, E6z1  ! exp(z1) decomposition
    complex(kind=real64) :: E0z2, E1z2, E2z2, E3z2, E4z2, E5z2, E6z2  ! exp(z2) decomposition
    complex(kind=real64) :: E0z3, E1z3, E2z3, E3z3, E4z3, E5z3, E6z3  ! exp(z3) decomposition
    complex(kind=real64) :: W0, T01, T02, W1, W2, TT1, TT2, TT3       ! Components of the fundamental solution D*
    complex(kind=real64) :: Q1, Q2, S01, S02, S03, S1, S2, S3, S4, S5 ! Components of the fundamental solution S*
    complex(kind=real64) :: fs_d(0:3,0:3), fs_s(0:3,0:3), auxs        ! Fundamental solutions matrices
    complex(kind=real64) :: cte_d(0:3,0:3), cte_s(0:3,0:3)            ! Constants of fundamental solutions matrices
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
      ngp_rho=32!15
      ngp_theta=32!5+nint(25.0d0*(thetaf-thetai)/c_pi_2)
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
          d1r1=1.d0/r
          d1r2=d1r1**2
          d1r3=d1r2*d1r1
          d1r4=d1r3*d1r1
          d1r5=d1r4*d1r1
          ! r_{,k}
          drdx=rv*d1r1
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
          !
          ! ONLY AT MOST WEAKLY SINGULAR PARTS
          !
          ! Calculation of the COMPONENTS of the FUNDAMENTAL SOLUTIONS
          z1=-c_im*p%k1*r
          z2=-c_im*p%k2*r
          z3=-c_im*p%k3*r
          call fbem_decomposed_zexp_6(z1,E0z1,E1z1,E2z1,E3z1,E4z1,E5z1,E6z1)
          call fbem_decomposed_zexp_6(z2,E0z2,E1z2,E2z2,E3z2,E4z2,E5z2,E6z2)
          call fbem_decomposed_zexp_6(z3,E0z3,E1z3,E2z3,E3z3,E4z3,E5z3,E6z3)
          ! Components of D*
          W0=p%W0(2)+p%W0(3)*d1r1*E2z1+p%W0(4)*d1r1*E2z2+p%W0(5)*d1r2*E3z1+p%W0(6)*d1r2*E3z2
          T01=p%T01(1)*d1r1+p%T01(2)*d1r1*E2z1+p%T01(3)*d1r1*E2z2+p%T01(4)*d1r2*E3z1+p%T01(5)*d1r2*E3z2+p%T01(6)*d1r3*E4z1&
             +p%T01(7)*d1r3*E4z2
          T02=p%T02(1)*d1r1+p%T02(2)+p%T02(3)*d1r1*E2z1+p%T02(4)*d1r1*E2z2+p%T02(5)*d1r2*E3z1+p%T02(6)*d1r2*E3z2&
             +p%T02(7)*d1r3*E4z1+p%T02(8)*d1r3*E4z2
          W1=p%W1(1)*d1r1+p%W1(2)*d1r1*E2z1+p%W1(3)*d1r1*E2z2+p%W1(4)*d1r1*E2z3+p%W1(5)*d1r2*E3z1+p%W1(6)*d1r2*E3z2&
            +p%W1(7)*d1r2*E3z3+p%W1(8)*d1r3*E4z1+p%W1(9)*d1r3*E4z2+p%W1(10)*d1r3*E4z3
          W2=p%W2(1)*d1r1+p%W2(2)+p%W2(3)*d1r1*E2z3+p%W2(4)*d1r2*E3z1+p%W2(5)*d1r2*E3z2+p%W2(6)*d1r2*E3z3+p%W2(7)*d1r3*E4z1&
            +p%W2(8)*d1r3*E4z2+p%W2(9)*d1r3*E4z3
          TT1=p%T1(2)+p%T1(3)*d1r1*E2z1+p%T1(4)*d1r1*E2z2+p%T1(5)*d1r1*E2z3+p%T1(6)*d1r2*E3z1+p%T1(7)*d1r2*E3z2&
             +p%T1(8)*d1r2*E3z3+p%T1(9)*d1r3*E4z1+p%T1(10)*d1r3*E4z2+p%T1(11)*d1r3*E4z3+p%T1(12)*d1r4*E5z1+p%T1(13)*d1r4*E5z2&
             +p%T1(14)*d1r4*E5z3
          TT2=p%T2(2)+p%T2(3)*d1r1*E2z3+p%T2(4)*d1r2*E3z1+p%T2(5)*d1r2*E3z2+p%T2(6)*d1r2*E3z3+p%T2(7)*d1r3*E4z1&
             +p%T2(8)*d1r3*E4z2+p%T2(9)*d1r3*E4z3+p%T2(10)*d1r4*E5z1+p%T2(11)*d1r4*E5z2+p%T2(12)*d1r4*E5z3
          TT3=p%T3(2)+p%T3(3)*d1r1*E2z1+p%T3(4)*d1r1*E2z2+p%T3(5)*d1r2*E3z1+p%T3(6)*d1r2*E3z2+p%T3(7)*d1r2*E3z3&
             +p%T3(8)*d1r3*E4z1+p%T3(9)*d1r3*E4z2+p%T3(10)*d1r3*E4z3+p%T3(11)*d1r4*E5z1+p%T3(12)*d1r4*E5z2+p%T3(13)*d1r4*E5z3
          ! Components of S*
          Q1=p%Q1(2)*d1r1+p%Q1(3)*d1r1*E2z1+p%Q1(4)*d1r1*E2z2+p%Q1(5)*d1r1*E2z3+p%Q1(6)*d1r2*E3z1&
            +p%Q1(7)*d1r2*E3z2+p%Q1(8)*d1r2*E3z3+p%Q1(9)*d1r3*E4z1+p%Q1(10)*d1r3*E4z2+p%Q1(11)*d1r3*E4z3
          Q2=p%Q2(2)*d1r1+p%Q2(3)+p%Q2(4)*d1r1*E2z3+p%Q2(5)*d1r2*E3z1+p%Q2(6)*d1r2*E3z2+p%Q2(7)*d1r2*E3z3&
            +p%Q2(8)*d1r3*E4z1+p%Q2(9)*d1r3*E4z2+p%Q2(10)*d1r3*E4z3
          S01=p%S01(2)+p%S01(3)*d1r1*E2z1+p%S01(4)*d1r1*E2z2+p%S01(5)*d1r1*E2z3+p%S01(6)*d1r2*E3z1&
             +p%S01(7)*d1r2*E3z2+p%S01(8)*d1r2*E3z3+p%S01(9)*d1r3*E4z1+p%S01(10)*d1r3*E4z2+p%S01(11)*d1r3*E4z3&
             +p%S01(12)*d1r4*E5z1+p%S01(13)*d1r4*E5z2+p%S01(14)*d1r4*E5z3
          S02=p%S02(2)+p%S02(3)*d1r1*E2z1+p%S02(4)*d1r1*E2z2+p%S02(5)*d1r2*E3z1+p%S02(6)*d1r2*E3z2&
             +p%S02(7)*d1r2*E3z3+p%S02(8)*d1r3*E4z1+p%S02(9)*d1r3*E4z2+p%S02(10)*d1r3*E4z3+p%S02(11)*d1r4*E5z1&
             +p%S02(12)*d1r4*E5z2+p%S02(13)*d1r4*E5z3
          S03=p%S03(2)+p%S03(3)*d1r1*E2z3+p%S03(4)*d1r2*E3z1+p%S03(5)*d1r2*E3z2+p%S03(6)*d1r2*E3z3&
             +p%S03(7)*d1r3*E4z1+p%S03(8)*d1r3*E4z2+p%S03(9)*d1r3*E4z3+p%S03(10)*d1r4*E5z1+p%S03(11)*d1r4*E5z2&
             +p%S03(12)*d1r4*E5z3
          S1=p%S1(2)*d1r1+p%S1(3)*d1r1*E2z3+p%S1(4)*d1r2*E3z1+p%S1(5)*d1r2*E3z2+p%S1(6)*d1r2*E3z3&
            +p%S1(7)*d1r3*E4z1+p%S1(8)*d1r3*E4z2+p%S1(9)*d1r3*E4z3+p%S1(10)*d1r4*E5z1+p%S1(11)*d1r4*E5z2+p%S1(12)*d1r4*E5z3&
            +p%S1(13)*d1r5*E6z1+p%S1(14)*d1r5*E6z2+p%S1(15)*d1r5*E6z3
          S2=p%S2(2)*d1r1+p%S2(3)*d1r1*E2z1+p%S2(4)*d1r1*E2z2+p%S2(5)*d1r2*E3z1+p%S2(6)*d1r2*E3z2&
            +p%S2(7)*d1r2*E3z3+p%S2(8)*d1r3*E4z1+p%S2(9)*d1r3*E4z2+p%S2(10)*d1r3*E4z3+p%S2(11)*d1r4*E5z1+p%S2(12)*d1r4*E5z2&
            +p%S2(13)*d1r4*E5z3+p%S2(14)*d1r5*E6z1+p%S2(15)*d1r5*E6z2+p%S2(16)*d1r5*E6z3
          S3=p%S3(2)*d1r1+p%S3(3)*d1r1*E2z1+p%S3(4)*d1r1*E2z2+p%S3(5)*d1r1*E2z3+p%S3(6)*d1r2*E3z1&
            +p%S3(7)*d1r2*E3z2+p%S3(8)*d1r2*E3z3+p%S3(9)*d1r3*E4z1+p%S3(10)*d1r3*E4z2+p%S3(11)*d1r3*E4z3+p%S3(12)*d1r4*E5z1&
            +p%S3(13)*d1r4*E5z2+p%S3(14)*d1r4*E5z3+p%S3(15)*d1r5*E6z1+p%S3(16)*d1r5*E6z2+p%S3(17)*d1r5*E6z3
          S4=p%S4(2)*d1r1+p%S4(3)+p%S4(4)*d1r2*E3z3+p%S4(5)*d1r3*E4z1+p%S4(6)*d1r3*E4z2+p%S4(7)*d1r3*E4z3&
            +p%S4(8)*d1r4*E5z1+p%S4(9)*d1r4*E5z2+p%S4(10)*d1r4*E5z3+p%S4(11)*d1r5*E6z1+p%S4(12)*d1r5*E6z2+p%S4(13)*d1r5*E6z3
          S5=p%S5(2)*d1r1+p%S5(3)+p%S5(4)*d1r1*E2z1+p%S5(5)*d1r1*E2z2+p%S5(6)*d1r2*E3z1+p%S5(7)*d1r2*E3z2&
            +p%S5(8)*d1r3*E4z1+p%S5(9)*d1r3*E4z2+p%S5(10)*d1r3*E4z3+p%S5(11)*d1r4*E5z1+p%S5(12)*d1r4*E5z2+p%S5(13)*d1r4*E5z3&
            +p%S5(14)*d1r5*E6z1+p%S5(15)*d1r5*E6z2+p%S5(16)*d1r5*E6z3
          ! FUNDAMENTAL SOLUTIONS
          ! Fluid load / Fluid response
          fs_d(0,0)=(p%W0(1)*d1r2+W0)*drdni
          fs_s(0,0)=Q1*drdn*drdni+Q2*n_dot_ni
          ! Fluid load / Solid response
          fs_d(0,1:3)=-W1*drdx*drdni+W2*n_i
          fs_s(0,1:3)=(p%S01(1)*d1r2+S01)*drdx*drdn*drdni+(p%S02(1)*d1r2+S02)*n*drdni+(p%S03(1)*d1r2+S03)*n_i*drdn+S03*drdx*n_dot_ni
          ! Solid load / Fluid response
          fs_d(1:3,0)=-T01*drdx*drdni+T02*n_i
          fs_s(1:3,0)=-(p%S01(1)*d1r2+S01)*drdx*drdn*drdni+(p%S02(1)*d1r2+S02)*n_i*drdn+(p%S03(1)*d1r2+S03)*n*drdni-S03*drdx*n_dot_ni
          ! Solid load / Solid response
          do il=1,3
            do ik=1,3
              fs_d(il,ik)=(p%T1(1)*d1r2+TT1)*drdx(il)*drdx(ik)*drdni+(p%T2(1)*d1r2+TT2)*drdni*c_dkr(il,ik)-TT2*drdx(il)*n_i(ik)&
                         -TT3*drdx(ik)*n_i(il)+p%T2(1)*d1r2*((n_i(il)-n(il))*drdx(ik)-(n_i(ik)-n(ik))*drdx(il))
              fs_s(il,ik)=S1*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni-c_dkr(il,ik)*drdn*drdni+drdx(il)*drdx(ik)*n_dot_ni)&
                         +S2*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni)+S3*drdx(il)*drdx(ik)*drdn*drdni&
                         +S4*(c_dkr(il,ik)*n_dot_ni+n_i(ik)*n(il))+S5*n(ik)*n_i(il)-c_dkr(il,ik)*p%S1(1)*d1r3*drdn*drdni&
                         +p%S3(1)*d1r3*drdx(il)*drdx(ik)*drdn*drdni
            end do
          end do
          ! ADD INTEGRAND EVALUATION
          do il=0,3
            do ik=0,3
              m(:,il,ik)=m(:,il,ik)+fs_s(il,ik)*phif1jw
              l(:,il,ik)=l(:,il,ik)+fs_d(il,ik)*phif2jw
            end do
          end do
          !
          ! REGULAR PARTS OF CPV AND HFP INTEGRALS
          !
          ! Fluid load / Fluid response
          m(:,0,0)=m(:,0,0)+p%Q1(1)*d1r3*drdn*drdni*(phi_f1-phi_f1_i)*jw
          m(:,0,0)=m(:,0,0)+p%Q2(1)*d1r3*n_dot_ni*(phi_f1-phi_f1_i-psi_i(1,:)*rv(1)-psi_i(2,:)*rv(2)-psi_i(3,:)*rv(3))*jw
          m(:,0,0)=m(:,0,0)-p%Q2(1)*(psi_i(1,:)*n(1)+psi_i(2,:)*n(2)+psi_i(3,:)*n(3))*d1r2*drdni*jw
          ! Fluid load / Solid response
          do ik=1,3
            m(:,0,ik)=m(:,0,ik)+p%S03(1)*d1r2*drdx(ik)*n_dot_ni*(phi_f1-phi_f1_i)*jw
            m(:,0,ik)=m(:,0,ik)-phi_f1_i*p%S03(1)*d1r2*n(ik)*drdni*jw
          end do
          ! Solid load / Fluid response
          do il=1,3
            m(:,il,0)=m(:,il,0)-p%S03(1)*d1r2*drdx(il)*n_dot_ni*(phi_f1-phi_f1_i)*jw
            m(:,il,0)=m(:,il,0)+phi_f1_i*p%S03(1)*d1r2*n(il)*drdni*jw
          end do
          ! Solid load / Solid response
          do il=1,3
            do ik=1,3
              ! L21
              l(:,il,ik)=l(:,il,ik)+p%T2(1)*d1r2*(n(il)*drdx(ik)-n(ik)*drdx(il))*(phi_f2-phi_f2_i)*jw
              ! Regular parts of M2
              if (il.eq.ik) then
                ! M21
                auxs=p%S4(1)*d1r3*n_dot_ni
                m(:,il,ik)=m(:,il,ik)+auxs*(phi_f1-phi_f1_i-psi_i(1,:)*rv(1)-psi_i(2,:)*rv(2)-psi_i(3,:)*rv(3))*jw
                ! M221
                auxs=-p%S4(1)*3.d0*d1r3*drdn*drdni
                m(:,il,ik)=m(:,il,ik)+phi_f1_i*auxs*jw
                ! M231
                auxs=-p%S4(1)*d1r2*drdni
                m(:,il,ik)=m(:,il,ik)+auxs*(psi_i(1,:)*n(1)+psi_i(2,:)*n(2)+psi_i(3,:)*n(3))*jw
              end if
              ! Regular parts of M3
              ! M31
              auxs=d1r3*(p%S2(1)*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni)&
                        +p%S1(1)*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni)&
                        +p%S1(1)*drdx(il)*drdx(ik)*n_dot_ni&
                        +p%S4(1)*n_i(ik)*n(il)+p%S5(1)*n_i(il)*n(ik))
              m(:,il,ik)=m(:,il,ik)+auxs*(phi_f1-phi_f1_i-psi_i(1,:)*rv(1)-psi_i(2,:)*rv(2)-psi_i(3,:)*rv(3))*jw
              ! M321
              auxs=d1r3*(-5.d0*p%S1(1)*drdx(il)*drdx(ik)*drdn*drdni+p%S1(1)*drdx(ik)*(n_i(il)*drdn-n(il)*drdni)&
                         +p%S2(1)*drdx(il)*(n_i(ik)*drdn-n(ik)*drdni))
              m(:,il,ik)=m(:,il,ik)+phi_f1_i*auxs*jw
              ! M331
              do ij=1,3
                auxs=d1r2*(drdx(ij)*(p%S2(1)*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni)&
                                    +p%S1(1)*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni)&
                                    -p%S1(1)*drdx(il)*n(ik)*drdni)&
                          -((p%S1(1)/3.d0-p%S5(1))*n_i(il)*drdx(ij)+p%S1(1)/3.d0*n_i(ij)*drdx(il))*(n(ik)-n_i(ik)*n_dot_ni)&
                          +p%S4(1)*n_i(ik)*drdx(ij)*(n(il)-n_i(il)*n_dot_ni)&
                          -p%S1(1)/3.d0*(c_dkr(il,ik)+n_i(il)*n_i(ik))*n(ij)*drdni&
                          -p%S1(1)/3.d0*(c_dkr(ij,ik)-n_i(ij)*n_i(ik))*n(il)*drdni)
                m(:,il,ik)=m(:,il,ik)+psi_i(ij,:)*auxs*jw
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
    ! Fluid load / Fluid response
    m(:,0,0)=m(:,0,0)+phi_f1_i*p%Q2(1)*mli1
    do ij=1,3
      m(:,0,0)=m(:,0,0)+psi_i(ij,:)*p%Q2(1)*mli2(ij)
    end do
    ! Fluid load / Solid response
    do ik=1,3
      m(:,0,ik)=m(:,0,ik)+phi_f1_i*p%S03(1)*mli2(ik)
    end do
    ! Solid load / Fluid response
    do il=1,3
      m(:,il,0)=m(:,il,0)-phi_f1_i*p%S03(1)*mli2(il)
    end do
    ! Solid load / Solid response
    do il=1,3
      do ik=1,3
        !
        ! Line integrals (L)
        !
        ! L22
        l(:,il,ik)=l(:,il,ik)+phi_f2_i*p%T2(1)*lli(il,ik)
        !
        ! Line integrals (M)
        !
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
    !
    ! Calculate cte_s and cte_d, and multiply them to m and l, respectively
    !
    ! Fluid load / Fluid response constants
    cte_d(0,0)=-c_1_4pi/p%J
    cte_s(0,0)=-c_1_4pi
    ! Fluid load / Solid response constants
    cte_d(0,1:3)=c_1_4pi/p%mu
    cte_s(0,1:3)=c_1_4pi
    ! Solid load / Fluid response constants
    cte_d(1:3,0)=-c_1_4pi/p%J
    cte_s(1:3,0)=-c_1_4pi
    ! Solid load / Solid response constants
    cte_d(1:3,1:3)=c_1_4pi
    cte_s(1:3,1:3)=c_1_4pi*p%mu
    ! Multiply by constants
    do il=0,3
      do ik=0,3
        m(:,il,ik)=cte_s(il,ik)*m(:,il,ik)
        l(:,il,ik)=cte_d(il,ik)*l(:,il,ik)
      end do
    end do
    ! If the normal has to be reversed, then l=-l
    if (reverse) l=-l
  end subroutine fbem_bem_harpor3d_hbie_int

  subroutine fbem_bem_harpor3d_hbie_auto(e,reverse,x_i,n_i,p,qsp,ns,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                     !! Integration element
    logical                            :: reverse               !! Reverse orientation
    real(kind=real64)                  :: x_i(3)                !! Collocation point
    real(kind=real64)                  :: n_i(3)                !! Unit normal at the collocation point
    type(fbem_bem_harpor3d_parameters) :: p                     !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp                   !! Quasi-singular integration parameters
    integer                            :: ns                    !! Maximum level of subdivisions
    complex(kind=real64)               :: m(e%n_pnodes,0:3,0:3) !! m integration kernel
    complex(kind=real64)               :: l(e%n_snodes,0:3,0:3) !! l integration kernel
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
        call fbem_bem_harpor3d_hbie_int(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi,p,m,l)
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
          call fbem_bem_harpor3d_hbie_ext_pre(ps,e,reverse,x_i,n_i,p,m,l)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_harpor3d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,p,qsp,1,ns,m,l)
        end if
    end select
  end subroutine fbem_bem_harpor3d_hbie_auto

  ! ================================================================================================================================

  ! ================================================================================================================================
  ! BOTH SBIE AND HBIE

  subroutine fbem_bem_harpor3d_shbie_ext_pre(ps,e,reverse,x_i,n_i,p,h,g,m,l)
    implicit none
    ! I/O
    integer                            :: ps                    !! Selected precalculated dataset
    type(fbem_bem_element)             :: e                     !! Element
    logical                            :: reverse               !! Reverse normal vector
    real(kind=real64)                  :: x_i(3)                !! Collocation point position vector
    real(kind=real64)                  :: n_i(3)                !! Collocation point unit normal vector
    type(fbem_bem_harpor3d_parameters) :: p                     !! Parameters of the region
    complex(kind=real64)               :: h(e%n_pnodes,0:3,0:3) !! h integration kernels vector
    complex(kind=real64)               :: g(e%n_snodes,0:3,0:3) !! g integration kernels vector
    complex(kind=real64)               :: m(e%n_pnodes,0:3,0:3) !! m integration kernels vector
    complex(kind=real64)               :: l(e%n_snodes,0:3,0:3) !! l integration kernels vector
    ! Local
    integer              :: kip                                       ! Counter variable for integration points loop
    real(kind=real64)    :: x(3)                                      ! Position vector at integration point
    real(kind=real64)    :: n(3)                                      ! Unit normal vector at integration point
    real(kind=real64)    :: pphijw(e%n_pnodes)                        ! phi^p * jacobian * weight at integration point
    real(kind=real64)    :: sphijw(e%n_snodes)                        ! phi^s * jacobian * weight at integration point
    real(kind=real64)    :: rv(3)                                     ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, d1r1, d1r2, d1r3, d1r4, d1r5           ! Distance vector module and its inverse
    real(kind=real64)    :: drdx(3)                                   ! Distance vector derivatives with respect to x_k
    real(kind=real64)    :: drdn                                      ! Partial derivative of r respect to unit normal
    real(kind=real64)    :: drdni                                     ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)    :: n_dot_ni                                  ! Dot product of n and n_i
    integer              :: il, ik                                    ! Counter for load / observation components
    complex(kind=real64) :: z(3)                                      ! Arguments z=ikr
    complex(kind=real64) :: EnR(0:6,3)                                ! Exponential function decomposition for each wavenumber
    complex(kind=real64) :: eta, vartheta, psi, chi                   ! Components of the fundamental solution
    complex(kind=real64) :: W0, T01, T02, W1, W2, TT1, TT2, TT3       ! Components of the fundamental solution
    complex(kind=real64) :: Q1, Q2, S01, S02, S03, S1, S2, S3, S4, S5 ! Components of the fundamental solution
    complex(kind=real64) :: fs_u(0:3,0:3), fs_t(0:3,0:3)              ! Fundamental solutions matrices
    complex(kind=real64) :: fs_d(0:3,0:3), fs_s(0:3,0:3)              ! Fundamental solutions matrices
    complex(kind=real64) :: cte_u(0:3,0:3), cte_t(0:3,0:3)            ! Constants of fundamental solutions matrices
    complex(kind=real64) :: cte_d(0:3,0:3), cte_s(0:3,0:3)            ! Constants of fundamental solutions matrices
    ! Initialize
    h=(0.d0,0.d0)
    g=(0.d0,0.d0)
    m=(0.d0,0.d0)
    l=(0.d0,0.d0)
    ! Numerical integration
    do kip=1,e%ps_ngp(ps)
      ! Functional shape functions * jacobian * weight
      pphijw=e%ps_pphijw(:,kip,ps)
      sphijw=e%ps_sphijw(:,kip,ps)
      ! Calculation of the components of the fundamental solution
      x=e%ps_x(:,kip,ps)
      n=e%ps_n(:,kip,ps)
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
      z(3)=-c_im*p%k3*r
      call fbem_zexp_decomposed(3,z,EnR)
      EnR(2,:)=EnR(2,:)*d1r1
      EnR(3,:)=EnR(3,:)*d1r2
      EnR(4,:)=EnR(4,:)*d1r3
      EnR(5,:)=EnR(5,:)*d1r4
      EnR(6,:)=EnR(6,:)*d1r5
      eta=d1r1+p%eta(1)+p%eta(2)*EnR(2,1)+p%eta(3)*EnR(2,2)
      vartheta=p%vartheta(1)+p%vartheta(2)*EnR(2,1)+p%vartheta(3)*EnR(2,2)+p%vartheta(4)*EnR(3,1)+p%vartheta(5)*EnR(3,2)
      psi=p%psi(1)*d1r1+p%psi(2)+EnR(2,3)+p%psi(3)*EnR(3,1)+p%psi(4)*EnR(3,2)+p%psi(5)*EnR(3,3)+p%psi(6)*EnR(4,1)&
         +p%psi(7)*EnR(4,2)+p%psi(8)*EnR(4,3)
      chi=p%chi(1)*d1r1+p%chi(2)*EnR(2,1)+p%chi(3)*EnR(2,2)+EnR(2,3)+p%chi(4)*EnR(3,1)+p%chi(5)*EnR(3,2)+p%chi(6)*EnR(3,3)&
         +p%chi(7)*EnR(4,1)+p%chi(8)*EnR(4,2)+p%chi(9)*EnR(4,3)
      W0=p%W0(1)*d1r2+p%W0(2)+p%W0(3)*EnR(2,1)+p%W0(4)*EnR(2,2)+p%W0(5)*EnR(3,1)+p%W0(6)*EnR(3,2)
      T01=p%T01(1)*d1r1+p%T01(2)*EnR(2,1)+p%T01(3)*EnR(2,2)+p%T01(4)*EnR(3,1)+p%T01(5)*EnR(3,2)+p%T01(6)*EnR(4,1)&
         +p%T01(7)*EnR(4,2)
      T02=p%T02(1)*d1r1+p%T02(2)+p%T02(3)*EnR(2,1)+p%T02(4)*EnR(2,2)+p%T02(5)*EnR(3,1)+p%T02(6)*EnR(3,2)+p%T02(7)*EnR(4,1)&
         +p%T02(8)*EnR(4,2)
      W1 =p%W1(1)*d1r1+p%W1(2)*EnR(2,1)+p%W1(3)*EnR(2,2)+p%W1(4)*EnR(2,3)+p%W1(5)*EnR(3,1)+p%W1(6)*EnR(3,2)+p%W1(7)*EnR(3,3)&
         +p%W1(8)*EnR(4,1)+p%W1(9)*EnR(4,2)+p%W1(10)*EnR(4,3)
      W2 =p%W2(1)*d1r1+p%W2(2)+p%W2(3)*EnR(2,3)+p%W2(4)*EnR(3,1)+p%W2(5)*EnR(3,2)+p%W2(6)*EnR(3,3)+p%W2(7)*EnR(4,1)&
         +p%W2(8)*EnR(4,2)+p%W2(9)*EnR(4,3)
      TT1=p%T1(1)*d1r2+p%T1(2)+p%T1(3)*EnR(2,1)+p%T1(4)*EnR(2,2)+p%T1(5)*EnR(2,3)+p%T1(6)*EnR(3,1)+p%T1(7)*EnR(3,2)&
         +p%T1(8)*EnR(3,3)+p%T1(9)*EnR(4,1)+p%T1(10)*EnR(4,2)+p%T1(11)*EnR(4,3)+p%T1(12)*EnR(5,1)+p%T1(13)*EnR(5,2)&
         +p%T1(14)*EnR(5,3)
      TT2=p%T2(1)*d1r2+p%T2(2)+p%T2(3)*EnR(2,3)+p%T2(4)*EnR(3,1)+p%T2(5)*EnR(3,2)+p%T2(6)*EnR(3,3)+p%T2(7)*EnR(4,1)&
         +p%T2(8)*EnR(4,2)+p%T2(9)*EnR(4,3)+p%T2(10)*EnR(5,1)+p%T2(11)*EnR(5,2)+p%T2(12)*EnR(5,3)
      TT3=p%T3(1)*d1r2+p%T3(2)+p%T3(3)*EnR(2,1)+p%T3(4)*EnR(2,2)+p%T3(5)*EnR(3,1)+p%T3(6)*EnR(3,2)+p%T3(7)*EnR(3,3)&
         +p%T3(8)*EnR(4,1)+p%T3(9)*EnR(4,2)+p%T3(10)*EnR(4,3)+p%T3(11)*EnR(5,1)+p%T3(12)*EnR(5,2)+p%T3(13)*EnR(5,3)
      Q1=p%Q1(1)*d1r3+p%Q1(2)*d1r1+p%Q1(3)*EnR(2,1)+p%Q1(4)*EnR(2,2)+p%Q1(5)*EnR(2,3)+p%Q1(6)*EnR(3,1)&
        +p%Q1(7)*EnR(3,2)+p%Q1(8)*EnR(3,3)+p%Q1(9)*EnR(4,1)+p%Q1(10)*EnR(4,2)+p%Q1(11)*EnR(4,3)
      Q2=p%Q2(1)*d1r3+p%Q2(2)*d1r1+p%Q2(3)+p%Q2(4)*EnR(2,3)+p%Q2(5)*EnR(3,1)+p%Q2(6)*EnR(3,2)+p%Q2(7)*EnR(3,3)&
        +p%Q2(8)*EnR(4,1)+p%Q2(9)*EnR(4,2)+p%Q2(10)*EnR(4,3)
      S01=p%S01(1)*d1r2+p%S01(2)+p%S01(3)*EnR(2,1)+p%S01(4)*EnR(2,2)+p%S01(5)*EnR(2,3)+p%S01(6)*EnR(3,1)&
         +p%S01(7)*EnR(3,2)+p%S01(8)*EnR(3,3)+p%S01(9)*EnR(4,1)+p%S01(10)*EnR(4,2)+p%S01(11)*EnR(4,3)&
         +p%S01(12)*EnR(5,1)+p%S01(13)*EnR(5,2)+p%S01(14)*EnR(5,3)
      S02=p%S02(1)*d1r2+p%S02(2)+p%S02(3)*EnR(2,1)+p%S02(4)*EnR(2,2)+p%S02(5)*EnR(3,1)+p%S02(6)*EnR(3,2)&
         +p%S02(7)*EnR(3,3)+p%S02(8)*EnR(4,1)+p%S02(9)*EnR(4,2)+p%S02(10)*EnR(4,3)+p%S02(11)*EnR(5,1)&
         +p%S02(12)*EnR(5,2)+p%S02(13)*EnR(5,3)
      S03=p%S03(1)*d1r2+p%S03(2)+p%S03(3)*EnR(2,3)+p%S03(4)*EnR(3,1)+p%S03(5)*EnR(3,2)+p%S03(6)*EnR(3,3)&
         +p%S03(7)*EnR(4,1)+p%S03(8)*EnR(4,2)+p%S03(9)*EnR(4,3)+p%S03(10)*EnR(5,1)+p%S03(11)*EnR(5,2)&
         +p%S03(12)*EnR(5,3)
      S1=p%S1(1)*d1r3+p%S1(2)*d1r1+p%S1(3)*EnR(2,3)+p%S1(4)*EnR(3,1)+p%S1(5)*EnR(3,2)+p%S1(6)*EnR(3,3)&
        +p%S1(7)*EnR(4,1)+p%S1(8)*EnR(4,2)+p%S1(9)*EnR(4,3)+p%S1(10)*EnR(5,1)+p%S1(11)*EnR(5,2)+p%S1(12)*EnR(5,3)&
        +p%S1(13)*EnR(6,1)+p%S1(14)*EnR(6,2)+p%S1(15)*EnR(6,3)
      S2=p%S2(1)*d1r3+p%S2(2)*d1r1+p%S2(3)*EnR(2,1)+p%S2(4)*EnR(2,2)+p%S2(5)*EnR(3,1)+p%S2(6)*EnR(3,2)&
        +p%S2(7)*EnR(3,3)+p%S2(8)*EnR(4,1)+p%S2(9)*EnR(4,2)+p%S2(10)*EnR(4,3)+p%S2(11)*EnR(5,1)+p%S2(12)*EnR(5,2)&
        +p%S2(13)*EnR(5,3)+p%S2(14)*EnR(6,1)+p%S2(15)*EnR(6,2)+p%S2(16)*EnR(6,3)
      S3=p%S3(1)*d1r3+p%S3(2)*d1r1+p%S3(3)*EnR(2,1)+p%S3(4)*EnR(2,2)+p%S3(5)*EnR(2,3)+p%S3(6)*EnR(3,1)&
        +p%S3(7)*EnR(3,2)+p%S3(8)*EnR(3,3)+p%S3(9)*EnR(4,1)+p%S3(10)*EnR(4,2)+p%S3(11)*EnR(4,3)+p%S3(12)*EnR(5,1)&
        +p%S3(13)*EnR(5,2)+p%S3(14)*EnR(5,3)+p%S3(15)*EnR(6,1)+p%S3(16)*EnR(6,2)+p%S3(17)*EnR(6,3)
      S4=p%S4(1)*d1r3+p%S4(2)*d1r1+p%S4(3)+p%S4(4)*EnR(3,3)+p%S4(5)*EnR(4,1)+p%S4(6)*EnR(4,2)+p%S4(7)*EnR(4,3)&
        +p%S4(8)*EnR(5,1)+p%S4(9)*EnR(5,2)+p%S4(10)*EnR(5,3)+p%S4(11)*EnR(6,1)+p%S4(12)*EnR(6,2)+p%S4(13)*EnR(6,3)
      S5=p%S5(1)*d1r3+p%S5(2)*d1r1+p%S5(3)+p%S5(4)*EnR(2,1)+p%S5(5)*EnR(2,2)+p%S5(6)*EnR(3,1)+p%S5(7)*EnR(3,2)&
        +p%S5(8)*EnR(4,1)+p%S5(9)*EnR(4,2)+p%S5(10)*EnR(4,3)+p%S5(11)*EnR(5,1)+p%S5(12)*EnR(5,2)+p%S5(13)*EnR(5,3)&
        +p%S5(14)*EnR(6,1)+p%S5(15)*EnR(6,2)+p%S5(16)*EnR(6,3)
      ! Fundamental solutions
      fs_u(0,0)=eta
      fs_t(0,0)=W0*drdn
      fs_d(0,0)=W0*drdni
      fs_s(0,0)=Q1*drdn*drdni+Q2*n_dot_ni
      fs_u(0,1:3)=vartheta*drdx
      fs_t(0,1:3)=T01*drdx*drdn+T02*n
      fs_d(0,1:3)=-W1*drdx*drdni+W2*n_i
      fs_s(0,1:3)=S01*drdx*drdn*drdni+S02*n*drdni+S03*(n_i*drdn+drdx*n_dot_ni)
      fs_u(1:3,0)=vartheta*drdx
      fs_t(1:3,0)=W1*drdx*drdn+W2*n
      fs_d(1:3,0)=-T01*drdx*drdni+T02*n_i
      fs_s(1:3,0)=-S01*drdx*drdn*drdni+S02*n_i*drdn-S03*(-n*drdni+drdx*n_dot_ni)
      do il=1,3
        do ik=1,3
          fs_u(il,ik)=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
          fs_t(il,ik)=TT1*drdx(il)*drdx(ik)*drdn+TT2*(drdn*c_dkr(il,ik)+drdx(ik)*n(il))+TT3*drdx(il)*n(ik)
          fs_d(il,ik)=TT1*drdx(il)*drdx(ik)*drdni-TT2*(drdx(il)*n_i(ik)-c_dkr(il,ik)*drdni)-TT3*drdx(ik)*n_i(il)
          fs_s(il,ik)=S1*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni-c_dkr(il,ik)*drdn*drdni+drdx(il)*drdx(ik)*n_dot_ni)&
                     +S2*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni)+S3*drdx(il)*drdx(ik)*drdn*drdni&
                     +S4*(c_dkr(il,ik)*n_dot_ni+n_i(ik)*n(il))+S5*n(ik)*n_i(il)
        end do
      end do
      ! Add the term of the numerical integration summation (fundamental solutions * functional shape functions * jacobian * weight)
      do il=0,3
        do ik=0,3
          h(:,il,ik)=h(:,il,ik)+fs_t(il,ik)*pphijw
          g(:,il,ik)=g(:,il,ik)+fs_u(il,ik)*sphijw
          m(:,il,ik)=m(:,il,ik)+fs_s(il,ik)*pphijw
          l(:,il,ik)=l(:,il,ik)+fs_d(il,ik)*sphijw
        end do
      end do
    end do
    ! Multiply by constants
    cte_u(0,0)=-c_1_4pi
    cte_t(0,0)=-c_1_4pi
    cte_d(0,0)=-c_1_4pi/p%J
    cte_s(0,0)=-c_1_4pi
    cte_u(0,1:3)=-c_1_4pi
    cte_t(0,1:3)= c_1_4pi
    cte_d(0,1:3)=c_1_4pi/p%mu
    cte_s(0,1:3)=c_1_4pi
    cte_u(1:3,0)=-c_1_4pi/p%J
    cte_t(1:3,0)=-c_1_4pi/p%mu
    cte_d(1:3,0)=-c_1_4pi/p%J
    cte_s(1:3,0)=-c_1_4pi
    cte_u(1:3,1:3)=c_1_4pi/p%mu
    cte_t(1:3,1:3)=c_1_4pi
    cte_d(1:3,1:3)=c_1_4pi
    cte_s(1:3,1:3)=c_1_4pi*p%mu
    do il=0,3
      do ik=0,3
        h(:,il,ik)=cte_t(il,ik)*h(:,il,ik)
        g(:,il,ik)=cte_u(il,ik)*g(:,il,ik)
        m(:,il,ik)=cte_s(il,ik)*m(:,il,ik)
        l(:,il,ik)=cte_d(il,ik)*l(:,il,ik)
      end do
    end do
    ! Reverse if needed
    if (reverse) then
      h=-h
      m=-m
    end if
  end subroutine fbem_bem_harpor3d_shbie_ext_pre

  subroutine fbem_bem_harpor3d_shbie_auto(e,reverse,x_i,n_i,p,qsp,ns,h,g,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                     !! Integration element
    logical                            :: reverse               !! Reverse orientation
    real(kind=real64)                  :: x_i(3)                !! Collocation point
    real(kind=real64)                  :: n_i(3)                !! Unit normal at the collocation point
    type(fbem_bem_harpor3d_parameters) :: p                     !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp                   !! Quasi-singular integration parameters
    integer                            :: ns                    !! Maximum level of subdivisions
    complex(kind=real64)               :: h(e%n_pnodes,0:3,0:3) !! h integration kernel
    complex(kind=real64)               :: g(e%n_snodes,0:3,0:3) !! g integration kernel
    complex(kind=real64)               :: m(e%n_pnodes,0:3,0:3) !! m integration kernel
    complex(kind=real64)               :: l(e%n_snodes,0:3,0:3) !! l integration kernel
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
        call fbem_bem_harpor3d_sbie_int(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi,p,h,g)
        call fbem_bem_harpor3d_hbie_int(e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi,p,m,l)
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
          call fbem_bem_harpor3d_shbie_ext_pre(ps,e,reverse,x_i,n_i,p,h,g,m,l)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_harpor3d_sbie_ext_adp(e,reverse,xi_s,x_i,p,qsp,1,ns,h,g)
          call fbem_bem_harpor3d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,p,qsp,1,ns,m,l)
        end if
    end select
  end subroutine fbem_bem_harpor3d_shbie_auto

  ! ================================================================================================================================

end module fbem_bem_harpor3d
