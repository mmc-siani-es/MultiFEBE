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
!! <b> This module implements the calculation of element-wise integrals of Boundary Integral Equations of the 2D Biot's
!! poroelastodynamics.</b>
module fbem_bem_harpor2d

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
  ! INITIAL SETUP
  public :: fbem_bem_harpor2d_parameters
  public :: fbem_bem_harpor2d_calculate_basic_parameters
  public :: fbem_bem_harpor2d_calculate_parameters
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE)
  ! Exterior integration
  public :: fbem_bem_harpor2d_sbie_ext_pre
  public :: fbem_bem_harpor2d_sbie_ext_st
  public :: fbem_bem_harpor2d_sbie_ext_adp
  ! Interior integration
  public :: fbem_bem_harpor2d_sbie_int
  ! Automatic integration
  public :: fbem_bem_harpor2d_sbie_auto
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (HBIE)
  ! Exterior integration
  public :: fbem_bem_harpor2d_hbie_ext_pre
  public :: fbem_bem_harpor2d_hbie_ext_st
  public :: fbem_bem_harpor2d_hbie_ext_adp
  ! Interior integration
  public :: fbem_bem_harpor2d_hbie_int
  ! Automatic integration
  public :: fbem_bem_harpor2d_hbie_auto
  ! ================================================================================================================================

  ! ================================================================================================================================
  ! BOTH SBIE AND HBIE SIMULTANEOUSLY
  ! Exterior integration
  public :: fbem_bem_harpor2d_shbie_ext_pre
  ! Automatic integration
  public :: fbem_bem_harpor2d_shbie_auto
  ! ================================================================================================================================

  ! ================================================================================================================================
  !! Data structure for region parameters (material properties and frequency dependant parameters)
  type fbem_bem_harpor2d_parameters
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
    ! Coefficients of the components of the fundamental solutions
    complex(kind=real64) :: eta(4)
    complex(kind=real64) :: vartheta(4)
    complex(kind=real64) :: psi(5), chi(5)
    complex(kind=real64) :: W0(5)
    complex(kind=real64) :: T01(5), T02(6)
    complex(kind=real64) :: W1(10), W2(6)
    complex(kind=real64) :: T1(8), T2(8), T3(7)
    complex(kind=real64) :: Q01(11), Q02(7)
    complex(kind=real64) :: S01(13), S02(12), S03(11)
    complex(kind=real64) :: S1(9), S2(10), S3(11), S4(7), S5(10)
  end type fbem_bem_harpor2d_parameters
  ! ================================================================================================================================


contains

  ! ================================================================================================================================
  ! INITIAL SETUP

  !! Subroutine that calculate basic parameters for a given region properties and a given frequency
  subroutine fbem_bem_harpor2d_calculate_basic_parameters(lambda,mu,rho1,rho2,rhoa,R,Q,b,omega,pars)
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
    type(fbem_bem_harpor2d_parameters) :: pars
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
  end subroutine fbem_bem_harpor2d_calculate_basic_parameters

  !! Subroutine that calculate all parameters for a given region properties and a given frequency
  subroutine fbem_bem_harpor2d_calculate_parameters(lambda,mu,rho1,rho2,rhoa,R,Q,b,omega,pars)
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
    type(fbem_bem_harpor2d_parameters) :: pars
    ! Local
    complex(kind=real64)               :: rhohat11, rhohat22, rhohat12
    complex(kind=real64)               :: Z, J
    complex(kind=real64)               :: ca, cb, cc
    complex(kind=real64)               :: k1, k2, k3, tmpk
    complex(kind=real64)               :: alpha1, alpha2, beta1, beta2, vartheta_c
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
    !
    ! Coefficients of equations of fundamental solution's components
    !
    ! Constants
    alpha1=k1**2-mu/(lambda+2.0*mu)*k3**2
    alpha2=k2**2-mu/(lambda+2.0*mu)*k3**2
    beta1=mu/(lambda+2.0*mu)*k1**2-k1**2*k2**2/k3**2
    beta2=mu/(lambda+2.0*mu)*k2**2-k1**2*k2**2/k3**2
    vartheta_c=(Q/R-Z)/(lambda+2.0d0*mu)/(k1**2-k2**2)
    ! eta
    pars%eta(1)=-1.0d0
    pars%eta(2)=-((alpha1*zlog(k1)-alpha2*zlog(k2))/(k1**2-k2**2)+zlog(0.5d0*c_im)+c_gamma)
    pars%eta(3)= alpha1/(k1**2-k2**2)
    pars%eta(4)=-alpha2/(k1**2-k2**2)
    ! vartheta
    pars%vartheta(1)=-0.5d0*vartheta_c*(k1**2-k2**2)
    pars%vartheta(2)=-0.5d0*vartheta_c*(k1**2*zlog(k1)-k2**2*zlog(k2)+(k1**2-k2**2)*(zlog(0.5d0*c_im)+c_gamma-0.5d0))
    pars%vartheta(3)= vartheta_c*c_im*k1
    pars%vartheta(4)=-vartheta_c*c_im*k2
    ! psi
    pars%psi(1)=-0.5d0*(lambda+3.0d0*mu)/(lambda+2.0d0*mu)
    pars%psi(2)=-0.5d0*(zlog(k3)+(beta1*zlog(k1)-beta2*zlog(k2))/(k1**2-k2**2)&
                +(lambda+3.0d0*mu)/(lambda+2.0d0*mu)*(zlog(0.5d0*c_im)+c_gamma)+0.5d0*(lambda+mu)/(lambda+2.0d0*mu))
    pars%psi(3)= 1.0d0/(c_im*k3)
    pars%psi(4)=-beta1/(k1**2-k2**2)/(c_im*k1)
    pars%psi(5)= beta2/(k1**2-k2**2)/(c_im*k2)
    ! chi
    pars%chi(1)=-0.5d0*(lambda+mu)/(lambda+2.0d0*mu)
    pars%chi(2)=0.125d0*(k3**2-(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))
    pars%chi(3)=0.125d0*(k3**2*zlog(k3)-(k1**2*beta1*zlog(k1)-k2**2*beta2*zlog(k2))/(k1**2-k2**2)&
               +(k3**2-(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))*(zlog(0.5d0*c_im)+c_gamma-0.75d0))
    pars%chi(4)=-beta1/(k1**2-k2**2)
    pars%chi(5)= beta2/(k1**2-k2**2)
    ! W_0
    pars%W0(1)= J
    pars%W0(2)=-0.5d0*(Z*vartheta_c*(k1**2-k2**2)+J*(k1**2*alpha1-k2**2*alpha2)/(k1**2-k2**2))
    pars%W0(3)=-0.5d0*(Z*vartheta_c*(k1**2*zlog(k1)-k2**2*zlog(k2)+(k1**2-k2**2)*(zlog(0.5d0*c_im)+c_gamma-0.5d0))&
                +J/(k1**2-k2**2)*(k1**2*alpha1*zlog(k1)-k2**2*alpha2*zlog(k2)&
                +(k1**2*alpha1-k2**2*alpha2)*(zlog(0.5d0*c_im)+c_gamma-0.5d0)))
    pars%W0(4)= (Z*vartheta_c+J*alpha1/(k1**2-k2**2))*c_im*k1
    pars%W0(5)=-(Z*vartheta_c+J*alpha2/(k1**2-k2**2))*c_im*k2
    ! T_01
    pars%T01(1)=mu*vartheta_c*(k1**2-k2**2)
    pars%T01(2)=-2.0d0*mu*vartheta_c*k1**2
    pars%T01(3)= 2.0d0*mu*vartheta_c*k2**2
    pars%T01(4)= 4.0d0*mu*vartheta_c*c_im*k1
    pars%T01(5)=-4.0d0*mu*vartheta_c*c_im*k2
    ! T_02
    pars%T02(1)=(Q/R-Z)*(lambda+mu)/(lambda+2.0d0*mu)-Q/R
    pars%T02(2)=-0.5d0*(Q/R-Z)*mu/(lambda+2.0d0*mu)&
                +(Q/R-Z)*(lambda+mu)/(lambda+2.0d0*mu)*((k1**2*zlog(k1)-k2**2*zlog(k2))/(k1**2-k2**2)+zlog(0.5d0*c_im)+c_gamma)&
                -Q/R*((alpha1*zlog(k1)-alpha2*zlog(k2))/(k1**2-k2**2)+zlog(0.5d0*c_im)+c_gamma)
    pars%T02(3)=-lambda*vartheta_c*k1**2+Q/R*alpha1/(k1**2-k2**2)
    pars%T02(4)= lambda*vartheta_c*k2**2-Q/R*alpha2/(k1**2-k2**2)
    pars%T02(5)=-2.0d0*mu*vartheta_c*c_im*k1
    pars%T02(6)= 2.0d0*mu*vartheta_c*c_im*k2
    ! W_1
    pars%W1(1)=0.5d0*(Q/R*mu/(lambda+2.0d0*mu)-Z)
    pars%W1(2)=0.125d0*Z*(k3**2-(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))
    pars%W1(3)=0.125d0*Z*(k3**2*zlog(k3)-(k1**2*beta1*zlog(k1)-k2**2*beta2*zlog(k2))/(k1**2-k2**2)&
               +(k3**2-(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))*(zlog(0.5d0*c_im)+c_gamma-0.75d0))
    pars%W1(4)=-mu*vartheta_c*k1**2
    pars%W1(5)= mu*vartheta_c*k2**2
    pars%W1(6)= 2.0d0*mu*vartheta_c*c_im*k1
    pars%W1(7)=-2.0d0*mu*vartheta_c*c_im*k2
    pars%W1(8)= Z
    pars%W1(9)=-Z*beta1/(k1**2-k2**2)
    pars%W1(10)=Z*beta2/(k1**2-k2**2)
    ! W_2
    pars%W2(1)=0.5d0*(Q/R*mu/(lambda+2.0d0*mu)+Z)
    pars%W2(2)=0.5d0*(Z*(zlog(k3)-k1**2*k2**2/k3**2*(zlog(k1)-zlog(k2))/(k1**2-k2**2)+zlog(0.5d0*c_im)+c_gamma+0.5d0)&
               +Q/R*mu/(lambda+2.0d0*mu)*((k1**2*zlog(k1)-k2**2*zlog(k2))/(k1**2-k2**2)+zlog(0.5d0*c_im)+c_gamma-0.5d0))
    pars%W2(3)=-Z
    pars%W2(4)=-Z/(c_im*k3)
    pars%W2(5)= Z*beta1/(k1**2-k2**2)/(c_im*k1)-mu*vartheta_c*c_im*k1
    pars%W2(6)=-Z*beta2/(k1**2-k2**2)/(c_im*k2)+mu*vartheta_c*c_im*k2
    ! T_1
    pars%T1(1)=mu/(lambda+2.0d0*mu)
    pars%T1(2)=-0.5d0*(0.5d0*(k3**2-(2.0d0*lambda+mu)/mu*(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))&
                     +Q/R/J*(Q/R-Z)/(lambda+2.0d0*mu))
    pars%T1(3)=-0.5d0*(0.5d0*(k3**2*zlog(k3)-(2.0d0*lambda+mu)/mu*(k1**2*beta1*zlog(k1)-k2**2*beta2*zlog(k2))/(k1**2-k2**2)&
               +(k3**2-(2.0d0*lambda+mu)/mu*(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))*(zlog(0.5d0*c_im)+c_gamma))&
               -0.125d0*(3.0d0*k3**2-(4.0d0*lambda+3.0d0*mu)/mu*(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))&
               +Q/R/J*(Q/R-Z)/(lambda+2.0d0*mu)*((k1**2*zlog(k1)-k2**2*zlog(k2))/(k1**2-k2**2)&
               +zlog(0.5d0*c_im)+c_gamma-0.5d0))
    pars%T1(4)= (Q/R/J*vartheta_c-lambda/mu*beta1/(k1**2-k2**2))*c_im*k1
    pars%T1(5)=-(Q/R/J*vartheta_c-lambda/mu*beta2/(k1**2-k2**2))*c_im*k2
    pars%T1(6)=-2.0d0
    pars%T1(7)= 2.0d0*beta1/(k1**2-k2**2)
    pars%T1(8)=-2.0d0*beta2/(k1**2-k2**2)
    ! T_2
    pars%T2(1)=-2.0d0*(lambda+mu)/(lambda+2.0d0*mu)
    pars%T2(2)=-0.25d0*(k3**2-(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))
    pars%T2(3)= 2.0d0*c_im*k3
    pars%T2(4)=-2.0d0*c_im*k1*beta1/(k1**2-k2**2)
    pars%T2(5)= 2.0d0*c_im*k2*beta2/(k1**2-k2**2)
    pars%T2(6)= 8.0d0
    pars%T2(7)=-8.0d0*beta1/(k1**2-k2**2)
    pars%T2(8)= 8.0d0*beta2/(k1**2-k2**2)
    ! T_3
    pars%T3(1)=-mu/(lambda+2.0d0*mu)
    pars%T3(2)=0.25d0*(k3**2+(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))
    pars%T3(3)=0.25d0*(k3**2*zlog(k3)+(k1**2*beta1*zlog(k1)-k2**2*beta2*zlog(k2))/(k1**2-k2**2)&
               +(k3**2+(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))*(zlog(0.5d0*c_im)+c_gamma)&
               -0.25d0*(k3**2+3.0d0*(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2)))
    pars%T3(4)=-c_im*k3
    pars%T3(5)=-2.0d0
    pars%T3(6)= 2.0d0*beta1/(k1**2-k2**2)
    pars%T3(7)=-2.0d0*beta2/(k1**2-k2**2)
    ! Q_01
    pars%Q01(1) =2.0d0*J
    pars%Q01(2) =0.5d0/mu*(-Z**2*(lambda+3.0d0*mu)/(lambda+2.0d0*mu)&
                 +Z*Q/R*2.0d0*mu/(lambda+2.0d0*mu)+mu*J*(k1**2*alpha1-k2**2*alpha2)/(k1**2-k2**2))
    pars%Q01(3) =0.125d0*Z**2/mu*(k3**2-(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))
    pars%Q01(4) =0.125d0*Z**2/mu*(k3**2*zlog(k3)-(k1**2*beta1*zlog(k1)-k2**2*beta2*zlog(k2))/(k1**2-k2**2)&
                 +(k3**2-(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))*(zlog(0.5d0*c_im)+c_gamma-0.75d0))
    pars%Q01(5) =-(2.0d0*Z*vartheta_c+J*alpha1/(k1**2-k2**2))*k1**2
    pars%Q01(6) = (2.0d0*Z*vartheta_c+J*alpha2/(k1**2-k2**2))*k2**2
    pars%Q01(7) = 2.0d0*(2.0d0*Z*vartheta_c+J*alpha1/(k1**2-k2**2))*c_im*k1
    pars%Q01(8) =-2.0d0*(2.0d0*Z*vartheta_c+J*alpha2/(k1**2-k2**2))*c_im*k2
    pars%Q01(9) = Z**2/mu
    pars%Q01(10)=-Z**2/mu*beta1/(k1**2-k2**2)
    pars%Q01(11)= Z**2/mu*beta2/(k1**2-k2**2)
    ! Q_02

    pars%Q02(1)=J
    pars%Q02(2)=-0.5d0*(Z**2/mu*(lambda+3.0d0*mu)/(lambda+2.0d0*mu)+(Q/R-Z)*2.0d0*Z/(lambda+2.0d0*mu)&
               +J*(k1**2*alpha1-k2**2*alpha2)/(k1**2-k2**2))
    pars%Q02(3)=-0.5d0*(Z**2/mu*(zlog(k3)+(beta1*zlog(k1)-beta2*zlog(k2))/(k1**2-k2**2)&
                +(lambda+3.0d0*mu)/(lambda+2.0d0*mu)*(zlog(0.5d0*c_im)+c_gamma)+0.5d0*(lambda+mu)/(lambda+2.0d0*mu))&
                +(Q/R-Z)*2.0d0*Z/(lambda+2.0d0*mu)*((k1**2*zlog(k1)-k2**2*zlog(k2))/(k1**2-k2**2)+zlog(0.5d0*c_im)+c_gamma-0.5d0)&
                +J*(k1**2*alpha1*zlog(k1)-k2**2*alpha2*zlog(k2)+(k1**2*alpha1-k2**2*alpha2)*(zlog(0.5d0*c_im)+c_gamma-0.5d0))&
                /(k1**2-k2**2))
    pars%Q02(4)= Z**2/mu
    pars%Q02(5)= Z**2/mu/(c_im*k3)
    pars%Q02(6)=-Z**2/mu*beta1/(k1**2-k2**2)/(c_im*k1)+2.0d0*Z*vartheta_c*c_im*k1+J*c_im*k1*alpha1/(k1**2-k2**2)
    pars%Q02(7)= Z**2/mu*beta2/(k1**2-k2**2)/(c_im*k2)-2.0d0*Z*vartheta_c*c_im*k2-J*c_im*k2*alpha2/(k1**2-k2**2)
    ! S_01
    pars%S01(1) =2.0d0*(Q/R*mu/(lambda+2.0d0*mu)-Z)
    pars%S01(2) =vartheta_c*mu*(k1**4-k2**4)
    pars%S01(3) =-0.25d0*Z*(k3**2-(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))&
                 +vartheta_c*mu*(k1**4*zlog(k1)-k2**4*zlog(k2)+(k1**4-k2**4)*(zlog(0.5d0*c_im)+c_gamma-0.5d0))
    pars%S01(4) = -8.0d0*mu*vartheta_c*k1**2
    pars%S01(5) =  8.0d0*mu*vartheta_c*k2**2
    pars%S01(6) = 16.0d0*mu*vartheta_c*c_im*k1
    pars%S01(7) =-16.0d0*mu*vartheta_c*c_im*k2
    pars%S01(8) =  2.0d0*Z*c_im*k3
    pars%S01(9) = (2.0d0*Z*k2**2/k3**2-Q/R*2.0d0*mu/(lambda+2.0d0*mu))*c_im*k1**3/(k1**2-k2**2)
    pars%S01(10)=-(2.0d0*Z*k1**2/k3**2-Q/R*2.0d0*mu/(lambda+2.0d0*mu))*c_im*k2**3/(k1**2-k2**2)
    pars%S01(11)= 8.0d0*Z
    pars%S01(12)=-8.0d0*Z*beta1/(k1**2-k2**2)
    pars%S01(13)= 8.0d0*Z*beta2/(k1**2-k2**2)
    ! S_02
    pars%S02(1) =Q/R*mu/(lambda+2.0d0*mu)+Z
    pars%S02(2) =0.5d0*((Q/R*mu/(lambda+2.0d0*mu)-0.5d0*Z)*k3**2&
                 -(Q/R-0.25d0*Z)*2.0d0*mu/(lambda+2.0d0*mu)*(k1**4-k2**4)/(k1**2-k2**2)&
                 -Z/J/R/(lambda+2.0d0*mu)*(lambda+0.5d0*mu+Q*(Q/R-Z)))
    pars%S02(3) =0.5d0*(-0.5d0*Z*(k3**2*zlog(k3)-(2.0d0*lambda+mu)/mu*(k1**2*beta1*zlog(k1)-k2**2*beta2*zlog(k2))/(k1**2-k2**2)&
                 +(k3**2-(2.0d0*lambda+mu)/mu*(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))*(zlog(0.5d0*c_im)+c_gamma))&
                 +0.125d0*Z*(3.0d0*k3**2-(4.0d0*lambda+3.0d0*mu)/mu*(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))&
                 -Q*Z/R/J*(Q/R-Z)/(lambda+2.0d0*mu)*((k1**2*zlog(k1)-k2**2*zlog(k2))/(k1**2-k2**2)+zlog(0.5d0*c_im)+c_gamma-0.5d0)&
                 +(Q/R-Z)*lambda/(lambda+2.0d0*mu)*(k1**4*zlog(k1)-k2**4*zlog(k2)&
                 +(k1**4-k2**4)*(zlog(0.5d0*c_im)+c_gamma-0.5d0))/(k1**2-k2**2)&
                 -Q/R*(k1**2*alpha1*zlog(k1)-k2**2*alpha2*zlog(k2)+(k1**2*alpha1-k2**2*alpha2)*(zlog(0.5d0*c_im)&
                 +c_gamma-0.5d0))/(k1**2-k2**2))
    pars%S02(4) = 2.0d0*mu*vartheta_c*k1**2
    pars%S02(5) =-2.0d0*mu*vartheta_c*k2**2
    pars%S02(6) =-4.0d0*mu*vartheta_c*c_im*k1
    pars%S02(7) = 4.0d0*mu*vartheta_c*c_im*k2
    pars%S02(8) = (vartheta_c*(Q/R*Z/J-lambda*k1**2)+Q/R*alpha1/(k1**2-k2**2)-Z*lambda/mu*beta1/(k1**2-k2**2))*c_im*k1
    pars%S02(9) =-(vartheta_c*(Q/R*Z/J-lambda*k2**2)+Q/R*alpha2/(k1**2-k2**2)-Z*lambda/mu*beta2/(k1**2-k2**2))*c_im*k2
    pars%S02(10)=-2.0d0*Z
    pars%S02(11)= 2.0d0*Z*beta1/(k1**2-k2**2)
    pars%S02(12)=-2.0d0*Z*beta2/(k1**2-k2**2)
    ! S_03
    pars%S03(1) =Q/R*mu/(lambda+2.0d0*mu)
    pars%S03(2) =-0.25d0*Z*(k3**2+(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))
    pars%S03(3) = 0.25d0*Z*(0.25d0*(k3**2+3.0d0*(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))&
                  -(k3**2*zlog(k3)+(k1**2*beta1*zlog(k1)-k2**2*beta2*zlog(k2))/(k1**2-k2**2)&
                  +(k3**2+(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))*(zlog(0.5d0*c_im)+c_gamma)))
    pars%S03(4) =-2.0d0*mu*vartheta_c*k1**2
    pars%S03(5) = 2.0d0*mu*vartheta_c*k2**2
    pars%S03(6) = 4.0d0*mu*vartheta_c*c_im*k1
    pars%S03(7) =-4.0d0*mu*vartheta_c*c_im*k2
    pars%S03(8) = Z*c_im*k3
    pars%S03(9) = 2.0d0*Z
    pars%S03(10)=-2.0d0*Z*beta1/(k1**2-k2**2)
    pars%S03(11)= 2.0d0*Z*beta2/(k1**2-k2**2)
    ! S_1
    pars%S1(1)= 2.0d0*lambda*mu/(lambda+2.0d0*mu)
    pars%S1(2)=-0.5d0*mu*(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2)
    pars%S1(3)= mu*k3**2
    pars%S1(4)=-6.0d0*mu*c_im*k3
    pars%S1(5)= 4.0d0*mu*c_im*k1*beta1/(k1**2-k2**2)
    pars%S1(6)=-4.0d0*mu*c_im*k2*beta2/(k1**2-k2**2)
    pars%S1(7)=-16.0d0*mu
    pars%S1(8)= 16.0d0*mu*beta1/(k1**2-k2**2)
    pars%S1(9)=-16.0d0*mu*beta2/(k1**2-k2**2)
    ! S_2
    pars%S2(1) =4.0d0*mu**2/(lambda+2.0d0*mu)
    pars%S2(2) =0.5d0*(mu*k3**2-(2.0d0*lambda+mu)*(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))+mu*Q/R/J*vartheta_c*(k1**2-k2**2)
    pars%S2(3) = 2.0d0*(lambda*beta1/(k1**2-k2**2)-mu*Q/R/J*vartheta_c)*k1**2
    pars%S2(4) =-2.0d0*(lambda*beta2/(k1**2-k2**2)-mu*Q/R/J*vartheta_c)*k2**2
    pars%S2(5) =-4.0d0*mu*c_im*k3
    pars%S2(6) = 4.0d0*((mu-lambda)*beta1/(k1**2-k2**2)+mu*Q/R/J*vartheta_c)*c_im*k1
    pars%S2(7) =-4.0d0*((mu-lambda)*beta2/(k1**2-k2**2)+mu*Q/R/J*vartheta_c)*c_im*k2
    pars%S2(8) =-16.0d0*mu
    pars%S2(9) = 16.0d0*mu*beta1/(k1**2-k2**2)
    pars%S2(10)=-16.0d0*mu*beta2/(k1**2-k2**2)
    ! S_3
    pars%S3(1)=16.0d0*mu*(lambda+mu)/(lambda+2.0d0*mu)
    pars%S3(2) =mu*(k3**2-(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))
    pars%S3(3) =4.0d0*mu*k3**2
    pars%S3(4) =-4.0d0*mu*k1**2*beta1/(k1**2-k2**2)
    pars%S3(5) = 4.0d0*mu*k2**2*beta2/(k1**2-k2**2)
    pars%S3(6) =-32.0d0*mu*c_im*k3
    pars%S3(7) = 32.0d0*mu*c_im*k1*beta1/(k1**2-k2**2)
    pars%S3(8) =-32.0d0*mu*c_im*k2*beta2/(k1**2-k2**2)
    pars%S3(9) =-96.0d0*mu
    pars%S3(10)= 96.0d0*mu*beta1/(k1**2-k2**2)
    pars%S3(11)=-96.0d0*mu*beta2/(k1**2-k2**2)
    ! S_4
    pars%S4(1)=2.0d0*mu**2/(lambda+2.0d0*mu)
    pars%S4(2)=-0.5d0*mu*(k3**2+(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))
    pars%S4(3)=-0.5d0*mu*(k3**2*zlog(k3)+(k1**2*beta1*zlog(k1)-k2**2*beta2*zlog(k2))/(k1**2-k2**2)&
               +(k3**2+(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))*(zlog(0.5d0*c_im)+c_gamma)&
               -0.25d0*(k3**2+3.0d0*(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2)))
    pars%S4(4)= 2.0d0*mu*c_im*k3
    pars%S4(5)= 4.0d0*mu
    pars%S4(6)=-4.0d0*mu*beta1/(k1**2-k2**2)
    pars%S4(7)= 4.0d0*mu*beta2/(k1**2-k2**2)
    ! S_5
    pars%S5(1) =2.0d0*mu*(lambda-mu)/(lambda+2.0d0*mu)
    pars%S5(2) =0.5d0*mu*(k3**2-(2.0d0*lambda/mu*(lambda/mu+2.0d0)+1.0d0)*(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))&
                +Q/R/J*(2.0d0*(Q/R-Z)*(lambda+mu)/(lambda+2.0d0*mu)-Q/R)
    pars%S5(3) =0.5d0*mu*(k3**2*zlog(k3)&
                -(2.0d0*lambda/mu*(lambda/mu+2.0d0)+1.0d0)*(k1**2*beta1*zlog(k1)-k2**2*beta2*zlog(k2))/(k1**2-k2**2)&
                +(k3**2-(2.0d0*lambda/mu*(lambda/mu+2.0d0)+1.0d0)*(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))*(zlog(0.5d0*c_im)&
                +c_gamma))&
                -0.125d0*mu*(3.0d0*k3**2-(8.0d0*lambda+3.0d0*mu)/mu*(k1**2*beta1-k2**2*beta2)/(k1**2-k2**2))&
                +Q/R/J*((Q/R-Z)*2.0d0*(lambda+mu)/(lambda+2.0d0*mu)*((k1**2*zlog(k1)-k2**2*zlog(k2))/(k1**2-k2**2)&
                +zlog(0.5d0*c_im)+c_gamma-0.5d0*mu/(lambda+mu))&
                -Q/R*((alpha1*zlog(k1)-alpha2*zlog(k2))/(k1**2-k2**2)+zlog(0.5d0*c_im)+c_gamma))
    pars%S5(4) = lambda**2/mu*k1**2*beta1/(k1**2-k2**2)+Q/R/J*(Q/R*alpha1/(k1**2-k2**2)-2.0d0*lambda*vartheta_c*k1**2)
    pars%S5(5) =-lambda**2/mu*k2**2*beta2/(k1**2-k2**2)-Q/R/J*(Q/R*alpha2/(k1**2-k2**2)-2.0d0*lambda*vartheta_c*k2**2)
    pars%S5(6) =-4.0d0*(mu*Q/R/J*vartheta_c-lambda*beta1/(k1**2-k2**2))*c_im*k1
    pars%S5(7) = 4.0d0*(mu*Q/R/J*vartheta_c-lambda*beta2/(k1**2-k2**2))*c_im*k2
    pars%S5(8) = 4.0d0*mu
    pars%S5(9) =-4.0d0*mu*beta1/(k1**2-k2**2)
    pars%S5(10)= 4.0d0*mu*beta2/(k1**2-k2**2)
  end subroutine fbem_bem_harpor2d_calculate_parameters

  ! ================================================================================================================================

  ! ================================================================================================================================
  ! SINGULAR BOUNDARY INTEGRAL EQUATION (SBIE)

  subroutine fbem_bem_harpor2d_sbie_ext_pre(ps,e,reverse,x_i,p,h,g)
    implicit none
    ! I/O
    integer                            :: ps                    !! Selected precalculated dataset
    type(fbem_bem_element)             :: e                     !! Element
    logical                            :: reverse               !! Reverse normal vector
    real(kind=real64)                  :: x_i(2)                !! Collocation point position vector
    type(fbem_bem_harpor2d_parameters) :: p                     !! Parameters of the region
    complex(kind=real64)               :: h(e%n_pnodes,0:2,0:2) !! h integration kernels matrix
    complex(kind=real64)               :: g(e%n_snodes,0:2,0:2) !! g integration kernels matrix
    ! Local
    integer              :: kip                            ! Counter variable for integration points loop
    real(kind=real64)    :: x(2)                           ! Position vector at integration point
    real(kind=real64)    :: n(2)                           ! Unit normal vector at integration point
    real(kind=real64)    :: pphijw(e%n_pnodes)             ! phi^p * jacobian * weight at integration point
    real(kind=real64)    :: sphijw(e%n_snodes)             ! phi^s * jacobian * weight at integration point
    real(kind=real64)    :: rv(2)                          ! Radiovector between collocation point and integration point (x-x_i)
    real(kind=real64)    :: r, r2, d1r1, logr              ! Radiovector module, squared, inverses and log(r)
    real(kind=real64)    :: drdx(2)                        ! Radiovector derivatives with respect to x_k
    real(kind=real64)    :: drdn                           ! Partial derivative of r respect to unit normal
    integer              :: il, ik                         ! Counter variables for load direction and observation direction
    complex(kind=real64) :: z(3)                           ! Bessel functions arguments: z=ikr
    complex(kind=real64) :: KnR(0:2,3)                     ! Bessel functions decomposition
    complex(kind=real64) :: eta, vartheta, psi, chi        ! Fundamental solution components
    complex(kind=real64) :: W0, T01, T02, W1, W2           ! Fundamental solution components
    complex(kind=real64) :: T1, T2, T3                     ! Fundamental solution components
    complex(kind=real64) :: fs_u(0:2,0:2), fs_t(0:2,0:2)   ! Fundamental solutions
    complex(kind=real64) :: cte_u(0:2,0:2), cte_t(0:2,0:2) ! Fundamental solution constants
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
      z(3)=c_im*p%k3*r
      call fbem_BesselKnR_decomposed(3,z,KnR)
      eta=p%eta(1)*logr+p%eta(2)+p%eta(3)*KnR(0,1)+p%eta(4)*KnR(0,2)
      vartheta=(p%vartheta(1)*logr+p%vartheta(2))*r+p%vartheta(3)*KnR(1,1)+p%vartheta(4)*KnR(1,2)
      psi=p%psi(1)*logr+p%psi(2)+KnR(0,3)+(p%psi(3)*KnR(1,3)+p%psi(4)*KnR(1,1)+p%psi(5)*KnR(1,2))*d1r1
      chi=p%chi(1)+(p%chi(2)*logr+p%chi(3))*r2+KnR(2,3)+p%chi(4)*KnR(2,1)+p%chi(5)*KnR(2,2)
      W0=p%W0(1)*d1r1+(p%W0(2)*logr+p%W0(3))*r+p%W0(4)*KnR(1,1)+p%W0(5)*KnR(1,2)
      T01=p%T01(1)+p%T01(2)*KnR(0,1)+p%T01(3)*KnR(0,2)+(p%T01(4)*KnR(1,1)+p%T01(5)*KnR(1,2))*d1r1
      T02=p%T02(1)*logr+p%T02(2)+p%T02(3)*KnR(0,1)+p%T02(4)*KnR(0,2)+(p%T02(5)*KnR(1,1)+p%T02(6)*KnR(1,2))*d1r1
      W1=p%W1(1)+(p%W1(2)*logr+p%W1(3))*r2+p%W1(4)*KnR(0,1)+p%W1(5)*KnR(0,2)+(p%W1(6)*KnR(1,1)+p%W1(7)*KnR(1,2))*d1r1&
         +p%W1(8)*KnR(2,3)+p%W1(9)*KnR(2,1)+p%W1(10)*KnR(2,2)
      W2=p%W2(1)*logr+p%W2(2)+p%W2(3)*KnR(0,3)+(p%W2(4)*KnR(1,3)+p%W2(5)*KnR(1,1)+p%W2(6)*KnR(1,2))*d1r1
      T1=p%T1(1)*d1r1+(p%T1(2)*logr+p%T1(3))*r+p%T1(4)*KnR(1,1)+p%T1(5)*KnR(1,2)+(p%T1(6)*KnR(2,3)+p%T1(7)*KnR(2,1)&
        +p%T1(8)*KnR(2,2))*d1r1
      T2=p%T2(1)*d1r1+p%T2(2)*r+p%T2(3)*KnR(1,3)+p%T2(4)*KnR(1,1)+p%T2(5)*KnR(1,2)+(p%T2(6)*KnR(2,3)+p%T2(7)*KnR(2,1)&
        +p%T2(8)*KnR(2,2))*d1r1
      T3=p%T3(1)*d1r1+(p%T3(2)*logr+p%T3(3))*r+p%T3(4)*KnR(1,3)+(p%T3(5)*KnR(2,3)+p%T3(6)*KnR(2,1)+p%T3(7)*KnR(2,2))*d1r1
      ! Fundamental solutions
      fs_u(0,0)=eta
      fs_t(0,0)=W0*drdn
      fs_u(0,1:2)=vartheta*drdx
      fs_t(0,1:2)=T01*drdx*drdn+T02*n
      fs_u(1:2,0)=vartheta*drdx
      fs_t(1:2,0)=W1*drdx*drdn+W2*n
      do il=1,2
        do ik=1,2
          fs_u(il,ik)=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
          fs_t(il,ik)=T1*drdx(il)*n(ik)+T2*drdx(il)*drdx(ik)*drdn+T3*(drdx(ik)*n(il)+drdn*c_dkr(il,ik))
        end do
      end do
      ! Add the term of the numerical integration summation (fundamental solutions * functional shape functions * jacobian * weight)
      do il=0,2
        do ik=0,2
          h(:,il,ik)=h(:,il,ik)+fs_t(il,ik)*pphijw
          g(:,il,ik)=g(:,il,ik)+fs_u(il,ik)*sphijw
        end do
      end do
    end do
    ! Multiply by constants
    cte_u(0,0)=-c_1_2pi
    cte_t(0,0)=-c_1_2pi
    cte_u(0,1:2)=-c_1_2pi
    cte_t(0,1:2)= c_1_2pi
    cte_u(1:2,0)=-c_1_2pi/p%J
    cte_t(1:2,0)=-c_1_2pi/p%mu
    cte_u(1:2,1:2)=c_1_2pi/p%mu
    cte_t(1:2,1:2)=c_1_2pi
    do il=0,2
      do ik=0,2
        h(:,il,ik)=cte_t(il,ik)*h(:,il,ik)
        g(:,il,ik)=cte_u(il,ik)*g(:,il,ik)
      end do
    end do
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_harpor2d_sbie_ext_pre

  subroutine fbem_bem_harpor2d_sbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,p,gln,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                     !! Integration element
    logical                            :: reverse               !! Reverse normal vector
    real(kind=real64)                  :: xi_s(1,2)             !! Coordinates of the subdivision of the parent element (xi space [-1,1])
    real(kind=real64)                  :: x_i(2)                !! Collocation point position vector
    real(kind=real64)                  :: barxip(1)             !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)                  :: barr                  !! Telles jacobian at barxip
    type(fbem_bem_harpor2d_parameters) :: p                     !! Parameters of the region
    integer                            :: gln                   !! 1D Gauss-Legendre number of integration points (<=32)
    complex(kind=real64)               :: h(e%n_pnodes,0:2,0:2) !! h kernel vector
    complex(kind=real64)               :: g(e%n_snodes,0:2,0:2) !! g kernel vector
    ! Local
    integer                      :: kphi                           ! Counter variable for shape functions loops
    integer                      :: kip                            ! Counter variable of integration points
    real(kind=real64)            :: gamma                          ! Coordinate gamma (Telles transformation space [0,1])
    real(kind=real64)            :: w                              ! Weights of an integration point
    type(fbem_telles_parameters) :: telles_parameters              ! Telles parameters
    real(kind=real64)            :: jt                             ! Telles jacobian: xip->gamma
    real(kind=real64)            :: xip                            ! Coordinate xip (subdivision space [0,1])
    real(kind=real64)            :: js                             ! Subdivision jacobian: xi->xip
    real(kind=real64)            :: xi                             ! Coordinate xi [xi1,xi2]
    real(kind=real64)            :: aux(10)                        ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)               ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi(e%n_gnodes)           ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: pphi(e%n_pnodes)               ! Functional shape functions values
    real(kind=real64)            :: sphi(e%n_snodes)               ! Functional shape functions values
    real(kind=real64)            :: x(2)                           ! Position vector at xi
    real(kind=real64)            :: T(2)                           ! Tangent vector at xi
    real(kind=real64)            :: N(2)                           ! Normal vector at xi
    real(kind=real64)            :: rv(2)                          ! Radiovector between collocation point and integration point (x-x_i)
    real(kind=real64)            :: r, r2, d1r1, logr              ! Radiovector module, squared, its inverses and log(r)
    real(kind=real64)            :: drdx(2)                        ! Radiovector derivatives with respect to x_k
    real(kind=real64)            :: jg                             ! Geometric jacobian
    real(kind=real64)            :: drdn                           ! Partial derivative of r respect to unit normal
    real(kind=real64)            :: jw                             ! Jacobian * weight
    real(kind=real64)            :: pphijw(e%n_pnodes)             ! Auxiliary variables for integrand evaluation
    real(kind=real64)            :: sphijw(e%n_snodes)             ! Auxiliary variables for integrand evaluation
    integer                      :: il, ik                         ! Counter variables for load direction and observation direction
    complex(kind=real64)         :: z(3)                           ! Bessel functions arguments: z=ikr
    complex(kind=real64)         :: KnR(0:2,3)                     ! Bessel functions decomposition
    complex(kind=real64)         :: eta, vartheta, psi, chi        ! Fundamental solution components
    complex(kind=real64)         :: W0, T01, T02, W1, W2           ! Fundamental solution components
    complex(kind=real64)         :: T1, T2, T3                     ! Fundamental solution components
    complex(kind=real64)         :: fs_u(0:2,0:2), fs_t(0:2,0:2)   ! Fundamental solutions
    complex(kind=real64)         :: cte_u(0:2,0:2), cte_t(0:2,0:2) ! Fundamental solution constants
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
      d1r1=1.0d0/r
      logr=log(r)
      drdx=rv*d1r1
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
      z(3)=c_im*p%k3*r
      call fbem_BesselKnR_decomposed(3,z,KnR)
      eta=p%eta(1)*logr+p%eta(2)+p%eta(3)*KnR(0,1)+p%eta(4)*KnR(0,2)
      vartheta=(p%vartheta(1)*logr+p%vartheta(2))*r+p%vartheta(3)*KnR(1,1)+p%vartheta(4)*KnR(1,2)
      psi=p%psi(1)*logr+p%psi(2)+KnR(0,3)+(p%psi(3)*KnR(1,3)+p%psi(4)*KnR(1,1)+p%psi(5)*KnR(1,2))*d1r1
      chi=p%chi(1)+(p%chi(2)*logr+p%chi(3))*r2+KnR(2,3)+p%chi(4)*KnR(2,1)+p%chi(5)*KnR(2,2)
      W0=p%W0(1)*d1r1+(p%W0(2)*logr+p%W0(3))*r+p%W0(4)*KnR(1,1)+p%W0(5)*KnR(1,2)
      T01=p%T01(1)+p%T01(2)*KnR(0,1)+p%T01(3)*KnR(0,2)+(p%T01(4)*KnR(1,1)+p%T01(5)*KnR(1,2))*d1r1
      T02=p%T02(1)*logr+p%T02(2)+p%T02(3)*KnR(0,1)+p%T02(4)*KnR(0,2)+(p%T02(5)*KnR(1,1)+p%T02(6)*KnR(1,2))*d1r1
      W1=p%W1(1)+(p%W1(2)*logr+p%W1(3))*r2+p%W1(4)*KnR(0,1)+p%W1(5)*KnR(0,2)+(p%W1(6)*KnR(1,1)+p%W1(7)*KnR(1,2))*d1r1&
         +p%W1(8)*KnR(2,3)+p%W1(9)*KnR(2,1)+p%W1(10)*KnR(2,2)
      W2=p%W2(1)*logr+p%W2(2)+p%W2(3)*KnR(0,3)+(p%W2(4)*KnR(1,3)+p%W2(5)*KnR(1,1)+p%W2(6)*KnR(1,2))*d1r1
      T1=p%T1(1)*d1r1+(p%T1(2)*logr+p%T1(3))*r+p%T1(4)*KnR(1,1)+p%T1(5)*KnR(1,2)+(p%T1(6)*KnR(2,3)+p%T1(7)*KnR(2,1)&
        +p%T1(8)*KnR(2,2))*d1r1
      T2=p%T2(1)*d1r1+p%T2(2)*r+p%T2(3)*KnR(1,3)+p%T2(4)*KnR(1,1)+p%T2(5)*KnR(1,2)+(p%T2(6)*KnR(2,3)+p%T2(7)*KnR(2,1)&
        +p%T2(8)*KnR(2,2))*d1r1
      T3=p%T3(1)*d1r1+(p%T3(2)*logr+p%T3(3))*r+p%T3(4)*KnR(1,3)+(p%T3(5)*KnR(2,3)+p%T3(6)*KnR(2,1)+p%T3(7)*KnR(2,2))*d1r1
      ! Fundamental solutions
      fs_u(0,0)=eta
      fs_t(0,0)=W0*drdn
      fs_u(0,1:2)=vartheta*drdx
      fs_t(0,1:2)=T01*drdx*drdn+T02*n
      fs_u(1:2,0)=vartheta*drdx
      fs_t(1:2,0)=W1*drdx*drdn+W2*n
      do il=1,2
        do ik=1,2
          fs_u(il,ik)=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
          fs_t(il,ik)=T1*drdx(il)*n(ik)+T2*drdx(il)*drdx(ik)*drdn+T3*(drdx(ik)*n(il)+drdn*c_dkr(il,ik))
        end do
      end do
      ! Add the term of the numerical integration summation (fundamental solutions * functional shape functions * jacobian * weight)
      do il=0,2
        do ik=0,2
          h(:,il,ik)=h(:,il,ik)+fs_t(il,ik)*pphijw
          g(:,il,ik)=g(:,il,ik)+fs_u(il,ik)*sphijw
        end do
      end do
    end do
    ! Multiply by constants
    cte_u(0,0)=-c_1_2pi
    cte_t(0,0)=-c_1_2pi
    cte_u(0,1:2)=-c_1_2pi
    cte_t(0,1:2)= c_1_2pi
    cte_u(1:2,0)=-c_1_2pi/p%J
    cte_t(1:2,0)=-c_1_2pi/p%mu
    cte_u(1:2,1:2)=c_1_2pi/p%mu
    cte_t(1:2,1:2)=c_1_2pi
    do il=0,2
      do ik=0,2
        h(:,il,ik)=cte_t(il,ik)*h(:,il,ik)
        g(:,il,ik)=cte_u(il,ik)*g(:,il,ik)
      end do
    end do
    ! Reverse if needed
    if (reverse) h=-h
  end subroutine fbem_bem_harpor2d_sbie_ext_st

  recursive subroutine fbem_bem_harpor2d_sbie_ext_adp(e,reverse,xi_s,x_i,p,qsp,ks,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                     !! Element
    logical                            :: reverse               !! Reverse orientation
    real(kind=real64)                  :: xi_s(1,2)             !! Subdivision of the parent element
    real(kind=real64)                  :: x_i(2)                !! Collocation point position vector
    type(fbem_bem_harpor2d_parameters) :: p                     !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp                   !! Quasi-singular integration parameters
    integer                            :: ks                    !! Current level of subdivisions
    integer                            :: ns                    !! Maximum level of subdivision
    complex(kind=real64)               :: h(e%n_pnodes,0:2,0:2) !! h integration kernels matrix
    complex(kind=real64)               :: g(e%n_snodes,0:2,0:2) !! g integration kernels matrix
    ! Local
    integer              :: gln_near                  ! 1D Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer              :: gln                       ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical              :: subdivide                 ! True if subdivision has to be performed
    real(kind=real64)    :: barxi(1)                  ! Nearest element coordinate with respect to collocation point
    real(kind=real64)    :: barxip(1)                 ! Nearest element subdivision local coordinate with respect to collocation point
    real(kind=real64)    :: rmin                      ! Minimum distance between collocation point and barxi on the element
    real(kind=real64)    :: barr                      ! Telles jacobian at barxi
    real(kind=real64)    :: cl                        ! Characteristic length
    real(kind=real64)    :: d                         ! Normalized distance between collocation point and element subdivision
    integer              :: method                    ! Method used in nearest point algorithm
    real(kind=real64)    :: tmp_xi_s(1,2)             ! Subdivision
    real(kind=real64)    :: x_s(2,e%n_gnodes)         ! Coordinates of the element subdivision
    complex(kind=real64) :: h_tmp(e%n_pnodes,0:2,0:2) ! h integration kernels matrix (temporary)
    complex(kind=real64) :: g_tmp(e%n_snodes,0:2,0:2) ! g integration kernels matrix (temporary)
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
        call fbem_warning_message(error_unit,0,'fbem_bem_harpor2d_sbie_ext_adp',ns,'maximum number of subdivisions reached')
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
      call fbem_bem_harpor2d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_harpor2d_sbie_ext_adp(e,reverse,tmp_xi_s,x_i,p,qsp,ks+1,ns,h,g)
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_harpor2d_sbie_ext_st(e,reverse,xi_s,x_i,barxip,barr,p,gln,h_tmp,g_tmp)
      h=h+h_tmp
      g=g+g_tmp
    end if
  end subroutine fbem_bem_harpor2d_sbie_ext_adp

  !! This subroutine calculates the kernels for SBIE interior integration, needing only basic data.
  subroutine fbem_bem_harpor2d_sbie_int(ngp,type_g,type_f1,type_f2,delta_f,x_nodes,reverse,xi_i,p,h,g)
    implicit none
    ! I/O
    integer                            :: ngp                              !! Number of Gauss point to be used (<=32)
    integer                            :: type_g                           !! Geometrial interpolation
    integer                            :: type_f1                          !! Functional interpolation (primary variables)
    integer                            :: type_f2                          !! Functional interpolation (secondary variables)
    real(kind=real64)                  :: delta_f                          !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64)                  :: x_nodes(2,fbem_n_nodes(type_g))  !! Position vectors of geometrical nodes
    logical                            :: reverse                          !! Reverse normal vector
    real(kind=real64)                  :: xi_i                             !! Reference coordinate of the singular point.
    type(fbem_bem_harpor2d_parameters) :: p                                !! Parameters of the region
    complex(kind=real64)               :: h(fbem_n_nodes(type_f1),0:2,0:2) !! h integration kernels matrix
    complex(kind=real64)               :: g(fbem_n_nodes(type_f2),0:2,0:2) !! g integration kernels matrix
    ! Local
    integer                            :: l, k                             ! Counter variables for load direction and observation direction
    integer                            :: nnodes_g                         ! Number of nodes of the geometrical interpolation
    integer                            :: kphi                             ! Counter variable for shape functions loops
    integer                            :: kip                              ! Counter of integration points
    real(kind=real64)                  :: x_i(2)                           ! Real coordinates of collocation point
    real(kind=real64)                  :: phi_f1(fbem_n_nodes(type_f1))    ! Functional shape functions (primary variables) values at xi
    real(kind=real64)                  :: phi_f1_i(fbem_n_nodes(type_f1))  ! Functional shape functions (primary variables) values at xi_i
    real(kind=real64)                  :: phi_f2(fbem_n_nodes(type_f2))    ! Functional shape functions (secondary variables) values at xi
    real(kind=real64)                  :: phi_g(fbem_n_nodes(type_g))      ! Geometrical shape functions values at xi
    real(kind=real64)                  :: dphidxi_g(fbem_n_nodes(type_g))  ! Geometrical shape functions derivatives values at xi
    integer                            :: nsub                             ! Number of subdivision of the element
    integer                            :: ksub                             ! Counter of subdivision
    real(kind=real64)                  :: gamma                            ! Coordinate gamma
    real(kind=real64)                  :: w                                ! Weights of each integration point
    type(fbem_telles_parameters)       :: telles_parameters                ! Telles parameters
    real(kind=real64)                  :: jt                               ! Telles jacobian
    real(kind=real64)                  :: xip                              ! Coordinate  xip of subdivided element [0,1]
    real(kind=real64)                  :: js                               ! Jacobian of the xi [xisub(1,:),xisub[2,:]] -> xip [0,1] transformation
    real(kind=real64)                  :: xip_i(2)                         ! Singular point in xip space
    real(kind=real64)                  :: xi                               ! Coordinate xi
    real(kind=real64)                  :: xisub(2,2)                       ! Coordinates of subdivisions in xi space
    real(kind=real64)                  :: aux(10)                          ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)                  :: x(2)                             ! Position vector at xi
    real(kind=real64)                  :: T(2)                             ! Tangent vector at xi
    real(kind=real64)                  :: N(2)                             ! Normal vector at xi
    real(kind=real64)                  :: rv(2)                            ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)                  :: r, logr, d1r, d1r2               ! Distance vector module, log(r), 1/r and 1/r^2
    real(kind=real64)                  :: ra, rb                           ! Distance vector from collocation point to element vertices
    real(kind=real64)                  :: drdx(2)                          ! Distance vector derivatives with respect to x_k
    real(kind=real64)                  :: jg                               ! Geometric jacobian
    real(kind=real64)                  :: drdn                             ! Partial derivative of r respect to unit normal
    real(kind=real64)                  :: drdt                             ! Partial derivative of r respect to unit tangent
    real(kind=real64)                  :: jw                               ! Jacobians * weight
    real(kind=real64)                  :: phif1jw(fbem_n_nodes(type_f1))   ! Functional shape functions (primary variables) * jw
    real(kind=real64)                  :: phif2jw(fbem_n_nodes(type_f2))   ! Functional shape functions (secondary variables) * jw
    complex(kind=real64)               :: z1, z2, z3
    complex(kind=real64)               :: K0r_z1, K1r_z1, K2r_z1
    complex(kind=real64)               :: K0r_z2, K1r_z2, K2r_z2
    complex(kind=real64)               :: K0r_z3, K1r_z3, K2r_z3
    complex(kind=real64)               :: eta_r, vartheta_r, psi_r, chi_r
    complex(kind=real64)               :: W0_r, T01_r, T02_r, W1_r, W2_r
    complex(kind=real64)               :: T1_r, T2_r, T3_r
    complex(kind=real64)               :: fs_u(0:2,0:2), fs_t(0:2,0:2)
    complex(kind=real64)               :: cte_u(0:2,0:2), cte_t(0:2,0:2)
    !
    ! Initialize
    !
    ! The kernel matrices are set to zero.
    h=(0.d0,0.d0)
    g=(0.d0,0.d0)
    ! Number of geometrical nodes
    nnodes_g=fbem_n_nodes(type_g)
    ! If xi_i belongs to one of the vertices, subdivision is not needed.
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
        xisub(1,1)=-1.d0
        xisub(2,1)= 1.d0
        ! Coordinate xip of the collocation point
        ! If xi_i is at xi=-1
        if (xi_i.lt.0.d0) xip_i(1)=0.d0
        ! If xi_i is at xi=1
        if (xi_i.gt.0.d0) xip_i(1)=1.d0
      ! If 2 subdivisions
      case (2)
        ! Coordinates xi of the subdivision 1
        xisub(1,1)=-1.d0
        xisub(2,1)=xi_i
        ! Coordinate xip of the collocation point
        xip_i(1)=1.d0
        ! Coordinates xi of the subdivision 2
        xisub(1,2)=xi_i
        xisub(2,2)=1.d0
        ! Coordinate xip of the collocation point
        xip_i(2)=0.d0
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
    x_i=0.d0
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
    ! Integrate WEAKLY SINGULAR integrals
    !
    ! Loop through SUBDIVISIONS
    do ksub=1,nsub
      ! Jacobian of xip->xi transformation (is constant)
      js=xisub(2,ksub)-xisub(1,ksub)
      ! Telles transformation parameters (barr=0 at the collocation point)
      telles_parameters=fbem_telles01_calculate_parameters(xip_i(ksub),0.d0)
      ! Loop through INTEGRATION POINTS
      do kip=1,gl01_n(ngp)
        ! GAMMA->XIP->XI COORDINATE TRANSFORMATION
        ! Coordinate and weight in gamma [0,1]
        gamma=gl01_xi(kip,ngp)
        w=gl01_w(kip,ngp)
        ! Calculate xip and gamma->xip jacobian
        call fbem_telles_xi_and_jacobian(telles_parameters,gamma,xip,jt)
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
        jg=sqrt(T(1)**2+T(2)**2)
        ! Unit normal vector
        n=N/jg
        ! Distance vector
        rv=x-x_i
        ! Distance vector norm, log(r), 1/r and 1/r^2
        r=sqrt(rv(1)**2+rv(2)**2)
        logr=log(r)
        ! FUNDAMENTAL SOLUTIONS
        ! Fluid load / Fluid response
        fs_u(0,0)=p%eta(1)*logr
        fs_t(0,0)=0.d0
        ! Fluid load / Solid response
        fs_u(0,1:2)=0.d0
        fs_t(0,1:2)=p%T02(1)*n*logr
        ! Solid load / Fluid response
        fs_u(1:2,0)=0.d0
        fs_t(1:2,0)=p%W2(1)*n*logr
        ! Solid load / Solid response
        do l=1,2
          do k=1,2
            fs_u(l,k)=c_dkr(l,k)*p%psi(1)*logr
            fs_t(l,k)=0.d0
          end do
        end do
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
        ! FUNCTIONAL SHAPE FUNCTION * JACOBIANS * WEIGHT
        jw=jg*js*jt*w
        phif1jw=phi_f1*jw
        phif2jw=phi_f2*jw
        ! ADD INTEGRAND EVALUATION
        do l=0,2
          do k=0,2
            h(:,l,k)=h(:,l,k)+fs_t(l,k)*phif1jw
            g(:,l,k)=g(:,l,k)+fs_u(l,k)*phif2jw
          end do
        end do
      end do ! Loop through INTEGRATION POINTS
    end do ! Loop through SUBDIVISIONS
    !
    ! Integrate REGULAR integrals of the SINGULAR part
    !
    ! Loop through SUBDIVISIONS
    do ksub=1,nsub
      ! Jacobian of xip->xi transformation (is constant)
      js=xisub(2,ksub)-xisub(1,ksub)
      ! Loop through INTEGRATION POINTS
      do kip=1,gl01_n(ngp)
        ! XIP->XI COORDINATE TRANSFORMATION
        ! Coordinate and weight in xip [0,1]
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
        jg=dsqrt(T(1)**2+T(2)**2)
        ! Unit normal vector
        n=N/jg
        ! Unit tangent vector
        t=T/jg
        ! Distance vector
        rv=x-x_i
        ! Distance vector norm, log(r), 1/r and 1/r^2
        r=sqrt(rv(1)**2+rv(2)**2)
        logr=log(r)
        d1r=1.d0/r
        ! r_{,1} and r_{,2}
        drdx=rv*d1r
        ! dr/dn
        drdn=dot_product(drdx,n)
        ! dr/dGamma
        drdt=dot_product(drdx,t)
        ! FUNDAMENTAL SOLUTIONS
        ! Fluid load / Fluid response
        fs_u(0,0)=0.d0
        fs_t(0,0)=p%W0(1)*d1r*drdn
        ! Fluid load / Solid response
        fs_u(0,1:2)=0.d0
        fs_t(0,1:2)=0.d0
        ! Solid load / Fluid response
        fs_u(1:2,0)=0.d0
        fs_t(1:2,0)=0.d0
        ! Solid load / Solid response
        do l=1,2
          do k=1,2
            fs_u(l,k)=0.d0
            fs_t(l,k)=p%T2(1)*d1r*drdx(l)*drdx(k)*drdn+p%T3(1)*d1r*drdn*c_dkr(l,k)
          end do
        end do
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
        ! FUNCTIONAL SHAPE FUNCTION * JACOBIANS * WEIGHT
        jw=jg*js*w
        phif1jw=phi_f1*jw
        phif2jw=phi_f2*jw
        ! ADD INTEGRAND EVALUATION
        do l=0,2
          do k=0,2
            h(:,l,k)=h(:,l,k)+fs_t(l,k)*phif1jw
            g(:,l,k)=g(:,l,k)+fs_u(l,k)*phif2jw
          end do
        end do
        ! REGULAR PART of the CPV integral
        h(:,1,2)=h(:,1,2)-p%T1(1)*d1r*drdt*(phi_f1-phi_f1_i)*jw
        h(:,2,1)=h(:,2,1)+p%T1(1)*d1r*drdt*(phi_f1-phi_f1_i)*jw
      end do ! Loop through INTEGRATION POINTS
    end do ! Loop through SUBDIVISIONS
    !
    ! Integrate REGULAR integrals of the REGULAR part
    !
    ! Loop through SUBDIVISIONS
    do ksub=1,nsub
      ! Jacobian of xip->xi transformation (is constant)
      js=xisub(2,ksub)-xisub(1,ksub)
      ! Loop through INTEGRATION POINTS
      do kip=1,gl01_n(ngp)
        ! XIP->XI COORDINATE TRANSFORMATION
        ! Coordinate and weight in xip [0,1]
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
        jg=sqrt(T(1)**2+T(2)**2)
        ! Unit normal vector
        n=N/jg
        ! Unit tangent vector
        t=T/jg
        ! Distance vector
        rv=x-x_i
        ! Distance vector norm, log(r), 1/r and 1/r^2
        r=sqrt(rv(1)**2+rv(2)**2)
        logr=log(r)
        d1r=1.0d0/r
        d1r2=d1r**2
        ! r_{,k}
        drdx=rv*d1r
        ! dr/dn
        drdn=dot_product(drdx,n)
        ! dr/dGamma
        drdt=dot_product(drdx,t)
        ! z_m = i k_m r
        z1=c_im*p%k1*r
        z2=c_im*p%k2*r
        z3=c_im*p%k3*r
        ! Calculate K_n^R for z1, z2 and z3
        call fbem_modified_bessel_K0r_K1r_K2r_2(z1,K0r_z1,K1r_z1,K2r_z1)
        call fbem_modified_bessel_K0r_K1r_K2r_2(z2,K0r_z2,K1r_z2,K2r_z2)
        call fbem_modified_bessel_K0r_K1r_K2r_2(z3,K0r_z3,K1r_z3,K2r_z3)
        ! Calculate components (regular parts)
        eta_r=p%eta(2)+p%eta(3)*K0r_z1+p%eta(4)*K0r_z2
        vartheta_r=(p%vartheta(1)*logr+p%vartheta(2))*r+p%vartheta(3)*K1r_z1+p%vartheta(4)*K1r_z2
        psi_r=p%psi(2)+K0r_z3+(p%psi(3)*K1r_z3+p%psi(4)*K1r_z1+p%psi(5)*K1r_z2)*d1r
        chi_r=p%chi(1)+(p%chi(2)*logr+p%chi(3))*r**2+K2r_z3+p%chi(4)*K2r_z1+p%chi(5)*K2r_z2
        W0_r=(p%W0(2)*logr+p%W0(3))*r+p%W0(4)*K1r_z1+p%W0(5)*K1r_z2
        T01_r=p%T01(1)+p%T01(2)*K0r_z1+p%T01(3)*K0r_z2+(p%T01(4)*K1r_z1+p%T01(5)*K1r_z2)*d1r
        T02_r=p%T02(2)+p%T02(3)*K0r_z1+p%T02(4)*K0r_z2+(p%T02(5)*K1r_z1+p%T02(6)*K1r_z2)*d1r
        W1_r=p%W1(1)+(p%W1(2)*logr+p%W1(3))*r**2+p%W1(4)*K0r_z1+p%W1(5)*K0r_z2+(p%W1(6)*K1r_z1+p%W1(7)*K1r_z2)*d1r+p%W1(8)*K2r_z3&
           +p%W1(9)*K2r_z1+p%W1(10)*K2r_z2
        W2_r=p%W2(2)+p%W2(3)*K0r_z3+(p%W2(4)*K1r_z3+p%W2(5)*K1r_z1+p%W2(6)*K1r_z2)*d1r
        T1_r=(p%T1(2)*logr+p%T1(3))*r+p%T1(4)*K1r_z1+p%T1(5)*K1r_z2+(p%T1(6)*K2r_z3+p%T1(7)*K2r_z1+p%T1(8)*K2r_z2)*d1r
        T2_r=p%T2(2)*r+p%T2(3)*K1r_z3+p%T2(4)*K1r_z1+p%T2(5)*K1r_z2+(p%T2(6)*K2r_z3+p%T2(7)*K2r_z1+p%T2(8)*K2r_z2)*d1r
        T3_r=(p%T3(2)*logr+p%T3(3))*r+p%T3(4)*K1r_z3+(p%T3(5)*K2r_z3+p%T3(6)*K2r_z1+p%T3(7)*K2r_z2)*d1r
        ! FUNDAMENTAL SOLUTIONS
        ! Fluid load / Fluid response
        fs_u(0,0)=eta_r
        fs_t(0,0)=W0_r*drdn
        ! Fluid load / Solid response
        fs_u(0,1:2)=vartheta_r*drdx
        fs_t(0,1:2)=T01_r*drdx*drdn+T02_r*n
        ! Solid load / Fluid response
        fs_u(1:2,0)=vartheta_r*drdx
        fs_t(1:2,0)=W1_r*drdx*drdn+W2_r*n
        ! Solid load / Solid response
        do l=1,2
          do k=1,2
            fs_u(l,k)=psi_r*c_dkr(l,k)-chi_r*drdx(l)*drdx(k)
            fs_t(l,k)=T1_r*drdx(l)*n(k)+T2_r*drdx(l)*drdx(k)*drdn+T3_r*(drdx(k)*n(l)+drdn*c_dkr(l,k))
          end do
        end do
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
        ! FUNCTIONAL SHAPE FUNCTION * JACOBIANS * WEIGHT
        jw=jg*js*w
        phif1jw=phi_f1*jw
        phif2jw=phi_f2*jw
        ! ADD INTEGRAND EVALUATION
        do l=0,2
          do k=0,2
            h(:,l,k)=h(:,l,k)+fs_t(l,k)*phif1jw
            g(:,l,k)=g(:,l,k)+fs_u(l,k)*phif2jw
          end do
        end do
      end do ! Loop through INTEGRATION POINTS
    end do ! Loop through SUBDIVISIONS
    !
    ! Add analytical part of the CPV integral (H12 and H21)
    !
    ! Calculate ra and rb
    ra=sqrt((x_nodes(1,1)-x_i(1))**2+(x_nodes(2,1)-x_i(2))**2)
    rb=sqrt((x_nodes(1,2)-x_i(1))**2+(x_nodes(2,2)-x_i(2))**2)
    ! The addition depends on the position of xi_i
    select case (nsub)
      ! If just 1 subdivision, xi_i=-1 or xi_i=1
      case (1)
        ! If xi_i is at xi=-1
        if (xi_i.lt.0.d0) then
          h(:,1,2)=h(:,1,2)-p%T1(1)*phi_f1_i*log(rb)
          h(:,2,1)=h(:,2,1)+p%T1(1)*phi_f1_i*log(rb)
        ! If xi_i is at xi=1
        else
          h(:,1,2)=h(:,1,2)+p%T1(1)*phi_f1_i*log(ra)
          h(:,2,1)=h(:,2,1)-p%T1(1)*phi_f1_i*log(ra)
        end if
      ! If 2 subdivisions, xi_i is inside de element -1<xi_i<1
      case (2)
        h(:,1,2)=h(:,1,2)-p%T1(1)*phi_f1_i*(log(rb)-log(ra))
        h(:,2,1)=h(:,2,1)+p%T1(1)*phi_f1_i*(log(rb)-log(ra))
    end select
    !
    ! Calculate cte_t and cte_u, and multiply them to h and g, respectively
    !
    ! Fluid load / Fluid response constants
    cte_u(0,0)=-c_1_2pi
    cte_t(0,0)=-c_1_2pi
    ! Fluid load / Solid response constants
    cte_u(0,1:2)=-c_1_2pi
    cte_t(0,1:2)= c_1_2pi
    ! Solid load / Fluid response constants
    cte_u(1:2,0)=-c_1_2pi/p%J
    cte_t(1:2,0)=-c_1_2pi/p%mu
    ! Solid load / Solid response constants
    cte_u(1:2,1:2)=c_1_2pi/p%mu
    cte_t(1:2,1:2)=c_1_2pi
    ! Multiply by constants
    do l=0,2
      do k=0,2
        h(:,l,k)=cte_t(l,k)*h(:,l,k)
        g(:,l,k)=cte_u(l,k)*g(:,l,k)
      end do
    end do
    !
    ! If the element orientation is reversed
    !
    if (reverse) h=-h
  end subroutine fbem_bem_harpor2d_sbie_int

  subroutine fbem_bem_harpor2d_sbie_auto(e,reverse,x_i,p,qsp,ns,h,g)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                     !! Integration element
    logical                            :: reverse               !! Reverse orientation
    real(kind=real64)                  :: x_i(2)                !! Collocation point
    type(fbem_bem_harpor2d_parameters) :: p                     !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp                   !! Quasi-singular integration parameters
    integer                            :: ns                    !! Maximum level of subdivisions
    complex(kind=real64)               :: h(e%n_pnodes,0:2,0:2) !! h integration kernel
    complex(kind=real64)               :: g(e%n_snodes,0:2,0:2) !! g integration kernel
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
        call fbem_bem_harpor2d_sbie_int(30,e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi(1),p,h,g)
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
          call fbem_bem_harpor2d_sbie_ext_pre(ps,e,reverse,x_i,p,h,g)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_harpor2d_sbie_ext_adp(e,reverse,xi_s,x_i,p,qsp,1,ns,h,g)
        end if
    end select
  end subroutine fbem_bem_harpor2d_sbie_auto

  ! ================================================================================================================================

  ! ================================================================================================================================
  ! HYPERSINGULAR BOUNDARY INTEGRAL EQUATION (HBIE)

  subroutine fbem_bem_harpor2d_hbie_ext_pre(ps,e,reverse,x_i,n_i,p,m,l)
    implicit none
    ! I/O
    integer                            :: ps                    !! Selected precalculated dataset
    type(fbem_bem_element)             :: e                     !! Element
    logical                            :: reverse               !! Reverse normal vector
    real(kind=real64)                  :: x_i(2)                !! Collocation point position vector
    real(kind=real64)                  :: n_i(2)                !! Collocation point unit normal vector
    type(fbem_bem_harpor2d_parameters) :: p                     !! Parameters of the region
    complex(kind=real64)               :: m(e%n_pnodes,0:2,0:2) !! m kernels matrix
    complex(kind=real64)               :: l(e%n_snodes,0:2,0:2) !! l kernels matrix
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
    complex(kind=real64) :: z(3)                           ! Arguments z=ikr
    complex(kind=real64) :: KnR(0:2,3)                     ! Bessel functions decomposition
    complex(kind=real64) :: W0, T01, T02, W1, W2           ! Fundamental solution components
    complex(kind=real64) :: T1, T2, T3                     ! Fundamental solution components
    complex(kind=real64) :: Q01, Q02, S01, S02, S03        ! Fundamental solution components
    complex(kind=real64) :: S1, S2, S3, S4, S5             ! Fundamental solution components
    complex(kind=real64) :: fs_d(0:2,0:2), fs_s(0:2,0:2)   ! Fundamental solutions
    complex(kind=real64) :: cte_d(0:2,0:2), cte_s(0:2,0:2) ! Fundamental solution constants
    ! Initialisation
    m=(0.d0,0.d0)
    l=(0.d0,0.d0)
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
      z(3)=c_im*p%k3*r
      call fbem_BesselKnR_decomposed(3,z,KnR)
      W0=p%W0(1)*d1r1+(p%W0(2)*logr+p%W0(3))*r+p%W0(4)*KnR(1,1)+p%W0(5)*KnR(1,2)
      T01=p%T01(1)+p%T01(2)*KnR(0,1)+p%T01(3)*KnR(0,2)+(p%T01(4)*KnR(1,1)+p%T01(5)*KnR(1,2))*d1r1
      T02=p%T02(1)*logr+p%T02(2)+p%T02(3)*KnR(0,1)+p%T02(4)*KnR(0,2)+(p%T02(5)*KnR(1,1)+p%T02(6)*KnR(1,2))*d1r1
      W1=p%W1(1)+(p%W1(2)*logr+p%W1(3))*r2+p%W1(4)*KnR(0,1)+p%W1(5)*KnR(0,2)+(p%W1(6)*KnR(1,1)+p%W1(7)*KnR(1,2))*d1r1&
        +p%W1(8)*KnR(2,3)+p%W1(9)*KnR(2,1)+p%W1(10)*KnR(2,2)
      W2=p%W2(1)*logr+p%W2(2)+p%W2(3)*KnR(0,3)+(p%W2(4)*KnR(1,3)+p%W2(5)*KnR(1,1)+p%W2(6)*KnR(1,2))*d1r1
      T1=p%T1(1)*d1r1+(p%T1(2)*logr+p%T1(3))*r+p%T1(4)*KnR(1,1)+p%T1(5)*KnR(1,2)&
        +(p%T1(6)*KnR(2,3)+p%T1(7)*KnR(2,1)+p%T1(8)*KnR(2,2))*d1r1
      T2=p%T2(1)*d1r1+p%T2(2)*r+p%T2(3)*KnR(1,3)+p%T2(4)*KnR(1,1)+p%T2(5)*KnR(1,2)&
        +(p%T2(6)*KnR(2,3)+p%T2(7)*KnR(2,1)+p%T2(8)*KnR(2,2))*d1r1
      T3=p%T3(1)*d1r1+(p%T3(2)*logr+p%T3(3))*r+p%T3(4)*KnR(1,3)+(p%T3(5)*KnR(2,3)+p%T3(6)*KnR(2,1)+p%T3(7)*KnR(2,2))*d1r1
      Q01=p%Q01(1)*d1r2+p%Q01(2)+(p%Q01(3)*logr+p%Q01(4))*r2+p%Q01(5)*KnR(0,1)+p%Q01(6)*KnR(0,2)&
          +(p%Q01(7)*KnR(1,1)+p%Q01(8)*KnR(1,2))*d1r1+p%Q01(9)*KnR(2,3)+p%Q01(10)*KnR(2,1)+p%Q01(11)*KnR(2,2)
      Q02=p%Q02(1)*d1r2+p%Q02(2)*logr+p%Q02(3)+p%Q02(4)*KnR(0,3)+(p%Q02(5)*KnR(1,3)+p%Q02(6)*KnR(1,1)+p%Q02(7)*KnR(1,2))*d1r1
      S01=p%S01(1)*d1r1+(p%S01(2)*logr+p%S01(3))*r+(p%S01(4)*KnR(0,1)+p%S01(5)*KnR(0,2))*d1r1+(p%S01(6)*KnR(1,1)+p%S01(7)*KnR(1,2))*d1r2&
          +p%S01(8)*KnR(1,3)+p%S01(9)*KnR(1,1)+p%S01(10)*KnR(1,2)+(p%S01(11)*KnR(2,3)+p%S01(12)*KnR(2,1)+p%S01(13)*KnR(2,2))*d1r1
      S02=p%S02(1)*d1r1+(p%S02(2)*logr+p%S02(3))*r+(p%S02(4)*KnR(0,1)+p%S02(5)*KnR(0,2))*d1r1+(p%S02(6)*KnR(1,1)+p%S02(7)*KnR(1,2))*d1r2&
          +p%S02(8)*KnR(1,1)+p%S02(9)*KnR(1,2)+(p%S02(10)*KnR(2,3)+p%S02(11)*KnR(2,1)+p%S02(12)*KnR(2,2))*d1r1
      S03=p%S03(1)*d1r1+(p%S03(2)*logr+p%S03(3))*r+(p%S03(4)*KnR(0,1)+p%S03(5)*KnR(0,2))*d1r1+(p%S03(6)*KnR(1,1)+p%S03(7)*KnR(1,2))*d1r2&
          +p%S03(8)*KnR(1,3)+(p%S03(9)*KnR(2,3)+p%S03(10)*KnR(2,1)+p%S03(11)*KnR(2,2))*d1r1
      S1=p%S1(1)*d1r2+p%S1(2)+p%S1(3)*KnR(0,3)+(p%S1(4)*KnR(1,3)+p%S1(5)*KnR(1,1)+p%S1(6)*KnR(1,2))*d1r1&
         +(p%S1(7)*KnR(2,3)+p%S1(8)*KnR(2,1)+p%S1(9)*KnR(2,2))*d1r2
      S2=p%S2(1)*d1r2+p%S2(2)+p%S2(3)*KnR(0,1)+p%S2(4)*KnR(0,2)+(p%S2(5)*KnR(1,3)+p%S2(6)*KnR(1,1)+p%S2(7)*KnR(1,2))*d1r1&
         +(p%S2(8)*KnR(2,3)+p%S2(9)*KnR(2,1)+p%S2(10)*KnR(2,2))*d1r2
      S3=p%S3(1)*d1r2+p%S3(2)+p%S3(3)*KnR(0,3)+p%S3(4)*KnR(0,1)+p%S3(5)*KnR(0,2)+(p%S3(6)*KnR(1,3)+p%S3(7)*KnR(1,1)+p%S3(8)*KnR(1,2))*d1r1&
         +(p%S3(9)*KnR(2,3)+p%S3(10)*KnR(2,1)+p%S3(11)*KnR(2,2))*d1r2
      S4=p%S4(1)*d1r2+p%S4(2)*logr+p%S4(3)+p%S4(4)*KnR(1,3)*d1r1+(p%S4(5)*KnR(2,3)+p%S4(6)*KnR(2,1)+p%S4(7)*KnR(2,2))*d1r2
      S5=p%S5(1)*d1r2+p%S5(2)*logr+p%S5(3)+p%S5(4)*KnR(0,1)+p%S5(5)*KnR(0,2)+(p%S5(6)*KnR(1,1)+p%S5(7)*KnR(1,2))*d1r1&
         +(p%S5(8)*KnR(2,3)+p%S5(9)*KnR(2,1)+p%S5(10)*KnR(2,2))*d1r2
      ! Fundamental solutions
      fs_d(0,0)=W0*drdni
      fs_s(0,0)=Q01*drdn*drdni+Q02*n_dot_ni
      fs_d(0,1:2)=-W1*drdx*drdni+W2*n_i
      fs_s(0,1:2)=S01*drdx*drdn*drdni+S02*n*drdni+S03*(n_i*drdn+drdx*n_dot_ni)
      fs_d(1:2,0)=-T01*drdx*drdni+T02*n_i
      fs_s(1:2,0)=-S01*drdx*drdn*drdni+S02*n_i*drdn-S03*(-n*drdni+drdx*n_dot_ni)
      do il=1,2
        do ik=1,2
          fs_d(il,ik)=-T1*drdx(ik)*n_i(il)+T2*drdx(il)*drdx(ik)*drdni-T3*(drdx(il)*n_i(ik)-drdni*c_dkr(il,ik))
          fs_s(il,ik)=S1*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni-drdn*drdni*c_dkr(il,ik)+drdx(il)*drdx(ik)*n_dot_ni)&
                     +S2*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni)+S3*drdx(il)*drdx(ik)*drdn*drdni&
                     +S4*(n_dot_ni*c_dkr(il,ik)+n(il)*n_i(ik))+S5*n(ik)*n_i(il)
        end do
      end do
      ! Add the term of the numerical integration summation (fundamental solutions * functional shape functions * jacobian * weight)
      do il=0,2
        do ik=0,2
          m(:,il,ik)=m(:,il,ik)+fs_s(il,ik)*pphijw
          l(:,il,ik)=l(:,il,ik)+fs_d(il,ik)*sphijw
        end do
      end do
    end do
    ! Multiply by constants
    cte_d(0,0)=-c_1_2pi/p%J
    cte_s(0,0)=-c_1_2pi
    cte_d(0,1:2)=c_1_2pi/p%mu
    cte_s(0,1:2)=c_1_2pi
    cte_d(1:2,0)=-c_1_2pi/p%J
    cte_s(1:2,0)=-c_1_2pi
    cte_d(1:2,1:2)=c_1_2pi
    cte_s(1:2,1:2)=c_1_2pi
    do il=0,2
      do ik=0,2
        m(:,il,ik)=cte_s(il,ik)*m(:,il,ik)
        l(:,il,ik)=cte_d(il,ik)*l(:,il,ik)
      end do
    end do
    ! Reverse if needed
    if (reverse) m=-m
  end subroutine fbem_bem_harpor2d_hbie_ext_pre

  !! This subroutine calculates the kernels for HBIE exterior integration (near collocation points) using Telles transformation
  !! within a subdivision of the element, needing only needs basic data.
  subroutine fbem_bem_harpor2d_hbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,p,gln,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                     !! Integration element
    logical                            :: reverse               !! Reverse normal vector
    real(kind=real64)                  :: xi_s(1,2)             !! Coordinates of the subdivision of the parent element (xi space [-1,1])
    real(kind=real64)                  :: x_i(2)                !! Collocation point position vector
    real(kind=real64)                  :: n_i(2)                !! Collocation point unit normal vector
    real(kind=real64)                  :: barxip(1)             !! Nearest local coordinate of the subdivision with respect to x_i
    real(kind=real64)                  :: barr                  !! Telles jacobian at barxip
    type(fbem_bem_harpor2d_parameters) :: p                     !! Parameters of the region
    integer                            :: gln                   !! 1D Gauss-Legendre number of integration points (<=32)
    complex(kind=real64)               :: m(e%n_pnodes,0:2,0:2) !! m kernels matrix
    complex(kind=real64)               :: l(e%n_snodes,0:2,0:2) !! l kernels matrix
    ! Local
    integer                      :: kphi                           ! Counter variable for shape functions loops
    integer                      :: kip                            ! Counter variable of integration points
    real(kind=real64)            :: gamma                          ! Coordinate gamma (Telles transformation space [0,1])
    real(kind=real64)            :: w                              ! Weights of an integration point
    type(fbem_telles_parameters) :: telles_parameters              ! Telles parameters
    real(kind=real64)            :: jt                             ! Telles jacobian: xip->gamma
    real(kind=real64)            :: xip                            ! Coordinate xip (subdivision space [0,1])
    real(kind=real64)            :: js                             ! Subdivision jacobian: xi->xip
    real(kind=real64)            :: xi                             ! Coordinate xi [xi1,xi2]
    real(kind=real64)            :: aux(10)                        ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: gphi(e%n_gnodes)               ! Geometrical shape functions values
    real(kind=real64)            :: dgphidxi(e%n_gnodes)           ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: pphi(e%n_pnodes)               ! Functional shape functions values
    real(kind=real64)            :: sphi(e%n_snodes)               ! Functional shape functions values
    real(kind=real64)            :: x(2)                           ! Position vector at xi
    real(kind=real64)            :: T(2)                           ! Tangent vector at xi
    real(kind=real64)            :: N(2)                           ! Normal vector at xi
    real(kind=real64)            :: rv(2)                          ! Radiovector between collocation point and integration point (x-x_i)
    real(kind=real64)            :: r, r2, d1r1, d1r2, logr        ! Radiovector module, squared, inverses and log(r)
    real(kind=real64)            :: drdx(2)                        ! Radiovector derivatives with respect to x_k
    real(kind=real64)            :: jg                             ! Geometric jacobian
    real(kind=real64)            :: drdn                           ! Partial derivative of r respect to unit normal
    real(kind=real64)            :: drdni                          ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)            :: n_dot_ni                       ! Dot product of unit normals
    real(kind=real64)            :: jw                             ! Jacobians * weights
    real(kind=real64)            :: pphijw(e%n_pnodes)             ! Auxiliary variables for integrand evaluation
    real(kind=real64)            :: sphijw(e%n_snodes)             ! Auxiliary variables for integrand evaluation
    integer                      :: il, ik                         ! Counter variables for load direction and observation direction
    complex(kind=real64)         :: z(3)                           ! Bessel functions arguments: z=ikr
    complex(kind=real64)         :: KnR(0:2,3)                     ! Bessel functions decomposition
    complex(kind=real64)         :: W0, T01, T02, W1, W2           ! Fundamental solution components
    complex(kind=real64)         :: T1, T2, T3                     ! Fundamental solution components
    complex(kind=real64)         :: Q01, Q02, S01, S02, S03        ! Fundamental solution components
    complex(kind=real64)         :: S1, S2, S3, S4, S5             ! Fundamental solution components
    complex(kind=real64)         :: fs_d(0:2,0:2), fs_s(0:2,0:2)   ! Fundamental solutions
    complex(kind=real64)         :: cte_d(0:2,0:2), cte_s(0:2,0:2) ! Fundamental solution constants
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
      z(3)=c_im*p%k3*r
      call fbem_BesselKnR_decomposed(3,z,KnR)
      W0=p%W0(1)*d1r1+(p%W0(2)*logr+p%W0(3))*r+p%W0(4)*KnR(1,1)+p%W0(5)*KnR(1,2)
      T01=p%T01(1)+p%T01(2)*KnR(0,1)+p%T01(3)*KnR(0,2)+(p%T01(4)*KnR(1,1)+p%T01(5)*KnR(1,2))*d1r1
      T02=p%T02(1)*logr+p%T02(2)+p%T02(3)*KnR(0,1)+p%T02(4)*KnR(0,2)+(p%T02(5)*KnR(1,1)+p%T02(6)*KnR(1,2))*d1r1
      W1=p%W1(1)+(p%W1(2)*logr+p%W1(3))*r2+p%W1(4)*KnR(0,1)+p%W1(5)*KnR(0,2)+(p%W1(6)*KnR(1,1)+p%W1(7)*KnR(1,2))*d1r1&
        +p%W1(8)*KnR(2,3)+p%W1(9)*KnR(2,1)+p%W1(10)*KnR(2,2)
      W2=p%W2(1)*logr+p%W2(2)+p%W2(3)*KnR(0,3)+(p%W2(4)*KnR(1,3)+p%W2(5)*KnR(1,1)+p%W2(6)*KnR(1,2))*d1r1
      T1=p%T1(1)*d1r1+(p%T1(2)*logr+p%T1(3))*r+p%T1(4)*KnR(1,1)+p%T1(5)*KnR(1,2)&
        +(p%T1(6)*KnR(2,3)+p%T1(7)*KnR(2,1)+p%T1(8)*KnR(2,2))*d1r1
      T2=p%T2(1)*d1r1+p%T2(2)*r+p%T2(3)*KnR(1,3)+p%T2(4)*KnR(1,1)+p%T2(5)*KnR(1,2)&
        +(p%T2(6)*KnR(2,3)+p%T2(7)*KnR(2,1)+p%T2(8)*KnR(2,2))*d1r1
      T3=p%T3(1)*d1r1+(p%T3(2)*logr+p%T3(3))*r+p%T3(4)*KnR(1,3)+(p%T3(5)*KnR(2,3)+p%T3(6)*KnR(2,1)+p%T3(7)*KnR(2,2))*d1r1
      Q01=p%Q01(1)*d1r2+p%Q01(2)+(p%Q01(3)*logr+p%Q01(4))*r2+p%Q01(5)*KnR(0,1)+p%Q01(6)*KnR(0,2)&
          +(p%Q01(7)*KnR(1,1)+p%Q01(8)*KnR(1,2))*d1r1+p%Q01(9)*KnR(2,3)+p%Q01(10)*KnR(2,1)+p%Q01(11)*KnR(2,2)
      Q02=p%Q02(1)*d1r2+p%Q02(2)*logr+p%Q02(3)+p%Q02(4)*KnR(0,3)+(p%Q02(5)*KnR(1,3)+p%Q02(6)*KnR(1,1)+p%Q02(7)*KnR(1,2))*d1r1
      S01=p%S01(1)*d1r1+(p%S01(2)*logr+p%S01(3))*r+(p%S01(4)*KnR(0,1)+p%S01(5)*KnR(0,2))*d1r1+(p%S01(6)*KnR(1,1)+p%S01(7)*KnR(1,2))*d1r2&
          +p%S01(8)*KnR(1,3)+p%S01(9)*KnR(1,1)+p%S01(10)*KnR(1,2)+(p%S01(11)*KnR(2,3)+p%S01(12)*KnR(2,1)+p%S01(13)*KnR(2,2))*d1r1
      S02=p%S02(1)*d1r1+(p%S02(2)*logr+p%S02(3))*r+(p%S02(4)*KnR(0,1)+p%S02(5)*KnR(0,2))*d1r1+(p%S02(6)*KnR(1,1)+p%S02(7)*KnR(1,2))*d1r2&
          +p%S02(8)*KnR(1,1)+p%S02(9)*KnR(1,2)+(p%S02(10)*KnR(2,3)+p%S02(11)*KnR(2,1)+p%S02(12)*KnR(2,2))*d1r1
      S03=p%S03(1)*d1r1+(p%S03(2)*logr+p%S03(3))*r+(p%S03(4)*KnR(0,1)+p%S03(5)*KnR(0,2))*d1r1+(p%S03(6)*KnR(1,1)+p%S03(7)*KnR(1,2))*d1r2&
          +p%S03(8)*KnR(1,3)+(p%S03(9)*KnR(2,3)+p%S03(10)*KnR(2,1)+p%S03(11)*KnR(2,2))*d1r1
      S1=p%S1(1)*d1r2+p%S1(2)+p%S1(3)*KnR(0,3)+(p%S1(4)*KnR(1,3)+p%S1(5)*KnR(1,1)+p%S1(6)*KnR(1,2))*d1r1&
         +(p%S1(7)*KnR(2,3)+p%S1(8)*KnR(2,1)+p%S1(9)*KnR(2,2))*d1r2
      S2=p%S2(1)*d1r2+p%S2(2)+p%S2(3)*KnR(0,1)+p%S2(4)*KnR(0,2)+(p%S2(5)*KnR(1,3)+p%S2(6)*KnR(1,1)+p%S2(7)*KnR(1,2))*d1r1&
         +(p%S2(8)*KnR(2,3)+p%S2(9)*KnR(2,1)+p%S2(10)*KnR(2,2))*d1r2
      S3=p%S3(1)*d1r2+p%S3(2)+p%S3(3)*KnR(0,3)+p%S3(4)*KnR(0,1)+p%S3(5)*KnR(0,2)+(p%S3(6)*KnR(1,3)+p%S3(7)*KnR(1,1)+p%S3(8)*KnR(1,2))*d1r1&
         +(p%S3(9)*KnR(2,3)+p%S3(10)*KnR(2,1)+p%S3(11)*KnR(2,2))*d1r2
      S4=p%S4(1)*d1r2+p%S4(2)*logr+p%S4(3)+p%S4(4)*KnR(1,3)*d1r1+(p%S4(5)*KnR(2,3)+p%S4(6)*KnR(2,1)+p%S4(7)*KnR(2,2))*d1r2
      S5=p%S5(1)*d1r2+p%S5(2)*logr+p%S5(3)+p%S5(4)*KnR(0,1)+p%S5(5)*KnR(0,2)+(p%S5(6)*KnR(1,1)+p%S5(7)*KnR(1,2))*d1r1&
         +(p%S5(8)*KnR(2,3)+p%S5(9)*KnR(2,1)+p%S5(10)*KnR(2,2))*d1r2
      ! Fundamental solutions
      fs_d(0,0)=W0*drdni
      fs_s(0,0)=Q01*drdn*drdni+Q02*n_dot_ni
      fs_d(0,1:2)=-W1*drdx*drdni+W2*n_i
      fs_s(0,1:2)=S01*drdx*drdn*drdni+S02*n*drdni+S03*(n_i*drdn+drdx*n_dot_ni)
      fs_d(1:2,0)=-T01*drdx*drdni+T02*n_i
      fs_s(1:2,0)=-S01*drdx*drdn*drdni+S02*n_i*drdn-S03*(-n*drdni+drdx*n_dot_ni)
      do il=1,2
        do ik=1,2
          fs_d(il,ik)=-T1*drdx(ik)*n_i(il)+T2*drdx(il)*drdx(ik)*drdni-T3*(drdx(il)*n_i(ik)-drdni*c_dkr(il,ik))
          fs_s(il,ik)=S1*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni-drdn*drdni*c_dkr(il,ik)+drdx(il)*drdx(ik)*n_dot_ni)&
                     +S2*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni)+S3*drdx(il)*drdx(ik)*drdn*drdni&
                     +S4*(n_dot_ni*c_dkr(il,ik)+n(il)*n_i(ik))+S5*n(ik)*n_i(il)
        end do
      end do
      ! Add the term of the numerical integration summation (fundamental solutions * functional shape functions * jacobian * weight)
      do il=0,2
        do ik=0,2
          m(:,il,ik)=m(:,il,ik)+fs_s(il,ik)*pphijw
          l(:,il,ik)=l(:,il,ik)+fs_d(il,ik)*sphijw
        end do
      end do
    end do
    ! Multiply by constants
    cte_d(0,0)=-c_1_2pi/p%J
    cte_s(0,0)=-c_1_2pi
    cte_d(0,1:2)=c_1_2pi/p%mu
    cte_s(0,1:2)=c_1_2pi
    cte_d(1:2,0)=-c_1_2pi/p%J
    cte_s(1:2,0)=-c_1_2pi
    cte_d(1:2,1:2)=c_1_2pi
    cte_s(1:2,1:2)=c_1_2pi
    do il=0,2
      do ik=0,2
        m(:,il,ik)=cte_s(il,ik)*m(:,il,ik)
        l(:,il,ik)=cte_d(il,ik)*l(:,il,ik)
      end do
    end do
    ! Reverse if needed
    if (reverse) m=-m
  end subroutine fbem_bem_harpor2d_hbie_ext_st

  recursive subroutine fbem_bem_harpor2d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,p,qsp,ks,ns,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                     !! Element
    logical                            :: reverse               !! Reverse orientation
    real(kind=real64)                  :: xi_s(1,2)             !! Subdivision of the parent element
    real(kind=real64)                  :: x_i(2)                !! Collocation point position vector
    real(kind=real64)                  :: n_i(2)                !! Collocation point unit normal vector
    type(fbem_bem_harpor2d_parameters) :: p                     !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp                   !! Quasi-singular integration parameters
    integer                            :: ks                    !! Current level of subdivisions
    integer                            :: ns                    !! Maximum level of subdivision
    complex(kind=real64)               :: m(e%n_pnodes,0:2,0:2) !! m integration kernels matrix
    complex(kind=real64)               :: l(e%n_snodes,0:2,0:2) !! l integration kernels matrix
    ! Local
    integer              :: gln_near                  ! 1D Gauss-Legendre integ. points required to integrate only the quasi-singular integrand
    integer              :: gln                       ! 1D Gauss-Legendre integ. points required to integrate the whole integrand
    logical              :: subdivide                 ! True if subdivision has to be performed
    real(kind=real64)    :: barxi(1)                  ! Nearest element coordinate with respect to collocation point
    real(kind=real64)    :: barxip(1)                 ! Nearest element subdivision local coordinate with respect to collocation point
    real(kind=real64)    :: rmin                      ! Minimum distance between collocation point and barxi on the element
    real(kind=real64)    :: barr                      ! Telles jacobian at barxi
    real(kind=real64)    :: cl                        ! Characteristic length
    real(kind=real64)    :: d                         ! Normalized distance between collocation point and element subdivision
    integer              :: method                    ! Method used in nearest point algorithm
    real(kind=real64)    :: tmp_xi_s(1,2)             ! Subdivision
    real(kind=real64)    :: x_s(2,e%n_gnodes)         ! Coordinates of the element subdivision
    complex(kind=real64) :: m_tmp(e%n_pnodes,0:2,0:2) ! m kernels matrix (temporary)
    complex(kind=real64) :: l_tmp(e%n_snodes,0:2,0:2) ! l kernels matrix (temporary)
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
        call fbem_warning_message(error_unit,0,'fbem_bem_harpor2d_hbie_ext_adp',ns,'maximum number of subdivisions reached')
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
      call fbem_bem_harpor2d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
      ! SUBLINE 2
      tmp_xi_s(1,1)=0.5d0*(xi_s(1,1)+xi_s(1,2))
      tmp_xi_s(1,2)=xi_s(1,2)
      call fbem_bem_harpor2d_hbie_ext_adp(e,reverse,tmp_xi_s,x_i,n_i,p,qsp,ks+1,ns,m,l)
    ! Calculate the subdivided element using Telles transformation
    else
      barr=fbem_telles_barr(d,fbem_f_any)
      gln=max(gln_near,e%gln_far)
      call fbem_bem_harpor2d_hbie_ext_st(e,reverse,xi_s,x_i,n_i,barxip,barr,p,gln,m_tmp,l_tmp)
      m=m+m_tmp
      l=l+l_tmp
    end if
  end subroutine fbem_bem_harpor2d_hbie_ext_adp

  !! This subroutine calculates the kernels for HBIE interior integration, needing only basic data.
  subroutine fbem_bem_harpor2d_hbie_int(ngp,type_g,type_f1,type_f2,delta_f,x_nodes,reverse,xi_i,p,m,l)
    implicit none
    ! I/O
    integer                            :: ngp                              !! Number of Gauss point to be used (<=32)
    integer                            :: type_g                           !! Geometrial interpolation
    integer                            :: type_f1                          !! Functional interpolation (primary variables)
    integer                            :: type_f2                          !! Functional interpolation (secondary variables)
    real(kind=real64)                  :: delta_f                          !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64)                  :: x_nodes(2,fbem_n_nodes(type_g))  !! Position vectors of geometrical nodes
    logical                            :: reverse                          !! Reverse normal vector
    real(kind=real64)                  :: xi_i                             !! Reference coordinate of the singular point.
    type(fbem_bem_harpor2d_parameters) :: p                                !! Parameters of the region
    complex(kind=real64)               :: m(fbem_n_nodes(type_f1),0:2,0:2) !! m integration kernels matrix
    complex(kind=real64)               :: l(fbem_n_nodes(type_f2),0:2,0:2) !! l integration kernels matrix
    ! Local

    integer                            :: il, ik                              ! Counter variables for load direction and observation direction
    integer                            :: nnodes_g                            ! Number of nodes of the geometrical interpolation
    integer                            :: kphi                                ! Counter variable for shape functions loops
    integer                            :: kip                                 ! Counter of integration points
    real(kind=real64)                  :: x_i(2)                              ! Real coordinates of collocation point
    real(kind=real64)                  :: t_i(2)                              ! Unit tangent at collocation point
    real(kind=real64)                  :: n_i(2)                              ! Unit normal at collocation point
    real(kind=real64)                  :: phi_f1(fbem_n_nodes(type_f1))       ! Functional shape functions (primary variables) values at xi
    real(kind=real64)                  :: phi_f1_i(fbem_n_nodes(type_f1))     ! Functional shape functions (primary variables) values at xi_i
    real(kind=real64)                  :: dphidxi_f1_i(fbem_n_nodes(type_f1)) ! Functional shape functions derivatives (primary variables) values at xi_i
    real(kind=real64)                  :: phi_f2(fbem_n_nodes(type_f2))       ! Functional shape functions (primary variables) values at xi
    real(kind=real64)                  :: phi_f2_i(fbem_n_nodes(type_f2))     ! Functional shape functions (primary variables) values at xi_i
    real(kind=real64)                  :: dphidxi_f2_i(fbem_n_nodes(type_f2)) ! Functional shape functions derivatives (primary variables) values at xi_i
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
    real(kind=real64)                  :: rv(2)                               ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)                  :: r, logr, d1r, d1r2                  ! Distance vector module
    real(kind=real64)                  :: ra, rb                              ! Distance vector from collocation point to element vertices
    real(kind=real64)                  :: drdx(2)                             ! Distance vector derivatives with respect to x_k
    real(kind=real64)                  :: jg                                  ! Geometric jacobian
    real(kind=real64)                  :: jgi                                 ! Geometric jacobian at the collocation point
    real(kind=real64)                  :: drdn                                ! Partial derivative of r respect to unit normal
    real(kind=real64)                  :: drdni                               ! Partial derivative of r respect to unit normal at collocation point
    real(kind=real64)                  :: n_dot_ni                            ! Dot product of unit normals
    real(kind=real64)                  :: drdt                                ! Partial derivative of r respect to unit tangent
    real(kind=real64)                  :: jw                                  ! Jacobians * weight
    real(kind=real64)                  :: phif1jw(fbem_n_nodes(type_f1))      ! Auxiliary variables for integrand evaluation
    real(kind=real64)                  :: phif2jw(fbem_n_nodes(type_f2))      ! Auxiliary variables for integrand evaluation
    type(fbem_telles_parameters)       :: telles_parameters                   ! Telles parameters
    real(kind=real64)                  :: jt                                  ! Telles jacobian
    complex(kind=real64)               :: mr(fbem_n_nodes(type_f1),0:2,0:2)   ! m integration kernels matrix (regular part)
    complex(kind=real64)               :: lr(fbem_n_nodes(type_f2),0:2,0:2)   ! l integration kernels matrix (regular part)
    complex(kind=real64)               :: z1, z2, z3
    complex(kind=real64)               :: K0r_z1, K1r_z1, K2r_z1
    complex(kind=real64)               :: K0r_z2, K1r_z2, K2r_z2
    complex(kind=real64)               :: K0r_z3, K1r_z3, K2r_z3
    complex(kind=real64)               :: W0r, T01r, T02r, W1r, W2r
    complex(kind=real64)               :: T1r, T2r, T3r
    complex(kind=real64)               :: Q01r, Q02r, S01r, S02r, S03r
    complex(kind=real64)               :: S1r, S2r, S3r, S4r, S5r
    complex(kind=real64)               :: fs_d(0:2,0:2), fs_s(0:2,0:2)
    complex(kind=real64)               :: cte_d(0:2,0:2), cte_s(0:2,0:2)
    !
    ! Initialize
    !
    ! The kernel matrices are set to zero
    m =(0.d0,0.d0)
    l =(0.d0,0.d0)
    mr=(0.d0,0.d0)
    lr=(0.d0,0.d0)
    ! Number of geometrical nodes
    nnodes_g=fbem_n_nodes(type_g)
    ! The HBIE can not be collocated at vertices
    if (fbem_check_xi_vertex(xi_i).eqv.(.true.)) then
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                              'the HBIE cannot be collocated at a vertex')
    ! If xi_i is inside the element, 2 subdivisions are needed
    else
      nsub=2
      ! Coordinates xi of the subdivision 1
      xisub(1,1)=-1.d0
      xisub(2,1)=xi_i
      ! Coordinate xip of the collocation point
      xip_i(1)=1.d0
      ! Coordinates xi of the subdivision 2
      xisub(1,2)=xi_i
      xisub(2,2)=1.d0
      ! Coordinate xip of the collocation point
      xip_i(2)=0.d0
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
    x_i=0.d0
    T_i=0.d0
    do kphi=1,nnodes_G
      x_i=x_i+phi_g(kphi)*x_nodes(:,kphi)
      T_i=T_i+dphidxi_g(kphi)*x_nodes(:,kphi)
    end do
    ! Normal vector
    N_i(1)=T_i(2)
    N_i(2)=-T_i(1)
    ! Geometric jacobian
    jgi=dsqrt(T_i(1)**2+T_i(2)**2)
    ! Unit normal vector at collocation point
    n_i=N_i/jgi
    ! Unit tangent vector at collocation point
    t_i=T_i/jgi
    !
    ! Functional shape functions and its first derivatives (primary variables) at xi_i
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
    ! Functional shape functions and its first derivatives (secondary variables) at xi_i
    !
#   define etype type_f2
#   define delta delta_f
#   define xi xi_i
#   define phi phi_f2_i
#   define dphidxi dphidxi_f2_i
#   include <phi_and_dphidxi_1d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
#   undef dphidxi
    !
    ! Integrate WEAKLY SINGULAR integrals
    !
    ! Loop through SUBDIVISIONS
    do ksub=1,nsub
      ! Telles transformation parameters (barr=0 at the collocation point)
      telles_parameters=fbem_telles01_calculate_parameters(xip_i(ksub),0.d0)
      ! Jacobian of xip->xi transformation (is constant)
      js=xisub(2,ksub)-xisub(1,ksub)
      ! Loop through INTEGRATION POINTS
      do kip=1,gl01_n(ngp)
        ! GAMMA->XIP->XI COORDINATE TRANSFORMATION
        ! Coordinate and weight in gamma [0,1]
        gamma=gl01_xi(kip,ngp)
        w=gl01_w(kip,ngp)
        ! Calculate xip and gamma->xip jacobian
        call fbem_telles_xi_and_jacobian(telles_parameters,gamma,xip,jt)
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
        jg=sqrt(T(1)**2+T(2)**2)
        ! Unit normal vector
        n=N/jg
        ! Unit tangent vector
        t=T/jg
        ! Distance vector
        rv=x-x_i
        ! Distance vector norm and its inverse
        r=sqrt(rv(1)**2+rv(2)**2)
        logr=log(r)
        d1r=1.0d0/r
        d1r2=1.0d0/r**2
        ! r_{,k}
        drdx=rv*d1r
        ! dr/dn
        drdn=dot_product(drdx,n)
        ! dr/dn_i
        drdni=-dot_product(drdx,n_i)
        ! n · n_i
        n_dot_ni=dot_product(n,n_i)
        ! dr/dGamma
        drdt=dot_product(drdx,t)
        ! FUNDAMENTAL SOLUTIONS
        ! Fluid load / Fluid response
        fs_d(0,0)=0.d0
        fs_s(0,0)=p%Q02(2)*n_dot_ni*logr
        ! Fluid load / Solid response
        fs_d(0,1:2)=p%W2(1)*n_i*logr
        fs_s(0,1:2)=0.d0
        ! Solid load / Fluid response
        fs_d(1:2,0)=p%T02(1)*n_i*logr
        fs_s(1:2,0)=0.d0
        ! Solid load / Solid response
        do il=1,2
          do ik=1,2
            fs_d(il,ik)=0.d0
            fs_s(il,ik)=(p%S4(2)*(n_dot_ni*c_dkr(il,ik)+n(il)*n_i(ik))+p%S5(2)*n(ik)*n_i(il))*logr
          end do
        end do
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
        ! FUNCTIONAL SHAPE FUNCTIONS * JACOBIAN * WEIGHT
        jw=jg*js*jt*w
        phif1jw=phi_f1*jw
        phif2jw=phi_f2*jw
        ! ADD INTEGRAND EVALUATION
        do il=0,2
          do ik=0,2
            m(:,il,ik)=m(:,il,ik)+fs_s(il,ik)*phif1jw
            l(:,il,ik)=l(:,il,ik)+fs_d(il,ik)*phif2jw
          end do
        end do
      end do ! Loop through INTEGRATION POINTS
    end do ! Loop through SUBDIVISIONS
    !
    ! Integrate REGULAR integrals of the SINGULAR parts
    !
    ! Loop through SUBDIVISIONS
    do ksub=1,nsub
      ! Jacobian of xip->xi transformation (is constant)
      js=xisub(2,ksub)-xisub(1,ksub)
      ! Loop through INTEGRATION POINTS
      do kip=1,gl01_n(ngp)
        ! XIP->XI COORDINATE TRANSFORMATION
        ! Coordinate and weight in xip [0,1]
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
        jg=sqrt(T(1)**2+T(2)**2)
        ! Unit normal vector
        n=N/jg
        ! Unit tangent vector
        t=T/jg
        ! Distance vector
        rv=x-x_i
        ! Distance vector norm and its inverse
        r=sqrt(rv(1)**2+rv(2)**2)
        logr=log(r)
        d1r=1.0d0/r
        d1r2=1.0d0/r**2
        ! r_{,k}
        drdx=rv*d1r
        ! dr/dn
        drdn=dot_product(drdx,n)
        ! dr/dn_i
        drdni=-dot_product(drdx,n_i)
        ! n · n_i
        n_dot_ni=dot_product(n,n_i)
        ! dr/dGamma
        drdt=dot_product(drdx,t)
        ! FUNDAMENTAL SOLUTIONS
        ! Fluid load / Fluid response
        fs_d(0,0)=p%W0(1)*d1r*drdni
        fs_s(0,0)=p%Q01(1)*d1r2*drdn*drdni
        ! Fluid load / Solid response
        fs_d(0,1:2)=0.d0
        fs_s(0,1:2)=p%S01(1)*d1r*drdx*drdn*drdni+p%S03(1)*d1r*n_i*drdn
        ! Solid load / Fluid response
        fs_d(1:2,0)=0.d0
        fs_s(1:2,0)=p%Z*d1r*n_i*drdn-p%S01(1)*d1r*drdx*drdn*drdni+p%S03(1)*d1r*n*drdni&
                   +p%S03(1)*(n(2)*c_dkr(1,il)-n(1)*c_dkr(2,1:2))*d1r*(drdx(2)*n_i(1)-drdx(1)*n_i(2)-drdt)&
                   +p%S03(1)*((n(2)-n_i(2))*c_dkr(1,1:2)-(n(1)-n_i(1))*c_dkr(2,1:2))*d1r*drdt
        ! Solid load / Solid response
        do il=1,2
          do ik=1,2
            fs_d(il,ik)=p%T2(1)*d1r*drdx(il)*drdx(ik)*drdni+p%T3(1)*d1r*drdni*c_dkr(il,ik)
            fs_s(il,ik)=p%S3(1)*d1r2*drdx(il)*drdx(ik)*drdn*drdni&
                       +(1.0d0-c_dkr(il,ik))*0.125d0*p%S3(1)*d1r2*(2.0d0*drdx(1)*drdx(2)*n_dot_ni+n(1)*n_i(2)+n(2)*n_i(1))&
                       +c_dkr(il,ik)*0.25d0*p%S3(1)*d1r2*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni)
          end do
        end do
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
        ! FUNCTIONAL SHAPE FUNCTIONS * JACOBIAN * WEIGHT
        jw=jg*js*w
        phif1jw=phi_f1*jw
        phif2jw=phi_f2*jw
        ! ADD INTEGRAND EVALUATION
        do il=0,2
          do ik=0,2
            m(:,il,ik)=m(:,il,ik)+fs_s(il,ik)*phif1jw
            l(:,il,ik)=l(:,il,ik)+fs_d(il,ik)*phif2jw
          end do
        end do
      end do ! Loop through INTEGRATION POINTS
    end do ! Loop through SUBDIVISIONS
    !
    ! Integrate REGULAR integrals of the CPV and HFP integrals
    !
    ! Loop through SUBDIVISIONS
    do ksub=1,nsub
      ! Jacobian of xip->xi transformation (is constant)
      js=xisub(2,ksub)-xisub(1,ksub)
      ! Loop through INTEGRATION POINTS
      do kip=1,gl01_n(ngp)
        ! XIP->XI COORDINATE TRANSFORMATION
        ! Coordinate and weight in xip [0,1]
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
        jg=sqrt(T(1)**2+T(2)**2)
        ! Unit normal vector
        n=N/jg
        ! Unit tangent vector
        t=T/jg
        ! Distance vector
        rv=x-x_i
        ! Distance vector norm and its inverse
        r=sqrt(rv(1)**2+rv(2)**2)
        logr=log(r)
        d1r=1.0d0/r
        d1r2=1.0d0/r**2
        ! r_{,k}
        drdx=rv*d1r
        ! dr/dn
        drdn=dot_product(drdx,n)
        ! dr/dn_i
        drdni=-dot_product(drdx,n_i)
        ! n · n_i
        n_dot_ni=dot_product(n,n_i)
        ! dr/dGamma
        drdt=dot_product(drdx,t)
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
        ! FUNCTIONAL SHAPE FUNCTIONS * JACOBIAN * WEIGHT
        jw=jg*js*w
        phif1jw=phi_f1*jw
        phif2jw=phi_f2*jw
        ! Regular parts of the HFP
        m(:,0,0)=m(:,0,0)+p%Q02(1)*d1r2*(n_dot_ni-abs(drdt))*phif1jw
        if (ksub.eq.1) then
          m(:,0,0)=m(:,0,0)+p%Q02(1)*d1r2*abs(drdt)*(phi_f1-phi_f1_i+dphidxi_f1_i/jgi*r)*jw
        else
          m(:,0,0)=m(:,0,0)+p%Q02(1)*d1r2*abs(drdt)*(phi_f1-phi_f1_i-dphidxi_f1_i/jgi*r)*jw
        end if
        ! Regular parts of the CPV
        do ik=1,2
          m(:,0,ik)=m(:,0,ik)+p%S03(1)*(n_i(1)*c_dkr(2,ik)-n_i(2)*c_dkr(1,ik))*d1r*drdt*(phi_f1-phi_f1_i)*jw
        end do
        ! Regular parts of the CPV
        do il=1,2
          m(:,il,0)=m(:,il,0)+p%S03(1)*(n_i(2)*c_dkr(1,il)-n_i(1)*c_dkr(2,il))*d1r*drdt*(phi_f1-phi_f1_i)*jw
        end do
        ! Regular parts of the CPV
        l(:,1,2)=l(:,1,2)-p%T1(1)*d1r*(drdx(2)*n_i(1)-drdx(1)*n_i(2)-drdt)*phif2jw
        l(:,2,1)=l(:,2,1)+p%T1(1)*d1r*(drdx(2)*n_i(1)-drdx(1)*n_i(2)-drdt)*phif2jw
        ! Other Regular parts of the CPV
        l(:,1,2)=l(:,1,2)-p%T1(1)*d1r*drdt*(phi_f2-phi_f2_i)*jw
        l(:,2,1)=l(:,2,1)+p%T1(1)*d1r*drdt*(phi_f2-phi_f2_i)*jw
        ! Regular parts of the HFP
        do il=1,2
          m(:,il,il)=m(:,il,il)+0.125d0*p%S3(1)*d1r2*(n_dot_ni-abs(drdt))*phif1jw
        end do
        do il=1,2
          if (ksub.eq.1) then
            m(:,il,il)=m(:,il,il)+0.125d0*p%S3(1)*d1r2*abs(drdt)*(phi_f1-phi_f1_i+dphidxi_f1_i/jgi*r)*jw
          else
            m(:,il,il)=m(:,il,il)+0.125d0*p%S3(1)*d1r2*abs(drdt)*(phi_f1-phi_f1_i-dphidxi_f1_i/jgi*r)*jw
          end if
        end do
      end do ! Loop through INTEGRATION POINTS
    end do ! Loop through SUBDIVISIONS
    !
    ! Add analytical terms of CPV and HFP integrals
    !
    ! Calculate ra and rb
    ra=sqrt((x_nodes(1,1)-x_i(1))**2+(x_nodes(2,1)-x_i(2))**2)
    rb=sqrt((x_nodes(1,2)-x_i(1))**2+(x_nodes(2,2)-x_i(2))**2)
    ! Fluid load / Fluid response
    ! HFP
    m(:,0,0)=m(:,0,0)+p%Q02(1)*(-phi_f1_i*(1.0d0/ra+1.0d0/rb)+dphidxi_f1_i/jgi*(log(rb)-log(ra)))
    ! Fluid load / Solid response
    ! CPV
    do ik=1,2
      m(:,0,ik)=m(:,0,ik)+p%S03(1)*(n_i(1)*c_dkr(2,ik)-n_i(2)*c_dkr(1,ik))*(log(rb)-log(ra))*phi_f1_i
    end do
    ! Solid load / Fluid response
    ! Regular parts of the CPV
    do il=1,2
      m(:,il,0)=m(:,il,0)+p%S03(1)*(n_i(2)*c_dkr(1,il)-n_i(1)*c_dkr(2,il))*(log(rb)-log(ra))*phi_f1_i
    end do
    ! Solid load / Solid response
    ! CPV
    l(:,1,2)=l(:,1,2)-p%T1(1)*(log(rb)-log(ra))*phi_f2_i
    l(:,2,1)=l(:,2,1)+p%T1(1)*(log(rb)-log(ra))*phi_f2_i
    ! HFP
    do il=1,2
      m(:,il,il)=m(:,il,il)+0.125d0*p%S3(1)*(-phi_f1_i*(1.0d0/ra+1.0d0/rb)+dphidxi_f1_i/jgi*(dlog(rb)-dlog(ra)))
    end do
    !
    ! Integrate REGULAR integrals of the REGULAR parts
    !
    ! Loop through subdivisions
    do ksub=1,nsub
      ! Jacobian of xip->xi transformation (is constant)
      js=xisub(2,ksub)-xisub(1,ksub)
      ! Loop through integration points
      do kip=1,gl01_n(ngp)
        ! XIP->XI COORDINATE TRANSFORMATION
        ! Coordinate and weight in xip [0,1]
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
        jg=sqrt(T(1)**2+T(2)**2)
        ! Unit normal vector
        n=N/jg
        ! Unit tangent vector
        t=T/jg
        ! Distance vector
        rv=x-x_i
        ! Distance vector norm and its inverse
        r=sqrt(rv(1)**2+rv(2)**2)
        logr=log(r)
        d1r=1.0d0/r
        d1r2=1.0d0/r**2
        ! r_{,k}
        drdx=rv*d1r
        ! dr/dn
        drdn=dot_product(drdx,n)
        ! dr/dn_i
        drdni=-dot_product(drdx,n_i)
        ! n · n_i
        n_dot_ni=dot_product(n,n_i)
        ! dr/dGamma
        drdt=dot_product(drdx,t)
        ! CALCULATE THE COMPONENTS OF THE REGULAR PARTS OF THE FUNDAMENTAL SOLUTIONS
        ! z_m = i k_m r
        z1=c_im*p%k1*r
        z2=c_im*p%k2*r
        z3=c_im*p%k3*r
        ! Calculate K_n^R for z1, z2 and z3
        call fbem_modified_bessel_K0r_K1r_K2r_2(z1,K0r_z1,K1r_z1,K2r_z1)
        call fbem_modified_bessel_K0r_K1r_K2r_2(z2,K0r_z2,K1r_z2,K2r_z2)
        call fbem_modified_bessel_K0r_K1r_K2r_2(z3,K0r_z3,K1r_z3,K2r_z3)
        ! Calculate components (regular parts)
        W0r=(p%W0(2)*logr+p%W0(3))*r+p%W0(4)*K1r_z1+p%W0(5)*K1r_z2
        T01r=p%T01(1)+p%T01(2)*K0r_z1+p%T01(3)*K0r_z2+(p%T01(4)*K1r_z1+p%T01(5)*K1r_z2)*d1r
        T02r=p%T02(2)+p%T02(3)*K0r_z1+p%T02(4)*K0r_z2+(p%T02(5)*K1r_z1+p%T02(6)*K1r_z2)*d1r
        W1r=p%W1(1)+(p%W1(2)*logr+p%W1(3))*r**2+p%W1(4)*K0r_z1+p%W1(5)*K0r_z2+(p%W1(6)*K1r_z1+p%W1(7)*K1r_z2)*d1r&
             +p%W1(8)*K2r_z3+p%W1(9)*K2r_z1+p%W1(10)*K2r_z2
        W2r=p%W2(2)+p%W2(3)*K0r_z3+(p%W2(4)*K1r_z3+p%W2(5)*K1r_z1+p%W2(6)*K1r_z2)*d1r
        T1r=(p%T1(2)*logr+p%T1(3))*r+p%T1(4)*K1r_z1+p%T1(5)*K1r_z2+(p%T1(6)*K2r_z3+p%T1(7)*K2r_z1+p%T1(8)*K2r_z2)*d1r
        T2r=p%T2(2)*r+p%T2(3)*K1r_z3+p%T2(4)*K1r_z1+p%T2(5)*K1r_z2+(p%T2(6)*K2r_z3+p%T2(7)*K2r_z1+p%T2(8)*K2r_z2)*d1r
        T3r=(p%T3(2)*logr+p%T3(3))*r+p%T3(4)*K1r_z3+(p%T3(5)*K2r_z3+p%T3(6)*K2r_z1+p%T3(7)*K2r_z2)*d1r
        Q01r=p%Q01(2)+(p%Q01(3)*logr+p%Q01(4))*r**2+p%Q01(5)*K0r_z1+p%Q01(6)*K0r_z2&
              +(p%Q01(7)*K1r_z1+p%Q01(8)*K1r_z2)*d1r+p%Q01(9)*K2r_z3+p%Q01(10)*K2r_z1+p%Q01(11)*K2r_z2
        Q02r=p%Q02(3)+p%Q02(4)*K0r_z3+(p%Q02(5)*K1r_z3+p%Q02(6)*K1r_z1+p%Q02(7)*K1r_z2)*d1r
        S01r=(p%S01(2)*logr+p%S01(3))*r+(p%S01(4)*K0r_z1+p%S01(5)*K0r_z2)*d1r+(p%S01(6)*K1r_z1+p%S01(7)*K1r_z2)*d1r2&
              +p%S01(8)*K1r_z3+p%S01(9)*K1r_z1+p%S01(10)*K1r_z2+(p%S01(11)*K2r_z3+p%S01(12)*K2r_z1+p%S01(13)*K2r_z2)*d1r
        S02r=(p%S02(2)*logr+p%S02(3))*r+(p%S02(4)*K0r_z1+p%S02(5)*K0r_z2)*d1r+(p%S02(6)*K1r_z1+p%S02(7)*K1r_z2)*d1r2&
              +p%S02(8)*K1r_z1+p%S02(9)*K1r_z2+(p%S02(10)*K2r_z3+p%S02(11)*K2r_z1+p%S02(12)*K2r_z2)*d1r
        S03r=(p%S03(2)*logr+p%S03(3))*r+(p%S03(4)*K0r_z1+p%S03(5)*K0r_z2)*d1r+(p%S03(6)*K1r_z1+p%S03(7)*K1r_z2)*d1r2&
              +p%S03(8)*K1r_z3+(p%S03(9)*K2r_z3+p%S03(10)*K2r_z1+p%S03(11)*K2r_z2)*d1r
        S1r=p%S1(2)+p%S1(3)*K0r_z3+(p%S1(4)*K1r_z3+p%S1(5)*K1r_z1+p%S1(6)*K1r_z2)*d1r&
             +(p%S1(7)*K2r_z3+p%S1(8)*K2r_z1+p%S1(9)*K2r_z2)*d1r2
        S2r=p%S2(2)+p%S2(3)*K0r_z1+p%S2(4)*K0r_z2+(p%S2(5)*K1r_z3+p%S2(6)*K1r_z1+p%S2(7)*K1r_z2)*d1r&
             +(p%S2(8)*K2r_z3+p%S2(9)*K2r_z1+p%S2(10)*K2r_z2)*d1r2
        S3r=p%S3(2)+p%S3(3)*K0r_z3+p%S3(4)*K0r_z1+p%S3(5)*K0r_z2+(p%S3(6)*K1r_z3+p%S3(7)*K1r_z1+p%S3(8)*K1r_z2)*d1r&
             +(p%S3(9)*K2r_z3+p%S3(10)*K2r_z1+p%S3(11)*K2r_z2)*d1r2
        S4r=p%S4(3)+p%S4(4)*K1r_z3*d1r+(p%S4(5)*K2r_z3+p%S4(6)*K2r_z1+p%S4(7)*K2r_z2)*d1r2
        S5r=p%S5(3)+p%S5(4)*K0r_z1+p%S5(5)*K0r_z2+(p%S5(6)*K1r_z1+p%S5(7)*K1r_z2)*d1r&
             +(p%S5(8)*K2r_z3+p%S5(9)*K2r_z1+p%S5(10)*K2r_z2)*d1r2
        ! FUNDAMENTAL SOLUTIONS
        ! Fluid load / Fluid response
        fs_d(0,0)=W0r*drdni
        fs_s(0,0)=Q01r*drdn*drdni+Q02r*n_dot_ni
        ! Fluid load / Solid response
        fs_d(0,1:2)=-W1r*drdx*drdni+W2r*n_i
        fs_s(0,1:2)=S01r*drdx*drdn*drdni+S02r*n*drdni+S03r*(n_i*drdn+drdx*n_dot_ni)
        ! Solid load / Fluid response
        fs_d(1:2,0)=-T01r*drdx*drdni+T02r*n_i
        fs_s(1:2,0)=-S01r*drdx*drdn*drdni+S02r*n_i*drdn-S03r*(-n*drdni+drdx*n_dot_ni)
        ! Solid load / Solid response
        do il=1,2
          do ik=1,2
            fs_d(il,ik)=-T1r*drdx(ik)*n_i(il)+T2r*drdx(il)*drdx(ik)*drdni-T3r*(drdx(il)*n_i(ik)-drdni*c_dkr(il,ik))
            fs_s(il,ik)=S1r*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni-drdn*drdni*c_dkr(il,ik)+drdx(il)*drdx(ik)*n_dot_ni)&
                       +S2r*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni)+S3r*drdx(il)*drdx(ik)*drdn*drdni&
                       +S4r*(n_dot_ni*c_dkr(il,ik)+n(il)*n_i(ik))+S5r*n(ik)*n_i(il)
          end do
        end do
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
        ! FUNCTIONAL SHAPE FUNCTIONS * JACOBIAN * WEIGHT
        jw=jg*js*w
        phif1jw=phi_f1*jw
        phif2jw=phi_f2*jw
        ! ADD INTEGRAND EVALUATION
        do il=0,2
          do ik=0,2
            mr(:,il,ik)=mr(:,il,ik)+fs_s(il,ik)*phif1jw
            lr(:,il,ik)=lr(:,il,ik)+fs_d(il,ik)*phif2jw
          end do
        end do
      end do ! Loop through integration points
    end do ! Loop through subdivisions
    !
    ! Add SINGULAR and REGULAR parts
    !
    m=m+mr
    l=l+lr
    !
    ! Calculate cte_s and cte_d, and multiply them to m and l, respectively
    !
    ! Fluid load / Fluid response constants
    cte_d(0,0)=-c_1_2pi/p%J
    cte_s(0,0)=-c_1_2pi
    ! Fluid load / Solid response constants
    cte_d(0,1:2)=c_1_2pi/p%mu
    cte_s(0,1:2)=c_1_2pi
    ! Solid load / Fluid response constants
    cte_d(1:2,0)=-c_1_2pi/p%J
    cte_s(1:2,0)=-c_1_2pi
    ! Solid load / Solid response constants
    cte_d(1:2,1:2)=c_1_2pi
    cte_s(1:2,1:2)=c_1_2pi
    ! Multiply by constants
    do il=0,2
      do ik=0,2
        m(:,il,ik)=cte_s(il,ik)*m(:,il,ik)
        l(:,il,ik)=cte_d(il,ik)*l(:,il,ik)
      end do
    end do
    !
    ! If the element orientation is reversed
    !
    if (reverse) l=-l
  end subroutine fbem_bem_harpor2d_hbie_int

  subroutine fbem_bem_harpor2d_hbie_auto(e,reverse,x_i,n_i,p,qsp,ns,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                     !! Integration element
    logical                            :: reverse               !! Reverse orientation
    real(kind=real64)                  :: x_i(2)                !! Collocation point
    real(kind=real64)                  :: n_i(2)                !! Unit normal at the collocation point
    type(fbem_bem_harpor2d_parameters) :: p                     !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp                   !! Quasi-singular integration parameters
    integer                            :: ns                    !! Maximum level of subdivisions
    complex(kind=real64)               :: m(e%n_pnodes,0:2,0:2) !! m integration kernel
    complex(kind=real64)               :: l(e%n_snodes,0:2,0:2) !! l integration kernel
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
        call fbem_bem_harpor2d_hbie_int(30,e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi(1),p,m,l)
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
          call fbem_bem_harpor2d_hbie_ext_pre(ps,e,reverse,x_i,n_i,p,m,l)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_harpor2d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,p,qsp,1,ns,m,l)
        end if
    end select
  end subroutine fbem_bem_harpor2d_hbie_auto

  ! ================================================================================================================================
  ! BOTH SBIE AND HBIE SIMULTANEOUSLY

  subroutine fbem_bem_harpor2d_shbie_ext_pre(ps,e,reverse,x_i,n_i,p,h,g,m,l)
    implicit none
    ! I/O
    integer                            :: ps                    !! Selected precalculated dataset
    type(fbem_bem_element)             :: e                     !! Element
    logical                            :: reverse               !! Reverse normal vector
    real(kind=real64)                  :: x_i(2)                !! Collocation point position vector
    real(kind=real64)                  :: n_i(2)                !! Collocation point unit normal vector
    type(fbem_bem_harpor2d_parameters) :: p                     !! Parameters of the region
    complex(kind=real64)               :: h(e%n_pnodes,0:2,0:2) !! h kernels matrix
    complex(kind=real64)               :: g(e%n_snodes,0:2,0:2) !! g kernels matrix
    complex(kind=real64)               :: m(e%n_pnodes,0:2,0:2) !! m kernels matrix
    complex(kind=real64)               :: l(e%n_snodes,0:2,0:2) !! l kernels matrix
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
    complex(kind=real64) :: z(3)                           ! Arguments z=ikr
    complex(kind=real64) :: KnR(0:2,3)                     ! Bessel functions decomposition
    complex(kind=real64) :: eta, vartheta, psi, chi        ! Fundamental solution components
    complex(kind=real64) :: W0, T01, T02, W1, W2           ! Fundamental solution components
    complex(kind=real64) :: T1, T2, T3                     ! Fundamental solution components
    complex(kind=real64) :: Q01, Q02, S01, S02, S03        ! Fundamental solution components
    complex(kind=real64) :: S1, S2, S3, S4, S5             ! Fundamental solution components
    complex(kind=real64) :: fs_u(0:2,0:2), fs_t(0:2,0:2)   ! Fundamental solutions
    complex(kind=real64) :: fs_d(0:2,0:2), fs_s(0:2,0:2)   ! Fundamental solutions
    complex(kind=real64) :: cte_u(0:2,0:2), cte_t(0:2,0:2) ! Fundamental solution constants
    complex(kind=real64) :: cte_d(0:2,0:2), cte_s(0:2,0:2) ! Fundamental solution constants
    ! Initialisation
    h=(0.0d0,0.0d0)
    g=(0.0d0,0.0d0)
    m=(0.d0,0.d0)
    l=(0.d0,0.d0)
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
      z(3)=c_im*p%k3*r
      call fbem_BesselKnR_decomposed(3,z,KnR)
      eta=p%eta(1)*logr+p%eta(2)+p%eta(3)*KnR(0,1)+p%eta(4)*KnR(0,2)
      vartheta=(p%vartheta(1)*logr+p%vartheta(2))*r+p%vartheta(3)*KnR(1,1)+p%vartheta(4)*KnR(1,2)
      psi=p%psi(1)*logr+p%psi(2)+KnR(0,3)+(p%psi(3)*KnR(1,3)+p%psi(4)*KnR(1,1)+p%psi(5)*KnR(1,2))*d1r1
      chi=p%chi(1)+(p%chi(2)*logr+p%chi(3))*r2+KnR(2,3)+p%chi(4)*KnR(2,1)+p%chi(5)*KnR(2,2)
      W0=p%W0(1)*d1r1+(p%W0(2)*logr+p%W0(3))*r+p%W0(4)*KnR(1,1)+p%W0(5)*KnR(1,2)
      T01=p%T01(1)+p%T01(2)*KnR(0,1)+p%T01(3)*KnR(0,2)+(p%T01(4)*KnR(1,1)+p%T01(5)*KnR(1,2))*d1r1
      T02=p%T02(1)*logr+p%T02(2)+p%T02(3)*KnR(0,1)+p%T02(4)*KnR(0,2)+(p%T02(5)*KnR(1,1)+p%T02(6)*KnR(1,2))*d1r1
      W1=p%W1(1)+(p%W1(2)*logr+p%W1(3))*r2+p%W1(4)*KnR(0,1)+p%W1(5)*KnR(0,2)+(p%W1(6)*KnR(1,1)+p%W1(7)*KnR(1,2))*d1r1&
        +p%W1(8)*KnR(2,3)+p%W1(9)*KnR(2,1)+p%W1(10)*KnR(2,2)
      W2=p%W2(1)*logr+p%W2(2)+p%W2(3)*KnR(0,3)+(p%W2(4)*KnR(1,3)+p%W2(5)*KnR(1,1)+p%W2(6)*KnR(1,2))*d1r1
      T1=p%T1(1)*d1r1+(p%T1(2)*logr+p%T1(3))*r+p%T1(4)*KnR(1,1)+p%T1(5)*KnR(1,2)&
        +(p%T1(6)*KnR(2,3)+p%T1(7)*KnR(2,1)+p%T1(8)*KnR(2,2))*d1r1
      T2=p%T2(1)*d1r1+p%T2(2)*r+p%T2(3)*KnR(1,3)+p%T2(4)*KnR(1,1)+p%T2(5)*KnR(1,2)&
        +(p%T2(6)*KnR(2,3)+p%T2(7)*KnR(2,1)+p%T2(8)*KnR(2,2))*d1r1
      T3=p%T3(1)*d1r1+(p%T3(2)*logr+p%T3(3))*r+p%T3(4)*KnR(1,3)+(p%T3(5)*KnR(2,3)+p%T3(6)*KnR(2,1)+p%T3(7)*KnR(2,2))*d1r1
      Q01=p%Q01(1)*d1r2+p%Q01(2)+(p%Q01(3)*logr+p%Q01(4))*r2+p%Q01(5)*KnR(0,1)+p%Q01(6)*KnR(0,2)&
          +(p%Q01(7)*KnR(1,1)+p%Q01(8)*KnR(1,2))*d1r1+p%Q01(9)*KnR(2,3)+p%Q01(10)*KnR(2,1)+p%Q01(11)*KnR(2,2)
      Q02=p%Q02(1)*d1r2+p%Q02(2)*logr+p%Q02(3)+p%Q02(4)*KnR(0,3)+(p%Q02(5)*KnR(1,3)+p%Q02(6)*KnR(1,1)+p%Q02(7)*KnR(1,2))*d1r1
      S01=p%S01(1)*d1r1+(p%S01(2)*logr+p%S01(3))*r+(p%S01(4)*KnR(0,1)+p%S01(5)*KnR(0,2))*d1r1+(p%S01(6)*KnR(1,1)+p%S01(7)*KnR(1,2))*d1r2&
          +p%S01(8)*KnR(1,3)+p%S01(9)*KnR(1,1)+p%S01(10)*KnR(1,2)+(p%S01(11)*KnR(2,3)+p%S01(12)*KnR(2,1)+p%S01(13)*KnR(2,2))*d1r1
      S02=p%S02(1)*d1r1+(p%S02(2)*logr+p%S02(3))*r+(p%S02(4)*KnR(0,1)+p%S02(5)*KnR(0,2))*d1r1+(p%S02(6)*KnR(1,1)+p%S02(7)*KnR(1,2))*d1r2&
          +p%S02(8)*KnR(1,1)+p%S02(9)*KnR(1,2)+(p%S02(10)*KnR(2,3)+p%S02(11)*KnR(2,1)+p%S02(12)*KnR(2,2))*d1r1
      S03=p%S03(1)*d1r1+(p%S03(2)*logr+p%S03(3))*r+(p%S03(4)*KnR(0,1)+p%S03(5)*KnR(0,2))*d1r1+(p%S03(6)*KnR(1,1)+p%S03(7)*KnR(1,2))*d1r2&
          +p%S03(8)*KnR(1,3)+(p%S03(9)*KnR(2,3)+p%S03(10)*KnR(2,1)+p%S03(11)*KnR(2,2))*d1r1
      S1=p%S1(1)*d1r2+p%S1(2)+p%S1(3)*KnR(0,3)+(p%S1(4)*KnR(1,3)+p%S1(5)*KnR(1,1)+p%S1(6)*KnR(1,2))*d1r1&
         +(p%S1(7)*KnR(2,3)+p%S1(8)*KnR(2,1)+p%S1(9)*KnR(2,2))*d1r2
      S2=p%S2(1)*d1r2+p%S2(2)+p%S2(3)*KnR(0,1)+p%S2(4)*KnR(0,2)+(p%S2(5)*KnR(1,3)+p%S2(6)*KnR(1,1)+p%S2(7)*KnR(1,2))*d1r1&
         +(p%S2(8)*KnR(2,3)+p%S2(9)*KnR(2,1)+p%S2(10)*KnR(2,2))*d1r2
      S3=p%S3(1)*d1r2+p%S3(2)+p%S3(3)*KnR(0,3)+p%S3(4)*KnR(0,1)+p%S3(5)*KnR(0,2)+(p%S3(6)*KnR(1,3)+p%S3(7)*KnR(1,1)+p%S3(8)*KnR(1,2))*d1r1&
         +(p%S3(9)*KnR(2,3)+p%S3(10)*KnR(2,1)+p%S3(11)*KnR(2,2))*d1r2
      S4=p%S4(1)*d1r2+p%S4(2)*logr+p%S4(3)+p%S4(4)*KnR(1,3)*d1r1+(p%S4(5)*KnR(2,3)+p%S4(6)*KnR(2,1)+p%S4(7)*KnR(2,2))*d1r2
      S5=p%S5(1)*d1r2+p%S5(2)*logr+p%S5(3)+p%S5(4)*KnR(0,1)+p%S5(5)*KnR(0,2)+(p%S5(6)*KnR(1,1)+p%S5(7)*KnR(1,2))*d1r1&
         +(p%S5(8)*KnR(2,3)+p%S5(9)*KnR(2,1)+p%S5(10)*KnR(2,2))*d1r2
      ! Fundamental solutions
      fs_u(0,0)=eta
      fs_t(0,0)=W0*drdn
      fs_d(0,0)=W0*drdni
      fs_s(0,0)=Q01*drdn*drdni+Q02*n_dot_ni
      fs_u(0,1:2)=vartheta*drdx
      fs_t(0,1:2)=T01*drdx*drdn+T02*n
      fs_d(0,1:2)=-W1*drdx*drdni+W2*n_i
      fs_s(0,1:2)=S01*drdx*drdn*drdni+S02*n*drdni+S03*(n_i*drdn+drdx*n_dot_ni)
      fs_u(1:2,0)=vartheta*drdx
      fs_t(1:2,0)=W1*drdx*drdn+W2*n
      fs_d(1:2,0)=-T01*drdx*drdni+T02*n_i
      fs_s(1:2,0)=-S01*drdx*drdn*drdni+S02*n_i*drdn-S03*(-n*drdni+drdx*n_dot_ni)
      do il=1,2
        do ik=1,2
          fs_u(il,ik)=psi*c_dkr(il,ik)-chi*drdx(il)*drdx(ik)
          fs_t(il,ik)=T1*drdx(il)*n(ik)+T2*drdx(il)*drdx(ik)*drdn+T3*(drdx(ik)*n(il)+drdn*c_dkr(il,ik))
          fs_d(il,ik)=-T1*drdx(ik)*n_i(il)+T2*drdx(il)*drdx(ik)*drdni-T3*(drdx(il)*n_i(ik)-drdni*c_dkr(il,ik))
          fs_s(il,ik)=S1*(drdx(il)*n_i(ik)*drdn-drdx(ik)*n(il)*drdni-drdn*drdni*c_dkr(il,ik)+drdx(il)*drdx(ik)*n_dot_ni)&
                     +S2*(drdx(ik)*n_i(il)*drdn-drdx(il)*n(ik)*drdni)+S3*drdx(il)*drdx(ik)*drdn*drdni&
                     +S4*(n_dot_ni*c_dkr(il,ik)+n(il)*n_i(ik))+S5*n(ik)*n_i(il)
        end do
      end do
      ! Add the term of the numerical integration summation (fundamental solutions * functional shape functions * jacobian * weight)
      do il=0,2
        do ik=0,2
          h(:,il,ik)=h(:,il,ik)+fs_t(il,ik)*pphijw
          g(:,il,ik)=g(:,il,ik)+fs_u(il,ik)*sphijw
          m(:,il,ik)=m(:,il,ik)+fs_s(il,ik)*pphijw
          l(:,il,ik)=l(:,il,ik)+fs_d(il,ik)*sphijw
        end do
      end do
    end do
    ! Multiply by constants
    cte_u(0,0)=-c_1_2pi
    cte_t(0,0)=-c_1_2pi
    cte_d(0,0)=-c_1_2pi/p%J
    cte_s(0,0)=-c_1_2pi
    cte_u(0,1:2)=-c_1_2pi
    cte_t(0,1:2)= c_1_2pi
    cte_d(0,1:2)=c_1_2pi/p%mu
    cte_s(0,1:2)=c_1_2pi
    cte_u(1:2,0)=-c_1_2pi/p%J
    cte_t(1:2,0)=-c_1_2pi/p%mu
    cte_d(1:2,0)=-c_1_2pi/p%J
    cte_s(1:2,0)=-c_1_2pi
    cte_u(1:2,1:2)=c_1_2pi/p%mu
    cte_t(1:2,1:2)=c_1_2pi
    cte_d(1:2,1:2)=c_1_2pi
    cte_s(1:2,1:2)=c_1_2pi
    do il=0,2
      do ik=0,2
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
  end subroutine fbem_bem_harpor2d_shbie_ext_pre

  subroutine fbem_bem_harpor2d_shbie_auto(e,reverse,x_i,n_i,p,qsp,ns,h,g,m,l)
    implicit none
    ! I/O
    type(fbem_bem_element)             :: e                     !! Integration element
    logical                            :: reverse               !! Reverse orientation
    real(kind=real64)                  :: x_i(2)                !! Collocation point
    real(kind=real64)                  :: n_i(2)                !! Unit normal at the collocation point
    type(fbem_bem_harpor2d_parameters) :: p                     !! Parameters of the region
    type(fbem_qs_parameters)           :: qsp                   !! Quasi-singular integration parameters
    integer                            :: ns                    !! Maximum level of subdivisions
    complex(kind=real64)               :: h(e%n_pnodes,0:2,0:2) !! h integration kernel
    complex(kind=real64)               :: g(e%n_snodes,0:2,0:2) !! g integration kernel
    complex(kind=real64)               :: m(e%n_pnodes,0:2,0:2) !! m integration kernel
    complex(kind=real64)               :: l(e%n_snodes,0:2,0:2) !! l integration kernel
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
        call fbem_bem_harpor2d_sbie_int(30,e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi(1),p,h,g)
        call fbem_bem_harpor2d_hbie_int(30,e%gtype,e%ptype,e%stype,e%ptype_delta,e%x,reverse,barxi(1),p,m,l)
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
          call fbem_bem_harpor2d_shbie_ext_pre(ps,e,reverse,x_i,n_i,p,h,g,m,l)
        ! Integrate using an adaptative algorithm
        else
          call fbem_bem_harpor2d_sbie_ext_adp(e,reverse,xi_s,x_i,p,qsp,1,ns,h,g)
          call fbem_bem_harpor2d_hbie_ext_adp(e,reverse,xi_s,x_i,n_i,p,qsp,1,ns,m,l)
        end if
    end select
  end subroutine fbem_bem_harpor2d_shbie_auto

end module fbem_bem_harpor2d

