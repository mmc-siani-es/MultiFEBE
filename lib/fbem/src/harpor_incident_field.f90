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
!! <b> This module implements some incident fields for Biot's poroelastic media. </b>
module fbem_harpor_incident_field

  ! Fortran 2003 standard
  use iso_fortran_env
  use fbem_numerical
  use fbem_harela_incident_field

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all are private
  private

  ! Subroutines
  public :: fbem_harpor_vertical_plane_wave_layered_halfspace_setup
  public :: fbem_harpor_vertical_plane_wave_layered_halfspace
  public :: fbem_harpor_incident_plane_wave
  public :: harpor_calculate_permeable_rayleigh_wavenumber
  public :: fbem_harpor_incident_rayleigh_wave
  public :: fbem_phurkhao
  ! Parameters
  integer, parameter, public :: fbem_harpor_p1_wave         =1
  integer, parameter, public :: fbem_harpor_p2_wave         =2
  integer, parameter, public :: fbem_harpor_shx_wave        =3
  integer, parameter, public :: fbem_harpor_shy_wave        =4
  integer, parameter, public :: fbem_harpor_R_wave_permeable=5

contains

  subroutine fbem_harpor_vertical_plane_wave_layered_halfspace_setup(nl,top,properties,ic_k,omega,VA,k)
    implicit none
    ! I/O
    integer              :: nl               !! Number of layers
    real(kind=real64)    :: top(nl)          !! Top coordinate of each layer
    complex(kind=real64) :: properties(9,nl) !! Properties (phi,lambda,mu,Q,R,rhof,rhos,rhoa,b) of each layer
    complex(kind=real64) :: ic_k(nl)         !! Interface dissipation constant k
    real(kind=real64)    :: omega            !! Frequency w [rad/s]
    complex(kind=real64) :: VA(3,2,8,nl)     !! Amplitudes of each wave for each variable: u1,u2,U1,U2,tau11,tau12,tau22,tau
    complex(kind=real64) :: k(3,nl)          !! Wavenumbers for each layer
    ! Local
    real(kind=real64)    :: phi(nl)
    complex(kind=real64) :: lambda(nl)
    complex(kind=real64) :: mu(nl)
    complex(kind=real64) :: Q(nl)
    complex(kind=real64) :: R(nl)
    real(kind=real64)    :: rhof(nl)
    real(kind=real64)    :: rhos(nl)
    real(kind=real64)    :: rhoa(nl)
    real(kind=real64)    :: b(nl)
    complex(kind=real64) :: A(3,2,nl)
    integer              :: i
    complex(kind=real64) :: kp1, kp2, ks
    complex(kind=real64) :: hatrho11(nl), hatrho12(nl), hatrho22(nl)
    complex(kind=real64) :: phif(2,nl)
    real(kind=real64)    :: topr(nl)
    complex(kind=real64) :: TR, expmk1z, exppk1z, expmk2z, exppk2z, expmk1m1z, exppk1m1z, expmk2m1z, exppk2m1z
    complex(kind=real64) :: MA(4*nl,4*nl), VB(4*nl,1)
    !
    ! Initialization
    !
    A=(0.d0,0.d0)
    VA=(0.d0,0.d0)
    k=(0.d0,0.d0)
    do i=1,nl
      phi(i)=dreal(properties(1,i))
      lambda(i)=properties(2,i)
      mu(i)=properties(3,i)
      Q(i)=properties(4,i)
      R(i)=properties(5,i)
      rhof(i)=dreal(properties(6,i))
      rhos(i)=dreal(properties(7,i))
      rhoa(i)=dreal(properties(8,i))
      b(i)=dreal(properties(9,i))
      call harpor_calculate_bulk_wavenumbers(lambda(i),mu(i),Q(i),R(i),phi(i),rhos(i),rhof(i),rhoa(i),b(i),omega,kp1,kp2,ks)
      k(1,i)=kp1
      k(2,i)=kp2
      k(3,i)=ks
      hatrho11(i)=(1.d0-phi(i))*rhos(i)+rhoa(i)-c_im*b(i)/omega
      hatrho22(i)=       phi(i)*rhof(i)+rhoa(i)-c_im*b(i)/omega
      hatrho12(i)=                     -rhoa(i)+c_im*b(i)/omega
      phif(1,i)=(k(1,i)**2*(lambda(i)+2.d0*mu(i)+Q(i)**2/R(i))-omega**2*hatrho11(i))/(omega**2*hatrho12(i)-k(1,i)**2*Q(i))
      phif(2,i)=(k(2,i)**2*(lambda(i)+2.d0*mu(i)+Q(i)**2/R(i))-omega**2*hatrho11(i))/(omega**2*hatrho12(i)-k(2,i)**2*Q(i))
      topr(i)=top(i)-top(1)
    end do
    !
    ! P1-wave
    !
    ! Initialize
    MA=(0.d0,0.d0)
    VB=(0.d0,0.d0)
    !
    ! Boundary conditions at free-surface (z=x_3^(1))
    !
    ! Equation: u_3^(1)=1
    !
    ! Unknown: A_I1^(1)
    MA(1,1)=1.d0
    ! Unknown: A_R1^(1)
    MA(1,2)=1.d0
    ! Unknown: A_I2^(1)
    MA(1,3)=1.d0
    ! Unknown: A_R2^(1)
    MA(1,4)=1.d0
    ! Constant
    VB(1,1)=1.d0
    !
    ! Equation: tau_33^(1)=0
    !
    ! Unknown: A_I1^(1)
    MA(2,1)=-(lambda(1)+Q(1)**2/R(1)+2.d0*mu(1)+phif(1,1)*Q(1))*c_im*k(1,1)
    ! Unknown: A_R1^(1)
    MA(2,2)= (lambda(1)+Q(1)**2/R(1)+2.d0*mu(1)+phif(1,1)*Q(1))*c_im*k(1,1)
    ! Unknown: A_I2^(1)
    MA(2,3)=-(lambda(1)+Q(1)**2/R(1)+2.d0*mu(1)+phif(2,1)*Q(1))*c_im*k(2,1)
    ! Unknown: A_R2^(1)
    MA(2,4)= (lambda(1)+Q(1)**2/R(1)+2.d0*mu(1)+phif(2,1)*Q(1))*c_im*k(2,1)
    !
    ! Equation: tau^(1)=0
    !
    ! Unknown: A_I1^(1)
    MA(3,1)=-(Q(1)+phif(1,1)*R(1))*c_im*k(1,1)
    ! Unknown: A_R1^(1)
    MA(3,2)= (Q(1)+phif(1,1)*R(1))*c_im*k(1,1)
    ! Unknown: A_I2^(1)
    MA(3,3)=-(Q(1)+phif(2,1)*R(1))*c_im*k(2,1)
    ! Unknown: A_R2^(1)
    MA(3,4)= (Q(1)+phif(2,1)*R(1))*c_im*k(2,1)
    !
    ! Equation: A_I2^(nl)=0
    !
    ! Unknown: A_I2^(nl)
    MA(4,4*nl-1)=1.d0
    !
    ! Interface conditions between layer i and layer i-1 (z=x_3^(i))
    !
    do i=2,nl
      expmk1z  =exp(-c_im*k(1,i)*topr(i))
      exppk1z  =exp( c_im*k(1,i)*topr(i))
      expmk2z  =exp(-c_im*k(2,i)*topr(i))
      exppk2z  =exp( c_im*k(2,i)*topr(i))
      expmk1m1z=exp(-c_im*k(1,i-1)*topr(i))
      exppk1m1z=exp( c_im*k(1,i-1)*topr(i))
      expmk2m1z=exp(-c_im*k(2,i-1)*topr(i))
      exppk2m1z=exp( c_im*k(2,i-1)*topr(i))
      !
      ! Equation: u_3^(i) = u_3^(i-1)
      !
      ! Unknown: A_I1^(i)
      MA(4*i-3,4*i-3)= expmk1z
      ! Unknown: A_R1^(i)
      MA(4*i-3,4*i-2)= exppk1z
      ! Unknown: A_I2^(i)
      MA(4*i-3,4*i-1)= expmk2z
      ! Unknown: A_R2^(i)
      MA(4*i-3,4*i  )= exppk2z
      ! Unknown: A_I1^(i-1)
      MA(4*i-3,4*(i-1)-3)=-expmk1m1z
      ! Unknown: A_R1^(i-1)
      MA(4*i-3,4*(i-1)-2)=-exppk1m1z
      ! Unknown: A_I2^(i-1)
      MA(4*i-3,4*(i-1)-1)=-expmk2m1z
      ! Unknown: A_R2^(i-1)
      MA(4*i-3,4*(i-1)  )=-exppk2m1z
      !
      ! Equation: phi^(i)*(U_3^(i)-u_3^(i)) = phi^(i-1)*(U_3^(i-1)-u_3^(i-1))
      !
      ! Unknown: A_I1^(i)
      MA(4*i-2,4*i-3)=phi(i)*(phif(1,i)-1.d0)*expmk1z
      ! Unknown: A_R1^(i)
      MA(4*i-2,4*i-2)=phi(i)*(phif(1,i)-1.d0)*exppk1z
      ! Unknown: A_I2^(i)
      MA(4*i-2,4*i-1)=phi(i)*(phif(2,i)-1.d0)*expmk2z
      ! Unknown: A_R2^(i)
      MA(4*i-2,4*i  )=phi(i)*(phif(2,i)-1.d0)*exppk2z
      ! Unknown: A_I1^(i-1)
      MA(4*i-2,4*(i-1)-3)=-phi(i-1)*(phif(1,i-1)-1.d0)*expmk1m1z
      ! Unknown: A_R1^(i-1)
      MA(4*i-2,4*(i-1)-2)=-phi(i-1)*(phif(1,i-1)-1.d0)*exppk1m1z
      ! Unknown: A_I2^(i-1)
      MA(4*i-2,4*(i-1)-1)=-phi(i-1)*(phif(2,i-1)-1.d0)*expmk2m1z
      ! Unknown: A_R2^(i-1)
      MA(4*i-2,4*(i-1)  )=-phi(i-1)*(phif(2,i-1)-1.d0)*exppk2m1z
      !
      ! Equation: tau_33^(i)+tau^(i)=tau_33^(i-1)+tau^(i-1)
      !
      ! Unknown: A_I1^(i)
      MA(4*i-1,4*i-3)=-(lambda(i)+Q(i)**2/R(i)+2.d0*mu(i)+phif(1,i)*(Q(i)+R(i))+Q(i))*c_im*k(1,i)*expmk1z
      ! Unknown: A_R1^(i)
      MA(4*i-1,4*i-2)= (lambda(i)+Q(i)**2/R(i)+2.d0*mu(i)+phif(1,i)*(Q(i)+R(i))+Q(i))*c_im*k(1,i)*exppk1z
      ! Unknown: A_I2^(i)
      MA(4*i-1,4*i-1)=-(lambda(i)+Q(i)**2/R(i)+2.d0*mu(i)+phif(2,i)*(Q(i)+R(i))+Q(i))*c_im*k(2,i)*expmk2z
      ! Unknown: A_R2^(i)
      MA(4*i-1,4*i  )= (lambda(i)+Q(i)**2/R(i)+2.d0*mu(i)+phif(2,i)*(Q(i)+R(i))+Q(i))*c_im*k(2,i)*exppk2z
      ! Unknown: A_I1^(i-1)
      MA(4*i-1,4*(i-1)-3)= (lambda(i-1)+Q(i-1)**2/R(i-1)+2.d0*mu(i-1)+phif(1,i-1)*(Q(i-1)+R(i-1))+Q(i-1))*c_im*k(1,i-1)*expmk1m1z
      ! Unknown: A_R1^(i-1)
      MA(4*i-1,4*(i-1)-2)=-(lambda(i-1)+Q(i-1)**2/R(i-1)+2.d0*mu(i-1)+phif(1,i-1)*(Q(i-1)+R(i-1))+Q(i-1))*c_im*k(1,i-1)*exppk1m1z
      ! Unknown: A_I2^(i-1)
      MA(4*i-1,4*(i-1)-1)= (lambda(i-1)+Q(i-1)**2/R(i-1)+2.d0*mu(i-1)+phif(2,i-1)*(Q(i-1)+R(i-1))+Q(i-1))*c_im*k(2,i-1)*expmk2m1z
      ! Unknown: A_R2^(i-1)
      MA(4*i-1,4*(i-1)  )=-(lambda(i-1)+Q(i-1)**2/R(i-1)+2.d0*mu(i-1)+phif(2,i-1)*(Q(i-1)+R(i-1))+Q(i-1))*c_im*k(2,i-1)*exppk2m1z
      !
      ! Equation: -tau^(i)/phi(i)+tau^(i-1)/phi(i-1)=k^(i)*i*omega*phi^(i)*(U_3^(i)-u_3^(i))
      !
      ! Unknown: A_I1^(i)
      MA(4*i,4*i-3)=( (Q(i)+phif(1,i)*R(i))*c_im*k(1,i)/phi(i)-ic_k(i)*c_im*omega*phi(i)*(phif(1,i)-1.d0))*expmk1z
      ! Unknown: A_R1^(i)
      MA(4*i,4*i-2)=(-(Q(i)+phif(1,i)*R(i))*c_im*k(1,i)/phi(i)-ic_k(i)*c_im*omega*phi(i)*(phif(1,i)-1.d0))*exppk1z
      ! Unknown: A_I2^(i)
      MA(4*i,4*i-1)=( (Q(i)+phif(2,i)*R(i))*c_im*k(2,i)/phi(i)-ic_k(i)*c_im*omega*phi(i)*(phif(2,i)-1.d0))*expmk2z
      ! Unknown: A_R2^(i)
      MA(4*i,4*i  )=(-(Q(i)+phif(2,i)*R(i))*c_im*k(2,i)/phi(i)-ic_k(i)*c_im*omega*phi(i)*(phif(2,i)-1.d0))*exppk2z
      ! Unknown: A_I1^(i-1)
      MA(4*i,4*(i-1)-3)=-(Q(i-1)+phif(1,i-1)*R(i-1))*c_im*k(1,i-1)/phi(i-1)*expmk1m1z
      ! Unknown: A_R1^(i-1)
      MA(4*i,4*(i-1)-2)= (Q(i-1)+phif(1,i-1)*R(i-1))*c_im*k(1,i-1)/phi(i-1)*exppk1m1z
      ! Unknown: A_I2^(i-1)
      MA(4*i,4*(i-1)-1)=-(Q(i-1)+phif(2,i-1)*R(i-1))*c_im*k(2,i-1)/phi(i-1)*expmk2m1z
      ! Unknown: A_R2^(i-1)
      MA(4*i,4*(i-1)  )= (Q(i-1)+phif(2,i-1)*R(i-1))*c_im*k(2,i-1)/phi(i-1)*exppk2m1z
    end do
    !
    ! Solve LSE using LAPACK routines with scaling and refining
    !
    call fbem_solve_lse_c(4*nl,MA,1,VB)
    !
    ! Unmap solution
    !
    do i=1,nl
      A(1,1,i)=VB(4*i-3,1)
      A(1,2,i)=VB(4*i-2,1)
      A(2,1,i)=VB(4*i-1,1)
      A(2,2,i)=VB(4*i  ,1)
    end do
    ! Save amplitudes of each wave for each variable
    do i=1,nl
      ! u1
      VA(1,1,1,i)=(0.d0,0.d0)
      VA(1,2,1,i)=(0.d0,0.d0)
      VA(2,1,1,i)=(0.d0,0.d0)
      VA(2,2,1,i)=(0.d0,0.d0)
      ! u2
      VA(1,1,2,i)=A(1,1,i)
      VA(1,2,2,i)=A(1,2,i)
      VA(2,1,2,i)=A(2,1,i)
      VA(2,2,2,i)=A(2,2,i)
      ! U1
      VA(1,1,3,i)=(0.d0,0.d0)
      VA(1,2,3,i)=(0.d0,0.d0)
      VA(2,1,3,i)=(0.d0,0.d0)
      VA(2,2,3,i)=(0.d0,0.d0)
      ! U2
      VA(1,1,4,i)=phif(1,i)*A(1,1,i)
      VA(1,2,4,i)=phif(1,i)*A(1,2,i)
      VA(2,1,4,i)=phif(2,i)*A(2,1,i)
      VA(2,2,4,i)=phif(2,i)*A(2,2,i)
      ! tau11
      VA(1,1,5,i)=-(lambda(i)+Q(i)**2/R(i)+phif(1,i)*Q(i))*c_im*k(1,i)*A(1,1,i)
      VA(1,2,5,i)= (lambda(i)+Q(i)**2/R(i)+phif(1,i)*Q(i))*c_im*k(1,i)*A(1,2,i)
      VA(2,1,5,i)=-(lambda(i)+Q(i)**2/R(i)+phif(2,i)*Q(i))*c_im*k(2,i)*A(2,1,i)
      VA(2,2,5,i)= (lambda(i)+Q(i)**2/R(i)+phif(2,i)*Q(i))*c_im*k(2,i)*A(2,2,i)
      ! tau12
      VA(1,1,6,i)=(0.d0,0.d0)
      VA(1,2,6,i)=(0.d0,0.d0)
      VA(2,1,6,i)=(0.d0,0.d0)
      VA(2,2,6,i)=(0.d0,0.d0)
      ! tau22
      VA(1,1,7,i)=-(lambda(i)+Q(i)**2/R(i)+2.d0*mu(i)+phif(1,i)*Q(i))*c_im*k(1,i)*A(1,1,i)
      VA(1,2,7,i)= (lambda(i)+Q(i)**2/R(i)+2.d0*mu(i)+phif(1,i)*Q(i))*c_im*k(1,i)*A(1,2,i)
      VA(2,1,7,i)=-(lambda(i)+Q(i)**2/R(i)+2.d0*mu(i)+phif(2,i)*Q(i))*c_im*k(2,i)*A(2,1,i)
      VA(2,2,7,i)= (lambda(i)+Q(i)**2/R(i)+2.d0*mu(i)+phif(2,i)*Q(i))*c_im*k(2,i)*A(2,2,i)
      ! tau
      VA(1,1,8,i)=-(Q(i)+phif(1,i)*R(i))*c_im*k(1,i)*A(1,1,i)
      VA(1,2,8,i)= (Q(i)+phif(1,i)*R(i))*c_im*k(1,i)*A(1,2,i)
      VA(2,1,8,i)=-(Q(i)+phif(2,i)*R(i))*c_im*k(2,i)*A(2,1,i)
      VA(2,2,8,i)= (Q(i)+phif(2,i)*R(i))*c_im*k(2,i)*A(2,2,i)
    end do
    !
    ! S-wave
    !
    ! Incident wave amplitude (free-surface with unitary displacements: topr=0, u=1, tau=0)
    A(3,1,1)=0.5d0
    ! Reflected wave amplitude (free-surface with unitary: topr=0, u=1, tau=0)
    A(3,2,1)=0.5d0
    do i=2,nl
      TR=(mu(i-1)*k(3,i-1))/(mu(i)*k(3,i))
      A(3,1,i)=0.5d0*exp( c_im*k(3,i)*topr(i))*((1.d0+TR)*exp(-c_im*k(3,i-1)*topr(i))*A(3,1,i-1)&
                                               +(1.d0-TR)*exp( c_im*k(3,i-1)*topr(i))*A(3,2,i-1))
      A(3,2,i)=0.5d0*exp(-c_im*k(3,i)*topr(i))*((1.d0-TR)*exp(-c_im*k(3,i-1)*topr(i))*A(3,1,i-1)&
                                               +(1.d0+TR)*exp( c_im*k(3,i-1)*topr(i))*A(3,2,i-1))
    end do
    ! Save amplitudes of each wave for each variable
    do i=1,nl
      ! u1
      VA(3,1,1,i)=A(3,1,i)
      VA(3,2,1,i)=A(3,2,i)
      ! u2
      VA(3,1,2,i)=(0.d0,0.d0)
      VA(3,2,2,i)=(0.d0,0.d0)
      ! U1
      VA(3,1,3,i)=-hatrho12(i)/hatrho22(i)*A(3,1,i)
      VA(3,2,3,i)=-hatrho12(i)/hatrho22(i)*A(3,2,i)
      ! U2
      VA(3,1,4,i)=(0.d0,0.d0)
      VA(3,2,4,i)=(0.d0,0.d0)
      ! tau11
      VA(3,1,5,i)=(0.d0,0.d0)
      VA(3,2,5,i)=(0.d0,0.d0)
      ! tau12
      VA(3,1,6,i)=-mu(i)*c_im*k(3,i)*A(3,1,i)
      VA(3,2,6,i)= mu(i)*c_im*k(3,i)*A(3,2,i)
      ! tau22
      VA(3,1,7,i)=(0.d0,0.d0)
      VA(3,2,7,i)=(0.d0,0.d0)
      ! tau
      VA(3,1,8,i)=(0.d0,0.d0)
      VA(3,2,8,i)=(0.d0,0.d0)
    end do
  end subroutine fbem_harpor_vertical_plane_wave_layered_halfspace_setup


  subroutine fbem_harpor_vertical_plane_wave_layered_halfspace(Rn,nl,top,VA,k,type_of_wave,tol,layer,x,n,ui,ti)
    implicit none
    ! I/O
    integer              :: Rn
    integer              :: nl           !! Number of layers
    real(kind=real64)    :: top(nl)      !! Top coordinate of each layer
    complex(kind=real64) :: VA(3,2,8,nl) !! Amplitudes of each wave for each variable: u1,u2,U1,U2,tau11,tau12,tau22,tau
    complex(kind=real64) :: k(3,nl)      !! Wavenumbers for each layer
    integer              :: type_of_wave
    real(kind=real64)    :: tol
    integer              :: layer
    real(kind=real64)    :: x(rn)
    real(kind=real64)    :: n(rn)
    complex(kind=real64) :: ui(0:rn)
    complex(kind=real64) :: ti(0:rn)
    ! Local
    integer :: h,v
    complex(kind=real64) :: eikx(3,2)
    !
    ! Initialize
    !
    ui=(0.d0,0.d0)
    ti=(0.d0,0.d0)
    select case (Rn)
      case (2)
        h=1
        v=2
      case (3)
        select case (type_of_wave)
          case (fbem_harpor_shx_wave)
            h=1
          case (fbem_harpor_shy_wave)
            h=2
        end select
        v=3
      case  default
        stop 'fbem_harpor_vertical_plane_wave_layered_halfspace: Rn not valid'
    end select
    !
    ! Check
    !
    if (.not.((layer.ge.1).and.(layer.le.nl))) then
      stop 'fbem_harpor_vertical_plane_wave_layered_halfspace: Invalid selected layer'
    end if
    if (layer.eq.nl) then
      if (.not.(x(v).lt.(top(layer)+tol))) then
        stop 'fbem_harpor_vertical_plane_wave_layered_halfspace: Coordinate does not belong to the selected layer'
      end if
    else
      if (.not.((x(v).lt.(top(layer)+tol)).and.(x(v).gt.(top(layer+1)-tol)))) then
        stop 'fbem_harpor_vertical_plane_wave_layered_halfspace: Coordinate does not belong to the selected layer'
      end if
    end if
    ! exp(+-ikz)
    eikx(1,1)=exp(-c_im*k(1,layer)*(x(v)-top(1)))
    eikx(1,2)=exp( c_im*k(1,layer)*(x(v)-top(1)))
    eikx(2,1)=exp(-c_im*k(2,layer)*(x(v)-top(1)))
    eikx(2,2)=exp( c_im*k(2,layer)*(x(v)-top(1)))
    eikx(3,1)=exp(-c_im*k(3,layer)*(x(v)-top(1)))
    eikx(3,2)=exp( c_im*k(3,layer)*(x(v)-top(1)))
    ! Calculate
    select case (type_of_wave)
    case (fbem_harpor_p1_wave)
      ui(0)= VA(1,1,8,layer)*eikx(1,1)+VA(1,2,8,layer)*eikx(1,2)+VA(2,1,8,layer)*eikx(2,1)+VA(2,2,8,layer)*eikx(2,2)
      ui(v)= VA(1,1,2,layer)*eikx(1,1)+VA(1,2,2,layer)*eikx(1,2)+VA(2,1,2,layer)*eikx(2,1)+VA(2,2,2,layer)*eikx(2,2)
      ti(0)=(VA(1,1,4,layer)*eikx(1,1)+VA(1,2,4,layer)*eikx(1,2)+VA(2,1,4,layer)*eikx(2,1)+VA(2,2,4,layer)*eikx(2,2))*n(v)
      ti(1)=(VA(1,1,5,layer)*eikx(1,1)+VA(1,2,5,layer)*eikx(1,2)+VA(2,1,5,layer)*eikx(2,1)+VA(2,2,5,layer)*eikx(2,2))*n(1)
      if (Rn.eq.3) then
      ti(2)=(VA(1,1,5,layer)*eikx(1,1)+VA(1,2,5,layer)*eikx(1,2)+VA(2,1,5,layer)*eikx(2,1)+VA(2,2,5,layer)*eikx(2,2))*n(2)
      end if
      ti(v)=(VA(1,1,7,layer)*eikx(1,1)+VA(1,2,7,layer)*eikx(1,2)+VA(2,1,7,layer)*eikx(2,1)+VA(2,2,7,layer)*eikx(2,2))*n(v)
    case (fbem_harpor_shx_wave,fbem_harpor_shy_wave)
      ui(h)= VA(3,1,1,layer)*eikx(3,1)+VA(3,2,1,layer)*eikx(3,2)
      ti(0)=(VA(3,1,3,layer)*eikx(3,1)+VA(3,2,3,layer)*eikx(3,2))*n(h)
      ti(h)=(VA(3,1,6,layer)*eikx(3,1)+VA(3,2,6,layer)*eikx(3,2))*n(v)
      ti(v)=(VA(3,1,6,layer)*eikx(3,1)+VA(3,2,6,layer)*eikx(3,2))*n(h)
    end select
  end subroutine fbem_harpor_vertical_plane_wave_layered_halfspace

  !
  !   Calcula campo de movimientos y tensores de tension del
  !   campo incidente en un semiespacio POROELASTICO.
  !   (solo incidencia vertical)
  !
  subroutine fbem_harpor_incident_plane_wave(d,omega,kp1,kp2,ks,lambda,mu,cq,cr,ro11bar,ro22bar,ro12bar,&
                                             space,type_of_wave,x,z_fs,n,ui,ti)
    implicit none
    ! I/O variables
    integer    :: d
    real*8     :: omega
    complex*16 :: kp1
    complex*16 :: kp2
    complex*16 :: ks
    complex*16 :: lambda
    complex*16 :: mu
    complex*16 :: cq
    complex*16 :: cr
    complex*16 :: ro11bar
    complex*16 :: ro22bar
    complex*16 :: ro12bar
    integer    :: space
    integer    :: type_of_wave
    real*8     :: x(3)
    real*8     :: z_fs
    real*8     :: n(3)
    complex*16 :: ui(0:3)
    complex*16 :: ti(0:3)
    ! Local variables
    integer    :: i
    complex*16 :: tau
    complex*16 :: uis(3)
    complex*16 :: uif(3)
    complex*16 :: sigma(6)
    complex*16 :: epss(6)
    complex*16 :: epsf(6)
    complex*16 :: cimag
    complex*16 :: zp1,zp2,zs
    complex*16 :: cte
    complex*16 :: fhi1s,fhi1f,fhi2s,fhi2f  ! autovalores de onda p porosa
    complex*16 :: ai,ar
    real*8     :: xr,yr,zr



    !
    ! Inicializacion
    !
    ! Si el problema es 2d:
    !   x2 (y) en 2d == x3 (z) en 3d
    !   x1 (x) en 2d == x2 (y) en 3d
    if (d.eq.2) then
      x(3)=x(2)
      n(3)=n(2)
      x(2)=x(1)
      n(2)=n(1)
      x(1)=0.d0
      n(1)=0.d0
      select case (type_of_wave)
        case (fbem_harpor_shx_wave)
          type_of_wave=fbem_harpor_shy_wave
        case (fbem_harpor_shy_wave)
          stop 'harpor_incident_field: it does not make sense for a 2D problem having a shy wave'
      end select
    end if
    ! Unidad imaginaria
    cimag=(0.d0,1.d0)
    ! Inicializacion de desplazamientos, deformaciones y tensiones
    uis  =(0.d0,0.d0)
    uif  =(0.d0,0.d0)
    tau  =(0.d0,0.d0)
    sigma=(0.d0,0.d0)
    epss =(0.d0,0.d0)
    epsf =(0.d0,0.d0)
    ! Coordenadas relativas
    xr=x(1)
    yr=x(2)
    zr=x(3)-z_fs
    ! Numeros de onda multiplicados por i
    zp1=cimag*kp1
    zp2=cimag*kp2
    zs =cimag*ks
    ! Algunas constantes
    cte=ro12bar/ro22bar
    fhi1s=(1.d0,0.d0)
    fhi2s=(1.d0,0.d0)
    fhi1f=kp1**2*(lambda+2.d0*mu+cq**2/cr)-omega**2*ro11bar
    fhi1f=fhi1f/(omega**2*ro12bar-kp1**2*cq)
    fhi2f=kp2**2*(lambda+2.d0*mu+cq**2/cr)-omega**2*ro11bar
    fhi2f=fhi2f/(omega**2*ro12bar-kp2**2*cq)
    ! Amplitudes de las ondas
    select case (type_of_wave)
      case (fbem_harpor_shx_wave,fbem_harpor_shy_wave)
        !ai=0.5*(1./(1.-beta*(1.+cte)))  ! desplazamiento promedio = 1
        ai=(0.5d0,0.0d0)                     ! desplazamiento esqueleto solido = 1
        ar=ai
      case (fbem_harpor_p1_wave)
        !ai=0.5*(1./((1.-beta)*fhi1s+beta*fhi1f))  ! desplazamiento promedio = 1
        ai=(0.5d0,0.0d0)                                   ! desplazamiento esqueleto solido = 1
        ar=ai
    end select
    if (space.eq.1) ar=(0.d0,0.d0)
    ! Solucion
    select case (type_of_wave)
      case (fbem_harpor_shx_wave)
        uis(1)=ai*cdexp(-zs*zr)+ar*cdexp(zs*zr)
        uif(1)=-cte*uis(1)
        epss(3)=(-ai*cdexp(-zs*zr)+ar*cdexp(zs*zr))*zs
        sigma(3)=mu*epss(3)
      case (fbem_harpor_shy_wave)
        uis(2)=ai*cdexp(-zs*zr)+ar*cdexp(zs*zr)
        uif(2)=-cte*uis(2)
        epss(5)=(-ai*cdexp(-zs*zr)+ar*cdexp(zs*zr))*zs
        sigma(5)=mu*epss(5)
      case(fbem_harpor_p1_wave)
        uis(3)=fhi1s*(ai*cdexp(-zp1*zr)+ar*cdexp(zp1*zr))
        uif(3)=fhi1f*(ai*cdexp(-zp1*zr)+ar*cdexp(zp1*zr))
        epss(6)=zp1*fhi1s*(-ai*cdexp(-zp1*zr)+ar*cdexp(zp1*zr))
        epsf(6)=zp1*fhi1f*(-ai*cdexp(-zp1*zr)+ar*cdexp(zp1*zr))
        sigma(6)=(lambda+2.d0*mu+cq**2/cr)*epss(6)+cq*epsf(6)
        sigma(1)=(lambda+cq**2/cr)*epss(6)+cq*epsf(6)
        sigma(4)=sigma(1)
        tau=cq*epss(6)+cr*epsf(6)
      case(fbem_harpor_p2_wave)
        uis(3)=fhi2s*(ai*cdexp(-zp2*zr)+ar*cdexp(zp2*zr))
        uif(3)=fhi2f*(ai*cdexp(-zp2*zr)+ar*cdexp(zp2*zr))
        epss(6)=zp1*fhi2s*(-ai*cdexp(-zp2*zr)+ar*cdexp(zp2*zr))
        epsf(6)=zp1*fhi2f*(-ai*cdexp(-zp2*zr)+ar*cdexp(zp2*zr))
        sigma(6)=(lambda+2.d0*mu+cq**2/cr)*epss(6)+cq*epsf(6)
        sigma(1)=(lambda+cq**2/cr)*epss(6)+cq*epsf(6)
        sigma(4)=sigma(1)
        tau=cq*epss(6)+cr*epsf(6)
      case default
        stop 'fbem_harpor_incident_plane_wave: not valid type_of_wave'
    end select
    ! Componentes finales, teniendo en cuenta que los tensores simetricos estan almacenados en vectores con:
    ! 1: xx
    ! 2: xy
    ! 3: xz
    ! 4: yy
    ! 5: yz
    ! 6: zz
    ui(0)=tau
    ui(1)=uis(1)
    ui(2)=uis(2)
    ui(3)=uis(3)
    ti(0)=uif  (1)*n(1)+uif  (2)*n(2)+uif  (3)*n(3)
    ti(1)=sigma(1)*n(1)+sigma(2)*n(2)+sigma(3)*n(3)
    ti(2)=sigma(2)*n(1)+sigma(4)*n(2)+sigma(5)*n(3)
    ti(3)=sigma(3)*n(1)+sigma(5)*n(2)+sigma(6)*n(3)
    !
    !  Si el problema es 2D:
    !    x2 [y] en 2d == x3 [z] en 3d
    !    x1 [x] en 2d == x2 [y] en 3d
    !
    if (d.eq.2) then
      ui(1)=ui(2)
      ti(1)=ti(2)
      x(1)=x(2)
      n(1)=n(2)
      ui(2)=ui(3)
      ti(2)=ti(3)
      x(2)=x(3)
      n(2)=n(3)
    end if
  end subroutine fbem_harpor_incident_plane_wave

  !! Phurkhao incident field
  subroutine fbem_phurkhao(a,rhof,rhos,lambda,mu,phi,rhoa,R,Q,b,omega,omega1,x,n,wavetype,u,t)
    implicit none
    ! I/O
    real(kind=real64)    :: a
    real(kind=real64)    :: rhof
    real(kind=real64)    :: rhos
    real(kind=real64)    :: lambda
    real(kind=real64)    :: mu
    real(kind=real64)    :: phi
    real(kind=real64)    :: rhoa
    real(kind=real64)    :: R
    real(kind=real64)    :: Q
    real(kind=real64)    :: b
    real(kind=real64)    :: omega
    real(kind=real64)    :: omega1
    real(kind=real64)    :: x(3)
    real(kind=real64)    :: n(3)
    integer              :: wavetype
    complex(kind=real64) :: u(0:3)
    complex(kind=real64) :: t(0:3)
    ! Local
    real(kind=real64)    :: rho
    real(kind=real64)    :: alpha
    real(kind=real64)    :: M
    real(kind=real64)    :: f
    real(kind=real64)    :: qa
    real(kind=real64)    :: alpha_f
    real(kind=real64)    :: ba
    real(kind=real64)    :: beta
    real(kind=real64)    :: omegaa
    real(kind=real64)    :: vc, vs
    real(kind=real64)    :: n1, n2, delta
    complex(kind=real64) :: z1, z2, tmpz
    complex(kind=real64) :: gamma1, gamma2, gamma1p, gamma1m, gamma2p, gamma2m
    complex(kind=real64) :: ta, tb, tc, phi0, phi0f, expy, expmy
    ! Initialize
    u=(0.d0,0.d0)
    t=(0.d0,0.d0)
    ! Calculate constants
    f=phi
    rho=phi*rhof+(1.d0-phi)*rhos
    M=R/phi**2
    alpha=phi*(1.d0+Q/R)
    alpha_f=M/mu
    ba=b*a/(phi**2*sqrt(mu*rho))
    qa=rhof/rho
    vc=sqrt((lambda+2.d0*mu+alpha**2*M)/rho)
    vs=sqrt(mu/rho)
    beta=vc/vs
    omegaa=omega*a/vs
  !  ! Wavenumber calculation (my way)
  !  ta=alpha_f*(beta**2-alpha**2*alpha_f)
  !  tb=((2.d0*alpha*alpha_f*qa-alpha_f-beta**2*qa/f)*omegaa+c_im*ba*beta**2)*omegaa
  !  tc=(-c_im*ba+qa*omegaa*(1.0d0/f-qa))*omegaa**3
  !  gamma1=zsqrt(0.5d0*(-tb-zsqrt(tb*tb-4.d0*ta*tc))/ta)
  !  gamma2=zsqrt(0.5d0*(-tb+zsqrt(tb*tb-4.d0*ta*tc))/ta)
  !  gamma1=dble(gamma1)-c_im*dimag(gamma1)
  !  gamma2=dble(gamma2)-c_im*dimag(gamma2)
    ! Wavenumber calculation (as Phurkhao)
    n1=(2.d0*alpha*alpha_f*qa-alpha_f-qa*beta**2/f)/(alpha_f*(1.d0-alpha**2*alpha_f/beta**2))
    n2=(beta**2*(qa/f-qa**2))/(alpha_f*(1.d0-alpha**2*alpha_f/beta**2))
    delta=ba*beta**2/(alpha_f*(1.d0-alpha**2*alpha_f/beta**2))/omegaa
    z1=0.5d0*(-n1+c_im*delta-sqrt((n1-c_im*delta)**2-4.d0*(n2+c_im*delta)))
    z2=0.5d0*(-n1+c_im*delta+sqrt((n1-c_im*delta)**2-4.d0*(n2+c_im*delta)))
    if (real(z1).gt.real(z2)) then
      tmpz=z1
      z1=z2
      z2=tmpz
    end if
    gamma1=omegaa/beta*sqrt(z1)
    gamma2=omegaa/beta*sqrt(z2)
    if (real(gamma1).lt.0.d0) gamma1=-gamma1
    if (real(gamma2).lt.0.d0) gamma2=-gamma2
    !
    ! Campo obtenido al sumar una onda que sube y otra que baja (tanto parte real como imaginaria de los desplazamientos simétricas)
    !
    ! Depending on the wave type
    select case (wavetype)
      !
      ! P1-wave
      !
      case (fbem_harpor_p1_wave)
        phi0=1.d0
        !phi0f=(beta**2*gamma1**2-omegaa**2)/(qa*omegaa**2-alpha*alpha_f*gamma1**2) ! (my way)
        phi0f=(1.0d0-z1)/(alpha*alpha_f*z1/beta**2-qa) ! (as Phurkhao)
        expy=zexp(-c_im*gamma1*x(2))
        expmy=zexp(c_im*gamma1*x(2))
        u(0)=-f*alpha_f*(alpha*phi0+phi0f)*gamma1**2*(expy+expmy)*mu
        u(1)=(0.d0,0.d0)
        u(2)=-c_im*gamma1*phi0*(expy-expmy)*a
        t(0)=-c_im*gamma1*(phi0f/f+phi0)*(expy-expmy)*a*n(2)
        t(1)=(0.d0,0.d0)
        t(2)=(-(beta**2*phi0+alpha*alpha_f*phi0f)*gamma1**2*(expy+expmy)*mu-u(0))*n(2)
        ! Dimensionless wave field
        u=u/(-omegaa**2*mu*2.d0)
        t=t/(-omegaa**2*mu*2.d0)
      !
      ! P2-wave
      !
      case (fbem_harpor_p2_wave)
        !phi0=(alpha*alpha_f*gamma2**2-qa*omegaa**2)/(omegaa**2-beta**2*gamma2**2) ! (my way)
        phi0=(alpha*alpha_f*z2/beta**2-qa)/(1.0d0-z2) ! (as Phurkhao)
        phi0f=1.0d0
        expy=zexp(-c_im*gamma2*x(2))
        expmy=zexp(c_im*gamma2*x(2))
        u(0)=-f*alpha_f*(alpha*phi0+phi0f)*gamma2**2*(expy+expmy)*mu
        u(1)=(0.d0,0.d0)
        u(2)=-c_im*gamma2*phi0*(expy-expmy)*a
        t(0)=-c_im*gamma2*(phi0f/f+phi0)*(expy-expmy)*a*n(2)
        t(1)=(0.d0,0.d0)
        t(2)=(-(beta**2*phi0+alpha*alpha_f*phi0f)*gamma2**2*(expy+expmy)*mu-u(0))*n(2)
        ! Dimensionless wave field
        u=u/(-mu*2.d0)
        t=t/(-mu*2.d0)
    end select
  end subroutine fbem_phurkhao

  !! Rayleigh wave in poroelastic half-space. Propagation along the x positive direction
  subroutine fbem_harpor_incident_rayleigh_wave(Rn,lambda,mu,Q,R,phi,rhos,rhof,rhoa,b,omega,xn,nn,ui,ti)
    implicit none
    ! I/O
    integer              :: Rn
    complex(kind=real64) :: lambda, mu, Q, R
    real(kind=real64)    :: phi, rhos, rhof, rhoa, b
    real(kind=real64)    :: omega
    real(kind=real64)    :: xn(3)
    real(kind=real64)    :: nn(3)
    complex(kind=real64) :: ui(0:3)
    complex(kind=real64) :: ti(0:3)
    ! Local
    complex(kind=real64) :: kp1, kp2, ks, kr, krp1, krp2, krs
    complex(kind=real64) :: hatrho11, hatrho12, hatrho22, D(2)
    complex(kind=real64) :: AP1, AP2, AS, Auy(4,3), Asy(4,3), ekry(3), ekrx
    complex(kind=real64) :: up(4), sp(4)
    complex(kind=real64) :: up0(4), sp0(4)
    real(kind=real64)    :: x(3)
    !
    ! Initialize
    !
    ui=0.d0
    ti=0.d0
    !
    ! CALCULATE WAVENUMBERS
    !
    call harpor_calculate_bulk_wavenumbers(lambda,mu,Q,R,phi,rhos,rhof,rhoa,b,omega,kp1,kp2,ks)
    call harpor_calculate_permeable_rayleigh_wavenumber(lambda,mu,Q,R,phi,rhos,rhof,rhoa,b,omega,1.d-14,kr)
    krp1=sqrt(kr**2-kp1**2)
    krp2=sqrt(kr**2-kp2**2)
    krs =sqrt(kr**2-ks**2 )
    if (dreal(krp1).lt.0.d0) krp1=-krp1
    if (dreal(krp2).lt.0.d0) krp2=-krp2
    if (dreal(krs ).lt.0.d0) krs =-krs
    !
    ! MATRICES OF AMPLITUDES
    ! Rows: variables
    ! Columns: RP1, RP2, RS
    !
    ! Auxiliary constants
    hatrho11=(1.d0-phi)*rhos + rhoa - c_im*b/omega
    hatrho22=       phi*rhof + rhoa - c_im*b/omega
    hatrho12=                - rhoa + c_im*b/omega
    D(1)=((lambda+2.d0*mu)*kp1**2-omega**2*(hatrho11-Q/R*hatrho12))/(omega**2*(hatrho12-Q/R*hatrho22))
    D(2)=((lambda+2.d0*mu)*kp2**2-omega**2*(hatrho11-Q/R*hatrho12))/(omega**2*(hatrho12-Q/R*hatrho22))
    !
    ! Potentials amplitudes
    !
    AP1=1.d0
    AP2=-kp1**2/kp2**2*(Q+R*D(1))/(Q+R*D(2))
    AS=2.d0*c_im*kr*(krp1+AP2*krp2)/(2.d0*kr**2-ks**2)
    !
    ! Displacements amplitudes. Variables: u_x, u_y, U_x, U_y
    !
    ! u_x
    Auy(1,1)=-c_im*kr
    Auy(1,2)=-c_im*kr
    Auy(1,3)= krs
    ! u_y
    Auy(2,1)= krp1
    Auy(2,2)= krp2
    Auy(2,3)= c_im*kr
    ! U_x
    Auy(3,1)=-c_im*kr*D(1)
    Auy(3,2)=-c_im*kr*D(2)
    Auy(3,3)=-krs*hatrho12/hatrho22
    ! U_y
    Auy(4,1)= krp1*D(1)
    Auy(4,2)= krp2*D(2)
    Auy(4,3)=-c_im*kr*hatrho12/hatrho22
    ! Final matrix
    Auy(:,1)=Auy(:,1)*AP1
    Auy(:,2)=Auy(:,2)*AP2
    Auy(:,3)=Auy(:,3)*AS
    !
    ! Stresses amplitudes. Variables: tau_xx, tau_xy, tau_yy, tau
    !
    ! tau_xx
    Asy(1,1)=-(lambda+Q**2/R+Q*D(1))*kp1**2-2.d0*mu*kr**2
    Asy(1,2)=-(lambda+Q**2/R+Q*D(2))*kp2**2-2.d0*mu*kr**2
    Asy(1,3)=-2.d0*mu*c_im*kr*krs
    ! tau_xy
    Asy(2,1)=-2.d0*mu*c_im*kr*krp1
    Asy(2,2)=-2.d0*mu*c_im*kr*krp2
    Asy(2,3)= mu*(2.d0*kr**2-ks**2)
    ! tau_yy
    Asy(3,1)=-(lambda+2.d0*mu+Q**2/R+Q*D(1))*kp1**2+2.d0*mu*kr**2
    Asy(3,2)=-(lambda+2.d0*mu+Q**2/R+Q*D(2))*kp2**2+2.d0*mu*kr**2
    Asy(3,3)= 2.d0*mu*c_im*kr*krs
    ! tau
    Asy(4,1)=-(Q+R*D(1))*kp1**2
    Asy(4,2)=-(Q+R*D(2))*kp2**2
    Asy(4,3)= 0.d0
    ! Final matrix
    Asy(:,1)=Asy(:,1)*AP1
    Asy(:,2)=Asy(:,2)*AP2
    Asy(:,3)=Asy(:,3)*AS
    !
    ! INCIDENT FIELD IN TERMS OF POTENTIALS at (0,0)
    !
    ! Point
    x=0.d0
    ! Depth exponentials
    ekry(1)=exp(krp1*x(2))
    ekry(2)=exp(krp2*x(2))
    ekry(3)=exp( krs*x(2))
    ! Exponential in the propagation direction
    ekrx=exp(-c_im*kr*x(1))
    ! Calculate displacements and stresses
    up0=matmul(Auy,ekry)*ekrx
    sp0=matmul(Asy,ekry)*ekrx
    !
    ! INCIDENT FIELD IN TERMS OF POTENTIALS at xn
    !
    ! Point
    select case (Rn)
      case (2)
        x=xn
      case (3)
        x(1)=xn(1)
        x(2)=xn(3)
        x(3)=0.d0
    end select
    ! Depth exponentials
    ekry(1)=exp(krp1*x(2))
    ekry(2)=exp(krp2*x(2))
    ekry(3)=exp( krs*x(2))
    ! Exponential in the propagation direction
    ekrx=exp(-c_im*kr*x(1))
    ! Calculate displacements and stresses
    up=matmul(Auy,ekry)*ekrx
    sp=matmul(Asy,ekry)*ekrx
    !
    ! NORMALIZE
    !
    ! Unitary vertical displacements
    up=up/up0(2)
    sp=sp/up0(2)
    !
    ! OUTPUT VARIABLES
    !
    select case (Rn)
      case (2)
        ui(0)=sp(4)
        ui(1)=up(1)
        ui(2)=up(2)
        ui(3)=0.d0
        ti(0)=up(3)*nn(1)+up(4)*nn(2)
        ti(1)=sp(1)*nn(1)+sp(2)*nn(2)
        ti(2)=sp(2)*nn(1)+sp(3)*nn(2)
        ti(3)=0.d0
      case (3)
        ui(0)=sp(4)
        ui(1)=up(1)
        ui(2)=0.d0
        ui(3)=up(2)
        ti(0)=up(3)*nn(1)+up(4)*nn(3)
        ti(1)=sp(1)*nn(1)+sp(2)*nn(3)
        ti(2)=0.d0
        ti(3)=sp(2)*nn(1)+sp(3)*nn(3)
    end select
  end subroutine fbem_harpor_incident_rayleigh_wave

  !! Bulk wavenumbers of poroelastic media
  subroutine harpor_calculate_bulk_wavenumbers(lambda,mu,Q,R,phi,rhos,rhof,rhoa,b,omega,kp1,kp2,ks)
    implicit none
    ! I/O
    complex(kind=real64) :: lambda, mu, Q, R
    real(kind=real64)    :: phi, rhos, rhof, rhoa, b
    real(kind=real64)    :: omega
    complex(kind=real64) :: kp1, kp2, ks
    ! Local
    integer              :: i, n
    complex(kind=real64) :: hatrho11, hatrho22, hatrho12
    complex(kind=real64) :: a0, a1, a2, x1, x2
    complex(kind=real64) :: k(4)
    ! Auxiliary variables
    hatrho11=(1.d0-phi)*rhos + rhoa - c_im*b/omega
    hatrho22=       phi*rhof + rhoa - c_im*b/omega
    hatrho12=                - rhoa + c_im*b/omega
    ! Calculate kp1 and kp2 by solving the biquadratic equation
    a0=omega**4*(hatrho11*hatrho22-hatrho12**2)/(R*(lambda+2.d0*mu))
    a1=omega**2*(hatrho22/R+(hatrho11+(Q/R)**2*hatrho22-2.d0*Q/R*hatrho12)/(lambda+2.d0*mu))
    a2=1.d0
    call fbem_solve_quadratic_equation(a2,-a1,a0,x1,x2)
    x1=sqrt(x1)
    x2=sqrt(x2)
    k(1)= x1
    k(2)= x2
    k(3)=-x1
    k(4)=-x2
    n=0
    do i=1,4
      if (dreal(k(i)).gt.0) then
        n=n+1
        k(n)=k(i)
      end if
    end do
    if (n.gt.2) stop 'Error: this should not happen'
    if (dreal(k(1)).lt.dreal(k(2))) then
      kp1=k(1)
      kp2=k(2)
    else
      kp1=k(2)
      kp2=k(1)
    end if
    ! Calculate ks
    k(1)= omega*sqrt((hatrho11-hatrho12**2/hatrho22)/mu)
    k(2)=-omega*sqrt((hatrho11-hatrho12**2/hatrho22)/mu)
    do i=1,2
      if (real(k(i)).gt.0) then
        ks=k(i)
        exit
      end if
    end do
  end subroutine harpor_calculate_bulk_wavenumbers

  !! Rayleigh surface wavenumber of poroelastic media with permeable half-plane
  subroutine harpor_calculate_permeable_rayleigh_wavenumber(lambda,mu,Q,R,phi,rhos,rhof,rhoa,b,omega,relative_error,kr)
    implicit none
    ! I/O
    complex(kind=real64) :: lambda, mu, Q, R
    real(kind=real64)    :: phi, rhos, rhof, rhoa, b
    real(kind=real64)    :: omega
    real(kind=real64)    :: relative_error
    complex(kind=real64) :: kr
    ! Local
    complex(kind=real64) :: kp1, kp2, ks
    complex(kind=real64) :: r_1, f_1, r_2, f_2, r_n, f_n, dfdr
    integer              :: iteration, mode
    logical              :: next
    real(kind=real64)    :: relax
    !
    ! Mode
    !
    mode=0
    if (abs(b).eq.0.d0) mode=1
    !
    ! Calculate kp1, kp2, ks
    !
    call harpor_calculate_bulk_wavenumbers(lambda,mu,Q,R,phi,rhos,rhof,rhoa,b,omega,kp1,kp2,ks)
    !
    ! Calculate a good approximation based on the Bardet's model (ignoring kp2)
    !
    call harela_calculate_rayleigh_wavenumber(kp1,ks,relative_error,kr)
    !
    ! Refine it by using the characteristic equation solved by the Newton-Raphson Method or Secant Method
    ! Note: r=ks/kr=cr/cs
    !
    !
    ! NEWTON-RAPHSON METHOD
    !
    relax=0.25d0
    r_1=ks/kr
    f_1=harpor_permeable_rayleigh_equation(mode,mu,lambda,kp1,kp2,ks,ks/r_1)
    iteration=0
    next=.true.
    do while (next)
      dfdr=harpor_permeable_rayleigh_equation_d(mode,mu,lambda,kp1,kp2,ks,ks/r_1)
      r_n=r_1-relax*f_1/dfdr
      f_n=harpor_permeable_rayleigh_equation(mode,mu,lambda,kp1,kp2,ks,ks/r_n)
      if ((abs((r_1-r_n)/r_n).lt.relative_error).and.(abs(f_n).lt.relative_error)) then
        next=.false.
      else
        if (iteration.gt.1000) then
          write(error_unit,*) 'Last 2 iterations [Re(r), Im(r), Re(f(r)), Im(f(r))]:'
          write(error_unit,'(4e25.16)') dreal(r_1), dimag(r_1), dreal(f_1), dimag(f_1)
          write(error_unit,'(4e25.16)') dreal(r_n), dimag(r_n), dreal(f_n), dimag(f_n)
          write(*,*) 'WARNING: Rayleigh wavenumber calculation of a permeable poroelastic half-plane has reached too many iterations.'
          exit
        end if
        r_1=r_n
        f_1=f_n
        iteration=iteration+1
      end if
    end do
    kr=ks/r_n
!    !
!    ! SECANT METHOD
!    !
!    r_1=0.95d0*ks/kr
!    f_1=harpor_permeable_rayleigh_equation(mode,mu,lambda,kp1,kp2,ks,ks/r_1)
!    r_2=1.05d0*ks/kr
!    f_2=harpor_permeable_rayleigh_equation(mode,mu,lambda,kp1,kp2,ks,ks/r_2)
!    iteration=0
!    next=.true.
!    do while (next)
!      dfdr=(f_2-f_1)/(r_2-r_1)
!      r_n=r_1-f_1/dfdr
!      f_n=harpor_permeable_rayleigh_equation(mode,mu,lambda,kp1,kp2,ks,ks/r_n)
!      r_1=r_2
!      f_1=f_2
!      r_2=r_n
!      f_2=f_n
!      if ((abs((r_1-r_2)/r_2).lt.relative_error).and.(abs(f_n).lt.relative_error)) then
!        next=.false.
!      else
!        iteration=iteration+1
!        if (iteration.gt.100) then
!          write(error_unit,*) 'Last 2 iterations [Re(r), Im(r), Re(f(r)), Im(f(r))]:'
!          write(error_unit,'(4e25.16)') dreal(r_1), dimag(r_1), dreal(f_1), dimag(f_1)
!          write(error_unit,'(4e25.16)') dreal(r_2), dimag(r_2), dreal(f_2), dimag(f_2)
!          write(*,*) 'WARNING: Rayleigh wavenumber calculation of a permeable poroelastic half-plane has reached too many iterations.'
!          exit
!        end if
!      end if
!    end do
!    kr=ks/r_n
!    write(*,*) r_n, abs((r_1-r_n)/r_n), harpor_permeable_rayleigh_equation(mode,mu,lambda,kp1,kp2,ks,ks/r_n)
!    stop
  end subroutine harpor_calculate_permeable_rayleigh_wavenumber

  function harpor_permeable_rayleigh_equation(mode,mu,lambda,kp1,kp2,ks,kr)
    implicit none
    ! I/O
    complex(kind=real64) :: harpor_permeable_rayleigh_equation
    integer              :: mode !! Mode: 0 (complex), 1 (real part), 0 (imaginary part)
    complex(kind=real64) :: mu, lambda, kp1, kp2, ks, kr
    ! Local
    complex(kind=real64) :: r, r2, H1, H2, G1, G2
    ! Compute
    r =ks/kr
    r2=r**2
    H1=(kp1**2-mu/(lambda+2.d0*mu)*ks**2)/(kp2**2-kp1**2)
    H2=(kp2**2-mu/(lambda+2.d0*mu)*ks**2)/(kp2**2-kp1**2)
    G1=(kp1/ks)**2
    G2=(kp2/ks)**2
    harpor_permeable_rayleigh_equation=(2.d0-r2)**2-4.d0*sqrt(1.d0-r2)*(H2*sqrt(1.d0-G1*r2)-H1*sqrt(1.d0-G2*r2))
    select case (mode)
      case (1)
        harpor_permeable_rayleigh_equation=dreal(harpor_permeable_rayleigh_equation)
      case (2)
        harpor_permeable_rayleigh_equation=dimag(harpor_permeable_rayleigh_equation)
    end select
  end function harpor_permeable_rayleigh_equation

  function harpor_permeable_rayleigh_equation_d(mode,mu,lambda,kp1,kp2,ks,kr)
    implicit none
    ! I/O
    complex(kind=real64) :: harpor_permeable_rayleigh_equation_d
    integer              :: mode !! Mode: 0 (complex), 1 (real part), 0 (imaginary part)
    complex(kind=real64) :: mu, lambda, kp1, kp2, ks, kr
    ! Local
    complex(kind=real64) :: r, r2, H1, H2, G1, G2
    ! Compute
    r =ks/kr
    r2=r**2
    H1=(kp1**2-mu/(lambda+2.d0*mu)*ks**2)/(kp2**2-kp1**2)
    H2=(kp2**2-mu/(lambda+2.d0*mu)*ks**2)/(kp2**2-kp1**2)
    G1=(kp1/ks)**2
    G2=(kp2/ks)**2
    harpor_permeable_rayleigh_equation_d=4.d0*r*(-(2.d0-r2)+1.d0/sqrt(1.d0-r2)*(H2*sqrt(1.d0-G1*r2)-H1*sqrt(1.d0-G2*r2))&
                                                 +sqrt(1.d0-r2)*(H2*G1/sqrt(1.d0-G1*r2)-H1*G2/sqrt(1.d0-G2*r2)))
    select case (mode)
      case (1)
        harpor_permeable_rayleigh_equation_d=dreal(harpor_permeable_rayleigh_equation_d)
      case (2)
        harpor_permeable_rayleigh_equation_d=dimag(harpor_permeable_rayleigh_equation_d)
    end select
  end function harpor_permeable_rayleigh_equation_d

end module fbem_harpor_incident_field
