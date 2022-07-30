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
!! <b> This module implements some incident fields for inviscid fluids. </b>
module fbem_harpot_incident_field

  ! Fortran 2003 standard
  use iso_fortran_env

  ! FBEM library
  use fbem_numerical
  use fbem_string_handling

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all are private
  private

  ! Subroutines
  public :: fbem_harpot_pointwave
  public :: fbem_harpot_planewave
  public :: fbem_harpot_properties

contains

  !! Full-space or half-space point wave for inviscid fluids. For the half-space solution, the point can be
  !! at any coordinate, and any boundary condition  for the (p=0 or Un=0).
  subroutine fbem_harpot_pointwave(rn,omega,rho,c,A,x0,space,np,xp,bc,symconf,xs,x,n,p_inc,Un_inc)
    implicit none
    ! Input/Output
    integer              :: rn            !! Dimensional space: 2 or 3.
    real(kind=real64)    :: omega         !! Angular frequency.
    real(kind=real64)    :: rho           !! Density.
    complex(kind=real64) :: c             !! Wave propagation speed.
    complex(kind=real64) :: A             !! Amplitude of the point source at r=1
    real(kind=real64)    :: x0(rn)        !! Point coordinates.
    integer              :: space         !! Type of space: 1 full-space, 2 half-space.
    integer              :: np            !! Axis corresponding to the normal of the plane (x_np=xp).
    real(kind=real64)    :: xp            !! Position of the plane.
    integer              :: bc            !! Boundary condition at the plane: 0 p=0, 1 Un=0.
    !
    real(kind=real64)    :: symconf(rn)   !! Symmetry configuration vector. It must be meet that symconf(np)=0.
    real(kind=real64)    :: xs(rn)        !! Origin of symmetry/antisymmetry decomposition
    !
    real(kind=real64)    :: x(rn)         !! Observation point
    real(kind=real64)    :: n(rn)         !! Unit normal of the observation point
    !
    complex(kind=real64) :: p_inc         !! Pressure of the incident wave
    complex(kind=real64) :: Un_inc        !! Normal displacement of the incident wave
    ! Local
    complex(kind=real64) :: k           ! Wavenumber
    complex(kind=real64) :: p1          ! p a r=1 de la solucion fundamental
    real(kind=real64)    :: rv(rn), r
    real(kind=real64)    :: drdx(rn)
    real(kind=real64)    :: drdn
    k=omega/c
    p1=c_1_4pi*exp(-c_im*k*1.d0)/1.d0
    rv=x-x0
    r=sqrt(dot_product(rv,rv))
    drdx=rv/r
    drdn=dot_product(drdx,n)
    if (space.eq.2) stop 'not yet'
    select case (rn)
      case (2)
        stop 'not yet'
      case (3)
        p_inc=A/p1*c_1_4pi*exp(-c_im*k*r)/r
        Un_inc=-A/p1/(rho*omega**2)*c_1_4pi*(1.d0/r+c_im*k)/r*drdn*exp(-c_im*k*r)
    end select
  end subroutine fbem_harpot_pointwave

  !! Full-space or half-space incident plane wave for inviscid fluids. For the half-space solution, the plane can have the
  !! orientation of any axis, any coordinate, and any boundary condition (p=0 or Un=0).
  subroutine fbem_harpot_planewave(rn,omega,rho,c,A,x0,varphi,theta,space,np,xp,bc,symconf,xs,x,n,p_inc,Un_inc)
    implicit none
    ! Input/Output
    integer              :: rn            !! Dimensional space: 2 or 3
    real(kind=real64)    :: omega         !! Angular frequency
    real(kind=real64)    :: rho           !! Density
    complex(kind=real64) :: c             !! Wave propagation speed
    complex(kind=real64) :: A             !! Amplitude of the incident wave
    real(kind=real64)    :: x0(rn)        !! Origin of incident plane wave
    real(kind=real64)    :: varphi        !! Angle of incidence
    real(kind=real64)    :: theta         !! Angle of incidence
    integer              :: space         !! Type of space: 1 full-space, 2 half-space.
    integer              :: np            !! Axis corresponding to the normal of the plane (x_np=xp)
    real(kind=real64)    :: xp            !! Position of the plane
    integer              :: bc            !! Boundary condition at the plane: 0 p=0, 1 Un=0.
    real(kind=real64)    :: symconf(rn)   !! Symmetry configuration vector. It must be meet that symconf(np)=0.
    real(kind=real64)    :: xs(rn)        !! Origin of symmetry/antisymmetry decomposition
    real(kind=real64)    :: x(rn)         !! Observation point
    real(kind=real64)    :: n(rn)         !! Unit normal of the observation point
    complex(kind=real64) :: p_inc         !! Pressure of the incident wave
    complex(kind=real64) :: Un_inc        !! Normal displacement of the incident wave
    ! Local
    complex(kind=real64) :: k           ! Wavenumber
    integer              :: ki, kj      ! Counters
    real(kind=real64)    :: q_dot_x0xs  ! Dot product q·(x0-xs)
    real(kind=real64)    :: q(rn)       ! Propagation direction vector
    complex(kind=real64) :: AC, AE      ! Effective amplitudes
    complex(kind=real64) :: half_exp_m  ! Exponential term -
    complex(kind=real64) :: half_exp_p  ! Exponential term +
    complex(kind=real64) :: Ws, Wa      ! Symmetric and antisymmetric part of an W
    complex(kind=real64) :: W(rn)       ! W for each coordinate
    complex(kind=real64) :: Zs, Za      ! Symmetric and antisymmetric part of an Z
    complex(kind=real64) :: Z(rn)       ! Z for each coordinate
    complex(kind=real64) :: WZqn        ! W*Z*q*n products
    complex(kind=real64) :: Ar          ! Amplitude of the reflected wave
    complex(kind=real64) :: p_ref       ! Pressure of the reflected wave
    complex(kind=real64) :: Un_ref      ! Normal displacement of the reflected wave
    ! Wavenumber
    k=omega/c
    !
    ! Incident wave
    !
    ! Propagation direction vector
    select case (rn)
      case (2)
        q(1)=cos(theta)
        q(2)=sin(theta)
      case (3)
        q(1)=cos(theta)*sin(varphi)
        q(2)=cos(theta)*cos(varphi)
        q(3)=sin(theta)
    end select
    ! q·(x0-xs)
    q_dot_x0xs=0.0d0
    do ki=1,rn
      q_dot_x0xs=q_dot_x0xs+q(ki)*(x0(ki)-xs(ki))
    end do
    ! Components effective amplitudes
    AC=A*exp(c_im*k*q_dot_x0xs)
    AE=-c_im*k/(rho*omega**2)*AC
    ! Calculate W_j^s, W_j^a, Z_j^s, Z_j^a, and the final W_j and Z_j considering the symmetry configuration
    do ki=1,rn
      half_exp_m=0.5d0*exp(-c_im*k*q(ki)*(x(ki)-xs(ki)))
      half_exp_p=0.5d0*exp( c_im*k*q(ki)*(x(ki)-xs(ki)))
      Ws=half_exp_m+half_exp_p
      Wa=half_exp_m-half_exp_p
      Zs=Wa
      Za=Ws
      ! Depending on the symmetry configuration
      select case (nint(symconf(ki)))
        case (-1)
          W(ki)=Wa
          Z(ki)=Za
        case ( 0)
          W(ki)=Ws+Wa
          Z(ki)=Zs+Za
        case ( 1)
          W(ki)=Ws
          Z(ki)=Zs
      end select
    end do
    ! Calculate incident p and Un
    p_inc=AC
    Un_inc=(0.0d0,0.0d0)
    do ki=1,rn
      p_inc=p_inc*w(ki)
      WZqn=q(ki)*n(ki)
      do kj=1,rn
        if (kj.eq.ki) then
          WZqn=WZqn*Z(kj)
        else
          WZqn=WZqn*W(kj)
        end if
      end do
      Un_inc=Un_inc+WZqn
    end do
    Un_inc=Un_inc*AE
    !
    ! Reflected wave for the half-space solution
    !
    if (space.eq.2) then
      ! Check restriction
      if (nint(symconf(np)).ne.0) then
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                'The symmetry configuration is not valid because the half-space plane cannot be a symmetry plane')
      end if
      ! Propagation direction vector of the reflected wave
      q(np)=-q(np)
      ! q·(x0-xs)
      q_dot_x0xs=0.0d0
      do ki=1,rn
        q_dot_x0xs=q_dot_x0xs+q(ki)*(x0(ki)-xs(ki))
      end do
      ! Amplitude of the reflected wave
      select case (bc)
        ! p=0 at x_np=xp
        case (0)
          Ar=-A*exp(-2.0d0*c_im*k*q(np)*(xp-x0(np)))
        ! Un=0 at x_np=xp
        case (1)
          Ar=A*exp(-2.0d0*c_im*k*q(np)*(xp-x0(np)))
      end select
      ! Components effective amplitudes
      AC=Ar*exp(c_im*k*q_dot_x0xs)
      AE=-c_im*k/(rho*omega**2)*AC
      ! Calculate W_j^s, W_j^a, Z_j^s, Z_j^a, and the final W_j and Z_j considering the symmetry configuration
      do ki=1,rn
        half_exp_m=0.5d0*exp(-c_im*k*q(ki)*(x(ki)-xs(ki)))
        half_exp_p=0.5d0*exp( c_im*k*q(ki)*(x(ki)-xs(ki)))
        Ws=half_exp_m+half_exp_p
        Wa=half_exp_m-half_exp_p
        Zs=Wa
        Za=Ws
        ! Depending on the symmetry configuration
        select case (nint(symconf(ki)))
          case (-1)
            W(ki)=Wa
            Z(ki)=Za
          case ( 0)
            W(ki)=Ws+Wa
            Z(ki)=Zs+Za
          case ( 1)
            W(ki)=Ws
            Z(ki)=Zs
        end select
      end do
      ! Calculate incident p and Un
      p_ref=AC
      Un_ref=(0.0d0,0.0d0)
      do ki=1,rn
        p_ref=p_ref*w(ki)
        WZqn=q(ki)*n(ki)
        do kj=1,rn
          if (kj.eq.ki) then
            WZqn=WZqn*Z(kj)
          else
            WZqn=WZqn*W(kj)
          end if
        end do
        Un_ref=Un_ref+WZqn
      end do
      Un_ref=Un_ref*AE
      !
      ! Add incident and reflected waves
      !
      p_inc=p_inc+p_ref
      Un_inc=Un_inc+Un_ref
    end if
  end subroutine fbem_harpot_planewave

  !! Calculate all properties by giving two of them
  subroutine fbem_harpot_properties(pK,K,prho,rho,pc,c)
    implicit none
    ! I/O
    logical           :: pK
    real(kind=real64) :: K
    logical           :: prho
    real(kind=real64) :: rho
    logical           :: pc
    real(kind=real64) :: c
    ! Local
    integer :: n
    ! Check input constants
    n=0
    if (pK  ) n=n+1
    if (prho) n=n+1
    if (pc  ) n=n+1
    if (n.ne.2) stop 'only 2 properties are needed'

    if (pK.and.prho) then
      c=sqrt(K/rho)
    end if
    if (pK.and.pc) then
      rho=K/c**2
    end if
    if (prho.and.pc) then
      K=rho*c**2
    end if
    ! Check if constants are valid
    if (  K.le.0.d0) stop 'K must be >0'
    if (  c.le.0.d0) stop 'c must be >0'
    if (rho.le.0.d0) stop 'rho must be >0'
  end subroutine fbem_harpot_properties

end module fbem_harpot_incident_field
