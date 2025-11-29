! ---------------------------------------------------------------------
! Copyright (C) 2014-2025 Universidad de Las Palmas de Gran Canaria:
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
!! @version 1.0
!!
!! <b> This module has sound/acoustics related function and subroutines.</b>
!!
!! See: Maeso and Aznarez, Estrategias para la reduccion del impacto acustico en el entorno de carreteras. Una aplicacion
!!      del Metodo de Elementos de Contorno, 2005. http://hdl.handle.net/10553/1500
!!
!! TODO:

module fbem_acoustics

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_string_handling

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all are private
  private

  ! Public
  public :: fbem_octave_f
  public :: fbem_onethird_octave_f
  public :: fbem_pressure_spl_db
  public :: fbem_pressure_dBA_correction
  public :: fbem_IL_i
  public :: fbem_SIL
  public :: fbem_beta_st
  public :: fbem_delanybazley1970_thick_beta_n
  public :: fbem_delanybazley1970_thick_beta_st
  public :: fbem_delanybazley1970_thin_beta_n
  public :: fbem_delanybazley1970_thin_beta_st

  real(kind=real64), parameter, dimension(10) :: fbem_octave_f=(/31.5d0,&
                                                                 63.0d0,&
                                                                  125.0d0,&
                                                                  250.0d0,&
                                                                  500.0d0,&
                                                                 1000.0d0,&
                                                                 2000.0d0,&
                                                                 4000.0d0,&
                                                                 8000.0d0,&
                                                                 16000.0d0/)
  real(kind=real64), parameter, dimension(30) :: fbem_onethird_octave_f=(/25.0d0,   31.5d0,   40.0d0,&
                                                                          50.0d0,   63.0d0,   80.0d0,&
                                                                         100.0d0,  125.0d0,  160.0d0,&
                                                                         200.0d0,  250.0d0,  315.0d0,&
                                                                         400.0d0,  500.0d0,  630.0d0,&
                                                                         800.0d0, 1000.0d0, 1250.0d0,&
                                                                        1600.0d0, 2000.0d0, 2500.0d0,&
                                                                        3150.0d0, 4000.0d0, 5000.0d0,&
                                                                        6300.0d0, 8000.0d0,10000.0d0,&
                                                                       12500.0d0,16000.0d0,20000.0d0/)

contains

  !! Calculate pressure in SPL (Sound Pressure Level) by giving pressure in Pa.
  function fbem_pressure_spl_db(p)
    implicit none
    !! Pressure in Pa
    real(kind=real64) :: p
    !! Pressure in SPL
    real(kind=real64) :: fbem_pressure_spl_db
    fbem_pressure_spl_db=10.0d0*dlog10(0.5d0*(p**2))+94.0d0
  end function fbem_pressure_spl_db

  !! A-weighting correction in dB in order to calculate in pressure in dBA.
  !! ANSI Standards S1.4-1983 and S1.42-2001
  function fbem_pressure_dBA_correction(f)
    implicit none
    !! Frequency, it must be between 16 and 20000 Hz
    real(kind=real64) :: f
    !! Correction in dB
    real(kind=real64) :: fbem_pressure_dBA_correction
    ! Correction factor R_A
    real(kind=real64) :: ra
    ! Check i
    if ((f.lt.16.0d0).or.(f.gt.20000.0d0)) then
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                              'frequency is out of range to calculate A-weighted dB correction.')
    end if
    ra=(12200.0d0**2)*f**4/(f**2+20.6d0**2)/(f**2+12200.0d0**2)/dsqrt((f**2+107.7d0**2)*(f**2+737.9d0**2))
    fbem_pressure_dBA_correction=2.0d0+20.0d0*dlog10(ra)
  end function fbem_pressure_dBA_correction

  !! Insertion Loss coefficient
  function fbem_IL_i(p_before,p_after)
    implicit none
    !! Before correction measure pressure in Pa
    real(kind=real64) :: p_before
    !! After correction measure pressure in Pa
    real(kind=real64) :: p_after
    !! Insertion Loss
    real(kind=real64) :: fbem_IL_i
    fbem_IL_i=-20.0d0*dlog10(p_after/p_before)
  end function fbem_IL_i

  !! Mean spectral Insertion Loss coeficient
  function fbem_SIL(n,IL)
    implicit none
    !! Spectral IL
    real(kind=real64) :: fbem_SIL
    !! Number of frequencies
    integer :: n
    !! Before correction measure pressure in Pa
    real(kind=real64) :: IL(n)
    ! Auxiliary variables
    real(kind=real64) :: sump
    ! Counter
    integer :: i
    ! Calculate sum
    sump=0.0d0
    do i=1,n
      sump=sump+10.0d0**(-IL(i)/10.0d0)
    end do
    ! Spectral IL
    fbem_SIL=10.0d0*dlog10(sump)
  end function fbem_SIL

  !! Delany & Bazley (1970) model of absorbing materials.
  !! M.E. Delany, E.N. Bazley, Acoustical properties of fibrous absorbent materials, Applied Acoustics, Volume 3, Issue 2, 1970.
  !! See also: Sagartzazu, "Review in Sound Absorbing Materials"
  !! TODO: Review the range of validity of this model
  function fbem_delanybazley1970_thick_beta_n(sigma,f)
    implicit none
    ! I/O
    complex(kind=real64) :: fbem_delanybazley1970_thick_beta_n !! Admittance according to Delany & Bazley (1970)
    real(kind=real64)    :: sigma                              !! Absorbing material flow resistivity (N·s/m^4)
    real(kind=real64)    :: f                                  !! Frequency (Hz)
    ! Local
    real(kind=real64)    :: fs
    complex(kind=real64) :: beta
    fs=1000.d0*f/sigma
    beta=1.d0/(1.d0+9.08d0*fs**(-0.75d0)-c_im*11.9d0*fs**(-0.73d0))
    fbem_delanybazley1970_thick_beta_n = beta
  end function fbem_delanybazley1970_thick_beta_n

  !! Delany & Bazley (1970) model of absorbing materials (thin layer).
  !! M.E. Delany, E.N. Bazley, Acoustical properties of fibrous absorbent materials, Applied Acoustics, Volume 3, Issue 2, 1970.
  !! See also: Sagartzazu, "Review in Sound Absorbing Materials"
  !! TODO: Review the range of validity of this model
  function fbem_delanybazley1970_thin_beta_n(rho,c,sigma,t,f)
    implicit none
    ! I/O
    complex(kind=real64) :: fbem_delanybazley1970_thin_beta_n !! Admittance according to Delany & Bazley (1970)
    real(kind=real64)    :: rho                               !! Fluid density (kg/m^3)
    real(kind=real64)    :: c                                 !! Fluid speed propagation (m/s)
    real(kind=real64)    :: sigma                             !! Absorbing material flow resistivity (N·s/m^4)
    real(kind=real64)    :: t                                 !! Thickness (m)
    real(kind=real64)    :: f                                 !! Frequency (Hz)
    ! Local
    real(kind=real64)    :: fs
    real(kind=real64)    :: k
    complex(kind=real64) :: k_layer
    complex(kind=real64) :: beta
    complex(kind=real64) :: z
    fs=1000.d0*f/sigma
    k = 2.d0*c_pi*f/c
    k_layer = k*(10.3d0*fs**(-0.59d0)+c_im*(1.d0+10.8d0*fs**(-0.7d0)))
    beta = 1.d0/(1.d0+9.08d0*fs**(-0.75d0)-c_im*11.9d0*fs**(-0.73d0))
    z = k_layer*t
    !fbem_delanybazley1970_thin_beta_n = beta*(cdexp(z)+cdexp(-z))/(cdexp(z)-cdexp(-z))
    fbem_delanybazley1970_thin_beta_n = beta*tanh(k_layer*t)
  end function fbem_delanybazley1970_thin_beta_n

  !! Admittance obtained from statistical absorption coefficient, from admittance (normal incidence)
  !! Morse & Bolt 1944, Bies & Hansen (A4.6).
  function fbem_beta_st(beta)
    implicit none
    ! I/O
    real(kind=real64)    :: fbem_beta_st !! Admittance from statistical absorption coefficient
    complex(kind=real64) :: beta         !! Admittance (normal incidence)
    ! Local
    complex(kind=real64) :: z
    real(kind=real64)    :: xi, theta
    real(kind=real64)    :: alpha_st
    z = 1.d0/beta
    xi = cdabs(z)
    theta = datan(dimag(z)/dreal(z))
    alpha_st = 1.d0-dcos(theta)/xi*dlog(1.d0+2.d0*xi*dcos(theta)+xi**2)&
             + dcos(2.d0*theta)/xi/dsin(theta)*datan(xi*dsin(theta)/(1.d0+xi*dcos(theta)))
    alpha_st = 8.d0*dcos(theta)/xi*alpha_st
    fbem_beta_st = (1.d0-sqrt(1.d0-alpha_st))/(1.d0+sqrt(1.d0-alpha_st))
  end function fbem_beta_st

  !! Delany & Bazley (1970) model of absorbing materials. Admittance (statistical).
  !! M.E. Delany, E.N. Bazley, Acoustical properties of fibrous absorbent materials, Applied Acoustics, Volume 3, Issue 2, 1970.
  !! See also: Sagartzazu, "Review in Sound Absorbing Materials"
  !! TODO: Review the range of validity of this model
  function fbem_delanybazley1970_thick_beta_st(sigma,f)
    implicit none
    ! I/O
    real(kind=real64) :: fbem_delanybazley1970_thick_beta_st !! Admittance according to Delany & Bazley (1970)
    real(kind=real64) :: sigma                               !! Absorbing material flow resistivity (N·s/m^4)
    real(kind=real64) :: f                                   !! Frequency (Hz)
    ! Local
    complex(kind=real64) :: beta_n
    beta_n = fbem_delanybazley1970_thick_beta_n(sigma,f)
    fbem_delanybazley1970_thick_beta_st = fbem_beta_st(beta_n)
  end function fbem_delanybazley1970_thick_beta_st

  !! Delany & Bazley (1970) model of absorbing materials (thin layer). Admittance (statistical).
  !! M.E. Delany, E.N. Bazley, Acoustical properties of fibrous absorbent materials, Applied Acoustics, Volume 3, Issue 2, 1970.
  !! See also: Sagartzazu, "Review in Sound Absorbing Materials"
  !! TODO: Review the range of validity of this model
  function fbem_delanybazley1970_thin_beta_st(rho,c,sigma,t,f)
    implicit none
    ! I/O
    complex(kind=real64) :: fbem_delanybazley1970_thin_beta_st !! Admittance according to Delany & Bazley (1970)
    real(kind=real64)    :: rho                                !! Fluid density (kg/m^3)
    real(kind=real64)    :: c                                  !! Fluid speed propagation (m/s)
    real(kind=real64)    :: sigma                              !! Absorbing material flow resistivity (N·s/m^4)
    real(kind=real64)    :: t                                  !! Thickness (m)
    real(kind=real64)    :: f                                  !! Frequency (Hz)
    ! Local
    complex(kind=real64) :: beta_n
    beta_n = fbem_delanybazley1970_thin_beta_n(rho,c,sigma,t,f)
    fbem_delanybazley1970_thin_beta_st = fbem_beta_st(beta_n)
  end function fbem_delanybazley1970_thin_beta_st

end module fbem_acoustics
