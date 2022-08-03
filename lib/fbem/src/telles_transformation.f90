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
!! <b> This module implements the Telles transformation:
!!     - J.C.F. Telles. A self-adaptive co-ordinate transformation for efficient numerical evaluation of general boundary element
!!       integrals. International Journal for Numerical Methods in Engineering, 24:959-973, 1987.
!!     - J.C.F. Telles and R.F. Oliveira. Third degree polynomial transformation for boundary element integrals: further
!!       improvements. Engineering Analysis with Boundary Elements, 13:135-141, 1994.
!!</b>
!!
!! <h2>DESCRIPTION</h2>
!!
!! The Telles transformation is a third degree polynomial transformation that regularizes the integrand of quasi-singular
!! integrals. The key point of the Telles transformation is finding the best value of the jacobian at the nearest point of the
!! element with respect to the singular point. Telles found that is enough to find a relationship between the jacobian at nearest
!! point and the normalized distance for each type of quasi-singular integrand. He studied optimum relationships for
!! <tt>f={ln(r),1/r,1/r^2,1/r^3}</tt> and his results are very useful.
!!
!! However, the numerical experiments carried out by him didn't take into account all possible locations of singular points with
!! respect to the element. So, new numerical experiments exploring singular points at all possible locations have been designed and
!! executed by me. New jacobian vs distance functions have been obtained for <tt>f={ln(r),1/r,1/r^2,1/r^3}</tt>. Because all the
!! relationships follow the same patterns for all integrands, a mean relationship has been built too.
module fbem_telles_transformation

  ! Fortran 2003 standard
  use iso_fortran_env

  ! fbem module common
  use fbem_numerical
  use fbem_string_handling

  ! No implicit variables are allowed
  implicit none

  ! By default all are private
  private

  ! Public type and functions
  public :: fbem_telles_parameters
  public :: fbem_telles_barr
  public :: fbem_telles11_calculate_parameters
  public :: fbem_telles01_calculate_parameters
  public :: fbem_telles_xi
  public :: fbem_telles_jacobian
  public :: fbem_telles_xi_and_jacobian

  !! Telles parameters data structure (a,b,c,d).
  type fbem_telles_parameters
    real(kind=real64) :: c(4) !! <tt>c(1)=a, c(2)=b, c(3)=c, c(4)=d</tt>
  end type fbem_telles_parameters

contains

  !! Obtained optimum Telles jacobian at nearest point of the element with respect the singular point
  function fbem_telles_barr(d, f)
    !! Estimated optimum value of the jacobian at the nearest point .
    real(kind=real64) :: fbem_telles_barr
    !! Normalized distance between nearest point and collocation point.
    real(kind=real64) :: d
    !! Integrand type: <tt>f={fbem_f_lnr,fbem_f_1r1,fbem_f_1r2,fbem_f_1r3}</tt>.
    integer :: f
    if (d.lt.0.0d0) then
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,'d<0 is invalid.')
    end if
    select case (f)
      case (fbem_f_any)
        ! Mean relationship: lnr, 1/r, 1/r^2 and 1/r^3
        if (d.lt.3.0d0) then
          fbem_telles_barr=d/(0.89039d0*d+0.32883d0)
        else
          fbem_telles_barr=1.0d0
        end if
      case (fbem_f_lnr)
        ! For lnr
        if (d.lt.0.5d0) then
          fbem_telles_barr=d/(0.964859d0*d+0.017571d0)
        else
          fbem_telles_barr=1.0d0
        end if
      case (fbem_f_1r1)
        ! For 1/r
        if (d.lt.2.0d0) then
          fbem_telles_barr=d/(0.683416d0*d+0.633168d0)
        else
          fbem_telles_barr=1.0d0
        end if
      case (fbem_f_1r2)
        ! For 1/r^2
        if (d.lt.3.0d0) then
          fbem_telles_barr=d/(0.864422d0*d+0.406734d0)
        else
          fbem_telles_barr=1.0d0
        end if
      case (fbem_f_1r3)
        ! For 1/r^3
        if (d.lt.4.0d0) then
          fbem_telles_barr=d/(0.88073d0*d+0.477080d0)
        else
          fbem_telles_barr=1.0d0
        end if
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'given integrand is invalid')
    end select
    if (fbem_telles_barr.gt.1.0d0) fbem_telles_barr=1.0d0
  end function fbem_telles_barr

  !! Calculation of parameters for Telles transformation in [-1,1] domain
  function fbem_telles11_calculate_parameters(bar_xi, bar_r)
    !! Telles transformation in [-1,1] parameters (a, b, c, d)
    type(fbem_telles_parameters) :: fbem_telles11_calculate_parameters
    !! Nearest point reference coordinate
    real(kind=real64) :: bar_xi
    !! Jacobian value at the nearest point
    real(kind=real64) :: bar_r
    real(kind=real64) :: w, p, q, R2, R31, R32, bar_gamma, Qaux, d13
    w=bar_xi/(1.0d0+2.0d0*bar_r)
    p=1.0d0/(3.0d0*(1.0d0+2.0d0*bar_r))*(3.0d0-2.0d0*bar_r-3.0d0*w*bar_xi)
    q=w/2.0d0*((3.0d0-2.0d0*bar_r)/(1.0d0+2.0d0*bar_r)-2.0d0*w**2-1.0d0)
    R2=dsqrt(q**2+p**3)
    R31=-q+R2
    R32=-q-R2
    d13=1.0d0/3.0d0
    if (R31.ge.0.0d0) then
      R31=R31**d13
    else
      R31=-(dabs(R31))**d13
    end if
    if (R32.ge.0.0d0) then
      R32=R32**d13
    else
      R32=-(dabs(R32))**d13
    end if
    bar_gamma=R31+R32+w
    Qaux=1.0d0+3.0d0*bar_gamma**2
    fbem_telles11_calculate_parameters%c(1)=(1.0d0-bar_r)/Qaux
    fbem_telles11_calculate_parameters%c(2)=-3.0d0*bar_gamma*fbem_telles11_calculate_parameters%c(1)
    fbem_telles11_calculate_parameters%c(3)=(bar_r+3.0d0*bar_gamma**2)/Qaux
    fbem_telles11_calculate_parameters%c(4)=-fbem_telles11_calculate_parameters%c(2)
  end function fbem_telles11_calculate_parameters

  !! Calculation of parameters for Telles transformation in [0,1] domain
  function fbem_telles01_calculate_parameters(bar_xi, bar_r)
    !! Telles transformation in [-1,1] parameters (a, b, c, d)
    type(fbem_telles_parameters) :: fbem_telles01_calculate_parameters
    !! Nearest point reference coordinate
    real(kind=real64) :: bar_xi
    !! Jacobian value at the nearest point
    real(kind=real64) :: bar_r
    real(kind=real64) :: w, p, q, R2, R31, R32, bar_gamma, Qaux, waux, d13
    waux=1.0d0+2.0d0*bar_r
    w=(bar_xi+bar_r)/waux
    p=(3.0d0*bar_xi+bar_r)/(3.0d0*waux)-w**2
    q=w*((3.0d0*bar_xi+bar_r)/(2.0d0*waux)-w**2)-bar_xi/(2.0d0*waux)
    R2=dsqrt(q**2+p**3)
    R31=-q+R2
    R32=-q-R2
    d13=1.0d0/3.0d0
    if (R31.ge.0.0d0) then
      R31=R31**d13
    else
      R31=-(dabs(R31))**d13
    end if
    if (R32.ge.0.0d0) then
      R32=R32**d13
    else
      R32=-(dabs(R32))**d13
    end if
    bar_gamma=R31+R32+w
    Qaux=3.0d0*bar_gamma*(bar_gamma-1.0d0)+1.0d0
    fbem_telles01_calculate_parameters%c(1)=(1.0d0-bar_r)/Qaux
    fbem_telles01_calculate_parameters%c(2)=-3.0d0*bar_gamma*fbem_telles01_calculate_parameters%c(1)
    fbem_telles01_calculate_parameters%c(3)=(3.0d0*bar_gamma*(bar_gamma-bar_r)+bar_r)/Qaux
    fbem_telles01_calculate_parameters%c(4)=0.d0
  end function fbem_telles01_calculate_parameters

  !! Calculation of the xi coordinate <tt>xi(gamma)</tt>
  function fbem_telles_xi(params, vgamma)
    !! xi coordinate
    real(kind=real64) :: fbem_telles_xi
    !! Telles transformation parameters
    type(fbem_telles_parameters) :: params
    !! gamma coordinate
    real(kind=real64) :: vgamma
    fbem_telles_xi=params%c(1)*vgamma**3+params%c(2)*vgamma**2+params%c(3)*vgamma+params%c(4)
  end function fbem_telles_xi

  !! Calculation of the jacobian of the transformation <tt>xi(gamma)</tt>
  function fbem_telles_jacobian(params, vgamma)
    !! jacobian
    real(kind=real64) :: fbem_telles_jacobian
    !! Telles transformation parameters
    type(fbem_telles_parameters) :: params
    !! gamma coordinate
    real(kind=real64) :: vgamma
    fbem_telles_jacobian=3.0d0*params%c(1)*vgamma**2+2.0d0*params%c(2)*vgamma+params%c(3)
  end function fbem_telles_jacobian

  !! Calculation of the xi coordinate and jacobian
  subroutine fbem_telles_xi_and_jacobian(params, vgamma, xi, jacobian)
    !! Telles transformation parameters
    type(fbem_telles_parameters) :: params
    !! gamma coordinate
    real(kind=real64) :: vgamma
    !! xi coordinate
    real(kind=real64) :: xi
    !! jacobian
    real(kind=real64) :: jacobian
    xi=params%c(1)*vgamma**3+params%c(2)*vgamma**2+params%c(3)*vgamma+params%c(4)
    jacobian=3.0d0*params%c(1)*vgamma**2+2.0d0*params%c(2)*vgamma+params%c(3)
  end subroutine

end module fbem_telles_transformation
