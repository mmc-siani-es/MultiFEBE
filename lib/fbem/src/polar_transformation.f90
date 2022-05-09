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
!! <b> This module implements the polar transformation for 3D BEM weakly singular integrals. </b>
module fbem_polar_transformation

  ! Fortran 2003 standard
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_string_handling
  use fbem_shape_functions

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all are private
  private

  public :: fbem_polar_transformation_setup
  public :: fbem_polar_transformation_angular
  public :: fbem_degsubquad_transformation_setup

contains

  !! This subroutine do the setup previous to the polar transformation
  subroutine fbem_degsubquad_transformation_setup(type_g,xi_i,ndsq,dsq,xi_dsq)
    implicit none
    ! I/O
    integer           :: type_g             !! Geometrical interpolation
    real(kind=real64) :: xi_i(2)            !! Reference coordinates of the singular point.
    integer           :: ndsq               !! Number of degenerated quadrilateral subdivisions
    integer           :: dsq(8)             !! Vector that contains what degenerated quadrilaterals need to be integrated
    real(kind=real64) :: xi_dsq(2,4,8)      !! Matrix that contains the angles theta of the subtriangles to be integrated
    ! Local
    logical           :: i_in_edge          ! Location of collocation point
    integer           :: kdsq
    !
    ! Check if xi is on the border of the element or not, depending on this the integration scheme is different.
    !
    i_in_edge=fbem_check_xi1xi2_edge(type_g,xi_i)
    !
    ! Selection of degenerated quadrilaterals that need to be integrated
    !
    ! If the collocation point is not on an edge, all have to be integrated.
    if (i_in_edge.eqv.(.false.)) then
      select case (fbem_n_edges(type_g))
        case (3)
          ndsq=6
        case (4)
          ndsq=8
      end select
      do kdsq=1,ndsq
        dsq(kdsq)=kdsq
      end do
    ! If the collocation point is on an edge, only some of them need to be integrated.
    else
      ndsq=0
      select case (fbem_n_edges(type_g))
        ! Quadrilateral
        case (4)
          ! If collocation point is on node 1
          if ((xi_i(1).le.-1.0d0+check_xi_tolerance).and.(xi_i(2).le.-1.0d0+check_xi_tolerance)) then
            ndsq=2
            dsq(1)=4
            dsq(2)=5
          end if
          ! If collocation point is on node 2
          if ((xi_i(1).ge. 1.0d0-check_xi_tolerance).and.(xi_i(2).le.-1.0d0+check_xi_tolerance)) then
            ndsq=2
            dsq(1)=6
            dsq(2)=7
          end if
          ! If collocation point is on node 3
          if ((xi_i(1).ge. 1.0d0-check_xi_tolerance).and.(xi_i(2).ge. 1.0d0-check_xi_tolerance)) then
            ndsq=2
            dsq(1)=8
            dsq(2)=1
          end if
          ! If collocation point is on node 4
          if ((xi_i(1).le.-1.0d0+check_xi_tolerance).and.(xi_i(2).ge. 1.0d0-check_xi_tolerance)) then
            ndsq=2
            dsq(1)=2
            dsq(2)=3
          end if
          ! If it is not on a vertex node
          if (ndsq.eq.0) then
            ! If collocation point is on edge 1
            if (xi_i(2).le.-1.0d0+check_xi_tolerance) then
              ndsq=4
              dsq(1)=4
              dsq(2)=5
              dsq(3)=6
              dsq(4)=7
            end if
            ! If collocation point is on edge 2
            if (xi_i(1).ge. 1.0d0-check_xi_tolerance) then
              ndsq=4
              dsq(1)=6
              dsq(2)=7
              dsq(3)=8
              dsq(4)=1
            end if
            ! If collocation point is on edge 3
            if (xi_i(2).ge. 1.0d0-check_xi_tolerance) then
              ndsq=4
              dsq(1)=8
              dsq(2)=1
              dsq(3)=2
              dsq(4)=3
            end if
            ! If collocation point is on edge 4
            if (xi_i(1).le.-1.0d0+check_xi_tolerance) then
              ndsq=4
              dsq(1)=2
              dsq(2)=3
              dsq(3)=4
              dsq(4)=5
            end if
          end if
        ! Triangular
        case (3)
          ! If collocation point is on node 1
          if (xi_i(1).ge.1.0d0-check_xi_tolerance) then
            ndsq=1
            dsq(1)=3
          end if
          ! If collocation point is on node 2
          if (xi_i(2).ge.1.0d0-check_xi_tolerance) then
            ndsq=1
            dsq(1)=6
          end if
          ! If collocation point is on node 3
          if ((xi_i(1).le.check_xi_tolerance).and.(xi_i(2).le.check_xi_tolerance)) then
            ndsq=2
            dsq(1)=1
            dsq(2)=2
          end if
          ! If it is not on a vertex node
          if (ndsq.eq.0) then
            ! If collocation point is on edge 1
            if ((xi_i(1)+xi_i(2)).ge.1.0d0-check_xi_tolerance) then
              ndsq=4
              dsq(1)=3
              dsq(2)=4
              dsq(3)=5
              dsq(4)=6
            end if
            ! If collocation point is on edge 2
            if (xi_i(1).le.check_xi_tolerance) then
              ndsq=3
              dsq(1)=6
              dsq(2)=1
              dsq(3)=2
            end if
            ! If collocation point is on edge 3
            if (xi_i(2).le.check_xi_tolerance) then
              ndsq=3
              dsq(1)=1
              dsq(2)=2
              dsq(3)=3
            end if
          end if
      end select
    end if
    !
    ! Coordinates of each degenerated quadrilateral
    !
    do kdsq=1,ndsq
      xi_dsq(:,1,kdsq)=xi_i
      xi_dsq(:,2,kdsq)=xi_dsq(:,1,kdsq)
      select case (fbem_n_edges(type_g))
        ! Quadrilateral
        case (4)
          select case (dsq(kdsq))
            case (1)
              xi_dsq(:,3,kdsq)=[  -1.d0,  -1.d0]
              xi_dsq(:,4,kdsq)=[xi_i(1),  -1.d0]
            case (2)
              xi_dsq(:,3,kdsq)=[xi_i(1),  -1.d0]
              xi_dsq(:,4,kdsq)=[   1.d0,  -1.d0]
            case (3)
              xi_dsq(:,3,kdsq)=[   1.d0,  -1.d0]
              xi_dsq(:,4,kdsq)=[   1.d0,xi_i(2)]
            case (4)
              xi_dsq(:,3,kdsq)=[   1.d0,xi_i(2)]
              xi_dsq(:,4,kdsq)=[   1.d0,   1.d0]
            case (5)
              xi_dsq(:,3,kdsq)=[   1.d0,   1.d0]
              xi_dsq(:,4,kdsq)=[xi_i(1),   1.d0]
            case (6)
              xi_dsq(:,3,kdsq)=[xi_i(1),   1.d0]
              xi_dsq(:,4,kdsq)=[  -1.d0,   1.d0]
            case (7)
              xi_dsq(:,3,kdsq)=[  -1.d0,   1.d0]
              xi_dsq(:,4,kdsq)=[  -1.d0,xi_i(2)]
            case (8)
              xi_dsq(:,3,kdsq)=[  -1.d0,xi_i(2)]
              xi_dsq(:,4,kdsq)=[  -1.d0,  -1.d0]
          end select
        ! Triangular
        case (3)
          select case (dsq(kdsq))
            case (1)
              xi_dsq(:,3,kdsq)=[   1.d0,   0.d0]
              xi_dsq(:,4,kdsq)=0.5d0*[1.d0-xi_i(2)+xi_i(1),1.d0+xi_i(2)-xi_i(1)]
            case (2)
              xi_dsq(:,3,kdsq)=0.5d0*[1.d0-xi_i(2)+xi_i(1),1.d0+xi_i(2)-xi_i(1)]
              xi_dsq(:,4,kdsq)=[   0.d0,   1.d0]
            case (3)
              xi_dsq(:,3,kdsq)=[   0.d0,   1.d0]
              xi_dsq(:,4,kdsq)=[   0.d0,xi_i(2)]
            case (4)
              xi_dsq(:,3,kdsq)=[   0.d0,xi_i(2)]
              xi_dsq(:,4,kdsq)=[   0.d0,   0.d0]
            case (5)
              xi_dsq(:,3,kdsq)=[   0.d0,   0.d0]
              xi_dsq(:,4,kdsq)=[xi_i(1),   0.d0]
            case (6)
              xi_dsq(:,3,kdsq)=[xi_i(1),   0.d0]
              xi_dsq(:,4,kdsq)=[   1.d0,   0.d0]
          end select
      end select
    end do
  end subroutine fbem_degsubquad_transformation_setup

  !! This subroutine do the setup previous to the polar transformation
  subroutine fbem_polar_transformation_setup(type_g,xi_i,nsubtri,subtriangle,theta_subtri,thetap_subtri)
    implicit none
    ! I/O
    integer           :: type_g             !! Geometrical interpolation
    real(kind=real64) :: xi_i(2)            !! Reference coordinates of the singular point.
    integer           :: nsubtri            !! Number of subtriangles
    integer           :: subtriangle(8)     !! Vector that contains what subtriangles need to be integrated
    real(kind=real64) :: theta_subtri(2,8)  !! Matrix that contains the angles theta of the subtriangles to be integrated
    real(kind=real64) :: thetap_subtri(2,8) !! Matrix that contains the angles thetap of the subtriangles to be integrated
    ! Local
    logical           :: i_in_edge          ! Location of collocation point
    integer           :: ksubtri            ! Counter variable for subtriangles loop
    real(kind=real64) :: theta1             ! Angle of the vertex 1
    real(kind=real64) :: theta2             ! Angle of the vertex 2
    real(kind=real64) :: theta3             ! Angle of the vertex 3
    real(kind=real64) :: theta4             ! Angle of the vertex 4
    !
    ! Check if xi is on the border of the element or not, depending on this the integration scheme is different.
    !
    i_in_edge=fbem_check_xi1xi2_edge(type_g,xi_i)
    ! Selection of subtriangles to integrate
    ! If the collocation point is not on an edge, all triangles must be integrated
    if (i_in_edge.eqv.(.false.)) then
      select case (type_g)
        ! Quadrilaterals
        case (fbem_quad4,fbem_quad8,fbem_quad9)
          ! Number of subtriangles
          nsubtri=8
        case (fbem_tri3,fbem_tri6)
          ! Number of subtriangles
          nsubtri=6
      end select
      ! Mark subtriangles
      do ksubtri=1,nsubtri
        subtriangle(ksubtri)=ksubtri
      end do
    ! If the collocation point is on an edge.
    else
      ! Find which triangles must be integrated
      ! Initialization of ntri
      nsubtri=0
      select case (type_g)
        ! Quadrilaterals
        case (fbem_quad4,fbem_quad8,fbem_quad9)
          ! If collocation point is on node 1
          if ((xi_i(1).le.-1.0d0+check_xi_tolerance).and.(xi_i(2).le.-1.0d0+check_xi_tolerance)) then
            nsubtri=2
            subtriangle(1)=4
            subtriangle(2)=5
          end if
          ! If collocation point is on node 2
          if ((xi_i(1).ge. 1.0d0-check_xi_tolerance).and.(xi_i(2).le.-1.0d0+check_xi_tolerance)) then
            nsubtri=2
            subtriangle(1)=6
            subtriangle(2)=7
          end if
          ! If collocation point is on node 3
          if ((xi_i(1).ge. 1.0d0-check_xi_tolerance).and.(xi_i(2).ge. 1.0d0-check_xi_tolerance)) then
            nsubtri=2
            subtriangle(1)=8
            subtriangle(2)=1
          end if
          ! If collocation point is on node 4
          if ((xi_i(1).le.-1.0d0+check_xi_tolerance).and.(xi_i(2).ge. 1.0d0-check_xi_tolerance)) then
            nsubtri=2
            subtriangle(1)=2
            subtriangle(2)=3
          end if
          ! If it is not on a vertex node
          if (nsubtri.eq.0) then
            ! If collocation point is on edge 1
            if (xi_i(2).le.-1.0d0+check_xi_tolerance) then
              nsubtri=4
              subtriangle(1)=4
              subtriangle(2)=5
              subtriangle(3)=6
              subtriangle(4)=7
            end if
            ! If collocation point is on edge 2
            if (xi_i(1).ge. 1.0d0-check_xi_tolerance) then
              nsubtri=4
              subtriangle(1)=6
              subtriangle(2)=7
              subtriangle(3)=8
              subtriangle(4)=1
            end if
            ! If collocation point is on edge 3
            if (xi_i(2).ge. 1.0d0-check_xi_tolerance) then
              nsubtri=4
              subtriangle(1)=8
              subtriangle(2)=1
              subtriangle(3)=2
              subtriangle(4)=3
            end if
            ! If collocation point is on edge 4
            if (xi_i(1).le.-1.0d0+check_xi_tolerance) then
              nsubtri=4
              subtriangle(1)=2
              subtriangle(2)=3
              subtriangle(3)=4
              subtriangle(4)=5
            end if
          end if
        ! Triangulars
        case (fbem_tri3,fbem_tri6)
          ! If collocation point is on node 1
          if (xi_i(1).ge.1.0d0-check_xi_tolerance) then
            nsubtri=1
            subtriangle(1)=3
          end if
          ! If collocation point is on node 2
          if (xi_i(2).ge.1.0d0-check_xi_tolerance) then
            nsubtri=1
            subtriangle(1)=6
          end if
          ! If collocation point is on node 3
          if ((xi_i(1).le.check_xi_tolerance).and.(xi_i(2).le.check_xi_tolerance)) then
            nsubtri=2
            subtriangle(1)=1
            subtriangle(2)=2
          end if
          ! If it is not on a vertex node
          if (nsubtri.eq.0) then
            ! If collocation point is on edge 1
            if ((xi_i(1)+xi_i(2)).ge.1.0d0-check_xi_tolerance) then
              nsubtri=4
              subtriangle(1)=3
              subtriangle(2)=4
              subtriangle(3)=5
              subtriangle(4)=6
            end if
            ! If collocation point is on edge 2
            if (xi_i(1).le.check_xi_tolerance) then
              nsubtri=3
              subtriangle(1)=6
              subtriangle(2)=1
              subtriangle(3)=2
            end if
            ! If collocation point is on edge 3
            if (xi_i(2).le.check_xi_tolerance) then
              nsubtri=3
              subtriangle(1)=1
              subtriangle(2)=2
              subtriangle(3)=3
            end if
          end if
      end select
    end if
    !
    ! Initial (thetai) and final (thetaf) angles for each subtriangle to be integrated, and its tranformation to thetap space.
    !
    do ksubtri=1,nsubtri
      select case (type_g)
        ! Quadrilaterals
        case (fbem_quad4,fbem_quad8,fbem_quad9)
          select case (subtriangle(ksubtri))
            case (1)
              theta1=c_pi -asin((-1.0d0-xi_i(2))/sqrt((-1.0d0-xi_i(1))**2+(-1.0d0-xi_i(2))**2))
              theta_subtri(1,ksubtri) =theta1
              theta_subtri(2,ksubtri) =1.5d0*c_pi
              thetap_subtri(1,ksubtri)=( 1.0d0+xi_i(2))*log(tan(0.5d0*(theta_subtri(1,ksubtri)-c_pi)))
              thetap_subtri(2,ksubtri)=( 1.0d0+xi_i(2))*log(tan(0.5d0*(theta_subtri(2,ksubtri)-c_pi)))
            case (2)
              theta2=c_2pi+asin((-1.0d0-xi_i(2))/sqrt(( 1.0d0-xi_i(1))**2+(-1.0d0-xi_i(2))**2))
              theta_subtri(1,ksubtri) =1.5d0*c_pi
              theta_subtri(2,ksubtri) =theta2
              thetap_subtri(1,ksubtri)=( 1.0d0+xi_i(2))*log(tan(0.5d0*(theta_subtri(1,ksubtri)-c_pi)))
              thetap_subtri(2,ksubtri)=( 1.0d0+xi_i(2))*log(tan(0.5d0*(theta_subtri(2,ksubtri)-c_pi)))
            case (3)
              theta2=c_2pi+asin((-1.0d0-xi_i(2))/sqrt(( 1.0d0-xi_i(1))**2+(-1.0d0-xi_i(2))**2))
              theta_subtri(1,ksubtri) =theta2
              theta_subtri(2,ksubtri) =c_2pi
              thetap_subtri(1,ksubtri)=( 1.0d0-xi_i(1))*log(tan(0.5d0*(theta_subtri(1,ksubtri)+c_pi_2)))
              thetap_subtri(2,ksubtri)=( 1.0d0-xi_i(1))*log(tan(0.5d0*(theta_subtri(2,ksubtri)+c_pi_2)))
            case (4)
              theta3=      asin(( 1.0d0-xi_i(2))/sqrt(( 1.0d0-xi_i(1))**2+( 1.0d0-xi_i(2))**2))
              theta_subtri(1,ksubtri) =0.0d0
              theta_subtri(2,ksubtri) =theta3
              thetap_subtri(1,ksubtri)=( 1.0d0-xi_i(1))*log(tan(0.5d0*(theta_subtri(1,ksubtri)+c_pi_2)))
              thetap_subtri(2,ksubtri)=( 1.0d0-xi_i(1))*log(tan(0.5d0*(theta_subtri(2,ksubtri)+c_pi_2)))
            case (5)
              theta3=      asin(( 1.0d0-xi_i(2))/sqrt(( 1.0d0-xi_i(1))**2+( 1.0d0-xi_i(2))**2))
              theta_subtri(1,ksubtri) =theta3
              theta_subtri(2,ksubtri) =c_pi_2
              thetap_subtri(1,ksubtri)=( 1.0d0-xi_i(2))*log(tan(0.5d0*theta_subtri(1,ksubtri)))
              thetap_subtri(2,ksubtri)=( 1.0d0-xi_i(2))*log(tan(0.5d0*theta_subtri(2,ksubtri)))
            case (6)
              theta4=c_pi -asin(( 1.0d0-xi_i(2))/sqrt((-1.0d0-xi_i(1))**2+( 1.0d0-xi_i(2))**2))
              theta_subtri(1,ksubtri) =c_pi_2
              theta_subtri(2,ksubtri) =theta4
              thetap_subtri(1,ksubtri)=( 1.0d0-xi_i(2))*log(tan(0.5d0*theta_subtri(1,ksubtri)))
              thetap_subtri(2,ksubtri)=( 1.0d0-xi_i(2))*log(tan(0.5d0*theta_subtri(2,ksubtri)))
            case (7)
              theta4=c_pi -asin(( 1.0d0-xi_i(2))/sqrt((-1.0d0-xi_i(1))**2+( 1.0d0-xi_i(2))**2))
              theta_subtri(1,ksubtri) =theta4
              theta_subtri(2,ksubtri) =c_pi
              thetap_subtri(1,ksubtri)=( 1.0d0+xi_i(1))*log(tan(0.5d0*(theta_subtri(1,ksubtri)-c_pi_2)))
              thetap_subtri(2,ksubtri)=( 1.0d0+xi_i(1))*log(tan(0.5d0*(theta_subtri(2,ksubtri)-c_pi_2)))
            case (8)
              theta1=c_pi -asin((-1.0d0-xi_i(2))/sqrt((-1.0d0-xi_i(1))**2+(-1.0d0-xi_i(2))**2))
              theta_subtri(1,ksubtri) =c_pi
              theta_subtri(2,ksubtri) =theta1
              thetap_subtri(1,ksubtri)=( 1.0d0+xi_i(1))*log(tan(0.5d0*(theta_subtri(1,ksubtri)-c_pi_2)))
              thetap_subtri(2,ksubtri)=( 1.0d0+xi_i(1))*log(tan(0.5d0*(theta_subtri(2,ksubtri)-c_pi_2)))
          end select
        ! Triangulars
        case (fbem_tri3,fbem_tri6)
          select case (subtriangle(ksubtri))
            case (1)
              theta1=c_2pi-asin(       xi_i(2) /sqrt((1.0d0-xi_i(1))**2+       xi_i(2)**2 ))
              theta_subtri(1,ksubtri) =theta1
              theta_subtri(2,ksubtri) =c_pi_4+c_2pi
              thetap_subtri(1,ksubtri)=(1.0d0-xi_i(1)-xi_i(2))/c_sqrt2*log(tan(0.5d0*(theta_subtri(1,ksubtri)+c_pi_4)))
              thetap_subtri(2,ksubtri)=(1.0d0-xi_i(1)-xi_i(2))/c_sqrt2*log(tan(0.5d0*(theta_subtri(2,ksubtri)+c_pi_4)))
            case (2)
              theta2=c_pi -asin((1.0d0-xi_i(2))/sqrt(       xi_i(1)**2 +(1.0d0-xi_i(2))**2))
              theta_subtri(1,ksubtri) =c_pi_4
              theta_subtri(2,ksubtri) =theta2
              thetap_subtri(1,ksubtri)=(1.0d0-xi_i(1)-xi_i(2))/c_sqrt2*log(tan(0.5d0*(theta_subtri(1,ksubtri)+c_pi_4)))
              thetap_subtri(2,ksubtri)=(1.0d0-xi_i(1)-xi_i(2))/c_sqrt2*log(tan(0.5d0*(theta_subtri(2,ksubtri)+c_pi_4)))
            case (3)
              theta2=c_pi -asin((1.0d0-xi_i(2))/sqrt(       xi_i(1)**2 +(1.0d0-xi_i(2))**2))
              theta_subtri(1,ksubtri) =theta2
              theta_subtri(2,ksubtri) =c_pi
              thetap_subtri(1,ksubtri)=(      xi_i(1)                )*log(tan(0.5d0*(theta_subtri(1,ksubtri)-c_pi_2)))
              thetap_subtri(2,ksubtri)=(      xi_i(1)                )*log(tan(0.5d0*(theta_subtri(2,ksubtri)-c_pi_2)))
            case (4)
              theta3=c_pi +asin(       xi_i(2) /sqrt(       xi_i(1)**2 +       xi_i(2)**2 ))
              theta_subtri(1,ksubtri) =c_pi
              theta_subtri(2,ksubtri) =theta3
              thetap_subtri(1,ksubtri)=(      xi_i(1)                )*log(tan(0.5d0*(theta_subtri(1,ksubtri)-c_pi_2)))
              thetap_subtri(2,ksubtri)=(      xi_i(1)                )*log(tan(0.5d0*(theta_subtri(2,ksubtri)-c_pi_2)))
            case (5)
              theta3=c_pi +asin(       xi_i(2) /sqrt(       xi_i(1)**2 +       xi_i(2)**2 ))
              theta_subtri(1,ksubtri) =theta3
              theta_subtri(2,ksubtri) =1.5d0*c_pi
              thetap_subtri(1,ksubtri)=(              xi_i(2)        )*log(tan(0.5d0*(theta_subtri(1,ksubtri)-c_pi)))
              thetap_subtri(2,ksubtri)=(              xi_i(2)        )*log(tan(0.5d0*(theta_subtri(2,ksubtri)-c_pi)))
            case (6)
              theta1=c_2pi-asin(       xi_i(2) /sqrt((1.0d0-xi_i(1))**2+       xi_i(2)**2 ))
              theta_subtri(1,ksubtri) =1.5d0*c_pi
              theta_subtri(2,ksubtri) =theta1
              thetap_subtri(1,ksubtri)=(              xi_i(2)        )*log(tan(0.5d0*(theta_subtri(1,ksubtri)-c_pi)))
              thetap_subtri(2,ksubtri)=(              xi_i(2)        )*log(tan(0.5d0*(theta_subtri(2,ksubtri)-c_pi)))
          end select
      end select
    end do
  end subroutine

  !! This subroutine performs the angular transformation
  subroutine fbem_polar_transformation_angular(type_g,xi_i,subtriangle,thetap,theta,rhoij)
    implicit none
    ! I/O
    integer           :: type_g                          !! Geometrial interpolation
    real(kind=real64) :: xi_i(2)                         !! Reference coordinates of the singular point.
    integer           :: subtriangle                     !! Selected subtriangle
    real(kind=real64) :: thetap                          !! Angle coordinate thetap
    real(kind=real64) :: theta                           !! Angle coordinate theta
    real(kind=real64) :: rhoij                           !! Maximum rho (radial)
    !
    ! Transformation thetap->theta, and calculation of rhoij (edge maximum rho)
    !
    select case (type_g)
      ! Quadrilaterals
      case (fbem_quad4,fbem_quad8,fbem_quad9)
        select case (subtriangle)
          case (1,2)
            theta=2.0d0*atan(exp(thetap/( 1.0d0+xi_i(2))))+c_pi
            rhoij=(-1.0d0-xi_i(2))/sin(theta)
          case (3,4)
            theta=2.0d0*atan(exp(thetap/( 1.0d0-xi_i(1))))-c_pi_2
            rhoij=( 1.0d0-xi_i(1))/cos(theta)
          case (5,6)
            theta=2.0d0*atan(exp(thetap/( 1.0d0-xi_i(2))))
            rhoij=( 1.0d0-xi_i(2))/sin(theta)
          case (7,8)
            theta=2.0d0*atan(exp(thetap/( 1.0d0+xi_i(1))))+c_pi_2
            rhoij=(-1.0d0-xi_i(1))/cos(theta)
        end select
      ! Triangulars
      case (fbem_tri3,fbem_tri6)
        select case (subtriangle)
          case (1,2)
            theta=2.0d0*atan(exp(c_sqrt2*thetap/(1.0d0-xi_i(1)-xi_i(2))))-c_pi_4
            rhoij=(1.0d0-xi_i(1)-xi_i(2))/(cos(theta)+sin(theta))
          case (3,4)
            theta=2.0d0*atan(exp(thetap/xi_i(1)))+c_pi_2
            rhoij=-xi_i(1)/cos(theta)
          case (5,6)
            theta=2.0d0*atan(exp(thetap/xi_i(2)))+c_pi
            rhoij=-xi_i(2)/sin(theta)
        end select
    end select
  end subroutine

end module fbem_polar_transformation
