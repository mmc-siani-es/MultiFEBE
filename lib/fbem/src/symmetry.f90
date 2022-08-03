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
!! <b> This module implements symmetry related routines. </b>
!!
!! TODO: make full use of the module
module fbem_symmetry

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! needed fbem module
  use fbem_numerical
  use fbem_string_handling

  ! No implicit variables are allowed
  implicit none

  ! By default all are private
  private

  ! Public routines
  public :: fbem_symmetry_multipliers
  public :: fbem_reflect_position_2d
  public :: fbem_reflect_position_3d
  public :: fbem_reflect_vector_2d
  public :: fbem_reflect_vector_3d
  public :: fbem_rotate_vector_2d
  public :: fbem_rotate_vector_3d
  public :: fbem_check_compatible_splanes_configuration_2d
  public :: fbem_check_compatible_splanes_configuration_3d
  public :: fbem_set_of_points_reflection_2d
  public :: fbem_set_of_points_reflection_3d

contains

  !! This subroutine calculates the multipliers for each symmetrical element for a configuration of symmetry planes.
  subroutine fbem_symmetry_multipliers(step,rn,n_symplanes,symplane_m,symplane_s,symplane_t,symplane_r,&
                                                           symconf_m,symconf_s,symconf_t,symconf_r,reversed)
    implicit none
    ! I/O variables
    integer           :: step                       !! Step of the process, each symmetrical element
    integer           :: rn                         !! Dimensional space
    integer           :: n_symplanes                !! Number of symmetry planes
    real(kind=real64) :: symplane_m(3,3)            !! Geometric multipliers for each symmetry plane
    real(kind=real64) :: symplane_s(3)              !! Scalar variable multiplier for each symmetry plane
    real(kind=real64) :: symplane_t(3,3)            !! Translation variable multiplier for each symmetry plane
    real(kind=real64) :: symplane_r(3,3)            !! Rotation variable multiplier for each symmetry plane
    real(kind=real64) :: symconf_m(rn)              !! Geometric multipliers for the selected step
    real(kind=real64) :: symconf_s                  !! Scalar variable multipliers for the selected step
    real(kind=real64) :: symconf_t(rn)              !! Translation variable multipliers for the selected step
    real(kind=real64) :: symconf_r(rn)              !! Rotation variable multipliers for the selected step
    logical           :: reversed                   !! True if the element is reversed
    ! Local variables
    integer           :: kc
    !
    ! Multipliers for each step
    !
    select case (step)
      !
      ! Root element
      !
      case (1)
        do kc=1,rn
          symconf_m(kc)=1.0d0
          symconf_s    =1.0d0
          symconf_t(kc)=1.0d0
          symconf_r(kc)=1.0d0
        end do
        reversed=.false.
      !
      ! Element reflected with respect to the SP1
      !
      case (2)
        do kc=1,rn
          symconf_m(kc)=symplane_m(kc,1)
          symconf_s    =symplane_s(1)
          symconf_t(kc)=symplane_t(kc,1)
          symconf_r(kc)=symplane_r(kc,1)
        end do
        reversed=.true.
      !
      ! Element reflected with respect to the SP1 and SP2
      !
      case (3)
        do kc=1,rn
          symconf_m(kc)=symplane_m(kc,1)*symplane_m(kc,2)
          symconf_s    =symplane_s(1)*symplane_s(2)
          symconf_t(kc)=symplane_t(kc,1)*symplane_t(kc,2)
          symconf_r(kc)=symplane_r(kc,1)*symplane_r(kc,2)
        end do
        reversed=.false.
      !
      ! Element reflected with respect to the SP2
      !
      case (4)
        do kc=1,rn
          symconf_m(kc)=symplane_m(kc,2)
          symconf_s    =symplane_s(2)
          symconf_t(kc)=symplane_t(kc,2)
          symconf_r(kc)=symplane_r(kc,2)
        end do
        reversed=.true.
      !
      ! Element reflected with respect to the SP3
      !
      case (5)
        do kc=1,rn
          symconf_m(kc)=symplane_m(kc,3)
          symconf_s=symplane_s(3)
          symconf_t(kc)=symplane_t(kc,3)
          symconf_r(kc)=symplane_r(kc,3)
        end do
        reversed=.true.
      !
      ! Element reflected with respect to the SP1 and SP3
      !
      case (6)
        do kc=1,rn
          symconf_m(kc)=symplane_m(kc,1)*symplane_m(kc,3)
          symconf_s=symplane_s(1)*symplane_s(3)
          symconf_t(kc)=symplane_t(kc,1)*symplane_t(kc,3)
          symconf_r(kc)=symplane_r(kc,1)*symplane_r(kc,3)
        end do
        reversed=.false.
      !
      ! Element reflected with respect to the SP1, SP2 and SP3
      !
      case (7)
        do kc=1,rn
          symconf_m(kc)=symplane_m(kc,1)*symplane_m(kc,2)*symplane_m(kc,3)
          symconf_s=symplane_s(1)*symplane_s(2)*symplane_s(3)
          symconf_t(kc)=symplane_t(kc,1)*symplane_t(kc,2)*symplane_t(kc,3)
          symconf_r(kc)=symplane_r(kc,1)*symplane_r(kc,2)*symplane_r(kc,3)
        end do
        reversed=.true.
      !
      ! Element reflected with respect to the SP2 and SP3
      !
      case (8)
        do kc=1,rn
          symconf_m(kc)=symplane_m(kc,2)*symplane_m(kc,3)
          symconf_s=symplane_s(2)*symplane_s(3)
          symconf_t(kc)=symplane_t(kc,2)*symplane_t(kc,3)
          symconf_r(kc)=symplane_r(kc,2)*symplane_r(kc,3)
        end do
        reversed=.false.
    end select
  end subroutine fbem_symmetry_multipliers

  !! Calculation of the reflected position vector with respect to a plane in 2D space
  function fbem_reflect_position_2d(r,x,n)
    implicit none
    !! Reflected position vector
    real(kind=real64) :: fbem_reflect_position_2d(2)
    !! Original position vector
    real(kind=real64) :: r(2)
    !! Position vector of a point of the plane
    real(kind=real64) :: x(2)
    !! Unit normal vector of the plane
    real(kind=real64) :: n(2)
    real(kind=real64) :: aux
    ! Formula for point reflection with respect to a line: rr=r-2·n·((r-x)·n)
    aux=(r(1)-x(1))*n(1)+(r(2)-x(2))*n(2)
    fbem_reflect_position_2d(1)=r(1)-2.0d0*n(1)*aux
    fbem_reflect_position_2d(2)=r(2)-2.0d0*n(2)*aux
  end function fbem_reflect_position_2d

  !! Calculation of the reflected position vector with respect to a plane in 3D space
  function fbem_reflect_position_3d(r,x,n)
    implicit none
    !! Reflected position vector
    real(kind=real64) :: fbem_reflect_position_3d(3)
    !! Original position vector
    real(kind=real64) :: r(3)
    !! Position vector of a point of the plane
    real(kind=real64) :: x(3)
    !! Unit normal vector of the plane
    real(kind=real64) :: n(3)
    real(kind=real64) :: aux
    ! Formula for point reflection with respect to a plane: rr=r-2·n·((r-x)·n)
    aux=(r(1)-x(1))*n(1)+(r(2)-x(2))*n(2)+(r(3)-x(3))*n(3)
    fbem_reflect_position_3d(1)=r(1)-2.0d0*n(1)*aux
    fbem_reflect_position_3d(2)=r(2)-2.0d0*n(2)*aux
    fbem_reflect_position_3d(3)=r(3)-2.0d0*n(3)*aux
  end function fbem_reflect_position_3d

  !! Calculation of the reflected vector with respect to a plane in 2D space
  function fbem_reflect_vector_2d(v,n)
    implicit none
    !! Reflected vector
    real(kind=real64) :: fbem_reflect_vector_2d(2)
    !! Original vector
    real(kind=real64) :: v(2)
    !! Unit normal vector of the plane
    real(kind=real64) :: n(2)
    real(kind=real64) :: aux
    ! Formula for vector reflection with respect to a plane: rr=v-2·n·(v·n)
    aux=v(1)*n(1)+v(2)*n(2)
    fbem_reflect_vector_2d(1)=v(1)-2.0d0*n(1)*aux
    fbem_reflect_vector_2d(2)=v(2)-2.0d0*n(2)*aux
  end function fbem_reflect_vector_2d

  !! Calculation of the reflected vector with respect to a plane in 3D space
  function fbem_reflect_vector_3d(v,n)
    implicit none
    !! Reflected vector
    real(kind=real64) :: fbem_reflect_vector_3d(3)
    !! Original vector
    real(kind=real64) :: v(3)
    !! Unit normal vector of the plane
    real(kind=real64) :: n(3)
    real(kind=real64) :: aux
    ! Formula for vector reflection with respect to a plane: rr=v-2·n·(v·n)
    aux=v(1)*n(1)+v(2)*n(2)+v(3)*n(3)
    fbem_reflect_vector_3d(1)=v(1)-2.0d0*n(1)*aux
    fbem_reflect_vector_3d(2)=v(2)-2.0d0*n(2)*aux
    fbem_reflect_vector_3d(3)=v(3)-2.0d0*n(3)*aux
  end function fbem_reflect_vector_3d

  !! Calculation of the rotated vector in 2D space
  function fbem_rotate_vector_2d(v,angle)
    implicit none
    !! Rotated vector
    real(kind=real64) :: fbem_rotate_vector_2d(2)
    !! Original vector
    real(kind=real64) :: v(2)
    !! Rotation angle
    real(kind=real64) :: angle
    real(kind=real64) :: rot(2,2), cos_angle, sin_angle
    integer :: i, j
    ! Save this calculation
    cos_angle=cos(angle)
    sin_angle=sin(angle)
    ! Rotation matrix
    rot(1,1)=cos_angle
    rot(1,2)=-sin_angle
    rot(2,1)=sin_angle
    rot(2,2)=cos_angle
    ! Rotation matrix multiplied by original vector
    do i=1,2
      fbem_rotate_vector_2d(i)=0.0d0
      do j=1,2
        fbem_rotate_vector_2d(i)=fbem_rotate_vector_2d(i)+rot(i,j)*v(j)
      end do
    end do
  end function fbem_rotate_vector_2d

  !! Calculation of the rotated vector with respect to an axis in 3D space
  function fbem_rotate_vector_3d(v,axis,angle)
    implicit none
    !! Rotated vector
    real(kind=real64) :: fbem_rotate_vector_3d(3)
    !! Original vector
    real(kind=real64) :: v(3)
    !! Rotation axis as a unit vector
    real(kind=real64) :: axis(3)
    !! Rotation angle
    real(kind=real64) :: angle
    real(kind=real64) :: rot(3,3), cos_angle, sin_angle, auxcos
    integer :: i, j
    ! Save this calculation
    cos_angle=cos(angle)
    sin_angle=sin(angle)
    auxcos=1.0d0-cos_angle
    ! Rotation matrix
    rot(1,1)=cos_angle+(axis(1)**2)*auxcos
    rot(1,2)=axis(1)*axis(2)*auxcos-axis(3)*sin_angle
    rot(1,3)=axis(1)*axis(3)*auxcos+axis(2)*sin_angle
    rot(2,1)=axis(2)*axis(1)*auxcos+axis(3)*sin_angle
    rot(2,2)=cos_angle+(axis(2)**2)*auxcos
    rot(2,3)=axis(2)*axis(3)*auxcos-axis(1)*sin_angle
    rot(3,1)=axis(3)*axis(1)*auxcos-axis(2)*sin_angle
    rot(3,2)=axis(3)*axis(2)*auxcos+axis(1)*sin_angle
    rot(3,3)=cos_angle+(axis(3)**2)*auxcos
    ! Rotation matrix multiplied by original vector
    do i=1,3
      fbem_rotate_vector_3d(i)=0.0d0
      do j=1,3
        fbem_rotate_vector_3d(i)=fbem_rotate_vector_3d(i)+rot(i,j)*v(j)
      end do
    end do
  end function fbem_rotate_vector_3d

  !! Determine if the symmetry planes configuration in 2D is valid. The valid cases are:
  !!   -n_splanes=0,1,2
  !!   -With n_splanes=2:
  !!     -If both are of the same type (q=0 or u=0 in the symmetry planes), angle alpha between them must be alpha=pi/n, n=2,3,4,...
  !!     -If both are of different type, angle alpha between them must be alpha=pi/n, n=2,4,6,....
  subroutine fbem_check_compatible_splanes_configuration_2d(n_splanes,splane_n,splane_s)
    implicit none
    !! Number of planes
    integer :: n_splanes
    !! Unit normal of each plane
    real(kind=real64) :: splane_n(2,n_splanes)
    !! Plane nature: 1 (symmetric), -1 (antisymmetric)
    integer :: splane_s(n_splanes)
    ! Counters
    integer :: i, j
    real(kind=real64) :: alpha, cosalpha
    ! If n_splanes>1, compatibility must be checked
    !
    if (n_splanes.gt.1) then
      ! Check that n_splanes<=3
      if (n_splanes.gt.2) then
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                'the number of symmetry planes must be <= 2')
      end if
      ! If alpha is the angle between planes: n=pi/alpha, n=2,3,4,...;
      ! If the pair of planes has different nature, n=2,4,6....
      do i=1,n_splanes
        do j=i+1,n_splanes
          ! Because vectors are normalized, scalar product of both normals are the cosine of the angle
          cosalpha=splane_n(1,i)*splane_n(1,j)+splane_n(2,i)*splane_n(2,j)
          ! Angle is degrees (the result will be between 0º and 180º)
          alpha=dacos(cosalpha)*180.0d0/c_pi
          ! Because it is needed the small angle (between 0º and 90º)
          if (alpha.ge.90.0d0) alpha=180.0d0-alpha
          ! If angle is close to 0º (it is used alpha<1º), both planes are nearly parallel
          if (alpha.lt.1.0d0) then
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                    'two symmetry planes are nearly parallel')
          end if
          ! Check the condition of two planes to generate a finite number of source points in fundamental solution:
          ! if alpha is the angle between two planes, n=180/alpha must be integer and >2
          if (((180.0d0-dble(nint(180.0d0/alpha)*alpha)).lt.1.0d-6).and.((180.0d0/alpha).gt.(2.0d0-(1.0d-6)))) then
            ! Both planes are valid
          else
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                    'the angle between two symmetry planes is invalid')
          end if
          ! If both planes has the different nature (one u=0 and the other q=0), besides, n must be even.
          if (splane_s(i).ne.splane_s(j)) then
            if ((180.0d0/alpha-nint(180.0d0/alpha/2.0d0)*2.0d0).lt.1.0d-6) then
              ! Both planes are compatible
            else
              call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                      'the nature of the two symmetry planes are incompatible')
            end if
          end if
        end do
      end do
    end if ! If n_splanes>1
  end subroutine fbem_check_compatible_splanes_configuration_2d

  !! Determine if the symmetry planes configuration in 3D is valid. The valid cases are:
  !!   -n_splanes=0,1,2,3
  !!   -With n_splanes=2:
  !!     -If both are of the same type (q=0 or u=0 in the symmetry planes), angle alpha between them must be alpha=pi/n, n=2,3,4,...
  !!     -If both are of different type, angle alpha between them must be alpha=pi/n, n=2,4,6,....
  !!   -With n_splanes=3:
  !!     -Only combinations of planes verify conditions of n_splanes=2, and one plane is ortoghonal to the other two.
  subroutine fbem_check_compatible_splanes_configuration_3d(n_splanes,splane_n,splane_s)
    implicit none
    !! Number of planes
    integer :: n_splanes
    !! Unit normal of each plane
    real(kind=real64) :: splane_n(3,n_splanes)
    !! Plane nature: 1 (symmetric), -1 (antisymmetric)
    integer :: splane_s(n_splanes)
    ! Counters
    integer :: i, j, k
    real(kind=real64) :: alpha, cosalpha
    ! If n_splanes>1, compatibility must be checked
    !
    if (n_splanes.gt.1) then
      ! Check that n_splanes<=3
      if (n_splanes.gt.3) then
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                'the number of symmetry planes must be <= 3')
      end if
      ! All plane pairs must verify that if alpha is the angle between planes: n=pi/alpha, n=2,3,4,...; if the pair of planes
      ! has different nature, n=2,4,6....
      ! Counter of orthogonal angles between planes
      k=0
      do i=1,n_splanes
        do j=i+1,n_splanes
          ! Because vectors are normalized, scalar product of both normals are the cosine of the angle
          cosalpha=splane_n(1,i)*splane_n(1,j)+splane_n(2,i)*splane_n(2,j)+splane_n(3,i)*splane_n(3,j)
          ! Angle is degrees (the result will be between 0º and 180º)
          alpha=dacos(cosalpha)*180.0d0/c_pi
          ! Because it is needed the small angle (between 0º and 90º)
          if (alpha.ge.90.0d0) alpha=180.0d0-alpha
          ! Increment the counter of orthogonal angles if alpha=90
          if ((alpha.gt.(90.0d0-(1.0d-6))).and.(alpha.lt.(90.0d0+(1.0d-6)))) k=k+1
          ! If angle is close to 0º (it is used alpha<1º), both planes are nearly parallel
          if (alpha.lt.1.0d0) then
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                    'two symmetry planes are nearly parallel')
          end if
          ! Check the condition of two planes to generate a finite number of source points in fundamental solution:
          ! if alpha is the angle between two planes, n=180/alpha must be integer and >2
          if (((180.0d0-dble(nint(180.0d0/alpha)*alpha)).lt.1.0d-6).and.((180.0d0/alpha).gt.(2.0d0-(1.0d-6)))) then
            ! Both planes are valid
          else
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                    'the angle between two symmetry planes is invalid')
          end if
          ! If both planes has the different nature (one u=0 and the other q=0), besides, n must be even.
          if (splane_s(i).ne.splane_s(j)) then
            if ((180.0d0/alpha-nint(180.0d0/alpha/2.0d0)*2.0d0).lt.1.0d-6) then
              ! Both planes are compatible
            else
              call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                      'the nature of the two symmetry planes are incompatible')
            end if
          end if
        end do
      end do
      ! If there are 3 planes, one of them must be orthogonal to the other two.
      if (n_splanes.eq.3) then
        if (k.lt.2) then
          call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                  'one plane must be orthogonal to the other two')
        end if
      end if
    end if ! If n_splanes>1
  end subroutine fbem_check_compatible_splanes_configuration_3d

  !! Calculate the finite set of points reflected on a set of planes in 2D space
  subroutine fbem_set_of_points_reflection_2d(x_i,n_i,n_splanes,splane_x,splane_n,splane_s,tolerance_splane,&
                                              max_n_points,x_list,n_list,s_list,n_points)
    implicit none
    !! Position vector of the root point
    real(kind=real64) :: x_i(2)
    !! Unit normal vector of the root point
    real(kind=real64) :: n_i(2)
    !! Number of planes
    integer :: n_splanes
    !! Position vectors of a point of each plane
    real(kind=real64) :: splane_x(2,n_splanes)
    !! Unit normal vectors of each plane
    real(kind=real64) :: splane_n(2,n_splanes)
    !! Vector that contains (anti)symmetry plane behaviour:
    !! -( 1): same sign
    !! -(-1): opposite sign
    integer :: splane_s(n_splanes)
    !! Tolerance used to determine if a point is similar to other
    real(kind=real64) :: tolerance_splane
    !! Maximum number of points
    integer :: max_n_points
    !! List of position vectors of points
    real(kind=real64) :: x_list(2,max_n_points)
    !! List of unit normal vectors of points
    real(kind=real64) :: n_list(2,max_n_points)
    !! List of point signs
    integer :: s_list(max_n_points)
    !! Number of points generated
    integer :: n_points
    ! Loop counter variables
    integer :: i, j
    ! Matrix to check if a point has been reflected
    logical :: check_x(max_n_points,n_splanes)
    ! Counter of the analysis point
    integer :: actual_point
    ! Counter of the point of maximum index
    integer :: max_point
    ! Control variables
    logical :: new_reflection, in_list
    ! Distance and temporal point and normal
    real(kind=real64) :: r, x_tmp(2), n_tmp(2)
    ! Initialize reflection control variable
    do i=1,max_n_points
      do j=1,n_splanes
        check_x(i,j)=.false.
      end do
    end do
    ! Initialize counters
    actual_point=1
    max_point=1
    ! Save root point
    ! Position vector
    x_list(1,actual_point)=x_i(1)
    x_list(2,actual_point)=x_i(2)
    ! Unit normal vector
    n_list(1,actual_point)=n_i(1)
    n_list(2,actual_point)=n_i(2)
    ! Save root point sign
    s_list(actual_point)=1
    ! Initialize loop control variable
    new_reflection=.true.
    ! Do while
    do while (new_reflection.eqv.(.true.))
      ! Loop through symmetry planes for the actual point
      do i=1,n_splanes
        ! If the actual point has not been reflected with respect to plane i
        if (check_x(actual_point,i).eqv.(.false.)) then
          ! Reflect it
          x_tmp=fbem_reflect_position_2d(x_list(:,actual_point),splane_x(:,i),splane_n(:,i))
          n_tmp=fbem_reflect_vector_2d(n_list(:,actual_point),splane_n(:,i))
          ! Check if the reflected point is in the list
          in_list=.false.
          do j=1,max_point
            ! Calculate distance between the reflect point and the point j in the list
            r=dsqrt((x_tmp(1)-x_list(1,j))**2+(x_tmp(2)-x_list(2,j))**2)
            ! If the distance is less or equal than the tolerance
            if (r.le.tolerance_splane) then
              ! Control variable that indicate if x_tmp is already in the list
              in_list=.true.
              ! This point is yet a reflection
              check_x(actual_point,i)=.true.
              exit
            end if
          end do
          ! If the reflected point is not in the list
          if (in_list.eqv.(.false.)) then
            ! Increment the index of the number of points
            max_point=max_point+1
            ! Register that the point has been reflected
            check_x(actual_point,i)=.true.
            check_x(max_point,i)=.true.
            ! Save position vector in the list
            x_list(1,max_point)=x_tmp(1)
            x_list(2,max_point)=x_tmp(2)
            ! Save unit normal vector in the list
            n_list(1,max_point)=n_tmp(1)
            n_list(2,max_point)=n_tmp(2)
            ! Save point sign
            s_list(max_point)=splane_s(i)*s_list(actual_point)
          end if
        end if
      end do
      ! For the root point
      if (actual_point.eq.1) then
        ! Increment actual point counter if new points has been generated
        if (max_point.gt.actual_point) then
          actual_point=actual_point+1
        ! If new points hasn't been generated, terminate, because it is on the symmetry planes intersection
        else
          new_reflection=.false.
        end if
      else
        ! Increment actual point counter if new points has been generated
        if (max_point.gt.actual_point) actual_point=actual_point+1
        ! Check if new actual point has been reflected to all planes, if so, terminate.
        new_reflection=.false.
        do i=1,n_splanes
          if (check_x(actual_point,i).eqv.(.false.)) then
            new_reflection=.true.
            exit
          end if
        end do
      end if
    end do
    ! The number of generated points
    n_points=max_point
  end subroutine fbem_set_of_points_reflection_2d

  !! Calculate the finite set of points reflected on a set of planes in 3D space
  subroutine fbem_set_of_points_reflection_3d(x_i,n_i,n_splanes,splane_x,splane_n,splane_s,tolerance_splane,&
                                              max_n_points,x_list,n_list,s_list,n_points)
    implicit none
    !! Position vector of the root point
    real(kind=real64) :: x_i(3)
    !! Unit normal vector of the root point
    real(kind=real64) :: n_i(3)
    !! Number of planes
    integer :: n_splanes
    !! Position vectors of a point of each plane
    real(kind=real64) :: splane_x(3,n_splanes)
    !! Unit normal vectors of each plane
    real(kind=real64) :: splane_n(3,n_splanes)
    !! Vector that contains (anti)symmetry plane behaviour:
    !! -( 1): same sign
    !! -(-1): opposite sign
    integer :: splane_s(n_splanes)
    !! Tolerance used to determine if a point is similar to other
    real(kind=real64) :: tolerance_splane
    !! Maximum number of points
    integer :: max_n_points
    !! List of position vectors of points
    real(kind=real64) :: x_list(3,max_n_points)
    !! List of unit normal vectors of points
    real(kind=real64) :: n_list(3,max_n_points)
    !! List of point signs
    integer :: s_list(max_n_points)
    !! Number of points generated
    integer :: n_points
    ! Loop counter variables
    integer :: i, j
    ! Matrix to check if a point has been reflected
    logical :: check_x(max_n_points,n_splanes)
    ! Counter of the analysis point
    integer :: actual_point
    ! Counter of the point of maximum index
    integer :: max_point
    ! Control variables
    logical :: new_reflection, in_list
    ! Distance and temporal point and normal
    real(kind=real64) :: r, x_tmp(3), n_tmp(3)
    ! Initialize reflection control variable
    do i=1,max_n_points
      do j=1,n_splanes
        check_x(i,j)=.false.
      end do
    end do
    ! Initialize counters
    actual_point=1
    max_point=1
    ! Save root point
    ! Position vector
    x_list(1,actual_point)=x_i(1)
    x_list(2,actual_point)=x_i(2)
    x_list(3,actual_point)=x_i(3)
    ! Unit normal vector
    n_list(1,actual_point)=n_i(1)
    n_list(2,actual_point)=n_i(2)
    n_list(3,actual_point)=n_i(3)
    ! Save root point sign
    s_list(actual_point)=1
    ! Initialize loop control variable
    new_reflection=.true.
    ! Do while
    do while (new_reflection.eqv.(.true.))
      ! Loop through symmetry planes for the actual point
      do i=1,n_splanes
        ! If the actual point has not been reflected with respect to plane i
        if (check_x(actual_point,i).eqv.(.false.)) then
          ! Reflect it
          x_tmp=fbem_reflect_position_3d(x_list(:,actual_point),splane_x(:,i),splane_n(:,i))
          n_tmp=fbem_reflect_vector_3d(n_list(:,actual_point),splane_n(:,i))
          ! Check if the reflected point is in the list
          in_list=.false.
          do j=1,max_point
            ! Calculate distance between the reflect point and the point j in the list
            r=dsqrt((x_tmp(1)-x_list(1,j))**2+(x_tmp(2)-x_list(2,j))**2+(x_tmp(3)-x_list(3,j))**2)
            ! If the distance is less or equal than the tolerance
            if (r.le.tolerance_splane) then
              ! Control variable that indicate if x_tmp is already in the list
              in_list=.true.
              ! This point is yet a reflection
              check_x(actual_point,i)=.true.
              exit
            end if
          end do
          ! If the reflected point is not in the list
          if (in_list.eqv.(.false.)) then
            ! Increment the index of the number of points
            max_point=max_point+1
            ! Register that the point has been reflected
            check_x(actual_point,i)=.true.
            check_x(max_point,i)=.true.
            ! Save position vector in the list
            x_list(1,max_point)=x_tmp(1)
            x_list(2,max_point)=x_tmp(2)
            x_list(3,max_point)=x_tmp(3)
            ! Save unit normal vector in the list
            n_list(1,max_point)=n_tmp(1)
            n_list(2,max_point)=n_tmp(2)
            n_list(3,max_point)=n_tmp(3)
            ! Save point sign
            s_list(max_point)=splane_s(i)*s_list(actual_point)
          end if
        end if
      end do
      ! For the root point
      if (actual_point.eq.1) then
        ! Increment actual point counter if new points has been generated
        if (max_point.gt.actual_point) then
          actual_point=actual_point+1
        ! If new points hasn't been generated, terminate, because it is on the symmetry planes intersection
        else
          new_reflection=.false.
        end if
      else
        ! Increment actual point counter if new points has been generated
        if (max_point.gt.actual_point) actual_point=actual_point+1
        ! Check if new actual point has been reflected to all planes, if so, terminate.
        new_reflection=.false.
        do i=1,n_splanes
          if (check_x(actual_point,i).eqv.(.false.)) then
            new_reflection=.true.
            exit
          end if
        end do
      end if
    end do
    ! The number of generated points
    n_points=max_point
  end subroutine fbem_set_of_points_reflection_3d

end module fbem_symmetry
