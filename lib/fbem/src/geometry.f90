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

! For safety reasons, is possible to compile the code with a xi or xi1,xi2 domain checking, e.g., in every function call
! the xi or xi1,xi2 are checked to see if they are in the domain. If CHECK_XIS == 1, checking is enabled, and a warning message
! is displayed. If CHECK_XIS == 2, checking is enabled, an error message is displayed and execution stops.Other values of disable
! any checking.
#define CHECK_XIS 0

!! @author Jacob David Rodriguez Bordon (jacobdavid.rodriguezbordon@ulpgc.es)
!!
!! @version 2.0
!!
!! <b> This module implements geometrical functions and subroutines about elements. </b>
module fbem_geometry

  ! Fortran 2003 standard
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_numerical
  use fbem_string_handling
  use fbem_quad_rules
  use fbem_telles_transformation
  use fbem_shape_functions

  ! No implicit variables are allowed
  implicit none

  ! By default all are private
  private

  ! Coordinate transformation
  public :: fbem_coordinate_transformation_L

  ! Element basic vectors
  public :: fbem_position
  public :: fbem_position2d
  public :: fbem_position3d
  public :: fbem_tangent_xi
  public :: fbem_tangent_xi1
  public :: fbem_tangent_xi2
  public :: fbem_utangent_xi
  public :: fbem_utangent_xi1
  public :: fbem_utangent_xi2
  public :: fbem_tangent_derivatives_2d
  public :: fbem_utangents_at_boundary
  public :: fbem_normal2d
  public :: fbem_unormal2d
  public :: fbem_normal3d
  public :: fbem_unormal3d
  public :: fbem_normal2d_and_tangent_to_center   ! to be deprecated
  public :: fbem_normal3d_and_bis_tangent_at_node ! to be deprecated
  public :: fbem_local_tangential_cartesian_axes

  ! Element jacobians
  public :: fbem_jacobian2d
  public :: fbem_jacobian3d

  ! Dataset values at integration points
  public :: fbem_dataset_at_integration_points
  public :: fbem_dataset_at_integration_points_dif

  ! Element subdivision
  public :: fbem_obtain_element_subdivision_coordinates

  ! Element properties
  public :: fbem_length2d
  public :: fbem_length3d
  public :: fbem_length2d_subdivision
  public :: fbem_characteristic_length
  public :: fbem_characteristic_length_subdivision
  public :: fbem_area2d
  public :: fbem_area3d
  public :: fbem_element_size
  public :: fbem_element_size_ngp
  public :: fbem_element_centroid
  public :: fbem_element_centroid_ngp
  public :: fbem_geometry_element_ball

  ! Element-point functions
  public :: fbem_nearest_xi_nodes              ! For 1D and 2D elements (using only their nodes)
  public :: fbem_nearest_xi                    ! For 1D element using a sampling algorithm
  public :: fbem_nearest_xi_subdivision        ! For a part of a 1D element using a sampling algorithm
  public :: fbem_nearest_minimization_1d       ! For 2D element using a minimization algorithm
  public :: fbem_nearest_xi1xi2                ! For 2D element using a sampling algorithm
  public :: fbem_nearest_xi1xi2_subdivision    ! For a part of a 2D element using a sampling algorithm
  public :: fbem_nearest_minimization_2d       ! For 2D element using a minimization algorithm
  public :: fbem_nearest_element_point_bem     ! For 1D and 2D elements (adaptative algorithm using all methods efficiently)
  public :: fbem_nearest_minimization_3d_shell ! For degenerated shell elements
  public :: fbem_local_coordinates

  ! Transformations
  public :: fbem_rigid_solid_transformation_matrix
  public :: fbem_rotation_infinitesimal_displacement
  public :: fbem_transformation_collapse_nodal_positions
  public :: fbem_move_xi_from_vertex
  public :: fbem_move_xi1xi2_from_edge

  !! Common name for <tt>position2d_1d</tt> and <tt>position2d_2d</tt>.
  interface fbem_position2d
    module procedure position2d_1d
    module procedure position2d_2d
  end interface fbem_position2d

  !! Common name for <tt>position3d_1d</tt> and <tt>position3d_nd</tt>.
  interface fbem_position3d
    module procedure position3d_1d
    module procedure position3d_nd
  end interface fbem_position3d

  !! Common name for <tt>fbem_jacobian2d_1d</tt> and <tt>fbem_jacobian2d_2d</tt>.
  interface fbem_jacobian2d
    module procedure fbem_jacobian2d_1d
    module procedure fbem_jacobian2d_2d
  end interface fbem_jacobian2d

  !! Common name for <tt>fbem_jacobian3d_1d</tt> and <tt>fbem_jacobian3d_nd</tt>.
  interface fbem_jacobian3d
    module procedure fbem_jacobian3d_1d
    module procedure fbem_jacobian3d_nd
  end interface fbem_jacobian3d

contains

  !! Transformation matrix for displacements/rotations "a" and forces/moments "p" in a rigid solid:
  !! a_s = T  · a_m
  !! p_m = T' · p_s
  subroutine fbem_rigid_solid_transformation_matrix(rn,x_s,x_m,n_dof_s,T)
    integer           :: rn                  !! Ambient dimension
    real(kind=real64) :: x_s(rn)             !! Slave node position
    real(kind=real64) :: x_m(rn)             !! Master node position
    integer           :: n_dof_s             !! Number of degrees of freedom of the slave node
    real(kind=real64) :: T(n_dof_s,3*(rn-1)) !! Transformation matrix
    real(kind=real64) :: r(rn)               !! Distance vector
    T=0.d0
    r=x_s-x_m
    select case (rn)
      case (2)
        T(1,1)=1.d0
        T(2,2)=1.d0
        T(1,3)=-r(2)
        T(2,3)= r(1)
        if (n_dof_s.eq.3) then
          T(3,3)=1.d0
        end if
      case (3)
        T(1,1)=1.d0
        T(2,2)=1.d0
        T(3,3)=1.d0
        T(2,4)=-r(3)
        T(3,4)= r(2)
        T(1,5)= r(3)
        T(3,5)=-r(1)
        T(1,6)=-r(2)
        T(2,6)= r(1)
        if (n_dof_s.eq.6) then
          T(4,4)=1.d0
          T(5,5)=1.d0
          T(6,6)=1.d0
        end if
    end select
  end subroutine fbem_rigid_solid_transformation_matrix

  !! Infinitesimal displacement due to a rotation
  function fbem_rotation_infinitesimal_displacement(axis_point,axis_vector,angle,point)
    implicit none
    ! I/O
    real(kind=real64) :: fbem_rotation_infinitesimal_displacement(3) !! Displacement vector
    real(kind=real64) :: axis_point(3)                               !! A point of the axis
    real(kind=real64) :: axis_vector(3)                              !! Axis
    real(kind=real64) :: angle                                       !! Rotation angle (radians)
    real(kind=real64) :: point(3)                                    !! Observation point
    fbem_rotation_infinitesimal_displacement=fbem_cross_product(axis_vector,point-axis_point)*angle
  end function fbem_rotation_infinitesimal_displacement

  ! =====================================
  ! Coordinates transformation matrices L
  ! =====================================

  !! <tt>v = L·v'</tt>
  !!
  !! <tt>v'</tt> is the vector in the old coordinate system.
  !! <tt>v</tt> is the vector in the new coordinate system.
  subroutine fbem_coordinate_transformation_L(n,ep,e,L)
    implicit none
    !! Number of dimensions: n=2 in 2D, n=3 in 3D
    integer :: n
    !! Unit axis vectors of x' system of coordinates: ep(i,j): component i of axis j
    real(kind=real64) :: ep(n,n)
    !! Unit axis vectors of x system of coordinates: e(i,j): component i of axis j
    real(kind=real64) :: e(n,n)
    !! Transformation matrix L
    real(kind=real64) :: L(n,n)
    !
    integer :: i, j, k
    ! Row
    do i=1,n
      ! Column
      do j=1,n
        ! Initialize
        L(i,j)=0.0d0
        ! Dot product
        do k=1,n
          L(i,j)=L(i,j)+e(k,i)*ep(k,j)
        end do
      end do
    end do
  end subroutine fbem_coordinate_transformation_L

  ! ================
  ! Position vectors
  ! ================

  !! Calculation of 2D position vector of an 1D element.
  function position2d_1d(etype,x_nodes,xi)
    implicit none
    ! I/O
    integer            :: etype                          !! Element geometrical interpolation type: <tt>etype={fbem_line2,fbem_line3}</tt>.
    real(kind=real64)  :: x_nodes(2,fbem_n_nodes(etype)) !! 2D position of element nodes
    real(kind=real64)  :: xi                             !! Reference space coordinate
    real(kind=real64)  :: position2d_1d(2)               !! 2D position vector of an 1D element.
    ! Local
    integer           :: k
    real(kind=real64) :: phi(fbem_n_nodes(etype))
    real(kind=real64) :: aux(10)
    real(kind=real64) :: x(2)
#   if CHECK_XIS==1
    call fbem_check_xi_warning(xi)
#   endif
#   if CHECK_XIS==2
    call fbem_check_xi_error(xi)
#   endif
#   define delta 0.0d0
#   include <phi_1d.rc>
#   undef delta
    ! Components calculation
    x=0.
    do k=1,fbem_n_nodes(etype)
        x=x+phi(k)*x_nodes(:,k)
    end do
    position2d_1d=x
  end function position2d_1d

  !! Calculation of 3D position vector of an 1D element.
  function position3d_1d(etype,x_nodes,xi)
    implicit none
    ! I/O
    integer            :: etype                          !! Element geometrical interpolation type: <tt>etype={fbem_line2,fbem_line3}</tt>.
    real(kind=real64)  :: x_nodes(3,fbem_n_nodes(etype)) !! 2D position of element nodes
    real(kind=real64)  :: xi                             !! Reference space coordinate
    real(kind=real64)  :: position3d_1d(3)               !! 3D position vector of an 1D element.
    ! Local
    integer           :: k
    real(kind=real64) :: phi(fbem_n_nodes(etype))
    real(kind=real64) :: aux(10)
    real(kind=real64) :: x(3)
#   if CHECK_XIS==1
    call fbem_check_xi_warning(xi)
#   endif
#   if CHECK_XIS==2
    call fbem_check_xi_error(xi)
#   endif
#   define delta 0.0d0
#   include <phi_1d.rc>
#   undef delta
    ! Components calculation
    x=0.
    do k=1,fbem_n_nodes(etype)
        x=x+phi(k)*x_nodes(:,k)
    end do
    position3d_1d=x
  end function position3d_1d

  !! Calculation of 2D position vector of an 2D element.
  function position2d_2d(etype,x_nodes,xi)
    implicit none
    ! I/O
    integer            :: etype                          !! Element geometrical interpolation type: <tt>etype={fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    real(kind=real64)  :: x_nodes(2,fbem_n_nodes(etype)) !! 2D position of element nodes
    real(kind=real64)  :: xi(2)                          !! Reference space coordinates
    real(kind=real64)  :: position2d_2d(2)               !! 2D position vector of an 2D element.
    ! Local
    integer           :: k
    real(kind=real64) :: phi(fbem_n_nodes(etype))
    real(kind=real64) :: aux(10)
    real(kind=real64) :: x(2)
#   if CHECK_XIS==1
    call fbem_check_xi1xi2_warning(etype,xi)
#   endif
#   if CHECK_XIS==2
    call fbem_check_xi1xi2_error(etype,xi)
#   endif
#   define delta 0.0d0
#   include <phi_2d.rc>
#   undef delta
    ! Components calculation
    x=0.
    do k=1,fbem_n_nodes(etype)
        x=x+phi(k)*x_nodes(:,k)
    end do
    position2d_2d=x
  end function position2d_2d

  !! Calculation of 3D position vector of an 2D or 3D element (only implemented for 2D elements).
  function position3d_nd(etype,x_nodes,xi)
    implicit none
    integer            :: etype                          !! Element geometrical interpolation type: <tt>etype={fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    real(kind=real64)  :: x_nodes(3,fbem_n_nodes(etype)) !! 2D position of element nodes
    real(kind=real64)  :: xi(fbem_n_dimension(etype))    !! Reference space coordinates
    real(kind=real64)  :: position3d_nd(3)               !! 3D position vector of an nD element.
    ! Local
    integer           :: k
    real(kind=real64) :: phi(fbem_n_nodes(etype))
    real(kind=real64) :: aux(10)
    real(kind=real64) :: x(3)
#   if CHECK_XIS==1
    call fbem_check_xi1xi2_warning(etype,xi)
#   endif
#   if CHECK_XIS==2
    call fbem_check_xi1xi2_error(etype,xi)
#   endif
#   define delta 0.0d0
#   include <phi_2d.rc>
#   undef delta
    ! Components calculation
    x=0.
    do k=1,fbem_n_nodes(etype)
      x=x+phi(k)*x_nodes(:,k)
    end do
    position3d_nd=x
  end function position3d_nd

  !! Calculation of position vector of an 1D or 2D element.
  function fbem_position(rn,etype,x_nodes,xi_p)
    implicit none
    ! I/O
    integer           :: rn                              !! Dimensional space
    integer           :: etype                           !! Element geometrical interpolation type: <tt>etype={line2,line3,tri3,tri6,quad4,quad8,quad9}</tt>.
    real(kind=real64) :: x_nodes(rn,fbem_n_nodes(etype)) !! Position of element geometrical nodes
    real(kind=real64) :: xi_p(fbem_n_dimension(etype))   !! Reference space coordinates
    real(kind=real64) :: fbem_position(rn)
    ! Local
    integer           :: k
    real(kind=real64) :: phi(fbem_n_nodes(etype))
    real(kind=real64) :: aux(10)
    real(kind=real64) :: x(rn)
#   if CHECK_XIS==1
    call fbem_check_xi1xi2_warning(etype,xi)
#   endif
#   if CHECK_XIS==2
    call fbem_check_xi1xi2_error(etype,xi)
#   endif
    select case (fbem_n_dimension(etype))
      case (1)
#       define xi xi_p(1)
#       define delta 0.d0
#       include <phi_1d.rc>
#       undef xi
#       undef delta
      case (2)
#       define xi xi_p
#       define delta 0.d0
#       include <phi_2d.rc>
#       undef xi
#       undef delta
      case (3)
#       define xi xi_p
#       define delta 0.d0
#       include <phi_3d.rc>
#       undef xi
#       undef delta
    end select
    ! Components calculation
    x=0.d0
    do k=1,fbem_n_nodes(etype)
        x=x+phi(k)*x_nodes(:,k)
    end do
    fbem_position=x
  end function fbem_position

  ! ===============
  ! Tangent vectors
  ! ===============

  !! Calculation of the tangent vector of a 1D element
  function fbem_tangent_xi(rn,etype,x_nodes,xi)
    implicit none
    integer           :: rn                             !! Dimensional space
    integer           :: etype                          !! Element geometrical interpolation type: <tt>etype={fbem_line2,fbem_line3}</tt>.
    real(kind=real64) :: x_nodes(rn,fbem_n_nodes(etype)) !! 3D position of element nodes
    real(kind=real64) :: xi                             !! Reference space coordinate
    real(kind=real64) :: fbem_tangent_xi(rn)
    real(kind=real64) :: dphidxi_vector(fbem_n_nodes(etype))
    real(kind=real64) :: aux(10)
    integer           :: i,j
#   if CHECK_XIS==1
    call fbem_check_xi_warning(xi)
#   endif
#   if CHECK_XIS==2
    call fbem_check_xi_error(xi)
#   endif
#   define dphidxi dphidxi_vector
#   define delta 0.0d0
#   include <dphidxi_1d.rc>
#   undef delta
#   undef dphidxi
    ! Components calculation
    do j=1,rn
      fbem_tangent_xi(j)=0.0d0
    end do
    do i=1,fbem_n_nodes(etype)
      do j=1,rn
        fbem_tangent_xi(j)=fbem_tangent_xi(j)+dphidxi_vector(i)*x_nodes(j,i)
      end do
    end do
  end function fbem_tangent_xi

  !! Calculation of the tangent vector of a 2D element in <tt>xi_1</tt> direction
  function fbem_tangent_xi1(rn,etype,x_nodes,xi)
    implicit none
    integer           :: rn                             !! Dimensional space
    integer           :: etype                          !! Element geometrical interpolation type: <tt>etype={fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    real(kind=real64) :: x_nodes(rn,fbem_n_nodes(etype)) !! 3D position of element nodes
    real(kind=real64) :: xi(2)                          !! Reference space coordinates
    real(kind=real64) :: fbem_tangent_xi1(rn)
    real(kind=real64) :: dphidxi1_vector(fbem_n_nodes(etype))
    integer           :: i,j
    real(kind=real64) :: aux(10)
#   if CHECK_XIS==1
    call fbem_check_xi1xi2_warning(etype,xi)
#   endif
#   if CHECK_XIS==2
    call fbem_check_xi1xi2_error(etype,xi)
#   endif
#   define dphidxi1 dphidxi1_vector
#   define delta 0.0d0
#   include <dphidxi1_2d.rc>
#   undef delta
#   undef dphidxi1
    ! Components calculation
    do j=1,rn
      fbem_tangent_xi1(j)=0.0d0
    end do
    do i=1,fbem_n_nodes(etype)
      do j=1,rn
        fbem_tangent_xi1(j)=fbem_tangent_xi1(j)+dphidxi1_vector(i)*x_nodes(j,i)
      end do
    end do
  end function fbem_tangent_xi1

  !! Calculation of the tangent vector of a 2D element in <tt>xi_2</tt> direction
  function fbem_tangent_xi2(rn,etype,x_nodes,xi)
    implicit none
    integer           :: rn                             !! Dimensional space
    integer           :: etype                          !! Element geometrical interpolation type: <tt>etype={fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    real(kind=real64) :: x_nodes(rn,fbem_n_nodes(etype)) !! 3D position of element nodes
    real(kind=real64) :: xi(2)                          !! Reference space coordinates
    real(kind=real64) :: fbem_tangent_xi2(rn)
    real(kind=real64) :: dphidxi2_vector(fbem_n_nodes(etype))
    integer           :: i,j
    real(kind=real64) :: aux(10)
#   if CHECK_XIS==1
    call fbem_check_xi1xi2_warning(etype,xi)
#   endif
#   if CHECK_XIS==2
    call fbem_check_xi1xi2_error(etype,xi)
#   endif
#   define dphidxi2 dphidxi2_vector
#   define delta 0.0d0
#   include <dphidxi2_2d.rc>
#   undef delta
#   undef dphidxi2
    ! Components calculation
    do j=1,rn
      fbem_tangent_xi2(j)=0.0d0
    end do
    do i=1,fbem_n_nodes(etype)
      do j=1,rn
        fbem_tangent_xi2(j)=fbem_tangent_xi2(j)+dphidxi2_vector(i)*x_nodes(j,i)
      end do
    end do
  end function fbem_tangent_xi2

  !! Calculation of the unit tangent vector of a 1D element
  function fbem_utangent_xi(rn,etype,x_nodes,xi)
    implicit none
    integer           :: rn                             !! Dimensional space
    integer           :: etype                          !! Element geometrical interpolation type: <tt>etype={fbem_line2,fbem_line3}</tt>.
    real(kind=real64) :: x_nodes(rn,fbem_n_nodes(etype)) !! 3D position of element nodes
    real(kind=real64) :: xi                             !! Reference space coordinate
    real(kind=real64) :: fbem_utangent_xi(rn)
    real(kind=real64) :: dphidxi_vector(fbem_n_nodes(etype))
    real(kind=real64) :: aux(10), norm
    integer           :: i,j
#   if CHECK_XIS==1
    call fbem_check_xi_warning(xi)
#   endif
#   if CHECK_XIS==2
    call fbem_check_xi_error(xi)
#   endif
#   define dphidxi dphidxi_vector
#   define delta 0.0d0
#   include <dphidxi_1d.rc>
#   undef delta
#   undef dphidxi
    ! Components calculation
    do j=1,rn
      fbem_utangent_xi(j)=0.0d0
    end do
    do i=1,fbem_n_nodes(etype)
      do j=1,rn
        fbem_utangent_xi(j)=fbem_utangent_xi(j)+dphidxi_vector(i)*x_nodes(j,i)
      end do
    end do
    ! Norm
    norm=0.0d0
    do i=1,rn
      norm=norm+fbem_utangent_xi(i)**2
    end do
    norm=dsqrt(norm)
    ! Normalization
    do i=1,rn
      fbem_utangent_xi(i)=fbem_utangent_xi(i)/norm
    end do
  end function fbem_utangent_xi

  !! Calculation of the unit tangent vector of a 2D element in <tt>xi_1</tt> direction
  function fbem_utangent_xi1(rn,etype,x_nodes,xi)
    implicit none
    integer           :: rn                             !! Dimensional space
    integer           :: etype                          !! Element geometrical interpolation type: <tt>etype={fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    real(kind=real64) :: x_nodes(rn,fbem_n_nodes(etype)) !! 3D position of element nodes
    real(kind=real64) :: xi(2)                          !! Reference space coordinates
    real(kind=real64) :: fbem_utangent_xi1(rn)
    real(kind=real64) :: dphidxi1_vector(fbem_n_nodes(etype))
    integer           :: i,j
    real(kind=real64) :: aux(10), norm
#   if CHECK_XIS==1
    call fbem_check_xi1xi2_warning(etype,xi)
#   endif
#   if CHECK_XIS==2
    call fbem_check_xi1xi2_error(etype,xi)
#   endif
#   define dphidxi1 dphidxi1_vector
#   define delta 0.0d0
#   include <dphidxi1_2d.rc>
#   undef delta
#   undef dphidxi1
    ! Components calculation
    do j=1,rn
      fbem_utangent_xi1(j)=0.0d0
    end do
    do i=1,fbem_n_nodes(etype)
      do j=1,rn
        fbem_utangent_xi1(j)=fbem_utangent_xi1(j)+dphidxi1_vector(i)*x_nodes(j,i)
      end do
    end do
    ! Norm
    norm=0.0d0
    do i=1,rn
      norm=norm+fbem_utangent_xi1(i)**2
    end do
    norm=dsqrt(norm)
    ! Normalization
    do i=1,rn
      fbem_utangent_xi1(i)=fbem_utangent_xi1(i)/norm
    end do
  end function fbem_utangent_xi1

  !! Calculation of the unit tangent vector of a 2D element in <tt>xi_2</tt> direction
  function fbem_utangent_xi2(rn,etype,x_nodes,xi)
    implicit none
    integer           :: rn                             !! Dimensional space
    integer           :: etype                          !! Element geometrical interpolation type: <tt>etype={fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    real(kind=real64) :: x_nodes(rn,fbem_n_nodes(etype)) !! 3D position of element nodes
    real(kind=real64) :: xi(2)                          !! Reference space coordinates
    real(kind=real64) :: fbem_utangent_xi2(rn)
    real(kind=real64) :: dphidxi2_vector(fbem_n_nodes(etype))
    integer           :: i,j
    real(kind=real64) :: aux(10), norm
#   if CHECK_XIS==1
    call fbem_check_xi1xi2_warning(etype,xi)
#   endif
#   if CHECK_XIS==2
    call fbem_check_xi1xi2_error(etype,xi)
#   endif
#   define dphidxi2 dphidxi2_vector
#   define delta 0.0d0
#   include <dphidxi2_2d.rc>
#   undef delta
#   undef dphidxi2
    ! Components calculation
    do j=1,rn
      fbem_utangent_xi2(j)=0.0d0
    end do
    do i=1,fbem_n_nodes(etype)
      do j=1,rn
        fbem_utangent_xi2(j)=fbem_utangent_xi2(j)+dphidxi2_vector(i)*x_nodes(j,i)
      end do
    end do
    ! Norm
    norm=0.0d0
    do i=1,rn
      norm=norm+fbem_utangent_xi2(i)**2
    end do
    norm=dsqrt(norm)
    ! Normalization
    do i=1,rn
      fbem_utangent_xi2(i)=fbem_utangent_xi2(i)/norm
    end do
  end function fbem_utangent_xi2

  subroutine fbem_tangent_derivatives_2d(rn,etype,x_nodes,xi,dTidxij)
    implicit none
    ! I/O
    integer           :: rn                              !! Dimensional space
    integer           :: etype                           !! Element geometrical interpolation type: <tt>etype={fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    real(kind=real64) :: x_nodes(rn,fbem_n_nodes(etype)) !! 3D position of element nodes
    real(kind=real64) :: xi(2)                           !! Reference space coordinates
    real(kind=real64) :: dTidxij(rn,2,2)
    ! Local
    real(kind=real64) :: d2phidxi(fbem_n_nodes(etype),2,2)
    integer           :: i, j, n
    real(kind=real64) :: aux(10)
#   define delta 0.0d0
#   include <d2phidxi_2d.rc>
#   undef delta
    n=fbem_n_nodes(etype)
    dTidxij=0.d0
    do i=1,n
      dTidxij(:,1,1)=dTidxij(:,1,1)+d2phidxi(i,1,1)*x_nodes(:,i)
      dTidxij(:,1,2)=dTidxij(:,1,2)+d2phidxi(i,1,2)*x_nodes(:,i)
      dTidxij(:,2,2)=dTidxij(:,2,2)+d2phidxi(i,2,2)*x_nodes(:,i)
    end do
    dTidxij(:,2,1)=dTidxij(:,1,2)
  end subroutine fbem_tangent_derivatives_2d

  subroutine fbem_utangents_at_boundary(rn,etype,x_nodes,node,tbp,tbm)
    implicit none
    ! I/O variables
    integer           :: rn
    integer           :: etype
    real(kind=real64) :: x_nodes(rn,fbem_n_nodes(etype))
    integer           :: node
    real(kind=real64) :: tbp(rn)
    real(kind=real64) :: tbm(rn)
    ! Local variables
    integer           :: kn, kc, nnodes
    real(kind=real64) :: xiaux(fbem_n_dimension(etype))
    real(kind=real64) :: dphidxi(fbem_n_nodes(etype))
    real(kind=real64) :: dphidxi1(fbem_n_nodes(etype))
    real(kind=real64) :: dphidxi2(fbem_n_nodes(etype))
    real(kind=real64) :: dphidxi3(fbem_n_nodes(etype))
    real(kind=real64) :: aux(10)
    real(kind=real64) :: t(rn)
    real(kind=real64) :: t1(rn), t2(rn), t3(rn)
    real(kind=real64) :: norm
    !
    nnodes=fbem_n_nodes(etype)
    !
    select case (fbem_n_dimension(etype))
      case (1)
        ! Save xi
#       define xi xiaux(1)
#       define delta 0.0d0
#       include <xi_1d_at_node.rc>
#       undef xi
#       undef delta
        ! Calculate the shape functions derivatives at xi
#       define xi xiaux(1)
#       define delta 0.0d0
#       include <dphidxi_1d.rc>
#       undef xi
#       undef delta
        ! Calculate T
        do kc=1,rn
          T(kc)=0.0d0
          do kn=1,nnodes
            T(kc)=T(kc)+dphidxi(kn)*x_nodes(kc,kn)
          end do
        end do
        ! Normalize T to t
        norm=0.0d0
        do kc=1,rn
          norm=norm+T(kc)**2
        end do
        norm=dsqrt(norm)
        do kc=1,rn
          t(kc)=T(kc)/norm
        end do
        !
        ! The tb depend on the node
        !
        select case (node)
          ! Node 1 (vertex node)
          case (1)
            do kc=1,rn
              tbp(kc)=t(kc)
              tbm(kc)=t(kc)
            end do
          ! Node 2 (vertex node)
          case (2)
            do kc=1,rn
              tbp(kc)=-t(kc)
              tbm(kc)=-t(kc)
            end do
          ! For the other nodes they are not defined
          case default
            do kc=1,rn
              tbp(kc)=0.0d0
              tbm(kc)=0.0d0
            end do
        end select
      case (2)
        ! Save xi1 and xi2
#       define xi xiaux
#       define delta 0.0d0
#       include <xi_2d_at_node.rc>
#       undef xi
#       undef delta
        ! Calculate the shape functions derivatives at xi
#       define xi xiaux
#       define delta 0.0d0
#       include <dphidxik_2d.rc>
#       undef delta
#       undef xi
        ! Calculate T1 and T2
        do kc=1,rn
          T1(kc)=0.0d0
          T2(kc)=0.0d0
          do kn=1,nnodes
            T1(kc)=T1(kc)+dphidxi1(kn)*x_nodes(kc,kn)
            T2(kc)=T2(kc)+dphidxi2(kn)*x_nodes(kc,kn)
          end do
        end do
        ! Normalize T1 to t1
        norm=0.0d0
        do kc=1,rn
          norm=norm+T1(kc)**2
        end do
        norm=dsqrt(norm)
        do kc=1,rn
          t1(kc)=T1(kc)/norm
        end do
        ! Normalize T2 to t2
        norm=0.0d0
        do kc=1,rn
          norm=norm+T2(kc)**2
        end do
        norm=dsqrt(norm)
        do kc=1,rn
          t2(kc)=T2(kc)/norm
        end do
        ! The tb depend on the type of element
        select case (etype)
          !
          ! Triangular
          !
          case (fbem_tri3,fbem_tri6)
            !
            ! For nodes at the edge 1, i.e. edge between nodes 1 and 2. it is necessary an auxiliary T3 tangent vector.
            !
            ! It is used 1D interpolation of line2 (for tri3) and line3 (for tri6) in the [0,1] domain.
            select case (etype)
              case (fbem_tri3)
                dphidxi3(1)= 1.0d0
                dphidxi3(2)=-1.0d0
                dphidxi3(3)= 0.0d0
              case (fbem_tri6)
                dphidxi3(1)=4.0d0*xiaux(1)-1.0d0
                dphidxi3(2)=4.0d0*xiaux(1)-3.0d0
                dphidxi3(3)=0.0d0
                dphidxi3(4)=4.0d0*(1.0d0-2.0d0*xiaux(1))
                dphidxi3(5)=0.0d0
                dphidxi3(6)=0.0d0
            end select
            ! Calculate T3
            do kc=1,rn
              T3(kc)=0.0d0
              do kn=1,nnodes
                T3(kc)=T3(kc)+dphidxi3(kn)*x_nodes(kc,kn)
              end do
            end do
            ! Normalize T3 to t3
            norm=0.0d0
            do kc=1,rn
              norm=norm+T3(kc)**2
            end do
            norm=dsqrt(norm)
            do kc=1,rn
              t3(kc)=T3(kc)/norm
            end do
            ! The tb depend on the node of element
            select case (node)
              ! Node 1 (vertex node)
              case (1)
                do kc=1,rn
                  tbp(kc)=-t3(kc)
                  tbm(kc)=-t1(kc)
                end do
              ! Node 2 (vertex node)
              case (2)
                do kc=1,rn
                  tbp(kc)=-t2(kc)
                  tbm(kc)= t3(kc)
                end do
              ! Node 3 (vertex node)
              case (3)
                do kc=1,rn
                  tbp(kc)= t1(kc)
                  tbm(kc)= t2(kc)
                end do
              ! Node 4 (edge node)
              case (4)
                do kc=1,rn
                  tbp(kc)=-t3(kc)
                  tbm(kc)= t3(kc)
                end do
              ! Node 5 (edge node)
              case (5)
                do kc=1,rn
                  tbp(kc)=-t2(kc)
                  tbm(kc)= t2(kc)
                end do
              ! Node 6 (edge node)
              case (6)
                do kc=1,rn
                  tbp(kc)= t1(kc)
                  tbm(kc)=-t1(kc)
                end do
            end select
          !
          ! Quadrilateral
          !
          case (fbem_quad4,fbem_quad8,fbem_quad9)
            ! The tb depend on the node of element
            select case (node)
              ! Node 1 (vertex node)
              case (1)
                do kc=1,rn
                  tbp(kc)= t1(kc)
                  tbm(kc)= t2(kc)
                end do
              ! Node 2 (vertex node)
              case (2)
                do kc=1,rn
                  tbp(kc)= t2(kc)
                  tbm(kc)=-t1(kc)
                end do
              ! Node 3 (vertex node)
              case (3)
                do kc=1,rn
                  tbp(kc)=-t1(kc)
                  tbm(kc)=-t2(kc)
                end do
              ! Node 4 (vertex node)
              case (4)
                do kc=1,rn
                  tbp(kc)=-t2(kc)
                  tbm(kc)= t1(kc)
                end do
             ! Node 5 (edge node)
              case (5)
                do kc=1,rn
                  tbp(kc)= t1(kc)
                  tbm(kc)=-t1(kc)
                end do
             ! Node 6 (edge node)
              case (6)
                do kc=1,rn
                  tbp(kc)= t2(kc)
                  tbm(kc)=-t2(kc)
                end do
             ! Node 7 (edge node)
              case (7)
                do kc=1,rn
                  tbp(kc)=-t1(kc)
                  tbm(kc)= t1(kc)
                end do
             ! Node 8 (edge node)
              case (8)
                do kc=1,rn
                  tbp(kc)=-t2(kc)
                  tbm(kc)= t2(kc)
                end do
              ! Node 9 (interior node), not defined
              case (9)
                do kc=1,rn
                  tbp(kc)=0.0d0
                  tbm(kc)=0.0d0
                end do
            end select
        end select
    end select
  end subroutine

  ! ==============
  ! Normal vectors
  ! ==============

  !! Calculation 2D normal vector of a 1D element (<tt> N = T x e_3</tt>)
  function fbem_normal2d(etype,x_nodes,xi)
    implicit none
    !! Element geometrical interpolation type: <tt>etype={fbem_line2,fbem_line3}</tt>.
    integer :: etype
    !! 2D position of element nodes
    real(kind=real64) :: x_nodes(2,fbem_n_nodes(etype))
    !! Reference space coordinate
    real(kind=real64) :: xi
    real(kind=real64) :: fbem_normal2d(2)
    real(kind=real64) :: dphidxi_vector(fbem_n_nodes(etype))
    real(kind=real64) :: aux(10)
    real(kind=real64) :: tmp
    integer :: i,j
#   if CHECK_XIS==1
    call fbem_check_xi_warning(xi)
#   endif
#   if CHECK_XIS==2
    call fbem_check_xi_error(xi)
#   endif
#   define dphidxi dphidxi_vector
#   define delta 0.0d0
#   include <dphidxi_1d.rc>
#   undef delta
#   undef dphidxi
    ! Components of tangent vector calculation
    do j=1,2
      fbem_normal2d(j)=0.0d0
    end do
    do i=1,fbem_n_nodes(etype)
      do j=1,2
        fbem_normal2d(j)=fbem_normal2d(j)+dphidxi_vector(i)*x_nodes(j,i)
      end do
    end do
    ! Normal vector
    tmp=fbem_normal2d(1)
    fbem_normal2d(1)=fbem_normal2d(2)
    fbem_normal2d(2)=-tmp
  end function fbem_normal2d

  !! Calculation 2D unit normal vector of a 1D element (<tt> n = T x e_3/|T x e_3|</tt>)
  function fbem_unormal2d(etype,x_nodes,xi)
    implicit none
    !! Element geometrical interpolation type: <tt>etype={fbem_line2,fbem_line3}</tt>.
    integer :: etype
    !! 2D position of element nodes
    real(kind=real64) :: x_nodes(2,fbem_n_nodes(etype))
    !! Reference space coordinate
    real(kind=real64) :: xi
    real(kind=real64) :: fbem_unormal2d(2)
    real(kind=real64) :: dphidxi_vector(fbem_n_nodes(etype))
    real(kind=real64) :: aux(10)
    real(kind=real64) :: tmp
    integer :: i,j
#   if CHECK_XIS==1
    call fbem_check_xi_warning(xi)
#   endif
#   if CHECK_XIS==2
    call fbem_check_xi_error(xi)
#   endif
#   define dphidxi dphidxi_vector
#   define delta 0.0d0
#   include <dphidxi_1d.rc>
#   undef delta
#   undef dphidxi
    ! Components of tangent vector calculation
    fbem_unormal2d=0.0d0
    do i=1,fbem_n_nodes(etype)
      fbem_unormal2d=fbem_unormal2d+dphidxi_vector(i)*x_nodes(:,i)
    end do
    ! Normal vector
    tmp=fbem_unormal2d(1)
    fbem_unormal2d(1)=fbem_unormal2d(2)
    fbem_unormal2d(2)=-tmp
    ! Normalize
    fbem_unormal2d=fbem_unormal2d/sqrt(dot_product(fbem_unormal2d,fbem_unormal2d))
  end function fbem_unormal2d

  !! Calculation 3D normal vector of a 2D element
  function fbem_normal3d(etype,x_nodes,xi)
    implicit none
    !! Element geometrical interpolation type: <tt>etype={fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    integer :: etype
    !! 3D position of element nodes
    real(kind=real64) :: x_nodes(3,fbem_n_nodes(etype))
    !! Reference space coordinates
    real(kind=real64) :: xi(2)
    real(kind=real64) :: fbem_normal3d(3)
    real(kind=real64) :: tangent_xi1(3)
    real(kind=real64) :: tangent_xi2(3)
    real(kind=real64) :: dphidxi1_vector(fbem_n_nodes(etype)), dphidxi2_vector(fbem_n_nodes(etype))
    integer :: i,j
    real(kind=real64) :: aux(10)
#   if CHECK_XIS==1
    call fbem_check_xi1xi2_warning(etype,xi)
#   endif
#   if CHECK_XIS==2
    call fbem_check_xi1xi2_error(etype,xi)
#   endif
#   define dphidxi1 dphidxi1_vector
#   define dphidxi2 dphidxi2_vector
#   define delta 0.0d0
#   include <dphidxik_2d.rc>
#   undef dphidxi1
#   undef dphidxi2
#   undef delta
    ! Components calculation of T1 and T2
    do j=1,3
      tangent_xi1(j)=0.0d0
      tangent_xi2(j)=0.0d0
    end do
    do i=1,fbem_n_nodes(etype)
      do j=1,3
        tangent_xi1(j)=tangent_xi1(j)+dphidxi1_vector(i)*x_nodes(j,i)
        tangent_xi2(j)=tangent_xi2(j)+dphidxi2_vector(i)*x_nodes(j,i)
      end do
    end do
    ! Normal vector as T1 x T2
    fbem_normal3d(1)=tangent_xi1(2)*tangent_xi2(3)-tangent_xi1(3)*tangent_xi2(2)
    fbem_normal3d(2)=tangent_xi1(3)*tangent_xi2(1)-tangent_xi1(1)*tangent_xi2(3)
    fbem_normal3d(3)=tangent_xi1(1)*tangent_xi2(2)-tangent_xi1(2)*tangent_xi2(1)
  end function fbem_normal3d

  !! Calculation 3D normal vector of a 2D element
  function fbem_unormal3d(etype,x_nodes,xi)
    implicit none
    !! Element geometrical interpolation type: <tt>etype={fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    integer :: etype
    !! 3D position of element nodes
    real(kind=real64) :: x_nodes(3,fbem_n_nodes(etype))
    !! Reference space coordinates
    real(kind=real64) :: xi(2)
    real(kind=real64) :: fbem_unormal3d(3)
    real(kind=real64) :: tangent_xi1(3)
    real(kind=real64) :: tangent_xi2(3)
    real(kind=real64) :: dphidxi1_vector(fbem_n_nodes(etype)), dphidxi2_vector(fbem_n_nodes(etype))
    integer :: i,j
    real(kind=real64) :: aux(10)
    real(kind=real64) :: tmp
#   if CHECK_XIS==1
    call fbem_check_xi1xi2_warning(etype,xi)
#   endif
#   if CHECK_XIS==2
    call fbem_check_xi1xi2_error(etype,xi)
#   endif
#   define dphidxi1 dphidxi1_vector
#   define dphidxi2 dphidxi2_vector
#   define delta 0.0d0
#   include <dphidxik_2d.rc>
#   undef dphidxi1
#   undef dphidxi2
#   undef delta
    ! Components calculation of T1 and T2
    do j=1,3
      tangent_xi1(j)=0.0d0
      tangent_xi2(j)=0.0d0
    end do
    do i=1,fbem_n_nodes(etype)
      do j=1,3
        tangent_xi1(j)=tangent_xi1(j)+dphidxi1_vector(i)*x_nodes(j,i)
        tangent_xi2(j)=tangent_xi2(j)+dphidxi2_vector(i)*x_nodes(j,i)
      end do
    end do
    ! Normal vector as T1 x T2
    fbem_unormal3d(1)=tangent_xi1(2)*tangent_xi2(3)-tangent_xi1(3)*tangent_xi2(2)
    fbem_unormal3d(2)=tangent_xi1(3)*tangent_xi2(1)-tangent_xi1(1)*tangent_xi2(3)
    fbem_unormal3d(3)=tangent_xi1(1)*tangent_xi2(2)-tangent_xi1(2)*tangent_xi2(1)
    ! Normalize
    tmp=dsqrt(fbem_unormal3d(1)**2+fbem_unormal3d(2)**2+fbem_unormal3d(3)**2)
    fbem_unormal3d(1)=fbem_unormal3d(1)/tmp
    fbem_unormal3d(2)=fbem_unormal3d(2)/tmp
    fbem_unormal3d(3)=fbem_unormal3d(3)/tmp
  end function fbem_unormal3d

  !! Calculation of 2D normal vector of a 1D element and tangent to center (only valid for vertex nodes)
  subroutine fbem_normal2d_and_tangent_to_center(etype,x_nodes,node,n,tb)
    implicit none
    !! Element geometrical interpolation type: <tt>etype={fbem_line2,fbem_line3}</tt>.
    integer :: etype
    !! 2D position of element nodes
    real(kind=real64) :: x_nodes(2,fbem_n_nodes(etype))
    !! Local node
    integer :: node
    !! Normal vector
    real(kind=real64) :: n(2)
    !! Tangent vector at the node with angle bisector direction
    real(kind=real64) :: tb(2)
    real(kind=real64) :: xi
    real(kind=real64) :: T(2)
    real(kind=real64) :: dphidxi_vector(fbem_n_nodes(etype))
    integer :: i,j
    real(kind=real64) :: norm, aux(10)
    ! If the node is interior tb is null
    if (node.gt.2) then
      do i=1,2
        tb(i)=0.0d0
      end do
    else
#     if CHECK_XIS==1
      call fbem_check_xi_warning(xi)
#     endif
#     if CHECK_XIS==2
      call fbem_check_xi_error(xi)
#     endif
#     define delta 0.0d0
#     include <xi_1d_at_node.rc>
#     define dphidxi dphidxi_vector
#     include <dphidxi_1d.rc>
#     undef dphidxi
#     undef delta
      ! Components calculation of T
      do j=1,2
        T(j)=0.0d0
      end do
      do i=1,fbem_n_nodes(etype)
        do j=1,2
          T(j)=T(j)+dphidxi_vector(i)*x_nodes(j,i)
        end do
      end do
      ! Normal vector
      n(1)=T(2)
      n(2)=-T(1)
      ! Tangent vector
      tb(1)=T(1)
      tb(2)=T(2)
      ! Normalization of all vectors
      ! Norm of normal
      norm=dsqrt(n(1)**2+n(2)**2)
      ! Unit normal
      do i=1,2
        n(i)=n(i)/norm
      end do
      ! Norm of tangent
      norm=dsqrt(tb(1)**2+tb(2)**2)
      ! Unit tangent
      do i=1,2
        tb(i)=tb(i)/norm
      end do
      ! If collocation point is at node 2, the tangent has to be inverted
      if (node.eq.2) then
        do i=1,2
          tb(i)=-tb(i)
        end do
      end if
    end if
  end subroutine fbem_normal2d_and_tangent_to_center

  !! Calculation of 3D normal vector of a 2D element and bisector tangent (only valid for edge or vertex nodes)
  subroutine fbem_normal3d_and_bis_tangent_at_node(etype,x_nodes,node,n,tb)
    implicit none
    !! Element geometrical interpolation type: <tt>etype={fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    integer :: etype
    !! 3D position of element nodes
    real(kind=real64) :: x_nodes(3,fbem_n_nodes(etype))
    !! Local node
    integer :: node
    !! Normal vector
    real(kind=real64) :: n(3)
    !! Tangent vector at the node with angle bisector direction
    real(kind=real64) :: tb(3)
    real(kind=real64) :: xi(2)
    real(kind=real64) :: T1(3)
    real(kind=real64) :: T2(3)
    real(kind=real64) :: T3(3)
    real(kind=real64) :: dphidxi1_vector(fbem_n_nodes(etype))
    real(kind=real64) :: dphidxi2_vector(fbem_n_nodes(etype))
    real(kind=real64) :: dphidxi3_vector(fbem_n_nodes(etype))
    integer :: i,j
    real(kind=real64) :: norm, aux(10)
#   if CHECK_XIS==1
    call fbem_check_xi1xi2_warning(etype,xi)
#   endif
#   if CHECK_XIS==2
    call fbem_check_xi1xi2_error(etype,xi)
#   endif
#   define delta 0.0d0
#   include <xi_2d_at_node.rc>
#   define dphidxi1 dphidxi1_vector
#   define dphidxi2 dphidxi2_vector
#   include <dphidxik_2d.rc>
#   undef dphidxi1
#   undef dphidxi2
#   undef delta
    ! Components calculation of T1 and T2
    do j=1,3
      T1(j)=0.0d0
      T2(j)=0.0d0
    end do
    do i=1,fbem_n_nodes(etype)
      do j=1,3
        T1(j)=T1(j)+dphidxi1_vector(i)*x_nodes(j,i)
        T2(j)=T2(j)+dphidxi2_vector(i)*x_nodes(j,i)
      end do
    end do
    ! Normal vector as T1 x T2
    N(1)=T1(2)*T2(3)-T1(3)*T2(2)
    N(2)=T1(3)*T2(1)-T1(1)*T2(3)
    N(3)=T1(1)*T2(2)-T1(2)*T2(1)
    ! Normalization of all vectors
    ! Norm of normal
    norm=dsqrt(N(1)**2+N(2)**2+N(3)**2)
    ! Unit normal
    do i=1,3
      n(i)=N(i)/norm
    end do
    ! Norm of tangent 1
    norm=dsqrt(T1(1)**2+T1(2)**2+T1(3)**2)
    ! Unit tangent 1
    do i=1,3
      t1(i)=T1(i)/norm
    end do
    ! Norm of tangent 2
    norm=dsqrt(T2(1)**2+T2(2)**2+T2(3)**2)
    ! Unit tangent 2
    do i=1,3
      t2(i)=T2(i)/norm
    end do
    ! Tangent vector at the node with angle bisector direction
    select case (etype)
      case (fbem_tri3,fbem_tri6)
        ! Calculation of T3 (tangent vector of edge 1 of triangles)
        ! It is used 1D interpolation of line2 (for tri3) and line3 (for tri6) in the [0,1] domain
        select case (etype)
          case (fbem_tri3)
            dphidxi3_vector(1)=1.0d0
            dphidxi3_vector(2)=-1.0d0
            dphidxi3_vector(3)=0.0d0
          case (fbem_tri6)
            dphidxi3_vector(1)=4.0d0*xi(1)-1.0d0
            dphidxi3_vector(2)=4.0d0*xi(1)-3.0d0
            dphidxi3_vector(3)=0.0d0
            dphidxi3_vector(4)=4.0d0*(1.0d0-2.0d0*xi(1))
            dphidxi3_vector(5)=0.0d0
            dphidxi3_vector(6)=0.0d0
        end select
        ! Components calculation of T3
        do j=1,3
          T3(j)=0.0d0
        end do
        do i=1,fbem_n_nodes(etype)
          do j=1,3
            T3(j)=T3(j)+dphidxi3_vector(i)*x_nodes(j,i)
          end do
        end do
        ! Norm of tangent 3
        norm=dsqrt(T3(1)**2+T3(2)**2+T3(3)**2)
        ! Unit tangent 3
        do i=1,3
          t3(i)=T3(i)/norm
        end do
        ! Now, we can calculate angle bisector tangents
        select case (node)
          case (1)
            do i=1,3
              TB(i)=-t1(i)-t3(i)
            end do
          case (2)
            do i=1,3
              TB(i)=t3(i)-t2(i)
            end do
          case (3)
            do i=1,3
              TB(i)=t1(i)+t2(i)
            end do
          case (4)
            ! modo 1
            !do i=1,3
            !  TB(i)=-t1(i)-t2(i)
            !end do
            ! modo 2 (t3 x n)
            TB(1)=t3(2)*n(3)-t3(3)*n(2)
            TB(2)=t3(3)*n(1)-t3(1)*n(3)
            TB(3)=t3(1)*n(2)-t3(2)*n(1)
          case (5)
            ! modo 1
            !do i=1,3
            !  tb(i)=T1(i)
            !end do
            ! modo 2 (t2 x n)
            TB(1)=t2(2)*n(3)-t2(3)*n(2)
            TB(2)=t2(3)*n(1)-t2(1)*n(3)
            TB(3)=t2(1)*n(2)-t2(2)*n(1)
          case (6)
            ! modo 1
            do i=1,3
              TB(i)=t2(i)
            end do
            ! modo 2 (n x t1)
            TB(1)=n(2)*t1(3)-n(3)*t1(2)
            TB(2)=n(3)*t1(1)-n(1)*t1(3)
            TB(3)=n(1)*t1(2)-n(2)*t1(1)
        end select
      case (fbem_quad4,fbem_quad8,fbem_quad9)
        select case (node)
          case (1)
            do i=1,3
              TB(i)=t1(i)+t2(i)
            end do
          case (2)
            do i=1,3
              TB(i)=-t1(i)+t2(i)
            end do
          case (3)
            do i=1,3
              TB(i)=-t1(i)-t2(i)
            end do
          case (4)
            do i=1,3
              TB(i)=t1(i)-t2(i)
            end do
          case (5)
            ! modo 1
            !do i=1,3
            !  TB(i)=t2(i)
            !end do
            ! modo 2 (n x t1)
            TB(1)=n(2)*t1(3)-n(3)*t1(2)
            TB(2)=n(3)*t1(1)-n(1)*t1(3)
            TB(3)=n(1)*t1(2)-n(2)*t1(1)
          case (6)
            ! modo 1
            !do i=1,3
            !  TB(i)=-t1(i)
            !end do
            ! modo 2 (n x t2)
            TB(1)=n(2)*t2(3)-n(3)*t2(2)
            TB(2)=n(3)*t2(1)-n(1)*t2(3)
            TB(3)=n(1)*t2(2)-n(2)*t2(1)
          case (7)
            ! modo 1
            !do i=1,3
            !  TB(i)=-t2(i)
            !end do
            ! modo 2 (t1 x n)
            TB(1)=t1(2)*n(3)-t1(3)*n(2)
            TB(2)=t1(3)*n(1)-t1(1)*n(3)
            TB(3)=t1(1)*n(2)-t1(2)*n(1)
          case (8)
            ! modo 1
            !do i=1,3
            !  TB(i)=t1(i)
            !end do
            ! modo 2 (t2 x n)
            TB(1)=t2(2)*n(3)-t2(3)*n(2)
            TB(2)=t2(3)*n(1)-t2(1)*n(3)
            TB(3)=t2(1)*n(2)-t2(2)*n(1)
          ! In this case, the node is interior, and has no bisector tangent, 0.0d0 is given.
          case (9)
            do i=1,3
              TB(i)=0.0d0
            end do
        end select
    end select
    ! Norm of bisector tangent
    norm=dsqrt(TB(1)**2+TB(2)**2+TB(3)**2)
    ! Unit bisector tangent
    do i=1,3
      tb(i)=TB(i)/norm
    end do
  end subroutine fbem_normal3d_and_bis_tangent_at_node


  ! =========
  ! Jacobians
  ! =========

  !! Calculation 2D jacobian of a 1D element (module of the tangent vector)
  function fbem_jacobian2d_1d(etype,x_nodes,xi)
    implicit none
    !! Element geometrical interpolation type: <tt>etype={fbem_line2,fbem_line3}</tt>.
    integer :: etype
    !! 2D position of element nodes
    real(kind=real64) :: x_nodes(2,fbem_n_nodes(etype))
    !! Reference space coordinate
    real(kind=real64) :: xi
    real(kind=real64) :: fbem_jacobian2d_1d
    real(kind=real64) :: dphidxi_vector(fbem_n_nodes(etype))
    real(kind=real64) :: tangent(2)
    real(kind=real64) :: aux(10)
    integer :: i,j
#   if CHECK_XIS==1
    call fbem_check_xi_warning(xi)
#   endif
#   if CHECK_XIS==2
    call fbem_check_xi_error(xi)
#   endif
#   define dphidxi dphidxi_vector
#   define delta 0.0d0
#   include <dphidxi_1d.rc>
#   undef delta
#   undef dphidxi
    ! Components of tangent vector calculation
    do j=1,2
      tangent(j)=0.0d0
    end do
    do i=1,fbem_n_nodes(etype)
      do j=1,2
        tangent(j)=tangent(j)+dphidxi_vector(i)*x_nodes(j,i)
      end do
    end do
    ! Jacobian as the module of the tangent vector
    fbem_jacobian2d_1d=dsqrt(tangent(1)**2+tangent(2)**2)
  end function fbem_jacobian2d_1d

  !! Calculation of 2D jacobian of a 2D element (usual jacobian)
  function fbem_jacobian2d_2d(etype,x_nodes,xi)
    implicit none
    !! Element geometrical interpolation type: <tt>etype={fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    integer :: etype
    !! 2D position of element nodes
    real(kind=real64) :: x_nodes(2,fbem_n_nodes(etype))
    !! Reference space coordinates
    real(kind=real64) :: xi(2)
    real(kind=real64) :: fbem_jacobian2d_2d
    real(kind=real64) :: dxdxi1(2)
    real(kind=real64) :: dxdxi2(2)
    real(kind=real64) :: dphidxi1_vector(fbem_n_nodes(etype)), dphidxi2_vector(fbem_n_nodes(etype))
    integer :: i,j
    real(kind=real64) :: aux(10)
#   if CHECK_XIS==1
    call fbem_check_xi1xi2_warning(etype,xi)
#   endif
#   if CHECK_XIS==2
    call fbem_check_xi1xi2_error(etype,xi)
#   endif
#   define dphidxi1 dphidxi1_vector
#   define dphidxi2 dphidxi2_vector
#   define delta 0.0d0
#   include <dphidxik_2d.rc>
#   undef dphidxi1
#   undef dphidxi2
#   undef delta
    ! dxdxi1(1) = dx_1/dxi_1, dxdxi1(2) = dx_2/dxi_1
    ! dxdxi2(1) = dx_1/dxi_2, dxdxi2(2) = dx_2/dxi_2
    do j=1,2
      dxdxi1(j)=0.0d0
      dxdxi2(j)=0.0d0
    end do
    do i=1,fbem_n_nodes(etype)
      do j=1,2
        dxdxi1(j)=dxdxi1(j)+dphidxi1_vector(i)*x_nodes(j,i)
        dxdxi2(j)=dxdxi2(j)+dphidxi2_vector(i)*x_nodes(j,i)
      end do
    end do
    ! Jacobian
    fbem_jacobian2d_2d=dxdxi1(1)*dxdxi2(2)-dxdxi1(2)*dxdxi2(1)
  end function fbem_jacobian2d_2d

  !! Calculation of 3D jacobian of a 1D element (module of the tangent vector)
  function fbem_jacobian3d_1d(etype,x_nodes,xi)
    implicit none
    !! Element geometrical interpolation type: <tt>etype={fbem_line2,fbem_line3}</tt>.
    integer :: etype
    !! 3D position of element nodes
    real(kind=real64) :: x_nodes(3,fbem_n_nodes(etype))
    !! Reference space coordinate
    real(kind=real64) :: xi
    real(kind=real64) :: fbem_jacobian3d_1d
    real(kind=real64) :: dphidxi_vector(fbem_n_nodes(etype))
    real(kind=real64) :: tangent(3)
    real(kind=real64) :: aux(10)
    integer :: i,j
#   if CHECK_XIS==1
    call fbem_check_xi_warning(xi)
#   endif
#   if CHECK_XIS==2
    call fbem_check_xi_error(xi)
#   endif
#   define delta 0.0d0
#   define dphidxi dphidxi_vector
#   include <dphidxi_1d.rc>
#   undef dphidxi
#   undef delta
    ! Components of tangent vector calculation
    do j=1,3
      tangent(j)=0.0d0
    end do
    do i=1,fbem_n_nodes(etype)
      do j=1,3
        tangent(j)=tangent(j)+dphidxi_vector(i)*x_nodes(j,i)
      end do
    end do
    ! Jacobian as the module of the tangent vector
    fbem_jacobian3d_1d=dsqrt(tangent(1)**2+tangent(2)**2+tangent(3)**2)
  end function fbem_jacobian3d_1d

  !! Calculation of 3D jacobian of a 2D (module of the tangent vector) or 3D element (usual jacobian). Only implemented for
  !! 2D elements.
  function fbem_jacobian3d_nd(etype,x_nodes,xi)
    implicit none
    !! Element geometrical interpolation type: <tt>etype={fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    integer :: etype
    !! 3D position of element nodes
    real(kind=real64) :: x_nodes(3,fbem_n_nodes(etype))
    !! Reference space coordinates
    real(kind=real64) :: xi(fbem_n_dimension(etype))
    real(kind=real64) :: fbem_jacobian3d_nd
    real(kind=real64) :: tangent_xi1(3)
    real(kind=real64) :: tangent_xi2(3)
    real(kind=real64) :: dphidxi1_vector(fbem_n_nodes(etype)), dphidxi2_vector(fbem_n_nodes(etype))
    real(kind=real64) :: normal(3)
    integer :: i,j
    real(kind=real64) :: aux(10)
    if (fbem_n_dimension(etype).eq.2) then
#     if CHECK_XIS==1
      call fbem_check_xi1xi2_warning(etype,xi)
#     endif
#     if CHECK_XIS==2
      call fbem_check_xi1xi2_error(etype,xi)
#     endif
#     define dphidxi1 dphidxi1_vector
#     define dphidxi2 dphidxi2_vector
#     define delta 0.0d0
#     include <dphidxik_2d.rc>
#     undef dphidxi1
#     undef dphidxi2
#     undef delta
      ! Components of tangent vectors calculation
      do j=1,3
        tangent_xi1(j)=0.0d0
        tangent_xi2(j)=0.0d0
      end do
      do i=1,fbem_n_nodes(etype)
        do j=1,3
          tangent_xi1(j)=tangent_xi1(j)+dphidxi1_vector(i)*x_nodes(j,i)
          tangent_xi2(j)=tangent_xi2(j)+dphidxi2_vector(i)*x_nodes(j,i)
        end do
      end do
      ! Normal vector as T1 x T2
      normal(1)=tangent_xi1(2)*tangent_xi2(3)-tangent_xi1(3)*tangent_xi2(2)
      normal(2)=tangent_xi1(3)*tangent_xi2(1)-tangent_xi1(1)*tangent_xi2(3)
      normal(3)=tangent_xi1(1)*tangent_xi2(2)-tangent_xi1(2)*tangent_xi2(1)
      ! Jacobian as the module of the normal vector
      fbem_jacobian3d_nd=dsqrt(normal(1)**2+normal(2)**2+normal(3)**2)
    end if
  end function fbem_jacobian3d_nd

  !
  ! 2D line:
  ! M=(v1_1 n_1')
  !   (v1_2 n_2')
  !
  ! 3D surface:
  ! M=(v1_1 v2_1 n_1)
  !   (v1_2 v2_2 n_2)
  !   (v1_3 v2_3 n_3)
  !
  function fbem_local_tangential_cartesian_axes(rn,gtype,x,reverse,xi)
    implicit none
    ! I/O variables
    integer           :: rn
    integer           :: gtype
    real(kind=real64) :: x(rn,fbem_n_nodes(gtype))
    real(kind=real64) :: xi(fbem_n_dimension(gtype))
    logical           :: reverse
    real(kind=real64) :: fbem_local_tangential_cartesian_axes(rn,rn)
    ! Local variables
    integer           :: kphi
    real(kind=real64) :: xi1d
    real(kind=real64) :: aux(10)
    real(kind=real64) :: dphidxi(fbem_n_nodes(gtype))
    real(kind=real64) :: dphidxi1(fbem_n_nodes(gtype))
    real(kind=real64) :: dphidxi2(fbem_n_nodes(gtype))
    real(kind=real64) :: dphidxi3(fbem_n_nodes(gtype))
    real(kind=real64) :: J(rn,rn), invJ(rn,rn), detJ
    real(kind=real64) :: T(rn), T1(3), T2(3), N(3)
    real(kind=real64) :: jg, jxi1, jxi2
    real(kind=real64) :: v1(3), v2(3), Mvt(2,2), detMvt, dxidv(2,2)
    select case (fbem_n_dimension(gtype))
      ! Line gradient
      case (1)
        select case (rn)
          case (2)
            xi1d=xi(1)
            ! Geometrical first derivatives at xi
#           define etype gtype
#           define delta 0.0d0
#           define xi xi1d
#           include <dphidxi_1d.rc>
#           undef etype
#           undef delta
#           undef xi
            ! Calculate x_i, T1 and T2
            T=0.d0
            do kphi=1,fbem_n_nodes(gtype)
              T=T+dphidxi(kphi)*x(:,kphi)
            end do
            ! Geometric jacobian
            jg=sqrt(dot_product(T,T))
            ! Unit tangent and normal
            t=T/jg
            n(1)=t(2)
            n(2)=-t(1)
            ! Local cartesian axes: v1==-t, v2==n
            fbem_local_tangential_cartesian_axes(1,1)=-t(1)
            fbem_local_tangential_cartesian_axes(2,1)=-t(2)
            fbem_local_tangential_cartesian_axes(1,2)=n(1)
            fbem_local_tangential_cartesian_axes(2,2)=n(2)
            if (reverse) fbem_local_tangential_cartesian_axes=-fbem_local_tangential_cartesian_axes
          case default
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'element not valid for tangential cartesian axes calculation')
        end select
      case (2)
        select case (rn)
          ! Surface tangential gradient
          case (3)
            ! Geometrical first derivatives at xi
#           define etype gtype
#           define delta 0.0d0
#           include <dphidxik_2d.rc>
#           undef etype
#           undef delta
            ! Calculate x_i, T1 and T2
            T1=0.d0
            T2=0.d0
            do kphi=1,fbem_n_nodes(gtype)
              T1=T1+dphidxi1(kphi)*x(:,kphi)
              T2=T2+dphidxi2(kphi)*x(:,kphi)
            end do
            ! Normal vector as T1 x T2 at xi_i
            N(1)=T1(2)*T2(3)-T1(3)*T2(2)
            N(2)=T1(3)*T2(1)-T1(1)*T2(3)
            N(3)=T1(1)*T2(2)-T1(2)*T2(1)
            ! Geometric jacobian
            jg=sqrt(dot_product(N,N))
            ! Unit normal
            n=N/jg
            ! In jt1 and jt2 it is saved the tangent vectors norm
            jxi1=sqrt(dot_product(T1,T1))
            t1=T1/jxi1
            jxi2=sqrt(dot_product(T2,T2))
            t2=T2/jxi2
            ! The orthogonal system at xi
            ! v1 = t1
            v1=t1
            ! v2 = n x t1
            v2(1)=n(2)*t1(3)-n(3)*t1(2)
            v2(2)=n(3)*t1(1)-n(1)*t1(3)
            v2(3)=n(1)*t1(2)-n(2)*t1(1)
            ! Local cartesian axes: v1, v2, n
            fbem_local_tangential_cartesian_axes(1,1)=v1(1)
            fbem_local_tangential_cartesian_axes(2,1)=v1(2)
            fbem_local_tangential_cartesian_axes(3,1)=v1(3)
            fbem_local_tangential_cartesian_axes(1,2)=v2(1)
            fbem_local_tangential_cartesian_axes(2,2)=v2(2)
            fbem_local_tangential_cartesian_axes(3,2)=v2(3)
            fbem_local_tangential_cartesian_axes(1,3)=n(1)
            fbem_local_tangential_cartesian_axes(2,3)=n(2)
            fbem_local_tangential_cartesian_axes(3,3)=n(3)
            if (reverse) then
              fbem_local_tangential_cartesian_axes(:,1)=-fbem_local_tangential_cartesian_axes(:,1)
              fbem_local_tangential_cartesian_axes(:,3)=-fbem_local_tangential_cartesian_axes(:,3)
            end if
        end select
      ! Volume gradient
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'element not valid for tangential cartesian axes calculation')
    end select
  end function fbem_local_tangential_cartesian_axes







  ! =====================================
  ! Dataset at element integration points
  ! =====================================

  !! Calculate functional shape functions, position vectors, normals, geometric jacobian*weight on integration points for
  !! given rules.
  subroutine fbem_dataset_at_integration_points(rn,etype,deltaf,x_nodes,nrules,glrule,ngp,fphi,x,n,jacw)
    implicit none
    ! I/O
    integer                        :: rn                              !! Dimensional space
    integer                        :: etype                           !! Element functional interpolation type: <tt>etype={fbem_line2,fbem_line3,fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    real(kind=real64)              :: deltaf                          !! Displacement for discontinuous elements
    real(kind=real64)              :: x_nodes(rn,fbem_n_nodes(etype)) !! Position vectors of geometrical nodes
    integer                        :: nrules                          !! Number of rules where evaluate
    integer                        :: glrule(nrules)                  !! Number of Gauss points in each direction of Gauss-Legendre quadrature. For triangular elements, if rule(i)<=15 then Wandzura symmetrical rule is used.
    integer                        :: ngp(nrules)                     !! Number of integration points for each rule
    real(kind=real64), allocatable :: fphi(:,:,:)                     !! Functional shape functions on integration points
    real(kind=real64), allocatable :: x(:,:,:)                        !! Position vectors on integration points
    real(kind=real64), allocatable :: n(:,:,:)                        !! Unit normal vector on integration points
    real(kind=real64), allocatable :: jacw(:,:)                       !! Geometric jacobian multiplied by weight on integration points
    ! Local
    integer           :: i, kphi, k, kxi, k1, k2, kt,nnodes, order
    real(kind=real64) :: aux(10)
    real(kind=real64) :: gphi_vector(fbem_n_nodes(etype))
    real(kind=real64) :: dphidxi_vector(fbem_n_nodes(etype))
    real(kind=real64) :: dphidxi1_vector(fbem_n_nodes(etype))
    real(kind=real64) :: dphidxi2_vector(fbem_n_nodes(etype))
    real(kind=real64) :: fphi_vector(fbem_n_nodes(etype))
    real(kind=real64) :: xi_p(2), wgp(2), xp(rn), T(rn), T1(rn), T2(rn), Np(rn), jac
    !
    ! Number of nodes
    !
    nnodes=fbem_n_nodes(etype)
    !
    ! Calculate the number of integration points for each rule
    !
    select case (fbem_n_dimension(etype))
      ! Line elements
      case (1)
        do i=1,nrules
          ngp(i)=glrule(i)
        end do
      ! Surface elements
      case (2)
        select case(fbem_n_edges(etype))
          ! Triangular elements
          case (3)
            do i=1,nrules
              ! If ngp>15 then Gauss-Legendre*Gauss-Jacobi
              if (glrule(i).gt.15) then
                ngp(i)=glrule(i)**2
              ! If ngp<=15 Wandzura rules are used
              else
                ! Conversion between ngp and quadrature order
                order=2*glrule(i)-1
                ! Number of integration points
                ngp(i)=wantri_n(order)
              end if
            end do
          ! Quadrilateral elements
          case (4)
            do i=1,nrules
              ngp(i)=glrule(i)**2
            end do
        case default
          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'only triangular and quadrilateral elements')
        end select
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'only line and surface elements')
    end select
    !
    ! Deallocate and allocate the dataset
    !
    if (allocated(fphi).eqv.(.true.)) deallocate(fphi)
    if (allocated(x).eqv.(.true.)) deallocate(x)
    if (allocated(n).eqv.(.true.)) deallocate(n)
    if (allocated(jacw).eqv.(.true.)) deallocate(jacw)
    allocate (fphi(nnodes,maxval(ngp),nrules))
    allocate (x(rn,maxval(ngp),nrules))
    allocate (n(rn,maxval(ngp),nrules))
    allocate (jacw(maxval(ngp),nrules))
    !
    ! Calculate the dataset
    !
    select case (fbem_n_dimension(etype))
      !
      ! Line elements
      !
      case (1)
       ! Loop through rules
        do i=1,nrules
          ! Loop through xi coordinate
          do kxi=1,gl11_n(glrule(i))
            ! Index of gaussian point
            kt=kxi
            ! xi coordinate and weight
            xi_p(1)=gl11_xi(kxi,glrule(i))
            wgp(1)=gl11_w(kxi,glrule(i))
            ! Geometrical shape functions and first derivatives at xi
#           define xi xi_p(1)
#           define delta 0.0d0
#           define phi gphi_vector
#           define dphidxi dphidxi_vector
#           include <phi_and_dphidxi_1d.rc>
#           undef xi
#           undef phi
#           undef dphidxi
#           undef delta
            ! Functional shape functions at xi
#           define xi xi_p(1)
#           define delta deltaf
#           define phi fphi_vector
#           include <phi_1d.rc>
#           undef xi
#           undef delta
#           undef phi
            ! Components calculation of x and T at xi
            do k=1,rn
              xp(k)=0.0d0
              T(k)=0.0d0
            end do
            do kphi=1,nnodes
              do k=1,rn
                xp(k)=xp(k)+gphi_vector(kphi)*x_nodes(k,kphi)
                T(k)=T(k)+dphidxi_vector(kphi)*x_nodes(k,kphi)
              end do
            end do
            ! Geometric jacobian
            select case (rn)
              case (2)
                jac=dsqrt(T(1)**2+T(2)**2)
              case (3)
                jac=dsqrt(T(1)**2+T(2)**2+T(3)**2)
            end select
            ! Position vector
            do k=1,rn
              x(k,kt,i)=xp(k)
            end do
            ! Unit normal vector (note that it is not valid for 3D, it does not make sense)
            select case (rn)
              ! 2D: the unit normal is saved in n
              case (2)
                Np(1)=T(2)
                Np(2)=-T(1)
                n(1,kt,i)=Np(1)/jac
                n(2,kt,i)=Np(2)/jac
              ! 3D: the unit normal is undefined, so the unit tangent is saved in n
              case (3)
                n(1,kt,i)=T(1)/jac
                n(2,kt,i)=T(2)/jac
                n(3,kt,i)=T(3)/jac
            end select
            ! Jacobian multiply by weights
            jacw(kt,i)=jac*wgp(1)
            ! Functional shape functions
            do kphi=1,nnodes
              fphi(kphi,kt,i)=fphi_vector(kphi)
            end do
          end do ! Loop through xi coordinate
        end do ! Loop through rules
      !
      ! Surface elements
      !
      case (2)
        select case (fbem_n_edges(etype))
          !
          ! Quadrilateral elements
          !
          case (4)
            ! Loop through rules
            do i=1,nrules
              ! Loop through xi_1 direction
              do k1=1,gl11_n(glrule(i))
                ! xi_1 coordinate and weight
                xi_p(1)=gl11_xi(k1,glrule(i))
                wgp(1)=gl11_w(k1,glrule(i))
                ! Loop through xi_2 direction
                do k2=1,gl11_n(glrule(i))
                  ! Total index of gaussian point
                  kt=k2+(k1-1)*gl11_n(glrule(i))
                  ! xi_2 coordinate and weight
                  xi_p(2)=gl11_xi(k2,glrule(i))
                  wgp(2)=gl11_w(k2,glrule(i))
                  ! Geometrical shape functions and first derivatives at xi
#                 define xi xi_p
#                 define delta 0.0d0
#                 define phi gphi_vector
#                 define dphidxi1 dphidxi1_vector
#                 define dphidxi2 dphidxi2_vector
#                 include <phi_and_dphidxik_2d.rc>
#                 undef xi
#                 undef phi
#                 undef dphidxi1
#                 undef dphidxi2
#                 undef delta
                  ! Functional shape functions at xi
#                 define xi xi_p
#                 define delta deltaf
#                 define phi fphi_vector
#                 include <phi_2d.rc>
#                 undef xi
#                 undef delta
#                 undef phi
                  ! Components calculation of x, T1 and T2 at xi
                  do k=1,rn
                    xp(k)=0.0d0
                    T1(k)=0.0d0
                    T2(k)=0.0d0
                  end do
                  do kphi=1,nnodes
                    do k=1,rn
                      xp(k)=xp(k)+gphi_vector(kphi)*x_nodes(k,kphi)
                      T1(k)=T1(k)+dphidxi1_vector(kphi)*x_nodes(k,kphi)
                      T2(k)=T2(k)+dphidxi2_vector(kphi)*x_nodes(k,kphi)
                    end do
                  end do
                  ! Jacobian and normal
                  select case(rn)
                    ! 2D: the jacobian is the classic jacobian, but the normal is in 3D space, so it is undefined
                    case (2)
                      ! Normal (undefined, set to 0.)
                      Np(1)=0.0d0
                      Np(2)=0.0d0
                      ! Jacobian
                      jac=T1(1)*T2(2)-T1(2)*T2(1)
                    ! 3D: the jacobian is |T1 x T2|, which is a non-unit normal
                    case (3)
                      ! Normal vector as T1 x T2 at xi
                      Np(1)=T1(2)*T2(3)-T1(3)*T2(2)
                      Np(2)=T1(3)*T2(1)-T1(1)*T2(3)
                      Np(3)=T1(1)*T2(2)-T1(2)*T2(1)
                      ! Jacobian
                      jac=dsqrt(Np(1)**2+Np(2)**2+Np(3)**2)
                  end select
                  ! Position vector
                  do k=1,rn
                    x(k,kt,i)=xp(k)
                  end do
                  ! Unit normal vector
                  do k=1,rn
                    n(k,kt,i)=Np(k)/jac
                  end do
                  ! Jacobian multiply by weights
                  jacw(kt,i)=jac*wgp(1)*wgp(2)
                  ! Functional shape functions
                  do kphi=1,nnodes
                    fphi(kphi,kt,i)=fphi_vector(kphi)
                  end do
                end do ! Loop through xi_2 direction
              end do ! Loop through xi_1 direction
            end do ! Loop through rules
          !
          ! Triangular elements
          !
          case (3)
            ! Loop through rules
            do i=1,nrules
              ! If ngp>15 then Gauss-Legendre*Gauss-Jacobi
              if (glrule(i).gt.15) then
                ! Loop through rule points
                do k1=1,gl01_n(glrule(i))
                  do k2=1,gj01_n(glrule(i))
                    ! Total index of gaussian point
                    kt=k2+(k1-1)*gl01_n(glrule(i))
                    ! xi_1 coordinate and xi_2 coordinate
                    xi_p(1)=(1.0d0-gj01_xi(k2,glrule(i)))*gl01_xi(k1,glrule(i))
                    xi_p(2)=gj01_xi(k2,glrule(i))
                    ! Geometrical shape functions and first derivatives at xi
#                   define xi xi_p
#                   define delta 0.0d0
#                   define phi gphi_vector
#                   define dphidxi1 dphidxi1_vector
#                   define dphidxi2 dphidxi2_vector
#                   include <phi_and_dphidxik_2d.rc>
#                   undef xi
#                   undef delta
#                   undef phi
#                   undef dphidxi1
#                   undef dphidxi2
                    ! Functional shape functions at xi
#                   define xi xi_p
#                   define delta deltaf
#                   define phi fphi_vector
#                   include <phi_2d.rc>
#                   undef xi
#                   undef delta
#                   undef phi
                    ! Components calculation of x, T1 and T2 at xi
                    do k=1,rn
                      xp(k)=0.0d0
                      T1(k)=0.0d0
                      T2(k)=0.0d0
                    end do
                    do kphi=1,nnodes
                      do k=1,rn
                        xp(k)=xp(k)+gphi_vector(kphi)*x_nodes(k,kphi)
                        T1(k)=T1(k)+dphidxi1_vector(kphi)*x_nodes(k,kphi)
                        T2(k)=T2(k)+dphidxi2_vector(kphi)*x_nodes(k,kphi)
                      end do
                    end do
                    ! Jacobian and normal
                    select case(rn)
                      ! 2D: the jacobian is the classic jacobian, but the normal is in 3D space, so it is undefined
                      case (2)
                        ! Normal (undefined, set to 0.)
                        Np(1)=0.0d0
                        Np(2)=0.0d0
                        ! Jacobian
                        jac=T1(1)*T2(2)-T1(2)*T2(1)
                      ! 3D: the jacobian is |T1 x T2|, which is a non-unit normal
                      case (3)
                        ! Normal vector as T1 x T2 at xi
                        Np(1)=T1(2)*T2(3)-T1(3)*T2(2)
                        Np(2)=T1(3)*T2(1)-T1(1)*T2(3)
                        Np(3)=T1(1)*T2(2)-T1(2)*T2(1)
                        ! Jacobian
                        jac=dsqrt(Np(1)**2+Np(2)**2+Np(3)**2)
                    end select
                    ! Position vector
                    do k=1,rn
                      x(k,kt,i)=xp(k)
                    end do
                    ! Unit normal vector
                    do k=1,rn
                      n(k,kt,i)=Np(k)/jac
                    end do
                    ! Jacobian * weights
                    jacw(kt,i)=jac*gl01_w(k1,glrule(i))*gj01_w(k2,glrule(i))
                    ! Functional shape functions
                    do kphi=1,nnodes
                      fphi(kphi,kt,i)=fphi_vector(kphi)
                    end do
                  end do
                end do
              ! If ngp<=15 Wandzura rules are used
              else
                ! Conversion between ngp and quadrature order
                order=2*glrule(i)-1
                ! Loop through rule points
                do kt=1,wantri_n(order)
                  ! xi_1 coordinate
                  xi_p(1)=wantri_xi1(kt,order)
                  ! xi_2 coordinate
                  xi_p(2)=wantri_xi2(kt,order)
                  ! Weight
                  wgp(1)=wantri_w(kt,order)
                    ! Geometrical shape functions and first derivatives at xi
#                   define xi xi_p
#                   define delta 0.0d0
#                   define phi gphi_vector
#                   define dphidxi1 dphidxi1_vector
#                   define dphidxi2 dphidxi2_vector
#                   include <phi_and_dphidxik_2d.rc>
#                   undef xi
#                   undef delta
#                   undef phi
#                   undef dphidxi1
#                   undef dphidxi2
                    ! Functional shape functions at xi
#                   define xi xi_p
#                   define delta deltaf
#                   define phi fphi_vector
#                   include <phi_2d.rc>
#                   undef xi
#                   undef delta
#                   undef phi
                  ! Components calculation of x, T1 and T2 at xi
                  do k=1,rn
                    xp(k)=0.0d0
                    T1(k)=0.0d0
                    T2(k)=0.0d0
                  end do
                  do kphi=1,nnodes
                    do k=1,rn
                      xp(k)=xp(k)+gphi_vector(kphi)*x_nodes(k,kphi)
                      T1(k)=T1(k)+dphidxi1_vector(kphi)*x_nodes(k,kphi)
                      T2(k)=T2(k)+dphidxi2_vector(kphi)*x_nodes(k,kphi)
                    end do
                  end do
                  ! Jacobian and normal
                  select case(rn)
                    ! 2D: the jacobian is the classic jacobian, but the normal is in 3D space, so it is undefined
                    case (2)
                      ! Normal (undefined, set to 0.)
                      Np(1)=0.0d0
                      Np(2)=0.0d0
                      ! Jacobian
                      jac=T1(1)*T2(2)-T1(2)*T2(1)
                    ! 3D: the jacobian is |T1 x T2|, which is a non-unit normal
                    case (3)
                      ! Normal vector as T1 x T2 at xi
                      Np(1)=T1(2)*T2(3)-T1(3)*T2(2)
                      Np(2)=T1(3)*T2(1)-T1(1)*T2(3)
                      Np(3)=T1(1)*T2(2)-T1(2)*T2(1)
                      ! Jacobian
                      jac=dsqrt(Np(1)**2+Np(2)**2+Np(3)**2)
                  end select
                  ! Position vector
                  do k=1,rn
                    x(k,kt,i)=xp(k)
                  end do
                  ! Unit normal vector
                  do k=1,rn
                    n(k,kt,i)=Np(k)/jac
                  end do
                  ! Jacobian multiply by weight
                  jacw(kt,i)=jac*wgp(1)
                  ! Functional shape functions
                  do kphi=1,nnodes
                    fphi(kphi,kt,i)=fphi_vector(kphi)
                  end do
                end do
              end if
            end do ! Loop through rules
        end select
      case default
          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={line2,line3,tri3,tri6,quad4,quad8,quad9}')
    end select
  end subroutine fbem_dataset_at_integration_points

  !! Calculate functional shape functions, position vectors, normals, geometric jacobian*weight on integration points for
  !! given rules. This subroutine lets to introduce different shape functions for primary variables, secondary variables and
  !! geometry.
  subroutine fbem_dataset_at_integration_points_dif(rn,typeg,typef1,typef2,deltaf,x_nodes,nrules,glrule,ngp,f1phi,f2phi,x,n,jacw)
    implicit none
    ! I/O
    integer                        :: rn                              !! Dimensional space
    integer                        :: typeg                           !! Geometrical interpolation
    integer                        :: typef1                          !! Functional interpolation (primary variables)
    integer                        :: typef2                          !! Functional interpolation (secondary variables)
    real(kind=real64)              :: deltaf                          !! Displacement of discontinuous elements (0 for continuous)
    real(kind=real64)              :: x_nodes(rn,fbem_n_nodes(typeg)) !! Position vectors of geometrical nodes
    integer                        :: nrules                          !! Number of rules where evaluate
    integer                        :: glrule(nrules)                  !! Number of Gauss points in each direction of Gauss-Legendre quadrature. For triangular elements, if rule(i)<=15 then Wandzura symmetrical rule is used.
    integer                        :: ngp(nrules)                     !! Number of integration points for each rule
    real(kind=real64), allocatable :: f1phi(:,:,:)                    !! Functional shape functions (primary variables) on integration points
    real(kind=real64), allocatable :: f2phi(:,:,:)                    !! Functional shape functions (secondary variables) on integration points
    real(kind=real64), allocatable :: x(:,:,:)                        !! Position vectors on integration points
    real(kind=real64), allocatable :: n(:,:,:)                        !! Unit normal vector on integration points
    real(kind=real64), allocatable :: jacw(:,:)                       !! Geometric jacobian multiplied by weight on integration points
    ! Local
    integer           :: i, kphi, k, kxi, k1, k2, kt, order
    integer           :: nnodesg, nnodesf1, nnodesf2
    real(kind=real64) :: aux(10)
    real(kind=real64) :: gphi_vector(fbem_n_nodes(typeg))
    real(kind=real64) :: dphidxi_vector(fbem_n_nodes(typeg))
    real(kind=real64) :: dphidxi1_vector(fbem_n_nodes(typeg))
    real(kind=real64) :: dphidxi2_vector(fbem_n_nodes(typeg))
    real(kind=real64) :: f1phi_vector(fbem_n_nodes(typef1))
    real(kind=real64) :: f2phi_vector(fbem_n_nodes(typef2))
    real(kind=real64) :: xi_p(2), wgp(2), xp(rn), T(rn), T1(rn), T2(rn), Np(rn), jac
    !
    ! Number of nodes
    !
    nnodesg=fbem_n_nodes(typeg)
    nnodesf1=fbem_n_nodes(typef1)
    nnodesf2=fbem_n_nodes(typef2)
    !
    ! Calculate the number of integration points for each rule
    !
    select case (typeg)
      ! Line elements
      case (fbem_line2,fbem_line3,fbem_line4)
        do i=1,nrules
          ngp(i)=glrule(i)
        end do
      ! Quadrilateral elements
      case (fbem_quad4,fbem_quad8,fbem_quad9)
        do i=1,nrules
          ngp(i)=glrule(i)**2
        end do
      ! Triangular elements
      case (fbem_tri3,fbem_tri6)
        do i=1,nrules
          ! If ngp>15 then Gauss-Legendre*Gauss-Jacobi
          if (glrule(i).gt.15) then
            ngp(i)=glrule(i)**2
          ! If ngp<=15 Wandzura rules are used
          else
            ! Conversion between ngp and quadrature order
            order=2*glrule(i)-1
            ! Number of integration points
            ngp(i)=wantri_n(order)
          end if
        end do
    end select
    !
    ! Deallocate and allocate the dataset
    !
    if (allocated(f1phi).eqv.(.true.)) deallocate(f1phi)
    if (allocated(f2phi).eqv.(.true.)) deallocate(f2phi)
    if (allocated(x    ).eqv.(.true.)) deallocate(x    )
    if (allocated(n    ).eqv.(.true.)) deallocate(n    )
    if (allocated(jacw ).eqv.(.true.)) deallocate(jacw )
    allocate (f1phi(nnodesf1,maxval(ngp),nrules))
    allocate (f2phi(nnodesf2,maxval(ngp),nrules))
    allocate (x(rn,maxval(ngp),nrules))
    allocate (n(rn,maxval(ngp),nrules))
    allocate (jacw(maxval(ngp),nrules))
    !
    ! Calculate the dataset
    !
    select case (typeg)
      !
      ! Line elements
      !
      case (fbem_line2,fbem_line3,fbem_line4)
       ! Loop through rules
        do i=1,nrules
          ! Loop through xi coordinate
          do kxi=1,gl11_n(glrule(i))
            ! Index of gaussian point
            kt=kxi
            ! xi coordinate and weight
            xi_p(1)=gl11_xi(kxi,glrule(i))
            wgp(1)=gl11_w(kxi,glrule(i))
            ! Geometrical shape functions and first derivatives at xi
#           define etype typeg
#           define xi xi_p(1)
#           define delta 0.0d0
#           define phi gphi_vector
#           define dphidxi dphidxi_vector
#           include <phi_and_dphidxi_1d.rc>
#           undef etype
#           undef xi
#           undef phi
#           undef dphidxi
#           undef delta
            ! Functional shape functions (primary variables) at xi
#           define etype typef1
#           define xi xi_p(1)
#           define delta deltaf
#           define phi f1phi_vector
#           include <phi_1d.rc>
#           undef etype
#           undef xi
#           undef delta
#           undef phi
            ! Functional shape functions (secondary variables) at xi
#           define etype typef2
#           define xi xi_p(1)
#           define delta deltaf
#           define phi f2phi_vector
#           include <phi_1d.rc>
#           undef etype
#           undef xi
#           undef delta
#           undef phi
            ! Components calculation of x and T at xi
            do k=1,rn
              xp(k)=0.0d0
              T(k)=0.0d0
            end do
            do kphi=1,nnodesg
              do k=1,rn
                xp(k)=xp(k)+gphi_vector(kphi)*x_nodes(k,kphi)
                T(k)=T(k)+dphidxi_vector(kphi)*x_nodes(k,kphi)
              end do
            end do
            ! Geometric jacobian
            select case (rn)
              case (2)
                jac=dsqrt(T(1)**2+T(2)**2)
              case (3)
                jac=dsqrt(T(1)**2+T(2)**2+T(3)**2)
            end select
            ! Position vector
            do k=1,rn
              x(k,kt,i)=xp(k)
            end do
            ! Unit normal vector (note that it is not valid for 3D, it does not make sense)
            select case (rn)
              ! 2D: the unit normal is saved in n
              case (2)
                Np(1)=T(2)
                Np(2)=-T(1)
                n(1,kt,i)=Np(1)/jac
                n(2,kt,i)=Np(2)/jac
              ! 3D: the unit normal is undefined, so the unit tangent is saved in n
              case (3)
                n(1,kt,i)=T(1)/jac
                n(2,kt,i)=T(2)/jac
                n(3,kt,i)=T(3)/jac
            end select
            ! Jacobian multiply by weights
            jacw(kt,i)=jac*wgp(1)
            ! Functional shape functions
            do kphi=1,nnodesf1
              f1phi(kphi,kt,i)=f1phi_vector(kphi)
            end do
            do kphi=1,nnodesf2
              f2phi(kphi,kt,i)=f2phi_vector(kphi)
            end do
          end do ! Loop through xi coordinate
        end do ! Loop through rules
      !
      ! Quadrilateral elements
      !
      case (fbem_quad4,fbem_quad8,fbem_quad9)
        ! Loop through rules
        do i=1,nrules
          ! Loop through xi_1 direction
          do k1=1,gl11_n(glrule(i))
            ! xi_1 coordinate and weight
            xi_p(1)=gl11_xi(k1,glrule(i))
            wgp(1)=gl11_w(k1,glrule(i))
            ! Loop through xi_2 direction
            do k2=1,gl11_n(glrule(i))
              ! Total index of gaussian point
              kt=k2+(k1-1)*gl11_n(glrule(i))
              ! xi_2 coordinate and weight
              xi_p(2)=gl11_xi(k2,glrule(i))
              wgp(2)=gl11_w(k2,glrule(i))
              ! Geometrical shape functions and first derivatives at xi
#             define etype typeg
#             define xi xi_p
#             define delta 0.0d0
#             define phi gphi_vector
#             define dphidxi1 dphidxi1_vector
#             define dphidxi2 dphidxi2_vector
#             include <phi_and_dphidxik_2d.rc>
#             undef etype
#             undef xi
#             undef phi
#             undef dphidxi1
#             undef dphidxi2
#             undef delta
              ! Functional shape functions (primary variables) at xi
#             define etype typef1
#             define xi xi_p
#             define delta deltaf
#             define phi f1phi_vector
#             include <phi_2d.rc>
#             undef etype
#             undef xi
#             undef delta
#             undef phi
              ! Functional shape functions (secondary variables) at xi
#             define etype typef2
#             define xi xi_p
#             define delta deltaf
#             define phi f2phi_vector
#             include <phi_2d.rc>
#             undef etype
#             undef xi
#             undef delta
#             undef phi
              ! Components calculation of x, T1 and T2 at xi
              do k=1,rn
                xp(k)=0.0d0
                T1(k)=0.0d0
                T2(k)=0.0d0
              end do
              do kphi=1,nnodesg
                do k=1,rn
                  xp(k)=xp(k)+gphi_vector(kphi)*x_nodes(k,kphi)
                  T1(k)=T1(k)+dphidxi1_vector(kphi)*x_nodes(k,kphi)
                  T2(k)=T2(k)+dphidxi2_vector(kphi)*x_nodes(k,kphi)
                end do
              end do
              ! Jacobian and normal
              select case(rn)
                ! 2D: the jacobian is the classic jacobian, but the normal is in 3D space, so it is undefined
                case (2)
                  ! Normal (undefined, set to 0.)
                  Np(1)=0.0d0
                  Np(2)=0.0d0
                  ! Jacobian
                  jac=T1(1)*T2(2)-T1(2)*T2(1)
                ! 3D: the jacobian is |T1 x T2|, which is a non-unit normal
                case (3)
                  ! Normal vector as T1 x T2 at xi
                  Np(1)=T1(2)*T2(3)-T1(3)*T2(2)
                  Np(2)=T1(3)*T2(1)-T1(1)*T2(3)
                  Np(3)=T1(1)*T2(2)-T1(2)*T2(1)
                  ! Jacobian
                  jac=dsqrt(Np(1)**2+Np(2)**2+Np(3)**2)
              end select
              ! Position vector
              do k=1,rn
                x(k,kt,i)=xp(k)
              end do
              ! Unit normal vector
              do k=1,rn
                n(k,kt,i)=Np(k)/jac
              end do
              ! Jacobian multiply by weights
              jacw(kt,i)=jac*wgp(1)*wgp(2)
              ! Functional shape functions
              do kphi=1,nnodesf1
                f1phi(kphi,kt,i)=f1phi_vector(kphi)
              end do
              do kphi=1,nnodesf2
                f2phi(kphi,kt,i)=f2phi_vector(kphi)
              end do
            end do ! Loop through xi_2 direction
          end do ! Loop through xi_1 direction
        end do ! Loop through rules
      !
      ! Triangular elements
      !
      case (fbem_tri3,fbem_tri6)
        ! Loop through rules
        do i=1,nrules
          ! If ngp>15 then Gauss-Legendre*Gauss-Jacobi
          if (glrule(i).gt.15) then
            ! Loop through rule points
            do k1=1,gl01_n(glrule(i))
              do k2=1,gj01_n(glrule(i))
                ! Total index of gaussian point
                kt=k2+(k1-1)*gl01_n(glrule(i))
                ! xi_1 coordinate and xi_2 coordinate
                xi_p(1)=(1.0d0-gj01_xi(k2,glrule(i)))*gl01_xi(k1,glrule(i))
                xi_p(2)=gj01_xi(k2,glrule(i))
                ! Geometrical shape functions and first derivatives at xi
#               define etype typeg
#               define xi xi_p
#               define delta 0.0d0
#               define phi gphi_vector
#               define dphidxi1 dphidxi1_vector
#               define dphidxi2 dphidxi2_vector
#               include <phi_and_dphidxik_2d.rc>
#               undef etype
#               undef xi
#               undef delta
#               undef phi
#               undef dphidxi1
#               undef dphidxi2
                ! Functional shape functions (primary variables) at xi
#               define etype typef1
#               define xi xi_p
#               define delta deltaf
#               define phi f1phi_vector
#               include <phi_2d.rc>
#               undef etype
#               undef xi
#               undef delta
#               undef phi
                ! Functional shape functions (secondary variables) at xi
#               define etype typef2
#               define xi xi_p
#               define delta deltaf
#               define phi f2phi_vector
#               include <phi_2d.rc>
#               undef etype
#               undef xi
#               undef delta
#               undef phi
                ! Components calculation of x, T1 and T2 at xi
                do k=1,rn
                  xp(k)=0.0d0
                  T1(k)=0.0d0
                  T2(k)=0.0d0
                end do
                do kphi=1,nnodesg
                  do k=1,rn
                    xp(k)=xp(k)+gphi_vector(kphi)*x_nodes(k,kphi)
                    T1(k)=T1(k)+dphidxi1_vector(kphi)*x_nodes(k,kphi)
                    T2(k)=T2(k)+dphidxi2_vector(kphi)*x_nodes(k,kphi)
                  end do
                end do
                ! Jacobian and normal
                select case(rn)
                  ! 2D: the jacobian is the classic jacobian, but the normal is in 3D space, so it is undefined
                  case (2)
                    ! Normal (undefined, set to 0.)
                    Np(1)=0.0d0
                    Np(2)=0.0d0
                    ! Jacobian
                    jac=T1(1)*T2(2)-T1(2)*T2(1)
                  ! 3D: the jacobian is |T1 x T2|, which is a non-unit normal
                  case (3)
                    ! Normal vector as T1 x T2 at xi
                    Np(1)=T1(2)*T2(3)-T1(3)*T2(2)
                    Np(2)=T1(3)*T2(1)-T1(1)*T2(3)
                    Np(3)=T1(1)*T2(2)-T1(2)*T2(1)
                    ! Jacobian
                    jac=dsqrt(Np(1)**2+Np(2)**2+Np(3)**2)
                end select
                ! Position vector
                do k=1,rn
                  x(k,kt,i)=xp(k)
                end do
                ! Unit normal vector
                do k=1,rn
                  n(k,kt,i)=Np(k)/jac
                end do
                ! Jacobian * weights
                jacw(kt,i)=jac*gl01_w(k1,glrule(i))*gj01_w(k2,glrule(i))
                ! Functional shape functions
                do kphi=1,nnodesf1
                  f1phi(kphi,kt,i)=f1phi_vector(kphi)
                end do
                do kphi=1,nnodesf2
                  f2phi(kphi,kt,i)=f2phi_vector(kphi)
                end do
              end do
            end do
          ! If ngp<=15 Wandzura rules are used
          else
            ! Conversion between ngp and quadrature order
            order=2*glrule(i)-1
            ! Loop through rule points
            do kt=1,wantri_n(order)
              ! xi_1 coordinate
              xi_p(1)=wantri_xi1(kt,order)
              ! xi_2 coordinate
              xi_p(2)=wantri_xi2(kt,order)
              ! Weight
              wgp(1)=wantri_w(kt,order)
              ! Geometrical shape functions and first derivatives at xi
#             define etype typeg
#             define xi xi_p
#             define delta 0.0d0
#             define phi gphi_vector
#             define dphidxi1 dphidxi1_vector
#             define dphidxi2 dphidxi2_vector
#             include <phi_and_dphidxik_2d.rc>
#             undef etype
#             undef xi
#             undef delta
#             undef phi
#             undef dphidxi1
#             undef dphidxi2
              ! Functional shape functions (primary variables) at xi
#             define etype typef1
#             define xi xi_p
#             define delta deltaf
#             define phi f1phi_vector
#             include <phi_2d.rc>
#             undef etype
#             undef xi
#             undef delta
#             undef phi
              ! Functional shape functions (secondary variables) at xi
#             define etype typef2
#             define xi xi_p
#             define delta deltaf
#             define phi f2phi_vector
#             include <phi_2d.rc>
#             undef etype
#             undef xi
#             undef delta
#             undef phi
              ! Components calculation of x, T1 and T2 at xi
              do k=1,rn
                xp(k)=0.0d0
                T1(k)=0.0d0
                T2(k)=0.0d0
              end do
              do kphi=1,nnodesg
                do k=1,rn
                  xp(k)=xp(k)+gphi_vector(kphi)*x_nodes(k,kphi)
                  T1(k)=T1(k)+dphidxi1_vector(kphi)*x_nodes(k,kphi)
                  T2(k)=T2(k)+dphidxi2_vector(kphi)*x_nodes(k,kphi)
                end do
              end do
              ! Jacobian and normal
              select case(rn)
                ! 2D: the jacobian is the classic jacobian, but the normal is in 3D space, so it is undefined
                case (2)
                  ! Normal (undefined, set to 0.)
                  Np(1)=0.0d0
                  Np(2)=0.0d0
                  ! Jacobian
                  jac=T1(1)*T2(2)-T1(2)*T2(1)
                ! 3D: the jacobian is |T1 x T2|, which is a non-unit normal
                case (3)
                  ! Normal vector as T1 x T2 at xi
                  Np(1)=T1(2)*T2(3)-T1(3)*T2(2)
                  Np(2)=T1(3)*T2(1)-T1(1)*T2(3)
                  Np(3)=T1(1)*T2(2)-T1(2)*T2(1)
                  ! Jacobian
                  jac=dsqrt(Np(1)**2+Np(2)**2+Np(3)**2)
              end select
              ! Position vector
              do k=1,rn
                x(k,kt,i)=xp(k)
              end do
              ! Unit normal vector
              do k=1,rn
                n(k,kt,i)=Np(k)/jac
              end do
              ! Jacobian multiply by weight
              jacw(kt,i)=jac*wgp(1)
              ! Functional shape functions
              do kphi=1,nnodesf1
                f1phi(kphi,kt,i)=f1phi_vector(kphi)
              end do
              do kphi=1,nnodesf2
                f2phi(kphi,kt,i)=f2phi_vector(kphi)
              end do
            end do
          end if
        end do ! Loop through rules
      case default
          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'typeg={line2,line3,tri3,tri6,quad4,quad8,quad9}')
    end select
  end subroutine fbem_dataset_at_integration_points_dif

  ! ================================================================================================================================
  ! ELEMENT SUBDIVISION
  ! ================================================================================================================================

  !! Obtain the coordinates of an element subdivision
  subroutine fbem_obtain_element_subdivision_coordinates(rn,gtype,x,xi_s,x_s)
    implicit none
    ! I/O
    integer           :: rn                                                     !! Spatial dimension
    integer           :: gtype                                                  !! Type of geometrical interpolation
    real(kind=real64) :: x(rn,fbem_n_nodes(gtype))                              !! Coordinates of the nodes of the parent element
    real(kind=real64) :: xi_s(fbem_n_dimension(gtype),fbem_n_vertices(gtype))   !! Local coordinates of the subdivision in the parent element
    real(kind=real64) :: x_s(rn,fbem_n_nodes(gtype))                            !! Coordinates of the nodes of the element subdivision
    ! Local
    integer           :: k, i
    real(kind=real64) :: xis_subdivision_nodes(fbem_n_dimension(gtype),fbem_n_nodes(gtype))
    real(kind=real64) :: xis(fbem_n_dimension(gtype))
    real(kind=real64) :: xi(fbem_n_dimension(gtype))
    real(kind=real64) :: xi1d
    real(kind=real64) :: phis(fbem_n_vertices(gtype))
    real(kind=real64) :: phi(fbem_n_nodes(gtype))
    real(kind=real64) :: aux(10)
    ! Obtain the local coordinates of the nodes in the subdivision space
    xis_subdivision_nodes=fbem_xi_hybrid(gtype,0.d0)
    ! Obtain x
    select case (fbem_n_dimension(gtype))
      case (1)
        do k=1,fbem_n_nodes(gtype)
          xis=xis_subdivision_nodes(:,k)
          ! xis -> xi
#         define delta 0.0d0
#         define xi xis(1)
#         define phi phis
#         include <phi_line2.rc>
#         undef delta
#         undef xi
#         undef phi
          xi=0.d0
          do i=1,2
            xi=xi+phis(i)*xi_s(:,i)
          end do
          ! xi -> x
          xi1d=xi(1)
#         define etype gtype
#         define delta 0.0d0
#         define xi xi1d
#         include <phi_1d.rc>
#         undef etype
#         undef delta
#         undef xi
          x_s(:,k)=0.d0
          do i=1,fbem_n_nodes(gtype)
            x_s(:,k)=x_s(:,k)+phi(i)*x(:,i)
          end do
        end do
      case (2)
        select case (fbem_n_vertices(gtype))
          case (3)
            do k=1,fbem_n_nodes(gtype)
              xis=xis_subdivision_nodes(:,k)
              ! xis -> xi
#             define delta 0.0d0
#             define xi xis
#             define phi phis
#             include <phi_tri3.rc>
#             undef delta
#             undef xi
#             undef phi
              xi=0.d0
              do i=1,3
                xi=xi+phis(i)*xi_s(:,i)
              end do
              ! xi -> x
#             define etype gtype
#             define delta 0.0d0
#             include <phi_2d.rc>
#             undef etype
#             undef delta
              x_s(:,k)=0.d0
              do i=1,fbem_n_nodes(gtype)
                x_s(:,k)=x_s(:,k)+phi(i)*x(:,i)
              end do
            end do
          case (4)
            do k=1,fbem_n_nodes(gtype)
              xis=xis_subdivision_nodes(:,k)
              ! xis -> xi
#             define delta 0.0d0
#             define xi xis
#             define phi phis
#             include <phi_quad4.rc>
#             undef delta
#             undef xi
#             undef phi
              xi=0.d0
              do i=1,4
                xi=xi+phis(i)*xi_s(:,i)
              end do
              ! xi -> x
#             define etype gtype
#             define delta 0.0d0
#             include <phi_2d.rc>
#             undef etype
#             undef delta
              x_s(:,k)=0.d0
              do i=1,fbem_n_nodes(gtype)
                x_s(:,k)=x_s(:,k)+phi(i)*x(:,i)
              end do
            end do
        end select
    end select
  end subroutine fbem_obtain_element_subdivision_coordinates

  ! ================================================================================================================================
  ! ELEMENT PROPERTIES
  ! ================================================================================================================================

  !! Calculation of the length of a 1D element in a 2D space
  function fbem_length2d(etype,x_nodes,rule,tol,step)
    implicit none
    !! Element geometrical interpolation type: <tt>etype={fbem_line2,fbem_line3}</tt>.
    integer :: etype
    !! 2D position of element nodes
    real(kind=real64) :: x_nodes(2,fbem_n_nodes(etype))
    !! Selected rule (number of points of Gauss-Legendre quadrature). If <tt>rule<=0</tt> then integration stop criterion is
    !! tolerance. If <tt>rule>0</tt> then integrate with specified number of Gauss-Legendre quadrature points.
    integer :: rule
    !! Integration relative error height. If <tt>tol<=0</tt> then integrate with default tolerance <tt>tol=1.0E-6</tt>.
    real(kind=real64) :: tol
    !! Iterative procedure starts at <tt>rule=1</tt> and does <tt>rule=rule+step</tt> until <tt>tol</tt> is satisfied. If
    !! <tt>step<=0</tt>, default step is <tt>step<=3</tt>.
    integer :: step
    real(kind=real64) :: fbem_length2d
    real(kind=real64) :: xi, length, local_tol, length_old, length_error
    integer :: i, j, local_step
    ! Integration controlled by specified rule
    if ((rule.gt.0).and.(rule.le.gl11_nr)) then
      length=0.0d0
      do i=1,gl11_n(rule)
        xi=gl11_xi(i,rule)
        length=length+fbem_jacobian2d_1d(etype,x_nodes,xi)*gl11_w(i,rule)
      end do
      fbem_length2d=length
    ! Integration controlled by tolerance
    else
      ! If tolerance is zero or less, default tolerance is used
      if (tol.le.0.0d0) then
        local_tol=1.0d-6
      else
        local_tol=tol
      end if
      ! If step is zero or less, or greater than the number of rules-1, default step is used
      if ((step.le.0).or.(step.gt.gl11_nr-1)) then
        local_step=3
      else
        local_step=step
      end if
      ! Initial error greater than tolerance
      length_error=local_tol+1.0d0
      ! First evaluation of length with rule=1
      j=1
      xi=gl11_xi(1,j)
      length=fbem_jacobian2d_1d(etype,x_nodes,xi)*gl11_w(1,j)
      ! While tolerance is not reached increment rule with step
      j=1+local_step
      do while (length_error.ge.local_tol)
        ! save old length
        length_old=length
        ! calculate new length with incremented rule
        length=0.0d0
        do i=1,gl11_n(j)
          xi=gl11_xi(i,j)
          length=length+fbem_jacobian2d_1d(etype,x_nodes,xi)*gl11_w(i,j)
        end do
        ! calculate error of previous rule
        length_error=dabs((length_old-length)/length)
        j=j+local_step
        ! if next rule is greater or equal to the maximum rule, estimate error and length with two last rules
        if (j.ge.gl11_nr) then
          ! calculate length with previous to last rule
          length=0.0d0
          do i=1,gl11_n(gl11_nr-1)
            xi=gl11_xi(i,gl11_nr-1)
            length=length+fbem_jacobian2d_1d(etype,x_nodes,xi)*gl11_w(i,gl11_nr-1)
          end do
          ! save length
          length_old=length
          ! calculate length with last rule
          length=0.0d0
          do i=1,gl11_n(gl11_nr)
            xi=gl11_xi(i,gl11_nr)
            length=length+fbem_jacobian2d_1d(etype,x_nodes,xi)*gl11_w(i,gl11_nr)
          end do
          ! estimate error of previous rule
          length_error=dabs((length_old-length)/length)
          ! if error is less than specified, print a warning message and exit the loop
          if (length_error.ge.local_tol) then
            call fbem_warning_message(error_unit,3,__FILE__,__LINE__,&
                                      'length calculated with less than specified error height.')
          end if
          exit
        end if
      end do
      fbem_length2d=length
    end if
  end function fbem_length2d

  !! Calculation of the length of a 1D element in a 3D space
  function fbem_length3d(etype,x_nodes,rule,tol,step)
    implicit none
    !! Element geometrical interpolation type: <tt>etype={fbem_line2,fbem_line3}</tt>.
    integer :: etype
    !! 3D position of element nodes
    real(kind=real64) :: x_nodes(3,fbem_n_nodes(etype))
    !! Selected rule (number of points of Gauss-Legendre quadrature). If <tt>rule<=0</tt> then integration stop criterion is
    !! tolerance. If <tt>rule>0</tt> then integrate with specified number of Gauss-Legendre quadrature points.
    integer :: rule
    !! Integration relative error height. If <tt>tol<=0</tt> then integrate with default tolerance <tt>tol=1.0E-6</tt>.
    real(kind=real64) :: tol
    !! Iterative procedure starts at <tt>rule=1</tt> and does <tt>rule=rule+step</tt> until <tt>tol</tt> is satisfied. If
    !! <tt>step<=0</tt>, default step is <tt>step<=3</tt>.
    integer :: step
    !! Length of a 1D element in a 3D space
    real(kind=real64) :: fbem_length3d
    real(kind=real64) :: xi, length, local_tol, length_old, length_error
    integer :: i, j, local_step
    ! Integration controlled by specified rule
    if ((rule.gt.0).and.(rule.le.gl11_nr)) then
      length=0.0d0
      do i=1,gl11_n(rule)
        xi=gl11_xi(i,rule)
        length=length+fbem_jacobian3d_1d(etype,x_nodes,xi)*gl11_w(i,rule)
      end do
      fbem_length3d=length
    ! Integration controlled by tolerance
    else
      ! If tolerance is zero or less, default tolerance is used
      if (tol.le.0.0d0) then
        local_tol=1.0d-6
      else
        local_tol=tol
      end if
      ! If step is zero or less, or greater than the number of rules-1, default step is used
      if ((step.le.0).or.(step.gt.gl11_nr-1)) then
        local_step=3
      else
        local_step=step
      end if
      ! Initial error greater than tolerance
      length_error=local_tol+1.0d0
      ! First evaluation of length with rule=1
      j=1
      xi=gl11_xi(1,1)
      length=fbem_jacobian3d_1d(etype,x_nodes,xi)*gl11_w(1,1)
      ! While tolerance is not reached increment rule with step
      j=1+local_step
      do while (length_error.ge.local_tol)
        ! save old length
        length_old=length
        ! calculate new length with incremented rule
        length=0.0d0
        do i=1,gl11_n(j)
          xi=gl11_xi(i,j)
          length=length+fbem_jacobian3d_1d(etype,x_nodes,xi)*gl11_w(i,j)
        end do
        ! calculate error of previous rule
        length_error=dabs((length_old-length)/length)
        j=j+local_step
        ! if next rule is greater or equal to the maximum rule, estimate error and length with two last rules
        if (j.ge.gl11_nr) then
          ! calculate length with previous to last rule
          length=0.0d0
          do i=1,gl11_n(gl11_nr-1)
            xi=gl11_xi(i,gl11_nr-1)
            length=length+fbem_jacobian3d_1d(etype,x_nodes,xi)*gl11_w(i,gl11_nr-1)
          end do
          ! save length
          length_old=length
          ! calculate length with last rule
          length=0.0d0
          do i=1,gl11_n(gl11_nr)
            xi=gl11_xi(i,gl11_nr)
            length=length+fbem_jacobian3d_1d(etype,x_nodes,xi)*gl11_w(i,gl11_nr)
          end do
          ! estimate error of previous rule
          length_error=dabs((length_old-length)/length)
          ! if error is less than specified, print a warning message and exit the loop
          if (length_error.ge.local_tol) then
            call fbem_warning_message(error_unit,3,__FILE__,__LINE__,&
                                      'length calculated with less than specified error height.')
          end if
          exit
        end if
      end do
      fbem_length3d=length
    end if
  end function fbem_length3d

  !! Calculation of the length of a subdivision of a 1D element in a 2D space
  function fbem_length2d_subdivision(etype,x_nodes,xi1,xi2,rule,tol,step)
    implicit none
    ! I/O variables
    real(kind=real64) :: fbem_length2d_subdivision
    integer           :: etype                          !! Element geometrical interpolation type: <tt>etype={fbem_line2,fbem_line3}</tt>.
    real(kind=real64) :: x_nodes(2,fbem_n_nodes(etype)) !! 2D position of element nodes
    real(kind=real64) :: xi1                            !! Minimum xi coordinate
    real(kind=real64) :: xi2                            !! Maximum xi coordinate
    !! Selected rule (number of points of Gauss-Legendre quadrature). If <tt>rule<=0</tt> then integration stop criterion is
    !! tolerance. If <tt>rule>0</tt> then integrate with specified number of Gauss-Legendre quadrature points.
    integer :: rule
    !! Integration relative error height. If <tt>tol<=0</tt> then integrate with default tolerance <tt>tol=1.0E-6</tt>.
    real(kind=real64) :: tol
    !! Iterative procedure starts at <tt>rule=1</tt> and does <tt>rule=rule+step</tt> until <tt>tol</tt> is satisfied. If
    !! <tt>step<=0</tt>, default step is <tt>step<=3</tt>.
    integer :: step
    ! Local variables
    real(kind=real64) :: xi, length, local_tol, length_old, length_error
    integer :: i, j, local_step
    ! Integration controlled by specified rule
    if ((rule.gt.0).and.(rule.le.gl11_nr)) then
      length=0.0d0
      do i=1,gl11_n(rule)
        xi=xi1+0.5d0*(gl11_xi(i,rule)+1.0d0)*(xi2-xi1)
        length=length+fbem_jacobian2d_1d(etype,x_nodes,xi)*0.5d0*(xi2-xi1)*gl11_w(i,rule)
      end do
      fbem_length2d_subdivision=length
    ! Integration controlled by tolerance
    else
      ! If tolerance is zero or less, default tolerance is used
      if (tol.le.0.0d0) then
        local_tol=1.0d-6
      else
        local_tol=tol
      end if
      ! If step is zero or less, or greater than the number of rules-1, default step is used
      if ((step.le.0).or.(step.gt.gl11_nr-1)) then
        local_step=3
      else
        local_step=step
      end if
      ! Initial error greater than tolerance
      length_error=local_tol+1.0d0
      ! First evaluation of length with rule=1
      j=1
      xi=xi1+0.5d0*(gl11_xi(1,j)+1.0d0)*(xi2-xi1)
      length=fbem_jacobian2d_1d(etype,x_nodes,xi)*0.5d0*(xi2-xi1)*gl11_w(1,j)
      ! While tolerance is not reached increment rule with step
      j=1+local_step
      do while (length_error.ge.local_tol)
        ! save old length
        length_old=length
        ! calculate new length with incremented rule
        length=0.0d0
        do i=1,gl11_n(j)
          xi=xi1+0.5d0*(gl11_xi(i,j)+1.0d0)*(xi2-xi1)
          length=length+fbem_jacobian2d_1d(etype,x_nodes,xi)*0.5d0*(xi2-xi1)*gl11_w(i,j)
        end do
        ! calculate error of previous rule
        length_error=dabs((length_old-length)/length)
        j=j+local_step
        ! if next rule is greater or equal to the maximum rule, estimate error and length with two last rules
        if (j.ge.gl11_nr) then
          ! calculate length with previous to last rule
          length=0.0d0
          do i=1,gl11_n(gl11_nr-1)
            xi=xi1+0.5d0*(gl11_xi(i,gl11_nr-1)+1.0d0)*(xi2-xi1)
            length=length+fbem_jacobian2d_1d(etype,x_nodes,xi)*0.5d0*(xi2-xi1)*gl11_w(i,gl11_nr-1)
          end do
          ! save length
          length_old=length
          ! calculate length with last rule
          length=0.0d0
          do i=1,gl11_n(gl11_nr)
            xi=xi1+0.5d0*(gl11_xi(i,gl11_nr)+1.0d0)*(xi2-xi1)
            length=length+fbem_jacobian2d_1d(etype,x_nodes,xi)*0.5d0*(xi2-xi1)*gl11_w(i,gl11_nr)
          end do
          ! estimate error of previous rule
          length_error=dabs((length_old-length)/length)
          ! if error is less than specified, print a warning message and exit the loop
          if (length_error.ge.local_tol) then
            call fbem_warning_message(error_unit,3,__FILE__,__LINE__,&
                                      'length calculated with less than specified error height.')
          end if
          exit
        end if
      end do
      fbem_length2d_subdivision=length
    end if
  end function fbem_length2d_subdivision

  !! Calculation of the area of a 2D element in 2D space.
  function fbem_area2d(etype,x_nodes,rule,tol,step)
    implicit none
    !! Element geometrical interpolation type: <tt>etype={fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    integer :: etype
    !! 2D position of element nodes
    real(kind=real64) :: x_nodes(2,fbem_n_nodes(etype))
    !! Selected rule. If the element is a quadrilateral, it corresponds to the number of points of Gauss-Legendre quadrature (up to
    !! 32). If the element is a triangular, it corresponds to the order of a Wandzura symmetrical quadrature rule (up to 30). If
    !! <tt>rule<=0</tt> then integration stop criterion is tolerance. If <tt>rule>0</tt> then integrate with specified the
    !! specified rule.
    integer :: rule
    !! Integration relative error height. If <tt>tol<=0</tt> then integrate with default tolerance <tt>tol=1.0E-6</tt>.
    real(kind=real64) :: tol
    !! Iterative procedure starts at <tt>rule=1</tt> and does <tt>rule=rule+step</tt> until <tt>tol</tt> is satisfied. If
    !! <tt>step<=0</tt>, default step is <tt>step<=3</tt>.
    integer :: step
    !! Area of a 2D element in 2D space.
    real(kind=real64) :: fbem_area2d
    real(kind=real64) :: xi(2), area, local_tol, area_old, area_error
    integer :: i, j, k, local_step
    ! Switch between triangular (wandzura) or quadrilateral (gauss-legendre) elements
    select case (etype)
      ! Triangular
      case (fbem_tri3,fbem_tri6)
        ! Integration controlled by specified rule
        if ((rule.gt.0).and.(rule.le.wantri_nr)) then
          area=0.0d0
          do i=1,wantri_n(rule)
            xi(1)=wantri_xi1(i,rule)
            xi(2)=wantri_xi2(i,rule)
            area=area+fbem_jacobian2d_2d(etype,x_nodes,xi)*wantri_w(i,rule)
          end do
          fbem_area2d=area
        ! Integration controlled by tolerance
        else
          ! If tolerance is zero or less, default tolerance is used
          if (tol.le.0.0d0) then
            local_tol=1.0d-6
          else
            local_tol=tol
          end if
          ! If step is zero or less, or greater than the number of rules-1, default step is used
          if ((step.le.0).or.(step.gt.wantri_nr-1)) then
            local_step=3
          else
            local_step=step
          end if
          ! Initial error greater than tolerance
          area_error=local_tol+1.0d0
          ! First evaluation of area with rule=1
          j=1
          xi(1)=wantri_xi1(1,1)
          xi(2)=wantri_xi2(1,1)
          area=fbem_jacobian2d_2d(etype,x_nodes,xi)*wantri_w(1,1)
          ! While tolerance is not reached increment rule with step
          j=1+local_step
          do while (area_error.ge.local_tol)
            ! save old area
            area_old=area
            ! calculate new area with incremented rule
            area=0.0d0
            do i=1,wantri_n(j)
              xi(1)=wantri_xi1(i,j)
              xi(2)=wantri_xi2(i,j)
              area=area+fbem_jacobian2d_2d(etype,x_nodes,xi)*wantri_w(i,j)
            end do
            ! calculate error of previous rule
            area_error=dabs((area_old-area)/area)
            j=j+local_step
            ! if next rule is greater or equal to the maximum rule, estimate error and area with two last rules
            if (j.ge.wantri_nr) then
              ! calculate area with previous to last rule
              area=0.0d0
              do i=1,wantri_n(wantri_nr-1)
                xi(1)=wantri_xi1(i,wantri_nr-1)
                xi(2)=wantri_xi2(i,wantri_nr-1)
                area=area+fbem_jacobian2d_2d(etype,x_nodes,xi)*wantri_w(i,wantri_nr-1)
              end do
              ! save area
              area_old=area
              ! calculate area with last rule
              area=0.0d0
              do i=1,wantri_n(wantri_nr-1)
                xi(1)=wantri_xi1(i,wantri_nr-1)
                xi(2)=wantri_xi2(i,wantri_nr-1)
                area=area+fbem_jacobian2d_2d(etype,x_nodes,xi)*wantri_w(i,wantri_nr-1)
              end do
              ! estimate error of previous rule
              area_error=dabs((area_old-area)/area)
              ! if error is less than specified, print a warning message and exit the loop
              if (area_error.ge.local_tol) then
                call fbem_warning_message(error_unit,3,__FILE__,__LINE__,&
                                          'area calculated with less than specified error height.')
              end if
              exit
            end if
          end do
          fbem_area2d=area
        end if
      ! Quadrilateral
      case (fbem_quad4,fbem_quad8,fbem_quad9)
        ! Integration controlled by specified rule
        if ((rule.gt.0).and.(rule.le.gl11_nr)) then
          area=0.0d0
          do i=1,gl11_n(rule)
            xi(1)=gl11_xi(i,rule)
            do j=1,gl11_n(rule)
              xi(2)=gl11_xi(j,rule)
              area=area+fbem_jacobian2d_2d(etype,x_nodes,xi)*gl11_w(i,rule)*gl11_w(j,rule)
            end do
          end do
          fbem_area2d=area
        ! Integration controlled by tolerance
        else
          ! If tolerance is zero or less, default tolerance is used
          if (tol.le.0.0d0) then
            local_tol=1.0d-6
          else
            local_tol=tol
          end if
          ! If step is zero or less, or greater than the number of rules-1, default step is used
          if ((step.le.0).or.(step.gt.gl11_nr-1)) then
            local_step=3
          else
            local_step=step
          end if
          ! Initial error greater than tolerance
          area_error=local_tol+1.0d0
          ! First evaluation of area with rule=1
          k=1
          xi(1)=gl11_xi(1,1)
          xi(2)=gl11_xi(1,1)
          area=fbem_jacobian2d_2d(etype,x_nodes,xi)*gl11_w(1,1)*gl11_w(1,1)
          ! While tolerance is not reached increment rule with step
          k=1+local_step
          do while (area_error.ge.local_tol)
            ! save old area
            area_old=area
            ! calculate new area with incremented rule
            area=0.0d0
            do i=1,gl11_n(k)
              xi(1)=gl11_xi(i,k)
              do j=1,gl11_n(k)
                xi(2)=gl11_xi(j,k)
                area=area+fbem_jacobian2d_2d(etype,x_nodes,xi)*gl11_w(i,k)*gl11_w(j,k)
              end do
            end do
            ! calculate error of previous rule
            area_error=dabs((area_old-area)/area)
            k=k+local_step
            ! if next rule is greater or equal to the maximum rule, estimate error and area with two last rules
            if (k.ge.gl11_nr) then
              ! calculate area with previous to last rule
              area=0.0d0
              do i=1,gl11_n(gl11_nr-1)
                xi(1)=gl11_xi(i,gl11_nr-1)
                do j=1,gl11_n(gl11_nr-1)
                  xi(2)=gl11_xi(j,gl11_nr-1)
                  area=area+fbem_jacobian2d_2d(etype,x_nodes,xi)*gl11_w(i,gl11_nr-1)*gl11_w(j,gl11_nr-1)
                end do
              end do
              ! save previous area
              area_old=area
              ! calculate area with last rule
              area=0.0d0
              do i=1,gl11_n(gl11_nr)
                xi(1)=gl11_xi(i,gl11_nr)
                do j=1,gl11_n(gl11_nr)
                  xi(2)=gl11_xi(j,gl11_nr)
                  area=area+fbem_jacobian2d_2d(etype,x_nodes,xi)*gl11_w(i,gl11_nr)*gl11_w(j,gl11_nr)
                end do
              end do
              ! estimate error of previous rule
              area_error=dabs((area_old-area)/area)
              ! if error is less than specified, print a warning message and exit the loop
              if (area_error.ge.local_tol) then
                call fbem_warning_message(error_unit,3,__FILE__,__LINE__,&
                                          'area calculated with less than specified error height.')
              end if
              exit
            end if
          end do
          fbem_area2d=area
        end if
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={tri3,tri6,quad4,quad8,quad9}')
    end select
  end function fbem_area2d

  !! Calculation of the area of a 2D element in 2D space.
  function fbem_area3d(etype,x_nodes,rule,tol,step)
    implicit none
    !! Element geometrical interpolation type: <tt>etype={fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    integer :: etype
    !! 3D position of element nodes
    real(kind=real64) :: x_nodes(3,fbem_n_nodes(etype))
    !! Selected rule. If the element is a quadrilateral, it corresponds to the number of points of Gauss-Legendre quadrature (up to
    !! 32). If the element is a triangular, it corresponds to the order of a Wandzura symmetrical quadrature rule (up to 30). If
    !! <tt>rule<=0</tt> then integration stop criterion is tolerance. If <tt>rule>0</tt> then integrate with specified the
    !! specified rule.
    integer :: rule
    !! Integration relative error height. If <tt>tol<=0</tt> then integrate with default tolerance <tt>tol=1.0E-6</tt>.
    real(kind=real64) :: tol
    !! Iterative procedure starts at <tt>rule=1</tt> and does <tt>rule=rule+step</tt> until <tt>tol</tt> is satisfied. If
    !! <tt>step<=0</tt>, default step is <tt>step<=3</tt>.
    integer :: step
    !! Area of a 2D element in 2D space.
    real(kind=real64) :: fbem_area3d
    real(kind=real64) :: xi(2), area, local_tol, area_old, area_error
    integer :: i, j, k, local_step
    ! Switch between triangular (dunavant) or quadrilateral (gauss-legendre) elements
    select case (etype)
      ! Triangular
      case (fbem_tri3,fbem_tri6)
        ! Integration controlled by specified rule
        if ((rule.gt.0).and.(rule.le.wantri_nr)) then
          area=0.0d0
          do i=1,wantri_n(rule)
            xi(1)=wantri_xi1(i,rule)
            xi(2)=wantri_xi2(i,rule)
            area=area+fbem_jacobian3d_nd(etype,x_nodes,xi)*wantri_w(i,rule)
          end do
          fbem_area3d=0.5d0*area
        ! Integration controlled by tolerance
        else
          ! If tolerance is zero or less, default tolerance is used
          if (tol.le.0.0d0) then
            local_tol=1.0d-6
          else
            local_tol=tol
          end if
          ! If step is zero or less, or greater than the number of rules-1, default step is used
          if ((step.le.0).or.(step.gt.wantri_nr-1)) then
            local_step=3
          else
            local_step=step
          end if
          ! Initial error greater than tolerance
          area_error=local_tol+1.0d0
          ! First evaluation of area with rule=1
          j=1
          xi(1)=wantri_xi1(1,1)
          xi(2)=wantri_xi2(1,1)
          area=fbem_jacobian3d_nd(etype,x_nodes,xi)*wantri_w(1,1)
          ! While tolerance is not reached increment rule with step
          j=1+local_step
          do while (area_error.ge.local_tol)
            ! save old area
            area_old=area
            ! calculate new area with incremented rule
            area=0.0d0
            do i=1,wantri_n(j)
              xi(1)=wantri_xi1(i,j)
              xi(2)=wantri_xi2(i,j)
              area=area+fbem_jacobian3d_nd(etype,x_nodes,xi)*wantri_w(i,j)
            end do
            ! calculate error of previous rule
            area_error=dabs((area_old-area)/area)
            j=j+local_step
            ! if next rule is greater or equal to the maximum rule, estimate error and area with two last rules
            if (j.ge.wantri_nr) then
              ! calculate area with previous to last rule
              area=0.0d0
              do i=1,wantri_n(wantri_nr-1)
                xi(1)=wantri_xi1(i,wantri_nr-1)
                xi(2)=wantri_xi2(i,wantri_nr-1)
                area=area+fbem_jacobian3d_nd(etype,x_nodes,xi)*wantri_w(i,wantri_nr-1)
              end do
              ! save area
              area_old=area
              ! calculate area with last rule
              area=0.0d0
              do i=1,wantri_n(wantri_nr-1)
                xi(1)=wantri_xi1(i,wantri_nr-1)
                xi(2)=wantri_xi2(i,wantri_nr-1)
                area=area+fbem_jacobian3d_nd(etype,x_nodes,xi)*wantri_w(i,wantri_nr-1)
              end do
              ! estimate error of previous rule
              area_error=dabs(area-area_old)
              ! if error is less than specified, print a warning message and exit the loop
              if (area_error.ge.local_tol) then
                call fbem_warning_message(error_unit,3,__FILE__,__LINE__,&
                                          'area calculated with less than specified error height.')
              end if
              exit
            end if
          end do
          fbem_area3d=area
        end if
      ! Quadrilateral
      case (fbem_quad4,fbem_quad8,fbem_quad9)
        ! Integration controlled by specified rule
        if ((rule.gt.0).and.(rule.le.gl11_nr)) then
          area=0.0d0
          do i=1,gl11_n(rule)
            xi(1)=gl11_xi(i,rule)
            do j=1,gl11_n(rule)
              xi(2)=gl11_xi(j,rule)
              area=area+fbem_jacobian3d_nd(etype,x_nodes,xi)*gl11_w(i,rule)*gl11_w(j,rule)
            end do
          end do
          fbem_area3d=area
        ! Integration controlled by tolerance
        else
          ! If tolerance is zero or less, default tolerance is used
          if (tol.le.0.0d0) then
            local_tol=1.0d-6
          else
            local_tol=tol
          end if
          ! If step is zero or less, or greater than the number of rules-1, default step is used
          if ((step.le.0).or.(step.gt.gl11_nr-1)) then
            local_step=3
          else
            local_step=step
          end if
          ! Initial error greater than tolerance
          area_error=local_tol+1.0d0
          ! First evaluation of area with rule=1
          k=1
          xi(1)=gl11_xi(1,1)
          xi(2)=gl11_xi(1,1)
          area=fbem_jacobian3d_nd(etype,x_nodes,xi)*gl11_w(1,1)*gl11_w(1,1)
          ! While tolerance is not reached increment rule with step
          k=1+local_step
          do while (area_error.ge.local_tol)
            ! save old area
            area_old=area
            ! calculate new area with incremented rule
            area=0.0d0
            do i=1,gl11_n(k)
              xi(1)=gl11_xi(i,k)
              do j=1,gl11_n(k)
                xi(2)=gl11_xi(j,k)
                area=area+fbem_jacobian3d_nd(etype,x_nodes,xi)*gl11_w(i,k)*gl11_w(j,k)
              end do
            end do
            ! calculate error of previous rule
            area_error=dabs((area_old-area)/area)
            k=k+local_step
            ! if next rule is greater or equal to the maximum rule, estimate error and area with two last rules
            if (k.ge.gl11_nr) then
              ! calculate area with previous to last rule
              area=0.0d0
              do i=1,gl11_n(gl11_nr-1)
                xi(1)=gl11_xi(i,gl11_nr-1)
                do j=1,gl11_n(gl11_nr-1)
                  xi(2)=gl11_xi(j,gl11_nr-1)
                  area=area+fbem_jacobian3d_nd(etype,x_nodes,xi)*gl11_w(i,gl11_nr-1)*gl11_w(j,gl11_nr-1)
                end do
              end do
              ! save previous area
              area_old=area
              ! calculate area with last rule
              area=0.0d0
              do i=1,gl11_n(gl11_nr)
                xi(1)=gl11_xi(i,gl11_nr)
                do j=1,gl11_n(gl11_nr)
                  xi(2)=gl11_xi(j,gl11_nr)
                  area=area+fbem_jacobian3d_nd(etype,x_nodes,xi)*gl11_w(i,gl11_nr)*gl11_w(j,gl11_nr)
                end do
              end do
              ! estimate error of previous rule
              area_error=dabs((area_old-area)/area)
              ! if error is less than specified, print a warning message and exit the loop
              if (area_error.ge.local_tol) then
                call fbem_warning_message(error_unit,3,__FILE__,__LINE__,&
                                          'area calculated with less than specified error height.')
              end if
              exit
            end if
          end do
          fbem_area3d=area
        end if
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={tri3,tri6,quad4,quad8,quad9}')
    end select
  end function fbem_area3d

  !! Function that returns the characteristic size an element (approx. maximum length of its edges)
  function fbem_characteristic_length(rn,etype,x_nodes,tol)
    implicit none
    ! I/O
    integer                        :: rn                              !! Spatial dimension
    real(kind=real64)              :: fbem_characteristic_length      !! Characteristic width
    integer                        :: etype                           !! Element type: <tt>etype={fbem_line2,fbem_line3,fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    real(kind=real64)              :: x_nodes(rn,fbem_n_nodes(etype)) !! Position of element nodes
    real(kind=real64)              :: tol                             !! Tolerance
    ! Local
    real(kind=real64)              :: l(fbem_n_edges(etype))
    real(kind=real64), allocatable :: x_nodes_aux(:,:)
    integer                        :: edge, n
    ! Loop through the edges of the element
    do edge=1,fbem_n_edges(etype)
      allocate (x_nodes_aux(rn,fbem_n_nodes(fbem_edge_type(edge,etype))))
      do n=1,fbem_n_nodes(fbem_edge_type(edge,etype))
        x_nodes_aux(:,n)=x_nodes(:,fbem_edge_node(n,edge,etype))
      end do
      select case (rn)
        case (2)
          l(edge)=fbem_length2d(fbem_edge_type(edge,etype),x_nodes_aux,0,tol,0)
        case (3)
          l(edge)=fbem_length3d(fbem_edge_type(edge,etype),x_nodes_aux,0,tol,0)
      end select
      deallocate (x_nodes_aux)
    end do
    fbem_characteristic_length=maxval(l)
  end function fbem_characteristic_length

  !! Function that returns the characteristic size of an element subdivision (approx. maximum length of its edges)
  function fbem_characteristic_length_subdivision(rn,etype,x_nodes,xise,tol)
    implicit none
    ! I/O
    integer                        :: rn                                     !! Spatial dimension
    real(kind=real64)              :: fbem_characteristic_length_subdivision !! Characteristic width
    integer                        :: etype                                  !! Element type: <tt>etype={fbem_line2,fbem_line3,fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    real(kind=real64)              :: x_nodes(rn,fbem_n_nodes(etype))        !! Position of element nodes
    real(kind=real64)              :: xise(2,4)                              !! xi coordinates of the subdivision
    real(kind=real64)              :: tol                                    !! Tolerance
    ! Local
    real(kind=real64)              :: x_nodes_aux(rn,fbem_n_nodes(etype))
    real(kind=real64)              :: xie(fbem_n_dimension(etype),fbem_n_nodes(etype))
    real(kind=real64)              :: xis(fbem_n_dimension(etype))
    real(kind=real64)              :: xise1(1,2)                             !! xi coordinates of the subdivision (temporary)
    real(kind=real64)              :: xise3(2,3)                             !! xi coordinates of the subdivision (temporary)
    integer                        :: n
    ! Obtain the local coordinates of the nodes of the original type of element
    xie=fbem_xi_hybrid(etype,0.d0)
    ! Depending on the type of element
    select case (etype)
      case (fbem_line2,fbem_line3,fbem_line4)
        xise1=xise(1:1,1:2)
        ! Obtain the real coordinates using the element subdivision as a support
        do n=1,fbem_n_nodes(etype)
          xis=fbem_position(1,fbem_line2,xise1,xie(1,n))
          x_nodes_aux(:,n)=fbem_position(rn,etype,x_nodes,xis)
        end do
      case (fbem_tri3,fbem_tri6)
        xise3=xise(:,1:3)
        ! Obtain the real coordinates using the element subdivision as a support
        do n=1,fbem_n_nodes(etype)
          xis=fbem_position(2,fbem_tri3,xise3,xie(:,n))
          x_nodes_aux(:,n)=fbem_position(rn,etype,x_nodes,xis)
        end do
      case (fbem_quad4,fbem_quad8,fbem_quad9)
        ! Obtain the real coordinates using the element subdivision as a support
        do n=1,fbem_n_nodes(etype)
          xis=fbem_position(2,fbem_quad4,xise,xie(:,n))
          x_nodes_aux(:,n)=fbem_position(rn,etype,x_nodes,xis)
        end do
    end select
    fbem_characteristic_length_subdivision=fbem_characteristic_length(rn,etype,x_nodes_aux,tol)
  end function fbem_characteristic_length_subdivision

  !! Calculate the centroid of the element with a given number of integration points
  function fbem_element_centroid_ngp(rn,type_g,x_nodes,esize,glp)
    implicit none
    ! I/O
    integer                      :: rn                               !! R^n
    integer                      :: type_g                           !! Geometrical interpolation
    real(kind=real64)            :: x_nodes(rn,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
    real(kind=real64)            :: esize                            !! Size: length for a line element, area for a surface element, ...
    integer                      :: glp                              !! Number of Gauss points for the integration
    real(kind=real64)            :: fbem_element_centroid_ngp(rn)    !! Centroid
    ! Local
    real(kind=real64)            :: xmoment(rn)                      ! First moment
    integer                      :: kphi                             ! Counter variable for shape functions loops
    integer                      :: kc                               ! Counter for coordinate loops
    integer                      :: nnodes_g                         ! Number of nodes of the element
    integer                      :: kip, kip1, kip2                  ! Counter variable of integration points
    integer                      :: rule                             ! Integration rule
    real(kind=real64)            :: xi, xi2d(2)                      ! Coordinate xi
    real(kind=real64)            :: w, w2d(2)                        ! Weights of an integration point
    real(kind=real64)            :: aux(10)                          ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: phi_g(fbem_n_nodes(type_g))      ! Geometrical shape functions values
    real(kind=real64)            :: dphidxi_g(fbem_n_nodes(type_g))  ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: dphidxi1_g(fbem_n_nodes(type_g)) ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: dphidxi2_g(fbem_n_nodes(type_g)) ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: x(rn)                            ! Position vector of the Gauss point
    real(kind=real64)            :: T(rn)                            ! Tangent vector at xi
    real(kind=real64)            :: T1(rn), T2(rn), T1f(3), T2f(3)   ! Tangent vectors at xi
    real(kind=real64)            :: N(3)                             ! Normal vector at xi
    real(kind=real64)            :: jg                               ! Geometric jacobian
    real(kind=real64)            :: jw                               ! Jacobian * weight
    ! Number of nodes of the element
    nnodes_g=fbem_n_nodes(type_g)
    ! Initialize integrals
    xmoment=0.d0
    ! Switch depending on the number of dimensions of the element
    select case (fbem_n_dimension(type_g))
      !
      ! LINE ELEMENTS
      !
      case (1)
        ! Integration rule
        rule=glp
        ! Loop through integrations points
        do kip=1,gl11_n(rule)
          ! xi coordinate and weight
          xi=gl11_xi(kip,rule)
          w=gl11_w(kip,rule)
          ! Geometrical shape functions and first derivatives at xi
#         define etype type_g
#         define delta 0.0d0
#         define phi phi_g
#         define dphidxi dphidxi_g
#         include <phi_and_dphidxi_1d.rc>
#         undef etype
#         undef delta
#         undef phi
#         undef dphidxi
          ! Components calculation of x and T at xi
          x=0.d0
          T=0.d0
          do kphi=1,nnodes_g
            x=x+phi_g(kphi)*x_nodes(:,kphi)
            T=T+dphidxi_g(kphi)*x_nodes(:,kphi)
          end do
          ! Geometric jacobian
          jg=sqrt(dot_product(T,T))
          ! Add integration points
          jw=jg*w
          xmoment(:)=xmoment(:)+x(:)*jw
        end do ! Loop through integrations points
      !
      ! AREA ELEMENTS
      !
      case (2)
        select case (fbem_n_edges(type_g))
          !
          ! TRIANGULAR ELEMENTS
          !
          case (3)
            ! Integration rule
            rule=2*glp-1
            ! Loop through integrations points
            do kip=1,wantri_n(rule)
              ! xi1, xi2 coordinates and weight
              xi2d(1)=wantri_xi1(kip,rule)
              xi2d(2)=wantri_xi2(kip,rule)
              w=wantri_w(kip,rule)
              ! Geometrical shape functions and first derivatives at xi
#             define etype type_g
#             define delta 0.0d0
#             define xi xi2d
#             define phi phi_g
#             define dphidxi1 dphidxi1_g
#             define dphidxi2 dphidxi2_g
#             include <phi_and_dphidxik_2d.rc>
#             undef etype
#             undef delta
#             undef xi
#             undef phi
#             undef dphidxi1
#             undef dphidxi2
              ! Components calculation of x and T at xi
              x=0.d0
              T1=0.d0
              T2=0.d0
              do kphi=1,nnodes_g
                x=x+phi_g(kphi)*x_nodes(:,kphi)
                T1=T1+dphidxi1_g(kphi)*x_nodes(:,kphi)
                T2=T2+dphidxi2_g(kphi)*x_nodes(:,kphi)
              end do
              T1f=0.d0
              T2f=0.d0
              do kc=1,rn
                T1f(kc)=T1(kc)
                T2f(kc)=T2(kc)
              end do
              ! Normal vector as T1 x T2 at xi
              N(1)=T1f(2)*T2f(3)-T1f(3)*T2f(2)
              N(2)=T1f(3)*T2f(1)-T1f(1)*T2f(3)
              N(3)=T1f(1)*T2f(2)-T1f(2)*T2f(1)
              ! Geometric jacobian
              jg=sqrt(dot_product(N,N))
              ! Add integration point
              jw=jg*w
              xmoment(:)=xmoment(:)+x(:)*jw
            end do
          !
          ! QUADRILATERAL ELEMENTS
          !
          case (4)
            ! Integration rule
            rule=glp
            ! Loop through integrations points
            do kip1=1,gl11_n(rule)
              ! xi1 coordinate and weight
              xi2d(1)=gl11_xi(kip1,rule)
              w2d(1)=gl11_w(kip1,rule)
              do kip2=1,gl11_n(rule)
                ! xi2 coordinate and weight
                xi2d(2)=gl11_xi(kip2,rule)
                w2d(2)=gl11_w(kip2,rule)
                ! Geometrical shape functions and first derivatives at xi
#               define etype type_g
#               define delta 0.0d0
#               define xi xi2d
#               define phi phi_g
#               define dphidxi1 dphidxi1_g
#               define dphidxi2 dphidxi2_g
#               include <phi_and_dphidxik_2d.rc>
#               undef etype
#               undef delta
#               undef xi
#               undef phi
#               undef dphidxi1
#               undef dphidxi2
                ! Components calculation of x and T at xi
                x=0.d0
                T1=0.d0
                T2=0.d0
                do kphi=1,nnodes_g
                  x=x+phi_g(kphi)*x_nodes(:,kphi)
                  T1=T1+dphidxi1_g(kphi)*x_nodes(:,kphi)
                  T2=T2+dphidxi2_g(kphi)*x_nodes(:,kphi)
                end do
                T1f=0.d0
                T2f=0.d0
                do kc=1,rn
                  T1f(kc)=T1(kc)
                  T2f(kc)=T2(kc)
                end do
                ! Normal vector as T1 x T2 at xi
                N(1)=T1f(2)*T2f(3)-T1f(3)*T2f(2)
                N(2)=T1f(3)*T2f(1)-T1f(1)*T2f(3)
                N(3)=T1f(1)*T2f(2)-T1f(2)*T2f(1)
                ! Geometric jacobian
                jg=sqrt(dot_product(N,N))
                ! Add integration point
                jw=jg*w2d(1)*w2d(2)
                xmoment(:)=xmoment(:)+x(:)*jw
              end do
            end do
        end select
      !
      ! VOLUME ELEMENTS
      !
      case (3)
        xmoment=0.d0
        write(*,*) 'warning: function fbem_element_centroid_ngp : volume elements not implemented yet'
    end select
    ! Calculate centroid
    fbem_element_centroid_ngp=xmoment/esize
  end function fbem_element_centroid_ngp

  !! Calculate the size of the element with a given relative error
  function fbem_element_centroid(rn,type_g,x_nodes,esize,rerror)
    implicit none
    ! I/O
    integer                      :: rn                               !! R^n
    integer                      :: type_g                           !! Geometrical interpolation
    real(kind=real64)            :: x_nodes(rn,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
    real(kind=real64)            :: esize                            !! Size: length for a line element, area for a surface element, volume for volume element
    real(kind=real64)            :: rerror                           !! Relative error
    real(kind=real64)            :: fbem_element_centroid(rn)        !! Centroid
    ! Local
    integer                      :: glp                              ! Number of Gauss points for the integration
    real(kind=real64)            :: cl                               ! Characteristic length (diameter of the element)
    real(kind=real64)            :: centroid_old(rn)
    real(kind=real64)            :: centroid_new(rn)
    logical                      :: ok
    real(kind=real64)            :: rerror_new
    ! Check
    if (rerror.le.0.d0) then
      write(error_unit,*) 'rerror = ',rerror
      call fbem_error_message(error_unit,0,'fbem_element_centroid',0,'relative error must be >0')
    end if
    ! Iterative procedure
    glp=1
    centroid_old=fbem_element_centroid_ngp(rn,type_g,x_nodes,esize,glp)
    cl=fbem_characteristic_length(rn,type_g,x_nodes,rerror)
    ok=.false.
    do while (ok.eqv.(.false.))
      glp=glp+1
      centroid_new=fbem_element_centroid_ngp(rn,type_g,x_nodes,esize,glp)
      rerror_new=sqrt(dot_product(centroid_new-centroid_old,centroid_new-centroid_old))/cl
      if (rerror_new.lt.rerror) then
        ok=.true.
      else
        if (glp.eq.15) then
          write(error_unit,*) 'Required relative error = ', rerror
          write(error_unit,*) 'Achieved relative error = ', rerror_new
          call fbem_warning_message(error_unit,0,'fbem_element_centroid',0,'>30th polynomial order required to calculate the centroid.')
          ok=.true.
        end if
      end if
      centroid_old=centroid_new
    end do
    fbem_element_centroid=centroid_new
  end function fbem_element_centroid

  !! Calculate the size of the element with a given relative error
  function fbem_element_size(rn,type_g,x_nodes,rerror)
    implicit none
    ! I/O
    integer                      :: rn                               !! R^n
    integer                      :: type_g                           !! Geometrical interpolation
    real(kind=real64)            :: x_nodes(rn,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
    real(kind=real64)            :: rerror                           !! Relative error
    real(kind=real64)            :: fbem_element_size                !! Size: length for a line element, area for a surface element, volume...
    ! Local
    integer                      :: glp                              ! Number of Gauss points for the integration
    real(kind=real64)            :: size_old
    real(kind=real64)            :: size_new
    real(kind=real64)            :: rerror_new
    logical                      :: ok
    ! Check
    if (rerror.le.0.d0) then
      write(error_unit,*) 'rerror = ',rerror
      call fbem_error_message(error_unit,0,'fbem_element_centroid',0,'relative error must be >0')
    end if
    ! Iterative procedure
    glp=1
    size_old=fbem_element_size_ngp(rn,type_g,x_nodes,glp)
    ok=.false.
    do while (ok.eqv.(.false.))
      glp=glp+1
      size_new=fbem_element_size_ngp(rn,type_g,x_nodes,glp)
      rerror_new=abs(size_old-size_new)/size_new
      if (rerror_new.lt.rerror) then
        ok=.true.
      else
        if (glp.eq.15) then
          write(error_unit,*) 'Required relative error = ', rerror
          write(error_unit,*) 'Achieved relative error = ', rerror_new
          call fbem_warning_message(error_unit,0,'fbem_element_size',0,'>27th polynomial order required to calculate the size.')
          ok=.true.
        end if
      end if
      size_old=size_new
    end do
    fbem_element_size=size_new
  end function fbem_element_size

  !! Calculate the size of the element with a given number of integration points
  function fbem_element_size_ngp(rn,type_g,x_nodes,glp)
    implicit none
    ! I/O
    integer                      :: rn                               !! R^n
    integer                      :: type_g                           !! Geometrical interpolation
    real(kind=real64)            :: x_nodes(rn,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
    integer                      :: glp                              !! Number of Gauss points for the integration
    real(kind=real64)            :: fbem_element_size_ngp            !! Size: length for a line element, area for a surface element, ...
    ! Local
    integer                      :: kphi                             ! Counter variable for shape functions loops
    integer                      :: kc                               ! Counter for coordinate loops
    integer                      :: nnodes_g                         ! Number of nodes of the element
    integer                      :: kip, kip1, kip2                  ! Counter variable of integration points
    integer                      :: rule                             ! Integration rule
    real(kind=real64)            :: xi, xi2d(2)                      ! Coordinate xi
    real(kind=real64)            :: w, w2d(2)                        ! Weights of an integration point
    real(kind=real64)            :: aux(10)                          ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: dphidxi_g(fbem_n_nodes(type_g))  ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: dphidxi1_g(fbem_n_nodes(type_g)) ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: dphidxi2_g(fbem_n_nodes(type_g)) ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: T(rn)                            ! Tangent vector at xi
    real(kind=real64)            :: T1(rn), T2(rn), T1f(3), T2f(3)   ! Tangent vectors at xi
    real(kind=real64)            :: N(3)                             ! Normal vector at xi
    real(kind=real64)            :: jg                               ! Geometric jacobian
    ! Number of nodes of the element
    nnodes_g=fbem_n_nodes(type_g)
    ! Initialize integrals
    fbem_element_size_ngp=0.d0
    ! Switch depending on the number of dimensions of the element
    select case (fbem_n_dimension(type_g))
      !
      ! LINE ELEMENTS
      !
      case (1)
        ! Integration rule
        rule=glp
        ! Loop through integrations points
        do kip=1,gl11_n(rule)
          ! xi coordinate and weight
          xi=gl11_xi(kip,rule)
          w=gl11_w(kip,rule)
          ! Geometrical shape functions and first derivatives at xi
#         define etype type_g
#         define delta 0.0d0
#         define dphidxi dphidxi_g
#         include <dphidxi_1d.rc>
#         undef etype
#         undef delta
#         undef dphidxi
          ! Components calculation of x and T at xi
          T=0.d0
          do kphi=1,nnodes_g
            T=T+dphidxi_g(kphi)*x_nodes(:,kphi)
          end do
          ! Geometric jacobian
          jg=sqrt(dot_product(T,T))
          ! Add integration points
          fbem_element_size_ngp=fbem_element_size_ngp+jg*w
        end do ! Loop through integrations points
      !
      ! AREA ELEMENTS
      !
      case (2)
        select case (fbem_n_edges(type_g))
          !
          ! TRIANGULAR ELEMENTS
          !
          case (3)
            ! Integration rule
            rule=2*glp-1
            ! Loop through integrations points
            do kip=1,wantri_n(rule)
              ! xi1, xi2 coordinates and weight
              xi2d(1)=wantri_xi1(kip,rule)
              xi2d(2)=wantri_xi2(kip,rule)
              w=wantri_w(kip,rule)
              ! Geometrical shape functions and first derivatives at xi
#             define etype type_g
#             define delta 0.0d0
#             define xi xi2d
#             define dphidxi1 dphidxi1_g
#             define dphidxi2 dphidxi2_g
#             include <dphidxik_2d.rc>
#             undef etype
#             undef delta
#             undef xi
#             undef dphidxi1
#             undef dphidxi2
              ! Components calculation of x and T at xi
              T1=0.d0
              T2=0.d0
              do kphi=1,nnodes_g
                T1=T1+dphidxi1_g(kphi)*x_nodes(:,kphi)
                T2=T2+dphidxi2_g(kphi)*x_nodes(:,kphi)
              end do
              T1f=0.d0
              T2f=0.d0
              do kc=1,rn
                T1f(kc)=T1(kc)
                T2f(kc)=T2(kc)
              end do
              ! Normal vector as T1 x T2 at xi
              N(1)=T1f(2)*T2f(3)-T1f(3)*T2f(2)
              N(2)=T1f(3)*T2f(1)-T1f(1)*T2f(3)
              N(3)=T1f(1)*T2f(2)-T1f(2)*T2f(1)
              ! Geometric jacobian
              jg=sqrt(dot_product(N,N))
              ! Add integration point
              fbem_element_size_ngp=fbem_element_size_ngp+jg*w
            end do
          !
          ! QUADRILATERAL ELEMENTS
          !
          case (4)
            ! Integration rule
            rule=glp
            ! Loop through integrations points
            do kip1=1,gl11_n(rule)
              ! xi1 coordinate and weight
              xi2d(1)=gl11_xi(kip1,rule)
              w2d(1)=gl11_w(kip1,rule)
              do kip2=1,gl11_n(rule)
                ! xi2 coordinate and weight
                xi2d(2)=gl11_xi(kip2,rule)
                w2d(2)=gl11_w(kip2,rule)
                ! Geometrical first derivatives at xi
#               define etype type_g
#               define delta 0.0d0
#               define xi xi2d
#               define dphidxi1 dphidxi1_g
#               define dphidxi2 dphidxi2_g
#               include <dphidxik_2d.rc>
#               undef etype
#               undef delta
#               undef xi
#               undef dphidxi1
#               undef dphidxi2
                ! Components calculation of x and T at xi
                T1=0.d0
                T2=0.d0
                do kphi=1,nnodes_g
                  T1=T1+dphidxi1_g(kphi)*x_nodes(:,kphi)
                  T2=T2+dphidxi2_g(kphi)*x_nodes(:,kphi)
                end do
                T1f=0.d0
                T2f=0.d0
                do kc=1,rn
                  T1f(kc)=T1(kc)
                  T2f(kc)=T2(kc)
                end do
                ! Normal vector as T1 x T2 at xi
                N(1)=T1f(2)*T2f(3)-T1f(3)*T2f(2)
                N(2)=T1f(3)*T2f(1)-T1f(1)*T2f(3)
                N(3)=T1f(1)*T2f(2)-T1f(2)*T2f(1)
                ! Geometric jacobian
                jg=sqrt(dot_product(N,N))
                ! Add integration point
                fbem_element_size_ngp=fbem_element_size_ngp+jg*w2d(1)*w2d(2)
              end do
            end do
        end select
      !
      ! VOLUME ELEMENTS
      !
      case (3)
        write(*,*) 'warning: function fbem_element_size_ngp : volume elements not implemented yet'
        fbem_element_size_ngp=1.d0
    end select
  end function fbem_element_size_ngp

  !! Calculate a simplified version of an element by definining its centre (centroid) and radius (radius of a ball that contains the
  !! element). The element covering ball is approximate.
  subroutine fbem_geometry_element_ball(rn,type_g,x_nodes,glp,centre,radius)
    implicit none
    ! I/O
    integer                      :: rn                               !! R^n
    integer                      :: type_g                           !! Geometrical interpolation
    real(kind=real64)            :: x_nodes(rn,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
    integer                      :: glp                              !! Number of Gauss points for the integration
    real(kind=real64)            :: centre(rn)                       !! Ball centre
    real(kind=real64)            :: radius                           !! Ball radius
    ! Local
    real(kind=real64)            :: esize                            ! Element size: length for a line element, area for a area element, volume for a volume element
    real(kind=real64)            :: xmoment(rn)                      ! First moment
    integer                      :: kphi                             ! Counter variable for shape functions loops
    integer                      :: kc                               ! Counter for coordinate loops
    integer                      :: nnodes_g                         ! Number of nodes of the element
    integer                      :: kip, kip1, kip2                  ! Counter variable of integration points
    integer                      :: rule                             ! Integration rule
    real(kind=real64)            :: xi, xi2d(2)                      ! Coordinate xi
    real(kind=real64)            :: w, w2d(2)                        ! Weights of an integration point
    real(kind=real64)            :: aux(10)                          ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: phi_g(fbem_n_nodes(type_g))      ! Geometrical shape functions values
    real(kind=real64)            :: dphidxi_g(fbem_n_nodes(type_g))  ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: dphidxi1_g(fbem_n_nodes(type_g)) ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: dphidxi2_g(fbem_n_nodes(type_g)) ! Geometrical shape functions first derivatives values
    real(kind=real64)            :: x(rn)                            ! Position vector of the Gauss point
    real(kind=real64)            :: r(rn)                            ! Distance vector
    real(kind=real64)            :: T(rn)                            ! Tangent vector at xi
    real(kind=real64)            :: T1(rn), T2(rn), T1f(3), T2f(3)   ! Tangent vectors at xi
    real(kind=real64)            :: N(3)                             ! Normal vector at xi
    real(kind=real64)            :: jg                               ! Geometric jacobian
    real(kind=real64)            :: jw                               ! Jacobian * weight
    real(kind=real64)            :: tmp_radius                       ! Temporary radius
    ! Number of nodes of the element
    nnodes_g=fbem_n_nodes(type_g)
    ! Initialize integrals
    esize=0.d0
    xmoment=0.d0
    ! Switch depending on the number of dimensions of the element
    select case (fbem_n_dimension(type_g))
      !
      ! LINE ELEMENTS
      !
      case (1)
        ! Integration rule
        rule=glp
        ! Loop through integrations points
        do kip=1,gl11_n(rule)
          ! xi coordinate and weight
          xi=gl11_xi(kip,rule)
          w=gl11_w(kip,rule)
          ! Geometrical shape functions and first derivatives at xi
#         define etype type_g
#         define delta 0.0d0
#         define phi phi_g
#         define dphidxi dphidxi_g
#         include <phi_and_dphidxi_1d.rc>
#         undef etype
#         undef delta
#         undef phi
#         undef dphidxi
          ! Components calculation of x and T at xi
          x=0.d0
          T=0.d0
          do kphi=1,nnodes_g
            x=x+phi_g(kphi)*x_nodes(:,kphi)
            T=T+dphidxi_g(kphi)*x_nodes(:,kphi)
          end do
          ! Geometric jacobian
          jg=sqrt(dot_product(T,T))
          ! Add integration points
          jw=jg*w
          esize=esize+jw
          xmoment(:)=xmoment(:)+x(:)*jw
        end do ! Loop through integrations points
      !
      ! AREA ELEMENTS
      !
      case (2)
        select case (fbem_n_edges(type_g))
          !
          ! TRIANGULAR ELEMENTS
          !
          case (3)
            ! Integration rule
            rule=2*glp-1
            ! Loop through integrations points
            do kip=1,wantri_n(rule)
              ! xi1, xi2 coordinates and weight
              xi2d(1)=wantri_xi1(kip,rule)
              xi2d(2)=wantri_xi2(kip,rule)
              w=wantri_w(kip,rule)
              ! Geometrical shape functions and first derivatives at xi
#             define etype type_g
#             define delta 0.0d0
#             define xi xi2d
#             define phi phi_g
#             define dphidxi1 dphidxi1_g
#             define dphidxi2 dphidxi2_g
#             include <phi_and_dphidxik_2d.rc>
#             undef etype
#             undef delta
#             undef xi
#             undef phi
#             undef dphidxi1
#             undef dphidxi2
              ! Components calculation of x and T at xi
              x=0.d0
              T1=0.d0
              T2=0.d0
              do kphi=1,nnodes_g
                x=x+phi_g(kphi)*x_nodes(:,kphi)
                T1=T1+dphidxi1_g(kphi)*x_nodes(:,kphi)
                T2=T2+dphidxi2_g(kphi)*x_nodes(:,kphi)
              end do
              T1f=0.d0
              T2f=0.d0
              do kc=1,rn
                T1f(kc)=T1(kc)
                T2f(kc)=T2(kc)
              end do
              ! Normal vector as T1 x T2 at xi
              N(1)=T1f(2)*T2f(3)-T1f(3)*T2f(2)
              N(2)=T1f(3)*T2f(1)-T1f(1)*T2f(3)
              N(3)=T1f(1)*T2f(2)-T1f(2)*T2f(1)
              ! Geometric jacobian
              jg=sqrt(dot_product(N,N))
              ! Add integration point
              jw=jg*w
              esize=esize+jw
              xmoment(:)=xmoment(:)+x(:)*jw
            end do
          !
          ! QUADRILATERAL ELEMENTS
          !
          case (4)
            ! Integration rule
            rule=glp
            ! Loop through integrations points
            do kip1=1,gl11_n(rule)
              ! xi1 coordinate and weight
              xi2d(1)=gl11_xi(kip1,rule)
              w2d(1)=gl11_w(kip1,rule)
              do kip2=1,gl11_n(rule)
                ! xi2 coordinate and weight
                xi2d(2)=gl11_xi(kip2,rule)
                w2d(2)=gl11_w(kip2,rule)
                ! Geometrical shape functions and first derivatives at xi
#               define etype type_g
#               define delta 0.0d0
#               define xi xi2d
#               define phi phi_g
#               define dphidxi1 dphidxi1_g
#               define dphidxi2 dphidxi2_g
#               include <phi_and_dphidxik_2d.rc>
#               undef etype
#               undef delta
#               undef xi
#               undef phi
#               undef dphidxi1
#               undef dphidxi2
                ! Components calculation of x and T at xi
                x=0.d0
                T1=0.d0
                T2=0.d0
                do kphi=1,nnodes_g
                  x=x+phi_g(kphi)*x_nodes(:,kphi)
                  T1=T1+dphidxi1_g(kphi)*x_nodes(:,kphi)
                  T2=T2+dphidxi2_g(kphi)*x_nodes(:,kphi)
                end do
                T1f=0.d0
                T2f=0.d0
                do kc=1,rn
                  T1f(kc)=T1(kc)
                  T2f(kc)=T2(kc)
                end do
                ! Normal vector as T1 x T2 at xi
                N(1)=T1f(2)*T2f(3)-T1f(3)*T2f(2)
                N(2)=T1f(3)*T2f(1)-T1f(1)*T2f(3)
                N(3)=T1f(1)*T2f(2)-T1f(2)*T2f(1)
                ! Geometric jacobian
                jg=sqrt(dot_product(N,N))
                ! Add integration point
                jw=jg*w2d(1)*w2d(2)
                esize=esize+jw
                xmoment(:)=xmoment(:)+x(:)*jw
              end do
            end do
        end select
      !
      ! VOLUME ELEMENTS
      !
      case (3)
        stop 'subroutine fbem_geometry_element_ball : volume elements not implemented yet'
    end select
    ! Calculate centroid
    centre=xmoment/esize
    !
    ! Calculate an approximate radius using the nodes
    !
    radius=0.d0
    do kphi=1,nnodes_g
      r=x_nodes(:,kphi)-centre
      tmp_radius=sqrt(dot_product(r,r))
      if (radius.lt.tmp_radius) then
        radius=tmp_radius
      end if
    end do
    !
    ! Calculate an approximate radius using gauss points
    !
    ! Switch depending on the number of dimensions of the element
    select case (fbem_n_dimension(type_g))
      !
      ! LINE ELEMENTS
      !
      case (1)
        ! Integration rule
        rule=glp
        ! Loop through integrations points
        do kip=1,gl11_n(rule)
          ! xi coordinate and weight
          xi=gl11_xi(kip,rule)
          w=gl11_w(kip,rule)
          ! Geometrical shape functions and first derivatives at xi
#         define etype type_g
#         define delta 0.0d0
#         define phi phi_g
#         include <phi_1d.rc>
#         undef etype
#         undef delta
#         undef phi
          ! Components calculation of x and T at xi
          x=0.d0
          do kphi=1,nnodes_g
            x=x+phi_g(kphi)*x_nodes(:,kphi)
          end do
          r=x-centre
          tmp_radius=sqrt(dot_product(r,r))
          if (radius.lt.tmp_radius) then
            radius=tmp_radius
          end if
        end do ! Loop through integrations points
      !
      ! AREA ELEMENTS
      !
      case (2)
        select case (fbem_n_edges(type_g))
          !
          ! TRIANGULAR ELEMENTS
          !
          case (3)
            ! Integration rule
            rule=2*glp-1
            ! Loop through integrations points
            do kip=1,wantri_n(rule)
              ! xi1, xi2 coordinates and weight
              xi2d(1)=wantri_xi1(kip,rule)
              xi2d(2)=wantri_xi2(kip,rule)
              w=wantri_w(kip,rule)
              ! Geometrical shape functions and first derivatives at xi
#             define etype type_g
#             define delta 0.0d0
#             define xi xi2d
#             define phi phi_g
#             include <phi_2d.rc>
#             undef etype
#             undef delta
#             undef xi
#             undef phi
              ! Components calculation of x and T at xi
              x=0.d0
              do kphi=1,nnodes_g
                x=x+phi_g(kphi)*x_nodes(:,kphi)
              end do
              r=x-centre
              tmp_radius=sqrt(dot_product(r,r))
              if (radius.lt.tmp_radius) then
                radius=tmp_radius
              end if
            end do
          !
          ! QUADRILATERAL ELEMENTS
          !
          case (4)
            ! Integration rule
            rule=glp
            ! Loop through integrations points
            do kip1=1,gl11_n(rule)
              ! xi1 coordinate and weight
              xi2d(1)=gl11_xi(kip1,rule)
              w2d(1)=gl11_w(kip1,rule)
              do kip2=1,gl11_n(rule)
                ! xi2 coordinate and weight
                xi2d(2)=gl11_xi(kip2,rule)
                w2d(2)=gl11_w(kip2,rule)
                ! Geometrical shape functions and first derivatives at xi
#               define etype type_g
#               define delta 0.0d0
#               define xi xi2d
#               define phi phi_g
#               include <phi_2d.rc>
#               undef etype
#               undef delta
#               undef xi
#               undef phi
                ! Components calculation of x and T at xi
                x=0.d0
                do kphi=1,nnodes_g
                  x=x+phi_g(kphi)*x_nodes(:,kphi)
                end do
                r=x-centre
                tmp_radius=sqrt(dot_product(r,r))
                if (radius.lt.tmp_radius) then
                  radius=tmp_radius
                end if
              end do
            end do
        end select
      !
      ! VOLUME ELEMENTS
      !
      case (3)
        stop 'subroutine fbem_geometry_element_ball : volume elements not implemented yet'
    end select
  end subroutine fbem_geometry_element_ball

  ! ==================================================================
  ! NEAREST REFERENCE COORDINATE OF AN ELEMENT WITH RESPECT TO A POINT
  ! ==================================================================
  !
  ! There are 4 methods:
  !   (1) Using the element ball. Computationally cheap, but with poor precision.
  !   (2) Using the element nodes. Computationally cheap (less than (1)), and with poor precision.
  !   (3) Using a sampling algorithm. Computationally expensive, but robust and with good precision.
  !   (4) Using a minimization algorithm (using a geometry linearization). Cheap with excellent precision, but do not always converge.
  !     It converges always when the point belongs to the element and when the element plane and straight. For curved elements,
  !     it converges near the element.
  !
  ! An adaptative algorithm for BEM requirements that uses (2), (3) and (4) is also implemented. Being cl the characteristic length,
  ! the algorithm does the following:
  !   Use the element nodes to compute a first approximation to barxi and rmin.
  !   Calculate d=rmin/cl
  !   If d<=1
  !     Use the minimization algorithm to compute a numerically exact value of barxi and rmin.
  !     If minimization algorithm does not converge
  !       Use the sampling algorithm
  !     End if
  !   End if


  !! Calculation of <tt>xi</tt> or <tt>xi_1,xi_2</tt> nearest reference coordinates of an element with respect to a point, using
  !! only the element nodes.
  subroutine fbem_nearest_xi_nodes(rn,etype,x_nodes,p,nearest_xi,r_nearest_xi)
    implicit none
    ! I/O
    integer           :: rn                                   !! Dimensional space
    integer           :: etype                                !! Element type: <tt>etype={fbem_line2,fbem_line3,fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    real(kind=real64) :: x_nodes(rn,fbem_n_nodes(etype))      !! Position vectors of the element nodes
    real(kind=real64) :: p(rn)                                !! Point
    real(kind=real64) :: nearest_xi(fbem_n_dimension(etype))  !! Nearest <tt>xi</tt> coordinates (reference space)
    real(kind=real64) :: r_nearest_xi                         !! Minimum distance
    ! Local
    real(kind=real64) :: rmod, r(rn)
    integer           :: kn, nearest_node
    !
    ! First guess with first node
    !
    ! Calculate distance vector
    r=x_nodes(:,1)-p
    rmod=dsqrt(dot_product(r,r))
    ! Assign the node 1 as the nearest
    r_nearest_xi=rmod
    nearest_node=1
    !
    ! Compare with the other nodes
    !
    do kn=2,fbem_n_nodes(etype)
      ! Calculate distance vector
      r=x_nodes(:,kn)-p
      rmod=dsqrt(dot_product(r,r))
      ! Compare
      if (rmod.lt.r_nearest_xi) then
        r_nearest_xi=rmod
        nearest_node=kn
      end if
    end do
    ! Obtain the xi reference coordinate
    select case (fbem_n_dimension(etype))
      case (1)
#       define node nearest_node
#       define xi nearest_xi(1)
#       define delta 0.0d0
#       include <xi_1d_at_node.rc>
#       undef node
#       undef xi
#       undef delta
      case (2)
#       define node nearest_node
#       define xi nearest_xi
#       define delta 0.0d0
#       include <xi_2d_at_node.rc>
#       undef node
#       undef xi
#       undef delta
    end select
  end subroutine fbem_nearest_xi_nodes

  ! ===================
  ! SAMPLING ALGORITHMS
  ! ===================

  !! Calculation of the <tt>xi</tt> nearest reference coordinate of an 1D element with respect to a point.
  !!
  !! The algorithm skeleton is:
  !!   - Sampling stage. <tt>n_samples-1</tt> interior points are sampled, and nearest <tt>xi</tt> is selected from all sampled points
  !!     including extreme points which are always evaluated.
  !!   - Resampling stage. <tt>n_resamples</tt> points are evaluated around previous minimum found. This stage is repeated
  !!     <tt>n_iterations</tt> times.
  subroutine fbem_nearest_xi(rn,etype,x_nodes,p,n_samples,n_resamples,n_iterations,nearest_xi,r_nearest_xi)
    implicit none
    ! I/O
    integer           :: rn                              !! Dimensional space
    integer           :: etype                           !! Element geometrical interpolation type: <tt>etype={fbem_line2,fbem_line3}</tt>.
    real(kind=real64) :: x_nodes(rn,fbem_n_nodes(etype)) !! Position of element nodes
    real(kind=real64) :: p(rn)                           !! Point
    integer           :: n_samples                       !! Number of samples, if <tt>n_samples<1</tt> default is used <tt>n_samples=25</tt>.
    integer           :: n_resamples                     !! Number of resamples, if <tt>n_resamples<1</tt> default is used <tt>n_samples=14</tt>.
    integer           :: n_iterations                    !! Number of iterations of resampling, if <tt>n_iterations<0</tt> default is used <tt>n_samples=2</tt>.
    real(kind=real64) :: nearest_xi                      !! Nearest <tt>xi</tt> coordinate (reference space)
    real(kind=real64) :: r_nearest_xi                    !! Minimum distance (minimum module of distance vector)
    ! Local
    integer           :: local_n_samples
    integer           :: local_n_resamples
    integer           :: local_n_iterations
    real(kind=real64) :: phi_vector(fbem_n_nodes(etype)), aux(10)
    real(kind=real64) :: x(rn), rv(rn), r
    real(kind=real64) :: xi, xi_min, r_min, xi_min_old, r_min_old, width
    integer           :: j, k, l
    ! Copy procedure parameters
    local_n_samples=n_samples
    local_n_resamples=n_resamples
    local_n_iterations=n_iterations
    ! SAMPLING
    ! if n_samples<1 then default is used (25)
    if (local_n_samples.lt.1) then
      local_n_samples=25
    end if
    ! loop through sample points
    do k=0,local_n_samples+1
      xi=-1.0d0+2.0d0/dble(local_n_samples+1)*dble(k)
      ! Shape functions at xi
#     define phi phi_vector
#     define delta 0.0d0
#     include <phi_1d.rc>
#     undef delta
#     undef phi
      ! Calculate position vector and distance vector
      x=0.0d0
      do j=1,fbem_n_nodes(etype)
        x=x+phi_vector(j)*x_nodes(:,j)
      end do
      ! Distance vector
      rv=x-p
      r=dot_product(rv,rv)
      r=sqrt(r)
      ! In the first iteration assign the minimum
      if (k.eq.0) then
        xi_min=xi
        r_min=r
      end if
      ! Check if is less than minimum
      if (r.lt.r_min) then
        xi_min=xi
        r_min=r
      end if
    end do
    ! Save sampling results
    xi_min_old=xi_min
    r_min_old=r_min
    ! RESAMPLING
    ! if n_resamples<1 then default is used (14)
    if (local_n_resamples.lt.1) then
      local_n_resamples=14
    end if
    ! if n_iterations<0 then default is used (2)
    if (local_n_iterations.lt.0) then
      local_n_iterations=2
    end if
    ! Search width
    width=2.0d0/dble(local_n_samples+1)
    ! If iterations > 0 do
    if (local_n_iterations.gt.0) then
      ! loop through iterations
      do l=1,local_n_iterations
        ! loop through resampling points
        do k=1,local_n_resamples
          xi=xi_min_old-width+(2.0d0*width)/dble(local_n_resamples+1)*dble(k)
          ! Only sample points in the element
          if ((xi.gt.-1.0d0).and.(xi.lt.1.0d0)) then
            ! Shape functions at xi
#           define phi phi_vector
#           define delta 0.0d0
#           include <phi_1d.rc>
#           undef delta
#           undef phi
            ! Calculate position vector and distance vector
            x=0.0d0
            do j=1,fbem_n_nodes(etype)
              x=x+phi_vector(j)*x_nodes(:,j)
            end do
            ! Distance vector
            rv=x-p
            r=dot_product(rv,rv)
            r=sqrt(r)
            ! Check if is less than minimum
            if (r.lt.r_min) then
              xi_min=xi
              r_min=r
            end if
          end if
        end do
        ! Stablish new center and width
        xi_min_old=xi_min
        r_min_old=r_min
        width=2.0d0*width/dble(local_n_resamples+1)
      end do
    end if
    ! Save results to output variables
    nearest_xi=xi_min
    r_nearest_xi=r_min
  end subroutine fbem_nearest_xi

  !! Calculation of the <tt>xi</tt> nearest reference coordinate of a line element with respect to a point, limiting the searching
  !! to a subdivision of the element.
  !!
  !! The algorithm skeleton is:
  !!   - Sampling stage. <tt>n_samples-1</tt> interior points are sampled, and nearest <tt>xi</tt> is selected from all sampled
  !!     points including extreme points which are always evaluated.
  !!   - Resampling stage. <tt>n_resamples</tt> points are evaluated around previous minimum found. This stage is repeated
  !!     <tt>n_iterations</tt> times.
  !!
  !!  Future modification: change nearest xi to nearest xis
  subroutine fbem_nearest_xi_subdivision(rn,etype,x_nodes,xis_min,xis_max,p,n_samples,n_resamples,n_iterations,nearest_xi,r_nearest_xi)
    implicit none
    integer           :: rn
    integer           :: etype                           !! Element type: <tt>etype={fbem_line2,fbem_line3}</tt>.
    real(kind=real64) :: x_nodes(rn,fbem_n_nodes(etype)) !! Position of element nodes
    real(kind=real64) :: xis_min, xis_max                !! Subdivision of the element
    real(kind=real64) :: p(rn)                           !! Point
    integer           :: n_samples                       !! Number of samples, if <tt>n_samples<1</tt> default is used <tt>n_samples=25</tt>.
    integer           :: n_resamples                     !! Number of resamples, if <tt>n_resamples<1</tt> default is used <tt>n_resamples=14</tt>.
    integer           :: n_iterations                    !! Number of iterations of resampling, if <tt>n_iterations<0</tt> default is used <tt>n_samples=2</tt>.
    real(kind=real64) :: nearest_xi                      !! Nearest <tt>xi</tt> coordinate (reference space)
    real(kind=real64) :: r_nearest_xi                    !! Minimum distance (minimum module of distance vector)
    ! Local
    integer           :: local_n_samples
    integer           :: local_n_resamples
    integer           :: local_n_iterations
    real(kind=real64) :: phi_vector(fbem_n_nodes(etype)), aux(10)
    real(kind=real64) :: x(rn), rv(rn), r
    real(kind=real64) :: xi, xi_min, r_min, xi_min_old, r_min_old, width
    integer           :: j, k, l
    ! Copy procedure parameters
    local_n_samples=n_samples
    local_n_resamples=n_resamples
    local_n_iterations=n_iterations
    ! SAMPLING
    ! if n_samples<1 then default is used (25)
    if (local_n_samples.lt.1) then
      local_n_samples=25
    end if
    ! loop through sample points
    do k=0,local_n_samples+1
      xi=xis_min+(xis_max-xis_min)/dble(local_n_samples+1)*dble(k)
      ! Shape functions at xi
#     define phi phi_vector
#     define delta 0.0d0
#     include <phi_1d.rc>
#     undef delta
#     undef phi
      ! Calculate position vector and distance vector
      x=0.0d0
      do j=1,fbem_n_nodes(etype)
        x=x+phi_vector(j)*x_nodes(:,j)
      end do
      ! Distance vector
      rv=x-p
      r=dot_product(rv,rv)
      r=sqrt(r)
      ! In the first iteration assign minimum
      if (k.eq.0) then
        xi_min=xi
        r_min=r
      end if
      ! Check if is less than minimum
      if (r.lt.r_min) then
        xi_min=xi
        r_min=r
      end if
    end do
    ! Save sampling results
    xi_min_old=xi_min
    r_min_old=r_min
    ! RESAMPLING
    ! if n_resamples<1 then default is used (14)
    if (local_n_resamples.lt.1) then
      local_n_resamples=14
    end if
    ! if n_iterations<0 then default is used (2)
    if (local_n_iterations.lt.0) then
      local_n_iterations=2
    end if
    ! Search width
    width=(xis_max-xis_min)/dble(local_n_samples+1)
    ! If iterations > 0 do
    if (local_n_iterations.gt.0) then
      ! loop through iterations
      do l=1,local_n_iterations
        ! loop through resampling points
        do k=1,local_n_resamples
          xi=xi_min_old-width+(2.0d0*width)/dble(local_n_resamples+1)*dble(k)
          ! Only sample points in the element
          if ((xi.gt.xis_min).and.(xi.lt.xis_max)) then
            ! Shape functions at xi
#           define phi phi_vector
#           define delta 0.0d0
#           include <phi_1d.rc>
#           undef delta
#           undef phi
            ! Calculate position vector and distance vector
            x=0.0d0
            do j=1,fbem_n_nodes(etype)
              x=x+phi_vector(j)*x_nodes(:,j)
            end do
            ! Distance vector
            rv=x-p
            r=dot_product(rv,rv)
            r=sqrt(r)
            ! Check if is less than minimum
            if (r.lt.r_min) then
              xi_min=xi
              r_min=r
            end if
          end if
        end do
        ! Stablish new center and width
        xi_min_old=xi_min
        r_min_old=r_min
        width=2.0d0*width/dble(local_n_resamples+1)
      end do
    end if
    ! Save results to output variables
    nearest_xi=xi_min
    r_nearest_xi=r_min
  end subroutine fbem_nearest_xi_subdivision

  !! Calculation of <tt>xi_1,xi_2</tt> nearest reference coordinates of an element with respect to a point
  !!
  !! The algorithm skeleton is:
  !!   - Sampling stage. <tt>n_samples-1</tt> interior points are sampled, and nearest <tt>xi</tt> is selected from all sampled
  !!     points including extreme points which are always evaluated.
  !!   - Resampling stage. <tt>n_resamples</tt> points are evaluated around previous minimum found. This stage is repeated
  !!     <tt>n_iterations</tt> times.
  subroutine fbem_nearest_xi1xi2(rn,etype,x_nodes,p,n_samples,n_resamples,n_iterations,nearest_xi,r_nearest_xi)
    implicit none
    ! I/O
    integer           :: rn                              !! Dimensional space
    integer           :: etype                           !! Element type: <tt>etype={fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    real(kind=real64) :: x_nodes(rn,fbem_n_nodes(etype)) !! Position of the element nodes
    real(kind=real64) :: p(rn)                           !! Point
    integer           :: n_samples                       !! Number of samples, if <tt>n_samples<1</tt> default is used <tt>n_samples=25</tt>.
    integer           :: n_resamples                     !! Number of resamples, if <tt>n_resamples<1</tt> default is used <tt>n_samples=14</tt>.
    integer           :: n_iterations                    !! Number of iterations of resampling, if <tt>n_iterations<0</tt> default is used <tt>n_samples=2</tt>.
    real(kind=real64) :: nearest_xi(2)                   !! Nearest <tt>xi_1,xi_2</tt> coordinates (reference space)
    real(kind=real64) :: r_nearest_xi                    !! Minimum distance (minimum module of distance vector)
    ! Local
    integer           :: local_n_samples
    integer           :: local_n_resamples
    integer           :: local_n_iterations
    real(kind=real64) :: phi_g(fbem_n_nodes(etype)), aux(10)
    real(kind=real64) :: x(rn), rv(rn), r
    real(kind=real64) :: xi(2), xi_min(2), r_min, xi_min_old(2), r_min_old, width
    integer           :: j, k1, k2, l
    ! Copy procedure parameters
    local_n_samples=n_samples
    local_n_resamples=n_resamples
    local_n_iterations=n_iterations
    ! SAMPLING
    ! if n_samples<1 then default is used (25)
    if (local_n_samples.lt.1) then
      local_n_samples=25
    end if
    select case (etype)
      ! loop through sample points for triangular elements
      case(fbem_tri3,fbem_tri6)
        do k1=0,local_n_samples+1
          xi(1)=dble(k1)/dble(local_n_samples+1)
          do k2=0,local_n_samples+1-k1
            xi(2)=dble(k2)/dble(local_n_samples+1)
            ! Shape functions at xi
#           define phi phi_g
#           define delta 0.0d0
#           include <phi_2d.rc>
#           undef phi
#           undef delta
            ! Calculate position vector and distance vector
            x=0.d0
            do j=1,fbem_n_nodes(etype)
              x=x+phi_g(j)*x_nodes(:,j)
            end do
            ! Distance vector
            rv=x-p
            r=dot_product(rv,rv)
            r=sqrt(r)
            ! The first sampling point is assumed initially as the minimum
            if ((k1.eq.0).and.(k2.eq.0)) then
              xi_min(1)=xi(1)
              xi_min(2)=xi(2)
              r_min=r
            end if
            ! Check if is less than minimum
            if (r.lt.r_min) then
              xi_min(1)=xi(1)
              xi_min(2)=xi(2)
              r_min=r
            end if
          end do
        end do
      ! loop through sample points for quadrilateral elements
      case(fbem_quad4,fbem_quad8,fbem_quad9)
        ! loop through sample points for quadrilateral elements
        do k1=0,local_n_samples+1
          xi(1)=-1.0d0+2.0d0*dble(k1)/dble(local_n_samples+1)
          do k2=0,local_n_samples+1
            xi(2)=-1.0d0+2.0d0*dble(k2)/dble(local_n_samples+1)
            ! Shape functions at xi
#           define phi phi_g
#           define delta 0.0d0
#           include <phi_2d.rc>
#           undef phi
#           undef delta
            ! Calculate position vector and distance vector
            x=0.d0
            do j=1,fbem_n_nodes(etype)
              x=x+phi_g(j)*x_nodes(:,j)
            end do
            ! Distance vector
            rv=x-p
            r=dot_product(rv,rv)
            r=sqrt(r)
            ! The first sampling point is assumed initially as the minimum
            if ((k1.eq.0).and.(k2.eq.0)) then
              xi_min(1)=xi(1)
              xi_min(2)=xi(2)
              r_min=r
            end if
            ! Check if is less than minimum
            if (r.lt.r_min) then
              xi_min(1)=xi(1)
              xi_min(2)=xi(2)
              r_min=r
            end if
          end do
        end do
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={tri3,tri6,quad4,quad8,quad9}')
    end select
    ! Save sampling results
    xi_min_old(1)=xi_min(1)
    xi_min_old(2)=xi_min(2)
    r_min_old=r_min
    ! RESAMPLING
    ! if n_resamples<1 then default is used (14)
    if (local_n_resamples.lt.1) then
      local_n_resamples=14
    end if
    ! if n_iterations<0 then default is used (2)
    if (local_n_iterations.lt.0) then
      local_n_iterations=2
    end if
    ! If iterations>0 do
    if (local_n_iterations.gt.0) then
      select case (etype)
        ! loop through resampling points for triangular elements
        case(fbem_tri3,fbem_tri6)
          ! Resampling search width
          width=1.0d0/dble(local_n_samples+1)
          ! loop through iterations
          do l=1,local_n_iterations
            do k1=1,local_n_resamples
              xi(1)=xi_min_old(1)+2.0d0*(dble(k1)-dble(local_n_resamples+1)/2.0d0)/dble(local_n_resamples+1)*width
              do k2=1,local_n_resamples
                xi(2)=xi_min_old(2)+2.0d0*(dble(k2)-dble(local_n_resamples+1)/2.0d0)/dble(local_n_resamples+1)*width
                ! Only sample points in the element
                if ((xi(1).gt.0.d0).and.(xi(2).gt.0.d0).and.((xi(1)+xi(2)).lt.1.d0)) then
                  ! Shape functions at xi
#                 define delta 0.0d0
#                 define phi phi_g
#                 include <phi_2d.rc>
#                 undef phi
#                 undef delta
                  ! Calculate position vector and distance vector
                  x=0.d0
                  do j=1,fbem_n_nodes(etype)
                    x=x+phi_g(j)*x_nodes(:,j)
                  end do
                  ! Distance vector
                  rv=x-p
                  r=dot_product(rv,rv)
                  r=sqrt(r)
                  ! Check if is less than minimum
                  if (r.lt.r_min) then
                    xi_min(1)=xi(1)
                    xi_min(2)=xi(2)
                    r_min=r
                  end if
                end if
              end do
            end do
            ! Stablish new center and width
            xi_min_old(1)=xi_min(1)
            xi_min_old(2)=xi_min(2)
            r_min_old=r_min
            width=2.0d0*width/dble(local_n_resamples+1)
          end do
        ! loop through resampling points for quadrilateral elements
        case(fbem_quad4,fbem_quad8,fbem_quad9)
          ! Resampling search width
          width=2.0d0/dble(n_samples+1)
          ! loop through iterations
          do l=1,local_n_iterations
            do k1=1,local_n_resamples
              xi(1)=xi_min_old(1)+2.0d0*(dble(k1)-dble(local_n_resamples+1)/2.0d0)/dble(local_n_resamples+1)*width
              if ((xi(1).gt.-1.0d0).and.(xi(1).lt.1.0d0)) then
                do k2=1,local_n_resamples
                  xi(2)=xi_min_old(2)+2.0d0*(dble(k2)-dble(local_n_resamples+1)/2.0d0)/dble(local_n_resamples+1)*width
                  ! Only sample points in the element
                  if ((xi(2).gt.-1.0d0).and.(xi(2).lt.1.0d0)) then
                    ! Shape functions at xi
#                   define phi phi_g
#                   define delta 0.0d0
#                   include <phi_2d.rc>
#                   undef phi
#                   undef delta
                    ! Calculate position vector and distance vector
                    x=0.d0
                    do j=1,fbem_n_nodes(etype)
                      x=x+phi_g(j)*x_nodes(:,j)
                    end do
                    ! Distance vector
                    rv=x-p
                    r=dot_product(rv,rv)
                    r=sqrt(r)
                    ! Check if is less than minimum
                    if (r.lt.r_min) then
                      xi_min(1)=xi(1)
                      xi_min(2)=xi(2)
                      r_min=r
                    end if
                  end if
                end do
              end if
            end do
            ! Stablish new center and width
            xi_min_old(1)=xi_min(1)
            xi_min_old(2)=xi_min(2)
            r_min_old=r_min
            width=width/dble(local_n_resamples+1)
            width=2.0d0*width/dble(local_n_resamples+1)
          end do
      end select
    end if
    ! Save results to output variables
    nearest_xi(1)=xi_min(1)
    nearest_xi(2)=xi_min(2)
    r_nearest_xi=r_min
  end subroutine fbem_nearest_xi1xi2

  !! Calculation of <tt>xis_1,xis_2</tt> nearest reference coordinates of a subdivision of an element with respect to a point.
  !!
  !! The algorithm skeleton is:
  !!   - Sampling stage. <tt>n_samples-1</tt> interior points are sampled, and nearest <tt>xi</tt> is selected from all sampled
  !!     points including extreme points which are always evaluated.
  !!   - Resampling stage. <tt>n_resamples</tt> points are evaluated around previous minimum found. This stage is repeated
  !!     <tt>n_iterations</tt> times.
  subroutine fbem_nearest_xi1xi2_subdivision(rn,etype,x_nodes,xise,p,n_samples,n_resamples,n_iterations,nearest_xis,r_nearest_xis)
    implicit none
    ! I/O
    integer           :: rn                              !! Dimensional space
    integer           :: etype                           !! Element type: <tt>etype={fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    real(kind=real64) :: x_nodes(rn,fbem_n_nodes(etype)) !! Position of the element nodes
    real(kind=real64) :: xise(2,4)                       !! Element coordinates of the subdivision
    real(kind=real64) :: p(rn)                           !! Point
    integer           :: n_samples                       !! Number of samples, if <tt>n_samples<1</tt> default is used <tt>n_samples=25</tt>.
    integer           :: n_resamples                     !! Number of resamples, if <tt>n_resamples<1</tt> default is used <tt>n_samples=14</tt>.
    integer           :: n_iterations                    !! Number of iterations of resampling, if <tt>n_iterations<0</tt> default is used <tt>n_samples=2</tt>.
    real(kind=real64) :: nearest_xis(2)                  !! Nearest <tt>xis_1,xis_2</tt> coordinates (reference space)
    real(kind=real64) :: r_nearest_xis                   !! Minimum distance (minimum module of distance vector)
    ! Local
    integer           :: local_n_samples
    integer           :: local_n_resamples
    integer           :: local_n_iterations
    real(kind=real64) :: phi_g(fbem_n_nodes(etype)), aux(10)
    real(kind=real64) :: x(rn), rv(rn), r
    real(kind=real64) :: xi(2), xis(2), xis_min(2), r_min, xis_min_old(2), r_min_old, width
    integer           :: j, k1, k2, l
    ! Copy procedure parameters
    local_n_samples=n_samples
    local_n_resamples=n_resamples
    local_n_iterations=n_iterations
    ! SAMPLING
    ! if n_samples<1 then default is used (25)
    if (local_n_samples.lt.1) then
      local_n_samples=25
    end if
    select case (etype)
      ! loop through sample points for triangular elements
      case(fbem_tri3,fbem_tri6)
        do k1=0,local_n_samples+1
          xis(1)=dble(k1)/dble(local_n_samples+1)
          do k2=0,local_n_samples+1-k1
            xis(2)=dble(k2)/dble(local_n_samples+1)
            ! Transformation: xis->xi
            ! Shape functions at xis
#           define xi xis
#           define phi phi_g
#           define delta 0.0d0
#           include <phi_tri3.rc>
#           undef xi
#           undef phi
#           undef delta
            ! Calculate xi
            xi=0.d0
            do j=1,3
              xi=xi+phi_g(j)*xise(:,j)
            end do
            ! Shape functions at xi
#           define phi phi_g
#           define delta 0.0d0
#           include <phi_2d.rc>
#           undef phi
#           undef delta
            ! Calculate position vector and distance vector
            x=0.d0
            do j=1,fbem_n_nodes(etype)
              x=x+phi_g(j)*x_nodes(:,j)
            end do
            ! Distance vector
            rv=x-p
            r=dot_product(rv,rv)
            r=sqrt(r)
            ! The first sampling point is assumed initially as the minimum
            if ((k1.eq.0).and.(k2.eq.0)) then
              xis_min(1)=xis(1)
              xis_min(2)=xis(2)
              r_min=r
            end if
            ! Check if is less than minimum
            if (r.lt.r_min) then
              xis_min(1)=xis(1)
              xis_min(2)=xis(2)
              r_min=r
            end if
          end do
        end do
      ! loop through sample points for quadrilateral elements
      case(fbem_quad4,fbem_quad8,fbem_quad9)
        ! loop through sample points for quadrilateral elements
        do k1=0,local_n_samples+1
          xis(1)=-1.0d0+2.0d0*dble(k1)/dble(local_n_samples+1)
          do k2=0,local_n_samples+1
            xis(2)=-1.0d0+2.0d0*dble(k2)/dble(local_n_samples+1)
            ! Transformation: xis->xi
            ! Shape functions at xis
#           define xi xis
#           define phi phi_g
#           define delta 0.0d0
#           include <phi_quad4.rc>
#           undef xi
#           undef phi
#           undef delta
            ! Calculate xi
            xi=0.d0
            do j=1,4
              xi=xi+phi_g(j)*xise(:,j)
            end do
            ! Shape functions at xi
#           define phi phi_g
#           define delta 0.0d0
#           include <phi_2d.rc>
#           undef phi
#           undef delta
            ! Calculate position vector and distance vector
            x=0.d0
            do j=1,fbem_n_nodes(etype)
              x=x+phi_g(j)*x_nodes(:,j)
            end do
            ! Distance vector
            rv=x-p
            r=dot_product(rv,rv)
            r=sqrt(r)
            ! The first sampling point is assumed initially as the minimum
            if ((k1.eq.0).and.(k2.eq.0)) then
              xis_min(1)=xis(1)
              xis_min(2)=xis(2)
              r_min=r
            end if
            ! Check if is less than minimum
            if (r.lt.r_min) then
              xis_min(1)=xis(1)
              xis_min(2)=xis(2)
              r_min=r
            end if
          end do
        end do
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={tri3,tri6,quad4,quad8,quad9}')
    end select
    ! Save sampling results
    xis_min_old(1)=xis_min(1)
    xis_min_old(2)=xis_min(2)
    r_min_old=r_min
    ! RESAMPLING
    ! if n_resamples<1 then default is used (14)
    if (local_n_resamples.lt.1) then
      local_n_resamples=14
    end if
    ! if n_iterations<0 then default is used (2)
    if (local_n_iterations.lt.0) then
      local_n_iterations=2
    end if
    ! If iterations>0 do
    if (local_n_iterations.gt.0) then
      select case (etype)
        ! loop through resampling points for triangular elements
        case(fbem_tri3,fbem_tri6)
          ! Resampling search width
          width=1.0d0/dble(local_n_samples+1)
          ! loop through iterations
          do l=1,local_n_iterations
            do k1=1,local_n_resamples
              xis(1)=xis_min_old(1)+2.0d0*(dble(k1)-dble(local_n_resamples+1)/2.0d0)/dble(local_n_resamples+1)*width
              do k2=1,local_n_resamples
                xis(2)=xis_min_old(2)+2.0d0*(dble(k2)-dble(local_n_resamples+1)/2.0d0)/dble(local_n_resamples+1)*width
                ! Only sample points in the element
                if ((xis(1).gt.0.d0).and.(xis(2).gt.0.d0).and.((xis(1)+xis(2)).lt.1.d0)) then
                  ! Transformation: xis->xi
                  ! Shape functions at xis
#                 define xi xis
#                 define phi phi_g
#                 define delta 0.0d0
#                 include <phi_tri3.rc>
#                 undef xi
#                 undef phi
#                 undef delta
                  ! Calculate xi
                  xi=0.d0
                  do j=1,3
                    xi=xi+phi_g(j)*xise(:,j)
                  end do
                  ! Shape functions at xi
#                 define delta 0.0d0
#                 define phi phi_g
#                 include <phi_2d.rc>
#                 undef phi
#                 undef delta
                  ! Calculate position vector and distance vector
                  x=0.d0
                  do j=1,fbem_n_nodes(etype)
                    x=x+phi_g(j)*x_nodes(:,j)
                  end do
                  ! Distance vector
                  rv=x-p
                  r=dot_product(rv,rv)
                  r=sqrt(r)
                  ! Check if is less than minimum
                  if (r.lt.r_min) then
                    xis_min(1)=xis(1)
                    xis_min(2)=xis(2)
                    r_min=r
                  end if
                end if
              end do
            end do
            ! Stablish new center and width
            xis_min_old(1)=xis_min(1)
            xis_min_old(2)=xis_min(2)
            r_min_old=r_min
            width=2.0d0*width/dble(local_n_resamples+1)
          end do
        ! loop through resampling points for quadrilateral elements
        case(fbem_quad4,fbem_quad8,fbem_quad9)
          ! Resampling search width
          width=2.0d0/dble(local_n_samples+1)
          ! loop through iterations
          do l=1,local_n_iterations
            do k1=1,local_n_resamples
              xis(1)=xis_min_old(1)+2.0d0*(dble(k1)-dble(local_n_resamples+1)/2.0d0)/dble(local_n_resamples+1)*width
              if ((xi(1).gt.-1.0d0).and.(xi(1).lt.1.0d0)) then
                do k2=1,local_n_resamples
                  xis(2)=xis_min_old(2)+2.0d0*(dble(k2)-dble(local_n_resamples+1)/2.0d0)/dble(local_n_resamples+1)*width
                  ! Only sample points in the element
                  if ((xis(2).gt.-1.0d0).and.(xis(2).lt.1.0d0)) then
                    ! Transformation: xis->xi
                    ! Shape functions at xis
#                   define xi xis
#                   define phi phi_g
#                   define delta 0.0d0
#                   include <phi_quad4.rc>
#                   undef xi
#                   undef phi
#                   undef delta
                    ! Calculate xi
                    xi=0.d0
                    do j=1,4
                      xi=xi+phi_g(j)*xise(:,j)
                    end do
                    ! Shape functions at xi
#                   define phi phi_g
#                   define delta 0.0d0
#                   include <phi_2d.rc>
#                   undef phi
#                   undef delta
                    ! Calculate position vector and distance vector
                    x=0.d0
                    do j=1,fbem_n_nodes(etype)
                      x=x+phi_g(j)*x_nodes(:,j)
                    end do
                    ! Distance vector
                    rv=x-p
                    r=dot_product(rv,rv)
                    r=sqrt(r)
                    ! Check if is less than minimum
                    if (r.lt.r_min) then
                      xis_min(1)=xis(1)
                      xis_min(2)=xis(2)
                      r_min=r
                    end if
                  end if
                end do
              end if
            end do
            ! Stablish new center and width
            xis_min_old(1)=xis_min(1)
            xis_min_old(2)=xis_min(2)
            r_min_old=r_min
            width=width/dble(local_n_resamples+1)
            width=2.0d0*width/dble(local_n_resamples+1)
          end do
      end select
    end if
    ! Save results to output variables
    nearest_xis(1)=xis_min(1)
    nearest_xis(2)=xis_min(2)
    r_nearest_xis=r_min
  end subroutine fbem_nearest_xi1xi2_subdivision

  ! =======================
  ! MINIMIZATION ALGORITHMS
  ! =======================

  recursive subroutine nearest_minimization_iteration_1d(rn,gtype,x,x_i,barxi,new_barxi)
    implicit none
    ! I/O
    integer            :: rn
    integer            :: gtype
    real(kind=real64)  :: x(rn,fbem_n_nodes(gtype))
    real(kind=real64)  :: x_i(rn)
    real(kind=real64)  :: barxi
    real(kind=real64)  :: new_barxi
    ! Local
    integer            :: k
    real(kind=real128) :: phi(fbem_n_nodes(gtype))
    real(kind=real128) :: dphidxi(fbem_n_nodes(gtype))
    real(kind=real128) :: aux(10)
    real(kind=real128) :: x_barxi(rn), dxdxi_barxi(rn)
    real(kind=real128) :: a(rn), b(rn)
    ! Calculate
#   define etype gtype
#   define xi barxi
#   define delta 0.d0
#   include <phi_and_dphidxi_1d.rc>
#   undef etype
#   undef xi
#   undef delta
    x_barxi=0.d0
    dxdxi_barxi=0.d0
    do k=1,fbem_n_nodes(gtype)
      x_barxi=x_barxi+phi(k)*x(:,k)
      dxdxi_barxi=dxdxi_barxi+dphidxi(k)*x(:,k)
    end do
    a=x_barxi-dxdxi_barxi*barxi-x_i
    b=dxdxi_barxi
    new_barxi=-dot_product(a,b)/dot_product(b,b)
  end subroutine nearest_minimization_iteration_1d

  recursive subroutine nearest_minimization_iteration_2d(rn,gtype,x,x_i,barxi,new_barxi)
    implicit none
    ! I/O
    integer            :: rn
    integer            :: gtype
    real(kind=real64)  :: x(rn,fbem_n_nodes(gtype))
    real(kind=real64)  :: x_i(rn)
    real(kind=real64)  :: barxi(2)
    real(kind=real64)  :: new_barxi(2)
    ! Local
    integer            :: i
    real(kind=real128) :: phi(fbem_n_nodes(gtype))
    real(kind=real128) :: dphidxi1(fbem_n_nodes(gtype))
    real(kind=real128) :: dphidxi2(fbem_n_nodes(gtype))
    real(kind=real128) :: aux(10)
    real(kind=real128) :: x_barxi(rn), dxdxi1_barxi(rn), dxdxi2_barxi(rn)
    real(kind=real128) :: a(rn), b(rn), c(rn), bb, cc, bc, ab, ac, detA
#   define etype gtype
#   define xi barxi
#   define delta 0.d0
#   include <phi_and_dphidxik_2d.rc>
#   undef etype
#   undef xi
#   undef delta
    x_barxi=0.d0
    dxdxi1_barxi=0.d0
    dxdxi2_barxi=0.d0
    do i=1,fbem_n_nodes(gtype)
      x_barxi     =x_barxi     +phi(i)*x(:,i)
      dxdxi1_barxi=dxdxi1_barxi+dphidxi1(i)*x(:,i)
      dxdxi2_barxi=dxdxi2_barxi+dphidxi2(i)*x(:,i)
    end do
    a=x_barxi-dxdxi1_barxi*barxi(1)-dxdxi2_barxi*barxi(2)-x_i
    b=dxdxi1_barxi
    c=dxdxi2_barxi
    bb=dot_product(b,b)
    cc=dot_product(c,c)
    bc=dot_product(b,c)
    ab=dot_product(a,b)
    ac=dot_product(a,c)
    detA=bb*cc-bc**2
    new_barxi(1)=-(ab*cc-bc*ac)/detA
    new_barxi(2)=-(bb*ac-ab*bc)/detA
  end subroutine nearest_minimization_iteration_2d

  recursive subroutine nearest_minimization_iteration_3d_shell(gtype,x,v3_md,tv3_md,x_i,barxi,new_barxi)
    implicit none
    ! I/O
    integer            :: gtype
    real(kind=real64)  :: x(3,fbem_n_nodes(gtype))
    real(kind=real64)  :: v3_md(3,fbem_n_nodes(gtype))
    real(kind=real64)  :: tv3_md(fbem_n_nodes(gtype))
    real(kind=real64)  :: x_i(3)
    real(kind=real64)  :: barxi(3)
    real(kind=real64)  :: new_barxi(3)
    ! Local
    integer            :: i, j
    real(kind=real128) :: phi(fbem_n_nodes(gtype)), varphi(fbem_n_nodes(gtype))
    real(kind=real128) :: dphidxi1(fbem_n_nodes(gtype)), dvarphidxi1(fbem_n_nodes(gtype))
    real(kind=real128) :: dphidxi2(fbem_n_nodes(gtype)), dvarphidxi2(fbem_n_nodes(gtype))
    real(kind=real128) :: dphidxi3(fbem_n_nodes(gtype)), dvarphidxi3(fbem_n_nodes(gtype))
    real(kind=real128) :: aux(10)
    real(kind=real128) :: x_barxi(3), dxdxi1_barxi(3), dxdxi2_barxi(3), dxdxi3_barxi(3)
    real(kind=real128) :: a(3), b(3,3), K(3,3), f(3), InvK(3,3), DetK
    ! Shape functions and its derivatives
#   define etype gtype
#   define xi barxi
#   define delta 0.d0
#   include <phi_and_dphidxik_2d.rc>
#   undef etype
#   undef xi
#   undef delta
    dphidxi3=0.d0
    varphi=phi*0.5d0*barxi(3)*tv3_md
    dvarphidxi1=dphidxi1*0.5d0*barxi(3)*tv3_md
    dvarphidxi2=dphidxi2*0.5d0*barxi(3)*tv3_md
    dvarphidxi3=phi*0.5d0*tv3_md
    ! x, dxdxik
    x_barxi=0.d0
    dxdxi1_barxi=0.d0
    dxdxi2_barxi=0.d0
    dxdxi3_barxi=0.d0
    do i=1,fbem_n_nodes(gtype)
      x_barxi     =x_barxi     +phi(i)*x(:,i)     +varphi(i)*v3_md(:,i)
      dxdxi1_barxi=dxdxi1_barxi+dphidxi1(i)*x(:,i)+dvarphidxi1(i)*v3_md(:,i)
      dxdxi2_barxi=dxdxi2_barxi+dphidxi2(i)*x(:,i)+dvarphidxi2(i)*v3_md(:,i)
      dxdxi3_barxi=dxdxi3_barxi+dphidxi3(i)*x(:,i)+dvarphidxi3(i)*v3_md(:,i)
    end do
    ! Build LSE
    a=x_barxi-dxdxi1_barxi*barxi(1)-dxdxi2_barxi*barxi(2)-dxdxi3_barxi*barxi(3)-x_i
    b(:,1)=dxdxi1_barxi
    b(:,2)=dxdxi2_barxi
    b(:,3)=dxdxi3_barxi
    do i=1,3
      do j=i,3
        K(i,j)=dot_product(b(:,i),b(:,j))
        if (i.ne.j) K(j,i)=K(i,j)
      end do
      f(i)=-dot_product(a,b(:,i))
    end do
    call fbem_invert_3x3_matrix(K,InvK,DetK)
    new_barxi=matmul(InvK,f)
  end subroutine nearest_minimization_iteration_3d_shell

  recursive subroutine fbem_nearest_minimization_1d(rn,gtype,x,x_i,error,nmax,barxi,rmin,info)
    implicit none
    ! I/O
    integer           :: rn
    integer           :: gtype
    real(kind=real64) :: x(rn,fbem_n_nodes(gtype))
    real(kind=real64) :: x_i(rn)
    real(kind=real64) :: error
    integer           :: nmax
    real(kind=real64) :: barxi !! Warning: it should be initialized
    real(kind=real64) :: rmin
    integer info
    ! Local
    real(kind=real64)  :: barxi_history(nmax+1)
    real(kind=real64)  :: error_history(nmax)
    integer            :: k
    real(kind=real128) :: r(rn)
    real(kind=real128) :: phi(fbem_n_nodes(gtype))
    real(kind=real128) :: aux(10)
    ! Initialize barxi
    if ((barxi.lt.-1.d0).or.(barxi.gt. 1.d0)) then
      barxi_history(1)=0.d0
    else
      barxi_history(1)=barxi
    end if
    ! Iterate
    k=1
    info=0
    do while (info.eq.0)
      call nearest_minimization_iteration_1d(rn,gtype,x,x_i,barxi_history(k),barxi_history(k+1))
      error_history(k)=abs(barxi_history(k+1)-barxi_history(k))*0.5d0
      if (error_history(k).le.error) then
        barxi=barxi_history(k+1)
        info=1
      else
        if (k.eq.nmax) then
          barxi=barxi_history(k+1)
          info=2
        else
          if (k.gt.5) then
            if (error_history(k).gt.error_history(k-2)) then
              barxi=barxi_history(k+1)
              info=3
            else
              k=k+1
            end if
          else
            k=k+1
          end if
        end if
      end if
    end do
    if (info.eq.1) then
      ! If the nearest is outside the element, use the vertices to check which is nearest.
      if ((barxi.lt.-1.d0).or.(barxi.gt.1.d0)) then
        barxi=-1.d0
        r=x(:,1)-x_i
        rmin=sqrt(dot_product(r,r))
        r=x(:,2)-x_i
        if (sqrt(dot_product(r,r)).lt.rmin) barxi=1.d0
      endif
      ! Calculate rmin
#     define etype gtype
#     define xi barxi
#     define delta 0.d0
#     include <phi_1d.rc>
#     undef etype
#     undef xi
#     undef delta
      r=0.d0
      do k=1,fbem_n_nodes(gtype)
        r=r+phi(k)*x(:,k)
      end do
      r=r-x_i
      rmin=sqrt(dot_product(r,r))
    end if
  end subroutine fbem_nearest_minimization_1d

  recursive subroutine fbem_nearest_minimization_2d(rn,gtype,x,x_i,error,nmax,barxi,rmin,info)
    implicit none
    ! I/O
    integer           :: rn
    integer           :: gtype
    real(kind=real64) :: x(rn,fbem_n_nodes(gtype))
    real(kind=real64) :: x_i(rn)
    real(kind=real64) :: error
    integer           :: nmax
    real(kind=real64) :: barxi(2) !! Warning: it should be initialized
    real(kind=real64) :: rmin
    integer           :: info
    ! Local
    real(kind=real64)               :: barxi_history(2,nmax+1)
    real(kind=real64)               :: error_history(nmax)
    integer                         :: k
    real(kind=real64)               :: xi(2,fbem_n_nodes(gtype))
    integer                         :: kedge, edge_gtype, methode
    real(kind=real64), allocatable  :: x_edge(:,:), rmin_edge(:), barxi_edge(:,:)
    real(kind=real64)               :: d, cl
    real(kind=real128), allocatable :: barx_edge(:,:)
    real(kind=real128)              :: r(rn)
    real(kind=real128)              :: phi(fbem_n_nodes(gtype))
    real(kind=real128)              :: aux(10)
    ! If barxi is not initially in the element, assign to it its center.
    if (.not.fbem_check_xi1xi2(gtype,barxi)) then
      select case (fbem_n_edges(gtype))
        case (3)
          barxi_history(:,1)=1.d0/3.d0
        case (4)
          barxi_history(:,1)=0.d0
      end select
    else
      barxi_history(:,1)=barxi
    end if
    ! Iterate
    k=1
    info=0
    do while (info.eq.0)
      call nearest_minimization_iteration_2d(rn,gtype,x,x_i,barxi_history(:,k),barxi_history(:,k+1))
      error_history(k)=sqrt((barxi_history(1,k+1)-barxi_history(1,k))**2+(barxi_history(2,k+1)-barxi_history(2,k))**2)*0.5d0
      if (error_history(k).le.error) then
        barxi=barxi_history(:,k+1)
        info=1
      else
        if (k.eq.nmax) then
          barxi=barxi_history(:,k+1)
          info=2
        else
          if (k.gt.5) then
            if (error_history(k).ge.error_history(k-2)) then
              barxi=barxi_history(:,k+1)
              info=3
            else
              k=k+1
            end if
          else
            k=k+1
          end if
        end if
      end if
    end do
    if (info.eq.1) then
      ! If the nearest is outside the element, use the edges to check which is nearest.
      if (.not.fbem_check_xi1xi2(gtype,barxi)) then
        allocate (rmin_edge(fbem_n_edges(gtype)),barx_edge(rn,fbem_n_edges(gtype)),barxi_edge(1,fbem_n_edges(gtype)))
        ! Check all edges
        do kedge=1,fbem_n_edges(gtype)
          edge_gtype=fbem_edge_type(kedge,gtype)
          allocate (x_edge(rn,fbem_n_nodes(edge_gtype)))
          do k=1,fbem_n_nodes(edge_gtype)
            x_edge(:,k)=x(:,fbem_edge_node(k,kedge,gtype))
          end do
          cl=fbem_characteristic_length(rn,edge_gtype,x_edge,1.d-12)
          call fbem_nearest_element_point_bem(rn,edge_gtype,x_edge,cl,x_i,barxi_edge(:,kedge),rmin_edge(kedge),d,methode)
          deallocate (x_edge)
        end do
        ! Choose the edge with the smaller rmin and calculate the nearest barxi in the element space
        kedge=minloc(rmin_edge,1)
        edge_gtype=fbem_edge_type(kedge,gtype)
        allocate (x_edge(2,fbem_n_nodes(edge_gtype)))
        ! xi coordinates of each node of the element
#       define etype gtype
#       define delta 0.d0
#       include <xi_2d.rc>
#       undef etype
#       undef delta
        do k=1,fbem_n_nodes(edge_gtype)
          x_edge(:,k)=xi(:,fbem_edge_node(k,kedge,gtype))
        end do
        ! xi coordinate (element space) of barxi_edge (edge space)
#       define etype edge_gtype
#       define xi barxi_edge(1,kedge)
#       define delta 0.0d0
#       include <phi_1d.rc>
#       undef etype
#       undef xi
#       undef delta
        barxi=0.d0
        do k=1,fbem_n_nodes(edge_gtype)
          barxi=barxi+phi(k)*x_edge(:,k)
        end do
        deallocate (x_edge,rmin_edge,barx_edge,barxi_edge)
      end if
      ! Calculate rmin
#     define etype gtype
#     define xi barxi
#     define delta 0.0d0
#     include <phi_2d.rc>
#     undef etype
#     undef xi
#     undef delta
      r=0.d0
      do k=1,fbem_n_nodes(gtype)
        r=r+phi(k)*x(:,k)
      end do
      r=r-x_i
      rmin=sqrt(dot_product(r,r))
    end if
  end subroutine fbem_nearest_minimization_2d

  subroutine fbem_nearest_minimization_3d_shell(gtype,x,v_md,tv_md,x_i,error,nmax,barxi,info)
    implicit none
    ! I/O
    integer           :: gtype
    real(kind=real64) :: x(3,fbem_n_nodes(gtype))
    real(kind=real64) :: v_md(3,3,fbem_n_nodes(gtype))
    real(kind=real64) :: tv_md(3,fbem_n_nodes(gtype))
    real(kind=real64) :: x_i(3)
    real(kind=real64) :: error
    integer           :: nmax
    real(kind=real64) :: barxi(3) !! Warning: it should be initialized
    integer           :: info ! 1: converged, 2: n>nmax, 3: does not converge
    ! Local
    real(kind=real64) :: v3_md(3,fbem_n_nodes(gtype))
    real(kind=real64) :: tv3_md(fbem_n_nodes(gtype))
    real(kind=real64) :: barxi_history(3,nmax+1)
    real(kind=real64) :: error_history(nmax)
    integer           :: k
    v3_md=v_md(:,3,:)
    tv3_md=tv_md(3,:)
    ! If barxi is not initially in the element, assign to it its center.
    if ((.not.fbem_check_xi1xi2(gtype,barxi(1:2))).or.(.not.fbem_check_xi(barxi(3)))) then
      select case (fbem_n_edges(gtype))
        case (3)
          barxi_history(1:2,1)=1.d0/3.d0
        case (4)
          barxi_history(1:2,1)=0.d0
      end select
      barxi_history(3,1)=0.d0
    else
      barxi_history(:,1)=barxi
    end if
    ! Iterate
    k=1
    info=0
    do while (info.eq.0)
      call nearest_minimization_iteration_3d_shell(gtype,x,v3_md,tv3_md,x_i,barxi_history(:,k),barxi_history(:,k+1))
      error_history(k)=sqrt(dot_product(barxi_history(:,k+1)-barxi_history(:,k),barxi_history(:,k+1)-barxi_history(:,k)))*0.5d0
      if (error_history(k).le.error) then
        barxi=barxi_history(:,k+1)
        info=1
      else
        if (k.eq.nmax) then
          barxi=barxi_history(:,k+1)
          info=2
        else
          if (k.gt.5) then
            if (error_history(k).ge.error_history(k-2)) then
              barxi=barxi_history(:,k+1)
              info=3
            else
              k=k+1
            end if
          else
            k=k+1
          end if
        end if
      end if
    end do
  end subroutine fbem_nearest_minimization_3d_shell

  !! Calculate the local coordinates of a point x belonging to an element.
  subroutine fbem_local_coordinates(rn,gtype,x,cl,x_i,xi_i,d)
    implicit none
    ! I/O
    integer           :: rn                            !! Spatial dimension
    integer           :: gtype                         !! Geometrical interpolation type
    real(kind=real64) :: x(rn,fbem_n_nodes(gtype))     !! Coordinates of the nodes
    real(kind=real64) :: x_i(rn)                       !! Collocation point
    real(kind=real64) :: cl                            !! Characteristic length
    real(kind=real64) :: xi_i(fbem_n_dimension(gtype)) !! Element reference coordinates of the point to x_i (Warning: it should be initialized)
    real(kind=real64) :: d                             !! Dimensionless distance d=rmin/cl
    ! Local
    real(kind=real64) :: error                         ! Absolute error for barxi calculation
    integer           :: nmax, info                    ! Maximum number of iterations and information
    real(kind=real64) :: rmin                          ! Distance between x(xi_i) and x_i
    ! Calculate
    select case (fbem_n_dimension(gtype))
      case (1)
        error=1.d-14
        nmax=25
        call fbem_nearest_minimization_1d(rn,gtype,x,x_i,error,nmax,xi_i(1),rmin,info)
      case (2)
        error=1.d-14
        nmax=25
        call fbem_nearest_minimization_2d(rn,gtype,x,x_i,error,nmax,xi_i,rmin,info)
      case (3)
        stop 'fbem_local_coordinates: n_dimension=3 not implemented yet'
      case default
        stop 'fbem_local_coordinates: n_dimension not valid'
    end select
    d=rmin/cl
    if (info.ne.1) then
      stop 'fbem_local_coordinates: fbem_local_coordinates has not converged'
    end if
  end subroutine fbem_local_coordinates

  ! ============================
  ! ADAPTATIVE ALGORITHM FOR BEM
  ! ============================

  !! Use the element nodes to compute a first approximation to barxi and rmin.
  !! Calculate d=rmin/cl
  !! If d<=1
  !!   Use the minimization algorithm to compute a numerically exact value of barxi and rmin.
  !!   If minimization algorithm does not converge
  !!     Use the sampling algorithm
  !!   End if
  !! End if
  recursive subroutine fbem_nearest_element_point_bem(rn,gtype,x,cl,x_i,barxi,rmin,d,method)
    implicit none
    ! I/O
    integer           :: rn                             !! Spatial dimension
    integer           :: gtype                          !! Geometrical interpolation type
    real(kind=real64) :: x(rn,fbem_n_nodes(gtype))      !! Coordinates of the nodes
    real(kind=real64) :: x_i(rn)                        !! Collocation point
    real(kind=real64) :: cl                             !! Characteristic length
    real(kind=real64) :: barxi(fbem_n_dimension(gtype)) !! Element reference coordinates of the element nearest point to x_i
    real(kind=real64) :: rmin                           !! Distance between x(barxi) and x_i
    real(kind=real64) :: d                              !! Dimensionless distance d=rmin/cl
    integer           :: method                         !! Method that has been used: 1 only nodes, 2 minimization, 3 sampling.
    ! Local
    real(kind=real64) :: error                          ! Absolute error for barxi calculation
    integer           :: nmax, info                     ! Maximum number of iterations and information
    ! Calculate
    call fbem_nearest_xi_nodes(rn,gtype,x,x_i,barxi,rmin)
    d=rmin/cl
    if (d.lt.1.d0) then
      select case (fbem_n_dimension(gtype))
        case (1)
          error=1.d-14
          nmax=20
          call fbem_nearest_minimization_1d(rn,gtype,x,x_i,error,nmax,barxi(1),rmin,info)
          if (info.ne.1) then
            call fbem_nearest_xi(rn,gtype,x,x_i,25,14,2,barxi(1),rmin)
            method=3
          else
            method=2
          end if
        case (2)
          error=1.d-14
          nmax=20
          call fbem_nearest_minimization_2d(rn,gtype,x,x_i,error,nmax,barxi,rmin,info)
          if (info.ne.1) then
            call fbem_nearest_xi1xi2(rn,gtype,x,x_i,25,14,2,barxi,rmin)
            method=3
          else
            method=2
          end if
        case default
          stop 'fbem_nearest_element_point_bem: n_dimension not valid'
      end select
      d=rmin/cl
    else
      method=1
    end if
  end subroutine fbem_nearest_element_point_bem

  ! ==============
  ! Transformation
  ! ==============

  !! This subroutine moves each nodal position to the average position of all its common nodes and itself, and then all nodes
  !! belonging to a symmetry plane are positioned exactly at the symmetry plane.
  subroutine fbem_transformation_collapse_nodal_positions(rn,n_nodes,node,n_symplanes,symplane_eid)
    implicit none
    ! I/O
    integer           :: rn
    integer           :: n_nodes
    type(fbem_node)   :: node(n_nodes)
    integer           :: n_symplanes
    integer           :: symplane_eid(3)
    ! Local
    integer           :: kni, knj
    real(kind=real64) :: x(rn)
    logical           :: processed(n_nodes)
    processed=.false.
    do kni=1,n_nodes
      if ((node(kni)%n_nodes.gt.0).and.(processed(kni).eqv.(.false.))) then
        x=node(kni)%x
        do knj=1,node(kni)%n_nodes
          x=x+node(node(kni)%node(knj))%x
        end do
        x=x/dble(1+node(kni)%n_nodes)
        node(kni)%x=x
        processed(kni)=.true.
        do knj=1,node(kni)%n_nodes
          node(node(kni)%node(knj))%x=x
          processed(node(kni)%node(knj))=.true.
        end do
      end if
    end do
    if (n_symplanes.gt.0) then
      do kni=1,n_nodes
        if (node(kni)%n_symplanes.gt.0) then
          do knj=1,node(kni)%n_symplanes
            node(kni)%x(symplane_eid(node(kni)%symplane(knj)))=0.d0
          end do
        end if
      end do
    end if
  end subroutine fbem_transformation_collapse_nodal_positions

  !! Function that moves the given 1D reference coordinate to the incenter if it is in an vertex
  function fbem_move_xi_from_vertex(xi,delta)
    implicit none
    !! xi moved
    real(kind=real64) :: fbem_move_xi_from_vertex
    !! xi coordinate
    real(kind=real64) :: xi
    !! Displacement from the xi coordinate to the element centroid.
    real(kind=real64) :: delta
    ! Only displace the coordinate if it is in an vertex
    if (fbem_check_xi_vertex(xi).eqv.(.true.)) then
      fbem_move_xi_from_vertex=xi*(1.0d0-delta)
    else
      fbem_move_xi_from_vertex=xi
    end if
  end function fbem_move_xi_from_vertex

  !! Function that moves the given 2D reference coordinate to the incenter if it is in an edge
  function fbem_move_xi1xi2_from_edge(etype,xi,delta)
    implicit none
    !! xi1,xi2 moved
    real(kind=real64) :: fbem_move_xi1xi2_from_edge(2)
    !! Element interpolation type: <tt>etype={fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    integer :: etype
    !! xi1,xi2 coordinate
    real(kind=real64) :: xi(2)
    !! Displacement from the xi1,xi2 coordinate to the element incenter. Displacement is measured as a per unit radius with respect
    !! to the incenter and vertex.
    real(kind=real64) :: delta
    ! Only displace the coordinate if it is in an edge
    if (fbem_check_xi1xi2_edge(etype,xi).eqv.(.true.)) then
      select case (etype)
        case (fbem_tri3,fbem_tri6)
          fbem_move_xi1xi2_from_edge(1)=xi(1)*(1.0d0-delta)+0.333333333333333333d0*delta
          fbem_move_xi1xi2_from_edge(2)=xi(2)*(1.0d0-delta)+0.333333333333333333d0*delta
        case (fbem_quad4,fbem_quad8,fbem_quad9)
          fbem_move_xi1xi2_from_edge(1)=xi(1)*(1.0d0-delta)
          fbem_move_xi1xi2_from_edge(2)=xi(2)*(1.0d0-delta)
        case default
          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={tri3,tri6,quad4,quad8,quad9}')
      end select
    else
      fbem_move_xi1xi2_from_edge(1)=xi(1)
      fbem_move_xi1xi2_from_edge(2)=xi(2)
    end if
  end function fbem_move_xi1xi2_from_edge

  !! Function that moves the given 3D reference coordinate to the incenter
  function fbem_move_xi1xi2xi3_from_face(etype,xi,delta)
    implicit none
    !! xi1,xi2 moved
    real(kind=real64) :: fbem_move_xi1xi2xi3_from_face(3)
    !! Element interpolation type: <tt>etype={fbem_tet4,fbem_tet10,fbem_hex8,fbem_hex20,fbem_hex27}</tt>.
    integer :: etype
    !! xi1,xi2,xi3 coordinate
    real(kind=real64) :: xi(3)
    !! Displacement from the xi1,xi2,xi3 coordinate to the element incenter. Displacement is measured as a per unit radius with respect
    !! to the incenter and vertex.
    real(kind=real64) :: delta
    ! Movimiento de todo hacia el incentro
    select case (etype)
      case (fbem_tet4,fbem_tet10)
        fbem_move_xi1xi2xi3_from_face=xi*(1.0d0-delta)+0.25d0*delta
      case (fbem_hex8,fbem_hex20,fbem_hex27)
        fbem_move_xi1xi2xi3_from_face=xi*(1.0d0-delta)
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={fbem_tet4,fbem_tet10,fbem_hex8,fbem_hex20,fbem_hex27}')
    end select
  end function fbem_move_xi1xi2xi3_from_face

end module fbem_geometry
