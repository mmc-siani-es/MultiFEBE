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

! For safety reasons, is possible to compile the code with a xi or xi1,xi2 domain checking, e.g., in every function call
! the xi or xi1,xi2 are checked to see if they are in the domain. If CHECK_XIS == 1, checking is enabled, and a warning message
! is displayed. If CHECK_XIS == 2, checking is enabled, an error message is displayed and execution stops. Other values disable
! any checking.
#define CHECK_XIS 0

!! @author Jacob David Rodriguez Bordon (jacobdavid.rodriguezbordon@ulpgc.es)
!!
!! @version 2.0
!!
!! <b> This module implements parameters and functions related to typical elements used by FEM and BEM codes. </b>
!!
!! <h2>DESCRIPTION</h2>
!!
!! The implemented elements are:
!! -0D element:
!!   -<tt>point1</tt>: 1 node Lagrange constant element.
!!
!! -1D elements:
!!   -<tt>line2</tt>: 2 nodes Lagrange linear element.
!!>
!!       |
!! 1 --------- 2 ---> xi
!!       |
!!
!!     line2
!!<
!!   -<tt>line3</tt>: 3 nodes Lagrange quadratic element.
!!>
!!       |
!! 1 --- 3 --- 2 ---> xi
!!       |
!!
!!     line3
!!<
!!   -<tt>line4</tt>: 4 nodes Lagrange cubic element.
!!>
!!       |
!! 1 - 3 - 4 - 2 ---> xi
!!       |
!!
!!     line4
!!<
!! -2D elements:
!!   -<tt>tri3</tt>: Triangular 3 nodes Lagrange linear element.
!!>
!!         xi2
!!         ^
!!        /
!!       2
!!      / \
!!     /   \
!!    /     \
!!   /       \
!!  /         \
!! 3 --------- 1 --- > xi1
!!
!!     tri3
!!<
!!   -<tt>tri6</tt>: Triangular 6 nodes Lagrange quadratic element.
!!>
!!         xi2
!!         ^
!!        /
!!       2
!!      / \
!!     /   \
!!    5     4
!!   /       \
!!  /         \
!! 3 --- 6 --- 1 --- > xi1
!!
!!     tri6
!!<
!!   -<tt>quad4</tt>: Quadrilateral 4 nodes Lagrange bi-linear element.
!!>
!!      xi2
!!       ^
!!       |
!!       |
!! 4 --------- 3
!! |     |     |
!! |     |     |
!! |     |-----|---- > xi1
!! |           |
!! |           |
!! 1 --------- 2
!!
!!     quad4
!!<
!!   -<tt>quad8</tt>: Quadrilateral 8 nodes serendipity quadratic element.
!!>
!!      xi2
!!       ^
!!       |
!!       |
!! 4 --- 7 --- 3
!! |     |     |
!! |     |     |
!! 8     |---- 6 --- > xi1
!! |           |
!! |           |
!! 1 --- 5 --- 2
!!
!!     quad8
!!<
!!   -<tt>quad9</tt>: Quadrilateral 9 nodes Lagrange bi-quadratic element.
!!>
!!      xi2
!!       ^
!!       |
!!       |
!! 4 --- 7 --- 3
!! |     |     |
!! |     |     |
!! 8     9 --- 6 --- > xi1
!! |           |
!! |           |
!! 1 --- 5 --- 2
!!
!!     quad9
!!<
!!
!! Element shape functions and derivatives have a parameter named <tt>delta</tt>. It is the distance from edges that nodes
!! are moved to the element inside, i.e. if <tt>delta>0.0</tt> the element is discontinuous.
!!
!! Module uses resource files that contains reference coordinates and calculations of shape functions and its derivatives. It
!! lets the module easily maintainable, and programmer can use the inline code instead of functions, which save computational
!! cost (function call costs). For example, use:
!!
!!>
!!#define dphidxi dphidxi_v
!!#include <dphidxi_1d.rc>
!!#undef dphidxi
!!<
!!
!! instead of:
!!
!!>
!!dphidxi_v=fbem_dphidxi(etype,xi)
!!<
!!
!! Note: it is necessary to have several variables declared, check resource files. The variable name used in resource files can be
!! changed by the pre-processor by using <tt>#define</tt> directive.
!!
!! <b>DATE:</b>
!! January, 2013
!!
!! <b>BUGS:</b>
!! Not found yet
!!
!! <b>TODO:</b>
!! -Implement 3D elements.
!! -Modify implementation of interfaces and functions to integrate 2D and 3D elements in the same functions.
!! -Improve documentation.

module fbem_shape_functions

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! needed fbem module
  use fbem_numerical
  use fbem_string_handling
  use fbem_quad_rules

  ! No implicit variables are allowed
  implicit none

  ! By default all are private
  private

  ! Initialization
  public :: fbem_init_shape_functions_module

  ! Standard shape functions
  ! Internal identifier parameters for implementing the elements. Functions should use only the parameter name and not the numbers.
  integer, parameter, public :: fbem_point1=1  !! Internal identifier for <tt>point1</tt> Lagrange element (constant).
  integer, parameter, public :: fbem_line2 =2  !! Internal identifier for <tt>line2</tt> Lagrange element (linear).
  integer, parameter, public :: fbem_line3 =3  !! Internal identifier for <tt>line3</tt> Lagrange element (quadratic).
  integer, parameter, public :: fbem_line4 =4  !! Internal identifier for <tt>line4</tt> Lagrange element (cubic).
  integer, parameter, public :: fbem_tri3  =5  !! Internal identifier for <tt>tri3</tt> Lagrange element (linear).
  integer, parameter, public :: fbem_tri6  =6  !! Internal identifier for <tt>tri6</tt> Lagrange element (quadratic).
  integer, parameter, public :: fbem_quad4 =7  !! Internal identifier for <tt>quad4</tt> Lagrange element (bilinear).
  integer, parameter, public :: fbem_quad8 =8  !! Internal identifier for <tt>quad8</tt> serendipity element (biquadratic incomplete).
  integer, parameter, public :: fbem_quad9 =9  !! Internal identifier for <tt>quad9</tt> Lagrange element (biquadratic complete).
  integer, parameter, public :: fbem_tet4  =10 !! Internal identifier for <tt>tet4</tt> Lagrange element (linear).
  integer, parameter, public :: fbem_tet10 =11 !! Internal identifier for <tt>tet10</tt> Lagrange element (quadratic).
  integer, parameter, public :: fbem_hex8  =12 !! Internal identifier for <tt>hex8</tt> Lagrange element (trilinear).
  integer, parameter, public :: fbem_hex20 =13 !! Internal identifier for <tt>hex20</tt> serendipity element (triquadratic incomplete).
  integer, parameter, public :: fbem_hex27 =14 !! Internal identifier for <tt>hex27</tt> Lagrange element (triquadratic).
  ! Special shape functions
  integer, parameter, public :: fbem_line2point1=15 !! Internal identifier for <tt>line2point1</tt> with line2 interpolation in xi1 and constant in xi2.
  integer, parameter, public :: fbem_line3point1=16 !! Internal identifier for <tt>line3point1</tt> with line3 interpolation in xi1 and constant in xi2.
  integer, parameter, public :: fbem_line3_qp1t =17 !! Internal identifier for <tt>line3</tt> for traction quarter-point elements (F. Chirino et al. (1989)).
  integer, parameter, public :: fbem_line3_qp2t =18 !! Internal identifier for <tt>line3</tt> for traction quarter-point elements (F. Chirino et al. (1989)).
  integer, parameter, public :: fbem_line3_mqp1u=19 !! Internal identifier for <tt>line3</tt> modified for displacement quarter-point elements (L.J. Gray et al. (2003)), with the crack-tip at node 1 (L.J. Gray et al. Improved quarter-point crack tip element (2003)).
  integer, parameter, public :: fbem_line3_mqp2u=20 !! Internal identifier for <tt>line3</tt> modified for displacement quarter-point elements (L.J. Gray et al. (2003)), with the crack-tip at node 2 (L.J. Gray et al. Improved quarter-point crack tip element (2003)).
  integer, parameter, public :: fbem_line3_mqp1t=21 !! Internal identifier for <tt>line3</tt> modified for traction singular-quarter-point elements (F. Chirino et al. (1998), L.J. Gray et al. (2003)), with the crack-tip at node 1 (L.J. Gray et al. Improved quarter-point crack tip element (2003)).
  integer, parameter, public :: fbem_line3_mqp2t=22 !! Internal identifier for <tt>line3</tt> modified for traction singular-quarter-point elements (F. Chirino et al. (1998), L.J. Gray et al. (2003)), with the crack-tip at node 2 (L.J. Gray et al. Improved quarter-point crack tip element (2003)).
  ! Some new special purpose
  integer, parameter, public :: fbem_tri4       =23 !! Internal identifier for <tt>tri4</tt> (linear triangle with cubic bubble node)
  integer, parameter, public :: fbem_tri1       =24 !! Internal identifier for <tt>tri1</tt> (constant triangle)
  integer, parameter, public :: fbem_quad1      =25 !! Internal identifier for <tt>quad1</tt> (constant quadrilateral)
  integer, parameter, public :: fbem_line1      =26 !! Internal identifier for <tt>line1</tt> (constant line)

  ! Parameters for each element type (only for the standard Lagrange elements)

  !! Element reference Spatial dimensions.
  !!
  !! Example: <tt>fbem_n_dimension(fbem_tri3)</tt> returns <tt>2</tt>, which means that <tt>fbem_tri3</tt> is a 2D element.
  !!
  !! It is useful when implementing functions/subroutines.
  integer, allocatable, public :: fbem_n_dimension(:)

  !! Element order
  integer, allocatable, public :: fbem_element_order(:)

  !! Element number of nodes.
  !!
  !! Example: <tt>fbem_n_nodes(fbem_quad9)</tt> returns <tt>9</tt>, which means that <tt>fbem_quad9</tt> has 9 nodes.
  !!
  !! It is useful when implementing functions/subroutines.
  integer, allocatable, public :: fbem_n_nodes(:)

  !! Element number of vertices.
  !!
  !! Example: <tt>fbem_n_vertices(fbem_quad4)</tt> returns <tt>4</tt>, which means that <tt>fbem_quad4</tt> is a quadrilateral
  !! element.
  !!
  !! It is useful when implementing functions/subroutines.
  integer, allocatable, public :: fbem_n_vertices(:)

  !! Element inversion vector (converts local numeration)
  integer, allocatable, public :: fbem_element_invert(:,:)

  !! Type of node location in each element:
  integer, parameter, public   :: fbem_loctype_vertex         =1
  integer, parameter, public   :: fbem_loctype_edge_interior  =2
  integer, parameter, public   :: fbem_loctype_face_interior  =3
  integer, parameter, public   :: fbem_loctype_volume_interior=4
  integer, allocatable, public :: fbem_node_loctype(:,:)

  !! Number of (sub-)edges of each type of element.
  integer, allocatable, public :: fbem_n_edges(:)
  !! Type of (sub-)edges of each type of element.
  integer, allocatable, public :: fbem_edge_type(:,:)
  !! Nodes of (sub-)edges of each type of element.
  integer, allocatable, public :: fbem_edge_node(:,:,:)

  !! Number of (sub-)faces of each type of element.
  integer, allocatable, public :: fbem_n_faces(:)
  !! Type of (sub-)faces of each type of element.
  integer, allocatable, public :: fbem_face_type(:,:)
  !! Nodes of (sub-)faces of each type of element.
  integer, allocatable, public :: fbem_face_node(:,:,:)

  !! Domain tolerance used by <tt>fbem_check_xi</tt> and <tt>fbem_check_xi1xi2</tt> functions.
  real(kind=real64), parameter, public :: check_xi_tolerance=0.5d-12 ! 1.d-6

  ! Reference coordinates
  public :: fbem_xi
  public :: fbem_xi1xi2
  public :: fbem_xi_hybrid

  ! Coordinates check
  public :: fbem_check_xi
  public :: fbem_check_xi_warning
  public :: fbem_check_xi_error
  public :: fbem_check_xi_vertex
  public :: fbem_check_xi1xi2
  public :: fbem_check_xi1xi2_warning
  public :: fbem_check_xi1xi2_error
  public :: fbem_check_xi1xi2_edge
  public :: fbem_check_xi_on_element_boundary

  ! Shape functions
  public :: fbem_phi
  public :: fbem_phi_hybrid

  ! Shape functions first derivatives
  public :: fbem_dphidxi
  public :: fbem_dphidxi1
  public :: fbem_dphidxi2
  public :: fbem_dphidx

  ! Tangential derivative
  public :: fbem_psi

  ! Shape function integrals
  public :: fbem_shape_functions_integrals
  public :: fbem_shape_functions_x_integrals

  !! This interface switchs between functions that evaluates shape functions for 1D elements and 2D/3D elements.
  interface fbem_phi
    module procedure fbem_phi_1d
    module procedure fbem_phi_nd
  end interface fbem_phi

contains

  ! Shape function derivatives with respect to global coordinates for the calculation of tangential or volume gradients of a field
  ! variable interpolated over an element.
  function fbem_psi(rn,gtype,x,ftype,ftype_delta,xi)
    implicit none
    ! I/O variables
    integer           :: rn
    integer           :: gtype
    real(kind=real64) :: x(rn,fbem_n_nodes(gtype))
    integer           :: ftype
    real(kind=real64) :: ftype_delta
    real(kind=real64) :: xi(fbem_n_dimension(gtype))
    real(kind=real64) :: fbem_psi(fbem_n_nodes(ftype),rn)
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
      case (0)
        fbem_psi=0
      ! Line gradient
      case (1)
        xi1d=xi(1)
        ! Geometrical first derivatives at xi
#       define etype gtype
#       define delta 0.0d0
#       define xi xi1d
#       include <dphidxi_1d.rc>
#       undef etype
#       undef delta
#       undef xi
        ! Calculate tangent (velocity) vector
        T=0.d0
        do kphi=1,fbem_n_nodes(gtype)
          T=T+dphidxi(kphi)*x(:,kphi)
        end do
        ! Geometric jacobian
        jg=sqrt(dot_product(T,T))
        ! Unit normal
        t=T/jg
        ! Functional first derivatives at xi
#       define etype ftype
#       define delta ftype_delta
#       define xi xi1d
#       include <dphidxi_1d.rc>
#       undef etype
#       undef delta
#       undef xi
        ! Build psi
        do kphi=1,fbem_n_nodes(ftype)
          fbem_psi(kphi,:)=dphidxi(kphi)/jg*t
        end do
      case (2)
        select case (rn)
          ! Surface gradient
          case (2)
            ! Geometrical first derivatives at xi
#           define etype gtype
#           define delta 0.0d0
#           include <dphidxik_2d.rc>
#           undef etype
#           undef delta
            ! Calculate jacobian matrix, its inverse and determinant
            J=0.d0
            do kphi=1,fbem_n_nodes(gtype)
              J(1,:)=J(1,:)+dphidxi1(kphi)*x(:,kphi)
              J(2,:)=J(2,:)+dphidxi2(kphi)*x(:,kphi)
            end do
            call fbem_invert_2x2_matrix(J,invJ,detJ)
            ! Functional first derivatives at xi
#           define etype ftype
#           define delta ftype_delta
#           include <dphidxik_2d.rc>
#           undef etype
#           undef delta
#           undef xi
            ! Build psi
            do kphi=1,fbem_n_nodes(ftype)
              fbem_psi(kphi,1)=invJ(1,1)*dphidxi1(kphi)+invJ(1,2)*dphidxi2(kphi)
              fbem_psi(kphi,2)=invJ(2,1)*dphidxi1(kphi)+invJ(2,2)*dphidxi2(kphi)
            end do
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
            ! Functional first derivatives at xi
#           define etype ftype
#           define delta ftype_delta
#           include <dphidxik_2d.rc>
#           undef etype
#           undef delta
#           undef xi
            ! The orthogonal system at xi
            ! v1 = t1
            v1=t1
            ! v2 = n x t1
            v2(1)=n(2)*t1(3)-n(3)*t1(2)
            v2(2)=n(3)*t1(1)-n(1)*t1(3)
            v2(3)=n(1)*t1(2)-n(2)*t1(1)
            ! Orthogonal coordinates transformation matrix
            Mvt(1,1)=1.d0
            Mvt(1,2)=dot_product(v1,t2)
            Mvt(2,1)=0.d0
            Mvt(2,2)=dot_product(v2,t2)
            ! xi derivatives with respect to surface orthogonal coordinates
            detMvt=Mvt(1,1)*Mvt(2,2)-Mvt(1,2)*Mvt(2,1)
            dxidv(1,1)= Mvt(2,2)/(jxi1*detMvt)
            dxidv(1,2)=-Mvt(1,2)/(jxi1*detMvt)
            dxidv(2,1)=-Mvt(2,1)/(jxi2*detMvt)
            dxidv(2,2)= Mvt(1,1)/(jxi2*detMvt)
            ! Build psi
            do kphi=1,fbem_n_nodes(ftype)
              fbem_psi(kphi,:)=(dphidxi1(kphi)*dxidv(1,1)+dphidxi2(kphi)*dxidv(2,1))*v1+&
                               (dphidxi1(kphi)*dxidv(1,2)+dphidxi2(kphi)*dxidv(2,2))*v2
            end do
        end select
      ! Volume gradient
      case (3)
        if (ftype_delta.ne.0.d0) stop '3D elements only with delta=0'
        ! Geometrical first derivatives at xi
#       define etype gtype
#       define delta 0.0d0
#       include <dphidxik_3d.rc>
#       undef etype
#       undef delta
        ! Calculate jacobian matrix, its inverse and determinant
        J=0.d0
        do kphi=1,fbem_n_nodes(gtype)
          J(1,:)=J(1,:)+dphidxi1(kphi)*x(:,kphi)
          J(2,:)=J(2,:)+dphidxi2(kphi)*x(:,kphi)
          J(3,:)=J(3,:)+dphidxi3(kphi)*x(:,kphi)
        end do
        call fbem_invert_3x3_matrix(J,invJ,detJ)
        ! Functional first derivatives at xi
#       define etype ftype
#       define delta ftype_delta
#       include <dphidxik_3d.rc>
#       undef etype
#       undef delta
#       undef xi
        ! Build psi
        do kphi=1,fbem_n_nodes(ftype)
          fbem_psi(kphi,1)=invJ(1,1)*dphidxi1(kphi)+invJ(1,2)*dphidxi2(kphi)+invJ(1,3)*dphidxi3(kphi)
          fbem_psi(kphi,2)=invJ(2,1)*dphidxi1(kphi)+invJ(2,2)*dphidxi2(kphi)+invJ(2,3)*dphidxi3(kphi)
          fbem_psi(kphi,3)=invJ(3,1)*dphidxi1(kphi)+invJ(3,2)*dphidxi2(kphi)+invJ(3,3)*dphidxi3(kphi)
        end do
    end select
  end function fbem_psi


  !! Initialization of this module
  subroutine fbem_init_shape_functions_module
    implicit none
    integer :: maxtype, maxnodes, maxedges, maxnodesedge
    ! Change by hand this parameter
    maxtype=26
    if (.not.allocated(fbem_n_dimension)) then
      allocate (fbem_n_dimension(maxtype))
      fbem_n_dimension=0
      fbem_n_dimension(fbem_point1     )=0
      fbem_n_dimension(fbem_line1      )=1
      fbem_n_dimension(fbem_line2      )=1
      fbem_n_dimension(fbem_line3      )=1
      fbem_n_dimension(fbem_line4      )=1
      fbem_n_dimension(fbem_tri3       )=2
      fbem_n_dimension(fbem_tri6       )=2
      fbem_n_dimension(fbem_quad4      )=2
      fbem_n_dimension(fbem_quad8      )=2
      fbem_n_dimension(fbem_quad9      )=2
      fbem_n_dimension(fbem_tet4       )=3
      fbem_n_dimension(fbem_tet10      )=3
      fbem_n_dimension(fbem_hex8       )=3
      fbem_n_dimension(fbem_hex20      )=3
      fbem_n_dimension(fbem_hex27      )=3
      fbem_n_dimension(fbem_line2point1)=2
      fbem_n_dimension(fbem_line3point1)=2
      fbem_n_dimension(fbem_line3_qp1t )=1
      fbem_n_dimension(fbem_line3_qp2t )=1
      fbem_n_dimension(fbem_line3_mqp1u)=1
      fbem_n_dimension(fbem_line3_mqp2u)=1
      fbem_n_dimension(fbem_line3_mqp1t)=1
      fbem_n_dimension(fbem_line3_mqp2t)=1
      fbem_n_dimension(fbem_tri4       )=2
      fbem_n_dimension(fbem_tri1       )=2
      fbem_n_dimension(fbem_quad1      )=2
    end if
    if (.not.allocated(fbem_element_order)) then
      allocate (fbem_element_order(maxtype))
      fbem_element_order=0
      fbem_element_order(fbem_point1     )=0
      fbem_element_order(fbem_line1      )=0
      fbem_element_order(fbem_line2      )=1
      fbem_element_order(fbem_line3      )=2
      fbem_element_order(fbem_line4      )=3
      fbem_element_order(fbem_tri3       )=1
      fbem_element_order(fbem_tri6       )=2
      fbem_element_order(fbem_quad4      )=1
      fbem_element_order(fbem_quad8      )=2
      fbem_element_order(fbem_quad9      )=2
      fbem_element_order(fbem_tet4       )=1
      fbem_element_order(fbem_tet10      )=2
      fbem_element_order(fbem_hex8       )=1
      fbem_element_order(fbem_hex20      )=2
      fbem_element_order(fbem_hex27      )=2
      fbem_element_order(fbem_line2point1)=1
      fbem_element_order(fbem_line3point1)=2
      fbem_element_order(fbem_line3_qp1t )=2
      fbem_element_order(fbem_line3_qp2t )=2
      fbem_element_order(fbem_line3_mqp1u)=2
      fbem_element_order(fbem_line3_mqp2u)=2
      fbem_element_order(fbem_line3_mqp1t)=2
      fbem_element_order(fbem_line3_mqp2t)=2
      fbem_element_order(fbem_tri4       )=1
      fbem_element_order(fbem_tri1       )=0
      fbem_element_order(fbem_quad1      )=0
    end if
    if (.not.allocated(fbem_n_nodes)) then
      allocate (fbem_n_nodes(maxtype))
      fbem_n_nodes=0
      fbem_n_nodes(fbem_point1     )=1
      fbem_n_nodes(fbem_line1      )=1
      fbem_n_nodes(fbem_line2      )=2
      fbem_n_nodes(fbem_line3      )=3
      fbem_n_nodes(fbem_line4      )=4
      fbem_n_nodes(fbem_tri3       )=3
      fbem_n_nodes(fbem_tri6       )=6
      fbem_n_nodes(fbem_quad4      )=4
      fbem_n_nodes(fbem_quad8      )=8
      fbem_n_nodes(fbem_quad9      )=9
      fbem_n_nodes(fbem_tet4       )=4
      fbem_n_nodes(fbem_tet10      )=10
      fbem_n_nodes(fbem_hex8       )=8
      fbem_n_nodes(fbem_hex20      )=20
      fbem_n_nodes(fbem_hex27      )=27
      fbem_n_nodes(fbem_line2point1)=2
      fbem_n_nodes(fbem_line3point1)=3
      fbem_n_nodes(fbem_line3_qp1t )=3
      fbem_n_nodes(fbem_line3_qp2t )=3
      fbem_n_nodes(fbem_line3_mqp1u)=3
      fbem_n_nodes(fbem_line3_mqp2u)=3
      fbem_n_nodes(fbem_line3_mqp1t)=3
      fbem_n_nodes(fbem_line3_mqp2t)=3
      fbem_n_nodes(fbem_tri4       )=4
      fbem_n_nodes(fbem_tri1       )=1
      fbem_n_nodes(fbem_quad1      )=1
    end if
    maxnodes=maxval(fbem_n_nodes)
    if (.not.allocated(fbem_n_vertices)) then
      allocate (fbem_n_vertices(maxtype))
      fbem_n_vertices=0
      fbem_n_vertices(fbem_point1     )=0
      fbem_n_vertices(fbem_line1      )=0
      fbem_n_vertices(fbem_line2      )=2
      fbem_n_vertices(fbem_line3      )=2
      fbem_n_vertices(fbem_line4      )=2
      fbem_n_vertices(fbem_tri3       )=3
      fbem_n_vertices(fbem_tri6       )=3
      fbem_n_vertices(fbem_quad4      )=4
      fbem_n_vertices(fbem_quad8      )=4
      fbem_n_vertices(fbem_quad9      )=4
      fbem_n_vertices(fbem_tet4       )=4
      fbem_n_vertices(fbem_tet10      )=4
      fbem_n_vertices(fbem_hex8       )=8
      fbem_n_vertices(fbem_hex20      )=8
      fbem_n_vertices(fbem_hex27      )=8
!      fbem_n_vertices(fbem_line2point1)=2
!      fbem_n_vertices(fbem_line3point1)=2
!      fbem_n_vertices(fbem_line3_qp1t )=2
!      fbem_n_vertices(fbem_line3_qp2t )=2
!      fbem_n_vertices(fbem_line3_mqp1u)=2
!      fbem_n_vertices(fbem_line3_mqp2u)=2
!      fbem_n_vertices(fbem_line3_mqp1t)=2
!      fbem_n_vertices(fbem_line3_mqp2t)=2
      fbem_n_vertices(fbem_tri4       )=3
      fbem_n_vertices(fbem_tri1       )=0
      fbem_n_vertices(fbem_quad1      )=0
    end if
    if (.not.allocated(fbem_element_invert)) then
      allocate (fbem_element_invert(maxnodes,maxtype))
      fbem_element_invert=0
      ! point1
      fbem_element_invert(1,fbem_point1)=1
      ! line1
      fbem_element_invert(1,fbem_line1)=1
      ! line2
      fbem_element_invert(1,fbem_line2)=2
      fbem_element_invert(2,fbem_line2)=1
      ! line3
      fbem_element_invert(:,fbem_line3)=fbem_element_invert(:,fbem_line2)
      fbem_element_invert(3,fbem_line3)=3
      ! line4
      fbem_element_invert(:,fbem_line4)=fbem_element_invert(:,fbem_line2)
      fbem_element_invert(3,fbem_line4)=4
      fbem_element_invert(4,fbem_line4)=3
      ! tri3
      fbem_element_invert(1,fbem_tri3)=1
      fbem_element_invert(2,fbem_tri3)=3
      fbem_element_invert(3,fbem_tri3)=2
      ! tri6
      fbem_element_invert(:,fbem_tri6)=fbem_element_invert(:,fbem_tri3)
      fbem_element_invert(4,fbem_tri6)=6
      fbem_element_invert(5,fbem_tri6)=5
      fbem_element_invert(6,fbem_tri6)=4
      ! quad4
      fbem_element_invert(1,fbem_quad4)=1
      fbem_element_invert(2,fbem_quad4)=4
      fbem_element_invert(3,fbem_quad4)=3
      fbem_element_invert(4,fbem_quad4)=2
      ! quad8
      fbem_element_invert(:,fbem_quad8)=fbem_element_invert(:,fbem_quad4)
      fbem_element_invert(5,fbem_quad8)=8
      fbem_element_invert(6,fbem_quad8)=7
      fbem_element_invert(7,fbem_quad8)=6
      fbem_element_invert(8,fbem_quad8)=5
      ! quad9
      fbem_element_invert(:,fbem_quad9)=fbem_element_invert(:,fbem_quad8)
      fbem_element_invert(9,fbem_quad9)=9
      ! Invert volume elements?
!      ! Crack-tip elements derived from line3
!      fbem_element_invert(:,fbem_line3_qp1t )=fbem_element_invert(:,fbem_line3)
!      fbem_element_invert(:,fbem_line3_qp2t )=fbem_element_invert(:,fbem_line3)
!      fbem_element_invert(:,fbem_line3_mqp1u)=fbem_element_invert(:,fbem_line3)
!      fbem_element_invert(:,fbem_line3_mqp2u)=fbem_element_invert(:,fbem_line3)
!      fbem_element_invert(:,fbem_line3_mqp1t)=fbem_element_invert(:,fbem_line3)
!      fbem_element_invert(:,fbem_line3_mqp2t)=fbem_element_invert(:,fbem_line3)
      ! tri4
      fbem_element_invert(1,fbem_tri4)=1
      fbem_element_invert(2,fbem_tri4)=3
      fbem_element_invert(3,fbem_tri4)=2
      ! tri1
      fbem_element_invert(1,fbem_tri1)=1
      ! quad1
      fbem_element_invert(1,fbem_quad1)=1
    end if
    if (.not.allocated(fbem_node_loctype)) then
      allocate (fbem_node_loctype(maxnodes,maxtype))
      fbem_node_loctype=0
      ! point1
      fbem_node_loctype(1,fbem_point1)=fbem_loctype_vertex
      ! line1
      fbem_node_loctype(1,fbem_line1)=fbem_loctype_edge_interior
      ! line2
      fbem_node_loctype(1,fbem_line2)=fbem_loctype_vertex
      fbem_node_loctype(2,fbem_line2)=fbem_loctype_vertex
      ! line3
      fbem_node_loctype(:,fbem_line3)=fbem_node_loctype(:,fbem_line2)
      fbem_node_loctype(3,fbem_line3)=fbem_loctype_edge_interior
      ! line4
      fbem_node_loctype(:,fbem_line4)=fbem_node_loctype(:,fbem_line2)
      fbem_node_loctype(3,fbem_line4)=fbem_loctype_edge_interior
      fbem_node_loctype(4,fbem_line4)=fbem_loctype_edge_interior
      ! tri3
      fbem_node_loctype(1,fbem_tri3)=fbem_loctype_vertex
      fbem_node_loctype(2,fbem_tri3)=fbem_loctype_vertex
      fbem_node_loctype(3,fbem_tri3)=fbem_loctype_vertex
      ! tri6
      fbem_node_loctype(:,fbem_tri6)=fbem_node_loctype(:,fbem_tri3)
      fbem_node_loctype(4,fbem_tri6)=fbem_loctype_edge_interior
      fbem_node_loctype(5,fbem_tri6)=fbem_loctype_edge_interior
      fbem_node_loctype(6,fbem_tri6)=fbem_loctype_edge_interior
      ! quad4
      fbem_node_loctype(1,fbem_quad4)=fbem_loctype_vertex
      fbem_node_loctype(2,fbem_quad4)=fbem_loctype_vertex
      fbem_node_loctype(3,fbem_quad4)=fbem_loctype_vertex
      fbem_node_loctype(4,fbem_quad4)=fbem_loctype_vertex
      ! quad8
      fbem_node_loctype(:,fbem_quad8)=fbem_node_loctype(:,fbem_quad4)
      fbem_node_loctype(5,fbem_quad8)=fbem_loctype_edge_interior
      fbem_node_loctype(6,fbem_quad8)=fbem_loctype_edge_interior
      fbem_node_loctype(7,fbem_quad8)=fbem_loctype_edge_interior
      fbem_node_loctype(8,fbem_quad8)=fbem_loctype_edge_interior
      ! quad9
      fbem_node_loctype(:,fbem_quad9)=fbem_node_loctype(:,fbem_quad8)
      fbem_node_loctype(9,fbem_quad9)=fbem_loctype_face_interior
      ! tet4
      fbem_node_loctype(1,fbem_tet4)=fbem_loctype_vertex
      fbem_node_loctype(2,fbem_tet4)=fbem_loctype_vertex
      fbem_node_loctype(3,fbem_tet4)=fbem_loctype_vertex
      fbem_node_loctype(4,fbem_tet4)=fbem_loctype_vertex
      ! tet10
      fbem_node_loctype( :,fbem_tet10)=fbem_node_loctype(:,fbem_tet4)
      fbem_node_loctype( 5,fbem_tet10)=fbem_loctype_edge_interior
      fbem_node_loctype( 6,fbem_tet10)=fbem_loctype_edge_interior
      fbem_node_loctype( 7,fbem_tet10)=fbem_loctype_edge_interior
      fbem_node_loctype( 8,fbem_tet10)=fbem_loctype_edge_interior
      fbem_node_loctype( 9,fbem_tet10)=fbem_loctype_edge_interior
      fbem_node_loctype(10,fbem_tet10)=fbem_loctype_edge_interior
!      ! Crack-tip elements derived from line3
!      fbem_node_loctype(:,fbem_line3_qp1t )=fbem_node_loctype(:,fbem_line3)
!      fbem_node_loctype(:,fbem_line3_qp2t )=fbem_node_loctype(:,fbem_line3)
!      fbem_node_loctype(:,fbem_line3_mqp1u)=fbem_node_loctype(:,fbem_line3)
!      fbem_node_loctype(:,fbem_line3_mqp2u)=fbem_node_loctype(:,fbem_line3)
!      fbem_node_loctype(:,fbem_line3_mqp1t)=fbem_node_loctype(:,fbem_line3)
!      fbem_node_loctype(:,fbem_line3_mqp2t)=fbem_node_loctype(:,fbem_line3)
      ! tri4
      fbem_node_loctype(1,fbem_tri4)=fbem_loctype_vertex
      fbem_node_loctype(2,fbem_tri4)=fbem_loctype_vertex
      fbem_node_loctype(3,fbem_tri4)=fbem_loctype_vertex
      fbem_node_loctype(4,fbem_tri4)=fbem_loctype_face_interior
      ! tri1
      fbem_node_loctype(1,fbem_tri1)=fbem_loctype_face_interior
      ! quad1
      fbem_node_loctype(1,fbem_quad1)=fbem_loctype_face_interior
    end if
    if (.not.allocated(fbem_n_edges)) then
      allocate (fbem_n_edges(maxtype))
      fbem_n_edges=0
      fbem_n_edges(fbem_point1     )=0
      fbem_n_edges(fbem_line1      )=0
      fbem_n_edges(fbem_line2      )=1
      fbem_n_edges(fbem_line3      )=1
      fbem_n_edges(fbem_line4      )=1
      fbem_n_edges(fbem_tri3       )=3
      fbem_n_edges(fbem_tri6       )=3
      fbem_n_edges(fbem_quad4      )=4
      fbem_n_edges(fbem_quad8      )=4
      fbem_n_edges(fbem_quad9      )=4
      fbem_n_edges(fbem_tet4       )=6
      fbem_n_edges(fbem_tet10      )=6
!      fbem_n_edges(fbem_line3_qp1t )=1
!      fbem_n_edges(fbem_line3_qp2t )=1
!      fbem_n_edges(fbem_line3_mqp1u)=1
!      fbem_n_edges(fbem_line3_mqp2u)=1
!      fbem_n_edges(fbem_line3_mqp1t)=1
!      fbem_n_edges(fbem_line3_mqp2t)=1
      fbem_n_edges(fbem_tri4       )=3
      fbem_n_edges(fbem_tri1       )=0
      fbem_n_edges(fbem_quad1      )=0
    end if
    maxedges=maxval(fbem_n_edges)
    if (.not.allocated(fbem_edge_type)) then
      allocate (fbem_edge_type(maxedges,maxtype))
      fbem_edge_type=0
      ! line2
      fbem_edge_type(1,fbem_line2)=fbem_line2
      ! line3
      fbem_edge_type(1,fbem_line3)=fbem_line3
      ! line4
      fbem_edge_type(1,fbem_line4)=fbem_line4
      ! tri3
      fbem_edge_type(1,fbem_tri3)=fbem_line2
      fbem_edge_type(2,fbem_tri3)=fbem_line2
      fbem_edge_type(3,fbem_tri3)=fbem_line2
      ! tri6
      fbem_edge_type(1,fbem_tri6)=fbem_line3
      fbem_edge_type(2,fbem_tri6)=fbem_line3
      fbem_edge_type(3,fbem_tri6)=fbem_line3
      ! quad4
      fbem_edge_type(1,fbem_quad4)=fbem_line2
      fbem_edge_type(2,fbem_quad4)=fbem_line2
      fbem_edge_type(3,fbem_quad4)=fbem_line2
      fbem_edge_type(4,fbem_quad4)=fbem_line2
      ! quad8
      fbem_edge_type(1,fbem_quad8)=fbem_line3
      fbem_edge_type(2,fbem_quad8)=fbem_line3
      fbem_edge_type(3,fbem_quad8)=fbem_line3
      fbem_edge_type(4,fbem_quad8)=fbem_line3
      ! quad9
      fbem_edge_type(:,fbem_quad9)=fbem_edge_type(:,fbem_quad8)
      ! tet4
      fbem_edge_type(1,fbem_tet4)=fbem_line2
      fbem_edge_type(2,fbem_tet4)=fbem_line2
      fbem_edge_type(3,fbem_tet4)=fbem_line2
      fbem_edge_type(4,fbem_tet4)=fbem_line2
      fbem_edge_type(5,fbem_tet4)=fbem_line2
      fbem_edge_type(6,fbem_tet4)=fbem_line2
      ! tet10
      fbem_edge_type(1,fbem_tet10)=fbem_line3
      fbem_edge_type(2,fbem_tet10)=fbem_line3
      fbem_edge_type(3,fbem_tet10)=fbem_line3
      fbem_edge_type(4,fbem_tet10)=fbem_line3
      fbem_edge_type(5,fbem_tet10)=fbem_line3
      fbem_edge_type(6,fbem_tet10)=fbem_line3
!      ! Crack-tip elements derived from line3
!      fbem_edge_type(:,fbem_line3_qp1t )=fbem_edge_type(:,fbem_line3)
!      fbem_edge_type(:,fbem_line3_qp2t )=fbem_edge_type(:,fbem_line3)
!      fbem_edge_type(:,fbem_line3_mqp1u)=fbem_edge_type(:,fbem_line3)
!      fbem_edge_type(:,fbem_line3_mqp2u)=fbem_edge_type(:,fbem_line3)
!      fbem_edge_type(:,fbem_line3_mqp1t)=fbem_edge_type(:,fbem_line3)
!      fbem_edge_type(:,fbem_line3_mqp2t)=fbem_edge_type(:,fbem_line3)
      ! tri3
      fbem_edge_type(1,fbem_tri4)=fbem_line2
      fbem_edge_type(2,fbem_tri4)=fbem_line2
      fbem_edge_type(3,fbem_tri4)=fbem_line2
    end if
    ! Change by hand
    maxnodesedge=4
    if (.not.allocated(fbem_edge_node)) then
      allocate (fbem_edge_node(maxnodesedge,maxedges,maxtype))
      fbem_edge_node=0
      ! line2
      fbem_edge_node(1,1,fbem_line2)=1
      fbem_edge_node(2,1,fbem_line2)=2
      ! line3
      fbem_edge_node(:,:,fbem_line3)=fbem_edge_node(:,:,fbem_line2)
      fbem_edge_node(3,1,fbem_line3)=3
      ! line4
      fbem_edge_node(:,:,fbem_line4)=fbem_edge_node(:,:,fbem_line2)
      fbem_edge_node(3,1,fbem_line4)=3
      fbem_edge_node(4,1,fbem_line4)=4
      ! tri3
      fbem_edge_node(1,1,fbem_tri3)=1
      fbem_edge_node(2,1,fbem_tri3)=2
      fbem_edge_node(1,2,fbem_tri3)=2
      fbem_edge_node(2,2,fbem_tri3)=3
      fbem_edge_node(1,3,fbem_tri3)=3
      fbem_edge_node(2,3,fbem_tri3)=1
      ! tri6
      fbem_edge_node(:,:,fbem_tri6)=fbem_edge_node(:,:,fbem_tri3)
      fbem_edge_node(3,1,fbem_tri6)=4
      fbem_edge_node(3,2,fbem_tri6)=5
      fbem_edge_node(3,3,fbem_tri6)=6
      ! quad4
      fbem_edge_node(1,1,fbem_quad4)=1
      fbem_edge_node(2,1,fbem_quad4)=2
      fbem_edge_node(1,2,fbem_quad4)=2
      fbem_edge_node(2,2,fbem_quad4)=3
      fbem_edge_node(1,3,fbem_quad4)=3
      fbem_edge_node(2,3,fbem_quad4)=4
      fbem_edge_node(1,4,fbem_quad4)=4
      fbem_edge_node(2,4,fbem_quad4)=1
      ! quad8
      fbem_edge_node(:,:,fbem_quad8)=fbem_edge_node(:,:,fbem_quad4)
      fbem_edge_node(3,1,fbem_quad8)=5
      fbem_edge_node(3,2,fbem_quad8)=6
      fbem_edge_node(3,3,fbem_quad8)=7
      fbem_edge_node(3,4,fbem_quad8)=8
      ! quad9
      fbem_edge_node(:,:,fbem_quad9)=fbem_edge_node(:,:,fbem_quad8)
      ! tet4
      fbem_edge_node(1,1,fbem_tet4)=1
      fbem_edge_node(2,1,fbem_tet4)=2
      fbem_edge_node(1,2,fbem_tet4)=2
      fbem_edge_node(2,2,fbem_tet4)=3
      fbem_edge_node(1,3,fbem_tet4)=3
      fbem_edge_node(2,3,fbem_tet4)=1
      fbem_edge_node(1,4,fbem_tet4)=1
      fbem_edge_node(2,4,fbem_tet4)=4
      fbem_edge_node(1,5,fbem_tet4)=2
      fbem_edge_node(2,5,fbem_tet4)=4
      fbem_edge_node(1,6,fbem_tet4)=3
      fbem_edge_node(2,6,fbem_tet4)=4
      ! tet10
      fbem_edge_node(:,:,fbem_tet10)=fbem_edge_node(:,:,fbem_tet4)
      fbem_edge_node(3,1,fbem_tet10)=5
      fbem_edge_node(3,2,fbem_tet10)=6
      fbem_edge_node(3,3,fbem_tet10)=7
      fbem_edge_node(3,4,fbem_tet10)=8
      fbem_edge_node(3,5,fbem_tet10)=9
      fbem_edge_node(3,6,fbem_tet10)=10
      !
!      ! Crack-tip elements derived from line3
!      fbem_edge_node(:,:,fbem_line3_qp1t )=fbem_edge_node(:,:,fbem_line3)
!      fbem_edge_node(:,:,fbem_line3_qp2t )=fbem_edge_node(:,:,fbem_line3)
!      fbem_edge_node(:,:,fbem_line3_mqp1u)=fbem_edge_node(:,:,fbem_line3)
!      fbem_edge_node(:,:,fbem_line3_mqp2u)=fbem_edge_node(:,:,fbem_line3)
!      fbem_edge_node(:,:,fbem_line3_mqp1t)=fbem_edge_node(:,:,fbem_line3)
!      fbem_edge_node(:,:,fbem_line3_mqp2t)=fbem_edge_node(:,:,fbem_line3)
      ! tri4
      fbem_edge_node(1,1,fbem_tri4)=1
      fbem_edge_node(2,1,fbem_tri4)=2
      fbem_edge_node(1,2,fbem_tri4)=2
      fbem_edge_node(2,2,fbem_tri4)=3
      fbem_edge_node(1,3,fbem_tri4)=3
      fbem_edge_node(2,3,fbem_tri4)=1
    end if

    ! Faces and so on

  end subroutine fbem_init_shape_functions_module

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Nodes coordinates in reference space
  !---------------------------------------------------------------------------------------------------------------------------------

  !! Nodes coordinates in reference space for 1D elements
  function fbem_xi(etype,delta)
    implicit none
    !! Element geometrical interpolation type: <tt>etype={fbem_line2,fbem_line3}</tt>.
    integer :: etype
    !! Displacement for discontinuous elements
    real(kind=real64) :: delta
    !! Reference coordinates of the element nodes of an <tt>etype</tt> element.
    real(kind=real64) :: fbem_xi(fbem_n_nodes(etype))
#   define xi fbem_xi
#   include <xi_1d.rc>
#   undef xi
  end function fbem_xi

  !! Nodes coordinates in reference space for 2D elements
  function fbem_xi1xi2(etype,delta)
    implicit none
    !! Element geometrical interpolation type: <tt>etype={fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    integer :: etype
    !! Displacement for discontinuous elements
    real(kind=real64) :: delta
    !! Matrix that contains reference coordinates of the element nodes of an <tt>etype</tt> element.
    real(kind=real64) :: fbem_xi1xi2(2,fbem_n_nodes(etype))
#   define xi fbem_xi1xi2
#   include <xi_2d.rc>
#   undef xi
  end function fbem_xi1xi2

  !! Nodes coordinates in reference space for 2D elements
  function fbem_xi_hybrid(etype,delta)
    implicit none
    integer           :: etype                                                        !! Element geometrical interpolation type: <tt>etype={fbem_line2,fbem_line3,fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    real(kind=real64) :: delta                                                        !! Displacement for discontinuous elements
    real(kind=real64) :: fbem_xi_hybrid(fbem_n_dimension(etype),fbem_n_nodes(etype)) !! Matrix that contains reference coordinates of the element nodes of an <tt>etype</tt> element.
    real(kind=real64) :: xiaux(fbem_n_nodes(etype)) ! Auxiliary variable to store the xi coordinates for a 1D element
    integer           :: kn                        ! Counter
    ! Select depending on the number of dimensions of the element
    select case (fbem_n_dimension(etype))
      case (1)
#       define xi xiaux
#       include <xi_1d.rc>
#       undef xi
        do kn=1,fbem_n_nodes(etype)
          fbem_xi_hybrid(1,kn)=xiaux(kn)
        end do
      case (2)
#       define xi fbem_xi_hybrid
#       include <xi_2d.rc>
#       undef xi
      case (3)
        if (delta.ne.0.d0) stop '3D elements only with delta=0'
#       define xi fbem_xi_hybrid
#       include <xi_3d.rc>
#       undef xi
    end select
  end function fbem_xi_hybrid

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Check if coordinates in reference space are in the domain
  !---------------------------------------------------------------------------------------------------------------------------------

  !! Check if given <tt>xi</tt> coordinate is in the element
  function fbem_check_xi(xi)
    implicit none
    !! <tt>.true.</tt> if <tt>xi</tt> is in the element, <tt>.false.</tt> otherwise.
    logical :: fbem_check_xi
    !! Coordinate in reference space.
    real(kind=real64) :: xi
    if ((xi.lt.(-1.0d0-check_xi_tolerance)).or.(xi.gt.(1.0d0+check_xi_tolerance))) then
      fbem_check_xi=.false.
    else
      fbem_check_xi=.true.
    end if
  end function fbem_check_xi

  !! Check if given <tt>xi</tt> coordinate is in the element and print a warning message if not.
  subroutine fbem_check_xi_warning(xi)
    implicit none
    !! Coordinate in reference space.
    real(kind=real64) :: xi
    if ((xi.lt.(-1.0d0-check_xi_tolerance)).or.(xi.gt.(1.0d0+check_xi_tolerance))) then
      call fbem_warning_message(error_unit,1,__FILE__,__LINE__,'xi not in element domain')
    end if
  end subroutine fbem_check_xi_warning

  !! Check if given <tt>xi</tt> coordinate is in the element, if not print an error message and execution stops.
  subroutine fbem_check_xi_error(xi)
    implicit none
    !! Coordinate in reference space.
    real(kind=real64) :: xi
    if ((xi.lt.(-1.0d0-check_xi_tolerance)).or.(xi.gt.(1.0d0+check_xi_tolerance))) then
      call fbem_error_message(output_unit,0,__FILE__,__LINE__,'xi not in element domain')
    end if
  end subroutine fbem_check_xi_error

  !! Check if given <tt>xi</tt> coordinate is at a vertex
  function fbem_check_xi_vertex(xi)
    implicit none
    logical :: fbem_check_xi_vertex
    !! Coordinate in reference space.
    real(kind=real64) :: xi
    fbem_check_xi_vertex=.false.
    if ((xi.lt.-1.0d0+check_xi_tolerance).or.(xi.gt.1.0d0-check_xi_tolerance)) then
      fbem_check_xi_vertex=.true.
    end if
  end function fbem_check_xi_vertex

  !! Check if given <tt>xi</tt> coordinate is in the element
  function fbem_check_xi1xi2(etype,xi)
    implicit none
    logical :: fbem_check_xi1xi2
    !! Element interpolation type: <tt>etype={fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    integer :: etype
    !! Coordinates in reference space.
    real(kind=real64) :: xi(fbem_n_dimension(etype))
    select case (etype)
      case (fbem_quad4,fbem_quad8,fbem_quad9)
        if ((xi(1).lt.-1.0d0-check_xi_tolerance).or.(xi(1).gt.1.0d0+check_xi_tolerance).or.&
            (xi(2).lt.-1.0d0-check_xi_tolerance).or.(xi(2).gt.1.0d0+check_xi_tolerance)) then
          fbem_check_xi1xi2=.false.
        else
          fbem_check_xi1xi2=.true.
        end if
      case (fbem_tri3,fbem_tri6)
        if ((xi(1).lt.0.0d0-check_xi_tolerance).or.(xi(2).lt.0.0d0-check_xi_tolerance).or.&
            ((xi(1)+xi(2)).gt.1.0d0+check_xi_tolerance)) then
          fbem_check_xi1xi2=.false.
        else
          fbem_check_xi1xi2=.true.
        end if
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={tri3,tri6,quad4,quad8,quad9}')
    end select
  end function fbem_check_xi1xi2

  !! Check if given <tt>xi</tt> coordinate is in the element and print a warning message if not.
  subroutine fbem_check_xi1xi2_warning(etype,xi)
    implicit none
    !! Element interpolation type: <tt>etype={fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    integer :: etype
    !! Coordinates in reference space.
    real(kind=real64) :: xi(fbem_n_dimension(etype))
    select case (etype)
      case (fbem_quad4,fbem_quad8,fbem_quad9)
        if ((xi(1).lt.-1.0d0-check_xi_tolerance).or.(xi(1).gt.1.0d0+check_xi_tolerance).or.&
            (xi(2).lt.-1.0d0-check_xi_tolerance).or.(xi(2).gt.1.0d0+check_xi_tolerance)) then
          call fbem_warning_message(error_unit,1,__FILE__,__LINE__,'(xi_1,xi_2) not in element domain')
        end if
      case (fbem_tri3,fbem_tri6)
        if ((xi(1).lt.0.0d0-check_xi_tolerance).or.(xi(2).lt.0.0d0-check_xi_tolerance).or.&
            ((xi(1)+xi(2)).gt.1.0d0+check_xi_tolerance)) then
          call fbem_warning_message(error_unit,1,__FILE__,__LINE__,'(xi_1,xi_2) not in element domain')
        end if
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={tri3,tri6,quad4,quad8,quad9}')
    end select
  end subroutine fbem_check_xi1xi2_warning

  !! Check if given <tt>xi</tt> coordinate is in the element, if not print an error message and execution stops.
  subroutine fbem_check_xi1xi2_error(etype,xi)
    implicit none
    !! Element interpolation type: <tt>etype={fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    integer :: etype
    !! Coordinates in reference space.
    real(kind=real64) :: xi(fbem_n_dimension(etype))
    select case (etype)
      case (fbem_quad4,fbem_quad8,fbem_quad9)
        if ((xi(1).lt.-1.0d0-check_xi_tolerance).or.(xi(1).gt.1.0d0+check_xi_tolerance).or.&
            (xi(2).lt.-1.0d0-check_xi_tolerance).or.(xi(2).gt.1.0d0+check_xi_tolerance)) then
          call fbem_error_message(output_unit,0,__FILE__,__LINE__,'(xi_1,xi_2) not in element domain')
        end if
      case (fbem_tri3,fbem_tri6)
        if ((xi(1).lt.0.0d0-check_xi_tolerance).or.(xi(2).lt.0.0d0-check_xi_tolerance).or.&
            ((xi(1)+xi(2)).gt.1.0d0+check_xi_tolerance)) then
          call fbem_error_message(output_unit,0,__FILE__,__LINE__,'(xi_1,xi_2) not in element domain')
        end if
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={tri3,tri6,quad4,quad8,quad9}')
    end select
  end subroutine fbem_check_xi1xi2_error

  !! Check if given <tt>xi</tt> coordinate is on the edge of the element
  function fbem_check_xi1xi2_edge(etype,xi)
    implicit none
    logical :: fbem_check_xi1xi2_edge
    !! Element interpolation type: <tt>etype={fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    integer :: etype
    !! Coordinates in reference space.
    real(kind=real64) :: xi(fbem_n_dimension(etype))
    select case (etype)
      case (fbem_quad4,fbem_quad8,fbem_quad9)
        if ((xi(1).lt.-1.0d0-check_xi_tolerance).or.(xi(1).gt.1.0d0+check_xi_tolerance).or.&
            (xi(2).lt.-1.0d0-check_xi_tolerance).or.(xi(2).gt.1.0d0+check_xi_tolerance)) then
          fbem_check_xi1xi2_edge=.false.
        else
          if ((xi(1).lt.-1.0d0+check_xi_tolerance).or.(xi(1).gt.1.0d0-check_xi_tolerance).or.&
              (xi(2).lt.-1.0d0+check_xi_tolerance).or.(xi(2).gt.1.0d0-check_xi_tolerance)) then
            fbem_check_xi1xi2_edge=.true.
          else
            fbem_check_xi1xi2_edge=.false.
          end if
        end if
      case (fbem_tri3,fbem_tri6)
        if ((xi(1).lt.0.0d0-check_xi_tolerance).or.(xi(2).lt.0.0d0-check_xi_tolerance).or.&
            ((xi(1)+xi(2)).gt.1.0d0+check_xi_tolerance)) then
          fbem_check_xi1xi2_edge=.false.
        else
          if ((xi(1).lt.check_xi_tolerance).or.(xi(2).lt.check_xi_tolerance).or.&
              ((xi(1)+xi(2)).gt.1.0d0-check_xi_tolerance)) then
            fbem_check_xi1xi2_edge=.true.
          else
            fbem_check_xi1xi2_edge=.false.
          end if
        end if
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={tri3,tri6,quad4,quad8,quad9}')
    end select
  end function fbem_check_xi1xi2_edge

  !! Check if a given <tt>xi</tt> coordinate is in the element boundary
  function fbem_check_xi_on_element_boundary(etype,xi)
    implicit none
    logical           :: fbem_check_xi_on_element_boundary
    integer           :: etype                       !! Element interpolation type: <tt>etype={fbem_tri3,fbem_tri6,fbem_quad4,fbem_quad8,fbem_quad9}</tt>.
    real(kind=real64) :: xi(fbem_n_dimension(etype)) !! Coordinates in reference space.
    select case (etype)
      case (fbem_line2,fbem_line3,fbem_line4)
        fbem_check_xi_on_element_boundary=.false.
        if ((xi(1).lt.-1.0d0+check_xi_tolerance).or.(xi(1).gt.1.0d0-check_xi_tolerance)) then
          fbem_check_xi_on_element_boundary=.true.
        end if
      case (fbem_quad4,fbem_quad8,fbem_quad9)
        if ((xi(1).lt.-1.0d0-check_xi_tolerance).or.(xi(1).gt.1.0d0+check_xi_tolerance).or.&
            (xi(2).lt.-1.0d0-check_xi_tolerance).or.(xi(2).gt.1.0d0+check_xi_tolerance)) then
          fbem_check_xi_on_element_boundary=.false.
        else
          if ((xi(1).lt.-1.0d0+check_xi_tolerance).or.(xi(1).gt.1.0d0-check_xi_tolerance).or.&
              (xi(2).lt.-1.0d0+check_xi_tolerance).or.(xi(2).gt.1.0d0-check_xi_tolerance)) then
            fbem_check_xi_on_element_boundary=.true.
          else
            fbem_check_xi_on_element_boundary=.false.
          end if
        end if
      case (fbem_tri3,fbem_tri6)
        if ((xi(1).lt.0.0d0-check_xi_tolerance).or.(xi(2).lt.0.0d0-check_xi_tolerance).or.&
            ((xi(1)+xi(2)).gt.1.0d0+check_xi_tolerance)) then
          fbem_check_xi_on_element_boundary=.false.
        else
          if ((xi(1).lt.check_xi_tolerance).or.(xi(2).lt.check_xi_tolerance).or.&
              ((xi(1)+xi(2)).gt.1.0d0-check_xi_tolerance)) then
            fbem_check_xi_on_element_boundary=.true.
          else
            fbem_check_xi_on_element_boundary=.false.
          end if
        end if
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={line2,line3,line4,tri3,tri6,quad4,quad8,quad9}')
    end select
  end function fbem_check_xi_on_element_boundary

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Shape functions
  !---------------------------------------------------------------------------------------------------------------------------------

  function fbem_phi_1d(etype,delta,xi)
    implicit none
    !! Element type
    integer :: etype
    !! Displacement for discontinuous elements
    real(kind=real64) :: delta
    !! Coordinate in reference space
    real(kind=real64) :: xi
    real(kind=real64) :: fbem_phi_1d(fbem_n_nodes(etype))
    real(kind=real64) :: aux(10)
#   if CHECK_XIS==1
    call fbem_check_xi_warning(xi)
#   endif
#   if CHECK_XIS==2
    call fbem_check_xi_error(xi)
#   endif
#   define phi fbem_phi_1d
#   include <phi_1d.rc>
#   undef phi
  end function fbem_phi_1d

  function fbem_phi_nd(etype,delta,xi)
    implicit none
    !! Element type
    integer :: etype
    !! Displacement for discontinuous elements (if delta=0.0d0, then continuous element)
    real(kind=real64) :: delta
    !! Coordinates in reference space
    real(kind=real64) :: xi(fbem_n_dimension(etype))
    real(kind=real64) :: fbem_phi_nd(fbem_n_nodes(etype))
    real(kind=real64) :: aux(10)
#   if CHECK_XIS==1
    call fbem_check_xi1xi2_warning(etype,xi)
#   endif
#   if CHECK_XIS==1
    call fbem_check_xi1xi2_error(etype,xi)
#   endif
#   define phi fbem_phi_nd
#   include <phi_2d.rc>
#   undef phi
  end function fbem_phi_nd

  function fbem_phi_hybrid(etype,delta,xi_p)
    implicit none
    ! I/O
    integer           :: etype                                !! Element type
    real(kind=real64) :: delta                                !! Displacement for discontinuous elements (if delta=0.0d0, then continuous element)
    real(kind=real64) :: xi_p(fbem_n_dimension(etype))        !! Coordinates in reference space
    real(kind=real64) :: fbem_phi_hybrid(fbem_n_nodes(etype)) !! Shape functions values
    ! Local
    real(kind=real64) :: aux(10)
#   if CHECK_XIS==1
    call fbem_check_xi1xi2_warning(etype,xi)
#   endif
#   if CHECK_XIS==1
    call fbem_check_xi1xi2_error(etype,xi)
#   endif
    select case (fbem_n_dimension(etype))
      case (1)
#       define xi xi_p(1)
#       define phi fbem_phi_hybrid
#       include <phi_1d.rc>
#       undef xi
#       undef phi
      case (2)
#       define xi xi_p
#       define phi fbem_phi_hybrid
#       include <phi_2d.rc>
#       undef xi
#       undef phi
      case default
        stop 'fbem_phi_hybrid: fbem_n_dimension(etype) not valid'
    end select
  end function fbem_phi_hybrid

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Shape functions first derivatives
  !---------------------------------------------------------------------------------------------------------------------------------

  function fbem_dphidxi(etype,delta,xi)
    implicit none
    !! Element type
    integer :: etype
    !! Displacement for discontinuous elements (if delta=0.0d0, then continuous element)
    real(kind=real64) :: delta
    !! Coordinate in reference space
    real(kind=real64) :: xi
    real(kind=real64) :: fbem_dphidxi(fbem_n_nodes(etype))
    real(kind=real64) :: aux(10)
#   if CHECK_XIS==1
    call fbem_check_xi_warning(xi)
#   endif
#   if CHECK_XIS==2
    call fbem_check_xi_error(xi)
#   endif
#   define dphidxi fbem_dphidxi
#   include <dphidxi_1d.rc>
#   undef dphidxi
  end function fbem_dphidxi

  function fbem_dphidxi1(etype,delta,xi)
    implicit none
    !! Element type
    integer :: etype
    !! Displacement for discontinuous elements (if delta=0.0d0, then continuous element)
    real(kind=real64) :: delta
    !! Coordinates in reference space
    real(kind=real64) :: xi(fbem_n_dimension(etype))
    real(kind=real64) :: fbem_dphidxi1(fbem_n_nodes(etype))
    real(kind=real64) :: aux(10)
#   if CHECK_XIS==1
    call fbem_check_xi1xi2_warning(etype,xi)
#   endif
#   if CHECK_XIS==1
    call fbem_check_xi1xi2_error(etype,xi)
#   endif
#   define dphidxi1 fbem_dphidxi1
#   include <dphidxi1_2d.rc>
#   undef dphidxi1
  end function fbem_dphidxi1

  function fbem_dphidxi2(etype,delta,xi)
    implicit none
    !! Element type
    integer :: etype
    !! Displacement for discontinuous elements (if delta=0.0d0, then continuous element)
    real(kind=real64) :: delta
    !! Coordinates in reference space
    real(kind=real64) :: xi(fbem_n_dimension(etype))
    real(kind=real64) :: fbem_dphidxi2(fbem_n_nodes(etype))
    real(kind=real64) :: aux(10)
#   if CHECK_XIS==1
    call fbem_check_xi1xi2_warning(etype,xi)
#   endif
#   if CHECK_XIS==1
    call fbem_check_xi1xi2_error(etype,xi)
#   endif
#   define dphidxi2 fbem_dphidxi2
#   include <dphidxi2_2d.rc>
#   undef dphidxi2
  end function fbem_dphidxi2

  !! Shape function derivatives with respect to cartesian coordinates.
  function fbem_dphidx(rn,etype,delta,x,xind)
    implicit none
    ! I/O
    integer           :: rn
    integer           :: etype
    real(kind=real64) :: delta
    real(kind=real64) :: x(rn,fbem_n_nodes(etype))
    real(kind=real64) :: xind(fbem_n_dimension(etype))
    real(kind=real64) :: fbem_dphidx(fbem_n_nodes(etype),rn)
    ! Local
    real(kind=real64) :: aux(10)
    real(kind=real64) :: dphidxi(fbem_n_nodes(etype))
    real(kind=real64) :: dphidxi1(fbem_n_nodes(etype))
    real(kind=real64) :: dphidxi2(fbem_n_nodes(etype))
    real(kind=real64) :: xi1d
    real(kind=real64) :: T(rn), dxdxi(rn,rn), dxidx(rn,rn), norm
    real(kind=real64) :: jxi1, jxi2
    real(kind=real64) :: N(3), v1(3), v2(3), t1(3), t2(3)
    real(kind=real64) :: M(2,2), J(2,2)
    integer           :: n_nodes, n_dimension, i, k
    n_nodes=fbem_n_nodes(etype)
    n_dimension=fbem_n_dimension(etype)
    select case (n_dimension)
      !
      ! LINE ELEMENT
      !
      case (1)
        xi1d=xind(1)
#       define xi xi1d
#       include <dphidxi_1d.rc>
#       undef xi
        T=0.d0
        do i=1,n_nodes
          T=T+dphidxi(i)*x(:,i)
        end do
        norm=sqrt(dot_product(T,T))
        do k=1,rn
          fbem_dphidx(:,k)=T(k)*dphidxi
        end do
        fbem_dphidx=fbem_dphidx/norm**2
      !
      ! SURFACE ELEMENT
      !
      case (2)
#       define xi xind
#       include <dphidxik_2d.rc>
#       undef xi
        dxdxi=0.d0
        do i=1,n_nodes
          dxdxi(1,:)=dxdxi(1,:)+dphidxi1(i)*x(:,i)
          dxdxi(2,:)=dxdxi(2,:)+dphidxi2(i)*x(:,i)
        end do
        select case (rn)
          case (2)
            call fbem_invert_2x2_matrix(dxdxi,dxidx,norm)
            fbem_dphidx(:,1)=dxidx(1,1)*dphidxi1(:)+dxidx(1,2)*dphidxi2(:)
            fbem_dphidx(:,2)=dxidx(2,1)*dphidxi1(:)+dxidx(2,2)*dphidxi2(:)
          case (3)
            ! Line jacobians: J_{xi_j} = |dx/dxi_j|
            jxi1=sqrt(dot_product(dxdxi(1,:),dxdxi(1,:)))
            jxi2=sqrt(dot_product(dxdxi(2,:),dxdxi(2,:)))
            ! Unit tangents
            t1=dxdxi(1,:)/jxi1
            t2=dxdxi(2,:)/jxi2
            ! v1 = t1
            v1=t1
            ! Unit normal vector: n = t1 x t2 / |t1 x t2|
            N(1)=t1(2)*t2(3)-t1(3)*t2(2)
            N(2)=t1(3)*t2(1)-t1(1)*t2(3)
            N(3)=t1(1)*t2(2)-t1(2)*t2(1)
            n=N/sqrt(dot_product(N,N))
            ! v2 = n x v1
            v2(1)=n(2)*v1(3)-n(3)*v1(2)
            v2(2)=n(3)*v1(1)-n(1)*v1(3)
            v2(3)=n(1)*v1(2)-n(2)*v1(1)
            ! M_{ij} = ti · vj
            M(1,1)=dot_product(t1,v1)
            M(1,2)=dot_product(t1,v2)
            M(2,1)=dot_product(t2,v1)
            M(2,2)=dot_product(t2,v2)
            ! Jacobian matrix: d/dv_i = J_{ij} d/dxi_j
            J(1,1)= M(2,2)/jxi1
            J(1,2)=-M(1,2)/jxi1
            J(2,1)=-M(2,1)/jxi2
            J(2,2)= M(1,1)/jxi2
            J=J/(M(1,1)*M(2,2)-M(1,2)*M(2,1))
            ! L · J
            do i=1,n_nodes
              fbem_dphidx(i,:)=(J(1,1)*dphidxi1(i)+J(1,2)*dphidxi2(i))*v1&
                              +(J(2,1)*dphidxi1(i)+J(2,2)*dphidxi2(i))*v2
            end do
        end select
      case (3)
        stop 'not yet : fbem_dphidx'
    end select
  end function fbem_dphidx


  ! --------------------------------------------------------------------------------------------------------------------------------
  ! Shape function integrals
  ! --------------------------------------------------------------------------------------------------------------------------------

  !! Calculate phi integrals
  subroutine fbem_shape_functions_integrals(rn,type_g,type_f,delta_f,x_nodes,glp,weight)
    implicit none
    ! I/O
    integer                      :: rn                               !! R^n
    integer                      :: type_g                           !! Geometrical interpolation: line2, line3, tri3, tri6, quad4, quad8, quad9
    integer                      :: type_f                           !! Functional interpolation
    real(kind=real64)            :: delta_f                          !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64)            :: x_nodes(rn,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
    integer                      :: glp                              !! Gauss-Legendre number of points
    real(kind=real64)            :: weight(fbem_n_nodes(type_f))     !! Nodal weight (S_Omega phi·dOmega)
    ! Local
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
    real(kind=real64)            :: T(rn)                            ! Tangent vector at xi
    real(kind=real64)            :: T1(rn), T2(rn), T1f(3), T2f(3)   ! Tangent vectors at xi
    real(kind=real64)            :: N(3)                             ! Normal vector at xi
    real(kind=real64)            :: jg                               ! Geometric jacobian
    real(kind=real64)            :: jw                               ! Jacobian * weight
    real(kind=real64)            :: phi_f(fbem_n_nodes(type_f))      ! Functional shape functions values
    ! Number of nodes of the element
    nnodes_g=fbem_n_nodes(type_g)
    ! Initialize integrals
    weight=0.d0
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
          ! Components calculation of T at xi
          T=0.d0
          do kphi=1,nnodes_g
            T=T+dphidxi_g(kphi)*x_nodes(:,kphi)
          end do
          ! Geometric jacobian
          jg=sqrt(dot_product(T,T))
          ! Functional shape functions at xi
#         define etype type_f
#         define delta delta_f
#         define phi phi_f
#         include <phi_1d.rc>
#         undef etype
#         undef delta
#         undef phi
          ! Add integration point
          jw=jg*w
          weight=weight+phi_f*jw
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
#             undef dphidxi
              ! Components calculation of T1 and T2 at xi
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
              ! Functional shape functions at xi
#             define etype type_f
#             define delta delta_f
#             define xi xi2d
#             define phi phi_f
#             include <phi_2d.rc>
#             undef etype
#             undef delta
#             undef xi
#             undef phi
              ! Add integration point
              jw=jg*w
              weight=weight+phi_f*jw
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
#               undef dphidxi
                ! Components calculation of T1 and T2 at xi
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
                ! Functional shape functions at xi
#               define etype type_f
#               define delta delta_f
#               define xi xi2d
#               define phi phi_f
#               include <phi_2d.rc>
#               undef etype
#               undef delta
#               undef xi
#               undef phi
                ! Add integration point
                jw=jg*w2d(1)*w2d(2)
                weight=weight+phi_f*jw
              end do
            end do
        end select
      !
      ! VOLUME ELEMENTS
      !
      case (3)
        stop 'not implemented yet'
    end select
  end subroutine fbem_shape_functions_integrals

  !! Calculate phi*x integrals
  subroutine fbem_shape_functions_x_integrals(rn,type_g,type_f,delta_f,x_nodes,p,glp,weight)
    implicit none
    ! I/O
    integer                      :: rn                               !! R^n
    integer                      :: type_g                           !! Geometrical interpolation
    integer                      :: type_f                           !! Functional interpolation
    real(kind=real64)            :: delta_f                          !! Displacement for discontinuous functional elements (if delta=0.0d0, then continuous element)
    real(kind=real64)            :: x_nodes(rn,fbem_n_nodes(type_g)) !! Position vectors of geometrical nodes
    real(kind=real64)            :: p(rn)                            !! Point with respect moments are calculated
    integer                      :: glp                              !! Gauss-Legendre number of points
    real(kind=real64)            :: weight(fbem_n_nodes(type_f),rn)  !! Nodal weights (S_Omega phi·dOmega)
    ! Local
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
    real(kind=real64)            :: phi_f(fbem_n_nodes(type_f))      ! Functional shape functions values
    ! Number of nodes of the element
    nnodes_g=fbem_n_nodes(type_g)
    ! Initialize integrals
    weight=0.d0
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
          ! Functional shape functions at xi
#         define etype type_f
#         define delta delta_f
#         define phi phi_f
#         include <phi_1d.rc>
#         undef etype
#         undef delta
#         undef phi
          ! Add integration point
          jw=jg*w
          do kc=1,rn
            weight(:,kc)=weight(:,kc)+(x(kc)-p(kc))*phi_f*jw
          end do
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
#             undef dphidxi
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
              ! Functional shape functions at xi
#             define etype type_f
#             define delta delta_f
#             define xi xi2d
#             define phi phi_f
#             include <phi_2d.rc>
#             undef etype
#             undef delta
#             undef xi
#             undef phi
              ! Add integration point
              jw=jg*w
              do kc=1,rn
                weight(:,kc)=weight(:,kc)+(x(kc)-p(kc))*phi_f*jw
              end do
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
#               undef dphidxi
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
                ! Functional shape functions at xi
#               define etype type_f
#               define delta delta_f
#               define xi xi2d
#               define phi phi_f
#               include <phi_2d.rc>
#               undef etype
#               undef delta
#               undef xi
#               undef phi
                ! Add integration point
                jw=jg*w2d(1)*w2d(2)
                do kc=1,rn
                  weight(:,kc)=weight(:,kc)+(x(kc)-p(kc))*phi_f*jw
                end do
              end do
            end do
        end select
      !
      ! VOLUME ELEMENTS
      !
      case (3)
        stop 'not implemented yet'
    end select
  end subroutine fbem_shape_functions_x_integrals

end module fbem_shape_functions
