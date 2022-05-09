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
!! <b> This module implements several general derived types and routines associated with the BEM.</b>
module fbem_bem_general

  ! Fortran 2003 standard
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_quad_rules
  use fbem_telles_transformation
  use fbem_geometry

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all are private
  private

  public :: fbem_bem_element
  public :: fbem_bem_logd1rphi_1d_int
  public :: fbem_bem_int_logr

  !! BEM element prepared for calculation
  type fbem_bem_element
    ! Geometry
    integer                        :: n                !! Element ambient space dimension (R^n)
    integer                        :: gtype            !! Type of interpolation
    integer                        :: d                !! Element intrinsic space dimension
    integer                        :: n_gnodes         !! Number of nodes
    real(kind=real64), allocatable :: x(:,:)           !! Coordinates of each node
    ! Primary field
    integer                        :: ptype            !! Type of interpolation
    real(kind=real64)              :: ptype_delta      !! Delta parameter for the interpolation (>0 discontinuous)
    integer                        :: n_pnodes         !! Number of nodes
    ! Secondary field (load field for body loads)
    integer                        :: stype            !! Type of interpolation
    real(kind=real64)              :: stype_delta      !! Delta parameter for the interpolation (>0 discontinuous)
    integer                        :: n_snodes         !! Number of nodes
    ! Auxiliary variables
    real(kind=real64)              :: cl               !! Characteristic length
    integer                        :: gln_far          !! 1D Gauss-Legendre integ. points required to integrate any shape function
    real(kind=real64), allocatable :: bball_centre(:)  !! Centre of a ball that contains the element
    real(kind=real64)              :: bball_radius     !! Radius of a ball that contains the element
    ! Design macro-element
    integer                        :: dmetype          !! Type of interpolation
    integer                        :: dme_d            !! Local dimension
    integer                        :: dme_n_gnodes     !! Number of nodes
    real(kind=real64), allocatable :: dme_x(:,:)       !! Coordinates of each node in the ambient space
    real(kind=real64)              :: dme_cl           !! Characteristic length
    ! Precalculated dataset setup
    integer                        :: n_ps                 !! Number of precalculated datasets
    integer          , allocatable :: ps_gln(:)            !! 1D Gauss-Legendre integ. points for each precalculated dataset
    integer                        :: ps_gln_max           !! maxval(ps_gln)
    integer          , allocatable :: ps_ngp(:)            !! Number of integration points of each precalculated dataset
    integer                        :: ps_ngp_max           !! maxval(ps_ngp)
    real(kind=real64), allocatable :: ps_x(:,:,:)          !! Coordinates at each integration point of each precalculated dataset
    real(kind=real64), allocatable :: ps_t(:,:,:,:)        !! Unit tangents at each integration point of each precalculated dataset ps_dme_t(kc,tang,ip,rule)
    real(kind=real64), allocatable :: ps_n(:,:,:)          !! Unit normal at each integration point of each precalculated dataset
    real(kind=real64), allocatable :: ps_j(:,:)            !! Geometric jacobian at each integration point of each precalculated dataset
    real(kind=real64), allocatable :: ps_w(:,:)            !! Weight at each integration point of each precalculated dataset
    real(kind=real64), allocatable :: ps_gphi(:,:,:)       !! phi^g at each integration point of each precalculated dataset
    real(kind=real64), allocatable :: ps_dgphidxi(:,:,:,:) !! dphi^g/dxi_k at each integration point of each precalculated dataset (ps_dgphidxi(node,k,ip,rule))
    real(kind=real64), allocatable :: ps_dxdxi(:,:,:,:)    !! dx_i/dxi_j (jacobian matrix) at each integration point of each precalculated dataset (ps_dxdxi(i,j,ip,rule))
    real(kind=real64), allocatable :: ps_dxidx(:,:,:,:)    !! dxi_i/dx_j (inverse of the jacobian matrix) at each integration point of each precalculated dataset (ps_dxidx(i,j,ip,rule))
    real(kind=real64), allocatable :: ps_pphi(:,:,:)       !! phi^p at each integration point of each precalculated dataset
    real(kind=real64), allocatable :: ps_sphi(:,:,:)       !! phi^s at each integration point of each precalculated dataset
    real(kind=real64), allocatable :: ps_pphijw(:,:,:)     !! phi^p*j*weight at each integration point of each precalculated dataset
    real(kind=real64), allocatable :: ps_sphijw(:,:,:)     !! phi^s*j*weight at each integration point of each precalculated dataset
    ! On the design macro-element in which the element is embedded
    real(kind=real64), allocatable :: ps_dme_gphi(:,:,:)       !! phi^{dme-g} at each integration point of each precalculated dataset
    real(kind=real64), allocatable :: ps_dme_dgphidxi(:,:,:,:) !! dphi^{dme-g}/dxik at each integration point of each precalculated dataset (ps_dme_gdphidxi(node,k,ip,rule))
    real(kind=real64), allocatable :: ps_dme_dxdxi(:,:,:,:)    !! dx_i/dxi_j (jacobian matrix) at each integration point of each precalculated dataset (ps_dme_dxdxi(i,j,ip,rule))
    real(kind=real64), allocatable :: ps_dme_dxidx(:,:,:,:)    !! dxi_i/dx_j (inverse of the jacobian matrix) at each integration point of each precalculated dataset (ps_dme_dxidx(i,j,ip,rule))
    real(kind=real64), allocatable :: ps_dme_dgphidx(:,:,:,:)  !! dphi^{dme-g}/dxj (inverse of the jacobian matrix) at each integration point of each precalculated dataset (ps_dme_dxidx(node,j,ip,rule))
    real(kind=real64), allocatable :: ps_dme_t(:,:,:,:)        !! Unit tangents at each integration point of each precalculated dataset ps_dme_t(kc,tang,ip,rule)
    real(kind=real64), allocatable :: ps_dme_j(:,:)            !! Geometric jacobian at each integration point of each precalculated dataset
  contains
    procedure, public :: init
    procedure, public :: init_precalculated_datasets
  end type fbem_bem_element

contains

  ! Initialization of fbem_bem_element
  subroutine init(e)
    implicit none
    class(fbem_bem_element) :: e
    e%gtype=0
    e%d=0
    e%n_gnodes=0
    e%n=0
    if (allocated(e%x)) deallocate(e%x)
    e%ptype=0
    e%ptype_delta=0.d0
    e%n_pnodes=0
    e%stype=0
    e%stype_delta=0
    e%n_snodes=0
    e%dmetype=0
    e%dme_d=0
    e%dme_n_gnodes=0
    if (allocated(e%dme_x)) deallocate(e%dme_x)
    e%dme_cl=0.d0
    e%cl=0.d0
    e%gln_far=0
    if (allocated(e%bball_centre)) deallocate(e%bball_centre)
    e%bball_radius=0.d0
    e%n_ps=0
    if (allocated(e%ps_gln)) deallocate(e%ps_gln)
    e%ps_gln_max=0
    if (allocated(e%ps_ngp)) deallocate(e%ps_ngp)
    e%ps_ngp_max=0
    if (allocated(e%ps_x)) deallocate(e%ps_x)
    if (allocated(e%ps_t)) deallocate(e%ps_t)
    if (allocated(e%ps_n)) deallocate(e%ps_n)
    if (allocated(e%ps_j)) deallocate(e%ps_j)
    if (allocated(e%ps_w)) deallocate(e%ps_w)
    if (allocated(e%ps_gphi)) deallocate(e%ps_gphi)
    if (allocated(e%ps_dgphidxi)) deallocate(e%ps_dgphidxi)
    if (allocated(e%ps_dxdxi)) deallocate(e%ps_dxdxi)
    if (allocated(e%ps_dxidx)) deallocate(e%ps_dxidx)
    if (allocated(e%ps_pphi)) deallocate(e%ps_pphi)
    if (allocated(e%ps_sphi)) deallocate(e%ps_sphi)
    if (allocated(e%ps_pphijw)) deallocate(e%ps_pphijw)
    if (allocated(e%ps_sphijw)) deallocate(e%ps_sphijw)
    if (allocated(e%ps_dme_gphi)) deallocate(e%ps_dme_gphi)
    if (allocated(e%ps_dme_dgphidxi)) deallocate(e%ps_dme_dgphidxi)
    if (allocated(e%ps_dme_dxdxi)) deallocate(e%ps_dme_dxdxi)
    if (allocated(e%ps_dme_dxidx)) deallocate(e%ps_dme_dxidx)
    if (allocated(e%ps_dme_dgphidx)) deallocate(e%ps_dme_dgphidx)
    if (allocated(e%ps_dme_t)) deallocate(e%ps_dme_t)
    if (allocated(e%ps_dme_j)) deallocate(e%ps_dme_j)
  end subroutine init

  !! Initialization of precalculated datasets for fbem_bem_element
  subroutine init_precalculated_datasets(e,n_ps,ps_gln)
    implicit none
    ! I/O
    class(fbem_bem_element) :: e            !! Element
    integer                 :: n_ps         !! Number of precalculated datasets
    integer                 :: ps_gln(n_ps) !! 1D Gauss-Legendre integ. points for each precalculated dataset
    ! Local
    integer           :: i, kphi, kxi, k1, k2, kt, order, kc
    real(kind=real64) :: aux(10)
    real(kind=real64) :: gphi(e%n_gnodes)
    real(kind=real64) :: dgphidxi(e%n_gnodes)
    real(kind=real64) :: dgphidxi1(e%n_gnodes)
    real(kind=real64) :: dgphidxi2(e%n_gnodes)
    real(kind=real64) :: pphi(e%n_pnodes)
    real(kind=real64) :: sphi(e%n_snodes)
    real(kind=real64) :: xi_p(e%d), wgp(e%d), j
    real(kind=real64) :: xp(e%n), T(e%n), T1(e%n), T2(e%n), N(e%n)
    real(kind=real64) :: dme_d
    real(kind=real64), allocatable :: dme_xi(:), dme_gphi(:)
    real(kind=real64), allocatable :: dme_dgphidxi(:), dme_dgphidxi1(:), dme_dgphidxi2(:)
    !
    ! Initialization
    !
    if (allocated(e%ps_gln)) deallocate(e%ps_gln)
    if (allocated(e%ps_ngp)) deallocate(e%ps_ngp)
    if (allocated(e%ps_x)) deallocate(e%ps_x)
    if (allocated(e%ps_t)) deallocate(e%ps_t)
    if (allocated(e%ps_n)) deallocate(e%ps_n)
    if (allocated(e%ps_j)) deallocate(e%ps_j)
    if (allocated(e%ps_w)) deallocate(e%ps_w)
    if (allocated(e%ps_gphi)) deallocate(e%ps_gphi)
    if (allocated(e%ps_dgphidxi)) deallocate(e%ps_dgphidxi)
    if (allocated(e%ps_dxdxi)) deallocate(e%ps_dxdxi)
    if (allocated(e%ps_dxidx)) deallocate(e%ps_dxidx)
    if (allocated(e%ps_pphi)) deallocate(e%ps_pphi)
    if (allocated(e%ps_sphi)) deallocate(e%ps_sphi)
    if (allocated(e%ps_pphijw)) deallocate(e%ps_pphijw)
    if (allocated(e%ps_sphijw)) deallocate(e%ps_sphijw)
    if (allocated(e%ps_dme_gphi)) deallocate(e%ps_dme_gphi)
    if (allocated(e%ps_dme_dgphidxi)) deallocate(e%ps_dme_dgphidxi)
    if (allocated(e%ps_dme_dxdxi)) deallocate(e%ps_dme_dxdxi)
    if (allocated(e%ps_dme_dxidx)) deallocate(e%ps_dme_dxidx)
    if (allocated(e%ps_dme_dgphidx)) deallocate(e%ps_dme_dgphidx)
    if (allocated(e%ps_dme_t)) deallocate(e%ps_dme_t)
    if (allocated(e%ps_dme_j)) deallocate(e%ps_dme_j)
    !
    ! Calculate the number of integration points for each rule
    !
    ! Initialize
    allocate (e%ps_gln(n_ps))
    allocate (e%ps_ngp(n_ps))
    e%n_ps=n_ps
    e%ps_gln=ps_gln
    select case (e%gtype)
      ! Line elements
      case (fbem_line2,fbem_line3,fbem_line4)
        do i=1,e%n_ps
          e%ps_ngp(i)=e%ps_gln(i)
        end do
      ! Quadrilateral elements
      case (fbem_quad4,fbem_quad8,fbem_quad9)
        do i=1,e%n_ps
          e%ps_ngp(i)=e%ps_gln(i)**2
        end do
      ! Triangular elements
      case (fbem_tri3,fbem_tri6)
        do i=1,e%n_ps
          ! If ngp>15 then Gauss-Legendre*Gauss-Jacobi
          if (ps_gln(i).gt.15) then
            e%ps_ngp(i)=e%ps_gln(i)**2
          ! If ngp<=15 Wandzura rules are used
          else
            ! Conversion between ngp and quadrature order
            order=2*e%ps_gln(i)-1
            ! Number of integration points
            e%ps_ngp(i)=wantri_n(order)
          end if
        end do
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'typeg={line2,line3,line4,tri3,tri6,quad4,quad8,quad9}')
    end select
    e%ps_gln_max=maxval(e%ps_gln)
    e%ps_ngp_max=maxval(e%ps_ngp)
    !
    ! Calculate the dataset
    !
    ! Initialization
    allocate(e%ps_x       (e%n             ,e%ps_ngp_max,e%ps_gln_max))
    allocate(e%ps_t       (e%n       ,e%n-1,e%ps_ngp_max,e%ps_gln_max))
    allocate(e%ps_n       (e%n             ,e%ps_ngp_max,e%ps_gln_max))
    allocate(e%ps_j       (                 e%ps_ngp_max,e%ps_gln_max))
    allocate(e%ps_w       (                 e%ps_ngp_max,e%ps_gln_max))
    allocate(e%ps_gphi    (e%n_gnodes,      e%ps_ngp_max,e%ps_gln_max))
    allocate(e%ps_dgphidxi(e%n_gnodes,  e%d,e%ps_ngp_max,e%ps_gln_max))
    allocate(e%ps_dxdxi   (e%n       ,  e%d,e%ps_ngp_max,e%ps_gln_max))
    allocate(e%ps_dxidx   (e%d       ,  e%n,e%ps_ngp_max,e%ps_gln_max))
    allocate(e%ps_pphi    (e%n_pnodes,      e%ps_ngp_max,e%ps_gln_max))
    allocate(e%ps_sphi    (e%n_snodes,      e%ps_ngp_max,e%ps_gln_max))
    allocate(e%ps_pphijw  (e%n_pnodes,      e%ps_ngp_max,e%ps_gln_max))
    allocate(e%ps_sphijw  (e%n_snodes,      e%ps_ngp_max,e%ps_gln_max))
    if (e%dmetype.ne.0) then
      allocate(e%ps_dme_gphi    (e%dme_n_gnodes        ,e%ps_ngp_max,e%ps_gln_max))
      allocate(e%ps_dme_dgphidxi(e%dme_n_gnodes,e%dme_d,e%ps_ngp_max,e%ps_gln_max))
      allocate(e%ps_dme_dxdxi   (e%n           ,e%dme_d,e%ps_ngp_max,e%ps_gln_max))
      allocate(e%ps_dme_dxidx   (e%dme_d       ,e%n    ,e%ps_ngp_max,e%ps_gln_max))
      allocate(e%ps_dme_dgphidx (e%dme_n_gnodes,e%n    ,e%ps_ngp_max,e%ps_gln_max))
      allocate(e%ps_dme_t       (           e%n,e%dme_d,e%ps_ngp_max,e%ps_gln_max))
      allocate(e%ps_dme_j       (                       e%ps_ngp_max,e%ps_gln_max))
      allocate (dme_xi(e%dme_d))
      dme_xi=0.d0
      allocate (dme_gphi(e%dme_n_gnodes))
      allocate (dme_dgphidxi(e%dme_n_gnodes),dme_dgphidxi1(e%dme_n_gnodes),dme_dgphidxi2(e%dme_n_gnodes))
    end if
    ! Calculate
    select case (e%gtype)
      !
      ! Line elements
      !
      case (fbem_line2,fbem_line3,fbem_line4)
       ! Loop through precalculated datasets
        do i=1,e%n_ps
          ! Loop through xi coordinate
          do kxi=1,gl11_n(e%ps_gln(i))
            ! Index of gaussian point
            kt=kxi
            ! xi coordinate and weight
            xi_p(1)=gl11_xi(kxi,e%ps_gln(i))
            wgp(1)=gl11_w(kxi,e%ps_gln(i))
            e%ps_w(kt,i)=wgp(1)
            ! Geometrical shape functions and first derivatives at xi
#           define etype e%gtype
#           define xi xi_p(1)
#           define delta 0.d0
#           define phi gphi
#           define dphidxi dgphidxi
#           include <phi_and_dphidxi_1d.rc>
#           undef etype
#           undef xi
#           undef phi
#           undef dphidxi
#           undef delta
            ! Functional shape functions (primary variables) at xi
#           define etype e%ptype
#           define xi xi_p(1)
#           define delta e%ptype_delta
#           define phi pphi
#           include <phi_1d.rc>
#           undef etype
#           undef xi
#           undef delta
#           undef phi
            ! Functional shape functions (secondary variables) at xi
#           define etype e%stype
#           define xi xi_p(1)
#           define delta e%stype_delta
#           define phi sphi
#           include <phi_1d.rc>
#           undef etype
#           undef xi
#           undef delta
#           undef phi
            e%ps_gphi(:,kt,i)=gphi
            e%ps_dgphidxi(:,1,kt,i)=dgphidxi
            e%ps_pphi(:,kt,i)=pphi
            e%ps_sphi(:,kt,i)=sphi
            ! Components calculation of x and T at xi
            xp=0.d0
            T=0.d0
            do kphi=1,e%n_gnodes
              xp=xp+gphi(kphi)*e%x(:,kphi)
              T=T+dgphidxi(kphi)*e%x(:,kphi)
            end do
            ! Position vector
            e%ps_x(:,kt,i)=xp
            ! Tangent vector (jacobian matrix)
            e%ps_dxdxi(:,1,kt,i)=T
            ! Geometric jacobian
            j=sqrt(dot_product(T,T))
            e%ps_j(kt,i)=j
            ! Inverse of the jacobian matrix
            e%ps_dxidx(1,:,kt,i)=T/j**2
            ! Unit tangent
            e%ps_t(:,1,kt,i)=T/j
            ! Unit normal vector
            select case (e%n)
              ! 2D: the unit normal is saved in n
              case (2)
                e%ps_n(1,kt,i)= e%ps_t(2,1,kt,i)
                e%ps_n(2,kt,i)=-e%ps_t(1,1,kt,i)
              ! 3D: the unit normal is undefined (set to 0.)
              case (3)
                e%ps_n(:,kt,i)=0.d0
            end select
            ! Functional shape functions * jacobian * weight
            e%ps_pphijw(:,kt,i)=pphi(:)*j*wgp(1)
            e%ps_sphijw(:,kt,i)=sphi(:)*j*wgp(1)
            ! Design macro-element in which the element is embedded
            if (e%dmetype.ne.0) then
              ! Calculate the local coordinates of xp in the design macro-element
              call fbem_local_coordinates(e%n,e%dmetype,e%dme_x,e%dme_cl,xp,dme_xi,dme_d)
              if (dme_d.gt.1.d-12) then
                write(error_unit,*) 'd=',dme_d
                write(error_unit,*) 'dme_x=',e%dme_x
                write(error_unit,*) 'x=',xp
                write(error_unit,*) 'xi=',dme_xi
                stop 'the integration point is not in the design macro-element'
              end if
              ! Shape functions and first derivatives at dme_xi
              select case (e%dme_d)
                ! 1D DESIGN MACRO-ELEMENT
                case (1)
#                 define etype e%dmetype
#                 define xi dme_xi(1)
#                 define delta 0.d0
#                 define phi dme_gphi
#                 define dphidxi dme_dgphidxi
#                 include <phi_and_dphidxi_1d.rc>
#                 undef etype
#                 undef xi
#                 undef phi
#                 undef dphidxi
#                 undef delta
                  e%ps_dme_gphi(:,kt,i)=dme_gphi
                  e%ps_dme_dgphidxi(:,1,kt,i)=dme_dgphidxi
                  ! Calculation of T
                  T=0.d0
                  do kphi=1,e%dme_n_gnodes
                    T=T+dme_dgphidxi(kphi)*e%dme_x(:,kphi)
                  end do
                  ! Jacobian matrix
                  e%ps_dme_dxdxi(:,1,kt,i)=T
                  ! Geometric jacobian
                  j=sqrt(dot_product(T,T))
                  e%ps_dme_j(kt,i)=j
                  ! Unit tangent
                  e%ps_dme_t(:,1,kt,i)=T/j
                  ! Inverse of the jacobian matrix
                  e%ps_dme_dxidx(1,:,kt,i)=T/j**2
                  ! Shape functions derivatives with respect to cartesian coordinates
                  do kc=1,e%n
                    e%ps_dme_dgphidx(:,kc,kt,i)=T(kc)*dme_dgphidxi(:)
                  end do
                  e%ps_dme_dgphidx(:,:,kt,i)=e%ps_dme_dgphidx(:,:,kt,i)/j**2
                ! 2D DESIGN MACRO-ELEMENT
                case (2)
#                 define etype e%dmetype
#                 define xi dme_xi
#                 define delta 0.d0
#                 define phi dme_gphi
#                 define dphidxi1 dme_dgphidxi1
#                 define dphidxi2 dme_dgphidxi2
#                 include <phi_and_dphidxik_2d.rc>
#                 undef etype
#                 undef xi
#                 undef phi
#                 undef dphidxi1
#                 undef dphidxi2
#                 undef delta
                  e%ps_dme_gphi(:,kt,i)=dme_gphi
                  e%ps_dme_dgphidxi(:,1,kt,i)=dme_dgphidxi1
                  e%ps_dme_dgphidxi(:,2,kt,i)=dme_dgphidxi2
                  ! Calculation of T
                  T1=0.d0
                  T2=0.d0
                  do kphi=1,e%dme_n_gnodes
                    T1=T1+dme_dgphidxi1(kphi)*e%dme_x(:,kphi)
                    T2=T2+dme_dgphidxi2(kphi)*e%dme_x(:,kphi)
                  end do
                  ! Jacobian matrix
                  e%ps_dme_dxdxi(:,1,kt,i)=T1
                  e%ps_dme_dxdxi(:,2,kt,i)=T2
                  ! Unit tangents
                  e%ps_dme_t(:,1,kt,i)=T1/sqrt(dot_product(T1,T1))
                  e%ps_dme_t(:,2,kt,i)=T2/sqrt(dot_product(T2,T2))
                  ! Inverse of the jacobian matrix and jacobian determinant
                  select case (e%n)
                    case (2)
                      call fbem_invert_2x2_matrix(e%ps_dme_dxdxi(:,:,kt,i),e%ps_dme_dxidx(:,:,kt,i),j)
                      e%ps_dme_j(kt,i)=j
                      e%ps_dme_dgphidx(:,1,kt,i)=e%ps_dme_dxidx(1,1,kt,i)*dme_dgphidxi1(:)+e%ps_dme_dxidx(2,1,kt,i)*dme_dgphidxi2(:)
                      e%ps_dme_dgphidx(:,2,kt,i)=e%ps_dme_dxidx(1,2,kt,i)*dme_dgphidxi1(:)+e%ps_dme_dxidx(2,2,kt,i)*dme_dgphidxi2(:)
                    case (3)
                      stop 'not yet init_precalculated_datasets'
                  end select
                case default
                  stop 'Invalid case: e%dme_d'
              end select
            end if
          end do ! Loop through xi coordinate
        end do ! Loop through rules
      !
      ! Quadrilateral elements
      !
      case (fbem_quad4,fbem_quad8,fbem_quad9)
        ! Loop through precalculated datasets
        do i=1,e%n_ps
          ! Loop through xi_1 direction
          do k1=1,gl11_n(e%ps_gln(i))
            ! xi_1 coordinate and weight
            xi_p(1)=gl11_xi(k1,e%ps_gln(i))
            wgp(1)=gl11_w(k1,e%ps_gln(i))
            ! Loop through xi_2 direction
            do k2=1,gl11_n(e%ps_gln(i))
              ! Total index of gaussian point
              kt=k2+(k1-1)*gl11_n(e%ps_gln(i))
              ! xi_2 coordinate and weight
              xi_p(2)=gl11_xi(k2,e%ps_gln(i))
              wgp(2)=gl11_w(k2,e%ps_gln(i))
              e%ps_w(kt,i)=wgp(1)*wgp(2)
              ! Geometrical shape functions and first derivatives at xi
#             define etype e%gtype
#             define xi xi_p
#             define delta 0.d0
#             define phi gphi
#             define dphidxi1 dgphidxi1
#             define dphidxi2 dgphidxi2
#             include <phi_and_dphidxik_2d.rc>
#             undef etype
#             undef xi
#             undef phi
#             undef dphidxi1
#             undef dphidxi2
#             undef delta
              ! Functional shape functions (primary variables) at xi
#             define etype e%ptype
#             define xi xi_p
#             define delta e%ptype_delta
#             define phi pphi
#             include <phi_2d.rc>
#             undef etype
#             undef xi
#             undef delta
#             undef phi
              ! Functional shape functions (secondary variables) at xi
#             define etype e%stype
#             define xi xi_p
#             define delta e%stype_delta
#             define phi sphi
#             include <phi_2d.rc>
#             undef etype
#             undef xi
#             undef delta
#             undef phi
              e%ps_gphi(:,kt,i)=gphi
              e%ps_dgphidxi(:,1,kt,i)=dgphidxi1
              e%ps_dgphidxi(:,2,kt,i)=dgphidxi2
              e%ps_pphi(:,kt,i)=pphi
              e%ps_sphi(:,kt,i)=sphi
              ! Components calculation of x, T1 and T2 at xi
              xp=0.d0
              T1=0.d0
              T2=0.d0
              do kphi=1,e%n_gnodes
                xp=xp+gphi(kphi)*e%x(:,kphi)
                T1=T1+dgphidxi1(kphi)*e%x(:,kphi)
                T2=T2+dgphidxi2(kphi)*e%x(:,kphi)
              end do
              ! Position vector
              e%ps_x(:,kt,i)=xp
              ! Unit tangents
              e%ps_t(:,1,kt,i)=T1/sqrt(dot_product(T1,T1))
              e%ps_t(:,2,kt,i)=T2/sqrt(dot_product(T2,T2))
              ! Jacobian and unit normal vector
              select case(e%n)
                ! 2D: the jacobian is the classic jacobian, but the normal is in 3D space, so it is undefined
                case (2)
                  ! Normal (undefined, set to 0.)
                  e%ps_n(:,kt,i)=0.d0
                  ! Jacobian
                  j=T1(1)*T2(2)-T1(2)*T2(1)
                ! 3D: the jacobian is |T1 x T2|, which is a non-unit normal
                case (3)
                  ! Normal vector as T1 x T2 at xi
                  N(1)=T1(2)*T2(3)-T1(3)*T2(2)
                  N(2)=T1(3)*T2(1)-T1(1)*T2(3)
                  N(3)=T1(1)*T2(2)-T1(2)*T2(1)
                  ! Jacobian
                  j=sqrt(dot_product(N,N))
                  ! Unit normal
                  e%ps_n(:,kt,i)=N/j
              end select
              e%ps_j(kt,i)=j
              ! Functional shape functions * jacobian * weights
              e%ps_pphijw(:,kt,i)=pphi(:)*j*wgp(1)*wgp(2)
              e%ps_sphijw(:,kt,i)=sphi(:)*j*wgp(1)*wgp(2)
              ! Design macro-element in which the element is embedded
              if (e%dmetype.ne.0) then
                stop 'dme not yet for 2D BE'
              end if
            end do ! Loop through xi_2 direction
          end do ! Loop through xi_1 direction
        end do ! Loop through rules
      !
      ! Triangular elements
      !
      case (fbem_tri3,fbem_tri6)
        ! Loop through precalculated datasets
        do i=1,e%n_ps
          ! If ngp>15 then Gauss-Legendre*Gauss-Jacobi
          if (e%ps_gln(i).gt.15) then
            ! Loop through rule points
            do k1=1,gl01_n(e%ps_gln(i))
              do k2=1,gj01_n(e%ps_gln(i))
                ! Total index of gaussian point
                kt=k2+(k1-1)*gl01_n(e%ps_gln(i))
                ! xi_1 coordinate and xi_2 coordinate and weights
                xi_p(1)=(1.0d0-gj01_xi(k2,e%ps_gln(i)))*gl01_xi(k1,e%ps_gln(i))
                xi_p(2)=gj01_xi(k2,e%ps_gln(i))
                wgp(1)=gl01_w(k1,e%ps_gln(i))
                wgp(2)=gj01_w(k2,e%ps_gln(i))
                e%ps_w(kt,i)=wgp(1)*wgp(2)
                ! Geometrical shape functions and first derivatives at xi
#               define etype e%gtype
#               define xi xi_p
#               define delta 0.d0
#               define phi gphi
#               define dphidxi1 dgphidxi1
#               define dphidxi2 dgphidxi2
#               include <phi_and_dphidxik_2d.rc>
#               undef etype
#               undef xi
#               undef delta
#               undef phi
#               undef dphidxi1
#               undef dphidxi2
                ! Functional shape functions (primary variables) at xi
#               define etype e%ptype
#               define xi xi_p
#               define delta e%ptype_delta
#               define phi pphi
#               include <phi_2d.rc>
#               undef etype
#               undef xi
#               undef delta
#               undef phi
                ! Functional shape functions (secondary variables) at xi
#               define etype e%stype
#               define xi xi_p
#               define delta e%stype_delta
#               define phi sphi
#               include <phi_2d.rc>
#               undef etype
#               undef xi
#               undef delta
#               undef phi
                e%ps_gphi(:,kt,i)=gphi
                e%ps_dgphidxi(:,1,kt,i)=dgphidxi1
                e%ps_dgphidxi(:,2,kt,i)=dgphidxi2
                e%ps_pphi(:,kt,i)=pphi
                e%ps_sphi(:,kt,i)=sphi
                ! Components calculation of x, T1 and T2 at xi
                xp=0.d0
                T1=0.d0
                T2=0.d0
                do kphi=1,e%n_gnodes
                  xp=xp+gphi(kphi)*e%x(:,kphi)
                  T1=T1+dgphidxi1(kphi)*e%x(:,kphi)
                  T2=T2+dgphidxi2(kphi)*e%x(:,kphi)
                end do
                ! Position vector
                e%ps_x(:,kt,i)=xp
                ! Unit tangents
                e%ps_t(:,1,kt,i)=T1/sqrt(dot_product(T1,T1))
                e%ps_t(:,2,kt,i)=T2/sqrt(dot_product(T2,T2))
                ! Jacobian and unit normal vector
                select case(e%n)
                  ! 2D: the jacobian is the classic jacobian, but the normal is in 3D space, so it is undefined
                  case (2)
                    ! Normal (undefined, set to 0.)
                    e%ps_n(:,kt,i)=0.d0
                    ! Jacobian
                    j=T1(1)*T2(2)-T1(2)*T2(1)
                  ! 3D: the jacobian is |T1 x T2|, which is a non-unit normal
                  case (3)
                    ! Normal vector as T1 x T2 at xi
                    N(1)=T1(2)*T2(3)-T1(3)*T2(2)
                    N(2)=T1(3)*T2(1)-T1(1)*T2(3)
                    N(3)=T1(1)*T2(2)-T1(2)*T2(1)
                    ! Jacobian
                    j=sqrt(dot_product(N,N))
                    ! Unit normal
                    e%ps_n(:,kt,i)=N/j
                end select
                e%ps_j(kt,i)=j
                ! Functional shape functions * jacobian * weights
                e%ps_pphijw(:,kt,i)=pphi(:)*j*wgp(1)*wgp(2)
                e%ps_sphijw(:,kt,i)=sphi(:)*j*wgp(1)*wgp(2)
                if (e%dmetype.ne.0) then
                  stop 'dme not yet for 2D BE'
                end if
              end do
            end do
          ! If ngp<=15 Wandzura rules are used
          else
            ! Conversion between ngp and quadrature order
            order=2*e%ps_gln(i)-1
            ! Loop through rule points
            do kt=1,wantri_n(order)
              ! xi_1 coordinate
              xi_p(1)=wantri_xi1(kt,order)
              ! xi_2 coordinate
              xi_p(2)=wantri_xi2(kt,order)
              ! Weight
              wgp(1)=wantri_w(kt,order)
              e%ps_w(kt,i)=wgp(1)
              ! Geometrical shape functions and first derivatives at xi
#             define etype e%gtype
#             define xi xi_p
#             define delta 0.d0
#             define phi gphi
#             define dphidxi1 dgphidxi1
#             define dphidxi2 dgphidxi2
#             include <phi_and_dphidxik_2d.rc>
#             undef etype
#             undef xi
#             undef delta
#             undef phi
#             undef dphidxi1
#             undef dphidxi2
              ! Functional shape functions (primary variables) at xi
#             define etype e%ptype
#             define xi xi_p
#             define delta e%ptype_delta
#             define phi pphi
#             include <phi_2d.rc>
#             undef etype
#             undef xi
#             undef delta
#             undef phi
              ! Functional shape functions (secondary variables) at xi
#             define etype e%stype
#             define xi xi_p
#             define delta e%stype_delta
#             define phi sphi
#             include <phi_2d.rc>
#             undef etype
#             undef xi
#             undef delta
#             undef phi
              e%ps_gphi(:,kt,i)=gphi
              e%ps_dgphidxi(:,1,kt,i)=dgphidxi1
              e%ps_dgphidxi(:,2,kt,i)=dgphidxi2
              e%ps_pphi(:,kt,i)=pphi
              e%ps_sphi(:,kt,i)=sphi
              ! Components calculation of x, T1 and T2 at xi
              xp=0.d0
              T1=0.d0
              T2=0.d0
              do kphi=1,e%n_gnodes
                xp=xp+gphi(kphi)*e%x(:,kphi)
                T1=T1+dgphidxi1(kphi)*e%x(:,kphi)
                T2=T2+dgphidxi2(kphi)*e%x(:,kphi)
              end do
              ! Position vector
              e%ps_x(:,kt,i)=xp
              ! Unit tangents
              e%ps_t(:,1,kt,i)=T1/sqrt(dot_product(T1,T1))
              e%ps_t(:,2,kt,i)=T2/sqrt(dot_product(T2,T2))
              ! Jacobian and unit normal vector
              select case(e%n)
                ! 2D: the jacobian is the classic jacobian, but the normal is in 3D space, so it is undefined
                case (2)
                  ! Normal (undefined, set to 0.)
                  e%ps_n(:,kt,i)=0.d0
                  ! Jacobian
                  j=T1(1)*T2(2)-T1(2)*T2(1)
                ! 3D: the jacobian is |T1 x T2|, which is a non-unit normal
                case (3)
                  ! Normal vector as T1 x T2 at xi
                  N(1)=T1(2)*T2(3)-T1(3)*T2(2)
                  N(2)=T1(3)*T2(1)-T1(1)*T2(3)
                  N(3)=T1(1)*T2(2)-T1(2)*T2(1)
                  ! Jacobian
                  j=sqrt(dot_product(N,N))
                  ! Unit normal
                  e%ps_n(:,kt,i)=N/j
              end select
              e%ps_j(kt,i)=j
              ! Functional shape functions * jacobian * weights
              e%ps_pphijw(:,kt,i)=pphi(:)*j*wgp(1)
              e%ps_sphijw(:,kt,i)=sphi(:)*j*wgp(1)
              if (e%dmetype.ne.0) then
                stop 'dme not yet for 2D BE'
              end if
            end do
          end if
        end do
      case default
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'etype={line2,line3,line4,tri3,tri6,quad4,quad8,quad9}')
    end select
  end subroutine init_precalculated_datasets

  !! This subroutine calculates the integrals log(1/r)*phi^f over a line when singular (weakly singular integral).
  subroutine fbem_bem_logd1rphi_1d_int(rn,e_gtype,e_x,e_ftype,e_ftype_delta,xi_i,lws)
    implicit none
    ! I/O
    integer           :: rn                            !! Spatial dimension
    integer           :: e_gtype                       !! Geometrical interpolation
    real(kind=real64) :: e_x(rn,fbem_n_nodes(e_gtype)) !! Coordinates of the nodes
    integer           :: e_ftype                       !! Functional interpolation (phi)
    real(kind=real64) :: e_ftype_delta                 !! Functional interpolation displacement towards inside
    real(kind=real64) :: xi_i(1)                       !! Local coordinate of the singular point.
    real(kind=real64) :: lws(fbem_n_nodes(e_ftype))    !! log(1/r)*phi^f kernel
    ! Local
    integer                      :: gln                             ! 1D Gauss-Legendre number of integration points (<=32)
    integer                      :: kphi                            ! Counter variable for shape functions loops
    integer                      :: kip                             ! Counter of integration points
    real(kind=real64)            :: x_i(2)                          ! Real coordinates of collocation point
    real(kind=real64)            :: gphi(fbem_n_nodes(e_gtype))     ! Geometrical shape functions values at xi
    real(kind=real64)            :: dgphidxi(fbem_n_nodes(e_gtype)) ! Geometrical shape functions derivatives values at xi
    real(kind=real64)            :: fphi(fbem_n_nodes(e_ftype))     ! Functional shape functions values at xi
    integer                      :: nsub                            ! Number of subdivision of the element
    integer                      :: ksub                            ! Counter of subdivision
    real(kind=real64)            :: gamma                           ! Coordinate gamma
    real(kind=real64)            :: w                               ! Weights of each integration point
    real(kind=real64)            :: xip                             ! Coordinate  xip of subdivided element [0,1]
    real(kind=real64)            :: js                              ! Jacobian of the xi [xisub(1,:),xisub[2,:]] -> xip [0,1] transformation
    real(kind=real64)            :: xip_i(2)                        ! Singular point in xip space
    real(kind=real64)            :: xi                              ! Coordinate xi
    real(kind=real64)            :: xisub(2,2)                      ! Coordinates of element subdivisions
    real(kind=real64)            :: aux(10)                         ! Auxiliary variable needed for shape_functions module resources
    real(kind=real64)            :: x(2)                            ! Position vector at xi
    real(kind=real64)            :: T(2)                            ! Tangent vector at xi
    real(kind=real64)            :: rv(2)                           ! Distance vector between collocation point and integration point (x-x_i)
    real(kind=real64)            :: r, d1r, logd1r                  ! Distance vector module, its inverse and log(1/r)
    real(kind=real64)            :: jg                              ! Geometric jacobian
    real(kind=real64)            :: jw                              ! Jacobians * weight
    type(fbem_telles_parameters) :: telles_parameters               ! Telles parameters
    real(kind=real64)            :: jt                              ! Telles jacobian
    ! Integration points to be used
    gln=30
    ! Initialize kernel
    lws=0.d0
    ! Setup the subdivisions for transformation xip -> xi: [0,1] -> [xisub_1,xisub_2]
    if (fbem_check_xi_vertex(xi_i(1))) then
      nsub=1
      xisub(1,1)=-1.0d0
      xisub(2,1)= 1.0d0
      if (xi_i(1).lt.0.0d0) xip_i(1)=0.0d0
      if (xi_i(1).gt.0.0d0) xip_i(1)=1.0d0
    else
      nsub=2
      xisub(1,1)=-1.0d0
      xisub(2,1)=xi_i(1)
      xip_i(1)=1.0d0
      xisub(1,2)=xi_i(1)
      xisub(2,2)=1.0d0
      xip_i(2)=0.0d0
    end if
    ! Calculate spatial coordinates at xi_i
#   define etype e_gtype
#   define delta 0.d0
#   define xi xi_i(1)
#   define phi gphi
#   include <phi_1d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
    x_i=0.d0
    do kphi=1,fbem_n_nodes(e_gtype)
      x_i=x_i+gphi(kphi)*e_x(:,kphi)
    end do
    ! Loop through subdivisions
    do ksub=1,nsub
      ! Telles transformation (gamma [0,1] -> xip [0,1]) setup (null jacobian at the collocation point)
      telles_parameters=fbem_telles01_calculate_parameters(xip_i(ksub),0.0d0)
      ! Add integration points
      do kip=1,gl01_n(gln)
        ! GAMMA COORDINATE
        gamma=gl01_xi(kip,gln)
        w=gl01_w(kip,gln)
        ! GAMMA->XIP TRANSFORMATION
        call fbem_telles_xi_and_jacobian(telles_parameters,gamma,xip,jt)
        ! XIP-> XI TRANSFORMATION
        js=xisub(2,ksub)-xisub(1,ksub)
        xi=js*xip+xisub(1,ksub)
        ! XI->X TRANSFORMATION
#       define etype e_gtype
#       define delta 0.d0
#       define phi gphi
#       define dphidxi dgphidxi
#       include <phi_and_dphidxi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
#       undef dphidxi
        x=0.d0
        T=0.d0
        do kphi=1,fbem_n_nodes(e_gtype)
          x=x+gphi(kphi)*e_x(:,kphi)
          T=T+dgphidxi(kphi)*e_x(:,kphi)
        end do
        jg=sqrt(dot_product(T,T))
        ! Distance vector and other derived terms
        rv=x-x_i
        r=sqrt(dot_product(rv,rv))
        d1r=1.d0/r
        logd1r=log(d1r)
        ! Jacobians * weight
        jw=jg*js*jt*w
        ! FUNCTIONAL SHAPE FUNCTIONS
        ! Functional shape functions (secondary variables) at xi
#       define etype e_ftype
#       define delta e_ftype_delta
#       define phi fphi
#       include <phi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
        ! Add integration points
        lws(:)=lws(:)+logd1r*fphi*jw
      end do
    end do
  end subroutine fbem_bem_logd1rphi_1d_int

  !! Calculation of the integral log(r) over a line when it is singular (weakly singular integral).
  function fbem_bem_int_logr(rn,e_gtype,e_x,xi_i)
    implicit none
    ! I/O
    real(kind=real64) :: fbem_bem_int_logr             !! log(1/r)*phi^f kernel
    integer           :: rn                            !! Spatial dimension
    integer           :: e_gtype                       !! Geometrical interpolation
    real(kind=real64) :: e_x(rn,fbem_n_nodes(e_gtype)) !! Position vector of the nodes
    real(kind=real64) :: xi_i(1)                       !! Local coordinate of the singular point.
    ! Local
    integer                      :: gln                             ! Gauss-Legendre number of integration points (<=32)
    integer                      :: kphi                            ! Counter variable for shape functions loops
    integer                      :: kip                             ! Counter of integration points
    real(kind=real64)            :: x_i(2)                          ! Real coordinates of collocation point
    real(kind=real64)            :: gphi(fbem_n_nodes(e_gtype))     ! Geometrical shape functions values at xi
    real(kind=real64)            :: dgphidxi(fbem_n_nodes(e_gtype)) ! Geometrical shape functions derivatives values at xi
    integer                      :: nsub                            ! Number of subdivision of the element
    integer                      :: ksub                            ! Counter of subdivision
    real(kind=real64)            :: gamma                           ! Coordinate gamma
    real(kind=real64)            :: w                               ! Weights of each integration point
    type(fbem_telles_parameters) :: telles_parameters               ! Telles parameters
    real(kind=real64)            :: jt                              ! Telles jacobian
    real(kind=real64)            :: xip                             ! Coordinate  xip of subdivided element [0,1]
    real(kind=real64)            :: js                              ! Jacobian of the xi [xisub(1,:),xisub[2,:]] -> xip [0,1] transformation
    real(kind=real64)            :: xip_i(2)                        ! Singular point in xip space
    real(kind=real64)            :: xi                              ! xi coordinate (element space [-1,1])
    real(kind=real64)            :: xisub(2,2)                      ! Coordinates of element subdivisions
    real(kind=real64)            :: aux(10)                         ! Auxiliary variable
    real(kind=real64)            :: x(2), T(2)                      ! Position and tangent vectors
    real(kind=real64)            :: jg                              ! Geometric jacobian: xi [-1,1] -> x
    real(kind=real64)            :: rv(2), r, logr                  ! Distance vector, and its modulus and natural logarithm
    ! Initialization
    gln=30
    fbem_bem_int_logr=0.d0
    ! Setup the subdivisions for transformation xip -> xi: [0,1] -> [xisub_1,xisub_2]
    if (fbem_check_xi_vertex(xi_i(1))) then
      nsub=1
      xisub(1,1)=-1.0d0
      xisub(2,1)= 1.0d0
      if (xi_i(1).lt.0.0d0) xip_i(1)=0.0d0
      if (xi_i(1).gt.0.0d0) xip_i(1)=1.0d0
    else
      nsub=2
      xisub(1,1)=-1.0d0
      xisub(2,1)=xi_i(1)
      xip_i(1)=1.0d0
      xisub(1,2)=xi_i(1)
      xisub(2,2)=1.0d0
      xip_i(2)=0.0d0
    end if
    ! Calculate position vector at xi_i
#   define etype e_gtype
#   define delta 0.d0
#   define xi xi_i(1)
#   define phi gphi
#   include <phi_1d.rc>
#   undef etype
#   undef delta
#   undef xi
#   undef phi
    x_i=0.d0
    do kphi=1,fbem_n_nodes(e_gtype)
      x_i=x_i+gphi(kphi)*e_x(:,kphi)
    end do
    ! Loop through subdivisions
    do ksub=1,nsub
      ! Telles transformation (gamma [0,1] -> xip [0,1]) setup (null jacobian at the collocation point)
      telles_parameters=fbem_telles01_calculate_parameters(xip_i(ksub),0.0d0)
      ! Add integration points
      do kip=1,gl01_n(gln)
        ! gamma [0,1]
        gamma=gl01_xi(kip,gln)
        w=gl01_w(kip,gln)
        ! gamma [0,1] -> xip [0,1]
        call fbem_telles_xi_and_jacobian(telles_parameters,gamma,xip,jt)
        ! xip [0,1] -> xi [-1,1]
        js=xisub(2,ksub)-xisub(1,ksub)
        xi=js*xip+xisub(1,ksub)
        ! xi [-1,1] -> x, T
#       define etype e_gtype
#       define delta 0.d0
#       define phi gphi
#       define dphidxi dgphidxi
#       include <phi_and_dphidxi_1d.rc>
#       undef etype
#       undef delta
#       undef phi
#       undef dphidxi
        ! Position vector x and tangent vector T
        x=0.d0
        T=0.d0
        do kphi=1,fbem_n_nodes(e_gtype)
          x=x+gphi(kphi)*e_x(:,kphi)
          T=T+dgphidxi(kphi)*e_x(:,kphi)
        end do
        ! Jacobian: xi [-1,1] -> x
        jg=sqrt(dot_product(T,T))
        ! Distance vector, and its modulus and natural logarithm
        rv=x-x_i
        r=sqrt(dot_product(rv,rv))
        logr=log(r)
        ! Add integration points
        fbem_bem_int_logr=fbem_bem_int_logr+logr*jg*js*jt*w
      end do
    end do
  end function fbem_bem_int_logr

end module fbem_bem_general
