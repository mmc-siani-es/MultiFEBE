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
! is displayed. If CHECK_XIS == 2, checking is enabled, an error message is displayed and execution stops.Other values of disable
! any checking.
#define CHECK_XIS 0

!! @author Jacob David Rodriguez Bordon (jacobdavid.rodriguezbordon@ulpgc.es)
!!
!! @version 2.0
!!
!! <b> This module implements the fundamental part of a methodology that estimates the number of Gauss-Legendre points needed for
!! quasi-singular integration in BEM. See <tt>/doc/complementary/quasi-singular/</tt> folder. </b>
!!
!! <h2>DESCRIPTION</h2>
!!
!! The methodology has the following characteristics:
!!  -It uses Gauss-Legendre quadrature with or without Telles transformation applied from <tt>N=2</tt> to <tt>N=30</tt>
!!   (order <tt>O=2*2-1=3</tt> to <tt>O=2*30-1=59</tt>). When using Telles transformation the <tt>jacobian(d)</tt> relationship
!!   is <tt>fbem_telles_barr(d,fbem_f_any)</tt>.
!!  -Geometric interpolation elements are: 1D (<tt>fbem_line2,fbem_line3</tt>) and 2D (<tt>fbem_tri3,fbem_tri6,fbem_quad4,
!!   fbem_quad8,fbem_quad9</tt>).
!!  -Quasi-singular terms <tt>f</tt> are: <tt>ln(r)</tt>, <tt>1/r</tt>, <tt>1/r^2</tt> and <tt>1/r^3</tt>.
!!  -Studied error heights are: <tt>1.0E-01</tt>, <tt>1.0E-02</tt>, ..., <tt>1.0E-11</tt> and <tt>1.0E-12</tt>.
!!  -It is needed an external subroutine that calculates nearest element point (in reference space) to singular point, and the
!!   corresponding normalized distance.
!!  -In order to be conservative, the normalized distance to be used should be (<tt>real_distance/max(real_element_width))</tt>),
!!   and the number of Gauss points in each coordinate direction for 2D elements should be the same. A more adaptative alternative
!!   is use different normalized distances for each coordinate direction, and calculate number of Gauss points required for each
!!   normalized distance.
!!  -When normalized distance between singular point and integration domain is too small, the number of Gauss points tends to
!!   infinity in same cases. Because the study is truncated at <tt>N=30</tt>, a small distance limit exists for studied error
!!   heights. The routines warn about these situations.
!!  -It takes in account element distortion and curvature by using the maximum between number of Gauss points needed by <tt>f</tt>
!!   and number of Gauss points needed by <tt>phi*jacobian</tt>.
!!
!! @see fbem_telles_transformation

module fbem_quasisingular_integration

  ! Fortran 2003 standard
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_quad_rules

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all are private
  private

  ! Public functions
  public :: fbem_qs_parameters
  public :: fbem_qs_calculate_parameters
  public :: fbem_qs_n_estimation_standard
  public :: fbem_qs_n_estimation_telles

  public :: fbem_qs_phijac_ngp_1d
  public :: fbem_qs_phijac_ngp_2d

!  public :: fbem_qs_test_standard_line
!  public :: fbem_qs_test_disttrans_line
!  public :: fbem_qs_test_standard_quad
!  public :: fbem_qs_test_standard_tri

  ! ================================================================================================================================
  ! COEFFICIENTS OF A CONTINUOUS APPROXIMATION OF THE NUMBER OF GAUSS-LEGENDRE POINTS (VALID FOR LINE, QUADS AND TRIS ELEMENTS)
  ! ================================================================================================================================
  !
  ! The estimation is built by using a proper combination of approximated experimental curves N(d) along three characteristic paths
  ! around a line element:
  !   - Curve 1, N1(d): path 1 (center / perpendicular path), very demanding for small d.
  !   - Curve 2, N2(d): path 2 (vertex / perpendicular path), less demanding than curve 1 for small d, and similar for large d.
  !   - Curve 3, N3(d): path 3 (vertex / coplanar path     ), less demanding than curve 2 for small d, but a little bit more
  !     demanding for large d.
  !
  ! The approximation of each particular curve is of the type:
  !   N(d)=10**(alpha0+alpha1*log10(d))/(1+beta1*log10(d)), d=[dmin,dmax]
  !
  ! Eight different integrands are considered ln(r), 1/r, 1/r^2, ..., 1/r^7.
  !
  ! d=r/L is the dimensionless distance, where r is the real minimum distance between collocation point and element, and L is the
  ! length of the element (for 2D elements is the longest edge length).
  !
  ! LINE ELEMENT
  !
  ! When d<=2, it is required a precise value of d, the position of the element nearest to the collocation point, and the relative
  ! orientation between collocation point and element. The path/curve pairs should be:
  !   - Path 1: curve 1
  !   - Path 2: curve 2
  !   - Path 3: curve 3
  !
  !       2          1          2
  !       2          1          2
  !       2          1          2
  ! 3-3-3 v==========c==========v 3-3-3
  !       2          1          2
  !       2          1          2
  !       2          1          2
  !
  ! When d>2, only matters the distance d (an approximate value is enough), and the curve 3 is used as N(d).
  !
  ! QUADRILATERAL ELEMENT (TRIANGULAR ELEMENT IS ANALOGOUS)
  !
  ! Five characteristic paths can be defined:
  !   - Path 1: center      / perpendicular path
  !   - Path 2: edge center / perpendicular path
  !   - Path 3: vertex      / perpendicular path
  !   - Path 4: edge center / coplanar path
  !   - Path 5: vertex      / coplanar path
  !
  ! Assignment of curves built for line element applied to 2D elements:
  !
  ! When d<=2, it is required a precise value of d and the position of the element nearest to the collocation point. The path/curve
  ! pair should be:
  !   - Path 1: curve 1
  !   - Path 2: curve 1
  !   - Path 3: curve 2
  !   - Path 4: curve 1
  !   - Path 5: curve 2
  !
  !             2 2        1 1        2 2
  !             22         11         22
  !       2-2-2 v==========c==========v 2-2-2
  !          1 /2       1 /1       1//2
  !          1//2       1/ 1       1/ 2
  !    1-1-1 c----------c----------c 1-1-1
  !       2 /1       1 /1       2//1
  !       2//1       1/ 1       2/ 1
  ! 2-2-2 v==========c==========v 2-2-2
  !      22         11         22
  !     2 2        1 1        2 2
  !
  ! When d>2, only matters the distance d (an approximate value is enough), and the curve 3 is used as N(d).
  !
  ! General notes:
  !   - When d<=2, characteristic path 1 is the most restrictive, while when d>2 the most restrictive is characteristic path 3. For
  !     d>2, the difference increases for stronger integrands, being their curves almost similar for 1/r, and having a difference
  !     <=1 integration point for 1/r^7. Therefore, for d>2 the curve associated with the characteristic path 3.
  !   - For Standard integration, when d<=2, the difference between characteristic paths 2 and 3 is small, these differences
  !     decreases as the degree of singularity increase. In the case of Telles integration, there exist considerable differences
  !     between both.
  !   - Although the approximation has been built around the line element, it is a good estimation if used directly to
  !     quadrilaterals and triangles, taking into acount which curve have to be used for each characteristic path (see previous
  !     figure).

  ! ====================
  ! STANDARD INTEGRATION
  ! ====================

  ! -----
  ! ln(r)
  ! -----

  ! alpha0(x)=c0+c1*x+c2*x**2, x=[log10(relative error)]
  ! lnr_alpha0(curve,coefficient)
  real(kind=real64), parameter :: lnr_alpha0(3,0:2)=reshape((/ 2.511371433253553d-1, 1.588915580092415d-1, 1.303338290205561d-1,&
                                                              -8.318131267946077d-2,-9.469173675969373d-2,-9.041069352366727d-2,&
                                                              -1.973643360180158d-3,-2.500979779162619d-3,-2.302123487272636d-3/),&
                                                            (/3,3/))
  ! alpha1(x)=c0+c1*x+c2*x**2, x=[log10(relative error)]
  ! lnr_alpha1(curve,coefficient)
  real(kind=real64), parameter :: lnr_alpha1(3,0:2)=reshape((/-8.301769338383449d-1,-8.575272011921337d-1,-5.785210418300009d-1,&
                                                              -7.458331284021328d-2,-6.866562987114574d-2,-3.107089548772586d-2,&
                                                              -1.992535165956408d-3,-1.735469217262563d-3,-4.204818399993008d-4/),&
                                                            (/3,3/))
  ! beta1(x)=c0+c1*x+c2*x**2, x=[log10(relative error)]
  ! lnr_beta1(curve,coefficient)
  real(kind=real64), parameter :: lnr_beta1(3,0:2)=reshape((/-1.114783338380396d-1,-3.918234478283774d-1,-1.933844506192423d-1,&
                                                             -6.169581618278328d-2,-6.478096580939043d-2,-3.303326826061181d-2,&
                                                             -1.873698813959166d-3,-1.657459305149584d-3,-6.497990057394490d-4/),&
                                                           (/3,3/))

  ! ---------------
  ! 1/r^n, n=1,..,7
  ! ---------------

  ! alpha0(x,y)=c0+c1*x+c2*y+c3*x*y+c4*x**2+c5*y**2, x=[log10(relative error)], y=[integrand (1,2,...,7)]
  ! d1rn_alpha0(curve,coefficient)
  real(kind=real64), parameter :: d1rn_alpha0(3,0:5)=reshape((/ 2.443943976922257d-1, 1.855869272054077d-1, 1.964994120809251d-1,&
                                                               -7.845125775746743d-2,-8.441784317192081d-2,-7.533682025056084d-2,&
                                                                3.895326800087784d-2, 3.992519918023225d-2, 6.300967117113785d-2,&
                                                                8.783124822346978d-4, 1.003034542932314d-3, 1.453521059340246d-3,&
                                                               -1.719172638314362d-3,-2.012460287158036d-3,-1.653933461489553d-3,&
                                                               -1.575912106080523d-3,-1.596684136739970d-3,-2.613368441217468d-3/),&
                                                             (/3,6/))
  ! alpha1(x,y)=c0+c1*x+c2*y+c3*x*y+c4*x**2+c5*y**2, x=[log10(relative error)], y=[integrand (1,2,...,7)]
  ! d1rn_alpha1(curve,coefficient)
  real(kind=real64), parameter :: d1rn_alpha1(3,0:5)=reshape((/-8.421216556191472d-1,-8.012899104725332d-1,-5.891778953073477d-1,&
                                                               -8.778262996280104d-2,-7.060561516342915d-2,-4.515354064860566d-2,&
                                                                1.679496012553309d-2, 5.713299243166406d-3, 5.728946548799128d-3,&
                                                                6.522695414940627d-4,-1.952251227204736d-4, 1.106131193536092d-6,&
                                                               -2.649807213627071d-3,-1.984455335909311d-3,-1.179100966502172d-3,&
                                                               -8.472009602597869d-4,-6.680823196695776d-4,-5.355782068751523d-4/),&
                                                             (/3,6/))
  ! beta1(x,y)=c0+c1*x+c2*y+c3*x*y+c4*x**2+c5*y**2, x=[log10(relative error)], y=[integrand (1,2,...,7)]
  ! d1rn_beta1(curve,coefficient)
  real(kind=real64), parameter :: d1rn_beta1(3,0:5)=reshape((/-2.693032219216048d-2,-2.353515879579821d-1,-9.828561879516765d-2,&
                                                              -6.138311391493901d-2,-5.141856801187955d-2,-3.074365737570111d-2,&
                                                               2.273921394673689d-2, 1.693100431819354d-2, 1.306905121742232d-2,&
                                                               9.789745450155518d-4,-1.491381226735667d-5, 6.065437465834727d-5,&
                                                              -2.030377872751224d-3,-1.330369292483563d-3,-7.853234287066991d-4,&
                                                              -9.059289455187508d-4,-1.264704521827383d-3,-1.025889829482388d-3/),&
                                                            (/3,6/))

  ! ==================
  ! TELLES INTEGRATION
  ! ==================

  ! -----
  ! ln(r)
  ! -----

  ! alpha0(x)=c0+c1*x+c2*x**2, x=[log10(relative error)]
  ! lnr_alpha0(curve,coefficient)
  real(kind=real64), parameter :: telles_lnr_alpha0(3,0:2)=reshape((/ 2.5834148732821033d-01, 5.1283007739426911d-01, 4.1061632629507161d-01,&
                                                                     -8.8865167808209675d-02,-1.1479138136155583d-02,-3.0669024691563598d-02,&
                                                                     -2.3640194116454241d-03, 1.7167160355773211d-03, 5.1920314153891956d-04/),&
                                                                   (/3,3/))
  ! alpha1(x)=c0+c1*x+c2*x**2, x=[log10(relative error)]
  ! lnr_alpha1(curve,coefficient)
  real(kind=real64), parameter :: telles_lnr_alpha1(3,0:2)=reshape((/-9.1721308433764848d-01,-1.0632414684959128d+00,-1.1740312565911180d+00,&
                                                                     -6.2837979237008845d-02,-1.1633428696469847d-01,-1.3946434974037192d-01,&
                                                                     -1.1580362636310136d-03,-4.4173578665990879d-03,-5.4264705520737193d-03/),&
                                                                   (/3,3/))
  ! beta1(x)=c0+c1*x+c2*x**2, x=[log10(relative error)]
  ! lnr_beta1(curve,coefficient)
  real(kind=real64), parameter :: telles_lnr_beta1(3,0:2)=reshape((/-4.9431525769444484d-01,-9.1742055835524161d-01,-1.1668939380105636d+00,&
                                                                    -6.5628053690049357d-02,-1.4448614058834844d-01,-1.7686541251294283d-01,&
                                                                    -1.2683068366630545d-03,-5.5849846784841173d-03,-6.9562836980345301d-03/),&
                                                                  (/3,3/))

  ! ---------------
  ! 1/r^n, n=1,..,7
  ! ---------------

  ! alpha0(x,y)=c0+c1*x+c2*y+c3*x*y+c4*x**2+c5*y**2, x=[log10(relative error)], y=[integrand (1,2,...,7)]
  ! d1rn_alpha0(curve,coefficient)
  real(kind=real64), parameter :: telles_d1rn_alpha0(3,0:5)=reshape((/ 2.9379042486857926d-01, 1.1155167367380846d-01, 9.5587534756324991d-02,&
                                                                      -7.6518438657183094d-02,-9.0043220465500617d-02,-8.5841370957250140d-02,&
                                                                       3.0011332045426935d-02, 3.9605615471752335d-02, 5.9724590614136340d-02,&
                                                                       5.8900790667561719d-04, 1.0669688965738889d-03, 1.4101322853742803d-03,&
                                                                      -1.7642732320416377d-03,-2.2428131063082419d-03,-2.0424217927592597d-03,&
                                                                      -1.3634192833061319d-03,-1.6030309471606837d-03,-2.4264102214015559d-03/),&
                                                                    (/3,6/))
  ! alpha1(x,y)=c0+c1*x+c2*y+c3*x*y+c4*x**2+c5*y**2, x=[log10(relative error)], y=[integrand (1,2,...,7)]
  ! d1rn_alpha1(curve,coefficient)
  real(kind=real64), parameter :: telles_d1rn_alpha1(3,0:5)=reshape((/-9.8339088502232042d-01,-6.5856321275500784d-01,-2.4634605754204883d-01,&
                                                                      -8.6761356304300027d-02,-4.7354561095619759d-02, 1.0051411705899942d-02,&
                                                                       2.6028760963988956d-02, 1.9423558756522371d-03,-1.1994978255365170d-02,&
                                                                       7.7677507397256337d-04,-2.9087437598165618d-04,-1.0022065157916291d-03,&
                                                                      -2.3619290417702245d-03,-1.1859755938049318d-03, 7.9639034155470554d-04,&
                                                                      -1.2891431252005461d-03,-4.5988734937654318d-04,-1.8103634920482468d-04/),&
                                                                    (/3,6/))
  ! beta1(x,y)=c0+c1*x+c2*y+c3*x*y+c4*x**2+c5*y**2, x=[log10(relative error)], y=[integrand (1,2,...,7)]
  ! d1rn_beta1(curve,coefficient)
  real(kind=real64), parameter :: telles_d1rn_beta1(3,0:5)=reshape((/-4.5453529921363467d-01,-2.6610069546799336d-01, 7.5605400079843565d-02,&
                                                                     -7.3941529478921839d-02,-3.6048387668689889d-02, 2.4021121564322980d-02,&
                                                                      3.7672240427426076d-02, 1.3288129149913631d-02,-9.3693044556669332d-03,&
                                                                      9.8276802658567980d-04,-1.8099231920766130d-04,-1.3532246088221670d-03,&
                                                                     -1.8897770518175436d-03,-7.6015457149682889d-04, 1.3564674989692248d-03,&
                                                                     -1.8079562352345093d-03,-1.2247286683263207d-03,-7.4467072775070487d-04/),&
                                                                   (/3,6/))

  ! -------------------------------------
  ! PARAMETERS FOR A GIVEN RELATIVE ERROR
  ! -------------------------------------

  !! Quasi-singular parameters for the estimation of the number of Gausss-Legendre points for a given relative error.
  type fbem_qs_parameters
    ! Relative error
    real(kind=real64) :: re
    ! Standard integration
    real(kind=real64) :: dmin(3,0:7)
    real(kind=real64) :: dmax(3,0:7)
    real(kind=real64) :: alpha0(3,0:7)
    real(kind=real64) :: alpha1(3,0:7)
    real(kind=real64) :: beta1(3,0:7)
    ! Telles integration
    real(kind=real64) :: telles_dmin(3,0:7)
    real(kind=real64) :: telles_dmax(3,0:7)
    real(kind=real64) :: telles_alpha0(3,0:7)
    real(kind=real64) :: telles_alpha1(3,0:7)
    real(kind=real64) :: telles_beta1(3,0:7)
  end type fbem_qs_parameters

contains

  ! ----- !
  ! SETUP !
  ! ------!

  !! Calculate quasi-singular integration parameters
  subroutine fbem_qs_calculate_parameters(relative_error,p)
    implicit none
    real(kind=real64)        :: relative_error
    type(fbem_qs_parameters) :: p
    real(kind=real64)        :: log10re
    integer                  :: c, i
    real(kind=real64)        :: c1d(0:2), c2d(0:5)
    ! Adjust the relative error within valid ranges
    log10re=log10(abs(relative_error))
    if (log10re.gt.-3.d0) log10re=-3.d0
    if (log10re.lt.-15.d0) log10re=-15.d0
    p%re=10.d0**log10re
    !
    ! Standard integration parameters
    !
    ! Integrand ln(r)
    do c=1,3
      ! Calculate alpha0, alpha1 and beta1
      c1d=lnr_alpha0(c,:)
      p%alpha0(c,0)=f_poly_1d_quadratic(c1d,log10re)
      c1d=lnr_alpha1(c,:)
      p%alpha1(c,0)=f_poly_1d_quadratic(c1d,log10re)
      c1d=lnr_beta1(c,:)
      p%beta1(c,0)=f_poly_1d_quadratic(c1d,log10re)
      ! Calculate dmin for N=30, dmax for N=2 => x=log10(d), y=log10(N), x=(y-alpha0)/(alpha1-y*beta1)
      p%dmin(c,0)=10.d0**((log10(30.)-p%alpha0(c,0))/(p%alpha1(c,0)-log10(30.)*p%beta1(c,0)))
      p%dmax(c,0)=10.d0**((log10( 2.)-p%alpha0(c,0))/(p%alpha1(c,0)-log10( 2.)*p%beta1(c,0)))
    end do
    ! Integrand 1/r^n, n=1,..,7
    do i=1,7
      do c=1,3
        ! Calculate alpha0, alpha1 and beta1
        c2d=d1rn_alpha0(c,:)
        p%alpha0(c,i)=f_poly_2d_quadratic(c2d,log10re,dble(i))
        c2d=d1rn_alpha1(c,:)
        p%alpha1(c,i)=f_poly_2d_quadratic(c2d,log10re,dble(i))
        c2d=d1rn_beta1(c,:)
        p%beta1(c,i)=f_poly_2d_quadratic(c2d,log10re,dble(i))
        ! Calculate dmin for N=30, dmax for N=2 => x=log10(d), y=log10(N), x=(y-alpha0)/(alpha1-y*beta1)
        p%dmin(c,i)=10.d0**((log10(30.)-p%alpha0(c,i))/(p%alpha1(c,i)-log10(30.)*p%beta1(c,i)))
        p%dmax(c,i)=10.d0**((log10( 2.)-p%alpha0(c,i))/(p%alpha1(c,i)-log10( 2.)*p%beta1(c,i)))
      end do
    end do
    !
    ! Telles integration parameters
    !
    ! Integrand ln(r)
    do c=1,3
      ! Calculate alpha0, alpha1 and beta1
      c1d=telles_lnr_alpha0(c,:)
      p%telles_alpha0(c,0)=f_poly_1d_quadratic(c1d,log10re)
      c1d=telles_lnr_alpha1(c,:)
      p%telles_alpha1(c,0)=f_poly_1d_quadratic(c1d,log10re)
      c1d=telles_lnr_beta1(c,:)
      p%telles_beta1(c,0)=f_poly_1d_quadratic(c1d,log10re)
      ! Calculate dmin for N=30, dmax for N=2 => x=log10(d), y=log10(N), x=(y-alpha0)/(alpha1-y*beta1)
      p%telles_dmin(c,0)=10.d0**((log10(30.)-p%telles_alpha0(c,0))/(p%telles_alpha1(c,0)-log10(30.)*p%telles_beta1(c,0)))
      p%telles_dmax(c,0)=10.d0**((log10( 2.)-p%telles_alpha0(c,0))/(p%telles_alpha1(c,0)-log10( 2.)*p%telles_beta1(c,0)))
      ! Correct ranges
      if (p%telles_dmin(c,0).lt.1.d-6) p%telles_dmin(c,0)=1.d-6
      if (p%telles_dmin(c,0).gt.p%telles_dmax(c,0)) p%telles_dmin(c,0)=1.d-6
    end do
    ! Integrand 1/r^n, n=1,..,7
    do i=1,7
      do c=1,3
        ! Calculate alpha0, alpha1 and beta1
        c2d=telles_d1rn_alpha0(c,:)
        p%telles_alpha0(c,i)=f_poly_2d_quadratic(c2d,log10re,dble(i))
        c2d=telles_d1rn_alpha1(c,:)
        p%telles_alpha1(c,i)=f_poly_2d_quadratic(c2d,log10re,dble(i))
        c2d=telles_d1rn_beta1(c,:)
        p%telles_beta1(c,i)=f_poly_2d_quadratic(c2d,log10re,dble(i))
        ! Calculate dmin for N=30, dmax for N=2 => x=log10(d), y=log10(N), x=(y-alpha0)/(alpha1-y*beta1)
        p%telles_dmin(c,i)=10.d0**((log10(30.)-p%telles_alpha0(c,i))/(p%telles_alpha1(c,i)-log10(30.)*p%telles_beta1(c,i)))
        p%telles_dmax(c,i)=10.d0**((log10( 2.)-p%telles_alpha0(c,i))/(p%telles_alpha1(c,i)-log10( 2.)*p%telles_beta1(c,i)))
        ! Correct ranges
        if (p%telles_dmin(c,i).lt.1.d-6) p%telles_dmin(c,i)=1.d-6
        if (p%telles_dmin(c,i).gt.p%telles_dmax(c,i)) p%telles_dmin(c,i)=1.d-6
      end do
    end do
  end subroutine fbem_qs_calculate_parameters

  !! Function that calculates an estimation of the number of Gauss-Legendre points needed to integrate a quasi-singular function
  !! in a standard way.
  function fbem_qs_n_estimation_standard(rn,etype,f,p,d,barxi,r,q)
    implicit none
    ! I/O
    integer                     :: fbem_qs_n_estimation_standard   !! If 0, d is out of range. If 2<=N<=30 a correct estimation is given.
    integer                     :: rn                              !! Spatial dimension
    integer                     :: etype                           !! Element type: line2, line3, quad4, quad8, quad9
    integer                     :: f                               !! 0: ln(r); 1-7: 1/r^f
    type(fbem_qs_parameters)    :: p                               !! Quasi-singular integration parameters
    real(kind=real64)           :: d                               !! Dimensionless distance d=r/L
    real(kind=real64)           :: barxi(fbem_n_dimension(etype)) !! Local coordinate of the nearest point of the element to the collocation point. -1<=barxi<=1 for lineX and quadX, and 0<=barxi<=1 for triX
    real(kind=real64), optional :: r(rn)                           !! Dimensionless distance vector: r=(x-x_i)/|x-x_i|
    real(kind=real64), optional :: q(rn)                           !! Orientation vector: the average unit tangent vector for a lineX, and the average unit normal vector for triX and quadX.
    ! Local
    real(kind=real64)           :: gamma
    real(kind=real64)           :: N1, N2, N3
    real(kind=real64)           :: log10d
    !
    ! d<=2 (a precise estimation of d is required)
    !
    if (d.le.2.d0) then
      ! log10(d)
      log10d=log10(d)
      ! Interpolation between both curves
      select case (fbem_n_dimension(etype))
        ! LINE ELEMENT
        case (1)
          ! if -1<barxi<1 (ZONE INSIDE), interpolate between N1 and N2
          if (abs(barxi(1)).lt.1.d0) then
            ! N1(d)
            if (d.le.p%dmin(1,f)) then
              fbem_qs_n_estimation_standard=0
              return
            else
              if (d.ge.p%dmax(1,f)) then
                N1=2.d0
              else
                N1=10.d0**((p%alpha0(1,f)+p%alpha1(1,f)*log10d)/(1.d0+p%beta1(1,f)*log10d))
              end if
            end if
            ! N2(d)
            if (d.le.p%dmin(2,f)) then
              fbem_qs_n_estimation_standard=0
              return
            else
              if (d.ge.p%dmax(2,f)) then
                N2=2.d0
              else
                N2=10.d0**((p%alpha0(2,f)+p%alpha1(2,f)*log10d)/(1.d0+p%beta1(2,f)*log10d))
              end if
            end if
            ! Interpolate between N1 and N2
            fbem_qs_n_estimation_standard=ceiling(barxi(1)**2*(N2-N1)+N1)
          ! if barxi==1 or barxi==-1 (ZONE OUTSIDE), interpolate between N2 and N3
          else
            if (present(r).and.present(q)) then
              ! Relative orientation (gamma in [0,1])
              gamma=abs(dot_product(r,q))
              if (gamma.gt.1.d0) then
                call fbem_error_message(output_unit,0,__FILE__,__LINE__,'r and q must be unit vectors')
              end if
              ! N2(d)
              if (d.le.p%dmin(2,f)) then
                fbem_qs_n_estimation_standard=0
                return
              else
                if (d.ge.p%dmax(2,f)) then
                  N2=2.d0
                else
                  N2=10.d0**((p%alpha0(2,f)+p%alpha1(2,f)*log10d)/(1.d0+p%beta1(2,f)*log10d))
                end if
              end if
              ! N3(d)
              if (d.le.p%dmin(3,f)) then
                fbem_qs_n_estimation_standard=0
                return
              else
                if (d.ge.p%dmax(3,f)) then
                  N3=2.d0
                else
                  N3=10.d0**((p%alpha0(3,f)+p%alpha1(3,f)*log10d)/(1.d0+p%beta1(3,f)*log10d))
                end if
              end if
              ! Interpolate between N2 and N3
              fbem_qs_n_estimation_standard=ceiling((1.d0-gamma)*N2+gamma*N3)
            else
              ! N2(d)
              if (d.le.p%dmin(2,f)) then
                fbem_qs_n_estimation_standard=0
                return
              else
                if (d.ge.p%dmax(2,f)) then
                  N2=2.d0
                else
                  N2=10.d0**((p%alpha0(2,f)+p%alpha1(2,f)*log10d)/(1.d0+p%beta1(2,f)*log10d))
                end if
              end if
              fbem_qs_n_estimation_standard=ceiling(N2)
            endif
          end if
        ! TRIANGULAR OR QUADRILATERAL ELEMENTS
        case (2)
          ! N1(d)
          if (d.le.p%dmin(1,f)) then
            fbem_qs_n_estimation_standard=0
            return
          else
            if (d.ge.p%dmax(1,f)) then
              N1=2.d0
            else
              N1=10.d0**((p%alpha0(1,f)+p%alpha1(1,f)*log10d)/(1.d0+p%beta1(1,f)*log10d))
            end if
          end if
          ! N2(d)
          if (d.le.p%dmin(2,f)) then
            fbem_qs_n_estimation_standard=0
            return
          else
            if (d.ge.p%dmax(2,f)) then
              N2=2.d0
            else
              N2=10.d0**((p%alpha0(2,f)+p%alpha1(2,f)*log10d)/(1.d0+p%beta1(2,f)*log10d))
            end if
          end if
          ! Interpolate between N1 and N2 by using the barxi coordinate
          select case (fbem_n_edges(etype))
            ! TRIANGULAR
            case (3)
              fbem_qs_n_estimation_standard=ceiling(4.0d0*(barxi(2)*(barxi(2)+barxi(1)-1.0d0)+(barxi(1)-1.0d0)*barxi(1))*(N2-N1)+N2)
            ! QUADRILATERAL
            case (4)
              fbem_qs_n_estimation_standard=ceiling((barxi(1)**2)*(barxi(2)**2)*(N2-N1)+N1)
            case default
              call fbem_error_message(output_unit,0,__FILE__,__LINE__,'incorrect element type')
          end select
        case default
          call fbem_error_message(output_unit,0,__FILE__,__LINE__,'incorrect element type')
      end select
    !
    ! d>2 (a rough estimation is enough)
    !
    else
      ! N3(d), but using nint() instead of ceiling(), this produces approximately an average between N2(d) and N3(d).
      if (d.ge.p%dmax(3,f)) then
        fbem_qs_n_estimation_standard=2
      else
        log10d=log10(d)
        fbem_qs_n_estimation_standard=nint(10.d0**((p%alpha0(3,f)+p%alpha1(3,f)*log10d)/(1.d0+p%beta1(3,f)*log10d)))
      end if
    end if
    if (fbem_qs_n_estimation_standard.gt.30) then
      fbem_qs_n_estimation_standard=0
    end if
  end function fbem_qs_n_estimation_standard

  !! Function that calculates an estimation of the number of Gauss-Legendre points needed to integrate a quasi-singular function
  !! using the Telles transformation.
  function fbem_qs_n_estimation_telles(rn,etype,f,p,d,barxi,r,q)
    implicit none
    ! I/O
    integer                     :: fbem_qs_n_estimation_telles     !! If 0, d is out of range. If 2<=N<=30 a correct estimation is given.
    integer                     :: rn                              !! Spatial dimension
    integer                     :: etype                           !! Element type: line2, line3, quad4, quad8, quad9
    integer                     :: f                               !! 0: ln(r); 1-7: 1/r^f
    type(fbem_qs_parameters)    :: p                               !! Quasi-singular integration parameters
    real(kind=real64)           :: d                               !! Dimensionless distance d=r/L
    real(kind=real64)           :: barxi(fbem_n_dimension(etype)) !! Local coordinate of the nearest point of the element to the collocation point. -1<=barxi<=1 for lineX and quadX, and 0<=barxi<=1 for triX
    real(kind=real64), optional :: r(rn)                           !! Dimensionless distance vector: r=(x-x_i)/|x-x_i|
    real(kind=real64), optional :: q(rn)                           !! Orientation vector: the average unit tangent vector for a lineX, and the average unit normal vector for triX and quadX.
    ! Local
    real(kind=real64)           :: gamma
    real(kind=real64)           :: N1, N2, N3
    real(kind=real64)           :: log10d
    !
    ! d<=2 (a precise estimation of d is required)
    !
    if (d.le.2.d0) then
      ! log10(d)
      log10d=log10(d)
      ! Interpolation between both curves
      select case (fbem_n_dimension(etype))
        ! LINE ELEMENT
        case (1)
          ! if -1<barxi<1 (ZONE INSIDE), interpolate between N1 and N2
          if (abs(barxi(1)).lt.1.d0) then
            ! N1(d)
            if (d.le.p%telles_dmin(1,f)) then
              fbem_qs_n_estimation_telles=0
              return
            else
              if (d.ge.p%telles_dmax(1,f)) then
                N1=2.d0
              else
                N1=10.d0**((p%telles_alpha0(1,f)+p%telles_alpha1(1,f)*log10d)/(1.d0+p%telles_beta1(1,f)*log10d))
              end if
            end if
            ! N2(d)
            if (d.le.p%telles_dmin(2,f)) then
              fbem_qs_n_estimation_telles=0
              return
            else
              if (d.ge.p%telles_dmax(2,f)) then
                N2=2.d0
              else
                N2=10.d0**((p%telles_alpha0(2,f)+p%telles_alpha1(2,f)*log10d)/(1.d0+p%telles_beta1(2,f)*log10d))
              end if
            end if
            ! Interpolate between N1 and N2
            fbem_qs_n_estimation_telles=ceiling(barxi(1)**2*(N2-N1)+N1)
          ! if barxi==1 or barxi==-1 (ZONE OUTSIDE), interpolate between N2 and N3
          else
            if (present(r).and.present(q)) then
              ! Relative orientation (gamma in [0,1])
              gamma=abs(dot_product(r,q))
              if (gamma.gt.1.d0) then
                call fbem_error_message(output_unit,0,__FILE__,__LINE__,'r and q must be unit vectors')
              end if
              ! N2(d)
              if (d.le.p%telles_dmin(2,f)) then
                fbem_qs_n_estimation_telles=0
                return
              else
                if (d.ge.p%telles_dmax(2,f)) then
                  N2=2.d0
                else
                  N2=10.d0**((p%telles_alpha0(2,f)+p%telles_alpha1(2,f)*log10d)/(1.d0+p%telles_beta1(2,f)*log10d))
                end if
              end if
              ! N3(d)
              if (d.le.p%telles_dmin(3,f)) then
                fbem_qs_n_estimation_telles=0
                return
              else
                if (d.ge.p%telles_dmax(3,f)) then
                  N3=2.d0
                else
                  N3=10.d0**((p%telles_alpha0(3,f)+p%telles_alpha1(3,f)*log10d)/(1.d0+p%telles_beta1(3,f)*log10d))
                end if
              end if
              ! Interpolate between N2 and N3
              fbem_qs_n_estimation_telles=ceiling((1.d0-gamma)*N2+gamma*N3)
            else
              ! N2(d)
              if (d.le.p%telles_dmin(2,f)) then
                fbem_qs_n_estimation_telles=0
                return
              else
                if (d.ge.p%telles_dmax(2,f)) then
                  N2=2.d0
                else
                  N2=10.d0**((p%telles_alpha0(2,f)+p%telles_alpha1(2,f)*log10d)/(1.d0+p%telles_beta1(2,f)*log10d))
                end if
              end if
              fbem_qs_n_estimation_telles=ceiling(N2)
            end if
          end if
        ! TRIANGULAR OR QUADRILATERAL ELEMENTS
        case (2)
          ! N1(d)
          if (d.le.p%telles_dmin(1,f)) then
            fbem_qs_n_estimation_telles=0
            return
          else
            if (d.ge.p%telles_dmax(1,f)) then
              N1=2.d0
            else
              N1=10.d0**((p%telles_alpha0(1,f)+p%telles_alpha1(1,f)*log10d)/(1.d0+p%telles_beta1(1,f)*log10d))
            end if
          end if
          ! N2(d)
          if (d.le.p%telles_dmin(2,f)) then
            fbem_qs_n_estimation_telles=0
            return
          else
            if (d.ge.p%telles_dmax(2,f)) then
              N2=2.d0
            else
              N2=10.d0**((p%telles_alpha0(2,f)+p%telles_alpha1(2,f)*log10d)/(1.d0+p%telles_beta1(2,f)*log10d))
            end if
          end if
          ! Interpolate between N1 and N2 by using the barxi coordinate
          select case (fbem_n_edges(etype))
            ! TRIANGULAR
            case (3)
              fbem_qs_n_estimation_telles=ceiling(4.0d0*(barxi(2)*(barxi(2)+barxi(1)-1.0d0)+(barxi(1)-1.0d0)*barxi(1))*(N2-N1)+N2)
            ! QUADRILATERAL
            case (4)
              fbem_qs_n_estimation_telles=ceiling((barxi(1)**2)*(barxi(2)**2)*(N2-N1)+N1)
            case default
              call fbem_error_message(output_unit,0,__FILE__,__LINE__,'incorrect element type')
          end select
        case default
          call fbem_error_message(output_unit,0,__FILE__,__LINE__,'incorrect element type')
      end select
    !
    ! d>2 (a rough estimation is enough)
    ! Telles transformation does not improve the integrand, but for completeness an estimation is given.
    !
    else
      ! N3(d), but using nint() instead of ceiling(), this produces approximately an average between N2(d) and N3(d).
      if (d.ge.p%dmax(3,f)) then
        fbem_qs_n_estimation_telles=2
      else
        log10d=log10(d)
        fbem_qs_n_estimation_telles=nint(10.d0**((p%alpha0(3,f)+p%alpha1(3,f)*log10d)/(1.d0+p%beta1(3,f)*log10d)))
      end if
    end if
    if (fbem_qs_n_estimation_telles.gt.30) then
      fbem_qs_n_estimation_telles=0
    end if
  end function fbem_qs_n_estimation_telles

  ! Number of Gauss-Legendre points needed to integrate phi*jacobian integrand
  ! --------------------------------------------------------------------------

  !! Calculate number of Gauss points needed for shape function integration for 1D elements in 2D space
  subroutine fbem_qs_phijac_ngp_1d(rn,type_g,type_f,x_nodes,error_height,ngp,achieve)
    implicit none
    ! I/O
    integer           :: rn                               !! Dimensional space
    integer           :: type_g                           !! Geometrical interpolation
    integer           :: type_f                           !! Functional interpolation
    real(kind=real64) :: x_nodes(rn,fbem_n_nodes(type_g)) !! Position of element nodes
    real(kind=real64) :: error_height                     !! Error height [1.0E-12,1.0E-01]. If the error height is <=0, 1.0E-6 is used.
    integer           :: ngp                              !! Number of Gauss-Legendre points needed to integrate <tt>phi*jacobian</tt> integrand for 1D elements.
    logical           :: achieve                          !! True if the number of Gauss points achieves the error height, false otherwise.
    ! Local
    real(kind=real64) :: local_error_height
    real(kind=real64) :: integral(fbem_n_nodes(type_f))
    real(kind=real64) :: integral_old(fbem_n_nodes(type_f))
    real(kind=real64) :: phi_f(fbem_n_nodes(type_f))
    real(kind=real64) :: dphidxi_g(fbem_n_nodes(type_g))
    real(kind=real64) :: aux(10)
    real(kind=real64) :: xi, t(rn), jacobian
    logical           :: check_continue
    integer           :: j, k, l, nnodes_g, nnodes_f
    ! Save the number of nodes
    nnodes_g=fbem_n_nodes(type_g)
    nnodes_f=fbem_n_nodes(type_f)
    ! If error_height is zero or less, default error_height of 1.0E6 is used, if less than 1.0E-12, 1.0E-12 i used
    if (error_height.le.0.0d0) then
      local_error_height=1.0d-6
    else
      if (error_height.lt.1.0d-12) then
        local_error_height=1.0d-12
      else
        local_error_height=error_height
      end if
    end if
    ! Initialize loop that increments number of Gauss points until error_height is reached
    check_continue=.true.
    achieve=.true.
    integral=0.d0
    j=1
    do while (check_continue.eqv.(.true.))
      ! Save old value of integrals and initialize for a new calculation
      integral_old=integral
      integral=0.0d0
      ! Calculate integrals
      do k=1,gl11_n(j)
        xi=gl11_xi(k,j)
        ! Get phi and dphidxi
#       define etype type_g
#       define delta 0.0d0
#       define dphidxi dphidxi_g
#       include <dphidxi_1d.rc>
#       undef etype
#       undef delta
#       undef dphidxi
        ! Calculate tangent vector
        t=0.d0
        do l=1,nnodes_g
          t=t+dphidxi_g(l)*x_nodes(:,l)
        end do
        ! Jacobian as the module of the tangent vector
        jacobian=sqrt(dot_product(t,t))
        ! Get phi_f
#       define etype type_f
#       define phi phi_f
#       define delta 0.0d0
#       include <phi_1d.rc>
#       undef etype
#       undef phi
#       undef delta
        ! Calculate integrals for each phi_f (note: 1 is added to phi_f in order to ensure that the integral is >0)
        integral=integral+(phi_f+1.0d0)*jacobian*gl11_w(k,j)
      end do
      ! Check if all integrals have error less than error_height, except for j=1 (first integral evaluation)
      if (j.eq.1) then
        j=j+1
      else
        check_continue=.false.
        do l=1,nnodes_f
          ! If relative error > error height, is necessary to increase the number of Gauss points.
          if (dabs((integral_old(l)-integral(l))/integral(l)).gt.local_error_height) then
            check_continue=.true.
            j=j+1
            exit
          end if
        end do
        if (j.eq.(gl11_nr+1)) then
          check_continue=.false.
          achieve=.false.
        end if
      end if
    end do
    ! Return the number of Gauss points that integrate all integrals with error_height
    ngp=j-1
  end subroutine fbem_qs_phijac_ngp_1d

  !! Calculate number of Gauss-Legendre points needed for shape function integration for 2D elements in 3D space
  subroutine fbem_qs_phijac_ngp_2d(rn,type_g,type_f,x_nodes,error_height,ngp,achieve)
    implicit none
    ! I/O
    integer           :: rn                               !! Dimensional space
    integer           :: type_g                           !! Geometrical interpolation
    integer           :: type_f                           !! Functional interpolation
    real(kind=real64) :: x_nodes(rn,fbem_n_nodes(type_g)) !! Position of element nodes
    real(kind=real64) :: error_height                     !! Error height [1.0E-12,1.0E-01]. If the error height is <=0, 1.0E-6 is used.
    integer           :: ngp                              !! Number of Gauss-Legendre points needed to integrate <tt>phi*jacobian</tt> integrand for 2D elements (for each direction).
    logical           :: achieve                          !! True if the number of Gauss points achieves the error height, false otherwise.
    ! Local
    real(kind=real64) :: local_error_height
    real(kind=real64) :: integral(fbem_n_nodes(type_f))
    real(kind=real64) :: integral_old(fbem_n_nodes(type_f))
    real(kind=real64) :: phi_f(fbem_n_nodes(type_f))
    real(kind=real64) :: dphidxi1(fbem_n_nodes(type_g))
    real(kind=real64) :: dphidxi2(fbem_n_nodes(type_g))
    real(kind=real64) :: aux(10),xi(2)
    real(kind=real64) :: T1(rn),T2(rn), N(rn), jacobian, jw
    logical           :: check_continue
    integer           :: j, k1, k2, l, nnodes_g, nnodes_f
    ! Number of nodes of the element
    nnodes_g=fbem_n_nodes(type_g)
    nnodes_f=fbem_n_nodes(type_f)
    ! If error_height is zero or less, default error_height of 1.0E6 is used, if less than 1.0E-12, 1.0E-12 i used
    if (error_height.le.0.0d0) then
      local_error_height=1.0d-6
    else
      if (error_height.lt.1.0d-12) then
        local_error_height=1.0d-12
      else
        local_error_height=error_height
      end if
    end if
    ! Initialize loop that increments number of Gauss points until error_height is reached
    check_continue=.true.
    achieve=.true.
    integral=0.d0
    j=1
    do while (check_continue.eqv.(.true.))
      ! Save old value of integrals and initialize for a new calculation
      integral_old=integral
      integral=0.0d0
      ! Calculate integrals depending if quadrilateral or triangular
      select case (type_g)
        ! Triangular elements
        case (fbem_tri3,fbem_tri6)
          do k1=1,gl01_n(j)
            do k2=1,gj01_n(j)
              xi(1)=(1.0d0-gj01_xi(k2,j))*gl01_xi(k1,j)
              xi(2)=gj01_xi(k2,j)
              ! Get phi, dphidxi1 and dphidxi2
#             define etype type_g
#             define delta 0.0d0
#             include <dphidxik_2d.rc>
#             undef etype
#             undef delta
              ! Calculate tangents vectors
              T1=0.0d0
              T2=0.0d0
              do l=1,nnodes_g
                  T1=T1+dphidxi1(l)*x_nodes(:,l)
                  T2=T2+dphidxi2(l)*x_nodes(:,l)
              end do
              ! Jacobian as the module of the normal vector
              select case (rn)
                case (2)
                  ! Jacobian
                  jacobian=T1(1)*T2(2)-T1(2)*T2(1)
                case (3)
                  ! Normal vector
                  N(1)=T1(2)*T2(3)-T1(3)*T2(2)
                  N(2)=T1(3)*T2(1)-T1(1)*T2(3)
                  N(3)=T1(1)*T2(2)-T1(2)*T2(1)
                  ! Jacobian
                  jacobian=sqrt(dot_product(N,N))
              end select
              ! Jacobian * weights
              jw=jacobian*gl01_w(k1,j)*gj01_w(k2,j)
              ! Get phi_f
#             define etype type_f
#             define phi phi_f
#             define delta 0.0d0
#             include <phi_2d.rc>
#             undef etype
#             undef phi
#             undef delta
              ! Calculate integrals for each phi_f (note: 1 is added to fphi in order to ensure that the integral is >0)
              integral=integral+(phi_f+1.0d0)*jw
            end do
          end do
        ! Quadrilateral elements
        case (fbem_quad4,fbem_quad8,fbem_quad9)
          do k1=1,gl11_n(j)
            xi(1)=gl11_xi(k1,j)
            do k2=1,gl11_n(j)
              xi(2)=gl11_xi(k2,j)
              ! Get phi, dphidxi1 and dphidxi2
#             define etype type_g
#             define delta 0.0d0
#             include <dphidxik_2d.rc>
#             undef etype
#             undef delta
              ! Calculate tangents vectors
              T1=0.0d0
              T2=0.0d0
              do l=1,nnodes_g
                  T1=T1+dphidxi1(l)*x_nodes(:,l)
                  T2=T2+dphidxi2(l)*x_nodes(:,l)
              end do
              ! Jacobian as the module of the normal vector
              select case (rn)
                case (2)
                  ! Jacobian
                  jacobian=T1(1)*T2(2)-T1(2)*T2(1)
                case (3)
                  ! Normal vector
                  N(1)=T1(2)*T2(3)-T1(3)*T2(2)
                  N(2)=T1(3)*T2(1)-T1(1)*T2(3)
                  N(3)=T1(1)*T2(2)-T1(2)*T2(1)
                  ! Jacobian
                  jacobian=sqrt(dot_product(N,N))
              end select
              ! Jacobian * weights
              jw=jacobian*gl11_w(k1,j)*gl11_w(k2,j)
              ! Get phi_f
#             define etype type_f
#             define phi phi_f
#             define delta 0.0d0
#             include <phi_2d.rc>
#             undef etype
#             undef phi
#             undef delta
              ! Calculate integrals for each phi_f (note: 1 is added to fphi in order to ensure that the integral is >0)
              integral=integral+(phi_f+1.0d0)*jw
            end do
          end do
        case default
          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'element type undefined')
      end select
      ! Check if all integrals have error less than error_height, except for j=1 (first integral evaluation)
      if (j.eq.1) then
        j=j+1
      else
        ! By default, there's no next iteration
        check_continue=.false.
        ! Loop through errors for each shape function * jacobian integral
        do l=1,nnodes_f
          ! If relative error > error height, is necessary to increase the number of Gauss points.
          if (dabs((integral_old(l)-integral(l))/integral(l)).gt.local_error_height) then
            check_continue=.true.
            j=j+1
            exit
          end if
        end do
        if (j.eq.(gl11_nr+1)) then
          check_continue=.false.
          achieve=.false.
        end if
      end if
    end do
    ! Return the number of Gauss points that integrate all integrals with error_height
    ngp=j-1
  end subroutine fbem_qs_phijac_ngp_2d


  subroutine fbem_qs_test_standard_line

    implicit none

    real(kind=real128) :: xnodes(2,2), x_i(2), w, phi(2), x(2), t(2), rv(2), r, jg
    real(kind=real128) :: xi
    real(kind=real128) :: eps
    real(kind=real128) :: d, d_min, d_max, delta_d
    integer            :: e_min, delta_e, e_max
    integer            :: i, j, k, e, npoints
    real(kind=real128) :: ref_integral, integral, pdensity
    integer            :: alpha
    integer            :: nd
    integer            :: idistance, ini_nd, ini_nd_tmp
    integer            :: ngp_min, ngp_max, delta_ngp, ngp
    integer, allocatable :: function_n(:)
    real(kind=real128), allocatable :: function_d(:)
    character(len=32)  :: fmtstr

    fmtstr='(3i11,e25.16,i11)'

    ! Reference element (line2)
    xnodes(1,1)=0.0d0
    xnodes(1,2)=0.0d0
    xnodes(2,1)=1.0d0
    xnodes(2,2)=0.0d0

    ! Distance discretization
    d_min=1.d-6
    d_max=1.d6
    delta_d=0.001d0
    nd=nint((log10(d_max)-log10(d_min))/delta_d)+1

    ! Error height discretization
    e_min=-15
    e_max=-3
    delta_e=1

    ! Number of Gauss Points used
    ngp_min=2
    ngp_max=30
    delta_ngp=1
    allocate (function_n(ngp_min:ngp_max))
    allocate (function_d(ngp_min:ngp_max))

    ! File configuration
    open(unit=60,file='n_d_standard_line.dat',recl=1024)

    ! Header
    write(60,*) '# alpha, log10(eps), ch. path, d, N(d)'
    write(60,*) '# alpha: integrand 1/r^alpha'
    write(60,*) '# ch. path (characteristic path) for line element at 0<=x<=1, y=0:'
    write(60,*) '#   1 x=0.5 , y=d (out-of-plane center)'
    write(60,*) '#   2 x=1   , y=d (out-of-plane vertex)'
    write(60,*) '#   3 x=1+d , y=0 (in-plane vertex)'

    ! Loop through alpha
    do alpha=1,9
      write(*,*) 'integrand type = ', alpha

      ! Loop through error heights
      do e=e_min,e_max,delta_e
        write(*,*) ' log10(relative_error) = ', e

        ! ------------------------------------------
        ! Collocation point in center (out-of-plane)
        ! ------------------------------------------

        write(*,*) '  Collocation point in center (out-of-plane)'

        ! Initialize
        function_n=0
        function_d=0.d0
        ini_nd_tmp=nd

        ! Loop though necessary number of gauss points
        do ngp=ngp_min,ngp_max,delta_ngp

          ini_nd=ini_nd_tmp
          ! Loop through the distance d
          do idistance=ini_nd,1,-1
            ! Generate the collocation point
            d=10.0d0**(dble(idistance-1)*delta_d)*d_min
            x_i(1)=0.5d0
            x_i(2)=d

!            if ((ngp.eq.ngp_min).and.(e.eq.e_min)) then
!              ref_integral=0.
!              do k=1,gl11_n(32)
!                xi=gl11_xi(k,32)
!                w=gl11_w(k,32)
!                phi(1)=0.5*(1.-xi)
!                phi(2)=0.5*(1.+xi)
!                x=phi(1)*xnodes(1,:)+phi(2)*xnodes(2,:)
!                rv=x-x_i
!                r=sqrt(dot_product(rv,rv))
!                t=-0.5*xnodes(1,:)+0.5*xnodes(2,:)
!                jg=sqrt(dot_product(t,t))
!                if (alpha.eq.0) then
!                  ref_integral=ref_integral+log(r)*jg*w
!                else
!                  ref_integral=ref_integral+(1./r**alpha)*jg*w
!                end if
!              end do
!              write(21,*) alpha, d, ref_integral
!            end if

            !
            ! Calculate the reference numerical solution with ngp+2
            !
            ref_integral=0.
            do k=1,gl11_n(ngp+2)
              xi=gl11_xi(k,ngp+2)
              w=gl11_w(k,ngp+2)
              phi(1)=0.5*(1.-xi)
              phi(2)=0.5*(1.+xi)
              x=phi(1)*xnodes(1,:)+phi(2)*xnodes(2,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              t=-0.5*xnodes(1,:)+0.5*xnodes(2,:)
              jg=sqrt(dot_product(t,t))
              if (alpha.eq.0) then
                ref_integral=ref_integral+log(r)*jg*w
              else
                ref_integral=ref_integral+(1./r**alpha)*jg*w
              end if
            end do
            !
            ! Calculate the approximate integral solution with ngp
            !
            integral=0.0d0
            do k=1,gl11_n(ngp)
              xi=gl11_xi(k,ngp)
              w=gl11_w(k,ngp)
              phi(1)=0.5*(1.-xi)
              phi(2)=0.5*(1.+xi)
              x=phi(1)*xnodes(1,:)+phi(2)*xnodes(2,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              t=-0.5*xnodes(1,:)+0.5*xnodes(2,:)
              jg=sqrt(dot_product(t,t))
              if (alpha.eq.0) then
                integral=integral+log(r)*jg*w
              else
                integral=integral+(1./r**alpha)*jg*w
              end if
            end do

            !
            ! Stop criterion
            !
            eps=abs((integral-ref_integral)/ref_integral)
            if (eps.gt.10.0d0**e) then
              function_d(ngp)=d
              function_n(ngp)=ngp
              ini_nd_tmp=idistance
              exit
            end if

          end do ! Loop through the distance d

        end do ! Loop though necessary number of gauss points

        ! Write
        do ngp=ngp_max,ngp_min,-delta_ngp
          if (function_n(ngp).ne.0) then
            write(60,fmtstr) alpha, e, 1, function_d(ngp), function_n(ngp)
          end if
        end do
        write(60,*)
        write(60,*)

        ! ------------------------------------------
        ! Collocation point in vertex (out-of-plane)
        ! ------------------------------------------

        write(*,*) '  Collocation point in vertex y'

        ! Initialize
        function_n=0
        function_d=0.d0
        ini_nd_tmp=nd

        ! Loop though necessary number of gauss points
        do ngp=ngp_min,ngp_max,delta_ngp

          ini_nd=ini_nd_tmp
          ! Loop through the distance d
          do idistance=ini_nd,1,-1
            ! Generate the collocation point
            d=10.0d0**(dble(idistance-1)*delta_d)*d_min
            x_i(1)=0.0d0
            x_i(2)=d

            !
            ! Calculate the reference numerical solution with ngp+2
            !
            ref_integral=0.0d0
            do k=1,gl11_n(ngp+2)
              xi=gl11_xi(k,ngp+2)
              w=gl11_w(k,ngp+2)
              phi(1)=0.5*(1.-xi)
              phi(2)=0.5*(1.+xi)
              x=phi(1)*xnodes(1,:)+phi(2)*xnodes(2,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              t=-0.5*xnodes(1,:)+0.5*xnodes(2,:)
              jg=sqrt(dot_product(t,t))
              if (alpha.eq.0) then
                ref_integral=ref_integral+log(r)*jg*w
              else
                ref_integral=ref_integral+(1./r**alpha)*jg*w
              end if
            end do
            !
            ! Calculate the approximate integral solution with ngp
            !
            integral=0.0d0
            do k=1,gl11_n(ngp)
              xi=gl11_xi(k,ngp)
              w=gl11_w(k,ngp)
              phi(1)=0.5*(1.-xi)
              phi(2)=0.5*(1.+xi)
              x=phi(1)*xnodes(1,:)+phi(2)*xnodes(2,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              t=-0.5*xnodes(1,:)+0.5*xnodes(2,:)
              jg=sqrt(dot_product(t,t))
              if (alpha.eq.0) then
                integral=integral+log(r)*jg*w
              else
                integral=integral+(1./r**alpha)*jg*w
              end if
            end do

            !
            ! Stop criterion
            !
            eps=abs((integral-ref_integral)/ref_integral)
            if (eps.gt.10.0d0**e) then
              function_d(ngp)=d
              function_n(ngp)=ngp
              ini_nd_tmp=idistance
              exit
            end if

          end do ! Loop through the distance d

        end do ! Loop though necessary number of gauss points

        ! Write
        do ngp=ngp_max,ngp_min,-delta_ngp
          if (function_n(ngp).ne.0) then
            write(60,fmtstr) alpha, e, 2, function_d(ngp), function_n(ngp)
          end if
        end do
        write(60,*)
        write(60,*)

        ! --------------------------------------
        ! Collocation point in vertex (in-plane)
        ! --------------------------------------

        write(*,*) '  Collocation point in vertex (in-plane)'

        ! Initialize
        function_n=0
        function_d=0.d0
        ini_nd_tmp=nd

        ! Loop though necessary number of gauss points
        do ngp=ngp_min,ngp_max,delta_ngp

          ini_nd=ini_nd_tmp
          ! Loop through the distance d
          do idistance=ini_nd,1,-1
            ! Generate the collocation point
            d=10.0d0**(dble(idistance-1)*delta_d)*d_min
            x_i(1)=-d
            x_i(2)=0.0d0

            !
            ! Calculate the reference numerical solution with ngp+2
            !
            ref_integral=0.0d0
            do k=1,gl11_n(ngp+2)
              xi=gl11_xi(k,ngp+2)
              w=gl11_w(k,ngp+2)
              phi(1)=0.5*(1.-xi)
              phi(2)=0.5*(1.+xi)
              x=phi(1)*xnodes(1,:)+phi(2)*xnodes(2,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              t=-0.5*xnodes(1,:)+0.5*xnodes(2,:)
              jg=sqrt(dot_product(t,t))
              if (alpha.eq.0) then
                ref_integral=ref_integral+log(r)*jg*w
              else
                ref_integral=ref_integral+(1./r**alpha)*jg*w
              end if
            end do
            !
            ! Calculate the approximate integral solution with ngp
            !
            integral=0.0d0
            do k=1,gl11_n(ngp)
              xi=gl11_xi(k,ngp)
              w=gl11_w(k,ngp)
              phi(1)=0.5*(1.-xi)
              phi(2)=0.5*(1.+xi)
              x=phi(1)*xnodes(1,:)+phi(2)*xnodes(2,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              t=-0.5*xnodes(1,:)+0.5*xnodes(2,:)
              jg=sqrt(dot_product(t,t))
              if (alpha.eq.0) then
                integral=integral+log(r)*jg*w
              else
                integral=integral+(1./r**alpha)*jg*w
              end if
            end do

            !
            ! Stop criterion
            !
            eps=abs((integral-ref_integral)/ref_integral)
            if (eps.gt.10.0d0**e) then
              function_d(ngp)=d
              function_n(ngp)=ngp
              ini_nd_tmp=idistance
              exit
            end if

          end do ! Loop through the distance d

        end do ! Loop though necessary number of gauss points

        ! Write
        do ngp=ngp_max,ngp_min,-delta_ngp
          if (function_n(ngp).ne.0) then
            write(60,fmtstr) alpha, e, 3, function_d(ngp), function_n(ngp)
          end if
        end do
        write(60,*)
        write(60,*)

      end do ! Loop through error heights

    end do ! Loop through alpha

    close(60)

  end subroutine fbem_qs_test_standard_line

  subroutine fbem_qs_test_disttrans_line

    implicit none

    real(kind=real128) :: xnodes(2,2), x_i(2), w, phi(2), x(2), t(2), rv(2), r, jg
    real(kind=real128) :: xi
    real(kind=real128) :: eps
    real(kind=real128) :: d, d_min, d_max, delta_d
    real(kind=real128) :: a, c, eta, etap, eta1, eta2, dxideta, detadetap
    integer            :: e_min, delta_e, e_max
    integer            :: i, j, k, e, npoints
    real(kind=real128) :: ref_integral, integral, pdensity
    integer            :: alpha
    integer            :: nd
    integer            :: idistance, ini_nd, ini_nd_tmp
    integer            :: ngp_min, ngp_max, delta_ngp, ngp
    integer, allocatable :: function_n(:)
    real(kind=real128), allocatable :: function_d(:)
    character(len=32)  :: fmtstr

    fmtstr='(3i11,e25.16,i11)'

    ! Reference element (line2)
    xnodes(1,1)=0.0d0
    xnodes(1,2)=0.0d0
    xnodes(2,1)=1.0d0
    xnodes(2,2)=0.0d0

    ! Distance discretization
    d_min=1.d-6
    d_max=1.d6
    delta_d=0.001d0
    nd=nint((log10(d_max)-log10(d_min))/delta_d)+1

    ! Error height discretization
    e_min=-15
    e_max=-3
    delta_e=1

    ! Number of Gauss Points used
    ngp_min=2
    ngp_max=30
    delta_ngp=1
    allocate (function_n(ngp_min:ngp_max))
    allocate (function_d(ngp_min:ngp_max))

    ! File configuration
    open(unit=60,file='n_d_disttrans_line.dat',recl=1024)

    ! Header
    write(60,*) '# alpha, log10(eps), ch. path, d, N(d)'
    write(60,*) '# alpha: integrand 1/r^alpha'
    write(60,*) '# ch. path (characteristic path) for line element at 0<=x<=1, y=0:'
    write(60,*) '#   1 x=0.5 , y=d (out-of-plane center)'
    write(60,*) '#   2 x=1   , y=d (out-of-plane vertex)'
    write(60,*) '#   3 x=1+d , y=0 (in-plane vertex)'

    ! Loop through alpha
    do alpha=1,9
      write(*,*) 'integrand type = ', alpha

      ! Loop through error heights
      do e=e_min,e_max,delta_e
        write(*,*) ' log10(relative_error) = ', e

        ! ------------------------------------------
        ! Collocation point in center (out-of-plane)
        ! ------------------------------------------

        write(*,*) '  Collocation point in center (out-of-plane)'

        ! Initialize
        function_n=0
        function_d=0.d0
        ini_nd_tmp=nd

        ! Loop though necessary number of gauss points
        do ngp=ngp_min,ngp_max,delta_ngp

          ini_nd=ini_nd_tmp
          ! Loop through the distance d
          do idistance=ini_nd,1,-1
            ! Generate the collocation point
            d=10.0d0**(dble(idistance-1)*delta_d)*d_min
            x_i(1)=0.5d0
            x_i(2)=d

            c=0.d0
            a=(1.d0+1.d-3)*d/0.5d0

            ! 1/r
            !eta1=log(sqrt(a**2+(1.d0+c)**2)-(1.d0+c))
            !eta2=log(sqrt(a**2+(1.d0-c)**2)+(1.d0-c))

            ! 1/r^2
            !eta1=atan((-1.d0-c)/a)/a
            !eta2=atan(( 1.d0-c)/a)/a

            ! 1/r^3 (new)
            eta1=(-1.d0-c)/sqrt(a**2+(-1.d0-c)**2)/a**2
            eta2=( 1.d0-c)/sqrt(a**2+( 1.d0-c)**2)/a**2


            !
            ! Calculate the reference numerical solution with ngp+2
            !
            ref_integral=0.
            do k=1,gl11_n(ngp+2)
              etap=gl11_xi(k,ngp+2)
              w=gl11_w(k,ngp+2)
              eta=0.5d0*(1.d0-etap)*eta1+0.5d0*(1.d0+etap)*eta2
              detadetap=0.5d0*(eta2-eta1)

              ! 1/r
              !xi=0.5d0*(exp(eta)-a**2*exp(-eta))+c
              !dxideta=sqrt(a**2+(xi-c)**2)

              ! 1/r^2
              !xi=a*tan(a*eta)+c
              !dxideta=a**2+(xi-c)**2

              ! 1/r^3
              if ((1.d0-(a**2*eta)**2).gt.0) then
                xi=c+a**3*eta/sqrt(1.d0-(a**2*eta)**2)
                dxideta=sqrt(a**2+(xi-c)**2)**3
              else
                stop 'Houston: we have a problem'
              end if


              phi(1)=0.5*(1.-xi)
              phi(2)=0.5*(1.+xi)
              x=phi(1)*xnodes(1,:)+phi(2)*xnodes(2,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              t=-0.5*xnodes(1,:)+0.5*xnodes(2,:)
              jg=sqrt(dot_product(t,t))
              if (alpha.eq.0) then
                ref_integral=ref_integral+log(r)*jg*dxideta*detadetap*w
              else
                ref_integral=ref_integral+(dxideta/r)/r**(alpha-1)*jg*detadetap*w
              end if
            end do
            !
            ! Calculate the approximate integral solution with ngp
            !
            integral=0.
            do k=1,gl11_n(ngp)
              etap=gl11_xi(k,ngp)
              w=gl11_w(k,ngp)
              eta=0.5d0*(1.d0-etap)*eta1+0.5d0*(1.d0+etap)*eta2
              detadetap=0.5d0*(eta2-eta1)

              ! 1/r
              !xi=0.5d0*(exp(eta)-a**2*exp(-eta))+c
              !dxideta=sqrt(a**2+(xi-c)**2)

              ! 1/r^2
              !xi=a*tan(a*eta)+c
              !dxideta=a**2+(xi-c)**2

              ! 1/r^3
              if ((1.d0-(a**2*eta)**2).gt.0) then
                xi=c+a**3*eta/sqrt(1.d0-(a**2*eta)**2)
                dxideta=sqrt(a**2+(xi-c)**2)**3
              else
                stop 'Houston: we have a problem'
              end if

              phi(1)=0.5*(1.-xi)
              phi(2)=0.5*(1.+xi)
              x=phi(1)*xnodes(1,:)+phi(2)*xnodes(2,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              t=-0.5*xnodes(1,:)+0.5*xnodes(2,:)
              jg=sqrt(dot_product(t,t))
              if (alpha.eq.0) then
                integral=integral+log(r)*jg*dxideta*detadetap*w
              else
                integral=integral+(dxideta/r)/r**(alpha-1)*jg*detadetap*w
              end if
            end do

            !
            ! Stop criterion
            !
            eps=abs((integral-ref_integral)/ref_integral)
            if (eps.gt.10.0d0**e) then
              function_d(ngp)=d
              function_n(ngp)=ngp
              ini_nd_tmp=idistance
              exit
            end if

          end do ! Loop through the distance d

        end do ! Loop though necessary number of gauss points

        ! Write
        do ngp=ngp_max,ngp_min,-delta_ngp
          if (function_n(ngp).ne.0) then
            write(60,fmtstr) alpha, e, 1, function_d(ngp), function_n(ngp)
          end if
        end do
        write(60,*)
        write(60,*)

        ! ------------------------------------------
        ! Collocation point in vertex (out-of-plane)
        ! ------------------------------------------

        write(*,*) '  Collocation point in vertex (out-of-plane)'

        ! Initialize
        function_n=0
        function_d=0.d0
        ini_nd_tmp=nd

        ! Loop though necessary number of gauss points
        do ngp=ngp_min,ngp_max,delta_ngp

          ini_nd=ini_nd_tmp
          ! Loop through the distance d
          do idistance=ini_nd,1,-1
            ! Generate the collocation point
            d=10.0d0**(dble(idistance-1)*delta_d)*d_min
            x_i(1)=1.d0
            x_i(2)=d

            c=1.d0
            a=(1.d0+1.d-3)*d/0.5d0

            ! 1/r
            !eta1=log(sqrt(a**2+(1.d0+c)**2)-(1.d0+c))
            !eta2=log(sqrt(a**2+(1.d0-c)**2)+(1.d0-c))

            ! 1/r^2
            !eta1=atan((-1.d0-c)/a)/a
            !eta2=atan(( 1.d0-c)/a)/a

            ! 1/r^3 (new)
            eta1=(-1.d0-c)/sqrt(a**2+(-1.d0-c)**2)/a**2
            eta2=( 1.d0-c)/sqrt(a**2+( 1.d0-c)**2)/a**2


            !
            ! Calculate the reference numerical solution with ngp+2
            !
            ref_integral=0.
            do k=1,gl11_n(32)
              etap=gl11_xi(k,32)
              w=gl11_w(k,32)
              eta=0.5d0*(1.d0-etap)*eta1+0.5d0*(1.d0+etap)*eta2
              detadetap=0.5d0*(eta2-eta1)

              ! 1/r
              !xi=0.5d0*(exp(eta)-a**2*exp(-eta))+c
              !dxideta=sqrt(a**2+(xi-c)**2)

              ! 1/r^2
              !xi=a*tan(a*eta)+c
              !dxideta=a**2+(xi-c)**2

              ! 1/r^3
              if ((1.d0-(a**2*eta)**2).gt.0) then
                xi=c+a**3*eta/sqrt(1.d0-(a**2*eta)**2)
                dxideta=sqrt(a**2+(xi-c)**2)**3
              else
                stop 'Houston: we have a problem'
              end if

              phi(1)=0.5*(1.-xi)
              phi(2)=0.5*(1.+xi)
              x=phi(1)*xnodes(1,:)+phi(2)*xnodes(2,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              t=-0.5*xnodes(1,:)+0.5*xnodes(2,:)
              jg=sqrt(dot_product(t,t))
              if (alpha.eq.0) then
                ref_integral=ref_integral+log(r)*jg*dxideta*detadetap*w
              else
                ref_integral=ref_integral+(dxideta/r)/r**(alpha-1)*jg*detadetap*w
              end if
            end do
            !
            ! Calculate the approximate integral solution with ngp
            !
            integral=0.
            do k=1,gl11_n(ngp)
              etap=gl11_xi(k,ngp)
              w=gl11_w(k,ngp)
              eta=0.5d0*(1.d0-etap)*eta1+0.5d0*(1.d0+etap)*eta2
              detadetap=0.5d0*(eta2-eta1)

              ! 1/r
              !xi=0.5d0*(exp(eta)-a**2*exp(-eta))+c
              !dxideta=sqrt(a**2+(xi-c)**2)

              ! 1/r^2
              !xi=a*tan(a*eta)+c
              !dxideta=a**2+(xi-c)**2

              ! 1/r^3
              if ((1.d0-(a**2*eta)**2).gt.0) then
                xi=c+a**3*eta/sqrt(1.d0-(a**2*eta)**2)
                dxideta=sqrt(a**2+(xi-c)**2)**3
              else
                stop 'Houston: we have a problem'
              end if

              phi(1)=0.5*(1.-xi)
              phi(2)=0.5*(1.+xi)
              x=phi(1)*xnodes(1,:)+phi(2)*xnodes(2,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              t=-0.5*xnodes(1,:)+0.5*xnodes(2,:)
              jg=sqrt(dot_product(t,t))
              if (alpha.eq.0) then
                integral=integral+log(r)*jg*dxideta*detadetap*w
              else
                integral=integral+(dxideta/r)/r**(alpha-1)*jg*detadetap*w
              end if
            end do

            !
            ! Stop criterion
            !
            eps=abs((integral-ref_integral)/ref_integral)
            if (eps.gt.10.0d0**e) then
              function_d(ngp)=d
              function_n(ngp)=ngp
              ini_nd_tmp=idistance
              exit
            end if

          end do ! Loop through the distance d

        end do ! Loop though necessary number of gauss points

        ! Write
        do ngp=ngp_max,ngp_min,-delta_ngp
          if (function_n(ngp).ne.0) then
            write(60,fmtstr) alpha, e, 2, function_d(ngp), function_n(ngp)
          end if
        end do
        write(60,*)
        write(60,*)

        ! ------------------------------------------
        ! Collocation point in vertex (in-plane)
        ! ------------------------------------------

        write(*,*) '  Collocation point in vertex (in-plane)'

        ! Initialize
        function_n=0
        function_d=0.d0
        ini_nd_tmp=nd

        ! Loop though necessary number of gauss points
        do ngp=ngp_min,ngp_max,delta_ngp

          ini_nd=ini_nd_tmp
          ! Loop through the distance d
          do idistance=ini_nd,1,-1
            ! Generate the collocation point
            d=10.0d0**(dble(idistance-1)*delta_d)*d_min
            x_i(1)=1.d0+d
            x_i(2)=0.d0

            c=1.d0
            a=(1.d0+1.d-3)*d/0.5d0

            ! 1/r
            !eta1=log(sqrt(a**2+(1.d0+c)**2)-(1.d0+c))
            !eta2=log(sqrt(a**2+(1.d0-c)**2)+(1.d0-c))

            ! 1/r^2
            !eta1=atan((-1.d0-c)/a)/a
            !eta2=atan(( 1.d0-c)/a)/a

            ! 1/r^3 (new)
            eta1=(-1.d0-c)/sqrt(a**2+(-1.d0-c)**2)/a**2
            eta2=( 1.d0-c)/sqrt(a**2+( 1.d0-c)**2)/a**2


            !
            ! Calculate the reference numerical solution with ngp+2
            !
            ref_integral=0.
            do k=1,gl11_n(ngp+2)
              etap=gl11_xi(k,ngp+2)
              w=gl11_w(k,ngp+2)
              eta=0.5d0*(1.d0-etap)*eta1+0.5d0*(1.d0+etap)*eta2
              detadetap=0.5d0*(eta2-eta1)

              ! 1/r
              !xi=0.5d0*(exp(eta)-a**2*exp(-eta))+c
              !dxideta=sqrt(a**2+(xi-c)**2)

              ! 1/r^2
              !xi=a*tan(a*eta)+c
              !dxideta=a**2+(xi-c)**2

              ! 1/r^3
              if ((1.d0-(a**2*eta)**2).gt.0) then
                xi=c+a**3*eta/sqrt(1.d0-(a**2*eta)**2)
                dxideta=sqrt(a**2+(xi-c)**2)**3
              else
                stop 'Houston: we have a problem'
              end if

              phi(1)=0.5*(1.-xi)
              phi(2)=0.5*(1.+xi)
              x=phi(1)*xnodes(1,:)+phi(2)*xnodes(2,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              t=-0.5*xnodes(1,:)+0.5*xnodes(2,:)
              jg=sqrt(dot_product(t,t))
              if (alpha.eq.0) then
                ref_integral=ref_integral+log(r)*jg*dxideta*detadetap*w
              else
                ref_integral=ref_integral+(dxideta/r)/r**(alpha-1)*jg*detadetap*w
              end if
            end do
            !
            ! Calculate the approximate integral solution with ngp
            !
            integral=0.
            do k=1,gl11_n(ngp)
              etap=gl11_xi(k,ngp)
              w=gl11_w(k,ngp)
              eta=0.5d0*(1.d0-etap)*eta1+0.5d0*(1.d0+etap)*eta2
              detadetap=0.5d0*(eta2-eta1)

              ! 1/r
              !xi=0.5d0*(exp(eta)-a**2*exp(-eta))+c
              !dxideta=sqrt(a**2+(xi-c)**2)

              ! 1/r^2
              !xi=a*tan(a*eta)+c
              !dxideta=a**2+(xi-c)**2

              ! 1/r^3
              if ((1.d0-(a**2*eta)**2).gt.0) then
                xi=c+a**3*eta/sqrt(1.d0-(a**2*eta)**2)
                dxideta=sqrt(a**2+(xi-c)**2)**3
              else
                stop 'Houston: we have a problem'
              end if

              phi(1)=0.5*(1.-xi)
              phi(2)=0.5*(1.+xi)
              x=phi(1)*xnodes(1,:)+phi(2)*xnodes(2,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              t=-0.5*xnodes(1,:)+0.5*xnodes(2,:)
              jg=sqrt(dot_product(t,t))
              if (alpha.eq.0) then
                integral=integral+log(r)*jg*dxideta*detadetap*w
              else
                integral=integral+(dxideta/r)/r**(alpha-1)*jg*detadetap*w
              end if
            end do

            !
            ! Stop criterion
            !
            eps=abs((integral-ref_integral)/ref_integral)
            if (eps.gt.10.0d0**e) then
              function_d(ngp)=d
              function_n(ngp)=ngp
              ini_nd_tmp=idistance
              exit
            end if

          end do ! Loop through the distance d

        end do ! Loop though necessary number of gauss points

        ! Write
        do ngp=ngp_max,ngp_min,-delta_ngp
          if (function_n(ngp).ne.0) then
            write(60,fmtstr) alpha, e, 3, function_d(ngp), function_n(ngp)
          end if
        end do
        write(60,*)
        write(60,*)

      end do ! Loop through error heights

    end do ! Loop through alpha

    close(60)

  end subroutine fbem_qs_test_disttrans_line

  subroutine fbem_qs_test_standard_quad

    implicit none

    real(kind=real128) :: xnodes(4,3), x_i(3), w, x(3), t1(3), t2(3), n(3), rv(3), r, jg
    real(kind=real128) :: phi(4), dphidxi(4,2)
    real(kind=real128) :: xi(2)
    real(kind=real128) :: eps
    real(kind=real128) :: d, d_min, d_max, delta_d
    integer            :: e_min, delta_e, e_max
    integer            :: i, j, k1, k2, e, npoints
    real(kind=real128) :: ref_integral, integral, pdensity
    integer            :: alpha
    integer            :: nd
    integer            :: idistance, ini_nd, ini_nd_tmp
    integer            :: ngp_min, ngp_max, delta_ngp, ngp
    integer, allocatable :: function_n(:)
    real(kind=real128), allocatable :: function_d(:)
    character(len=32)  :: fmtstr

    fmtstr='(3i11,e25.16,i11)'

    ! Reference element (quad4)
    xnodes=0.
    xnodes(1,1)=0.
    xnodes(1,2)=0.
    xnodes(2,1)=1.
    xnodes(2,2)=0.
    xnodes(3,1)=1.
    xnodes(3,2)=1.
    xnodes(4,1)=0.
    xnodes(4,2)=1.

    ! Distance discretization
    d_min=1.d-6
    d_max=1.d6
    delta_d=0.001d0
    nd=nint((log10(d_max)-log10(d_min))/delta_d)+1

    ! Error height discretization
    e_min=-15
    e_max=-3
    delta_e=1

    ! Number of Gauss Points used
    ngp_min=2
    ngp_max=30
    delta_ngp=1
    allocate (function_n(ngp_min:ngp_max))
    allocate (function_d(ngp_min:ngp_max))

    ! File configuration
    open(unit=60,file='n_d_standard_quad.dat',recl=1024)

    ! Header
    write(60,*) '# alpha, log10(eps), ch. path, d, N(d)'
    write(60,*) '# alpha: integrand 1/r^alpha'
    write(60,*) '# ch. path (characteristic path) for quad element at 0<=x<=1, 0<=y<=1, z=0:'
    write(60,*) '#   1 x=0.5 , y=0.5 , z=d'
    write(60,*) '#   2 x=0.5 , y=1.0 , z=d'
    write(60,*) '#   3 x=1   , y=1   , z=d'
    write(60,*) '#   4 x=0.5 , y=1+d , z=0'
    write(60,*) '#   5 x=1   , y=1+d , z=0'

    ! Loop through alpha
    do alpha=1,9
      write(*,*) 'integrand type = ', alpha

      ! Loop through error heights
      do e=e_min,e_max,delta_e
        write(*,*) ' log10(relative_error) = ', e

        ! ------------------------------------------
        ! Collocation point in center (out-of-plane)
        ! ------------------------------------------

        write(*,*) '  Collocation point in center (out-of-plane)'

        ! Initialize
        function_n=0
        function_d=0.d0
        ini_nd_tmp=nd

        ! Loop though necessary number of gauss points
        do ngp=ngp_min,ngp_max,delta_ngp

          ini_nd=ini_nd_tmp
          ! Loop through the distance d
          do idistance=ini_nd,1,-1
            ! Generate the collocation point
            d=10.0d0**(dble(idistance-1)*delta_d)*d_min
            x_i(1)=0.5
            x_i(2)=0.5
            x_i(3)=d

            !
            ! Calculate the reference numerical solution with ngp+2
            !
            ref_integral=0.
            do k1=1,gl11_n(ngp+2)
            do k2=1,gl11_n(ngp+2)
              xi(1)=gl11_xi(k1,ngp+2)
              xi(2)=gl11_xi(k2,ngp+2)
              w=gl11_w(k1,ngp+2)*gl11_w(k2,ngp+2)
              phi(1)=0.25*(1.-xi(1))*(1.-xi(2))
              phi(2)=0.25*(1.+xi(1))*(1.-xi(2))
              phi(3)=0.25*(1.+xi(1))*(1.+xi(2))
              phi(4)=0.25*(1.-xi(1))*(1.+xi(2))
              dphidxi(1,1)=-0.25*(1.-xi(2))
              dphidxi(2,1)= 0.25*(1.-xi(2))
              dphidxi(3,1)= 0.25*(1.+xi(2))
              dphidxi(4,1)=-0.25*(1.+xi(2))
              dphidxi(1,2)=-0.25*(1.-xi(1))
              dphidxi(2,2)=-0.25*(1.+xi(1))
              dphidxi(3,2)= 0.25*(1.+xi(1))
              dphidxi(4,2)= 0.25*(1.-xi(1))
              x=  phi(1)*xnodes(1,:)
              x=x+phi(2)*xnodes(2,:)
              x=x+phi(3)*xnodes(3,:)
              x=x+phi(4)*xnodes(4,:)
              t1=   dphidxi(1,1)*xnodes(1,:)
              t1=t1+dphidxi(2,1)*xnodes(2,:)
              t1=t1+dphidxi(3,1)*xnodes(3,:)
              t1=t1+dphidxi(4,1)*xnodes(4,:)
              t2=   dphidxi(1,2)*xnodes(1,:)
              t2=t2+dphidxi(2,2)*xnodes(2,:)
              t2=t2+dphidxi(3,2)*xnodes(3,:)
              t2=t2+dphidxi(4,2)*xnodes(4,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              n=fbem_cross_product(t1,t2)
              jg=sqrt(dot_product(n,n))
              if (alpha.eq.0) then
                ref_integral=ref_integral+log(r)*jg*w
              else
                ref_integral=ref_integral+(1./r**alpha)*jg*w
              end if
            end do
            end do
            !
            ! Calculate the approximate integral solution with ngp
            !
            integral=0.
            do k1=1,gl11_n(ngp)
            do k2=1,gl11_n(ngp)
              xi(1)=gl11_xi(k1,ngp)
              xi(2)=gl11_xi(k2,ngp)
              w=gl11_w(k1,ngp)*gl11_w(k2,ngp)
              phi(1)=0.25*(1.-xi(1))*(1.-xi(2))
              phi(2)=0.25*(1.+xi(1))*(1.-xi(2))
              phi(3)=0.25*(1.+xi(1))*(1.+xi(2))
              phi(4)=0.25*(1.-xi(1))*(1.+xi(2))
              dphidxi(1,1)=-0.25*(1.-xi(2))
              dphidxi(2,1)= 0.25*(1.-xi(2))
              dphidxi(3,1)= 0.25*(1.+xi(2))
              dphidxi(4,1)=-0.25*(1.+xi(2))
              dphidxi(1,2)=-0.25*(1.-xi(1))
              dphidxi(2,2)=-0.25*(1.+xi(1))
              dphidxi(3,2)= 0.25*(1.+xi(1))
              dphidxi(4,2)= 0.25*(1.-xi(1))
              x=  phi(1)*xnodes(1,:)
              x=x+phi(2)*xnodes(2,:)
              x=x+phi(3)*xnodes(3,:)
              x=x+phi(4)*xnodes(4,:)
              t1=   dphidxi(1,1)*xnodes(1,:)
              t1=t1+dphidxi(2,1)*xnodes(2,:)
              t1=t1+dphidxi(3,1)*xnodes(3,:)
              t1=t1+dphidxi(4,1)*xnodes(4,:)
              t2=   dphidxi(1,2)*xnodes(1,:)
              t2=t2+dphidxi(2,2)*xnodes(2,:)
              t2=t2+dphidxi(3,2)*xnodes(3,:)
              t2=t2+dphidxi(4,2)*xnodes(4,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              n=fbem_cross_product(t1,t2)
              jg=sqrt(dot_product(n,n))
              if (alpha.eq.0) then
                integral=integral+log(r)*jg*w
              else
                integral=integral+(1./r**alpha)*jg*w
              end if
            end do
            end do

            !
            ! Stop criterion
            !
            eps=abs((integral-ref_integral)/ref_integral)
            if (eps.gt.10.0d0**e) then
              function_d(ngp)=d
              function_n(ngp)=ngp
              ini_nd_tmp=idistance
              exit
            end if

          end do ! Loop through the distance d

        end do ! Loop though necessary number of gauss points

        ! Write
        do ngp=ngp_max,ngp_min,-delta_ngp
          if (function_n(ngp).ne.0) then
            write(60,fmtstr) alpha, e, 1, function_d(ngp), function_n(ngp)
          end if
        end do
        write(60,*)
        write(60,*)

        ! -----------------------------------------------
        ! Collocation point in edge center (out-of-plane)
        ! -----------------------------------------------

        write(*,*) '  Collocation point in edge center (out-of-plane)'

        ! Initialize
        function_n=0
        function_d=0.d0
        ini_nd_tmp=nd

        ! Loop though necessary number of gauss points
        do ngp=ngp_min,ngp_max,delta_ngp

          ini_nd=ini_nd_tmp
          ! Loop through the distance d
          do idistance=ini_nd,1,-1
            ! Generate the collocation point
            d=10.0d0**(dble(idistance-1)*delta_d)*d_min
            x_i(1)=0.5
            x_i(2)=1.
            x_i(3)=d

            !
            ! Calculate the reference numerical solution with ngp+2
            !
            ref_integral=0.
            do k1=1,gl11_n(ngp+2)
            do k2=1,gl11_n(ngp+2)
              xi(1)=gl11_xi(k1,ngp+2)
              xi(2)=gl11_xi(k2,ngp+2)
              w=gl11_w(k1,ngp+2)*gl11_w(k2,ngp+2)
              phi(1)=0.25*(1.-xi(1))*(1.-xi(2))
              phi(2)=0.25*(1.+xi(1))*(1.-xi(2))
              phi(3)=0.25*(1.+xi(1))*(1.+xi(2))
              phi(4)=0.25*(1.-xi(1))*(1.+xi(2))
              dphidxi(1,1)=-0.25*(1.-xi(2))
              dphidxi(2,1)= 0.25*(1.-xi(2))
              dphidxi(3,1)= 0.25*(1.+xi(2))
              dphidxi(4,1)=-0.25*(1.+xi(2))
              dphidxi(1,2)=-0.25*(1.-xi(1))
              dphidxi(2,2)=-0.25*(1.+xi(1))
              dphidxi(3,2)= 0.25*(1.+xi(1))
              dphidxi(4,2)= 0.25*(1.-xi(1))
              x=  phi(1)*xnodes(1,:)
              x=x+phi(2)*xnodes(2,:)
              x=x+phi(3)*xnodes(3,:)
              x=x+phi(4)*xnodes(4,:)
              t1=   dphidxi(1,1)*xnodes(1,:)
              t1=t1+dphidxi(2,1)*xnodes(2,:)
              t1=t1+dphidxi(3,1)*xnodes(3,:)
              t1=t1+dphidxi(4,1)*xnodes(4,:)
              t2=   dphidxi(1,2)*xnodes(1,:)
              t2=t2+dphidxi(2,2)*xnodes(2,:)
              t2=t2+dphidxi(3,2)*xnodes(3,:)
              t2=t2+dphidxi(4,2)*xnodes(4,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              n=fbem_cross_product(t1,t2)
              jg=sqrt(dot_product(n,n))
              if (alpha.eq.0) then
                ref_integral=ref_integral+log(r)*jg*w
              else
                ref_integral=ref_integral+(1./r**alpha)*jg*w
              end if
            end do
            end do
            !
            ! Calculate the approximate integral solution with ngp
            !
            integral=0.
            do k1=1,gl11_n(ngp)
            do k2=1,gl11_n(ngp)
              xi(1)=gl11_xi(k1,ngp)
              xi(2)=gl11_xi(k2,ngp)
              w=gl11_w(k1,ngp)*gl11_w(k2,ngp)
              phi(1)=0.25*(1.-xi(1))*(1.-xi(2))
              phi(2)=0.25*(1.+xi(1))*(1.-xi(2))
              phi(3)=0.25*(1.+xi(1))*(1.+xi(2))
              phi(4)=0.25*(1.-xi(1))*(1.+xi(2))
              dphidxi(1,1)=-0.25*(1.-xi(2))
              dphidxi(2,1)= 0.25*(1.-xi(2))
              dphidxi(3,1)= 0.25*(1.+xi(2))
              dphidxi(4,1)=-0.25*(1.+xi(2))
              dphidxi(1,2)=-0.25*(1.-xi(1))
              dphidxi(2,2)=-0.25*(1.+xi(1))
              dphidxi(3,2)= 0.25*(1.+xi(1))
              dphidxi(4,2)= 0.25*(1.-xi(1))
              x=  phi(1)*xnodes(1,:)
              x=x+phi(2)*xnodes(2,:)
              x=x+phi(3)*xnodes(3,:)
              x=x+phi(4)*xnodes(4,:)
              t1=   dphidxi(1,1)*xnodes(1,:)
              t1=t1+dphidxi(2,1)*xnodes(2,:)
              t1=t1+dphidxi(3,1)*xnodes(3,:)
              t1=t1+dphidxi(4,1)*xnodes(4,:)
              t2=   dphidxi(1,2)*xnodes(1,:)
              t2=t2+dphidxi(2,2)*xnodes(2,:)
              t2=t2+dphidxi(3,2)*xnodes(3,:)
              t2=t2+dphidxi(4,2)*xnodes(4,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              n=fbem_cross_product(t1,t2)
              jg=sqrt(dot_product(n,n))
              if (alpha.eq.0) then
                integral=integral+log(r)*jg*w
              else
                integral=integral+(1./r**alpha)*jg*w
              end if
            end do
            end do

            !
            ! Stop criterion
            !
            eps=abs((integral-ref_integral)/ref_integral)
            if (eps.gt.10.0d0**e) then
              function_d(ngp)=d
              function_n(ngp)=ngp
              ini_nd_tmp=idistance
              exit
            end if

          end do ! Loop through the distance d

        end do ! Loop though necessary number of gauss points

        ! Write
        do ngp=ngp_max,ngp_min,-delta_ngp
          if (function_n(ngp).ne.0) then
            write(60,fmtstr) alpha, e, 2, function_d(ngp), function_n(ngp)
          end if
        end do
        write(60,*)
        write(60,*)

        ! ------------------------------------------
        ! Collocation point in vertex (out-of-plane)
        ! ------------------------------------------

        write(*,*) '  Collocation point in vertex (out-of-plane)'

        ! Initialize
        function_n=0
        function_d=0.d0
        ini_nd_tmp=nd

        ! Loop though necessary number of gauss points
        do ngp=ngp_min,ngp_max,delta_ngp

          ini_nd=ini_nd_tmp
          ! Loop through the distance d
          do idistance=ini_nd,1,-1
            ! Generate the collocation point
            d=10.0d0**(dble(idistance-1)*delta_d)*d_min
            x_i(1)=1.
            x_i(2)=1.
            x_i(3)=d

            !
            ! Calculate the reference numerical solution with ngp+2
            !
            ref_integral=0.
            do k1=1,gl11_n(ngp+2)
            do k2=1,gl11_n(ngp+2)
              xi(1)=gl11_xi(k1,ngp+2)
              xi(2)=gl11_xi(k2,ngp+2)
              w=gl11_w(k1,ngp+2)*gl11_w(k2,ngp+2)
              phi(1)=0.25*(1.-xi(1))*(1.-xi(2))
              phi(2)=0.25*(1.+xi(1))*(1.-xi(2))
              phi(3)=0.25*(1.+xi(1))*(1.+xi(2))
              phi(4)=0.25*(1.-xi(1))*(1.+xi(2))
              dphidxi(1,1)=-0.25*(1.-xi(2))
              dphidxi(2,1)= 0.25*(1.-xi(2))
              dphidxi(3,1)= 0.25*(1.+xi(2))
              dphidxi(4,1)=-0.25*(1.+xi(2))
              dphidxi(1,2)=-0.25*(1.-xi(1))
              dphidxi(2,2)=-0.25*(1.+xi(1))
              dphidxi(3,2)= 0.25*(1.+xi(1))
              dphidxi(4,2)= 0.25*(1.-xi(1))
              x=  phi(1)*xnodes(1,:)
              x=x+phi(2)*xnodes(2,:)
              x=x+phi(3)*xnodes(3,:)
              x=x+phi(4)*xnodes(4,:)
              t1=   dphidxi(1,1)*xnodes(1,:)
              t1=t1+dphidxi(2,1)*xnodes(2,:)
              t1=t1+dphidxi(3,1)*xnodes(3,:)
              t1=t1+dphidxi(4,1)*xnodes(4,:)
              t2=   dphidxi(1,2)*xnodes(1,:)
              t2=t2+dphidxi(2,2)*xnodes(2,:)
              t2=t2+dphidxi(3,2)*xnodes(3,:)
              t2=t2+dphidxi(4,2)*xnodes(4,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              n=fbem_cross_product(t1,t2)
              jg=sqrt(dot_product(n,n))
              if (alpha.eq.0) then
                ref_integral=ref_integral+log(r)*jg*w
              else
                ref_integral=ref_integral+(1./r**alpha)*jg*w
              end if
            end do
            end do
            !
            ! Calculate the approximate integral solution with ngp
            !
            integral=0.
            do k1=1,gl11_n(ngp)
            do k2=1,gl11_n(ngp)
              xi(1)=gl11_xi(k1,ngp)
              xi(2)=gl11_xi(k2,ngp)
              w=gl11_w(k1,ngp)*gl11_w(k2,ngp)
              phi(1)=0.25*(1.-xi(1))*(1.-xi(2))
              phi(2)=0.25*(1.+xi(1))*(1.-xi(2))
              phi(3)=0.25*(1.+xi(1))*(1.+xi(2))
              phi(4)=0.25*(1.-xi(1))*(1.+xi(2))
              dphidxi(1,1)=-0.25*(1.-xi(2))
              dphidxi(2,1)= 0.25*(1.-xi(2))
              dphidxi(3,1)= 0.25*(1.+xi(2))
              dphidxi(4,1)=-0.25*(1.+xi(2))
              dphidxi(1,2)=-0.25*(1.-xi(1))
              dphidxi(2,2)=-0.25*(1.+xi(1))
              dphidxi(3,2)= 0.25*(1.+xi(1))
              dphidxi(4,2)= 0.25*(1.-xi(1))
              x=  phi(1)*xnodes(1,:)
              x=x+phi(2)*xnodes(2,:)
              x=x+phi(3)*xnodes(3,:)
              x=x+phi(4)*xnodes(4,:)
              t1=   dphidxi(1,1)*xnodes(1,:)
              t1=t1+dphidxi(2,1)*xnodes(2,:)
              t1=t1+dphidxi(3,1)*xnodes(3,:)
              t1=t1+dphidxi(4,1)*xnodes(4,:)
              t2=   dphidxi(1,2)*xnodes(1,:)
              t2=t2+dphidxi(2,2)*xnodes(2,:)
              t2=t2+dphidxi(3,2)*xnodes(3,:)
              t2=t2+dphidxi(4,2)*xnodes(4,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              n=fbem_cross_product(t1,t2)
              jg=sqrt(dot_product(n,n))
              if (alpha.eq.0) then
                integral=integral+log(r)*jg*w
              else
                integral=integral+(1./r**alpha)*jg*w
              end if
            end do
            end do

            !
            ! Stop criterion
            !
            eps=abs((integral-ref_integral)/ref_integral)
            if (eps.gt.10.0d0**e) then
              function_d(ngp)=d
              function_n(ngp)=ngp
              ini_nd_tmp=idistance
              exit
            end if

          end do ! Loop through the distance d

        end do ! Loop though necessary number of gauss points

        ! Write
        do ngp=ngp_max,ngp_min,-delta_ngp
          if (function_n(ngp).ne.0) then
            write(60,fmtstr) alpha, e, 3, function_d(ngp), function_n(ngp)
          end if
        end do
        write(60,*)
        write(60,*)

        ! -------------------------------------------
        ! Collocation point in edge center (in-plane)
        ! -------------------------------------------

        write(*,*) '  Collocation point in edge center (in-plane)'

        ! Initialize
        function_n=0
        function_d=0.d0
        ini_nd_tmp=nd

        ! Loop though necessary number of gauss points
        do ngp=ngp_min,ngp_max,delta_ngp

          ini_nd=ini_nd_tmp
          ! Loop through the distance d
          do idistance=ini_nd,1,-1
            ! Generate the collocation point
            d=10.0d0**(dble(idistance-1)*delta_d)*d_min
            x_i(1)=0.5
            x_i(2)=1.+d
            x_i(3)=0.

            !
            ! Calculate the reference numerical solution with ngp+2
            !
            ref_integral=0.
            do k1=1,gl11_n(ngp+2)
            do k2=1,gl11_n(ngp+2)
              xi(1)=gl11_xi(k1,ngp+2)
              xi(2)=gl11_xi(k2,ngp+2)
              w=gl11_w(k1,ngp+2)*gl11_w(k2,ngp+2)
              phi(1)=0.25*(1.-xi(1))*(1.-xi(2))
              phi(2)=0.25*(1.+xi(1))*(1.-xi(2))
              phi(3)=0.25*(1.+xi(1))*(1.+xi(2))
              phi(4)=0.25*(1.-xi(1))*(1.+xi(2))
              dphidxi(1,1)=-0.25*(1.-xi(2))
              dphidxi(2,1)= 0.25*(1.-xi(2))
              dphidxi(3,1)= 0.25*(1.+xi(2))
              dphidxi(4,1)=-0.25*(1.+xi(2))
              dphidxi(1,2)=-0.25*(1.-xi(1))
              dphidxi(2,2)=-0.25*(1.+xi(1))
              dphidxi(3,2)= 0.25*(1.+xi(1))
              dphidxi(4,2)= 0.25*(1.-xi(1))
              x=  phi(1)*xnodes(1,:)
              x=x+phi(2)*xnodes(2,:)
              x=x+phi(3)*xnodes(3,:)
              x=x+phi(4)*xnodes(4,:)
              t1=   dphidxi(1,1)*xnodes(1,:)
              t1=t1+dphidxi(2,1)*xnodes(2,:)
              t1=t1+dphidxi(3,1)*xnodes(3,:)
              t1=t1+dphidxi(4,1)*xnodes(4,:)
              t2=   dphidxi(1,2)*xnodes(1,:)
              t2=t2+dphidxi(2,2)*xnodes(2,:)
              t2=t2+dphidxi(3,2)*xnodes(3,:)
              t2=t2+dphidxi(4,2)*xnodes(4,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              n=fbem_cross_product(t1,t2)
              jg=sqrt(dot_product(n,n))
              if (alpha.eq.0) then
                ref_integral=ref_integral+log(r)*jg*w
              else
                ref_integral=ref_integral+(1./r**alpha)*jg*w
              end if
            end do
            end do
            !
            ! Calculate the approximate integral solution with ngp
            !
            integral=0.
            do k1=1,gl11_n(ngp)
            do k2=1,gl11_n(ngp)
              xi(1)=gl11_xi(k1,ngp)
              xi(2)=gl11_xi(k2,ngp)
              w=gl11_w(k1,ngp)*gl11_w(k2,ngp)
              phi(1)=0.25*(1.-xi(1))*(1.-xi(2))
              phi(2)=0.25*(1.+xi(1))*(1.-xi(2))
              phi(3)=0.25*(1.+xi(1))*(1.+xi(2))
              phi(4)=0.25*(1.-xi(1))*(1.+xi(2))
              dphidxi(1,1)=-0.25*(1.-xi(2))
              dphidxi(2,1)= 0.25*(1.-xi(2))
              dphidxi(3,1)= 0.25*(1.+xi(2))
              dphidxi(4,1)=-0.25*(1.+xi(2))
              dphidxi(1,2)=-0.25*(1.-xi(1))
              dphidxi(2,2)=-0.25*(1.+xi(1))
              dphidxi(3,2)= 0.25*(1.+xi(1))
              dphidxi(4,2)= 0.25*(1.-xi(1))
              x=  phi(1)*xnodes(1,:)
              x=x+phi(2)*xnodes(2,:)
              x=x+phi(3)*xnodes(3,:)
              x=x+phi(4)*xnodes(4,:)
              t1=   dphidxi(1,1)*xnodes(1,:)
              t1=t1+dphidxi(2,1)*xnodes(2,:)
              t1=t1+dphidxi(3,1)*xnodes(3,:)
              t1=t1+dphidxi(4,1)*xnodes(4,:)
              t2=   dphidxi(1,2)*xnodes(1,:)
              t2=t2+dphidxi(2,2)*xnodes(2,:)
              t2=t2+dphidxi(3,2)*xnodes(3,:)
              t2=t2+dphidxi(4,2)*xnodes(4,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              n=fbem_cross_product(t1,t2)
              jg=sqrt(dot_product(n,n))
              if (alpha.eq.0) then
                integral=integral+log(r)*jg*w
              else
                integral=integral+(1./r**alpha)*jg*w
              end if
            end do
            end do

            !
            ! Stop criterion
            !
            eps=abs((integral-ref_integral)/ref_integral)
            if (eps.gt.10.0d0**e) then
              function_d(ngp)=d
              function_n(ngp)=ngp
              ini_nd_tmp=idistance
              exit
            end if

          end do ! Loop through the distance d

        end do ! Loop though necessary number of gauss points

        ! Write
        do ngp=ngp_max,ngp_min,-delta_ngp
          if (function_n(ngp).ne.0) then
            write(60,fmtstr) alpha, e, 4, function_d(ngp), function_n(ngp)
          end if
        end do
        write(60,*)
        write(60,*)

        ! --------------------------------------
        ! Collocation point in vertex (in-plane)
        ! --------------------------------------

        write(*,*) '  Collocation point in vertex (in-plane)'

        ! Initialize
        function_n=0
        function_d=0.d0
        ini_nd_tmp=nd

        ! Loop though necessary number of gauss points
        do ngp=ngp_min,ngp_max,delta_ngp

          ini_nd=ini_nd_tmp
          ! Loop through the distance d
          do idistance=ini_nd,1,-1
            ! Generate the collocation point
            d=10.0d0**(dble(idistance-1)*delta_d)*d_min
            x_i(1)=1.
            x_i(2)=1.+d
            x_i(3)=0.

            !
            ! Calculate the reference numerical solution with ngp+2
            !
            ref_integral=0.
            do k1=1,gl11_n(ngp+2)
            do k2=1,gl11_n(ngp+2)
              xi(1)=gl11_xi(k1,ngp+2)
              xi(2)=gl11_xi(k2,ngp+2)
              w=gl11_w(k1,ngp+2)*gl11_w(k2,ngp+2)
              phi(1)=0.25*(1.-xi(1))*(1.-xi(2))
              phi(2)=0.25*(1.+xi(1))*(1.-xi(2))
              phi(3)=0.25*(1.+xi(1))*(1.+xi(2))
              phi(4)=0.25*(1.-xi(1))*(1.+xi(2))
              dphidxi(1,1)=-0.25*(1.-xi(2))
              dphidxi(2,1)= 0.25*(1.-xi(2))
              dphidxi(3,1)= 0.25*(1.+xi(2))
              dphidxi(4,1)=-0.25*(1.+xi(2))
              dphidxi(1,2)=-0.25*(1.-xi(1))
              dphidxi(2,2)=-0.25*(1.+xi(1))
              dphidxi(3,2)= 0.25*(1.+xi(1))
              dphidxi(4,2)= 0.25*(1.-xi(1))
              x=  phi(1)*xnodes(1,:)
              x=x+phi(2)*xnodes(2,:)
              x=x+phi(3)*xnodes(3,:)
              x=x+phi(4)*xnodes(4,:)
              t1=   dphidxi(1,1)*xnodes(1,:)
              t1=t1+dphidxi(2,1)*xnodes(2,:)
              t1=t1+dphidxi(3,1)*xnodes(3,:)
              t1=t1+dphidxi(4,1)*xnodes(4,:)
              t2=   dphidxi(1,2)*xnodes(1,:)
              t2=t2+dphidxi(2,2)*xnodes(2,:)
              t2=t2+dphidxi(3,2)*xnodes(3,:)
              t2=t2+dphidxi(4,2)*xnodes(4,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              n=fbem_cross_product(t1,t2)
              jg=sqrt(dot_product(n,n))
              if (alpha.eq.0) then
                ref_integral=ref_integral+log(r)*jg*w
              else
                ref_integral=ref_integral+(1./r**alpha)*jg*w
              end if
            end do
            end do
            !
            ! Calculate the approximate integral solution with ngp
            !
            integral=0.
            do k1=1,gl11_n(ngp)
            do k2=1,gl11_n(ngp)
              xi(1)=gl11_xi(k1,ngp)
              xi(2)=gl11_xi(k2,ngp)
              w=gl11_w(k1,ngp)*gl11_w(k2,ngp)
              phi(1)=0.25*(1.-xi(1))*(1.-xi(2))
              phi(2)=0.25*(1.+xi(1))*(1.-xi(2))
              phi(3)=0.25*(1.+xi(1))*(1.+xi(2))
              phi(4)=0.25*(1.-xi(1))*(1.+xi(2))
              dphidxi(1,1)=-0.25*(1.-xi(2))
              dphidxi(2,1)= 0.25*(1.-xi(2))
              dphidxi(3,1)= 0.25*(1.+xi(2))
              dphidxi(4,1)=-0.25*(1.+xi(2))
              dphidxi(1,2)=-0.25*(1.-xi(1))
              dphidxi(2,2)=-0.25*(1.+xi(1))
              dphidxi(3,2)= 0.25*(1.+xi(1))
              dphidxi(4,2)= 0.25*(1.-xi(1))
              x=  phi(1)*xnodes(1,:)
              x=x+phi(2)*xnodes(2,:)
              x=x+phi(3)*xnodes(3,:)
              x=x+phi(4)*xnodes(4,:)
              t1=   dphidxi(1,1)*xnodes(1,:)
              t1=t1+dphidxi(2,1)*xnodes(2,:)
              t1=t1+dphidxi(3,1)*xnodes(3,:)
              t1=t1+dphidxi(4,1)*xnodes(4,:)
              t2=   dphidxi(1,2)*xnodes(1,:)
              t2=t2+dphidxi(2,2)*xnodes(2,:)
              t2=t2+dphidxi(3,2)*xnodes(3,:)
              t2=t2+dphidxi(4,2)*xnodes(4,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              n=fbem_cross_product(t1,t2)
              jg=sqrt(dot_product(n,n))
              if (alpha.eq.0) then
                integral=integral+log(r)*jg*w
              else
                integral=integral+(1./r**alpha)*jg*w
              end if
            end do
            end do

            !
            ! Stop criterion
            !
            eps=abs((integral-ref_integral)/ref_integral)
            if (eps.gt.10.0d0**e) then
              function_d(ngp)=d
              function_n(ngp)=ngp
              ini_nd_tmp=idistance
              exit
            end if

          end do ! Loop through the distance d

        end do ! Loop though necessary number of gauss points

        ! Write
        do ngp=ngp_max,ngp_min,-delta_ngp
          if (function_n(ngp).ne.0) then
            write(60,fmtstr) alpha, e, 5, function_d(ngp), function_n(ngp)
          end if
        end do
        write(60,*)
        write(60,*)

      end do ! Loop through error heights

    end do ! Loop through alpha

    close(60)

  end subroutine fbem_qs_test_standard_quad







  subroutine fbem_qs_test_standard_tri

    implicit none

    real(kind=real128) :: xnodes(3,3), xi(2), x_i(3), w, x(3), t1(3), t2(3), n(3), rv(3), r, jg
    real(kind=real128) :: phi(3), dphidxi(3,2)
    real(kind=real128) :: xip(2)
    real(kind=real128) :: eps
    real(kind=real128) :: d, d_min, d_max, delta_d
    integer            :: e_min, delta_e, e_max
    integer            :: i, j, k1, k2, e, npoints
    real(kind=real128) :: ref_integral, integral, pdensity
    integer            :: alpha
    integer            :: nd
    integer            :: idistance, ini_nd, ini_nd_tmp
    integer            :: ngp_min, ngp_max, delta_ngp, ngp
    integer, allocatable :: function_n(:)
    real(kind=real128), allocatable :: function_d(:)
    character(len=32)  :: fmtstr

    fmtstr='(3i11,e25.16,i11)'

    ! Reference element (quad4)
    xnodes=0.
    xnodes(1,1)=1.
    xnodes(1,2)=0.
    xnodes(2,1)=0.5
    xnodes(2,2)=0.5*sqrt(3.d0)
    xnodes(3,1)=0.
    xnodes(3,2)=0.

    ! Distance discretization
    d_min=1.d-6
    d_max=1.d6
    delta_d=0.001d0
    nd=nint((log10(d_max)-log10(d_min))/delta_d)+1

    ! Error height discretization
    e_min=-15
    e_max=-3
    delta_e=1

    ! Number of Gauss Points used
    ngp_min=2
    ngp_max=30
    delta_ngp=1
    allocate (function_n(ngp_min:ngp_max))
    allocate (function_d(ngp_min:ngp_max))

    ! File configuration
    open(unit=60,file='n_d_standard_tri.dat',recl=1024)

    ! Header
    write(60,*) '# alpha, log10(eps), ch. path, d, N(d)'
    write(60,*) '# alpha: integrand 1/r^alpha'
    write(60,*) '# ch. path (characteristic path) for tri element with vertices at (0,0), (1,0) and (0.5,0.5*sqrt(3)):'
    write(60,*) '#   1 x=0.5 , y=sqrt(3)/6 , z=d'
    write(60,*) '#   2 x=0.5 , y=0         , z=d'
    write(60,*) '#   3 x=1   , y=0         , z=d'
    write(60,*) '#   4 x=0.5 , y=-d        , z=0'
    write(60,*) '#   5 x=1   , y=-d        , z=0'

    ! Loop through alpha
    do alpha=1,9
      write(*,*) 'integrand type = ', alpha

      ! Loop through error heights
      do e=e_min,e_max,delta_e
        write(*,*) ' log10(relative_error) = ', e

        ! ------------------------------------------
        ! Collocation point in center (out-of-plane)
        ! ------------------------------------------

        write(*,*) '  Collocation point in center (out-of-plane)'

        ! Initialize
        function_n=0
        function_d=0.d0
        ini_nd_tmp=nd

        ! Loop though necessary number of gauss points
        do ngp=ngp_min,ngp_max,delta_ngp

          ini_nd=ini_nd_tmp
          ! Loop through the distance d
          do idistance=ini_nd,1,-1
            ! Generate the collocation point
            d=10.0d0**(dble(idistance-1)*delta_d)*d_min
            x_i(1)=0.5d0
            x_i(2)=sqrt(3.d0)/6.d0
            x_i(3)=d

            !
            ! Calculate the reference numerical solution with ngp+2
            !
            ref_integral=0.
            do k1=1,gl01_n(ngp+2)
            do k2=1,gl01_n(ngp+2)
              xip(1)=gl01_xi(k1,ngp+2)
              xip(2)=gl01_xi(k2,ngp+2)
              xi(1)=(1.d0-xip(2))*xip(1)
              xi(2)=xip(2)
              w=gl01_w(k1,ngp+2)*gl01_w(k2,ngp+2)*(1.d0-xip(2))
              phi(1)=xi(1)
              phi(2)=xi(2)
              phi(3)=1.-xi(1)-xi(2)
              dphidxi(1,1)= 1.
              dphidxi(2,1)= 0.
              dphidxi(3,1)=-1.
              dphidxi(1,2)= 0.
              dphidxi(2,2)= 1.
              dphidxi(3,2)=-1.
              x=  phi(1)*xnodes(1,:)
              x=x+phi(2)*xnodes(2,:)
              x=x+phi(3)*xnodes(3,:)
              t1=   dphidxi(1,1)*xnodes(1,:)
              t1=t1+dphidxi(2,1)*xnodes(2,:)
              t1=t1+dphidxi(3,1)*xnodes(3,:)
              t2=   dphidxi(1,2)*xnodes(1,:)
              t2=t2+dphidxi(2,2)*xnodes(2,:)
              t2=t2+dphidxi(3,2)*xnodes(3,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              n=fbem_cross_product(t1,t2)
              jg=sqrt(dot_product(n,n))
              if (alpha.eq.0) then
                ref_integral=ref_integral+log(r)*jg*w
              else
                ref_integral=ref_integral+(1./r**alpha)*jg*w
              end if
            end do
            end do
            !
            ! Calculate the approximate integral solution with ngp
            !
            integral=0.
            do k1=1,gl01_n(ngp)
            do k2=1,gl01_n(ngp)
              xip(1)=gl01_xi(k1,ngp)
              xip(2)=gl01_xi(k2,ngp)
              xi(1)=(1.d0-xip(2))*xip(1)
              xi(2)=xip(2)
              w=gl01_w(k1,ngp)*gl01_w(k2,ngp)*(1.d0-xip(2))
              phi(1)=xi(1)
              phi(2)=xi(2)
              phi(3)=1.-xi(1)-xi(2)
              dphidxi(1,1)= 1.
              dphidxi(2,1)= 0.
              dphidxi(3,1)=-1.
              dphidxi(1,2)= 0.
              dphidxi(2,2)= 1.
              dphidxi(3,2)=-1.
              x=  phi(1)*xnodes(1,:)
              x=x+phi(2)*xnodes(2,:)
              x=x+phi(3)*xnodes(3,:)
              t1=   dphidxi(1,1)*xnodes(1,:)
              t1=t1+dphidxi(2,1)*xnodes(2,:)
              t1=t1+dphidxi(3,1)*xnodes(3,:)
              t2=   dphidxi(1,2)*xnodes(1,:)
              t2=t2+dphidxi(2,2)*xnodes(2,:)
              t2=t2+dphidxi(3,2)*xnodes(3,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              n=fbem_cross_product(t1,t2)
              jg=sqrt(dot_product(n,n))
              if (alpha.eq.0) then
                integral=integral+log(r)*jg*w
              else
                integral=integral+(1./r**alpha)*jg*w
              end if
            end do
            end do

            !
            ! Stop criterion
            !
            eps=abs((integral-ref_integral)/ref_integral)
            if (eps.gt.10.0d0**e) then
              function_d(ngp)=d
              function_n(ngp)=ngp
              ini_nd_tmp=idistance
              exit
            end if

          end do ! Loop through the distance d

        end do ! Loop though necessary number of gauss points

        ! Write
        do ngp=ngp_max,ngp_min,-delta_ngp
          if (function_n(ngp).ne.0) then
            write(60,fmtstr) alpha, e, 1, function_d(ngp), function_n(ngp)
          end if
        end do
        write(60,*)
        write(60,*)

        ! -----------------------------------------------
        ! Collocation point in edge center (out-of-plane)
        ! -----------------------------------------------

        write(*,*) '  Collocation point in edge center (out-of-plane)'

        ! Initialize
        function_n=0
        function_d=0.d0
        ini_nd_tmp=nd

        ! Loop though necessary number of gauss points
        do ngp=ngp_min,ngp_max,delta_ngp

          ini_nd=ini_nd_tmp
          ! Loop through the distance d
          do idistance=ini_nd,1,-1
            ! Generate the collocation point
            d=10.0d0**(dble(idistance-1)*delta_d)*d_min
            x_i(1)=0.5d0
            x_i(2)=0.d0
            x_i(3)=d

            !
            ! Calculate the reference numerical solution with ngp+2
            !
            ref_integral=0.
            do k1=1,gl01_n(ngp+2)
            do k2=1,gl01_n(ngp+2)
              xip(1)=gl01_xi(k1,ngp+2)
              xip(2)=gl01_xi(k2,ngp+2)
              xi(1)=(1.d0-xip(2))*xip(1)
              xi(2)=xip(2)
              w=gl01_w(k1,ngp+2)*gl01_w(k2,ngp+2)*(1.d0-xip(2))
              phi(1)=xi(1)
              phi(2)=xi(2)
              phi(3)=1.-xi(1)-xi(2)
              dphidxi(1,1)= 1.
              dphidxi(2,1)= 0.
              dphidxi(3,1)=-1.
              dphidxi(1,2)= 0.
              dphidxi(2,2)= 1.
              dphidxi(3,2)=-1.
              x=  phi(1)*xnodes(1,:)
              x=x+phi(2)*xnodes(2,:)
              x=x+phi(3)*xnodes(3,:)
              t1=   dphidxi(1,1)*xnodes(1,:)
              t1=t1+dphidxi(2,1)*xnodes(2,:)
              t1=t1+dphidxi(3,1)*xnodes(3,:)
              t2=   dphidxi(1,2)*xnodes(1,:)
              t2=t2+dphidxi(2,2)*xnodes(2,:)
              t2=t2+dphidxi(3,2)*xnodes(3,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              n=fbem_cross_product(t1,t2)
              jg=sqrt(dot_product(n,n))
              if (alpha.eq.0) then
                ref_integral=ref_integral+log(r)*jg*w
              else
                ref_integral=ref_integral+(1./r**alpha)*jg*w
              end if
            end do
            end do
            !
            ! Calculate the approximate integral solution with ngp
            !
            integral=0.
            do k1=1,gl01_n(ngp)
            do k2=1,gl01_n(ngp)
              xip(1)=gl01_xi(k1,ngp)
              xip(2)=gl01_xi(k2,ngp)
              xi(1)=(1.d0-xip(2))*xip(1)
              xi(2)=xip(2)
              w=gl01_w(k1,ngp)*gl01_w(k2,ngp)*(1.d0-xip(2))
              phi(1)=xi(1)
              phi(2)=xi(2)
              phi(3)=1.-xi(1)-xi(2)
              dphidxi(1,1)= 1.
              dphidxi(2,1)= 0.
              dphidxi(3,1)=-1.
              dphidxi(1,2)= 0.
              dphidxi(2,2)= 1.
              dphidxi(3,2)=-1.
              x=  phi(1)*xnodes(1,:)
              x=x+phi(2)*xnodes(2,:)
              x=x+phi(3)*xnodes(3,:)
              t1=   dphidxi(1,1)*xnodes(1,:)
              t1=t1+dphidxi(2,1)*xnodes(2,:)
              t1=t1+dphidxi(3,1)*xnodes(3,:)
              t2=   dphidxi(1,2)*xnodes(1,:)
              t2=t2+dphidxi(2,2)*xnodes(2,:)
              t2=t2+dphidxi(3,2)*xnodes(3,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              n=fbem_cross_product(t1,t2)
              jg=sqrt(dot_product(n,n))
              if (alpha.eq.0) then
                integral=integral+log(r)*jg*w
              else
                integral=integral+(1./r**alpha)*jg*w
              end if
            end do
            end do

            !
            ! Stop criterion
            !
            eps=abs((integral-ref_integral)/ref_integral)
            if (eps.gt.10.0d0**e) then
              function_d(ngp)=d
              function_n(ngp)=ngp
              ini_nd_tmp=idistance
              exit
            end if

          end do ! Loop through the distance d

        end do ! Loop though necessary number of gauss points

        ! Write
        do ngp=ngp_max,ngp_min,-delta_ngp
          if (function_n(ngp).ne.0) then
            write(60,fmtstr) alpha, e, 2, function_d(ngp), function_n(ngp)
          end if
        end do
        write(60,*)
        write(60,*)

        ! ------------------------------------------
        ! Collocation point in vertex (out-of-plane)
        ! ------------------------------------------

        write(*,*) '  Collocation point in vertex (out-of-plane)'

        ! Initialize
        function_n=0
        function_d=0.d0
        ini_nd_tmp=nd

        ! Loop though necessary number of gauss points
        do ngp=ngp_min,ngp_max,delta_ngp

          ini_nd=ini_nd_tmp
          ! Loop through the distance d
          do idistance=ini_nd,1,-1
            ! Generate the collocation point
            d=10.0d0**(dble(idistance-1)*delta_d)*d_min
            x_i(1)=1.d0
            x_i(2)=0.d0
            x_i(3)=d

            !
            ! Calculate the reference numerical solution with ngp+2
            !
            ref_integral=0.
            do k1=1,gl01_n(ngp+2)
            do k2=1,gl01_n(ngp+2)
              xip(1)=gl01_xi(k1,ngp+2)
              xip(2)=gl01_xi(k2,ngp+2)
              xi(1)=(1.d0-xip(2))*xip(1)
              xi(2)=xip(2)
              w=gl01_w(k1,ngp+2)*gl01_w(k2,ngp+2)*(1.d0-xip(2))
              phi(1)=xi(1)
              phi(2)=xi(2)
              phi(3)=1.-xi(1)-xi(2)
              dphidxi(1,1)= 1.
              dphidxi(2,1)= 0.
              dphidxi(3,1)=-1.
              dphidxi(1,2)= 0.
              dphidxi(2,2)= 1.
              dphidxi(3,2)=-1.
              x=  phi(1)*xnodes(1,:)
              x=x+phi(2)*xnodes(2,:)
              x=x+phi(3)*xnodes(3,:)
              t1=   dphidxi(1,1)*xnodes(1,:)
              t1=t1+dphidxi(2,1)*xnodes(2,:)
              t1=t1+dphidxi(3,1)*xnodes(3,:)
              t2=   dphidxi(1,2)*xnodes(1,:)
              t2=t2+dphidxi(2,2)*xnodes(2,:)
              t2=t2+dphidxi(3,2)*xnodes(3,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              n=fbem_cross_product(t1,t2)
              jg=sqrt(dot_product(n,n))
              if (alpha.eq.0) then
                ref_integral=ref_integral+log(r)*jg*w
              else
                ref_integral=ref_integral+(1./r**alpha)*jg*w
              end if
            end do
            end do
            !
            ! Calculate the approximate integral solution with ngp
            !
            integral=0.
            do k1=1,gl01_n(ngp)
            do k2=1,gl01_n(ngp)
              xip(1)=gl01_xi(k1,ngp)
              xip(2)=gl01_xi(k2,ngp)
              xi(1)=(1.d0-xip(2))*xip(1)
              xi(2)=xip(2)
              w=gl01_w(k1,ngp)*gl01_w(k2,ngp)*(1.d0-xip(2))
              phi(1)=xi(1)
              phi(2)=xi(2)
              phi(3)=1.-xi(1)-xi(2)
              dphidxi(1,1)= 1.
              dphidxi(2,1)= 0.
              dphidxi(3,1)=-1.
              dphidxi(1,2)= 0.
              dphidxi(2,2)= 1.
              dphidxi(3,2)=-1.
              x=  phi(1)*xnodes(1,:)
              x=x+phi(2)*xnodes(2,:)
              x=x+phi(3)*xnodes(3,:)
              t1=   dphidxi(1,1)*xnodes(1,:)
              t1=t1+dphidxi(2,1)*xnodes(2,:)
              t1=t1+dphidxi(3,1)*xnodes(3,:)
              t2=   dphidxi(1,2)*xnodes(1,:)
              t2=t2+dphidxi(2,2)*xnodes(2,:)
              t2=t2+dphidxi(3,2)*xnodes(3,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              n=fbem_cross_product(t1,t2)
              jg=sqrt(dot_product(n,n))
              if (alpha.eq.0) then
                integral=integral+log(r)*jg*w
              else
                integral=integral+(1./r**alpha)*jg*w
              end if
            end do
            end do

            !
            ! Stop criterion
            !
            eps=abs((integral-ref_integral)/ref_integral)
            if (eps.gt.10.0d0**e) then
              function_d(ngp)=d
              function_n(ngp)=ngp
              ini_nd_tmp=idistance
              exit
            end if

          end do ! Loop through the distance d

        end do ! Loop though necessary number of gauss points

        ! Write
        do ngp=ngp_max,ngp_min,-delta_ngp
          if (function_n(ngp).ne.0) then
            write(60,fmtstr) alpha, e, 3, function_d(ngp), function_n(ngp)
          end if
        end do
        write(60,*)
        write(60,*)

        ! -------------------------------------------
        ! Collocation point in edge center (in-plane)
        ! -------------------------------------------

        write(*,*) '  Collocation point in edge center (in-plane)'

        ! Initialize
        function_n=0
        function_d=0.d0
        ini_nd_tmp=nd

        ! Loop though necessary number of gauss points
        do ngp=ngp_min,ngp_max,delta_ngp

          ini_nd=ini_nd_tmp
          ! Loop through the distance d
          do idistance=ini_nd,1,-1
            ! Generate the collocation point
            d=10.0d0**(dble(idistance-1)*delta_d)*d_min
            x_i(1)=0.5d0
            x_i(2)=-d
            x_i(3)=0.

            !
            ! Calculate the reference numerical solution with ngp+2
            !
            ref_integral=0.
            do k1=1,gl01_n(ngp+2)
            do k2=1,gl01_n(ngp+2)
              xip(1)=gl01_xi(k1,ngp+2)
              xip(2)=gl01_xi(k2,ngp+2)
              xi(1)=(1.d0-xip(2))*xip(1)
              xi(2)=xip(2)
              w=gl01_w(k1,ngp+2)*gl01_w(k2,ngp+2)*(1.d0-xip(2))
              phi(1)=xi(1)
              phi(2)=xi(2)
              phi(3)=1.-xi(1)-xi(2)
              dphidxi(1,1)= 1.
              dphidxi(2,1)= 0.
              dphidxi(3,1)=-1.
              dphidxi(1,2)= 0.
              dphidxi(2,2)= 1.
              dphidxi(3,2)=-1.
              x=  phi(1)*xnodes(1,:)
              x=x+phi(2)*xnodes(2,:)
              x=x+phi(3)*xnodes(3,:)
              t1=   dphidxi(1,1)*xnodes(1,:)
              t1=t1+dphidxi(2,1)*xnodes(2,:)
              t1=t1+dphidxi(3,1)*xnodes(3,:)
              t2=   dphidxi(1,2)*xnodes(1,:)
              t2=t2+dphidxi(2,2)*xnodes(2,:)
              t2=t2+dphidxi(3,2)*xnodes(3,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              n=fbem_cross_product(t1,t2)
              jg=sqrt(dot_product(n,n))
              if (alpha.eq.0) then
                ref_integral=ref_integral+log(r)*jg*w
              else
                ref_integral=ref_integral+(1./r**alpha)*jg*w
              end if
            end do
            end do
            !
            ! Calculate the approximate integral solution with ngp
            !
            integral=0.
            do k1=1,gl01_n(ngp)
            do k2=1,gl01_n(ngp)
              xip(1)=gl01_xi(k1,ngp)
              xip(2)=gl01_xi(k2,ngp)
              xi(1)=(1.d0-xip(2))*xip(1)
              xi(2)=xip(2)
              w=gl01_w(k1,ngp)*gl01_w(k2,ngp)*(1.d0-xip(2))
              phi(1)=xi(1)
              phi(2)=xi(2)
              phi(3)=1.-xi(1)-xi(2)
              dphidxi(1,1)= 1.
              dphidxi(2,1)= 0.
              dphidxi(3,1)=-1.
              dphidxi(1,2)= 0.
              dphidxi(2,2)= 1.
              dphidxi(3,2)=-1.
              x=  phi(1)*xnodes(1,:)
              x=x+phi(2)*xnodes(2,:)
              x=x+phi(3)*xnodes(3,:)
              t1=   dphidxi(1,1)*xnodes(1,:)
              t1=t1+dphidxi(2,1)*xnodes(2,:)
              t1=t1+dphidxi(3,1)*xnodes(3,:)
              t2=   dphidxi(1,2)*xnodes(1,:)
              t2=t2+dphidxi(2,2)*xnodes(2,:)
              t2=t2+dphidxi(3,2)*xnodes(3,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              n=fbem_cross_product(t1,t2)
              jg=sqrt(dot_product(n,n))
              if (alpha.eq.0) then
                integral=integral+log(r)*jg*w
              else
                integral=integral+(1./r**alpha)*jg*w
              end if
            end do
            end do

            !
            ! Stop criterion
            !
            eps=abs((integral-ref_integral)/ref_integral)
            if (eps.gt.10.0d0**e) then
              function_d(ngp)=d
              function_n(ngp)=ngp
              ini_nd_tmp=idistance
              exit
            end if

          end do ! Loop through the distance d

        end do ! Loop though necessary number of gauss points

        ! Write
        do ngp=ngp_max,ngp_min,-delta_ngp
          if (function_n(ngp).ne.0) then
            write(60,fmtstr) alpha, e, 4, function_d(ngp), function_n(ngp)
          end if
        end do
        write(60,*)
        write(60,*)

        ! --------------------------------------
        ! Collocation point in vertex (in-plane)
        ! --------------------------------------

        write(*,*) '  Collocation point in vertex (in-plane)'

        ! Initialize
        function_n=0
        function_d=0.d0
        ini_nd_tmp=nd

        ! Loop though necessary number of gauss points
        do ngp=ngp_min,ngp_max,delta_ngp

          ini_nd=ini_nd_tmp
          ! Loop through the distance d
          do idistance=ini_nd,1,-1
            ! Generate the collocation point
            d=10.0d0**(dble(idistance-1)*delta_d)*d_min
            x_i(1)=1.d0
            x_i(2)=-d
            x_i(3)=0.

            !
            ! Calculate the reference numerical solution with ngp+2
            !
            ref_integral=0.
            do k1=1,gl01_n(ngp+2)
            do k2=1,gl01_n(ngp+2)
              xip(1)=gl01_xi(k1,ngp+2)
              xip(2)=gl01_xi(k2,ngp+2)
              xi(1)=(1.d0-xip(2))*xip(1)
              xi(2)=xip(2)
              w=gl01_w(k1,ngp+2)*gl01_w(k2,ngp+2)*(1.d0-xip(2))
              phi(1)=xi(1)
              phi(2)=xi(2)
              phi(3)=1.-xi(1)-xi(2)
              dphidxi(1,1)= 1.
              dphidxi(2,1)= 0.
              dphidxi(3,1)=-1.
              dphidxi(1,2)= 0.
              dphidxi(2,2)= 1.
              dphidxi(3,2)=-1.
              x=  phi(1)*xnodes(1,:)
              x=x+phi(2)*xnodes(2,:)
              x=x+phi(3)*xnodes(3,:)
              t1=   dphidxi(1,1)*xnodes(1,:)
              t1=t1+dphidxi(2,1)*xnodes(2,:)
              t1=t1+dphidxi(3,1)*xnodes(3,:)
              t2=   dphidxi(1,2)*xnodes(1,:)
              t2=t2+dphidxi(2,2)*xnodes(2,:)
              t2=t2+dphidxi(3,2)*xnodes(3,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              n=fbem_cross_product(t1,t2)
              jg=sqrt(dot_product(n,n))
              if (alpha.eq.0) then
                ref_integral=ref_integral+log(r)*jg*w
              else
                ref_integral=ref_integral+(1./r**alpha)*jg*w
              end if
            end do
            end do
            !
            ! Calculate the approximate integral solution with ngp
            !
            integral=0.
            do k1=1,gl01_n(ngp)
            do k2=1,gl01_n(ngp)
              xip(1)=gl01_xi(k1,ngp)
              xip(2)=gl01_xi(k2,ngp)
              xi(1)=(1.d0-xip(2))*xip(1)
              xi(2)=xip(2)
              w=gl01_w(k1,ngp)*gl01_w(k2,ngp)*(1.d0-xip(2))
              phi(1)=xi(1)
              phi(2)=xi(2)
              phi(3)=1.-xi(1)-xi(2)
              dphidxi(1,1)= 1.
              dphidxi(2,1)= 0.
              dphidxi(3,1)=-1.
              dphidxi(1,2)= 0.
              dphidxi(2,2)= 1.
              dphidxi(3,2)=-1.
              x=  phi(1)*xnodes(1,:)
              x=x+phi(2)*xnodes(2,:)
              x=x+phi(3)*xnodes(3,:)
              t1=   dphidxi(1,1)*xnodes(1,:)
              t1=t1+dphidxi(2,1)*xnodes(2,:)
              t1=t1+dphidxi(3,1)*xnodes(3,:)
              t2=   dphidxi(1,2)*xnodes(1,:)
              t2=t2+dphidxi(2,2)*xnodes(2,:)
              t2=t2+dphidxi(3,2)*xnodes(3,:)
              rv=x-x_i
              r=sqrt(dot_product(rv,rv))
              n=fbem_cross_product(t1,t2)
              jg=sqrt(dot_product(n,n))
              if (alpha.eq.0) then
                integral=integral+log(r)*jg*w
              else
                integral=integral+(1./r**alpha)*jg*w
              end if
            end do
            end do

            !
            ! Stop criterion
            !
            eps=abs((integral-ref_integral)/ref_integral)
            if (eps.gt.10.0d0**e) then
              function_d(ngp)=d
              function_n(ngp)=ngp
              ini_nd_tmp=idistance
              exit
            end if

          end do ! Loop through the distance d

        end do ! Loop though necessary number of gauss points

        ! Write
        do ngp=ngp_max,ngp_min,-delta_ngp
          if (function_n(ngp).ne.0) then
            write(60,fmtstr) alpha, e, 5, function_d(ngp), function_n(ngp)
          end if
        end do
        write(60,*)
        write(60,*)

      end do ! Loop through error heights

    end do ! Loop through alpha

    close(60)

  end subroutine fbem_qs_test_standard_tri

end module fbem_quasisingular_integration
