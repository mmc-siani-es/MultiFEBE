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
!! <b>This module has numerical constants, functions and subroutines.</b>
module fbem_numerical

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all are private
  private

  ! Numerical routines
  public :: fbem_zarg
  public :: fbem_solve_quadratic_equation
  public :: fbem_invert_2x2_matrix
  public :: fbem_invert_3x3_matrix
  public :: fbem_invert_4x4_matrix
  public :: fbem_cross_product
  public :: fbem_vector_normalization
  public :: fbem_vector_norm
  public :: f_poly_1d_linear, f_poly_1d_quadratic, f_poly_2d_linear, f_poly_2d_quadratic
  public :: fbem_modified_bessel_K0_and_K1
  public :: fbem_modified_bessel_K0r_K1r_K2r
  public :: fbem_modified_bessel_K0r_K1r_K2r_2
  public :: fbem_BesselKnR_decomposed
  public :: fbem_decomposed_zexp
  public :: fbem_decomposed_zexp_6
  public :: fbem_zexp_decomposed
  public :: fbem_quicksort
  public :: fbem_quicksort_matrix, fbem_quicksort_matrix_r64
  public :: fbem_binary_tree
  public :: fbem_binary_tree_matrix2
  public :: fbem_unique
  public :: fbem_solve_lse_c
  public :: dsyevj3

  ! ================================================================================================================================
  !type fbem_table
  !  integer                        :: nc, nr
  !  real(kind=real64), allocatable :: cell(:,:)
  !end type fbem_table
  ! ================================================================================================================================

  ! Numerical constants
  real(kind=real64),    parameter, public :: c_e       =2.71828182845904523536028747135d0 !! <tt>e</tt> number
  real(kind=real64),    parameter, public :: c_log2e   =1.44269504088896340735992468100d0 !! <tt>log2(e)</tt>
  real(kind=real64),    parameter, public :: c_log10e  =0.43429448190325182765112891892d0 !! <tt>log10(e)</tt>
  real(kind=real64),    parameter, public :: c_sqrt2   =1.41421356237309504880168872421d0 !! <tt>sqrt(2)</tt>
  real(kind=real64),    parameter, public :: c_sqrt1_2 =0.70710678118654752440084436210d0 !! <tt>sqrt(1/2)</tt>
  real(kind=real64),    parameter, public :: c_sqrt3   =1.73205080756887729352744634151d0 !! <tt>sqrt(3)</tt>
  real(kind=real64),    parameter, public :: c_pi      =3.14159265358979323846264338328d0 !! <tt>pi</tt> number
  real(kind=real64),    parameter, public :: c_2pi     =6.28318530717958623199592693709d0 !! <tt>2*pi</tt>
  real(kind=real64),    parameter, public :: c_pi_2    =1.57079632679489661923132169164d0 !! <tt>pi/2</tt>
  real(kind=real64),    parameter, public :: c_pi_4    =0.78539816339744830961566084582d0 !! <tt>pi/4</tt>
  real(kind=real64),    parameter, public :: c_sqrtpi  =1.77245385090551602729816748334d0 !! <tt>sqrt(pi)</tt>
  real(kind=real64),    parameter, public :: c_2_sqrtpi=1.12837916709551257389615890312d0 !! <tt>2/sqrt(pi)</tt>
  real(kind=real64),    parameter, public :: c_1_pi    =0.31830988618379067153776752675d0 !! <tt>1/pi</tt>
  real(kind=real64),    parameter, public :: c_2_pi    =0.63661977236758134307553505349d0 !! <tt>2/pi</tt>
  real(kind=real64),    parameter, public :: c_ln10    =2.30258509299404568401799145468d0 !! <tt>ln(10)</tt>
  real(kind=real64),    parameter, public :: c_ln2     =0.69314718055994530941723212146d0 !! <tt>ln(2)</tt>
  real(kind=real64),    parameter, public :: c_lnpi    =1.14472988584940017414342735135d0 !! <tt>ln(pi)</tt>
  real(kind=real64),    parameter, public :: c_1_2pi   =0.15915494309189533576888376337d0 !! <tt>1/(2*pi)</tt>
  real(kind=real64),    parameter, public :: c_1_4pi   =0.07957747154594767280411105048d0 !! <tt>1/(4*pi)</tt>
  real(kind=real64),    parameter, public :: c_gamma   =0.57721566490153286060651209008d0 !! <tt>gamma</tt> number (Euler constant)
  complex(kind=real64), parameter, public :: c_im      =(0.0d0,1.0d0)                     !! Imaginary unit
  !! Kronecker delta
  real(kind=real64),    parameter, public, dimension(3,3) :: c_dkr=reshape((/1.0d0,0.0d0,0.0d0,&
                                                                             0.0d0,1.0d0,0.0d0,&
                                                                             0.0d0,0.0d0,1.0d0/),(/3,3/))
  real(kind=real64),    parameter, public, dimension(2,2) :: c_lc2=reshape((/ 0.0d0,-1.0d0,&
                                                                              1.0d0, 0.0d0/),(/2,2/))
  ! Identifier for functions
  integer, parameter, public :: fbem_f_any=-1 !! Identifier for a generic function.
  integer, parameter, public :: fbem_f_lnr= 0 !! Identifier for <tt>f=ln(r)</tt>.
  integer, parameter, public :: fbem_f_1r1= 1 !! Identifier for <tt>f=1/r</tt>.
  integer, parameter, public :: fbem_f_1r2= 2 !! Identifier for <tt>f=1/r**2</tt>.
  integer, parameter, public :: fbem_f_1r3= 3 !! Identifier for <tt>f=1/r**3</tt>.
  integer, parameter, public :: fbem_f_1r4= 4 !! Identifier for <tt>f=1/r**3</tt>.
  integer, parameter, public :: fbem_f_1r5= 5 !! Identifier for <tt>f=1/r**3</tt>.
  integer, parameter, public :: fbem_f_1r6= 6 !! Identifier for <tt>f=1/r**3</tt>.
  integer, parameter, public :: fbem_f_1r7= 7 !! Identifier for <tt>f=1/r**3</tt>.

!  !! Relative error used for <tt>fbem_modified_bessel_k0_and_k1</tt>
!  real(kind=real64), parameter, public                 :: fbem_modified_bessel_k0_and_k1_re=1.0d-6

  ! Bessel series constants
  real(kind=real64), dimension(25), parameter :: c_s_I0s=(/0.2500000000000000d+00,0.1562500000000000d-01,0.4340277777777778d-03,&
                                                           0.6781684027777777d-05,0.6781684027777778d-07,0.4709502797067901d-09,&
                                                           0.2402807549524440d-11,0.9385966990329842d-14,0.2896903392077112d-16,&
                                                           0.7242258480192779d-19,0.1496334396734045d-21,0.2597802772107717d-24,&
                                                           0.3842903509035085d-27,0.4901662639075364d-30,0.5446291821194849d-33,&
                                                           0.5318644356635594d-36,0.4600903422695150d-39,0.3550079801462307d-42,&
                                                           0.2458504017633177d-45,0.1536565011020735d-48,0.8710686003518909d-52,&
                                                           0.4499321282809353d-55,0.2126333309456217d-58,0.9228877211181494d-62,&
                                                           0.3691550884472597d-65/)
  real(kind=real64), dimension(25), parameter :: c_s_I1s=(/0.1250000000000000d+00,0.5208333333333333d-02,0.1085069444444444d-03,&
                                                           0.1356336805555556d-05,0.1130280671296296d-07,0.6727861138668430d-10,&
                                                           0.3003509436905549d-12,0.1042885221147760d-14,0.2896903392077112d-17,&
                                                           0.6583871345629799d-20,0.1246945330611704d-22,0.1998309824698244d-25,&
                                                           0.2744931077882203d-28,0.3267775092716909d-31,0.3403932388246780d-34,&
                                                           0.3128614327432703d-37,0.2556057457052861d-40,0.1868463053401214d-43,&
                                                           0.1229252008816588d-46,0.7316976242955882d-50,0.3959402728872232d-53,&
                                                           0.1956226644699719d-56,0.8859722122734235d-60,0.3691550884472597d-63,&
                                                           0.1419827263258691d-66/)
  real(kind=real64), dimension(25), parameter :: c_s_I2s=(/0.4166666666666666d-01,0.1302083333333333d-02, 0.2170138888888889d-04,&
                                                           0.2260561342592593d-06, 0.1614686673280423d-08,0.8409826423335537d-11,&
                                                           0.3337232707672833d-13,0.1042885221147760d-15, 0.2633548538251920d-18,&
                                                           0.5486559454691499d-21, 0.9591887158551572d-24,0.1427364160498746d-26,&
                                                           0.1829954051921469d-29,0.2042359432948068d-32, 0.2002313169556930d-35,&
                                                           0.1738119070795946d-38, 0.1345293398448874d-41,0.9342315267006071d-45,&
                                                           0.5853580994364706d-48,0.3325898292252674d-51, 0.1721479447335753d-54,&
                                                           0.8150944352915494d-58, 0.3543888849093694d-61,0.1419827263258691d-64,&
                                                           0.5258619493550709d-68/)
  real(kind=real64), dimension(25), parameter :: c_s_K0s=(/0.1056960837746168d+00,0.1441850523591355d-01, 0.5451899602568578d-03,&
                                                           0.1021401413595785d-04, 0.1157035094151340d-06,0.8819883064451181d-09,&
                                                           0.4843198560366339d-11,0.2009199025022224d-13, 0.6523109713385804d-16,&
                                                           0.1703200013148379d-18, 0.3655038691331977d-21,0.6562036847904769d-24,&
                                                           0.1000276306268431d-26,0.1310874511539864d-29, 0.1492835847185592d-32,&
                                                           0.1491089034246152d-35, 0.1316933544567889d-38,0.1035875091927791d-41,&
                                                           0.7303045169403070d-45,0.4641231481427954d-48, 0.2672563063411289d-51,&
                                                           0.1400907588171493d-54, 0.6712995532223491d-58,0.2952080188129701d-61,&
                                                           0.1195598278789771d-64/)
  real(kind=real64), dimension(25), parameter :: c_s_K1s=(/0.1681960837746168d+00,0.1134844793505348d-01, 0.2997217162395400d-03,&
                                                           0.4356873015494250d-05, 0.4045163759053851d-07,0.2616078891824172d-09,&
                                                           0.1248343508052904d-11,0.4580762857954694d-14, 0.1333590976597932d-16,&
                                                           0.3156580672502777d-19, 0.6195643263104267d-22,0.1024915744483027d-24,&
                                                           0.1448572802368345d-27,0.1769617849337929d-30, 0.1887319386408531d-33,&
                                                           0.1772626006921546d-36, 0.1477459813170169d-39,0.1100228849678733d-42,&
                                                           0.7364507769843892d-46,0.4455063202516885d-49, 0.2447600070050589d-52,&
                                                           0.1226685844691296d-55, 0.5631078452364297d-59,0.2376430354041650d-62,&
                                                           0.9251518577738951d-66/)
  real(kind=real64), dimension(25), parameter :: c_s_K2s=(/0.6995425014709447d-01,0.3162632817096702d-02, 0.6428462102568579d-04,&
                                                           0.7638215249589183d-06, 0.6009474894831275d-08,0.3375221445071909d-10,&
                                                           0.1424128705699592d-12,0.4685051380069469d-15, 0.1236296783618592d-17,&
                                                           0.2676205222541409d-20, 0.4839663180530602d-23,0.7422781329200106d-26,&
                                                           0.9779148952583729d-29,0.1118775902292131d-31, 0.1121966187120059d-34,&
                                                           0.9944484431275031d-38, 0.7846909195550826d-41,0.5547855824728694d-44,&
                                                           0.3534782657041685d-47,0.2040146447927005d-50, 0.1071658636749542d-53,&
                                                           0.5145153287684215d-57, 0.2266606936342094d-60,0.9194725487208603d-64,&
                                                           0.3445964730620170d-67/)
  real(kind=real64), dimension(25), parameter :: c_l_K0s=(/-0.1250000000000000d+00,0.7031250000000000d-01, -0.7324218750000000d-01,&
                                                            0.1121520996093750d+00, -0.2271080017089844d+00,0.5725014209747314d+00,&
                                                           -0.1727727502584457d+01,0.6074042001273483d+01, -0.2438052969955606d+02,&
                                                            0.1100171402692467d+03, -0.5513358961220206d+03,0.3038090510922384d+04,&
                                                           -0.1825775547429317d+05,0.1188384262567832d+06, -0.8328593040162892d+06,&
                                                            0.6252951493434796d+07, -0.5006958953198892d+08,0.4259392165047668d+09,&
                                                           -0.3836255180230433d+10,0.3646840080706556d+11, -0.3649010818849834d+12,&
                                                            0.3833534661393945d+13, -0.4218971570284096d+14,0.4854014686852901d+15,&
                                                           -0.5827244631566907d+16/)
  real(kind=real64), dimension(25), parameter :: c_l_K1s=(/0.3750000000000000d+00,-0.1171875000000000d+00, 0.1025390625000000d+00,&
                                                          -0.1441955566406250d+00, 0.2775764465332031d+00,-0.6765925884246826d+00,&
                                                           0.1993531733751297d+01,-0.6883914268109947d+01, 0.2724882731126855d+02,&
                                                          -0.1215978918765359d+03, 0.6038440767050702d+03,-0.3302272294480852d+04,&
                                                           0.1971837591223662d+05,-0.1276412726461746d+06, 0.8902978767070677d+06,&
                                                          -0.6656367718817687d+07, 0.5310411010968521d+08,-0.4502786003050392d+09,&
                                                           0.4043620325107754d+10,-0.3833857520742789d+11, 0.3827011346598605d+12,&
                                                          -0.4011838599133197d+13, 0.4406481417852277d+14,-0.5060568503314725d+15,&
                                                           0.6065091351222699d+16/)

  interface fbem_invert_2x2_matrix
    module procedure fbem_invert_2x2_matrix_r64
    module procedure fbem_invert_2x2_matrix_c64
    module procedure fbem_invert_2x2_matrix_r128
    module procedure fbem_invert_2x2_matrix_c128
  end interface fbem_invert_2x2_matrix

  interface fbem_invert_3x3_matrix
    module procedure fbem_invert_3x3_matrix_r64
    module procedure fbem_invert_3x3_matrix_c64
    module procedure fbem_invert_3x3_matrix_r128
    module procedure fbem_invert_3x3_matrix_c128
  end interface fbem_invert_3x3_matrix

  interface fbem_invert_4x4_matrix
    module procedure fbem_invert_4x4_matrix_r64
    module procedure fbem_invert_4x4_matrix_c64
    module procedure fbem_invert_4x4_matrix_r128
    module procedure fbem_invert_4x4_matrix_c128
  end interface fbem_invert_4x4_matrix

  interface fbem_cross_product
    module procedure fbem_cross_product_r64
    module procedure fbem_cross_product_r128
  end interface fbem_cross_product

contains

  !! Function that returns the argument (angle) of a complex number
  function fbem_zarg(z)
    implicit none
    real(kind=real64)    :: fbem_zarg
    complex(kind=real64) :: z
    fbem_zarg=atan2(imag(z),real(z))
  end function fbem_zarg

  !! Safe solution of the quadratic equation for complex coefficients. From:
  !!   - Numerical Recipes in Fortran 77: The Art of Scientific Computing
  subroutine fbem_solve_quadratic_equation(a,b,c,x1,x2)
    implicit none
    ! I/O
    complex(kind=real64) :: a, b, c
    complex(kind=real64) :: x1, x2
    ! Local
    complex(kind=real64) :: bp, q, sqrtbp2mac
    ! Compute
    bp=0.5d0*b
    sqrtbp2mac=sqrt(bp**2-a*c)
    if (dble(conjg(bp)*sqrtbp2mac).ge.0.d0) then
      q=-bp-sqrtbp2mac
    else
      q=-bp+sqrtbp2mac
    end if
    x1=q/a
    x2=c/q
  end subroutine fbem_solve_quadratic_equation

  !! Inverse and determinant of a 2x2 matrix
  subroutine fbem_invert_2x2_matrix_r64(A,InvA,DetA,info)
    implicit none
    real(kind=real64) :: A(2,2)
    real(kind=real64) :: InvA(2,2)
    real(kind=real64) :: DetA
    integer, optional :: info
    InvA(1,1)= A(2,2)
    InvA(1,2)=-A(1,2)
    InvA(2,1)=-A(2,1)
    InvA(2,2)= A(1,1)
    DetA=A(1,1)*A(2,2)-A(1,2)*A(2,1)
    if (present(info)) then
      if (abs(DetA).le.1.d-15) then
        InvA=0.d0
        info=-1
        return
      else
        info=0
      end if
    end if
    InvA=InvA/DetA
  end subroutine fbem_invert_2x2_matrix_r64

  !! Inverse and determinant of a 2x2 matrix
  subroutine fbem_invert_2x2_matrix_r128(A,InvA,DetA,info)
    implicit none
    real(kind=real128) :: A(2,2)
    real(kind=real128) :: InvA(2,2)
    real(kind=real128) :: DetA
    integer, optional  :: info
    InvA(1,1)= A(2,2)
    InvA(1,2)=-A(1,2)
    InvA(2,1)=-A(2,1)
    InvA(2,2)= A(1,1)
    DetA=A(1,1)*A(2,2)-A(1,2)*A(2,1)
    if (present(info)) then
      if (abs(DetA).le.1.d-33) then
        InvA=0.d0
        info=-1
        return
      else
        info=0
      end if
    end if
    InvA=InvA/DetA
  end subroutine fbem_invert_2x2_matrix_r128

  !! Inverse and determinant of a 2x2 matrix
  subroutine fbem_invert_2x2_matrix_c64(A,InvA,DetA,info)
    implicit none
    complex(kind=real64) :: A(2,2)
    complex(kind=real64) :: InvA(2,2)
    complex(kind=real64) :: DetA
    integer, optional    :: info
    InvA(1,1)= A(2,2)
    InvA(1,2)=-A(1,2)
    InvA(2,1)=-A(2,1)
    InvA(2,2)= A(1,1)
    DetA=A(1,1)*A(2,2)-A(1,2)*A(2,1)
    if (present(info)) then
      if (abs(DetA).le.1.d-15) then
        InvA=0.d0
        info=-1
        return
      else
        info=0
      end if
    end if
    InvA=InvA/DetA
  end subroutine fbem_invert_2x2_matrix_c64

  !! Inverse and determinant of a 2x2 matrix
  subroutine fbem_invert_2x2_matrix_c128(A,InvA,DetA,info)
    implicit none
    complex(kind=real128) :: A(2,2)
    complex(kind=real128) :: InvA(2,2)
    complex(kind=real128) :: DetA
    integer, optional     :: info
    InvA(1,1)= A(2,2)
    InvA(1,2)=-A(1,2)
    InvA(2,1)=-A(2,1)
    InvA(2,2)= A(1,1)
    DetA=A(1,1)*A(2,2)-A(1,2)*A(2,1)
    if (present(info)) then
      if (abs(DetA).le.1.q-33) then
        InvA=0.d0
        info=-1
        return
      else
        info=0
      end if
    end if
    InvA=InvA/DetA
  end subroutine fbem_invert_2x2_matrix_c128

  !! Inverse and determinant of a 3x3 matrix
  subroutine fbem_invert_3x3_matrix_r64(A,InvA,DetA,info)
    implicit none
    real(kind=real64) :: A(3,3)
    real(kind=real64) :: InvA(3,3)
    real(kind=real64) :: DetA
    integer, optional  :: info
    InvA=0.d0
    InvA(1,1)=A(2,2)*A(3,3)-A(3,2)*A(2,3)
    InvA(1,2)=A(3,2)*A(1,3)-A(1,2)*A(3,3)
    InvA(1,3)=A(1,2)*A(2,3)-A(1,3)*A(2,2)
    DetA=InvA(1,1)*A(1,1)+InvA(1,2)*A(2,1)+InvA(1,3)*A(3,1)
    if (present(info)) then
      if (abs(DetA).le.1.d-15) then
        InvA=0.d0
        info=-1
        return
      else
        info=0
      end if
    end if
    InvA(1,1)=InvA(1,1)/DetA
    InvA(1,2)=InvA(1,2)/DetA
    InvA(1,3)=InvA(1,3)/DetA
    InvA(2,1)=(A(2,3)*A(3,1)-A(2,1)*A(3,3))/DetA
    InvA(2,2)=(A(1,1)*A(3,3)-A(3,1)*A(1,3))/DetA
    InvA(2,3)=(A(2,1)*A(1,3)-A(1,1)*A(2,3))/DetA
    InvA(3,1)=(A(2,1)*A(3,2)-A(2,2)*A(3,1))/DetA
    InvA(3,2)=(A(3,1)*A(1,2)-A(1,1)*A(3,2))/DetA
    InvA(3,3)=(A(1,1)*A(2,2)-A(1,2)*A(2,1))/DetA
  end subroutine fbem_invert_3x3_matrix_r64

  !! Inverse and determinant of a 3x3 matrix
  subroutine fbem_invert_3x3_matrix_r128(A,InvA,DetA,info)
    implicit none
    real(kind=real128) :: A(3,3)
    real(kind=real128) :: InvA(3,3)
    real(kind=real128) :: DetA
    integer, optional  :: info
    InvA=0.d0
    InvA(1,1)=A(2,2)*A(3,3)-A(3,2)*A(2,3)
    InvA(1,2)=A(3,2)*A(1,3)-A(1,2)*A(3,3)
    InvA(1,3)=A(1,2)*A(2,3)-A(1,3)*A(2,2)
    DetA=InvA(1,1)*A(1,1)+InvA(1,2)*A(2,1)+InvA(1,3)*A(3,1)
    if (present(info)) then
      if (abs(DetA).le.1.q-33) then
        InvA=0.d0
        info=-1
        return
      else
        info=0
      end if
    end if
    InvA(1,1)=InvA(1,1)/DetA
    InvA(1,2)=InvA(1,2)/DetA
    InvA(1,3)=InvA(1,3)/DetA
    InvA(2,1)=(A(2,3)*A(3,1)-A(2,1)*A(3,3))/DetA
    InvA(2,2)=(A(1,1)*A(3,3)-A(3,1)*A(1,3))/DetA
    InvA(2,3)=(A(2,1)*A(1,3)-A(1,1)*A(2,3))/DetA
    InvA(3,1)=(A(2,1)*A(3,2)-A(2,2)*A(3,1))/DetA
    InvA(3,2)=(A(3,1)*A(1,2)-A(1,1)*A(3,2))/DetA
    InvA(3,3)=(A(1,1)*A(2,2)-A(1,2)*A(2,1))/DetA
  end subroutine fbem_invert_3x3_matrix_r128

  !! Inverse and determinant of a 3x3 matrix
  subroutine fbem_invert_3x3_matrix_c64(A,InvA,DetA,info)
    implicit none
    complex(kind=real64) :: A(3,3)
    complex(kind=real64) :: InvA(3,3)
    complex(kind=real64) :: DetA
    integer, optional    :: info
    InvA=0.d0
    InvA(1,1)=A(2,2)*A(3,3)-A(3,2)*A(2,3)
    InvA(1,2)=A(3,2)*A(1,3)-A(1,2)*A(3,3)
    InvA(1,3)=A(1,2)*A(2,3)-A(1,3)*A(2,2)
    DetA=InvA(1,1)*A(1,1)+InvA(1,2)*A(2,1)+InvA(1,3)*A(3,1)
    if (present(info)) then
      if (abs(DetA).le.1.d-15) then
        InvA=0.d0
        info=-1
        return
      else
        info=0
      end if
    end if
    InvA(1,1)=InvA(1,1)/DetA
    InvA(1,2)=InvA(1,2)/DetA
    InvA(1,3)=InvA(1,3)/DetA
    InvA(2,1)=(A(2,3)*A(3,1)-A(2,1)*A(3,3))/DetA
    InvA(2,2)=(A(1,1)*A(3,3)-A(3,1)*A(1,3))/DetA
    InvA(2,3)=(A(2,1)*A(1,3)-A(1,1)*A(2,3))/DetA
    InvA(3,1)=(A(2,1)*A(3,2)-A(2,2)*A(3,1))/DetA
    InvA(3,2)=(A(3,1)*A(1,2)-A(1,1)*A(3,2))/DetA
    InvA(3,3)=(A(1,1)*A(2,2)-A(1,2)*A(2,1))/DetA
  end subroutine fbem_invert_3x3_matrix_c64

  !! Inverse and determinant of a 3x3 matrix
  subroutine fbem_invert_3x3_matrix_c128(A,InvA,DetA,info)
    implicit none
    complex(kind=real128) :: A(3,3)
    complex(kind=real128) :: InvA(3,3)
    complex(kind=real128) :: DetA
    integer, optional     :: info
    InvA=0.d0
    InvA(1,1)=A(2,2)*A(3,3)-A(3,2)*A(2,3)
    InvA(1,2)=A(3,2)*A(1,3)-A(1,2)*A(3,3)
    InvA(1,3)=A(1,2)*A(2,3)-A(1,3)*A(2,2)
    DetA=InvA(1,1)*A(1,1)+InvA(1,2)*A(2,1)+InvA(1,3)*A(3,1)
    if (present(info)) then
      if (abs(DetA).le.1.q-33) then
        InvA=0.d0
        info=-1
        return
      else
        info=0
      end if
    end if
    InvA(1,1)=InvA(1,1)/DetA
    InvA(1,2)=InvA(1,2)/DetA
    InvA(1,3)=InvA(1,3)/DetA
    InvA(2,1)=(A(2,3)*A(3,1)-A(2,1)*A(3,3))/DetA
    InvA(2,2)=(A(1,1)*A(3,3)-A(3,1)*A(1,3))/DetA
    InvA(2,3)=(A(2,1)*A(1,3)-A(1,1)*A(2,3))/DetA
    InvA(3,1)=(A(2,1)*A(3,2)-A(2,2)*A(3,1))/DetA
    InvA(3,2)=(A(3,1)*A(1,2)-A(1,1)*A(3,2))/DetA
    InvA(3,3)=(A(1,1)*A(2,2)-A(1,2)*A(2,1))/DetA
  end subroutine fbem_invert_3x3_matrix_c128

  subroutine fbem_invert_4x4_matrix_r64(A,InvA,DetA,info)
    implicit none
    real(kind=real64) :: A(4,4)
    real(kind=real64) :: InvA(4,4)
    real(kind=real64) :: DetA
    integer, optional :: info
    real(kind=real64) :: CofA(4,4)
    DetA=A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-&
         A(3,3)*A(4,2)))-A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+&
         A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))+A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-&
         A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))-A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+&
         A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    if (present(info)) then
      if (abs(DetA).le.1.q-15) then
        InvA=0.d0
        info=-1
        return
      else
        info=0
      end if
    end if
    CofA(1,1)=A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))
    CofA(1,2)=A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))
    CofA(1,3)=A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
    CofA(1,4)=A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
    CofA(2,1)=A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))
    CofA(2,2)=A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))
    CofA(2,3)=A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
    CofA(2,4)=A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
    CofA(3,1)=A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2))
    CofA(3,2)=A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3))
    CofA(3,3)=A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1))
    CofA(3,4)=A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2))
    CofA(4,1)=A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3))
    CofA(4,2)=A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
    CofA(4,3)=A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2))
    CofA(4,4)=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
    InvA=transpose(CofA)/DetA
  end subroutine fbem_invert_4x4_matrix_r64

  subroutine fbem_invert_4x4_matrix_r128(A,InvA,DetA,info)
    implicit none
    real(kind=real128) :: A(4,4)
    real(kind=real128) :: InvA(4,4)
    real(kind=real128) :: DetA
    integer, optional  :: info
    real(kind=real128) :: CofA(4,4)
    DetA=A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-&
         A(3,3)*A(4,2)))-A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+&
         A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))+A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-&
         A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))-A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+&
         A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    if (present(info)) then
      if (abs(DetA).le.1.q-33) then
        InvA=0.d0
        info=-1
        return
      else
        info=0
      end if
    end if
    CofA(1,1)=A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))
    CofA(1,2)=A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))
    CofA(1,3)=A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
    CofA(1,4)=A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
    CofA(2,1)=A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))
    CofA(2,2)=A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))
    CofA(2,3)=A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
    CofA(2,4)=A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
    CofA(3,1)=A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2))
    CofA(3,2)=A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3))
    CofA(3,3)=A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1))
    CofA(3,4)=A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2))
    CofA(4,1)=A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3))
    CofA(4,2)=A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
    CofA(4,3)=A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2))
    CofA(4,4)=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
    InvA=transpose(CofA)/DetA
  end subroutine fbem_invert_4x4_matrix_r128

  subroutine fbem_invert_4x4_matrix_c64(A,InvA,DetA,info)
    implicit none
    complex(kind=real64) :: A(4,4)
    complex(kind=real64) :: InvA(4,4)
    complex(kind=real64) :: DetA
    integer, optional    :: info
    complex(kind=real64) :: CofA(4,4)
    DetA=A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-&
         A(3,3)*A(4,2)))-A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+&
         A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))+A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-&
         A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))-A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+&
         A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    if (present(info)) then
      if (abs(DetA).le.1.q-15) then
        InvA=0.d0
        info=-1
        return
      else
        info=0
      end if
    end if
    CofA(1,1)=A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))
    CofA(1,2)=A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))
    CofA(1,3)=A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
    CofA(1,4)=A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
    CofA(2,1)=A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))
    CofA(2,2)=A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))
    CofA(2,3)=A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
    CofA(2,4)=A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
    CofA(3,1)=A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2))
    CofA(3,2)=A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3))
    CofA(3,3)=A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1))
    CofA(3,4)=A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2))
    CofA(4,1)=A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3))
    CofA(4,2)=A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
    CofA(4,3)=A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2))
    CofA(4,4)=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
    InvA=transpose(CofA)/DetA
  end subroutine fbem_invert_4x4_matrix_c64

  subroutine fbem_invert_4x4_matrix_c128(A,InvA,DetA,info)
    implicit none
    complex(kind=real128) :: A(4,4)
    complex(kind=real128) :: InvA(4,4)
    complex(kind=real128) :: DetA
    integer, optional     :: info
    complex(kind=real128) :: CofA(4,4)
    DetA=A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-&
         A(3,3)*A(4,2)))-A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+&
         A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))+A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-&
         A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))-A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+&
         A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    if (present(info)) then
      if (abs(DetA).le.1.q-33) then
        InvA=0.d0
        info=-1
        return
      else
        info=0
      end if
    end if
    CofA(1,1)=A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))
    CofA(1,2)=A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))
    CofA(1,3)=A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
    CofA(1,4)=A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
    CofA(2,1)=A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))
    CofA(2,2)=A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))
    CofA(2,3)=A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
    CofA(2,4)=A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
    CofA(3,1)=A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2))
    CofA(3,2)=A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3))
    CofA(3,3)=A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1))
    CofA(3,4)=A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2))
    CofA(4,1)=A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3))
    CofA(4,2)=A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
    CofA(4,3)=A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2))
    CofA(4,4)=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
    InvA=transpose(CofA)/DetA
  end subroutine fbem_invert_4x4_matrix_c128

  function fbem_cross_product_r64(a,b)
    implicit none
    real(kind=real64) :: fbem_cross_product_r64(3)
    real(kind=real64) :: a(3)
    real(kind=real64) :: b(3)
    fbem_cross_product_r64(1)=a(2)*b(3)-a(3)*b(2)
    fbem_cross_product_r64(2)=a(3)*b(1)-a(1)*b(3)
    fbem_cross_product_r64(3)=a(1)*b(2)-a(2)*b(1)
  end function

  function fbem_cross_product_r128(a,b)
    implicit none
    real(kind=real128) :: fbem_cross_product_r128(3)
    real(kind=real128) :: a(3)
    real(kind=real128) :: b(3)
    fbem_cross_product_r128(1)=a(2)*b(3)-a(3)*b(2)
    fbem_cross_product_r128(2)=a(3)*b(1)-a(1)*b(3)
    fbem_cross_product_r128(3)=a(1)*b(2)-a(2)*b(1)
  end function

  function fbem_vector_normalization(n,a)
    implicit none
    integer           :: n
    real(kind=real64) :: fbem_vector_normalization(n)
    real(kind=real64) :: a(n)
    fbem_vector_normalization=a/sqrt(dot_product(a,a))
  end function fbem_vector_normalization

  function fbem_vector_norm(n,a)
    implicit none
    integer           :: n
    real(kind=real64) :: fbem_vector_norm
    real(kind=real64) :: a(n)
    fbem_vector_norm=sqrt(dot_product(a,a))
  end function fbem_vector_norm

  ! -------------------- !
  ! POLYNOMIAL FUNCTIONS !
  ! ---------------------!

  function f_poly_1d_linear(c,x)
    implicit none
    real(kind=real64) :: f_poly_1d_linear
    real(kind=real64) :: c(0:1)
    real(kind=real64) :: x
    f_poly_1d_linear=c(0)+c(1)*x
  end function f_poly_1d_linear

  function f_poly_1d_quadratic(c,x)
    implicit none
    real(kind=real64) :: f_poly_1d_quadratic
    real(kind=real64) :: c(0:2)
    real(kind=real64) :: x
    f_poly_1d_quadratic=c(0)+c(1)*x+c(2)*x**2
  end function f_poly_1d_quadratic

  function f_poly_2d_linear(c,x,y)
    implicit none
    real(kind=real64) :: f_poly_2d_linear
    real(kind=real64) :: c(0:2)
    real(kind=real64) :: x, y
    f_poly_2d_linear=c(0)+c(1)*x+c(2)*y
  end function f_poly_2d_linear

  function f_poly_2d_quadratic(c,x,y)
    implicit none
    real(kind=real64) :: f_poly_2d_quadratic
    real(kind=real64) :: c(0:5)
    real(kind=real64) :: x, y
    f_poly_2d_quadratic=c(0)+c(1)*x+c(2)*y+c(3)*x*y+c(4)*x**2+c(5)*y**2
  end function f_poly_2d_quadratic

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! Bessel series study for argument z (pure imaginary) from 1.0d-12 to 1.0d12.
  ! ===========================================================================
  !
  ! 64 bits floats are used.
  !
  ! The number of terms are truncated when the relative error is 1.0d-6.
  ! --------------------------------------------------------------------
  !
  ! It is found that:
  !  - The series for small arguments breaks when |z|>30 due finite precision arithmetics
  !  - The series for large arguments breaks when |z|<6.5 due finite precision arithmetics
  !  - The relative error between both series is approximately 1.0d-6 in the range 6.5<|z|<20
  !  - The number of terms of both series are approximately equal at |z|=7. This value seems to be the limit for both series
  !    in order to minimize the computational cost.
  !
  ! Conclusion:
  !
  !  - When |z|< 7, the series for small arguments are used.
  !    - If |z|<=1e-3, then N=2.
  !    - If |z|> 1e-3, then N=nint(1./(-0.107*log10(|z|)+0.163))
  !  - When |z|>=7, the series for large arguments are used. The number of terms is N=2 for |z|>1e3.
  !    - If |z|< 1e3, then N=nint(1./(0.21*log10(|z|)-0.065))
  !    - If |z|>=1e3, then N=2
  !
  ! The number of terms are truncated when the relative error is 1.0d-9.
  ! --------------------------------------------------------------------
  !
  ! It is found that:
  !  - The series for small arguments breaks when |z|>30 due finite precision arithmetics
  !  - The series for large arguments breaks when |z|<9.75 due finite precision arithmetics
  !  - The relative error between both series is approximately 1.0d-9 in the range 10<|z|<20
  !  - The number of terms of both series are approximately equal at |z|=10. This value seems to be the limit for both series
  !    in order to minimize the computational cost.
  !
  ! Conclusion:
  !
  !  - When |z|< 10, the series for small arguments are used.
  !    - If |z|<=1e-4, then N=2.
  !    - If |z|> 1e-4, then N=nint(1./(-0.08*log10(|z|)+0.1277))
  !  - When |z|>=10, the series for large arguments are used. The number of terms is N=2 for |z|>1e3.
  !    - If |z|< 1e4, then N=nint(1./(0.135*log10(|z|)-0.0592))
  !    - If |z|>=1e4, then N=2
  !
  !
  ! The subroutines above are programmed with relative error 1.0d-6
  ! --------------------------------------------------------------------------------------------------------------------------------

  !! Calculation of K0 and K1 (modified Bessel functions of the second kind of order 0 and 1) with complex argument.
  !! From:
  !!   -Abramowitz & Stegun, Handbook of Mathematical Functions with Formulas, Graphs and Tables, 1972.
  subroutine fbem_modified_bessel_K0_and_K1(z, K0, K1)
    implicit none
    complex(kind=real64) :: z              !! Argument
    complex(kind=real64) :: K0             !! K0(z)
    complex(kind=real64) :: K1             !! K1(z)
    integer              :: k              ! Terms counter
    integer              :: k_max          ! Maximum number of terms of the series
    ! Calculations variables
    real(kind=real64)    :: absz
    complex(kind=real64) :: z2, z_2, powz2k, logz_2, powzk, pow1zk, sqrtpi2zemz
    ! Series values
    complex(kind=real64) :: I0s, I1s, K0s, K1s
    ! abs(z)
    absz=zabs(z)
    ! z^2
    z2=z**2
    ! For "small" arguments |z|<7
    if (absz.lt.7.0d0) then
      ! Number of terms in the series
      if (absz.le.1.0d-3) then
        k_max=2
      else
        k_max=nint(1.0d0/(-0.107d0*dlog10(absz)+0.163d0))
      end if
      ! Calculation of series IOs, I1s, K0s, K1s
      ! Initialization
      I0s=(0.0d0,0.0d0)
      I1s=(0.0d0,0.0d0)
      K0s=(0.0d0,0.0d0)
      K1s=(0.0d0,0.0d0)
      powz2k=1.0d0
      ! Iterative procedure from k=1 to k_max or when relative error is reached
      do k=1,k_max
        powz2k=powz2k*z2
        I0s=I0s+c_s_I0s(k)*powz2k
        I1s=I1s+c_s_I1s(k)*powz2k
        K0s=K0s+c_s_K0s(k)*powz2k
        K1s=K1s+c_s_K1s(k)*powz2k
      end do
      ! Calculate definitive mlndz2, K0R, d1z and K1R
      z_2=0.5d0*z
      logz_2=zlog(z_2)
      K0=-logz_2-c_gamma-logz_2*I0s+K0s
      K1=1.0d0/z+z_2*(logz_2+c_gamma-0.5d0)+z_2*(logz_2*I1s-0.5d0*K1s)
    ! For "large" arguments |z|>=7
    else
      ! Number of terms in the series
      if (absz.ge.1.0d3) then
        k_max=2
      else
        k_max=nint(1.0d0/(0.21d0*dlog10(absz)-0.065d0))
      end if
      ! Calculation of series K0s, K1s
      ! Initialization
      K0s=1.0d0
      K1s=1.0d0
      powzk=1
      ! Iterative procedure from k=1 to k_max or when relative error is reached
      do k=1,k_max
        powzk=powzk*z
        pow1zk=1.0d0/powzk
        K0s=K0s+c_l_K0s(k)*pow1zk
        K1s=K1s+c_l_K1s(k)*pow1zk
      end do
      sqrtpi2zemz=zsqrt(0.5d0*c_pi/z)*zexp(-z)
      ! Total K0, K1 and K2
      K0=sqrtpi2zemz*K0s
      K1=sqrtpi2zemz*K1s
    end if
  end subroutine fbem_modified_bessel_K0_and_K1

  !! Calculation of K0R, K1R and K2R (parts of modified Bessel functions of the second kind of order 0, 1 and 2) with complex
  !! argument. The routine returns K0R, K1R and K2R, and in order to get K0, K1 and K2:
  !!   - K0(z) = -log(z/2) - gamma + K0R(z)
  !!   - K1(z) = 1/z + z/2*(log(z/2) + gamma - 1/2) + K1R(z)
  !!   - K2(z) = 2/z^2 - 1/2 + K2R(z)
  !!
  !! K0R, K1R and K2R -> 0 as |z|->0, are at least of order O(z*log(z)).
  !!
  !! From:
  !!   -Abramowitz & Stegun, Handbook of Mathematical Functions with Formulas, Graphs and Tables, 1972.
  subroutine fbem_modified_bessel_K0R_K1R_K2R(z,K0R,K1R,K2R)
    implicit none
    complex(kind=real64) :: z              !! Argument
    complex(kind=real64) :: K0R            !! K0R(z)
    complex(kind=real64) :: K1R            !! K1R(z)
    complex(kind=real64) :: K2R            !! K2R(z)
    integer              :: k              ! Terms counter
    integer              :: k_max          ! Maximum number of terms of the series
    ! Calculations variables
    real(kind=real64)    :: absz
    complex(kind=real64) :: z2, powz2k, logz_2, z_2, powzk, pow1zk, sqrtpi2zemz, div1z
    ! Series values
    complex(kind=real64) :: I0s, I1s, K0s, K1s
    ! abs(z)
    absz=zabs(z)
    ! z^2
    z2=z**2
    ! For "small" arguments |z|<7
    if (absz.lt.7.0d0) then
      ! Number of terms in the series
      if (absz.le.1.0d-3) then
        k_max=2
      else
        k_max=nint(1.0d0/(-0.107d0*dlog10(absz)+0.163d0))
      end if
      ! Calculation of series IOs, I1s, K0s, K1s
      ! Initialization
      I0s=(0.0d0,0.0d0)
      I1s=(0.0d0,0.0d0)
      K0s=(0.0d0,0.0d0)
      K1s=(0.0d0,0.0d0)
      powz2k=1.0d0
      ! Iterative procedure from k=1 to k_max or when relative error is reached
      do k=1,k_max
        powz2k=powz2k*z2
        I0s=I0s+c_s_I0s(k)*powz2k
        I1s=I1s+c_s_I1s(k)*powz2k
        K0s=K0s+c_s_K0s(k)*powz2k
        K1s=K1s+c_s_K1s(k)*powz2k
      end do
      ! Calculate definitive mlndz2, K0R, d1z and K1R
      z_2=0.5d0*z
      logz_2=zlog(z_2)
      K0R=-logz_2*I0s+K0s
      K1R=z_2*(logz_2*I1s-0.5d0*K1s)
      K2R=K0R+2.0d0/z*K1R
    ! For "large" arguments |z|>=7
    else
      ! Number of terms in the series
      if (absz.ge.1.0d3) then
        k_max=2
      else
        k_max=nint(1.0d0/(0.21d0*dlog10(absz)-0.065d0))
      end if
      ! Calculation of series K0s, K1s
      ! Initialization
      K0s=1.0d0
      K1s=1.0d0
      powzk=1
      ! Iterative procedure from k=1 to k_max or when relative error is reached
      do k=1,k_max
        powzk=powzk*z
        pow1zk=1.0d0/powzk
        K0s=K0s+c_l_K0s(k)*pow1zk
        K1s=K1s+c_l_K1s(k)*pow1zk
      end do
      sqrtpi2zemz=zsqrt(0.5d0*c_pi/z)*zexp(-z)
      div1z=1.0d0/z
      z_2=0.5d0*z
      ! Total K0, K1 and K2
      K0R=sqrtpi2zemz*K0s
      K1R=sqrtpi2zemz*K1s
      K2R=K0R+2.0d0*div1z*K1R
      ! Substract the singularities
      logz_2=zlog(0.5d0*z)
      K0R=K0R+logz_2+c_gamma
      K1R=K1R-div1z-z_2*(logz_2+c_gamma-0.5d0)
      K2R=K2R-2.0d0/(z2)+0.5d0
    end if
  end subroutine fbem_modified_bessel_K0R_K1R_K2R

  !! Calculation of K0R, K1R and K2R (parts of modified Bessel functions of the second kind of order 0, 1 and 2) with complex
  !! argument. The routine returns K0R, K1R and K2R, and in order to get K0, K1 and K2:
  !!   - K0(z) = -log(z/2) - gamma + K0R(z)
  !!   - K1(z) = 1/z + z/2*(log(z/2) + gamma - 1/2) + K1R(z)
  !!   - K2(z) = 2/z^2 - 1/2 -z^2/8*(log(z/2) + gamma - 3/4) + K2R(z)
  !!
  !! K0R is of order O(z)
  !! K1R is of order O(z^2)
  !! K2R is of order O(z^3)
  !!
  !! This routine uses the precalculated coefficients of the series.
  !!
  !! From:
  !!   -Abramowitz & Stegun, Handbook of Mathematical Functions with Formulas, Graphs and Tables, 1972.
  subroutine fbem_modified_bessel_K0R_K1R_K2R_2(z,K0R,K1R,K2R)
    implicit none
    complex(kind=real64) :: z              !! Argument
    complex(kind=real64) :: K0R            !! K0R(z)
    complex(kind=real64) :: K1R            !! K1R(z)
    complex(kind=real64) :: K2R            !! K2R(z)
    integer              :: k              ! Terms counter
    integer              :: k_max          ! Maximum number of terms of the series
    ! Calculations variables
    real(kind=real64)    :: absz
    complex(kind=real64) :: z2, powz2k, logz_2, z_2, powzk, pow1zk, sqrtpi2zemz, div1z
    ! Series values
    complex(kind=real64) :: I0s, I1s, I2s, K0s, K1s, K2s
    ! abs(z)
    absz=abs(z)
    ! z^2
    z2=z**2
    ! For "small" arguments |z|<7
    if (absz.lt.7.0d0) then
      ! Number of terms in the series
      if (absz.le.1.0d-3) then
        k_max=2
      else
        k_max=nint(1.0d0/(-0.107d0*log10(absz)+0.163d0))
      end if
      ! Calculation of series IOs, I1s, K0s, K1s
      ! Initialization
      I0s=(0.0d0,0.0d0)
      I1s=(0.0d0,0.0d0)
      I2s=(0.0d0,0.0d0)
      K0s=(0.0d0,0.0d0)
      K1s=(0.0d0,0.0d0)
      K2s=(0.0d0,0.0d0)
      powz2k=1.0d0
      ! Iterative procedure from k=1 to k_max or when relative error is reached
      do k=1,k_max
        powz2k=powz2k*z2
        I0s=I0s+c_s_I0s(k)*powz2k
        I1s=I1s+c_s_I1s(k)*powz2k
        I2s=I2s+c_s_I2s(k)*powz2k
        K0s=K0s+c_s_K0s(k)*powz2k
        K1s=K1s+c_s_K1s(k)*powz2k
        K2s=K2s+c_s_K2s(k)*powz2k
      end do
      ! Calculate definitive mlndz2, K0R, d1z and K1R
      z_2=0.5d0*z
      logz_2=log(z_2)
      K0R=-logz_2*I0s+K0s
      K1R=z_2*(logz_2*I1s-0.5d0*K1s)
      K2R=-0.25d0*z2*(logz_2*I2s-0.5d0*K2s)
    ! For "large" arguments |z|>=7
    else
      ! Number of terms in the series
      if (absz.ge.1.0d3) then
        k_max=2
      else
        k_max=nint(1.0d0/(0.21d0*log10(absz)-0.065d0))
      end if
      ! Calculation of series K0s, K1s
      ! Initialization
      K0s=1.0d0
      K1s=1.0d0
      powzk=1
      ! Iterative procedure from k=1 to k_max or when relative error is reached
      do k=1,k_max
        powzk=powzk*z
        pow1zk=1.0d0/powzk
        K0s=K0s+c_l_K0s(k)*pow1zk
        K1s=K1s+c_l_K1s(k)*pow1zk
      end do
      sqrtpi2zemz=sqrt(0.5d0*c_pi/z)*exp(-z)
      div1z=1.0d0/z
      z_2=0.5d0*z
      ! Total K0, K1 and K2
      K0R=sqrtpi2zemz*K0s
      K1R=sqrtpi2zemz*K1s
      K2R=K0R+2.0d0*div1z*K1R
      ! Substract the singularities
      logz_2=log(0.5d0*z)
      K0R=K0R+logz_2+c_gamma
      K1R=K1R-div1z-z_2*(logz_2+c_gamma-0.5d0)
      K2R=K2R-2.0d0/(z2)+0.5d0+0.125d0*z2*(logz_2+c_gamma-0.75d0)
    end if
  end subroutine fbem_modified_bessel_K0R_K1R_K2R_2

  !! Calculation of K0R, K1R and K2R (parts of modified Bessel functions of the second kind of order 0, 1 and 2) with complex
  !! argument. The routine returns K0R, K1R and K2R, and in order to get K0, K1 and K2:
  !!   - K0(z) = -log(z/2) - gamma + K0R(z)
  !!   - K1(z) = 1/z + z/2*(log(z/2) + gamma - 1/2) + K1R(z)
  !!   - K2(z) = 2/z^2 - 1/2 -z^2/8*(log(z/2) + gamma - 3/4) + K2R(z)
  !!
  !! K0R is of order O(z)
  !! K1R is of order O(z^2)
  !! K2R is of order O(z^3)
  !!
  !! This routine uses the precalculated coefficients of the series.
  !! This routine calculates them for several arguments.
  !!
  !! From:
  !!   -Abramowitz & Stegun, Handbook of Mathematical Functions with Formulas, Graphs and Tables, 1972.
  !  - When |z|< 10, the series for small arguments are used.
  !    - If |z|<=1e-4, then N=2.
  !    - If |z|> 1e-4, then N=nint(1./(-0.08*log10(|z|)+0.1277))
  !  - When |z|>=10, the series for large arguments are used. The number of terms is N=2 for |z|>1e3.
  !    - If |z|< 1e4, then N=nint(1./(0.135*log10(|z|)-0.0592))
  !    - If |z|>=1e4, then N=2

  subroutine fbem_BesselKnR_decomposed(nz,z,KnR)
    implicit none
    ! I/O
    integer              :: nz             !! Number of arguments
    complex(kind=real64) :: z(nz)          !! Arguments
    complex(kind=real64) :: KnR(0:2,nz)    !! KnR(0,kz) == K0R(z(kz)), KnR(1,kz) == K1R(z(kz)), KnR(2,kz) == K2R(z(kz))
    ! Local
    integer              :: kz, k
    integer              :: k_max
    real(kind=real64)    :: absz
    complex(kind=real64) :: zkz, z2, powz2k, logz_2, z_2, powzk, pow1zk, sqrtpi2zemz, div1z
    complex(kind=real64) :: I0s, I1s, I2s, K0s, K1s, K2s
    ! Loop through the number of arguments
    do kz=1,nz
      zkz=z(kz)
      ! abs(z)
      absz=abs(zkz)
      ! z^2
      z2=zkz**2
!      ! For "small" arguments |z|<7
!      if (absz.lt.7.0d0) then
!        ! Number of terms in the series
!        if (absz.le.1.0d-3) then
!          k_max=2
!        else
!          k_max=nint(1.0d0/(-0.107d0*log10(absz)+0.163d0))
!        end if
      ! For "small" arguments |z|<10
      if (absz.lt.10.0d0) then
        ! Number of terms in the series
        if (absz.le.1.0d-4) then
          k_max=2
        else
          k_max=nint(1.0d0/(-0.08d0*log10(absz)+0.1277d0))
        end if
        ! Calculation of series IOs, I1s, K0s, K1s
        ! Initialization
        I0s=(0.0d0,0.0d0)
        I1s=(0.0d0,0.0d0)
        I2s=(0.0d0,0.0d0)
        K0s=(0.0d0,0.0d0)
        K1s=(0.0d0,0.0d0)
        K2s=(0.0d0,0.0d0)
        powz2k=1.0d0
        ! Iterative procedure from k=1 to k_max or when relative error is reached
        do k=1,k_max
          powz2k=powz2k*z2
          I0s=I0s+c_s_I0s(k)*powz2k
          I1s=I1s+c_s_I1s(k)*powz2k
          I2s=I2s+c_s_I2s(k)*powz2k
          K0s=K0s+c_s_K0s(k)*powz2k
          K1s=K1s+c_s_K1s(k)*powz2k
          K2s=K2s+c_s_K2s(k)*powz2k
        end do
        ! Calculate definitive mlndz2, K0R, d1z and K1R
        z_2=0.5d0*zkz
        logz_2=log(z_2)
        KnR(0,kz)=-logz_2*I0s+K0s
        KnR(1,kz)=z_2*(logz_2*I1s-0.5d0*K1s)
        KnR(2,kz)=-0.25d0*z2*(logz_2*I2s-0.5d0*K2s)
!      ! For "large" arguments |z|>=7
!      else
!        ! Number of terms in the series
!        if (absz.ge.1.0d3) then
!          k_max=2
!        else
!          k_max=nint(1.0d0/(0.21d0*log10(absz)-0.065d0))
!        end if
      ! For "large" arguments |z|>=10
      else
        ! Number of terms in the series
        if (absz.ge.1.0d4) then
          k_max=2
        else
          k_max=nint(1.0d0/(0.135d0*log10(absz)-0.0592d0))
        end if
        ! Calculation of series K0s, K1s
        ! Initialization
        K0s=1.0d0
        K1s=1.0d0
        powzk=1.0d0
        ! Iterative procedure from k=1 to k_max or when relative error is reached
        do k=1,k_max
          powzk=powzk*zkz
          pow1zk=1.0d0/powzk
          K0s=K0s+c_l_K0s(k)*pow1zk
          K1s=K1s+c_l_K1s(k)*pow1zk
        end do
        sqrtpi2zemz=sqrt(0.5d0*c_pi/zkz)*exp(-zkz)
        div1z=1.0d0/zkz
        z_2=0.5d0*zkz
        ! Total K0, K1 and K2
        KnR(0,kz)=sqrtpi2zemz*K0s
        KnR(1,kz)=sqrtpi2zemz*K1s
        KnR(2,kz)=KnR(0,kz)+2.0d0*div1z*KnR(1,kz)
        ! Substract the singularities
        logz_2=log(0.5d0*zkz)
        KnR(0,kz)=KnR(0,kz)+logz_2+c_gamma
        KnR(1,kz)=KnR(1,kz)-div1z-z_2*(logz_2+c_gamma-0.5d0)
        KnR(2,kz)=KnR(2,kz)-2.0d0/(z2)+0.5d0+0.125d0*z2*(logz_2+c_gamma-0.75d0)
      end if
    end do
  end subroutine fbem_BesselKnR_decomposed

  !! Calculation of terms of the expansion of zexp(z)
  subroutine fbem_decomposed_zexp(z,E0,E1,E2,E3,E4)
    implicit none
    ! I/O
    complex*16 :: z  !! Argument z
    complex*16 :: E0 !! zexp(z)
    complex*16 :: E1 !! zexp(z)-1
    complex*16 :: E2 !! zexp(z)-1-z
    complex*16 :: E3 !! zexp(z)-1-z-z**2/2
    complex*16 :: E4 !! zexp(z)-1-z-z**2/2-z**3/6
    ! Local
    real*8     :: absz
    complex*16 :: term(0:19)
    real*8     :: factorial
    complex*16 :: zpowk
    integer    :: k
    integer    :: n
    complex*16 :: expz, z2, z3
    absz=abs(z)
    ! If |z|<=1
    if (absz.le.1.d0) then
      if (absz.le.1.d-6) then
        n=1
      else
        n=nint(10**(1.176d0+0.175d0*log10(absz)))
      end if
      term(0)=(1.d0,0.d0)
      term(1)=z
      factorial=1.d0
      zpowk=z
      do k=2,n+4
        factorial=factorial*dble(k)
        zpowk=zpowk*z
        term(k)=zpowk/factorial
      end do
      E4=(0.d0,0.d0)
      do k=n+4,4,-1
        E4=E4+term(k)
      end do
      E3=E4+term(3)
      E2=E3+term(2)
      E1=E2+term(1)
      E0=E1+term(0)
    ! If |z|>1
    else
      expz=exp(z)
      z2=z*z
      z3=z2*z
      E0=expz
      E1=expz-1.d0
      E2=E1-z
      E3=E2-0.5d0*z2
      E4=E3-0.166666666666666667d0*z3
    end if
  end subroutine fbem_decomposed_zexp

  !! Calculation of terms of the expansion of zexp(z) up to 6
  subroutine fbem_decomposed_zexp_6(z,E0,E1,E2,E3,E4,E5,E6)
    implicit none
    ! I/O
    complex*16 :: z  !! Argument z
    complex*16 :: E0 !! zexp(z)
    complex*16 :: E1 !! zexp(z)-1
    complex*16 :: E2 !! zexp(z)-1-z
    complex*16 :: E3 !! zexp(z)-1-z-z**2/2
    complex*16 :: E4 !! zexp(z)-1-z-z**2/2-z**3/6
    complex*16 :: E5 !! zexp(z)-1-z-z**2/2-z**3/6-z**4/24
    complex*16 :: E6 !! zexp(z)-1-z-z**2/2-z**3/6-z**4/24-z**5/120
    ! Local
    real*8     :: absz
    complex*16 :: term(0:21)
    real*8     :: factorial
    complex*16 :: zpowk
    integer    :: k
    integer    :: n
    complex*16 :: expz, z2, z3, z4, z5
    absz=abs(z)
    ! If |z|<=1
    if (absz.le.1.d0) then
      if (absz.le.1.d-6) then
        n=1
      else
        n=nint(10**(1.176d0+0.175d0*log10(absz)))
      end if
      term(0)=(1.d0,0.d0)
      term(1)=z
      factorial=1.d0
      zpowk=z
      do k=2,n+6
        factorial=factorial*dble(k)
        zpowk=zpowk*z
        term(k)=zpowk/factorial
      end do
      E6=(0.d0,0.d0)
      do k=n+6,6,-1
        E6=E6+term(k)
      end do
      E5=E6+term(5)
      E4=E5+term(4)
      E3=E4+term(3)
      E2=E3+term(2)
      E1=E2+term(1)
      E0=E1+term(0)
    ! If |z|>1
    else
      expz=exp(z)
      z2=z*z
      z3=z2*z
      z4=z3*z
      z5=z4*z
      E0=expz
      E1=expz-1.d0
      E2=E1-z
      E3=E2-0.5d0*z2
      E4=E3-(1.d0/6.d0)*z3
      E5=E4-(1.d0/24.d0)*z4
      E6=E5-(1.d0/120.d0)*z5
    end if
  end subroutine fbem_decomposed_zexp_6

  !! Calculation of terms of the expansion of zexp(z) up to the sixth term for nz complex numbers.
  subroutine fbem_zexp_decomposed(nz,z,E)
    implicit none
    ! I/O
    !! Number of arguments
    integer    :: nz
    !! Arguments z
    complex*16 :: z(nz)
    !! Terms of the expansions:
    !! E(0,kz)=zexp(z(kz))
    !! E(1,kz)=zexp(z(kz))-1
    !! E(2,kz)=zexp(z(kz))-1-z(kz)
    !! E(3,kz)=zexp(z(kz))-1-z(kz)-z(kz)**2/2
    !! E(4,kz)=zexp(z(kz))-1-z(kz)-z(kz)**2/2-z(kz)**3/6
    !! E(5,kz)=zexp(z(kz))-1-z(kz)-z(kz)**2/2-z(kz)**3/6-z(kz)**4/24
    !! E(6,kz)=zexp(z(kz))-1-z(kz)-z(kz)**2/2-z(kz)**3/6-z(kz)**4/24-z(kz)**5/120
    complex*16 :: E(0:6,nz)
    ! Local
    integer    :: kz
    real*8     :: absz
    complex*16 :: zkz
    complex*16 :: term(0:21)
    real*8     :: factorial
    complex*16 :: zpowk
    integer    :: k
    integer    :: n
    complex*16 :: expz, z2, z3, z4, z5
    ! For each z
    do kz=1,nz
      zkz=z(kz)
      absz=abs(zkz)
      ! If |z|<=1
      if (absz.le.1.d0) then
        if (absz.le.1.d-6) then
          n=1
        else
          n=nint(10**(1.176d0+0.175d0*log10(absz)))
        end if
        term(0)=(1.d0,0.d0)
        term(1)=zkz
        factorial=1.d0
        zpowk=zkz
        do k=2,n+6
          factorial=factorial*dble(k)
          zpowk=zpowk*zkz
          term(k)=zpowk/factorial
        end do
        E(6,kz)=(0.d0,0.d0)
        do k=n+6,6,-1
          E(6,kz)=E(6,kz)+term(k)
        end do
        E(5,kz)=E(6,kz)+term(5)
        E(4,kz)=E(5,kz)+term(4)
        E(3,kz)=E(4,kz)+term(3)
        E(2,kz)=E(3,kz)+term(2)
        E(1,kz)=E(2,kz)+term(1)
        E(0,kz)=E(1,kz)+term(0)
      ! If |z|>1
      else
        expz=exp(zkz)
        z2=zkz*zkz
        z3=z2*zkz
        z4=z3*zkz
        z5=z4*zkz
        E(0,kz)=expz
        E(1,kz)=expz-1.d0
        E(2,kz)=E(1,kz)-zkz
        E(3,kz)=E(2,kz)-0.5d0*z2
        E(4,kz)=E(3,kz)-(1.d0/6.d0)*z3
        E(5,kz)=E(4,kz)-(1.d0/24.d0)*z4
        E(6,kz)=E(5,kz)-(1.d0/120.d0)*z5
      end if
    end do
  end subroutine fbem_zexp_decomposed

  !! Standard quicksort routine for sorting (ascend) a vector of integers.
  recursive subroutine fbem_quicksort(start,finish,length,vector)
    implicit none
    integer :: start
    integer :: finish
    integer :: length
    integer :: vector(length)
    integer :: i
    integer :: j
    integer :: x
    integer :: tmp
    i=start
    j=finish
    x=vector((i+j)/2)
    do
      do while (vector(i).lt.x)
        i=i+1
      end do
      do while (vector(j).gt.x)
        j=j-1
      end do
      if (i.le.j) then
        if (vector(i).ne.vector(j)) then
          tmp=vector(i)
          vector(i)=vector(j)
          vector(j)=tmp
        end if
        i=i+1
        j=j-1
      end if
      if (i.gt.j) exit
    end do
    if (start.lt.j) call fbem_quicksort(start,j,length,vector)
    if (i.lt.finish) call fbem_quicksort(i,finish,length,vector)
  end subroutine fbem_quicksort

  !! Standard quicksort routine for sorting each column in ascending order with preference
  !! in the column order: 1, 2, 3, ...
  recursive subroutine fbem_quicksort_matrix(start,finish,length,columns,matrix)
    implicit none
    integer :: start
    integer :: finish
    integer :: length
    integer :: columns
    integer :: matrix(length,columns)
    integer :: x(columns)
    integer :: l
    integer :: i
    integer :: j
    integer :: p
    integer :: k
    integer :: tmp
    i=start
    j=finish
    x=matrix((i+j)/2,:)
    do
      do
        l=i
        do k=1,columns
          if (matrix(i,k).lt.x(k)) then
            i=i+1
            exit
          else if (matrix(i,k).eq.x(k)) then
            cycle
          else
            exit
          end if
        end do
        if (l.eq.i) exit
      end do
      do
        l=j
        do k=1,columns
          if (matrix(j,k).gt.x(k)) then
            j=j-1
            exit
          else if (matrix(j,k).eq.x(k)) then
            cycle
          else
            exit
          end if
        end do
        if (l.eq.j) exit
      end do
      if (i.le.j) then
        if (sum(abs(matrix(i,:)-matrix(j,:))).ne.0) then
          do k=1,columns
            tmp=matrix(i,k)
            matrix(i,k)=matrix(j,k)
            matrix(j,k)=tmp
          end do
        end if
        i=i+1
        j=j-1
      end if
      if (i.gt.j) exit
    end do
    if (start.lt.j)  call fbem_quicksort_matrix(start,j,length,columns,matrix)
    if (i.lt.finish) call fbem_quicksort_matrix(i,finish,length,columns,matrix)
  end subroutine fbem_quicksort_matrix

  !! Standard quicksort routine for sorting each column in ascending order with preference
  !! in the column order: 1, 2, 3, ...
  recursive subroutine fbem_quicksort_matrix_r64(start,finish,length,columns,matrix)
    implicit none
    integer(kind=int64) :: start
    integer(kind=int64) :: finish
    integer(kind=int64) :: length
    integer             :: columns
    real(kind=real64)   :: matrix(length,columns)
    real(kind=real64)   :: x(columns)
    integer(kind=int64) :: l
    integer(kind=int64) :: i
    integer(kind=int64) :: j
    integer(kind=int64) :: p
    integer(kind=int64) :: k
    real(kind=real64)   :: tmp
    i=start
    j=finish
    x=matrix((i+j)/2,:)
    do
      do
        l=i
        do k=1,columns
          if (matrix(i,k).lt.x(k)) then
            i=i+1
            exit
          else if (matrix(i,k).eq.x(k)) then
            cycle
          else
            exit
          end if
        end do
        if (l.eq.i) exit
      end do
      do
        l=j
        do k=1,columns
          if (matrix(j,k).gt.x(k)) then
            j=j-1
            exit
          else if (matrix(j,k).eq.x(k)) then
            cycle
          else
            exit
          end if
        end do
        if (l.eq.j) exit
      end do
      if (i.le.j) then
        if (sum(abs(matrix(i,:)-matrix(j,:))).ne.0.d0) then
          do k=1,columns
            tmp=matrix(i,k)
            matrix(i,k)=matrix(j,k)
            matrix(j,k)=tmp
          end do
        end if
        i=i+1
        j=j-1
      end if
      if (i.gt.j) exit
    end do
    if (start.lt.j)  call fbem_quicksort_matrix_r64(start,j,length,columns,matrix)
    if (i.lt.finish) call fbem_quicksort_matrix_r64(i,finish,length,columns,matrix)
  end subroutine fbem_quicksort_matrix_r64

  ! Search an integer in a sorted vector, if it finds it, it gives the position
  ! if it is not there, it gives the insertion position.
  recursive subroutine fbem_binary_tree(start,finish,x,length,vector,kins,kp)
    implicit none
    integer :: start
    integer :: finish
    integer :: x
    integer :: length
    integer :: vector(length)
    integer :: kins
    integer :: kp
    integer :: k, mid
    if (start.le.finish) then
      mid=start+(finish-start)/2;
      if (vector(mid).eq.x) then
        kins=0
        kp=mid
        return
      else
        if (vector(mid).gt.x) then
          call fbem_binary_tree(start,mid-1,x,length,vector,kins,kp)
        else
          call fbem_binary_tree(mid+1,finish,x,length,vector,kins,kp)
        end if
      end if
    else
      kins=start
      kp=0
      return
    end if
  end subroutine fbem_binary_tree

  ! Search an integer in a sorted matrix, if it finds it, it gives the position
  ! if it is not there, it gives the insertion position.
  recursive subroutine fbem_binary_tree_matrix2(start,finish,x,length,matrix,kins,kp)
    implicit none
    integer :: start
    integer :: finish
    integer :: x(2)
    integer :: length
    integer :: matrix(length,2)
    integer :: kins
    integer :: kp
    integer :: k, mid
    if (start.le.finish) then
      mid=start+(finish-start)/2;
      if (matrix(mid,1).eq.x(1)) then
        if (matrix(mid,2).eq.x(2)) then
          kins=0
          kp=mid
          return
        else
          if (matrix(mid,2).gt.x(2)) then
            call fbem_binary_tree_matrix2(start,mid-1,x,length,matrix,kins,kp)
          else
            call fbem_binary_tree_matrix2(mid+1,finish,x,length,matrix,kins,kp)
          end if
        end if
      else
        if (matrix(mid,1).gt.x(1)) then
          call fbem_binary_tree_matrix2(start,mid-1,x,length,matrix,kins,kp)
        else
          call fbem_binary_tree_matrix2(mid+1,finish,x,length,matrix,kins,kp)
        end if
      end if
    else
      kins=start
      kp=0
      return
    end if
  end subroutine fbem_binary_tree_matrix2

  subroutine fbem_unique(v,w)
   implicit none
   integer, intent(in) :: v(:)
   integer, allocatable, intent(out) :: w(:)
   integer :: i,c
   logical :: mask(size(v))
   mask=.false.
   do i=1,size(v)
     c=count(v(i).eq.v)
     if (c.eq.1) then
       mask(i)=.true.
     else
       if (.not.any(v(i).eq.v.and.mask)) mask(i)=.true.
     end if
   end do
   allocate(w(count(mask)))
   w=pack(v,mask)
  end subroutine fbem_unique

  !! This subroutine solves a Linear System of Equations (LSE) using several subroutines from the Lapack library. It performs
  !! scaling and refines the solution. It requires a complete copy of A, so be careful with memory requirements. For complex
  !! numbers of <tt>kind=real64</tt> and multiple right hand sides. As usual, A is modified (scaled and factorized) at exit,
  !! and b contains the solution.
  subroutine fbem_solve_lse_c(n_dof,A,n_rhs,b)
    ! No implicit variables are allowed
    implicit none
    ! I/O
    integer              :: n_dof           !! Number of degrees of freedom
    complex(kind=real64) :: A(n_dof,n_dof)  !! Matrix A
    integer              :: n_rhs           !! Number of right hand sides
    complex(kind=real64) :: b(n_dof,n_rhs)  !! Vector b
    ! Local
    integer              :: ipiv(n_dof)         ! Row order
    character(len=1)     :: equed               ! Scaling applied
    real(kind=real64)    :: r(n_dof)            ! Row scaling factors
    real(kind=real64)    :: c(n_dof)            ! Column scaling factors
    complex(kind=real64) :: Ao(n_dof,n_dof)     ! Matrix A unfactorized
    complex(kind=real64) :: bcopy(n_dof,n_rhs)  ! Copy of vector b required by the refine stage
    real(kind=real64)    :: rowcnd
    real(kind=real64)    :: colcnd
    real(kind=real64)    :: amax
    real(kind=real64)    :: summ
    character(len=1)     :: norm
    real(kind=real64)    :: anorm
    real(kind=real64)    :: rcond
    real(kind=real64)    :: rwork(n_dof)
    complex(kind=real64) :: work(4*n_dof)
    integer              :: info
    character(len=1)     :: trans
    integer              :: i, j
    real(kind=real64)    :: ferr(n_rhs)
    real(kind=real64)    :: berr(n_rhs)
    !
    ! Scaling stage
    !
    call zgeequ(n_dof,n_dof,A,n_dof,r,c,rowcnd,colcnd,amax,info)
    if (info.gt.0) then
      write(error_unit,*) 'zgeequ info: ', info
      stop
    end if
    call zlaqge(n_dof,n_dof,A,n_dof,r,c,rowcnd,colcnd,amax,equed)
    Ao=A
    if ((equed.eq.'R').or.(equed.eq.'B')) then
      do i=1,n_rhs
        b(:,i)=r(:)*b(:,i)
      end do
    end if
    !
    ! Factorization stage
    !
    call zgetrf(n_dof,n_dof,A,n_dof,ipiv,info)
    if (info.gt.0) then
      write(error_unit,*) 'zgetrf info: ', info
      stop
    end if
    !
    ! Solve stage
    !
    bcopy=b
    trans='N'
    call zgetrs(trans,n_dof,n_rhs,A,n_dof,ipiv,b,n_dof,info)
    if (info.gt.0) then
      write(error_unit,*) 'zgetrs info: ', info
      stop
    end if
    !
    ! Refine stage
    !
    call zgerfs(trans,n_dof,n_rhs,Ao,n_dof,A,n_dof,ipiv,bcopy,n_dof,b,n_dof,ferr,berr,work,rwork,info)
    if (info.gt.0) then
      write(error_unit,*) 'zgerfs info: ', info
      stop
    end if
    !
    ! Unscaling of b (solution)
    !
    if ((equed.eq.'C').or.(equed.eq.'B')) then
      do i=1,n_rhs
        b(:,i)=c(:)*b(:,i)
      end do
    end if
  end subroutine fbem_solve_lse_c

!* ----------------------------------------------------------------------------
!* Numerical diagonalization of 3x3 matrcies
!* Copyright (C) 2006  Joachim Kopp
!* ----------------------------------------------------------------------------
!* This library is free software; you can redistribute it and/or
!* modify it under the terms of the GNU Lesser General Public
!* License as published by the Free Software Foundation; either
!* version 2.1 of the License, or (at your option) any later version.
!*
!* This library is distributed in the hope that it will be useful,
!* but WITHOUT ANY WARRANTY; without even the implied warranty of
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!* Lesser General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with this library; if not, write to the Free Software
!* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
!* ----------------------------------------------------------------------------

! ----------------------------------------------------------------------------
      SUBROUTINE DSYEVJ3(A, Q, W)
! ----------------------------------------------------------------------------
! Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
! matrix A using the Jacobi algorithm.
! The upper triangular part of A is destroyed during the calculation,
! the diagonal elements are read but not destroyed, and the lower
! triangular elements are not referenced at all.
! ----------------------------------------------------------------------------
! Parameters:
!   A: The symmetric input matrix
!   Q: Storage buffer for eigenvectors
!   W: Storage buffer for eigenvalues
! ----------------------------------------------------------------------------
!     .. Arguments ..
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION Q(3,3)
      DOUBLE PRECISION W(3)

!     .. Parameters ..
      INTEGER          N
      PARAMETER        ( N = 3 )

!     .. Local Variables ..
      DOUBLE PRECISION SD, SO
      DOUBLE PRECISION S, C, T
      DOUBLE PRECISION G, H, Z, THETA
      DOUBLE PRECISION THRESH
      INTEGER          I, X, Y, R

!     Initialize Q to the identitity matrix
!     --- This loop can be omitted if only the eigenvalues are desired ---
      DO 10 X = 1, N
        Q(X,X) = 1.0D0
        DO 11, Y = 1, X-1
          Q(X, Y) = 0.0D0
          Q(Y, X) = 0.0D0
   11   CONTINUE
   10 CONTINUE

!     Initialize W to diag(A)
      DO 20 X = 1, N
        W(X) = A(X, X)
   20 CONTINUE

!     Calculate SQR(tr(A))
      SD = 0.0D0
      DO 30 X = 1, N
        SD = SD + ABS(W(X))
   30 CONTINUE
      SD = SD**2

!     Main iteration loop
      DO 40 I = 1, 50
!       Test for convergence
        SO = 0.0D0
        DO 50 X = 1, N
          DO 51 Y = X+1, N
            SO = SO + ABS(A(X, Y))
   51     CONTINUE
   50   CONTINUE
        IF (SO .EQ. 0.0D0) THEN
          RETURN
        END IF

        IF (I .LT. 4) THEN
          THRESH = 0.2D0 * SO / N**2
        ELSE
          THRESH = 0.0D0
        END IF

!       Do sweep
        DO 60 X = 1, N
          DO 61 Y = X+1, N
            G = 100.0D0 * ( ABS(A(X, Y)) )
            IF ( I .GT. 4 .AND. ABS(W(X)) + G .EQ. ABS(W(X)).AND. ABS(W(Y)) + G .EQ. ABS(W(Y)) ) THEN
              A(X, Y) = 0.0D0
            ELSE IF (ABS(A(X, Y)) .GT. THRESH) THEN
!             Calculate Jacobi transformation
              H = W(Y) - W(X)
              IF ( ABS(H) + G .EQ. ABS(H) ) THEN
                T = A(X, Y) / H
              ELSE
                THETA = 0.5D0 * H / A(X, Y)
                IF (THETA .LT. 0.0D0) THEN
                  T = -1.0D0 / (SQRT(1.0D0 + THETA**2) - THETA)
                ELSE
                  T = 1.0D0 / (SQRT(1.0D0 + THETA**2) + THETA)
                END IF
              END IF

              C = 1.0D0 / SQRT( 1.0D0 + T**2 )
              S = T * C
              Z = T * A(X, Y)

!             Apply Jacobi transformation
              A(X, Y) = 0.0D0
              W(X)    = W(X) - Z
              W(Y)    = W(Y) + Z
              DO 70 R = 1, X-1
                T       = A(R, X)
                A(R, X) = C * T - S * A(R, Y)
                A(R, Y) = S * T + C * A(R, Y)
   70         CONTINUE
              DO 80, R = X+1, Y-1
                T       = A(X, R)
                A(X, R) = C * T - S * A(R, Y)
                A(R, Y) = S * T + C * A(R, Y)
   80         CONTINUE
              DO 90, R = Y+1, N
                T       = A(X, R)
                A(X, R) = C * T - S * A(Y, R)
                A(Y, R) = S * T + C * A(Y, R)
   90         CONTINUE

!             Update eigenvectors
!             --- This loop can be omitted if only the eigenvalues are desired ---
              DO 100, R = 1, N
                T       = Q(R, X)
                Q(R, X) = C * T - S * Q(R, Y)
                Q(R, Y) = S * T + C * Q(R, Y)
  100         CONTINUE
            END IF
   61     CONTINUE
   60   CONTINUE
   40 CONTINUE

      PRINT *, "DSYEVJ3: No convergence."

      END SUBROUTINE
! End of subroutine DSYEVJ3

end module fbem_numerical
