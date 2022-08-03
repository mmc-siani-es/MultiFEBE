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
!! TO DO:
!!   - The way quadrature rules are stored is not efficient in the sense they are arranged as matrices, having the first index as
!!     the integration point index, and the second index as the quadrature rule index. Hence, almost half of the matrix is filled
!!     with zeros. However, this way of storing them is very easy to use. This would be changed to a vector-like storage, where
!!     an offset vector have to be defined for each set of quadrature rules. This would change the way other routines have access
!!     to the quadrature rules.
!!
!! <b>This module implements quadrature rules for 1D and 2D integration.</b>
!!
!! <h2>DESCRIPTION</h2>
!!
!! For each quadrature there are:
!!   -A scalar indicating the number of rules: <tt>?_nr</tt>.
!!   -A vector indicating the number of gauss points for each rule: <tt>?_n</tt>.
!!   -Matrices indicating the coordinates for each point and rule: <tt>?_xi</tt> (1D), or <tt>?_xi1</tt> and  <tt>?_xi2</tt> (2D).
!!   -A matrix indicating the weights for each point and rule: <tt>?_w</tt>.
!!
!! Coordinates and weights of a rule are stored in columns in their matrices in order to take the most of column-major order of
!! Fortran arrays. Points and weights are declared in <tt>quad_rules.f90</tt>, but initialized in resource files
!! <tt>quad_rules_rc/?.rc</tt> via <tt>data</tt> statements.
!!
!! There are 7 quadratures rules defined:
!!
!! - Gauss-Legendre (weight function is <tt>1</tt>) in <tt>[-1,1]</tt> domain (<tt>gl11</tt>):
!!   - 32 rules, from <tt>n=1</tt> to <tt>n=32</tt>, from order 1 to order 63.
!!   - http://processingjs.nihongoresources.com/bezierinfo/legendre-gauss-values.php
!!
!! - Gauss-Legendre (weight function is <tt>1</tt>) in <tt>[0,1]</tt> domain (<tt>gl01</tt>):
!!   - 32 rules, from <tt>n=1</tt> to <tt>n=32</tt>, from order 1 to order 63.
!!   - http://processingjs.nihongoresources.com/bezierinfo/legendre-gauss-values.php
!!
!! - Gauss-Jacobi (weight function is <tt>1-x</tt>) in <tt>[0,1]</tt> domain (<tt>gj01</tt>):
!!   - 32 rules, from <tt>n=1</tt> to <tt>n=32</tt>, from order 1 to order 63.
!!   - John Burkardt's Website (http://people.sc.fsu.edu/~jburkardt/)
!!
!! - Gauss-Ln (weight function is <tt>ln(1/x)</tt>) in <tt>[0,1]</tt> domain (<tt>gln01</tt>):
!!   - 32 rules, from <tt>n=1</tt> to <tt>n=32</tt>, from order 1 to order 63.
!!   - W. Gautschi, "OPQ: A MATLAB suite of programs for generating orthogonal polynomials and related quadrature rules".
!!   - http://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
!!
!! - Gauss-LnPoly (weight function is <tt>ln(1/x)</tt>) in <tt>[0,1]</tt> domain (<tt>gln01</tt>):
!!   - 7 rules, from <tt>n=1</tt> to <tt>n=7</tt>, from order 1 to order 63.
!!   - Crow, Quadrature of integrands with a logarithmic singularity, 1993.
!!
!! - Wandzura rules for triangle region (<tt>wantri</tt>):
!!   - Weights sum is 0.5, so the reference triangle area assumed is 0.5.
!!   - 30 rules, from <tt>n=1</tt> to <tt>n=175</tt>, from order 1 to order 30.
!!   - S. Wandzura, "Symmetric quadrature rules on a triangle". Computers and Mathematics with Applications, Vol. 45, 1829-1840,
!!     2003.
!!   - S. Wandzura, private communications October 2012 (Thank you very much).
!!   - John Burkardt's Website (http://people.sc.fsu.edu/~jburkardt/)
!!
!! WARNING: Be careful when compiling this module, some resource files contain a lot of continuation lines. Also, the Fortran 2003
!!          module <tt>iso_fortran_env</tt> is used, but not all compilers support it. At least GNU Fortran (gfortran) 4.6 and
!!          Intel Fortran (ifort) 13.0 can compile the module.
module fbem_quad_rules

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all are private
  private

  ! Gauss-Legendre [-1,1]: 16516 bytes
  integer, public           :: gl11_nr            !! Gauss-Legendre [-1,1]: number of rules defined
  integer, public           :: gl11_n(32)         !! Gauss-Legendre [-1,1]: number of points for each rule
  real(kind=real64), public :: gl11_xi(32,32)     !! Gauss-Legendre [-1,1]: coordinates for each rule
  real(kind=real64), public :: gl11_w(32,32)      !! Gauss-Legendre [-1,1]: weights for each rule
  ! Gauss-Legendre [0,1]: 16516 bytes
  integer, public           :: gl01_nr            !! Gauss-Legendre [0,1]: number of rules defined
  integer, public           :: gl01_n(32)         !! Gauss-Legendre [0,1]: number of points for each rule
  real(kind=real64), public :: gl01_xi(32,32)     !! Gauss-Legendre [0,1]: coordinates for each rule
  real(kind=real64), public :: gl01_w(32,32)      !! Gauss-Legendre [0,1]: weights for each rule
  ! Gauss-Jacobi [0,1]: 16516 bytes
  integer, public           :: gj01_nr            !! Gauss-Jacobi [0,1]: number of rules defined
  integer, public           :: gj01_n(32)         !! Gauss-Jacobi [0,1]: number of points for each rule
  real(kind=real64), public :: gj01_xi(32,32)     !! Gauss-Jacobi [0,1]: coordinates for each rule
  real(kind=real64), public :: gj01_w(32,32)      !! Gauss-Jacobi [0,1]: weights for each rule
  ! Gauss-Ln [0,1]: 16516 bytes
  integer, public           :: gln01_nr           !! Gauss-Ln [0,1]: number of rules defined
  integer, public           :: gln01_n(32)        !! Gauss-Ln [0,1]: number of points for each rule
  real(kind=real64), public :: gln01_xi(32,32)    !! Gauss-Ln [0,1]: coordinates for each rule
  real(kind=real64), public :: gln01_w(32,32)     !! Gauss-Ln [0,1]: weights for each rule
  ! Gauss-Crow [0,1]: 816 bytes
  integer, public           :: gc01_nr            !! Gauss-Crow [0,1]: number of rules defined
  integer, public           :: gc01_n(7)          !! Gauss-Crow [0,1]: number of points for each rule
  real(kind=real64), public :: gc01_xi(7,7)       !! Gauss-Crow [0,1]: coordinates for each rule
  real(kind=real64), public :: gc01_w(7,7)        !! Gauss-Crow [0,1]: weights for each rule
  ! Wandzura triangle: 126844 bytes
  integer, public           :: wantri_nr          !! Wandzura triangle: number of rules defined
  integer, public           :: wantri_n(30)       !! Wandzura triangle: number of points for each rule
  real(kind=real64), public :: wantri_xi1(176,30) !! Wandzura triangle: coordinates 1 for each rule
  real(kind=real64), public :: wantri_xi2(176,30) !! Wandzura triangle: coordinates 2 for each rule
  real(kind=real64), public :: wantri_w(176,30)   !! Wandzura triangle: weights for each rule

  ! Total memory required: 192908 bytes

  ! Data initialization
  ! Gauss-Legendre [-1,1]
  include 'resources_quad_rules/gl11.rc'
  ! Gauss-Legendre [0,1]
  include 'resources_quad_rules/gl01.rc'
  ! Gauss-Jacobi [0,1]
  include 'resources_quad_rules/gj01.rc'
  ! Gauss-Ln [0,1]
  include 'resources_quad_rules/gln01.rc'
  ! Gauss-Crow [0,1]
  include 'resources_quad_rules/gc01.rc'
  ! Wandzura triangle
  include 'resources_quad_rules/wantri.rc'

end module fbem_quad_rules
