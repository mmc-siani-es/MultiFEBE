! ---------------------------------------------------------------------
! Copyright (C) 2014-2023 Universidad de Las Palmas de Gran Canaria:
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
!! <b> Subroutine that build the symmetry planes for a given bodyload element, taking into
!! account the problem symmetry planes and deactivating those symmetry planes that
!! the element belongs to. </b>

subroutine build_symplane_bodyload_elements(se_int,se_n_symplanes,se_n_symelements,se_symplane_m,se_symplane_s,se_symplane_t,se_symplane_r)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures

  ! Module of problem variables
  use problem_variables

  ! No implicit variables
  implicit none

  ! I/O variables
  integer                :: se_int                   !! Selected element
  integer                :: se_n_symplanes           !! Number of symmetry planes: 0 ,1, 2 or 3.
  integer                :: se_n_symelements         !! Number of symmetrical elements
  real(kind=real64)      :: se_symplane_m(3,3)       !! Multiplier for each sym. plane for geometrical reflection: 1. or -1.
  real(kind=real64)      :: se_symplane_s(3)         !! Multiplier for each sym. plane for scalar variables reflection: 1. or -1.
  real(kind=real64)      :: se_symplane_t(3,3)       !! Multiplier for each sym. plane for translational vectorial variables reflection: 1. or -1.
  real(kind=real64)      :: se_symplane_r(3,3)       !! Multiplier for each sym. plane for rotational vectorial variables reflection: 1. or -1.
  ! Local variables
  integer                :: ks, ksj, se_ks
  logical                :: deactivate

  !
  ! Build the symmetry plane configuration for the current element
  !
  !
  ! For elements not lying on a symmetry plane, mirroring with respect to all of them are needed, so copy
  ! use the problem symmetry plane configuration. This is the situation in most cases.
  !
  if (element(se_int)%n_symplanes.eq.0) then
    se_n_symplanes  =n_symplanes
    se_n_symelements=2**n_symplanes
    se_symplane_m   =symplane_m
    se_symplane_s   =symplane_s
    se_symplane_t   =symplane_t
    se_symplane_r   =symplane_r
  !
  ! For elements lying on any symmetry plane, mirroring with respecto to these is deactivated.
  !
  else
    se_n_symplanes  = n_symplanes-element(se_int)%n_symplanes
    se_n_symelements=2**se_n_symplanes
    se_ks=0
    do ks=1,n_symplanes
      deactivate=.false.
      do ksj=1,element(se_int)%n_symplanes
        if (element(se_int)%symplane(ksj).eq.ks) deactivate=.true.
      end do
      if (.not.deactivate) then
          se_ks=se_ks+1
          se_symplane_m(:,se_ks)=symplane_m(:,ks)
          se_symplane_s(  se_ks)=symplane_s(  ks)
          se_symplane_t(:,se_ks)=symplane_t(:,ks)
          se_symplane_r(:,se_ks)=symplane_r(:,ks)
      end if
    end do
  end if

end subroutine build_symplane_bodyload_elements
