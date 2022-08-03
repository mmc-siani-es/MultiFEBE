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

subroutine transformation_element_to_quarter_point(sn)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_symmetry
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! No implicit variables are allowed
  implicit none

  ! I/O
  integer                        :: sn
  ! Local
  integer                        :: knl
  integer                        :: ke, se
  integer                        :: kc, knc, snc
  integer                        :: midnode
  real(kind=real64)              :: ratio, length
  real(kind=real64), allocatable :: xnodes(:,:)
  real(kind=real64)              :: r(problem%n)
  character(len=fbem_fmtstr)     :: fmtstr            ! String used for write format string

  ! Length between the tip and the tip node with respect to the total length
  ratio=0.25d0+1.d-6

  ! Loop through the elements of the node
  do ke=1,node(sn)%n_elements
    se=node(sn)%element(ke)
    knl=node(sn)%element_node_iid(ke)
    allocate (xnodes(problem%n,element(se)%n_nodes))
    do kc=1,problem%n
      do knc=1,element(se)%n_nodes
        snc=element(se)%node(knc)
        xnodes(kc,knc)=node(snc)%x(kc)
      end do
    end do
    select case (element(se)%type)
      !
      ! LINE3
      !
      case (fbem_line3)
        if (problem%n.eq.3) then
          call fbem_error_message(error_unit,0,'element',se,&
                                  'a line3 element cannot be transformed into a quarter-point element for 3D problems')
        end if
        ! Calculate the vector between node 1 and 2 and length
        length=0.0d0
        do kc=1,problem%n
          r(kc)=xnodes(kc,2)-xnodes(kc,1)
          length=length+r(kc)**2
        end do
        length=dsqrt(length)
        midnode=element(se)%node(3)
        select case (knl)
          case (1)
            do kc=1,problem%n
              node(midnode)%x(kc)=node(sn)%x(kc)+ratio*r(kc)
            end do
          case (2)
            do kc=1,problem%n
              node(midnode)%x(kc)=node(sn)%x(kc)-ratio*r(kc)
            end do
          case default
            call fbem_error_message(error_unit,0,'node',node(midnode)%id,&
                                  'this node (mid-node) cannot be in the crack front')
        end select
        if (verbose_level.ge.3) then
          write(fmtstr,*) '(3x,a,i8,a,',problem%n,'e25.16)'
          call fbem_trimall(fmtstr)
          write(output_unit,fmtstr) 'New position for node', node(midnode)%id, ':', (node(midnode)%x(kc),kc=1,problem%n)
        end if
      case default
        stop 'this element is not available for quarter-point transformation'
    end select
    deallocate (xnodes)
  end do

end subroutine transformation_element_to_quarter_point
