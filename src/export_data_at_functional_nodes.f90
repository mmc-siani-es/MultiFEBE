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
!! <b> Subroutine that export the data at geometrical nodes. </b>

subroutine export_data_at_functional_nodes

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_string_handling
  use fbem_data_structures

  ! Module of problem variables
  use problem_variables

  ! Local variables
  implicit none
  integer                                 :: ke, kn, i
  integer                                 :: tmp_unit(7)
  character(len=fbem_filename_max_length) :: tmp_filename

  ! Starting message
  if (verbose_level.ge.1)  write(output_unit,'(a)') 'Exporting element data at functional nodes ...'

  ! Open files
  do i=1,6
    tmp_unit(i)=fbem_get_valid_unit()
    select case (i)
      case (1)
        tmp_filename=trim(output_filename)//'.x_fn_element.dat'
      case (2)
        tmp_filename=trim(output_filename)//'.n_fn_element.dat'
      case (3)
        tmp_filename=trim(output_filename)//'.x_fn_node.dat'
      case (4)
        tmp_filename=trim(output_filename)//'.n_fn_node.dat'
      case (5)
        tmp_filename=trim(output_filename)//'.t1_fn_node.dat'
      case (6)
        tmp_filename=trim(output_filename)//'.t2_fn_node.dat'
    end select
    call fbem_trim2b(tmp_filename)
    open(unit=tmp_unit(i),file=trim(tmp_filename),action='write',recl=fbem_file_record_length)
    if (verbose_level.ge.2) write(output_unit,'(a,a,a)') ' Opening "', trim(tmp_filename), '"'
  end do

  ! ============================================================================================================================== !
  ! ELEMENT-WISE VECTORS                                                                                                           !
  ! ============================================================================================================================== !

  do ke=1,n_elements
    do kn=1,element(ke)%n_nodes
      select case (problem%n)
        case (2)
          ! Position vectors of functional nodes
          write(tmp_unit(1),'(2e25.16)') element(ke)%x_fn(1,kn), element(ke)%x_fn(2,kn)
          ! Unit normals at functional nodes
          if (element(ke)%n_dimension.eq.1) then
            write(tmp_unit(2),'(4e25.16)') element(ke)%x_fn(1,kn), element(ke)%x_fn(2,kn),&
                                           vectors_scale_factor*(element(ke)%n_fn(1,kn)),&
                                           vectors_scale_factor*(element(ke)%n_fn(2,kn))
          end if
        case (3)
          ! Position vectors of functional nodes
          write(tmp_unit(1),'(3e25.16)') element(ke)%x_fn(1,kn), element(ke)%x_fn(2,kn), element(ke)%x_fn(3,kn)
          ! Unit normals at functional nodes
          if (element(ke)%n_dimension.eq.2) then
            write(tmp_unit(2),'(6e25.16)') element(ke)%x_fn(1,kn), element(ke)%x_fn(2,kn), element(ke)%x_fn(3,kn),&
                                           vectors_scale_factor*(element(ke)%n_fn(1,kn)),&
                                           vectors_scale_factor*(element(ke)%n_fn(2,kn)),&
                                           vectors_scale_factor*(element(ke)%n_fn(3,kn))
          end if
      end select
    end do
  end do

  ! ============================================================================================================================== !
  ! NODE-WISE VECTORS                                                                                                              !
  ! ============================================================================================================================== !

  do kn=1,n_nodes
    select case (problem%n)
      case (2)
        ! Position vectors of geometrical nodes
        write(tmp_unit(3),'(2e25.16)') node(kn)%x_fn(1), node(kn)%x_fn(2)
        ! Unit normals at functional nodes
        if (element(node(kn)%element(1))%n_dimension.eq.1) then
          write(tmp_unit(4),'(4e25.16)') node(kn)%x_fn(1), node(kn)%x_fn(2),&
                                         (node(kn)%n_fn(1)),&
                                         (node(kn)%n_fn(2))
        end if
        write(tmp_unit(5),'(4e25.16)') node(kn)%x_fn(1), node(kn)%x_fn(2),&
                                       (node(kn)%t1_fn(1)),&
                                       (node(kn)%t1_fn(2))
        if (element(node(kn)%element(1))%n_dimension.eq.2) then
          write(tmp_unit(6),'(4e25.16)') node(kn)%x_fn(1), node(kn)%x_fn(2),&
                                         (node(kn)%t2_fn(1)),&
                                         (node(kn)%t2_fn(2))
        end if
      case (3)
        ! Position vectors of geometrical nodes
        write(tmp_unit(3),'(3e25.16)') node(kn)%x_fn(1), node(kn)%x_fn(2), node(kn)%x_fn(3)
        ! Unit normals at functional nodes
        if (element(node(kn)%element(1))%n_dimension.eq.2) then
          write(tmp_unit(4),'(6e25.16)') node(kn)%x_fn(1), node(kn)%x_fn(2), node(kn)%x_fn(3),&
                                         (node(kn)%n_fn(1)),&
                                         (node(kn)%n_fn(2)),&
                                         (node(kn)%n_fn(3))
        end if
        write(tmp_unit(5),'(6e25.16)') node(kn)%x_fn(1), node(kn)%x_fn(2), node(kn)%x_fn(3),&
                                       (node(kn)%t1_fn(1)),&
                                       (node(kn)%t1_fn(2)),&
                                       (node(kn)%t1_fn(3))
        if (element(node(kn)%element(1))%n_dimension.eq.2) then
          write(tmp_unit(6),'(6e25.16)') node(kn)%x_fn(1), node(kn)%x_fn(2), node(kn)%x_fn(3),&
                                         (node(kn)%t2_fn(1)),&
                                         (node(kn)%t2_fn(2)),&
                                         (node(kn)%t2_fn(3))
        end if
    end select
  end do

  ! Close files
  do i=1,6
    close(tmp_unit(i))
  end do

  ! Ending message
  if (verbose_level.ge.1) write(output_unit,'(a)') 'done.'

end subroutine export_data_at_functional_nodes
