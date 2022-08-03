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

subroutine export_data_at_geometrical_nodes

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_string_handling
  use fbem_data_structures
  use fbem_shape_functions

  ! Module of problem variables
  use problem_variables

  ! Local variables
  implicit none
  integer                                 :: ke, se, kn, sn, knj, i, k
  integer                                 :: tmp_unit(3)
  character(len=fbem_filename_max_length) :: tmp_filename
  character(len=fbem_fmtstr)              :: fmtstr
  integer, allocatable                    :: drawline_point(:)
  integer                                 :: n_drawlines
  real(kind=real64)                       :: n(problem%n)

  ! Starting message
  if (verbose_level.ge.1)  write(output_unit,'(a)') 'Exporting geometrical of data mesh ...'

  ! Open files
  do i=1,3
    tmp_unit(i)=fbem_get_valid_unit()
    select case (i)
      case (1)
        tmp_filename=trim(output_filename)//'.egn.dat'
      case (2)
        tmp_filename=trim(output_filename)//'.egplines.dat'
      case (3)
        tmp_filename=trim(output_filename)//'.ngn.dat'
    end select
    call fbem_trim2b(tmp_filename)
    open(unit=tmp_unit(i),file=trim(tmp_filename),action='write',recl=fbem_file_record_length)
    if (verbose_level.ge.2) write(output_unit,'(a,a,a)') ' Opening "', trim(tmp_filename), '"'
  end do

  ! ======================================= !
  ! ELEMENT-WISE GEOMETRIC VECTORS AT NODES !
  ! ======================================= !

  ! File format:
  !
  ! - 2D:
  ! <element id> <element dimension> <part id> <node id> <node index> <x> <t1> <t2> <n>
  !
  ! - 3D:
  ! <element id> <element dimension> <part id> <node id> <node index> <x> <t1> <t2> <t3> <n>
  !
  ! Note: vectors (x,t1,t2,t3,n) have n (2 or 3) components, and non-existant vectors are null.

  select case (problem%n)
    case (2)
      write(fmtstr,*) '(i',fbem_nchar_int(element_eid_max)+1,',i2,i',fbem_nchar_int(part_eid_max)+1,&
                      ',i',fbem_nchar_int(node_eid_max)+1,',i3,8',fmt_real,')'
      call fbem_trimall(fmtstr)
      do ke=1,n_elements
        select case (element(ke)%n_dimension)
          case (0)
          case (1)
            do kn=1,element(ke)%n_nodes
              sn=element(ke)%node(kn)
                write(tmp_unit(1),fmtstr) element(ke)%id, element(ke)%n_dimension, element(ke)%part, node(sn)%id, kn, &
                                          element(ke)%x_gn(:,kn), &
                                          element(ke)%t1_gn(:,kn), &
                                          [0.,0.], &
                                          element(ke)%n_gn(:,kn)
            end do
          case (2)
            do kn=1,element(ke)%n_nodes
              sn=element(ke)%node(kn)
                write(tmp_unit(1),fmtstr) element(ke)%id, element(ke)%n_dimension, element(ke)%part, node(sn)%id, kn, &
                                          element(ke)%x_gn(:,kn),&
                                          element(ke)%t1_gn(:,kn),&
                                          element(ke)%t2_gn(:,kn),&
                                          [0.,0.]
            end do
        end select
      end do
    case (3)
      write(fmtstr,*) '(i',fbem_nchar_int(element_eid_max)+1,',i2,i',fbem_nchar_int(part_eid_max)+1,&
                      ',i',fbem_nchar_int(node_eid_max)+1,',i3,15',fmt_real,')'
      call fbem_trimall(fmtstr)
      do ke=1,n_elements
        select case (element(ke)%n_dimension)
          case (0)
          case (1)
            do kn=1,element(ke)%n_nodes
              sn=element(ke)%node(kn)
              write(tmp_unit(1),fmtstr) element(ke)%id, element(ke)%n_dimension, element(ke)%part, node(sn)%id, kn, &
                                        element(ke)%x_gn(:,kn), &
                                        element(ke)%t1_gn(:,kn), &
                                        [0.,0.,0.], &
                                        [0.,0.,0.], &
                                        [0.,0.,0.]
            end do
          case (2)
            do kn=1,element(ke)%n_nodes
              sn=element(ke)%node(kn)
              write(tmp_unit(1),fmtstr) element(ke)%id, element(ke)%n_dimension, element(ke)%part, node(sn)%id, kn, &
                                        element(ke)%x_gn(:,kn), &
                                        element(ke)%t1_gn(:,kn), &
                                        element(ke)%t2_gn(:,kn), &
                                        [0.,0.,0.], &
                                        element(ke)%n_gn(:,kn)
            end do
          case (3)
            do kn=1,element(ke)%n_nodes
              sn=element(ke)%node(kn)
              write(tmp_unit(1),fmtstr) element(ke)%id, element(ke)%n_dimension, element(ke)%part, node(sn)%id, kn, &
                                        element(ke)%x_gn(:,kn), &
                                        element(ke)%t1_gn(:,kn), &
                                        element(ke)%t2_gn(:,kn), &
                                        element(ke)%t3_gn(:,kn), &
                                        [0.,0.,0.]
            end do
        end select
      end do
  end select

  ! ========================= !
  ! ELEMENT LINES USING NODES !
  ! ========================= !

  write(fmtstr,*) '(i',fbem_nchar_int(element_eid_max)+1,',i',fbem_nchar_int(part_eid_max)+1,',',problem%n,fmt_real,')'
  call fbem_trim2b(fmtstr)
  do ke=1,n_elements
    select case (element(ke)%type)
      case (fbem_line2)
        n_drawlines=2
        allocate (drawline_point(n_drawlines))
        drawline_point=(/1,2/)
      case (fbem_line3)
        n_drawlines=3
        allocate (drawline_point(n_drawlines))
        drawline_point=(/1,3,2/)
      case (fbem_line4)
        n_drawlines=4
        allocate (drawline_point(n_drawlines))
        drawline_point=(/1,3,4,2/)
      case (fbem_tri3)
        n_drawlines=4
        allocate (drawline_point(n_drawlines))
        drawline_point=(/1,2,3,1/)
      case (fbem_tri6)
        n_drawlines=7
        allocate (drawline_point(n_drawlines))
        drawline_point=(/1,4,2,5,3,6,1/)
      case (fbem_quad4)
        n_drawlines=5
        allocate (drawline_point(n_drawlines))
        drawline_point=(/1,2,3,4,1/)
      case (fbem_quad8,fbem_quad9)
        n_drawlines=9
        allocate (drawline_point(n_drawlines))
        drawline_point=(/1,5,2,6,3,7,4,8,1/)
      case default
        stop 'not implemented yet'
    end select
    do kn=1,n_drawlines
      sn=element(ke)%node(drawline_point(kn))
      write(tmp_unit(2),fmtstr) element(ke)%id, element(ke)%part, (node(sn)%x(k),k=1,problem%n)
    end do
    write(tmp_unit(2),*)
    write(tmp_unit(2),*)
    deallocate (drawline_point)
  end do

  ! =========================== !
  ! NODE-WISE GEOMETRIC VECTORS !
  ! =========================== !

  ! Format:
  !
  ! <node id> <x> <n>
  !
  ! Note: vectors (x,n) have n (2 or 3) components, and non-existant vectors are null.

  write(fmtstr,*) '(i',fbem_nchar_int(node_eid_max)+1,',',2*problem%n,fmt_real,')'
  call fbem_trim2b(fmtstr)
  do kn=1,n_nodes
    n=0.
    if (((node(kn)%dimensional_degree.eq.1).and.(problem%n.eq.2)).or.((node(kn)%dimensional_degree.eq.2).and.(problem%n.eq.3))) then
      do ke=1,node(kn)%n_elements
        se=node(kn)%element(ke)
        knj=node(kn)%element_node_iid(ke)
        n=n+element(se)%n_gn(:,knj)
      end do
      n=n/sqrt(dot_product(n,n))
    end if
    write(tmp_unit(3),fmtstr) node(kn)%id, node(kn)%x, n
  end do

  ! Close files
  do i=1,3
    close(tmp_unit(i))
  end do

  ! Ending message
  if (verbose_level.ge.1) write(output_unit,'(a)') 'done.'

end subroutine export_data_at_geometrical_nodes
