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
  use fbem_gmsh

  ! Module of problem variables
  use problem_variables

  ! Local variables
  implicit none
  integer                                 :: ke, se, kn, sn, knj, i, k, ss, sp, ks, kr
  integer                                 :: global_sn, global_se
  integer                                 :: tmp_unit(3)
  character(len=fbem_filename_max_length) :: tmp_filename
  character(len=fbem_fmtstr)              :: fmtstr
  integer, allocatable                    :: drawline_point(:)
  integer                                 :: n_drawlines
  real(kind=real64)                       :: n(problem%n)
  type(fbem_gmsh_mesh)                    :: tmp_mesh

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
  ! - 3D:d
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

  ! ============ !
  ! SHELLS IN 3D !
  ! ============ !


  if ((n_fe_subregions.gt.0).and.(problem%n.eq.3)) then

    ! TO-DO: include all the shell elements, for now only quad9 and tri6

    ! Initialization
    call tmp_mesh%init
    tmp_mesh%version = '2.2'
    tmp_mesh%filetype = 0
    tmp_mesh%datasize = 8

    ! Count the number of mesh entities to be exported
    ! The physicalnames are exactly the parts of the original mesh, even though
    ! there will be parts with no elements.
    tmp_mesh%n_physicalnames = n_parts
    tmp_mesh%n_nodes = 0
    tmp_mesh%n_elements = 0
    do kr=1,n_regions
      if (region(kr)%class.ne.fbem_fe) cycle
      do ks=1,region(kr)%n_fe_subregions
        ss=region(kr)%fe_subregion(ks)
        sp=fe_subregion(ss)%part
        do ke=1,part(sp)%n_elements
          se=part(sp)%element(ke)
          if (element(se)%n_dimension.eq.2) then
            !
            ! Solo vamos a exportar los TRI6 y QUAD9
            !
            !
            ! QUAD9 => HEX27
            !
            if (element(se)%type.eq.fbem_quad9) then
              tmp_mesh%n_elements = tmp_mesh%n_elements + 1
              tmp_mesh%n_nodes = tmp_mesh%n_nodes + 27
            !
            ! TRI6 => PRISM18
            !
            else if (element(se)%type.eq.fbem_tri6) then
              tmp_mesh%n_elements = tmp_mesh%n_elements + 1
              tmp_mesh%n_nodes = tmp_mesh%n_nodes + 18
            else
              call fbem_warning_message(error_unit,0,'element',element(se)%id,'only tri6 and quad9 shells are exported to a volumetric mesh')
            end if
          end if
        end do
      end do
    end do

    ! Allocation
    allocate(tmp_mesh%physicalname_dim (tmp_mesh%n_physicalnames))
    allocate(tmp_mesh%physicalname_eid (tmp_mesh%n_physicalnames))
    allocate(tmp_mesh%physicalname_name(tmp_mesh%n_physicalnames))
    allocate(tmp_mesh%node_eid(tmp_mesh%n_nodes))
    allocate(tmp_mesh%node_x(3,tmp_mesh%n_nodes))
    allocate(tmp_mesh%element_eid     (tmp_mesh%n_elements))
    allocate(tmp_mesh%element_type    (tmp_mesh%n_elements))
    allocate(tmp_mesh%element_physical(tmp_mesh%n_elements))
    allocate(tmp_mesh%element_node (27,tmp_mesh%n_elements)) ! Ojo, suponemos un # nodos maximo de 27

    global_sn = 1
    global_se = 1

    do kr=1,n_regions
      if (region(kr)%class.ne.fbem_fe) cycle
      do ks=1,region(kr)%n_fe_subregions
        ss=region(kr)%fe_subregion(ks)
        sp=fe_subregion(ss)%part

        tmp_mesh%physicalname_dim (sp)=element(part(sp)%element(1))%n_dimension
        tmp_mesh%physicalname_eid (sp)=part(sp)%id
        tmp_mesh%physicalname_name(sp)=trim(part(sp)%name)

        do ke=1,part(sp)%n_elements
          se=part(sp)%element(ke)
          if (element(se)%n_dimension.eq.2) then
            !
            ! Solo vamos a exportar los TRI6 y QUAD9
            !
            !
            ! QUAD9 => HEX8----------HEX27
            !
            if (element(se)%type.eq.fbem_quad9) then
              ! Bottom plane vertices
              do kn=1,4
                tmp_mesh%node_eid(global_sn+(kn-1))=global_sn+(kn-1)
                tmp_mesh%node_x(:,global_sn+(kn-1))=element(se)%x_gn(:,kn)-0.5d0*element(se)%tv_midnode(3,kn)*element(se)%v_midnode(:,3,kn)
              end do
              ! Top plane vertices
              do kn=1,4
                tmp_mesh%node_eid(global_sn+4+(kn-1))=global_sn+4+(kn-1)
                tmp_mesh%node_x(:,global_sn+4+(kn-1))=element(se)%x_gn(:,kn)+0.5d0*element(se)%tv_midnode(3,kn)*element(se)%v_midnode(:,3,kn)
              end do
              ! Bottom plane edge nodes
              kn = 5
              tmp_mesh%node_eid(global_sn+8)=global_sn+8
              tmp_mesh%node_x(:,global_sn+8)=element(se)%x_gn(:,kn)-0.5d0*element(se)%tv_midnode(3,kn)*element(se)%v_midnode(:,3,kn)
              kn = 6
              tmp_mesh%node_eid(global_sn+11)=global_sn+11
              tmp_mesh%node_x(:,global_sn+11)=element(se)%x_gn(:,kn)-0.5d0*element(se)%tv_midnode(3,kn)*element(se)%v_midnode(:,3,kn)
              kn = 7
              tmp_mesh%node_eid(global_sn+13)=global_sn+13
              tmp_mesh%node_x(:,global_sn+13)=element(se)%x_gn(:,kn)-0.5d0*element(se)%tv_midnode(3,kn)*element(se)%v_midnode(:,3,kn)
              kn = 8
              tmp_mesh%node_eid(global_sn+9)=global_sn+9
              tmp_mesh%node_x(:,global_sn+9)=element(se)%x_gn(:,kn)-0.5d0*element(se)%tv_midnode(3,kn)*element(se)%v_midnode(:,3,kn)
              ! Top plane edge nodes
              kn = 5
              tmp_mesh%node_eid(global_sn+16)=global_sn+16
              tmp_mesh%node_x(:,global_sn+16)=element(se)%x_gn(:,kn)+0.5d0*element(se)%tv_midnode(3,kn)*element(se)%v_midnode(:,3,kn)
              kn = 6
              tmp_mesh%node_eid(global_sn+18)=global_sn+18
              tmp_mesh%node_x(:,global_sn+18)=element(se)%x_gn(:,kn)+0.5d0*element(se)%tv_midnode(3,kn)*element(se)%v_midnode(:,3,kn)
              kn = 7
              tmp_mesh%node_eid(global_sn+19)=global_sn+19
              tmp_mesh%node_x(:,global_sn+19)=element(se)%x_gn(:,kn)+0.5d0*element(se)%tv_midnode(3,kn)*element(se)%v_midnode(:,3,kn)
              kn = 8
              tmp_mesh%node_eid(global_sn+17)=global_sn+17
              tmp_mesh%node_x(:,global_sn+17)=element(se)%x_gn(:,kn)+0.5d0*element(se)%tv_midnode(3,kn)*element(se)%v_midnode(:,3,kn)
              ! Mid-plane vertices
              kn = 1
              tmp_mesh%node_eid(global_sn+10)=global_sn+10
              tmp_mesh%node_x(:,global_sn+10)=element(se)%x_gn(:,kn)
              kn = 2
              tmp_mesh%node_eid(global_sn+12)=global_sn+12
              tmp_mesh%node_x(:,global_sn+12)=element(se)%x_gn(:,kn)
              kn = 3
              tmp_mesh%node_eid(global_sn+14)=global_sn+14
              tmp_mesh%node_x(:,global_sn+14)=element(se)%x_gn(:,kn)
              kn = 4
              tmp_mesh%node_eid(global_sn+15)=global_sn+15
              tmp_mesh%node_x(:,global_sn+15)=element(se)%x_gn(:,kn)
              ! Bottom plane mid node
              kn = 9
              tmp_mesh%node_eid(global_sn+20)=global_sn+20
              tmp_mesh%node_x(:,global_sn+20)=element(se)%x_gn(:,kn)-0.5d0*element(se)%tv_midnode(3,kn)*element(se)%v_midnode(:,3,kn)
              ! Mid-plane edge nodes and interior
              kn = 5
              tmp_mesh%node_eid(global_sn+21)=global_sn+21
              tmp_mesh%node_x(:,global_sn+21)=element(se)%x_gn(:,kn)
              kn = 6
              tmp_mesh%node_eid(global_sn+23)=global_sn+23
              tmp_mesh%node_x(:,global_sn+23)=element(se)%x_gn(:,kn)
              kn = 7
              tmp_mesh%node_eid(global_sn+24)=global_sn+24
              tmp_mesh%node_x(:,global_sn+24)=element(se)%x_gn(:,kn)
              kn = 8
              tmp_mesh%node_eid(global_sn+22)=global_sn+22
              tmp_mesh%node_x(:,global_sn+22)=element(se)%x_gn(:,kn)
              kn = 9
              tmp_mesh%node_eid(global_sn+26)=global_sn+26
              tmp_mesh%node_x(:,global_sn+26)=element(se)%x_gn(:,kn)
              ! Top plane mid node
              kn = 9
              tmp_mesh%node_eid(global_sn+25)=global_sn+25
              tmp_mesh%node_x(:,global_sn+25)=element(se)%x_gn(:,kn)+0.5d0*element(se)%tv_midnode(3,kn)*element(se)%v_midnode(:,3,kn)
              ! Element construction
              tmp_mesh%element_eid     (global_se) = global_se
              tmp_mesh%element_type    (global_se) = 12
              tmp_mesh%element_physical(global_se) = part(sp)%id
              do kn=1,27
                tmp_mesh%element_node    (kn,global_se) = global_sn+(kn-1)
              end do
              global_sn = global_sn + 27
              global_se = global_se + 1
            !
            ! TRI6 => PRISM18
            !
            else if (element(se)%type.eq.fbem_tri6) then
              ! Bottom plane vertices
              do kn=1,3
                tmp_mesh%node_eid(global_sn+(kn-1))=global_sn+(kn-1)
                tmp_mesh%node_x(:,global_sn+(kn-1))=element(se)%x_gn(:,kn)-0.5d0*element(se)%tv_midnode(3,kn)*element(se)%v_midnode(:,3,kn)
              end do
              ! Top plane vertices
              do kn=1,3
                tmp_mesh%node_eid(global_sn+3+(kn-1))=global_sn+3+(kn-1)
                tmp_mesh%node_x(:,global_sn+3+(kn-1))=element(se)%x_gn(:,kn)+0.5d0*element(se)%tv_midnode(3,kn)*element(se)%v_midnode(:,3,kn)
              end do
              ! Bottom plane edge nodes
              kn = 4
              tmp_mesh%node_eid(global_sn+6)=global_sn+6
              tmp_mesh%node_x(:,global_sn+6)=element(se)%x_gn(:,kn)-0.5d0*element(se)%tv_midnode(3,kn)*element(se)%v_midnode(:,3,kn)
              kn = 5
              tmp_mesh%node_eid(global_sn+9)=global_sn+9
              tmp_mesh%node_x(:,global_sn+9)=element(se)%x_gn(:,kn)-0.5d0*element(se)%tv_midnode(3,kn)*element(se)%v_midnode(:,3,kn)
              kn = 6
              tmp_mesh%node_eid(global_sn+7)=global_sn+7
              tmp_mesh%node_x(:,global_sn+7)=element(se)%x_gn(:,kn)-0.5d0*element(se)%tv_midnode(3,kn)*element(se)%v_midnode(:,3,kn)
              ! Top plane edge nodes
              kn = 4
              tmp_mesh%node_eid(global_sn+12)=global_sn+12
              tmp_mesh%node_x(:,global_sn+12)=element(se)%x_gn(:,kn)+0.5d0*element(se)%tv_midnode(3,kn)*element(se)%v_midnode(:,3,kn)
              kn = 5
              tmp_mesh%node_eid(global_sn+14)=global_sn+14
              tmp_mesh%node_x(:,global_sn+14)=element(se)%x_gn(:,kn)+0.5d0*element(se)%tv_midnode(3,kn)*element(se)%v_midnode(:,3,kn)
              kn = 6
              tmp_mesh%node_eid(global_sn+13)=global_sn+13
              tmp_mesh%node_x(:,global_sn+13)=element(se)%x_gn(:,kn)+0.5d0*element(se)%tv_midnode(3,kn)*element(se)%v_midnode(:,3,kn)
              ! Mid-plane vertices
              kn = 1
              tmp_mesh%node_eid(global_sn+8)=global_sn+8
              tmp_mesh%node_x(:,global_sn+8)=element(se)%x_gn(:,kn)
              kn = 2
              tmp_mesh%node_eid(global_sn+10)=global_sn+10
              tmp_mesh%node_x(:,global_sn+10)=element(se)%x_gn(:,kn)
              kn = 3
              tmp_mesh%node_eid(global_sn+11)=global_sn+11
              tmp_mesh%node_x(:,global_sn+11)=element(se)%x_gn(:,kn)
              kn = 4
              tmp_mesh%node_eid(global_sn+15)=global_sn+15
              tmp_mesh%node_x(:,global_sn+15)=element(se)%x_gn(:,kn)
              kn = 5
              tmp_mesh%node_eid(global_sn+17)=global_sn+17
              tmp_mesh%node_x(:,global_sn+17)=element(se)%x_gn(:,kn)
              kn = 6
              tmp_mesh%node_eid(global_sn+16)=global_sn+16
              tmp_mesh%node_x(:,global_sn+16)=element(se)%x_gn(:,kn)
              ! Element construction
              tmp_mesh%element_eid     (global_se) = global_se
              tmp_mesh%element_type    (global_se) = 13
              tmp_mesh%element_physical(global_se) = part(sp)%id
              do kn=1,18
                tmp_mesh%element_node    (kn,global_se) = global_sn+(kn-1)
              end do
              global_sn = global_sn + 18
              global_se = global_se + 1
            end if
          end if
        end do
      end do
    end do

    ! Write
    tmp_filename=trim(output_filename)//'.shell3d.msh'
    call tmp_mesh%write(tmp_filename)

  end if
















  ! Close files
  do i=1,3
    close(tmp_unit(i))
  end do

  ! Ending message
  if (verbose_level.ge.1) write(output_unit,'(a)') 'done.'

end subroutine export_data_at_geometrical_nodes
