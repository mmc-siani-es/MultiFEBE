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
!! <b> This module implements mesh derived type and procedures (as a class) </b>
module fbem_mesh_module

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_data_structures

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all are private
  private

  ! Public
  public :: fbem_mesh
  public :: fbem_convert_mesh_file_format

  !! Mesh: parts > elements (and auxiliary subelements: subedges and subfaces) > nodes
  type fbem_mesh
    character(len=fbem_stdcharlen)  :: name
    integer                         :: n

    integer                         :: n_parts
    type(fbem_part), allocatable    :: part(:)
    integer                         :: part_eid_min
    integer                         :: part_eid_max
    integer, allocatable            :: part_iid(:)

    integer                         :: n_elements
    type(fbem_element), allocatable :: element(:)
    integer                         :: element_eid_min
    integer                         :: element_eid_max
    integer, allocatable            :: element_iid(:)

    integer                         :: n_nodes
    type(fbem_node), allocatable    :: node(:)
    integer                         :: node_eid_min
    integer                         :: node_eid_max
    integer, allocatable            :: node_iid(:)

    integer                         :: n_subedges
    type(fbem_element), allocatable :: subedge(:)

    integer                         :: n_subfaces
    type(fbem_element), allocatable :: subface(:)

  contains
    procedure, pass(mesh)           :: read_from_file
    procedure, pass(mesh)           :: write_to_file
    procedure, pass(mesh)           :: build_node_nodes_connectivity
    procedure, pass(mesh)           :: build_node_elements_connectivity
    procedure, pass(mesh)           :: read_node_values_from_file
    procedure, pass(mesh)           :: destroy
  end type fbem_mesh

contains

  subroutine destroy(mesh)
    implicit none
    class(fbem_mesh) :: mesh
    if (allocated(mesh%part       )) deallocate(mesh%part       )
    if (allocated(mesh%part_iid   )) deallocate(mesh%part_iid   )
    if (allocated(mesh%element    )) deallocate(mesh%element    )
    if (allocated(mesh%element_iid)) deallocate(mesh%element_iid)
    if (allocated(mesh%node       )) deallocate(mesh%node       )
    if (allocated(mesh%node_iid   )) deallocate(mesh%node_iid   )
  end subroutine destroy

  subroutine read_from_file(mesh,n,tol,filename,format)
    implicit none
    ! I/O
    class(fbem_mesh)                        :: mesh
    integer                                 :: n
    real(kind=real64)                       :: tol               !! Tolerance (radius of the ball)
    character(len=*)                        :: filename
    character(len=*)                        :: format
    ! Local
    integer                                 :: fileunit          ! Unit of the file to read from
    integer                                 :: iostat_var
    character(len=fbem_string_max_length)   :: iomsg_var
    integer                                 :: format_tag        ! Mode of reading: 1 (native format), 2 (gmsh format)
    character(len=fbem_stdcharlen)          :: section_name      ! Name of the section
    logical                                 :: found             ! Logical variable for sections and keywords
    integer                                 :: i, j              ! Counters
    integer                                 :: tmp_int           ! Temporary integer
    integer                                 :: tmp_n_tags        ! Temporary variable to read the number of tags of the element
    integer                                 :: tmp_tags(8)       ! Temporary variable to read the number of tags of the element
    character(len=fbem_stdcharlen)          :: tmp_type          ! Temporary variable to read the type of the element
    character(len=fbem_file_record_length)  :: linestr           ! Line
    character(len=fbem_file_record_length)  :: tmp_str              ! Temporary string
    integer                                 :: nwords, nc
    integer                                 :: tmp_id
    real(kind=real64)                       :: tmp_real
    ! Initialize mesh
    call mesh%destroy
    !
    ! AMBIENT SPACE OF THE MESH
    !
    if ((n.ge.1).and.(n.le.3)) then
      mesh%n=n
    else
      call fbem_error_message(error_unit,0,'n',n,'invalid value of n (R^n).')
    end if
    !
    ! OPEN FILE
    !
    fileunit=fbem_get_valid_unit()
    open(unit=fileunit,file=trim(filename),action='read',status='old',recl=fbem_file_record_length,iostat=iostat_var,iomsg=iomsg_var)
    if (iostat_var.ne.0) call fbem_error_message(error_unit,0,'iostat',abs(iostat_var),trim(iomsg_var))
    format_tag=0
    if (trim(format).eq.'multifebe') format_tag=1
    if (trim(format).eq.'gmsh'     ) format_tag=2
    if (format_tag.eq.0) then
      call fbem_error_message(error_unit,0,trim(format),0,'is an unknown format.')
    end if
    !
    ! Mesh format version
    !
    if (format_tag.eq.2) then
      section_name='MeshFormat'
      call fbem_search_section_gmsh(fileunit,section_name,found)
      if (found) then
        read(fileunit,'(a)') linestr
        nwords=fbem_count_words(linestr)
        if (nwords.eq.3) then
          tmp_str=fbem_extract_word(linestr,1)
          if (trim(tmp_str).eq.'2.2') then
            write(output_unit,*) 'This is a gmsh *.msh file version 2.2'
          else
            call fbem_error_message(error_unit,0,trim(filename),0,'multifebe can read only gmsh *.msh file version 2.2')
          end if
        else
          call fbem_error_message(error_unit,0,trim(filename),0,'This line in $MeshFormat must have 3 numbers.')
        end if
      else
        call fbem_error_message(error_unit,0,trim(filename),0,'The section $MeshFormat is required.')
      end if
    end if
    !
    ! PARTS
    !
    select case (format_tag)
      case (1)
        section_name='parts'
        call fbem_search_section(fileunit,section_name,found)
      case (2)
        section_name='PhysicalNames'
        call fbem_search_section_gmsh(fileunit,section_name,found)
    end select
    if (found) then
      ! READ THE NUMBER OF PARTS
      read(fileunit,*) mesh%n_parts
      if (mesh%n_parts.gt.0) then
        allocate (mesh%part(mesh%n_parts))
        do i=1,mesh%n_parts
          mesh%part(i)%id=0
          mesh%part(i)%name=''
          mesh%part(i)%type=0
          mesh%part(i)%entity=0
          mesh%part(i)%n_nodes=0
          mesh%part(i)%n_elements=0
          mesh%part(i)%export=.true.
        end do
      else
        call fbem_error_message(error_unit,0,trim(filename),mesh%n_parts,'the number of parts must be >0')
      end if
      ! READ EACH PART
      do i=1,mesh%n_parts
        select case (format_tag)
          case (1)
            read(fileunit,*) mesh%part(i)%id, mesh%part(i)%name
          case (2)
            read(fileunit,*) tmp_int, mesh%part(i)%id, mesh%part(i)%name
        end select
        if (mesh%part(i)%id.le.0) then
          call fbem_error_message(error_unit,0,'part',mesh%part(i)%id,'identifiers must be >0')
        end if
        call fbem_trim(mesh%part(i)%name)
      end do
      ! CHECK AND BUILD PARTS IDENTIFIERS
      mesh%part_eid_min=mesh%part(1)%id
      mesh%part_eid_max=mesh%part(1)%id
      do i=2,mesh%n_parts
        if (mesh%part(i)%id.lt.mesh%part_eid_min) mesh%part_eid_min=mesh%part(i)%id
        if (mesh%part(i)%id.gt.mesh%part_eid_max) mesh%part_eid_max=mesh%part(i)%id
      end do
      allocate (mesh%part_iid(mesh%part_eid_min:mesh%part_eid_max))
      mesh%part_iid=0
      do i=1,mesh%n_parts
        if (mesh%part_iid(mesh%part(i)%id).ne.0) then
          call fbem_error_message(error_unit,0,'part',mesh%part(i)%id,'is repeated')
        else
          mesh%part_iid(mesh%part(i)%id)=i
        end if
      end do
    else
      select case (format_tag)
        case (1)
          call fbem_error_message(error_unit,0,trim(filename),0,'The section [parts] is required.')
        case (2)
          call fbem_error_message(error_unit,0,trim(filename),0,'The section PhysicalNames is required.')
      end select
    end if
    !
    ! ELEMENTS
    !
    select case (format_tag)
      case (1)
        section_name='elements'
        call fbem_search_section(fileunit,section_name,found)
      case (2)
        section_name='Elements'
        call fbem_search_section_gmsh(fileunit,section_name,found)
    end select
    if (found) then
      ! READ THE NUMBER OF ELEMENTS
      read(fileunit,*) mesh%n_elements
      if (mesh%n_elements.gt.0) then
        allocate (mesh%element(mesh%n_elements))
        do i=1,mesh%n_elements
          mesh%element(i)%id=0
          mesh%element(i)%type=0
          mesh%element(i)%n_dimension=0
          mesh%element(i)%size=0.d0
          mesh%element(i)%part=0
          mesh%element(i)%n_nodes=0
          mesh%element(i)%element=0
          mesh%element(i)%dm_n_elements=0
          mesh%element(i)%dm_mode=0
          mesh%element(i)%plane=.false.
          mesh%element(i)%bball_radius=0.d0
          mesh%element(i)%n_subedges=0
          mesh%element(i)%n_subfaces=0
          mesh%element(i)%n_supelements=0
          mesh%element(i)%type_g=0
          mesh%element(i)%type_f1=0
          mesh%element(i)%type_f2=0
          mesh%element(i)%type_f1f=0
          mesh%element(i)%type_f2f=0
          mesh%element(i)%discontinuous=.false.
          mesh%element(i)%delta_f=0.d0
          mesh%element(i)%csize=0.d0
          mesh%element(i)%n_phi=0
          mesh%element(i)%fe_type=0
          mesh%element(i)%K_intmode=0
          mesh%element(i)%K_intngp=0
          mesh%element(i)%A=0.d0
          mesh%element(i)%I=0.d0
          mesh%element(i)%orthotropic_shell_fd1=0.d0
          mesh%element(i)%export=.true.
        end do
      else
        call fbem_error_message(error_unit,0,trim(filename),mesh%n_elements,'the number of elements must be >0')
      end if
      ! READ EACH ELEMENT
      do i=1,mesh%n_elements
        read(fileunit,*) mesh%element(i)%id, tmp_type, tmp_n_tags, (tmp_tags(j),j=1,tmp_n_tags)
        if (mesh%element(i)%id.le.0) then
          call fbem_error_message(error_unit,0,'element',mesh%element(i)%id,'identifiers must be >0')
        end if
        if ((tmp_tags(1).ge.mesh%part_eid_min).and.(tmp_tags(1).le.mesh%part_eid_max)) then
          if (mesh%part_iid(tmp_tags(1)).eq.0) then
            call fbem_warning_message(error_unit,0,'part',tmp_tags(1),'does not exist')
            call fbem_error_message(error_unit,0,'element',mesh%element(i)%id,'is assigned to the previous invalid part')
          else
            mesh%element(i)%part=mesh%part_iid(tmp_tags(1))
          end if
        else
          call fbem_warning_message(error_unit,0,'part',tmp_tags(1),'does not exist')
          call fbem_error_message(error_unit,0,'element',mesh%element(i)%id,'is assigned to the previous invalid part')
        end if
        if ((trim(tmp_type).eq.'point1').or.(trim(tmp_type).eq.'15')) mesh%element(i)%type=fbem_point1
        if ((trim(tmp_type).eq.'line2' ).or.(trim(tmp_type).eq.'1' )) mesh%element(i)%type=fbem_line2
        if ((trim(tmp_type).eq.'line3' ).or.(trim(tmp_type).eq.'8' )) mesh%element(i)%type=fbem_line3
        if ((trim(tmp_type).eq.'line4' ).or.(trim(tmp_type).eq.'26')) mesh%element(i)%type=fbem_line4
        if ((trim(tmp_type).eq.'tri3'  ).or.(trim(tmp_type).eq.'2' )) mesh%element(i)%type=fbem_tri3
        if ((trim(tmp_type).eq.'tri6'  ).or.(trim(tmp_type).eq.'9' )) mesh%element(i)%type=fbem_tri6
        if ((trim(tmp_type).eq.'quad4' ).or.(trim(tmp_type).eq.'3' )) mesh%element(i)%type=fbem_quad4
        if ((trim(tmp_type).eq.'quad8' ).or.(trim(tmp_type).eq.'16')) mesh%element(i)%type=fbem_quad8
        if ((trim(tmp_type).eq.'quad9' ).or.(trim(tmp_type).eq.'10')) mesh%element(i)%type=fbem_quad9
        if ((trim(tmp_type).eq.'tet4'  ).or.(trim(tmp_type).eq.'4' )) mesh%element(i)%type=fbem_tet4
        if ((trim(tmp_type).eq.'tet10' ).or.(trim(tmp_type).eq.'11')) mesh%element(i)%type=fbem_tet10
        if ((trim(tmp_type).eq.'hex8'  ).or.(trim(tmp_type).eq.'5' )) mesh%element(i)%type=fbem_hex8
        if ((trim(tmp_type).eq.'hex20' ).or.(trim(tmp_type).eq.'17')) mesh%element(i)%type=fbem_hex20
        if ((trim(tmp_type).eq.'hex27' ).or.(trim(tmp_type).eq.'12')) mesh%element(i)%type=fbem_hex27
        if (mesh%element(i)%type.eq.0) then
          call fbem_error_message(error_unit,0,'element',mesh%element(i)%id,'unrecognized type of element')
        end if
        mesh%element(i)%n_nodes=fbem_n_nodes(mesh%element(i)%type)
        mesh%element(i)%n_dimension=fbem_n_dimension(mesh%element(i)%type)
        ! Back to the line in order to read the nodes of the element
        backspace(fileunit)
        allocate(mesh%element(i)%node(mesh%element(i)%n_nodes))
        read(fileunit,*)tmp_id,tmp_type,tmp_n_tags,(tmp_tags(j),j=1,tmp_n_tags),(mesh%element(i)%node(j),j=1,mesh%element(i)%n_nodes)
      end do
      ! CHECK AND BUILD ELEMENT IDENTIFIERS
      mesh%element_eid_min=mesh%element(1)%id
      mesh%element_eid_max=mesh%element(1)%id
      do i=2,mesh%n_elements
        if (mesh%element(i)%id.lt.mesh%element_eid_min) mesh%element_eid_min=mesh%element(i)%id
        if (mesh%element(i)%id.gt.mesh%element_eid_max) mesh%element_eid_max=mesh%element(i)%id
      end do
      allocate (mesh%element_iid(mesh%element_eid_min:mesh%element_eid_max))
      mesh%element_iid=0
      do i=1,mesh%n_elements
        if (mesh%element_iid(mesh%element(i)%id).ne.0) then
          call fbem_error_message(error_unit,0,'element',mesh%element(i)%id,'is repeated')
        else
          mesh%element_iid(mesh%element(i)%id)=i
        end if
      end do
      ! BUILD CONNECTIVITY PART->ELEMENTS CONNECTIVITY
      do i=1,mesh%n_parts
        mesh%part(i)%n_elements=0
      end do
      do j=1,mesh%n_elements
        i=mesh%element(j)%part
        mesh%part(i)%n_elements=mesh%part(i)%n_elements+1
      end do
      do i=1,mesh%n_parts
        allocate (mesh%part(i)%element(mesh%part(i)%n_elements))
        mesh%part(i)%n_elements=0
      end do
      do j=1,mesh%n_elements
        i=mesh%element(j)%part
        mesh%part(i)%n_elements=mesh%part(i)%n_elements+1
        mesh%part(i)%element(mesh%part(i)%n_elements)=j
      end do
      ! CHECK THAT EACH PART CONTAINS ELEMENTS OF THE SAME DIMENSION
      do i=1,mesh%n_parts
        mesh%part(i)%n_dimension=mesh%element(mesh%part(i)%element(1))%n_dimension
        do j=2,mesh%part(i)%n_elements
          if (mesh%element(mesh%part(i)%element(j))%n_dimension.ne.mesh%part(i)%n_dimension) then
            call fbem_error_message(error_unit,0,'part',mesh%part(i)%id,'contain elements of different dimensions')
          end if
        end do
      end do
    else
      select case (format_tag)
        case (1)
          call fbem_error_message(error_unit,0,trim(filename),0,'The section [elements] is required.')
        case (2)
          call fbem_error_message(error_unit,0,trim(filename),0,'The section Elements is required.')
      end select
    end if
    !
    ! READ NODES
    !
    select case (format_tag)
      case (1)
        section_name='nodes'
        call fbem_search_section(fileunit,section_name,found)
      case (2)
        section_name='Nodes'
        call fbem_search_section_gmsh(fileunit,section_name,found)
    end select
    if (found) then
      ! READ THE NUMBER OF NODES
      read(fileunit,*) mesh%n_nodes
      if (mesh%n_nodes.gt.0) then
        allocate (mesh%node(mesh%n_nodes))
        do i=1,mesh%n_nodes
          allocate(mesh%node(i)%x(mesh%n))
          mesh%node(i)%x=0.d0
          mesh%node(i)%export=.true.
        end do
      else
        call fbem_error_message(error_unit,0,trim(filename),mesh%n_nodes,'the number of nodes must be >0')
      end if
      ! READ EACH NODE
      do i=1,mesh%n_nodes
        read(fileunit,'(a)') linestr
        nwords=fbem_count_words(linestr)
        if (nwords.gt.4) then
          call fbem_error_message(error_unit,0,trim(filename),nwords,'invalid number of arguments in the list of nodes')
        end if
        nwords=min0(nwords-1,mesh%n)+1
        read(linestr,*) mesh%node(i)%id
        do j=2,nwords
          tmp_str=fbem_extract_word(linestr,j)
          call fbem_trimall(tmp_str)
          nc=scan(tmp_str,'eEdDqQ')
          if (nc.eq.0) then
            tmp_str(len_trim(tmp_str)+1:len_trim(tmp_str)+2)='D0'
          else
            nc=scan(tmp_str,'eE')
            if (nc.ne.0) tmp_str(nc:nc)='D'
          end if
          read(tmp_str,*) mesh%node(i)%x(j-1)
        end do
        if (mesh%node(i)%id.le.0) then
          call fbem_error_message(error_unit,0,'node',mesh%node(i)%id,'identifiers must be >0')
        end if
      end do
      ! CHECK AND BUILD NODES IDENTIFIERS
      mesh%node_eid_min=mesh%node(1)%id
      mesh%node_eid_max=mesh%node(1)%id
      do i=2,mesh%n_nodes
        if (mesh%node(i)%id.lt.mesh%node_eid_min) mesh%node_eid_min=mesh%node(i)%id
        if (mesh%node(i)%id.gt.mesh%node_eid_max) mesh%node_eid_max=mesh%node(i)%id
      end do
      allocate (mesh%node_iid(mesh%node_eid_min:mesh%node_eid_max))
      mesh%node_iid=0
      do i=1,mesh%n_nodes
        if (mesh%node_iid(mesh%node(i)%id).ne.0) then
          call fbem_error_message(error_unit,0,'node',mesh%node(i)%id,'is repeated')
        else
          mesh%node_iid(mesh%node(i)%id)=i
        end if
      end do
      ! CHANGE EID TO IID IN THE ELEMENT->NODES CONNECTIVITY
      do i=1,mesh%n_elements
        do j=1,mesh%element(i)%n_nodes
          if ((mesh%element(i)%node(j).ge.mesh%node_eid_min).and.(mesh%element(i)%node(j).le.mesh%node_eid_max)) then
            if (mesh%node_iid(mesh%element(i)%node(j)).eq.0) then
              call fbem_warning_message(error_unit,0,'node',mesh%element(i)%node(j),'does not exist')
              call fbem_error_message(error_unit,0,'element',mesh%element(i)%id,'contain the previous node')
            else
              mesh%element(i)%node(j)=mesh%node_iid(mesh%element(i)%node(j))
            end if
          else
            call fbem_warning_message(error_unit,0,'node',mesh%element(i)%node(j),'does not exist')
            call fbem_error_message(error_unit,0,'element',mesh%element(i)%id,'contain the previous node')
          end if
        end do
      end do
    else
      select case (format_tag)
        case (1)
          call fbem_error_message(error_unit,0,trim(filename),0,'The section [nodes] is required.')
        case (2)
          call fbem_error_message(error_unit,0,trim(filename),0,'The section Nodes is required.')
      end select
    end if
    !
    ! CLOSE FILE
    !
    close(unit=fileunit)
    !
    ! NODE->NODES CONNECTIVITY
    !
    call mesh%build_node_nodes_connectivity(tol)
    call mesh%build_node_elements_connectivity
    !
    ! OTHER CONNECTIVITIES
    !
!    call fbem_node_parts_connectivity(n_nodes,node,n_elements,element,n_parts,part)
!    call fbem_part_nodes_connectivity(n_nodes,node,n_parts,part)
!    call fbem_build_mesh_subelements(n_nodes,node,n_elements,element,n_subedges,subedge,n_subfaces,subface)
  end subroutine read_from_file

  subroutine write_to_file(mesh,filename,format)
    implicit none
    ! I/O
    class(fbem_mesh)                        :: mesh
    character(len=*)                        :: filename
    character(len=*)                        :: format
    ! Local
    integer                                 :: fileunit          ! Unit of the file to read from
    integer                                 :: iostat_var
    character(len=fbem_string_max_length)   :: iomsg_var
    integer                                 :: format_tag        ! Mode of reading: 1 (multifebe format), 2 (gmsh format)
    character(len=fbem_stdcharlen)          :: section_name      ! Name of the section
    logical                                 :: found             ! Logical variable for sections and keywords
    integer                                 :: i, j              ! Counters
    integer                                 :: tmp_int           ! Temporary integer
    integer                                 :: tmp_n_tags        ! Temporary variable to read the number of tags of the element
    integer                                 :: tmp_tags(8)       ! Temporary variable to read the number of tags of the element
    character(len=fbem_stdcharlen)          :: tmp_type          ! Temporary variable to read the type of the element
    character(len=fbem_file_record_length)  :: linestr           ! Line
    integer                                 :: nwords
    integer                                 :: tmp_id
    character(len=fbem_fmtstr)              :: fmtstr
    !
    ! CHECK FORMAT
    !
    format_tag=0
    if (trim(format).eq.'multifebe') format_tag=1
    if (trim(format).eq.'gmsh'     ) format_tag=2
    if (format_tag.eq.0) then
      call fbem_error_message(error_unit,0,trim(format),0,'is an unknown format.')
    end if
    !
    ! OPEN FILE
    !
    fileunit=fbem_get_valid_unit()
    ! status='new' ???
    open(unit=fileunit,file=trim(filename),action='write',recl=fbem_file_record_length,iostat=iostat_var,iomsg=iomsg_var)
    if (iostat_var.ne.0) call fbem_error_message(error_unit,0,'iostat',iostat_var,trim(iomsg_var))
    !
    ! HEADER
    !
    select case (format_tag)
      case (1)
        stop 'not yet : write_to_file'
      case (2)
        write(fileunit,'(a)') '$MeshFormat'
        write(fileunit,'(a)') '2.2 0 8'
        write(fileunit,'(a)') '$EndMeshFormat'
    end select
    !
    ! PARTS
    !
    select case (format_tag)
      case (1)
        stop 'not yet : write_to_file'
      case (2)
        write(fileunit,'(a)') '$PhysicalNames'
        write(fmtstr,*) '(i',fbem_nchar_int(mesh%n_parts),')'
        call fbem_trimall(fmtstr)
        write(fileunit,fmtstr) mesh%n_parts
        do i=1,mesh%n_parts
          write(fmtstr,*) '(i1,1x,i',fbem_nchar_int(mesh%part(i)%id),',1x,a1,a',len_trim(mesh%part(i)%name),',a1)'
          call fbem_trimall(fmtstr)
          write(fileunit,fmtstr) mesh%part(i)%n_dimension, mesh%part(i)%id, '"',trim(mesh%part(i)%name),'"'
        end do
        write(fileunit,'(a)') '$EndPhysicalNames'
    end select
    !
    ! NODES
    !
    select case (format_tag)
      case (1)
        stop 'not yet : write_to_file'
      case (2)
        write(fileunit,'(a)') '$Nodes'
        write(fmtstr,*) '(i',fbem_nchar_int(mesh%n_nodes),')'
        call fbem_trimall(fmtstr)
        write(fileunit,fmtstr) mesh%n_nodes
        select case (mesh%n)
          case (2)
            do i=1,mesh%n_nodes
              write(fmtstr,*) '(i',fbem_nchar_int(mesh%node(i)%id),',3e25.16)'
              call fbem_trimall(fmtstr)
              write(fileunit,fmtstr) mesh%node(i)%id, mesh%node(i)%x(1:2), 0.
            end do
          case (3)
            do i=1,mesh%n_nodes
              write(fmtstr,*) '(i',fbem_nchar_int(mesh%node(i)%id),',3e25.16)'
              call fbem_trimall(fmtstr)
              write(fileunit,fmtstr) mesh%node(i)%id, mesh%node(i)%x(1:3)
            end do
          case default
            stop 'mesh%n not valid'
        end select
        write(fileunit,'(a)') '$EndNodes'
    end select
    !
    ! ELEMENTS
    !
    select case (format_tag)
      case (1)
        stop 'not yet : write_to_file'
      case (2)
        write(fileunit,'(a)') '$Elements'
        write(fmtstr,*) '(i',fbem_nchar_int(mesh%n_elements),')'
        call fbem_trimall(fmtstr)
        write(fileunit,fmtstr) mesh%n_elements
        do i=1,mesh%n_elements
          write(fmtstr,*) '(i',fbem_nchar_int(mesh%element(i)%id),')'
          call fbem_trimall(fmtstr)
          write(fileunit,fmtstr,advance='no') mesh%element(i)%id
          select case (mesh%element(i)%type)
            case (fbem_point1)
              tmp_int=15
            case (fbem_line2)
              tmp_int=1
            case (fbem_line3)
              tmp_int=8
            case (fbem_line4)
              tmp_int=26
            case (fbem_tri3)
              tmp_int=2
            case (fbem_tri6)
              tmp_int=9
            case (fbem_quad4)
              tmp_int=3
            case (fbem_quad8)
              tmp_int=16
            case (fbem_quad9)
              tmp_int=10
            case (fbem_tet4)
              tmp_int=4
            case (fbem_tet10)
              tmp_int=11
            case (fbem_hex8)
              tmp_int=5
            case (fbem_hex20)
              tmp_int=17
            case (fbem_hex27)
              tmp_int=12
            case default
              stop 'unknown element()%type'
          end select
          write(fmtstr,*) '(1x,i',fbem_nchar_int(tmp_int),')'
          call fbem_trimall(fmtstr)
          write(fileunit,fmtstr,advance='no') tmp_int
          tmp_int=mesh%part(mesh%element(i)%part)%id
          write(fmtstr,*) '(1x,i1,1x,i',fbem_nchar_int(tmp_int),',1x,i',fbem_nchar_int(tmp_int),')'
          call fbem_trimall(fmtstr)
          write(fileunit,fmtstr,advance='no') 2, tmp_int, tmp_int
          tmp_int=0
          do j=1,mesh%element(i)%n_nodes
            if (tmp_int.lt.mesh%node(mesh%element(i)%node(j))%id) then
              tmp_int=mesh%node(mesh%element(i)%node(j))%id
            end if
          end do
          write(fmtstr,*) '(',mesh%element(i)%n_nodes,'i',fbem_nchar_int(tmp_int)+1,')'
          call fbem_trimall(fmtstr)
          write(fileunit,fmtstr) (mesh%node(mesh%element(i)%node(j))%id,j=1,mesh%element(i)%n_nodes)
        end do
        write(fileunit,'(a)') '$EndElements'
    end select
    close(unit=fileunit)
  end subroutine write_to_file

  !! Build node->nodes connectivity
  subroutine build_node_nodes_connectivity(mesh,tol)
    implicit none
    ! I/O
    class(fbem_mesh)     :: mesh           !! Mesh
    real(kind=real64)    :: tol            !! Tolerance (radius of the ball)
    ! Local
    integer              :: kni, knj       ! Counters
    integer              :: snj            ! Identifier
    real(kind=real64)    :: rv(mesh%n)     ! Distance vector
    real(kind=real64)    :: r              ! Distance
    integer, allocatable :: node_nodes(:)  ! node -> nodes connectivity by distance
    !
    ! Initialise
    !
    allocate (node_nodes(mesh%n_nodes))
    do kni=1,mesh%n_nodes
      mesh%node(kni)%n_nodes=0
    end do
    !
    ! Build the connectivity vector of nodes for each node (distance symmetry property is fully exploited)
    !
    do kni=1,mesh%n_nodes
      if (mesh%node(kni)%n_nodes.eq.0) then
        node_nodes=0
        do knj=kni+1,mesh%n_nodes
          rv=mesh%node(kni)%x-mesh%node(knj)%x
          r=sqrt(dot_product(rv,rv))
          if (r.le.tol) then
            mesh%node(kni)%n_nodes=mesh%node(kni)%n_nodes+1
            node_nodes(mesh%node(kni)%n_nodes)=knj
          end if
        end do
        if (mesh%node(kni)%n_nodes.gt.0) then
          ! Allocate and copy the connectivity vector node i -> node j to the data structure
          allocate (mesh%node(kni)%node(mesh%node(kni)%n_nodes))
          mesh%node(kni)%node=node_nodes(1:mesh%node(kni)%n_nodes)
          ! Copy the connectivity node j -> node i
          do knj=1,mesh%node(kni)%n_nodes
            snj=node_nodes(knj)
            if (mesh%node(snj)%n_nodes.eq.0) then
              mesh%node(snj)%n_nodes=mesh%node(kni)%n_nodes
              allocate (mesh%node(snj)%node(mesh%node(snj)%n_nodes))
              mesh%node(snj)%node=node_nodes(1:mesh%node(snj)%n_nodes)
              ! Replace the appearance of node j in node_nodes by node i
              mesh%node(snj)%node(knj)=kni
            else
              call fbem_error_message(error_unit,0,'node',mesh%node(kni)%id,&
              'the set of nodes such that dist(x(node)-x)<=gtol must be the same for all nodes in the set: check geometric tolerance and mesh')
            end if
          end do
        end if
      end if
    end do
    !
    ! Finalise
    !
    deallocate (node_nodes)
!    ! Write
!    allocate (node_node_i(mesh%n_nodes))
!    do kni=1,mesh%n_nodes
!      node_node_i(1:mesh%node(kni)%n_nodes)=mesh%node(kni)%node
!      call fbem_quicksort(1,mesh%node(kni)%n_nodes,mesh%node(kni)%n_nodes,node_node_i(1:mesh%node(kni)%n_nodes))
!      write(22,'(i11)') kni
!      do knj=1,mesh%node(kni)%n_nodes
!        write(22,'(i11)') node_node_i(knj)
!      end do
!    end do
!    pause
  end subroutine build_node_nodes_connectivity

  !! Build node->elements connectivity
  subroutine build_node_elements_connectivity(mesh)
    implicit none
    ! I/O
    class(fbem_mesh)   :: mesh           !! Mesh
    ! Local
    integer            :: kni, ke, se, knj
    integer            :: nd                  ! Number of dimensions of the element
    integer            :: n_nd_elements(3)    ! Number of 1d, 2d and 3d elements connected to a node
    ! Loop through nodes
    do kni=1,mesh%n_nodes
      !
      ! Calculate the number of elements that are connected with the node
      !
      ! Initialize
      mesh%node(kni)%n_elements=0
      ! Loop through elements
      do ke=1,mesh%n_elements
        ! Loop through the nodes of the element
        do knj=1,mesh%element(ke)%n_nodes
          if (mesh%element(ke)%node(knj).eq.kni) mesh%node(kni)%n_elements=mesh%node(kni)%n_elements+1
        end do
      end do
      !
      ! An unconnected node is not allowed
      !
      if (mesh%node(kni)%n_elements.eq.0) then
        call fbem_error_message(error_unit,0,'node',mesh%node(kni)%id,'this node is not connected any element')
      end if
      !
      ! Build the connectivity vectors
      !
      ! Allocate the connectivity vectors
      allocate (mesh%node(kni)%element(mesh%node(kni)%n_elements))
      allocate (mesh%node(kni)%element_node_iid(mesh%node(kni)%n_elements))
      allocate (mesh%node(kni)%element_node_loctype(mesh%node(kni)%n_elements))
      ! Initialize
      mesh%node(kni)%n_elements=0
      ! Loop through elements
      do ke=1,mesh%n_elements
        do knj=1,mesh%element(ke)%n_nodes
          if (mesh%element(ke)%node(knj).eq.kni) then
            mesh%node(kni)%n_elements=mesh%node(kni)%n_elements+1
            mesh%node(kni)%element(mesh%node(kni)%n_elements)=ke
            mesh%node(kni)%element_node_iid(mesh%node(kni)%n_elements)=knj
            if (mesh%element(ke)%type.ne.fbem_point1) then
              mesh%node(kni)%element_node_loctype(mesh%node(kni)%n_elements)=fbem_node_loctype(knj,mesh%element(ke)%type)
            else
              mesh%node(kni)%element_node_loctype(mesh%node(kni)%n_elements)=fbem_loctype_vertex
            end if
          end if
        end do
      end do
      !
      ! Determine what classes (line, surface, volume) of elements are connected to the node
      !
      ! Count the number of 1D (line), 2D (surface) and 3D (volume) elements connected to a node
      n_nd_elements=0
      do ke=1,mesh%node(kni)%n_elements
        se=mesh%node(kni)%element(ke)
        nd=mesh%element(se)%n_dimension
        if (mesh%element(se)%type.ne.fbem_point1) n_nd_elements(nd)=n_nd_elements(nd)+1
      end do
      ! Possible dimensional situations of the node connected with its elements:
      ! F F F = 0: Unconnected node.
      ! T F F = 1: Only connected to 1D elements.
      ! F T F = 2: Only connected to 2D elements.
      ! T T F = 3: Connected to 1D and 2D elements.
      ! F F T = 4: Only connected to 3D elements.
      ! T F T = 5: Connected to 1D and 3D elements.
      ! F T T = 6: Connected to 2D and 3D elements.
      ! T T T = 7: Connected to 1D, 2D and 3D elements.
      mesh%node(kni)%dimensional_degree=0
      do nd=1,3
        if (n_nd_elements(nd).gt.0) mesh%node(kni)%dimensional_degree=mesh%node(kni)%dimensional_degree+2**(nd-1)
      end do
      ! We interchange the code 3 and 4, to have more coherent codes.
      if (mesh%node(kni)%dimensional_degree.eq.3) then
        mesh%node(kni)%dimensional_degree=4
      else
        if (mesh%node(kni)%dimensional_degree.eq.4) mesh%node(kni)%dimensional_degree=3
      end if
      ! Definitive code for dimensional_degree:
      ! 0: Unconnected node
      ! 1: Only connected to 1D elements.
      ! 2: Only connected to 2D elements.
      ! 3: Only connected to 3D elements.
      ! 4: Connected to 1D and 2D elements, i.e. the node is connected with >=2 elements which are 1D and 2D elements.
      ! 5: Connected to 1D and 3D elements, i.e. the node is connected with >=2 elements which are 1D and 3D elements.
      ! 6: Connected to 2D and 3D elements, i.e. the node is connected with >=2 elements which are 2D and 3D elements.
      ! 7: Connected to 1D, 2D and 3D elements, i.e. the node is connected with >=3 elements which are 1D, 2D and 3D elements.
    end do
  end subroutine build_node_elements_connectivity

  !! Read node values from file
  subroutine read_node_values_from_file(mesh,filename)
    implicit none
    ! I/O
    class(fbem_mesh)                        :: mesh
    character(len=*)                        :: filename
    ! Local
    integer                                 :: fileunit          ! Unit of the file to read from
    integer                                 :: i, j, k, l        ! Counters
    character(len=fbem_stdcharlen)          :: tmp_type          ! Temporary variable to read the type of the element
    integer                                 :: tmp_type_tag
    integer                                 :: tmp_id, kdof_min, kdof_max, kgroup_min, kgroup_max
    real(kind=real64), allocatable          :: value_r(:)
    complex(kind=real64), allocatable       :: value_c(:)
    integer                                 :: n
    integer                                 :: n_values
    integer                                 :: sn
    ! OPEN FILE
    fileunit=fbem_get_valid_unit()
    open(unit=fileunit,file=trim(filename),action='read',recl=fbem_file_record_length)
    ! READ
    ! <number of nodes>
    read(fileunit,*) n
    if ((n.le.0).or.(n.gt.mesh%n_nodes)) then
      call fbem_error_message(error_unit,0,trim(filename),n,'n does not meet 0<n<=n_nodes')
    end if
    do i=1,n
      ! <id> <type: real or complex> <kdof_min> <kdof_max> <kgroup_min> <kgroup_max> <value(kdof_min,kgroup_min)> ... <value(kdof_max,kgroup_min)> ... <value(kdof_min,kgroup_max)> <value(kdof_max,kgroup_max)>
      read(fileunit,*) tmp_id, tmp_type, kdof_min, kdof_max, kgroup_min, kgroup_max
      if ((tmp_id.ge.mesh%node_eid_min).and.(tmp_id.le.mesh%node_eid_max)) then
        if (mesh%node_iid(tmp_id).eq.0) then
          call fbem_error_message(error_unit,0,trim(filename),tmp_id,'this node does not exist in the mesh.')
        else
          sn=mesh%node_iid(tmp_id)
        end if
      else
        call fbem_error_message(error_unit,0,trim(filename),tmp_id,'this node does not exist in the mesh.')
      end if
      tmp_type_tag=0
      if (trim(tmp_type).eq.'real'   ) tmp_type_tag=1
      if (trim(tmp_type).eq.'complex') tmp_type_tag=2
      if (tmp_type_tag.eq.0) then
        call fbem_error_message(error_unit,0,trim(filename),tmp_id,'the node value data type is not valid.')
      end if
      if ((kdof_max.lt.kdof_min).or.(kgroup_max.lt.kgroup_min)) then
        call fbem_error_message(error_unit,0,trim(filename),0,'kdof_max<kdof_min or kgroup_max<kgroup_min.')
      end if
      n_values=(kdof_max-kdof_min+1)*(kgroup_max-kgroup_min+1)
      backspace(fileunit)
      select case (tmp_type_tag)
        case (1)
          allocate(value_r(n_values))
          read(fileunit,*) tmp_id, tmp_type, kdof_min, kdof_max, kgroup_min, kgroup_max, (value_r(j),j=1,n_values)
          if (allocated(mesh%node(sn)%value_r)) then
            call fbem_error_message(error_unit,0,trim(filename),mesh%node(sn)%id,'this node is repeated.')
          end if
          allocate(mesh%node(sn)%value_r(kdof_min:kdof_max,kgroup_min:kgroup_max))
          l=1
          do j=kgroup_min,kgroup_max
            do k=kdof_min,kdof_max
              mesh%node(sn)%value_r(k,j)=value_r(l)
              l=l+1
            end do
          end do
          deallocate(value_r)
        case (2)
          allocate(value_c(n_values))
          read(fileunit,*) tmp_id, tmp_type, kdof_min, kdof_max, kgroup_min, kgroup_max, (value_c(j),j=1,n_values)
          if (allocated(mesh%node(sn)%value_r)) then
            call fbem_error_message(error_unit,0,trim(filename),mesh%node(sn)%id,'this node is repeated.')
          end if
          allocate(mesh%node(sn)%value_c(kdof_min:kdof_max,kgroup_min:kgroup_max))
          l=1
          do j=kgroup_min,kgroup_max
            do k=kdof_min,kdof_max
              mesh%node(sn)%value_c(k,j)=value_c(l)
              l=l+1
            end do
          end do
          deallocate(value_c)
      end select
    end do
    ! CLOSE FILE
    close(unit=fileunit)
  end subroutine read_node_values_from_file

  ! Read a mesh from a file in a given format, and write it in another given format
  subroutine fbem_convert_mesh_file_format(n,infilename,informat,outfilename,outformat)
    implicit none
    ! I/O
    integer          :: n
    character(len=*) :: infilename
    character(len=*) :: informat
    character(len=*) :: outfilename
    character(len=*) :: outformat
    ! Local
    type(fbem_mesh)  :: mesh
    ! Read and write mesh
    write(output_unit,*) 'Reading mesh from "', trim(infilename), '"'
    call mesh%read_from_file(n,-1.d0,trim(infilename),trim(informat))
    write(output_unit,*) 'Done'
    write(output_unit,*) 'Writing mesh to "', trim(outfilename), '"'
    call mesh%write_to_file(trim(outfilename),trim(outformat))
    write(output_unit,*) 'Done'
    call mesh%destroy
  end subroutine fbem_convert_mesh_file_format

end module fbem_mesh_module
