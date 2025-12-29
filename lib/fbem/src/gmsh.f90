! ---------------------------------------------------------------------
! Copyright (C) 2014-2025 Universidad de Las Palmas de Gran Canaria:
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
!! <b> This module implements data types and routines for handle Gmsh files. </b>

! TODO: leer y escribir version posterior de malla msh. Falta rutinas de escritura.
! LA IDEA ES PONER AQUI LECTURA Y ESCRITURA DE FICHEROS A TRAVES DE ESTRUCTURAS DE DATOS ESPECIFICAS
! DE GMSH, Y LUEGO EN OTRO MODULO (POR EJEMPLO EL MODULO MESH) LAS TRANSFORMACIONES ENTRE ESTRUCTURAS DE DATOS.


module fbem_gmsh

  ! Fortran 2003 intrinsic module
  use iso_fortran_env
  ! fbem modules
  use fbem_string_handling

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all is private
  private

  ! Public
  public :: fbem_gmsh_mesh
  public :: fbem_gmsh_search_section
  public :: fbem_export_gmsh_fmt_real
  public :: fbem_export_gmsh_NodeData
  public :: fbem_export_gmsh_ElementData
  public :: fbem_export_gmsh_ElementNodeData

  ! Element types in Gmsh:
  !
  ! 1   2-node line.
  ! 2   3-node triangle.
  ! 3   4-node quadrangle.
  ! 4   4-node tetrahedron.
  ! 5   8-node hexahedron.
  ! 6   6-node prism.
  ! 7   5-node pyramid.
  ! 8   3-node second order line (2 nodes associated with the vertices and 1 with the edge).
  ! 9   6-node second order triangle (3 nodes associated with the vertices and 3 with the edges).
  ! 10  9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face).
  ! 11  10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges).
  ! 12  27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume).
  ! 13  18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces).
  ! 14  14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face).
  ! 15  1-node point.
  ! 16  8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges).
  ! 17  20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges).
  !
  ! Limited to second order
  !
  integer, dimension(17), parameter :: fbem_gmsh_type_n_nodes=[2,3,4,4,8,6,5,3,6,9,10,27,18,14,1,8,20]

  ! Exporting precision
  character(len=fbem_fmtstr) :: fmt_real='e16.8e2'

  ! MSH 2.2 FILE FORMAT: https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-2-_0028Legacy_0029
  ! TODO: MSH 4.1 FILE FORMAT: https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format

  type fbem_gmsh_mesh
    ! $MeshFormat (valid 2.2)
    character(len=3)                :: version  ! Only 2.2
    integer                         :: filetype ! Only ASCII (filetype=0) is valid.
    integer                         :: datasize ! Only 8 is valid
    ! $PhysicalNames (valid 2.2)
    integer                         :: n_physicalnames
    integer, allocatable            :: physicalname_dim(:)
    integer, allocatable            :: physicalname_eid(:)
    integer                         :: physicalname_eid_min
    integer                         :: physicalname_eid_max
    integer, allocatable            :: physicalname_iid(:)
    character(len=127), allocatable :: physicalname_name(:)
    ! $Nodes
    integer                         :: n_nodes
    integer, allocatable            :: node_eid(:)
    integer                         :: node_eid_min
    integer                         :: node_eid_max
    integer, allocatable            :: node_iid(:)
    real(kind=real64), allocatable  :: node_x(:,:) ! (index1,index2) : index1 is 1==x, 2==y, 3==z; index2 is point id.
    ! $Elements
    integer                         :: n_elements
    integer, allocatable            :: element_eid(:)
    integer                         :: element_eid_min
    integer                         :: element_eid_max
    integer, allocatable            :: element_iid(:)
    integer, allocatable            :: element_type(:)
    integer, allocatable            :: element_physical(:) ! eid of the physical
    integer, allocatable            :: element_node(:,:)   ! eid of nodes
    ! $Periodic (not used, not available)
  contains
    procedure, pass(mesh)           :: init
    procedure, pass(mesh)           :: read
    !procedure, pass(mesh)           :: check
    !procedure, pass(mesh)           :: write
  end type fbem_gmsh_mesh

contains

  subroutine init(mesh)
    implicit none
    class(fbem_gmsh_mesh) :: mesh
    mesh%version=''
    mesh%filetype=0
    mesh%datasize=0
    mesh%n_physicalnames=0
    mesh%n_nodes=0
    mesh%n_elements=0
    if (allocated(mesh%physicalname_dim  )) deallocate(mesh%physicalname_dim  )
    if (allocated(mesh%physicalname_eid  )) deallocate(mesh%physicalname_eid  )
    if (allocated(mesh%physicalname_iid  )) deallocate(mesh%physicalname_iid  )
    if (allocated(mesh%physicalname_name )) deallocate(mesh%physicalname_name )
    if (allocated(mesh%node_eid          )) deallocate(mesh%node_eid          )
    if (allocated(mesh%node_iid          )) deallocate(mesh%node_iid          )
    if (allocated(mesh%node_x            )) deallocate(mesh%node_x            )
    if (allocated(mesh%element_eid       )) deallocate(mesh%element_eid       )
    if (allocated(mesh%element_iid       )) deallocate(mesh%element_iid       )
    if (allocated(mesh%element_type      )) deallocate(mesh%element_type      )
    if (allocated(mesh%element_physical  )) deallocate(mesh%element_physical  )
    if (allocated(mesh%element_node      )) deallocate(mesh%element_node      )
  end subroutine init

  subroutine set_eid_to_iid(entity_str,n_eid,eid,eid_min,eid_max,iid)
    implicit none
    character(len=32)                 :: entity_str
    integer, intent(in)               :: n_eid
    integer, intent(in)               :: eid(n_eid)
    integer, intent(out)              :: eid_min
    integer, intent(out)              :: eid_max
    integer, intent(out), allocatable :: iid(:)
    integer :: i
    eid_min=minval(eid)
    eid_max=maxval(eid)
    if (allocated(iid)) deallocate(iid)
    allocate (iid(eid_min:eid_max))
    iid=0
    do i=1,n_eid
      if (iid(eid(i)).ne.0) then
        call fbem_error_message(error_unit,0,trim(entity_str),eid(i),'is repeated.')
      else
        iid(eid(i))=i
      end if
    end do
  end subroutine set_eid_to_iid

  !! It finds the first line in a file that matches "$section". The file position remains in the next line.
  subroutine fbem_gmsh_search_section(selected_unit,section,found)
    implicit none
    ! I/O
    integer                                :: selected_unit !! Unit of the file
    character(len=*)                       :: section       !! Section name to be found. It can't contain blanks.
    logical                                :: found         !! True if the section has been found, false otherwise.
    ! Local
    character(len=fbem_file_record_length) :: line          ! Line
    integer                                :: word_length   ! Word length
    integer                                :: file_line     ! Current line, at the end, it takes the line where the section line is located.
    integer                                :: ios           ! Error flag
    ! Rewind to the begin of file
    rewind(selected_unit)
    file_line=0
    ! Initialize
    found=.false.
    ! Reading process
    do
      ! Read line
      read(selected_unit,'(a)',iostat=ios) line
      if (is_iostat_end(ios)) then
        exit
      end if
      ! If at the beginning of the file, check if BOM is present.
      if (file_line.eq.0) then
        if ((ichar(line(1:1)).eq.239).and.(ichar(line(2:2)).eq.187).and.(ichar(line(3:3)).eq.191)) line(1:3)=''
      end if
      file_line=file_line+1
      ! Clean the line from spaces at the start and at the end
      call fbem_trim(line)
      ! Check if the line starts with '$'
      if (line(1:1).eq.'$') then
        word_length=len_trim(line)
        ! Pick up the string after '$'
        line=line(2:word_length)
        ! Check if the section name is the same
        if (trim(line).eq.trim(section)) then
          found=.true.
          exit
        end if
      end if
    end do
  end subroutine fbem_gmsh_search_section

  subroutine read(mesh,filename)
    implicit none
    ! I/O
    class(fbem_gmsh_mesh), intent(out)      :: mesh
    character(len=*), intent(in)            :: filename
    ! Local
    integer                                 :: fileunit          ! Unit of the file to read from
    integer                                 :: iostat_var
    character(len=fbem_string_max_length)   :: iomsg_var
    character(len=fbem_stdcharlen)          :: section_name      ! Name of the section
    logical                                 :: found             ! Logical variable for sections and keywords
    integer                                 :: i, j              ! Counters
    integer                                 :: tmp_int           ! Temporary integer
    integer                                 :: tmp_n_tags        ! Temporary variable to read the number of tags of the element
    integer                                 :: tmp_tags(8)       ! Temporary variable to read the number of tags of the element
    character(len=fbem_stdcharlen)          :: tmp_type          ! Temporary variable to read the type of the element
    character(len=fbem_file_record_length)  :: linestr           ! Line
    character(len=fbem_file_record_length)  :: tmp_str           ! Temporary string
    integer                                 :: nwords, nc
    integer                                 :: tmp_id
    real(kind=real64)                       :: tmp_real
    character(len=32)                       :: entity_str
    !
    ! Initialize mesh
    !
    call mesh%init
    !
    ! Open file
    !
    call fbem_open_file_to_read(filename,'This is a Gmsh file.',fileunit)
    !
    ! Read $MeshFormat
    !
    section_name='MeshFormat'
    call fbem_gmsh_search_section(fileunit,section_name,found)
    if (found) then
      read(fileunit,'(a)') linestr
      nwords=fbem_count_words(linestr)
      if (nwords.eq.3) then
        tmp_str=fbem_extract_word(linestr,1)
        if (trim(tmp_str).eq.'2.2') then
          mesh%version='2.2'
          call fbem_timestamp_w_message(output_unit,2,'This is a gmsh *.msh file version 2.2')
        else
          call fbem_error_message(error_unit,0,trim(filename),0,'multifebe can read only gmsh *.msh file version 2.2.')
        end if
      else
        call fbem_error_message(error_unit,0,trim(filename),0,'This line in $MeshFormat must have 3 arguments.')
      end if
    else
      call fbem_error_message(error_unit,0,trim(filename),0,'The section $MeshFormat is required.')
    end if
    !
    ! Read $PhysicalNames
    !
    section_name='PhysicalNames'
    call fbem_gmsh_search_section(fileunit,section_name,found)
    if (found) then
      ! Read numPhysicalNames
      read(fileunit,'(a)') linestr
      nwords=fbem_count_words(linestr)
      if (nwords.eq.1) then
        tmp_str=fbem_extract_word(linestr,1)
        read(tmp_str,*) mesh%n_physicalnames
      else
        call fbem_error_message(error_unit,0,trim(filename),0,'This line in $PhysicalNames must have 1 argument.')
      end if
      ! Check and initialize
      if (mesh%n_physicalnames.gt.0) then
        allocate (mesh%physicalname_dim (mesh%n_physicalnames))
        allocate (mesh%physicalname_eid (mesh%n_physicalnames))
        allocate (mesh%physicalname_name(mesh%n_physicalnames))
      else
        call fbem_error_message(error_unit,0,trim(filename),mesh%n_physicalnames,'the number of PhysicalNames must be >0')
      end if
      ! Read each PhysicalName: dimension, physicalTag, and "name"
      do i=1,mesh%n_physicalnames
        read(fileunit,'(a)') linestr
        nwords=fbem_count_words(linestr)
        if (nwords.eq.3) then
          ! dimension
          tmp_str=fbem_extract_word(linestr,1)
          read(tmp_str,*) mesh%physicalname_dim(i)
          if ((mesh%physicalname_dim(i).lt.0).or.(mesh%physicalname_dim(i).gt.3)) then
            call fbem_error_message(error_unit,0,trim(filename),i+1,'The dimension of a PhysicalName must be 0, 1, 2 or 3.')
          end if
          ! physicalTag
          tmp_str=fbem_extract_word(linestr,2)
          read(tmp_str,*) mesh%physicalname_eid(i)
          if (mesh%physicalname_eid(i).le.0) then
            call fbem_error_message(error_unit,0,trim(filename),i+1,'The tag of a PhysicalName must be >0.')
          end if
          ! name
          tmp_str=fbem_extract_word(linestr,3)
          mesh%physicalname_name(i)=trim(tmp_str)
        else
          call fbem_error_message(error_unit,0,trim(filename),i+1,'This line in $PhysicalNames must have 3 arguments.')
        end if
      end do
      ! Set and check external ids (tags)
      entity_str='physicalTag'
      call set_eid_to_iid(entity_str,mesh%n_physicalnames,mesh%physicalname_eid,&
                          mesh%physicalname_eid_min,mesh%physicalname_eid_max,mesh%physicalname_iid)
    else
      call fbem_error_message(error_unit,0,trim(filename),0,'The section $PhysicalNames is required.')
    end if
    !
    ! Read $Nodes
    !
    section_name='Nodes'
    call fbem_gmsh_search_section(fileunit,section_name,found)
    if (found) then
      ! Read numNodes
      read(fileunit,'(a)') linestr
      nwords=fbem_count_words(linestr)
      if (nwords.eq.1) then
        tmp_str=fbem_extract_word(linestr,1)
        read(tmp_str,*) mesh%n_nodes
      else
        call fbem_error_message(error_unit,0,trim(filename),1,'This line in $Nodes must have 1 argument.')
      end if
      ! Check and initialize
      if (mesh%n_nodes.gt.0) then
        allocate (mesh%node_eid(mesh%n_nodes))
        allocate (mesh%node_x(3,mesh%n_nodes))
      else
        call fbem_error_message(error_unit,0,trim(filename),1,'the number of Nodes must be >0.')
      end if
      ! Read each Node: nodeTag, and x, y, z
      do i=1,mesh%n_nodes
        read(fileunit,'(a)') linestr
        nwords=fbem_count_words(linestr)
        if (nwords.eq.4) then
          ! nodeTag
          tmp_str=fbem_extract_word(linestr,1)
          read(tmp_str,*) mesh%node_eid(i)
          if (mesh%node_eid(i).le.0) then
            call fbem_error_message(error_unit,0,trim(filename),i+1,'The tag of a Node must be >0.')
          end if
          ! x, y, z
          do j=2,4
            tmp_str=fbem_extract_word(linestr,j)
            call fbem_trimall(tmp_str)
            ! Transform the string to a double, so is read as double always.
            nc=scan(tmp_str,'eEdDqQ')
            if (nc.eq.0) then
              tmp_str(len_trim(tmp_str)+1:len_trim(tmp_str)+2)='D0'
            else
              nc=scan(tmp_str,'eE')
              if (nc.ne.0) tmp_str(nc:nc)='D'
            end if
            read(tmp_str,*) mesh%node_x(j-1,i)
          end do
        else
          call fbem_error_message(error_unit,0,trim(filename),i+1,'This line in $Nodes must have 4 arguments.')
        end if
      end do
      ! Set and check external ids (tags)
      entity_str='nodeTag'
      call set_eid_to_iid(entity_str,mesh%n_nodes,mesh%node_eid,mesh%node_eid_min,mesh%node_eid_max,mesh%node_iid)
    else
      call fbem_error_message(error_unit,0,trim(filename),0,'The section $Nodes is required.')
    end if
    !
    ! ELEMENTS
    !
    section_name='Elements'
    call fbem_gmsh_search_section(fileunit,section_name,found)
    if (found) then
      ! Read numElements
      read(fileunit,'(a)') linestr
      nwords=fbem_count_words(linestr)
      if (nwords.eq.1) then
        tmp_str=fbem_extract_word(linestr,1)
        read(tmp_str,*) mesh%n_elements
      else
        call fbem_error_message(error_unit,0,trim(filename),1,'This line in $Elements must have 1 argument.')
      end if
      ! Check and initialize
      if (mesh%n_elements.gt.0) then
        allocate (mesh%element_eid     (mesh%n_elements))
        allocate (mesh%element_type    (mesh%n_elements))
        allocate (mesh%element_physical(mesh%n_elements))
        allocate (mesh%element_node    (maxval(fbem_gmsh_type_n_nodes),mesh%n_elements))
      else
        call fbem_error_message(error_unit,0,trim(filename),mesh%n_nodes,'the number of Elements must be >0.')
      end if
      ! Read each element
      do i=1,mesh%n_elements
        read(fileunit,'(a)') linestr
        nwords=fbem_count_words(linestr)
        if (nwords.ge.4) then
          ! elementTag
          tmp_str=fbem_extract_word(linestr,1)
          read(tmp_str,*) mesh%element_eid(i)
          if (mesh%element_eid(i).le.0) then
            call fbem_error_message(error_unit,0,trim(filename),mesh%element_eid(i),'The tag of an Element must be >0.')
          end if
          ! elementType
          tmp_str=fbem_extract_word(linestr,2)
          read(tmp_str,*) mesh%element_type(i)
          if ((mesh%element_type(i).le.0).or.(mesh%element_type(i).gt.size(fbem_gmsh_type_n_nodes))) then
            call fbem_error_message(error_unit,0,trim(filename),mesh%element_eid(i),'Incorrect type of element.')
          end if
          ! numberTags
          tmp_str=fbem_extract_word(linestr,3)
          read(tmp_str,*) tmp_n_tags
          if ((tmp_n_tags.lt.1).or.(tmp_n_tags.gt.8)) then
            call fbem_error_message(error_unit,0,trim(filename),mesh%element_eid(i),'Incorrect number of tags for an element.')
          end if
          ! physicalTag
          tmp_str=fbem_extract_word(linestr,4)
          read(tmp_str,*) mesh%element_physical(i)
          if (findloc(mesh%physicalname_eid,mesh%element_physical(i),1).eq.0) then
            call fbem_error_message(error_unit,0,trim(filename),mesh%element_eid(i),'The indicated physical entity for this element does not exist.')
          end if
          ! Nodes
          if (nwords.ne.(3+tmp_n_tags+fbem_gmsh_type_n_nodes(mesh%element_type(i)))) then
            call fbem_error_message(error_unit,0,trim(filename),mesh%element_eid(i),'The number of nodes for this element is incorrect.')
          end if
          do j=1,fbem_gmsh_type_n_nodes(mesh%element_type(i))
            tmp_str=fbem_extract_word(linestr,3+tmp_n_tags+j)
            read(tmp_str,*) mesh%element_node(j,i)
            if (findloc(mesh%node_eid,mesh%element_node(j,i),1).eq.0) then
              call fbem_error_message(error_unit,0,trim(filename),mesh%element_eid(i),'An indicated element node does not exist for this element.')
            end if
          end do
        else
          call fbem_error_message(error_unit,0,trim(filename),i+1,'In this line of $Elements, the number of arguments must be >=4')
        end if
      end do
      ! Set and check external ids (tags)
      entity_str='elementTag'
      call set_eid_to_iid(entity_str,mesh%n_elements,mesh%element_eid,mesh%element_eid_min,mesh%element_eid_max,mesh%element_iid)
    else
      call fbem_error_message(error_unit,0,trim(filename),0,'The section Elements is required.')
    end if
    !
    ! CLOSE FILE
    !
    call fbem_close_file(filename,fileunit)
  end subroutine read

  subroutine fbem_export_gmsh_fmt_real(fmt_real_val)
    implicit none
    character(len=*) :: fmt_real_val
    fmt_real=fmt_real_val
  end subroutine

  subroutine fbem_export_gmsh_NodeData(fileunit,viewname,timeval,timestep,ncomp,nmax,nnode,nodeid,nodeval)
    implicit none
    ! I/O
    integer           :: fileunit
    character(len=*)  :: viewname
    real(kind=real64) :: timeval
    integer           :: timestep
    integer           :: ncomp
    integer           :: nmax
    integer           :: nnode
    integer           :: nodeid(nmax)
    real(kind=real64) :: nodeval(9,nmax)
    ! Local
    integer                    :: k
    character(len=fbem_fmtstr) :: fmt1, fmt2, fmt_integer
    !
    ! From gmsh documentation:
    !
    !  $NodeData
    !    numStringTags(ASCII int)
    !    stringTag(string) ...
    !    numRealTags(ASCII int)
    !    realTag(ASCII double) ...
    !    numIntegerTags(ASCII int)
    !    integerTag(ASCII int) ...
    !    nodeTag(size_t) value(double) ...
    !    ...
    !  $EndNodeData
    !
    ! check
    if ((ncomp.ne.1).and.(ncomp.ne.3).and.(ncomp.ne.9)) then
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid value of ncomp')
    end if
    write(fmt_integer,*) 'i', fbem_nchar_int(maxval(nodeid))+1
    call fbem_trimall(fmt_integer)
    ! write
    write(fileunit,'(a9)') '$NodeData'
    write(fileunit,'(a1)' ) '1'
    write(fmt1,*) '(a1,a',len_trim(viewname),',a1)'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) '"',viewname,'"'
    write(fileunit,'(a1)' ) '1'
    write(fmt1,*) '(',trim(fmt_real),')'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) timeval
    write(fileunit,'(a1)') '3'
    write(fmt1,*) '(i',fbem_nchar_int(timestep),')'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) timestep
    write(fileunit,'(i1)') ncomp
    write(fmt1,*) '(i',fbem_nchar_int(nnode),')'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) nnode
    do k=1,nnode
      write(fmt1,*) '(i',fbem_nchar_int(nodeid(k)),',',ncomp,trim(fmt_real),')'
      call fbem_trimall(fmt1)
      write(fileunit,fmt1) nodeid(k), nodeval(1:ncomp,k)
    end do
    write(fileunit,'(a12)') '$EndNodeData'
  end subroutine fbem_export_gmsh_NodeData

  subroutine fbem_export_gmsh_ElementData(fileunit,viewname,timeval,timestep,ncomp,nmax,nelement,elementid,elementval)
    implicit none
    ! I/O
    integer           :: fileunit
    character(len=*)  :: viewname
    real(kind=real64) :: timeval
    integer           :: timestep
    integer           :: ncomp
    integer           :: nmax
    integer           :: nelement
    integer           :: elementid(nmax)
    real(kind=real64) :: elementval(9,nmax)
    ! Local
    integer                    :: k
    character(len=fbem_fmtstr) :: fmt1, fmt2, fmt_integer
    !
    ! From gmsh documentation:
    !
    !  $ElementData
    !    numStringTags(ASCII int)
    !    stringTag(string) ...
    !    numRealTags(ASCII int)
    !    realTag(ASCII double) ...
    !    numIntegerTags(ASCII int)
    !    integerTag(ASCII int) ...
    !    elementTag(size_t) value(double) ...
    !    ...
    !  $EndElementData
    !
    ! check
    if ((ncomp.ne.1).and.(ncomp.ne.3).and.(ncomp.ne.9)) then
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid value of ncomp')
    end if
    write(fmt_integer,*) 'i', fbem_nchar_int(maxval(elementid))+1
    call fbem_trim2b(fmt_integer)
    ! write
    write(fileunit,'(a12)') '$ElementData'
    write(fileunit,'(a1)' ) '1'
    write(fmt1,*) '(a1,a',len_trim(viewname),',a1)'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) '"',viewname,'"'
    write(fileunit,'(a1)' ) '1'
    write(fmt1,*) '(',trim(fmt_real),')'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) timeval
    write(fileunit,'(a1)') '3'
    write(fmt1,*) '(i',fbem_nchar_int(timestep),')'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) timestep
    write(fileunit,'(i1)') ncomp
    write(fmt1,*) '(i',fbem_nchar_int(nelement),')'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) nelement
    do k=1,nelement
      write(fmt1,*) '(i',fbem_nchar_int(elementid(k)),',',ncomp,trim(fmt_real),')'
      call fbem_trimall(fmt1)
      write(fileunit,fmt1) elementid(k), elementval(1:ncomp,k)
    end do
    write(fileunit,'(a15)') '$EndElementData'
  end subroutine fbem_export_gmsh_ElementData

  subroutine fbem_export_gmsh_ElementNodeData(fileunit,viewname,timeval,timestep,ncomp,nemax,nnemax,nelem,elemid,nnodeselem,elemval)
    implicit none
    ! I/O
    integer           :: fileunit
    character(len=*)  :: viewname
    real(kind=real64) :: timeval
    integer           :: timestep
    integer           :: ncomp
    integer           :: nemax
    integer           :: nnemax
    integer           :: nelem
    integer           :: elemid(nemax)
    integer           :: nnodeselem(nemax)
    real(kind=real64) :: elemval(9,nnemax,nemax)
    ! Local
    integer                    :: k
    character(len=fbem_fmtstr) :: fmt1, fmt2, fmt_integer
    !
    ! From gmsh documentation:
    !
    ! $ElementNodeData
    !   numStringTags(ASCII int)
    !   stringTag(string) ...
    !   numRealTags(ASCII int)
    !   realTag(ASCII double) ...
    !   numIntegerTags(ASCII int)
    !   integerTag(ASCII int) ...
    !   elementTag(size_t) numNodesPerElement(int) value(double) ...
    !   ...
    ! $EndElementNodeData
    !
    ! check
    if ((ncomp.ne.1).and.(ncomp.ne.3).and.(ncomp.ne.9)) then
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid value of ncomp')
    end if
    write(fmt_integer,*) 'i', fbem_nchar_int(maxval(elemid))+1
    call fbem_trim2b(fmt_integer)
    ! write
    write(fileunit,'(a16)') '$ElementNodeData'
    write(fileunit,'(a1)' ) '1'
    write(fmt1,*) '(a1,a',len_trim(viewname),',a1)'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) '"',viewname,'"'
    write(fileunit,'(a1)' ) '1'
    write(fmt1,*) '(',trim(fmt_real),')'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) timeval
    write(fileunit,'(a1)') '3'
    write(fmt1,*) '(i',fbem_nchar_int(timestep),')'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) timestep
    write(fileunit,'(i2)') ncomp
    write(fmt1,*) '(i',fbem_nchar_int(nelem),')'
    call fbem_trimall(fmt1)
    write(fileunit,fmt1) nelem
    do k=1,nelem
      write(fmt1,*) '(i',fbem_nchar_int(elemid(k)),',1x,i',fbem_nchar_int(nnodeselem(k)),',',ncomp*nnodeselem(k),trim(fmt_real),')'
      call fbem_trimall(fmt1)
      write(fileunit,fmt1) elemid(k), nnodeselem(k), elemval(1:ncomp,1:nnodeselem(k),k)
    end do
    write(fileunit,'(a19)') '$EndElementNodeData'
  end subroutine fbem_export_gmsh_ElementNodeData

end module fbem_gmsh
