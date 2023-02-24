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
!! <b> Read and setup geometrical sensitivity analysis. </b>
subroutine read_sensitivity(fileunit)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_geometry
  use fbem_data_structures
  use fbem_mesh_module

  ! Problem variables module
  use problem_variables

  implicit none
  ! I/O
  integer                                 :: fileunit          !! Unit of the file to read from
  ! Local
  character(len=fbem_stdcharlen)          :: section_name      ! Name of the section
  logical                                 :: found             ! Logical variable for sections and keywords
  character(len=fbem_stdcharlen)          :: mode              ! Input mode (string)
  integer                                 :: mode_tag          ! Input mode tag
  integer                                 :: i, j, k, l, knj, snj, ka, kc, ki, kj
  character(len=fbem_filename_max_length) :: auxmesh_filename, auxdata_filename
  character(len=fbem_stdcharlen)          :: mesh_format
  real(kind=real64)                       :: delta_a
  integer                                 :: n_affectedparts
  integer, allocatable                    :: affectedpart(:)
  type(fbem_mesh)                         :: tmpmesh
  logical                                 :: good
  integer                                 :: sp, sn, method, max_d
  real(kind=real64), allocatable          :: xi(:), phi(:), dphidx(:,:)
  real(kind=real64)                       :: x(problem%n), centre(problem%n), radius, r, rv(problem%n), cl, d
  integer, allocatable                    :: node_n_elements(:)
  integer, allocatable                    :: node_element(:,:)
  integer                                 :: min_n, max_n
  logical                                 :: ok
  real(kind=real64)                       :: max_dxda, dxda

  ! <input mode>
  ! <n_designvariables>
  ! <read data for each design variable according to input mode>

  section_name='sensitivity'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! Read input mode
    read(fileunit,*) mode
    mode_tag=0
    if (trim(mode).eq.'mesh_movement' ) mode_tag=1
    if (trim(mode).eq.'auxiliary_mesh') mode_tag=2
    if (trim(mode).eq.'design_mesh'   ) mode_tag=3
    if (mode_tag.eq.0) then
      call fbem_error_message(error_unit,0,section_name,0,'the indicated input mode is not valid.')
    end if

    ! Read the number of design variables
    read(fileunit,*) problem%n_designvariables
     ! Allocate and initialize
    if (problem%n_designvariables.gt.0) then
      do i=1,n_nodes
        allocate (node(i)%dxda(problem%n,problem%n_designvariables))
        node(i)%dxda=0.d0
      end do
      do i=1,n_internalpoints
        allocate (internalpoint(i)%dxda(problem%n,problem%n_designvariables))
        allocate (internalpoint(i)%dnda(problem%n,problem%n,problem%n_designvariables))
        allocate (internalpoint(i)%d2xdadx(problem%n,problem%n,problem%n_designvariables))
        internalpoint(i)%dxda=0.d0
        internalpoint(i)%dnda=0.d0
        internalpoint(i)%d2xdadx=0.d0
      end do
    else
      call fbem_error_message(error_unit,0,section_name,0,'the number of design variables must be >0')
    end if

    select case (mode_tag)

      ! ===================================================================
      ! Assign dxda by a displaced mesh with respect to the design variable
      ! ===================================================================
      !
      !  <auxmesh format> <auxmesh filename> <Delta a> <n_affectedparts> <affectedpart(1:n)>

      case (1)
        ! Read and process the data for each design variable
        do i=1,problem%n_designvariables
          ! Read
          read(fileunit,*) mesh_format, auxmesh_filename, delta_a, n_affectedparts
          allocate (affectedpart(n_affectedparts))
          backspace(fileunit)
          read(fileunit,*) mesh_format, auxmesh_filename, delta_a, n_affectedparts, (affectedpart(j),j=1,n_affectedparts)
          ! Check mesh file path
          call fbem_trim(auxmesh_filename)
          if (.not.fbem_file_exists(auxmesh_filename)) then
            write(output_unit,'(a82)') 'Mesh file does not exist ([sensitivity] : auxmesh_filename), check the given path.'
            write(output_unit,*)
          end if
          if (fbem_path_is_relative(auxmesh_filename)) then
            auxmesh_filename=trim(input_filedir)//trim(auxmesh_filename)
          end if
          ! Build the mesh
          call tmpmesh%read_from_file(problem%n,geometric_tolerance,auxmesh_filename,mesh_format)
          ! Check of this auxiliary mesh is topologically identical to the base mesh
          good=.true.
          if (n_parts.ne.tmpmesh%n_parts) good=.false.
          if (n_elements.ne.tmpmesh%n_elements) good=.false.
          if (n_nodes.ne.tmpmesh%n_nodes) good=.false.
          ! Aqui faltarian chequear las connectividades, los tipos base de elementos, etc...
          if (.not.good) then
            call fbem_error_message(error_unit,0,auxmesh_filename,0,'the auxiliary mesh is not topologically equal to the base mesh')
          end if
          ! Chequear que las partes asignadas existen
          ! Chequear que los nodos tiene misma numeracion interna tambien ...
          do j=1,n_affectedparts
            sp=part_iid(affectedpart(j))
            do k=1,part(sp)%n_nodes
              sn=part(sp)%node(k)
              node(sn)%dxda(:,i)=(tmpmesh%node(sn)%x-node(sn)%x)/delta_a
            end do
          end do
          deallocate(affectedpart)
        end do

        write(error_unit,*) 'mesh_movement mode does not allow to calculate internal point sensitivities correctly'

      ! ====================================================================================
      ! Assign dxda according to an interpolated velocity field defined by an auxiliary mesh
      ! ====================================================================================
      !
      !  <auxmesh format> <auxmesh filename> <auxmesh dxda data> <n_affectedparts> <affectedpart(1:n)>

      case (2)
        ! Read and process the data for each design variable
        do i=1,problem%n_designvariables
          ! Read
          read(fileunit,*) mesh_format, auxmesh_filename, auxdata_filename, n_affectedparts
          allocate (affectedpart(n_affectedparts))
          backspace(fileunit)
          read(fileunit,*) mesh_format, auxmesh_filename, auxdata_filename, n_affectedparts, (affectedpart(j),j=1,n_affectedparts)
          ! Read and build the auxiliary mesh
          call tmpmesh%read_from_file(problem%n,geometric_tolerance,auxmesh_filename,mesh_format)
          ! Calculate some data on elements that is needed.
          do j=1,tmpmesh%n_elements
            allocate (tmpmesh%element(j)%x_gn(problem%n,tmpmesh%element(j)%n_nodes))
            do k=1,tmpmesh%element(j)%n_nodes
              sn=tmpmesh%element(j)%node(k)
              tmpmesh%element(j)%x_gn(:,k)=tmpmesh%node(sn)%x
            end do
            call fbem_geometry_element_ball(problem%n,tmpmesh%element(j)%type,tmpmesh%element(j)%x_gn,3,centre,radius)
            cl=fbem_characteristic_length(problem%n,tmpmesh%element(j)%type,tmpmesh%element(j)%x_gn,1.d-6)
            allocate (tmpmesh%element(j)%bball_centre(problem%n))
            tmpmesh%element(j)%bball_centre=centre
            tmpmesh%element(j)%bball_radius=radius
            tmpmesh%element(j)%csize=cl
          end do
          ! Read the velocity of design variables at nodes
          call tmpmesh%read_node_values_from_file(auxdata_filename)
          ! Assign interpolated values from the auxiliary mesh to the nodes of the affected parts
          do j=1,n_affectedparts
            sp=part_iid(affectedpart(j))
            do k=1,part(sp)%n_nodes
              sn=part(sp)%node(k)
              x=node(sn)%x
              do l=1,tmpmesh%n_elements
                rv=x-tmpmesh%element(l)%bball_centre
                r=sqrt(dot_product(rv,rv))
                if (r.lt.1.5d0*tmpmesh%element(l)%bball_radius) then
                  allocate (xi(tmpmesh%element(l)%n_dimension))
                  cl=tmpmesh%element(l)%csize
                  call fbem_nearest_element_point_bem(problem%n,tmpmesh%element(l)%type,tmpmesh%element(l)%x_gn,cl,x,xi,r,d,method)
                  if (d.lt.1.d-12) then
                    allocate (phi(tmpmesh%element(l)%n_nodes))
                    phi=fbem_phi_hybrid(tmpmesh%element(l)%type,0.d0,xi)
                    node(sn)%dxda(:,i)=0.d0
                    do knj=1,tmpmesh%element(l)%n_nodes
                      snj=tmpmesh%element(l)%node(knj)
                      node(sn)%dxda(:,i)=node(sn)%dxda(:,i)+phi(knj)*tmpmesh%node(snj)%value_r(:,1)
                    end do
                    deallocate (phi)
                  end if
                  deallocate (xi)
                end if
              end do
              !write(99,'(2i11,4e25.16)') i, node(sn)%id, node(sn)%x, node(sn)%dxda(:,i)
              !write(99,'(2i11,6e25.16)') i, node(sn)%id, node(sn)%x, node(sn)%dxda(:,i)
            end do
          end do
          deallocate(affectedpart)
        end do
        write(error_unit,*) 'auxiliary_mesh mode does not allow to calculate internal point sensitivities correctly'

      ! =========================================================================================================================
      ! Design mesh where dxda is not transfered to the problem mesh (isoparametric), but the shape functions of the design macro
      ! elements are integrated in the elements.
      ! =========================================================================================================================
      !
      !  <design mesh format> <design mesh filename> <design mesh data>
      case (3)
        ! Read design mesh and data
        read(fileunit,*) mesh_format, auxmesh_filename, auxdata_filename
        ! Read and build the design mesh
        call design_mesh%read_from_file(problem%n,geometric_tolerance,auxmesh_filename,mesh_format)
        ! Read the velocity field of design variables at nodes
        call design_mesh%read_node_values_from_file(auxdata_filename)
        ! Transfer the velocity field from value_r to dxda
        do j=1,design_mesh%n_nodes
          ! Check
          if (.not.allocated(design_mesh%node(j)%value_r)) then
            call fbem_error_message(error_unit,0,auxdata_filename,0,'the mesh data must the real')
          end if
          if (lbound(design_mesh%node(j)%value_r,1).ne.1) then
            call fbem_error_message(error_unit,0,auxdata_filename,0,'kdof_min must be 1')
          end if
          if (ubound(design_mesh%node(j)%value_r,1).ne.problem%n) then
            call fbem_error_message(error_unit,0,auxdata_filename,0,'kdof_max must be equal to the problem dimension')
          end if
          if (lbound(design_mesh%node(j)%value_r,2).ne.1) then
            call fbem_error_message(error_unit,0,auxdata_filename,0,'kgroup_min must be 1')
          end if
          if (ubound(design_mesh%node(j)%value_r,2).ne.problem%n_designvariables) then
            call fbem_error_message(error_unit,0,auxdata_filename,0,'kgroup_max mu be equal to the number of design variables')
          end if
          ! Allocate and copy
          allocate (design_mesh%node(j)%dxda(problem%n,problem%n_designvariables))
          design_mesh%node(j)%dxda=design_mesh%node(j)%value_r
        end do
        ! Calculate some needed geometrical data of the design mesh
        do j=1,design_mesh%n_elements
          allocate (design_mesh%element(j)%x_gn(problem%n,design_mesh%element(j)%n_nodes))
          do k=1,design_mesh%element(j)%n_nodes
            sn=design_mesh%element(j)%node(k)
            design_mesh%element(j)%x_gn(:,k)=design_mesh%node(sn)%x
          end do
          call fbem_geometry_element_ball(problem%n,design_mesh%element(j)%type,design_mesh%element(j)%x_gn,3,centre,radius)
          cl=fbem_characteristic_length(problem%n,design_mesh%element(j)%type,design_mesh%element(j)%x_gn,1.d-6)
          allocate (design_mesh%element(j)%bball_centre(problem%n))
          design_mesh%element(j)%bball_centre=centre
          design_mesh%element(j)%bball_radius=radius
          design_mesh%element(j)%csize=cl
        end do
        ! Build connectivity of the mesh with the design mesh
        ! Nodes -> Design Mesh
        do j=1,n_nodes
          x=node(j)%x
          ! Calculate node(j)%dm_n_elements
          node(j)%dm_n_elements=0
          max_d=0
          k=0
          do l=1,design_mesh%n_elements
            rv=x-design_mesh%element(l)%bball_centre
            r=sqrt(dot_product(rv,rv))
            if (r.lt.1.5d0*design_mesh%element(l)%bball_radius) then
              allocate (xi(design_mesh%element(l)%n_dimension))
              cl=design_mesh%element(l)%csize
              call fbem_nearest_element_point_bem(problem%n,design_mesh%element(l)%type,design_mesh%element(l)%x_gn,cl,x,xi,r,d,method)
              if (d.lt.1.d-12) then
                k=k+1
                if (design_mesh%element(l)%n_dimension.gt.max_d) max_d=design_mesh%element(l)%n_dimension
              end if
              deallocate (xi)
            end if
          end do
          node(j)%dm_n_elements=k
          if (node(j)%dm_n_elements.gt.0) then
            ! Build the connectivity
            allocate (node(j)%dm_element(k))
            allocate (node(j)%dm_element_d(k))
            allocate (node(j)%dm_element_xi(max_d,k))
            k=0
            do l=1,design_mesh%n_elements
              rv=x-design_mesh%element(l)%bball_centre
              r=sqrt(dot_product(rv,rv))
              if (r.lt.1.5d0*design_mesh%element(l)%bball_radius) then
                allocate (xi(design_mesh%element(l)%n_dimension))
                cl=design_mesh%element(l)%csize
                call fbem_nearest_element_point_bem(problem%n,design_mesh%element(l)%type,design_mesh%element(l)%x_gn,cl,x,xi,r,d,method)
                if (d.lt.1.d-12) then
                  k=k+1
                  max_d=design_mesh%element(l)%n_dimension
                  ! Save data
                  node(j)%dm_element(k)=l
                  node(j)%dm_element_d(k)=max_d
                  node(j)%dm_element_xi=0.d0
                  node(j)%dm_element_xi(1:max_d,k)=xi
                  ! Copy the velocity to the node dxda variable
                  allocate (phi(design_mesh%element(l)%n_nodes))
                  phi=fbem_phi_hybrid(design_mesh%element(l)%type,0.d0,xi)
                  do i=1,problem%n_designvariables
                    node(j)%dxda(:,i)=0.d0
                    do knj=1,design_mesh%element(l)%n_nodes
                      snj=design_mesh%element(l)%node(knj)
                      node(j)%dxda(:,i)=node(j)%dxda(:,i)+phi(knj)*design_mesh%node(snj)%dxda(:,i)
                    end do
                  end do
                  deallocate (phi)
                end if
                deallocate (xi)
              end if
            end do
          end if
!          ! Write
!          do i=1,problem%n_designvariables
!            write(99,'(2i11,4e25.16)') i, node(j)%id, node(j)%x, node(j)%dxda(:,i)
!            !write(99,'(2i11,6e25.16)') i, node(j)%id, node(j)%x, node(j)%dxda(:,i)
!          end do
        end do
        ! Elements -> Design Mesh
        do j=1,n_elements
          ! Initial check
          max_n=0
          min_n=1
          do k=1,element(j)%n_nodes
            snj=element(j)%node(k)
            if (node(snj)%dm_n_elements.gt.max_n) max_n=node(snj)%dm_n_elements
            if (node(snj)%dm_n_elements.lt.min_n) min_n=node(snj)%dm_n_elements
          end do
          ! All nodes are connected to at least one element of the design mesh
          if (min_n.gt.0) then
            ! Build auxiliary variables
            allocate (node_n_elements(element(j)%n_nodes),node_element(max_n,element(j)%n_nodes))
            do k=1,element(j)%n_nodes
              snj=element(j)%node(k)
              node_n_elements(k)=node(snj)%dm_n_elements
              do l=1,node_n_elements(k)
                node_element(l,k)=node(snj)%dm_element(l)
              end do
            end do
            ! Calculate common elements of the design mesh in the nodes of the element
            do i=1,node_n_elements(1)
              ! Selected element of the node 1
              snj=node_element(i,1)
              ! Explore if element snj is present in the rest of the nodes, if not set zero node_element(i,1)
              do k=2,element(j)%n_nodes
                ok=.false.
                do l=1,node_n_elements(k)
                  if (node_element(l,k).eq.snj) then
                    ok=.true.
                    exit
                  end if
                end do
                if (ok.eqv.(.false.)) then
                  node_element(i,1)=0
                  exit
                end if
              end do
            end do
            ! Extract the common elements from node_element(i,1)
            element(j)%dm_n_elements=0
            do i=1,node_n_elements(1)
              if (node_element(i,1).gt.0) element(j)%dm_n_elements=element(j)%dm_n_elements+1
            end do
            if (element(j)%dm_n_elements.gt.0) then
              allocate (element(j)%dm_element(element(j)%dm_n_elements))
              k=0
              do i=1,node_n_elements(1)
                if (node_element(i,1).gt.0) then
                  k=k+1
                  element(j)%dm_element(k)=node_element(i,1)
                end if
              end do
              ! The element is correctly embedded in an element of the design mesh
              if (element(j)%dm_n_elements.eq.1) then
                element(j)%dm_mode=1
              ! The element is the connected to >=2 design elements. It is set to be an isoparametric with respect to
              ! the design variables (dxda transferred to its nodes)
              else
                element(j)%dm_mode=0
              end if
              ! Save xi of each node
              max_d=0
              do i=1,element(j)%dm_n_elements
                snj=element(j)%dm_element(i)
                if (design_mesh%element(snj)%n_dimension.gt.max_d) max_d=design_mesh%element(snj)%n_dimension
              end do
              allocate (element(j)%dm_xi_gn(max_d,element(j)%n_nodes,element(j)%dm_n_elements))
              do i=1,element(j)%dm_n_elements
                snj=element(j)%dm_element(i)
                do k=1,element(j)%n_nodes
                  sn=element(j)%node(k)
                  do l=1,node(sn)%dm_n_elements
                    if (node(sn)%dm_element(l).eq.snj) then
                      element(j)%dm_xi_gn(1:node(sn)%dm_element_d(l),k,i)=node(sn)%dm_element_xi(1:node(sn)%dm_element_d(l),l)
                      exit
                    end if
                  end do
                end do
              end do
            else
              element(j)%dm_mode=0
              element(j)%dm_n_elements=0
            end if
            deallocate (node_n_elements,node_element)
          ! Some node is not connected to an element of the design mesh
          else
            element(j)%dm_mode=0
            element(j)%dm_n_elements=0
          end if
        end do
        ! Internal points -> Design Mesh
        do j=1,n_internalpoints
          x=internalpoint(j)%x
          ! Calculate internalpoint(j)%dm_n_elements
          internalpoint(j)%dm_n_elements=0
          max_d=0
          k=0
          do l=1,design_mesh%n_elements
            rv=x-design_mesh%element(l)%bball_centre
            r=sqrt(dot_product(rv,rv))
            if (r.lt.1.5d0*design_mesh%element(l)%bball_radius) then
              allocate (xi(design_mesh%element(l)%n_dimension))
              cl=design_mesh%element(l)%csize
              call fbem_nearest_element_point_bem(problem%n,design_mesh%element(l)%type,design_mesh%element(l)%x_gn,cl,x,xi,r,d,method)
              if (d.lt.1.d-12) then
                k=k+1
                if (design_mesh%element(l)%n_dimension.gt.max_d) max_d=design_mesh%element(l)%n_dimension
              end if
              deallocate (xi)
            end if
          end do
          internalpoint(j)%dm_n_elements=k
          if (internalpoint(j)%dm_n_elements.gt.0) then
            ! Build the connectivity
            allocate (internalpoint(j)%dm_element(k))
            allocate (internalpoint(j)%dm_element_d(k))
            allocate (internalpoint(j)%dm_element_xi(max_d,k))
            k=0
            do l=1,design_mesh%n_elements
              rv=x-design_mesh%element(l)%bball_centre
              r=sqrt(dot_product(rv,rv))
              if (r.lt.1.5d0*design_mesh%element(l)%bball_radius) then
                allocate (xi(design_mesh%element(l)%n_dimension))
                cl=design_mesh%element(l)%csize
                call fbem_nearest_element_point_bem(problem%n,design_mesh%element(l)%type,design_mesh%element(l)%x_gn,cl,x,xi,r,d,method)
                if (d.lt.1.d-12) then
                  k=k+1
                  max_d=design_mesh%element(l)%n_dimension
                  ! Save data
                  internalpoint(j)%dm_element(k)=l
                  internalpoint(j)%dm_element_d(k)=max_d
                  internalpoint(j)%dm_element_xi=0.d0
                  internalpoint(j)%dm_element_xi(1:max_d,k)=xi
                  ! Copy the velocity to the node dxda variable
                  allocate (phi(design_mesh%element(l)%n_nodes),dphidx(design_mesh%element(l)%n_nodes,problem%n))
                  phi=fbem_phi_hybrid(design_mesh%element(l)%type,0.d0,xi)
                  dphidx=fbem_dphidx(problem%n,design_mesh%element(l)%type,0.d0,design_mesh%element(l)%x_gn,xi)
                  do i=1,problem%n_designvariables
                    ! dxda
                    internalpoint(j)%dxda(:,i)=0.d0
                    do knj=1,design_mesh%element(l)%n_nodes
                      snj=design_mesh%element(l)%node(knj)
                      internalpoint(j)%dxda(:,i)=internalpoint(j)%dxda(:,i)+phi(knj)*design_mesh%node(snj)%dxda(:,i)
                    end do
                    ! d(dxda)/dxj
                    internalpoint(j)%d2xdadx(:,:,i)=0.d0
                    do ki=1,problem%n
                      do kj=1,problem%n
                        do knj=1,design_mesh%element(l)%n_nodes
                          snj=design_mesh%element(l)%node(knj)
                          internalpoint(j)%d2xdadx(ki,kj,i)=internalpoint(j)%d2xdadx(ki,kj,i)+dphidx(knj,kj)*design_mesh%node(snj)%dxda(ki,i)
                        end do
                      end do
                    end do
                  end do
                  deallocate (phi,dphidx)
                end if
                deallocate (xi)
              end if
            end do
          end if

          ! Build dnda
          do i=1,problem%n_designvariables
            internalpoint(j)%dnda(1,1,i)=0.d0
            internalpoint(j)%dnda(2,1,i)=-internalpoint(j)%d2xdadx(1,2,i)
            internalpoint(j)%dnda(1,2,i)=-internalpoint(j)%d2xdadx(2,1,i)
            internalpoint(j)%dnda(2,2,i)=0.d0
          end do
          if (problem%n.eq.3) stop 'internal points in sensitivity analysis not yet'
!          ! Write
!          do i=1,problem%n_designvariables
!            write(98,'(2i11,8e25.16)') i, internalpoint(j)%id, internalpoint(j)%x, internalpoint(j)%dxda(:,i), internalpoint(j)%dnda(:,1,i), internalpoint(j)%dnda(:,2,i)
!            !write(98,'(2i11,15e25.16)') i, internalpoint(j)%id, internalpoint(j)%x, internalpoint(j)%dxda(:,i), internalpoint(j)%dnda(:,1,i), internalpoint(j)%dnda(:,2,i), internalpoint(j)%dnda(:,3,i)
!          end do
        end do

    end select

    ! Calculate de maximum value of dxda
    max_dxda=0.d0
    do i=1,n_nodes
      do ka=1,problem%n_designvariables
        dxda=sqrt(dot_product(node(i)%dxda(:,ka),node(i)%dxda(:,ka)))
        if (dxda.gt.max_dxda) max_dxda=dxda
      end do
    end do

    ! Set the "sensitivity" flag to each element: false if dxda is zero in it, true otherwise.
    do i=1,n_elements
      element(i)%sensitivity=.false.
      do j=1,element(i)%n_nodes
        do ka=1,problem%n_designvariables
          do kc=1,problem%n
            if (abs(node(element(i)%node(j))%dxda(kc,ka)).gt.1.d-15*max_dxda) then
              element(i)%sensitivity=.true.
              exit
            end if
          end do
          if (element(i)%sensitivity) exit
        end do
        if (element(i)%sensitivity) exit
      end do
    end do

    ! Set the "sensitivity" flag to each region: false if dxda is zero in it, true otherwise.
    do i=1,n_regions
      region(i)%sensitivity=.false.
      select case (region(i)%class)
        case (fbem_be)
          do j=1,region(i)%n_boundaries
            do k=1,part(boundary(region(i)%boundary(j))%part)%n_elements
              if (element(part(boundary(region(i)%boundary(j))%part)%element(k))%sensitivity) then
                region(i)%sensitivity=.true.
                exit
              end if
            end do
            if (region(i)%sensitivity) exit
          end do
        case (fbem_fe)
          do j=1,region(i)%n_fe_subregions
            do k=1,part(fe_subregion(region(i)%fe_subregion(j))%part)%n_elements
              if (element(part(fe_subregion(region(i)%fe_subregion(j))%part)%element(k))%sensitivity) then
                region(i)%sensitivity=.true.
                exit
              end if
            end do
            if (region(i)%sensitivity) exit
          end do
      end select
    end do

    ! Set the "sensitivity" flag to each internal point (for now all, it should be false if dxda and dnda are null)
    do i=1,n_internalpoints
      internalpoint(i)%sensitivity=.true.
    end do

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  else
    call fbem_error_message(error_unit,0,'['//trim(section_name)//']',0,'this section is required')
  end if

end subroutine read_sensitivity
