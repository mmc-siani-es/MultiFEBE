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

!!
!! This connectivity is done through their subelements
!!   - [BEAM TIP LOAD] BE body point load element (se) <-> FE subregion line element (se_fe). It is done through their nodes:
!!     - sn=element(se)%node(1) : the node of the point element
!!     - node(sn)%n_nodes, node(sn)%node(:) : connectivity of the node of the point element with other nodes that share its
!!                                            position.
!!     - sn_fe=node(sn)%node(kn_fe), se_fe=node(sn_fe)%element(1) : only one FE node must be present in the previous connectivity,
!!                                                                  and it must be a vertex node of a line (beam) finite element
!!                                                                  without other connections.
!!     - node(sn)%coupled_node=sn_fe, node(sn_fe)%coupled_node=sn : connectivity between both nodes
!!   - [BEAM/SHELL EDGE LOAD] BE body line load element (se) <-> FE subregion line/surface element (se_fe) (only 3D problems)
!!     For a beam, it is done like 2D DBEM/FEM coupling, i.e. element by element (like build_connectivity_be_fe).
!!     For a shell, it is done through their subedges:
!!     - sd=element(se)%subedge(1) : the edge of the line element.
!!     - subedge(sd)%n_nodes, subedge(sd)%node(:) : connectivity of the subedge of the line element with its nodes.
!!     - sd_fe is found by exploring the subedges of FE subregions of shell elements.
!!     - subedge(sd_fe)%n_supelements should be 1: if not, error.
!!     - subedge(sd_fe)%element=sd, subedge(sd)%element=sd_fe : connectivity between both subedges
!!
subroutine build_connectivity_be_bodyload_element_fe(se)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! No implicit allowed
  implicit none

  ! I/O variables
  integer              :: se                  !! Selected BE body load element
  ! Local variables
  integer              :: j, k, l, jf, lf, n  ! Counters
  integer              :: sn                  ! Selected element and node
  integer              :: se_fe, ke_fe        ! Selected finite element
  integer              :: kn_be, sn_be        ! Selected boundary element node
  integer              :: kn_fe, sn_fe        ! Selected boundary element node
  integer              :: kn_cn, sn_cn        ! Selected common node
  integer              :: sse, kse_fe, sse_fe ! Selected subedges
  integer, allocatable :: n_fe_nodes(:)      ! Auxiliary variable to count the number of fe nodes in the same position
  integer, allocatable :: iid_fe_nodes(:,:)  ! Auxiliary variable to save the id of the fe nodes in the same position
  logical              :: uncoupled          ! Auxiliary variable that is True if an element is totally uncoupled with a FE element
  integer              :: n_coupled_fe       ! Auxiliary variable that stores the number of coupled FE elements
  integer              :: tmp_iid_fe         ! Temporary variable that stores a FE element iid
  integer, allocatable :: iid_coupled_fe(:)  ! Auxiliary variable that stores the id of each couple FE element
  logical              :: exists             ! Auxiliary variable to know if a given FE element exists in the common FE nodes


  !
  ! Chequear que pasa cuando los elementos no son del mismo tipo, avisarlo como error......
  !


  ! ==================
  ! BE BODY POINT LOAD
  ! ==================
  !
  ! [BEAM TIP LOAD] BE body point load element (se) <-> FE subregion line element (se_fe)

  if (element(se)%n_dimension.eq.0) then
    ! Node of the BE point load
    sn=element(se)%node(1)
    ! Calculate the number of FE nodes connected with the node of the BE point load element
    n=0
    do j=1,node(sn)%n_nodes
      if (part(node(node(sn)%node(j))%part(1))%type.eq.fbem_part_fe_subregion) then
        n=n+1
        sn_fe=node(sn)%node(j)
        se_fe=node(sn_fe)%element(1)
      end if
    end do
    if (n.eq.0) then
      uncoupled=.false.
      node(sn)%coupled_node=0
      return
    else
      if ((n.eq.1).and.(node(sn_fe)%n_elements.eq.1).and.(element(se_fe)%n_dimension.eq.1).and.&
          (node(sn_fe)%element_node_loctype(1).eq.fbem_loctype_vertex)) then
        node(sn)%coupled_node=sn_fe
        node(sn_fe)%coupled_node=sn
      else
        stop 'this BE body point load is not allowed'
      end if
    end if
  end if

  ! =================
  ! BE BODY LINE LOAD
  ! =================
  !
  ! [BEAM/SHELL EDGE LOAD] BE body line load element (se) <-> FE subregion line/surface element (se_fe).

  if (element(se)%n_dimension.eq.1) then

    !
    ! Este trozo de codigo casi todo es igual que en build_connectivity_be_fe. Para el acople SHELL EDGE / LINE LOAD, hay ademas
    ! que identificar si el acople se produce solo sobre una de las aristas del shell.
    !

    ! Calculate the number of FE nodes of the nodes of the BE element
    allocate (n_fe_nodes(element(se)%n_nodes))
    do k=1,element(se)%n_nodes
      sn=element(se)%node(k)
      n_fe_nodes(k)=0
      do l=1,node(sn)%n_nodes
        if (part(node(node(sn)%node(l))%part(1))%type.eq.fbem_part_fe_subregion) then
          n_fe_nodes(k)=n_fe_nodes(k)+1
        end if
      end do
    end do

    ! If any BE node is not connected to a FE node, then the BE element is uncoupled.
    uncoupled=.false.
    do k=1,element(se)%n_nodes
      if (n_fe_nodes(k).eq.0) then
        uncoupled=.true.
        element(se)%element=0
        return
      end if
    end do

    ! If all BE nodes of the BE element have common FE nodes, it is necessary to continue.
    ! Save the id of each common FE node of the BE nodes of the BE element
    allocate (iid_fe_nodes(maxval(n_fe_nodes),element(se)%n_nodes))
    do k=1,element(se)%n_nodes
      sn=element(se)%node(k)
      n_fe_nodes(k)=0
      do l=1,node(sn)%n_nodes
        if (part(node(node(sn)%node(l))%part(1))%type.eq.fbem_part_fe_subregion) then
          n_fe_nodes(k)=n_fe_nodes(k)+1
          iid_fe_nodes(n_fe_nodes(k),k)=node(sn)%node(l)
        end if
      end do
    end do
    ! Calculate the number of common FE elements of all nodes of the BE element.
    n_coupled_fe=0
    ! The first node of the selected BE element is taken as the reference node to calculate the common FE elements
    sn=element(se)%node(1)
    do j=1,n_fe_nodes(1)
      ! Loop through the elements of the FE node of the selected BE node
      do l=1,node(iid_fe_nodes(j,1))%n_elements
        ! Selected FE element
        tmp_iid_fe=node(iid_fe_nodes(j,1))%element(l)
        uncoupled=.false.
        ! Find if the other FE nodes of the other nodes of the BE element belong to the same element.
        do k=2,element(se)%n_nodes
          exists=.false.
          do jf=1,n_fe_nodes(k)
            do lf=1,node(iid_fe_nodes(jf,k))%n_elements
              if (tmp_iid_fe.eq.node(iid_fe_nodes(jf,k))%element(lf)) then
                exists=.true.
                exit
              end if
            end do
          end do
          if (exists.eqv.(.false.)) then
            uncoupled=.true.
            exit
          end if
        end do
        if (uncoupled.eqv.(.false.)) then
          n_coupled_fe=n_coupled_fe+1
        end if
      end do
    end do

    ! If there is no coupled FE elements
    if (n_coupled_fe.eq.0) then
      element(se)%element=0
      return
    end if

    ! If there is coupled FE elements, it is necessary to continue.
    ! Save the id of each coupled FE element
    allocate (iid_coupled_fe(n_coupled_fe))
    n_coupled_fe=0
    ! The first node of the selected BE element is taken as the reference node to calculate the common FE elements
    sn=element(se)%node(1)
    do j=1,n_fe_nodes(1)
      do l=1,node(iid_fe_nodes(j,1))%n_elements
        ! Selected FE element
        tmp_iid_fe=node(iid_fe_nodes(j,1))%element(l)
        uncoupled=.false.
        ! Find if the other FE nodes of the other nodes of the BE element belong to the same element.
        do k=2,element(se)%n_nodes
          exists=.false.
          do jf=1,n_fe_nodes(k)
            do lf=1,node(iid_fe_nodes(jf,k))%n_elements
              if (tmp_iid_fe.eq.node(iid_fe_nodes(jf,k))%element(lf)) then
                exists=.true.
                exit
              end if
            end do
          end do
          if (exists.eqv.(.false.)) then
            uncoupled=.true.
            exit
          end if
        end do
        if (uncoupled.eqv.(.false.)) then
          n_coupled_fe=n_coupled_fe+1
          iid_coupled_fe(n_coupled_fe)=tmp_iid_fe
        end if
      end do
    end do
    ! Connect both elements
    select case (n_coupled_fe)
      !
      ! If it is uncoupled
      !
      case (0)
        stop 'Fatal error: the program flow should not pass through here'
      !
      ! If it is coupled with 1 FE element
      !
      case (1)
        ! Save its iid
        se_fe=iid_coupled_fe(1)
        select case (element(se_fe)%n_dimension)
          !
          ! LINE LOAD <---> BEAM (ELEMENT - ELEMENT)
          !
          case (1)
            !
            ! Save the iid of the FE element in the BE element, and viceversa.
            !
            element(se)%element=se_fe
            element(se_fe)%element=se
            ! Check
            if (element(se)%type.ne.element(se)%type) then
              stop 'Error: coupled BE_bodyload-FE elements must have the same base type'
            end if
            !
            ! Find the node connection between both elements
            !
            ! Allocate the node connection between both
            allocate (element(se)%element_node(element(se)%n_nodes))
            allocate (element(se_fe)%element_node(element(se_fe)%n_nodes))
            ! Initialize the connectivity vector
            element(se)%element_node=0
            element(se_fe)%element_node=0
            ! Loop through the nodes of the BE element
            do kn_be=1,element(se)%n_nodes
              ! Selected BE node
              sn_be=element(se)%node(kn_be)
              ! Loop through the common FE nodes of sn_be
              do kn_cn=1,node(sn_be)%n_nodes
                ! Selected common node
                sn_cn=node(sn_be)%node(kn_cn)
                ! Only if it is a FE node
                if (part(node(sn_cn)%part(1))%type.eq.fbem_part_fe_subregion) then
                  ! Loop through the nodes of the coupled FE element
                  do kn_fe=1,element(se_fe)%n_nodes
                    sn_fe=element(se_fe)%node(kn_fe)
                    ! If the nodes are the same, then they are connected
                    if (sn_fe.eq.sn_cn) then
                      element(se)%element_node(kn_be)=sn_fe
                      element(se_fe)%element_node(kn_fe)=sn_be
                      exit
                    end if
                  end do
                end if
              end do
            end do

          !
          ! LINE LOAD <---> SHELL EDGE (SUBEDGE - SUBEDGE)
          !
          case (2)
            ! Subedge of the line load element
            sse=element(se)%subedge(1)
            ! Allocate and initialize data structures
            allocate (subedge(sse)%element_node(subedge(sse)%n_nodes))
            subedge(sse)%element_node=0
            do kse_fe=1,element(se_fe)%n_subedges
              sse_fe=element(se_fe)%subedge(kse_fe)
              if (allocated(subedge(sse_fe)%element_node).eqv.(.false.)) then
                allocate (subedge(sse_fe)%element_node(subedge(sse_fe)%n_nodes))
                subedge(sse_fe)%element=0
                subedge(sse_fe)%element_node=0
              end if
            end do
            !
            ! Find the edge of the shell that connects both elements (see_fe)
            !
            ! Build shell edges nodes to line load nodes connectivity.
            ! Loop through the nodes of the BE element
            do kn_be=1,subedge(sse)%n_nodes
              ! Selected BE node
              sn_be=subedge(sse)%node(kn_be)
              ! Loop through the common FE nodes of sn_be
              do kn_cn=1,node(sn_be)%n_nodes
                ! Selected common node
                sn_cn=node(sn_be)%node(kn_cn)
                ! Only if it is a FE node
                if (part(node(sn_cn)%part(1))%type.eq.fbem_part_fe_subregion) then
                  ! Find if sn_cn belongs to se_fe
                  exists=.false.
                  do ke_fe=1,node(sn_cn)%n_elements
                    if (node(sn_cn)%element(ke_fe).eq.se_fe) then
                      exists=.true.
                      exit
                    end if
                  end do
                  ! If so, build shell edges nodes to line load nodes connectivity.
                  if (exists) then
                    ! Loop through the edges of se_fe
                    do kse_fe=1,element(se_fe)%n_subedges
                      sse_fe=element(se_fe)%subedge(kse_fe)
                      ! Only those that are not already connected are candidates
                      if (subedge(sse_fe)%element.eq.0) then
                        ! Loop through the nodes of the edge
                        do kn_fe=1,subedge(sse_fe)%n_nodes
                          sn_fe=subedge(sse_fe)%node(kn_fe)
                          ! If the nodes are the same, then they are connected
                          if (sn_fe.eq.sn_cn) then
                            subedge(sse_fe)%element_node(kn_fe)=sn_be
                            exit
                          end if
                        end do
                      end if
                    end do
                  end if
                end if
              end do
            end do
            ! Now, loop through the edges of se_fe, only that which all the nodes are connected is the possible one.
            uncoupled=.true.
            do kse_fe=1,element(se_fe)%n_subedges
              sse_fe=element(se_fe)%subedge(kse_fe)
              if (subedge(sse_fe)%element.eq.0) then
                n=0
                do kn_fe=1,subedge(sse_fe)%n_nodes
                  if (subedge(sse_fe)%element_node(kn_fe).gt.0) n=n+1
                end do
                if (n.eq.subedge(sse_fe)%n_nodes) then
                  uncoupled=.false.
                  exit
                end if
              end if
            end do
            ! Return if all edges are uncoupled
            if (uncoupled) return
            ! Check that the edge type is the same
            if (subedge(sse)%type.ne.subedge(sse_fe)%type) then
              stop 'Error: a coupled BE_bodyload-FE edge must have the same base type'
            end if
            ! Check that sse_fe is a boundary edge
            if (subedge(sse_fe)%n_supelements.ne.1) then
              stop 'Error: a coupled BE_bodyload-FE edge must be a boundary edge'
            end if
            !
            ! Save the iid of the coupling edges
            !
            subedge(sse)%element=sse_fe
            subedge(sse_fe)%element=sse
            ! Build connectivities between both
            ! Loop through the nodes of the BE element
            do kn_be=1,subedge(sse)%n_nodes
              ! Selected BE node
              sn_be=subedge(sse)%node(kn_be)
              ! Loop through the common FE nodes of sn_be
              do kn_cn=1,node(sn_be)%n_nodes
                ! Selected common node
                sn_cn=node(sn_be)%node(kn_cn)
                ! Only if it is a FE node
                if (part(node(sn_cn)%part(1))%type.eq.fbem_part_fe_subregion) then
                  ! Find if sn_cn belongs to se_fe
                  exists=.false.
                  do ke_fe=1,node(sn_cn)%n_elements
                    if (node(sn_cn)%element(ke_fe).eq.se_fe) then
                      exists=.true.
                      exit
                    end if
                  end do
                  ! If so
                  if (exists) then
                    ! Loop through the edges of se_fe
                    do kse_fe=1,element(se_fe)%n_subedges
                      ! Only the selected sse_fe is the coupled one
                      if (element(se_fe)%subedge(kse_fe).eq.sse_fe) then
                        do kn_fe=1,subedge(sse_fe)%n_nodes
                          sn_fe=subedge(sse_fe)%node(kn_fe)
                          ! If the nodes are the same, then they are connected
                          if (sn_fe.eq.sn_cn) then
                            subedge(sse)%element_node(kn_be)=sn_fe
                            subedge(sse_fe)%element_node(kn_fe)=sn_be
                            exit
                          end if
                        end do
                      end if
                    end do
                  end if
                end if
              end do
            end do

        end select
      !
      ! If 2
      !
      case default
        stop ' >= 2 coupled finite elements to a boundary element are not allowed'
        ! ?¿?¿?¿? En principio bastaria con acoplar con el BE solo uno de los fe, p.ej. si es una lamina al borde de un solido fem
        ! acoplar solo con la lamina ?¿?¿?¿?¿?¿?
    end select

  end if

  ! ====================
  ! BE BODY SURFACE LOAD
  ! ====================
  !
  ! [SHELL SURFACE LOAD] BE body surface load element (se) <-> FE subregion surface element (se_fe).

  if (element(se)%n_dimension.eq.2) then

    ! Calculate the number of FE nodes at the same position of each node of the BE
    allocate (n_fe_nodes(element(se)%n_nodes))
    do k=1,element(se)%n_nodes
      sn=element(se)%node(k)
      n_fe_nodes(k)=0
      do l=1,node(sn)%n_nodes
        if (part(node(node(sn)%node(l))%part(1))%type.eq.fbem_part_fe_subregion) then
          n_fe_nodes(k)=n_fe_nodes(k)+1
        end if
      end do
    end do

    ! If any BE node is not connected to a FE node, then the BE is uncoupled.
    uncoupled=.false.
    do k=1,element(se)%n_nodes
      if (n_fe_nodes(k).eq.0) then
        uncoupled=.true.
        element(se)%element=0
        return
      end if
    end do

    ! If all BE nodes of the BE share position with FE nodes, it is necessary to continue.
    ! Save the iid of each FE node sharing position with the BE nodes of the BE.
    allocate (iid_fe_nodes(maxval(n_fe_nodes),element(se)%n_nodes))
    do k=1,element(se)%n_nodes
      sn=element(se)%node(k)
      n_fe_nodes(k)=0
      do l=1,node(sn)%n_nodes
        if (part(node(node(sn)%node(l))%part(1))%type.eq.fbem_part_fe_subregion) then
          n_fe_nodes(k)=n_fe_nodes(k)+1
          iid_fe_nodes(n_fe_nodes(k),k)=node(sn)%node(l)
        end if
      end do
    end do

    ! Calculate the number of common FE of all nodes of the BE element.
    n_coupled_fe=0
    ! The first node of the selected BE element is taken as the reference node to calculate the common FE.
    sn=element(se)%node(1)
    do j=1,n_fe_nodes(1)
      ! Loop through the elements of the FE node of the selected BE node
      do l=1,node(iid_fe_nodes(j,1))%n_elements
        ! Selected FE element
        tmp_iid_fe=node(iid_fe_nodes(j,1))%element(l)
        uncoupled=.false.
        ! Find if the other FE nodes of the other nodes of the BE element belong to the same element.
        do k=2,element(se)%n_nodes
          exists=.false.
          do jf=1,n_fe_nodes(k)
            do lf=1,node(iid_fe_nodes(jf,k))%n_elements
              if (tmp_iid_fe.eq.node(iid_fe_nodes(jf,k))%element(lf)) then
                exists=.true.
                exit
              end if
            end do
          end do
          if (exists.eqv.(.false.)) then
            uncoupled=.true.
            exit
          end if
        end do
        if (uncoupled.eqv.(.false.)) then
          n_coupled_fe=n_coupled_fe+1
        end if
      end do
    end do

    ! If there is no coupled FE elements
    if (n_coupled_fe.eq.0) then
      element(se)%element=0
      return
    end if

    ! If there is coupled FE elements, it is necessary to continue.
    ! Save the id of each coupled FE element
    allocate (iid_coupled_fe(n_coupled_fe))
    n_coupled_fe=0
    ! The first node of the selected BE element is taken as the reference node to calculate the common FE elements
    sn=element(se)%node(1)
    do j=1,n_fe_nodes(1)
      do l=1,node(iid_fe_nodes(j,1))%n_elements
        ! Selected FE element
        tmp_iid_fe=node(iid_fe_nodes(j,1))%element(l)
        uncoupled=.false.
        ! Find if the other FE nodes of the other nodes of the BE element belong to the same element.
        do k=2,element(se)%n_nodes
          exists=.false.
          do jf=1,n_fe_nodes(k)
            do lf=1,node(iid_fe_nodes(jf,k))%n_elements
              if (tmp_iid_fe.eq.node(iid_fe_nodes(jf,k))%element(lf)) then
                exists=.true.
                exit
              end if
            end do
          end do
          if (exists.eqv.(.false.)) then
            uncoupled=.true.
            exit
          end if
        end do
        if (uncoupled.eqv.(.false.)) then
          n_coupled_fe=n_coupled_fe+1
          iid_coupled_fe(n_coupled_fe)=tmp_iid_fe
        end if
      end do
    end do


    ! Connect both elements
    select case (n_coupled_fe)

      ! =========
      ! UNCOUPLED
      ! =========

      case (0)
        stop 'Fatal error: the program flow should not pass through here'

      ! =================
      ! COUPLED WITH 1 FE
      ! =================

      case (1)
        ! Save the coupled FE iid
        se_fe=iid_coupled_fe(1)

        if (element(se)%n_dimension.eq.element(se_fe)%n_dimension) then
          if (problem%n.eq.2) then
            stop 'Error: be surface load superimposed to solid element'
          end if
          ! Save the iid of the FE element in the BE element, and viceversa.
          element(se)%element=se_fe
          element(se_fe)%element=se
          ! Check type
          if (element(se)%type.ne.element(se)%type) then
            stop 'Error: coupled BE load - structural FE must have the same base type'
          end if
          !
          ! Find the node connection between both elements
          !
          ! Allocate the node connection between both
          allocate (element(se)%element_node(element(se)%n_nodes))
          allocate (element(se_fe)%element_node(element(se_fe)%n_nodes))
          ! Initialize the connectivity vector
          element(se)%element_node=0
          element(se_fe)%element_node=0
          ! Loop through the nodes of the BE element
          do kn_be=1,element(se)%n_nodes
            ! Selected BE node
            sn_be=element(se)%node(kn_be)
            ! Loop through the common FE nodes of sn_be
            do kn_cn=1,node(sn_be)%n_nodes
              ! Selected common node
              sn_cn=node(sn_be)%node(kn_cn)
              ! Only if it is a FE node
              if (part(node(sn_cn)%part(1))%type.eq.fbem_part_fe_subregion) then
                ! Loop through the nodes of the coupled FE element
                do kn_fe=1,element(se_fe)%n_nodes
                  sn_fe=element(se_fe)%node(kn_fe)
                  ! If the nodes are the same, then they are connected
                  if (sn_fe.eq.sn_cn) then
                    element(se)%element_node(kn_be)=sn_fe
                    element(se_fe)%element_node(kn_fe)=sn_be
                    exit
                  end if
                end do
              end if
            end do
          end do
        end if

      ! ==================
      ! COUPLED WITH 2 FEs
      ! ==================
      case default
        stop ' >= 2 coupled finite elements to a boundary element are not allowed'

    end select

  end if

end subroutine build_connectivity_be_bodyload_element_fe
