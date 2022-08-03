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
!! BE WITH SOLID (N)D FE
!!   2D PROBLEM: SOLID FE (2D) - ORDINARY BE (1D)
!!   3D PROBLEM: SOLID FE (3D) - ORDINARY BE (2D)
!! BE WITH STRUCTURAL (N-1)D FE
!!   2D PROBLEM: BEAM FE (1D)  - ORDINARY OR CRACK-LIKE BE (1D)
!!   3D PROBLEM: SHELL FE (2D) - ORDINARY OR CRACK-LIKE BE (2D)
!!
!! BE - STRUCTURAL (N-1)D FE COUPLING: the BE is coupled with one or the two sides of the structure.
!! element - element
!! BE - SOLID FE COUPLING: the BE is coupled with a part of the boundary of a FE.
!! element - subedge/subface:
!! se_be: ordinary BE
!! se_fe: solid FE
!! sd_fe: edge of the solid FE (if 2D)
!! sf_fe: face of the solid FE (if 3D)
!!     - BE (1D or 2D):
!!       - se_fe=element(se_be)%element: coupled FE (0 if uncoupled).
!!       - sn_fe=element(se_be)%element_node(kn_be): node to node coupling from the BE point of view.
!!     - FE (2D):
!!       - sd_fe=element(se_fe)%subedge(kd_fe)
!!       - se_be=subedge(sd_fe)%element: coupled BE (0 if uncoupled).
!!       - sn_be=subedge(sd_fe)%element_node(kn_sd_fe): node to node coupling from the FE subedge point of view.
!!     - FE (3D):
!!       - sf_fe=element(se_fe)%subface(kf_fe)
!!       - se_be=subface(sf_fe)%element: coupled BE (0 if uncoupled).
!!       - sn_be=subface(sf_fe)%element_node(kn_sf_fe): node to node coupling from the FE subface point of view.


subroutine build_connectivity_be_fe(se)

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
  integer              :: se                        !! Selected BE
  ! Local variables
  integer              :: j, k, l, jf, lf           ! Counters
  integer              :: sn                        ! Selected element and node
  integer              :: se_fe                     ! Selected finite element
  integer              :: kn_be, sn_be              ! Selected boundary element node
  integer              :: kn_fe, sn_fe              ! Selected boundary element node
  integer              :: kn_cn, sn_cn              ! Selected common node
  integer              :: kse_fe, sse_fe            ! Selected subedges
  integer, allocatable :: n_fe_nodes(:)             ! Auxiliary variable to count the number of fe nodes in the same position
  integer, allocatable :: iid_fe_nodes(:,:)         ! Auxiliary variable to save the id of the fe nodes in the same position
  logical              :: uncoupled                 ! Auxiliary variable that is True if an element is totally uncoupled with a FE element
  integer              :: n_coupled_fe              ! Auxiliary variable that stores the number of coupled FE elements
  integer              :: tmp_iid_fe                ! Temporary variable that stores a FE element iid
  integer, allocatable :: iid_coupled_fe(:)         ! Auxiliary variable that stores the id of each couple FE element
  logical              :: exists                    ! Auxiliary variable to know if a given FE element exists in the common FE nodes


  ! Only for continuum BE
  if (element(se)%n_dimension.eq.(problem%n-1)) then

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
        !
        ! COUPLED WITH A SOLID FE
        !
        !   2D: BE - SOLID FE SUBEDGE
        !   3D: BE - SOLID FE SUBFACE
        !
        if (element(se)%n_dimension.eq.(element(se_fe)%n_dimension-1)) then
          ! Save to the BE data structure
          element(se)%element=se_fe
          allocate (element(se)%element_node(element(se)%n_nodes))
          element(se)%element_node=0
          do k=1,element(se)%n_nodes
            do j=1,n_fe_nodes(k)
              do l=1,node(iid_fe_nodes(j,k))%n_elements
                if (node(iid_fe_nodes(j,k))%element(l).eq.se_fe) then
                  element(se)%element_node(k)=iid_fe_nodes(j,k)
                  exit
                end if
              end do
              if (element(se)%element_node(k).ne.0) exit
            end do
            if (element(se)%element_node(k).eq.0) stop 'Fatal error: the program flow should not pass through here'
          end do
          ! Save to the FE subedge data structure
          select case (element(se_fe)%n_dimension)
            ! 2D SOLID FE SUBEDGE
            case (2)
              ! Allocate and initialize data structures
              do kse_fe=1,element(se_fe)%n_subedges
                sse_fe=element(se_fe)%subedge(kse_fe)
                if (allocated(subedge(sse_fe)%element_node).eqv.(.false.)) then
                  allocate (subedge(sse_fe)%element_node(subedge(sse_fe)%n_nodes))
                  subedge(sse_fe)%element=0
                  subedge(sse_fe)%element_node=0
                end if
              end do
              ! Find the connection between BE and FE subedges
              do kse_fe=1,element(se_fe)%n_subedges
                sse_fe=element(se_fe)%subedge(kse_fe)
                ! Only subedges that were unconnected to BE
                if (subedge(sse_fe)%element.eq.0) then
                  ! Node to node connection
                  do k=1,subedge(sse_fe)%n_nodes
                    do l=1,element(se)%n_nodes
                      if (subedge(sse_fe)%node(k).eq.element(se)%element_node(l)) then
                        subedge(sse_fe)%element_node(k)=element(se)%node(l)
                        exit
                      end if
                    end do
                  end do
                  ! Connect subedge if only all its nodes are connected
                  subedge(sse_fe)%element=se
                  do k=1,subedge(sse_fe)%n_nodes
                    if (subedge(sse_fe)%element_node(k).eq.0) then
                      subedge(sse_fe)%element=0
                      deallocate (subedge(sse_fe)%element_node)
                      exit
                    end if
                  end do
                  ! Check types
                  if (subedge(sse_fe)%element.ne.0) then
                    if (subedge(sse_fe)%type.ne.element(se)%type) then
                      stop 'Error: coupled BE - FE subedge must have the same base type'
                    end if
                  end if
                end if
              end do
              ! Check that the selected BE is connected to only one solid FE subedge
              k=0
              do kse_fe=1,element(se_fe)%n_subedges
                sse_fe=element(se_fe)%subedge(kse_fe)
                if (subedge(sse_fe)%element.eq.se) then
                  k=k+1
                end if
              end do
              if (k.ne.1) then
                call fbem_error_message(error_unit,0,'element',element(se)%id,&
                                        'something is wrong with the BE-FE coupling of this BE')
              end if
            ! 3D SOLID FE SUBFACE
            case (3)
              stop 'not yet 45' ! deberia ser exactamente igual que 2D pero cambiar subedge por subface
          end select
        end if
        !
        ! COUPLED WITH A STRUCTURAL FE
        !
        !   2D: BE - BEAM FE
        !   3D: BE - SHELL FE
        !
        if (element(se)%n_dimension.eq.element(se_fe)%n_dimension) then
          ! Save the iid of the FE element in the BE element, and viceversa.
          element(se)%element=se_fe
          element(se_fe)%element=se
          ! Check type
          if (element(se)%type.ne.element(se)%type) then
            stop 'Error: coupled BE - structural FE must have the same mesh type'
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
        ! ?¿?¿?¿? En principio bastaria con acoplar con el BE solo uno de los fe, p.ej. si es una lamina al borde de un solido fem
        ! acoplar solo con la lamina ?¿?¿?¿?¿?¿? pero entonces hay que acoplar la lamina y el solido. Mirar Mixed-dimensional
        ! coupling.

    end select

  end if

end subroutine build_connectivity_be_fe
