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
!! <b> Subroutine that builds the data at functional nodes. </b>

subroutine build_data_at_functional_nodes

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_data_structures
  use fbem_geometry

  ! Module of problem variables
  use problem_variables

  ! Local variables
  implicit none
  integer                        :: kb, sp, ke, se, kn, sn, kln, kc, kej, sej
  integer                        :: tmp_n_gnodes, tmp_n_f2nodes, tmp_n_dimension
  integer                        :: ss1, ss2
  real(kind=real64)              :: xi_1d, xi_2d(2)
  real(kind=real64)              :: v1(problem%n), v2(problem%n), v3(problem%n)
  real(kind=real64)              :: tbp(3), esp(3)
  real(kind=real64), allocatable :: x_nodes(:,:)
  real(kind=real64), allocatable :: x_fn(:)
  real(kind=real64), allocatable :: n_fn(:)
  real(kind=real64), allocatable :: xi_nd(:)
  real(kind=real64)              :: norm

  ! Starting message
  if (verbose_level.ge.1)  write(output_unit,'(a)') 'Building element data at functional nodes ...'

  ! Allocate
  allocate (x_fn(problem%n),n_fn(problem%n))

  ! ================================================================================================================================
  ! ELEMENT-WISE VECTORS AT BE FUNCTIONAL NODES
  ! ================================================================================================================================

  ! Loop through BOUNDARY ELEMENTS
  do ke=1,n_elements
    if (part(element(ke)%part)%type.eq.fbem_part_be_boundary) then

      ! Copy x_gn to x_nodes
      allocate (x_nodes(problem%n,element(ke)%n_nodes))
      do kn=1,element(ke)%n_nodes
        do kc=1,problem%n
          x_nodes(kc,kn)=element(ke)%x_gn(kc,kn)
        end do
      end do
      ! Allocate xi vector
      allocate (xi_nd(element(ke)%n_dimension))

      ! -----
      ! xi_fn
      ! -----

      allocate (element(ke)%xi_fn(element(ke)%n_dimension,element(ke)%n_nodes))
      element(ke)%xi_fn=fbem_xi_hybrid(element(ke)%type,element(ke)%delta_f)

      ! ----
      ! x_fn
      ! ----

      allocate (element(ke)%x_fn(problem%n,element(ke)%n_nodes))
      do kn=1,element(ke)%n_nodes
        do kc=1,element(ke)%n_dimension
          xi_nd(kc)=element(ke)%xi_fn(kc,kn)
        end do
        x_fn=fbem_position(problem%n,element(ke)%type,x_nodes,xi_nd)
        do kc=1,problem%n
          element(ke)%x_fn(kc,kn)=x_fn(kc)
        end do
      end do

      ! ----
      ! n_fn
      ! ----

      ! It is only done for 1D elements in 2D and 2D elements in 3D, otherwise n_gn is a zero vector.
      allocate (element(ke)%n_fn(problem%n,element(ke)%n_nodes))
      do kn=1,element(ke)%n_nodes
        do kc=1,problem%n
          element(ke)%n_fn(kc,kn)=0.0d0
        end do
      end do
      ! 1D element in 2D
      if ((element(ke)%n_dimension.eq.1).and.(problem%n.eq.2)) then
        do kn=1,element(ke)%n_nodes
          xi_1d=element(ke)%xi_fn(1,kn)
          n_fn=fbem_unormal2d(element(ke)%type,x_nodes,xi_1d)
          do kc=1,2
            element(ke)%n_fn(kc,kn)=n_fn(kc)
          end do
        end do
      end if
      ! 2D elements in 3D
      if ((element(ke)%n_dimension.eq.2).and.(problem%n.eq.3)) then
        do kn=1,element(ke)%n_nodes
          xi_2d(1)=element(ke)%xi_fn(1,kn)
          xi_2d(2)=element(ke)%xi_fn(2,kn)
          n_fn=fbem_unormal3d(element(ke)%type,x_nodes,xi_2d)
          do kc=1,3
            element(ke)%n_fn(kc,kn)=n_fn(kc)
          end do
        end do
      end if

      ! Deallocate
      deallocate (x_nodes,xi_nd)

    end if
  end do ! Loop through BOUNDARY ELEMENTS

  ! ================================================================================================================================
  ! ELEMENT-WISE VECTORS AT COUPLED BE BODY LOAD ELEMENT FUNCTIONAL NODES
  ! ================================================================================================================================

  ! Loop through COUPLED BE BODY LOAD ELEMENTS
  do kb=1,n_be_bodyloads
    sp=be_bodyload(kb)%part
    select case (be_bodyload(kb)%coupling)

      ! ----------------------------------
      ! FE BEAM TIP - BE LINE/SURFACE LOAD
      ! ----------------------------------

      case (fbem_bl_coupling_beam_tip)
        select case (problem%n)
          !
          ! 2D: LINE LOAD AT THE BEAM TIP
          !
          case (2)
            do ke=1,part(sp)%n_elements
              se=part(sp)%element(ke)

              ! Geometrical description
              ! Copy x_gn to x_nodes
              tmp_n_gnodes=fbem_n_nodes(element(se)%type_g)
              tmp_n_dimension=fbem_n_dimension(element(se)%type_g)
              allocate (x_nodes(problem%n,tmp_n_gnodes))
              allocate (xi_nd(tmp_n_dimension))
              x_nodes=element(se)%x_gn

              ! Functional interpolation of the load
              element(se)%type_f2=fbem_point1
              tmp_n_f2nodes=fbem_n_nodes(element(se)%type_f2)
              allocate (element(se)%xi_fn(tmp_n_dimension,tmp_n_f2nodes)) ! la dimension de un elemento constante es igual a su soporte geometrico
              allocate (element(se)%x_fn(problem%n,tmp_n_f2nodes))
              if (tmp_n_f2nodes.eq.1) then
                element(se)%xi_fn(:,1)=0.d0 ! si fuese un triangulo (2D) seria (1/3,1/3)
                element(se)%x_fn(:,1)=node(element(se)%node(1))%x
              else
                element(se)%xi_fn=fbem_xi_hybrid(element(se)%type_f2,element(se)%delta_f)
                do kn=1,tmp_n_f2nodes
                  xi_nd=element(se)%xi_fn(:,kn)
                  x_fn=fbem_position(problem%n,element(se)%type_g,x_nodes,xi_nd)
                  element(se)%x_fn(:,kn)=x_fn
                end do
              end if

              ! Deallocate
              deallocate (x_nodes,xi_nd)

            end do
          !
          ! 3D: SURFACE LOAD ATH THE BEAM TIP
          !
          case (3)
            ! Modelo de Luis (en principio solo habria que crear un cuadrilatero que emule a una circunferencia, se puede hacer..)
            stop 'not yet 117'
        end select

      ! -----------------------------------------------------
      ! FE BEAM - BE LINE LOAD AND FE SHELL - BE SURFACE LOAD
      ! -----------------------------------------------------

      case (fbem_bl_coupling_beam_line,fbem_bl_coupling_shell_surface)

        do ke=1,part(sp)%n_elements
          se=part(sp)%element(ke)

          ! Copy x_gn to x_nodes
          allocate (x_nodes(problem%n,element(se)%n_nodes))
          do kn=1,element(se)%n_nodes
            do kc=1,problem%n
              x_nodes(kc,kn)=element(se)%x_gn(kc,kn)
            end do
          end do
          ! Allocate xi vector
          allocate (xi_nd(element(se)%n_dimension))

          ! -----
          ! xi_fn
          ! -----

          allocate (element(se)%xi_fn(element(se)%n_dimension,element(se)%n_nodes))
          element(se)%xi_fn=fbem_xi_hybrid(element(se)%type,element(se)%delta_f)

          ! ----
          ! x_fn
          ! ----

          allocate (element(se)%x_fn(problem%n,element(se)%n_nodes))
          do kn=1,element(se)%n_nodes
            do kc=1,element(se)%n_dimension
              xi_nd(kc)=element(se)%xi_fn(kc,kn)
            end do
            x_fn=fbem_position(problem%n,element(se)%type,x_nodes,xi_nd)
            do kc=1,problem%n
              element(se)%x_fn(kc,kn)=x_fn(kc)
            end do
          end do

          ! ----
          ! n_fn
          ! ----

          ! It is only done for 1D elements in 2D and 2D elements in 3D, otherwise n_gn is a zero vector.
          allocate (element(se)%n_fn(problem%n,element(se)%n_nodes))
          do kn=1,element(se)%n_nodes
            do kc=1,problem%n
              element(se)%n_fn(kc,kn)=0.0d0
            end do
          end do
          ! 1D element in 2D
          if ((element(se)%n_dimension.eq.1).and.(problem%n.eq.2)) then
            do kn=1,element(se)%n_nodes
              xi_1d=element(se)%xi_fn(1,kn)
              n_fn=fbem_unormal2d(element(se)%type,x_nodes,xi_1d)
              do kc=1,2
                element(se)%n_fn(kc,kn)=n_fn(kc)
              end do
            end do
          end if
          ! 2D elements in 3D
          if ((element(se)%n_dimension.eq.2).and.(problem%n.eq.3)) then
            do kn=1,element(se)%n_nodes
              xi_2d(1)=element(se)%xi_fn(1,kn)
              xi_2d(2)=element(se)%xi_fn(2,kn)
              n_fn=fbem_unormal3d(element(se)%type,x_nodes,xi_2d)
              do kc=1,3
                element(se)%n_fn(kc,kn)=n_fn(kc)
              end do
            end do
          end if

          ! Deallocate
          deallocate (x_nodes,xi_nd)

        end do

      ! -----------------------
      ! FE SHELL - BE EDGE LOAD
      ! -----------------------

      case (fbem_bl_coupling_shell_edge)
        do ke=1,part(sp)%n_elements
          se=part(sp)%element(ke)

          ! Geometrical description of the element
          ! Copy x_gn to x_nodes
          tmp_n_gnodes=fbem_n_nodes(element(se)%type_g)
          tmp_n_dimension=fbem_n_dimension(element(se)%type_g)
          allocate (x_nodes(problem%n,tmp_n_gnodes))
          allocate (xi_nd(tmp_n_dimension))
          x_nodes=element(se)%x_gn

          ! Functional interpolation of the load
          select case (element(se)%type)
            case (fbem_line2)
              element(se)%type_f2=fbem_line2point1
            case (fbem_line3)
              element(se)%type_f2=fbem_line3point1
          end select


          tmp_n_f2nodes=fbem_n_nodes(element(se)%type_f2)
          allocate (element(se)%xi_fn(tmp_n_dimension,tmp_n_f2nodes)) ! la dimension de un elemento constante es igual a su soporte geometrico
          allocate (element(se)%x_fn(problem%n,tmp_n_f2nodes))
          if (tmp_n_f2nodes.eq.1) then
            element(se)%xi_fn(:,1)=0.d0 ! si fuese un triangulo (2D) seria (1/3,1/3)
            element(se)%x_fn(:,1)=node(element(se)%node(1))%x
          else
            element(se)%xi_fn=fbem_xi_hybrid(element(se)%type_f2,element(se)%delta_f)
            do kn=1,tmp_n_f2nodes
              xi_nd=element(se)%xi_fn(:,kn)
              x_fn=fbem_position(problem%n,element(se)%type_g,x_nodes,xi_nd)
              element(se)%x_fn(:,kn)=x_fn
            end do
          end if

          ! Deallocate
          deallocate (x_nodes,xi_nd)

        end do

    end select
  end do ! Loop through COUPLED BE BODY LOAD ELEMENTS

  ! ================================================================================================================================
  ! VECTORS AT BE FUNCTIONAL NODES
  ! ================================================================================================================================

  ! Loop through BE NODES
  do kn=1,n_nodes
    if (part(node(kn)%part(1))%type.eq.fbem_part_be_boundary) then

      ! -----
      ! x_fn
      ! ----

      ! Save to the structure
      se=node(kn)%element(1)
      kln=node(kn)%element_node_iid(1)
      node(kn)%x_fn=element(se)%x_fn(:,kln)

      ! -----
      ! n_fn
      ! ----
      !
      ! The unit normal at the functional node is the unit average normal of all elements that contain the node.

      ! If an element containing the node is continuous (then all elements are continuous)
      if (element(node(kn)%element(1))%discontinuous.eqv.(.false.)) then
        ! Initialize
        n_fn=0.0d0
        ! Loop through the ELEMENTS of the NODE
        do ke=1,node(kn)%n_elements
          se=node(kn)%element(ke)
          kln=node(kn)%element_node_iid(ke)
          n_fn=n_fn+element(se)%n_fn(:,kln)
        end do
        ! If the functional node belongs to any symmetry plane, it is necessary to include the normals of symmetrical elements.
        select case (node(kn)%n_symplanes)
          ! If it belongs to 1 symmetry plane (it can happen in 2D and 3D)
          case (1)
            ! Symmetry plane of the node
            ss1=node(kn)%symplane(1)
            ! Reflect and add the normal with reversed orientation.
            do kc=1,problem%n
              n_fn(kc)=n_fn(kc)+symplane_m(kc,ss1)*n_fn(kc)
            end do
          ! If it belongs to 2 symmetry planes (it can happen only in 3D)
          case (2)
            ! Symmetry planes of the node
            ss1=node(kn)%symplane(1)
            ss2=node(kn)%symplane(2)
            ! Reflect and add the normal with reversed orientation with respect to the first symmetry plane.
            do kc=1,problem%n
              n_fn(kc)=n_fn(kc)+symplane_m(kc,ss1)*n_fn(kc)
            end do
            ! Reflect and add the normal with reversed orientation with respect to the second symmetry plane.
            do kc=1,problem%n
              n_fn(kc)=n_fn(kc)+symplane_m(kc,ss2)*n_fn(kc)
            end do
        end select
        ! Normalize
        n_fn=n_fn/dsqrt(dot_product(n_fn,n_fn))
        ! Save to the structure
        node(kn)%n_fn=n_fn
      ! If the element is discontinuous
      else
        se=node(kn)%element(1)
        kln=node(kn)%element_node_iid(1)
        node(kn)%n_fn=element(se)%n_fn(:,kln)
      end if

      ! ----------------
      ! t1_fn and t2_fn
      ! ---------------

      ! For 2D, the unit tangent 1 at the functional node is t1 = e3 x n.
      if (problem%n.eq.2) then
        node(kn)%t1_fn(1) =-node(kn)%n_fn(2)
        node(kn)%t1_fn(2) = node(kn)%n_fn(1)
        node(kn)%t1_fn=node(kn)%t1_fn/sqrt(dot_product(node(kn)%t1_fn,node(kn)%t1_fn))
      end if
      ! For 3D, a reference vector can be introduced by the user (it is stored in t1). This reference vector together with n let to
      ! define a plane where t1 and n lies, being t2 its normal. If this reference vector is not introduced, then e1, e2 and e3 are
      ! used in this order, respectively.
      if (problem%n.eq.3) then
        ! Build t2 = n x reference vector (if any)
        node(kn)%t2_fn(1) = node(kn)%n_fn(2)*node(kn)%t1_fn(3)-node(kn)%n_fn(3)*node(kn)%t1_fn(2)
        node(kn)%t2_fn(2) = node(kn)%n_fn(3)*node(kn)%t1_fn(1)-node(kn)%n_fn(1)*node(kn)%t1_fn(3)
        node(kn)%t2_fn(3) = node(kn)%n_fn(1)*node(kn)%t1_fn(2)-node(kn)%n_fn(2)*node(kn)%t1_fn(1)
        norm=dsqrt(dot_product(node(kn)%t2_fn,node(kn)%t2_fn))
        if (norm.le.geometric_tolerance) then
          ! e1 as reference vector
          node(kn)%t2_fn(1) = 0.0d0
          node(kn)%t2_fn(2) = node(kn)%n_fn(3)
          node(kn)%t2_fn(3) =-node(kn)%n_fn(2)
          norm=dsqrt(dot_product(node(kn)%t2_fn,node(kn)%t2_fn))
          if (norm.le.geometric_tolerance) then
            ! e2 as reference vector
            node(kn)%t2_fn(1) =-node(kn)%n_fn(3)
            node(kn)%t2_fn(2) = 0.0d0
            node(kn)%t2_fn(3) = node(kn)%n_fn(1)
            norm=dsqrt(dot_product(node(kn)%t2_fn,node(kn)%t2_fn))
            if (norm.le.geometric_tolerance) then
              ! e3 as reference vector
              node(kn)%t2_fn(1) = node(kn)%n_fn(2)
              node(kn)%t2_fn(2) =-node(kn)%n_fn(1)
              node(kn)%t2_fn(3) = 0.0d0
              norm=dsqrt(dot_product(node(kn)%t2_fn,node(kn)%t2_fn))
              if (norm.le.geometric_tolerance) then
                stop 'Fatal error. Check subroutine build_data_at_functional_nodes() for BE.'
              end if
            end if
          end if
        end if
        node(kn)%t2_fn=node(kn)%t2_fn/sqrt(dot_product(node(kn)%t2_fn,node(kn)%t2_fn))
        ! Build t1 = t2 x n
        node(kn)%t1_fn(1) = node(kn)%t2_fn(2)*node(kn)%n_fn(3)-node(kn)%t2_fn(3)*node(kn)%n_fn(2)
        node(kn)%t1_fn(2) = node(kn)%t2_fn(3)*node(kn)%n_fn(1)-node(kn)%t2_fn(1)*node(kn)%n_fn(3)
        node(kn)%t1_fn(3) = node(kn)%t2_fn(1)*node(kn)%n_fn(2)-node(kn)%t2_fn(2)*node(kn)%n_fn(1)
        node(kn)%t1_fn=node(kn)%t1_fn/sqrt(dot_product(node(kn)%t1_fn,node(kn)%t1_fn))
      end if

    end if
  end do ! Loop through BE NODES

  ! ================================================================================================================================
  ! DIRECTOR VECTORS (LOCAL BENDING ROTATIONS) OF SHELL NODES (v_1 and v_2) DEGENERATED FINITE ELEMENTS
  ! ================================================================================================================================

  ! Loop through FINITE ELEMENTS
  do ke=1,n_elements
    if (part(element(ke)%part)%type.eq.fbem_part_fe_subregion) then

      ! Not required for rigid regions/subregions
      if (fe_subregion(part(element(ke)%part)%entity)%master_node.ne.0) cycle

      select case (problem%n)
        case (2)
          select case (element(ke)%n_dimension)

            ! -------------
            ! POINT ELEMENTS
            ! -------------

            case (0)
              ! Nothing to do

            ! -------------
            ! BEAM ELEMENTS
            ! -------------

            case (1)
              ! Nothing to do

            ! --------------
            ! SOLID ELEMENTS
            ! --------------

            case (2)
              ! Nothing to do
          end select



        case (3)
          select case (element(ke)%n_dimension)

            ! -------------
            ! POINT ELEMENTS
            ! -------------

            case (0)
              ! Nothing to do

            ! -------------
            ! BEAM ELEMENTS
            ! -------------

            case (1)
              ! Nothing to do


            ! --------------
            ! SHELL ELEMENTS
            ! --------------

            case (2)
              do kn=1,element(ke)%n_nodes
                sn=element(ke)%node(kn)
                !
                ! SINGULAR NODE
                !
                ! If the node belongs to a joint (a singular point whose edge belongs to >2 elements), the v1 and v2 vectors
                ! are obtained using a simple global criteria starting from v3. Since this nodes are converted to 6 DOF with
                ! global rotations, this node at different elements does not have to have equal v1 and v2. Nothing more must be
                ! done for this nodes when they belong to symmetry planes.
                if (node(sn)%is_singular) then
                  ! e1 as reference vector: v2 = v3 x e1
                  v2(1) = 0.0d0
                  v2(2) = element(ke)%v_midnode(3,3,kn)
                  v2(3) =-element(ke)%v_midnode(2,3,kn)
                  norm=dsqrt(dot_product(v2,v2))
                  if (norm.le.geometric_tolerance) then
                    ! e2 as reference vector: v2 = v3 x e2
                    v2(1) =-element(ke)%v_midnode(3,3,kn)
                    v2(2) = 0.0d0
                    v2(3) = element(ke)%v_midnode(1,3,kn)
                    norm=dsqrt(dot_product(v2,v2))
                    if (norm.le.geometric_tolerance) then
                      ! e3 as reference vector: v2 = v3 x e3
                      v2(1) = element(ke)%v_midnode(2,3,kn)
                      v2(2) =-element(ke)%v_midnode(1,3,kn)
                      v2(3) = 0.0d0
                      norm=dsqrt(dot_product(v2,v2))
                      if (norm.le.geometric_tolerance) then
                        stop 'Fatal error. Check subroutine build_data_at_functional_nodes() for FE.'
                      end if
                    end if
                  end if
                  v2=v2/dsqrt(dot_product(v2,v2))
                  ! v1 = v2 x v3
                  v1=fbem_cross_product(v2,element(ke)%v_midnode(:,3,kn))
                  v1=v1/sqrt(dot_product(v1,v1))
                !
                ! REGULAR NODE
                !
                ! If the node belong to a regular point (not a singular point), then it is a 5 DOF node with local rotations.
                !
                else
                  !
                  ! The node is at the boundary. Rotation alpha (direction v2) is tangent with the boundary, and rotation beta
                  ! (direction v1) is hence perpendicular to it and v3.
                  !
                  if (node(sn)%in_boundary) then
                    !
                    ! If the node belong to any symmetry plane
                    !
                    if (node(sn)%n_symplanes.gt.0) then
                      !
                      ! Check that v3 is orthogonal to the normals of the symmetry planes.
                      !
                      do kej=1,node(sn)%n_symplanes
                        ! Normal vector of the symmetry plane
                        esp=0.d0
                        esp(symplane_eid(node(sn)%symplane(kej)))=1.d0
                        if (abs(dot_product(element(ke)%v_midnode(:,3,kn),esp)).gt.geometric_tolerance) then
                          stop 'v3 vector should be orthogonal to the normal of all symmetry planes containing the node'
                        end if
                      end do
                      !
                      ! Whether the node belong to one or two symmetry planes (three is not allowed), the normal (esp) of the
                      ! first symmetry plane present in the list: x,  y, z; is used to build the direction
                      ! of the alpha rotation (v2 direction) and beta rotation (v1 direction) in such a way that:
                      !
                      !   - beta (v1 direction) is coincident with the esp vector: v1=esp; i.e. beta is the rotation with
                      !     direction normal with the symmetry plane.
                      !   - alpha (v2 direction) is v2 = v3 x esp; i.e. alpha is the rotation with direction tangent to the
                      !     symmetry plane.
                      !
                      esp=0.d0
                      esp(symplane_eid(node(sn)%symplane(1)))=1.d0
                      ! v1 = esp
                      v1=esp
                      v1=v1/sqrt(dot_product(v1,v1))
                      ! v2 = v3 x esp
                      v2=fbem_cross_product(element(ke)%v_midnode(:,3,kn),esp)
                      v2=v2/sqrt(dot_product(v2,v2))

!                      ! ¿?¿¿?¿?¿?
!                      ! Warning message for older files with regard to the local definition of rotations
!                      ! ¿?¿?¿?¿?¿?
!                      ! To be deprecated
!                      write(*,*) '-----'
!                      write(*,*) 'Node ', node(sn)%id
!                      write(*,*) 'Warning: if the file is old, check the boundary condition of this node, '
!                      write(*,*) 'it has been produced a change of criteria when selecting the local rotations.'
!                      if (node(sn)%n_symplanes.eq.1) write(*,*) 'Now, alpha rotation is the rotation with direction tangent to the symmetry plane, and beta with the normal.'
!                      if (node(sn)%n_symplanes.eq.2) write(*,*) 'Now, alpha rotation is the rotation with tangent direction, and beta with the normal, with respect to the first symmetry plane present in this list: x, y or z.'
!                      write(*,*) '-----'

                    !
                    ! If the node does not belong to any symmetry plane
                    !
                    else
                      ! Summatory of tbp-tbm of all elements at the selected node, and normalization
                      tbp=0.d0
                      do kej=1,node(sn)%n_elements
                        sej=node(sn)%element(kej)
                        tbp=tbp+element(sej)%tbp_gn(:,node(sn)%element_node_iid(kej))
                        tbp=tbp-element(sej)%tbm_gn(:,node(sn)%element_node_iid(kej))
                      end do
                      tbp=tbp/sqrt(dot_product(tbp,tbp))
                      ! v1 = tbp x v3
                      v1=fbem_cross_product(tbp,element(ke)%v_midnode(:,3,kn))
                      v1=v1/sqrt(dot_product(v1,v1))
                      ! v2 = v3 x v1
                      v2=fbem_cross_product(element(ke)%v_midnode(:,3,kn),v1)
                      v2=v2/sqrt(dot_product(v2,v2))
                    end if
                  !
                  ! The node is not at the boundary. A global criteria is followed in such a way that v1, v2 and v3 are the same
                  ! at the node for all the elements that contain it.
                  !
                  else
                    ! e1 as reference vector: v2 = v3 x e1
                    v2(1) = 0.0d0
                    v2(2) = element(ke)%v_midnode(3,3,kn)
                    v2(3) =-element(ke)%v_midnode(2,3,kn)
                    norm=dsqrt(dot_product(v2,v2))
                    if (norm.le.geometric_tolerance) then
                      ! e2 as reference vector: v2 = v3 x e2
                      v2(1) =-element(ke)%v_midnode(3,3,kn)
                      v2(2) = 0.0d0
                      v2(3) = element(ke)%v_midnode(1,3,kn)
                      norm=dsqrt(dot_product(v2,v2))
                      if (norm.le.geometric_tolerance) then
                        ! e3 as reference vector: v2 = v3 x e3
                        v2(1) = element(ke)%v_midnode(2,3,kn)
                        v2(2) =-element(ke)%v_midnode(1,3,kn)
                        v2(3) = 0.0d0
                        norm=dsqrt(dot_product(v2,v2))
                        if (norm.le.geometric_tolerance) then
                          stop 'Fatal error. Check subroutine build_data_at_functional_nodes() for FE.'
                        end if
                      end if
                    end if
                    v2=v2/norm
                    ! v1 = v2 x v3
                    v1=fbem_cross_product(v2,element(ke)%v_midnode(:,3,kn))
                    v1=v1/sqrt(dot_product(v1,v1))
                  end if
                end if
                ! Save results
                element(ke)%v_midnode(:,1,kn)=v1
                element(ke)%v_midnode(:,2,kn)=v2

                ! Write v1, v2, v3
!                write(22,'(2i11,6e25.16)') element(ke)%id, node(element(ke)%node(kn))%id, element(ke)%x_gn(:,kn), element(ke)%v_midnode(:,1,kn)
!                write(23,'(2i11,6e25.16)') element(ke)%id, node(element(ke)%node(kn))%id, element(ke)%x_gn(:,kn), element(ke)%v_midnode(:,2,kn)
!                write(24,'(2i11,6e25.16)') element(ke)%id, node(element(ke)%node(kn))%id, element(ke)%x_gn(:,kn), element(ke)%v_midnode(:,3,kn)

              end do

            ! --------------
            ! SOLID ELEMENTS
            ! --------------

            case (3)
              ! Nothing to do

          end select
      end select
    end if
  end do ! Loop through FINITE ELEMENTS

  ! Ending message
  if (verbose_level.ge.1) write(output_unit,'(a)') 'done.'

end subroutine build_data_at_functional_nodes
