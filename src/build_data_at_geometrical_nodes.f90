! ---------------------------------------------------------------------
! Copyright (C) 2014-2022 Universidad de Las Palmas de Gran Canaria:
!                         Jacob D.R. Bordon
!                         Guillermo M. Alamo
!                         Luis A. Padron
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
!! <b> Subroutine that builds the data at geometrical nodes. </b>

subroutine build_data_at_geometrical_nodes

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_data_structures
  use fbem_geometry
  use fbem_fem_beams
  use fbem_fem_shells

  ! Module of problem variables
  use problem_variables

  ! Local variables
  implicit none
  integer                        :: kb, sp, ke, se, kn, sn, kc, kej, sej, sp1, sp2
  integer                        :: sse, sse_fe, pi, pj, pk
  integer                        :: sn_fe, se_fe, kn_fe
  integer                        :: tmp_n_nodes, tmp_n_dimension
  real(kind=real64)              :: xi_1d, xi_2d(2), xi_3d(3), tmp
  real(kind=real64)              :: v1(problem%n), v2(problem%n), v3(problem%n)
  real(kind=real64)              :: t1(problem%n), t2(problem%n), t3(problem%n)
  real(kind=real64), allocatable :: x_nodes(:,:)
  real(kind=real64), allocatable :: x_gn(:)
  real(kind=real64), allocatable :: n_gn(:)
  real(kind=real64), allocatable :: t1_gn(:)
  real(kind=real64), allocatable :: t2_gn(:)
  real(kind=real64), allocatable :: tbp_gn(:)
  real(kind=real64), allocatable :: tbm_gn(:)

  ! Starting message
  if (verbose_level.ge.1)  write(output_unit,'(a)') 'Building element data at geometrical nodes ...'

  ! Allocate
  allocate (x_gn(problem%n),n_gn(problem%n),t1_gn(problem%n),t2_gn(problem%n),tbp_gn(problem%n),tbm_gn(problem%n))

  ! ================================================================================================================================
  ! GENERAL DATA AT GEOMETRICAL NODES OF THE ELEMENTS
  ! ================================================================================================================================

  ! Loop through the elements
  do ke=1,n_elements

    ! Tip/edge coupled BE body loads are treated later
    if (part(element(ke)%part)%type.eq.fbem_part_be_bodyload) then
      if (be_bodyload(part(element(ke)%part)%entity)%coupling.ne.0) then
        cycle
      end if
    end if


    ! Aqui habra problemas con las cargas puntuales no acopladas...?¿?¿

    !
    ! xi_gn
    !
    if (element(ke)%n_dimension.gt.0) then
      allocate (element(ke)%xi_gn(element(ke)%n_dimension,element(ke)%n_nodes))
      element(ke)%xi_gn=fbem_xi_hybrid(element(ke)%type,0.0d0)
    end if
    !
    ! x_gn
    !
    allocate (x_nodes(problem%n,element(ke)%n_nodes))
    do kn=1,element(ke)%n_nodes
      sn=element(ke)%node(kn)
      do kc=1,problem%n
        x_nodes(kc,kn)=node(sn)%x(kc)
      end do
    end do
    allocate (element(ke)%x_gn(problem%n,element(ke)%n_nodes))
    do kn=1,element(ke)%n_nodes
      do kc=1,problem%n
        element(ke)%x_gn(kc,kn)=x_nodes(kc,kn)
      end do
    end do
    !
    ! jacobian_gn
    !
    if (element(ke)%n_dimension.gt.0) then
      allocate (element(ke)%jacobian_gn(element(ke)%n_nodes))
      select case (problem%n)
        case (2)
          select case (element(ke)%n_dimension)
            case (1)
              do kn=1,element(ke)%n_nodes
                xi_1d=element(ke)%xi_gn(1,kn)
                element(ke)%jacobian_gn(kn)=fbem_jacobian2d(element(ke)%type,x_nodes,xi_1d)
              end do
            case (2)
              do kn=1,element(ke)%n_nodes
                xi_2d(1)=element(ke)%xi_gn(1,kn)
                xi_2d(2)=element(ke)%xi_gn(2,kn)
                element(ke)%jacobian_gn(kn)=fbem_jacobian2d(element(ke)%type,x_nodes,xi_2d)
              end do
          end select
        case (3)
          select case (element(ke)%n_dimension)
            case (1)
              do kn=1,element(ke)%n_nodes
                xi_1d=element(ke)%xi_gn(1,kn)
                element(ke)%jacobian_gn(kn)=fbem_jacobian3d(element(ke)%type,x_nodes,xi_1d)
              end do
            case (2)
              do kn=1,element(ke)%n_nodes
                xi_2d(1)=element(ke)%xi_gn(1,kn)
                xi_2d(2)=element(ke)%xi_gn(2,kn)
                element(ke)%jacobian_gn(kn)=fbem_jacobian3d(element(ke)%type,x_nodes,xi_2d)
              end do
          end select
      end select
    end if
    !
    ! t1_gn and t2_gn
    !
    if (element(ke)%n_dimension.gt.0) then
      select case (element(ke)%n_dimension)
        ! 1D elements
        case (1)
          allocate (element(ke)%t1_gn(problem%n,element(ke)%n_nodes))
          do kn=1,element(ke)%n_nodes
            xi_1d=element(ke)%xi_gn(1,kn)
            t1_gn=fbem_utangent_xi(problem%n,element(ke)%type,x_nodes,xi_1d)
            do kc=1,problem%n
              element(ke)%t1_gn(kc,kn)=t1_gn(kc)
            end do
          end do
        ! 2D elements
        case (2)
          allocate (element(ke)%t1_gn(problem%n,element(ke)%n_nodes))
          allocate (element(ke)%t2_gn(problem%n,element(ke)%n_nodes))
          do kn=1,element(ke)%n_nodes
            xi_2d(1)=element(ke)%xi_gn(1,kn)
            xi_2d(2)=element(ke)%xi_gn(2,kn)
            t1_gn=fbem_utangent_xi1(problem%n,element(ke)%type,x_nodes,xi_2d)
            t2_gn=fbem_utangent_xi2(problem%n,element(ke)%type,x_nodes,xi_2d)
            do kc=1,problem%n
              element(ke)%t1_gn(kc,kn)=t1_gn(kc)
              element(ke)%t2_gn(kc,kn)=t2_gn(kc)
            end do
          end do
      end select
    end if
    !
    ! n_gn
    !
    ! It is only valid done for 1D elements in 2D and 2D elements in 3D, otherwise n_gn is a zero vector.
    !
    if (element(ke)%n_dimension.gt.0) then
      allocate (element(ke)%n_gn(problem%n,element(ke)%n_nodes))
      do kn=1,element(ke)%n_nodes
        do kc=1,problem%n
          element(ke)%n_gn(kc,kn)=0.0d0
        end do
      end do
      ! 1D element in 2D
      if ((element(ke)%n_dimension.eq.1).and.(problem%n.eq.2)) then
        do kn=1,element(ke)%n_nodes
          xi_1d=element(ke)%xi_gn(1,kn)
          n_gn=fbem_unormal2d(element(ke)%type,x_nodes,xi_1d)
          do kc=1,2
            element(ke)%n_gn(kc,kn)=n_gn(kc)
          end do
        end do
      end if
      ! 2D elements in 3D
      if ((element(ke)%n_dimension.eq.2).and.(problem%n.eq.3)) then
        do kn=1,element(ke)%n_nodes
          xi_2d(1)=element(ke)%xi_gn(1,kn)
          xi_2d(2)=element(ke)%xi_gn(2,kn)
          n_gn=fbem_unormal3d(element(ke)%type,x_nodes,xi_2d)
          do kc=1,3
            element(ke)%n_gn(kc,kn)=n_gn(kc)
          end do
        end do
      end if
    end if
    !
    ! tbp_gn and tbm_gn
    !
    if (element(ke)%n_dimension.gt.0) then
      allocate (element(ke)%tbp_gn(problem%n,element(ke)%n_nodes))
      allocate (element(ke)%tbm_gn(problem%n,element(ke)%n_nodes))
      do kn=1,element(ke)%n_nodes
        call fbem_utangents_at_boundary(problem%n,element(ke)%type,x_nodes,kn,tbp_gn,tbm_gn)
        do kc=1,problem%n
          element(ke)%tbp_gn(kc,kn)=tbp_gn(kc)
          element(ke)%tbm_gn(kc,kn)=tbm_gn(kc)
        end do
      end do
    end if

    ! Deallocate
    deallocate (x_nodes)
  end do ! Loop through the elements

  ! ================================================================================================================================
  ! DIRECTOR VECTORS (v_1,v_2,v_3) OF DEGENERATED BEAMS AND (v_3) OF SHELL DEGENERATED FINITE ELEMENTS
  ! ================================================================================================================================

  ! Build: element()%v_midnode, element()%tv_midnode

  ! Loop through FINITE ELEMENTS
  do ke=1,n_elements

  ! aqui nos saltamos esto para las regiones rigidas, pero igual es necesario para cuando queremos introducir la masa/inercia
  ! en estas regiones si son vigas o laminas "rigidas".
  ! ademas, si queremos elementos rigidos de union viga-lamina, cada elemento rigido debe tener su master node

    if (part(element(ke)%part)%type.eq.fbem_part_fe_subregion) then

      ! no es necesario para regiones rigidas
      if (fe_subregion(part(element(ke)%part)%entity)%master_node.ne.0) cycle
      if (element(ke)%fe_type.ne.0) cycle

      select case (problem%n)

        ! ==
        ! 2D
        ! ==

        case (2)
          select case (element(ke)%n_dimension)

            ! ---------------
            ! POINT ELEMENTS
            ! ---------------

            case (0)
              ! Nothing to do

            ! -------------
            ! BEAM ELEMENTS
            ! -------------

            case (1)

              if (allocated(element(ke)%tn_midnode).eqv.(.false.)) then
                call fbem_error_message(error_unit,0,'element',element(ke)%id,'thicknesses of this beam element are undefined.')
              end if
              if (allocated(element(ke)%v_midnode).eqv.(.false.)) then
                call fbem_error_message(error_unit,0,'element',element(ke)%id,'director vector of this beam element is undefined.')
              end if
              allocate (element(ke)%tv_midnode(3,element(ke)%n_nodes))
              element(ke)%tv_midnode=0.d0
              do kn=1,element(ke)%n_nodes
                sn=element(ke)%node(kn)
                !
                ! Calculate v1, v2 and v3
                !
                ! For simplicity, it is assumed that v1 == t1 (cross section is always orthogonal to the mid-line).                
                v1=element(ke)%t1_gn(:,kn)
                ! Save
                element(ke)%v_midnode(1,1,kn)=v1(1)
                element(ke)%v_midnode(2,1,kn)=v1(2)
                element(ke)%v_midnode(3,1,kn)=0.
                ! v2 = e3 x v1
                element(ke)%v_midnode(1,2,kn)=-v1(2)
                element(ke)%v_midnode(2,2,kn)= v1(1)
                element(ke)%v_midnode(3,2,kn)= 0.d0
                ! v3 == e3
                element(ke)%v_midnode(1,3,kn)=0.d0
                element(ke)%v_midnode(2,3,kn)=0.d0
                element(ke)%v_midnode(3,3,kn)=1.d0
                !
                ! Transform thicknesses with respect to orthogonal directions into thicknesses with respect to v2 and v3
                !
                element(ke)%tv_midnode(2,kn)=element(ke)%tn_midnode(2,kn)/dot_product(v1,element(ke)%t1_gn(:,kn))
                element(ke)%tv_midnode(3,kn)=element(ke)%tn_midnode(3,kn)
                
!                ! Write v1, v2, v3
!                write(22,'(2i11,6e25.16)') element(ke)%id, node(element(ke)%node(kn))%id, element(ke)%x_gn(:,kn), element(ke)%v_midnode(:,1,kn)
!                write(23,'(2i11,6e25.16)') element(ke)%id, node(element(ke)%node(kn))%id, element(ke)%x_gn(:,kn), element(ke)%v_midnode(:,2,kn)
!                write(24,'(2i11,6e25.16)') element(ke)%id, node(element(ke)%node(kn))%id, element(ke)%x_gn(:,kn), element(ke)%v_midnode(:,3,kn)

              end do

            ! --------------
            ! SOLID ELEMENTS
            ! --------------

            case (2)
              select case (problem%subtype)
                case (fbem_mechanics_plane_strain)
                  if (allocated(element(ke)%tn_midnode)) then
                    call fbem_error_message(error_unit,0,__FILE__,__LINE__,'fatal error, this should not happen')
                  else
                    allocate (element(ke)%tn_midnode(3,element(ke)%n_nodes))
                    element(ke)%tn_midnode=1
                  end if
                case (fbem_mechanics_plane_stress)
                  if (.not.allocated(element(ke)%tn_midnode)) then
                    call fbem_error_message(error_unit,0,'element',element(ke)%id,'for plane stress, you must enter thickness at [cross sections]')
                  end if
              end select

          end select

        ! ==
        ! 3D
        ! ==

        case (3)
          select case (element(ke)%n_dimension)

            ! ---------------
            ! POINT ELEMENTS
            ! ---------------

            case (0)
              ! Nothing to do

            ! -------------
            ! BEAM ELEMENTS
            ! -------------

            case (1)
              if (allocated(element(ke)%tn_midnode).eqv.(.false.)) then
                call fbem_error_message(error_unit,0,'element',element(ke)%id,'thicknesses of this beam element are undefined.')
              end if
              if (allocated(element(ke)%v_midnode).eqv.(.false.)) then
                call fbem_error_message(error_unit,0,'element',element(ke)%id,'director vector of this beam element is undefined.')
              end if
              allocate (element(ke)%tv_midnode(3,element(ke)%n_nodes))
              element(ke)%tv_midnode=0.d0
              do kn=1,element(ke)%n_nodes
                sn=element(ke)%node(kn)
                !
                ! Calculate v1, v2 and v3
                !
                ! For simplicity, it is assumed that v1 == t1 (cross section is always orthogonal to the mid-line).
                v1=element(ke)%t1_gn(:,kn)
                ! Save
                element(ke)%v_midnode(:,1,kn)=v1
                ! Copy provisional v2
                v2=element(ke)%v_midnode(:,2,kn)
                ! v3 = (v1 x v2^{indicated})/|v1 x v2^{indicated}|
                v3=fbem_cross_product(v1,v2)
                v3=v3/sqrt(dot_product(v3,v3))
                ! v2 = v3 x v1
                v2=fbem_cross_product(v3,v1)
                v2=v2/sqrt(dot_product(v2,v2))
                ! Save
                element(ke)%v_midnode(:,2,kn)=v2
                element(ke)%v_midnode(:,3,kn)=v3
                !
                ! Transform thicknesses with respect to orthogonal directions into thicknesses with respect to v2 and v3
                !
                ! Build t1, t2 and t3, which are orthogonal to t1
                t1=element(ke)%t1_gn(:,kn)
                t2=element(ke)%v_midnode(:,2,kn)
                t3=fbem_cross_product(t1,t2)
                t3=t3/sqrt(dot_product(t3,t3))
                t2=fbem_cross_product(t3,t1)
                t2=t2/sqrt(dot_product(t2,t2))
                ! Calculate
                element(ke)%tv_midnode(2,kn)=element(ke)%tn_midnode(2,kn)/dot_product(v2,t2)
                element(ke)%tv_midnode(3,kn)=element(ke)%tn_midnode(3,kn)/dot_product(v3,t3)

                ! Write v1, v2, v3
                ! gnuplot
                !set view equal xyz
                !splot "fort.22" u 3:4:5:6:7:8 w vect, "fort.23" u 3:4:5:6:7:8 w vect, "fort.24" u 3:4:5:6:7:8 w vect
                !write(22,'(2i11,6e25.16)') element(ke)%id, node(element(ke)%node(kn))%id, element(ke)%x_gn(:,kn), element(ke)%v_midnode(:,1,kn)
                !write(23,'(2i11,8e25.16)') element(ke)%id, node(element(ke)%node(kn))%id, element(ke)%x_gn(:,kn), element(ke)%v_midnode(:,2,kn), element(ke)%tv_midnode(2,kn)
                !write(24,'(2i11,8e25.16)') element(ke)%id, node(element(ke)%node(kn))%id, element(ke)%x_gn(:,kn), element(ke)%v_midnode(:,3,kn), element(ke)%tv_midnode(3,kn)
                
              end do

!              do pi=0,100
!                do pj=0,10
!                  do pk=0,10
!                    xi_3d(1)=-1.d0+2.d0*dble(pi)/dble(100)
!                    xi_3d(2)=-1.d0+2.d0*dble(pj)/dble(10)
!                    xi_3d(3)=-1.d0+2.d0*dble(pk)/dble(10)
!                    v1=fbem_fem_degbeam3d_x(element(ke)%type,element(ke)%x_gn,element(ke)%v_midnode,element(ke)%tv_midnode,xi_3d)
!                    write(63,'(3e25.16)') v1
!                  end do
!                end do
!              end do

            ! --------------
            ! SHELL ELEMENTS
            ! --------------

            case (2)
              if (.not.allocated(element(ke)%tn_midnode)) then
                call fbem_error_message(error_unit,0,'element',element(ke)%id,'the thickness of this shell element is undefined.')
              end if
              allocate (element(ke)%v_midnode(3,3,element(ke)%n_nodes),element(ke)%tv_midnode(3,element(ke)%n_nodes))
              element(ke)%v_midnode=0.d0
              element(ke)%tv_midnode=0.d0
              do kn=1,element(ke)%n_nodes
                sn=element(ke)%node(kn)
                v3=0.d0
                if (node(sn)%is_singular) then
                  ! En uniones con elementos inclinados, podría existir diferencias importantes entre las v3 de este
                  ! nodo singular y las de un nodo interior de la arista que lo contiene
                  v3=element(ke)%n_gn(:,kn)
                else
                  do kej=1,node(sn)%n_elements
                    sej=node(sn)%element(kej)
                    v3=v3+element(sej)%n_gn(:,node(sn)%element_node_iid(kej))
                  end do
                end if
                ! If the node belongs to any symmetry plane, it is necessary to include the v_3 vector of symmetrical elements.
                select case (node(sn)%n_symplanes)
                  ! If it belongs to 1 symmetry plane (it can happen in 2D and 3D)
                  case (1)
                    ! Symmetry plane of the node
                    sp1=node(sn)%symplane(1)
                    ! Reflect and add the normal with reversed orientation.
                    do kc=1,problem%n
                      v3(kc)=v3(kc)+symplane_m(kc,sp1)*v3(kc)
                    end do
                  ! If it belongs to 2 symmetry planes (it can happen only in 3D)
                  case (2)
                    ! Symmetry planes of the node
                    sp1=node(sn)%symplane(1)
                    sp2=node(sn)%symplane(2)
                    ! Reflect and add the normal with reversed orientation with respect to the first symmetry plane.
                    do kc=1,problem%n
                      v3(kc)=v3(kc)+symplane_m(kc,sp1)*v3(kc)
                    end do
                    ! Reflect and add the normal with reversed orientation with respect to the second symmetry plane.
                    do kc=1,problem%n
                      v3(kc)=v3(kc)+symplane_m(kc,sp2)*v3(kc)
                    end do
                  case (3)
                    stop 'weird situation, a shell node belongs to three symmetry planes'
                end select
                if (sqrt(dot_product(v3,v3)).le.geometric_tolerance) then
                  call fbem_error_message(error_unit,0,'element',element(ke)%id,'a v3 vector of this element is close to zero. This shell is joined with other shell with an angle very close to 0º or 180º,  or it is coplanar with respect to a symmetry plane')
                end if
                v3=v3/sqrt(dot_product(v3,v3))
                element(ke)%v_midnode(:,3,kn)=v3
                element(ke)%tv_midnode(3,kn)=element(ke)%tn_midnode(3,kn)/dot_product(v3,element(ke)%n_gn(:,kn))
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

  ! ================================================================================================================================
  ! ELEMENT LOCAL AXIS (ep_1,ep_2,ep_3) FOR STRAIGHT BEAM FINITE ELEMENTS AND DEGENERATED SHELL FINITE ELEMENTS
  ! ================================================================================================================================

  ! Build: element()%ep

  ! Loop through FINITE ELEMENTS
  do ke=1,n_elements

    if (part(element(ke)%part)%type.eq.fbem_part_fe_subregion) then

      if (fe_subregion(part(element(ke)%part)%entity)%master_node.ne.0) cycle

      select case (element(ke)%n_dimension)

        ! ==========================================================================================================================
        ! ZERO-DIMENSIONAL ELEMENTS
        ! ==========================================================================================================================

        case (0)
          ! Nothing to do

        ! ==========================================================================================================================
        ! ONE-DIMENSIONAL ELEMENTS
        ! ==========================================================================================================================

        case (1)

          select case (element(ke)%fe_type)

            !-----------------------------------------------------------------------------------------------------------------------
            ! DEGENERATED BEAM FINITE ELEMENT
            !
            case (0)

              if (.not.allocated(element(ke)%ep)) then
                allocate(element(ke)%ep(3,3))
                element(ke)%ep=0
              end if

            !-----------------------------------------------------------------------------------------------------------------------

            !-----------------------------------------------------------------------------------------------------------------------
            ! STRAIGHT EULER-BERNOULLI AND TIMOSHENKO BEAM FINITE ELEMENTS
            !
            case (1,2)

              select case (problem%n)
                case (2)
                  element(ke)%ep(:,1)= element(ke)%t1_gn(:,1)
                  element(ke)%ep(1,2)=-element(ke)%ep(2,1)
                  element(ke)%ep(2,2)= element(ke)%ep(1,1)
                case (3)
                  element(ke)%ep(:,1)=element(ke)%t1_gn(:,1)
                  if (fbem_vector_norm(3,element(ke)%ep(:,2)).eq.0.d0) then
                    call fbem_error_message(error_unit,0,'element',element(ke)%id,'strbeam local axis ep2 is not correctly defined')
                  end if
                  element(ke)%ep(:,3)=fbem_cross_product(element(ke)%ep(:,1),element(ke)%ep(:,2))
                  if (fbem_vector_norm(3,element(ke)%ep(:,3)).lt.1.d-12) then
                    call fbem_error_message(error_unit,0,'element',element(ke)%id,'strbeam local axes ambiguous: ep1 x ep2 is approx. null')
                  end if
                  element(ke)%ep(:,3)=fbem_vector_normalization(3,element(ke)%ep(:,3))
                  element(ke)%ep(:,2)=fbem_cross_product(element(ke)%ep(:,3),element(ke)%ep(:,1))
              end select

            !-----------------------------------------------------------------------------------------------------------------------

            !-----------------------------------------------------------------------------------------------------------------------
            ! BAR FINITE ELEMENT
            !
            case (3)
              ! Nothing to do
            !-----------------------------------------------------------------------------------------------------------------------

            !-----------------------------------------------------------------------------------------------------------------------
            ! DISCRETE TRANSLATIONAL AND ROTATIONAL/TRANSLATIONAL SPRING FINITE ELEMENTS
            !
            case (4,5)

              if (abs(element(ke)%ep(1,1)-1.d0).lt.1.d-12) then
                do kc=1,problem%n
                  element(ke)%ep(kc,kc)=1
                end do
              else
                select case (problem%n)
                  case (2)
                    element(ke)%ep(:,1)= element(ke)%t1_gn(:,1)
                    element(ke)%ep(1,2)=-element(ke)%ep(2,1)
                    element(ke)%ep(2,2)= element(ke)%ep(1,1)
                  case (3)
                    element(ke)%ep(:,1)=element(ke)%t1_gn(:,1)
                    if (fbem_vector_norm(3,element(ke)%ep(:,2)).eq.0.d0) then
                      call fbem_error_message(error_unit,0,'element',element(ke)%id,'distra/disrotra local axis ep2 is not correctly defined')
                    end if
                    element(ke)%ep(:,3)=fbem_cross_product(element(ke)%ep(:,1),element(ke)%ep(:,2))
                    if (fbem_vector_norm(3,element(ke)%ep(:,3)).lt.1.d-12) then
                      call fbem_error_message(error_unit,0,'element',element(ke)%id,'distra/disrotra local axes ambiguous: ep1 x ep2 is approx. null')
                    end if
                    element(ke)%ep(:,3)=fbem_vector_normalization(3,element(ke)%ep(:,3))
                    element(ke)%ep(:,2)=fbem_cross_product(element(ke)%ep(:,3),element(ke)%ep(:,1))
                end select
              end if
            !-----------------------------------------------------------------------------------------------------------------------

            !-----------------------------------------------------------------------------------------------------------------------
            ! OTHER TYPES
            !
            case default
              !
              ! AQUI FALTA METER LOS DISROTRA
              !
              call fbem_error_message(error_unit,0,'element',element(ke)%id,'invalid type of 1D element')
            !-----------------------------------------------------------------------------------------------------------------------

          end select

        ! ==========================================================================================================================
        ! TWO-DIMENSIONAL ELEMENTS
        ! ==========================================================================================================================

        case (2)

          select case (problem%n)

            !-----------------------------------------------------------------------------------------------------------------------
            ! SOLID / CONTINUUM ELEMENTS
            !
            case (2)

              ! Nothing to do

            !-----------------------------------------------------------------------------------------------------------------------

            !-----------------------------------------------------------------------------------------------------------------------
            ! DEGENERATED SHELL FINITE ELEMENT
            !
            case (3)

              if (.not.allocated(element(ke)%ep)) then
                allocate(element(ke)%ep(3,3))
                element(ke)%ep=0
              end if

            !-----------------------------------------------------------------------------------------------------------------------

          end select

        ! ==========================================================================================================================
        ! THREE-DIMENSIONAL ELEMENTS
        ! ==========================================================================================================================

        case (3)

          ! Nothing to do

        !===========================================================================================================================
        ! OTHER XD
        !===========================================================================================================================

        case default

          call fbem_error_message(error_unit,0,'element',element(ke)%id,'invalid type dimension')

      end select

    end if

  end do ! Loop through FINITE ELEMENTS

  ! ================================================================================================================================
  ! GEOMETRICAL ELEMENT ASSOCIATED WITH THE COUPLED BE BODY LOADS FOR BE-FE INTERACTION
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
              ! Geometrical description of the line load
              element(se)%type_g=fbem_line2
              tmp_n_nodes=fbem_n_nodes(element(se)%type_g)
              tmp_n_dimension=fbem_n_dimension(element(se)%type_g)
              allocate (element(se)%x_gn(problem%n,tmp_n_nodes))
              allocate (element(se)%xi_gn(tmp_n_dimension,tmp_n_nodes))
              ! xi coordinates of the geometrical nodes
              element(se)%xi_gn=fbem_xi_hybrid(element(se)%type_g,0.0d0)
              ! Thickness vector (normal)
              sn_fe=node(element(se)%node(1))%coupled_node
              se_fe=node(sn_fe)%element(1)
              kn_fe=node(sn_fe)%element_node_iid(1)
              v2=-element(se_fe)%t1_gn(2,kn_fe)*0.5d0*element(se_fe)%A
              v2= element(se_fe)%t1_gn(1,kn_fe)*0.5d0*element(se_fe)%A
              ! Position of the nodes
              element(se)%x_gn(:,1)=node(sn_fe)%x+v2
              element(se)%x_gn(:,2)=node(sn_fe)%x-v2
            end do
          !
          ! 3D: SURFACE LOAD ATH THE BEAM TIP
          !
          case (3)

            ! Carga en punta
            ! Modelo de Luis (en principio solo habria que crear un cuadrilatero que emule a una circunferencia, se puede hacer..)
            stop 'not yet 47'


        end select

      ! ----------------------
      ! FE BEAM - BE LINE LOAD
      ! ----------------------

      ! Modelo de Luis

      case (fbem_bl_coupling_beam_line)

        do ke=1,part(sp)%n_elements
          se=part(sp)%element(ke)

          ! lo normal para todos los elementos

          !
          ! xi_gn
          !
          allocate (element(se)%xi_gn(element(se)%n_dimension,element(se)%n_nodes))
          element(se)%xi_gn=fbem_xi_hybrid(element(se)%type,0.0d0)
          !
          ! x_gn
          !
          allocate (x_nodes(problem%n,element(se)%n_nodes))
          do kn=1,element(se)%n_nodes
            sn=element(se)%node(kn)
            do kc=1,problem%n
              x_nodes(kc,kn)=node(sn)%x(kc)
            end do
          end do
          allocate (element(se)%x_gn(problem%n,element(se)%n_nodes))
          do kn=1,element(se)%n_nodes
            do kc=1,problem%n
              element(se)%x_gn(kc,kn)=x_nodes(kc,kn)
            end do
          end do
          !
          ! jacobian_gn
          !
          allocate (element(se)%jacobian_gn(element(se)%n_nodes))
          select case (problem%n)
            case (2)
              select case (element(se)%n_dimension)
                case (1)
                  do kn=1,element(se)%n_nodes
                    xi_1d=element(se)%xi_gn(1,kn)
                    element(se)%jacobian_gn(kn)=fbem_jacobian2d(element(se)%type,x_nodes,xi_1d)
                  end do
                case (2)
                  do kn=1,element(se)%n_nodes
                    xi_2d(1)=element(se)%xi_gn(1,kn)
                    xi_2d(2)=element(se)%xi_gn(2,kn)
                    element(se)%jacobian_gn(kn)=fbem_jacobian2d(element(se)%type,x_nodes,xi_2d)
                  end do
              end select
            case (3)
              select case (element(se)%n_dimension)
                case (1)
                  do kn=1,element(se)%n_nodes
                    xi_1d=element(se)%xi_gn(1,kn)
                    element(se)%jacobian_gn(kn)=fbem_jacobian3d(element(se)%type,x_nodes,xi_1d)
                  end do
                case (2)
                  do kn=1,element(se)%n_nodes
                    xi_2d(1)=element(se)%xi_gn(1,kn)
                    xi_2d(2)=element(se)%xi_gn(2,kn)
                    element(se)%jacobian_gn(kn)=fbem_jacobian3d(element(se)%type,x_nodes,xi_2d)
                  end do
              end select
          end select
          !
          ! t1_gn and t2_gn
          !
          select case (element(se)%n_dimension)
            ! 1D elements
            case (1)
              allocate (element(se)%t1_gn(problem%n,element(se)%n_nodes))
              do kn=1,element(se)%n_nodes
                xi_1d=element(se)%xi_gn(1,kn)
                t1_gn=fbem_utangent_xi(problem%n,element(se)%type,x_nodes,xi_1d)
                do kc=1,problem%n
                  element(se)%t1_gn(kc,kn)=t1_gn(kc)
                end do
              end do
            ! 2D elements
            case (2)
              allocate (element(se)%t1_gn(problem%n,element(se)%n_nodes))
              allocate (element(se)%t2_gn(problem%n,element(se)%n_nodes))
              do kn=1,element(se)%n_nodes
                xi_2d(1)=element(se)%xi_gn(1,kn)
                xi_2d(2)=element(se)%xi_gn(2,kn)
                t1_gn=fbem_utangent_xi1(problem%n,element(se)%type,x_nodes,xi_2d)
                t2_gn=fbem_utangent_xi2(problem%n,element(se)%type,x_nodes,xi_2d)
                do kc=1,problem%n
                  element(se)%t1_gn(kc,kn)=t1_gn(kc)
                  element(se)%t2_gn(kc,kn)=t2_gn(kc)
                end do
              end do
          end select
          !
          ! n_gn
          !
          ! It is only valid done for 1D elements in 2D and 2D elements in 3D, otherwise n_gn is a zero vector.
          !
          allocate (element(se)%n_gn(problem%n,element(se)%n_nodes))
          do kn=1,element(se)%n_nodes
            do kc=1,problem%n
              element(se)%n_gn(kc,kn)=0.0d0
            end do
          end do
          ! 1D element in 2D
          if ((element(se)%n_dimension.eq.1).and.(problem%n.eq.2)) then
            do kn=1,element(se)%n_nodes
              xi_1d=element(se)%xi_gn(1,kn)
              n_gn=fbem_unormal2d(element(se)%type,x_nodes,xi_1d)
              do kc=1,2
                element(se)%n_gn(kc,kn)=n_gn(kc)
              end do
            end do
          end if
          ! 2D elements in 3D
          if ((element(se)%n_dimension.eq.2).and.(problem%n.eq.3)) then
            do kn=1,element(se)%n_nodes
              xi_2d(1)=element(se)%xi_gn(1,kn)
              xi_2d(2)=element(se)%xi_gn(2,kn)
              n_gn=fbem_unormal3d(element(se)%type,x_nodes,xi_2d)
              do kc=1,3
                element(se)%n_gn(kc,kn)=n_gn(kc)
              end do
            end do
          end if
          !
          ! tbp_gn and tbm_gn
          !
          allocate (element(se)%tbp_gn(problem%n,element(se)%n_nodes))
          allocate (element(se)%tbm_gn(problem%n,element(se)%n_nodes))
          do kn=1,element(se)%n_nodes
            call fbem_utangents_at_boundary(problem%n,element(se)%type,x_nodes,kn,tbp_gn,tbm_gn)
            do kc=1,problem%n
              element(se)%tbp_gn(kc,kn)=tbp_gn(kc)
              element(se)%tbm_gn(kc,kn)=tbm_gn(kc)
            end do
          end do
          ! Deallocate
          deallocate (x_nodes)

        end do

      ! -----------------------
      ! FE SHELL - BE EDGE LOAD
      ! -----------------------

      case (fbem_bl_coupling_shell_edge)
        do ke=1,part(sp)%n_elements
          se=part(sp)%element(ke)
          ! Geometrical description of the surface load
          select case (element(se)%type)
            case (fbem_line2)
              element(se)%type_g=fbem_quad4
            case (fbem_line3)
              element(se)%type_g=fbem_quad9
          end select
          tmp_n_nodes=fbem_n_nodes(element(se)%type_g)
          tmp_n_dimension=fbem_n_dimension(element(se)%type_g)
          allocate (element(se)%x_gn(problem%n,tmp_n_nodes))
          allocate (element(se)%xi_gn(tmp_n_dimension,tmp_n_nodes))
          ! xi coordinates of the geometrical nodes
          element(se)%xi_gn=fbem_xi_hybrid(element(se)%type_g,0.0d0)
          ! Thickness vector for each node is taken from the Shell FE
          sse=element(se)%subedge(1)
          sse_fe=subedge(sse)%element
          se_fe=subedge(sse_fe)%supelement(1)
          do kn=1,subedge(sse)%n_nodes
            sn=subedge(sse)%node(kn)
            sn_fe=subedge(sse)%element_node(kn)
            do kn_fe=1,element(se_fe)%n_nodes
              if (element(se_fe)%node(kn_fe).eq.sn_fe) then
                v3=element(se_fe)%v_midnode(:,3,kn_fe)*element(se_fe)%tv_midnode(3,kn_fe)
                exit
              end if
            end do
            select case (element(se)%type)
              case (fbem_line2)
                ! QUAD4 SURFACE LOAD GEOMETRY
                select case (kn)
                  case (1)
                    element(se)%x_gn(:,4)=node(sn_fe)%x+0.5d0*v3
                    element(se)%x_gn(:,1)=node(sn_fe)%x-0.5d0*v3
                  case (2)
                    element(se)%x_gn(:,3)=node(sn_fe)%x+0.5d0*v3
                    element(se)%x_gn(:,2)=node(sn_fe)%x-0.5d0*v3
                end select
              case (fbem_line3)
                ! QUAD9 SURFACE LOAD GEOMETRY
                select case (kn)
                  case (1)
                    element(se)%x_gn(:,4)=node(sn_fe)%x+0.5d0*v3
                    element(se)%x_gn(:,8)=node(sn_fe)%x
                    element(se)%x_gn(:,1)=node(sn_fe)%x-0.5d0*v3
                  case (2)
                    element(se)%x_gn(:,3)=node(sn_fe)%x+0.5d0*v3
                    element(se)%x_gn(:,6)=node(sn_fe)%x
                    element(se)%x_gn(:,2)=node(sn_fe)%x-0.5d0*v3
                  case (3)
                    element(se)%x_gn(:,7)=node(sn_fe)%x+0.5d0*v3
                    element(se)%x_gn(:,9)=node(sn_fe)%x
                    element(se)%x_gn(:,5)=node(sn_fe)%x-0.5d0*v3
                end select
            end select
          end do
!          ! Write
!          do kn=1,fbem_n_nodes(element(se)%type_g)
!            write(*,'(2i11,3e25.16)') element(se)%id, kn, element(se)%x_gn(:,kn)
!          end do
        end do

      ! --------------------------
      ! FE SHELL - BE SURFACE LOAD
      ! --------------------------

      case (fbem_bl_coupling_shell_surface)
        do ke=1,part(sp)%n_elements
          se=part(sp)%element(ke)

          ! lo normal para todos los elementos

          !
          ! xi_gn
          !
          allocate (element(se)%xi_gn(element(se)%n_dimension,element(se)%n_nodes))
          element(se)%xi_gn=fbem_xi_hybrid(element(se)%type,0.0d0)
          !
          ! x_gn
          !
          allocate (x_nodes(problem%n,element(se)%n_nodes))
          do kn=1,element(se)%n_nodes
            sn=element(se)%node(kn)
            do kc=1,problem%n
              x_nodes(kc,kn)=node(sn)%x(kc)
            end do
          end do
          allocate (element(se)%x_gn(problem%n,element(se)%n_nodes))
          do kn=1,element(se)%n_nodes
            do kc=1,problem%n
              element(se)%x_gn(kc,kn)=x_nodes(kc,kn)
            end do
          end do
          !
          ! jacobian_gn
          !
          allocate (element(se)%jacobian_gn(element(se)%n_nodes))
          select case (problem%n)
            case (2)
              select case (element(se)%n_dimension)
                case (1)
                  do kn=1,element(se)%n_nodes
                    xi_1d=element(se)%xi_gn(1,kn)
                    element(se)%jacobian_gn(kn)=fbem_jacobian2d(element(se)%type,x_nodes,xi_1d)
                  end do
                case (2)
                  do kn=1,element(se)%n_nodes
                    xi_2d(1)=element(se)%xi_gn(1,kn)
                    xi_2d(2)=element(se)%xi_gn(2,kn)
                    element(se)%jacobian_gn(kn)=fbem_jacobian2d(element(se)%type,x_nodes,xi_2d)
                  end do
              end select
            case (3)
              select case (element(se)%n_dimension)
                case (1)
                  do kn=1,element(se)%n_nodes
                    xi_1d=element(se)%xi_gn(1,kn)
                    element(se)%jacobian_gn(kn)=fbem_jacobian3d(element(se)%type,x_nodes,xi_1d)
                  end do
                case (2)
                  do kn=1,element(se)%n_nodes
                    xi_2d(1)=element(se)%xi_gn(1,kn)
                    xi_2d(2)=element(se)%xi_gn(2,kn)
                    element(se)%jacobian_gn(kn)=fbem_jacobian3d(element(se)%type,x_nodes,xi_2d)
                  end do
              end select
          end select
          !
          ! t1_gn and t2_gn
          !
          select case (element(se)%n_dimension)
            ! 1D elements
            case (1)
              allocate (element(se)%t1_gn(problem%n,element(se)%n_nodes))
              do kn=1,element(se)%n_nodes
                xi_1d=element(se)%xi_gn(1,kn)
                t1_gn=fbem_utangent_xi(problem%n,element(se)%type,x_nodes,xi_1d)
                do kc=1,problem%n
                  element(se)%t1_gn(kc,kn)=t1_gn(kc)
                end do
              end do
            ! 2D elements
            case (2)
              allocate (element(se)%t1_gn(problem%n,element(se)%n_nodes))
              allocate (element(se)%t2_gn(problem%n,element(se)%n_nodes))
              do kn=1,element(se)%n_nodes
                xi_2d(1)=element(se)%xi_gn(1,kn)
                xi_2d(2)=element(se)%xi_gn(2,kn)
                t1_gn=fbem_utangent_xi1(problem%n,element(se)%type,x_nodes,xi_2d)
                t2_gn=fbem_utangent_xi2(problem%n,element(se)%type,x_nodes,xi_2d)
                do kc=1,problem%n
                  element(se)%t1_gn(kc,kn)=t1_gn(kc)
                  element(se)%t2_gn(kc,kn)=t2_gn(kc)
                end do
              end do
          end select
          !
          ! n_gn
          !
          ! It is only valid done for 1D elements in 2D and 2D elements in 3D, otherwise n_gn is a zero vector.
          !
          allocate (element(se)%n_gn(problem%n,element(se)%n_nodes))
          do kn=1,element(se)%n_nodes
            do kc=1,problem%n
              element(se)%n_gn(kc,kn)=0.0d0
            end do
          end do
          ! 1D element in 2D
          if ((element(se)%n_dimension.eq.1).and.(problem%n.eq.2)) then
            do kn=1,element(se)%n_nodes
              xi_1d=element(se)%xi_gn(1,kn)
              n_gn=fbem_unormal2d(element(se)%type,x_nodes,xi_1d)
              do kc=1,2
                element(se)%n_gn(kc,kn)=n_gn(kc)
              end do
            end do
          end if
          ! 2D elements in 3D
          if ((element(se)%n_dimension.eq.2).and.(problem%n.eq.3)) then
            do kn=1,element(se)%n_nodes
              xi_2d(1)=element(se)%xi_gn(1,kn)
              xi_2d(2)=element(se)%xi_gn(2,kn)
              n_gn=fbem_unormal3d(element(se)%type,x_nodes,xi_2d)
              do kc=1,3
                element(se)%n_gn(kc,kn)=n_gn(kc)
              end do
            end do
          end if
          !
          ! tbp_gn and tbm_gn
          !
          allocate (element(se)%tbp_gn(problem%n,element(se)%n_nodes))
          allocate (element(se)%tbm_gn(problem%n,element(se)%n_nodes))
          do kn=1,element(se)%n_nodes
            call fbem_utangents_at_boundary(problem%n,element(se)%type,x_nodes,kn,tbp_gn,tbm_gn)
            do kc=1,problem%n
              element(se)%tbp_gn(kc,kn)=tbp_gn(kc)
              element(se)%tbm_gn(kc,kn)=tbm_gn(kc)
            end do
          end do
          ! Deallocate
          deallocate (x_nodes)

        end do

    end select
  end do ! Loop through COUPLED BE BODY LOAD ELEMENTS

  ! Ending message
  if (verbose_level.ge.1) write(output_unit,'(a)') 'done.'

end subroutine build_data_at_geometrical_nodes
