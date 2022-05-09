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
!! <b> Subroutine that builds element data at collocation points. </b>

subroutine build_data_at_collocation_points

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_string_handling
  use fbem_numerical
  use fbem_shape_functions
  use fbem_data_structures
  use fbem_geometry

  ! Module of problem variables
  use problem_variables

  ! Local variables
  implicit none
  integer                        :: kb, sp, ke, kn, kc, kip, ka
  integer                        :: se, sn, kn1, kn2, dme, se_fe
  integer                        :: tmp_n_gnodes, tmp_n_dimension, tmp_f2nodes
  real(kind=real64)              :: x_i(problem%n), n_i(problem%n), t1_i(problem%n), t2_i(problem%n)
  real(kind=real64)              :: ep1(3), ep2(3), ep3(3), r_integration, A, Iy, Iz
  real(kind=real64), allocatable :: xi(:)
  real(kind=real64)              :: xi_1d
  real(kind=real64)              :: xi_2d(2)
  real(kind=real64), allocatable :: x_nodes(:,:), dme_x_nodes(:,:)
  real(kind=real64), allocatable :: dxda_nodes(:,:,:)
  real(kind=real64), allocatable :: xi_nodes(:,:)
  real(kind=real64), allocatable :: phi(:), dphidx(:,:)
  real(kind=real64)              :: d
  integer                        :: tmp_unit

  ! Starting message
  if (verbose_level.ge.1)  write(output_unit,'(a)') 'Building element data at collocation points ...'

  ! ====================
  ! BE BOUNDARY ELEMENTS
  ! ====================

  ! Loop through BE BOUNDARIES
  do kb=1,n_boundaries
    ! Loop through the ELEMENTS of the BE BOUNDARY
    do ke=1,part(boundary(kb)%part)%n_elements
      ! Selected element
      se=part(boundary(kb)%part)%element(ke)
      ! Copy the coordinates of the geometrical nodes
      allocate (x_nodes(problem%n,element(se)%n_nodes))
      x_nodes=element(se)%x_gn
      ! Copy the reference coordinates of the geometrical nodes
      allocate (xi_nodes(element(se)%n_dimension,element(se)%n_nodes))
      xi_nodes=element(se)%xi_gn
      !
      ! CALCULATE XI_I FOR THE FORMULATIONS OF THE NODES
      !
      ! If the element is continuous
      if (element(se)%discontinuous.eqv.(.false.)) then
        ! Switch between element dimensions
        select case (element(se)%n_dimension)
          ! 1D
          case (1)
            ! Allocate "xi_i" members
            allocate (element(se)%xi_i_sbie(1,element(se)%n_nodes))
            allocate (element(se)%xi_i_sbie_mca(1,element(se)%n_nodes))
            allocate (element(se)%xi_i_hbie(1,element(se)%n_nodes))
            ! Loop through the nodes of the element
            do kn=1,element(se)%n_nodes
              ! Selected node
              sn=element(se)%node(kn)
              ! If the node has SBIE, then copy directly the xi coordinate of the geometrical node to the "xi_i_sbie" member
              if (node(sn)%sbie.eq.fbem_sbie) then
                element(se)%xi_i_sbie(1,kn)=xi_nodes(1,kn)
              end if
              ! If the node has SBIE MCA, then move towards inside the element the xi coordinate of the geometrical node,
              ! and copy it to the "xi_i_sbie_mca" member
              if (node(sn)%sbie.eq.fbem_sbie_mca) then
                ! Save xi coordinate of the geometrical node
                xi_1d=xi_nodes(1,kn)
                ! Move the xi coordinate
                xi_1d=fbem_move_xi_from_vertex(xi_1d,element(se)%delta_sbie_mca(kn))
                ! Copy it
                element(se)%xi_i_sbie_mca(1,kn)=xi_1d
              end if
              ! If the node has HBIE formulation, then move towards inside the element the xi coordinate of
              ! the geometrical node, and copy it to the "xi_i_hbie" member
              if (node(sn)%hbie.eq.fbem_hbie) then
                ! Save xi coordinate of the geometrical node
                xi_1d=xi_nodes(1,kn)
                ! Move the xi coordinate
                xi_1d=fbem_move_xi_from_vertex(xi_1d,element(se)%delta_hbie(kn))
                ! Copy it
                element(se)%xi_i_hbie(1,kn)=xi_1d
              end if
            end do
          ! 2D
          case (2)
            ! Allocate "xi_i" member
            allocate (element(se)%xi_i_sbie(2,element(se)%n_nodes))
            allocate (element(se)%xi_i_sbie_mca(2,element(se)%n_nodes))
            allocate (element(se)%xi_i_hbie(2,element(se)%n_nodes))
            ! Loop through the nodes of the element
            do kn=1,element(se)%n_nodes
              ! Selected node
              sn=element(se)%node(kn)
              ! If the node has SBIE or SBIE + SBIE MCA formulation, then copy directly the xi coordinate of the geometrical node
              ! to the "xi_i_sbie" member
              if (node(sn)%sbie.eq.fbem_sbie) then
                element(se)%xi_i_sbie(1,kn)=xi_nodes(1,kn)
                element(se)%xi_i_sbie(2,kn)=xi_nodes(2,kn)
              end if
              ! If the node has SBIE MCA, then move towards inside the element the xi coordinate of the geometrical node,
              ! and copy it to the "xi_i_sbie_mca" member
              if (node(sn)%sbie.eq.fbem_sbie_mca) then
                ! Save xi coordinate of the geometrical node
                xi_2d(1)=xi_nodes(1,kn)
                xi_2d(2)=xi_nodes(2,kn)
                ! Move the xi coordinate
                xi_2d=fbem_move_xi1xi2_from_edge(element(se)%type,xi_2d,element(se)%delta_sbie_mca(kn))
                ! Copy it
                element(se)%xi_i_sbie_mca(1,kn)=xi_2d(1)
                element(se)%xi_i_sbie_mca(2,kn)=xi_2d(2)
              end if
              ! If the node has HBIE formulation, then move towards inside the element the xi coordinate of
              ! the geometrical node, and copy it to the "xi_i_hbie" member
              if (node(sn)%hbie.eq.fbem_hbie) then
                ! Save xi coordinate of the geometrical node
                xi_2d(1)=xi_nodes(1,kn)
                xi_2d(2)=xi_nodes(2,kn)
                ! Move the xi coordinate
                xi_2d=fbem_move_xi1xi2_from_edge(element(se)%type,xi_2d,element(se)%delta_hbie(kn))
                ! Copy it
                element(se)%xi_i_hbie(1,kn)=xi_2d(1)
                element(se)%xi_i_hbie(2,kn)=xi_2d(2)
              end if
            end do
        end select
      ! If the element is discontinuous
      else
        allocate (element(se)%xi_i_sbie(element(se)%n_dimension,element(se)%n_nodes))
        element(se)%xi_i_sbie=fbem_xi_hybrid(element(se)%type,element(se)%delta_f)
        allocate (element(se)%xi_i_hbie(element(se)%n_dimension,element(se)%n_nodes))
        element(se)%xi_i_hbie=fbem_xi_hybrid(element(se)%type,element(se)%delta_f)
      end if
      !
      ! Calculate x_i
      !
      ! Allocate "x_i" members
      allocate (element(se)%x_i_sbie(problem%n,element(se)%n_nodes))
      allocate (element(se)%x_i_sbie_mca(problem%n,element(se)%n_nodes))
      allocate (element(se)%x_i_hbie(problem%n,element(se)%n_nodes))
      ! Switch element dimensions
      select case (element(se)%n_dimension)
        ! 1D
        case (1)
          ! Switch problem dimensions
          select case (problem%n)
            ! 2D
            case (2)
              do kn=1,element(se)%n_nodes
                ! Selected node
                sn=element(se)%node(kn)
                ! If the node has SBIE formulation
                if (node(sn)%sbie.eq.fbem_sbie) then
                  xi_1d=element(se)%xi_i_sbie(1,kn)
                  x_i=fbem_position2d(element(se)%type,x_nodes,xi_1d)
                  do kc=1,2
                    element(se)%x_i_sbie(kc,kn)=x_i(kc)
                  end do
                end if
                ! If the node has SBIE MCA formulation
                if (node(sn)%sbie.eq.fbem_sbie_mca) then
                  xi_1d=element(se)%xi_i_sbie_mca(1,kn)
                  x_i=fbem_position2d(element(se)%type,x_nodes,xi_1d)
                  do kc=1,2
                    element(se)%x_i_sbie_mca(kc,kn)=x_i(kc)
                  end do
                end if
                ! If the node has HBIE formulation
                if (node(sn)%hbie.eq.fbem_hbie) then
                  xi_1d=element(se)%xi_i_hbie(1,kn)
                  x_i=fbem_position2d(element(se)%type,x_nodes,xi_1d)
                  do kc=1,2
                    element(se)%x_i_hbie(kc,kn)=x_i(kc)
                  end do
                end if
              end do
            ! 3D
            case (3)
              do kn=1,element(se)%n_nodes
                ! Selected node
                sn=element(se)%node(kn)
                ! If the node has SBIE formulation
                if (node(sn)%sbie.eq.fbem_sbie) then
                  xi_1d=element(se)%xi_i_sbie(1,kn)
                  x_i=fbem_position3d(element(se)%type,x_nodes,xi_1d)
                  do kc=1,3
                    element(se)%x_i_sbie(kc,kn)=x_i(kc)
                  end do
                end if
                ! If the node has SBIE MCA formulation
                if (node(sn)%sbie.eq.fbem_sbie_mca) then
                  xi_1d=element(se)%xi_i_sbie_mca(1,kn)
                  x_i=fbem_position3d(element(se)%type,x_nodes,xi_1d)
                  do kc=1,3
                    element(se)%x_i_sbie_mca(kc,kn)=x_i(kc)
                  end do
                end if
                ! If the node has HBIE formulation
                if (node(sn)%hbie.eq.fbem_hbie) then
                  xi_1d=element(se)%xi_i_hbie(1,kn)
                  x_i=fbem_position3d(element(se)%type,x_nodes,xi_1d)
                  do kc=1,3
                    element(se)%x_i_hbie(kc,kn)=x_i(kc)
                  end do
                end if
              end do
          end select
        ! 2D
        case (2)
          ! Switch problem dimensions
          select case (problem%n)
            ! 2D
            case (2)
              do kn=1,element(se)%n_nodes
                ! Selected node
                sn=element(se)%node(kn)
                ! If the node has SBIE formulation
                if (node(sn)%sbie.eq.fbem_sbie) then
                  xi_2d(1)=element(se)%xi_i_sbie(1,kn)
                  xi_2d(2)=element(se)%xi_i_sbie(2,kn)
                  x_i=fbem_position2d(element(se)%type,x_nodes,xi_2d)
                  do kc=1,2
                    element(se)%x_i_sbie(kc,kn)=x_i(kc)
                  end do
                end if
                ! If the node has SBIE MCA formulation
                if (node(sn)%sbie.eq.fbem_sbie_mca) then
                  xi_2d(1)=element(se)%xi_i_sbie_mca(1,kn)
                  xi_2d(2)=element(se)%xi_i_sbie_mca(2,kn)
                  x_i=fbem_position2d(element(se)%type,x_nodes,xi_2d)
                  do kc=1,2
                    element(se)%x_i_sbie_mca(kc,kn)=x_i(kc)
                  end do
                end if
                ! If the node has HBIE formulation
                if (node(sn)%hbie.eq.fbem_hbie) then
                  xi_2d(1)=element(se)%xi_i_hbie(1,kn)
                  xi_2d(2)=element(se)%xi_i_hbie(2,kn)
                  x_i=fbem_position2d(element(se)%type,x_nodes,xi_2d)
                  do kc=1,2
                    element(se)%x_i_hbie(kc,kn)=x_i(kc)
                  end do
                end if
              end do
            ! 3D
            case (3)
              do kn=1,element(se)%n_nodes
                ! Selected node
                sn=element(se)%node(kn)
                ! If the node has SBIE formulation
                if (node(sn)%sbie.eq.fbem_sbie) then
                  xi_2d(1)=element(se)%xi_i_sbie(1,kn)
                  xi_2d(2)=element(se)%xi_i_sbie(2,kn)
                  x_i=fbem_position3d(element(se)%type,x_nodes,xi_2d)
                  do kc=1,3
                    element(se)%x_i_sbie(kc,kn)=x_i(kc)
                  end do
                end if
                ! If the node has SBIE MCA formulation
                if (node(sn)%sbie.eq.fbem_sbie_mca) then
                  xi_2d(1)=element(se)%xi_i_sbie_mca(1,kn)
                  xi_2d(2)=element(se)%xi_i_sbie_mca(2,kn)
                  x_i=fbem_position3d(element(se)%type,x_nodes,xi_2d)
                  do kc=1,3
                    element(se)%x_i_sbie_mca(kc,kn)=x_i(kc)
                  end do
                end if
                ! If the node has HBIE formulation
                if (node(sn)%hbie.eq.fbem_hbie) then
                  xi_2d(1)=element(se)%xi_i_hbie(1,kn)
                  xi_2d(2)=element(se)%xi_i_hbie(2,kn)
                  x_i=fbem_position3d(element(se)%type,x_nodes,xi_2d)
                  do kc=1,3
                    element(se)%x_i_hbie(kc,kn)=x_i(kc)
                  end do
                end if
              end do
          end select
      end select
      !
      ! Calculate n_i
      !
      ! It is only valid done for 1D elements in 2D and 2D elements in 3D, otherwise a warning message is displayed.
      !
      ! Allocate "n_i" member
      allocate (element(se)%n_i_sbie(problem%n,element(se)%n_nodes))
      allocate (element(se)%n_i_sbie_mca(problem%n,element(se)%n_nodes))
      allocate (element(se)%n_i_hbie(problem%n,element(se)%n_nodes))
      ! 1D element in 2D
      if ((element(se)%n_dimension.eq.1).and.(problem%n.eq.2)) then
        do kn=1,element(se)%n_nodes
          ! Selected node
          sn=element(se)%node(kn)
          ! If the node has SBIE formulation
          if (node(sn)%sbie.eq.fbem_sbie) then
            xi_1d=element(se)%xi_i_sbie(1,kn)
            n_i=fbem_unormal2d(element(se)%type,x_nodes,xi_1d)
            do kc=1,2
              element(se)%n_i_sbie(kc,kn)=n_i(kc)
            end do
          end if
          ! If the node has SBIE MCA formulation
          if (node(sn)%sbie.eq.fbem_sbie_mca) then
            xi_1d=element(se)%xi_i_sbie_mca(1,kn)
            n_i=fbem_unormal2d(element(se)%type,x_nodes,xi_1d)
            do kc=1,2
              element(se)%n_i_sbie_mca(kc,kn)=n_i(kc)
            end do
          end if
          ! If the node has HBIE formulation
          if (node(sn)%hbie.eq.fbem_hbie) then
            xi_1d=element(se)%xi_i_hbie(1,kn)
            n_i=fbem_unormal2d(element(se)%type,x_nodes,xi_1d)
            do kc=1,2
              element(se)%n_i_hbie(kc,kn)=n_i(kc)
            end do
          end if
        end do
      end if
      ! 2D elements in 3D
      if ((element(se)%n_dimension.eq.2).and.(problem%n.eq.3)) then
        do kn=1,element(se)%n_nodes
          ! Selected node
          sn=element(se)%node(kn)
          ! If the node has SBIE formulation
          if (node(sn)%sbie.eq.fbem_sbie) then
            xi_2d(1)=element(se)%xi_i_sbie(1,kn)
            xi_2d(2)=element(se)%xi_i_sbie(2,kn)
            n_i=fbem_unormal3d(element(se)%type,x_nodes,xi_2d)
            do kc=1,3
              element(se)%n_i_sbie(kc,kn)=n_i(kc)
            end do
          end if
          ! If the node has SBIE MCA formulation
          if (node(sn)%sbie.eq.fbem_sbie_mca) then
            xi_2d(1)=element(se)%xi_i_sbie_mca(1,kn)
            xi_2d(2)=element(se)%xi_i_sbie_mca(2,kn)
            n_i=fbem_unormal3d(element(se)%type,x_nodes,xi_2d)
            do kc=1,3
              element(se)%n_i_sbie_mca(kc,kn)=n_i(kc)
            end do
          end if
          ! If the node has HBIE formulation
          if (node(sn)%hbie.eq.fbem_hbie) then
            xi_2d(1)=element(se)%xi_i_hbie(1,kn)
            xi_2d(2)=element(se)%xi_i_hbie(2,kn)
            n_i=fbem_unormal3d(element(se)%type,x_nodes,xi_2d)
            do kc=1,3
              element(se)%n_i_hbie(kc,kn)=n_i(kc)
            end do
          end if
        end do
      end if
      !
      ! Calculate xi of the design element and dxda (for sensitivity analysis) at collocation points
      !
      if (problem%sensitivity) then
        !
        ! Isoparametric parametrization
        !
        if (element(se)%dm_mode.eq.0) then
          ! Allocate data structures
          allocate (element(se)%dm_xi_i_sbie(element(se)%n_dimension,element(se)%n_nodes))
          allocate (element(se)%dm_xi_i_sbie_mca(element(se)%n_dimension,element(se)%n_nodes))
          allocate (element(se)%dm_xi_i_hbie(element(se)%n_dimension,element(se)%n_nodes))
          allocate (element(se)%dxda_i_sbie(problem%n,element(se)%n_nodes,problem%n_designvariables))
          allocate (element(se)%dxda_i_sbie_mca(problem%n,element(se)%n_nodes,problem%n_designvariables))
          allocate (element(se)%dxda_i_hbie(problem%n,element(se)%n_nodes,problem%n_designvariables))
          allocate (element(se)%d2xdadx_i_sbie(problem%n,problem%n,element(se)%n_nodes,problem%n_designvariables))
          allocate (element(se)%d2xdadx_i_sbie_mca(problem%n,problem%n,element(se)%n_nodes,problem%n_designvariables))
          allocate (element(se)%d2xdadx_i_hbie(problem%n,problem%n,element(se)%n_nodes,problem%n_designvariables))
          allocate (element(se)%dnda_i_sbie(problem%n,element(se)%n_nodes,problem%n_designvariables))
          allocate (element(se)%dnda_i_sbie_mca(problem%n,element(se)%n_nodes,problem%n_designvariables))
          allocate (element(se)%dnda_i_hbie(problem%n,element(se)%n_nodes,problem%n_designvariables))
          allocate (dxda_nodes(problem%n,element(se)%n_nodes,problem%n_designvariables))
          allocate (phi(element(se)%n_nodes))
          allocate (dphidx(element(se)%n_nodes,problem%n))
          ! xi_i
          element(se)%dm_xi_i_sbie=element(se)%xi_i_sbie
          element(se)%dm_xi_i_sbie_mca=element(se)%xi_i_sbie_mca
          element(se)%dm_xi_i_hbie=element(se)%xi_i_hbie
          ! dxda of the nodes of the element
          do kn=1,element(se)%n_nodes
            sn=element(se)%node(kn)
            dxda_nodes(:,kn,:)=node(sn)%dxda
          end do
          ! dxda_i
          do kn=1,element(se)%n_nodes
            sn=element(se)%node(kn)
            ! If the node has SBIE formulation
            if (node(sn)%sbie.eq.fbem_sbie) then
              phi=fbem_phi_hybrid(element(se)%type_g,0.d0,element(se)%dm_xi_i_sbie(:,kn))
              element(se)%dxda_i_sbie(:,kn,:)=0.d0
              do kn2=1,element(se)%n_nodes
                element(se)%dxda_i_sbie(:,kn,:)=element(se)%dxda_i_sbie(:,kn,:)+phi(kn2)*dxda_nodes(:,kn2,:)
              end do
            end if
            ! If the node has SBIE MCA formulation
            if (node(sn)%sbie.eq.fbem_sbie_mca) then
              phi=fbem_phi_hybrid(element(se)%type_g,0.d0,element(se)%dm_xi_i_sbie_mca(:,kn))
              element(se)%dxda_i_sbie_mca(:,kn,:)=0.d0
              do kn2=1,element(se)%n_nodes
                element(se)%dxda_i_sbie_mca(:,kn,:)=element(se)%dxda_i_sbie_mca(:,kn,:)+phi(kn2)*dxda_nodes(:,kn2,:)
              end do
            end if
            ! If the node has HBIE formulation
            if (node(sn)%hbie.eq.fbem_hbie) then
              phi=fbem_phi_hybrid(element(se)%type_g,0.d0,element(se)%dm_xi_i_hbie(:,kn))
              element(se)%dxda_i_hbie(:,kn,:)=0.d0
              do kn2=1,element(se)%n_nodes
                element(se)%dxda_i_hbie(:,kn,:)=element(se)%dxda_i_hbie(:,kn,:)+phi(kn2)*dxda_nodes(:,kn2,:)
              end do
            end if
          end do
          ! d2xdadx_i
          do kn=1,element(se)%n_nodes
            sn=element(se)%node(kn)
            ! If the node has SBIE formulation
            if (node(sn)%sbie.eq.fbem_sbie) then
              dphidx=fbem_dphidx(problem%n,element(se)%type_g,0.d0,x_nodes,element(se)%dm_xi_i_sbie(:,kn))
              element(se)%d2xdadx_i_sbie(:,:,kn,:)=0.d0
              do kc=1,problem%n
                do kn2=1,element(se)%n_nodes
                  element(se)%d2xdadx_i_sbie(:,kc,kn,:)=element(se)%d2xdadx_i_sbie(:,kc,kn,:)+dphidx(kn2,kc)*dxda_nodes(:,kn2,:)
                end do
              end do
            end if
            ! If the node has SBIE MCA formulation
            if (node(sn)%sbie.eq.fbem_sbie_mca) then
              dphidx=fbem_dphidx(problem%n,element(se)%type_g,0.d0,x_nodes,element(se)%dm_xi_i_sbie_mca(:,kn))
              element(se)%d2xdadx_i_sbie_mca(:,:,kn,:)=0.d0
              do kc=1,problem%n
                do kn2=1,element(se)%n_nodes
                  element(se)%d2xdadx_i_sbie_mca(:,kc,kn,:)=element(se)%d2xdadx_i_sbie_mca(:,kc,kn,:)+dphidx(kn2,kc)*dxda_nodes(:,kn2,:)
                end do
              end do
            end if
            ! If the node has HBIE formulation
            if (node(sn)%hbie.eq.fbem_hbie) then
              dphidx=fbem_dphidx(problem%n,element(se)%type_g,0.d0,x_nodes,element(se)%dm_xi_i_hbie(:,kn))
              element(se)%d2xdadx_i_hbie(:,:,kn,:)=0.d0
              do kc=1,problem%n
                do kn2=1,element(se)%n_nodes
                  element(se)%d2xdadx_i_hbie(:,kc,kn,:)=element(se)%d2xdadx_i_hbie(:,kc,kn,:)+dphidx(kn2,kc)*dxda_nodes(:,kn2,:)
                end do
              end do
            end if
          end do
          ! dnda_i
          if (problem%n.ne.2) stop 'dnda_i: only for 2d'
          do kn=1,element(se)%n_nodes
            sn=element(se)%node(kn)
            ! If the node has SBIE formulation
            if (node(sn)%sbie.eq.fbem_sbie) then
              n_i=element(se)%n_i_sbie(:,kn)
              t1_i(1)=-n_i(2)
              t1_i(2)=n_i(1)
              do ka=1,problem%n_designvariables
                d=n_i(1)*dot_product(element(se)%d2xdadx_i_sbie(1,:,kn,ka),t1_i)&
                 +n_i(2)*dot_product(element(se)%d2xdadx_i_sbie(2,:,kn,ka),t1_i)
                element(se)%dnda_i_sbie(:,kn,ka)=-t1_i(:)*d
              end do
            end if
            ! If the node has SBIE MCA formulation
            if (node(sn)%sbie.eq.fbem_sbie_mca) then
              n_i=element(se)%n_i_sbie_mca(:,kn)
              t1_i(1)=-n_i(2)
              t1_i(2)=n_i(1)
              do ka=1,problem%n_designvariables
                d=n_i(1)*dot_product(element(se)%d2xdadx_i_sbie_mca(1,:,kn,ka),t1_i)&
                 +n_i(2)*dot_product(element(se)%d2xdadx_i_sbie_mca(2,:,kn,ka),t1_i)
                element(se)%dnda_i_sbie_mca(:,kn,ka)=-t1_i(:)*d
              end do
            end if
            ! If the node has HBIE formulation
            if (node(sn)%hbie.eq.fbem_hbie) then
              n_i=element(se)%n_i_hbie(:,kn)
              t1_i(1)=-n_i(2)
              t1_i(2)=n_i(1)
              do ka=1,problem%n_designvariables
                d=n_i(1)*dot_product(element(se)%d2xdadx_i_hbie(1,:,kn,ka),t1_i)&
                 +n_i(2)*dot_product(element(se)%d2xdadx_i_hbie(2,:,kn,ka),t1_i)
                element(se)%dnda_i_hbie(:,kn,ka)=-t1_i(:)*d
              end do
            end if
          end do
          ! Deallocate
          deallocate (dxda_nodes,phi,dphidx)
        !
        ! Macro-Element parametrization
        !
        else
          ! Design Macro-Element
          dme=element(se)%dm_element(1)
          ! Allocate data structures
          allocate (element(se)%dm_xi_i_sbie(design_mesh%element(dme)%n_dimension,element(se)%n_nodes))
          allocate (element(se)%dm_xi_i_sbie_mca(design_mesh%element(dme)%n_dimension,element(se)%n_nodes))
          allocate (element(se)%dm_xi_i_hbie(design_mesh%element(dme)%n_dimension,element(se)%n_nodes))
          allocate (element(se)%dxda_i_sbie(problem%n,element(se)%n_nodes,problem%n_designvariables))
          allocate (element(se)%dxda_i_sbie_mca(problem%n,element(se)%n_nodes,problem%n_designvariables))
          allocate (element(se)%dxda_i_hbie(problem%n,element(se)%n_nodes,problem%n_designvariables))
          allocate (element(se)%d2xdadx_i_sbie(problem%n,problem%n,element(se)%n_nodes,problem%n_designvariables))
          allocate (element(se)%d2xdadx_i_sbie_mca(problem%n,problem%n,element(se)%n_nodes,problem%n_designvariables))
          allocate (element(se)%d2xdadx_i_hbie(problem%n,problem%n,element(se)%n_nodes,problem%n_designvariables))
          allocate (element(se)%dnda_i_sbie(problem%n,element(se)%n_nodes,problem%n_designvariables))
          allocate (element(se)%dnda_i_sbie_mca(problem%n,element(se)%n_nodes,problem%n_designvariables))
          allocate (element(se)%dnda_i_hbie(problem%n,element(se)%n_nodes,problem%n_designvariables))
          allocate (dxda_nodes(problem%n,design_mesh%element(dme)%n_nodes,problem%n_designvariables))
          allocate (dme_x_nodes(problem%n,design_mesh%element(dme)%n_nodes))
          allocate (xi(design_mesh%element(dme)%n_dimension))
          allocate (phi(design_mesh%element(dme)%n_nodes))
          allocate (dphidx(design_mesh%element(dme)%n_nodes,problem%n))
          ! x and dxda of the design element
          do kn=1,design_mesh%element(dme)%n_nodes
            sn=design_mesh%element(dme)%node(kn)
            dxda_nodes(:,kn,:)=design_mesh%node(sn)%dxda
            dme_x_nodes(:,kn)=design_mesh%node(sn)%x
          end do
          ! xi_i
          do kn=1,element(se)%n_nodes
            sn=element(se)%node(kn)
            ! If the node has SBIE formulation
            if (node(sn)%sbie.eq.fbem_sbie) then
              x_i=element(se)%x_i_sbie(:,kn)
              call fbem_local_coordinates(problem%n,design_mesh%element(dme)%type,design_mesh%element(dme)%x_gn,&
                                          design_mesh%element(dme)%csize,x_i,xi,d)
              if (d.gt.1.d-12) then
                call fbem_error_message(error_unit,0,'element',element(se)%id,&
                                        'a collocation point is not completely embedded in its design macro-element.')
              end if
              element(se)%dm_xi_i_sbie(:,kn)=xi
            end if
            ! If the node has SBIE MCA formulation
            if (node(sn)%sbie.eq.fbem_sbie_mca) then
              x_i=element(se)%x_i_sbie_mca(:,kn)
              call fbem_local_coordinates(problem%n,design_mesh%element(dme)%type,design_mesh%element(dme)%x_gn,&
                                          design_mesh%element(dme)%csize,x_i,xi,d)
              if (d.gt.1.d-12) then
                call fbem_error_message(error_unit,0,'element',element(se)%id,&
                                        'a collocation point is not completely embedded in its design macro-element.')
              end if
              element(se)%dm_xi_i_sbie_mca(:,kn)=xi
            end if
            ! If the node has HBIE formulation
            if (node(sn)%hbie.eq.fbem_hbie) then
              x_i=element(se)%x_i_hbie(:,kn)
              call fbem_local_coordinates(problem%n,design_mesh%element(dme)%type,design_mesh%element(dme)%x_gn,&
                                          design_mesh%element(dme)%csize,x_i,xi,d)
              if (d.gt.1.d-12) then
                call fbem_error_message(error_unit,0,'element',element(se)%id,&
                                        'a collocation point is not completely embedded in its design macro-element.')
              end if
              element(se)%dm_xi_i_hbie(:,kn)=xi
            end if
          end do
          ! dxda_i
          do kn=1,element(se)%n_nodes
            sn=element(se)%node(kn)
            ! If the node has SBIE formulation
            if (node(sn)%sbie.eq.fbem_sbie) then
              phi=fbem_phi_hybrid(design_mesh%element(dme)%type,0.d0,element(se)%dm_xi_i_sbie(:,kn))
              element(se)%dxda_i_sbie(:,kn,:)=0.d0
              do kn2=1,design_mesh%element(dme)%n_nodes
                element(se)%dxda_i_sbie(:,kn,:)=element(se)%dxda_i_sbie(:,kn,:)+phi(kn2)*dxda_nodes(:,kn2,:)
              end do
            end if
            ! If the node has SBIE MCA formulation
            if (node(sn)%sbie.eq.fbem_sbie_mca) then
              phi=fbem_phi_hybrid(design_mesh%element(dme)%type,0.d0,element(se)%dm_xi_i_sbie_mca(:,kn))
              element(se)%dxda_i_sbie_mca(:,kn,:)=0.d0
              do kn2=1,design_mesh%element(dme)%n_nodes
                element(se)%dxda_i_sbie_mca(:,kn,:)=element(se)%dxda_i_sbie_mca(:,kn,:)+phi(kn2)*dxda_nodes(:,kn2,:)
              end do
            end if
            ! If the node has HBIE formulation
            if (node(sn)%hbie.eq.fbem_hbie) then
              phi=fbem_phi_hybrid(design_mesh%element(dme)%type,0.d0,element(se)%dm_xi_i_hbie(:,kn))
              element(se)%dxda_i_hbie(:,kn,:)=0.d0
              do kn2=1,design_mesh%element(dme)%n_nodes
                element(se)%dxda_i_hbie(:,kn,:)=element(se)%dxda_i_hbie(:,kn,:)+phi(kn2)*dxda_nodes(:,kn2,:)
              end do
            end if
          end do
          ! d2xdadx_i
          do kn=1,element(se)%n_nodes
            sn=element(se)%node(kn)
            ! If the node has SBIE formulation
            if (node(sn)%sbie.eq.fbem_sbie) then
              dphidx=fbem_dphidx(problem%n,design_mesh%element(dme)%type,0.d0,dme_x_nodes,element(se)%dm_xi_i_sbie(:,kn))
              element(se)%d2xdadx_i_sbie(:,:,kn,:)=0.d0
              do kc=1,problem%n
                do kn2=1,design_mesh%element(dme)%n_nodes
                  element(se)%d2xdadx_i_sbie(:,kc,kn,:)=element(se)%d2xdadx_i_sbie(:,kc,kn,:)+dphidx(kn2,kc)*dxda_nodes(:,kn2,:)
                end do
              end do
            end if
            ! If the node has SBIE MCA formulation
            if (node(sn)%sbie.eq.fbem_sbie_mca) then
              dphidx=fbem_dphidx(problem%n,design_mesh%element(dme)%type,0.d0,dme_x_nodes,element(se)%dm_xi_i_sbie_mca(:,kn))
              element(se)%d2xdadx_i_sbie_mca(:,:,kn,:)=0.d0
              do kc=1,problem%n
                do kn2=1,design_mesh%element(dme)%n_nodes
                  element(se)%d2xdadx_i_sbie_mca(:,kc,kn,:)=element(se)%d2xdadx_i_sbie_mca(:,kc,kn,:)+dphidx(kn2,kc)*dxda_nodes(:,kn2,:)
                end do
              end do
            end if
            ! If the node has HBIE formulation
            if (node(sn)%hbie.eq.fbem_hbie) then
              dphidx=fbem_dphidx(problem%n,design_mesh%element(dme)%type,0.d0,dme_x_nodes,element(se)%dm_xi_i_hbie(:,kn))
              element(se)%d2xdadx_i_hbie(:,:,kn,:)=0.d0
              do kc=1,problem%n
                do kn2=1,design_mesh%element(dme)%n_nodes
                  element(se)%d2xdadx_i_hbie(:,kc,kn,:)=element(se)%d2xdadx_i_hbie(:,kc,kn,:)+dphidx(kn2,kc)*dxda_nodes(:,kn2,:)
                end do
              end do
            end if
          end do
          ! dnda_i
          if (problem%n.ne.2) stop 'dnda_i: only for 2d'
          do kn=1,element(se)%n_nodes
            sn=element(se)%node(kn)
            ! If the node has SBIE formulation
            if (node(sn)%sbie.eq.fbem_sbie) then
              n_i=element(se)%n_i_sbie(:,kn)
              t1_i(1)=-n_i(2)
              t1_i(2)=n_i(1)
              do ka=1,problem%n_designvariables
                d=n_i(1)*dot_product(element(se)%d2xdadx_i_sbie(1,:,kn,ka),t1_i)&
                 +n_i(2)*dot_product(element(se)%d2xdadx_i_sbie(2,:,kn,ka),t1_i)
                element(se)%dnda_i_sbie(:,kn,ka)=-t1_i(:)*d
              end do
            end if
            ! If the node has SBIE MCA formulation
            if (node(sn)%sbie.eq.fbem_sbie_mca) then
              n_i=element(se)%n_i_sbie_mca(:,kn)
              t1_i(1)=-n_i(2)
              t1_i(2)=n_i(1)
              do ka=1,problem%n_designvariables
                d=n_i(1)*dot_product(element(se)%d2xdadx_i_sbie_mca(1,:,kn,ka),t1_i)&
                 +n_i(2)*dot_product(element(se)%d2xdadx_i_sbie_mca(2,:,kn,ka),t1_i)
                element(se)%dnda_i_sbie_mca(:,kn,ka)=-t1_i(:)*d
              end do
            end if
            ! If the node has HBIE formulation
            if (node(sn)%hbie.eq.fbem_hbie) then
              n_i=element(se)%n_i_hbie(:,kn)
              t1_i(1)=-n_i(2)
              t1_i(2)=n_i(1)
              do ka=1,problem%n_designvariables
                d=n_i(1)*dot_product(element(se)%d2xdadx_i_hbie(1,:,kn,ka),t1_i)&
                 +n_i(2)*dot_product(element(se)%d2xdadx_i_hbie(2,:,kn,ka),t1_i)
                element(se)%dnda_i_hbie(:,kn,ka)=-t1_i(:)*d
              end do
              !write(69,'(6e25.16)') element(se)%x_i_hbie(:,kn), element(se)%n_i_hbie(:,kn), element(se)%dnda_i_hbie(:,kn,1)
            end if
          end do
          ! Deallocate
          deallocate (dxda_nodes,dme_x_nodes,phi,dphidx,xi)
        end if
      end if
      ! Deallocate
      deallocate (x_nodes,xi_nodes)
    end do ! Loop through the ELEMENTS of the BE BOUNDARY
  end do ! Loop through BE BOUNDARIES

!  !
!  ! Export collocation points
!  !
!  tmp_unit=fbem_get_valid_unit()
!  open(unit=tmp_unit,file='x_i.dat',action='write',recl=fbem_file_record_length)
!  do ke=1,n_elements
!    do kn=1,element(ke)%n_nodes
!      sn=element(ke)%node(kn)
!      if (node(sn)%sbie.eq.fbem_sbie) then
!        select case (problem%n)
!          case (2)
!            write(tmp_unit,'(2i11,2e25.16)') boundary(part(element(ke)%part)%entity)%id, 1, element(ke)%x_i_sbie(1,kn), element(ke)%x_i_sbie(2,kn)
!          case (3)
!            write(tmp_unit,'(2i11,3e25.16)') boundary(part(element(ke)%part)%entity)%id, 1, element(ke)%x_i_sbie(1,kn), element(ke)%x_i_sbie(2,kn), element(ke)%x_i_sbie(3,kn)
!        end select
!      end if
!      if (node(sn)%sbie.eq.fbem_sbie_mca) then
!        select case (problem%n)
!          case (2)
!            write(tmp_unit,'(2i11,2e25.16)') boundary(part(element(ke)%part)%entity)%id, 2, element(ke)%x_i_sbie_mca(1,kn), element(ke)%x_i_sbie_mca(2,kn)
!          case (3)
!            write(tmp_unit,'(2i11,3e25.16)') boundary(part(element(ke)%part)%entity)%id, 2, element(ke)%x_i_sbie_mca(1,kn), element(ke)%x_i_sbie_mca(2,kn), element(ke)%x_i_sbie_mca(3,kn)
!        end select
!      end if
!      if (node(sn)%hbie.eq.fbem_hbie) then
!        select case (problem%n)
!          case (2)
!            write(tmp_unit,'(2i11,2e25.16)') boundary(part(element(ke)%part)%entity)%id, 3, element(ke)%x_i_hbie(1,kn), element(ke)%x_i_hbie(2,kn)
!          case (3)
!            write(tmp_unit,'(2i11,3e25.16)') boundary(part(element(ke)%part)%entity)%id, 3, element(ke)%x_i_hbie(1,kn), element(ke)%x_i_hbie(2,kn), element(ke)%x_i_hbie(3,kn)
!        end select
!      end if
!    end do
!  end do
!  close(tmp_unit)
!  stop

  ! =============================
  ! COUPLED BE BODY LOAD ELEMENTS
  ! =============================

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

              ! Remember that this elements have in reality 1 higher dimension than indicated in %n_dimension, one must look
              ! at type_g instead of type or type_f (which are constants by now...)

              tmp_n_gnodes=fbem_n_nodes(element(se)%type_g)
              tmp_n_dimension=fbem_n_dimension(element(se)%type_g)
              tmp_f2nodes=fbem_n_nodes(element(se)%type_f2)
              ! Copy the coordinates of the geometrical nodes
              allocate (x_nodes(problem%n,tmp_n_gnodes))
              x_nodes=element(se)%x_gn
              ! Copy the reference coordinates of the geometrical nodes
              allocate (xi_nodes(tmp_n_dimension,tmp_n_gnodes))
              xi_nodes=element(se)%xi_gn
              !
              ! CALCULATE XI_I FOR THE FORMULATIONS OF THE NODES
              !
              select case (element(se)%type_f2)
                case (fbem_point1)
                  allocate (element(se)%cbl_n_cp(1))
                  element(se)%cbl_n_cp(1)=1
                  allocate (element(se)%cbl_xi_i(1,1,1))
                  element(se)%cbl_xi_i(1,1,1)= 0.d0
                case (fbem_line2)
                  allocate (element(se)%cbl_n_cp(2))
                  element(se)%cbl_n_cp(1)=1
                  element(se)%cbl_n_cp(2)=1
                  allocate (element(se)%cbl_xi_i(1,1,2))
                  element(se)%cbl_xi_i(1,1,1)=-1.d0
                  element(se)%cbl_xi_i(1,1,2)= 1.d0
                case (fbem_line3)
                  allocate (element(se)%cbl_n_cp(3))
                  element(se)%cbl_n_cp(1)=1
                  element(se)%cbl_n_cp(2)=1
                  element(se)%cbl_n_cp(3)=1
                  allocate (element(se)%cbl_xi_i(1,1,3))
                  element(se)%cbl_xi_i(1,1,1)=-1.d0
                  element(se)%cbl_xi_i(1,1,2)= 1.d0
                  element(se)%cbl_xi_i(1,1,3)= 0.d0
                case (fbem_line4)
                  allocate (element(se)%cbl_n_cp(4))
                  element(se)%cbl_n_cp(1)=1
                  element(se)%cbl_n_cp(2)=1
                  element(se)%cbl_n_cp(3)=1
                  element(se)%cbl_n_cp(4)=1
                  allocate (element(se)%cbl_xi_i(1,1,4))
                  element(se)%cbl_xi_i(1,1,1)=-1.d0
                  element(se)%cbl_xi_i(1,1,2)= 1.d0
                  element(se)%cbl_xi_i(1,1,3)=-1.d0/3.d0
                  element(se)%cbl_xi_i(1,1,4)= 1.d0/3.d0
              end select
              !
              ! CALCULATE x^i
              !
              allocate (element(se)%cbl_x_i(problem%n,maxval(element(se)%cbl_n_cp),tmp_f2nodes))
              do kn=1,tmp_f2nodes
                do kc=1,element(se)%cbl_n_cp(kn)
                  xi_1d=element(se)%cbl_xi_i(1,kc,kn)
                  x_i=fbem_position2d(element(se)%type_g,x_nodes,xi_1d)
                  element(se)%cbl_x_i(:,kc,kn)=x_i
                end do
              end do

              ! Deallocate
              deallocate (x_nodes,xi_nodes)

            end do
          !
          ! 3D: SURFACE LOAD ATH THE BEAM TIP
          !
          case (3)
            ! Modelo de Luis (en principio solo habria que crear un cuadrilatero que emule a una circunferencia, se puede hacer..)
            stop 'not yet 46'
        end select

      ! -------------------------- !
      ! 3D: FE BEAM - BE LINE LOAD !
      ! -------------------------- !

      case (fbem_bl_coupling_beam_line)
        do ke=1,part(be_bodyload(kb)%part)%n_elements
          se=part(be_bodyload(kb)%part)%element(ke)
          !
          ! Save the radius of integration
          !
          se_fe=element(se)%element
          r_integration=element(se_fe)%r_integration
          !
          ! Save the basic data of the element
          !
          allocate (x_nodes(3,element(se)%n_nodes))
          x_nodes=element(se)%x_gn
          allocate (xi_nodes(1,element(se)%n_nodes))
          xi_nodes=element(se)%xi_gn
          !
          ! CALCULATE XI_I (IN THE ELEMENT) FOR THE FORMULATIONS OF THE NODES
          !
          if (element(se)%discontinuous) then
            stop 'Error: discontinuous elements not allowed for beam-line load coupling'
          else
            allocate (element(se)%xi_i_sbie(1,element(se)%n_nodes))
            allocate (element(se)%xi_i_sbie_mca(1,element(se)%n_nodes))
            allocate (element(se)%xi_i_hbie(1,element(se)%n_nodes))
            do kn=1,element(se)%n_nodes
              sn=element(se)%node(kn)
              if (node(sn)%sbie.eq.fbem_sbie) then
                element(se)%xi_i_sbie(1,kn)=xi_nodes(1,kn)
              end if
              if (node(sn)%sbie.eq.fbem_sbie_mca) then
                xi_1d=xi_nodes(1,kn)
                xi_1d=fbem_move_xi_from_vertex(xi_1d,element(se)%delta_sbie_mca(kn))
                element(se)%xi_i_sbie_mca(1,kn)=xi_1d
              end if
              if (node(sn)%hbie.eq.fbem_hbie) then
                xi_1d=xi_nodes(1,kn)
                xi_1d=fbem_move_xi_from_vertex(xi_1d,element(se)%delta_hbie(kn))
                element(se)%xi_i_hbie(1,kn)=xi_1d
              end if
            end do
          end if
          !
          ! CALCULATE X_I
          !
          allocate (element(se)%x_i_sbie    (3,element(se)%n_nodes))
          allocate (element(se)%x_i_sbie_mca(3,element(se)%n_nodes))
          allocate (element(se)%x_i_hbie    (3,element(se)%n_nodes))
          do kn=1,element(se)%n_nodes
            sn=element(se)%node(kn)
            if (node(sn)%sbie.eq.fbem_sbie) then
              xi_1d=element(se)%xi_i_sbie(1,kn)
              x_i=fbem_position3d(element(se)%type,x_nodes,xi_1d)
              do kc=1,3
                element(se)%x_i_sbie(kc,kn)=x_i(kc)
              end do
            end if
            if (node(sn)%sbie.eq.fbem_sbie_mca) then
              xi_1d=element(se)%xi_i_sbie_mca(1,kn)
              x_i=fbem_position3d(element(se)%type,x_nodes,xi_1d)
              do kc=1,3
                element(se)%x_i_sbie_mca(kc,kn)=x_i(kc)
              end do
            end if
          end do
          !
          ! CALCULATE X_I FOR FINAL COLLOCATION POINTS
          !
          allocate (element(se)%cbl_n_cp(element(se)%n_nodes))
          allocate (element(se)%cbl_x_i(3,4,element(se)%n_nodes))
          element(se)%cbl_n_cp=4
          do kn=1,element(se)%n_nodes
            sn=element(se)%node(kn)
            if (node(sn)%sbie.eq.fbem_sbie) then
              xi_1d=element(se)%xi_i_sbie(1,kn)
            end if
            if (node(sn)%sbie.eq.fbem_sbie_mca) then
              xi_1d=element(se)%xi_i_sbie_mca(1,kn)
            end if
            ep1=fbem_utangent_xi(3,element(se)%type_g,element(se)%x_gn,xi_1d)
            ! e1 as reference vector: ep2 = ep1 x e1
            ep2(1) = 0.0d0
            ep2(2) = ep1(3)
            ep2(3) =-ep1(2)
            if (sqrt(dot_product(ep2,ep2)).le.1.d-14) then
              ! e2 as reference vector: ep2 = ep1 x e2
              ep2(1) =-ep1(3)
              ep2(2) = 0.0d0
              ep2(3) = ep1(1)
              if (sqrt(dot_product(ep2,ep2)).le.1.d-14) then
                ! e3 as reference vector: ep2 = ep1 x e3
                ep2(1) = ep1(2)
                ep2(2) =-ep1(1)
                ep2(3) = 0.0d0
                if (sqrt(dot_product(ep2,ep2)).le.1.d-14) then
                  stop 'Fatal error 666'
                end if
              end if
            end if
            ep2=ep2/dsqrt(dot_product(ep2,ep2))
            ! ep3 = ep1 x ep2
            ep3=fbem_cross_product(ep1,ep2)
            ep3=ep3/sqrt(dot_product(ep3,ep3))
            ! Save position of collocation points
            if (node(sn)%sbie.eq.fbem_sbie) then
              element(se)%cbl_x_i(:,1,kn)=element(se)%x_i_sbie(:,kn)+ep2*r_integration
              element(se)%cbl_x_i(:,2,kn)=element(se)%x_i_sbie(:,kn)-ep2*r_integration
              element(se)%cbl_x_i(:,3,kn)=element(se)%x_i_sbie(:,kn)+ep3*r_integration
              element(se)%cbl_x_i(:,4,kn)=element(se)%x_i_sbie(:,kn)-ep3*r_integration
            end if
            if (node(sn)%sbie.eq.fbem_sbie_mca) then
              element(se)%cbl_x_i(:,1,kn)=element(se)%x_i_sbie_mca(:,kn)+ep2*r_integration
              element(se)%cbl_x_i(:,2,kn)=element(se)%x_i_sbie_mca(:,kn)-ep2*r_integration
              element(se)%cbl_x_i(:,3,kn)=element(se)%x_i_sbie_mca(:,kn)+ep3*r_integration
              element(se)%cbl_x_i(:,4,kn)=element(se)%x_i_sbie_mca(:,kn)-ep3*r_integration
            end if
          end do
          deallocate (x_nodes,xi_nodes)

        end do

      ! ------------------------------ !
      ! 3D: FE SHELL - BE SURFACE LOAD !
      ! ------------------------------ !

      case (fbem_bl_coupling_shell_surface)

        ! Loop through the ELEMENTS of the BE BODY LOAD
        do ke=1,part(be_bodyload(kb)%part)%n_elements
          ! Selected element
          se=part(be_bodyload(kb)%part)%element(ke)
          ! Copy the coordinates of the geometrical nodes
          allocate (x_nodes(problem%n,element(se)%n_nodes))
          x_nodes=element(se)%x_gn
          ! Copy the reference coordinates of the geometrical nodes
          allocate (xi_nodes(element(se)%n_dimension,element(se)%n_nodes))
          xi_nodes=element(se)%xi_gn
          !
          ! CALCULATE XI_I FOR THE FORMULATIONS OF THE NODES
          !
          ! If the element is continuous
          if (element(se)%discontinuous.eqv.(.false.)) then
            ! Switch between element dimensions
            select case (element(se)%n_dimension)
              ! 1D
              case (1)
                ! Allocate "xi_i" members
                allocate (element(se)%xi_i_sbie(1,element(se)%n_nodes))
                allocate (element(se)%xi_i_sbie_mca(1,element(se)%n_nodes))
                allocate (element(se)%xi_i_hbie(1,element(se)%n_nodes))
                ! Loop through the nodes of the element
                do kn=1,element(se)%n_nodes
                  ! Selected node
                  sn=element(se)%node(kn)
                  ! If the node has SBIE, then copy directly the xi coordinate of the geometrical node to the "xi_i_sbie" member
                  if (node(sn)%sbie.eq.fbem_sbie) then
                    element(se)%xi_i_sbie(1,kn)=xi_nodes(1,kn)
                  end if
                  ! If the node has SBIE MCA, then move towards inside the element the xi coordinate of the geometrical node,
                  ! and copy it to the "xi_i_sbie_mca" member
                  if (node(sn)%sbie.eq.fbem_sbie_mca) then
                    ! Save xi coordinate of the geometrical node
                    xi_1d=xi_nodes(1,kn)
                    ! Move the xi coordinate
                    xi_1d=fbem_move_xi_from_vertex(xi_1d,element(se)%delta_sbie_mca(kn))
                    ! Copy it
                    element(se)%xi_i_sbie_mca(1,kn)=xi_1d
                  end if
                  ! If the node has HBIE formulation, then move towards inside the element the xi coordinate of
                  ! the geometrical node, and copy it to the "xi_i_hbie" member
                  if (node(sn)%hbie.eq.fbem_hbie) then
                    ! Save xi coordinate of the geometrical node
                    xi_1d=xi_nodes(1,kn)
                    ! Move the xi coordinate
                    xi_1d=fbem_move_xi_from_vertex(xi_1d,element(se)%delta_hbie(kn))
                    ! Copy it
                    element(se)%xi_i_hbie(1,kn)=xi_1d
                  end if
                end do
              ! 2D
              case (2)
                ! Allocate "xi_i" member
                allocate (element(se)%xi_i_sbie(2,element(se)%n_nodes))
                allocate (element(se)%xi_i_sbie_mca(2,element(se)%n_nodes))
                allocate (element(se)%xi_i_hbie(2,element(se)%n_nodes))
                ! Loop through the nodes of the element
                do kn=1,element(se)%n_nodes
                  ! Selected node
                  sn=element(se)%node(kn)
                  ! If the node has SBIE or SBIE + SBIE MCA formulation, then copy directly the xi coordinate of the geometrical node
                  ! to the "xi_i_sbie" member
                  if (node(sn)%sbie.eq.fbem_sbie) then
                    element(se)%xi_i_sbie(1,kn)=xi_nodes(1,kn)
                    element(se)%xi_i_sbie(2,kn)=xi_nodes(2,kn)
                  end if
                  ! If the node has SBIE MCA, then move towards inside the element the xi coordinate of the geometrical node,
                  ! and copy it to the "xi_i_sbie_mca" member
                  if (node(sn)%sbie.eq.fbem_sbie_mca) then
                    ! Save xi coordinate of the geometrical node
                    xi_2d(1)=xi_nodes(1,kn)
                    xi_2d(2)=xi_nodes(2,kn)
                    ! Move the xi coordinate
                    xi_2d=fbem_move_xi1xi2_from_edge(element(se)%type,xi_2d,element(se)%delta_sbie_mca(kn))
                    ! Copy it
                    element(se)%xi_i_sbie_mca(1,kn)=xi_2d(1)
                    element(se)%xi_i_sbie_mca(2,kn)=xi_2d(2)
                  end if
                  ! If the node has HBIE formulation, then move towards inside the element the xi coordinate of
                  ! the geometrical node, and copy it to the "xi_i_hbie" member
                  if (node(sn)%hbie.eq.fbem_hbie) then
                    ! Save xi coordinate of the geometrical node
                    xi_2d(1)=xi_nodes(1,kn)
                    xi_2d(2)=xi_nodes(2,kn)
                    ! Move the xi coordinate
                    xi_2d=fbem_move_xi1xi2_from_edge(element(se)%type,xi_2d,element(se)%delta_hbie(kn))
                    ! Copy it
                    element(se)%xi_i_hbie(1,kn)=xi_2d(1)
                    element(se)%xi_i_hbie(2,kn)=xi_2d(2)
                  end if
                end do
            end select
          ! If the element is discontinuous
          else
            allocate (element(se)%xi_i_sbie(element(se)%n_dimension,element(se)%n_nodes))
            element(se)%xi_i_sbie=fbem_xi_hybrid(element(se)%type,element(se)%delta_f)
            allocate (element(se)%xi_i_hbie(element(se)%n_dimension,element(se)%n_nodes))
            element(se)%xi_i_hbie=fbem_xi_hybrid(element(se)%type,element(se)%delta_f)
          end if
          !
          ! Calculate x_i
          !
          ! Allocate "x_i" members
          allocate (element(se)%x_i_sbie(problem%n,element(se)%n_nodes))
          allocate (element(se)%x_i_sbie_mca(problem%n,element(se)%n_nodes))
          allocate (element(se)%x_i_hbie(problem%n,element(se)%n_nodes))
          ! Switch element dimensions
          select case (element(se)%n_dimension)
            ! 1D
            case (1)
              ! Switch problem dimensions
              select case (problem%n)
                ! 2D
                case (2)
                  do kn=1,element(se)%n_nodes
                    ! Selected node
                    sn=element(se)%node(kn)
                    ! If the node has SBIE formulation
                    if (node(sn)%sbie.eq.fbem_sbie) then
                      xi_1d=element(se)%xi_i_sbie(1,kn)
                      x_i=fbem_position2d(element(se)%type,x_nodes,xi_1d)
                      do kc=1,2
                        element(se)%x_i_sbie(kc,kn)=x_i(kc)
                      end do
                    end if
                    ! If the node has SBIE MCA formulation
                    if (node(sn)%sbie.eq.fbem_sbie_mca) then
                      xi_1d=element(se)%xi_i_sbie_mca(1,kn)
                      x_i=fbem_position2d(element(se)%type,x_nodes,xi_1d)
                      do kc=1,2
                        element(se)%x_i_sbie_mca(kc,kn)=x_i(kc)
                      end do
                    end if
                    ! If the node has HBIE formulation
                    if (node(sn)%hbie.eq.fbem_hbie) then
                      xi_1d=element(se)%xi_i_hbie(1,kn)
                      x_i=fbem_position2d(element(se)%type,x_nodes,xi_1d)
                      do kc=1,2
                        element(se)%x_i_hbie(kc,kn)=x_i(kc)
                      end do
                    end if
                  end do
                ! 3D
                case (3)
                  do kn=1,element(se)%n_nodes
                    ! Selected node
                    sn=element(se)%node(kn)
                    ! If the node has SBIE formulation
                    if (node(sn)%sbie.eq.fbem_sbie) then
                      xi_1d=element(se)%xi_i_sbie(1,kn)
                      x_i=fbem_position3d(element(se)%type,x_nodes,xi_1d)
                      do kc=1,3
                        element(se)%x_i_sbie(kc,kn)=x_i(kc)
                      end do
                    end if
                    ! If the node has SBIE MCA formulation
                    if (node(sn)%sbie.eq.fbem_sbie_mca) then
                      xi_1d=element(se)%xi_i_sbie_mca(1,kn)
                      x_i=fbem_position3d(element(se)%type,x_nodes,xi_1d)
                      do kc=1,3
                        element(se)%x_i_sbie_mca(kc,kn)=x_i(kc)
                      end do
                    end if
                    ! If the node has HBIE formulation
                    if (node(sn)%hbie.eq.fbem_hbie) then
                      xi_1d=element(se)%xi_i_hbie(1,kn)
                      x_i=fbem_position3d(element(se)%type,x_nodes,xi_1d)
                      do kc=1,3
                        element(se)%x_i_hbie(kc,kn)=x_i(kc)
                      end do
                    end if
                  end do
              end select
            ! 2D
            case (2)
              ! Switch problem dimensions
              select case (problem%n)
                ! 2D
                case (2)
                  do kn=1,element(se)%n_nodes
                    ! Selected node
                    sn=element(se)%node(kn)
                    ! If the node has SBIE formulation
                    if (node(sn)%sbie.eq.fbem_sbie) then
                      xi_2d(1)=element(se)%xi_i_sbie(1,kn)
                      xi_2d(2)=element(se)%xi_i_sbie(2,kn)
                      x_i=fbem_position2d(element(se)%type,x_nodes,xi_2d)
                      do kc=1,2
                        element(se)%x_i_sbie(kc,kn)=x_i(kc)
                      end do
                    end if
                    ! If the node has SBIE MCA formulation
                    if (node(sn)%sbie.eq.fbem_sbie_mca) then
                      xi_2d(1)=element(se)%xi_i_sbie_mca(1,kn)
                      xi_2d(2)=element(se)%xi_i_sbie_mca(2,kn)
                      x_i=fbem_position2d(element(se)%type,x_nodes,xi_2d)
                      do kc=1,2
                        element(se)%x_i_sbie_mca(kc,kn)=x_i(kc)
                      end do
                    end if
                    ! If the node has HBIE formulation
                    if (node(sn)%hbie.eq.fbem_hbie) then
                      xi_2d(1)=element(se)%xi_i_hbie(1,kn)
                      xi_2d(2)=element(se)%xi_i_hbie(2,kn)
                      x_i=fbem_position2d(element(se)%type,x_nodes,xi_2d)
                      do kc=1,2
                        element(se)%x_i_hbie(kc,kn)=x_i(kc)
                      end do
                    end if
                  end do
                ! 3D
                case (3)
                  do kn=1,element(se)%n_nodes
                    ! Selected node
                    sn=element(se)%node(kn)
                    ! If the node has SBIE formulation
                    if (node(sn)%sbie.eq.fbem_sbie) then
                      xi_2d(1)=element(se)%xi_i_sbie(1,kn)
                      xi_2d(2)=element(se)%xi_i_sbie(2,kn)
                      x_i=fbem_position3d(element(se)%type,x_nodes,xi_2d)
                      do kc=1,3
                        element(se)%x_i_sbie(kc,kn)=x_i(kc)
                      end do
                    end if
                    ! If the node has SBIE MCA formulation
                    if (node(sn)%sbie.eq.fbem_sbie_mca) then
                      xi_2d(1)=element(se)%xi_i_sbie_mca(1,kn)
                      xi_2d(2)=element(se)%xi_i_sbie_mca(2,kn)
                      x_i=fbem_position3d(element(se)%type,x_nodes,xi_2d)
                      do kc=1,3
                        element(se)%x_i_sbie_mca(kc,kn)=x_i(kc)
                      end do
                    end if
                    ! If the node has HBIE formulation
                    if (node(sn)%hbie.eq.fbem_hbie) then
                      xi_2d(1)=element(se)%xi_i_hbie(1,kn)
                      xi_2d(2)=element(se)%xi_i_hbie(2,kn)
                      x_i=fbem_position3d(element(se)%type,x_nodes,xi_2d)
                      do kc=1,3
                        element(se)%x_i_hbie(kc,kn)=x_i(kc)
                      end do
                    end if
                  end do
              end select
          end select
          !
          ! Calculate n_i
          !
          ! It is only valid done for 1D elements in 2D and 2D elements in 3D, otherwise a warning message is displayed.
          !
          ! Allocate "n_i" member
          allocate (element(se)%n_i_sbie(problem%n,element(se)%n_nodes))
          allocate (element(se)%n_i_sbie_mca(problem%n,element(se)%n_nodes))
          allocate (element(se)%n_i_hbie(problem%n,element(se)%n_nodes))
          ! 1D element in 2D
          if ((element(se)%n_dimension.eq.1).and.(problem%n.eq.2)) then
            do kn=1,element(se)%n_nodes
              ! Selected node
              sn=element(se)%node(kn)
              ! If the node has SBIE formulation
              if (node(sn)%sbie.eq.fbem_sbie) then
                xi_1d=element(se)%xi_i_sbie(1,kn)
                n_i=fbem_unormal2d(element(se)%type,x_nodes,xi_1d)
                do kc=1,2
                  element(se)%n_i_sbie(kc,kn)=n_i(kc)
                end do
              end if
              ! If the node has SBIE MCA formulation
              if (node(sn)%sbie.eq.fbem_sbie_mca) then
                xi_1d=element(se)%xi_i_sbie_mca(1,kn)
                n_i=fbem_unormal2d(element(se)%type,x_nodes,xi_1d)
                do kc=1,2
                  element(se)%n_i_sbie_mca(kc,kn)=n_i(kc)
                end do
              end if
              ! If the node has HBIE formulation
              if (node(sn)%hbie.eq.fbem_hbie) then
                xi_1d=element(se)%xi_i_hbie(1,kn)
                n_i=fbem_unormal2d(element(se)%type,x_nodes,xi_1d)
                do kc=1,2
                  element(se)%n_i_hbie(kc,kn)=n_i(kc)
                end do
              end if
            end do
          end if
          ! 2D elements in 3D
          if ((element(se)%n_dimension.eq.2).and.(problem%n.eq.3)) then
            do kn=1,element(se)%n_nodes
              ! Selected node
              sn=element(se)%node(kn)
              ! If the node has SBIE formulation
              if (node(sn)%sbie.eq.fbem_sbie) then
                xi_2d(1)=element(se)%xi_i_sbie(1,kn)
                xi_2d(2)=element(se)%xi_i_sbie(2,kn)
                n_i=fbem_unormal3d(element(se)%type,x_nodes,xi_2d)
                do kc=1,3
                  element(se)%n_i_sbie(kc,kn)=n_i(kc)
                end do
              end if
              ! If the node has SBIE MCA formulation
              if (node(sn)%sbie.eq.fbem_sbie_mca) then
                xi_2d(1)=element(se)%xi_i_sbie_mca(1,kn)
                xi_2d(2)=element(se)%xi_i_sbie_mca(2,kn)
                n_i=fbem_unormal3d(element(se)%type,x_nodes,xi_2d)
                do kc=1,3
                  element(se)%n_i_sbie_mca(kc,kn)=n_i(kc)
                end do
              end if
              ! If the node has HBIE formulation
              if (node(sn)%hbie.eq.fbem_hbie) then
                xi_2d(1)=element(se)%xi_i_hbie(1,kn)
                xi_2d(2)=element(se)%xi_i_hbie(2,kn)
                n_i=fbem_unormal3d(element(se)%type,x_nodes,xi_2d)
                do kc=1,3
                  element(se)%n_i_hbie(kc,kn)=n_i(kc)
                end do
              end if
            end do
          end if
          ! Deallocate
          deallocate (x_nodes,xi_nodes)
        end do

      ! -----------------------
      ! FE SHELL - BE EDGE LOAD
      ! -----------------------

      case (fbem_bl_coupling_shell_edge)
        do ke=1,part(sp)%n_elements
          se=part(sp)%element(ke)

          !
          ! Remember that this elements have in reality 1 higher dimension than indicated in %n_dimension, one must look
          ! at type_g for geometrical aspects, and type_f2 for functional aspects
          !
          ! Geometry
          tmp_n_gnodes=fbem_n_nodes(element(se)%type_g)
          tmp_n_dimension=fbem_n_dimension(element(se)%type_g)
          allocate (x_nodes(problem%n,tmp_n_gnodes))
          x_nodes=element(se)%x_gn
          allocate (xi_nodes(tmp_n_dimension,tmp_n_gnodes))
          xi_nodes=element(se)%xi_gn
          ! Number of nodes of functional interpolation
          tmp_f2nodes=fbem_n_nodes(element(se)%type_f2)

          !
          ! Calculate xi_i for the formulations of the nodes
          !
          ! Allocate
          allocate (element(se)%cbl_n_cp(tmp_f2nodes))
          allocate (element(se)%cbl_xi_i(2,1,tmp_f2nodes))
          ! Loop through the nodes of the root element (longitudinal direction xi1 / k1 of the derived element)
          do kn1=1,element(se)%n_nodes
            ! Selected root node
            sn=element(se)%node(kn1)
            ! Loop through the nodes in the thickness direction (xi2 / k2)
            do kn2=1,tmp_f2nodes/element(se)%n_nodes
              ! Node index of the derived element
              kn=(kn2-1)*element(se)%n_nodes+kn1
              ! Only one collocation point
              element(se)%cbl_n_cp(kn)=1
              ! Root node with SBIE
              if (node(sn)%sbie.eq.fbem_sbie) then
                element(se)%cbl_xi_i(:,1,kn)=element(se)%xi_fn(:,kn)
              end if
              ! Root node with SBIE MCA (only move the xi1 coordinate)
              if (node(sn)%sbie.eq.fbem_sbie_mca) then
                ! Coordinate of the functional node
                xi_2d=element(se)%xi_fn(:,kn)
                ! Move the xi1 towards inside coordinate
                xi_1d=xi_2d(1)
                xi_1d=fbem_move_xi_from_vertex(xi_1d,element(se)%delta_sbie_mca(kn1))
                xi_2d(1)=xi_1d
                ! Copy it
                element(se)%cbl_xi_i(:,1,kn)=xi_2d
              end if
            end do
          end do
          !
          ! Calculate x^i
          !
          allocate (element(se)%cbl_x_i(problem%n,maxval(element(se)%cbl_n_cp),tmp_f2nodes))
          do kn=1,tmp_f2nodes
            do kc=1,element(se)%cbl_n_cp(kn)
              xi_2d=element(se)%cbl_xi_i(:,kc,kn)
              x_i=fbem_position3d(element(se)%type_g,x_nodes,xi_2d)
              element(se)%cbl_x_i(:,kc,kn)=x_i
            end do
          end do
          ! Deallocate
          deallocate (x_nodes,xi_nodes)

        end do

    end select
  end do ! Loop through COUPLED BE BODY LOAD ELEMENTS

  ! Ending message
  if (verbose_level.ge.1) write(output_unit,'(a)') 'done.'

end subroutine build_data_at_collocation_points
