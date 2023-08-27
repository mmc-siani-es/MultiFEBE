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

subroutine assign_solution_mechanics_static

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_geometry
  use fbem_symmetry

  ! Module of problem variables
  use problem_variables

  ! No implicit variables are allowed
  implicit none

  ! Local variables
  integer                        :: kr, kb, ke, kn, ks
  integer                        :: sb, se, sn, ss, sp, snm
  integer                        :: k, kt, km
  logical, allocatable           :: node_used(:)
  integer                        :: sn_fe
  real(kind=real64), allocatable :: T(:,:)

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'START assigning the LSE solution to variables')

  allocate (node_used(n_nodes))
  node_used=.false.

  ! =============
  ! RIGID REGIONS
  ! =============

  ! Loop through regions
  do kr=1,n_regions
    if ((region(kr)%class.eq.fbem_fe).and.(region(kr)%type.eq.fbem_rigid)) then
      !
      ! Solution of the master node
      !
      sn=region(kr)%master_node
      node_used(sn)=.true.
      do k=1,node(sn)%n_dof
        if (node(sn)%ctype(k,1).eq.1) then
          node(sn)%value_r(k,1)=b_r(node(sn)%col(k,1),1)
        else
          node(sn)%value_r(k,1)=node(sn)%cvalue_r(k,1,1)
        end if
      end do
      !
      ! Solution of the slave nodes
      !
      snm=sn
      do kn=1,node(snm)%n_slaves
        sn=node(snm)%slave(kn)
        node_used(sn)=.true.
        allocate (T(node(sn)%n_dof,node(snm)%n_dof))
        call fbem_rigid_solid_transformation_matrix(problem%n,node(sn)%x,node(snm)%x,node(sn)%n_dof,T)
        do k=1,node(sn)%n_dof
          node(sn)%value_r(k,1)=0
          do km=1,node(snm)%n_dof
            node(sn)%value_r(k,1)=node(sn)%value_r(k,1)+T(k,km)*node(snm)%value_r(km,1)
          end do
        end do
        deallocate (T)
      end do
    end if
  end do

  ! ================
  ! FLEXIBLE REGIONS
  ! ================

  ! Loop through regions
  do kr=1,n_regions
    select case (region(kr)%class)

      ! =========
      ! BE REGION
      ! =========

      case (fbem_be)

        ! =============
        ! BE BOUNDARIES
        ! =============

        do kb=1,region(kr)%n_boundaries
          sb=region(kr)%boundary(kb)
          do ke=1,part(boundary(sb)%part)%n_elements
            se=part(boundary(sb)%part)%element(ke)
            do kn=1,element(se)%n_nodes
              sn=element(se)%node(kn)
              if (.not.node_used(sn)) then
                node_used(sn)=.true.
                select case (boundary(sb)%coupling)

                  ! ===========
                  ! BE BOUNDARY
                  ! ===========

                  case (fbem_boundary_coupling_be)
                    select case (boundary(sb)%class)

                      ! =================
                      ! ORDINARY BOUNDARY
                      ! =================

                      case (fbem_boundary_class_ordinary)
                        ! Index of variables for each coordinate k:
                        ! node(sn)%col(          k,1): u_k
                        ! node(sn)%col(problem%n+k,1): p_k
                        do k=1,problem%n
                          select case (node(sn)%ctype(k,1))
                            case (0)
                              node(sn)%value_r(          k,1)=node(sn)%cvalue_r(k,1,1)
                              node(sn)%value_r(problem%n+k,1)=b_r(node(sn)%col(problem%n+k,1),1)
                            case (1)
                              node(sn)%value_r(          k,1)=b_r(node(sn)%col(k,1),1)
                              node(sn)%value_r(problem%n+k,1)=node(sn)%cvalue_r(k,1,1)
                            case (10)
                              node(sn)%value_r(          k,1)=b_r(node(sn)%col(k,1),1)
                              if (.not.region(kr)%boundary_reversion(kb)) then
                                node(sn)%value_r(problem%n+k,1)=node(sn)%cvalue_r(k,1,1)*node(sn)%n_fn(k)
                              else
                                node(sn)%value_r(problem%n+k,1)=-node(sn)%cvalue_r(k,1,1)*node(sn)%n_fn(k)
                              end if
                          end select
                        end do

                      ! ===================
                      ! CRACK-LIKE BOUNDARY
                      ! ===================

                      case (fbem_boundary_class_cracklike)
                        ! Index of variables for each coordinate k:
                        ! node(sn)%col(          k,1): u_k for face +
                        ! node(sn)%col(problem%n+k,1): p_k for face +
                        ! node(sn)%col(          k,2): u_k for face -
                        ! node(sn)%col(problem%n+k,2): p_k for face -
                        ! Face +
                        do k=1,problem%n
                          select case (node(sn)%ctype(k,1))
                            case (0)
                              node(sn)%value_r(          k,1)=node(sn)%cvalue_r(k,1,1)
                              node(sn)%value_r(problem%n+k,1)=b_r(node(sn)%col(problem%n+k,1),1)
                            case (1)
                              node(sn)%value_r(          k,1)=b_r(node(sn)%col(k,1),1)
                              node(sn)%value_r(problem%n+k,1)=node(sn)%cvalue_r(k,1,1)
                          end select
                        end do
                        ! Face -
                        do k=1,problem%n
                          select case (node(sn)%ctype(k,2))
                            case (0)
                              node(sn)%value_r(          k,2)=node(sn)%cvalue_r(k,1,2)
                              node(sn)%value_r(problem%n+k,2)=b_r(node(sn)%col(problem%n+k,2),1)
                            case (1)
                              node(sn)%value_r(          k,2)=b_r(node(sn)%col(k,2),1)
                              node(sn)%value_r(problem%n+k,2)=node(sn)%cvalue_r(k,1,2)
                          end select
                        end do
                    end select

                  ! ==============
                  ! BE-BE BOUNDARY
                  ! ==============

                  case (fbem_boundary_coupling_be_be)
                    ! Index of variables for each coordinate k:
                    ! node(sn)%col(          k,1): u_k for region 1
                    ! node(sn)%col(problem%n+k,1): p_k for region 1
                    ! node(sn)%col(          k,2): u_k for region 2 (not used since u_k_{region 1}=u_k_{region 2})
                    ! node(sn)%col(problem%n+k,2): p_k for region 2 (not used since p_k_{region 1}=-p_k_{region 2})
                    do k=1,problem%n
                      node(sn)%value_r(          k,1)=b_r(node(sn)%col(k,1),1)
                      node(sn)%value_r(problem%n+k,1)=b_r(node(sn)%col(problem%n+k,1),1)
                      node(sn)%value_r(          k,2)= node(sn)%value_r(          k,1)
                      node(sn)%value_r(problem%n+k,2)=-node(sn)%value_r(problem%n+k,1)
                    end do

                  ! ==============
                  ! BE-FE BOUNDARY
                  ! ==============

                  case (fbem_boundary_coupling_be_fe)
                    ! The fe node connected with the be node
                    sn_fe=element(se)%element_node(kn)
                    select case (boundary(sb)%class)

                      ! =================
                      ! ORDINARY BOUNDARY
                      ! =================

                      case (fbem_boundary_class_ordinary)
                        ! Index of variables for each coordinate k:
                        ! node(sn)%col(          k,1): u_k
                        ! node(sn)%col(problem%n+k,1): t_k
                        do k=1,problem%n
                          if (node(sn_fe)%ctype(k,1).eq.1) then
                            node(sn)%value_r(k,1)=b_r(node(sn_fe)%col(k,1),1)
                          else
                            node(sn)%value_r(k,1)=node(sn_fe)%cvalue_r(k,1,1)
                          end if
                          node(sn)%value_r(problem%n+k,1)=b_r(node(sn)%col(problem%n+k,1),1)
                        end do

                      ! ===================
                      ! CRACK-LIKE BOUNDARY
                      ! ===================

                      case (fbem_boundary_class_cracklike)
!                        !
!                        ! SBIE (solve for Σt = t^+ + t^-)
!                        !
!                        ! Index of variable values for each coordinate k:
!                        ! node(sn)%value_r(          k,1): u_k^+
!                        ! node(sn)%value_r(problem%n+k,1): t_k^+ (temporarily before calculating t_k^+ and t_k^- it will be Σt)
!                        ! node(sn)%value_r(          k,2): u_k^-
!                        ! node(sn)%value_r(problem%n+k,2): t_k^- (temporarily before calculating t_k^+ and t_k^- it will be 0)
!                        do k=1,problem%n
!                          if (node(sn_fe)%ctype(k,1).eq.1) then
!                            node(sn)%value_r(k,1)=b_r(node(sn_fe)%col(k,1),1)
!                          else
!                            node(sn)%value_r(k,1)=node(sn_fe)%cvalue_r(k,1,1)
!                          end if
!                          node(sn)%value_r(problem%n+k,1)=b_r(node(sn)%col(problem%n+k,1),1)
!                          node(sn)%value_r(          k,2)=node(sn)%value_r(k,1)
!                          node(sn)%value_r(problem%n+k,2)=0.d0
!                        end do

                        !
                        ! Usando Dual BEM (para CC mas generales)
                        !
                        ! Index of variables for each coordinate k:
                        ! node(sn)%col(          k,1): u_k for face +
                        ! node(sn)%col(problem%n+k,1): t_k for face +
                        ! node(sn)%col(          k,2): u_k for face - (not used since u_l_{+} = u_k_{-})
                        ! node(sn)%col(problem%n+k,2): t_k for face -
                        do k=1,problem%n
                          if (node(sn_fe)%ctype(k,1).eq.1) then
                            node(sn)%value_r(k,1)=b_r(node(sn_fe)%col(k,1),1)
                          else
                            node(sn)%value_r(k,1)=node(sn_fe)%cvalue_r(k,1,1)
                          end if
                          node(sn)%value_r(problem%n+k,1)=b_r(node(sn)%col(problem%n+k,1),1)
                          node(sn)%value_r(          k,2)=node(sn)%value_r(k,1)
                          node(sn)%value_r(problem%n+k,2)=b_r(node(sn)%col(problem%n+k,2),1)
                        end do

                    end select

                  ! =================
                  ! BE-FE-BE BOUNDARY
                  ! =================

                  case (fbem_boundary_coupling_be_fe_be)
                    ! The fe node connected with the be node
                    sn_fe=element(se)%element_node(kn)
                    ! Index of variables for each coordinate k:
                    ! node(sn)%col(          k,1): u_k for region 1
                    ! node(sn)%col(problem%n+k,1): t_k for region 1
                    ! node(sn)%col(          k,2): u_k for region 2 (not used since u_k_{region 1}=u_k_{region 2})
                    ! node(sn)%col(problem%n+k,2): t_k for region 2
                    do k=1,problem%n
                      if (node(sn_fe)%ctype(k,1).eq.1) then
                        node(sn)%value_r(k,1)=b_r(node(sn_fe)%col(k,1),1)
                      else
                        node(sn)%value_r(k,1)=node(sn_fe)%cvalue_r(k,1,1)
                      end if
                      node(sn)%value_r(problem%n+k,1)=b_r(node(sn)%col(problem%n+k,1),1)
                      node(sn)%value_r(          k,2)=node(sn)%value_r(k,1)
                      node(sn)%value_r(problem%n+k,2)=b_r(node(sn)%col(problem%n+k,2),1)
                    end do

                end select
              end if
            end do
          end do
        end do

        ! =====================
        ! COUPLED BE BODY LOADS
        ! =====================

        do kb=1,region(kr)%n_be_bodyloads
          sb=region(kr)%be_bodyload(kb)
          sp=be_bodyload(sb)%part
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
                    do kn=1,element(se)%n_nodes
                      sn=element(se)%node(kn)
                      if (.not.node_used(sn)) then
                        node_used(sn)=.true.

                        ! Index of variables:
                        ! - node(sn)%col(          k,nf): u_k for each coordinate k and each functional node (thickness) associated with node sn (not active since equal to tip/edge displacement field)
                        ! - node(sn)%col(problem%n+k,nf): b_k for each coordinate k and each functional node (thickness) associated with node sn
                        ! Assign row to the equation and column to the variables
                        do kt=1,fbem_n_nodes(element(se)%type_f2)/element(se)%n_nodes
                          do k=1,problem%n
                            node(sn)%value_r(problem%n+k,kt)=b_r(node(sn)%col(problem%n+k,kt),1)

                            ! node(sn)%col(          k,nf): u_k tiene que sacarse por puntos internos si se desea

                          end do
                        end do

                      end if
                    end do
                  end do
                !
                ! 3D: SURFACE LOAD ATH THE BEAM TIP
                !
                case (3)
                  ! Modelo de Luis (en principio solo habria que crear un cuadrilatero que emule a una circunferencia, se puede hacer..)
                  stop 'not yet 37'
              end select

            ! -----------------------
            ! FE SHELL - BE EDGE LOAD
            ! -----------------------

            case (fbem_bl_coupling_shell_edge)
              do ke=1,part(sp)%n_elements
                se=part(sp)%element(ke)
                do kn=1,element(se)%n_nodes
                  sn=element(se)%node(kn)
                  if (.not.node_used(sn)) then
                    node_used(sn)=.true.




                    ! Index of variables:
                    ! - node(sn)%col(          k,nf): u_k for each coordinate k and each functional node (thickness) associated with node sn (not active since equal to tip/edge displacement field)
                    ! - node(sn)%col(problem%n+k,nf): b_k for each coordinate k and each functional node (thickness) associated with node sn
                    ! Assign row to the equation and column to the variables
                    do kt=1,fbem_n_nodes(element(se)%type_f2)/element(se)%n_nodes
                      do k=1,problem%n
                        node(sn)%value_r(problem%n+k,kt)=b_r(node(sn)%col(problem%n+k,kt),1)

                        ! node(sn)%col(          k,nf): u_k tiene que sacarse por puntos internos si se desea

                      end do
                    end do




                  end if
                end do
              end do

            ! -----------------------------------------------------
            ! FE BEAM - BE LINE LOAD AND FE SHELL - BE SURFACE LOAD
            ! -----------------------------------------------------

            case (fbem_bl_coupling_beam_line,fbem_bl_coupling_shell_surface)
              do ke=1,part(sp)%n_elements
                se=part(sp)%element(ke)
                do kn=1,element(se)%n_nodes
                  sn=element(se)%node(kn)
                  if (.not.node_used(sn)) then
                    node_used(sn)=.true.
                    sn_fe=element(se)%element_node(kn)

                    if (node(sn_fe)%n_dof.ge.problem%n) then

                      do k=1,problem%n
                        if (node(sn_fe)%ctype(k,1).eq.1) then
                          node(sn)%value_r(k,1)=b_r(node(sn_fe)%col(k,1),1)
                        else
                          node(sn)%value_r(k,1)=node(sn_fe)%cvalue_r(k,1,1)
                        end if
                      end do

                    else

                      call fbem_warning_message(error_unit,0,'node',node(sn_fe)%id,'has <n dof, not assigning to be bodyload')

                    end if

                    do k=1,problem%n
                      node(sn)%value_r(problem%n+k,1)=b_r(node(sn)%col(problem%n+k,1),1)
                    end do

                  end if
                end do
              end do

          end select

        end do

      ! =========
      ! FE REGION
      ! =========

      case (fbem_fe)
        !
        ! FLEXIBLE REGIONS
        !
        if (region(kr)%type.ne.fbem_rigid) then
          do ks=1,region(kr)%n_fe_subregions
            ss=region(kr)%fe_subregion(ks)
            do ke=1,part(fe_subregion(ss)%part)%n_elements
              se=part(fe_subregion(ss)%part)%element(ke)
              do kn=1,element(se)%n_nodes
                sn=element(se)%node(kn)
                if (.not.node_used(sn)) then
                  node_used(sn)=.true.
                  do k=1,node(sn)%n_dof
                    if (node(sn)%ctype(k,1).eq.1) then
                      node(sn)%value_r(k,1)=b_r(node(sn)%col(k,1),1)
                    else
                      node(sn)%value_r(k,1)=node(sn)%cvalue_r(k,1,1)
                    end if
                  end do
                end if
              end do
            end do
          end do
        end if

    end select
  end do

  if (verbose_level.ge.1) call fbem_timestamp_w_message(output_unit,2,'END assigning the LSE solution to variables')

end subroutine assign_solution_mechanics_static
