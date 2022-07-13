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

subroutine assemble_fem_staela_coupling(se_int,mode)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! Problem variables module
  use problem_variables
  use fbem_numerical
  use fbem_shape_functions
  use fbem_fem_solids
  use fbem_fem_beams
  use fbem_fem_shells

  ! No implicit variables are allowed
  implicit none

  ! I/O variables
  integer                        :: se_int                                     !! FE
  integer                        :: mode                                       !! mode=0 (to Ax=b), mode=1 (to q=Ka-f_bem-f_ext)

  ! Local variables
  integer                        :: kn_row, sn_row, kn, sn, kc, sr
  integer                        :: kdof_row, kdof_col, kdofe
  integer                        :: row
  integer                        :: row_local
  integer                        :: col
  integer                        :: col_local
  integer                        :: se_be, sb_be, sn_be
  integer                        :: kd_fe, sd_fe, sd_type, sd_fe_n_nodes
  integer                        :: se_be_type, se_be_n_nodes
  integer                        :: ndof_se, ndof_load
  real(kind=real64)              :: nu
  real(kind=real64), allocatable :: nodal_axes(:,:,:)
  real(kind=real64), allocatable :: Qedge(:,:)
  real(kind=real64), allocatable :: Qmid(:,:)
  integer                        :: se_int_type, se_int_n_nodes, se_int_n_dof
  integer, allocatable           :: se_int_n_dof_node(:)
  real(kind=real64), allocatable :: se_int_x(:,:), se_int_t(:)
  real(kind=real64), allocatable :: t_e(:), f_e(:)

  !
  ! Falta meter el tratamiento para regiones rigidas
  !

  !
  ! Region of the element (in order to take material properties)
  !
  sr=fe_subregion(part(element(se_int)%part)%entity)%region

  ! ===================================
  ! nD CONTINUUM (SOLID) FINITE ELEMENT
  ! ===================================

  if (element(se_int)%n_dimension.eq.problem%n) then
    select case (problem%n)

      ! -----------
      ! 2D SOLID FE
      ! -----------

      case (2)


        ! ------------------
        ! BOUNDARY/EDGE LOAD
        ! ------------------

        ! Nota: hay que alocatar ctype, para los subedges....
        !
        ! Cargas externas¿?¿?¿? - revisar ¿?¿??????????????????????????????????????????????????????????????????

        do kd_fe=1,element(se_int)%n_subedges
          sd_fe=element(se_int)%subedge(kd_fe)
          !
          ! Solo en coord. globales por ahora .......
          !
          if (subedge(sd_fe)%ctype(2,1).gt.0) then
            ! LOAD DATA
            sd_type=subedge(sd_fe)%type
            sd_fe_n_nodes=subedge(sd_fe)%n_nodes
            ! FE DATA
            se_int_type=element(se_int)%type
            se_int_n_nodes=element(se_int)%n_nodes
            allocate (se_int_x(2,se_int_n_nodes))
            allocate (se_int_t(se_int_n_nodes))
            allocate (se_int_n_dof_node(se_int_n_nodes))
            se_int_x=element(se_int)%x_gn
            select case (problem%subtype)
              case (fbem_mechanics_plane_strain)
                se_int_t=1.d0
              case (fbem_mechanics_plane_stress)
                se_int_t=element(se_int)%tn_midnode(3,:)
            end select
            do kn=1,se_int_n_nodes
              sn=element(se_int)%node(kn)
              se_int_n_dof_node(kn)=node(sn)%n_dof
            end do
            se_int_n_dof=sum(se_int_n_dof_node)
            ! CALCULATE Qedge
            allocate (Qedge(2*se_int_n_nodes,2*sd_fe_n_nodes))
            call fbem_fem_ela2d_solid_Q_edge(se_int_type,se_int_x,se_int_t,kd_fe,sd_type,4,Qedge)
            !
            ! ASSEMBLE to Ax=b
            !
            ! Initialize local row counter
            row_local=1
            ! Loop through the FE DOFs
            do kn_row=1,se_int_n_nodes
              sn_row=element(se_int)%node(kn_row)
              do kdof_row=1,se_int_n_dof_node(kn_row)
                if (node(sn_row)%ctype(kdof_row,1).eq.1) then
                  row=node(sn_row)%row(kdof_row,1)
                  ! Initialize local column counter
                  col_local=1
                  ! Loop through FE edge nodes
                  do kn=1,sd_fe_n_nodes
                    ! Note that Qedge is built using the element parent orientation for the edge. Because of that, the reversion
                    ! is needed in general, although boundary load are only allowed in theory on boundaries, and there the subedge
                    ! is not reversed.
                    if (element(se_int)%subedge_reversion(kd_fe)) then
                      do kdof_col=1,problem%n
                        b_r(row,1)=b_r(row,1)+Qedge(row_local,col_local)*subedge(sd_fe)%cvalue_r(kdof_col,fbem_element_invert(kn,sd_type),2)
                        col_local=col_local+1
                      end do
                    else
                      do kdof_col=1,problem%n
                        b_r(row,1)=b_r(row,1)+Qedge(row_local,col_local)*subedge(sd_fe)%cvalue_r(kdof_col,kn,2)
                        col_local=col_local+1
                      end do
                    end if
                  end do
                end if
                row_local=row_local+1
              end do
            end do
            deallocate (se_int_x,se_int_t,se_int_n_dof_node,Qedge)
          end if
        end do

        ! ----------------
        ! BEM-FEM COUPLING
        ! ----------------

        do kd_fe=1,element(se_int)%n_subedges
          sd_fe=element(se_int)%subedge(kd_fe)
          if (subedge(sd_fe)%element.gt.0) then
            ! CONNECTIVITY
            se_be=subedge(sd_fe)%element
            ! FE DATA
            sd_type=subedge(sd_fe)%type
            se_int_type=element(se_int)%type
            se_int_n_nodes=element(se_int)%n_nodes
            allocate (se_int_x(2,se_int_n_nodes))
            allocate (se_int_t(se_int_n_nodes))
            allocate (se_int_n_dof_node(se_int_n_nodes))
            se_int_x=element(se_int)%x_gn
            select case (problem%subtype)
              case (fbem_mechanics_plane_strain)
                se_int_t=1.d0
              case (fbem_mechanics_plane_stress)
                se_int_t=element(se_int)%tn_midnode(3,:)
            end select
            do kn=1,se_int_n_nodes
              sn=element(se_int)%node(kn)
              se_int_n_dof_node(kn)=node(sn)%n_dof
            end do
            se_int_n_dof=sum(se_int_n_dof_node)
            sd_fe_n_nodes=subedge(sd_fe)%n_nodes
            ! BE DATA
            se_be_n_nodes=element(se_be)%n_nodes
            se_be_type=element(se_be)%type
            se_be_n_nodes=element(se_be)%n_nodes
            ! CALCULATE Qedge
            allocate (Qedge(2*se_int_n_nodes,2*se_be_n_nodes))
            call fbem_fem_ela2d_solid_Q_edge(se_int_type,se_int_x,se_int_t,kd_fe,se_be_type,4,Qedge)
            !
            ! ASSEMBLE
            !
            select case (mode)
              !
              ! ASSEMBLE to Ax=b
              !
              case (0)

                ! Initialize local row counter
                row_local=1
                ! Loop through the FE DOF
                do kn_row=1,se_int_n_nodes
                  sn_row=element(se_int)%node(kn_row)
                  do kdof_row=1,se_int_n_dof_node(kn_row)
                    if (node(sn_row)%ctype(kdof_row,1).eq.1) then
                      row=node(sn_row)%row(kdof_row,1)
                      ! Initialize local column counter
                      col_local=1
                      ! Loop through FE edge DOF
                      do kn=1,sd_fe_n_nodes
                        ! Note that Qedge is built using the element parent orientation for the edge. Because of that the reversion
                        ! is needed in general. Although here if the edge is a boundary edge (it should be for coupling), the subedge
                        ! must have the same orientation
                        if (element(se_int)%subedge_reversion(kd_fe)) then
                          sn_be=subedge(sd_fe)%element_node(fbem_element_invert(kn,sd_type))
                        else
                          sn_be=subedge(sd_fe)%element_node(kn)
                        end if
                        do kdof_col=1,problem%n
                          col=node(sn_be)%col(problem%n+kdof_col,1)
                          A_r(row,col)=A_r(row,col)+Qedge(row_local,col_local)
                          col_local=col_local+1
                        end do
                      end do
                    end if
                    row_local=row_local+1
                  end do
                end do

              !
              ! ASSEMBLE to q=Ka-f_bem-f_ext
              !
              case (1)

                ! Initialization
                allocate (t_e(problem%n*sd_fe_n_nodes))
                allocate (f_e(problem%n*se_int_n_nodes))
                f_e=0
                t_e=0
                ! Build the BE load vector (t_e)
                do kn=1,sd_fe_n_nodes
                  if (element(se_int)%subedge_reversion(kd_fe)) then
                    sn_be=subedge(sd_fe)%element_node(fbem_element_invert(kn,sd_type))
                  else
                    sn_be=subedge(sd_fe)%element_node(kn)
                  end if
                  t_e((problem%n*(kn-1)+1):(problem%n*kn))=node(sn_be)%value_r((problem%n+1):2*problem%n,1)
                end do
                ! Build the FE load vector (f_e)
                f_e=matmul(Qedge,t_e)
                ! Save to q
                kdofe=1
                do kn=1,element(se_int)%n_nodes
                  do kc=1,element(se_int)%node_n_dof(kn)
                    element(se_int)%value_r(kc,kn,2)=element(se_int)%value_r(kc,kn,2)+f_e(kdofe)
                    kdofe=kdofe+1
                  end do
                end do
                ! Finalization
                deallocate (t_e,f_e)

              case default

                stop 'invalid mode for coupling assemble'

            end select
            deallocate (se_int_x,se_int_t,se_int_n_dof_node,Qedge)
          end if
        end do

      ! -----------
      ! 3D SOLID FE
      ! -----------

      case (3)
        stop 'not yet 23' ! Debería ser lo mismo que 2D pero con subfaces....

    end select
  end if

  ! ===================================================================================
  ! (n-1)D STRUCTURAL FINITE ELEMENT (SHELL FE IN 3D AMBIENT, OR BEAM FE IN 2D AMBIENT)
  ! ===================================================================================

  if (element(se_int)%n_dimension.eq.(problem%n-1)) then

    ! ----------------------------------------------------------
    ! FE MID-FACE TO BE FACE COUPLING (FE ELEMENT TO BE ELEMENT)
    ! ----------------------------------------------------------

    if (element(se_int)%element.gt.0) then
      !
      ! Calculate Qmid
      !
      ! Initialize
      ndof_se=3*(problem%n-1)*element(se_int)%n_nodes
      ndof_load=problem%n*element(se_int)%n_nodes
      allocate (Qmid(ndof_se,ndof_load))
      ! Calculate
      select case (element(se_int)%n_dimension)

        ! =====================
        ! 1D STRUCTURAL ELEMENT
        ! =====================

        case (1)
          select case (element(se_int)%fe_type)

            ! ------------------------
            ! DEGENERATED BEAM ELEMENT
            ! ------------------------

            case (0)
              call fbem_fem_degbeam_Q_midline(problem%n,element(se_int)%type,element(se_int)%x_gn,&
                                              element(se_int)%Q_intmode,element(se_int)%Q_intngp,Qmid)

            ! ---------------------------------------------------
            ! STRAIGHT EULER-BERNOULLI OR TIMOSHENKO BEAM ELEMENT
            ! ---------------------------------------------------

            case (1,2)
              !
              ! TO-DO: Los ejes nodales habría que guardarlos para cada nodo y construirlos para cada elemento
              !
              allocate(nodal_axes(problem%n,problem%n,element(se_int)%n_nodes))
              nodal_axes=0
              do kn=1,element(se_int)%n_nodes
                do kc=1,problem%n
                  nodal_axes(kc,kc,kn)=1
                end do
              end do
              !
              ! TODO: Usar el paradigma de region()%material(1)
              ! Use material of the FE REGION
              !
              nu=region(sr)%property_r(3)
              !
              ! Calculate
              !
              call fbem_fem_strbeam_Q_midline(problem%n,element(se_int)%type,element(se_int)%fe_options(1),element(se_int)%fe_type,&
                                              element(se_int)%x_gn,element(se_int)%ep,&
                                              element(se_int)%A,element(se_int)%I(1:3),element(se_int)%ksh,&
                                              1.0d0,nu,nodal_axes,Qmid)
            case default
              call fbem_error_message(error_unit,0,'element',element(se_int)%id,'this type of element cannot be coupled to BE')

          end select

        ! =====================
        ! 2D STRUCTURAL ELEMENT
        ! =====================

        case (2)

          ! --------------------------------
          ! DEGENERATED SHELL FINITE ELEMENT
          ! --------------------------------

          call fbem_fem_degshell_Q_midsurface(element(se_int)%type,&
                                              element(se_int)%x_gn,element(se_int)%v_midnode,element(se_int)%tv_midnode,&
                                              element(se_int)%node_n_dof,element(se_int)%Q_intmode,element(se_int)%Q_intngp,Qmid)

      end select
      !
      ! Assemble Qmid
      !
      ! FE -> BE CONNECTIVITY
      se_be=element(se_int)%element
      sb_be=part(element(se_be)%part)%entity
      select case (part(element(se_be)%part)%type)

        ! ===========
        ! BE BOUNDARY
        ! ===========

        case (fbem_part_be_boundary)
          !
          ! ASSEMBLE
          !
          select case (mode)
            !
            ! ASSEMBLE to Ax=b
            !
            case (0)
              ! ASSEMBLE LOAD MATRIX
              ! Local row
              row_local=1
              do kn_row=1,element(se_int)%n_nodes
                sn_row=element(se_int)%node(kn_row)
                do kdof_row=1,element(se_int)%node_n_dof(kn_row)
                  if (node(sn_row)%ctype(kdof_row,1).eq.1) then
                    row=node(sn_row)%row(kdof_row,1)
                    ! Local column
                    col_local=1
                    do kn=1,element(se_int)%n_nodes
                      ! Selected FE node and its corresponding BE node
                      sn=element(se_int)%node(kn)
                      sn_be=element(se_int)%element_node(kn)
                      ! Loop through coordinates
                      do kdof_col=1,problem%n
                        select case (boundary(sb_be)%coupling)

                          ! ==============
                          ! BE-FE BOUNDARY
                          ! ==============

                          case (fbem_boundary_coupling_be_fe)
                            select case (boundary(sb_be)%class)

                              ! -----------------
                              ! ORDINARY BOUNDARY
                              ! -----------------

                              case (fbem_boundary_class_ordinary)
                                ! node(sn_be)%col(problem%n+kdof_col,1): t_k
                                ! q_k = - t_k
                                col=node(sn_be)%col(problem%n+kdof_col,1)
                                A_r(row,col)=A_r(row,col)+Qmid(row_local,col_local)

                              ! -------------------
                              ! CRACK-LIKE BOUNDARY
                              ! -------------------

                              case (fbem_boundary_class_cracklike)
    !                            !
    !                            ! SBIE (Σt_k unknown)
    !                            !
    !                            ! node(sn_be)%col(problem%n+kdof_col,1): Σt_k = t_k^+ + t_k^-
    !                            ! q_k = - t_k^+ - t_k^-
    !                            col=node(sn_be)%col(problem%n+kdof_col,1)
    !                            A_r(row,col)=A_r(row,col)+Qmid(row_local,col_local)
                                !
                                ! DBEM coupling (general interface conditions and individual tractions)
                                !
                                ! node(sn_be)%col(problem%n+kdof_col,1): t_k^+
                                ! node(sn_be)%col(problem%n+kdof_col,2): t_k^-
                                ! q_k = - t_k^+ - t_k^-
                                col=node(sn_be)%col(problem%n+kdof_col,1)
                                A_r(row,col)=A_r(row,col)+Qmid(row_local,col_local)
                                col=node(sn_be)%col(problem%n+kdof_col,2)
                                A_r(row,col)=A_r(row,col)+Qmid(row_local,col_local)

                            end select

                          ! =================
                          ! BE-FE-BE BOUNDARY
                          ! =================

                          case (fbem_boundary_coupling_be_fe_be)
                            ! node(sn_be)%col(problem%n+kdof_col,1): t_k^(1) (region 1)
                            ! node(sn_be)%col(problem%n+kdof_col,2): t_k^(2) (region 2)
                            ! q_k = - t_k^(1) - t_k^(2)
                            col=node(sn_be)%col(problem%n+kdof_col,1)
                            A_r(row,col)=A_r(row,col)+Qmid(row_local,col_local)
                            col=node(sn_be)%col(problem%n+kdof_col,2)
                            A_r(row,col)=A_r(row,col)+Qmid(row_local,col_local)

                        end select
                        col_local=col_local+1
                      end do
                    end do
                  end if
                  row_local=row_local+1
                end do
              end do
            !
            ! ASSEMBLE to q=Ka-f_bem-f_ext
            !
            case (1)


              ! Initialization
              allocate (t_e(ndof_load))
              allocate (f_e(ndof_se))
              f_e=0
              t_e=0
              ! Build the BE load vector (t_e)
              do kn=1,element(se_int)%n_nodes
                ! Selected FE node and its corresponding BE node
                sn=element(se_int)%node(kn)
                sn_be=element(se_int)%element_node(kn)
                select case (boundary(sb_be)%coupling)

                  ! ==============
                  ! BE-FE BOUNDARY
                  ! ==============

                  case (fbem_boundary_coupling_be_fe)
                    select case (boundary(sb_be)%class)

                      ! -----------------
                      ! ORDINARY BOUNDARY
                      ! -----------------

                      case (fbem_boundary_class_ordinary)
                        ! node(sn_be)%col(problem%n+kdof_col,1): t_k
                        ! q_k = - t_k
                        t_e((problem%n*(kn-1)+1):(problem%n*kn))=node(sn_be)%value_r((problem%n+1):2*problem%n,1)

                      ! -------------------
                      ! CRACK-LIKE BOUNDARY
                      ! -------------------

                      case (fbem_boundary_class_cracklike)
!~                         !
!~                         ! SBIE (Σt_k unknown)
!~                         !
!~                         ! node(sn_be)%col(problem%n+kdof_col,1): Σt_k = t_k^+ + t_k^-
!~                         ! q_k = - t_k^+ - t_k^-
!~                         t_e((problem%n*(kn-1)+1):(problem%n*kn))=node(sn_be)%value_r((problem%n+1):2*problem%n,1)
                        !
                        ! DBEM coupling (potential for general interface conditions and individual tractions)
                        !
                        ! node(sn_be)%col(problem%n+kdof_col,1): t_k^+
                        ! node(sn_be)%col(problem%n+kdof_col,2): t_k^-
                        ! q_k = - t_k^+ - t_k^-
                        t_e((problem%n*(kn-1)+1):(problem%n*kn))=node(sn_be)%value_r((problem%n+1):2*problem%n,1)&
                                                                +node(sn_be)%value_r((problem%n+1):2*problem%n,2)


                    end select

                  ! =================
                  ! BE-FE-BE BOUNDARY
                  ! =================

                  case (fbem_boundary_coupling_be_fe_be)
                    ! node(sn_be)%col(problem%n+kdof_col,1): t_k^(1) (region 1)
                    ! node(sn_be)%col(problem%n+kdof_col,2): t_k^(2) (region 2)
                    ! q_k = - t_k^(1) - t_k^(2)
                    t_e((problem%n*(kn-1)+1):(problem%n*kn))=node(sn_be)%value_r((problem%n+1):2*problem%n,1)&
                                                            +node(sn_be)%value_r((problem%n+1):2*problem%n,2)

                end select

              end do
              ! Build the FE load vector (f_e)
              f_e=matmul(Qmid,t_e)
              ! Save to q
              kdofe=1
              do kn=1,element(se_int)%n_nodes
                do kc=1,element(se_int)%node_n_dof(kn)
                  element(se_int)%value_r(kc,kn,2)=element(se_int)%value_r(kc,kn,2)+f_e(kdofe)
                  kdofe=kdofe+1
                end do
              end do
              ! Finalization
              deallocate (t_e,f_e)

            case default
              stop 'invalid mode for coupling assemble'
          end select

        ! ============
        ! BE BODY LOAD
        ! ============

        case (fbem_part_be_bodyload)
          !
          ! ASSEMBLE
          !
          select case (mode)
            !
            ! ASSEMBLE to Ax=b
            !
            case (0)
              ! Local row
              row_local=1
              do kn_row=1,element(se_int)%n_nodes
                sn_row=element(se_int)%node(kn_row)
                do kdof_row=1,element(se_int)%node_n_dof(kn_row)
                  if (node(sn_row)%ctype(kdof_row,1).eq.1) then
                    row=node(sn_row)%row(kdof_row,1)
                    ! Local column
                    col_local=1
                    do kn=1,element(se_int)%n_nodes
                      ! Selected FE node and its corresponding BE node
                      sn=element(se_int)%node(kn)
                      sn_be=element(se_int)%element_node(kn)
                      ! Loop through coordinates
                      do kdof_col=1,problem%n
                        col=node(sn_be)%col(problem%n+kdof_col,1)
                        A_r(row,col)=A_r(row,col)+Qmid(row_local,col_local)
                        col_local=col_local+1
                      end do
                    end do
                  end if
                  row_local=row_local+1
                end do
              end do

            !
            ! ASSEMBLE to q=Ka-f_bem-f_ext
            !
            case (1)

              ! Initialization
              allocate (t_e(ndof_load))
              allocate (f_e(ndof_se))
              f_e=0
              t_e=0
              ! Build the BE load vector (t_e)
              do kn=1,element(se_int)%n_nodes
                ! Selected FE node and its corresponding BE node
                sn=element(se_int)%node(kn)
                sn_be=element(se_int)%element_node(kn)
                t_e((problem%n*(kn-1)+1):(problem%n*kn))=node(sn_be)%value_r((problem%n+1):2*problem%n,1)
              end do
              ! Build the FE load vector (f_e)
              f_e=matmul(Qmid,t_e)
              ! Save to q
              kdofe=1
              do kn=1,element(se_int)%n_nodes
                do kc=1,element(se_int)%node_n_dof(kn)
                  element(se_int)%value_r(kc,kn,2)=element(se_int)%value_r(kc,kn,2)+f_e(kdofe)
                  kdofe=kdofe+1
                end do
              end do
              ! Finalization
              deallocate (t_e,f_e)

            case default
              stop 'invalid mode for coupling assemble'
          end select




      end select
      !
      ! Finalization
      !
      deallocate (Qmid)
    end if

  end if

  ! ========================================================
  ! (n-2)D STRUCTURAL FINITE ELEMENT (BEAM FE IN 3D AMBIENT)
  ! ========================================================

  if ((element(se_int)%n_dimension.eq.1).and.(problem%n.eq.3)) then

    if (element(se_int)%element.gt.0) then
      !
      ! Build Qmid
      !
      ! Initialize
      ndof_se=6*element(se_int)%n_nodes
      ndof_load=3*element(se_int)%n_nodes
      allocate (Qmid(ndof_se,ndof_load))
      ! Calculate
      select case (element(se_int)%fe_type)

        ! ------------------------
        ! DEGENERATED BEAM ELEMENT
        ! ------------------------

        case (0)
          call fbem_fem_degbeam_Q_midline(problem%n,element(se_int)%type,element(se_int)%x_gn,&
                                          element(se_int)%Q_intmode,element(se_int)%Q_intngp,Qmid)

        ! ---------------------------------------------------
        ! STRAIGHT EULER-BERNOULLI OR TIMOSHENKO BEAM ELEMENT
        ! ---------------------------------------------------

        case (1,2)
          !
          ! TO-DO: Los ejes nodales habría que guardarlos para cada nodo y construirlos para cada elemento
          !
          allocate(nodal_axes(problem%n,problem%n,element(se_int)%n_nodes))
          nodal_axes=0
          do kn=1,element(se_int)%n_nodes
            do kc=1,problem%n
              nodal_axes(kc,kc,kn)=1
            end do
          end do
          !
          ! TODO: Usar el paradigma de region()%material(1)
          ! Use material of the FE REGION
          !
          nu=region(sr)%property_r(3)
          !
          ! Calculate
          !
          call fbem_fem_strbeam_Q_midline(problem%n,element(se_int)%type,element(se_int)%fe_options(1),element(se_int)%fe_type,&
                                          element(se_int)%x_gn,element(se_int)%ep,&
                                          element(se_int)%A,element(se_int)%I(1:3),element(se_int)%ksh,&
                                          1.0d0,nu,nodal_axes,Qmid)

        case default
          ! AQUI FALTA METER LOS DISTRA, DISROTRA, Y BAR ELEMENTS (habría que ver qué y cómo)
          call fbem_error_message(error_unit,0,'element',element(se_int)%id,'this type of element cannot be coupled to BE body loads')

      end select
      !
      ! Assemble Qmid
      !
      ! FE -> BE CONNECTIVITY
      se_be=element(se_int)%element
      sb_be=part(element(se_be)%part)%entity
      select case (part(element(se_be)%part)%type)

        ! ===========
        ! BE BOUNDARY
        ! ===========

        case (fbem_part_be_boundary)
          stop 'Fatal error: this should not happen'

        ! ============
        ! BE BODY LOAD
        ! ============

        case (fbem_part_be_bodyload)
          select case (mode)
            !
            ! ASSEMBLE to Ax=b
            !
            case (0)
              ! Local row
              row_local=1
              do kn_row=1,element(se_int)%n_nodes
                sn_row=element(se_int)%node(kn_row)
                do kdof_row=1,node(sn_row)%n_dof
                  if (node(sn_row)%ctype(kdof_row,1).eq.1) then
                    row=node(sn_row)%row(kdof_row,1)
                    ! Local column
                    col_local=1
                    do kn=1,element(se_int)%n_nodes
                      ! Selected Shell FE node and its corresponding BE node
                      sn=element(se_int)%node(kn)
                      sn_be=element(se_int)%element_node(kn)
                      ! Loop through coordinates
                      do kdof_col=1,problem%n
                        col=node(sn_be)%col(problem%n+kdof_col,1)
                        A_r(row,col)=A_r(row,col)+Qmid(row_local,col_local)
                        col_local=col_local+1
                      end do
                    end do
                  end if
                  row_local=row_local+1
                end do
              end do
            !
            ! ASSEMBLE to q=Ka-f_bem-f_ext
            !
            case (1)
              ! Initialization
              allocate (t_e(ndof_load))
              allocate (f_e(ndof_se))
              f_e=0
              t_e=0
              ! Build the BE load vector (t_e)
              do kn=1,element(se_int)%n_nodes
                ! Selected FE node and its corresponding BE node
                sn=element(se_int)%node(kn)
                sn_be=element(se_int)%element_node(kn)
                t_e((problem%n*(kn-1)+1):(problem%n*kn))=node(sn_be)%value_r((problem%n+1):2*problem%n,1)
              end do
              ! Build the FE load vector (f_e)
              f_e=matmul(Qmid,t_e)
              ! Save to q
              kdofe=1
              do kn=1,element(se_int)%n_nodes
                do kc=1,element(se_int)%node_n_dof(kn)
                  element(se_int)%value_r(kc,kn,2)=element(se_int)%value_r(kc,kn,2)+f_e(kdofe)
                  kdofe=kdofe+1
                end do
              end do
              ! Finalization
              deallocate (t_e,f_e)

            case default
              stop 'invalid mode for coupling assemble'

          end select

      end select

      deallocate (Qmid)

    end if
  end if

end subroutine assemble_fem_staela_coupling
