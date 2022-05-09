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

subroutine assemble_fem_harela_dfdx(se_int)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! Problem variables module
  use problem_variables
  use fbem_shape_functions
  use fbem_fem_solids
  use fbem_fem_beams
  use fbem_fem_shells

  ! No implicit variables are allowed
  implicit none

  ! I/O variables
  integer                        :: se_int                           !! Selected finite element

  ! Local variables
  integer                        :: ka
  integer                        :: kn_row, sn_row, kn, sn, kc
  integer                        :: kdof_row, kdof_col
  integer                        :: row
  integer                        :: row_local
  integer                        :: col
  integer                        :: col_local
  integer                        :: se_be, sb_be, kn_be, sn_be
  integer                        :: kd_fe, sd_fe, sd_type, sd_fe_n_nodes
  integer                        :: sd_be
  integer                        :: se_be_type, se_be_stype, se_be_n_nodes, se_be_n_snodes
  integer                        :: knt
  real(kind=real64)              :: length
  real(kind=real64), allocatable :: dQdx(:,:,:,:), dQda(:,:)
  integer                        :: se_int_type, se_int_n_nodes, se_int_n_dof
  integer, allocatable           :: se_int_n_dof_node(:)
  real(kind=real64), allocatable :: se_int_x(:,:), se_int_v(:,:,:), se_int_t(:)

  ! ===================================
  ! nD CONTINUUM (SOLID) FINITE ELEMENT
  ! ===================================

  if (element(se_int)%n_dimension.eq.problem%n) then
    select case (problem%n)

      ! -----------
      ! 2D SOLID FE
      ! -----------

      case (2)
        !
        ! BEM-FEM COUPLING
        !
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
            sb_be=part(element(se_be)%part)%entity
            se_be_type=element(se_be)%type
            se_be_n_nodes=element(se_be)%n_nodes
            ! CALCULATE Qedge
            allocate (dQdx(2*se_int_n_nodes,2*se_be_n_nodes,2,se_int_n_nodes),dQda(2*se_int_n_nodes,2*se_be_n_nodes))
            call fbem_fem_ela2d_solid_dQdx_edge(se_int_type,se_int_x,se_int_t,kd_fe,se_be_type,4,dQdx)
            ! ASSEMBLE
            do ka=1,problem%n_designvariables
              ! Build dQ/da=dQ/dx_k^(i)*dx_k^(i)/da
              dQda=(0.d0,0.d0)
              do kn=1,se_int_n_nodes
                sn=element(se_int)%node(kn)
                do kc=1,problem%n
                  dQda=dQda+dQdx(:,:,kc,kn)*node(sn)%dxda(kc,ka)
                end do
              end do
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
                    ! Loop through FE edge nodes
                    do kn=1,sd_fe_n_nodes
                      ! Note that Qedge is built using the element parent orientation for the edge. Because of that the reversion
                      ! is needed in general. Although here if the edge is a boundary edge (it should be for coupling), the subedge
                      ! must have the same orientation
                      if (element(se_int)%subedge_reversion(kd_fe)) then
                        sn_be=subedge(sd_fe)%element_node(fbem_element_invert(kn,sd_type))
                      else
                        sn_be=subedge(sd_fe)%element_node(kn)
                      end if
                      ! Save the index of the BE node in the BE element
                      do kn_be=1,element(se_be)%n_nodes
                        if (element(se_be)%node(kn_be).eq.sn_be) exit
                      end do
                      do kdof_col=1,problem%n
                        select case (region(boundary(sb_be)%region(1))%type)

                          ! --------------
                          ! INVISCID FLUID
                          ! --------------

                          case (fbem_potential)
                            ! node(sn_be)%col(1,1): p
                            ! q_k = p·n_k
                            if (boundary(sb_be)%region_boundary_reversion(1)) then
                              bsa_c(row,ka)=bsa_c(row,ka)-dQda(row_local,col_local)*node(sn_be)%value_c(1,1)*element(se_be)%n_fn(kdof_col,kn_be)
                            else
                              bsa_c(row,ka)=bsa_c(row,ka)+dQda(row_local,col_local)*node(sn_be)%value_c(1,1)*element(se_be)%n_fn(kdof_col,kn_be)
                            end if

                          ! ------------------
                          ! VISCOELASTIC SOLID
                          ! ------------------

                          case (fbem_viscoelastic)
                            ! node(sn_be)%col(problem%n+kdof_col,1): t_k
                            ! q_k = - t_k
                            bsa_c(row,ka)=bsa_c(row,ka)-dQda(row_local,col_local)*node(sn_be)%value_c(problem%n+kdof_col,1)

                          ! ------------------
                          ! POROELASTIC MEDIUM
                          ! ------------------

                          case (fbem_poroelastic)
                            stop 'stop: not harpor_sa with perturbed boundary'

                        end select
                        col_local=col_local+1
                      end do
                    end do
                  end if
                  row_local=row_local+1
                end do
              end do
            end do
            deallocate (se_int_x,se_int_t,se_int_n_dof_node,dQdx,dQda)
          end if
        end do

      ! -----------
      ! 3D SOLID FE
      ! -----------

      case (3)
        stop 'not yet 22' ! Debería ser lo mismo que 2D pero con subfaces....

    end select
  end if

  ! etc.... copiar y adaptar de assemble_fem_harela_coupling

end subroutine assemble_fem_harela_dfdx
