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

subroutine assemble_fem_staela_stiffness_matrix(se,se_n_dof,K)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! Problem variables module
  use problem_variables
  use fbem_numerical
  use fbem_geometry

  ! No implicit variables are allowed
  implicit none

  ! I/O variables
  integer                        :: se                   !! Element to assemble
  integer                        :: se_n_dof             !! Number of degrees of freedom of the element
  real(kind=real64)              :: K(se_n_dof,se_n_dof) !! Stiffness matrix

  ! Local variables
  integer                        :: se_n_nodes           ! Number of nodes of the element
  integer                        :: row                  ! Row of the LSE
  integer                        :: row_local            ! Row of the element matrix
  integer                        :: col                  ! Column of the LSE
  integer                        :: col_local            ! Column of the element matrix
  integer                        :: row_node             ! Row of the element matrix
  integer                        :: col_node             ! Column of the element matrix
  integer                        :: knr, snr, kcr        ! Node counter, selected node, and DOF counter for row loop
  integer                        :: knc, snc, kcc        ! Node counter, selected node, and DOF counter for column loop
  integer                        :: snrm, sncm
  real(kind=real64), allocatable :: TR(:,:), TC(:,:)
  real(kind=real64), allocatable :: Kt(:,:)

  ! Initialize
  se_n_nodes=element(se)%n_nodes
  row_node=1
  ! Loop through the NODES of the ELEMENT to assemble the EQUATION of each degree of freedom
  do knr=1,se_n_nodes
    snr=element(se)%node(knr)
    ! Initialize
    col_node=1
    ! Loop through the NODES of the ELEMENT to assemble the STIFFNESS of each DOF
    do knc=1,se_n_nodes
      snc=element(se)%node(knc)
      select case (node(snr)%rigid_link)
        !
        ! (ROW) Normal or master node
        !
        case (0,1)
          select case (node(snc)%rigid_link)
            !
            ! (COL) Normal or master node
            !
            case (0,1)
              do kcr=1,element(se)%node_n_dof(knr)
                row_local=row_node+kcr-1
                if (node(snr)%ctype(kcr,1).eq.1) then
                  row=node(snr)%row(kcr,1)
                  do kcc=1,element(se)%node_n_dof(knc)
                    col_local=col_node+kcc-1
                    if (node(snc)%ctype(kcc,1).eq.1) then
                      col=node(snc)%col(kcc,1)
                      A_r(row,col)=A_r(row,col)+K(row_local,col_local)
                    else
                      b_r(row,1)=b_r(row,1)-K(row_local,col_local)*node(snc)%cvalue_r(kcc,1,1)
                    end if
                  end do
                end if
              end do
            !
            ! (COL) Slave node
            !
            case (2)
              sncm=node(snc)%master
              allocate (TC(element(se)%node_n_dof(knc),node(sncm)%n_dof))
              call fbem_rigid_solid_transformation_matrix(problem%n,node(snc)%x,node(sncm)%x,element(se)%node_n_dof(knc),TC)
              allocate (Kt(element(se)%node_n_dof(knr),node(sncm)%n_dof))
              Kt=matmul(K(row_node:(row_node+element(se)%node_n_dof(knr)-1),col_node:(col_node+element(se)%node_n_dof(knc)-1)),TC)
              do kcr=1,element(se)%node_n_dof(knr)
                if (node(snr)%ctype(kcr,1).eq.1) then
                  row=node(snr)%row(kcr,1)
                  do kcc=1,node(sncm)%n_dof
                    if (node(sncm)%ctype(kcc,1).eq.1) then
                      col=node(sncm)%col(kcc,1)
                      A_r(row,col)=A_r(row,col)+Kt(kcr,kcc)
                    else
                      b_r(row,1)=b_r(row,1)-Kt(kcr,kcc)*node(sncm)%cvalue_r(kcc,1,1)
                    end if
                  end do
                end if
              end do
              deallocate (Kt,TC)
          end select

        !
        ! (ROW) Slave node
        !
        case (2)
          select case (node(snc)%rigid_link)
            !
            ! (COL) Normal or master node
            !
            case (0,1)
              snrm=node(snr)%master
              allocate (TR(element(se)%node_n_dof(knr),node(snrm)%n_dof))
              call fbem_rigid_solid_transformation_matrix(problem%n,node(snr)%x,node(snrm)%x,element(se)%node_n_dof(knr),TR)
              allocate (Kt(node(snrm)%n_dof,element(se)%node_n_dof(knc)))
              Kt=matmul(transpose(TR),K(row_node:(row_node+element(se)%node_n_dof(knr)-1),col_node:(col_node+element(se)%node_n_dof(knc)-1)))
              do kcr=1,node(snrm)%n_dof
                if (node(snrm)%ctype(kcr,1).eq.1) then
                  row=node(snrm)%row(kcr,1)
                  do kcc=1,element(se)%node_n_dof(knc)
                    if (node(snc)%ctype(kcc,1).eq.1) then
                      col=node(snc)%col(kcc,1)
                      A_r(row,col)=A_r(row,col)+Kt(kcr,kcc)
                    else
                      b_r(row,1)=b_r(row,1)-Kt(kcr,kcc)*node(snc)%cvalue_r(kcc,1,1)
                    end if
                  end do
                end if
              end do
              deallocate (Kt,TR)
            !
            ! (COL) Slave node
            !
            case (2)
              snrm=node(snr)%master
              sncm=node(snc)%master
              allocate (TR(element(se)%node_n_dof(knr),node(snrm)%n_dof))
              allocate (TC(element(se)%node_n_dof(knc),node(sncm)%n_dof))
              call fbem_rigid_solid_transformation_matrix(problem%n,node(snr)%x,node(snrm)%x,element(se)%node_n_dof(knr),TR)
              call fbem_rigid_solid_transformation_matrix(problem%n,node(snc)%x,node(sncm)%x,element(se)%node_n_dof(knc),TC)
              allocate (Kt(node(snrm)%n_dof,node(sncm)%n_dof))
              Kt=matmul(matmul(transpose(TR),K(row_node:(row_node+element(se)%node_n_dof(knr)-1),col_node:(col_node+element(se)%node_n_dof(knc)-1))),TC)
              do kcr=1,node(snrm)%n_dof
                if (node(snrm)%ctype(kcr,1).eq.1) then
                  row=node(snrm)%row(kcr,1)
                  do kcc=1,node(sncm)%n_dof
                    if (node(sncm)%ctype(kcc,1).eq.1) then
                      col=node(sncm)%col(kcc,1)
                      A_r(row,col)=A_r(row,col)+Kt(kcr,kcc)
                    else
                      b_r(row,1)=b_r(row,1)-Kt(kcr,kcc)*node(sncm)%cvalue_r(kcc,1,1)
                    end if
                  end do
                end if
              end do
              deallocate (Kt,TC,TR)

          end select

      end select
      col_node=col_node+element(se)%node_n_dof(knc)
    end do
    row_node=row_node+element(se)%node_n_dof(knr)
  end do

end subroutine assemble_fem_staela_stiffness_matrix
