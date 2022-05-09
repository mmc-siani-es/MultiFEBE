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

subroutine assemble_fem_harela_dKdx(se,se_n_dof,dKdx)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! Problem variables module
  use problem_variables

  ! No implicit variables are allowed
  implicit none

  ! I/O variables
  integer              :: se                                                    !! Element to assemble
  integer              :: se_n_dof                                              !! Number of degrees of freedom of the element
  complex(kind=real64) :: dKdx(se_n_dof,se_n_dof,problem%n,element(se)%n_nodes) !! dKdx matrix

  ! Local variables
  integer              :: se_n_nodes              ! Number of nodes of the element
  integer              :: row                     ! Row of the LSE
  integer              :: row_local               ! Row of the element matrix
  integer              :: col                     ! Column of the LSE
  integer              :: col_local               ! Column of the element matrix
  integer              :: knr, snr, kcr           ! Node counter, selected node, and DOF counter for row loop
  integer              :: knc, snc, kcc           ! Node counter, selected node, and DOF counter for column loop
  integer              :: ka, kn, sn, kc
  complex(kind=real64) :: dKda(se_n_dof,se_n_dof)

  ! Initialize
  se_n_nodes=element(se)%n_nodes

  ! Loop through design variables
  do ka=1,problem%n_designvariables

    ! Build dK/da=dK/dx_k^(i)*dx_k^(i)/da
    dKda=(0.d0,0.d0)
    do kn=1,se_n_nodes
      sn=element(se)%node(kn)
      do kc=1,problem%n
        dKda=dKda+dKdx(:,:,kc,kn)*node(sn)%dxda(kc,ka)
      end do
    end do
    ! Assemble dKda
    row_local=1
    ! Loop through the NODES of the ELEMENT to assemble the EQUATION of each DOF
    do knr=1,se_n_nodes
      snr=element(se)%node(knr)
      ! Loop through the DOFs of the node
      do kcr=1,node(snr)%n_dof
        ! The equation is ACTIVE is the associated DOF is ACTIVE
        if (node(snr)%ctype(kcr,1).eq.1) then
          ! ROW of the EQUATION in the LSE
          row=node(snr)%row(kcr,1)
          ! Initialize
          col_local=1
          ! Loop through the NODES of the ELEMENT to assemble dKda for each node
          do knc=1,se_n_nodes
            snc=element(se)%node(knc)
            ! Loop through the DOFs of the node
            do kcc=1,node(snc)%n_dof
              bsa_c(row,ka)=bsa_c(row,ka)-dKda(row_local,col_local)*node(snc)%value_c(kcc,1)
              col_local=col_local+1
            end do
          end do
        end if
        row_local=row_local+1
      end do
    end do

  end do ! Loop through design variables

end subroutine assemble_fem_harela_dKdx
