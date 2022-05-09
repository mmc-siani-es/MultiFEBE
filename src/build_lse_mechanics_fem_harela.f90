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

subroutine build_lse_mechanics_fem_harela(kf,kr)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_numerical
  use fbem_shape_functions
  use fbem_geometry
  use fbem_symmetry
  use fbem_quasisingular_integration
  use fbem_telles_transformation
  use fbem_bem_general
  use fbem_bem_staela2d
  use fbem_bem_staela3d
  use fbem_fem_beams
  use fbem_fem_shells
  use fbem_fem_solids

  ! Module of problem variables
  use problem_variables

  ! No implicit variables
  implicit none

  ! I/O variables
  integer                           :: kf
  integer                           :: kr
  ! Local variables
  integer                           :: ks_int, ss_int
  integer                           :: ke_int, se_int
  integer                           :: kn, sn
  integer                           :: kc
  integer                           :: ndof_mat
  complex(kind=real64), allocatable :: K(:,:), f_e(:)
  integer                           :: row, kdofe
  character(len=fbem_fmtstr)        :: fmtstr

  if (verbose_level.ge.2) then
    write(fmtstr,*) '(a26,1x,i',fbem_nchar_int(region(kr)%id),')'
    call fbem_trimall(fmtstr)
    call fbem_timestamp_message(output_unit,2)
    write(output_unit,fmtstr) 'START assembling FE region',region(kr)%id
  end if




  !
  ! Ojo, que en el acople BEM-FEM si habria que ensamblar las Q ?¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿?????????????????????????
  !
  if (region(kr)%type.eq.fbem_rigid) return





  ! Loop through the FE SUBREGIONS of the FE REGION
  do ks_int=1,region(kr)%n_fe_subregions
    ss_int=region(kr)%fe_subregion(ks_int)

    ! Loop through the FINITE ELEMENTS of the FE SUBREGION
    do ke_int=1,part(fe_subregion(ss_int)%part)%n_elements
      se_int=part(fe_subregion(ss_int)%part)%element(ke_int)

      if (verbose_level.ge.3) then
        write(fmtstr,*) '(a24,1x,i',fbem_nchar_int(element(se_int)%id),')'
        call fbem_trimall(fmtstr)
        call fbem_timestamp_message(output_unit,2)
        write(output_unit,fmtstr) 'START assembling element',element(se_int)%id
      end if

      !
      ! Build and assemble the stiffness matrix
      !
      ndof_mat=3*(problem%n-1)*element(se_int)%n_nodes
      allocate (K(ndof_mat,ndof_mat))
      call build_fem_K_harmonic(kf,se_int,ndof_mat,K)
      call assemble_fem_harela_stiffness_matrix(se_int,ndof_mat,K)
      deallocate(K)

      !
      ! Build and assemble coupling terms of finite elements
      !
      ! TODO: Add -f_bem to q (q=K·a-f_bem-f_ext)
      !
      call assemble_fem_harela_coupling(se_int,0)

      !
      ! Build and assemble distributed loads
      !
      ! TODO: Pasar al bucle de ensamblaje de cargas distribuidas independiente o no?
      !
      ! Initialize and build f_e
      ndof_mat=3*(problem%n-1)*element(se_int)%n_nodes
      allocate (f_e(ndof_mat))
      call build_fem_f_distributed_harmonic(kf,se_int,ndof_mat,f_e)

      !
      ! TODO: ojo no preparado para regiones rigidas (el ensamblaje)¿¿¿¿¿¿¿¿¿¿¿?????????????????????????????????????????????
      !

      !
      ! Assemble f_e to b (A·x=b)
      !
      kdofe=1
      do kn=1,element(se_int)%n_nodes
        sn=element(se_int)%node(kn)
        do kc=1,element(se_int)%node_n_dof(kn)
          if (node(sn)%ctype(kc,1).eq.1) then
            row=node(sn)%row(kc,1)
            b_c(row,1)=b_c(row,1)+f_e(kdofe)
          end if
          kdofe=kdofe+1
        end do
      end do
      deallocate(f_e)

    end do ! Loop through the FINITE ELEMENTS of the FE SUBREGION

  end do ! Loop through the FE SUBREGION of the FE REGION

  if (verbose_level.ge.2) then
    write(fmtstr,*) '(a24,1x,i',fbem_nchar_int(region(kr)%id),')'
    call fbem_trimall(fmtstr)
    call fbem_timestamp_message(output_unit,2)
    write(output_unit,fmtstr) 'END assembling FE region',region(kr)%id
  end if

end subroutine build_lse_mechanics_fem_harela
