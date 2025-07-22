! ---------------------------------------------------------------------
! Copyright (C) 2014-2024 Universidad de Las Palmas de Gran Canaria:
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

subroutine transfer_conditions_bem_bodyloads_mechanics_harmonic

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_string_handling
  use fbem_geometry
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! No implicit variables are allowed
  implicit none

  ! Local variables
  integer           :: kb, sp
  integer           :: kn, sn

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START transferring conditions (BE bodyloads -> BE bodyload nodes)')

  ! Loop through BE BODY LOADS
  do kb=1,n_be_bodyloads
    sp=be_bodyload(kb)%part
    select case (be_bodyload(kb)%coupling)

      case (fbem_bl_uncoupled)
        select case (region(be_bodyload(kb)%region)%type)

          ! --------------
          ! INVISCID FLUID
          ! --------------

          case (fbem_potential)
            do kn=1,part(sp)%n_nodes
              sn=part(sp)%node(kn)
              allocate (node(sn)%ctype(1,1))
              allocate (node(sn)%cvalue_c(1,1,1))
              node(sn)%ctype=be_bodyload(kb)%ctype
              node(sn)%cvalue_c=be_bodyload(kb)%cvalue_c
            end do

          ! ------------------
          ! VISCOELASTIC SOLID
          ! ------------------

          case (fbem_viscoelastic)
            do kn=1,part(sp)%n_nodes
              sn=part(sp)%node(kn)
              allocate (node(sn)%ctype(problem%n,1))
              allocate (node(sn)%cvalue_c(problem%n,1,1))
              node(sn)%ctype=be_bodyload(kb)%ctype
              node(sn)%cvalue_c=be_bodyload(kb)%cvalue_c
            end do

          ! -----------------
          ! POROELASTIC MEDIUM
          ! -----------------

          case (fbem_poroelastic)
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,'body loads not available for poroelastic media')

        end select

      case (fbem_bl_coupling_beam_tip)
        ! N/A

      case (fbem_bl_coupling_beam_line)
        ! N/A

      case (fbem_bl_coupling_shell_edge)
        ! N/A

      case (fbem_bl_coupling_shell_surface)
        ! N/A
    end select

  end do ! Loop through BE BODY LOADS

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END transferring conditions (BE bodyloads -> BE bodyload nodes)')

end subroutine transfer_conditions_bem_bodyloads_mechanics_harmonic
