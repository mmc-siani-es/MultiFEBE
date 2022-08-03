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


subroutine build_data_of_be_elements

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_geometry
  use fbem_quasisingular_integration
  use fbem_shape_functions

  ! Module of problem variables
  use problem_variables

  ! Local variables
  implicit none
  integer                        :: ke, se, sp
  integer                        :: se_n_nodes
  integer                        :: se_type_g, se_type_f1, se_type_f2
  integer                        :: n_phi_f1, n_phi_f2
  real(kind=real64), allocatable :: x_gn(:,:)
  real(kind=real64)              :: csize, centre(problem%n), radius
  character(len=fbem_fmtstr)     :: fmtstr
  logical                        :: ok

  ! Starting message
  if (verbose_level.ge.1)  write(output_unit,'(a)') 'Building data of boundary elements ...'

  !
  ! Loop through BE ELEMENTS and BE BODY LOAD ELEMENTS
  !
  do ke=1,n_elements
    ! Selected element
    se=ke
    sp=element(se)%part
    if ((part(sp)%type.eq.fbem_part_be_boundary).or.(part(sp)%type.eq.fbem_part_be_bodyload)) then
      se_type_g=element(se)%type_g
      se_type_f1=element(se)%type_f1
      se_type_f2=element(se)%type_f2
      se_n_nodes=element(se)%n_nodes
      allocate (x_gn(problem%n,fbem_n_nodes(se_type_g)))
      ! Coordinates of each (geometrical) node of the element
      x_gn=element(se)%x_gn
      ! Size of the element
      element(se)%size=fbem_element_size(problem%n,se_type_g,x_gn,1.d-9)
      ! Characteristic size of the element
      csize=fbem_characteristic_length(problem%n,se_type_g,x_gn,1.d-9)
      ! Estimate number of 1D Gauss-Legendre quadrature points to integrate any shape function integration
      select case (part(sp)%type)
        case (fbem_part_be_boundary)
          select case (fbem_n_dimension(se_type_g))
            case (1)
              call fbem_qs_phijac_ngp_1d(problem%n,se_type_g,se_type_f1,x_gn,qsi_relative_error,n_phi_f1,ok)
              call fbem_qs_phijac_ngp_1d(problem%n,se_type_g,se_type_f2,x_gn,qsi_relative_error,n_phi_f2,ok)
            case (2)
              call fbem_qs_phijac_ngp_2d(problem%n,se_type_g,se_type_f1,x_gn,qsi_relative_error,n_phi_f1,ok)
              call fbem_qs_phijac_ngp_2d(problem%n,se_type_g,se_type_f2,x_gn,qsi_relative_error,n_phi_f2,ok)
          end select
        case (fbem_part_be_bodyload)
          n_phi_f1=0
          select case (fbem_n_dimension(se_type_g))
            case (1)
              call fbem_qs_phijac_ngp_1d(problem%n,se_type_g,se_type_f2,x_gn,qsi_relative_error,n_phi_f2,ok)
            case (2)
              call fbem_qs_phijac_ngp_2d(problem%n,se_type_g,se_type_f2,x_gn,qsi_relative_error,n_phi_f2,ok)
          end select
      end select
      ! If ok is false, then it means that even with the maximum number of integration points, the integration cannot
      ! be performed with the specified relative error, probably because the element is too distorted or invalid.
      if (ok.eqv.(.false.)) then
        write(error_unit,'(a,i3)') 'N_phif1 =',n_phi_f1
        write(error_unit,'(a,i3)') 'N_phif2 =',n_phi_f2
        call fbem_warning_message(error_unit,0,'element',element(se)%id,'requires too many integration points, check its validity.')
      end if
      ! Save to the element data structure
      element(se)%csize=csize
      element(se)%n_phi=max(n_phi_f1,n_phi_f2)
      ! Calculate the centroid
      allocate (element(se)%centroid(problem%n))
      element(se)%centroid=fbem_element_centroid(problem%n,se_type_g,x_gn,element(se)%size,1.d-9)
      ! Calculate a ball that contains the element
      call fbem_geometry_element_ball(problem%n,se_type_g,x_gn,element(se)%n_phi,centre,radius)
      ! Save to the element data structure
      allocate (element(se)%bball_centre(problem%n))
      element(se)%bball_centre=centre
      element(se)%bball_radius=radius
      !  Message
      if (verbose_level.ge.2) then
        write(fmtstr,*) '(a,i',fbem_nchar_int(element(se)%id),')'
        call fbem_trimall(fmtstr)
        write(output_unit,fmtstr)         ' Element: ', element(se)%id
        if (verbose_level.ge.3) then
          write(output_unit,'(a,e25.16)') '  Characteristic size: ', element(se)%csize
          write(fmtstr,*) '(a,i',fbem_nchar_int(element(se)%n_phi),')'
          call fbem_trimall(fmtstr)
          write(output_unit,fmtstr)       '  1D Gauss-Legendre N to integrate shape functions: ', element(se)%n_phi
        end if
      end if
      deallocate (x_gn)
    end if
  end do

  ! Ending message
  if (verbose_level.ge.1) write(output_unit,'(a)') 'done.'

end subroutine build_data_of_be_elements
