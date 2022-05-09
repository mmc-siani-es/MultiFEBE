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

subroutine read_special_boundary_elements(input_fileunit)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_symmetry
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! No implicit variables are allowed
  implicit none

  ! I/O
  integer                        :: input_fileunit              !! Input file unit
  ! Local
  character(len=fbem_stdcharlen) :: section_name                ! Name of the section
  logical                        :: found                       ! Logical variable for sections and keywords
  integer                        :: n_special_elements          !
  character(len=fbem_stdcharlen) :: tmp_type                    ! Temporary variable to read the type of the element
  integer                        :: ke, se_eid, se

  ! Return if not needed
  if (n_be_regions.eq.0) return

  ! Locate the section
  section_name='special be elements'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(input_fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! Read the number of special elements to be set
    read(input_fileunit,*) n_special_elements
    ! Read each element
    do ke=1,n_special_elements
      read(input_fileunit,*) se_eid, tmp_type
      se=element_iid(se_eid)
      ! Switch depending on the general type of element
      select case (element(se)%type)
        case (fbem_line3)
          element(se)%type_f1 =0
          element(se)%type_f2 =0
          element(se)%type_f1f=0
          element(se)%type_f2f=0
          if (trim(tmp_type).eq.'line3_qp1_normal') then
            element(se)%type_f1=fbem_line3
            element(se)%type_f2=fbem_line3
            element(se)%type_f1f=fbem_line3
            element(se)%type_f2f=fbem_line3
            call transformation_element_to_quarter_point(element(se)%node(1))
          end if
          if (trim(tmp_type).eq.'line3_qp2_normal') then
            element(se)%type_f1=fbem_line3
            element(se)%type_f2=fbem_line3
            element(se)%type_f1f=fbem_line3
            element(se)%type_f2f=fbem_line3
            call transformation_element_to_quarter_point(element(se)%node(2))
          end if
          if (trim(tmp_type).eq.'line3_mqp1_normal') then
            element(se)%type_f1=fbem_line3_mqp1u
            element(se)%type_f2=fbem_line3
            element(se)%type_f1f=fbem_line3_mqp1u
            element(se)%type_f2f=fbem_line3
            call transformation_element_to_quarter_point(element(se)%node(1))
          end if
          if (trim(tmp_type).eq.'line3_mqp2_normal') then
            element(se)%type_f1=fbem_line3_mqp2u
            element(se)%type_f2=fbem_line3
            element(se)%type_f1f=fbem_line3_mqp2u
            element(se)%type_f2f=fbem_line3
            call transformation_element_to_quarter_point(element(se)%node(2))
          end if
          if (trim(tmp_type).eq.'line3_qp1') then
            element(se)%type_f1=fbem_line3
            element(se)%type_f2=fbem_line3_qp1t
            element(se)%type_f1f=fbem_line3
            element(se)%type_f2f=fbem_line3_qp1t
            call transformation_element_to_quarter_point(element(se)%node(1))
          end if
          if (trim(tmp_type).eq.'line3_qp2') then
            element(se)%type_f1=fbem_line3
            element(se)%type_f2=fbem_line3_qp2t
            element(se)%type_f1f=fbem_line3
            element(se)%type_f2f=fbem_line3_qp2t
            call transformation_element_to_quarter_point(element(se)%node(2))
          end if
          if (trim(tmp_type).eq.'line3_mqp1') then
            element(se)%type_f1=fbem_line3_mqp1u
            element(se)%type_f2=fbem_line3_mqp1t
            element(se)%type_f1f=fbem_line3_mqp1u
            element(se)%type_f2f=fbem_line3_mqp1t
            call transformation_element_to_quarter_point(element(se)%node(1))
          end if
          if (trim(tmp_type).eq.'line3_mqp2') then
            element(se)%type_f1=fbem_line3_mqp2u
            element(se)%type_f2=fbem_line3_mqp2t
            element(se)%type_f1f=fbem_line3_mqp2u
            element(se)%type_f2f=fbem_line3_mqp2t
            call transformation_element_to_quarter_point(element(se)%node(2))
          end if
          !
          ! Quarter-point element configurations for poroelastic crack problems
          !
          !
          ! All variables have sqrt(r) behaviour, except t_k and Un
          if (trim(tmp_type).eq.'line3_qp1_poro_conf2') then
            element(se)%type_f1 =fbem_line3_mqp1u
            element(se)%type_f2 =fbem_line3_mqp1t
            element(se)%type_f1f=fbem_line3_mqp1u
            element(se)%type_f2f=fbem_line3_mqp1t
            call transformation_element_to_quarter_point(element(se)%node(1))
          end if
          if (trim(tmp_type).eq.'line3_qp2_poro_conf2') then
            element(se)%type_f1 =fbem_line3_mqp2u
            element(se)%type_f2 =fbem_line3_mqp2t
            element(se)%type_f1f=fbem_line3_mqp2u
            element(se)%type_f2f=fbem_line3_mqp2t
            call transformation_element_to_quarter_point(element(se)%node(2))
          end if
          ! All variables have sqrt(r) behaviour, except tau
          if (trim(tmp_type).eq.'line3_qp1_poro_conf3') then
            element(se)%type_f1 =fbem_line3_mqp1u
            element(se)%type_f2 =fbem_line3_mqp1u
            element(se)%type_f1f=fbem_line3_mqp1t
            element(se)%type_f2f=fbem_line3_mqp1u
            call transformation_element_to_quarter_point(element(se)%node(1))
          end if
          if (trim(tmp_type).eq.'line3_qp2_poro_conf3') then
            element(se)%type_f1 =fbem_line3_mqp2u
            element(se)%type_f2 =fbem_line3_mqp2u
            element(se)%type_f1f=fbem_line3_mqp2t
            element(se)%type_f2f=fbem_line3_mqp2u
            call transformation_element_to_quarter_point(element(se)%node(2))
          end if
          ! All variables have sqrt(r) behaviour, except t_k
          if (trim(tmp_type).eq.'line3_qp1_poro_conf4') then
            element(se)%type_f1 =fbem_line3_mqp1u
            element(se)%type_f2 =fbem_line3_mqp1t
            element(se)%type_f1f=fbem_line3_mqp1u
            element(se)%type_f2f=fbem_line3_mqp1u
            call transformation_element_to_quarter_point(element(se)%node(1))
          end if
          if (trim(tmp_type).eq.'line3_qp2_poro_conf4') then
            element(se)%type_f1 =fbem_line3_mqp2u
            element(se)%type_f2 =fbem_line3_mqp2t
            element(se)%type_f1f=fbem_line3_mqp2u
            element(se)%type_f2f=fbem_line3_mqp2u
            call transformation_element_to_quarter_point(element(se)%node(2))
          end if
          ! All variables have sqrt(r) behaviour, except tau and t_k
          if (trim(tmp_type).eq.'line3_qp1_poro_conf5') then
            element(se)%type_f1 =fbem_line3_mqp1u
            element(se)%type_f2 =fbem_line3_mqp1t
            element(se)%type_f1f=fbem_line3_mqp1t
            element(se)%type_f2f=fbem_line3_mqp1u
            call transformation_element_to_quarter_point(element(se)%node(1))
          end if
          if (trim(tmp_type).eq.'line3_qp2_poro_conf5') then
            element(se)%type_f1 =fbem_line3_mqp2u
            element(se)%type_f2 =fbem_line3_mqp2t
            element(se)%type_f1f=fbem_line3_mqp2t
            element(se)%type_f2f=fbem_line3_mqp2u
            call transformation_element_to_quarter_point(element(se)%node(2))
          end if
          if (element(se)%type_f1.eq.0) then
              call fbem_error_message(error_unit,0,'element',element(se)%id,'the special element type does not exist')
          end if
        case default
          call fbem_error_message(error_unit,0,'element',element(se)%id,'this type of element does not admit any speciality')
      end select

    end do

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  end if


end subroutine read_special_boundary_elements
