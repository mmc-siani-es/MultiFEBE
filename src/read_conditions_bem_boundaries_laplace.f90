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


subroutine read_conditions_bem_boundaries_laplace(input_fileunit)

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

  ! I/O variables
  integer                                :: input_fileunit    !! Input file unit
  ! Local variables
  character(len=fbem_stdcharlen)         :: section_name      ! Name of the section
  logical                                :: found
  integer                                :: i, j, k
  integer                                :: sn
  character(len=fbem_stdcharlen)         :: keyword
  character(len=fbem_file_record_length) :: line          ! Line
  character(len=fbem_file_record_length) :: word          ! Line
  integer                                :: nw, nc, nl
  integer                                :: ios
  integer                                :: n
  integer                                :: eid, sb
  real(kind=real64)                      :: tmp_value
  logical                                :: bc

!  ! Section name
!  section_name='conditions over be boundaries'
!  ! Initialize line counter
!  nl=1
!  do
!    ! Read the line and exit if the end of the file
!    read(input_fileunit,'(a)',iostat=ios) line
!    if (is_iostat_end(ios)) exit
!    ! Remove beginning, ending and duplicated blanks of a string
!    call fbem_trim2b(line)
!    ! Count the number of characters
!    nc=len_trim(line)
!    ! If the line is not empty, continue processing the line
!    if (nc.gt.0) then
!      ! If the line starts with "[", then it is the end of the section
!      if (line(1:1).eq.'[') exit
!      ! If the line starts with "#", then it is a comment
!      if (line(1:1).eq.'#') cycle
!      ! Count the number of words
!      nw=fbem_count_words(line)
!      ! Check for errors
!      if (nw.lt.4) then
!        call fbem_error_message(error_unit,0,section_name,nl,'not enough arguments in this line.')
!      end if
!      !
!      ! The first word must be "boundary"
!      !
!      word=fbem_extract_word(line,1)
!      if (trim(word).eq.'boundary') then
!        ! The second word must be a valid boundary id
!        word=fbem_extract_word(line,2)
!        ! Read an integer from word
!        read(word,*) eid
!        ! Check if the written boundary exists. If so, the eid is transformed to iid (sb).
!        if ((eid.ge.boundary_eid_min).and.(eid.le.boundary_eid_max)) then
!          if (boundary_iid(eid).eq.0) then
!            call fbem_error_message(error_unit,0,section_name,nl,'the boundary identifier given in this line does not exist.')
!          else
!            sb=boundary_iid(eid)
!          end if
!        else
!          call fbem_error_message(error_unit,0,section_name,nl,'the boundary identifier given in this line does not exist.')
!        end if
!        ! The second word must be a valid type of boundary condition:
!        !   - u               : Dirichlet u=U
!        !   - j               : Neumann (physical flux) j=k·du/dn=J
!        !   - u(j)-linear     : linear Robin expressed as u(j)
!        !   - u(j)-non-linear : non-linear Robin expressed as u(j)
!        word=fbem_extract_word(line,3)
!        bc=.false.
!        if (trim(word).eq.'u') then
!          bc=.true.
!          read(word,*) tmp_value
!          boundary(i)%ctype(1,1)=0
!          boundary(i)%cvalue_r(1,1,1)=tmp_value
!        end if
!        if (trim(word).eq.'j') then
!          bc=.true.
!          read(word,*) tmp_value
!          boundary(i)%ctype(1,1)=1
!          boundary(i)%cvalue_r(1,1,1)=tmp_value
!        end if
!        if (bc.eqv.(.false.)) call fbem_error_message(error_unit,0,section_name,nl,'invalid boundary condition in this line.')
!      else
!        call fbem_error_message(error_unit,0,section_name,nl,'invalid syntax in this line of the section.')
!      end if
!    end if
!    ! Increment line counter
!    nl=nl+1
!  end do



!
! Old way, ¿to be deprecated?
!

  ! Locate the section
  section_name='conditions over be boundaries'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(input_fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    do i=1,n_boundaries
      ! Find the boundary
      call fbem_search_section(input_fileunit,section_name,found)
      write(keyword,*) 'boundary ', boundary(i)%id
      call fbem_trim2b(keyword)
      call fbem_search_keyword(input_fileunit,keyword,':',found)
      if (found) then
        select case (boundary(i)%coupling)

          ! ===========
          ! BE BOUNDARY
          ! ===========

          case (fbem_boundary_coupling_be)
            select case (boundary(i)%class)

              ! -----------------
              ! ORDINARY BOUNDARY
              ! -----------------

              case (fbem_boundary_class_ordinary)
                ! Read the type of boundary condition
                read(input_fileunit,*) boundary(i)%ctype(1,1)
                call fbem_search_section(input_fileunit,section_name,found)
                call fbem_search_keyword(input_fileunit,keyword,':',found)
                ! Switch depending on the type of boundary condition
                select case (boundary(i)%ctype(1,1))
                  ! 0: p=P
                  ! 1: j=J=k·dp/dn
                  case (0,1)
                    read(input_fileunit,*) boundary(i)%ctype(1,1), boundary(i)%cvalue_r(1,1,1)
                  case default
                    call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                end select

              ! -------------------
              ! CRACK-LIKE BOUNDARY
              ! -------------------

              case (fbem_boundary_class_cracklike)
                ! Face +
                read(input_fileunit,*) boundary(i)%ctype(1,1)
                call fbem_search_section(input_fileunit,section_name,found)
                call fbem_search_keyword(input_fileunit,keyword,':',found)
                ! Switch depending on the type of boundary condition
                select case (boundary(i)%ctype(1,1))
                  ! 0: p^+=P^+
                  ! 1: j^+=J^+=k·dp^+/dn^+
                  case (0,1)
                    read(input_fileunit,*) boundary(i)%ctype(1,1), boundary(i)%cvalue_r(1,1,1)
                  case default
                    call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                end select
                ! Face -
                read(input_fileunit,*) boundary(i)%ctype(1,2)
                backspace(input_fileunit)
                ! Switch depending on the type of boundary condition
                select case (boundary(i)%ctype(1,2))
                  ! 0: p^-=P^-
                  ! 1: j^-=J^-=k·dp^-/dn^-
                  case (0,1)
                    read(input_fileunit,*) boundary(i)%ctype(1,2), boundary(i)%cvalue_r(1,1,2)
                  case default
                    call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'invalid type of boundary condition')
                end select

            end select

          ! ==============
          ! BE-BE BOUNDARY
          ! ==============

          case (fbem_boundary_coupling_be_be)
            ! It does not need B.C.

          ! ==============
          ! BE-FE BOUNDARY
          ! ==============

          case (fbem_boundary_coupling_be_fe)
            stop 'not implemented yet'

          ! =================
          ! BE-FE-BE BOUNDARY
          ! =================

          case (fbem_boundary_coupling_be_fe_be)
            stop 'not implemented yet'

        end select
      end if
    end do

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  end if

end subroutine read_conditions_bem_boundaries_laplace
