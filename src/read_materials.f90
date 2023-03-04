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

!! @author Jacob David Rodriguez Bordon (jacobdavid.rodriguezbordon@ulpgc.es)
!!
!! @version 2.0
!!
!! <b> Subroutine that reads the materials from a file. </b>

subroutine read_materials(fileunit)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_data_structures
  use fbem_harela_incident_field
  use fbem_harpot_incident_field

  ! Problem variables module
  use problem_variables

  implicit none
  ! I/O
  integer                                :: fileunit
  ! Local
  character(len=fbem_stdcharlen)         :: section_name
  logical                                :: found
  character(len=fbem_file_record_length) :: line, word
  integer                                :: n_words
  logical                                :: valid_material, valid_property
  integer                                :: i, j, k, np

  n_materials=0
  section_name='materials'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! Read the number of materials
    read(fileunit,*) n_materials

    ! Allocate and initialize
    if (n_materials.gt.0) then
      allocate (material(n_materials))
    else
      call fbem_error_message(error_unit,0,section_name,0,'the number of materials must be >0')
    end if

    ! Read each material
    do i=1,n_materials
      read(fileunit,'(a)') line
      n_words=fbem_count_words(line)
      if (n_words.lt.2) then
        call fbem_error_message(error_unit,0,section_name,i,'wrong number of arguments in this line')
      end if
      ! Identifier
      word=fbem_extract_word(line,1)
      read(word,*) material(i)%id
      ! Check id
      if (material(i)%id.le.0) then
        call fbem_error_message(error_unit,0,'material',material(i)%id,'identifiers must be greater than 0')
      end if
      ! Type of material
      word=fbem_extract_word(line,2)
      read(word,*) material(i)%type
      call fbem_trimall(material(i)%type)
      valid_material=.false.

      ! =====
      ! FLUID
      ! =====

      !
      ! Pensar en meter distintos tipos de amortiguamiento
      !

      if (trim(material(i)%type).eq.'fluid') then
        !
        ! List of properties:
        !
        !  1 K
        !  2 rho
        !  3 xi
        !  4 c
        !
        valid_material=.true.
        allocate(material(i)%property_defined(4))
        allocate(material(i)%property(4,1))
        material(i)%property_defined=.false.
        material(i)%property=0.d0
        ! Calculate the number of properties introduced by the user
        np=n_words-2
        if (mod(np,2).ne.0) then
          call fbem_error_message(error_unit,0,'material',material(i)%id,'wrong number of arguments')
        end if
        np=np/2
        ! Read each property
        do j=1,np
          valid_property=.false.
          word=fbem_extract_word(line,2+2*j-1)
          if (trim(word).eq.'K') then
            valid_property=.true.
            if (material(i)%property_defined(1)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of K')
            end if
            material(i)%property_defined(1)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(1,1)
          end if
          if (trim(word).eq.'rho') then
            valid_property=.true.
            if (material(i)%property_defined(2)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of rho')
            end if
            material(i)%property_defined(2)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(2,1)
          end if
          if (trim(word).eq.'xi') then
            valid_property=.true.
            if (material(i)%property_defined(3)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of xi')
            end if
            material(i)%property_defined(3)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(3,1)
          end if
          if (trim(word).eq.'c') then
            valid_property=.true.
            if (material(i)%property_defined(4)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of c')
            end if
            material(i)%property_defined(4)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(4,1)
          end if
          if (.not.valid_property) then
            call fbem_error_message(error_unit,0,'material',material(i)%id,'unknown property')
          end if
        end do
        call fbem_harpot_properties(material(i)%property_defined(1),material(i)%property(1,1),&
                                    material(i)%property_defined(2),material(i)%property(2,1),&
                                    material(i)%property_defined(4),material(i)%property(4,1))
        material(i)%property_defined([1,2,4])=.true.
        if (material(i)%property_defined(3)) then
          if (material(i)%property(3,1).lt.0.d0) then
            call fbem_error_message(error_unit,0,'material',material(i)%id,'xi must be >=0')
          end if
        end if
      end if

      ! ===========
      ! RIGID SOLID
      ! ===========

      if (trim(material(i)%type).eq.'rigid_solid') then
        !
        ! List of properties:
        !
        !  1 rho
        !
        valid_material=.true.
        allocate(material(i)%property_defined(1))
        allocate(material(i)%property(1,1))
        material(i)%property_defined=.false.
        material(i)%property=0.d0
        ! Calculate the number of properties introduced by the user
        np=n_words-2
        if (mod(np,2).ne.0) then
          call fbem_error_message(error_unit,0,'material',material(i)%id,'wrong number of arguments')
        end if
        np=np/2
        ! Read each property
        do j=1,np
          valid_property=.false.
          word=fbem_extract_word(line,2+2*j-1)
          if (trim(word).eq.'rho') then
            valid_property=.true.
            if (material(i)%property_defined(1)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of rho')
            end if
            material(i)%property_defined(1)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(1,1)
            if (material(i)%property(1,1).le.0.d0) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'rho must be >0')
            end if
          end if
          if (.not.valid_property) then
            call fbem_error_message(error_unit,0,'material',material(i)%id,'unknown property')
          end if
        end do

        ! TODO: la funcionalidad de la masa del solido rigido no esta todavia metida

      end if

      ! =============
      ! ELASTIC SOLID
      ! =============


      !
      ! TODO: Pensar en meter distintos tipos de amortiguamiento
      !

      if (trim(material(i)%type).eq.'elastic_solid') then
        !
        ! List of properties:
        !
        !  1 E
        !  2 nu
        !  3 lambda
        !  4 mu or G
        !  5 K
        !  6 rho
        !  7 xi
        !
        valid_material=.true.
        allocate(material(i)%property_defined(7))
        allocate(material(i)%property(7,1))
        material(i)%property_defined=.false.
        material(i)%property=0.d0
        ! Calculate the number of properties introduced by the user
        np=n_words-2
        if (mod(np,2).ne.0) then
          call fbem_error_message(error_unit,0,'material',material(i)%id,'wrong number of arguments')
        end if
        np=np/2
        ! Read each property
        do j=1,np
          valid_property=.false.
          word=fbem_extract_word(line,2+2*j-1)
          if (trim(word).eq.'E') then
            valid_property=.true.
            if (material(i)%property_defined(1)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of E')
            end if
            material(i)%property_defined(1)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(1,1)
          end if
          if (trim(word).eq.'nu') then
            valid_property=.true.
            if (material(i)%property_defined(2)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of nu')
            end if
            material(i)%property_defined(2)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(2,1)
          end if
          if (trim(word).eq.'lambda') then
            valid_property=.true.
            if (material(i)%property_defined(3)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of lambda')
            end if
            material(i)%property_defined(3)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(3,1)
          end if
          if ((trim(word).eq.'mu').or.(trim(word).eq.'G')) then
            valid_property=.true.
            if (material(i)%property_defined(4)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of mu')
            end if
            material(i)%property_defined(4)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(4,1)
          end if
          if (trim(word).eq.'K') then
            valid_property=.true.
            if (material(i)%property_defined(5)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of K')
            end if
            material(i)%property_defined(5)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(5,1)
          end if
          if (trim(word).eq.'rho') then
            valid_property=.true.
            if (material(i)%property_defined(6)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of rho')
            end if
            material(i)%property_defined(6)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(6,1)
          end if
          if (trim(word).eq.'xi') then
            valid_property=.true.
            if (material(i)%property_defined(7)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of xi')
            end if
            material(i)%property_defined(7)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(7,1)
          end if
          if (.not.valid_property) then
            call fbem_error_message(error_unit,0,'material',material(i)%id,'unknown property')
          end if
        end do
        ! Calculate the other properties (elastic constants) and check them
        call fbem_ela_properties(material(i)%property_defined(1),material(i)%property(1,1),&
                               material(i)%property_defined(2),material(i)%property(2,1),&
                               material(i)%property_defined(3),material(i)%property(3,1),&
                               material(i)%property_defined(4),material(i)%property(4,1),&
                               material(i)%property_defined(5),material(i)%property(5,1))
        material(i)%property_defined(1:5)=.true.
        if (material(i)%property_defined(6)) then
          if (material(i)%property(6,1).le.0.d0) then
            call fbem_error_message(error_unit,0,'material',material(i)%id,'rho must be >0')
          end if
        end if
        if (material(i)%property_defined(7)) then
          if (material(i)%property(7,1).lt.0.d0) then
            call fbem_error_message(error_unit,0,'material',material(i)%id,'xi must be >=0')
          end if
        end if
      end if

      ! =======================
      ! BIOT POROELASTIC MEDIUM
      ! =======================

      !
      ! Pensar en meter distintos tipos de amortiguamiento
      !

      if (trim(material(i)%type).eq.'biot_poroelastic_medium') then
        !
        ! List of properties:
        !
        !  1 phi
        !  2 E
        !  3 nu
        !  4 lambda
        !  5 mu
        !  6 K
        !  7 Q
        !  8 R
        !  9 rho_f
        ! 10 rho_s
        ! 11 rho_a
        ! 12 xi
        ! 13 b
        !
        valid_material=.true.
        allocate(material(i)%property_defined(13))
        allocate(material(i)%property(13,1))
        material(i)%property_defined=.false.
        material(i)%property=0.d0
        ! Calculate the number of properties introduced by the user
        np=n_words-2
        if (mod(np,2).ne.0) then
          call fbem_error_message(error_unit,0,'material',material(i)%id,'wrong number of arguments')
        end if
        np=np/2
        ! Read each property
        do j=1,np
          valid_property=.false.
          word=fbem_extract_word(line,2+2*j-1)
          if (trim(word).eq.'phi') then
            valid_property=.true.
            if (material(i)%property_defined(1)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of phi')
            end if
            material(i)%property_defined(1)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(1,1)
          end if
          if (trim(word).eq.'E') then
            valid_property=.true.
            if (material(i)%property_defined(2)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of E')
            end if
            material(i)%property_defined(2)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(2,1)
          end if
          if (trim(word).eq.'nu') then
            valid_property=.true.
            if (material(i)%property_defined(3)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of nu')
            end if
            material(i)%property_defined(3)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(3,1)
          end if
          if (trim(word).eq.'lambda') then
            valid_property=.true.
            if (material(i)%property_defined(4)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of lambda')
            end if
            material(i)%property_defined(4)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(4,1)
          end if
          if (trim(word).eq.'mu') then
            valid_property=.true.
            if (material(i)%property_defined(5)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of mu')
            end if
            material(i)%property_defined(5)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(5,1)
          end if
          if (trim(word).eq.'K') then
            valid_property=.true.
            if (material(i)%property_defined(6)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of K')
            end if
            material(i)%property_defined(6)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(6,1)
          end if
          if (trim(word).eq.'Q') then
            valid_property=.true.
            if (material(i)%property_defined(7)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of Q')
            end if
            material(i)%property_defined(7)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(7,1)
          end if
          if (trim(word).eq.'R') then
            valid_property=.true.
            if (material(i)%property_defined(8)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of R')
            end if
            material(i)%property_defined(8)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(8,1)
          end if
          if (trim(word).eq.'rho_f') then
            valid_property=.true.
            if (material(i)%property_defined(9)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of rho_f')
            end if
            material(i)%property_defined(9)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(9,1)
          end if
          if (trim(word).eq.'rho_s') then
            valid_property=.true.
            if (material(i)%property_defined(10)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of rho_s')
            end if
            material(i)%property_defined(10)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(10,1)
          end if
          if (trim(word).eq.'rho_a') then
            valid_property=.true.
            if (material(i)%property_defined(11)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of rho_a')
            end if
            material(i)%property_defined(11)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(11,1)
          end if
          if (trim(word).eq.'xi') then
            valid_property=.true.
            if (material(i)%property_defined(12)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of xi')
            end if
            material(i)%property_defined(12)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(12,1)
          end if
          if (trim(word).eq.'b') then
            valid_property=.true.
            if (material(i)%property_defined(13)) then
              call fbem_error_message(error_unit,0,'material',material(i)%id,'repeated value of b')
            end if
            material(i)%property_defined(13)=.true.
            word=fbem_extract_word(line,2+2*j)
            read(word,*) material(i)%property(13,1)
          end if
          if (.not.valid_property) then
            call fbem_error_message(error_unit,0,'material',material(i)%id,'unknown property')
          end if
        end do
        ! Calculate the other properties (drained elastic constants) and check them
        call fbem_ela_properties(material(i)%property_defined(2),material(i)%property(2,1),&
                               material(i)%property_defined(3),material(i)%property(3,1),&
                               material(i)%property_defined(4),material(i)%property(4,1),&
                               material(i)%property_defined(5),material(i)%property(5,1),&
                               material(i)%property_defined(6),material(i)%property(6,1))
        material(i)%property_defined(2:6)=.true.
        if (material(i)%property_defined(1)) then
          if ((material(i)%property(1,1).le.0.d0).or.(material(i)%property(1,1).ge.1.d0)) then
            call fbem_error_message(error_unit,0,'material',material(i)%id,'phi must be in ]0,1[')
          end if
        end if
        if (material(i)%property_defined(7)) then
          if (material(i)%property(7,1).le.0.d0) then
            call fbem_error_message(error_unit,0,'material',material(i)%id,'Q must be >0')
          end if
        end if
        if (material(i)%property_defined(8)) then
          if (material(i)%property(8,1).le.0.d0) then
            call fbem_error_message(error_unit,0,'material',material(i)%id,'R must be >0')
          end if
        end if
        if (material(i)%property_defined(9)) then
          if (material(i)%property(9,1).le.0.d0) then
            call fbem_error_message(error_unit,0,'material',material(i)%id,'rho_f must be >0')
          end if
        end if
        if (material(i)%property_defined(10)) then
          if (material(i)%property(10,1).le.0.d0) then
            call fbem_error_message(error_unit,0,'material',material(i)%id,'rho_s must be >0')
          end if
        end if
        if (material(i)%property_defined(11)) then
          if (material(i)%property(11,1).le.0.d0) then
            call fbem_error_message(error_unit,0,'material',material(i)%id,'rho_a must be >0')
          end if
        end if
        if (material(i)%property_defined(13)) then
          if (material(i)%property(13,1).le.0.d0) then
            call fbem_error_message(error_unit,0,'material',material(i)%id,'b must be >0')
          end if
        end if
      end if

      ! ========================
      ! UNKNOWN TYPE OF MATERIAL
      ! ========================

      if (.not.valid_material) then
        call fbem_error_message(error_unit,0,'material',material(i)%id,'unknown type of material')
      end if

    end do

    ! ====================================
    ! CHECK AND BUILD MATERIAL IDENTIFIERS
    ! ====================================

    material_eid_min=material(1)%id
    material_eid_max=material(1)%id
    do i=2,n_materials
      if (material(i)%id.lt.material_eid_min) material_eid_min=material(i)%id
      if (material(i)%id.gt.material_eid_max) material_eid_max=material(i)%id
    end do
    allocate (material_iid(material_eid_min:material_eid_max))
    material_iid=0
    do i=1,n_materials
      if (material_iid(material(i)%id).ne.0) then
        call fbem_error_message(error_unit,0,'material',material(i)%id,'is repeated')
      else
        material_iid(material(i)%id)=i
      end if
    end do

    ! Ending message
    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  ! Descomentar las siguientes 2 lineas cuando se termine de implementar los materiales para las regiones y sea obligatorio
  !else

  !  call fbem_error_message(error_unit,0,'['//trim(section_name)//']',0,'this section is required')

  end if

end subroutine read_materials
