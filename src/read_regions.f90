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


!! @author Jacob David Rodriguez Bordon (jacobdavid.rodriguezbordon@ulpgc.es)
!!
!! @version 2.0
!!
!! <b> Subroutine that reads the regions from a file. </b>

subroutine read_regions(fileunit)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_string_handling
  use fbem_numerical
  use fbem_shape_functions
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! Local variables
  implicit none
  ! I/O
  integer                                :: fileunit               ! Unit of the file to read from
  ! Local
  character(len=fbem_stdcharlen)         :: section_name           ! Name of the section
  character(len=fbem_fmtstr)             :: fmtstr                 ! String used for write format string
  logical                                :: found
  character(len=fbem_stdcharlen)         :: tmp_type               ! Temporary variable to read the type of the region
  character(len=fbem_stdcharlen)         :: tmp_variable           ! Temporary variable to read the type of the region
  integer                                :: i, j, si               ! Counters
  complex(kind=real64)                   :: cte_r, cte_q, cte_b, cte_p, cte_trho, cte_phi, cte_rhof
  character(len=fbem_file_record_length) :: line, word
  integer                                :: nw
  integer                                :: km

  section_name='regions'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section [regions]')

    ! Read the number of regions
    read(fileunit,*) n_regions

    ! Allocate and initialize
    if (n_regions.gt.0) then
      allocate (region(n_regions))
      do i=1,n_regions
        region(i)%id=0
        region(i)%name=''
        region(i)%class=0
        region(i)%type=0
        region(i)%subtype=0
        region(i)%model=0
        region(i)%space=0
        region(i)%n_materials=0
        region(i)%n_fe_subregions=0
        region(i)%n_boundaries=0
        region(i)%n_be_bodyloads=0
        region(i)%n_internalpoints=0
        region(i)%n_incidentfields=0

        ! rigid feature
        region(i)%master_node=0

        region(i)%sensitivity=.false.
        region(i)%export=.false.
      end do
    else
      call fbem_error_message(error_unit,0,section_name,0,'the number of regions must be >0')
    end if

    ! Loop through regions
    do i=1,n_regions

      ! ========================
      ! READ REGION ID AND CLASS
      ! ========================

      ! Read
      nw=0
      do while (nw.eq.0)
        read(fileunit,'(a)') line
        nw=fbem_count_words(line)
      end do
      if (nw.lt.2) then
        call fbem_error_message(error_unit,0,section_name,0,'wrong number of arguments in a line of this section')
      end if
      ! Extract and check id
      word=fbem_extract_word(line,1)
      read(word,*) region(i)%id
      if (region(i)%id.le.0) then
        call fbem_error_message(error_unit,0,'region',region(i)%id,'identifiers must be greater than 0')
      end if
      ! Extract and check the region class
      word=fbem_extract_word(line,2)
      region(i)%class=0
      if (trim(word).eq.'be') region(i)%class=fbem_be
      if (trim(word).eq.'fe') region(i)%class=fbem_fe
      if (region(i)%class.eq.0) then
        call fbem_error_message(error_unit,0,'region',region(i)%id,'the region class must be "be" or "fe"')
      end if
      ! Extract and check the region space
      if ((nw.gt.2).and.(region(i)%class.eq.fbem_be)) then
        word=fbem_extract_word(line,3)
        region(i)%space=0
        if (trim(word).eq.'full-space'             ) region(i)%space=fbem_full_space
        if (trim(word).eq.'half-space'             ) region(i)%space=fbem_half_space
        if (trim(word).eq.'multilayered_half-space') region(i)%space=fbem_multilayered_half_space
        if (region(i)%space.eq.0) then
          call fbem_error_message(error_unit,0,'region',region(i)%id,'the region space must be "full-space", "half-space" or "multilayered_half-space"')
        end if
        if (region(i)%space.eq.fbem_half_space) then
          ! Default values
          if (nw.lt.6) then
            select case (problem%n)
              case (2)
                region(i)%halfspace_n=2
              case (3)
                region(i)%halfspace_n=3
            end select
            region(i)%halfspace_x=0.d0
            region(i)%halfspace_bc=1
          else
            ! Hyperplane normal
            word=fbem_extract_word(line,4)
            read(word,*) region(i)%halfspace_n
            if ((abs(region(i)%halfspace_n).lt.1).or.(abs(region(i)%halfspace_n).gt.problem%n)) then
              call fbem_error_message(error_unit,0,'region',region(i)%id,'half-space hyperplane normal specification is invalid')
            end if
            ! Hyperplane coordinate
            word=fbem_extract_word(line,5)
            read(word,*) region(i)%halfspace_x
            ! Hyperplane boundary condition
            word=fbem_extract_word(line,6)
            read(word,*) region(i)%halfspace_bc
          end if
        else if (region(i)%space.eq.fbem_multilayered_half_space) then
          ! Hyperplane normal
          word=fbem_extract_word(line,4)
          read(word,*) region(i)%halfspace_n
          if ((abs(region(i)%halfspace_n).lt.1).or.(abs(region(i)%halfspace_n).gt.problem%n)) then
            call fbem_error_message(error_unit,0,'region',region(i)%id,'half-space hyperplane normal specification is invalid')
          end if
          ! Number of layers/materials
          word=fbem_extract_word(line,5)
          read(word,*) region(i)%n_materials
          allocate (region(i)%material(region(i)%n_materials))
          allocate (region(i)%multilayered_x(region(i)%n_materials+1))
          ! Hyperplane coordinates (layer top) (from free-surface to infinity)
          do j=1,region(i)%n_materials+1
            word=fbem_extract_word(line,5+j)
            read(word,*) region(i)%multilayered_x(j)
          end do
          ! Check layering coordinates
          if (region(i)%halfspace_n.gt.0) then
            do j=2,region(i)%n_materials+1
              if ((region(i)%multilayered_x(j)-region(i)%multilayered_x(j-1)).ge.0.d0) then
                call fbem_error_message(error_unit,0,'region',region(i)%id,'invalid layered half-space coordinates')
              end if
            end do
          else
            do j=2,region(i)%n_materials+1
              if ((region(i)%multilayered_x(j)-region(i)%multilayered_x(j-1)).le.0.d0) then
                call fbem_error_message(error_unit,0,'region',region(i)%id,'invalid layered half-space coordinates')
              end if
            end do
          end if
          ! Hyperplane boundary condition
          word=fbem_extract_word(line,7+region(i)%n_materials)
          read(word,*) region(i)%halfspace_bc
        end if
      ! By default, full-space.
      else
        region(i)%space=fbem_full_space
      end if
      ! Write message
      if ((verbose_level.ge.3).and.(region(i)%class.eq.fbem_be)) then
        select case (region(i)%space)
          case (fbem_full_space)
            write(output_unit,'(2x,a10)') 'Full-space'
          case (fbem_half_space)
            write(output_unit,'(2x,a40)') 'Half-space with the following properties'
            select case (region(i)%halfspace_n)
              case (-1)
                write(output_unit,'(2x,a29)'      ) 'Hyperplane normal: n=(-1,0,0)'
                write(output_unit,'(2x,a26,f15.7)') 'Hyperplane coordinate: x1=',region(i)%halfspace_x
              case ( 1)
                write(output_unit,'(2x,a29)'      ) 'Hyperplane normal: n=(+1,0,0)'
                write(output_unit,'(2x,a26,f15.7)') 'Hyperplane coordinate: x1=',region(i)%halfspace_x
              case (-2)
                write(output_unit,'(2x,a29)'      ) 'Hyperplane normal: n=(0,-1,0)'
                write(output_unit,'(2x,a26,f15.7)') 'Hyperplane coordinate: x2=',region(i)%halfspace_x
              case ( 2)
                write(output_unit,'(2x,a29)'      ) 'Hyperplane normal: n=(0,+1,0)'
                write(output_unit,'(2x,a26,f15.7)') 'Hyperplane coordinate: x2=',region(i)%halfspace_x
              case (-3)
                write(output_unit,'(2x,a29)'      ) 'Hyperplane normal: n=(0,0,-1)'
                write(output_unit,'(2x,a26,f15.7)') 'Hyperplane coordinate: x3=',region(i)%halfspace_x
              case ( 3)
                write(output_unit,'(2x,a29)'      ) 'Hyperplane normal: n=(0,0,+1)'
                write(output_unit,'(2x,a26,f15.7)') 'Hyperplane coordinate: x3=',region(i)%halfspace_x
              case default
                stop 'Fatal error: this should not happen'
            end select
            write(output_unit,'(2x,a34,i2)') 'Hyperplane boundary condition: bc=',region(i)%halfspace_bc
          case (fbem_multilayered_half_space)
            write(output_unit,'(2x,a40)') 'Multilayered half-space with the following properties'
            select case (region(i)%halfspace_n)
              case (-1)
                write(output_unit,'(2x,a29)'      ) 'Hyperplane normal: n=(-1,0,0)'
                write(output_unit,'(2x,a26,f15.7)') 'Hyperplane coordinate: x1=',region(i)%halfspace_x
              case ( 1)
                write(output_unit,'(2x,a29)'      ) 'Hyperplane normal: n=(+1,0,0)'
                write(output_unit,'(2x,a26,f15.7)') 'Hyperplane coordinate: x1=',region(i)%halfspace_x
              case (-2)
                write(output_unit,'(2x,a29)'      ) 'Hyperplane normal: n=(0,-1,0)'
                write(output_unit,'(2x,a26,f15.7)') 'Hyperplane coordinate: x2=',region(i)%halfspace_x
              case ( 2)
                write(output_unit,'(2x,a29)'      ) 'Hyperplane normal: n=(0,+1,0)'
                write(output_unit,'(2x,a26,f15.7)') 'Hyperplane coordinate: x2=',region(i)%halfspace_x
              case (-3)
                write(output_unit,'(2x,a29)'      ) 'Hyperplane normal: n=(0,0,-1)'
                write(output_unit,'(2x,a26,f15.7)') 'Hyperplane coordinate: x3=',region(i)%halfspace_x
              case ( 3)
                write(output_unit,'(2x,a29)'      ) 'Hyperplane normal: n=(0,0,+1)'
                write(output_unit,'(2x,a26,f15.7)') 'Hyperplane coordinate: x3=',region(i)%halfspace_x
              case default
                stop 'Fatal error: this should not happen'
            end select
            write(output_unit,'(2x,a34,i2)') 'Hyperplane boundary condition: bc=',region(i)%halfspace_bc
          case default
            stop 'Fatal error: this should not happen'
        end select
      end if

      ! =====================================================================
      ! READ BE BOUNDARIES (IF A BE REGION) OR FE SUBREGIONS (IF A FE REGION)
      ! =====================================================================

      select case (region(i)%class)

        ! ---------
        ! BE REGION
        ! ---------

        case (fbem_be)
          read(fileunit,*) region(i)%n_boundaries
          if (region(i)%n_boundaries.gt.0) then
            allocate (region(i)%boundary(region(i)%n_boundaries))
            allocate (region(i)%boundary_reversion(region(i)%n_boundaries))
            backspace(fileunit)
            read(fileunit,*) region(i)%n_boundaries, (region(i)%boundary(j),j=1,region(i)%n_boundaries)
            do j=1,region(i)%n_boundaries
              if (region(i)%boundary(j).lt.0) then
                region(i)%boundary_reversion(j)=.true.
                region(i)%boundary(j)=-region(i)%boundary(j)
              else
                region(i)%boundary_reversion(j)=.false.
              end if
            end do
          end if

        ! ---------
        ! FE REGION
        ! ---------

        case (fbem_fe)
          read(fileunit,*) region(i)%n_fe_subregions
          allocate (region(i)%fe_subregion(region(i)%n_fe_subregions))
          backspace(fileunit)
          read(fileunit,*) region(i)%n_fe_subregions, (region(i)%fe_subregion(j),j=1,region(i)%n_fe_subregions)

      end select

      ! =======================
      ! READ MATERIAL AND MODEL
      ! =======================

      ! Continue reading depending on the problem type
      select case (problem%type)

        ! ==========================================================================================================================
        ! MECHANICS
        ! ==========================================================================================================================

        case (fbem_mechanics)
          ! Switch between types of analysis
          select case (problem%analysis)

            ! ----------------------------------------------------------------------------------------------------------------------
            ! STATIC
            ! ----------------------------------------------------------------------------------------------------------------------

            case (fbem_static)

              if (region(i)%space.eq.fbem_multilayered_half_space) stop 'not multilayered halfspace for elastostatics yet'

              read(fileunit,'(a)') line
              nw=fbem_count_words(line)
              word=fbem_extract_word(line,1)
              if (trim(word).eq.'material') then

                ! Read the material id
                word=fbem_extract_word(line,2)
                read(word,*) km
                km=material_iid(km)

                ! -----------
                ! RIGID SOLID
                ! -----------

                if (trim(material(km)%type).eq.'rigid_solid') then

                  if (region(i)%class.ne.fbem_fe) stop 'Rigid solid materials only for FE regions'
                  region(i)%type=fbem_rigid
                  region(i)%subtype=0
                  if (trim(fbem_extract_word(line,3)).eq.'master_node_x') then
                    allocate(region(i)%master_node_x(problem%n))
                    do j=1,problem%n
                      word=fbem_extract_word(line,3+j)
                      read(word,*) region(i)%master_node_x(j)
                    end do
                    write(*,*) region(i)%master_node_x
                  else if (trim(fbem_extract_word(line,3)).eq.'master_node') then
                    word=fbem_extract_word(line,4)
                    read(word,*) region(i)%master_node
                    write(*,*) region(i)%master_node
                  else
                    call fbem_error_message(error_unit,0,'region',region(i)%id,'requires master_node_x or master_node')
                  end if

                ! ------------- !
                ! ELASTIC SOLID !
                ! ------------- !

                else if (trim(material(km)%type).eq.'elastic_solid') then
                  !
                  ! Complete list of region properties:
                  !
                  !   - property_r(1): -
                  !   - property_r(2): mu, shear modulus
                  !   - property_r(3): nu, Poisson's ratio
                  !   - property_r(4): -
                  !   - property_r(5): E, Young modulus
                  !   - property_r(6): lambda, Lame constant
                  !
                  region(i)%type=fbem_viscoelastic
                  region(i)%subtype=0
                  ! Allocate region properties
                  allocate(region(i)%property_r(6))
                  !
                  ! Read from "material" data structure
                  !
                  if (material(km)%property_defined(4)) then
                    region(i)%property_r(2)=material(km)%property(4,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'mu is required for the material of this region')
                  end if
                  if (material(km)%property_defined(2)) then
                    region(i)%property_r(3)=material(km)%property(2,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'nu is required for the material of this region')
                  end if
                  if (material(km)%property_defined(1)) then
                    region(i)%property_r(5)=material(km)%property(1,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'E is required for the material of this region')
                  end if
                  if (material(km)%property_defined(3)) then
                    region(i)%property_r(6)=material(km)%property(3,1)
                   else
                    call fbem_error_message(error_unit,0,section_name,i,'lambda is required for the material of this region')
                  end if

                ! -------------------------------------------------- !
                ! BIOT POROELASTIC MEDIUM (UNDER DRAINED CONDITIONS) !
                ! -------------------------------------------------- !

                else if (trim(material(km)%type).eq.'biot_poroelastic_medium') then
                  !
                  ! Complete list of region properties:
                  !
                  !   - property_r(1): -
                  !   - property_r(2): mu, shear modulus
                  !   - property_r(3): nu, Poisson's ratio
                  !   - property_r(4): -
                  !   - property_r(5): E, Young modulus
                  !   - property_r(6): lambda, Lame constant
                  !
                  region(i)%type=fbem_viscoelastic
                  region(i)%subtype=0
                  ! Allocate region properties
                  allocate(region(i)%property_r(6))
                  !
                  ! Read from "material" data structure
                  !
                  if (material(km)%property_defined(5)) then
                    region(i)%property_r(2)=material(km)%property(5,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'mu is required for the material of this region')
                  end if
                  if (material(km)%property_defined(3)) then
                    region(i)%property_r(3)=material(km)%property(3,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'nu is required for the material of this region')
                  end if
                  if (material(km)%property_defined(2)) then
                    region(i)%property_r(5)=material(km)%property(2,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'E is required for the material of this region')
                  end if
                  if (material(km)%property_defined(4)) then
                    region(i)%property_r(6)=material(km)%property(4,1)
                   else
                    call fbem_error_message(error_unit,0,section_name,i,'lambda is required for the material of this region')
                  end if

                else
                  call fbem_error_message(error_unit,0,section_name,i,'invalid or unknown type of material')
                end if

              else if (trim(word).eq.'orthotropic') then
                !
                ! ORTHOTROPIC (3D)
                !
                if (problem%n.ne.3) stop 'Sorry, orthotropic only for 3D'
                if (region(i)%class.ne.fbem_fe) stop 'Only for FEM shells'
                region(i)%type=fbem_viscoelastic
                region(i)%subtype=fbem_elastic_orthotropic
                ! Check if correct number of parameters
                if (nw.ne.10) then
                  call fbem_error_message(error_unit,0,section_name,nw,'the number of words in this line must be 10')
                end if
                ! Complete list of properties
                !   - property_r( 1): E11
                !   - property_r( 2): E22
                !   - property_r( 3): E33
                !   - property_r( 4): G12 (=G21)
                !   - property_r( 5): G13 (=G31)
                !   - property_r( 6): G23 (=G32)
                !   - property_r( 7): nu12
                !   - property_r( 8): nu21
                !   - property_r( 9): nu13
                !   - property_r(10): nu31
                !   - property_r(11): nu23
                !   - property_r(12): nu32
                ! Read region properties:
                ! 9 constants are needed: E11, E22, E33, G12, G13, G23, nu12, nu13, nu23
                allocate(region(i)%property_r(12))
                word=fbem_extract_word(line, 2)
                read(word,*) region(i)%property_r( 1)
                word=fbem_extract_word(line, 3)
                read(word,*) region(i)%property_r( 2)
                word=fbem_extract_word(line, 4)
                read(word,*) region(i)%property_r( 3)
                word=fbem_extract_word(line, 5)
                read(word,*) region(i)%property_r( 4)
                word=fbem_extract_word(line, 6)
                read(word,*) region(i)%property_r( 5)
                word=fbem_extract_word(line, 7)
                read(word,*) region(i)%property_r( 6)
                word=fbem_extract_word(line, 8)
                read(word,*) region(i)%property_r( 7)
                word=fbem_extract_word(line, 9)
                read(word,*) region(i)%property_r( 9)
                word=fbem_extract_word(line,10)
                read(word,*) region(i)%property_r(11)
                ! Reciprocity relations:
                ! E11*nu12=E22*nu21 => nu21=E11/E22*nu12
                ! E11*nu13=E33*nu31 => nu31=E11/E33*nu13
                ! E22*nu23=E33*nu32 => nu32=E22/E33*nu23
                region(i)%property_r( 8)=region(i)%property_r(1)/region(i)%property_r(2)*region(i)%property_r( 7)
                region(i)%property_r(10)=region(i)%property_r(1)/region(i)%property_r(3)*region(i)%property_r( 9)
                region(i)%property_r(12)=region(i)%property_r(2)/region(i)%property_r(3)*region(i)%property_r(11)
              else
                !
                ! ISOTROPIC
                !

                !
                ! En teoria para BEM tension plana, hay que poner nu'=nu/(1+nu)
                !

                region(i)%type=fbem_viscoelastic
                region(i)%subtype=0
                ! Allocate region properties
                allocate(region(i)%property_r(6))
                ! Read region properties (real 3D values):
                !   - property_r(2): mu, shear modulus
                !   - property_r(3): nu, Poisson's ratio
                read(line,*) region(i)%property_r(2), region(i)%property_r(3)
                ! Complete list of properties
                !   - property_r(1): -
                !   - property_r(2): mu, shear modulus
                !   - property_r(3): nu, Poisson's ratio
                !   - property_r(4): -
                !   - property_r(5): E, Young modulus
                !   - property_r(6): lambda, Lame constant
                ! Convert
                region(i)%property_r(5)=2.0d0*region(i)%property_r(2)*(1.0d0+region(i)%property_r(3))
                region(i)%property_r(6)=2.0d0*region(i)%property_r(2)*region(i)%property_r(3)/(1.0d0-2.0d0*region(i)%property_r(3))
              end if

            ! ----------------------------------------------------------------------------------------------------------------------
            ! HARMONIC
            ! ----------------------------------------------------------------------------------------------------------------------

            case (fbem_harmonic)

              ! Read the type of region
              read(fileunit,*) tmp_type
              ! Back to the previous line
              backspace(fileunit)
              ! Assign the default invalid value to the type
              region(i)%type=0
              if (trim(tmp_type).eq.'material') then
                region(i)%type   =0
                region(i)%subtype=0
              else if ((trim(tmp_type).eq.'fluid').or.(trim(tmp_type).eq.'inviscid_fluid')) then
                region(i)%type   =fbem_potential
                region(i)%subtype=0
              else if (trim(tmp_type).eq.'hysteretic_fluid') then
                region(i)%type   =fbem_potential
                region(i)%subtype=1
              else if (trim(tmp_type).eq.'viscoelastic') then
                region(i)%type   =fbem_viscoelastic
                region(i)%subtype=0
              else if (trim(tmp_type).eq.'viscoelastic2') then
                region(i)%type   =fbem_viscoelastic
                region(i)%subtype=0
              else if (trim(tmp_type).eq.'elastic') then
                region(i)%type   =fbem_viscoelastic
                region(i)%subtype=fbem_elastic
              else if (trim(tmp_type).eq.'viscoelastic_orthotropic') then
                region(i)%type   =fbem_viscoelastic
                region(i)%subtype=fbem_elastic_orthotropic
              else if (trim(tmp_type).eq.'bardet') then
                region(i)%type   =fbem_viscoelastic
                region(i)%subtype=fbem_bardet
              else if (trim(tmp_type).eq.'brt_cp1') then
                region(i)%type   =fbem_viscoelastic
                region(i)%subtype=fbem_brt_cp1
              else if (trim(tmp_type).eq.'brt_cp2') then
                region(i)%type   =fbem_viscoelastic
                region(i)%subtype=fbem_brt_cp2
              else if (trim(tmp_type).eq.'brt_cpm') then
                region(i)%type   =fbem_viscoelastic
                region(i)%subtype=fbem_brt_cpm
              else if (trim(tmp_type).eq.'poroelastic') then
                region(i)%type   =fbem_poroelastic
                region(i)%subtype=0
              else if (trim(tmp_type).eq.'poroelastic2') then
                region(i)%type   =fbem_poroelastic
                region(i)%subtype=0
              else
                call fbem_error_message(error_unit,0,'region',region(i)%id,'the region type is invalid')
              end if

              !
              ! Read from "material" data structure
              !
              if (region(i)%type.eq.0) then

                ! Read the material id
                read(fileunit,'(a)') line
                word=fbem_extract_word(line,2)
                read(word,*) km
                km=material_iid(km)

                ! ----- !
                ! FLUID !
                ! ----- !

                if (trim(material(km)%type).eq.'fluid') then
                  region(i)%type   =fbem_potential
                  region(i)%subtype=0
                  if (region(i)%space.eq.fbem_multilayered_half_space) stop 'not multilayered halfspace for fluids yet'
                  !
                  ! region()%property_*
                  ! -------------------
                  !
                  !   - property_r(1): rho, density (only real)
                  !   - property_?(2): K, bulk modulus (real or complex)
                  !   - property_r(3): xi, hysteretic damping ratio (only real)
                  !   - property_?(4): c, wave propagation speed (real or complex)
                  !
                  ! Allocate region properties
                  allocate(region(i)%property_r(4))
                  allocate(region(i)%property_c(4))
                  region(i)%property_r=0.
                  region(i)%property_c=0.
                  !
                  ! Read from "material" data structure
                  !
                  !  1 K
                  !  2 rho
                  !  3 xi
                  !  4 c
                  !
                  if (material(km)%property_defined(2)) then
                    region(i)%property_r(1)=material(km)%property(2,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'rho is required for the material of this region')
                  end if
                  if (material(km)%property_defined(1)) then
                    region(i)%property_r(2)=material(km)%property(1,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'K is required for the material of this region')
                  end if
                  if (material(km)%property_defined(3)) then
                    region(i)%property_r(3)=material(km)%property(3,1)
                  else
                    region(i)%property_r(3)=0.d0
                  end if
                  if (material(km)%property_defined(4)) then
                    region(i)%property_r(4)=material(km)%property(4,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'c is required for the material of this region')
                  end if
                  region(i)%property_c(2)=region(i)%property_r(2)*(1+c_im*2*region(i)%property_r(3))
                  region(i)%property_c(4)=sqrt(region(i)%property_c(2)/region(i)%property_r(1))

                ! ------------- !
                ! ELASTIC SOLID !
                ! ------------- !

                else if (trim(material(km)%type).eq.'elastic_solid') then
                  region(i)%type   =fbem_viscoelastic
                  region(i)%subtype=0
                  if (region(i)%space.eq.fbem_half_space) stop 'not halfspace for elastodynamics yet'
                  if (region(i)%space.eq.fbem_multilayered_half_space) stop 'not multilayered halfspace for elastodynamics yet'
                  !
                  ! Complete list of region properties:
                  !
                  !   - property_r(1): rho, density
                  !   - property_r(2): mu, shear modulus
                  !   - property_r(3): nu, Poisson's ratio
                  !   - property_r(4): xi, damping coefficient
                  !   - property_r(5): E, Young modulus
                  !   - property_r(6): lambda, Lame constant
                  !   - property_r(7): c1, primary wave (P-wave) propagation speed (c_p)
                  !   - property_r(8): c2, secondary wave (S-wave) propagation speed (c_s)
                  !   - property_c(2): mu, complex shear modulus, includes damping
                  !   - property_c(5): E, complex Young modulus, includes damping
                  !   - property_c(6): lambda, complex Lame constant, includes damping
                  !   - property_c(7): c1, complex primary wave (P-wave) propagation speed (c_p), includes damping
                  !   - property_c(8): c2, complex secondary wave (S-wave) propagation speed (c_s), includes damping
                  !
                  ! Allocate region properties
                  allocate(region(i)%property_r(8))
                  allocate(region(i)%property_c(2:8))
                  region(i)%property_r=0.
                  region(i)%property_c=0.
                  !
                  ! Read from "material" data structure
                  !
                  if (material(km)%property_defined(6)) then
                    region(i)%property_r(1)=material(km)%property(6,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'rho is required for the material of this region')
                  end if
                  if (material(km)%property_defined(4)) then
                    region(i)%property_r(2)=material(km)%property(4,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'mu is required for the material of this region')
                  end if
                  if (material(km)%property_defined(2)) then
                    region(i)%property_r(3)=material(km)%property(2,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'nu is required for the material of this region')
                  end if
                  if (material(km)%property_defined(7)) then
                    region(i)%property_r(4)=material(km)%property(7,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'xi is required for the material of this region')
                  end if
                  if (material(km)%property_defined(1)) then
                    region(i)%property_r(5)=material(km)%property(1,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'E is required for the material of this region')
                  end if
                  if (material(km)%property_defined(3)) then
                    region(i)%property_r(6)=material(km)%property(3,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'lambda is required for the material of this region')
                  end if
                  region(i)%property_r(7)=dsqrt((region(i)%property_r(6)+2.0d0*region(i)%property_r(2))/region(i)%property_r(1))
                  region(i)%property_r(8)=dsqrt(region(i)%property_r(2)/region(i)%property_r(1))
                  region(i)%property_c(2)=region(i)%property_r(2)*(1.0d0+c_im*2.0d0*region(i)%property_r(4))
                  region(i)%property_c(5)=2.0d0*region(i)%property_c(2)*(1.0d0+region(i)%property_r(3))
                  region(i)%property_c(6)=2.0d0*region(i)%property_c(2)*region(i)%property_r(3)/(1.0d0-2.0d0*region(i)%property_r(3))
                  region(i)%property_c(7)=zsqrt((region(i)%property_c(6)+2.0d0*region(i)%property_c(2))/region(i)%property_r(1))
                  region(i)%property_c(8)=zsqrt(region(i)%property_c(2)/region(i)%property_r(1))

                ! ----------------------- !
                ! BIOT POROELASTIC MEDIUM !
                ! ----------------------- !

                else if (trim(material(km)%type).eq.'biot_poroelastic_medium') then
                  region(i)%type   =fbem_poroelastic
                  region(i)%subtype=0
                  if (region(i)%space.eq.fbem_half_space) stop 'not halfspace for poroelasticity yet'
                  if (region(i)%space.eq.fbem_multilayered_half_space) stop 'not multilayered halfspace for poroelasticity yet'
                  !
                  ! List of properties
                  ! ------------------
                  !
                  !   - property_r(1) : rhof, density of the fluid
                  !   - property_r(2) : rhos, density of the solid skeleton
                  !   - property_r(3) : lambda, Lame constant
                  !   - property_r(4) : mu, shear modulus
                  !   - property_r(5) : E, Young modulus
                  !   - property_r(6) : nu, Poisson's ratio
                  !   - property_r(7) : xi, damping coefficient
                  !   - property_r(8) : phi, porosity
                  !   - property_r(9) : rhoa, additional density
                  !   - property_r(10): R coupling parameter
                  !   - property_r(11): Q coupling parameter
                  !   - property_r(12): b coupling parameter
                  !   - property_r(13): rho1, density of the solid skeleton with respect to the total volume
                  !   - property_r(14): rho2, density of the fluid with respect to the total volume
                  !   - property_c(3) : lambda, Lame constant, includes damping
                  !   - property_c(4) : mu, shear modulus, includes damping
                  !   - property_c(5) : E, Young modulus, includes damping
                  !   - property_c(6) : nu, Poisson's ratio, includes damping
                  !   - property_c(10): R coupling parameter, includes damping
                  !   - property_c(11): Q coupling parameter, includes damping
                  !
                  ! Allocate region properties
                  allocate(region(i)%property_r(14))
                  allocate(region(i)%property_c(3:11))
                  region(i)%property_r=0.d0
                  region(i)%property_c=0.d0
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
                  if (material(km)%property_defined(1)) then
                    region(i)%property_r(8)=material(km)%property(1,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'phi is required for the material of this region')
                  end if
                  if (material(km)%property_defined(2)) then
                    region(i)%property_r(5)=material(km)%property(2,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'E is required for the material of this region')
                  end if
                  if (material(km)%property_defined(3)) then
                    region(i)%property_r(6)=material(km)%property(3,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'nu is required for the material of this region')
                  end if
                  if (material(km)%property_defined(4)) then
                    region(i)%property_r(3)=material(km)%property(4,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'lambda is required for the material of this region')
                  end if
                  if (material(km)%property_defined(5)) then
                    region(i)%property_r(4)=material(km)%property(5,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'mu is required for the material of this region')
                  end if
                  if (material(km)%property_defined(7)) then
                    region(i)%property_r(11)=material(km)%property(7,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'Q is required for the material of this region')
                  end if
                  if (material(km)%property_defined(8)) then
                    region(i)%property_r(10)=material(km)%property(8,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'R is required for the material of this region')
                  end if
                  if (material(km)%property_defined(9)) then
                    region(i)%property_r(1)=material(km)%property(9,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'rhof is required for the material of this region')
                  end if
                  if (material(km)%property_defined(10)) then
                    region(i)%property_r(2)=material(km)%property(10,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'rhos is required for the material of this region')
                  end if
                  if (material(km)%property_defined(11)) then
                    region(i)%property_r(9)=material(km)%property(11,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'rhoa is required for the material of this region')
                  end if
                  if (material(km)%property_defined(12)) then
                    region(i)%property_r(7)=material(km)%property(12,1)
                  else
                    region(i)%property_r(7)=0
                  end if
                  if (material(km)%property_defined(13)) then
                    region(i)%property_r(12)=material(km)%property(13,1)
                  else
                    call fbem_error_message(error_unit,0,section_name,i,'b is required for the material of this region')
                  end if

                  region(i)%property_c( 3)=region(i)%property_r( 3)*(1+c_im*2*region(i)%property_r(7))
                  region(i)%property_c( 4)=region(i)%property_r( 4)*(1+c_im*2*region(i)%property_r(7))
                  region(i)%property_c( 5)=region(i)%property_r( 5)*(1+c_im*2*region(i)%property_r(7))
                  region(i)%property_c( 6)=0.5d0*region(i)%property_c(3)/(region(i)%property_c(3)+region(i)%property_c(4))
                  ! Duda: Q y R se ven afectadas por el amort. histeretico. Por lo que pone el libro de
                  ! Detournay y Cheng sí.
                  region(i)%property_c(10)=region(i)%property_r(10)*(1+c_im*2*region(i)%property_r(7))
                  region(i)%property_c(11)=region(i)%property_r(11)*(1+c_im*2*region(i)%property_r(7))

                else

                  call fbem_error_message(error_unit,0,section_name,i,'invalid or unknown type of material')

                end if

              !
              ! Read properties directly
              !
              else

                !
                ! Read region properties
                !
                select case (region(i)%type)

                  ! -----
                  ! FLUID
                  ! -----

                  ! List of properties
                  ! ------------------
                  !
                  !   - property_r(1): rho, density (only real)
                  !   - property_?(2): K, bulk modulus (real or complex)
                  !   - property_r(3): xi, hysteretic damping ratio (only real)
                  !   - property_?(4): c, wave propagation speed (real or complex)
                  !
                  case (fbem_potential)
                    if (region(i)%space.eq.fbem_multilayered_half_space) stop 'not multilayered halfspace for fluid yet'
                    ! Allocate region properties
                    allocate(region(i)%property_r(4))
                    allocate(region(i)%property_c(4))
                    region(i)%property_r=0.
                    region(i)%property_c=0.
                    ! Read region properties:
                    select case (region(i)%subtype)
                      ! Read rho (real), c (real)
                      case (0)
                        read(fileunit,*) tmp_type, region(i)%property_r(1), region(i)%property_r(4)
                        region(i)%property_r(2)=region(i)%property_r(1)*(region(i)%property_r(4))**2
                        region(i)%property_c(2)=region(i)%property_r(2)
                        region(i)%property_c(4)=region(i)%property_r(4)
                      ! Read rho (real), K (real), xi (real)
                      case (1)
                        read(fileunit,*) tmp_type, region(i)%property_r(1), region(i)%property_r(2), region(i)%property_r(3)
                        region(i)%property_c(2)=region(i)%property_r(2)*(1.d0+c_im*2.d0*region(i)%property_r(3))
                        region(i)%property_r(4)=sqrt(region(i)%property_r(2)/region(i)%property_r(1))
                        region(i)%property_c(4)=sqrt(region(i)%property_c(2)/region(i)%property_r(1))
                    end select

                  ! ------------------
                  ! VISCOELASTIC SOLID
                  ! ------------------

                  case (fbem_viscoelastic)

                    if (region(i)%space.eq.fbem_half_space) stop 'not halfspace for elasticity yet'

                    select case (region(i)%subtype)
                      !
                      ! VISCOELASTIC
                      !
                      ! List of properties
                      ! ------------------
                      !
                      !   - property_r(1): rho, density
                      !   - property_r(2): mu, shear modulus
                      !   - property_r(3): nu, Poisson's ratio
                      !   - property_r(4): xi, damping coefficient
                      !   - property_r(5): E, Young modulus
                      !   - property_r(6): lambda, Lame constant
                      !   - property_r(7): c1, primary wave (P-wave) propagation speed (c_p)
                      !   - property_r(8): c2, secondary wave (S-wave) propagation speed (c_s)
                      !   - property_c(2): mu, complex shear modulus, includes damping
                      !   - property_c(3): nu, complex (imaginary part different from 0 only if damping different in mu and lambda)
                      !   - property_c(5): E, complex Young modulus, includes damping
                      !   - property_c(6): lambda, complex Lame constant, includes damping
                      !   - property_c(7): c1, complex primary wave (P-wave) propagation speed (c_p), includes damping
                      !   - property_c(8): c2, complex secondary wave (S-wave) propagation speed (c_s), includes damping
                      !
                      case (0)
                        ! Allocate region properties
                        allocate(region(i)%property_r(8))
                        allocate(region(i)%property_c(2:8))

                        if (trim(tmp_type).eq.'viscoelastic') then
                          ! Read region properties:
                          !   - property_r(1): rho, density
                          !   - property_r(2): mu, shear modulus
                          !   - property_r(3): nu, Poisson's ratio
                          !   - property_r(4): xi, damping coefficient
                          read(fileunit,*) tmp_type, region(i)%property_r(1), region(i)%property_r(2), region(i)%property_r(3),&
                                                     region(i)%property_r(4)
                          ! Conversion
                          region(i)%property_r(5)=2.0d0*region(i)%property_r(2)*(1.0d0+region(i)%property_r(3))
                          region(i)%property_r(6)=2.0d0*region(i)%property_r(2)*region(i)%property_r(3)/(1.0d0-2.0d0*region(i)%property_r(3))
                          region(i)%property_r(7)=dsqrt((region(i)%property_r(6)+2.0d0*region(i)%property_r(2))/region(i)%property_r(1))
                          region(i)%property_r(8)=dsqrt(region(i)%property_r(2)/region(i)%property_r(1))
                          region(i)%property_c(2)=region(i)%property_r(2)*(1.0d0+c_im*2.0d0*region(i)%property_r(4))
                          region(i)%property_c(3)=region(i)%property_r(3)
                          region(i)%property_c(5)=2.0d0*region(i)%property_c(2)*(1.0d0+region(i)%property_r(3))
                          region(i)%property_c(6)=2.0d0*region(i)%property_c(2)*region(i)%property_r(3)/(1.0d0-2.0d0*region(i)%property_r(3))
                          region(i)%property_c(7)=zsqrt((region(i)%property_c(6)+2.0d0*region(i)%property_c(2))/region(i)%property_r(1))
                          region(i)%property_c(8)=zsqrt(region(i)%property_c(2)/region(i)%property_r(1))

                        else

                          ! Read region properties:
                          !   - property_r(1): rho, density
                          !   - property_c(6): K, bulk modulus (complex)
                          !   - property_c(2): mu, shear modulus (complex)
                          read(fileunit,*) tmp_type, region(i)%property_r(1), region(i)%property_c(6), region(i)%property_c(2)

                          ! Conversion
                          region(i)%property_r(2)=dble(region(i)%property_c(2))
                          region(i)%property_c(3)=0.5d0*(3.d0*region(i)%property_c(6)-2.d0*region(i)%property_c(2))/(3.d0*region(i)%property_c(6)+region(i)%property_c(2))
                          region(i)%property_r(3)=dble(region(i)%property_c(3))

                          ! Amortiguamiento medio entre el de mu y K
                          region(i)%property_r(4)=0.5d0*(imag(region(i)%property_c(2))/(2.d0*dble(region(i)%property_c(2)))+imag(region(i)%property_c(6))/(2.d0*dble(region(i)%property_c(6))))

                          region(i)%property_r(5)=dble(9.d0*region(i)%property_c(6)*region(i)%property_c(2)/(3.d0*region(i)%property_c(6)+region(i)%property_c(2)))
                          region(i)%property_r(6)=dble(region(i)%property_c(6)-2.d0/3.d0*region(i)%property_c(2))
                          region(i)%property_r(7)=dsqrt((region(i)%property_r(6)+2.0d0*region(i)%property_r(2))/region(i)%property_r(1))
                          region(i)%property_r(8)=dsqrt(region(i)%property_r(2)/region(i)%property_r(1))
                          region(i)%property_c(5)=9.d0*region(i)%property_c(6)*region(i)%property_c(2)/(3.d0*region(i)%property_c(6)+region(i)%property_c(2))
                          region(i)%property_c(6)=region(i)%property_c(6)-2.d0/3.d0*region(i)%property_c(2)
                          region(i)%property_c(7)=zsqrt((region(i)%property_c(6)+2.0d0*region(i)%property_c(2))/region(i)%property_r(1))
                          region(i)%property_c(8)=zsqrt(region(i)%property_c(2)/region(i)%property_r(1))



                        end if

                      !
                      ! ELASTIC
                      !
                      ! List of properties
                      ! ------------------
                      !
                      !   - property_r(1): rho, density
                      !   - property_r(2): mu, shear modulus
                      !   - property_r(3): nu, Poisson's ratio
                      !   - property_r(4): =0.0d0
                      !   - property_r(5): E, Young modulus
                      !   - property_r(6): lambda, Lame constant
                      !   - property_r(7): c1, primary wave (P-wave) propagation speed (c_p)
                      !   - property_r(8): c2, secondary wave (S-wave) propagation speed (c_s)
                      !   - property_c(2): =property_r(2)
                      !   - property_c(5): =property_r(5)
                      !   - property_c(6): =property_r(6)
                      !   - property_c(7): =property_r(7)
                      !   - property_c(8): =property_r(8)
                      case (fbem_elastic)
                        ! Allocate region properties
                        allocate(region(i)%property_r(8))
                        allocate(region(i)%property_c(2:8))
                        ! Read region properties:
                        !   - property_r(1): rho, density
                        !   - property_r(2): mu, shear modulus
                        !   - property_r(3): nu, Poisson's ratio
                        read(fileunit,*) tmp_type, region(i)%property_r(1), region(i)%property_r(2), region(i)%property_r(3)
                        ! Conversion
                        region(i)%property_r(4)=0.0d0
                        region(i)%property_r(5)=2.0d0*region(i)%property_r(2)*(1.0d0+region(i)%property_r(3))
                        region(i)%property_r(6)=2.0d0*region(i)%property_r(2)*region(i)%property_r(3)/(1.0d0-2.0d0*region(i)%property_r(3))
                        region(i)%property_r(7)=dsqrt((region(i)%property_r(6)+2.0d0*region(i)%property_r(2))/region(i)%property_r(1))
                        region(i)%property_r(8)=dsqrt(region(i)%property_r(2)/region(i)%property_r(1))
                        region(i)%property_c(2)=region(i)%property_r(2)
                        region(i)%property_c(3)=region(i)%property_r(3)
                        region(i)%property_c(5)=region(i)%property_r(5)
                        region(i)%property_c(6)=region(i)%property_r(6)
                        region(i)%property_c(7)=region(i)%property_r(7)
                        region(i)%property_c(8)=region(i)%property_r(8)
                      !
                      ! VISCOELASTIC ORTHOTROPIC
                      !
                      !   - property_c( 1): E11
                      !   - property_c( 2): E22
                      !   - property_c( 3): E33
                      !   - property_c( 4): G12 (=G21)
                      !   - property_c( 5): G13 (=G31)
                      !   - property_c( 6): G23 (=G32)
                      !   - property_r( 7): nu12
                      !   - property_r( 8): nu21
                      !   - property_r( 9): nu13
                      !   - property_r(10): nu31
                      !   - property_r(11): nu23
                      !   - property_r(12): nu32
                      !   - property_r(13): rho
                      case (fbem_elastic_orthotropic)
                        if (problem%n.ne.3) stop 'Sorry, orthotropic only for 3D'
                        if (region(i)%class.ne.fbem_fe) stop 'Only for FEM shells'
                        allocate(region(i)%property_r(7:13))
                        allocate(region(i)%property_c(1:6))
                        ! Read region properties:
                        ! E11, E22, E33, G12, G13, G23, nu12, nu13, nu23, rho
                        read(fileunit,*) tmp_type, region(i)%property_c( 1), region(i)%property_c( 2), region(i)%property_c( 3), &
                                                   region(i)%property_c( 4), region(i)%property_c( 5), region(i)%property_c( 6), &
                                                   region(i)%property_r( 7), region(i)%property_r( 9), region(i)%property_r(11), &
                                                   region(i)%property_r(13)
                        ! Reciprocity relations:
                        ! E11*nu12=E22*nu21 => nu21=E11/E22*nu12
                        ! E11*nu13=E33*nu31 => nu31=E11/E33*nu13
                        ! E22*nu23=E33*nu32 => nu32=E22/E33*nu23
                        region(i)%property_r( 8)=dreal(region(i)%property_c(1)/region(i)%property_c(2)*region(i)%property_r( 7))
                        region(i)%property_r(10)=dreal(region(i)%property_c(1)/region(i)%property_c(3)*region(i)%property_r( 9))
                        region(i)%property_r(12)=dreal(region(i)%property_c(2)/region(i)%property_c(3)*region(i)%property_r(11))
                      !
                      ! Bardet's viscoelasticity model of poroelasticity
                      !
                      ! List of properties
                      ! ------------------
                      !
                      !   - property_r(1) : rhof, density of the fluid
                      !   - property_r(2) : rhos, density of the solid skeleton
                      !   - property_r(3) : lambda, Lame constant
                      !   - property_r(4) : mu, shear modulus
                      !   - property_r(5) : E, Young modulus
                      !   - property_r(6) : nu, Poisson's ratio
                      !   - property_r(7) : xi, damping coefficient
                      !   - property_r(8) : phi, porosity
                      !   - property_r(9) : rhoa, additional density
                      !   - property_r(10): R coupling parameter
                      !   - property_r(11): Q coupling parameter
                      !   - property_r(12): b coupling parameter
                      !   - property_r(13): rho1, density of the solid skeleton with respect to the total volume
                      !   - property_r(14): rho2, density of the fluid with respect to the total volume
                      !   - property_c(3) : lambda, Lame constant, includes damping
                      !   - property_c(4) : mu, shear modulus, includes damping
                      !   - property_c(5) : E, Young modulus, includes damping
                      !   - property_c(7) : c1, complex primary wave (P-wave) propagation speed (c_p), includes damping (non attenuated)
                      !   - property_c(8) : c2, complex primary wave (S-wave) propagation speed (c_s), includes damping (non attenuated)
                      !   - property_c(9) : xi1, attenuation factor of P-wave
                      !   - property_c(10): xi2, attenuation factor of S-wave
                      case (fbem_bardet)
                        ! Allocate region properties
                        allocate(region(i)%property_r(14))
                        allocate(region(i)%property_c(3:10))
                        ! Read region properties:
                        !   - property_r(1) : rhof, density of the fluid
                        !   - property_r(2) : rhos, density of the solid skeleton
                        !   - property_r(3) : lambda, Lame constant
                        !   - property_r(4) : mu, shear modulus
                        !   - property_r(7) : xi, damping coefficient
                        !   - property_r(8) : phi, porosity
                        !   - property_r(9) : rhoa, additional density (it is not used)
                        !   - property_r(10): R coupling parameter
                        !   - property_r(11): Q coupling parameter
                        !   - property_r(12): b coupling parameter
                        read(fileunit,*) tmp_type, region(i)%property_r( 1), region(i)%property_r( 2), region(i)%property_r( 3),&
                                                   region(i)%property_r( 4), region(i)%property_r( 7), region(i)%property_r( 8),&
                                                   region(i)%property_r( 9), region(i)%property_r(10), region(i)%property_r(11),&
                                                   region(i)%property_r(12)
                        ! Calculate rho1 and rho2
                        region(i)%property_r(13)=(1.0d0-region(i)%property_r(8))*region(i)%property_r(2)
                        region(i)%property_r(14)=region(i)%property_r(8)*region(i)%property_r(1)
                        ! Conversion
                        region(i)%property_r(6)=0.5d0*region(i)%property_r(3)/(region(i)%property_r(3)+region(i)%property_r(4))
                        region(i)%property_c(3)=region(i)%property_r(3)*(1.0d0+c_im*2.0d0*region(i)%property_r(7))
                        region(i)%property_c(4)=region(i)%property_r(4)*(1.0d0+c_im*2.0d0*region(i)%property_r(7))
                        region(i)%property_c(5)=region(i)%property_r(5)*(1.0d0+c_im*2.0d0*region(i)%property_r(7))
                        ! Calculate constant parts of c1 and c2, i.e. c1 and c2 and attenuation factors
                        cte_r   =region(i)%property_r(10)
                        cte_q   =region(i)%property_r(11)
                        cte_b   =region(i)%property_r(12)
                        cte_p   =region(i)%property_c(3)+2.0d0*region(i)%property_c(4)+cte_q**2/cte_r
                        cte_trho=region(i)%property_r(13)+region(i)%property_r(14)
                        cte_phi =region(i)%property_r(8)
                        cte_rhof=region(i)%property_r(1)
                        region(i)%property_c(7)=zsqrt((cte_p+2.0d0*cte_q+cte_r)/cte_trho)
                        region(i)%property_c(8)=zsqrt(region(i)%property_c(4)/cte_trho)
                        region(i)%property_c(9)=cte_trho/cte_b*((cte_q+cte_r)/(cte_p+2.0d0*cte_q+cte_r)*cte_phi*cte_rhof/cte_trho)**2
                        region(i)%property_c(10)=cte_trho/cte_b*(cte_phi*cte_rhof/cte_trho)**2
                      !
                      ! Bougacha-Roesset-Tassoulas viscoelasticity model of poroelasticity
                      !
                      ! List of properties
                      ! ------------------
                      !
                      !   - property_r(1) : rhof, density of the fluid
                      !   - property_r(2) : rhos, density of the solid skeleton
                      !   - property_r(3) : lambda, Lame constant
                      !   - property_r(4) : mu, shear modulus
                      !   - property_r(5) : E, Young modulus
                      !   - property_r(6) : nu, Poisson's ratio
                      !   - property_r(7) : xi, damping coefficient
                      !   - property_r(8) : phi, porosity
                      !   - property_r(9) : rhoa, additional density
                      !   - property_r(10): R coupling parameter
                      !   - property_r(11): Q coupling parameter
                      !   - property_r(12): b coupling parameter
                      !   - property_r(13): rho1, density of the solid skeleton with respect to the total volume
                      !   - property_r(14): rho2, density of the fluid with respect to the total volume
                      !   - property_c(3) : lambda, Lame constant, includes damping
                      !   - property_c(4) : mu, shear modulus, includes damping
                      !   - property_c(5) : E, Young modulus, includes damping
                      case (fbem_brt_cp1,fbem_brt_cp2,fbem_brt_cpm)
                        ! Allocate region properties
                        allocate(region(i)%property_r(14))
                        allocate(region(i)%property_c(3:5))
                        ! Read region properties:
                        !   - property_r(1) : rhof, density of the fluid
                        !   - property_r(2) : rhos, density of the solid skeleton
                        !   - property_r(3) : lambda, Lame constant
                        !   - property_r(4) : mu, shear modulus
                        !   - property_r(7) : xi, damping coefficient
                        !   - property_r(8) : phi, porosity
                        !   - property_r(9) : rhoa, additional density (it is not used)
                        !   - property_r(10): R coupling parameter
                        !   - property_r(11): Q coupling parameter
                        !   - property_r(12): b coupling parameter
                        read(fileunit,*) tmp_type, region(i)%property_r( 1), region(i)%property_r( 2), region(i)%property_r( 3),&
                                                   region(i)%property_r( 4), region(i)%property_r( 7), region(i)%property_r( 8),&
                                                   region(i)%property_r( 9), region(i)%property_r(10), region(i)%property_r(11),&
                                                   region(i)%property_r(12)
                        ! Calculate rho1 and rho2
                        region(i)%property_r(13)=(1.0d0-region(i)%property_r(8))*region(i)%property_r(2)
                        region(i)%property_r(14)=region(i)%property_r(8)*region(i)%property_r(1)
                        ! Conversion
                        region(i)%property_r(6)=0.5d0*region(i)%property_r(3)/(region(i)%property_r(3)+region(i)%property_r(4))
                        region(i)%property_c(3)=region(i)%property_r(3)*(1.0d0+c_im*2.0d0*region(i)%property_r(7))
                        region(i)%property_c(4)=region(i)%property_r(4)*(1.0d0+c_im*2.0d0*region(i)%property_r(7))
                        region(i)%property_c(5)=region(i)%property_r(5)*(1.0d0+c_im*2.0d0*region(i)%property_r(7))

                    end select

                  ! -----------------
                  ! POROELASTIC MEDIA
                  ! -----------------
                  !
                  ! List of properties
                  ! ------------------
                  !
                  !   - property_r(1) : rhof, density of the fluid
                  !   - property_r(2) : rhos, density of the solid skeleton
                  !   - property_r(3) : lambda, Lame constant
                  !   - property_r(4) : mu, shear modulus
                  !   - property_r(5) : E, Young modulus
                  !   - property_r(6) : nu, Poisson's ratio
                  !   - property_r(7) : xi, damping coefficient
                  !   - property_r(8) : phi, porosity
                  !   - property_r(9) : rhoa, additional density
                  !   - property_r(10): R coupling parameter
                  !   - property_r(11): Q coupling parameter
                  !   - property_r(12): b coupling parameter
                  !   - property_r(13): rho1, density of the solid skeleton with respect to the total volume
                  !   - property_r(14): rho2, density of the fluid with respect to the total volume
                  !   - property_c(3) : lambda, Lame constant, includes damping
                  !   - property_c(4) : mu, shear modulus, includes damping
                  !   - property_c(5) : E, Young modulus, includes damping
                  !   - property_c(6) : nu, Poisson's ratio, includes damping
                  !   - property_c(10): R coupling parameter, includes damping
                  !   - property_c(11): Q coupling parameter, includes damping
                  !
                  case (fbem_poroelastic)

                    if (region(i)%space.eq.fbem_half_space) stop 'not halfspace for poroelasticity yet'
                    if (region(i)%space.eq.fbem_multilayered_half_space) stop 'not multilayered halfspace for poroelasticity yet'

                    ! Allocate region properties
                    allocate(region(i)%property_r(14))
                    allocate(region(i)%property_c(3:11))
                    region(i)%property_r=0.d0
                    region(i)%property_c=0.d0

                    if (trim(tmp_type).eq.'poroelastic') then
                      !
                      ! Formato de entrada 1
                      !
                      ! Read region properties:
                      !   - property_r(1) : rhof, density of the fluid
                      !   - property_r(2) : rhos, density of the solid skeleton
                      !   - property_r(3) : lambda, Lame constant
                      !   - property_r(4) : mu, shear modulus
                      !   - property_r(7) : xi, damping coefficient
                      !   - property_r(8) : phi, porosity
                      !   - property_r(9) : rhoa, additional density
                      !   - property_r(10): R coupling parameter
                      !   - property_r(11): Q coupling parameter
                      !   - property_r(12): b coupling parameter
                      read(fileunit,*) tmp_type, tmp_variable, &
                                       region(i)%property_r( 1), region(i)%property_r( 2), region(i)%property_r( 3),&
                                       region(i)%property_r( 4), region(i)%property_r( 7), region(i)%property_r( 8),&
                                       region(i)%property_r( 9), region(i)%property_r(10), region(i)%property_r(11),&
                                       region(i)%property_r(12)
                      ! Calculate rho1 and rho2
                      region(i)%property_r(13)=(1.0d0-region(i)%property_r(8))*region(i)%property_r(2)
                      region(i)%property_r(14)=region(i)%property_r(8)*region(i)%property_r(1)
                      ! Conversion
                      region(i)%property_r( 6)=0.5d0*region(i)%property_r(3)/(region(i)%property_r(3)+region(i)%property_r(4))
                      region(i)%property_c( 3)=region(i)%property_r( 3)*(1.0d0+c_im*2.0d0*region(i)%property_r(7))
                      region(i)%property_c( 4)=region(i)%property_r( 4)*(1.0d0+c_im*2.0d0*region(i)%property_r(7))
                      region(i)%property_c( 5)=region(i)%property_r( 5)*(1.0d0+c_im*2.0d0*region(i)%property_r(7))
                      region(i)%property_c( 6)=0.5d0*region(i)%property_c(3)/(region(i)%property_c(3)+region(i)%property_c(4))
                      region(i)%property_c(10)=region(i)%property_r(10)
                      region(i)%property_c(11)=region(i)%property_r(11)
                    else
                      !
                      ! Formato entrada 2: rhof, rhos, K, mu, phi, rhoa, R, Q, b
                      !
                      ! Read region properties:
                      !   - property_r(1) : rhof, density of the fluid
                      !   - property_r(2) : rhos, density of the solid skeleton
                      !   - property_c(3) : K, bulk modulus (complex) (luego se transforma a lambda)
                      !   - property_c(4) : mu, shear modulus (complex)
                      !   - property_r(8) : phi, porosity
                      !   - property_r(9) : rhoa, additional density
                      !   - property_r(10): R coupling parameter
                      !   - property_r(11): Q coupling parameter
                      !   - property_r(12): b coupling parameter
                      read(fileunit,*) tmp_type, tmp_variable, &
                                       region(i)%property_r( 1), region(i)%property_r( 2), region(i)%property_c( 3),&
                                       region(i)%property_c( 4), region(i)%property_r( 8), region(i)%property_r( 9), &
                                       region(i)%property_r(10), region(i)%property_r(11), region(i)%property_r(12)
                      ! Calculate rho1 and rho2
                      region(i)%property_r(13)=(1.0d0-region(i)%property_r(8))*region(i)%property_r(2)
                      region(i)%property_r(14)=region(i)%property_r(8)*region(i)%property_r(1)
                      ! Real part of K and mu
                      region(i)%property_r( 3)=dreal(region(i)%property_c(3))
                      region(i)%property_r( 4)=dreal(region(i)%property_c(4))
                      ! nu
                      region(i)%property_c( 6)=0.5d0*(3.d0*region(i)%property_c(3)-2.d0*region(i)%property_c(4))/&
                                                     (3.d0*region(i)%property_c(3)+region(i)%property_c(4))
                      region(i)%property_r( 6)=dreal(region(i)%property_c(6))
                      ! Conversion of K to lambda
                      region(i)%property_r( 3)=3.d0*region(i)%property_r(3)*region(i)%property_r(6)/(1.d0+region(i)%property_r(6))
                      region(i)%property_c( 3)=3.d0*region(i)%property_c(3)*region(i)%property_c(6)/(1.d0+region(i)%property_c(6))
                      ! Complex E, R, Q
                      ! Take de dilatational damping
                      region(i)%property_r( 7)=dimag(region(i)%property_c(3))/(2.d0*region(i)%property_r(3))
                      region(i)%property_c( 5)=region(i)%property_r( 5)*(1.0d0+c_im*2.0d0*region(i)%property_r(7))
                      region(i)%property_c(10)=region(i)%property_r(10)*(1.0d0+c_im*2.0d0*region(i)%property_r(7))
                      region(i)%property_c(11)=region(i)%property_r(11)*(1.0d0+c_im*2.0d0*region(i)%property_r(7))
                    end if



                end select

              end if

          end select

        ! ==========================================================================================================================
        ! LAPLACE
        ! ==========================================================================================================================

        case (fbem_laplace)
          ! Allocate region properties
          allocate(region(i)%property_r(1))
          ! Read region properties (real 3D values):
          !   - property_r(1): k, conductivity
          read(fileunit,*) region(i)%property_r(1)

      end select

      ! ==================
      ! READ BE BODY LOADS
      ! ==================

      if (region(i)%class.eq.fbem_be) then
        read(fileunit,*) region(i)%n_be_bodyloads
        allocate(region(i)%be_bodyload(region(i)%n_be_bodyloads))
        backspace(fileunit)
        read(fileunit,*) region(i)%n_be_bodyloads, (region(i)%be_bodyload(j),j=1,region(i)%n_be_bodyloads)
      end if

      ! ====================
      ! READ INCIDENT FIELDS
      ! ====================

      if ((problem%analysis.eq.fbem_harmonic).and.(region(i)%class.eq.fbem_be)) then
        ! Read
        read(fileunit,*) region(i)%n_incidentfields
        allocate(region(i)%incidentfield(region(i)%n_incidentfields))
        backspace(fileunit)
        read(fileunit,*) region(i)%n_incidentfields, (region(i)%incidentfield(j),j=1,region(i)%n_incidentfields)
        ! Convert eid to iid in region(i)%incidentfield(j)
        do j=1,region(i)%n_incidentfields
          ! Check if the boundary is within id limits
          if ((region(i)%incidentfield(j).ge.incidentfield_eid_min).and.&
              (region(i)%incidentfield(j).le.incidentfield_eid_max)) then
            if (incidentfield_iid(region(i)%incidentfield(j)).eq.0) then
              call fbem_error_message(error_unit,0,'region',region(i)%id,'there is an incident field that does not exist.')
            else
              ! Transform eid to iid
              region(i)%incidentfield(j)=incidentfield_iid(region(i)%incidentfield(j))
            end if
          else
            call fbem_error_message(error_unit,0,'region',region(i)%id,'there is an incident field that does not exist.')
          end if
        end do
        ! Check that the type of region of the incident wave field is the same as the region
        do j=1,region(i)%n_incidentfields
          si=region(i)%incidentfield(j)
          if (incidentfield(si)%region_type.ne.region(i)%type) then
            call fbem_error_message(error_unit,0,'region',region(i)%id,'there is an incompatible incident field.')
          end if
        end do
      end if

    end do

    ! ==================================
    ! CHECK AND BUILD REGION IDENTIFIERS
    ! ==================================

    region_eid_min=region(1)%id
    region_eid_max=region(1)%id
    do i=2,n_regions
      if (region(i)%id.lt.region_eid_min) region_eid_min=region(i)%id
      if (region(i)%id.gt.region_eid_max) region_eid_max=region(i)%id
    end do
    allocate (region_iid(region_eid_min:region_eid_max))
    region_iid=0
    do i=1,n_regions
      if (region_iid(region(i)%id).ne.0) then
        call fbem_error_message(error_unit,0,'region',region(i)%id,'is repeated')
      else
        region_iid(region(i)%id)=i
      end if
    end do

    ! =========================================
    ! CALCULATE THE NUMBER OF BE AND FE REGIONS
    ! =========================================

    ! Initialize
    n_be_regions=0
    n_fe_regions=0
    ! Calculate
    do i=1,n_regions
      if (region(i)%class.eq.fbem_be) n_be_regions=n_be_regions+1
      if (region(i)%class.eq.fbem_fe) n_fe_regions=n_fe_regions+1
    end do
    ! Write results
    if (verbose_level.ge.4) then
      write(fmtstr,*) '(3x,a,i',fbem_nchar_int(n_be_regions),')'
      call fbem_trimall(fmtstr)
      write(output_unit,fmtstr) 'Number of BE regions: ', n_be_regions
      write(fmtstr,*) '(3x,a,i',fbem_nchar_int(n_fe_regions),')'
      call fbem_trimall(fmtstr)
      write(output_unit,fmtstr) 'Number of FE regions: ', n_fe_regions
    end if

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section [regions]')

  else

    call fbem_error_message(error_unit,0,'[regions]',0,'this section is required')

  end if

end subroutine read_regions
