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
!! <b> Subroutine that reads the cross sections from a file. </b>
!!
subroutine read_cross_sections(fileunit)

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
  integer                        :: fileunit               ! Unit of the file to read from
  ! Local
  character(len=fbem_stdcharlen) :: section_name           ! Name of the section
  logical                        :: found
  integer                        :: i, j                   ! Counters
  integer                        :: n_cross_sections
  character(len=fbem_stdcharlen) :: tmp_class, tmp_cross_section    ! Class of structural element
  integer                        :: tmp_n_fe_subregions
  integer, allocatable           :: tmp_fe_subregion(:)
  real(kind=real64)              :: area, inertia(6), tmp_fd1(3), ksh(3)
  real(kind=real64)              :: thickness(3), v2ref(3), auxr(10), nu, m, n
  real(kind=real64)              :: mass_matrix(6,6), rho_add, zmin, zmax, zc
  complex(kind=real64)           :: auxc(10)
  complex(kind=real64)           :: Usource
  real(kind=real64)              :: em_pars(3)
  integer                        :: ks, ss
  integer                        :: ke, se
  integer                        :: kn, sn
  logical                        :: valid

  ! Return if not needed
  if (n_fe_regions.eq.0) return

  section_name='cross sections'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! Read the number of cross sections
    read(fileunit,*) n_cross_sections
    if (n_cross_sections.le.0) then
      call fbem_error_message(error_unit,0,section_name,n_cross_sections,'the number of cross sections must be >0')
    end if

    ! Loop through CROSS SECTIONS
    do i=1,n_cross_sections

      ! Initialization
      valid=.false.
      mass_matrix=0
      area=0
      inertia=0
      ksh=0
      v2ref=0
      thickness=0
      tmp_fd1=0

      ! Read
      read(fileunit,*) tmp_class, tmp_n_fe_subregions
      allocate (tmp_fe_subregion(tmp_n_fe_subregions))
      backspace (fileunit)


      !
      ! Masa puntual como elemento
      !
      !
      ! PMASS
      !
      if ((trim(tmp_class).eq.'point_mass').or.(trim(tmp_class).eq.'pmass')) then
        valid=.true.
        read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section
        backspace (fileunit)
        if (trim(tmp_cross_section).eq.'balanced') then
          read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, (mass_matrix(kn,kn),kn=1,3*(problem%n-1))
        else if (trim(tmp_cross_section).eq.'unbalanced') then
          read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, mass_matrix(1:3*(problem%n-1),1:3*(problem%n-1))
        else
          call fbem_error_message(error_unit,0,section_name,0,'invalid type of mass matrix')
        end if
        if (problem%analysis.eq.fbem_static) then
          call fbem_warning_message(error_unit,0,section_name,0,'point_mass elements are ignored in static analysis')
        end if
      end if

      !
      ! DEGBEAM
      !
      !
      ! Falta leer factor de correccion de cortante customizado
      !
      if ((trim(tmp_class).eq.'degenerated_beam').or.(trim(tmp_class).eq.'degbeam')) then
        valid=.true.
        select case (problem%n)
          case (2)
            read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), thickness(2)
            thickness(3)=1. ! Plane strain assumption, if plane stress the other option must be activated
          case (3)
            read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), thickness(2:3), v2ref
        end select
      end if

      !
      ! STRBEAM
      !
      if ((trim(tmp_class).eq.'strbeam_eb').or.(trim(tmp_class).eq.'strbeam_t')) then
        valid=.true.
        ! By default
        select case (problem%n)
          case (2)
            ! area and inertia must take into account if plane strain or plane stress
            !
            ! lo suyo seria asumir rectangle, y introducir canto, y luego el espesor ponerlo automaticamente 1, o el espesor
            ! de plain stress
            !
            read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), area, inertia(3), ksh(2)
            tmp_cross_section='generic'
          case (3)
            read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section
            backspace (fileunit)
            if (trim(tmp_cross_section).eq.'circle') then
              read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, auxr(  1), v2ref ! Diameter
            else if (trim(tmp_cross_section).eq.'hollow_circle') then
              read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, auxr(1:2), v2ref ! Outer and inner diameter
            else if (trim(tmp_cross_section).eq.'rectangle') then
              read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, auxr(1:2), v2ref ! Width (y') and height (z')
            else if (trim(tmp_cross_section).eq.'generic') then
              read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, area, (inertia(j),j=1,3), (ksh(j),j=1,3), v2ref
            else
              call fbem_error_message(error_unit,0,section_name,0,'invalid type of cross section for a 3D strbeam_eb')
            end if
        end select
      end if

      !
      ! STRBEAM_FLOODED
      !
      if (trim(tmp_class).eq.'strbeam_flooded') then
        valid=.true.
        ! en 2D se filtran los elementos segun y, y en 3D segun z
        read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), zmin, zmax, rho_add
        if (zmin.ge.zmax) then
          call fbem_error_message(error_unit,0,section_name,0,'strbeam_flood requires zmin<zmax')
        end if
      end if

      !
      ! STRBEAM_SUBMERGED
      !
      if (trim(tmp_class).eq.'strbeam_submerged') then
        valid=.true.
        ! en 2D se filtran los elementos segun y, y en 3D segun z
        read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), zmin, zmax, rho_add
        if (zmin.ge.zmax) then
          call fbem_error_message(error_unit,0,section_name,0,'strbeam_submerged requires zmin<zmax')
        end if
      end if

      !
      ! BAR
      !
      if (trim(tmp_class).eq.'bar') then
        valid=.true.
        read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section
        backspace (fileunit)
        if (trim(tmp_cross_section).eq.'circle') then
          read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, auxr(  1) ! Diameter
        else if (trim(tmp_cross_section).eq.'hollow_circle') then
          read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, auxr(1:2) ! Outer and inner diameter
        else if (trim(tmp_cross_section).eq.'rectangle') then
          read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, auxr(1:2) ! Width (y') and height (z')
        else if (trim(tmp_cross_section).eq.'generic') then
          read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, area
        else
          call fbem_error_message(error_unit,0,section_name,0,'invalid type of cross section for a 3D strbeam_eb')
        end if
      end if

      !
      ! DISTRA
      !
      !   -distra 2D  :  kx' ky'
      !   -distra 3D  :  kx' ky' kz'
      if (trim(tmp_class).eq.'distra') then
        valid=.true.
        select case (problem%analysis)
          case (fbem_static)
            read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section
            backspace (fileunit)
            if (trim(tmp_cross_section).eq.'global') then
              read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, auxr(1:problem%n)
            else if (trim(tmp_cross_section).eq.'local') then
              select case (problem%n)
                case (2)
                  read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, auxr(1:2)
                case (3)
                  read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, auxr(1:3), v2ref
              end select
            else
              call fbem_error_message(error_unit,0,section_name,0,'invalid type of axes for distra')
            end if
          case (fbem_harmonic)
            read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section
            backspace (fileunit)
            if (trim(tmp_cross_section).eq.'global') then
              read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, auxc(1:problem%n), auxr(1:problem%n)
            else if (trim(tmp_cross_section).eq.'local') then
              select case (problem%n)
                case (2)
                  read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, auxc(1:2), auxr(1:2)
                case (3)
                  read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, auxc(1:3), auxr(1:3), v2ref
              end select
            else
              call fbem_error_message(error_unit,0,section_name,0,'invalid type of axes for distra')
            end if
          case default
            stop 'not yet 123'
        end select
      end if

      !
      ! DISTRA_EM
      !
      !   -distra_em 2D  :  kx' ky' Bl R L Usource
      !   -distra_em 3D  :  kx' ky' kz' Bl R L Usource
      if (trim(tmp_class).eq.'distra_em') then
        valid=.true.
        select case (problem%analysis)
          case (fbem_static)
            stop 'Error : static : distra_em : not available'
          case (fbem_harmonic)
            select case (problem%n)
              case (2)
                read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), auxc(1:2), auxr(1:2), em_pars, Usource
              case (3)
                read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), auxc(1:3), auxr(1:3), v2ref, em_pars, Usource
            end select
          case default
            stop 'Error : unknown analysis'
        end select
      end if

      !
      ! DISROTRA
      !
      !   -disrotra 2D:  kx' ky' krz' ky'rz' -    -    -      -
      !   -disrotra 3D:  kx' ky' kz'  krx'   kry' krz' ky'rz' kz'ry'
      if (trim(tmp_class).eq.'disrotra') then
        valid=.true.
        select case (problem%analysis)
          case (fbem_static)
            read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section
            backspace (fileunit)
            if (trim(tmp_cross_section).eq.'global') then
              select case (problem%n)
                case (2)
                  read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, auxr(1:4)
                case (3)
                  read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, auxr(1:8)
              end select
            else if (trim(tmp_cross_section).eq.'local') then
              select case (problem%n)
                case (2)
                  read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, auxr(1:4)
                case (3)
                  read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, auxr(1:8), v2ref
              end select
            else
              call fbem_error_message(error_unit,0,section_name,0,'invalid type of axes for distra')
            end if
          case (fbem_harmonic)
            read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section
            backspace (fileunit)
            if (trim(tmp_cross_section).eq.'global') then
              select case (problem%n)
                case (2)
                  read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, auxc(1:4), auxr(1:4)
                case (3)
                  read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, auxc(1:8), auxr(1:8)
              end select
            else if (trim(tmp_cross_section).eq.'local') then
              select case (problem%n)
                case (2)
                  read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, auxc(1:4), auxr(1:4)
                case (3)
                  read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_cross_section, auxc(1:8), auxr(1:8), v2ref
              end select
            else
              call fbem_error_message(error_unit,0,section_name,0,'invalid type of axes for disrotra')
            end if
          case default
            stop 'Error : unknown analysis'
        end select
      end if

      !
      ! DISROTRA_EM
      !
      !   -disrotra 2D:  kx' ky' krz' ky'rz' -    -    -      - Bl R L Usource
      !   -disrotra 3D:  kx' ky' kz'  krx'   kry' krz' ky'rz' kz'ry' Bl R L Usource
      if (trim(tmp_class).eq.'disrotra_em') then
        valid=.true.
        select case (problem%analysis)
          case (fbem_static)
            stop 'Error : static : disrotra_em : not available'
          case (fbem_harmonic)
            select case (problem%n)
              case (2)
                read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), auxc(1:4), auxr(1:4), em_pars, Usource
              case (3)
                read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), auxc(1:8), auxr(1:8), v2ref, em_pars, Usource
            end select
          case default
            stop 'Error : unknown analysis'
        end select
      end if

      !
      ! DEGSHELL
      !
      !
      ! Falta leer factor de correccion de cortante customizado
      !
      if ((trim(tmp_class).eq.'shell').or.(trim(tmp_class).eq.'degshell')) then
        valid=.true.
        if (problem%n.eq.2) then
          call fbem_error_message(error_unit,0,section_name,0,'invalid class of cross section for a 2D problem')
        end if
        read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), thickness(3)
      end if

      !
      ! ORTHOTROPIC SHELL ORIENTATION
      !
      if (trim(tmp_class).eq.'orthotropic_shell_orientation') then
        valid=.true.
        if (problem%n.eq.2) then
          call fbem_error_message(error_unit,0,section_name,0,'invalid class of cross section for a 2D problem')
        end if
        read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), tmp_fd1
        tmp_fd1=tmp_fd1/sqrt(dot_product(tmp_fd1,tmp_fd1))
      end if

      !
      ! 2D PLANE STRESS
      !
      if (trim(tmp_class).eq.'plane_stress_thickness') then
        valid=.true.
        if (problem%n.eq.3) then
          call fbem_error_message(error_unit,0,section_name,0,'invalid class of cross section for a 3D problem')
        end if
        if ((problem%n.eq.2).and.(problem%subtype.ne.fbem_mechanics_plane_strain)) then
          call fbem_error_message(error_unit,0,section_name,0,'invalid class of cross section for a 2D problem')
        end if
        read(fileunit,*) tmp_class, tmp_n_fe_subregions, (tmp_fe_subregion(ks),ks=1,tmp_n_fe_subregions), thickness(3)
      end if

      !
      ! CHECK VALIDITY
      !
      if (.not.valid) then
        call fbem_error_message(error_unit,0,section_name,0,'invalid cross section input')
      end if
      !
      ! Check if the indicated FE subregions exist and convert eid to iid in tmp_fe_subregion
      !
      do ks=1,tmp_n_fe_subregions
        if ((tmp_fe_subregion(ks).ge.fe_subregion_eid_min).and.(tmp_fe_subregion(ks).le.fe_subregion_eid_max)) then
          if (fe_subregion_iid(tmp_fe_subregion(ks)).eq.0) then
            call fbem_warning_message(error_unit,0,section_name,0,'wrong definition of a FE subregion in this section')
            call fbem_error_message(error_unit,0,'FE subregion',tmp_fe_subregion(ks),'does not exist')
          else
            tmp_fe_subregion(ks)=fe_subregion_iid(tmp_fe_subregion(ks))
          end if
        else
          call fbem_warning_message(error_unit,0,section_name,0,'wrong definition of a FE subregion in this section')
          call fbem_error_message(error_unit,0,'FE subregion',tmp_fe_subregion(ks),'does not exist')
        end if
      end do

      !
      ! Assign the cross sections to the elements of each FE subregion where applicable
      !
      do ks=1,tmp_n_fe_subregions
        ss=tmp_fe_subregion(ks)
        nu=region(fe_subregion(ss)%region)%property_r(3)
        do ke=1,part(fe_subregion(ss)%part)%n_elements
          se=part(fe_subregion(ss)%part)%element(ke)

          !
          ! PMASS
          !
          if (((trim(tmp_class).eq.'point_mass').or.(trim(tmp_class).eq.'pmass')).and.(element(se)%n_dimension.eq.0)) then
            element(se)%fe_type=0
            if (.not.allocated(element(se)%mass_matrix)) then
              allocate (element(se)%mass_matrix(1:3*(problem%n-1),1:3*(problem%n-1)))
            end if
            element(se)%mass_matrix=mass_matrix(1:3*(problem%n-1),1:3*(problem%n-1))
          end if

          !
          ! DEGBEAM
          !
          if (((trim(tmp_class).eq.'degenerated_beam').or.(trim(tmp_class).eq.'degbeam')).and.(element(se)%n_dimension.eq.1)) then
            element(se)%fe_type=0
            if (.not.allocated(element(se)%v_midnode)) then
              allocate (element(se)%v_midnode(3,3,element(se)%n_nodes))
              element(se)%v_midnode=0
            end if
            if (.not.allocated(element(se)%tn_midnode)) then
              allocate (element(se)%tn_midnode(3,element(se)%n_nodes))
              element(se)%tn_midnode=0
            end if
            element(se)%v_midnode(1,2,:)=v2ref(1)
            element(se)%v_midnode(2,2,:)=v2ref(2)
            element(se)%v_midnode(3,2,:)=v2ref(3)
            element(se)%tn_midnode(2,:)=thickness(2)
            element(se)%tn_midnode(3,:)=thickness(3)
            element(se)%cs_type=3
            element(se)%cs_param(1)=thickness(2)
            element(se)%cs_param(2)=thickness(3)
            element(se)%A=thickness(2)*thickness(3)
            element(se)%I(2)=thickness(3)**3*thickness(2)/12.d0
            element(se)%I(3)=thickness(2)**3*thickness(3)/12.d0
            element(se)%I(1)=element(se)%I(2)+element(se)%I(3)
            element(se)%I(4)=0
            element(se)%I(5)=0
            element(se)%I(6)=0
            element(se)%ksh(1)=0.d0
            element(se)%ksh(2)=10.d0*(1.d0+nu)/(12.d0+11.d0*nu)
            element(se)%ksh(3)=element(se)%ksh(2)
            element(se)%K_intmode=2
            element(se)%K_intngp =0
            element(se)%M_intmode=0
            element(se)%M_intngp =0
            element(se)%Q_intmode=0
            element(se)%Q_intngp =0
            element(se)%r_integration=(thickness(2)+thickness(3))/c_pi
          end if

          !
          ! STRBEAM
          !
          if (((trim(tmp_class).eq.'strbeam_eb').or.(trim(tmp_class).eq.'strbeam_t')).and.(element(se)%n_dimension.eq.1)) then
            if (.not.allocated(element(se)%ep)) then
              allocate (element(se)%ep(problem%n,problem%n))
              element(se)%ep=0.d0
            end if
            if (trim(tmp_class).eq.'strbeam_eb') then
              element(se)%fe_type=1
            else
              element(se)%fe_type=2
            end if
            if (trim(tmp_cross_section).eq.'circle') then
              element(se)%cs_type=1
              element(se)%cs_param(1)=auxr(1)
              element(se)%A=c_pi*(0.5d0*auxr(1))**2
              element(se)%I(1)=0.5d0*c_pi*(0.5d0*auxr(1))**4
              element(se)%I(2)=0.25d0*c_pi*(0.5d0*auxr(1))**4
              element(se)%I(3)=element(se)%I(2)
              element(se)%I(4)=0.d0
              element(se)%I(5)=0.d0
              element(se)%I(6)=0.d0
              element(se)%ksh(1)=1
              element(se)%ksh(2)=6.d0*(1.d0+nu)/(7.d0+6.d0*nu)
              element(se)%ksh(3)=element(se)%ksh(2)
              element(se)%r_integration=0.5d0*auxr(1)
            else if (trim(tmp_cross_section).eq.'hollow_circle') then
              element(se)%cs_type=2
              element(se)%cs_param(1)=auxr(1)
              element(se)%cs_param(2)=auxr(2)
              element(se)%A=c_pi*(0.5d0*auxr(1))**2-c_pi*(0.5d0*auxr(2))**2
              element(se)%I(1)=0.5d0*c_pi*(0.5d0*auxr(1))**4-0.5d0*c_pi*(0.5d0*auxr(2))**4
              element(se)%I(2)=0.25d0*c_pi*(0.5d0*auxr(1))**4-0.25d0*c_pi*(0.5d0*auxr(2))**4
              element(se)%I(3)=element(se)%I(2)
              element(se)%I(4)=0.d0
              element(se)%I(5)=0.d0
              element(se)%I(6)=0.d0
              element(se)%ksh(1)=1
              m=auxr(2)/auxr(1)
              element(se)%ksh(2)=6*(1+nu)*(1+m**2)**2/((7+6*nu)*(1+m**2)**2+(20+12*nu)*m**2)
              element(se)%ksh(3)=element(se)%ksh(2)
              element(se)%r_integration=0.5d0*auxr(1)
            else if (trim(tmp_cross_section).eq.'rectangle') then
              element(se)%cs_type=3
              element(se)%cs_param(1)=auxr(1)
              element(se)%cs_param(2)=auxr(2)
              element(se)%A=auxr(1)*auxr(2)
              element(se)%I(2)=auxr(1)*auxr(2)**3/12.d0
              element(se)%I(3)=auxr(1)**3*auxr(2)/12.d0
              element(se)%I(1)=element(se)%I(2)+element(se)%I(3)
              element(se)%I(4)=0.d0
              element(se)%I(5)=0.d0
              element(se)%I(6)=0.d0
              if (auxr(1).ge.auxr(2)) then
                m=0.5d0*auxr(1)
                n=0.5d0*auxr(2)
              else
                m=0.5d0*auxr(2)
                n=0.5d0*auxr(1)
              end if
              element(se)%ksh(1)=m*n**3*(16.d0/3.d0-3.36d0*n/m*(1.d0-(n/m)**4/12.d0))/element(se)%I(1)
              element(se)%ksh(2)=10.d0*(1.d0+nu)/(12.d0+11.d0*nu)
              element(se)%ksh(3)=element(se)%ksh(2)
              element(se)%r_integration=0.5d0*(sqrt(element(se)%I(2)/element(se)%A)+sqrt(element(se)%I(3)/element(se)%A)) ! radius of gyration
            else if (trim(tmp_cross_section).eq.'generic') then
              element(se)%cs_type=0
              element(se)%A=area
              element(se)%I=inertia
              element(se)%ksh=ksh
              element(se)%r_integration=0.5d0*(sqrt(inertia(2)/area)+sqrt(inertia(3)/area)) ! radius of gyration
            end if
            ! Reference vector for y' axis
            element(se)%ep(:,2)=v2ref(1:problem%n)
            ! Force a straight element
            if (element(se)%type.eq.fbem_line3) then
              sn=element(se)%node(3)
              node(sn)%x=0.5d0*(node(element(se)%node(1))%x+node(element(se)%node(2))%x)
            end if
          end if

          !
          ! STRBEAM_FLOODED
          !
          if ((trim(tmp_class).eq.'strbeam_flooded').and.(element(se)%n_dimension.eq.1)) then
            if (element(se)%cs_type.ne.2) then
              call fbem_error_message(error_unit,0,section_name,0,'strbeam_flooded is available only for hollow_circle cross-sections')
            end if
            ! Centroid of the element
            zc=0.5d0*(node(element(se)%node(1))%x(problem%n)+node(element(se)%node(2))%x(problem%n))
            ! Check if inside de defined region
            if ((zc.le.zmax).and.(zc.ge.zmin)) then
              element(se)%flooded=.true.
              !
              ! Nota: se considera que longitudinalmente el fluido tambien aporta masa (esto es discutible segun en que casos)
              !
              ! Masa dentro del tubo
              !
              element(se)%A_flooded(1)=0.25d0*c_pi*element(se)%cs_param(2)**2
              element(se)%A_flooded(2)=0.25d0*c_pi*element(se)%cs_param(2)**2
              element(se)%A_flooded(3)=0.25d0*c_pi*element(se)%cs_param(2)**2
              element(se)%rho_flooded=rho_add
            end if
          end if

          !
          ! STRBEAM_SUBMERGED
          !
          if ((trim(tmp_class).eq.'strbeam_submerged').and.(element(se)%n_dimension.eq.1)) then
            if (element(se)%cs_type.ne.2) then
              call fbem_error_message(error_unit,0,section_name,0,'strbeam_submerged is available only for hollow_circle cross-sections')
            end if
            ! Centroid of the element
            zc=0.5d0*(node(element(se)%node(1))%x(problem%n)+node(element(se)%node(2))%x(problem%n))
            ! Check if inside de defined region
            if ((zc.le.zmax).and.(zc.ge.zmin)) then
              element(se)%submerged=.true.
              !
              ! Nota: se considera que longitudinalmente no hay masa añadida debido a las fuerzas hidrodinamicas
              !
              element(se)%A_submerged(1)=0.d0
              element(se)%A_submerged(2)=0.25d0*c_pi*element(se)%cs_param(1)**2
              element(se)%A_submerged(3)=0.25d0*c_pi*element(se)%cs_param(1)**2
              element(se)%rho_submerged=rho_add
            end if
          end if

          !
          ! BAR
          !
          if (trim(tmp_class).eq.'bar'.and.(element(se)%n_dimension.eq.1)) then
            element(se)%fe_type=3
            if (trim(tmp_cross_section).eq.'circle') then
              element(se)%cs_type=1
              element(se)%cs_param(1)=auxr(1)
              element(se)%A=c_pi*(0.5d0*auxr(1))**2
            else if (trim(tmp_cross_section).eq.'hollow_circle') then
              element(se)%cs_type=2
              element(se)%cs_param(1)=auxr(1)
              element(se)%cs_param(2)=auxr(2)
              element(se)%A=c_pi*(0.5d0*auxr(1))**2-c_pi*(0.5d0*auxr(2))**2
            else if (trim(tmp_cross_section).eq.'rectangle') then
              element(se)%cs_type=3
              element(se)%cs_param(1)=auxr(1)
              element(se)%cs_param(2)=auxr(2)
              element(se)%A=auxr(1)*auxr(2)
            else if (trim(tmp_cross_section).eq.'generic') then
              element(se)%cs_type=0
              element(se)%A=area
            end if
            if (element(se)%type.ne.fbem_line2) then
              call fbem_error_message(error_unit,0,section_name,0,'bar elements can be only of mesh element type line2')
            end if
          end if

          !
          ! DISTRA
          !
          if ((trim(tmp_class).eq.'distra').and.(element(se)%n_dimension.eq.1)) then
            element(se)%fe_type=4
            if (.not.allocated(element(se)%ep)) then
              allocate (element(se)%ep(problem%n,problem%n))
              element(se)%ep=0.d0
            end if
            element(se)%k_r=0
            element(se)%k_c=0
            element(se)%c=0
            select case (problem%analysis)
              case (fbem_static)
                element(se)%k_r(1:problem%n)=auxr(1:problem%n)
              case (fbem_harmonic)
                element(se)%k_c(1:problem%n)=auxc(1:problem%n)
                element(se)%c  (1:problem%n)=auxr(1:problem%n)
              case default
                stop 'not yet 123'
            end select
            if (trim(tmp_cross_section).eq.'global') then
              element(se)%ep(1,1)=1
            else
              ! Reference vector for y' axis
              element(se)%ep(:,2)=v2ref(1:problem%n)
            end if
            !if (element(se)%type.ne.fbem_line2) then
            !  call fbem_error_message(error_unit,0,section_name,0,'distra elements can be only of mesh element type line2')
            !end if
          end if

          !
          ! DISTRA_EM
          !
          if ((trim(tmp_class).eq.'distra_em').and.(element(se)%n_dimension.eq.1)) then
            element(se)%fe_type=4
            element(se)%fe_options(1)=1
            if (.not.allocated(element(se)%ep)) then
              allocate (element(se)%ep(problem%n,problem%n))
              element(se)%ep=0.d0
            end if
            element(se)%k_r=0
            element(se)%k_c=0
            element(se)%c=0
            select case (problem%analysis)
              case (fbem_static)
                element(se)%k_r(1:problem%n)=auxr(1:problem%n)
              case (fbem_harmonic)
                element(se)%k_c(1:problem%n)=auxc(1:problem%n)
                element(se)%c  (1:problem%n)=auxr(1:problem%n)
              case default
                stop 'not yet 123'
            end select
            ! Reference vector for y' axis
            element(se)%ep(:,2)=v2ref(1:problem%n)
            ! Electromagnetic parameters
            element(se)%em_U=Usource
            element(se)%em_Bl=em_pars(1)
            element(se)%em_R=em_pars(2)
            element(se)%em_L=em_pars(3)
            !if (element(se)%type.ne.fbem_line2) then
            !  call fbem_error_message(error_unit,0,section_name,0,'distra elements can be only of mesh element type line2')
            !end if
          end if

          !
          ! DISROTRA
          !
          if ((trim(tmp_class).eq.'disrotra').and.(element(se)%n_dimension.eq.1)) then
            element(se)%fe_type=5
            if (.not.allocated(element(se)%ep)) then
              allocate (element(se)%ep(problem%n,problem%n))
              element(se)%ep=0.d0
            end if
            element(se)%k_r=0
            element(se)%k_c=0
            element(se)%c=0
            select case (problem%analysis)
              case (fbem_static)
                select case (problem%n)
                  case (2)
                    element(se)%k_r(1:4)=auxr(1:4)
                  case (3)
                    element(se)%k_r(1:8)=auxr(1:8)
                end select
              case (fbem_harmonic)
                select case (problem%n)
                  case (2)
                    element(se)%k_c(1:4)=auxc(1:4)
                    element(se)%c  (1:4)=auxr(1:4)
                  case (3)
                    element(se)%k_c(1:8)=auxc(1:8)
                    element(se)%c  (1:8)=auxr(1:8)
                end select

              case default
                stop 'not yet 123'
            end select
            if (trim(tmp_cross_section).eq.'global') then
              element(se)%ep(1,1)=1
            else
              ! Reference vector for y' axis
              element(se)%ep(:,2)=v2ref(1:problem%n)
            end if
            !if (element(se)%type.ne.fbem_line2) then
            !  call fbem_error_message(error_unit,0,section_name,0,'disrotra elements can be only of mesh element type line2')
            !end if
          end if


          !
          ! DISROTRA_EM
          !
          if ((trim(tmp_class).eq.'disrotra_em').and.(element(se)%n_dimension.eq.1)) then
            element(se)%fe_type=5
            element(se)%fe_options(1)=1
            if (.not.allocated(element(se)%ep)) then
              allocate (element(se)%ep(problem%n,problem%n))
              element(se)%ep=0.d0
            end if
            element(se)%k_r=0
            element(se)%k_c=0
            element(se)%c=0
            select case (problem%analysis)
              case (fbem_static)
                select case (problem%n)
                  case (2)
                    element(se)%k_r(1:4)=auxr(1:4)
                  case (3)
                    element(se)%k_r(1:8)=auxr(1:8)
                end select
              case (fbem_harmonic)
                select case (problem%n)
                  case (2)
                    element(se)%k_c(1:4)=auxc(1:4)
                    element(se)%c  (1:4)=auxr(1:4)
                  case (3)
                    element(se)%k_c(1:8)=auxc(1:8)
                    element(se)%c  (1:8)=auxr(1:8)
                end select

              case default
                stop 'not yet 123'
            end select
            ! Reference vector for y' axis
            element(se)%ep(:,2)=v2ref(1:problem%n)
            ! Electromagnetic parameters
            element(se)%em_U=Usource
            element(se)%em_Bl=em_pars(1)
            element(se)%em_R=em_pars(2)
            element(se)%em_L=em_pars(3)
            !if (element(se)%type.ne.fbem_line2) then
            !  call fbem_error_message(error_unit,0,section_name,0,'disrotra elements can be only of mesh element type line2')
            !end if
          end if


          !
          ! DEGSHELL
          !
          if (((trim(tmp_class).eq.'shell').or.(trim(tmp_class).eq.'degshell')).and.(element(se)%n_dimension.eq.2)) then
            element(se)%fe_type=0
            if (allocated(element(se)%tn_midnode).eqv.(.false.)) then
              allocate (element(se)%tn_midnode(3,element(se)%n_nodes))
              element(se)%tn_midnode=0.
            end if
            element(se)%tn_midnode(3,:)=thickness(3)
            element(se)%ksh(1)=5.d0/6.d0
            element(se)%ksh(2)=5.d0/6.d0
            element(se)%ksh(3)=0.d0
            element(se)%mitc=.true.
            element(se)%K_intmode=0
            element(se)%K_intngp =0
            element(se)%M_intmode=0
            element(se)%M_intngp =0
            element(se)%Q_intmode=0
            element(se)%Q_intngp =0
          end if

          !
          ! ORTHOTROPIC SHELL ORIENTATION
          !
          if ((trim(tmp_class).eq.'orthotropic_shell_orientation').and.(element(se)%n_dimension.eq.2)) then
            element(se)%orthotropic_shell_fd1=tmp_fd1
          end if

          !
          ! 2D PLANE STRESS THICKNESS
          !
          if ((trim(tmp_class).eq.'plane_stress_thickness').and.(element(se)%n_dimension.eq.2)) then
            if (allocated(element(se)%tn_midnode).eqv.(.false.)) then
              allocate (element(se)%tn_midnode(3,element(se)%n_nodes))
              element(se)%tn_midnode=0.
            end if
            element(se)%tn_midnode(3,:)=thickness(3)
          end if

        end do
      end do

      deallocate (tmp_fe_subregion)

    end do ! Loop through CROSS SECTIONS

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  end if

end subroutine read_cross_sections
