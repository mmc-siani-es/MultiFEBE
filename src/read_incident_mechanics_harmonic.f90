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

!! Section format:
!! >
!! <n>
!! [for each incident field:
!! <id>
!! <class>
!! <space> [if space=half-space: <np> <xp> <bc>]
!!
!! <variable> <amplitude> <x0(1)> <x0(2)> [if 3D: <x0(3)>] [if 3D: <varphi>] <theta>
!!
!! Propuesta para siguiente version
!! la nueva linea anterior se divide en dos:
!! <variable> <component> <x0(1)> <x0(2)> [if 3D: <x0(3)>] <n0(1)> <n0(2)> [if 3D: <n0(3)>]
!! [if 3D: <varphi>] <theta>
!!
!! <xs(1)> <xs(2)> [if 3D: <xs(3)>] <symconf(1)> <symconf(2)> [if 3D: <symconf(3)>]
!! <region_type> <wave_type>
!! ]
!! >

subroutine read_incident_mechanics_harmonic(fileunit)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_symmetry
  use fbem_data_structures
  use fbem_harela_incident_field
  use fbem_harpor_incident_field

  ! Problem variables module
  use problem_variables

  ! No implicit variables are allowed
  implicit none

  ! I/O
  integer                        :: fileunit    !! Input file unit
  ! Local
  integer                        :: i, j, kc, km
  character(len=fbem_stdcharlen) :: tmp_class
  character(len=fbem_stdcharlen) :: tmp_space
  character(len=fbem_stdcharlen) :: tmp_regiontype
  character(len=fbem_stdcharlen) :: tmp_wavetype
  logical                        :: found

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section [incident waves]')
  call fbem_search_section(fileunit,'incident waves',found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section [incident waves]')

    ! Read the number of incident fields
    read(fileunit,*) n_incidentfields
    if (n_incidentfields.gt.0) then
      ! Allocate the vector
      allocate (incidentfield(n_incidentfields))
      do i=1,n_incidentfields
        allocate (incidentfield(i)%x0(problem%n),incidentfield(i)%xs(problem%n),incidentfield(i)%symconf(problem%n))
      end do
    end if

    ! Loop through the incident fields
    do i=1,n_incidentfields
      !
      ! Read external identifier
      !
      read(fileunit,*) incidentfield(i)%id
      if (incidentfield(i)%id.le.0) then
        call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,'the identifier must be greater than 0')
      end if
      !
      ! Read the class of incident field
      !
      read(fileunit,*) tmp_class
      ! Assign value
      incidentfield(i)%class=0
      if (trim(tmp_class).eq.'point') incidentfield(i)%class=fbem_point
      if (trim(tmp_class).eq.'line' ) incidentfield(i)%class=fbem_line
      if (trim(tmp_class).eq.'plane') incidentfield(i)%class=fbem_plane
      if (incidentfield(i)%class.eq.0) then
        call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,'the class of incident wave field is invalid.')
      end if
      !
      ! Read the type of space
      !
      read(fileunit,*) tmp_space
      ! Assign value
      incidentfield(i)%space=0
      if (trim(tmp_space).eq.'full-space') incidentfield(i)%space=fbem_full_space
      if (trim(tmp_space).eq.'half-space') incidentfield(i)%space=fbem_half_space
      if (trim(tmp_space).eq.'multilayered_half-space') incidentfield(i)%space=fbem_multilayered_half_space
      if (incidentfield(i)%space.eq.0) then
        call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,'the type of space is invalid.')
      end if
      !
      ! Homogeneous half-space
      !
      if (incidentfield(i)%space.eq.fbem_half_space) then
        ! Read the plane axis, coordinate, and boundary condition
        backspace(fileunit)
        read(fileunit,*) tmp_space, incidentfield(i)%np, incidentfield(i)%xp, incidentfield(i)%bc
        if ((incidentfield(i)%np.lt.1).or.(incidentfield(i)%np.gt.problem%n)) then
          call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,'1<=np<=dimension of the problem')
        end if
      end if
      !
      ! Multilayered half-space
      !
      if (incidentfield(i)%space.eq.fbem_multilayered_half_space) then
        ! Read layers
        backspace(fileunit)
        read(fileunit,*) tmp_space, incidentfield(i)%np, incidentfield(i)%n_layers
        if ((incidentfield(i)%np.eq.0).or.(abs(incidentfield(i)%np).gt.3)) then
          call fbem_error_message(error_unit,0,'incidentfield',incidentfield(i)%id,'np must be +1, -1, +2, -2, +3 or -3')
        end if
        if (incidentfield(i)%n_layers.lt.1) then
          call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,'n_layers must be >=1')
        end if
        backspace(fileunit)
        allocate (incidentfield(i)%layer_ztop(incidentfield(i)%n_layers))
        allocate (incidentfield(i)%layer_material(incidentfield(i)%n_layers))
        allocate (incidentfield(i)%layer_ctype(incidentfield(i)%n_layers))
        allocate (incidentfield(i)%layer_cvalue(incidentfield(i)%n_layers))
        read(fileunit,*) tmp_space, incidentfield(i)%np, incidentfield(i)%n_layers, &
                         (incidentfield(i)%layer_ztop(j),j=1,incidentfield(i)%n_layers), &
                         (incidentfield(i)%layer_material(j),j=1,incidentfield(i)%n_layers), &
                         (incidentfield(i)%layer_ctype(j),j=1,incidentfield(i)%n_layers), &
                         (incidentfield(i)%layer_cvalue(j),j=1,incidentfield(i)%n_layers)
        ! Check layer top coordinates and their order (general for np=+-1,+-2,+-3)
        if (incidentfield(i)%np.gt.0) then
          do j=2,incidentfield(i)%n_layers
            if ((incidentfield(i)%layer_ztop(j)-incidentfield(i)%layer_ztop(j-1)).ge.0.d0) then
              call fbem_error_message(error_unit,0,'incidentfield',incidentfield(i)%id,&
              'the top coordinates of layers must be monotically decreasing (starting from free-surface coordinate inwards)')
            end if
          end do
        else if (incidentfield(i)%np.lt.0) then
          do j=2,incidentfield(i)%n_layers
            if ((incidentfield(i)%layer_ztop(j)-incidentfield(i)%layer_ztop(j-1)).le.0.d0) then
              call fbem_error_message(error_unit,0,'incidentfield',incidentfield(i)%id,&
              'the top coordinates of layers must be monotically increasing (starting from free-surface coordinate inwards)')
            end if
          end do
        end if
      end if
      !
      ! Read the term in which the field is written (0 (primary), 1 (secondary)), amplitude, origin of the incident wave,
      ! angle varphi (if 3D), and angle theta (in degrees)
      !
      ! Depending on the problem dimension
      select case (problem%n)
        case (2)
          read(fileunit,*) incidentfield(i)%variable, incidentfield(i)%amplitude, (incidentfield(i)%x0(kc),kc=1,2),&
                                 incidentfield(i)%theta
          incidentfield(i)%varphi=0.d0
        case (3)
          read(fileunit,*) incidentfield(i)%variable, incidentfield(i)%amplitude, (incidentfield(i)%x0(kc),kc=1,3),&
                                 incidentfield(i)%varphi, incidentfield(i)%theta
      end select
      ! Convert degrees to rad
      incidentfield(i)%varphi=incidentfield(i)%varphi*c_pi/180.d0
      incidentfield(i)%theta=incidentfield(i)%theta*c_pi/180.d0
      !
      ! Read the symmetry/anti-symmetry decomposition configuration
      !
      read(fileunit,*) (incidentfield(i)%xs(kc),kc=1,problem%n), (incidentfield(i)%symconf(kc),kc=1,problem%n)
      ! Check the symmetry configuration
      if (incidentfield(i)%space.eq.fbem_half_space) then
        if (incidentfield(i)%symconf(incidentfield(i)%np).ne.0) then
          call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,&
                                  'the symmetry/anti-symmetry decomposition is not valid since symconf(np) must be 0')
        end if
      end if
      !
      ! Read the type of region and the type of wave
      !
      read(fileunit,*) tmp_regiontype, tmp_wavetype
      ! Assign values
      incidentfield(i)%region_type=0
      if (trim(tmp_regiontype).eq.'fluid'       ) incidentfield(i)%region_type=fbem_potential
      if (trim(tmp_regiontype).eq.'elastic'     ) incidentfield(i)%region_type=fbem_viscoelastic
      if (trim(tmp_regiontype).eq.'viscoelastic') incidentfield(i)%region_type=fbem_viscoelastic
      if (trim(tmp_regiontype).eq.'poroelastic' ) incidentfield(i)%region_type=fbem_poroelastic
      if (incidentfield(i)%region_type.eq.0) then
        call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,&
        'the region type of the incident wave field must be "fluid", "elastic", "viscoelastic" or "poroelastic"')
      end if
      incidentfield(i)%wave_type=0
      ! Switch depending on the type of region
      select case (incidentfield(i)%region_type)
        !
        ! Inviscid fluid
        !
        case (fbem_potential)
          if (trim(tmp_wavetype).eq.'p') incidentfield(i)%wave_type=1
          if (incidentfield(i)%wave_type.eq.0) then
            call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,&
                                    'the wave type for a fluid can be only "p".')
          end if
        !
        ! Viscoelastic solid
        !
        case (fbem_viscoelastic)
          if (trim(tmp_wavetype).eq.'p'       ) incidentfield(i)%wave_type=fbem_harela_p_wave
          if (trim(tmp_wavetype).eq.'sv'      ) incidentfield(i)%wave_type=fbem_harela_sv_wave
          if (trim(tmp_wavetype).eq.'sh'      ) incidentfield(i)%wave_type=fbem_harela_sh_wave
          if (trim(tmp_wavetype).eq.'rayleigh') incidentfield(i)%wave_type=fbem_harela_rayleigh_wave
          if (incidentfield(i)%wave_type.eq.0) then
            call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,&
                                    'the wave type for a viscoelastic solid can be "p", "sv", "sh" or "rayleigh".')
          end if
        !
        ! Poroelastic medium
        !
        case (fbem_poroelastic)
          if (trim(tmp_wavetype).eq.'p1'      ) incidentfield(i)%wave_type=fbem_harpor_p1_wave
          if (trim(tmp_wavetype).eq.'p2'      ) incidentfield(i)%wave_type=fbem_harpor_p2_wave
          if (trim(tmp_wavetype).eq.'shx'     ) incidentfield(i)%wave_type=fbem_harpor_shx_wave
          if (trim(tmp_wavetype).eq.'shy'     ) incidentfield(i)%wave_type=fbem_harpor_shy_wave
          if (trim(tmp_wavetype).eq.'rayleigh') incidentfield(i)%wave_type=fbem_harpor_R_wave_permeable
          if (incidentfield(i)%wave_type.eq.0) then
            call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,&
                                    'the wave type for a poroelastic medium can be "p1", "shx", "shy" or "rayleigh".')
          end if
      end select
      !
      ! Check if implemented
      !
      ! Switch depending on the class of source
      select case (incidentfield(i)%class)

        ! ---------- !
        ! Point wave !
        ! ---------- !

        case (fbem_point)
          ! nothing to do

        ! --------- !
        ! Line wave !
        ! --------- !

        case (fbem_line)
          stop 'not implemented yet'

        ! ---------- !
        ! Plane wave !
        ! ---------- !

        case (fbem_plane)
          ! Switch depending on the type of region
          select case (incidentfield(i)%region_type)
            !
            ! Inviscid fluid
            !
            case (fbem_potential)
              ! The multilayered half-space has not been implemented
              if (incidentfield(i)%space.eq.fbem_multilayered_half_space) then
                call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,&
                                        'multilayered half-space for inviscid fluids has not been implemented')
              end if
            !
            ! Viscoelastic solid
            !
            case (fbem_viscoelastic)
              ! Some limitations exist:
              !   - x0 and xs must be 0.,0.,0.
              !   - if 2D: np=2, bc=1, wavetype=sh=>use the fluid model and make a variable transformation
              !   - if 3D: np=3, bc=1
              !   - symconf
              do kc=1,problem%n
                if (incidentfield(i)%x0(kc).ne.0.d0) then
                  call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,'all components of x0 can be only 0.')
                end if
                if (incidentfield(i)%xs(kc).ne.0.d0) then
                  call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,'all components of xs can be only 0.')
                end if
              end do
              ! Switch depending on the dimension of the problem
              select case (problem%n)
                ! 2D
                case (2)
                  if (incidentfield(i)%space.eq.fbem_half_space) then
                    if (incidentfield(i)%np.ne.2) then
                      call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,'np can be only 2.')
                    end if
                    if (incidentfield(i)%bc.ne.1) then
                      call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,'bc must be 1.')
                    end if
                  end if
                  if (incidentfield(i)%wave_type.eq.fbem_harela_sh_wave) then
                    call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,&
                                            'for anti-plane problems, use and equivalent fluid model.')
                  end if
                  if (nint(incidentfield(i)%symconf(2)).ne.0) then
                    call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,&
                                            'the symmetry/anti-symmetry decomposition can not be done for the y direction')
                  end if
                  incidentfield(i)%symconf(2)=incidentfield(i)%symconf(1)

                  ! Dudas con respecto a lo de la simetria del campo incidente....

                ! 3D
                case (3)
                  if (incidentfield(i)%space.eq.fbem_half_space) then
                    if (incidentfield(i)%np.ne.3) then
                      call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,'np can be only 2.')
                    end if
                    if (incidentfield(i)%bc.ne.1) then
                      call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,'bc must be 1.')
                    end if
                  end if
                  if (nint(incidentfield(i)%symconf(1)).ne.0) then
                    call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,&
                                            'the symmetry/anti-symmetry decomposition can not be done for the x direction')
                  end if
                  if (nint(incidentfield(i)%symconf(3)).ne.0) then
                    call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,&
                                            'the symmetry/anti-symmetry decomposition can not be done for the z direction')
                  end if
              end select

              ! =======================
              ! Multilayered half-space
              ! =======================

              if (incidentfield(i)%space.eq.fbem_multilayered_half_space) then
                ! solo en dirección z positiva o negativa
                if (abs(incidentfield(i)%np).ne.3) then
                  call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,'np must be 3 or -3, i.e. free-surface with n=(0,0,1) or n=(0,0,-1)')
                end if
                ! Solo incidencia vertical
                if (abs(incidentfield(i)%theta-c_pi_2).gt.1.d-12) then
                  call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,'only vertical incidence is implemented')
                end if
                ! Solo sh
                if (abs(incidentfield(i)%varphi).gt.1.d-12) then
                  call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,'only varphi=0 is implemented')
                end if

                ! =================================================================
                ! CHANGE EID TO IID IN THE INCIDENT FIELD -> MATERIALS CONNECTIVITY
                ! =================================================================

                do j=1,incidentfield(i)%n_layers
                  if ((incidentfield(i)%layer_material(j).ge.material_eid_min).and.(incidentfield(i)%layer_material(j).le.material_eid_max)) then
                    if (material_iid(incidentfield(i)%layer_material(j)).eq.0) then
                      call fbem_warning_message(error_unit,0,'material',incidentfield(i)%layer_material(j),'does not exist')
                      call fbem_error_message(error_unit,0,'incidentfield',incidentfield(i)%id,'contain the previous material')
                    else
                      incidentfield(i)%layer_material(j)=material_iid(incidentfield(i)%layer_material(j))
                    end if
                  else
                    call fbem_warning_message(error_unit,0,'material',incidentfield(i)%layer_material(j),'does not exist')
                    call fbem_error_message(error_unit,0,'incidentfield',incidentfield(i)%id,'contain the previous material')
                  end if
                end do
                ! Only isotropic_linear_elastic_solid are allowed
                do j=1,incidentfield(i)%n_layers
                  if (trim(material(incidentfield(i)%layer_material(j))%type).ne.'elastic_solid') then
                    call fbem_warning_message(error_unit,0,'material',material(incidentfield(i)%layer_material(j))%id,&
                                              'is not a elastic_solid.')
                    call fbem_error_message(error_unit,0,'incidentfield',incidentfield(i)%id,'contain the previous material')
                  end if
                end do

                !
                ! Allocate and copy material properties
                !

                allocate (incidentfield(i)%layer_lambda(incidentfield(i)%n_layers))
                allocate (incidentfield(i)%layer_mu(incidentfield(i)%n_layers))
                allocate (incidentfield(i)%layer_rho(incidentfield(i)%n_layers))
                do j=1,incidentfield(i)%n_layers
                  km=incidentfield(i)%layer_material(j)
                  ! Copy lambda
                  if (material(km)%property_defined(3)) then
                    if (material(km)%property_defined(7)) then
                      incidentfield(i)%layer_lambda(j)=material(km)%property(3,1)*(1.d0+c_im*2.d0*material(km)%property(7,1))
                    else
                      incidentfield(i)%layer_lambda(j)=material(km)%property(3,1)
                    end if
                  else
                    call fbem_warning_message(error_unit,0,'material',material(km)%id,'lambda is required')
                    call fbem_error_message(error_unit,0,'incidentfield',incidentfield(i)%id,'contain the previous material')
                  end if
                  ! Copy mu
                  if (material(km)%property_defined(4)) then
                    if (material(km)%property_defined(7)) then
                      incidentfield(i)%layer_mu(j)=material(km)%property(4,1)*(1.d0+c_im*2.d0*material(km)%property(7,1))
                    else
                      incidentfield(i)%layer_mu(j)=material(km)%property(4,1)
                    end if
                  else
                    call fbem_warning_message(error_unit,0,'material',material(km)%id,'mu is required')
                    call fbem_error_message(error_unit,0,'incidentfield',incidentfield(i)%id,'contain the previous material')
                  end if
                  ! Copy rho
                  if (material(km)%property_defined(6)) then
                    incidentfield(i)%layer_rho(j)=material(km)%property(6,1)
                  else
                    call fbem_warning_message(error_unit,0,'material',material(km)%id,'rho is required')
                    call fbem_error_message(error_unit,0,'incidentfield',incidentfield(i)%id,'contain the previous material')
                  end if
                end do

              end if

            ! =======================
            ! BIOT POROELASTIC MEDIUM
            ! =======================

            case (fbem_poroelastic)
              ! Some limitations exist:
              !   - x0 and xs must be 0.,0.,0.
              !   - if 2D: np=2, bc=1, shy does not make sense
              !   - if 3D: np=3, bc=1
              !   - Only theta=90, varphi=0
              !   - symconf must be 0
              ! Switch depending on the dimension of the problem
              select case (problem%n)
                ! 2D
                case (2)
                  if (incidentfield(i)%space.eq.fbem_half_space) then
                    if (incidentfield(i)%np.ne.2) then
                      call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,'np can be only 2.')
                    end if
                    if (incidentfield(i)%bc.ne.1) then
                      call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,'bc must be 1.')
                    end if
                  end if
                  if (incidentfield(i)%wave_type.eq.fbem_harpor_shy_wave) then
                    call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,&
                                            'a shy wave does not make sense.')
                  end if
                ! 3D
                case (3)
                  if (incidentfield(i)%space.eq.fbem_half_space) then
                    if (incidentfield(i)%np.ne.3) then
                      call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,'np can be only 3.')
                    end if
                    if (incidentfield(i)%bc.ne.1) then
                      call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,'bc must be 1.')
                    end if
                  end if
              end select

              do kc=1,problem%n
                if (nint(incidentfield(i)%symconf(kc)).ne.0) then
                  call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,'feature not available')
                end if
              end do

              ! =============================
              ! Check multilayered half-space
              ! =============================

              if (incidentfield(i)%space.eq.fbem_multilayered_half_space) then
                ! Solo incidencia vertical
                if (abs(incidentfield(i)%theta-c_pi_2).gt.1.d-12) then
                  call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,'only vertical incidence is implemented')
                end if

                ! =================================================================
                ! CHANGE EID TO IID IN THE INCIDENT FIELD -> MATERIALS CONNECTIVITY
                ! =================================================================

                do j=1,incidentfield(i)%n_layers
                  if ((incidentfield(i)%layer_material(j).ge.material_eid_min).and.(incidentfield(i)%layer_material(j).le.material_eid_max)) then
                    if (material_iid(incidentfield(i)%layer_material(j)).eq.0) then
                      call fbem_warning_message(error_unit,0,'material',incidentfield(i)%layer_material(j),'does not exist')
                      call fbem_error_message(error_unit,0,'incidentfield',incidentfield(i)%id,'contain the previous material')
                    else
                      incidentfield(i)%layer_material(j)=material_iid(incidentfield(i)%layer_material(j))
                    end if
                  else
                    call fbem_warning_message(error_unit,0,'material',incidentfield(i)%layer_material(j),'does not exist')
                    call fbem_error_message(error_unit,0,'incidentfield',incidentfield(i)%id,'contain the previous material')
                  end if
                end do
                ! Only biot_poroelastic_medium are allowed
                do j=1,incidentfield(i)%n_layers
                  if (trim(material(incidentfield(i)%layer_material(j))%type).ne.'biot_poroelastic_medium') then
                    call fbem_warning_message(error_unit,0,'material',material(incidentfield(i)%layer_material(j))%id,&
                                              'is not a biot_poroelastic_medium.')
                    call fbem_error_message(error_unit,0,'incidentfield',incidentfield(i)%id,'contain the previous material')
                  end if
                end do


                !
                ! Allocate and copy material properties
                !

                allocate (incidentfield(i)%layer_properties(9,incidentfield(i)%n_layers))
                do j=1,incidentfield(i)%n_layers
                  km=incidentfield(i)%layer_material(j)
                  !
                  ! 1 phi
                  ! 2 lambda
                  ! 3 mu
                  ! 4 Q
                  ! 5 R
                  ! 6 rhof
                  ! 7 rhos
                  ! 8 rhoa
                  ! 9 b
                  !
                  ! Copy phi
                  if (material(km)%property_defined(1)) then
                    incidentfield(i)%layer_properties(1,j)=material(km)%property(1,1)
                  else
                    call fbem_warning_message(error_unit,0,'material',material(km)%id,'phi is required')
                    call fbem_error_message(error_unit,0,'incidentfield',incidentfield(i)%id,'contain the previous material')
                  end if
                  ! Copy lambda
                  if (material(km)%property_defined(4)) then
                    if (material(km)%property_defined(12)) then
                      incidentfield(i)%layer_properties(2,j)=material(km)%property(4,1)*(1.d0+c_im*2.d0*material(km)%property(12,1))
                    else
                      incidentfield(i)%layer_properties(2,j)=material(km)%property(4,1)
                    end if
                  else
                    call fbem_warning_message(error_unit,0,'material',material(km)%id,'lambda is required')
                    call fbem_error_message(error_unit,0,'incidentfield',incidentfield(i)%id,'contain the previous material')
                  end if
                  ! Copy mu
                  if (material(km)%property_defined(5)) then
                    if (material(km)%property_defined(12)) then
                      incidentfield(i)%layer_properties(3,j)=material(km)%property(5,1)*(1.d0+c_im*2.d0*material(km)%property(12,1))
                    else
                      incidentfield(i)%layer_properties(3,j)=material(km)%property(5,1)
                    end if
                  else
                    call fbem_warning_message(error_unit,0,'material',material(km)%id,'mu is required')
                    call fbem_error_message(error_unit,0,'incidentfield',incidentfield(i)%id,'contain the previous material')
                  end if
                  ! Copy Q
                  if (material(km)%property_defined(7)) then
                    incidentfield(i)%layer_properties(4,j)=material(km)%property(7,1)
                  else
                    call fbem_warning_message(error_unit,0,'material',material(km)%id,'Q is required')
                    call fbem_error_message(error_unit,0,'incidentfield',incidentfield(i)%id,'contain the previous material')
                  end if
                  ! Copy R
                  if (material(km)%property_defined(8)) then
                    incidentfield(i)%layer_properties(5,j)=material(km)%property(8,1)
                  else
                    call fbem_warning_message(error_unit,0,'material',material(km)%id,'R is required')
                    call fbem_error_message(error_unit,0,'incidentfield',incidentfield(i)%id,'contain the previous material')
                  end if
                  ! Copy rhof
                  if (material(km)%property_defined(9)) then
                    incidentfield(i)%layer_properties(6,j)=material(km)%property(9,1)
                  else
                    call fbem_warning_message(error_unit,0,'material',material(km)%id,'rho_f is required')
                    call fbem_error_message(error_unit,0,'incidentfield',incidentfield(i)%id,'contain the previous material')
                  end if
                  ! Copy rhos
                  if (material(km)%property_defined(10)) then
                    incidentfield(i)%layer_properties(7,j)=material(km)%property(10,1)
                  else
                    call fbem_warning_message(error_unit,0,'material',material(km)%id,'rho_s is required')
                    call fbem_error_message(error_unit,0,'incidentfield',incidentfield(i)%id,'contain the previous material')
                  end if
                  ! Copy rhoa
                  if (material(km)%property_defined(11)) then
                    incidentfield(i)%layer_properties(8,j)=material(km)%property(11,1)
                  else
                    call fbem_warning_message(error_unit,0,'material',material(km)%id,'rho_a is required')
                    call fbem_error_message(error_unit,0,'incidentfield',incidentfield(i)%id,'contain the previous material')
                  end if
                  ! Copy b
                  if (material(km)%property_defined(13)) then
                    incidentfield(i)%layer_properties(9,j)=material(km)%property(13,1)
                  else
                    call fbem_warning_message(error_unit,0,'material',material(km)%id,'rhof is required')
                    call fbem_error_message(error_unit,0,'incidentfield',incidentfield(i)%id,'contain the previous material')
                  end if


                end do

              end if

          end select

      end select

    end do

    ! ==========================================
    ! CHECK AND BUILD INCIDENT FIELD IDENTIFIERS
    ! ==========================================

    incidentfield_eid_min=incidentfield(1)%id
    incidentfield_eid_max=incidentfield(1)%id
    do i=2,n_incidentfields
      if (incidentfield(i)%id.lt.incidentfield_eid_min) incidentfield_eid_min=incidentfield(i)%id
      if (incidentfield(i)%id.gt.incidentfield_eid_max) incidentfield_eid_max=incidentfield(i)%id
    end do
    allocate (incidentfield_iid(incidentfield_eid_min:incidentfield_eid_max))
    incidentfield_iid=0
    do i=1,n_incidentfields
      if (incidentfield_iid(incidentfield(i)%id).ne.0) then
        call fbem_error_message(error_unit,0,'incident wave',incidentfield(i)%id,'this incident wave field is repeated.')
      else
        incidentfield_iid(incidentfield(i)%id)=i
      end if
    end do

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section [incident waves]')

  else

    n_incidentfields=0

  end if

end subroutine read_incident_mechanics_harmonic
