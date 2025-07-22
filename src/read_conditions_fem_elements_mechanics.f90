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

subroutine read_conditions_fem_elements_mechanics(input_fileunit)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_geometry
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_symmetry
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! No implicit variables are allowed
  implicit none

  ! Local variables
  integer                               :: input_fileunit    ! Input file unit
  logical                               :: found             ! Logical variable for sections and keywords
  integer                               :: kn, kg
  integer                               :: ke, se
  integer                               :: kc
  integer                               :: kp, ks
  character(len=fbem_stdcharlen)        :: keyword           ! Auxiliary variable to save a keyword
  character(len=fbem_stdcharlen)        :: tmp_ctype
  character(len=fbem_stdcharlen)        :: section_name

  section_name='conditions over fe elements'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(input_fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    !
    ! ELEMENTS
    !
    ! Loop through the "fe" element
    do ke=1,n_elements
      if (part(element(ke)%part)%type.eq.fbem_part_fe_subregion) then
        ! Find the element
        call fbem_search_section(input_fileunit,'conditions over fe elements',found)
        write(keyword,*) 'element ', element(ke)%id
        call fbem_trim2b(keyword)
        call fbem_search_keyword(input_fileunit,keyword,':',found)
        if (found) then
          !
          ! Esto debera ser refinado y mejorado....... introduciendo incluso el vector director para la carga local....
          ! En principio esto vale bastante bien, pero habria que enriquecerlo para las shells y las vigas degeneradas, para poder
          ! poner cargas en aristas y demas tambien....
          !
          ! Read the type of distributed load
          read(input_fileunit,*) tmp_ctype
          element(ke)%ctype(1,1)=0
          element(ke)%ctype(2,1)=0
          if (trim(tmp_ctype).eq.'local' ) element(ke)%ctype(1,1)=1
          if (trim(tmp_ctype).eq.'global') element(ke)%ctype(2,1)=1
          if (trim(tmp_ctype).eq.'local_and_global') then
            element(ke)%ctype(1,1)=1
            element(ke)%ctype(2,1)=1
          end if
          if ((element(ke)%ctype(1,1).eq.0).and.(element(ke)%ctype(2,1).eq.0)) then
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                    'invalid type of distributed load')
          end if
          !
          ! Depending on the problem type and analysis
          !
          select case (problem%type)
            case (fbem_mechanics)
              select case (problem%analysis)
                case (fbem_static)
                  ! Read the values for local distributed load
                  if (element(ke)%ctype(1,1).eq.1) then
                    do kc=1,problem%n
                      read(input_fileunit,*) (element(ke)%cvalue_r(kc,kn,1),kn=1,element(ke)%n_nodes)
                    end do
                  end if
                  ! Read the values for global distributed load
                  if (element(ke)%ctype(2,1).eq.1) then
                    do kc=1,problem%n
                      read(input_fileunit,*) (element(ke)%cvalue_r(kc,kn,2),kn=1,element(ke)%n_nodes)
                    end do
                  end if
                case (fbem_harmonic)
                  ! Read the values for local distributed load
                  if (element(ke)%ctype(1,1).eq.1) then
                    do kc=1,problem%n
                      read(input_fileunit,*) (element(ke)%cvalue_c(kc,kn,1),kn=1,element(ke)%n_nodes)
                    end do
                  end if
                  ! Read the values for global distributed load
                  if (element(ke)%ctype(2,1).eq.1) then
                    do kc=1,problem%n
                      read(input_fileunit,*) (element(ke)%cvalue_c(kc,kn,2),kn=1,element(ke)%n_nodes)
                    end do
                  end if
              end select
          end select
        end if
      end if
    end do
    !
    ! PART
    !
    ! Loop through each "fe" part
    do kp=1,n_parts
      if (part(kp)%type.eq.fbem_part_fe_subregion) then
        ! Find the part
        call fbem_search_section(input_fileunit,'conditions over fe elements',found)
        write(keyword,*) 'part ', part(kp)%id
        call fbem_trim2b(keyword)
        call fbem_search_keyword(input_fileunit,keyword,':',found)
        if (found) then
          ! Read the type of distributed load
          read(input_fileunit,*) tmp_ctype
          ! First element
          ke=part(kp)%element(1)
          element(ke)%ctype(1,1)=0
          element(ke)%ctype(2,1)=0
          if (trim(tmp_ctype).eq.'local' ) element(ke)%ctype(1,1)=1
          if (trim(tmp_ctype).eq.'global') element(ke)%ctype(2,1)=1
          if (trim(tmp_ctype).eq.'local_and_global') then
            element(ke)%ctype(1,1)=1
            element(ke)%ctype(2,1)=1
          end if
          if ((element(ke)%ctype(1,1).eq.0).and.(element(ke)%ctype(2,1).eq.0)) then
            call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                    'invalid type of distributed load')
          end if
          !
          ! Depending on the problem type and analysis
          !
          select case (problem%type)
            case (fbem_mechanics)
              select case (problem%analysis)
                case (fbem_static)
                  ! Read the values for local distributed load
                  if (element(ke)%ctype(1,1).eq.1) then
                    read(input_fileunit,*) (element(ke)%cvalue_r(kc,1,1),kc=1,problem%n)
                  end if
                  ! Read the values for global distributed load
                  if (element(ke)%ctype(2,1).eq.1) then
                    read(input_fileunit,*) (element(ke)%cvalue_r(kc,1,2),kc=1,problem%n)
                  end if
                  ! Copy to all nodes
                  do kc=1,problem%n
                    element(ke)%cvalue_r(kc,:,1)=element(ke)%cvalue_r(kc,1,1)
                    element(ke)%cvalue_r(kc,:,2)=element(ke)%cvalue_r(kc,1,2)
                  end do
                  ! Copy to the rest of the group
                  do ke=2,part(kp)%n_elements
                    element(part(kp)%element(ke))%ctype(1,1)=element(part(kp)%element(1))%ctype(1,1)
                    element(part(kp)%element(ke))%ctype(2,1)=element(part(kp)%element(1))%ctype(2,1)
                    do kc=1,problem%n
                      element(part(kp)%element(ke))%cvalue_r(kc,:,1)=element(part(kp)%element(1))%cvalue_r(kc,1,1)
                      element(part(kp)%element(ke))%cvalue_r(kc,:,2)=element(part(kp)%element(1))%cvalue_r(kc,1,2)
                    end do
                  end do
                case (fbem_harmonic)
                  ! Read the values for local distributed load
                  if (element(ke)%ctype(1,1).eq.1) then
                    read(input_fileunit,*) (element(ke)%cvalue_c(kc,1,1),kc=1,problem%n)
                  end if
                  ! Read the values for global distributed load
                  if (element(ke)%ctype(2,1).eq.1) then
                    read(input_fileunit,*) (element(ke)%cvalue_c(kc,1,2),kc=1,problem%n)
                  end if
                  ! Copy to all nodes
                  do kc=1,problem%n
                    element(ke)%cvalue_c(kc,:,1)=element(ke)%cvalue_c(kc,1,1)
                    element(ke)%cvalue_c(kc,:,2)=element(ke)%cvalue_c(kc,1,2)
                  end do
                  ! Copy to the rest of the group
                  do ke=2,part(kp)%n_elements
                    element(part(kp)%element(ke))%ctype(1,1)=element(part(kp)%element(1))%ctype(1,1)
                    element(part(kp)%element(ke))%ctype(2,1)=element(part(kp)%element(1))%ctype(2,1)
                    do kc=1,problem%n
                      element(part(kp)%element(ke))%cvalue_c(kc,:,1)=element(part(kp)%element(1))%cvalue_c(kc,1,1)
                      element(part(kp)%element(ke))%cvalue_c(kc,:,2)=element(part(kp)%element(1))%cvalue_c(kc,1,2)
                    end do
                  end do
              end select
          end select
        end if
      end if
    end do
    !
    ! FE SUBREGION
    !
    ! Loop through each FE subregion
    do ks=1,n_fe_subregions
      kp=fe_subregion(ks)%part
      ! Find the fe subregion
      call fbem_search_section(input_fileunit,'conditions over fe elements',found)
      write(keyword,*) 'fe_subregion ', fe_subregion(ks)%id
      call fbem_trim2b(keyword)
      call fbem_search_keyword(input_fileunit,keyword,':',found)

      !
      ! Code similar to "part" finding
      !
      if (found) then
        ! Read the type of distributed load
        read(input_fileunit,*) tmp_ctype
        ! First element
        ke=part(kp)%element(1)
        element(ke)%ctype(1,1)=0
        element(ke)%ctype(2,1)=0
        if (trim(tmp_ctype).eq.'local' ) element(ke)%ctype(1,1)=1
        if (trim(tmp_ctype).eq.'global') element(ke)%ctype(2,1)=1
        if (trim(tmp_ctype).eq.'local_and_global') then
          element(ke)%ctype(1,1)=1
          element(ke)%ctype(2,1)=1
        end if
        if ((element(ke)%ctype(1,1).eq.0).and.(element(ke)%ctype(2,1).eq.0)) then
          call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                  'invalid type of distributed load')
        end if
        !
        ! Depending on the problem type and analysis
        !
        select case (problem%type)
          case (fbem_mechanics)
            select case (problem%analysis)
              case (fbem_static)
                ! Read the values for local distributed load
                if (element(ke)%ctype(1,1).eq.1) then
                  read(input_fileunit,*) (element(ke)%cvalue_r(kc,1,1),kc=1,problem%n)
                end if
                ! Read the values for global distributed load
                if (element(ke)%ctype(2,1).eq.1) then
                  read(input_fileunit,*) (element(ke)%cvalue_r(kc,1,2),kc=1,problem%n)
                end if
                ! Copy to all nodes
                do kc=1,problem%n
                  element(ke)%cvalue_r(kc,:,1)=element(ke)%cvalue_r(kc,1,1)
                  element(ke)%cvalue_r(kc,:,2)=element(ke)%cvalue_r(kc,1,2)
                end do
                ! Copy to the rest of the group
                do ke=2,part(kp)%n_elements
                  element(part(kp)%element(ke))%ctype(1,1)=element(part(kp)%element(1))%ctype(1,1)
                  element(part(kp)%element(ke))%ctype(2,1)=element(part(kp)%element(1))%ctype(2,1)
                  do kc=1,problem%n
                    element(part(kp)%element(ke))%cvalue_r(kc,:,1)=element(part(kp)%element(1))%cvalue_r(kc,1,1)
                    element(part(kp)%element(ke))%cvalue_r(kc,:,2)=element(part(kp)%element(1))%cvalue_r(kc,1,2)
                  end do
                end do
              case (fbem_harmonic)
                ! Read the values for local distributed load
                if (element(ke)%ctype(1,1).eq.1) then
                  read(input_fileunit,*) (element(ke)%cvalue_c(kc,1,1),kc=1,problem%n)
                end if
                ! Read the values for global distributed load
                if (element(ke)%ctype(2,1).eq.1) then
                  read(input_fileunit,*) (element(ke)%cvalue_c(kc,1,2),kc=1,problem%n)
                end if
                ! Copy to all nodes
                do kc=1,problem%n
                  element(ke)%cvalue_c(kc,:,1)=element(ke)%cvalue_c(kc,1,1)
                  element(ke)%cvalue_c(kc,:,2)=element(ke)%cvalue_c(kc,1,2)
                end do
                ! Copy to the rest of the group
                do ke=2,part(kp)%n_elements
                  element(part(kp)%element(ke))%ctype(1,1)=element(part(kp)%element(1))%ctype(1,1)
                  element(part(kp)%element(ke))%ctype(2,1)=element(part(kp)%element(1))%ctype(2,1)
                  do kc=1,problem%n
                    element(part(kp)%element(ke))%cvalue_c(kc,:,1)=element(part(kp)%element(1))%cvalue_c(kc,1,1)
                    element(part(kp)%element(ke))%cvalue_c(kc,:,2)=element(part(kp)%element(1))%cvalue_c(kc,1,2)
                  end do
                end do
            end select
        end select
      end if

    end do

    !
    ! GROUP
    !
    do kg=1,n_groups
      select case (group(kg)%type)
        !
        ! ELEMENTS
        !
        ! Cuidado que el formato de entrada es diferente a cuando se leen por elementos...¿?¿¿?¿?
        !
        !
        case (fbem_group_type_elements)
          call fbem_search_section(input_fileunit,'conditions over fe elements',found)
          write(keyword,*) 'group ', group(kg)%id
          call fbem_trim2b(keyword)
          call fbem_search_keyword(input_fileunit,keyword,':',found)
          if (found) then
            ! Read the type of distributed load
            read(input_fileunit,*) tmp_ctype
            ! First element
            ke=group(kg)%object(1)
            element(ke)%ctype(1,1)=0
            element(ke)%ctype(2,1)=0
            if (trim(tmp_ctype).eq.'local' ) element(ke)%ctype(1,1)=1
            if (trim(tmp_ctype).eq.'global') element(ke)%ctype(2,1)=1
            if (trim(tmp_ctype).eq.'local_and_global') then
              element(ke)%ctype(1,1)=1
              element(ke)%ctype(2,1)=1
            end if
            if ((element(ke)%ctype(1,1).eq.0).and.(element(ke)%ctype(2,1).eq.0)) then
              call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                      'invalid type of distributed load')
            end if
            !
            ! Depending on the problem type and analysis
            !
            select case (problem%type)
              case (fbem_mechanics)
                select case (problem%analysis)
                  case (fbem_static)
                    ! Read the values for local distributed load
                    if (element(ke)%ctype(1,1).eq.1) then
                      read(input_fileunit,*) (element(ke)%cvalue_r(kc,1,1),kc=1,problem%n)
                    end if
                    ! Read the values for global distributed load
                    if (element(ke)%ctype(2,1).eq.1) then
                      read(input_fileunit,*) (element(ke)%cvalue_r(kc,1,2),kc=1,problem%n)
                    end if
                    ! Copy to all nodes
                    do kc=1,problem%n
                      element(ke)%cvalue_r(kc,:,1)=element(ke)%cvalue_r(kc,1,1)
                      element(ke)%cvalue_r(kc,:,2)=element(ke)%cvalue_r(kc,1,2)
                    end do
                    ! Copy to the rest of the group
                    do ke=2,group(kg)%n_objects
                      element(group(kg)%object(ke))%ctype(1,1)=element(group(kg)%object(1))%ctype(1,1)
                      element(group(kg)%object(ke))%ctype(2,1)=element(group(kg)%object(1))%ctype(2,1)
                      do kc=1,problem%n
                        element(group(kg)%object(ke))%cvalue_r(kc,:,1)=element(group(kg)%object(1))%cvalue_r(kc,1,1)
                        element(group(kg)%object(ke))%cvalue_r(kc,:,2)=element(group(kg)%object(1))%cvalue_r(kc,1,2)
                      end do
                    end do
                  case (fbem_harmonic)
                    ! Read the values for local distributed load
                    if (element(ke)%ctype(1,1).eq.1) then
                      read(input_fileunit,*) (element(ke)%cvalue_c(kc,1,1),kc=1,problem%n)
                    end if
                    ! Read the values for global distributed load
                    if (element(ke)%ctype(2,1).eq.1) then
                      read(input_fileunit,*) (element(ke)%cvalue_c(kc,1,2),kc=1,problem%n)
                    end if
                    ! Copy to all nodes
                    do kc=1,problem%n
                      element(ke)%cvalue_c(kc,:,1)=element(ke)%cvalue_c(kc,1,1)
                      element(ke)%cvalue_c(kc,:,2)=element(ke)%cvalue_c(kc,1,2)
                    end do
                    ! Copy to the rest of the group
                    do ke=2,group(kg)%n_objects
                      element(group(kg)%object(ke))%ctype(1,1)=element(group(kg)%object(1))%ctype(1,1)
                      element(group(kg)%object(ke))%ctype(2,1)=element(group(kg)%object(1))%ctype(2,1)
                      do kc=1,problem%n
                        element(group(kg)%object(ke))%cvalue_c(kc,:,1)=element(group(kg)%object(1))%cvalue_c(kc,1,1)
                        element(group(kg)%object(ke))%cvalue_c(kc,:,2)=element(group(kg)%object(1))%cvalue_c(kc,1,2)
                      end do
                    end do
                end select
            end select
          end if
        !
        ! SUBEDGES
        !
        case (fbem_group_type_subedges)
          call fbem_search_section(input_fileunit,'conditions over fe elements',found)
          write(keyword,*) 'group ', group(kg)%id
          call fbem_trim2b(keyword)
          call fbem_search_keyword(input_fileunit,keyword,':',found)
          if (found) then
            ! Read the type of distributed load
            read(input_fileunit,*) tmp_ctype
            ! First element
            ke=group(kg)%object(1)
            subedge(ke)%ctype(1,1)=0
            subedge(ke)%ctype(2,1)=0
            if (trim(tmp_ctype).eq.'local' ) subedge(ke)%ctype(1,1)=1
            if (trim(tmp_ctype).eq.'global') subedge(ke)%ctype(2,1)=1
            if (trim(tmp_ctype).eq.'local_and_global') then
              subedge(ke)%ctype(1,1)=1
              subedge(ke)%ctype(2,1)=1
            end if
            if ((subedge(ke)%ctype(1,1).eq.0).and.(subedge(ke)%ctype(2,1).eq.0)) then
              call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                                      'invalid type of distributed load')
            end if
            !
            ! Depending on the problem type and analysis
            !
            select case (problem%type)
              case (fbem_mechanics)
                select case (problem%analysis)
                  case (fbem_static)
                    ! Read the values for local distributed load
                    if (subedge(ke)%ctype(1,1).eq.1) then
                      read(input_fileunit,*) (subedge(ke)%cvalue_r(kc,1,1),kc=1,problem%n)
                    end if
                    ! Read the values for global distributed load
                    if (subedge(ke)%ctype(2,1).eq.1) then
                      read(input_fileunit,*) (subedge(ke)%cvalue_r(kc,1,2),kc=1,problem%n)
                    end if
                    ! Copy to all nodes
                    do kc=1,problem%n
                      subedge(ke)%cvalue_r(kc,:,1)=subedge(ke)%cvalue_r(kc,1,1)
                      subedge(ke)%cvalue_r(kc,:,2)=subedge(ke)%cvalue_r(kc,1,2)
                    end do
                    ! Copy to the rest of the group
                    do ke=2,group(kg)%n_objects
                      subedge(group(kg)%object(ke))%ctype(1,1)=subedge(group(kg)%object(1))%ctype(1,1)
                      subedge(group(kg)%object(ke))%ctype(2,1)=subedge(group(kg)%object(1))%ctype(2,1)
                      do kc=1,problem%n
                        subedge(group(kg)%object(ke))%cvalue_r(kc,:,1)=subedge(group(kg)%object(1))%cvalue_r(kc,1,1)
                        subedge(group(kg)%object(ke))%cvalue_r(kc,:,2)=subedge(group(kg)%object(1))%cvalue_r(kc,1,2)
                      end do
                    end do
                  case (fbem_harmonic)
                    ! Read the values for local distributed load
                    if (subedge(ke)%ctype(1,1).eq.1) then
                      read(input_fileunit,*) (subedge(ke)%cvalue_c(kc,1,1),kc=1,problem%n)
                    end if
                    ! Read the values for global distributed load
                    if (subedge(ke)%ctype(2,1).eq.1) then
                      read(input_fileunit,*) (subedge(ke)%cvalue_c(kc,1,2),kc=1,problem%n)
                    end if
                    ! Copy to all nodes
                    do kc=1,problem%n
                      subedge(ke)%cvalue_c(kc,:,1)=subedge(ke)%cvalue_c(kc,1,1)
                      subedge(ke)%cvalue_c(kc,:,2)=subedge(ke)%cvalue_c(kc,1,2)
                    end do
                    ! Copy to the rest of the group
                    do ke=2,group(kg)%n_objects
                      subedge(group(kg)%object(ke))%ctype(1,1)=subedge(group(kg)%object(1))%ctype(1,1)
                      subedge(group(kg)%object(ke))%ctype(2,1)=subedge(group(kg)%object(1))%ctype(2,1)
                      do kc=1,problem%n
                        subedge(group(kg)%object(ke))%cvalue_c(kc,:,1)=subedge(group(kg)%object(1))%cvalue_c(kc,1,1)
                        subedge(group(kg)%object(ke))%cvalue_c(kc,:,2)=subedge(group(kg)%object(1))%cvalue_c(kc,1,2)
                      end do
                    end do
                end select
            end select
          end if

      end select

    end do

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  end if

end subroutine read_conditions_fem_elements_mechanics
