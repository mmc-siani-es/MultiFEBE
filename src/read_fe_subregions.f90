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
!! <b> Subroutine that reads the FE subregions from a file. </b>
subroutine read_fe_subregions(fileunit)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_data_structures

  ! Problem variables module
  use problem_variables

  ! Local variables
  implicit none
  ! I/O
  integer                        :: fileunit          ! Unit of the file to read from
  ! Local
  character(len=fbem_stdcharlen) :: section_name      ! Name of the section
  logical                        :: found
  integer                        :: i, j, k           ! Counters
  integer                        :: kr, ksr, ssr

  ! Return if not needed
  if (n_fe_regions.eq.0) then
    n_fe_subregions=0
    return
  end if

  ! Locate the section
  section_name='fe subregions'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! Read the number of FE regions
    read(fileunit,*) n_fe_subregions

    ! Allocate and initialize
    if (n_fe_subregions.gt.0) then
      allocate (fe_subregion(n_fe_subregions))
      do i=1,n_fe_subregions
        fe_subregion(i)%id=0
        fe_subregion(i)%name=''
        fe_subregion(i)%class=0
        fe_subregion(i)%type=0
        fe_subregion(i)%subtype=0
        fe_subregion(i)%model=0
        fe_subregion(i)%n_materials=1
        allocate (fe_subregion(i)%material(1))
        fe_subregion(i)%part=0
        fe_subregion(i)%region=0

        !
        ! Experimental...
        !
        fe_subregion(i)%master_node=0


        fe_subregion(i)%sensitivity=.false.
        fe_subregion(i)%export=.false.
      end do
    else
      call fbem_error_message(error_unit,0,section_name,0,'the number of FE subregions must be >0')
    end if

    ! Read FE subregion id, its part, material and model.
    do i=1,n_fe_subregions
      read(fileunit,*) fe_subregion(i)%id, fe_subregion(i)%part, fe_subregion(i)%material(1), fe_subregion(i)%model
      ! Check id
      if (fe_subregion(i)%id.le.0) then
        call fbem_error_message(error_unit,0,'FE subregion',fe_subregion(i)%id,'identifiers must be greater than 0')
      end if

      ! Check material and model (that should be defined previously)

    end do

    ! =========================================
    ! CHECK AND BUILD FE SUBREGIONS IDENTIFIERS
    ! =========================================

    fe_subregion_eid_min=fe_subregion(1)%id
    fe_subregion_eid_max=fe_subregion(1)%id
    do i=2,n_fe_subregions
      if (fe_subregion(i)%id.lt.fe_subregion_eid_min) fe_subregion_eid_min=fe_subregion(i)%id
      if (fe_subregion(i)%id.gt.fe_subregion_eid_max) fe_subregion_eid_max=fe_subregion(i)%id
    end do
    allocate (fe_subregion_iid(fe_subregion_eid_min:fe_subregion_eid_max))
    fe_subregion_iid=0
    do i=1,n_fe_subregions
      if (fe_subregion_iid(fe_subregion(i)%id).ne.0) then
        call fbem_error_message(error_unit,0,'FE subregion',fe_subregion(i)%id,'is repeated')
      else
        fe_subregion_iid(fe_subregion(i)%id)=i
      end if
    end do

    ! ===================================================================
    ! CHECK AND BUILD CONNECTIVITIES BETWEEN FE SUBREGIONS AND FE REGIONS
    ! ===================================================================

    ! Convert eid to iid in FE region -> FE subregion connectivity, and check if the indicated FE subregion actually exist.
    do kr=1,n_regions
      if (region(kr)%class.eq.fbem_fe) then
        do ksr=1,region(kr)%n_fe_subregions
          ssr=region(kr)%fe_subregion(ksr)
          if ((ssr.ge.fe_subregion_eid_min).and.(ssr.le.fe_subregion_eid_max)) then
            if (fe_subregion_iid(ssr).eq.0) then
              call fbem_warning_message(error_unit,0,'FE region',region(kr)%id,'wrong definition of its subregions')
              call fbem_error_message(error_unit,0,'FE subregion',ssr,'does not exist')
            else
              region(kr)%fe_subregion(ksr)=fe_subregion_iid(ssr)
            end if
          else
            call fbem_warning_message(error_unit,0,'FE region',region(kr)%id,'wrong definition of its subregions')
            call fbem_error_message(error_unit,0,'FE subregion',ssr,'does not exist')
          end if
        end do
      end if
    end do

    ! Build FE subregion -> FE region connectivities.
    do i=1,n_fe_subregions
      do j=1,n_regions
        if (region(j)%class.eq.fbem_fe) then
          do k=1,region(j)%n_fe_subregions
            if (region(j)%fe_subregion(k).eq.i) then
              if (fe_subregion(i)%region.eq.0) then
                fe_subregion(i)%region=j
              else
                call fbem_warning_message(error_unit,0,'FE region',region(fe_subregion(i)%region)%id,&
                                          'is connected with the following FE subregion')
                call fbem_warning_message(error_unit,0,'FE region',region(j)%id,&
                                          'is connected with the following FE subregion')
                call fbem_error_message(error_unit,0,'FE subregion',fe_subregion(i)%id,&
                                        'is connected to more than one FE region')
              end if
            end if
          end do
        end if
      end do
    end do

    ! Ending message
    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  else

    call fbem_error_message(error_unit,0,'['//trim(section_name)//']',0,'this section is required')

  end if

end subroutine read_fe_subregions
