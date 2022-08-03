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
!! <b> Subroutine that reads the parts from a file. </b>
subroutine read_parts(fileunit,mode)

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
  integer                        :: fileunit          !! Unit of the file to read from
  integer                        :: mode              !! Mode of reading: 1 (native format), 2 (gmsh format)
  ! Local
  character(len=fbem_stdcharlen) :: section_name      ! Name of the section
  logical                        :: found             ! Logical variable for sections and keywords
  integer                        :: i, j              ! Counters
  integer                        :: tmp_int           ! Temporary integer

  ! Locate the section
  select case (mode)
    case (1)
      section_name='parts'
      call fbem_search_section(fileunit,section_name,found)
    case (2)
      section_name='PhysicalNames'
      call fbem_search_section_gmsh(fileunit,section_name,found)
  end select
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! Read the number of parts
    read(fileunit,*) n_parts

    ! Allocate and initialize
    if (n_parts.gt.0) then
      allocate (part(n_parts))
      do i=1,n_parts
        part(i)%id=0
        part(i)%name=''
        part(i)%type=0
        part(i)%entity=0
        part(i)%n_nodes=0
        part(i)%n_elements=0
        part(i)%export=.false.
      end do
    else
      call fbem_error_message(error_unit,0,trim(section_name),0,'the number of parts must be >0')
    end if

    ! Read each part line
    do i=1,n_parts
      ! Read part id and name
      select case (mode)
        case (1)
          read(fileunit,*) part(i)%id, part(i)%name
        case (2)
          read(fileunit,*) tmp_int, part(i)%id, part(i)%name
      end select
      ! Check id
      if (part(i)%id.le.0) then
        call fbem_error_message(error_unit,0,'part',part(i)%id,'identifiers must be greater than 0')
      end if
      ! Remove starting and ending blanks
      call fbem_trim(part(i)%name)
    end do

    ! =================================
    ! CHECK AND BUILD PARTS IDENTIFIERS
    ! =================================

    part_eid_min=part(1)%id
    part_eid_max=part(1)%id
    do i=2,n_parts
      if (part(i)%id.lt.part_eid_min) part_eid_min=part(i)%id
      if (part(i)%id.gt.part_eid_max) part_eid_max=part(i)%id
    end do
    allocate (part_iid(part_eid_min:part_eid_max))
    part_iid=0
    do i=1,n_parts
      if (part_iid(part(i)%id).ne.0) then
        call fbem_error_message(error_unit,0,'part',part(i)%id,'is repeated')
      else
        part_iid(part(i)%id)=i
      end if
    end do

    ! =========================================================
    ! CHANGE EID TO IID IN THE BE BOUNDARY -> PART CONNECTIVITY
    ! =========================================================

    do i=1,n_boundaries
      if ((boundary(i)%part.ge.part_eid_min).and.(boundary(i)%part.le.part_eid_max)) then
        if (part_iid(boundary(i)%part).eq.0) then
          call fbem_warning_message(error_unit,0,'part',boundary(i)%part,'does not exist')
          call fbem_error_message(error_unit,0,'BE boundary',boundary(i)%id,'is connected to the non-existent previously indicated part')
        else
          boundary(i)%part=part_iid(boundary(i)%part)
        end if
      else
        call fbem_warning_message(error_unit,0,'part',boundary(i)%part,'does not exist')
        call fbem_error_message(error_unit,0,'BE boundary',boundary(i)%id,'is connected to the non-existent previously indicated part')
      end if
    end do

    ! ==========================================================
    ! CHANGE EID TO IID IN THE FE SUBREGION -> PART CONNECTIVITY
    ! ==========================================================

    do i=1,n_fe_subregions
      if ((fe_subregion(i)%part.ge.part_eid_min).and.(fe_subregion(i)%part.le.part_eid_max)) then
        if (part_iid(fe_subregion(i)%part).eq.0) then
          call fbem_warning_message(error_unit,0,'part',fe_subregion(i)%part,'does not exist')
          call fbem_error_message(error_unit,0,'FE subregion',fe_subregion(i)%id,'is connected to the non-existent previously indicated part')
        else
          fe_subregion(i)%part=part_iid(fe_subregion(i)%part)
        end if
      else
        call fbem_warning_message(error_unit,0,'part',fe_subregion(i)%part,'does not exist')
        call fbem_error_message(error_unit,0,'FE subregion',fe_subregion(i)%id,'is connected to the non-existent previously indicated part')
      end if
    end do

    ! ==========================================================
    ! CHANGE EID TO IID IN THE BE BODY LOAD -> PART CONNECTIVITY
    ! ==========================================================

    do i=1,n_be_bodyloads
      if ((be_bodyload(i)%part.ge.part_eid_min).and.(be_bodyload(i)%part.le.part_eid_max)) then
        if (part_iid(be_bodyload(i)%part).eq.0) then
          call fbem_warning_message(error_unit,0,'part',be_bodyload(i)%part,'does not exist')
          call fbem_error_message(error_unit,0,'BE body load',be_bodyload(i)%id,'is connected to the non-existent previously indicated part')
        else
          be_bodyload(i)%part=part_iid(be_bodyload(i)%part)
        end if
      else
        call fbem_warning_message(error_unit,0,'part',be_bodyload(i)%part,'does not exist')
        call fbem_error_message(error_unit,0,'BE body load',be_bodyload(i)%id,'is connected to the non-existent previously indicated part')
      end if
    end do

    ! ======================================================================
    ! BUILD PART -> BE BOUNDARY OR FE SUBREGION OR BE BODY LOAD CONNECTIVITY
    ! ======================================================================

    do i=1,n_parts

      ! -----------
      ! BE BOUNDARY
      ! -----------

      do j=1,n_boundaries
        if (boundary(j)%part.eq.i) then
          if (part(i)%entity.eq.0) then
            part(i)%type=fbem_part_be_boundary
            part(i)%entity=j
          else
            call fbem_warning_message(error_unit,0,'BE boundary',boundary(j)%id,'is connected to the following part')
            select case (part(i)%type)
              case (fbem_part_be_boundary)
                call fbem_warning_message(error_unit,0,'BE boundary',boundary(part(i)%entity)%id,'is connected to the following part')
              case (fbem_part_fe_subregion)
                call fbem_warning_message(error_unit,0,'FE subregion',fe_subregion(part(i)%entity)%id,'is connected to the following part')
              case (fbem_part_be_bodyload)
                call fbem_warning_message(error_unit,0,'BE body load',be_bodyload(part(i)%entity)%id,'is connected to the following part')
            end select
            call fbem_error_message(error_unit,0,'part',part(i)%id,'connected to the two previously indicated physical entities')
          end if
        end if
      end do

      ! ------------
      ! FE SUBREGION
      ! ------------

      do j=1,n_fe_subregions
        if (fe_subregion(j)%part.eq.i) then
          if (part(i)%entity.eq.0) then
            part(i)%type=fbem_part_fe_subregion
            part(i)%entity=j
          else
            call fbem_warning_message(error_unit,0,'FE subregion',fe_subregion(j)%id,'is connected to the following part')
            select case (part(i)%type)
              case (fbem_part_be_boundary)
                call fbem_warning_message(error_unit,0,'BE boundary',boundary(part(i)%entity)%id,'is connected to the following part')
              case (fbem_part_fe_subregion)
                call fbem_warning_message(error_unit,0,'FE subregion',fe_subregion(part(i)%entity)%id,'is connected to the following part')
              case (fbem_part_be_bodyload)
                call fbem_warning_message(error_unit,0,'BE body load',be_bodyload(part(i)%entity)%id,'is connected to the following part')
            end select
            call fbem_error_message(error_unit,0,'part',part(i)%id,'connected to the two previously indicated physical entities')
          end if
        end if
      end do

      ! ------------
      ! BE BODY LOAD
      ! ------------

      do j=1,n_be_bodyloads
        if (be_bodyload(j)%part.eq.i) then
          if (part(i)%entity.eq.0) then
            part(i)%type=fbem_part_be_bodyload
            part(i)%entity=j
          else
            call fbem_warning_message(error_unit,0,'BE body load',be_bodyload(j)%id,'is connected to the following part')
            select case (part(i)%type)
              case (fbem_part_be_boundary)
                call fbem_warning_message(error_unit,0,'BE boundary',boundary(part(i)%entity)%id,'is connected to the following part')
              case (fbem_part_fe_subregion)
                call fbem_warning_message(error_unit,0,'FE subregion',fe_subregion(part(i)%entity)%id,'is connected to the following part')
              case (fbem_part_be_bodyload)
                call fbem_warning_message(error_unit,0,'BE body load',be_bodyload(part(i)%entity)%id,'is connected to the following part')
            end select
            call fbem_error_message(error_unit,0,'part',part(i)%id,'connected to the two previously indicated physical entities')
          end if
        end if
      end do

      if (part(i)%type.eq.0) then
        call fbem_warning_message(error_unit,0,'part',part(i)%id,'is not connected to any physical entity.')
      end if

    end do

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  else

    call fbem_error_message(error_unit,0,trim(section_name),0,'this section is required')

  end if

end subroutine read_parts
