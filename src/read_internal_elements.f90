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
!! <b> Read and create internal elements (resultants calculations) </b>

subroutine read_internal_elements(input_fileunit)
  ! Fortran 2003 intrinsic module
  use iso_fortran_env
  ! fbem modules
  use fbem_numerical
  use fbem_geometry
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_symmetry
  use fbem_data_structures
  use fbem_quad_rules
  ! Problem variables module
  use problem_variables
  ! No implicit variables are allowed
  implicit none
  ! I/O
  integer                                 :: input_fileunit    ! Input file unit
  ! Local
  character(len=fbem_stdcharlen)          :: section_name      ! Name of the section
  logical                                 :: found             ! Logical variable for sections and keywords
  integer                                 :: i, j, j1, j2, k, sn, ke, se
  integer                                 :: sp, sr, kp, kr
  integer                                 :: n_additional_internalpoints
  type(fbem_internalpoint), allocatable   :: internalpoint_tmp(:)
  integer, allocatable                    :: internalpoint_iid_tmp(:)
  integer                                 :: iid, eid
  integer                                 :: etype
  real(kind=real64), allocatable          :: xi(:), phi(:)
  real(kind=real64)                       :: centroid_part(problem%n), area_part, delta_f, xi1d, xi2d(2), xi3d(3)
  real(kind=real64)                       :: aux(10)
  character(len=fbem_stdcharlen)          :: fileformat

  ! Locate the section
  section_name='internal elements'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(input_fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    internalelements=.true.

    ! ---------------------------
    ! READ INTERNAL ELEMENTS FILE
    ! ---------------------------

    ! From specific mesh file
    call fbem_search_section(input_fileunit,section_name,found)
    if (found) call fbem_search_keyword(input_fileunit,'mesh_file','=',found)
    if (found) then
      read(input_fileunit,*) internalelements_fileformat, internalelements_filename
      call fbem_trimall(internalelements_filename)
      ! If path to the same directory, remove it.
      if (internalelements_filename(1:2).eq.'./') then
        internalelements_filename=trim(internalelements_filename(3:len_trim(internalelements_filename)))
      end if
      ! Not yet allowed the full syntax to navigate between folders
      if (internalelements_filename(1:2).eq.'..') then
        call fbem_error_message(error_unit,0,'mesh_file',0,'wrong path to the file')
      end if
      ! Add path from the input file
      if (internalelements_filename(1:1).ne.'/') then
        internalelements_filename=trim(pwd)//trim(internalelements_filename)
      end if
      call fbem_trimall(internalelements_filename)
      if (verbose_level.ge.3) then
        write(output_unit,'(2x,a,i1)') 'internalelements_fileformat = ', internalelements_fileformat
        write(output_unit,'(2x,a,a)') 'internalelements_filename = ', trim(internalelements_filename)
      end if
    ! From main mesh file
    else
      internalelements_fileformat=mesh_file_mode
      if (internalelements_fileformat.eq.0) then
        internalelements_filename=input_filename
      else
        internalelements_filename=mesh_filename
      end if
    end if

    ! Read the mesh from the file
    select case (internalelements_fileformat)
      case (0,1)
        fileformat='multifebe'
      case (2)
        fileformat='gmsh'
      case default
        call fbem_error_message(error_unit,0,'mesh_file',0,'invalid file format')
    end select
    call internalelements_mesh%read_from_file(problem%n,geometric_tolerance,internalelements_filename,fileformat)
    ! Characterize and check the elements
    do j=1,internalelements_mesh%n_elements
      ! Check the dimension
      !if (internalelements_mesh%element(j)%n_dimension.ne.(problem%n-1)) then
      !  call fbem_error_message(error_unit,0,'mesh_file',internalelements_mesh%element(j)%id,&
      !                          'this element does not have the correct dimension')
      !end if
      ! Copy the coordinates of the nodes of each element to the element x_gn
      allocate (internalelements_mesh%element(j)%x_gn(problem%n,internalelements_mesh%element(j)%n_nodes))
      etype=internalelements_mesh%element(j)%type
      do k=1,internalelements_mesh%element(j)%n_nodes
        sn=internalelements_mesh%element(j)%node(k)
        internalelements_mesh%element(j)%x_gn(:,k)=internalelements_mesh%node(sn)%x
      end do
      ! Calculate the characteristic length
      internalelements_mesh%element(j)%csize=fbem_characteristic_length(problem%n,etype,internalelements_mesh%element(j)%x_gn,1.d-9)
      ! Calculate the size
      internalelements_mesh%element(j)%size=fbem_element_size(problem%n,etype,internalelements_mesh%element(j)%x_gn,1.d-9)
      ! Calculate the centroid of the element
      allocate (internalelements_mesh%element(j)%centroid(problem%n))
      internalelements_mesh%element(j)%centroid=fbem_element_centroid(problem%n,etype,internalelements_mesh%element(j)%x_gn,&
                                                                      internalelements_mesh%element(j)%size,1.d-9)
    end do

    ! ----------------------------------------------------------------
    ! READ THE CONNECTIVITY BETWEEN INTERNAL ELEMENTS PARTS TO REGIONS
    ! ----------------------------------------------------------------
    !
    ! region_association : <n_parts>
    ! <internal element mesh part id.> <region id.>
    ! ...
    call fbem_search_section(input_fileunit,section_name,found)
    if (found) call fbem_search_keyword(input_fileunit,'region_association',':',found)
    if (found) then
      ! For each part of the internal elements mesh, read the region where they are.
      read(input_fileunit,*) j
      do j=1,j
        read(input_fileunit,*) sp, sr
        ! Check that sp is a valid part in the internal elements mesh, and then obtain its iid and save it to kp
        if ((sp.ge.internalelements_mesh%part_eid_min).and.(sp.le.internalelements_mesh%part_eid_max)) then
          if (internalelements_mesh%part_iid(sp).eq.0) then
            call fbem_error_message(error_unit,0,section_name,sp,'the indicated part in the internal elements mesh does not exist')
          else
            kp=internalelements_mesh%part_iid(sp)
          end if
        else
          call fbem_error_message(error_unit,0,section_name,sp,'the indicated part in the internal elements mesh does not exist')
        end if
        ! Check that sr is a valid region, obtain its iid, and save it to the internal element part member "entity"
        if ((sr.ge.region_eid_min).and.(sr.le.region_eid_max)) then
          if (region_iid(sr).eq.0) then
            call fbem_error_message(error_unit,0,section_name,sr,'the indicated region does not exist')
          else
            kr=region_iid(sr)
            internalelements_mesh%part(kp)%entity=kr
          end if
        else
          call fbem_error_message(error_unit,0,section_name,sr,'the indicated region does not exist')
        end if
      end do
    else
      call fbem_error_message(error_unit,0,section_name,0,'the list "region_association" must be defined')
    end if

    ! ----------------------------------------------------
    ! READ THE STRESS APPROXIMATION ORDER OVER THE ELEMENT
    ! ----------------------------------------------------

    call fbem_search_section(input_fileunit,section_name,found)
    if (found) call fbem_search_keyword(input_fileunit,'order','=',found)
    if (found) then
      read(input_fileunit,*) internalelements_order
      if ((internalelements_order.lt.1).or.(internalelements_order.gt.30)) then
        call fbem_error_message(error_unit,0,section_name,internalelements_order,'the variable "order" must be between 1 and 30')
      end if
    else
      internalelements_order=2
    end if

    ! -----------------------------------
    ! GENERATE ADDITIONAL INTERNAL POINTS
    ! -----------------------------------

    if (n_internalpoints.gt.0) then
      !
      ! Copy the internal points to temporary variables
      !
      allocate (internalpoint_tmp(n_internalpoints))
      allocate (internalpoint_iid_tmp(internalpoint_eid_min:internalpoint_eid_max))
      do i=1,n_internalpoints
        allocate (internalpoint_tmp(i)%x(problem%n))
        internalpoint_tmp(i)%id=internalpoint(i)%id
        internalpoint_tmp(i)%x=internalpoint(i)%x
        internalpoint_tmp(i)%region=internalpoint(i)%region
        internalpoint_tmp(i)%export=internalpoint(i)%export
      end do
      internalpoint_iid_tmp=internalpoint_iid
      !
      ! Deallocate the internal points
      !
      do i=1,n_internalpoints
        deallocate (internalpoint(i)%x)
      end do
      deallocate (internalpoint)
      deallocate (internalpoint_iid)
    end if
    !
    ! Count the number of additional internal points necessary
    !
    n_additional_internalpoints=0
    do kp=1,internalelements_mesh%n_parts
      if (internalelements_mesh%part(kp)%entity.eq.0) cycle
      if (region(internalelements_mesh%part(kp)%entity)%class.eq.fbem_be) then
        do ke=1,internalelements_mesh%part(kp)%n_elements
          se=internalelements_mesh%part(kp)%element(ke)
          !
          ! Internal points at nodes of discontinuous elements
          !
          n_additional_internalpoints=n_additional_internalpoints+internalelements_mesh%element(se)%n_nodes
!            select case (internalelements_mesh%element(se)%n_dimension)
!              case (1)
!                n_additional_internalpoints=n_additional_internalpoints+
!                if (problem%n.ne.2) then
!                  call fbem_error_message(error_unit,0,section_name,internalelements_mesh%element(se)%id,&
!                                        'line internal elements are allowed only in 2D')
!                end if
!              case (2)
!                select case (fbem_n_edges(internalelements_mesh%element(se)%type))
!                  case (3)
!                    n_additional_internalpoints=n_additional_internalpoints+wantri_n(internalelements_order)
!                  case (4)
!                    n_additional_internalpoints=n_additional_internalpoints+gl11_n(internalelements_order)**2
!                end select
!                if (problem%n.ne.3) then
!                  call fbem_error_message(error_unit,0,section_name,internalelements_mesh%element(se)%id,&
!                                        'surface internal elements are allowed only in 3D')
!                end if
!              case default
!                call fbem_error_message(error_unit,0,section_name,internalelements_mesh%element(se)%id,&
!                                        'the indicated element has an invalid dimension')
!            end select
        end do
      end if
    end do
    !
    ! Allocate the internal points list including the additional internal points, and save the internal points in it.
    !
    allocate (internalpoint(n_internalpoints+n_additional_internalpoints))
    do i=1,n_internalpoints
      allocate (internalpoint(i)%x(problem%n))
      internalpoint(i)%id=internalpoint_tmp(i)%id
      internalpoint(i)%x=internalpoint_tmp(i)%x
      internalpoint(i)%region=internalpoint_tmp(i)%region
      internalpoint(i)%export=internalpoint_tmp(i)%export
    end do
    !
    ! Initialize external and internal identifiers
    !
    if (n_internalpoints.gt.0) then
      allocate (internalpoint_iid(internalpoint_eid_min:internalpoint_eid_max+n_additional_internalpoints))
      internalpoint_iid(internalpoint_eid_min:internalpoint_eid_max)=internalpoint_iid_tmp
      eid=internalpoint_eid_max+1
      iid=n_internalpoints+1
    else
      internalpoint_eid_min=1
      internalpoint_eid_max=n_additional_internalpoints
      allocate (internalpoint_iid(internalpoint_eid_min:internalpoint_eid_max))
      eid=1
      iid=1
    end if
    !
    ! Generate additional internal points
    !
    do kp=1,internalelements_mesh%n_parts
      if (internalelements_mesh%part(kp)%entity.eq.0) cycle
      if (region(internalelements_mesh%part(kp)%entity)%class.eq.fbem_be) then
        do ke=1,internalelements_mesh%part(kp)%n_elements
          se=internalelements_mesh%part(kp)%element(ke)
          allocate (xi(internalelements_mesh%element(se)%n_dimension))
          etype=internalelements_mesh%element(se)%type
          !
          ! MODE: Internal points at nodes of discontinuous elements (internal element order neglected)
          !
          internalelements_mesh%element(se)%type_g=etype
          internalelements_mesh%element(se)%type_f1=etype
          internalelements_mesh%element(se)%type_f2=etype
          internalelements_mesh%element(se)%discontinuous=.true.
          internalelements_mesh%element(se)%n_internalpoints=internalelements_mesh%element(se)%n_nodes
          allocate (internalelements_mesh%element(se)%internalpoint(internalelements_mesh%element(se)%n_nodes))
          allocate (phi(internalelements_mesh%element(se)%n_nodes))
          ! Used delta parameter for discontinuous elements
          select case (etype)
            ! Linear elements
            case (fbem_line2,fbem_tri3,fbem_quad4,fbem_tet4,fbem_hex8)
              delta_f=0.42264973d0
            ! Quadratic elements
            case (fbem_line3,fbem_tri6,fbem_quad8,fbem_quad9,fbem_tet10,fbem_hex20,fbem_hex27)
              delta_f=0.22540333d0
            ! Cubic elements
            case (fbem_line4)
              delta_f=0.138863688d0
          end select
          internalelements_mesh%element(se)%delta_f=delta_f
          do j=1,internalelements_mesh%element(se)%n_nodes
            internalpoint(iid)%id=eid
            internalpoint_iid(eid)=iid
            allocate (internalpoint(iid)%x(problem%n),internalpoint(iid)%n(problem%n))
            select case (internalelements_mesh%element(se)%n_dimension)
              case (1)
#               define delta delta_f
#               define node j
#               define xi xi1d
#               include <xi_1d_at_node.rc>
#               undef node
#               undef delta
#               define delta 0.d0
#               include <phi_1d.rc>
#               undef delta
#               undef xi
                internalpoint(iid)%x=0.d0
                do i=1,internalelements_mesh%element(se)%n_nodes
                  internalpoint(iid)%x=internalpoint(iid)%x+phi(i)*internalelements_mesh%element(se)%x_gn(:,i)
                end do
                internalpoint(iid)%n=0.
                if (problem%n.eq.2) internalpoint(iid)%n=fbem_unormal2d(etype,internalelements_mesh%element(se)%x_gn,xi(1))
              case (2)
#               define delta delta_f
#               define node j
#               include <xi_2d_at_node.rc>
#               undef node
#               undef delta
#               define delta 0.d0
#               include <phi_2d.rc>
#               undef delta
                internalpoint(iid)%x=0.d0
                do i=1,internalelements_mesh%element(se)%n_nodes
                  internalpoint(iid)%x=internalpoint(iid)%x+phi(i)*internalelements_mesh%element(se)%x_gn(:,i)
                end do
                internalpoint(iid)%n=0.
                if (problem%n.eq.3) internalpoint(iid)%n=fbem_unormal3d(etype,internalelements_mesh%element(se)%x_gn,xi)
              case (3)
#               define delta delta_f
#               define node j
#               include <xi_3d_at_node.rc>
#               undef node
#               undef delta
#               define delta 0.d0
#               include <phi_3d.rc>
#               undef delta
                internalpoint(iid)%x=0.d0
                do i=1,internalelements_mesh%element(se)%n_nodes
                  internalpoint(iid)%x=internalpoint(iid)%x+phi(i)*internalelements_mesh%element(se)%x_gn(:,i)
                end do
                internalpoint(iid)%n=0.
            end select
            internalpoint(iid)%region=internalelements_mesh%part(internalelements_mesh%element(se)%part)%entity
            internalpoint(iid)%export=.true.
            internalelements_mesh%element(se)%internalpoint(j)=iid
            eid=eid+1
            iid=iid+1
          end do

!            !
!            ! Internal points directly at Gauss points
!            !
!            select case (internalelements_mesh%element(se)%n_dimension)
!              case (1)
!                internalelements_mesh%element(se)%n_internalpoints=gl11_n(internalelements_order)
!                allocate (internalelements_mesh%element(se)%internalpoint(gl11_n(internalelements_order)))
!                do j=1,gl11_n(internalelements_order)
!                  internalpoint(iid)%id=eid
!                  internalpoint_iid(eid)=iid
!                  allocate (internalpoint(iid)%x(problem%n),internalpoint(iid)%n(problem%n))
!                  xi=gl11_xi(j,internalelements_order)
!                  internalpoint(iid)%x=fbem_position(problem%n,etype,internalelements_mesh%element(se)%x_gn,xi)
!                  internalpoint(iid)%n=fbem_unormal2d(etype,internalelements_mesh%element(se)%x_gn,xi(1))
!                  internalpoint(iid)%region=internalelements_mesh%part(internalelements_mesh%element(se)%part)%entity
!                  internalpoint(iid)%export=.true.
!                  internalelements_mesh%element(se)%internalpoint(j)=iid
!                  eid=eid+1
!                  iid=iid+1
!                end do
!              case (2)
!                select case (fbem_n_edges(internalelements_mesh%element(se)%type))
!                  case (3)
!                    internalelements_mesh%element(se)%n_internalpoints=wantri_n(internalelements_order)
!                    allocate (internalelements_mesh%element(se)%internalpoint(wantri_n(internalelements_order)))
!                    do j=1,wantri_n(internalelements_order)
!                      internalpoint(iid)%id=eid
!                      internalpoint_iid(eid)=iid
!                      allocate (internalpoint(iid)%x(problem%n),internalpoint(iid)%n(problem%n))
!                      xi(1)=wantri_xi1(j,internalelements_order)
!                      xi(2)=wantri_xi2(j,internalelements_order)
!                      internalpoint(iid)%x=fbem_position(problem%n,etype,internalelements_mesh%element(se)%x_gn,xi)
!                      internalpoint(iid)%n=fbem_unormal3d(etype,internalelements_mesh%element(se)%x_gn,xi)
!                      internalpoint(iid)%region=internalelements_mesh%part(internalelements_mesh%element(se)%part)%entity
!                      internalpoint(iid)%export=.true.
!                      internalelements_mesh%element(se)%internalpoint(j)=iid
!                      eid=eid+1
!                      iid=iid+1
!                    end do
!                  case (4)
!                    internalelements_mesh%element(se)%n_internalpoints=gl11_n(internalelements_order)**2
!                    allocate (internalelements_mesh%element(se)%internalpoint(gl11_n(internalelements_order)**2))
!                    do j1=1,gl11_n(internalelements_order)
!                      do j2=1,gl11_n(internalelements_order)
!                        internalpoint(iid)%id=eid
!                        internalpoint_iid(eid)=iid
!                        allocate (internalpoint(iid)%x(problem%n),internalpoint(iid)%n(problem%n))
!                        xi(1)=gl11_xi(j1,internalelements_order)
!                        xi(2)=gl11_xi(j2,internalelements_order)
!                        internalpoint(iid)%x=fbem_position(problem%n,etype,internalelements_mesh%element(se)%x_gn,xi)
!                        internalpoint(iid)%n=fbem_unormal3d(etype,internalelements_mesh%element(se)%x_gn,xi)
!                        internalpoint(iid)%region=internalelements_mesh%part(internalelements_mesh%element(se)%part)%entity
!                        internalpoint(iid)%export=.true.
!                        k=j2+(j1-1)*gl11_n(internalelements_order)
!                        internalelements_mesh%element(se)%internalpoint(k)=iid
!                        eid=eid+1
!                        iid=iid+1
!                      end do
!                    end do
!                end select
!              case (3)
!                stop 'volume internal post-processing elements not yet allowed'
!            end select
          deallocate (phi,xi)
        end do
      end if
    end do
    internalpoint_eid_max=internalpoint_eid_max+n_additional_internalpoints
    n_internalpoints=n_internalpoints+n_additional_internalpoints

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  else
    internalelements=.false.
  end if

end subroutine read_internal_elements
