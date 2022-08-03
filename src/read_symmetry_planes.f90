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
!! <b> Subroutine that reads the symmetry planes configuration from a file. </b>
!!
!!>
!!
!! This section has only 1 variable:
!!
!!   ----------------------- -------- -------- --------------------------------- ------------------------ -------------------------
!!        User variable name Required     Type                       User values            Default value     Program variable name
!!   ----------------------- -------- -------- --------------------------------- ------------------------ -------------------------
!!                         x       no     ints                    s1 t11 t12 t13                        -        symplane variables
!!                         y       no     ints                    s2 t21 t22 t23                        -        symplane variables
!!                         z       no     ints                    s3 t31 t32 t33                        -        symplane variables
!!   ----------------------- -------- -------- --------------------------------- ------------------------ -------------------------
!!
!!<
subroutine read_symmetry_planes(fileunit)

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
  character(len=fbem_fmtstr)     :: fmtstr            ! String used for write format string
  logical                        :: found             ! Logical variable for sections and keywords
  integer                        :: i, j              ! Counters
  integer                        :: tmp_s, tmp_t(3)   ! Auxiliary variables for symmetry planes reading
  character(len=fbem_stdcharlen) :: tmp_str           ! Auxiliary variables for symmetry planes reading

  section_name='symmetry planes'
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'SEARCHING section ['//trim(section_name)//']')
  call fbem_search_section(fileunit,section_name,found)
  if (found) then

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START reading section ['//trim(section_name)//']')

    ! Initialize
    n_symplanes=0

    ! Find the "symmetry planes" section in the file
    call fbem_search_section(fileunit,section_name,found)
    ! Find the variable "x". The symmetry plane associated with the normal n=(1,0,0)
    call fbem_search_keyword(fileunit,'x','=',found)
    ! If "x" is found, read from the indicated file
    if (found) then
      ! Increment the counter of symmetry planes
      n_symplanes=n_symplanes+1
      ! Save the eid of the symmetry plane
      symplane_eid(n_symplanes)=1
      ! Read the signs variables at both sides of the symmetry plane
      read(fileunit,*) tmp_s, (tmp_t(i),i=1,problem%n)
      ! Copy to real data
      symplane_s(n_symplanes)=dble(tmp_s)
      do i=1,problem%n
        symplane_t(i,n_symplanes)=dble(tmp_t(i))
      end do
    end if

    ! Find the "symmetry planes" section in the file
    call fbem_search_section(fileunit,section_name,found)
    ! Find the variable "y". The symmetry plane associated with the normal n=(0,1,0)
    call fbem_search_keyword(fileunit,'y','=',found)
    ! If "y" is present, read from the indicated file
    if (found) then
      ! Increment the counter of symmetry planes
      n_symplanes=n_symplanes+1
      ! Save the eid of the symmetry plane
      symplane_eid(n_symplanes)=2
      ! Read the signs variables at both sides of the symmetry plane
      read(fileunit,*) tmp_s, (tmp_t(i),i=1,problem%n)
      ! Copy to real data
      symplane_s(n_symplanes)=dble(tmp_s)
      do i=1,problem%n
        symplane_t(i,n_symplanes)=dble(tmp_t(i))
      end do
    end if

    ! Find the "symmetry planes" section in the file
    call fbem_search_section(fileunit,section_name,found)
    ! Find the variable "z". The symmetry plane associated with the normal n=(0,0,1)
    call fbem_search_keyword(fileunit,'z','=',found)
    ! If "z" is present, read from the indicated file
    if (found) then
      ! If the problem is 3D, read
      if (problem%n.eq.3) then
        ! Increment the counter of symmetry planes
        n_symplanes=n_symplanes+1
        ! Save the eid of the symmetry plane
        symplane_eid(n_symplanes)=3
        ! Read the signs variables at both sides of the symmetry plane
        read(fileunit,*) tmp_s, (tmp_t(i),i=1,problem%n)
        ! Copy to real data
        symplane_s(n_symplanes)=dble(tmp_s)
        do i=1,problem%n
          symplane_t(i,n_symplanes)=dble(tmp_t(i))
        end do
      ! Otherwise, it is not valid
      else
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
                              'the "z" symmetry plane is valid only for 3D problems')
      end if
    end if

    ! Read "plane_n1"
    call fbem_search_section(fileunit,section_name,found)
    call fbem_search_keyword(fileunit,'plane_n1',':',found)
    if (.not.found) call fbem_search_keyword(fileunit,'plane_yz',':',found)
    if (found) then
      ! Increment the counter of symmetry planes
      n_symplanes=n_symplanes+1
      ! Save the eid of the symmetry plane
      symplane_eid(n_symplanes)=1
      read(fileunit,*) tmp_str
      if (trim(tmp_str).eq.'symmetry') then
        symplane_s(n_symplanes)  = 1.
        symplane_t(1,n_symplanes)=-1.
        symplane_t(2,n_symplanes)= 1.
        symplane_t(3,n_symplanes)= 1.
      else if (trim(tmp_str).eq.'antisymmetry') then
        symplane_s(n_symplanes)  =-1.
        symplane_t(1,n_symplanes)= 1.
        symplane_t(2,n_symplanes)=-1.
        symplane_t(3,n_symplanes)=-1.
      else
        call fbem_error_message(error_unit,0,section_name,0,'invalid argument for plane_n1')
      end if
    end if

    ! Read "plane_n2"
    call fbem_search_section(fileunit,section_name,found)
    call fbem_search_keyword(fileunit,'plane_n2',':',found)
    if (.not.found) call fbem_search_keyword(fileunit,'plane_zx',':',found)
    if (found) then
      ! Increment the counter of symmetry planes
      n_symplanes=n_symplanes+1
      ! Save the eid of the symmetry plane
      symplane_eid(n_symplanes)=2
      read(fileunit,*) tmp_str
      if (trim(tmp_str).eq.'symmetry') then
        symplane_s(n_symplanes)  = 1.
        symplane_t(1,n_symplanes)= 1.
        symplane_t(2,n_symplanes)=-1.
        symplane_t(3,n_symplanes)= 1.
      else if (trim(tmp_str).eq.'antisymmetry') then
        symplane_s(n_symplanes)  =-1.
        symplane_t(1,n_symplanes)=-1.
        symplane_t(2,n_symplanes)= 1.
        symplane_t(3,n_symplanes)=-1.
      else
        call fbem_error_message(error_unit,0,section_name,0,'invalid argument for plane_n2')
      end if
    end if

    ! If the problem is 3D, read
    if (problem%n.eq.3) then
      ! Read "plane_n3"
      call fbem_search_section(fileunit,section_name,found)
      call fbem_search_keyword(fileunit,'plane_n3',':',found)
      if (.not.found) call fbem_search_keyword(fileunit,'plane_xy',':',found)
      if (found) then
        ! Increment the counter of symmetry planes
        n_symplanes=n_symplanes+1
        ! Save the eid of the symmetry plane
        symplane_eid(n_symplanes)=3
        read(fileunit,*) tmp_str
        if (trim(tmp_str).eq.'symmetry') then
          symplane_s(n_symplanes)  = 1.
          symplane_t(1,n_symplanes)= 1.
          symplane_t(2,n_symplanes)= 1.
          symplane_t(3,n_symplanes)=-1.
        else if (trim(tmp_str).eq.'antisymmetry') then
          symplane_s(n_symplanes)  =-1.
          symplane_t(1,n_symplanes)=-1.
          symplane_t(2,n_symplanes)=-1.
          symplane_t(3,n_symplanes)= 1.
        else
          call fbem_error_message(error_unit,0,section_name,0,'invalid argument for plane_n1')
        end if
      end if
    end if

    ! Initialize symplane_iid
    symplane_iid=0
    ! Build eid to iid vector for symmetry planes
    do i=1,n_symplanes
      symplane_iid(symplane_eid(i))=i
    end do

    ! Internal order of symplane_eid(i), always from x to z, i.e. the following possibilities:
    ! symplane_eid(1) == x, symplane_eid(2) == z
    ! symplane_eid(1) == y, symplane_eid(2) == z
    ! But never:
    ! symplane_eid(1) == z, symplane_eid(2) == y
    ! symplane_eid(1) == y, symplane_eid(2) == x

    ! Number of symmetrical elements
    n_symelements=2**n_symplanes

    ! Write
    if (verbose_level.ge.3) then
      write(output_unit,'(a,i1)') '   Number of symmetry planes = ', n_symplanes
      if (n_symplanes.gt.0) then
        select case (problem%n)
          case (2)
            write(output_unit,'(3x,a,1x,a,1x,a,1x,a)') 'symplane', '_s', 't1', 't2'
          case (3)
            write(output_unit,'(3x,a,1x,a,1x,a,1x,a,1x,a)') 'symplane', '_s', 't1', 't2', 't3'
        end select
        write(fmtstr,*) '(3x,i8,1x,i2,',problem%n,'i3)'
        call fbem_trim2b(fmtstr)
        do i=1,n_symplanes
          write(output_unit,fmtstr) symplane_eid(i), idint(symplane_s(i)), (idint(symplane_t(j,i)),j=1,problem%n)
        end do
      end if
    end if

    ! Build symmetry multipliers
    ! Loop through symmetry planes
    do i=1,n_symplanes
      ! Geometrical symmetry, it depends on the symmetry plane
      select case (symplane_eid(i))
        ! Plane 1 (x2-x3): normal vector is e1, and pass through the point (0,0,0)
        case (1)
          symplane_m(1,i)=-1.d0
          symplane_m(2,i)= 1.d0
          symplane_m(3,i)= 1.d0
        ! Plane 2 (x3-x1): normal vector is e2, and pass through the point (0,0,0)
        case (2)
          symplane_m(1,i)= 1.d0
          symplane_m(2,i)=-1.d0
          symplane_m(3,i)= 1.d0
        ! Plane 3 (x1-x2): normal vector is e3, and pass through the point (0,0,0)
        case (3)
          symplane_m(1,i)= 1.d0
          symplane_m(2,i)= 1.d0
          symplane_m(3,i)=-1.d0
      end select
      ! Symmetry of rotations
      symplane_r(:,i)=-symplane_t(:,i)
    end do

    if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END reading section ['//trim(section_name)//']')

  else
    n_symplanes=0
    n_symelements=1
  end if

end subroutine read_symmetry_planes
