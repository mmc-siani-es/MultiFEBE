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


subroutine export_solution_laplace

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_numerical
  use fbem_shape_functions
  use fbem_geometry
  use fbem_symmetry

  ! Module of problem variables
  use problem_variables

  ! No implicit variables allowed
  implicit none

  ! Local variables
  integer                                 :: kr, kb, ke, kn, kc
  integer                                 :: sb, se, sn
  integer                                 :: kip, sip
  integer                                 :: k
  integer                                 :: output_fileunit
  character(len=fbem_filename_max_length) :: tmp_filename
  logical                                 :: node_used(n_nodes)
  character(len=fbem_fmtstr)              :: fmt1, fmt2
  integer                                 :: ncint, ncreal, nc, ncmax
  integer                                 :: face
  character(len=fbem_string_max_length)   :: tmp_string

  ! Starting message
  if (verbose_level.ge.1)  write(output_unit,'(a)') 'Exporting solution ...'

  ! ================================================================================================================================
  ! EXPORT NODAL SOLUTIONS
  ! ================================================================================================================================

  if (export_nso) then

    ! Open ".nso" file
    ! ----------------

    output_fileunit=fbem_get_valid_unit()
    tmp_filename=trim(output_filename)//'.nso'
    call fbem_trim2b(tmp_filename)
    open(unit=output_fileunit,file=trim(tmp_filename),action='write',recl=fbem_file_record_length)
    if (verbose_level.ge.2) write(output_unit,'(a,a,a)') ' Opening "', trim(tmp_filename), '"'

    ! Header of the file
    ! ------------------
    write(output_fileunit,'(a)'  ) '# Program      : multifebe'
    write(output_fileunit,'(a,a5)')'# Version      : ', multifebe_version
    write(output_fileunit,'(a)'  ) '# File_format  : nso'
    write(output_fileunit,'(a,i1)')'# Problem_dim  : ', problem%n
    write(output_fileunit,'(a,a)') '# Input_file   : ', trim(input_filename)
    write(output_fileunit,'(a,a)') '# Description  : ', trim(problem%description)
    write(output_fileunit,'(a)',advance='no') '# Timestamp    : '
    call fbem_timestamp(output_fileunit,timestamp_date_start,timestamp_time_start,1)
    write(output_fileunit,*)
    write(output_fileunit,'(a)'  ) '# Columns  Description'
    write(output_fileunit,'(a)'  ) '# C1-C2      Region id and class.'
    write(output_fileunit,'(a)'  ) '# C3-C5      Boundary id, class and face.'
    select case (problem%n)
    case (2)
    write(output_fileunit,'(a)'  ) '# C6-C8      Node id, x1 and x2.'
    write(output_fileunit,'(a)'  ) '# >=C9       See manual.'

    case (3)
    write(output_fileunit,'(a)'  ) '# C6-C9      Node id, x1, x2 and x3.'
    write(output_fileunit,'(a)'  ) '# >=C10      See manual.'
    end select
    write(output_fileunit,'(a)') '#'
    ! Column header
    ! Number of characters of an integer column
    write(fmt1,*) '(',fmt_integer,')'
    call fbem_trimall(fmt1)
    write(tmp_string,fmt1) 0
    ncint=len_trim(tmp_string)
    ! Number of characters of a real column
    write(fmt1,*) '(',fmt_real,')'
    call fbem_trimall(fmt1)
    write(tmp_string,fmt1) 0.
    ncreal=len_trim(tmp_string)
    ! Comment character
    write(output_fileunit,'(a)',advance='no') '#'
    ! First column
    do k=1,ncint-3
      write(output_fileunit,'(a)',advance='no') '_'
    end do
    write(output_fileunit,'(a2)',advance='no') 'C1'
    ! Rest of the columns
    select case (problem%n)
      case (2)
        ncmax=11
      case (3)
        ncmax=13
    end select
    do kc=2,ncmax
      ! Depending if integer or real column
      if ((kc.ge.2).and.(kc.le.6)) then
        nc=ncint
      else
        nc=ncreal
      end if
      ! Write
      do k=1,nc-fbem_nchar_int(kc)-1
        write(output_fileunit,'(a)',advance='no') '_'
      end do
      write(fmt1,*) '(a1,i',fbem_nchar_int(kc),')'
      call fbem_trimall(fmt1)
      write(output_fileunit,fmt1,advance='no') 'C',kc
    end do
    write(output_fileunit,*)

    ! Export formats
    ! --------------

    ! Export format for columns 1-6
    write(fmt1,*) '(6',fmt_integer,',',problem%n,fmt_real,')'
    call fbem_trim2b(fmt1)
    ! Export format for columns >=10
    write(fmt2,*) '(1',fmt_real,')'
    call fbem_trim2b(fmt2)

    ! Write
    ! -----

    ! Initialize
    node_used=.false.
    ! Loop through regions
    do kr=1,n_regions
      ! Switch between region classes
      select case (region(kr)%class)

        ! --------------------------------------------------------------------------------------------------------------------------
        ! BE REGION
        !
        case (fbem_be)

          ! +----------------+
          ! | BOUNDARY NODES |
          ! +----------------+

          ! Loop through entities
          do kb=1,region(kr)%n_boundaries
            sb=region(kr)%boundary(kb)
            ! Initialize node usage
            do kn=1,part(boundary(sb)%part)%n_nodes
              sn=part(boundary(sb)%part)%node(kn)
              node_used(sn)=.false.
            end do
            ! Loop through elements and nodes
            do ke=1,part(boundary(sb)%part)%n_elements
              se=part(boundary(sb)%part)%element(ke)
              do kn=1,element(se)%n_nodes
                sn=element(se)%node(kn)
                if (node_used(sn).eqv.(.false.)) then
                  node_used(sn)=.true.
                  select case (boundary(sb)%coupling)
                    !
                    ! Uncoupled boundary of BE-FE coupled boundary
                    !
                    case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
                      select case (boundary(sb)%class)
                        !
                        ! Ordinary boundary
                        !
                        case (fbem_boundary_class_ordinary)
                          face=1
                          ! Write columns 1-9
                          write(output_fileunit,fmt1,advance='no') region(kr)%id, region(kr)%class, &
                                                                   boundary(sb)%id, boundary(sb)%class, face, &
                                                                   node(sn)%id, element(se)%x_fn(:,kn)
                          ! Write columns >=10
                          write(output_fileunit,fmt2,advance='no') node(sn)%value_r(1,face)
                          write(output_fileunit,fmt2,advance='no') node(sn)%value_r(2,face)
                          write(output_fileunit,*)
                        !
                        ! Crack-like boundaries
                        !
                        case (fbem_boundary_class_cracklike)
                          ! Face + and face -
                          do face=1,2
                            ! Write columns 1-9
                            write(output_fileunit,fmt1,advance='no') region(kr)%id, region(kr)%class, &
                                                                     boundary(sb)%id, boundary(sb)%class, face, &
                                                                     node(sn)%id, element(se)%x_fn(:,kn)
                            ! Write columns >=10
                            write(output_fileunit,fmt2,advance='no') node(sn)%value_r(1,face)
                            write(output_fileunit,fmt2,advance='no') node(sn)%value_r(2,face)
                            write(output_fileunit,*)
                          end do
                      end select
                    !
                    ! BE-BE coupled boundary of BE-FE-BE coupled boundary
                    !
                    case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)
                      ! The region is the region 1 of the boundary
                      if (region(kr)%boundary_reversion(kb).eqv.(.false.)) then
                        face=1
                      ! The region is the region 2 of the boundary
                      else
                        face=2
                      end if
                      ! Write columns 1-9
                      write(output_fileunit,fmt1,advance='no') region(kr)%id, region(kr)%class, &
                                                               boundary(sb)%id, boundary(sb)%class, face, &
                                                               node(sn)%id, element(se)%x_fn(:,kn)
                      ! Write columns >=10
                      write(output_fileunit,fmt2,advance='no') node(sn)%value_r(1,face)
                      write(output_fileunit,fmt2,advance='no') node(sn)%value_r(2,face)
                      write(output_fileunit,*)
                  end select
                end if
              end do
            end do
          end do

          ! +----------------+
          ! | INTERNAL NODES |
          ! +----------------+

          ! Loop through the internal points of the region
          do kip=1,region(kr)%n_internalpoints
            sip=region(kr)%internalpoint(kip)
            ! Write columns 1-11 (2D) or 1-12 (3D)
            write(output_fileunit,fmt1,advance='no') region(kr)%id, region(kr)%class, &
                                                     0, 0, 0, &
                                                     internalpoint(sip)%id, internalpoint(sip)%x
            ! Write columns >=12 (2D) or >=13 (3D)
            ! PRIMARY VARIABLES
            write(output_fileunit,fmt2,advance='no') internalpoint(sip)%value_r(1,0)
            ! SECONDARY VARIABLES
            do kc=1,problem%n
              write(output_fileunit,fmt2,advance='no') internalpoint(sip)%value_r(1,kc)
            end do
            write(output_fileunit,*)
          end do

        ! --------------------------------------------------------------------------------------------------------------------------

        ! --------------------------------------------------------------------------------------------------------------------------
        ! FE region
        !
        case (fbem_fe)
          stop 'not yet 85'
        ! --------------------------------------------------------------------------------------------------------------------------

      end select
    end do

    ! Close file
    close(unit=output_fileunit)
    if (verbose_level.ge.2) write(output_unit,'(a,a,a)') ' Closing "', trim(tmp_filename), '"'

  end if

  ! Ending message
  if (verbose_level.ge.1) write(output_unit,'(a)') 'done.'

end subroutine export_solution_laplace
