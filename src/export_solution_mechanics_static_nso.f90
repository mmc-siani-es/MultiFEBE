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

subroutine export_solution_mechanics_static_nso(output_fileunit)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_shape_functions
  use fbem_geometry
  use fbem_string_handling
  use fbem_numerical

  ! Module of problem variables
  use problem_variables

  ! No implicit variables allowed
  implicit none

  ! I/O variables
  integer                                 :: output_fileunit
 ! Local variables
  integer                                 :: k
  integer                                 :: kr, kb, ke, kn, kip, sip, kc
  integer                                 :: sb, se, sn, ks, ss, sp
  logical                                 :: node_used(n_nodes)
  character(len=fbem_fmtstr)              :: fmt1, fmt2
  integer                                 :: ncint, ncreal, nc, ncmax
  character(len=fbem_string_max_length)   :: tmp_string
  integer                                 :: face, k_start, k_end

  ! ====== !
  ! HEADER !
  ! ====== !

  if (export_overwrite) then
    write(output_fileunit,'(a)'  ) '# Program      : multifebe'
    write(output_fileunit,'(a)'  ) '# Version      : 1.0'
    write(output_fileunit,'(a)'  ) '# File_format  : nso'
    write(output_fileunit,'(a)'  ) '# Specification: 1'
    write(output_fileunit,'(a,a)') '# Input_file   : ', trim(input_filename)
    write(output_fileunit,'(a,i1)')'# Problem n    : ', problem%n
    write(output_fileunit,'(a,a)') '# Description  : ', trim(problem%description)
    write(tmp_string,*) timestamp_date_start(1:4),'-',timestamp_date_start(5:6),'-',timestamp_date_start(7:8),' ',&
                        timestamp_time_start(1:2),':',timestamp_time_start(3:4),':',timestamp_time_start(5:10)
    call fbem_trim2b(tmp_string)
    write(output_fileunit,'(a,a)') '# Timestamp    : ', trim(tmp_string)
    write(output_fileunit,*)
    ! Column description
    write(output_fileunit,'(a)'  ) '# Columns  Description'
    write(output_fileunit,'(a)'  ) '# 1-2      Step index and value.'
    write(output_fileunit,'(a)'  ) '# 3-5      Region id, class and type.'
    write(output_fileunit,'(a)'  ) '# 6-8      (if col 4 == 1) Boundary id, class and face.'
    write(output_fileunit,'(a)'  ) '# 6-8      (if col 4 == 2) Subregion id, number of DOF and 0.'
    select case (problem%n)
    case (2)
    write(output_fileunit,'(a)'  ) '# 9-11     Node id, x1 and x2.'
    write(output_fileunit,'(a)'  ) '# >=12     Node variables. Depend on the region class and type (see documentation).'
    case (3)
    write(output_fileunit,'(a)'  ) '# 9-12     Node id, x1, x2 and x3.'
    write(output_fileunit,'(a)'  ) '# >=13     Node variables. Depend on the region class and type (see documentation).'
    end select
    write(output_fileunit,'(a)') '#'
  end if

  !
  ! Notation:
  !
  ! Elastic:
  !   - uk: displacement (k component)
  !   - tk: traction (k component)
  !   -
  !
  !
  ! 2D PROBLEM
  ! ----------
  !
  ! BE region ($4 == 1):
  !
  ! |-----------------|     |-----|-----|-----|-----|
  ! | 5               |     |12,13|14,15|16,17|18,19|
  ! |-----------------| ... |-----|-----|-----|-----|
  ! | 2 (elastic)     |     |   u1|   u2|   t1|   t2|
  ! |-----------------|     |-----|-----|-----|-----|
  !
  ! FE region ($4 == 2):
  !
  ! |-----------------|     |-----|-----|-----|-----|
  ! | 5               |     |12,13|14,15|16,17|18,19|
  ! |-----------------| ... |-----|-----|-----|-----|
  ! | 2 (elastic)     |     |   u1|   u2|   F1|   F2|
  ! |-----------------|     |-----|-----|-----|-----|
  !
  ! 3D PROBLEM
  ! ----------
  !
  ! BE region ($4 == 1):
  !
  ! |-----------------|     |-----|-----|-----|-----|-----|-----|
  ! | 5               |     |13,14|15,16|17,18|19,20|21,22|23,24|
  ! |-----------------| ... |-----|-----|-----|-----|-----|-----|
  ! | 2 (elastic)     |     |   u1|   u2|   u3|   t1|   t2|   t3|
  ! |-----------------|     |-----|-----|-----|-----|-----|-----|

  ! --------------
  ! Export formats
  ! --------------

  ! Export format for columns 1-11 (2D) or 1-12 (3D)
  write(fmt1,*) '(1',fmt_integer,',1',fmt_real,',7',fmt_integer,',',problem%n,fmt_real,')'
  call fbem_trim2b(fmt1)
  ! Export format for columns >=12 (2D) or >=13 (3D)
  write(fmt2,*) '(',fmt_real,')'
  call fbem_trim2b(fmt2)

  ! Loop through regions
  do kr=1,n_regions
    select case (region(kr)%class)

      ! ======================================================================================================================
      ! BE region
      !
      case (fbem_be)

        ! +----------------+
        ! | BOUNDARY NODES |
        ! +----------------+

        do kb=1,region(kr)%n_boundaries
          sb=region(kr)%boundary(kb)
          node_used=.false.
          do ke=1,part(boundary(sb)%part)%n_elements
            se=part(boundary(sb)%part)%element(ke)
            do kn=1,element(se)%n_nodes
              sn=element(se)%node(kn)
              if ((node_used(sn).eqv.(.false.)).and.(node(sn)%export)) then
                node_used(sn)=.true.
                !
                ! Boundary coupling
                !
                select case (boundary(sb)%coupling)
                !
                ! Uncoupled boundary of BE-FE coupled boundary
                !
                case (fbem_boundary_coupling_be,fbem_boundary_coupling_be_fe)
                    !
                    ! Boundary class
                    !
                    select case (boundary(sb)%class)
                      !
                      ! Ordinary boundary
                      !
                      case (fbem_boundary_class_ordinary)

                        face=1
                        write(output_fileunit,fmt1,advance='no') 0, 0.d0, &
                                                                 region(kr)%id, region(kr)%class, region(kr)%type, &
                                                                 boundary(sb)%id, boundary(sb)%class, face, &
                                                                 node(sn)%id, element(se)%x_fn(:,kn)
                        do k=1,2*problem%n
                          write(output_fileunit,fmt2,advance='no') node(sn)%value_r(k,face)
                        end do
                        write(output_fileunit,*)

                      !
                      ! Crack-like boundaries
                      !
                      case (fbem_boundary_class_cracklike)

                        do face=1,2
                          write(output_fileunit,fmt1,advance='no') 0, 0.d0, &
                                                                   region(kr)%id, region(kr)%class, region(kr)%type, &
                                                                   boundary(sb)%id, boundary(sb)%class, face, &
                                                                   node(sn)%id, element(se)%x_fn(:,kn)
                          do k=1,2*problem%n
                            write(output_fileunit,fmt2,advance='no') node(sn)%value_r(k,face)
                          end do
                          write(output_fileunit,*)
                        end do

                    end select
                !
                ! BE-BE coupled boundary of BE-FE-BE coupled boundary
                !
                case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)

                  if (region(kr)%boundary_reversion(kb).eqv.(.false.)) then
                    face=1
                  else
                    face=2
                  end if
                  write(output_fileunit,fmt1,advance='no') 0, 0.d0, &
                                                           region(kr)%id, region(kr)%class, region(kr)%type, &
                                                           boundary(sb)%id, boundary(sb)%class, face, &
                                                           node(sn)%id, element(se)%x_fn(:,kn)
                  do k=1,2*problem%n
                    write(output_fileunit,fmt2,advance='no') node(sn)%value_r(k,face)
                  end do
                  write(output_fileunit,*)

                end select
              end if
            end do
          end do
        end do


      ! +-----------------+
      ! | INTERNAL POINTS |
      ! +-----------------+

      ! Loop through the internal points of the region
      do kip=1,region(kr)%n_internalpoints
        sip=region(kr)%internalpoint(kip)
        if (internalpoint(sip)%export) then
          ! Write columns 1-11 (2D) or 1-12 (3D)
          write(output_fileunit,fmt1,advance='no') 0, 0.d0, &
                                                   region(kr)%id, region(kr)%class, region(kr)%type, &
                                                   0, 0, 0, &
                                                   internalpoint(sip)%id, internalpoint(sip)%x

          ! Displacements
          do k=1,problem%n
            write(output_fileunit,fmt2,advance='no') internalpoint(sip)%value_r(k,0)
          end do
          ! Stresses
          do kc=1,problem%n
            do k=1,problem%n
              write(output_fileunit,fmt2,advance='no') internalpoint(sip)%value_r(k,kc)
            end do
          end do

          write(output_fileunit,*)
        end if
      end do

      ! +------------------------+
      ! | ICOUPLED BE BODY LOADS |
      ! +------------------------+

      face=1
      node_used=.false.
      do kb=1,region(kr)%n_be_bodyloads
        sb=region(kr)%be_bodyload(kb)
        sp=be_bodyload(sb)%part
        select case (be_bodyload(sb)%coupling)

          ! ----------------------------------
          ! FE BEAM TIP - BE LINE/SURFACE LOAD
          ! ----------------------------------

          case (fbem_bl_coupling_beam_tip)
            stop 'not yet'

          ! -----------------------
          ! FE SHELL - BE EDGE LOAD
          ! -----------------------

          case (fbem_bl_coupling_shell_edge)
            stop 'not yet'

          ! -----------------------------------------------------
          ! FE BEAM - BE LINE LOAD AND FE SHELL - BE SURFACE LOAD
          ! -----------------------------------------------------

          case (fbem_bl_coupling_beam_line,fbem_bl_coupling_shell_surface)
            do ke=1,part(sp)%n_elements
              se=part(sp)%element(ke)
              do kn=1,element(se)%n_nodes
                sn=element(se)%node(kn)
                if (node_used(sn).eqv.(.false.).and.(node(sn)%export)) then
                  node_used(sn)=.true.
                  ! Write columns 1-11 (2D) or 1-12 (3D)
                  write(output_fileunit,fmt1,advance='no') 0, 0.d0, &
                                                           region(kr)%id, region(kr)%class, region(kr)%type, &
                                                           be_bodyload(sb)%id, 0, 0, &
                                                           node(sn)%id, element(se)%x_fn(:,kn)
                  ! Write columns >=12 (2D) or >=13 (3D)
                  do k=1,2*problem%n
                    write(output_fileunit,fmt2,advance='no') node(sn)%value_r(k,face)
                  end do
                  write(output_fileunit,*)
                end if
              end do
            end do

        end select

      end do



      ! ======================================================================================================================

      ! ======================================================================================================================
      ! FE region
      !
      case (fbem_fe)
        node_used=.false.
        do ks=1,region(kr)%n_fe_subregions
          ss=region(kr)%fe_subregion(ks)
          do ke=1,part(fe_subregion(ss)%part)%n_elements
            se=part(fe_subregion(ss)%part)%element(ke)
            do kn=1,element(se)%n_nodes
              sn=element(se)%node(kn)
              if ((.not.node_used(sn)).and.(node(sn)%export)) then
                node_used(sn)=.true.
                ! Write columns 1-11 (2D) or 1-12 (3D)
                write(output_fileunit,fmt1,advance='no') 0, 0.d0, region(kr)%id, region(kr)%class, region(kr)%type, &
                                                         fe_subregion(ss)%id, node(sn)%n_dof, 0, &
                                                         node(sn)%id, node(sn)%x
                ! Write columns >=12 (2D) or >=13 (3D)
                do k=1,node(sn)%n_dof
                  write(output_fileunit,fmt2,advance='no') node(sn)%value_r(k,1)
                end do
                do k=1,node(sn)%n_dof
                  write(output_fileunit,fmt2,advance='no') node(sn)%value_r(k,2)
                end do
                write(output_fileunit,*)
              end if
            end do
          end do
        end do


      ! ======================================================================================================================

    end select
    write(output_fileunit,*)
    write(output_fileunit,*)
  end do

end subroutine export_solution_mechanics_static_nso
