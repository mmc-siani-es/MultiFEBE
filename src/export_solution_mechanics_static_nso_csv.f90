! ---------------------------------------------------------------------
! Copyright (C) 2014-2026 Universidad de Las Palmas de Gran Canaria:
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

subroutine export_solution_mechanics_static_nso_csv

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_shape_functions
  use fbem_geometry
  use fbem_string_handling
  use fbem_numerical
  use fbem_fem_beams
  use fbem_fem_shells


  ! Module of problem variables
  use problem_variables

  ! other modules
  use csv_module

  ! No implicit variables allowed
  implicit none

  ! Local variables
  integer                                 :: k, kr, kb, ke, kn, kip, sip, kc
  integer                                 :: sb, se, se1, sn, ks, ss, sp
  logical                                 :: node_used(n_nodes)
  integer                                 :: face
  character(len=fbem_filename_max_length) :: tmp_filename
  type(csv_file)                          :: file_csv
  logical                                 :: status_ok, do_append
  integer                                 :: n_cols
  character(len=64)                       :: header(64)
  real(kind=real64)                       :: v1(3), v2(3)
  real(kind=real64)                       :: xi
  real(kind=real64), allocatable          :: N(:,:), nodal_axes(:,:,:)
  real(kind=real64), allocatable          :: ae(:)
  real(kind=real64)                       :: E, nu
  real(kind=real64)                       :: a6(6), f6(6)

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

          ! OPEN FILE HANDLING
          call file_csv%initialize
          write(tmp_filename,*) trim(output_filename),'.region_',region(kr)%id,'.boundary_',boundary(sb)%id,'.nso.csv'
          call fbem_trimall(tmp_filename)
          do_append=.false.

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

                ! CONFIGURE TABLE
                select case (problem%n)
                  case (2)
                    header(1) = 'Node'
                    header(2) = 'x1'
                    header(3) = 'x2'
                    header(4) = 'u1'
                    header(5) = 'u2'
                    header(6) = 't1'
                    header(7) = 't2'
                    n_cols=7
                  case (3)
                    header( 1) = 'Node'
                    header( 2) = 'x1'
                    header( 3) = 'x2'
                    header( 4) = 'x3'
                    header( 5) = 'u1'
                    header( 6) = 'u2'
                    header( 7) = 'u3'
                    header( 8) = 't1'
                    header( 9) = 't2'
                    header(10) = 't3'
                    n_cols=10
                end select
                ! OPEN FILE
                call file_csv%open(trim(tmp_filename),n_cols=n_cols,status_ok=status_ok,append=do_append)
                if (.not.status_ok) call fbem_error_message(error_unit,0,trim(tmp_filename),0,'error when opening the file.')
                ! HEADER
                do k=1,n_cols
                  call file_csv%add(header(k),trim_str=.true.)
                end do
                call file_csv%next_row()
                ! DATA
                node_used=.false.
                do ke=1,part(boundary(sb)%part)%n_elements
                  se=part(boundary(sb)%part)%element(ke)
                  do kn=1,element(se)%n_nodes
                    sn=element(se)%node(kn)
                    if ((.not.node_used(sn)).and.(node(sn)%export)) then
                      node_used(sn)=.true.
                      face=1
                      call file_csv%add(node(sn)%id)
                      do k=1,problem%n
                        call file_csv%add(element(se)%x_fn(k,kn))
                      end do
                      do k=1,2*problem%n
                        call file_csv%add(node(sn)%value_r(k,face))
                      end do
                      call file_csv%next_row()
                    end if
                  end do
                end do

              !
              ! Crack-like boundaries
              !
              case (fbem_boundary_class_cracklike)

                ! CONFIGURE TABLE
                select case (problem%n)
                  case (2)
                    header(1) = 'Face'
                    header(2) = 'Node'
                    header(3) = 'x1'
                    header(4) = 'x2'
                    header(5) = 'u1'
                    header(6) = 'u2'
                    header(7) = 't1'
                    header(8) = 't2'
                    n_cols=8
                  case (3)
                    header( 1) = 'Face'
                    header( 2) = 'Node'
                    header( 3) = 'x1'
                    header( 4) = 'x2'
                    header( 5) = 'x3'
                    header( 6) = 'u1'
                    header( 7) = 'u2'
                    header( 8) = 'u3'
                    header( 9) = 't1'
                    header(10) = 't2'
                    header(11) = 't3'
                    n_cols=11
                end select
                ! OPEN FILE
                call file_csv%open(trim(tmp_filename),n_cols=n_cols,status_ok=status_ok,append=do_append)
                if (.not.status_ok) call fbem_error_message(error_unit,0,trim(tmp_filename),0,'error when opening the file.')
                ! HEADER
                do k=1,n_cols
                  call file_csv%add(header(k),trim_str=.true.)
                end do
                call file_csv%next_row()
                ! DATA
                node_used=.false.
                do ke=1,part(boundary(sb)%part)%n_elements
                  se=part(boundary(sb)%part)%element(ke)
                  do kn=1,element(se)%n_nodes
                    sn=element(se)%node(kn)
                    if ((.not.node_used(sn)).and.(node(sn)%export)) then
                      node_used(sn)=.true.
                      do face=1,2
                        call file_csv%add(face)
                        call file_csv%add(node(sn)%id)
                        do k=1,problem%n
                          call file_csv%add(element(se)%x_fn(k,kn))
                        end do
                        do k=1,2*problem%n
                          call file_csv%add(node(sn)%value_r(k,face))
                        end do
                        call file_csv%next_row()
                      end do
                    end if
                  end do
                end do

            end select

          !
          ! BE-BE coupled boundary of BE-FE-BE coupled boundary
          !
          case (fbem_boundary_coupling_be_be,fbem_boundary_coupling_be_fe_be)

            ! CONFIGURE TABLE
            select case (problem%n)
              case (2)
                header(1) = 'Node'
                header(2) = 'x1'
                header(3) = 'x2'
                header(4) = 'u1'
                header(5) = 'u2'
                header(6) = 't1'
                header(7) = 't2'
                n_cols=7
              case (3)
                header( 1) = 'Node'
                header( 2) = 'x1'
                header( 3) = 'x2'
                header( 4) = 'x3'
                header( 5) = 'u1'
                header( 6) = 'u2'
                header( 7) = 'u3'
                header( 8) = 't1'
                header( 9) = 't2'
                header(10) = 't3'
                n_cols=10
            end select
            ! OPEN FILE
            call file_csv%open(trim(tmp_filename),n_cols=n_cols,status_ok=status_ok,append=do_append)
            if (.not.status_ok) call fbem_error_message(error_unit,0,trim(tmp_filename),0,'error when opening the file.')
            ! HEADER
            do k=1,n_cols
              call file_csv%add(header(k),trim_str=.true.)
            end do
            call file_csv%next_row()
            ! DATA
            node_used=.false.
            do ke=1,part(boundary(sb)%part)%n_elements
              se=part(boundary(sb)%part)%element(ke)
              do kn=1,element(se)%n_nodes
                sn=element(se)%node(kn)
                if ((.not.node_used(sn)).and.(node(sn)%export)) then
                  node_used(sn)=.true.
                  face=1
                  if (region(kr)%boundary_reversion(kb)) face=2
                  call file_csv%add(node(sn)%id)
                  do k=1,problem%n
                    call file_csv%add(element(se)%x_fn(k,kn))
                  end do
                  do k=1,2*problem%n
                    call file_csv%add(node(sn)%value_r(k,face))
                  end do
                  call file_csv%next_row()
                end if
              end do
            end do

          end select

          ! CLOSE FILE
          call file_csv%close(status_ok)
          if (.not.status_ok) call fbem_error_message(error_unit,0,trim(tmp_filename),0,'error when closing the file.')

        end do

      ! +-----------------+
      ! | INTERNAL POINTS |
      ! +-----------------+

      if (region(kr)%n_internalpoints > 0) then
        ! INITIALIZE FILE HANDLING
        call file_csv%initialize
        write(tmp_filename,*) trim(output_filename),'.region_',region(kr)%id,'.internal_points.nso.csv'
        call fbem_trimall(tmp_filename)
        do_append=.false.
        ! CONFIGURE TABLE
        select case (problem%n)
          case (2)
            header(1) = 'Internal point'
            header(2) = 'x1'
            header(3) = 'x2'
            header(4) = 'u1'
            header(5) = 'u2'
            header(6) = 's11'
            header(7) = 's12'
            header(8) = 's21'
            header(9) = 's22'
            n_cols=9
          case (3)
            header( 1) = 'Internal point'
            header( 2) = 'x1'
            header( 3) = 'x2'
            header( 4) = 'x3'
            header( 5) = 'u1'
            header( 6) = 'u2'
            header( 7) = 'u3'
            header( 8) = 's11'
            header( 9) = 's12'
            header(10) = 's13'
            header(11) = 's21'
            header(12) = 's22'
            header(13) = 's23'
            header(14) = 's31'
            header(15) = 's32'
            header(16) = 's33'
            n_cols=16
        end select
        ! OPEN FILE
        call file_csv%open(trim(tmp_filename),n_cols=n_cols,status_ok=status_ok,append=do_append)
        if (.not.status_ok) call fbem_error_message(error_unit,0,trim(tmp_filename),0,'error when opening the file.')
        ! HEADER
        do k=1,n_cols
          call file_csv%add(header(k),trim_str=.true.)
        end do
        call file_csv%next_row()
        ! DATA
        do kip=1,region(kr)%n_internalpoints
          sip=region(kr)%internalpoint(kip)
          if (internalpoint(sip)%export) then
            call file_csv%add(internalpoint(sip)%id)
            do k=1,problem%n
              call file_csv%add(internalpoint(sip)%x(k))
            end do
            do k=1,problem%n
              call file_csv%add(internalpoint(sip)%value_r(k,0))
            end do
            do kc=1,problem%n
              do k=1,problem%n
                call file_csv%add(internalpoint(sip)%value_r(k,kc))
              end do
            end do
            call file_csv%next_row()
          end if
        end do
        ! CLOSE FILE
        call file_csv%close(status_ok)
        if (.not.status_ok) call fbem_error_message(error_unit,0,trim(tmp_filename),0,'error when closing the file.')
      end if

      ! +------------------------+
      ! | ICOUPLED BE BODY LOADS |
      ! +------------------------+

      node_used=.false.
      do kb=1,region(kr)%n_be_bodyloads
        sb=region(kr)%be_bodyload(kb)
        sp=be_bodyload(sb)%part

        ! INITIALIZE FILE HANDLING
        call file_csv%initialize
        write(tmp_filename,*) trim(output_filename),'.region_',region(kr)%id,'.be_bodyload_',be_bodyload(sb)%id,'.nso.csv'
        call fbem_trimall(tmp_filename)
        do_append=.false.

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

            ! CONFIGURE TABLE
            select case (problem%n)
              case (2)
                header(1) = 'Node'
                header(2) = 'x1'
                header(3) = 'x2'
                header(4) = 'u1'
                header(5) = 'u2'
                header(6) = 'b1'
                header(7) = 'b2'
                n_cols=7
              case (3)
                header( 1) = 'Node'
                header( 2) = 'x1'
                header( 3) = 'x2'
                header( 4) = 'x3'
                header( 5) = 'u1'
                header( 6) = 'u2'
                header( 7) = 'u3'
                header( 8) = 'b1'
                header( 9) = 'b2'
                header(10) = 'b3'
                n_cols=10
            end select
            ! OPEN FILE
            call file_csv%open(trim(tmp_filename),n_cols=n_cols,status_ok=status_ok,append=do_append)
            if (.not.status_ok) call fbem_error_message(error_unit,0,trim(tmp_filename),0,'error when opening the file.')
            ! HEADER
            do k=1,n_cols
              call file_csv%add(header(k),trim_str=.true.)
            end do
            call file_csv%next_row()
            ! DATA
            do ke=1,part(sp)%n_elements
              se=part(sp)%element(ke)
              do kn=1,element(se)%n_nodes
                sn=element(se)%node(kn)
                if (.not.(node_used(sn)).and.(node(sn)%export)) then
                  node_used(sn)=.true.
                  face=1
                  call file_csv%add(node(sn)%id)
                  do k=1,problem%n
                    call file_csv%add(element(se)%x_fn(k,kn))
                  end do
                  do k=1,2*problem%n
                    call file_csv%add(node(sn)%value_r(k,face))
                  end do
                  call file_csv%next_row()
                end if
              end do
            end do

        end select

        ! CLOSE FILE
        call file_csv%close(status_ok)
        if (.not.status_ok) call fbem_error_message(error_unit,0,trim(tmp_filename),0,'error when closing the file.')

      end do

      ! ======================================================================================================================

      ! ======================================================================================================================
      ! FE region
      !
      case (fbem_fe)

        ! Mirar si interesa esto aqui o dentro del bucle, aunque con ello las distintas fe subregions que
        ! comparten nodos repitan resultados
        node_used=.false.

        do ks=1,region(kr)%n_fe_subregions
          ss=region(kr)%fe_subregion(ks)
          sp=fe_subregion(ss)%part
          se1=part(sp)%element(1)

          ! INITIALIZE FILE HANDLING
          call file_csv%initialize
          write(tmp_filename,*) trim(output_filename),'.region_',region(kr)%id,'.fe_subregion_',fe_subregion(ss)%id,'.nso.csv'
          call fbem_trimall(tmp_filename)
          do_append=.false.

          select case (element(se1)%n_dimension)

            ! ====================================================================================================================
            ! ONE-DIMENSIONAL ELEMENTS
            ! ====================================================================================================================

            case (1)

              select case (element(se1)%fe_type)

                case (0,5)

                  !---------------------------------------------------------------------------------------------------------------
                  ! DEGENERATED BEAM FINITE ELEMENT AND DISCRETE ROTATIONAL/TRANSLATIONAL SPRINGS/DASHPOTS (DISROTRA)
                  !
                  ! CONFIGURE TABLE
                  select case (problem%n)
                    case (2)
                      header(1) = 'Node'
                      header(2) = 'x1'
                      header(3) = 'x2'
                      header(4) = 'u1'
                      header(5) = 'u2'
                      header(6) = 'r3'
                      header(7) = 'f1'
                      header(8) = 'f2'
                      header(9) = 'm3'
                      n_cols=9
                    case (3)
                      header( 1) = 'Node'
                      header( 2) = 'x1'
                      header( 3) = 'x2'
                      header( 4) = 'x3'
                      header( 5) = 'u1'
                      header( 6) = 'u2'
                      header( 7) = 'u3'
                      header( 8) = 'r1'
                      header( 9) = 'r2'
                      header(10) = 'r3'
                      header(11) = 'f1'
                      header(12) = 'f2'
                      header(13) = 'f3'
                      header(14) = 'm1'
                      header(15) = 'm2'
                      header(16) = 'm3'
                      n_cols=16
                  end select
                  ! OPEN FILE
                  call file_csv%open(trim(tmp_filename),n_cols=n_cols,status_ok=status_ok,append=do_append)
                  if (.not.status_ok) call fbem_error_message(error_unit,0,trim(tmp_filename),0,'error when opening the file.')
                  ! HEADER
                  do k=1,n_cols
                    call file_csv%add(header(k),trim_str=.true.)
                  end do
                  call file_csv%next_row()
                  ! DATA
                  do ke=1,part(sp)%n_elements
                    se=part(sp)%element(ke)
                    do kn=1,element(se)%n_nodes
                      sn=element(se)%node(kn)
                      if ((.not.node_used(sn)).and.(node(sn)%export)) then
                        node_used(sn)=.true.
                        call file_csv%add(node(sn)%id)
                        do k=1,problem%n
                          call file_csv%add(node(sn)%x(k))
                        end do
                        do kc=1,2
                          do k=1,3*(problem%n-1)
                            call file_csv%add(node(sn)%value_r(k,kc))
                          end do
                        end do
                        call file_csv%next_row()
                      end if
                    end do
                  end do
                  !
                  !---------------------------------------------------------------------------------------------------------------

                case (1,2)

                  !---------------------------------------------------------------------------------------------------------------
                  ! STRAIGHT EULER-BERNOULLI AND TIMOSHENKO BEAM FINITE ELEMENTS
                  !
                  ! CONFIGURE TABLE
                  select case (problem%n)
                    case (2)
                      header(1) = 'Node'
                      header(2) = 'x1'
                      header(3) = 'x2'
                      header(4) = 'u1'
                      header(5) = 'u2'
                      header(6) = 'r3'
                      header(7) = 'f1'
                      header(8) = 'f2'
                      header(9) = 'm3'
                      n_cols=9
                    case (3)
                      header( 1) = 'Node'
                      header( 2) = 'x1'
                      header( 3) = 'x2'
                      header( 4) = 'x3'
                      header( 5) = 'u1'
                      header( 6) = 'u2'
                      header( 7) = 'u3'
                      header( 8) = 'r1'
                      header( 9) = 'r2'
                      header(10) = 'r3'
                      header(11) = 'f1'
                      header(12) = 'f2'
                      header(13) = 'f3'
                      header(14) = 'm1'
                      header(15) = 'm2'
                      header(16) = 'm3'
                      n_cols=16
                  end select
                  ! OPEN FILE
                  call file_csv%open(trim(tmp_filename),n_cols=n_cols,status_ok=status_ok,append=do_append)
                  if (.not.status_ok) call fbem_error_message(error_unit,0,trim(tmp_filename),0,'error when opening the file.')
                  ! HEADER
                  do k=1,n_cols
                    call file_csv%add(header(k),trim_str=.true.)
                  end do
                  call file_csv%next_row()
                  ! DATA
                  do ke=1,part(sp)%n_elements
                    se=part(sp)%element(ke)
                    ! If the element has nodes without rotations
                    if (element(se)%fe_options(1) == 0) then
                      ! Build the vector of element DOF and nodal axes
                      if (allocated(ae)) deallocate (ae)
                      allocate (ae(3*(problem%n-1)*element(se)%n_nodes))
                      ae=0
                      kc=0
                      do kn=1,element(se)%n_nodes
                        sn=element(se)%node(kn)
                        do k=1,node(sn)%n_dof
                          kc=kc+1
                          ae(kc)=node(sn)%value_r(k,1)
                        end do
                      end do
                      !
                      ! TO-DO: Los ejes nodales habría que guardarlos para cada nodo y construirlos para cada elemento
                      !
                      if (allocated(nodal_axes)) deallocate (nodal_axes)
                      allocate(nodal_axes(problem%n,problem%n,element(se)%n_nodes))
                      nodal_axes=0
                      do kn=1,element(se)%n_nodes
                        do kc=1,problem%n
                          nodal_axes(kc,kc,kn)=1
                        end do
                      end do
                      ! Material properties
                      E =region(kr)%property_r(5)
                      nu=region(kr)%property_r(3)
                    end if
                    do kn=1,element(se)%n_nodes
                      sn=element(se)%node(kn)
                      if ((.not.node_used(sn)).and.(node(sn)%export)) then
                        node_used(sn)=.true.
                        ! If the element has nodes without rotations, obtain displacements and rotations at the node
                        if (element(se)%fe_options(1) == 0) then
                          xi=element(se)%xi_gn(1,kn)
                          allocate (N(3*(problem%n-1),3*(problem%n-1)*element(se)%n_nodes))
                          N=fbem_fem_strbeam_N_w_rotations(problem%n,element(se)%type,element(se)%fe_options(1),element(se)%fe_type,&
                                                  element(se)%x_gn,element(se)%ep,element(se)%A,element(se)%I(1:3),element(se)%ksh,&
                                                  E,nu,nodal_axes,xi)
                          a6=0
                          a6(1:3*(problem%n-1)) = matmul(N,ae)
                          f6=0
                          f6=node(sn)%value_r(1:node(sn)%n_dof,2)
                          deallocate (N)
                        else
                          a6(1:3*(problem%n-1)) = node(sn)%value_r(1:3*(problem%n-1),1)
                          f6(1:3*(problem%n-1)) = node(sn)%value_r(1:3*(problem%n-1),2)
                        end if
                        ! EXPORT
                        call file_csv%add(node(sn)%id)
                        do k=1,problem%n
                          call file_csv%add(node(sn)%x(k))
                        end do
                        do k=1,3*(problem%n-1)
                          call file_csv%add(a6(k))
                        end do
                        do k=1,3*(problem%n-1)
                          call file_csv%add(f6(k))
                        end do
                        call file_csv%next_row()
                      end if
                    end do

                  end do
                  !
                  !---------------------------------------------------------------------------------------------------------------

                case (3,4,6)

                  !---------------------------------------------------------------------------------------------------------------
                  ! BAR FINITE ELEMENT, DISCRETE TRANSLATIONAL SPRINGS/DASHPOTS (DISTRA), AND DISCRETE SPRING-DASHPOT
                  !
                  ! CONFIGURE TABLE
                  select case (problem%n)
                    case (2)
                      header(1) = 'Node'
                      header(2) = 'x1'
                      header(3) = 'x2'
                      header(4) = 'u1'
                      header(5) = 'u2'
                      header(7) = 'f1'
                      header(8) = 'f2'
                      n_cols=7
                    case (3)
                      header( 1) = 'Node'
                      header( 2) = 'x1'
                      header( 3) = 'x2'
                      header( 4) = 'x3'
                      header( 5) = 'u1'
                      header( 6) = 'u2'
                      header( 7) = 'u3'
                      header( 8) = 'f1'
                      header( 9) = 'f2'
                      header(10) = 'f3'
                      n_cols=10
                  end select
                  ! OPEN FILE
                  call file_csv%open(trim(tmp_filename),n_cols=n_cols,status_ok=status_ok,append=do_append)
                  if (.not.status_ok) call fbem_error_message(error_unit,0,trim(tmp_filename),0,'error when opening the file.')
                  ! HEADER
                  do k=1,n_cols
                    call file_csv%add(header(k),trim_str=.true.)
                  end do
                  call file_csv%next_row()
                  ! DATA
                  do ke=1,part(sp)%n_elements
                    se=part(sp)%element(ke)
                    do kn=1,element(se)%n_nodes
                      sn=element(se)%node(kn)
                      if ((.not.node_used(sn)).and.(node(sn)%export)) then
                        node_used(sn)=.true.
                        call file_csv%add(node(sn)%id)
                        do k=1,problem%n
                          call file_csv%add(node(sn)%x(k))
                        end do
                        do kc=1,2
                          do k=1,problem%n
                            call file_csv%add(node(sn)%value_r(k,kc))
                          end do
                        end do
                        call file_csv%next_row()
                      end if
                    end do
                  end do
                  !
                  !---------------------------------------------------------------------------------------------------------------

                case default

                  !---------------------------------------------------------------------------------------------------------------
                  ! OTHER TYPES
                  !
                  call fbem_error_message(error_unit,0,'fe_subregion',fe_subregion(ss)%id,'invalid type of element')
                  !
                  !---------------------------------------------------------------------------------------------------------------

              end select

            ! ====================================================================================================================
            ! TWO-DIMENSIONAL ELEMENTS
            ! ====================================================================================================================

            case (2)

              select case (problem%n)

                case (2)

                  !---------------------------------------------------------------------------------------------------------------
                  ! SOLID / CONTINUUM ELEMENTS
                  !
                  ! CONFIGURE TABLE
                  header(1) = 'Node'
                  header(2) = 'x1'
                  header(3) = 'x2'
                  header(4) = 'u1'
                  header(5) = 'u2'
                  header(6) = 'f1'
                  header(7) = 'f2'
                  n_cols=7
                  ! OPEN FILE
                  call file_csv%open(trim(tmp_filename),n_cols=n_cols,status_ok=status_ok,append=do_append)
                  if (.not.status_ok) call fbem_error_message(error_unit,0,trim(tmp_filename),0,'error when opening the file.')
                  ! HEADER
                  do k=1,n_cols
                    call file_csv%add(header(k),trim_str=.true.)
                  end do
                  call file_csv%next_row()
                  ! DATA
                  do ke=1,part(sp)%n_elements
                    se=part(sp)%element(ke)
                    do kn=1,element(se)%n_nodes
                      sn=element(se)%node(kn)
                      if ((.not.node_used(sn)).and.(node(sn)%export)) then
                        node_used(sn)=.true.
                        call file_csv%add(node(sn)%id)
                        call file_csv%add(node(sn)%x(1))
                        call file_csv%add(node(sn)%x(2))
                        call file_csv%add(node(sn)%value_r(1,1))
                        call file_csv%add(node(sn)%value_r(2,1))
                        call file_csv%add(node(sn)%value_r(1,2))
                        call file_csv%add(node(sn)%value_r(2,2))
                        call file_csv%next_row()
                      end if
                    end do
                  end do
                  !
                  !---------------------------------------------------------------------------------------------------------------

                case (3)

                  !---------------------------------------------------------------------------------------------------------------
                  ! DEGENERATED SHELL FINITE ELEMENT
                  !
                  ! CONFIGURE TABLE
                  header( 1) = 'Node'
                  header( 2) = 'x1'
                  header( 3) = 'x2'
                  header( 4) = 'x3'
                  header( 5) = 'u1'
                  header( 6) = 'u2'
                  header( 7) = 'u3'
                  header( 8) = 'r1'
                  header( 9) = 'r2'
                  header(10) = 'r3'
                  header(11) = 'f1'
                  header(12) = 'f2'
                  header(13) = 'f3'
                  header(14) = 'm1'
                  header(15) = 'm2'
                  header(16) = 'm3'
                  n_cols=16
                  ! OPEN FILE
                  call file_csv%open(trim(tmp_filename),n_cols=n_cols,status_ok=status_ok,append=do_append)
                  if (.not.status_ok) call fbem_error_message(error_unit,0,trim(tmp_filename),0,'error when opening the file.')
                  ! HEADER
                  do k=1,n_cols
                    call file_csv%add(header(k),trim_str=.true.)
                  end do
                  call file_csv%next_row()
                  ! DATA
                  do ke=1,part(sp)%n_elements
                    se=part(sp)%element(ke)
                    do kn=1,element(se)%n_nodes
                      sn=element(se)%node(kn)
                      if ((.not.node_used(sn)).and.(node(sn)%export)) then
                        node_used(sn)=.true.
                        ! For 5 DOF nodes, transform to 6 DOF node for export.
                        if (node(sn)%n_dof == 5) then
                          v1 = element(se)%v_midnode(:,1,kn)
                          v2 = element(se)%v_midnode(:,2,kn)
                          a6 = fbem_fem_degshell_XDOF_to_6DOF(node(sn)%value_r(1:5,1),v1,v2,5)
                          f6 = fbem_fem_degshell_XDOF_to_6DOF(node(sn)%value_r(1:5,2),v1,v2,5)
                        else
                          a6 = node(sn)%value_r(1:6,1)
                          f6 = node(sn)%value_r(1:6,2)
                        end if
                        ! EXPORT
                        call file_csv%add(node(sn)%id)
                        do k=1,3
                          call file_csv%add(node(sn)%x(k))
                        end do
                        do k=1,6
                          call file_csv%add(a6)
                        end do
                        do k=1,6
                          call file_csv%add(f6)
                        end do
                        call file_csv%next_row()
                      end if
                    end do
                  end do
                  !
                  !---------------------------------------------------------------------------------------------------------------

              end select

            ! ====================================================================================================================
            ! THREE-DIMENSIONAL ELEMENTS
            ! ====================================================================================================================

            case (3)

              !-------------------------------------------------------------------------------------------------------------------
              ! SOLID / CONTINUUM ELEMENTS
              !
              ! CONFIGURE TABLE
              header( 1) = 'Node'
              header( 2) = 'x1'
              header( 3) = 'x2'
              header( 4) = 'x3'
              header( 5) = 'u1'
              header( 6) = 'u2'
              header( 7) = 'u3'
              header( 8) = 'f1'
              header( 9) = 'f2'
              header(10) = 'f3'
              n_cols=10
              ! OPEN FILE
              call file_csv%open(trim(tmp_filename),n_cols=n_cols,status_ok=status_ok,append=do_append)
              if (.not.status_ok) call fbem_error_message(error_unit,0,trim(tmp_filename),0,'error when opening the file.')
              ! HEADER
              do k=1,n_cols
                call file_csv%add(header(k),trim_str=.true.)
              end do
              call file_csv%next_row()
              ! DATA
              do ke=1,part(sp)%n_elements
                se=part(sp)%element(ke)
                do kn=1,element(se)%n_nodes
                  sn=element(se)%node(kn)
                  if ((.not.node_used(sn)).and.(node(sn)%export)) then
                    node_used(sn)=.true.
                    call file_csv%add(node(sn)%id)
                    call file_csv%add(node(sn)%x(1))
                    call file_csv%add(node(sn)%x(2))
                    call file_csv%add(node(sn)%x(3))
                    call file_csv%add(node(sn)%value_r(1,1))
                    call file_csv%add(node(sn)%value_r(2,1))
                    call file_csv%add(node(sn)%value_r(3,1))
                    call file_csv%add(node(sn)%value_r(1,2))
                    call file_csv%add(node(sn)%value_r(2,2))
                    call file_csv%add(node(sn)%value_r(3,2))
                    call file_csv%next_row()
                  end if
                end do
              end do
              !
              !-------------------------------------------------------------------------------------------------------------------

          end select

          ! CLOSE FILE
          call file_csv%close(status_ok)
          if (.not.status_ok) call fbem_error_message(error_unit,0,trim(tmp_filename),0,'error when closing the file.')

        end do

      ! ======================================================================================================================

    end select

  end do

end subroutine export_solution_mechanics_static_nso_csv
