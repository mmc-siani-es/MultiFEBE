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

!
!   2D PROBLEM: BEAM FE (1D)  - BE BODY LINE LOAD BE (1D)
!   3D PROBLEM: BEAM FE (1D)  - BE BODY SURFACE LOAD BE (1D)
!   3D PROBLEM: SHELL FE (2D) - BE BODY SURFACE LOAD BE (2D)
!

subroutine build_connectivities_be_bodyload_elements_and_fe

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
  integer                               :: i, j, k
  integer                               :: sn, se, sse
  logical                               :: coupled

  ! Starting message
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START building connectivities between BE bodyloads and FE')

  ! Loop through BE body loads
  do i=1,n_be_bodyloads
    ! Build connectivities of the first element of the BE body load
    se=part(be_bodyload(i)%part)%element(1)
    call build_connectivity_be_bodyload_element_fe(se)
    select case (element(se)%n_dimension)
      !
      ! POINT LOAD
      !
      case (0)
        ! True if coupled
        if (node(element(se)%node(1))%coupled_node.eq.0) then
          coupled=.false.
          be_bodyload(i)%coupling=0
        else
          ! BEAM TIP -- POINT LOAD
          coupled=.true.
          be_bodyload(i)%coupling=fbem_bl_coupling_beam_tip
        end if
      !
      ! LINE LOAD
      !
      case (1)
        if (element(se)%element.eq.0) then
          if (subedge(element(se)%subedge(1))%element.eq.0) then
            coupled=.false.
            be_bodyload(i)%coupling=0
          else
            ! SHELL EDGE -- LINE LOAD
            coupled=.true.
            be_bodyload(i)%coupling=fbem_bl_coupling_shell_edge
          end if
        else
          ! BEAM -- LINE LOAD
          coupled=.true.
          be_bodyload(i)%coupling=fbem_bl_coupling_beam_line
        end if
      !
      ! SURFACE LOAD
      !
      case (2)
        if (element(se)%element.eq.0) then
          coupled=.false.
          be_bodyload(i)%coupling=0
        else
          ! SHELL SURFACE -- SURFACE LOAD
          coupled=.true.
          be_bodyload(i)%coupling=fbem_bl_coupling_shell_surface
        end if
      case default
        stop 'BE BODY LOAD - FEM: not yet 1 ...'
    end select
    ! Loop through the rest of elements of the BE body load i.
    do j=2,part(be_bodyload(i)%part)%n_elements
      se=part(be_bodyload(i)%part)%element(j)
      call build_connectivity_be_bodyload_element_fe(se)
      select case (element(se)%n_dimension)
        !
        ! POINT LOAD
        !
        case (0)
          if (node(element(se)%node(1))%coupled_node.eq.0) then
            if (coupled.eqv.(.true.)) then
              call fbem_error_message(error_unit,0,'BE body load',be_bodyload(i)%id,&
                                  'all the loads in a BE body load must be or coupled or uncoupled with a finite element')
            end if
          else
            if (coupled.eqv.(.false.)) then
              call fbem_error_message(error_unit,0,'BE body load',be_bodyload(i)%id,&
                                  'all the loads in a BE body load must be or coupled or uncoupled with a finite element')
            end if
          end if
        !
        ! LINE LOAD
        !
        case (1)
          if (element(se)%element.eq.0) then
            if (subedge(element(se)%subedge(1))%element.eq.0) then
              if (coupled.eqv.(.true.)) then
                call fbem_error_message(error_unit,0,'BE body load',be_bodyload(i)%id,&
                                    'all the loads in a BE body load must be or coupled or uncoupled with a finite element')
              end if
            else
              ! SHELL EDGE -- LINE LOAD
              if (coupled.eqv.(.false.).or.(be_bodyload(i)%coupling.ne.fbem_bl_coupling_shell_edge)) then
                call fbem_error_message(error_unit,0,'BE body load',be_bodyload(i)%id,&
                                    'all the loads in a BE body load must be or coupled or uncoupled with a finite element')
              end if
            end if
          else
            ! BEAM -- LINE LOAD
            if (coupled.eqv.(.false.).or.(be_bodyload(i)%coupling.ne.fbem_bl_coupling_beam_line)) then
              call fbem_error_message(error_unit,0,'BE body load',be_bodyload(i)%id,&
                                  'all the loads in a BE body load must be or coupled or uncoupled with a finite element')
            end if
          end if
        !
        ! SURFACE LOAD
        !
        case (2)
          if (element(se)%element.eq.0) then
            if (coupled) then
              call fbem_error_message(error_unit,0,'BE body load',be_bodyload(i)%id,&
                                  'all the loads in a BE body load must be or coupled or uncoupled with a finite element')
            end if
          else
            ! BEAM -- LINE LOAD
            if (coupled.eqv.(.false.).or.(be_bodyload(i)%coupling.ne.fbem_bl_coupling_shell_surface)) then
              call fbem_error_message(error_unit,0,'BE body load',be_bodyload(i)%id,&
                                  'all the loads in a BE body load must be or coupled or uncoupled with a finite element')
            end if
          end if
        case default
          stop 'BE BODY LOAD - FEM: not yet 1 ...'
      end select
    end do
  end do
  if (verbose_level.ge.3) then
    write(output_unit,*)
    write(output_unit,'(2x,a)') '================'
    write(output_unit,'(2x,a)') 'BE LOAD COUPLING'
    write(output_unit,'(2x,a)') '================'
    write(output_unit,*)
    do i=1,n_be_bodyloads
      select case (be_bodyload(i)%coupling)
        case (0)
          write(*,'(a,i11,a)') 'BE body load: ', be_bodyload(i)%id, ' (UNCOUPLED)'
        case (fbem_bl_coupling_beam_tip)
          write(*,'(a,2i11)') 'BE body load: ', be_bodyload(i)%id, be_bodyload(i)%part
          do j=1,part(be_bodyload(i)%part)%n_elements
            se=part(be_bodyload(i)%part)%element(j)
            write(*,'(a,i11)') 'BE body point load element: ', element(se)%id
            sn=element(se)%node(1) ! Node of the BE point load
            write(*,'(a,2i11)') 'BE point of view: ', node(sn)%id, node(node(sn)%coupled_node)%id
            sn=node(element(se)%node(1))%coupled_node ! Node of the FE beam tip
            write(*,'(a,2i11)') 'FE point of view: ', node(sn)%id, node(node(sn)%coupled_node)%id
          end do
        case (fbem_bl_coupling_beam_line,fbem_bl_coupling_shell_surface)
          write(*,'(a,2i11)') 'BE body load: ', be_bodyload(i)%id, be_bodyload(i)%part
          do j=1,part(be_bodyload(i)%part)%n_elements
            se=part(be_bodyload(i)%part)%element(j)
            if (be_bodyload(i)%coupling.eq.fbem_bl_coupling_beam_line) then
              write(*,'(a,i11)') 'BE body line load element: ', element(se)%id
            else
              write(*,'(a,i11)') 'BE body surface load element: ', element(se)%id
            end if
            write(*,'(a,8i11)') 'BE point of view: ', element(se)%id, (node(element(se)%node(k))%id,k=1,element(se)%n_nodes), &
                                                      element(element(se)%element)%id, (node(element(se)%element_node(k))%id,k=1,element(se)%n_nodes)
            se=element(se)%element ! FE element coupled with the BE line load
            write(*,'(a,8i11)') 'FE point of view: ', element(se)%id, (node(element(se)%node(k))%id,k=1,element(se)%n_nodes), &
                                                      element(element(se)%element)%id, (node(element(se)%element_node(k))%id,k=1,element(se)%n_nodes)
          end do
        case (fbem_bl_coupling_shell_edge)
          write(*,'(a,2i11)') 'BE body load: ', be_bodyload(i)%id, be_bodyload(i)%part
          do j=1,part(be_bodyload(i)%part)%n_elements
            se=part(be_bodyload(i)%part)%element(j)
            sse=element(se)%subedge(1)
            write(*,'(a,i11)') 'BE body line load element: ', element(se)%id
            write(*,'(a,8i11)') 'BE point of view: ', element(se)%id, (node(subedge(sse)%node(k))%id,k=1,subedge(sse)%n_nodes), &
                                                      element(subedge(subedge(sse)%element)%supelement(1))%id, (node(subedge(sse)%element_node(k))%id,k=1,subedge(sse)%n_nodes)
            sse=subedge(sse)%element ! FE element coupled with the BE line load
            write(*,'(a,8i11)') 'FE point of view: ', element(subedge(sse)%supelement(1))%id, (node(subedge(sse)%node(k))%id,k=1,subedge(sse)%n_nodes), &
                                                      element(subedge(subedge(sse)%element)%supelement(1))%id, (node(subedge(sse)%element_node(k))%id,k=1,subedge(sse)%n_nodes)
          end do
      end select
    end do

  end if

  ! Ending message
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END building connectivities between BE bodyloads and FE')

  ! Meter por aqui las opciones de interpolacion funcional de la carga en punta, asi como la geometria exacta (si no se quiere)
  ! la simple recta o circulo. Habra que modificar cosas para esto mas adelante.

end subroutine build_connectivities_be_bodyload_elements_and_fe
