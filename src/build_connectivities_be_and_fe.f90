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

! BE WITH SOLID (N)D FE
!   2D PROBLEM: SOLID FE (2D) - ORDINARY BE (1D)
!   3D PROBLEM: SOLID FE (3D) - ORDINARY BE (2D)
! BE WITH STRUCTURAL (N-1)D FE
!   2D PROBLEM: BEAM FE (1D)  - ORDINARY OR CRACK-LIKE BE (1D)
!   3D PROBLEM: SHELL FE (2D) - ORDINARY OR CRACK-LIKE BE (2D)
subroutine build_connectivities_be_and_fe

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
  integer                               :: se
  logical                               :: coupled
  integer                               :: se_be, se_fe, sn_be, sn_fe, kd_fe, sd_fe
  character(len=fbem_fmtstr)            :: fmtstr

  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'START building connectivities between BE and FE')

  !
  ! Loop through boundaries
  !
  do i=1,n_boundaries
    !
    ! Build BE-FE connectivities
    !
    ! First element of the boundary i
    se=part(boundary(i)%part)%element(1)
    call build_connectivity_be_fe(se)
    ! True if coupled
    if (element(se)%element.eq.0) then
      coupled=.false.
    else
      coupled=.true.
    end if
    ! Rest of elements of the boundary i
    do j=2,part(boundary(i)%part)%n_elements
      se=part(boundary(i)%part)%element(j)
      call build_connectivity_be_fe(se)
      ! Each element must have the same coupling of the first element.
      if (element(se)%element.eq.0) then
        if (coupled) then
          call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'all its BEs must be coupled or uncoupled')
        end if
      else
        if (coupled.eqv.(.false.)) then
          call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,'all its BEs must be coupled or uncoupled')
        end if
      end if
    end do
    !
    ! Assign the class of boundary and check if FE-BE are compatible
    !
    select case (boundary(i)%class)
      case (fbem_boundary_class_ordinary)
        select case (boundary(i)%n_regions)
          ! EXTERIOR BOUNDARY
          case (1)
            if (coupled) then
              boundary(i)%coupling=fbem_boundary_coupling_be_fe
              ! (N-1)D BE - (N-1)D FE
              se=part(boundary(i)%part)%element(1)
              if (element(element(se)%element)%n_dimension.eq.(problem%n-1)) then
                do j=2,part(boundary(i)%part)%n_elements
                  se=part(boundary(i)%part)%element(j)
                  if (element(element(se)%element)%n_dimension.ne.(problem%n-1)) then
                    call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,&
                                            'all FE connected with this boundary must have the same dimension')
                  end if
                end do
              end if
              ! (N-1)D BE - ND FE
              se=part(boundary(i)%part)%element(1)
              if (element(element(se)%element)%n_dimension.eq.problem%n) then
                do j=2,part(boundary(i)%part)%n_elements
                  se=part(boundary(i)%part)%element(j)
                  if (element(element(se)%element)%n_dimension.ne.problem%n) then
                    call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,&
                                            'all FE connected with this boundary must have the same dimension')
                  end if
                end do
              end if
            else
              boundary(i)%coupling=fbem_boundary_coupling_be
            end if
          ! INTERFACE BOUNDARY
          case (2)
            if (coupled) then
              boundary(i)%coupling=fbem_boundary_coupling_be_fe_be
              do j=1,part(boundary(i)%part)%n_elements
                se=part(boundary(i)%part)%element(j)
                if (element(element(se)%element)%n_dimension.ne.(problem%n-1)) then
                  call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,&
                                          'all FE connected with this boundary must have the same dimension')
                end if
              end do
            else
              boundary(i)%coupling=fbem_boundary_coupling_be_be
            end if
        end select
      ! CRACK-LIKE BOUNDARY
      case (fbem_boundary_class_cracklike)
        if (coupled) then
          boundary(i)%coupling=fbem_boundary_coupling_be_fe
          do j=1,part(boundary(i)%part)%n_elements
            se=part(boundary(i)%part)%element(j)
            if (element(element(se)%element)%n_dimension.ne.(problem%n-1)) then
              call fbem_error_message(error_unit,0,'boundary',boundary(i)%id,&
                                      'all FE connected with this boundary must have the same dimension')
            end if
          end do
        else
          boundary(i)%coupling=fbem_boundary_coupling_be
        end if
    end select
  end do
  !
  ! Write results
  !
  if (verbose_level.ge.3) then
    write(output_unit,*)
    write(output_unit,'(2x,a)') '=============='
    write(output_unit,'(2x,a)') 'BOUNDARY CLASS'
    write(output_unit,'(2x,a)') '=============='
    write(output_unit,*)
    write(output_unit,'(2x,a11,1x,a8,1x,a11,1x,a11,1x,a6)') '___boundary','___class','BE_region_1','BE_region_2','FE_dim'
    do i=1,n_boundaries
      se=part(boundary(i)%part)%element(1)
      select case (boundary(i)%coupling)
        case (fbem_boundary_coupling_be)
          write(output_unit,'(2x,i11,1x,a8)',advance='no') boundary(i)%id, '      BE'
        case (fbem_boundary_coupling_be_fe)
          write(output_unit,'(2x,i11,1x,a8)',advance='no') boundary(i)%id, '   BE-FE'
        case (fbem_boundary_coupling_be_be)
          write(output_unit,'(2x,i11,1x,a8)',advance='no') boundary(i)%id, '   BE-BE'
        case (fbem_boundary_coupling_be_fe_be)
          write(output_unit,'(2x,i11,1x,a8)',advance='no') boundary(i)%id, 'BE-FE-BE'
      end select
      select case (boundary(i)%n_regions)
        case (1)
          write(output_unit,'(1x,i11,1x,11x)',advance='no') region(boundary(i)%region(1))%id
        case (2)
          write(output_unit,'(1x,i11,1x,i11)',advance='no') region(boundary(i)%region(1))%id, region(boundary(i)%region(2))%id
      end select
      if (element(se)%element.ne.0) then
        write(output_unit,'(1x,i6)') element(element(se)%element)%n_dimension
      else
        write(output_unit,*)
      end if
    end do
    write(output_unit,*)
    write(output_unit,'(2x,a)') '=================='
    write(output_unit,'(2x,a)') 'BE-FE CONNECTIVITY'
    write(output_unit,'(2x,a)') '=================='
    write(output_unit,*)
    do i=1,n_boundaries
      if ((boundary(i)%coupling.eq.fbem_boundary_coupling_be_fe).or.(boundary(i)%coupling.eq.fbem_boundary_coupling_be_fe_be)) then
        do j=1,part(boundary(i)%part)%n_elements
          se_be=part(boundary(i)%part)%element(j)
          se_fe=element(se_be)%element
          !
          ! BE->FE
          !
          write(fmtstr,*) '(2x,a3,i',fbem_nchar_int(element(se_be)%id),',a7,i',fbem_nchar_int(element(se_fe)%id),',a20)'
          call fbem_trimall(fmtstr)
          write(output_unit,fmtstr,advance='no') 'BE ',element(se_be)%id,' -> FE ',element(se_fe)%id,'; BE_node->FE_node={'
          do k=1,element(se_be)%n_nodes
            sn_be=element(se_be)%node(k)
            sn_fe=element(se_be)%element_node(k)
            write(fmtstr,*) '(i',fbem_nchar_int(node(sn_be)%id),',a2,i',fbem_nchar_int(node(sn_fe)%id),',a1)'
            call fbem_trimall(fmtstr)
            if (k.lt.element(se_be)%n_nodes) then
              write(output_unit,fmtstr,advance='no') node(sn_be)%id,'->',node(sn_fe)%id,','
            else
              write(output_unit,fmtstr) node(sn_be)%id,'->',node(sn_fe)%id,'}'
            end if
          end do
          !
          ! 2D SOLID -> BE
          !
          ! FE->BE
          if ((problem%n.eq.2).and.(element(se_fe)%n_dimension.eq.2)) then
            do kd_fe=1,element(se_fe)%n_subedges
              sd_fe=element(se_fe)%subedge(kd_fe)
              if (subedge(sd_fe)%element.ne.0) then
                se_be=subedge(sd_fe)%element
                write(fmtstr,*) '(2x,a3,i',fbem_nchar_int(element(se_fe)%id),',a11,i',fbem_nchar_int(kd_fe),',a7,i',fbem_nchar_int(element(se_be)%id),',a20)'
                call fbem_trimall(fmtstr)
                write(output_unit,fmtstr,advance='no') 'FE ',element(se_fe)%id,' / subedge ',kd_fe,' -> BE ',element(se_be)%id,'; BE_node->FE_node={'
                do k=1,subedge(sd_fe)%n_nodes
                  sn_be=subedge(sd_fe)%node(k)
                  sn_fe=subedge(sd_fe)%element_node(k)
                  write(fmtstr,*) '(i',fbem_nchar_int(node(sn_be)%id),',a2,i',fbem_nchar_int(node(sn_fe)%id),',a1)'
                  call fbem_trimall(fmtstr)
                  if (k.lt.subedge(sd_fe)%n_nodes) then
                    write(output_unit,fmtstr,advance='no') node(sn_be)%id,'->',node(sn_fe)%id,','
                  else
                    write(output_unit,fmtstr) node(sn_be)%id,'->',node(sn_fe)%id,'}'
                  end if
                end do
                exit
              end if
            end do
          end if



          !
          ! Falta... para 3D etc....
          !



        end do
      end if
    end do
  end if



  ! Chequear las normales de los BE con respecto a solid FE, para que sea coherente...


  ! Ending message
  if (verbose_level.ge.2) call fbem_timestamp_w_message(output_unit,2,'END building connectivities between BE and FE')

end subroutine build_connectivities_be_and_fe
