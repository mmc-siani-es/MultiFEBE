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


subroutine export_solution_mechanics_static_sa

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

  integer                                 :: ka, kr, kb, ke, kn, kip, kc, kfn, kncbl
  integer                                 :: sb, se, sn, ks, ss, sp
  integer                                 :: k, ndof_fe_node
  integer                                 :: output_fileunit
  character(len=fbem_filename_max_length) :: tmp_filename      ! Temporary file name
  logical                                 :: node_used(n_nodes)
  character(len=fbem_fmtstr)              :: fmtstr            ! String used for write format string

  ! Starting message
  if (verbose_level.ge.1)  write(output_unit,'(a)') 'Exporting sensitivities solution ...'

  ! ================================================================================================================================
  ! EXPORT NODAL SOLUTIONS
  ! ================================================================================================================================

  if (export_nso) then
    ! Loop through design variables
    do ka=1,problem%n_designvariables

      ! --------------------
      ! Open ".sa*.nso" file
      ! --------------------

      output_fileunit=fbem_get_valid_unit()
      write(tmp_filename,*) trim(output_filename),'.sa',ka,'.nso'
      call fbem_trimall(tmp_filename)
      open(unit=output_fileunit,file=trim(tmp_filename),action='write',recl=fbem_file_record_length)
      if (verbose_level.ge.2) write(output_unit,'(a,a,a)') ' Opening "', trim(tmp_filename), '"'

      ! Loop through regions
      do kr=1,n_regions
        write(fmtstr,*) '(a,i',fbem_nchar_int(region(kr)%id),',a)'
        call fbem_trim2b(fmtstr)
        write(output_fileunit,fmtstr,advance='no') '# Region ', region(kr)%id, ', '
        select case (region(kr)%class)

          ! ============================================================================================================================
          ! BE region
          !
          case (fbem_be)
            write(output_fileunit,'(a)') 'viscoelastic BE'
            write(output_fileunit,'(a)') '# Columns: region, node, x_k, u_k, p_k, face'
            write(fmtstr,*) '(2',fmt_integer,',',3*(problem%n),fmt_real,',i3)'
            call fbem_trim2b(fmtstr)
            do kb=1,region(kr)%n_boundaries
              sb=region(kr)%boundary(kb)
              node_used=.false.
              do ke=1,part(boundary(sb)%part)%n_elements
                se=part(boundary(sb)%part)%element(ke)
                do kn=1,element(se)%n_nodes
                  sn=element(se)%node(kn)
                  if (node_used(sn).eqv.(.false.)) then
                    node_used(sn)=.true.
                    !
                    ! Boundary coupling
                    !
                    select case (boundary(sb)%coupling)
                      !
                      ! Uncoupled boundary
                      !
                      case (fbem_boundary_coupling_be)
                        !
                        ! Boundary class
                        !
                        select case (boundary(sb)%class)
                          !
                          ! Ordinary boundary
                          !
                          case (fbem_boundary_class_ordinary)
                            write(output_fileunit,fmtstr) region(kr)%id, node(sn)%id, (node(sn)%x(k),k=1,problem%n),&
                                                                                      (node(sn)%dvda_r(k,1,ka),k=1,2*problem%n), 1
                          !
                          ! Crack-like boundaries
                          !
                          case (fbem_boundary_class_cracklike)
                            write(output_fileunit,fmtstr) region(kr)%id, node(sn)%id, (node(sn)%x(k),k=1,problem%n),&
                                                                                      (node(sn)%dvda_r(k,1,ka),k=1,2*problem%n), 1
                            write(output_fileunit,fmtstr) region(kr)%id, node(sn)%id, (node(sn)%x(k),k=1,problem%n),&
                                                                                      (node(sn)%dvda_r(k,2,ka),k=1,2*problem%n), 2
                        end select
                      !
                      ! BE-BE coupled boundary
                      !
                      case (fbem_boundary_coupling_be_be)
                        if (region(kr)%boundary_reversion(kb).eqv.(.false.)) then
                          write(output_fileunit,fmtstr) region(kr)%id, node(sn)%id, (node(sn)%x(k),k=1,problem%n),&
                                                                                    (node(sn)%dvda_r(k,1,ka),k=1,2*problem%n), 1
                        else
                          write(output_fileunit,fmtstr) region(kr)%id, node(sn)%id, (node(sn)%x(k),k=1,problem%n),&
                                                                                    (node(sn)%dvda_r(k,2,ka),k=1,2*problem%n), 2
                        end if
                      !
                      ! BE-FE coupled boundary
                      !
                      case (fbem_boundary_coupling_be_fe)
                        !
                        ! Boundary class
                        !
                        select case (boundary(sb)%class)
                          !
                          ! Ordinary boundary
                          !
                          case (fbem_boundary_class_ordinary)
                            write(output_fileunit,fmtstr) region(kr)%id, node(sn)%id, (node(sn)%x(k),k=1,problem%n),&
                                                                                      (node(sn)%dvda_r(k,1,ka),k=1,2*problem%n), 1
                          !
                          ! Crack-like boundaries
                          !
                          case (fbem_boundary_class_cracklike)
                            write(output_fileunit,fmtstr) region(kr)%id, node(sn)%id, (node(sn)%x(k),k=1,problem%n),&
                                                                                      (node(sn)%dvda_r(k,1,ka),k=1,2*problem%n), 1
                            write(output_fileunit,fmtstr) region(kr)%id, node(sn)%id, (node(sn)%x(k),k=1,problem%n),&
                                                                                      (node(sn)%dvda_r(k,2,ka),k=1,2*problem%n), 2
                        end select
                      !
                      ! BE-FE-BE coupled boundary
                      !
                      case (fbem_boundary_coupling_be_fe_be)
                        if (region(kr)%boundary_reversion(kb).eqv.(.false.)) then
                          write(output_fileunit,fmtstr) region(kr)%id, node(sn)%id, (node(sn)%x(k),k=1,problem%n),&
                                                                                    (node(sn)%dvda_r(k,1,ka),k=1,2*problem%n), 1
                        else
                          write(output_fileunit,fmtstr) region(kr)%id, node(sn)%id, (node(sn)%x(k),k=1,problem%n),&
                                                                                    (node(sn)%dvda_r(k,2,ka),k=1,2*problem%n), 2
                        end if
                    end select
                  end if
                end do
              end do
            end do
          ! ============================================================================================================================

          ! ============================================================================================================================
          ! FE region
          !
          case (fbem_fe)
            write(output_fileunit,'(a)') 'viscoelastic FE'
            write(output_fileunit,'(a)') '# Columns: region, node, x_k, u_k[, theta_k]'
            node_used=.false.
            do ks=1,region(kr)%n_fe_subregions
              ss=region(kr)%fe_subregion(ks)
              do ke=1,part(fe_subregion(ss)%part)%n_elements
                se=part(fe_subregion(ss)%part)%element(ke)
                do kn=1,element(se)%n_nodes
                  sn=element(se)%node(kn)
                  if (node_used(sn).eqv.(.false.)) then
                    node_used(sn)=.true.
                    ndof_fe_node=node(sn)%n_dof
                    write(fmtstr,*) '(2',fmt_integer,',',problem%n+node(sn)%n_dof,fmt_real,')'
                    call fbem_trim2b(fmtstr)
                    write(output_fileunit,fmtstr) region(kr)%id, node(sn)%id, (node(sn)%x(k),k=1,problem%n),&
                                                                              (node(sn)%dvda_r(k,1,ka),k=1,node(sn)%n_dof)
                  end if
                end do
              end do
            end do
          ! ============================================================================================================================

        end select
        write(output_fileunit,*)
        write(output_fileunit,*)
      end do

      ! ----------
      ! Close file
      ! ----------

      close(unit=output_fileunit)
      if (verbose_level.ge.2) write(output_unit,'(a,a,a)') ' Closing "', trim(tmp_filename), '"'

    end do
  end if

  ! Ending message
  if (verbose_level.ge.1) write(output_unit,'(a)') 'done.'

end subroutine export_solution_mechanics_static_sa
