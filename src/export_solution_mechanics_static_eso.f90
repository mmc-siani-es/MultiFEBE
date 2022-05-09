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

subroutine export_solution_mechanics_static_eso(output_fileunit)

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

  !
  ! Falta meter los ejes de los esfuerzos internos
  !

  if (export_overwrite) then
    write(output_fileunit,'(a)'  ) '# Program      : multifebe'
    write(output_fileunit,'(a)'  ) '# Version      : 1.0'
    write(output_fileunit,'(a)'  ) '# File_format  : eso'
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
    write(output_fileunit,'(a)'  ) '# 6-8      (if col 4 == 2) Subregion id, finite element dimension and finite element type.'
    ! Falta if col4==1
    write(output_fileunit,'(a)'  ) '# 9-11     (if col 4 == 2) Element id, element dimension, element type.'
    write(output_fileunit,'(a)'  ) '# 12-13    (if col 4 == 2) element node index, element node number of DOFs.'

    select case (problem%n)
    case (2)
    write(output_fileunit,'(a)'  ) '# 14-16    Node id, x1 and x2.'
    write(output_fileunit,'(a)'  ) '# >=17     Node variables (see documentation).'
    case (3)
    write(output_fileunit,'(a)'  ) '# 14-17    Node id, x1, x2 and x3.'
    write(output_fileunit,'(a)'  ) '# >=18     Node variables (see documentation).'
    end select
    write(output_fileunit,'(a)') '#'
  end if
  ! --------------
  ! Export formats
  ! --------------

  ! Export format for columns 1-16 (2D) or 1-17 (3D)
  write(fmt1,*) '(1',fmt_integer,',1',fmt_real,',12',fmt_integer,',',problem%n,fmt_real,')'
  call fbem_trim2b(fmt1)
  ! Export format for columns >=17 (2D) or >=18 (3D)
  write(fmt2,*) '(',fmt_real,')'
  call fbem_trim2b(fmt2)

  ! -----
  ! Write
  ! -----

  ! Loop through regions
  do kr=1,n_regions
    select case (region(kr)%class)

      ! ==========================================================================================================================
      ! BE region
      !
      case (fbem_be)

      ! ==========================================================================================================================
      ! FE region
      !
      case (fbem_fe)
        do ks=1,region(kr)%n_fe_subregions
          ss=region(kr)%fe_subregion(ks)
          do ke=1,part(fe_subregion(ss)%part)%n_elements
            se=part(fe_subregion(ss)%part)%element(ke)
            select case (element(se)%n_dimension)

              ! ==================================================================================================================
              ! ONE-DIMENSIONAL ELEMENTS
              ! ==================================================================================================================

              case (1)

                select case (element(se)%fe_type)

                  case (0,1,2)

                    !-------------------------------------------------------------------------------------------------------------
                    ! DEGENERATED BEAM AND STRAIGHT EULER-BERNOULLI AND TIMOSHENKO BEAM FINITE ELEMENTS
                    !
                    do kn=1,element(se)%n_nodes
                      sn=element(se)%node(kn)
                      write(output_fileunit,fmt1,advance='no') 0, 0.d0, &
                                                               region(kr)%id, region(kr)%class, region(kr)%type, &
                                                               fe_subregion(ss)%id, element(se)%n_dimension, element(se)%fe_type, &
                                                               element(se)%id, element(se)%n_dimension, element(se)%fe_type,&
                                                               kn, element(se)%node_n_dof(kn), node(sn)%id, node(sn)%x
                      !
                      ! Stress resultants (local axes v_k/ep_k): value_r(k,n,1): NX,VY,BZ OR NX,VY,VZ,BX,BY,BZ
                      ! Equilibrating loads (global axes)      : value_r(k,n,2): FX,FY[,MZ] OR FX,FY,FZ[,MX,MY,MZ]
                      !
                      do k=1,3*(problem%n-1)
                        write(output_fileunit,fmt2,advance='no') element(se)%value_r(k,kn,1)
                      end do
!~                       do k=1,element(se)%node_n_dof(kn)
!~                         write(output_fileunit,fmt2,advance='no') element(se)%value_r(k,kn,2)
!~                       end do
                      write(output_fileunit,*)
                    end do
                    !
                    !-------------------------------------------------------------------------------------------------------------

                  case (3)

                    !-------------------------------------------------------------------------------------------------------------
                    ! BAR FINITE ELEMENTS
                    !
                    do kn=1,element(se)%n_nodes
                      sn=element(se)%node(kn)
                      write(output_fileunit,fmt1,advance='no') 0, 0.d0, &
                                                               region(kr)%id, region(kr)%class, region(kr)%type, &
                                                               fe_subregion(ss)%id, element(se)%n_dimension, element(se)%fe_type, &
                                                               element(se)%id, element(se)%n_dimension, element(se)%fe_type,&
                                                               kn, element(se)%node_n_dof(kn), node(sn)%id, node(sn)%x
                      !
                      ! Stress resultants (local axes v_k/ep_k): value_r(k,n,1): NX
                      ! Equilibrating loads (global axes)      : value_r(k,n,2): FX,FY OR FX,FY,FZ
                      !
                      write(output_fileunit,fmt2,advance='no') element(se)%value_r(1,kn,1)
!~                       do k=1,element(se)%node_n_dof(kn)
!~                         write(output_fileunit,fmt2,advance='no') element(se)%value_r(k,kn,2)
!~                       end do
                      write(output_fileunit,*)
                    end do
                    !
                    !-------------------------------------------------------------------------------------------------------------

                  case (4)

                    !-------------------------------------------------------------------------------------------------------------
                    ! DISCRETE TRANSLATIONAL SPRING FINITE ELEMENTS
                    !
                    do kn=1,element(se)%n_nodes
                      sn=element(se)%node(kn)
                      write(output_fileunit,fmt1,advance='no') 0, 0.d0, &
                                                               region(kr)%id, region(kr)%class, region(kr)%type, &
                                                               fe_subregion(ss)%id, element(se)%n_dimension, element(se)%fe_type, &
                                                               element(se)%id, element(se)%n_dimension, element(se)%fe_type,&
                                                               kn, element(se)%node_n_dof(kn), node(sn)%id, node(sn)%x
                      !
                      ! Spring forces                     : value_r(k,n,1): NX,NY,NZ
                      ! Equilibrating loads (nodal axes)  : value_r(k,n,2): FX,FY,FZ
                      !
                      do k=1,problem%n
                        write(output_fileunit,fmt2,advance='no') element(se)%value_r(k,kn,1)
                      end do
!~                       do k=1,problem%n
!~                         write(output_fileunit,fmt2,advance='no') element(se)%value_r(k,kn,2)
!~                       end do
                      write(output_fileunit,*)
                    end do
                    !
                    !-------------------------------------------------------------------------------------------------------------

                  case (5)

                    !-------------------------------------------------------------------------------------------------------------
                    ! DISCRETE ROTATIONAL/TRANSLATIONAL SPRING FINITE ELEMENTS
                    !
                    do kn=1,element(se)%n_nodes
                      sn=element(se)%node(kn)
                      write(output_fileunit,fmt1,advance='no') 0, 0.d0, &
                                                               region(kr)%id, region(kr)%class, region(kr)%type, &
                                                               fe_subregion(ss)%id, element(se)%n_dimension, element(se)%fe_type, &
                                                               element(se)%id, element(se)%n_dimension, element(se)%fe_type,&
                                                               kn, element(se)%node_n_dof(kn), node(sn)%id, node(sn)%x
                      !
                      ! Spring forces/moments             : value_r(k,n,1): NX,NY,NZ,BX,BY,BZ
                      ! Equilibrating loads (nodal axes)  : value_r(k,n,2): FX,FY,FZ,MX,MY,MZ
                      !
                      do k=1,3*(problem%n-1)
                        write(output_fileunit,fmt2,advance='no') element(se)%value_r(k,kn,1)
                      end do
!~                       do k=1,3*(problem%n-1)
!~                         write(output_fileunit,fmt2,advance='no') element(se)%value_r(k,kn,2)
!~                       end do
                      write(output_fileunit,*)
                    end do
                    !
                    !-------------------------------------------------------------------------------------------------------------

                  case default

                    !-------------------------------------------------------------------------------------------------------------
                    ! OTHER TYPES
                    !
                    call fbem_error_message(error_unit,0,'element',element(se)%id,'invalid type of 1D element')
                    !
                    !-------------------------------------------------------------------------------------------------------------

                end select

              ! ==================================================================================================================
              ! TWO-DIMENSIONAL ELEMENTS
              ! ==================================================================================================================

              case (2)

                select case (problem%n)

                  case (2)

                    !-------------------------------------------------------------------------------------------------------------
                    ! SOLID / CONTINUUM ELEMENTS
                    !

                    !
                    !
                    !-------------------------------------------------------------------------------------------------------------

                  case (3)

                    !-------------------------------------------------------------------------------------------------------------
                    ! DEGENERATED SHELL FINITE ELEMENT
                    !
                    do kn=1,element(se)%n_nodes
                      sn=element(se)%node(kn)
                      write(output_fileunit,fmt1,advance='no') 0, 0.d0, &
                                                               region(kr)%id, region(kr)%class, region(kr)%type, &
                                                               fe_subregion(ss)%id, element(se)%n_dimension, element(se)%fe_type, &
                                                               element(se)%id, element(se)%n_dimension, element(se)%fe_type,&
                                                               kn, element(se)%node_n_dof(kn), node(sn)%id, node(sn)%x
                      !
                      ! Stress resultants (local axes ep_k)               : value_r(k,n,1): NX,NY,NXY,MX,MY,MXY,VX,VY
                      ! Equilibrating loads (global u, global/local theta): value_r(k,n,2): FX,FY,FZ[,MA,MB][,MX,MY,MZ]
                      !
                      do k=1,8
                        write(output_fileunit,fmt2,advance='no') element(se)%value_r(k,kn,1)
                      end do
!~                       do k=1,6
!~                         write(output_fileunit,fmt2,advance='no') element(se)%value_r(k,kn,2)
!~                       end do
                      write(output_fileunit,*)
                    end do
                    !
                    !-------------------------------------------------------------------------------------------------------------

                end select

              ! ==================================================================================================================
              ! THREE-DIMENSIONAL ELEMENTS
              ! ==================================================================================================================

              case (3)

                !-----------------------------------------------------------------------------------------------------------------
                ! SOLID / CONTINUUM ELEMENTS
                !
                call fbem_error_message(error_unit,0,'element',element(se)%id,'3D elements not available yet')
                !
                !-----------------------------------------------------------------------------------------------------------------

            end select

          end do
        end do
      ! ==========================================================================================================================

    end select
    write(output_fileunit,*)
    write(output_fileunit,*)
  end do

end subroutine export_solution_mechanics_static_eso
