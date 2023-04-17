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


subroutine export_solution_mechanics_harmonic_sa(kf)

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

  ! I/O variables
  integer                                 :: kf
  ! Local variables
  real(kind=real64)                       :: omega
  integer                                 :: kr, kb, ke, kn, kc, ka
  integer                                 :: sb, se, sn, ks, ss
  integer                                 :: k
  integer                                 :: output_fileunit
  character(len=fbem_filename_max_length) :: tmp_filename
  logical                                 :: node_used(n_nodes)
  character(len=fbem_fmtstr)              :: fmt1, fmt2
  integer                                 :: ncint, ncreal, nc, ncmax
  character(len=fbem_string_max_length)   :: tmp_string
  integer                                 :: face, k_start, k_end

  ! Starting message
  if (verbose_level.ge.1)  write(output_unit,'(a)') 'Exporting sensitivities solution ...'

  ! Frequency
  omega=frequency(kf)
  if (frequency_units.eq.'f') omega=omega*c_1_2pi

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
      if (kf.eq.1) then
        open(unit=output_fileunit,file=trim(tmp_filename),action='write',recl=fbem_file_record_length)
      else
        open(unit=output_fileunit,file=trim(tmp_filename),access='append',recl=fbem_file_record_length)
      end if
      if (verbose_level.ge.2) write(output_unit,'(a,a,a)') ' Opening "', trim(tmp_filename), '"'

      ! ------------------
      ! Header of the file
      ! ------------------

      if (kf.eq.1) then
        write(output_fileunit,'(a)'  ) '# Program      : multifebe'
        write(output_fileunit,'(a,a5)')'# Version      : ', multifebe_version
        write(output_fileunit,'(a)'  ) '# File_format  : sa*.nso'
        write(output_fileunit,'(a,i1)')'# Problem_dim  : ', problem%n
        write(output_fileunit,'(a,a)') '# Input_file   : ', trim(input_filename)
        write(output_fileunit,'(a,a)') '# Description  : ', trim(problem%description)
        write(output_fileunit,'(a)',advance='no') '# Timestamp    : '
        call fbem_timestamp(output_fileunit,timestamp_date_start,timestamp_time_start,1)
        write(output_fileunit,*)
        ! Column description
        write(output_fileunit,'(a)'  ) '# Columns  Description'
        write(output_fileunit,'(a)'  ) '# C1-C2    Frequency index and value.'
        write(output_fileunit,'(a)'  ) '# C3-C5    Region id, class and type.'
        write(output_fileunit,'(a)'  ) '# C6-C8    Boundary id, class and face.'
        select case (problem%n)
        case (2)
        write(output_fileunit,'(a)'  ) '# C9-C11   Node id, x1 and x2.'
        write(output_fileunit,'(a)'  ) '# >=C12    Node variables. Depend on the region class and type.'

        case (3)
        write(output_fileunit,'(a)'  ) '# C9-C12   Node id, x1, x2 and x3.'
        write(output_fileunit,'(a)'  ) '# >=C13    Node variables. Depend on the region class and type.'
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
            ncmax=35
          case (3)
            ncmax=44
        end select
        do kc=2,ncmax
          ! Depending if integer or real column
          if ((kc.ge.3).and.(kc.le.9)) then
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
      end if
      !
      ! Each variable is expressed in cartesian notation (real,imaginary) or in polar notation (absolute value,argument), depending
      ! on the export settings. Each line contains different variables at different columns depending on the problem dimension and
      ! on the type of region (fluid, elastic, poroelastic). The variables are divided into two groups: total field and incident field
      ! (*_i) variables. Within each group, the variables are divided into two subgroups: primary and secondary variables.
      !
      ! Notation:
      !
      ! Fluid:
      !   - p : pressure
      !   - Un: normal displacement
      ! Elastic:
      !   - uk: displacement (k component)
      !   - tk: traction (k component)
      ! Poroelastic:
      !   - tau: fluid stress
      !   - uk : solid displacement (k component)
      !   - Un : fluid normal displacement
      !   - tk : solid traction (k component)
      !
      !
      ! 2D PROBLEM
      ! ----------
      !
      ! BE region:
      !
      ! |-----------------|     |-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|------|-----|
      ! | 5               |     |12,13|14,15|16,17|18,19|20,21|22,23|24,25|26,27|28,29|30,31|32,33 |34,35|
      ! |-----------------| ... |-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|------|-----|
      ! | 1 (fluid)       |     |    p|   Un|  p_i| Un_i|     |     |     |     |     |     |      |     |
      ! | 2 (elastic)     |     |   u1|   u2|   t1|   t2| u1_i| u2_i| t1_i| t2_i|     |     |      |     |
      ! | 3 (poroelastic) |     |  tau|   u1|   u2|   Un|   t1|   t2|tau_i| u1_i| u2_i| Un_i| t1_i | t2_i|
      ! |-----------------|     |-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|------|-----|
      !
      ! 3D PROBLEM
      ! ----------
      !
      ! BE region:
      !
      ! |-----------------|     |-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
      ! | 5               |     |13,14|15,16|17,18|19,20|21,22|23,24|25,26|27,28|29,30|31,32|33,34|35,36|37,38|39,40|41,42|43,44|
      ! |-----------------| ... |-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
      ! | 1 (fluid)       |     |    p|   Un|  p_i| Un_i|     |     |     |     |     |     |     |     |     |     |     |     |
      ! | 2 (elastic)     |     |   u1|   u2|   u3|   t1|   t2|   t3| u1_i| u2_i| u3_i| t1_i| t2_i| t3_i|     |     |     |     |
      ! | 3 (poroelastic) |     |  tau|   u1|   u2|   u3|   Un|   t1|   t2|   t3|tau_i| u1_i| u2_i| u3_i| Un_i| t1_i| t2_i| t3_i|
      ! |-----------------|     |-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|

      ! --------------
      ! Export formats
      ! --------------

      ! Export format for columns 1-11 (2D) or 1-12 (3D)
      write(fmt1,*) '(1',fmt_integer,',1',fmt_real,',7',fmt_integer,',',problem%n,fmt_real,')'
      call fbem_trim2b(fmt1)
      ! Export format for columns >=12 (2D) or >=13 (3D)
      write(fmt2,*) '(2',fmt_real,')'
      call fbem_trim2b(fmt2)

      ! -----
      ! Write
      ! -----

      ! Loop through REGIONS
      do kr=1,n_regions
        select case (region(kr)%class)

          ! ==========================================================================================================================
          ! BE REGION
          !
          case (fbem_be)

            ! +----------------+
            ! | BOUNDARY NODES |
            ! +----------------+

            ! Select the start and end indices of variables
            select case (region(kr)%type)
              ! INVISCID FLUID
              case (fbem_potential)
                k_start=1
                k_end  =2
              ! VISCOELASTIC SOLID
              case (fbem_viscoelastic)
                k_start=1
                k_end  =2*problem%n
              ! POROELASTIC MEDIA
              case (fbem_poroelastic)
                k_start=0
                k_end  =1+2*problem%n
            end select
            ! Loop through boundaries
            do kb=1,region(kr)%n_boundaries
              sb=region(kr)%boundary(kb)
              ! Initialize
              node_used=.false.
              ! Loop through elements and nodes
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
                            ! Write columns 1-11 (2D) or 1-12 (3D)
                            write(output_fileunit,fmt1,advance='no') kf, omega, &
                                                                     region(kr)%id, region(kr)%class, region(kr)%type, &
                                                                     boundary(sb)%id, boundary(sb)%class, face, &
                                                                     node(sn)%id, element(se)%x_fn(:,kn)
                            ! Write columns >=12 (2D) or >=13 (3D)
                            ! TOTAL FIELD
                            do k=k_start,k_end
                              select case (complex_notation)
                                case (1)
                                  write(output_fileunit,fmt2,advance='no') abs(node(sn)%dvda_c(k,face,ka)), fbem_zarg(node(sn)%dvda_c(k,face,ka))
                                case (2)
                                  write(output_fileunit,fmt2,advance='no') real(node(sn)%dvda_c(k,face,ka)), imag(node(sn)%dvda_c(k,face,ka))
                              end select
                            end do
                            write(output_fileunit,*)
                          !
                          ! Crack-like boundaries
                          !
                          case (fbem_boundary_class_cracklike)
                            ! Face + and face -
                            do face=1,2
                              ! Write columns 1-11 (2D) or 1-12 (3D)
                              write(output_fileunit,fmt1,advance='no') kf, omega, &
                                                                       region(kr)%id, region(kr)%class, region(kr)%type, &
                                                                       boundary(sb)%id, boundary(sb)%class, face, &
                                                                       node(sn)%id, element(se)%x_fn(:,kn)
                              ! Write columns >=12 (2D) or >=13 (3D)
                              ! TOTAL FIELD
                              do k=k_start,k_end
                                select case (complex_notation)
                                  case (1)
                                    write(output_fileunit,fmt2,advance='no') abs(node(sn)%dvda_c(k,face,ka)), fbem_zarg(node(sn)%dvda_c(k,face,ka))
                                  case (2)
                                    write(output_fileunit,fmt2,advance='no') real(node(sn)%dvda_c(k,face,ka)), imag(node(sn)%dvda_c(k,face,ka))
                                end select
                              end do
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
                        ! Write columns 1-11 (2D) or 1-12 (3D)
                        write(output_fileunit,fmt1,advance='no') kf, omega, &
                                                                 region(kr)%id, region(kr)%class, region(kr)%type, &
                                                                 boundary(sb)%id, boundary(sb)%class, face, &
                                                                 node(sn)%id, element(se)%x_fn(:,kn)
                        ! Write columns >=12 (2D) or >=13 (3D)
                        ! TOTAL FIELD
                        do k=k_start,k_end
                          select case (complex_notation)
                            case (1)
                              write(output_fileunit,fmt2,advance='no') abs(node(sn)%dvda_c(k,face,ka)), fbem_zarg(node(sn)%dvda_c(k,face,ka))
                            case (2)
                              write(output_fileunit,fmt2,advance='no') real(node(sn)%dvda_c(k,face,ka)), imag(node(sn)%dvda_c(k,face,ka))
                          end select
                        end do
                        write(output_fileunit,*)
                    end select
                  end if
                end do
              end do
            end do
          ! ==========================================================================================================================

          ! ==========================================================================================================================
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
                  if ((node_used(sn).eqv.(.false.)).and.(node(sn)%export.eqv.(.true.))) then
                    node_used(sn)=.true.
                    ! Write columns 1-11 (2D) or 1-12 (3D)
                    write(output_fileunit,fmt1,advance='no') kf, omega, region(kr)%id, region(kr)%class, region(kr)%type, 0, 0, 0, node(sn)%id, node(sn)%x
                    ! Write columns >=12 (2D) or >=13 (3D)
                    do k=1,node(sn)%n_dof
                      select case (complex_notation)
                        case (1)
                          write(output_fileunit,fmt2,advance='no') abs(node(sn)%dvda_c(k,1,ka)), fbem_zarg(node(sn)%dvda_c(k,1,ka))
                        case (2)
                          write(output_fileunit,fmt2,advance='no') real(node(sn)%dvda_c(k,1,ka)), imag(node(sn)%dvda_c(k,1,ka))
                      end select
                    end do
                    write(output_fileunit,*)
                  end if
                end do
              end do
            end do
          ! ==========================================================================================================================

        end select

      end do ! Loop through REGIONS

      ! ----------
      ! Close file
      ! ----------

      close(unit=output_fileunit)
      if (verbose_level.ge.2) write(output_unit,'(a,a,a)') ' Closing "', trim(tmp_filename), '"'

    end do

  end if

  ! Ending message
  if (verbose_level.ge.1) write(output_unit,'(a)') 'done.'

end subroutine export_solution_mechanics_harmonic_sa
