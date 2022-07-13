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

subroutine export_solution_mechanics_static_tot(output_fileunit)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_quad_rules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_shape_functions
  use fbem_geometry
  use fbem_symmetry

  ! Module of problem variables
  use problem_variables

  ! No implicit variables allowed
  implicit none

  ! I/O variables
  integer                                 :: output_fileunit
  ! Local variables
  integer                                 :: kr, kb, ke, kn, kip, kc, kfn, kncbl
  integer                                 :: sb, se, sn, ks, ss, sp
  integer                                 :: k, ndof_fe_node
  logical                                 :: node_used(n_nodes)
  character(len=fbem_fmtstr)              :: fmt1, fmt2, fmt3, fmtstr    ! String used for write format string
  integer                                 :: ncint, ncreal, nc, ncmax
  character(len=fbem_string_max_length)   :: tmp_string
  integer                                 :: face, kface, nfaces
  ! RESULTANTS
  integer                                 :: j, j1, j2, rule
  real(kind=real64), allocatable          :: xi(:), pphi(:), sphi(:)
  real(kind=real64)                       :: x(problem%n), n(problem%n), xm(problem%n), xmt(problem%n)
  real(kind=real64)                       :: w, jg, r(problem%n), r1, r2, r3
  real(kind=real64)                       :: A
  real(kind=real64)                       :: u(problem%n), t(problem%n)
  real(kind=real64)                       :: D (problem%n,2), G (2*problem%n-3,2) ! Kinematic resultants (displacements and rotation)
  real(kind=real64)                       :: F (problem%n,2), M (2*problem%n-3,2) ! Stress resultants (forces and moments)
  real(kind=real64)                       :: DS(problem%n,2), GS(2*problem%n-3,2) ! Kinematic resultants (displacements and rotation)
  real(kind=real64)                       :: FS(problem%n,2), MS(2*problem%n-3,2) ! Stress resultants (forces and moments)
  real(kind=real64)                       :: DT(problem%n,2), GT(2*problem%n-3,2) ! Kinematic resultants (displacements and rotation)
  real(kind=real64)                       :: FT(problem%n,2), MT(2*problem%n-3,2) ! Stress resultants (forces and moments)
  logical                                 :: sb_reversion, reversed
  real(kind=real64)                       :: symconf_m(3), symconf_s, symconf_t(3), symconf_r(3)

  ! ------------------
  ! Header of the file
  ! ------------------

  if (export_overwrite) then
    ! Info about the program that generated the file and its general description
    write(output_fileunit,'(a)'  ) '# Program      : multifebe'
    write(output_fileunit,'(a)'  ) '# Version      : 0.6'
    write(output_fileunit,'(a)'  ) '# File_format  : tot'
    write(output_fileunit,'(a)'  ) '# Specification: 1'
    write(output_fileunit,'(a,a)') '# Input_file   : ', trim(input_filename)
    write(output_fileunit,'(a,a)') '# Description  : ', trim(problem%description)
    write(tmp_string,*) timestamp_date_start(1:4),'-',timestamp_date_start(5:6),'-',timestamp_date_start(7:8),' ',&
                        timestamp_time_start(1:2),':',timestamp_time_start(3:4),':',timestamp_time_start(5:10)
    call fbem_trim2b(tmp_string)
    write(output_fileunit,'(a,a)') '# Timestamp    : ', trim(tmp_string)
    write(output_fileunit,*)
    ! Column description
    write(output_fileunit,'(a)'  ) '# Columns  Description'
    write(output_fileunit,'(a)'  ) '# 1-3      Region id, class and type.'
    write(output_fileunit,'(a)'  ) '# 4-6      Boundary id, class and face.'
    select case (problem%n)
    case (2)
    write(output_fileunit,'(a)'  ) '# 7-8      Reference point'
    write(output_fileunit,'(a)'  ) '# 9-10     Resultant displacements along x1 and x2'
    write(output_fileunit,'(a)'  ) '# 11       Resultant rotation in x3 direction with respect to the reference point'
    write(output_fileunit,'(a)'  ) '# 12-13    Resultant forces along x1 and x2'
    write(output_fileunit,'(a)'  ) '# 14       Resultant moment in x3 direction with respect to the reference point'
    case (3)
    write(output_fileunit,'(a)'  ) '# 7-9      Reference point'
    write(output_fileunit,'(a)'  ) '# 10-12    Resultant displacements along x1, x2 and x3'
    write(output_fileunit,'(a)'  ) '# 13-15    Resultant rotations in x1, x2 and x3 directions with respect to the reference point'
    write(output_fileunit,'(a)'  ) '# 16-18    Resultant forces along x1, x2 and x3'
    write(output_fileunit,'(a)'  ) '# 19-21    Resultant moments in x1, x2 and x3 directions with respect to the reference point'
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
    do k=1,ncint-2
      write(output_fileunit,'(a)',advance='no') '_'
    end do
    write(output_fileunit,'(a)',advance='no') '1'
    ! Rest of the columns
    select case (problem%n)
      case (2)
        ncmax=14
      case (3)
        ncmax=21
    end select
    do kc=2,ncmax
      ! Depending if integer or real column
      if (kc.le.6) then
        nc=ncint
      else
        nc=ncreal
      end if
      ! Write
      do k=1,nc-fbem_nchar_int(kc)
        write(output_fileunit,'(a)',advance='no') '_'
      end do
      write(fmt1,*) '(i',fbem_nchar_int(kc),')'
      call fbem_trimall(fmt1)
      write(output_fileunit,fmt1,advance='no') kc
    end do
    write(output_fileunit,*)
  end if

  ! Export formats
  ! --------------

  ! Export format for columns 1-6
  write(fmt1,*) '(6',fmt_integer,')'
  call fbem_trim2b(fmt1)
  ! Export format for columns >=7
  write(fmt2,*) '(',fmt_real,')'
  call fbem_trim2b(fmt2)
  ! Export format for vectors
  write(fmt3,*) '(',problem%n,fmt_real,')'
  call fbem_trim2b(fmt3)

  ! -------------------
  ! Calculate and write
  ! -------------------

  ! Loop through regions
  do kr=1,n_regions
    if (region(kr)%class.eq.fbem_be) then
      ! Loop through the BOUNDARIES of each REGION
      do kb=1,region(kr)%n_boundaries
        D=0.
        G=0.
        F=0.
        M=0.
        !
        ! General setup
        !
        sb=region(kr)%boundary(kb)
        sb_reversion=region(kr)%boundary_reversion(kb)
        select case (boundary(sb)%class)
          case (fbem_boundary_class_ordinary)
            nfaces=1
            select case (boundary(sb)%n_regions)
              case (1)
                face=1
              case (2)
                if (sb_reversion) then
                  face=2
                else
                  face=1
                end if
            end select
          case (fbem_boundary_class_cracklike)
            nfaces=2
        end select
        !
        ! Point with respect moments are calculated, xm.
        !
        xm=0.d0
        A=0.d0
        do ke=1,part(boundary(sb)%part)%n_elements
          se=part(boundary(sb)%part)%element(ke)
          xm=xm+element(se)%centroid*element(se)%size
          A=A+element(se)%size
        end do
        xm=xm/A
        if (tot_apply_symmetry.and.(n_symplanes.gt.0)) then
          xmt=0.d0
          do ks=1,n_symelements
            call fbem_symmetry_multipliers(ks,3,n_symplanes,symplane_m,symplane_s,symplane_t,symplane_r,&
                                           symconf_m,symconf_s,symconf_t,symconf_r,reversed)
            xmt=xmt+symconf_m(1:problem%n)*xm
          end do
          xm=xmt/dble(n_symelements)
        end if
        ! If tot_xm=0, then the origin is used.
        if (tot_xm.eq.0) xm=0.d0
        !
        ! Loop through the ELEMENTS of each BOUNDARY
        !
        do ke=1,part(boundary(sb)%part)%n_elements
          se=part(boundary(sb)%part)%element(ke)
          allocate (xi(element(se)%n_dimension),pphi(fbem_n_nodes(element(se)%type_f1)),sphi(fbem_n_nodes(element(se)%type_f2)))
          !
          ! Integrate
          !
          select case (boundary(sb)%class)

            ! =================
            ! ORDINARY BOUNDARY
            ! =================

            case (fbem_boundary_class_ordinary)
              select case (element(se)%n_dimension)
                case (1)
                  rule=element(se)%n_phi
                  do j=1,gl11_n(rule)
                    xi=gl11_xi(j,rule)
                    w=gl11_w(j,rule)
                    x=fbem_position(problem%n,element(se)%type_g,element(se)%x_gn,xi)
                    r=x-xm
                    jg=fbem_jacobian2d(element(se)%type_g,element(se)%x_gn,xi(1))
                    pphi=fbem_phi_hybrid(element(se)%type_f1,element(se)%delta_f,xi)
                    sphi=fbem_phi_hybrid(element(se)%type_f2,element(se)%delta_f,xi)
                    u=0.
                    do kn=1,fbem_n_nodes(element(se)%type_f1)
                      sn=element(se)%node(kn)
                      u=u+pphi(kn)*node(sn)%value_r(1:problem%n,face)
                    end do
                    t=0.
                    do kn=1,fbem_n_nodes(element(se)%type_f2)
                      sn=element(se)%node(kn)
                      t=t+sphi(kn)*node(sn)%value_r((problem%n+1):(2*problem%n),face)
                    end do
                    D(:,face)=D(:,face)+u*w*jg
                    r3=sqrt(r(1)**2+r(2)**2)
                    if (r3/element(se)%csize.gt.1.d-12) G(1,face)=G(1,face)+(r(1)*u(2)-r(2)*u(1))/r3**2*w*jg
                    F(:,face)=F(:,face)+t*w*jg
                    M(1,face)=M(1,face)+(r(1)*t(2)-r(2)*t(1))*w*jg
                  end do
                case (2)
                  select case (fbem_n_edges(element(se)%type_g))
                    case (3)
                      rule=2*element(se)%n_phi-1
                      do j=1,wantri_n(rule)
                        xi(1)=wantri_xi1(j,rule)
                        xi(2)=wantri_xi2(j,rule)
                        w=wantri_w(j,rule)
                        x=fbem_position(problem%n,element(se)%type_g,element(se)%x_gn,xi)
                        r=x-xm
                        jg=fbem_jacobian3d(element(se)%type_g,element(se)%x_gn,xi)
                        pphi=fbem_phi_hybrid(element(se)%type_f1,element(se)%delta_f,xi)
                        sphi=fbem_phi_hybrid(element(se)%type_f2,element(se)%delta_f,xi)
                        u=0.
                        do kn=1,fbem_n_nodes(element(se)%type_f1)
                          sn=element(se)%node(kn)
                          u=u+pphi(kn)*node(sn)%value_r(1:problem%n,face)
                        end do
                        t=0.
                        do kn=1,fbem_n_nodes(element(se)%type_f2)
                          sn=element(se)%node(kn)
                          t=t+sphi(kn)*node(sn)%value_r((problem%n+1):(2*problem%n),face)
                        end do
                        D(:,face)=D(:,face)+u*w*jg
                        r1=sqrt(r(2)**2+r(3)**2)
                        r2=sqrt(r(3)**2+r(1)**2)
                        r3=sqrt(r(1)**2+r(2)**2)
                        if (r1/element(se)%csize.gt.1.d-12) G(1,face)=G(1,face)+(r(2)*u(3)-r(3)*u(2))/r1**2*w*jg
                        if (r2/element(se)%csize.gt.1.d-12) G(2,face)=G(2,face)+(r(3)*u(1)-r(1)*u(3))/r2**2*w*jg
                        if (r3/element(se)%csize.gt.1.d-12) G(3,face)=G(3,face)+(r(1)*u(2)-r(2)*u(1))/r3**2*w*jg
                        F(:,face)=F(:,face)+t*w*jg
                        M(1,face)=M(1,face)+(r(2)*t(3)-r(3)*t(2))*w*jg
                        M(2,face)=M(2,face)+(r(3)*t(1)-r(1)*t(3))*w*jg
                        M(3,face)=M(3,face)+(r(1)*t(2)-r(2)*t(1))*w*jg
                      end do
                    case (4)
                      rule=element(se)%n_phi
                      do j1=1,gl11_n(rule)
                        do j2=1,gl11_n(rule)
                          xi(1)=gl11_xi(j1,rule)
                          xi(2)=gl11_xi(j2,rule)
                          w=gl11_w(j1,rule)*gl11_w(j2,rule)
                          x=fbem_position(problem%n,element(se)%type_g,element(se)%x_gn,xi)
                          r=x-xm
                          jg=fbem_jacobian3d(element(se)%type_g,element(se)%x_gn,xi)
                          pphi=fbem_phi_hybrid(element(se)%type_f1,element(se)%delta_f,xi)
                          sphi=fbem_phi_hybrid(element(se)%type_f2,element(se)%delta_f,xi)
                          u=0.
                          do kn=1,fbem_n_nodes(element(se)%type_f1)
                            sn=element(se)%node(kn)
                            u=u+pphi(kn)*node(sn)%value_r(1:problem%n,face)
                          end do
                          t=0.
                          do kn=1,fbem_n_nodes(element(se)%type_f2)
                            sn=element(se)%node(kn)
                            t=t+sphi(kn)*node(sn)%value_r((problem%n+1):(2*problem%n),face)
                          end do
                          D(:,face)=D(:,face)+u*w*jg
                          r1=sqrt(r(2)**2+r(3)**2)
                          r2=sqrt(r(3)**2+r(1)**2)
                          r3=sqrt(r(1)**2+r(2)**2)
                          if (r1/element(se)%csize.gt.1.d-12) G(1,face)=G(1,face)+(r(2)*u(3)-r(3)*u(2))/r1**2*w*jg
                          if (r2/element(se)%csize.gt.1.d-12) G(2,face)=G(2,face)+(r(3)*u(1)-r(1)*u(3))/r2**2*w*jg
                          if (r3/element(se)%csize.gt.1.d-12) G(3,face)=G(3,face)+(r(1)*u(2)-r(2)*u(1))/r3**2*w*jg
                          F(:,face)=F(:,face)+t*w*jg
                          M(1,face)=M(1,face)+(r(2)*t(3)-r(3)*t(2))*w*jg
                          M(2,face)=M(2,face)+(r(3)*t(1)-r(1)*t(3))*w*jg
                          M(3,face)=M(3,face)+(r(1)*t(2)-r(2)*t(1))*w*jg
                        end do
                      end do
                  end select
              end select

            ! ===================
            ! CRACK-LIKE BOUNDARY
            ! ===================

            case (fbem_boundary_class_cracklike)
              select case (element(se)%n_dimension)
                case (1)
                  rule=element(se)%n_phi
                  do j=1,gl11_n(rule)
                    xi=gl11_xi(j,rule)
                    w=gl11_w(j,rule)
                    x=fbem_position(problem%n,element(se)%type_g,element(se)%x_gn,xi(1))
                    r=x-xm
                    n=fbem_unormal2d(element(se)%type_g,element(se)%x_gn,xi(1))
                    if (sb_reversion) n=-n
                    jg=fbem_jacobian2d(element(se)%type_g,element(se)%x_gn,xi(1))
                    pphi=fbem_phi_hybrid(element(se)%type_f1,element(se)%delta_f,xi)
                    sphi=fbem_phi_hybrid(element(se)%type_f2,element(se)%delta_f,xi)
                    do face=1,2
                      if (face.eq.2) n=-n
                      ! u_k
                      u=0.
                      do kn=1,fbem_n_nodes(element(se)%type_f1)
                        sn=element(se)%node(kn)
                        u=u+pphi(kn)*node(sn)%value_r(1:problem%n,face)
                      end do
                      ! t_k
                      t=0.
                      do kn=1,fbem_n_nodes(element(se)%type_f2)
                        sn=element(se)%node(kn)
                        t=t+sphi(kn)*node(sn)%value_r((problem%n+1):(2*problem%n),face)
                      end do
                      D(:,face)=D(:,face)+u*w*jg
                      r3=sqrt(r(1)**2+r(2)**2)
                      if (r3/element(se)%csize.gt.1.d-12) G(1,face)=G(1,face)+(r(1)*u(2)-r(2)*u(1))/r3**2*w*jg
                      F(:,face)=F(:,face)+t*w*jg
                      M(1,face)=M(1,face)+(r(1)*t(2)-r(2)*t(1))*w*jg
                    end do
                  end do
                case (2)
                  select case (fbem_n_edges(element(se)%type_g))
                    case (3)
                      rule=2*element(se)%n_phi-1
                      do j=1,wantri_n(rule)
                        xi(1)=wantri_xi1(j,rule)
                        xi(2)=wantri_xi2(j,rule)
                        w=wantri_w(j,rule)
                        x=fbem_position(problem%n,element(se)%type_g,element(se)%x_gn,xi)
                        r=x-xm
                        n=fbem_unormal3d(element(se)%type_g,element(se)%x_gn,xi)
                        if (sb_reversion) n=-n
                        jg=fbem_jacobian3d(element(se)%type_g,element(se)%x_gn,xi)
                        pphi=fbem_phi_hybrid(element(se)%type_f1,element(se)%delta_f,xi)
                        sphi=fbem_phi_hybrid(element(se)%type_f2,element(se)%delta_f,xi)
                        do face=1,2
                          if (face.eq.2) n=-n
                          ! u_k
                          u=0.
                          do kn=1,fbem_n_nodes(element(se)%type_f1)
                            sn=element(se)%node(kn)
                            u=u+pphi(kn)*node(sn)%value_r(1:problem%n,face)
                          end do
                          ! t_k
                          t=0.
                          do kn=1,fbem_n_nodes(element(se)%type_f2)
                            sn=element(se)%node(kn)
                            t=t+sphi(kn)*node(sn)%value_r((problem%n+1):(2*problem%n),face)
                          end do
                          D(:,face)=D(:,face)+u*w*jg
                          r1=sqrt(r(2)**2+r(3)**2)
                          r2=sqrt(r(3)**2+r(1)**2)
                          r3=sqrt(r(1)**2+r(2)**2)
                          if (r1/element(se)%csize.gt.1.d-12) G(1,face)=G(1,face)+(r(2)*u(3)-r(3)*u(2))/r1**2*w*jg
                          if (r2/element(se)%csize.gt.1.d-12) G(2,face)=G(2,face)+(r(3)*u(1)-r(1)*u(3))/r2**2*w*jg
                          if (r3/element(se)%csize.gt.1.d-12) G(3,face)=G(3,face)+(r(1)*u(2)-r(2)*u(1))/r3**2*w*jg
                          F(:,face)=F(:,face)+t*w*jg
                          M(1,face)=M(1,face)+(r(2)*t(3)-r(3)*t(2))*w*jg
                          M(2,face)=M(2,face)+(r(3)*t(1)-r(1)*t(3))*w*jg
                          M(3,face)=M(3,face)+(r(1)*t(2)-r(2)*t(1))*w*jg
                        end do
                      end do
                    case (4)
                      rule=element(se)%n_phi
                      do j1=1,gl11_n(rule)
                        do j2=1,gl11_n(rule)
                          xi(1)=gl11_xi(j1,rule)
                          xi(2)=gl11_xi(j2,rule)
                          w=gl11_w(j1,rule)*gl11_w(j2,rule)
                          x=fbem_position(problem%n,element(se)%type_g,element(se)%x_gn,xi)
                          r=x-xm
                          n=fbem_unormal3d(element(se)%type_g,element(se)%x_gn,xi)
                          if (sb_reversion) n=-n
                          jg=fbem_jacobian3d(element(se)%type_g,element(se)%x_gn,xi)
                          pphi=fbem_phi_hybrid(element(se)%type_f1,element(se)%delta_f,xi)
                          sphi=fbem_phi_hybrid(element(se)%type_f2,element(se)%delta_f,xi)
                          do face=1,2
                            if (face.eq.2) n=-n
                            ! u_k
                            u=0.
                            do kn=1,fbem_n_nodes(element(se)%type_f1)
                              sn=element(se)%node(kn)
                              u=u+pphi(kn)*node(sn)%value_r(1:problem%n,face)
                            end do
                            ! t_k
                            t=0.
                            do kn=1,fbem_n_nodes(element(se)%type_f2)
                              sn=element(se)%node(kn)
                              t=t+sphi(kn)*node(sn)%value_r((problem%n+1):(2*problem%n),face)
                            end do
                            D(:,face)=D(:,face)+u*w*jg
                            r1=sqrt(r(2)**2+r(3)**2)
                            r2=sqrt(r(3)**2+r(1)**2)
                            r3=sqrt(r(1)**2+r(2)**2)
                            if (r1/element(se)%csize.gt.1.d-12) G(1,face)=G(1,face)+(r(2)*u(3)-r(3)*u(2))/r1**2*w*jg
                            if (r2/element(se)%csize.gt.1.d-12) G(2,face)=G(2,face)+(r(3)*u(1)-r(1)*u(3))/r2**2*w*jg
                            if (r3/element(se)%csize.gt.1.d-12) G(3,face)=G(3,face)+(r(1)*u(2)-r(2)*u(1))/r3**2*w*jg
                            F(:,face)=F(:,face)+t*w*jg
                            M(1,face)=M(1,face)+(r(2)*t(3)-r(3)*t(2))*w*jg
                            M(2,face)=M(2,face)+(r(3)*t(1)-r(1)*t(3))*w*jg
                            M(3,face)=M(3,face)+(r(1)*t(2)-r(2)*t(1))*w*jg
                          end do
                        end do
                      end do
                  end select
              end select
          end select
          deallocate (xi,pphi,sphi)
        end do
        ! Average kinematics
        D=D/A
        G=G/A
        !
        ! Apply symmetry conditions
        !
        if (tot_apply_symmetry.and.(n_symplanes.gt.0)) then
          DT=0.
          GT=0.
          FT=0.
          MT=0.
          ! Variables
          do ks=1,n_symelements
            ! Setup
            call fbem_symmetry_multipliers(ks,3,n_symplanes,symplane_m,symplane_s,symplane_t,symplane_r,&
                                           symconf_m,symconf_s,symconf_t,symconf_r,reversed)
            ! Values of variables at each symmetry step
            ! DISPLACEMENTS
            do k=1,problem%n
              DS(k,:)=symconf_t(k)*D(k,:)
            end do
            ! ROTATIONS
            select case (problem%n)
              case (2)
                GS(1,:)=symconf_r(3)*G(1,:)
              case (3)
                do k=1,problem%n
                  GS(k,:)=symconf_r(k)*G(k,:)
                end do
            end select
            ! FORCES
            do k=1,problem%n
              FS(k,:)=symconf_t(k)*F(k,:)
            end do
            ! MOMENTS
            select case (problem%n)
              case (2)
                MS(1,:)=symconf_r(3)*M(1,:)
              case (3)
                do k=1,problem%n
                  MS(k,:)=symconf_r(k)*M(k,:)
                end do
            end select
            ! Add all contributions
            ! DISPLACEMENTS
            do k=1,problem%n
              DT(k,:)=DT(k,:)+DS(k,:)
            end do
            ! ROTATIONS
            select case (problem%n)
              case (2)
                GT(1,:)=GT(1,:)+GS(1,:)
              case (3)
                GT(1,:)=GT(1,:)+GS(1,:)
                GT(2,:)=GT(2,:)+GS(2,:)
                GT(3,:)=GT(3,:)+GS(3,:)
            end select
            ! FORCES
            do k=1,problem%n
              FT(k,:)=FT(k,:)+FS(k,:)
            end do
            ! MOMENTS
            select case (problem%n)
              case (2)
                MT(1,:)=MT(1,:)+MS(1,:)
              case (3)
                MT(1,:)=MT(1,:)+MS(1,:)
                MT(2,:)=MT(2,:)+MS(2,:)
                MT(3,:)=MT(3,:)+MS(3,:)
            end select
          end do
          DT=DT/dble(n_symelements)
          GT=GT/dble(n_symelements)
        else
          DT=D
          GT=G
          FT=F
          MT=M
        end if
        !
        ! Write the resultants of each face of the boundary
        !
        do kface=1,2
          if ((nfaces.eq.2).or.((nfaces.eq.1).and.(kface.eq.face))) then
            ! Write columns 1-8
            write(output_fileunit,fmt1,advance='no') region(kr)%id, region(kr)%class, region(kr)%type, &
                                                     boundary(sb)%id, boundary(sb)%class, kface
            ! Write reference point
            write(output_fileunit,fmt3,advance='no') xm
            ! Write columns (Resultant of displacements)
            do k=1,problem%n
              write(output_fileunit,fmt2,advance='no') DT(k,kface)
            end do
            ! Write columns (Resultant of rotations)
            do k=1,2*problem%n-3
              write(output_fileunit,fmt2,advance='no') GT(k,kface)
            end do
            ! Write columns (Resultant of forces)
            do k=1,problem%n
              write(output_fileunit,fmt2,advance='no') FT(k,kface)
            end do
            ! Write columns (Resultant of moments)
            do k=1,2*problem%n-3
              write(output_fileunit,fmt2,advance='no') MT(k,kface)
            end do
            write(output_fileunit,*)
          end if
        end do

      end do ! Loop through the BOUNDARIES of the REGION

      ! ============
      ! BE BODY LOAD
      ! ============

      ! Loop through the BE LOADS of each REGION
      do kb=1,region(kr)%n_be_bodyloads
        D=0.
        G=0.
        F=0.
        M=0.
        !
        ! General setup
        !
        sb=region(kr)%be_bodyload(kb)
        nfaces=1
        face=1
        !
        ! Point with respect moments are calculated, xm.
        !
        xm=0.d0
        A=0.d0
        do ke=1,part(be_bodyload(sb)%part)%n_elements
          se=part(be_bodyload(sb)%part)%element(ke)
          xm=xm+element(se)%centroid*element(se)%size
          A=A+element(se)%size
        end do
        xm=xm/A
        if (tot_apply_symmetry.and.(n_symplanes.gt.0)) then
          xmt=0.d0
          do ks=1,n_symelements
            call fbem_symmetry_multipliers(ks,3,n_symplanes,symplane_m,symplane_s,symplane_t,symplane_r,&
                                           symconf_m,symconf_s,symconf_t,symconf_r,reversed)
            xmt=xmt+symconf_m(1:problem%n)*xm
          end do
          xm=xmt/dble(n_symelements)
        end if
        ! If tot_xm=0, then the origin is used.
        if (tot_xm.eq.0) xm=0.d0
        !
        ! Loop through the ELEMENTS of each BE LOAD
        !
        do ke=1,part(be_bodyload(sb)%part)%n_elements
          se=part(be_bodyload(sb)%part)%element(ke)
          allocate (xi(element(se)%n_dimension),pphi(fbem_n_nodes(element(se)%type_f1)),sphi(fbem_n_nodes(element(se)%type_f2)))
          select case (element(se)%n_dimension)
            case (1)
              rule=element(se)%n_phi
              do j=1,gl11_n(rule)
                xi=gl11_xi(j,rule)
                w=gl11_w(j,rule)
                x=fbem_position(problem%n,element(se)%type_g,element(se)%x_gn,xi)
                r=x-xm
                select case (problem%n)
                  case (2)
                    jg=fbem_jacobian2d(element(se)%type_g,element(se)%x_gn,xi(1))
                  case (3)
                    jg=fbem_jacobian3d(element(se)%type_g,element(se)%x_gn,xi(1))
                  case default
                    stop 'error'
                end select
                pphi=fbem_phi_hybrid(element(se)%type_f1,element(se)%delta_f,xi)
                sphi=fbem_phi_hybrid(element(se)%type_f2,element(se)%delta_f,xi)
                u=0.
                do kn=1,fbem_n_nodes(element(se)%type_f1)
                  sn=element(se)%node(kn)
                  u=u+pphi(kn)*node(sn)%value_r(1:problem%n,face)
                end do
                t=0.
                do kn=1,fbem_n_nodes(element(se)%type_f2)
                  sn=element(se)%node(kn)
                  t=t+sphi(kn)*node(sn)%value_r((problem%n+1):(2*problem%n),face)
                end do
                D(:,face)=D(:,face)+u*w*jg
                select case (problem%n)
                  case (2)
                    r3=sqrt(r(1)**2+r(2)**2)
                    if (r3/element(se)%csize.gt.1.d-12) G(1,face)=G(1,face)+(r(1)*u(2)-r(2)*u(1))/r3**2*w*jg
                    F(:,face)=F(:,face)+t*w*jg
                    M(1,face)=M(1,face)+(r(1)*t(2)-r(2)*t(1))*w*jg
                  case (3)
                    r1=sqrt(r(2)**2+r(3)**2)
                    r2=sqrt(r(3)**2+r(1)**2)
                    r3=sqrt(r(1)**2+r(2)**2)
                    if (r1/element(se)%csize.gt.1.d-12) G(1,face)=G(1,face)+(r(2)*u(3)-r(3)*u(2))/r1**2*w*jg
                    if (r2/element(se)%csize.gt.1.d-12) G(2,face)=G(2,face)+(r(3)*u(1)-r(1)*u(3))/r2**2*w*jg
                    if (r3/element(se)%csize.gt.1.d-12) G(3,face)=G(3,face)+(r(1)*u(2)-r(2)*u(1))/r3**2*w*jg
                    F(:,face)=F(:,face)+t*w*jg
                    M(1,face)=M(1,face)+(r(2)*t(3)-r(3)*t(2))*w*jg
                    M(2,face)=M(2,face)+(r(3)*t(1)-r(1)*t(3))*w*jg
                    M(3,face)=M(3,face)+(r(1)*t(2)-r(2)*t(1))*w*jg
                  case default
                    stop 'error'
                end select
              end do
            case (2)
              select case (fbem_n_edges(element(se)%type_g))
                case (3)
                  rule=2*element(se)%n_phi-1
                  do j=1,wantri_n(rule)
                    xi(1)=wantri_xi1(j,rule)
                    xi(2)=wantri_xi2(j,rule)
                    w=wantri_w(j,rule)
                    x=fbem_position(problem%n,element(se)%type_g,element(se)%x_gn,xi)
                    r=x-xm
                    select case (problem%n)
                      case (2)
                        jg=fbem_jacobian2d(element(se)%type_g,element(se)%x_gn,xi)
                      case (3)
                        jg=fbem_jacobian3d(element(se)%type_g,element(se)%x_gn,xi)
                      case default
                        stop 'error'
                    end select
                    pphi=fbem_phi_hybrid(element(se)%type_f1,element(se)%delta_f,xi)
                    sphi=fbem_phi_hybrid(element(se)%type_f2,element(se)%delta_f,xi)
                    u=0.
                    do kn=1,fbem_n_nodes(element(se)%type_f1)
                      sn=element(se)%node(kn)
                      u=u+pphi(kn)*node(sn)%value_r(1:problem%n,face)
                    end do
                    t=0.
                    do kn=1,fbem_n_nodes(element(se)%type_f2)
                      sn=element(se)%node(kn)
                      t=t+sphi(kn)*node(sn)%value_r((problem%n+1):(2*problem%n),face)
                    end do
                    D(:,face)=D(:,face)+u*w*jg
                    select case (problem%n)
                      case (2)
                        r3=sqrt(r(1)**2+r(2)**2)
                        if (r3/element(se)%csize.gt.1.d-12) G(1,face)=G(1,face)+(r(1)*u(2)-r(2)*u(1))/r3**2*w*jg
                        F(:,face)=F(:,face)+t*w*jg
                        M(1,face)=M(1,face)+(r(1)*t(2)-r(2)*t(1))*w*jg
                      case (3)
                        r1=sqrt(r(2)**2+r(3)**2)
                        r2=sqrt(r(3)**2+r(1)**2)
                        r3=sqrt(r(1)**2+r(2)**2)
                        if (r1/element(se)%csize.gt.1.d-12) G(1,face)=G(1,face)+(r(2)*u(3)-r(3)*u(2))/r1**2*w*jg
                        if (r2/element(se)%csize.gt.1.d-12) G(2,face)=G(2,face)+(r(3)*u(1)-r(1)*u(3))/r2**2*w*jg
                        if (r3/element(se)%csize.gt.1.d-12) G(3,face)=G(3,face)+(r(1)*u(2)-r(2)*u(1))/r3**2*w*jg
                        F(:,face)=F(:,face)+t*w*jg
                        M(1,face)=M(1,face)+(r(2)*t(3)-r(3)*t(2))*w*jg
                        M(2,face)=M(2,face)+(r(3)*t(1)-r(1)*t(3))*w*jg
                        M(3,face)=M(3,face)+(r(1)*t(2)-r(2)*t(1))*w*jg
                      case default
                        stop 'error'
                    end select
                  end do
                case (4)
                  rule=element(se)%n_phi
                  do j1=1,gl11_n(rule)
                    do j2=1,gl11_n(rule)
                      xi(1)=gl11_xi(j1,rule)
                      xi(2)=gl11_xi(j2,rule)
                      w=gl11_w(j1,rule)*gl11_w(j2,rule)
                      x=fbem_position(problem%n,element(se)%type_g,element(se)%x_gn,xi)
                      r=x-xm
                      select case (problem%n)
                        case (2)
                          jg=fbem_jacobian2d(element(se)%type_g,element(se)%x_gn,xi)
                        case (3)
                          jg=fbem_jacobian3d(element(se)%type_g,element(se)%x_gn,xi)
                        case default
                          stop 'error'
                      end select
                      pphi=fbem_phi_hybrid(element(se)%type_f1,element(se)%delta_f,xi)
                      sphi=fbem_phi_hybrid(element(se)%type_f2,element(se)%delta_f,xi)
                      u=0.
                      do kn=1,fbem_n_nodes(element(se)%type_f1)
                        sn=element(se)%node(kn)
                        u=u+pphi(kn)*node(sn)%value_r(1:problem%n,face)
                      end do
                      t=0.
                      do kn=1,fbem_n_nodes(element(se)%type_f2)
                        sn=element(se)%node(kn)
                        t=t+sphi(kn)*node(sn)%value_r((problem%n+1):(2*problem%n),face)
                      end do
                      D(:,face)=D(:,face)+u*w*jg
                      select case (problem%n)
                        case (2)
                          r3=sqrt(r(1)**2+r(2)**2)
                          if (r3/element(se)%csize.gt.1.d-12) G(1,face)=G(1,face)+(r(1)*u(2)-r(2)*u(1))/r3**2*w*jg
                          F(:,face)=F(:,face)+t*w*jg
                          M(1,face)=M(1,face)+(r(1)*t(2)-r(2)*t(1))*w*jg
                        case (3)
                          r1=sqrt(r(2)**2+r(3)**2)
                          r2=sqrt(r(3)**2+r(1)**2)
                          r3=sqrt(r(1)**2+r(2)**2)
                          if (r1/element(se)%csize.gt.1.d-12) G(1,face)=G(1,face)+(r(2)*u(3)-r(3)*u(2))/r1**2*w*jg
                          if (r2/element(se)%csize.gt.1.d-12) G(2,face)=G(2,face)+(r(3)*u(1)-r(1)*u(3))/r2**2*w*jg
                          if (r3/element(se)%csize.gt.1.d-12) G(3,face)=G(3,face)+(r(1)*u(2)-r(2)*u(1))/r3**2*w*jg
                          F(:,face)=F(:,face)+t*w*jg
                          M(1,face)=M(1,face)+(r(2)*t(3)-r(3)*t(2))*w*jg
                          M(2,face)=M(2,face)+(r(3)*t(1)-r(1)*t(3))*w*jg
                          M(3,face)=M(3,face)+(r(1)*t(2)-r(2)*t(1))*w*jg
                        case default
                          stop 'error'
                      end select
                    end do
                  end do
              end select
          end select
          deallocate (xi,pphi,sphi)
        end do
        ! Average kinematics
        D=D/A
        G=G/A
        !
        ! Apply symmetry conditions
        !
        if (tot_apply_symmetry.and.(n_symplanes.gt.0)) then
          DT=0.
          GT=0.
          FT=0.
          MT=0.
          ! Variables
          do ks=1,n_symelements
            ! Setup
            call fbem_symmetry_multipliers(ks,3,n_symplanes,symplane_m,symplane_s,symplane_t,symplane_r,&
                                           symconf_m,symconf_s,symconf_t,symconf_r,reversed)
            ! Values of variables at each symmetry step
            ! DISPLACEMENTS
            do k=1,problem%n
              DS(k,:)=symconf_t(k)*D(k,:)
            end do
            ! ROTATIONS
            select case (problem%n)
              case (2)
                GS(1,:)=symconf_r(3)*G(1,:)
              case (3)
                do k=1,problem%n
                  GS(k,:)=symconf_r(k)*G(k,:)
                end do
            end select
            ! FORCES
            do k=1,problem%n
              FS(k,:)=symconf_t(k)*F(k,:)
            end do
            ! MOMENTS
            select case (problem%n)
              case (2)
                MS(1,:)=symconf_r(3)*M(1,:)
              case (3)
                do k=1,problem%n
                  MS(k,:)=symconf_r(k)*M(k,:)
                end do
            end select
            ! Add all contributions
            ! DISPLACEMENTS
            do k=1,problem%n
              DT(k,:)=DT(k,:)+DS(k,:)
            end do
            ! ROTATIONS
            select case (problem%n)
              case (2)
                GT(1,:)=GT(1,:)+GS(1,:)
              case (3)
                GT(1,:)=GT(1,:)+GS(1,:)
                GT(2,:)=GT(2,:)+GS(2,:)
                GT(3,:)=GT(3,:)+GS(3,:)
            end select
            ! FORCES
            do k=1,problem%n
              FT(k,:)=FT(k,:)+FS(k,:)
            end do
            ! MOMENTS
            select case (problem%n)
              case (2)
                MT(1,:)=MT(1,:)+MS(1,:)
              case (3)
                MT(1,:)=MT(1,:)+MS(1,:)
                MT(2,:)=MT(2,:)+MS(2,:)
                MT(3,:)=MT(3,:)+MS(3,:)
            end select
          end do
          DT=DT/dble(n_symelements)
          GT=GT/dble(n_symelements)
        else
          DT=D
          GT=G
          FT=F
          MT=M
        end if
        !
        ! Write the resultants of each face of the BE LOAD
        !
        ! Export format for columns 1-6
        read(fmt_integer(2:len_trim(fmt_integer)),*) kface
        write(fmt1,*) '(6',fmt_integer,')'
        call fbem_trim2b(fmt1)
        ! Write columns 1-8
        write(output_fileunit,fmt1,advance='no') region(kr)%id, region(kr)%class, region(kr)%type, &
                                                 be_bodyload(sb)%id, 0, 0
        ! Write reference point
        write(output_fileunit,fmt3,advance='no') xm
        ! Write columns (Resultant of displacements)
        do k=1,problem%n
          write(output_fileunit,fmt2,advance='no') DT(k,1)
        end do
        ! Write columns (Resultant of rotations)
        do k=1,2*problem%n-3
          write(output_fileunit,fmt2,advance='no') GT(k,1)
        end do
        ! Write columns (Resultant of forces)
        do k=1,problem%n
          write(output_fileunit,fmt2,advance='no') FT(k,1)
        end do
        ! Write columns (Resultant of moments)
        do k=1,2*problem%n-3
          write(output_fileunit,fmt2,advance='no') MT(k,1)
        end do
        write(output_fileunit,*)


      end do ! Loop through the BE LOADS of the REGION

    end if
  end do ! Loop through the regions

end subroutine export_solution_mechanics_static_tot
