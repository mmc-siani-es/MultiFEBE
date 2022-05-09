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

subroutine export_solution_mechanics_harmonic_ier(kf,output_fileunit)

  ! Fortran 2003 intrinsic module
  use iso_fortran_env

  ! fbem modules
  use fbem_data_structures
  use fbem_string_handling
  use fbem_numerical
  use fbem_geometry
  use fbem_shape_functions
  use fbem_symmetry
  use fbem_quad_rules
  use fbem_fem_shells

  ! Module of problem variables
  use problem_variables

  ! No implicit variables allowed
  implicit none

  ! I/O variables
  integer                                 :: kf
  integer                                 :: output_fileunit
  ! Local variables
  real(kind=real64)                       :: omega
  integer                                 :: j, j1, j2, k
  integer                                 :: kr, kb, ke, kn, knl, sip, kc, kp, sp, kdof
  integer                                 :: sr, sb, se, sn, ks, ss, ks2, ss2, ke2, se2
  integer                                 :: etype
  character(len=fbem_filename_max_length) :: tmp_filename                ! Temporary file name
  character(len=fbem_fmtstr)              :: fmt1, fmt2, fmt3, fmtstr    ! String used for write format string
  integer                                 :: ncint, ncreal, nc, ncmax
  character(len=fbem_string_max_length)   :: tmp_string
  integer                                 :: face, kface, nfaces
  real(kind=real64)                       :: nu
  complex(kind=real64)                    :: mu, kappa, E
  ! SIF
  integer                                 :: midnode
  complex(kind=real64)                    :: delta_u_I_midnode, delta_u_II_midnode
  real(kind=real64)                       :: r_midnode
  integer                                 :: endnode
  complex(kind=real64)                    :: delta_u_I_endnode, delta_u_II_endnode
  real(kind=real64)                       :: r_endnode
  real(kind=real64)                       :: ratio
  real(kind=real64)                       :: n_tip(problem%n), t1_tip(problem%n), t2_tip(problem%n)
  complex(kind=real64)                    :: t_I, t_II
  complex(kind=real64)                    :: K_I, K_II
  ! RESULTANTS
  integer                                 :: rule, method, info, nex, se_int_n_dof, se_int_n_nodes
  integer, allocatable                    :: se_int_n_dof_node(:)
  real(kind=real64), allocatable          :: xi(:), pphi(:), sphi(:), Nsf(:,:), Ssigma_i(:,:), Lc(:,:)
  real(kind=real64)                       :: x(problem%n), n(problem%n), xm(problem%n), xmt(problem%n), xi2d(2), xi3d(3)
  real(kind=real64)                       :: w, jg, r(problem%n), r1, r2, r3, rmin, dmin
  real(kind=real64)                       :: A
  complex(kind=real64), allocatable       :: a_g(:), a_l(:)
  complex(kind=real64)                    :: u(problem%n), t(problem%n), sigma(6)
  complex(kind=real64)                    :: D (problem%n,2), G (2*problem%n-3,2) ! Kinematic resultants (displacements and rotation)
  complex(kind=real64)                    :: F (problem%n,2), M (2*problem%n-3,2) ! Stress resultants (forces and moments)
  complex(kind=real64)                    :: DS(problem%n,2), GS(2*problem%n-3,2) ! Kinematic resultants (displacements and rotation)
  complex(kind=real64)                    :: FS(problem%n,2), MS(2*problem%n-3,2) ! Stress resultants (forces and moments)
  complex(kind=real64)                    :: DT(problem%n,2), GT(2*problem%n-3,2) ! Kinematic resultants (displacements and rotation)
  complex(kind=real64)                    :: FT(problem%n,2), MT(2*problem%n-3,2) ! Stress resultants (forces and moments)
  logical                                 :: sb_reversion, reversed
  real(kind=real64)                       :: symconf_m(3), symconf_s, symconf_t(3), symconf_r(3)

  ! Frequency
  omega=frequency(kf)
  if (frequency_units.eq.'f') omega=omega*c_1_2pi

  ! Falta el calculo de esto para fluidos y medio poroso

  !
  ! Meter una opcion para que considere o no la simetria en este calculo
  !


  !
  ! Ojo!: esto no valido porque ahora
  !

    ! 'ojo: esto no valido porque ahora se utilizan elementos internos BEM discontinuos para....'


  ! ------------------
  ! Header of the file
  ! ------------------

  if ((kf.eq.1).and.export_overwrite) then
    ! Info about the program that generated the file and its general description
    write(output_fileunit,'(a)'  ) '# Program      : multifebe'
    write(output_fileunit,'(a)'  ) '# Version      : 0.6'
    write(output_fileunit,'(a)'  ) '# File_format  : ier'
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
    write(output_fileunit,'(a)'  ) '# 1-2      Frequency index and value.'
    write(output_fileunit,'(a)'  ) '# 3        Part id'
    select case (problem%n)
    case (2)
    write(output_fileunit,'(a)'  ) '# 4-5      Centroid'
    write(output_fileunit,'(a)'  ) '# 6-9      Resultant displacements along x1 and x2 (complex)'
    write(output_fileunit,'(a)'  ) '# 10-11    Resultant rotation in x3 direction with respect to the centroid (complex)'
    write(output_fileunit,'(a)'  ) '# 12-15    Resultant forces along x1 and x2 (complex)'
    write(output_fileunit,'(a)'  ) '# 16-17    Resultant moment in x3 direction with respect to the centroid (complex)'
    case (3)
    write(output_fileunit,'(a)'  ) '# 4-6      Centroid'
    write(output_fileunit,'(a)'  ) '# 7-12     Resultant displacements along x1, x2 and x3 (complex)'
    write(output_fileunit,'(a)'  ) '# 13-18    Resultant rotations in x1, x2 and x3 directions with respect to the centroid (complex)'
    write(output_fileunit,'(a)'  ) '# 19-24    Resultant forces along x1, x2 and x3 (complex)'
    write(output_fileunit,'(a)'  ) '# 25-30    Resultant moments in x1, x2 and x3 directions with respect to the centroid (complex)'
    end select
    write(output_fileunit,'(a)') '#'
    select case (complex_notation)
      case (1)
        write(output_fileunit,'(a)'  ) '# Complex notation: polar'
      case (2)
        write(output_fileunit,'(a)'  ) '# Complex notation: cartesian'
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
        ncmax=17
      case (3)
        ncmax=30
    end select
    do kc=2,ncmax
      ! Depending if integer or real column
      if (kc.eq.3) then
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

  ! Export format for columns 1-3
  write(fmt1,*) '(',fmt_integer,',',fmt_real,',',fmt_integer,')'
  call fbem_trim2b(fmt1)
  ! Export format for columns >=4
  write(fmt2,*) '(2',fmt_real,')'
  call fbem_trim2b(fmt2)
  ! Export format for vectors
  write(fmt3,*) '(',problem%n,fmt_real,')'
  call fbem_trim2b(fmt3)

  ! -------------------
  ! Calculate and write
  ! -------------------

  ! Loop through parts
  do kp=1,internalelements_mesh%n_parts
    ! Initialize
    D=(0.d0,0.d0)
    G=(0.d0,0.d0)
    F=(0.d0,0.d0)
    M=(0.d0,0.d0)
    !
    ! Calculate the centroid of the part and its area
    !
    xm=0.d0
    A=0.d0
    do ke=1,internalelements_mesh%part(kp)%n_elements
      se=internalelements_mesh%part(kp)%element(ke)
      xm=xm+internalelements_mesh%element(se)%centroid*internalelements_mesh%element(se)%size
      A=A+internalelements_mesh%element(se)%size
    end do
    xm=xm/A
    if (tot_apply_symmetry.and.(n_symplanes.gt.0)) then
      xmt=0.d0
      ! Reference point applying symmetry
      do ks=1,n_symelements
        call fbem_symmetry_multipliers(ks,3,n_symplanes,symplane_m,symplane_s,symplane_t,symplane_r,&
                                       symconf_m,symconf_s,symconf_t,symconf_r,reversed)
        xmt=xmt+symconf_m(1:problem%n)*xm
      end do
      xm=xmt/dble(n_symelements)
    end if
    !
    ! Calculate stress and kinematic resultants over all elements of this part
    !


    select case (region(internalelements_mesh%part(kp)%entity)%class)

      ! ==================================
      ! BOUNDARY ELEMENT INTERNAL ELEMENTS
      ! ==================================

      case (fbem_be)

        select case (region(internalelements_mesh%part(kp)%entity)%type)

          ! --------------
          ! INVISCID FLUID
          ! --------------

          case (fbem_potential)
            stop 'not yet 87'

          ! ------------------
          ! VISCOELASTIC SOLID
          ! ------------------

          case (fbem_viscoelastic)
            ! Loop through the elements of each part
            do ke=1,internalelements_mesh%part(kp)%n_elements
              se=internalelements_mesh%part(kp)%element(ke)
              etype=internalelements_mesh%element(se)%type
              allocate (xi(internalelements_mesh%element(se)%n_dimension))
              select case (internalelements_mesh%element(se)%n_dimension)
                case (1)
                  do j=1,gl11_n(internalelements_order)
                    xi=gl11_xi(j,internalelements_order)
                    w=gl11_w(j,internalelements_order)
                    jg=fbem_jacobian2d(etype,internalelements_mesh%element(se)%x_gn,xi(1))
                    sip=internalelements_mesh%element(se)%internalpoint(j)
                    u=internalpoint(sip)%value_c(:,0)
                    t(1)=sum(internalpoint(sip)%value_c(1,1:problem%n)*internalpoint(sip)%n)
                    t(2)=sum(internalpoint(sip)%value_c(2,1:problem%n)*internalpoint(sip)%n)
                    r=internalpoint(sip)%x-xm
                    D(:,1)=D(:,1)+u*w*jg
                    r3=sqrt(r(1)**2+r(2)**2)
                    if (r3/internalelements_mesh%element(se)%csize.gt.1.d-12) G(1,1)=G(1,1)+(r(1)*u(2)-r(2)*u(1))/r3**2*w*jg
                    F(:,1)=F(:,1)+t(:)*w*jg
                    M(1,1)=M(1,1)+(r(1)*t(2)-r(2)*t(1))*w*jg
                  end do
                case (2)
                  select case (fbem_n_edges(internalelements_mesh%element(se)%type))
                    case (3)
                      do j=1,wantri_n(internalelements_order)
                        xi(1)=wantri_xi1(j,internalelements_order)
                        xi(2)=wantri_xi2(j,internalelements_order)
                        w=wantri_w(j,internalelements_order)
                        jg=fbem_jacobian3d(etype,internalelements_mesh%element(se)%x_gn,xi)
                        sip=internalelements_mesh%element(se)%internalpoint(j)
                        u=internalpoint(sip)%value_c(:,0)
                        t(1)=sum(internalpoint(sip)%value_c(1,1:problem%n)*internalpoint(sip)%n)
                        t(2)=sum(internalpoint(sip)%value_c(2,1:problem%n)*internalpoint(sip)%n)
                        t(3)=sum(internalpoint(sip)%value_c(3,1:problem%n)*internalpoint(sip)%n)
                        r=internalpoint(sip)%x-xm
                        D(:,1)=D(:,1)+u*w*jg
                        r1=sqrt(r(2)**2+r(3)**2)
                        r2=sqrt(r(3)**2+r(1)**2)
                        r3=sqrt(r(1)**2+r(2)**2)
                        if (r1/internalelements_mesh%element(se)%csize.gt.1.d-12) G(1,1)=G(1,1)+(r(2)*u(3)-r(3)*u(2))/r1**2*w*jg
                        if (r2/internalelements_mesh%element(se)%csize.gt.1.d-12) G(2,1)=G(2,1)+(r(3)*u(1)-r(1)*u(3))/r2**2*w*jg
                        if (r3/internalelements_mesh%element(se)%csize.gt.1.d-12) G(3,1)=G(3,1)+(r(1)*u(2)-r(2)*u(1))/r3**2*w*jg
                        F(:,1)=F(:,1)+t(:)*w*jg
                        M(1,1)=M(1,1)+(r(2)*t(3)-r(3)*t(2))*w*jg
                        M(2,1)=M(2,1)+(r(3)*t(1)-r(1)*t(3))*w*jg
                        M(3,1)=M(3,1)+(r(1)*t(2)-r(2)*t(1))*w*jg
                      end do
                    case (4)
                      do j1=1,gl11_n(internalelements_order)
                        do j2=1,gl11_n(internalelements_order)
                          xi(1)=gl11_xi(j1,internalelements_order)
                          xi(2)=gl11_xi(j2,internalelements_order)
                          w=gl11_w(j1,internalelements_order)*gl11_w(j2,internalelements_order)
                          jg=fbem_jacobian3d(etype,internalelements_mesh%element(se)%x_gn,xi)
                          k=j2+(j1-1)*gl11_n(internalelements_order)
                          sip=internalelements_mesh%element(se)%internalpoint(k)
                          u=internalpoint(sip)%value_c(:,0)
                          t(1)=sum(internalpoint(sip)%value_c(1,1:problem%n)*internalpoint(sip)%n)
                          t(2)=sum(internalpoint(sip)%value_c(2,1:problem%n)*internalpoint(sip)%n)
                          t(3)=sum(internalpoint(sip)%value_c(3,1:problem%n)*internalpoint(sip)%n)
                          r=internalpoint(sip)%x-xm
                          D(:,1)=D(:,1)+u*w*jg
                          r1=sqrt(r(2)**2+r(3)**2)
                          r2=sqrt(r(3)**2+r(1)**2)
                          r3=sqrt(r(1)**2+r(2)**2)
                          if (r1/internalelements_mesh%element(se)%csize.gt.1.d-12) G(1,1)=G(1,1)+(r(2)*u(3)-r(3)*u(2))/r1**2*w*jg
                          if (r2/internalelements_mesh%element(se)%csize.gt.1.d-12) G(2,1)=G(2,1)+(r(3)*u(1)-r(1)*u(3))/r2**2*w*jg
                          if (r3/internalelements_mesh%element(se)%csize.gt.1.d-12) G(3,1)=G(3,1)+(r(1)*u(2)-r(2)*u(1))/r3**2*w*jg
                          F(:,1)=F(:,1)+t(:)*w*jg
                          M(1,1)=M(1,1)+(r(2)*t(3)-r(3)*t(2))*w*jg
                          M(2,1)=M(2,1)+(r(3)*t(1)-r(1)*t(3))*w*jg
                          M(3,1)=M(3,1)+(r(1)*t(2)-r(2)*t(1))*w*jg
                        end do
                      end do
                  end select
              end select
              deallocate (xi)
            end do

          ! -----------------
          ! POROELASTIC MEDIA
          ! -----------------

          case (fbem_poroelastic)
            stop 'not yet 88'

        end select





      ! ================================
      ! FINITE ELEMENT INTERNAL ELEMENTS
      ! ================================

      case (fbem_fe)

        ! Loop through the elements of each part
        do ke=1,internalelements_mesh%part(kp)%n_elements
          se=internalelements_mesh%part(kp)%element(ke)
          etype=internalelements_mesh%element(se)%type
          allocate (xi(internalelements_mesh%element(se)%n_dimension))
          select case (internalelements_mesh%element(se)%n_dimension)
            case (1)
              stop 'not yet 89'
              do j=1,gl11_n(internalelements_order)
                xi=gl11_xi(j,internalelements_order)
                w=gl11_w(j,internalelements_order)
                jg=fbem_jacobian2d(etype,internalelements_mesh%element(se)%x_gn,xi(1))

                ! ...

              end do
            case (2)
              select case (fbem_n_edges(internalelements_mesh%element(se)%type))
                case (3)

                  stop 'not yet 90'

                  do j=1,wantri_n(internalelements_order)
                    xi(1)=wantri_xi1(j,internalelements_order)
                    xi(2)=wantri_xi2(j,internalelements_order)
                    w=wantri_w(j,internalelements_order)
                    jg=fbem_jacobian3d(etype,internalelements_mesh%element(se)%x_gn,xi)
                    sip=internalelements_mesh%element(se)%internalpoint(j)
                    u=internalpoint(sip)%value_c(:,0)
                    t(1)=sum(internalpoint(sip)%value_c(1,1:problem%n)*internalpoint(sip)%n)
                    t(2)=sum(internalpoint(sip)%value_c(2,1:problem%n)*internalpoint(sip)%n)
                    t(3)=sum(internalpoint(sip)%value_c(3,1:problem%n)*internalpoint(sip)%n)
                    r=internalpoint(sip)%x-xm
                    D(:,1)=D(:,1)+u*w*jg
                    r1=sqrt(r(2)**2+r(3)**2)
                    r2=sqrt(r(3)**2+r(1)**2)
                    r3=sqrt(r(1)**2+r(2)**2)
                    if (r1/internalelements_mesh%element(se)%csize.gt.1.d-12) G(1,1)=G(1,1)+(r(2)*u(3)-r(3)*u(2))/r1**2*w*jg
                    if (r2/internalelements_mesh%element(se)%csize.gt.1.d-12) G(2,1)=G(2,1)+(r(3)*u(1)-r(1)*u(3))/r2**2*w*jg
                    if (r3/internalelements_mesh%element(se)%csize.gt.1.d-12) G(3,1)=G(3,1)+(r(1)*u(2)-r(2)*u(1))/r3**2*w*jg
                    F(:,1)=F(:,1)+t(:)*w*jg
                    M(1,1)=M(1,1)+(r(2)*t(3)-r(3)*t(2))*w*jg
                    M(2,1)=M(2,1)+(r(3)*t(1)-r(1)*t(3))*w*jg
                    M(3,1)=M(3,1)+(r(1)*t(2)-r(2)*t(1))*w*jg
                  end do
                case (4)
                  do j1=1,gl11_n(internalelements_order)
                    do j2=1,gl11_n(internalelements_order)
                      xi(1)=gl11_xi(j1,internalelements_order)
                      xi(2)=gl11_xi(j2,internalelements_order)
                      w=gl11_w(j1,internalelements_order)*gl11_w(j2,internalelements_order)
                      jg=fbem_jacobian3d(etype,internalelements_mesh%element(se)%x_gn,xi)
                      x=fbem_position(problem%n,etype,internalelements_mesh%element(se)%x_gn,xi)
                      n=fbem_unormal3d(etype,internalelements_mesh%element(se)%x_gn,xi)


                      ! Buscar los elementos de la region que contengan al punto, y hallar u y sigma promedio

                      nex=0
                      u=0
                      sigma=0
                      kr=internalelements_mesh%part(kp)%entity
                      if (region(kr)%subtype.eq.fbem_elastic_orthotropic) then
                        stop 'orthotropic posprocessing not yet'
                      end if
                      E=region(kr)%property_c(5)
                      nu=region(kr)%property_r(3)
                      do ks2=1,region(kr)%n_fe_subregions
                        ss2=region(kr)%fe_subregion(ks2)
                        do ke2=1,part(fe_subregion(ss2)%part)%n_elements
                          se2=part(fe_subregion(ss2)%part)%element(ke2)

                          ! Check if the point is in the element
                          if (.not.allocated(element(se2)%bball_centre)) then
                            allocate(element(se2)%bball_centre(problem%n))
                          end if
                          ! Calculate the element ball
                          call fbem_geometry_element_ball(problem%n,element(se2)%type,element(se2)%x_gn,5,element(se2)%bball_centre,element(se2)%bball_radius)
                          r=element(se2)%bball_centre-x
                          if ((sqrt(dot_product(r,r))-element(se2)%bball_radius).ge.(4.d0*element(se2)%bball_radius)) then
                            cycle
                          end if
                          call fbem_nearest_element_point_bem(problem%n,element(se2)%type,element(se2)%x_gn,2.d0*element(se2)%bball_radius,x,xi2d,rmin,dmin,method)



                          ! solo para shells por ahora
                          if ((problem%n.eq.3).and.(element(se2)%n_dimension.eq.2)) then
                            if (rmin.le.0.51d0*maxval(element(se2)%tv_midnode(3,:))) then
                              xi3d(1)=xi2d(1)
                              xi3d(2)=xi2d(2)
                              xi3d(3)=0.d0
                              call fbem_nearest_minimization_3d_shell(element(se2)%type,element(se2)%x_gn,element(se2)%v_midnode,element(se2)%tv_midnode,x,1.d-13,30,xi3d,info)
                              ! Only valid for points x_i in the element
                              if ((info.ne.1).or.(.not.fbem_check_xi1xi2(element(se2)%type,xi3d(1:2))).or.(.not.fbem_check_xi(xi3d(3)))) then
                                cycle
                              else
                                ! One element found
                                nex=nex+1
                                se_int_n_nodes=fbem_n_nodes(element(se2)%type)
                                allocate (Nsf(3,5*se_int_n_nodes),Ssigma_i(6,5*se_int_n_nodes))
                                call fbem_fem_degshell_N(element(se2)%type,element(se2)%v_midnode,element(se2)%tv_midnode(3,:),xi3d,Nsf)
                                call fbem_fem_degshell_stress_tensor(element(se2)%type,element(se2)%x_gn,element(se2)%v_midnode,element(se2)%tv_midnode(3,:),nu,element(se2)%ksh,xi3d,Ssigma_i)
                                ! Number of DOF of each node of the element
                                allocate (se_int_n_dof_node(se_int_n_nodes))
                                do kn=1,se_int_n_nodes
                                  sn=element(se2)%node(kn)
                                  se_int_n_dof_node(kn)=node(sn)%n_dof
                                end do
                                se_int_n_dof=sum(se_int_n_dof_node)
                                ! Calculate COORDINATE TRANSFORMATION MATRIX (Lc)
                                allocate (Lc(se_int_n_dof,5*se_int_n_nodes))
                                call fbem_fem_degshell_Lc(element(se2)%type,se_int_n_dof,se_int_n_dof_node,element(se2)%v_midnode,Lc)
                                ! Build the element DOF vector with displacements and local rotations
                                allocate (a_g(se_int_n_dof),a_l(5*se_int_n_nodes))
                                kdof=1
                                do kn=1,se_int_n_nodes
                                  sn=element(se2)%node(kn)
                                  a_g((kdof):(kdof-1+se_int_n_dof_node(kn)))=node(sn)%value_c(:,1)
                                  kdof=kdof+se_int_n_dof_node(kn)
                                end do
                                a_l=matmul(transpose(Lc),a_g)
                                u=u+matmul(Nsf,a_l)
                                sigma=sigma+E*matmul(Ssigma_i,a_l)
                                deallocate (Nsf,Ssigma_i,se_int_n_dof_node,Lc,a_g,a_l)
                              end if
                            else

                              cycle
                            end if
                          else
                            stop 'not yet 91'
                          end if
                        end do
                      end do
                      if (nex.eq.0) then
                        stop 'internal element out of any finite element'
                      end if
                      u=u/nex
                      sigma=sigma/nex

                      t(1)=sigma(1)*n(1)+sigma(4)*n(2)+sigma(6)*n(3)
                      t(2)=sigma(4)*n(1)+sigma(2)*n(2)+sigma(5)*n(3)
                      t(3)=sigma(6)*n(1)+sigma(5)*n(2)+sigma(3)*n(3)

                      r=x-xm
                      D(:,1)=D(:,1)+u*w*jg
                      r1=sqrt(r(2)**2+r(3)**2)
                      r2=sqrt(r(3)**2+r(1)**2)
                      r3=sqrt(r(1)**2+r(2)**2)
                      if (r1/internalelements_mesh%element(se)%csize.gt.1.d-12) G(1,1)=G(1,1)+(r(2)*u(3)-r(3)*u(2))/r1**2*w*jg
                      if (r2/internalelements_mesh%element(se)%csize.gt.1.d-12) G(2,1)=G(2,1)+(r(3)*u(1)-r(1)*u(3))/r2**2*w*jg
                      if (r3/internalelements_mesh%element(se)%csize.gt.1.d-12) G(3,1)=G(3,1)+(r(1)*u(2)-r(2)*u(1))/r3**2*w*jg
                      F(:,1)=F(:,1)+t(:)*w*jg
                      M(1,1)=M(1,1)+(r(2)*t(3)-r(3)*t(2))*w*jg
                      M(2,1)=M(2,1)+(r(3)*t(1)-r(1)*t(3))*w*jg
                      M(3,1)=M(3,1)+(r(1)*t(2)-r(2)*t(1))*w*jg
                    end do
                  end do
              end select
          end select
          deallocate (xi)
        end do









    end select








    ! Average kinematics
    D=D/A
    G=G/A
    !
    ! Apply symmetry conditions
    !
    if (tot_apply_symmetry.and.(n_symplanes.gt.0)) then
      DT=(0.d0,0.d0)
      GT=(0.d0,0.d0)
      FT=(0.d0,0.d0)
      MT=(0.d0,0.d0)
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
    ! Write the resultant of each part
    !
    ! Write columns 1-3
    write(output_fileunit,fmt1,advance='no') kf, omega, internalelements_mesh%part(kp)%id
    ! Write centroid
    write(output_fileunit,fmt3,advance='no') xm
    ! Write columns (Resultant of displacements)
    do k=1,problem%n
      select case (complex_notation)
        case (1)
          write(output_fileunit,fmt2,advance='no') abs(DT(k,1)), fbem_zarg(DT(k,1))
        case (2)
          write(output_fileunit,fmt2,advance='no') real(DT(k,1)), imag(DT(k,1))
      end select
    end do
    ! Write columns (Resultant of rotations)
    do k=1,2*problem%n-3
      select case (complex_notation)
        case (1)
          write(output_fileunit,fmt2,advance='no') abs(GT(k,1)), fbem_zarg(GT(k,1))
        case (2)
          write(output_fileunit,fmt2,advance='no') real(GT(k,1)), imag(GT(k,1))
      end select
    end do
    ! Write columns (Resultant of forces)
    do k=1,problem%n
      select case (complex_notation)
        case (1)
          write(output_fileunit,fmt2,advance='no') abs(FT(k,1)), fbem_zarg(FT(k,1))
        case (2)
          write(output_fileunit,fmt2,advance='no') real(FT(k,1)), imag(FT(k,1))
      end select
    end do
    ! Write columns (Resultant of moments)
    do k=1,2*problem%n-3
      select case (complex_notation)
        case (1)
          write(output_fileunit,fmt2,advance='no') abs(MT(k,1)), fbem_zarg(MT(k,1))
        case (2)
          write(output_fileunit,fmt2,advance='no') real(MT(k,1)), imag(MT(k,1))
      end select
    end do
    write(output_fileunit,*)

  end do ! Loop through the parts of the internal elements mesh


end subroutine export_solution_mechanics_harmonic_ier
